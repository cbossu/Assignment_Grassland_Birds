#!/usr/bin/env python
"""
FEEMS on breeding birds + spatial assignment of nonbreeding birds to breeding demes.

Workflow
  1. Read a PLINK bed (Beagle-imputed, FST-outlier subset) holding BOTH breeding and
     nonbreeding birds, plus a metadata CSV saying which is which.
  2. Keep SNPs polymorphic in the breeding set (SpatialGraph silently drops invariant
     SNPs, which would misalign the nonbreeding genotype columns).
  3. Build the grid from breeding coordinates only and fit FEEMS on breeding birds only
     (lambda / lambda_q by cross-validation unless supplied).
  4. Compute FEEMS posterior allele frequencies at every deme (sampled or not) and the
     posterior probability that each nonbreeding bird came from each deme
     (feems.spatial_prediction; point_mu or truncated-normal likelihood).
  5. Optional leave-one-deme-out validation on breeding birds (great-circle error of the
     MAP deme) so you know how much to trust step 4.

Metadata CSV columns (header required):
  sample_id  - must match the IID column of the .fam
  group      - 'breeding' or 'nonbreeding'
  lon, lat   - decimal degrees. Required for breeding birds; for nonbreeding birds these are
               the capture locations and are only carried through to the output.

Coordinates are (lon, lat) everywhere, as FEEMS expects.
"""
import argparse
import json
import os
import sys

import numpy as np
import pandas as pd
from pandas_plink import read_plink

from feems import SpatialGraph
from feems.utils import prepare_graph_inputs
from feems.spatial_graph import query_node_attributes
from feems.cross_validation import run_cv_joint, train_test_split
from feems.spatial_prediction import (
    predict_deme_point_mu,
    predict_deme_trunc_normal_mu,
)


# --------------------------------------------------------------------------- helpers
def haversine_km(lon1, lat1, lon2, lat2):
    lon1, lat1, lon2, lat2 = map(np.radians, (lon1, lat1, lon2, lat2))
    a = (np.sin((lat2 - lat1) / 2) ** 2
         + np.cos(lat1) * np.cos(lat2) * np.sin((lon2 - lon1) / 2) ** 2)
    return 2 * 6371.0 * np.arcsin(np.sqrt(a))


def credible_set_size(p, level=0.95):
    """Number of demes needed to reach `level` posterior mass."""
    s = np.sort(p)[::-1]
    return int(np.searchsorted(np.cumsum(s), level) + 1)


#def predict(g, sp_graph, method):
#    if method == "trunc":
#        z, _ = predict_deme_trunc_normal_mu(g, sp_graph)
#    else:
#        z, _ = predict_deme_point_mu(g, sp_graph)
#    return np.exp(z)  # [samples x demes], columns in permuted node order
from scipy.special import softmax

def predict_ll(g, sp_graph, method):
    if method == "trunc":
        z, _ = predict_deme_trunc_normal_mu(g, sp_graph)
    else:
        z, _ = predict_deme_point_mu(g, sp_graph)
    return np.asarray(z)  # log scale, [samples x demes], permuted node order

def predict(g, sp_graph, method, T=1.0):
    return softmax(predict_ll(g, sp_graph, method) / T, axis=1)

    # ---- calibration-step -------------------------------------------
def calibrate_T(LL, TC, out_csv, Ts=np.geomspace(1, 1e5, 51)):
    cov = []
    for T in Ts:
        P_ = softmax(LL / T, axis=1)
        pt = P_[np.arange(len(TC)), TC]
        cov.append(np.mean([p[p > q].sum() < 0.95 for p, q in zip(P_, pt)]))
    cov = np.asarray(cov)
    pd.DataFrame({"T": Ts, "coverage": cov}).to_csv(out_csv, index=False)
    ok = cov >= 0.95
    return float(Ts[ok.argmax()]) if ok.any() else float("nan")

def deme_table(sp_graph):
    perm = query_node_attributes(sp_graph, "permuted_idx")
    n_samp = query_node_attributes(sp_graph, "n_samples")[perm]
    xy = sp_graph.node_pos[perm]
    return pd.DataFrame({
        "deme": perm,               # original node id (matches nodepos.csv row)
        "deme_lon": xy[:, 0],
        "deme_lat": xy[:, 1],
        "n_breeding_samples": n_samp,
    })


def summarize(P, demes, ids, cap_lon=None, cap_lat=None):
    rows = []
    for i, sid in enumerate(ids):
        p = P[i]
        k = int(np.argmax(p))
        rows.append({
            "sample_id": sid,
            "capture_lon": np.nan if cap_lon is None else cap_lon[i],
            "capture_lat": np.nan if cap_lat is None else cap_lat[i],
            "map_deme": int(demes.deme.iloc[k]),
            "map_lon": demes.deme_lon.iloc[k],
            "map_lat": demes.deme_lat.iloc[k],
            "map_prob": p[k],
            "map_deme_sampled": bool(demes.n_breeding_samples.iloc[k] > 0),
            "n_demes_95cs": credible_set_size(p, 0.95),
            "postmean_lon": np.sum(p * demes.deme_lon.values),
            "postmean_lat": np.sum(p * demes.deme_lat.values),
        })
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bfile", required=True, help="PLINK prefix (bed/bim/fam)")
    ap.add_argument("--meta", required=True, help="metadata CSV (see header)")
    ap.add_argument("--grid", required=True,
                    help="discrete global grid .shp (grid_250 / grid_500 ship with feems)")
    ap.add_argument("--outer", default=None,
                    help="outer boundary file, two columns lon lat (whitespace). "
                         "Default: convex hull of breeding sites + --buffer")
    ap.add_argument("--buffer", type=float, default=2.0,
                    help="degrees added around convex hull when --outer not given")
    ap.add_argument("--translated", action="store_true",
                    help="date-line wrap (only needed if range reaches western Alaska)")
    ap.add_argument("--lamb", type=float, default=None, help="skip CV and use this lambda")
    ap.add_argument("--lamb_q", type=float, default=None, help="skip CV and use this lambda_q")
    ap.add_argument("--cv_folds", type=int, default=None,
                    help="CV folds over sampled demes (default: leave-one-deme-out)")
    ap.add_argument("--method", choices=["point_mu", "trunc"], default="point_mu",
                    help="trunc integrates over posterior frequency uncertainty "
                         "(slower: dense pseudo-inverse of the Laplacian)")
    ap.add_argument("--loo", action="store_true",
                    help="leave-one-deme-out assignment validation on breeding birds")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--temp", type=float, default=1.0,
                    help="temperature for assignment posteriors (calibrate with --loo)")
    ap.add_argument("--loo_ind", action="store_true",
                    help="leave-individual-out validation (stratified K-fold by bird)")
    ap.add_argument("--ind_folds", type=int, default=10,
                    help="folds for --loo_ind; set to the number of breeding birds for true LOO")
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    out = lambda f: os.path.join(args.outdir, f)

    # ---- genotypes -------------------------------------------------------------
    bim, fam, G = read_plink(args.bfile, verbose=False)
    X = np.asarray(G.compute()).T  # samples x snps
    iid = fam["iid"].astype(str).values
    print(f"[read] {X.shape[0]} samples x {X.shape[1]} SNPs from {args.bfile}")

    meta = pd.read_csv(args.meta)
    meta["sample_id"] = meta["sample_id"].astype(str)
    meta["group"] = meta["group"].str.strip().str.lower()
    missing = set(meta.sample_id) - set(iid)
    if missing:
        sys.exit(f"[error] {len(missing)} meta samples not in .fam, e.g. {sorted(missing)[:5]}")
    idx = {s: i for i, s in enumerate(iid)}

    br = meta[meta.group == "breeding"].reset_index(drop=True)
    nb = meta[meta.group == "nonbreeding"].reset_index(drop=True)
    if br[["lon", "lat"]].isna().any().any():
        sys.exit("[error] every breeding bird needs lon/lat")
    Xb = X[[idx[s] for s in br.sample_id]]
    Xn = X[[idx[s] for s in nb.sample_id]]
    print(f"[groups] breeding={len(br)}  nonbreeding={len(nb)}")

    # imputed data should be hard calls; flag dosages / missingness
    for name, M in (("breeding", Xb), ("nonbreeding", Xn)):
        nmiss = np.isnan(M).sum()
        nonint = (~np.isnan(M) & ~np.isclose(M, np.round(M))).sum()
        print(f"[check] {name}: missing={nmiss}  non-integer={nonint}")
    if np.isnan(Xb).any():
        print("[warn] missing breeding genotypes -> SNP-mean imputation (FEEMS needs none)")
        col_mean = np.nanmean(Xb, axis=0)
        r, c = np.where(np.isnan(Xb))
        Xb[r, c] = col_mean[c]
    # nonbreeding NaNs are fine: the assignment likelihood skips them

    # ---- SNPs polymorphic in breeding set (keeps columns aligned) -------------
    ac = Xb.sum(0)
    keep = (ac > 0) & (ac < 2 * Xb.shape[0])
    if (~keep).any():
        print(f"[snps] dropping {(~keep).sum()} SNPs monomorphic among breeding birds")
    Xb, Xn = Xb[:, keep], Xn[:, keep]
    bim.loc[keep].to_csv(out("snps_used.tsv"), sep="\t", index=False)

    # ---- graph -----------------------------------------------------------------
    coord_b = br[["lon", "lat"]].to_numpy(float)
    outer = np.loadtxt(args.outer) if args.outer else None
    outer, edges, grid, _ = prepare_graph_inputs(
        coord=coord_b, ggrid=args.grid, translated=args.translated,
        buffer=args.buffer, outer=outer)
    sp = SpatialGraph(Xb, coord_b, grid, edges, scale_snps=True)
    assert sp.genotypes.shape[1] == Xn.shape[1], "SNP columns misaligned"
    print(f"[graph] {len(sp)} demes, {sp.n_observed_nodes} with breeding samples")

    from math import radians, sin, cos, asin, sqrt
    # ---- indexing check 1: are breeding birds attached to nearby nodes? ----
    def hav(lon1, lat1, lon2, lat2):
        lon1, lat1, lon2, lat2 = map(np.radians, [lon1, lat1, lon2, lat2])
        a = np.sin((lat2-lat1)/2)**2 + np.cos(lat1)*np.cos(lat2)*np.sin((lon2-lon1)/2)**2
        return 6371 * 2 * np.arcsin(np.sqrt(a))

    node_idx = sp.assned_node_idx                     # node each breeding bird is attached to
    npos = sp.node_pos                                # (n_nodes, 2) lon, lat
    d = hav(coord_b[:, 0], coord_b[:, 1], npos[node_idx, 0], npos[node_idx, 1])
    print("[check1] bird-to-own-node km (50/90/100%):", np.percentile(d, [50, 90, 100]))

    e = np.array(sp.edges)  # or list(sp.edges()) depending on version
    print("[grid] median edge km:",
          np.median(hav(sp.node_pos[e[:,0],0], sp.node_pos[e[:,0],1],
                        sp.node_pos[e[:,1],0], sp.node_pos[e[:,1],1])))

    # ---- lambda by CV -----------------------------------------------------------
    if args.lamb is None or args.lamb_q is None:
        lamb_grid = np.geomspace(1e-3, 1e2, 10)[::-1]
        lamb_q_grid = np.geomspace(1, 1e5, 9)[::-1]
        cv_err = run_cv_joint(sp, lamb_grid, lamb_q_grid,
                              n_folds=args.cv_folds, factr=1e10)
        m = np.mean(cv_err, axis=0)              # [lamb_q x lamb]
        iq, il = np.unravel_index(np.nanargmin(m), m.shape)
        lamb, lamb_q = float(lamb_grid[il]), float(lamb_q_grid[iq])
        pd.DataFrame(m, index=lamb_q_grid, columns=lamb_grid).to_csv(out("cv_error.csv"))
        if il in (0, len(lamb_grid) - 1) or iq in (0, len(lamb_q_grid) - 1):
            print("\n[warn] CV optimum on the edge of the grid; widen the grid", flush=True)
    else:
        lamb, lamb_q = args.lamb, args.lamb_q
    print(f"[fit] lambda={lamb:.4g}  lambda_q={lamb_q:.4g}")

    sp.fit_null_model(verbose=False)
    sp.fit(lamb=lamb, lamb_q=lamb_q, verbose=False)

    # exports for feems' ggplot_feems.R
    np.savetxt(out("edgew.csv"),
               np.vstack((np.array(sp.edges).T, sp.w)).T, delimiter=",")
    np.savetxt(out("nodepos.csv"),
               np.vstack((sp.node_pos.T,
                          [sp.nodes[n]["n_samples"] for n in range(len(sp.nodes))])).T,
               delimiter=",")

    # ---- export migration surface --------------------------------------------
    edges = np.array(list(sp.edges()))
    pos = np.array([sp.nodes[n]["pos"] for n in range(len(sp))])
    w = np.asarray(sp.w)
    assert len(w) == len(edges), "edge weights and edge list differ in length"
    lw = np.log10(w)
    pd.DataFrame({
        "node1": edges[:, 0], "node2": edges[:, 1],
        "lon1": pos[edges[:, 0], 0], "lat1": pos[edges[:, 0], 1],
        "lon2": pos[edges[:, 1], 0], "lat2": pos[edges[:, 1], 1],
        "w": w, "log10w_centered": lw - lw.mean(),   # FEEMS's usual display scale
    }).to_csv(out("edge_weights.csv"), index=False)
    print(f"[surface] {len(w)} edges; log10 w range {lw.min():.2f} to {lw.max():.2f}")

    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt, cartopy.crs as ccrs
    from feems import Viz
    proj = ccrs.EquidistantConic(central_longitude=-100, central_latitude=40)
    fig = plt.figure(dpi=300); ax = fig.add_subplot(1, 1, 1, projection=proj)
    v = Viz(ax, sp, projection=proj, edge_width=0.5, edge_alpha=1, edge_zorder=100,
            sample_pt_size=10, obs_node_size=7.5, sample_pt_color="black", cbar_font_size=10)
    v.draw_map(); v.draw_edges(use_weights=True); v.draw_obs_nodes(use_ids=False)
    v.draw_edge_colorbar()
    fig.savefig(out("migration_surface_feems.png"), bbox_inches="tight")
    # ---- assign nonbreeding birds -------------------------------------------------
    demes = deme_table(sp)
    demes.to_csv(out("demes.csv"), index=False)
    if len(nb):
        P = predict(Xn, sp, args.method,T=args.temp)
        summ = summarize(P, demes, nb.sample_id.values,
                         nb.get("lon", pd.Series([np.nan] * len(nb))).values,
                         nb.get("lat", pd.Series([np.nan] * len(nb))).values)
        summ.to_csv(out("nonbreeding_assignment_summary.csv"), index=False)
        long = pd.DataFrame(P, columns=demes.deme.values)
        node_pos = sp.node_pos            # (n_nodes, 2) lon/lat, FEEMS's order
        print("[check2] max |demes - sp.node_pos|:",
              np.abs(demes[["deme_lon", "deme_lat"]].to_numpy(float) - sp.node_pos).max())  # should be ~0
        long.insert(0, "sample_id", nb.sample_id.values)
        long = long.melt(id_vars="sample_id", var_name="deme", value_name="post_prob")
        long = long.merge(demes, on="deme")
        long.to_csv(out("nonbreeding_assignment_long.csv"), index=False)
        print(f"[assign] {len(nb)} birds; median MAP prob={summ.map_prob.median():.3f}; "
              f"median 95% set={int(summ.n_demes_95cs.median())} demes; "
              f"{(~summ.map_deme_sampled).sum()} MAP to unsampled demes")

    # ---- self-assignment test of full fit -------------------------------------------
    Pfull = predict(Xb, sp, args.method)
    dtf = deme_table(sp)
    top = Pfull.argmax(axis=1)
    own = sp.assned_node_idx
    own_col = np.array([np.where(dtf.deme.values == n)[0][0] for n in own])
    print("[self] top deme == own deme:", (top == own_col).mean())
    print("[self] median km to own deme:",
          np.median(haversine_km(coord_b[:, 0], coord_b[:, 1],
                                 dtf.deme_lon.values[top], dtf.deme_lat.values[top])))

    # ---- sink diagnostics ---------------------------------------------------
    deg = np.array([sp.degree(n) for n in range(len(sp))])
    q_full = np.full(len(sp), np.nan)
    q_full[np.asarray(sp.perm_idx)[:sp.n_observed_nodes]] = np.asarray(sp.q)   # sampled nodes only
    pos = sp.node_pos
    obs = ~np.isnan(q_full)
    print(f"[sink] degree median {np.median(deg):.0f} (range {deg.min()}-{deg.max()}); "
          f"q median {np.nanmedian(q_full):.3g} (range {np.nanmin(q_full):.3g}-{np.nanmax(q_full):.3g})")
    for lon, lat in [(-99.025, 19.597), (-104.882, 51.728), (-101.009, 20.871)]:
        n = int(np.argmin(np.hypot(pos[:, 0] - lon, pos[:, 1] - lat)))
        q_rank = (np.nansum(q_full > q_full[n]) + 1) if obs[n] else None
        print(f"[sink] node {n} ({lon},{lat}): sampled={obs[n]}, degree {deg[n]}, "
              f"q {q_full[n]:.3g}" + (f" (rank {q_rank} of {obs.sum()} by q)" if q_rank else ""))


    # ---- which demes look like the global average? --------------------------
    node_of = np.asarray(sp.assned_node_idx)
    g = np.nanmean(Xb, axis=0)                                   # global mean genotype per SNP
    rows = []
    for n in np.unique(node_of):
        Xn_ = Xb[node_of == n]
        f = np.nanmean(Xn_, axis=0)
        het = np.nanmean((Xn_ == 1), axis=0).mean()
        rows.append((n, len(Xn_), np.sqrt(np.nanmean((f - g) ** 2)), het, *sp.node_pos[n]))
    dd = pd.DataFrame(rows, columns=["node", "n", "rms_dist_to_global", "obs_het", "lon", "lat"])
    dd["rank_closest"] = dd.rms_dist_to_global.rank().astype(int)
    print(dd.sort_values("rms_dist_to_global").head(8).to_string(index=False))
    print(dd[dd.node.isin([229, 17, 222])].to_string(index=False))

    # ---- leave-one-deme-out validation ---------------------------------------------
    if args.loo:
        obs = [(n, list(sp.nodes[n]["sample_idx"])) for n in range(len(sp))
               if len(sp.nodes[n]["sample_idx"]) > 0]
        rows = []
        LL, TC = [], []
        for k, (node, samples) in enumerate(obs, 1):
            sp.factor = None  # deepcopy can't copy the cholmod factor
            is_train = ~np.isin(np.arange(Xb.shape[0]), samples)
            tr, _ = train_test_split(sp, is_train)
            tr.fit_null_model(verbose=False)
            tr.fit(lamb=lamb, lamb_q=lamb_q, verbose=False)
           # Ph = predict(Xb[samples], tr, args.method)
            Zh = predict_ll(Xb[samples], tr, args.method)
            Ph = softmax(Zh / args.temp, axis=1)
            dt = deme_table(tr)
            true_col = int(np.where(dt.deme.values == node)[0][0])
            LL.append(Zh)
            TC.extend([true_col] * len(samples))
            if k == 1:
                # held-out deme really is unsampled in the training graph
                print("[check] held-out n_samples in tr:", len(tr.nodes[node]["sample_idx"]))
                # true_col points at the same location as the original node
                print("[check] sp node pos:", sp.node_pos[node],
                      " dt row:", dt.deme_lon.iloc[true_col], dt.deme_lat.iloc[true_col])
                # dt order matches graph order overall
                print("[check3] max |dt - tr.node_pos|:",
                      np.abs(dt[["deme_lon", "deme_lat"]].to_numpy(float) - tr.node_pos).max())
            for j, s in enumerate(samples):
                p = Ph[j]
                kmap = int(np.argmax(p))
                rows.append({
                    "sample_id": br.sample_id.iloc[s], "true_deme": node,
                    "true_lon": coord_b[s, 0], "true_lat": coord_b[s, 1],
                    "map_lon": dt.deme_lon.iloc[kmap], "map_lat": dt.deme_lat.iloc[kmap],
                    "err_km": haversine_km(coord_b[s, 0], coord_b[s, 1],
                                           dt.deme_lon.iloc[kmap], dt.deme_lat.iloc[kmap]),
                    "prob_true_deme": p[true_col],
                    "rank_true_deme": int((p > p[true_col]).sum() + 1),
                    # in the 95% set if mass of strictly-better demes is < 0.95
                    "true_in_95cs": bool(p[p > p[true_col]].sum() < 0.95),
                })
            print(f"[loo] {k}/{len(obs)} demes done", flush=True)
        loo = pd.DataFrame(rows)
        loo.to_csv(out("loo_validation.csv"), index=False)
        print(f"[loo] median error {loo.err_km.median():.0f} km; "
              f"true deme in 95% set {100 * loo.true_in_95cs.mean():.0f}%")
        LL = np.vstack(LL); TC = np.asarray(TC)
        np.save(out("loo_loglik.npy"), LL); np.save(out("loo_truecol.npy"), TC)
        Ts = np.geomspace(1, 1e5, 26)
        cov = []
        for T in Ts:
            P_ = softmax(LL / T, axis=1)
            pt = P_[np.arange(len(TC)), TC]
            cov.append(np.mean([p[p > q].sum() < 0.95 for p, q in zip(P_, pt)]))
        cov = np.asarray(cov)
        pd.DataFrame({"T": Ts, "coverage": cov}).to_csv(
            out("loo_temperature_calibration.csv"), index=False)
        ok = cov >= 0.95
        T_cal = float(Ts[ok.argmax()]) if ok.any() else float("nan")
        print(f"[calib] smallest T with >=95% LOO coverage: {T_cal:.1f}")

    # ---- leave-individual-out validation (stratified K-fold) ------------------
    if args.loo_ind:
        rng_ = np.random.default_rng(args.seed)
        K = min(args.ind_folds, Xb.shape[0])
        node_of = np.full(Xb.shape[0], -1)
        fold = np.full(Xb.shape[0], -1)
        for n in range(len(sp)):
            s = np.asarray(list(sp.nodes[n]["sample_idx"]), dtype=int)
            if len(s) == 0:
                continue
            node_of[s] = n
            off = rng_.integers(K)
            fold[rng_.permutation(s)] = (off + np.arange(len(s))) % K   # spread deme across folds
        assert (fold >= 0).all() and (node_of >= 0).all()

        rows, LL, TC, CO = [], [], [], []
        for k in range(K):
            test = np.where(fold == k)[0]
            if len(test) == 0:
                continue
            sp.factor = None
            tr, _ = train_test_split(sp, fold != k)
            tr.fit_null_model(verbose=False)
            tr.fit(lamb=lamb, lamb_q=lamb_q, verbose=False)
            Zh = predict_ll(Xb[test], tr, args.method)
            Ph = softmax(Zh / args.temp, axis=1)
            dt = deme_table(tr)
            for j, s in enumerate(test):
                node = node_of[s]
                tcol = int(np.where(dt.deme.values == node)[0][0])
                p = Ph[j]; kmap = int(np.argmax(p))
                LL.append(Zh[j]); TC.append(tcol); CO.append(dt.deme.values)
                rows.append({
                    "sample_id": br.sample_id.iloc[s], "fold": k, "true_deme": node,
                    "deme_sampled_in_train": len(tr.nodes[node]["sample_idx"]) > 0,
                    "true_lon": coord_b[s, 0], "true_lat": coord_b[s, 1],
                    "map_lon": dt.deme_lon.iloc[kmap], "map_lat": dt.deme_lat.iloc[kmap],
                    "err_km": haversine_km(coord_b[s, 0], coord_b[s, 1],
                                           dt.deme_lon.iloc[kmap], dt.deme_lat.iloc[kmap]),
                    "prob_true_deme": p[tcol],
                    "rank_true_deme": int((p > p[tcol]).sum() + 1),
                    "true_in_95cs": bool(p[p > p[tcol]].sum() < 0.95),
                })
            print(f"[loo_ind] fold {k + 1}/{K} done ({len(test)} birds)", flush=True)

        li = pd.DataFrame(rows)
        li.to_csv(out("loo_ind_validation.csv"), index=False)
        LL = np.vstack(LL); TC = np.asarray(TC)
        np.save(out("loo_ind_loglik.npy"), LL); np.save(out("loo_ind_truecol.npy"), TC);np.save(out("loo_ind_colorder.npy"), np.vstack(CO))
        T_ind = calibrate_T(LL, TC, out("loo_ind_temperature_calibration.csv"))

        for lab, d in [("all", li),
                       ("deme still sampled", li[li.deme_sampled_in_train]),
                       ("deme unsampled (singleton)", li[~li.deme_sampled_in_train])]:
            if len(d):
                print(f"[loo_ind] {lab:28s} n={len(d):3d}  median error {d.err_km.median():5.0f} km  "
                      f"top-deme hit {100 * (d.rank_true_deme == 1).mean():4.1f}%  "
                      f"in 95% set {100 * d.true_in_95cs.mean():4.1f}%")
        print(f"[calib_ind] smallest T with >=95% coverage: {T_ind:.1f}")


    json.dump({"lamb": lamb, "lamb_q": lamb_q, "method": args.method,
               "n_breeding": len(br), "n_nonbreeding": len(nb),
               "n_snps": int(Xb.shape[1]), "n_demes": len(sp),
               "n_sampled_demes": int(sp.n_observed_nodes)},
              open(out("run_info.json"), "w"), indent=2)
    print(f"[done] outputs in {args.outdir}")


if __name__ == "__main__":
    main()

