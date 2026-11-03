#!/usr/bin/env python
"""cluster_assign.py — score FEEMS assignment at the genoscape-cluster level.

Sums each bird's tempered deme posterior within genoscape clusters, then
  * scores the leave-individual-out breeding birds (accuracy, coverage, confusion),
  * writes cluster posteriors for the nonbreeding birds.

Needs from feems_assign_nonbreeding.py --loo_ind:
  loo_ind_loglik.npy, loo_ind_truecol.npy, loo_ind_colorder.npy, loo_ind_validation.csv
and from deme_clusters.R:  deme_clusters.csv (deme, cluster)

Usage:
  python cluster_assign.py --run_dir feems_LOSH_outer --clusters deme_clusters.csv --T 200
"""
import argparse
import os

import numpy as np
import pandas as pd
from scipy.special import softmax


def cluster_probs(P, col_demes, lab, clusters):
    """P: [n x d] deme probs; col_demes: [n x d] deme id of each column."""
    idx = {c: i for i, c in enumerate(clusters)}
    C = np.zeros((P.shape[0], len(clusters)))
    for i in range(P.shape[0]):
        ci = np.array([idx[lab.get(int(d), "none")] for d in col_demes[i]])
        np.add.at(C[i], ci, P[i])
    return C


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run_dir", required=True)
    ap.add_argument("--clusters", required=True, help="deme_clusters.csv from deme_clusters.R")
    ap.add_argument("--T", type=float, required=True, help="calibrated temperature")
    ap.add_argument("--equal_prior", action="store_true",
                    help="also score with equal prior weight per cluster "
                         "(deme probs divided by number of demes in their cluster)")
    a = ap.parse_args()
    out = lambda f: os.path.join(a.run_dir, f)

    dc = pd.read_csv(a.clusters)
    lab = dict(zip(dc.deme.astype(int), dc.cluster.astype(str)))
    clusters = sorted(set(lab.values()) | {"none"})
    n_dem = pd.Series(list(lab.values())).value_counts().reindex(clusters).fillna(0).values
    print("[demes per cluster] " + ", ".join(f"{c}={int(n)}" for c, n in zip(clusters, n_dem)))

    def reweight(C):
        """equal prior per named cluster: divide summed deme mass by n demes, renormalise.
        'none' demes (no genoscape colour) are dropped here; with few demes they would
        otherwise get an inflated share."""
        keep = (n_dem > 0) & (np.array(clusters) != "none")
        W = np.where(keep, C / np.maximum(n_dem, 1), 0.0)
        return W / W.sum(1, keepdims=True)

    versions = [("", "deme-uniform prior")]
    if a.equal_prior:
        versions.append(("_eqprior", "equal prior per cluster"))

    # ---- leave-individual-out breeding birds ----------------------------------
    LL = np.load(out("loo_ind_loglik.npy"))
    CO = np.load(out("loo_ind_colorder.npy"))
    val = pd.read_csv(out("loo_ind_validation.csv"))
    assert len(LL) == len(CO) == len(val), "LOO arrays and table differ in length"

    P = softmax(LL / a.T, axis=1)
    C0 = cluster_probs(P, CO, lab, clusters)
    true_c = val.true_deme.astype(int).map(lambda d: lab.get(d, "none")).values
    ti = np.array([clusters.index(c) for c in true_c])
    summary = []

    for suf, desc in versions:
        C = reweight(C0) if suf else C0
        map_i = C.argmax(1)
        p_true = C[np.arange(len(C)), ti]
        in95 = np.array([C[i][C[i] > p_true[i]].sum() < 0.95 for i in range(len(C))])
        set95 = np.array([((np.cumsum(np.sort(c)[::-1]) - np.sort(c)[::-1]) < 0.95).sum() for c in C])

        res = val[["sample_id", "true_deme", "deme_sampled_in_train"]].copy()
        res["true_cluster"] = true_c
        res["map_cluster"] = [clusters[i] for i in map_i]
        res["map_cluster_prob"] = C.max(1)
        res["prob_true_cluster"] = p_true
        res["true_in_95cs"] = in95
        res["n_clusters_95cs"] = set95
        for j, c in enumerate(clusters):
            res[f"p_{c}"] = C[:, j]
        res.to_csv(out(f"loo_ind_cluster_assignment{suf}.csv"), index=False)

        conf = pd.crosstab(res.true_cluster, res.map_cluster, rownames=["true"], colnames=["assigned"])
        conf.to_csv(out(f"loo_ind_cluster_confusion{suf}.csv"))

        acc = (res.true_cluster == res.map_cluster).mean()
        base = res.true_cluster.value_counts(normalize=True).max()
        bal = np.mean([(d.map_cluster == c).mean() for c, d in res.groupby("true_cluster") if c != "none"])
        lines = [f"Cluster-level scoring of leave-individual-out birds (T = {a.T:g}; {desc})",
                 f"  birds: {len(res)}",
                 f"  MAP cluster correct: {100 * acc:.1f}%"
                 f"   (baseline, always guess commonest cluster: {100 * base:.1f}%)",
                 f"  balanced accuracy (mean of per-cluster accuracy, named clusters): {100 * bal:.1f}%",
                 f"  mean prob on true cluster: {res.prob_true_cluster.mean():.2f}",
                 f"  true cluster in 95% set: {100 * res.true_in_95cs.mean():.1f}%",
                 f"  median clusters in 95% set: {int(res.n_clusters_95cs.median())} of {len(clusters) - 1}",
                 "", "  by true cluster:"]
        for c, d in res.groupby("true_cluster"):
            lines.append(f"    {c:12s} n={len(d):3d}  correct {100 * (d.map_cluster == c).mean():5.1f}%  "
                         f"mean p(true) {d.prob_true_cluster.mean():.2f}")
        lines += ["", "  confusion (rows = true, cols = assigned):", conf.to_string()]
        if (res.true_cluster == "none").any():
            lines.append("\n  note: some true demes have no genoscape colour ('none'); "
                         "increase buffer_km in deme_clusters.R")
        txt = "\n".join(lines)
        print(txt + "\n")
        open(out(f"loo_ind_cluster_metrics{suf}.txt"), "w").write(txt + "\n")
        summary.append((desc, acc, bal))

    if len(summary) > 1:
        print("[compare] " + " | ".join(f"{d}: accuracy {100 * x:.1f}%, balanced {100 * b:.1f}%"
                                        for d, x, b in summary) + "\n")

    # ---- nonbreeding birds -----------------------------------------------------
    f_long = out("nonbreeding_assignment_long.csv")
    if os.path.exists(f_long):
        lg = pd.read_csv(f_long)
        med_deme = lg.groupby("sample_id").post_prob.max().median()
        if med_deme > 0.5:
            print(f"[warn] median MAP deme prob in {f_long} is {med_deme:.2f}: "
                  "this run looks UNTEMPERED. Rerun it with --temp before using these.")
        lg["cluster"] = lg.deme.astype(int).map(lambda d: lab.get(d, "none"))
        nb0 = lg.pivot_table(index="sample_id", columns="cluster", values="post_prob",
                             aggfunc="sum", fill_value=0.0).reindex(columns=clusters, fill_value=0.0)
        for suf, desc in versions:
            pc = reweight(nb0.values) if suf else nb0.values / nb0.values.sum(1, keepdims=True)
            nb = pd.DataFrame(pc, index=nb0.index, columns=[f"p_{c}" for c in clusters])
            nb.insert(0, "map_cluster", [clusters[i] for i in pc.argmax(1)])
            nb.insert(1, "map_cluster_prob", pc.max(1))
            nb.reset_index().to_csv(out(f"nonbreeding_cluster_assignment{suf}.csv"), index=False)
            print(f"[nonbreeding, {desc}] {len(nb)} birds; median MAP deme prob {med_deme:.3f}; "
                  f"median MAP cluster prob {nb.map_cluster_prob.median():.2f}\n"
                  f"{nb.map_cluster.value_counts().to_string()}\n")


if __name__ == "__main__":
    main()
