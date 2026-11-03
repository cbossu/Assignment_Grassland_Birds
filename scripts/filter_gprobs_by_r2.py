#!/usr/bin/env python3
"""
Post-process Beagle 3.3.2 imputation output from an ANGSD beagle-format run.

Beagle writes, for `out=<prefix> like=<infile>`:
    <prefix>.<infile>.log        run summary
    <prefix>.<infile>.gprobs.gz  imputed genotype PROBABILITIES, same layout as
                                 the ANGSD input: marker allele1 allele2, then
                                 P(AA) P(Aa) P(aa) per individual
    <prefix>.<infile>.r2         marker  allelicR2   (imputation quality, 0-1,
                                 NaN when uninformative)
    <prefix>.<infile>.phased.gz  best-guess phased haplotypes (not used here)

This script does three independent things you can mix and match:

1. Filter markers by allelic R2 (drop poorly-imputed sites before downstream
   population-genetic analyses). The filtered .gprobs.gz stays in the exact
   ANGSD/Beagle probability format, so it can be fed straight into tools that
   accept a -beagle genotype-likelihood/probability file (e.g. PCAngsd,
   NGSadmix, ngsLD) — no VCF conversion needed for those.

2. Call hard genotypes (argmax of the 3 probabilities per individual, with an
   optional minimum-posterior-probability cutoff below which a call is set
   to missing) and write them as a plain tab-delimited genotype table.

3. Write a VCF (VCFv4.2) of the hard-called genotypes, e.g. for Fst outlier
   analysis (vcftools --weir-fst-pop, PLINK, OutFLANK, pcadapt, BayeScan via
   a converter, etc).

   IMPORTANT CAVEATS for the VCF:
   - ANGSD's "allele1"/"allele2" (from -doMajorMinor) are the *major* and
     *minor* allele inferred from the data, not necessarily the reference
     genome's base. They're written to VCF REF/ALT as-is. This is fine for
     Fst — it's computed from allele frequencies and is symmetric in which
     allele is called REF vs ALT — but don't treat REF as "the reference
     genome allele" for anything that cares (e.g. annotation lookups).
   - There's no genotype quality/depth info from Beagle, only the posterior
     genotype probability, which is written as FORMAT/GP and used (via
     --min-call-prob) to decide GT and to fill FORMAT/DS (dosage).
   - Marker IDs are parsed as "<chrom>_<pos>" by splitting on the LAST
     underscore, which matches how ANGSD names markers. Override with
     --chrom-split if your naming differs.
   - Records are written in the order they appear in the .gprobs.gz file.
     ANGSD writes markers in chromosomal order, but if you concatenated
     multiple per-chromosome runs out of order, sort afterward, e.g.:
         bcftools sort out.vcf -Oz -o out.sorted.vcf.gz && tabix -p vcf out.sorted.vcf.gz
   - Output is plain (unindexed) text .vcf on purpose, not .vcf.gz, because
     plain gzip is NOT the same as bgzip and tools like tabix/bcftools will
     reject a plain-gzipped VCF. To compress + index for bcftools/tabix:
         bgzip out.vcf && tabix -p vcf out.vcf.gz
     vcftools and PLINK are both fine reading the plain .vcf directly.

Usage:
    python filter_gprobs_by_r2.py \\
        --gprobs impute.angsd_beagle.test.combined.subsetted.beagle.gz.gprobs.gz \\
        --r2 impute.angsd_beagle.test.combined.subsetted.beagle.gz.r2 \\
        --min-r2 0.3 \\
        --out-gprobs filtered.gprobs.gz \\
        --hard-calls filtered.geno.tsv.gz \\
        --vcf filtered.vcf \\
        --min-call-prob 0.9 \\
        --sample-names samples.txt          # optional, one ID per column, in order,
                                             # same order as your ANGSD -bam list
"""
import argparse
import gzip
import sys
from datetime import date


def open_maybe_gz(path, mode="rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def open_maybe_gz_write(path):
    if path.endswith(".gz"):
        return gzip.open(path, "wt")
    return open(path, "w")


def load_r2(r2_path, min_r2):
    """Return the set of marker IDs passing the R2 threshold.

    NaN R2 (uninformative marker) is dropped by default, same as missing.
    """
    keep = set()
    n_total = 0
    n_nan = 0
    with open_maybe_gz(r2_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            marker, r2 = parts[0], parts[1]
            n_total += 1
            if r2.lower() == "nan":
                n_nan += 1
                continue
            try:
                if float(r2) >= min_r2:
                    keep.add((marker, r2))
            except ValueError:
                continue
    sys.stderr.write(
        f"[r2] {n_total} markers in .r2 file, {n_nan} NaN, "
        f"{len(keep)} pass R2 >= {min_r2}\n"
    )
    # store r2 value alongside marker so the VCF can annotate INFO/R2
    return {m: r for m, r in keep}


def call_genotype(pAA, pAa, paa):
    probs = (pAA, pAa, paa)
    best = max(range(3), key=lambda i: probs[i])
    return best, probs[best]  # 0=hom allele1, 1=het, 2=hom allele2


def split_marker(marker, split_char):
    idx = marker.rfind(split_char)
    if idx == -1:
        return None, None
    return marker[:idx], marker[idx + 1:]


def collect_chroms(gprobs_path, keep, split_char):
    """First pass: collect distinct chromosome names (in order of first
    appearance) among markers that will pass the R2 filter, so the VCF
    header can declare ##contig lines and avoid bcftools/tabix warnings."""
    chroms = []
    seen = set()
    with open_maybe_gz(gprobs_path) as fh:
        for line in fh:
            fields = line.split()
            if not fields or fields[0] == "marker":
                continue
            marker = fields[0]
            if marker not in keep:
                continue
            chrom, _pos = split_marker(marker, split_char)
            if chrom and chrom not in seen:
                seen.add(chrom)
                chroms.append(chrom)
    return chroms


def vcf_header(sample_names, source_gprobs, chroms):
    lines = [
        "##fileformat=VCFv4.2",
        f"##fileDate={date.today().strftime('%Y%m%d')}",
        "##source=filter_gprobs_by_r2.py (Beagle 3.3.2 gprobs -> VCF)",
        f"##sourceFile={source_gprobs}",
    ]
    lines += [f"##contig=<ID={c}>" for c in chroms]
    lines += [
        '##INFO=<ID=R2,Number=1,Type=Float,Description="Beagle 3.3.2 allelic R2 (imputation quality)">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype (argmax of Beagle posterior genotype probabilities)">',
        '##FORMAT=<ID=GP,Number=3,Type=Float,Description="Estimated genotype probabilities (hom-allele1,het,hom-allele2) from Beagle 3.3.2">',
        '##FORMAT=<ID=DS,Number=1,Type=Float,Description="Estimated allele2 dosage: 0*P(hom-allele1)+1*P(het)+2*P(hom-allele2)">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(sample_names),
    ]
    return "\n".join(lines) + "\n"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gprobs", required=True, help="Beagle .gprobs.gz output")
    ap.add_argument("--r2", required=True, help="Beagle .r2 output")
    ap.add_argument("--min-r2", type=float, default=0.3,
                     help="minimum allelic R2 to keep a marker (default 0.3)")
    ap.add_argument("--out-gprobs", default=None,
                     help="write the R2-filtered .gprobs.gz here (probability format, "
                          "ready for PCAngsd/NGSadmix/ngsLD -beagle input)")
    ap.add_argument("--hard-calls", default=None,
                     help="also write argmax hard-called genotypes (0/0, 0/1, 1/1) as a "
                          "tab-delimited table")
    ap.add_argument("--vcf", default=None,
                     help="also write a VCFv4.2 file of hard-called genotypes (e.g. for "
                          "Fst outlier analysis). Plain text, not bgzipped — see caveats "
                          "in the module docstring for compressing/indexing afterward.")
    ap.add_argument("--min-call-prob", type=float, default=0.0,
                     help="hard calls (in --hard-calls and --vcf) below this posterior "
                          "probability are written as missing (./. in VCF), default 0.0 "
                          "= always call the argmax. For Fst work, 0.8-0.9 is typical.")
    ap.add_argument("--sample-names", default=None,
                     help="optional file, one sample ID per line in column order, "
                          "same order as your ANGSD -bam list. Used as the header for "
                          "--hard-calls and as VCF sample columns; defaults to ind0, ind1, ...")
    ap.add_argument("--chrom-split", default="_",
                     help="character used to split marker IDs into chrom/pos for the VCF "
                          "(splits on the LAST occurrence). Default '_' matches ANGSD's "
                          "default marker naming.")
    args = ap.parse_args()

    r2_by_marker = load_r2(args.r2, args.min_r2)
    keep = set(r2_by_marker)

    sample_names = None
    if args.sample_names:
        with open(args.sample_names) as fh:
            sample_names = [l.strip() for l in fh if l.strip()]

    out_gprobs = open_maybe_gz_write(args.out_gprobs) if args.out_gprobs else None
    hard_calls = open_maybe_gz_write(args.hard_calls) if args.hard_calls else None
    vcf_out = open_maybe_gz_write(args.vcf) if args.vcf else None

    vcf_chroms = collect_chroms(args.gprobs, keep, args.chrom_split) if vcf_out else None

    n_in = 0
    n_kept = 0
    header_written = False
    n_marker_parse_fail = 0
    vcf_header_written = False

    with open_maybe_gz(args.gprobs) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split()

            # Beagle's own header line (starts with "marker") - detect and handle once
            if fields[0] == "marker":
                n_ind = (len(fields) - 3) // 3
                if not sample_names:
                    sample_names = [f"ind{i}" for i in range(n_ind)]
                elif len(sample_names) != n_ind:
                    sys.exit(f"ERROR: --sample-names has {len(sample_names)} names but "
                             f"gprobs file has {n_ind} individuals")

                if out_gprobs:
                    out_gprobs.write(line + "\n")
                if hard_calls:
                    hard_calls.write("marker\tallele1\tallele2\t" + "\t".join(sample_names) + "\n")
                if vcf_out:
                    vcf_out.write(vcf_header(sample_names, args.gprobs, vcf_chroms))
                    vcf_header_written = True
                header_written = True
                continue

            n_in += 1
            marker = fields[0]
            if marker not in keep:
                continue
            n_kept += 1

            if out_gprobs:
                out_gprobs.write(line + "\n")

            if hard_calls or vcf_out:
                a1, a2 = fields[1], fields[2]
                probs = [float(x) for x in fields[3:]]
                calls = []
                gp_dose = []
                for i in range(0, len(probs), 3):
                    pAA, pAa, paa = probs[i:i + 3]
                    idx, p = call_genotype(pAA, pAa, paa)
                    dose = 0 * pAA + 1 * pAa + 2 * paa
                    if p < args.min_call_prob:
                        calls.append("./.")
                    else:
                        calls.append(["0/0", "0/1", "1/1"][idx])
                    gp_dose.append((pAA, pAa, paa, dose))

                if hard_calls:
                    hard_calls.write(f"{marker}\t{a1}\t{a2}\t" + "\t".join(calls) + "\n")

                if vcf_out:
                    chrom, pos = split_marker(marker, args.chrom_split)
                    if chrom is None:
                        n_marker_parse_fail += 1
                        continue
                    samp_fields = []
                    for gt, (pAA, pAa, paa, dose) in zip(calls, gp_dose):
                        samp_fields.append(f"{gt}:{pAA:.4g},{pAa:.4g},{paa:.4g}:{dose:.3f}")
                    r2_val = r2_by_marker.get(marker, "NaN")
                    vcf_out.write(
                        f"{chrom}\t{pos}\t{marker}\t{a1}\t{a2}\t.\tPASS\tR2={r2_val}\t"
                        f"GT:GP:DS\t" + "\t".join(samp_fields) + "\n"
                    )

    if out_gprobs:
        out_gprobs.close()
    if hard_calls:
        hard_calls.close()
    if vcf_out:
        vcf_out.close()

    sys.stderr.write(f"[gprobs] {n_in} markers read, {n_kept} kept after R2 filter"
                      f"{' (no header line found in gprobs file!)' if not header_written else ''}\n")
    if args.vcf:
        if not vcf_header_written:
            sys.stderr.write("[vcf] WARNING: no VCF written - gprobs file had no header line\n")
        elif n_marker_parse_fail:
            sys.stderr.write(f"[vcf] WARNING: {n_marker_parse_fail} markers could not be split "
                              f"into chrom/pos with --chrom-split='{args.chrom_split}' and were "
                              f"skipped from the VCF\n")


if __name__ == "__main__":
    main()
