# NAT assay parameters

*Eduard Grebe — Vitalant Research Institute*

The Estimator can load the published limits of detection (LoD) of seven
blood-screening nucleic-acid test (NAT) assays instead of requiring manual entry.
This page documents where those numbers come from. It is a **condensed** summary
of the companion analysis `residualrisk_analysis/assays/ASSAYS.qmd`, which holds
the full per-analyte tables, fiducial limits, probit-fitting code and complete
reference list.

Two things to keep in mind throughout:

- **Only HIV-1 Group M is relevant here.** This tool models HIV transfusion
  transmission, so the canned values are the HIV-1 Group M screening
  (multiplex) LoDs. `ASSAYS.qmd` additionally covers HIV-1 Group O, HIV-2, HCV
  and HBV for each assay.
- **The model works in copies/mL.** Manufacturers report LoDs in international
  units (IU/mL). Each value below is the manufacturer's IU/mL LoD multiplied by
  an assay-specific IU→copies conversion factor (see *IU/mL → copies/mL
  conversion*). The conversion factor and WHO standard are surfaced in the app
  for transparency; they are **provenance only** and do not themselves enter the
  calculation — the operative inputs are the copies/mL 50% LoD, its standard
  deviation, and the 95%:50% LoD ratio.

## Canned assay limits of detection

HIV-1 Group M, screening (multiplex) assay, in **copies/mL**. `RSE` is the
relative standard error of the 50% LoD (`SD ÷ 50% LoD`); for every assay except
Bio-Manguinhos it derives from a 95% CI of the 50% LoD (the coefficient of
variation `CoV = (CI_hi − CI_lo) / 3.92 / PE`, computed in IU/mL and invariant
under the conversion). That CI is the manufacturer's where the insert publishes a
50% LoD; for the cobas TaqScreen MPX and MPX v2.0 — whose inserts give only the 95%
LoD — both the 50% LoD and its CI come from our probit fit of the insert's
per-concentration reactivity data (see *Probit-fitted 50% LoDs*). `cp/IU` is the
conversion factor applied.

| Assay | 50% LoD | 95% LoD | 50% LoD SD | RSE | WHO IS (HIV-1) | cp/IU |
|---|---:|---:|---:|---:|---|---:|
| Procleix Ultrio (Tigris) | 5.0 | 12.2 | 0.364 | 7.28% | 1st IS 97/656 | 0.60 |
| Procleix Ultrio Plus (Tigris) | 2.7 | 12.3 | 0.191 | 7.07% | 2nd IS 97/650 | 0.58 |
| Procleix Ultrio Elite (Panther) | 3.1 | 10.4 | 0.234 | 7.55% | 2nd IS 97/650 | 0.58 |
| cobas TaqScreen MPX (s 201) | 5.5 | 29.4 | 0.385 | 7.00% | 1st IS 97/656 | 0.60 |
| cobas TaqScreen MPX v2.0 (s 201) | 5.3 | 26.8 | 0.250 | 4.72% | 2nd IS 97/650 | 0.58 |
| cobas MPX (5800/6800/8800) | 1.3 | 9.0 | 0.0785 | 6.04% | 3rd IS 10/152 | 0.35 |
| Brazilian NAT Platform (Bio-Manguinhos) | 27.13 | 55.6 | 3.527 | 13.00% † | 2nd IS 97/650 | 0.58 |

The 95%:50% LoD ratio used by the simulation is simply `95% LoD ÷ 50% LoD` of the
values above (e.g. 2.44 for Ultrio, 6.92 for cobas MPX); it is dimensionless and
unaffected by the conversion. **†** Bio-Manguinhos RSE is *assumed* — see
*Bio-Manguinhos: a provisional standard deviation*.

## Sources for the limits of detection

Every value is traceable to a specific table, page and revision of the
manufacturer's package insert (or, for Bio-Manguinhos, a peer-reviewed paper).
The IU/mL → copies/mL arithmetic is shown so each copies/mL figure can be
reproduced.

- **Procleix Ultrio (Tigris)** — Grifols Diagnostic Solutions, IFU 502623
  Rev. 007 (May 2020), Tables 36a/36b, p.82. HIV-1 50%/95% = **8.4 / 20.3 IU/mL**
  (discriminatory dHIV-1 assay vs 1st IS 97/656 — the only HIV-1 IU value, since
  the multiplex screen is calibrated in copies/mL against a subtype-B virus)
  × 0.60 → 5.0 / 12.2.
- **Procleix Ultrio Plus (Tigris)** — Grifols, IFU AW-12765 Rev. 006 (June 2020),
  Table 16, p.43. HIV-1 50%/95% = **4.7 / 21.2 IU/mL** × 0.58 → 2.7 / 12.3.
- **Procleix Ultrio Elite (Panther)** — Grifols, Master IFU v7.0 (15 May 2026),
  Table 21, p.42. HIV-1 50%/95% = **5.4 / 18.0 IU/mL** × 0.58 → 3.1 / 10.4.
- **cobas TaqScreen MPX (s 201)** — Roche Molecular Systems, P/N 04584252 190
  (© 2009), Table 8, p.28. 95% LoD = **49 IU/mL**; the insert publishes no 50%
  LoD, so the **50% (9.1 IU/mL) is probit-fitted** from the insert's
  per-concentration reactivity data (see *Probit-fitted 50% LoDs*). × 0.60 →
  5.5 / 29.4.
- **cobas TaqScreen MPX v2.0 (s 201)** — Roche, IFU 06457258001-06EN, Doc Rev. 5.0
  (Oct 2021), Table 1, p.21. 95% LoD = **46.2 IU/mL** (vs 2nd IS 97/650); **50%
  (9.2 IU/mL) probit-fitted**. × 0.58 → 5.3 / 26.8.
- **cobas MPX (5800/6800/8800)** — Roche, IFU 07237278001, Doc Rev. 6.0 (Jan 2020),
  Table 11, p.24, for the 95% LoD = **25.7 IU/mL** (EDTA plasma, vs 3rd IS
  10/152). The 50% LoD = **3.8 IU/mL** is taken from the **European CE/IVD v3.0
  insert** (IFU 09199659001-03EN, Doc Rev. 3.0, ~Mar 2025, Table 16, p.31), which
  the U.S. copy omits. × 0.35 → 1.3 / 9.0.
- **Brazilian NAT Platform (Bio-Manguinhos)** — Rocha D, Andrade E, Godoy DT,
  et al. *The Brazilian experience of nucleic acid testing…* **Transfusion**
  2018;58(4), doi:[10.1111/trf.14478](https://doi.org/10.1111/trf.14478),
  Table 2 (Probit; 24 replicates/dilution for HIV), p.865. HIV-1 50%/95% =
  **46.77 / 95.86 IU/mL** (vs 2nd IS 97/650) × 0.58 → 27.13 / 55.60.

## WHO International Standards

IU/mL values are only strictly comparable when calibrated against the **same** WHO
International Standard (IS). The HIV-1 RNA IS has gone through three generations,
and the assays here span all three — so the IU figures (and the IU→copies factor)
are **not interchangeable** across assays:

| HIV-1 RNA IS | NIBSC code | Used by |
|---|---|---|
| 1st IS | 97/656 (genotype B; cross-reacts with HBV DNA) | Procleix Ultrio (dHIV-1), cobas TaqScreen MPX |
| 2nd IS | 97/650 (Group M, subtype B) | Ultrio Plus, Ultrio Elite, cobas TaqScreen MPX v2.0, Bio-Manguinhos |
| 3rd IS | 10/152 | cobas MPX |

NIBSC codes are taken verbatim from the inserts; generations are confirmed against
the NIBSC catalogue / WHO ECBS reports. A re-issued IS is re-calibrated, so the
apparent sensitivity differences between assays partly reflect standard changes
rather than analytical performance alone.

## IU/mL → copies/mL conversion

**No WHO International Standard assigns a genome-copies-per-IU value.** Every
standard in scope is value-assigned in **IU only** — deliberately, because
copy-number measurements vary by orders of magnitude across NAT methods (WHO's
stated rationale: copy number "is not a robust measure that can be compared
readily between laboratories"). Every factor is therefore a **method-anchored
convention, not a physical constant.**

For HIV-1 Group M the factor also drifts across IS generations:

| HIV-1 IS | cp/IU | Basis |
|---|---:|---|
| 1st IS 97/656 | 0.60 | Manufacturer/assay convention (Roche cobas, Grifols inserts) |
| 2nd IS 97/650 | 0.58 | IS-specific cross-calibration — Lelie & van Drimmelen, *J Med Virol* 2020 (PMID 32285945): 0.58 (0.51–0.66) |
| 3rd IS 10/152 | 0.35 | Hologic Aptima HIV-1 Quant Dx insert (12 cp/mL = 35 IU/mL); Lelie & van Drimmelen 2020 |

**Direction.** All factors are **copies per 1 IU** (copies/mL = IU/mL × factor).
The equivalent "1 copy ≈ 1.7 IU" is the inverse relationship and must **not** be
used as a multiplier (it would inflate copy-based LoDs ~2.9-fold).

**Internal cross-check (3rd IS).** cobas MPX reports HIV-1 Group O directly in
copies/mL (95% LoD 8.2 cp/mL) and Group M in IU/mL (95% LoD 25.7). Since both are
the same assay, the implied factor is 8.2 / 25.7 ≈ **0.32 cp/IU** — corroborating
**0.35** for the 3rd IS and contradicting 0.6 (which would give a ~2× mismatch).

## Probit-fitted 50% LoDs

The three **Grifols Procleix** assays publish both a 50% and a 95% LoD directly.
Among the Roche assays:

- **cobas MPX** — the U.S. insert gives only the 95% LoD; its **50% LoD comes from
  the European CE/IVD v3.0 insert** (the same assay, identical 95% values).
- **cobas TaqScreen MPX** and **MPX v2.0** — the inserts publish only the 95% LoD,
  but tabulate the full per-concentration reactivity data. Their **50% LoDs are
  fitted here** by probit regression (a binomial GLM with a probit link in
  `log10(concentration)`, dose-estimated with `MASS::dose.p`). As validation, the
  same fit reproduces the manufacturer's reported 95% LoD; the fitted 50% is then
  carried into the table above.

## Bio-Manguinhos: a provisional standard deviation

Rocha et al. (2018) report only **point** 50%/95% LoDs (Probit on just 24
replicates/dilution for HIV), with no confidence interval, fiducial limits or
per-dilution hit-rate table. Unlike every other assay — whose 50% LoD SD is
derived from a manufacturer 95% CI — its SD cannot be derived and must be
**assumed**.

The value used is an assumed relative SD (RSE) of **13%** (6.08 IU/mL ≈ 3.527
cp/mL). A 24-replicate study warrants a larger relative SD than the well-powered
assays (RSE ~5–8%); anchoring on Procleix Ultrio (a similarly steep curve, ~120
replicates/dilution, RSE ~7.3%) and scaling RSE ∝ √(1/N) gives ~13%. An earlier
ballpark of **4.95 IU/mL (RSE 10.6%)** was used in prior analyses. This should be
revisited if the authors' per-dilution hit-rate table becomes available (it could
then be probit-refit for a proper delta-method standard error).

---

**Authoritative source.** This page is a condensed summary. The full reference —
with every per-virus and per-analyte LoD table, 95% fiducial/confidence limits,
the IU→copies conversion research (per virus × WHO IS, source-graded), the probit
code, and the complete citation list — is `residualrisk_analysis/assays/ASSAYS.qmd`
in the companion analysis repository.
