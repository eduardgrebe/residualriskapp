# Credits

## Developers

The model, Python package and web application were designed and developed by Eduard Grebe (Vitalant Research Institute) (<egrebe@vitalant.org>).

Additional assistance by:

* Brian Custer, Vitalant Research Institute (<bcuster@vitalant.org>) – Conceptualization, supervision, oversight, guidance and financial support
* Vivian I. Avelino-Silva: Conceptualization
* Marjorie D. Bravo: Collaboration and data curation
* Michael P. Busch: Conceptualization and guidance
* Artur Belov, U.S. Food and Drug Administration (<artur.belov@fda.hhs.gov>) – Infectivity model development and `k` parameter posterior distributions – human and animal data

## References

Original framework for estimating risk-day equivalents:

> Weusten, JJAM, Van Drimmelen HAJ, Lelie NP. Mathematic Modeling of the Risk of HBV, HCV, and HIV Transmission by Window‐phase Donations Not Detected by NAT. *Transfusion.* 2002;42(5):537-548. doi:[10.1046/j.1537-2995.2002.00099.x](https://doi.org/10.1046/j.1537-2995.2002.00099.x).

> Weusten, J, Vermeulen M, Van Drimmelen, Lelie N. Refinement of a Viral Transmission Risk Model for Blood Donations in Seroconversion Window Phase Screened by Nucleic Acid Testing in Different Pool Sizes and Repeat Test Algorithms. *Transfusion.* 2011;51(1):203-215. doi:[10.1111/j.1537-2995.2010.02804.x](https://doi.org/10.1111/j.1537-2995.2010.02804.x).

Baseline residual risk model and PrEP risk model framework:

> Grebe E, Busch MP, Notari EP, et al. HIV incidence in US first-time blood donors and transfusion risk with a 12-month deferral for men who have sex with men. *Blood.* 2020;136(11):1359-1367. doi:[10.1182/blood.2020007003](https://doi.org/10.1182/blood.2020007003).

> Grebe E, Avelino-Silva VI, Bravo MD, Busch MP, Custer B. Development of a risk assessment model of HIV transfusion transmission associated with undisclosed use of pre-exposure prophylaxis (PrEP) by blood donors. [ISBT Abstract PA28-L04. Oral presentation; 35th Regional Congress of the ISBT, Milan, Italy.] Vox Sang. 2025;120(Suppl. 1):110.

Transmissibility dose-response model:

> Belov A, Yang H, Forshee RA, et al. Modeling the risk of HIV transfusion transmission. *J Acquir Immune Defic Syndr.* 2023;92(2):173-179. doi:[10.1097/QAI.0000000000003115](https://doi.org/10.1097/QAI.0000000000003115).

## Funding

Development of this tool was primarily sponsored by [Vitalant Research Institute](https://research.vitalant.org), with additional support from [Eduard Grebe Consulting](https://grebe.consulting).

## Usage terms and license

Source code and text copyright (c) 2025-2026 Vitalant, with components (c) Eduard Grebe Consulting. All code released under the AGPL. 

You are free to use, copy and host instances of the app, for noncommercial and commercial applications, as long as the creators are credited, and the terms of the AGPL are complied with.

[![License: AGPL v3](https://img.shields.io/badge/License-AGPL%20v3-blue.svg)](https://www.gnu.org/licenses/agpl-3.0)

## Citation

We request that you cite the creators if you use this software in analysis, reports or publications, with the following suggested citation:

> Grebe E. (2026). Residual HIV Transfusion Transmission Risk Estimator (Version 1.1.0) [Computer software]. Vitalant Research Institute. [https://codeberg.org/eduardgrebe/residualriskapp](https://codeberg.org/eduardgrebe/residualriskapp).

BibTeX Entry:

```
@software{grebe2026,
  author       = {Grebe, Eduard},
  title        = {Residual HIV Transfusion Transmission Risk Estimator},
  version      = {1.1.0},
  year         = {2026},
  publisher    = {Vitalant Research Institute},
  doi          = {},
  url          = {https://codeberg.org/eduardgrebe/residualriskapp}
}
```