# Artificial Selection of Multispecies Microbial Inoculants

Code and data to reproduce the figures of:

> Benítez Domínguez B, *et al.* Artificial Selection of Multispecies Microbial Inoculants. *bioRxiv* (2026). DOI: [to be added]

## Abstract

<p align="justify">
Evolutionary engineering can optimize microbial consortia through iterative ecological perturbation and selection. It remains unclear whether selection improves community function simply by enriching communities with the strongest individual contributors or whether it can also enrich communities with synergistic interactions among species. Here, we address this question by optimizing the initial composition of multispecies wine yeast inoculants. Each generation we independently perturb the starting abundance of each of six yeast strains, toggling it between two of six possible doses (from 0 to 10<sup>7</sup> CFU/mL), and select for two functions of contrasting complexity. Unlike previous studies focused on serially passaged communities, our approach optimizes the initial inoculum composition for single-batch fermentation, a mode of operation common to many biotechnological processes. Selecting for a function that all individual strains can independently contribute to results in inoculants that are enriched with the strongest contributors. In contrast, selecting for a function that no strain can perform alone yields more complex communities whose performance emerges from synergistic interactions. These findings show that evolutionary engineering can exploit both individual contributions and interspecies interactions to optimize microbial communities. Our approach can be adapted to optimize other microbial inoculants of low-to-medium diversity.
</p>

## Repository structure

Each figure has its own folder with a self-contained R script, the data it needs and its output.

| Folder | Script | Figure |
|---|---|---|
| `Fig2/` | `Fig2_monocultures.R` | Fig 2. Two target functions at contrasting levels of functional complexity, and the monoculture performance that motivates selecting communities (bottom panels) |
| `Fig3/` | `Fig3_HS_FR.R` | Fig 3. Sugar consumption in the HS-FR must is largely additive and dominated by a single strain |
| `Fig4/` | `Fig4_lactic_acid.R` | Fig 4. Lactic acid production is interactive |
| `Fig5/` | `Fig5_epistasis.R` | Fig 5. Pairwise interactions are context-dependent |
| `S1_Fig/` | `S1_Fig_null.R` | S1 Fig. Contribution of pairwise interactions compared with experimental noise |
| `S2_Fig/` | `S2_Fig_Td10_SC2.R` | S2 Fig. Interaction between *T. delbrueckii* (Td10) and *S. cerevisiae* SC2 near the high-sugar optimum (third round) |
| `S1_Methods/` | `S1_Methods_FigA_WH.R` | Fig A in S1 Methods. Contribution of each interaction order compared with experimental noise |

```
FigX/
├── FigX_name.R     script
├── data/
│   ├── it1/        CSV files of round 1
│   ├── it2/        CSV files of round 2
│   └── it3/        CSV files of round 3
└── output/         figures (.svg, .pdf, .png), created when the script runs
```

## Data

One CSV file per function and round, stored in the folder of its round (`data/it1/`, `data/it2/`, `data/it3/`):

| File | Function | Round |
|---|---|---|
| `hsfr_it1_data.csv`, `hsfr_it2_data.csv`, `hsfr_it3_data.csv` | Sugar consumed in the HS-FR must (g/L) | 1, 2, 3 |
| `lactic_it1_data.csv`, `lactic_it2_data.csv`, `lactic_it3_data.csv` | L-lactic acid produced in SGM (g/L) | 1, 2, 3 |

Each row is one community position on the plate. The lactic acid files contain only the communities that finished the fermentation (≤ 5 g/L of residual sugar), so rounds 1 and 2 have 33 and 44 of the 64 combinations. Where a figure needs all 64 (Fig 5), a missing community counts as zero.

Columns:

| Column | Description |
|---|---|
| `Well_template` | Plate and position of the community (e.g. `plate2_H03`) |
| `Environment` | Code of the inoculum, one digit per strain in the order Mp21, Lt24, Td10, SC3, SC2, SC1 (0 = low, 1 = high dose) |
| `Mp21`, `Lt24`, `Td10`, `SC3`, `SC2`, `SC1` | Dose of each strain, 0 = low and 1 = high dose of that round |
| `*_R1`, `*_R2`, `*_R3` | Function measured in the three replicate wells |
| `Mean_*`, `SD_*` | Mean and standard deviation of the three replicate wells |

Doses (CFU/mL) of each strain in each round, as low / high. A dose of 0 means the strain was not inoculated. The doses of round 1 were the same for every strain; those of rounds 2 and 3 follow from the winning consortium of the previous round (Fig 1).

HS-FR (sugar consumption):

| Strain | it1 | it2 | it3 |
|---|---|---|---|
| Mp21 | 10³ / 10⁵ | 10⁵ / 10⁶ | 10⁴ / 10⁵ |
| Lt24 | 10³ / 10⁵ | 10⁵ / 10⁶ | 10⁴ / 10⁵ |
| Td10 | 10³ / 10⁵ | 10⁵ / 10⁶ | 10⁶ / 10⁷ |
| SC3  | 10³ / 10⁵ | 10⁵ / 10⁶ | 10⁶ / 10⁷ |
| SC2  | 10³ / 10⁵ | 0 / 10³ | 10³ / 10⁴ |
| SC1  | 10³ / 10⁵ | 0 / 10³ | 10³ / 10⁴ |

SGM (L-lactic acid production):

| Strain | it1 | it2 | it3 |
|---|---|---|---|
| Mp21 | 10³ / 10⁵ | 0 / 10³ | 10³ / 10⁴ |
| Lt24 | 10³ / 10⁵ | 10⁵ / 10⁶ | 10⁶ / 10⁷ |
| Td10 | 10³ / 10⁵ | 0 / 10³ | 0 / 0 |
| SC3  | 10³ / 10⁵ | 10⁵ / 10⁶ | 10⁶ / 10⁷ |
| SC2  | 10³ / 10⁵ | 0 / 10³ | 10³ / 10⁴ |
| SC1  | 10³ / 10⁵ | 0 / 10³ | 10³ / 10⁴ |

Fig 2 also uses three monoculture files in `data/monocultures/`, one row per strain and dose (CFU/mL):

| File | Measure |
|---|---|
| `hsfr_sugar_consumed.csv` | Sugar consumed in the HS-FR must (g/L) |
| `sgm_lactic_acid.csv` | L-lactic acid produced in SGM (g/L) |
| `sgm_residual_sugar.csv` | Residual sugar in SGM (g/L); a fermentation counts as finished with ≤ 5 g/L |

Columns: `strain`, `dose`, `R1`, `R2`, `R3` (the three replicate wells), `Mean` and `SD`.

Strains: *Metschnikowia pulcherrima* (Mp21), *Lachancea thermotolerans* (Lt24), *Torulaspora delbrueckii* (Td10) and *Saccharomyces cerevisiae* (SC3, SC2, SC1).

## Requirements

R (≥ 4.1) with the packages `ggplot2` (≥ 3.5.0), `patchwork`, `scales` and `svglite` (Fig 2 also uses `ggh4x`). Each script lists the packages it needs at the top.

To run a script, open it, set `fig_dir` to the folder of that figure (or run it from inside that folder) and source it. The figures are written to `output/` inside the same folder.
