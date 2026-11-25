# HCV Recombination Patterns from a Global Whole-Genome Sequence Analysis

Code, phylogenetic trees, reference protein structures, and supplementary tables accompanying the following study, submitted to *Virus Evolution* (VEVOLU-2026-103).

> Tshiabuila, D., James, S. E., Maponga, T. G., Delphin, M., Matthews, P. C., Wilkinson, E., Dor, G., Baxter, C., Martin, D. P., & de Oliveira, T. (2026). *HCV Recombination Patterns from a Global Whole-Genome Sequence Analysis*.

## Data availability

- [GitHub repository](https://github.com/CERI-KRISP/HCV-Recombination-Patterns)
- [Zenodo archive DOI](https://doi.org/10.5281/zenodo.24467237)
- [Archive contents and verification status](data_archive_note.md)

## Dataset and quality control

The starting inventory contained 5,938 genomes, including 143 curated references. The ORF screen retained 5,776 genomes. Manual alignment review removed one sequence before the 5,775-genome preliminary FastTree screen. Removal of nine outliers left 5,766 aligned genomes, including 5,623 queries and 143 references. The query set comprised 5,614 genotype-assigned genomes and nine unassigned genomes. The reference panel comprised 134 non-recombinant references and nine known recombinant controls.

## Repository contents

| Path | Contents |
|---|---|
| [scripts](scripts/) | ORF extraction, TranslatorX and MAFFT alignment, stop-codon QC, AlphaFold3 CIF-to-PDB conversion, and GARD coordinate mapping |
| [reference_proteins](reference_proteins/) | Ten H77 protein sequences recovered from the archived model-chain records and the ten corresponding AlphaFold3 PDB models used for the structural analysis |
| [trees](trees/) | Preliminary FastTree, final IQ-TREE 3 global tree, and minor, major and available component trees for the 54 supported events |
| [supplementary_tables](supplementary_tables/) | The five manuscript supplementary tables |

The ORF and alignment scripts are the archived study utilities. The ORF helper extracts the longest ATG-initiated ORF and flags records below the length threshold. The optional `--trim-incomplete-stop` applies the documented partial-stop rule. This option was added for the release and has not been benchmarked against the 907 recorded trims. Use the archived final alignment to reproduce downstream results, which also reflect manual QC edits. The alignment script expects a local TranslatorX installation supplied through its documented environment variable. The CIF conversion script reads the CIF files from the working directory. The GARD coordinate script reads the archived retained-column map and reports both master-alignment flanks for a one-based GARD boundary.

The [regional tree manifest](trees/regional_validation_trees/tree_manifest.csv) maps published event IDs to both primary trees. The [component manifest](trees/regional_validation_trees/component_tree_manifest.csv) lists available disjoint-segment trees. Trees are archived as inferred, and genotype 7 display rooting is applied in the figures.

The reference protein directory preserves the archived PDB filenames. The model named `ns1` in those files is catalogued as `p7` in the FASTA directory.

## Supplementary tables

| Table | Contents |
|---|---|
| S1 | Accession inventory and quality-control status for all 5,938 genomes |
| S2 | 129 RDP5 candidate events included in the study. E090 is excluded because the export lists two supporting methods. |
| S3 | RDP5 and regional phylogenetic comparison for 62 LANL-annotated recombinant accessions |
| S4 | 76 supported event–accession records across 54 events and 59 distinct genomes, including five reference controls. Includes alignment and H77 breakpoints, validation spans and lengths, and branch support. The unassigned screen uses its own source-alignment columns for regional spans. |
| S5 | Concordance between ten internal RDP5 endpoints and mapped GARD boundary flanks for seven representatives |

## Software

The manuscript reports MAFFT 7.526, TranslatorX 1.1, FastTree 2.2.0, IQ-TREE 3.1.3, RDP5.93, HyPhy 2.5.103 with GARD analysis 0.2, and AlphaFold3. The Python structure conversion uses Biopython.

## Ethics and funding

Stellenbosch University Health Research Ethics Committee approval N22/08/089. Funding is described in the manuscript.
