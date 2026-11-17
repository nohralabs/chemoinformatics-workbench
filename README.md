# Chemoinformatics Workbench
[![Tests](https://github.com/nohralabs/chemoinformatics-workbench/actions/workflows/tests.yml/badge.svg)](https://github.com/nohralabs/chemoinformatics-workbench/actions/workflows/tests.yml)

A reproducible molecular-data platform for structure validation, chemical standardization, descriptor calculation, similarity searching, chemical-space exploration and QSAR modelling using Python and RDKit.

> **Project status:** Under active development. Features and case studies will be released progressively with tested code, documented data provenance and reproducible results.

## Scientific Objective

Chemical datasets frequently contain invalid structures, inconsistent molecular representations, salts, duplicates, missing stereochemistry and incompatible property records. These issues can distort similarity analyses, molecular descriptors and predictive models.

This project develops a transparent workflow that converts heterogeneous molecular records into standardized, traceable and analysis-ready chemical data.

## Planned Workflow

1. Import molecular records from SMILES and tabular files
2. Validate molecular structures and record parsing errors
3. Standardize structures, charges, fragments and tautomers
4. Generate canonical SMILES, molecular formulae and exact masses
5. Calculate physicochemical and structural descriptors
6. Generate molecular fingerprints
7. Perform similarity and substructure searches
8. Explore chemical space using dimensionality reduction and clustering
9. Store curated structures and calculated properties in SQL
10. Develop validated QSAR baselines with applicability-domain and uncertainty analysis

## Scientific Case Studies

### Fragrance Molecular Space

Analysis of public fragrance molecules using structural descriptors, similarity maps and molecular-property relationships.

### Cyclic Carbonates and PHU Building Blocks

Chemoinformatics analysis of cyclic carbonates, amines and hydroxyurethane-related structures connected to sustainable non-isocyanate polyurethane chemistry.

### Polymer and Adhesive Building Blocks

Comparison of monomers and functional building blocks relevant to polymers, coatings and adhesive technologies.

All datasets will be public, openly licensed or independently generated. No confidential employer or client information will be included.

## Technology Stack

- Python
- RDKit
- pandas and NumPy
- SciPy and scikit-learn
- Matplotlib and Seaborn
- SQLite and SQL
- Jupyter
- pytest
- GitHub Actions

## Repository Structure

```text
chemoinformatics-workbench/
├── data/
│   ├── raw/
│   ├── interim/
│   └── processed/
├── notebooks/
├── src/
│   └── chemoinformatics_workbench/
├── tests/
├── reports/
│   └── figures/
├── docs/
├── README.md
├── pyproject.toml
└── environment.yml
```

## Development Roadmap

- [ ] Establish the reproducible Python and RDKit environment
- [ ] Create the molecular-data schema
- [ ] Implement SMILES validation
- [ ] Implement molecular standardization
- [ ] Calculate descriptors and fingerprints
- [ ] Build similarity and substructure search
- [ ] Create the SQL molecular database
- [ ] Develop scientific case studies
- [ ] Add QSAR modelling and uncertainty analysis
- [ ] Add automated tests and continuous integration
- [ ] Publish a documented release

## Reproducibility and Quality

Every transformation will preserve the original molecular record, the standardized result, processing status and provenance. Core functions will be tested against chemically meaningful edge cases, including invalid SMILES, mixtures, charged species, stereoisomers and duplicate structures.

## Author

**Bassam Nohra, PhD**  
Chemistry, Materials Science, Chemoinformatics and Scientific AI
