# BiRAGAS CRISPR Complete v3.0 — DNA + RNA Analysis Platform

**Ayass Bioscience LLC** — The first platform to unify DNA editing and RNA targeting with causal inference for drug discovery.

## What's New in v3.0 (vs v2.0 DNA-only)

| Feature | v2.0 (DNA only) | v3.0 (DNA + RNA) |
|---------|-----------------|-------------------|
| **Cas13 RNA targeting** | No | **Cas13a, Cas13b, Cas13d (CasRx), dCas13** |
| **RNA base editing** | No | **A-to-I (ADAR2) + C-to-U (APOBEC)** |
| **CRISPRi/CRISPRa** | No | **Transcriptional repression/activation** |
| **Perturb-seq/CROP-seq** | No | **Single-cell CRISPR screen analysis** |
| **ncRNA targeting** | No | **lncRNA, miRNA, siRNA, circRNA** |
| **Spatial transcriptomics** | No | **CRISPR-TO support** |
| DNA knockout | 210,859 configs | **210,859 configs** (unchanged) |
| Combinations | 22.2B | **22.2B** (unchanged) |

## Quick Start

### Double-click the app
Open **`BiRAGAS_CRISPR_Complete_App.html`** in any browser. Works offline, no installation needed.

### Full backend
```bash
pip install numpy scipy networkx fastapi uvicorn pydantic
python -m BiRAGAS_CRISPR_Complete --serve
# Open http://localhost:8000
```

### CLI — DNA + RNA
```bash
# DNA guide design
python -m BiRAGAS_CRISPR_Complete --design BRCA1 --nuclease NGG

# RNA guide design (Cas13d)
python -m BiRAGAS_CRISPR_Complete --design BRCA1 --nuclease Cas13d

# DNA knockout strategy
python -m BiRAGAS_CRISPR_Complete --knockout TP53 --nuclease NGG

# RNA knockdown strategy
python -m BiRAGAS_CRISPR_Complete --knockout TP53 --nuclease Cas13d

# CRISPRi / CRISPRa
python -m BiRAGAS_CRISPR_Complete --crispri BRAF
python -m BiRAGAS_CRISPR_Complete --crispra BRAF

# Full pipeline (DNA + RNA)
python -m BiRAGAS_CRISPR_Complete --analyze /path/to/CRISPR --disease Melanoma
```

## Architecture

```
BiRAGAS_CRISPR_Complete/
├── core/                          # DNA Engines
│   ├── editing_engine.py          # Unified DNA+RNA guide design
│   ├── screening_engine.py        # MAGeCK + BAGEL2 + DrugZ
│   ├── knockout_engine.py         # 7-method ensemble (210K)
│   ├── mega_scale_engine.py       # Sparse matrix (22.2B)
│   ├── combination_engine.py      # 6-model synergy
│   └── ace_scoring_engine.py      # 15-stream ACE
├── rna/                           # RNA Engines (NEW)
│   ├── rna_base_edit_engine.py    # A-to-I + C-to-U editing
│   ├── transcriptome_engine.py    # CRISPRi/CRISPRa/Perturb-seq
│   └── noncoding_engine.py        # lncRNA/miRNA/siRNA
├── autonomous/                    # Self-correction
├── pipeline/unified_orchestrator.py  # DNA+RNA master pipeline
├── api/server.py                  # FastAPI (DNA+RNA endpoints)
└── BiRAGAS_CRISPR_Complete_App.html  # Standalone web app
```

## License
Copyright (c) Ayass Bioscience LLC. All Rights Reserved.
