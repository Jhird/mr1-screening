# mr1-screening

## In silico screening of cancer metabolites that may bind MR1 receptor

### Instructions

* Start with the `notebooks` folder. Each notebook is numbered (0-indexed) in the order that they should be run. Briefly:

    * **nb0** - Downloads all metabolic compounds found in the KEGG database and tries to add as many external identifiers as possible (PubChem, CHEMBL, CAS numbers). Uses [BioPython](https://biopython.org/)
    * **nb1** - Downloads SMILES and InChI string representations for all PubChem ids found in nb0. Uses [PubchemPy](https://pubchempy.readthedocs.io/en/latest/guide/gettingstarted.html)
    * **nb2** - Downloads SMILES and InChI string representations for all CAS numbers found `tables/nature_supplementary_cas_no.tsv`. These CAS numbers correspond to compounds predicted/found to bind MR1 and downloaded from the supplementary material of this [Nature Immunology paper](https://www.nature.com/articles/ni.3679)
    * **nb3** - Contains examples of functions that can help with calculation of SMILES similarities. RDKit is so far the best library for computing and comparing molecular fingerprints. However, RDKit does not have the [LINGO algorithm](https://www.ncbi.nlm.nih.gov/pubmed/15807504) implemented. The original publication can be found under `references`

* The `data` folder contains all text/string representations of the **reference** molecules (from the Nature paper) and the metabolites to **investigate**