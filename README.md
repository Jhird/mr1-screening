# mr1-screening

## In silico screening of cancer metabolites that may bind MR1 receptor

### Instructions

* Start with the `notebooks` folder. Each notebook is numbered (0-indexed) in the order that they should be run. Briefly:

    * **nb0** - Downloads all metabolic compounds found in the KEGG database and tries to add as many external identifiers as possible (PubChem, CHEMBL, CAS numbers). Uses [BioPython](https://biopython.org/)
    * **nb1** - Downloads SMILES and InChI string representations for all PubChem ids found in nb0. Uses [PubchemPy](https://pubchempy.readthedocs.io/en/latest/guide/gettingstarted.html)
    * **nb2** - Downloads SMILES and InChI string representations for all CAS numbers found `tables/nature_supplementary_cas_no.tsv`. These CAS numbers correspond to compounds predicted/found to bind MR1 and downloaded from the supplementary material of this [Nature Immunology paper](https://www.nature.com/articles/ni.3679)