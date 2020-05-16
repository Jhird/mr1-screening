# mr1-screening

## In silico screening of cancer metabolites that may bind MR1 receptor

### Instructions

* Start with the `notebooks` folder. Each notebook is numbered (0-indexed) in the order that they should be run. Briefly:

    * **nb0** - Downloads all metabolic compounds found in the KEGG database and tries to add as many external identifiers as possible (PubChem, CHEMBL, CAS numbers). Uses [BioPython](https://biopython.org/)
    * **nb1** - Downloads SMILES and InChI string representations for all PubChem ids found in nb0. Uses [PubchemPy](https://pubchempy.readthedocs.io/en/latest/guide/gettingstarted.html)
    * **nb2** - Downloads SMILES and InChI string representations for all CAS numbers found `tables/nature_supplementary_cas_no.tsv`. These CAS numbers correspond to compounds predicted/found to bind MR1 and downloaded from the supplementary material of this [Nature Immunology paper](https://www.nature.com/articles/ni.3679)
    * **nb3** - Contains examples of functions that can help with calculation of SMILES similarities. RDKit is so far the best library for computing and comparing molecular fingerprints. However, RDKit does not have the [LINGO algorithm](https://www.ncbi.nlm.nih.gov/pubmed/15807504) implemented. The original publication can be found under `references`

* The `data` folder contains all text/string representations of the **reference** molecules (from the Nature paper) and the metabolites to **investigate** (from KEGG metabolic compounds database)

## `Preliminary results`

### [C00579](https://www.kegg.jp/entry/C00579)
* This compound is consumed/produced by [dehydrolipoamide dehydrogenase](https://www.kegg.jp/dbget-bin/www_bget?ec:1.8.1.4).
* It turns out that this enzyme has been targeted for treating cancer in multiple studies:
    * [Targeting the Achilles' heel of cancer cells via integrin-mediated delivery of ROS-generating dihydrolipoamide dehydrogenase.](https://www.ncbi.nlm.nih.gov/pubmed/30872792)
    * [Dihydrolipoamide dehydrogenase regulates cystine deprivation-induced ferroptosis in head and neck cancer.](https://www.ncbi.nlm.nih.gov/pubmed/31931284)
    * [Proteomic identification of dihydrolipoamide dehydrogenase as a target of autoantibodies in patients with endometrial cancer.](https://www.ncbi.nlm.nih.gov/pubmed/25202086)
    
### [C15519](https://www.kegg.jp/entry/C15519)
* This compound is 25-Hydroxycholesterol, which has been implicated in cancer progression as the studies below indicate. `NOTE` We might want to exclude this from analysis since it is a lipid, which is most likely presented by CD1, not MR1
    * [The Cholesterol Metabolite 25-Hydroxycholesterol Activates Estrogen Receptor α-Mediated Signaling in Cancer Cells and in Cardiomyocytes](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0016631)
    * [25-HC Decreases the Sensitivity of Human Gastric Cancer Cells to 5-fluorouracil and Promotes Cells Invasion via the TLR2/NF-κB Signaling Pathway](https://pubmed.ncbi.nlm.nih.gov/30664194/?from_term=25-hydroxycholesterol+cancer&from_pos=10)
    * [25-Hydroxycholesterol Promotes Migration and Invasion of Lung Adenocarcinoma Cells](https://pubmed.ncbi.nlm.nih.gov/30664194/?from_term=25-hydroxycholesterol+cancer&from_pos=10)
    * [25-Hydroxycholesterol Secreted by Macrophages in Response to Toll-like Receptor Activation Suppresses Immunoglobulin A Production](https://pubmed.ncbi.nlm.nih.gov/19805370/?from_term=25-hydroxycholesterol+cancer&from_pos=8)
    * [Oncogenic Roles of the Cholesterol Metabolite 25-hydroxycholesterol in Bladder Cancer](https://pubmed.ncbi.nlm.nih.gov/32382321/?from_term=25-hydroxycholesterol+cancer&from_pos=3)

### [C05295](https://www.kegg.jp/entry/C05295)
* This compound is an androgen hormone produced from testosterone by the enzyme [aromatase](https://www.kegg.jp/dbget-bin/www_bget?ec:1.14.14.14)
* There are several hormone-dependent cancers that are worth exploring
    * [Circulating Steroid Hormone Variations Throughout Different Stages of Prostate Cancer](https://pubmed.ncbi.nlm.nih.gov/28924064/?from_term=steroid+hormone+cancer&from_pos=1)
    * [Steroid Hormone Profiling in Human Breast Adipose Tissue Using Semi-Automated Purification and Highly Sensitive Determination of Estrogens by GC-APCI-MS/MS](https://pubmed.ncbi.nlm.nih.gov/29147745/?from_term=steroid+hormone+cancer&from_pos=6)

### [C15605](https://www.kegg.jp/entry/C15605)
* This metabolite is produced by the [indoleamine 2,3-dioxygenase (IDO) enzyme](https://www.kegg.jp/dbget-bin/www_bget?ec:1.13.11.52) which is part of the [Kynurenine pathway](https://www.kegg.jp/kegg-bin/show_pathway?ec00380+1.13.11.52). It turns out that the IDO enzyme is involved in cancer and other diseases
    * [Role of IDO and TDO in Cancers and Related Diseases and the Therapeutic Implications](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6584917/)
    * [The Kynurenine Pathway: A Primary Resistance Mechanism in Patients with Glioblastoma](http://ar.iiarjournals.org/content/37/5/2159.full)

### [C0018](https://www.kegg.jp/entry/C00118)
* This is Glyceraldehyde 3-phosphate, an obvious target since it is part of the Glycolytic Pathway, which has a long history of being implicated in cancer
    * [Deregulation of glycolysis in cancer: glyceraldehyde-3-phosphate dehydrogenase as a therapeutic target.](https://www.ncbi.nlm.nih.gov/pubmed/23445303)
    * [Glyceraldehyde-3-phosphate dehydrogenase promotes liver tumorigenesis by modulating phosphoglycerate dehydrogenase](https://www.ncbi.nlm.nih.gov/pubmed/28387968)
    * [Glyceraldehyde-3-phosphate dehydrogenase promotes cancer growth and metastasis through upregulation of SNAIL expression](https://www.ncbi.nlm.nih.gov/pubmed/27878251)
    * [Glyceraldehyde-3-phosphate dehydrogenase: a promising target for molecular therapy in hepatocellular carcinoma.](https://www.ncbi.nlm.nih.gov/pubmed/22964488)

### [C00169](https://www.kegg.jp/entry/C00169)
* This is carbamoyl phosphate and there is some evidence that the enzyme that produces it is implicated in cancer progression
    * [Caspase Recruitment Domain Family Member 10 Regulates Carbamoyl Phosphate Synthase 1 and Promotes Cancer Growth in Bladder Cancer Cells ](https://pubmed.ncbi.nlm.nih.gov/31565867/?from_term=carbamoyl+phosphate+cancer&from_pos=1)

### [C00085](https://www.kegg.jp/entry/C00085)
* This is F6P (Fructose 6-phosphate), a mejor metabolite in the Glycolytic pathway, which is very important in the Warburg effect in cancer. There is a universe of publications around this:
    * [NF-κB Upregulates glutamine-fructose-6-phosphate Transaminase 2 to Promote Migration in Non-Small Cell Lung Cancer](https://pubmed.ncbi.nlm.nih.gov/30885209/?from_term=fructose+6-phosphate+cancer&from_pos=2)
    * [Fructose 2,6-Bisphosphate in Cancer Cell Metabolism](https://pubmed.ncbi.nlm.nih.gov/30234009/?from_term=fructose+6-phosphate+cancer&from_pos=1)
    * [Decreased fructose-1,6-bisphosphatase-2 Expression Promotes Glycolysis and Growth in Gastric Cancer Cells](https://pubmed.ncbi.nlm.nih.gov/24063558/?from_term=fructose+6-phosphate+cancer&from_pos=9)
    * [Higher Autocrine Motility factor/glucose-6-phosphate Isomerase Expression Is Associated With Tumorigenesis and Poorer Prognosis in Gastric Cancer](https://pubmed.ncbi.nlm.nih.gov/30464597/?from_term=fructose+6-phosphate+cancer&from_pos=10)

### [C05296](https://www.kegg.jp/entry/C05296)
* Same as above, this is an androgen steroid

### [C05299](https://www.kegg.jp/entry/C05299)
* Similar to the previous two compounds above, this is an estrogen steroid

### [C04295](https://www.kegg.jp/entry/C04295)
* Similar as above, this is an androgen steroid

### [C00951](https://www.kegg.jp/entry/C00951)
* Estradiol, another steroid involved in prostate and breast cancer signaling pathways