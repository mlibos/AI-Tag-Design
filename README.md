# AI-Tag-Design
AI Design Project @IPI for generating new epitope tags and nanobodies that bind to them. 

Frameworks used by our pipeline:
* [**ProteinMPNN**](https://github.com/dauparas/ProteinMPNN/tree/main)
* [**LigandMPNN**](https://github.com/dauparas/LigandMPNN)
* [**BoltzDesign1**](https://github.com/yehlincho/BoltzDesign1)

**De Novo Generated Binders for Generated Epitope Tags**
![alt text](https://github.com/mlibos/AI-Tag-Design/blob/main/boltzdesign1_models/best_models_pymol.png)

## Current Pipeline:
* Redesign of epitope tags using existing tags as templates
* Use ProteinMPNN as the initial redesign framework 
* Tested different initial inputs consisting of tags/binders
* Generation of 10,000 candidates from best initial input
* Filtering of candidate tags by biophysical properties
* Highest diversity candidate tags selected for binder design
* Binders were designed via two pathways: LigandMPNN and BoltzDesign1
* LigandMPNN was used to redesign binding residues in original binder
* BoltzDesign1 was used for de novo binder design
* Highest quality designs were selected based on iptm and confidence scores
