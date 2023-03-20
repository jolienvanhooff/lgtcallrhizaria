# lgtcallrhizaria

Calling of lateral gene transfer (LGT) events in Rhizaria through detection (phylogeny-based) and prediction (species distribution- or classifier-based).

[**Lateral gene transfer leaves lasting traces in Rhizaria**](https://www.biorxiv.org/content/10.1101/2023.01.27.525846v1)

----

## LGT detection (detect_gene_origin.py)

### Usage

```bash
scripts/detect_gene_origin.py -l ORTHOGROUPS_MEMBERS_LIST -i IQTREES_FOLDER -o OUTPUT -c CONFIG
```
Replace the following: 
* ORTHOGROUPS_MEMBERS_LIST: list containing names of orthogroups, followed by a colon, followed by the identifiers of sequence members, such as: 
    * OG0208575: SorSpe009682_OG0208575 SorSpe014733_OG0208575 SorSpe028706_OG0208575
    * NB: this is the Orthogroups.txt format outputted by OrthoFinder. 
* IQTREES_FOLDER: directory containing the newick-formatted trees that are to be analyzed
* OUTPUT: name of output directory (default: output_detect)
* CONFIG: configuration file in which settings and additional input files are specified (default: config_files/detect_config.ini)

It is possible to restart detect_gene_origin.py by simply entering the same command. For a complete restart, the output directory should be removed before.

### Input and additional data
The data used for this project can be obtained online:
* ORTHOGROUPS_MEMBERS_LIST: [Figshare project: Lateral gene transfers in Rhizaria](https://figshare.com/projects/Lateral_gene_transfers_LGTs_in_Rhizaria/158240) - orthogroups_members_list.txt  
* IQTREES_FOLDER: [Figshare project: Lateral gene transfers in Rhizaria](https://figshare.com/projects/Lateral_gene_transfers_LGTs_in_Rhizaria/158240) - SupplementaryDataset1.tar.gz
* CONFIG: in this repository, config_files/detect_config.ini
All input required in CONFIG is found in this repository (folder: data), except: 
* notung_path: here one has to specify where the notung algorithm can be found

### Tools used 
* Python packages
  * ETE3 (3.1.2)
  * BioPython (1.78) 
* Notung (2.9.1.5) *path to Notung can be specified in config file, under 'notung_path'*

----

## LGT prediction (predict_gene_origin.py)

### Usage & input

```bash
detect_gene_origin.py -l ORTHOGROUPS -f FASTA_FILES -i IQTREES_FOLDER -o OUTPUT -c CONFIG
```
Replace the following: 
* ORTHOGROUPS: list containing names of orthogroups
* IQTREES_FOLDER: directory containing the newick-formatted trees that are to be analyzed
* OUTPUT: name of output directory (default: output_detect)
* CONFIG: configuration file in which settings and additional input files are specified (default: config_files/detect_config.ini)

It is possible to restart detect_gene_origin.py by simply entering the same command. For a complete restart, the output directory should be removed before.

### Additional data

### Tools used
* Python packages
  * ETE3 (3.1.2)
  * BioPython (1.78) 
  * scikit-learn (0.23.2)
  * joblib (1.0.1)
* MCL (14.137) *for efficiency reasons, this can be execulted before before predict_gene_origin.py, provided that the correct file organization is being used (in the output folder, in subfolder 'OGXXXX_clusters', where OGXXXX is the name of the orthogroup, for example OG1253601)*
* DIAMOND (2.1.4) *path can be specified in config file*

----

## Manuscript
[Preprint: Lateral gene transfer leaves lasting traces in Rhizaria](https://www.biorxiv.org/content/10.1101/2023.01.27.525846v1)

## Questions & comments
[jolien.vanhooff@wur.nl](mailto:jolien.vanhooff@wur.nl)



