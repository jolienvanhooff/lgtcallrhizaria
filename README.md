# lgtcallrhizaria

Calling of lateral gene transfer (LGT) events in Rhizaria through detection (phylogeny-based) and prediction (species distribution- or classifier-based).

## LGT detection (detect_gene_origin.py)

### Usage

```bash
scripts/detect_gene_origin.py -l <orthogroups_members> -i iqtrees -o output_detect -c config_files/detect_config.ini
```

### Input

### Additional data

### Tools used 
* Python packages
  * ETE3 (3.1.2)
  * ETE Toolchain (3.0.0)
  * BioPython (1.78) 
* Notung (2.9.1.5) *path can be specified in config file*

## LGT prediction (predict_gene_origin.py)

### Usage

```bash
detect_gene_origin.py -i orthogroups.txt -i fasta_files -o output_predict -c config_files/predict_config.ini
```
 
### Input

### Additional data

### Tools used
* Python packages
  * ETE3 (3.1.2)
  * ETE Toolchain (3.0.0)
  * BioPython (1.78) 
  * scikit-learn (0.23.2)
  * joblib (1.0.1)
* MCL (14.137) *was applied before predict_gene_origin.py*
* DIAMOND (2.1.4) *path can be specified in config file*

## Reference
[Preprint: Lateral gene transfer leaves lasting traces in Rhizaria](https://www.biorxiv.org/content/10.1101/2023.01.27.525846v1)