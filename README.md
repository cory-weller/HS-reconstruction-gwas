# HS-reconstruction-gwas
 
This repository contains the scripts used to generate and process data, as well as generate figures, for the manuscript:

*Accurate, ultra-low coverage genome reconstruction and association studies in Hybrid Swarm mapping populations*

Cory A. Weller (caw5cv@virginia.edu) & Alan O. Bergland (aob2x@virginia.edu)

## Downloading singularity containers
This workflow allows for Singularity containers to process data in a reproducible manner without installing required programs and libraries. The singularity recipe files `Singularity.harp` and `Singularity.utils` automatically push pre-built image files to Singularity Hub. Singularity must be installed on your system. Once singularity is installed, the pre-built image files can be pulled to this repository using the commands:

```bash
singularity pull --name harp.sif shub://cory-weller/HS-reconstruction-gwas:harp
singularity pull --name utils.sif shub://cory-weller/HS-reconstruction-gwas:utils
```