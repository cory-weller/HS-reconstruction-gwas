# HS-reconstruction-gwas
 
This repository contains the scripts used to generate and process data, as well as generate figures, for the manuscript:

*Accurate, ultra-low coverage genome reconstruction and association studies in Hybrid Swarm mapping populations*

Cory A. Weller (caw5cv@virginia.edu) & Alan O. Bergland (aob2x@virginia.edu)

## Singularity

This workflow allows for Singularity containers to process data in a reproducible manner without installing required programs and libraries. You will first need to install singularity on your system, if it is not already available. Many HPC systems already have pre-loaded `singularity` that can be loaded as a module.

Otherwise, install singularity 3.x following the instructions from [sylabs.io](https://sylabs.io/guides/3.6/user-guide/quick_start.html#quick-installation-steps).

Then, you can retrieve the pre-built singularity image files from Singularity Hub. 

```bash
singularity pull --name harp.sif shub://cory-weller/HS-reconstruction-gwas:harp
singularity pull --name utils.sif shub://cory-weller/HS-reconstruction-gwas:utils
```