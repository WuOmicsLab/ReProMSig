# ReProMSig

## Description

[Reproducible Prognosis Molecular Signature (ReProMSig)] (https://omics.bjcancer.org/prognosis/) platform could help develop and validate a multivariable prognostic/predictive biomarker in a transparent and reproducible way, with the following advanced features:

It streamlines the analysis process in development of a multivariable prediction model using molecular profiles and/or clinicopathological factors, as well as evaluation of its prognostic and/or predictive value.

The full detail of modelling procedures and results can be provided as a signature report file (example report), which is well-designed following the TRIPOD statement.

Long-term storage and management of user datasets and signatures are supported for registered users.

Risk assessment for single patient using a developed biomarker is provided for research purpose.

This repository hosts the source codes of ReProMSig analysis pipeline, which can be used for local analysis. It includes scripts ($(pwd)/scripts/), and an example folder $(pwd)/example/ including demo input files and output files.   


## Installation
<b>R</b> (>= 3.6.1) and <b>Python</b> are required to be installed. 

Run the installiation script to install the required packages automatically before using the pipeline. $(pwd) referes to the ReProMSig direcotry path.

```bash
Rscript $(pwd)/package.install.R
```

To process the user provided config files (YAML format), yaml package needs to be installed.
```bash
pip install shyaml
```
