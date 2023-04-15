# ReProMSig

## Description

Reproducible Prognosis Molecular Signature (ReProMSig) (https://omics.bjcancer.org/prognosis/) platform could help develop and validate a multivariable prognostic/predictive biomarker in a transparent and reproducible way, with the following advanced features:

- It streamlines the analysis process in development of a multivariable prediction model using molecular profiles and/or clinicopathological factors, as well as evaluation of its prognostic and/or predictive value.

- The full detail of modelling procedures and results can be provided as a signature report file (example report), which is well-designed following the TRIPOD statement.

- Long-term storage and management of user datasets and signatures are supported for registered users.

- Risk assessment for single patient using a developed biomarker is provided for research purpose.

This repository hosts source codes of ReProMSig analysis pipeline, which can be used for local analysis. The generated signature models could also be uploaded to the [ReProMSig web server](https://omics.bjcancer.org/prognosis/) for sharing to public.


## Installation
System requirements: <b>R >= 3.6.1</b> and <b>Python</b>.

1) Run the installiation script to install the required packages automatically before using the pipeline. $(pwd) referes to the ReProMSig direcotry path.

```bash
Rscript $(pwd)/package.install.R
```

2) Install [shyaml](https://github.com/0k/shyaml) package for processing user provided config files (YAML format).
```bash
pip install shyaml
```

## Prepare data before running

The main function $(pwd)/scripts/repromsig.sh takes two config files in YAML format as input. And user also should provide clinicopathological and/or molecular datasets that will be taken as training and validation cohort(s).

###  1) Config file for analysis (YAML format)
This YAML file consists of information for data path and analysis parameters. Please see$(pwd)/example/analysis.yaml for an example file and xx for a complete list of configurations.

###  2) Config file for reporting (YAML format)
This YAML file consists of information for generating the reporting file of the developed siganture, including generanl descriptions and the structured items according to the TRIPOD guideline. Please see$(pwd)/example/reporting.yaml for an example.

### 3) Clinicopathological / Molecular profile files
Please visit [ReProMsig tutorial](https://omics.bjcancer.org/prognosis/) (section '1.1 Private datasets') for details.
#### Patient annotation
Patient annotation file consists of three groups of columns including fixed columns, endpoint columns and custom columns. The names of fixed columns should be identical to the table template file. These clinicopathological parameters will be used for query samples suitable for analysis, as well as for prognosis model development. Please note that missing values in all columns should be provided as "NA". You can generate the formatted patient annotation file by modifying the template file.

#### Molecular profiles
Molecular profile file is a feature-by-sample matrix in TXT/CSV/Excel format. Feature IDs should be provided at the first column with a column name "ID". The additional columns are molecular features for each sample with sample identifiers as column names. Molecular features could be gene mutation status, mRNA/non-coding RNA/protein quantification levels, methylation levels, etc. 
