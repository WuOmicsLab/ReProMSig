# ReProMSig

## Description

<b>R</b>eproducible <b>P</b>rognosis <b>M</b>olecular <b>Sig</b>nature (ReProMSig) (https://omics.bjcancer.org/prognosis/) platform could help develop and validate a multivariable prognostic/predictive biomarker in a transparent and reproducible way, with the following advanced features:

- It streamlines the analysis process in development of a multivariable prediction model using molecular profiles and/or clinicopathological factors, as well as evaluation of its prognostic and/or predictive value.

- The full detail of modelling procedures and results can be provided as a signature report file (example report), which is well-designed following the TRIPOD statement.

- Long-term storage and management of user datasets and signatures are supported for registered users.

- Risk assessment for single patient using a developed biomarker is provided for research purpose.

This repository hosts source code of ReProMSig analysis pipeline, which can be used for local analysis. The generated signature models could also be uploaded to the [ReProMSig web server](https://omics.bjcancer.org/prognosis/) for sharing to public.


## Installation
System requirements: <b>R >= 3.6.1</b> and <b>Python</b>.

1) Run the installiation script to install the required packages automatically before using the pipeline.

```bash
Rscript scripts/package.install.R
```
  <b>Note</b>: if automatic installation fails for some packages, please try manual installation with the failed package as below:
```R
# Linux users: specify the specific version, e.g.
remotes::install_version("glmnet", version = "3.0-2", repos = "https://cran.us.r-project.org")

# Mac/Windows users: install the package in binary format (type = "binary"), e.g.
install.packages('glmnet', type='binary')
```

2) Install [shyaml](https://github.com/0k/shyaml) package for processing user provided config files (YAML format).
```bash
pip install shyaml
```

3) Install [pandoc](https://www.pandoc.org/](https://www.pandoc.org/installing.html) for converting the RMarkdown document to a HTML-format reporting file
pandoc shoud be installed and version 1.12.3 or higher is required.

## Prepare data before running

The portal script `repromsig.sh` takes two config files in YAML format as input. User need to provide clinicopathological and/or molecular profiles that will be used as training and validation cohort(s).

####  1) Config file for analysis (YAML format)
This YAML file consists of data path and analysis parameters. Please see `ColoGuide_Stage_II_local/input/analysis.yaml` for an example and `config/analysis.default.setting.yaml` for a complete list of configurations.

####  2) Config file for reporting (YAML format)
This YAML file consists of structured information needed for generating the reporting file of a developed siganture, according to the [TRIPOD guideline](https://www.tripod-statement.org/). Please see `ColoGuide_Stage_II_local/input/reporting.yaml` for an example.

#### 3) Clinicopathological / Molecular profile files
Please visit [ReProMsig tutorial](https://omics.bjcancer.org/prognosis/) (section '1.1 Private datasets') for file format details.

- <b>Patient annotation</b>
Patient annotation file consists of three groups of columns including fixed columns, endpoint columns and custom columns. The names of fixed columns should be identical to the table template file. These clinicopathological parameters will be used for query samples suitable for analysis, as well as for prognosis model development. Please note that missing values in all columns should be provided as "NA". You can generate the formatted patient annotation file by modifying the template file.

- <b>Molecular profiles</b>
Molecular profile file is a feature-by-sample matrix in TXT/CSV/Excel format. Feature IDs should be provided at the first column with a column name "ID". The additional columns are molecular features for each sample with sample identifiers as column names. Molecular features could be gene mutation status, mRNA/non-coding RNA/protein quantification levels, methylation levels, etc. 

## Usages 
```bash
# The portal script is repromsig.sh, which takes two aforementioned yaml files as input. 
bash scripts/repromsig.sh [analysis.yaml] [reporting.yaml]

# Running an example project
bash scripts/repromsig.sh ColoGuide_Stage_II_local/input/analysis.yaml  ColoGuide_Stage_II_local/input/reporting.yaml
```

This analysis will create multiple output sub-folders, including information extracted from analysis.yaml file (`rda` dir), output files from signature modelling (`model`dir), independence test (`independence` dir), discrimination and calibration  evaluation (`performance` dir), survival differences inspection between risk groups (`external_evaluate` dir), summary tables and figures  for TRIPOD reporting (`tripod` dir) and the RData file and reporting html file (<b>`upload`</b> dir) that could be uploaded to "My signature" module of [ReProMSig web server](https://omics.bjcancer.org/prognosis/) for sharing.

## Script details
`repromsig.sh` utilizes multiple scripts to perform data processing, extracting, modelling and reporting, as shown below:

Script |Description
:-|:-
ymal.process.R | Data processing and extract analysis parameters from the user-provided `analysis.yaml`.
model.analysis.R | Perform predictor selection, multivariable prediction model building, signature score calculation and patient risk group stratification.
independence.analysis.R |	Perform univariate and multivariate Cox regression analyses, to test whether the signature is an independent prognostic or predictive factor. 
performance.analysis.R | Perform model evaluation, including time dependent receiver operating characteristic (ROC), prediction error (PE) , and calibration analysis.
KM.evaluate.analysis.R | Inspect the survival differences between risk groups by Kaplan-Meier analysis and log-rank test. 
tripod.report.input.R | Generate the variables, tables and graphs for a signature that will be shown in the reporting file.
tripod.report.html.R | Export the reporting html file by extract data from analysis output and used-provided "tripod yaml file".
local.rdata.generate.R | Export the RData file that could be uploaded to the ReProSig website, for displaying and sharing.

## License
ReProMSig is free for academic users of non-commercial purposes. Commercial use of ReProMSig requires a license. If ReProMSig package was used for your analysis, please cite our package.

## Contact information
Lihua Cao (lihuacao@bjcancer.org), Tingting Zhao(zhaott@bjmu.edu.cn), Jianmin Wu (wujm@bjmu.edu.cn).



