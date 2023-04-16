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

####  1) Config file for analysis (YAML format)
This YAML file consists of information for data path and analysis parameters. Please see$(pwd)/example/analysis.yaml for an example file and xx for a complete list of configurations.

####  2) Config file for reporting (YAML format)
This YAML file consists of information for generating the reporting file of the developed siganture, including generanl descriptions and the structured items according to the TRIPOD guideline. Please see$(pwd)/example/reporting.yaml for an example.

#### 3) Clinicopathological / Molecular profile files
Please visit [ReProMsig tutorial](https://omics.bjcancer.org/prognosis/) (section '1.1 Private datasets') for details.

- <b>Patient annotation</b>
Patient annotation file consists of three groups of columns including fixed columns, endpoint columns and custom columns. The names of fixed columns should be identical to the table template file. These clinicopathological parameters will be used for query samples suitable for analysis, as well as for prognosis model development. Please note that missing values in all columns should be provided as "NA". You can generate the formatted patient annotation file by modifying the template file.

- <b>Molecular profiles</b>
Molecular profile file is a feature-by-sample matrix in TXT/CSV/Excel format. Feature IDs should be provided at the first column with a column name "ID". The additional columns are molecular features for each sample with sample identifiers as column names. Molecular features could be gene mutation status, mRNA/non-coding RNA/protein quantification levels, methylation levels, etc. 

## Usages 
```bash
# The main function is repromsig.sh, which takes two aforementioned yaml files as input. 
bash scripts/repromsig.sh $analysis.yaml $reporting.yaml

# running example project
bash scripts/repromsig.sh ColoGuide_Stage_II_local/input/analysis.yaml  ColoGuide_Stage_II_local/input/reporting.yaml
```

This analysis will create multiple result folders containing output RData files, tables and plots described here.

<b>Note</b>: the RData file and Reporting html file in the <b>upload</b> sub-directory are the  core output files that could be uploaded to "My signature" module of ReProMSig (only for registered users) for displaying and sharing.

## Script details
`repromsig.sh` utilizes multiple scripts to perform data processing, extracting, modeling and reporting, see below for details:

Script |Description
:-|:-
ymal.process.R | Data processing, analysis parameters extraction from the user-provided "analysis yaml file".
model.analysis.R | Predictors selection, multivariable prediction model building, signature score calculation and patient risk group stratification
independence.analysis.R |	Perform univariate and multivariate Cox regression analyses, to test whether the signature is an independent prognostic or predictive factor. 
performance.analysis.R | Model evaluation, including time dependent receiver operating characteristic (ROC), prediction error (PE) , and calibration analysis.
KM.evaluate.analysis.R | Inspect the survival differences between risk groups by Kaplan-Meier analysis and log-rank test. 
tripod.report.input.R | Generate the variables, tables and graphs for a signature that will be displayed on the web and in the reporting file.
tripod.report.html.R | Save out the reporting html file by extract data from analysis output and used-provided "tripod yaml file".
local.rdata.generate.R | Save out the RData file that could be uploaded to the ReProSig website, for displaying and sharing.

## License
ReProMSig is free for academic users of non-commercial purposes. Commercial use of ReProMSig requires a license. If ReProMSig package was used for your analysis, please cite our package.

## Contact information
Lihua Cao (lihuacao@bjcancer.org), Tingting Zhao, Jianmin Wu (wujm@bjmu.edu.cn).


