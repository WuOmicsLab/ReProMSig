# bash scripts/repromsig.sh ColoGuide_Stage_II_local/input/analysis.yaml ColoGuide_Stage_II_local/input/reporting.yaml

#!/usr/bin/bash
# parameters
if [ $# != 2 ] ; then 
    echo "Parameters: repromsig.sh analysis.yaml reporting.yaml"
    exit 1; 
fi

config_yaml_file=$1
config_report_yaml_file=$2

## get file/path from config.yaml, need to pip install shyaml------------------
scirpt_dir=$(cd `dirname $0`; pwd)
repromsig_dir=$(dirname $scirpt_dir)
echo $repromsig_dir

##
buiding_sig_dir=$(cat $config_yaml_file | shyaml get-value output_dir)
if [ ! -d $buiding_sig_dir ];then
    mkdir -p $buiding_sig_dir
fi
##
basename=$(basename $buiding_sig_dir)
buiding_sig_dir=$(cd `dirname $buiding_sig_dir`; pwd)
buiding_sig_dir=$buiding_sig_dir/$basename/
echo $buiding_sig_dir


## run the analysis codes ------------------
cd $repromsig_dir
## step0 
Rscript $scirpt_dir/yaml.process.R $config_yaml_file $repromsig_dir

## step1 Construct model to obtain signature ---
Rscript $scirpt_dir/model.analysis.R $buiding_sig_dir/sig.ini

## step2 Independence of signature from clinical features ---
Rscript $scirpt_dir/independence.analysis.R $buiding_sig_dir/sig.ini

## step3 input for Time-dependent ROC and prediction error curves ---
Rscript $scirpt_dir/performance.analysis.R $buiding_sig_dir/sig.ini

## step4 KM evaluation ---
Rscript $scirpt_dir/KM.evaluate.analysis.R $buiding_sig_dir/sig.ini

## step5 generate files for tripod report ---
Rscript $scirpt_dir/tripod.report.input.R $buiding_sig_dir/sig.ini

## step6 generate the tripod report html file ---
Rscript $scirpt_dir/tripod.report.html.R $config_report_yaml_file $buiding_sig_dir/sig.ini $buiding_sig_dir/tripod.ini

## step7 combine the rdatas and extract necessary datas ---
Rscript $scirpt_dir/local.rdata.generate.R $buiding_sig_dir/sig.ini $buiding_sig_dir/tripod.ini
