#!/bin/bash

eval "$(micromamba shell hook --shell=bash)"
micromamba activate get_exp_permutation_data_env;

# 入力ファイル名とループ回数の取得
input_csv=$1 # GRIN2B_exp_original/ctrl_edited_fltrexpr_contig_tpm_onratio.csv
start_count=$2 # 1
loop_count=$3 # 100
test_name=$4 # GRIN2B_exp

# シード値を1からループ回数までループ
for seed in $(seq $start_count $loop_count)
do
  echo "Generating Permutation Data ${seed}"
  # 出力ファイル名
  output_dir="${test_name}_p${seed}"
  mkdir ${output_dir};

  # Pythonスクリプトの実行
  python /tmp/get_exp_permutation_data.py ${seed} \
  ${input_csv} \
  ${output_dir}/ctrl_edited_fltrexpr_contig_tpm_onratio.csv
done

micromamba deactivate;
echo "Done."
exit 0