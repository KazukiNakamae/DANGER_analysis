
#!/usr/local/bin/python3.9.1
# -*- coding: utf-8 -*-

# GOtermごとのt分布信頼区間を求める

__author__ = "Kazuki Nakamae <kazkinakamae@gmail.com>"
__version__ = "1.00"
__date__ = "24 May 2023"

import sys
import argparse
import os
from os import makedirs
import pandas as pd
import scipy.stats as stats
import matplotlib.pyplot as plt
import seaborn as sns

#########################################
def main():
    parser = argparse.ArgumentParser(description='collect_MaChIAto_data',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # arguments
    parser.add_argument('-t','--target_result_dir', type=str,  help='target directory', required=True)
    parser.add_argument('-n','--neg_collection_dir', type=str,  help='negative control directory', required=True)
    parser.add_argument('-c','--confidence_level', type=float,  help='threshold of confidence level', default=0.95, required=True)
    parser.add_argument('-o','--output_dir', type=str,  help='output directory', required=True)
    args = parser.parse_args()
    if not os.path.isdir(args.output_dir):
        os.makedirs(args.output_dir)
        print("Create Output Foloder: " + args.output_dir)

    def get_dindex_list(file):
        print("Read " + file)
        df = pd.read_csv(file, sep="\t")
        return df

    def calculate_confidence_interval(data_dir, input_fn, confidence_level=0.95):
        # データフレームを格納するリストを作成します
        df_list = []

        # 指定ディレクトリ内のすべてのCSVファイルを読み込みます
        for dirname in os.listdir(data_dir):
            bp_dindex_file = os.path.join(data_dir, dirname, "DANGER_analysis_result/" + input_fn)
            if os.path.isfile(bp_dindex_file):
                df_list.append(get_dindex_list(bp_dindex_file))

        # すべてのデータフレームを結合します
        combined_df = pd.concat(df_list)

        # データフレームを1番目の列でグループ化します
        grouped = combined_df.groupby(combined_df.columns[0])
        # 各グループに対して、2番目の列の信頼区間を計算します
        confidence_intervals = {}
        for label, group_df in grouped:
            # t分布を使用して信頼区間を計算します
            mean = group_df.iloc[:, 2].mean()
            sem = stats.sem(group_df.iloc[:, 2])
            ci = stats.t.interval(confidence_level, len(group_df.iloc[:, 2]) - 1, loc=mean, scale=sem) # 信頼区間
            confidence_intervals[label] = ci

        return confidence_intervals

    def plot_bar(df, output_fn):
        # データフレームを数値情報でソート（降順）
        df_sorted = df.sort_values(by=df.columns[2], ascending=False)
        # グラデーションの色パレットを作成
        palette = reversed(sns.color_palette("coolwarm", len(df_sorted)))
        # 棒グラフの作成
        if len(df_sorted) > 500:
            print("Sorry. Significant GO terms is too large. The image cannot be automatically made.")
            return 0
        if len(df_sorted) < 66:
            hight = 10
        else:
            hight = 0.15*len(df_sorted)
        plt.figure(figsize=(10, hight))
        sns.barplot(x=df.columns[2], y=df.columns[0], data=df_sorted, palette=palette)
        # tiffファイルとして保存
        plt.savefig(output_fn, format="tiff", dpi=350)

    print("D-index on Biological Process...")
    # サンプルデータのD-indexデータ取得
    target_bp_dindex_file = os.path.join(args.target_result_dir, "DANGER_analysis_result/DANGER_index_on_Biological_Process.txt")
    if os.path.isfile(target_bp_dindex_file):
        print("Read " + target_bp_dindex_file)
        target_bp_dindex_df = get_dindex_list(target_bp_dindex_file)
    else:
        print("The target file is not found. You shoud check the input paramters.")
        print("Done.")
        sys.exit(0)
    # ネガコンデータのD-indexの信頼区間を算出
    bp_negative_borders = calculate_confidence_interval(args.neg_collection_dir, "DANGER_index_on_Biological_Process.txt", args.confidence_level)
    target_bp_dindex_df['Confidence_Interval_Lower_on_Null_Distribution'] = \
        target_bp_dindex_df[target_bp_dindex_df.columns[0]].to_frame().applymap(lambda x, bp_negative_borders: bp_negative_borders.get(x, [None])[0], bp_negative_borders=bp_negative_borders)
    target_bp_dindex_df['Confidence_Interval_Upper_on_Null_Distribution'] = \
        target_bp_dindex_df[target_bp_dindex_df.columns[0]].to_frame().applymap(lambda x, bp_negative_borders: bp_negative_borders.get(x, [None])[1], bp_negative_borders=bp_negative_borders)
    # import pdb; pdb.set_trace()
    target_bp_dindex_df.to_csv(os.path.join(args.output_dir, "All_DANGER_index_on_Biological_Process.txt"), index=False, header=True, sep='\t')
    target_bp_dindex_df.query('Dindex > Confidence_Interval_Upper_on_Null_Distribution').to_csv(os.path.join(args.output_dir, "Significant_DANGER_index_on_Biological_Process.txt"), index=False, header=True, sep='\t')
    plot_bar(target_bp_dindex_df.query('Dindex > Confidence_Interval_Upper_on_Null_Distribution'), os.path.join(args.output_dir, "Significant_DANGER_index_on_Biological_Process.tiff"))

    print("D-index on Cellular Component...")
    # サンプルデータのD-indexデータ取得
    target_cc_dindex_file = os.path.join(args.target_result_dir, "DANGER_analysis_result/DANGER_index_on_Cellular_Component.txt")
    if os.path.isfile(target_cc_dindex_file):
        print("Read " + target_cc_dindex_file)
        target_cc_dindex_df = get_dindex_list(target_cc_dindex_file)
    else:
        print("The target file is not found. You shoud check the input paramters.")
        print("Done.")
        sys.exit(0)
    # ネガコンデータのD-indexの信頼区間を算出
    cc_negative_borders = calculate_confidence_interval(args.neg_collection_dir, "DANGER_index_on_Cellular_Component.txt", args.confidence_level)
    target_cc_dindex_df['Confidence_Interval_Lower_on_Null_Distribution'] = \
        target_cc_dindex_df[target_cc_dindex_df.columns[0]].to_frame().applymap(lambda x, cc_negative_borders: cc_negative_borders.get(x, [None])[0], cc_negative_borders=cc_negative_borders)
    target_cc_dindex_df['Confidence_Interval_Upper_on_Null_Distribution'] = \
        target_cc_dindex_df[target_cc_dindex_df.columns[0]].to_frame().applymap(lambda x, cc_negative_borders: cc_negative_borders.get(x, [None])[1], cc_negative_borders=cc_negative_borders)
    # import pdb; pdb.set_trace()
    target_cc_dindex_df.to_csv(os.path.join(args.output_dir, "All_DANGER_index_on_Cellular_Component.txt"), index=False, header=True, sep='\t')
    target_cc_dindex_df.query('Dindex > Confidence_Interval_Upper_on_Null_Distribution').to_csv(os.path.join(args.output_dir, "Significant_DANGER_index_on_Cellular_Component.txt"), index=False, header=True, sep='\t')
    plot_bar(target_cc_dindex_df.query('Dindex > Confidence_Interval_Upper_on_Null_Distribution'),  os.path.join(args.output_dir, "Significant_DANGER_index_on_Cellular_Component.tiff"))
    
    print("D-index on Molecular Function...")
    target_mf_dindex_file = os.path.join(args.target_result_dir, "DANGER_analysis_result/DANGER_index_on_Molecular_Function.txt")
    if os.path.isfile(target_mf_dindex_file):
        print("Read " + target_mf_dindex_file)
        target_mf_dindex_df = get_dindex_list(target_mf_dindex_file)
    else:
        print("The target file is not found. You shoud check the input paramters.")
        print("Done.")
        sys.exit(0)
    # ネガコンデータのD-indexの信頼区間を算出
    mf_negative_borders = calculate_confidence_interval(args.neg_collection_dir, "DANGER_index_on_Molecular_Function.txt", args.confidence_level)
    target_mf_dindex_df['Confidence_Interval_Lower_on_Null_Distribution'] = \
        target_mf_dindex_df[target_mf_dindex_df.columns[0]].to_frame().applymap(lambda x, mf_negative_borders: mf_negative_borders.get(x, [None])[0], mf_negative_borders=mf_negative_borders)
    target_mf_dindex_df['Confidence_Interval_Upper_on_Null_Distribution'] = \
        target_mf_dindex_df[target_mf_dindex_df.columns[0]].to_frame().applymap(lambda x, mf_negative_borders: mf_negative_borders.get(x, [None])[1], mf_negative_borders=mf_negative_borders)
    # import pdb; pdb.set_trace()
    target_mf_dindex_df.to_csv(os.path.join(args.output_dir, "All_DANGER_index_on_Molecular_Function.txt"), index=False, header=True, sep='\t')
    target_mf_dindex_df.query('Dindex > Confidence_Interval_Upper_on_Null_Distribution').to_csv(os.path.join(args.output_dir, "Significant_DANGER_index_on_Molecular_Function.txt"), index=False, header=True, sep='\t')
    plot_bar(target_mf_dindex_df.query('Dindex > Confidence_Interval_Upper_on_Null_Distribution'),  os.path.join(args.output_dir, "Significant_DANGER_index_on_Molecular_Function.tiff"))

    print("Done.")
    sys.exit(0)

if __name__ == '__main__':
    main()
quit()