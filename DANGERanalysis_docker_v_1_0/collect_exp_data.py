#!/usr/local/bin/python3.9.1
# -*- coding: utf-8 -*-

# WT/Edited発現定量結果をまとめる

__author__ = "Kazuki Nakamae <kazkinakamae@gmail.com>"
__version__ = "1.00"
__date__ = "1 March 2023"

import argparse
import sys
import math
from os import makedirs
from os.path import join, abspath, isdir
import pandas as pd

#########################################

def main():
    parser = argparse.ArgumentParser(description='collect_MaChIAto_data',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # arguments
    parser.add_argument('-o','--outdir', type=str,  help='output directory', required=True)
    parser.add_argument('-w','--wt_sample_files', type=str,  help='files of WT sample', default=[], required=True, nargs='+')
    parser.add_argument('-e','--edited_sample_files', type=str,  help='files of Edited sample', default=[], required=True, nargs='+')
    parser.add_argument('-t','--threshold', type=float,  help='threshold of Edited/Wt ratio', required=True)
    args = parser.parse_args()

    # set input/output directory path
    OUTPUT_DIRECTORY = abspath(args.outdir)
    if not isdir(OUTPUT_DIRECTORY):
        print("Create " + OUTPUT_DIRECTORY + "...")
        makedirs(OUTPUT_DIRECTORY)
    # import pdb; pdb.set_trace() #---BREAKPOINT----------------------------------------------------------------------------------------------#

    id_pddf = pd.DataFrame()
    wt_expected_count_pddf = pd.DataFrame()
    wt_TPM_pddf = pd.DataFrame()
    wt_FPKM_pddf = pd.DataFrame()
    edited_expected_count_pddf = pd.DataFrame()
    edited_TPM_pddf = pd.DataFrame()
    edited_FPKM_pddf = pd.DataFrame()
    wt_sum_expected_count_pddf = pd.DataFrame()
    wt_sum_TPM_pddf = pd.DataFrame()
    wt_sum_FPKM_pddf = pd.DataFrame()
    edited_sum_expected_count_pddf = pd.DataFrame()
    edited_sum_TPM_pddf = pd.DataFrame()
    edited_sum_FPKM_pddf = pd.DataFrame()
    summary_expected_count_pddf = pd.DataFrame()
    summary_TPM_pddf = pd.DataFrame()
    summary_FPKM_pddf = pd.DataFrame()

    id_pddf["id"] = pd.read_csv(args.wt_sample_files[0], sep='\t')["transcript_id"]

    print("Collect data...")
    for i,wt_sample_file in enumerate(args.wt_sample_files):
      temp_pddf = pd.read_csv(wt_sample_file, sep='\t')
      wt_expected_count_pddf["ctrl_" + str(i+1)] = temp_pddf["expected_count"]
      wt_TPM_pddf["ctrl_" + str(i+1)] = temp_pddf["TPM"]
      wt_FPKM_pddf["ctrl_" + str(i+1)] = temp_pddf["FPKM"]

    wt_sum_expected_count_pddf["ctrl_ave"] = wt_expected_count_pddf.mean(axis=1)
    wt_sum_TPM_pddf["ctrl_ave"] = wt_TPM_pddf.mean(axis=1)
    wt_sum_FPKM_pddf["ctrl_ave"] = wt_FPKM_pddf.mean(axis=1)

    for i,edited_sample_file in enumerate(args.edited_sample_files):
      temp_pddf = pd.read_csv(edited_sample_file, sep='\t')
      edited_expected_count_pddf["edited_" + str(i+1)] = temp_pddf["expected_count"]
      edited_TPM_pddf["edited_" + str(i+1)] = temp_pddf["TPM"]
      edited_FPKM_pddf["edited_" + str(i+1)] = temp_pddf["FPKM"]

    edited_sum_expected_count_pddf["treated_ave"] = edited_expected_count_pddf.mean(axis=1)
    edited_sum_TPM_pddf["treated_ave"] = edited_TPM_pddf.mean(axis=1)
    edited_sum_FPKM_pddf["treated_ave"] = edited_FPKM_pddf.mean(axis=1)

    print("Calcutate Edited/WT ratio...")
    summary_expected_count_pddf["ON-ratio"] = edited_sum_expected_count_pddf["treated_ave"]/wt_sum_expected_count_pddf["ctrl_ave"]
    summary_TPM_pddf["ON-ratio"] = edited_sum_TPM_pddf["treated_ave"]/wt_sum_TPM_pddf["ctrl_ave"]
    summary_FPKM_pddf["ON-ratio"] = edited_sum_FPKM_pddf["treated_ave"]/wt_sum_FPKM_pddf["ctrl_ave"]

    def _put_labels(ratio, threshold):
        if math.isnan(ratio):
            return "unknown"
        elif ratio > threshold:
            return "upregulated"
        elif ratio < 1/threshold:
            return "downregulated"
        else:
            return "unchanged"

    print("Threshold: " + str(args.threshold) + "...")
    summary_expected_count_pddf["Exp"] = summary_expected_count_pddf["ON-ratio"].apply(_put_labels, threshold=float(args.threshold))
    summary_TPM_pddf["Exp"] = summary_TPM_pddf["ON-ratio"].apply(_put_labels, threshold=float(args.threshold))
    summary_FPKM_pddf["Exp"] = summary_FPKM_pddf["ON-ratio"].apply(_put_labels, threshold=float(args.threshold))

    id_pddf.reset_index()
    wt_expected_count_pddf.reset_index()
    wt_TPM_pddf.reset_index()
    wt_FPKM_pddf.reset_index()
    edited_expected_count_pddf.reset_index()
    edited_TPM_pddf.reset_index()
    edited_FPKM_pddf.reset_index()
    wt_sum_expected_count_pddf.reset_index()
    wt_sum_TPM_pddf.reset_index()
    wt_sum_FPKM_pddf.reset_index()
    edited_sum_expected_count_pddf.reset_index()
    edited_sum_TPM_pddf.reset_index()
    edited_sum_FPKM_pddf.reset_index()
    summary_expected_count_pddf.reset_index()
    summary_TPM_pddf.reset_index()
    summary_FPKM_pddf.reset_index()

    result_cnt_pddf = pd.concat([id_pddf, wt_expected_count_pddf, edited_expected_count_pddf, wt_sum_expected_count_pddf, edited_sum_expected_count_pddf, summary_expected_count_pddf], axis=1)
    result_TPM_pddf = pd.concat([id_pddf, wt_TPM_pddf, edited_TPM_pddf, wt_sum_TPM_pddf, edited_sum_TPM_pddf, summary_TPM_pddf], axis=1)
    result_FPKM_pddf = pd.concat([id_pddf, wt_FPKM_pddf, edited_FPKM_pddf, wt_sum_FPKM_pddf, edited_sum_FPKM_pddf, summary_FPKM_pddf], axis=1)

    print("Save results in " + OUTPUT_DIRECTORY + "...")

    # save
    result_cnt_pddf.to_csv(join(OUTPUT_DIRECTORY, 'ctrl_edited_fltrexpr_contig_expected_count_onratio.csv'), index=False, sep=",")
    result_TPM_pddf.to_csv(join(OUTPUT_DIRECTORY, 'ctrl_edited_fltrexpr_contig_tpm_onratio.csv'), index=False, sep=",")
    result_FPKM_pddf.to_csv(join(OUTPUT_DIRECTORY, 'ctrl_edited_fltrexpr_contig_fpkm_onratio.csv'), index=False, sep=",")

    print("Done.")
    sys.exit(0)

if __name__ == '__main__':
    main()
quit()