from os.path import join
from os import mkdir
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np 
import sys
import math


table_fn_list = sys.argv[1:]

# import pdb; pdb.set_trace()

# テーブル読み取り
def mergeEnrichmentTable(input_fp, table_fn_list):
  enriment_table_cname = ["GO.ID", "Term", "Annotated", "Significant", "Expected", "pvalue", "mm"]
  enriment_table = pd.DataFrame(columns = enriment_table_cname)
  for table_fn in table_fn_list:
    temp_enriment_table = pd.read_table(join(table_fn, input_fp), sep= "\t", names=enriment_table_cname)
    enriment_table = pd.concat([enriment_table, temp_enriment_table], axis=0)
    temp_enriment_table = pd.DataFrame(columns = enriment_table_cname)
  return enriment_table

bp_enriment_table = mergeEnrichmentTable("BP_enrichment_table.txt", table_fn_list)
cc_enriment_table = mergeEnrichmentTable("CC_enrichment_table.txt", table_fn_list)
mf_enriment_table = mergeEnrichmentTable("MF_enrichment_table.txt", table_fn_list)

def calcDindexTable(enriment_table):
  # DANGER index算出
  danger_table_cname = ["GO","Term", "Dindex"]
  danger_table = pd.DataFrame(columns=danger_table_cname)
  def row_calcDindex(row):
    d_index_ele = float(row["Annotated"])*pow(math.e, 4-float(row["mm"]))
    return d_index_ele
  for GO in enriment_table["GO.ID"].unique().tolist():
    # 特定のGOをもつエントリー一覧取得
    go_table = enriment_table[enriment_table['GO.ID'] == GO]
    # ミスマッチ数ごとにD indexの子要素を計算
    d_index_ele_list = go_table.apply(row_calcDindex, axis=1)
    description = enriment_table[enriment_table['GO.ID'] == GO]["Term"].iloc[0]
    d_index = d_index_ele_list.sum() # D-indexを計算
    # 集計テーブルに追加
    temp_danger_table = pd.DataFrame(dict(GO=[GO], Term=[description], Dindex=[d_index]))
    danger_table = pd.concat([danger_table, temp_danger_table], axis=0)
    temp_danger_table = pd.DataFrame(columns=danger_table_cname)
  sorted_danger_table = danger_table.sort_values('Dindex', ascending=False)
  return sorted_danger_table

bp_danger_table = calcDindexTable(bp_enriment_table)
cc_danger_table = calcDindexTable(cc_enriment_table)
mf_danger_table = calcDindexTable(mf_enriment_table)

mkdir('DANGER_analysis_result')
bp_danger_table.to_csv(join('DANGER_analysis_result', "DANGER_index_on_Biological_Process.txt"), header=True, index=False, sep='\t')
cc_danger_table.to_csv(join('DANGER_analysis_result', "DANGER_index_on_Cellular_Component.txt"), header=True, index=False, sep='\t')
mf_danger_table.to_csv(join('DANGER_analysis_result', "DANGER_index_on_Molecular_Function.txt"), header=True, index=False, sep='\t')

def makeDangerHistgram(danger_table, title, label):
  plt.hist(danger_table['Dindex'].to_list(), bins=int(danger_table['Dindex'].max()+1))
  plt.title(title)
  plt.xlabel('D-index')
  plt.ylabel('Gene Counts')
  plt.savefig(join('DANGER_analysis_result', label + "histogram.tiff"), format="tiff", dpi=350)
  plt.clf()

makeDangerHistgram(bp_danger_table, "Biological Process", "BP_GO_")
makeDangerHistgram(cc_danger_table, "Cellular Component", "CC_GO_")
makeDangerHistgram(mf_danger_table, "Molecular Function", "MF_GO_")

def makeDangerHistgram(danger_table):
  dindex_q75, dindex_q25 = np.percentile(danger_table['Dindex'].to_list(), [75 ,25])
  dindex_iqr = dindex_q75 - dindex_q25
  upper_bound = dindex_q75 + (dindex_iqr * 1.5)
  return [upper_bound, danger_table[danger_table['Dindex'].to_list()>upper_bound]]

bp_upper_bound, bp_high_danger_table = makeDangerHistgram(bp_danger_table)
cc_upper_bound, cc_high_danger_table = makeDangerHistgram(cc_danger_table)
mf_upper_bound, mf_high_danger_table = makeDangerHistgram(mf_danger_table)

bp_high_danger_table.to_csv(join('DANGER_analysis_result', "High_risk_DANGER_index_on_Biological_Process.txt"), header=True, index=False, sep='\t')
cc_high_danger_table.to_csv(join('DANGER_analysis_result', "High_risk_DANGER_index_on_Cellular_Component.txt"), header=True, index=False, sep='\t')
mf_high_danger_table.to_csv(join('DANGER_analysis_result', "High_risk_DANGER_index_on_Molecular_Function.txt"), header=True, index=False, sep='\t')

print("Biological Process")
print("Sum of DANGER indexes")
print(bp_danger_table['Dindex'].sum())
print("Upper boundary")
print(bp_upper_bound)

print("Cellular Component")
print("Sum of DANGER indexes")
print(cc_danger_table['Dindex'].sum())
print("Upper boundary")
print(cc_upper_bound)

print("Molecular Function")
print("Sum of DANGER indexes")
print(mf_danger_table['Dindex'].sum())
print("Upper boundary")
print(mf_upper_bound)


# import pdb; pdb.set_trace()

