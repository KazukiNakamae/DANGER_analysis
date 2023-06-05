import pandas as pd
import numpy as np
import sys

# コマンドライン引数の取得
seed = int(sys.argv[1])
input_csv = sys.argv[2]
output_csv = sys.argv[3]

# CSVファイルの読み込み
df = pd.read_csv(input_csv, header=None, names=['seq', 'contig', 'pos', 'match', 'strand', 'mm'], sep='\t')

# mm列目をシャッフル
df['mm'] = df.copy().sample(frac=1, random_state=seed).reset_index(drop=True)['mm']

# 結果を新しいCSVファイルに保存
df.to_csv(output_csv, index=False, header=False, sep='\t')