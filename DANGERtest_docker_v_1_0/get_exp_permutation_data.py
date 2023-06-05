import pandas as pd
import numpy as np
import sys

# コマンドライン引数の取得
seed = int(sys.argv[1])
input_csv = sys.argv[2]
output_csv = sys.argv[3]

# CSVファイルの読み込み
df = pd.read_csv(input_csv)

# 1列目をシャッフル
df['id'] = df.copy().sample(frac=1, random_state=seed).reset_index(drop=True).values

# 結果を新しいCSVファイルに保存
df.to_csv(output_csv, index=False)