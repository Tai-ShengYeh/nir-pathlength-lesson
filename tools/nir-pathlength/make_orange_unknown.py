# 第 5 章 Orange 練習用的小資料：2258 nm 二階導數訊號（×10⁻³）的 4 個校正點 + 1 個假想未知樣品。
# 寫成 data/unknown_2258_orange.csv（Orange 欄名旗標格式）。
#   python tools/nir-pathlength/make_orange_unknown.py
# 校正點：data/spectra_wide_orange.csv 的 1 mm、重複 1 光譜，Savitzky-Golay 9/2 二階導數（同第 1、5 章），
#         取 2258.067 nm，乘以 10³（與 site/nir-pathlength/lib/ch05data.js 的 signal 相同）。
# 未知樣品：訊號 y0 = 4.20（×10⁻³），與 make_data_ch05.py 的 Y0 相同；acid_pct 留空（未知）。
# 為什麼不只放一列未知樣品、欄名取成 2258.067：Orange 會把模型的 2258.067 視為「由整條光譜
# 經 SG 前處理算出來」的變數，遇到同名的普通欄位時會重新從（不存在的）光譜計算，得到 0，
# 預測值變成截距、而且不報錯。所以校正點與未知樣品放在同一個檔案、同一個普通欄位裡。
import os, csv
import numpy as np, pandas as pd
from scipy.signal import savgol_filter

here = os.path.dirname(os.path.abspath(__file__))
root = os.path.join(here, '..', '..')
Y0 = 4.20
COL = 'd2_2258_x1e3'

df = pd.read_csv(os.path.join(root, 'data', 'spectra_wide_orange.csv'))
wc = list(df.columns[4:])
s = df[(df['mC#path_mm'] == 1.0) & (df['mC#replicate'] == 1)].sort_values('cC#acid_pct')
D = savgol_filter(s[wc].values, 9, 2, deriv=2, axis=1)
sig = D[:, wc.index('2258.067')] * 1e3

out = os.path.join(root, 'data', 'unknown_2258_orange.csv')
with open(out, 'w', newline='', encoding='utf-8') as fh:
    wr = csv.writer(fh)
    wr.writerow(['mS#set', 'mS#sample', 'cC#acid_pct', COL])
    for name, c, y in zip(s['mS#sample'], s['cC#acid_pct'], sig):
        wr.writerow(['calibration', name, '%g' % c, '%.6f' % y])
    wr.writerow(['unknown', 'unknown', '', '%.2f' % Y0])
print('wrote', os.path.normpath(out))
print(open(out, encoding='utf-8').read())

# 交叉檢查：與 ch05data.js 的 signal（4 位小數）一致
ref = [8.5087, 6.3073, 5.1725, 0.9509]
assert np.allclose(np.round(sig, 4), ref), sig
print('signal matches ch05data.js:', ref)
