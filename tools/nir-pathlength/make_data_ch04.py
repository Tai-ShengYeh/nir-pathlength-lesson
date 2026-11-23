# 第 4 章「PCA 與 PLS」用的數字：寫成 site/nir-pathlength/lib/ch04data.js，並產生 Orange 用的寬表 data/spectra_wide_orange.csv
# 設定全部事先固定（不依 RMSECV 挑）：
#   二階導數 = Savitzky–Golay 9 點、二次（同第 1 章）；波段 = 報告 1 mm 可用波段 1480–1871 與 2116–2284 nm；
#   PLS 不縮放（光譜同單位）；1 個成分（只有一種成分在變：醋酸取代水）。
import os, json, itertools
import numpy as np, pandas as pd
from scipy.signal import savgol_filter
from sklearn.decomposition import PCA
from sklearn.cross_decomposition import PLSRegression

here = os.path.dirname(os.path.abspath(__file__)); root = os.path.join(here, '..', '..')
d = pd.read_csv(os.path.join(root, 'data', 'spectra_long.csv'))

# ---- Orange 寬表：每列一條光譜，欄名為波長；表頭旗標 c=目標、m=附加資訊（C 連續、S 文字）
w = d.pivot_table(index=['path_mm', 'sample', 'acid_pct', 'replicate'], columns='wavelength_nm', values='absorbance').reset_index()
w.columns = [f'{c:.3f}' if isinstance(c, float) else c for c in w.columns]
w = w.rename(columns={'path_mm': 'mC#path_mm', 'sample': 'mS#sample', 'acid_pct': 'cC#acid_pct', 'replicate': 'mC#replicate'})
w.to_csv(os.path.join(root, 'data', 'spectra_wide_orange.csv'), index=False)

d = d[d.replicate == 1]
g = d.pivot_table(index=['path_mm', 'acid_pct'], columns='wavelength_nm', values='absorbance').sort_index()
wl = g.columns.values; X = g.values; idx = g.index.to_frame().values
path, conc = idx[:, 0], idx[:, 1]

# ---- 1. PCA：20 條原始光譜（5 光徑 × 4 樣品），只置中、不標準化
p20 = PCA(2).fit(X); s20 = p20.transform(X)
r = lambda a, b: float(np.corrcoef(a, b)[0, 1])

# ---- 2. 1 mm：二階導數 + 兩段可用波段
m1 = (path == 1.0); y = conc[m1]; X1 = X[m1]
D1 = savgol_filter(X1, 9, 2, deriv=2, axis=1)
band = ((wl >= 1480) & (wl <= 1872)) | ((wl >= 2116) & (wl <= 2285))
Db, Rb = D1[:, band], X1[:, band]
p1 = PCA(2).fit(Db); s1 = p1.transform(Db)

def loo(M, yy, nc):
    e = []
    for k in range(4):
        tr = [i for i in range(4) if i != k]
        e.append(float(PLSRegression(nc, scale=False).fit(M[tr], yy[tr]).predict(M[[k]]).ravel()[0] - yy[k]))
    return float(np.sqrt(np.mean(np.square(e)))), e

cv1, e1 = loo(Db, y, 1); cv2, e2 = loo(Db, y, 2); cvr, er = loo(Rb, y, 1)
cvs, _ = loo(Db, y, 1) if False else (None, None)
from sklearn.preprocessing import StandardScaler
def loo_scaled(M, yy, nc):                                  # Orange PLS 小工具的預設：縮放特徵與目標
    e = []
    for k in range(4):
        tr = [i for i in range(4) if i != k]
        e.append(float(PLSRegression(nc, scale=True).fit(M[tr], yy[tr]).predict(M[[k]]).ravel()[0] - yy[k]))
    return float(np.sqrt(np.mean(np.square(e))))
cv_sc = loo_scaled(Db, y, 1)

# 打亂標籤（第 3 章的檢查）：24 種排列的 PLS-1 RMSECV
perms = sorted(set(itertools.permutations(range(4))))
sh = [{'perm': list(pp), 'true': list(pp) == [0, 1, 2, 3], 'rmsecv': round(loo(Db, y[list(pp)], 1)[0], 3)} for pp in perms]
wrong_better = sum(s['rmsecv'] <= cv1 for s in sh if not s['true'])

# 迴歸係數（全部 4 點、1 個成分）
pl = PLSRegression(1, scale=False).fit(Db, y); coef = pl.coef_.ravel()
wlb = wl[band]
top = np.argsort(-np.abs(coef))[:5]
fitR2 = {nc: float(PLSRegression(nc, scale=False).fit(Db, y).score(Db, y)) for nc in (1, 2, 3)}

out = {
    'pca20': {'evr': [round(float(v), 3) for v in p20.explained_variance_ratio_],
              'pc1': [round(float(v), 3) for v in s20[:, 0]], 'pc2': [round(float(v), 3) for v in s20[:, 1]],
              'path': [float(v) for v in path], 'conc': [float(v) for v in conc],
              'r_pc1_path': round(abs(r(s20[:, 0], path)), 3), 'r_pc1_conc': round(abs(r(s20[:, 0], conc)), 3)},
    'band_wl': [round(float(v), 1) for v in wlb], 'n_band': int(band.sum()),
    'pca1': {'evr': [round(float(v), 3) for v in p1.explained_variance_ratio_], 'pc1': [round(float(v), 4) for v in s1[:, 0]],
             'r_pc1_conc': round(abs(r(s1[:, 0], y)), 4)},
    'conc': [float(v) for v in y],
    'rmsecv': {'pls1': round(cv1, 3), 'pls2': round(cv2, 3), 'raw1': round(cvr, 3), 'scaled1': round(cv_sc, 3), 'band2258': 2.182},
    'err': {'pls1': [round(v, 3) for v in e1], 'pls2': [round(v, 3) for v in e2], 'raw1': [round(v, 3) for v in er]},
    'fitR2': {str(k): round(v, 4) for k, v in fitR2.items()},
    'coef': [round(float(v), 3) for v in coef], 'coef_top_wl': [round(float(wlb[i]), 1) for i in top],
    'shuffle': sh, 'wrong_better': wrong_better,
}
js = '// 由 tools/nir-pathlength/make_data_ch04.py 產生（PCA 與 PLS）\nconst CH4 = ' + json.dumps(out, ensure_ascii=False) + ';\n'
open(os.path.join(root, 'site', 'nir-pathlength', 'lib', 'ch04data.js'), 'w', encoding='utf-8').write(js)
for k in ['pca1', 'n_band', 'rmsecv', 'err', 'fitR2', 'coef_top_wl', 'wrong_better']: print(k, out[k])
print('pca20', {k: out['pca20'][k] for k in ['evr', 'r_pc1_path', 'r_pc1_conc']})
print('shuffle', sorted(s['rmsecv'] for s in sh))
# 係數集中在哪？各區段 |coef| 占比
for lo, hi in [(1480, 1700), (1700, 1840), (1840, 1872), (2116, 2200), (2200, 2285)]:
    mm = (wlb >= lo) & (wlb <= hi); print(lo, hi, round(float(np.abs(coef[mm]).sum() / np.abs(coef).sum()), 3))
