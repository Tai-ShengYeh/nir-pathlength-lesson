# -*- coding: utf-8 -*-
# 第 4 章：由「同一份」cell 定義生成本機版與 Colab 版筆記本。用法：python build_notebooks_ch04.py
import base64, gzip, os
import nbformat as nbf

HERE = os.path.dirname(os.path.abspath(__file__))
DATA_CSV = os.path.join(HERE, '..', 'data', 'spectra_long.csv')
MD, CODE = 'md', 'code'
COMMON = []


def md(s):
    COMMON.append((MD, s.strip('\n')))


def code(s):
    COMMON.append((CODE, s.strip('\n')))


md(r'''
# 第 4 章 PCA 與 PLS：程式練習

**目標**：用同一組只有 4 個樣品的光譜，先用 **PCA**「看」資料，再用只有 1 個成分的 **PLS** 做校正，並檢查它是不是真的有預測力。

| 節 | 內容 |
|---|---|
| 0 | 載入資料、二階導數、波段；PCA 與 PLS 工具函式 |
| 1 | PCA：20 條原始光譜（5 光徑 × 4 樣品），PC1 其實是「光徑」 |
| 2 | PCA：1 mm 二階導數，PC1 是「濃度」 |
| 3 | PLS 留一交叉驗證（LOOCV）：1 / 2 個成分、原始光譜、縮放、與單一波段比較 |
| 4 | 洗牌檢驗：把濃度亂排 24 種 |
| 5 | PLS-1 迴歸係數：模型看了哪些波長？ |
| 6 | 學生練習：換設定比較；Orange 版指引 |

只用 numpy、pandas、matplotlib、scipy（scipy 只用 `savgol_filter` 當二階導數）。PCA 用 SVD、PLS-1 用手寫的 NIPALS，不需要 scikit-learn。
圖上文字一律用英文（Colab 沒有中文字型）。每節結尾有「動手改改看」，並附該節的**自我驗證 PASS/FAIL**。

**所有設定都是事先固定的，不是看 RMSECV 調出來的**：二階導數 SG 9 點二次；波段 1480–1872 與 2116–2285 nm（73 個波長）；
PLS 不縮放；1 個成分。
''')

md(r'''
## 0. 載入資料、二階導數、波段

- **資料**：重複 1 的全部 5 個光徑 × 4 個樣品 = 20 條光譜；濃度用 Atago 11/02 讀值：0、3.32、4.87、7.45 %。
- **二階導數**：Savitzky–Golay 9 點、二次多項式、點距為 1（同第 1、3 章）。
- **波段**：報告中 1 mm 可用的兩段 1480–1872 與 2116–2285 nm，共 73 個波長。這是看圖與化學決定的，不是挑 RMSECV 最低的。
- **PLS 不縮放**：所有波長都是吸光度、同單位；縮放會把雜訊多的波長放大到與訊號一樣重要。
- **1 個成分**：只有一種成分在變（醋酸取代水）。

**PCA**：只置中、不標準化，用 SVD。PC 的正負號是任意的，所以和別的變數比相關時取絕對值。
**PLS1（NIPALS）**：X、y 先置中；每個成分 w = X′y（正規化）、t = Xw、p = X′t / t′t、q = y′t / t′t，然後從 X、y 扣掉 t p′、t q。
''')

code(r'''
import itertools
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter

np.set_printoptions(precision=5, linewidth=130)
plt.rcParams['figure.dpi'] = 90

# ---- 事先固定的設定（先決定、不依 RMSECV 挑選）----
SG_WIN = 9                                  # 二階導數：Savitzky-Golay 9 點、二次多項式
BANDS = [(1480, 1872), (2116, 2285)]        # 兩段波段（nm），共 73 個波長
N_COMP = 1                                  # PLS 成分數
SCALE = False                               # PLS 不縮放

# ---- 自我驗證的小工具：每節結尾呼叫 show_checks() ----
checks = []
_shown = 0
def chk(name, got, key, tol):
    checks.append((bool(abs(got - key) <= tol), name, float(got), key))

def show_checks(title):
    global _shown
    rows = checks[_shown:]
    print('\n[自我驗證] %s' % title)
    print('%-6s | %-52s | %11s | %11s' % ('result', 'item', 'yours', 'key')); print('-' * 88)
    for ok, name, got, key in rows:
        print('%-6s | %-52s | %11.5g | %11.5g' % ('PASS' if ok else 'FAIL', name, got, key))
    print('本節：%d / %d PASS' % (sum(r[0] for r in rows), len(rows)))
    _shown = len(checks)
''')

code(r'''
df = pd.read_csv(DATA_PATH)
d_r1 = df[df['replicate'] == 1]
g20 = d_r1.pivot_table(index=['path_mm', 'acid_pct'], columns='wavelength_nm', values='absorbance').sort_index()
wl = g20.columns.values; n_wl = len(wl)
X20 = g20.values                                                    # 20 x 142，列依「光徑、濃度」排序
path20 = g20.index.get_level_values(0).values; conc20 = g20.index.get_level_values(1).values
paths = np.unique(path20)
m1 = np.isclose(path20, 1.0)
X1 = X20[m1]; conc = conc20[m1]                                     # 1 mm：4 x 142；conc = 0, 3.32, 4.87, 7.45
print('光徑（mm）：', paths, '\n濃度（Atago 11/02, %）：', conc)
print('20 條光譜矩陣：', X20.shape, '；1 mm 光譜矩陣：', X1.shape)

band = np.zeros(n_wl, bool)
for lo, hi in BANDS:
    band |= (wl >= lo) & (wl <= hi)
wlb = wl[band]
D1 = savgol_filter(X1, SG_WIN, 2, deriv=2, axis=1)                  # 二階導數 4 x 142
Db = D1[:, band]                                                    # 4 x 73：二階導數、只取兩段波段
Rb = X1[:, band]                                                    # 4 x 73：原始吸光度、同樣兩段波段（比較用）
print('波段內波長數：%d（%s）' % (band.sum(), ' 與 '.join('%g-%g nm' % b for b in BANDS)))
''')

code(r'''
# ---- PCA：只置中、不標準化（用 SVD）----
def pca_fit(M):
    Mc = M - M.mean(axis=0)
    U, S, Vt = np.linalg.svd(Mc, full_matrices=False)
    return {'scores': U * S, 'evr': S ** 2 / np.sum(S ** 2), 'loadings': Vt}

def abs_r(a, b):
    """PC 的正負號是任意的，所以取絕對值。"""
    return float(abs(np.corrcoef(a, b)[0, 1]))

# ---- PLS1（NIPALS，X 與 y 置中、不縮放）----
def pls1_fit(X, y, nc):
    xm = X.mean(axis=0); ym = y.mean()
    Xc = X - xm; yc = y - ym
    W = np.zeros((X.shape[1], nc)); P = np.zeros((X.shape[1], nc)); Q = np.zeros(nc)
    for a in range(nc):
        w = Xc.T @ yc; w = w / np.sqrt(w @ w)
        t = Xc @ w; tt = t @ t
        p = Xc.T @ t / tt; q = (yc @ t) / tt
        Xc = Xc - np.outer(t, p); yc = yc - t * q                    # 扣除（deflation）
        W[:, a] = w; P[:, a] = p; Q[a] = q
    B = W @ np.linalg.solve(P.T @ W, Q)                               # 迴歸係數（對置中後的 X）
    return {'B': B, 'b0': ym - xm @ B}

def pls1_predict(fit, Xnew):
    return np.atleast_2d(Xnew) @ fit['B'] + fit['b0']

def rmse(e):
    return float(np.sqrt(np.mean(np.square(e))))

def pls_loocv(M, y, nc=N_COMP, scale=SCALE):
    """留一交叉驗證：回傳 (RMSECV, 誤差 = 預測 - 真值)。scale=True 時每折用訓練折的平均與標準差標準化 X 與 y。"""
    n = len(y); e = []
    for k in range(n):
        tr = [i for i in range(n) if i != k]
        if scale:
            mu = M[tr].mean(axis=0); sd = M[tr].std(axis=0, ddof=1); ym = y[tr].mean(); ys = y[tr].std(ddof=1)
            fit = pls1_fit((M[tr] - mu) / sd, (y[tr] - ym) / ys, nc)
            pr = pls1_predict(fit, (M[[k]] - mu) / sd)[0] * ys + ym
        else:
            pr = pls1_predict(pls1_fit(M[tr], y[tr], nc), M[[k]])[0]
        e.append(float(pr - y[k]))
    return rmse(e), np.array(e)

chk('number of wavelengths', n_wl, 142, 0)
chk('number of spectra (5 paths x 4 samples)', X20.shape[0], 20, 0)
chk('wavelengths in the two bands', band.sum(), 73, 0)
show_checks('第 0 節')
''')

md(r'''
### 動手改改看（第 0 節）
1. 為什麼 PLS 不縮放、只取 1 個成分、波段取 1480–1872 與 2116–2285？這些設定是在看到 RMSECV **之前**就決定的，為什麼這點重要？
2. 把 `BANDS` 改成只留第一段，波段內波長數是多少？（改完記得改回來再看自我驗證）
''')

md(r'''
## 1. PCA：20 條原始光譜

5 個光徑 × 4 個樣品 = 20 條光譜，只置中、不標準化，看第 1、2 個主成分（PC）解釋多少變異，又和什麼有關。
''')

code(r'''
p20 = pca_fit(X20); s20 = p20['scores']; ev20 = p20['evr']
r_path = abs_r(s20[:, 0], path20); r_conc = abs_r(s20[:, 0], conc20)
print('解釋變異：PC1 = %.3f，PC2 = %.3f' % (ev20[0], ev20[1]))
print('|r(PC1, 光徑)| = %.3f；|r(PC1, 濃度)| = %.3f' % (r_path, r_conc))
print('PC1 幾乎就是「光徑」：光徑把整條光譜等比例放大，這個「大小」變化蓋過了濃度造成的細微差異。')

fig, ax = plt.subplots(figsize=(8, 4.8))
sc = ax.scatter(s20[:, 0], s20[:, 1], c=path20, cmap='viridis', s=60)
fig.colorbar(sc, label='path length (mm)')
ax.set_xlabel('PC1 (%.1f%%)' % (100 * ev20[0])); ax.set_ylabel('PC2 (%.1f%%)' % (100 * ev20[1]))
ax.set_title('PCA of 20 raw spectra (colored by path length)')
plt.tight_layout(); plt.show()

chk('PCA 20 raw: explained variance PC1', ev20[0], 0.917, 0.001)
chk('PCA 20 raw: explained variance PC2', ev20[1], 0.069, 0.001)
chk('PCA 20 raw: |r(PC1, path)|', r_path, 0.91, 0.002)
chk('PCA 20 raw: |r(PC1, conc)|', r_conc, 0.022, 0.002)
show_checks('第 1 節')
''')

md(r'''
### 動手改改看（第 1 節）
1. 散佈圖改成「依濃度上色」（`c=conc20`），看到什麼？濃度的資訊藏在哪個 PC？（試試 PC2、PC3）
2. 為什麼 PC1 的正負號每次可能不同，卻不影響結論？
3. PCA 是「描述」資料的工具（沒用到濃度）：它能告訴你最大的變異來自哪裡，卻不能直接預測濃度。
''')

md(r'''
## 2. PCA：1 mm、二階導數、兩段波段

換成 1 mm 的 4 條二階導數光譜（4 × 73）。二階導數去掉了基線與光徑造成的緩慢變化，只剩吸收峰的形狀。
''')

code(r'''
p1 = pca_fit(Db); s1 = p1['scores']; ev1 = p1['evr']
r_conc1 = abs_r(s1[:, 0], conc)
print('解釋變異：PC1 = %.3f，PC2 = %.3f' % (ev1[0], ev1[1]))
print('|r(PC1, 濃度)| = %.4f' % r_conc1)
print('二階導數之後，PC1 幾乎就是濃度。')

fig, ax = plt.subplots(figsize=(7, 4.2))
ax.plot(conc, s1[:, 0], 'o', color='steelblue', ms=9)
for c, s in zip(conc, s1[:, 0]):
    ax.annotate('%g' % c, (c, s), textcoords='offset points', xytext=(0, 8), ha='center', fontsize=8)
b, a = np.polyfit(conc, s1[:, 0], 1); xx = np.linspace(conc.min(), conc.max(), 20)
ax.plot(xx, a + b * xx, 'k--', lw=1)
ax.set_xlabel('acid (%)'); ax.set_ylabel('PC1 score (%.1f%%)' % (100 * ev1[0]))
ax.set_title('PC1 of 2nd-derivative spectra vs concentration (1 mm)')
plt.tight_layout(); plt.show()

chk('PCA 1 mm d2: explained variance PC1', ev1[0], 0.969, 0.001)
chk('PCA 1 mm d2: |r(PC1, conc)|', r_conc1, 0.9984, 0.002)
show_checks('第 2 節')
''')

md(r'''
### 動手改改看（第 2 節）
1. 對 `Rb`（原始吸光度、不微分）做 PCA，PC1 的解釋變異與 |r(PC1, conc)| 變成多少？
2. 為什麼只有 4 個樣品時，PCA 最多只有 3 個有意義的 PC？（看 `ev1`）
''')

md(r'''
## 3. PLS 留一交叉驗證（LOOCV）

每一折：留下 1 個樣品，用其餘 3 個擬合 PLS，預測被留下的樣品；誤差 = 預測 − 真值（%），**RMSECV = √(mean(err²))**。
同時比較：2 個成分、原始光譜（不微分）、縮放 X 與 y（Orange 的 PLS 小工具預設）、第 1 章的單一波段 2258 nm（RMSECV = 2.182 %）。
''')

code(r'''
cv1, e1 = pls_loocv(Db, conc, 1, False)
cv2, e2 = pls_loocv(Db, conc, 2, False)
cvr, er = pls_loocv(Rb, conc, 1, False)                             # 原始光譜（不微分）、同樣波段
cvs, es = pls_loocv(Db, conc, 1, True)                              # 縮放 X 與 y
RM2258 = 2.182                                                      # 第 1 章：固定單一波段 2258 nm 的 RMSECV

print(pd.DataFrame({'left_out': conc, 'PLS1_err': e1, 'PLS2_err': e2, 'raw_PLS1_err': er}).to_string(index=False, float_format=lambda v: '%.3f' % v))
print('RMSECV（%%）：PLS-1 = %.3f；PLS-2 = %.3f；原始光譜 PLS-1 = %.3f；縮放 PLS-1 = %.3f；單一波段 2258 nm（第 1 章）= %.3f' %
      (cv1, cv2, cvr, cvs, RM2258))
print('1 個成分最好：只有一種成分在變，多加的成分只是在擬合雜訊。二階導數把 RMSECV 從 2.5 降到 0.3。')

def fit_r2(nc):
    pr = pls1_predict(pls1_fit(Db, conc, nc), Db)
    return float(1 - np.sum((conc - pr) ** 2) / np.sum((conc - conc.mean()) ** 2))
fitR2 = [fit_r2(nc) for nc in (1, 2, 3)]
print('全部 4 點擬合的 R2：1 個成分 = %.4f；2 個 = %.4f；3 個 = %.4f' % tuple(fitR2))
print('3 個成分的 R2 = 1.0000：4 個樣品、3 個成分，模型可以把 4 點完全吃掉。擬合 R2 愈高不代表愈好（過度擬合），要看 RMSECV。')

fig, ax = plt.subplots(figsize=(9, 4))
xx = np.arange(4)
ax.bar(xx - 0.2, e1, 0.4, color='red', label='1 component (%.3f)' % cv1)
ax.bar(xx + 0.2, e2, 0.4, color='steelblue', label='2 components (%.3f)' % cv2)
ax.axhline(0, color='k'); ax.set_xticks(xx); ax.set_xticklabels(['leave out %g' % c for c in conc])
ax.set_ylabel('prediction - true (% acid)'); ax.set_title('LOOCV errors of PLS (2nd derivative, 73 wavelengths)'); ax.legend(loc='upper left', fontsize=8)
plt.tight_layout(); plt.show()

chk('RMSECV PLS-1 (%)', cv1, 0.317, 0.002)
chk('PLS-1 errors: max abs diff to key', np.abs(e1 - np.array([0.46, 0.252, -0.341, -0.104])).max(), 0, 0.002)
chk('RMSECV PLS-2 (%)', cv2, 0.461, 0.002)
chk('PLS-2 errors: max abs diff to key', np.abs(e2 - np.array([-0.209, 0.43, -0.406, 0.676])).max(), 0, 0.002)
chk('RMSECV raw spectra PLS-1 (%)', cvr, 2.547, 0.002)
chk('RMSECV scaled PLS-1 (%)', cvs, 2.458, 0.002)
chk('fit R2, 1 component', fitR2[0], 0.9968, 0.002)
chk('fit R2, 2 components', fitR2[1], 0.9986, 0.002)
chk('fit R2, 3 components', fitR2[2], 1.0, 0.002)
show_checks('第 3 節')
''')

md(r'''
### 動手改改看（第 3 節）
1. 為什麼 2 個成分的 RMSECV（0.461）比 1 個（0.317）差，擬合 R² 卻更高？
2. 原始光譜（2.547）、縮放（2.458）與單一波段 2258 nm（2.182）都遠不如二階導數 + 不縮放（0.317）。各自的原因是什麼？
3. 4 個樣品置中後只有 3 個自由度，為什麼 3 個成分一定 R² = 1？
''')

md(r'''
## 4. 洗牌檢驗：把濃度亂排，PLS 還有沒有預測力？

把 4 個濃度的順序打亂（4! = 24 種排列，只有 1 種是真的），對每種排列做同一套 PLS-1 LOOCV（波段、成分數都不變）。
若 0.317 % 真的反映光譜與濃度的關係，亂排後的 RMSECV 應該都比較差。
''')

code(r'''
perms = sorted(itertools.permutations(range(4)))                    # 24 種，第 1 個是原本的順序
sh = pd.DataFrame({'perm': [''.join(str(i + 1) for i in p) for p in perms],
                   'true': [p == (0, 1, 2, 3) for p in perms],
                   'rmsecv': [pls_loocv(Db, conc[list(p)], 1, False)[0] for p in perms]})
wrong = ~sh['true']
n_wrong_better = int((sh.loc[wrong, 'rmsecv'] <= cv1).sum())
srt = np.sort(sh['rmsecv'].values)
print(sh.sort_values('rmsecv').head(5).to_string(index=False, float_format=lambda v: '%.4f' % v))
print('真正順序的 RMSECV = %.3f，在 24 種裡排第 %d（最低）；錯誤排列中 RMSECV <= %.3f 的有 %d 種；第二低 = %.3f' %
      (sh.loc[sh['true'], 'rmsecv'].iloc[0], int(sh['rmsecv'].rank().loc[sh['true']].iloc[0]), cv1, n_wrong_better, srt[1]))
print('對比第 3 章：挑波長時亂排的濃度仍能湊出高 R2；這裡波段與成分數事先固定，亂排就湊不出好的預測。')

s = sh.sort_values('rmsecv')
fig, ax = plt.subplots(figsize=(9, 4))
ax.bar(range(24), s['rmsecv'], color=np.where(s['true'], 'red', 'lightgray'))
ax.set_xlabel('24 orderings (sorted; red = true order)'); ax.set_ylabel('RMSECV (% acid)')
ax.set_title('PLS-1 RMSECV for all 24 concentration orderings')
plt.tight_layout(); plt.show()

chk('shuffle: number of orderings', len(sh), 24, 0)
chk('shuffle: true ordering is the lowest (1 = yes)', float(sh['rmsecv'].idxmin() == sh.index[sh['true']][0]), 1, 0)
chk('shuffle: true ordering RMSECV', sh.loc[sh['true'], 'rmsecv'].iloc[0], 0.317, 0.002)
chk('shuffle: wrong orderings with RMSECV <= 0.317', n_wrong_better, 0, 0)
chk('shuffle: second lowest RMSECV', srt[1], 0.921, 0.002)
show_checks('第 4 節')
''')

md(r'''
### 動手改改看（第 4 節）
1. 把 `pls_loocv(..., 1, False)` 的成分數改成 2，真正順序還是最低嗎？
2. 這個檢驗的邏輯與第 3 章洗牌檢驗相同，結論有什麼不同？為什麼？（提示：本章沒有「從 73 個波長挑最大」這一步）
''')

md(r'''
## 5. PLS-1 迴歸係數：模型看了哪些波長？

用全部 4 個樣品、1 個成分擬合，看迴歸係數對波長的分布。
''')

code(r'''
fit_all = pls1_fit(Db, conc, 1)
coef = fit_all['B']
top = np.argsort(-np.abs(coef))[:5]
print('|係數| 最大的 5 個波長（nm）：', np.round(wlb[top], 1))
def share(lo, hi):
    m = (wlb >= lo) & (wlb <= hi)
    return float(np.abs(coef[m]).sum() / np.abs(coef).sum())
sh_a = share(1840, 1872); sh_b = share(2200, 2285)
print('|係數| 占比：1840-1872 nm = %.3f（水的 1900 nm 帶的短波側）；2200-2285 nm = %.3f（含 2258 nm 的醋酸吸收帶）' % (sh_a, sh_b))
print('模型自己找到了與化學相符的區域（不是我們告訴它的），而且用了很多波長一起平均雜訊。')

fig, ax = plt.subplots(figsize=(10, 4))
ax.vlines(wlb, 0, coef, color='steelblue', lw=1.5, label='coefficient')
ax.plot(wlb[top], coef[top], 'r^', ms=8, label='top 5 |coef|')
ax.axhline(0, color='k', lw=0.8)
ax.set_xlabel('wavelength (nm)'); ax.set_ylabel('PLS-1 regression coefficient'); ax.set_title('PLS-1 coefficients (2nd derivative, 73 wavelengths)')
ax.legend(loc='upper left')
plt.tight_layout(); plt.show()

key_top = np.array([1862.4, 1871.0, 1853.9, 2258.1, 2245.5])
chk('top-5 |coef| wavelengths match (of 5)', (np.abs(wlb[top] - key_top) < 0.05).sum(), 5, 0)
chk('top-1 wavelength (nm)', wlb[top[0]], 1862.4, 0.05)
chk('top-2 wavelength (nm)', wlb[top[1]], 1871.0, 0.05)
chk('top-3 wavelength (nm)', wlb[top[2]], 1853.9, 0.05)
chk('share of |coef| in 1840-1872 nm', sh_a, 0.391, 0.002)
chk('share of |coef| in 2200-2285 nm', sh_b, 0.276, 0.002)
show_checks('第 5 節')
''')

md(r'''
### 動手改改看（第 5 節）
1. 係數有正有負、形狀像二階導數的吸收峰。為什麼？（提示：二階導數把吸收峰變成「負的尖峰」）
2. 只留 |係數| 最大的 5 個波長重做 LOOCV（這是在迴圈外挑波長！），RMSECV 會變多少？這犯了第 3 章的哪個錯？
''')

md(r'''
## 6. 學生練習：換設定、比較

本章的設定（二階導數、兩段波段、不縮放、1 個成分）是事先固定的。請你自己換設定，看 RMSECV 怎麼變。
Orange 版的操作說明在 `exercises/orange/ch04_orange_guide.md`。
''')

code(r'''
def try_variant(name, M, nc, scale):
    print('%-34s nc = %d, scale = %-5s -> RMSECV = %.3f %%' % (name, nc, scale, pls_loocv(M, conc, nc, scale)[0]))

# 請改下面的設定再執行：
try_variant('baseline (d2, no scale, 1 comp)', Db, 1, False)
try_variant('scale = True', Db, 1, True)
try_variant('raw spectra (no derivative)', Rb, 1, False)
try_variant('2 components', Db, 2, False)
try_variant('raw + scale = True', Rb, 1, True)
try_variant('scale = True, 2 components', Db, 2, True)

labels = ['d2 + 1 comp\n(fixed)', 'd2 + 2 comp', 'd2 + scale', 'raw + 1 comp', 'single band\n2258 nm']
vals = [cv1, cv2, cvs, cvr, RM2258]
fig, ax = plt.subplots(figsize=(9, 4.2))
ax.bar(labels, vals, color=['red', 'steelblue', 'steelblue', 'steelblue', 'gray'])
ax.set_ylabel('RMSECV (% acid)'); ax.set_title('RMSECV of different settings (n = 4, LOOCV)')
plt.tight_layout(); plt.show()
''')

md(r'''
### 動手改改看（第 6 節）
1. `scale = True` 為什麼變差？（提示：縮放把每個波長的變異都放成 1，雜訊多的波長被放大到與訊號一樣重要）
2. 為什麼不能把 `nc` 設成 3 再做 LOOCV？（提示：每折只剩 3 個訓練樣品，置中後只剩 2 個自由度；試著執行看看會怎樣）
3. 在這 4 個樣品上，你「試了很多組設定、挑 RMSECV 最低的」報出來，會有什麼問題？（回想第 3 章）
4. 到 Orange 照 `exercises/orange/ch04_orange_guide.md` 的步驟重做，結果應與本頁一致（PLS 小工具的預設是縮放，會是 2.458）。
''')

code(r'''
print('全部自我驗證：%d / %d PASS' % (sum(c[0] for c in checks), len(checks)))
for c in checks:
    if not c[0]:
        print('FAIL:', c[1], c[2], c[3])
''')

md(r'''
### 總結
1. PCA 描述資料：原始光譜的 PC1 = 光徑（0.917），二階導數的 PC1 = 濃度（0.969）。
2. PLS 1 個成分、二階導數、兩段波段、不縮放：RMSECV = 0.317 %；這些設定都是事先固定的。
3. 洗牌檢驗：24 種排列只有真正順序最低，沒有錯誤排列比它好。
4. 係數集中在水與醋酸的吸收帶，與化學相符。
5. 只有 4 個樣品，0.317 % 仍是「初步」結果：需要更多獨立樣品驗證。
''')


def local_head():
    return [
        (MD, '> **本機 Jupyter 版**：請在 `teaching/exercises/` 目錄下開啟，資料讀自 `../data/spectra_long.csv`。\n'
             '> 需要套件：numpy、pandas、matplotlib、scipy。Colab 請改開 `ch04_pca_pls_colab.ipynb`。'),
        (CODE, "DATA_PATH = '../data/spectra_long.csv'   # 本機版：讀上一層 data 資料夾"),
    ]


def colab_head():
    with open(DATA_CSV, 'rb') as fh:
        raw = fh.read()
    b64 = base64.b64encode(gzip.compress(raw, 9, mtime=0)).decode('ascii')
    wrapped = '\n'.join(b64[i:i + 100] for i in range(0, len(b64), 100))
    embed = (
        "# 資料內嵌：把 spectra_long.csv 以 gzip + base64 存成字串，執行時解出成檔案\n"
        "import base64, gzip\n"
        "DATA_B64 = \"\"\"\n" + wrapped + "\n\"\"\"\n"
        "with open('spectra_long.csv', 'wb') as _f:\n"
        "    _f.write(gzip.decompress(base64.b64decode(DATA_B64)))\n"
        "DATA_PATH = 'spectra_long.csv'\n"
        "print('已解出 spectra_long.csv：', len(gzip.decompress(base64.b64decode(DATA_B64))), 'bytes')"
    )
    upload = (
        "# 【可選】想改用自己上傳的 spectra_long.csv，執行這一格並選檔；非 Colab 環境會自動略過。\n"
        "try:\n"
        "    from google.colab import files\n"
        "    _up = files.upload()\n"
        "    if 'spectra_long.csv' in _up:\n"
        "        DATA_PATH = 'spectra_long.csv'\n"
        "        print('已改用你上傳的 spectra_long.csv')\n"
        "    else:\n"
        "        print('沒有上傳 spectra_long.csv，繼續使用內嵌資料')\n"
        "except Exception as _e:\n"
        "    print('略過上傳（非 Colab 環境或未選檔）：', type(_e).__name__)"
    )
    return [
        (MD, '> **Google Colab 版**\n>\n'
             '> 1. 開啟 [colab.research.google.com](https://colab.research.google.com)。\n'
             '> 2. 選單「**檔案 → 上傳筆記本**」，選擇 `ch04_pca_pls_colab.ipynb`。\n'
             '> 3. 依序執行每一格。\n>\n'
             '> 資料已**內嵌**在筆記本中，下一格會自動解出 `spectra_long.csv`。Colab 已內建所需套件。'),
        (CODE, embed),
        (MD, '若想改用自己上傳的檔案，可執行下面這一格（可選）。'),
        (CODE, upload),
    ]


def build(name, head):
    nb = nbf.v4.new_notebook()
    cells = [COMMON[0]] + head + COMMON[1:]
    nb['cells'] = [nbf.v4.new_markdown_cell(s) if k == MD else nbf.v4.new_code_cell(s) for k, s in cells]
    nb['metadata'] = {'kernelspec': {'display_name': 'Python 3', 'language': 'python', 'name': 'python3'},
                      'language_info': {'name': 'python'}}
    path = os.path.join(HERE, name)
    with open(path, 'w', encoding='utf-8') as fh:
        nbf.write(nb, fh)
    print('wrote', path, '(%d cells)' % len(cells))


if __name__ == '__main__':
    build('ch04_pca_pls.ipynb', local_head())
    build('ch04_pca_pls_colab.ipynb', colab_head())
