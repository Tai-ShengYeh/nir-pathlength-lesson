# -*- coding: utf-8 -*-
# 第 3 章：由「同一份」cell 定義生成本機版與 Colab 版筆記本。用法：python build_notebooks_ch03.py
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
# 第 3 章 波長挑選的陷阱：程式練習

**目標**：只有 4 個樣品時，從 142 個波長挑 R² 最高的，會發生什麼事？本章用第 1、2 章同一組資料（1 mm、重複 1）：

| 節 | 內容 |
|---|---|
| 0 | 載入資料、二階導數、可用波段 |
| 1 | 每個波長的 R²，挑最高的（含飽和的陷阱） |
| 2 | 只在可用波段挑 vs 固定 2258 nm |
| 3 | 洗牌檢驗：把濃度亂排 24 種 |
| 4 | 留一交叉驗證（LOOCV）三種做法：挑波長在迴圈外／內、固定波長 |
| 5 | 學生練習：改變挑選池；選做：142 個隨機波長 |
| 6 | 自我驗證（PASS/FAIL 表） |

只用 numpy、pandas、matplotlib、scipy（scipy 只用 `savgol_filter` 當二階導數）。圖上文字一律用英文（Colab 沒有中文字型）。
每節結尾有「動手改改看」。
''')

md(r'''
## 0. 載入資料

1 mm 光徑、重複 1，4 個樣品（水＋三瓶醋），濃度用 Atago 11/02 讀值：0、3.32、4.87、7.45 %，共 142 個波長。

- **二階導數**：Savitzky–Golay 9 點、二次多項式、點距為 1（與第 1、2 章相同）。
- **可用波段**：4 條光譜吸光度**都小於 1.5** 的波長。吸光度太高時偵測器接近飽和、訊號不再線性，不可信。
''')

code(r'''
import itertools
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter

np.set_printoptions(precision=5, linewidth=130)
plt.rcParams['figure.dpi'] = 90

# ---- 可調整的參數 ----
SG_WIN   = 9        # 二階導數 SG 視窗（二次多項式）
A_MAX    = 1.5      # 可用波段：4 條光譜吸光度都低於此值
FIXED_WL = 2258     # 由化學選定的固定吸收帶（nm）
N_RAND   = 1000     # 第 5 節選做：隨機資料重複次數
SEED     = 2026
''')

code(r'''
df = pd.read_csv(DATA_PATH)
d1 = df[np.isclose(df['path_mm'], 1.0) & (df['replicate'] == 1)]
g = d1.pivot_table(index='acid_pct', columns='wavelength_nm', values='absorbance').sort_index()
wl = g.columns.values; X = g.values; conc = g.index.values          # X: 4 x 142；conc = 0, 3.32, 4.87, 7.45
n_wl = len(wl)

D = savgol_filter(X, SG_WIN, 2, deriv=2, axis=1)                    # 二階導數，單位：每點^2
usable = (X < A_MAX).all(axis=0)                                    # 4 條光譜都 A < 1.5

def r2(a, b):
    return float(np.corrcoef(a, b)[0, 1] ** 2)

def r2_all(rows, yy):
    """每個波長（欄）的二階導數對 yy 的 R²。rows 選出要用的樣品列。"""
    return np.array([r2(D[rows, j], yy) for j in range(n_wl)])

def best_idx(R, ok=None):
    """在 ok 為 True 的波長中挑 R² 最高者（平手取第一個）。"""
    ok = np.ones(len(R), bool) if ok is None else ok
    return int(np.argmax(np.where(ok, R, -1.0)))

i_fixed = int(np.argmin(np.abs(wl - FIXED_WL)))
print('光譜矩陣：', X.shape, ' 濃度：', conc)
print('可用波長數（4 條光譜 A < %.1f）：%d / %d' % (A_MAX, usable.sum(), n_wl))
''')

md(r'''
### 動手改改看（第 0 節）
1. 印出 `X.max()`，哪個波長區的吸光度最高？為什麼 1 mm 光徑在那裡會「飽和」？
2. 把 `A_MAX` 改成 1.2 或 1.8，可用波長數怎麼變？
''')

md(r'''
## 1. 每個波長的 R²，挑最高的

對 142 個波長各做一次「二階導數 vs 濃度」的直線迴歸，看 R²。n = 4 時，R² 很容易高；把**全部**波長裡 R² 最高的挑出來，
看起來最漂亮——但它是不是可信？
''')

code(r'''
R = r2_all(slice(None), conc)
j_all = best_idx(R)
n_over95 = int((R > 0.95).sum())
print('142 個波長全挑：%.1f nm，R2 = %.4f，4 條光譜吸光度 = %s' % (wl[j_all], R[j_all], np.round(X[:, j_all], 3)))
print('該波長有 %d 條光譜吸光度 >= %.1f（偵測器飽和／非線性），不是可信的訊號。' % ((X[:, j_all] >= A_MAX).sum(), A_MAX))
print('R2 > 0.95 的波長數：', n_over95, '（n = 4 時很容易有一堆「看起來很好」的波長）')

fig, ax = plt.subplots(figsize=(10, 4))
ax.plot(wl, R, 'k', lw=1, label='R2')
bad = ~usable                                                       # 飽和區（灰色）
ax.fill_between(wl, 0, 1, where=bad, color='gray', alpha=0.25, step='mid', label='saturated (A >= 1.5)')
ax.plot(wl[j_all], R[j_all], 'r^', ms=10, label='best of all 142')
ax.axhline(0.95, color='gray', ls=':')
ax.set_ylim(0, 1.02); ax.set_xlabel('wavelength (nm)'); ax.set_ylabel('R2 (2nd derivative vs concentration)')
ax.set_title('R2 at every wavelength (n = 4)'); ax.legend(loc='lower left')
plt.tight_layout(); plt.show()
''')

md(r'''
### 動手改改看（第 1 節）
1. 把門檻 0.95 改成 0.99，還剩幾個波長？
2. 畫出 `j_all` 那個波長的 4 個點（`plt.plot(conc, D[:, j_all], 'o')`），它們真的沿著濃度排成一直線嗎？
3. 飽和波段（灰色）的 R² 為什麼可以很高？這是化學還是儀器？
''')

md(r'''
## 2. 只在可用波段挑，以及固定用化學帶

排除飽和波段後，再挑一次；並與「事先由化學決定」的 2258 nm 比較。
''')

code(r'''
j_use = best_idx(R, usable)
print('可用波段最佳：%.1f nm，R2 = %.4f' % (wl[j_use], R[j_use]))
print('固定化學帶：  %.1f nm，R2 = %.4f' % (wl[i_fixed], R[i_fixed]))
print('挑出來的 R2 更高，但那是「從 %d 個裡挑最大」的結果，不代表它真的更會預測。' % usable.sum())
''')

md(r'''
### 動手改改看（第 2 節）
1. 在可用波段裡，R² 第 2、3 高的是哪些波長？和最佳差多少？
2. 為什麼說 2258 nm 的 R² 比較「誠實」？（提示：它是挑之前就由化學決定的）
''')

md(r'''
## 3. 洗牌檢驗：把濃度亂排，看挑波長還能不能「找到」好 R²

把 4 個濃度的順序打亂（共 4! = 24 種排列，只有 1 種是真的）。亂排之後濃度與光譜**沒有任何真實關係**。
對每種排列，記錄：

- **best**：在可用波段挑最佳波長的 R²；
- **fixed**：固定 2258 nm 的 R²。
''')

code(r'''
perms = sorted(itertools.permutations(range(4)))                    # 24 種，第 1 個是原本的順序
rows = []
for p in perms:
    yp = conc[list(p)]
    Rp = r2_all(slice(None), yp)
    rows.append({'perm': ''.join(str(i + 1) for i in p), 'true': p == (0, 1, 2, 3),
                 'best': float(Rp[usable].max()), 'fixed': float(Rp[i_fixed])})
shuf = pd.DataFrame(rows)
wrong = ~shuf['true']
n_wrong_best = int((shuf.loc[wrong, 'best'] > 0.9).sum())
n_wrong_fixed = int((shuf.loc[wrong, 'fixed'] > 0.9).sum())
print(shuf.head(6).to_string(index=False, float_format=lambda v: '%.4f' % v))
print('錯誤排列共 %d 種；其中「可用波段最佳 R2」> 0.9 的有 %d 種；「固定 2258 nm R2」> 0.9 的有 %d 種' %
      (wrong.sum(), n_wrong_best, n_wrong_fixed))
print('濃度明明是亂排的，挑最佳波長卻幾乎總能湊出 R2 > 0.9：高 R2 在 n = 4 時沒有證明力。')

fig, ax = plt.subplots(1, 2, figsize=(10.5, 3.8))
for a, (nm, title) in zip(ax, [('best', 'best of usable wavelengths'), ('fixed', 'fixed 2258 nm')]):
    s = shuf.sort_values(nm)
    a.bar(range(24), s[nm], color=np.where(s['true'], 'red', 'lightgray'))
    a.axhline(0.9, color='k', ls='--'); a.set_ylim(0, 1.02)
    a.set_xlabel('24 orderings (sorted; red = true order)'); a.set_ylabel('R2'); a.set_title(title)
plt.tight_layout(); plt.show()
''')

md(r'''
### 動手改改看（第 3 節）
1. 紅色是真正的順序。「挑最佳波長」時，真正順序在 24 種裡排第幾？「固定 2258 nm」呢？
2. 把門檻 0.9 改成 0.95 或 0.8，兩個計數怎麼變？
3. 這個檢驗的邏輯是什麼？（亂排的濃度不可能有真實關係，卻能得到高 R² → 「高 R²」本身不能證明什麼）
''')

md(r'''
## 4. 留一交叉驗證（LOOCV）的三種做法

每一折：留下 1 個樣品，用其餘 3 點擬合 **d2 = a + b·c**，再反推被留下樣品的濃度：pred = (d2 − a) / b，
誤差 err = pred − true，**RMSECV = √(mean(err²))**。三種做法只差「波長什麼時候挑」：

| 做法 | 波長怎麼來 | 問題 |
|---|---|---|
| (a) 挑在迴圈**外** | 先用全部 4 點挑好 1862.4 nm，再做 LOOCV | 被留下的樣品已參與挑選，**資訊洩漏** |
| (b) 挑在迴圈**內** | 每折只用 3 個訓練點重新挑最佳可用波長 | **誠實**，但每折挑到的波長都不同 |
| (c) 固定 2258 nm | 挑之前就由化學決定 | 沒有挑選偏差 |
''')

code(r'''
def inv_pred(tr, k, j):
    """用 tr 列的 3 點擬合 d2 = a + b*c，反推第 k 列的濃度。"""
    b, a = np.polyfit(conc[tr], D[tr, j], 1)
    return float((D[k, j] - a) / b)

def rmse(e):
    return float(np.sqrt(np.mean(np.square(e))))

e_out, e_in, e_fix, pick_in = [], [], [], []
for k in range(4):
    tr = [i for i in range(4) if i != k]
    e_out.append(inv_pred(tr, k, j_use) - conc[k])                          # (a) 迴圈外挑的波長
    jr = best_idx(r2_all(tr, conc[tr]), usable)                             # (b) 每折只用訓練點重新挑
    pick_in.append(float(wl[jr])); e_in.append(inv_pred(tr, k, jr) - conc[k])
    e_fix.append(inv_pred(tr, k, i_fixed) - conc[k])                        # (c) 固定 2258 nm
e_out, e_in, e_fix, pick_in = map(np.array, (e_out, e_in, e_fix, pick_in))

cv = pd.DataFrame({'left_out': conc, 'outside_err': e_out, 'inside_pick_nm': pick_in, 'inside_err': e_in, 'fixed_err': e_fix})
print(cv.to_string(index=False, float_format=lambda v: '%.3f' % v))
RM = {'outside': rmse(e_out), 'inside': rmse(e_in), 'fixed': rmse(e_fix)}
print('RMSECV：挑波長在迴圈外 %.3f %%（偷看了答案）｜在迴圈內 %.3f %%（誠實）｜固定 2258 nm %.3f %%' %
      (RM['outside'], RM['inside'], RM['fixed']))
print('在迴圈外挑波長，被留下的樣品其實已經參與了「挑」的過程 -> 資訊洩漏，誤差被低估約 6 倍。')
print('留下 4.87 %% 時，挑到的波長 (%.1f nm) 預測差最多（誤差約 %.1f %%），主導了誠實的 RMSECV。' % (pick_in[2], e_in[2]))

fig, ax = plt.subplots(figsize=(9, 4))
xx = np.arange(4)
ax.bar(xx - 0.27, e_out, 0.27, color='lightgray', label='select outside loop (%.2f)' % RM['outside'])
ax.bar(xx, e_in, 0.27, color='red', label='select inside loop (%.2f)' % RM['inside'])
ax.bar(xx + 0.27, e_fix, 0.27, color='steelblue', label='fixed 2258 nm (%.2f)' % RM['fixed'])
ax.axhline(0, color='k'); ax.set_xticks(xx); ax.set_xticklabels(['leave out %g' % c for c in conc])
ax.set_ylabel('prediction - true (% acid)'); ax.set_title('LOOCV errors'); ax.legend(loc='upper left', fontsize=8)
plt.tight_layout(); plt.show()
''')

md(r'''
**重點**：同一組資料、同一種迴歸，只因為「波長在迴圈外挑」，RMSECV 就從 2.6 % 變成 0.43 %。
固定 2258 nm 的 2.18 % 與第 1 章相同；它與誠實版（2.58 %）同量級，才是這組資料真正的預測力。

### 動手改改看（第 4 節）
1. 為什麼「誠實」版的每折會挑到不同的波長？這說明挑出來的波長有多穩定？
2. 固定 2258 nm 的 RMSECV（2.18 %）與誠實的 2.58 % 接近，這說明了什麼？
3. 若你在迴圈外挑波長，報告 RMSECV = 0.43 %，審稿人該不該相信？
''')

md(r'''
## 5. 學生練習：自己寫誠實 LOOCV、改變挑選池

**練習 A**：把挑選池從「可用波段」改成**全部 142 個波長**（含飽和），重做「迴圈內挑波長」的 LOOCV，看 RMSECV 變好還是變差。
下面已替你包成函式 `loocv_inside(ok)`，其中 `ok` 就是挑選池（布林向量）；請你自己改池、試試看，並讀懂函式每一行。

**選做**：把 4 × 142 的光譜換成 **142 個互相獨立的隨機 N(0,1) 欄位**（「142 個毫無關係的隨機波長」），重複 1000 次，
看「挑最高 R²」的中位數。
''')

code(r'''
def loocv_inside(ok):
    """誠實 LOOCV：每折只用 3 個訓練點，在挑選池 ok 裡重新挑最佳波長。回傳 RMSECV。"""
    e = []
    for k in range(4):
        tr = [i for i in range(4) if i != k]
        jr = best_idx(r2_all(tr, conc[tr]), ok)
        e.append(inv_pred(tr, k, jr) - conc[k])
    return rmse(e)

RM_usable = loocv_inside(usable)
RM_all = loocv_inside(np.ones(n_wl, bool))
print('迴圈內挑波長：只用可用波段 RMSECV = %.3f %%；用全部 142 個 RMSECV = %.3f %%' % (RM_usable, RM_all))

# ---- 選做：142 個毫無關係的隨機波長 ----
rng = np.random.default_rng(SEED)
def max_r2_random(n_cols=n_wl):
    Z = rng.standard_normal((4, n_cols)); Zc = Z - Z.mean(axis=0); yc = conc - conc.mean()
    return float(np.max((yc @ Zc) ** 2 / ((yc @ yc) * (Zc ** 2).sum(axis=0))))   # 每欄對 conc 的 R2，取最大

max_r2_rand = np.array([max_r2_random() for _ in range(N_RAND)])
med_rand = float(np.median(max_r2_rand))
print('隨機資料「最大 R2」的中位數 = %.3f（亂數不同，數字會有些微差異）' % med_rand)
print('這些欄位與濃度完全無關，卻挑得到 R2 約 0.99：n = 4 時，「最大 R2」只反映挑了多少個候選。')
''')

md(r'''
### 動手改改看（第 5 節）
1. 把 142 改成 10 個隨機欄位，「最大 R²」的中位數變多少？
2. 把樣品數 4 改成 12（自己造 12 個濃度），隨機資料的最大 R² 還有 0.99 嗎？
3. 寫一個函式：給定挑選池與樣品數，回傳誠實 LOOCV 的 RMSECV，並比較不同池的結果。
''')

md(r'''
## 6. 自我驗證

與答案鍵逐項比對（R² 容許 ±0.001，RMSECV 與誤差容許 ±0.002，計數須完全相同；選做的隨機實驗只檢查大致落在 0.99 附近）。
若改過上面的參數，部分項目 FAIL 是正常的，改回預設、整本重跑即可。
''')

code(r'''
checks = []
def chk(name, got, key, tol):
    checks.append((abs(got - key) <= tol, name, float(got), key))

chk('usable wavelengths (A < 1.5)', usable.sum(), 88, 0)
chk('best of all 142: wavelength (nm)', wl[j_all], 2309.7, 0.05)
chk('best of all 142: R2', R[j_all], 0.9985, 0.001)
chk('best of all 142: spectra with A >= 1.5', (X[:, j_all] >= A_MAX).sum(), 2, 0)
chk('best of all 142: A of 10% bottle', X[3, j_all], 1.575, 0.001)
chk('wavelengths with R2 > 0.95', n_over95, 26, 0)
chk('best of usable: wavelength (nm)', wl[j_use], 1862.4, 0.05)
chk('best of usable: R2', R[j_use], 0.9957, 0.001)
chk('fixed 2258: wavelength (nm)', wl[i_fixed], 2258.1, 0.05)
chk('fixed 2258: R2', R[i_fixed], 0.9342, 0.001)
chk('shuffle: number of orderings', len(shuf), 24, 0)
chk('shuffle: wrong orderings, best R2 > 0.9', n_wrong_best, 18, 0)
chk('shuffle: wrong orderings, fixed R2 > 0.9', n_wrong_fixed, 1, 0)
chk('RMSECV, selection outside loop (%)', RM['outside'], 0.432, 0.002)
chk('RMSECV, selection inside loop (%)', RM['inside'], 2.579, 0.002)
chk('RMSECV, fixed 2258 nm (%)', RM['fixed'], 2.182, 0.002)
chk('inside-loop picks match (of 4 folds)', (np.abs(pick_in - np.array([1853.9, 1871.0, 1796.1, 2233.1])) < 0.05).sum(), 4, 0)
chk('outside errors: max abs diff to key', np.abs(e_out - np.array([-0.704, 0.409, -0.098, -0.271])).max(), 0, 0.002)
chk('inside errors: max abs diff to key', np.abs(e_in - np.array([-1.147, 0.694, -4.169, 2.725])).max(), 0, 0.002)
chk('fixed errors: max abs diff to key', np.abs(e_fix - np.array([1.935, -0.685, -1.199, 3.66])).max(), 0, 0.002)
chk('optional: random median max R2 (approx 0.99)', med_rand, 0.99, 0.01)

print('%-6s | %-48s | %12s | %12s' % ('result', 'item', 'yours', 'key'))
print('-' * 88)
for ok, name, got, key in checks:
    print('%-6s | %-48s | %12.5g | %12.5g' % ('PASS' if ok else 'FAIL', name, got, key))
print('-' * 88)
print('TOTAL: %d / %d PASS' % (sum(c[0] for c in checks), len(checks)))
''')

md(r'''
### 總結
1. n = 4 時，從 142 個波長挑 R² 最高的，總能挑到 R² > 0.99；連 142 個**隨機**波長也辦得到。
2. 全部波長一起挑，會挑到**飽和**波段（2309.7 nm，A > 1.5）；要先排除不可信的波段。
3. 洗牌檢驗：濃度亂排時，「挑最佳波長」仍有 18/23 種得到 R² > 0.9，固定 2258 nm 只有 1/23。
4. LOOCV 必須把「挑波長」放進迴圈：迴圈外 0.43 %，迴圈內 2.58 %，固定 2258 nm 2.18 %。
5. 結論：樣品數很少時，**用化學知識事先決定波長**，比用資料挑選更可信。
''')


def local_head():
    return [
        (MD, '> **本機 Jupyter 版**：請在 `teaching/exercises/` 目錄下開啟，資料讀自 `../data/spectra_long.csv`。\n'
             '> 需要套件：numpy、pandas、matplotlib、scipy。Colab 請改開 `ch03_wavelength_pitfall_colab.ipynb`。'),
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
             '> 2. 選單「**檔案 → 上傳筆記本**」，選擇 `ch03_wavelength_pitfall_colab.ipynb`。\n'
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
    build('ch03_wavelength_pitfall.ipynb', local_head())
    build('ch03_wavelength_pitfall_colab.ipynb', colab_head())
