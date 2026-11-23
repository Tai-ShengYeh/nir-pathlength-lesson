# -*- coding: utf-8 -*-
# 第 8 章：由「同一份」cell 定義生成本機版與 Colab 版筆記本。用法：python build_notebooks_ch08.py
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
# 第 8 章 參考方法的誤差：滴定 vs Atago，程式練習

**目標**：前面幾章都把「參考值」當成真值。這一章反過來問：**參考值本身有多準？**
算出四組參考值、滴定的誤差預算、參考值彼此差多少、NIR 能不能當旁證、NaOH 力價的比例誤差、
滴定量怎麼設計、各種誤差大小比一比，最後產生一份**滴定紀錄表**（`ch08_titration_sheet.csv`）。
所有計算都是確定性的，**沒有任何亂數**，每個人都應該算出同一個答案。

| 節 | 內容 |
|---|---|
| 0 | 載入資料與二階導數（第 4、7 節用到光譜） |
| 1 | 四組參考值：標示、Atago 11/02、Atago 11/03、滴定 |
| 2 | 誤差預算 |
| 3 | 參考值彼此差多少 |
| 4 | NIR 旁證很弱 |
| 5 | 力價的比例誤差 |
| 6 | 滴定量的設計 |
| 7 | 誤差比一比 |
| 8 | 產生滴定紀錄表 |

只用 numpy、pandas、matplotlib、scipy。圖上文字一律用英文（Colab 沒有中文字型）。
每節結尾有「動手改改看」，並附該節的**自我驗證 PASS/FAIL**（相對誤差 1e-3；答案為 0 時絕對誤差 1e-3）。

**事先固定的設定**：滴定資料（另一組學生，11/03）：取 10 mL 醋稀釋到 100 mL，再取 25 mL，以 NaOH 滴定；
NaOH 由稱重配製（0.4315 g ÷ 40 ÷ 0.1 L → 力價 f = 1.07875，未標定）。醋酸當量 0.006 g/mL（0.1 N × 60.05 g/mol ÷ 1000），結果單位是 w/v %。
器材允差採 A 級玻璃器皿的典型值（示意），標準不確定度 = 允差 ÷ √3。
NIR：1 mm、重複 1 的 4 瓶，二階導數 Savitzky–Golay 9 點、二次多項式，取 2258 nm，乘 10³。
''')

md(r'''
## 0. 載入資料與二階導數
''')

code(r'''
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from scipy.signal import savgol_filter
from scipy import stats

np.set_printoptions(precision=5, linewidth=130)
plt.rcParams['figure.dpi'] = 90

# ---- 事先固定的設定 ----
SG_WIN = 9
WL0 = 2258

# ---- 自我驗證的小工具：每節結尾呼叫 show_checks() ----
checks = []
_shown = 0
def chk(name, got, key):
    """容許誤差：相對 1e-3；答案是 0 時用絕對 1e-3（答案只有四位有效數字）。"""
    tol = 1e-3 if key == 0 else 1e-3 * abs(key)
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
g = (df[(df['replicate'] == 1) & np.isclose(df['path_mm'], 1.0)]
     .pivot_table(index='acid_pct', columns='wavelength_nm', values='absorbance').sort_index())
wl = g.columns.values; X1 = g.values                                               # 4 瓶 x 波長（Atago 11/02：0, 3.32, 4.87, 7.45 %）
D1 = savgol_filter(X1, SG_WIN, 2, deriv=2, axis=1)
j = int(np.argmin(np.abs(wl - WL0)))
sig = D1[:, j] * 1e3                                                              # 2258 nm 訊號 y（x 10^-3）
print('1 mm、重複 1：%d 瓶 x %d 個波長；2258 nm 取第 %d 點（%.1f nm）' % (X1.shape[0], X1.shape[1], j, wl[j]))
print('訊號 y（x1e-3）：', np.round(sig, 4))

def fit_line(xs, ys):
    n = len(xs); b, a = np.polyfit(xs, ys, 1)
    s = float(np.sqrt(((ys - a - b * xs) ** 2).sum() / (n - 2)))
    return dict(a=float(a), b=float(b), s=s, s_conc=s / abs(b), sxx=float(((xs - xs.mean()) ** 2).sum()))
''')

md(r'''
### 動手改改看（第 0 節）
1. 這一章的第 1–3、5、6 節完全不需要光譜，只用實驗常數。為什麼參考方法的誤差可以不靠 NIR 就估出來？
2. 為什麼第 4 節只用 3 個點、而且自由度只有 1？
''')

md(r'''
## 1. 四組參考值

- **Atago** 讀值來自 ATAGO PAL-Easy ACID：這是以**導電度**為基礎的酸度計，**不是折射計**；結果是 w/v %。
- **滴定**：酸度 (w/v %) = V × f × 0.006 × (100/25) × (100/10)。
  V = NaOH 用量 (mL)；f = NaOH 力價（以稱重配製、未標定）；0.006 g/mL = 每 mL 0.1 N NaOH 對應的醋酸克數；
  100/25 = 100 mL 瓶裡只取 25 mL 滴定；100/10 = 10 mL 醋稀釋到 100 mL。
- 講義用 0.006 g/mL；精確值 0.1 N × 60.05 ÷ 1000 = 0.006005，差 0.08 %，可忽略。
''')

code(r'''
V = np.array([8.8, 17.0, 44.0])                           # NaOH 用量（mL）
f = 0.4315 / 40 / 0.1 * 10                                # 1.07875
K = f * 0.006 * (100 / 25) * (100 / 10)                   # 每 mL NaOH 對應的酸度 %
titr = V * K
label = np.array([3.0, 5.0, 10.0]); at1 = np.array([3.32, 4.87, 7.45]); at2 = np.array([3.01, 4.12, 6.86])
ratio = titr / at2
daydiff = (at1 - at2) / at1 * 100
print('f = 0.4315 / 40 / 0.1 * 10 = %.5f；K = f x 0.006 x 4 x 10 = %.5f %%/mL' % (f, K))
print('%-6s %7s %10s %10s %9s %9s %11s %12s' % ('bottle', 'label', 'Atago11/02', 'Atago11/03', 'V (mL)', 'titration', 'titr/Atago03', 'day diff (%)'))
for i in range(3):
    print('%-6d %7.2f %10.2f %10.2f %9.1f %9.3f %11.4f %12.3f' % (i + 1, label[i], at1[i], at2[i], V[i], titr[i], ratio[i], daydiff[i]))
eq_exact = 0.1 * 60.05 / 1000
print('講義用 0.006 g/mL；精確值 0.1 N x 60.05 / 1000 = %.6f g/mL，差 %.2f %%（可忽略）。' % (eq_exact, (eq_exact - 0.006) / eq_exact * 100))
print('同一批醋，隔一天用同一支 Atago 量，讀值就差 8-15 %；而滴定和標示、Atago 三者又彼此對不上。第 3 節會量化。')

chk('f (NaOH factor, unstandardized)', f, 1.079)
chk('K (% per mL NaOH)', K, 0.2589)
for i in range(3): chk('titration acidity, bottle %d (V=%.1f mL)' % (i + 1, V[i]), titr[i], [2.278, 4.401, 11.39][i])
for i in range(3): chk('titration / Atago 11/03, bottle %d' % (i + 1), ratio[i], [0.7569, 1.068, 1.661][i])
for i in range(3): chk('Atago day difference (%%), bottle %d' % (i + 1), daydiff[i], [9.337, 15.4, 7.919][i])
chk('0.1 N x 60.05 / 1000 (g/mL)', eq_exact, 0.006005)
chk('aliquot factor 100/25', 100 / 25, 4)
chk('dilution factor 100/10', 100 / 10, 10)
show_checks('第 1 節')
''')

md(r'''
### 動手改改看（第 1 節）
1. 為什麼 3 個瓶子的 titr/Atago 比值是 0.76、1.07、1.66，而不是同一個常數？（如果只是力價不準，比值應該怎樣？）
2. Atago 兩天的差 8–15 %，是儀器不穩，還是樣品變了？要怎麼設計實驗分辨？
3. 把 0.006 改成 0.006005，所有滴定酸度改變多少 %？
''')

md(r'''
## 2. 誤差預算

單次滴定的相對標準不確定度。器材允差採 **A 級玻璃器皿的典型值（示意，請以你手上器材標示為準）**；標準不確定度 = 允差 ÷ √3（矩形分布）。
滴定體積不確定度 u_V = √(2 × 0.02² + u_滴定管² + u_終點²)：兩次讀數各 ±0.02 mL、滴定管允差、終點一滴。

**這是「下限」，沒有算進力價。** NaOH 由稱重配製、沒有標定；NaOH 會吸濕、吸 CO₂，真正的力價很可能比 1.07875 **低 2–5 %**，
算出的酸度就**同比例偏高**。這是「比例偏差」（所有瓶子同一比例），不是隨機誤差。以 KHP 標定、三重複後，力價的相對不確定度約 0.2 %（典型值，示意）。
''')

code(r'''
u_read = 0.02                                             # 每次讀滴定管 +/- 0.02 mL
tol = dict(buret=0.05, endpoint=0.05, pip10=0.02, flask100=0.08, pip25=0.03)
u = {k: v / np.sqrt(3) for k, v in tol.items()}
uV = float(np.sqrt(2 * u_read ** 2 + u['buret'] ** 2 + u['endpoint'] ** 2))
print('滴定體積不確定度 uV = sqrt(2 x %.2f^2 + u_buret^2 + u_endpoint^2) = %.5f mL（兩次讀數各 %.2f mL；滴定管允差 %.2f、終點一滴 %.2f mL）'
      % (u_read, uV, u_read, tol['buret'], tol['endpoint']))
rel_titre = uV / V * 100
rel_pip10 = np.repeat(u['pip10'] / 10 * 100, 3); rel_flask = np.repeat(u['flask100'] / 100 * 100, 3); rel_pip25 = np.repeat(u['pip25'] / 25 * 100, 3)
total_rel = np.sqrt(rel_titre ** 2 + rel_pip10 ** 2 + rel_flask ** 2 + rel_pip25 ** 2)
total_abs = titr * total_rel / 100
print('%-8s %8s %8s %9s %9s %9s %10s %10s' % ('V (mL)', 'titre %', 'pip10 %', 'flask100 %', 'pip25 %', 'total %', 'acidity %', 'total abs'))
for i in range(3):
    print('%-8.1f %8.4f %8.4f %9.4f %9.4f %9.4f %10.3f %10.5f' % (V[i], rel_titre[i], rel_pip10[i], rel_flask[i], rel_pip25[i], total_rel[i], titr[i], total_abs[i]))
print('結論：器材造成的不確定度只有 0.2-0.6 %（絕對值約 0.013-0.021 % 酸度），體積小時由滴定管讀數主導，體積大時由 10 mL 吸量管主導。')
print('這是『下限』：它沒有算進力價。NaOH 由稱重配製、沒有標定；NaOH 會吸濕、吸 CO2，真正的力價很可能比 1.07875 低 2-5 %。')
print('這是另一種誤差：它是『比例偏差』（所有瓶子同一比例），不是隨機誤差；算出的酸度會同比例偏高。以 KHP 標定、三重複後，力價的相對不確定度約 0.2 %（典型值，示意）。')

fig, ax = plt.subplots(figsize=(8.6, 4.8))
rows_ = [('titre', rel_titre, 'steelblue'), ('pip10', rel_pip10, 'silver'), ('flask100', rel_flask, 'grey'), ('pip25', rel_pip25, 'dimgrey'), ('total', total_rel, 'firebrick')]
w = 0.16
for k, (nm, vals, col) in enumerate(rows_):
    ax.bar(np.arange(3) + (k - 2) * w, vals, w, label=nm, color=col)
ax.set_xticks(range(3)); ax.set_xticklabels(['V = %.1f mL' % v for v in V]); ax.set_ylim(0, 0.7)
ax.set_ylabel('relative standard uncertainty (%)'); ax.set_title('Titration error budget (Class A tolerances, illustrative)'); ax.legend(fontsize=8)
plt.tight_layout(); plt.show()

chk('uV: titre volume uncertainty (mL)', uV, 0.04967)
for i in range(3): chk('rel. u titre (%%), V=%.1f' % V[i], rel_titre[i], [0.5644, 0.2922, 0.1129][i])
for i in range(3): chk('rel. u 10 mL pipette (%%), V=%.1f' % V[i], rel_pip10[i], [0.1155, 0.1155, 0.1155][i])
for i in range(3): chk('rel. u 100 mL flask (%%), V=%.1f' % V[i], rel_flask[i], [0.04619, 0.04619, 0.04619][i])
for i in range(3): chk('rel. u 25 mL pipette (%%), V=%.1f' % V[i], rel_pip25[i], [0.06928, 0.06928, 0.06928][i])
for i in range(3): chk('total relative uncertainty (%%), V=%.1f' % V[i], total_rel[i], [0.5821, 0.325, 0.1817][i])
for i in range(3): chk('total absolute uncertainty (%% acid), V=%.1f' % V[i], total_abs[i], [0.01326, 0.0143, 0.0207][i])
show_checks('第 2 節')
''')

md(r'''
### 動手改改看（第 2 節）
1. 哪一項在 V = 8.8 mL 時最大？在 V = 44 mL 時呢？為什麼？
2. 如果真正的 NaOH 力價比 1.07875 低 3 %，3 個瓶子算出的酸度各偏多少（% 酸度）？這比器材誤差大幾倍？
3. 把 0.006 g/mL 的 0.08 % 差異也當成誤差，總不確定度會明顯變大嗎？
''')

md(r'''
## 3. 參考值彼此差多少

做法：如果某個參考值是真的，滴定時應該用掉 `need = 參考值 / K` mL 的 NaOH；實際用量是 V。
`gap = need − V`；`gap/uV` 是差了幾個讀數不確定度。
''')

code(r'''
refs = {'label': label, 'atago_1103': at2, 'atago_1102': at1}
need = {k: r / K for k, r in refs.items()}
gap = {k: need[k] - V for k in refs}
gap_u = {k: gap[k] / uV for k in refs}
for nm in refs:
    print('[%s] need (mL)：%s' % (nm, ' '.join('%7.2f' % x for x in need[nm])))
    print('%s gap (mL)：%s；gap/uV：%s' % (' ' * (len(nm) + 2), ' '.join('%7.3f' % x for x in gap[nm]), ' '.join('%8.2f' % x for x in gap_u[nm])))
all_gu = np.abs(np.concatenate([gap_u[k] for k in refs]))
print('9 個差距的大小介於 %.0f 到 %.0f 倍 uV。' % (all_gu.min(), all_gu.max()))
print('結論：差距是讀數不確定度的 22-352 倍，不可能是讀滴定管的誤差。')
print('第 1 瓶（標示 3 %）的滴定值（2.28 %）低於標示、Atago 兩天的讀值，三個都比它高：它是最可疑的一瓶，應該重測；但這不等於『證明滴定錯了』。')

fig, ax = plt.subplots(figsize=(8.6, 4.8))
w = 0.26
for k, (nm, col) in enumerate([('label', 'grey'), ('atago_1103', 'steelblue'), ('atago_1102', 'darkorange')]):
    ax.bar(np.arange(3) + (k - 1) * w, gap_u[nm], w, label={'label': 'label', 'atago_1103': 'Atago 11/03', 'atago_1102': 'Atago 11/02'}[nm], color=col)
ax.axhline(0, color='k', lw=0.8); ax.set_xticks(range(3)); ax.set_xticklabels(['bottle %d (V = %.1f mL)' % (i + 1, V[i]) for i in range(3)])
ax.set_ylabel('gap / uV  (titre needed - titre observed)'); ax.set_title('Gap of each reference to the titration, in reading uncertainties', fontsize=10); ax.legend(fontsize=8)
plt.tight_layout(); plt.show()

key_need = {'label': [11.59, 19.31, 38.62], 'atago_1103': [11.63, 15.91, 26.5], 'atago_1102': [12.82, 18.81, 28.78]}
key_gap = {'label': [2.787, 2.312, -5.375], 'atago_1103': [2.826, -1.087, -17.5], 'atago_1102': [4.023, 1.81, -15.22]}
key_gu = {'label': [56.13, 46.56, -108.2], 'atago_1103': [56.9, -21.88, -352.4], 'atago_1102': [81.01, 36.45, -306.5]}
for nm in refs:
    for i in range(3): chk('%s: need (mL), bottle %d' % (nm, i + 1), need[nm][i], key_need[nm][i])
    for i in range(3): chk('%s: gap (mL), bottle %d' % (nm, i + 1), gap[nm][i], key_gap[nm][i])
    for i in range(3): chk('%s: gap / uV, bottle %d' % (nm, i + 1), gap_u[nm][i], key_gu[nm][i])
chk('smallest |gap| / uV', all_gu.min(), 21.88)
chk('largest |gap| / uV', all_gu.max(), 352.4)
chk('bottle 1: references above the titration (of 3)', sum(int(gap[k][0] > 0) for k in refs), 3)
show_checks('第 3 節')
''')

md(r'''
### 動手改改看（第 3 節）
1. 如果只看 gap_in_uV，哪一個瓶子、哪一個參考的差最大？這能說明是哪一種方法錯嗎？
2. 第 3 瓶（V = 44 mL）三個 gap 都是負的。這個瓶子的滴定值（11.39 %）比標示 10 % 還高，可能的原因有哪些？
3. 為什麼說『最可疑』不等於『證明錯了』？還缺什麼證據？
''')

md(r'''
## 4. NIR 旁證很弱

1 mm、2258 nm 二階導數 ×10³（重複 1）。用**水＋第 2、3 瓶**（以滴定值當 x）建 3 點直線，預測第 1 瓶，並附 Fieller 95 % 區間（殘差自由度只有 1）。
''')

code(r'''
ref4 = np.concatenate([[0.0], titr])                      # 水 + 3 瓶的滴定值
keep = [0, 2, 3]                                          # 水、第 2 瓶、第 3 瓶；留下第 1 瓶當未知
xs, ys = ref4[keep], sig[keep]; n3 = 3
f3 = fit_line(xs, ys); b3, a3, s3 = f3['b'], f3['a'], f3['s']
sxx3 = float(((xs - xs.mean()) ** 2).sum()); tq = float(stats.t.ppf(.975, n3 - 2)); kk = tq * tq * s3 * s3; y0 = sig[1]
g3 = float(kk / (b3 * b3 * sxx3))
A = b3 * b3 - kk / sxx3; B = -2 * b3 * (y0 - a3) + 2 * kk * xs.mean() / sxx3
C = (y0 - a3) ** 2 - kk * (1 + 1 / n3) - kk * xs.mean() ** 2 / sxx3
pred3 = float((y0 - a3) / b3)
lo3 = hi3 = float('nan')
if A > 0:
    r = np.sqrt(B * B - 4 * A * C); lo3, hi3 = sorted([float((-B - r) / (2 * A)), float((-B + r) / (2 * A))])
print('3 點：水 + 第 2、3 瓶（用滴定值當 x）；殘差自由度 = 1，t(0.975, 1) = %.3f' % tq)
print('預測第 1 瓶：%.3f %%（滴定值是 %.3f %%）；g = %.4f（<1，Fieller 區間有界）' % (pred3, titr[0], g3))
print('Fieller 95 %% 區間：%.2f 到 %.2f %%' % (lo3, hi3))
print('結論：區間從 -14 到 13 %，什麼都證明不了；它同時包含『滴定對』與『標示對』。')
print('而且這是循環論證：用 NIR 來評斷參考方法，可是 NIR 的校正又是靠這些參考方法建立的。NIR 不能當裁判。')

chk('NIR prediction of bottle 1 (%)', pred3, 3.091)
chk('g statistic', g3, 0.6414)
chk('Fieller lower limit (%)', lo3, -14.25)
chk('Fieller upper limit (%)', hi3, 12.66)
show_checks('第 4 節')
''')

md(r'''
### 動手改改看（第 4 節）
1. 預測值 3.09 % 剛好接近標示 3 %、遠離滴定 2.28 %。這算不算『NIR 支持標示』？看區間再回答。
2. 自由度 1 的 t(0.975) 是多少？如果有 8 個校正點（自由度 6），區間會縮小到多少量級？
3. 為什麼用 NIR 評斷參考方法是循環論證？
''')

md(r'''
## 5. 力價的比例誤差

把參考值全部乘 (1 + e)，e = −3 %：R² 不變，斜率與預測同比例偏。解法是用 KHP（鄰苯二甲酸氫鉀，M = 204.22 g/mol）標定：
f = m / (204.22 × 0.1 × V/1000)。
''')

code(r'''
e = -0.03
r2f = lambda x, y: float(np.corrcoef(x, y)[0, 1] ** 2)
fa0 = fit_line(ref4, sig); fa1 = fit_line(ref4 * (1 + e), sig)
y0f = 4.20
r2_0 = r2f(ref4, sig); r2_1 = r2f(ref4 * (1 + e), sig)
pred0 = (y0f - fa0['a']) / fa0['b']; pred1 = (y0f - fa1['a']) / fa1['b']
khp_g = 0.1 * 0.020 * 204.22                              # 20 mL 的 0.1 N NaOH 需要的 KHP 克數
print('e = %.2f：R^2 %.4f -> %.4f（完全不變）；斜率 %.4f -> %.4f；訊號 y0 = %.2f 的預測 %.3f -> %.3f %%' % (e, r2_0, r2_1, fa0['b'], fa1['b'], y0f, pred0, pred1))
print('結論：比例誤差不會讓 R^2 變差，所以光看 R^2 完全看不出力價錯了；它只會讓斜率與預測同比例偏（預測 = 原來 x (1 + e)）。')
print('解法：用 KHP（鄰苯二甲酸氫鉀，M = 204.22 g/mol，一級標準）標定 NaOH，f = m / (204.22 x 0.1 x V/1000)。')
print('標定 20 mL 的 0.1 N NaOH 需要 KHP：0.1 x 0.020 x 204.22 = %.4f g。' % khp_g)

chk('R^2 (original)', r2_0, 0.9898)
chk('R^2 (references x 0.97)', r2_1, 0.9898)
chk('slope (original)', fa0['b'], -0.6415)
chk('slope (references x 0.97)', fa1['b'], -0.6614)
chk('prediction at y0 = 4.20 (original, %)', pred0, 6.131)
chk('prediction at y0 = 4.20 (x 0.97, %)', pred1, 5.947)
chk('prediction ratio scaled / original', pred1 / pred0, 0.97)
chk('KHP mass for 20 mL of 0.1 N NaOH (g)', khp_g, 0.4084)
show_checks('第 5 節')
''')

md(r'''
### 動手改改看（第 5 節）
1. 為什麼 R² 不變？（提示：R² 只看點排成直線的程度，與 x 軸整體伸縮無關）
2. 如果真正的力價低 5 %，預測會偏多少 %？
3. 若只有 R² 能看，要怎麼發現力價誤差？（提示：需要獨立的標準，例如 KHP）
''')

md(r'''
## 6. 滴定量的設計

D mL 醋稀釋到 100 mL，取 25 mL 滴定。滴定體積 V = 酸度 ÷ (f × 0.006 × 4 × (100/D))。
目標 15–40 mL（50 mL 滴定管）：太小，讀數誤差佔比大；超過 50 mL 要重新裝滿滴定管。**D 越大，滴定量越大。**
''')

code(r'''
concs = [3, 5, 7, 11.39]; dil = [5, 10, 20]
titre_fn = lambda c, D: c / (f * 0.006 * (100 / 25) * (100 / D))
Vm = np.array([[titre_fn(c, D) for D in dil] for c in concs])
Rm = uV / Vm * 100
inrange = (Vm >= 15) & (Vm <= 40)
over50 = Vm > 50
print('%-10s %s' % ('acidity %', '      '.join('D = %2d mL' % d for d in dil)))
for i, c in enumerate(concs):
    print('%-10.2f' % c, end='')
    for jj in range(len(dil)):
        print(' %7.2f mL %s' % (Vm[i, jj], '(>50)' if over50[i, jj] else ('(OK) ' if inrange[i, jj] else '     ')), end='')
    print()
print('15-40 mL 的格子有 %d 個；超過 50 mL 的有 %d 個。' % (inrange.sum(), over50.sum()))
print('相對讀數不確定度 uV / V（%）：體積越大越小，D = 10 mL、V 約 20 mL 時約 0.25 %。')
print('D 越大，滴定量越大。D = 10 mL 對 5-7 % 的醋很合適（19-27 mL）；3 % 的醋 D = 10 mL 只有 11.6 mL（太小），改用 D = 20 mL；11.39 % 的醋 D = 10 mL 要 44 mL，改用 D = 5 mL 才有 22 mL。')

cls = np.where(Vm > 50, 4, np.where(Vm > 40, 3, np.where(Vm >= 15, 2, 1)))
fig, ax = plt.subplots(figsize=(8, 5))
ax.imshow(cls, cmap=ListedColormap(['gainsboro', 'mediumseagreen', 'khaki', 'firebrick']), vmin=0.5, vmax=4.5, aspect='auto')
for i in range(len(concs)):
    for jj in range(len(dil)):
        ax.text(jj, i, '%.1f' % Vm[i, jj], ha='center', va='center')
ax.set_xticks(range(3)); ax.set_xticklabels(dil); ax.set_yticks(range(4)); ax.set_yticklabels(concs)
ax.set_xlabel('vinegar volume D (mL) diluted to 100 mL'); ax.set_ylabel('acidity (%)')
ax.set_title('Titre volume V (mL): green = 15-40 mL target; red = over 50 mL')
plt.tight_layout(); plt.show()

key_V = [[5.794, 11.59, 23.17], [9.656, 19.31, 38.62], [13.52, 27.04, 54.07], [22.0, 43.99, 87.99]]
key_R = [[0.8572, 0.4286, 0.2143], [0.5143, 0.2572, 0.1286], [0.3674, 0.1837, 0.09185], [0.2258, 0.1129, 0.05645]]
for i, c in enumerate(concs):
    for jj, D in enumerate(dil): chk('titre V (mL): %.2f %%, D=%d' % (c, D), Vm[i, jj], key_V[i][jj])
for i, c in enumerate(concs):
    for jj, D in enumerate(dil): chk('rel. reading u (%%): %.2f %%, D=%d' % (c, D), Rm[i, jj], key_R[i][jj])
chk('cells within 15-40 mL (of 12)', int(inrange.sum()), 5)
chk('cells above 50 mL (of 12)', int(over50.sum()), 2)
show_checks('第 6 節')
''')

md(r'''
### 動手改改看（第 6 節）
1. 11.39 % 的醋用 D = 20 mL 會得到幾 mL？為什麼不好？（提示：滴定管只有 50 mL）
2. 要滴 3 % 的醋，D 取多少最合適？（提示：目標 15–40 mL）
3. 滴定量 V 加倍，讀數誤差佔比減半；那為什麼不一直把 V 做得很大？
''')

md(r'''
## 7. 誤差比一比（單位都是 % 酸度）

**注意每一項的依據不同**：
- **實測**：NIR 校正殘差（4 點、自由度 2）；Atago 兩天差異的 SD（3 瓶各 2 次）；滴定與標示的差（第 1 瓶，但標示本身準不準不知道）。
- **假設**：滴定器材下限（用典型允差推估，不是量到的）。
''')

code(r'''
res = sig - (fa0['a'] + fa0['b'] * ref4)
nir_resid_titr = float(np.sqrt((res ** 2).sum() / 2) / abs(fa0['b']))     # 4 點、自由度 2
titr_floor_abs = total_abs
atago_day_sd_abs = float(np.sqrt(np.mean((at1 - at2) ** 2) / 2))
titr_3pct_gap_abs = float(3.0 - titr[0])
triplicate_factor = float(1 / np.sqrt(3))
print('NIR 校正殘差（4 點，滴定當參考，自由度 2）  ：%.4f %%   [實測]' % nir_resid_titr)
print('Atago 兩天的差異 SD（3 瓶、各 2 次）        ：%.4f %%   [實測]' % atago_day_sd_abs)
print('滴定 vs 標示：第 1 瓶（3 %% 標示）差        ：%.4f %%   [實測的差，但標示本身是否準不知道]' % titr_3pct_gap_abs)
print('滴定器材下限（3 瓶）                        ：%s %%   [假設：典型允差，不是量到的]' % ' / '.join('%.4f' % x for x in titr_floor_abs))
print('三重複可把隨機誤差縮小為 1/sqrt(3) = %.4f 倍' % triplicate_factor)
print('閱讀：器材下限（約 0.02 %）比其他幾項小 20-50 倍，但它是假設值；Atago 兩天差 0.41 % 和 NIR 殘差 0.61 % 是實測，量級一樣。')
print('所以 NIR 校正殘差裡，很大一部分可能是參考值自己的誤差，而不是光譜不好。')

ev = [titr_floor_abs[2], titr_floor_abs[1], titr_floor_abs[0], titr_3pct_gap_abs, atago_day_sd_abs, nir_resid_titr]
en = ['titration floor, V=44 mL (assumed)', 'titration floor, V=17 mL (assumed)', 'titration floor, V=8.8 mL (assumed)',
      'titration vs label, bottle 1 (measured)', 'Atago day-to-day SD (measured)', 'NIR residual, 4 pts (measured)']
ec = ['grey'] * 3 + ['steelblue'] * 3
fig, ax = plt.subplots(figsize=(9.5, 4.8))
ax.hlines(range(6), 0.005, ev, color=ec, lw=2); ax.plot(ev, range(6), 'o', ms=9, color='none')
for k in range(6):
    ax.plot(ev[k], k, 'o', ms=9, color=ec[k]); ax.text(ev[k], k + 0.2, '%.3f' % ev[k], ha='center', fontsize=8)
ax.set_xscale('log'); ax.set_xlim(0.005, 2); ax.set_yticks(range(6)); ax.set_yticklabels(en, fontsize=8)
ax.set_xlabel('error size (% acid, log scale)'); ax.set_title('Reference-method errors vs NIR residual')
plt.tight_layout(); plt.show()

chk('NIR residual, titration reference (% acid)', nir_resid_titr, 0.612)
for i in range(3): chk('titration floor, V=%.1f mL (%% acid)' % V[i], titr_floor_abs[i], [0.01326, 0.0143, 0.0207][i])
chk('Atago day-to-day SD (% acid)', atago_day_sd_abs, 0.4096)
chk('titration vs label 3 % (% acid)', titr_3pct_gap_abs, 0.7217)
chk('triplicate factor 1/sqrt(3)', triplicate_factor, 0.5774)
show_checks('第 7 節')
''')

md(r'''
### 動手改改看（第 7 節）
1. 哪些數字是『實測』，哪些是『假設』？為什麼一定要分清楚？
2. 如果滴定做三重複，器材下限變成多少？這會不會改變結論？
3. 力價的比例偏差（2–5 %）折算成酸度，對第 2 瓶（4.4 %）是多少 % 酸度？它和上面哪一項同量級？
''')

md(r'''
## 8. 產生滴定紀錄表

`ch08_titration_sheet.csv`（UTF-8 含 BOM，Excel 可直接開）共 20 列、8 欄（步驟、樣品、重複、秤重_g、起點_mL、終點_mL、用量_mL、備註）：
KHP 標定 3 次 → 空白 2 次 → 原液 3 次 → 盲樣 A、B、C、D 各 3 次。秤重、起點、終點、用量、備註 5 欄**留白**，實驗時手填。
「用量」＝終點 − 起點；樣品的酸度要先扣掉空白用量再計算。列印版流程、公式與稀釋建議見 `ch08_titration_sop.md`。
''')

code(r'''
import csv
cols = ['步驟', '樣品', '重複', '秤重_g', '起點_mL', '終點_mL', '用量_mL', '備註']
blank5 = [''] * 5
rows8 = [['KHP 標定', 'KHP', i] + blank5 for i in (1, 2, 3)]
rows8 += [['空白', '空白', i] + blank5 for i in (1, 2)]
rows8 += [['樣品滴定', '原液', i] + blank5 for i in (1, 2, 3)]
for L in 'ABCD':
    rows8 += [['樣品滴定', '盲樣' + L, i] + blank5 for i in (1, 2, 3)]
with open('ch08_titration_sheet.csv', 'w', encoding='utf-8-sig', newline='') as fh:
    w = csv.writer(fh, lineterminator='\n'); w.writerow(cols); w.writerows(rows8)
ts = pd.read_csv('ch08_titration_sheet.csv', encoding='utf-8-sig')
print('已寫出 ch08_titration_sheet.csv：%d 列，%d 欄' % (len(ts), ts.shape[1]))
print('順序：KHP 標定 3 次 -> 空白 2 次 -> 原液 3 次 -> 盲樣 A、B、C、D 各 3 次；秤重、起點、終點、用量、備註 5 欄留白，實驗時手填。')
print("『用量』= 終點 - 起點；樣品的酸度要先扣掉空白用量再計算。流程、公式與稀釋建議見 ch08_titration_sop.md。")
print(ts.iloc[[0, 3, 5, 8, 11, 14, 19], :3].to_string())

chk('sheet: number of rows (3 + 2 + 3 + 12)', len(ts), 20)
chk('sheet: number of columns', ts.shape[1], 8)
chk('sheet: KHP rows', int((ts['步驟'] == 'KHP 標定').sum()), 3)
chk('sheet: blank rows', int((ts['步驟'] == '空白').sum()), 2)
chk('sheet: stock vinegar rows', int((ts['樣品'] == '原液').sum()), 3)
chk('sheet: blind-sample rows (A-D x 3)', int(ts['樣品'].str.startswith('盲樣').sum()), 12)
chk('sheet: replicate-number mismatches', int((ts['重複'].values != np.array([1, 2, 3, 1, 2, 1, 2, 3] + [1, 2, 3] * 4)).sum()), 0)
chk('sheet: numeric cells not blank (weight/start/end/titre)', int(ts.iloc[:, 3:7].notna().sum().sum()), 0)
chk('sheet: remark cells not blank', int(ts['備註'].notna().sum()), 0)
show_checks('第 8 節')
''')

md(r'''
### 動手改改看（第 8 節）
1. 用 Excel 開 `ch08_titration_sheet.csv`，中文有亂碼嗎？（檔案含 BOM，應該不會）
2. 為什麼 KHP 要稱重 3 份、空白要做 2 次？
3. 為什麼盲樣的順序最好由別人編碼、不要依濃度排？
''')

code(r'''
print('全部自我驗證：%d / %d PASS' % (sum(c[0] for c in checks), len(checks)))
for c in checks:
    if not c[0]:
        print('FAIL:', c[1], c[2], c[3])
''')

md(r'''
### 總結
1. 四組參考值（標示、Atago 兩天、滴定）互相對不上；Atago 兩天就差 8–15 %。Atago 是電導式酸度計，不是折射計。
2. 滴定器材的不確定度只有 0.2–0.6 %（約 0.013–0.021 % 酸度）；未標定 NaOH 的力價偏差（可能偏低 2–5 %）是另一種比例偏差。
3. 參考值之間的差距是讀數不確定度的 22–352 倍，不是讀數誤差；第 1 瓶（標示 3 %）的滴定值最可疑，應重測。
4. NIR 旁證太弱（區間 −14 到 13 %），而且循環論證，不能當裁判。
5. 力價的比例誤差不影響 R²，只讓斜率與預測同比例偏；用 KHP 標定（20 mL 需約 0.41 g）。
6. 滴定量要設計到 15–40 mL；低酸度樣品加大 D（20 mL），高酸度樣品減小 D（5 mL）。
7. 比誤差時要分清楚實測與假設；NIR 殘差有一大部分可能來自參考值。
8. 已產生滴定紀錄表（20 列）；列印版流程見 `ch08_titration_sop.md`。
''')


def local_head():
    return [
        (MD, '> **本機 Jupyter 版**：請在 `teaching/exercises/` 目錄下開啟，資料讀自 `../data/spectra_long.csv`。\n'
             '> 需要套件：numpy、pandas、matplotlib、scipy。Colab 請改開 `ch08_reference_methods_colab.ipynb`。'),
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
             '> 2. 選單「**檔案 → 上傳筆記本**」，選擇 `ch08_reference_methods_colab.ipynb`。\n'
             '> 3. 依序執行每一格。\n>\n'
             '> 資料已**內嵌**在筆記本中，下一格會自動解出 `spectra_long.csv`。Colab 已內建所需套件。\n'
             '> 第 8 節產生的 `ch08_titration_sheet.csv` 會存在 Colab 的檔案面板（左側資料夾圖示），請自行下載。'),
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
    print('wrote', name, '(%d cells)' % len(cells))


if __name__ == '__main__':
    build('ch08_reference_methods.ipynb', local_head())
    build('ch08_reference_methods_colab.ipynb', colab_head())
