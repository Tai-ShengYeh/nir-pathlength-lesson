# -*- coding: utf-8 -*-
# 第 5 章：由「同一份」cell 定義生成本機版與 Colab 版筆記本。用法：python build_notebooks_ch05.py
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
# 第 5 章 預測未知樣品與不確定度：程式練習

**目標**：用單一波長 2258 nm 的二階導數訊號建立「反向校正」y = a + b·c，預測未知樣品的酸度 c₀，並且誠實地算出它的 **95 % 區間**；
看懂為什麼校正點很少（n = 3、4）時區間會大到沒有意義，甚至**無界**。

| 節 | 內容 |
|---|---|
| 0 | 載入資料、取 2258 nm 二階導數訊號、工具函式 |
| 1 | 用 Atago 4 點建立校正：斜率、殘差標準差 s、t 值、g |
| 2 | 藏起 4.87 % 那瓶，只用 3 點校正：預測值、區間無界 |
| 3 | 假想未知樣品 y₀ = 4.20：簡易區間 vs Fieller 精確區間；量 3 次 |
| 4 | 改用滴定參考值（4 點）校正 |
| 5 | 由殘差換算的偵測／定量下限（和第 1 章的 LOD 不同） |
| 6 | 設計示意：第 2 章的設計（8 濃度 × 3 次）區間會縮到多小 |
| 7 | 怎麼報告結果；學生練習 |

只用 numpy、pandas、matplotlib、scipy（`savgol_filter` 做二階導數、`stats.t` 查 t 值）。圖上文字一律用英文（Colab 沒有中文字型）。
每節結尾有「動手改改看」，並附該節的**自我驗證 PASS/FAIL**。

**方法**
- 訊號 y：1 mm、重複 1，SG 9 點二次二階導數（點距為 1）在 2258 nm 的值 × 10³。
- 反向校正：以已知濃度 c 對訊號 y 做直線 y = a + b·c；殘差標準差 s = √(Σ殘差² / (n − 2))，t = t(0.975, n − 2)，Sxx = Σ(c − c̄)²。
- 預測：c₀ = (y₀ − a) / b。
- 簡易區間（半寬）：t·s/|b| · √(1/m + 1/n + (y₀ − ȳ)² / (b²·Sxx))，m 是未知樣品重複量測的次數。
- **Fieller 精確區間**：解所有滿足 |y₀ − a − b·c| ≤ t·s·√(1/m + 1/n + (c − c̄)²/Sxx) 的 c，這是一個一元二次不等式。
- g = t²s² / (b²·Sxx)。**g ≥ 1 表示斜率在 95 % 水準下和 0 沒有顯著差異，區間無界**（任何濃度都「相容」於這個訊號）。
''')

md(r'''
## 0. 載入資料、2258 nm 二階導數訊號

- **資料**：1 mm、重複 1 的 4 個樣品（Atago 11/02 讀值 0、3.32、4.87、7.45 %）。
- **訊號**：Savitzky–Golay 9 點、二次、二階導數、點距為 1，取 2258 nm 的值再乘 10³。
- **滴定參考值**：滴定用去 NaOH 體積 V = 8.8、17、44 mL，酸度 = V × 1.07875 × 0.24（%）；純水是 0。
''')

code(r'''
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from scipy import stats

np.set_printoptions(precision=5, linewidth=130)
plt.rcParams['figure.dpi'] = 90

# ---- 事先固定的設定 ----
SG_WIN = 9            # 二階導數：Savitzky-Golay 9 點、二次多項式
WL0 = 2258            # 單一波長（nm）
Y0 = 4.20             # 假想未知樣品的訊號（x 10^-3）

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
d1 = df[(df['replicate'] == 1) & (np.isclose(df['path_mm'], 1.0))]
g1 = d1.pivot_table(index='acid_pct', columns='wavelength_nm', values='absorbance').sort_index()
wl = g1.columns.values; atago = g1.index.values                     # 4 個濃度：0, 3.32, 4.87, 7.45
D = savgol_filter(g1.values, SG_WIN, 2, deriv=2, axis=1)            # 4 x 142 二階導數
j = int(np.argmin(np.abs(wl - WL0)))                                # 最接近 2258 nm 的波長
sig = D[:, j] * 1e3                                                 # 訊號 y（x 10^-3）
titr = np.array([0.0] + [v * 1.07875 * 0.24 for v in (8.8, 17, 44)])  # 滴定參考值（%）
print('波長 = %.1f nm' % wl[j])
print('Atago 酸度（%）：', atago)
print('滴定酸度（%）  ：', np.round(titr, 3))
print('2258 nm 訊號 y ：', np.round(sig, 4))
''')

code(r'''
def calib(xs, ys):
    """直線校正 y = a + b*c。回傳字典：n、a、b、s（殘差標準差）、Sxx、c 與 y 的平均、t、斜率的 t 值 tb、g。"""
    xs = np.asarray(xs, float); ys = np.asarray(ys, float); n = len(xs)
    xbar = xs.mean(); ybar = ys.mean(); sxx = float(((xs - xbar) ** 2).sum())
    b = float(((xs - xbar) * (ys - ybar)).sum() / sxx); a = float(ybar - b * xbar)
    s = float(np.sqrt(((ys - (a + b * xs)) ** 2).sum() / (n - 2)))
    t = float(stats.t.ppf(0.975, n - 2))
    return dict(n=n, a=a, b=b, s=s, sxx=sxx, xbar=float(xbar), ybar=float(ybar), t=t,
                tb=float(abs(b) / (s / np.sqrt(sxx))), g=float(t ** 2 * s ** 2 / (b ** 2 * sxx)))

def predict(c, y0, m=1):
    """預測 c0 = (y0 - a)/b，簡易半寬，與 Fieller 精確區間（g >= 1 時 lo = hi = None，表示無界）。"""
    x0 = (y0 - c['a']) / c['b']
    simple = c['t'] * c['s'] / abs(c['b']) * np.sqrt(1 / m + 1 / c['n'] + (y0 - c['ybar']) ** 2 / (c['b'] ** 2 * c['sxx']))
    # Fieller：解 (y0 - a - b*x)^2 = t^2 s^2 (1/m + 1/n + (x - xbar)^2 / Sxx) 的二次方程 A x^2 + B x + C = 0
    k = c['t'] ** 2 * c['s'] ** 2
    A = c['b'] ** 2 - k / c['sxx']
    B = -2 * c['b'] * (y0 - c['a']) + 2 * k * c['xbar'] / c['sxx']
    C = (y0 - c['a']) ** 2 - k * (1 / m + 1 / c['n']) - k * c['xbar'] ** 2 / c['sxx']
    if A <= 0:                                   # A <= 0 等價於 g >= 1：斜率不顯著，區間無界
        lo = hi = None
    else:
        r = np.sqrt(B * B - 4 * A * C); lo, hi = sorted([(-B - r) / (2 * A), (-B + r) / (2 * A)])
    return dict(x0=float(x0), simple=float(simple), lo=None if lo is None else float(lo), hi=None if hi is None else float(hi))

def plot_calib(ax, c, xs, ys, y0=None, m=1, xlim=(-1, 9), hide=None):
    """校正直線、點、95 % 預測帶（單次量測半寬 t s sqrt(1/m + 1/n + (c-xbar)^2/Sxx)）；可加 y0 水平線。"""
    xx = np.linspace(xlim[0], xlim[1], 200); fit = c['a'] + c['b'] * xx
    hw = c['t'] * c['s'] * np.sqrt(1 / m + 1 / c['n'] + (xx - c['xbar']) ** 2 / c['sxx'])
    ax.fill_between(xx, fit - hw, fit + hw, color='steelblue', alpha=0.15, label='95%% prediction band (m=%d)' % m)
    ax.plot(xx, fit, color='steelblue', lw=1.5, label='calibration line')
    ax.plot(xs, ys, 'o', color='k', ms=7, label='calibration points')
    if hide is not None:
        ax.plot(*hide, 's', color='red', ms=9, label='hidden bottle')
    if y0 is not None:
        ax.axhline(y0, color='red', ls='--', lw=1, label='y0 = %.2f' % y0)
    ax.set_xlim(*xlim); ax.set_xlabel('reference acid (%)'); ax.set_ylabel('signal y (x 1e-3)')
    ax.legend(fontsize=8)

chk('signal y at 0 %', sig[0], 8.5087, 0.001)
chk('signal y at 3.32 %', sig[1], 6.3073, 0.001)
chk('signal y at 4.87 %', sig[2], 5.1725, 0.001)
chk('signal y at 7.45 %', sig[3], 0.9509, 0.001)
show_checks('第 0 節')
''')

md(r'''
### 動手改改看（第 0 節）
1. 訊號隨酸度**增加而變小**（斜率為負）。為什麼二階導數的吸收峰是「負的尖峰」？這會不會影響預測？
2. 把 `WL0` 改成 2250，訊號變多少？（改完記得改回 2258 再看自我驗證）
''')

md(r'''
## 1. 用 Atago 4 點建立校正

用全部 4 個樣品（0、3.32、4.87、7.45 %）對訊號做直線。n = 4，自由度 n − 2 = 2，所以 t(0.975, 2) = 4.3027 很大：
點少時，t 值本身就是不確定度的重要來源。
**斜率的 t 值** = |b| / (s/√Sxx) 檢查斜率是不是和 0 有顯著差異；**g = t²s²/(b²Sxx)** 是它的另一種寫法：g < 1 才有有限區間。
''')

code(r'''
c4 = calib(atago, sig)
print('斜率 b = %.4f（訊號/%%）；截距 a = %.4f；殘差標準差 s = %.4f' % (c4['b'], c4['a'], c4['s']))
print('t(0.975, %d) = %.4f；斜率的 t 值 = %.4f；g = %.4f' % (c4['n'] - 2, c4['t'], c4['tb'], c4['g']))
print('g < 1：斜率顯著不為 0，Fieller 區間有界。' if c4['g'] < 1 else 'g >= 1：區間無界。')

fig, ax = plt.subplots(figsize=(7, 4.4))
plot_calib(ax, c4, atago, sig, xlim=(-1, 9))
ax.set_title('Calibration with Atago reference (n = 4)')
plt.tight_layout(); plt.show()

chk('Atago 4-point: slope b', c4['b'], -0.9852, 0.002)
chk('Atago 4-point: s', c4['s'], 0.9975, 0.002)
chk('Atago 4-point: t(0.975, 2)', c4['t'], 4.3027, 0.002)
chk('Atago 4-point: slope t-statistic', c4['tb'], 5.3267, 0.002)
chk('Atago 4-point: g', c4['g'], 0.6525, 0.002)
show_checks('第 1 節')
''')

md(r'''
### 動手改改看（第 1 節）
1. 斜率的 t 值 = 5.33，聽起來「很顯著」，為什麼 g 還高達 0.65？（提示：g = t²/tb²，即 t(0.975, 2) = 4.30 與斜率 t 值的比）
2. 若有 8 個校正點、s 和 b 不變，t(0.975, 6) 是多少？g 會變多少？
''')

md(r'''
## 2. 藏起 4.87 % 那瓶：只用 3 點校正

模擬「這一瓶是未知樣品」：把 4.87 % 的樣品從校正中拿掉，只用 0、3.32、7.45 % 三點校正，再預測它。
n = 3 時自由度只剩 1，t(0.975, 1) = 12.71，區間會大得離譜。
''')

code(r'''
tr = [0, 1, 3]                                   # 校正用 0、3.32、7.45 %；藏起第 3 個（4.87 %）
c3 = calib(atago[tr], sig[tr]); p3 = predict(c3, sig[2])
print('3 點校正：b = %.4f；s = %.4f；t(0.975, 1) = %.4f；g = %.4f' % (c3['b'], c3['s'], c3['t'], c3['g']))
print('藏起來的那瓶：訊號 y = %.4f -> 預測酸度 = %.4f %%' % (sig[2], p3['x0']))
print('對照：Atago 讀值 = %.2f %%；滴定 = %.3f %%' % (atago[2], titr[2]))
if p3['lo'] is None:
    print('g = %.3f >= 1：斜率與 0 沒有顯著差異，Fieller 95 %% 區間無界（任何濃度都相容）。' % c3['g'])
else:
    print('Fieller 95 %% 區間 [%.2f, %.2f] %%' % (p3['lo'], p3['hi']))
print('點預測 3.67 % 與 Atago 4.87 %、滴定 4.40 % 差約 1 %，而且我們「沒有辦法」宣稱它有多準。')

fig, ax = plt.subplots(figsize=(7, 4.4))
plot_calib(ax, c3, atago[tr], sig[tr], xlim=(-1, 9), hide=(atago[2], sig[2]))
ax.set_title('3-point calibration, 4.87% bottle hidden')
plt.tight_layout(); plt.show()

chk('3-point: slope b', c3['b'], -1.0258, 0.002)
chk('3-point: s', c3['s'], 0.9507, 0.002)
chk('3-point: t(0.975, 1)', c3['t'], 12.7062, 0.002)
chk('3-point: g', c3['g'], 4.9775, 0.002)
chk('3-point: interval unbounded (1 = yes)', float(p3['lo'] is None), 1, 0)
chk('3-point: predicted acid of hidden bottle', p3['x0'], 3.6711, 0.002)
show_checks('第 2 節')
''')

md(r'''
### 動手改改看（第 2 節）
1. g = 4.98 遠大於 1，這在圖上看起來是什麼樣子？（預測帶與水平線 y₀ 的交會）
2. 點預測 3.67 % 與 Atago 4.87 %、滴定 4.40 % 相差 1 % 以上，卻不能說「錯了」。為什麼？（提示：區間無界）
3. 改藏第 2 個樣品（`tr = [0, 2, 3]`），預測值與 g 會怎麼變？
''')

md(r'''
## 3. 假想未知樣品 y₀ = 4.20（4 點 Atago 校正）

假設量到一個未知樣品，訊號 y₀ = 4.20（× 10⁻³）。比較：
- **簡易區間**：c₀ ± t·s/|b| · √(1/m + 1/n + (y₀ − ȳ)²/(b²Sxx))（把「y 的不確定度」直接除以斜率）；
- **Fieller 精確區間**：解二次不等式，g 較大時**不對稱**，也能偵測無界。

再看「同一個未知樣品量 3 次（m = 3，取平均）」區間縮多少。
''')

code(r'''
p4 = predict(c4, Y0); p4m3 = predict(c4, Y0, 3)
print('預測 x0 = %.4f %%' % p4['x0'])
print('簡易區間：%.4f ± %.4f -> [%.4f, %.4f]' % (p4['x0'], p4['simple'], p4['x0'] - p4['simple'], p4['x0'] + p4['simple']))
print('Fieller 區間 (m = 1)：[%.4f, %.4f]' % (p4['lo'], p4['hi']))
print('Fieller 區間 (m = 3)：[%.4f, %.4f]' % (p4m3['lo'], p4m3['hi']))
print('Fieller 區間不對稱：x0 到下限 %.2f，到上限 %.2f。' % (p4['x0'] - p4['lo'], p4['hi'] - p4['x0']))
print('下限是負的酸度（−1.68 %）沒有物理意義：區間太寬，代表資訊不足，不是「可能是負的」。')

fig, ax = plt.subplots(figsize=(8, 4.6))
plot_calib(ax, c4, atago, sig, y0=Y0, xlim=(-3, 17))
ax.axvspan(p4['lo'], p4['hi'], ymin=0, ymax=0.04, color='red', alpha=0.6, label='Fieller 95% interval (m=1)')
ax.plot([p4['x0']], [Y0], 'r*', ms=12)
ax.legend(fontsize=8, loc='upper right')
ax.set_title('Prediction of y0 = 4.20 with Fieller interval (n = 4)')
plt.tight_layout(); plt.show()

chk('y0 = 4.20: predicted x0', p4['x0'], 4.9604, 0.002)
chk('y0 = 4.20: simple half-width', p4['simple'], 4.9442, 0.002)
chk('y0 = 4.20: Fieller lower', p4['lo'], -1.6831, 0.002)
chk('y0 = 4.20: Fieller upper', p4['hi'], 15.548, 0.002)
chk('y0 = 4.20, m = 3: Fieller lower', p4m3['lo'], 0.7828, 0.002)
chk('y0 = 4.20, m = 3: Fieller upper', p4m3['hi'], 13.0822, 0.002)
show_checks('第 3 節')
''')

md(r'''
### 動手改改看（第 3 節）
1. 簡易區間 ±4.94 對稱，Fieller [−1.68, 15.55] 不對稱。哪一個比較誠實？g = 0.65 愈接近 1 差異愈大，為什麼？
2. 量 3 次（m = 3）只把下限從 −1.68 拉到 0.78，不是縮成 1/√3。為什麼？（提示：1/m 只是根號裡三項之一，另兩項來自校正本身）
3. 把 `Y0` 改成 8.0（接近 0 % 的訊號）或 1.0（超出校正範圍），區間怎麼變？
''')

md(r'''
## 4. 改用滴定參考值校正（4 點）

Atago 讀值來自 PAL-Easy ACID（以導電度換算酸度，不是折射計），本身有讀值誤差；滴定是較準的參考方法。把校正的橫軸從 Atago 換成滴定酸度（0、2.278、4.401、11.392 %），看 s、g 與預測的變化。
''')

code(r'''
c4t = calib(titr, sig); p4t = predict(c4t, Y0)
print('滴定 4 點：b = %.4f；s = %.4f；g = %.4f' % (c4t['b'], c4t['s'], c4t['g']))
print('y0 = 4.20：x0 = %.4f %%；Fieller 區間 [%.4f, %.4f]；簡易 ±%.4f' % (p4t['x0'], p4t['lo'], p4t['hi'], p4t['simple']))
print('殘差 s 從 0.9975 降到 0.3927，g 從 0.65 降到 0.095，區間縮窄到 [3.16, 9.45]（寬度約 6.3 %；Atago 是 17.2 %）。')
print('代價：滴定的最高點 11.39 % 拉長了範圍（Sxx 變大），也代表參考值的準確度直接決定預測的上限。')

fig, ax = plt.subplots(figsize=(8, 4.6))
plot_calib(ax, c4t, titr, sig, y0=Y0, xlim=(-1, 13))
ax.axvspan(p4t['lo'], p4t['hi'], ymin=0, ymax=0.04, color='red', alpha=0.6, label='Fieller 95% interval (m=1)')
ax.plot([p4t['x0']], [Y0], 'r*', ms=12)
ax.legend(fontsize=8, loc='upper right')
ax.set_xlabel('titration acid (%)'); ax.set_title('Calibration against titration reference (n = 4)')
plt.tight_layout(); plt.show()

chk('titration 4-point: slope b', c4t['b'], -0.6415, 0.002)
chk('titration 4-point: s', c4t['s'], 0.3927, 0.002)
chk('titration 4-point: g', c4t['g'], 0.0954, 0.002)
chk('titration, y0 = 4.20: predicted x0', p4t['x0'], 6.1309, 0.002)
chk('titration, y0 = 4.20: Fieller lower', p4t['lo'], 3.1568, 0.002)
chk('titration, y0 = 4.20: Fieller upper', p4t['hi'], 9.4453, 0.002)
show_checks('第 4 節')
''')

md(r'''
### 動手改改看（第 4 節）
1. 同一個訊號 y₀ = 4.20，用 Atago 校正預測 4.96 %，用滴定校正預測 6.13 %。這兩個數字差在哪裡？報告時必須說明用的是哪一種參考值。
2. s 變小的原因可能是什麼？（提示：Atago 的讀值解析度與滴定的不確定度）
''')

md(r'''
## 5. 由殘差換算的偵測／定量下限

常見的經驗式：偵測下限 LOD = 3 s/|b|，定量下限 LOQ = 10 s/|b|，s 是校正的殘差標準差，b 是斜率。
這裡的 s 是**殘差**，包含了裝樣（fill）的差異與參考值本身的誤差，所以和第 1 章用「儀器雜訊」算出的 LOD **不同**，通常大很多。
''')

code(r'''
lod_a, loq_a = 3 * c4['s'] / abs(c4['b']), 10 * c4['s'] / abs(c4['b'])
lod_t, loq_t = 3 * c4t['s'] / abs(c4t['b']), 10 * c4t['s'] / abs(c4t['b'])
print('Atago 參考：LOD = %.3f %%；LOQ = %.3f %%' % (lod_a, loq_a))
print('滴定參考  ：LOD = %.3f %%；LOQ = %.3f %%' % (lod_t, loq_t))
print('第 1 章的 LOD 只含儀器雜訊（重複掃描同一裝樣）；這裡的 s 還含裝樣、參考值與校正的誤差，才是「實際使用時」能看到的下限。')

chk('Atago residual LOD (3s/|b|)', lod_a, 3.038, 0.005)
chk('Atago residual LOQ (10s/|b|)', loq_a, 10.125, 0.005)
chk('titration residual LOD (3s/|b|)', lod_t, 1.836, 0.005)
chk('titration residual LOQ (10s/|b|)', loq_t, 6.120, 0.005)
show_checks('第 5 節')
''')

md(r'''
### 動手改改看（第 5 節）
1. Atago 的 LOQ = 10.1 % 比校正範圍的上限 7.45 % 還大。這代表什麼？（我們的 4 個校正點，沒有一個真的「定量」得夠準）
2. 這個 LOD 與第 1 章的儀器雜訊 LOD 差在哪些誤差來源？
''')

md(r'''
## 6. 設計示意：用第 2 章的設計，區間會縮到多小？

第 2 章的設計：0–7 % 共 8 個濃度（0, 1, …, 7 %），每個濃度裝樣 3 次（n = 24）。下面**假設 s、b 不變**（只換 n 與 Sxx），預測同一個 y₀ = 4.20。
這只是示意（假設較樂觀）：真實資料的 s、b 會因參考值、裝樣而不同，但它顯示了「多放幾個校正點、拉開範圍」的威力。
''')

code(r'''
def design(c, m):
    xs = np.repeat(np.arange(8.0), 3); n = len(xs); sxx = float(((xs - xs.mean()) ** 2).sum()); t = float(stats.t.ppf(0.975, n - 2))
    cc = dict(c, n=n, sxx=sxx, xbar=float(xs.mean()), ybar=float(c['a'] + c['b'] * xs.mean()), t=t,
              g=float(t ** 2 * c['s'] ** 2 / (c['b'] ** 2 * sxx)))
    return dict(g=cc['g'], **predict(cc, Y0, m))

d_a1 = design(c4, 1); d_a3 = design(c4, 3); d_t3 = design(c4t, 3)
print('Atago，設計 8 濃度 x 3 次，m = 1：g = %.4f；x0 = %.3f；區間 [%.4f, %.4f]' % (d_a1['g'], d_a1['x0'], d_a1['lo'], d_a1['hi']))
print('Atago，設計 8 濃度 x 3 次，m = 3：區間 [%.4f, %.4f]' % (d_a3['lo'], d_a3['hi']))
print('滴定，設計 8 濃度 x 3 次，m = 3：g = %.4f；x0 = %.3f；區間 [%.4f, %.4f]' % (d_t3['g'], d_t3['x0'], d_t3['lo'], d_t3['hi']))

rows = [('Atago n=4, m=1', p4), ('Atago n=4, m=3', p4m3), ('Atago design, m=1', d_a1), ('Atago design, m=3', d_a3),
        ('titration n=4, m=1', p4t), ('titration design, m=3', d_t3)]
fig, ax = plt.subplots(figsize=(9, 4.2))
for i, (lab, p) in enumerate(rows):
    ax.plot([p['lo'], p['hi']], [i, i], '-', lw=5, color='steelblue' if 'Atago' in lab else 'darkorange', solid_capstyle='butt')
    ax.plot([p['x0']], [i], 'k|', ms=14)
ax.set_yticks(range(len(rows))); ax.set_yticklabels([r[0] for r in rows]); ax.invert_yaxis()
ax.set_xlabel('predicted acid (%) with 95% Fieller interval'); ax.set_title('Interval width for y0 = 4.20')
plt.tight_layout(); plt.show()

chk('design (24 points), Atago: g', d_a1['g'], 0.035, 0.002)
chk('design, Atago m = 1: Fieller lower', d_a1['lo'], 2.8134, 0.002)
chk('design, Atago m = 1: Fieller upper', d_a1['hi'], 7.2134, 0.002)
chk('design, Atago m = 3: Fieller lower', d_a3['lo'], 3.6741, 0.002)
chk('design, Atago m = 3: Fieller upper', d_a3['hi'], 6.3526, 0.002)
chk('design, titration: g', d_t3['g'], 0.0128, 0.002)
chk('design, titration m = 3: Fieller lower', d_t3['lo'], 5.3266, 0.002)
chk('design, titration m = 3: Fieller upper', d_t3['hi'], 7.0033, 0.002)
show_checks('第 6 節')
''')

md(r'''
### 動手改改看（第 6 節）
1. 區間從 [−1.68, 15.55]（n = 4）縮到 [2.81, 7.21]（24 點）。縮窄的最主要原因是 t 值、n，還是 Sxx？（分別改 `design` 裡的 `t`、`n`、`sxx` 試試）
2. 設計若只有 8 個濃度、每個 1 次（n = 8），區間是多少？
3. 為什麼說這個示意「偏樂觀」？（提示：s 與 b 是從 4 點估出來的；真實的 24 點會發現新的誤差來源）
''')

md(r'''
## 7. 怎麼報告結果

一個預測值如果沒有**區間**、**n**、**參考方法**、**波長**，就不算完整的報告。下面把 Atago 4 點的預測整理成一行。
''')

code(r'''
def fmt(v):
    return ('%.2f' % v).replace('-', '−')

def report(x0, p, n, ref, wl_nm=2258):
    """格式化一行報告；區間無界時直接寫無界。"""
    if p['lo'] is None:
        return '酸度 %s %%（95 %% 區間無界，n = %d，%s 參考值，單一波長 %d nm）' % (fmt(x0), n, ref, wl_nm)
    return '酸度 %s %%（95 %% 區間 %s 至 %s %%，n = %d，%s 參考值，單一波長 %d nm）' % (fmt(x0), fmt(p['lo']), fmt(p['hi']), n, ref, wl_nm)

print(report(p4['x0'], p4, c4['n'], 'Atago'))
print(report(p3['x0'], p3, c3['n'], 'Atago'))

# 學生練習：把上面的報告改成「滴定參考值、4 點」的版本，並再寫一行「設計 8 濃度 x 3 次、量 3 次」的版本。
# print(report(p4t['x0'], p4t, c4t['n'], '滴定'))
''')

md(r'''
### 動手改改看（第 7 節）
1. 取消註解最後一行，改成滴定參考值的報告；比較 Atago 與滴定兩行，哪個比較能讓讀者信任？
2. 加上第 6 節的設計示意，寫出一行「設計 8 濃度 × 3 次、量 3 次」的報告（注意：這是示意，不是實測）。
3. 一份報告若只寫「酸度 4.96 %」，讀者無法判斷什麼？
''')

code(r'''
print('全部自我驗證：%d / %d PASS' % (sum(c[0] for c in checks), len(checks)))
for c in checks:
    if not c[0]:
        print('FAIL:', c[1], c[2], c[3])
''')

md(r'''
### 總結
1. 反向校正 y = a + b·c，預測 c₀ = (y₀ − a)/b；4 點 Atago：b = −0.985、s = 1.00、g = 0.65。
2. 藏起一瓶、只用 3 點：g = 4.98 ≥ 1，區間無界；點預測 3.67 % 不能宣稱準或不準。
3. y₀ = 4.20：簡易 ±4.94 對稱；Fieller [−1.68, 15.55] 不對稱且誠實；量 3 次只縮到 [0.78, 13.08]。
4. 滴定參考值（s = 0.39）區間縮到 [3.16, 9.45]；參考值的品質決定預測的品質。
5. 殘差換算的 LOD / LOQ（Atago 3.04 / 10.1 %；滴定 1.84 / 6.12 %）與第 1 章的儀器雜訊 LOD 不同。
6. 區間縮窄靠「更多、範圍更廣的校正點」（8 濃度 × 3 次：[2.81, 7.21]），而不是只靠重複量測。
''')


def local_head():
    return [
        (MD, '> **本機 Jupyter 版**：請在 `teaching/exercises/` 目錄下開啟，資料讀自 `../data/spectra_long.csv`。\n'
             '> 需要套件：numpy、pandas、matplotlib、scipy。Colab 請改開 `ch05_prediction_uncertainty_colab.ipynb`。'),
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
             '> 2. 選單「**檔案 → 上傳筆記本**」，選擇 `ch05_prediction_uncertainty_colab.ipynb`。\n'
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
    build('ch05_prediction_uncertainty.ipynb', local_head())
    build('ch05_prediction_uncertainty_colab.ipynb', colab_head())
