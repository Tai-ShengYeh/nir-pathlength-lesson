# -*- coding: utf-8 -*-
# 第 9 章：由「同一份」cell 定義生成本機版與 Colab 版筆記本。用法：python build_notebooks_ch09.py
import base64, gzip, os
import nbformat as nbf

HERE = os.path.dirname(os.path.abspath(__file__))
DATA_CSV = os.path.join(HERE, '..', 'data', 'spectra_long.csv')
TEMPLATE_MD = os.path.join(HERE, 'ch09_validation_report_template.md')
MD, CODE = 'md', 'code'
COMMON = []


def md(s):
    COMMON.append((MD, s.strip('\n')))


def code(s):
    COMMON.append((CODE, s.strip('\n')))


md(r'''
# 第 9 章 方法確效報告：現有證據夠不夠？程式練習

**目標**：把前面八章的結果收在一起，回答：「這個 1 mm、2258 nm 的二階導數 NIR 方法，現在能不能說『已確效』？」
依序檢查用途與目標、選擇性與範圍、線性、準確度、精密度、LOD/LOQ、預測區間（現在 vs 下一次實驗 vs 目標），
最後產生一份**確效總表**（`ch09_validation_scorecard.csv`，後兩欄留給學生做完下一次實驗後填）。
配套的報告範本是 `ch09_validation_report_template.md`。所有計算都是確定性的，**沒有任何亂數**。

| 節 | 內容 |
|---|---|
| 0 | 載入資料、二階導數與共用函式 |
| 1 | 用途與目標（±0.5 % 是示意） |
| 2 | 選擇性與範圍：兩個醋酸帶、四種參考值的 R² |
| 3 | 線性：4 點、沒有重複 |
| 4 | 準確度：留一法 RMSECV |
| 5 | 精密度：量過什麼、沒量什麼 |
| 6 | LOD / LOQ：殘差法 vs 儀器雜訊法 |
| 7 | 預測區間：現在 vs 下一次實驗 vs 目標 |
| 8 | 確效總表與報告範本 |

只用 numpy、pandas、matplotlib、scipy。圖上文字一律用英文（Colab 沒有中文字型）。
每節結尾有「動手改改看」，並附該節的**自我驗證 PASS/FAIL**（相對誤差 1e-3；答案為 0 時絕對誤差 1e-3）。

**數字的來源（和第 9 章講義一致）**：
- **重算**（從 `spectra_long.csv`）：二階導數 R²、單波長留一法 RMSECV、1 mm 2258 nm 最大吸光度、殘差法與儀器雜訊法的 LOD/LOQ、純誤差與缺適性的自由度。
- **取自前面各章**（在程式裡寫死並註明出處）：PLS RMSECV（第 4 章）、Atago 兩天差異 SD（第 8 章）、現在的預測區間（第 5 章）、下一次實驗的區間半寬（第 7 章）、掃描 RMS 雜訊（第 6 章）。
''')

md(r'''
## 0. 載入資料、二階導數與共用函式
''')

code(r'''
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter

np.set_printoptions(precision=5, linewidth=130)
plt.rcParams['figure.dpi'] = 90

# ---- 事先固定的設定 ----
SG_WIN = 9
WL0, WL1 = 2258, 1670
TARGET_HW = 0.5          # 示意的用途目標：預測區間半寬 +/- 0.5 %（不是法規或客戶要求）

# ---- 取自前面各章的數字（寫死；括號內是出處，講義 CH9 用的是同一批）----
PLS_RMSECV = 0.317       # 第 4 章：1 mm 二階導數、PLS 1 成分、4 瓶留一法 RMSECV
ATAGO_DAY_SD = 0.4096    # 第 8 章：Atago 同 3 瓶隔天重測的差異 SD（% 酸度）
NOW_LO, NOW_HI, NOW_X0 = 3.157, 9.445, 6.131     # 第 5 章：4 瓶滴定值校正、未知樣品 y0 = 4.20 的 95 % Fieller 區間（填 1 次）
PLAN_S = np.array([0.1, 0.2, 0.4, 0.6, 1.0])                         # 第 7 章：預試驗可能量到的填樣標準差（% 酸度）
PLAN_HW_M3 = np.array([0.1274, 0.2548, 0.5107, 0.7688, 1.296])       # 第 7 章：8 濃度 x 3 次、未知樣品填 3 次、4 % 的區間半寬
PLAN_S_REF = np.array([0.0, 0.2, 0.4])                               # 第 7 章：再加上參考值（標準液）誤差，填樣 SD 固定 0.4 %
PLAN_HW_LV = np.array([0.6033, 0.6307, 0.7086])                      # 第 7 章：對應的區間半寬
SCAN_RMS_AU = 0.001236   # 第 6 章：同一瓶重複掃描的 RMS 雜訊（吸光度）
SIGMA_CH01 = 2.7e-4      # 第 1 章：平坦區 1550-1650 nm 二階導數的雜訊 sigma（1 mm）

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
j2258 = int(np.argmin(np.abs(wl - WL0))); j1670 = int(np.argmin(np.abs(wl - WL1)))
sig2258 = D1[:, j2258] * 1e3                                                       # 2258 nm 訊號 y（x 10^-3）
sig1670 = D1[:, j1670] * 1e3
print('1 mm、重複 1：%d 瓶 x %d 個波長；2258 nm 取 %.1f nm、1670 nm 取 %.1f nm' % (X1.shape[0], X1.shape[1], wl[j2258], wl[j1670]))
print('2258 nm 訊號 y（x1e-3）：', np.round(sig2258, 4))
print('1670 nm 訊號 y（x1e-3）：', np.round(sig1670, 4))

# 四組參考值（含水 = 0）
K = 0.4315 / 40 / 0.1 * 10 * 0.006 * (100 / 25) * (100 / 10)      # 滴定：每 mL NaOH 對應的酸度 %（第 8 章，f = 1.07875）
refs = {'atago': np.array([0, 3.32, 4.87, 7.45]),                  # ATAGO PAL-Easy ACID 181（導電度式酸度計）11/02
        'titr': np.array([0, 8.8, 17, 44]) * K,                    # 滴定（NaOH 稱重配製、未標定、只做一次）
        'label': np.array([0, 3, 5, 10.0]),                        # 瓶身標示
        'pal30s': np.array([0, 2.8, 4.8, 10.4])}                   # PAL-30S 折射計，11/17 量，直接顯示酸度 %
ref_names = {'atago': 'Atago', 'titr': 'titration', 'label': 'label', 'pal30s': 'PAL-30S'}

def fit_line(xs, ys):
    n = len(xs); b, a = np.polyfit(xs, ys, 1)
    s = float(np.sqrt(((ys - a - b * xs) ** 2).sum() / (n - 2)))
    return dict(a=float(a), b=float(b), s=s, s_conc=s / abs(b), sxx=float(((xs - xs.mean()) ** 2).sum()))

def r2f(x, y):
    return float(np.corrcoef(x, y)[0, 1] ** 2)

def loo_rmse(x, y):
    """反向校正的留一法：每次留一瓶，用其餘 3 瓶擬合 訊號 = a + b x，再把被留下瓶的訊號換算回濃度，與它的參考值比。"""
    e = []
    for i in range(len(x)):
        m = np.arange(len(x)) != i
        ft = fit_line(x[m], y[m]); e.append((y[i] - ft['a']) / ft['b'] - x[i])
    return float(np.sqrt(np.mean(np.square(e))))

print('滴定係數 K = %.5f %%/mL；滴定酸度（%%）：%s' % (K, ' / '.join('%.3f' % v for v in refs['titr'])))
''')

md(r'''
### 動手改改看（第 0 節）
1. 為什麼要把 1670 nm 也算進來？（提示：第 2 章用兩個醋酸帶互相對照）
2. 把 `SG_WIN` 改成 11，R² 會變嗎？
''')

md(r'''
## 1. 用途與目標

- **用途（這門課的假設）**：用 1 mm 光徑、2258 nm 二階導數，由 NIR 光譜估計醋的酸度（w/v %）。
- **目標（示意）**：預測區間半寬 ±0.5 % 酸度。這個 0.5 % 是為了讓後面能對照而設的**示意目標**，不是法規、不是客戶要求。
  真正的確效要先由用途決定可接受的誤差，再回頭看證據夠不夠；目標改了，結論也會改。
- 確效 = 事先寫好目標與判斷準則，再用實驗證據逐項檢查。現在只有 4 瓶（水 + 3 瓶醋）、每瓶 1 次，所以多半只能說「初步」。
''')

code(r'''
print('用途（這門課的假設）：用 1 mm 光徑、2258 nm 二階導數，由 NIR 光譜估計醋的酸度（w/v %）。')
print('目標（示意）：預測區間半寬 +/- %.1f %% 酸度。這個 %.1f %% 是為了讓後面能對照而設的『示意目標』，不是法規、不是客戶要求；' % (TARGET_HW, TARGET_HW))
print('  真正的確效要先由用途決定可接受的誤差，再回頭看證據夠不夠。目標改了，結論也會改。')
print('確效 = 事先寫好目標與判斷準則，再用實驗證據逐項檢查。現在只有 4 瓶（水 + 3 瓶醋）、每瓶 1 次，所以多半只能說『初步』。')
print('各項要檢查：選擇性、範圍、線性、準確度、精密度、LOD/LOQ、穩健性、系統適用性（第 8 節的總表）。')
chk('number of calibration samples (water + 3 bottles)', X1.shape[0], 4)
show_checks('第 1 節')
''')

md(r'''
### 動手改改看（第 1 節）
1. 如果用途是「判斷醋是否在 4–6 %」，目標半寬 0.5 % 合理嗎？如果用途是「區分 3 % 與 10 %」呢？
2. 為什麼目標要在看數據**之前**訂？
''')

md(r'''
## 2. 選擇性與範圍

- **Atago** = ATAGO PAL-Easy ACID 181：**導電度式**酸度計，不是折射計。**PAL-30S** = **折射計**（11/17 量），直接顯示酸度 %。
- 兩個醋酸帶一致是**一致性**，不是**選擇性**。
- 標示與 PAL-30S 的 R² 都比滴定高，所以 n = 4 的 R² **不能用來選哪個參考方法比較對**。
''')

code(r'''
print('Atago = ATAGO PAL-Easy ACID 181：導電度式酸度計，不是折射計。PAL-30S = 折射計（11/17 量），直接顯示酸度 %。')
sigs = {'2258': sig2258, '1670': sig1670}
r2_tab = {w: {k: r2f(x, y) for k, x in refs.items()} for w, y in sigs.items()}
print('%-10s %10s %10s' % ('reference', '2258 nm', '1670 nm'))
for k in refs: print('%-10s %10.4f %10.4f' % (ref_names[k], r2_tab['2258'][k], r2_tab['1670'][k]))
print('\n閱讀：')
print('  1. 兩個醋酸帶（2258、1670 nm）對同一批參考值的 R^2 都很高，這只是『一致性』（兩個譜帶都跟濃度一起變），')
print('     不是『選擇性』。選擇性是指加入糖、其他酸、色素等干擾物後，訊號仍只反映醋酸；這件事沒有測過。')
print('  2. 標示與 PAL-30S 的 R^2（0.996、0.998）都比滴定（0.990）高；Atago（0.934）最低。')
print('     只有 4 點（含水）的 R^2 無法用來選『哪個參考方法比較對』：R^2 高只代表 4 點排得像直線，')
print('     不代表參考值準（第 8 章：比例誤差不改變 R^2）。同一組光譜，換參考值 R^2 就從 0.93 變到 0.998，')
print('     差異來自參考值本身的誤差與非線性，不是 NIR 變好或變壞。')
Amax = float(X1[:, j2258].max()); A_LIMIT = 1.5
print('\n範圍：1 mm、2258 nm 的最大原始吸光度 %.3f（< %.1f，沒有飽和）；已量到的範圍 Atago 0-%.2f %%、滴定 0-%.2f %%。'
      % (Amax, A_LIMIT, refs['atago'].max(), refs['titr'].max()))
print('  滴定的 11.39 % 只靠一次有爭議的滴定（第 8 章），不宜當作已證明的範圍上限；下一次實驗計畫 0-7 %、8 個等間隔濃度。')

fig, ax = plt.subplots(figsize=(9.5, 5.2))
xk = np.arange(4); wbar = 0.38
for i, (w, c) in enumerate((('2258', 'steelblue'), ('1670', 'silver'))):
    v = [r2_tab[w][k] for k in refs]
    ax.bar(xk + (i - 0.5) * wbar, v, wbar, color=c, label='%s nm' % w)
    for xx, vv in zip(xk + (i - 0.5) * wbar, v): ax.text(xx, vv + 0.0008, '%.3f' % vv, ha='center', fontsize=8)
ax.set_ylim(0.9, 1.012); ax.set_xticks(xk); ax.set_xticklabels([ref_names[k] for k in refs])
ax.set_ylabel('R^2 (2nd derivative vs reference)'); ax.legend(loc='upper left')
ax.set_title('Same spectra, different reference: R^2 with n = 4 cannot pick a reference', fontsize=10)
plt.tight_layout(); plt.show()

key_r2 = {'2258': {'atago': 0.9342, 'titr': 0.9898, 'label': 0.9955, 'pal30s': 0.9984},
          '1670': {'atago': 0.9448, 'titr': 0.9866, 'label': 0.9983, 'pal30s': 0.999}}
for w in ('2258', '1670'):
    for k in refs: chk('R^2 %s nm vs %s' % (w, ref_names[k]), r2_tab[w][k], key_r2[w][k])
chk('label and PAL-30S R^2 both above titration (2258 nm)', float(r2_tab['2258']['label'] > r2_tab['2258']['titr'] and r2_tab['2258']['pal30s'] > r2_tab['2258']['titr']), 1)
chk('max raw absorbance, 1 mm at 2258 nm', Amax, 1.441)
chk('calibration range, titration upper end (%)', round(float(refs['titr'].max()), 2), 11.39)
show_checks('第 2 節')
''')

md(r'''
### 動手改改看（第 2 節）
1. 「兩個譜帶都跟濃度一起變」為什麼不能證明選擇性？要加什麼樣品才能檢驗？
2. 標示的 R² 比滴定高，能不能推論「標示比滴定準」？
3. 11.39 % 這個點拿掉後，滴定的範圍變成多少？
''')

md(r'''
## 3. 線性

4 點（含水）、每濃度 1 次，所以**純誤差自由度 = 0**，無法做缺適性檢定。
同一濃度的重複如果只是「同一瓶換個填樣」，會共用配製誤差而低估變異；現在 4 瓶各只配一次，R² 同時包含配製、填樣、參考值與光譜誤差，無法拆開。
''')

code(r'''
fa = fit_line(refs['atago'], sig2258); ft = fit_line(refs['titr'], sig2258)
n_obs = 4; n_lvl = 4
pure_df_now = n_obs - n_lvl; resid_df_now = n_obs - 2
print('現在：%d 瓶、%d 個濃度、每濃度 1 次。殘差自由度 = %d；純誤差自由度 = n - 濃度數 = %d；缺適性自由度 = %d。'
      % (n_obs, n_lvl, resid_df_now, pure_df_now, resid_df_now - pure_df_now))
print('R^2：Atago %.4f、滴定 %.4f（2258 nm）。沒有重複，就沒有『純誤差』可以拿來檢定直線是否夠用（缺適性檢定）。' % (r2f(refs['atago'], sig2258), r2f(refs['titr'], sig2258)))
n_plan_lvl, n_plan_rep = 8, 3; n_plan = n_plan_lvl * n_plan_rep
pure_df_plan = n_plan - n_plan_lvl; lof_df_plan = n_plan_lvl - 2
print('下一次（第 7 章）：%d 濃度 x %d 次 = %d 觀測；純誤差自由度 %d；缺適性自由度 %d；才有條件檢定直線與曲線。' % (n_plan_lvl, n_plan_rep, n_plan, pure_df_plan, lof_df_plan))
print('提醒：同一濃度的重複，如果是『同一次配製、同一瓶、換個填樣』，它們共用配製誤差；')
print('      這樣得到的純誤差只含填樣與儀器誤差，會低估真正的樣品間變異。重複要盡量包含獨立配製。')
print('      同樣地，現在 4 瓶各只配一次，R^2 同時包含配製、填樣、參考值與光譜誤差，無法拆開。')
chk('R^2 Atago (lin)', r2f(refs['atago'], sig2258), 0.9342)
chk('R^2 titration (lin)', r2f(refs['titr'], sig2258), 0.9898)
chk('pure-error df, now', pure_df_now, 0)
chk('pure-error df, plan (8 x 3)', pure_df_plan, 16)
chk('lack-of-fit df, plan (8 levels)', lof_df_plan, 6)
show_checks('第 3 節')
''')

md(r'''
### 動手改改看（第 3 節）
1. 如果 8 個濃度只各做 1 次（8 觀測），純誤差自由度是多少？缺適性自由度呢？
2. 「同瓶重複填 3 次」與「3 瓶獨立配製」的純誤差，哪個比較大？為什麼？
''')

md(r'''
## 4. 準確度

留一法 RMSECV：4 瓶輪流留下 1 瓶，用其餘 3 瓶擬合，再把被留下瓶的訊號換算回濃度。
這**不是準確度**：沒有獨立樣品（n_test = 0），而且每次只用 3 點擬合（自由度 1）。
- PLS 的 RMSECV（0.32 %）**小於** Atago 自己兩天的差異 SD（0.41 %）：過度配適的徵兆。
- PAL-30S 的 0.26 % 只是**一致性**，不是準確度。
''')

code(r'''
print('做法：4 瓶輪流留下 1 瓶，用其餘 3 瓶擬合 訊號 = a + b x，再把被留下瓶的訊號換算回濃度；RMSECV = 這 4 個誤差的均方根。')
print('注意：這不是『準確度』：沒有獨立樣品（n_test = 0），而且每次只用 3 點擬合（自由度 1）。它是『現有資料能做到的內部檢查』。')
rm_loo = {k: loo_rmse(x, sig2258) for k, x in refs.items()}
print('%-10s %9s' % ('reference', 'RMSECV'))
for k in refs: print('%-10s %9.4f' % (ref_names[k], rm_loo[k]))
print('PLS（第 4 章，1 mm 二階導數、1 成分）RMSECV = %.3f %%' % PLS_RMSECV)
print('Atago 兩天差異 SD（第 8 章，同樣 3 瓶隔天重測）= %.4f %%' % ATAGO_DAY_SD)
print('\n閱讀：')
print('  1. PLS 的 RMSECV（%.2f %%）比 Atago 自己兩天的差異（%.2f %%）還小：NIR 預測得比參考方法自己重複的結果還一致，' % (PLS_RMSECV, ATAGO_DAY_SD))
print('     這不合理；是過度配適（用 4 個點、整段光譜擬合）的徵兆。RMSECV 不能小於參考值本身的誤差還說它『準』。')
print('  2. PAL-30S 的 0.26 % 很小，但只代表 NIR 與 PAL-30S 兩個方法『彼此一致』（consistency），不代表準確（accuracy）：')
print('     兩者有可能一起偏。準確度需要獨立的盲樣加上有標定過的滴定（第 7、8 章）。')
print('  3. 同一組光譜，RMSECV 隨參考值從 2.18（Atago）到 0.26（PAL-30S）變化 8 倍，所以 RMSECV 的大小主要反映參考值的品質。')

names = ['Atago', 'titration', 'label', 'PLS (ch04)', 'PAL-30S']
vals = [rm_loo['atago'], rm_loo['titr'], rm_loo['label'], PLS_RMSECV, rm_loo['pal30s']]
cols = ['grey', 'grey', 'grey', 'steelblue', 'grey']
fig, ax = plt.subplots(figsize=(9.5, 5))
ax.bar(range(5), vals, color=cols)
for i, v in enumerate(vals): ax.text(i, v / 2, '%.2f' % v, ha='center', color='white', fontsize=9, fontweight='bold')
ax.axhline(ATAGO_DAY_SD, ls='--', color='firebrick', lw=2)
ax.set_xticks(range(5)); ax.set_xticklabels(names); ax.set_ylim(0, 2.6); ax.set_ylabel('RMSECV (% acid)')
ax.set_title('Leave-one-out RMSECV by reference (n = 4)')
ax.legend(handles=[plt.Rectangle((0, 0), 1, 1, color='grey'), plt.Rectangle((0, 0), 1, 1, color='steelblue'),
                   plt.Line2D([0], [0], color='firebrick', ls='--', lw=2)],
          labels=['single wavelength 2258 nm', 'PLS, all wavelengths (ch04)', 'Atago day-to-day SD = %.2f' % ATAGO_DAY_SD], loc='upper right', fontsize=8)
plt.tight_layout(); plt.show()

chk('RMSECV Atago, 2258 nm (% acid)', rm_loo['atago'], 2.182)
chk('RMSECV titration, 2258 nm (% acid)', rm_loo['titr'], 1.152)
chk('RMSECV PAL-30S, 2258 nm (% acid)', rm_loo['pal30s'], 0.2559)
chk('RMSECV label, 2258 nm (% acid; ch08 nirfit)', rm_loo['label'], 0.6769)
chk('PLS RMSECV (ch04, hard-coded)', PLS_RMSECV, 0.317)
chk('Atago day-to-day SD (ch08, hard-coded)', ATAGO_DAY_SD, 0.4096)
chk('PLS RMSECV below Atago day SD (1 = yes)', float(PLS_RMSECV < ATAGO_DAY_SD), 1)
chk('independent test samples (n_test)', 0, 0)
show_checks('第 4 節')
''')

md(r'''
### 動手改改看（第 4 節）
1. 如果 Atago 兩天差異是 0.41 %，NIR 的 RMSECV 理論上最好能做到多少？（提示：參考值誤差進了校正）
2. 為什麼 PAL-30S 的 RMSECV 小，還不能說 NIR「準確」？
3. 要得到真正的準確度，需要什麼樣的樣品？
''')

md(r'''
## 5. 精密度

只量過**掃描雜訊**（0.001236 AU，第 6 章）與參考方法 Atago 的日間差異；**填樣、配製、日間、操作者、樣品槽**對 NIR 都沒量。
校正殘差 SD 混了光譜、配製、填樣與參考值誤差，不是精密度。
''')

code(r'''
resid_titr = float(np.sqrt(((sig2258 - (ft['a'] + ft['b'] * refs['titr'])) ** 2).sum() / 2) / abs(ft['b']))   # 校正殘差 SD（% 酸度，4 點、自由度 2）
print('【量過的】')
print('  - 掃描雜訊：同一瓶重複掃描的 RMS = %.6f AU（第 6 章）。這只是儀器的短期重複性。' % SCAN_RMS_AU)
print('  - 校正殘差 SD（滴定當參考、4 點、自由度 2）= %.3f %%；它混了光譜、配製、填樣與參考值誤差，不是『精密度』。' % resid_titr)
print('  - 參考方法 Atago 隔天差異 SD = %.4f %%（第 8 章）；這是參考方法的日間精密度，不是 NIR 的。' % ATAGO_DAY_SD)
print('【沒量的】')
print('  - 填樣重複性（同一瓶重新填 10 次）')
print('  - 樣品配製的重複性（獨立配 3 份）')
print('  - 日間（不同天）、不同操作者、不同樣品槽的中間精密度')
print('下一次（第 7 章預試驗）先填樣 10 次量出填樣 SD，因為它直接決定預測區間多寬（第 7 節）。')
chk('scan RMS noise (AU; ch06, hard-coded)', SCAN_RMS_AU, 0.001236)
chk('NIR calibration residual SD, titration ref (% acid)', resid_titr, 0.612)
show_checks('第 5 節')
''')

md(r'''
### 動手改改看（第 5 節）
1. 掃描雜訊 0.0012 AU 換算成酸度大概是多少 %？為什麼它不等於方法的精密度？
2. 填樣、日間、操作者三種變異，哪一種你最擔心？怎麼用最少次實驗量出來？
''')

md(r'''
## 6. LOD / LOQ

- **殘差法**：LOD = 3 s / |b|、LOQ = 10 s / |b|；s = 校正線的殘差 SD，b = 斜率。s 包含了參考值與配製的誤差。
- **儀器雜訊法**（第 1 章）：σ = 2.7×10⁻⁴（平坦區 1550–1650 nm 的二階導數雜訊，1 mm）；斜率用**原始**（沒乘 10³）的 1 mm、2258 nm 二階導數對 Atago 濃度。
''')

code(r'''
def lod_res(x, y):
    f = fit_line(x, y); return np.array([3 * f['s_conc'], 10 * f['s_conc']])
res_titr = lod_res(refs['titr'], sig2258); res_atago = lod_res(refs['atago'], sig2258)
print('殘差法：LOD = 3 s / |b|、LOQ = 10 s / |b|；s = 校正線的殘差 SD（訊號單位），b = 斜率。s 包含了參考值與配製的誤差。')
print('  滴定當參考  ：LOD %.3f %%、LOQ %.3f %%' % tuple(res_titr))
print('  Atago 當參考：LOD %.3f %%、LOQ %.3f %%（LOQ 已超過這 4 瓶的最大濃度 %.2f %%）' % (res_atago[0], res_atago[1], refs['atago'].max()))
slope_raw = fit_line(refs['atago'], D1[:, j2258])['b']      # 第 1 章用的斜率：原始（沒乘 1e3）二階導數對 Atago 濃度
lod_ch01 = 3 * SIGMA_CH01 / abs(slope_raw); loq_ch01 = 10 * SIGMA_CH01 / abs(slope_raw)
print('儀器雜訊法（第 1 章）：sigma = %.1e（平坦區 1550-1650 nm 的二階導數雜訊，1 mm）；斜率（原始、未乘 1e3）= %.7f 每 %%' % (SIGMA_CH01, slope_raw))
print('  LOD = 3 sigma / |b| = %.3f %%、LOQ = 10 sigma / |b| = %.3f %%' % (lod_ch01, loq_ch01))
print('閱讀：兩種做法差 2-4 倍。儀器雜訊法只考慮光譜雜訊，是下限；殘差法把配製、參考值的誤差都算進去，比較保守但也被 4 點、自由度 2 拖累。')
print('所以報告時要寫清楚『用哪一種方法、s 是什麼』，而且新校正做完要重算。')

lm = np.array([[res_titr[0], res_atago[0], lod_ch01], [res_titr[1], res_atago[1], loq_ch01]])
fig, ax = plt.subplots(figsize=(9.5, 5.2)); xk = np.arange(3); wbar = 0.38
for i, (lab, c) in enumerate((('LOD (3 s / |b|)', 'steelblue'), ('LOQ (10 s / |b|)', 'firebrick'))):
    ax.bar(xk + (i - 0.5) * wbar, lm[i], wbar, color=c, label=lab)
    for xx, vv in zip(xk + (i - 0.5) * wbar, lm[i]): ax.text(xx, vv + 0.15, '%.2f' % vv, ha='center', fontsize=8)
ax.set_xticks(xk); ax.set_xticklabels(['residual, titration ref', 'residual, Atago ref', 'instrument noise (ch01)'])
ax.set_ylim(0, 11.5); ax.set_ylabel('% acid'); ax.legend(loc='upper right', fontsize=8)
ax.set_title('LOD and LOQ: three ways to estimate (preliminary)')
plt.tight_layout(); plt.show()

chk('residual LOD, titration (%)', res_titr[0], 1.836)
chk('residual LOQ, titration (%)', res_titr[1], 6.12)
chk('residual LOD, Atago (%)', res_atago[0], 3.038)
chk('residual LOQ, Atago (%)', res_atago[1], 10.13)
chk('ch01 slope, 1 mm 2258 nm, raw units (AU/nm^2 per %)', slope_raw, -0.0009852)
chk('instrument-noise LOD, ch01 (%)', lod_ch01, 0.8222)
chk('instrument-noise LOQ, ch01 (%)', loq_ch01, 2.741)
show_checks('第 6 節')
''')

md(r'''
### 動手改改看（第 6 節）
1. 如果校正殘差 SD 減半，殘差法的 LOD 與 LOQ 變多少？
2. LOQ 大於已量到的最大濃度，這個 LOQ 還可信嗎？
''')

md(r'''
## 7. 預測區間：現在 vs 下一次實驗 vs 目標

現在（第 5 章）半寬約 3.14 %；下一次（第 7 章，8 濃度 × 3 次、未知樣品填 3 次）半寬取決於預試驗量到的填樣 SD；目標（示意）是 0.5 %。
''')

code(r'''
now_hw = (NOW_HI - NOW_LO) / 2
print('現在（第 5 章）：4 瓶、未知樣品填 1 次，預測 %.3f %%，95 %% 區間 %.3f 到 %.3f %%，半寬 %.3f %%。' % (NOW_X0, NOW_LO, NOW_HI, now_hw))
print('目標（示意）：半寬 %.1f %%；現在是目標的 %.1f 倍，離目標很遠。' % (TARGET_HW, now_hw / TARGET_HW))
print('下一次（第 7 章）：8 濃度 x 3 次、未知樣品填 3 次、4 %；區間半寬取決於預試驗量到的填樣 SD：')
print('  %-26s %s' % ('填樣 SD（%）', ' '.join('%6.1f' % v for v in PLAN_S)))
print('  %-26s %s' % ('區間半寬（%）', ' '.join('%6.3f' % v for v in PLAN_HW_M3)))
print('  %-26s %s' % ('達到目標 0.5 %？', ' '.join('%6s' % ('yes' if v <= TARGET_HW else 'no') for v in PLAN_HW_M3)))
print('  再加上參考值（標準液）誤差（填樣 SD 固定 0.4 %）：')
print('  %-26s %s' % ('參考值 SD（%）', ' '.join('%6.1f' % v for v in PLAN_S_REF)))
print('  %-26s %s' % ('區間半寬（%）', ' '.join('%6.3f' % v for v in PLAN_HW_LV)))
print('閱讀：')
print('  1. 8 x 3 的設計能把半寬從 %.2f %% 壓到 0.13-1.3 %%，目標 0.5 %% 只有在填樣 SD 約 0.4 %% 以下才有機會（0.4 %% 時 0.51，略超過）。' % now_hw)
print('  2. 只要把參考值誤差算進去（標準液 SD 0.2 %、0.4 %），半寬就從 0.60 變成 0.63、0.71，全都超過 0.5 %：')
print('     到那時，瓶頸是參考方法，不是光譜。')

iv = np.concatenate([[now_hw], PLAN_HW_M3, PLAN_HW_LV])
inm = ['now: 4 bottles, 1 fill (ch05)'] + ['plan 8x3, fill SD %.1f' % s for s in PLAN_S] + ['plan 8x3, fill 0.4 + ref SD %.1f' % s for s in PLAN_S_REF]
icol = ['firebrick'] + ['steelblue'] * 5 + ['darkorange'] * 3
ypos = np.arange(len(iv))[::-1]
fig, ax = plt.subplots(figsize=(10, 5.4))
ax.hlines(ypos, 0.08, iv, color=icol, lw=2)
for y_, v_, c_ in zip(ypos, iv, icol):
    ax.plot(v_, y_, 'o', ms=9, color=c_); ax.text(v_ * 1.07, y_, '%.2f' % v_, ha='left', va='center', fontsize=8)
ax.axvline(TARGET_HW, ls='--', lw=2, color='darkgreen'); ax.text(TARGET_HW * 1.05, len(iv) - 0.6, 'illustrative target 0.5', color='darkgreen', fontsize=8)
ax.set_xscale('log'); ax.set_xlim(0.08, 4.5); ax.set_ylim(-0.6, len(iv) - 0.2)
ax.set_yticks(ypos); ax.set_yticklabels(inm, fontsize=8)
ax.set_xlabel('95 % prediction interval half-width (% acid, log scale)'); ax.set_title('Prediction interval: now vs next experiment vs target')
plt.tight_layout(); plt.show()

chk('now: interval lower limit (%; ch05, hard-coded)', NOW_LO, 3.157)
chk('now: interval upper limit (%; ch05, hard-coded)', NOW_HI, 9.445)
chk('now: prediction x0 (%; ch05, hard-coded)', NOW_X0, 6.131)
chk('now: half-width (%)', now_hw, 3.144)
for i, s in enumerate(PLAN_S): chk('plan half-width, fill SD %.1f (%%)' % s, PLAN_HW_M3[i], [0.1274, 0.2548, 0.5107, 0.7688, 1.296][i])
for i, s in enumerate(PLAN_S_REF): chk('plan half-width, ref SD %.1f + fill 0.4 (%%)' % s, PLAN_HW_LV[i], [0.6033, 0.6307, 0.7086][i])
chk('plan designs within the 0.5 % target (of 8)', int((np.concatenate([PLAN_HW_M3, PLAN_HW_LV]) <= TARGET_HW).sum()), 2)
chk('now half-width / target', now_hw / TARGET_HW, 6.288)
show_checks('第 7 節')
''')

md(r'''
### 動手改改看（第 7 節）
1. 如果把目標放寬到 ±1 %，下一次的 8 種情況有幾個達標？
2. 為什麼說「加入參考值誤差後，瓶頸是參考方法」？要改善要做什麼？（提示：第 8 章）
3. 現在的區間是 3.16–9.45 %，覆蓋了 3 % 到 10 % 的大部分標示值；這樣的方法有用途嗎？
''')

md(r'''
## 8. 確效總表與報告範本

`ch09_validation_scorecard.csv`（UTF-8 含 BOM，Excel 可直接開）共 8 列、6 欄：項目、狀態、現在的證據、怎麼補、下一次實驗結果、結論。
狀態用文字：**沒有** / **部分或初步** / **有證據**。**下一次實驗結果**與**結論**兩欄留白，做完第 7 章的實驗後由學生填寫。
報告範本是 `ch09_validation_report_template.md`（一個驗證項目一個小節，「現有證據」已預填）。

**能說**：1 mm、2258 nm 的二階導數與滴定值（n = 4）呈線性關係；這是**可行性研究**的初步證據。
**不能說**：「方法已確效」；「準確度是 x %」；「可取代滴定」。
''')

code(r'''
status_txt = ['沒有', '部分或初步', '有證據']                    # 0 / 1 / 2
card = [
    ('選擇性', 0, '兩個醋酸帶一致只是一致性；沒有測過糖等干擾物', '市售盲樣、加入干擾物（第 7 章）'),
    ('線性', 1, '4 點、沒有重複，無法做缺適性檢定', '8 濃度 × 3 次（第 7 章）'),
    ('範圍', 1, '2258 nm 在 1 mm 未飽和；11.4 % 只靠一次有爭議的滴定', '0–7 %、8 個等間隔濃度'),
    ('準確度', 0, '沒有獨立樣品；參考值彼此對不上', '盲樣＋標定過的滴定（第 7、8 章）'),
    ('精密度', 0, '只有掃描雜訊；填樣、日間誤差沒量', '預試驗填樣 10 次（第 7 章）'),
    ('LOD／LOQ', 1, '殘差法 1.8 / 6.1 %、儀器雜訊法 0.8 %（初步）', '新校正重算'),
    ('穩健性', 0, '溫度、操作者、樣品槽都沒測', '刻意改變條件的小實驗'),
    ('系統適用性', 0, '水帶規則在 1 mm 只能當參考（水帶飽和）', '每次實驗量共同對照（第 7 章）'),
]
cols9 = ['項目', '狀態', '現在的證據', '怎麼補', '下一次實驗結果', '結論']
lines = [','.join(cols9)]
for item, c, now, how in card:
    fields = [item, status_txt[c], now, how, '', '']
    assert not any(ch in fl for fl in fields for ch in ',"\n')           # 欄位內不含逗號、引號、換行，不需要加引號（三個語言輸出才會逐位元一致）
    lines.append(','.join(fields))
with open('ch09_validation_scorecard.csv', 'w', encoding='utf-8-sig', newline='') as fh:
    fh.write('\n'.join(lines) + '\n')
sc = pd.read_csv('ch09_validation_scorecard.csv', encoding='utf-8-sig', dtype=str, keep_default_na=False)
print('已寫出 ch09_validation_scorecard.csv：%d 列，%d 欄' % sc.shape)
print('狀態：沒有 %d 項、部分或初步 %d 項、有證據 %d 項（8 項裡沒有任何一項是『有證據』）。'
      % ((sc['狀態'] == '沒有').sum(), (sc['狀態'] == '部分或初步').sum(), (sc['狀態'] == '有證據').sum()))
print('『下一次實驗結果』『結論』兩欄留白，做完第 7 章的實驗後由學生填寫。')
print('報告範本：ch09_validation_report_template.md（一個驗證項目一個小節，現有證據已預填；實驗後把『結果』『結論』補上）。')
with open('ch09_validation_report_template.md', encoding='utf-8') as fh:
    tpl = fh.read()
print('\n能說什麼：')
print('  - 1 mm、2258 nm 的二階導數與滴定值（n = 4）呈線性關係；這是一個可行性研究（feasibility study）的初步證據。')
print('不能說什麼：')
print('  - 『方法已確效』；『準確度是 x %』（沒有獨立樣品）；『可取代滴定』（參考方法彼此對不上，NIR 本身依賴它們校正）。')

chk('scorecard: number of rows', sc.shape[0], 8)
chk('scorecard: number of columns', sc.shape[1], 6)
chk('scorecard: status 沒有 (CH9.card 狀態 0)', int((sc['狀態'] == '沒有').sum()), 5)
chk('scorecard: status 部分或初步 (CH9.card 狀態 1)', int((sc['狀態'] == '部分或初步').sum()), 3)
chk('scorecard: status 有證據 (CH9.card 狀態 2)', int((sc['狀態'] == '有證據').sum()), 0)
chk('scorecard: per-item status mismatches vs CH9.card',
    int((sc['狀態'].values != np.array(['沒有', '部分或初步', '部分或初步', '沒有', '沒有', '部分或初步', '沒有', '沒有'])).sum()), 0)
chk('scorecard: blank 下一次實驗結果 cells', int((sc['下一次實驗結果'] == '').sum()), 8)
chk('scorecard: blank 結論 cells', int((sc['結論'] == '').sum()), 8)
chk('template: item names found in report template (of 8)', sum(1 for c_ in card if c_[0] in tpl), 8)
show_checks('第 8 節')
''')

md(r'''
### 動手改改看（第 8 節）
1. 用 Excel 開 `ch09_validation_scorecard.csv`，中文有亂碼嗎？（檔案含 BOM，應該不會）
2. 做完第 7 章的實驗後，哪幾項的狀態可以從「沒有」升到「部分或初步」？哪一項仍然升不上去？
3. 為什麼「有證據」一項都沒有，卻還是值得做這個方法？
''')

code(r'''
print('全部自我驗證：%d / %d PASS' % (sum(c[0] for c in checks), len(checks)))
for c in checks:
    if not c[0]:
        print('FAIL:', c[1], c[2], c[3])
''')

md(r'''
### 總結
1. 目標 ±0.5 % 是示意；先有用途與目標，才能判斷證據夠不夠。
2. 兩個醋酸帶一致是一致性不是選擇性；n = 4 的 R² 不能用來選參考方法。
3. 線性：4 點、沒有重複，純誤差自由度 0；下一次 8 × 3 才有 16 與 6。
4. 準確度：沒有獨立樣品；PLS RMSECV 0.32 < Atago 日間 SD 0.41 是過度配適的徵兆；PAL-30S 0.26 是一致性不是準確度。
5. 精密度：只量了掃描雜訊；填樣、配製、日間都沒量。
6. LOD/LOQ：殘差法 1.8 / 6.1 %（滴定）、儀器雜訊法 0.8 / 2.7 %，皆為初步。
7. 預測區間：現在半寬 3.14 %，下一次 0.13–1.3 %；示意目標 0.5 % 只有填樣 SD 夠小、且參考值誤差很小時才可能達到。
8. 確效總表 8 項：沒有 5、部分或初步 3、有證據 0；可說「可行性研究、初步線性」，不能說「已確效」。
''')


def local_head():
    return [
        (MD, '> **本機 Jupyter 版**：請在 `teaching/exercises/` 目錄下開啟，資料讀自 `../data/spectra_long.csv`，報告範本讀自同一資料夾的 `ch09_validation_report_template.md`。\n'
             '> 需要套件：numpy、pandas、matplotlib、scipy。Colab 請改開 `ch09_validation_report_colab.ipynb`。'),
        (CODE, "DATA_PATH = '../data/spectra_long.csv'   # 本機版：讀上一層 data 資料夾"),
    ]


def _embed(path):
    with open(path, 'rb') as fh:
        raw = fh.read()
    b64 = base64.b64encode(gzip.compress(raw, 9, mtime=0)).decode('ascii')
    return '\n'.join(b64[i:i + 100] for i in range(0, len(b64), 100))


def colab_head():
    embed = (
        "# 資料內嵌：把 spectra_long.csv 與報告範本以 gzip + base64 存成字串，執行時解出成檔案\n"
        "import base64, gzip\n"
        "DATA_B64 = \"\"\"\n" + _embed(DATA_CSV) + "\n\"\"\"\n"
        "TEMPLATE_B64 = \"\"\"\n" + _embed(TEMPLATE_MD) + "\n\"\"\"\n"
        "with open('spectra_long.csv', 'wb') as _f:\n"
        "    _f.write(gzip.decompress(base64.b64decode(DATA_B64)))\n"
        "with open('ch09_validation_report_template.md', 'wb') as _f:\n"
        "    _f.write(gzip.decompress(base64.b64decode(TEMPLATE_B64)))\n"
        "DATA_PATH = 'spectra_long.csv'\n"
        "print('已解出 spectra_long.csv：', len(gzip.decompress(base64.b64decode(DATA_B64))), 'bytes；ch09_validation_report_template.md：',\n"
        "      len(gzip.decompress(base64.b64decode(TEMPLATE_B64))), 'bytes')"
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
             '> 2. 選單「**檔案 → 上傳筆記本**」，選擇 `ch09_validation_report_colab.ipynb`。\n'
             '> 3. 依序執行每一格。\n>\n'
             '> 資料與報告範本已**內嵌**在筆記本中，下一格會自動解出 `spectra_long.csv` 與 `ch09_validation_report_template.md`。Colab 已內建所需套件。\n'
             '> 第 8 節產生的 `ch09_validation_scorecard.csv` 會存在 Colab 的檔案面板（左側資料夾圖示），請自行下載。'),
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
    build('ch09_validation_report.ipynb', local_head())
    build('ch09_validation_report_colab.ipynb', colab_head())
