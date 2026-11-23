# 第 9 章「方法確效報告」用的數字：寫成 site/nir-pathlength/lib/ch09data.js
# 大部分數字直接取自前面各章的資料檔（nirdata.js、ch03–ch08data.js），保證和各章一致；
# 1670 nm 的滴定 R² 在這裡重算（同第 2 章、報告的方法）。確定性計算，無亂數。
import os, json
import numpy as np, pandas as pd
from scipy.signal import savgol_filter

here = os.path.dirname(os.path.abspath(__file__)); root = os.path.join(here, '..', '..')
lib = os.path.join(root, 'site', 'nir-pathlength', 'lib')
def load(f):
    t = open(os.path.join(lib, f), encoding='utf-8').read(); return json.loads(t[t.index('{'):t.rindex('}') + 1])
nir, c3, c4, c5, c6, c7, c8 = (load(f) for f in ('nirdata.js', 'ch03data.js', 'ch04data.js', 'ch05data.js', 'ch06data.js', 'ch07data.js', 'ch08data.js'))

# 1 mm、二階導數（同前面各章）；兩個醋酸帶對三種參考值的 R²
d = pd.read_csv(os.path.join(root, 'data', 'spectra_long.csv')); d = d[(d.replicate == 1) & (d.path_mm == 1.0)]
g = d.pivot_table(index='acid_pct', columns='wavelength_nm', values='absorbance').sort_index()
wl = g.columns.values; D = savgol_filter(g.values, 9, 2, deriv=2, axis=1)
refs = {'atago': np.array([0, 3.32, 4.87, 7.45]), 'titr': np.array([0, *c8['refs']['titr']]), 'label': np.array([0, 3.0, 5.0, 10.0]),
        'pal30s': np.array([0, *c8['refs']['pal30s']])}   # Atago = PAL-Easy ACID 181（導電度式）；PAL-30S 折射計 11/17
r2 = {str(w0): {k: float(np.corrcoef(v, D[:, int(np.argmin(abs(wl - w0)))])[0, 1] ** 2) for k, v in refs.items()} for w0 in (2258, 1670)}
A2258 = g.values[:, int(np.argmin(abs(wl - 2258)))]

out = {
  'target_hw': 0.5,                                   # 示意的用途目標：預測區間 ±0.5 %
  'spec': {'r2': r2, 'n': 4},
  'range': {'atago': [0, 7.45], 'titr': [0, round(float(refs['titr'][-1]), 2)], 'plan': [0, 7], 'A2258_1mm_max': float(A2258.max()), 'A_limit': 1.5},
  'lin': {'r2_atago': c6['one']['r2_atago'], 'r2_titr': c6['one']['r2_titr'], 'pure_df_now': c7['spacing']['now']['pure_df'],
          'lof_df_plan': c7['spacing']['eight']['lof_df'], 'pure_df_plan': c7['spacing']['eight']['pure_df']},
  'acc': {'rmsecv_atago': c6['one']['rmse_atago'], 'rmsecv_titr': c6['one']['rmse_titr'], 'rmsecv_pls': c4['rmsecv']['pls1'], 'rmsecv_pal30s': c8['nirfit']['pal30s']['rmsecv'],
          'atago_day_sd': c8['tonir']['atago_day_sd_abs'], 'n_test': 0, 'ref_gap_max_uV': max(abs(v) for k in c8['agree'] for v in c8['agree'][k]['gap_in_uV'])},
  'prec': {'scan_rms_AU': c6['rep']['rms'], 'scan_conc_025': c7['scan']['conc025'], 'scan_conc_1mm_est': c7['scan']['conc1_est'],
           'fill_measured': False, 'resid_titr': c7['scan']['resid_titr'], 'atago_day_sd': c8['tonir']['atago_day_sd_abs']},
  # 第 1 章的儀器雜訊法：平坦區 1550–1650 nm 的二階導數雜訊 σ(1 mm) = 2.7e-4（ch01/index.html），斜率取 nirdata 的 1 mm 2258 nm
  'lod': {'res_titr': c5['lod_res']['titr'], 'res_atago': c5['lod_res']['atago'],
          'ch01_sigma_1mm': 2.7e-4, 'ch01_slope_1mm': nir['band']['2258']['1']['slope'],
          'ch01_lod_1mm': 3 * 2.7e-4 / abs(nir['band']['2258']['1']['slope']), 'ch01_loq_1mm': 10 * 2.7e-4 / abs(nir['band']['2258']['1']['slope'])},
  'suit': {'wb_band': 0.10, 'flags': 2, 'wb_peak_1mm': c6['wb']['1.0']['peak'], 'wb_usable_1mm': c6['wb']['1.0']['usable']},
  'pred': {'now_lo': c5['four_titr']['pred']['lo'], 'now_hi': c5['four_titr']['pred']['hi'],
           'now_hw': (c5['four_titr']['pred']['hi'] - c5['four_titr']['pred']['lo']) / 2,
           'now_x0': c5['four_titr']['pred']['x0'], 'now_m': 1,
           # 下一次實驗（8 濃度 × 3 次、未知樣品填 3 次、4 %）：寬度取決於預試驗量到的填樣誤差 0.1–1.0 %（第 7 章）
           'plan_s': c7['sens']['s'], 'plan_hw_m3': c7['sens']['m3'],
           'plan_hw_lv': c7['sens_ref']['hw'], 'plan_s_ref': c7['sens_ref']['s_ref'], 'plan_s_fill': c7['sens_ref']['s_fill']},
}
# 總表：每一項的狀態（2 有證據、1 部分／初步、0 沒有）與補的方式
out['card'] = [
  ['選擇性', 0, '兩個醋酸帶一致只是一致性；沒有測過糖等干擾物', '市售盲樣、加入干擾物（第 7 章）'],
  ['線性', 1, '4 點、沒有重複，無法做缺適性檢定', '8 濃度 × 3 次（第 7 章）'],
  ['範圍', 1, '2258 nm 在 1 mm 未飽和；11.4 % 只靠一次有爭議的滴定', '0–7 %、8 個等間隔濃度'],
  ['準確度', 0, '沒有獨立樣品；參考值彼此對不上', '盲樣＋標定過的滴定（第 7、8 章）'],
  ['精密度', 0, '只有掃描雜訊；填樣、日間誤差沒量', '預試驗填樣 10 次（第 7 章）'],
  ['LOD／LOQ', 1, '殘差法 1.8 / 6.1 %、儀器雜訊法 0.8 %（初步）', '新校正重算'],
  ['穩健性', 0, '溫度、操作者、樣品槽都沒測', '刻意改變條件的小實驗'],
  ['系統適用性', 0, '水帶規則在 1 mm 只能當參考（水帶飽和）', '每次實驗量共同對照（第 7 章）'],
]

def rnd(o):
    if isinstance(o, float): return float('%.4g' % o)
    if isinstance(o, dict): return {k: rnd(v) for k, v in o.items()}
    if isinstance(o, list): return [rnd(v) for v in o]
    return o
out = rnd(out)
open(os.path.join(lib, 'ch09data.js'), 'w', encoding='utf-8').write(
    '// 由 tools/nir-pathlength/make_data_ch09.py 產生（方法確效報告）\nconst CH9 = ' + json.dumps(out, ensure_ascii=False) + ';\n')
print(json.dumps({k: v for k, v in out.items() if k != 'card'}, ensure_ascii=False, indent=1))
