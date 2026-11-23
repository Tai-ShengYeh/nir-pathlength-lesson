# 第 2 章 Orange 練習用的小資料：1 mm、重複 1 的 4 條光譜（水 + 3 瓶醋）＋ 各種參考值。
# 寫成 data/spectra_1mm_refs_orange.csv（Orange 欄名旗標格式，與 spectra_wide_orange.csv 相同）。
#   python tools/nir-pathlength/make_orange_refs.py
# 吸光度直接取自 data/spectra_long.csv（字串原樣複製，不重新格式化）；參考值與
# site/nir-pathlength/ch02/index.html 的常數相同：
#   Atago 11/02 = 0, 3.32, 4.87, 7.45；Atago 11/03 = 0, 3.01, 4.12, 6.86；標示 = 0, 3, 5, 10
#   滴定 NaOH 體積 V = 8.8, 17, 44 mL（水 = 0）；Orange 裡再用 Formula 算 酸度 = V × 1.07875 × 0.24
# 所有參考值都放在 meta（mC#），沒有預設 target：學生在 Select Columns 自己指定。
# 列順序 = 濃度由小到大（水、3 %、5 %、10 % 瓶），所以 Orange 的 LOO 順序與 R / Python 相同。
import os, csv

here = os.path.dirname(os.path.abspath(__file__))
root = os.path.join(here, '..', '..')
SAMPLES = ['water', 'vinegar_3%', 'vinegar_5%', 'vinegar_10%']
REF = {  # sample: (label_pct, atago_1102, atago_1103, naoh_ml)
    'water':       ('0', '0', '0', '0'),
    'vinegar_3%':  ('3', '3.32', '3.01', '8.8'),
    'vinegar_5%':  ('5', '4.87', '4.12', '17'),
    'vinegar_10%': ('10', '7.45', '6.86', '44'),
}

spec, wls = {}, []
with open(os.path.join(root, 'data', 'spectra_long.csv'), newline='', encoding='utf-8') as fh:
    for r in csv.DictReader(fh):
        if float(r['path_mm']) != 1.0 or r['replicate'] != '1':
            continue
        assert r['sample'] in REF and float(r['acid_pct']) == float(REF[r['sample']][1])
        assert float(r['label_pct']) == float(REF[r['sample']][0])
        spec.setdefault(r['sample'], []).append((float(r['wavelength_nm']), r['wavelength_nm'], r['absorbance']))
for s in SAMPLES:
    spec[s].sort()
    assert len(spec[s]) == 142
wl_names = ['%.3f' % w for w, _, _ in spec[SAMPLES[0]]]          # 與 spectra_wide_orange.csv 欄名相同
assert all([w for w, _, _ in spec[s]] == [w for w, _, _ in spec[SAMPLES[0]]] for s in SAMPLES)

out = os.path.join(root, 'data', 'spectra_1mm_refs_orange.csv')
with open(out, 'w', newline='', encoding='utf-8') as fh:
    wr = csv.writer(fh)
    wr.writerow(['mS#sample', 'mC#label_pct', 'mC#atago_1102', 'mC#atago_1103', 'mC#naoh_ml'] + wl_names)
    for s in SAMPLES:
        wr.writerow([s, *REF[s]] + [a for _, _, a in spec[s]])
print('wrote', os.path.normpath(out), '| 4 rows x', len(wl_names), 'wavelengths')

# 交叉檢查：與 spectra_wide_orange.csv 的 1 mm 列逐值相同
with open(os.path.join(root, 'data', 'spectra_wide_orange.csv'), newline='', encoding='utf-8') as fh:
    rows = list(csv.reader(fh))
hdr = rows[0]
assert hdr[4:] == wl_names
for r in rows[1:]:
    if float(r[0]) == 1.0 and r[3] == '1':
        assert [float(v) for v in r[4:]] == [float(a) for _, _, a in spec[r[1]]], r[1]
print('identical to the 1 mm rows of spectra_wide_orange.csv')
