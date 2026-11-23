# 打包每章的教學包（zip），並寫出 portal/packs.js 給課程入口頁顯示檔案大小。
# 用法（在 teaching/ 底下）：python tools/build_packs.py
# 教學包內的資料夾結構與 repo 相同（exercises/、exercises/orange/、data/），R 與 Jupyter 的 ../data 路徑、Orange 的相對路徑都能直接用。
import os, glob, json, zipfile, io

ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
SITE = 'https://tai-shengyeh.github.io/nir-pathlength-lesson/'
COLAB = 'https://colab.research.google.com/github/Tai-ShengYeh/nir-pathlength-lesson/blob/main/'
CH = [  # (編號, 檔名代稱, 標題, 練習檔前綴)
    ('01', 'pathlength', '最佳光徑選擇', 'exercises/pathlength_exercise'),
    ('02', 'references', '參考值與實驗設計', 'exercises/ch02_reference_design'),
    ('03', 'wavelength_pitfall', '波長挑選的陷阱', 'exercises/ch03_wavelength_pitfall'),
    ('04', 'pca_pls', 'PCA 與 PLS', 'exercises/ch04_pca_pls'),
    ('05', 'prediction', '預測未知樣品與不確定度', 'exercises/ch05_prediction_uncertainty'),
    ('06', 'data_quality', '資料品質與離群值', 'exercises/ch06_data_quality'),
    ('07', 'experiment_design', '設計下一次實驗', 'exercises/ch07_experiment_design'),
]

def files_for(n, prefix):
    fs = [prefix + '.R', prefix + '.ipynb', prefix + '_colab.ipynb']
    fs += sorted(glob.glob(f'exercises/ch{n}_*.csv') + glob.glob(f'exercises/ch{n}_*.md'))   # 第 7 章的順序表與紀錄表
    fs += sorted(glob.glob(f'exercises/orange/ch{n}_*'))                      # .ows 與操作指南
    fs += sorted(glob.glob(f'exercises/orange/img/ch{n}_*'))                  # 指南截圖
    return [f.replace('\\', '/') for f in fs if os.path.isfile(f)]

def readme(n, title, fs):
    has_orange = any(f.endswith('.ows') for f in fs)
    lines = [f'近紅外光譜：醋酸定量與最佳光徑　第 {int(n)} 章　{title}', '=' * 40, '',
             f'動畫課程（有旁白與小測驗）：{SITE}site/nir-pathlength/' + ('' if n == '00' else f'ch{n}/'), '',
             '這個教學包裡有：']
    lines += ['  ' + f for f in fs]
    lines += ['  data/            光譜資料（所有練習共用）', '', '怎麼開始：',
              '1. 先看動畫課程，再做練習。',
              '2. R：用 RStudio 開啟 exercises/ 裡的 .R 檔，按 Source；最後會列出自我驗證 PASS / FAIL。',
              '3. Python：用 Jupyter 開啟 exercises/ 裡的 .ipynb（不是 _colab 那本），選「Run All」。',
              '   或用 Google Colab：上傳 _colab.ipynb（資料已內嵌），選「執行階段 → 全部執行」。']
    if has_orange:
        lines += ['4. Orange：先安裝 Orange Data Mining 與 Spectroscopy 外掛，開啟 exercises/orange/ 裡的 .ows 檔，',
                  '   照同資料夾的 _guide.md 操作（請保持資料夾結構，不要只拿出 .ows 檔）。']
    lines += ['', '全部檢查都 PASS，才算完成。解壓縮時請保留資料夾結構。']
    return ('\r\n'.join(lines) + '\r\n').encode('utf-8-sig')

def build(name, items):
    path = os.path.join(ROOT, 'packs', name)
    with zipfile.ZipFile(path, 'w', zipfile.ZIP_DEFLATED) as z:
        for arc, src in items:
            if isinstance(src, str): z.write(os.path.join(ROOT, src), arc)
            else: z.writestr(arc, src)
    return os.path.getsize(path)

os.chdir(ROOT); os.makedirs('packs', exist_ok=True)
data = sorted(f.replace('\\', '/') for f in glob.glob('data/*.csv'))
top = 'nir-pathlength'
info, everything = {}, set()
for n, slug, title, prefix in CH:
    fs = files_for(n, prefix); everything.update(fs)
    items = [(f'{top}/{f}', f) for f in fs + data] + [(f'{top}/README.txt', readme(n, title, fs))]
    name = f'nir_ch{n}_{slug}.zip'; size = build(name, items)
    ows = next((f for f in fs if f.endswith('.ows')), None)
    guide = next((f for f in fs if f.endswith('_guide.md')), None)
    info['ch' + n] = {'zip': 'packs/' + name, 'size': size, 'ows': ows, 'guide': guide,
                      'r': prefix + '.R', 'ipynb': prefix + '.ipynb', 'colab': COLAB + prefix + '_colab.ipynb'}
    print(name, round(size / 1024), 'KB', len(fs), 'files')
allf = sorted(everything) + data + sorted(glob.glob('exercises/orange/README.md')) + ['最佳光徑選擇_講義.md']
items = [(f'{top}/{f}', f) for f in allf] + [(f'{top}/README.txt', readme('00', '全部章節', sorted(everything)).replace('第 0 章　'.encode(), b''))]
size = build('nir_all_chapters.zip', items); info['all'] = {'zip': 'packs/nir_all_chapters.zip', 'size': size}
print('nir_all_chapters.zip', round(size / 1024), 'KB')
open('portal/packs.js', 'w', encoding='utf-8').write('// 由 tools/build_packs.py 產生\nconst PACKS = ' + json.dumps(info, ensure_ascii=False, indent=1) + ';\n')
