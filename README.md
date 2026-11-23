# 近紅外光譜：最佳光徑選擇

以醋酸的近紅外（NIR）實測光譜為例，教專題研究學生如何選擇最佳光徑。課程會重新教一次二階導數（Savitzky–Golay）、線性回歸、偵測極限（LOD）與留一交叉驗證（LOOCV）。

| 教材 | 連結 |
|---|---|
| 🎬 動畫課程（第 1 章 9 分鐘、第 2、3 章各 7 分鐘，有旁白與小測驗） | https://tai-shengyeh.github.io/nir-pathlength-lesson/ |
| [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Tai-ShengYeh/nir-pathlength-lesson/blob/main/exercises/pathlength_exercise_colab.ipynb) | Python 練習，一鍵在 Colab 開啟（資料已內嵌） |
| Jupyter 練習 | [exercises/pathlength_exercise.ipynb](exercises/pathlength_exercise.ipynb) |
| R 練習 | [exercises/pathlength_exercise.R](exercises/pathlength_exercise.R)（只用 base R） |
| 第 2 章練習：參考值與實驗設計 | [R](exercises/ch02_reference_design.R) · [Jupyter](exercises/ch02_reference_design.ipynb) · [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Tai-ShengYeh/nir-pathlength-lesson/blob/main/exercises/ch02_reference_design_colab.ipynb) |
| 第 3 章練習：波長挑選的陷阱 | [R](exercises/ch03_wavelength_pitfall.R) · [Jupyter](exercises/ch03_wavelength_pitfall.ipynb) · [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Tai-ShengYeh/nir-pathlength-lesson/blob/main/exercises/ch03_wavelength_pitfall_colab.ipynb) |
| 學生講義 | [最佳光徑選擇_講義.md](最佳光徑選擇_講義.md) |
| 資料 | [data/spectra_long.csv](data/spectra_long.csv)（21 條光譜 × 142 波長） |

**結論**：在靈敏度與光徑仍成正比（比爾–朗伯定律成立）的範圍內，取最長的光徑，也就是 1 mm。

練習的最後一節會自動比對答案鍵，R 與 Python 版各有 87 項檢查，全部 PASS 才算完成。

## 本機播放動畫課程

```
python tools/serve.py 8765 site
```

然後開啟 http://localhost:8765/nir-pathlength/

## 資料夾說明

- `site/nir-pathlength/`：動畫課程（靜態網頁）
- `exercises/`：R 與 Python 練習；`build_notebooks.py` 由同一份內容產生 Jupyter 與 Colab 兩本筆記本
- `data/`：整理後的光譜資料
- `content/`、`books/`、`tools/`：課程的旁白稿、設定與製作工具

## 授權與致謝

動畫課程以 [Papermorph](https://github.com/DozenTwelve/Papermorph) 製作，播放器與繪圖引擎（`site/nir-pathlength/lib/` 中的 engine、contents、units、volume、math 檔案）與 `tools/serve.py` 依 MIT 授權使用，詳見 [THIRD_PARTY.md](THIRD_PARTY.md)。旁白由 Microsoft Edge TTS（zh-TW-HsiaoChenNeural）產生。
