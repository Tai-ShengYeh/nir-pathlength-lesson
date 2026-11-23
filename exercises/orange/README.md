# Orange 練習（第 1–6 章）

用 [Orange Data Mining](https://orangedatamining.com/) 拖拉 widget 完成與 R / Python 練習相同的分析，不必寫程式。每章一個工作流程檔（`.ows`）加一份逐步講義。

> 測試版本：**Orange 3.40.0** + **Spectroscopy 外掛 0.9.2**（Windows 11）。其他版本的介面名稱或預設值可能略有不同。

## 安裝（只要做一次）

1. 到 Orange 官網下載 Windows 版安裝程式並安裝（本講義用 3.40.0 測試）。
2. 開啟 Orange，選單 **Options → Add-ons…**（附加元件）。
3. 勾選 **Spectroscopy**，按 **OK**。安裝完成後依提示**重新啟動 Orange**。
4. 重開後，左側 widget 面板會多出 **Spectroscopy** 類別（**Spectra**、**Preprocess Spectra** 都在裡面）。

> Add-ons 視窗一直轉圈或失敗，通常是網路或防火牆擋住，換個網路再試；學校電腦可能需要系統管理員權限。

## 各章內容

| 章 | 講義 | 工作流程 | 資料 | 做什麼 | Orange 讀出的關鍵數字 |
|---|---|---|---|---|---|
| 1 | [ch01_orange_guide.md](ch01_orange_guide.md) | `ch01_pathlength.ows`（8 個 widget） | `data/spectra_wide_orange.csv` | 看 5 種光徑的飽和平台；1 mm、二階導數、只用 2258 nm 做線性迴歸與留一法；換光徑、換 1670 nm | 1 mm / 2258 nm：R² 0.934，LOO RMSE 2.198 |
| 2 | [ch02_orange_guide.md](ch02_orange_guide.md) | `ch02_references.ows`（8 個 widget） | `data/spectra_1mm_refs_orange.csv` | Formula 算滴定酸度與滴定 / Atago 比值；Select Columns 輪流換 4 種參考值看 R² | R² 0.996 / 0.990 / 0.956 / 0.934（標示 / 滴定 / Atago 11/03 / 11/02，2258 nm）；比值 0.757 / 1.068 / 1.661 |
| 3 | [ch03_orange_guide.md](ch03_orange_guide.md) | `ch03_wavelength_pitfall.ows`（10 個 widget） | `data/spectra_wide_orange.csv` | 88 個可用波長中挑 R² 最高的：在交叉驗證外挑（Rank）vs 在交叉驗證內挑（Preprocess → Preprocessor） | 挑到 1862.4 nm；LOO RMSE 外 0.427 vs 內 2.578 |
| 4 | [ch04_orange_guide.md](ch04_orange_guide.md) | `ch04_pca_pls.ows`（9 個 widget） | `data/spectra_wide_orange.csv` | PCA 看 PC1 抓到光徑；1 mm 二階導數 PLS（1 成分、不縮放）留一法 | PC1 0.917；PLS RMSECV 0.317 |
| 5 | [ch05_orange_guide.md](ch05_orange_guide.md) | `ch05_prediction.ows`（12 個 widget） | `data/spectra_wide_orange.csv`、`data/unknown_2258_orange.csv` | 預測未知樣品（**只有點預測**）：藏起 4.87 % 瓶用 3 瓶預測；4 瓶校正預測 y0 = 4.20 × 10⁻³。區間（Fieller）Orange 做不到，請到 R / Python | 3.669（講義 3.671，真值 4.87）；4.891（講義 4.960） |
| 6 | [ch06_orange_guide.md](ch06_orange_guide.md) | `ch06_data_quality.ows`（20 個 widget） | `data/spectra_wide_orange.csv` | 資料品質與離群值：水帶高度 ÷ 同光徑中位數（±10 % 旗標）、0.25 mm 光譜目視、Outliers / PCA 的陷阱（光徑混在一起）。**沒有 Q 殘差**；積分是梯形法（與講義比值最多差 0.0003） | 旗標 0.25 mm 4.87 % 瓶（0.7725，講義 0.7723）與 0.5 mm 3.32 % 瓶（1.1701，講義 1.17）；Outliers 預設（LOF）抓 **0** 條 |

截圖都在 `img/`（`chNN_orange_*.png`）。

## 開啟工作流程

**File → Open** 選 `.ows`。資料路徑以相對路徑存在 `.ows` 裡（`../../data/…`），只要保持這個結構，整個資料夾搬到哪裡都能開：

```
（任意位置）/
├── data/
│   ├── spectra_wide_orange.csv
│   ├── spectra_1mm_refs_orange.csv
│   └── unknown_2258_orange.csv
└── exercises/orange/
    ├── ch01_pathlength.ows
    ├── ch02_references.ows
    ├── ch03_wavelength_pitfall.ows
    ├── ch04_pca_pls.ows
    ├── ch05_prediction.ows
    └── ch06_data_quality.ows
```

如果 **File** widget 顯示紅色錯誤（找不到檔案），雙擊它重新選資料檔即可。在 Orange 裡**另存**後，資料路徑會變成你電腦上的絕對路徑；要發給別人請用原始的 `.ows`。

## 一定要先知道的兩件事

1. **Orange 的 Linear Regression 是直接迴歸，講義是反向校正。** 講義（與 R / Python 練習）擬合 `訊號 = a + b·濃度` 再反推濃度；Orange 擬合 `濃度 = a + b·訊號`。4 點全用時 **R² 完全相同**，但**留一法的 RMSE 不同**（例如 1 mm / 2258 nm：Orange 2.198，講義 2.182；5 mm / 2258 nm：Orange 4.967，講義 17.16）。第 1–3 章講義列出兩邊的數字，並在第 1 章第 5 節解釋差別。第 4 章的 PLS 兩邊算法相同，數字一致（0.317）。第 5 章的預測值也因此不同（y0 = 4.20：Orange 4.891，講義 4.960），而且 Orange 沒有預測區間，見第 5 章講義第 4、5 節。
2. **Orange 的留一法依檔案的列順序**（`spectra_wide_orange.csv` 是樣品名稱排序：10 %、3 %、5 %、水），R / Python 依濃度由小到大。逐筆比較預測值時，請用濃度對照，不要用列號。（第 2 章的資料檔本來就照濃度排列。）

另外：Test and Score 在 **Leave one out** 下顯示的 R2 是「留一預測值」的 R²，不是校正 R²；校正 R² 要切到 **Test on train data**。

## 給老師：重新產生與驗證

所有 `.ows` 都由 Orange 自己的 scheme API 產生（不是手寫 XML），並以 headless 方式重新載入、直接讀 widget 上的數字驗證：

```
python tools/nir-pathlength/make_orange_refs.py          # 產生 data/spectra_1mm_refs_orange.csv（第 2 章）
python tools/nir-pathlength/make_orange_unknown.py       # 產生 data/unknown_2258_orange.csv（第 5 章）
python exercises/orange/tools/build_ch01_ows.py          # 產生 .ows（ch01 – ch06 同理）
python exercises/orange/tools/verify_ch01_ows.py --img exercises/orange/img   # 重新載入、印出數字、存截圖
```

- `tools/owsutil.py`、`tools/owsbuild.py`：共用的載入、等待、讀表、存檔與路徑檢查工具。
- 每次存檔都會解碼 `.ows` 裡所有設定（含 context 設定的內容），確認沒有任何本機絕對路徑。
- 驗證腳本另以 numpy 重做同一套直接迴歸，並列出講義的反向校正數字對照。
