# =====================================================================
# NIR 醋酸最佳光徑選擇：R 程式練習（只用 base R：stats、graphics、utils）
#
# 執行方式：
#   Rscript pathlength_exercise.R          （在任何目錄都可以；程式會自動切到本檔所在資料夾）
#   或在 RStudio 開啟後，先把工作目錄設到 teaching/exercises，再逐段執行。
#
# 資料：../data/spectra_long.csv     圖：figures_R/*.png（圖上文字一律用英文）
#
# 每一節的節奏：小例子 -> 手刻 -> 套到真實資料 -> 用 R 內建函數對照。
# 每節結尾有「動手改改看」，請改參數後重跑。
# =====================================================================

options(width = 130, digits = 7)

# ---- 切換到本檔所在資料夾（Rscript 才有 --file= 參數）----
.args <- commandArgs(trailingOnly = FALSE)
.fa <- sub("^--file=", "", .args[grep("^--file=", .args)])
if (length(.fa) == 1 && nzchar(.fa)) setwd(dirname(normalizePath(.fa)))
DATA_PATH <- "../data/spectra_long.csv"
FIG_DIR <- "figures_R"
dir.create(FIG_DIR, showWarnings = FALSE)

hr <- function(title) cat("\n", strrep("=", 70), "\n", title, "\n", strrep("=", 70), "\n", sep = "")
fmt <- function(x, d = 4) formatC(x, digits = d, format = "g")

# ---- 可調整的參數（練習題會請你改這幾個值）----
A_LIM      <- 1.5   # 保守的線性吸光度上限
A_SAT      <- 2.0   # 飽和判斷值
SG_WIN     <- 9     # 二階導數的 SG 視窗點數（奇數）
SMOOTH_WIN <- 7     # 估計雜訊用的 SG 平滑視窗
RATIO_TOL  <- 0.2   # 第 7 節：靈敏度比與理想值的容許偏差
BANDS <- c("2258" = 2258, "1670" = 1670)   # 醋酸特徵帶 (nm)：C-H 組合帶、C-H 一倍頻
PATHS <- c(0.25, 0.5, 1, 2, 5)
ACID  <- c(0, 3.32, 4.87, 7.45)            # Atago 實測酸度 (%)
LABEL <- c("water 0%", "3.32%", "4.87%", "7.45%")
COLS  <- c("black", "blue", "darkgreen", "red")
pk <- function(p) as.character(p)          # 光徑 -> list 的鍵

# =====================================================================
# 0. 載入資料、畫出 5 個光徑的光譜
# =====================================================================
hr("0. 載入資料與光譜圖")

df <- read.csv(DATA_PATH)
cat("資料維度：", nrow(df), "列 x", ncol(df), "欄\n")
print(head(df, 3))

main <- df[df$replicate == 1, ]                 # 主分析只用 replicate == 1
wl <- sort(unique(main$wavelength_nm))
n_wl <- length(wl)

spec <- list()                                  # spec[["1"]] -> 4 x 142 矩陣（列：水, 3.32, 4.87, 7.45）
for (p in PATHS) {
  S <- NULL
  for (a in ACID) {
    h <- main[abs(main$path_mm - p) < 1e-9 & abs(main$acid_pct - a) < 1e-9, ]
    h <- h[order(h$wavelength_nm), ]
    stopifnot(nrow(h) == n_wl, all(abs(h$wavelength_nm - wl) < 1e-6))
    S <- rbind(S, h$absorbance)
  }
  spec[[pk(p)]] <- S
}
cat(sprintf("波長範圍 %.1f - %.1f nm，共 %d 點；spec[['1']] 維度 = %s\n",
            wl[1], wl[n_wl], n_wl, paste(dim(spec[["1"]]), collapse = " x ")))

png(file.path(FIG_DIR, "fig01_spectra.png"), width = 1800, height = 380, res = 110)
par(mfrow = c(1, 5), mar = c(4, 4, 2, 0.5))
for (p in PATHS) {
  S <- spec[[pk(p)]]
  plot(wl, S[1, ], type = "n", ylim = range(S), xlab = "wavelength (nm)",
       ylab = if (p == PATHS[1]) "absorbance (AU)" else "", main = paste0(p, " mm"))
  for (i in 1:4) lines(wl, S[i, ], col = COLS[i])
  abline(h = 1.5, lty = 2, col = "gray40")
  if (p == PATHS[1]) legend("topright", legend = LABEL, col = COLS, lty = 1, cex = 0.7)
}
dev.off()

# ---- 動手改改看（第 0 節）----
# 1. 把 ACID 的 3.32 改成 3.0，哪裡會出錯？為什麼比較浮點數要用 abs(x - a) < 1e-9？
# 2. 只畫 spec[["1"]] 的 7.45% 那條，xlim 限制在 2100-2300 nm。

# =====================================================================
# 1. 基礎①：Beer-Lambert 與飽和
# =====================================================================
hr("1. 基礎①：Beer-Lambert 與飽和  (A = epsilon * b * c)")

# 玩具例子：每 mm 吸 2.0 AU 的樣品，儀器上限 2.3 AU
b <- PATHS
A_true <- 2.0 * b                  # 理想：正比於光徑
A_obs  <- pmin(A_true, 2.3)        # 儀器：超過 2.3 就削平
print(data.frame(path_mm = b, A_true = A_true, A_obs = A_obs), row.names = FALSE)

frac_below <- function(S, thr) mean(apply(S, 2, max) < thr)   # 4 條光譜在該波長「全部」< thr

tab1 <- data.frame(path_mm = PATHS,
                   max_A = sapply(PATHS, function(p) max(spec[[pk(p)]])),
                   frac_lim = sapply(PATHS, function(p) frac_below(spec[[pk(p)]], A_LIM)),
                   frac_sat = sapply(PATHS, function(p) frac_below(spec[[pk(p)]], A_SAT)))
names(tab1)[3:4] <- c(paste0("frac_A<", A_LIM), paste0("frac_A<", A_SAT))
print(round(tab1, 4), row.names = FALSE, check.names = FALSE)

png(file.path(FIG_DIR, "fig02_saturation.png"), width = 1100, height = 420, res = 110)
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
plot(PATHS, A_true, log = "x", type = "l", lty = 3, ylim = c(0, 4), xlab = "path length (mm)", ylab = "absorbance")
lines(PATHS, A_obs, type = "b", col = "red", pch = 15)
lines(PATHS, sapply(PATHS, function(p) max(spec[[pk(p)]])), type = "b", col = "blue", pch = 16)
legend("topleft", c("toy: ideal", "toy: instrument", "real data: max A"), lty = c(3, 1, 1),
       col = c("black", "red", "blue"), pch = c(NA, 15, 16), cex = 0.8)
plot(PATHS, 100 * sapply(PATHS, function(p) frac_below(spec[[pk(p)]], A_LIM)), type = "b", pch = 16,
     ylim = c(0, 100), xlab = "path length (mm)", ylab = "% of wavelengths")
lines(PATHS, 100 * sapply(PATHS, function(p) frac_below(spec[[pk(p)]], A_SAT)), type = "b", pch = 15, col = "red")
legend("topright", c(paste("A <", A_LIM), paste("A <", A_SAT)), col = c("black", "red"), pch = c(16, 15), lty = 1, cex = 0.8)
dev.off()

# ---- 動手改改看（第 1 節）----
# 1. 把 A_LIM 改成 2.0 重跑，哪個光徑的比例變化最大？為什麼？
# 2. 玩具例子把儀器上限改成 1.5，哪幾個光徑被削平？
# 3. 為什麼要求「4 條光譜在該波長都 < 門檻」，而不是只看 7.45% 那一條？

# =====================================================================
# 2. 基礎②：二階導數與 Savitzky-Golay
# =====================================================================
hr("2. 基礎②：二階導數與 Savitzky-Golay")

# ---- 2.1 玩具例子：二階導數消除基線偏移與斜率（有限差分，點間距 = 1）----
x <- -10:10; a_true <- 0.02
y_plain  <- a_true * x^2
y_offset <- a_true * x^2 + 0.5
y_tilt   <- a_true * x^2 + 0.5 + 0.03 * x
fd2 <- function(y) { n <- length(y); y[1:(n - 2)] - 2 * y[2:(n - 1)] + y[3:n] }   # y[i-1] - 2 y[i] + y[i+1]
cat("前 3 個二階差分：\n")
cat("  plain         :", round(fd2(y_plain)[1:3], 6), "\n")
cat("  +offset       :", round(fd2(y_offset)[1:3], 6), "\n")
cat("  +offset+slope :", round(fd2(y_tilt)[1:3], 6), "\n")
cat("理論值 2a =", 2 * a_true, "\n")

png(file.path(FIG_DIR, "fig03_toy_second_derivative.png"), width = 1000, height = 400, res = 110)
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
matplot(x, cbind(y_plain, y_offset, y_tilt), type = "b", pch = 16, cex = 0.5, lty = 1, xlab = "index", ylab = "y", main = "toy signal")
matplot(x[2:20], cbind(fd2(y_plain), fd2(y_offset), fd2(y_tilt)), type = "b", pch = 16, cex = 0.5, lty = 1,
        ylim = c(0, 0.06), xlab = "index", ylab = "2nd difference", main = "2nd difference (all identical)")
dev.off()

# ---- 2.2 SG = 視窗內擬合二次式，二階導數 = 2*a2 ----
set.seed(0)
y_noisy <- y_tilt + rnorm(length(y_tilt), 0, 0.002)
w <- 9; h <- (w - 1) / 2
xs <- -h:h                                       # 視窗內索引 -4..4（點間距 = 1）
center <- 11                                     # R 從 1 開始；第 11 點是 x = 0
seg <- y_noisy[(center - h):(center + h)]
fit_w <- lm(seg ~ xs + I(xs^2))
cat(sprintf("方法 A：lm 擬合二次式，二階導數 = 2*a2 = %.6f\n", 2 * coef(fit_w)[3]))

sg_coeffs <- function(window, deriv) {
  # 設計矩陣 X = [1, x, x^2]，最小平方解 = (X'X)^-1 X' y，第 k 列 = 係數 a_k 的權重
  hh <- (window - 1) / 2
  xx <- -hh:hh
  X <- cbind(1, xx, xx^2)
  P <- solve(t(X) %*% X) %*% t(X)               # 3 x window
  if (deriv == 0) return(as.numeric(P[1, ]))    # 平滑：a0
  if (deriv == 2) return(as.numeric(2 * P[3, ]))# 二階導數：2*a2
  stop("只實作 deriv = 0 或 2")
}
c2 <- sg_coeffs(9, 2)
cat("9 點二次二階導數係數 c2 =", round(c2, 5), "\n")
cat(sprintf("方法 B：c2 與視窗 9 點相乘加總 = %.6f\n", sum(c2 * seg)))
cat("c2 總和 =", format(sum(c2), digits = 3), "（為 0，所以常數偏移被消掉）\n")

# ---- 2.3 套到整條光譜（手刻 SG；以索引為單位、點間距 = 1，與 scipy savgol_filter 同慣例）----
# 單位陷阱：波長網格不等距（間距約 4.5-15.9 nm），所以這裡的 d2A/dlambda2 其實是「每點^2」，
# 不是真正的「每 nm^2」。整份分析都用同一個索引單位，比較光徑沒問題，但數值不能直接當 AU/nm^2。
sg_apply <- function(S, coef) {
  half <- (length(coef) - 1) / 2
  out <- matrix(NA_real_, nrow(S), ncol(S))     # 視窗不足的頭尾給 NA
  for (i in (half + 1):(ncol(S) - half)) out[, i] <- S[, (i - half):(i + half), drop = FALSE] %*% coef
  out
}
c2 <- sg_coeffs(SG_WIN, 2)
cat(sprintf("波長間距 (nm)：最小 %.2f，最大 %.2f（不等距）\n", min(diff(wl)), max(diff(wl))))
D2 <- lapply(spec, sg_apply, coef = c2)

# R 沒有內建 SG：用 lm 在每個 9 點視窗擬合二次式取 2*a2 來對照（只對 1 mm、水那一條）
half <- (SG_WIN - 1) / 2
xs_w <- -half:half
d2_lm <- rep(NA_real_, n_wl)
for (i in (half + 1):(n_wl - half)) {
  yy <- spec[["1"]][1, (i - half):(i + half)]
  d2_lm[i] <- 2 * coef(lm(yy ~ xs_w + I(xs_w^2)))[3]
}
cat(sprintf("手刻 SG 與「lm 逐視窗擬合」最大差異 = %.2e\n", maxdiff <- max(abs(D2[["1"]][1, ] - d2_lm), na.rm = TRUE)))

band_index <- function(w0) which.min(abs(wl - w0))           # 1-based 索引
IDX <- sapply(BANDS, band_index)
for (k in names(BANDS)) cat(sprintf("特徵帶 %s nm -> 索引 %d (實際波長 %.2f nm)\n", k, IDX[[k]], wl[IDX[[k]]]))

png(file.path(FIG_DIR, "fig04_bands_d2_1mm.png"), width = 1100, height = 800, res = 110)
par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
for (k in names(BANDS)) {
  i <- IDX[[k]]; m <- wl > wl[i] - 150 & wl < wl[i] + 150
  matplot(wl[m], t(spec[["1"]][, m]), type = "l", lty = 1, col = COLS, xlab = "wavelength (nm)", ylab = "absorbance",
          main = paste("1 mm absorbance near", k, "nm")); abline(v = wl[i], col = "gold", lwd = 3)
  matplot(wl[m], t(D2[["1"]][, m]), type = "b", pch = 16, cex = 0.5, lty = 1, col = COLS, xlab = "wavelength (nm)",
          ylab = "2nd derivative", main = paste("1 mm 2nd derivative near", k, "nm")); abline(v = wl[i], col = "gold", lwd = 3)
}
dev.off()
cat("1 mm 在 2258 nm 帶的二階導數值（水, 3.32, 4.87, 7.45）：", fmt(D2[["1"]][, IDX[["2258"]]], 5), "\n")

# ---- 動手改改看（第 2 節）----
# 1. 把 SG_WIN 改成 15 重跑。二階導數曲線更平滑還是更雜？最後的自我驗證有些項目會 FAIL，為什麼？
# 2. 玩具例子把 y_tilt 再加上 0.01*x^3（三次項），二階導數還會相同嗎？
# 3. 為什麼 SG 比有限差分「抗雜訊」？比較 sum(c2^2) 與有限差分 [1,-2,1] 的 6。

# =====================================================================
# 3. 基礎③：線性迴歸
# =====================================================================
hr("3. 基礎③：線性迴歸（x = acid_pct，y = 特徵帶二階導數，n = 4）")

fit_line <- function(x, y) {
  sxx <- sum((x - mean(x))^2)
  sxy <- sum((x - mean(x)) * (y - mean(y)))
  slope <- sxy / sxx
  intercept <- mean(y) - slope * mean(x)
  pred <- intercept + slope * x
  r2 <- 1 - sum((y - pred)^2) / sum((y - mean(y))^2)
  list(slope = slope, intercept = intercept, r2 = r2)
}

xt <- c(1, 2, 3, 4); yt <- c(2.1, 3.9, 6.2, 7.8)
ft <- fit_line(xt, yt)
cat(sprintf("手刻：斜率 %.4f 截距 %.4f R2 %.4f\n", ft$slope, ft$intercept, ft$r2))
lt <- lm(yt ~ xt)
cat(sprintf("lm()：斜率 %.4f 截距 %.4f R2 %.4f\n", coef(lt)[2], coef(lt)[1], summary(lt)$r.squared))

CAL <- list()
for (p in PATHS) for (k in names(BANDS)) {
  sig <- D2[[pk(p)]][, IDX[[k]]]                 # 4 個二階導數值
  f <- fit_line(ACID, sig)
  CAL[[paste(p, k)]] <- list(sig = sig, slope = f$slope, intercept = f$intercept, r2 = f$r2)
}
getc <- function(p, k) CAL[[paste(p, k)]]
tab3 <- do.call(rbind, lapply(names(BANDS), function(k) data.frame(
  path_mm = PATHS, band_nm = k,
  slope = sapply(PATHS, function(p) getc(p, k)$slope),
  intercept = sapply(PATHS, function(p) getc(p, k)$intercept),
  R2 = sapply(PATHS, function(p) getc(p, k)$r2))))
print(tab3, row.names = FALSE, digits = 4)
lm1 <- lm(getc(1, "2258")$sig ~ ACID)
cat(sprintf("\n1 mm / 2258 nm：手刻斜率 %.6e，lm() 斜率 %.6e\n", getc(1, "2258")$slope, coef(lm1)[2]))

png(file.path(FIG_DIR, "fig05_calibration.png"), width = 1100, height = 440, res = 110)
par(mfrow = c(1, 2), mar = c(4, 5, 2, 1))
for (k in names(BANDS)) {
  yy <- sapply(PATHS, function(p) getc(p, k)$sig)
  matplot(ACID, yy, type = "b", pch = 16, lty = 1, xlab = "Atago acidity (%)", ylab = "2nd derivative (per index^2)",
          main = paste(k, "nm band"))
  legend("bottomleft", paste(PATHS, "mm"), col = 1:5, lty = 1, pch = 16, cex = 0.7)
}
dev.off()

# ---- 動手改改看（第 3 節）----
# 1. 刪掉 0.25 mm 的 4.87%（ACID 與 sig 各拿掉第 3 個元素，n = 3）：斜率與 R2 怎麼變？報告說它疑似整體縮小約 15%。
# 2. 把 x 改成標示值 c(0,3,5,10)，1 mm 的 R2 是否更高？（第 8 節會正式做）
# 3. 為什麼 n = 4 時 R2 接近 1 不代表模型好？

# =====================================================================
# 4. Beer-Lambert 正比檢查
# =====================================================================
hr("4. Beer-Lambert 正比檢查：靈敏度比 slope(p)/slope(0.25 mm) vs 理想 p/0.25")

ratio_of <- function(p, k) getc(p, k)$slope / getc(0.25, k)$slope
tab4 <- do.call(rbind, lapply(names(BANDS), function(k) data.frame(
  band_nm = k, path_mm = PATHS, ratio = sapply(PATHS, ratio_of, k = k), ideal = PATHS / 0.25,
  ratio_over_ideal = sapply(PATHS, ratio_of, k = k) / (PATHS / 0.25))))
print(tab4, row.names = FALSE, digits = 3)

png(file.path(FIG_DIR, "fig06_sensitivity_ratio.png"), width = 600, height = 450, res = 110)
par(mar = c(4, 4, 2, 1))
plot(PATHS, PATHS / 0.25, log = "x", type = "l", lty = 3, ylim = c(-2, 22), xlab = "path length (mm)", ylab = "sensitivity ratio")
lines(PATHS, sapply(PATHS, ratio_of, k = "2258"), type = "b", pch = 16, col = "blue")
lines(PATHS, sapply(PATHS, ratio_of, k = "1670"), type = "b", pch = 15, col = "red")
legend("topleft", c("ideal (p / 0.25)", "2258 nm", "1670 nm"), lty = c(3, 1, 1), pch = c(NA, 16, 15), col = c("black", "blue", "red"), cex = 0.8)
dev.off()

# ---- 動手改改看（第 4 節）----
# 1. 若以 0.5 mm 當基準，理想值要怎麼改？結論一樣嗎？
# 2. 為什麼 5 mm 的比值是負的？（回到第 0 節圖看 2258 nm 附近的光譜形狀）

# =====================================================================
# 5. 基礎④：LOD = 3 sigma / |斜率|
# =====================================================================
hr("5. 基礎④：LOD = 3σ/|斜率|")

# ---- 5.1 雜訊如何傳遞到二階導數：sigma_d2 = sigma_A * sqrt(sum(c2^2)) ----
set.seed(1)
noise <- rnorm(200000, 0, 1)
d2_noise <- as.numeric(stats::filter(noise, rev(c2), sides = 2))   # 對稱係數，rev 與否相同
cat(sprintf("模擬：二階導數雜訊 sd = %.4f\n", sd(d2_noise, na.rm = TRUE)))
cat(sprintf("公式：sigma_A * sqrt(sum(c2^2)) = %.4f\n", sqrt(sum(c2^2))))
cat(sprintf("對照：有限差分 [1,-2,1] 的放大倍率 = %.4f\n", sqrt(6)))

# ---- 5.2 照原程式的算法：局部雜訊（特徵帶索引 +/- 8 點）----
# 母體標準差（除以 N，等同 numpy 的 np.std；R 的 sd() 是除以 N-1，不要用）
pop_sd <- function(v) { v <- as.numeric(v); sqrt(mean((v - mean(v))^2)) }
c_smooth <- sg_coeffs(SMOOTH_WIN, 0)                 # 7 點二次「平滑」係數
SMOOTH <- lapply(spec, sg_apply, coef = c_smooth)
sum_c2 <- sqrt(sum(c2^2))

NOISE_LOC <- list(); LOD_LOC <- list()
for (p in PATHS) {
  resid <- spec[[pk(p)]] - SMOOTH[[pk(p)]]            # 原始 - 平滑
  for (k in names(BANDS)) {
    i <- IDX[[k]]
    loc <- max(i - 8, 1):min(i + 8, n_wl)
    r_loc <- resid[, loc]
    stopifnot(!anyNA(r_loc))
    sA <- pop_sd(r_loc)
    NOISE_LOC[[paste(p, k)]] <- sA
    LOD_LOC[[paste(p, k)]] <- 3 * sA * sum_c2 / abs(getc(p, k)$slope)
  }
}
tab5 <- do.call(rbind, lapply(names(BANDS), function(k) data.frame(
  path_mm = PATHS, band_nm = k,
  noise_A = sapply(PATHS, function(p) NOISE_LOC[[paste(p, k)]]),
  LOD_pct = sapply(PATHS, function(p) LOD_LOC[[paste(p, k)]]))))
print(tab5, row.names = FALSE, digits = 4)

# ---- 5.3 陷阱：吸光度越高光譜曲率越大，殘差把真實譜形也算成雜訊 ----
# 改在平坦且未飽和的 1550-1650 nm 估算（4 條光譜、7 點二次平滑殘差、母體標準差，共 16 點）
flat <- which(wl >= 1550 & wl <= 1650)
cat("\n平坦區 1550-1650 nm 共", length(flat), "個點\n")
NOISE_FLAT <- list(); LOD_FLAT <- list()
for (p in PATHS) {
  r_flat <- (spec[[pk(p)]] - SMOOTH[[pk(p)]])[, flat]
  stopifnot(!anyNA(r_flat))
  NOISE_FLAT[[pk(p)]] <- pop_sd(r_flat)
  for (k in names(BANDS)) LOD_FLAT[[paste(p, k)]] <- 3 * NOISE_FLAT[[pk(p)]] * sum_c2 / abs(getc(p, k)$slope)
}
print(data.frame(path_mm = PATHS,
                 noise_flat_AU = sapply(PATHS, function(p) NOISE_FLAT[[pk(p)]]),
                 noise_local_2258_AU = sapply(PATHS, function(p) NOISE_LOC[[paste(p, "2258")]]),
                 LOD_flat_2258_pct = sapply(PATHS, function(p) LOD_FLAT[[paste(p, "2258")]]),
                 LOD_local_2258_pct = sapply(PATHS, function(p) LOD_LOC[[paste(p, "2258")]])), row.names = FALSE, digits = 4)
cat("\n報告的平坦區雜訊：0.25 mm 0.00018、0.5 mm 0.00026、1 mm 0.00027 AU\n")
cat(sprintf("本程式算出來：      0.25 mm %.5f、0.5 mm %.5f、1 mm %.5f AU\n",
            NOISE_FLAT[["0.25"]], NOISE_FLAT[["0.5"]], NOISE_FLAT[["1"]]))
cat(sprintf("2258 nm 帶 LOD(0.25 mm) / LOD(1 mm) = %.2f  (報告：1 mm 約為 0.25 mm 的 1/2.6)\n",
            LOD_FLAT[["0.25 2258"]] / LOD_FLAT[["1 2258"]]))

png(file.path(FIG_DIR, "fig07_noise_lod.png"), width = 1100, height = 420, res = 110)
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
nl <- sapply(PATHS, function(p) NOISE_LOC[[paste(p, "2258")]]); nf <- sapply(PATHS, function(p) NOISE_FLAT[[pk(p)]])
plot(PATHS, nl, log = "xy", type = "b", pch = 16, ylim = range(c(nl, nf)), xlab = "path length (mm)", ylab = "noise (AU)")
lines(PATHS, nf, type = "b", pch = 15, col = "red")
legend("topleft", c("local (2258 nm +/- 8 pts)", "flat region 1550-1650 nm"), pch = c(16, 15), col = c("black", "red"), lty = 1, cex = 0.8)
ll <- sapply(PATHS, function(p) LOD_LOC[[paste(p, "2258")]]); lf <- sapply(PATHS, function(p) LOD_FLAT[[paste(p, "2258")]])
plot(PATHS, ll, log = "xy", type = "b", pch = 16, ylim = range(c(ll, lf)), xlab = "path length (mm)", ylab = "LOD (% acid)")
lines(PATHS, lf, type = "b", pch = 15, col = "red")
legend("topleft", c("LOD (local noise)", "LOD (flat-region noise)"), pch = c(16, 15), col = c("black", "red"), lty = 1, cex = 0.8)
dev.off()

# ---- 動手改改看（第 5 節）----
# 1. 把 SMOOTH_WIN 改成 9 或 5，平坦區雜訊變大還是變小？為什麼？
# 2. 把平坦區改成 1900-2050 nm（水的強吸收帶），雜訊與 LOD 會怎樣？
# 3. 用 replicate == 2 的 0.25 mm 水與 replicate == 1 的水相減，算 pop_sd(a - b)：
#    這是「填樣再現性」的下限，和平坦區雜訊（0.00018 AU）相比是大還是小？
#    提示：df[df$replicate == 2, "absorbance"] 與 spec[["0.25"]][1, ]
# 4. 局部雜訊範圍從 +/-8 點改成 +/-3 點，LOD 怎麼變？

# =====================================================================
# 6. 基礎⑤：留一交叉驗證 (LOOCV)
# =====================================================================
hr("6. 基礎⑤：留一交叉驗證與 RMSECV")

loocv <- function(x, y) {
  e <- numeric(length(x))
  for (j in seq_along(x)) {
    f <- fit_line(x[-j], y[-j])                    # 留下第 j 點，用其他點擬合
    x_hat <- (y[j] - f$intercept) / f$slope        # 預測被留下點的濃度
    e[j] <- x[j] - x_hat
  }
  list(rmsecv = sqrt(mean(e^2)), e = e)
}
lt2 <- loocv(xt, yt)
cat("玩具例子：每次的誤差 =", round(lt2$e, 4), " RMSECV =", round(lt2$rmsecv, 4), "\n")
for (p in PATHS) for (k in names(BANDS)) CAL[[paste(p, k)]]$rmsecv <- loocv(ACID, getc(p, k)$sig)$rmsecv
tab6 <- do.call(rbind, lapply(names(BANDS), function(k) data.frame(
  path_mm = PATHS, band_nm = k, R2 = sapply(PATHS, function(p) getc(p, k)$r2),
  RMSECV_pct = sapply(PATHS, function(p) getc(p, k)$rmsecv))))
print(tab6, row.names = FALSE, digits = 4)

# ---- 動手改改看（第 6 節）----
# 1. 只留 n = 3 點再做 LOOCV，還能做嗎？誤差怎麼變？
# 2. 把玩具例子 yt[4] 從 7.8 改成 9.0（離群點），RMSECV 與 R2 哪個反應比較大？
# 3. 為什麼 n = 4 的交叉驗證結果本身就很不穩定？

# =====================================================================
# 7. 決策規則
# =====================================================================
hr("7. 決策規則：在 Beer-Lambert 成立的範圍內取最長光徑")

band_amax <- function(p, k) {
  i <- IDX[[k]]; hh <- (SG_WIN - 1) / 2
  max(spec[[pk(p)]][, max(i - hh, 1):min(i + hh, n_wl)])
}
dec <- data.frame(path_mm = PATHS, ratio = NA_real_, ideal = PATHS / 0.25, ratio_over_ideal = NA_real_,
                  Amax_2258 = NA_real_, BL_ok = NA, unsaturated = NA, pass = NA)
for (r in seq_along(PATHS)) {
  p <- PATHS[r]
  dec$ratio[r] <- ratio_of(p, "2258")
  dec$ratio_over_ideal[r] <- dec$ratio[r] / dec$ideal[r]
  dec$Amax_2258[r] <- band_amax(p, "2258")
  dec$BL_ok[r] <- abs(dec$ratio_over_ideal[r] - 1) <= RATIO_TOL
  dec$unsaturated[r] <- dec$Amax_2258[r] < A_SAT
  dec$pass[r] <- dec$BL_ok[r] && dec$unsaturated[r]
}
print(dec, row.names = FALSE, digits = 3)
SELECTED <- max(dec$path_mm[dec$pass])
cat(sprintf("\n選擇的光徑 = %g mm\n", SELECTED))

usable_segments <- function(S, thr) {
  ok <- apply(S, 2, max) < thr
  segs <- character(0); i <- 1
  while (i <= n_wl) {
    if (ok[i]) {
      j <- i
      while (j + 1 <= n_wl && ok[j + 1]) j <- j + 1
      if (j > i) segs <- c(segs, sprintf("%.0f-%.0f", wl[i], wl[j]))   # 至少 2 點
      i <- j + 1
    } else i <- i + 1
  }
  segs
}
USABLE <- lapply(spec, usable_segments, thr = A_LIM)
cat(sprintf("\n各光徑可用波段（4 條光譜皆 A < %g）：\n", A_LIM))
for (p in PATHS) cat(sprintf("  %5g mm: %s\n", p, if (length(USABLE[[pk(p)]])) paste(USABLE[[pk(p)]], collapse = ", ") else "(無)"))
cat("（1 mm 多出的 1347-1413 是光譜最短波長端、被 1450 nm 水帶隔開的孤立小段，報告沒有採用。）\n")

# ---- 動手改改看（第 7 節）----
# 1. RATIO_TOL 改成 0.02 會選到哪個光徑？改成 0.6 呢？
# 2. A_LIM 改成 2.0，1 mm 的可用波段變成什麼？
# 3. 改用 1670 nm 帶判斷，結論會一樣嗎？

# =====================================================================
# 8. 延伸：改用其他參考值
# =====================================================================
hr("8. 延伸：標示值與手動滴定值作參考（1 mm）")

REFS <- list("Atago (0/3.32/4.87/7.45)" = c(0, 3.32, 4.87, 7.45),
             "label (0/3/5/10)" = c(0, 3, 5, 10),
             "titration (0/2.28/4.40/11.39)" = c(0, 2.28, 4.40, 11.39))
R2_EXT <- list()
for (nm in names(REFS)) for (k in names(BANDS)) R2_EXT[[paste(nm, k)]] <- fit_line(REFS[[nm]], getc(1, k)$sig)$r2
print(data.frame(reference = names(REFS),
                 R2_2258 = sapply(names(REFS), function(nm) R2_EXT[[paste(nm, "2258")]]),
                 R2_1670 = sapply(names(REFS), function(nm) R2_EXT[[paste(nm, "1670")]]), row.names = NULL), digits = 3)

png(file.path(FIG_DIR, "fig08_reference_compare.png"), width = 1100, height = 440, res = 110)
par(mfrow = c(1, 2), mar = c(4, 5, 2, 1))
for (k in names(BANDS)) {
  plot(NA, xlim = c(0, 11.5), ylim = range(getc(1, k)$sig), xlab = "reference acidity (%)", ylab = "2nd derivative (1 mm)",
       main = paste(k, "nm"))
  for (i in seq_along(REFS)) lines(REFS[[i]], getc(1, k)$sig, type = "b", pch = 16, col = i + 1)
  legend("bottomleft", c("Atago", "label", "titration"), col = 2:4, pch = 16, lty = 1, cex = 0.8)
}
dev.off()

# ---- 動手改改看（第 8 節）----
# 1. 把 0.25 mm 與 2 mm 的 R2 也用標示值算一次，短光徑的改善與 1 mm 一樣明顯嗎？
# 2. 滴定/Atago 的比值（0.76 -> 1.07 -> 1.66）不是常數，為什麼這排除了「NaOH 力價誤差」的解釋？
# 3. 要建立真正的校正，你會怎麼改實驗設計？

# =====================================================================
# 9. 自我驗證
# =====================================================================
hr("9. 自我驗證（比值與 R2 容許 +/-0.01，其他相對誤差 2%）")
# 如果改過上面的參數（A_LIM、SG_WIN 等），許多項目 FAIL 是正常的；改回預設、整支重跑即可。

checks <- data.frame(ok = logical(0), item = character(0), got = numeric(0), key = numeric(0), stringsAsFactors = FALSE)
chk <- function(name, got, key, kind, tol = NA) {
  ok <- switch(kind,
               abs = abs(got - key) <= 0.01 + 1e-12,
               rel = abs(got - key) <= 0.02 * abs(key),
               tol = abs(got - key) <= tol)
  checks[nrow(checks) + 1, ] <<- list(ok, name, got, key)
}

K2258 <- list(  # A<1.5, A<2.0, Amax, slope, ratio, R2, RMSECV, noise_A, LOD
  "0.25" = c(1.000, 1.000, 0.5877, -2.497e-4, 1.00, 0.9589, 1.076, 4.353e-4, 0.5961),
  "0.5"  = c(0.8662, 0.9648, 0.8827, -4.811e-4, 1.93, 0.8533, 3.344, 9.618e-4, 0.6835),
  "1"    = c(0.6197, 0.8099, 1.575, -9.852e-4, 3.95, 0.9342, 2.182, 2.218e-3, 0.7698),
  "2"    = c(0.3662, 0.4437, 2.225, -8.427e-4, 3.38, 0.9498, 1.223, 2.068e-3, 0.8391),
  "5"    = c(0.01408, 0.169, 2.245, 8.893e-5, -0.36, 0.3229, 17.16, 1.396e-3, 5.367))
K1670 <- list(  # Amax, slope, R2, RMSECV, noise_A, LOD  （取自 pathlength_metrics.csv）
  "0.25" = c(0.4232, -2.859e-05, 0.855, 3.712, 2.997e-4, 3.584),
  "0.5"  = c(0.4706, -7.549e-05, 0.7873, 4.173, 4.357e-4, 1.973),
  "1"    = c(0.7422, -1.535e-04, 0.9448, 1.945, 5.066e-4, 1.128),
  "2"    = c(1.073, -3.010e-04, 0.8508, 3.493, 8.118e-4, 0.9219),
  "5"    = c(1.909, -4.342e-04, 0.999, 0.1577, 1.917e-3, 1.509))
for (p in PATHS) {
  K <- K2258[[pk(p)]]; t <- sprintf("2258 %g mm ", p); S <- spec[[pk(p)]]
  chk(paste0(t, "frac A<1.5"), frac_below(S, 1.5), K[1], "abs")
  chk(paste0(t, "frac A<2.0"), frac_below(S, 2.0), K[2], "abs")
  chk(paste0(t, "Amax(band)"), band_amax(p, "2258"), K[3], "rel")
  chk(paste0(t, "slope"), getc(p, "2258")$slope, K[4], "rel")
  chk(paste0(t, "ratio"), ratio_of(p, "2258"), K[5], "abs")
  chk(paste0(t, "R2"), getc(p, "2258")$r2, K[6], "abs")
  chk(paste0(t, "RMSECV%"), getc(p, "2258")$rmsecv, K[7], "rel")
  chk(paste0(t, "noise_A"), NOISE_LOC[[paste(p, "2258")]], K[8], "rel")
  chk(paste0(t, "LOD%"), LOD_LOC[[paste(p, "2258")]], K[9], "rel")
  K <- K1670[[pk(p)]]; t <- sprintf("1670 %g mm ", p)
  chk(paste0(t, "Amax(band)"), band_amax(p, "1670"), K[1], "rel")
  chk(paste0(t, "slope"), getc(p, "1670")$slope, K[2], "rel")
  chk(paste0(t, "R2"), getc(p, "1670")$r2, K[3], "abs")
  chk(paste0(t, "RMSECV%"), getc(p, "1670")$rmsecv, K[4], "rel")
  chk(paste0(t, "noise_A"), NOISE_LOC[[paste(p, "1670")]], K[5], "rel")
  chk(paste0(t, "LOD%"), LOD_LOC[[paste(p, "1670")]], K[6], "rel")
}
chk("hand-made SG vs lm window fit (max diff < 1e-9)", maxdiff, 0, "tol", 1e-9)
for (pp in list(c(0.25, 0.00018), c(0.5, 0.00026), c(1, 0.00027)))
  chk(sprintf("flat-region noise %g mm", pp[1]), NOISE_FLAT[[pk(pp[1])]], pp[2], "tol", 6e-6)
chk("LOD ratio 0.25mm/1mm (2258, flat noise)", LOD_FLAT[["0.25 2258"]] / LOD_FLAT[["1 2258"]], 2.6, "tol", 0.1)
# 1 mm 的程式輸出多一段 1347-1413（孤立小段），報告沒列入；所以用「答案鍵的區段都出現」來比對
chk("usable bands 1 mm (key segs present)", as.numeric(all(c("1480-1871", "2116-2284") %in% USABLE[["1"]])), 1, "tol", 0)
chk("usable bands 0.5 mm (key segs present)", as.numeric(all(c("1347-1889", "2011-2435") %in% USABLE[["0.5"]])), 1, "tol", 0)
chk("selected path = 1 mm", SELECTED, 1, "tol", 0)
chk("ext 2258 label R2", R2_EXT[["label (0/3/5/10) 2258"]], 0.996, "abs")
chk("ext 1670 label R2", R2_EXT[["label (0/3/5/10) 1670"]], 0.998, "abs")
chk("ext 2258 titration R2", R2_EXT[["titration (0/2.28/4.40/11.39) 2258"]], 0.990, "abs")
chk("ext 1670 titration R2", R2_EXT[["titration (0/2.28/4.40/11.39) 1670"]], 0.987, "abs")

cat(sprintf("%-6s | %-50s | %14s | %14s\n", "result", "item", "yours", "key"))
cat(strrep("-", 92), "\n", sep = "")
for (r in seq_len(nrow(checks)))
  cat(sprintf("%-6s | %-50s | %14s | %14s\n", if (checks$ok[r]) "PASS" else "FAIL", checks$item[r],
              formatC(checks$got[r], digits = 6, format = "g"), formatC(checks$key[r], digits = 6, format = "g")))
cat(strrep("-", 92), "\n", sep = "")
cat(sprintf("TOTAL: %d / %d PASS\n", sum(checks$ok), nrow(checks)))
cat("圖已存到", normalizePath(FIG_DIR), "\n")
