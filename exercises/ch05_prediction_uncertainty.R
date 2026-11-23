# =====================================================================
# 第 5 章 預測未知樣品與不確定度：R 程式練習（只用 base R）
# 用單一波長 2258 nm 的二階導數訊號做「反向校正」y = a + b*c，預測未知樣品的酸度 c0，並算出誠實的 95 % 區間
#
# 執行方式：Rscript ch05_prediction_uncertainty.R （程式會自動切到本檔所在資料夾）
# 資料：../data/spectra_long.csv      圖：figures_R/ch05_*.png（圖上文字一律用英文）
# 每節結尾有「動手改改看」，每節的最後有該節的自我驗證 PASS/FAIL；第 7 節是報告格式與學生練習。
#
# 方法：訊號 y = 1 mm、重複 1、SG 9 點二次二階導數（點距 1）在 2258 nm 的值 x 10^3。
#   校正：y = a + b*c；s = sqrt(殘差平方和/(n-2))；t = qt(0.975, n-2)；Sxx = sum((c - mean(c))^2)
#   預測：c0 = (y0 - a)/b
#   簡易半寬：t*s/|b| * sqrt(1/m + 1/n + (y0 - ybar)^2/(b^2*Sxx))；m = 未知樣品重複量測次數
#   Fieller 精確區間：解 |y0 - a - b*c| <= t*s*sqrt(1/m + 1/n + (c - cbar)^2/Sxx) 的二次不等式
#   g = t^2*s^2/(b^2*Sxx)；g >= 1 表示斜率不顯著不為 0，區間無界。
# =====================================================================

options(width = 130, digits = 7)
.args <- commandArgs(trailingOnly = FALSE)
.fa <- sub("^--file=", "", .args[grep("^--file=", .args)])
if (length(.fa) == 1 && nzchar(.fa)) setwd(dirname(normalizePath(.fa)))
DATA_PATH <- "../data/spectra_long.csv"
FIG_DIR <- "figures_R"
dir.create(FIG_DIR, showWarnings = FALSE)
hr <- function(title) cat("\n", strrep("=", 70), "\n", title, "\n", strrep("=", 70), "\n", sep = "")

# ---- 事先固定的設定 ----
SG_WIN <- 9            # 二階導數：Savitzky-Golay 9 點、二次多項式
WL0    <- 2258         # 單一波長（nm）
Y0     <- 4.20         # 假想未知樣品的訊號（x 10^-3）

# ---- 自我驗證的小工具：每節結尾呼叫 show_checks() ----
checks <- data.frame(ok = logical(0), item = character(0), got = numeric(0), key = numeric(0), stringsAsFactors = FALSE)
.shown <- 0
chk <- function(name, got, key, tol) {
  checks[nrow(checks) + 1, ] <<- list(abs(got - key) <= tol, name, got, key)
}
show_checks <- function(title) {
  rows <- seq_len(nrow(checks)); rows <- rows[rows > .shown]
  cat(sprintf("\n[自我驗證] %s\n", title))
  cat(sprintf("%-6s | %-52s | %11s | %11s\n", "result", "item", "yours", "key")); cat(strrep("-", 88), "\n", sep = "")
  for (r in rows)
    cat(sprintf("%-6s | %-52s | %11s | %11s\n", if (checks$ok[r]) "PASS" else "FAIL", checks$item[r],
                formatC(checks$got[r], digits = 5, format = "g"), formatC(checks$key[r], digits = 5, format = "g")))
  cat(sprintf("本節：%d / %d PASS\n", sum(checks$ok[rows]), length(rows)))
  .shown <<- nrow(checks)
}

# =====================================================================
# 0. 載入資料、2258 nm 二階導數訊號
# =====================================================================
hr("0. 載入資料：1 mm、重複 1 的 4 個樣品，2258 nm 二階導數訊號")
df <- read.csv(DATA_PATH)
d1 <- df[df$replicate == 1 & abs(df$path_mm - 1) < 1e-9, ]
atago <- sort(unique(d1$acid_pct))                       # 0, 3.32, 4.87, 7.45
wl <- sort(unique(d1$wavelength_nm)); n_wl <- length(wl)
X1 <- t(sapply(atago, function(a) {
  h <- d1[abs(d1$acid_pct - a) < 1e-9, ]; h <- h[order(h$wavelength_nm), ]
  stopifnot(nrow(h) == n_wl); h$absorbance
}))                                                      # 4 x 142

# 二階導數：SG 9 點、二次、點距為 1（與第 3、4 章相同）。2258 nm 離邊緣很遠，邊緣怎麼處理不影響結果。
sg_coeffs2 <- function(window) {
  hh <- (window - 1) / 2; xx <- -hh:hh
  M <- cbind(1, xx, xx^2)
  as.numeric(2 * (solve(t(M) %*% M) %*% t(M))[3, ])
}
sg_d2 <- function(S, window = SG_WIN) {
  cf <- sg_coeffs2(window); half <- (window - 1) / 2
  out <- matrix(NA_real_, nrow(S), ncol(S))
  for (i in (half + 1):(ncol(S) - half)) out[, i] <- S[, (i - half):(i + half), drop = FALSE] %*% cf
  out[, 1:half] <- out[, half + 1]
  out[, (ncol(S) - half + 1):ncol(S)] <- out[, ncol(S) - half]
  out
}
D1 <- sg_d2(X1)                                          # 4 x 142 二階導數
j <- which.min(abs(wl - WL0))                            # 最接近 2258 nm 的波長
sig <- D1[, j] * 1e3                                     # 訊號 y（x 10^-3）
titr <- c(0, 8.8, 17, 44) * c(0, rep(1.07875 * 0.24, 3)) # 滴定參考值（%）：V x 1.07875 x 0.24；純水 0
cat(sprintf("波長 = %.1f nm\n", wl[j]))
cat("Atago 酸度（%）：", atago, "\n")
cat("滴定酸度（%）  ：", round(titr, 3), "\n")
cat("2258 nm 訊號 y ：", round(sig, 4), "\n")

# ---- 工具函式 ----
# 直線校正 y = a + b*c。回傳 list：n、a、b、s（殘差標準差）、sxx、xbar、ybar、t、tb（斜率的 t 值）、g
calib <- function(xs, ys) {
  n <- length(xs); xbar <- mean(xs); ybar <- mean(ys); sxx <- sum((xs - xbar)^2)
  b <- sum((xs - xbar) * (ys - ybar)) / sxx; a <- ybar - b * xbar
  s <- sqrt(sum((ys - (a + b * xs))^2) / (n - 2))
  t <- qt(0.975, n - 2)
  list(n = n, a = a, b = b, s = s, sxx = sxx, xbar = xbar, ybar = ybar, t = t,
       tb = abs(b) / (s / sqrt(sxx)), g = t^2 * s^2 / (b^2 * sxx))
}
# 預測 c0 = (y0 - a)/b、簡易半寬、Fieller 精確區間（g >= 1 時 lo = hi = NA，表示無界）
predict_c <- function(c, y0, m = 1) {
  x0 <- (y0 - c$a) / c$b
  simple <- c$t * c$s / abs(c$b) * sqrt(1 / m + 1 / c$n + (y0 - c$ybar)^2 / (c$b^2 * c$sxx))
  # Fieller：解 (y0 - a - b*x)^2 = t^2 s^2 (1/m + 1/n + (x - xbar)^2/Sxx) 的二次方程 A x^2 + B x + C = 0
  k <- c$t^2 * c$s^2
  A <- c$b^2 - k / c$sxx
  B <- -2 * c$b * (y0 - c$a) + 2 * k * c$xbar / c$sxx
  C <- (y0 - c$a)^2 - k * (1 / m + 1 / c$n) - k * c$xbar^2 / c$sxx
  if (A <= 0) { lo <- NA_real_; hi <- NA_real_ }          # A <= 0 等價於 g >= 1：斜率不顯著，區間無界
  else { r <- sqrt(B * B - 4 * A * C); rt <- sort(c((-B - r) / (2 * A), (-B + r) / (2 * A))); lo <- rt[1]; hi <- rt[2] }
  list(x0 = x0, simple = simple, lo = lo, hi = hi)
}
# 校正直線、點、95 % 預測帶（半寬 t s sqrt(1/m + 1/n + (c-cbar)^2/Sxx)）；可加 y0 水平線與藏起來的點
plot_calib <- function(c, xs, ys, xlim, y0 = NULL, m = 1, hide = NULL, xlab = "reference acid (%)") {
  xx <- seq(xlim[1], xlim[2], length.out = 200); fit <- c$a + c$b * xx
  hw <- c$t * c$s * sqrt(1 / m + 1 / c$n + (xx - c$xbar)^2 / c$sxx)
  yr <- range(c(fit - hw, fit + hw, ys, y0, hide[2]))
  yr <- c(max(min(yr), -25), min(max(yr), 35))
  plot(xx, fit, type = "n", xlim = xlim, ylim = yr, xlab = xlab, ylab = "signal y (x 1e-3)")
  polygon(c(xx, rev(xx)), c(fit - hw, rev(fit + hw)), col = adjustcolor("steelblue", 0.15), border = NA)
  lines(xx, fit - hw, col = "steelblue", lty = 3); lines(xx, fit + hw, col = "steelblue", lty = 3)
  lines(xx, fit, col = "steelblue", lwd = 2)
  points(xs, ys, pch = 19, cex = 1.3)
  if (!is.null(hide)) points(hide[1], hide[2], pch = 15, col = "red", cex = 1.5)
  if (!is.null(y0)) abline(h = y0, col = "red", lty = 2)
  leg <- c("calibration line", sprintf("95%% prediction band (m=%d)", m), "calibration points")
  lt <- c(1, 3, NA); pc <- c(NA, NA, 19); cl <- c("steelblue", "steelblue", "black")
  if (!is.null(hide)) { leg <- c(leg, "hidden bottle"); lt <- c(lt, NA); pc <- c(pc, 15); cl <- c(cl, "red") }
  if (!is.null(y0)) { leg <- c(leg, sprintf("y0 = %.2f", y0)); lt <- c(lt, 2); pc <- c(pc, NA); cl <- c(cl, "red") }
  legend("topright", leg, lty = lt, pch = pc, col = cl, cex = 0.75)
}

chk("signal y at 0 %", sig[1], 8.5087, 0.001)
chk("signal y at 3.32 %", sig[2], 6.3073, 0.001)
chk("signal y at 4.87 %", sig[3], 5.1725, 0.001)
chk("signal y at 7.45 %", sig[4], 0.9509, 0.001)
show_checks("第 0 節")

# ---- 動手改改看（第 0 節）----
# 1. 訊號隨酸度增加而變小（斜率為負）。為什麼二階導數的吸收峰是「負的尖峰」？這會不會影響預測？
# 2. 把 WL0 改成 2250，訊號變多少？（改完記得改回 2258 再看自我驗證）

# =====================================================================
# 1. 用 Atago 4 點建立校正
# =====================================================================
hr("1. Atago 4 點校正")
c4 <- calib(atago, sig)
cat(sprintf("斜率 b = %.4f（訊號/%%）；截距 a = %.4f；殘差標準差 s = %.4f\n", c4$b, c4$a, c4$s))
cat(sprintf("t(0.975, %d) = %.4f；斜率的 t 值 = %.4f；g = %.4f\n", c4$n - 2, c4$t, c4$tb, c4$g))
cat(if (c4$g < 1) "g < 1：斜率顯著不為 0，Fieller 區間有界。\n" else "g >= 1：區間無界。\n")
cat("n = 4，自由度只有 2，t(0.975, 2) = 4.30 很大：點少時，t 值本身就是不確定度的重要來源。\n")

png(file.path(FIG_DIR, "ch05_calib_atago.png"), width = 800, height = 500, res = 110)
par(mar = c(4, 4, 2, 1))
plot_calib(c4, atago, sig, c(-1, 9))
title("Calibration with Atago reference (n = 4)")
dev.off()

chk("Atago 4-point: slope b", c4$b, -0.9852, 0.002)
chk("Atago 4-point: s", c4$s, 0.9975, 0.002)
chk("Atago 4-point: t(0.975, 2)", c4$t, 4.3027, 0.002)
chk("Atago 4-point: slope t-statistic", c4$tb, 5.3267, 0.002)
chk("Atago 4-point: g", c4$g, 0.6525, 0.002)
show_checks("第 1 節")

# ---- 動手改改看（第 1 節）----
# 1. 斜率的 t 值 = 5.33，聽起來「很顯著」，為什麼 g 還高達 0.65？（提示：g = t^2/tb^2，即 t(0.975, 2) = 4.30 與斜率 t 值的比）
# 2. 若有 8 個校正點、s 和 b 不變，t(0.975, 6) 是多少？g 會變多少？

# =====================================================================
# 2. 藏起 4.87 % 那瓶：只用 3 點校正
# =====================================================================
hr("2. 藏起 4.87 % 的樣品，3 點校正")
tr <- c(1, 2, 4)                                         # 校正用 0、3.32、7.45 %；藏起第 3 個（4.87 %）
c3 <- calib(atago[tr], sig[tr]); p3 <- predict_c(c3, sig[3])
cat(sprintf("3 點校正：b = %.4f；s = %.4f；t(0.975, 1) = %.4f；g = %.4f\n", c3$b, c3$s, c3$t, c3$g))
cat(sprintf("藏起來的那瓶：訊號 y = %.4f -> 預測酸度 = %.4f %%\n", sig[3], p3$x0))
cat(sprintf("對照：Atago 讀值 = %.2f %%；滴定 = %.3f %%\n", atago[3], titr[3]))
if (is.na(p3$lo)) {
  cat(sprintf("g = %.3f >= 1：斜率與 0 沒有顯著差異，Fieller 95 %% 區間無界（任何濃度都相容）。\n", c3$g))
} else cat(sprintf("Fieller 95 %% 區間 [%.2f, %.2f] %%\n", p3$lo, p3$hi))
cat("點預測 3.67 % 與 Atago 4.87 %、滴定 4.40 % 差約 1 %，而且我們「沒有辦法」宣稱它有多準。\n")

png(file.path(FIG_DIR, "ch05_hide_bottle.png"), width = 800, height = 500, res = 110)
par(mar = c(4, 4, 2, 1))
plot_calib(c3, atago[tr], sig[tr], c(-1, 9), hide = c(atago[3], sig[3]))
title("3-point calibration, 4.87% bottle hidden")
dev.off()

chk("3-point: slope b", c3$b, -1.0258, 0.002)
chk("3-point: s", c3$s, 0.9507, 0.002)
chk("3-point: t(0.975, 1)", c3$t, 12.7062, 0.002)
chk("3-point: g", c3$g, 4.9775, 0.002)
chk("3-point: interval unbounded (1 = yes)", as.numeric(is.na(p3$lo)), 1, 0)
chk("3-point: predicted acid of hidden bottle", p3$x0, 3.6711, 0.002)
show_checks("第 2 節")

# ---- 動手改改看（第 2 節）----
# 1. g = 4.98 遠大於 1，這在圖上看起來是什麼樣子？（預測帶與水平線 y0 的交會）
# 2. 點預測 3.67 % 與 Atago 4.87 %、滴定 4.40 % 相差 1 % 以上，卻不能說「錯了」。為什麼？（提示：區間無界）
# 3. 改藏第 2 個樣品（tr <- c(1, 3, 4)），預測值與 g 會怎麼變？

# =====================================================================
# 3. 假想未知樣品 y0 = 4.20（Atago 4 點）
# =====================================================================
hr("3. 假想未知樣品 y0 = 4.20：簡易區間 vs Fieller 區間")
p4 <- predict_c(c4, Y0); p4m3 <- predict_c(c4, Y0, 3)
cat(sprintf("預測 x0 = %.4f %%\n", p4$x0))
cat(sprintf("簡易區間：%.4f +/- %.4f -> [%.4f, %.4f]\n", p4$x0, p4$simple, p4$x0 - p4$simple, p4$x0 + p4$simple))
cat(sprintf("Fieller 區間 (m = 1)：[%.4f, %.4f]\n", p4$lo, p4$hi))
cat(sprintf("Fieller 區間 (m = 3)：[%.4f, %.4f]\n", p4m3$lo, p4m3$hi))
cat(sprintf("Fieller 區間不對稱：x0 到下限 %.2f，到上限 %.2f。\n", p4$x0 - p4$lo, p4$hi - p4$x0))
cat("下限是負的酸度（-1.68 %）沒有物理意義：區間太寬，代表資訊不足，不是「可能是負的」。\n")

png(file.path(FIG_DIR, "ch05_fieller_band.png"), width = 900, height = 520, res = 110)
par(mar = c(4, 4, 2, 1))
plot_calib(c4, atago, sig, c(-3, 17), y0 = Y0)
yl <- par("usr")[3] + 0.03 * diff(par("usr")[3:4])
segments(p4$lo, yl, p4$hi, yl, col = "red", lwd = 6)
abline(v = c(p4$lo, p4$hi), col = "red", lty = 3)
points(p4$x0, Y0, pch = 8, col = "red", cex = 1.6)
title("Prediction of y0 = 4.20 with Fieller interval (n = 4)")
legend("bottomleft", "Fieller 95% interval (m=1)", lty = 1, lwd = 4, col = "red", cex = 0.75)
dev.off()

chk("y0 = 4.20: predicted x0", p4$x0, 4.9604, 0.002)
chk("y0 = 4.20: simple half-width", p4$simple, 4.9442, 0.002)
chk("y0 = 4.20: Fieller lower", p4$lo, -1.6831, 0.002)
chk("y0 = 4.20: Fieller upper", p4$hi, 15.548, 0.002)
chk("y0 = 4.20, m = 3: Fieller lower", p4m3$lo, 0.7828, 0.002)
chk("y0 = 4.20, m = 3: Fieller upper", p4m3$hi, 13.0822, 0.002)
show_checks("第 3 節")

# ---- 動手改改看（第 3 節）----
# 1. 簡易區間 +/-4.94 對稱，Fieller [-1.68, 15.55] 不對稱。哪一個比較誠實？g = 0.65 愈接近 1 差異愈大，為什麼？
# 2. 量 3 次（m = 3）只把下限從 -1.68 拉到 0.78，不是縮成 1/sqrt(3)。為什麼？（提示：1/m 只是根號裡三項之一，另兩項來自校正本身）
# 3. 把 Y0 改成 8.0（接近 0 % 的訊號）或 1.0（超出校正範圍），區間怎麼變？

# =====================================================================
# 4. 改用滴定參考值校正（4 點）
# =====================================================================
hr("4. 滴定參考值 4 點校正")
c4t <- calib(titr, sig); p4t <- predict_c(c4t, Y0)
cat(sprintf("滴定 4 點：b = %.4f；s = %.4f；g = %.4f\n", c4t$b, c4t$s, c4t$g))
cat(sprintf("y0 = 4.20：x0 = %.4f %%；Fieller 區間 [%.4f, %.4f]；簡易 +/-%.4f\n", p4t$x0, p4t$lo, p4t$hi, p4t$simple))
cat("殘差 s 從 0.9975 降到 0.3927，g 從 0.65 降到 0.095，區間縮窄到 [3.16, 9.45]（寬度約 6.3 %；Atago 是 17.2 %）。\n")
cat("代價：滴定的最高點 11.39 % 拉長了範圍（Sxx 變大），也代表參考值的準確度直接決定預測的上限。\n")

png(file.path(FIG_DIR, "ch05_calib_titration.png"), width = 900, height = 520, res = 110)
par(mar = c(4, 4, 2, 1))
plot_calib(c4t, titr, sig, c(-1, 13), y0 = Y0, xlab = "titration acid (%)")
yl <- par("usr")[3] + 0.03 * diff(par("usr")[3:4])
segments(p4t$lo, yl, p4t$hi, yl, col = "red", lwd = 6)
points(p4t$x0, Y0, pch = 8, col = "red", cex = 1.6)
title("Calibration against titration reference (n = 4)")
legend("bottomleft", "Fieller 95% interval (m=1)", lty = 1, lwd = 4, col = "red", cex = 0.75)
dev.off()

chk("titration 4-point: slope b", c4t$b, -0.6415, 0.002)
chk("titration 4-point: s", c4t$s, 0.3927, 0.002)
chk("titration 4-point: g", c4t$g, 0.0954, 0.002)
chk("titration, y0 = 4.20: predicted x0", p4t$x0, 6.1309, 0.002)
chk("titration, y0 = 4.20: Fieller lower", p4t$lo, 3.1568, 0.002)
chk("titration, y0 = 4.20: Fieller upper", p4t$hi, 9.4453, 0.002)
show_checks("第 4 節")

# ---- 動手改改看（第 4 節）----
# 1. 同一個訊號 y0 = 4.20，用 Atago 校正預測 4.96 %，用滴定校正預測 6.13 %。這兩個數字差在哪裡？報告時必須說明用的是哪一種參考值。
# 2. s 變小的原因可能是什麼？（提示：Atago 的讀值解析度與滴定的不確定度）

# =====================================================================
# 5. 由殘差換算的偵測／定量下限
# =====================================================================
hr("5. 殘差換算的 LOD = 3s/|b|、LOQ = 10s/|b|")
lod_a <- 3 * c4$s / abs(c4$b);  loq_a <- 10 * c4$s / abs(c4$b)
lod_t <- 3 * c4t$s / abs(c4t$b); loq_t <- 10 * c4t$s / abs(c4t$b)
cat(sprintf("Atago 參考：LOD = %.3f %%；LOQ = %.3f %%\n", lod_a, loq_a))
cat(sprintf("滴定參考  ：LOD = %.3f %%；LOQ = %.3f %%\n", lod_t, loq_t))
cat("第 1 章的 LOD 只含儀器雜訊（重複掃描同一裝樣）；這裡的 s 還含裝樣、參考值與校正的誤差，才是「實際使用時」能看到的下限。\n")

chk("Atago residual LOD (3s/|b|)", lod_a, 3.038, 0.005)
chk("Atago residual LOQ (10s/|b|)", loq_a, 10.125, 0.005)
chk("titration residual LOD (3s/|b|)", lod_t, 1.836, 0.005)
chk("titration residual LOQ (10s/|b|)", loq_t, 6.120, 0.005)
show_checks("第 5 節")

# ---- 動手改改看（第 5 節）----
# 1. Atago 的 LOQ = 10.1 % 比校正範圍的上限 7.45 % 還大。這代表什麼？（我們的 4 個校正點，沒有一個真的「定量」得夠準）
# 2. 這個 LOD 與第 1 章的儀器雜訊 LOD 差在哪些誤差來源？

# =====================================================================
# 6. 設計示意：第 2 章的設計（0-7 %、8 濃度 x 3 次）
# =====================================================================
hr("6. 設計示意：8 濃度 x 3 次（n = 24），假設 s、b 不變")
design <- function(c, m) {
  xs <- rep(0:7, each = 3); n <- length(xs); sxx <- sum((xs - mean(xs))^2); t <- qt(0.975, n - 2)
  cc <- c; cc$n <- n; cc$sxx <- sxx; cc$xbar <- mean(xs); cc$ybar <- c$a + c$b * mean(xs); cc$t <- t
  cc$g <- t^2 * c$s^2 / (c$b^2 * sxx)
  r <- predict_c(cc, Y0, m); r$g <- cc$g; r
}
d_a1 <- design(c4, 1); d_a3 <- design(c4, 3); d_t3 <- design(c4t, 3)
cat(sprintf("Atago，設計 8 濃度 x 3 次，m = 1：g = %.4f；x0 = %.3f；區間 [%.4f, %.4f]\n", d_a1$g, d_a1$x0, d_a1$lo, d_a1$hi))
cat(sprintf("Atago，設計 8 濃度 x 3 次，m = 3：區間 [%.4f, %.4f]\n", d_a3$lo, d_a3$hi))
cat(sprintf("滴定，設計 8 濃度 x 3 次，m = 3：g = %.4f；x0 = %.3f；區間 [%.4f, %.4f]\n", d_t3$g, d_t3$x0, d_t3$lo, d_t3$hi))

rows <- list("Atago n=4, m=1" = p4, "Atago n=4, m=3" = p4m3, "Atago design, m=1" = d_a1, "Atago design, m=3" = d_a3,
             "titration n=4, m=1" = p4t, "titration design, m=3" = d_t3)
png(file.path(FIG_DIR, "ch05_intervals.png"), width = 900, height = 450, res = 110)
par(mar = c(4, 11, 2, 1))
k <- length(rows)
plot(NA, xlim = range(sapply(rows, function(p) c(p$lo, p$hi))), ylim = c(k + 0.5, 0.5), yaxt = "n", ylab = "",
     xlab = "predicted acid (%) with 95% Fieller interval", main = "Interval width for y0 = 4.20")
axis(2, at = seq_len(k), labels = names(rows), las = 1, cex.axis = 0.8)
for (i in seq_len(k)) {
  p <- rows[[i]]; col <- if (grepl("Atago", names(rows)[i])) "steelblue" else "darkorange"
  segments(p$lo, i, p$hi, i, lwd = 6, col = col, lend = 1); points(p$x0, i, pch = "|", cex = 1.5)
}
dev.off()

chk("design (24 points), Atago: g", d_a1$g, 0.035, 0.002)
chk("design, Atago m = 1: Fieller lower", d_a1$lo, 2.8134, 0.002)
chk("design, Atago m = 1: Fieller upper", d_a1$hi, 7.2134, 0.002)
chk("design, Atago m = 3: Fieller lower", d_a3$lo, 3.6741, 0.002)
chk("design, Atago m = 3: Fieller upper", d_a3$hi, 6.3526, 0.002)
chk("design, titration: g", d_t3$g, 0.0128, 0.002)
chk("design, titration m = 3: Fieller lower", d_t3$lo, 5.3266, 0.002)
chk("design, titration m = 3: Fieller upper", d_t3$hi, 7.0033, 0.002)
show_checks("第 6 節")

# ---- 動手改改看（第 6 節）----
# 1. 區間從 [-1.68, 15.55]（n = 4）縮到 [2.81, 7.21]（24 點）。縮窄的最主要原因是 t 值、n，還是 Sxx？（分別改 design 裡的 t、n、sxx 試試）
# 2. 設計若只有 8 個濃度、每個 1 次（n = 8），區間是多少？
# 3. 為什麼說這個示意「偏樂觀」？（提示：s 與 b 是從 4 點估出來的；真實的 24 點會發現新的誤差來源）

# =====================================================================
# 7. 怎麼報告結果
# =====================================================================
hr("7. 報告格式與學生練習")
fmt <- function(v) sub("-", "−", sprintf("%.2f", v), fixed = TRUE)   # 負號用 Unicode 減號
report <- function(x0, p, n, ref, wl_nm = 2258) {
  if (is.na(p$lo)) return(sprintf("酸度 %s %%（95 %% 區間無界，n = %d，%s 參考值，單一波長 %d nm）", fmt(x0), n, ref, wl_nm))
  sprintf("酸度 %s %%（95 %% 區間 %s 至 %s %%，n = %d，%s 參考值，單一波長 %d nm）", fmt(x0), fmt(p$lo), fmt(p$hi), n, ref, wl_nm)
}
cat(report(p4$x0, p4, c4$n, "Atago"), "\n")
cat(report(p3$x0, p3, c3$n, "Atago"), "\n")
# 學生練習：把上面的報告改成「滴定參考值、4 點」的版本，並再寫一行「設計 8 濃度 x 3 次、量 3 次」的版本。
# cat(report(p4t$x0, p4t, c4t$n, "滴定"), "\n")

# ---- 動手改改看（第 7 節）----
# 1. 取消註解最後一行，改成滴定參考值的報告；比較 Atago 與滴定兩行，哪個比較能讓讀者信任？
# 2. 加上第 6 節的設計示意，寫出一行「設計 8 濃度 x 3 次、量 3 次」的報告（注意：這是示意，不是實測）。
# 3. 一份報告若只寫「酸度 4.96 %」，讀者無法判斷什麼？

# =====================================================================
# 總結
# =====================================================================
hr("總結")
cat(sprintf("全部自我驗證：%d / %d PASS\n", sum(checks$ok), nrow(checks)))
if (!all(checks$ok)) { cat("FAIL 項目：\n"); print(checks[!checks$ok, ], row.names = FALSE) }
cat("1. 反向校正 y = a + b*c，預測 c0 = (y0 - a)/b；4 點 Atago：b = -0.985、s = 1.00、g = 0.65。\n")
cat("2. 藏起一瓶、只用 3 點：g = 4.98 >= 1，區間無界；點預測 3.67 % 不能宣稱準或不準。\n")
cat("3. y0 = 4.20：簡易 +/-4.94 對稱；Fieller [-1.68, 15.55] 不對稱且誠實；量 3 次只縮到 [0.78, 13.08]。\n")
cat("4. 滴定參考值（s = 0.39）區間縮到 [3.16, 9.45]；參考值的品質決定預測的品質。\n")
cat("5. 殘差換算的 LOD / LOQ（Atago 3.04 / 10.1 %；滴定 1.84 / 6.12 %）與第 1 章的儀器雜訊 LOD 不同。\n")
cat("6. 區間縮窄靠「更多、範圍更廣的校正點」（8 濃度 x 3 次：[2.81, 7.21]），而不是只靠重複量測。\n")
cat("圖已存到", normalizePath(FIG_DIR), "\n")
