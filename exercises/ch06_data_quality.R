# =====================================================================
# 第 6 章 資料品質與離群值：R 程式練習（只用 base R）
# 看到「可疑」的光譜時怎麼辦？先定規則、同類相比、找原因、重測；不要因為 R² 好看就刪資料。
#
# 執行方式：Rscript ch06_data_quality.R （程式會自動切到本檔所在資料夾）
# 資料：../data/spectra_long.csv      圖：figures_R/ch06_*.png（圖上文字一律用英文）
# 每節結尾有「動手改改看」，每節的最後有該節的自我驗證 PASS/FAIL；第 7 節是決策流程總結。
#
# 事先固定的規則（看結果之前就定好，事後不准改）：
#   水帶高度 = 1440-1460 nm 平均吸光度 - 1580-1650 nm（平坦區）平均吸光度
#   同一光徑 4 瓶的水帶高度，各自除以該光徑 4 瓶的中位數 = 比值；|比值 - 1| > 0.10 就標「可疑」
#   （醋裡的水最多少約 7 %，同光徑的水帶高度應該相差不到 10 %）
#   該光徑 1440-1460 nm 的最大吸光度 >= 1.5 時水帶已飽和，這個檢查不能用
#   二階導數：Savitzky-Golay 9 點、二次多項式、點距 1（與第 3、4、5 章相同），取 2258 nm，乘 1e3
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
SG_WIN  <- 9            # 二階導數：Savitzky-Golay 9 點、二次多項式
WL0     <- 2258         # 單一波長（nm）
WB_LO   <- c(1440, 1460) # 水帶區（nm）
FLAT_LO <- c(1580, 1650) # 平坦區（nm），同第 1 章估雜訊的區域
TOL_WB  <- 0.10         # 水帶比值容許範圍 +/-10 %
SAT     <- 1.5          # 最大吸光度 >= 1.5：水帶飽和，不能用這個檢查
CONC    <- c(0, 3.32, 4.87, 7.45)   # Atago 讀值（%）
PATHS   <- c(0.25, 0.5, 1, 2, 5)    # 光徑（mm）

# ---- 自我驗證的小工具：每節結尾呼叫 show_checks() ----
checks <- data.frame(ok = logical(0), item = character(0), got = numeric(0), key = numeric(0), stringsAsFactors = FALSE)
.shown <- 0
chk <- function(name, got, key, tol) {
  checks[nrow(checks) + 1, ] <<- list(abs(got - key) <= tol, name, got, key)
}
# 答案是四位有效數字（第 4 位四捨五入），容許誤差 = 第 4 位的半個單位
chk4 <- function(name, got, key) {
  tol <- if (key == 0) 1e-4 else 0.5 * 10^(floor(log10(abs(key))) - 3) * 1.001
  chk(name, got, key, tol)
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
# 0. 載入資料
# =====================================================================
hr("0. 載入資料：5 種光徑 x 4 瓶（重複 1），共 20 條光譜")
df <- read.csv(DATA_PATH)
r1 <- df[df$replicate == 1, ]
wl <- sort(unique(r1$wavelength_nm)); n_wl <- length(wl)
# 20 x n_wl 的光譜矩陣，列依「光徑、酸度」排序（先 0.25 mm 的 4 瓶，再 0.5 mm ...）
Pv <- rep(PATHS, each = 4); Cv <- rep(CONC, times = length(PATHS))
X <- t(sapply(seq_along(Pv), function(i) {
  h <- r1[abs(r1$path_mm - Pv[i]) < 1e-9 & abs(r1$acid_pct - Cv[i]) < 1e-9, ]; h <- h[order(h$wavelength_nm), ]
  stopifnot(nrow(h) == n_wl); h$absorbance
}))
row_of <- function(p, c) which(abs(Pv - p) < 1e-9 & abs(Cv - c) < 1e-9)
cat(sprintf("光譜矩陣 X：%d 條光譜 x %d 個波長（%.0f-%.0f nm）\n", nrow(X), ncol(X), min(wl), max(wl)))
cat("光徑（mm）：", PATHS, "；每個光徑的 4 瓶 Atago 酸度（%）：", CONC, "\n")
cat("本章只有一對「真正的重複」：0.25 mm 的純水量了兩次（replicate 1、2）；其餘都只量 1 次。\n")

# ---- 動手改改看（第 0 節）----
# 1. 為什麼要把判斷規則寫在檔案開頭（看結果之前），而不是看完圖再決定門檻？
# 2. 20 條光譜裡有幾對重複？這會限制我們能回答哪些問題？

# =====================================================================
# 1. 量尺：同一個光徑的重複有多像
# =====================================================================
hr("1. 量尺：0.25 mm 純水的兩次重複")
rep_w <- df[abs(df$path_mm - 0.25) < 1e-9 & df$sample == "water", ]
rep1 <- rep_w[rep_w$replicate == 1, ]; rep1 <- rep1[order(rep1$wavelength_nm), "absorbance"]
rep2 <- rep_w[rep_w$replicate == 2, ]; rep2 <- rep2[order(rep2$wavelength_nm), "absorbance"]
dif <- rep2 - rep1
rms <- sqrt(mean(dif^2)); mx <- max(abs(dif))
cat(sprintf("兩次重複相減：RMS = %.6f AU；最大絕對差 = %.6f AU\n", rms, mx))
cat("這是本資料裡「同一瓶、同一光徑」唯一的一對重複，代表儀器雜訊加上這一次的操作差異。\n")
cat("注意：只有這一對，而且是純水。重新裝樣（fill）造成的誤差、酸液樣品的重複性，都量不到。\n")

png(file.path(FIG_DIR, "ch06_replicate_overlay.png"), width = 900, height = 620, res = 110)
par(mfrow = c(2, 1), mar = c(3.5, 4, 2, 1))
plot(wl, rep1, type = "l", col = "steelblue", lwd = 2, xlab = "", ylab = "absorbance (AU)", main = "0.25 mm water, two replicates")
lines(wl, rep2, col = "darkorange", lwd = 1.5, lty = 2)
legend("topright", c("replicate 1", "replicate 2"), col = c("steelblue", "darkorange"), lty = c(1, 2), lwd = 2, cex = 0.8)
plot(wl, dif * 1e3, type = "l", col = "firebrick", xlab = "wavelength (nm)", ylab = "difference (x 1e-3 AU)",
     main = sprintf("replicate 2 - replicate 1  (RMS = %.2e AU, max = %.2e AU)", rms, mx))
abline(h = 0, col = "grey50", lty = 3)
dev.off()

chk4("replicate difference RMS (AU)", rms, 0.001236)
chk4("replicate difference max |diff| (AU)", mx, 0.004527)
show_checks("第 1 節")

# ---- 動手改改看（第 1 節）----
# 1. RMS 差 0.0012 AU 與純水在 0.25 mm 的訊號（約 0.4 AU）相比有多小？這個數字能當成「裝樣誤差」嗎？為什麼不行？
# 2. 要量得到「裝樣誤差」，實驗上至少要多做什麼？

# =====================================================================
# 2. 水帶高度檢查：規則先定好，同光徑 4 瓶互相比
# =====================================================================
hr("2. 水帶高度：同一光徑的 4 瓶，誰偏離中位數超過 10 %？")
WB <- wl >= WB_LO[1] & wl <= WB_LO[2]; FLAT <- wl >= FLAT_LO[1] & wl <= FLAT_LO[2]
wb <- list()
for (p in PATHS) {
  rr <- sapply(CONC, function(c) row_of(p, c))
  h <- sapply(rr, function(r) mean(X[r, WB]) - mean(X[r, FLAT]))
  peak <- max(sapply(rr, function(r) max(X[r, WB])))
  ratio <- h / median(h); usable <- peak < SAT
  flag <- usable & (abs(ratio - 1) > TOL_WB)
  wb[[as.character(p)]] <- list(h = h, ratio = ratio, peak = peak, usable = usable, flag = flag)
  cat(sprintf("光徑 %-4g mm：水帶最大吸光度 %.3f -> %s\n", p, peak, if (usable) "可用" else "水帶已飽和，不能用這個檢查"))
  cat(sprintf("   比值（0, 3.32, 4.87, 7.45 %%）：%s\n", paste(sprintf("%.4g", ratio), collapse = ", ")))
  if (any(flag)) cat(sprintf("   >>> 可疑：%s %% 那瓶（偏離中位數 %.2f %%）\n", paste(CONC[flag], collapse = ", "), 100 * abs(ratio[flag] - 1)))
}
cat("\n0.25 mm 的 4.87 % 那瓶比值 0.7723（水帶低 22.8 %），0.5 mm 的 3.32 % 那瓶比值 1.17（高 17 %）。\n")
cat("1 mm 以上水帶吸光度 >= 1.5：這個檢查會失靈（飽和讓水帶高度變得不隨裝樣量變），所以「沒有被標記」不等於「沒有問題」。\n")

png(file.path(FIG_DIR, "ch06_waterband_ratio.png"), width = 900, height = 520, res = 110)
par(mar = c(4, 4, 2, 1))
plot(NA, xlim = c(0.5, length(PATHS) + 0.5), ylim = c(0.7, 1.25), xaxt = "n", xlab = "path length (mm)",
     ylab = "water-band height / median of 4 bottles", main = "Water-band ratio per path (rule: |ratio - 1| <= 0.10)")
axis(1, at = seq_along(PATHS), labels = PATHS)
rect(0.5, 0.9, length(PATHS) + 0.5, 1.1, col = adjustcolor("steelblue", 0.12), border = NA)
abline(h = c(0.9, 1, 1.1), lty = c(3, 1, 3), col = "steelblue")
for (i in seq_along(PATHS)) {
  w <- wb[[as.character(PATHS[i])]]; xj <- i + c(-0.15, -0.05, 0.05, 0.15)
  pc <- ifelse(w$flag, 17, ifelse(w$usable, 19, 1)); cl <- ifelse(w$flag, "red", ifelse(w$usable, "black", "grey50"))
  points(xj, w$ratio, pch = pc, col = cl, cex = 1.4)
}
legend("topright", c("ok", "flagged (> 10 %)", "path not usable (saturated)"), pch = c(19, 17, 1), col = c("black", "red", "grey50"), cex = 0.75)
dev.off()

r25 <- wb[["0.25"]]; r5 <- wb[["0.5"]]; r1m <- wb[["1"]]
chk4("0.25 mm ratio, 0 %", r25$ratio[1], 1.045)
chk4("0.25 mm ratio, 3.32 %", r25$ratio[2], 1.014)
chk4("0.25 mm ratio, 4.87 %", r25$ratio[3], 0.7723)
chk4("0.25 mm ratio, 7.45 %", r25$ratio[4], 0.9862)
chk4("0.5 mm ratio, 0 %", r5$ratio[1], 0.9988)
chk4("0.5 mm ratio, 3.32 %", r5$ratio[2], 1.17)
chk4("0.5 mm ratio, 4.87 %", r5$ratio[3], 1.001)
chk4("0.5 mm ratio, 7.45 %", r5$ratio[4], 0.9762)
chk4("1 mm ratio, 0 %", r1m$ratio[1], 1.02)
chk4("1 mm ratio, 3.32 %", r1m$ratio[2], 1.005)
chk4("1 mm ratio, 4.87 %", r1m$ratio[3], 0.9952)
chk4("1 mm ratio, 7.45 %", r1m$ratio[4], 0.9711)
chk("0.25 mm: flagged bottle is 4.87 % (index 3)", which(r25$flag), 3, 0)
chk("0.5 mm: flagged bottle is 3.32 % (index 2)", which(r5$flag), 2, 0)
chk("1 mm: number of flags (rule not usable)", sum(r1m$flag), 0, 0)
chk4("1 mm: peak A in 1440-1460 nm", r1m$peak, 1.645)
chk("1 mm: usable (1 = yes)", as.numeric(r1m$usable), 0, 0)
show_checks("第 2 節")

# ---- 動手改改看（第 2 節）----
# 1. 把 TOL_WB 改成 0.05 或 0.20，被標記的瓶子會怎麼變？這個門檻能「看完結果再調」嗎？
# 2. 為什麼用「同光徑 4 瓶的中位數」當基準，而不是平均？（提示：中位數對一個離群值不敏感）
# 3. 1 mm 比值全都在 +/-3 % 以內，卻不能說這 4 瓶都沒問題。為什麼？

# =====================================================================
# 3. PCA 陷阱：20 條不同光徑的光譜一起做，離群值是「光徑」而不是「壞瓶」
# =====================================================================
hr("3. 自動離群偵測的陷阱：20 條原始光譜 PCA（2 個成分）+ Q 殘差")
pc <- prcomp(X, center = TRUE, scale. = FALSE)
Xc <- scale(X, center = TRUE, scale = FALSE)
Xhat <- pc$x[, 1:2, drop = FALSE] %*% t(pc$rotation[, 1:2, drop = FALSE])
Q <- rowSums((Xc - Xhat)^2)                               # Q 殘差：兩個成分解釋不了的平方和
ord <- order(-Q)                                          # Q 由大到小
lab <- sprintf("%g mm / %g %%", Pv, Cv)
ev <- pc$sdev^2 / sum(pc$sdev^2)
cat(sprintf("前兩個成分解釋的變異：%.2f %%、%.2f %%\n", 100 * ev[1], 100 * ev[2]))
cat("Q 由大到小前 8 名：\n")
for (k in 1:8) cat(sprintf("  第 %d 名：%-16s Q = %.4f\n", k, lab[ord[k]], Q[ord[k]]))
rk <- function(p, c) which(ord == row_of(p, c))
cat(sprintf("前一節標記的兩瓶：0.25 mm / 4.87 %% 排第 %d，0.5 mm / 3.32 %% 排第 %d。\n", rk(0.25, 4.87), rk(0.5, 3.32)))
cat("前三名是 2 mm 純水、2 mm 3.32 %、0.25 mm 4.87 %：2 mm 的兩條只是「光徑不同、飽和更嚴重」，不是壞瓶。\n")
cat("教訓：把不同光徑混在一起做 PCA，模型主要在描述光徑；Q 排名會被光徑的非線性主導，抓不到真正的問題瓶。\n")

png(file.path(FIG_DIR, "ch06_pca_q.png"), width = 950, height = 520, res = 110)
par(mar = c(7, 4, 2, 1))
cols <- ifelse(seq_along(Q) %in% c(row_of(0.25, 4.87), row_of(0.5, 3.32)), "red", "steelblue")
barplot(Q, names.arg = lab, las = 2, col = cols, cex.names = 0.65, ylab = "Q residual (2 PCs)", main = "PCA Q residual of 20 raw spectra (red = flagged by water-band rule)")
dev.off()

chk4("Q of 0.25 mm / 4.87 %", Q[row_of(0.25, 4.87)], 1.296)
chk("rank 1 by Q is 2 mm water (row 13)", ord[1], 13, 0)
chk("rank 2 by Q is 2 mm 3.32 % (row 14)", ord[2], 14, 0)
chk("rank 3 by Q is 0.25 mm 4.87 % (row 3)", ord[3], 3, 0)
chk("rank of 0.25 mm 4.87 % (flagged)", rk(0.25, 4.87), 3, 0)
chk("rank of 0.5 mm 3.32 % (flagged)", rk(0.5, 3.32), 8, 0)
show_checks("第 3 節")

# ---- 動手改改看（第 3 節）----
# 1. 只取 0.25 mm 的 4 條光譜做 PCA（rows 1-4），4.87 % 那條還是 Q 最大嗎？（比較 q_in_path_025：1.296 是否最高？）
# 2. 0.5 mm 的 3.32 % 為什麼只排第 8？PCA 偵測不到它，但水帶規則抓得到，差別在哪？
# 3. 為什麼說「離群值偵測的結果，必須先問『相對於誰』」？

# =====================================================================
# 4. 1 mm、2258 nm 二階導數：大殘差不等於壞光譜
# =====================================================================
hr("4. 1 mm 的 2258 nm 二階導數校正：留一法（LOO）預測誤差")
# 二階導數：SG 9 點、二次、點距為 1（與第 3、4、5 章相同）。2258 nm 離邊緣很遠，邊緣怎麼處理不影響結果。
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
D <- sg_d2(X)                                             # 20 x n_wl 二階導數
j <- which.min(abs(wl - WL0))                             # 最接近 2258 nm 的波長
m1 <- sapply(CONC, function(c) row_of(1, c))
y1 <- D[m1, j] * 1e3                                      # 1 mm 的 4 個訊號 y（x 10^-3）
atago <- CONC
titr <- c(0, 8.8, 17, 44) * c(0, rep(1.07875 * 0.24, 3))  # 滴定參考值（%）：V x 1.07875 x 0.24；純水 0

# 反向校正 y = a + b*c；LOO：每次拿掉一個點，用其餘點建線，預測被拿掉那瓶的酸度，誤差 = 預測 - 參考
fit_line <- function(cs, ys) { b <- sum((cs - mean(cs)) * (ys - mean(ys))) / sum((cs - mean(cs))^2); c(a = mean(ys) - b * mean(cs), b = b) }
loo <- function(ys, cs) {
  e <- sapply(seq_along(cs), function(i) { f <- fit_line(cs[-i], ys[-i]); (ys[i] - f["a"]) / f["b"] - cs[i] })
  list(err = as.numeric(e), rmse = sqrt(mean(e^2)))
}
r2 <- function(ys, cs) cor(cs, ys)^2

la <- loo(y1, atago); lt <- loo(y1, titr)
cat("2258 nm 訊號 y（1 mm）：", round(y1, 4), "\n")
cat("滴定酸度（%）         ：", round(titr, 3), "\n")
cat(sprintf("Atago 參考：LOO 誤差 = %s；RMSE = %.4f %%；R2 = %.4f\n", paste(sprintf("%.4g", la$err), collapse = ", "), la$rmse, r2(y1, atago)))
cat(sprintf("滴定參考  ：LOO 誤差 = %s；RMSE = %.4f %%；R2 = %.4f\n", paste(sprintf("%.4g", lt$err), collapse = ", "), lt$rmse, r2(y1, titr)))
cat("同樣 4 條光譜，換成滴定參考值，RMSE 從 2.18 降到 1.15、R2 從 0.934 升到 0.990：\n")
cat("Atago 那個最大的誤差（7.45 % 那瓶，+3.66）與其說是光譜壞了，不如說是參考值（Atago 讀值）本身比較粗。\n")

png(file.path(FIG_DIR, "ch06_loo_errors.png"), width = 900, height = 520, res = 110)
par(mar = c(4, 4, 2, 1))
ym <- range(c(la$err, lt$err)) * 1.15
plot(NA, xlim = c(0.5, 4.5), ylim = ym, xaxt = "n", xlab = "bottle (Atago acid, %)", ylab = "LOO prediction error (%)",
     main = "1 mm, 2258 nm: LOO errors, Atago vs titration reference")
axis(1, at = 1:4, labels = CONC); abline(h = 0, col = "grey50", lty = 3)
points((1:4) - 0.08, la$err, pch = 19, col = "steelblue", cex = 1.5)
points((1:4) + 0.08, lt$err, pch = 17, col = "darkorange", cex = 1.5)
legend("topleft", c(sprintf("Atago ref (RMSE %.2f)", la$rmse), sprintf("titration ref (RMSE %.2f)", lt$rmse)),
       pch = c(19, 17), col = c("steelblue", "darkorange"), cex = 0.8)
dev.off()

chk4("titration ref: 11.39 % = 44 mL value", titr[4], 11.39)
chk4("Atago LOO error, 0 %", la$err[1], 1.935)
chk4("Atago LOO error, 3.32 %", la$err[2], -0.6849)
chk4("Atago LOO error, 4.87 %", la$err[3], -1.199)
chk4("Atago LOO error, 7.45 %", la$err[4], 3.66)
chk4("Atago LOO RMSE", la$rmse, 2.182)
chk4("titration LOO error, 0 %", lt$err[1], -1.352)
chk4("titration LOO error, 2.278 %", lt$err[2], 0.8128)
chk4("titration LOO error, 4.401 %", lt$err[3], 0.285)
chk4("titration LOO error, 11.39 %", lt$err[4], -1.655)
chk4("titration LOO RMSE", lt$rmse, 1.152)
chk4("R2 Atago", r2(y1, atago), 0.9342)
chk4("R2 titration", r2(y1, titr), 0.9898)
show_checks("第 4 節")

# ---- 動手改改看（第 4 節）----
# 1. Atago 的 LOO RMSE 2.18 % 比 4 個校正點的酸度範圍（0-7.45 %）的 1/3 還大。這個校正在實用上夠用嗎？
# 2. 兩種參考值用的是同一批光譜，為什麼結果差這麼多？這說明「誤差」可能來自哪裡？
# 3. 為什麼 LOO 比單看 R2 更適合在點很少時評估預測能力？

# =====================================================================
# 5. 只拿掉一個點，R2 就會變：刪點不是修正
# =====================================================================
hr("5. Drop-one：依序拿掉 1 mm 的一個樣品，R2（Atago）怎麼變？")
drop_r2 <- sapply(1:4, function(i) r2(y1[-i], atago[-i]))
cat(sprintf("全部 4 點：R2 = %.4f\n", r2(y1, atago)))
for (i in 1:4) cat(sprintf("拿掉 %.2f %%：剩下 3 點 R2 = %.4f（%s）\n", atago[i], drop_r2[i], if (drop_r2[i] > r2(y1, atago)) "升高" else "降低"))
cat("拿掉 7.45 %（殘差最大的點）R2 升到 0.9994。\n")
cat("結論：拿掉殘差最大的點，R2 幾乎一定會升高 —— 這只是因為 3 個點更容易排成一條線，不是資料變乾淨了。\n")
cat("n = 4 時只拿掉一點就剩 3 點、自由度更少；R2 變好看不是刪資料的理由。\n")

chk4("drop-one R2: without 0 %", drop_r2[1], 0.9686)
chk4("drop-one R2: without 3.32 %", drop_r2[2], 0.942)
chk4("drop-one R2: without 4.87 %", drop_r2[3], 0.9701)
chk4("drop-one R2: without 7.45 %", drop_r2[4], 0.9994)
show_checks("第 5 節")

# ---- 動手改改看（第 5 節）----
# 1. 為什麼 3 點的 R2 幾乎總是比 4 點高？（提示：2 點的 R2 一定是 1）
# 2. 如果有人說「我把 7.45 % 那瓶刪掉後 R2 升到 0.9994，所以校正更好了」，你要怎麼回應？
# 3. 四種拿法，R2 全都比 0.9342 高（拿掉 3.32 % 只升到 0.942，拿掉 7.45 % 升到 0.9994）。「升高」本身能告訴你哪個點該刪嗎？

# =====================================================================
# 6. 0.25 mm：修正還是刪掉？三種處理，三種分數
# =====================================================================
hr("6. 0.25 mm 的 4.87 % 被水帶規則標記：三種處理方式")
m025 <- sapply(CONC, function(c) row_of(0.25, c))
y25 <- D[m025, j] * 1e3
fix <- y25; fix[3] <- fix[3] / r25$ratio[3]               # 規則事先定好：訊號除以它的水帶比值
rm_orig <- loo(y25, atago)$rmse; rm_fix <- loo(fix, atago)$rmse; rm_drop <- loo(y25[-3], atago[-3])$rmse
r2_orig <- r2(y25, atago); r2_fix <- r2(fix, atago); r2_drop <- r2(y25[-3], atago[-3])
cat("0.25 mm 訊號 y        ：", round(y25, 4), "\n")
cat("4.87 % 除以比值 0.7723 後：", round(fix, 4), "\n")
cat(sprintf("原始     ：LOO RMSE = %.4f %%；R2 = %.4f\n", rm_orig, r2_orig))
cat(sprintf("除以比值 ：LOO RMSE = %.4f %%；R2 = %.4f\n", rm_fix, r2_fix))
cat(sprintf("刪掉     ：LOO RMSE = %.4f %%；R2 = %.4f（只剩 3 點）\n", rm_drop, r2_drop))
cat("三種處理的 R2 都在 0.958-0.962，差不多；但 LOO RMSE 差很多：原始最小（1.08），修正後變差（1.53），刪掉最差（3.09）。\n")
cat("原始最好看不代表它對：標記只是「需要查原因」，不是「該刪」。這一瓶可能是裝樣量少，但我們沒有重測，不能下結論。\n")

chk4("0.25 mm LOO RMSE, original", rm_orig, 1.076)
chk4("0.25 mm LOO RMSE, divided by ratio", rm_fix, 1.527)
chk4("0.25 mm LOO RMSE, dropped", rm_drop, 3.088)
chk4("0.25 mm R2, original", r2_orig, 0.9589)
chk4("0.25 mm R2, divided by ratio", r2_fix, 0.9578)
chk4("0.25 mm R2, dropped", r2_drop, 0.9616)
show_checks("第 6 節")

# ---- 動手改改看（第 6 節）----
# 1. 三種處理的 R2 都差不多（0.958-0.962），為什麼 LOO RMSE 卻差到 3 倍？
# 2. 「除以水帶比值」這個修正規則，是在看完 LOO 結果之前還是之後決定的？如果是之後決定的，會有什麼問題？
# 3. 若這瓶真的是裝樣量不足，最正確的做法是什麼？（提示：重測）

# =====================================================================
# 7. 決策流程總結
# =====================================================================
hr("7. 決策流程：遇到可疑光譜怎麼辦")
cat("1. 先定規則：看結果之前，把判斷標準（水帶 +/-10 %、飽和門檻 1.5）寫下來。\n")
cat("2. 同類相比：只跟同光徑、同批的樣品比；不同光徑混在一起（如全部做 PCA）會被光徑主導。\n")
cat("3. 找原因：可疑的光譜先查裝樣、氣泡、光徑、飽和，不要先想刪或修。\n")
cat("4. 重測三次：同一瓶重新裝樣測 3 次，才分得出是裝樣誤差、儀器雜訊還是樣品本身。\n")
cat("5. 紀錄兩種結果：不管你刪不刪，報告「含」與「不含」可疑點的兩種結果，讓讀者自己判斷。\n")

# ---- 動手改改看（第 7 節）----
# 1. 對 0.25 mm / 4.87 % 與 0.5 mm / 3.32 % 這兩瓶，照上面五步，你會怎麼寫實驗紀錄？
# 2. 為什麼「兩種結果都要紀錄」比「只留較好看的那一種」更有說服力？

# =====================================================================
# 總結
# =====================================================================
hr("總結")
cat(sprintf("全部自我驗證：%d / %d PASS\n", sum(checks$ok), nrow(checks)))
if (!all(checks$ok)) { cat("FAIL 項目：\n"); print(checks[!checks$ok, ], row.names = FALSE) }
cat("1. 0.25 mm 純水兩次重複：RMS 0.0012 AU、最大 0.0045 AU；這是唯一的重複，量不到裝樣誤差。\n")
cat("2. 水帶規則（事先固定 +/-10 %）標出 0.25 mm / 4.87 %（比值 0.772）與 0.5 mm / 3.32 %（1.17）；1 mm 以上水帶飽和，不能用。\n")
cat("3. 20 條混光徑的 PCA Q 殘差前三名是 2 mm 兩條與 0.25 mm / 4.87 %；0.5 mm / 3.32 % 只排第 8：自動離群偵測要先問相對於誰。\n")
cat("4. 1 mm 的 2258 nm：Atago LOO RMSE 2.18、滴定 1.15；大殘差更像是參考值粗，不一定是光譜壞。\n")
cat("5. Drop-one：拿掉最大殘差，R2 升到 0.9994；刪點讓 R2 變好看，不等於資料變好。\n")
cat("6. 0.25 mm / 4.87 %：原始 / 除以比值 / 刪掉的 LOO RMSE = 1.08 / 1.53 / 3.09，R2 幾乎相同；標記不等於該刪。\n")
cat("7. 流程：先定規則、同類相比、找原因、重測三次、紀錄兩種結果。\n")
cat("圖已存到", normalizePath(FIG_DIR), "\n")
