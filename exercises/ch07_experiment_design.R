# =====================================================================
# 第 7 章 設計下一次實驗：R 程式練習（只用 base R）
# 這一章不收新資料，而是用現有資料算出「下一次實驗該怎麼設計」的各種數字：
#   現在的校正有多不確定、雜訊多大、預試驗要做幾次、濃度怎麼排、測量順序怎麼排、課堂怎麼分組。
# 所有計算都是確定性的（沒有任何亂數）；最後一節產生實驗用的「順序表」與「紀錄表」。
#
# 執行方式：Rscript ch07_experiment_design.R （程式會自動切到本檔所在資料夾）
# 資料：../data/spectra_long.csv      圖：figures_R/ch07_*.png（圖上文字一律用英文）
# 輸出：ch07_run_sheet.csv（實驗順序表，UTF-8 含 BOM，Excel 可直接開）
# 每節結尾有「動手改改看」，每節最後有該節的自我驗證 PASS/FAIL。
#
# 事先固定的設定（與第 3、4、5 章相同）：
#   二階導數：Savitzky-Golay 9 點、二次多項式、點距 1，取 2258 nm，乘 1e3
#   校正：1 mm、重複 1 的 4 瓶，反向校正 y = a + b*c
#   參考值：Atago 讀值 (0, 3.32, 4.87, 7.45) %；滴定值 = 0 與 V x 1.07875 x 0.24 (V = 8.8, 17, 44 mL)
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
SG_WIN <- 9
WL0    <- 2258
CONC   <- c(0, 3.32, 4.87, 7.45)       # Atago 讀值（%）

# ---- 自我驗證的小工具：每節結尾呼叫 show_checks() ----
checks <- data.frame(ok = logical(0), item = character(0), got = numeric(0), key = numeric(0), stringsAsFactors = FALSE)
.shown <- 0
# 容許誤差：相對 1e-3；答案是 0 時用絕對 1e-3（答案只有四位有效數字）
chk <- function(name, got, key) {
  tol <- if (key == 0) 1e-3 else 1e-3 * abs(key)
  checks[nrow(checks) + 1, ] <<- list(abs(got - key) <= tol, name, as.numeric(got), key)
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
# 0. 載入資料與二階導數
# =====================================================================
hr("0. 載入資料與二階導數")
df <- read.csv(DATA_PATH)
r1 <- df[df$replicate == 1 & abs(df$path_mm - 1) < 1e-9, ]
wl <- sort(unique(r1$wavelength_nm))
X1 <- t(sapply(CONC, function(c) {                     # 4 x n_wl：1 mm、重複 1 的 4 瓶原始光譜
  h <- r1[abs(r1$acid_pct - c) < 1e-9, ]; h <- h[order(h$wavelength_nm), ]
  stopifnot(nrow(h) == length(wl)); h$absorbance
}))
sg_coeffs2 <- function(window) {
  hh <- (window - 1) / 2; xx <- -hh:hh
  M <- cbind(1, xx, xx^2)
  as.numeric(2 * (solve(t(M) %*% M) %*% t(M))[3, ])
}
sg_d2 <- function(S, window = SG_WIN) {                # 二階導數；2200-2300 nm 離邊緣很遠，邊緣怎麼處理不影響結果
  cf <- sg_coeffs2(window); half <- (window - 1) / 2
  out <- matrix(NA_real_, nrow(S), ncol(S))
  for (i in (half + 1):(ncol(S) - half)) out[, i] <- S[, (i - half):(i + half), drop = FALSE] %*% cf
  out[, 1:half] <- out[, half + 1]
  out[, (ncol(S) - half + 1):ncol(S)] <- out[, ncol(S) - half]
  out
}
D1 <- sg_d2(X1)
j <- which.min(abs(wl - WL0))
sig <- D1[, j] * 1e3                                    # 1 mm 4 瓶的 2258 nm 訊號 y（x 10^-3）
atago <- CONC
titr <- c(0, 8.8, 17, 44) * c(0, rep(1.07875 * 0.24, 3))   # 滴定參考值（%）
cat(sprintf("1 mm、重複 1：%d 瓶 x %d 個波長；2258 nm 取第 %d 點（%.1f nm）\n", nrow(X1), length(wl), j, wl[j]))
cat("訊號 y（x1e-3）：", round(sig, 4), "\n")
cat("Atago 酸度（%）：", atago, "；滴定酸度（%）：", round(titr, 3), "\n")

# ---- 動手改改看（第 0 節）----
# 1. 為什麼這章所有計算都不用亂數？（提示：讓每個人算出同一個答案，才能互相核對）
# 2. 這章用的資料只有 4 瓶、各量 1 次。下面哪些數字會因此「估不準」？

# =====================================================================
# 1. 現在的 1 mm 校正：斜率 b、殘差 s、換算成濃度的 s/|b|
# =====================================================================
hr("1. 現在的 1 mm 校正（4 點、自由度只有 2）")
fit_line <- function(xs, ys) {
  n <- length(xs); sxx <- sum((xs - mean(xs))^2)
  b <- sum((xs - mean(xs)) * (ys - mean(ys))) / sxx; a <- mean(ys) - b * mean(xs)
  s <- sqrt(sum((ys - a - b * xs)^2) / (n - 2))
  list(a = a, b = b, s = s, s_conc = s / abs(b), sxx = sxx)
}
fa <- fit_line(atago, sig); ft <- fit_line(titr, sig)
cat(sprintf("Atago 參考：b = %.4f，s = %.4f，s/|b| = %.4f %%（Sxx = %.2f）\n", fa$b, fa$s, fa$s_conc, fa$sxx))
cat(sprintf("滴定參考  ：b = %.4f，s = %.4f，s/|b| = %.4f %%（Sxx = %.2f）\n", ft$b, ft$s, ft$s_conc, ft$sxx))
cat("s/|b| 是「殘差換算成酸度」的大小：訊號每差 s，濃度就差 s/|b| %。\n")
cat("兩種參考用的是同一批光譜，s/|b| 卻差很多，所以這個殘差裡混了參考值本身的誤差。\n")

chk("Atago slope b", fa$b, -0.9852)
chk("Atago residual s", fa$s, 0.9975)
chk("Atago s/|b| (%)", fa$s_conc, 1.013)
chk("titration slope b", ft$b, -0.6415)
chk("titration residual s", ft$s, 0.3927)
chk("titration s/|b| (%)", ft$s_conc, 0.612)
show_checks("第 1 節")

# ---- 動手改改看（第 1 節）----
# 1. 為什麼斜率是負的？（二階導數在吸收峰是負的，酸越多越負）
# 2. 4 點直線的殘差自由度只有 2，s 的估計有多可靠？第 3 節會回答。

# =====================================================================
# 2. 掃描雜訊：唯一的一對重複（0.25 mm 水）
# =====================================================================
hr("2. 掃描雜訊：0.25 mm 水的兩次重複")
rw <- df[abs(df$path_mm - 0.25) < 1e-9 & df$sample == "water", ]
rp1 <- rw[rw$replicate == 1, ]; rp1 <- rp1[order(rp1$wavelength_nm), "absorbance"]
rp2 <- rw[rw$replicate == 2, ]; rp2 <- rp2[order(rp2$wavelength_nm), "absorbance"]
Rr <- sg_d2(rbind(rp1, rp2)) * 1e3
win <- wl >= 2200 & wl <= 2300
dd <- (Rr[2, ] - Rr[1, ])[win]
sd025 <- sqrt(mean((dd - mean(dd))^2)) / sqrt(2)         # 兩次相減的標準差 / sqrt(2) = 單次掃描雜訊
A1 <- mean(X1[, j]); A025 <- rp1[j]
scale <- 10^(A1 - A025)
conc025 <- sd025 / abs(fa$b); conc1_est <- sd025 * scale / abs(fa$b)
cat(sprintf("0.25 mm 單次掃描雜訊（二階導數 x1e3，2200-2300 nm）：%.4f -> 換算成酸度 %.4f %%\n", sd025, conc025))
cat(sprintf("2258 nm 的吸光度：1 mm 平均 %.4f、0.25 mm 水 %.4f -> 放大倍數 10^(A1 - A025) = %.3f\n", A1, A025, scale))
cat(sprintf("粗估 1 mm 的掃描雜訊：約 %.4f %%（數量級估計）；校正殘差 s/|b| 是 %.3f %%（Atago）/ %.3f %%（滴定）\n",
            conc1_est, fa$s_conc, ft$s_conc))
cat("注意：1 mm 的數字是「假設透光率的雜訊固定」才成立的數量級估計，不是實測。\n")
cat("結論：即使放大了，掃描雜訊（約 0.17 %）仍遠小於校正殘差（0.6-1.0 %）。\n")
cat("所以現在校正的殘差主要來自別處：裝樣、配製、參考值。這正是下一次實驗要量的。\n")

chk("scan noise sd (0.25 mm, x1e-3)", sd025, 0.02458)
chk("scan noise in acid units (0.25 mm, %)", conc025, 0.02495)
chk("A1: mean 1 mm A at 2258 nm", A1, 1.339)
chk("A025: 0.25 mm water A at 2258 nm", A025, 0.5122)
chk("scale = 10^(A1 - A025)", scale, 6.715)
chk("1 mm scan noise estimate (%)", conc1_est, 0.1675)
show_checks("第 2 節")

# ---- 動手改改看（第 2 節）----
# 1. 為什麼說 1 mm 的 0.17 % 只是「數量級」？它依賴哪個假設？
# 2. 掃描雜訊比校正殘差小很多，這對「同一次裝樣重複掃描很多次」划不划算有什麼意思？

# =====================================================================
# 3. 預試驗：只重複 k 次，標準差 s 有多不可靠
# =====================================================================
hr("3. 預試驗：重新填樣 k 次得到標準差 s，真正的 sigma 落在哪裡？")
cat("95 % 區間：sigma 介於 lo * s 與 hi * s 之間，lo = sqrt((k-1)/chi2_0.975)，hi = sqrt((k-1)/chi2_0.025)，自由度 k-1。\n")
KS <- c(2, 3, 5, 10, 20)
pl <- data.frame(k = KS, lo = sqrt((KS - 1) / qchisq(.975, KS - 1)), hi = sqrt((KS - 1) / qchisq(.025, KS - 1)))
for (i in seq_along(KS)) cat(sprintf("k = %2d：sigma 在 %.3f s 到 %.3f s 之間\n", pl$k[i], pl$lo[i], pl$hi[i]))
cat("現在的校正殘差自由度只有 2，相當於 k = 3 那一列：真正的 sigma 可能是 s 的 0.52 倍到 6.3 倍。\n")
cat("預試驗（約 10 次）可以把這個範圍縮到 0.69 倍到 1.8 倍。\n")

png(file.path(FIG_DIR, "ch07_pilot_factors.png"), width = 900, height = 520, res = 110)
par(mar = c(4, 4, 2, 1))
plot(pl$k, pl$hi, type = "b", pch = 19, col = "firebrick", log = "y", ylim = c(0.3, 40), xlab = "number of repeats k",
     ylab = "factor on observed s (95 %)", main = "How far true sigma can be from observed s")
lines(pl$k, pl$lo, type = "b", pch = 17, col = "steelblue"); abline(h = 1, lty = 3, col = "grey40")
legend("topright", c("upper factor (hi)", "lower factor (lo)"), pch = c(19, 17), col = c("firebrick", "steelblue"), cex = 0.8)
dev.off()

for (i in seq_along(KS)) {
  chk(sprintf("pilot k=%d lo", KS[i]), pl$lo[i], c(0.4461, 0.5207, 0.5991, 0.6878, 0.7605)[i])
  chk(sprintf("pilot k=%d hi", KS[i]), pl$hi[i], c(31.91, 6.285, 2.874, 1.826, 1.461)[i])
}
show_checks("第 3 節")

# ---- 動手改改看（第 3 節）----
# 1. 只重複 2 次（k = 2）的 s，真 sigma 可能大到 s 的幾倍？這樣的預試驗有用嗎？
# 2. k 從 10 加到 20，hi 從 1.83 降到 1.46；再多做 10 次值得嗎？

# =====================================================================
# 4. 區間寬度敏感度：寬度取決於還不知道的誤差
# =====================================================================
hr("4. 預測區間的半寬（Fieller）與誤差大小")
half_width <- function(s_conc, m, xs, b = ft$b, x0 = 4) {   # 截距 a = 0；未知樣品 4 %，量 m 次
  n <- length(xs); s <- s_conc * abs(b); sxx <- sum((xs - mean(xs))^2); t <- qt(.975, n - 2)
  k <- t * t * s * s; xb <- mean(xs); y0 <- b * x0
  A <- b * b - k / sxx; B <- -2 * b * y0 + 2 * k * xb / sxx; C <- y0 * y0 - k * (1 / m + 1 / n) - k * xb * xb / sxx
  r <- sqrt(B * B - 4 * A * C); rt <- sort(c((-B - r) / (2 * A), (-B + r) / (2 * A)))
  (rt[2] - rt[1]) / 2
}
X83 <- rep(0:7, each = 3); X8 <- 0:7
SS <- c(0.1, 0.2, 0.4, 0.6, 1.0)
m1 <- sapply(SS, function(s) half_width(s, 1, X83)); m3 <- sapply(SS, function(s) half_width(s, 3, X83))
cat("(a) 把 24 次填樣當成彼此獨立（只有配製誤差為 0 時才成立）；半寬（%）\n")
for (i in seq_along(SS)) cat(sprintf("   填樣誤差 s = %.1f %%：未知樣品量 1 次 -> %.4f；量 3 次 -> %.4f\n", SS[i], m1[i], m3[i]))
hw_levels <- function(s_fill, s_ref, m = 3) {              # 同一濃度 3 次填樣共用配製誤差 s_ref；用 8 個濃度平均（自由度 6）
  s_pt <- sqrt(s_ref^2 + s_fill^2 / 3)
  n <- 8; b <- ft$b; xs <- X8; s <- s_pt * abs(b); sxx <- sum((xs - mean(xs))^2); t <- qt(.975, n - 2)
  k <- t * t * s * s; xb <- mean(xs); y0 <- b * 4; w_unk <- (s_fill^2 / m) / s_pt^2
  A <- b * b - k / sxx; B <- -2 * b * y0 + 2 * k * xb / sxx; C <- y0 * y0 - k * (w_unk + 1 / n) - k * xb * xb / sxx
  r <- sqrt(B * B - 4 * A * C); rt <- sort(c((-B - r) / (2 * A), (-B + r) / (2 * A)))
  (rt[2] - rt[1]) / 2
}
SF <- 0.4; SR <- c(0, 0.2, 0.4)
hw <- sapply(SR, function(r) hw_levels(SF, r))
cat("(b) 分析 8 個濃度平均（自由度 6），不論配製誤差是否為 0 都成立；同一濃度的 3 次填樣共用配製誤差 s_ref；填樣誤差 s_fill = 0.4 %\n")
for (i in seq_along(SR)) cat(sprintf("   s_ref = %.1f %% -> 半寬 %.4f %%\n", SR[i], hw[i]))
fills <- c(1, 3, 6, 12); pt_sd <- sqrt(0.4^2 + SF^2 / fills)
cat("對照：把 24 次填樣當獨立點，s = 0.4 % 時半寬 +/- 0.51 %，只在配製誤差為 0 時才對；分析 8 個濃度平均則 s_ref = 0 時 +/- 0.60 %、s_ref = 0.4 % 時 +/- 0.71 %，兩種情況都有效。\n")
cat("(c) 填樣次數加倍能不能壓過配製誤差？s_ref = 0.4 %，每個濃度平均 f 次填樣的點標準差：\n")
for (i in seq_along(fills)) cat(sprintf("   f = %2d -> %.4f %%\n", fills[i], pt_sd[i]))
cat("多填樣只能壓低填樣那一項；配製誤差 0.4 % 是下限：0.5657 -> 0.4163，越來越慢。\n")

png(file.path(FIG_DIR, "ch07_point_sd_vs_fills.png"), width = 900, height = 520, res = 110)
par(mar = c(4, 4, 2, 1))
plot(fills, pt_sd, type = "b", pch = 19, col = "steelblue", log = "x", ylim = c(0.3, 0.62), xlab = "fills per level f (log scale)",
     ylab = "point SD (%)", main = "Point SD vs number of fills (s_fill = s_ref = 0.4 %)")
abline(h = 0.4, lty = 2, col = "firebrick"); text(8, 0.41, "floor = s_ref", col = "firebrick", cex = 0.8)
dev.off()

for (i in seq_along(SS)) chk(sprintf("half-width m=1, s=%.1f", SS[i]), m1[i], c(0.2119, 0.424, 0.8498, 1.279, 2.156)[i])
for (i in seq_along(SS)) chk(sprintf("half-width m=3, s=%.1f", SS[i]), m3[i], c(0.1274, 0.2548, 0.5107, 0.7688, 1.296)[i])
for (i in seq_along(SR)) chk(sprintf("half-width 8 levels, s_ref=%.1f", SR[i]), hw[i], c(0.6033, 0.6307, 0.7086)[i])
for (i in seq_along(fills)) chk(sprintf("point SD, fills=%d", fills[i]), pt_sd[i], c(0.5657, 0.4619, 0.432, 0.4163)[i])
show_checks("第 4 節")

# ---- 動手改改看（第 4 節）----
# 1. 填樣誤差 s 加倍，半寬大約也加倍嗎？（看 m1 那一列）
# 2. 未知樣品量 3 次比量 1 次，半寬少多少？為什麼不是少到 1/3？
# 3. 配製誤差 s_ref 等於填樣誤差時，半寬比 s_ref = 0 大多少？這提醒我們實驗要先量什麼？

# =====================================================================
# 5. 濃度怎麼排：同樣 24 次測量
# =====================================================================
hr("5. 濃度排法：斜率的標準誤與線性檢查")
designs <- list(now = atago, ends = rep(c(0, 7), each = 12), three = rep(c(0, 3.5, 7), each = 8), eight = X83)
sxx_d <- sapply(designs, function(v) sum((v - mean(v))^2))
sp <- lapply(names(designs), function(k) {
  v <- designs[[k]]; lv <- length(unique(v))
  c(n = length(v), sxx = sxx_d[[k]], se_rel = sqrt(sxx_d[["now"]] / sxx_d[[k]]), levels = lv, lof_df = lv - 2, pure_df = length(v) - lv)
}); names(sp) <- names(designs)
cat(sprintf("%-6s %4s %8s %8s %7s %7s %8s\n", "design", "n", "Sxx", "se_rel", "levels", "lof_df", "pure_df"))
for (k in names(sp)) cat(sprintf("%-6s %4d %8.2f %8.4f %7d %7d %8d\n", k, as.integer(sp[[k]]["n"]), sp[[k]]["sxx"], sp[[k]]["se_rel"],
                                 as.integer(sp[[k]]["levels"]), as.integer(sp[[k]]["lof_df"]), as.integer(sp[[k]]["pure_df"])))
cat("se_rel：斜率標準誤相對於現在 4 點各 1 次（越小越好）。ends 最小，但只有 2 個濃度，沒辦法檢查線性（lof_df = 0）。\n")
cat("線性檢查（缺適性檢定）需要：濃度數 - 2 > 0，而且純誤差自由度 > 0。eight 兩項都夠（6 與 16）。\n")
cx <- seq(0, 7, length.out = 29); curve <- cx - 0.04 * cx * (7 - cx)     # 彎曲的真實反應（示意）
mid_bias <- 0.04 * 3.5 * 3.5
cat(sprintf("示意：真實反應有彎曲時，中間最多偏 %.2f %%；只量兩端的人完全看不到。\n", mid_bias))

png(file.path(FIG_DIR, "ch07_spacing.png"), width = 1000, height = 480, res = 110)
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
se <- sapply(sp, function(z) z[["se_rel"]]); lofd <- sapply(sp, function(z) z[["lof_df"]])
bp <- barplot(se, col = c("grey60", "steelblue", "steelblue", "steelblue"), ylim = c(0, 1.15), ylab = "slope SE relative to 'now'", main = "Design vs slope SE")
text(bp, se + 0.05, sprintf("%.2f", se), cex = 0.8); text(bp, 0.08, paste0("lof df=", lofd), cex = 0.7)
plot(cx, curve, type = "l", lwd = 2, col = "firebrick", xlab = "true acid (%)", ylab = "response (acid units)", main = "Curved response, 2 end points only")
abline(0, 1, lty = 2, col = "grey40"); points(c(0, 7), c(0, 7), pch = 19, cex = 1.3)
legend("topleft", c("curved truth", "straight line through ends"), lty = c(1, 2), col = c("firebrick", "grey40"), cex = 0.75)
dev.off()

for (k in names(sp)) {
  chk(sprintf("%s: Sxx", k), sp[[k]][["sxx"]], c(now = 29.09, ends = 294, three = 196, eight = 126)[[k]])
  chk(sprintf("%s: slope SE relative", k), sp[[k]][["se_rel"]], c(now = 1, ends = 0.3146, three = 0.3852, eight = 0.4805)[[k]])
  chk(sprintf("%s: lack-of-fit df", k), sp[[k]][["lof_df"]], c(now = 2, ends = 0, three = 1, eight = 6)[[k]])
  chk(sprintf("%s: pure-error df", k), sp[[k]][["pure_df"]], c(now = 0, ends = 22, three = 21, eight = 16)[[k]])
}
chk("curve: mid-point bias (%)", mid_bias, 0.49)
show_checks("第 5 節")

# ---- 動手改改看（第 5 節）----
# 1. ends 的斜率最準，為什麼我們仍然不選它？
# 2. 現在的 4 點（now）lof_df = 2、pure_df = 0：能檢查線性嗎？能估純誤差嗎？
# 3. 8 濃度 x 3 次的 se_rel = 0.48，只比 now 好 2 倍，卻多出 16 個純誤差自由度。這值得嗎？

# =====================================================================
# 6. 測量順序：從低到高 vs 隨機 vs 分區塊
# =====================================================================
hr("6. 測量順序：儀器漂移（整段相當於 0.5 % 濃度）對斜率的影響")
n <- 24; drift_total <- 0.5; dr <- drift_total * (0:(n - 1)) / (n - 1)
Sxx <- sum((X83 - mean(X83))^2); Sdd <- sum((dr - mean(dr))^2)
asc <- sort(X83); asc_bias <- sum((asc - mean(asc)) * dr) / Sxx * 100
rand_sd <- sqrt(Sdd / ((n - 1) * Sxx)) * 100
sxx_b <- sum((X8 - mean(X8))^2)
block_sd <- sqrt(sum(sapply(1:3, function(i) { b_ <- dr[((i - 1) * 8 + 1):(i * 8)]; sxx_b * sum((b_ - mean(b_))^2) / 7 }))) / Sxx * 100
base <- c(3, 6, 0, 5, 2, 7, 1, 4)
example <- unlist(lapply(0:2, function(blk) base[((0:7 + 3 * blk) %% 8) + 1]))
example_bias <- sum((example - 3.5) * dr) / Sxx * 100
cat(sprintf("由低到高排：斜率偏差 = %.3f %%（系統性偏差，每次都一樣）\n", asc_bias))
cat(sprintf("完全隨機排：偏差平均為 0、標準差 = %.3f %%（95 %%：+/- %.3f %%）\n", rand_sd, 1.96 * rand_sd))
cat(sprintf("分 3 區塊（每塊 8 個濃度各 1 次、塊內隨機）：標準差 = %.4f %%（95 %%：+/- %.4f %%）\n", block_sd, 1.96 * block_sd))
cat("範例順序（3 區塊，每塊 8 濃度，拉丁方式輪替，不靠亂數）：\n")
for (b in 0:2) cat(sprintf("   區塊 %d：%s\n", b + 1, paste(example[(b * 8 + 1):(b * 8 + 8)], collapse = " ")))
cat(sprintf("這個範例順序在上述漂移下的斜率偏差 = %.4f %%\n", example_bias))
cat("結論：由低到高排會把漂移變成斜率誤差；分區塊把漂移抵銷掉大部分；範例順序的偏差只有 0.069 %。\n")

png(file.path(FIG_DIR, "ch07_order_bias.png"), width = 900, height = 520, res = 110)
par(mar = c(5, 4, 2, 1))
v <- c(asc_bias, 1.96 * rand_sd, 1.96 * block_sd, abs(example_bias))
bp <- barplot(v, names.arg = c("low to high\n(bias)", "random\n(95 % range)", "3 blocks\n(95 % range)", "example order\n(|bias|)"),
              col = c("firebrick", "grey60", "steelblue", "darkorange"), ylim = c(0, 7.5), ylab = "slope error (% acid)",
              main = "Effect of run order under drift (0.5 % over 24 runs)")
text(bp, v + 0.25, sprintf("%.2f", v), cex = 0.85)
dev.off()

chk("ascending order: slope bias (%)", asc_bias, 6.522)
chk("random order: sd (%)", rand_sd, 1.369)
chk("random order: 95 % range (%)", 1.96 * rand_sd, 2.684)
chk("3 blocks: sd (%)", block_sd, 0.4744)
chk("3 blocks: 95 % range (%)", 1.96 * block_sd, 0.9298)
chk("example order: slope bias (%)", example_bias, 0.06901)
chk("example order: mismatches with reference list", sum(example != c(3, 6, 0, 5, 2, 7, 1, 4, 5, 2, 7, 1, 4, 3, 6, 0, 1, 4, 3, 6, 0, 5, 2, 7)), 0)
show_checks("第 6 節")

# ---- 動手改改看（第 6 節）----
# 1. 由低到高排，漂移 0.5 % 就產生 6.5 % 的斜率偏差，為什麼這麼大？（提示：漂移和濃度同方向上升）
# 2. 完全隨機的偏差平均是 0，為什麼還是不如分區塊？
# 3. 把 base 換成另一種排法，example_bias 會變多少？

# =====================================================================
# 7. 課堂分組：8 組，每組 4 個濃度
# =====================================================================
hr("7. 課堂分組")
groups <- lapply(0:7, function(gi) sort(unique((gi + c(0, 1, 3, 6)) %% 8)))
per_level <- sapply(0:7, function(lv) sum(sapply(groups, function(g) lv %in% g)))
for (gi in 0:7) cat(sprintf("第 %d 組：濃度 %s %%\n", gi + 1, paste(groups[[gi + 1]], collapse = ", ")))
cat("每個濃度剛好 4 組做；每組都有高有低（位移 0, 1, 3, 6）。全班合起來等於 8 濃度 x 4 次。\n")

ref_groups <- list(c(0, 1, 3, 6), c(1, 2, 4, 7), c(0, 2, 3, 5), c(1, 3, 4, 6), c(2, 4, 5, 7), c(0, 3, 5, 6), c(1, 4, 6, 7), c(0, 2, 5, 7))
chk("groups: mismatches with reference", sum(!mapply(function(a, b) identical(as.numeric(a), as.numeric(b)), groups, ref_groups)), 0)
chk("per-level group counts: mismatches with 4", sum(per_level != 4), 0)
show_checks("第 7 節")

# ---- 動手改改看（第 7 節）----
# 1. 如果有 6 組，要怎麼改位移讓每個濃度仍然被做到相同次數？
# 2. 為什麼每組都要有低有高，而不是第 1 組做 0-3、第 2 組做 4-7？

# =====================================================================
# 8. 產生順序表與紀錄表
# =====================================================================
hr("8. 產生 ch07_run_sheet.csv（順序表）")
cols <- c("階段", "區塊", "樣品", "順序", "濃度_pct", "第幾次填樣", "第幾次掃描", "日期", "操作者", "光徑_mm", "瓶號", "溫度_C", "檔名", "備註")
pilot_rows <- sprintf("預試驗,預試驗,預試驗瓶,%d,4,%d,%d", 1:13, c(1:10, 11, 11, 11), c(rep(1, 10), 1:3))
fill_no <- ave(example, example, FUN = seq_along)
main_rows <- sprintf("主實驗,%d,標準品,%d,%d,%d,1", rep(1:3, each = 8), 1:24, as.integer(example), as.integer(fill_no))
val_seq <- c("A1", "B1", "C1", "D1", "B2", "C2", "D2", "A2", "C3", "D3", "A3", "B3")      # 字母 = 盲樣，數字 = 第幾次填樣
val_rows <- sprintf("驗證,盲樣,盲樣%s,%d,,%s,1", substr(val_seq, 1, 1), 1:12, substr(val_seq, 2, 2))
body <- c(pilot_rows, main_rows, val_rows)
txt <- paste0(paste(cols, collapse = ","), "\n", paste0(body, strrep(",", 7), "\n", collapse = ""))
con <- file("ch07_run_sheet.csv", "wb")
writeBin(c(as.raw(c(0xEF, 0xBB, 0xBF)), charToRaw(enc2utf8(txt))), con); close(con)
rs <- read.csv("ch07_run_sheet.csv", fileEncoding = "UTF-8-BOM", check.names = FALSE)
cat(sprintf("已寫出 ch07_run_sheet.csv：%d 列（預試驗 13 + 主實驗 24 + 驗證 12），%d 欄\n", nrow(rs), ncol(rs)))
cat("預試驗：同一瓶約 4 % 的醋重新填樣 10 次（填樣 1-10），再對第 11 次填樣連續掃描 3 次（掃描雜訊檢查）；主實驗：8 濃度 x 3 區塊，依範例順序。\n")
cat("驗證：4 個盲樣（盲樣A-D）各填樣 3 次，輪替順序 A1 B1 C1 D1 / B2 C2 D2 A2 / C3 D3 A3 B3；濃度欄留空（操作者不知道）。\n")
cat("日期、操作者、光徑、瓶號、溫度、檔名、備註 7 欄留白，實驗時手填。列印版的計畫與檢查表見 ch07_record_sheet.md。\n")

chk("run sheet: number of rows", nrow(rs), 49)
chk("run sheet: number of columns", ncol(rs), 14)
chk("run sheet: pilot rows", sum(rs[["階段"]] == "預試驗"), 13)
chk("run sheet: main rows", sum(rs[["階段"]] == "主實驗"), 24)
chk("run sheet: validation rows", sum(rs[["階段"]] == "驗證"), 12)
chk("run sheet: main levels mismatch with example order", sum(rs[["濃度_pct"]][rs[["階段"]] == "主實驗"] != example), 0)
pl <- rs[rs[["階段"]] == "預試驗", ]
chk("run sheet: pilot fill/scan pattern mismatches (fills 1-11, fill 11 scanned 3 times)",
    sum(pl[["第幾次填樣"]] != c(1:10, 11, 11, 11)) + sum(pl[["第幾次掃描"]] != c(rep(1, 10), 1:3)), 0)
vl <- rs[rs[["階段"]] == "驗證", ]
chk("run sheet: validation order mismatches (A1 B1 C1 D1 B2 C2 D2 A2 C3 D3 A3 B3)",
    sum(paste0(sub("盲樣", "", vl[["樣品"]]), vl[["第幾次填樣"]]) != val_seq), 0)
chk("run sheet: validation concentration cells not blank", sum(!is.na(vl[["濃度_pct"]])), 0)
show_checks("第 8 節")

# ---- 動手改改看（第 8 節）----
# 1. 用 Excel 開 ch07_run_sheet.csv，中文有亂碼嗎？（檔案含 BOM，應該不會）
# 2. 預試驗為什麼要同時有「重新填樣 10 次」與「第 11 次填樣掃描 3 次」？兩者分別量到什麼誤差？
# 3. 主實驗第 k 區塊的「第幾次填樣」等於 k，為什麼？
# 4. 驗證階段為什麼要由別人編碼盲樣、濃度欄留空？

# =====================================================================
# 總結
# =====================================================================
hr("總結")
cat(sprintf("全部自我驗證：%d / %d PASS\n", sum(checks$ok), nrow(checks)))
if (!all(checks$ok)) { cat("FAIL 項目：\n"); print(checks[!checks$ok, ], row.names = FALSE) }
cat("1. 現在的 1 mm 校正：殘差換算成酸度是 1.01 %（Atago）/ 0.61 %（滴定），自由度只有 2。\n")
cat("2. 掃描雜訊很小（0.25 mm 約 0.025 %，1 mm 粗估 0.17 %），殘差主要來自裝樣、配製與參考值。\n")
cat("3. 預試驗重複 10 次，真正的 sigma 在 0.69-1.83 倍 s 之間；只重複 3 次是 0.52-6.3 倍。\n")
cat("4. 區間半寬隨填樣誤差成比例變大；配製誤差是填樣次數壓不掉的下限。\n")
cat("5. 同樣 24 次測量：8 濃度 x 3 次可以做線性檢查（lof_df 6、pure_df 16）；只量兩端看不到彎曲。\n")
cat("6. 測量順序：由低到高偏 6.5 %，隨機 +/- 2.7 %，分 3 區塊 +/- 0.93 %，範例順序只偏 0.069 %。\n")
cat("7. 課堂 8 組各做 4 個濃度，每個濃度剛好 4 組。\n")
cat("8. 已產生順序表（49 列：預試驗 13 + 主實驗 24 + 驗證 12）；紀錄表見 ch07_record_sheet.md。\n")
cat("圖已存到", basename(FIG_DIR), "\n")
