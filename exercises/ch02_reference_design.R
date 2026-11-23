# =====================================================================
# 第 2 章 參考值與實驗設計：R 程式練習（只用 base R）
#
# 執行方式：Rscript ch02_reference_design.R （程式會自動切到本檔所在資料夾）
# 資料：../data/spectra_long.csv      圖：figures_R/ch02_*.png（圖上文字一律用英文）
# 每節結尾有「動手改改看」。第 7 節是自我驗證 PASS/FAIL 表。
# =====================================================================

options(width = 130, digits = 7)
.args <- commandArgs(trailingOnly = FALSE)
.fa <- sub("^--file=", "", .args[grep("^--file=", .args)])
if (length(.fa) == 1 && nzchar(.fa)) setwd(dirname(normalizePath(.fa)))
DATA_PATH <- "../data/spectra_long.csv"
FIG_DIR <- "figures_R"
dir.create(FIG_DIR, showWarnings = FALSE)
hr <- function(title) cat("\n", strrep("=", 70), "\n", title, "\n", strrep("=", 70), "\n", sep = "")

# ---- 可調整的參數 ----
SG_WIN     <- 9        # 二階導數 SG 視窗
F_TITER    <- 1.07875  # NaOH 力價 f
SIG_FILL   <- 0.05     # 第 5 節：填樣/光徑相對誤差（標準差）
SIG_NOISE  <- 1.2e-4   # 第 5 節：加成性雜訊（二階導數單位）
N_SIM      <- 2000     # 蒙地卡羅次數
SEED       <- 2026
TRUE_SLOPE <- -9.852e-4  # 1 mm、2258 nm 校正線（每 % 酸度）
TRUE_ICPT  <-  9.087e-3

# =====================================================================
# 0. 載入資料、定義參考值
# =====================================================================
hr("0. 載入資料與參考值")
df <- read.csv(DATA_PATH)
main <- df[df$replicate == 1, ]
wl <- sort(unique(main$wavelength_nm)); n_wl <- length(wl)
ACID_ATAGO_1102 <- c(0, 3.32, 4.87, 7.45)    # Atago 11/02
ACID_ATAGO_1103 <- c(0, 3.01, 4.12, 6.86)    # Atago 11/03
ACID_LABEL      <- c(0, 3, 5, 10)            # 瓶身標示
SAMPLES <- c("water", "vinegar_3%", "vinegar_5%", "vinegar_10%")

get_spec <- function(p) {                    # 回傳 4 x 142 矩陣（水, 3%, 5%, 10% 瓶）
  t(sapply(SAMPLES, function(s) {
    h <- main[abs(main$path_mm - p) < 1e-9 & main$sample == s, ]
    h <- h[order(h$wavelength_nm), ]
    stopifnot(nrow(h) == n_wl)
    h$absorbance
  }))
}
spec1 <- get_spec(1); spec025 <- get_spec(0.25)
cat("1 mm 光譜矩陣維度：", paste(dim(spec1), collapse = " x "), "\n")

# 手刻 SG（與第 1 章相同）
sg_coeffs2 <- function(window) {
  hh <- (window - 1) / 2; xx <- -hh:hh
  X <- cbind(1, xx, xx^2)
  as.numeric(2 * (solve(t(X) %*% X) %*% t(X))[3, ])
}
sg_apply <- function(S, coef) {
  half <- (length(coef) - 1) / 2
  out <- matrix(NA_real_, nrow(S), ncol(S))
  for (i in (half + 1):(ncol(S) - half)) out[, i] <- S[, (i - half):(i + half), drop = FALSE] %*% coef
  out
}
fit_line <- function(x, y) {
  sxx <- sum((x - mean(x))^2); sxy <- sum((x - mean(x)) * (y - mean(y)))
  slope <- sxy / sxx; icpt <- mean(y) - slope * mean(x)
  pred <- icpt + slope * x
  list(slope = slope, intercept = icpt, r2 = 1 - sum((y - pred)^2) / sum((y - mean(y))^2))
}
D2_1 <- sg_apply(spec1, sg_coeffs2(SG_WIN))
i2258 <- which.min(abs(wl - 2258)); i1670 <- which.min(abs(wl - 1670))
sig2258 <- D2_1[, i2258]; sig1670 <- D2_1[, i1670]
cat(sprintf("1 mm 2258 nm 二階導數（水,3%%,5%%,10%%）：%s\n", paste(formatC(sig2258, digits = 4, format = "e"), collapse = ", ")))
# 滴定值稍後在第 1 節算出，這裡先佔位

# ---- 動手改改看（第 0 節）----
# 1. 為什麼這章把 Atago 分成 11/02 與 11/03 兩組？兩組同一瓶的差距有多大？
# 2. 把 spec1 的列順序印出來，確認第 4 列是 10% 標示的瓶子。

# =====================================================================
# 1. 滴定計算
# =====================================================================
hr("1. 滴定計算")
# 步驟：取 10 mL 醋稀釋到 100 mL，再取 25 mL 滴定；每 mL 0.1 N NaOH 相當 0.006 g 醋酸
# 酸度 (g/100 mL) = V * f * 0.006 * (100/25) / 10 * 100
V1 <- 8.8
a1 <- V1 * F_TITER * 0.006 * (100 / 25) / 10 * 100
cat(sprintf("手算 3%% 瓶：8.8 mL * %.5f * 0.006 = %.6f g（25 mL 中）\n", F_TITER, V1 * F_TITER * 0.006))
cat(sprintf("  x 4（還原到 100 mL 稀釋液）= %.6f g；/ 10 mL 原液 *100 = %.4f g/100 mL\n",
            V1 * F_TITER * 0.006 * 4, a1))

titration <- function(V, f = F_TITER, aliquot = 25, dil_vol = 100, sample_mL = 10) {
  V * f * 0.006 * (dil_vol / aliquot) / sample_mL * 100
}
V_NAOH <- c(8.8, 17, 44)
ACID_TITR <- titration(V_NAOH)
print(data.frame(label = c("3%", "5%", "10%"), V_mL = V_NAOH, titration_pct = round(ACID_TITR, 3)), row.names = FALSE)
ACID_TITR_X <- c(0, ACID_TITR)

# ---- 動手改改看（第 1 節）----
# 1. 把 F_TITER 改成 1.0，三個結果各變成多少？比例是多少？
# 2. 若 3% 瓶實際用量是 8.08 mL（手稿筆誤），酸度會差多少？
# 3. 把 aliquot 改成 10（取 10 mL 稀釋液）時，同樣 10% 瓶要用掉多少 mL？（第 6 節會正式算）

# =====================================================================
# 2. 比值檢驗：比例誤差 vs 非比例誤差
# =====================================================================
hr("2. 滴定 / Atago(11/03) 比值：比例誤差 vs 非比例誤差")
ratio_real <- ACID_TITR / ACID_ATAGO_1103[-1]
cat("實際比值 滴定/Atago(11/03)：", round(ratio_real, 3), "\n")
cat("兩次 Atago 相差（11/02 相對 11/03）：", round(100 * (ACID_ATAGO_1102[-1] / ACID_ATAGO_1103[-1] - 1), 1), "%\n")

# 示範：力價 f 錯 5%，三瓶都乘 1.05，比值「一起」放大 1.05 倍，但彼此的相對關係不變
ACID_F5 <- titration(V_NAOH, f = F_TITER * 1.05)
cat("f 錯 5% 後的滴定值：", round(ACID_F5, 3), "\n")
cat("f 錯 5% 後的比值    ：", round(ACID_F5 / ACID_ATAGO_1103[-1], 3), "\n")
rr5 <- ACID_F5 / ACID_ATAGO_1103[-1]
cat("比值除以第一個（相對形狀）：實際", round(ratio_real / ratio_real[1], 3), "/ f 錯 5%：", round(rr5 / rr5[1], 3), "（完全相同）
")
cat("結論：比例誤差讓所有比值同乘一個常數，不會讓比值從 0.76 變到 1.66；\n")
cat("      實際比值不是常數 -> 不能用力價誤差解釋，Atago 在高濃度讀低。\n")
cat("（這是比例誤差只影響絕對值、不影響線性；R2 完全不變。）\n")
cat(sprintf("R2 檢查：真實滴定 R2 = %.4f，f 錯 5%% R2 = %.4f\n",
            fit_line(ACID_TITR_X, sig2258)$r2, fit_line(c(0, ACID_F5), sig2258)$r2))

png(file.path(FIG_DIR, "ch02_titration_vs_atago.png"), width = 1000, height = 450, res = 110)
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
plot(ACID_ATAGO_1103[-1], ACID_TITR, pch = 16, col = "red", xlim = c(0, 12), ylim = c(0, 12),
     xlab = "Atago 11/03 (%)", ylab = "titration (%)", main = "titration vs Atago")
points(ACID_ATAGO_1102[-1], ACID_TITR, pch = 17, col = "blue")
abline(0, 1, lty = 2)
legend("topleft", c("y = x", "Atago 11/03", "Atago 11/02"), lty = c(2, NA, NA), pch = c(NA, 16, 17), col = c("black", "red", "blue"), cex = 0.8)
barplot(rbind(ratio_real, ratio_real * 1.05), beside = TRUE, names.arg = c("3%", "5%", "10%"), ylab = "titration / Atago(11/03)",
        main = "ratio (red: f off by 5%)", col = c("gray70", "salmon"))
abline(h = 1, lty = 2)
dev.off()

cat("\n說明：Atago 以導電度換算酸度。醋酸是弱酸，導電度不隨濃度線性增加，且受基質影響，\n")
cat("      『可能』是高濃度讀低的原因，但不同型號的量程與換算方式不同，需查對儀器型號說明書確認。\n")

# ---- 動手改改看（第 2 節）----
# 1. 把 f 錯誤改成 -8%，比值還是不是常數？
# 2. 若 Atago 只是單純讀值「整體偏低 20%」（比例誤差），比值會長怎樣？和實際觀察一致嗎？
# 3. 兩次 Atago 讀值相差 8-15%（例如 4.87 vs 4.12），這對「用 Atago 當校正參考值」有什麼影響？

# =====================================================================
# 3. NIR 對各參考值的線性（1 mm，n = 4）
# =====================================================================
hr("3. NIR 對各參考值的線性（1 mm；2258 與 1670 nm 的二階導數）")
REFS <- list("label (0/3/5/10)" = ACID_LABEL,
             "titration" = ACID_TITR_X,
             "Atago 11/03" = ACID_ATAGO_1103,
             "Atago 11/02" = ACID_ATAGO_1102)
R2T <- data.frame(reference = names(REFS), R2_2258 = NA_real_, R2_1670 = NA_real_)
for (r in seq_along(REFS)) {
  R2T$R2_2258[r] <- fit_line(REFS[[r]], sig2258)$r2
  R2T$R2_1670[r] <- fit_line(REFS[[r]], sig1670)$r2
}
print(R2T, digits = 4, row.names = FALSE)
cat("兩個獨立吸收帶 (2258、1670 nm) 都與滴定、標示較一致；n=4 只能當『一致的跡象』，不是統計證明。\n")

png(file.path(FIG_DIR, "ch02_nir_vs_refs.png"), width = 1100, height = 440, res = 110)
par(mfrow = c(1, 2), mar = c(4, 5, 2, 1))
for (b in c("2258", "1670")) {
  sg <- if (b == "2258") sig2258 else sig1670
  plot(NA, xlim = c(0, 11.5), ylim = range(sg), xlab = "reference acidity (%)", ylab = "2nd derivative (1 mm)", main = paste(b, "nm"))
  for (i in seq_along(REFS)) lines(REFS[[i]], sg, type = "b", pch = 16, col = i)
  legend("bottomleft", c("label", "titration", "Atago 11/03", "Atago 11/02"), col = 1:4, pch = 16, lty = 1, cex = 0.8)
}
dev.off()

# ---- 動手改改看（第 3 節）----
# 1. 把 4 個參考值各自去掉水 (x = 0)，只用三瓶醋重算 R2（n = 3），排名會變嗎？
# 2. 為什麼 n = 4 時 R2 差 0.99 vs 0.93 仍不能說「統計上顯著不同」？
# 3. 兩個吸收帶獨立測量，卻都偏好滴定／標示，這支持了什麼推論？還缺什麼？

# =====================================================================
# 4. 重複性與填樣誤差
# =====================================================================
hr("4. 重複性與填樣誤差")
pop_sd <- function(v) { v <- as.numeric(v); sqrt(mean((v - mean(v))^2)) }   # 母體標準差 (np.std)
w1 <- main[abs(main$path_mm - 0.25) < 1e-9 & main$sample == "water", ]
w2 <- df[df$replicate == 2 & abs(df$path_mm - 0.25) < 1e-9 & df$sample == "water", ]
w1 <- w1[order(w1$wavelength_nm), ]; w2 <- w2[order(w2$wavelength_nm), ]
stopifnot(nrow(w1) == n_wl, nrow(w2) == n_wl)
REP_SD <- pop_sd(w1$absorbance - w2$absorbance)
cat(sprintf("0.25 mm 水兩次量測差值的母體標準差 = %.6f AU -> 儀器很穩\n", REP_SD))

i1940 <- which.min(abs(wl - 1940)); i1450 <- which.min(abs(wl - 1450))
water <- spec025[1, ]; v5 <- spec025[3, ]
R1940 <- v5[i1940] / water[i1940]; R1450 <- v5[i1450] / water[i1450]
cat(sprintf("0.25 mm：1940 nm 帶 (%.1f nm) 4.87%% 醋 %.3f vs 水 %.3f -> 比值 %.3f（%.0f%%）\n",
            wl[i1940], v5[i1940], water[i1940], R1940, 100 * (R1940 - 1)))
cat(sprintf("          1450 nm 帶 (%.1f nm) 4.87%% 醋 %.3f vs 水 %.3f -> 比值 %.3f（%.0f%%）\n",
            wl[i1450], v5[i1450], water[i1450], R1450, 100 * (R1450 - 1)))
cat("醋含水較少，水帶本來就會稍低；但報告認為『整體縮小約 15%』比預期多，疑似填樣或光徑誤差。\n")

# 填樣相對誤差約與 1/光徑 成正比：同一個絕對光徑誤差 db，相對誤差 = db / b
db <- 0.04   # 假設絕對光徑（或填樣）誤差 0.04 mm
print(data.frame(path_mm = c(0.25, 0.5, 1, 2, 5), rel_err_pct = round(100 * db / c(0.25, 0.5, 1, 2, 5), 2)), row.names = FALSE)

png(file.path(FIG_DIR, "ch02_repeat_fill.png"), width = 1100, height = 420, res = 110)
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
plot(wl, water, type = "l", ylim = range(c(water, v5)), xlab = "wavelength (nm)", ylab = "absorbance (0.25 mm)", main = "0.25 mm: water vs 4.87%")
lines(wl, v5, col = "red"); legend("topright", c("water", "4.87%"), col = c("black", "red"), lty = 1, cex = 0.8)
plot(wl, w1$absorbance - w2$absorbance, type = "l", xlab = "wavelength (nm)", ylab = "rep1 - rep2 (AU)", main = "water replicate difference")
abline(h = 0, lty = 2)
dev.off()

# ---- 動手改改看（第 4 節）----
# 1. 把 0.25 mm 水的差值換成 1 mm（沒有重複量測，想想為什麼做不到），需要補什麼實驗？
# 2. db 改成 0.01 mm，0.25 mm 與 1 mm 的相對誤差各是多少？
# 3. 4.87% 的 1450 nm 與 1940 nm 縮小幅度不同 (約 10% vs 17%)，若是純填樣誤差你預期會一樣還是不同？

# =====================================================================
# 5. 實驗設計模擬（蒙地卡羅）
# =====================================================================
hr("5. 實驗設計模擬：A 現況 vs B 稀釋系列 x 3 次填樣")
# 真實值：signal = (截距 + 斜率*c) * (1 + eps_fill) + eps_noise
cv_group <- function(x, y, grp) {                 # 留一組交叉驗證，回傳 RMSECV（% 酸度）
  e <- numeric(0)
  for (g in unique(grp)) {
    tr <- grp != g
    f <- fit_line(x[tr], y[tr])
    e <- c(e, x[!tr] - (y[!tr] - f$intercept) / f$slope)
  }
  sqrt(mean(e^2))
}
sim_design <- function(conc, nrep, nsim = N_SIM, seed = SEED) {
  set.seed(seed)
  x <- rep(conc, each = nrep); grp <- rep(seq_along(conc), each = nrep)
  slopes <- numeric(nsim); rm <- numeric(nsim)
  for (s in seq_len(nsim)) {
    y <- (TRUE_ICPT + TRUE_SLOPE * x) * (1 + rnorm(length(x), 0, SIG_FILL)) + rnorm(length(x), 0, SIG_NOISE)
    slopes[s] <- fit_line(x, y)$slope
    rm[s] <- cv_group(x, y, grp)
  }
  list(slopes = slopes, rmsecv = rm)
}
simA <- sim_design(c(0, 3.32, 4.87, 7.45), 1)
simB <- sim_design(0:7, 3)
rsdA <- sd(simA$slopes) / abs(TRUE_SLOPE); rsdB <- sd(simB$slopes) / abs(TRUE_SLOPE)
medA <- median(simA$rmsecv); medB <- median(simB$rmsecv)
print(data.frame(design = c("A: 4 conc x 1 fill", "B: 8 conc x 3 fills"),
                 n = c(4, 24), slope_RSD_pct = round(100 * c(rsdA, rsdB), 2),
                 RMSECV_median_pct = round(c(medA, medB), 3)), row.names = FALSE)
cat(sprintf("B / A 斜率相對標準差比值 = %.3f（理論約 0.48：B 的 Sxx 約為 A 的 4.3 倍）\n", rsdB / rsdA))
cat("B 的斜率不確定度明顯較小；RMSECV 的 B 用『留一濃度』（同濃度 3 次一起留下），較誠實。\n")

png(file.path(FIG_DIR, "ch02_simulation.png"), width = 1100, height = 420, res = 110)
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
rg <- range(c(simA$slopes, simB$slopes))
hist(simA$slopes, breaks = 40, xlim = rg, col = rgb(1, 0, 0, 0.4), main = "slope distribution", xlab = "fitted slope (per %)")
hist(simB$slopes, breaks = 40, col = rgb(0, 0, 1, 0.4), add = TRUE)
abline(v = TRUE_SLOPE, lty = 2); legend("topleft", c("A", "B"), fill = c(rgb(1, 0, 0, 0.4), rgb(0, 0, 1, 0.4)), cex = 0.8)
boxplot(list(A = simA$rmsecv, B = simB$rmsecv), ylab = "RMSECV (% acid)", main = "RMSECV", outline = FALSE)
dev.off()

# ---- 動手改改看（第 5 節）----
# 1. 把 SIG_FILL 改成 0.02 或 0.10，A 與 B 的斜率相對標準差怎麼變？比值 B/A 變不變？
# 2. 設計 C：8 個濃度 x 1 次填樣；與 A、B 比較，多重複與多濃度哪個貢獻較大？
# 3. 把 SIG_NOISE 調成 1e-3（雜訊大 10 倍），結論還成立嗎？

# =====================================================================
# 6. 滴定的改進清單
# =====================================================================
hr("6. 滴定的改進清單")
cat("1. 用 KHP（鄰苯二甲酸氫鉀）標定 NaOH：NaOH 吸濕、吸 CO2，不是一級標準品。力價錯屬比例誤差。\n")
cat("2. 三重複 + 空白試驗。\n")
cat("3. 10% 瓶用掉 44 mL，接近 50 mL 滴定管：改取 10 mL 稀釋液（或稀釋更多倍）。\n")
V_10 <- ACID_TITR[3] / titration(1, aliquot = 10)   # 每 mL 對應的酸度，用 10 mL 稀釋液時
cat(sprintf("   取 10 mL 稀釋液時，10%% 瓶約需 %.1f mL NaOH（原本 44 mL）\n", V_10))
cat("4. 量體積卻當成 10 g -> 結果是 w/v %；醋密度約 1.01-1.02，差異 < 2%。\n")
DENS <- 1.015
cat(sprintf("   w/w = w/v / 密度：10%% 瓶 %.2f / %.3f = %.2f（差 %.1f%%）\n",
            ACID_TITR[3], DENS, ACID_TITR[3] / DENS, 100 * (1 - 1 / DENS)))

# ---- 動手改改看（第 6 節）----
# 1. 密度改成 1.02，w/w 差多少？
# 2. 若稀釋改為 10 mL -> 250 mL、取 25 mL 滴定，10% 瓶需幾 mL？
# 3. 比例誤差 vs 非比例誤差：以上哪些改進屬於前者、哪些屬於後者？

# =====================================================================
# 7. 自我驗證
# =====================================================================
hr("7. 自我驗證")
checks <- data.frame(ok = logical(0), item = character(0), got = numeric(0), key = numeric(0), stringsAsFactors = FALSE)
chk <- function(name, got, key, tol) {
  checks[nrow(checks) + 1, ] <<- list(abs(got - key) <= tol, name, got, key)
}
chk("titration 3% (8.8 mL)", ACID_TITR[1], 2.278, 0.005)
chk("titration 5% (17 mL)", ACID_TITR[2], 4.401, 0.005)
chk("titration 10% (44 mL)", ACID_TITR[3], 11.39, 0.005)
chk("ratio titr/Atago 3%", ratio_real[1], 0.757, 0.002)
chk("ratio titr/Atago 5%", ratio_real[2], 1.068, 0.002)
chk("ratio titr/Atago 10%", ratio_real[3], 1.661, 0.002)
KEYR2 <- list("label (0/3/5/10)" = c(0.9955, 0.9983), "titration" = c(0.9898, 0.9866),
              "Atago 11/03" = c(0.9555, 0.9633), "Atago 11/02" = c(0.9342, 0.9448))
for (nm in names(KEYR2)) {
  chk(paste("R2 2258", nm), R2T$R2_2258[R2T$reference == nm], KEYR2[[nm]][1], 0.001)
  chk(paste("R2 1670", nm), R2T$R2_1670[R2T$reference == nm], KEYR2[[nm]][2], 0.001)
}
chk("replicate RMS 0.25 mm water (AU)", REP_SD, 0.001225, 1e-6)
chk("1940 nm ratio 4.87%/water", R1940, 1.052 / 1.270, 0.003)
chk("1450 nm ratio 4.87%/water", R1450, 0.557 / 0.621, 0.003)
chk("sim: B slope RSD < 0.7 x A (1 = true)", as.numeric(rsdB < 0.7 * rsdA), 1, 0)
chk("sim: B median RMSECV < A median RMSECV (1 = true)", as.numeric(medB < medA), 1, 0)
chk("10 mL aliquot: NaOH for 10% bottle (mL)", V_10, 17.6, 0.05)

cat(sprintf("%-6s | %-48s | %12s | %12s\n", "result", "item", "yours", "key")); cat(strrep("-", 88), "\n", sep = "")
for (r in seq_len(nrow(checks)))
  cat(sprintf("%-6s | %-48s | %12s | %12s\n", if (checks$ok[r]) "PASS" else "FAIL", checks$item[r],
              formatC(checks$got[r], digits = 5, format = "g"), formatC(checks$key[r], digits = 5, format = "g")))
cat(strrep("-", 88), "\n", sep = "")
cat(sprintf("TOTAL: %d / %d PASS\n", sum(checks$ok), nrow(checks)))
cat("圖已存到", normalizePath(FIG_DIR), "\n")
