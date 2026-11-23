# =====================================================================
# 第 3 章 波長挑選的陷阱：R 程式練習（只用 base R）
# 只有 4 個樣品時，從 142 個波長挑 R² 最高的會發生什麼事
#
# 執行方式：Rscript ch03_wavelength_pitfall.R （程式會自動切到本檔所在資料夾）
# 資料：../data/spectra_long.csv      圖：figures_R/ch03_*.png（圖上文字一律用英文）
# 每節結尾有「動手改改看」。第 6 節是自我驗證 PASS/FAIL 表。
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
SG_WIN   <- 9        # 二階導數 SG 視窗（二次多項式）
A_MAX    <- 1.5      # 可用波段：4 條光譜吸光度都低於此值
FIXED_WL <- 2258     # 由化學選定的固定吸收帶（nm）
N_RAND   <- 1000     # 第 5 節選做：隨機資料重複次數
SEED     <- 2026

# =====================================================================
# 0. 載入資料
# =====================================================================
hr("0. 載入資料：1 mm、重複 1、4 個樣品、142 個波長")
df <- read.csv(DATA_PATH)
d1 <- df[abs(df$path_mm - 1) < 1e-9 & df$replicate == 1, ]
conc <- sort(unique(d1$acid_pct))                       # 0, 3.32, 4.87, 7.45
wl <- sort(unique(d1$wavelength_nm)); n_wl <- length(wl)
X <- t(sapply(conc, function(a) {
  h <- d1[abs(d1$acid_pct - a) < 1e-9, ]; h <- h[order(h$wavelength_nm), ]
  stopifnot(nrow(h) == n_wl); h$absorbance
}))                                                      # 4 x 142 吸光度
cat("濃度（Atago 11/02, %）：", conc, "\n光譜矩陣維度：", paste(dim(X), collapse = " x "), "\n")

# 二階導數：SG 9 點、二次、點距為 1。邊緣與 scipy savgol_filter 預設 mode='interp' 相同：
# 邊緣用第一個（最後一個）視窗的二次擬合外推；二階導數是常數，所以等於該視窗的中心值。
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
D <- sg_d2(X)
usable <- apply(X < A_MAX, 2, all)                      # 4 條光譜都 A < 1.5
cat("可用波長數（4 條光譜 A < 1.5）：", sum(usable), "/", n_wl, "\n")

r2 <- function(a, b) cor(a, b)^2
r2_all <- function(M, yy) as.numeric(cor(M, yy))^2     # 每個波長（欄）對 yy 的 R²
best_idx <- function(R, ok = rep(TRUE, length(R))) which.max(ifelse(ok, R, -1))   # 平手取第一個
i_fixed <- which.min(abs(wl - FIXED_WL))

# ---- 動手改改看（第 0 節）----
# 1. 印出 range(X)，哪個波長區的吸光度最高？為什麼 1 mm 光徑在那裡會「飽和」？
# 2. 把 A_MAX 改成 1.2 或 1.8，可用波長數怎麼變？

# =====================================================================
# 1. 每個波長的 R²，挑最高的
# =====================================================================
hr("1. 每個波長：二階導數 vs 濃度的 R²")
R <- r2_all(D, conc)
j_all <- best_idx(R)
n_over95 <- sum(R > 0.95)
cat(sprintf("142 個波長全挑：%.1f nm，R2 = %.4f，4 條光譜吸光度 = %s\n", wl[j_all], R[j_all],
            paste(formatC(X[, j_all], digits = 3, format = "f"), collapse = ", ")))
cat(sprintf("該波長有 %d 條光譜吸光度 >= %.1f（偵測器飽和／非線性），不是可信的訊號。\n", sum(X[, j_all] >= A_MAX), A_MAX))
cat("R2 > 0.95 的波長數：", n_over95, "（n = 4 時很容易有一堆『看起來很好』的波長）\n")

png(file.path(FIG_DIR, "ch03_r2_vs_wavelength.png"), width = 1000, height = 420, res = 110)
par(mar = c(4, 4, 2, 1))
plot(wl, R, type = "l", ylim = c(0, 1), xlab = "wavelength (nm)", ylab = "R2 (2nd derivative vs concentration)",
     main = "R2 at every wavelength (n = 4)")
runs <- rle(!usable); ends <- cumsum(runs$lengths); starts <- ends - runs$lengths + 1   # 飽和區（灰色）
for (k in which(runs$values)) rect(wl[max(starts[k] - 1, 1)], -1, wl[min(ends[k] + 1, n_wl)], 2, col = rgb(0.5, 0.5, 0.5, 0.25), border = NA)
lines(wl, R)
points(wl[j_all], R[j_all], pch = 17, col = "red", cex = 1.4)
abline(h = 0.95, lty = 3)
legend("bottomleft", c("R2", "saturated (A >= 1.5)", "best of all 142"), lty = c(1, NA, NA), pch = c(NA, 15, 17),
       col = c("black", rgb(0.5, 0.5, 0.5, 0.4), "red"), cex = 0.8)
dev.off()

# ---- 動手改改看（第 1 節）----
# 1. 把門檻 0.95 改成 0.99，還剩幾個波長？
# 2. 畫出 j_all 那個波長的 4 條光譜點（plot），它們真的沿著濃度排成一直線嗎？
# 3. 飽和波段（灰色）的 R2 為什麼可以很高？這是化學還是儀器？

# =====================================================================
# 2. 只在可用波段挑，以及固定用化學帶
# =====================================================================
hr("2. 只在可用波段挑 vs 固定 2258 nm")
j_use <- best_idx(R, usable)
cat(sprintf("可用波段最佳：%.1f nm，R2 = %.4f\n", wl[j_use], R[j_use]))
cat(sprintf("固定化學帶：  %.1f nm，R2 = %.4f\n", wl[i_fixed], R[i_fixed]))
cat("挑出來的 R2 更高，但那是『從 88 個裡挑最大』的結果，不代表它真的更會預測。\n")

# ---- 動手改改看（第 2 節）----
# 1. 在可用波段裡，R2 第 2、3 高的是哪些波長？和最佳差多少？
# 2. 為什麼說 2258 nm 的 R2 比較「誠實」？（提示：它是挑之前就由化學決定的）

# =====================================================================
# 3. 洗牌檢驗：把濃度亂排，看挑波長還能不能「找到」好 R²
# =====================================================================
hr("3. 洗牌檢驗：24 種濃度排列")
perm_all <- function(v) {                               # 全部排列（4! = 24）
  if (length(v) == 1) return(list(v))
  out <- list()
  for (i in seq_along(v)) for (p in perm_all(v[-i])) out[[length(out) + 1]] <- c(v[i], p)
  out
}
perms <- perm_all(1:4)                                  # 依字典序產生，第 1 個是原本的順序
shuf <- data.frame(perm = sapply(perms, paste, collapse = ""), true = FALSE, best = NA_real_, fixed = NA_real_)
for (k in seq_along(perms)) {
  yp <- conc[perms[[k]]]
  Rp <- r2_all(D, yp)
  shuf$true[k] <- identical(perms[[k]], 1:4)
  shuf$best[k] <- max(Rp[usable])                      # 可用波段最佳
  shuf$fixed[k] <- Rp[i_fixed]                         # 固定 2258 nm
}
wrong <- !shuf$true
n_wrong_best <- sum(shuf$best[wrong] > 0.9)
n_wrong_fixed <- sum(shuf$fixed[wrong] > 0.9)
print(head(shuf, 6), digits = 4, row.names = FALSE)
cat(sprintf("錯誤排列共 %d 種；其中『可用波段最佳 R2』> 0.9 的有 %d 種；『固定 2258 nm R2』> 0.9 的有 %d 種\n",
            sum(wrong), n_wrong_best, n_wrong_fixed))
cat("濃度明明是亂排的，挑最佳波長卻幾乎總能湊出 R2 > 0.9：高 R2 在 n = 4 時沒有證明力。\n")

png(file.path(FIG_DIR, "ch03_shuffle.png"), width = 1100, height = 420, res = 110)
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
for (nm in c("best", "fixed")) {
  v <- shuf[[nm]]; o <- order(v)
  barplot(v[o], col = ifelse(shuf$true[o], "red", "gray70"), ylim = c(0, 1), ylab = "R2", xlab = "24 orderings (sorted)",
          main = if (nm == "best") "best of usable wavelengths" else "fixed 2258 nm")
  abline(h = 0.9, lty = 2)
}
dev.off()

# ---- 動手改改看（第 3 節）----
# 1. 紅色是真正的順序。「挑最佳波長」時，真正順序在 24 種裡排第幾？「固定 2258 nm」呢？
# 2. 把門檻 0.9 改成 0.95 或 0.8，兩個計數怎麼變？
# 3. 這個檢驗的邏輯是什麼？（亂排的濃度不可能有真實關係，卻能得到高 R2 -> 『高 R2』本身不能證明什麼）

# =====================================================================
# 4. 留一交叉驗證（LOOCV）的三種做法
# =====================================================================
hr("4. LOOCV 三種做法")
# 每折：留下 1 個樣品，用其餘 3 點擬合 d2 = a + b*c，再反推被留下樣品的濃度
fit3 <- function(tr, j) {
  cx <- conc[tr]; dy <- D[tr, j]
  b <- sum((cx - mean(cx)) * (dy - mean(dy))) / sum((cx - mean(cx))^2)
  c(a = mean(dy) - b * mean(cx), b = b)
}
inv_pred <- function(tr, k, j) { p <- fit3(tr, j); unname((D[k, j] - p["a"]) / p["b"]) }
rmse <- function(e) sqrt(mean(e^2))

e_out <- e_in <- e_fix <- numeric(4); pick_in <- numeric(4)
for (k in 1:4) {
  tr <- setdiff(1:4, k)
  e_out[k] <- inv_pred(tr, k, j_use) - conc[k]                     # (a) 挑波長用了全部 4 點（在迴圈外）
  jr <- best_idx(r2_all(D[tr, , drop = FALSE], conc[tr]), usable)  # (b) 每折只用 3 個訓練點重新挑
  pick_in[k] <- wl[jr]
  e_in[k] <- inv_pred(tr, k, jr) - conc[k]
  e_fix[k] <- inv_pred(tr, k, i_fixed) - conc[k]                   # (c) 固定 2258 nm
}
cv <- data.frame(left_out = conc, outside_err = e_out, inside_pick_nm = pick_in, inside_err = e_in, fixed_err = e_fix)
print(cv, digits = 4, row.names = FALSE)
RM <- c(outside = rmse(e_out), inside = rmse(e_in), fixed = rmse(e_fix))
cat(sprintf("RMSECV：挑波長在迴圈外 %.3f %%（偷看了答案）｜挑波長在迴圈內 %.3f %%（誠實）｜固定 2258 nm %.3f %%\n",
            RM["outside"], RM["inside"], RM["fixed"]))
cat("在迴圈外挑波長，被留下的樣品其實已經參與了『挑』的過程 -> 資訊洩漏，誤差被低估約 6 倍。\n")
cat("留下 4.87 % 時，挑到的波長 (1796.1 nm) 預測差最多（誤差約 -4.2 %），主導了誠實的 RMSECV。\n")

png(file.path(FIG_DIR, "ch03_loocv.png"), width = 900, height = 420, res = 110)
par(mar = c(4, 4, 2, 1))
barplot(rbind(e_out, e_in, e_fix), beside = TRUE, names.arg = paste0("leave out ", conc), col = c("gray70", "red", "steelblue"),
        ylab = "prediction - true (% acid)", main = "LOOCV errors")
abline(h = 0)
legend("topleft", c(sprintf("select outside loop (%.2f)", RM["outside"]), sprintf("select inside loop (%.2f)", RM["inside"]),
                    sprintf("fixed 2258 nm (%.2f)", RM["fixed"])), fill = c("gray70", "red", "steelblue"), cex = 0.8)
dev.off()

# ---- 動手改改看（第 4 節）----
# 1. 為什麼「誠實」版的每折會挑到不同的波長？這說明挑出來的波長有多穩定？
# 2. 固定 2258 nm 的 RMSECV (2.18 %) 與誠實的 2.58 % 接近。這說明了什麼？
# 3. 若你在迴圈外挑波長，報告 RMSECV = 0.43 %，審稿人該不該相信？

# =====================================================================
# 5. 學生練習：自己寫誠實 LOOCV、改變挑選池
# =====================================================================
hr("5. 學生練習")
cat("練習 A：把挑選池改成全部 142 個波長（含飽和），重做『迴圈內挑波長』的 LOOCV。\n")
loocv_inside <- function(ok) {
  e <- numeric(4)
  for (k in 1:4) {
    tr <- setdiff(1:4, k)
    jr <- best_idx(r2_all(D[tr, , drop = FALSE], conc[tr]), ok)
    e[k] <- inv_pred(tr, k, jr) - conc[k]
  }
  rmse(e)
}
RM_usable <- loocv_inside(usable)
RM_all <- loocv_inside(rep(TRUE, n_wl))
cat(sprintf("迴圈內挑波長：只用可用波段 RMSECV = %.3f %%；用全部 142 個 RMSECV = %.3f %%\n", RM_usable, RM_all))

cat("\n選做：把 4 x 142 的光譜換成 142 個互相獨立的隨機 N(0,1) 欄位（『142 個毫無關係的隨機波長』），\n")
cat("      重複", N_RAND, "次，看『挑最高 R2』的中位數。\n")
set.seed(SEED)
max_r2_rand <- replicate(N_RAND, max(r2_all(matrix(rnorm(4 * n_wl), 4, n_wl), conc)))
med_rand <- median(max_r2_rand)
cat(sprintf("隨機資料「最大 R2」的中位數 = %.3f（亂數不同，數字會有些微差異）\n", med_rand))
cat("這些欄位與濃度完全無關，卻挑得到 R2 約 0.99：n = 4 時，『最大 R2』只反映挑了多少個候選。\n")

# ---- 動手改改看（第 5 節）----
# 1. 把 142 改成 10 個隨機欄位，「最大 R2」的中位數變多少？
# 2. 把樣品數 4 改成 12（自己造 12 個濃度），隨機資料的最大 R2 還有 0.99 嗎？
# 3. 寫一個函式：給定挑選池與樣品數，回傳誠實 LOOCV 的 RMSECV，並比較不同池的結果。

# =====================================================================
# 6. 自我驗證
# =====================================================================
hr("6. 自我驗證")
checks <- data.frame(ok = logical(0), item = character(0), got = numeric(0), key = numeric(0), stringsAsFactors = FALSE)
chk <- function(name, got, key, tol) {
  checks[nrow(checks) + 1, ] <<- list(abs(got - key) <= tol, name, got, key)
}
chk("usable wavelengths (A < 1.5)", sum(usable), 88, 0)
chk("best of all 142: wavelength (nm)", wl[j_all], 2309.7, 0.05)
chk("best of all 142: R2", R[j_all], 0.9985, 0.001)
chk("best of all 142: spectra with A >= 1.5", sum(X[, j_all] >= A_MAX), 2, 0)
chk("best of all 142: A of 10% bottle", X[4, j_all], 1.575, 0.001)
chk("wavelengths with R2 > 0.95", n_over95, 26, 0)
chk("best of usable: wavelength (nm)", wl[j_use], 1862.4, 0.05)
chk("best of usable: R2", R[j_use], 0.9957, 0.001)
chk("fixed 2258: wavelength (nm)", wl[i_fixed], 2258.1, 0.05)
chk("fixed 2258: R2", R[i_fixed], 0.9342, 0.001)
chk("shuffle: number of orderings", nrow(shuf), 24, 0)
chk("shuffle: wrong orderings, best R2 > 0.9", n_wrong_best, 18, 0)
chk("shuffle: wrong orderings, fixed R2 > 0.9", n_wrong_fixed, 1, 0)
chk("RMSECV, selection outside loop (%)", RM["outside"], 0.432, 0.002)
chk("RMSECV, selection inside loop (%)", RM["inside"], 2.579, 0.002)
chk("RMSECV, fixed 2258 nm (%)", RM["fixed"], 2.182, 0.002)
chk("inside-loop picks match (of 4 folds)", sum(abs(pick_in - c(1853.9, 1871.0, 1796.1, 2233.1)) < 0.05), 4, 0)
chk("outside errors: max abs diff to key", max(abs(e_out - c(-0.704, 0.409, -0.098, -0.271))), 0, 0.002)
chk("inside errors: max abs diff to key", max(abs(e_in - c(-1.147, 0.694, -4.169, 2.725))), 0, 0.002)
chk("fixed errors: max abs diff to key", max(abs(e_fix - c(1.935, -0.685, -1.199, 3.66))), 0, 0.002)
chk("optional: random median max R2 (approx 0.99)", med_rand, 0.99, 0.01)

cat(sprintf("%-6s | %-48s | %12s | %12s\n", "result", "item", "yours", "key")); cat(strrep("-", 88), "\n", sep = "")
for (r in seq_len(nrow(checks)))
  cat(sprintf("%-6s | %-48s | %12s | %12s\n", if (checks$ok[r]) "PASS" else "FAIL", checks$item[r],
              formatC(checks$got[r], digits = 5, format = "g"), formatC(checks$key[r], digits = 5, format = "g")))
cat(strrep("-", 88), "\n", sep = "")
cat(sprintf("TOTAL: %d / %d PASS\n", sum(checks$ok), nrow(checks)))
cat("圖已存到", normalizePath(FIG_DIR), "\n")
