# =====================================================================
# 第 4 章 PCA 與 PLS：R 程式練習（只用 base R）
# 用同一組 4 個樣品的光譜：先用 PCA「看」資料，再用只有 1 個成分的 PLS 做校正
#
# 執行方式：Rscript ch04_pca_pls.R （程式會自動切到本檔所在資料夾）
# 資料：../data/spectra_long.csv      圖：figures_R/ch04_*.png（圖上文字一律用英文）
# 每節結尾有「動手改改看」，每節的最後有該節的自我驗證 PASS/FAIL；第 6 節是學生練習與總結。
# 所有設定都是「事先固定」的（見下方參數），不是看 RMSECV 調出來的。
# =====================================================================

options(width = 130, digits = 7)
.args <- commandArgs(trailingOnly = FALSE)
.fa <- sub("^--file=", "", .args[grep("^--file=", .args)])
if (length(.fa) == 1 && nzchar(.fa)) setwd(dirname(normalizePath(.fa)))
DATA_PATH <- "../data/spectra_long.csv"
FIG_DIR <- "figures_R"
dir.create(FIG_DIR, showWarnings = FALSE)
hr <- function(title) cat("\n", strrep("=", 70), "\n", title, "\n", strrep("=", 70), "\n", sep = "")

# ---- 事先固定的設定（先決定、不依 RMSECV 挑選）----
SG_WIN  <- 9                          # 二階導數：Savitzky-Golay 9 點、二次多項式（同第 1 章）
BANDS   <- list(c(1480, 1872), c(2116, 2285))   # 報告中 1 mm 可用的兩段波段（nm），共 73 個波長
N_COMP  <- 1                          # PLS 成分數：只有一種成分在變（醋酸取代水），所以 1 個
SCALE   <- FALSE                      # PLS 不縮放：所有波長同單位（吸光度），縮放會放大雜訊

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
# 0. 載入資料、二階導數、波段
# =====================================================================
hr("0. 載入資料：5 個光徑 x 4 個樣品（重複 1）")
df <- read.csv(DATA_PATH)
d_r1 <- df[df$replicate == 1, ]
paths <- sort(unique(d_r1$path_mm))                      # 5 個光徑（mm）
conc <- sort(unique(d_r1$acid_pct))                      # 0, 3.32, 4.87, 7.45
wl <- sort(unique(d_r1$wavelength_nm)); n_wl <- length(wl)
get_spec <- function(p, a) {
  h <- d_r1[abs(d_r1$path_mm - p) < 1e-9 & abs(d_r1$acid_pct - a) < 1e-9, ]; h <- h[order(h$wavelength_nm), ]
  stopifnot(nrow(h) == n_wl); h$absorbance
}
# 20 條原始光譜：列依「光徑、濃度」排序
grid <- expand.grid(conc = conc, path = paths)           # conc 變動最快
path20 <- grid$path; conc20 <- grid$conc
X20 <- t(mapply(get_spec, path20, conc20))               # 20 x 142
# 1 mm 的 4 條光譜
X1 <- t(sapply(conc, function(a) get_spec(1, a)))        # 4 x 142
cat("光徑（mm）：", paths, "\n濃度（Atago 11/02, %）：", conc, "\n")
cat("20 條光譜矩陣：", paste(dim(X20), collapse = " x "), "；1 mm 光譜矩陣：", paste(dim(X1), collapse = " x "), "\n")

# 二階導數：SG 9 點、二次、點距為 1（與第 3 章相同）。邊緣用第一個／最後一個視窗的值；
# 本章只用波段內（離邊緣很遠）的波長，邊緣怎麼處理不影響結果。
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
band <- Reduce(`|`, lapply(BANDS, function(b) wl >= b[1] & wl <= b[2]))   # 兩段波段
wlb <- wl[band]
D1 <- sg_d2(X1)                                          # 4 x 142 二階導數
Db <- D1[, band]                                         # 4 x 73：二階導數、只取兩段波段
Rb <- X1[, band]                                         # 4 x 73：原始吸光度、同樣兩段波段（比較用）
cat(sprintf("波段內波長數：%d（%s）\n", sum(band), paste(sapply(BANDS, function(b) sprintf("%g-%g nm", b[1], b[2])), collapse = " 與 ")))

# ---- 工具函式：PCA（prcomp，只置中）與 PLS1（NIPALS，手寫）----
pca_fit <- function(M) prcomp(M, center = TRUE, scale. = FALSE)
evr <- function(p) p$sdev^2 / sum(p$sdev^2)               # 解釋變異比例
abs_r <- function(a, b) abs(cor(a, b))                    # PC 的正負號是任意的，所以取絕對值

# PLS1：X、y 先置中（不縮放），每個成分：w = X'y 正規化 -> t = Xw -> p = X't/t't, q = y't/t't -> 扣除（deflation）
pls1_fit <- function(X, y, nc) {
  xm <- colMeans(X); ym <- mean(y)
  Xc <- sweep(X, 2, xm); yc <- y - ym
  W <- P <- matrix(0, ncol(X), nc); Q <- numeric(nc)
  for (a in seq_len(nc)) {
    w <- crossprod(Xc, yc); w <- w / sqrt(sum(w^2))
    t <- as.numeric(Xc %*% w); tt <- sum(t^2)
    p <- as.numeric(crossprod(Xc, t)) / tt; q <- sum(yc * t) / tt
    Xc <- Xc - outer(t, p); yc <- yc - t * q
    W[, a] <- w; P[, a] <- p; Q[a] <- q
  }
  B <- as.numeric(W %*% solve(crossprod(P, W), Q))       # 迴歸係數（對置中後的 X）
  list(B = B, b0 = ym - sum(xm * B), nc = nc)
}
pls1_predict <- function(fit, Xnew) as.numeric(as.matrix(Xnew) %*% fit$B + fit$b0)

# 留一交叉驗證：回傳誤差（預測 - 真值）與 RMSECV。scale = TRUE 時每折用訓練折的平均與標準差標準化 X 與 y
rmse <- function(e) sqrt(mean(e^2))
pls_loocv <- function(M, y, nc = N_COMP, scale = SCALE) {
  n <- nrow(M); e <- numeric(n)
  for (k in seq_len(n)) {
    tr <- setdiff(seq_len(n), k)
    if (scale) {
      mu <- colMeans(M[tr, , drop = FALSE]); sdv <- apply(M[tr, , drop = FALSE], 2, sd); ym <- mean(y[tr]); ys <- sd(y[tr])
      fit <- pls1_fit(sweep(sweep(M[tr, , drop = FALSE], 2, mu), 2, sdv, "/"), (y[tr] - ym) / ys, nc)
      pr <- pls1_predict(fit, sweep(sweep(M[k, , drop = FALSE], 2, mu), 2, sdv, "/")) * ys + ym
    } else {
      pr <- pls1_predict(pls1_fit(M[tr, , drop = FALSE], y[tr], nc), M[k, , drop = FALSE])
    }
    e[k] <- pr - y[k]
  }
  list(rmsecv = rmse(e), err = e)
}

chk("number of wavelengths", n_wl, 142, 0)
chk("number of spectra (5 paths x 4 samples)", nrow(X20), 20, 0)
chk("wavelengths in the two bands", sum(band), 73, 0)
show_checks("第 0 節")

# ---- 動手改改看（第 0 節）----
# 1. 為什麼 PLS 不縮放、只取 1 個成分、波段取 1480-1872 與 2116-2285？這些設定是在看到 RMSECV 之前就決定的，為什麼這點重要？
# 2. 把 BANDS 改成只留第一段，band 內波長數是多少？（改完記得改回來再看自我驗證）

# =====================================================================
# 1. PCA：20 條原始光譜（5 光徑 x 4 樣品）
# =====================================================================
hr("1. PCA：20 條原始光譜")
p20 <- pca_fit(X20); s20 <- p20$x
ev20 <- evr(p20)
r_path <- abs_r(s20[, 1], path20); r_conc <- abs_r(s20[, 1], conc20)
cat(sprintf("解釋變異：PC1 = %.3f，PC2 = %.3f\n", ev20[1], ev20[2]))
cat(sprintf("|r(PC1, 光徑)| = %.3f；|r(PC1, 濃度)| = %.3f\n", r_path, r_conc))
cat("PC1 幾乎就是『光徑』：光徑把整條光譜等比例放大，這個「大小」變化蓋過了濃度造成的細微差異。\n")

path_cols <- hcl.colors(length(paths), "Dark 3")
png(file.path(FIG_DIR, "ch04_pca20.png"), width = 900, height = 520, res = 110)
par(mar = c(4, 4, 2, 1))
plot(s20[, 1], s20[, 2], pch = 19, cex = 1.3, col = path_cols[match(path20, paths)],
     xlab = sprintf("PC1 (%.1f%%)", 100 * ev20[1]), ylab = sprintf("PC2 (%.1f%%)", 100 * ev20[2]),
     main = "PCA of 20 raw spectra (colored by path length)")
legend("topright", legend = paste0(paths, " mm"), pch = 19, col = path_cols, cex = 0.8, title = "path")
dev.off()

chk("PCA 20 raw: explained variance PC1", ev20[1], 0.917, 0.001)
chk("PCA 20 raw: explained variance PC2", ev20[2], 0.069, 0.001)
chk("PCA 20 raw: |r(PC1, path)|", r_path, 0.91, 0.002)
chk("PCA 20 raw: |r(PC1, conc)|", r_conc, 0.022, 0.002)
show_checks("第 1 節")

# ---- 動手改改看（第 1 節）----
# 1. 畫散佈圖時改成「依濃度上色」（col = 依 conc20），看到什麼？濃度的資訊藏在哪個 PC？（試試 PC2、PC3）
# 2. 為什麼 PC1 的正負號每次執行可能不同，卻不影響結論？
# 3. 這是「描述」資料的工具（沒用到濃度），它能告訴你最大的變異來自哪裡，卻不能直接預測濃度。

# =====================================================================
# 2. PCA：1 mm、二階導數、兩段波段
# =====================================================================
hr("2. PCA：1 mm 二階導數（4 x 73）")
p1 <- pca_fit(Db); s1 <- p1$x
ev1 <- evr(p1)
r_conc1 <- abs_r(s1[, 1], conc)
cat(sprintf("解釋變異：PC1 = %.3f，PC2 = %.3f\n", ev1[1], ev1[2]))
cat(sprintf("|r(PC1, 濃度)| = %.4f\n", r_conc1))
cat("二階導數去掉了基線與光徑造成的緩慢變化，只剩吸收峰的形狀；此時 PC1 幾乎就是濃度。\n")

png(file.path(FIG_DIR, "ch04_pca_1mm.png"), width = 800, height = 480, res = 110)
par(mar = c(4, 4, 2, 1))
plot(conc, s1[, 1], pch = 19, cex = 1.4, col = "steelblue", xlab = "acid (%)", ylab = sprintf("PC1 score (%.1f%%)", 100 * ev1[1]),
     main = "PC1 of 2nd-derivative spectra vs concentration (1 mm)")
text(conc, s1[, 1], labels = conc, pos = 3, cex = 0.8)
abline(lm(s1[, 1] ~ conc), lty = 2)
dev.off()

chk("PCA 1 mm d2: explained variance PC1", ev1[1], 0.969, 0.001)
chk("PCA 1 mm d2: |r(PC1, conc)|", r_conc1, 0.9984, 0.002)
show_checks("第 2 節")

# ---- 動手改改看（第 2 節）----
# 1. 對 X1[, band]（原始吸光度、不微分）做 PCA，PC1 的解釋變異與 |r(PC1, conc)| 變成多少？
# 2. 為什麼只有 4 個樣品時，PCA 最多只有 3 個有意義的 PC？（看 p1$sdev）

# =====================================================================
# 3. PLS 留一交叉驗證（LOOCV）
# =====================================================================
hr("3. PLS LOOCV：1 mm 二階導數 4 x 73")
cat("每折：留下 1 個樣品，用其餘 3 個擬合 PLS，預測被留下的樣品；誤差 = 預測 - 真值（%）。\n")
cv1 <- pls_loocv(Db, conc, 1, FALSE)
cv2 <- pls_loocv(Db, conc, 2, FALSE)
cvr <- pls_loocv(Rb, conc, 1, FALSE)                     # 原始光譜（不微分）、同樣波段
cvs <- pls_loocv(Db, conc, 1, TRUE)                      # 縮放 X 與 y（Orange 的預設）
RM2258 <- 2.182                                          # 第 1 章：固定單一波段 2258 nm 的 RMSECV
tab <- data.frame(left_out = conc, PLS1_err = cv1$err, PLS2_err = cv2$err, raw_PLS1_err = cvr$err)
print(tab, digits = 3, row.names = FALSE)
cat(sprintf("RMSECV（%%）：PLS-1 = %.3f；PLS-2 = %.3f；原始光譜 PLS-1 = %.3f；縮放 PLS-1 = %.3f；單一波段 2258 nm（第 1 章）= %.3f\n",
            cv1$rmsecv, cv2$rmsecv, cvr$rmsecv, cvs$rmsecv, RM2258))
cat("1 個成分最好：只有一種成分在變，多加的成分只是在擬合雜訊。二階導數把 RMSECV 從 2.5 降到 0.3。\n")

fitR2 <- sapply(1:3, function(nc) { f <- pls1_fit(Db, conc, nc); pr <- pls1_predict(f, Db); 1 - sum((conc - pr)^2) / sum((conc - mean(conc))^2) })
cat(sprintf("全部 4 點擬合的 R2：1 個成分 = %.4f；2 個 = %.4f；3 個 = %.4f\n", fitR2[1], fitR2[2], fitR2[3]))
cat("3 個成分的 R2 = 1.0000：4 個樣品、3 個成分，模型可以把 4 點完全吃掉。擬合 R2 愈高不代表愈好（過度擬合），要看 RMSECV。\n")

png(file.path(FIG_DIR, "ch04_loocv.png"), width = 900, height = 420, res = 110)
par(mar = c(4, 4, 2, 1))
barplot(rbind(cv1$err, cv2$err), beside = TRUE, names.arg = paste0("leave out ", conc), col = c("red", "steelblue"),
        ylab = "prediction - true (% acid)", main = "LOOCV errors of PLS (2nd derivative, 73 wavelengths)")
abline(h = 0)
legend("topleft", c(sprintf("1 component (%.3f)", cv1$rmsecv), sprintf("2 components (%.3f)", cv2$rmsecv)), fill = c("red", "steelblue"), cex = 0.8)
dev.off()

chk("RMSECV PLS-1 (%)", cv1$rmsecv, 0.317, 0.002)
chk("PLS-1 errors: max abs diff to key", max(abs(cv1$err - c(0.46, 0.252, -0.341, -0.104))), 0, 0.002)
chk("RMSECV PLS-2 (%)", cv2$rmsecv, 0.461, 0.002)
chk("PLS-2 errors: max abs diff to key", max(abs(cv2$err - c(-0.209, 0.43, -0.406, 0.676))), 0, 0.002)
chk("RMSECV raw spectra PLS-1 (%)", cvr$rmsecv, 2.547, 0.002)
chk("RMSECV scaled PLS-1 (%)", cvs$rmsecv, 2.458, 0.002)
chk("fit R2, 1 component", fitR2[1], 0.9968, 0.002)
chk("fit R2, 2 components", fitR2[2], 0.9986, 0.002)
chk("fit R2, 3 components", fitR2[3], 1.0, 0.002)
show_checks("第 3 節")

# ---- 動手改改看（第 3 節）----
# 1. 為什麼 2 個成分的 RMSECV (0.461) 比 1 個 (0.317) 差，擬合 R2 卻更高？
# 2. 原始光譜 (2.547)、縮放 (2.458) 與單一波段 2258 nm (2.182) 都遠不如二階導數 + 不縮放 (0.317)。各自的原因是什麼？
# 3. 4 個樣品只有 3 個自由度（置中後），為什麼 3 個成分一定 R2 = 1？

# =====================================================================
# 4. 洗牌檢驗：把濃度亂排，PLS 還有沒有預測力？
# =====================================================================
hr("4. 洗牌檢驗：24 種濃度排列")
perm_all <- function(v) {                                # 全部排列（4! = 24），字典序，第 1 個是原本的順序
  if (length(v) == 1) return(list(v))
  out <- list()
  for (i in seq_along(v)) for (p in perm_all(v[-i])) out[[length(out) + 1]] <- c(v[i], p)
  out
}
perms <- perm_all(1:4)
sh <- data.frame(perm = sapply(perms, paste, collapse = ""), true = sapply(perms, identical, 1:4),
                 rmsecv = sapply(perms, function(p) pls_loocv(Db, conc[p], 1, FALSE)$rmsecv))
wrong <- !sh$true
n_wrong_better <- sum(sh$rmsecv[wrong] <= cv1$rmsecv)
srt <- sort(sh$rmsecv)
print(head(sh[order(sh$rmsecv), ], 5), digits = 4, row.names = FALSE)
cat(sprintf("真正順序的 RMSECV = %.3f，在 24 種裡排第 %d（最低）；錯誤排列中 RMSECV <= %.3f 的有 %d 種；第二低 = %.3f\n",
            sh$rmsecv[sh$true], rank(sh$rmsecv)[sh$true], cv1$rmsecv, n_wrong_better, srt[2]))
cat("對比第 3 章：挑波長時亂排的濃度仍能湊出高 R2；這裡波段與成分數事先固定，亂排就湊不出好的預測。\n")

png(file.path(FIG_DIR, "ch04_shuffle.png"), width = 900, height = 420, res = 110)
par(mar = c(4, 4, 2, 1))
o <- order(sh$rmsecv)
barplot(sh$rmsecv[o], col = ifelse(sh$true[o], "red", "gray70"), ylab = "RMSECV (% acid)", xlab = "24 orderings (sorted; red = true order)",
        main = "PLS-1 RMSECV for all 24 concentration orderings")
dev.off()

chk("shuffle: number of orderings", nrow(sh), 24, 0)
chk("shuffle: true ordering is the lowest (1 = yes)", as.numeric(which.min(sh$rmsecv) == which(sh$true)), 1, 0)
chk("shuffle: true ordering RMSECV", sh$rmsecv[sh$true], 0.317, 0.002)
chk("shuffle: wrong orderings with RMSECV <= 0.317", n_wrong_better, 0, 0)
chk("shuffle: second lowest RMSECV", srt[2], 0.921, 0.002)
show_checks("第 4 節")

# ---- 動手改改看（第 4 節）----
# 1. 把 N_COMP 的比較改成 2 個成分，真正順序還是最低嗎？
# 2. 這個檢驗的邏輯與第 3 章洗牌檢驗相同，結論有什麼不同？為什麼？（提示：本章沒有「從 73 個波長挑最大」這一步）

# =====================================================================
# 5. PLS-1 迴歸係數：模型看了哪些波長？
# =====================================================================
hr("5. 迴歸係數（全部 4 點、1 個成分）")
fit_all <- pls1_fit(Db, conc, 1)
coef <- fit_all$B
top <- order(-abs(coef))[1:5]
cat("|係數| 最大的 5 個波長（nm）：", sprintf("%.1f", wlb[top]), "\n")
share <- function(lo, hi) { m <- wlb >= lo & wlb <= hi; sum(abs(coef[m])) / sum(abs(coef)) }
sh_a <- share(1840, 1872); sh_b <- share(2200, 2285)
cat(sprintf("|係數| 占比：1840-1872 nm = %.3f（水的 1900 nm 帶的短波側）；2200-2285 nm = %.3f（含 2258 nm 的醋酸吸收帶）\n", sh_a, sh_b))
cat("模型自己找到了與化學相符的區域（不是我們告訴它的），而且用了很多波長一起平均雜訊。\n")

png(file.path(FIG_DIR, "ch04_coef.png"), width = 1000, height = 420, res = 110)
par(mar = c(4, 4, 2, 1))
plot(wlb, coef, type = "h", lwd = 2, col = "steelblue", xlab = "wavelength (nm)", ylab = "PLS-1 regression coefficient",
     main = "PLS-1 coefficients (2nd derivative, 73 wavelengths)")
abline(h = 0)
points(wlb[top], coef[top], pch = 17, col = "red", cex = 1.2)
legend("topleft", c("coefficient", "top 5 |coef|"), lty = c(1, NA), pch = c(NA, 17), col = c("steelblue", "red"), cex = 0.8)
dev.off()

key_top <- c(1862.4, 1871.0, 1853.9, 2258.1, 2245.5)
chk("top-5 |coef| wavelengths match (of 5)", sum(abs(wlb[top] - key_top) < 0.05), 5, 0)
chk("top-1 wavelength (nm)", wlb[top[1]], 1862.4, 0.05)
chk("top-2 wavelength (nm)", wlb[top[2]], 1871.0, 0.05)
chk("top-3 wavelength (nm)", wlb[top[3]], 1853.9, 0.05)
chk("share of |coef| in 1840-1872 nm", sh_a, 0.391, 0.002)
chk("share of |coef| in 2200-2285 nm", sh_b, 0.276, 0.002)
show_checks("第 5 節")

# ---- 動手改改看（第 5 節）----
# 1. 係數有正有負、形狀像二階導數的吸收峰。為什麼？（提示：二階導數把吸收峰變成「負的尖峰」）
# 2. 只留 |係數| 最大的 5 個波長重做 LOOCV（這在迴圈外挑波長！）RMSECV 會變多少？這犯了第 3 章的哪個錯？

# =====================================================================
# 6. 學生練習：換設定、比較
# =====================================================================
hr("6. 學生練習：scale、原始光譜、成分數")
cat("本章的設定（二階導數、兩段波段、不縮放、1 個成分）是事先固定的。下面請你自己換設定，看 RMSECV 怎麼變。\n")
cat("Orange 版的操作說明在 exercises/orange/ch04_orange_guide.md。\n\n")
try_variant <- function(name, M, nc, scale) {
  cat(sprintf("%-34s nc = %d, scale = %-5s -> RMSECV = %.3f %%\n", name, nc, scale, pls_loocv(M, conc, nc, scale)$rmsecv))
}
# 請改下面的設定再執行：
try_variant("baseline (d2, no scale, 1 comp)", Db, 1, FALSE)
try_variant("scale = TRUE", Db, 1, TRUE)
try_variant("raw spectra (no derivative)", Rb, 1, FALSE)
try_variant("2 components", Db, 2, FALSE)
try_variant("raw + scale = TRUE", Rb, 1, TRUE)
try_variant("scale = TRUE, 2 components", Db, 2, TRUE)

variants <- c("d2 + 1 comp\n(fixed)" = cv1$rmsecv, "d2 + 2 comp" = cv2$rmsecv, "d2 + scale" = cvs$rmsecv,
              "raw + 1 comp" = cvr$rmsecv, "single band\n2258 nm" = RM2258)
png(file.path(FIG_DIR, "ch04_compare.png"), width = 900, height = 450, res = 110)
par(mar = c(5, 4, 2, 1))
barplot(variants, col = c("red", "steelblue", "steelblue", "steelblue", "gray70"), ylab = "RMSECV (% acid)",
        main = "RMSECV of different settings (n = 4, LOOCV)", cex.names = 0.75)
dev.off()

# ---- 動手改改看（第 6 節）----
# 1. scale = TRUE 為什麼變差？（提示：縮放把每個波長的變異都放成 1，雜訊多的波長被放大到與訊號一樣重要）
# 2. 為什麼不能把 nc 設成 3 再做 LOOCV？（提示：每折只剩 3 個訓練樣品，置中後只剩 2 個自由度；試著執行看看會怎樣）
# 3. 在這 4 個樣品上，你『試了很多組設定、挑 RMSECV 最低的』報出來，會有什麼問題？（回想第 3 章）
# 4. 到 Orange 照 exercises/orange/ch04_orange_guide.md 的步驟重做，結果應與本頁一致（PLS 小工具的預設是縮放，會是 2.458）。

# =====================================================================
# 總結
# =====================================================================
hr("總結")
cat(sprintf("全部自我驗證：%d / %d PASS\n", sum(checks$ok), nrow(checks)))
if (!all(checks$ok)) { cat("FAIL 項目：\n"); print(checks[!checks$ok, ], row.names = FALSE) }
cat("1. PCA 描述資料：原始光譜的 PC1 = 光徑（0.917），二階導數的 PC1 = 濃度（0.969）。\n")
cat("2. PLS 1 個成分、二階導數、兩段波段、不縮放：RMSECV = 0.317 %；設定都是事先固定的。\n")
cat("3. 洗牌檢驗：24 種排列只有真正順序最低，沒有錯誤排列比它好。\n")
cat("4. 係數集中在水與醋酸的吸收帶，與化學相符。\n")
cat("5. 只有 4 個樣品，0.317 % 仍是『初步』結果：需要更多獨立樣品驗證。\n")
cat("圖已存到", normalizePath(FIG_DIR), "\n")
