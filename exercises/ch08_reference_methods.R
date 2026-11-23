# =====================================================================
# 第 8 章 參考方法的誤差：滴定 vs Atago，R 程式練習（只用 base R）
# 前面幾章都把「參考值」當成真值。這一章反過來問：參考值本身有多準？
#   四組參考值（標示、Atago 兩天、滴定）怎麼算、誤差預算有多大、彼此差多少、
#   NIR 能不能當旁證、NaOH 力價的比例誤差、滴定量怎麼設計、各種誤差大小比一比，
#   最後產生一份「滴定紀錄表」ch08_titration_sheet.csv。
# 所有計算都是確定性的（沒有任何亂數）。
#
# 執行方式：Rscript ch08_reference_methods.R （程式會自動切到本檔所在資料夾）
# 資料：../data/spectra_long.csv（只有第 4 節、第 7 節用到光譜）  圖：figures_R/ch08_*.png（圖上文字一律用英文）
# 輸出：ch08_titration_sheet.csv（滴定紀錄表，UTF-8 含 BOM，Excel 可直接開）；列印版 SOP 見 ch08_titration_sop.md
# 每節結尾有「動手改改看」，每節最後有該節的自我驗證 PASS/FAIL。
#
# 事先固定的設定：
#   滴定資料（另一組學生，11/03）：取 10 mL 醋稀釋到 100 mL，再取 25 mL，以 NaOH 滴定；
#   NaOH 由稱重配製（0.4315 g / 40 / 0.1 L -> 力價 f = 1.07875，未標定）。
#   醋酸當量 0.006 g/mL（0.1 N x 60.05 g/mol / 1000）；結果單位是 w/v %。
#   器材允差採 A 級玻璃器皿的典型值（示意），標準不確定度 = 允差 / sqrt(3)（矩形分布）。
#   NIR：1 mm、重複 1 的 4 瓶，二階導數 Savitzky-Golay 9 點、二次多項式，取 2258 nm，乘 1e3
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
# 0. 載入資料與二階導數（第 4、7 節會用到）
# =====================================================================
hr("0. 載入資料與二階導數")
df <- read.csv(DATA_PATH)
r1 <- df[df$replicate == 1 & abs(df$path_mm - 1) < 1e-9, ]
wl <- sort(unique(r1$wavelength_nm))
ACID4 <- c(0, 3.32, 4.87, 7.45)                          # 4 瓶的 Atago 11/02 讀值（資料檔裡的 acid_pct）
X1 <- t(sapply(ACID4, function(c) {                      # 4 x n_wl：1 mm、重複 1 的 4 瓶原始光譜
  h <- r1[abs(r1$acid_pct - c) < 1e-9, ]; h <- h[order(h$wavelength_nm), ]
  stopifnot(nrow(h) == length(wl)); h$absorbance
}))
sg_coeffs2 <- function(window) {
  hh <- (window - 1) / 2; xx <- -hh:hh
  M <- cbind(1, xx, xx^2)
  as.numeric(2 * (solve(t(M) %*% M) %*% t(M))[3, ])
}
sg_d2 <- function(S, window = SG_WIN) {                  # 二階導數；2200-2300 nm 離邊緣很遠，邊緣怎麼處理不影響結果
  cf <- sg_coeffs2(window); half <- (window - 1) / 2
  out <- matrix(NA_real_, nrow(S), ncol(S))
  for (i in (half + 1):(ncol(S) - half)) out[, i] <- S[, (i - half):(i + half), drop = FALSE] %*% cf
  out[, 1:half] <- out[, half + 1]
  out[, (ncol(S) - half + 1):ncol(S)] <- out[, ncol(S) - half]
  out
}
D1 <- sg_d2(X1)
j <- which.min(abs(wl - WL0))
sig <- D1[, j] * 1e3                                      # 1 mm 4 瓶的 2258 nm 訊號 y（x 10^-3）
cat(sprintf("1 mm、重複 1：%d 瓶 x %d 個波長；2258 nm 取第 %d 點（%.1f nm）\n", nrow(X1), length(wl), j, wl[j]))
cat("訊號 y（x1e-3）：", round(sig, 4), "\n")
fit_line <- function(xs, ys) {
  n <- length(xs); sxx <- sum((xs - mean(xs))^2)
  b <- sum((xs - mean(xs)) * (ys - mean(ys))) / sxx; a <- mean(ys) - b * mean(xs)
  s <- sqrt(sum((ys - a - b * xs)^2) / (n - 2))
  list(a = a, b = b, s = s, s_conc = s / abs(b), sxx = sxx)
}

# ---- 動手改改看（第 0 節）----
# 1. 這一章的第 1-3、5、6 節完全不需要光譜，只用實驗常數。為什麼參考方法的誤差可以不靠 NIR 就估出來？
# 2. 為什麼第 4 節只用 3 個點、而且自由度只有 1？

# =====================================================================
# 1. 四組參考值
# =====================================================================
hr("1. 四組參考值：標示、Atago 11/02、Atago 11/03、滴定")
cat("Atago 讀值來自 ATAGO PAL-Easy ACID：這是以電導度為基礎的酸度計，不是折射計；結果是 w/v %。\n")
cat("滴定：酸度 (w/v %) = V x f x 0.006 x (100/25) x (100/10)\n")
cat("  V = NaOH 用量 (mL)；f = NaOH 力價（以稱重配製、未標定）；0.006 g/mL = 每 mL 0.1 N NaOH 對應的醋酸克數；\n")
cat("  100/25 = 100 mL 瓶裡只取 25 mL 滴定；100/10 = 10 mL 醋稀釋到 100 mL。\n")
V <- c(8.8, 17, 44)                                       # NaOH 用量（mL）
f <- 0.4315 / 40 / 0.1 * 10                               # 1.07875
K <- f * 0.006 * (100 / 25) * (100 / 10)                  # 每 mL NaOH 對應的酸度 %
titr <- V * K
label <- c(3, 5, 10); at1 <- c(3.32, 4.87, 7.45); at2 <- c(3.01, 4.12, 6.86)
ratio <- titr / at2
daydiff <- (at1 - at2) / at1 * 100
cat(sprintf("f = 0.4315 / 40 / 0.1 * 10 = %.5f；K = f x 0.006 x 4 x 10 = %.5f %%/mL\n", f, K))
cat(sprintf("%-6s %7s %10s %10s %9s %9s %11s %12s\n", "bottle", "label", "Atago11/02", "Atago11/03", "V (mL)", "titration", "titr/Atago03", "day diff (%)"))
for (i in 1:3) cat(sprintf("%-6d %7.2f %10.2f %10.2f %9.1f %9.3f %11.4f %12.3f\n", i, label[i], at1[i], at2[i], V[i], titr[i], ratio[i], daydiff[i]))
eq_exact <- 0.1 * 60.05 / 1000
cat(sprintf("講義用 0.006 g/mL；精確值 0.1 N x 60.05 / 1000 = %.6f g/mL，差 %.2f %%（可忽略）。\n", eq_exact, (eq_exact - 0.006) / eq_exact * 100))
cat("同一批醋，隔一天用同一支 Atago 量，讀值就差 8-15 %；而滴定和標示、Atago 三者又彼此對不上。第 3 節會量化。\n")

chk("f (NaOH factor, unstandardized)", f, 1.079)
chk("K (% per mL NaOH)", K, 0.2589)
for (i in 1:3) chk(sprintf("titration acidity, bottle %d (V=%.1f mL)", i, V[i]), titr[i], c(2.278, 4.401, 11.39)[i])
for (i in 1:3) chk(sprintf("titration / Atago 11/03, bottle %d", i), ratio[i], c(0.7569, 1.068, 1.661)[i])
for (i in 1:3) chk(sprintf("Atago day difference (%%), bottle %d", i), daydiff[i], c(9.337, 15.4, 7.919)[i])
chk("0.1 N x 60.05 / 1000 (g/mL)", eq_exact, 0.006005)
chk("aliquot factor 100/25", 100 / 25, 4)
chk("dilution factor 100/10", 100 / 10, 10)
show_checks("第 1 節")

# ---- 動手改改看（第 1 節）----
# 1. 為什麼 3 個瓶子的 titr/Atago 比值是 0.76、1.07、1.66，而不是同一個常數？（如果只是力價不準，比值應該怎樣？）
# 2. Atago 兩天的差 8-15 %，是儀器不穩，還是樣品變了？要怎麼設計實驗分辨？
# 3. 把 0.006 改成 0.006005，所有滴定酸度改變多少 %？

# =====================================================================
# 2. 誤差預算
# =====================================================================
hr("2. 滴定的誤差預算（單次滴定，相對標準不確定度）")
cat("器材允差採 A 級玻璃器皿的典型值（示意，請以你手上器材標示為準）；標準不確定度 = 允差 / sqrt(3)（矩形分布）。\n")
u_read <- 0.02                                            # 每次讀滴定管 +/- 0.02 mL
tol <- c(buret = 0.05, endpoint = 0.05, pip10 = 0.02, flask100 = 0.08, pip25 = 0.03)
u <- tol / sqrt(3)
uV <- sqrt(2 * u_read^2 + u[["buret"]]^2 + u[["endpoint"]]^2)   # 兩次讀數 + 滴定管允差 + 終點一滴
cat(sprintf("滴定體積不確定度 uV = sqrt(2 x %.2f^2 + u_buret^2 + u_endpoint^2) = %.5f mL（兩次讀數各 %.2f mL；滴定管允差 %.2f、終點一滴 %.2f mL）\n",
            u_read, uV, u_read, tol[["buret"]], tol[["endpoint"]]))
rel_titre <- uV / V * 100
rel_pip10 <- rep(u[["pip10"]] / 10 * 100, 3); rel_flask <- rep(u[["flask100"]] / 100 * 100, 3); rel_pip25 <- rep(u[["pip25"]] / 25 * 100, 3)
total_rel <- sqrt(rel_titre^2 + rel_pip10^2 + rel_flask^2 + rel_pip25^2)
total_abs <- titr * total_rel / 100
cat(sprintf("%-8s %8s %8s %9s %9s %9s %10s %10s\n", "V (mL)", "titre %", "pip10 %", "flask100 %", "pip25 %", "total %", "acidity %", "total abs"))
for (i in 1:3) cat(sprintf("%-8.1f %8.4f %8.4f %9.4f %9.4f %9.4f %10.3f %10.5f\n", V[i], rel_titre[i], rel_pip10[i], rel_flask[i], rel_pip25[i], total_rel[i], titr[i], total_abs[i]))
cat("結論：器材造成的不確定度只有 0.2-0.6 %（絕對值約 0.013-0.021 % 酸度），體積小時由滴定管讀數主導，體積大時由 10 mL 吸量管主導。\n")
cat("這是『下限』：它沒有算進力價。NaOH 由稱重配製、沒有標定；NaOH 會吸濕、吸 CO2，真正的力價很可能比 1.07875 低 2-5 %。\n")
cat("這是另一種誤差：它是『比例偏差』（所有瓶子同一比例），不是隨機誤差；算出的酸度會同比例偏高。以 KHP 標定、三重複後，力價的相對不確定度約 0.2 %（典型值，示意）。\n")

png(file.path(FIG_DIR, "ch08_budget_bars.png"), width = 950, height = 520, res = 110)
par(mar = c(4, 4, 2, 1))
bm <- rbind(titre = rel_titre, pip10 = rel_pip10, flask100 = rel_flask, pip25 = rel_pip25, total = total_rel)
bp <- barplot(bm, beside = TRUE, names.arg = sprintf("V = %.1f mL", V), col = c("steelblue", "grey70", "grey50", "grey35", "firebrick"),
              ylab = "relative standard uncertainty (%)", ylim = c(0, 0.7), main = "Titration error budget (Class A tolerances, illustrative)")
legend("topright", rownames(bm), fill = c("steelblue", "grey70", "grey50", "grey35", "firebrick"), cex = 0.8)
dev.off()

chk("uV: titre volume uncertainty (mL)", uV, 0.04967)
for (i in 1:3) chk(sprintf("rel. u titre (%%), V=%.1f", V[i]), rel_titre[i], c(0.5644, 0.2922, 0.1129)[i])
for (i in 1:3) chk(sprintf("rel. u 10 mL pipette (%%), V=%.1f", V[i]), rel_pip10[i], c(0.1155, 0.1155, 0.1155)[i])
for (i in 1:3) chk(sprintf("rel. u 100 mL flask (%%), V=%.1f", V[i]), rel_flask[i], c(0.04619, 0.04619, 0.04619)[i])
for (i in 1:3) chk(sprintf("rel. u 25 mL pipette (%%), V=%.1f", V[i]), rel_pip25[i], c(0.06928, 0.06928, 0.06928)[i])
for (i in 1:3) chk(sprintf("total relative uncertainty (%%), V=%.1f", V[i]), total_rel[i], c(0.5821, 0.325, 0.1817)[i])
for (i in 1:3) chk(sprintf("total absolute uncertainty (%% acid), V=%.1f", V[i]), total_abs[i], c(0.01326, 0.0143, 0.0207)[i])
show_checks("第 2 節")

# ---- 動手改改看（第 2 節）----
# 1. 哪一項在 V = 8.8 mL 時最大？在 V = 44 mL 時呢？為什麼？
# 2. 如果真正的 NaOH 力價比 1.07875 低 3 %，3 個瓶子算出的酸度各偏多少（% 酸度）？這比器材誤差大幾倍？
# 3. 把 0.006 g/mL 的 0.08 % 差異也當成誤差，總不確定度會明顯變大嗎？

# =====================================================================
# 3. 參考值彼此差多少
# =====================================================================
hr("3. 參考值彼此差多少：換成滴定體積（mL）與讀數不確定度 uV 比")
cat("做法：如果某個參考值是真的，滴定時應該用掉 need = 參考值 / K mL 的 NaOH；實際用量是 V。gap = need - V；gap/uV 是差了幾個讀數不確定度。\n")
refs <- list(label = label, atago_1103 = at2, atago_1102 = at1)
need <- lapply(refs, function(r) r / K)
gap <- lapply(need, function(n) n - V)
gap_u <- lapply(gap, function(g) g / uV)
for (nm in names(refs)) {
  cat(sprintf("[%s] need (mL)：%s\n", nm, paste(sprintf("%7.2f", need[[nm]]), collapse = " ")))
  cat(sprintf("%*s gap (mL)：%s；gap/uV：%s\n", nchar(nm) + 2, "", paste(sprintf("%7.3f", gap[[nm]]), collapse = " "), paste(sprintf("%8.2f", gap_u[[nm]]), collapse = " ")))
}
all_gu <- abs(unlist(gap_u))
cat(sprintf("9 個差距的大小介於 %.0f 到 %.0f 倍 uV。\n", min(all_gu), max(all_gu)))
cat("結論：差距是讀數不確定度的 22-352 倍，不可能是讀滴定管的誤差。\n")
cat("第 1 瓶（標示 3 %）的滴定值（2.28 %）低於標示、Atago 兩天的讀值，三個都比它高：它是最可疑的一瓶，應該重測；但這不等於『證明滴定錯了』。\n")

png(file.path(FIG_DIR, "ch08_gaps.png"), width = 950, height = 520, res = 110)
par(mar = c(4, 4, 2, 1))
gm <- rbind(label = gap_u$label, "Atago 11/03" = gap_u$atago_1103, "Atago 11/02" = gap_u$atago_1102)
barplot(gm, beside = TRUE, names.arg = sprintf("bottle %d (V = %.1f mL)", 1:3, V), col = c("grey60", "steelblue", "darkorange"),
        ylab = "gap / uV  (titre needed - titre observed)", main = "How far each reference is from the titration, in units of reading uncertainty")
abline(h = 0)
legend("bottomleft", rownames(gm), fill = c("grey60", "steelblue", "darkorange"), cex = 0.8)
dev.off()

key_need <- list(label = c(11.59, 19.31, 38.62), atago_1103 = c(11.63, 15.91, 26.5), atago_1102 = c(12.82, 18.81, 28.78))
key_gap <- list(label = c(2.787, 2.312, -5.375), atago_1103 = c(2.826, -1.087, -17.5), atago_1102 = c(4.023, 1.81, -15.22))
key_gu <- list(label = c(56.13, 46.56, -108.2), atago_1103 = c(56.9, -21.88, -352.4), atago_1102 = c(81.01, 36.45, -306.5))
for (nm in names(refs)) {
  for (i in 1:3) chk(sprintf("%s: need (mL), bottle %d", nm, i), need[[nm]][i], key_need[[nm]][i])
  for (i in 1:3) chk(sprintf("%s: gap (mL), bottle %d", nm, i), gap[[nm]][i], key_gap[[nm]][i])
  for (i in 1:3) chk(sprintf("%s: gap / uV, bottle %d", nm, i), gap_u[[nm]][i], key_gu[[nm]][i])
}
chk("smallest |gap| / uV", min(all_gu), 21.88)
chk("largest |gap| / uV", max(all_gu), 352.4)
chk("bottle 1: references above the titration (of 3)", sum(sapply(gap, function(g) g[1] > 0)), 3)
show_checks("第 3 節")

# ---- 動手改改看（第 3 節）----
# 1. 如果只看 gap_in_uV，哪一個瓶子、哪一個參考的差最大？這能說明是哪一種方法錯嗎？
# 2. 第 3 瓶（V = 44 mL）三個 gap 都是負的。這個瓶子的滴定值（11.39 %）比標示 10 % 還高，可能的原因有哪些？
# 3. 為什麼說『最可疑』不等於『證明錯了』？還缺什麼證據？

# =====================================================================
# 4. NIR 旁證很弱
# =====================================================================
hr("4. NIR 旁證：3 點直線預測第 1 瓶（1 mm、2258 nm 二階導數 x1e3）")
ref4 <- c(0, titr)                                        # 水 + 3 瓶的滴定值
keep <- c(1, 3, 4)                                        # 水、第 2 瓶、第 3 瓶；留下第 1 瓶當未知
xs <- ref4[keep]; ys <- sig[keep]; n3 <- 3
f3 <- fit_line(xs, ys); b3 <- f3$b; a3 <- f3$a; s3 <- f3$s
sxx3 <- sum((xs - mean(xs))^2); tq <- qt(.975, n3 - 2); kk <- tq * tq * s3 * s3; y0 <- sig[2]
g3 <- kk / (b3 * b3 * sxx3)
A <- b3 * b3 - kk / sxx3; B <- -2 * b3 * (y0 - a3) + 2 * kk * mean(xs) / sxx3
C <- (y0 - a3)^2 - kk * (1 + 1 / n3) - kk * mean(xs)^2 / sxx3
pred3 <- (y0 - a3) / b3
lo3 <- hi3 <- NA_real_
if (A > 0) { r <- sqrt(B * B - 4 * A * C); rt <- sort(c((-B - r) / (2 * A), (-B + r) / (2 * A))); lo3 <- rt[1]; hi3 <- rt[2] }
cat(sprintf("3 點：水 + 第 2、3 瓶（用滴定值當 x）；殘差自由度 = 1，t(0.975, 1) = %.3f\n", tq))
cat(sprintf("預測第 1 瓶：%.3f %%（滴定值是 %.3f %%）；g = %.4f（<1，Fieller 區間有界）\n", pred3, titr[1], g3))
cat(sprintf("Fieller 95 %% 區間：%.2f 到 %.2f %%\n", lo3, hi3))
cat("結論：區間從 -14 到 13 %，什麼都證明不了；它同時包含『滴定對』與『標示對』。\n")
cat("而且這是循環論證：用 NIR 來評斷參考方法，可是 NIR 的校正又是靠這些參考方法建立的。NIR 不能當裁判。\n")

chk("NIR prediction of bottle 1 (%)", pred3, 3.091)
chk("g statistic", g3, 0.6414)
chk("Fieller lower limit (%)", lo3, -14.25)
chk("Fieller upper limit (%)", hi3, 12.66)
show_checks("第 4 節")

# ---- 動手改改看（第 4 節）----
# 1. 預測值 3.09 % 剛好接近標示 3 %、遠離滴定 2.28 %。這算不算『NIR 支持標示』？看區間再回答。
# 2. 自由度 1 的 t(0.975) 是多少？如果有 8 個校正點（自由度 6），區間會縮小到多少量級？
# 3. 為什麼用 NIR 評斷參考方法是循環論證？

# =====================================================================
# 5. 力價的比例誤差
# =====================================================================
hr("5. 力價的比例誤差：參考值全部乘 (1 + e)")
e <- -0.03
r2f <- function(x, y) cor(x, y)^2
fa0 <- fit_line(ref4, sig); fa1 <- fit_line(ref4 * (1 + e), sig)
y0f <- 4.20
r2_0 <- r2f(ref4, sig); r2_1 <- r2f(ref4 * (1 + e), sig)
pred0 <- (y0f - fa0$a) / fa0$b; pred1 <- (y0f - fa1$a) / fa1$b
khp_g <- 0.1 * 0.020 * 204.22                             # 20 mL 的 0.1 N NaOH 需要的 KHP 克數
cat(sprintf("e = %.2f：R^2 %.4f -> %.4f（完全不變）；斜率 %.4f -> %.4f；訊號 y0 = %.2f 的預測 %.3f -> %.3f %%\n",
            e, r2_0, r2_1, fa0$b, fa1$b, y0f, pred0, pred1))
cat("結論：比例誤差不會讓 R^2 變差，所以光看 R^2 完全看不出力價錯了；它只會讓斜率與預測同比例偏（預測 = 原來 x (1 + e)）。\n")
cat("解法：用 KHP（鄰苯二甲酸氫鉀，M = 204.22 g/mol，一級標準）標定 NaOH，f = m / (204.22 x 0.1 x V/1000)。\n")
cat(sprintf("標定 20 mL 的 0.1 N NaOH 需要 KHP：0.1 x 0.020 x 204.22 = %.4f g。\n", khp_g))

chk("R^2 (original)", r2_0, 0.9898)
chk("R^2 (references x 0.97)", r2_1, 0.9898)
chk("slope (original)", fa0$b, -0.6415)
chk("slope (references x 0.97)", fa1$b, -0.6614)
chk("prediction at y0 = 4.20 (original, %)", pred0, 6.131)
chk("prediction at y0 = 4.20 (x 0.97, %)", pred1, 5.947)
chk("prediction ratio scaled / original", pred1 / pred0, 0.97)
chk("KHP mass for 20 mL of 0.1 N NaOH (g)", khp_g, 0.4084)
show_checks("第 5 節")

# ---- 動手改改看（第 5 節）----
# 1. 為什麼 R^2 不變？（提示：R^2 只看點排成直線的程度，與 x 軸整體伸縮無關）
# 2. 如果真正的力價低 5 %，預測會偏多少 %？
# 3. 若只有 R^2 能看，要怎麼發現力價誤差？（提示：需要獨立的標準，例如 KHP）

# =====================================================================
# 6. 滴定量的設計
# =====================================================================
hr("6. 滴定體積的選擇：D mL 醋稀釋到 100 mL，取 25 mL 滴定")
cat("滴定體積 V = 酸度 / (f x 0.006 x 4 x (100 / D))；目標 15-40 mL（50 mL 滴定管）；太小讀數誤差佔比大，超過 50 mL 要重新裝滿滴定管。\n")
concs <- c(3, 5, 7, 11.39); dil <- c(5, 10, 20)
titre_fn <- function(c, D) c / (f * 0.006 * (100 / 25) * (100 / D))
Vm <- outer(concs, dil, Vectorize(titre_fn))
Rm <- uV / Vm * 100
inrange <- Vm >= 15 & Vm <= 40
over50 <- Vm > 50
cat(sprintf("%-10s %s\n", "acidity %", paste(sprintf("D = %2d mL", dil), collapse = "      ")))
for (i in seq_along(concs)) {
  cat(sprintf("%-10.2f", concs[i]))
  for (jj in seq_along(dil)) cat(sprintf(" %7.2f mL %s", Vm[i, jj], if (over50[i, jj]) "(>50)" else if (inrange[i, jj]) "(OK) " else "     "))
  cat("\n")
}
cat(sprintf("15-40 mL 的格子有 %d 個；超過 50 mL 的有 %d 個。\n", sum(inrange), sum(over50)))
cat("相對讀數不確定度 uV / V（%）：體積越大越小，D = 10 mL、V 約 20 mL 時約 0.25 %。\n")
cat("D 越大，滴定量越大。D = 10 mL 對 5-7 % 的醋很合適（19-27 mL）；3 % 的醋 D = 10 mL 只有 11.6 mL（太小），改用 D = 20 mL；11.39 % 的醋 D = 10 mL 要 44 mL，改用 D = 5 mL 才有 22 mL。\n")

png(file.path(FIG_DIR, "ch08_aliquot_heatmap.png"), width = 900, height = 560, res = 110)
par(mar = c(7, 5, 3, 1))
cls <- ifelse(Vm > 50, 4, ifelse(Vm > 40, 3, ifelse(Vm >= 15, 2, 1)))
z <- matrix(0, length(dil), length(concs))
for (i in seq_along(concs)) for (jj in seq_along(dil)) z[jj, length(concs) + 1 - i] <- cls[i, jj]
image(seq_along(dil), seq_along(concs), z, breaks = 0.5:4.5, col = c("grey85", "palegreen3", "khaki", "firebrick"), axes = FALSE,
      xlab = "vinegar volume D (mL) diluted to 100 mL", ylab = "acidity (%)", main = "Titre volume V (mL): green = 15-40 mL target")
axis(1, at = seq_along(dil), labels = dil); axis(2, at = seq_along(concs), labels = rev(concs), las = 1); box()
for (i in seq_along(concs)) for (jj in seq_along(dil)) text(jj, length(concs) + 1 - i, sprintf("%.1f", Vm[i, jj]), cex = 1.1)
legend("bottom", c("< 15", "15-40 (target)", "40-50", "> 50 (refill)"), fill = c("grey85", "palegreen3", "khaki", "firebrick"), cex = 0.75, horiz = TRUE, inset = c(0, -0.42), xpd = TRUE)
dev.off()

key_V <- rbind(c(5.794, 11.59, 23.17), c(9.656, 19.31, 38.62), c(13.52, 27.04, 54.07), c(22.0, 43.99, 87.99))
key_R <- rbind(c(0.8572, 0.4286, 0.2143), c(0.5143, 0.2572, 0.1286), c(0.3674, 0.1837, 0.09185), c(0.2258, 0.1129, 0.05645))
for (i in seq_along(concs)) for (jj in seq_along(dil)) chk(sprintf("titre V (mL): %.2f %%, D=%d", concs[i], dil[jj]), Vm[i, jj], key_V[i, jj])
for (i in seq_along(concs)) for (jj in seq_along(dil)) chk(sprintf("rel. reading u (%%): %.2f %%, D=%d", concs[i], dil[jj]), Rm[i, jj], key_R[i, jj])
chk("cells within 15-40 mL (of 12)", sum(inrange), 5)
chk("cells above 50 mL (of 12)", sum(over50), 2)
show_checks("第 6 節")

# ---- 動手改改看（第 6 節）----
# 1. 11.39 % 的醋用 D = 20 mL 會得到幾 mL？為什麼不好？（提示：滴定管只有 50 mL）
# 2. 要滴 3 % 的醋，D 取多少最合適？（提示：目標 15-40 mL）
# 3. 滴定量 V 加倍，讀數誤差佔比減半；那為什麼不一直把 V 做得很大？

# =====================================================================
# 7. 誤差比一比
# =====================================================================
hr("7. 誤差比一比（單位都是 % 酸度）")
res <- sig - (fa0$a + fa0$b * ref4)
nir_resid_titr <- sqrt(sum(res^2) / 2) / abs(fa0$b)       # 4 點、自由度 2
titr_floor_abs <- total_abs
atago_day_sd_abs <- sqrt(mean((at1 - at2)^2) / 2)
titr_3pct_gap_abs <- 3.0 - titr[1]
triplicate_factor <- 1 / sqrt(3)
cat(sprintf("NIR 校正殘差（4 點，滴定當參考，自由度 2）  ：%.4f %%   [實測]\n", nir_resid_titr))
cat(sprintf("Atago 兩天的差異 SD（3 瓶、各 2 次）        ：%.4f %%   [實測]\n", atago_day_sd_abs))
cat(sprintf("滴定 vs 標示：第 1 瓶（3 %% 標示）差        ：%.4f %%   [實測的差，但標示本身是否準不知道]\n", titr_3pct_gap_abs))
cat(sprintf("滴定器材下限（3 瓶）                        ：%s %%   [假設：典型允差，不是量到的]\n", paste(sprintf("%.4f", titr_floor_abs), collapse = " / ")))
cat(sprintf("三重複可把隨機誤差縮小為 1/sqrt(3) = %.4f 倍\n", triplicate_factor))
cat("閱讀：器材下限（約 0.02 %）比其他幾項小 20-50 倍，但它是假設值；Atago 兩天差 0.41 % 和 NIR 殘差 0.61 % 是實測，量級一樣。\n")
cat("所以 NIR 校正殘差裡，很大一部分可能是參考值自己的誤差，而不是光譜不好。\n")

png(file.path(FIG_DIR, "ch08_error_compare.png"), width = 1000, height = 520, res = 110)
par(mar = c(4, 17, 2, 1))
ev <- c(titr_floor_abs[3], titr_floor_abs[2], titr_floor_abs[1], titr_3pct_gap_abs, atago_day_sd_abs, nir_resid_titr)
en <- c("titration floor, V=44 mL (assumed)", "titration floor, V=17 mL (assumed)", "titration floor, V=8.8 mL (assumed)",
        "titration vs label, bottle 1 (measured)", "Atago day-to-day SD (measured)", "NIR residual, 4 pts (measured)")
ec <- c(rep("grey55", 3), rep("steelblue", 3))
plot(ev, seq_along(ev), log = "x", xlim = c(0.005, 2), ylim = c(0.5, 6.8), pch = 19, col = ec, cex = 1.6, yaxt = "n", ylab = "", xlab = "error size (% acid, log scale)",
     main = "Reference-method errors vs NIR residual")
segments(0.005, seq_along(ev), ev, seq_along(ev), col = ec, lwd = 2)
axis(2, at = seq_along(ev), labels = en, las = 1, cex.axis = 0.75)
text(ev, seq_along(ev) + 0.3, sprintf("%.3f", ev), cex = 0.8)
dev.off()

chk("NIR residual, titration reference (% acid)", nir_resid_titr, 0.612)
for (i in 1:3) chk(sprintf("titration floor, V=%.1f mL (%% acid)", V[i]), titr_floor_abs[i], c(0.01326, 0.0143, 0.0207)[i])
chk("Atago day-to-day SD (% acid)", atago_day_sd_abs, 0.4096)
chk("titration vs label 3 % (% acid)", titr_3pct_gap_abs, 0.7217)
chk("triplicate factor 1/sqrt(3)", triplicate_factor, 0.5774)
show_checks("第 7 節")

# ---- 動手改改看（第 7 節）----
# 1. 哪些數字是『實測』，哪些是『假設』？為什麼一定要分清楚？
# 2. 如果滴定做三重複，器材下限變成多少？這會不會改變結論？
# 3. 力價的比例偏差（2-5 %）折算成酸度，對第 2 瓶（4.4 %）是多少 % 酸度？它和上面哪一項同量級？

# =====================================================================
# 8. 產生滴定紀錄表
# =====================================================================
hr("8. 產生 ch08_titration_sheet.csv（滴定紀錄表）")
cols <- c("步驟", "樣品", "重複", "秤重_g", "起點_mL", "終點_mL", "用量_mL", "備註")
rows8 <- c(sprintf("KHP 標定,KHP,%d", 1:3), sprintf("空白,空白,%d", 1:2), sprintf("樣品滴定,原液,%d", 1:3),
           unlist(lapply(c("A", "B", "C", "D"), function(L) sprintf("樣品滴定,盲樣%s,%d", L, 1:3))))
txt <- paste0(paste(cols, collapse = ","), "\n", paste0(rows8, strrep(",", 5), "\n", collapse = ""))
con <- file("ch08_titration_sheet.csv", "wb")
writeBin(c(as.raw(c(0xEF, 0xBB, 0xBF)), charToRaw(enc2utf8(txt))), con); close(con)
ts <- read.csv("ch08_titration_sheet.csv", fileEncoding = "UTF-8-BOM", check.names = FALSE)
cat(sprintf("已寫出 ch08_titration_sheet.csv：%d 列，%d 欄\n", nrow(ts), ncol(ts)))
cat("順序：KHP 標定 3 次 -> 空白 2 次 -> 原液 3 次 -> 盲樣 A、B、C、D 各 3 次；秤重、起點、終點、用量、備註 5 欄留白，實驗時手填。\n")
cat("『用量』= 終點 - 起點；樣品的酸度要先扣掉空白用量再計算。流程、公式與稀釋建議見 ch08_titration_sop.md。\n")

chk("sheet: number of rows (3 + 2 + 3 + 12)", nrow(ts), 20)
chk("sheet: number of columns", ncol(ts), 8)
chk("sheet: KHP rows", sum(ts[["步驟"]] == "KHP 標定"), 3)
chk("sheet: blank rows", sum(ts[["步驟"]] == "空白"), 2)
chk("sheet: stock vinegar rows", sum(ts[["樣品"]] == "原液"), 3)
chk("sheet: blind-sample rows (A-D x 3)", sum(grepl("^盲樣", ts[["樣品"]])), 12)
chk("sheet: replicate-number mismatches", sum(ts[["重複"]] != c(1:3, 1:2, 1:3, rep(1:3, 4))), 0)
chk("sheet: numeric cells not blank (weight/start/end/titre)", sum(!is.na(ts[, 4:7])), 0)
chk("sheet: remark cells not blank", sum(!is.na(ts[["備註"]])), 0)
show_checks("第 8 節")

# ---- 動手改改看（第 8 節）----
# 1. 用 Excel 開 ch08_titration_sheet.csv，中文有亂碼嗎？（檔案含 BOM，應該不會）
# 2. 為什麼 KHP 要稱重 3 份、空白要做 2 次？
# 3. 為什麼盲樣的順序最好由別人編碼、不要依濃度排？

# =====================================================================
# 總結
# =====================================================================
hr("總結")
cat(sprintf("全部自我驗證：%d / %d PASS\n", sum(checks$ok), nrow(checks)))
if (!all(checks$ok)) { cat("FAIL 項目：\n"); print(checks[!checks$ok, ], row.names = FALSE) }
cat("1. 四組參考值（標示、Atago 兩天、滴定）互相對不上；Atago 兩天就差 8-15 %。Atago 是電導式酸度計，不是折射計。\n")
cat("2. 滴定器材的不確定度只有 0.2-0.6 %（約 0.013-0.021 % 酸度）；未標定 NaOH 的力價偏差（可能偏低 2-5 %）是另一種比例偏差。\n")
cat("3. 參考值之間的差距是讀數不確定度的 22-352 倍，不是讀數誤差；第 1 瓶（標示 3 %）的滴定值最可疑，應重測。\n")
cat("4. NIR 旁證太弱（區間 -14 到 13 %），而且循環論證，不能當裁判。\n")
cat("5. 力價的比例誤差不影響 R^2，只讓斜率與預測同比例偏；用 KHP 標定（20 mL 需約 0.41 g）。\n")
cat("6. 滴定量要設計到 15-40 mL；低酸度樣品加大 D（20 mL），高酸度樣品減小 D（5 mL）。\n")
cat("7. 比誤差時要分清楚實測與假設；NIR 殘差有一大部分可能來自參考值。\n")
cat("8. 已產生滴定紀錄表（20 列）；列印版流程見 ch08_titration_sop.md。\n")
cat("圖已存到", basename(FIG_DIR), "\n")
