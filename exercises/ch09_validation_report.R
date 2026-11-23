# =====================================================================
# 第 9 章 方法確效報告：現有證據夠不夠？R 程式練習（只用 base R）
# 前面八章的結果收在一起，回答：「這個 1 mm、2258 nm 的二階導數 NIR 方法，現在能不能說『已確效』？」
#   用途與目標、選擇性與範圍、線性、準確度、精密度、LOD/LOQ、預測區間（現在 vs 下一次實驗 vs 目標），
#   最後產生一份「確效總表」ch09_validation_scorecard.csv（後兩欄留給學生做完下一次實驗後填）。
#   配套的報告範本：ch09_validation_report_template.md。
# 所有計算都是確定性的（沒有任何亂數）。
#
# 執行方式：Rscript ch09_validation_report.R （程式會自動切到本檔所在資料夾）
# 資料：../data/spectra_long.csv  圖：figures_R/ch09_*.png（圖上文字一律用英文）
# 輸出：ch09_validation_scorecard.csv（UTF-8 含 BOM，Excel 可直接開）
# 每節結尾有「動手改改看」，每節最後有該節的自我驗證 PASS/FAIL。
#
# 數字的來源（和第 9 章講義 CH9 一致）：
#   【重算】從 spectra_long.csv 重算：二階導數 R²、單波長留一法 RMSECV、1 mm 2258 nm 最大吸光度、
#           殘差法與儀器雜訊法的 LOD/LOQ、純誤差與缺適性的自由度。
#   【取自前面各章】PLS RMSECV（第 4 章）、Atago 兩天差異 SD（第 8 章）、現在的預測區間（第 5 章）、
#           下一次實驗的區間半寬（第 7 章）、掃描 RMS 雜訊（第 6 章）；在程式裡寫死並註明出處。
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
WL1    <- 1670
TARGET_HW <- 0.5          # 示意的用途目標：預測區間半寬 +/- 0.5 %（不是法規或客戶要求）

# ---- 取自前面各章的數字（寫死；括號內是出處，講義 CH9 用的是同一批）----
PLS_RMSECV   <- 0.317     # 第 4 章：1 mm 二階導數、PLS 1 成分、4 瓶留一法 RMSECV
ATAGO_DAY_SD <- 0.4096    # 第 8 章：Atago 同 3 瓶隔天重測的差異 SD（% 酸度）
NOW_LO <- 3.157; NOW_HI <- 9.445; NOW_X0 <- 6.131   # 第 5 章：4 瓶滴定值校正、未知樣品 y0 = 4.20 的 95 % Fieller 區間（填 1 次）
PLAN_S     <- c(0.1, 0.2, 0.4, 0.6, 1.0)                           # 第 7 章：預試驗可能量到的填樣標準差（% 酸度）
PLAN_HW_M3 <- c(0.1274, 0.2548, 0.5107, 0.7688, 1.296)             # 第 7 章：8 濃度 x 3 次、未知樣品填 3 次、4 % 的區間半寬
PLAN_S_REF <- c(0.0, 0.2, 0.4)                                     # 第 7 章：再加上參考值（標準液）誤差，填樣 SD 固定 0.4 %
PLAN_HW_LV <- c(0.6033, 0.6307, 0.7086)                            # 第 7 章：對應的區間半寬
SCAN_RMS_AU <- 0.001236   # 第 6 章：同一瓶重複掃描的 RMS 雜訊（吸光度）
SIGMA_CH01  <- 2.7e-4     # 第 1 章：平坦區 1550-1650 nm 二階導數的雜訊 sigma（1 mm）

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
# 0. 載入資料、二階導數與共用函式
# =====================================================================
hr("0. 載入資料、二階導數與共用函式")
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
sg_d2 <- function(S, window = SG_WIN) {                  # 二階導數；1670、2258 nm 離邊緣很遠，邊緣怎麼處理不影響結果
  cf <- sg_coeffs2(window); half <- (window - 1) / 2
  out <- matrix(NA_real_, nrow(S), ncol(S))
  for (i in (half + 1):(ncol(S) - half)) out[, i] <- S[, (i - half):(i + half), drop = FALSE] %*% cf
  out[, 1:half] <- out[, half + 1]
  out[, (ncol(S) - half + 1):ncol(S)] <- out[, ncol(S) - half]
  out
}
D1 <- sg_d2(X1)
j2258 <- which.min(abs(wl - WL0)); j1670 <- which.min(abs(wl - WL1))
sig2258 <- D1[, j2258] * 1e3                              # 2258 nm 訊號 y（x 10^-3）
sig1670 <- D1[, j1670] * 1e3
cat(sprintf("1 mm、重複 1：%d 瓶 x %d 個波長；2258 nm 取 %.1f nm、1670 nm 取 %.1f nm\n", nrow(X1), length(wl), wl[j2258], wl[j1670]))
cat("2258 nm 訊號 y（x1e-3）：", round(sig2258, 4), "\n")
cat("1670 nm 訊號 y（x1e-3）：", round(sig1670, 4), "\n")

# 四組參考值（含水 = 0）
K <- 0.4315 / 40 / 0.1 * 10 * 0.006 * (100 / 25) * (100 / 10)   # 滴定：每 mL NaOH 對應的酸度 %（第 8 章，f = 1.07875）
refs <- list(atago  = c(0, 3.32, 4.87, 7.45),            # ATAGO PAL-Easy ACID 181（導電度式酸度計）11/02
             titr   = c(0, 8.8, 17, 44) * K,             # 滴定（NaOH 稱重配製、未標定、只做一次）
             label  = c(0, 3, 5, 10),                    # 瓶身標示
             pal30s = c(0, 2.8, 4.8, 10.4))              # PAL-30S 折射計，11/17 量，直接顯示酸度 %
ref_names <- c(atago = "Atago", titr = "titration", label = "label", pal30s = "PAL-30S")
fit_line <- function(xs, ys) {
  n <- length(xs); sxx <- sum((xs - mean(xs))^2)
  b <- sum((xs - mean(xs)) * (ys - mean(ys))) / sxx; a <- mean(ys) - b * mean(xs)
  s <- sqrt(sum((ys - a - b * xs)^2) / (n - 2))
  list(a = a, b = b, s = s, s_conc = s / abs(b), sxx = sxx)
}
r2f <- function(x, y) cor(x, y)^2
# 反向校正的留一法：每次留一瓶，用其餘 3 瓶擬合 訊號 = a + b x，再把被留下瓶的訊號換算回濃度，與它的參考值比
loo_rmse <- function(x, y) {
  e <- sapply(seq_along(x), function(i) { ft <- fit_line(x[-i], y[-i]); (y[i] - ft$a) / ft$b - x[i] })
  sqrt(mean(e^2))
}
cat(sprintf("滴定係數 K = %.5f %%/mL；滴定酸度（%%）：%s\n", K, paste(sprintf("%.3f", refs$titr), collapse = " / ")))

# ---- 動手改改看（第 0 節）----
# 1. 為什麼要把 1670 nm 也算進來？（提示：第 2 章用兩個醋酸帶互相對照）
# 2. 把 SG_WIN 改成 11，R^2 會變嗎？

# =====================================================================
# 1. 用途與目標
# =====================================================================
hr("1. 用途與目標")
cat("用途（這門課的假設）：用 1 mm 光徑、2258 nm 二階導數，由 NIR 光譜估計醋的酸度（w/v %）。\n")
cat(sprintf("目標（示意）：預測區間半寬 +/- %.1f %% 酸度。這個 %.1f %% 是為了讓後面能對照而設的『示意目標』，不是法規、不是客戶要求；\n", TARGET_HW, TARGET_HW))
cat("  真正的確效要先由用途決定可接受的誤差，再回頭看證據夠不夠。目標改了，結論也會改。\n")
cat("確效 = 事先寫好目標與判斷準則，再用實驗證據逐項檢查。現在只有 4 瓶（水 + 3 瓶醋）、每瓶 1 次，所以多半只能說『初步』。\n")
cat("各項要檢查：選擇性、範圍、線性、準確度、精密度、LOD/LOQ、穩健性、系統適用性（第 8 節的總表）。\n")
chk("number of calibration samples (water + 3 bottles)", nrow(X1), 4)
show_checks("第 1 節")

# ---- 動手改改看（第 1 節）----
# 1. 如果用途是『判斷醋是否在 4-6 %』，目標半寬 0.5 % 合理嗎？如果用途是『區分 3 % 與 10 %』呢？
# 2. 為什麼目標要在看數據『之前』訂？

# =====================================================================
# 2. 選擇性與範圍
# =====================================================================
hr("2. 選擇性與範圍：兩個醋酸帶、四種參考值的 R^2")
cat("Atago = ATAGO PAL-Easy ACID 181：導電度式酸度計，不是折射計。PAL-30S = 折射計（11/17 量），直接顯示酸度 %。\n")
r2_tab <- sapply(list("2258" = sig2258, "1670" = sig1670), function(y) sapply(refs, function(x) r2f(x, y)))
cat(sprintf("%-10s %10s %10s\n", "reference", "2258 nm", "1670 nm"))
for (k in names(refs)) cat(sprintf("%-10s %10.4f %10.4f\n", ref_names[[k]], r2_tab[k, "2258"], r2_tab[k, "1670"]))
cat("\n閱讀：\n")
cat("  1. 兩個醋酸帶（2258、1670 nm）對同一批參考值的 R^2 都很高，這只是『一致性』（兩個譜帶都跟濃度一起變），\n")
cat("     不是『選擇性』。選擇性是指加入糖、其他酸、色素等干擾物後，訊號仍只反映醋酸；這件事沒有測過。\n")
cat("  2. 標示與 PAL-30S 的 R^2（0.996、0.998）都比滴定（0.990）高；Atago（0.934）最低。\n")
cat("     只有 4 點（含水）的 R^2 無法用來選『哪個參考方法比較對』：R^2 高只代表 4 點排得像直線，\n")
cat("     不代表參考值準（第 8 章：比例誤差不改變 R^2）。同一組光譜，換參考值 R^2 就從 0.93 變到 0.998，\n")
cat("     差異來自參考值本身的誤差與非線性，不是 NIR 變好或變壞。\n")
Amax <- max(X1[, j2258]); A_LIMIT <- 1.5
cat(sprintf("\n範圍：1 mm、2258 nm 的最大原始吸光度 %.3f（< %.1f，沒有飽和）；已量到的範圍 Atago 0-%.2f %%、滴定 0-%.2f %%。\n",
            Amax, A_LIMIT, max(refs$atago), max(refs$titr)))
cat("  滴定的 11.39 % 只靠一次有爭議的滴定（第 8 章），不宜當作已證明的範圍上限；下一次實驗計畫 0-7 %、8 個等間隔濃度。\n")

png(file.path(FIG_DIR, "ch09_r2_by_reference.png"), width = 950, height = 520, res = 110)
par(mar = c(4, 4, 3, 1))
rm <- t(r2_tab[, c("2258", "1670")])
bp <- barplot(rm, beside = TRUE, names.arg = unname(ref_names[colnames(rm)]), col = c("steelblue", "grey70"), ylim = c(0.9, 1.012), xpd = FALSE,
              ylab = "R^2 (2nd derivative vs reference)", main = "Same spectra, different reference: R^2 with n = 4 cannot pick a reference")
text(bp, rm, sprintf("%.3f", rm), pos = 3, cex = 0.8)
legend("topleft", c("2258 nm", "1670 nm"), fill = c("steelblue", "grey70"), cex = 0.9)
dev.off()

key_r2 <- list("2258" = c(atago = 0.9342, titr = 0.9898, label = 0.9955, pal30s = 0.9984),
               "1670" = c(atago = 0.9448, titr = 0.9866, label = 0.9983, pal30s = 0.999))
for (w in c("2258", "1670")) for (k in names(refs)) chk(sprintf("R^2 %s nm vs %s", w, ref_names[[k]]), r2_tab[k, w], key_r2[[w]][[k]])
chk("label and PAL-30S R^2 both above titration (2258 nm)", as.numeric(r2_tab["label", "2258"] > r2_tab["titr", "2258"] && r2_tab["pal30s", "2258"] > r2_tab["titr", "2258"]), 1)
chk("max raw absorbance, 1 mm at 2258 nm", Amax, 1.441)
chk("calibration range, titration upper end (%)", round(max(refs$titr), 2), 11.39)
show_checks("第 2 節")

# ---- 動手改改看（第 2 節）----
# 1. 『兩個譜帶都跟濃度一起變』為什麼不能證明選擇性？要加什麼樣品才能檢驗？
# 2. 標示的 R^2 比滴定高，能不能推論『標示比滴定準』？
# 3. 11.39 % 這個點拿掉後，滴定的範圍變成多少？

# =====================================================================
# 3. 線性
# =====================================================================
hr("3. 線性：4 點、沒有重複")
fa <- fit_line(refs$atago, sig2258); ft <- fit_line(refs$titr, sig2258)
n_obs <- 4; n_lvl <- 4
pure_df_now <- n_obs - n_lvl; resid_df_now <- n_obs - 2
cat(sprintf("現在：%d 瓶、%d 個濃度、每濃度 1 次。殘差自由度 = %d；純誤差自由度 = n - 濃度數 = %d；缺適性自由度 = %d。\n",
            n_obs, n_lvl, resid_df_now, pure_df_now, resid_df_now - pure_df_now))
cat(sprintf("R^2：Atago %.4f、滴定 %.4f（2258 nm）。沒有重複，就沒有『純誤差』可以拿來檢定直線是否夠用（缺適性檢定）。\n", r2f(refs$atago, sig2258), r2f(refs$titr, sig2258)))
n_plan_lvl <- 8; n_plan_rep <- 3; n_plan <- n_plan_lvl * n_plan_rep
pure_df_plan <- n_plan - n_plan_lvl; lof_df_plan <- n_plan_lvl - 2
cat(sprintf("下一次（第 7 章）：%d 濃度 x %d 次 = %d 觀測；純誤差自由度 %d；缺適性自由度 %d；才有條件檢定直線與曲線。\n",
            n_plan_lvl, n_plan_rep, n_plan, pure_df_plan, lof_df_plan))
cat("提醒：同一濃度的重複，如果是『同一次配製、同一瓶、換個填樣』，它們共用配製誤差；\n")
cat("      這樣得到的純誤差只含填樣與儀器誤差，會低估真正的樣品間變異。重複要盡量包含獨立配製。\n")
cat("      同樣地，現在 4 瓶各只配一次，R^2 同時包含配製、填樣、參考值與光譜誤差，無法拆開。\n")
chk("R^2 Atago (lin)", r2f(refs$atago, sig2258), 0.9342)
chk("R^2 titration (lin)", r2f(refs$titr, sig2258), 0.9898)
chk("pure-error df, now", pure_df_now, 0)
chk("pure-error df, plan (8 x 3)", pure_df_plan, 16)
chk("lack-of-fit df, plan (8 levels)", lof_df_plan, 6)
show_checks("第 3 節")

# ---- 動手改改看（第 3 節）----
# 1. 如果 8 個濃度只各做 1 次（8 觀測），純誤差自由度是多少？缺適性自由度呢？
# 2. 『同瓶重複填 3 次』與『3 瓶獨立配製』的純誤差，哪個比較大？為什麼？

# =====================================================================
# 4. 準確度
# =====================================================================
hr("4. 準確度：留一法 RMSECV（% 酸度）")
cat("做法：4 瓶輪流留下 1 瓶，用其餘 3 瓶擬合 訊號 = a + b x，再把被留下瓶的訊號換算回濃度；RMSECV = 這 4 個誤差的均方根。\n")
cat("注意：這不是『準確度』：沒有獨立樣品（n_test = 0），而且每次只用 3 點擬合（自由度 1）。它是『現有資料能做到的內部檢查』。\n")
rm_loo <- sapply(refs, function(x) loo_rmse(x, sig2258))
cat(sprintf("%-10s %9s\n", "reference", "RMSECV"))
for (k in names(refs)) cat(sprintf("%-10s %9.4f\n", ref_names[[k]], rm_loo[[k]]))
cat(sprintf("PLS（第 4 章，1 mm 二階導數、1 成分）RMSECV = %.3f %%\n", PLS_RMSECV))
cat(sprintf("Atago 兩天差異 SD（第 8 章，同樣 3 瓶隔天重測）= %.4f %%\n", ATAGO_DAY_SD))
cat("\n閱讀：\n")
cat(sprintf("  1. PLS 的 RMSECV（%.2f %%）比 Atago 自己兩天的差異（%.2f %%）還小：NIR 預測得比參考方法自己重複的結果還一致，\n", PLS_RMSECV, ATAGO_DAY_SD))
cat("     這不合理；是過度配適（用 4 個點、整段光譜擬合）的徵兆。RMSECV 不能小於參考值本身的誤差還說它『準』。\n")
cat("  2. PAL-30S 的 0.26 % 很小，但只代表 NIR 與 PAL-30S 兩個方法『彼此一致』（consistency），不代表準確（accuracy）：\n")
cat("     兩者有可能一起偏。準確度需要獨立的盲樣加上有標定過的滴定（第 7、8 章）。\n")
cat("  3. 同一組光譜，RMSECV 隨參考值從 2.18（Atago）到 0.26（PAL-30S）變化 8 倍，所以 RMSECV 的大小主要反映參考值的品質。\n")

png(file.path(FIG_DIR, "ch09_rmsecv_bars.png"), width = 1000, height = 540, res = 110)
par(mar = c(5, 4, 3, 1))
bv <- c(Atago = rm_loo[["atago"]], titration = rm_loo[["titr"]], label = rm_loo[["label"]], "PLS (ch04)" = PLS_RMSECV, "PAL-30S" = rm_loo[["pal30s"]])
bp <- barplot(bv, col = c("grey55", "grey55", "grey55", "steelblue", "grey55"), ylab = "RMSECV (% acid)", ylim = c(0, 2.6),
              main = "Leave-one-out RMSECV by reference (n = 4)", cex.names = 0.9)
text(bp, bv / 2, sprintf("%.2f", bv), col = "white", cex = 1.0, font = 2)
abline(h = ATAGO_DAY_SD, lty = 2, col = "firebrick", lwd = 2)
legend("topright", c("single wavelength 2258 nm", "PLS, all wavelengths (ch04)", sprintf("Atago day-to-day SD = %.2f", ATAGO_DAY_SD)),
       fill = c("grey55", "steelblue", NA), border = c("black", "black", NA), lty = c(NA, NA, 2), lwd = c(NA, NA, 2), col = c(NA, NA, "firebrick"), cex = 0.8)
dev.off()

chk("RMSECV Atago, 2258 nm (% acid)", rm_loo[["atago"]], 2.182)
chk("RMSECV titration, 2258 nm (% acid)", rm_loo[["titr"]], 1.152)
chk("RMSECV PAL-30S, 2258 nm (% acid)", rm_loo[["pal30s"]], 0.2559)
chk("RMSECV label, 2258 nm (% acid; ch08 nirfit)", rm_loo[["label"]], 0.6769)
chk("PLS RMSECV (ch04, hard-coded)", PLS_RMSECV, 0.317)
chk("Atago day-to-day SD (ch08, hard-coded)", ATAGO_DAY_SD, 0.4096)
chk("PLS RMSECV below Atago day SD (1 = yes)", as.numeric(PLS_RMSECV < ATAGO_DAY_SD), 1)
chk("independent test samples (n_test)", 0, 0)
show_checks("第 4 節")

# ---- 動手改改看（第 4 節）----
# 1. 如果 Atago 兩天差異是 0.41 %，NIR 的 RMSECV 理論上最好能做到多少？（提示：參考值誤差進了校正）
# 2. 為什麼 PAL-30S 的 RMSECV 小，還不能說 NIR『準確』？
# 3. 要得到真正的準確度，需要什麼樣的樣品？

# =====================================================================
# 5. 精密度
# =====================================================================
hr("5. 精密度：量過什麼、沒量什麼")
resid_titr <- sqrt(sum((sig2258 - (ft$a + ft$b * refs$titr))^2) / 2) / abs(ft$b)   # 校正殘差 SD（% 酸度，4 點、自由度 2）
cat("【量過的】\n")
cat(sprintf("  - 掃描雜訊：同一瓶重複掃描的 RMS = %.6f AU（第 6 章）。這只是儀器的短期重複性。\n", SCAN_RMS_AU))
cat(sprintf("  - 校正殘差 SD（滴定當參考、4 點、自由度 2）= %.3f %%；它混了光譜、配製、填樣與參考值誤差，不是『精密度』。\n", resid_titr))
cat(sprintf("  - 參考方法 Atago 隔天差異 SD = %.4f %%（第 8 章）；這是參考方法的日間精密度，不是 NIR 的。\n", ATAGO_DAY_SD))
cat("【沒量的】\n")
cat("  - 填樣重複性（同一瓶重新填 10 次）\n")
cat("  - 樣品配製的重複性（獨立配 3 份）\n")
cat("  - 日間（不同天）、不同操作者、不同樣品槽的中間精密度\n")
cat("下一次（第 7 章預試驗）先填樣 10 次量出填樣 SD，因為它直接決定預測區間多寬（第 7 節）。\n")
chk("scan RMS noise (AU; ch06, hard-coded)", SCAN_RMS_AU, 0.001236)
chk("NIR calibration residual SD, titration ref (% acid)", resid_titr, 0.612)
show_checks("第 5 節")

# ---- 動手改改看（第 5 節）----
# 1. 掃描雜訊 0.0012 AU 換算成酸度大概是多少 %？為什麼它不等於方法的精密度？
# 2. 填樣、日間、操作者三種變異，哪一種你最擔心？怎麼用最少次實驗量出來？

# =====================================================================
# 6. LOD / LOQ
# =====================================================================
hr("6. LOD / LOQ：殘差法 vs 儀器雜訊法")
cat("殘差法：LOD = 3 s / |b|、LOQ = 10 s / |b|；s = 校正線的殘差 SD（訊號單位），b = 斜率。s 包含了參考值與配製的誤差。\n")
lod_res <- function(x, y) { f <- fit_line(x, y); c(3 * f$s_conc, 10 * f$s_conc) }
res_titr <- lod_res(refs$titr, sig2258); res_atago <- lod_res(refs$atago, sig2258)
cat(sprintf("  滴定當參考  ：LOD %.3f %%、LOQ %.3f %%\n", res_titr[1], res_titr[2]))
cat(sprintf("  Atago 當參考：LOD %.3f %%、LOQ %.3f %%（LOQ 已超過這 4 瓶的最大濃度 %.2f %%）\n", res_atago[1], res_atago[2], max(refs$atago)))
slope_raw <- fit_line(refs$atago, D1[, j2258])$b           # 第 1 章用的斜率：原始（沒乘 1e3）二階導數對 Atago 濃度
lod_ch01 <- 3 * SIGMA_CH01 / abs(slope_raw); loq_ch01 <- 10 * SIGMA_CH01 / abs(slope_raw)
cat(sprintf("儀器雜訊法（第 1 章）：sigma = %.1e（平坦區 1550-1650 nm 的二階導數雜訊，1 mm）；斜率（原始、未乘 1e3）= %.7f 每 %%\n", SIGMA_CH01, slope_raw))
cat(sprintf("  LOD = 3 sigma / |b| = %.3f %%、LOQ = 10 sigma / |b| = %.3f %%\n", lod_ch01, loq_ch01))
cat("閱讀：兩種做法差 2-4 倍。儀器雜訊法只考慮光譜雜訊，是下限；殘差法把配製、參考值的誤差都算進去，比較保守但也被 4 點、自由度 2 拖累。\n")
cat("所以報告時要寫清楚『用哪一種方法、s 是什麼』，而且新校正做完要重算。\n")

png(file.path(FIG_DIR, "ch09_lod_loq.png"), width = 950, height = 520, res = 110)
par(mar = c(4, 4, 3, 1))
lm_ <- rbind(LOD = c(res_titr[1], res_atago[1], lod_ch01), LOQ = c(res_titr[2], res_atago[2], loq_ch01))
bp <- barplot(lm_, beside = TRUE, names.arg = c("residual, titration ref", "residual, Atago ref", "instrument noise (ch01)"),
              col = c("steelblue", "firebrick"), ylab = "% acid", ylim = c(0, 11.5), main = "LOD and LOQ: three ways to estimate (preliminary)")
text(bp, lm_, sprintf("%.2f", lm_), pos = 3, cex = 0.85)
legend("topright", c("LOD (3 s / |b|)", "LOQ (10 s / |b|)"), fill = c("steelblue", "firebrick"), cex = 0.85)
dev.off()

chk("residual LOD, titration (%)", res_titr[1], 1.836)
chk("residual LOQ, titration (%)", res_titr[2], 6.12)
chk("residual LOD, Atago (%)", res_atago[1], 3.038)
chk("residual LOQ, Atago (%)", res_atago[2], 10.13)
chk("ch01 slope, 1 mm 2258 nm, raw units (AU/nm^2 per %)", slope_raw, -0.0009852)
chk("instrument-noise LOD, ch01 (%)", lod_ch01, 0.8222)
chk("instrument-noise LOQ, ch01 (%)", loq_ch01, 2.741)
show_checks("第 6 節")

# ---- 動手改改看（第 6 節）----
# 1. 如果校正殘差 SD 減半，殘差法的 LOD 與 LOQ 變多少？
# 2. LOQ 大於已量到的最大濃度，這個 LOQ 還可信嗎？

# =====================================================================
# 7. 預測區間：現在 vs 下一次實驗 vs 目標
# =====================================================================
hr("7. 預測區間：現在 vs 下一次實驗 vs 目標（半寬，% 酸度）")
now_hw <- (NOW_HI - NOW_LO) / 2
cat(sprintf("現在（第 5 章）：4 瓶、未知樣品填 1 次，預測 %.3f %%，95 %% 區間 %.3f 到 %.3f %%，半寬 %.3f %%。\n", NOW_X0, NOW_LO, NOW_HI, now_hw))
cat(sprintf("目標（示意）：半寬 %.1f %%；現在是目標的 %.1f 倍，離目標很遠。\n", TARGET_HW, now_hw / TARGET_HW))
cat("下一次（第 7 章）：8 濃度 x 3 次、未知樣品填 3 次、4 %；區間半寬取決於預試驗量到的填樣 SD：\n")
cat(sprintf("  %-26s %s\n", "填樣 SD（%）", paste(sprintf("%6.1f", PLAN_S), collapse = " ")))
cat(sprintf("  %-26s %s\n", "區間半寬（%）", paste(sprintf("%6.3f", PLAN_HW_M3), collapse = " ")))
cat(sprintf("  %-26s %s\n", "達到目標 0.5 %？", paste(sprintf("%6s", ifelse(PLAN_HW_M3 <= TARGET_HW, "yes", "no")), collapse = " ")))
cat("  再加上參考值（標準液）誤差（填樣 SD 固定 0.4 %）：\n")
cat(sprintf("  %-26s %s\n", "參考值 SD（%）", paste(sprintf("%6.1f", PLAN_S_REF), collapse = " ")))
cat(sprintf("  %-26s %s\n", "區間半寬（%）", paste(sprintf("%6.3f", PLAN_HW_LV), collapse = " ")))
cat("閱讀：\n")
cat(sprintf("  1. 8 x 3 的設計能把半寬從 %.2f %% 壓到 0.13-1.3 %%，目標 0.5 %% 只有在填樣 SD 約 0.4 %% 以下才有機會（0.4 %% 時 0.51，略超過）。\n", now_hw))
cat("  2. 只要把參考值誤差算進去（標準液 SD 0.2 %、0.4 %），半寬就從 0.60 變成 0.63、0.71，全都超過 0.5 %：\n")
cat("     到那時，瓶頸是參考方法，不是光譜。\n")

png(file.path(FIG_DIR, "ch09_interval_compare.png"), width = 1000, height = 560, res = 110)
par(mar = c(4, 14, 3, 1))
iv <- c(now_hw, PLAN_HW_M3, PLAN_HW_LV)
inm <- c("now: 4 bottles, 1 fill (ch05)", sprintf("plan 8x3, fill SD %.1f", PLAN_S), sprintf("plan 8x3, fill 0.4 + ref SD %.1f", PLAN_S_REF))
icol <- c("firebrick", rep("steelblue", 5), rep("darkorange", 3))
plot(iv, rev(seq_along(iv)), log = "x", xlim = c(0.08, 4.5), ylim = c(0.5, length(iv) + 0.5), pch = 19, col = icol, cex = 1.6, yaxt = "n", ylab = "",
     xlab = "95 % prediction interval half-width (% acid, log scale)", main = "Prediction interval: now vs next experiment vs target")
segments(0.08, rev(seq_along(iv)), iv, rev(seq_along(iv)), col = icol, lwd = 2)
axis(2, at = rev(seq_along(iv)), labels = inm, las = 1, cex.axis = 0.78)
text(iv, rev(seq_along(iv)), sprintf("%.2f", iv), cex = 0.8, pos = 4)
abline(v = TARGET_HW, lty = 2, lwd = 2, col = "darkgreen")
text(TARGET_HW, length(iv) + 0.45, "illustrative target 0.5", col = "darkgreen", cex = 0.8, pos = 4)
dev.off()

chk("now: interval lower limit (%; ch05, hard-coded)", NOW_LO, 3.157)
chk("now: interval upper limit (%; ch05, hard-coded)", NOW_HI, 9.445)
chk("now: prediction x0 (%; ch05, hard-coded)", NOW_X0, 6.131)
chk("now: half-width (%)", now_hw, 3.144)
for (i in seq_along(PLAN_S)) chk(sprintf("plan half-width, fill SD %.1f (%%)", PLAN_S[i]), PLAN_HW_M3[i], c(0.1274, 0.2548, 0.5107, 0.7688, 1.296)[i])
for (i in seq_along(PLAN_S_REF)) chk(sprintf("plan half-width, ref SD %.1f + fill 0.4 (%%)", PLAN_S_REF[i]), PLAN_HW_LV[i], c(0.6033, 0.6307, 0.7086)[i])
chk("plan designs within the 0.5 % target (of 8)", sum(c(PLAN_HW_M3, PLAN_HW_LV) <= TARGET_HW), 2)
chk("now half-width / target", now_hw / TARGET_HW, 6.288)
show_checks("第 7 節")

# ---- 動手改改看（第 7 節）----
# 1. 如果把目標放寬到 +/- 1 %，下一次的 8 種情況有幾個達標？
# 2. 為什麼說『加入參考值誤差後，瓶頸是參考方法』？要改善要做什麼（提示：第 8 章）？
# 3. 現在的區間是 3.16-9.45 %，覆蓋了 3 % 到 10 % 的大部分標示值；這樣的方法有用途嗎？

# =====================================================================
# 8. 確效總表與報告範本
# =====================================================================
hr("8. 產生 ch09_validation_scorecard.csv（確效總表）")
status_txt <- c("沒有", "部分或初步", "有證據")             # 0 / 1 / 2
card <- data.frame(
  item = c("選擇性", "線性", "範圍", "準確度", "精密度", "LOD／LOQ", "穩健性", "系統適用性"),
  code = c(0, 1, 1, 0, 0, 1, 0, 0),
  now  = c("兩個醋酸帶一致只是一致性；沒有測過糖等干擾物", "4 點、沒有重複，無法做缺適性檢定",
           "2258 nm 在 1 mm 未飽和；11.4 % 只靠一次有爭議的滴定", "沒有獨立樣品；參考值彼此對不上",
           "只有掃描雜訊；填樣、日間誤差沒量", "殘差法 1.8 / 6.1 %、儀器雜訊法 0.8 %（初步）",
           "溫度、操作者、樣品槽都沒測", "水帶規則在 1 mm 只能當參考（水帶飽和）"),
  how  = c("市售盲樣、加入干擾物（第 7 章）", "8 濃度 × 3 次（第 7 章）", "0–7 %、8 個等間隔濃度", "盲樣＋標定過的滴定（第 7、8 章）",
           "預試驗填樣 10 次（第 7 章）", "新校正重算", "刻意改變條件的小實驗", "每次實驗量共同對照（第 7 章）"),
  stringsAsFactors = FALSE)
cols9 <- c("項目", "狀態", "現在的證據", "怎麼補", "下一次實驗結果", "結論")
body <- unlist(lapply(seq_len(nrow(card)), function(i) {
  fields <- c(card$item[i], status_txt[card$code[i] + 1], card$now[i], card$how[i], "", "")
  stopifnot(!any(grepl("[,\"\n]", fields)))              # 欄位內不含逗號、引號、換行，不需要加引號（三個語言輸出才會逐位元一致）
  paste0(paste(fields, collapse = ","), "\n")
}))
txt <- paste0(paste(cols9, collapse = ","), "\n", paste(body, collapse = ""))
con <- file("ch09_validation_scorecard.csv", "wb")
writeBin(c(as.raw(c(0xEF, 0xBB, 0xBF)), charToRaw(enc2utf8(txt))), con); close(con)
sc <- read.csv("ch09_validation_scorecard.csv", fileEncoding = "UTF-8-BOM", check.names = FALSE, colClasses = "character")
cat(sprintf("已寫出 ch09_validation_scorecard.csv：%d 列，%d 欄\n", nrow(sc), ncol(sc)))
cat(sprintf("狀態：沒有 %d 項、部分或初步 %d 項、有證據 %d 項（8 項裡沒有任何一項是『有證據』）。\n",
            sum(sc[["狀態"]] == "沒有"), sum(sc[["狀態"]] == "部分或初步"), sum(sc[["狀態"]] == "有證據")))
cat("『下一次實驗結果』『結論』兩欄留白，做完第 7 章的實驗後由學生填寫。\n")
cat("報告範本：ch09_validation_report_template.md（一個驗證項目一個小節，現有證據已預填；實驗後把『結果』『結論』補上）。\n")
tpl <- readLines("ch09_validation_report_template.md", encoding = "UTF-8", warn = FALSE)
cat("\n能說什麼：\n")
cat("  - 1 mm、2258 nm 的二階導數與滴定值（n = 4）呈線性關係；這是一個可行性研究（feasibility study）的初步證據。\n")
cat("不能說什麼：\n")
cat("  - 『方法已確效』；『準確度是 x %』（沒有獨立樣品）；『可取代滴定』（參考方法彼此對不上，NIR 本身依賴它們校正）。\n")

chk("scorecard: number of rows", nrow(sc), 8)
chk("scorecard: number of columns", ncol(sc), 6)
chk("scorecard: status 沒有 (CH9.card 狀態 0)", sum(sc[["狀態"]] == "沒有"), 5)
chk("scorecard: status 部分或初步 (CH9.card 狀態 1)", sum(sc[["狀態"]] == "部分或初步"), 3)
chk("scorecard: status 有證據 (CH9.card 狀態 2)", sum(sc[["狀態"]] == "有證據"), 0)
chk("scorecard: per-item status mismatches vs CH9.card", sum(sc[["狀態"]] != c("沒有", "部分或初步", "部分或初步", "沒有", "沒有", "部分或初步", "沒有", "沒有")), 0)
chk("scorecard: blank 下一次實驗結果 cells", sum(is.na(sc[["下一次實驗結果"]]) | sc[["下一次實驗結果"]] == ""), 8)
chk("scorecard: blank 結論 cells", sum(is.na(sc[["結論"]]) | sc[["結論"]] == ""), 8)
chk("template: item names found in report template (of 8)", sum(sapply(card$item, function(s) any(grepl(s, tpl, fixed = TRUE)))), 8)
show_checks("第 8 節")

# ---- 動手改改看（第 8 節）----
# 1. 用 Excel 開 ch09_validation_scorecard.csv，中文有亂碼嗎？（檔案含 BOM，應該不會）
# 2. 做完第 7 章的實驗後，哪幾項的狀態可以從『沒有』升到『部分或初步』？哪一項仍然升不上去？
# 3. 為什麼『有證據』一項都沒有，卻還是值得做這個方法？

# =====================================================================
# 總結
# =====================================================================
hr("總結")
cat(sprintf("全部自我驗證：%d / %d PASS\n", sum(checks$ok), nrow(checks)))
if (!all(checks$ok)) { cat("FAIL 項目：\n"); print(checks[!checks$ok, ], row.names = FALSE) }
cat("1. 目標 +/- 0.5 % 是示意；先有用途與目標，才能判斷證據夠不夠。\n")
cat("2. 兩個醋酸帶一致是一致性不是選擇性；n = 4 的 R^2 不能用來選參考方法。\n")
cat("3. 線性：4 點、沒有重複，純誤差自由度 0；下一次 8 x 3 才有 16 與 6。\n")
cat("4. 準確度：沒有獨立樣品；PLS RMSECV 0.32 < Atago 日間 SD 0.41 是過度配適的徵兆；PAL-30S 0.26 是一致性不是準確度。\n")
cat("5. 精密度：只量了掃描雜訊；填樣、配製、日間都沒量。\n")
cat("6. LOD/LOQ：殘差法 1.8 / 6.1 %（滴定）、儀器雜訊法 0.8 / 2.7 %，皆為初步。\n")
cat("7. 預測區間：現在半寬 3.14 %，下一次 0.13-1.3 %；示意目標 0.5 % 只有填樣 SD 夠小、且參考值誤差很小時才可能達到。\n")
cat("8. 確效總表 8 項：沒有 5、部分或初步 3、有證據 0；可說『可行性研究、初步線性』，不能說『已確效』。\n")
cat("圖已存到", basename(FIG_DIR), "\n")
