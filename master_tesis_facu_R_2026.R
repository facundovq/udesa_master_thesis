##############################
###        Tesis           ###
##############################

# ============================================================================
# SETUP
# ============================================================================

# install.packages(c("tidyverse","zoo","urca","vars","sandwich","lmtest","stargazer"))

library(tidyverse)
library(zoo)
library(urca)
library(vars)
library(sandwich)
library(lmtest)
library(stargazer)

# Force dplyr functions after vars/MASS overwrite them
select <- dplyr::select
filter <- dplyr::filter
lag    <- dplyr::lag
lead   <- dplyr::lead

data_path   <- "C:/Users/cufa9/Documents/Master/UdeSA/Tesis/input/"
output_path <- "C:/Users/cufa9/Documents/Master/UdeSA/Tesis/output/"


# ============================================================================
# 1. DATA — LOAD & PREPARE
# ============================================================================

cpi <- read.csv(paste0(data_path, "CPIAUCSL.csv"))
psr <- read.csv(paste0(data_path, "PSAVERT.csv"))
yc  <- read.csv(paste0(data_path, "T10Y3M.csv"))
rec <- read.csv(paste0(data_path, "USREC.csv"))

cpi$observation_date <- as.Date(cpi$observation_date)
psr$observation_date <- as.Date(psr$observation_date)
yc$observation_date  <- as.Date(yc$observation_date)
rec$observation_date <- as.Date(rec$observation_date)

yc$T10Y3M    <- as.numeric(as.character(yc$T10Y3M))
psr$PSAVERT  <- as.numeric(as.character(psr$PSAVERT))
cpi$CPIAUCSL <- as.numeric(as.character(cpi$CPIAUCSL))

cat("NA count — T10Y3M:", sum(is.na(yc$T10Y3M)),
    "| PSAVERT:", sum(is.na(psr$PSAVERT)),
    "| CPIAUCSL:", sum(is.na(cpi$CPIAUCSL)), "\n")

# --- Aggregate YC from daily to monthly ---
yc_monthly <- yc %>%
  mutate(observation_date = as.Date(paste0(format(observation_date, "%Y-%m"), "-01"))) %>%
  group_by(observation_date) %>%
  summarise(T10Y3M = mean(T10Y3M, na.rm = TRUE), .groups = "drop")

cat("Monthly YC observations:", nrow(yc_monthly), "\n")

# --- Filter & merge ---
start_date <- as.Date("1982-01-01")
end_date   <- as.Date("2026-02-01")

df_final <- cpi %>%
  filter(observation_date >= start_date & observation_date <= end_date) %>%
  inner_join(psr %>% filter(observation_date >= start_date & observation_date <= end_date),
             by = "observation_date") %>%
  inner_join(yc_monthly %>% filter(observation_date >= start_date & observation_date <= end_date),
             by = "observation_date")

cat("Final dataset:", nrow(df_final), "x", ncol(df_final), "\n")

# --- Plot raw data ---
ts_raw <- zoo(df_final[, c("CPIAUCSL","PSAVERT","T10Y3M")],
              order.by = df_final$observation_date)
plot(ts_raw, main = "Macroeconomic Variables (1982-2026)", xlab = "Time",
     col = c("steelblue","darkorange","forestgreen"))


# ============================================================================
# 2. STATIONARITY & AUTOCORRELATION
# ============================================================================

df_final <- df_final %>%
  arrange(observation_date) %>%
  mutate(Inflation = (log(CPIAUCSL) - log(lag(CPIAUCSL))) * 1200) %>%
  drop_na()

ts_trans <- zoo(df_final[, c("Inflation","PSAVERT","T10Y3M")],
                order.by = df_final$observation_date)
plot(ts_trans, main = "Transformed Macroeconomic Variables", xlab = "Time",
     col = c("steelblue","darkorange","forestgreen"))

# --- ADF Tests ---
cat("\n=== ADF Test: Inflation ===\n")
summary(ur.df(df_final$Inflation, type = "drift", selectlags = "AIC"))

cat("\n=== ADF Test: Personal Saving Rate ===\n")
summary(ur.df(df_final$PSAVERT, type = "drift", selectlags = "AIC"))

cat("\n=== ADF Test: Yield Curve Spread ===\n")
summary(ur.df(df_final$T10Y3M, type = "none", selectlags = "AIC"))

# --- ACF & PACF ---
par(mfrow = c(3, 2))
acf(df_final$Inflation,  main = "ACF: Inflation")
pacf(df_final$Inflation, main = "PACF: Inflation")
acf(df_final$PSAVERT,    main = "ACF: Personal Saving Rate")
pacf(df_final$PSAVERT,   main = "PACF: Personal Saving Rate")
acf(df_final$T10Y3M,     main = "ACF: Yield Curve Spread")
pacf(df_final$T10Y3M,    main = "PACF: Yield Curve Spread")
par(mfrow = c(1, 1))

# --- Diagnostic plot: saving rate vs inversions ---
df_final$inversion <- as.integer(df_final$T10Y3M < 0)

par(mar = c(4, 4, 3, 4))
plot(df_final$observation_date, df_final$PSAVERT,
     type = "l", col = "darkorange", lwd = 1.5,
     ylab = "Personal Saving Rate (%)", xlab = "",
     main = "Personal Saving Rate vs Yield Curve (shaded = inversion)")

inv_runs <- rle(df_final$inversion)
idx <- cumsum(c(1, inv_runs$lengths))
for (i in seq_along(inv_runs$values)) {
  if (inv_runs$values[i] == 1) {
    rect(df_final$observation_date[idx[i]], par("usr")[3],
         df_final$observation_date[min(idx[i+1], nrow(df_final))], par("usr")[4],
         col = rgb(1,0,0,0.15), border = NA)
  }
}
par(new = TRUE)
plot(df_final$observation_date, df_final$T10Y3M,
     type = "l", col = "steelblue", lwd = 1, axes = FALSE, xlab = "", ylab = "")
abline(h = 0, lty = 2, col = "black", lwd = 0.8)
axis(4, col.axis = "steelblue")
mtext("T10Y3M Spread (pp)", side = 4, line = 3, col = "steelblue")

cat("Inversion episodes:", sum(diff(c(0, df_final$inversion)) == 1),
    "| Total inversion months:", sum(df_final$inversion), "\n")
par(mar = c(5, 4, 4, 2))


# ============================================================================
# 3. LOCAL PROJECTIONS (Jorda 2005)
# ============================================================================

# --- BIC lag selection via VAR ---
lp_vars   <- df_final[, c("Inflation","PSAVERT","T10Y3M")]
var_model <- VAR(lp_vars, type = "const", lag.max = 12, ic = "SC")
p         <- var_model$p
cat("BIC-selected lag order: p =", p, "\n")

# --- Build LP dataset ---
H <- 36

df_lp <- df_final %>%
  select(observation_date, Inflation, PSAVERT, T10Y3M) %>%
  mutate(
    D_inversion = as.integer(T10Y3M < 0),
    D_covid     = as.integer(observation_date >= as.Date("2020-03-01") &
                               observation_date <= as.Date("2021-06-01")),
    trend       = seq_len(n())
  ) %>%
  left_join(rec %>% rename(D_recession = USREC), by = "observation_date")

# Lags — base R shift to avoid namespace issues inside loops
for (l in seq_len(p)) {
  n <- nrow(df_lp)
  df_lp[[paste0("Inflation_L", l)]] <- c(rep(NA, l), df_lp$Inflation[1:(n-l)])
  df_lp[[paste0("PSAVERT_L",   l)]] <- c(rep(NA, l), df_lp$PSAVERT[1:(n-l)])
  df_lp[[paste0("T10Y3M_L",    l)]] <- c(rep(NA, l), df_lp$T10Y3M[1:(n-l)])
}

# Outcome: PSAVERT_{t+h} - PSAVERT_{t-1}
df_lp$PSAVERT_Lminus1 <- df_lp$PSAVERT_L1

for (h in 0:H) {
  n <- nrow(df_lp)
  future_psr <- c(df_lp$PSAVERT[(h+1):n], rep(NA, h))
  df_lp[[paste0("y_h", h)]] <- future_psr - df_lp$PSAVERT_Lminus1
}

df_lp <- df_lp %>% drop_na()
cat("LP dataset:", nrow(df_lp), "observations\n")

df_lp_pre2020 <- df_lp %>%
  filter(observation_date < as.Date("2020-01-01"))

cat("Full sample:", nrow(df_lp), "| Pre-2020:", nrow(df_lp_pre2020), "\n")

# --- Control column definitions ---
lag_cols     <- c(paste0("Inflation_L", 1:p),
                  paste0("PSAVERT_L",   1:p),
                  paste0("T10Y3M_L",    1:p))
contemp_ctrl <- c("PSAVERT", "Inflation")
outlier_ctrl <- c("D_covid", "trend", "D_recession")

# --- LP runner ---
run_lp <- function(data, se_type = "HC1", H = 36, alpha = 0.10) {
  res1 <- res2 <- vector("list", H + 1)
  avail_outlier <- intersect(outlier_ctrl, names(data))
  
  for (h in 0:H) {
    y <- data[[paste0("y_h", h)]]
    
    # Experiment 1: Inversion dummy
    rhs1 <- intersect(c("D_inversion", contemp_ctrl, avail_outlier, lag_cols), names(data))
    m1   <- lm(y ~ ., data = as.data.frame(cbind(y = y, data[, rhs1])))
    vcov1 <- if (se_type == "HC1") vcovHC(m1, type = "HC1") else
      NeweyWest(m1, lag = h + 1, prewhite = FALSE)
    ct1  <- coeftest(m1, vcov1)
    ci1  <- coefci(m1, vcov = vcov1, level = 1 - alpha)
    res1[[h + 1]] <- data.frame(
      h       = h,
      coef    = ct1["D_inversion", "Estimate"],
      ci_low  = ci1["D_inversion", 1],
      ci_high = ci1["D_inversion", 2],
      pval    = ct1["D_inversion", "Pr(>|t|)"]
    )
    
    # Experiment 2: YC level (1pp decline)
    rhs2 <- intersect(c("T10Y3M", contemp_ctrl, avail_outlier, lag_cols), names(data))
    m2   <- lm(y ~ ., data = as.data.frame(cbind(y = y, data[, rhs2])))
    vcov2 <- if (se_type == "HC1") vcovHC(m2, type = "HC1") else
      NeweyWest(m2, lag = h + 1, prewhite = FALSE)
    ct2  <- coeftest(m2, vcov2)
    ci2  <- coefci(m2, vcov = vcov2, level = 1 - alpha)
    shock <- -1
    res2[[h + 1]] <- data.frame(
      h       = h,
      coef    = ct2["T10Y3M", "Estimate"] * shock,
      ci_low  = min(ci2["T10Y3M", ] * shock),
      ci_high = max(ci2["T10Y3M", ] * shock),
      pval    = ct2["T10Y3M", "Pr(>|t|)"]
    )
  }
  
  list(irf1 = do.call(rbind, res1),
       irf2 = do.call(rbind, res2))
}

# --- Run all four models ---
cat("Running Model 1 (HC1, Full)...\n");     r1 <- run_lp(df_lp,         "HC1")
cat("Running Model 2 (HAC, Full)...\n");     r2 <- run_lp(df_lp,         "HAC")
cat("Running Model 3 (HAC, Pre-2020)...\n"); r3 <- run_lp(df_lp_pre2020, "HAC")
cat("Running Model 4 (HC1, Pre-2020)...\n"); r4 <- run_lp(df_lp_pre2020, "HC1")
cat("All models estimated.\n")

irf1_m1 <- r1$irf1; irf2_m1 <- r1$irf2
irf1_m2 <- r2$irf1; irf2_m2 <- r2$irf2
irf1_m3 <- r3$irf1; irf2_m3 <- r3$irf2
irf1_m4 <- r4$irf1; irf2_m4 <- r4$irf2


# ============================================================================
# 4. IRF PLOTS
# ============================================================================

plot_irf <- function(irfs, colors, labels, title, ylab) {
  ylim <- range(sapply(irfs, function(x) c(x$ci_low, x$ci_high)), na.rm = TRUE)
  plot(NULL, xlim = c(0, H), ylim = ylim,
       xlab = "Horizon (months)", ylab = ylab, main = title)
  abline(h = 0, lty = 2, col = "black", lwd = 0.8)
  for (i in seq_along(irfs)) {
    irf <- irfs[[i]]
    polygon(c(irf$h, rev(irf$h)), c(irf$ci_low, rev(irf$ci_high)),
            col = adjustcolor(colors[i], alpha.f = 0.12), border = NA)
    lines(irf$h, irf$coef, col = colors[i], lwd = 2)
  }
  legend("topleft", legend = labels, col = colors, lwd = 2, bty = "n", cex = 0.8)
}

irf_colors <- c("steelblue","darkorange","forestgreen","crimson")
irf_labels <- c("M1: HC1, Full","M2: HAC, Full","M3: HAC, Pre-2020","M4: HC1, Pre-2020")

par(mfrow = c(1, 2))
plot_irf(list(irf1_m1, irf1_m2, irf1_m3, irf1_m4), irf_colors, irf_labels,
         "Exp 1: Yield Curve Inversion (Dummy)",
         "Cumulative change in Saving Rate (pp)")
plot_irf(list(irf2_m1, irf2_m2, irf2_m3, irf2_m4), irf_colors, irf_labels,
         "Exp 2: 1pp Decline in YC Spread",
         "Cumulative change in Saving Rate (pp)")
par(mfrow = c(1, 1))

# Model 1 & 4 side-by-side (headline)
par(mfrow = c(2, 2))
for (info in list(
  list(irf1_m1, "Exp1 - M1: HC1, Full Sample",  "steelblue"),
  list(irf1_m4, "Exp1 - M4: HC1, Pre-2020",      "crimson"),
  list(irf2_m1, "Exp2 - M1: HC1, Full Sample",   "steelblue"),
  list(irf2_m4, "Exp2 - M4: HC1, Pre-2020",      "crimson")
)) {
  irf <- info[[1]]; ttl <- info[[2]]; col <- info[[3]]
  ylim <- range(c(irf$ci_low, irf$ci_high), na.rm = TRUE)
  plot(irf$h, irf$coef, type = "l", col = col, lwd = 2,
       ylim = ylim, xlab = "Horizon (months)",
       ylab = "Cumul. change in Saving Rate (pp)", main = ttl)
  polygon(c(irf$h, rev(irf$h)), c(irf$ci_low, rev(irf$ci_high)),
          col = adjustcolor(col, 0.15), border = NA)
  abline(h = 0, lty = 2, lwd = 0.8)
  legend("topleft", legend = c("IRF","90% CI"),
         col = c(col, adjustcolor(col, 0.4)), lwd = c(2,8), bty = "n", cex = 0.8)
}
par(mfrow = c(1, 1))


# ============================================================================
# 5. SIGNIFICANCE TESTS
# ============================================================================

pval_stars <- function(p) {
  ifelse(p < 0.01, "***", ifelse(p < 0.05, "**", ifelse(p < 0.10, "*", "")))
}

make_pval_table <- function(irf_m1, irf_m4, label) {
  df <- data.frame(
    Horizon = irf_m1$h,
    M1_coef = round(irf_m1$coef, 3),
    M1_pval = round(irf_m1$pval, 3),
    M1_sig  = pval_stars(irf_m1$pval),
    M4_coef = round(irf_m4$coef, 3),
    M4_pval = round(irf_m4$pval, 3),
    M4_sig  = pval_stars(irf_m4$pval)
  )
  cat("\n", paste(rep("=", 65), collapse=""), "\n")
  cat("  ", label, "\n")
  cat("  Significance: * p<0.10  ** p<0.05  *** p<0.01\n")
  cat(paste(rep("=", 65), collapse=""), "\n")
  print(df, row.names = FALSE)
  invisible(df)
}

make_pval_table(irf1_m1, irf1_m4, "Experiment 1: Yield Curve Inversion")
make_pval_table(irf2_m1, irf2_m4, "Experiment 2: 1pp YC Spread Decline")

compute_joint <- function(irf, exclude_h0 = TRUE) {
  df   <- if (exclude_h0) irf[irf$h > 0, ] else irf
  z2   <- qnorm(df$pval / 2, lower.tail = FALSE)^2
  chi2 <- sum(z2)
  n_h  <- nrow(df)
  list(chi2 = chi2, df = n_h,
       p    = pchisq(chi2, df = n_h, lower.tail = FALSE))
}

joint_test <- function(irf, label) {
  jt <- compute_joint(irf)
  cat("\n", paste(rep("=", 65), collapse=""), "\n")
  cat("  Joint Test:", label, "\n")
  cat(paste(rep("=", 65), collapse=""), "\n")
  cat("  Horizons tested:", jt$df, "\n")
  cat("  Chi2 statistic: ", round(jt$chi2, 3), "\n")
  cat("  Joint p-value:  ", round(jt$p, 4), pval_stars(jt$p), "\n")
  sub  <- irf[irf$h > 0, ]
  top3 <- head(sub[order(sub$pval), c("h","pval")], 3)
  cat("  Top 3 horizons:\n"); print(top3, row.names = FALSE)
}

joint_test(irf1_m1, "Exp1 - Inversion Dummy, Model 1 (Full Sample)")
joint_test(irf2_m1, "Exp2 - YC Spread, Model 1 (Full Sample)")
joint_test(irf1_m4, "Exp1 - Inversion Dummy, Model 4 (Pre-2020)")
joint_test(irf2_m4, "Exp2 - YC Spread, Model 4 (Pre-2020)")


# ============================================================================
# 6. LATEX TABLE EXPORT
# ============================================================================

fmt_coef <- function(coef, pval, dec = 3) {
  stars <- ifelse(pval < 0.01, "***", ifelse(pval < 0.05, "**",
                                             ifelse(pval < 0.10, "*", "")))
  paste0(formatC(coef, digits = dec, format = "f"), stars)
}

fmt_ci <- function(lo, hi, dec = 3) {
  paste0("[", formatC(lo, digits = dec, format = "f"), ", ",
         formatC(hi, digits = dec, format = "f"), "]")
}

make_irf_table <- function(irf_m1, irf_m4, caption, label, filename) {
  rows <- lapply(which(irf_m1$h > 0), function(i) {
    r1  <- irf_m1[i, ]
    r4  <- irf_m4[i, ]
    sig <- r1$pval < 0.10 | r4$pval < 0.10
    h_str <- if (sig) paste0("\\textbf{", r1$h, "}") else as.character(r1$h)
    paste(h_str,
          fmt_coef(r1$coef, r1$pval),
          fmt_ci(r1$ci_low, r1$ci_high),
          fmt_coef(r4$coef, r4$pval),
          fmt_ci(r4$ci_low, r4$ci_high),
          sep = " & ")
  })
  body <- paste(sapply(rows, function(r) paste0(r, " \\\\")), collapse = "\n")
  latex <- paste0(
    "\\begin{table}[htbp]\n\\centering\\small\n",
    "\\caption{", caption, "}\n\\label{", label, "}\n",
    "\\begin{tabular}{c cccc}\n\\toprule\n",
    " & \\multicolumn{2}{c}{\\textbf{Model 1}} & \\multicolumn{2}{c}{\\textbf{Model 4}} \\\\\n",
    " & \\multicolumn{2}{c}{HC1 SEs, Full Sample} & \\multicolumn{2}{c}{HC1 SEs, Pre-2020} \\\\\n",
    "\\cmidrule(lr){2-3} \\cmidrule(lr){4-5}\n",
    "\\textbf{h} & \\textbf{Coef.} & \\textbf{90\\% CI} & \\textbf{Coef.} & \\textbf{90\\% CI} \\\\\n",
    "\\midrule\n", body, "\n\\midrule\n",
    "\\multicolumn{5}{l}{\\textit{Notes:} Dep. var.: cumulative change in PSAVERT (pp).} \\\\\n",
    "\\multicolumn{5}{l}{h=0 excluded (Cholesky). Controls: contemp. PSAVERT, Inflation,} \\\\\n",
    "\\multicolumn{5}{l}{COVID dummy, NBER recession dummy, trend, ", p, " BIC lags.} \\\\\n",
    "\\multicolumn{5}{l}{* p$<$0.10, ** p$<$0.05, *** p$<$0.01. 90\\% CI in brackets.} \\\\\n",
    "\\bottomrule\n\\end{tabular}\n\\end{table}"
  )
  writeLines(latex, paste0(output_path, filename))
  cat("Saved:", paste0(output_path, filename), "\n")
  invisible(latex)
}

make_irf_table(irf1_m1, irf1_m4,
               "IRF: Effect of Yield Curve Inversion on Personal Saving Rate",
               "tab:irf_exp1", "table_irf_exp1.tex")

make_irf_table(irf2_m1, irf2_m4,
               "IRF: Effect of 1pp Yield Curve Decline on Personal Saving Rate",
               "tab:irf_exp2", "table_irf_exp2.tex")

# --- Summary table ---
jt_m1e1 <- compute_joint(irf1_m1)
jt_m4e1 <- compute_joint(irf1_m4)
jt_m1e2 <- compute_joint(irf2_m1)
jt_m4e2 <- compute_joint(irf2_m4)

latex_summary <- paste0(
  "\\begin{table}[htbp]\n\\centering\\small\n",
  "\\caption{Local Projection Summary Statistics and Joint Significance Tests}\n",
  "\\label{tab:lp_summary}\n",
  "\\begin{tabular}{l cccc}\n\\toprule\n",
  " & \\multicolumn{2}{c}{\\textbf{Experiment 1}} & \\multicolumn{2}{c}{\\textbf{Experiment 2}} \\\\\n",
  " & \\multicolumn{2}{c}{Inversion Dummy} & \\multicolumn{2}{c}{1pp YC Decline} \\\\\n",
  "\\cmidrule(lr){2-3} \\cmidrule(lr){4-5}\n",
  " & \\textbf{Model 1} & \\textbf{Model 4} & \\textbf{Model 1} & \\textbf{Model 4} \\\\\n",
  " & Full Sample & Pre-2020 & Full Sample & Pre-2020 \\\\\n\\midrule\n",
  "Observations & ", nrow(df_lp), " & ", nrow(df_lp_pre2020),
  " & ", nrow(df_lp), " & ", nrow(df_lp_pre2020), " \\\\\n",
  "Horizons (h=1--36) & 36 & 36 & 36 & 36 \\\\\n",
  "BIC-selected lags & ", p, " & ", p, " & ", p, " & ", p, " \\\\\n\\midrule\n",
  "\\textit{Joint significance test:} & & & & \\\\\n",
  "$\\chi^2$ statistic & ", round(jt_m1e1$chi2,3), " & ", round(jt_m4e1$chi2,3),
  " & ", round(jt_m1e2$chi2,3), " & ", round(jt_m4e2$chi2,3), " \\\\\n",
  "Degrees of freedom & ", jt_m1e1$df, " & ", jt_m4e1$df,
  " & ", jt_m1e2$df, " & ", jt_m4e2$df, " \\\\\n",
  "Joint p-value & ", formatC(jt_m1e1$p, digits=3, format="f"),
  " & ", formatC(jt_m4e1$p, digits=3, format="f"), pval_stars(jt_m4e1$p),
  " & ", formatC(jt_m1e2$p, digits=3, format="f"),
  " & ", formatC(jt_m4e2$p, digits=3, format="f"), pval_stars(jt_m4e2$p), " \\\\\n",
  "\\midrule\n",
  "\\multicolumn{5}{l}{\\textit{Notes:} HC1 robust SEs. Model 1: full sample with COVID dummy.} \\\\\n",
  "\\multicolumn{5}{l}{Model 4: pre-2020 sample. Joint test: $\\chi^2=\\sum z_h^2$, $\\chi^2(36)$ under H$_0$.} \\\\\n",
  "\\multicolumn{5}{l}{h=0 excluded from joint test (Cholesky impact artifact).} \\\\\n",
  "\\bottomrule\n\\end{tabular}\n\\end{table}"
)

writeLines(latex_summary, paste0(output_path, "table_lp_summary.tex"))
cat("Saved:", paste0(output_path, "table_lp_summary.tex"), "\n")
cat("\n--- LaTeX Preview ---\n")
cat(latex_summary, "\n")
