suppressPackageStartupMessages({
  library(meta)
  library(ggplot2)
  library(ggrepel)
})

options(stringsAsFactors = FALSE, warn = -1)

# ============================================================
# [0] 基础路径与全局参数配置
# ============================================================
BASE      <- "C:/Users/TianhaoLi/Desktop/TreME_abundance"
META_CSV  <- file.path(BASE, "sample_metadata_0414_precise_corrected.csv")
ABUND_CSV <- file.path(BASE, "sample_celltype_abundance.csv")

# 输出目录：修改为包含 MCC 的独立结果文件夹
OUT_DIR   <- file.path(BASE, "results_total_immune_wilcox_withMCC——")

MIN_SAMPLES           <- 2
MIN_NONZERO           <- 0.03
FDR_THRESH            <- 0.10
P_THRESH              <- 0.05
WITHIN_FDR_THRESH     <- 0.10
EXCLUDE_CTS           <- c("ILC_Unknown")
T_CELL_PREFIXES       <- c("CD4_", "CD8_")
T_CELL_ONLY_THRESHOLD <- 0.999
TIMEPOINTS            <- c("Pre", "Post")

IMMUNE_PREFIXES <- c("CD4_", "CD8_", "ILC_", "Bplasma_", "Myeloid_")
is_immune <- function(ct) any(startsWith(ct, IMMUNE_PREFIXES))

# ============================================================
# [1] 数据加载与预处理
# ============================================================
cat("[1] Loading and cleaning dataset...\n")
meta_df_all  <- read.csv(META_CSV, check.names = FALSE)

# 【修改】保留 MCC 癌症类型，不再执行 MCC 过滤

abund_df_all <- read.csv(ABUND_CSV, row.names = 1, check.names = FALSE)
if ("majorCluster" %in% rownames(abund_df_all)) {
  abund_df_all <- abund_df_all[setdiff(rownames(abund_df_all), "majorCluster"), , drop = FALSE]
}

abund_df_all[] <- lapply(abund_df_all, function(x) as.numeric(as.character(x)))
abund_df_all[is.na(abund_df_all)] <- 0

# 处理基因/细胞类型命名重叠
if ("Fib_Fib-SOX6+" %in% colnames(abund_df_all)) {
  if ("Fib_SOX6+Fib" %in% colnames(abund_df_all)) {
    abund_df_all[["Fib_SOX6+Fib"]] <- abund_df_all[["Fib_SOX6+Fib"]] + abund_df_all[["Fib_Fib-SOX6+"]]
    abund_df_all <- abund_df_all[, setdiff(colnames(abund_df_all), "Fib_Fib-SOX6+"), drop = FALSE]
  } else {
    colnames(abund_df_all)[colnames(abund_df_all) == "Fib_Fib-SOX6+"] <- "Fib_SOX6+Fib"
  }
}

# 过滤仅保留免疫细胞
immune_cols  <- colnames(abund_df_all)[sapply(colnames(abund_df_all), is_immune)]
abund_df_all <- abund_df_all[, immune_cols, drop = FALSE]

# 映射响应（R/NR）
RESP_MAP <- c(
  CR = "R", `CR*` = "R", PR = "R", R = 'EXCLUDE', pCR = "EXCLUDE", MPR = "EXCLUDE", OR = "EXCLUDE",
  SD = "NR", PD = "NR", `NR/SD` = "NR", NR = "EXCLUDE", `non-MPR` = "EXCLUDE", `non-pCR` = "EXCLUDE",
  `Post-ICI (resistant)` = "EXCLUDE", EXCLUDE = "EXCLUDE", unknown = "EXCLUDE", unknowm = "EXCLUDE", Untreated = "EXCLUDE"
)
meta_df_all$Response_Binary <- unname(RESP_MAP[meta_df_all$Response_Unified])
meta_df_all$Response_Binary[is.na(meta_df_all$Response_Binary)] <- "EXCLUDE"

prepare_dataset <- function(meta_all, abund_all, tp) {
  m <- meta_all[meta_all$Timepoint_Standardized == tp & meta_all$Response_Binary != "EXCLUDE", , drop = FALSE]
  keep_samps <- intersect(as.character(m$Sample_new), rownames(abund_all))
  m <- m[as.character(m$Sample_new) %in% keep_samps, , drop = FALSE]
  a <- abund_all[keep_samps, , drop = FALSE]
  
  # 剔除仅有 T 细胞富集的离群样本
  t_cols <- colnames(a)[sapply(colnames(a), function(ct) any(startsWith(ct, T_CELL_PREFIXES)))]
  if (length(t_cols) > 0 && nrow(a) > 0) {
    sample_total <- rowSums(a, na.rm = TRUE)
    sample_tcell <- rowSums(a[, t_cols, drop = FALSE], na.rm = TRUE)
    tcell_ratio  <- ifelse(sample_total > 0, sample_tcell / sample_total, NA_real_)
    remove_samps <- names(tcell_ratio)[is.finite(tcell_ratio) & tcell_ratio >= T_CELL_ONLY_THRESHOLD]
    if (length(remove_samps) > 0) {
      a <- a[setdiff(rownames(a), remove_samps), , drop = FALSE]
      m <- m[!(as.character(m$Sample_new) %in% remove_samps), , drop = FALSE]
    }
  }
  
  valid_cts <- character(0)
  if (ncol(a) > 0 && nrow(a) > 0) {
    nonzero_ratio <- colMeans(a > 0, na.rm = TRUE)
    valid_cts     <- setdiff(names(nonzero_ratio[nonzero_ratio > MIN_NONZERO]), EXCLUDE_CTS)
  }
  list(meta = m, abund = a, valid_cts = valid_cts)
}

datasets <- lapply(setNames(TIMEPOINTS, TIMEPOINTS), function(tp) prepare_dataset(meta_df_all, abund_df_all, tp))

# ============================================================
# [2] Meta 分析引擎 (Total Immune Denominator + Continuous SMD)
# ============================================================
run_meta_total_wilcox <- function(ds) {
  meta_sub  <- ds$meta
  abund_sub <- ds$abund
  cts       <- ds$valid_cts
  out       <- list()
  cancer_groups <- split(meta_sub, meta_sub$Cancer_Type)
  
  # 分母：总免疫细胞绝对数
  total_immune_denom <- rowSums(abund_sub, na.rm = TRUE)
  
  for (ct in cts) {
    if (!(ct %in% colnames(abund_sub))) next
    
    # 计算目标亚群在总免疫细胞中的相对比例 (Fraction)
    pct_vec <- abund_sub[, ct] / total_immune_denom
    
    rows <- list()
    for (cancer in names(cancer_groups)) {
      c_meta   <- cancer_groups[[cancer]]
      r_samps  <- intersect(as.character(c_meta$Sample_new[c_meta$Response_Binary == "R"]), rownames(abund_sub))
      nr_samps <- intersect(as.character(c_meta$Sample_new[c_meta$Response_Binary == "NR"]), rownames(abund_sub))
      if (length(r_samps) < MIN_SAMPLES || length(nr_samps) < MIN_SAMPLES) next
      
      r_vals  <- pct_vec[r_samps]
      nr_vals <- pct_vec[nr_samps]
      if (length(r_vals) < 2 || length(nr_vals) < 2) next
      
      wx_p <- tryCatch(wilcox.test(r_vals, nr_vals, exact = FALSE)$p.value, error = function(e) NA_real_)
      
      rows[[length(rows) + 1]] <- data.frame(
        Cancer = cancer, 
        n_e = length(r_vals),  mean_e = mean(r_vals, na.rm = TRUE),  sd_e = sd(r_vals, na.rm = TRUE),
        n_c = length(nr_vals), mean_c = mean(nr_vals, na.rm = TRUE), sd_c = sd(nr_vals, na.rm = TRUE),
        wilcox_p = wx_p, stringsAsFactors = FALSE
      )
    }
    
    if (length(rows) < 2) next
    df_ct <- do.call(rbind, rows)
    
    # 拟合连续型标准化均值差 (SMD) 随机效应 Meta 模型 (REML)
    meta_fit <- tryCatch(
      metacont(n.e = df_ct$n_e, mean.e = df_ct$mean_e, sd.e = df_ct$sd_e,
               n.c = df_ct$n_c, mean.c = df_ct$mean_c, sd.c = df_ct$sd_c,
               studlab = df_ct$Cancer, sm = "SMD", method.tau = "REML", random = TRUE, fixed = FALSE),
      error = function(e) NULL
    )
    
    if (is.null(meta_fit)) next
    
    out[[ct]] <- list(
      df_ct = df_ct, meta_fit = meta_fit, 
      pooled_p = meta_fit$pval.random, lineage_prefix = "TotalImmune"
    )
  }
  
  # 全局 FDR (Benjamini-Hochberg)
  if (length(out) > 0) {
    all_p <- sapply(names(out), function(ct) out[[ct]]$pooled_p)
    all_fdr <- p.adjust(all_p, method = "BH")
    for (ct in names(out)) out[[ct]]$fdr <- all_fdr[ct]
  }
  
  out
}

# ============================================================
# [3] 单细胞类型森林图绘制函数 (Continuous SMD 专用)
# ============================================================
draw_individual_forest_smd <- function(ct, obj, tp, out_file) {
  df_ct    <- obj$df_ct
  meta_fit <- obj$meta_fit
  
  # 提取每个研究的 SMD 及 95% CI
  df_ct$TE <- meta_fit$TE
  df_ct$lo <- meta_fit$lower
  df_ct$hi <- meta_fit$upper
  
  # 提取 Pooled 结果
  pooled_g  <- meta_fit$TE.random
  pooled_lo <- meta_fit$lower.random
  pooled_hi <- meta_fit$upper.random
  pooled_p  <- meta_fit$pval.random
  i2_val    <- round(meta_fit$I2 * 100, 0)
  
  p_sig    <- is.finite(pooled_p) && pooled_p < P_THRESH
  col_pool <- if (p_sig) "#d62728" else "#555555"
  
  p_label    <- ifelse(pooled_p < 0.001, "p<0.001", sprintf("p=%.3f", pooled_p))
  reml_label <- sprintf("Continuous SMD pooled (I\u00b2=%d%%, %s)", i2_val, p_label)
  
  n_studies <- nrow(df_ct)
  y_studies <- rev(seq_len(n_studies))
  y_pooled  <- -0.6
  
  # 计算图内 (Within-plot) 的 Wilcoxon FDR
  df_ct$wilcox_fdr <- p.adjust(df_ct$wilcox_p, method = "BH")
  df_ct$col <- ifelse(!is.na(df_ct$wilcox_fdr) & df_ct$wilcox_fdr < WITHIN_FDR_THRESH, "#d62728", "#2166ac")
  
  all_x   <- c(df_ct$lo, df_ct$hi, pooled_lo, pooled_hi)
  all_x   <- all_x[is.finite(all_x)]
  x_range <- range(all_x)
  x_pad   <- max((x_range[2] - x_range[1]) * 0.15, 0.05)
  xlim    <- c(x_range[1] - x_pad, x_range[2] + x_pad)
  
  fmt_ci_q <- function(g, lo, hi, q) {
    base <- sprintf("%.2f [%.2f, %.2f]", g, lo, hi)
    if (is.na(q) || !is.finite(q)) return(base)
    if (q < 0.001) return(paste0(base, "  q<0.001"))
    paste0(base, sprintf("  q=%.3f", q))
  }
  df_ct$label <- mapply(fmt_ci_q, df_ct$TE, df_ct$lo, df_ct$hi, df_ct$wilcox_fdr)
  pool_label  <- sprintf("%.2f [%.2f, %.2f]", pooled_g, pooled_lo, pooled_hi)
  
  fig_h <- max(3.5, (n_studies + 3) * 0.42 + 1.5)
  png(out_file, width = 1500, height = ceiling(fig_h * 155), res = 160)
  par(mar = c(4.5, 9, 3, 11), mgp = c(2, 0.6, 0), xpd = FALSE)
  
  y_range <- c(y_pooled - 0.7, n_studies + 0.6)
  
  plot(NA, xlim = xlim, ylim = y_range,
       xlab = "Standardized Mean Difference (SMD) in Total Immune (R vs NR)", ylab = "", 
       yaxt = "n", bty = "n",
       main = sprintf("%s (Denominator: Total Immune)  |  %s", ct, tp), 
       cex.main = 1.05, font.main = 2, cex.lab = 0.85)
  
  abline(v = 0, lty = 2, col = "grey55", lwd = 1)
  
  # 背景灰色斑马纹
  for (i in seq_len(n_studies)) {
    yi <- y_studies[i]
    if (i %% 2 == 0) rect(xlim[1], yi - 0.45, xlim[2], yi + 0.45, col = "#f2f2f2", border = NA)
  }
  abline(h = 0, lty = 1, col = "grey75", lwd = 0.7)
  
  # 各研究 SMD 误差棒与数据点
  for (i in seq_len(n_studies)) {
    yi  <- y_studies[i]
    col <- df_ct$col[i]
    if(!is.finite(df_ct$lo[i]) || !is.finite(df_ct$hi[i])) next
    segments(df_ct$lo[i], yi, df_ct$hi[i], yi, col = col, lwd = 1.7)
    points(df_ct$TE[i], yi, pch = 15, cex = 0.9, col = col)
  }
  
  # Pooled 菱形
  h_dia <- 0.30
  polygon(x = c(pooled_lo, pooled_g, pooled_hi, pooled_g),
          y = c(y_pooled,  y_pooled + h_dia, y_pooled, y_pooled - h_dia),
          col = col_pool, border = col_pool)
  
  axis(2, at = y_studies, labels = df_ct$Cancer, las = 2, cex.axis = 0.78, tick = FALSE, hadj = 1)
  text(xlim[1], y_pooled, labels = reml_label, adj = c(0, 0.5), cex = 0.68, col = col_pool, font = if (p_sig) 2L else 1L, xpd = TRUE)
  
  mtext("SMD [95% CI]   q (within-plot FDR)", side = 3, line = 0.3, at = xlim[2] + x_pad * 0.2, adj = 0, cex = 0.68, font = 2, xpd = TRUE)
  for (i in seq_len(n_studies)) {
    mtext(df_ct$label[i], side = 4, at = y_studies[i], las = 2, cex = 0.62, col = df_ct$col[i], font = if (df_ct$col[i] == "#d62728") 2L else 1L, xpd = TRUE, line = 0.4)
  }
  mtext(pool_label, side = 4, at = y_pooled, las = 2, cex = 0.65, col = col_pool, font = 2L, xpd = TRUE, line = 0.4)
  
  legend("bottomright", legend = c(sprintf("Cancer FDR<%.2f (BH within plot)", WITHIN_FDR_THRESH), "Not significant", if (p_sig) sprintf("Pooled p<%.2f \u2713", P_THRESH) else sprintf("Pooled p\u2265%.2f", P_THRESH)),
         col = c("#d62728", "#2166ac", col_pool), pch = c(15, 15, 18), pt.cex = c(1.0, 1.0, 1.3), bty = "n", cex = 0.68)
  dev.off()
}

# ============================================================
# [4] 汇总森林图绘制函数 (对应图片样式)
# ============================================================
draw_overview_self_sorted <- function(tp_res, tp_label, out_path) {
  if (length(tp_res) == 0) return(invisible(NULL))
  
  rows <- lapply(names(tp_res), function(ct) {
    obj <- tp_res[[ct]]
    mf  <- obj$meta_fit
    data.frame(
      Cell_Type = ct, 
      TE        = mf$TE.random, 
      lo        = mf$lower.random, 
      hi        = mf$upper.random, 
      Meta_P    = mf$pval.random, 
      N_cancers = nrow(obj$df_ct), 
      stringsAsFactors = FALSE
    )
  })
  
  df <- do.call(rbind, rows)
  df <- df[is.finite(df$TE) & is.finite(df$lo) & is.finite(df$hi), , drop = FALSE]
  if (nrow(df) == 0) return(invisible(NULL))
  
  # 按效应量降序排列
  df <- df[order(-df$TE), , drop = FALSE]
  n_rows <- nrow(df)
  df$y   <- rev(seq_len(n_rows))
  
  # 判断显著性与颜色映射 (p < 0.05 为红字高亮)
  df$sig <- !is.na(df$Meta_P) & is.finite(df$Meta_P) & df$Meta_P < P_THRESH
  df$col <- ifelse(df$sig, "#d62728", "#333333")
  
  all_x <- c(df$lo, df$hi); all_x <- all_x[is.finite(all_x)]
  xrng  <- range(all_x)
  xpad  <- max((xrng[2] - xrng[1]) * 0.08, 0.05)
  xlim  <- c(xrng[1] - xpad, xrng[2] + xpad)
  
  fmt_p <- function(p) {
    if (is.na(p) || !is.finite(p)) return("n.s.")
    if (p < 1e-4) return(sprintf("%.2e", p))
    if (p < 0.001) return("p<0.001")
    sprintf("%.4f", p)
  }
  
  df$p_label  <- sapply(df$Meta_P, fmt_p)
  df$te_label <- sprintf("%.2f [%.2f, %.2f]", df$TE, df$lo, df$hi)
  
  fig_h <- max(8, n_rows * 0.30 + 2.5)
  png(out_path, width = 2400, height = ceiling(fig_h * 120), res = 130)
  par(mar = c(4, 14, 4, 14), mgp = c(2.5, 0.7, 0), xpd = FALSE)
  
  plot(NA, xlim = xlim, ylim = c(0.3, n_rows + 0.7), 
       xlab = "Pooled Log Odds Ratio / Effect Size (R vs NR)", ylab = "", yaxt = "n", bty = "n", 
       main = sprintf("Cross-cancer meta (self-sig, sorted by Log OR): %s [total_immune_wilcox_withMCC]", tp_label), 
       cex.main = 1.15, font.main = 2, cex.lab = 0.90)
  
  abline(v = 0, lty = 2, col = "grey55", lwd = 1.2)
  
  # 背景灰色斑马纹
  for (k in seq_len(n_rows)) {
    yi <- df$y[k]
    if (k %% 2 == 0) rect(xlim[1], yi - 0.45, xlim[2], yi + 0.45, col = "#f4f4f4", border = NA)
  }
  
  # CI 线段与中点
  for (i in seq_len(n_rows)) {
    yi  <- df$y[i]; col <- df$col[i]
    segments(df$lo[i], yi, df$hi[i], yi, col = col, lwd = 1.6)
    segments(df$lo[i], yi - 0.15, df$lo[i], yi + 0.15, col = col, lwd = 1.2)
    segments(df$hi[i], yi - 0.15, df$hi[i], yi + 0.15, col = col, lwd = 1.2)
    points(df$TE[i], yi, pch = if(df$sig[i]) 18 else 15, cex = if(df$sig[i]) 1.5 else 0.85, col = col)
  }
  
  # 左侧 Y 轴细胞亚群名称
  for (i in seq_len(n_rows)) {
    mtext(df$Cell_Type[i], side = 2, at = df$y[i], las = 2, cex = 0.72, xpd = TRUE, 
          col = df$col[i], font = if (df$sig[i]) 2L else 1L, line = 0.5)
  }
  
  # 右侧数据文本标注
  x_right <- xlim[2]
  mtext("Log OR [95% CI]", side = 3, at = x_right, adj = 0, cex = 0.75, font = 2, xpd = TRUE, line = 0.8)
  mtext("p-meta",          side = 3, at = x_right, adj = -2.2, cex = 0.75, font = 2, xpd = TRUE, line = 0.8)
  
  for (i in seq_len(n_rows)) {
    yi  <- df$y[i]; col <- df$col[i]; fn <- if (df$sig[i]) 2L else 1L
    mtext(df$te_label[i], side = 4, at = yi, las = 2, cex = 0.65, col = col, font = fn, xpd = TRUE, line = 0.5)
    mtext(df$p_label[i],  side = 4, at = yi, las = 2, cex = 0.65, col = col, font = fn, xpd = TRUE, line = 7.5)
  }
  
  legend("bottomleft", 
         legend = c(sprintf("Pooled p < %.2f (self)", P_THRESH), "Not significant"), 
         col = c("#d62728", "#333333"), pch = c(18, 15), pt.cex = c(1.5, 0.85), bty = "n", cex = 0.82)
  
  dev.off()
}

# ============================================================
# [4b] 火山-森林图 (Volcano Forest)
# ============================================================
draw_volcano_forest <- function(tp_res, tp_label, out_path) {
  if (length(tp_res) == 0) return(invisible(NULL))

  rows <- lapply(names(tp_res), function(ct) {
    obj <- tp_res[[ct]]
    mf  <- obj$meta_fit
    # 方向一致性：各癌种 SMD 符号的净投票 (正 - 负)
    te_studies <- mf$TE
    te_studies <- te_studies[is.finite(te_studies)]
    sign_consistency <- if (length(te_studies) > 0) sum(sign(te_studies)) else NA_real_

    data.frame(
      Cell_Type = ct,
      TE        = mf$TE.random,
      lo        = mf$lower.random,
      hi        = mf$upper.random,
      Meta_P    = mf$pval.random,
      SignCons  = sign_consistency,
      stringsAsFactors = FALSE
    )
  })

  df <- do.call(rbind, rows)
  df <- df[is.finite(df$TE) & is.finite(df$Meta_P), , drop = FALSE]
  if (nrow(df) == 0) return(invisible(NULL))

  df$neglog10p <- -log10(df$Meta_P)
  df$sig       <- is.finite(df$Meta_P) & df$Meta_P < P_THRESH

  # 颜色范围对称
  cons_max <- max(abs(df$SignCons), na.rm = TRUE)
  if (!is.finite(cons_max) || cons_max == 0) cons_max <- 1

  df_sig <- df[df$sig, , drop = FALSE]
  df_ns  <- df[!df$sig, , drop = FALSE]

  p <- ggplot(df, aes(x = TE, y = neglog10p)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.5) +
    # 非显著：极淡的 CI 误差棒 + 灰点（弱化，避免抢视觉）
    geom_errorbarh(data = df_ns, aes(xmin = lo, xmax = hi),
                   height = 0, linewidth = 0.3, color = "grey80", alpha = 0.35) +
    geom_point(data = df_ns, color = "grey65", alpha = 0.55, size = 1.6)

  # 显著：CI 误差棒 + 彩色点(白描边) + 标签
  if (nrow(df_sig) > 0) {
    p <- p +
      geom_errorbarh(data = df_sig, aes(xmin = lo, xmax = hi, color = SignCons),
                     height = 0, linewidth = 0.8, alpha = 0.85) +
      # 白色描边：底层稍大白点
      geom_point(data = df_sig, color = "white", size = 4.2) +
      geom_point(data = df_sig, aes(color = SignCons), size = 3) +
      ggrepel::geom_text_repel(
        data = df_sig, aes(label = Cell_Type, color = SignCons),
        size = 3.1, fontface = "bold", show.legend = FALSE,
        seed = 1, max.overlaps = Inf,
        box.padding = 0.6, point.padding = 0.4,
        min.segment.length = 0, segment.size = 0.3,
        segment.color = "grey55", segment.alpha = 0.7,
        force = 3, force_pull = 0.5
      )
  }

  p <- p +
    scale_color_gradientn(
      colours = c("#313695", "#4575B4", "#91BFDB", "#E0F3F8",
                  "#FFFFFF",
                  "#FEE090", "#FC8D59", "#D73027", "#A50026"),
      limits = c(-cons_max, cons_max),
      name = "Sign\nconsistency"
    ) +
    labs(
      x = "Pooled Effect Size (SMD, R vs NR)",
      y = expression(-log[10]("Summary p-value")),
      title = sprintf("Volcano Forest Meta-Analysis (%s)", tp_label)
    ) +
    theme_bw(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey93"),
      axis.title = element_text(face = "bold"),
      legend.position = "left",
      legend.title = element_text(size = 10, face = "bold"),
      legend.key.height = unit(0.9, "cm")
    )

  ggsave(out_path, plot = p, width = 9.5, height = 7.5, dpi = 300)
}

# ============================================================
# [5] 主流程执行
# ============================================================
cat("\n[3] Running total_immune_wilcox meta-analysis (with MCC)...\n")

for (tp in TIMEPOINTS) {
  cat(sprintf("\nProcessing timepoint: %s\n", tp))
  tp_out_dir <- file.path(OUT_DIR, tp)
  dir.create(tp_out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # 计算 Meta 结果
  res <- run_meta_total_wilcox(datasets[[tp]])
  
  # 1. 绘制全局汇总图
  cat("  -> Drawing overview forest plot...\n")
  draw_overview_self_sorted(res, tp, file.path(tp_out_dir, "overview_self_sig_sorted.png"))

  # 1b. 绘制火山-森林图
  cat("  -> Drawing volcano forest plot...\n")
  draw_volcano_forest(res, tp, file.path(tp_out_dir, "volcano_forest.png"))
  
  # 2. 循环绘制每个细胞类型的单独森林图
  cat(sprintf("  -> Drawing %d individual cell-type forest plots...\n", length(res)))
  for (ct in names(res)) {
    ct_filename <- paste0(gsub("[^A-Za-z0-9_]", "_", ct), ".png")
    out_file    <- file.path(tp_out_dir, ct_filename)
    draw_individual_forest_smd(ct, res[[ct]], tp, out_file)
  }
  
  # 3. 导出 CSV 汇总表格
  if (length(res) > 0) {
    sum_df <- do.call(rbind, lapply(names(res), function(ct) {
      mf <- res[[ct]]$meta_fit
      data.frame(
        Cell_Type  = ct,
        EffectSize = mf$TE.random,
        SE         = mf$seTE.random,
        P_val      = mf$pval.random,
        FDR_qval   = res[[ct]]$fdr,
        N_Cancers  = nrow(res[[ct]]$df_ct)
      )
    }))
    write.csv(sum_df, file.path(tp_out_dir, "summary_metrics.csv"), row.names = FALSE)
  }
}

cat("\n============================================================\n")
cat(" Analysis complete!\n")
cat(" Results saved in directory:\n ", OUT_DIR, "\n")
cat("============================================================\n")