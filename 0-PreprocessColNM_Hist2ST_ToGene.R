suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(purrr)
})

# === 路径 ===
base_dir  <- "/Users/a1-6/R/du_1130/hist2gene_20250417" # hist2st_20250506
test_dir  <- file.path(base_dir, "internal", "test")
spot_dir  <- file.path(base_dir, "spot_gene_pairs")
out_dir   <- file.path(test_dir, "analysis")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# === 数字与样本映射 ===
idx_to_case <- c(
  "0" = "A121573",
  "1" = "A798015",
  "2" = "A120838",
  "3" = "A938797",
  "4" = "A416371",
  "5" = "A595688"
)

# === 在 spot_gene_pairs 中按样本 ID + Rep 自动找文件（前缀不限） ===
find_rep_file <- function(case_id, rep = c("Rep1","Rep2"), dir = spot_dir) {
  rep <- match.arg(rep)
  pat <- paste0("_", case_id, "_", rep, "\\.csv$")
  hits <- list.files(dir, pattern = pat, full.names = TRUE)
  if (length(hits) == 0) stop("未找到文件：pattern=", pat, " in ", dir)
  if (length(hits) > 1)  stop("匹配到多个文件（期望唯一）：", paste(basename(hits), collapse = ", "))
  hits[[1]]
}

# === 兼容 gt/pre 两种命名：无下划线/有下划线 ===
# kind: "gt" 或 "pre"
find_gtpre_file <- function(kind = c("gt","pre"), idx, dir = test_dir) {
  kind <- match.arg(kind)
  prefix <- if (kind == "gt") "gt_with_gene" else "pre_with_gene"
  # 用单一 regex 同时接住两种：0 与 _0
  pat <- paste0("^", prefix, "_?", idx, "\\.csv$")
  hits <- list.files(dir, pattern = pat, full.names = TRUE)
  if (length(hits) == 0) stop("未找到文件：", prefix, "{_?}", idx, ".csv in ", dir)
  if (length(hits) > 1)  stop("匹配到多个文件（期望唯一）：", paste(basename(hits), collapse = ", "))
  hits[[1]]
}

# === 解析 Rep 表：保留首列；处理第2列 fn；拆第3列 gene 为数值多列 ===
# 输入列按位置：1=行数(保留)、2=fn、3=gene
parse_rep_csv <- function(rep_csv_path) {
  df_raw <- readr::read_csv(rep_csv_path, show_col_types = FALSE)
  if (ncol(df_raw) < 3) stop("文件列数不足 3：", rep_csv_path)
  
  df <- df_raw
  colnames(df)[2:3] <- c("fn", "gene")
  
  file_prefix <- tools::file_path_sans_ext(basename(rep_csv_path))
  extracted   <- str_extract(df$fn, "^[^_]+")
  extracted   <- str_remove(extracted, "\\.png$")
  df$fn       <- paste0(file_prefix, "_", extracted)
  
  gene_list_chr <- str_remove_all(df$gene, "\\[|\\]")
  gene_list     <- str_split(gene_list_chr, ",\\s*")
  gene_num      <- lapply(gene_list, as.numeric)
  
  lens <- lengths(gene_num)
  if (length(unique(lens)) != 1) {
    stop("gene 数组长度不一致：", rep_csv_path)
  }
  p <- unique(lens)
  
  mat <- do.call(rbind, gene_num)
  mat_df <- as.data.frame(mat, check.names = FALSE)
  colnames(mat_df) <- paste0("g", seq_len(p))
  
  # 返回：[1]=行数(原始), [2]=fn(已改), [3..]=g1..gp(数值)
  bind_cols(df[1], df["fn"], mat_df)
}

# === 浮点比较（带容差） ===
vec_equal_tol <- function(x, y, tol = 1e-8) {
  isTRUE(all.equal(as.numeric(x), as.numeric(y), tolerance = tol))
}

# === 主流程 ===
for (i in names(idx_to_case)) {
  case_id <- idx_to_case[[i]]
  message("处理索引 ", i, "（样本：", case_id, "）...")
  
  rep1_path <- find_rep_file(case_id, "Rep1")
  rep2_path <- find_rep_file(case_id, "Rep2")
  gt_path   <- find_gtpre_file("gt",  idx = i)
  pre_path  <- find_gtpre_file("pre", idx = i)
  
  rep1_df <- parse_rep_csv(rep1_path)  # 列：rownum, fn, g1..gp
  rep2_df <- parse_rep_csv(rep2_path)
  rep_all <- bind_rows(rep1_df, rep2_df)
  
  labels   <- rep_all[[2]]                        # 第2列 = fn（目标列名）
  spot_mat <- as.matrix(rep_all[, -(1:2), drop=FALSE])  # 行=spot，列=gene 序
  
  gt  <- readr::read_csv(gt_path,  show_col_types = FALSE)
  pre <- readr::read_csv(pre_path, show_col_types = FALSE)
  if (ncol(gt)  < 3) stop("GT 列数不足 3：",  gt_path)
  if (ncol(pre) < 3) stop("PRE 列数不足 3：", pre_path)
  
  # gt: [1]=行数(保留), [2]=gene, [3..]=各 spot
  gt_mat <- as.matrix(gt[, -(1:2), drop=FALSE])   # 行=gene，列=spot
  
  # 维度一致性
  if (nrow(spot_mat) != ncol(gt_mat))
    stop("spot 数与 GT 列数不一致：spot=", nrow(spot_mat), " vs GT列=", ncol(gt_mat))
  if (ncol(spot_mat) != nrow(gt_mat))
    stop("gene 维度不一致：Rep gene数=", ncol(spot_mat), " vs GT gene行数=", nrow(gt_mat))
  
  # 逐列比对（Rep 第k行 == GT 第k列）
  checks <- map_lgl(seq_len(nrow(spot_mat)), function(k) {
    vec_equal_tol(spot_mat[k, ], gt_mat[, k], tol = 1e-8)
  })
  if (!all(checks)) {
    bad <- which(!checks)[1]
    stop("数值校验失败（容差1e-8）。首个不一致 spot：", labels[bad],
         "；索引=", i, "；样本=", case_id)
  }
  message("数值校验通过：", basename(gt_path), " 对齐 ", basename(rep1_path), "/", basename(rep2_path))
  
  # 重命名 GT / PRE 的第3..末列为 labels（前两列：行数、gene 保持）
  new_gt  <- gt
  new_pre <- pre
  if ((ncol(new_gt) - 2) != length(labels))
    stop("重命名数量不匹配：GT样本列数=", ncol(new_gt) - 2, " vs labels=", length(labels))
  if ((ncol(new_pre) - 2) != length(labels))
    stop("重命名数量不匹配：PRE样本列数=", ncol(new_pre) - 2, " vs labels=", length(labels))
  
  colnames(new_gt)[3:ncol(new_gt)]   <- labels
  colnames(new_pre)[3:ncol(new_pre)] <- labels
  
  # 输出
  gt_out  <- file.path(out_dir, paste0("test_", i, "_gt_with_gene.csv"))
  pre_out <- file.path(out_dir, paste0("test_", i, "_pre_with_gene.csv"))
  readr::write_csv(new_gt,  gt_out)
  readr::write_csv(new_pre, pre_out)
  
  message("已写出：", gt_out)
  message("已写出：", pre_out)
}
message("全部处理完成。")


############# 二、针对outputs_external中的csv的列名处理，仅需运行一次 #############
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(purrr)
})

# ========== 路径 ==========
base_dir <- "/Users/a1-6/R/du_1130/hist2st_20250506"
spot_dir <- file.path(base_dir, "spot_gene_pairs")

# 数据集目录 ↔ spot 文件映射（左=数据集子目录名，右=spot文件名）
ds_map <- c(
  "18CRC" = "18crc.csv",
  "34CRC" = "34crc.csv",
  "human" = "human.csv",
  "rep1"  = "rep1.csv",
  "rep2"  = "rep2.csv"
)

# ========== 小工具 ==========
# 安全读取 csv，返回 NULL 表示失败；保留列，不改名
read_csv_safely <- function(path) {
  if (!file.exists(path)) {
    message("缺失文件：", path)
    return(NULL)
  }
  tryCatch(readr::read_csv(path, show_col_types = FALSE),
           error = function(e) { message("读取失败：", path, "；", e$message); NULL })
}

# 从数据集目录中查找 gt/pre 源文件（优先无编号文件，其次匹配唯一编号）
# 会返回唯一匹配，否则打印并返回 NULL
find_gt_pre_files <- function(ds_dir) {
  all_csv <- list.files(ds_dir, pattern = "\\.csv$", full.names = TRUE)
  if (length(all_csv) == 0) {
    message("目录下未找到任何 CSV：", ds_dir)
    return(list(gt = NULL, pre = NULL))
  }
  # 先找无编号标准名
  gt_std  <- file.path(ds_dir, "gt_with_gene.csv")
  pre_std <- file.path(ds_dir, "pre_with_gene.csv")
  
  gt_path <- if (file.exists(gt_std)) gt_std else {
    # 退而求其次，找 gt_with_gene*.csv（含编号），需唯一
    hits <- grep("^gt_with_gene(_?\\d+)?\\.csv$", basename(all_csv), value = TRUE)
    if (length(hits) == 0) {
      message("未找到 gt_with_gene*.csv：", ds_dir)
      NULL
    } else if (length(hits) > 1 && !("gt_with_gene2.csv" %in% hits)) {
      # 针对 18CRC 你特别提到是 2，这里做一个优先选择；否则报多匹配
      message("匹配到多个 gt_with_gene*.csv（期望唯一）：", paste(hits, collapse = ", "), " @ ", ds_dir)
      NULL
    } else {
      # 若包含 gt_with_gene2.csv，优先取它（兼容你给的 18CRC 情况）
      chosen <- if ("gt_with_gene2.csv" %in% hits) "gt_with_gene2.csv" else hits[[1]]
      file.path(ds_dir, chosen)
    }
  }
  
  pre_path <- if (file.exists(pre_std)) pre_std else {
    hits <- grep("^pre_with_gene(_?\\d+)?\\.csv$", basename(all_csv), value = TRUE)
    if (length(hits) == 0) {
      message("未找到 pre_with_gene*.csv：", ds_dir)
      NULL
    } else if (length(hits) > 1 && !("pre_with_gene2.csv" %in% hits)) {
      message("匹配到多个 pre_with_gene*.csv（期望唯一）：", paste(hits, collapse = ", "), " @ ", ds_dir)
      NULL
    } else {
      chosen <- if ("pre_with_gene2.csv" %in% hits) "pre_with_gene2.csv" else hits[[1]]
      file.path(ds_dir, chosen)
    }
  }
  
  list(gt = gt_path, pre = pre_path)
}

# 解析 spot csv：
# - 保留第1列（原始行数）；
# - 第2列标准化为 fn，并提取 "^[^_]+"，去掉 .png；
# - 第3列标准化为 gene：从 "[x, y, ...]" 转为数值多列 g1..gp；
# 返回列序： [1]=原始首列, [2]=fn(已抽取), [3..]=g1..gp(数值)
parse_spot_csv <- function(spot_csv_path) {
  df_raw <- read_csv_safely(spot_csv_path)
  if (is.null(df_raw)) return(NULL)
  if (ncol(df_raw) < 3) {
    message("spot 文件列数不足 3：", spot_csv_path)
    return(NULL)
  }
  
  df <- df_raw
  colnames(df)[2:3] <- c("fn", "gene")
  
  # 仅提取条码部分（下划线前），并去掉 .png；不加任何前缀
  extracted <- str_extract(df$fn, "^[^_]+")
  extracted <- str_remove(extracted, "\\.png$")
  df$fn <- extracted
  
  # gene 列从字符串数组 -> 数值多列
  gene_list_chr <- str_remove_all(df$gene, "\\[|\\]")
  gene_list     <- str_split(gene_list_chr, ",\\s*")
  gene_num      <- lapply(gene_list, as.numeric)
  
  lens <- lengths(gene_num)
  if (length(unique(lens)) != 1) {
    message("gene 数组长度不一致：", spot_csv_path)
    return(NULL)
  }
  p  <- unique(lens)
  m  <- do.call(rbind, gene_num)
  md <- as.data.frame(m, check.names = FALSE)
  colnames(md) <- paste0("g", seq_len(p))
  
  bind_cols(df[1], df["fn"], md)
}

# 把数据框的“数值区”（第3..末列）取为矩阵
num_mat <- function(df) as.matrix(df[, -(1:2), drop = FALSE])

# 浮点比较（带容差）
vec_equal_tol <- function(x, y, tol = 1e-8) {
  isTRUE(all.equal(as.numeric(x), as.numeric(y), tolerance = tol))
}

# 比对规则：先比第1行 3..末列，再比所有行 3..末列，均为 TRUE 才通过
compare_by_rule <- function(spot_df, gt_df, spot_name, gt_name, tol = 1e-8) {
  sm <- num_mat(spot_df)
  gm <- num_mat(gt_df)
  
  # 维度：spot 行数 = gm 列数；spot 列数(基因数) = gm 行数
  if (nrow(sm) != ncol(gm) || ncol(sm) != nrow(gm)) {
    message("维度不匹配：", spot_name, " vs ", gt_name,
            "；spot 行/列=", paste(dim(sm), collapse = "x"),
            "；gt 数值区 行/列=", paste(dim(gm), collapse = "x"))
    return(FALSE)
  }
  
  # 先比首行
  if (!vec_equal_tol(sm[1, ], gm[, 1], tol = tol)) {
    message("首行不一致：", spot_name, " [1,3:末]  vs  ", gt_name, " [列1]")
    return(FALSE)
  }
  
  # 再比所有行
  ok <- TRUE
  for (i in seq_len(nrow(sm))) {
    if (!vec_equal_tol(sm[i, ], gm[, i], tol = tol)) {
      message("存在不一致的行 i=", i, "：", spot_name, " vs ", gt_name)
      ok <- FALSE
      break
    }
  }
  ok
}

# 主流程：处理一个数据集目录及其对应 spot
process_one <- function(ds_name, spot_file) {
  ds_dir <- file.path(base_dir, ds_name)
  if (!dir.exists(ds_dir)) {
    message("数据集目录不存在：", ds_dir)
    return(invisible(NULL))
  }
  
  spot_path <- file.path(spot_dir, spot_file)
  spot_df   <- parse_spot_csv(spot_path)
  if (is.null(spot_df)) return(invisible(NULL))
  
  # 查找该目录下的 gt/pre 源文件
  paths <- find_gt_pre_files(ds_dir)
  gt_src <- paths$gt
  pre_src <- paths$pre
  if (is.null(gt_src)) { message("跳过（缺少 gt 源）：", ds_dir); return(invisible(NULL)) }
  if (is.null(pre_src)) { message("警告：缺少 pre 源：", ds_dir, "（将仅处理 gt）") }
  
  gt_df <- read_csv_safely(gt_src)
  if (is.null(gt_df) || ncol(gt_df) < 3) {
    message("gt 源无法使用：", gt_src)
    return(invisible(NULL))
  }
  pre_df <- if (!is.null(pre_src)) read_csv_safely(pre_src) else NULL
  
  # 比对
  ok <- compare_by_rule(spot_df, gt_df, basename(spot_path), basename(gt_src), tol = 1e-8)
  if (!isTRUE(ok)) {
    message("比对失败，未改名：", basename(spot_path), "  vs  ", basename(gt_src), "  @ ", ds_name)
    return(invisible(NULL))
  }
  
  # 用 spot 第2列（fn，已抽取）作为列名
  labels <- spot_df[[2]]
  # 列数校验
  if ((ncol(gt_df) - 2) != length(labels)) {
    message("GT 样本列数不匹配：", gt_src,
            "；GT样本列=", ncol(gt_df) - 2, "；labels=", length(labels))
    return(invisible(NULL))
  }
  if (!is.null(pre_df) && (ncol(pre_df) - 2) != length(labels)) {
    message("PRE 样本列数不匹配：", pre_src,
            "；PRE样本列=", ncol(pre_df) - 2, "；labels=", length(labels))
    return(invisible(NULL))
  }
  
  new_gt <- gt_df
  colnames(new_gt)[3:ncol(new_gt)] <- labels
  gt_out <- file.path(ds_dir, "gt_with_gene.csv")
  readr::write_csv(new_gt, gt_out)
  message("已写出：", gt_out)
  
  if (!is.null(pre_df)) {
    new_pre <- pre_df
    colnames(new_pre)[3:ncol(new_pre)] <- labels
    pre_out <- file.path(ds_dir, "pre_with_gene.csv")
    readr::write_csv(new_pre, pre_out)
    message("已写出：", pre_out)
  }
}

# ========== 执行 ==========
walk(names(ds_map), ~ process_one(.x, ds_map[[.x]]))

message("全部处理完成。")
