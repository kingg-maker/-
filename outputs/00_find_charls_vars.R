# ============================================================
# 00_find_charls_vars.R
# 扫描已下载的原始 CHARLS 文件，定位研究所需变量所在文件
# 用法：
#   1) 修改 raw_dir 为你存放下载文件的文件夹
#   2) 运行本脚本
#   3) 把生成的 charls_var_report.csv 发回给我
# ============================================================

library(haven)
library(readxl)
library(tidyverse)

# ---------- 配置 ----------
raw_dir <- "C:/Users/26301/Desktop/CHARLS_raw"   # 改为你的下载目录
out_csv <- "C:/Users/26301/Documents/Codex/2026-07-29/2026-2026-07-24-2026-10/outputs/charls_var_report.csv"

# 目标变量（前缀匹配，容忍不同波次的变量名变体）
targets <- c(
  "ba002", "bd001", "cb060",
  "ha065", "i024", "da056", "ed005", "ee003",
  "internet", "computer", "smart"
)

# ---------- 读取变量名（不加载数据，快） ----------
read_varnames <- function(path) {
  ext <- tolower(tools::file_ext(path))
  tryCatch({
    if (ext %in% c("dta")) return(names(haven::read_dta(path, n_max = 0)))
    if (ext %in% c("sav", "zsav")) return(names(haven::read_sav(path, n_max = 0)))
    if (ext == "csv") {
      hdr <- readLines(path, n = 1, encoding = "UTF-8")
      return(trimws(strsplit(hdr, ",")[[1]]))
    }
    if (ext %in% c("xls", "xlsx")) {
      d <- readxl::read_excel(path, n_max = 1)
      return(names(d))
    }
    NULL
  }, error = function(e) NULL)
}

# ---------- 扫描 ----------
files <- list.files(raw_dir, recursive = TRUE, full.names = TRUE)
files <- files[tolower(tools::file_ext(files)) %in%
                 c("dta", "sav", "zsav", "csv", "xls", "xlsx")]
cat("找到文件数:", length(files), "\n")

report <- list()
for (f in files) {
  vars <- read_varnames(f)
  if (is.null(vars) || length(vars) == 0) next
  hit <- targets[sapply(targets, function(t)
    any(grepl(t, vars, ignore.case = TRUE)))]
  if (length(hit) == 0) next
  report[[f]] <- data.frame(
    file = f,
    found_vars = paste(hit, collapse = "; "),
    has_ID = any(vars == "ID"),
    has_householdID = any(vars == "householdID"),
    n_vars = length(vars),
    stringsAsFactors = FALSE
  )
}

report_df <- bind_rows(report)
if (nrow(report_df)) {
  write_csv(report_df, out_csv)
  print(report_df %>% select(file, found_vars, has_ID, has_householdID))
  cat("\n报告已保存:", out_csv, "\n")
} else {
  cat("未找到任何目标变量，请检查 raw_dir 路径和文件格式。\n")
}
