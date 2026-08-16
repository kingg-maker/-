# ============================================================
# Build a harmonized 0-31 cognitive health score from raw CHARLS.
# Partial-missing version: an orientation item that is missing/refused is
# counted as 0; the orientation subscale is NA only if all five items are
# missing. This keeps more respondents than the strict version.
# Output: cognition_31_partial.csv.
# ============================================================

suppressPackageStartupMessages({
  library(haven)
  library(tidyverse)
})

raw_root <- "C:/Users/26301/Desktop/新增数据"
out_root <- "C:/Users/26301/Documents/Codex/2026-08-07/w/outputs/rq1_revised_v2"

missing_codes <- c(97, 98, 99, 997, 998, 999)

clean_id <- function(x, len = 12) {
  x <- as.character(x)
  ifelse(grepl("^[0-9]+$", x) & nchar(x) < len,
         str_pad(x, len, "left", "0"), x)
}

convert_id_2011 <- function(x) {
  x <- clean_id(x, 11)
  ifelse(nchar(x) == 11,
         paste0(substr(x, 1, 10), "0", substr(x, 11, 11)),
         x)
}

score_orient <- function(d, cols, correct) {
  mat <- as.matrix(d[cols])
  mode(mat) <- "numeric"
  mat[mat %in% missing_codes] <- NA
  score <- rowSums(
    sapply(seq_along(correct), function(j) {
      as.numeric(mat[, j] == correct[j])
    }),
    na.rm = TRUE
  )
  score[rowSums(is.na(mat)) == length(cols)] <- NA
  score
}

score_recall <- function(d, cols, has_zero = FALSE) {
  mat <- as.matrix(d[cols])
  mode(mat) <- "numeric"
  mat[mat %in% missing_codes] <- NA
  all_missing <- rowSums(is.na(mat)) == length(cols)
  if (has_zero) {
    score <- rowSums(mat > 0, na.rm = TRUE)
  } else {
    score <- rowSums(!is.na(mat))
  }
  score[all_missing] <- NA
  score
}

score_draw <- function(d, col) {
  x <- d[[col]]
  x[x %in% missing_codes] <- NA
  ifelse(is.na(x), NA, as.numeric(x == 1))
}

score_serial <- function(d, cols, correct = c(93, 86, 79, 72, 65)) {
  mat <- as.matrix(d[cols])
  mode(mat) <- "numeric"
  mat[mat %in% missing_codes] <- NA
  any_missing <- rowSums(is.na(mat)) > 0
  score <- numeric(nrow(mat))
  for (j in seq_along(correct)) {
    score <- score + as.numeric(mat[, j] == correct[j])
  }
  score[any_missing] <- NA
  score
}

out <- list()

# 2011
d <- read_dta(file.path(raw_root, "2011", "health_status_and_functioning.dta"))
out[["2011"]] <- data.frame(
  ID = convert_id_2011(d$ID),
  wave = 1L,
  orient = score_orient(d, c("dc001s1", "dc001s2", "dc001s3", "dc002", "dc003"), c(1, 2, 3, 1, 1)),
  recall_imm = score_recall(d, paste0("dc006s", 1:10), has_zero = FALSE),
  recall_del = score_recall(d, paste0("dc027s", 1:10), has_zero = FALSE),
  serial = score_serial(d, paste0("dc0", 19:23)),
  draw = score_draw(d, "dc025")
)

# 2013
d <- read_dta(file.path(raw_root, "2013", "Health_Status_and_Functioning.dta"))
out[["2013"]] <- data.frame(
  ID = clean_id(d$ID),
  wave = 2L,
  orient = score_orient(d, c("dc001s1", "dc001s2", "dc001s3", "dc002", "dc003"), c(1, 2, 3, 1, 1)),
  recall_imm = score_recall(d, paste0("dc006_1_s", 1:10), has_zero = FALSE),
  recall_del = score_recall(d, paste0("dc027s", 1:10), has_zero = FALSE),
  serial = score_serial(d, paste0("dc0", 19:23)),
  draw = score_draw(d, "dc025")
)

# 2015
d <- read_dta(file.path(raw_root, "2015", "Health_Status_and_Functioning.dta"))
out[["2015"]] <- data.frame(
  ID = clean_id(d$ID),
  wave = 3L,
  orient = score_orient(d, c("dc001s1", "dc001s2", "dc001s3", "dc002", "dc003"), c(1, 2, 3, 1, 1)),
  recall_imm = score_recall(d, paste0("dc006s", 1:10), has_zero = FALSE),
  recall_del = score_recall(d, paste0("dc027s", 1:10), has_zero = FALSE),
  serial = score_serial(d, paste0("dc0", 19:23)),
  draw = score_draw(d, "dc025")
)

# 2018
d <- read_dta(file.path(raw_root, "2018", "Cognition.dta"))
out[["2018"]] <- data.frame(
  ID = clean_id(d$ID),
  wave = 4L,
  orient = score_orient(d, c("dc001_w4", "dc006_w4", "dc003_w4", "dc005_w4", "dc002_w4"), rep(1, 5)),
  recall_imm = score_recall(d, paste0("dc028_w4_s", 1:10), has_zero = TRUE),
  recall_del = score_recall(d, paste0("dc047_w4_s", 1:10), has_zero = TRUE),
  serial = score_serial(d, paste0("dc014_w4_", 1:5, "_1")),
  draw = score_draw(d, "dc024_w4")
)

# 2020
d <- read_dta(file.path(raw_root, "2020", "Health_Status_and_Functioning.dta"))
out[["2020"]] <- data.frame(
  ID = clean_id(d$ID),
  wave = 5L,
  orient = score_orient(d, c("dc001", "dc005", "dc003", "dc004", "dc002"), rep(1, 5)),
  recall_imm = score_recall(d, paste0("dc012_s", 1:10), has_zero = TRUE),
  recall_del = score_recall(d, paste0("dc028_s", 1:10), has_zero = TRUE),
  serial = score_serial(d, paste0("dc007_", 1:5, "_1")),
  draw = score_draw(d, "dc009")
)

cognition <- bind_rows(out) %>%
  mutate(
    cog_cap_31 = ifelse(
      rowSums(is.na(cbind(orient, recall_imm, recall_del, serial, draw))) > 0,
      NA_real_,
      orient + recall_imm + recall_del + serial + draw
    )
  ) %>%
  select(ID, wave, orient, recall_imm, recall_del, serial, draw, cog_cap_31)

write_csv(cognition, file.path(out_root, "cognition_31_partial.csv"))

cat("cognition rows:", nrow(cognition), "\n")
cat("cog_cap_31 nonmissing by wave:\n")
print(cognition %>% group_by(wave) %>% summarise(n = n(), nonmissing = sum(!is.na(cog_cap_31)), .groups = "drop"))
