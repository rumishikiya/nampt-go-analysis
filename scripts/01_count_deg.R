# 01_count_deg.R
# DEG表の概要を確認する（データは手元のPCでのみ処理）

deg <- read.csv("data/DEG.csv", check.names = FALSE)

# 1列目を遺伝子名として扱う
colnames(deg)[1] <- "gene"

cat("全遺伝子数:", nrow(deg), "\n")
cat("列名:", paste(colnames(deg), collapse = ", "), "\n\n")

# 閾値の組み合わせごとにDEGの数を数える
summary_table <- expand.grid(FDR = c(0.05, 0.01), logFC = c(0, 0.58, 1))
summary_table$up   <- NA
summary_table$down <- NA

for (i in seq_len(nrow(summary_table))) {
  f <- summary_table$FDR[i]
  l <- summary_table$logFC[i]
  summary_table$up[i]   <- sum(deg$FDR < f & deg$logFC >  l, na.rm = TRUE)
  summary_table$down[i] <- sum(deg$FDR < f & deg$logFC < -l, na.rm = TRUE)
}
summary_table$total <- summary_table$up + summary_table$down

cat("閾値ごとのDEG数（logFC 0.58 ≒ 1.5倍、1 = 2倍）\n")
print(summary_table)

# ノックアウトの確認：Namptの行
cat("\nNamptの行\n")
print(deg[grepl("^Nampt$", deg$gene, ignore.case = TRUE), ])
