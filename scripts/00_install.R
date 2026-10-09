# 00_install.R
# GO解析に必要なパッケージをインストールする（最初に1回だけ実行）

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(
  c("clusterProfiler", "org.Mm.eg.db", "enrichplot"),
  update = FALSE, ask = FALSE
)

# 読み込めるか確認（エラーが出なければOK）
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
cat("\nインストール完了\n")
