# 04_go_all_ontologies.R
# BP（働き）・CC（場所）・MF（分子の機能）の3分類を、上昇・低下別にGO解析する
# 背景遺伝子 = DEG表に載っている全遺伝子

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(enrichplot)
})
options(width = 200)

deg <- read.csv("data/DEG.csv", check.names = FALSE)
colnames(deg)[1] <- "gene"

# ---- 設定（ここを変えれば条件を変えられる）----
fdr_cut   <- 0.05
logfc_cut <- 0      # 1.5倍以上に絞るなら 0.58
n_show    <- 15     # 表示する上位ターム数

# ---- 遺伝子リスト ----
universe  <- unique(deg$gene)
gene_sets <- list(
  up   = unique(deg$gene[deg$FDR < fdr_cut & deg$logFC >  logfc_cut]),
  down = unique(deg$gene[deg$FDR < fdr_cut & deg$logFC < -logfc_cut])
)
cat("背景:", length(universe), " 上昇:", length(gene_sets$up),
    " 低下:", length(gene_sets$down), "\n\n")

# ---- 関数 ----
has_terms <- function(ego) !is.null(ego) && nrow(as.data.frame(ego)) > 0
ratio <- function(x) sapply(strsplit(x, "/"), function(v) as.numeric(v[1]) / as.numeric(v[2]))

run_go <- function(genes, ont) {
  ego <- enrichGO(gene = genes, universe = universe, OrgDb = org.Mm.eg.db,
                  keyType = "SYMBOL", ont = ont, pAdjustMethod = "BH",
                  pvalueCutoff = 0.05, qvalueCutoff = 1,
                  minGSSize = 10, maxGSSize = 500)
  if (!has_terms(ego)) return(ego)
  simplify(ego, cutoff = 0.7, by = "p.adjust", select_fun = min)
}

show_top <- function(ego, label) {
  cat("====", label, "====\n")
  if (!has_terms(ego)) { cat("有意なタームなし\n\n"); return(invisible(NULL)) }
  res <- as.data.frame(ego)
  res$Fold <- round(ratio(res$GeneRatio) / ratio(res$BgRatio), 2)
  out <- head(res[order(res$p.adjust), c("Description", "Count", "Fold", "p.adjust")], n_show)
  out$p.adjust <- signif(out$p.adjust, 2)
  rownames(out) <- NULL
  print(out, right = FALSE)
  cat("有意なターム数（整理後）:", nrow(res), "\n\n")
}

# ---- 解析（3分類 × 上昇・低下 = 6通り）----
labels <- c(BP = "働き", CC = "場所", MF = "分子の機能")
dir_jp <- c(up = "KOで上昇", down = "KOで低下")

dir.create("results", showWarnings = FALSE)
summary_rows <- list()
pdf("results/go_dotplot_all.pdf", width = 9, height = 7)

for (ont in c("BP", "CC", "MF")) {
  for (d in c("up", "down")) {
    ego <- run_go(gene_sets[[d]], ont)
    show_top(ego, paste0(ont, "（", labels[ont], "）／", dir_jp[d]))

    n_terms <- if (has_terms(ego)) nrow(as.data.frame(ego)) else 0
    summary_rows[[length(summary_rows) + 1]] <- data.frame(分類 = ont, 方向 = dir_jp[d], 有意なターム数 = n_terms)

    if (has_terms(ego)) {
      write.csv(as.data.frame(ego),
                paste0("results/go_", tolower(ont), "_", d, ".csv"), row.names = FALSE)
      print(dotplot(ego, showCategory = n_show) + ggplot2::ggtitle(paste(ont, d)))
    }
  }
}
dev.off()

cat("==== まとめ ====\n")
summary_table <- do.call(rbind, summary_rows)
rownames(summary_table) <- NULL
print(summary_table, right = FALSE)
cat("\nresults/ に表（CSV）と図（go_dotplot_all.pdf）を保存しました\n")
