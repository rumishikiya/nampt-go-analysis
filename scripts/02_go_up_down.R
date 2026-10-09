# 02_go_up_down.R
# 上昇・低下遺伝子を分けてGO解析（生物学的プロセス, BP）を行う
# 背景遺伝子 = DEG表に載っている全遺伝子（発現が検出された遺伝子）

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

# ---- 遺伝子リストを作る ----
universe <- unique(deg$gene)
up   <- unique(deg$gene[deg$FDR < fdr_cut & deg$logFC >  logfc_cut])
down <- unique(deg$gene[deg$FDR < fdr_cut & deg$logFC < -logfc_cut])
cat("背景:", length(universe), " 上昇:", length(up), " 低下:", length(down), "\n\n")

# ---- GO解析（FDR < 0.05、タームサイズ10〜500、重複タームを整理）----
run_go <- function(genes) {
  ego <- enrichGO(gene = genes, universe = universe, OrgDb = org.Mm.eg.db,
                  keyType = "SYMBOL", ont = "BP", pAdjustMethod = "BH",
                  pvalueCutoff = 0.05, qvalueCutoff = 1,
                  minGSSize = 10, maxGSSize = 500)
  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) return(ego)
  simplify(ego, cutoff = 0.7, by = "p.adjust", select_fun = min)
}

has_terms <- function(ego) !is.null(ego) && nrow(as.data.frame(ego)) > 0
ratio <- function(x) sapply(strsplit(x, "/"), function(v) as.numeric(v[1]) / as.numeric(v[2]))

# ---- 上位タームを表示（遺伝子名は出さない＝共有しやすい形）----
show_top <- function(ego, label, n = 20) {
  cat("====", label, "====\n")
  if (!has_terms(ego)) { cat("有意なタームなし\n\n"); return(invisible(NULL)) }
  res <- as.data.frame(ego)
  res$Fold <- round(ratio(res$GeneRatio) / ratio(res$BgRatio), 2)
  out <- head(res[order(res$p.adjust), c("Description", "Count", "Fold", "p.adjust")], n)
  out$p.adjust <- signif(out$p.adjust, 2)
  rownames(out) <- NULL
  print(out, right = FALSE)
  cat("有意なターム数（整理後）:", nrow(res), "\n\n")
}

ego_up   <- run_go(up)
ego_down <- run_go(down)

show_top(ego_up,   "KOで上昇した遺伝子")
show_top(ego_down, "KOで低下した遺伝子")

# ---- 結果を保存（results/ はgitの管理外）----
dir.create("results", showWarnings = FALSE)
if (has_terms(ego_up))   write.csv(as.data.frame(ego_up),   "results/go_bp_up.csv",   row.names = FALSE)
if (has_terms(ego_down)) write.csv(as.data.frame(ego_down), "results/go_bp_down.csv", row.names = FALSE)

pdf("results/go_dotplot.pdf", width = 9, height = 7)
if (has_terms(ego_up))   print(dotplot(ego_up,   showCategory = 15) + ggplot2::ggtitle("Up in KO"))
if (has_terms(ego_down)) print(dotplot(ego_down, showCategory = 15) + ggplot2::ggtitle("Down in KO"))
dev.off()

cat("results/ に表（CSV）と図（PDF）を保存しました\n")
