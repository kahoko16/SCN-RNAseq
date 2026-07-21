# =========================================================
# SCN snRNA-seq analysis: Gpr176-positive cells
# -> Fraction of Avpr1a / Adcyap1r1(Pac1) / Mtnr1a / Vipr2 positive cells
#    within Gpr176+ cells (reverse of the previous analysis)
#
# 全細胞・サブタイプ別（VIP, AVP, CCK, Other）に、
# Gpr176陽性細胞のうち各受容体陽性細胞の割合を計算する。
#
# 受容体遺伝子:
#   Avpr1a    : AVP受容体 V1a
#   Adcyap1r1 : PACAP受容体 (Pac1)
#   Mtnr1a    : メラトニン受容体 1A ("Mtr1a" として依頼されたものに対応)
#   Vipr2     : VIP受容体 VPAC2
# =========================================================

# -------------------------
# 0. Packages
# -------------------------
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyverse)

# -------------------------
# 1. Load Seurat object (script1で保存したもの)
# -------------------------
scn_neu <- readRDS("scn_neu.rds")
Idents(scn_neu) <- "cluster2"

scn_dark  <- subset(scn_neu, subset = Condition == "Dark")
scn_light <- subset(scn_neu, subset = Condition == "Light")

# -------------------------
# 2. サブタイプ分類（既存スクリプトと同じ定義）
# -------------------------
assign_group <- function(cluster2) {
  group <- rep("Other", length(cluster2))
  group[cluster2 %in% c("Neu01", "Neu05", "Neu07")] <- "VIP"
  group[cluster2 %in% c("Neu02", "Neu11", "Neu18")] <- "AVP"
  group[cluster2 %in% c("Neu09")] <- "CCK"
  group
}

# -------------------------
# 3. Gpr176陽性細胞中の受容体陽性割合を計算する関数
#    - group = "All"   : SCN全体（Gpr176陽性細胞全て）
#    - group = 各サブタイプ : そのサブタイプのGpr176陽性細胞のみ
# -------------------------
calc_receptor_in_gpr176 <- function(seu, receptor_gene) {
  expr <- GetAssayData(seu, layer = "data")

  if (!(receptor_gene %in% rownames(expr))) {
    warning(paste0(receptor_gene, " not found in this object; skipped."))
    return(NULL)
  }

  df <- data.frame(
    cluster2      = as.character(seu$cluster2),
    gpr176_pos    = as.vector(expr["Gpr176", ] > 0),
    receptor_pos  = as.vector(expr[receptor_gene, ] > 0)
  )
  df$group <- assign_group(df$cluster2)

  gpr176_cells <- subset(df, gpr176_pos)

  # サブタイプ別
  by_group <- aggregate(receptor_pos ~ group, data = gpr176_cells, mean)

  # 全体（All）
  overall <- data.frame(
    group        = "All",
    receptor_pos = mean(gpr176_cells$receptor_pos)
  )

  result <- rbind(overall, by_group)
  result$gene <- receptor_gene
  result
}

# -------------------------
# 4. 対象遺伝子
# -------------------------
receptor_genes <- c(
  Avpr1a    = "Avpr1a",
  Pac1      = "Adcyap1r1",
  Mtnr1a    = "Mtnr1a",
  Vipr2     = "Vipr2"
)

# -------------------------
# 5. Dark / Light それぞれで実行し、結果をまとめる
# -------------------------
df_list <- list()

for (gene_label in names(receptor_genes)) {
  gene <- receptor_genes[[gene_label]]

  res_dark  <- calc_receptor_in_gpr176(scn_dark,  gene)
  res_light <- calc_receptor_in_gpr176(scn_light, gene)

  if (!is.null(res_dark))  res_dark$condition  <- "Dark"
  if (!is.null(res_light)) res_light$condition <- "Light"

  df_list[[gene_label]] <- bind_rows(res_dark, res_light) %>%
    mutate(gene_label = gene_label)
}

df_receptor_in_gpr176 <- bind_rows(df_list)

# group の順序を固定 (All -> VIP -> AVP -> CCK -> Other)
df_receptor_in_gpr176$group <- factor(
  df_receptor_in_gpr176$group,
  levels = c("All", "VIP", "AVP", "CCK", "Other")
)

print(df_receptor_in_gpr176)

# -------------------------
# 6. 共通プロット設定
# -------------------------
theme_all <- theme_classic() +
  theme(
    plot.title   = element_text(size = 22, face = "bold", hjust = 0.5),
    axis.title   = element_text(size = 18),
    axis.text    = element_text(size = 14),
    axis.text.x  = element_text(size = 14, face = "bold", angle = 45, hjust = 1),
    axis.text.y  = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text  = element_text(size = 14),
    strip.text   = element_text(size = 16, face = "bold")
  )

fill_condition <- scale_fill_manual(values = c("Dark" = "black", "Light" = "orange"))

# -------------------------
# 7. Barplot: 遺伝子ごとにfacet、x = サブタイプ、fill = 条件
# -------------------------
p_receptor_in_gpr176 <- ggplot(
  df_receptor_in_gpr176,
  aes(x = group, y = receptor_pos, fill = condition)
) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_wrap(~ gene_label, nrow = 1) +
  labs(
    title = "Fraction of receptor+ cells among Gpr176+ cells",
    y = "Fraction of receptor+ (in Gpr176+ cells)",
    x = "Cell type"
  ) +
  fill_condition +
  theme_all

p_receptor_in_gpr176

# 遺伝子ごとに個別プロットが欲しい場合
plots_by_gene <- lapply(names(receptor_genes), function(gene_label) {
  df_sub <- df_receptor_in_gpr176 %>% filter(gene_label == !!gene_label)

  ggplot(df_sub, aes(x = group, y = receptor_pos, fill = condition)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    labs(
      title = paste0(gene_label, " in Gpr176+ cells"),
      y = "Fraction of receptor+",
      x = "Cell type"
    ) +
    fill_condition +
    theme_all
})
names(plots_by_gene) <- names(receptor_genes)

# 例: Avpr1aだけ表示
# plots_by_gene[["Avpr1a"]]

# =========================================================
# 8. 受容体発現確認ドットプロット（Mtnr1a, Vipr2を追加）
#    ggplot_analysis.R の features_full を拡張したもの
# =========================================================

features_receptor <- c(
  "Vip", "Avp", "Nms", "Cck",
  "Gpr176", "Avpr1a", "Adcyap1r1", "Mtnr1a", "Vipr2"
)
features_receptor <- features_receptor[
  features_receptor %in% rownames(GetAssayData(scn_neu, layer = "data"))
]

# クラスター並び順をDarkの著者マーカー発現で決める（ggplot_analysis.R と同じ方法）
features_paper <- c(
  "Gad1", "Gad2", "Slc32a1",
  "Vip", "Avp", "Nms", "Cck",
  "Syt7", "Ttr", "Rorb", "Rgs16", "Penk", "Vgf"
)

avg_dark_for_order <- AverageExpression(
  scn_dark,
  features = features_paper,
  group.by = "cluster2"
)$RNA

avg_dark_for_order_scaled <- t(apply(avg_dark_for_order, 1, function(x) {
  rng <- range(x, na.rm = TRUE)
  if (diff(rng) == 0) rep(0, length(x)) else (x - rng[1]) / diff(rng)
}))

hc <- hclust(dist(t(avg_dark_for_order_scaled)))
cluster_order_auto <- hc$labels[hc$order]

make_dotplot_df <- function(seu, features, cluster_order) {
  avg_exp <- AverageExpression(seu, features = features, group.by = "cluster2")$RNA

  genes_present <- intersect(features, rownames(avg_exp))
  clusters_present <- colnames(avg_exp)
  avg_exp <- avg_exp[genes_present, clusters_present, drop = FALSE]

  avg_scaled <- t(apply(avg_exp, 1, function(x) {
    rng <- range(x, na.rm = TRUE)
    if (diff(rng) == 0) rep(0, length(x)) else (x - rng[1]) / diff(rng)
  }))
  colnames(avg_scaled) <- clusters_present
  rownames(avg_scaled) <- genes_present

  expr <- GetAssayData(seu, layer = "data")
  expr <- expr[genes_present, , drop = FALSE]

  pct_exp <- sapply(genes_present, function(g) {
    vals <- tapply(expr[g, ] > 0, as.character(seu$cluster2), mean)
    vals[clusters_present]
  })
  pct_exp <- t(pct_exp)
  colnames(pct_exp) <- clusters_present
  rownames(pct_exp) <- genes_present

  df_plot <- expand.grid(
    gene = genes_present,
    cluster = clusters_present,
    stringsAsFactors = FALSE
  )
  df_plot$avg <- as.vector(avg_scaled)
  df_plot$pct <- as.vector(pct_exp)

  cluster_order_use <- intersect(cluster_order, unique(df_plot$cluster))
  df_plot$cluster <- factor(df_plot$cluster, levels = cluster_order_use)
  df_plot$gene <- factor(df_plot$gene, levels = genes_present)
  df_plot <- df_plot[!is.na(df_plot$cluster), , drop = FALSE]

  df_plot
}

df_dot_dark_receptor <- make_dotplot_df(
  seu = scn_dark,
  features = features_receptor,
  cluster_order = cluster_order_auto
)

df_dot_light_receptor <- make_dotplot_df(
  seu = scn_light,
  features = features_receptor,
  cluster_order = cluster_order_auto
)

theme_dot <- theme_classic() +
  theme(
    plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

p_dark_receptor <- ggplot(df_dot_dark_receptor, aes(x = gene, y = cluster)) +
  geom_point(aes(size = pct, color = avg)) +
  scale_color_gradient(low = "white", high = "blue", limits = c(0, 1)) +
  scale_size(range = c(0, 5), limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75)) +
  labs(
    title = "Dark",
    x = "Features",
    y = "Cluster",
    color = "Average Expression",
    size = "Percent Expressed"
  ) +
  theme_dot

p_light_receptor <- ggplot(df_dot_light_receptor, aes(x = gene, y = cluster)) +
  geom_point(aes(size = pct, color = avg)) +
  scale_color_gradient(low = "white", high = "blue", limits = c(0, 1)) +
  scale_size(range = c(0, 5), limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75)) +
  labs(
    title = "Light",
    x = "Features",
    y = "Cluster",
    color = "Average Expression",
    size = "Percent Expressed"
  ) +
  theme_dot

p_dark_receptor
p_light_receptor
