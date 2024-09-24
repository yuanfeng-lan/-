# 载入所需 R 包
library(fgsea) # 或者: library(gage)

# 读取基因表达矩阵数据
exprs_file <- "your_expression_matrix.csv"
exprs_mat <- read.csv(exprs_file, row.names = 1) # 数据中基因名应为行名

# 加载针对特定样品组的基因列表或一些预定义的基因集（如Kegg, Biocarta等）
geneSet_file <- "your_gene_set.txt"
geneSets <- gmtPathways(geneSet_file)

# 定义gsea分析的参数
nperm <- 10^4
significance_level <- 0.05 # 显著性水平
min_size <- 15 # 基因集最小大小
max_size <- 500 # 基因集最大大小

# 对上述所有基因集执行gsea分析并指定分析算法，如fgsea算法
results <- fgsea(exprs_mat, geneSets,
                 nperm = nperm,
                 minSize = min_size,
                 maxSize = max_size,  
                 BPPARAM = SnowParam(workers = 4))

# 处理GSEA分析结果
esList <- results$hits %>% group_by(pathway) %>% summarize(enrichment_score = medianNES )
geneSet_df <- data.frame(id = names(esList),
                         esList,                     # 中位数标准化富集得分
                         gene_size = geneSetSizes(geneSets[names(esList)])   # 基因集大小
)

# 按富集得分排序结果并对其检验结果是否显著
geneSet_df <- geneSet_df %>%
  arrange(-enrichment_score, -gene_size) %>%
  mutate(padj = p.adjust(pvalue, method = "fdr", n = length(enrichment_score)))

# 按8-bit位图格式（.gmt）输出结果
result_gmt_file <- "gsea_results.gmt"
writeGeneSets(geneSets[names(esList)], geneSet_df$enrichment_score, geneSet_df$id, file = result_gmt_file, description = "GSEA analysis results by R-fgsea")