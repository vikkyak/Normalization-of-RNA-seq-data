require(recount)
set.seed(1)
param_ID <- "SRP001540"          #ID for Pickrell's data
n1 <- 40                         #number of replicates for group 1
n2 <- 29                         #number of replicates for group 2
FDR <- 0.1                       #FDR threshold
# download_study(param_ID, type="rse-gene", download=T)
load(file.path("/home/user/Desktop/Papers/Single_cell/New_TMM/SRP001540", 'rse_gene.Rdata'))
rse <- scale_counts(rse_gene)
x <- assays(rse)$counts
x <- as.data.frame(x)
###  Collapsing the data for technical replicates  ###
data <- cbind(x$SRR031822, x$SRR031953 + x$SRR031873,#Female1-2
              x$SRR031952 + x$SRR031871, x$SRR031868,#Female3-4
              x$SRR031819, x$SRR031897 + x$SRR031857,#Female5-6
              x$SRR031823, x$SRR031959, x$SRR031955,#Female7-9
              x$SRR031954, x$SRR031956, x$SRR031838,#Female10-12
              x$SRR031918, x$SRR031817, x$SRR031949 + x$SRR031852,#Female13-15
              x$SRR031841, x$SRR031865, x$SRR031896,#Female16-18
              x$SRR031853, x$SRR031820, x$SRR031874,#Female19-21
              x$SRR031895, x$SRR031870, x$SRR031839,#Female22-24
              x$SRR031958, x$SRR031867, x$SRR031848,#Female25-27
              x$SRR031847, x$SRR031818, x$SRR031919,#Female28-30
              x$SRR031866, x$SRR031849, x$SRR031877,#Female31-33
              x$SRR031814, x$SRR031914, x$SRR031812,#Female34-36
              x$SRR031842, x$SRR031843, x$SRR031860, x$SRR031837,#Female37-40
              x$SRR031917, x$SRR031821 + x$SRR031898,#Male1-2
              x$SRR031950 + x$SRR031850, x$SRR031876 + x$SRR031862,#Male3-4
              x$SRR031875, x$SRR031915, x$SRR031878 + x$SRR031863,#Male5-7
              x$SRR031869, x$SRR031864, x$SRR031845,#Male8-10
              x$SRR031951 + x$SRR031851, x$SRR031846,#Male11-12
              x$SRR031916, x$SRR031844, x$SRR031813,#Male13-15
              x$SRR031894, x$SRR031854, x$SRR031858,#Male16-18
              x$SRR031859, x$SRR031872, x$SRR031816,#Male19-21
              x$SRR031815, x$SRR031920 + x$SRR031899,#Male22-23
              x$SRR031957 + x$SRR031855, x$SRR031840,#Male24-25
              x$SRR031948, x$SRR031893, x$SRR031811, x$SRR031861)#Male26-29


colnames(data) <- c(paste("Female", 1:40, sep=""), paste("Male", 1:29, sep=""))
rownames(data) <- rownames(x)

###  Preparation for class labels  ###
data.cl <- c(rep(1, n1), rep(2, n2))
group <- data.cl
group <- factor(group)
counts <- data
require(ROC)
require(pROC)
require(edgeR)
y <- edgeR::DGEList(counts = counts, group = group)
###  Filtering low count genes  ##
keep = filterByExpr(y, group)
y <- y[keep, keep.lib.sizes = FALSE]
###  edgeR  ###
yG <- edgeR::calcNormFactors(y)
design <- model.matrix( ~ group)
yG <- edgeR::estimateDisp(yG, design)
yG <- edgeR::glmQLFit(yG, design)
yG <- edgeR::glmQLFTest(yG, coef = ncol(yG$design))
res <- edgeR::topTags(yG, n = nrow(counts), sort.by = "none")
res$padj <- p.adjust(res$table$PValue, method = "BH")
outcome = rep(1, dim(res$table)[1])
outcome[abs(res$table[, "logFC"]) <= 0.5] = 0
rocEdgeG <- roc(outcome, res$padj)
###  DESeq2  ###
require(DESeq2)
counts <- y$counts
colData <- data.frame(condition = as.factor(group))
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = colData,
                              design =  ~ condition)
dds <- DESeq(dds)
res <- results(dds,  independentFiltering = F)
outcome = rep(1, dim(res)[1])
outcome[abs(res$log2FoldChange) <= 0.5] = 0
rocDESeq <- roc(outcome, res$padj)

# Total Count
require(psych)
counts <- y$counts
totCounts <- colSums(counts)
sizeFactors(dds) <- totCounts / geometric.mean(totCounts)
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds)
res <- results(dds)
outcome = rep(1, dim(res)[1])
outcome[abs(res$log2FoldChange) <= 0.5] = 0
rocTOt <- roc(outcome, res$padj)

## TCC estimate
require(TCC)
counts <- y$counts
cond <- group
tcc <- new("TCC", counts, cond)
tcc <-
  TCC::calcNormFactors(
    tcc,
    norm.method = "tmm",
    test.method = "edger",
    iteration = 3,
    FDR = 0.1,
    floorPDEG = 0.05
  )
tccEsts <- tcc$norm.factors
tccEsts <-
  tccEsts * colSums(counts) / geometric.mean(colSums(counts))
sizeFactors(dds) <- tccEsts / geometric.mean(tccEsts)
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds)
res <- results(dds)
outcome = rep(1, dim(res)[1])
outcome[abs(res$log2FoldChange) <= 0.5] = 0
rocTCC <- roc(outcome, res$padj)

## PoissonSeq estimate
require(PoissonSeq)
counts <- y$counts
psests <- PS.Est.Depth(counts)
sizeFactors(dds) <- psests / geometric.mean(psests)
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds)
res <- results(dds)
outcome = rep(1, dim(res)[1])
outcome[abs(res$log2FoldChange) <= 0.5] = 0
rocPS <- roc(outcome, res$padj)
# proposed
source("~/New_TMM/Source.R", echo = TRUE)
require(compositions)
counts <- y$counts
pair <- c("1", "2")
dge <- DGEList(counts = counts, group = group)
Pro = DgEList(dge,
              optimal = TRUE,
              method = "TMM",
              refColumn = "1")
ProE <- edgeR::estimateCommonDisp(Pro)
ProE <- edgeR::estimateTagwiseDisp(ProE, trend = "movingave")
ProE <- exactTest(ProE, dispersion = "tagwise", pair = pair)
res <- topTags(ProE, n = nrow(ProE), sort.by = "none")
res$padj <- p.adjust(res$table$PValue, method = "BH")
outcome = rep(1, dim(res$table)[1])
outcome[abs(res$table[, "logFC"]) <= 0.5] = 0
rocProE <- roc(outcome, res$padj)
## ROC plot
plot(rocDESeq, main = "ROC Curves for Pickrell Data")
# lines(rocTOt, col = "blue", lty = 2)
lines(rocEdgeG, col = "green", lty = 3)
lines(rocProE, col = "red", lty = 4)
lines(rocPS, col = "pink", lty = 5)
lines(rocTCC, col = "yellow", lty = 6)
legend(
  "bottomright",
  c("Proposed",
    "TMM",
    "DESeq",
    # "Total Count",
    "PoissonSeq",
    "DEGES"),
  lty = c(4, 3, 1,  5, 6),
  col = c("red", "green", "black",   "pink", "yellow"),
  cex = c(0.7, 0.7, 0.7, 0.7, 0.7),
  lwd = 2
)

ResulT <-
  list(
    ProG = rocProE,
    EdgeG = rocEdgeG,
    PS = rocPS,
    TCC = rocTCC,
    # TOT = rocTOt,
    DESeq = rocDESeq
  )

# write.table(ResulT, file='~/New_TMM//Results/pickrell.txt')


