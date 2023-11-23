####  Parameters  ###

for (PDEG in seq(0.05, 0.75, 0.1)) {
  for (P1 in seq(0.2, 0.8, 0.1)) {
    N_trial <- 1            #number of trials
    G <- 10000                 #number of genes
    n1 <- 3                    #number of replicates for group 1
    n2 <- 3                    #number of replicates for group 2
    # PDEG <- 0.05               #proportion of DEG
    # P1 <- 0.40                  #proportion of up-regulated DEGs in group 1
    FC <- 4                    #degree of fold change (FC)
    
    ###  Load packages  ###
    library(edgeR)
    library(DESeq2)
    library(TCC)
    library(MBCluster.Seq)
    library(ROC)
    library(compositions)
    source("~/Desktop/Papers/Single_cell/New_TMM/Source.R", echo = TRUE)
    ###################
    ###  Main loop  ###
    ###################
    matome <- NULL
    for (i in 1:N_trial) {
      print(i)
      
      ###  Generation of simulated data  ###
      set.seed(i)
      fc.matrix <- makeFCMatrix(
        Ngene = G,
        PDEG = PDEG,
        DEG.assign = c(P1, 1 - P1),
        replicates = c(n1, n2)
      )
      tcc <- simulateReadCounts(
        Ngene = G,
        PDEG = PDEG,
        DEG.assign = c(P1, 1 - P1),
        replicates = c(n1, n2),
        fc.matrix = fc.matrix,
        DEG.foldchange = c(FC, FC)
      )
      data <- tcc$count
      data.cl <- c(rep(1, n1), rep(2, n2))
      counts <- data
      group <- data.cl
      
      ###  edgeR  ###
      y <- edgeR::DGEList(counts = counts, group = group)
      y <- edgeR::calcNormFactors(y)
      design <- model.matrix( ~ group)
      y <- edgeR::estimateDisp(y, design)
      fit <- edgeR::glmQLFit(y, design)
      qlf <- edgeR::glmQLFTest(fit, coef = 2)
      res <- edgeR::topTags(qlf, n = nrow(counts), sort.by = "none")
      
      #Calculation of AUC value
      p.value <- res$table$PValue
      ranking <- rank(p.value)
      obj <- as.numeric(tcc$simulation$trueDEG != 0)
      auc <- AUC(rocdemo.sca(truth = obj, data = -ranking))
      auc_edger <- auc
      
      ###  DESeq2  ###
      colData <- data.frame(condition = as.factor(data.cl))
      d <- DESeqDataSetFromMatrix(
        countData = data,
        colData = colData,
        design =  ~ condition
      )
      d <- DESeq(d)
      tmp <- results(d)
      
      #Calculation of AUC value
      p.value <- tmp$pvalue
      p.value[is.na(p.value)] <- 1
      ranking <- rank(p.value)
      obj <- as.numeric(tcc$simulation$trueDEG != 0)
      auc <- AUC(rocdemo.sca(truth = obj, data = -ranking))
      auc_deseq2 <- auc
      
      ## Median  normalization
      norm.med<-med(counts) 
      dge <- DGEList(counts = norm.med, group = group)
      Med <- edgeR::calcNormFactors(dge, method="none") 
      Med <- edgeR::estimateDisp(Med, design)
      Med <- edgeR::estimateTagwiseDisp(Med, trend = "movingave")
      fit <- edgeR::glmQLFit(Med, design)
      Qlf <- edgeR::glmQLFTest(fit, coef = 2)
      Res <- edgeR::topTags(Qlf, n = nrow(counts), sort.by = "none")
      #Calculation of AUC value
      P.value <- Res$table$PValue
      ranking <- rank(P.value)
      Obj <- as.numeric(tcc$simulation$trueDEG != 0)
      auc <- AUC(rocdemo.sca(truth = Obj, data = -ranking))
      auc_Med <- auc
      
      ##  upper quartile normalization
      norm.uq<-uq(counts)
      dge <- DGEList(counts = norm.uq, group = group)
      Uq <- edgeR::calcNormFactors(dge, method="none") 
      Uq <- edgeR::estimateDisp(Uq, design)
      Uq <- edgeR::estimateTagwiseDisp(Uq, trend = "movingave")
      fit <- edgeR::glmQLFit(Uq, design)
      Qlf <- edgeR::glmQLFTest(fit, coef = 2)
      Res <- edgeR::topTags(Qlf, n = nrow(counts), sort.by = "none")
      #Calculation of AUC value
      P.value <- Res$table$PValue
      ranking <- rank(P.value)
      Obj <- as.numeric(tcc$simulation$trueDEG != 0)
      auc <- AUC(rocdemo.sca(truth = Obj, data = -ranking))
      auc_Uq <- auc
      
      ## full quantile normalization (FQ)  
      norm.fq<-fq(counts)
      dge <- DGEList(counts = norm.fq, group = group)
      Fq <- edgeR::calcNormFactors(dge, method="none") 
      Fq <- edgeR::estimateDisp(Fq, design)
      Fq <- edgeR::estimateTagwiseDisp(Fq, trend = "movingave")
      fit <- edgeR::glmQLFit(Fq, design)
      Qlf <- edgeR::glmQLFTest(fit, coef = 2)
      Res <- edgeR::topTags(Qlf, n = nrow(counts), sort.by = "none")
      #Calculation of AUC value
      P.value <- Res$table$PValue
      ranking <- rank(P.value)
      Obj <- as.numeric(tcc$simulation$trueDEG != 0)
      auc <- AUC(rocdemo.sca(truth = Obj, data = -ranking))
      auc_Fq <- auc
      
      # proposed
      dge <- DGEList(counts = counts, group = group)
      Pro = DgEList(dge,
                    optimal = TRUE,
                    method = "TMM",
                    refColumn = "1")
      design <- model.matrix( ~ group)
      Pro <- edgeR::estimateDisp(Pro, design)
      Pro <- edgeR::estimateTagwiseDisp(Pro, trend = "movingave")
      fit <- edgeR::glmQLFit(Pro, design)
      Qlf <- edgeR::glmQLFTest(fit, coef = 2)
      Res <- edgeR::topTags(Qlf, n = nrow(counts), sort.by = "none")
      #Calculation of AUC value
      P.value <- Res$table$PValue
      ranking <- rank(P.value)
      Obj <- as.numeric(tcc$simulation$trueDEG != 0)
      auc <- AUC(rocdemo.sca(truth = Obj, data = -ranking))
      auc_Pro <- auc
  
      
      matome <-
        rbind(matome, c(auc_edger, auc_deseq2, auc_Med, auc_Uq, auc_Fq, auc_Pro))
    }
    
    #Output(AUC)
    colnames(matome) <- c("edgeR", "DESeq2", "Med", "Uq", "Fq", "Pro")
    # matome
    out_f <- paste("Add3_", PDEG, "_", P1, "_", n1, "_gamma.txt", sep = "")
    write.table(
      matome,
      out_f,
      sep = "\t",
      append = F,
      quote = F,
      row.names = F
    )
    # summary(matome)
    
  }
  
}

# par(mfrow=c(1,2))
# boxplot(summary(matome[,1]), summary(matome[,2]), summary(matome[,3]),  summary(matome[,4]),
#         main = "ROC Plot of Two-group Simulated Data",
#         names = c("TMM", "DESeq", "TCC", "Proposed"),
#         col=c("red", "orange", "green", "yellow")
# )
# boxplot(summary(matome[,1]), summary(matome[,3]), summary(matome[,4]),  summary(matome[,2]),
#         main = "ROC Plot of Two-group Simulated Data",
#         names = c("TMM", "DESeq", "TCC", "Proposed"),
#         col=c("red", "black", "green", "blue")
# )
# dev.off()
