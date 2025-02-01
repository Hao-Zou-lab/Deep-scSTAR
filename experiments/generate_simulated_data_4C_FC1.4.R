setwd("E:\\simulate_for_test")
library(R.matlab)
library(Seurat)

rm(list=ls())
graphics.off()

ncells <- 10000  # Number of cells

# Our simulation involves three cell types/components.
# Cells are distributed according to a bivariate normal in a 2-D biological subspace. 
# Each cell type has a different x/y center and a different SD.

xmus <- c(0,5,5)
xsds <- c(1,0.1,1)
ymus <- c(5,5,0)
ysds <- c(1,0.1,1)


# Note that the different centers should not lie on the same y=mx line; this represents populations that differ only in library size. 
# Such differences should not be present in normalized data, and will be eliminated by the cosine normalization step.
# The centers above are chosen so as to guarantee good separation between the different components.

#####################################
# Generating data for batch 1, with a given proportion of cells in each component. 

FC_DE <- 1.4
fileDir <- paste0("./4C_", as.character(FC_DE))
dir.create(file.path(fileDir), showWarnings = FALSE)
ITR <- 20
#####################################
for (itr in 1:ITR){
  RN <- ceiling(runif(4, 0, 50))
  prop1 <- c(1, 0, 0)
  comp1 <- sample(1:3, prob=prop1, size=ncells, replace=TRUE)
  
  # Sampling locations for cells in each component.
  set.seed(RN[1])
  samples1 <- cbind(rnorm(n=ncells, mean=xmus[comp1],sd=xsds[comp1]),
                    rnorm(n=ncells, mean=ymus[comp1],sd=ysds[comp1]))
  
  # Plot the true cluster locations for batch 1
  ref.cols <- c("blue", "brown1", "gold2")
  clust1 <- ref.cols[comp1]
  plot(samples1, pch=16, cex=1.5, col=clust1, main=paste0("Batch 1"))
  
  # Random projection to D dimensional space, to mimic high-dimensional expression data.
  ngenes <- 1000
  set.seed(RN[2])
  proj <- matrix(rnorm(ngenes*ncells), nrow=ngenes, ncol=2)
  A1 <- samples1 %*% t(proj)
  
  # Add normally distributed noise.
  A1 <- A1 + rnorm(ngenes*ncells)
  rownames(A1) <- paste0("Cell", seq_len(ncells), "-1")
  colnames(A1) <- paste0("Gene", seq_len(ngenes))
  
  
  
  outmat <- outcol <- list()
  
  prop2 <- prop1
  
  
  # Setting proportions of each of the three cell types in batch 2.
  comp2 <- sample(1:3, prob=prop2, size=ncells, replace=TRUE) 
  
  # Sampling locations for cells in each component.  
  set.seed(RN[3])
  samples2 <- cbind(rnorm(n=ncells, mean=xmus[comp2], sd=xsds[comp2]),
                    rnorm(n=ncells, mean=ymus[comp2], sd=ysds[comp2]))
  
  # Plot the true cluter labels for batch 2.
  clust2 <- ref.cols[comp2]
  plot(samples2, pch=16, cex=1.5, col=clust2, main=paste("Batch 2,", "comp."))
  
  # Random projection, followed by adding batch effects and random noise.
  #idx_DEgenes <- sample.int(ngenes, size = ngenes, replace = FALSE)
  #idx_DEgenes <- idx_DEgenes[1:floor(ngenes/4)]
  idx_DEgenes <- 1:175
  idx_responserCells <- sample.int(ncells, size = ncells, replace = FALSE)
  idx_responserCells <- idx_responserCells[1:floor(ncells/1.5)]
  A2 <- samples2 %*% t(proj) 
  A2 <- A2 + matrix(rep(rnorm(ngenes), each=ncells), ncol=ngenes) # gene-specific batch effect (genes are columns)
  A2 <- A2 + rnorm(ngenes*ncells) # noise
  rownames(A2) <- paste0("Cell", seq_len(ncells), "-2")
  colnames(A2) <- paste0("Gene", seq_len(ngenes))
  idx_responserCells_1 <- idx_responserCells[1:length(idx_responserCells)/3]
  idx_DEgenes_1 <- 1:60
  idx_responserCells_2 <- idx_responserCells[(length(idx_responserCells)/3+1):(2*length(idx_responserCells)/3)]
  idx_DEgenes_2 <- 61:120
  idx_responserCells_3 <- idx_responserCells[(2*length(idx_responserCells)/3+1):length(idx_responserCells)]
  idx_DEgenes_3 <- 121:175
  A2[idx_responserCells_1, idx_DEgenes_1] <- FC_DE*A2[idx_responserCells_1, idx_DEgenes_1] 
  A2[idx_responserCells_2, idx_DEgenes_2] <- FC_DE*A2[idx_responserCells_2, idx_DEgenes_2] 
  A2[idx_responserCells_3, idx_DEgenes_3] <- FC_DE*A2[idx_responserCells_3, idx_DEgenes_3]
  
  outmat$ii <- A2
  outcol$ii <- clust2
  B1 <- t(A1)
  B2ii <- t(outmat$ii)
  
  rowNames <- matrix(" ", dim(B1)[1], 1)
  for (i in 1:dim(rowNames)[1]){
    if (i %in% idx_DEgenes){
      rowNames[i] <- paste0("DE_gene", i)
    }
    else {
      rowNames[i] <- paste0("Gene", i)
    }
  }
  colNames1 <- matrix("Batch1", 1, dim(B1)[2])
  colNames2 <- matrix("Batch2", 1, dim(B2ii)[2])
  for (i in 1:dim(colNames2)[2]){
    if (i %in% idx_responserCells_1){
      colNames2[i] <- paste0("Res_1_", colNames2[i])
    }
  }
  for (i in 1:dim(colNames2)[2]){
    if (i %in% idx_responserCells_2){
      colNames2[i] <- paste0("Res_2_", colNames2[i])
    }
  }
  for (i in 1:dim(colNames2)[2]){
    if (i %in% idx_responserCells_3){
      colNames2[i] <- paste0("Res_3_", colNames2[i])
    }
  }
  colNames <- c(colNames1, colNames2)
  X <- cbind(B1, B2ii)
  rownames(X) <- rowNames
  colnames(X) <- colNames
  write.table(X, paste0(fileDir, paste0("/", paste0(paste0(itr, ".xls")))), sep="\t", row.names = TRUE, col.names = TRUE)
}
