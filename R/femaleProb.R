#' femaleProb: probability based sex demultiplexing tool
#'
#' Calculates the probability of cells/nuclei belonging to male or female subjects.
#' Returns a data frame with probabilities based on 4 separate models:
#' \itemize{
#'   \item Univariate model, dependent on Xist expression
#'   \item Multivariate model, dependent on Xist and the sum of Y chromosome genes
#'   \item Multivariate + nCount, same as above and including unique RNA counts
#'   \item Xchrom + Ychrom genes, dependent on sum of Xist and Tsix vs the sum of Y chromosome genes
#' }
#'
#' @param Seuratobj Seurat object
#' @param lognormalized boolean, whether input data is log-normalized
#' @param ONLINE boolean, whether to fetch Y chromosome genes online from EnsDb
#' @param xistplots boolean, whether to generate diagnostic plots
#' @param itmax integer, iterations for densityMclust fitting before skipping
#'
#' @return A data frame with probability estimates for each cell across the four models
#' @export
#' @author Nickolas C. Chu



femaleProb <- function(Seuratobj, lognormalized = TRUE, ONLINE = TRUE, itmax = 100 )
{
  set.seed(1234)
  nCount_RNA <- NULL
  proportionF <- NULL
  Xist <- NULL
  Tsix <- NULL
  #make empty list
  Clusters <- list()
  annotation<-NULL
  enoughxist <- TRUE
  seurat_version=as.numeric(substring(Seuratobj@version, 1,2))
  
  #for catching clusters with low numbers of Xist expressing cells
  badclusters<-0
  
  #get number of clusters and Xist expression for all cells
  levels <- as.data.frame(as.numeric(as.character(unique(Seuratobj@active.ident))))
  if (lognormalized == TRUE){
    xist.exp <-FetchData(Seuratobj, "Xist")
    tsix.exp <- FetchData(Seuratobj, "Tsix")
  }else{
    xist.exp <- expm1(FetchData(Seuratobj, "Xist"))
    tsix.exp <- expm1(FetchData(Seuratobj, "Tsix"))
  }
  
  #sum up expression of all y chromosome genes. Unless ONLINE == FALSE, in this
  #case select genes will be used.
  data2 <- xist.exp
  data2[["Tsix"]] <- tsix.exp$Tsix
  
  if(seurat_version >= 5){
    if(ONLINE == TRUE){
      annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
      seqlevelsStyle(annotation) <- "UCSC"
      Ychrom_genes <- unique(annotation[annotation@seqnames == "chrY",]$gene_name)
      Seuratobj@assays$RNA$data[1:5, 1:5]
      Ychr_count <- as.matrix(Seuratobj@assays$RNA$data[rownames(Seuratobj@assays$RNA$data) %in% Ychrom_genes, ])
      transposed <- t(Ychr_count)
      Ygenes <- as.data.frame(rowSums(transposed))
    }else{
      Ygenes <- FetchData(Seuratobj,"Uty") + FetchData(Seuratobj, "Eif2s3y") + FetchData(Seuratobj, "Ddx3y")
    }
  }else{
    if(ONLINE == TRUE){
      annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
      seqlevelsStyle(annotation) <- "UCSC"
      Ychrom_genes <- unique(annotation[annotation@seqnames == "chrY",]$gene_name)
      Seuratobj@assays$RNA@data[1:5, 1:5]
      Ychr_count <- as.matrix(Seuratobj@assays$RNA@data[rownames(Seuratobj@assays$RNA@data) %in% Ychrom_genes, ])
      transposed <- t(Ychr_count)
      Ygenes <- as.data.frame(rowSums(transposed))
    }else{
      Ygenes <- FetchData(Seuratobj,"Uty") + FetchData(Seuratobj, "Eif2s3y") + FetchData(Seuratobj, "Ddx3y")
    }
  }
  
  
  
  #make a dataframe that has Xist, Ygenes, and cluster number
  data2[,"Ygenes"] <- Ygenes
  data2$Xgenes<- data2$Tsix + data2$Xist
  data2[,"seurat_clusters"] <- Seuratobj@active.ident
  data2[['nCount_RNA']] <- Seuratobj$nCount_RNA
  
  #preparing lists for graphs
  truelabels <- vector("character", nrow(levels))
  plot_list0 = list()
  plot_list0_titles = list()
  
  #Main loop function loops for the number of unique identities in the dataframe
  
  previous <- 0
  current <- -1
  for (i in 1:nrow(levels)) {
    print(paste("loop", i))
    current <- 1 + current
    print (paste("checking for cluster", current))
    
    found <- FALSE
    invalid <- FALSE
    newdata <- NULL
    constant_x <-FALSE
    #Checks for presence of cells belonging to cluster being checked.
    #If found subsets that cluster from the dataframe, otherwise checks next
    for (m in 1:max(levels)) {
      
      if (any (data2$seurat_clusters == current) ){
        newdata <- subset(data2,
                          subset = data2$seurat_clusters == current, 
                          select = c("Xist", "Ygenes","Xgenes", "nCount_RNA")
        )
        print("found")
        truelabels[i] <- current
        found<-TRUE
        if (nrow(newdata) < 5) {
          message("Skipping cluster ", current, " (<5 cells)")
          invalid<-TRUE
          break
        }
        if (n_distinct(newdata$Xist) == 1) {
          message("Skipping cluster ", current, " (Xist constant)")
          constant_x <-TRUE
        }
        break
      }else{
        print('not found')
        current <- 1 + current
        print(paste("how about", current, "?"))
      }
      
    }

    if (!found) next
    if (invalid) { badclusters <- badclusters + 1; next }
    
    #Produce density estimate for uni-variate and multi-variate

    
    if (NROW(newdata$Ygenes[newdata$Ygenes > 0]) >= 2) {
      densY <- tryCatch({
        densityMclust(newdata$Ygenes, G = 2, plot = FALSE,
                      control = emControl(itmax = itmax))
      }, error = function(e) {
        message("Cluster ", current, ": densY failed - ", e$message)
        NULL
      }) 
      if (!is.null(densY)) {
        newdata[["Prob.UniY.1"]] <- densY$z[,1]
        newdata[["Prob.UniY.2"]] <- densY$z[,2]
      }
    }

    if (NROW(newdata$Xist[newdata$Xist > 0]) >= 2 & constant_x == FALSE) {
      dens1 <- tryCatch({
        densityMclust(newdata$Xist, G = 2, plot = FALSE,
                      control = emControl(itmax = itmax))
      }, error = function(e) {
        message("Cluster ", current, ": dens1 failed - ", e$message)
        NULL
      }) # only Xist probabilities
    
      dens3 <- tryCatch({
        densityMclust(newdata[,c("Xist","Ygenes","nCount_RNA")], G = 2, plot = FALSE,
                      control = emControl(itmax = itmax))
      }, error = function(e) {
        message("Cluster ", current, ": dens3 failed - ", e$message)
        NULL
      }) # Xist, ygenes, and rna count
    
      dens2 <- tryCatch({
        densityMclust(newdata[,c("Xist","Ygenes")], G = 2, plot = FALSE,
                      control = emControl(itmax = itmax))
      }, error = function(e) {
        message("Cluster ", current, ": dens2 failed - ", e$message)
        NULL
      }) # Xist and ygenes
    
      dens4 <- tryCatch({
        densityMclust(newdata[,c("Xgenes","Ygenes")], G = 2, plot = FALSE,
                      control = emControl(itmax = itmax))
      }, error = function(e) {
        message("Cluster ", current, ": dens4 failed - ", e$message)
        NULL
      })

    
      #Add density estimate object to a growing list and store cluster number
      # Only proceed if dens1 succeeded
      if (!is.null(dens1)) {
        plot_list0[[i-badclusters]] <- dens1
        plot_list0_titles[[i-badclusters]] <- current
        newdata[["Prob.Uni.1"]] <- dens1$z[,1]
        newdata[["Prob.Uni.2"]] <- dens1$z[,2]
      }
    }else{
      print(paste("less than 2 cells express Xist in cluster", current))
      badclusters <- badclusters + 1
    }
    enoughYgenes <- TRUE
    if(!is.null(densY)) {
      stat1Y <- tryCatch(
        {stat1<- cor.test(newdata$Ygenes, newdata$Prob.UniY.1)
        },
        error=function(e) {
          message('An Error Occurred, check number of Y chrom expressing cells')
          print(e)
          enoughYgenes <- FALSE
          return(enoughYgenes)
          
        }
      )
      stat2Y <- tryCatch(
        {stat2<- cor.test(newdata$Ygenes, newdata$Prob.UniY.2)
        },
        error=function(e) {
          message('An Error Occurred, check number of Y chrom expressing cells')
          print(e)
          enoughYgenes <- FALSE
          return(enoughYgenes)
          
        }
      )

    if (is.logical(stat1Y) || is.logical(stat2Y)){
      enoughYgenes <- FALSE
    }

    }
    #not clear how mclust chooses to lable each class as 1 or 2, they can switch.
    #So to label the probabilities I use the correlation between probability and
    #Xist expression
    enoughxist <- TRUE
    
    # Only run cor.test if the probability columns exist and are numeric
    stat1 <- if ("Prob.Uni.1" %in% names(newdata) && is.numeric(newdata$Prob.Uni.1)) {
      tryCatch(cor.test(newdata$Xist, newdata$Prob.Uni.1),
               error = function(e) {
                 message("Cluster ", current, ": cor.test on Prob.Uni.1 failed - ", e$message)
                 NULL
               })
    } else NULL
    
    stat2 <- if ("Prob.Uni.2" %in% names(newdata) && is.numeric(newdata$Prob.Uni.2)) {
      tryCatch(cor.test(newdata$Xist, newdata$Prob.Uni.2),
               error = function(e) {
                 message("Cluster ", current, ": cor.test on Prob.Uni.2 failed - ", e$message)
                 NULL
               })
    } else NULL
    
    # If either correlation failed, mark as not enough Xist
    if (is.null(stat1) || is.null(stat2)) {
      enoughxist <- FALSE
    }

    
    if (enoughYgenes == TRUE) {
      #A positive correlation indicates the probability is for being female,
      #negative is male
      print(paste("there is enough Y chrom expression in cluster:", current))
      #univariate probabilities
      if (stat1Y$estimate > 0){
        colnames(newdata)[which(names(newdata) == "Prob.UniY.1")] <- "ProbMaleUniY"
      }else{
        colnames(newdata)[which(names(newdata) == "Prob.UniY.1")] <- "ProbFemaleUniY"
      }
      
      if (stat2Y$estimate > 0){
        colnames(newdata)[which(names(newdata) == "Prob.UniY.2")] <- "ProbMaleUniY"
      }else{
        colnames(newdata)[which(names(newdata) == "Prob.UniY.2")] <- "ProbFemaleUniY"
      }
      
    }
    
    
    if (enoughxist == TRUE) {
      #A positive correlation indicates the probability is for being female,
      #negative is male
      print(paste("there is enough Xist in cluster:", current))
      #univariate probabilities
      if (stat1$estimate > 0){
        colnames(newdata)[which(names(newdata) == "Prob.Uni.1")] <- "ProbFemaleUni"
      }else{
        colnames(newdata)[which(names(newdata) == "Prob.Uni.1")] <- "ProbMaleUni"
      }
      
      if (stat2$estimate > 0){
        colnames(newdata)[which(names(newdata) == "Prob.Uni.2")] <- "ProbFemaleUni"
      }else{
        colnames(newdata)[which(names(newdata) == "Prob.Uni.2")] <- "ProbMaleUni"
      }
      
      if (!is.null(dens2)) {
        newdata[['Prob.Multi.1']] <- dens2$z[,1]
        newdata[['Prob.Multi.2']] <- dens2$z[,2]
      } else {
        newdata[['Prob.Multi.1']] <- 0
        newdata[['Prob.Multi.2']] <- 0
      }
      if (!is.null(dens3)) {
        newdata[['Prob.Multi.ncount.1']] <- dens3$z[,1]
        newdata[['Prob.Multi.ncount.2']] <- dens3$z[,2]
      } else {
        newdata[['Prob.Multi.ncount.1']] <- 0
        newdata[['Prob.Multi.ncount.2']] <- 0
      }
        if (!is.null(dens4)) {
        newdata[["Prob.XY.1"]] <- dens4$z[,1]
        newdata[["Prob.XY.2"]] <- dens4$z[,2]
      } else {
        newdata[["Prob.XY.1"]] <- 0
        newdata[["Prob.XY.2"]] <- 0
      }
      
      stat3 <- tryCatch(cor.test(newdata$Xist, newdata$Prob.Multi.1), error=function(e) NULL)
      stat4 <- tryCatch(cor.test(newdata$Xist, newdata$Prob.Multi.2), error=function(e) NULL)
      stat5 <- tryCatch(cor.test(newdata$Xist, newdata$Prob.Multi.ncount.1), error=function(e) NULL)
      stat6 <- tryCatch(cor.test(newdata$Xist, newdata$Prob.Multi.ncount.2), error=function(e) NULL)
      stat7 <- tryCatch(cor.test(newdata$Xist, newdata$Prob.XY.1), error=function(e) NULL)
      stat8 <- tryCatch(cor.test(newdata$Xist, newdata$Prob.XY.2), error=function(e) NULL)
      
      if (!is.null(stat3) && stat3$estimate > 0){
        colnames(newdata)[which(names(newdata) == "Prob.Multi.1")] <- "ProbFemaleMulti"
      } else {
        colnames(newdata)[which(names(newdata) == "Prob.Multi.1")] <- "ProbMaleMulti"
      }
      if (!is.null(stat4) && stat4$estimate > 0){
        colnames(newdata)[which(names(newdata) == "Prob.Multi.2")] <- "ProbFemaleMulti"
      } else {
        colnames(newdata)[which(names(newdata) == "Prob.Multi.2")] <- "ProbMaleMulti"
      }
      if (!is.null(stat5) && stat5$estimate > 0){
        colnames(newdata)[which(names(newdata) == "Prob.Multi.ncount.1")] <- "ProbFemaleMultinCount"
      } else {
        colnames(newdata)[which(names(newdata) == "Prob.Multi.ncount.1")] <- "ProbMaleMultinCount"
      }
      if (!is.null(stat6) && stat6$estimate > 0){
        colnames(newdata)[which(names(newdata) == "Prob.Multi.ncount.2")] <- "ProbFemaleMultinCount"
      } else {
        colnames(newdata)[which(names(newdata) == "Prob.Multi.ncount.2")] <- "ProbMaleMultinCount"
      }
      if (!is.null(stat7) && stat7$estimate > 0) {
        colnames(newdata)[which(names(newdata) == "Prob.XY.1")] <- "ProbFemaleXY"
      } else {
        colnames(newdata)[which(names(newdata) == "Prob.XY.1")] <- "ProbMaleXY"
      }
      if (!is.null(stat8) && stat8$estimate > 0) {
        colnames(newdata)[which(names(newdata) == "Prob.XY.2")] <- "ProbFemaleXY"
      } else {
        colnames(newdata)[which(names(newdata) == "Prob.XY.2")] <- "ProbMaleXY"
      }
    }

    if(enoughxist == TRUE |enoughYgenes == TRUE){
    }else{
      print('not enough cells have xist or fits failed')
      break
    }
    
    
    
    #Determine at which value of ncount, the proportion of female cells stops changing
    max = 0
    female = 0
    femalelist = list("femalenum")
    newdata[["proportionF"]]<-0
    # ordere the dataframe by nCount_RNA
    ordered <-arrange(newdata,nCount_RNA)
    #starting from the lowest ncount keep track of the proportion of female cells,
    #where any cell expressing Xist is considered female. Can change this later
    for (k in 1:length(newdata$Xist)){
      cell = ordered$Xist[k]
      if (cell > 0 ){
        female = female + 1
      }else{
      }
      proportion = female/k
      ordered$proportionF[[k]] = proportion
    }
    # Must be a valid cluster and have at least one probability column present
    valid_prob_cols <- c("ProbFemaleUni","ProbMaleUni","ProbFemaleMulti","ProbMaleMulti",
                         "ProbFemaleMultinCount","ProbMaleMultinCount","ProbFemaleXY","ProbMaleXY", "ProbMaleUniY", "ProbFemaleUniY")
    
    has_probs <- any(valid_prob_cols %in% names(ordered))
    newdata$cell_id <- rownames(newdata)
    ordered$cell_id <- rownames(ordered)
                        
    if (!invalid && (enoughxist || enoughYgenes) && has_probs) {
      Clusters[[as.character(current)]] <- ordered
    }
  }

                        
 

  # Remove NULLs
  Clusters <- Filter(Negate(is.null), Clusters)
  
  # Remove data frames that are effectively empty (all NA in required fields)
  required_core <- c("Xist","Ygenes","Xgenes","nCount_RNA")
  Clusters <- Filter(function(df) {
    is.data.frame(df) &&
      nrow(df) > 0 &&
      all(required_core %in% names(df)) &&
      # not all NA across core columns
      !all(rowSums(is.na(df[required_core])) == length(required_core))
  }, Clusters)

                        
  merged <- bind_rows(Clusters, .id = "cluster")
  
  
  return(merged)
}

