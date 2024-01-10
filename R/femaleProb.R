#' femaleProb: probability based sex demultiplexing tool
#'
#' Calculates the probability of cells/nuclei belonging to male or female subjects
#' Returns a data frame with probabilitites based on 3 seperate models
#' Univariate mode,l dependent on Xist expression
#' Multivariate model, dependent on Xist and the sum of Y chromosome genes
#' Multivariate + nCount, same as above and including unique RNA counts
#'
#' @param Seuratobj seurat object
#' @param lognormalized boolean
#' @param ONLINE boolean
#' @param xistplots boolean
#' @export
#' @author Nickolas C. Chu


femaleProb <- function(Seuratobj, lognormalized = TRUE, ONLINE = TRUE, xistplots = FALSE )
{
  set.seed(1234)
  nCount_RNA <- NULL
  proportionF <- NULL
  Xist <- NULL
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
  }else{
    xist.exp <- expm1(FetchData(Seuratobj, "Xist"))
  }
  
  #sum up expression of all y chromosome genes. Unless ONLINE == FALSE, in this
  #case select genes will be used.
  data2 <- xist.exp
  
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
    
    #Checks for presence of cells belonging to cluster being checked.
    #If found subsets that cluster from the dataframe, otherwise checks next
    for (m in 1:max(levels)) {
      
      if (any (data2$seurat_clusters == current) ){
        newdata <- subset(data2, select = c("Xist", "Ygenes", "nCount_RNA"), data2$seurat_clusters == current )
        print("found")
        truelabels[i] <- current
        
        break
      }else{
        print('not found')
        current <- 1 + current
        print(paste("how about", current, "?"))
      }
      
    }
    
    #Produce density estimate for uni-variate and multi-variate
    if (NROW(newdata$Xist[newdata$Xist >= 2])){
      dens1 <- densityMclust(newdata$Xist, G = 2, plot = FALSE)#only Xist
      #multivariate probabilities
      dens3 <- densityMclust(newdata, G = 2, plot = FALSE)#Xist, ygenes, and rna count
      dens2 <- densityMclust(newdata[,c('Xist','Ygenes')], G = 2, plot = FALSE)#Xist and ygenes
      #Add density estimate object to a growing list and store cluster number
      
      plot_list0[[i-badclusters]]<- dens1
      plot_list0_titles[[i-badclusters]]<- current
      #Add probability of each class per cell to the dataframe that will be output
      newdata[['Prob.Uni.1']] <- dens1$z[,1]
      newdata[['Prob.Uni.2']] <- dens1$z[,2]
    }else{
      print(paste("less than 2 cells express Xist in cluster", current))
    }
    #not clear how mclust chooses to lable each class as 1 or 2, they can switch.
    #So to label the probabilities I use the correlation between probability and
    #Xist expression
    enoughxist <- TRUE
    #determine correlation of Xist expression to probability
    stat1 <- tryCatch(
      {stat1<- cor.test(newdata$Xist, newdata$Prob.Uni.1)
      },
      error=function(e) {
        message('An Error Occurred, check number of Xist expressing cells')
        print(e)
        enoughxist <- FALSE
        return(enoughxist)
        
      }
    )
    stat2 <- tryCatch(
      {stat2<- cor.test(newdata$Xist, newdata$Prob.Uni.2)
      },
      error=function(e) {
        message('An Error Occurred, check number of Xist expressing cells')
        print(e)
        enoughxist <- FALSE
        return(enoughxist)
        
      }
    )
    if (is.logical(stat1) || is.logical(stat2)){
      enoughxist <- FALSE
    }
    print(enoughxist)
    
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
      
      newdata[['Prob.Multi.1']] <- dens2$z[,1]
      newdata[['Prob.Multi.2']] <- dens2$z[,2]
      newdata[['Prob.Multi.ncount.1']] <- dens3$z[,1]
      newdata[['Prob.Multi.ncount.2']] <- dens3$z[,2]
      
      stat3<- cor.test(newdata$Xist, newdata$Prob.Multi.1)
      stat4<- cor.test(newdata$Xist, newdata$Prob.Multi.2)
      stat5<- cor.test(newdata$Xist, newdata$Prob.Multi.ncount.1)
      stat6<- cor.test(newdata$Xist, newdata$Prob.Multi.ncount.2)
      
      if (stat3$estimate > 0){
        colnames(newdata)[which(names(newdata) == "Prob.Multi.1")] <- "ProbFemaleMulti"
      }else{
        colnames(newdata)[which(names(newdata) == "Prob.Multi.1")] <- "ProbMaleMulti"
      }
      
      if (stat4$estimate > 0){
        colnames(newdata)[which(names(newdata) == "Prob.Multi.2")] <- "ProbFemaleMulti"
      }else{
        colnames(newdata)[which(names(newdata) == "Prob.Multi.2")] <- "ProbMaleMulti"
      }
      
      if (stat5$estimate > 0){
        colnames(newdata)[which(names(newdata) == "Prob.Multi.ncount.1")] <- "ProbFemaleMultinCount"
      }else{
        colnames(newdata)[which(names(newdata) == "Prob.Multi.ncount.1")] <- "ProbMaleMultinCount"
      }
      
      if (stat6$estimate > 0){
        colnames(newdata)[which(names(newdata) == "Prob.Multi.ncount.2")] <- "ProbFemaleMultinCount"
      }else{
        colnames(newdata)[which(names(newdata) == "Prob.Multi.ncount.2")] <- "ProbMaleMultinCount"
      }
      
    }else{#fill all columns with 0
      newdata[['Prob.Uni.1']] <- 0
      newdata[['Prob.Uni.2']] <- 0
      colnames(newdata)[which(names(newdata) == "Prob.Uni.1")] <- "ProbFemaleUni"
      colnames(newdata)[which(names(newdata) == "Prob.Uni.2")] <- "ProbMaleUni"
      newdata[['Prob.Multi.1']] <- 0
      newdata[['Prob.Multi.2']] <- 0
      newdata[['Prob.Multi.ncount.1']] <- 0
      newdata[['Prob.Multi.ncount.2']] <- 0
      colnames(newdata)[which(names(newdata) == "Prob.Multi.1")] <- "ProbFemaleMulti"
      colnames(newdata)[which(names(newdata) == "Prob.Multi.2")] <- "ProbMaleMulti"
      colnames(newdata)[which(names(newdata) == "Prob.Multi.ncount.1")] <- "ProbFemaleMultinCount"
      colnames(newdata)[which(names(newdata) == "Prob.Multi.ncount.2")] <- "ProbFemaleMultinCount"
      badclusters<-badclusters + 1
      print('not enough cells have xist')
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
    #add ordered dataframe into a list of dataframes. A dataframe per cluster
    Clusters <- append(Clusters, list(ordered),i)
    
  }
  #Create PDF for the change in female proportion with increasing rna count, per cluster
  if (xistplots == TRUE){
    print(truelabels)
    pdf(file = 'Proportion_Expressing_Xist.pdf')
    plot_list1 = list()
    noxist = 0
    for (N in 1:length(Clusters)){
      tempclust = Clusters[[N]]
      if(NROW(tempclust$Xist[tempclust$Xist >= 2])){
        hasxist<- N - noxist
        colnames(tempclust)[10] <- 'proportionF'
        p1 = ggplot(tempclust, aes(x = nCount_RNA, y = proportionF))+
          geom_point()+
          ggtitle(paste0("Cluster", truelabels[N]))
        plot_list1[[hasxist]] = p1
      }else{
        noxist <- noxist + 1
      }
    }
    for (i in 1:length(plot_list1)) {
      print(plot_list1[[i]])
    }
    dev.off()
    #create a list of histograms for Xist expression per cluster
    plot_list2 = list()
    for (N in 1:length(Clusters)){
      tempclust = Clusters[[N]]
      colnames(tempclust)[1] <- 'Xist'
      p2 = ggplot(tempclust, aes(x=Xist)) +
        geom_histogram(binwidth = 0.05) +
        #xlim(-1,NA)+
        #scale_x_continuous(limits = c(0,NA))+
        #scale_y_continuous()+
        ggtitle(paste0("Cluster", truelabels[N]))
      plot_list2[[N]] = p2
      
    }
    #create a PDF for the histograms for each cluster
    pdf(file = 'Xist_hist.pdf')
    for (i in 1:length(plot_list2)) {
      print(plot_list2[[i]])
    }
    dev.off()
    #Get number of clusters minus clusters with 2 < Xist expressing cells
    goodclusters<-  length(plot_list2) - badclusters
    #print PDF combining the ordered Xist expression proportions and Xist histograms
    pdf(file = 'combined.pdf')
    for(i in 1:goodclusters){
      current = plot_list0_titles[[i]]
      VAR1=plot_list0[[i]]
      VAR2=plot_list1[[i]]
      
      if(lognormalized == TRUE){
        bin_width <-0.05
      }
      else{
        bin_width <- 0.5
      }
      
      nbins <- seq(min(VAR1$data) - bin_width,
                   max(VAR1$data) + bin_width,
                   by = bin_width)
      
      plot(VAR2)
      plot(VAR1, what = "density", data =  VAR1$data, breaks = nbins , main = paste("Cluster", current, sep = " "), xlab = 'Xist' )
      
    }
    dev.off()
  }
  names(Clusters) <- truelabels
  merged <- bind_rows(Clusters, .id = "cluster")
  Clusters <- merged[rownames(Seuratobj@meta.data), ]
  
  return(Clusters)
}

