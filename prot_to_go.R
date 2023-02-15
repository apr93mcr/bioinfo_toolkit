#This function requires clusterProfiler package org.Hs.eg.db, and ggplot2
#prot_list must be a dataframe column or a vector containing a list of proteins or gene symbols
#From_type indicates thy type of gene symbols (SYMBOL) or uniprots ID ("UNIPROT") included in your protein list
#from_Type: possible option for from_Type: "UNIPROT","ENTREZID","ENSEMBL"
#simply_score: reduce the redundancy of GO terms (descending order), lower simply_score remove more redundant GO terms
#Min size evaluate the minimun of overlaping genes between you set and a given GO term
#When plot=T generates a ggplot barchar version depicting the GO terms and -log10(adjusted pvalue)

prot_to_go <- function(prot_list, from_Type="SYMBOL", min_size=5, simply_score=0.5, GO_onto=c("BP","MF","CC"), plot=F){
  
  if(!from_Type %in% c("UNIPROT","ENTREZID")){stop("Invalid key identifyer, use SYMBOL, UNIPROT or ENTREZID")}
  
  if(is.data.frame(prot_list)){prot_list <- prot_list[,1]}
  
  if(from_Type!="SYMBOL"){
    prot_list <- clusterProfiler::bitr(geneID = prot_list, fromType = "UNIPROT", toType = "SYMBOL", OrgDb="org.Hs.eg.db",drop = T)[2]
    prot_list <- prot_list[,1]}
  
  if(length(prot_list)==0){
    stop("The list of proteins is empty, check protein list format and type of gene or protein ID submitted")}
  
  res <- as.data.frame(matrix(nrow = 1, ncol = 10))
  colnames(res) <- c("ID","Description","GeneRatio","BgRatio","pvalue","p.adjust","qvalue","geneID","Count","Database")
  
  firstup <- function(x) {substr(x, 1, 1) <- toupper(substr(x, 1, 1)); x}
  
  for(i in GO_onto){
    cat(paste0("Searching in GO:",i,"\n"))
    
    test <- clusterProfiler::enrichGO(gene = prot_list, OrgDb = 'org.Hs.eg.db', keyType = "SYMBOL", ont = i)
    
    #Delete redundant GO terms
    if(dim(test@result)[1]!=0 & simply_score!=0){
      set.seed(12345)
      test <- clusterProfiler::simplify(test, cutoff=simply_score, by="p.adjust")
      test <- test@result}
      
    if(simply_score==0){test <- test@result}
    
    if(dim(test)[1]>0){
      test[,"Database"] <- paste0("GO:",i)
      res <- rbind(res,test)
      res <- res[res$Count >= min_size,]
      res$Description <- firstup(res$Description)
      }
  }
  
  if(dim(res)[1]<2){return("None GO terms were identified associated with your protein list")}
  
  if(dim(res)[1]>1 & plot==F){res <- res[-1,];  return(res)}
  
  if(dim(res)[1]>1 & plot==T){

    res <- res[-1,]
    res[,"log10pvalue"] <- -log10(res$p.adjust)
    
    library(ggplot2)
    p <- ggplot2::ggplot(res, aes(x = reorder(Description,log10pvalue), y = log10pvalue)) + ggplot2::ylab(expression("-log"[10]*"(adj p-value)")) + ggplot2::xlab("") +
      ggplot2::geom_point( aes(size= Count,colour = Database)) + ggplot2::scale_size(range=c(3,10)) + 
      ggplot2::coord_flip() + ggplot2::theme_bw() +
      ggplot2::theme(legend.title = element_blank(),
            axis.title = element_text(size = 15, face = "bold") ,
            axis.text.x = element_text(size = 15, face = "bold"),
            axis.text.y = element_text(size = 15, face = "bold"),
            legend.text = element_text(color = "black", size = 15),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank())
    
    print(p)
    return(res)
    }
}

#Runing example
#res_go <- prot_to_go(prot_list = protein_list, from_Type = "UNIPROT", simply_score = 0.5, plot = T, min_size = 5)
