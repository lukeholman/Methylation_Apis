make_cor_plot <- function(df, diverging, reorder = TRUE){
  
  df <- as.data.frame(df) %>%
    select_if(is.numeric)
  
  # function to reorder a correlation matrix using hierarchical clustering
  reorder_cormat <- function(cormat){ 
    dd <- as.dist((1-cormat)/2)
    hc <- hclust(dd)
    cormat[hc$order, hc$order]}
  
  cormat <- cor(df, use = 'pairwise.complete.obs', method = "spearman") 
  if(reorder) cormat <- cormat %>% reorder_cormat()
  cormat[upper.tri(cormat)] <- NA
  diag(cormat) <- NA
  
  cor.matrix <-  (reshape2::melt(cormat) %>% filter(!is.na(value)))[,1:2] %>% mutate(cor=0,p=0) 
  for(i in 1:nrow(cor.matrix)) { 
    cor.matrix[i, 3:4] <- cor.test(df[, names(df) == cor.matrix$Var1[i]], 
                                   df[, names(df) == cor.matrix$Var2[i]],
                                   use = 'pairwise.complete.obs',
                                   method = "spearman")[c(4,3)] %>% unlist()
  }
  
  cor.matrix$p <- p.adjust(cor.matrix$p, method = "BH") # apply B-H p value correction
  cor.matrix$sig <- ""
  cor.matrix$sig[cor.matrix$p < 0.05] <- "*"   
  cor.matrix$sig[cor.matrix$p < 0.001] <- "**"  
  cor.matrix$sig[cor.matrix$p < 0.0001] <- "***"   
  cor.matrix$label <- paste(format(round(cor.matrix$cor,2), nSmall =2), cor.matrix$sig, sep = "")
  brew <- brewer.pal(7, "Purples")
  
  p <- cor.matrix %>% 
    ggplot(aes(Var1, Var2, fill=cor)) + 
    geom_tile(colour = "grey10", size = 0.4, linetype = 3) + 
    geom_text(aes(label = label), colour = "black", size = 3.5) + 
    xlab(NULL) + ylab(NULL) +
    scale_x_discrete(expand=c(0.01,0.01)) + scale_y_discrete(expand = c(0.01, 0.01)) + 
    theme_minimal() + theme(panel.grid = element_blank(), 
                            legend.position = "none") 
  
  if(!diverging) return(p + scale_fill_gradient(low = brew[1], high = brew[4]))
  p + scale_fill_gradient2(low = "steelblue", high = "tomato")
}