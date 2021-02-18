make_pathway_figure <- function(GO_output, ylab, xlab, alt_colours = FALSE){
  
  if(nrow(GO_output) == 0) {
    warning("No significant GO terms to plot")
    return(NULL)
  }
  
  format_GO_results <- function(GO_output){
    levels_order <- GO_output %>%
      mutate(term_wrapped = str_wrap(term, width = 60)) %>%
      group_by(term_wrapped) %>%
      #summarise(p = prod(padj)^0.25, .groups = "drop") %>%
      summarise(NES_mean = mean(NES), .groups = "drop") %>%
      arrange(NES_mean) %>% pull(term_wrapped) 
    
    GO_output %>% as_tibble() %>%
      mutate(term = str_wrap(term, width = 60)) %>%
      mutate(term = factor(term, levels_order),
             sig = ifelse(pval < 0.05, "*", ""),
             sig = replace(sig, padj < 0.05, "**"))
  }
  
  p <- format_GO_results(GO_output) %>%
    mutate(Time = as.character(Time)) %>%
    mutate(Test_type = replace(Test_type, Test_type == "GO: Biological process", "GO: BP"),
           Test_type = replace(Test_type, Test_type == "GO: Molecular function", "GO: MF"),
           Test_type = replace(Test_type, Test_type == "GO: Cellular component", "GO: CC")) %>%
    ggplot(aes(Time, term, fill = NES)) + 
    geom_tile(colour = "black", size = 0.4)  + 
    geom_text(aes(label = sig)) +
    theme_bw() + ylab(ylab) + xlab(xlab) +
    scale_x_discrete(expand = c(0,0)) + 
    scale_y_discrete(expand = c(0,0)) +
    facet_grid(Test_type ~ ., scales = "free", space = "free") + 
    theme(strip.background = element_blank())
  
  if(!alt_colours){
    col_low <- brewer.pal(9, "Blues")[7]
    col_high <-  brewer.pal(9, "Reds")[7]
    p <- p + scale_fill_gradient2(low = col_low, high = col_high) 
  } else {
    col_low <- brewer.pal(9, "Purples")[7]
    col_high <-  brewer.pal(9, "Oranges")[7]
    p <- p + scale_fill_gradient2(low = col_low, high = col_high) 
   # p <- p + scale_fill_distiller(palette = "Purples", direction = 1) 
  }
  p
}
