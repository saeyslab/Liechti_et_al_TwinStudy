compareSubset <- function(All_traits, subset, names){
  
  subset_selected <- ggplot(All_traits %>% dplyr::filter(Group != "QC")) +
    geom_point(aes(x = umap_1, y = umap_2, col = names[subset+1][All_traits$Group != "QC"])) +
    geom_point(aes(x = umap_1, y = umap_2), col = "black",
               data = All_traits %>% dplyr::filter(Group == "QC")) +
    theme_void() +
    scale_color_discrete(name = "") +
    theme(legend.position="bottom")
  
  diff <- sapply(grep("% *MC", colnames(All_traits), value = TRUE),
                 function(trait){
                   test_res <- wilcox.test(All_traits[subset,trait],
                                           All_traits[!subset,trait])
                   fold_change <- median(All_traits[subset,trait])/
                     median(All_traits[!subset, trait])
                   
                   return(c(p_value = test_res$p.value, 
                            #`-log10p` = -log10(test_res$p.value),
                            fold_change = fold_change,
                            logfoldchange = log10(fold_change)))
                 }) %>% 
    t() %>% 
    data.frame(check.names = FALSE) %>% 
    tibble::rownames_to_column("Trait")
  
  diff$p_corrected <- p.adjust(diff$p_value, method = "BH")
  diff$`-log10p` <- -log10(diff$p_corrected)
  
  # max_fold <- max(abs(diff[,"logfoldchange"]))
  # max_p <- max(diff[,"-log10p"][!is.infinite(diff[,"-log10p"])], na.rm = TRUE)
  # print(max_p)
  # subset_volcano <- ggplot(diff) +
  #   geom_text(aes(x = logfoldchange, 
  #                 y = `-log10p`, 
  #                 label = Trait),
  #             size = 2) +
  #   xlim(-1.5*max_fold,1.5*max_fold)+
  #   ylim(0, 1.05*max_p)+
  #   theme_minimal()
  # 
  # best_feature <- diff[which.min(diff[,"p_value"]), "Trait"]
  # 
  # subset_bestFeature <- ggplot(All_traits) +
  #   geom_point(aes_string(x = "subset", y = paste0("`", best_feature, "`*100")), 
  #              position = ggbeeswarm::position_quasirandom(width = 0.3),
  #              size = 1) +
  #   geom_boxplot(aes_string(x = "subset", 
  #                           y = paste0("`", best_feature, "`*100"),
  #                           col = "subset"),
  #                outlier.alpha = 0,
  #                alpha = 0) +
  #   theme_minimal() +
  #   ylab(best_feature)
  # 
  # subset_overview <- subset_selected + subset_volcano + subset_bestFeature
  # 
  best_features <- diff[head(order(diff$`-log10p`,
                                   decreasing = TRUE), n = 4), "Trait"]
  subset_bestFeature <- list()
  
  for(best_feature in best_features){
    All_traits[, best_feature] <- All_traits[, best_feature] * 100
    subset_bestFeature[[best_feature]] <- ggplot(All_traits) +
      geom_point(aes_string(x = "names[subset+1]", y = paste0("`", best_feature, "`")), 
                 position = ggbeeswarm::position_quasirandom(width = 0.3),
                 size = 1) +
      geom_boxplot(aes_string(x = "names[subset+1]", 
                              y = paste0("`", best_feature, "`"),
                              col = "names[subset+1]"),
                   outlier.alpha = 0,
                   alpha = 0) +
      xlab("") +
      theme_minimal() +
      ylab(best_feature) + 
      theme(legend.position = "none")
  }
  
  best_correlated <- diff[head(order(diff$`-log10p`,
                                    decreasing = TRUE), n  = 10),]

  tnk_sub <- best_correlated[grep("TNK", best_correlated$Trait),]
  rownames(tnk_sub) <- gsub("TNK %MC([0-9]*)", "TNK MC\\1", tnk_sub$Trait)
  tnk_phenotype <- tnk_mc_mfi[rownames(tnk_sub),]
  tnk_heatmap <- pheno_heatmap(tnk_phenotype, tnk_sub[,c("logfoldchange", "-log10p")])
  
  bdc_sub <- best_correlated[grep("BDC", best_correlated$Trait),]
  rownames(bdc_sub) <- gsub("BDC %MC([0-9]*)", "BDC MC\\1", bdc_sub$Trait)
  bdc_phenotype <- bdc_mc_mfi[rownames(bdc_sub),]
  bdc_heatmap <- pheno_heatmap(bdc_phenotype, bdc_sub[,c("logfoldchange", "-log10p")])
  
  subset_overview <- #subset_volcano + 
    subset_selected | 
    ((subset_bestFeature[[1]] + subset_bestFeature[[2]]) / 
       (subset_bestFeature[[3]] +  subset_bestFeature[[4]]) + plot_layout(tag_level = 'new'))  |
    wrap_plots(grid::grid.grabExpr(ComplexHeatmap::draw(tnk_heatmap)),
               #grid::grid.grabExpr(ComplexHeatmap::draw(bdc_heatmap)), # For actual used example no BDC features
               ncol = 1)  + plot_layout(widths = c(1, 1, 1)) + plot_layout(tag_level = 'new')
  return(list(comparison = diff,
              plot = subset_overview))
}


compareSubsetAge <- function(All_traits, subset){
  
  subset_selected <- ggplot(All_traits) +
    geom_point(aes(x = umap_1, y = umap_2, col = subset)) +
    theme_void()
  
  cor <- sapply(grep("% *MC", colnames(All_traits), value = TRUE),
                 function(trait){
                   test_res <- cor.test(All_traits[subset,trait],
                                        All_traits[subset,"Age"],
                                        use = "pairwise.complete.obs",
                                        method = "spearman")
                   return(c(test_res$estimate["rho"],
                            p_value = test_res$p.value))#,
                           # "-log10p" = -log10(test_res$p.value)
                 }) %>% 
    t() %>% 
    data.frame(check.names = FALSE) %>% 
    tibble::rownames_to_column("Trait")
  
  cor$p_corrected <- p.adjust(cor$p_value, method = "BH")
  cor$`-log10p` <- -log10(cor$p_corrected)
  
  max_rho <- max(abs(cor[,"rho"]))
  max_p <- max(cor[,"-log10p"][!is.infinite(cor[,"-log10p"])], na.rm = TRUE)
  print(max_p)
  subset_volcano <- ggplot(cor) +
    geom_text(aes(x = rho, 
                  y = `-log10p`, 
                  label = Trait),
              size = 2) +
    xlim(-1.5*max_rho,1.5*max_rho)+
    ylim(0, 1.1* max_p)+
    theme_minimal()
  
  best_features <- cor[head(order(cor$`-log10p`,
                                  decreasing = TRUE), n = 4), "Trait"]
  print(best_features)
  subset_bestFeature <- list()
  
  for(best_feature in best_features){
    All_traits[, best_feature] <- All_traits[, best_feature] * 100
    subset_bestFeature[[best_feature]] <- ggpubr::ggscatter(All_traits[subset,],
                                                          x = "Age",
                                                          y = best_feature,
                                                          add = "reg.line", 
                                                          add.params = list(color = "#00BFC4"),
                                                          conf.int = TRUE, 
                                                          cor.coef = TRUE,
                                                          cor.method = "spearman") 
  }
  
  best_correlated <- cor[head(order(cor$`-log10p`,
                                    decreasing = TRUE), n  = 10),]
  
  age_tnk_sub <- best_correlated[grep("TNK", best_correlated$Trait),]
  rownames(age_tnk_sub) <- gsub("TNK %MC([0-9]*)", "TNK MC\\1", age_tnk_sub$Trait)
  age_tnk_phenotype <- tnk_mc_mfi[rownames(age_tnk_sub),]
  age_tnk_heatmap <- pheno_heatmap(age_tnk_phenotype, age_tnk_sub[,c("rho", "-log10p")])

  age_bdc_sub <- best_correlated[grep("BDC", best_correlated$Trait),]
  rownames(age_bdc_sub) <- gsub("BDC %MC([0-9]*)", "BDC MC\\1", age_bdc_sub$Trait)
  age_bdc_phenotype <- bdc_mc_mfi[rownames(age_bdc_sub),]
  age_bdc_heatmap <- pheno_heatmap(age_bdc_phenotype, age_bdc_sub[,c("rho", "-log10p")])
  
  subset_overview <- #subset_volcano + 
    ((subset_bestFeature[[best_features[1]]] + subset_bestFeature[[best_features[2]]]) / 
       (subset_bestFeature[[best_features[3]]] +  subset_bestFeature[[best_features[4]]]) + plot_layout(tag_level = 'new'))  |
    wrap_plots(grid::grid.grabExpr(ComplexHeatmap::draw(age_tnk_heatmap)),
               grid::grid.grabExpr(ComplexHeatmap::draw(age_bdc_heatmap)),
               ncol = 1) + plot_layout(tag_level = 'new')
  
  return(list(correlation = cor,
              plot = subset_overview))
}

pheno_heatmap <- function(phenotype, annotation){
  ComplexHeatmap::Heatmap(as.matrix(phenotype),
                          col = rev(RColorBrewer::brewer.pal(8,"RdBu")),
                          cluster_columns = FALSE,
                          cluster_rows = FALSE,
                          row_names_side = "left",
                          right_annotation = ComplexHeatmap::rowAnnotation(
                            df = annotation,
                            col = list("rho" = circlize::colorRamp2(c(-1, 0, 1),
                                                                    c("#7b3294", 
                                                                      "white", 
                                                                      "#008837")),
                                       "logfoldchange" = circlize::colorRamp2(c(-1, 0, 1),
                                                                    c("#7b3294", 
                                                                      "white", 
                                                                      "#008837")),
                                       "-log10p" = circlize::colorRamp2(c(30, 0),
                                                                        c("#008837",
                                                                          "white")))))
}