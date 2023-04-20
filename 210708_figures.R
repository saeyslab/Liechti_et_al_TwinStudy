# Load libraries and helper functions ------------------------------------------

library(flowCore)
library(FlowSOM)
library(dplyr)
library(ggplot2)
library(patchwork)
source("210708_loadDemographics.R")

# Set up variables -------------------------------------------------------------

recompute <- FALSE
date <- "210708"
#base_dir <- "/auto/net/fs4-Flowcytometry/Storage/u_ysa/Projects/twin_study/52e4c244-0e9b-4573-8e11-d10c249a871c-twinstudy2/"
base_dir <- "D:/Git_repos/twin_study"

demographics <- loadDemographics(base_dir)
demographics <- dplyr::filter(demographics, `ICS_Viable of CD3` > 25)

# Load results -----------------------------------------------------------------

## TNK 

fsom_TNK <- readRDS(file.path(base_dir, "RDS", paste0(date,"_TNK_FlowSOM.RDS")))

result_files_TNK <- list.files(file.path(base_dir, "RDS"), 
                               pattern = paste0(date, "_TNK_Results"),
                               full.names = TRUE) %>% 
  stringr::str_sort(numeric = TRUE)

results_TNK <- lapply(result_files_TNK, readRDS)

MC_pctgs_TNK <- do.call(rbind, lapply(results_TNK, function(x) x$metacluster_percentages))
colnames(MC_pctgs_TNK) <- paste0("TNK ", colnames(MC_pctgs_TNK))

MC_counts_TNK <- do.call(rbind, lapply(results_TNK, function(x) x$metacluster_counts))
colnames(MC_counts_TNK) <- paste0("TNK ", colnames(MC_counts_TNK))

MC_MFI_TNK <- do.call(rbind, lapply(results_TNK, function(x) x$metacluster_MFIs))
colnames(MC_MFI_TNK) <- paste0("TNK ", colnames(MC_MFI_TNK))

MC_outliers_TNK <- do.call(c, lapply(results_TNK, function(x){
  rowSums(attr(x$cluster_counts, "outliers"))
}))
MC_outliers_TNK <- MC_outliers_TNK / rowSums(MC_counts_TNK)

## BDC

fsom_BDC <- readRDS(file.path(base_dir, "RDS", paste0(date,"_BDC_FlowSOM.RDS")))

result_files_BDC <- list.files(file.path(base_dir, "RDS"), 
                               pattern = paste0(date, "_BDC_Results_[0-9]"),
                               full.names = TRUE) %>% 
  stringr::str_sort(numeric = TRUE)

results_BDC <- lapply(result_files_BDC, readRDS)

MC_pctgs_BDC <- do.call(rbind, lapply(results_BDC, function(x) x$metacluster_percentages))
colnames(MC_pctgs_BDC) <- paste0("BDC ", colnames(MC_pctgs_BDC))

MC_counts_BDC <- do.call(rbind, lapply(results_BDC, function(x) x$metacluster_counts))
colnames(MC_counts_BDC) <- paste0("BDC ", colnames(MC_counts_BDC))


MC_outliers_BDC <- do.call(c, lapply(results_BDC, function(x){
  rowSums(attr(x$cluster_counts, "outliers"))
}))
MC_outliers_BDC <- MC_outliers_BDC / rowSums(MC_counts_BDC)

## BDC normalised
result_files_BDC_norm <- list.files(file.path(base_dir, "RDS"), 
                               pattern = paste0(date, "_BDC_Results_normalized"),
                               full.names = TRUE) %>% 
  stringr::str_sort(numeric = TRUE)

results_BDC_norm <- lapply(result_files_BDC_norm, readRDS)

MC_pctgs_BDC_norm <- do.call(rbind, lapply(results_BDC_norm, function(x) x$metacluster_percentages))
colnames(MC_pctgs_BDC_norm) <- paste0("BDC ", colnames(MC_pctgs_BDC_norm))

MC_counts_BDC_norm <- do.call(rbind, lapply(results_BDC_norm, function(x) x$metacluster_counts))
colnames(MC_counts_BDC_norm) <- paste0("BDC ", colnames(MC_counts_BDC_norm))

MC_MFI_BDC_norm <- do.call(rbind, lapply(results_BDC_norm, function(x) x$metacluster_MFIs))
colnames(MC_MFI_BDC_norm) <- paste0("BDC ", colnames(MC_MFI_BDC_norm))


MC_outliers_BDC_norm <- do.call(c, lapply(results_BDC_norm, function(x){
  rowSums(attr(x$cluster_counts, "outliers"))
}))
MC_outliers_BDC_norm <- MC_outliers_BDC_norm / rowSums(MC_counts_BDC_norm)

# Combine into one data frame --------------------------------------------------

complete <- intersect(rownames(MC_pctgs_BDC_norm), rownames(MC_pctgs_TNK))
All_traits <- data.frame(demographics[complete,],
                         MC_pctgs_BDC_norm[complete,],
                         MC_pctgs_TNK[complete,],
                         check.names = FALSE)

counts <- data.frame(demographics[complete,],
                     MC_counts_BDC_norm[complete,],
                     MC_counts_TNK[complete,],
                     check.names = FALSE)
                     

date <- "210909"

# Figures FlowSOM model --------------------------------------------------------

for (panel in c("TNK", "BDC")) {
  if (panel == "TNK") {
    fsom <- fsom_TNK
  } else if (panel == "BDC"){
    fsom <- fsom_BDC
  }
  
  # FlowSOMmary
  FlowSOMmary(fsom, paste0(date, "_", panel, "_FlowSOM.pdf"))
  
  # Scatters
  channelpairs <- readxl::read_xlsx(
    file.path(base_dir,
              "Metadata",
              paste0("210420_", panel, "_lineage marker combinations for plotting.xlsx")),
    col_names = FALSE)
  channelpairs <- as.list(as.data.frame(t(channelpairs)))
  Plot2DScatters(fsom,
                 channelpairs = channelpairs,
                 metaclusters = c(1:40),
                 plotFile = paste0(date,"_", panel, "_Scatter_MC1-40.png"))
}

# Figures UMAP -----------------------------------------------------------------

# _ UMAP model -----------------------------------------------------------------

set.seed(1)
umap <- uwot::umap(scale(All_traits %>% dplyr::select(contains("%"))),
                   n_neighbors = 10)

All_traits <- data.frame(All_traits,
                         umap_1 = umap[,1],
                         umap_2 = umap[,2],
                         check.names = FALSE)
All_traits$Cohort[All_traits$Group == "QC"] = "QC"

ggplot(All_traits) +
  geom_point(aes(x = umap_1, y = umap_2, col = Age, shape = Cohort)) +
  theme_minimal() +
  scale_color_distiller(palette = "RdYlBu")


# _ Main overview colored by age -----------------------------------------------

umap_age <- ggplot(All_traits %>% dplyr::filter(Group != "QC")) +
  geom_point(aes(x = umap_1, y = umap_2, col = Age)) +
  geom_point(aes(x = umap_1, y = umap_2), col = "black",
             data = All_traits %>% dplyr::filter(Group == "QC")) +
  theme_void() +
  scale_color_distiller(palette = "RdYlBu") +
  theme(legend.position="bottom")


# _ Subset highlights ----------------------------------------------------------

tnk_mc_mfi <- GetMetaclusterMFIs(fsom_TNK, prettyColnames = TRUE, colsUsed = TRUE)
rownames(tnk_mc_mfi) <- paste0("TNK MC",seq_len(nrow(tnk_mc_mfi)))
tnk_mc_mfi <- tnk_mc_mfi[,stringr::str_sort(colnames(tnk_mc_mfi), numeric = TRUE)]
colnames(tnk_mc_mfi)[2] <- "CD1d" # Too long for plot !! make note in figure legend

bdc_mc_mfi <- GetMetaclusterMFIs(fsom_BDC, prettyColnames = TRUE, colsUsed = TRUE)
rownames(bdc_mc_mfi) <- paste0("BDC MC",seq_len(nrow(bdc_mc_mfi)))
bdc_mc_mfi <- bdc_mc_mfi[,stringr::str_sort(colnames(bdc_mc_mfi), numeric = TRUE)]


source("210708_compareSubset.R")


subset_A <- All_traits$umap_1 < -4
compare_A <- compareSubset(All_traits, subset_A, names = c("Main population", "Subpopulation"))
writexl::write_xlsx(compare_A$comparison, 
                    paste0(date,"_compare_A.xlsx"))
All_traits$subimmunotype <-subset_A

compare_A$plot
ggsave("compare_subimmunotype.pdf")
writexl::write_xlsx(All_traits,
                    paste0(date, "_All_traits.xlsx"))  

correlation_all <- compareSubsetAge(All_traits, rep(TRUE, nrow(All_traits)))
writexl::write_xlsx(correlation_all$correlation, 
                    paste0(date,"_correlation_all.xlsx"))

# Main figure ------------------------------------------------------------------
schematic <- png::readPNG("schematicOverview_simplified.png", native = TRUE)

wrap_elements(schematic) / 
(umap_age + correlation_all$plot  + plot_layout(widths = c(1, 2))) /
  (compare_A$plot) +  #  + plot_layout(tag_level = 'new')
  patchwork::plot_annotation(tag_levels = c("A"))
ggsave(paste0(date, "_main_overview.pdf"),
       width = 20, height = 15)

umap_age 

# Figure patient heatmap -------------------------------------------------------

dev.off()
pdf(paste0(date, "_patient_heatmap.pdf"), width = 20, height = 60)
pheatmap::pheatmap(All_traits[,grep("%", colnames(All_traits))],
                   scale = "column",
                   breaks = seq(-5,5, length.out = 100),
                   annotation_row = data.frame("QC" = (All_traits[,"Cohort"] == "QC") + 1,
                                               "Subset A" = as.numeric(subset_A),
                                               # "Subset C" = as.numeric(subset_C),
                                               # "Subset B" = as.numeric(subset_B),
                                               (All_traits[,"CMV", drop = FALSE] > 0) + 1,
                                               All_traits[,rev(c("Run#", "Age", "Gender"))],
                                               check.names = FALSE))
dev.off()

# Figure supplementary demographics --------------------------------------------

umap_gender <- ggplot(All_traits) +
  geom_point(aes(x = umap_1, y = umap_2, col = Gender)) +
  theme_void()

umap_CMV <- ggplot(All_traits[c(which(is.na(All_traits$CMV)),which(!is.na(All_traits$CMV))),]) +
  geom_point(aes(x = umap_1, y = umap_2, col = CMV > 0)) +
  theme_void() +
  scale_color_manual(values = c("TRUE" = "#d95f02", "FALSE" = "#7570b3"), 
                     na.value = "lightgrey")

umap_Run <- ggplot(All_traits) +
  geom_point(aes(x = umap_1, y = umap_2, col = `Run#`)) +
  geom_text(aes(x = umap_1, y = umap_2, label = `Run#`), col = "red", 
             data = dplyr::filter(All_traits, Cohort == "QC")) +
  theme_void()

(umap_age + umap_gender) / (umap_CMV + umap_Run) +
  plot_annotation(title = "Demographic information",
                  subtitle = paste0("Including ", nrow(All_traits), " samples,",
                                    " visualised with UMAP using 80 phenotypic features."))

ggsave(paste0(date,"_demographics.pdf"),
       width = 20, height = 15)

# Figures stability ------------------------------------------------------------

distances_longitudinal_HD <- All_traits[All_traits$Group == "Longitudinal",] %>% 
  group_by(`Sample ID`) %>% 
  summarise(Distance = mean(dist(cbind(across(contains("%"))))),
            Time = max(`Weeks from first timepoint`))

distances_longitudinal_ref <- All_traits[All_traits$Group == "Longitudinal",] %>% 
  mutate(`Sample ID` = sample(`Sample ID`)) %>% 
  group_by(`Sample ID`) %>% 
  summarise(Distance = mean(dist(cbind(across(contains("%"))))),
            Time = max(`Weeks from first timepoint`))

distances_twins <- All_traits[All_traits$Cohort == "UK",] %>% 
  dplyr::group_by(`New Family ID`, `Zygosity`) %>% 
  dplyr::summarise(Distance = mean(dist(cbind(across(contains("%"))))))

distances_twins_ref <- All_traits[All_traits$Cohort == "UK",] %>%  
  mutate(`New Family ID` = sample(`New Family ID`)) %>% 
  dplyr::group_by(`New Family ID`) %>% 
  dplyr::summarise(Distance = mean(dist(cbind(across(contains("%"))))))

pdf(paste0(date, "_Distances_longitudinal_HD.pdf"))
plot(density(na.omit(distances_longitudinal_ref$Distance)), 
     ylim = c(0,8),
     main = "Distances HD",
     xlab = "Distance")
lines(density(na.omit(distances_longitudinal_HD$Distance)), col = "red")
lines(density(na.omit(distances_twins$Distance)), col = "cyan")
lines(density(na.omit(distances_twins %>% 
                        dplyr::filter(Zygosity == "MZ") %>% 
                        dplyr::pull(Distance))), col = "blue")
lines(density(na.omit(distances_twins %>% 
                        dplyr::filter(Zygosity == "DZ") %>% 
                        dplyr::pull(Distance))), col = "green")
lines(density(na.omit(distances_twins_ref$Distance)), col = "grey")
legend("topright",
       legend = c("longitudinal","longitudinal_scrambled",
                  "twin", "twin_MZ", "twin_DZ", "twin_scrambled"),
       col = c("red", "black",
               "cyan", "blue", "green", "grey"),
       lty = 1)
dev.off()

distances_longitudinal_2D <- All_traits[All_traits$Group == "Longitudinal",] %>% 
  group_by(`Sample ID`) %>% 
  summarise(Distance = mean(dist(cbind(across(contains("umap"))))),
            Time = max(`Weeks from first timepoint`))

distances_longitudinal <- merge(distances_longitudinal_HD,
                                distances_longitudinal_2D,
                                by = "Sample ID",
                                suffixes = c("_HD", "_2D")) 

cor_dist_long <- cor(distances_longitudinal$Distance_HD, 
    distances_longitudinal$Distance_2D,
    use = "complete.obs")

ggplot(distances_longitudinal,
       aes(x = Distance_HD, y = Distance_2D)) +
  geom_point() + theme_minimal() +
  ggtitle(paste0("Correlation HD vs 2D: ", round(cor_dist_long, 2)))
ggsave(paste0(date,"_Distances_longitudinal_HDvs2D.pdf"))



distances_longitudinal <- All_traits[All_traits$Group == "Longitudinal",] %>% 
  group_by(`Sample ID`) %>% 
  summarise(Distance = mean(dist(cbind(.data[["umap_1"]], 
                                       .data[["umap_2"]]))),
            Time = max(.data[["Weeks from first timepoint"]]))

threshold <- 1

density_longitudinal <- ggplot(distances_longitudinal) + 
  geom_histogram(aes(Distance), binwidth = 0.2) +
  geom_vline(aes(xintercept = threshold)) +
  theme_minimal()
distance_count <- table(distances_longitudinal$Distance <= threshold)

to_plot <- merge(All_traits, distances_longitudinal,
                 by = "Sample ID", all.x = TRUE)
writexl::write_xlsx(to_plot,
                    "distances_longitudinal.xlsx")

umap_longitudinal <- 
  ggplot(to_plot[c(which(is.na(to_plot$Distance)),which(!is.na(to_plot$Distance))),] ) +
  geom_point(aes(x = umap_1, y = umap_2, col = Distance)) + # as.factor(`Run#`))) +
  geom_line(aes(x = umap_1, y = umap_2, group = as.factor(`Sample ID`), col = Distance),
            data = to_plot[!is.na(to_plot$Distance), ]) +
  scale_color_distiller(palette = "RdYlBu", na.value = "lightgrey") +
  theme_minimal() +
  ggtitle(paste0(distance_count["FALSE"], " distances larger than ", threshold))

umap_longitudinal_smallDistanceOnly <-
  ggplot(to_plot[c(which(is.na(to_plot$Distance)),which(!is.na(to_plot$Distance))),] ) +
  geom_point(aes(x = umap_1, y = umap_2, col = Distance)) +
  #geom_point(aes(x = umap_1, y = umap_2, col = as.factor(`Run#`))) +
  geom_line(aes(x = umap_1, y = umap_2, group = as.factor(`Sample ID`), col = Distance),
            data = to_plot[!is.na(to_plot$Distance) & to_plot$Distance <= threshold, ]) +
  scale_color_distiller(palette = "RdYlBu", na.value = "lightgrey") +
  theme_minimal() +
  ggtitle(paste0(distance_count["TRUE"], " distances smaller or equal to ", threshold))


density_longitudinal / (umap_longitudinal + umap_longitudinal_smallDistanceOnly )

ggsave(paste0(date, "_longitudinalChanges.pdf"),
       width = 20, height = 11.6)



distances_twins <- All_traits[All_traits$Cohort == "UK",] %>% 
  dplyr::group_by(`New Family ID`, `Zygosity`) %>% 
  dplyr::summarise(Distance = mean(dist(cbind(.data[["umap_1"]], .data[["umap_2"]]))))

density_twins <- ggplot(distances_twins) + 
  geom_histogram(aes(Distance, fill = Zygosity), binwidth = 0.1) +
  geom_vline(aes(xintercept = threshold)) +
  theme_minimal()
distance_twins_count <- table(distances_twins$Distance <= threshold)

umap_twins <- ggplot(All_traits ) +
  geom_point(aes(x = umap_1, y = umap_2, col = as.factor(`Zygosity`))) +
  geom_line(aes(x = umap_1, y = umap_2, col = as.factor(`Zygosity`), 
                group = as.factor(`New Family ID`)),
            data = All_traits[All_traits$Cohort == "UK",]) +
  theme_minimal()


density_twins / umap_twins

ggsave(paste0(date, "_twinChanges.pdf"),
       width = 20, height = 11.6)




# Figures CV -------------------------------------------------------------------

# _ Frequency ----

All_traits_stacked <- tidyr::pivot_longer(All_traits, 
                                          cols = grep("%", colnames(All_traits), value = TRUE),
                                          names_to = "Trait", 
                                          values_to = "Frequency")
All_traits_stacked$Panel <- stringr::str_sub(All_traits_stacked$Trait, 1, 3)

Frequency_plot <-    ggplot(All_traits_stacked, 
                            aes(x = reorder(Trait, Frequency, FUN = mean), 
                                y = log10(Frequency*100))) +
  geom_jitter(width = 0.2, size = 0.1, alpha = 0.3, color = "gray60") +
  geom_boxplot(fatten = NULL, outlier.shape = NA, coef = 0, fill = "white") + # fatten = NULL to make midline disappear, causes warnings which can be safely ignored
  # geom_boxplot(All_traits_stacked %>% dplyr::filter(Cohort == "QC"),
  #              mapping = aes(x = Trait,
  #                            y = log10(Frequency*100),
  #                            fill = Panel, color = Panel),
  #              width = 0.4,
  #              fatten = 0, outlier.shape = NA, coef = 0) +
  geom_point(All_traits_stacked %>% dplyr::filter(Cohort == "QC"),
               mapping = aes(x = Trait,
                             y = log10(Frequency*100),
                             fill = Panel, color = Panel),
             position = ggbeeswarm::position_quasirandom()) + 
  annotation_logticks(sides = c("l"), scaled = TRUE) +
  theme_classic() +
  xlab("Immune trait") + ylab("Frequency (log10 %)") +
  theme(axis.text.x = element_text(size = 7, angle = 90, vjust = 0.5, hjust=1),
        legend.position = "none") +
  scale_color_manual(values = c("BDC" = "#F8766D", 
                                "ICS" = "#00BA38",
                                "TNK" = "#619CFF")) +
  scale_fill_manual(values = c("BDC" = "#F8766D", 
                               "ICS" = "#00BA38",
                               "TNK" = "#619CFF"))
Frequency_plot

# _ QC CV vs count ----

CV_QC <- All_traits_stacked %>% 
  dplyr::filter(Cohort == "QC") %>% 
  group_by(Trait) %>% 
  dplyr::summarize(Panel = unique(Panel),
                   SD = sd(Frequency, na.rm=TRUE),
                   Mean = mean(Frequency, na.rm=TRUE)) %>% 
  mutate(CV = SD/Mean * 100)

CV_QC$Count <- sapply(CV_QC$Trait, 
                      function(trait){
                        mean(counts[counts$Group == "QC", 
                                    gsub("%", "", trait)])})

Poison_fun <- function(x) ((x^(1/2))/x)*100

CV_freq_vs_count <- ggplot(CV_QC[sample(nrow(CV_QC)), ], 
                           aes(y = `CV`, x = `Count`, color = `Panel`)) +
  stat_function(fun = Poison_fun, geom = "area", color = "gray", fill = "gray")  +
  geom_point(size = 2, alpha = 0.5) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)),
                expand = c(0, 0),
                #limits = c(0.05, 200000)) +
                limits = c(1, 200000)) +
  #scale_y_continuous(limits = c(0,450), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,150), expand = c(0, 0)) +
  annotation_logticks(sides = c("b"), scaled = TRUE) +
  coord_cartesian(clip = 'off') +
  geom_hline(yintercept = 30, color = "darkred", linetype="dashed") +
  ylab("%CV Frequency QC samples") + 
  xlab("Mean cell count QC samples (log)")  +
  theme_classic() +
  theme(legend.position = "none")+
  scale_color_manual(values = c("BDC" = "#F8766D", 
                                "ICS" = "#00BA38",
                                "TNK" = "#619CFF")) +
  scale_fill_manual(values = c("BDC" = "#F8766D", 
                               "ICS" = "#00BA38",
                               "TNK" = "#619CFF"))

CV_freq_vs_count

# _ CV all ----

All_traits_nomultiples_stacked <- All_traits_stacked %>% #Extract only weeks 0 and NA
  dplyr::filter(`Weeks from first timepoint` == 0 |
                  is.na(`Weeks from first timepoint`) &
                  duplicated == 0)

#Calculate Mean, SD and %CV for Frequency traits of all sample
CV_nomultiples_total <- All_traits_nomultiples_stacked %>% 
  group_by(Trait) %>%
  dplyr::summarize(Panel = unique(Panel),
                   SD_total = sd(Frequency, na.rm=TRUE),
                   Mean_total = mean(Frequency, na.rm=TRUE)) %>% 
  mutate(CV_total = SD_total/Mean_total * 100)

#Calculate CV values for just twinsUK  
CV_nomultiples_TwinsUK  <- All_traits_nomultiples_stacked  %>% 
  dplyr::filter(Cohort == "UK") %>%   
  group_by(Trait) %>%
  dplyr::summarize(Panel = unique(Panel),
                   SD_Twins = sd(Frequency, na.rm=TRUE),
                   Mean_Twins = mean(Frequency, na.rm=TRUE)) %>% 
  mutate(CV_Twins = SD_Twins/Mean_Twins * 100)

#Calculate CV values for just VRC  
CV_nomultiples_VRC  <- All_traits_nomultiples_stacked  %>% 
  dplyr::filter(Cohort == "VRC") %>%   
  group_by(Trait) %>%
  dplyr::summarize(Panel = unique(Panel),
                   SD_VRC = sd(Frequency, na.rm=TRUE),
                   Mean_VRC = mean(Frequency, na.rm=TRUE)) %>% 
  mutate(CV_VRC = SD_VRC/Mean_VRC * 100)


colnames(CV_QC)[3:7] <- paste0(colnames(CV_QC)[3:7], "_QC")

CV_nomultiples <- CV_nomultiples_total %>% 
  inner_join(CV_nomultiples_TwinsUK, by = c("Trait", "Panel")) %>% 
  inner_join(CV_nomultiples_VRC, by = c("Trait", "Panel")) %>% 
  inner_join(CV_QC, by = c("Trait", "Panel"))


PlotCVComparison <- function(CV_nomultiples,
                             var1,
                             var2,
                             lab1,
                             lab2,
                             title){
  ggplot(CV_nomultiples, 
         aes_string(x = var1, y = var2)) +
    geom_rect(aes(xmin=0, xmax=30, ymin=0, ymax=Inf), fill = "gray80") +
    geom_rect(aes(xmin=0, xmax=Inf, ymin=0, ymax=30), fill = "gray80") +
    geom_point(aes(color = Panel)) +
    scale_x_continuous(limits = c(0, 200), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0,200), expand = c(0, 0)) +
    scale_color_manual(values = c("BDC" = "#F8766D", 
                                  "ICS" = "#00BA38",
                                  "TNK" = "#619CFF")) +
    xlab(lab1) +
    ylab(lab2) +
    ggtitle(title) +
    geom_abline(intercept = 0, slope = 1, color = "darkred", linetype = "dashed") +
    theme_classic()
}

#Plot CV total vs QC
CV1 <- PlotCVComparison(CV_nomultiples,
                        "CV_QC", "CV_total",
                        "%CV - QC Sample", "%CV - All samples",
                        "QV vs all samples")

CV2 <- PlotCVComparison(CV_nomultiples,
                        "CV_Twins", "CV_VRC",
                        "%CV - Twins UK cohort", "%CV - VRC cohort",
                        "TwinsUK vs VRC")

# _ Overview figure ----

Frequency_plot / (CV_freq_vs_count + CV1 + CV2)
ggsave(paste0(date, "_CV_overview.pdf"),
       width = 20, height = 11.6)


# Figures batch effect ---------------------------------------------------------

set.seed(1)
umap_TNK <- uwot::umap(scale(All_traits %>% dplyr::select(contains("TNK %"))),
                       n_neighbors = 10)
umap_BDC_orig <- uwot::umap(scale(MC_pctgs_BDC[complete, ]),
                            n_neighbors = 10)
umap_BDC_norm <- uwot::umap(scale(All_traits %>% dplyr::select(contains("BDC %"))),
                            n_neighbors = 10)


All_traits[, c("umap_tnk_1", "umap_tnk_2",
               "umap_bdc_o_1", "umap_bdc_o_2",
               "umap_bdc_n_1", "umap_bdc_n_2")] <- cbind(umap_TNK,
                                                         umap_BDC_orig,
                                                         umap_BDC_norm)

set.seed(1)
p_umap_tnk <- ggplot(All_traits[sample(nrow(All_traits)),]) +
  geom_point(aes(x = umap_tnk_1, y = umap_tnk_2, col = `Run#`)) +
  theme_void() +
  scale_color_distiller(palette = "RdYlBu") +
  ggtitle("UMAP based on 40 TNK traits")

p_umap_bdc <- ggplot(All_traits[sample(nrow(All_traits)),]) +
  geom_point(aes(x = umap_bdc_o_1, y = umap_bdc_o_2, col = `Run#`)) +
  theme_void() +
  scale_color_distiller(palette = "RdYlBu") +
  ggtitle("UMAP based on 40 BDC traits - original")

p_umap_bdc_n <- ggplot(All_traits[sample(nrow(All_traits)),]) +
  geom_point(aes(x = umap_bdc_n_1, y = umap_bdc_n_2, col = `Run#`)) +
  theme_void() +
  scale_color_distiller(palette = "RdYlBu") +
  ggtitle("UMAP based on 40 BDC traits - normalized")

(p_umap_tnk + p_umap_bdc + p_umap_bdc_n)
ggsave(paste0(date, "_umap_separate.pdf"),
       width = 20, height = 5)

# demographics %>% 
#   filter(Group == "QC") %>% 
#   arrange("Run #") %>%
#   pull(`FlowJo ID`) %>% 
#   paste0("preprocessed_BDC/", . ,".fcs") %>% 
#   FlowSOM::PlotFileScatters(channels = GetChannels(fsom_BDC, 
#                                                    c("CD123", "CD11c")),
#                             names = paste0("Run ", 1:19),
#                             plotFile = "QC_CD123_CD11c.png",
#                             
#                             ncol = 2,
#                             nrow = 1)

norm_model <- readRDS("RDS/CytoNorm_model_BDC.RDS")

#plot(norm_model$clusterRes$`1`$splines$`1`$`B610-A`

plotSpline <- function(norm_model, cluster, batch, channel,
                       min = 0, max = 8){
     graphics::plot(norm_model$clusterRes[[cluster]]$quantiles[[batch]][,channel], 
                    norm_model$clusterRes[[cluster]]$refQuantiles[,channel],
                    xlim = c(min, max), ylim = c(min, max), pch = 19, 
                    bty = "n", xaxt = "n", yaxt = "n", 
                    xlab = "", ylab = "", main = FlowSOM::GetMarkers(norm_model$fsom, channel))
     graphics::lines(c(min, max), c(min, max), col = "#999999")
     x <- seq(min, max, 0.1)
     graphics::lines(x, norm_model$clusterRes[[cluster]]$splines[[batch]][[channel]](x), 
                     col = "#b30000")
}


png(paste0(date, "_splines_overview.png"),
    width = 10 * 200, height = 17 * 200)
graphics::layout(matrix(seq_len(10*17), ncol = 10))
graphics::par(mar = c(1, 1, 5, 0))
for(cluster in names(norm_model$clusterRes)){
  message(Sys.time(),": Cluster ",cluster)
  for(channel in names(norm_model$fsom$prettyColnames[1:17])){
    plotSpline(norm_model, 
               cluster = cluster, 
               batch = 2,
               channel = channel,
               min = -1, max = 5)
  }
}
dev.off()

# MFI --------------------------------------------------------------------------

MFI_traits <- data.frame(demographics[complete,],
                         MC_MFI_BDC_norm[complete,],
                         MC_MFI_TNK[complete,],
                         check.names = FALSE)

MFI_traits_stacked <- tidyr::pivot_longer(MFI_traits, 
                                          cols = grep("MC", colnames(MFI_traits), value = TRUE),
                                          names_to = "Trait", 
                                          values_to = "MFI")
MFI_traits_stacked$Panel <- stringr::str_sub(MFI_traits_stacked$Trait, 1, 3)

MFI_plot <- 
  MFI_traits_stacked %>% 
  mutate(Trait = forcats::fct_reorder(Trait, MFI, median, na.rm = TRUE)) %>%
  ggplot(aes(x = Trait, y = MFI)) + 
  #geom_jitter(width = 0.2, size = 0.1, alpha = 0.3, color = "gray60") +
  geom_boxplot(fatten = NULL, outlier.shape = NA, coef = 0, fill = "white") + # fatten = NULL to make midline disappear, causes warnings which can be safely ignored
  geom_point(MFI_traits_stacked %>% dplyr::filter(Group == "QC"),
             mapping = aes(x = Trait,
                           y = MFI, 
                           fill = Panel, color = Panel),
             position = ggbeeswarm::position_quasirandom(),
             size = 0.5) + 
  theme_classic() +
  xlab("Immune trait") + ylab("MFI (biexponentially transformed)") + 
  theme(axis.text.x = element_text(size = 7, angle = 90, vjust = 0.5, hjust=1),
        legend.position = "none") +
  scale_color_manual(values = c("BDC" = "#F8766D", 
                                "ICS" = "#00BA38",
                                "TNK" = "#619CFF")) +
  scale_fill_manual(values = c("BDC" = "#F8766D", 
                               "ICS" = "#00BA38",
                               "TNK" = "#619CFF"))
MFI_plot
ggsave(paste0(date, "_MFI_traits.pdf"),
       width = 65, height = 8, limitsize = FALSE)

# _ CV ----

CV_MFI_QC <- MFI_traits_stacked %>% 
  dplyr::filter(Group == "QC") %>% 
  group_by(Trait) %>% 
  dplyr::summarize(Panel = unique(Panel),
                   SD = sd(MFI, na.rm=TRUE),
                   Mean = mean(MFI, na.rm=TRUE)) %>% 
  mutate(CV = SD/Mean * 100)

CV_MFI_QC$Count <- sapply(CV_MFI_QC$Trait, 
                          function(trait){
                            mean(counts[counts$Group == "QC", 
                                        gsub("(.*) (.*) (.*)", "\\1 \\2", trait)])})
CV_MFI_freq_vs_count <- ggplot(CV_MFI_QC[sample(nrow(CV_MFI_QC)), ], 
                               aes(y = `CV`, x = `Count`, color = `Panel`)) +
  stat_function(fun = Poison_fun, geom = "area", color = "gray", fill = "gray")  +
  geom_point(size = 2, alpha = 0.5) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)),
                expand = c(0, 0),
                #limits = c(0.05, 200000)) +
                limits = c(1, 200000)) +
  #scale_y_continuous(limits = c(0,450), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,150), expand = c(0, 0)) +
  #annotation_logticks(sides = c("b"), scaled = TRUE) +
  coord_cartesian(clip = 'off') +
  geom_hline(yintercept = 30, color = "darkred", linetype="dashed") +
  ylab("%CV MFI QC samples") + 
  xlab("Mean cell count QC samples (log)")  +
  theme_classic() +
  theme(legend.position = "none")+
  scale_color_manual(values = c("BDC" = "#F8766D", 
                                "ICS" = "#00BA38",
                                "TNK" = "#619CFF")) +
  scale_fill_manual(values = c("BDC" = "#F8766D", 
                               "ICS" = "#00BA38",
                               "TNK" = "#619CFF"))

CV_MFI_freq_vs_count


MFI_traits_nomultiples_stacked <- MFI_traits_stacked %>% #Extract only weeks 0 and NA
  dplyr::filter(`Weeks from first timepoint` == 0 |
                  is.na(`Weeks from first timepoint`) &
                  duplicated == 0)

#Calculate Mean, SD and %CV for MFI traits of all sample
CV_MFI_nomultiples_total <- MFI_traits_nomultiples_stacked %>% 
  group_by(Trait) %>%
  dplyr::summarize(Panel = unique(Panel),
                   SD_total = sd(MFI, na.rm=TRUE),
                   Mean_total = mean(MFI, na.rm=TRUE)) %>% 
  mutate(CV_total = SD_total/Mean_total * 100)

#Calculate CV values for just twinsUK  
CV_MFI_nomultiples_TwinsUK  <- MFI_traits_nomultiples_stacked  %>% 
  dplyr::filter(Cohort == "UK") %>%   
  group_by(Trait) %>%
  dplyr::summarize(Panel = unique(Panel),
                   SD_Twins = sd(MFI, na.rm=TRUE),
                   Mean_Twins = mean(MFI, na.rm=TRUE)) %>% 
  mutate(CV_Twins = SD_Twins/Mean_Twins * 100)

#Calculate CV values for just VRC  
CV_MFI_nomultiples_VRC  <- MFI_traits_nomultiples_stacked  %>% 
  dplyr::filter(Cohort == "VRC") %>%   
  group_by(Trait) %>%
  dplyr::summarize(Panel = unique(Panel),
                   SD_VRC = sd(MFI, na.rm=TRUE),
                   Mean_VRC = mean(MFI, na.rm=TRUE)) %>% 
  mutate(CV_VRC = SD_VRC/Mean_VRC * 100)


colnames(CV_MFI_QC)[3:7] <- paste0(colnames(CV_MFI_QC)[3:7], "_QC")

CV_MFI_nomultiples <- CV_MFI_nomultiples_total %>% 
  inner_join(CV_MFI_nomultiples_TwinsUK, by = c("Trait", "Panel")) %>% 
  inner_join(CV_MFI_nomultiples_VRC, by = c("Trait", "Panel")) %>% 
  inner_join(CV_MFI_QC, by = c("Trait", "Panel"))

CV_MFI_1 <- PlotCVComparison(CV_MFI_nomultiples,
                        "CV_QC", "CV_total",
                        "%CV - QC Sample", "%CV - All samples",
                        "QV vs all samples")

CV_MFI_2 <- PlotCVComparison(CV_MFI_nomultiples,
                        "CV_Twins", "CV_VRC",
                        "%CV - Twins UK cohort", "%CV - VRC cohort",
                        "TwinsUK vs VRC")

CV_MFI_freq_vs_count + CV_MFI_1 + CV_MFI_2
ggsave(paste0(date, "_MFI_CV.pdf"),
       width = 20, height = 6)

# Outliers ----

outliers <- data.frame(demographics[complete,],
                       outliers_TNK = MC_outliers_TNK[complete],
                       outliers_BDC_original = MC_outliers_BDC[complete],
                       outliers_BDC = MC_outliers_BDC_norm[complete],
                       check.names = FALSE)
p_outliers_TNK <- ggplot(outliers) +
  geom_point(aes(x = reorder(`FlowJo ID`, `Run#`),
                 y = outliers_TNK * 100,
                 col = `Run#`)) +
  scale_color_distiller(palette = "RdYlBu") +
  ylim(0,30)

p_outliers_BDC_o <- ggplot(outliers) +
  geom_point(aes(x = reorder(`FlowJo ID`, `Run#`),
                 y = outliers_BDC_original  * 100,
                 col = `Run#`)) +
  scale_color_distiller(palette = "RdYlBu") +
  ylim(0,30)

p_outliers_BDC <- ggplot(outliers) +
  geom_point(aes(x = reorder(`FlowJo ID`, `Run#`),
                 y = outliers_BDC  * 100,
                 col = `Run#`)) +
  scale_color_distiller(palette = "RdYlBu") +
  ylim(0,30)

p_outliers_norm <- ggplot(outliers) +
  geom_point(aes(x = outliers_BDC_original * 100,
                 y = outliers_BDC  * 100,
                 col = `Run#`)) +
  scale_color_distiller(palette = "RdYlBu") +
  theme_minimal() +
  geom_abline(aes(slope = 1, intercept = 0))

(p_outliers_TNK / p_outliers_BDC_o / p_outliers_BDC) | p_outliers_norm
  
ggsave(paste0(date, "_outliers.pdf"),
       width = 20, height = 11.6)
  
