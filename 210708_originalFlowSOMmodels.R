devtools::install_github("saeyslab/FlowSOM")
devtools::install_github("saeyslab/CytoNorm")

# Load libraries and helper functions ------------------------------------------

library(flowCore)
library(FlowSOM)
source("twin_study/210708_loadDemographics.R")

# Set up variables -------------------------------------------------------------

recompute <- FALSE
date <- "210708"
base_dir <- "/auto/net/fs4-Flowcytometry/Storage/u_ysa/Projects/twin_study/52e4c244-0e9b-4573-8e11-d10c249a871c-twinstudy2/"

demographics <- loadDemographics(base_dir)
demographics <- dplyr::filter(demographics, `ICS_Viable of CD3` > 25)
training_set <- demographics[demographics[,"Run#"] %in% c(1,2), "FlowJo ID"]

# Panel specific variables -----------------------------------------------------

panel <- "BDC"
#panel <- "TNK"

# Raw files
fcs_path <- file.path(base_dir, panel, "fcs files_comped and pregated")

# New directories and files which will be created by the script
preprocessed_path <- file.path(base_dir, paste0(date, "_", panel, "_preprocessed"))
agg_path <- file.path(base_dir, "RDS", paste0(date, "_", panel, "_aggregate.fcs"))
fsom_path <- file.path(base_dir, "RDS", paste0(date, "_", panel, "_FlowSOM.RDS"))
output_path <- file.path(base_dir, "RDS", paste0(date, "_", panel, "_Results.RDS"))

# Markers to use
marker_files <- c("TNK" = "210203_TNK_Marker and Population for FlowSOM.xlsx",
                  "BDC" = "210420_BDC_Marker and subset list.xlsx")

marker_overview <- readxl::read_xlsx(file.path(base_dir, "Metadata", 
                                               marker_files[panel]),
                                     sheet = 1)

channels_to_use <- marker_overview$Channel[marker_overview[[4]] %in% c("Lineage", "Subset identification")]
channels_to_use <- paste0(channels_to_use, "-A")
channels_to_use <- gsub("U395-A", "U390-A", channels_to_use)
channels_to_use <- gsub("U670-A", "U660-A", channels_to_use)

channels_MFI <- marker_overview$Channel[marker_overview[[4]] %in% c("Functional")]
channels_MFI <- paste0(channels_MFI, "-A")
channels_MFI <- gsub("U395-A", "U390-A", channels_MFI)
channels_MFI <- gsub("U670-A", "U660-A", channels_MFI)

# Preprocessing preparation ----------------------------------------------------

# Manually defined BDC gate
BDC_gate <- flowCore::polygonGate(filterId = "CD14 or HLADR positive",
                                  .gate = matrix(c(-1, 2, 
                                                   1.8, 2,
                                                   1.8, -1,
                                                   4, -1,
                                                   4, 4,
                                                   -1, 4),
                                                 ncol = 2, byrow = TRUE,
                                                 dimnames = list(NULL, c("U660-A", "V510-A"))))

# Load in first file to estimate transformation 

ff <- read.FCS(file.path(fcs_path, "1.fcs"))
ff <- compensate(ff, keyword(ff)[["$SPILLOVER"]])
tf <- estimateLogicle(ff,
                      colnames(ff)[4:31])

# Preprocess the data ----------------------------------------------------------
if(!dir.exists(preprocessed_path)) dir.create(preprocessed_path)

set.seed(1)
for(file in paste0(demographics$`FlowJo ID`, ".fcs")){
  if(file.exists(file.path(fcs_path, file)) & 
     !file.exists(file.path(preprocessed_path, file))){
    message(Sys.time(),": Preprocessing ", file)
    ff <- read.FCS(file.path(fcs_path, file))
    ff <- compensate(ff, keyword(ff)[["$SPILLOVER"]])
    ff <- transform(ff, tf)
    if(panel == "BDC") ff <- ff[flowCore::filter(ff, BDC_gate)@subSet,]
    write.FCS(ff,
              file.path(preprocessed_path, file))
  }
}

# Make an aggregate of all training samples ------------------------------------

if(!recompute & file.exists(agg_path)){
  agg <- read.FCS(agg_path)
  message(Sys.time(),": Reloaded aggregate file from ", agg_path)
} else {
  message(Sys.time(), ": Aggregating ", panel)
  set.seed(1)
  agg <- AggregateFlowFrames(file.path(preprocessed_path, 
                                       paste0(training_set, ".fcs")),
                             cTotal = 50000*length(training_set),
                             channels = colnames(ff)[1:32])
  write.FCS(agg, agg_path)
  message(Sys.time(), ": Aggregating finished")
}

# Train the FlowSOM model ------------------------------------------------------

if(!recompute & file.exists(fsom_path)){
  fsom <- readRDS(fsom_path)
  message("Reloaded previously saved FlowSOM model from ",
          fsom_path,
          "\n Set recompute = TRUE to recompute.")
} else {
  
  n_meta <- 40
  t_fsom <- system.time(
    fsom <- FlowSOM(agg,
                    colsToUse = channels_to_use,
                    xdim = 12, ydim = 12,
                    scale = FALSE,
                    nClus = n_meta,
                    seed = 1)
  )
  
  # Remove fluorochrome and channel
  fsom$prettyColnames <- gsub(" [^ ]* <[^ ]*>$", "", fsom$prettyColnames)
  # Remove channel if no fluorochrome was specified
  fsom$prettyColnames <- gsub(" <[^ ]*>$", "", fsom$prettyColnames)
  fsom$time <- t_fsom
  
  saveRDS(fsom, fsom_path)
}

# Map all samples onto the model -----------------------------------------------

recompute <- TRUE

for(i_start in seq(1, nrow(demographics), by = 100)){
  i_end <- min(i_start+99, nrow(demographics))
  message(Sys.time(),": ", i_start,"-",i_end)
  
  output_file  <- gsub(".RDS$", paste0("_", i_start,"-", i_end, ".RDS"), output_path)
  if(!recompute & file.exists(output_file)){
    message("Skipping previously saved features.",
            "\n Set recompute = TRUE to recompute.")
  } else {
    files <- file.path(preprocessed_path, paste0(demographics$`FlowJo ID`[i_start:i_end], ".fcs"))
    exists <- file.exists(files)
    features <- GetFeatures(fsom,
                            files = files[exists],
                            level = c("clusters", "metaclusters"),
                            type = c("counts", "percentages", "MFIs"),
                            MFI = channels_MFI,
                            filenames = demographics$`FlowJo ID`[i_start:i_end][exists],
                            silent = FALSE)
    
    saveRDS(features, output_file)
  }
}

# Normalisation BDC panel ------------------------------------------------------

library(CytoNorm)

# Training
panel <- "BDC"

QC_samples <- demographics[demographics$Group == "QC",]
QC_files <- paste0(QC_samples$`FlowJo ID`,".fcs")

batch <- QC_samples$`Run#` %in% c(5:19)
batch <- as.character(1 + batch)
normalisation_path <- file.path(base_dir, paste0(date, "_", panel, "_normalized"))

nQ <- 50
if(file.exists(file.path(normalisation_path, "CytoNorm_model_BDC.RDS"))){
  norm_model <- readRDS(file.path(normalisation_path, "CytoNorm_model_BDC.RDS"))
} else {
  norm_model <- CytoNorm::CytoNorm.train(files = file.path(preprocessed_path, QC_files),
                                         labels = batch,
                                         channels = c(channels_to_use, channels_MFI),
                                         transformList = NULL,
                                         outputDir = normalisation_path,
                                         FlowSOM.params = list(nCells = 1000000,
                                                               xdim = 12, ydim = 12,
                                                               nClus = 10,
                                                               channels = channels_to_use),
                                         normMethod.train = QuantileNorm.train,
                                         normParams = list(quantileValues = ((1:nQ)/nQ)[-nQ],
                                                           goal = "1",
                                                           limit = c(-1,4)),
                                         seed = 1,
                                         clean = TRUE,
                                         plot = FALSE,
                                         verbose = TRUE)
  saveRDS(norm_model, 
          file.path(normalisation_path, "CytoNorm_model_BDC.RDS"))
}

# Normalising
for(i_start in seq(1, nrow(demographics), by = 100)){
  i_end <- min(i_start+99, nrow(demographics))
  message(Sys.time(),": ", i_start,"-",i_end)
  
  files <- file.path(preprocessed_path, paste0(demographics$`FlowJo ID`[i_start:i_end], ".fcs"))
  exists <- file.exists(files)
  run <- demographics$`Run#`[i_start:i_end][exists]
  
  CytoNorm::CytoNorm.normalize(model = norm_model,
                               files = files[exists],
                               labels = as.character(1+(run %in% c(5:19))),
                               transformList = NULL,
                               transformList.reverse = NULL,
                               outputDir = normalisation_path,
                               prefix = "",
                               normMethod.normalize = QuantileNorm.normalize,
                               clean = TRUE,
                               verbose = TRUE)
}

# Map all normalised samples onto the model 

for(i_start in seq(1, nrow(demographics), by = 100)){
  i_end <- min(i_start+99, nrow(demographics))
  message(Sys.time(),": ", i_start,"-",i_end)
  
  output_file  <- gsub(".RDS$", paste0("_normalized_", i_start,"-", i_end, ".RDS"), output_path)
  if(!recompute & file.exists(output_file)){
    message("Skipping previously saved features.",
            "\n Set recompute = TRUE to recompute.")
  } else {
    files <- file.path(normalisation_path, paste0(demographics$`FlowJo ID`[i_start:i_end], ".fcs"))
    exists <- file.exists(files)
    features <- GetFeatures(fsom,
                            files = files[exists],
                            level = c("clusters", "metaclusters"),
                            type = c("counts", "percentages", "MFIs"),
                            MFI = channels_MFI,
                            filenames = demographics$`FlowJo ID`[i_start:i_end][exists],
                            silent = FALSE)
    
    saveRDS(features, output_file)
  }
}
