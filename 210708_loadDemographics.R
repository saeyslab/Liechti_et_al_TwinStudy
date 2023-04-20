library(dplyr)

loadDemographics <- function(base_dir){
  demographics <- read.csv(paste0(base_dir, "/Metadata/200901_Twin_demographics_adjusted.csv"),
                           check.names = FALSE,
                           stringsAsFactors = FALSE)
  demographics$`FlowJo ID` <- as.character(demographics$`FlowJo ID`)
  demographics$`Sample ID` <- NULL # Remove, does not contain full information
  
  # Add run
  run <-  read.csv(paste0(base_dir, "/Metadata/201122_FlowJo ID and Run.csv"),
                   check.names = FALSE)
  rownames(run) <- run$`FlowJo ID` <- as.character(run$`FlowJo ID`)
  demographics <- merge(demographics, run, 
                        by = "FlowJo ID",
                        all = TRUE, sort = FALSE)
  
  # Add sample ID
  sampleID <-  read.csv(paste0(base_dir, "/Metadata/201122_FlowJo and Sample ID.csv"),
                        check.names = FALSE)
  rownames(sampleID) <- sampleID$`FlowJo ID` <- as.character(sampleID$`FlowJo ID`)
  demographics <- merge(demographics, sampleID, 
                        by = c("FlowJo ID"),
                        all = TRUE, sort = FALSE)
  
  # Add CMV status
  CMV_VRC <- read.csv(file.path(base_dir, "Metadata/2021-06-24 VRC serostatus for EBV CMV and HSV1_2.csv"),
                      check.names = FALSE,
                      stringsAsFactors = FALSE)
  CMV_VRC <- CMV_VRC %>% dplyr::filter(`Sample ID` != "S3-0860-01")
  CMV_twins <- read.csv(file.path(base_dir, "Metadata/2021-06-24 Twins serostatus for EBV CMV and HSV1_2.csv"),
                        check.names = FALSE,
                        stringsAsFactors = FALSE)
  CMV <- rbind(CMV_VRC, CMV_twins)
  CMV[is.na(CMV)] <- 0
  demographics <- merge(demographics, CMV, 
                        by = c("Sample ID"),
                        all.x = TRUE, sort = FALSE)
  
  # Add viability
  viability <- sapply(c("BDC", "TNK", "ICS"),
                        function(sheet) readxl::read_xlsx(file.path(base_dir, 
                                                                    "Metadata/210808_Viability summary_updated.xlsx"),
                                                          sheet = sheet) %>% as.data.frame())
  colnames(viability$ICS)[3:5] <- paste0("ICS_", colnames(viability$ICS)[3:5])
  viability_all <- merge(viability$BDC[,2:3], 
                         viability$TNK[,2:3],
                         by = "FlowJo ID", 
                         suffixes = c("_BDC", "_TNK"),
                         all = TRUE)
  viability_all <- merge(viability_all, 
                         viability$ICS[2:5], 
                         by = "FlowJo ID",
                         all = TRUE)
  colnames(viability_all)
  demographics <- merge(demographics, viability_all, 
                        by = c("FlowJo ID"),
                        all.x = TRUE, sort = FALSE)
  
  rownames(demographics) <- demographics$`FlowJo ID`
  demographics
}
