# Install and load required packages
# install.packages("remotes")
# library(remotes)
# remotes::install_github("MRCIEU/TwoSampleMR", force = TRUE)
library(TwoSampleMR)
library(dplyr)
library(stringr)
# Set working directory
setwd("D:/OneDrive - Hong Kong Metropolitan University/5.GWAS_Mendel Randomization Analysis")

### Load data ###
#load all available data
repeat {
  tryCatch({ao <- available_outcomes()
    if (!is.null(ao)) {break}},
    error = function(e) {
    message("Error in available_outcomes: ", e$message)})}

traits <- ao %>%
  mutate(batch = str_extract(id, "^\\w+-\\w+")) %>%
  filter(batch %in% c("ebi-a", "finn-b", "ieu-a", "ieu-b", "ukb-a", "ukb-b", "ukb-d", "ukb-e"),
         year > 2019,
         sample_size>10000)

# Loop through each trait in the filtered list
t0 <- Sys.time()
recorded_ids <- data.frame(xid= character(),
                           xname= character(),
                           yid = character(),
                           yname= character(),
                           pval = numeric(),
                           stringsAsFactors = FALSE)
recorded_file <- "20231226hk_year(2019)_samplesize(10000).csv"
write.table(recorded_ids, recorded_file, sep = ",", row.names = FALSE, col.names = TRUE)
batch_size <- 100
for (i in 1:nrow(traits)){
  x<-traits[i, , drop = FALSE]
  exposure <- extract_instruments(outcomes =x$id)
  if (is.null(exposure)) {exposure <- extract_instruments(outcomes =x$id)}
  for (j in 1:nrow(traits)) {
    y<- traits[j, , drop = FALSE]
    if(!is.null(exposure) && i!=j && x$population==y$population){
      tryCatch({
        mrresult <- mr(dat=harmonise_data(exposure_dat = exposure, 
                                          outcome_dat =extract_outcome_data(snps = exposure$SNP, 
                                                                            outcomes = y$id)), 
                       method_list = c("mr_ivw"))
        if (mrresult$pval < 0.01) {
          recorded_ids <- rbind(recorded_ids, data.frame(xid=x$id,
                                                         xname=x$trait,
                                                         yid=y$id,
                                                         yname=y$trait,
                                                         pval = mrresult$pval,stringsAsFactors = FALSE))}},
        error = function(e) {
          message("Error occurred: ", e$message,":",y$id)})}
    if (j %% batch_size == 0 || j == nrow(traits)) {
      write.table(recorded_ids, recorded_file, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
      recorded_ids <- data.frame(xid= character(),
                                 xname= character(),
                                 yid = character(),
                                 yname= character(),
                                 pval = numeric(),stringsAsFactors = FALSE)}}}