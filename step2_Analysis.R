# if (!require("devtools")) { install.packages("devtools") } else {}
# devtools::install_github("rondolab/MR-PRESSO",force = TRUE)
# devtools::install_github("qingyuanzhao/mr.raps",force = TRUE)
#install.packages("MendelianRandomization")
install.packages("PhenoScannerV2")
library(PhenoScannerV2)
library(TwoSampleMR)
library(data.table)
library(dplyr)
library(MRPRESSO)
library(mr.raps)
#data load
exposure_dat<- extract_instruments(outcomes ="ieu-b-4879")
exposure_dat2<- extract_instruments(outcomes ="finn-b-D3_SARCOIDOSIS")
exposure_dat <- clump_data(exposure_dat, clump_kb = 10000, clump_r2 = 0.001)
exposure_dat$F<-(exposure_dat$beta.exposure/exposure_dat$se.exposure)^2
exposure_dat2$F<-(exposure_dat2$beta.exposure/exposure_dat2$se.exposure)^2
filtered_dat <- exposure_dat %>% 
  filter(F > 10)
filtered_dat2 <- exposure_dat2 %>% 
  filter(F > 10)
exposure_dat2 <- clump_data(exposure_dat2, clump_kb = 10000, clump_r2 = 0.001)
outcome <- extract_outcome_data(snps=exposure_dat$SNP,
                                outcomes ="finn-b-D3_SARCOIDOSIS",
                                proxies = T, 
                                rsq = 0.8, 
                                align_alleles = 1, 
                                palindromes = 1, 
                                maf_threshold = 0.42)
outcome2 <- extract_outcome_data(snps=exposure_dat2$SNP,
                                outcomes ="ieu-b-4879",
                                proxies = T, 
                                rsq = 0.8, 
                                align_alleles = 1, 
                                palindromes = 1, 
                                maf_threshold = 0.42)
# # Data harmonization
data <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome)
data2 <- harmonise_data(exposure_dat = exposure_dat2, outcome_dat = outcome2)
### Mendelian Randomization Analysis ###
# Conduct MR analysis
mrresult <- mr(data, method_list=c('mr_ivw','mr_egger_regression','mr_weighted_median', 'mr_wald_ratio','mr_raps','mr_two_sample_ml'))
mrresult2 <- mr(data2, method_list=c('mr_ivw','mr_egger_regression','mr_weighted_median', 'mr_wald_ratio','mr_raps','mr_two_sample_ml'))
# Generating Odds Ratios
mrTab <- generate_odds_ratios(mrresult)
mrTab2 <- generate_odds_ratios(mrresult2)
#F statistics & R2
MR_F <- function(sample_size,num_IVs,r_square){numberator <- r_square * (sample_size - 1 - num_IVs)
                                               denominator <- (1 - r_square) * num_IVs
                                               f <- numberator / denominator
                                               return(f)}
my_PVE <- function(beta, se_beta, maf, N){numberator <- 2 * beta^2 * maf * (1 - maf)
                                          denominator <- 2 * beta^2 * maf * (1 - maf) + se_beta^2 * 2 * N * maf * (1 - maf)
                                          return(numberator/denominator)}
exposure_dat$pve <- my_PVE(beta = exposure_dat$beta.exposure, 
                           se_beta = exposure_dat$se.exposure,
                           maf = exposure_dat$eaf.exposure,
                           N = exposure_dat$samplesize.exposure)
exposure_dat2$pve <- my_PVE(beta = exposure_dat2$beta.exposure, 
                           se_beta = exposure_dat2$se.exposure,
                           maf = exposure_dat2$eaf.exposure,
                           N = 2046+215712)
IV_infor <- exposure_dat %>% group_by(exposure) %>%  summarise(samplesize=mean(samplesize.exposure),
                                                               nIV = n(),
                                                               R2=sum(pve))
IV_infor2 <- exposure_dat2 %>% group_by(exposure) %>%  summarise(samplesize=mean(samplesize.exposure),
                                                               nIV = n(),
                                                               R2=sum(pve))
IV_infor$F_stat <- MR_F(sample_size = IV_infor$samplesize, num_IVs = IV_infor$nIV, r_square = IV_infor$R2)
IV_infor2$F_stat <- MR_F(sample_size =2046+215712, num_IVs = IV_infor2$nIV, r_square = IV_infor2$R2)
# Single SNP MR analysis
res_single <- mr_singlesnp(data)
mr_funnel_plot(singlesnp_results = res_single)
res_single2 <- mr_singlesnp(data2)
mr_funnel_plot(singlesnp_results = res_single2)
# Heterogeneity Test
heterTab <- mr_heterogeneity(dat=data)
heterTab2 <- mr_heterogeneity(dat=data2)
# Pleiotropy Test
pleioTab <- mr_pleiotropy_test(data)
pleioTab2 <- mr_pleiotropy_test(data2)
#MR PRESSO test
set.seed(1234)
res_presso <- run_mr_presso(data, NbDistribution = 2000, SignifThreshold = 0.05)
set.seed(123)
res_presso2 <- run_mr_presso(data2, NbDistribution = 2000, SignifThreshold = 0.05)
#Re-run MR PRESSO after removing outliers if outlier exist
# names(res_presso) <- paste0(attributes(res_presso)$exposure," on ",attributes(res_presso)$outcome)
# res_presso_summary <- list()
# for (item in names(res_presso)){
#   global_P=res_presso[[item]]$`MR-PRESSO results`$`Global Test`$Pvalue
#   outlier_idx <- res_presso[[item]]$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`
#   snp_idx <- as.numeric(rownames(res_presso[[item]]$`MR-PRESSO results`$`Outlier Test`)[outlier_idx])
#   outlier_snp <- data$SNP[snp_idx]
#   res_presso_summary[[item]] <- list(global_P=global_P,
#                                      snp_idx=snp_idx,
#                                      outlier_snp=outlier_snp)
# }
# res_presso_summary <- as.data.frame(do.call(rbind, res_presso_summary))

# outliers <- unlist(res_presso_summary$snp_idx)
# outliers <- outliers[!is.na(outliers)]
# set.seed(1234)
# res_presso_2 <- run_mr_presso(data[-outliers,], NbDistribution = 4000, SignifThreshold = 0.05)
# ames(res_presso_2) <- paste0(attributes(res_presso_2)$exposure," on ",attributes(res_presso_2)$outcome)
out <- directionality_test(data)
out2 <- directionality_test(data2)
## Leave one out analysis
res_leaveoneout <- mr_leaveoneout(data)
mr_leaveoneout_plot(res_leaveoneout)
res_leaveoneout2 <- mr_leaveoneout(data2)
mr_leaveoneout_plot(res_leaveoneout2)
write.table(res_leaveoneout, file="forward.csv",sep = ",")
write.table(res_leaveoneout2, file="backward.csv",sep = ",")
#p_leaveoneout$`Asthma_childhood<12.bipolar`
#scatter plot 
mr_scatter_plot(mr_results=mrresult, dat=data)
mr_scatter_plot(mr_results=mrresult2, dat=data2)
#p_scatter$`Asthma_childhood<12.bipolar`
## Forest plot
library(forestploter)
df_forest <- mrTab[,c("outcome", "exposure", "method", "nsnp", "pval", "or", "or_lci95", "or_uci95")]
df_forest2 <- mrTab2[,c("outcome", "exposure", "method", "nsnp", "pval", "or", "or_lci95", "or_uci95")]
df_forest <- df_forest %>%
   mutate(outcome = ifelse(row_number() %% 5 == 1, "Sarcoidosis", ""),
          exposure = ifelse(row_number() %% 5 == 1,"TL", ""),
          nsnp = ifelse(row_number() %% 1 == 0, nsnp, ""),
          ` ` = paste(rep(" ", 20), collapse = " "),
          `OR (95% CI)` = sprintf("%.3f (%.3f to %.3f)", or, or_lci95, or_uci95),
          pval =sprintf("%.4f", pval))
df_forest2 <- df_forest2 %>%
  mutate(outcome = ifelse(row_number() %% 5 == 1,"TL", ""),
         exposure = ifelse(row_number() %% 5 == 1,"Sarcoidosis", ""),
         nsnp = ifelse(row_number() %% 1 == 0, nsnp, ""),
         ` ` = paste(rep(" ", 20), collapse = " "),
         `OR (95% CI)` = sprintf("%.3f (%.3f to %.3f)", or, or_lci95, or_uci95),
         pval =sprintf("%.4f", pval))
df<-rbind(df_forest, df_forest2)
tm <- forest_theme(base_size = 10,
                   # Confidence interval point shape, line type/color/width
                   ci_pch = 15,
                   ci_col = "blue",
                   ci_fill = "black",
                   ci_alpha = 0.8,
                   ci_lty = 1,
                   ci_lwd = 1.5,
                   ci_Theight = 0.2, # Set an T end at the end of CI 
                   # Reference line width/type/color
                   refline_lwd = 1,
                   refline_lty = "dashed",
                   refline_col = "red")
pt <- forest(df[,c("exposure","outcome","nsnp","method"," ","OR (95% CI)","pval")], 
             est = df$or,
             lower = df$or_lci95, 
             upper = df$or_uci95, 
             ci_column = 5,
             ref_line = 1, 
             xlim = c(0.35,1.02),
             ticks_at = c(0.35,0.7,0.9,1.02),
             theme = tm)
plot(pt)
write.table(data,"data_forward.csv",sep = ",", row.names = FALSE, col.names = TRUE)
write.table(data2,"data_reverse.csv",sep = ",", row.names = FALSE, col.names = TRUE)
write.table(exposure_dat,"data_exposure.csv",sep = ",", row.names = FALSE, col.names = TRUE)
write.table(exposure_dat2,"data_exposure2.csv",sep = ",", row.names = FALSE, col.names = TRUE)