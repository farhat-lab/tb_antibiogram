#read full and subsetted gentb output

gentb_output <- read.csv("gentb_binary_4.9.2021.txt", stringsAsFactors = FALSE, row.names = 1)

gentb_output_snp10 <- read.csv("gentb_subset_snp10_binary_4.9.21.txt", stringsAsFactors = FALSE, row.names = 1)

#Get geo metadata
geog <- read.delim("geog_cleaned_4.6.21.txt", sep = "\t", stringsAsFactors = FALSE)

#Get country names
countries <- sort(unique(geog$isolation_country))

#Get number for JAGS input files and save
jags_input_rr_df <- vector()
for (i in countries) {
    country_isolates <- na.omit(unlist(geog[which(geog$isolation_country == i),"BioSample"]))
    country_gentb <- na.omit(gentb_out[country_isolates,])
    abx_r_sum <- apply(country_gentb, 2, sum)
    country_gentb_rr <- country_gentb[which(country_gentb$rif == 1),]
    abx_r_rr_sum <- apply(country_gentb_rr, 2, sum)
    country_gentb_rs <- country_gentb[which(country_gentb$rif == 0),]
    abx_r_rs_sum <- apply(country_gentb_rs, 2, sum)
    country_isolates_snp10 <- na.omit(unlist(geog[which(geog$isolation_country == i),"BioSample"]))
    country_gentb_snp10 <- na.omit(gentb_out_snp10[country_isolates,])
    abx_r_sum_snp10 <- apply(country_gentb_snp10, 2, sum)
    country_gentb_snp10_rr <- country_gentb_snp10[which(country_gentb_snp10$rif == 1),]
    abx_r_rr_sum_snp10 <- apply(country_gentb_snp10_rr, 2, sum)
    country_gentb_snp10_rs <- country_gentb_snp10[which(country_gentb_snp10$rif == 0),]
    abx_r_rs_sum_snp10 <- apply(country_gentb_snp10_rs, 2, sum)
    jags_input_country_rr <- cbind(i, names(abx_r_sum), nrow(country_gentb), abx_r_sum, nrow(country_gentb_rr), abx_r_rr_sum, nrow(country_gentb_rs), abx_r_rs_sum, nrow(country_gentb_snp10), abx_r_sum_snp10, nrow(country_gentb_snp10_rr), abx_r_rr_sum_snp10, nrow(country_gentb_snp10_rs), abx_r_rs_sum_snp10)
    jags_input_rr_df <- rbind(jags_input_rr_df, jags_input_country_rr)
}

colnames(jags_input_rr_df) <- c("country", "drug", "gentb_n", "gentb_res","gentb_n_rr", "gentb_res_rr", "gentb_n_rs", "gentb_res_rs", "gentb_snp10_n", "gentb_snp10_res","gentb_snp10_n_rr", "gentb_snp10_res_rr", "gentb_snp10_n_rs", "gentb_snp10_res_rs")

write.table(jags_input_rr_df, "jags_model_input_w_rr_4.9.21.tsv", quote = FALSE, row.names = FALSE, sep = "\t")

jags_input_mdr_df <- vector()
for (i in countries) {
     country_isolates <- na.omit(unlist(geog[which(geog$isolation_country == i),"BioSample"]))
    country_gentb <- na.omit(gentb_out[country_isolates,])
    abx_r_sum <- apply(country_gentb, 2, sum)
    country_gentb_mdr <- country_gentb[which(country_gentb$rif == 1 & country_gentb$inh == 1),]
    abx_r_mdr_sum <- apply(country_gentb_mdr, 2, sum)
    country_gentb_nonmdr <- country_gentb[which(country_gentb$rif == 0 & country_gentb$inh == 0),]
    abx_r_nonmdr_sum <- apply(country_gentb_nonmdr, 2, sum)
    country_isolates_snp10 <- na.omit(unlist(geog[which(geog$isolation_country == i),"BioSample"]))
    country_gentb_snp10 <- na.omit(gentb_out_snp10[country_isolates,])
    abx_r_sum_snp10 <- apply(country_gentb_snp10, 2, sum)
    country_gentb_snp10_mdr <- country_gentb_snp10[which(country_gentb_snp10$rif == 1 & country_gentb_snp10$inh == 1),]
    abx_r_mdr_sum_snp10 <- apply(country_gentb_snp10_mdr, 2, sum)
    country_gentb_snp10_nonmdr <- country_gentb_snp10[which(country_gentb_snp10$rif == 0 & country_gentb_snp10$inh == 0),]
    abx_r_nonmdr_sum_snp10 <- apply(country_gentb_snp10_nonmdr, 2, sum)
    jags_input_country_mdr <- cbind(i, names(abx_r_sum), nrow(country_gentb), abx_r_sum, nrow(country_gentb_mdr), abx_r_mdr_sum, nrow(country_gentb_nonmdr), abx_r_nonmdr_sum, nrow(country_gentb_snp10), abx_r_sum_snp10, nrow(country_gentb_snp10_mdr), abx_r_mdr_sum_snp10, nrow(country_gentb_snp10_nonmdr), abx_r_nonmdr_sum_snp10)
    jags_input_mdr_df <- rbind(jags_input_mdr_df, jags_input_country_mdr)
}    

colnames(jags_input_mdr_df) <- c("country", "drug", "gentb_n", "gentb_res","gentb_n_mdr", "gentb_res_mdr", "gentb_n_nonmdr", "gentb_res_nonmdr", "gentb_snp10_n", "gentb_snp10_res","gentb_snp10_n_mdr", "gentb_snp10_res_mdr", "gentb_snp10_n_nonmdr", "gentb_snp10_res_nonmdr")

write.table(jags_input_mdr_df, "jags_model_input_w_mdr_4.9.21.tsv", quote = FALSE, row.names = FALSE, sep = "\t")