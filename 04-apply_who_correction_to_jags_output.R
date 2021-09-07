#Read JAGS output and subset to estimates with and without outbreak correction and exclude ciprofloxacin
jags_out <- read.csv("jags_model_output_w_inh_mono_7.20.2021.tsv", sep = "\t", stringsAsFactors = FALSE, header = TRUE)
jags_out <- jags_out[-which(jags_out$drug == "cip"),]

jags_out_snp10 <- jags_out[,c("country", "drug", colnames(jags_out)[grep("snp10", colnames(jags_out))])]

jags_out <- jags_out[,colnames(jags_out)[-grep("snp10", colnames(jags_out))]]

#Read WHO composite RR
who_composite_rr <- read.csv("who_tb_data/composite_rr_w_var_simplified_2018.txt", sep = "\t", header = TRUE, row.names = 1)

#Get geo metadata
geog <- read.csv("geog_cleaned_4.6.21.txt", sep = "\t", stringsAsFactors = FALSE)

#get antibiotic and country names
countries <- unique(jags_out$country)
antibiotics <- unique(jags_out$drug)

#Calculate the corrected estimates using the marginal probability formula described in Dixit et al.

who_corrected_abx_r_df <- vector()
for (country in countries) {
    jags_out_country <- jags_out[which(jags_out$country == country),]
    for (abx in antibiotics) {
        p_a_rr <- jags_out_country[which(jags_out_country$drug == abx), "mean_rr"]
        p_a_rs <- jags_out_country[which(jags_out_country$drug == abx), "mean_rs"]
        p_rr <- who_composite_rr[paste(country), "composite_rr_pct"]/100 #number is in percentage form
        p_rs <- 1-p_rr
        p_a_r <- (p_a_rr*p_rr)+(p_a_rs*p_rs)
        sd_a_rr <- jags_out_country[which(jags_out_country$drug == abx), "sd_rr"]
        var_a_rr <- sd_a_rr*sd_a_rr
        var_rr <- who_composite_rr[paste(country), "composite_rr_var"]
        var_a_r <- p_rr*p_rr*var_a_rr + p_a_rr*p_a_rr*var_rr
        who_corrected_abx_r <- cbind(jags_out_country[which(jags_out_country$drug == abx),], abx, p_a_r, var_a_r)
        who_corrected_abx_r_df <- rbind(who_corrected_abx_r_df, who_corrected_abx_r)
    }
}

write.csv(who_corrected_abx_r_df, "who_corrected_7.20.21.csv", row.names = FALSE)

who_corrected_snp10_abx_r_df <- vector()
for (country in countries) {
    jags_out_snp10_country <- jags_out_snp10[which(jags_out_snp10$country == country),]
    for (abx in antibiotics) {
        p_a_rr <- jags_out_snp10_country[which(jags_out_snp10_country$drug == abx), "mean_snp10_rr"]
        p_a_rs <- jags_out_snp10_country[which(jags_out_snp10_country$drug == abx), "mean_snp10_rs"]
        p_rr <- who_composite_rr[paste(country), "composite_rr_pct"]/100
        p_rs <- 1-p_rr
        p_a_r <- (p_a_rr*p_rr)+(p_a_rs*p_rs)
        sd_a_rr <- jags_out_snp10_country[which(jags_out_snp10_country$drug == abx), "sd_snp10_rr"]
        var_a_rr <- sd_a_rr*sd_a_rr
        var_rr <- who_composite_rr[paste(country), "composite_rr_var"]
        var_a_r <- p_rr*p_rr*var_a_rr + p_a_rr*p_a_rr*var_rr
        who_corrected_snp10_abx_r <- cbind(jags_out_snp10_country[which(jags_out_snp10_country$drug == abx),], abx, p_a_r, var_a_r)
        who_corrected_snp10_abx_r_df <- rbind(who_corrected_snp10_abx_r_df, who_corrected_snp10_abx_r)
    }
}

write.csv(who_corrected_snp10_abx_r_df, "who_corrected_snp10_7.20.21.csv", row.names = FALSE)

#Plot the antibiograms into PDF files. Note that this includes all countries not restricted by sample size. 

abx_order <- c("rif", "inh", "pza" ,"emb", "eth","str", "cap", "amk", "kan", "oflx", "levo", "pas")

pdf("Antibiogram_jags_who_corrected_8.25.21.pdf")
for (i in 1:length(countries)) {
    main_title <- paste(countries[i], "\n(n = ", unique(who_corrected_abx_r_df[which(who_corrected_abx_r_df$country == countries[i]), "gentb_n"]), ")")
   print(ggplot(who_corrected_abx_r_df, aes(x = factor(abx, level = abx_order ), y = p_a_r*100, group = country)) +
    geom_bar(stat="identity", color="maroon", fill = "maroon") + 
    geom_errorbar(aes(ymin=ifelse((p_a_r-1.96*sqrt(var_a_r))*100 >0,(p_a_r-1.96*sqrt(var_a_r))*100,0), ymax=ifelse((p_a_r+1.96*sqrt(var_a_r))*100 < 100, (p_a_r+1.96*sqrt(var_a_r))*100, 100)))  + 
    facet_wrap_paginate(~country, ncol = 1, nrow = 1, page = i) + ggtitle(main_title) +
    xlab("Antibiotic") + ylab("% Resistant") + theme_minimal())
}
dev.off()

countries <- sort(unique(who_corrected_snp10_abx_r_df$country))

pdf("Antibiogram_jags_who_corrected_snp10_8.25.21.pdf")
for (i in 1:length(countries)) {
    main_title <- paste(countries[i], "\n(n (>10 SNP apart)= ", unique(who_corrected_snp10_abx_r_df[which(who_corrected_snp10_abx_r_df$country == countries[i]), "gentb_snp10_n"]), ")")
   print(ggplot(who_corrected_snp10_abx_r_df, aes(x = factor(abx, level = abx_order ), y = p_a_r*100, group = country)) +
    geom_bar(stat="identity", color="maroon", fill = "maroon") + 
    geom_errorbar(aes(ymin=ifelse((p_a_r-1.96*sqrt(var_a_r))*100 >0,(p_a_r-1.96*sqrt(var_a_r))*100,0), ymax=ifelse((p_a_r+1.96*sqrt(var_a_r))*100 < 100, (p_a_r+1.96*sqrt(var_a_r))*100, 100)))  + 
    facet_wrap_paginate(~country, ncol = 1, nrow = 1, page = i) + ggtitle(main_title) +
    xlab("Antibiotic") + ylab("% Resistant") + theme_minimal())
}
dev.off()