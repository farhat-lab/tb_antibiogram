#install.packages("rjags")
library(rjags)

##alpha and beta calculation 
#alpha <- function(mean, var) {
#  if (var < (mean*(1-mean))) {
#    mean*(((mean*(1-mean))/var)-1)
#  } else {
#    "condition not met"
#  }
#}

#beta <- function(mean, var) {
#  if (var < (mean*(1-mean))) {
#    (1-mean)*(((mean*(1-mean))/var)-1)
#  } else {
#    "condition not met"
#  }
#}

se_sp <- read.csv("../bias_correction/se_sp_alpha_beta_w_inh_mono_shrt_crs.csv", row.names = 1, header = TRUE)

# 
# se_sp <- transform(se_sp, 
#                    se_alpha = apply(se_sp, 1, function(x) alpha(mean = x["TPR_bootstrapped"], var = x["sens_var"])), 
#                    se_beta = apply(se_sp, 1, function(x) beta(mean = x["TPR_bootstrapped"], var = x["sens_var"])),
#                    sp_alpha = apply(se_sp, 1, function(x) alpha(mean = x["TNR_bootstrapped"], var = x["spec_var"])),
#                    sp_beta = apply(se_sp, 1, function(x) beta(mean = x["TNR_bootstrapped"], var = x["spec_var"])))
# 
# write.csv(se_sp, "se_sp_alpha_beta.csv", quote = FALSE)
#       

#se_sp["inh_mono",] <- transform(se_sp["inh_mono",], sens_sd = (TPR_hiCI-TPR_lowCI)/3.92, spec_sd = (TNR_hiCI-TNR_lowCI)/3.92)
#se_sp <- transform(se_sp, sens_var = sens_sd*sens_sd, spec_var = spec_sd*spec_sd)

#se_sp <- transform(se_sp, 
#                  se_alpha = apply(se_sp, 1, function(x) alpha(mean = x["TPR_bootstrapped"], var = x["sens_var"])),
#                  se_beta = apply(se_sp, 1, function(x) beta(mean = x["TPR_bootstrapped"], var = x["sens_var"])),
#                  sp_alpha = apply(se_sp, 1, function(x) alpha(mean = x["TNR_bootstrapped"], var = x["spec_var"])),
#                  sp_beta = apply(se_sp, 1, function(x) beta(mean = x["TNR_bootstrapped"], var = x["spec_var"])))

#se_sp["shrt_crs",] <- transform(se_sp["shrt_crs",], sens_sd = (TPR_hiCI-TPR_lowCI)/3.92, spec_sd = (TNR_hiCI-TNR_lowCI)/3.92)
#se_sp["shrt_crs",] <- transform(se_sp["shrt_crs",], sens_var = sens_sd*sens_sd, spec_var = spec_sd*spec_sd)
#se_sp["shrt_crs",] <- transform(se_sp["shrt_crs",], 
#                   se_alpha = apply(se_sp["shrt_crs",], 1, function(x) alpha(mean = x["TPR"], var = x["sens_var"])),
#                   se_beta = apply(se_sp["shrt_crs",], 1, function(x) beta(mean = x["TPR"], var = x["sens_var"])),
#                   sp_alpha = apply(se_sp["shrt_crs",], 1, function(x) alpha(mean = x["TNR"], var = x["spec_var"])),
#                   sp_beta = apply(se_sp["shrt_crs",], 1, function(x) beta(mean = x["TNR"], var = x["spec_var"])))
#write.csv(se_sp, "se_sp_alpha_beta_w_inh_mono_shrt_crs.csv", quote = FALSE)


model_code = ' model {
  tpos ~ dbin(theta, sample_size)
  theta <- se*phi + (1-sp)*(1-phi)
  se ~ dbeta(se_alpha, se_beta) 
  sp ~ dbeta(sp_alpha, sp_beta)
  phi ~ dbeta(1, 1) 
}
'

input_file <- read.delim("jags_model_input_w_rr_4.9.21.tsv", stringsAsFactors = FALSE, header = TRUE, sep = "\t")

abx <- rownames(se_sp)[-14]

countries <- unique(input_file$country)


set.seed(123)
for(country in countries) {
  country_r <- input_file[which(input_file$country == country),]
  for(a in abx) {
    print(c("now calculating", paste(country), paste(a)))
    abx_se_sp <- se_sp[a,]
    country_abx_r <- country_r[which(country_r$drug == a), ]
    set.seed(123)
    jags_test <- jags.model(file = textConnection(model_code),
                            data = list(tpos = country_abx_r$gentb_res, 
                                        sample_size = country_abx_r$gentb_n, 
                                        se_alpha = abx_se_sp$se_alpha, 
                                        se_beta = abx_se_sp$se_beta, 
                                        sp_alpha = abx_se_sp$sp_alpha, 
                                        sp_beta = abx_se_sp$sp_beta), 
                            n.chains = 3, quiet = TRUE)
    update(jags_test, 10000)
    
    mcmc <- coda.samples(jags_test, variable.names=c("phi"), n.iter=20000)
    omcmc <- summary(mcmc)
    input_file[which(input_file$country == country & input_file$drug == a), "mean"] <- omcmc$statistics[1]
    input_file[which(input_file$country == country & input_file$drug == a), "sd"] <- omcmc$statistics[2]
    input_file[which(input_file$country == country & input_file$drug == a), "lo"] <- omcmc$quantiles[1]
    input_file[which(input_file$country == country & input_file$drug == a), "hi"] <- omcmc$quantiles[5]
    
    jags_test_snp10 <- jags.model(file = textConnection(model_code),
                            data = list(tpos = country_abx_r$gentb_snp10_res, 
                                        sample_size = country_abx_r$gentb_snp10_n, 
                                        se_alpha = abx_se_sp$se_alpha, 
                                        se_beta = abx_se_sp$se_beta, 
                                        sp_alpha = abx_se_sp$sp_alpha, 
                                        sp_beta = abx_se_sp$sp_beta), 
                            n.chains = 3, quiet = TRUE)
    update(jags_test_snp10, 10000)
    
    mcmc_snp10 <- coda.samples(jags_test_snp10, variable.names=c("phi"), n.iter=20000)
    omcmc_snp10 <- summary(mcmc_snp10)
    input_file[which(input_file$country == country & input_file$drug == a), "mean_snp10"] <- omcmc_snp10$statistics[1]
    input_file[which(input_file$country == country & input_file$drug == a), "sd_snp10"] <- omcmc_snp10$statistics[2]
    input_file[which(input_file$country == country & input_file$drug == a), "lo_snp10"] <- omcmc_snp10$quantiles[1]
    input_file[which(input_file$country == country & input_file$drug == a), "hi_snp10"] <- omcmc_snp10$quantiles[5]
  }
}

write.table(input_file, "jags_model_output.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

set.seed(123)
for(country in countries) {
  country_r <- input_file[which(input_file$country == country),]
  for(a in abx) {
    print(c("now calculating", paste(country), paste(a)))
    abx_se_sp <- se_sp[a,]
    country_abx_r <- country_r[which(country_r$drug == a), ]
    
    #JAGS model for all RR
    set.seed(123)
    jags_rr <- jags.model(file = textConnection(model_code),
                            data = list(tpos = country_abx_r$gentb_res_rr, 
                                        sample_size = country_abx_r$gentb_n_rr, 
                                        se_alpha = abx_se_sp$se_alpha, 
                                        se_beta = abx_se_sp$se_beta, 
                                        sp_alpha = abx_se_sp$sp_alpha, 
                                        sp_beta = abx_se_sp$sp_beta), 
                            n.chains = 3, quiet = TRUE)
    update(jags_rr, 10000)
    
    mcmc_rr <- coda.samples(jags_rr, variable.names=c("phi"), n.iter=20000)
    omcmc_rr <- summary(mcmc_rr)
    input_file[which(input_file$country == country & input_file$drug == a), "mean_rr"] <- omcmc_rr$statistics[1]
    input_file[which(input_file$country == country & input_file$drug == a), "sd_rr"] <- omcmc_rr$statistics[2]
    input_file[which(input_file$country == country & input_file$drug == a), "lo_rr"] <- omcmc_rr$quantiles[1]
    input_file[which(input_file$country == country & input_file$drug == a), "hi_rr"] <- omcmc_rr$quantiles[5]
    
    #JAGS model for all RS
    set.seed(123)
    jags_rs <- jags.model(file = textConnection(model_code),
                          data = list(tpos = country_abx_r$gentb_res_rs, 
                                      sample_size = country_abx_r$gentb_n_rs, 
                                      se_alpha = abx_se_sp$se_alpha, 
                                      se_beta = abx_se_sp$se_beta, 
                                      sp_alpha = abx_se_sp$sp_alpha, 
                                      sp_beta = abx_se_sp$sp_beta), 
                          n.chains = 3, quiet = TRUE)
    update(jags_rs, 10000)
    
    mcmc_rs <- coda.samples(jags_rs, variable.names=c("phi"), n.iter=20000)
    omcmc_rs <- summary(mcmc_rs)
    input_file[which(input_file$country == country & input_file$drug == a), "mean_rs"] <- omcmc_rs$statistics[1]
    input_file[which(input_file$country == country & input_file$drug == a), "sd_rs"] <- omcmc_rs$statistics[2]
    input_file[which(input_file$country == country & input_file$drug == a), "lo_rs"] <- omcmc_rs$quantiles[1]
    input_file[which(input_file$country == country & input_file$drug == a), "hi_rs"] <- omcmc_rs$quantiles[5]
    
    #JAGS model for snp10 subset RR
    set.seed(123)
    jags_snp10_rr <- jags.model(file = textConnection(model_code),
                                  data = list(tpos = country_abx_r$gentb_snp10_res_rr, 
                                              sample_size = country_abx_r$gentb_snp10_n_rr, 
                                              se_alpha = abx_se_sp$se_alpha, 
                                              se_beta = abx_se_sp$se_beta, 
                                              sp_alpha = abx_se_sp$sp_alpha, 
                                              sp_beta = abx_se_sp$sp_beta), 
                                  n.chains = 3, quiet = TRUE)
    update(jags_snp10_rr, 10000)
    
    mcmc_snp10_rr <- coda.samples(jags_snp10_rr, variable.names=c("phi"), n.iter=20000)
    omcmc_snp10_rr <- summary(mcmc_snp10_rr)
    input_file[which(input_file$country == country & input_file$drug == a), "mean_snp10_rr"] <- omcmc_snp10_rr$statistics[1]
    input_file[which(input_file$country == country & input_file$drug == a), "sd_snp10_rr"] <- omcmc_snp10_rr$statistics[2]
    input_file[which(input_file$country == country & input_file$drug == a), "lo_snp10_rr"] <- omcmc_snp10_rr$quantiles[1]
    input_file[which(input_file$country == country & input_file$drug == a), "hi_snp10_rr"] <- omcmc_snp10_rr$quantiles[5]
    
    #JAGS model for snp10 subset RS
    set.seed(123)
    jags_snp10_rs <- jags.model(file = textConnection(model_code),
                                data = list(tpos = country_abx_r$gentb_snp10_res_rs, 
                                            sample_size = country_abx_r$gentb_snp10_n_rs, 
                                            se_alpha = abx_se_sp$se_alpha, 
                                            se_beta = abx_se_sp$se_beta, 
                                            sp_alpha = abx_se_sp$sp_alpha, 
                                            sp_beta = abx_se_sp$sp_beta), 
                                n.chains = 3, quiet = TRUE)
    update(jags_snp10_rs, 10000)
    
    mcmc_snp10_rs <- coda.samples(jags_snp10_rs, variable.names=c("phi"), n.iter=20000)
    omcmc_snp10_rs <- summary(mcmc_snp10_rs)
    input_file[which(input_file$country == country & input_file$drug == a), "mean_snp10_rs"] <- omcmc_snp10_rs$statistics[1]
    input_file[which(input_file$country == country & input_file$drug == a), "sd_snp10_rs"] <- omcmc_snp10_rs$statistics[2]
    input_file[which(input_file$country == country & input_file$drug == a), "lo_snp10_rs"] <- omcmc_snp10_rs$quantiles[1]
    input_file[which(input_file$country == country & input_file$drug == a), "hi_snp10_rs"] <- omcmc_snp10_rs$quantiles[5]
    
  }
}

write.table(input_file, "jags_model_rr_output.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

#Change bias correction for INH MONO-R

set.seed(123)
for(country in countries) {
  country_r <- input_file[which(input_file$country == country),]
    print(c("now calculating", paste(country)))
    abx_se_sp <- se_sp["inh_mono",]
    country_abx_r <- country_r[which(country_r$drug == "inh"), ]
    
    #JAGS model for INH-R RS
    set.seed(123)
    jags_rs <- jags.model(file = textConnection(model_code),
                          data = list(tpos = country_abx_r$gentb_res_rs, 
                                      sample_size = country_abx_r$gentb_n_rs, 
                                      se_alpha = abx_se_sp$se_alpha, 
                                      se_beta = abx_se_sp$se_beta, 
                                      sp_alpha = abx_se_sp$sp_alpha, 
                                      sp_beta = abx_se_sp$sp_beta), 
                          n.chains = 3, quiet = TRUE)
    update(jags_rs, 10000)
    
    mcmc_rs <- coda.samples(jags_rs, variable.names=c("phi"), n.iter=20000)
    omcmc_rs <- summary(mcmc_rs)
    input_file[which(input_file$country == country & input_file$drug == "inh"), "mean_rs"] <- omcmc_rs$statistics[1]
    input_file[which(input_file$country == country & input_file$drug == "inh"), "sd_rs"] <- omcmc_rs$statistics[2]
    input_file[which(input_file$country == country & input_file$drug == "inh"), "lo_rs"] <- omcmc_rs$quantiles[1]
    input_file[which(input_file$country == country & input_file$drug == "inh"), "hi_rs"] <- omcmc_rs$quantiles[5]
    
    #JAGS model for snp10 subset RS
    set.seed(123)
    jags_snp10_rs <- jags.model(file = textConnection(model_code),
                                data = list(tpos = country_abx_r$gentb_snp10_res_rs, 
                                            sample_size = country_abx_r$gentb_snp10_n_rs, 
                                            se_alpha = abx_se_sp$se_alpha, 
                                            se_beta = abx_se_sp$se_beta, 
                                            sp_alpha = abx_se_sp$sp_alpha, 
                                            sp_beta = abx_se_sp$sp_beta), 
                                n.chains = 3, quiet = TRUE)
    update(jags_snp10_rs, 10000)
    
    mcmc_snp10_rs <- coda.samples(jags_snp10_rs, variable.names=c("phi"), n.iter=20000)
    omcmc_snp10_rs <- summary(mcmc_snp10_rs)
    input_file[which(input_file$country == country & input_file$drug == "inh"), "mean_snp10_rs"] <- omcmc_snp10_rs$statistics[1]
    input_file[which(input_file$country == country & input_file$drug == "inh"), "sd_snp10_rs"] <- omcmc_snp10_rs$statistics[2]
    input_file[which(input_file$country == country & input_file$drug == "inh"), "lo_snp10_rs"] <- omcmc_snp10_rs$quantiles[1]
    input_file[which(input_file$country == country & input_file$drug == "inh"), "hi_snp10_rs"] <- omcmc_snp10_rs$quantiles[5]
}

write.table(input_file, "jags_model_output_w_inh_mono_7.20.2021.tsv",  sep = "\t", quote = FALSE, row.names = FALSE)

input_file_mdr <- read.delim("jags_model_input_w_mdr_4.9.21.tsv", stringsAsFactors = FALSE, header = TRUE, sep = "\t")

for(country in countries) {
  country_mdr <- input_file_mdr[which(input_file_mdr$country == country),]
  for(a in abx) {
    print(c("now calculating", paste(country), paste(a)))
    abx_se_sp <- se_sp[a,]
    country_abx_mdr <- country_mdr[which(country_mdr$drug == a), ]
    
    #JAGS model for all MDR
    jags_mdr <- jags.model(file = textConnection(model_code),
                          data = list(tpos = country_abx_mdr$gentb_res_mdr, 
                                      sample_size = country_abx_mdr$gentb_n_mdr, 
                                      se_alpha = abx_se_sp$se_alpha, 
                                      se_beta = abx_se_sp$se_beta, 
                                      sp_alpha = abx_se_sp$sp_alpha, 
                                      sp_beta = abx_se_sp$sp_beta), 
                          n.chains = 3, quiet = TRUE)
    update(jags_mdr, 10000)
    
    mcmc_mdr <- coda.samples(jags_mdr, variable.names=c("phi"), n.iter=20000)
    omcmc_mdr <- summary(mcmc_mdr)
    input_file_mdr[which(input_file_mdr$country == country & input_file_mdr$drug == a), "mean_mdr"] <- omcmc_mdr$statistics[1]
    input_file_mdr[which(input_file_mdr$country == country & input_file_mdr$drug == a), "sd_mdr"] <- omcmc_mdr$statistics[2]
    input_file_mdr[which(input_file_mdr$country == country & input_file_mdr$drug == a), "lo_mdr"] <- omcmc_mdr$quantiles[1]
    input_file_mdr[which(input_file_mdr$country == country & input_file_mdr$drug == a), "hi_mdr"] <- omcmc_mdr$quantiles[5]
    
    # #JAGS model for all Non-MDR
    # jags_nonmdr <- jags.model(file = textConnection(model_code),
    #                       data = list(tpos = country_abx_mdr$gentb_res_nonmdr, 
    #                                   sample_size = country_abx_mdr$gentb_n_nonmdr, 
    #                                   se_alpha = abx_se_sp$se_alpha, 
    #                                   se_beta = abx_se_sp$se_beta, 
    #                                   sp_alpha = abx_se_sp$sp_alpha, 
    #                                   sp_beta = abx_se_sp$sp_beta), 
    #                       n.chains = 3, quiet = TRUE)
    # update(jags_nonmdr, 10000)
    # 
    # mcmc_nonmdr <- coda.samples(jags_nonmdr, variable.names=c("phi"), n.iter=20000)
    # omcmc_nonmdr <- summary(mcmc_nonmdr)
    # input_file_mdr[which(input_file_mdr$country == country & input_file_mdr$drug == a), "mean_nonmdr"] <- omcmc_nonmdr$statistics[1]
    # input_file_mdr[which(input_file_mdr$country == country & input_file_mdr$drug == a), "sd_nonmdr"] <- omcmc_nonmdr$statistics[2]
    # input_file_mdr[which(input_file_mdr$country == country & input_file_mdr$drug == a), "lo_nonmdr"] <- omcmc_nonmdr$quantiles[1]
    # input_file_mdr[which(input_file_mdr$country == country & input_file_mdr$drug == a), "hi_nonmdr"] <- omcmc_nonmdr$quantiles[5]
    
    #JAGS model for snp10 subset MDR
    jags_snp10_mdr <- jags.model(file = textConnection(model_code),
                                data = list(tpos = country_abx_mdr$gentb_snp10_res_mdr, 
                                            sample_size = country_abx_mdr$gentb_snp10_n_mdr, 
                                            se_alpha = abx_se_sp$se_alpha, 
                                            se_beta = abx_se_sp$se_beta, 
                                            sp_alpha = abx_se_sp$sp_alpha, 
                                            sp_beta = abx_se_sp$sp_beta), 
                                n.chains = 3, quiet = TRUE)
    update(jags_snp10_mdr, 10000)
    
    mcmc_snp10_mdr <- coda.samples(jags_snp10_mdr, variable.names=c("phi"), n.iter=20000)
    omcmc_snp10_mdr <- summary(mcmc_snp10_mdr)
    input_file_mdr[which(input_file_mdr$country == country & input_file_mdr$drug == a), "mean_snp10_mdr"] <- omcmc_snp10_mdr$statistics[1]
    input_file_mdr[which(input_file_mdr$country == country & input_file_mdr$drug == a), "sd_snp10_mdr"] <- omcmc_snp10_mdr$statistics[2]
    input_file_mdr[which(input_file_mdr$country == country & input_file_mdr$drug == a), "lo_snp10_mdr"] <- omcmc_snp10_mdr$quantiles[1]
    input_file_mdr[which(input_file_mdr$country == country & input_file_mdr$drug == a), "hi_snp10_mdr"] <- omcmc_snp10_mdr$quantiles[5]
    
    # #JAGS model for snp10 subset NON-MDR
    # jags_snp10_nonmdr <- jags.model(file = textConnection(model_code),
    #                             data = list(tpos = country_abx_mdr$gentb_snp10_res_nonmdr, 
    #                                         sample_size = country_abx_mdr$gentb_snp10_n_nonmdr, 
    #                                         se_alpha = abx_se_sp$se_alpha, 
    #                                         se_beta = abx_se_sp$se_beta, 
    #                                         sp_alpha = abx_se_sp$sp_alpha, 
    #                                         sp_beta = abx_se_sp$sp_beta), 
    #                             n.chains = 3, quiet = TRUE)
    # update(jags_snp10_nonmdr, 10000)
    # 
    # mcmc_snp10_nonmdr <- coda.samples(jags_snp10_nonmdr, variable.names=c("phi"), n.iter=20000)
    # omcmc_snp10_nonmdr <- summary(mcmc_snp10_nonmdr)
    # input_file_mdr[which(input_file_mdr$country == country & input_file_mdr$drug == a), "mean_snp10_nonmdr"] <- omcmc_snp10_nonmdr$statistics[1]
    # input_file_mdr[which(input_file_mdr$country == country & input_file_mdr$drug == a), "sd_snp10_nonmdr"] <- omcmc_snp10_nonmdr$statistics[2]
    # input_file_mdr[which(input_file_mdr$country == country & input_file_mdr$drug == a), "lo_snp10_nonmdr"] <- omcmc_snp10_nonmdr$quantiles[1]
    # input_file_mdr[which(input_file_mdr$country == country & input_file_mdr$drug == a), "hi_snp10_nonmdr"] <- omcmc_snp10_nonmdr$quantiles[5]
    
  }
}

write.table(input_file_mdr, "jags_model_mdr_output_4.9.21.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

#UK Levo test
uk_levo <- input_file[which(input_file$country == "United Kingdom of Great Britain and Northern Ireland" & input_file$drug == "levo"),]

uk_levo_jags <- jags.model(file = textConnection(model_code),
           data = list(tpos = uk_levo$gentb_res_rs, 
                       sample_size = uk_levo$gentb_snp10_n_rs, 
                       se_alpha = se_sp["levo",]$se_alpha, 
                       se_beta = se_sp["levo",]$se_beta, 
                       sp_alpha = se_sp["levo",]$sp_alpha, 
                       sp_beta = se_sp["levo",]$sp_beta), 
           n.chains = 3, quiet = TRUE)

update(uk_levo_jags, 10000)
mcmc_uk_levo_jags <- coda.samples(uk_levo_jags, variable.names=c("phi"), n.iter=20000)
omcmc_uk_levo_jags <- summary(mcmc_uk_levo_jags)

sa_input <- read.csv("south_africa_jags_input_4.9.21.csv", stringsAsFactors = FALSE)
sa_public <- sa_input[which(sa_input$isolate_type == "public"),]
sa_survey <- sa_input[which(sa_input$isolate_type == "survey"),]
sa_total <- sa_input[which(sa_input$isolate_type == "total"),]

set.seed(123)
for(isolate_type in c("public", "survey", "total")) {
  sa_subset<- sa_input[which(sa_input$isolate_type == isolate_type),]
  for(a in abx) {
    print(c("now calculating", paste(a)))
    abx_se_sp <- se_sp[a,]
    sa_subset_abx <- sa_subset[which(sa_subset$drug == a), ]
    
    #JAGS model for all RR
    sa_subset_jags_rr <- jags.model(file = textConnection(model_code),
                          data = list(tpos = sa_subset_abx$n_rr_abx_r, 
                                      sample_size = sa_subset_abx$n_rr, 
                                      se_alpha = abx_se_sp$se_alpha, 
                                      se_beta = abx_se_sp$se_beta, 
                                      sp_alpha = abx_se_sp$sp_alpha, 
                                      sp_beta = abx_se_sp$sp_beta), 
                          n.chains = 3, quiet = TRUE)
    update(sa_subset_jags_rr, 10000)
    
    mcmc_rr <- coda.samples(sa_subset_jags_rr, variable.names=c("phi"), n.iter=20000)
    omcmc_rr <- summary(mcmc_rr)
    sa_input[which(sa_input$isolate_type == isolate_type & sa_input$drug == a), "mean_rr"] <- omcmc_rr$statistics[1]
    sa_input[which(sa_input$isolate_type == isolate_type & sa_input$drug == a), "sd_rr"] <- omcmc_rr$statistics[2]
    sa_input[which(sa_input$isolate_type == isolate_type & sa_input$drug == a), "lo_rr"] <- omcmc_rr$quantiles[1]
    sa_input[which(sa_input$isolate_type == isolate_type & sa_input$drug == a), "hi_rr"] <- omcmc_rr$quantiles[5]
    
    #JAGS model for all RS
    sa_subset_jags_rs <- jags.model(file = textConnection(model_code),
                          data = list(tpos = sa_subset_abx$n_rs_abx_r, 
                                      sample_size = sa_subset_abx$n_rs, 
                                      se_alpha = abx_se_sp$se_alpha, 
                                      se_beta = abx_se_sp$se_beta, 
                                      sp_alpha = abx_se_sp$sp_alpha, 
                                      sp_beta = abx_se_sp$sp_beta), 
                          n.chains = 3, quiet = TRUE)
    update(sa_subset_jags_rs, 10000)
    
    mcmc_rs <- coda.samples(sa_subset_jags_rs, variable.names=c("phi"), n.iter=20000)
    omcmc_rs <- summary(mcmc_rs)
    sa_input[which(sa_input$isolate_type == isolate_type & sa_input$drug == a), "mean_rs"] <- omcmc_rs$statistics[1]
    sa_input[which(sa_input$isolate_type == isolate_type & sa_input$drug == a), "sd_rs"] <- omcmc_rs$statistics[2]
    sa_input[which(sa_input$isolate_type == isolate_type & sa_input$drug == a), "lo_rs"] <- omcmc_rs$quantiles[1]
    sa_input[which(sa_input$isolate_type == isolate_type & sa_input$drug == a), "hi_rs"] <- omcmc_rs$quantiles[5]

    #JAGS Model for MDR  
    sa_subset_jags_mdr <- jags.model(file = textConnection(model_code),
                                 data = list(tpos = sa_subset_abx$n_mdr_abx_r, 
                                             sample_size = sa_subset_abx$n_mdr, 
                                             se_alpha = abx_se_sp$se_alpha, 
                                             se_beta = abx_se_sp$se_beta, 
                                             sp_alpha = abx_se_sp$sp_alpha, 
                                             sp_beta = abx_se_sp$sp_beta), 
                                 n.chains = 3, quiet = TRUE)
    update(sa_subset_jags_mdr, 10000)
    
    mcmc_sa_subset_jags_mdr <- coda.samples(sa_subset_jags_mdr, variable.names=c("phi"), n.iter=20000)
    omcmc_sa_subset_jags_mdr <- summary(mcmc_sa_subset_jags_mdr)
    sa_input[which(sa_input$isolate_type == isolate_type & sa_input$drug == a), "mean_mdr"] <- omcmc_sa_subset_jags_mdr$statistics[1]
    sa_input[which(sa_input$isolate_type == isolate_type & sa_input$drug == a), "sd_mdr"] <- omcmc_sa_subset_jags_mdr$statistics[2]
    sa_input[which(sa_input$isolate_type == isolate_type & sa_input$drug == a), "lo_mdr"] <- omcmc_sa_subset_jags_mdr$quantiles[1]
    sa_input[which(sa_input$isolate_type == isolate_type & sa_input$drug == a), "hi_mdr"] <- omcmc_sa_subset_jags_mdr$quantiles[5]

    }
}

#Change bias correction for INH MONO-R

for(isolate_type in c("public", "survey", "total")) {
  sa_subset<- sa_input[which(sa_input$isolate_type == isolate_type),]
  print(c("now calculating", paste(sa_subset)))
  abx_se_sp <- se_sp["inh_mono",]
  sa_subset_abx <- sa_subset[which(sa_subset$drug == "inh"), ]

    #JAGS model for INH-R RS
  sa_subset_jags_rs <- jags.model(file = textConnection(model_code),
                                  data = list(tpos = sa_subset_abx$n_rs_abx_r, 
                                              sample_size = sa_subset_abx$n_rs, 
                                              se_alpha = abx_se_sp$se_alpha, 
                                              se_beta = abx_se_sp$se_beta, 
                                              sp_alpha = abx_se_sp$sp_alpha, 
                                              sp_beta = abx_se_sp$sp_beta), 
                                  n.chains = 3, quiet = TRUE)
  update(sa_subset_jags_rs, 10000)
  
  mcmc_rs <- coda.samples(sa_subset_jags_rs, variable.names=c("phi"), n.iter=20000)
  omcmc_rs <- summary(mcmc_rs)
  sa_input[which(sa_input$isolate_type == isolate_type & sa_input$drug == "inh"), "mean_rs"] <- omcmc_rs$statistics[1]
  sa_input[which(sa_input$isolate_type == isolate_type & sa_input$drug == "inh"), "sd_rs"] <- omcmc_rs$statistics[2]
  sa_input[which(sa_input$isolate_type == isolate_type & sa_input$drug == "inh"), "lo_rs"] <- omcmc_rs$quantiles[1]
  sa_input[which(sa_input$isolate_type == isolate_type & sa_input$drug == "inh"), "hi_rs"] <- omcmc_rs$quantiles[5]
  
}

write.csv(sa_input, "sa_jags_output_w_mono_inh_6.4.2021.csv", quote = FALSE, row.names = FALSE)

sa_test <- sa_input[which(sa_input$drug == "inh" & sa_input$isolate_type == "survey"),]

set.seed(123)
sa_inh_rs <- jags.model(file = textConnection(model_code),
           data = list(tpos = sa_test$n_abx_r, 
                       sample_size = sa_test$n_isolates, 
                       se_alpha = se_sp["inh", "se_alpha"], 
                       se_beta = se_sp["inh", "se_beta"], 
                       sp_alpha = se_sp["inh", "sp_alpha"], 
                       sp_beta = se_sp["inh", "sp_beta"]), 
           n.chains = 3, quiet = TRUE)

update(sa_inh_rs, 10000)

mcmc_sa_inh_rs <- coda.samples(sa_inh_rs, variable.names=c("phi"), n.iter=20000)
omcmc_sa_inh_rs <- summary(mcmc_sa_inh_rs)

omcmc_sa_inh_rs

jags_all_mono_inh <- jags.model(file = textConnection(model_code),
                      data = list(tpos = 1382, 
                                  sample_size = 14012, 
                                  se_alpha = se_sp["inh_mono","se_alpha"], 
                                  se_beta = se_sp["inh_mono","se_beta"], 
                                  sp_alpha = se_sp["inh_mono","sp_alpha"], 
                                  sp_beta = se_sp["inh_mono","sp_beta"], 
                      n.chains = 3, quiet = TRUE))

update(jags_all_mono_inh, 10000)

mcmc_jags_all_mono_inh <- coda.samples(jags_all_mono_inh, variable.names=c("phi"), n.iter=20000)
omcmc_jags_all_mono_inh <- summary(mcmc_jags_all_mono_inh)

jags_all_rr_inh <- jags.model(file = textConnection(model_code),
                                data = list(tpos = 4688, 
                                            sample_size = 5058, 
                                            se_alpha = se_sp["inh","se_alpha"], 
                                            se_beta = se_sp["inh","se_beta"], 
                                            sp_alpha = se_sp["inh","sp_alpha"], 
                                            sp_beta = se_sp["inh","sp_beta"], 
                                            n.chains = 3, quiet = TRUE))

update(jags_all_rr_inh, 10000)

mcmc_jags_all_rr_inh <- coda.samples(jags_all_rr_inh, variable.names=c("phi"), n.iter=20000)
omcmc_jags_all_rr_inh <- summary(mcmc_jags_all_rr_inh)

#WHO corrected global INH-R
(omcmc_jags_all_rr_inh$statistics["Mean"]*0.007031705) + (omcmc_jags_all_mono_inh$statistics["Mean"]*(1-0.007031705))

(0.007031705^2)*(omcmc_jags_all_rr_inh$statistics["SD"]^2) + (omcmc_jags_all_rr_inh$statistics["Mean"]^2)*0.0002384564


jags_all_mono_levo <- jags.model(file = textConnection(model_code),
                                data = list(tpos = 135, 
                                            sample_size = 14012, 
                                            se_alpha = se_sp["levo","se_alpha"], 
                                            se_beta = se_sp["levo","se_beta"], 
                                            sp_alpha = se_sp["levo","sp_alpha"], 
                                            sp_beta = se_sp["levo","sp_beta"], 
                                            n.chains = 3, quiet = TRUE))

update(jags_all_mono_levo, 10000)

mcmc_jags_all_mono_levo <- coda.samples(jags_all_mono_levo, variable.names=c("phi"), n.iter=20000)
omcmc_jags_all_mono_levo <- summary(mcmc_jags_all_mono_levo)

jags_mdr_pza <- jags.model(file = textConnection(model_code),
                                 data = list(tpos = 1551, 
                                             sample_size = 3964, 
                                             se_alpha = se_sp["pza","se_alpha"], 
                                             se_beta = se_sp["pza","se_beta"], 
                                             sp_alpha = se_sp["pza","sp_alpha"], 
                                             sp_beta = se_sp["pza","sp_beta"], 
                                             n.chains = 3, quiet = TRUE))

update(jags_mdr_pza, 10000)

mcmc_jags_mdr_pza <- coda.samples(jags_mdr_pza, variable.names=c("phi"), n.iter=20000)
omcmc_jags_mdr_pza <- summary(mcmc_jags_mdr_pza)


jags_mdr_levo <- jags.model(file = textConnection(model_code),
                                 data = list(tpos = 1042, 
                                             sample_size = 3964, 
                                             se_alpha = se_sp["levo","se_alpha"], 
                                             se_beta = se_sp["levo","se_beta"], 
                                             sp_alpha = se_sp["levo","sp_alpha"], 
                                             sp_beta = se_sp["levo","sp_beta"], 
                                             n.chains = 3, quiet = TRUE))

update(jags_mdr_levo, 10000)

mcmc_jags_mdr_levo <- coda.samples(jags_mdr_levo, variable.names=c("phi"), n.iter=20000)
omcmc_jags_mdr_levo <- summary(mcmc_jags_mdr_levo)


se_sp <- read.csv("../bias_correction/se_sp_alpha_beta_w_inh_mono_shrt_crs.csv", row.names = 1, header = TRUE)
se_sp["shrt_crs",] <- transform(se_sp["shrt_crs",], sens_sd = (TPR_hiCI-TPR_lowCI)/3.92, spec_sd = (TNR_hiCI-TNR_lowCI)/3.92)
se_sp["shrt_crs",] <- transform(se_sp["shrt_crs",], sens_var = sens_sd*sens_sd, spec_var = spec_sd*spec_sd)
se_sp["shrt_crs",] <- transform(se_sp["shrt_crs",], 
                   se_alpha = apply(se_sp["shrt_crs",], 1, function(x) alpha(mean = x["TPR"], var = x["sens_var"])),
                   se_beta = apply(se_sp["shrt_crs",], 1, function(x) beta(mean = x["TPR"], var = x["sens_var"])),
                   sp_alpha = apply(se_sp["shrt_crs",], 1, function(x) alpha(mean = x["TNR"], var = x["spec_var"])),
                   sp_beta = apply(se_sp["shrt_crs",], 1, function(x) beta(mean = x["TNR"], var = x["spec_var"])))
write.csv(se_sp, "se_sp_alpha_beta_w_inh_mono_shrt_crs.csv", quote = FALSE)

geno_counts <- read.csv("../geno_counts_6.29.21.csv", row.names = 1, header = TRUE)
shrt_crs_corected <- geno_counts[, c("n_tested_mdr_shrt_crs", "n_mdr_shrt_crs_r_2_drug")]
shrt_crs_corected$n_shrt_crs_elig <- shrt_crs_corected$n_tested_mdr_shrt_crs - shrt_crs_corected$n_mdr_shrt_crs_r_2_drug
shrt_crs_corected$mean_shrt_crs_elig <- NA
shrt_crs_corected$sd_shrt_crs_elig <- NA
shrt_crs_corected$lo_shrt_crs_elig <- NA
shrt_crs_corected$hi_shrt_crs_elig <- NA


for (i in 1:nrow(shrt_crs_corected)) {
  set.seed(123)
  jags_shrt_crs <- jags.model(file = textConnection(model_code),
                              data = list(tpos = shrt_crs_corected[i,"n_shrt_crs_elig"],
                                          sample_size = shrt_crs_corected[i,"n_tested_mdr_shrt_crs"], 
                                          se_alpha = se_sp["shrt_crs", "se_alpha"], 
                                          se_beta = se_sp["shrt_crs","se_beta"], 
                                          sp_alpha = se_sp["shrt_crs","sp_alpha"], 
                                          sp_beta = se_sp["shrt_crs","sp_beta"], 
                                          n.chains = 3, quiet = TRUE))
  update(jags_shrt_crs, 10000)
  
  mcmc_jags_shrt_crs <- coda.samples(jags_shrt_crs, variable.names=c("phi"), n.iter=20000)
  omcmc_jags_shrt_crs <- summary(mcmc_jags_shrt_crs)
  shrt_crs_corected[i,"mean_shrt_crs_elig"] <- omcmc_jags_shrt_crs$statistics[1]
  shrt_crs_corected[i,"sd_shrt_crs_elig"] <- omcmc_jags_shrt_crs$statistics[2]
  shrt_crs_corected[i,"lo_shrt_crs_elig"] <- omcmc_jags_shrt_crs$quantiles[1]
  shrt_crs_corected[i,"hi_shrt_crs_elig"] <- omcmc_jags_shrt_crs$quantiles[5]
}

shrt_crs_corected <- shrt_crs_corected[which(shrt_crs_corected$n_tested_mdr_shrt_crs > 99),]

write.csv(shrt_crs_corected, "shrt_crs_bias_corrected.csv")

jags_global_shrt_crs <- jags.model(file = textConnection(model_code),
                                   data = list(tpos = sum(shrt_crs_corected[which(shrt_crs_corected$n_tested_mdr_shrt_crs > 99),"n_shrt_crs_elig"]),
                                               sample_size = sum(shrt_crs_corected[which(shrt_crs_corected$n_tested_mdr_shrt_crs > 99),"n_tested_mdr_shrt_crs"]), 
                                               se_alpha = se_sp["shrt_crs", "se_alpha"], 
                                               se_beta = se_sp["shrt_crs","se_beta"], 
                                               sp_alpha = se_sp["shrt_crs","sp_alpha"], 
                                               sp_beta = se_sp["shrt_crs","sp_beta"], 
                                               n.chains = 3, quiet = TRUE))
update(jags_global_shrt_crs, 10000)

mcmc_jags_global_shrt_crs <- coda.samples(jags_global_shrt_crs, variable.names=c("phi"), n.iter=20000)
omcmc_jags_global_shrt_crs <- summary(mcmc_jags_global_shrt_crs)

omcmc_jags_global_shrt_crs$statistics
omcmc_jags_global_shrt_crs$quantiles
