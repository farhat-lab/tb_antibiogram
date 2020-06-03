
#R function that takes a country name and generates antibiogram corrected for MDR oversampling. Needs output of WDNN model with drugs as column names and converted into binary based on predetermined thersholds based on Chen et al. 2019,
# list of isolates corresponding to each row of the WDNN model output, geographic metadata with first column labelled 'BioSample' and second column labelled 'isolation_country'
# Also need to provide the prevalence of new MDR/RR for the country from WHO website in the form of probability 

create_antibiogram <- function(country, p_country_mdr, model_out = "model_output_binary.csv", isolates = "isolate_names.txt", geo_metadata = "geo_sampling.txt") {
  #read all required files
  model_output <- read.csv(paste(model_out))
  isolate_list <- readLines(paste(isolates))
  rownames(model_output) <- isolate_list
  geog <- read.csv(paste(geo_metadata), sep = "\t", stringsAsFactors = FALSE)
  
  country_isolates <- geog[which(geog$isolation_country == country),]
  country_resistant <- model_output[which(rownames(model_output) %in% country_isolates$BioSample),]
  
  country_resistant_count <- vector()
  
  for (i in 1:ncol(country_resistant)) {
    country_resistant_count[i]  <-  sum(country_resistant[,i])
  }
  
  names(country_resistant_count) <- colnames(country_resistant)
  
  #correct for mdr oversampling
  
  p_a_mdr <- function(a) {
    sum(country_resistant[which(country_resistant$rif == 1),paste(a)])/sum(country_resistant$rif)
  }
  
  p_mdr_a <- function(a) {
    sum(country_resistant[which(country_resistant[,paste(a)] ==1), "rif"])/sum(country_resistant[,paste(a)])
  }
  
  p_a_res <- function(a) {
    (p_a_mdr(a)*p_country_mdr)/p_mdr_a(a)
  }
  
  sum_r <- apply(country_resistant, 1, sum)
  mono_r <- which(sum_r == 1)
  a_mono_r <- function(a) {
    sum(as.numeric(row.names(country_resistant[which(country_resistant[,paste(a)] == 1),]) %in% names(mono_r)))
  }
  corrected_count <- vector()
  
  for (i in 1:ncol(country_resistant)) {
    if (colnames(country_resistant)[i] == "rif") { ##correction for mono_resistance does not apply for rifampicin
      corrected_count[i] <- p_a_res(colnames(country_resistant)[i])
    } else {
      corrected_count[i] <- p_a_res(colnames(country_resistant)[i]) + a_mono_r(colnames(country_resistant)[i])/nrow(country_resistant) 
    }
  }
  names(corrected_count) <- colnames(country_resistant)
  antibiogram_corrected <- barplot(corrected_count*100, ylim = c(0.0,max(corrected_count*100)+5), main = paste("Antibiogram for", country, "(corrected)", paste("(n = ", nrow(country_resistant), ")", sep = ""), sep = " "), xlab = "Antibiotic", ylab = "Percent of resistant isolates", col = "orange")
  text(antibiogram_corrected, y = 0, round(corrected_count*100, 1), cex = 1, pos = 3, col = "black")
}

pdf("Antibiograms.pdf")
create_antibiogram(country = "Peru", p_country_mdr = 0.063)
create_antibiogram(country = "India", p_country_mdr = 0.028)
create_antibiogram(country = "South Africa", p_country_mdr = 0.034)
create_antibiogram(country = "Russia", p_country_mdr = 0.32)
create_antibiogram(country = "Malawi", p_country_mdr = 0.0075)
create_antibiogram(country = "United Kingdom", p_country_mdr = 0.014)
create_antibiogram(country = "Belarus", p_country_mdr = 0.38)
create_antibiogram(country = "Uzbekistan", p_country_mdr = 0.15)
create_antibiogram(country = "Germany", p_country_mdr = 0.022)
create_antibiogram(country = "China", p_country_mdr = 0.071)
create_antibiogram(country = "Netherlands", p_country_mdr = 0.0092)
create_antibiogram(country = "United States", p_country_mdr = 0.015)
create_antibiogram(country = "Vietnam", p_country_mdr = 0.036)
dev.off()
