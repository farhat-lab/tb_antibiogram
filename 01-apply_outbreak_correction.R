##Load required package
#install.packages("igraph")
library(igraph)

##load RData file that contains the pairwise distance matrix, its annotations, list of isolates in the antibiogram dataset, 
load("pairwise_dist_geo.RData")

##Load the gentb model output in binary format
model_out <- read.csv("gentb_binary_4.9.2021.txt", stringsAsFactors = FALSE, row.names = 1)

##Ensure names match
malawi_names <- read.table("Malawi_BioSample_RunIDs.txt", stringsAsFactors = FALSE)
isolate_annot$isolate_ID <- unlist(lapply(isolate_annot$isolate_ID, function(x) {if(x %in% malawi_names$V5){x <- malawi_names[which(malawi_names$V5 == x), "V1"]} else {x <- x}}))
rownames(model_out) <- unlist(lapply(rownames(model_out), function(x) {if(x %in% malawi_names$V5){x <- malawi_names[which(malawi_names$V5 == x), "V1"]} else {x <- x}}))

##Get geometadata
geog <- read.csv("geog_cleaned_4.6.21.txt", sep = "\t", stringsAsFactors = FALSE)

##This function uses the igraph library to create a network of isolates then decompose them into separate networks if not interconnected by country. From each decomposed network, the function checks if all isolates are connected to each other... If yes, it randomly selects one isolate from the group. If no, it excludes the isolate with highest degree of connectedness and iterates previous steps until all connected isolates have been checked. 

subset_isolates_by_connectedness <- function(country,snp_threshold) {
    country_isolates <- na.omit(geog_complete[which(geog_complete$isolation_country == paste(country)),]$BioSample)
    country_isolates <- country_isolates[!is.na(country_isolates)]#subset to isolates belonging to country specified by user input 
    if (length(country_isolates) > 1) {
       pairwise_dist_country <- pairwise_dist[which(isolate_annot$isolate_ID %in% country_isolates),which(isolate_annot$isolate_ID %in% country_isolates)]
        country_isolate_annot <- isolate_annot$isolate_ID[as.numeric(row.names(pairwise_dist_country))]
        if (length(country_isolate_annot) > 0) {
            #identify indices of cells that are <= snp_threshold
        close_isolates <- which(pairwise_dist_country <= snp_threshold, arr.ind = TRUE)
        #remove isolate pairs that are duplicated 
        close_isolates <- close_isolates[apply(close_isolates, MARGIN =  1, FUN = function(x) !any(duplicated(x))), ]
        if(nrow(close_isolates) > 0) {
            close_pairs <- data.frame(t(apply(close_isolates,1,sort))) #sort the data frame before removing duplicates
            close_pairs <- close_pairs[!duplicated(close_pairs),]
            #Get the snp difference and add it to the data.frame
            snp_diff <- pairwise_dist_country[as.matrix(close_pairs)]
            close_pairs_values <- cbind(close_pairs, snp_diff)
            g <- graph_from_data_frame(close_pairs_values, directed = FALSE)
            #prune isolates by connectedness and randomly select from groups with 100% connectedness
            g1 <- g
            checked_isolates <- vector("character")
            selected_isolates <- vector("character")
            high_degree_isolates <- vector("character")
            max.degree <- 100
            while (length(max.degree) > 0) {
                dg <- decompose.graph(g1)
                for (i in 1:length(dg)) {
                    if(degree_distribution(dg[[i]])[length(degree_distribution(dg[[i]]))] == 1) {
                        checked_isolates <- append(checked_isolates, V(dg[[i]])$name)
                        set.seed(123)
                        dg1 <- delete_vertices(dg[[i]], sample(V(dg[[i]]), length(V(dg[[i]]))-1))
                        selected_isolates <- append(checked_isolates, V(dg1)$name)
                    } else { 
                        high_degree_isolates <- append(high_degree_isolates, names(which(degree(dg[[i]]) == degree(dg[[i]])[which.max(degree(dg[[i]]))])))
                        dg1 <- delete_vertices(dg[[i]], which(degree(dg[[i]]) == degree(dg[[i]])[which.max(degree(dg[[i]]))]))
        }
    }
                g1 <- delete_vertices(g1, V(g1)[which(V(g1)$name %in% c(checked_isolates, high_degree_isolates))])
                max.degree <- degree(g1)[which.max(degree(g1))]
            }

            #obtain isolate names
            isolate_groups_names <- country_isolate_annot[as.numeric(V(g)$name)]
            selected_isolates_names <- country_isolate_annot[as.numeric(selected_isolates)]
            #separate ungrouped and grouped isolates
            ungrouped_isolates <- country_isolate_annot[!country_isolate_annot %in% unlist(isolate_groups_names)]
            grouped_isolates <- country_isolate_annot[country_isolate_annot %in% unlist(isolate_groups_names)]
            grouped_isolate_subset <- country_isolate_annot[country_isolate_annot %in% unlist(selected_isolates_names)]

            #add subset from grouped isolates to ungrouped isolates to get final subset
            isolate_subset <- c(ungrouped_isolates, grouped_isolate_subset)
            } else {
            isolate_subset <- country_isolate_annot
            } 
        } else {
            isolate_subset <- country_isolate_annot
        }
        } else {
            isolate_subset <- country_isolates
        }
    return(isolate_subset)                                         
}

##Get country names
countries <- sort(unique(geog_complete$isolation_country))


#This function calculates number of total isolates available by country in the gentb model output
n_isolates <- function(country) {
  country_isolates <- geog_complete[which(geog_complete$isolation_country == paste(country)),]
  country_resistant <- model_out[which(rownames(model_out) %in% na.omit(country_isolates$BioSample)),]
  return(nrow(country_resistant))
}

##call the n_isolates function on all countries
all_country_count <- vector("list")
for (i in countries){
 all_country_count[[i]] <-  n_isolates(i)
}

##Call the above function on all countries and only keep isolates that are included in the pairwise distance matrix
set.seed(123)
all_country_subset <- vector("list")
for (i in countries){
 all_country_subset[[i]] <-  subset_isolates_by_connectedness(i, 10)
}

all_country_subset <- lapply(all_country_subset, function(x) x[x %in% isolate_annot$isolate_ID])

country_samples <- as.matrix(unlist(all_country_count))

country_samples_subset <- as.matrix(sapply(all_country_subset, length))

##subset the gentb model output to only selected isolates
model_subset <- model_out[unlist(all_country_subset), ]

##write the subsetted gentb model output
write.csv(model_subset, "gentb_subset_snp10_binary_4.9.21.txt", quote = FALSE)
