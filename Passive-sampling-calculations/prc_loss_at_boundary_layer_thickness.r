rm(list=ls()) # Reset all the workspace
source( 'prc_analysis_functions.r' ) #load functions - e.g. mass_in, mass_out, etc.
source( 'invlap.r' ) 
prc_data = non_dimensionalize(read.csv('12umPRC.csv')) #reads prc data from saved excel wksheet and add non-dimensional variables
target_data = non_dimensionalize(read.csv("target_data.csv")) #see above
prc_data$boundary_layer = apply( prc_data[,-c(1,2)],1,solve_boundary_layer ) # exclude the string rows that make the matrix coercian complicate things

good_prc_idx = prc_data$M < .85 & prc_data$M > .15 #define PRCs with acceptable losses
good_prc_data = prc_data[good_prc_idx,]
good_prc_data$station = as.numeric(gsub("[^0-9]","",good_prc_data$site)) # Finds the station number from the site string, e.g. NBH-2a -> 2



station_means = tapply(good_prc_data$boundary_layer,good_prc_data$station,mean)
good_prc_data$station_means = station_means[factor(good_prc_data$station)]
site_means = tapply(good_prc_data$boundary_layer,good_prc_data$site,mean)

get_fraction_in = function(data_row) {
  return(apply(site_means,1,function(d){ solve_diffusion_in(data_row,as.numeric(d)) } ))
}  #runs mass_in for each column of sites defined in prc analysis

feq_matrix = apply(target_data,1,get_fraction_in) #runs mass_in for each row or target compound
target_data = cbind(target_data,t(feq_matrix))

write.csv(prc_data,"12umPRC_with_boundary_layer.csv")
write.csv(good_prc_data,"12umPRC_GOOD_with_boundary_layer.csv")
write.csv(target_data,"target_data_fraction_to_eq.csv")
