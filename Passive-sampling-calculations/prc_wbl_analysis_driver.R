rm(list=ls()) # Reset all the workspace

prc_input_filename = '76umPRC.csv' #call data set you want to analyze
target_filename    = 'target_data.csv' #double check sheet thickness
(prc_BL_filename   = paste0('boundary_layer_',prc_input_filename))
(prc_good_BL_filename = paste0('good_boundary_layer_',prc_input_filename))
(target_feq_filename = paste0('target_feq_',prc_input_filename))

source( 'prc_analysis_functions.r' ) #load functions - e.g. mass_in, mass_out, etc.
source( 'invlap.r' ) 

prc_data = non_dimensionalize(read.csv(prc_input_filename)) #reads prc data from saved excel wksheet and add non-dimensional variables. 
target_data = non_dimensionalize(read.csv(target_filename)) #see above

#solve for boundary layer using M from PRCs using bisection method. starts guessing at .02cm will range from 0 - 0.1 cm
# UNCOMMENT TO change defaults
prc_data$boundary_layer = apply( prc_data[,-c(1,2)],1,function(row){ solve_boundary_layer(row,bl_guess=.01,bl_max=.2,tol=1e-8,maxIter=200)}) 
prc_data$boundary_layer = apply( prc_data[,-c(1,2)],1,solve_boundary_layer ) # exclude the string rows that make the matrix coercian complicate things

good_prc_idx = prc_data$M < .85 & prc_data$M > .15 #define PRCs with acceptable losses
good_prc_data = prc_data[good_prc_idx,]
good_prc_data$station = as.numeric(gsub("[^0-9]","",good_prc_data$site)) # Finds the station number from the site string, e.g. NBH-2a -> 2



station_means = tapply(good_prc_data$boundary_layer,good_prc_data$station,mean) #takes mean of "good" boundary layer by station
good_prc_data$station_means = station_means[factor(good_prc_data$station)] #adds averages to good_prc_data.csv

site_means = tapply(good_prc_data$boundary_layer,good_prc_data$site,mean)
good_prc_data$site_means = site_means[factor(good_prc_data$site)] #adds averages to good_prc_data.csv by site

get_fraction_in = function(data_row) {
  return(apply(site_means,1,function(d){ solve_diffusion_in(data_row,as.numeric(d)) } ))
}  #runs mass_in for each column of sites defined in prc analysis

feq_matrix = apply(target_data,1,get_fraction_in) #runs mass_in for each row or target compound
target_data = cbind(target_data,t(feq_matrix)) #combines target data w/ feq results

## UNCOMMENT TO SAVE WORK
write.csv(prc_data,prc_BL_filename)
write.csv(good_prc_data,prc_good_BL_filename)
write.csv(target_data,target_feq_filename)
