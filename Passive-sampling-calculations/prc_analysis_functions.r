# This non-dimensionalizes data as described in Thommpson et.al.
non_dimensionalize = function(data) {
  data$T = data$t..sec. * data$Dpe / data$l..cm.^2
  data$Y = data$Dw/data$Dpe
  return(data)
}

# Mass out analytic Laplace function
mass_out_function = function(s,kpew,d,Y){
  1/s - ( s^(3/2) * (tanh(sqrt(s))^(-1) + kpew*tanh(d*sqrt(s/Y)) / sqrt(Y)) )^(-1)
}

# Mass in analytic Laplace function
mass_in_function = function(s,kpew,d,Y){
  ( s^(3/2) * (tanh(sqrt(s))^(-1) + kpew*tanh(d*sqrt(s/Y)) / sqrt(Y)) )^(-1)
}

# This function solves for the fraction to diffused  IN after a given 
# time (in the prc_data_row) and boundary layer (BL)
solve_diffusion_in = function(prc_data_row,BL) {
  d = as.numeric(BL/prc_data_row['l..cm.'])
  kpew = as.numeric(prc_data_row['Kpew'])
  Y = as.numeric(prc_data_row['Y'])
  laplace_solution = function(s){ mass_in_function(s,kpew,d,Y) }
  return( invlap(laplace_solution, prc_data_row['T']) )
}

# This function solves for the fraction to diffused OUT after a given 
# time (in the prc_data_row) and boundary layer (BL)
solve_diffusion_out = function(prc_data_row,BL) {
  d = as.numeric(BL/prc_data_row['l..cm.'])
  kpew = as.numeric(prc_data_row['Kpew'])
  Y = as.numeric(prc_data_row['Y'])
  Time = as.numeric(prc_data_row['T'])
  laplace_solution = function(s){ mass_out_function(s,kpew,d,Y) }
  return( invlap(laplace_solution, Time) )
}

# Inverse solve for delta in Thompson, by inverting diffusion curve 
# in delta under the assumption that concentration on the PRC increases 
# as delta increases.
solve_boundary_layer = function(prc_data_row, bl_guess = .02, bl_max = .1, tol=1e-9, maxIter=100){
  l = 0        # lower bound (non-dimensional)
  u = bl_max   # upper bound (non-dimensional)
  d = bl_guess # initial guess (non-dimensional)
  # Bisection method for inverting C(delta) 
  # We use the assumption that d/d\delta C(\delta) < 0
  for( i in 1:maxIter) {
    s = solve_diffusion_out(prc_data_row,d)
    if( abs(prc_data_row['M'] - s) < tol ) {
      return(d)
    }
    else if( prc_data_row['M'] - s <= 0 ){
      u = d
      d = (d+l)/2
    }
    else if( prc_data_row['M'] - s >= 0 ){
      l = d
      d = (d+u)/2
    }
  }
  warning("Max Iterations reached")
  return(d*prc_data_row['l..cm.'])
}
