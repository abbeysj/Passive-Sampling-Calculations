# INVLAP  numerical inverse Laplace transform
#
# f = invlap(F, t, alpha, tol);
#         
# F       laplace-space function 
#           This is a handle to a function which takes laplace space vector p to F(p) 
# t       column vector of times for which real-space function values are
#           sought
# alpha   largest pole of F (default zero)
# tol     numerical tolerance of approaching pole (default 1e-9)
# f       vector of real-space values f(t)
#
# example: identity function in Laplace space:
#   function F = identity(s);                    # save these two lines
#            F = 1./(s.^2);                      # ...  as "identity.m"
#   invlap('identity', [1;2;3])                  # gives [1;2;3]
#
# algorithm: de Hoog et al's quotient difference method with accelerated 
#   convergence for the continued fraction expansion
#   [de Hoog, F. R., Knight, J. H., and Stokes, A. N. (1982). An improved 
#    method for numerical inversion of Laplace transforms. S.I.A.M. J. Sci. 
#    and Stat. Comput., 3, 357-366.]
# Modification: The time vector is split in segments of equal magnitude
#   which are inverted individually. This gives a better overall accuracy.   

#  details: de Hoog et al's algorithm f4 with modifications (T->2*T and 
#    introduction of tol). Corrected error in formulation of z.
#
#  Adapted from invlap.m by Kevin Joyce (kevin.t.joyce@gmail.com) on 28 April 2016
#
#  original copyright: Karl Hollenbeck
#             Department of Hydrodynamics and Water Resources
#             Technical University of Denmark, DK-2800 Lyngby
#             email: karl@isv16.isva.dtu.dk
#  22 Nov 1996, MATLAB 5 version 27 Jun 1997 updated 1 Oct 1998
#  IF YOU PUBLISH WORK BENEFITING FROM THIS M-FILE, PLEASE CITE IT AS:
#    Hollenbeck, K. J. (1998) INVLAP.M: A matlab function for numerical 
#    inversion of Laplace transforms by the de Hoog algorithm, 
#    http://www.isva.dtu.dk/staff/karl/invlap.htm 

invlap = function(F, t, alpha = 0, tol=1e-9){

  f = c();
  # split up t vector in pieces of same order of magnitude, invert one piece
  #   at a time. simultaneous inversion for times covering several orders of 
  #   magnitudes gives inaccurate results for the small times.
  allt = as.numeric(t) 				# save full times vector
  logallt = log10(allt)
  iminlogallt = floor(min(logallt))
  imaxlogallt = ceiling(max(logallt))
  for(ilogt in iminlogallt:imaxlogallt ) {	# loop through all pieces
    
    t = allt[which( (logallt>=ilogt) & (logallt<(ilogt+1)) )]
    if( length(t) > 0 ){                        # maybe no elements in that magnitude

      T = max(t)*2
      gamma = alpha-log(tol)/(2*T)
      # NOTE: The correction alpha -> alpha-log(tol)/(2*T) is not in de Hoog's
      #   paper, but in Mathematica's Mathsource (NLapInv.m) implementation of 
      #   inverse transforms
      nt = length(t)
      M = 20
      run = 0:(2*M);    # so there are 2M+1 terms in Fourier series expansion

      # find F argument, call F with it, get 'a' coefficients in power series
      s = gamma + 1i*pi*run/T
      #command = ['a = ' F '(s']
      #if nargin > 4,  			# pass on parameters
      #  for iarg = 1:nargin-4,
      #        command = [command ',P' int2str(iarg)]
      #  end
      #end
      #command = [command ');']
      #eval(command)
      a = F(s)
      a[1] = a[1]/2      			# zero term is halved

      # build up e and q tables. superscript is now row index, subscript column
      #   CAREFUL: paper uses null index, so all indeces are shifted by 1 here
      e = matrix(0,2*M+1, M+1)
      q = matrix(0,2*M  , M+1)   		# column 0 (here: 1) does not exist
      q[,2] = a[-1]/a[-length(a)]
      for(r in 2:(M+1)) {           		# step through columns (called r...)
        e[1:(2*(M-r+1)+1),r] = q[2:(2*(M-r+1)+2),r] - q[1:(2*(M-r+1)+1),r] + e[2:(2*(M-r+1)+2),r-1]
        if( r < M+1 ) {               		# one column fewer for q
          rq = r+1
          q[1:(2*(M-rq+1)+2),rq] = q[2:(2*(M-rq+1)+3),rq-1]*e[2:(2*(M-rq+1)+3),rq-1]/e[1:(2*(M-rq+1)+2),rq-1]
        }
      }

      # build up d vector (index shift: 1)
      d = 0*(1:(2*M+1))
      d[1] = a[1]
      d[2*(1:M)] = -q[1,2:(M+1)] # these 2 lines changed after niclas
      d[2*(1:M)+1] = -e[1,2:(M+1)] # ...

      # build up A and B vectors (index shift: 2) 
      #   - now make into matrices, one row for each time
      A = matrix(1,2*M+2,nt)
      B = matrix(1,2*M+2,nt)
      A[1,] = 0*A[1,];
      A[2,] = d[1]*A[2,];
      z = exp(1i*pi*t/T);
      # after niclas back to the paper (not: z = exp(-i*pi*t/T)) !!!
      for (n in 3:(2*M+2)) {
        A[n,] = A[n-1,] + d[n-1]*z*A[n-2,];  # different index 
        B[n,] = B[n-1,] + d[n-1]*z*B[n-2,];  #  shift for d!
      } 

      # double acceleration
      h2M = .5 * ( 1 + ( d[2*M]-d[2*M+1] )*z )
      R2Mz = -h2M*(1 - sqrt(1+d[2*M+1]*z/h2M^2))
      A[2*M+2,] = A[2*M+1,] + R2Mz * A[2*M,]
      B[2*M+2,] = B[2*M+1,] + R2Mz * B[2*M,]

      # inversion, vectorized for times, make result a column vector
      fpiece = 1/T * exp(gamma*t) * Re(A[2*M+2,]/B[2*M+2,]) 
      f = c(f, fpiece)			# put pieces together
    } # if not empty time piece
  } # loop through time vector pieces
  return(f)
}
