

g_series <- function(t, x, N_max=300) {
      if (any(t <= 0) || x <= 0) stop("Require t > 0 and x > 0")
      
      # --- 1. choose truncation point N so last retained term < tol ----
      #a_min  <- pi^2 * min(t) / (2 * x^2)          # slowest decay rate
      #N_max  <- ceiling( sqrt( max(1, -log(tol) / a_min) ) ) + 2  # safety margin
      
      n      <- seq_len(N_max)                      # 1,2,...,N_max
      coef   <- (-1)^(n + 1) * n^2                  # alternating n^2
      
      # --- 2. vectorised summation over n for all t in one shot --------
      a_vec  <- pi^2 * t / (2 * x^2)                # length = length(t)
      # outer returns an N_max × length(t) matrix; multiply by coef & sum rows
      g_vals <- drop( coef %*% exp( -outer(n^2, a_vec) ) )
      
      norm_const <- (pi^2)/x^2
      
      
      norm_const*g_vals
}


#-------------------------------------
#### define truncated densities
#### these all relate to the g function

dinv_gamma_trunc_g <-function(t,tstar,x){
      const <- x^3 / sqrt(2 * pi)              
      out <- numeric(length(t))                     # initialise to 0
      #idx <- (t < tstar) & (t > 0)                      # logical mask
      idx <- (t <= tstar) & (t > 0)
      N1=F_b(tstar, 3/2,x^2/2) # normalizing constant for the IG density
      
      out[idx] <- (const / N1) * t[idx]^(-5/2) * exp(-x^2 / (2 * t[idx]))
      out
}


dexp_trunc_g <- function(t,tstar,lam){
      out <- numeric(length(t))                      # initialise to 0
      idx <- t>= tstar
      N2 = exp(-lam*tstar) # normalizing constant for the Exp density
      
      out[idx] <- (lam / N2) * exp(-lam * t[idx])    # vectorised fill
      out
}

# normalizing constant for the inverse gamma
F_b <- function(b, alpha, beta){
      1 - pgamma(beta / b, shape = alpha, rate = 1, lower.tail = TRUE)
}

#--------------------------------------


#----------------------------------------------------
# random sampler from the truncated inverse gamma
rinv_gamma_trunc <- function(n, alpha, beta, b)
{
      if (n <= 0 || alpha <= 0 || beta <= 0 || b <= 0)
            stop("n, α, β, b must all be positive.")
      
      # tail probability to be kept
      p0 <- pgamma(beta / b, shape = alpha, rate = 1, lower.tail = TRUE)
      
      if (p0 >= 1)                       # truncation cuts off all mass
            stop("b is too small: P(T < b) = 0 for these parameters.")
      
      # 1. uniform on the retained tail
      u <- runif(n, min = p0, max = 1)
      
      # 2. right‑truncated Gamma  ⇒  qgamma of those uniforms
      g <- qgamma(u, shape = alpha, rate = 1, lower.tail = TRUE)
      
      # 3. back‑transform
      beta / g
}

# random sampler from the truncated exponential
rexp_trunc <- function(n, rate, a, b)
{
      if (n <= 0 || rate <= 0 || a < 0 || a >= b)
            stop("Need n>0, λ>0, 0≤a<b.")
      
      Fa <- 1 - exp(-rate * a)          # F(a)
      Fb <- 1 - exp(-rate * b)          # F(b)
      
      u  <- runif(n, min = 0, max = 1)  # U ~ Unif(0,1)
      v  <- Fa + u * (Fb - Fa)          # F(a) + U·[F(b)-F(a)]
      
      -log(1 - v) / rate                # inverse CDF
}
#------------------------------------------------------------



####################################################
# simulate one draw from the g density
####################################################

upper_bound_g <- function(t,x,N){
     
     m_tilde = ceiling(1/sqrt(pi^2 * t / 2 / x^2))
     N_max = max(2*m_tilde-1, 2*N-1)
     out = g_series(t,x,N_max)
     return(out)
}

lower_bound_g <- function(t,x,N){
     
     m_tilde = ceiling(1/sqrt(pi^2 * t / 2 / x^2))
     N_max = max(2*m_tilde, 2*N)
     out = g_series(t,x,N_max)
     return(out)
}

sample_from_g <- function(x){
      
      
      lam = pi^2/(2*x^2)
      
      c_g <- log(3) * log(4) / (3 * pi^2)
      tstar <-  c_g * 2*x^2 / log(3)
      gamma_star <- 9 * exp(-4 * x^2 / tstar)
      
      
      D1 <- 2  /  (1 - gamma_star)
      N1=F_b(tstar, 3/2,x^2/2) # normalizing constant for the IG density
      tildeD1 = D1 * N1 
      

      D2 = 2
      N2 = exp(-lam*tstar) # normalizing constant for the Exp density
      tildeD2 = D2*N2
      
      
      a = tildeD1/(tildeD1+tildeD2)
      
      
      # accept reject step
      sw=0
      while(!sw){
            
            # propose one draw
            coin = runif(1,0,1)
            if(coin < a ){
                  #simulate from inv-gamma
                  Y=rinv_gamma_trunc(1,3/2, x^2/2, tstar)
            }
            else{
                  #simulate from truncated exp
                  Y = rexp_trunc(1,lam,tstar,Inf)
            }
            
            dq_g = a*dinv_gamma_trunc_g(Y,tstar,x) + (1-a)*dexp_trunc_g(Y,tstar,lam)
            N=200
            sw1=0
            U=runif(1,0,1)
            while(!sw1){
                    N = N+50
                    UB = upper_bound_g(Y,x,N)
                    LB = lower_bound_g(Y,x,N)
                    if( U * ( (tildeD1 + tildeD2)*dq_g ) < LB ){
                         sw=1 # draw is accepted
                         out=Y
                         return(out)
                    }
                    if( U * ( (tildeD1 + tildeD2)*dq_g ) > UB ){
                         sw1=1 # decision done, reject
                    }
          } # end while sw1
      
      } # end while sw
}

sample_from_g_T <- function(x,T){
     if (any(x <= 0) || T <= 0) stop("Require x > 0 and T > 0")
     repeat{
          temp = sample_from_g(x)
          if(temp<T){
               out=temp
               break
          }          
     }
     return(out)
}
####################################################




#################################################
#################################################
#################################################
#################################################


f_series <- function(t_vec,x,r,N_max=300) {
     # t_vec: numeric vector of times at which to evaluate f_{T_a}
     # r, a: scalars with 0 < r < a
     # N: number of terms in the series to sum (the series converges rapidly)
     #
     # Returns: numeric vector of length length(t_vec) giving f_{T_a}(t_vec)
     
     
     # --- 1. choose truncation point N so last retained term < tol ----
     # a_min  <- pi^2 * min(t_vec) / (2 * x^2)          # slowest decay rate
     # N_max  <- ceiling( sqrt( max(1, -log(tol) / a_min) ) ) + 2  # safety margin
     
     # Ensure t_vec is a numeric vector
     t_vec <- as.numeric(t_vec)
     M <- length(t_vec)
     
     # Precompute n = 1,2,...,N as a column
     n <- seq_len(N_max)                    # 1,2,...,N
     
     n_pi_r <- n * pi * r / x                # vector of (n*pi*r/x)
     
     
     
     # The overall front constant - normalized
     norm_const <- pi / (x * r)             # scalar
     
     
     # Compute sin(n*pi*r/a) once:
     sin_npr   <- sin(n_pi_r)                # length-N vector
     sign_term <- (-1)^(n + 1)               # length-N vector of alternating signs
     
     # Now we want to form an N x M matrix of exp( - (n^2*pi^2)/(2*x^2) * t )
     # Let lambda_n = (n^2*pi^2)/(2*x^2). Precompute that:
     lambda_n <- (n^2 * pi^2) / (2 * x^2)     # length-N
     
     # We can use outer() to form a matrix of size N x M:
     # exp_mat[i,j] = exp( - lambda_n[i] * t_vec[j] )
     exp_mat <- exp( - outer(lambda_n, t_vec, FUN = function(lam, tv) lam * tv) )
     
     # Now each row i of exp_mat corresponds to e^{-lambda_n[i] * t_vec}
     # We need to multiply row i by (sign_term[i] * n[i] * sin_npr[i])
     # Then sum across rows.
     #
     # Construct an N x 1 weight vector: w[i] = sign_term[i] * n[i] * sin_npr[i]
     w <- sign_term * n * sin_npr              # length-N
     
     # Multiply each row i of exp_mat by w[i]:
     # This is equivalent to w %*% exp_mat in a row‐wise fashion.
     # In R, we can do:  weighted_mat <- exp_mat * w  (where w is recycled by row)
     weighted_mat <- exp_mat * w
     
     # Now sum over rows to get a length-M vector:
     series_sum <- colSums(weighted_mat)       # length M
     
     # Finally multiply by the front constant:
     f_vals <- norm_const * series_sum
     
     return(f_vals)
}

gamma_f <- function(t,r,x){
     y <- (x-r)/(2*x)
     out <- (1 + 1/y)*exp(-2 * x^2 * (2*y+1)/t)
     return(out)
}


#----------------------------------
dinv_gamma_trunc_f <- function(t, t_f, x, r) {             # truncated to (0, t_f)
     y=(x-r)/(2*x)
     
     const <- x*y*sqrt(2) / sqrt( pi)              
     N1_f=F_b(t_f, 1/2,2*x^2*y^2)     # normalizing constant for the IG density
     
     out <-numeric(length(t))
     idx <- (t <= t_f) & (t > 0)
     out[idx] <- (const/N1_f) * t[idx]^(-3/2) * exp( -2*x^2*y^2 / t[idx] )
     
     out
}
#----------------------------------


#----------------------------------
dexp_trunc_f <- function(t, t_f, lam){    # truncated to (t_f, Inf)
     out <- numeric(length(t))             # initialise to 0
     N2_f = exp(-lam*t_f)                # normalizing constant for the Exp density
     
     idx <- (t>= t_f)
     out[idx] <- (lam / N2_f) * exp(-lam * t[idx])    # vectorised fill
     out
}
#---------------------------------

lower_bound_f <- function(t,x,r,N){
     
     a = pi^2 * t / 2 / x^2
     m0 = ceiling(1/sqrt(2*a))
     N_max = max(200,m0+N)
     
     out = f_series(t,x,r,N_max)
     remainder_bound = (pi/x/r)*(1/2/a)*exp(-a * (m0+N)^2)
     
     out = out - remainder_bound
     return(out)
}


upper_bound_f <- function(t,x,r,N){
     
     a = pi^2 * t / 2 / x^2
     m0 = ceiling(1/sqrt(2*a))
     N_max = max(200,m0+N)
     
     out = f_series(t,x,r,N_max)
     remainder_bound = (pi/x/r)*(1/2/a)*exp(-a * (m0+N)^2)
     
     out = out + remainder_bound
     return(out)
}

sample_from_f <-function(x,r){
     
     
     
     #----------------------------
     y = (x-r)/(2*x)
     t1_f = 2 * x^2 * (2*y+1)/log(1 + 1/y)
     c_f = log(2)/pi^2
     t1_f = c_f*t1_f
     lam = pi^2/(2*x^2)
     t2_f = log(2)/lam
     #----------------------------
     
     if(t1_f<t2_f){
          option=1
     }else{
          option=2
     }
     
     
     
     if(option==2){ # need only 2 components for the envelope
          
          
          # chop off some inverse gamma values
          t_f = min(t1_f, t2_f) # this should be t2_f
          
          
          #----------------------
          gam1_star = gamma_f(t1_f,r,x)
          C1_f = (x/r) * (1/(1-gam1_star))
          N1_f=F_b(t_f, 1/2,2*x^2*y^2) # normalizing constant for the IG density
          tildeC1_f = C1_f * N1_f 
          #----------------------
          
          #----------------------
          C2_f = pi^3/(6*x*r*lam)
          N2_f = exp(-lam*t_f)     # normalizing constant for the Exp density
          tildeC2_f = C2_f*N2_f
          #---------------------
          
          pp = tildeC1_f / (tildeC1_f+ tildeC2_f)
          #-------------------------------------
          
          
          sw=0
          while(!sw){
               # first, simulate from q_f()
               coin = runif(1,0,1)
               if(coin < pp){
                    #simulate from inv-gamma
                    Y = rinv_gamma_trunc(1,1/2, 2*x^2*y^2, t_f)
               }
               else{
                    #simulate from truncated exp
                    Y = rexp_trunc(1,lam,t_f,Inf)
               }
          
               

               dq_g = pp*dinv_gamma_trunc_f(Y,t_f,x,r) + (1-pp)*dexp_trunc_f(Y,t_f,lam)
               U = runif(1,0,1)
               N=200
               sw1=0
               while(!sw1){
                    N = N+50
                    LB = lower_bound_f(Y,x,r,N)
                    UB = upper_bound_f(Y,x,r,N)
                    if( U * ((tildeC1_f + tildeC2_f)*dq_g ) < LB ){
                         sw1=1 # decision done; accept
                         out=Y
                         return(out)
                    }
                    if( U * ((tildeC1_f + tildeC2_f)*dq_g ) > UB ){
                         sw1=1 # decision done, reject
                    }
                    
               } # end while sw1
          } # end while sw
          
     }else{  # now option = 1 ... need 3 components 
          
          
          #-------------------------
          gam1_star = gamma_f(t1_f,r,x)
          C1_f = (x/r) * (1/(1-gam1_star))
          N1_f=F_b(t1_f, 1/2,2*x^2*y^2) # normalizing constant for the IG density
          tildeC1_f = C1_f * N1_f 
          #-------------------------
          
          #-------------------------
          aa=pi^2*t1_f/(2*x^2)
          Nf = ceiling(1/sqrt(2*aa))
          term1 = sum((1:Nf)*exp(-aa*(1:Nf)^2))
          term2 = (1/(2*aa))*exp(-aa*Nf^2)
          C_f=(pi/(x*r))*(term1+term2)
          NN_f= (t2_f-t1_f)
          tildeC_f = C_f * NN_f
          #------------------------
          
          #----------------------
          C2_f = pi^3/(6*x*r*lam)
          N2_f = exp(-lam*t2_f)    # normalizing constant for the Exp density
          tildeC2_f = C2_f*N2_f
          #---------------------
          
          pp1 = tildeC1_f / (tildeC1_f+ tildeC2_f+tildeC_f)
          pp = tildeC_f / (tildeC1_f+ tildeC2_f+tildeC_f)
          pp2 = tildeC2_f / (tildeC1_f+ tildeC2_f+tildeC_f)
          #-------------------------------------
          
       
          sw=0
          while(!sw){
               
               # first, simulate from q_f()
               coin = runif(1,0,1)
               if(coin < pp1){
                         #simulate from inv-gamma
                    Y=rinv_gamma_trunc(1,1/2, 2*x^2*y^2, t1_f)
               }
               if((coin>pp1) & (coin <pp1+pp)){
                    Y = runif(1, min=t1_f, max=t2_f)
               }
               if (coin > pp1+pp){
                    Y = rexp_trunc(1,lam,t2_f,Inf)
               }
                    
               
               dq_g = pp1*dinv_gamma_trunc_f(Y, t1_f,x,r) + pp*dunif(Y, min=t1_f, max=t2_f)+ pp2*dexp_trunc_f(Y, t2_f, lam)
               U = runif(1,0,1)
               sw1 = 0
               N=200
               while(!sw1){
                    N = N +50
                    UB = upper_bound_f(Y,x,r,N)
                    LB = lower_bound_f(Y,x,r,N)
                    if( U * ((tildeC1_f + tildeC_f + tildeC2_f)*dq_g ) < LB ){
                         sw1=1 # decision done, accept
                         out=Y
                         return(out)
                    }
                    if( U * ((tildeC1_f + tildeC_f + tildeC2_f)*dq_g ) > UB ){
                         sw1=1 # decision done, reject
                    }
                    
               } # end while sw1
          
          } # end while sw    
          
          
     } # end of else option =1
     
     
} # end sample_from_f

sample_from_f_T <- function(x,r,T){
     if (x <= 0 || T <= 0 || r <=0) stop("Require x > 0 and r>0 and T > 0")
     repeat{
          temp = sample_from_f(x,r)
          if(temp < T){
               out=temp
               break
          }
     }
     return(out)
}
############################################
############################################






#############################################
#############################################
#############################################
# Excursion densities
joint_pdf_excursion <- function(x_grid, t_grid, T, M = 100, N = 100) {
      # Create index vectors
      m <- 1:M
      n <- 1:N
      
      # Precompute constants
      pi_sq <- pi^2
      const_factor <- sqrt(2) * pi^(9/2) * T^(3/2)
      
      # Precompute index matrices
      sign_mat <- (-1)^(outer(m, n, "+"))
      m2n2_mat <- outer(m^2, n^2)
      
      # Create grid
      grid <- expand.grid(x = x_grid, t = t_grid)
      
      # Vectorized evaluation over the grid
      values <- mapply(function(x, t) {
            denom <- x^6
            exp_part <- exp(
                  - (outer(m^2, rep(1, N)) * pi_sq * (T - t) +
                           outer(rep(1, M), n^2) * pi_sq * t) / (2 * x^2)
            )
            sum_term <- sum(sign_mat * m2n2_mat * exp_part)
            const_factor / denom * sum_term
      }, grid$x, grid$t)
      
      # Return tidy data frame
      data.frame(x = grid$x, t = grid$t, value = values)
}

marginal_M_excursion <- function(x_vec, T, N = 100) {
      # Vector y_vec: vector of y > 0
      # N: number of terms in sum (truncate at ±N)
      
      n <- 1:N  # Only use positive n since sum is even
      n_mat <- matrix(n, nrow = length(x_vec), ncol = N, byrow = TRUE)
      y_mat <- matrix(x_vec, nrow = length(x_vec), ncol = N)
      
      n2 <- n_mat^2
      y2 <- y_mat^2
      
      term <- n2 * (3 - 4 * n2 * y2/T) * exp(-2 * n2 * y2/T)
      sum_terms <- rowSums(term)
      
      f_y <- -(8/T) * x_vec * sum_terms  # Multiply by 2 for negative n symmetry, total -4 × 2 = -8
      
      return(f_y)
}


marginal_tau_excursion <- function(t, T, max_terms = 200) {
      # Input validation
      if (any(t <= 0) || any(t >= T) || T <= 0) {
            stop("t must be in (0, T) and T > 0")
      }
      
      t <- as.vector(t)
      
      # Create sequences
      m_seq <- 1:max_terms
      n_seq <- 1:max_terms
      
      # Create all combinations using outer product
      m_matrix <- outer(m_seq, n_seq, function(x, y) x)
      n_matrix <- outer(m_seq, n_seq, function(x, y) y)
      
      # Pre-compute signs and squares
      sign_matrix <- (-1)^(m_matrix + n_matrix)
      m_sq_matrix <- m_matrix^2
      n_sq_matrix <- n_matrix^2
      
      # Vectorized computation for all t values
      result <- sapply(t, function(t_val) {
            # Compute denominators for this t
            denom_matrix <- (n_sq_matrix * t_val + m_sq_matrix * (T - t_val))^(5/2)
            
            # Compute all terms
            term_matrix <- sign_matrix * m_sq_matrix * n_sq_matrix / denom_matrix
            
            # Sum all terms
            series_sum <- sum(term_matrix)
            
            # Apply the 3T^{3/2} factor
            return(3 * T^(3/2) * series_sum)
      })
      
      return(result)
}


###############################################
###############################################
###############################################
## Restricted meander densities
#---------------------------------------------
joint_pdf_meander_r <- function(x_grid, t_grid, T, r, M = 200, N = 200) {
      # Indices
      m <- 1:M
      n <- 1:N
      
      # Constants
      pi2 <- pi^2
      
      # Outer products for precomputation
      sign_mat <- (-1)^(outer(m, n, "+"))     # M x N
      mn2_mat  <- outer(m, n^2)               # M x N
      m_mat    <- matrix(rep(m, N), M, N)     # M x N
      n2_mat   <- matrix(rep(n^2, each = M), M, N)
      m2_mat   <- matrix(rep(m^2, N), M, N)
      
      # Create grid
      grid <- expand.grid(x = x_grid, t = t_grid)
      
      # Evaluate function over grid
      grid$values <- mapply(function(x, t) {
            sin_term <- sin(m * pi * r / x)       # vector of length M
            sin_mat  <- matrix(rep(sin_term, N), M, N)
            
            exp_mat <- exp(
                  -pi2 / (2 * x^2) * (n2_mat * t + m2_mat * (T - t))
            )
            
            sum_term <- sum(sign_mat * mn2_mat * sin_mat * exp_mat)
            sum_term / x^5
      }, grid$x, grid$t)
      
      # Return tidy data frame
      data.frame(x = grid$x, t = grid$t, value = grid$values)
}


marginal_M_meander_r = function(x,r,T) {
      # x: vector of values
      # r, T: parameters
      # K: truncation parameter for sum over k from -K to K (excluding 0)
      K=100
      
      k_vals <- setdiff(-K:K, 0)  # exclude k = 0 since 2k = 0 contributes nothing
      
      # Vectorized outer product: each row is an x, each column is a k
      k_matrix <- outer(x, k_vals, function(xi, ki) 2 * ki * xi + r)
      
      # Compute the common squared terms
      z2 <- k_matrix^2
      
      # Evaluate the exponential and polynomial part
      exp_part <- exp(-z2 / (2 * T))
      poly_part <- 1 - z2 / T
      
      # Full summand: 2k * (...) * exp(...)
      summand <- sweep(exp_part * poly_part, 2, 2 * k_vals, FUN = "*")
      
      # Sum over k for each x
      sum_over_k <- rowSums(summand)
      
      # Final scaling
      result <- exp(r^2 / (2*T)) / r * sum_over_k
      return(result)
}





###############################################
###############################################
###############################################
## Brownaian meander densities
#---------------------------------------------
joint_pdf_meander <- function(x_vals, t_vals, T, M = 50, N = 50) {
      # Precompute constants
      const <- sqrt(2) * pi^(5/2) * sqrt(T)
      pi2 <- pi^2
      
      # Create m and n vectors
      m <- 1:M
      n <- 1:N
      
      # Precompute sign matrix: sgn[m, n] = (-1)^{m+n} - (-1)^n
      sign_mat <- outer(m, n, function(m, n) (-1)^(m + n) - (-1)^n)
      
      # Precompute n^2
      n_sq <- n^2
      
      # Initialize storage
      result <- data.frame(x = numeric(), t = numeric(), value = numeric())
      
      for (x in x_vals) {
            x2 <- x^2
            x4 <- x^4
            for (t in t_vals) {
                  # Compute exponential matrix
                  exp_mat <- outer(
                        m, n,
                        function(m, n) {
                              arg <- - (n^2 * pi2 * t + m^2 * pi2 * (T - t)) / (2 * x2)
                              exp(arg)
                        }
                  )
                  
                  # Multiply element-wise: sign * n^2 * exp(...)
                  sum_mat <- sign_mat * rep(n_sq, each = M) * exp_mat
                  
                  # Sum all values
                  total <- sum(sum_mat)
                  
                  # Final function value
                  f_xt <- (const / x4) * total
                  
                  # Store result
                  result <- rbind(result, data.frame(x = x, t = t, value = f_xt))
            }
      }
      
      return(result)
}



marginal_tau_meander = function(t_grid,T) {
      
      M=300
      N=300
      
      # m, n indices
      m <- 1:M
      n <- 1:N
      m_sq <- m^2
      n_sq <- n^2
      
      alpha=0.98
      
      # Precompute sign matrix
      sign_mat <- outer(m, n, function(m, n) (-1)^(m + n) - (-1)^n)
      alpha_mat <- outer(m,n, function(m,n) alpha^(m+n) )
      
      
      # Numerator matrix: n^2 for each entry
      num_mat <- outer(m, n, function(m, n) n^2)
      
      # Result vector
      result <- numeric(length(t_grid))
      
      for (i in seq_along(t_grid)) {
            t <- t_grid[i]
            denom_mat <- outer(m_sq, n_sq, function(m2, n2) n2 * t + m2 * (T - t))
            term_mat <- alpha_mat * sign_mat * num_mat / denom_mat^(3/2)
            result[i] <- sqrt(T) * sum(term_mat)
      }
      
      return(result)
}


marginal_M_meander <- function(x_vals,T) {
      
      N=100
      # x_vals: vector of x values
      # T: parameter
      # N: truncation index for the sum
      
      # Precompute constants
      const <- 2^(3/2) * sqrt(pi) * sqrt(T)
      pi2 <- pi^2
      k_vals <- (2 * (0:(N - 1)) + 1)^2  # (2n + 1)^2 for n = 0 to N-1
      
      # Initialize result vector
      result <- numeric(length(x_vals))
      
      for (i in seq_along(x_vals)) {
            x <- x_vals[i]
            x2 <- x^2
            x4 <- x^4
            exp_term <- exp(- k_vals * pi2 * T / (2 * x2))
            sum_term <- (-1 / x2 + k_vals * pi2 * T / x4) * exp_term
            result[i] <- const * sum(sum_term)
      }
      
      return(result)
}



######################################
######################################
#Devroye's MAXMEANDER algorithm for t=1, r is given
#Function to simulate the maximum value of BM|(BM>0,BM(t)=r) attains b/w [0,t]

MAXMEANDER <- function(r){
  Xi=6.8*exp(-9);Zeta=2.2*exp(-9);
  #Regime 1
  if(r>=1.5){
    c=5*r/((10*r^2)-8)
    #decision value 0=Reject,1=Accept,-1=Undecided
    DECISION=-1
    #REPEAT
    while(DECISION<1){
      E=rexp(n=1,rate=1);X=r+(c*E);
      V=runif(n=1,min=0,max=1);Y=10*r*V*exp(-1*E)/(1-Xi);
      k=2;S=f_k(r,1,X);
      DECISION=-1
      #REPEAT
      while(DECISION<0){
        #lower bound at step k
        L_k=S-(4*k*(1+(4*k*X*r))*exp((-2*(k*X)^2)+(2*k*X*r))/((1-Zeta)*r))
        #upper bound at step k
        U_k=S+((2*k*(r+((4*(k*X)^2)/r)))*exp((-2*(k*X)^2)+(2*k*X*r))/(1-Xi))
        if(Y<=L_k){DECISION=1;break;}
        if(Y>=U_k){DECISION=0;break;}
        S=S+f_k(r,k,X)
        k=k+1
      }#Until != Undecided
      
    }#Until Accept
    return(X)
  }#regime 1 ends
  
  #Regime 2
  if(r<1.5){
    mu=16*exp(-2*(pi^2)/3);
    nu=16*exp(-9);tau=4*exp(-9);
    p=(3*exp(9/8))/(1-mu); q=(123*exp((3*r)-(9/2)))/(2*(1-nu)*(4-(2*r)));
    DECISION=-1
    #REPEAT
    while(DECISION<1){
      U=runif(n=1,min=0,max=1);V=runif(n=1,min=0,max=1);
      if(U<=(p/(p+q))){
        N=rnorm(n=1,mean=0,sd=1);E1=rexp(n=1,rate=1);E2=rexp(n=1,rate=1);
        X=pi/sqrt((N^2)+(2*E1)+(2*E2));
        g_X=sqrt(2*pi)*exp(9/8)*(pi^4)*exp(-1*(pi^2)/(2*(X^2)))/((1-mu)*(X^6));
        Y=V*g_X;
        k=2;S=psi_k(r,1,X);
        DECISION=-1;
        #REPEAT
        while(DECISION<0){
          if(X>=1.5){DECISION=0;break;}
          T1=(1/(1-mu))*F_k(r,k,X)*(((k*pi)^3)*r/(X^4));
          if(Y<=(S-T1)){DECISION=1;break;}
          if(Y>=(S+T1)){DECISION=0;break;}
          S=S+psi_k(r,k,X)
          k=k+1
        }#Until != Undecided
      }
      else{
        E=rexp(n=1,rate=1);
        X=(1.5)+(E/(4-(2*r)));
        g1_X=q*(4-(2*r))*exp(-1*(4-(2*r))*(X-1.5));
        Y=V*g1_X;
        k=2;S=f_k(r,1,X);
        DECISION=-1;
        #REPEAT
        while(DECISION<0){
          L_k=S-(8*k*k*X*exp((2*k*X*r)-(2*k*k*X*X))/(1-tau));
          U_k=S+(164*(k^4)*(X^3)*exp((2*k*X*r)-(2*k*k*X*X))/(9*(1-nu)));
          if(Y<=L_k){DECISION=1;break;}
          if(Y>=U_k){DECISION=0;break;}
          S=S+f_k(r,k,X)
          k=k+1
        }#Until != Undecided
      }
    }#Until Accept
    
    return(X)
  }#regime 2 ends
}#MAXMEANDER ends

f_k<-function(r,k,x)
{
  f_k_x=2*k*exp(-2*((k*x)^2))*(((1-(r+(2*k*x))^2)*exp(-2*k*x*r))-((1-(r-(2*k*x))^2)*exp(2*k*x*r)))/r
  if(x<r){f_k_x=0}
  return(f_k_x)
}

psi_k<-function(r,k,x){
  psi_k_x=(F_k(r,k,x))*(((((k*k*pi*pi)-(2*x*x))/(x^3))*sin(pi*k*r/x))-((pi*k*r/(x^2))*cos(pi*k*r/x)))
  return(psi_k_x)
}

F_k<-function(r,k,x){
  F_k_x=((sqrt(2*pi)*exp(r*r/2)*pi*k*exp(-1*k*k*pi*pi/(2*x*x)))/(x*x*r))
  if(x<r){F_k_x=0}
  return(F_k_x)
}
#End of MAXMEANDER code

