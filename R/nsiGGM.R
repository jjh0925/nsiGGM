## Sunghwan' lasso ###
## Main algoirhtm

## " Ind.Path.gene " Should be included.

nsiGGM <- function(Y,
                   lambda1,
                   lambda2,
                   lambda3,
                   params.in_gr,
                   params.out_gr,
                   params.dummy,
                   params.verbose,
                   penalty="group",
                   rho=1,
                   penalize.diagonal=FALSE,
                   maxiter=50,
                   tol=1e-5,
                   warm=NULL,
                   screening="fast",
                   truncate=1e-5,
                   is.Naive.Group=TRUE
)
{
   ## initialize:
   p = dim(Y[[1]])[2]
   K = length(Y)
   n = rep(0, K)
   for(k in 1:K) {
      n[k] = dim(Y[[k]])[1]
   }

   ## assign feature names if none exist:
   if(length(dimnames(Y[[1]])[[2]])==0)
   {
      for(k in 1:K)
      {
         dimnames(Y[[k]])[[2]] = paste("V",1:p,sep="")
      }
   }

   # Centeralizaion over zero of Y: (i.e. mean values of each column of Y is equal to zero.)
   for(k in 1:K){
      Y[[k]] = (scale((Y[[k]]), scale=FALSE))
   }

   # set weights:
   weights = rep(1,K)

   ## if p is very high, identify connected nodes without storing full S matrices
   ### now get criterion over connected S:
   ## define S
   S = vector("list",length=K)

   for(k in 1:K)
   {
      ntemp = dim(Y[[k]])[1]
      S[[k]] = cov(Y[[k]])*(ntemp-1)/ntemp   ### <- This seems the estimation by MLE. S : Covariance matrix
   }

   ## define theta on all connected:   (so theta is really theta.connected).
   theta = list()
   for(k in 1:K)
   {
      theta[[k]] = matrix(0, p, p)
      dimnames(theta[[k]])[[1]] = dimnames(theta[[k]])[[2]] = dimnames(Y[[k]])[[2]]  ## by which we can name feature set to each data set.
   }

   ## now run JGL on each block of the "connected nodes" to fill in theta:
   Thetabl = admm.iters_in_out(Y, lambda1, lambda2, lambda3, params.in_gr, params.out_gr, params.dummy, params.verbose ,penalty=penalty, rho=rho, weights=weights, maxiter=maxiter, tol=tol, warm=warm.bl, is.Naive.Group=is.Naive.Group)

   # Update Theta with Thetabl's results:
   for(k in 1:K) {
      theta[[k]] = Thetabl$Z[[k]]
   }

   # round very small theta entries down to zero:
   if(dim(theta[[1]])[1]>0)
   {
      for(k in 1:K)
      {
         rounddown = abs(theta[[k]])<truncate; diag(rounddown)=FALSE
         theta[[k]]=theta[[k]]*(1-rounddown)
      }
   }

   out = list(theta = theta,  S = S)
   # class(out)="jgl"
   return(out)
}





admm.iters_in_out = function(Y,
                             lambda1,
                             lambda2,
                             lambda3,
                             params.in_gr,
                             params.out_gr,
                             params.dummy,
                             params.verbose,
                             penalty="group",
                             rho=1,
                             rho.increment=1,
                             weights,
                             maxiter = maxiter,
                             tol=1e-3,
                             warm=NULL,
                             is.Naive.Group=is.Naive.Group){
   K = length(Y)
   p = dim(Y[[1]])[2]
   n = weights

   ns = c()
   for(k in 1:K){ns[k] = dim(Y[[k]])[1]}

   S = list()
   for(k in 1:K){S[[k]] = cov(Y[[k]])*(ns[k]-1)/ns[k]}

   # initialize theta:  = Inverse covariance
   theta = list()
   for(k in 1:K){theta[[k]] = diag(1/diag(S[[k]]))}
   # initialize Z:
   Z = list(); for(k in 1:K){Z[[k]]=matrix(0,p,p)}
   # initialize W: (equal to U in the literature)
   W = list();  for(k in 1:K) {W[[k]] = matrix(0,p,p)}


   # iterations:
   iter = 0
   diff_value = 10

   while((iter==0) || (iter < maxiter && diff_value > tol))
   {

      theta_prev = theta
      # %--------------------------------------------------------------------------------%
      ## (1) Inference steps for Theta
      for(k in 1:K){
         edecomp = eigen(S[[k]] - rho*Z[[k]]/n[k] + rho*W[[k]]/n[k]) # In the slide #30, the equation at the bottom.
         D = edecomp$values
         V = edecomp$vectors
         D2 = n[k]/(2*rho) * ( -D + sqrt(D^2 + 4*rho/n[k]) )         # In the slide #30, the equation at the bottom.
         theta[[k]] = matrix(as.double(V %*% diag(D2) %*% t(V)), dim(V %*% diag(D2) %*% t(V))[1], dim(V %*% diag(D2) %*% t(V))[1])   ##  %*%  = matrix multiplication.
      }

      # %--------------------------------------------------------------------------------%
      ## (2) update Z:
      ## Define A matrices:
      A = list()

      for(k in 1:K){
         A[[k]] = theta[[k]] + W[[k]]
      }

      #if(penalty=="group")
      #{
      #if(is.Naive.Group){

      # Z = Structured_dsgl(A, rho, lambda1, lambda2, penalize.diagonal, Ind.Path.gene)

      #}else{

      #   minimize rho/2 ||Z-A||_F^2 + P(Z):
      #  Z = dsgl(A, rho, lam1, lam2, penalize.diagonal=TRUE)

      ##########################################
      ## Check up again the diagonal factors ###
      ##########################################

      tran.A <- foreach(i = 1:length(A), .combine=rbind) %do% { A[[i]][lower.tri(A[[i]], diag = TRUE)] }
      # tran.A <- t(do.call(cbind,lapply(A, as.vector)))
      ### Approximation step gets involved in that we draw an group effect
      X <- diag(dim(tran.A)[1]) # + matrix(rnorm(dim(tran.A)[1]^2,sd = 0.00001),dim(tran.A)[1],dim(tran.A)[1])

      dummy = 1 # % end of index for individual terms
      N = dim(X)[1]
      dims = dim(X)[2]
      T = dim(tran.A)[2]

      params.lambda1 = params.lambda4 = lambda1
      params.lambda2 = lambda2
      params.lambda3 = lambda3


      estimated_beta = struct_io_lasso(X, tran.A,  params.lambda1, params.lambda2, params.lambda3, params.lambda4, params.weight, params.in_gr_w
                                       ,params.out_gr_w, params.in_gr, params.out_gr, params.dummy,params.verbose, is.standard=FALSE)



      Z <- foreach(i=1:length(A)) %do% {
         tmp <- rep(NA, length(A[[1]]))
         tmp[as.vector(lower.tri(A[[1]], diag = TRUE))] <- estimated_beta[i,]
         tmpmat <- matrix(tmp, dim(A[[1]])[1], dim(A[[1]])[1])
         tmp[is.na(tmp)] <- t(tmpmat)[upper.tri(t(tmpmat))]
         tmpmat <- matrix(tmp, dim(A[[1]])[1], dim(A[[1]])[1])
         return(tmpmat)
      }
      # }
      #}

      # %--------------------------------------------------------------------------------% By this line, we complete the procedure for Z

      # update the dual variable W:
      for(k in 1:K){W[[k]] = W[[k]] + (theta[[k]]-Z[[k]])}

      # %--------------------------------------------------------------------------------% By this line, we complete the procedure for Z

      # bookkeeping:
      print(paste("ADMM",iter))
      iter = iter+1
      diff_value = 0
      for(k in 1:K) {diff_value = diff_value + sum(abs(theta[[k]] - theta_prev[[k]])) / sum(abs(theta_prev[[k]]))}
      # increment rho by a constant factor:
      rho = rho*rho.increment
   }
   diff = 0; for(k in 1:K){diff = diff + sum(abs(theta[[k]]-Z[[k]]))}
   out = list(theta=theta,Z=Z,diff=diff,iters=iter)
   return(out)
}


#################################
######   CMU lasso  #############
#################################

#estimated_beta = struct_io_lasso(X, tran.A,  params.lambda1, params.lambda2, params.lambda3, params.lambda4, params.weight, params.in_gr_w
#   ,params.out_gr_w, params.in_gr, params.out_gr, params.dummy, is.standard=FALSE)

struct_io_lasso <- function(D,
                            Y,
                            params.lambda1,
                            params.lambda2,
                            params.lambda3,
                            params.lambda4,
                            params.weight,
                            params.in_gr_w,
                            params.out_gr_w,
                            params.in_gr,
                            params.out_gr,
                            params.dummy,
                            params.verbose,
                            is.standard=TRUE){

   # Seunghak Lee (seunghak@cs.cmu.edu)
   # X: Design matrix (N by J matrix, N: # of samples, J: # of inputs)
   #Y: Response matrix (N by K matrix, N: # of samples, K: # of outputs)
   #params: parameters for our method
   #    <tuning parameters>
   #      params.lambda1
   #      params.lambda2
   #      params.lambda3
   #      params.lambda4
   #    <input group>
   #      params.in_gr: cell type
   #      (All inputs should belong to at least one group. It can belong to a group by itself.)
   #    <output group>
   #      params.out_gr: cell type
   #      (All inputs should belong to at least one group. It can belong to a group by itself.)
   # nargin tells us how many arguements are required.
   #                                          if nargin < 3
   #                                          error('Three arguments are needed; struct_io_lasso(D, Y, params)');
   #                                          end
   #if ~isfield(params,'lambda1') % tuning penalty for individual term
   #params.lambda1 = 0.01;
   #end
   #if ~isfield(params,'lambda2') % tuning penalty for input grouping
   #params.lambda2 = 0;
   #end
   #if ~isfield(params,'lambda3') % tuning penalty for output grouping
   #params.lambda3 = 0;
   #end
   #if ~isfield(params,'lambda4') % tuning penalty for interaction term
   #params.lambda4 = 0;
   #end
   #if ~isfield(params,'verbose')
   #params.verbose = 1;  % 1; verbose
   #end

   params.stopVal = 10e-8;
   params.MaxIterNum = 50;
   params.mode = 1
   params.init_beta = list()
   params.w = list()
   params.rho = list()
   params.nu = list()


   #standardize

   Y <- as.matrix(Y)   # tran.A
   D <- as.matrix(D)   # X

   if(is.standard){
      Y <- standardize_mat(Y)
      D <- standardize_mat(D)
   }

   dims = dim(D)[2] # % dimension
   N = dim(D)[1] #% samples
   T = dim(Y)[2] #% number of traits

   lambda1 = params.lambda1;
   lambda2 = params.lambda2;
   lambda3 = params.lambda3;
   lambda4 = params.lambda4;

   in_gr = params.in_gr;
   out_gr = params.out_gr;

   weight = matrix(rep(1,dims*T), dims,T)
   in_gr_w = matrix(rep(1,length(in_gr)), length(in_gr),1)
   out_gr_w = matrix(rep(1,length(out_gr)), length(out_gr),1)

   params.in_gr_w = in_gr_w;
   params.out_gr_w = out_gr_w;
   params.weight = weight;

   init_beta = params.init_beta;

   if(is.null(params.dummy)){
      dummy = dims
   } else {
      dummy = params.dummy
   }

   verbose = params.verbose;
   stopVal = params.stopVal;


   ######################################################
   in_gr_dic = list() ## It is supposed to contain the elements as the number of X dimension.
   tmp <- c()
   for (i in 1:length(in_gr)) {
      tmp <-c(tmp, rep(i, length(in_gr[[i]])))
   }

   n.tmp <- c()
   for(i in 1:length(unlist(in_gr))){
      n.tmp <- c(n.tmp, table(unlist(in_gr))[names(table(unlist(in_gr))) == unlist(in_gr)[i]])
   }

   for(i in 1:length(tmp)){
      k <- unlist(in_gr)[i]
      if(n.tmp[i] > 1){
         in_gr_dic[[k]] <-  c(tmp[i]-1, tmp[i])
      }else {
         in_gr_dic[[k]] <-  tmp[i]
      }
   }

   ###############################################################
   out_gr_dic = list() ## It is supposed to contain the elements as the number of X dimension.

   tmp <- c()
   for (i in 1:length(out_gr)) {
      tmp <-c(tmp, rep(i, length(out_gr[[i]])))
   }

   n.tmp <- c()
   for(i in 1:length(unlist(out_gr))){
      n.tmp <- c(n.tmp, table(unlist(out_gr))[names(table(unlist(out_gr))) == unlist(out_gr)[i]])
   }

   for(i in 1:length(tmp)){
      k <- unlist(out_gr)[i]
      if(n.tmp[i] > 1){
         out_gr_dic[[k]] <-  c(tmp[i]-1, tmp[i])
      }else {
         out_gr_dic[[k]] <-  tmp[i]
      }
   }

   ###############################################################


   # Solution for ridge-regression
   reg_param = (lambda1 + lambda2 + lambda3)*10;  ## These three parameter refers to three lambda
   beta = solve(t(D) %*% D + reg_param*diag(dims)) %*% (t(D)%*%Y) # This is the ridge regression solution.


   ## It seems to be an arguemnt for the objective function (5a), (5b), (5c).
   obj_value = objFnSIO(D, Y, beta,  params.lambda1, params.lambda2, params.lambda3, params.lambda4, params.weight, params.in_gr_w
                        ,params.out_gr_w, params.in_gr, params.out_gr, params.dummy)



   # Computing residual
   rsd_tild = Y - D %*% beta # N by K matrix

   ##########################################################################################

   for(outer_loop in 1:10000){

      for(m in 1:length(out_gr)){    ## out_gr has 1:10

         ogrp_size = length(out_gr[[m]])
         out_lambda = lambda3*out_gr_w[[m]]

         for(m2 in 1:length(in_gr)){

            igrp_size = length(in_gr[[m2]])
            in_lambda = lambda2*in_gr_w[m2, 1]  ##  Note that Original in_gr_w[m, 1]   here m seesm to be typo!!!!

            # Check if all entries are zero
            cond_1 = 0
            zero_chk = sum(sum(abs(beta[in_gr[[m2]], out_gr[[m]]])))

            ## This step seems very important in that at once we have all zero entries in a certain group, we skip the iteration to another.
            ## This step is involving the eita(n) and gamma(r) to get rid of beta that has all zero values at a certain group g and h. in eq(10)


            if(zero_chk != 0){
               #% The below process covers
               rsd_tild2 = rsd_tild;
               ## only the groupped features of output and input and its corresponding residual.

               ## Notice that this is the code for implementing  r = y - sum(bx)_(l!=j). See that by adding sum(bx)_(l!=j), we get what we want to see.
               rsd_tild2[ ,out_gr[[m]]] = rsd_tild2[ ,out_gr[[m]]] + D[ ,in_gr[[m2]], drop=FALSE] %*% beta[in_gr[[m2]], out_gr[[m]], drop=FALSE];


               ##############################################
               # For eq.(10)
               ##############################################
               if( in_lambda > 0 && out_lambda > 0 ){
                  for(in_gr_idx in 1:length(in_gr[[m2]])){

                     j = in_gr[[m2]][in_gr_idx]
                     for(out_gr_idx in 1:ogrp_size){   ##  ogrp_size is # of output group in our initial example dataset, output group size is 1.

                        if(j <= dummy){
                           lambda = lambda1 * weight[j, out_gr[[m]][out_gr_idx] ]
                        } else {
                           lambda = lambda4 * weight[j, out_gr[[m]][out_gr_idx] ]
                        }

                        if(beta[j, out_gr[[m]][out_gr_idx] ] != 0){
                           tmp2 = t(D[ ,j,drop=FALSE])  %*% rsd_tild2[ ,out_gr[[m]][out_gr_idx], drop=FALSE]    ## This is r*t(x) in the equation (10) in detailed paper.

                           if(abs(tmp2) <= lambda){
                           } else {
                              s_jk = sign(tmp2)
                              cond_1 = cond_1 + (tmp2 - lambda*s_jk)^2; ## This seems to have to do with the left hand side of eq.(10) in full-detailed manuscript.
                           }
                        }
                     }
                  }

                  block_thres = (in_lambda * sqrt(ogrp_size) + out_lambda * sqrt(igrp_size))^2  ## the right handside part in (10)

                  if(cond_1 <= block_thres){
                     beta[in_gr[[m2]], out_gr[[m]]] = 0
                     rsd_tild = rsd_tild2;
                  }
               }

               ##############################################
               # For eq.(11) Input
               ##############################################

               if(in_lambda > 0){
                  for(out_gr_idx in 1:ogrp_size){

                     cond_1 = 0
                     denom = sqrt(sum(beta[in_gr[[m2]],out_gr[[m]][out_gr_idx]]^2))

                     if(denom != 0){

                        rsd_tild2 = rsd_tild
                        rsd_tild2[ ,out_gr[[m]][out_gr_idx]] = rsd_tild2[ ,out_gr[[m]][out_gr_idx]] + D[ ,in_gr[[m2]], drop=FALSE] %*% beta[ in_gr[[m2]], out_gr[[m]][out_gr_idx], drop=FALSE]


                        for(in_gr_idx in 1:length(in_gr[[m2]])){
                           j = in_gr[[m2]][in_gr_idx]
                           if(j <= dummy){
                              lambda = lambda1 * weight[j, out_gr[[m]][out_gr_idx] ]
                           } else {
                              lambda = lambda4 * weight[j, out_gr[[m]][out_gr_idx] ]
                           }

                           if(beta[j, out_gr[[m]][out_gr_idx] ] != 0){
                              tmp2 = t(D[ ,j,drop=FALSE])  %*% rsd_tild2[ ,out_gr[[m]][out_gr_idx], drop=FALSE]    ## This is r*t(x) in the equation (11) in detailed paper.

                              if(abs(tmp2) <= lambda){
                              }else {
                                 s_jk = sign(tmp2)
                                 cond_1 = cond_1 + (tmp2 - lambda*s_jk)^2; ## This seems to have to do with eq.(11) in full-detailed manuscript.
                              }
                           } # This is for if(beta[j, out_gr[[m]][out_gr_idx] ] != 0){
                        }


                        if(cond_1 <= in_lambda^2) {
                           beta[in_gr[[m2]], out_gr[[m]][out_gr_idx]] = 0;
                           rsd_tild = rsd_tild2;

                        } else {
                           ################################################
                           #### From this non-zero index, we probably use another way of estimation scheme (direct group lasso or graph guided-lasso)

                           nonzero_idx = (beta[ in_gr[[m2]], out_gr[[m]][out_gr_idx]]  !=0 )  ##
                           num_nonzero = sum(nonzero_idx)  ##

                           for(in_gr_idx in 1:length(in_gr[[m2]])){
                              j = in_gr[[m2]][in_gr_idx]
                              zero_chk = abs(beta[j, out_gr[[m]][out_gr_idx]])

                              if(zero_chk == 0){
                              }else{

                                 if(num_nonzero == 1 && nonzero_idx[in_gr_idx] == 1){  ## I guess If only one feature is non-zero, nothing is to be updated more.
                                 }else{

                                    nonzero_idx2 = ( beta[j, out_gr[[m]]] != 0)
                                    num_nonzero2 = sum(nonzero_idx2);

                                    if(num_nonzero2 == 1 && nonzero_idx2[out_gr_idx] == 1){  ## I guess If only one feature is non-zero, nothing is to be updated more.
                                    }else{


                                       if(j <= dummy){
                                          lambda = lambda1 * weight[j, out_gr[[m]][out_gr_idx]]
                                       } else {
                                          lambda = lambda4 * weight[j, out_gr[[m]][out_gr_idx]]
                                       }

                                       rsd_tild2 = rsd_tild

                                       # Unlike the eq above (rsd_tild2[ ,out_gr[[m]][out_gr_idx]] = rsd_tild2[ ,out_gr[[m]][out_gr_idx]] + D[ ,in_gr[[m2]]] %*% beta[ in_gr[[m2]], out_gr[[m]][out_gr_idx]])
                                       # Here we only consider only one column of output index
                                       rsd_tild2[ ,out_gr[[m]][out_gr_idx]] = rsd_tild2[ ,out_gr[[m]][out_gr_idx]] + D[ ,j,drop=FALSE] %*% beta[j, out_gr[[m]][out_gr_idx], drop=FALSE]
                                       cond_3 = t(D[ ,j,drop=FALSE]) %*% rsd_tild2[ ,out_gr[[m]][out_gr_idx], drop=FALSE]

                                       if( abs(cond_3) <= lambda ) {  ## Probably, abs(cond_3) is only one I can use in JGL.

                                          beta[j,  out_gr[[m]][out_gr_idx]] = 0
                                          rsd_tild = rsd_tild2

                                       } else {
                                          rel_in_gr = in_gr_dic[[j]]  ## Group lasso technique involves
                                          rel_out_gr = out_gr_dic[[ out_gr[[m]][out_gr_idx]]]

                                          q2 = 0
                                          for(ridx in 1:length(rel_in_gr)){
                                             # In theory, in_gr[[rel_in_gr[ridx]]] function to bring an simulaneous zero value because it uses multiple groups which has any overlapping features.
                                             # beta[in_gr[[rel_in_gr[ridx]]], out_gr[[m]][out_gr_idx]]^2) contains involving all entries in this iteration.
                                             q2 = q2 + in_gr_w[rel_in_gr[ridx]] * lambda2 / sqrt(sum(beta[in_gr[[rel_in_gr[ridx]]], out_gr[[m]][out_gr_idx]]^2))
                                          }

                                          w3 = 0
                                          for(ridx in 1:length(rel_out_gr)){
                                             w3 = w3 + out_gr_w[rel_out_gr[ridx]] * lambda3 / sqrt(sum(beta[j, out_gr[[rel_out_gr[ridx]]]]^2))
                                          }

                                          beta_minus = min(0, 1/(1+q2+w3)*(cond_3 + lambda)) ## Note that cond_3 uses all entries in a correspoding column but it excludes the effect of beta_k^j when in computing its residual
                                          beta_plus  = max(0, 1/(1+q2+w3)*(cond_3 - lambda))

                                          new_beta = beta_minus + beta_plus
                                          rsd_tild[ ,out_gr[[m]][out_gr_idx]] = rsd_tild2[ ,out_gr[[m]][out_gr_idx]] - D[ ,j, drop=FALSE] %*% new_beta   ## Update residual : Since we didn't consider target features when computing residual, we subtract out the newly calculated residuals  ##
                                          beta[j, out_gr[[m]][out_gr_idx]] = new_beta


                                       } # This is for   if( abs(cond_3) <= lambda ) {


                                    } # This is for   if(num_nonzero2 == 1 && nonzero_idx2(out_gr_idx) == 1){
                                 } # This is for   if(num_nonzero == 1 && nonzero_idx(in_gr_idx) == 1)
                              } # This is for   if(zero_chk == 0){
                           } # This is for   for(in_gr_idx in 1:length(in_gr[[m2]]))
                        } # This is for   if(cond_1 <= in_lambda^2) {
                     } # This is for if(denom != 0){
                  } # This is for for(out_gr_idx in 1:ogrp_size)
               } # This is for if(in_lambda > 0){



               ##############################################
               # For eq.(12) Output
               ##############################################

               if(out_lambda > 0){
                  for(in_gr_idx in 1:length(in_gr[[m2]])) {   ## This is contray to one above.

                     j = in_gr[[m2]][in_gr_idx]
                     denom = sqrt(sum(beta[j, out_gr[[m]]]^2))
                     if(denom != 0){

                        rsd_tild2 = rsd_tild
                        rsd_tild2[ ,out_gr[[m]]] = rsd_tild2[ ,out_gr[[m]]] + D[ ,j, drop=FALSE] %*% beta[j, out_gr[[m]], drop=FALSE]
                        cond_1 = 0

                        for(out_gr_idx in 1:ogrp_size){
                           if(j <= dummy){
                              lambda = lambda1 * weight[j, out_gr[[m]][out_gr_idx] ]
                           } else {
                              lambda = lambda4 * weight[j, out_gr[[m]][out_gr_idx] ]
                           }

                           zero_chk = abs(beta[j, out_gr[[m]][out_gr_idx]])
                           if(zero_chk != 0){

                              tmp2 = t(D[ ,j,drop=FALSE])  %*% rsd_tild2[ ,out_gr[[m]][out_gr_idx], drop=FALSE]    ## This is r*t(x) in the equation

                              if (abs(tmp2) <= lambda){
                              } else {
                                 s_jk = sign(tmp2)
                                 cond_1 = cond_1 + (tmp2 - lambda*s_jk)^2
                              }
                           } # This is for if(zero_chk != 0){
                        }

                        if(cond_1 <= out_lambda^2){
                           beta[j, out_gr[[m]] ] = 0
                           rsd_tild = rsd_tild2
                        } else {
                           nonzero_idx = (beta[j, out_gr[[m]]] !=0  );
                           num_nonzero = sum(nonzero_idx);


                           for(out_gr_idx in 1:ogrp_size){
                              if(j <= dummy){
                                 lambda = lambda1 * weight[j, out_gr[[m]][out_gr_idx]]
                              } else{
                                 lambda = lambda4 * weight[j, out_gr[[m]][out_gr_idx]]
                              }

                              zero_chk = abs(beta[j, out_gr[[m]][out_gr_idx]])

                              if(zero_chk != 0){

                                 if(num_nonzero == 1 && nonzero_idx[out_gr_idx] == 1){
                                 }else{

                                    nonzero_idx2 = (beta[in_gr[[m2]], out_gr[[m]][out_gr_idx]]!=0)
                                    num_nonzero2 = sum(nonzero_idx2)
                                    if(num_nonzero2 == 1 && nonzero_idx2[in_gr_idx] == 1){
                                    }else{

                                       rsd_tild2 = rsd_tild
                                       rsd_tild2[ ,out_gr[[m]][out_gr_idx]] = rsd_tild2[ ,out_gr[[m]][out_gr_idx]] + D[ ,j,drop=FALSE] %*% beta[j, out_gr[[m]][out_gr_idx], drop=FALSE]
                                       cond_3 = t(D[, j, drop=FALSE]) %*% rsd_tild2[ ,out_gr[[m]][out_gr_idx], drop=FALSE]

                                       if(abs(cond_3) <= lambda) {
                                          beta[j, out_gr[[m]][out_gr_idx]] = 0
                                          rsd_tild = rsd_tild2;
                                       }else{
                                          rel_in_gr = in_gr_dic[[j]]
                                          rel_out_gr = out_gr_dic[[out_gr[[m]][out_gr_idx]]]

                                          q2 = 0
                                          for(ridx in 1:length(rel_in_gr)){
                                             # The three dots mean line continuation. MATLAB expressions normally end at the end of the line unless they are specifically continued (exception: within the [] list building operator.)
                                             q2 = q2 + in_gr_w[rel_in_gr[ridx]]  * lambda2 / sqrt(sum( beta[in_gr[[rel_in_gr[ridx]]] ,out_gr[[m]][out_gr_idx]]^2));
                                          }

                                          w3 = 0;
                                          for(ridx in 1:length(rel_out_gr)){
                                             w3 = w3 + out_gr_w[rel_out_gr[ridx]] * lambda3 / sqrt(sum(beta[j, out_gr[[rel_out_gr[ridx]]]]^2))
                                          }

                                          ## This is for the eq.(14)
                                          beta_minus = min(0, 1/(1+q2+w3)*(cond_3 + lambda))
                                          beta_plus  = max(0, 1/(1+q2+w3)*(cond_3 - lambda))
                                          new_beta = beta_minus + beta_plus

                                          rsd_tild[ ,out_gr[[m]][out_gr_idx]] = rsd_tild2[ ,out_gr[[m]][out_gr_idx]] - D[ ,j, drop=FALSE] %*% new_beta
                                          beta[j, out_gr[[m]][out_gr_idx]] = new_beta

                                       } # This is for if(abs(cond_3) <= lambda){
                                    } # This is for if(num_nonzero2 == 1 && nonzero_idx2(in_gr_idx) == 1){
                                 } # This is for if(num_nonzero == 1 && nonzero_idx(out_gr_idx) == 1){
                              } # This is for if(zero_chk != 0){
                           } # Thisis for  for(out_gr_idx in 1:ogrp_size){
                        } # This is for 'if(cond_1 <= out_lambda^2){'
                     } # This is for  'if(denom != 0){'
                  } # This is for  'for(in_gr_idx in 1:length(in_gr[[m2]]))'
               } # This is for  'if(out_lambda > 0){'


            } # for if(zero_chk != 0)
         }  # <- for "for(m2 in 1:length(in_gr))"
      } # for for(m in 1:length(out_gr)){



      prev_obj = obj_value;
      obj_value = objFnSIO(D, Y, beta,  params.lambda1, params.lambda2, params.lambda3, params.lambda4, params.weight, params.in_gr_w
                           ,params.out_gr_w, params.in_gr, params.out_gr, params.dummy)

      if(verbose == 1){
         cat( paste( "Obj_prev", prev_obj, "Obj_New", obj_value, "\n" ));
      }

      if((abs(obj_value-prev_obj)/abs(prev_obj)) <= stopVal) break

   } #for outer_loop=1:10000
   return(beta)
} # for function end


































# standardize matrices
standardize_mat <- function(x){
   m_list = apply(x,2,mean)
   for(i in 1:dim(x)[2]) {
      x[ ,i] = x[ ,i] - m_list[i]
      s = sqrt(sum(x[ ,i]^2))
      if(s!=0) x[ ,i] = x[ ,i]/s
   }
   return(x)
}

objFnSIO <- function(dictionary, Y, beta, params.lambda1, params.lambda2, params.lambda3, params.lambda4, params.weight, params.in_gr_w
                     ,params.out_gr_w, params.in_gr, params.out_gr, params.dummy){
   # used for struct inout_lasso
   lambda1 = params.lambda1;
   lambda2 = params.lambda2;
   lambda3 = params.lambda3;
   lambda4 = params.lambda4;
   weight = params.weight;
   in_gr_w = params.in_gr_w;
   out_gr_w = params.out_gr_w;
   in_gr = params.in_gr;
   out_gr = params.out_gr;
   dummy = params.dummy; # the last index of individual snps

   obj_val = (Y - (dictionary %*% beta))^2;
   obj_val = sum(sum(obj_val))/2;

   reg1 = sum(sum(lambda1 * weight[1:dummy, ] * abs(beta[1:dummy, ])  ));
   if(dummy != 0){
      reg2 = sum(sum(lambda4 * weight[((dummy+1):dim(weight)[1]) , ] * abs(beta[((dummy+1):dim(weight)[1]), ])   ));
   }

   reg3 = 0
   for(t in 1:dim(Y)[2]){
      for(m in 1:length(in_gr)){
         reg3 = reg3 + lambda2 * in_gr_w[[m]] * sqrt(sum(beta[in_gr[[m]], t]^2))
      }
   }

   reg4 = 0;
   for(j in 1:dim(dictionary)[2]){
      for(m in 1:length(out_gr)){
         reg4 = reg4 + lambda3  * out_gr_w[[m]] * sqrt(sum(beta[j, out_gr[[m]]]^2));
      }
   }

   obj_val = obj_val + reg1 + reg2 + reg3 + reg4;
   return(obj_val)
}







Structured_dsgl <- function(A, L=rho, lambda1, lambda2, penalize.diagonal, Ind.Path.gene)
{
   lam1 = lambda1*1/L
   lam2 = lambda2*1/L

   if(is.matrix(A[[1]])){
      p=dim(A[[1]])[1]
   }

   K = length(A)
   softA = A
   for(k in 1:K) {
      softA[[k]] = soft(A[[k]],lam1,penalize.diagonal=penalize.diagonal)
   } #if penalize.diagonal=FALSE was used in ggl(), then this will not penalize the diagonal.
   normsoftA = A[[1]]*0

   for(k in 1:K) {
      normsoftA = normsoftA + (softA[[k]])^2
   } ### This involves a "Group-lasso" part.


   for(i in 1:length(Ind.Path.gene)){
      normsoftA[lower.tri(normsoftA, diag = TRUE)][Ind.Path.gene[[i]]] <- sum(normsoftA[lower.tri(normsoftA, diag = TRUE)][Ind.Path.gene[[i]]])
      tmp <- rep(NA, length(normsoftA))
      tmp[as.vector(lower.tri(A[[1]], diag = TRUE))] <- normsoftA[lower.tri(normsoftA, diag = TRUE)]
      tmpmat <- matrix(tmp, dim(normsoftA)[1], dim(normsoftA)[1])
      tmp[is.na(tmp)] <- t(tmpmat)[upper.tri(t(tmpmat))]
      tmpmat <- matrix(tmp, dim(A[[1]])[1], dim(A[[1]])[1])
   }

   normsoftA <- tmpmat
   normsoftA = sqrt(normsoftA)  ## Newly define the norm2 matrix

   notshrunk = (normsoftA>lam2)*1
   # reset 0 elements of normsoftA to 1 so we don't get NAs later.

   normsoftA = normsoftA + (1-notshrunk)

   out = A
   for(k in 1:K)
   {
      out[[k]] = softA[[k]]*(1-lam2/normsoftA)
      out[[k]] = out[[k]]*notshrunk
   }
   return(out)
}








plot.jgl <-
   function(x,...)
   {
      .env = "environment: namespace:JGL"
      #UseMethod("plot")
      theta=x$theta
      library(igraph)
      K=length(theta)
      adj = make.adj.matrix(theta)
      diag(adj)=0
      gadj = graph.adjacency(adj,mode="upper",weighted=TRUE)
      #weight the edges according to the classes they belong to
      E(gadj)$color = 2^(K)-get.edge.attribute(gadj,"weight")
      #plot the net using igraph
      plot(gadj, vertex.frame.color="white",layout=layout.fruchterman.reingold,
           vertex.label=NA, vertex.label.cex=3, vertex.size=1)
   }


print.jgl <-
   function(x, ...)
   {
      #	.env = "environment: namespace:JGL"
      #	UseMethod("print")
      p = dim(x$theta[[1]])[1]
      K = length(x$theta)
      #cat("Call: ")
      #dput(x$call)
      cat("\n")

      # return number of connected nodes:
      cat("Number of connected nodes: ",sum(x$connected),"\n")
      subnetlengths = nedges = c()
      for(k in 1:K)
      {
         subnetlengths[k] = length(subnetworks(x$theta)[[k]])
         nedges[k] = (sum(x$theta[[k]]!=0)-p)/2
      }
      # number of subgraphs:
      cat("Number of subnetworks in each class: ", subnetlengths, "\n")
      # number of edges per class:
      cat("Number of edges in each class: ", nedges, "\n")

      # number of shared edges:
      sharededges = matrix(0,p,p)
      for(k in 1:K){sharededges = sharededges + (x$theta[[k]]!=0)}
      nsharededges = (sum(sharededges==K)-p)/2
      cat("Number of edges shared by all classes: ", nsharededges, "\n")

      # L1 norm of classes:
      norm = c()
      for(k in 1:K) {norm[k] = sum(abs(x$theta[[k]]))-sum(abs(diag(x$theta[[k]])))}
      cat("L1 norm of off-diagonal elements of classes' Thetas: ", norm, "\n")
   }



admm.iters.unconnected = function(Y, lambda1, lambda2,penalty="fused",rho=1,rho.increment=1,weights,maxiter = 1000,tol=1e-5,warm=NULL)

   ## Does "lambda1" only have to do with diagonal element?? Check if it is true.

{
   K = length(Y)
   for(k in 1:K){Y[[k]]=as.matrix(Y[[k]])}
   p = dim(Y[[1]])[2]
   n=weights

   # define S vectors: empirical covs of each element
   S = list()
   for(k in 1:K)
   {
      ntemp=dim(Y[[k]])[1]
      S[[k]]=rep(0,p)
      for(j in 1:p)
      {
         if(var(Y[[k]][,j])*(ntemp-1)/ntemp!=0){
            S[[k]][j] = var(Y[[k]][,j])*(ntemp-1)/ntemp
         }else{
            S[[k]][j] <- 1e-10
         }

      }
   }

   # initialize theta:
   theta = list()
   for(k in 1:K){
      theta[[k]] = 1/S[[k]]
   }



   # initialize Z:
   Z = list(); for(k in 1:K){Z[[k]] = rep(0,p)}
   # initialize W:
   W = list();	for(k in 1:K){W[[k]] = rep(0,p)}

   # initialize lambdas:
   lam1 = lambda1
   lam2 = lambda2

   # iterations:
   iter=0
   diff_value = 10
   while((iter==0) || (iter<maxiter && diff_value > tol))
   {
      print(iter)
      # update theta to minimize -logdet(theta) + <S,theta> + rho/2||theta - Z + W ||^2_2:
      theta.prev = theta
      for(k in 1:K)
      {
         B = n[k]*S[[k]] - rho*(Z[[k]] - W[[k]])
         theta[[k]] = 1/(2*rho) * ( -B + sqrt(B^2 + 4*rho*n[k]) )
         #		B = S[[k]] - rho/n[k]*(Z[[k]] - W[[k]])
         #		theta[[k]] = n[k]/(2*rho) * ( -B + sqrt(B^2 + 4*rho/n[k]) )  # written as in our paper
      }

      # update Z:
      # define A matrices:
      A = list()
      for(k in 1:K){ A[[k]] = theta[[k]] + W[[k]] }
      if(penalty=="fused")
      {
         # use flsa to minimize rho/2 ||Z-A||_F^2 + P(Z):
         if(K==2){Z = flsa2(A,rho,lam1,lam2,penalize.diagonal=TRUE)}
         if(K>2){Z = flsa.general(A,rho,lam1,lam2,penalize.diagonal=TRUE)}
      }
      if(penalty=="group")
      {
         # minimize rho/2 ||Z-A||_F^2 + P(Z):
         Z = dsgl(A,rho,lam1,lam2,penalize.diagonal=TRUE)
      }

      # update the dual variable W:
      for(k in 1:K){W[[k]] = W[[k]] + (theta[[k]]-Z[[k]])}

      # bookkeeping:
      iter = iter+1
      diff_value = 0
      for(k in 1:K) {

         diff_value = diff_value + sum(abs(theta[[k]] - theta.prev[[k]])) / sum(abs(theta.prev[[k]]))

      }
      # increment rho by a constant factor:
      rho = rho * rho.increment
   }

   diff = 0; for(k in 1:K){diff = diff + sum(abs(theta[[k]]-Z[[k]]))}
   out = list(theta=theta,Z=Z,diff=diff,iters=iter)
   return(out)
}



flsa2 <-
   function(A,L,lam1,lam2,penalize.diagonal)  #A is a list of 2 matrices from which we apply an L2 penalty to departures
   {
      # 1st apply fused penalty:
      S1 = abs(A[[1]]-A[[2]])<=2*lam2/L
      X1 = (A[[1]]+A[[2]])/2
      Y1 = X1

      S2 = (A[[1]] > A[[2]]+2*lam2/L)
      X2 = A[[1]] - lam2/L
      Y2 = A[[2]] + lam2/L

      S3 = (A[[2]] > A[[1]]+2*lam2/L)
      X3 = A[[1]] + lam2/L
      Y3 = A[[2]] - lam2/L

      X = soft(a = S1*X1 + S2*X2 + S3*X3, lam = lam1/L, penalize.diagonal=penalize.diagonal)
      Y = soft(a = S1*Y1 + S2*Y2 + S3*Y3, lam = lam1/L, penalize.diagonal=penalize.diagonal)

      return(list(X,Y))
   }

soft <-
   function(a,lam,penalize.diagonal){ # if last argument is FALSE, soft-threshold a matrix but don't penalize the diagonal
      out <- sign(a)*pmax(0, abs(a)-lam)
      if(!penalize.diagonal) diag(out) <- diag(a)
      return(out)
   }



### ADMM for FGL when Theta is known to be diagonal (unconnected nodes):

#  call:	admm.iters.unconnected(Yu,lambda1=lam1.unconnected,lambda2=lam2.unconnected,rho,weights,
#   penalize.diagonal=FALSE,maxiter,tol)
# lam1.unconnected, lam2.u...  vectors of penalties.

admm.iters.unconnected = function(Y,lambda1,lambda2,penalty="fused",rho=1,rho.increment=1,weights,maxiter = 1000,tol=1e-5,warm=NULL)
{
   K = length(Y)
   for(k in 1:K){Y[[k]]=as.matrix(Y[[k]])}
   p = dim(Y[[1]])[2]
   n=weights

   # define S vectors: empirical covs of each element
   S = list()
   for(k in 1:K)
   {
      ntemp=dim(Y[[k]])[1]
      S[[k]]=rep(0,p)
      for(j in 1:p)
      {
         if(var(Y[[k]][,j])*(ntemp-1)/ntemp!=0){
            S[[k]][j] = var(Y[[k]][,j])*(ntemp-1)/ntemp
         }else{
            S[[k]][j] <- 1e-10
         }
      }
   }

   # initialize theta:
   theta = list()
   for(k in 1:K){
      theta[[k]] = 1/S[[k]]
   }

   # initialize Z:
   Z = list(); for(k in 1:K){Z[[k]] = rep(0,p)}
   # initialize W:
   W = list();	for(k in 1:K){W[[k]] = rep(0,p)}

   # initialize lambdas:
   lam1 = lambda1
   lam2 = lambda2

   # iterations:
   iter=0
   diff_value = 10
   while((iter==0) || (iter<maxiter && diff_value > tol))
   {
      print(iter)
      # update theta to minimize -logdet(theta) + <S,theta> + rho/2||theta - Z + W ||^2_2:
      theta.prev = theta
      for(k in 1:K)
      {
         B = n[k]*S[[k]] - rho*(Z[[k]] - W[[k]])
         theta[[k]] = 1/(2*rho) * ( -B + sqrt(B^2 + 4*rho*n[k]) )
         #		B = S[[k]] - rho/n[k]*(Z[[k]] - W[[k]])
         #		theta[[k]] = n[k]/(2*rho) * ( -B + sqrt(B^2 + 4*rho/n[k]) )  # written as in our paper
      }

      # update Z:
      # define A matrices:
      A = list()
      for(k in 1:K){ A[[k]] = theta[[k]] + W[[k]] }
      if(penalty=="fused")
      {
         # use flsa to minimize rho/2 ||Z-A||_F^2 + P(Z):
         if(K==2){Z = flsa2(A,rho,lam1,lam2,penalize.diagonal=TRUE)}
         if(K>2){Z = flsa.general(A,rho,lam1,lam2,penalize.diagonal=TRUE)}
      }
      if(penalty=="group")
      {
         # minimize rho/2 ||Z-A||_F^2 + P(Z):
         Z = dsgl(A,rho,lam1,lam2,penalize.diagonal=TRUE)
      }

      # update the dual variable W:
      for(k in 1:K){W[[k]] = W[[k]] + (theta[[k]]-Z[[k]])}

      # bookkeeping:
      iter = iter+1
      diff_value = 0
      for(k in 1:K) {
         diff_value = diff_value + sum(abs(theta[[k]] - theta.prev[[k]])) / sum(abs(theta.prev[[k]]))
      }
      # increment rho by a constant factor:
      rho = rho * rho.increment
   }
   diff = 0; for(k in 1:K){diff = diff + sum(abs(theta[[k]]-Z[[k]]))}
   out = list(theta=theta,Z=Z,diff=diff,iters=iter)
   return(out)
}

penalty.as.matrix <-
   function(lambda, p, penalize.diagonal)
   {
      # for matrix penalties:  check dim and symmetry:
      if(is.matrix(lambda))
      {
         if(sum(lambda!= t(lambda))>0) {stop("error: penalty matrix is not symmetric")}
         if(sum(abs(dim(lambda)-p))!=0 ) {stop("error: penalty matrix has wrong dimension")}
      }
      # for scalar penalties: convert to matrix form:
      if(length(lambda)==1) {lambda=matrix(lambda,p,p)}
      # apply the penalize.diagonal argument:
      if(!penalize.diagonal) {diag(lambda)=0}
      return(lambda)
   }




### ADMM for FGL:
#admm.iters = function(Y,lambda1,lambda2,penalty="fused",rho=1,rho.increment=1,weights,penalize.diagonal,maxiter = 1000,tol=1e-5,warm=NULL){
#  K = length(Y)
#  p = dim(Y[[1]])[2]
#  n=weights

#  ns = c(); for(k in 1:K){ns[k] = dim(Y[[k]])[1]}
#  S = list(); for(k in 1:K){S[[k]] = cov(Y[[k]])*(ns[k]-1)/ns[k]}

# initialize theta:
#theta = list()
#  for(k in 1:K){theta[[k]] = diag(1/diag(S[[k]]))}
# initialize Z:
#  Z = list(); for(k in 1:K){Z[[k]]=matrix(0,p,p)}
#  # initialize W:
#  W = list();	for(k in 1:K) {W[[k]] = matrix(0,p,p) }

# initialize lambdas:  (shouldn't need to do this if the function is called from the main wrapper function, JGL)
#  lam1 = penalty.as.matrix(lambda1,p,penalize.diagonal=penalize.diagonal)
#  if(penalty=="fused") {lam2 = penalty.as.matrix(lambda2,p,penalize.diagonal=TRUE)}
#  if(penalty=="group") {lam2 = penalty.as.matrix(lambda2,p,penalize.diagonal=penalize.diagonal)}
# iterations:
#  iter=0
#  diff_value = 10
#  while((iter==0) || (iter<maxiter && diff_value > tol))
#  {
# reporting
#	if(iter%%10==0)
#    if(FALSE)
#    {
#      print(paste("iter=",iter))
#      if(penalty=="fused")
#      {
#        print(paste("crit=",crit(theta,S,n=rep(1,K),lam1,lam2,penalize.diagonal=penalize.diagonal)))
#        print(paste("crit=",crit(Z,S,n=rep(1,K),lam1,lam2,penalize.diagonal=penalize.diagonal)))
#      }
#      if(penalty=="group"){print(paste("crit=",gcrit(theta,S,n=rep(1,K),lam1,lam2,penalize.diagonal=penalize.diagonal)))}
#    }

# update theta:
#   theta.prev = theta
#    for(k in 1:K){
#      edecomp = eigen(S[[k]] - rho*Z[[k]]/n[k] + rho*W[[k]]/n[k])
#      D = edecomp$values
#      V = edecomp$vectors
#      D2 = n[k]/(2*rho) * ( -D + sqrt(D^2 + 4*rho/n[k]) )
#      theta[[k]] = V %*% diag(D2) %*% t(V)
#    }

# update Z:
# define A matrices:
#    A = list()
#    for(k in 1:K){ A[[k]] = theta[[k]] + W[[k]] }
#    if(penalty=="fused")
#    {
# use flsa to minimize rho/2 ||Z-A||_F^2 + P(Z):
#      if(K==2){Z = flsa2(A,rho,lam1,lam2,penalize.diagonal=TRUE)}
#      if(K>2){Z = flsa.general(A,rho,lam1,lam2,penalize.diagonal=TRUE)}  # the option to not penalize the diagonal is exercised when we initialize the lambda matrices
#    }
#    if(penalty=="group")
#    {
#  minimize rho/2 ||Z-A||_F^2 + P(Z):
#      Z = dsgl(A, rho, lam1, lam2, penalize.diagonal=TRUE)
#    }

# update the dual variable W:
#    for(k in 1:K){W[[k]] = W[[k]] + (theta[[k]]-Z[[k]])}

# bookkeeping:
#    iter = iter+1
#    diff_value = 0
#    for(k in 1:K) {diff_value = diff_value + sum(abs(theta[[k]] - theta.prev[[k]])) / sum(abs(theta.prev[[k]]))}
# increment rho by a constant factor:
#    rho = rho*rho.increment
#  }
#  diff = 0; for(k in 1:K){diff = diff + sum(abs(theta[[k]]-Z[[k]]))}
#  out = list(theta=theta,Z=Z,diff=diff,iters=iter)
#  return(out)
#}




dsgl <- function(A,L,lam1,lam2,penalize.diagonal)
{
   lam1 = lam1*1/L
   lam2 = lam2*1/L

   if(is.matrix(A[[1]])) {p=dim(A[[1]])[1]}
   if(is.vector(A[[1]])) {p=length(A[[1]])}
   K=length(A)
   softA = A
   for(k in 1:K) {softA[[k]] = soft(A[[k]],lam1,penalize.diagonal=penalize.diagonal) }   #if penalize.diagonal=FALSE was used in ggl(), then this will not penalize the diagonal.
   normsoftA = A[[1]]*0
   for(k in 1:K) {normsoftA = normsoftA + (softA[[k]])^2}

   normsoftA = sqrt(normsoftA)

   notshrunk = (normsoftA>lam2)*1
   # reset 0 elements of normsoftA to 1 so we don't get NAs later.
   normsoftA = normsoftA + (1-notshrunk)

   out = A
   for(k in 1:K)
   {
      out[[k]] = softA[[k]]*(1-lam2/normsoftA)
      out[[k]] = out[[k]]*notshrunk
   }
   return(out)
}



flsa.general <-
   function(A,L,lam1,lam2,penalize.diagonal)
   {
      trueA = A
      if(is.matrix(A[[1]])) {p=dim(A[[1]])[1]}
      if(is.vector(A[[1]])) {p=length(A[[1]])}
      K = length(A)
      # results matrices:
      X = list()
      #for(k in 1:K) {X[[k]] = matrix(NA,p,p)}
      for(k in 1:K) {X[[k]] = A[[1]]*NA}
      if(is.matrix(A[[1]])) {fusions = array(FALSE,dim=c(K,K,p,p))}
      if(is.vector(A[[1]])) {fusions = array(FALSE,dim=c(K,K,p,1))}

      # get starting newc: list of matrices.  newc[[k]][i,j] gives the (2*ordermats[[k]]-K-1) adjustment for lam2 at that element.  Basically, how many rank higher vs how many rank lower?
      newc = list()
      for(k in 1:K)
      {
         others = setdiff(1:K,k)
         others.smaller.k = 1:(k-1)
         newc[[k]] = A[[1]]*0
         for(o in others) {newc[[k]] = newc[[k]] + (A[[o]]-A[[k]]< -1e-4) - (A[[o]]-A[[k]]>1e-4)}
      }

      ######### start the loop here:
      for(iter in 1:(K-1))
      {

         # create order matrices:
         ordermats = list()
         for(k in 1:K)
         {
            others = setdiff(1:K,k)
            others.smaller.k = 1:(k-1)
            ordermats[[k]] = A[[1]]*0
            for(o in others) {ordermats[[k]] = ordermats[[k]] + (A[[k]]-A[[o]]>1e-4)}
            # to deal with ties, also add a unit to ordermat[[k]] if a class with a lower k has a tie at element i,j:
            if(k>1)
            {
               for(o in others.smaller.k) {ordermats[[k]] = ordermats[[k]] + (abs(A[[o]]-A[[k]])<1e-4)}
            }
            ordermats[[k]] = ordermats[[k]] + 1
         }

         # create beta.g matrices, holding the solution to Holger's "unconstrained problem"
         #  (prending we're not constraining the order of the solution to match the order of the A matrices)
         betas.g = list()
         for(k in 1:K)
         {
            betas.g[[k]] = A[[k]] - lam2/L*newc[[k]]
         }

         # identify and fuse all elements for which the betas.g are out of order:
         new.ordermats = list()
         for(k in 1:K)
         {
            others = setdiff(1:K,k)
            others.smaller.k = 1:(k-1)
            new.ordermats[[k]] = A[[1]]*0
            for(o in others) {new.ordermats[[k]] = new.ordermats[[k]] + (betas.g[[k]]-betas.g[[o]]>1e-4)}
            # to deal with ties, also add a unit to ordermat[[k]] if a class with a lower k has a tie at element i,j:
            if(k>1)
            {
               for(o in others.smaller.k) {new.ordermats[[k]] = new.ordermats[[k]] + (abs(betas.g[[o]]-betas.g[[k]])<1e-4)}
            }
            new.ordermats[[k]] = new.ordermats[[k]] + 1
         }

         # identify neighboring fusions:  "fusions": K x K x p x p array: K x K matrices, T/F for fusions
         for(k in 1:K){
            for(kp in 1:K){
               #given k,kp, declare a fusion when their ordermats entries are adjacent, and their new.ordermats entries have the opposite direction:
               fusions[k,kp,,] = fusions[k,kp,,]+
                  ((ordermats[[k]]-1==ordermats[[kp]])&(new.ordermats[[k]]<new.ordermats[[kp]]))+
                  ((ordermats[[k]]+1==ordermats[[kp]])&(new.ordermats[[k]]>new.ordermats[[kp]]))+
                  (abs(A[[k]]-A[[kp]])<1e-4)
               #(existing fusions, neighboring fusions, and ties)
               fusions = (fusions>0)*1
            }}


         # now we've noted fusions between all entries with adjacent ordermats entries and reversed new.ordermats entries
         # next: extend fusions to non-adjecent entries: if a-b and b-c, then connect a-c:
         for(k in 1:K){
            for(kp in 1:K){
               others = setdiff(1:K,c(k,kp))
               for(o in others)
               {
                  #identify elements in o which are fused with the same element in both k and kp, then add them to the list of k-kp fusions:
                  bothfused = fusions[k,o,,] & fusions[kp,o,,]
                  fusions[k,kp,,] = fusions[k,kp,,] | bothfused
               }
            }}

         # now recalculate A with the new fused entries:
         # to recalculate A, for each non-zero entry, identify the classes k with which it must be fused, and get their average:
         for(k in 1:K)
         {
            others = setdiff(1:K,k)
            #fusemean and denom: the mean value of all the trueA to be fused, and the number of values to be fused:
            fusemean = trueA[[k]]
            denom = A[[1]]*0+1
            for(o in others)
            {
               fusemean = fusemean+fusions[k,o,,]*trueA[[o]]  #add the values of the elements which must be fused to fusemean
               denom = denom+fusions[k,o,,]
            }
            # now redefine A[[k]]: unchanged from trueA if there's no fusion, and the mean of the fused elements when there is fusion:
            A[[k]] = fusemean/denom
         }

         #newc: list of matrices.  newc[[k]][i,j] gives the (2*ordermats[[k]]-K-1) adjustment for lam2 at that element.  Basically, how many rank higher vs how many rank lower?
         newc = list()
         for(k in 1:K)
         {
            others = setdiff(1:K,k)
            others.smaller.k = 1:(k-1)
            newc[[k]] = A[[1]]*0
            for(o in others) {newc[[k]] = newc[[k]] + (A[[o]]-A[[k]]< -1e-4) - (A[[o]]-A[[k]]>1e-4)}
         }

      } #end loop here

      # final version of betas.g:
      for(k in 1:K)
      {
         betas.g[[k]] = A[[k]] - lam2/L*newc[[k]]
      }
      # now soft-threshold the solution matices:
      for(k in 1:K)
      {
         X[[k]] = soft(betas.g[[k]],lam=lam1/L,penalize.diagonal=penalize.diagonal)
      }
      return(X)
   }

