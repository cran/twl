

just_choose <- function(x){
    if(!any(na.omit(x>0))){
	pr <- 1/length(x)
	warning("some unnormalized probabilities should be positive but are not")
    }
    rmultinom(1,1,x)
}

likeli_manip <- function(x,tune=0.02){
    tmp1 <- exp(x -max(x)+ 15)
    tmp1 <- (tmp1/sum(tmp1))+tune #+0.025
    return(tmp1)
}

ring <- function(x,start){
    if(start !=1) return(c(tail(x,-(start-1)),head(x,(start-1))))
    else return(x)
}

inf_func_mult <- function(x){
    pr <- exp(x)
    if(!any(na.omit(x>0))){
        pr <- 1/length(x)
	warning("some unnormalized probabilities should be positive but are not")
    }
    rmultinom(1,1,pr)
}

my_lognorm <- function(x,me,var){
    -(x-me)^2/(2*var) - log(2*pi*var)/2
}

mygroupcolVars <- function (x, ina, std = FALSE){
    m <- rowsum(x, ina)
    m2 <- rowsum(x^2, ina)
    ni <- c(table(ina))
    s <- (m2 - m^2/ni)/(ni - 1)
    if (std) 
        s <- sqrt(s)
    s
}


#' Main function to obtain posterior samples from a TWL model.
#' 
#' @param full_dat_mat list of matrices of the different data types.
#' @param full_dat list of data.tables with a single column labelled 'nam', denoting sample annotation.  A consistent naming convention of samples must be used across data types.
#' @param alpha_re Hyperparameter for the dirichlet prior model within each data type, influencing sparsity of clusterings.  A smaller number encourages fewer clusters.    Defaults to 7 and should be chosen as a function of sample size.
#' @param beta_re Hyperparameter for the dirichlet prior model across datatypes within each sample, influencing the degree to which  each data type's sample cluster labels affect those of the other data types.  Defaults to 0.4 and should be chosen as a function of the total number of data types being integrated in the analysis.
#' @param num_its Number of iterations.  Defaults to 5000.
#' @param num_all_clus Ceiling on the number of clusters.  Defaults to 30.  Should be chosen as some factor greater (for example, 5), than maximum number of hypothesized clusters in the data types.
#' @param output_every Frequency of sampling log statistics, reporting mixing, cluster distribution, and proportion of cluster sharing across data types.  Defaults to once every 20 iterations.  
#' @param manip TRUE/FALSE for whether likelihood manipulation should be used to increase mixing in situations where cluster means are far from one another in Euclidean distance.  This should not influence identified clusters nor parameters associated with them.  Defaults to TRUE.  
#' @param sav_inter A logical indicating whether a temporary file of the samples should be written out in the working directory every 50 iterations.  Allows for restarts when sampling is interrupted, and defaults to FALSE.
#' @return A list of lists of data.tables.  The list length is the number of iterations.  The length of each element is the number of data types.  The data.tables have 2 columns, sample annotation called `nam' and cluster assignment called 'clus'.
#' @examples
#' data(data_and_output)
#' \dontrun{clus_save <- TWLsample(misaligned_mat,misaligned,output_every=50,num_its=5000,manip=FALSE)
#' outpu_new <- pairwise_clus(clus_save,BURNIN=2000)
#' }
#' post_analy_cor(outpu_new,c("title1","title2","title3","title4","title5"),
#' tempfile(),ords='none') 
#' clus_labs <- post_analy_clus(outpu_new,clus_save,c(2:6),rep(0.6,5),c("title1","title2",
#' "title3","title4","title5"),tempfile())
#' output_nest <- cross_dat_analy(clus_save,4900)
#' @export
TWLsample <- function(full_dat_mat,full_dat,alpha_re=7,beta_re=0.4,num_its=5000,num_all_clus=30,output_every=20,manip=TRUE,sav_inter=FALSE){
    alpha <- alpha_re -1 # translated for use in code
    beta <- beta_re - 1 # translated for use in code
    t_full_dat_mat <- list()
    
    for(hh in seq_along(full_dat_mat)){
        t_full_dat_mat[[hh]] <- t(full_dat_mat[[hh]])
    }

    ALL_CLUS <- seq(num_all_clus)
    num_datasets <- length(full_dat_mat)

    num_row <- list()
    num_col <- list()

    def_lik <- list()
    likelihood_mat_pre <- list()
    def_lik_diag_var <- list()
    
    me_lis_prior <- list()
    cov_prior <- list()

    for(m in seq_along(full_dat_mat)){
        d_tm <- dim(full_dat_mat[[m]])
        num_row[[m]] <- d_tm[1]
        num_col[[m]] <- d_tm[2]

        full_dat[[m]][,clus:= sample(ALL_CLUS,num_row[[m]],replace=TRUE)]

        def_lik_diag_var[[m]] <- unlist(Rfast::colVars(full_dat_mat[[m]])) ## ASSUMES SAME SCALING OF VARIABLES WITHIN DATA SET -- THINK SHOULD BE COLUMN SPECIFIC

        def_lik[[m]] <- (log(1/sqrt(2*pi*def_lik_diag_var[[m]]))-3/4)
        
        # COMPLIANCE CHANGE
        likelihood_mat_pre[[m]] <- matrix(sum(def_lik[[m]]),nrow = num_row[[m]],ncol=num_all_clus,byrow=T)
        me_lis_prior[[m]] <-  colMeans(full_dat_mat[[m]])

        cov_prior[[m]] <- Rfast::colVars(full_dat_mat[[m]])*150/num_row[[m]] # 150 comes from simplications of factors found in the paper
    }

    ## WE GENERATE THIS ALTERNATIVE full_dat SO THAT THE MATCHING STILL WORKS
    match_lis <- list()
    for(bb in 1:num_datasets){

        ref <- full_dat[[bb]][,nam]
        tmp_dat_mat <- data.table::rbindlist(lapply(full_dat,function(x){x[nam%in%ref,list(clus,nam)]})) #RCPP?

        tmp_dat_mat_clus <- tmp_dat_mat[,table(c(clus,ALL_CLUS)),by=nam]
        mat_row_nams <- rle(tmp_dat_mat_clus[[1]])$values

        match_lis[[bb]] <- match(ref,mat_row_nams)
    }

    clus_save <- list()
    likelihood_mat <- list()

    for(j in 1:num_its){
        for(ii in seq_along(full_dat_mat)){
            count_na_cov <- 0
            likeli_lis <- list()
            freqs <- as.numeric(table(c(full_dat[[ii]][,clus],ALL_CLUS)))

            temp1 <- 1 # 5
            temp2 <- 1 # 0.7

            sparse_row <- 15 
            sparse_col <- 15 

            tab_clus_tmp <- table(full_dat[[ii]][,clus])
            uniq_clus <- unique(full_dat[[ii]][,clus])
            num_clus <- length(uniq_clus)
            tab_clus_tmp2 <- c(tab_clus_tmp)
            tab_clus_tmp3 <- c(tab_clus_tmp[match(uniq_clus,sort(uniq_clus))])

            reord_cov <- match(uniq_clus,sort(uniq_clus))
            cov_lis_glob_orig <- mygroupcolVars(full_dat_mat[[ii]],full_dat[[ii]][,clus])[reord_cov,]
            tmp_NA <- sum(is.nan(cov_lis_glob_orig[,1]))
            
            nan_inds <- which(is.nan(cov_lis_glob_orig),arr.ind=T)
            cov_lis_glob_orig[nan_inds] <- def_lik_diag_var[[ii]]
            arr_inds <- which(cov_lis_glob_orig==0,arr.ind=T)
            cov_lis_glob_orig[which(cov_lis_glob_orig==0)] <- def_lik_diag_var[[ii]][arr_inds[,2]]

            cov_lis_glob <-colMeans(cov_lis_glob_orig,na.rm=TRUE) #cov_lis <- 1/(1/(tab_clus_tmp2 + 1)* 1/def_lik_diag_var[[ii]] + tab_clus_tmp2/(tab_clus_tmp2+1) * 1/ cov_lis_glob)

            cov_lis <- 1/(outer(1/(tab_clus_tmp3 + 1),1/def_lik_diag_var[[ii]]) + outer(tab_clus_tmp3/(tab_clus_tmp3+1), 1/cov_lis_glob))
            me_lis_pre <- (rowsum(full_dat_mat[[ii]],full_dat[[ii]][,clus])/tab_clus_tmp2)[reord_cov,]

            me_lis_glob <- colMeans(me_lis_pre)

            me_lis_former <- (outer(1/(tab_clus_tmp3 + 11),me_lis_glob) + (tab_clus_tmp3+10)/(tab_clus_tmp3+11)*me_lis_pre )
            cov_lis_glob_me <- colSums(cov_lis_glob_orig*tab_clus_tmp3)/num_row[[ii]] # assumes column major

            w_l <- outer(tab_clus_tmp3,1/cov_lis_glob_me)
            precis_denom <- rep(1/cov_prior[[ii]],each=num_clus) + w_l

            me_lis_shrunk <- (w_l*me_lis_pre+ rep(me_lis_prior[[ii]]/cov_prior[[ii]],each=num_clus))/(precis_denom)
            me_lis <- matrix(rnorm(length(me_lis_shrunk), me_lis_shrunk, sqrt(1/precis_denom) ), nrow=num_clus,ncol=num_col[[ii]])

            likelihood_mat[[ii]] <- likelihood_mat_pre[[ii]]
            num_div <- num_col[[ii]]%%100

            for(gg in 1:num_clus){
                likeli_lis[[gg]] <- colSums(my_lognorm(t_full_dat_mat[[ii]], me_lis[gg,],cov_lis_glob_me))
                likelihood_mat[[ii]][,uniq_clus[gg]] <- likeli_lis[[gg]]
            }

            curr_nams <- full_dat[[ii]][,nam]
            tmp_dat_mat <- rbindlist(lapply(full_dat,function(x){x[nam%in%curr_nams,list(clus,nam)]})) #RCPP?
            
            tmp_dat_mat_clus <- tmp_dat_mat[,(table(c(ALL_CLUS,clus)) + beta),by=nam]
            mat_row_sparse <- matrix(tmp_dat_mat_clus[[2]],nrow=num_all_clus,ncol=num_row[[ii]])[,match_lis[[ii]]]

            mat_row_sparse_prob <- apply(mat_row_sparse,2,function(x){MCMCpack::rdirichlet(1,sparse_row * (x))})
            
            freqs_prob <- MCMCpack::rdirichlet(1,sparse_col*(freqs+alpha))[1,] ## 1 is the prior here, 0 plus the 'convenience' vector ALL_CLUS above, which makes table always give the same format # the 6 here can be smaller, as, say, 2

            if(!manip){
                probs_big_mat <- t(likelihood_mat[[ii]]) + temp1 * log(mat_row_sparse_prob) + temp2 * log(freqs_prob)
                probs_big_norm <- apply(probs_big_mat,2,function(x){x -max(x)+ 15}) # transpose to counter implicit transpose,  15 chosen according to what can exp sum to infinity when there are 40 clusters.  Would need all to have approx the same prob
                clus_samp <- apply(probs_big_norm,2,function(x){which(inf_func_mult(x)==1)}) # RCPP?

            } else {
                probs_part_mat <- exp(temp1 * log(freqs_prob) + temp2 * log(mat_row_sparse_prob))
                likeli_trans <- apply(t(likelihood_mat[[ii]]),2,likeli_manip)
                probs_big_norm <- likeli_trans*probs_part_mat
                clus_samp <- apply(probs_big_norm,2,function(x){which(just_choose(x)==1)})

            }
            full_dat[[ii]][,clus:=clus_samp]
        }
            clus_save[[j]] <- data.table::copy(full_dat)
            
            if ((j %% 500) == 0){
	       if(sav_inter){
	                       save(clus_save,file=paste0(getwd(),"/TWLclus_tmp_chain.RObject"),compress=F)
			       }
            }
            if ((j %% output_every) == 0){
                print(paste("ITERATION: ",j))
                print('')
                print('')
                for(bbb in 1:num_datasets){
                    print(paste0("Clust dist for ",bbb,"th dataset"))
                    print(((full_dat[[bbb]][,table(c(clus,ALL_CLUS))-1])))
                    print("")
                }
                for(bbb in 1:num_datasets){
                    print(paste0('proportion samples in same clus for ',bbb,"th dataset"))
                    print(mean(clus_save[[j-1]][[bbb]][,clus] == clus_save[[j]][[bbb]][,clus]))
                    print('')
                }
                print('')

                if(j>20){
                    for(vv in 1:(num_datasets-1)){
                        for(vvv in (vv+1):num_datasets){
                            print(paste0("average cluster overlap between ",vv,"th and ",vvv,"th datasets"))
                            print(mean(sapply(clus_save[(j-19):j],function(x){mean(merge(x[[vv]],x[[vvv]],'nam')[,clus.x==clus.y])})))
                        }
                    }
                }
            }
    }
    return(clus_save)
}

clus <- clus.x <- clus.y <- j <- nam <- NULL 
