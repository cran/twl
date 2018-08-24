
make_merge_by <- function(string){
    tmp <- function(x,y){
        data.table::setnames(x,names(x)[2],paste0("clus",abs(rnorm(1))))
        data.table::setnames(y,names(y)[2],paste0("clus",abs(rnorm(1))))
        merge(x,y,by=string)
    }
    return(tmp)
}

merge_on_nam <- make_merge_by("nam")

equal_without_nam <- function(dat_tab){
    dat_tab[,.SD,.SDcols = names(dat_tab) %like% 'clus'][,apply(.SD,1,function(x){all(x[1]==x[-1])})]
}

#' Create posterior similarity matrix from outputted list of clustering samples
#'
#' @param clus_save list of samples outputted from TWLsample function.
#' @param BURNIN number of samples devoted to burn-in.  Defaults to 2000.
#' @return outpu a list whose length is the number of datasets being integrated, and each elemnt of which is a posterior similarity matrix.  The dimension of each symmetric matrix is the number of samples in the respective dataset, and elements in the matrix are values between 0 and 1, and estimate of the probability 2 samples find themselves in the same clustering.
#' @examples
#' data(data_and_output)
#' \dontrun{clus_save <- TWLsample(misaligned_mat,misaligned,output_every=50,num_its=5000,manip=FALSE)
#' outpu_new <- pairwise_clus(clus_save,BURNIN=2000)
#' post_analy_cor(outpu_new,c("title1","title2","title3","title4","title5"),
#' tempfile(),ords='none') 
#' clus_labs <- post_analy_clus(outpu_new,clus_save,c(2:6),rep(0.6,5),c("title1","title2",
#' "title3","title4","title5"),tempfile())
#' output_nest <- cross_dat_analy(clus_save,4750)
#' }
#' @export
pairwise_clus <- function(clus_save,BURNIN=2000){
    outpu <- list()
    num_datasets <- length(clus_save[[1]])
    num_row <- sapply(clus_save[[1]],function(x){dim(x)[1]})
    num_its <- length(clus_save)
    for(whic in 1:num_datasets){
        k <- 0
        comp_lis <- list()
        for(ee in (BURNIN+1):num_its){
            k <- k+1
            me <- matrix(99,num_row[whic],num_row[whic])
            for (v in 1:(num_row[whic])){
                me[,v] <- clus_save[[ee]][[whic]][,clus[v]==clus]
            }
            comp_lis[[k]] <- me
        }
        outpu[[whic]] <- Reduce('+', comp_lis)/(num_its-BURNIN + 1)
        print(whic)
    }
    return(outpu)
}


#' Creates and saves correlation plots based on posterior similarity matrices
#'
#' @param outpu_new the output of the pairwise_clus function, and a list whose length is the number of datasets being integrated, and each elemnt of which is a posterior similarity matrix.  The dimension of each symmetric matrix is the number of samples in the respective dataset, and elements in the matrix are values between 0 and 1, and estimate of the probability 2 samples find themselves in the same clustering.
#' @param titles a vector of strings of length number of integrated datasets.  Elements of the vector are titles in the respective correlation plots
#' @param pdf_path file path where the plots will be saved as a pdf.
#' @param ords whether the correlation plots should be reordered according to that of hierarchical clustering for a more comprehensible plot.  Defaults to 'none'.  Passing any string apart from 'none' (i.e., 'yes') will result in the re-ordering.
#' @return dendro_ord regardless of whether correlation plots are reordered according to hierarchical clustering, a list of reorderings is returned of length the number of datasets on which analysis was performed.
#' @examples
#' data(data_and_output)
#' \dontrun{clus_save <- TWLsample(misaligned_mat,misaligned,output_every=50,num_its=5000,manip=FALSE)
#' outpu_new <- pairwise_clus(clus_save,BURNIN=2000)
#' post_analy_cor(outpu_new,c("title1","title2","title3","title4","title5"),
#' tempfile(),ords='none') 
#' clus_labs <- post_analy_clus(outpu_new,clus_save,c(2:6),rep(0.6,5),c("title1","title2",
#' "title3","title4","title5"),tempfile())
#' output_nest <- cross_dat_analy(clus_save,4750)
#' }
#' @export
post_analy_cor <- function(outpu_new,titles,pdf_path,ords='none'){
    dendro_ord <- list()
    for(b in seq_along(outpu_new)){
        if(ords!='none'){
            tmp <- hclust(as.dist(1-outpu_new[[b]]))
            ords1 <- tmp[[3]]
        } else {
            ords1 <- seq(dim(outpu_new[[b]])[1])
        }
        dendro_ord[[b]] <- ords1
    }

    pdf(pdf_path)
    for(b in seq_along(outpu_new)){
        corrplot::corrplot(outpu_new[[b]][dendro_ord[[b]],dendro_ord[[b]]],order='original',addgrid=NA,tl.pos='n',mar=rep(2,4))
    }
    dev.off()
    return(dendro_ord)
}




#' Assigns cluster labels by building dendrogram and thresholding at specified height
#'
#' @param outpu_new the output of the pairwise_clus function, and a list whose length is the number of datasets being integrated, and each elemnt of which is a posterior similarity matrix.  The dimension of each symmetric matrix is the number of samples in the respective dataset, and elements in the matrix are values between 0 and 1, and estimate of the probability 2 samples find themselves in the same clustering.
#' @param clus_sav_new list of samples outputted from TWLsample function.  See details for additional explanation of this parameter and height_clusts_vec.
#' @param num_clusts a vector of length the number of integrated datasets, specifying the number of cluster labels to be identified from the generated dendrogram for each dataset
#' @param height_clusts_vec vector of dendrogram heights of length the number of integrated datasets (if the analyst prefers manual inspection of outputted dendrograms and specification of the heights at which to threshold, thereby defining cluster membership).  Defaults to NULL.  See details for additional explanation of this parameter and num_clusts.
#' @param titles Vector of strings of length the number of datasets, used as prefixes in column labels of the outputted list of data.tables.
#' @param pdf_path  file path where the dendrogram figures will be saved as a pdf.
#' @return post_lab a list of data.tables of 2 columns each with names 'nam' and '*_clus', the nam specifying sample name annotation, and *_clus with the assigned cluster, where * is the corresponding element in the title argument vector.  
#' @examples
#' data(data_and_output)
#' \dontrun{clus_save <- TWLsample(misaligned_mat,misaligned,output_every=50,num_its=5000,manip=FALSE)
#' outpu_new <- pairwise_clus(clus_save,BURNIN=2000)
#' post_analy_cor(outpu_new,c("title1","title2","title3","title4","title5"),
#' tempfile(),ords='none') 
#' clus_labs <- post_analy_clus(outpu_new,clus_save,c(2:6),rep(0.6,5),c("title1","title2",
#' "title3","title4","title5"),tempfile())
#' output_nest <- cross_dat_analy(clus_save,4750)
#' }
#' @details At least one of either num_clusts or height_clusts_vec, or both, can be specified.  If both are specified, then heights is first used within the dendrogram for preliminary cluster assignment, then the X largest clusters of these receive final, outputted, assignment (the rest receiving a "clus_unknown" label), where X is the corresponding element in the num_clusts argument vector.
#' @export
post_analy_clus <- function(outpu_new,clus_sav_new,num_clusts,height_clusts_vec=NULL,titles,pdf_path){      ### EXAMINE THESE GRAPHS THEN MANUALLY SPECIFY THE NUMBER OF CLUSTERS IN HERE
    pdf(pdf_path)
    post_lab <- list()
    for(b in seq_along(outpu_new)){
        par(cex.lab=1.6)
        hclus1 <- hclust(as.dist(1-(outpu_new[[b]])))
        plot(hclus1,cex.lab=1.6,main="",labels=FALSE,xlab="Sample")

        if(is.null(height_clusts_vec[b])){
            hclus1_rec <- rect.hclust(hclus1,k=num_clusts[b])
        }
        if(!is.null(height_clusts_vec[b])){
            hclus1_rec <- rect.hclust(hclus1,h=height_clusts_vec[b])
            if(!is.null(num_clusts)){
                le <- sapply(hclus1_rec,length)
                thresh <- sort(le,decreasing=TRUE)[num_clusts[b]]
                hclus1_rec[which(le<thresh)] <- NULL
            }
        }
        post_lab[[b]] <- matrix(rep('clus_unknown',2*(dim(outpu_new[[b]])[1])),ncol=2)
        hclus1_rec <- hclus1_rec[!sapply(hclus1_rec,is.null)]

        for(kk in 1:length(hclus1_rec)){
            post_lab[[b]][hclus1_rec[[kk]],2] <- paste0('clus_',kk)
        }

        post_lab[[b]][,1] <- clus_sav_new[[1]][[b]][,nam]
        post_lab[[b]] <-  data.table::as.data.table(post_lab[[b]])
        names(post_lab[[b]]) <- c('nam',paste0(titles[b],'_clus'))
    }
    dev.off()
    return(post_lab)
}



#' Compares clustering across datasets using metrics described in associated TWL manuscript
#'
#' @param clus_save list of samples outputted from TWLsample function.
#' @param BURNIN number of samples devoted to burn-in.  Defaults to 2000.
#' @return outpu_lis a list of output metrics.  The first element is a list of lists of sample-specific pairwise cluster overlap.  The second element is an estimate of across all datasets cluster correspondence by averaging pairwise cluster overlap (the length is the vector therefore is the number of unique samples associated with at least 2 data sources.    
#' @examples
#' data(data_and_output)
#' \dontrun{clus_save <- TWLsample(misaligned_mat,misaligned,output_every=50,num_its=5000,manip=FALSE)
#' outpu_new <- pairwise_clus(clus_save,BURNIN=2000)
#' post_analy_cor(outpu_new,c("title1","title2","title3","title4","title5"),
#' tempfile(),ords='none') 
#' clus_labs <- post_analy_clus(outpu_new,clus_save,c(2:6),rep(0.6,5),c("title1","title2",
#' "title3","title4","title5"),tempfile())
#' output_nest <- cross_dat_analy(clus_save,4750)
#' }
#' @export
cross_dat_analy <- function(clus_save,BURNIN)
{
    num_its <- length(clus_save)
    INDS_NOT_INCLUDING_BURNIN <- (BURNIN+1):num_its
    num_datasets <- length(clus_save[[1]])
    sav_mat <- matrix(-99,num_datasets,num_datasets)
    sav_lis <- list()
    sav_lis_nams <- list()
    for(bb in 1:(num_datasets-1)){
	sav_lis[[bb]] <- list()
        sav_lis_nams[[bb]] <- list()
        for(vv in (bb+1):(num_datasets)){
            tmp <- sapply(clus_save[INDS_NOT_INCLUDING_BURNIN],function(x){(merge(x[[bb]],x[[vv]],'nam')[,clus.x==clus.y])})
            sav_mat[bb,vv] <- mean(tmp)
            row_me_tmp <- rowMeans(tmp)
            sav_lis_nams[[bb]][[vv]] <- sapply(clus_save[INDS_NOT_INCLUDING_BURNIN[1]],function(x){(merge(x[[bb]],x[[vv]],'nam')[,nam])})
            sav_lis[[bb]][[vv]] <- data.table::data.table(sav_lis_nams[[bb]][[vv]],row_me_tmp)
            names(sav_lis[[bb]][[vv]]) <- c('nam',paste0('cross_prec_',bb,'_',vv))
            print(paste0("comparing datasets ",bb," and ",vv))
        }
    }
    tmp <- lapply(clus_save[INDS_NOT_INCLUDING_BURNIN],function(x){Reduce(merge_on_nam,x)})
    cross_all <- data.table::data.table('nam'=tmp[[1]][,nam],'cross_prec_all' = rowMeans(sapply(tmp,equal_without_nam)))

    sav_mat2 <- MCMCpack::xpnd(MCMCpack::vech(t(sav_mat)))
    # Reseting the names the were temporily changed for the merge
    for(rr in seq_along(clus_save[[1]])){
        for(rrr in (BURNIN+1):length(clus_save)){
            data.table::setnames(clus_save[[rrr]][[rr]],c('nam','clus'))
        }
    }    
    outpu_list <- list(sav_lis,sav_mat2,cross_all)
    return(outpu_list)
}

clus <- clus.x <- clus.y <- j <- nam <- NULL 
