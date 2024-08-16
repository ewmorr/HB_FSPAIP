
############################################################
#Functions for multiple subsampling (rarefaction) of species counts tables
# and analysis of rarefied diversity metrics

#############################################################################
#perform rarefaction
multiple_subsamples = function(x = NULL, depth = NULL, iterations = NULL){
    # Performs multiple subsamples with vegan::rrarefy and returns 
    # a list of length iterations containing the resulting rarefied tables.
    # Attempts to convert non matrix objects (e.g., data.frame) to 
    # matrix and filters out columns with less than depth.
    # Input: x is a table of samples (rows) x species (columns)
    # depth is the desired sample depth per sample
    # iterations is the number of samples
    
    if(require(vegan) != T){
        install.packages(vegan)
    }
    library(vegan)
    
    if(is.matrix(x) == F){
        x = as.matrix(x)
    }
    
    x.min = x[rowSums(x) >= depth,]
    
    x_subsamples = list()
    
    for(i in 1:iterations){
        #rrarefy apparently throws a warnings if there are no counts of 1 in the data
        # annoying... but if you get that warning ignore
        x_subsamples[[i]] = rrarefy(
            x = x.min,
            sample = depth
        )
    }
    return(x_subsamples)
}

#############################################################################
#calculate richness of each sample (i.e., number of non zero entries per row)
richness_calc = function(x){
    x[x>0] = 1
    return(rowSums(x))
}
#############################################################################
#function to log+1 transform counts before dist
log_dist = function(x, method = "bray"){
    x_log = log(x+1)
    vegdist(x_log, method = method, binary = F, diag = T, upper = T)
}


#############################################################################
#Take avgs over calculated metrics and/or the subsampled count tables
#i.e., average a list of matrices
avg_matrix_list = function(x){
    #function to average a list of matrices element-wise
    # silently converts non matrix objects (e.g., dist or vector) 
    # to a matrix and returns matrix with original row and col names if any
    list_len = length(x)
    
    temp_mat = as.matrix(x[[list_len]])
    #get table size from last in list
    n_row = nrow(temp_mat)
    n_col = ncol(temp_mat)
    #get row and col names. 
    names_row = rownames(temp_mat)
    names_col = colnames(temp_mat)
    
    sum_vec = vector(length = n_row*n_col, mode = "numeric")
    
    for(i in 1:list_len){
        #make sure matrix
        x[[i]] = as.matrix(x[[i]])
        #1D the matrix and take rolling sum
        sum_vec = sum_vec + as.vector(x[[i]])
    }
    avg_vec = sum_vec/list_len
    
    #convert to original matrix dimensions and add names
    final_mat = matrix(avg_vec, nrow = n_row, ncol = n_col)
    rownames(final_mat) = names_row
    colnames(final_mat) = names_col
    
    return(final_mat)
}
#############################################################################
