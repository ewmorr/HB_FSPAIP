
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
        x_subsamples[[i]] = vegan::rrarefy(
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
    vegan::vegdist(x_log, method = method, binary = F, diag = T, upper = T)
}


#############################################################################
#Take avgs over calculated metrics and/or the subsampled count tables
#i.e., average a list of matrices
avg_matrix_list = function(x){
    #function to average a list of matrices element-wise
    # silently converts list elements to matrix (from e.g., dist or df)
    # returns matrix with original row and col names if any
    list_len = length(x)
    
    #convert to matrix
    x = lapply(x, as.matrix)
    
    #get row and col names.
    names_row = rownames(x[[list_len]])
    names_col = colnames(x[[list_len]])
    
    sum_mat = Reduce('+', x)
    avg_mat = sum_mat/list_len
    
    #add names
    rownames(avg_mat) = names_row
    colnames(avg_mat) = names_col
    
    return(avg_mat)
}
#############################################################################
