#' Collapses data for a given collapsing function
#' @param data A dataset with J columns and n rows.
#' @param collapse A list of length J indicating the scoring function to collapse from.
#' @param constraint Constraint is either "rasch" or "gpcm" depending on which parameter constraints should be run.
#' @return A list containing the collapsed data and a indicator vector for which parameters to test with the IMT. If no collapsing has occurred, the default indicator vector tests all parameters of the last item.
#' @examples
#' data(dataset)
#' collapse = split(rep(c(1:4), 10), rep(1:10, each = 4))
#' my_data = collapse_data(dataset, collapse, "rasch")
#' # See vignette("IMT-vignette") for more examples.
#' @export
collapse_data <-
function(data, collapse, constraint){

if(missing(constraint))
	constraint = "gpcm"

## Getting the constants needed
J = dim(data)[[2]]
N = dim(data)[[1]]

if(missing(collapse)){
	collapse = list()
	for(j in 1:J){
		collapse[[j]] = seq(1:(apply(data,2, function(x) length(unique(x))))[j])
		}
	}

params_per = unlist(lapply(collapse, function(x) {length(unique(x))-as.numeric(constraint=="rasch")}))
nparam = sum(params_per)
np = (nparam*(nparam+1)/2)

## code to create a vector of indicators for the IMT
ind_vec = NULL
for(k in 1:J){
	x = collapse[[k]]
	for(i in unique(x)[-length(unique(x))]){
		ind_vec = c(ind_vec, as.numeric((as.numeric(sum(x == i) > 1) + as.numeric(sum(x == i+1)>1)) >= 1) )
		}
	if(constraint == "gpcm")
		ind_vec = c(ind_vec,as.numeric(length(unique(x))!=length(x)))
}

if(sum(ind_vec) == 0){
	ind_vec = c(rep(0, sum(params_per[1:(J-1)])), rep(1, params_per[J]))
}
## If there is no collapsing occuring, the previous lines will set the IMT to test all categories of the last item specifically

newdata = matrix(c(0), nrow = N, ncol = J) ##collapsed data
for(j in 1:J){
	for(i in 1:N){
		newdata[i,j] = collapse[[j]][data[i,j]]
	}
}


list(data = newdata, ind = ind_vec)

}
