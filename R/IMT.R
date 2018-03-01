#' Runs information matrix test for an information matrix test GPCM model.
#' @param mod An IMT GPCM model.
#' @param constraint Constraint is either "rasch" or "gpcm".
#' @param R number of iterations for simulation of the variance-covariance matrix.
#' @param ind_vec Vector of 0's and 1's for item-level parameters to be tested in the information matrix test.
#' @return A list containing the information matrix test statistic and the associated degrees of freedom.
#' @importFrom MASS ginv
#' @examples
#' data(dataset)
#' collapse = split(rep(c(1:4), 10), rep(1:10, each = 4))
#' my_data = collapse_data(dataset, collapse, "rasch")
#' model = gpcm_IMT(my_data$data, constraint = "rasch")
#' \donttest{
#' test_fit = IMT(model, "rasch", R = 5000, my_data$ind)
#' #This line of code takes longer than 10 seconds to run
#' pvalue = pchisq(test_fit$Tstat, test_fit$df, lower.tail = FALSE)
#' }
#' # See vignette("IMT-vignette") for more examples
#' @export
IMT <-
function(mod, constraint, R, ind_vec){

if(missing(constraint))
	constraint = "gpcm"

if(missing(R))
	R = 100000

## Getting the constants needed
J = length(mod$coefficients)

nparam = (length(unlist(coef(mod)))) - (mod$constraint == "rasch")*J
params_per = unlist(lapply(mod$coefficients, function(x) {length(x) - as.numeric(mod$constraint == "rasch")}))
np = (nparam*(nparam+1)/2)
N = sum(mod$patterns$obs)


####


### Data version of dmat

dmat = matrix(c(0), nrow = N, ncol = np)

for(i in 1:N){   D = mod$score_IMT[i,]%*%t(mod$score_IMT[i,]) - (mod$Hess_IMT[i,,]+t(mod$Hess_IMT[i,,]))/2;
                 dmat[i,] = D[upper.tri(D,T)]}

### Data estimate of r -- use this everywhere!
r = colMeans(dmat)

#### MLE version -- This is what you could hope to use in practice using the MC approximation
mleparms=unlist(mod$coefficients)
randbetas = lapply(mod$coefficients, function(x)x[-length(x)])
randalphas = lapply(mod$coefficients, function(x)x[length(x)])

### Generate one set of thetas for everyone
randtheta=rnorm(R,0,1)

randdata=create.item.scores.list(randtheta, randbetas, randalphas) ##slow

print("Simulation Data generated")

if(mod$constraint == "rasch"){
	randmodel=gpcm_nopt_IMT(randdata$data.frame,start.val=unlist(randbetas),constraint=mod$constraint)
}
if(mod$constraint != "rasch"){
	randmodel=gpcm_nopt_IMT(randdata$data.frame,start.val=mleparms,constraint=mod$constraint)
	}

print("Hessian computed")


### Compute Dmat using MC simulated from MLE
d1mat.alt = matrix(c(0), nrow = nrow(randdata$data.frame), ncol = np)
d2mat.alt = matrix(c(0), nrow = nrow(randdata$data.frame), ncol = np)
dmat.alt = matrix(c(0), nrow=nrow(randdata$data.frame), ncol=np)

for(i in 1:nrow(randdata$data.frame)){   d1mat.alt = randmodel$score_IMT[i,]%*%t(randmodel$score_IMT[i,]);
                                           d2mat.alt = (randmodel$Hess_IMT[i,,]+t(randmodel$Hess_IMT[i,,]))/2;
                  D.alt = d1mat.alt - d2mat.alt
                  dmat.alt[i,] = D.alt[upper.tri(D.alt,T)]
                 # if(i%%25000 == 0){ print(i)}
}

print("Statistic being computed..")

### Compute test statistic based on MC from MLE

Yhat.alt =  cbind(dmat.alt,randmodel$score_IMT)
V2.alt=ginv(ginv(t(Yhat.alt)%*%Yhat.alt/nrow(randdata$data.frame))[1:ncol(dmat.alt),1:ncol(dmat.alt)])

nparam_itemset = sum(ind_vec)
np_new = (nparam_itemset*(nparam_itemset+1)/2)


tempmat = ind_vec%*%t(ind_vec)

tm = tempmat[upper.tri(tempmat, T)]
rtemp = subset(r, tm==1)


newT = t(rtemp)%*%ginv(V2.alt[tm == 1, tm == 1])%*%rtemp*N

list(Tstat = newT, df = np_new)

}
