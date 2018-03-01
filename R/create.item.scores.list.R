create.item.scores.list <-
function(theta, betas, alphas){
	J =  length(alphas)
	K = length(theta) # length theta
	
	itemscores <- matrix(c(0), nrow=K, ncol=J)
	P = list()
	
	for(m in 1:J){
		bt = betas[[m]]
		bt = c(0, bt)
		
		m_j = length(bt)
			
		inside <- matrix(c(0), nrow=K, ncol=m_j)
		for(i in 1:m_j)
			inside[,i] <- as.numeric(alphas[m])*(theta-bt[i])


		cumul <- matrix(c(0), nrow=K, ncol=m_j)
		for(i in 1:length(theta)){
			cumul[i,] <- cumsum(inside[i,])
		}

		nums <- exp(cumul)
		P[[m]] <- nums/rowSums(nums)

		for(k in 1:K)
			itemscores[k,m] <-sample(1:m_j, 1, replace=T, prob=P[[m]][k,])
		}
	
	data <- as.data.frame(itemscores)
	
	list(data.frame = data, P=P)
}
