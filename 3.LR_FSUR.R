####files for calculation of likelihood ratio (LR)
#file format#
#"allelefreq" is the (number of allele type)*(number of loci) table of allele frequency.
#	The first line is locus name and the first row is allele type label.
#	Each row corresponds to locus and must line up from the short telomere side to the long telomere side.
#"LD" is the table of haplotype frequency of two loci cluster. 
#	The first row is the allele of the former locus and the second row is the latter locus.
#	The third row is haplotype frequency.
#	If you want to use three loci cluster, please divide the file(e.g. haplotype "ABC" to haplotype "AB" and "BC").
#"LDlist" is the list of table of haplotype frequencies. 
#"rfreq" is the table of recombination rate.
#"PE" is the table of average probability of exclusion (APE).  
#	This needs when you use APE to calculate LR at a mutated locus.
#"genotype data" should be two rows for a female and one row for a male.
#	Each line corresponds to locus and must line up from the short telomere side to the long telomere side.
###
#exaple files;
allelefreq <- matrix(scan("27Xfreq.txt",skip=1),ncol=28,byrow=T) 
LD1_1 <- matrix(scan("LD1_1.txt",skip=1),ncol=3,byrow=T)
LD1_2 <- matrix(scan("LD1_2.txt",skip=1),ncol=3,byrow=T)
LD2 <- matrix(scan("LD2.txt",skip=1),ncol=3,byrow=T)
LD3 <- matrix(scan("LD3.txt",skip=1),ncol=3,byrow=T)
LD4 <- matrix(scan("LD4.txt",skip=1),ncol=3,byrow=T)
LDlist <- list(LD1_1,LD1_2,LD2,LD3,LD4)
rfreq <- matrix(scan("r.txt"),nrow=1)
PEsib <- matrix(scan("APEsib.txt"),nrow=1) 
####
#parameters
####
n <- 100000 #the number of pairs of sibling. 
N <- 425 #the number of haplotypic databese.
alnum <- dim(allelefreq)[1] #total number of kinds of allele
loci <- dim(allelefreq)[2]-1 #the number of marker
ld <- c(5,6,9,18,26) #list of marker positions of two loci cluster. 
#	Only first number of each two loci cluster should be listed.
#	If it is more than three loci cluster, list all position of loci except the last one.
mu <- 0.0015 #mutation rate. If you use APE to calculate LR at a mutated locus, mu should be 0.
newallele <- 0.001 #the frequency of unobserved allele. 

############
#change allele to its frequency
############

al_to_freq <- function(n,loci,alnum,allele,allelefreq){ 

	#"allele" <- "genotype data"

alfreq <- matrix(0,loci,2*n)
for (i in 1:(2*n)){
c1al <- matrix(0,loci,1)	
	for (al in 1:loci){
		for (k in 1:alnum){
			if (allele[al,i] == allelefreq[k,1]) {
			c1al[al] <- allelefreq[k,al+1]
			} 
		} #k
		if(c1al[al] == 0) c1al[al] <- newallele 
	}#al 
alfreq[,i] <- c1al
}#i
return(alfreq)
}

############
#likelihood of Full Sibling
############
#classification of IBS pattern

patternFS <- function(n,loci,mu,rfreq,child1,child2,child1freq,child2freq){

	#"child1" <- genotype of true child 
	#"child2" <- genotype of alleged child
	#"child1freq" <- result of "al_to_freq()" using genotype of true child
	#"child2freq" <- result of "al_to_freq()" using genotype of alleged child

pattern <- matrix(0,loci,n)  #IBS pattern
pIBD <- matrix(0,loci+1,n*5) #p(k0,k1,k2). p(k0) and p(k1) are devided into two lines which means Q and 1-Q.
ibs12 <- matrix(0,loci,n*6)　#line12=IBSalleles, line34=nonIBSalleles of child2,line56=nonIBSalleles of child1
ibs12freq <- matrix(0,loci,n*6)　#allele frequencies corresponding to ibs12

for(i in 1:n){
	p1 <- matrix(0,loci,1)　
	pIBD1 <- matrix(0,loci,5)　
	allele <- matrix(0,loci,6)　
	allelef <- matrix(0,loci,6)

	for(al in 1:loci){
		a1 <- 0
		a2 <- 0
		b1 <- 0
		b2 <- 0
		c1 <- 0
		c2 <- 0
		Q <- 0 　　　　
		
		if(child1[al,2*i-1] == child1[al,2*i])   a1 <- 1 　
		if(child2[al,2*i-1] == child2[al,2*i])   a2 <- 1 
		if(child1[al,2*i-1] == child2[al,2*i-1]) b1 <- 1 
		if(child1[al,2*i-1] == child2[al,2*i])   b2 <- 1 
		if(child1[al,2*i] == child2[al,2*i-1])   c1 <- 1 
		if(child1[al,2*i] == child2[al,2*i])     c2 <- 1 

		R1 <- (1-rfreq[al])^2+rfreq[al]^2
		R2 <-  2*rfreq[al]*(1-rfreq[al])

		if(a1*a2*b1*b2*c1*c2 == 1){ #AA AA
			p1[al] <- 1
			allele[al,c(1,3,5)] <- child2[al,2*i-1]
			allelef[al,c(1,3,5)] <- child2freq[al,2*i-1]
			Q <- (4-3*allelef[al,1])/(4-2*allelef[al,1])　
			pIBD1[al,1] <- Q*2*mu*R2
			pIBD1[al,2] <- (1-Q)*2*mu*R1
			pIBD1[al,3] <- Q*(2*mu*R1+(1-2*mu)*R2)
			pIBD1[al,4] <- (1-Q)*(2*mu*R2+(1-2*mu)*R1)
			pIBD1[al,5] <- Q*(1-2*mu)*R1+(1-Q)*(1-2*mu)*R2
		}else if(a1 == 1 && a2 == 0 && b1+b2 != 0){ #AA AB
			p1[al] <- 2
			allele[al,1] <- child1[al,2*i-1] 
			allele[al,3] <- child2[al,2*i-1]*b2+child2[al,2*i]*b1 
			allele[al,5] <- child1[al,2*i-1] 
			allelef[al,1] <- child1freq[al,2*i-1] 
			allelef[al,3] <- child2freq[al,2*i-1]*b2+child2freq[al,2*i]*b1 
			allelef[al,5] <- child1freq[al,2*i-1] 
			pIBD1[al,1] <- 2*mu*R1
			pIBD1[al,3] <- 2*mu*R2+(1-2*mu)*R1
			pIBD1[al,5] <- (1-2*mu)*R2
		}else if(a1 == 0 && a2 == 1 && b1+c1 != 0){ #AB AA
			p1[al] <- 3
			allele[al,1] <- child2[al,2*i-1] 
			allele[al,3] <- child2[al,2*i-1] 
			allele[al,5] <- child1[al,2*i-1]*c1+child1[al,2*i]*b1 
			allelef[al,1] <- child2freq[al,2*i-1] 
			allelef[al,3] <- child2freq[al,2*i-1] 
			allelef[al,5] <- child1freq[al,2*i-1]*c1+child1freq[al,2*i]*b1 
			pIBD1[al,1] <- 2*mu*R1
			pIBD1[al,3] <- 2*mu*R2+(1-2*mu)*R1
			pIBD1[al,5] <- (1-2*mu)*R2
		}else if(a1+a2 == 0 && b1+b2 == 1 && c1+c2 == 1){ #AB AB
			p1[al] <- 4
			allele[al,c(1,3,5)] <- child2[al,2*i-1]
			allele[al,c(2,4,6)] <- child2[al,2*i]
			allelef[al,c(1,3,5)] <- child2freq[al,2*i-1]
			allelef[al,c(2,4,6)] <- child2freq[al,2*i]
			Q <- ((4-3*allelef[al,1])/(4-2*allelef[al,1]) + (4-3*allelef[al,2])/(4-2*allelef[al,2])) /2
			pIBD1[al,1] <- Q*2*mu*R2
			pIBD1[al,2] <- (1-Q)*2*mu*R1
			pIBD1[al,3] <- Q*(2*mu*R1+(1-2*mu)*R2)
			pIBD1[al,4] <- (1-Q)*(2*mu*R2+(1-2*mu)*R1)
			pIBD1[al,5] <- Q*(1-2*mu)*R1+(1-Q)*(1-2*mu)*R2
		}else if(a1+a2 == 0 && b1+b2 == 1 || c1+c2 == 1){ #AB AC
			p1[al] <- 5
			allele[al,1] <- child2[al,2*i-1]*(b1+c1)+child2[al,2*i]*(b2+c2)
			allele[al,3] <- child2[al,2*i-1]*(b2+c2)+child2[al,2*i]*(b1+c1)
			allele[al,5] <- child1[al,2*i-1]*(c1+c2)+child1[al,2*i]*(b1+b2)
			allelef[al,1] <- child2freq[al,2*i-1]*(b1+c1)+child2freq[al,2*i]*(b2+c2)
			allelef[al,3] <- child2freq[al,2*i-1]*(b2+c2)+child2freq[al,2*i]*(b1+c1)
			allelef[al,5] <- child1freq[al,2*i-1]*(c1+c2)+child1freq[al,2*i]*(b1+b2)
			pIBD1[al,1] <- 2*mu*R1
			pIBD1[al,3] <- 2*mu*R2+(1-2*mu)*R1
			pIBD1[al,5] <- (1-2*mu)*R2
		}else { #AB CD || AB CC || AA BC || AA BB
			p1[al] <- 6 
			allele[al,3] <- child2[al,2*i-1]
			allele[al,4] <- child2[al,2*i]
			allele[al,5] <- child1[al,2*i-1]
			allele[al,6] <- child1[al,2*i]
			allelef[al,3] <- child2freq[al,2*i-1]
			allelef[al,4] <- child2freq[al,2*i]
			allelef[al,5] <- child1freq[al,2*i-1]
			allelef[al,6] <- child1freq[al,2*i]
			pIBD1[al,1] <- 2*mu*R1
			pIBD1[al,3] <- 2*mu*R2+(1-2*mu)*R1
			pIBD1[al,5] <- (1-2*mu)*R2
		}	
	}#al
	pattern[,i] <- p1
	ibs12[,(6*i-5):(6*i)] <- allele[,1:6]
	ibs12freq[,(6*i-5):(6*i)] <- allelef[,1:6]
	pIBD[1,5*i-4] <- 0.5*2*mu
	pIBD[1,5*i-2] <- 0.5*2*mu+0.5*(1-2*mu) 
	pIBD[1,5*i] <- 0.5*(1-2*mu) 
	pIBD[2:(loci+1),(5*i-4):(5*i)] <- pIBD1[,1:5]
}#i
list(pattern,pIBD,ibs12,ibs12freq)
}

###
#change independent loci to LD cluster
###

LDclusterFS <- function(n,LDlist,ld,pattern,ibs12){

	#pattern <- result of patternFS[[1]]
	#ibs12 <- result of patternFS[[3]]

LDfreq <- matrix(0,length(ld),2*n) #conditional frequencies of S2 non-IBD alleles
for (i in 1:n){
	ldfreq <- matrix(0,length(ld),2)
	for(j in 1:length(ld)){
		LD <- data.frame(LDlist[[j]])
		names(LD) <- c("L1","L2","HF")
		k <- ld[j]
		malef <- matrix(0,length(unique(LD[,1])),2)
		malef[,1]<- unique(LD[,1])
		for (a in 1:length(unique(LD[,1]))){
			for(h in 1:length(LD[,1])){
				if(malef[a,1] == LD[h,1]) malef[a,2] <- malef[a,2] + LD[h,3]
			}
		}
		Malef <- data.frame(L1=malef[,1],AF=malef[,2])

		if (pattern[k,i] == 1 || pattern[k,i] == 2 || pattern[k,i] == 3 || pattern[k,i] == 5){
			if (pattern[k+1,i] == 1 || pattern[k+1,i] == 2 || pattern[k+1,i] == 3 || pattern[k+1,i] == 5){
				locus1 <- c(ibs12[k,c(6*i-3,6*i-1)])
				locus2 <- c(rep(ibs12[k+1,6*i-3],2))
				id <- c(1:length(locus1))
				HAP <- data.frame(ID=id,L1=locus1,L2=locus2,CF=c(rep(0,length(locus1)))) 
				Htype <- merge(HAP,LD,all.x=T)
				Hapf <- merge(Htype,Malef,all.x=T)
				condfreq <- Hapf[order(Hapf$ID),] 
				condfreq[,4] <- condfreq[,5]/condfreq[,6] 
				x <- which(is.na(condfreq), arr.ind=TRUE) #new haplotype
				if (dim(x)[1] != 0) {
					for (b in 1:nrow(x)){
						if (x[b,2] == 6){#new allele
							condfreq[x[b,1],4] <- 1
						} else {
							condfreq[x[b,1],4] <- 1 / (condfreq[x[b,1],6]*N+1)
						}
					}
				}
				ldfreq[j,1] <- sum(condfreq[,4])/2

			} else if (pattern[k+1,i] == 4 || pattern[k+1,i] == 6){
				locus1 <- c(rep(ibs12[k,c(6*i-3,6*i-1)],2))
				locus2 <- c(rep(ibs12[k+1,6*i-3],2),rep(ibs12[k+1,6*i-2],2))
				id <- c(1:length(locus1))
				HAP <- data.frame(ID=id,L1=locus1,L2=locus2,CF=c(rep(0,length(locus1)))) 
				Htype <- merge(HAP,LD,all.x=T)
				Hapf <- merge(Htype,Malef,all.x=T)
				condfreq <- Hapf[order(Hapf$ID),] 
				condfreq[,4] <- condfreq[,5]/condfreq[,6] 
				x <- which(is.na(condfreq), arr.ind=TRUE) #new haplotype
				if (dim(x)[1] != 0) {
					for (b in 1:nrow(x)){
						if (x[b,2] == 6){#new allele
							condfreq[x[b,1],4] <- 1
						} else {
							condfreq[x[b,1],4] <- 1 / (condfreq[x[b,1],6]*N+1)
						}
					}
				}
				ldfreq[j,1] <- sum(condfreq[1:2,4])/2
				ldfreq[j,2] <- sum(condfreq[3:4,4])/2
			}
		} else if(pattern[k,i] == 4){
			if (pattern[k+1,i] == 1 || pattern[k+1,i] == 2 || pattern[k+1,i] == 3 || pattern[k+1,i] == 5){
				locus1 <- c(ibs12[k,c(6*i-3,6*i-2)])
				locus2 <- c(rep(ibs12[k+1,6*i-3],2))
				id <- c(1:length(locus1))
				HAP <- data.frame(ID=id,L1=locus1,L2=locus2,CF=c(rep(0,length(locus1)))) 
				Htype <- merge(HAP,LD,all.x=T)
				Hapf <- merge(Htype,Malef,all.x=T)
				condfreq <- Hapf[order(Hapf$ID),] 
				condfreq[,4] <- condfreq[,5]/condfreq[,6] 
				x <- which(is.na(condfreq), arr.ind=TRUE) #new haplotype
				if (dim(x)[1] != 0) {
					for (b in 1:nrow(x)){
						if (x[b,2] == 6){#new allele
							condfreq[x[b,1],4] <- 1
						} else {
							condfreq[x[b,1],4] <- 1 / (condfreq[x[b,1],6]*N+1)
						}
					}
				}
				ldfreq[j,1] <- sum(condfreq[,4])/2

			} else if (pattern[k+1,i] == 4 || pattern[k+1,i] == 6){
				locus1 <- c(rep(ibs12[k,c(6*i-3,6*i-2)],2))
				locus2 <- c(rep(ibs12[k+1,6*i-3],2),rep(ibs12[k+1,6*i-2],2))
				id <- c(1:length(locus1))
				HAP <- data.frame(ID=id,L1=locus1,L2=locus2,CF=c(rep(0,length(locus1)))) 
				Htype <- merge(HAP,LD,all.x=T)
				Hapf <- merge(Htype,Malef,all.x=T)
				condfreq <- Hapf[order(Hapf$ID),] 
				condfreq[,4] <- condfreq[,5]/condfreq[,6] 
				x <- which(is.na(condfreq), arr.ind=TRUE) #new haplotype
				if (dim(x)[1] != 0) {
					for (b in 1:nrow(x)){
						if (x[b,2] == 6){#new allele
							condfreq[x[b,1],4] <- 1
						} else {
							condfreq[x[b,1],4] <- 1 / (condfreq[x[b,1],6]*N+1)
						}
					}
				}
				ldfreq[j,1] <- sum(condfreq[1:2,4])/2
				ldfreq[j,2] <- sum(condfreq[3:4,4])/2
			}	
		} else if(pattern[k,i] == 6){
			if (pattern[k+1,i] == 1 || pattern[k+1,i] == 2 || pattern[k+1,i] == 3 || pattern[k+1,i] == 5){
				locus1 <- c(ibs12[k,c(6*i-3):(6*i)])
				locus2 <- c(rep(ibs12[k+1,6*i-3],4))
				id <- c(1:length(locus1))
				HAP <- data.frame(ID=id,L1=locus1,L2=locus2,CF=c(rep(0,length(locus1)))) 
				Htype <- merge(HAP,LD,all.x=T)
				Hapf <- merge(Htype,Malef,all.x=T)
				condfreq <- Hapf[order(Hapf$ID),] 
				condfreq[,4] <- condfreq[,5]/condfreq[,6] 
				x <- which(is.na(condfreq), arr.ind=TRUE) #new haplotype
				if (dim(x)[1] != 0) {
					for (b in 1:nrow(x)){
						if (x[b,2] == 6){#new allele
							condfreq[x[b,1],4] <- 1
						} else {
							condfreq[x[b,1],4] <- 1 / (condfreq[x[b,1],6]*N+1)
						}
					}
				}
				ldfreq[j,1] <- sum(condfreq[,4])/4

			} else if (pattern[k+1,i] == 4 || pattern[k+1,i] == 6){
				locus1 <- c(rep(ibs12[k,c(6*i-3):(6*i)],2))
				locus2 <- c(rep(ibs12[k+1,6*i-3],4),rep(ibs12[k+1,6*i-2],4))
				id <- c(1:length(locus1))
				HAP <- data.frame(ID=id,L1=locus1,L2=locus2,CF=c(rep(0,length(locus1)))) 
				Htype <- merge(HAP,LD,all.x=T)
				Hapf <- merge(Htype,Malef,all.x=T)
				condfreq <- Hapf[order(Hapf$ID),] 
				condfreq[,4] <- condfreq[,5]/condfreq[,6] 
				x <- which(is.na(condfreq), arr.ind=TRUE) #new haplotype
				if (dim(x)[1] != 0) {
					for (b in 1:nrow(x)){
						if (x[b,2] == 6){#new allele
							condfreq[x[b,1],4] <- 1
						} else {
							condfreq[x[b,1],4] <- 1 / (condfreq[x[b,1],6]*N+1)
						}
					}
				}
				ldfreq[j,1] <- sum(condfreq[1:4,4])/4
				ldfreq[j,2] <- sum(condfreq[5:8,4])/4
			}
		}
	LDfreq[,(2*i-1):(2*i)] <- ldfreq
	}#j
}#i
return(LDfreq)
}

###
#calculation of likelihood
###

LcalcFS <- function(n,loci,pattern,pIBD,freq,LDfreq){

	#pattern <- result of patternFS[[1]]
	#pIBD <- result of patternFS[[2]]
	#freq <- result of patternFS[[4]]
	#LDfreq <- result of LDclusterFS()

Lloci <- matrix(0,loci,n)
for(i in 1:n){
	LR1 <- matrix(0,loci,1)
	for(k in 1:loci){
		if (k %in% (ld+1)){
			s <- 1
			j <- which((ld+1) == k)
		}else s <- 2

		switch(s, 
		"1" = {if(pattern[k,i] == 1){
				LR1[k,] <- pIBD[k,5*i-2]*2*freq[k,6*i-3] + pIBD[k,5*i-1]*2*LDfreq[j,2*i-1] + pIBD[k,5*i]
			} else if(pattern[k,i] == 4){
				LR1[k,] <- pIBD[k,5*i-2]*(freq[k,6*i-3]+freq[k,6*i-2]) + pIBD[k,5*i-1]*(LDfreq[j,2*i-1]+LDfreq[j,2*i]) + pIBD[k,5*i]
			} else if(pattern[k,i] == 2 || pattern[k,i] == 3 || pattern[k,i] == 5) {
				LR1[k,] <- pIBD[k,5*i-2]*freq[k,6*i-3] + pIBD[k,5*i-1]*LDfreq[j,2*i-1] 			
			} else if(pattern[k,i] == 6){
				LR1[k,] <- pIBD[k,5*i-4]*(freq[k,6*i-3]+freq[k,6*i-2])/2 + pIBD[k,5*i-3]*(LDfreq[j,2*i-1]+LDfreq[j,2*i])/2
			}
		},
		"2" = {if(pattern[k,i] == 1){
				LR1[k,] <- (pIBD[k,5*i-2]+pIBD[k,5*i-1])*2*freq[k,6*i-3] + pIBD[k,5*i]				
			} else if(pattern[k,i] == 4){
				LR1[k,] <- (pIBD[k,5*i-2]+pIBD[k,5*i-1])*(freq[k,6*i-3]+freq[k,6*i-2]) + pIBD[k,5*i]				
			} else if(pattern[k,i] == 2 || pattern[k,i] == 3 || pattern[k,i] == 5) {
				LR1[k,] <- (pIBD[k,5*i-2]+pIBD[k,5*i-1])*freq[k,6*i-3] 
			} else if(pattern[k,i] == 6){
				LR1[k,] <- (pIBD[k,5*i-4]+pIBD[k,5*i-3])*(freq[k,6*i-3]+freq[k,6*i-2])/2
			}
		},
		)
	}#k
	Lloci[,i] <- LR1
}#i
return(Lloci)
}

############
#likelihood of unrelated female
############
#classification of genotype

patternUR <- function(n,loci,child2,child2freq){

	#"child2" <- genotype of alleged child
	#"child2freq" <- result of "al_to_freq()" using genotype of alleged child

pattern <- matrix(0,loci,n)  
alfreq <- matrix(0,loci,n) 
for(i in 1:n){
	p1 <- matrix(0,loci,1)
	freq <- matrix(0,loci,1)
	for(al in 1:loci){
		if(child2[al,2*i-1] == child2[al,2*i]) {#homo
			p1[al] <- 1
			freq[al] <- child2freq[al,2*i-1]^2 
		}else {
			p1[al] <- 2 #hetero
			freq[al] <- 2*child2freq[al,2*i-1]*child2freq[al,2*i]
		}
	}#al
	pattern[,i] <- p1
	alfreq[,i] <- freq
}#i
list(pattern,alfreq)
}

###
#change independent loci to LD cluster
###

LDclusterUR <- function(n,ld,LDlist,child2,pattern,ldfreq){

	#"child2" <- genotype of alleged child
	#pattern <- result of patternUR[[1]]
	#ldfreq <- result of patternUR[[2]]

for (i in 1:n){
	for(j in 1:length(ld)){
		LD <- data.frame(LDlist[[j]])
		names(LD) <- c("L1","L2","HF")
		k <- ld[j]
		malef <- matrix(0,length(unique(LD[,1])),2)
		malef[,1]<- unique(LD[,1])
		for (a in 1:length(unique(LD[,1]))){
			for(h in 1:length(LD[,1])){
				if(malef[a,1] == LD[h,1]) malef[a,2] <- malef[a,2] + LD[h,3]
			}
		}
		Malef <- data.frame(L1=malef[,1],AF=malef[,2])

	if(pattern[k,i] == 1 && pattern[k+1,i] == 1){
		locus1 <- c(child2[k,2*i-1])
		locus2 <- c(child2[k+1,2*i-1])
		id <- c(1:length(locus1))
		HAP <- data.frame(ID=id,L1=locus1,L2=locus2,CF=c(rep(0,length(locus1)))) 
		Htype <- merge(HAP,LD,all.x=T)
		Hapf <- merge(Htype,Malef,all.x=T)
		condfreq <- Hapf[order(Hapf$ID),] 
		condfreq[,4] <- condfreq[,5]/condfreq[,6] 
		x <- which(is.na(condfreq), arr.ind=TRUE) #new haplotype
		if (dim(x)[1] != 0) {
			for (b in 1:nrow(x)){
				if (x[b,2] == 6){#new allele
					condfreq[x[b,1],4] <- 1
				} else {
					condfreq[x[b,1],4] <- 1 / (condfreq[x[b,1],6]*N+1)
				}
			}
		}
		ldfreq[k+1,i] <- condfreq[1,4]^2
	} else if(pattern[k,i] == 1 && pattern[k+1,i] == 2){
		locus1 <- c(rep(child2[k,2*i-1],2))
		locus2 <- c(child2[k+1,(2*i-1):(2*i)])
		id <- c(1:length(locus1))
		HAP <- data.frame(ID=id,L1=locus1,L2=locus2,CF=c(rep(0,length(locus1)))) 
		Htype <- merge(HAP,LD,all.x=T)
		Hapf <- merge(Htype,Malef,all.x=T)
		condfreq <- Hapf[order(Hapf$ID),] 
		condfreq[,4] <- condfreq[,5]/condfreq[,6] 
		x <- which(is.na(condfreq), arr.ind=TRUE) #new haplotype
		if (dim(x)[1] != 0) {
			for (b in 1:nrow(x)){
				if (x[b,2] == 6){#new allele
					condfreq[x[b,1],4] <- 1
				} else {
					condfreq[x[b,1],4] <- 1 / (condfreq[x[b,1],6]*N+1)
				}
			}
		}
		ldfreq[k+1,i] <- 2*condfreq[1,4]*condfreq[2,4]
	} else if(pattern[k,i] == 2 && pattern[k+1,i] == 1){
		locus1 <- c(child2[k,(2*i-1):(2*i)])
		locus2 <- c(rep(child2[k+1,2*i-1],2))
		id <- c(1:length(locus1))
		HAP <- data.frame(ID=id,L1=locus1,L2=locus2,CF=c(rep(0,length(locus1)))) 
		Htype <- merge(HAP,LD,all.x=T)
		Hapf <- merge(Htype,Malef,all.x=T)
		condfreq <- Hapf[order(Hapf$ID),] 
		condfreq[,4] <- condfreq[,5]/condfreq[,6] 
		x <- which(is.na(condfreq), arr.ind=TRUE) #new haplotype
		if (dim(x)[1] != 0) {
			for (b in 1:nrow(x)){
				if (x[b,2] == 6){#new allele
					condfreq[x[b,1],4] <- 1
				} else {
					condfreq[x[b,1],4] <- 1 / (condfreq[x[b,1],6]*N+1)
				}
			}
		}
		ldfreq[k+1,i] <- condfreq[1,4]*condfreq[2,4]
	} else if(pattern[k,i] == 2 && pattern[k+1,i] == 2){
		locus1 <- c(rep(child2[k,(2*i-1):(2*i)],2))
		locus2 <- c(rep(child2[k+1,2*i-1],2),rep(child2[k+1,2*i],2))
		id <- c(1:length(locus1))
		HAP <- data.frame(ID=id,L1=locus1,L2=locus2,CF=c(rep(0,length(locus1)))) 
		Htype <- merge(HAP,LD,all.x=T)
		Hapf <- merge(Htype,Malef,all.x=T)
		condfreq <- Hapf[order(Hapf$ID),] 
		condfreq[,4] <- condfreq[,5]/condfreq[,6] 
		x <- which(is.na(condfreq), arr.ind=TRUE) #new haplotype
		if (dim(x)[1] != 0) {
			for (b in 1:nrow(x)){
				if (x[b,2] == 6){#new allele
					condfreq[x[b,1],4] <- 1
				} else {
					condfreq[x[b,1],4] <- 1 / (condfreq[x[b,1],6]*N+1)
				}
			}
		}
		ldfreq[k+1,i] <- condfreq[1,4]*condfreq[4,4]+condfreq[2,4]*condfreq[3,4]
	}
	}#j
}#i
return(ldfreq)
}

############
#calculation of likelihood ratio
############
LRcalc <- function(n,loci,Lx,Ly){

	#LRx <- likelihood of hypothesis questioned (i.e. result of LcalcFS())
	#LRy <- likelihood of alternative hypothesis (i.e. result of LDclusterUR())

LRloci <- matrix(0,loci,n)
for(i in 1:n){
	LR1 <- matrix(0,loci,1)
	for(k in 1:loci){
		LR1[k] <- Lx[k,i]/Ly[k,i]
	}
LRloci[,i] <- LR1	
}
LRtotal <- apply(LRloci, MARGIN=2, prod)
list(LRloci, LRtotal)
}

############
#calculation of likelihood ratio using APE
############
LRcalcPE <- function(n,APE,pattern,Lx,Ly){
	
	#APE <- PEsib
	#pattern <- result of patternFS[[1]]
	#Lx <- likelihood of hypothesis questioned (i.e. result of LcalcFS())
	#Ly <- likelihood of alternative hypothesis (i.e. result of LDclusterUR())

LRloci <- matrix(0,loci,n)
for(i in 1:n){
	LR1 <- matrix(0,loci,1)
	for(k in 1:loci){
		if(pattern[k,i] == 6){
			LR1[k] <- APE[k]
		} else {
			LR1[k] <- Lx[k,i]/Ly[k,i]
		}
	}
LRloci[,i] <- LR1	
}
LRtotal <- apply(LRloci, MARGIN=2, prod)
list(LRloci, LRtotal)
}

######
#how to use;
######
	allele <- daughter1m[[1]] #see "simulation_sample_generator"
child1freq <- al_to_freq(n,loci,alnum,allele,allelefreq)
	allele <- daughter2m[[1]] 
child2freq <- al_to_freq(n,loci,alnum,allele,allelefreq)
	child1 <- daughter1m[[1]] 
	child2 <- daughter2m[[1]]

patternfs <- patternFS(n,loci,mu,rfreq,child1,child2,child1freq,child2freq)
	pattern <- patternfs[[1]]
	ibs12 <- patternfs[[3]]
LDfreq <- LDclusterFS(n,LDlist,ld,pattern,ibs12)
	pIBD <- patternfs[[2]]
	freq <- patternfs[[4]]
Lx <- LcalcFS(n,loci,pattern,pIBD,freq,LDfreq)

patternur <- patternUR(n,loci,child2,child2freq)
	pattern <- patternur[[1]]
	ldfreq <-  patternur[[2]]
Ly <- LDclusterUR(n,ld,LDlist,child2,pattern,ldfreq)

LR <- LRcalc(n,loci,Lx,Ly) #Use 'LRcalcPE(n,APE,pattern,Lx,Ly)' when you want to use APE. 
