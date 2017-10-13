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
PEfd <- matrix(scan("APEfd.txt"),nrow=1)
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
#likelihood of Father-daughter
############
#classification of IBS pattern

patternFD <- function(n,loci,father,child,childfreq){

	#"father" <- genotype of alleged father
	#"child" <- genotype of daughter
	#"childfreq" <- result of "al_to_freq()" using genotype of daughter

pattern <- matrix(0,loci,n)  #IBS pattern
ibsfd <- matrix(0,loci,n*2)　#maternal alleles
ibsfdfreq <- matrix(0,loci,n)　#maternal alleles frequencies

for(i in 1:n){
	p1 <- matrix(0,loci,1)　
	allele <- matrix(0,loci,2)　
	allelef <- matrix(0,loci,1)

	for(al in 1:loci){
		if(child[al,2*i-1] == child[al,2*i]){  #homo
			if(father[al,i] == child[al,2*i-1]){ 
				p1[al] <- 1 #A AA
				allele[al,1] <- child[al,2*i-1]
				allelef[al,] <- childfreq[al,2*i-1]
			} else {p1[al] <- 3　#A BB
				allele[al,1] <- child[al,2*i-1] 
				allelef[al,] <- childfreq[al,2*i-1]
			}    
		}else if(child[al,2*i-1] != child[al,2*i]){ #hetero
			if (father[al,i] == child[al,2*i-1]){　
				p1[al] <- 2 #A AB
				allele[al,1] <- child[al,2*i] 
				allelef[al,] <- childfreq[al,2*i]
			} else if(father[al,i] == child[al,2*i]){
				p1[al] <- 2 
				allele[al,1] <- child[al,2*i-1] 
				allelef[al,] <- childfreq[al,2*i-1]
		      }else {p1[al] <- 4 #A BC
				allele[al,1:2] <- child[al,(2*i-1):(2*i)] 
				allelef[al,] <- sum(childfreq[al,(2*i-1):(2*i)])/2
			}
		}
	}#al
	pattern[,i] <- p1
	ibsfd[,(2*i-1):(2*i)] <- allele
	ibsfdfreq[,i] <- allelef
}#i
list(pattern,ibsfd,ibsfdfreq)
}

###
#change independent loci to LD cluster
###

LDclusterFD <- function(n,LDlist,ld,pattern,ibsfd,ldfreq){

	#"pattern" <- result of patternFD[[1]]
	#"ibsfd" <- result of patternFD[[2]]
	#"ldfreq" <- result of patternFD[[3]]

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

		if (pattern[k,i] == 1 || pattern[k,i] == 2){ #IBS(+) 
			if (pattern[k+1,i] == 1 || pattern[k+1,i] == 2){ #IBS(+)
				locus1 <- c(ibsfd[k,2*i-1])
				locus2 <- c(ibsfd[k+1,2*i-1])
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
				ldfreq[k+1,i] <- condfreq[,4]

			} else if (pattern[k+1,i] == 3 || pattern[k+1,i] == 4){ #IBS(-)
				locus1 <- c(rep(ibsfd[k,2*i-1],2))
				locus2 <- c(ibsfd[k+1,(2*i-1):(2*i)])
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
				ldfreq[k+1,i] <- sum(condfreq[,4])/2			
			}
		} else if(pattern[k,i] == 3 || pattern[k,i] == 4 ){ # IBS(-)
			if (pattern[k+1,i] == 1 || pattern[k+1,i] == 2){ # IBS(+)
				locus1 <- c(ibsfd[k,(2*i-1):(2*i)])
				locus2 <- c(rep(ibsfd[k+1,2*i-1],2))
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
				ldfreq[k+1,i] <- sum(condfreq[,4])/2

			} else if (pattern[k+1,i] == 3 || pattern[k+1,i] == 4){ # IBS(-)
				locus1 <- c(rep(ibsfd[k,(2*i-1):(2*i)],2))
				locus2 <- c(rep(ibsfd[k+1,2*i-1],2),rep(ibsfd[k+1,2*i],2))
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
				ldfreq[k+1,i] <- sum(condfreq[c(1,4),4])/2+sum(condfreq[c(2,3),4])/2
			}	
		}
	}#j
}#i
return(ldfreq)
}

###
#calculation of likelihood
###

LcalcFD <- function(n,loci,mu,pattern,ldfreq){

	#"pattern" <- result of patternFD[[1]]
	#"ldfreq" <- result of LDclusterFD

Lloci <- matrix(0,loci,n)
pIBD0 <- mu
pIBD1 <- 1-mu
for(i in 1:n){
	LR1 <- matrix(0,loci,1)
	for(k in 1:loci){
		if(pattern[k,i] == 1 || pattern[k,i] == 2){
				LR1[k,] <- pIBD1*ldfreq[k,i]
		} else if(pattern[k,i] == 3 || pattern[k,i] == 4){
				LR1[k,] <- pIBD0*ldfreq[k,i]
		}
	}#k
	Lloci[,i] <- LR1
}#i
return(Lloci)
}

############
#likelihood of unrelated male and female
############
#classification of genotype

patternUR <- function(n,loci,child,childfreq){

	#"child" <- genotype of alleged child
	#"childfreq" <- result of "al_to_freq()" using genotype of alleged child

pattern <- matrix(0,loci,n)  
alfreq <- matrix(0,loci,n) 
for(i in 1:n){
	p1 <- matrix(0,loci,1)
	freq <- matrix(0,loci,1)
	for(al in 1:loci){
		if(child[al,2*i-1] == child[al,2*i]) {#homo
			p1[al] <- 1
			freq[al] <- childfreq[al,2*i-1]^2 
		}else {
			p1[al] <- 2 #hetero
			freq[al] <- 2*childfreq[al,2*i-1]*childfreq[al,2*i]
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

LDclusterUR <- function(n,ld,LDlist,child,pattern,ldfreq){

	#"child" <- genotype of alleged child
	#"pattern" <- result of patternUR[[1]]
	#"ldfreq" <- result of patternUR[[2]]

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
		locus1 <- c(child[k,2*i-1])
		locus2 <- c(child[k+1,2*i-1])
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
		locus1 <- c(rep(child[k,2*i-1],2))
		locus2 <- c(child[k+1,(2*i-1):(2*i)])
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
		locus1 <- c(child[k,(2*i-1):(2*i)])
		locus2 <- c(rep(child[k+1,2*i-1],2))
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
		locus1 <- c(rep(child[k,(2*i-1):(2*i)],2))
		locus2 <- c(rep(child[k+1,2*i-1],2),rep(child[k+1,2*i],2))
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

	#"LRx" <- likelihood of hypothesis questioned (i.e. result of LcalcFD())
	#"LRy" <- likelihood of alternative hypothesis (i.e. result of LDclusterUR())

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
	
	#APE <- PEfd
	#pattern <- result of patternFD[[1]]
	#LRx <- likelihood of hypothesis questioned (i.e. result of LcalcFD())
	#LRy <- likelihood of alternative hypothesis (i.e. result of LDclusterUR())

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
childfreq <- al_to_freq(n,loci,alnum,allele,allelefreq)
	father <- father1
	child <- daughter1m
patternfd <- patternFD(n,loci,father,child,childfreq)
	pattern <- patternfd[[1]]
	ibsfd <- patternfd[[2]]
	ldfreq <- patternfd[[3]]
ldfreq <- LDclusterFD(n,LDlist,ld,pattern,ibsfd,ldfreq)
Lx <- LcalcFD(n,loci,mu,pattern,ldfreq)

patternur <- patternUR(n,loci,child,childfreq)
	pattern <- patternur[[1]]
	ldfreq <- patternur[[2]]
Ly <- LDclusterUR(n,ld,LDlist,child,pattern,ldfreq)

LR <- LRcalc(n,loci,Lx,Ly)  #Use 'LRcalcPE(n,loci,Lx,Ly)' when you want to use APE. 