

###file for generating simulation sample 
#Please set the number of lines (ncol) when scan() is used.
#"SUMallelfreq" is the (number of allele type)*(number of loci) table of cumulative allele frequency.
#	The total frequency par locus must be one.
#	The first line should be locus name and the first row is allele label.
#"sumLD" is the table of cumulative haplotype frequency. 
#	This simulation sample generator can receive two or three loci cluster.
#	The first line should be locus name.
#	From the first row to the second or third row should be the allele for each locus.
#	The last row should be cumulative haplotype frequency.
#	The total frequency par haplotype must be one.
#"LDlist" is the list of table of haplotype frequencies. 
#	Please divide into two groups by the number of loci par cluster.
#"rfreq" is the table of recombination rate.
###
#exaple files;
SUMallelefreq <- matrix(scan("sum27Xfreq.txt",skip=1,nlines=63),ncol=28,byrow=T)
sumLD1 <- matrix(scan("sumLD1.txt",skip=1),ncol=4,byrow=T)
sumLD2 <- matrix(scan("sumLD2.txt",skip=1),ncol=3,byrow=T)
sumLD3 <- matrix(scan("sumLD3.txt",skip=1),ncol=3,byrow=T)
sumLD4 <- matrix(scan("sumLD4.txt",skip=1),ncol=3,byrow=T)
LDlist2 <- list(sumLD2,sumLD3,sumLD4)
LDlist3 <- list(sumLD1)
rfreq <- matrix(scan("r.txt"),nrow=1)
###
#parameters
###
n <- 100000 #number of generating individuals
mu <- 0.0015 #mutation rate
alnum <- dim(SUMallelefreq)[1] #total number of kinds of allele
loci <- dim(SUMallelefreq)[2]-1 #the number of marker
ld2 <- c(9,18,26) #list of the first marker positions of two loci cluster. 
ld3 <- c(5) #list of the first marker positions of three loci cluster.

#####
#generating mother
#####
mother <- function(n,loci,alnum,SUMallelefreq,LDlist2,LDlist3){	
	moset <- matrix(0,loci,n*2)	 
	for(i in 1:(n*2)){
		mal <- matrix(0,loci,1)
		for(al in 1:loci){	
			x1 <- runif(1)			
			for (k1 in 1:alnum){	
				if (mal[al] != 0) break	
				if (x1 < SUMallelefreq[k1,al+1]) {
					mal[al,1] <- SUMallelefreq[k1,1]
				} else next
			}
		}
	
	#change independent loci to haplotype
		for (h in 1:length(LDlist2)){
		if(length(LDlist2) == 0) break
			hapf <- LDlist2[[h]]
			hapn <- dim(LDlist2[[h]])
			l <- ld2[h]
			x2 <- runif(1)					
			for (k2 in 1:hapn[1]){		
				if (mal[l,] != 0) break
				if (x2 < hapf[k2,3]){
					mal[l,] <- hapf[k2,1]
					mal[l+1,] <- hapf[k2,2]
				} else next
			}
		}

		for (g in 1:length(LDlist3)){
		if(length(LDlist3) == 0) break
			hapf <- LDlist3[[g]]
			hapn <- dim(LDlist3[[g]])
			l <- ld3[g]
			x3 <- runif(1)
			for (k3 in 1:hapn[1]){
				if (mal[l,] != 0) break
				if (x2 < hapf[k3,4]){
					mal[l,] <- hapf[k3,1]
					mal[l+1,] <- hapf[k3,2]
					mal[l+2,] <- hapf[k3,3]
				} else next
			}
		}
	moset[,i] <- mal
	}
	return(moset)
}

####
#generating father
####

father <- function(n,loci,alnum,SUMallelefreq,LDlist2,LDlist3){	
	faset <- matrix(0,loci,n)	 
	for(i in 1:(n)){
		fal <- matrix(0,loci,1)
		for(al in 1:loci){	
			x1 <-runif(1)			
			for (k1 in 1:alnum){	
				if (fal[al] != 0) break	
				if (x1 < SUMallelefreq[k1,al+1]) {
					fal[al,1] <- SUMallelefreq[k1,1]
				} else next
			}
		}
	
	#change independent loci to haplotype
		for (h in 1:length(LDlist2)){
		if(length(LDlist2) == 0) break
			hapf <- LDlist2[[h]]
			hapn <- dim(LDlist2[[h]])
			l <- ld2[h]
			x2 <- runif(1)					
			for (k2 in 1:hapn[1]){		
				if (fal[l,] != 0) break
				if (x2 < hapf[k2,3]){
					fal[l,] <- hapf[k2,1]
					fal[l+1,] <- hapf[k2,2]
				} else next
			}
		}

		for (g in 1:length(LDlist3)){
		if(length(LDlist3) == 0) break
			hapf <- LDlist3[[g]]
			hapn <- dim(LDlist3[[g]])
			l <- ld3[g]
			x3 <- runif(1)
			for (k3 in 1:hapn[1]){
				if (fal[l,] != 0) break
				if (x2 < hapf[k3,4]){
					fal[l,] <- hapf[k3,1]
					fal[l+1,] <- hapf[k3,2]
					fal[l+2,] <- hapf[k3,3]
				} else next
			}
		}
	faset[,i] <- fal
	}
	return(faset)
}

####
#generating daughter
####

daughter <- function(n,loci,rfreq,malset,faset){

#	"malset" <- result of mother(n)
#	"falset" <- result of father(n)

	dset <- matrix(0,loci,n*2) 
	drecomb <- matrix(0,loci,n) #record of recombinatin

for (i in 1:n){ 
	dre <- matrix(0,loci,1)
	dal<- matrix(0,loci,1)
	x1 <- runif(1) 　#the first locus
	if (x1 < 0.5){
		dal[1] <- malset[1,2*i-1]
		dre[1] <- 1
	}else {
		dal[1] <- malset[1,2*i]
		dre[1] <- 2
	}
	
	for (k in 2:loci){　#the following loci
		if (dre[k-1] == 1){
			s <- 1
		}else s <- 2
	
		x2 <- runif(1)

	switch(s,
		"1" = {
			if (x2 > rfreq[k-1]){
				dal[k] <- malset[k,2*i-1]
				dre[k] <- 1
			}else {dal[k] <- malset[k,2*i]
				dre[k] <- 2
			}
		},
		"2" = {
			if (x2 > rfreq[k-1]){
				dal[k] <- malset[k,2*i]
				dre[k] <- 2
			}else {dal[k] <- malset[k,2*i-1]
				dre[k] <- 1
			}
		}
	)　　
	}　
	dset[,2*i] <- dal
	dset[,2*i-1] <- falset[,i]
	drecomb[,i] <- dre
} 
list(dset,drecomb)
} 
 

####
#adding one step mutation
####

mutation <- function(n,mu,dsetm){

#	"dsetm" <- result of daughter[[1]]

dmut <- matrix(0,loci,n*2) #record of mutations
for (i in 1:n*2){
	for (a in 1:loci){
	x <- runif(1)
	if (x < mu*0.5){
		dsetm[a,i] <- dsetm[a,i] + 1
		dmut[a,i] <- 1
		} else if (x < mu){
		dsetm[a,i] <- dsetm[a,i] - 1
		dmut[a,i] <- 1
		} else next
	}
}
list(dsetm,dmut)
}



