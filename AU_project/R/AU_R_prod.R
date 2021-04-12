#AU test written in R
#
#Eric Hester
#2021-4-12

setwd("D:/Projects/AU_projec/")

# Functions ---------------------------------------------------------------
AU <- function(M){
  #contract matrix by finding only otus found in the specified group
  #pool samples
  sM <- apply(M,1,sum)
  #rows where the OTU is seen on the host
  rows <- which(sM>0)
  #contract matrix
  m <- M[rows,] 
  
  sm <- apply(m,1,sum)
  ssm <- sum(sm)
  abun <- sm/ssm #pooled relative abundance
  m0 <- m>0 #logical
  k <- apply(m0,1,sum) #number of samples where OTU shows up
  n <- ncol(M) #total number of samples
  ubiq <- k/n #ubiquity
  
  return(data.frame(rows, abun, ubiq))
}

STATS <- function(a, ss){
  temp <- meshgrid(ss,a)
  SS <- temp$X
  A <- temp$Y
  PI <- 1-(1-A)^SS
  E <- apply(PI,1,sum)
  s <- apply(PI^2,1,sum)
  V <- E-s
  p0 <- s/E
  m0 <- (E^2)/s
  
  return(data.frame(m0,p0,E,V))
}


# Code-----

#Load your count data from your OTU table
otus <- read.delim("otu_table.txt")
otus <- otus[,-c(1,202)]#for example OTU table, tidying is necessary

#Open a file that contains groupings for your samples
#In the example case, there are four groups:
#Two coral groups: Acropora and Proites
#Two algae groups: CCA, Turf
meta <- read.delim("meta.txt")

#All calculations are done per group
group <- "Acropora"
#Get the columns that belong to the group above
cols <- which(meta$group==group)

M <- otus[,cols]

#Calculate the sampling depth of each sample (i.e., sequencing depth)
ss <- apply(M,2,sum)

#Calculate number of samples
num <- ncol(M)

#Calculate total number of OTUs seen in group
numrows <- nnzero(apply(M,1,sum))

#Get abundance/ubiquity table from AU function
groupMat <- AU(M)

#Plot scatter plot
plot(abun~ubiq,data=groupMat,log='y')



#Calculation of confidence levels
conf <- array(data=NA,dim=numrows)
for(i in 1:numrows){
  ab = groupMat$abun[i]
  ciMat <- STATS(ab,ss)
  #
  k <- num*groupMat$ubiq[i]
  m0_minus_k <- abs(ciMat$m0-k)
  conf[i] <- pbeta(ciMat$p0,k+1,m0_minus_k) 
}


#AU Expectation curve
u <- seq(0,.995,length.out = 501)
log_a <- seq(-5, -0.5, length.out = 501)

temp <- meshgrid(u,log_a)
U <- temp$X
L <- temp$Y

A <- 10^L
k <- num*u

a <- t((10^log_a))

statsMat <- STATS(a,ss)
E <- statsMat$E

u <- E/num
lines(u,a,log='y',col='red')

#print a table of number of Sporadic verus Core OTUs
output <- data.frame(row.names = c(group))
output$Stable <- length(which(1-conf<0.01))
output$Sporadic <- length(which(1-conf>0.01))
output$PerCore <- signif(output$Stable/numrows*100,3)
output$PerSatellite <- signif(output$Sporadice/numrows*100,3)

print(output)



