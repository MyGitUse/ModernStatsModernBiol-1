#checking out the probability of seeing 3 mutations in a poisson distribution with lamda=5
dpois(x = 3, lambda = 5)

#what are the probabilities of see 0 to 12 mutations
dpois(x = 0:12, lambda = 5)



#make a barplot
barplot(dpois(0:12, 5), names.arg = 0:12, col = "red")


#how does R deal with categorical variables
genotype = c("AA","AO","BB","AO","OO","AO","AA","BO","BO",
             "AO","BB","AO","BO","AB","OO","AB","BB","AO","AO")
table(genotype)

#it automatically detects levels

genotypeF = factor(genotype)
levels(genotypeF)

#if levels are not in correct order, can reorder them:
#x = factor(x,levels(x)[c(4,5,1:3)])


#get the outcome of 15 coin tosses using Bernoulli model with probability equal to 0.5
rbinom(15, prob = 0.5, size = 1)

#this funciton does not have the same outcome everytime
#success and failure can have unequal outcomes in Bernoulli trial as long as the probabilities sum to 1

#uneven probability fo throwing a ball into right and left box
rbinom(12, prob = 2/3, size = 1)

#if we only care about how any times the ball goes into the right hand box we can set size equal to 
#the number of trials and x equal to 1 to tell the function to report a single number of observations

rbinom(1, prob = 2/3, size = 12) 

#in the above example we only need to show success because we know failure is equal to 1-p


#The number of successes in 15 Bernoulli trials with a probability of success of 0.3 is called a 
#binomial random variable or a random variable that follows the B(15,0.3) distribution

set.seed(235569515)
rbinom(1, prob = 0.3, size = 15)

#the complete probability mass distribution 
#set probabilities 
probabilities = dbinom(0:15, prob = 0.3, size = 15)
#display them to 2 decimal places
round(probabilities, 2)
#make a bar plot of probabilities
barplot(probabilities, names.arg = 0:15, col = "red")


#########find this out later



####1.3.4 A generative model for epitope detection

#ELISA false positive rate - 1%
#protein testing in 100 independent positions
#examine collection of 50 positions

#Verify by simulation that the sum of 50 independent Bernoulli variables with 
#p=0.01 is –to good enough approximation– the same as a Poisson(0.5) random variable.

set.seed(235569515)
sum(rbinom(50, prob = 0.01, size = 1))

dpois(x = 50, lambda = 0.5)

t.test(sum(rbinom(50, prob = 0.01, size = 1)),dpois(x = 50, lambda = 0.5))

load("data/e100.RData")

barplot(e100, ylim = c(0, 7), width = 0.7, xlim = c(-0.5, 100.5),
        names.arg = seq(along = e100), col = "darkolivegreen")

#If we look for the probability of seeing a number as big as 7 (or larger)
#This is, of course, the same as 1−P(  X≤6). 
#The probability P(X≤6) is the so-called cumulative distribution function at 6, 
#and R has the function ppois for computing it, which we can use in either of the 
#following two ways

1 - ppois(6, 0.5) #1st way
ppois(6, 0.5, lower.tail = FALSE) #if lower tail true prob are P[X≤x] ; 2nd way

#2nd way is better for small probabilities given that it doesnt have to deal with 
#the rounding of float integers

#hat are the chances that the maximum of 100 Poisson(0.5) trials is as large as 7? 
#We use extreme value analysis here:
#Order the data values x1,x2,..., x100 and rename them 
# x(1),x(2),...x(100), so that x(1) denotes the smallest and x(100) denotes the largest
#of the counts over the 100 positions
#together, x(1)...x(100) called the rank statistic of this sample of 100 values

maxes = replicate(100000, {
  max(rpois(100, 0.5))
})
table(maxes)

#Experiment with the random number generator that generates all possible numbers between 0
#and 1though the function called runif. Use it to generate a random variable with 
#4 levels (A, C, G, T) where a=1/8, c=3/8, g=3/8, and t=1/8

#not sure how to do this??
pvec = c(1/8, 3/8, 3/8,1/8)
runif(pvec, min=0, max=1)


#WWe often run simulation experiments to check whether the data we see are consistent 
#with the simplest possible four-box model where each box has the same probability
#Here we use a few R commands to generate such vectors of counts. First suppose 
#we have 8 characters of four different, equally likely types

pvec = rep(1/4, 4)
t(rmultinom(1, prob = pvec, size = 8)) #t transposes

#How do you interpret the difference between 
#rmultinom(n = 8, prob = pvec, size = 1) and rmultinom(n = 1, prob = pvec, size = 8)

rmultinom(n = 8, prob = pvec, size = 1)
#8 trials of throwing a ball into 1 of 4 boxes

rmultinom(n = 1, prob = pvec, size = 8)
#1 trial of throwing 8 balls into 1 of 4 boxes 

######1.4.1 Simulating for power########

#power has a special meaning in statistics. 
#It is the probability of detecting something if it is there
#(also called the true positive rate)

#lets say we have a DNA sequence and all nucleotides are equally likely
pvec=rep(1/4, 4)

#et’s determine if, by looking at a sequence of length n=20, we can 
#detect whether the original distribution of nucleotides is fair or 
#whether it comes from some other (“alternative”) process

#generate 1000 trials
obsunder0 = rmultinom(1000, prob = pvec, size = 20)
#show first 11 columns
obsunder0[, 1:11]
#Each column in the matrix obsunder0 is a simulated instance.

#this isnt enough to calculate power though
#We also need a measure of variability that will enable us 
#to describe how much variability is expected and how much is too much

#variability is computed computed as the sum of the squares of the 
#differences between the observed values and the expected values relative to the expected values

#therefore for the first column of data we would get
expected0 = pvec * 20
sum((obsunder0[, 1] - expected0)^2 / expected0)

#to calculate for all columns at once make a function for calculate stat:
stat = function(obsvd, exptd = 20 * pvec) {
  sum((obsvd - exptd)^2 / exptd)
}
stat(obsunder0[, 1])

#then compute for all instances and store in vector
S0 = apply(obsunder0, 2, stat) #the second variable 2 specifies to apply to columns (1 for rows)
summary(S0)

hist(S0, breaks = 25, col = "lavender", main = "")


#From the simulated data we can approximate, for instance, the 95% quantile 
#(a value that separates the smaller 95% of the values from the 5% largest values)

q95 = quantile(S0, probs = 0.95)
q95

#q95 represents our critical value for testing null hypothesis

#We must compute the probability that our test –based on the weighted sum-of-square differences– 
#will detect that the data in fact do not come from the null hypothesis. We compute the 
#probability of rejecting by simulation. We generate 1000 simulated instances from an 
#alternative process, parameterized by pvecA

pvecA = c(3/8, 1/4, 3/12, 1/8)
observed = rmultinom(1000, prob = pvecA, size = 20)
dim(observed)
observed[,1:7]

#see the mean for each nucleotide for all 1000 trials
apply(observed, 1, mean)

#determined what the expected values of A G T C in a 20 nucleotide sequence are
expectedA = pvecA * 20
expectedA

#determine variability
S1 = apply(observed, 2, stat)

#compare to the critical value we previously calculated
q95

#calculate how many trials are greater than the critical value
sum(S1 > q95)

#calculate the power 
power = mean(S1 > q95)
power

#With a sequence length of n=20 we have a power of about 20% 
#to detect the difference between the fair generating process and our alternative

#repeat the simulation experiments and suggest a new sequence length n
#that will ensure that the power is acceptable

pvecA = c(3/8, 1/4, 3/12, 1/8)
observed = rmultinom(1000, prob = pvecA, size = 28)
dim(observed)
observed[,1:7]

#see the mean for each nucleotide
apply(observed, 1, mean)

#determine expected
expectedA = pvecA * 28
expectedA

#determine variability
S2 = apply(observed, 2, stat)

#calculate critical value
q95

#calculate how many trials are greater than the critical value
sum(S2 > q95)

#calculate the power 
power = mean(S2 > q95)
power

#a length of ~28 nucleotides is needed for a power of 80%


#####classical statistics for classical data

#We didn’t need to simulate the data using Monte Carlo to compute the 95th percentiles; 
#there is an adequate theory to help us with the computations
#called the chi squared statistic (df=n-1)

#while chi squared is easier, the monte carlo method is helpful when many data point
#have low probabilities and their counts are mostly 0


##########CHAPTER SUMMARY#################


#Bernoulli distribution 
#used to represent a single binary trial such as a coin flip

#Binomial distribution
#used for the number of 1s in n binary trial and we generate the probabilities of seeing 
#k successes using the R function dbinom. 
#We also saw that we could simulate an n trial binomial using the function rbinom

#Poisson distribution
#most appropriate for small probabilities
#t has only one parameter λ
#the Poisson distribution for λ=np is approximately the same as the binomial distribution for 
#(n,p) if pis small

#Multinomial distribution 
#for discrete events that have more than two possible outcomes or levels


############EXERCISES################

#EXERCISE 1.1

#functions that start with an r as in rXXXX, where XXXX could be pois, binom, multinom
#are used to generate random discrete data 

#If we need a theoretical computation of a probability under one of these models, 
#we use the functions dXXXX, such as dbinom, which computes the probabilities of 
#events in the discrete binomial distribution (ie. density)

#dnorm computes the probability density function for the continuous normal distribution

#use pxxx functions to calculate cumulative tail probabilities (X>a)

#two other distributions
#(1) geom()  #geometric
#(2) hyper() #hypergeometic (not related to geometic; take samples without replacement)

####EXERCISE1.2

#How would you calculate the probability mass at the value X=2 for a 
#binomial B(10,0.3) with dbinom?

dbinom(2,10,0.3)

# Use dbinom to compute the cumulative distribution at the value 2, 
#corresponding to P(X<=2), and check with another R function/

sum(dbinom(c(2,1,0),10,0.3))

pbinom(2, 10, 0.3, lower.tail = TRUE, log.p = FALSE) #lower tail automatically set to true

###Exercise 1.3

prob_max <- function(n, lambda, m) {
 maxes <- rpois(n, lambda)
 mean(maxes<=m)
}


###Exercise 1.4
#Rewrite the function to have default values for its arguments 
#(i.e., values that are used by it if the argument is not specified in a call to the function).

prob_max <- function(n=1000, lambda=0.5, m=2) {
  maxes <- rpois(n, lambda)
  mean(maxes<=m)
}

####EXERCISE 1.5
#In the epitope example, use a simulation to find the probability of having a 
#maximum of 9 or larger in 100 trials. How many simulations do you need if 
#you would like to prove that "the probability is smaller than 0.000001"?

epitopes <- rpois(100, 0.5)
table(epitopes)
mean(epitopes>=9)

#number obs = 1/p-value = 1000000
maxes = replicate(1000000, {
  max(rpois(100, 0.5))
})
table(maxes)
mean(maxes>=9)


##################exercise 1.6

?Distributions

barplot(dcauchy(0:10), names.arg = 0:10, col = "red")

barplot(dchisq(0:10, df=2), names.arg = 0:10, col = "red")

barplot(dlnorm(0:10), names.arg = 0:10, col = "red")

barplot(dt(0:10, df=1), names.arg = 0:10, col = "red")

plot(df(0:10, 1, 10, log = FALSE), names.arg = 0:10, col = "red")

#############exercise 1/7
#Generate 100 instances of a Poisson(3) random variable. What is the mean? 
#What is the variance as computed by the R function var?

rep = replicate(100, {
  rpois(1,3)
})
mean(rep)
var(rep)


###########exercise1.8

#Test whether the C. elegans data is consistent with the uniform model 
#(all nucleotide fequencies the same) using a simulation.

#load the geonome library thing
library(BSgenome.Celegans.UCSC.ce2)

#read in celegens genome
genome <- BSgenome.Celegans.UCSC.ce2

#read dna seq from all chromosomes into a single object
#how to call on these more generally?
dnaseq <- DNAStringSet(c(genome$chrI,genome$chrII, genome$chrIII,genome$chrIV, genome$chrV, genome$chrX, genome$chrM))

#count bases in entire genome (nuclear and mitochondria)
#store in vector
basefreq <- alphabetFrequency(dnaseq, baseOnly=TRUE)

#find out number of bases in entire genome
genomesize <- width(dnaseq)

#calculate proportion of each observed nucleotide
baseprop <- basefreq/genomesize

#run monte carlo simulation with number of simulations equal to genome size
test<- rmultinom(genomesize, prob = rep(1/4, 4), size = 1)
testsize<- rowSums(test)
testprop <- testsize/genomesize

#run chi square test to determine if observed proportion differs from expected
chitest <- chisq.test(baseprop[,1:4], p = testprop)
chitest$p.value
#it does not differ
