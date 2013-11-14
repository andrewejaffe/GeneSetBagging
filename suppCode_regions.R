####
#### First just compare distributions
####

### Simulate data set
n = 100
x = rnorm(n)
s2 = mean( (x-mean(x))^2)


# Get 5000 samples from the posterior
# which is inverse-wishart if you assume 
# the jeffreys priors

pSamples = rep(NA,5000)
for(i in 1:5000){
  # Sample the sd from the inverse wishart (inverse gaussian in 1-d)
  tmpvar = rwish((n-1),(n*s2)^(-1))		
  # Draw posterior sample
  pSamples[i] = rnorm(1,mean(x),sd=sqrt(tmpvar/n))
}

# Get 5000 bootstrap samples from bootstrap 
# distribution of the mean

bSamples = rep(NA,5000)
for(i in 1:5000){
 bSamples[i] = mean(sample(x,replace=T))
}


# Plot the distributions, they are (approximately)
# the same. If you increase n they converge

plot(density(pSamples))
lines(density(bSamples),col="red")


####
#### Now compare tail probabilities
####


##Now let's do this over and over for a tail probability

bTailProb = pTailProb = rep(NA,100)

for(k in 1:100){
n = 100
x = rnorm(n)

# Get 5000 samples from the posterior
# which is inverse-wishart if you assume 
# the jeffreys priors

pSamples = rep(NA,5000)
for(i in 1:5000){
  s2 = mean( (x-mean(x))^2)
  tmpsd = rwish((n-1),(n*s2)^(-1))
  pSamples[i] = rnorm(1,mean(x),sd=sqrt(tmpsd/n))
}

# Get 5000 bootstrap samples from bootstrap 
# distribution of the mean

bSamples = rep(NA,5000)
for(i in 1:5000){
 bSamples[i] = mean(sample(x,replace=T))
}

bTailProb[k] = mean(bSamples > 0.2)
pTailProb[k] = mean(pSamples > 0.2)
cat(k)
}

plot(bTailProb,pTailProb)
abline(c(0,1))
