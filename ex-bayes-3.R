rm(list=ls())
d = read.csv("shuttle.txt", sep="\t")

y = d$damage
x = d$tempF

# Q1 calcul par max de vraisemblance
mle = summary(glm(y~x, family=binomial(link="probit")))
sigma = mle$cov.unscaled

# Q2
# Chacun sa prior !
# On peut par exemple prendre des priors N(0,10^2) indépendantes
prior = function(beta){
  return(dnorm(beta[1], 0, 10) * dnorm(beta[2], 0, 10))
}
logprior = function(beta){
  return(dnorm(beta[1], 0, 10, log=T) + dnorm(-beta[2], 0, 10, log=T))
}

# ou une prior impropre
prior = function(beta){return(1)}
logprior = function(beta){return(0)}

# Q3
# loi a posteriori
posterior = function(beta, y, x){
  p = pnorm(beta[1] + x*beta[2])
  lkd = prod(p^y) * prod((1-p)^(1-y))
  return(lkd * prior(beta))
}
logposterior = function(beta, y, x){
  p = pnorm(beta[1] + x*beta[2])
  loglkd = sum(y*log(p)) + sum((1-y)*log(1-p))
  if(!is.finite(loglkd)) return(-Inf)
  return(loglkd + logprior(beta))
}
# cette loi ne ressemble pas à une loi connue : on va avoir du mal à travailler avec directement


# Q4
# algorithme de Metropolis-Hastings

require(mvtnorm)

MH = function(beta0, niter, tau){
  beta = matrix(NA, nrow=niter, ncol=2)
  beta[1, ] = beta0
  
  for(i in 2:niter){
    proposal = rmvnorm(1, beta[i-1, ], tau^2*sigma)
    alpha = # probabilité d'acceptation
    if(runif(1) < alpha){
      # on accepte
    }
    else{
      # on rejette
    }
  }
  
  return(beta)
}


MH = function(beta0, niter, tau){
  beta = matrix(NA, nrow=niter, ncol=2)
  beta[1, ] = beta0
  acc = 0 # nombre d'acceptations
  
  for(i in 2:niter){
    proposal = rmvnorm(1, beta[i-1,], tau^2*sigma)
    logalpha = logposterior(proposal, y, x)-logposterior(beta[i-1,], y, x)
    if(log(runif(1)) < logalpha){
      beta[i,] = proposal
      acc = acc + 1
    }
    else{
      beta[i, ] = beta[i-1, ]
    }
  }
  print(acc/niter) #proportion d'acceptations
  return(beta)
}

niter = 2e3
b1 = MH(c(-3,0), niter, 1)
b2 = MH(c(2,0), niter, 1)
b3 = MH(c(2,-0.3), niter, 1)
b4 = MH(c(0, 0), niter, .1)
b5 = MH(c(0, -.01), niter, 1)


# étudions la sortie de l'algorithme
par(mfcol=c(3,3))
i = 2 # Changer en i=2 pour l'autre paramètre
# trace
plot(b1[, i], type="l")
plot(b2[, i], type="l")
plot(b3[, i], type="l")

# autocorrélations
acf(b1[100:niter, i])
acf(b2[100:niter, i])
acf(b3[100:niter, i])

# histogrammes
hist(b1[100:niter, i], breaks=50)
hist(b2[100:niter, i], breaks=50)
hist(b3[100:niter, i], breaks=50)

# Q6
# Effective Sample Size
niter/(2*sum(acf(b1[100:niter, 1], plot=F)$acf) - 1)
niter/(2*sum(acf(b2[100:niter, 1], plot=F)$acf) - 1)
niter/(2*sum(acf(b3[100:niter, 1], plot=F)$acf) - 1)

# Q7

(mu = colMeans(b2[100:niter, ]))
var(b2[100:niter, ])

# Q10
# Une estimation ponctuelle de la probabilité de défaillance en T=31 est
pnorm(mu[1] + mu[2]*1)
# on peut également faire un histogramme de cette probabilité
par(mfrow=c(1,1))
prob.pred = pnorm(b2[100:niter, 1] + 31*b2[100:niter, 2])
hist(prob.pred)
plot(density(prob.pred))

# Tout ceci peut bien sûr être fait via un package. Par exemple
library(mcmc)
library(expm) # pour la fonction expm qui permet d'optimiser la variance de la proposition
out = metrop(logposterior, c(2, 0), nbatch=1e4, 
             blen=1, scale=sqrtm(sigma), y=d$damage, x=d$tempF)
out$accept # taux d'acceptation
plot(ts(out$batch))

library(HDInterval)
hdi(out$batch[, 2]) # intervalle HPD pour beta2

# Q11 - non vue en cours
# Pour tirer selon une gaussienne N(mu,1) restreinte aux nombre positifs
xp = qnorm(runif(1) * pnorm(mu) + pnorm(-mu)) + mu
# Pour tirer selon une gaussienne N(mu,1) restreinte aux nombre positifs
xm = qnorm(runif(1) * pnorm(-mu)) + mu

# paramètres dans le Gibbs
X = cbind(1, x)
b = solve(t(X)%*%X)
a = b%*%t(X)
n = length(x)


# Gibbs à deux étapes : z sachant beta puis beta sachant z
GS = function(niter){
  beta = matrix(NA, nrow=niter, ncol=2)
  beta[1, ] = mle$coefficients[ , 1]
  
  for(i in 2:niter){
    mu = beta[i-1, 1] + beta[i-1, 2]*x
    zp = qnorm(runif(n) * pnorm(mu) + pnorm(-mu)) +mu
    zm = qnorm(runif(n) * pnorm(-mu)) + mu
    z = y*zp + (1-y)*zm
    beta[i, ] = rmvnorm(1, a%*%z, b)
  }
  return(beta)
}

bgs = GS(niter)
par(mfcol=c(3,1))
plot(bgs[ , 1], type="l")
acf(bgs[ , 1])
hist(bgs[ , 1])
colMeans(bgs)

# Q11f
system.time(MH(c(2, 0), niter,1))
system.time(GS(niter))
# Ces temps de calcul n'ont guère de sens pris seuls. Les deux méthodes prennent à peu près le même temps. 
# Pour les comparer, il faut voir quelle méthode produit le meilleur échantillon (par ex l'autocorrélation la plus faible).
sum(acf(b2[ , 1], plot=F)$acf)
sum(acf(bgs[ , 1], plot=F)$acf)
