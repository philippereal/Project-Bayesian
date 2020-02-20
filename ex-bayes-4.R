rm(list=ls())
data = read.csv("deathrate.csv")
y = data[, 17]
X = as.matrix(data[, 2:16])
n = length(y)

#Q1
reg.f = lm(y~X)
summary(reg.f)
betahat = reg.f$coefficients
residuals = reg.f$residuals
s2 = t(residuals)%*%residuals

X = cbind(1, X) # on ajoute une colonne de 1 pour beta_0

#Q2a
g = c(.1, 1, 10, 100, 1000)
# espérance a posteriori de beta0:
betahat[1] * g / (g + 1)

# espérance a posteriori de sigma^2
a = n/2
b = s2/2 + 1/(2*g+2) * ((t(betahat)%*%t(X)) %*% (X%*%betahat))
b / (a-1)

# Q2b
g = n
q = 2
X0 = X[,-(8:9)]
BF = (g+1)^(q/2) * 
  ((t(y)%*%y - g/(g+1) * t(y)%*%X0 %*% solve(t(X0)%*%X0) %*% t(X0)%*%y)/
  (t(y)%*%y - g/(g+1) * t(y)%*%X %*% solve(t(X)%*%X) %*% t(X)%*%y))^(n/2)
log10(BF)

#Q3
# fonction pour calculer la log-vraisemblance marginale
marglkd = function(gamma, X, g=n){
  q=sum(gamma)
  X1=X[ ,c(T,gamma)]
  if(q==0){return( -n/2 * log(t(y)%*%y))}
  m = -q/2*log(g+1) -
    n/2*log(t(y)%*%y - g/(g+1)* t(y)%*% X1 %*%
              solve(t(X1)%*%X1) %*%t(X1)%*%y)
return(m)
}

# calculons la log-vraisemblance marginale des 8 modèles
X_restreint = X[,1:4]
logprob3 = c(
marglkd(c(F,F,F), X_restreint),
marglkd(c(F,F,T), X_restreint),
marglkd(c(F,T,F), X_restreint),
marglkd(c(F,T,T), X_restreint),
marglkd(c(T,F,F), X_restreint),
marglkd(c(T,F,T), X_restreint),
marglkd(c(T,T,F), X_restreint),
marglkd(c(T,T,T), X_restreint))

# on peut ajouter une constante, qui évitera les erreurs numériques
logprob3 = logprob3-max(logprob3)
# les probabilités des modèles sont donc
prob3 = exp(logprob3)/sum(exp(logprob3))
round(prob3, 3)
# c'est le modèle (T, F, F) qui est de loin le plus probable a posteriori

#Q4

niter = 1e4 # nombre d'itérations
gamma = matrix(F, nrow = niter, ncol = 15)
gamma0 = sample(c(T, F), size = 15, replace = TRUE) # valeur initiale aléatoire
lkd = rep(0, niter)
modelnumber = rep(0, niter)

oldgamma = gamma0
for(i in 1:niter){
  newgamma = oldgamma
  for(j in 1:15){
    g1 = newgamma; g1[j]=TRUE
    g2 = newgamma; g2[j]=FALSE
    ml1 = marglkd(g1, X)
    ml2 = marglkd(g2, X)
    p = c(ml1,ml2)-min(ml1,ml2)
    # On souhaite tirer depuis une Bernoulli, avec probabilité de tirer TRUE égale à exp(p[1])/(exp(p[1])+exp(p[2])).
    # C'est ce que fait la ligne suivante. Notons que la fonction sample() calcule la constante de normalisation.
    newgamma[j] = sample(c(T,F), size=1, prob=exp(p)) 
  }
  gamma[i,] = newgamma
  lkd[i] = marglkd(newgamma, X )
  modelnumber[i] = sum(newgamma*2^(0:14))
  oldgamma = newgamma
}

colMeans(gamma)

# Vérifions le mélange de la chaîne de Markov à l'aide des autocorrélations.
par(mfrow=c(4,2))
for(i in 2:9) acf(as.numeric(gamma[,i]))
# Autocorrélation décroît rapidement. Pas besoin de sous-échantillonner.

# Vérifions la convergence + le mélange à l'aide de la trace (on utilise une moyenne glissante puisque les valeurs sont binaires).
library(zoo)
for(i in 2:15) plot(rollapply(gamma[,i], width=100, FUN=mean), type="l")

burnin = 500 # 500 itérations de burn-in
gammab = modelnumber[(burnin+1):niter] 
res = as.data.frame(table(gammab))
odo = order(res$Freq, decreasing=T)[1:50]
modcho = res$gammab[odo]
probtop50 = res$Freq[odo]/(niter-burnin)

indices = match(modcho, modelnumber)
cbind(probtop50, gamma[indices, ])


# Prédiction
Xnew = colMeans(X)
ynew.f = betahat %*% Xnew
hist(rnorm(niter, ynew.f, sqrt(s2/(n-p))))

ynew.b = rep(NA, niter)
for(i in 1:niter){
  X0 = X[, c(T, gamma[i,])]
  p0 = sum(gamma[i,])
  betahat0 = (lm(y~X0[,-1]))$coefficients
  s20 = sum((lm(y~X0[,-1]))$residuals^2)/(n-p0)
  sigma2=1/rgamma(1, n/2, s20/2 + .5/(g+1) * t(betahat0) %*% t(X0) %*% X0 %*% betahat0)
  beta = rmvnorm(1, g/(g+1) * betahat0, sigma2 *g/(g+1) * solve(t(X0)%*%X0))
  ynew.b[i] = beta %*% Xnew[c(T, gamma[i,])] + rnorm(1, 0, sqrt(sigma2))
}

hist(ynew.b)
