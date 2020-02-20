# Estimons l'intégrale de cos(x)*dnorm(x) sur R

n = 1e5
x = rnorm(n)
mean(cos(x))

plot(cumsum(x)/(1:n), type="l")

integrate(function(x) cos(x)*dnorm(x), -Inf, Inf)

# Estimons P[X > a] où X ~ N(0,1)

n = 1e5
x = rnorm(n)
a = 5
mean(x > a)
pnorm(a, lower=F)

# Importance sampling avec une loi exponentielle et une loi de Cauchy
lambda = 0.1
y = rexp(n, lambda)
mean((y>a) * dnorm(y) / dexp(y, lambda))

z = rcauchy(n)
mean((z>a) * dnorm(z) / dcauchy(z))
# Ces estimateurs sont bien plus efficaces.


# Feuille 2

deputes = read.csv2("deputes2019.csv")
attach(deputes)
par(mfrow=c(1, 1))

#Q1
#Exploration des données
summary(questions_orales)
barplot(table(questions_orales))
(n = length(questions_orales))
(nh = sum(sexe=="H"))
(nf = sum(sexe=="F"))
(qtot = sum(questions_orales))
(qh = sum(questions_orales[sexe=="H"]))
(qf = sum(questions_orales[sexe=="F"]))

#Le modèle de Poisson convient-il ?
lambdahat = qtot/n # EMV
points(0:14, dpois(0:14, lambdahat), col=2)
table(questions_orales) # nombre observé
round(n*dpois(0:15, lambdahat)) #nombre attendu
chisq.test(table(questions_orales), round(n*dpois(0:14, lambdahat))) # p-valeur large : on accepte H0 que les données sont Poisson

#Q2
#prior Gamma(2, 2)
par(mfrow=c(3, 1))
curve(dgamma(x, 2, 2), xlim=c(0, 4), main="Prior", ylab="density")
curve(dgamma(x, 2+qtot, 2+n), xlim=c(0, 4), main="Posterior model 1", ylab="density")
curve(dgamma(x, 2+qh, 2+nh), xlim=c(0, 4), main= "Posterior model 2", ylab="density")
curve(dgamma(x, 2+qf, 2+nf), col=2, add=T)
legend("topright", c("H", "F"), col=1:2, lty=1)

#Q3
# Intervalles de crédibilité
qgamma(c(.025, .975), 2+qtot, 2+n) # modèle 1
qgamma(c(.025, .975), 2+qh, 2+nh) # modèle 2 - hommes
qgamma(c(.025, .975), 2+qf, 2+nf) # modèle 2 - femmes

# Q4
# r=lambda1/lambda2
par(mfrow=c(1, 1))
niter = 10000
lambda1 = rgamma(niter, 2+qh, 2+nh)
lambda2 = rgamma(niter, 2+qf, 2+nf)
r = lambda1 / lambda2
mean(r)
sd(r)
hist(r)
quantile(r, c(0.025, 0.975))
# Convergence de notre estimateur
plot(1:niter, cumsum(r)/(1:niter), type="l")

#Q5: Monte Carlo standard
# Fonctions pour la vraisemblance
lkd.model1 = function(y, n, lambda) {
  return(exp(-n * lambda + y * log(lambda)))
}

lkd.model2 = function(y1, n1, y2, n2, lambda1, lambda2) {
  return(exp(-n1 * lambda1 + y1 * log(lambda1) - n2 * lambda2 + y2 * log(lambda2)))
  #return(lambda1^y1*exp(-n1*lambda1)*lambda2^y2*exp(-n2*lambda2))
}

BF_MC = function(a, b, y1, n1, y2, n2, M) {
  lambda1 = rgamma(M, a, b)
  m1 = sum(lkd.model1(y1 + y2, n1 + n2, lambda1)) / M
  lambda2.1 = rgamma(M, a, b)
  lambda2.2 = rgamma(M, a, b)
  m2 = sum(lkd.model2(y1, n1, y2, n2, lambda2.1, lambda2.2)) / M
  return(m2 / m1)
}


BF_MC = function(a, b, y1, n1, y2, n2, M) {
  lambda1 = rgamma(M, a, b)
  m1 = cumsum(lkd.model1(y1 + y2, n1 + n2, lambda1)) / (1:M)
  lambda2.1 = rgamma(M, a, b)
  lambda2.2 = rgamma(M, a, b)
  m2 = cumsum(lkd.model2(y1, n1, y2, n2, lambda2.1, lambda2.2)) / (1:M)
  return(m2 / m1)
}

M = 1e6
#a=2; b=2; y1=qh; n1=nh; y2=qf; n2=nf
resMC = BF_MC(2, 2, qh, nh, qf, nf, M)
resMC[M]
plot(100:M, resMC[100:M], type="l")
abline(h=trueBF, col=2)

#Q6: Importance sampling
BF_IS = function(a, b, y1, n1, y2, n2, M){
  mean1 = (a+y1+y2)/(b+n1+n2)
  mean2.1 = (a+y1)/(b+n1)
  mean2.2 = (a+y2)/(b+n2)
  
  sigma1 = sqrt((a+y1+y2)/(b+n1+n2)^2)
  sigma2.1 = sqrt((a+y1)/(b+n1)^2)
  sigma2.2 = sqrt((a+y2)/(b+n2)^2)
  
  lambda1 = rnorm(M, mean1, sigma1)
  m1 = cumsum( lkd.model1(y1+y2, n1+n2, lambda1) * 
                 dgamma(lambda1, a, b) / dnorm(lambda1, mean1, sigma1))/(1:M)
  
  lambda2.1 = rnorm(M, mean2.1, sigma2.1)
  lambda2.2 = rnorm(M, mean2.2, sigma2.2)
  m2 = cumsum(lkd.model2(y1, n1, y2, n2, lambda2.1, lambda2.2) * 
                dgamma(lambda2.1, a, b) * dgamma(lambda2.2, a, b) / 
                (dnorm(lambda2.1, mean2.1, sigma2.1) * 
                   dnorm(lambda2.2, mean2.2, sigma2.2))) / (1:M)
  
  return(m2/m1)
}

resIS = BF_IS(2, 2, qh, nh, qf, nf, M)
resIS[M]
plot(100:M, resIS[100:M], type="l")
abline(h=trueBF, col=2)

#Notons que l'estimateur par Importance Sampling converge bien plus vite que celui par Monte Carlo standard.

# Comparons le temps d'exécution :
system.time(BF_MC(2, 2, qh, nh, qf, nf, M))
system.time(BF_IS(2, 2, qh, nh, qf, nf, M))
# Importance sampling est un peu plus lent à nombre d'itérations égal, 
# donc bien plus efficace quand on prend en compte la vitesse de convergence.


#Q7 Valeur analytique
BF_analytical = function(a, b, y1, n1, y2, n2){
  m1 = b^a/gamma(a)*gamma(a+y1+y2)/(b+n1+n2)^(a+y1+y2)
  m2 = b^(2*a)/gamma(a)^2*gamma(a+y1)/(b+n1)^(a+y1)*gamma(a+y2)/(b+n2)^(a+y2)
  return(m2/m1)
}

# Ne marche pas parce que n est trop grand
BF_analytical(2, 2, qh, nh, qf, nf)

# Marche sur l'échelle log
BF_analytical2 = function(a, b, y1, n1, y2, n2){
  m1 = a*log(b)-lgamma(a)+lgamma(a+y1+y2)-(a+y1+y2)*log(b+n1+n2)
  m2 = 2*a*log(b)-2*lgamma(a)+lgamma(a+y1)-(a+y1)*log(b+n1)+lgamma(a+y2)-(a+y2)*log(b+n2)
  return(exp(m2-m1))
}
(trueBF = BF_analytical2(2, 2, qh, nh, qf, nf))
log10(trueBF) # évidence "substantielle" en faveur du modèle 1


#Q11
# la probabilité a posteriori du modèle 1 est
p = trueBF/(1+trueBF)

#Q12
nsim = 1e5
# cas a)
ech1 = rgamma(nsim, 2+qtot, 2+n)
mean(ech1)
sd(ech1)

# cas b)
ech2 = rgamma(nsim, 2+qf, 2+nf)
mean(ech2)
sd(ech2)

# cas c)
ech3 = rep(NA, nsim)
for(i in 1:nsim){
  modele = sample(c(1, 2), 1, prob=c(p, 1-p))
  if(modele==1){
    ech3[i] = rgamma(1, 2+qtot, 2+n)
  }
  else{
    ech3[i] = rgamma(1, 2+qf, 2+nf)
  }
}
mean(ech3)
sd(ech3)

par(mfrow=c(3, 1))
plot(density(ech1), xlim=c(.5, 4), main="Modèle 1")
plot(density(ech2), xlim=c(.5, 4), main="Modèle 2")
plot(density(ech3), xlim=c(.5, 4), main="Modèles pondérés")