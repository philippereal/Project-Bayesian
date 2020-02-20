# Feuille 1, exercice 2
rm(list = ls())

# Question 1 : loi beta
par(mfrow = c(3, 2))
curve(dbeta(x, 1, 1)) # la Beta(1,1) correspond à la loi uniforme
curve(dbeta(x, 2, 1))
curve(dbeta(x, 2, 4))
curve(dbeta(x, 43, 67)) # avec des valeurs plus élevées des paramètres, on a une loi plus piquée
curve(dbeta(x, 0.1, 2))
curve(dbeta(x, 0.1, 0.1))

# variante : toutes les courbes sur le même graphe
par(mfrow = c(1, 1))
curve(dbeta(x, 1, 1), ylim = c(0, 8)) # la Beta(1,1) correspond à la loi uniforme
curve(dbeta(x, 2, 1), add = T, col = 2)
curve(dbeta(x, 2, 4), add = T, col = 3)
curve(dbeta(x, 43, 67), add = T, col = 4) # avec des valeurs plus élevées des paramètres, on a une loi plus piquée
curve(dbeta(x, .1, 2), add = T, col = 5)
curve(dbeta(x, .1, .1), add = T, col = 6)

# variante : même graphique avec ggplot2
require(ggplot2)
x = seq(0, 1, len = 1e4)
p = qplot(x, geom = "blank")
p = p + stat_function(
  aes(x = x, y = ..y..),
  fun = dbeta,
  colour = "red",
  args = list(shape1 = 1, shape2 = 1)
)
p = p + stat_function(
  aes(x = x, y = ..y..),
  fun = dbeta,
  colour = "green",
  args = list(shape1 = 2, shape2 = 1)
)
p = p + stat_function(
  aes(x = x, y = ..y..),
  fun = dbeta,
  colour = "blue",
  args = list(shape1 = 2, shape2 = 4)
)
p = p + stat_function(
  aes(x = x, y = ..y..),
  fun = dbeta,
  colour = "black",
  args = list(shape1 = 43, shape2 = 67)
)
p = p + stat_function(
  aes(x = x, y = ..y..),
  fun = dbeta,
  colour = "orange",
  args = list(shape1 = .1, shape2 = 2)
)
p = p + stat_function(
  aes(x = x, y = ..y..),
  fun = dbeta,
  colour = "purple",
  args = list(shape1 = .1, shape2 = .1)
)
p
# à chaque fois, le paramètre args donne les paramètres à envoyer à la fonction dbeta

# Question 3
# Simulons quatre échantillons
theta = 0.6
y1 = rbinom(20, 1, theta)
y2 = rbinom(100, 1, theta)
y3 = rbinom(1000, 1, theta)
y4 = rbinom(1e4, 1, theta)

par(mfrow = c(3, 2))
curve(dbeta(x, 1, 1))
curve(dbeta(x, 1 + sum(y1), 1 + 20 - sum(y1)))
curve(dbeta(x, 1 + sum(y2), 1 + 100 - sum(y2)))
curve(dbeta(x, 1 + sum(y3), 1 + 1000 - sum(y3)))
curve(dbeta(x, 1 + sum(y4), 1 + 10000 - sum(y4)))
# lorsque la taille de l'échantillon augmente, la posterior se concentre autour de la vraie valeur

# variante :  même graphique avec ggplot2
x = seq(0, 1, len = 1e4)
p = qplot(x, geom = "blank")
p = p + stat_function(
  aes(x = x, y = ..y..),
  fun = dbeta,
  colour = "red",
  args = list(shape1 = 1, shape2 = 1)
)
p = p + stat_function(
  aes(x = x, y = ..y..),
  fun = dbeta,
  colour = "green",
  args = list(shape1 = 1 + sum(y1), shape2 = 1 + 20 - sum(y1))
)
p = p + stat_function(
  aes(x = x, y = ..y..),
  fun = dbeta,
  colour = "blue",
  args = list(shape1 = 1 + sum(y2), shape2 = 1 + 100 - sum(y2))
)
p = p + stat_function(
  aes(x = x, y = ..y..),
  fun = dbeta,
  colour = "black",
  args = list(shape1 = 1 + sum(y3), shape2 = 1 + 1000 - sum(y3))
)
p = p + stat_function(
  aes(x = x, y = ..y..),
  fun = dbeta,
  colour = "orange",
  args = list(shape1 = 1 + sum(y4), shape2 = 1 + 1e4 - sum(y4))
)
p

#avec la vraie valeur en plus
p = p + geom_vline(xintercept = theta, col = "purple")
p

p + xlim(c(.4, .8))


# Feuille 1, exercice 3

rm(list=ls())

# question 3

par(mfrow = c(2, 2))

# Lois a priori
curve(dgamma(x, 1, 0.5), xlim = c(0, 3), ylim = c(0, 4))
curve(dgamma(x, 1, 2), add = T, col = 2)
curve(dgamma(x, 1, 10), add = T, col = 3)
curve(dgamma(x, 2, 2), add = T, col = 4)

# question 4

lambda = 0.5
y1 = rpois(20, lambda)
y2 = rpois(100, lambda)
y3 = rpois(1000, lambda)

# Loi a posteriori après 20 observations
curve(dgamma(x, 1 + sum(y1), 0.5 + 20),
      xlim = c(0, 3),
      ylim = c(0, 4))
curve(dgamma(x, 1 + sum(y1), 2 + 20), add = T, col = 2)
curve(dgamma(x, 1 + sum(y1), 10 + 20), add = T, col = 3)
curve(dgamma(x, 2 + sum(y1), 2 + 20), add = T, col = 4)
curve(dgamma(x, .5 + sum(y1), 20), add = T, col = 5)

# Après 100 observations
curve(dgamma(x, 1 + sum(y2), 0.5 + 100),
      xlim = c(0, 3),
      ylim = c(0, 4))
curve(dgamma(x, 1 + sum(y2), 2 + 100), add = T, col = 2)
curve(dgamma(x, 1 + sum(y2), 10 + 100), add = T, col = 3)
curve(dgamma(x, 2 + sum(y2), 2 + 100), add = T, col = 4)
curve(dgamma(x, .5 + sum(y2), 100), add = T, col = 5)

# Après 1000 observations
curve(dgamma(x, 1 + sum(y3), 0.5 + 1000),
      xlim = c(0, 3),
      ylim = c(0, 10))
curve(dgamma(x, 1 + sum(y3), 2 + 1000), add = T, col = 2)
curve(dgamma(x, 1 + sum(y3), 10 + 1000), add = T, col = 3)
curve(dgamma(x, 2 + sum(y3), 2 + 1000), add = T, col = 4)
curve(dgamma(x, .5 + sum(y3), 1000), add = T, col = 5)


# question 5
# pour une loi Gamma(a, b), l'espérance vaut a/b et le mode (a-1)/b

# espérance a priori :
c(1, 1, 1, 2) / c(0.5, 2, 10, 2)

# espérance a posteriori :
(c(1, 1, 1, 2) + sum(y1)) / (c(0.5, 2, 10, 2) +  20)
(c(1, 1, 1, 2) + sum(y2)) / (c(0.5, 2, 10, 2) +  100)
(c(1, 1, 1, 2) + sum(y3)) / (c(0.5, 2, 10, 2) +  1000)

# MAP
(c(1, 1, 1, 2) + sum(y1) - 1) / (c(0.5, 2, 10, 2) +  20)
(c(1, 1, 1, 2) + sum(y2) - 1) / (c(0.5, 2, 10, 2) +  100)
(c(1, 1, 1, 2) + sum(y3) - 1) / (c(0.5, 2, 10, 2) +  1000)

# si on n'a pas de formule analytique, on peut tout de même trouver le mode
optimize(function(x) {
  dgamma(x, 1 + sum(y3), 1 + 1000)
}, c(0, 3), maximum = T)
# pour l'espérance sans formule analytique, on pourra faire du Monte-Carlo, cf TD 2

# question 6
# intervalles de crédibilité pour une prior gamma(1, 1)
qgamma(c(.025, .975), 1 + sum(y1), .5 + 20)
qgamma(c(.025, .975), 1 + sum(y2), .5 + 100)
qgamma(c(.025, .975), 1 + sum(y3), .5 + 1000)

# intervalles de crédibilité pour différentes prior, avec n=20
qgamma(c(.025, .975), 1 + sum(y1), .5 + 20)
qgamma(c(.025, .975), 1 + sum(y1), 2 + 20)
qgamma(c(.025, .975), 1 + sum(y1), 10 + 20)
qgamma(c(.025, .975), 2 + sum(y1), .2 + 20)

# intervalles de crédibilité pour différentes prior, avec n=1000
qgamma(c(.025, .975), 1 + sum(y3), .5 + 1000)
qgamma(c(.025, .975), 1 + sum(y3), 2 + 1000)
qgamma(c(.025, .975), 1 + sum(y3), 10 + 1000)
qgamma(c(.025, .975), 2 + sum(y3), .2 + 1000)
