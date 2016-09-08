library(esd)

it <- month.abb[5:9]

load('pr.aaca.rda')
load('t2m.aaca.rda')

Y <- aggregate(subset(pr.aaca,it=it),year,'wetmean')
Y <- PCA(pcafill(Y))
X <- aggregate(subset(t2m.aaca,it=it),year,'mean')
X <- PCA(pcafill(X))

cca <- CCA(X,Y)
plot(cca)



