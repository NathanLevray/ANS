#Nathan Marié L3 Info Groupe 2 Option ANS
rm(list=ls(all=TRUE)); 
graphics.off();

monRep="/home/nathan/Scolaritée/Licence Informatique/S6/ANS"
setwd(monRep)

#2 Jeu de données 

fichier = "datasets/notes.csv"
x = read.csv(fichier,header = T)
x

#2.2 Homogénéisation de la variable Anglais
x$Anglais = x$Anglais/100*20
x

#3 ACP
#3.1 Calcul sur le tableau de notes homogénéisé 

library(FactoMineR)
acp.res = PCA(x,scale.unit=F,graph=F,ncp=ncol(x))
acp.res

#3.2

acp.res$eig

barplot(acp.res$eig[,3],main="inertie expliquée cumulée",col="blue")

#3.3
acp.res$ind$coord

plot(acp.res,axes=c(1,2),choix="ind",title="Représentation des individus \n Plan des 2 premières composantes principales")

#3.4

acp.res$var$cor

plot(acp.res, axes=c(1,2),choix="varcor",title="Représentation des variables\n Corrélation entre les variables initiales et les 2 premières CP")



cor(x$Algo1,x$Algo2)
