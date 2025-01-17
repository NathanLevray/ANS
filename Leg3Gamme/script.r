#Nathan Marié L3 Info Groupe 2 Option ANS
rm(list=ls(all=TRUE)); 
graphics.off();

monRep="/home/nathan/Scolaritée/Licence Informatique/S6/ANS"
setwd(monRep)
#print(list.files())
#print(list.dirs())
#print(getwd())

# Importer les données
load("data/Leg3_GammeCorrectedData.FixeInclus.Rdata")
data= AOA.PFB.Fixe

summary(data)
class(data)
head(data)

# Les métadata

typage=sapply(data,class)
#print(typage)
indNumeric=which(typage=="numeric")
indNumeric=which(sapply(data,class)=="numeric")

df=data[,indNumeric]
class(df)
head(df)
summary(df)


dfNoClear = df #Une copie de df, car lorsque on va nettoyer df la maorité de paraFBD ne va plus être présent dans FBD

#for (i in 1:ncol(df)){
#  boxplot(df[,i],main=colnames(df)[i])
#}


#Nettoyer les données:
#############################

#On retire les colonnes constantes
col.constant <- sapply(df, FUN=function(x){sd(x,na.rm=T)==0})
#print(col.constant)
if (sum(col.constant)>0){
  cat("full constant column = ",colnames(df)[which(sum(col.constant)>0)])
  df=df[,!col.constant]
}


#On retire les élements vides sur les lignes
row.na = apply(df,1, function(x){any(is.na(x))})
if (sum(row.na)>0){
  cat("full NA line index= ",which(sum(row.na)>0))
  df = df[!row.na,]
  dfNoClear = dfNoClear[!row.na,]
}

boxplot(df,main="DF with no NA or constant")

#Corélation 
M= cor(df,use="pairwise.complete.obs")
corrplot::corrplot(M,tl.cex = 1)
boxplot(M,main="Correlation des valeurs",ylim=c(-1,1))

#Redondance
Rthreshod=0.75 #Modifiable 
idtoRemove=NULL
for (i in 1:ncol(df)){
  id=which(abs(M[i,])>Rthreshod)
  id=id[id>i]
  if(length(id)){
    idtoRemove=c(idtoRemove,id)
    #print(paste(colnames(df)[i],"correlated with", colnames(df)[id]))
  }
}
idtoRemove=unique(idtoRemove)
print(colnames(df)[idtoRemove])
df=df[,-idtoRemove] #retrait des colonnes redondantes

boxplot(df,main="Clean DF")
#On apperçoit que Course est toujours présente.
#On peut supposer que course n'a aucun rapport avec les autres données en effet le véhicule n'impacte que peu les algues.
df=subset(df, select = -Course)

#L'argumment peut aussi être fait avec Soundspeed ?
#Ou non peut-être que la vitesse du sound nous indique l'oppacité de l'eau ?? #Et donc la luminosité ????
df=df[ , !(names(df) %in% c("Soundspeed"))] #Equivalent à subset (df,select =-Soundspeed)

boxplot(df,main="Usable DF")


dfa = dfNoClear[,algue]
dfb = dfNoClear[,paraPFB]
dfc = dfNoClear[,c(algue,paraPFB)]
dfd = dfNoClear[,c(algue,paraPFB,"longitudeCorr","latitudeCorr")]

long = df$longitudeCorr
lat = df$latitudeCorr

save(df,dfa,dfb,dfc,dfd,long,lat,file="var.Rdata")#Sauvegarder les données nettoyées
#############################
load("var.Rdata")

#for (i in 1:ncol(df)){
#  boxplot(df[,i],main=colnames(df)[i])
#}

#Sur DF
#############################
#Matrice Similarité 
d=as.matrix(dist(scale(df))) #Scale de DF pour éviter problème d'échelle entre le x et y
similarite=1/(1+d)
#heatmap(similarite) #Trop de données
dd=as.dist(1-similarite)

#Hclust:
############
HC=hclust(dd,method="ward.D2")

#Silhouette:
########
vSil=NULL
vK=2:(300)
for (k in vK){
  sil=cluster::silhouette(cutree(HC,k),dd)
  vSil=c(vSil,mean(sil[,"sil_width"]))
}
plot(vK,vSil, main="silhouette selon K") 
KoptHierachic=vK[which.max(vSil)] #Permet de déterminer un K pour avoir des meilleurs clusters
print(KoptHierachic)
K=KoptHierachic

res.HC=cutree(HC,K)
plot(df,col=res.HC,main="Hierarchic clustering")
#Le Hierarchic clustering semble inutile avec ces données

#kmeans:
############

#Choix de K
KoptKmeanslist = list()
for(i in 1:100){
  KoptKmeanslist[[i]] = kmeans(df,i)  
}

betweenss_totss = list()
for (i in 1:100){
  betweenss_totss[[i]] = KoptKmeanslist[[i]]$betweenss/KoptKmeanslist[[i]]$totss
}

plot(1:100, betweenss_totss,ylab="Btw SS / Total SS",xlab ="Clusters K")

for (i in 2:10){
  x11();
  plot(df,col=KoptKmeanslist[[i]]$cluster,main=paste("kmedoid K=",i))
}

#A l'oeil il ne semble y avoir aucun K qui donne des clusters intéressant. #Peut-êre 4 ou 5
#K = Koptmeans = 4

res.KM=kmeans(scale(df),centers=3)
plot(df,col=res.KM$cluster, main="kmeans");
points(res.KM$centers,pch='+',cex=9)

#library("cluster")
#res=pam(1-similarite,K)
#res.PAM=res$clustering
#plot(df,col=res.PAM,main="PAM")

#res.PAM=pam(df,k=3)
#plot(df,col=res.PAM$clustering, main="kmedoid");
#points(res.PAM$medoids,pch='+',cex=6)


#Spectral Clustering:
########
# projection des données dans un autre espace -> laplacien=D^-1/2 Similarite D^-1/2
# espace des vecteurs propres de ce laplacien
sigma=0.05 #a modifier selon la dispersion des données
similarite=exp(-d*d/sigma)
diag(similarite)=0;
degre=rowSums(similarite)
ds <- degre^(-0.5)
Ds=diag(ds)
laplacien=Ds %*% similarite %*% Ds
eig=eigen(laplacien,symmetric=TRUE)
#KoptSpectral=which.max(eig$values) #Problème il y 1 cluster, donc le spectral clustering est inutile
plot(eig$values)
#Il y a plusieurs valeurs max
KoptSpectral = which(eig$values == max(eig$values))
K = tail(KoptSpectral,1)
print(K)

Z=eig$vectors[,1:K];
Zn=Z/apply(Z,MARGIN=1,FUN=function(x) norm(matrix(x),"f"));
cl=kmeans(Zn,centers=K,iter.max = 100, nstart = 10);
res.SC=cl$cluster;
plot(df,col=res.SC, main="spectral clustering");
plot(df$latitudeCorr,df$Cryptophyta,col=res.SC) #Ne montre rien 
#Spectral clusters ne montre rien de particulier


#Conclusion sur df il y a trop de données pour pouvoir effectuer un apprentisage
#Solution limiter les données

#Sur dfa (les algues)
#############################
#Matrice Similarité 
d=as.matrix(dist(scale(dfa))) 
similarite=1/(1+d)
dd=as.dist(1-similarite)

HC=hclust(dd,method="ward.D2")

#Silhouette:
########
#vSil=NULL
#vK=2:(300)
#for (k in vK){
#  sil=cluster::silhouette(cutree(HC,k),dd)
#  vSil=c(vSil,mean(sil[,"sil_width"]))
#}
#plot(vK,vSil, main="silhouette selon K") 
#KoptHierachic=vK[which.max(vSil)] #Permet de déterminer un K pour avoir des meilleurs clusters
#print(KoptHierachic)
#K=KoptHierachic
#Silhouette ne fonctionne pas le K augmente toujours

# Choix de K 

 k = list()
 for(i in 1:20){
   k[[i]] = kmeans(dfa,i)  
 }
 
 betweenss_totss = list()
 for (i in 1:20){
   betweenss_totss[[i]] = k[[i]]$betweenss/k[[i]]$totss
 }
 
 plot(1:20, betweenss_totss,type="b",ylab="Btw SS / Total SS",xlab ="Clusters K")

# Pour K = 5, on a la pente qui devient constante
K=5

res.HC=cutree(HC,K)
plot(long, lat, col=res.HC, pch = as.character(res.HC))

label = NULL
for (lab in unique(res.HC)){
  label = cbind(label, res.HC == lab)
}
dfa_lab = cbind(dfa, label)

Mcorr = cor(dfa_lab, use = "pairwise.complete.obs")
corrplot::corrplot(Mcorr,tl.cex = 1)
#On regarde les valeurs absolues
#Label 1 : Corrèle principalement avec Diatoms et Cryptophyta et légèrement avec le reste
#Label 2 : Corrèle principalement avec BlueGreen et Cryptophyta et légèrement avec Green 
#Label 3 : Corrèle principalement avec Diatoms
#Label 4 : Corrèle principalement avec Diatoms et légèrement Green
#Label 5 : Corrèle principalement avec Green et avec Cryptophyta et Diatoms


#Sur dfb (les paramètres physiques) 
#############################
#Matrice Similarité 
d=as.matrix(dist(scale(dfb))) 
similarite=1/(1+d)
dd=as.dist(1-similarite)

HC=hclust(dd,method="ward.D2")

#Silhouette:
########
#vSil=NULL
#vK=2:(300)
#for (k in vK){
#  sil=cluster::silhouette(cutree(HC,k),dd)
#  vSil=c(vSil,mean(sil[,"sil_width"]))
#}
#plot(vK,vSil, main="silhouette selon K") 
#KoptHierachic=vK[which.max(vSil)] #Permet de déterminer un K pour avoir des meilleurs clusters
#print(KoptHierachic)
#K=KoptHierachic
#Silhouette ne fonctionne pas le K augmente toujours, comme dans dfa

# Choix de K 

k = list()
for(i in 1:20){
  k[[i]] = kmeans(dfa,i)  
}

betweenss_totss = list()
for (i in 1:20){
  betweenss_totss[[i]] = k[[i]]$betweenss/k[[i]]$totss
}

plot(1:20, betweenss_totss,type="b",ylab="Btw SS / Total SS",xlab ="Clusters K")

# Pour K = 6 ?, on a la pente qui devient constante
K=6
res.HC=cutree(HC,K)
plot(long, lat, col=res.HC, pch = as.character(res.HC))

label = NULL
for (lab in unique(res.HC)){
  label = cbind(label, res.HC == lab)
}
dfb_lab = cbind(dfb, label)

Mcorr = cor(dfb_lab, use = "pairwise.complete.obs")
corrplot::corrplot(Mcorr,tl.cex = 1)
#Label 1 : Corrèle principalement avec Oxygène et Saturation et légèrement avec Température
#Label 2 : Corrèle principalement avec Température et avec Salinité
#Label 3 : Corrèle principalement avec Température et légèrement avec les autres 
#(Remarque: Le nombre de cluster est peut-être trop élévé ?)
#Label 4 : Corrèle principalement avec Oxygène et Saturation et un peu avec Salinité
#Label 5 : Corrèle principalement avec Température
#Label 6 : Corrèle principalement avec Salinité

#Sur dfc (les algues et les paramètres physiques)
#############################
#Matrice Similarité 
d=as.matrix(dist(scale(dfc))) 
similarite=1/(1+d)
dd=as.dist(1-similarite)


HC=hclust(dd,method="ward.D2")
#Silhouette non concluant, il donne K = 223 ce qui rend compliqué le fait de trouver les facteurs qui déterminent ces 223 groupes

# Choix de K 

k = list()
for(i in 1:20){
  k[[i]] = kmeans(dfa,i)  
}

betweenss_totss = list()
for (i in 1:20){
  betweenss_totss[[i]] = k[[i]]$betweenss/k[[i]]$totss
}

plot(1:20, betweenss_totss,type="b",ylab="Btw SS / Total SS",xlab ="Clusters K")

# Pour K = 6 ou 5, on a la pente qui devient constante
#La méthode est à revoir, elle est imprécise
K=6
res.HC=cutree(HC,K)
plot(long, lat, col=res.HC, pch = as.character(res.HC))

label = NULL
for (lab in unique(res.HC)){
  label = cbind(label, res.HC == lab)
}
dfc_lab = cbind(dfc, label)

Mcorr = cor(dfc_lab, use = "pairwise.complete.obs")
corrplot::corrplot(Mcorr,tl.cex = 1)

#Label 1 : Corrèle principalement avec Oxygène et Saturation et Salinité et légèrement avec Cryptophyta
#Label 2 : Corrèle principalement avec Diatoms et avec Cryptophyta et température
#Label 3 : Corrèle avec Diatoms, Cryptophyta, Température et Salinité
#Label 4 : Corrèle principalement avec Température et légèrement avec Diatoms et Salinité
#Label 5 : Corrèle principalement avec Oxygen et Saturation (Oxygen et Saturation sont fortement corrèlé entre eux) et moyennement avec Diatoms et Salinité
#Label 6 : Corrèle Fortement avec Green et moyennement avec Diatoms, Cryptophyta et Température

#Sur dfd (les algues et les paramètres physiques et la longitude latitude)
#############################
#Matrice Similarité 
d=as.matrix(dist(scale(dfc))) 
similarite=1/(1+d)
dd=as.dist(1-similarite)

HC=hclust(dd,method="ward.D2")
#Silhouette non concluant comme avec les autres, il donne K = 223

# Choix de K 

k = list()
for(i in 1:20){
  k[[i]] = kmeans(dfa,i)  
}

betweenss_totss = list()
for (i in 1:20){
  betweenss_totss[[i]] = k[[i]]$betweenss/k[[i]]$totss
}

plot(1:20, betweenss_totss,type="b",ylab="Btw SS / Total SS",xlab ="Clusters K")

# Pour K = 5, on a la pente qui devient constante
K=5


res.HC=cutree(HC,K)
plot(long, lat, col=res.HC, pch = as.character(res.HC))

label = NULL
for (lab in unique(res.HC)){
  label = cbind(label, res.HC == lab)
}
dfd_lab = cbind(dfd, label)

Mcorr = cor(dfd_lab, use = "pairwise.complete.obs")
corrplot::corrplot(Mcorr,tl.cex = 1)

#Label 1 : Corrèle principalement avec Cryptophyta, Oxygen, Saturation et légèrement en dessus des autres dans latitude
#Label 2 : Corrèle faublement dans Diatoms, Cryptophyta, Température et Salinité.
#Label 3 : COrrèle principalement dans Température avec un peu en Diatoms, salinité et latitude
#Label 4 : Corrèle principalement dans Oxygen et Saturation et un peu moins dans diatoms, salinité et latitude 
#Label 5 : Corrèle Fortement dans Green et moyennement dans Diatoms, Cryptophyta, température et longitude
#Remarque le label 5 est le seul qui corrèle avec longitude, en effet on voit que Green corrèle est celui qui corrèle le plus avec longitude
