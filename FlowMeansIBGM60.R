#################################################################################
#flowMeans óptimo (número de clusters prefijado) para el conjunto IBGM
################################################################################
library(flowCore)
library(flowMeans)
library(dplyr)
library(clue)
library(boot)

directorio<-"E:/datos/datos_ibgm/33_files_concat_especifico.csv"
datos<-read.csv(directorio)
tinic<-Sys.time()
salida<- flowMeans(datos[,c(8:19,21:34, 36:47)], NumC = 60, Standardize = FALSE)
tfin<-Sys.time()
tiempo<-tfin-tinic
etiquetas_ibgm<-data.frame(Manual=datos$OmiqFilter, Clustering=salida@Label)

#Calcular matrices F_socre, precision y recall
n_cells<-length(unique(etiquetas_ibgm[,"Manual"]))
n_clusters<-length(unique(etiquetas_ibgm[,"Clustering"]))
pr_mat <- re_mat <- F_mat <- matrix(NA, nrow = n_clusters, ncol = n_cells)
resultados<-data.frame(matrix(rep(0,8),nrow = 2, ncol = 4))
colnames(resultados)<-c("Pr","Re","F","NCelulas")
rownames(resultados)<-c("Maximo", "Hungaro")
for (i in 1:n_clusters) {
  for (j in 1:n_cells) {
    tp<-sum(etiquetas_ibgm[,"Clustering"]==i & etiquetas_ibgm[,"Manual"]==unique(etiquetas_ibgm[,"Manual"])[j])
    fp<-sum(etiquetas_ibgm[,"Clustering"]==i & etiquetas_ibgm[,"Manual"]!=unique(etiquetas_ibgm[,"Manual"])[j])
    fn<-sum(etiquetas_ibgm[,"Clustering"]!=i & etiquetas_ibgm[,"Manual"]==unique(etiquetas_ibgm[,"Manual"])[j])
    if(tp==0){
      pr_mat[i,j]<-0
      re_mat[i,j]<-0
    }else{
      pr_mat[i,j]<-tp/(tp+fp)
      re_mat[i,j]<-tp/(tp+fn)
    }
    if(pr_mat[i,j]==0 & re_mat[i,j] ==0){
      F_mat[i,j]=0
    }else{
      F_mat[i,j]<-(2*pr_mat[i,j]*re_mat[i,j])/(pr_mat[i,j]+re_mat[i,j])
    } 
  }
}

#Resolver con algoritmo del húngaro
if(n_clusters<n_cells){
  #Las poblaciones sin cluster se les asignará NA y tendrán valor 0 para pr, re y F
  asignacion<-solve_LSAP(F_mat, maximum = TRUE)
  cat("Numero clusters menor que el numero de grupos celulares\n")
  for (i in 1:n_clusters) {
    j_cell<-asignacion[i]
    resultados["Hungaro","Pr"]<-resultados["Hungaro","Pr"]+pr_mat[i,j_cell]
    resultados["Hungaro","Re"]<-resultados["Hungaro","Re"]+re_mat[i,j_cell]
    resultados["Hungaro","F"]<-resultados["Hungaro","F"]+F_mat[i,j_cell]
  }
  resultados["Hungaro","Pr"]<-resultados["Hungaro","Pr"]/n_cells
  resultados["Hungaro","Re"]<-resultados["Hungaro","Re"]/n_cells
  resultados["Hungaro","F"]<-resultados["Hungaro","F"]/n_cells
  resultados["Hungaro","NCelulas"]<-n_cells
}else{
  #Trabajo con la matriz traspuesta
  cat("Numero clusters mayor o igual que el numero de grupos celulares\n")
  asignacion<-solve_LSAP(t(F_mat), maximum = TRUE)
  for (j in 1:n_cells) {
    i_clus<-asignacion[j]
    resultados["Hungaro","Pr"]<-resultados["Hungaro","Pr"]+pr_mat[i_clus,j]
    resultados["Hungaro","Re"]<-resultados["Hungaro","Re"]+re_mat[i_clus,j]
    resultados["Hungaro","F"]<-resultados["Hungaro","F"]+F_mat[i_clus,j]
  }
  resultados["Hungaro","Pr"]<-resultados["Hungaro","Pr"]/n_cells
  resultados["Hungaro","Re"]<-resultados["Hungaro","Re"]/n_cells
  resultados["Hungaro","F"]<-resultados["Hungaro","F"]/n_cells
  resultados["Hungaro","NCelulas"]<-n_cells
}
#Resolución según el maximo
asignacion_maximo<-vector("list",length = n_cells)
for (j in 1:n_cells) {
  i_max<-which.max(F_mat[,j])
  asignacion_maximo[[j]]<-i_max
  resultados["Maximo","Pr"]<-resultados["Maximo","Pr"]+sum(etiquetas_ibgm[,"Manual"]==unique(etiquetas_ibgm[,"Manual"])[j])*pr_mat[i_max,j]
  resultados["Maximo","Re"]<-resultados["Maximo","Re"]+sum(etiquetas_ibgm[,"Manual"]==unique(etiquetas_ibgm[,"Manual"])[j])*re_mat[i_max,j]
  resultados["Maximo","F"]<-resultados["Maximo","F"]+sum(etiquetas_ibgm[,"Manual"]==unique(etiquetas_ibgm[,"Manual"])[j])*F_mat[i_max,j]
  resultados["Maximo","NCelulas"]<-resultados["Maximo","NCelulas"]+sum(etiquetas_ibgm[,"Manual"]==unique(etiquetas_ibgm[,"Manual"])[j])
}
resultados["Maximo","Pr"]<-resultados["Maximo","Pr"]/resultados["Maximo","NCelulas"]
resultados["Maximo","Re"]<-resultados["Maximo","Re"]/resultados["Maximo","NCelulas"]
resultados["Maximo","F"]<-resultados["Maximo","F"]/resultados["Maximo","NCelulas"]
resultados

################################################################################
#IC con cheap bootstrap
n_boot<-10

funcion_boot<-function(datos, indices){
  temp<-datos[indices,]
  n_cells<-length(unique(temp[,"Manual"]))
  n_clusters<-length(unique(temp[,"Clustering"]))
  pr_mat <- re_mat <- F_mat <- matrix(NA, nrow = n_clusters, ncol = n_cells)
  resultados<-rep(0,10)
  names(resultados)<-c("Pr_hungaro","Re_hungaro","F_hungaro","NCelulas_hungaro","Accuracy_hungaro","Pr_maximo","Re_maximo","F_maximo","NCelulas_maximo","Accuracy_maximo")
  for (i in 1:n_clusters) {
    for (j in 1:n_cells) {
      tp<-sum(temp[,"Clustering"]==i & temp[,"Manual"]==unique(temp[,"Manual"])[j])
      fp<-sum(temp[,"Clustering"]==i & temp[,"Manual"]!=unique(temp[,"Manual"])[j])
      fn<-sum(temp[,"Clustering"]!=i & temp[,"Manual"]==unique(temp[,"Manual"])[j])
      if(tp==0){
        pr_mat[i,j]<-0
        re_mat[i,j]<-0
      }else{
        pr_mat[i,j]<-tp/(tp+fp)
        re_mat[i,j]<-tp/(tp+fn)
      }
      if(pr_mat[i,j]==0 & re_mat[i,j] ==0){
        F_mat[i,j]=0
      }else{
        F_mat[i,j]<-(2*pr_mat[i,j]*re_mat[i,j])/(pr_mat[i,j]+re_mat[i,j])
      } 
    }
  }
  
  #Resuleve con el algoritmo del húngaro
  if(n_clusters<n_cells){
    #Las poblaciones sin cluster se les asignará NA y tendrán valor 0 para pr, re y F
    asignacion<-solve_LSAP(F_mat, maximum = TRUE)
    cat("Numero clusters menor que el numero de grupos celulares\n")
    asignacion_inversa<-rep(0,n_cells)
    for (i in 1:n_clusters) {
      j_cell<-asignacion[i]
      asignacion_inversa[j_cell]<-i
      resultados["Pr_hungaro"]<-resultados["Pr_hungaro"]+pr_mat[i,j_cell]
      resultados["Re_hungaro"]<-resultados["Re_hungaro"]+re_mat[i,j_cell]
      resultados["F_hungaro"]<-resultados["F_hungaro"]+F_mat[i,j_cell]
    }
    resultados["Accuracy_hungaro"]<-mean(temp$Clustering==asignacion_inversa[match(temp$Manual, unique(temp$Manual))])
    resultados["Pr_hungaro"]<-resultados["Pr_hungaro"]/n_cells
    resultados["Re_hungaro"]<-resultados["Re_hungaro"]/n_cells
    resultados["F_hungaro"]<-resultados["F_hungaro"]/n_cells
    resultados["NCelulas_hungaro"]<-n_cells
  }else{
    #Trabajo con la matriz traspuesta
    cat("Numero clusters mayor o igual que el numero de grupos celulares\n")
    asignacion<-solve_LSAP(t(F_mat), maximum = TRUE)
    for (j in 1:n_cells) {
      i_clus<-asignacion[j]
      resultados["Pr_hungaro"]<-resultados["Pr_hungaro"]+pr_mat[i_clus,j]
      resultados["Re_hungaro"]<-resultados["Re_hungaro"]+re_mat[i_clus,j]
      resultados["F_hungaro"]<-resultados["F_hungaro"]+F_mat[i_clus,j]
    }
    resultados["Accuracy_hungaro"]<-mean(temp$Clustering==asignacion[match(temp$Manual, unique(temp$Manual))])
    resultados["Pr_hungaro"]<-resultados["Pr_hungaro"]/n_cells
    resultados["Re_hungaro"]<-resultados["Re_hungaro"]/n_cells
    resultados["F_hungaro"]<-resultados["F_hungaro"]/n_cells
    resultados["NCelulas_hungaro"]<-n_cells
  }
  #Resolución con el máximo
  asignacion_maximo<-rep(0, n_cells)
  for (j in 1:n_cells) {
    i_max<-which.max(F_mat[,j])
    asignacion_maximo[j]<-i_max
    resultados["Pr_maximo"]<-resultados["Pr_maximo"]+sum(temp[,"Manual"]==unique(temp[,"Manual"])[j])*pr_mat[i_max,j]
    resultados["Re_maximo"]<-resultados["Re_maximo"]+sum(temp[,"Manual"]==unique(temp[,"Manual"])[j])*re_mat[i_max,j]
    resultados["F_maximo"]<-resultados["F_maximo"]+sum(temp[,"Manual"]==unique(temp[,"Manual"])[j])*F_mat[i_max,j]
    resultados["NCelulas_maximo"]<-resultados["NCelulas_maximo"]+sum(temp[,"Manual"]==unique(temp[,"Manual"])[j])
  }
  resultados["Accuracy_maximo"]<-mean(temp$Clustering==asignacion_maximo[match(temp$Manual, unique(temp$Manual))])
  resultados["Pr_maximo"]<-resultados["Pr_maximo"]/resultados["NCelulas_maximo"]
  resultados["Re_maximo"]<-resultados["Re_maximo"]/resultados["NCelulas_maximo"]
  resultados["F_maximo"]<-resultados["F_maximo"]/resultados["NCelulas_maximo"]
  return(resultados)
}

boot_etiquetas<-boot(etiquetas_ibgm,statistic = funcion_boot, R=n_boot)

df<-data.frame(Algoritmo="FlowMeans_Optimo",
               Asignacion=c(rep("Maximo",100),rep("Hungaro",100)),
               F_score=c(boot_etiquetas$t[,8],boot_etiquetas$t[,3]),
               Accuracy=c(boot_etiquetas$t[,10],boot_etiquetas$t[,5]),
               Pr=c(boot_etiquetas$t[,6],boot_etiquetas$t[,1]),
               Re=c(boot_etiquetas$t[,7],boot_etiquetas$t[,2]))
write.table(df, file = "boxplot_bgm_flowmeansopt.csv",
            sep = "\t", row.names = F)