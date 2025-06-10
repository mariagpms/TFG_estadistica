#################################################################################
#FlowSOM óptimo (número de clusters prefijado) para el conjunto IBGM, para 10 semillas diferentes
################################################################################
library(flowCore)
library(flowMeans)
library(dplyr)
library(clue)
library(boot)

directorio<-"E:/datos/33_files_concat_especifico.fcs"
datos<-read.FCS(directorio, transformation = FALSE, truncate_max_range = FALSE)

semillas<-c(1,123,1014,141079,1410,994,100,1001,29,59)
tiempo<-vector(mode="list",length = length(semillas))
etiquetas_ibgm<-vector(mode="list",length = length(semillas))

n_cells<-vector(length = length(semillas))
n_clusters<-vector(length = length(semillas))
pr_mat<-vector(mode="list",length = length(semillas))
re_mat<-vector(mode="list",length = length(semillas))
F_mat<-vector(mode="list",length = length(semillas))
resultados<-vector(mode="list",length = length(semillas))
asignacion<-vector(mode="list",length = length(semillas))
asignacion_maximo<-vector(mode="list",length = length(semillas))
#Mismo proceso que en FlowSOMIBGMdefault pero para las 10 semillas 
for (s in 1:length(semillas)) {
  tinic<-Sys.time()
  salida<-FlowSOM(input=datos,silent = FALSE, nClus = 60,colsToUse=c(8:19,21:34,36:47), scale = TRUE, seed = semillas[s])
  tfin<-Sys.time()
  tiempo[s]<-tfin-tinic
  etiquetas_ibgm[[s]]<-data.frame(Manual=exprs(datos)[,"OmiqFilter"], Clustering=GetMetaclusters(salida))
  
  n_cells[s]<-length(unique(etiquetas_ibgm[[s]][,"Manual"]))
  n_clusters[s]<-length(unique(etiquetas_ibgm[[s]][,"Clustering"]))
  pr_mat[[s]] <- re_mat[[s]] <- F_mat[[s]] <- matrix(NA, nrow = n_clusters[s], ncol = n_cells[s])
  resultados[[s]]<-data.frame(matrix(rep(0,8),nrow = 2, ncol = 4))
  colnames(resultados[[s]])<-c("Pr","Re","F","NCelulas")
  rownames(resultados[[s]])<-c("Maximo", "Hungaro")
  #Calcular matrices de métricas
  for (i in 1:n_clusters[s]) {
    for (j in 1:n_cells[s]) {
      tp<-sum(etiquetas_ibgm[[s]][,"Clustering"]==i & etiquetas_ibgm[[s]][,"Manual"]==unique(etiquetas_ibgm[[s]][,"Manual"])[j])
      fp<-sum(etiquetas_ibgm[[s]][,"Clustering"]==i & etiquetas_ibgm[[s]][,"Manual"]!=unique(etiquetas_ibgm[[s]][,"Manual"])[j])
      fn<-sum(etiquetas_ibgm[[s]][,"Clustering"]!=i & etiquetas_ibgm[[s]][,"Manual"]==unique(etiquetas_ibgm[[s]][,"Manual"])[j])
      if(tp==0){
        pr_mat[[s]][i,j]<-0
        re_mat[[s]][i,j]<-0
      }else{
        pr_mat[[s]][i,j]<-tp/(tp+fp)
        re_mat[[s]][i,j]<-tp/(tp+fn)
      }
      if(pr_mat[[s]][i,j]==0 & re_mat[[s]][i,j] ==0){
        F_mat[[s]][i,j]=0
      }else{
        F_mat[[s]][i,j]<-(2*pr_mat[[s]][i,j]*re_mat[[s]][i,j])/(pr_mat[[s]][i,j]+re_mat[[s]][i,j])
      } 
    }
  }
  
  #Resuleve con el algoritmo del húngaro
  if(n_clusters[s]<n_cells[s]){
    #Las poblaciones sin cluster se les asignará NA y tendrán valor 0 para pr, re y F
    asignacion[[s]]<-solve_LSAP(F_mat[[s]], maximum = TRUE)
    cat("Numero clusters menor que el numero de grupos celulares\n")
    for (i in 1:n_clusters[s]) {
      j_cell<-asignacion[[s]][i]
      resultados[[s]]["Hungaro","Pr"]<-resultados[[s]]["Hungaro","Pr"]+pr_mat[[s]][i,j_cell]
      resultados[[s]]["Hungaro","Re"]<-resultados[[s]]["Hungaro","Re"]+re_mat[[s]][i,j_cell]
      resultados[[s]]["Hungaro","F"]<-resultados[[s]]["Hungaro","F"]+F_mat[[s]][i,j_cell]
    }
    resultados[[s]]["Hungaro","Pr"]<-resultados[[s]]["Hungaro","Pr"]/n_cells[s]
    resultados[[s]]["Hungaro","Re"]<-resultados[[s]]["Hungaro","Re"]/n_cells[s]
    resultados[[s]]["Hungaro","F"]<-resultados[[s]]["Hungaro","F"]/n_cells[s]
    resultados[[s]]["Hungaro","NCelulas"]<-n_cells[s]
  }else{
    #Trabajo con la matriz traspuesta
    cat("Numero clusters mayor o igual que el numero de grupos celulares\n")
    asignacion[[s]]<-solve_LSAP(t(F_mat[[s]]), maximum = TRUE)
    for (j in 1:n_cells[s]) {
      i_clus<-asignacion[[s]][j]
      resultados[[s]]["Hungaro","Pr"]<-resultados[[s]]["Hungaro","Pr"]+pr_mat[[s]][i_clus,j]
      resultados[[s]]["Hungaro","Re"]<-resultados[[s]]["Hungaro","Re"]+re_mat[[s]][i_clus,j]
      resultados[[s]]["Hungaro","F"]<-resultados[[s]]["Hungaro","F"]+F_mat[[s]][i_clus,j]
    }
    resultados[[s]]["Hungaro","Pr"]<-resultados[[s]]["Hungaro","Pr"]/n_cells[s]
    resultados[[s]]["Hungaro","Re"]<-resultados[[s]]["Hungaro","Re"]/n_cells[s]
    resultados[[s]]["Hungaro","F"]<-resultados[[s]]["Hungaro","F"]/n_cells[s]
    resultados[[s]]["Hungaro","NCelulas"]<-n_cells[s]
  }
  #Resolución según el máximo
  asignacion_maximo[[s]]<-vector("list",length = n_cells[s])
  for (j in 1:n_cells[s]) {
    i_max<-which.max(F_mat[[s]][,j])
    asignacion_maximo[[s]][[j]]<-i_max
    resultados[[s]]["Maximo","Pr"]<-resultados[[s]]["Maximo","Pr"]+sum(etiquetas_ibgm[[s]][,"Manual"]==unique(etiquetas_ibgm[[s]][,"Manual"])[j])*pr_mat[[s]][i_max,j]
    resultados[[s]]["Maximo","Re"]<-resultados[[s]]["Maximo","Re"]+sum(etiquetas_ibgm[[s]][,"Manual"]==unique(etiquetas_ibgm[[s]][,"Manual"])[j])*re_mat[[s]][i_max,j]
    resultados[[s]]["Maximo","F"]<-resultados[[s]]["Maximo","F"]+sum(etiquetas_ibgm[[s]][,"Manual"]==unique(etiquetas_ibgm[[s]][,"Manual"])[j])*F_mat[[s]][i_max,j]
    resultados[[s]]["Maximo","NCelulas"]<-resultados[[s]]["Maximo","NCelulas"]+sum(etiquetas_ibgm[[s]][,"Manual"]==unique(etiquetas_ibgm[[s]][,"Manual"])[j])
  }
  resultados[[s]]["Maximo","Pr"]<-resultados[[s]]["Maximo","Pr"]/resultados[[s]]["Maximo","NCelulas"]
  resultados[[s]]["Maximo","Re"]<-resultados[[s]]["Maximo","Re"]/resultados[[s]]["Maximo","NCelulas"]
  resultados[[s]]["Maximo","F"]<-resultados[[s]]["Maximo","F"]/resultados[[s]]["Maximo","NCelulas"]
  resultados[[s]]
}
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

metricas_boot<-vector(mode="list",length = length(semillas))
for (s in 1:length(semillas)) {
  metricas_boot[[s]]<-boot(etiquetas_ibgm[[s]],statistic = funcion_boot, R=n_boot)$t
}


