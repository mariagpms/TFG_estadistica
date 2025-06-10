#################################################################################
#FlowSOM óptimo (número de clusters prefijado) para el conjunto IBGM, jerárquico con todos los marcadores
################################################################################
library(flowCore)
library(flowMeans)
library(dplyr)
library(clue)
library(boot)
library(readr)

directorio<-"D:/datos/33_files_concat_especifico.fcs"
datos<-read.FCS(directorio, transformation = FALSE, truncate_max_range = FALSE)

library(readr)
equivalencias <- read_csv("D:/equivalencias.csv")

#IC
n_boot<-10
library(boot)
#Funcion para calcular los IC de cheap bootstrap
funcion_boot<-function(datos, indices){
  temp<-datos[indices,]
  new_etiquetas<-unique(temp$Manual)[unique(temp$Manual)!="Otro"]
  n_cells<-length(new_etiquetas)
  n_clusters<-length(unique(temp[,"Clustering"]))
  pr_mat <- re_mat <- F_mat <- matrix(NA, nrow = n_clusters, ncol = n_cells)
  resultados<-rep(0,10)
  names(resultados)<-c("Pr_hungaro","Re_hungaro","F_hungaro","NCelulas_hungaro","Accuracy_hungaro","Pr_maximo","Re_maximo","F_maximo","NCelulas_maximo","Accuracy_maximo")
  #calcular matrices de métricas
  for (i in 1:n_clusters) {
    for (j in 1:n_cells) {
      tp<-sum(temp[,"Clustering"]==i & temp[,"Manual"]==new_etiquetas[j])
      fp<-sum(temp[,"Clustering"]==i & temp[,"Manual"]!=new_etiquetas[j])
      fn<-sum(temp[,"Clustering"]!=i & temp[,"Manual"]==new_etiquetas[j])
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
  
  library(clue)
  #Resuleve con el algoritmo del húngaro
  if(n_clusters<n_cells){
    #Las poblaciones sin cluster se les asignará NA y tendrán valor 0 para pr, re y F
    asignacion<-solve_LSAP(F_mat, maximum = TRUE)
    cat("Numero clusters menor que el numero de grupos celulares\n")
    numero<-0
    asignacion_inversa<-rep(0,n_cells)
    for (i in 1:n_clusters) {
      j_cell<-asignacion[i]
      asignacion_inversa[j_cell]<-i
      numero<-numero+sum(temp[,"Clustering"]==i & temp[,"Manual"]==new_etiquetas[j_cell])
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
    numero<-0
    for (j in 1:n_cells) {
      i_clus<-asignacion[j]
      numero<-numero+sum(temp[,"Clustering"]==i_clus & temp[,"Manual"]==new_etiquetas[j])
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
  #Resolución según el criterio del máximo
  asignacion_maximo<-rep(0, n_cells)
  numero<-0
  for (j in 1:n_cells) {
    i_max<-which.max(F_mat[,j])
    numero<-numero+sum(temp[,"Clustering"]==i_max & temp[,"Manual"]==new_etiquetas[j])
    asignacion_maximo[j]<-i_max
    resultados["Pr_maximo"]<-resultados["Pr_maximo"]+sum(temp[,"Manual"]==new_etiquetas[j])*pr_mat[i_max,j]
    resultados["Re_maximo"]<-resultados["Re_maximo"]+sum(temp[,"Manual"]==new_etiquetas[j])*re_mat[i_max,j]
    resultados["F_maximo"]<-resultados["F_maximo"]+sum(temp[,"Manual"]==new_etiquetas[j])*F_mat[i_max,j]
    resultados["NCelulas_maximo"]<-resultados["NCelulas_maximo"]+sum(temp[,"Manual"]==new_etiquetas[j])
  }
  resultados["Accuracy_maximo"]<-mean(temp$Clustering==asignacion_maximo[match(temp$Manual, unique(temp$Manual))])
  resultados["Pr_maximo"]<-resultados["Pr_maximo"]/resultados["NCelulas_maximo"]
  resultados["Re_maximo"]<-resultados["Re_maximo"]/resultados["NCelulas_maximo"]
  resultados["F_maximo"]<-resultados["F_maximo"]/resultados["NCelulas_maximo"]
  return(resultados)
}

#Evaluar resultados
evaluar_clustering<-function(etiquetas_ibgm){
  new_etiquetas<-unique(etiquetas_ibgm$Manual)[unique(etiquetas_ibgm$Manual)!="Otro"]
  n_cells<-length(new_etiquetas)
  n_clusters<-max(unique(as.numeric(etiquetas_ibgm[,"Clustering"])))
  pr_mat <- re_mat <- F_mat<-matrix(NA, nrow = n_clusters, ncol = n_cells)
  colnames(pr_mat)<-colnames(re_mat)<-colnames(F_mat)<-paste(new_etiquetas)
  rownames(pr_mat)<-rownames(re_mat)<-rownames(F_mat)<-paste("Clus",1:n_clusters)
  resultados<-data.frame(matrix(rep(0,10),nrow = 2, ncol = 5))
  colnames(resultados)<-c("Pr","Re","F","NCelulas","Accuracy")
  rownames(resultados)<-c("Maximo", "Hungaro")
  for (i in 1:n_clusters) {
    for (j in 1:n_cells) {
      tp<-sum(etiquetas_ibgm[,"Clustering"]==i & etiquetas_ibgm[,"Manual"]==new_etiquetas[j])
      fp<-sum(etiquetas_ibgm[,"Clustering"]==i & etiquetas_ibgm[,"Manual"]!=new_etiquetas[j])
      fn<-sum(etiquetas_ibgm[,"Clustering"]!=i & etiquetas_ibgm[,"Manual"]==new_etiquetas[j])
      tn<-sum(etiquetas_ibgm[,"Clustering"]!=i & etiquetas_ibgm[,"Manual"]!=new_etiquetas[j])
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
  
  library(clue)
  #Ahora se resuleve el problema de asignacion del húngaro
  if(n_clusters<n_cells){
    #Las poblaciones sin cluster se les asignará NA y tendrán valor 0 para pr, re y F
    asignacion<-solve_LSAP(F_mat, maximum = TRUE)
    cat("Numero clusters menor que el numero de grupos celulares\n")
    asignacion_inversa<-rep(0,n_cells)
    numero<-0
    for (i in 1:n_clusters) {
      j_cell<-asignacion[i]
      asignacion_inversa[j_cell]<-i
      numero<-numero+sum(etiquetas_ibgm[,"Clustering"]==i & etiquetas_ibgm[,"Manual"]==new_etiquetas[j_cell])
      resultados["Hungaro","Pr"]<-resultados["Hungaro","Pr"]+pr_mat[i,j_cell]
      resultados["Hungaro","Re"]<-resultados["Hungaro","Re"]+re_mat[i,j_cell]
      resultados["Hungaro","F"]<-resultados["Hungaro","F"]+F_mat[i,j_cell]
    }
    resultados["Hungaro","Accuracy"]<-numero/length(etiquetas_ibgm$Manual)
    resultados["Hungaro","Pr"]<-resultados["Hungaro","Pr"]/n_cells
    resultados["Hungaro","Re"]<-resultados["Hungaro","Re"]/n_cells
    resultados["Hungaro","F"]<-resultados["Hungaro","F"]/n_cells
    resultados["Hungaro","NCelulas"]<-n_cells
  }else{
    #Trabajo con la matriz traspuesta
    cat("Numero clusters mayor o igual que el numero de grupos celulares\n")
    asignacion<-solve_LSAP(t(F_mat), maximum = TRUE)
    numero<-0
    for (j in 1:n_cells) {
      i_clus<-asignacion[j]
      numero<-numero+sum(etiquetas_ibgm[,"Clustering"]==i_clus & etiquetas_ibgm[,"Manual"]==new_etiquetas[j])
      resultados["Hungaro","Pr"]<-resultados["Hungaro","Pr"]+pr_mat[i_clus,j]
      resultados["Hungaro","Re"]<-resultados["Hungaro","Re"]+re_mat[i_clus,j]
      resultados["Hungaro","F"]<-resultados["Hungaro","F"]+F_mat[i_clus,j]
    }
    resultados["Hungaro","Accuracy"]<-numero/length(etiquetas_ibgm$Manual)
    resultados["Hungaro","Pr"]<-resultados["Hungaro","Pr"]/n_cells
    resultados["Hungaro","Re"]<-resultados["Hungaro","Re"]/n_cells
    resultados["Hungaro","F"]<-resultados["Hungaro","F"]/n_cells
    resultados["Hungaro","NCelulas"]<-n_cells
  }
  asignacion_maximo<-rep(0,n_cells)
  numero<-0
  for (j in 1:n_cells) {
    i_max<-which.max(F_mat[,j])
    asignacion_maximo[j]<-i_max
    numero<-numero+sum(etiquetas_ibgm[,"Clustering"]==i_max & etiquetas_ibgm[,"Manual"]==new_etiquetas[j])
    resultados["Maximo","Pr"]<-resultados["Maximo","Pr"]+sum(etiquetas_ibgm[,"Manual"]==new_etiquetas[j])*pr_mat[i_max,j]
    resultados["Maximo","Re"]<-resultados["Maximo","Re"]+sum(etiquetas_ibgm[,"Manual"]==new_etiquetas[j])*re_mat[i_max,j]
    resultados["Maximo","F"]<-resultados["Maximo","F"]+sum(etiquetas_ibgm[,"Manual"]==unique(etiquetas_ibgm[,"Manual"])[j])*F_mat[i_max,j]
    resultados["Maximo","NCelulas"]<-resultados["Maximo","NCelulas"]+sum(etiquetas_ibgm[,"Manual"]==new_etiquetas[j])
  }
  resultados["Maximo","Accuracy"]<-numero/length(etiquetas_ibgm$Manual)
  resultados["Maximo","Pr"]<-resultados["Maximo","Pr"]/resultados["Maximo","NCelulas"]
  resultados["Maximo","Re"]<-resultados["Maximo","Re"]/resultados["Maximo","NCelulas"]
  resultados["Maximo","F"]<-resultados["Maximo","F"]/resultados["Maximo","NCelulas"]
  resultados  
  
  return(list(F_mat = F_mat, pr_mat = pr_mat, re_mat = re_mat, resultados = resultados, asignacion_hungaro=asignacion, asignacion_maximo=asignacion_maximo))
}
##########################################################################################
equivalencias$Clasif3<-"Otro"
equivalencias$Clasif3[c(35,47,42,32,23,48,28,17,46,33,4,39,51,24,8,41,14,30,16,15,26,3,25,5,13,21,18,50,36,40,19,9,52)]<-"T cells"
equivalencias$Clasif3[c(7,38,22,20,11,10,31,49,6,27)]<-"B cells"
equivalencias$Clasif3[c(44,29,1,37,2,12,43,53,45,34)]<-"No B & No T"

table(equivalencias$Clasif3[exprs(datos)[,"OmiqFilter"]])*100/length(exprs(datos)[,"OmiqFilter"])
#Clasificacion en las 3 grandes categorias
salida<-FlowSOM(input=datos,silent = FALSE, scale = TRUE,nClus = 3,colsToUse=c(8:19,21:34,36:47))
etiquetas_ibgm<-data.frame(Manual=equivalencias$Clasif3[exprs(datos)[,"OmiqFilter"]], Clustering=GetMetaclusters(salida))


result3<-evaluar_clustering(etiquetas_ibgm = etiquetas_ibgm)
result3$resultados
result3$asignacion_maximo
result3$asignacion_hungaro
unique(etiquetas_ibgm$Manual)
result3_boot<-boot(etiquetas_ibgm, statistic = funcion_boot, R = n_boot)

#B cells - cluster 3
#T cells - cluster 1
#No B & No T - cluster 2
##########################################################################################################################
#Para las T cells
datos_tcells<-Subset(datos,etiquetas_ibgm$Clustering==1)
salida_tcells<-FlowSOM(datos_tcells,silent = FALSE, scale=TRUE, nClus = 4, colsToUse=c(8:19,21:34,36:47))
equivalencias$ClasifTcells<-"Otro"
equivalencias$ClasifTcells[c(35,47,42,32,23,48)]<-"Tgd"
equivalencias$ClasifTcells[c(28,17,46,33,4)]<-"NKT"
equivalencias$ClasifTcells[c(39,51,24,8,41,14,30,16,15,26,3,25,5,13,21,18,50,36,40,19,9,52)]<-"Real T cells"

table(equivalencias$ClasifTcells[exprs(datos_tcells)[,"OmiqFilter"]])*100/length(exprs(datos_tcells)[,"OmiqFilter"])
etiquetas_tcells<-data.frame(Manual=equivalencias$ClasifTcells[exprs(datos_tcells)[,"OmiqFilter"]], Clustering=GetMetaclusters(salida_tcells))
result_tcells<-evaluar_clustering(etiquetas_tcells)
result_tcells$resultados
unique(etiquetas_tcells$Manual)
result_tcells$asignacion_maximo
result_tcells_boot<-boot(etiquetas_tcells, statistic = funcion_boot, R = n_boot)


#Real T cells - cluster 1
#NKT - cluster 1
#Tgd - cluster 4

##############################################################################################################################
#Para las Real T cells
datos_realtcells<-Subset(datos_tcells,etiquetas_tcells$Clustering==1)
salida_realtcells<-FlowSOM(datos_realtcells,silent = FALSE, scale=TRUE, nClus = 5, colsToUse=c(8:19,21:34,36:47))
equivalencias$ClasifRealTcells<-"Otro"
equivalencias$ClasifRealTcells[c(8)]<-"DP"
equivalencias$ClasifRealTcells[c(51,41,14,30,16,15)]<-"CD8"
equivalencias$ClasifRealTcells[c(24,26,3,25,5,13,21,18,50,36,40,19,9,52)]<-"CD4"
equivalencias$ClasifRealTcells[c(39)]<-"T cells"
table(equivalencias$ClasifRealTcells[exprs(datos_realtcells)[,"OmiqFilter"]])*100/length(exprs(datos_realtcells)[,"OmiqFilter"])

etiquetas_realtcells<-data.frame(Manual=equivalencias$ClasifRealTcells[exprs(datos_realtcells)[,"OmiqFilter"]], Clustering=GetMetaclusters(salida_realtcells))
result_realtcells<-evaluar_clustering(etiquetas_realtcells)
result_realtcells$resultados
result_realtcells$asignacion_maximo
unique(etiquetas_realtcells$Manual)
result_realtcells_boot<-boot(etiquetas_realtcells, statistic = funcion_boot, R = n_boot)
#CD8 - cluster 1
#CD4 - cluster 1
#DP - cluster 1
#T cells -cluster 3
##########################################################################################
#Para las CD8
datos_cd8<-Subset(datos_realtcells, etiquetas_realtcells$Clustering==1)
salida_cd8<-FlowSOM(datos_cd8,silent = FALSE, scale=TRUE, nClus = 7, colsToUse=c(8:19,21:34,36:47))
equivalencias$ClasifCd8<-"Otro"
equivalencias$ClasifCd8[51]<-"CD8"
equivalencias$ClasifCd8[41]<-"Trm"
equivalencias$ClasifCd8[14]<-"Memory"
equivalencias$ClasifCd8[30]<-"Naive"
equivalencias$ClasifCd8[16]<-"nrTEMRA"
equivalencias$ClasifCd8[15]<-"rTEMRA"
table(equivalencias$ClasifCd8[exprs(datos_cd8)[,"OmiqFilter"]])*100/length(exprs(datos_cd8)[,"OmiqFilter"])

etiquetas_cd8<-data.frame(Manual=equivalencias$ClasifCd8[exprs(datos_cd8)[,"OmiqFilter"]], Clustering=GetMetaclusters(salida_cd8))
result_cd8<-evaluar_clustering(etiquetas_cd8)
result_cd8$resultados
result_cd8$F_mat
result_cd8$asignacion_maximo
unique(etiquetas_cd8$Manual)
result_cd8_boot<-boot(etiquetas_cd8, statistic = funcion_boot, R = n_boot)

#Trm - 2
#Memory - 2
#Naive - 3
#rTEMRA - 2
#nrTEMRA - 2
#CD8 - 2
##############################################################################################
#Para las CD4
datos_cd4<-Subset(datos_realtcells, etiquetas_realtcells$Clustering==1)
salida_cd4<-FlowSOM(datos_cd4,silent = FALSE, scale=TRUE, nClus = 15, colsToUse=c(8:19,21:34,36:47))
equivalencias$ClasifCd4<-"Otro"
equivalencias$ClasifCd4[24]<-"CD4"
equivalencias$ClasifCd4[26]<-"Trm"
equivalencias$ClasifCd4[21]<-"Trm/Th2"
equivalencias$ClasifCd4[18]<-"Trm/Th1"
equivalencias$ClasifCd4[50]<-"Trm/Th17"
equivalencias$ClasifCd4[36]<-"Trm/Th1/Th17-like"
equivalencias$ClasifCd4[3]<-"Memory"
equivalencias$ClasifCd4[52]<-"Memory/Th2"
equivalencias$ClasifCd4[40]<-"Memory/Th1"
equivalencias$ClasifCd4[19]<-"Memory/Th17"
equivalencias$ClasifCd4[9]<-"Memory/Th1/Th17-like"
equivalencias$ClasifCd4[25]<-"Naive"
equivalencias$ClasifCd4[5]<-"nrTEMRA"
equivalencias$ClasifCd4[13]<-"rTEMRA"
table(equivalencias$ClasifCd4[exprs(datos_cd4)[,"OmiqFilter"]])*100/length(exprs(datos_cd4)[,"OmiqFilter"])

etiquetas_cd4<-data.frame(Manual=equivalencias$ClasifCd4[exprs(datos_cd4)[,"OmiqFilter"]], Clustering=GetMetaclusters(salida_cd4))
result_cd4<-evaluar_clustering(etiquetas_cd4)
result_cd4$resultados
result_cd4$F_mat
result_cd4$asignacion_maximo
unique(etiquetas_cd4$Manual)
result_cd4_boot<-boot(etiquetas_cd4, statistic = funcion_boot, R = n_boot)

#Trm/Th17 - 7
#Memory/Th1 - 7
#Trm/Th1/Th17-like - 7
#Trm/Th1 - 7
#Memory/Th1/Th17-like - 7

#Naive - 2
#Trm/Th2 - 12
#rTEMRA - 11
#Memory/Th2 - 13
#Memory/Th17 - 7
#nrTEMRA - 9
#CD4 - 11
####################################################################################
#para las NKT
datos_NKT<-Subset(datos_tcells,etiquetas_tcells$Clustering==1)
salida_NKT<-FlowSOM(datos_NKT,silent = FALSE, scale=TRUE, nClus = 6, colsToUse=c(8:19,21:34,36:47))
equivalencias$ClasifNKT<-"Otro"
equivalencias$ClasifNKT[28]<-"NKT"
equivalencias$ClasifNKT[17]<-"CD2-"
equivalencias$ClasifNKT[46]<-"CD2+CD8high"
equivalencias$ClasifNKT[33]<-"CD2+CD8dim"
equivalencias$ClasifNKT[4]<-"CD2+CD8-"

table(equivalencias$ClasifNKT[exprs(datos_NKT)[,"OmiqFilter"]])*100/length(exprs(datos_NKT)[,"OmiqFilter"])

etiquetas_NKT<-data.frame(Manual=equivalencias$ClasifNKT[exprs(datos_NKT)[,"OmiqFilter"]], Clustering=GetMetaclusters(salida_NKT))
result_NKT<-evaluar_clustering(etiquetas_NKT)
result_NKT$resultados
result_NKT$asignacion_maximo
unique(etiquetas_NKT$Manual)
result_NKT_boot<-boot(etiquetas_NKT, statistic = funcion_boot, R = n_boot)

#CD2- - 1
#NKT - 6
#CD2+CD8high -1
#CD2+CD8dim -1
#CD2+CD8- - 1
####################################################################################
#para las Tgd
datos_Tgd<-Subset(datos_tcells,etiquetas_tcells$Clustering==2)
salida_Tgd<-FlowSOM(datos_Tgd,silent = FALSE, scale=TRUE, nClus = 7, colsToUse=c(8:19,21:34,36:47))
equivalencias$ClasifTgd<-"Otro"
equivalencias$ClasifTgd[35]<-"Tgd"
equivalencias$ClasifTgd[47]<-"Trm"
equivalencias$ClasifTgd[42]<-"Memory"
equivalencias$ClasifTgd[32]<-"Naive"
equivalencias$ClasifTgd[23]<-"nrTEMRA"
equivalencias$ClasifTgd[48]<-"rTEMRA"


table(equivalencias$ClasifTgd[exprs(datos_Tgd)[,"OmiqFilter"]])*100/length(exprs(datos_Tgd)[,"OmiqFilter"])

etiquetas_Tgd<-data.frame(Manual=equivalencias$ClasifTgd[exprs(datos_Tgd)[,"OmiqFilter"]], Clustering=GetMetaclusters(salida_Tgd))
result_Tgd<-evaluar_clustering(etiquetas_Tgd)
result_Tgd$resultados
result_Tgd$asignacion_maximo
unique(etiquetas_Tgd$Manual)
result_Tgd_boot<-boot(etiquetas_Tgd, statistic = funcion_boot, R = n_boot)

#"Memory"  - 1
#"rTEMRA"  -3
#"Naive"  
#"Tgd"     -6
#"nrTEMRA" -1
#"Trm"- 4
######################################################################################
#Para las B cells
datos_bcells<-Subset(datos,etiquetas_ibgm$Clustering==1)
salida_bcells<-FlowSOM(datos_bcells,silent = FALSE, scale=TRUE, nClus = 4, colsToUse=c(8:19,21:34,36:47))
equivalencias$Clasifbcells<-"Otro"
equivalencias$Clasifbcells[7]<-"B cells"
equivalencias$Clasifbcells[38]<-"IgD+"
equivalencias$Clasifbcells[c(22,20,11,10,31,49,6,27)]<-"IgD-"

table(equivalencias$Clasifbcells[exprs(datos_bcells)[,"OmiqFilter"]])*100/length(exprs(datos_bcells)[,"OmiqFilter"])
etiquetas_bcells<-data.frame(Manual=equivalencias$Clasifbcells[exprs(datos_bcells)[,"OmiqFilter"]], Clustering=GetMetaclusters(salida_bcells))
result_bcells<-evaluar_clustering(etiquetas_bcells)
result_bcells$resultados
result_bcells$asignacion_maximo
unique(etiquetas_bcells$Manual)
result_bcells_boot<-boot(etiquetas_bcells, statistic = funcion_boot, R = n_boot)

#IgD- 1
#IgD+ 1
#B cells 2
######################################################################################
#Para las IgD-
datos_igd<-Subset(datos_bcells,etiquetas_bcells$Clustering==1)
salida_igd<-FlowSOM(datos_igd,silent = FALSE, scale=TRUE, nClus = 4, colsToUse=c(8:19,21:34,36:47))
equivalencias$ClasifIgd<-"Otro"
equivalencias$ClasifIgd[c(22,11,10,31)]<-"Plasmablasts"
equivalencias$ClasifIgd[c(20,49,6,27)]<-"Memory B cells"

table(equivalencias$ClasifIgd[exprs(datos_igd)[,"OmiqFilter"]])*100/length(exprs(datos_igd)[,"OmiqFilter"])
etiquetas_igd<-data.frame(Manual=equivalencias$ClasifIgd[exprs(datos_igd)[,"OmiqFilter"]], Clustering=GetMetaclusters(salida_igd))
result_igd<-evaluar_clustering(etiquetas_igd)
result_igd$resultados
result_igd$asignacion_maximo
unique(etiquetas_igd$Manual)
result_igd_boot<-boot(etiquetas_igd, statistic = funcion_boot, R = n_boot)

#Plasmablasts - 1
#Memory B cells - 3
######################################################################################
#Para las plasmablasts
datos_plasmablasts<-Subset(datos_igd,etiquetas_Igd$Clustering==1)
salida_plasmablasts<-FlowSOM(datos_plasmablasts,silent = FALSE, scale=TRUE, nClus = 5, colsToUse=c(8:19,21:34,36:47))
equivalencias$Clasifplasmablasts<-"Otro"
equivalencias$Clasifplasmablasts[22]<-"Plasmablasts"
equivalencias$Clasifplasmablasts[11]<-"Plasmablasts/IgA"
equivalencias$Clasifplasmablasts[10]<-"Plasmablasts/IgG"
equivalencias$Clasifplasmablasts[31]<-"Plasmablasts/IgM"


table(equivalencias$Clasifplasmablasts[exprs(datos_plasmablasts)[,"OmiqFilter"]])*100/length(exprs(datos_plasmablasts)[,"OmiqFilter"])
etiquetas_plasmablasts<-data.frame(Manual=equivalencias$Clasifplasmablasts[exprs(datos_plasmablasts)[,"OmiqFilter"]], Clustering=GetMetaclusters(salida_plasmablasts))
result_plasmablasts<-evaluar_clustering(etiquetas_plasmablasts)
result_plasmablasts$resultados
result_plasmablasts$asignacion_maximo
unique(etiquetas_plasmablasts$Manual)
result_plasmablasts_boot<-boot(etiquetas_plasmablasts, statistic = funcion_boot, R = n_boot)

#Plasmablasts - 1
#Plasmablasts/IgA  1
#Plasmablasts/IgG - 1
#Plasmablasts/IgM - 1
######################################################################################
#Para las memoryBcells
datos_memoryBcells<-Subset(datos_igd,etiquetas_Igd$Clustering==2)
salida_memoryBcells<-FlowSOM(datos_memoryBcells,silent = FALSE, scale=TRUE, nClus = 4, colsToUse=c(8:19,21:34,36:47))
equivalencias$ClasifmemoryBcells<-"Otro"
equivalencias$ClasifmemoryBcells[20]<-"memoryBcells"
equivalencias$ClasifmemoryBcells[49]<-"memoryBcells/IgA"
equivalencias$ClasifmemoryBcells[6]<-"memoryBcells/IgG"
equivalencias$ClasifmemoryBcells[27]<-"memoryBcells/IgM"


table(equivalencias$ClasifmemoryBcells[exprs(datos_memoryBcells)[,"OmiqFilter"]])*100/length(exprs(datos_memoryBcells)[,"OmiqFilter"])
etiquetas_memoryBcells<-data.frame(Manual=equivalencias$ClasifmemoryBcells[exprs(datos_memoryBcells)[,"OmiqFilter"]], Clustering=GetMetaclusters(salida_memoryBcells))
result_memoryBcells<-evaluar_clustering(etiquetas_memoryBcells)
result_memoryBcells$resultados
result_memoryBcells$asignacion_maximo
unique(etiquetas_memoryBcells$Manual)
result_memoryBcells_boot<-boot(etiquetas_memoryBcells, statistic = funcion_boot, R = n_boot)

#memoryBcells - 1
#memoryBcells/IgA - 1
#memoryBcells/IgG - 1
#memoryBcells/IgM - 1
##########################################################################################################################
#Para las no T no B cells
datos_nono<-Subset(datos,etiquetas_ibgm$Clustering==2)
salida_nono<-FlowSOM(datos_nono,silent = FALSE, scale=TRUE, nClus = 10 , colsToUse=c(8:19,21:34,36:47))
equivalencias$Clasifnono<-"Otro"
equivalencias$Clasifnono[12]<-"Basophils"
equivalencias$Clasifnono[2]<-"ILC"
equivalencias$Clasifnono[c(43,53,45,34)]<-"APC"
equivalencias$Clasifnono[c(44,29,1,37)]<-"DR-"

table(equivalencias$Clasifnono[exprs(datos_nono)[,"OmiqFilter"]])*100/length(exprs(datos_nono)[,"OmiqFilter"])
etiquetas_nono<-data.frame(Manual=equivalencias$Clasifnono[exprs(datos_nono)[,"OmiqFilter"]], Clustering=GetMetaclusters(salida_nono))
result_nono<-evaluar_clustering(etiquetas_nono)
result_nono$resultados
result_nono$asignacion_maximo
unique(etiquetas_nono$Manual)
result_nono_boot<-boot(etiquetas_nono, statistic = funcion_boot, R = n_boot)
#DR-
#APC-2
#ILC-
#basophils - 4
#############################################################################################################
#Para las DR-
datos_dr<-Subset(datos_nono,etiquetas_nono$Clustering==4)
salida_dr<-FlowSOM(datos_dr,silent = FALSE, scale=TRUE, nClus = 5 , colsToUse=c(8:19,21:34,36:47))
equivalencias$Clasifdr<-"Otro"
equivalencias$Clasifdr[44]<-"DR-"
equivalencias$Clasifdr[1]<-"Early NK"
equivalencias$Clasifdr[29]<-"Mature NK"
equivalencias$Clasifdr[37]<-"Terminal NK"

table(equivalencias$Clasifdr[exprs(datos_dr)[,"OmiqFilter"]])*100/length(exprs(datos_dr)[,"OmiqFilter"])
etiquetas_dr<-data.frame(Manual=equivalencias$Clasifdr[exprs(datos_dr)[,"OmiqFilter"]], Clustering=GetMetaclusters(salida_dr))
result_dr<-evaluar_clustering(etiquetas_dr)
result_dr$resultados
result_dr$asignacion_maximo
unique(etiquetas_dr$Manual)
result_dr_boot<-boot(etiquetas_dr, statistic = funcion_boot, R = n_boot)

#Teerminal NK - 2
# Early NK - 1
#Mature NK -3
#DR- -2

#############################################################################################################
#Para las APC
datos_apc<-Subset(datos_nono,etiquetas_nono$Clustering==2)
salida_apc<-FlowSOM(datos_apc,silent = FALSE, scale=TRUE, nClus = 5 , colsToUse=c(8:19,21:34,36:47))
equivalencias$Clasifapc<-"Otro"
equivalencias$Clasifapc[43]<-"APC"
equivalencias$Clasifapc[53]<-"Macrophages"
equivalencias$Clasifapc[45]<-"Monocytes"
equivalencias$Clasifapc[34]<-"cdc"

table(equivalencias$Clasifapc[exprs(datos_apc)[,"OmiqFilter"]])*100/length(exprs(datos_apc)[,"OmiqFilter"])
etiquetas_apc<-data.frame(Manual=equivalencias$Clasifapc[exprs(datos_apc)[,"OmiqFilter"]], Clustering=GetMetaclusters(salida_apc))
result_apc<-evaluar_clustering(etiquetas_apc)
result_apc$resultados
result_apc$asignacion_maximo
unique(etiquetas_apc$Manual)
result_apc_boot<-boot(etiquetas_apc, statistic = funcion_boot, R = n_boot)
#monocytes-1
#cdc-
#macrophages-1
#apc-5
##################################################################################################################################
