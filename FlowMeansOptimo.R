################################################################################
#FlowMeans óptimo (tal como se explica en el TFG) para los conjuntos del FlowCAP-I
################################################################################
#Tenemos 2 archivos de datos públicos del artículo de comparison
library(flowCore) #libreria para trabajar con los datos de citometria
library(flowMeans) #libreria que tiene el algoritmo
library(dplyr) #libreria para trabajar con los conjuntos de datos
library(boot) #libreria con la funcion bootstrap
library(clue) #libreria con el algoritmo del hungaro
ficheros<-c("FlowCAP_ND.fcs", "FlowCAP_WNV.fcs")
directorio<-"E:/datos/"

#Se han cargado los archivos para determinar las columnas que corresponden a marcadores ya que son los que se usan en los algoritmos de Clustering
cols_marcadores<-list(
  FlowCAP_ND.fcs = 3:12,
  FlowCAP_WNV.fcs = 3:8
)

tiempos<-list(
  FlowCAP_ND.fcs = NULL,
  FlowCAP_WNV.fcs = NULL)

etiquetas<-list(
  FlowCAP_ND.fcs = NULL,
  FlowCAP_WNV.fcs = NULL)

numClusters<-list(
  FlowCAP_ND.fcs = 7,
  FlowCAP_WNV.fcs = 4
)

#Las etquetas de los archivos de flowcap no son consistentes a lo largo de las muestras por lo que toca desglosar
#los archivos según las muestras, al hacer esto varias muestras dan error con un número de clusters automático.

for (archivo in ficheros) {
  print(archivo)
  datos<-read.FCS(paste(directorio,archivo, sep = ""), transformation = FALSE, truncate_max_range = FALSE)
  muestra<-exprs(datos)[,"sample"]
  datos<-flowCore::split(datos,muestra)
  etiquetas[[archivo]]<-data.frame(matrix(ncol = 3, nrow = 0))
  colnames(etiquetas[[archivo]])<-c("Muestra","Manual","Clustering")
  manual<-lapply(datos, function(d){d<-exprs(d)[,"label"]})
  datos<-lapply(datos, function(d){d<-exprs(d)[,cols_marcadores[[archivo]]]})
  for (j in 1:length(datos)) {
    cat("Muestra",j,"\n")
    if((archivo=="FlowCAP_ND.fcs")&(j %in% c(2,9,25))) next
    tiempos[[archivo]][j] <- system.time({
      salida<- flowMeans(datos[[j]], NumC = numClusters[[archivo]] ,Standardize = FALSE)
    })
    cat(tiempos[[archivo]][j],"\n")
    temp<-data.frame(Muestra=rep(j,length(salida@Label)), Clustering=salida@Label, Manual=manual[[j]])
    etiquetas[[archivo]]<-bind_rows(etiquetas[[archivo]],temp)
  }
}

metricas<-list(
  FlowCAP_ND.fcs = NULL,
  FlowCAP_WNV.fcs = NULL)

#################################################################################
#EVALUACION RESULTADOS

evaluando<-function(archivo, dataset){
  n_muestras<-max(dataset[,"Muestra"])
  resultados_hungaro<-data.frame(matrix(rep(0,n_muestras*5),nrow = n_muestras, ncol = 5))
  colnames(resultados_hungaro)<-c("Pr","Re","F","NCelulas","Accuracy")
  resultados_maximo<-data.frame(matrix(rep(0,n_muestras*5),nrow = n_muestras, ncol = 5))
  colnames(resultados_maximo)<-c("Pr","Re","F","NCelulas","Accuracy")
  for(m in 1:n_muestras){
    if((archivo=="FlowCAP_ND.fcs")&(m %in% c(2,9,25))) next
    temp<-dataset[dataset[,"Muestra"]==m,]
    temp<-temp[!is.na(temp[,"Manual"]),]
    temp<-temp[!is.na(temp[,"Clustering"]),]
    n_cells<-length(unique(temp[!is.na(temp[,"Manual"]),"Manual"]))
    n_clusters<-length(unique(temp[!is.na(temp[,"Clustering"]),"Clustering"]))
    pr_mat <- re_mat <- F_mat <- matrix(NA, nrow = n_clusters, ncol = n_cells)
    #Calculo matrices de F_score, precision y recall
    for (i in 1:n_clusters) {
      for (j in 1:n_cells) {
        tp<-sum(temp[,"Clustering"]==i & temp[,"Manual"]==j)
        fp<-sum(temp[,"Clustering"]==i & temp[,"Manual"]!=j)
        fn<-sum(temp[,"Clustering"]!=i & temp[,"Manual"]==j)
        if(tp==0){
          pr_mat[i,j]<-0
          re_mat[i,j]<-0
        }else{
          pr_mat[i,j]<-tp/(tp+fp)
          re_mat[i,j]<-tp/(tp+fn)
        }
        
        if(pr_mat[i,j]+re_mat[i,j] ==0){
          F_mat[i,j]=0
        }else{
          F_mat[i,j]<-(2*pr_mat[i,j]*re_mat[i,j])/(pr_mat[i,j]+re_mat[i,j])
        }
      }
    }
    #Resolver con el algoritmo del húngaro
    if(n_clusters<n_cells){
      #Las poblaciones sin cluster se les asignará NA y tendrán valor 0 para pr, re y F
      asignacion<-solve_LSAP(F_mat, maximum = TRUE)
      asignacion_inversa<-rep(0,n_cells)
      for (i in 1:n_clusters) {
        j_cell<-asignacion[i]
        asignacion_inversa[j_cell]<-i
        resultados_hungaro[m,"Pr"]<-resultados_hungaro[m,"Pr"]+pr_mat[i,j_cell]
        resultados_hungaro[m,"Re"]<-resultados_hungaro[m,"Re"]+re_mat[i,j_cell]
        resultados_hungaro[m,"F"]<-resultados_hungaro[m,"F"]+F_mat[i,j_cell]
      }
      resultados_hungaro[m,"Accuracy"]<-mean(temp$Clustering==asignacion_inversa[temp$Manual])
      resultados_hungaro[m,"Pr"]<-resultados_hungaro[m,"Pr"]/n_cells
      resultados_hungaro[m,"Re"]<-resultados_hungaro[m,"Re"]/n_cells
      resultados_hungaro[m,"F"]<-resultados_hungaro[m,"F"]/n_cells
      resultados_hungaro[m,"NCelulas"]<-n_cells
    }else{
      #Trabajo con la matriz traspuesta
      asignacion<-solve_LSAP(t(F_mat), maximum = TRUE)
      for (j in 1:n_cells) {
        i_clus<-asignacion[j]
        resultados_hungaro[m,"Pr"]<-resultados_hungaro[m,"Pr"]+pr_mat[i_clus,j]
        resultados_hungaro[m,"Re"]<-resultados_hungaro[m,"Re"]+re_mat[i_clus,j]
        resultados_hungaro[m,"F"]<-resultados_hungaro[m,"F"]+F_mat[i_clus,j]
      }
      resultados_hungaro[m,"Accuracy"]<-mean(temp$Clustering==asignacion[temp$Manual])
      resultados_hungaro[m,"Pr"]<-resultados_hungaro[m,"Pr"]/n_cells
      resultados_hungaro[m,"Re"]<-resultados_hungaro[m,"Re"]/n_cells
      resultados_hungaro[m,"F"]<-resultados_hungaro[m,"F"]/n_cells
      resultados_hungaro[m,"NCelulas"]<-n_cells
    }
    #Resolución según el máximo
    asignacion_maximo<-rep(0,n_cells)
    for (j in 1:n_cells) {
      i_max<-which.max(F_mat[,j])
      asignacion_maximo[j]<-i_max
      resultados_maximo[m,"Pr"]<-resultados_maximo[m,"Pr"]+sum(temp[,"Manual"]==j)*pr_mat[i_max,j]
      resultados_maximo[m,"Re"]<-resultados_maximo[m,"Re"]+sum(temp[,"Manual"]==j)*re_mat[i_max,j]
      resultados_maximo[m,"F"]<-resultados_maximo[m,"F"]+sum(temp[,"Manual"]==j)*F_mat[i_max,j]
      resultados_maximo[m,"NCelulas"]<-resultados_maximo[m,"NCelulas"]+sum(temp[,"Manual"]==j)
    }
    resultados_maximo[m,"Accuracy"]<-mean(temp$Clustering==asignacion_maximo[temp$Manual],na.rm=TRUE)
    resultados_maximo[m,"Pr"]<-resultados_maximo[m,"Pr"]/resultados_maximo[m,"NCelulas"]
    resultados_maximo[m,"Re"]<-resultados_maximo[m,"Re"]/resultados_maximo[m,"NCelulas"]
    resultados_maximo[m,"F"]<-resultados_maximo[m,"F"]/resultados_maximo[m,"NCelulas"]
  }
  return(list("Hungaro"=resultados_hungaro,"Maximo"=resultados_maximo))
}

for (archivo in ficheros) {
  print(archivo)
  temp<-etiquetas[[archivo]]
  metricas[[archivo]]<-evaluando(archivo = archivo, dataset = temp)
  metricas[[archivo]]$ResultadosHungaro<-apply(metricas[[archivo]]$Hungaro, MARGIN = 2, sum)/sum(metricas[[archivo]]$Hungaro$NCelulas!=0)
  metricas[[archivo]]$ResultadosMaximo<-apply(metricas[[archivo]]$Maximo*metricas[[archivo]]$Maximo$NCelulas, MARGIN =2, sum)/sum(metricas[[archivo]]$Maximo$NCelulas)
  print(sum(tiempos[[archivo]][!is.na(tiempos[[archivo]])]))
}
######################################################################################
#IC con cheap bootstrap
metricas_boot<-list(
  FlowCAP_ND.fcs = NULL,
  FlowCAP_WNV.fcs = NULL)
n_boot<-10

funcion_metricas<-function(temp, indices){
  temp<-temp[indices,]
  n_cells<-length(unique(temp[!is.na(temp[,"Manual"]),"Manual"]))
  n_clusters<-length(unique(temp[!is.na(temp[,"Clustering"]),"Clustering"]))
  pr_mat <- re_mat <- F_mat <- matrix(NA, nrow = n_clusters, ncol = n_cells)
  resultados<-rep(0,5*2)
  names(resultados)<-c("Pr_hungaro","Re_hungaro","F_hungaro","NCelulas_hungaro","Accuracy_hungaro","Pr_maximo","Re_maximo","F_maximo","NCelulas_maximo","Accuracy_maximo")
  #Calcular matrices F_score, precision y recall
  for (i in 1:n_clusters) {
    for (j in 1:n_cells) {
      tp<-sum(temp[,"Clustering"]==i & temp[,"Manual"]==j)
      fp<-sum(temp[,"Clustering"]==i & temp[,"Manual"]!=j)
      fn<-sum(temp[,"Clustering"]!=i & temp[,"Manual"]==j)
      if(tp==0){
        pr_mat[i,j]<-0
        re_mat[i,j]<-0
      }else{
        pr_mat[i,j]<-tp/(tp+fp)
        re_mat[i,j]<-tp/(tp+fn)
      }
      
      if(pr_mat[i,j]+re_mat[i,j] ==0){
        F_mat[i,j]=0
      }else{
        F_mat[i,j]<-(2*pr_mat[i,j]*re_mat[i,j])/(pr_mat[i,j]+re_mat[i,j])
      }
    }
  }
  #Resolver con el algoritmo del húngaro
  if(n_clusters<n_cells){
    #Las poblaciones sin cluster se les asignará NA y tendrán valor 0 para pr, re y F
    asignacion<-solve_LSAP(F_mat, maximum = TRUE)
    asignacion_inversa<-rep(0,n_cells)
    for (i in 1:n_clusters) {
      j_cell<-asignacion[i]
      asignacion_inversa[j_cell]<-i
      resultados["Pr_hungaro"]<-resultados["Pr_hungaro"]+pr_mat[i,j_cell]
      resultados["Re_hungaro"]<-resultados["Re_hungaro"]+re_mat[i,j_cell]
      resultados["F_hungaro"]<-resultados["F_hungaro"]+F_mat[i,j_cell]
    }
    resultados["Accuracy_hungaro"]<-mean(temp$Clustering==asignacion_inversa[temp$Manual])
    resultados["Pr_hungaro"]<-resultados["Pr_hungaro"]/n_cells
    resultados["Re_hungaro"]<-resultados["Re_hungaro"]/n_cells
    resultados["F_hungaro"]<-resultados["F_hungaro"]/n_cells
    resultados["NCelulas_hungaro"]<-n_cells
  }else{
    #Trabajo con la matriz traspuesta
    asignacion<-solve_LSAP(t(F_mat), maximum = TRUE)
    for (j in 1:n_cells) {
      i_clus<-asignacion[j]
      resultados["Pr_hungaro"]<-resultados["Pr_hungaro"]+pr_mat[i_clus,j]
      resultados["Re_hungaro"]<-resultados["Re_hungaro"]+re_mat[i_clus,j]
      resultados["F_hungaro"]<-resultados["F_hungaro"]+F_mat[i_clus,j]
    }
    resultados["Accuracy_hungaro"]<-mean(temp$Clustering==asignacion[temp$Manual])
    resultados["Pr_hungaro"]<-resultados["Pr_hungaro"]/n_cells
    resultados["Re_hungaro"]<-resultados["Re_hungaro"]/n_cells
    resultados["F_hungaro"]<-resultados["F_hungaro"]/n_cells
    resultados["NCelulas_hungaro"]<-n_cells
  }
  #Resolución según el máximo
  asignacion_maximo<-rep(0,n_cells)
  for (j in 1:n_cells) {
    i_max<-which.max(F_mat[,j])
    asignacion_maximo[j]<-i_max
    resultados["Pr_maximo"]<-resultados["Pr_maximo"]+sum(temp[,"Manual"]==j)*pr_mat[i_max,j]
    resultados["Re_maximo"]<-resultados["Re_maximo"]+sum(temp[,"Manual"]==j)*re_mat[i_max,j]
    resultados["F_maximo"]<-resultados["F_maximo"]+sum(temp[,"Manual"]==j)*F_mat[i_max,j]
    resultados["NCelulas_maximo"]<-resultados["NCelulas_maximo"]+sum(temp[,"Manual"]==j)
  }
  resultados["Accuracy_maximo"]<-mean(temp$Clustering==asignacion_maximo[temp$Manual])
  resultados["Pr_maximo"]<-resultados["Pr_maximo"]/resultados["NCelulas_maximo"]
  resultados["Re_maximo"]<-resultados["Re_maximo"]/resultados["NCelulas_maximo"]
  resultados["F_maximo"]<-resultados["F_maximo"]/resultados["NCelulas_maximo"]
  return(resultados)
}

for (archivo in ficheros) {
  cat(archivo, "\n")
  pr_boot<-re_boot<-F_boot<-N_boot<-Accuracy_boot<-data.frame(matrix(rep(0,n_boot*2),nrow = n_boot, ncol = 2))
  colnames(pr_boot)<-colnames(re_boot)<-colnames(F_boot)<-colnames(N_boot)<-colnames(Accuracy_boot)<-c("Hungaro","Maximo")
  dataset<-etiquetas[[archivo]]
  dataset<-dataset[!is.na(dataset[,"Manual"]),]
  dataset<-dataset[!is.na(dataset[,"Clustering"]),]
  muestras<-max(dataset[,"Muestra"])
  result_muestra<-vector("list",length = muestras)
  for (m in 1:muestras) {
    cat(m," ")
    if((archivo=="FlowCAP_ND.fcs")&(m %in% c(2,9,25))) next
    temp<-dataset[dataset[,"Muestra"]==m,]
    boot_temp<-boot(temp,statistic = funcion_metricas, R=n_boot)
    result_muestra[[m]]<-boot_temp$t
    colnames(result_muestra[[m]])<-c("Pr_hungaro","Re_hungaro","F_hungaro","NCelulas_hungaro","Accuracy_hungaro","Pr_maximo","Re_maximo","F_maximo","NCelulas_maximo","Accuracy_maximo")
  }
  #Se calcula el F_score, precision y recall final de cada archivo para cada una de las repeticiones de bootstrap
  for (r in 1:n_boot) {
    for(m in 1:muestras){
      if((archivo=="FlowCAP_ND.fcs")&(m %in% c(2,9,25))) next
      pr_boot[r,"Hungaro"]<-pr_boot[r,"Hungaro"]+result_muestra[[m]][r,"Pr_hungaro"]
      re_boot[r,"Hungaro"]<-re_boot[r,"Hungaro"]+result_muestra[[m]][r,"Re_hungaro"]
      F_boot[r,"Hungaro"]<-F_boot[r,"Hungaro"]+result_muestra[[m]][r,"F_hungaro"]
      N_boot[r,"Hungaro"]<-N_boot[r,"Hungaro"]+result_muestra[[m]][r,"NCelulas_hungaro"]
      Accuracy_boot[r,"Hungaro"]<-Accuracy_boot[r,"Hungaro"]+result_muestra[[m]][r,"Accuracy_hungaro"]
        
      pr_boot[r,"Maximo"]<-pr_boot[r,"Maximo"]+result_muestra[[m]][r,"Pr_maximo"]*result_muestra[[m]][r,"NCelulas_maximo"]
      re_boot[r,"Maximo"]<-re_boot[r,"Maximo"]+result_muestra[[m]][r,"Re_maximo"]*result_muestra[[m]][r,"NCelulas_maximo"]
      F_boot[r,"Maximo"]<-F_boot[r,"Maximo"]+result_muestra[[m]][r,"F_maximo"]*result_muestra[[m]][r,"NCelulas_maximo"]
      N_boot[r,"Maximo"]<-N_boot[r,"Maximo"]+result_muestra[[m]][r,"NCelulas_maximo"]
      Accuracy_boot[r,"Maximo"]<-Accuracy_boot[r,"Maximo"]+result_muestra[[m]][r,"Accuracy_maximo"]*result_muestra[[m]][r,"NCelulas_maximo"]
    }
    pr_boot[r,"Hungaro"]<-pr_boot[r,"Hungaro"]/sum(metricas[[archivo]]$Hungaro$NCelulas!=0)
    re_boot[r,"Hungaro"]<-re_boot[r,"Hungaro"]/sum(metricas[[archivo]]$Hungaro$NCelulas!=0)
    F_boot[r,"Hungaro"]<-F_boot[r,"Hungaro"]/sum(metricas[[archivo]]$Hungaro$NCelulas!=0)
    Accuracy_boot[r,"Hungaro"]<-Accuracy_boot[r,"Hungaro"]/sum(metricas[[archivo]]$Hungaro$NCelulas!=0)
    
    pr_boot[r,"Maximo"]<-pr_boot[r,"Maximo"]/N_boot[r,"Maximo"]
    re_boot[r,"Maximo"]<-re_boot[r,"Maximo"]/N_boot[r,"Maximo"]
    F_boot[r,"Maximo"]<-F_boot[r,"Maximo"]/N_boot[r,"Maximo"]
    Accuracy_boot[r,"Maximo"]<-Accuracy_boot[r,"Maximo"]/N_boot[r,"Maximo"]
  }
  metricas_boot[[archivo]]<-list(Pr=pr_boot, Re=re_boot, F_score=F_boot, N_cells=N_boot, Accuracy=Accuracy_boot)
}

#Almacenar los datos en dataframes
df_nd<-data.frame(Tipo=rep("ND_FlowMeansOpt",2*n_boot),
                  Asignacion=rep(c("Maximo","Hungaro"),each=n_boot),
                  F_score=c(metricas_boot[[1]][[3]]$Maximo,metricas_boot[[1]][[3]]$Hungaro),
                  pr=c(metricas_boot[[1]][[1]]$Maximo,metricas_boot[[1]][[1]]$Hungaro),
                  re=c(metricas_boot[[1]][[2]]$Maximo,metricas_boot[[1]][[2]]$Hungaro),
                  accuracy=c(metricas_boot[[1]][[5]]$Maximo,metricas_boot[[1]][[5]]$Hungaro))

df_wnv<-data.frame(Tipo=rep("WNV_FlowMeansOpt",2*n_boot),
                   Asignacion=rep(c("Maximo","Hungaro"),each=n_boot),
                   F_score=c(metricas_boot[[2]][[3]]$Maximo,metricas_boot[[2]][[3]]$Hungaro),
                   pr=c(metricas_boot[[2]][[1]]$Maximo,metricas_boot[[2]][[1]]$Hungaro),
                   re=c(metricas_boot[[2]][[2]]$Maximo,metricas_boot[[2]][[2]]$Hungaro),
                   accuracy=c(metricas_boot[[2]][[5]]$Maximo,metricas_boot[[2]][[5]]$Hungaro))
df<-rbind(df_nd, df_wnv)
nrow(df)

write.table(df, file = "boxplot_flowcap_flowMeansOpt.csv",
            sep = "\t", row.names = F)