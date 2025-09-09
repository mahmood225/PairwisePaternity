#----------------------------------------------------------------------------------
# function to calculate laboratory typing errors
#----------------------------------------------------------------------------------
Typing_Error_Calculation <- function(STR,
                                     Min_genotyped_STR= 7,
                                     Correct_Genotype_Probability = 0.8,
                                     Minimum_Typing_Error=0.001,
                                     Ensure_Correctness=FALSE,
                                     Mother_Offspring_Calulation=TRUE,
                                     Father_Offspring_Calulation=TRUE){
  STR[] <- sapply(STR[],as.character)
  STR=STR[,c("ID", "Sex", "DamID","SireID",colnames(STR)[which(!colnames(STR)%in%c("ID", "Sex", "DamID","SireID", "Date","Gen"))])]
  STR_Names=colnames(STR)[which(!colnames(STR)%in%c("ID","Sex", "DamID", "SireID","Gen"))]
  if(Father_Offspring_Calulation){
    STR=STR[rowSums(!is.na(STR[,STR_Names]))>=1,]
  }else{
    STR=STR[rowSums(!is.na(STR[,STR_Names]))>=Min_genotyped_STR,]
  }
  # mother offspring calculation
  if(Mother_Offspring_Calulation){
    STR_mother=STR[complete.cases(STR[,c("DamID")]),]
    STR_mother=STR_mother[STR_mother$DamID%in%STR$ID,]
    for (i in 1:nrow(STR_mother)) {
      STR_dam=STR[which(STR$ID==STR_mother$DamID[i]),STR_Names]
      STR_offs=STR_mother[i,STR_Names]
      similarity=stringdist::stringsim(toupper(STR_offs), toupper(STR_dam),method ="qgram")*2
      STR_mother[i,colnames(STR_offs)]=as.character(ifelse(similarity==0,1,0))
    }
    STR_mother[,STR_Names]=sapply(STR_mother[,which(!colnames(STR_mother)%in%c("ID","Sex", "DamID", "SireID"))],as.numeric)
    STR_mother=STR_mother[rowSums(!is.na(STR_mother[,STR_Names]))>=Min_genotyped_STR,]
    if(!Ensure_Correctness){
      STR_mother=STR_mother[1-(rowSums(STR_mother[,STR_Names],na.rm = T)/rowSums(!is.na(STR_mother[,STR_Names])))>=0.9,]
    }
  }
  # father offspring calculation
  if(Father_Offspring_Calulation){
    STR_father=STR[complete.cases(STR[,c("SireID")]),]
    STR_father=STR_father[STR_father$SireID%in%STR$ID,]
    for (i in 1:nrow(STR_father)) {
      STR_sire=STR[which(STR$ID==STR_father$SireID[i]),STR_Names]
      STR_offs=STR_father[i,STR_Names]
      similarity=stringdist::stringsim(toupper(STR_offs), toupper(STR_sire),method ="qgram")*2
      STR_father[i,colnames(STR_offs)]=as.character(ifelse(similarity==0,1,0))
    }
    STR_father[,STR_Names]=sapply(STR_father[,STR_Names],as.numeric)
    STR_father=STR_father[rowSums(!is.na(STR_father[,STR_Names]))>=Min_genotyped_STR,]
    if(!Ensure_Correctness){
      STR_father=STR_father[1-(rowSums(STR_father[,STR_Names],na.rm = T)/rowSums(!is.na(STR_father[,STR_Names])))>=0.9,]
    }
  }

  if(Father_Offspring_Calulation & Mother_Offspring_Calulation){
    Type_Error_mother=round(colSums(STR_mother[,STR_Names],na.rm = T)/
                              (colSums(!is.na(STR_mother[,STR_Names]),na.rm = T)*2),4)
    Type_Error_father=round(colSums(STR_father[,STR_Names],na.rm = T)/
                              (colSums(!is.na(STR_father[,STR_Names]),na.rm = T)*2),4)
    Type_Error=Type_Error_mother+Type_Error_father
  }
  if(Father_Offspring_Calulation & !Mother_Offspring_Calulation){
    Type_Error=round(colSums(STR_father[,STR_Names],na.rm = T)/
                       (colSums(!is.na(STR_father[,STR_Names]),na.rm = T)*2),4)
    Type_Error=Type_Error*2
  }
  if(!Father_Offspring_Calulation & Mother_Offspring_Calulation){
    Type_Error=round(colSums(STR_mother[,STR_Names],na.rm = T)/
                       (colSums(!is.na(STR_mother[,STR_Names]),na.rm = T)*2),4)
    Type_Error=Type_Error*2
  }
  if(min(Type_Error)<Minimum_Typing_Error){
    Minimum_Typing_Error=Minimum_Typing_Error-min(Type_Error)
    Type_Error=Type_Error+Minimum_Typing_Error
  }
  
  if(Father_Offspring_Calulation){
    # for father
    STR_father2=STR_father[rowSums(STR_father[,STR_Names],na.rm = T)>0,]
    STR_father2=STR_father2[STR_father2$ID%in%STR_mother$ID,]
    if(nrow(STR_father2)>0){
      gradient_father=as.data.frame(matrix(NA,ncol = length(STR_Names),nrow = nrow(STR_father2)))
      colnames(gradient_father)=STR_Names
      rownames(gradient_father)=STR_father2$ID
      for (i in 1:nrow(STR_father2)) {
        errors=colnames(STR_father2[i,STR_Names]>0)[which(STR_father2[i,STR_Names]>0)]
        for (j in errors) {
          STR_sire=tolower(STR[which(STR$ID==STR_father2$SireID[i]),j])
          STR_dam=tolower(STR[which(STR$ID==STR_father2$DamID[i]),j])
          STR_offs=tolower(STR[which(STR$ID==STR_father2$ID[i]),j])
          allel_off=unlist(strsplit(as.character(STR_offs),""))
          allel_mo=unlist(strsplit(as.character(STR_dam),""))
          allel_fa=unlist(strsplit(as.character(STR_sire),""))
          if(allel_off[1]==allel_off[2] & allel_fa[1]==allel_fa[2]){
            gradient_father[i,j]=abs(which(letters==allel_fa[1])-which(letters==allel_off[1]))
          }else{
            if(!is.na(STR_dam)){
              allel_off=allel_off[!allel_off%in%allel_mo]
              if(length(allel_off)==0)allel_off=unique(unlist(strsplit(as.character(STR_offs),"")))
              if(length(allel_off)==1){
                gradient_father[i,j]=mean(abs(which(letters%in%allel_off)-which(letters%in%allel_fa)))
              }else{
                gradient_father[i,j]=mean(c(mean(abs(which(letters%in%allel_off[2])-which(letters%in%allel_fa))),
                                            mean(abs(which(letters%in%allel_off[1])-which(letters%in%allel_fa)))))
              }
            }else{
              gradient_father[i,j]=mean(c(mean(abs(which(letters%in%allel_off[2])-which(letters%in%allel_fa))),
                                          mean(abs(which(letters%in%allel_off[1])-which(letters%in%allel_fa)))))
            }
          }
        }
      }
    }
  }
  if(Mother_Offspring_Calulation){
    STR_mother2=STR_mother[rowSums(STR_mother[,STR_Names],na.rm = T)>0,]
    STR_mother2=STR_mother2[STR_mother2$ID%in%STR_father$ID,]
    if(nrow(STR_mother2)>0){
      gradient_mother=as.data.frame(matrix(NA,ncol = length(STR_Names),nrow = nrow(STR_mother2)))
      colnames(gradient_mother)=STR_Names
      rownames(gradient_mother)=STR_mother2$ID
      for (i in 1:nrow(STR_mother2)) {
        errors=colnames(STR_mother2[i,STR_Names]>0)[which(STR_mother2[i,STR_Names]>0)]
        for (j in errors) {
          STR_sire=tolower(STR[which(STR$ID==STR_mother2$SireID[i]),j])
          STR_dam=tolower(STR[which(STR$ID==STR_mother2$DamID[i]),j])
          STR_offs=tolower(STR[which(STR$ID==STR_mother2$ID[i]),j])
          allel_off=unlist(strsplit(as.character(STR_offs),""))
          allel_mo=unlist(strsplit(as.character(STR_dam),""))
          allel_fa=unlist(strsplit(as.character(STR_sire),""))
          if(allel_off[1]==allel_off[2] & allel_mo[1]==allel_mo[2]){
            gradient_mother[i,j]=abs(which(letters==allel_mo[1])-which(letters==allel_off[1]))
          }else{
            if(!is.na(STR_sire)){
              allel_off=allel_off[!allel_off%in%allel_fa]
              if(length(allel_off)==0)allel_off=unique(unlist(strsplit(as.character(STR_offs),"")))
              if(length(allel_off)==1){
                gradient_mother[i,j]=mean(abs(which(letters%in%allel_off)-which(letters%in%allel_mo)))
              }else{
                gradient_mother[i,j]=mean(c(mean(abs(which(letters%in%allel_off[2])-which(letters%in%allel_mo))),
                                            mean(abs(which(letters%in%allel_off[1])-which(letters%in%allel_mo)))))
              }
            }else{
              gradient_mother[i,j]=mean(c(mean(abs(which(letters%in%allel_off[2])-which(letters%in%allel_mo))),
                                          mean(abs(which(letters%in%allel_off[1])-which(letters%in%allel_mo)))))
            }
          }
        }
      }
    }
  }
  return(Type_Error) 
}
