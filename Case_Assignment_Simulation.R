#----------------------------------------------------------------------------------
# Simulation function for each case separately
#----------------------------------------------------------------------------------
Case_Assignment_Simulation <- function(SRT,
                                      Cases=NULL,
                                      Typing_Error=NULL,
                                      Allele_freq=NULL,
                                      Num_Test=1000,
                                      Prop_Loci_Type=NULL,
                                      Min_genotyped_STR= 1,
                                      One_Parent_Unknown_Another=TRUE,
                                      Father_Known_Mother=TRUE,
                                      Mother_Known_Father=TRUE,
                                      n.cores = parallel::detectCores() - 1){

  STR[] <- sapply(STR[],as.character)
  rownames(STR)=STR$ID
  # Allele frequency
  if(is.null(Allele_freq)){
    Allele_freq <- list()
    for (i in which(!colnames(STR)%in%c("ID","Sex", "DamID", "SireID", "Date","Gen"))) {
      p=table(unlist(strsplit(STR[,i], ""), use.names=FALSE))/sum(table(unlist(strsplit(STR[,i], ""), use.names=FALSE)))
      Allele_freq[[length(Allele_freq)+1]] <- list(p)
    }
    names(Allele_freq)<-colnames(STR)[which(!colnames(STR)%in%c("ID","Sex", "DamID", "SireID", "Date","Gen"))]
  }
  # Proportion of loci typed
  if(is.null(Prop_Loci_Type)){
    Prop_Loci_Type=round(colSums(!is.na(STR[,which(!colnames(STR)%in%c("ID","Sex", "DamID", "SireID", "Date","Gen"))]),na.rm = T)/nrow(STR),4)
  }
  # Remove individuals with less than Min_genotyped_STR genotyped
  STR=STR[rowSums(!is.na(STR[,which(!colnames(STR)%in%c("ID","Sex", "DamID", "SireID", "Date","Gen"))]))>=Min_genotyped_STR,]
  # calculating Typing_Error
  if(is.null(Typing_Error)){
    # load Typing_Error_Calculation function
    source("./Typing_Error_Calculation.R")
    Typing_Error=Typing_Error_Calculation(STR=STR)
  }
  # Calculate the P-value distribution for the cases with information for only mother
  # Find cases when it is null
  if(is.null(Cases)){
    Cases=STR[,c("ID","DamID","SireID")]
  }else{
    Cases=STR[STR$ID%in%Cases,c("ID","DamID","SireID")]
    Cases$DamID[!Cases$DamID%in%STR$ID]=NA_character_
    Cases$SireID[!Cases$SireID%in%STR$ID]=NA_character_
  }
  All_result=list()
  
  #------------------------------------------------------------------------------------
  # For father cases when the mother is unknown
  if(nrow(Cases)>0 & One_Parent_Unknown_Another){
    #------------------------------------------------------------------------
    registerDoParallel(n.cores)
    results=foreach (1:Num_Test, .combine=rbind) %dopar% {
      # Load Functions
      source("./LOD_Both_Unknown.R")
      #------------------------------------------------
      #LOD father when the mother is unknown
      LOD_all=NULL
      for (i in 1:nrow(Cases)) {
        STR_Offs=STR[STR$ID==Cases$ID[i],colnames(STR)[which(!colnames(STR)%in%c("ID","Sex", "DamID", "SireID", "Date","Gen"))]]
        ind2=sapply(STR_Offs,function(v) intToUtf8(utf8ToInt(v)[sort(sample(nchar(v), 1))]),USE.NAMES = FALSE)
        # All genotype to missed ones
        for (Sim_allele in which(nchar(ind2)==0)) {
          ind2[Sim_allele]=names(unlist(Allele_freq[[Sim_allele]]))[min(which(cumsum(unlist(Allele_freq[[Sim_allele]]))>runif(1)))]
        }
        ind1=NULL
        for (Sim_allele in 1:length(ind2)) {
          ind1[Sim_allele]=names(unlist(Allele_freq[[Sim_allele]]))[min(which(cumsum(unlist(Allele_freq[[Sim_allele]]))>runif(1)))]
        }
        STR_Father=paste0(ind1,ind2)
        STR_Father=data.frame(matrix(STR_Father, 1))
        colnames(STR_Father)=colnames(STR_Offs)
        # missing the genotype in individuals based on parameters
        ind_mis=which(runif(length(STR_Father))>Prop_Loci_Type)
        if(length(ind_mis)>0)STR_Father[ind_mis]=NA_character_
        # Typing_Error for each allele
        Error_allel2=which(runif(length(STR_Father))<(Typing_Error/2) & !is.na(STR_Father))
        if(length(Error_allel2)>0){
          for (Sim_allele in Error_allel2) {
            alleles=unlist(Allele_freq[[Sim_allele]])[-which(names(unlist(Allele_freq[[Sim_allele]]))%in%unlist(strsplit(as.character(STR_Offs[Sim_allele]),"")))]
            if(length(alleles)==0)alleles=unlist(Allele_freq[[Sim_allele]])
            STR_Father[Sim_allele]=paste0(unlist(strsplit(as.character(STR_Father[Sim_allele]),""))[1],
                                sample(names(alleles),1,prob =alleles))
          }
        }
        LOD_all=c(LOD_all,LOD_Both_Unknown(STR_Offs,STR_Father,Allele_freq,Typing_Error))
      }
      return(LOD_all) 
    }
    stopImplicitCluster()
    Cases_Fa_Un_Mo=results
    colnames(Cases_Fa_Un_Mo)=Cases$ID
    rownames(Cases_Fa_Un_Mo)=paste0("Sim",1:nrow(Cases_Fa_Un_Mo))
    All_result[["One_Parent_Unknown_Another"]]=Cases_Fa_Un_Mo
    cat("\nCases with one parent when the other is unknown were done!\n")
  }
  #------------------------------------------------------------------------------------
  # For father when mother is known
  Cases_father=Cases[!is.na(Cases$DamID),]
  if(nrow(Cases_father)>0 & Father_Known_Mother){
    registerDoParallel(n.cores)
    results=foreach (1:Num_Test, .combine=rbind) %dopar% {
      # Load Functions
      source("./LOD_One_Parent_Known.R")
      #------------------------------------------------
      #LOD father when the mother is unknown
      LOD_all=NULL
      for (i in 1:nrow(Cases_father)) {
        STR_Known_Parent=STR[STR$ID==Cases_father$DamID[i],colnames(STR)[which(!colnames(STR)%in%c("ID","Sex", "DamID", "SireID", "Date","Gen"))]]
        STR_Offs=STR[STR$ID==Cases_father$ID[i],colnames(STR)[which(!colnames(STR)%in%c("ID","Sex", "DamID", "SireID", "Date","Gen"))]]
        STR_Father=as.data.frame(matrix(rep(NA_character_,ncol(STR_Offs)),1))
        colnames(STR_Father)=colnames(STR_Offs)
        for (j in 1:ncol(STR_Father)) {
          if(!is.na(STR_Known_Parent[j]) & !is.na(STR_Offs[j])){
            allel_known=unlist(strsplit(as.character(STR_Known_Parent[j]),""))
            allel_off=unlist(strsplit(as.character(STR_Offs[j]),""))
            allel_off=allel_off[!allel_off%in%allel_known]
            if(length(allel_off)==0)allel_off=unlist(strsplit(as.character(STR_Offs[j]),""))
            STR_Father[j]=sample(allel_off,1)
          }
          if(is.na(STR_Known_Parent[j]) & !is.na(STR_Offs[j])){
            allel_off=unlist(strsplit(as.character(STR_Offs[j]),""))
            STR_Father[j]=sample(allel_off,1)
          }
          if(is.na(STR_Known_Parent[j]) & is.na(STR_Offs[j])){
            STR_Father[j]=sample(names(unlist(Allele_freq[[j]])),1,prob =unlist(Allele_freq[[j]]) )
          }
        }
        ind1=NULL
        for (Sim_allele in 1:ncol(STR_Offs)) {
          ind1[Sim_allele]=sample(names(unlist(Allele_freq[[Sim_allele]])),1,prob = unlist(Allele_freq[[Sim_allele]]))
        }
        STR_Father[1,]=paste0(ind1,STR_Father)
        # missing the genotype in individuals based on parameters
        ind_mis=which(runif(length(STR_Father))>Prop_Loci_Type)
        if(length(ind_mis)>0)STR_Father[ind_mis]=NA_character_
        # Typing_Error for each allele
        Error_allel2=which(runif(length(STR_Father))<(Typing_Error/2) & !is.na(STR_Father))
        if(length(Error_allel2)>0){
          for (Sim_allele in Error_allel2) {
            alleles=unlist(Allele_freq[[Sim_allele]])[-which(names(unlist(Allele_freq[[Sim_allele]]))%in%unlist(strsplit(as.character(STR_Offs[Sim_allele]),"")))]
            if(length(alleles)==0)alleles=unlist(Allele_freq[[Sim_allele]])
            STR_Father[Sim_allele]=paste0(unlist(strsplit(as.character(STR_Father[Sim_allele]),""))[1],
                                          sample(names(alleles),1,prob =alleles))
          }
        }
        LOD_all=c(LOD_all,LOD_One_Parent_Known(STR_Offs,STR_Father,STR_Known_Parent,Allele_freq,Typing_Error))
      }
      return(LOD_all) 
    }
    stopImplicitCluster()
    Cases_Fa_Kn_Mo=results
    colnames(Cases_Fa_Kn_Mo)=Cases_father$ID
    rownames(Cases_Fa_Kn_Mo)=paste0("Sim",1:nrow(Cases_Fa_Kn_Mo))
    All_result[["Cases_Fa_Kn_Mo"]]=Cases_Fa_Kn_Mo
    cat("\nCases for father when the mother is known were done!\n")
  }
  #------------------------------------------------------------------------------------
  # For mother when father is known
  Cases_mother=Cases[!is.na(Cases$SireID),]
  if(nrow(Cases_mother)>0 & Mother_Known_Father){
    registerDoParallel(n.cores)
    results=foreach (1:Num_Test, .combine=rbind) %dopar% {
      # Load Functions
      source("./LOD_One_Parent_Known.R")
      #------------------------------------------------
      #LOD father when the mother is unknown
      LOD_all=NULL
      for (i in 1:nrow(Cases_mother)) {
        STR_Known_Parent=STR[STR$ID==Cases_mother$SireID[i],colnames(STR)[which(!colnames(STR)%in%c("ID","Sex", "DamID", "SireID", "Date","Gen"))]]
        STR_Offs=STR[STR$ID==Cases_mother$ID[i],colnames(STR)[which(!colnames(STR)%in%c("ID","Sex", "DamID", "SireID", "Date","Gen"))]]
        STR_Mother=as.data.frame(matrix(rep(NA_character_,ncol(STR_Offs)),1))
        colnames(STR_Mother)=colnames(STR_Offs)
        for (j in 1:ncol(STR_Mother)) {
          if(!is.na(STR_Known_Parent[j]) & !is.na(STR_Offs[j])){
            allel_known=unlist(strsplit(as.character(STR_Known_Parent[j]),""))
            allel_off=unlist(strsplit(as.character(STR_Offs[j]),""))
            allel_off=allel_off[!allel_off%in%allel_known]
            if(length(allel_off)==0)allel_off=unlist(strsplit(as.character(STR_Offs[j]),""))
            STR_Mother[j]=sample(allel_off,1)
          }
          if(is.na(STR_Known_Parent[j]) & !is.na(STR_Offs[j])){
            allel_off=unlist(strsplit(as.character(STR_Offs[j]),""))
            STR_Mother[j]=sample(allel_off,1)
          }
          if(is.na(STR_Known_Parent[j]) & is.na(STR_Offs[j])){
            STR_Mother[j]=sample(names(unlist(Allele_freq[[j]])),1,prob =unlist(Allele_freq[[j]]) )
          }
        }
        ind1=NULL
        for (Sim_allele in 1:ncol(STR_Offs)) {
          ind1[Sim_allele]=names(unlist(Allele_freq[[Sim_allele]]))[min(which(cumsum(unlist(Allele_freq[[Sim_allele]]))>runif(1)))]
        }
        STR_Mother=paste0(ind1,STR_Mother)
        
        STR_Mother=data.frame(matrix(STR_Mother, 1))
        colnames(STR_Mother)=colnames(STR_Offs)
        # missing the genotype in individuals based on parameters
        ind_mis=which(runif(length(STR_Mother))>=Prop_Loci_Type)
        if(length(ind_mis)>0)STR_Mother[ind_mis]=NA_character_
        # Typing_Error for each allele
        Error_allel2=which(runif(length(STR_Mother))<(Typing_Error/2) & !is.na(STR_Mother))
        if(length(Error_allel2)>0){
          for (Sim_allele in Error_allel2) {
            alleles=unlist(Allele_freq[[Sim_allele]])[-which(names(unlist(Allele_freq[[Sim_allele]]))%in%unlist(strsplit(as.character(STR_Offs[Sim_allele]),"")))]
            if(length(alleles)==0)alleles=unlist(Allele_freq[[Sim_allele]])
            STR_Mother[Sim_allele]=paste0(unlist(strsplit(as.character(STR_Mother[Sim_allele]),""))[1],
                                          sample(names(alleles),1,prob =alleles))
          }
        }
        LOD_all=c(LOD_all,LOD_One_Parent_Known(STR_Offs,STR_Mother,STR_Known_Parent,Allele_freq,Typing_Error))
      }
      return(LOD_all) 
    }
    stopImplicitCluster()
    Cases_Mo_Kn_Fa=results
    colnames(Cases_Mo_Kn_Fa)=Cases_mother$ID
    rownames(Cases_Mo_Kn_Fa)=paste0("Sim",1:nrow(Cases_Mo_Kn_Fa))
    All_result[["Cases_Mo_Kn_Fa"]]=Cases_Mo_Kn_Fa
    cat("\nCases for mother when the father is known were done!\n")
  }
  return(All_result)
}




