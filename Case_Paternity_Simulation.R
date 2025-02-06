#----------------------------------------------------------------------------------
# Simulation function for each case separately
#----------------------------------------------------------------------------------
Case_Paternity_Simulation <- function(SRT,
                                      Cases=NULL,
                                      Typing_Error=NULL,
                                      Allele_freq=NULL,
                                      Num_Test=1000,
                                      Prop_Loci_Type=NULL,
                                      Min_genotyped_allele= 1,
                                      Father_Unknown_Mother=TRUE,
                                      Mother_Unknown_Father=TRUE,
                                      Father_Known_Mother=TRUE,
                                      Mother_Known_Father=TRUE,
                                      Both_Parent_Jointly=TRUE,
                                      n.cores = parallel::detectCores() - 1){

  STR[] <- sapply(STR[],as.character)
  rownames(STR)=STR$ID
  # Allele frequency
  if(is.null(Allele_freq)){
    Allele_freq <- list()
    for (i in which(!colnames(STR)%in%c("ID","Sex", "DamID", "SireID", "Date"))) {
      p=table(unlist(strsplit(STR[,i], ""), use.names=FALSE))/sum(table(unlist(strsplit(STR[,i], ""), use.names=FALSE)))
      Allele_freq[[length(Allele_freq)+1]] <- list(p)
    }
    names(Allele_freq)<-colnames(STR)[which(!colnames(STR)%in%c("ID","Sex", "DamID", "SireID", "Date"))]
  }
  # Proportion of loci typed
  if(is.null(Prop_Loci_Type)){
    Prop_Loci_Type=round(colSums(!is.na(STR[,which(!colnames(STR)%in%c("ID","Sex", "DamID", "SireID", "Date"))]),na.rm = T)/nrow(STR),4)
  }
  # Remove individuals with less than Min_genotyped_allele genotyped
  STR=STR[rowSums(!is.na(STR[,which(!colnames(STR)%in%c("ID","Sex", "DamID", "SireID", "Date"))]))>=Min_genotyped_allele,]
  # calculating Typing_Error
  if(is.null(Typing_Error)){
    # load Typing_Error_Calculation function
    source("F:/My_Projects/Paternity_algorithm/Functions/Typing_Error_Calculation.R")
    Typing_Error=Typing_Error_Calculation(STR=STR)
  }
  # Calculate the P-value distribution for the cases with information for only mother
  # Find cases when it is null
  if(is.null(Cases)){
    Cases=STR[,c("ID","DamID","SireID")]
    Cases$DamID[!Cases$DamID%in%STR$ID]=NA_character_
    Cases$SireID[!Cases$SireID%in%STR$ID]=NA_character_
  }else{
    Cases=STR[STR$ID%in%Cases,c("ID","DamID","SireID")]
    Cases$DamID[!Cases$DamID%in%STR$ID]=NA_character_
    Cases$SireID[!Cases$SireID%in%STR$ID]=NA_character_
  }
  All_result=list()
  
  #------------------------------------------------------------------------------------
  # For father cases when the mother is unknown
  Cases_Father=Cases[!is.na(Cases$SireID),]
  if(nrow(Cases_Father)>0 & Father_Unknown_Mother){
    #------------------------------------------------------------------------
    registerDoParallel(n.cores)
    results=foreach (1:Num_Test, .combine=rbind) %dopar% {
      # Load Functions
      source("F:/My_Projects/Paternity_algorithm/Functions/LOD_Both_Unknown.R")
      #------------------------------------------------
      #LOD father when the mother is unknown
      LOD_all=NULL
      for (i in 1:nrow(Cases_Father)) {
        STR_Offs=STR[STR$ID==Cases_Father$ID[i],names(Typing_Error)]
        STR_father=STR[STR$ID==Cases_Father$SireID[i],names(Typing_Error)]
        STR_Candidate_offs=as.data.frame(matrix(rep(NA_character_,ncol(STR_Offs)),1))
        colnames(STR_Candidate_offs)=colnames(STR_Offs)
        for (j in 1:ncol(STR_father)) {
          if(!is.na(STR_father[j])){
            STR_Candidate_offs[j]=paste0(sample(unlist(strsplit(as.character(STR_father[j]),"")),1),
                                         sample(names(unlist(Allele_freq[[j]])),1,prob = unlist(Allele_freq[[j]])))
          }
          if(is.na(STR_father[j])){
            STR_Candidate_offs[j]=paste0(sample(names(unlist(Allele_freq[[j]])),1,prob = unlist(Allele_freq[[j]])),
                                         sample(names(unlist(Allele_freq[[j]])),1,prob = unlist(Allele_freq[[j]])))
          }
        }
        # Typing_Error for each allele
        Error_allel1=which(runif(length(STR_father))<(Typing_Error/2) & !is.na(STR_father))
        if(length(Error_allel1)>0){
          for (Sim_allele in Error_allel1) {
            alleles=unlist(Allele_freq[[Sim_allele]])[-which(names(unlist(Allele_freq[[Sim_allele]]))%in%unlist(strsplit(as.character(STR_father[Sim_allele]),"")))]
            if(length(alleles)==0)alleles=unlist(Allele_freq[[Sim_allele]])
            STR_Candidate_offs[Sim_allele]=paste0(sample(names(alleles),1,prob =alleles),
                                                  unlist(strsplit(as.character(STR_Candidate_offs[Sim_allele]),""))[2])
          }
        }
        STR_Candidate_offs[is.na(STR_Offs)]=NA_character_
        # calculating LOD
        LOD_all=c(LOD_all,LOD_Both_Unknown(STR_Candidate_offs,STR_father,Allele_freq,Typing_Error))
      }
      return(LOD_all) 
    }
    stopImplicitCluster()
    Cases_Fa_Un_Mo=results
    colnames(Cases_Fa_Un_Mo)=Cases_Father$ID
    rownames(Cases_Fa_Un_Mo)=paste0("Sim",1:nrow(Cases_Fa_Un_Mo))
    All_result[["Cases_Fa_Un_Mo"]]=Cases_Fa_Un_Mo
    cat("\nCases with fathers when the mother is unknown were done!\n")
  }
  #------------------------------------------------------------------------------------
  # For father cases when the mother is unknown
  Cases_Mother=Cases[!is.na(Cases$DamID),]
  if(nrow(Cases_Mother)>0 & Mother_Unknown_Father){
    #------------------------------------------------------------------------
    registerDoParallel(n.cores)
    results=foreach (1:Num_Test, .combine=rbind) %dopar% {
      # Load Functions
      source("F:/My_Projects/Paternity_algorithm/Functions/LOD_Both_Unknown.R")
      #------------------------------------------------
      #LOD father when the mother is unknown
      LOD_all=NULL
      for (i in 1:nrow(Cases_Mother)) {
        STR_Offs=STR[STR$ID==Cases_Mother$ID[i],names(Typing_Error)]
        STR_mother=STR[STR$ID==Cases_Mother$DamID[i],names(Typing_Error)]
        STR_Candidate_offs=as.data.frame(matrix(rep(NA_character_,ncol(STR_Offs)),1))
        colnames(STR_Candidate_offs)=colnames(STR_Offs)
        for (j in 1:ncol(STR_mother)) {
          if(!is.na(STR_mother[j])){
            STR_Candidate_offs[j]=paste0(sample(unlist(strsplit(as.character(STR_mother[j]),"")),1),
                                         sample(names(unlist(Allele_freq[[j]])),1,prob = unlist(Allele_freq[[j]])))
          }
          if(is.na(STR_mother[j])){
            STR_Candidate_offs[j]=paste0(sample(names(unlist(Allele_freq[[j]])),1,prob = unlist(Allele_freq[[j]])),
                                         sample(names(unlist(Allele_freq[[j]])),1,prob = unlist(Allele_freq[[j]])))
          }
        }
        # Typing_Error for each allele
        Error_allel1=which(runif(length(STR_mother))<(Typing_Error/2) & !is.na(STR_mother))
        if(length(Error_allel1)>0){
          for (Sim_allele in Error_allel1) {
            alleles=unlist(Allele_freq[[Sim_allele]])[-which(names(unlist(Allele_freq[[Sim_allele]]))%in%unlist(strsplit(as.character(STR_mother[Sim_allele]),"")))]
            if(length(alleles)==0)alleles=unlist(Allele_freq[[Sim_allele]])
            STR_Candidate_offs[Sim_allele]=paste0(sample(names(alleles),1,prob =alleles),
                                                  unlist(strsplit(as.character(STR_Candidate_offs[Sim_allele]),""))[2])
          }
        }
        STR_Candidate_offs[is.na(STR_Offs)]=NA_character_
        # calculating LOD
        LOD_all=c(LOD_all,LOD_Both_Unknown(STR_Candidate_offs,STR_mother,Allele_freq,Typing_Error))
      }
      return(LOD_all) 
    }
    stopImplicitCluster()
    Cases_Mo_Un_Fa=results
    colnames(Cases_Mo_Un_Fa)=Cases_Mother$ID
    rownames(Cases_Mo_Un_Fa)=paste0("Sim",1:nrow(Cases_Mo_Un_Fa))
    All_result[["Cases_Mo_Un_Fa"]]=Cases_Mo_Un_Fa
    cat("\nCases with mothers when the father is unknown were done!\n")
  }
  #------------------------------------------------------------------------------------
  # For parents when both parents are known
  Cases_Both=Cases[!is.na(Cases$DamID) & !is.na(Cases$SireID),]
  if(nrow(Cases_Both)>0 & any(c(Father_Known_Mother, Mother_Known_Father, Both_Parent_Jointly))){
    #------------------------------------------------------------------------
    registerDoParallel(n.cores)
    results=foreach (1:Num_Test, .combine=rbind) %dopar% {
      # Load Functions
      source("F:/My_Projects/Paternity_algorithm/Functions/LOD_One_Parent_Known.R")
      source("F:/My_Projects/Paternity_algorithm/Functions/LOD_Parents_Jointly.R")
      #------------------------------------------------
      #LOD father when the mother is unknown
      LOD_Father=NULL
      LOD_Mother=NULL
      LOD_Jointly=NULL
      LOD_all=NULL
      for (i in 1:nrow(Cases_Both)) {
        STR_Offs=STR[STR$ID==Cases_Both$ID[i],names(Typing_Error)]
        STR_mother=STR[STR$ID==Cases_Both$DamID[i],names(Typing_Error)]
        STR_father=STR[STR$ID==Cases_Both$SireID[i],names(Typing_Error)]
        STR_Candidate_offs=as.data.frame(matrix(rep(NA_character_,ncol(STR_Offs)),1))
        colnames(STR_Candidate_offs)=colnames(STR_Offs)
        for (j in 1:ncol(STR_Candidate_offs)) {
          if(!is.na(STR_father[j]) & !is.na(STR_mother[j])){
            STR_Candidate_offs[j]=paste0(sample(unlist(strsplit(as.character(STR_father[j]),"")),1),
                                         sample(unlist(strsplit(as.character(STR_mother[j]),"")),1))
          }
          if(!is.na(STR_father[j]) & is.na(STR_mother[j])){
            STR_Candidate_offs[j]=paste0(sample(unlist(strsplit(as.character(STR_father[j]),"")),1),
                                         sample(names(unlist(Allele_freq[[j]])),1,prob = unlist(Allele_freq[[j]])))
          }
          if(is.na(STR_father[j]) & !is.na(STR_mother[j])){
            STR_Candidate_offs[j]=paste0(sample(names(unlist(Allele_freq[[j]])),1,prob = unlist(Allele_freq[[j]])),
                                         sample(unlist(strsplit(as.character(STR_mother[j]),"")),1))
          }
          if(is.na(STR_father[j]) & is.na(STR_mother[j])){
            STR_Candidate_offs[j]=paste0(sample(names(unlist(Allele_freq[[j]])),1,prob = unlist(Allele_freq[[j]])),
                                         sample(names(unlist(Allele_freq[[j]])),1,prob = unlist(Allele_freq[[j]])))
          }
        }
        # Typing_Error for each allele
        Error_allel1=which(runif(length(STR_father))<(Typing_Error/2) & !is.na(STR_father))
        if(length(Error_allel1)>0){
          for (Sim_allele in Error_allel1) {
            alleles=unlist(Allele_freq[[Sim_allele]])[-which(names(unlist(Allele_freq[[Sim_allele]]))%in%unlist(strsplit(as.character(STR_father[Sim_allele]),"")))]
            if(length(alleles)==0)alleles=unlist(Allele_freq[[Sim_allele]])
            STR_Candidate_offs[Sim_allele]=paste0(sample(names(alleles),1,prob =alleles),
                                                  unlist(strsplit(as.character(STR_Candidate_offs[Sim_allele]),""))[2])
          }
        }
        Error_allel2=which(runif(length(STR_mother))<(Typing_Error/2) & !is.na(STR_mother))
        if(length(Error_allel2)>0){
          for (Sim_allele in Error_allel2) {
            alleles=unlist(Allele_freq[[Sim_allele]])[-which(names(unlist(Allele_freq[[Sim_allele]]))%in%unlist(strsplit(as.character(STR_mother[Sim_allele]),"")))]
            if(length(alleles)==0)alleles=unlist(Allele_freq[[Sim_allele]])
            STR_Candidate_offs[Sim_allele]=paste0(unlist(strsplit(as.character(STR_Candidate_offs[Sim_allele]),""))[1],
                                                  sample(names(alleles),1,prob =alleles))
          }
        }
        STR_Candidate_offs[is.na(STR_Offs)]=NA_character_
        # calculating LOD
        if(Father_Known_Mother)LOD_Father=c(LOD_Father,LOD_One_Parent_Known(STR_Candidate_offs,STR_father,STR_mother,Allele_freq,Typing_Error))
        if(Mother_Known_Father)LOD_Mother=c(LOD_Mother,LOD_One_Parent_Known(STR_Candidate_offs,STR_father,STR_mother,Allele_freq,Typing_Error))
        if(Both_Parent_Jointly)LOD_Jointly=c(LOD_Jointly,LOD_One_Parent_Known(STR_Candidate_offs,STR_father,STR_mother,Allele_freq,Typing_Error))
      }
      LOD_all=c(LOD_Father,LOD_Mother,LOD_Jointly)
      return(LOD_all) 
    }
    stopImplicitCluster()
    if(Father_Known_Mother){
      Cases_Fa_Kn_Mo=results[,1:nrow(Cases_Both)]
      colnames(Cases_Fa_Kn_Mo)=Cases_Both$ID
      rownames(Cases_Fa_Kn_Mo)=paste0("Sim",1:nrow(Cases_Fa_Kn_Mo))
      All_result[["Cases_Fa_Kn_Mo"]]=Cases_Fa_Kn_Mo
    }
    if(Mother_Known_Father){
      Cases_Mo_Kn_Fa=results[,(1+(sum(Father_Known_Mother)*nrow(Cases_Both))):(nrow(Cases_Both)*(sum(Father_Known_Mother)+1))]
      colnames(Cases_Mo_Kn_Fa)=Cases_Both$ID
      rownames(Cases_Mo_Kn_Fa)=paste0("Sim",1:nrow(Cases_Mo_Kn_Fa))
      All_result[["Cases_Mo_Kn_Fa"]]=Cases_Mo_Kn_Fa
    }
    if(Both_Parent_Jointly){
      Cases_Both_Jointly=results[,(1+(sum(c(Father_Known_Mother,Mother_Known_Father))*nrow(Cases_Both))):(nrow(Cases_Both)*(sum(c(Father_Known_Mother,Mother_Known_Father))+1))]
      colnames(Cases_Both_Jointly)=Cases_Both$ID
      rownames(Cases_Both_Jointly)=paste0("Sim",1:nrow(Cases_Both_Jointly))
      All_result[["Cases_Both_Jointly"]]=Cases_Both_Jointly
    }
    cat("\nCases with mothers when the father is unknown were done!\n")
  }
  return(All_result)
}



