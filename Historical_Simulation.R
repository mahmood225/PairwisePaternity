#----------------------------------------------------------------------------------
# Simulation function Using simulating of historical and current population
#----------------------------------------------------------------------------------
Historical_Simulation <- function(
  # STR parameters
  Num_SRT=15,
  Allele_num=c(rep(8,15)),
  Allele_num_random=F,
  Min_allele_num=4,
  Max_allele_num=15,
  Typing_Error=c(rep(0.01,15)),
  Mutation_Rate=0.00001,
  # general parameters
  Overlap_Gen_Num=2,
  Overlap_Gen_Prop=c(0.6,0.3,0.1),
  Sire_Prop=0.5,
  Dam_Prop=0.8,
  Num_Offs=c(1,2),
  Prop_Num_Offs=c(0.95,0.05),
  Sex_Prop=c(0.5,0.5),# first male and second female
  # historical parameters
  Num_Hist_Start=2000,
  Hist_Gen=50,
  Num_Hist_End=200,
  # Population Parameters
  Num_Pop_Start=200,
  Pop_Gen=10,
  Num_Pop_End=200,
  n.cores = parallel::detectCores() - 1){
  # check if the number of Num_Pop_Start isn't more than Num_Hist_End
  if(Num_Pop_Start>Num_Hist_End) { stop("Num_Pop_Start is more than Num_Hist_End\n")}
  # check if the number of Num_Pop_Start isn't more than Num_Hist_End
  if(!is.null(Typing_Error) & length(Typing_Error)!=Num_SRT) { stop("number of STR is not equal to the Typing_Error's length\n")}
  # calculating of population size in each generation for historical population
  Pop_Siz_Gen_His=seq(Num_Hist_Start,Num_Hist_End,length.out=Hist_Gen+1)
  ###---------------------------------------------------------------
  ### Historical population
  ###---------------------------------------------------------------
  ##-------------------------------------------
  cat("Simulating of the historical population started!!! \n\n")
  ## Simulating the first Historical population
  # ped
  ped=data.frame(ID=1:Pop_Siz_Gen_His[1],
                 DamID=0,
                 SireID=0,
                 Sex=if(Pop_Siz_Gen_His[1]%%2==0) {
                   rep(c("M","F"),each=Pop_Siz_Gen_His[1]/2)} else {
                     c(rep("M",each=floor(Pop_Siz_Gen_His[1]/2)),rep("F",each=ceiling(Pop_Siz_Gen_His[1]/2)))},
                 Gen=0)
  # STR
  if(Allele_num_random){Allele_num=round(runif(Num_SRT,min=Min_allele_num-0.49,max=Max_allele_num+0.49))}
  STR=as.data.frame(matrix(NA,nrow = nrow(ped),ncol = Num_SRT))
  rownames(STR)=1:Pop_Siz_Gen_His[1]
  for (i in 1:Num_SRT) {
    alleles=letters[1:Allele_num[i]]
    gamet1=rep(alleles,each=floor(nrow(STR)/Allele_num[i]))
    if(length(gamet1)<nrow(STR))gamet1=c(gamet1,sample(alleles,nrow(STR)-length(gamet1)))
    gamet2=gamet1[sample(1:length(gamet1),length(gamet1),replace = FALSE)]
    STR[,i]=paste0(gamet1,gamet2)
  }
  ##-------------------------------------------
  ## simulating the next generations
  for (i in 1:Hist_Gen) {
    # compute the number of sire and Dams needed to be sampled
    if(i>(Overlap_Gen_Num+1)) ggg=3 else ggg=i
    numb_selected_Dam=round(Overlap_Gen_Prop[1:ggg]/sum(Overlap_Gen_Prop[1:ggg])*((Pop_Siz_Gen_His[i]/2)*Dam_Prop))
    numb_selected_Sire=round(Overlap_Gen_Prop[1:ggg]/sum(Overlap_Gen_Prop[1:ggg])*((Pop_Siz_Gen_His[i]/2)*Sire_Prop))
    # sample the sire and Dams
    Dam_gen1=sample(ped$ID[ped$Sex=="F" & ped$Gen==i-1],numb_selected_Dam[1])
    if(length(numb_selected_Dam)>1){
      Dam_gen2=NULL
      for(j in 2:length(numb_selected_Dam)){
        Dam_ind=ped$ID[ped$Gen==i-j & ped$Sex=="F"]
        Dam_ind=Dam_ind[Dam_ind%in%ped$DamID[ped$Gen==(i-1)]]
        Dam_gen2=c(Dam_gen2,sample(Dam_ind,numb_selected_Dam[j]))
      }
      selected_Dam=c(Dam_gen1,Dam_gen2)
    }else selected_Dam=Dam_gen1
    Sire_gen1=sample(ped$ID[ped$Sex=="M"&ped$Gen==i-1],numb_selected_Sire[1])
    if(length(numb_selected_Sire)>1){
      Sire_gen2=NULL
      for(j in 2:length(numb_selected_Sire)){
        Sire_ind=ped$ID[ped$Gen==i-j & ped$Sex=="M"]
        Sire_ind=Sire_ind[Sire_ind%in%ped$SireID[ped$Gen==(i-1)]]
        Sire_gen2=c(Sire_gen2,sample(Sire_ind,numb_selected_Sire[j]))
      }
      selected_Sire=c(Sire_gen1,Sire_gen2)
    }else selected_Sire=Sire_gen1
    #compute the number of offspring number
    Num_Offs_Num=round(Pop_Siz_Gen_His[i+1]*Prop_Num_Offs)
    if(any(Num_Offs_Num%%Num_Offs!=0)){
      dd=which(Num_Offs_Num%%Num_Offs!=0)
      for (hh in dd) {
        Num_Offs_Num[1]=Num_Offs_Num[1]+Num_Offs_Num[hh]%%Num_Offs[hh]
        Num_Offs_Num[hh]=Num_Offs_Num[hh]-Num_Offs_Num[hh]%%Num_Offs[hh]
      }
    }
    if(sum(Num_Offs_Num)!=Pop_Siz_Gen_His[i+1]) Num_Offs_Num[1]=Num_Offs_Num[1]-(sum(Num_Offs_Num)-Pop_Siz_Gen_His[i+1])
    # samples mother and father
    ind=data.frame(DamID=rep(NA,Pop_Siz_Gen_His[i+1]),SireID=rep(NA,Pop_Siz_Gen_His[i+1]))
    for (ij in 1:length(Num_Offs_Num)) {
      # when ij is 1
      if(ij==1){
        # select for Dam
        if(length(selected_Dam)<Num_Offs_Num[ij]){
          ind$DamID[1:Num_Offs_Num[ij]]=c(selected_Dam,sample(selected_Dam,Num_Offs_Num[ij]-length(selected_Dam),replace =T))
        }else{
          ind$DamID[1:Num_Offs_Num[ij]]=sample(selected_Dam,Num_Offs_Num[ij],replace =F)
        }
        # select for sire
        if(length(selected_Sire)<Num_Offs_Num[ij]){
          ind$SireID[1:Num_Offs_Num[ij]]=c(selected_Sire,sample(selected_Sire,Num_Offs_Num[ij]-length(selected_Sire),replace =T))
        }else{
          ind$SireID[1:Num_Offs_Num[ij]]=sample(selected_Sire,Num_Offs_Num[ij],replace =F)
        }
        # when ij is more than 1
      }else{
        ind$DamID[(Num_Offs_Num[ij-1]+1):sum(Num_Offs_Num[1:ij])]=rep(sample(selected_Dam,Num_Offs_Num[ij]/Num_Offs[ij]),each=Num_Offs[ij])
        ind$SireID[(Num_Offs_Num[ij-1]+1):sum(Num_Offs_Num[1:ij])]=rep(sample(selected_Sire,Num_Offs_Num[ij]/Num_Offs[ij]),each=Num_Offs[ij])
      }
    }
    ped2=data.frame(ID=(max(ped$ID)+1):(Pop_Siz_Gen_His[i+1]+max(ped$ID)),
                    DamID=ind$DamID,
                    SireID=ind$SireID,
                    Sex=if(Pop_Siz_Gen_His[i+1]%%2==0) {
                      rep(c("M","F"),each=Pop_Siz_Gen_His[i+1]/2)} else {
                        c(rep("M",each=floor(Pop_Siz_Gen_His[i+1]/2)),rep("F",each=ceiling(Pop_Siz_Gen_His[i+1]/2)))},
                    Gen=i)
  ped=rbind(ped,ped2)
  ###---------------------------
  ### simulating the STR
  IDss=ped$ID[ped$Gen==i]
  registerDoParallel(n.cores)
  # Start simulation
  STR2=foreach (str_sim=1:length(IDss), .combine=rbind) %dopar% {
    STR_Dam=STR[rownames(STR)==as.character(ped$DamID[ped$ID==IDss[str_sim]]),]
    STR_Sire=STR[rownames(STR)==as.character(ped$SireID[ped$ID==IDss[str_sim]]),]
    Temp=paste0(sapply(STR_Sire,function(v) intToUtf8(utf8ToInt(v)[sort(sample(nchar(v), 1))]),USE.NAMES = FALSE),
                          sapply(STR_Dam,function(v) intToUtf8(utf8ToInt(v)[sort(sample(nchar(v), 1))]),USE.NAMES = FALSE))
    return(Temp)
  }
  stopImplicitCluster()
  rownames(STR2)=IDss
  STR2=as.data.frame(STR2)
  STR=rbind(STR,STR2)
  ped=ped[ped$Gen%in%c((i-Overlap_Gen_Num):i),]
  STR=STR[rownames(STR)%in%ped$ID,]
  cat(paste0("Historical generation ",i," is done! \n"))
  }
  
  rm(ped2,STR2)
  ###---------------------------------------------------------------
  ### Recent population
  ###---------------------------------------------------------------
  cat("Simulating of the recent population started!!! \n\n")
  ped=ped[ped$Gen==max(ped$Gen),]
  ped=ped[sample(1:nrow(ped),Num_Pop_Start, replace = F),]
  ped$Sex[1:round(Num_Pop_Start/2)]="M"
  ped$Sex[(round(Num_Pop_Start/2)+1):nrow(ped)]="F"
  ped$Gen=0
  STR=STR[rownames(STR)%in%as.character(ped$ID),]
  STR=STR[as.character(ped$ID),]
  ped$ID=1:nrow(ped)
  rownames(STR)=1:nrow(ped)
  ped$DamID=0;ped$SireID=0
  # calculating of population size in each generation for historical population
  Pop_Siz_Gen_pop=seq(Num_Pop_Start,Num_Pop_End,length.out=Pop_Gen+1)
  ##-------------------------------------------
  ## simulating the next generations
  for (i in 1:Pop_Gen) {
    # compute the number of sire and Dams needed to be sampled
    if(i>(Overlap_Gen_Num+1)) ggg=Overlap_Gen_Num+1 else ggg=i
    numb_selected_Dam=round(Overlap_Gen_Prop[1:ggg]/sum(Overlap_Gen_Prop[1:ggg])*((Pop_Siz_Gen_pop[i]/2)*Dam_Prop))
    numb_selected_Sire=round(Overlap_Gen_Prop[1:ggg]/sum(Overlap_Gen_Prop[1:ggg])*((Pop_Siz_Gen_pop[i]/2)*Sire_Prop))
    # sample the sire and Dams
    Dam_gen1=sample(ped$ID[ped$Sex=="F" & ped$Gen==i-1],numb_selected_Dam[1])
    if(length(numb_selected_Dam)>1){
      Dam_gen2=NULL
      for(j in 2:length(numb_selected_Dam)){
        Dam_ind=ped$ID[ped$Gen==i-j & ped$Sex=="F"]
        Dam_ind=Dam_ind[Dam_ind%in%ped$DamID[ped$Gen==(i-1)]]
        Dam_gen2=c(Dam_gen2,sample(Dam_ind,numb_selected_Dam[j]))
      }
      selected_Dam=c(Dam_gen1,Dam_gen2)
    }else selected_Dam=Dam_gen1
    Sire_gen1=sample(ped$ID[ped$Sex=="M"&ped$Gen==i-1],numb_selected_Sire[1])
    if(length(numb_selected_Sire)>1){
      Sire_gen2=NULL
      for(j in 2:length(numb_selected_Sire)){
        Sire_ind=ped$ID[ped$Gen==i-j & ped$Sex=="M"]
        Sire_ind=Sire_ind[Sire_ind%in%ped$SireID[ped$Gen==(i-1)]]
        Sire_gen2=c(Sire_gen2,sample(Sire_ind,numb_selected_Sire[j]))
      }
      selected_Sire=c(Sire_gen1,Sire_gen2)
    }else selected_Sire=Sire_gen1
    #compute the number of offspring number
    Num_Offs_Num=round(Pop_Siz_Gen_pop[i+1]*Prop_Num_Offs)
    if(any(Num_Offs_Num%%Num_Offs!=0)){
      dd=which(Num_Offs_Num%%Num_Offs!=0)
      for (hh in dd) {
        Num_Offs_Num[1]=Num_Offs_Num[1]+Num_Offs_Num[hh]%%Num_Offs[hh]
        Num_Offs_Num[hh]=Num_Offs_Num[hh]-Num_Offs_Num[hh]%%Num_Offs[hh]
      }
    }
    if(sum(Num_Offs_Num)!=Pop_Siz_Gen_pop[i+1]) Num_Offs_Num[1]=Num_Offs_Num[1]-(sum(Num_Offs_Num)-Pop_Siz_Gen_pop[i+1])
    # samples mother and father
    ind=data.frame(DamID=rep(NA,Pop_Siz_Gen_pop[i+1]),SireID=rep(NA,Pop_Siz_Gen_pop[i+1]))
    for (ij in 1:length(Num_Offs_Num)) {
      # when ij is 1
      if(ij==1){
        # select for Dam
        if(length(selected_Dam)<Num_Offs_Num[ij]){
          ind$DamID[1:Num_Offs_Num[ij]]=c(selected_Dam,sample(selected_Dam,Num_Offs_Num[ij]-length(selected_Dam),replace =T))
        }else{
          ind$DamID[1:Num_Offs_Num[ij]]=sample(selected_Dam,Num_Offs_Num[ij],replace =F)
        }
        # select for sire
        if(length(selected_Sire)<Num_Offs_Num[ij]){
          ind$SireID[1:Num_Offs_Num[ij]]=c(selected_Sire,sample(selected_Sire,Num_Offs_Num[ij]-length(selected_Sire),replace =T))
        }else{
          ind$SireID[1:Num_Offs_Num[ij]]=sample(selected_Sire,Num_Offs_Num[ij],replace =F)
        }
        # when ij is more than 1
      }else{
        ind$DamID[(Num_Offs_Num[ij-1]+1):sum(Num_Offs_Num[1:ij])]=rep(sample(selected_Dam,Num_Offs_Num[ij]/Num_Offs[ij]),each=Num_Offs[ij])
        ind$SireID[(Num_Offs_Num[ij-1]+1):sum(Num_Offs_Num[1:ij])]=rep(sample(selected_Sire,Num_Offs_Num[ij]/Num_Offs[ij]),each=Num_Offs[ij])
      }
    }
    ped2=data.frame(ID=(max(ped$ID)+1):(Pop_Siz_Gen_pop[i+1]+max(ped$ID)),
                    DamID=ind$DamID,
                    SireID=ind$SireID,
                    Sex=if(Pop_Siz_Gen_pop[i+1]%%2==0) {
                      rep(c("M","F"),each=Pop_Siz_Gen_pop[i+1]/2)} else {
                        c(rep("M",each=floor(Pop_Siz_Gen_pop[i+1]/2)),rep("F",each=ceiling(Pop_Siz_Gen_pop[i+1]/2)))},
                    Gen=i)
    ped=rbind(ped,ped2)
    ###---------------------------
    ### simulating the STR
    IDss=ped$ID[ped$Gen==i]
    registerDoParallel(n.cores)
    # Start simulation
    STR2=foreach (str_sim=1:length(IDss), .combine=rbind) %dopar% {
      STR_Dam=STR[rownames(STR)==as.character(ped$DamID[ped$ID==IDss[str_sim]]),]
      STR_Sire=STR[rownames(STR)==as.character(ped$SireID[ped$ID==IDss[str_sim]]),]
      Temp=paste0(sapply(STR_Sire,function(v) intToUtf8(utf8ToInt(v)[sort(sample(nchar(v), 1))]),USE.NAMES = FALSE),
                  sapply(STR_Dam,function(v) intToUtf8(utf8ToInt(v)[sort(sample(nchar(v), 1))]),USE.NAMES = FALSE))
      
      # Alter genotype randomly for offs based on the mutation rate
      # offs
      Error_allel1=which(runif(length(Temp))<(Mutation_Rate/2))
      Error_allel2=which(runif(length(Temp))<(Mutation_Rate/2))
      if(length(Error_allel1)>0){
        for (mut in Error_allel1) {
          alleles=letters[1:Allele_num[mut]][!letters[1:Allele_num[mut]]%in%unlist(strsplit(as.character(c(STR_Dam[mut],STR_Sire[mut])),""))]
          if(length(alleles)==0)alleles=letters[1:Allele_num[mut]][!letters[1:Allele_num[mut]]%in%unlist(strsplit(as.character(c(STR_Sire[mut])),""))]
          if(length(alleles)==0)alleles=letters[1:Allele_num[mut]]
          
          Temp[mut]=paste0(sample(alleles,1),
                               unlist(strsplit(as.character(Temp[mut]),""))[2])
        }
      }
      if(length(Error_allel2)>0){
        for (mut in Error_allel2) {
          alleles=letters[1:Allele_num[mut]][!letters[1:Allele_num[mut]]%in%unlist(strsplit(as.character(c(STR_Dam[mut],STR_Sire[mut])),""))]
          if(length(alleles)==0)alleles=letters[1:Allele_num[mut]][!letters[1:Allele_num[mut]]%in%unlist(strsplit(as.character(c(STR_Dam[mut])),""))]
          if(length(alleles)==0)alleles=letters[1:Allele_num[mut]]
          Temp[mut]=paste0(unlist(strsplit(as.character(Temp[mut]),""))[1],
                           sample(alleles,1))
        }
      }
      
      return(Temp)
    }
    stopImplicitCluster()
    rownames(STR2)=IDss
    STR2=as.data.frame(STR2)
    STR=rbind(STR,STR2)
    cat(paste0("Recent generation ",i," is done! \n"))
  }
  ###---------------------------------------------------------------
  ### Alter genotype randomly for the loci based on the typing error rate
  ###---------------------------------------------------------------
  if(!is.null(Typing_Error)){
    for (i in 1:ncol(STR)) {
      Error_allel1=sample(ped$ID[ped$Gen!=0],round(length(ped$ID[ped$Gen!=0])*(Typing_Error[i])))
      if(length(Error_allel1)>0){
        for (j in 1:length(Error_allel1)) {
          STR_Dam=STR[rownames(STR)==as.character(ped$DamID[ped$ID==Error_allel1[j]]),i]
          STR_Sire=STR[rownames(STR)==as.character(ped$SireID[ped$ID==Error_allel1[j]]),i]
          alleles=letters[1:Allele_num[i]][!letters[1:Allele_num[i]]%in%unlist(strsplit(as.character(c(STR_Dam,STR_Sire)),""))]
          if(length(alleles)==0)alleles=letters[1:Allele_num[i]][!letters[1:Allele_num[i]]%in%unlist(strsplit(as.character(c(STR_Sire)),""))]
          if(length(alleles)==0)alleles=letters[1:Allele_num[i]]
          STR[Error_allel1[j],i]=paste0(sample(alleles,1),
                           unlist(strsplit(as.character(STR[Error_allel1[j],i]),""))[2])
        }
      }
      Error_allel2=sample(ped$ID[ped$Gen!=0],round(length(ped$ID[ped$Gen!=0])*(Typing_Error[i])))
      if(length(Error_allel2)>0){
        for (j in 1:length(Error_allel2)) {
          STR_Dam=STR[rownames(STR)==as.character(ped$DamID[ped$ID==Error_allel2[j]]),i]
          STR_Sire=STR[rownames(STR)==as.character(ped$SireID[ped$ID==Error_allel2[j]]),i]
          alleles=letters[1:Allele_num[i]][!letters[1:Allele_num[i]]%in%unlist(strsplit(as.character(c(STR_Dam,STR_Sire)),""))]
          if(length(alleles)==0)alleles=letters[1:Allele_num[i]][!letters[1:Allele_num[i]]%in%unlist(strsplit(as.character(c(STR_Dam)),""))]
          if(length(alleles)==0)alleles=letters[1:Allele_num[i]]
          STR[Error_allel2[j],i]=paste0(unlist(strsplit(as.character(STR[Error_allel2[j],i]),""))[1],
                                        sample(alleles,1))
        }
      }
    }
  }
  ped[ped$Gen==0,c("DamID","SireID")]=NA
  Ped=ped
  All_result=list()
  All_result[["Simulated_Ped"]]=Ped
  All_result[["Simulated_STR"]]=STR
  return(All_result)
}



