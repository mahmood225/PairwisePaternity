#----------------------------------------------------------------------------------
# Calculating of LOD for both parents jointly
#----------------------------------------------------------------------------------
LOD_Parents_Jointly <- function(Offs_STR,
                                Mother_STR,
                                Father_STR,
                                Allele_freq,
                                Typing_Error){
  Allele_test=which(!is.na(Offs_STR) & !is.na(Mother_STR) & !is.na(Father_STR))
  Allele_Sim_Mo=stringdist::stringsim(toupper(Mother_STR), toupper(Offs_STR),method ="qgram")*2
  Allele_Sim_Fa=stringdist::stringsim(toupper(Father_STR), toupper(Offs_STR),method ="qgram")*2
  if(length(Allele_test)>0){
    Likelihood_all=NULL
    for (i in Allele_test) {
      Offs_alleles=unlist(strsplit(as.character(Offs_STR[i]),""))
      Father_alleles=unlist(strsplit(as.character(Father_STR[i]),""))
      Mother_alleles=unlist(strsplit(as.character(Mother_STR[i]),""))
      # calculating P(g_o)
      if(Offs_alleles[1]==Offs_alleles[2]){
        P_Offs=(Allele_freq[[names(Allele_freq)[i]]][[1]][which(names(Allele_freq[[i]][[1]])%in%Offs_alleles[1])])^2
      }else{
        P_Offs=(Allele_freq[[names(Allele_freq)[i]]][[1]][which(names(Allele_freq[[i]][[1]])%in%Offs_alleles[1])]*
                  Allele_freq[[names(Allele_freq)[i]]][[1]][which(names(Allele_freq[[i]][[1]])%in%Offs_alleles[2])])*2
      }
      # calculating P(g_m)
      if(Mother_alleles[1]==Mother_alleles[2]){
        P_mother=(Allele_freq[[names(Mother_STR)[i]]][[1]][which(names(Allele_freq[[i]][[1]])%in%Mother_alleles[1])])^2
      }else{
        P_mother=(Allele_freq[[names(Mother_STR)[i]]][[1]][which(names(Allele_freq[[i]][[1]])%in%Mother_alleles[1])]*
                    Allele_freq[[names(Mother_STR)[i]]][[1]][which(names(Allele_freq[[i]][[1]])%in%Mother_alleles[2])])*2
      }
      # calculating P(g_a)
      if(Father_alleles[1]==Father_alleles[2]){
        P_father=(Allele_freq[[names(Father_STR)[i]]][[1]][which(names(Allele_freq[[i]][[1]])%in%Father_alleles[1])])^2
      }else{
        P_father=(Allele_freq[[names(Father_STR)[i]]][[1]][which(names(Allele_freq[[i]][[1]])%in%Father_alleles[1])]*
                    Allele_freq[[names(Father_STR)[i]]][[1]][which(names(Allele_freq[[i]][[1]])%in%Father_alleles[2])])*2
      }
      # calculating T(g_o|g_m)
      if(Allele_Sim_Mo[i]==0){
        T_mother=0
      }else{
        if(Allele_Sim_Mo[i]==2){# if two alleles are similar
          T_mother=(Allele_freq[[names(Mother_STR)[i]]][[1]][which(names(Allele_freq[[i]][[1]])%in%Mother_alleles[1])]+
                      Allele_freq[[names(Mother_STR)[i]]][[1]][which(names(Allele_freq[[i]][[1]])%in%Mother_alleles[2])])/2
        }else{# if one alleles are similar
          if(Offs_alleles[1]==Offs_alleles[2]){
            T_mother=(Allele_freq[[names(Mother_STR)[i]]][[1]][which(names(Allele_freq[[i]][[1]])%in%Offs_alleles[1])])/2
          }else{
            if(Mother_alleles[1]==Mother_alleles[2]){
              T_mother=(Allele_freq[[names(Mother_STR)[i]]][[1]][which(names(Allele_freq[[i]][[1]])%in%Offs_alleles[!Offs_alleles%in%Mother_alleles])])
            }else{
              T_mother=(Allele_freq[[names(Mother_STR)[i]]][[1]][which(names(Allele_freq[[i]][[1]])%in%Offs_alleles[!Offs_alleles%in%Mother_alleles])])/2
            }
          }
        } 
      }
      # calculating T(g_o|g_a)
      if(Allele_Sim_Fa[i]==0){
        T_father=0
      }else{
        if(Allele_Sim_Fa[i]==2){# if two alleles are similar
          T_father=(Allele_freq[[names(Father_STR)[i]]][[1]][which(names(Allele_freq[[i]][[1]])%in%Father_alleles[1])]+
                      Allele_freq[[names(Father_STR)[i]]][[1]][which(names(Allele_freq[[i]][[1]])%in%Father_alleles[2])])/2
        }else{# if one alleles are similar
          if(Offs_alleles[1]==Offs_alleles[2]){
            T_father=(Allele_freq[[names(Father_STR)[i]]][[1]][which(names(Allele_freq[[i]][[1]])%in%Offs_alleles[1])])/2
          }else{
            if(Father_alleles[1]==Father_alleles[2]){
              T_father=(Allele_freq[[names(Father_STR)[i]]][[1]][which(names(Allele_freq[[i]][[1]])%in%Offs_alleles[!Offs_alleles%in%Father_alleles])])
            }else{
              T_father=(Allele_freq[[names(Father_STR)[i]]][[1]][which(names(Allele_freq[[i]][[1]])%in%Offs_alleles[!Offs_alleles%in%Father_alleles])])/2
            }
          }
        } 
      }
      # calculating T(g_o|g_m,g_a)
      Possible_geno=paste0(as.data.frame(expand.grid(Father_alleles,Mother_alleles))[,1],as.data.frame(expand.grid(Father_alleles,Mother_alleles))[,2])
      T_Offs_STR_ma_fa=sum((stringdist::stringsim(paste0(toupper(Offs_alleles),collapse = ""), toupper(Possible_geno),method ="qgram")*2)==2)/4
      # Calculating L(H1), L(H2) and L(H1,H2) or LOD
      Likelihood_all = c(Likelihood_all,as.numeric((P_mother*P_father*(
        (((1-Typing_Error[i])^3)*T_Offs_STR_ma_fa)+
          (Typing_Error[i]*((1-Typing_Error[i])^2)*(T_mother+T_father+P_Offs))+
          ((Typing_Error[i]^2)*((1-Typing_Error[i])^3)*P_Offs)+
          ((Typing_Error[i]^3)*P_Offs)
      ))/(
        P_mother*P_father*(
          (((1-Typing_Error[i])^3)*P_Offs)+
            ((Typing_Error[i]*((1-Typing_Error[i])^2))*P_Offs)+
            ((Typing_Error[i]^2)*((1-Typing_Error[i])^3)*P_Offs)+
            ((Typing_Error[i]^3)*P_Offs)
        )
      )))
    }
    Temp=round(log(prod(Likelihood_all)),3)
  }else Temp=NA
  return(Temp)
}
