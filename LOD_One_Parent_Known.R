#----------------------------------------------------------------------------------
# Calculating of LOD when one parent is known
#----------------------------------------------------------------------------------
LOD_One_Parent_Known <- function(Offs_STR,
                                    Alleged_Parent_STR,
                                    Known_Parent_STR,
                                    Allele_freq,
                                    Typing_Error){
  Allele_Sim=stringdist::stringsim(toupper(Alleged_Parent_STR), toupper(Offs_STR),method ="qgram")*2
  Allele_Sim_Mo=stringdist::stringsim(toupper(Known_Parent_STR), toupper(Offs_STR),method ="qgram")*2
  if(length(which(!is.na(Allele_Sim)))>0){
    Likelihood_all=NULL
    for (i in which(!is.na(Allele_Sim))) {
      Offs_alleles=unlist(strsplit(as.character(Offs_STR[i]),""))
      Alleged_alleles=unlist(strsplit(as.character(Alleged_Parent_STR[i]),""))
      Known_alleles=unlist(strsplit(as.character(Known_Parent_STR[i]),""))
      # calculating P(g_o)
      if(Offs_alleles[1]==Offs_alleles[2]){
        P_Offs=(Allele_freq[[names(Alleged_Parent_STR)[i]]][[1]][which(names(Allele_freq[[i]][[1]])%in%Offs_alleles[1])])^2
      }else{
        P_Offs=(Allele_freq[[names(Alleged_Parent_STR)[i]]][[1]][which(names(Allele_freq[[i]][[1]])%in%Offs_alleles[1])]*
                  Allele_freq[[names(Alleged_Parent_STR)[i]]][[1]][which(names(Allele_freq[[i]][[1]])%in%Offs_alleles[2])])*2
      }
      # calculating P(g_a)
      if(Alleged_alleles[1]==Alleged_alleles[2]){
        P_alleles=(Allele_freq[[names(Alleged_Parent_STR)[i]]][[1]][which(names(Allele_freq[[i]][[1]])%in%Alleged_alleles[1])])^2
      }else{
        P_alleles=(Allele_freq[[names(Alleged_Parent_STR)[i]]][[1]][which(names(Allele_freq[[i]][[1]])%in%Alleged_alleles[1])]*
                    Allele_freq[[names(Alleged_Parent_STR)[i]]][[1]][which(names(Allele_freq[[i]][[1]])%in%Alleged_alleles[2])])*2
      }
      # calculating T(g_o|g_a)
      if(Allele_Sim[i]==0){
        T_alleles=0
      }else{
        if(Allele_Sim[i]==2){# if two alleles are similar
          T_alleles=(Allele_freq[[names(Alleged_Parent_STR)[i]]][[1]][which(names(Allele_freq[[i]][[1]])%in%Alleged_alleles[1])]+
                      Allele_freq[[names(Alleged_Parent_STR)[i]]][[1]][which(names(Allele_freq[[i]][[1]])%in%Alleged_alleles[2])])/2
        }else{# if one alleles are similar
          if(Offs_alleles[1]==Offs_alleles[2]){
            T_alleles=(Allele_freq[[names(Alleged_Parent_STR)[i]]][[1]][which(names(Allele_freq[[i]][[1]])%in%Offs_alleles[1])])/2
          }else{
            if(Alleged_alleles[1]==Alleged_alleles[2]){
              T_alleles=(Allele_freq[[names(Alleged_Parent_STR)[i]]][[1]][which(names(Allele_freq[[i]][[1]])%in%Offs_alleles[!Offs_alleles%in%Alleged_alleles])])
            }else{
              T_alleles=(Allele_freq[[names(Alleged_Parent_STR)[i]]][[1]][which(names(Allele_freq[[i]][[1]])%in%Offs_alleles[!Offs_alleles%in%Alleged_alleles])])/2
            }
          }
        }
      }
      # when known genotype missed
      if(is.na(Known_Parent_STR[1,i])){
        # Calculating L(H1), L(H2) and L(H1,H2) or LOD
        Likelihood_all=c(Likelihood_all,as.numeric((P_alleles*(
          (((1-Typing_Error[i])^2)*T_alleles)+
            (Typing_Error[i]*((1-Typing_Error[i])^2)*P_Offs)+
            ((Typing_Error[i]^2)*P_Offs)
        ))/(P_alleles*(
          (((1-Typing_Error[i])^2)*P_Offs)+
            (Typing_Error[i]*((1-Typing_Error[i])^2)*P_Offs)+
            ((Typing_Error[i]^2)*P_Offs)
        ))))
      }else{
        # calculating P(g_m)
        if(Known_alleles[1]==Known_alleles[2]){
          P_known=(Allele_freq[[names(Alleged_Parent_STR)[i]]][[1]][which(names(Allele_freq[[i]][[1]])%in%Known_alleles[1])])^2
        }else{
          P_known=(Allele_freq[[names(Alleged_Parent_STR)[i]]][[1]][which(names(Allele_freq[[i]][[1]])%in%Known_alleles[1])]*
                      Allele_freq[[names(Alleged_Parent_STR)[i]]][[1]][which(names(Allele_freq[[i]][[1]])%in%Known_alleles[2])])*2
        }
        # calculating T(g_o|g_m,g_a)
        Possible_geno=paste0(as.data.frame(expand.grid(Known_alleles,Alleged_alleles))[,1],as.data.frame(expand.grid(Known_alleles,Alleged_alleles))[,2])
        T_Offs_STR_known_both=sum((stringdist::stringsim(paste0(toupper(Offs_alleles),collapse = ""), toupper(Possible_geno),method ="qgram")*2)==2)/4
        # calculating T(g_o|g_m)
        if(Allele_Sim_Mo[i]==0){
          T_known=0
        }else{
          if(Allele_Sim_Mo[i]==2){# if two alleles are similar
            T_known=(Allele_freq[[names(Known_Parent_STR)[i]]][[1]][which(names(Allele_freq[[i]][[1]])%in%Known_alleles[1])]+
                        Allele_freq[[names(Known_Parent_STR)[i]]][[1]][which(names(Allele_freq[[i]][[1]])%in%Known_alleles[2])])/2
          }else{# if one alleles are similar
            if(Offs_alleles[1]==Offs_alleles[2]){
              T_known=(Allele_freq[[names(Known_Parent_STR)[i]]][[1]][which(names(Allele_freq[[i]][[1]])%in%Offs_alleles[1])])/2
            }else{
              if(Known_alleles[1]==Known_alleles[2]){
                T_known=(Allele_freq[[names(Known_Parent_STR)[i]]][[1]][which(names(Allele_freq[[i]][[1]])%in%Offs_alleles[!Offs_alleles%in%Known_alleles])])
              }else{
                T_known=(Allele_freq[[names(Known_Parent_STR)[i]]][[1]][which(names(Allele_freq[[i]][[1]])%in%Offs_alleles[!Offs_alleles%in%Known_alleles])])/2
              }
            }
          }
        }
        # Calculating L(H1), L(H2) and L(H1,H2) or LOD
        Likelihood_all=c(Likelihood_all,as.numeric((P_known*P_alleles*(
          (((1-Typing_Error[i])^3)*T_Offs_STR_known_both)+
            (Typing_Error[i]*((1-Typing_Error[i])^2)*(T_known+T_alleles+P_Offs))+
            ((Typing_Error[i]^2)*((1-Typing_Error[i])^3)*P_Offs)+((Typing_Error[i]^3)*P_Offs)))/(P_known*P_alleles*(
              (((1-Typing_Error[i])^3)*T_known)+
                (Typing_Error[i]*((1-Typing_Error[i])^2)*(T_known+P_Offs))+
                ((Typing_Error[i]^2)*((1-Typing_Error[i])^3)*P_Offs)+((Typing_Error[i]^3)*P_Offs)))))
      }
    }
    Temp=round(log(prod(Likelihood_all)),3)
  }else Temp=NA
  return(Temp)
}
