#----------------------------------------------------------------------------------
# Calculating of LOD for one parents when the other is unknown
#----------------------------------------------------------------------------------
LOD_Both_Unknown <- function(Offs_STR,
                             Parent_STR,
                             Allele_freq,
                             Typing_Error){
  Allele_Sim=stringdist::stringsim(toupper(Parent_STR), toupper(Offs_STR),method ="qgram")*2
  if(length(which(!is.na(Allele_Sim)))>0){
    Likelihood_all=NULL
    for (i in which(!is.na(Allele_Sim))) {
      Offs_alleles=unlist(strsplit(as.character(Offs_STR[i]),""))
      Parent_alleles=unlist(strsplit(as.character(Parent_STR[i]),""))
      # calculating P(g_o)
      if(Offs_alleles[1]==Offs_alleles[2]){
        P_Offs=(Allele_freq[[names(Parent_STR)[i]]][[1]][which(names(Allele_freq[[i]][[1]])%in%Offs_alleles[1])])^2
      }else{
        P_Offs=(Allele_freq[[names(Parent_STR)[i]]][[1]][which(names(Allele_freq[[i]][[1]])%in%Offs_alleles[1])]*
                  Allele_freq[[names(Parent_STR)[i]]][[1]][which(names(Allele_freq[[i]][[1]])%in%Offs_alleles[2])])*2
      }
      # calculating P(g_m)
      if(Parent_alleles[1]==Parent_alleles[2]){
        P_parent=(Allele_freq[[names(Parent_STR)[i]]][[1]][which(names(Allele_freq[[i]][[1]])%in%Parent_alleles[1])])^2
      }else{
        P_parent=(Allele_freq[[names(Parent_STR)[i]]][[1]][which(names(Allele_freq[[i]][[1]])%in%Parent_alleles[1])]*
                    Allele_freq[[names(Parent_STR)[i]]][[1]][which(names(Allele_freq[[i]][[1]])%in%Parent_alleles[2])])*2
      }
      # calculating T(g_o|g_m)
      if(Allele_Sim[i]==0){
        T_parent=0
      }else{
        if(Allele_Sim[i]==2){# if two alleles are similar
          T_parent=(Allele_freq[[names(Parent_STR)[i]]][[1]][which(names(Allele_freq[[i]][[1]])%in%Parent_alleles[1])]+
                      Allele_freq[[names(Parent_STR)[i]]][[1]][which(names(Allele_freq[[i]][[1]])%in%Parent_alleles[2])])/2
        }else{# if one alleles are similar
          if(Offs_alleles[1]==Offs_alleles[2]){
            T_parent=(Allele_freq[[names(Parent_STR)[i]]][[1]][which(names(Allele_freq[[i]][[1]])%in%Offs_alleles[1])])/2
          }else{
            if(Parent_alleles[1]==Parent_alleles[2]){
              T_parent=(Allele_freq[[names(Parent_STR)[i]]][[1]][which(names(Allele_freq[[i]][[1]])%in%Offs_alleles[!Offs_alleles%in%Parent_alleles])])
            }else{
              T_parent=(Allele_freq[[names(Parent_STR)[i]]][[1]][which(names(Allele_freq[[i]][[1]])%in%Offs_alleles[!Offs_alleles%in%Parent_alleles])])/2
            }
          }
        }
      }
      # Calculating L(H1), L(H2) and L(H1,H2) or LOD
      Likelihood_all=c(Likelihood_all,as.numeric((P_parent*(
        (((1-Typing_Error[i])^2)*T_parent)+
          (Typing_Error[i]*((1-Typing_Error[i])^2)*P_Offs)+
          ((Typing_Error[i]^2)*P_Offs)
      ))/(P_parent*(
        (((1-Typing_Error[i])^2)*P_Offs)+
          (Typing_Error[i]*((1-Typing_Error[i])^2)*P_Offs)+
          ((Typing_Error[i]^2)*P_Offs)
      ))))
    }
    Temp=round(log(prod(Likelihood_all)),3)
  }else Temp=NA
  return(Temp)
}
