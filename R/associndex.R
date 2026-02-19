#' Calculate different association indices considering different interaction types.
#'
#' @description Calculates different metrics to characterize the interaction strength
#' for each pairwise interaction, based on the number of recruits found under
#' the canopy species and in the "Open" ground. Several indices have been proposed
#' to assess how much recruitment is affected by the canopy species relative
#' to the species performance in the absence of canopy plants (i.e. in open),
#' each of them with its own particularities. In the context of recruitment
#' interactions, the most frequently used index is the *Relative Interaction Index*
#' (**RII**; Armas et al., 2004), and some studies have recently used the
#' *additive and commutative symmetry intensity indices* (**NIntA** and **NIntC**)
#' y Diaz-Sierra et al. (2017). However, these indices are not well suited for
#' comparisons of interaction strength between pairs of species within a local
#' community, so the *Normalized Neighbour Suitability* index (**Ns**; Mingo, 2014)
#'  should be preferred (Alcantara et al. 2025). The function *`associndex`*
#'  calculates all the four indices besides providing the density of recruits
#'  under the canopy and in the open, and the data used to calculate it
#'  (relative cover and number of recruits in each microhabitat).
#'  **IMPORTANT NOTE: Data on recruitment in Open is required to calculate
#'  these indices.
#'  **Input**: The canopy-recruit interactions data set and the canopy cover
#'  data set
#'
#' @inheritParams check_interactions
#' @inheritParams check_cover
#'
#' @param expand Expands the data frame to include the non-detected interactions
#' between every pair of species that occur in the study site, assigning values
#' of 0 to the variables Fcr and Dcr. Explanation of its options:
#' - *yes*: Generates all possible interactions between canopy and recruit species
#' in the study site (i.e. detected and non-detected), adding values of interaction
#' frequency of 0 to the non-detected interactions.
#' - *no*: Considers only the pairwise interactions actually detected in the study site.
#'
#' @param rm_sp_no_cover **rm_sp_no_cover**: Defines the filtering behavior, specifying
#' which species should be removed based on the availability of cover data and the type
#' of interaction of interest. Explanation of its options:
#' - *allsp*: removes any species that lacks cover data (both canopy and recruit species).
#' - *onlycanopy*: removes only canopy species lacking cover data, retaining recruit
#' species even if they lack cover data.
#'
#' @param threshold_density **threshold_density**: indicates a threshold to retain only
#' those interactions in which the observed density is below the threshold
#' (i.e. Dcr < thr & Dro < thr). This is required to avoid the use of interactions
#' based on very few recruits which may have occurred by chance under a very scarce
#' canopy species. These cases will result in a large overestimation of Dcr, likely
#' being spurious outlayers. If no threshold value is provided or it is set to NULL,
#' the function estimates a threshold value based on recruitment density outliers
#' according to a Weibull distribution: a value is an outlier if it is above the limit
#' where less than 1 observation is expected. Explanation of its options:
#' - A positive value. High values provide more chance to include potentially spurious
#' interactions, unless such a high density is reliable based on the knowledge of the
#' study system. A way to decide the threshold value would consist in identifying possible
#' outlayers in the distribution of recruitment densities (*Dcr* and *Dro*).
#' - NULL: use package **extremevalues** to find the density threshold value (*K*)
#' above which, according to the Weibull distribution, less than 1 observation is
#' expected given the number of observations in the sample.
#'
#' @returns a data frame with the following information for each canopy-recruit
#' interaction with the following information:
#' - *Canopy*: **canopy species**.
#' - *Recruit*: **recruit species**.
#' - *Fcr*: **frequency of recruitment**, as the number of recruits found
#' under that canopy species.
#' - *Ac*: **Area occupied by the canopy species** in the study site,
#' in m^2^ (or m when relative cover is measured with transects).
#' - *Fro*: **frequency of recruitment**, as the number of recruits in
#' Open interspaces. - *Ao*: **Area  occupied by open interspaces** in the study
#' site in m^2^ (or m when relative cover is measured with transects).
#' - *Dcr*: **Density of recruits under a canopy species**: $Fcr/Ac$.
#' - *Dro*: **Density of recruits in open interspaces**: $Fro/Ao$.
#' - *max_Recr*: **maximum recruitment density** of a species in the study site
#' (i.e., under any canopy species or "Open").
#' - *Ns*: **Normalized Neighbour Suitability** index: $(Dcr - Dro) / maxRecr$.
#' - *NintC*: **Commutative symmetry intensity** index: $2*(Dcr - Dro)/((Dcr + Dro)+abs(Dcr-Dro))$
#' - *NintA*: **Additive symmetry intensity** index: $2*(Dcr - Dro)/((Dro) + abs(Dcr-Dro))$
#' - *RII*: **Relative Interaction Index**: $(Dcr - Dro)/(Dcr + Dro)$
#'
#' @export
#'
#' @examples associndex (Amoladeras_int, Amoladeras_cover, expand = "yes",
#' rm_sp_no_cover = "allsp", thresold_density=NULL)
#' associndex (Amoladeras_int, Amoladeras_cover, expand = "yes",
#' rm_sp_no_cover = "onlycanopy", thresold_density=NULL)
#' associndex (Amoladeras_int, Amoladeras_cover, expand = "no",
#' rm_sp_no_cover = "allsp", thresold_density=NULL)
#' associndex (Amoladeras_int, Amoladeras_cover, expand = "yes",
#' rm_sp_no_cover = "onlycanopy", thresold_density=NULL)

associndex<- function(int_data, cover_data, expand=c("yes","no"),
                      rm_sp_no_cover=c("allsp","onlycanopy"),
                      threshold_density=NULL)

{

  if (!"Open" %in% int_data$Canopy) stop(
    "ERROR: your data does not contain a node named Open or it is spelled
    differently. Data for recruitment in Open is required to calculate the indices."
    )

  # Suggest a threshold value based on recruitment density outliers according
  # to a Weibull distribution. Uses "extremevalues" package.
  int_data_0<-comm_to_RN_UNI(int_data, cover_data)
  db_inter_0 <- pre_associndex_UNISITE_UNI(int_data_0)
  db_inter_0$Dcr <- db_inter_0$Fcr/db_inter_0$Ac
  db_inter_0$Dro <- db_inter_0$Fro/db_inter_0$Ao
  y <- rbind(db_inter_0$Dcr,db_inter_0$Dro)
  y <- y[which(y>0)] # Zero densities are common, so we will check only for
  #suspiciously high density values.
  K <- extremevalues::getOutliers(y, method = "I", distribution = "weibull", rho=c(1,1)) # We
  #consider that a value is an outlier if it is above the limit where less
  #than 1 observation is expected.

  if(is.null(threshold_density)){
    threshold_density = K$limit[[2]]
    message("Based on Weibull distribution fitted to the observed values, the threshold density has been set to: ", threshold_density, ".")
  }


  if(expand=="yes" & rm_sp_no_cover=="allsp"){int_type ="rec"}
  if(expand=="yes" & rm_sp_no_cover=="onlycanopy"){int_type ="comp"}
  if(expand=="no" & rm_sp_no_cover=="onlycanopy"){int_type ="fac"}
  if(expand=="no" & rm_sp_no_cover=="allsp"){int_type ="unused"}

  if(int_type=="rec"){

    int_data<-comm_to_RN_UNI(int_data, cover_data)
    db_inter<-associndex_UNISITE_UNI(int_data, threshold_density)

  }

  if(int_type=="fac"){

    int_data<-comm_to_RN_BI(int_data, cover_data)
    db_inter<-associndex_UNISITE_BI(int_data, threshold_density)

  }

  if(int_type=="comp"){

    int_data<-comm_to_RN_UNI_COMP(int_data, cover_data)
    db_inter<-associndex_UNISITE_BI_COMP(int_data, threshold_density)

  }

  if(int_type=="unused"){

    db_inter<-NULL
    warning("The object is not provided because the combination of arguments
            is not commonly used in plant-plant interaction networks")
  }

  return(db_inter)
}

