\name{SoilHeterogeneity}
\alias{SoilHeterogeneity}
\title{ generate soil water potetnial and soil water isotopic profiles}
\description{ 
This function generates soil water potetnial and soil water isotopic profiles 
with values provided per considered soil layer. This is based on field data, 
with curves fit through randomly extracted values of soil water potential/isotope
compositions within the natural range of the observed soil profiles (here 
3x Standard deviation range from the mean obsereved values).

IMPORTANT: Note that for this function to properly work, the matrixes 
'SoilWaterPotential' and 'SoilWaterIsotopes' should be declared globally. 
These data structures provide matrixes of average and standard deviation of 
respectively soil water potentials with depth and for deuterium and oxygen with 
depth. The structure can be consulted in the rda-files add in the package 
structure
}
\usage{ SoilHeterogeneity(DataType=NULL)}
\arguments{
  \item{DataType}{
 Description: String which contains the datatype for which the construction of a
              soil profile is required; options are:'PSI' for soil water 
              potential profiles with detphs,'D2H' or 'D18O' when soil water
              isotope composition profiles are targetted}
\details{
}
\value{
Returning a matrix containing user defined number of soil water potential or isotope composition profiles, where every column represent a distinct profile, and the rows correspond to the discritized soil layers}
\references{}
\author{ Hannes De Deurwaerder, Marco D. Visser, Matteo Detto, Pascal Boeckx, 
Felicien Meunier, Pedro Herve Fernandez and Hans Verbeeck}
\note{}
\seealso{}
\examples{

SoilHeterogeneity <-function(DataType=NULL){
    #--------------------------------------------------------------------#
    # NOTE that for this function to properly work, all variables of     #
    # 'SoilWaterPotential' and 'SoilWaterIsotopes' objects should be     #
    # declared globally. see example provided                            #
    #--------------------------------------------------------------------#
    
  # declare empty matrix to fill up with new soil water potential profiles
  Profiles <- matrix(NA, length(Z), SWIFTitterations) 
  
  ### Situation A: for soil water potentials ### -------------------------------
  if(DataType=="PSI"){
      depths <- SoilWaterPotential['Depth'][,1] # considered monitoring depths
      SampledObs = matrix(NA,SWIFTitterations,nrow(SoilWaterPotential))
      colnames(SampledObs) <- depths
  
      # generate synthetic 'soil samples' given measured field data 
      for(ii in 1:length(depths)){
          SampledObs[,ii]<-rnorm(SWIFTitterations,
                                 SoilWaterPotential[ii, 'AvgPsi'], 
                                 abs((SoilWaterPotential[ii,'sdPsi'])))
      }
      
      # Curve fit new soil profiles using the synthetic soil samples 
      for(jj in 1:SWIFTitterations){
          # Curve fitting function  
          PsiProfileFun <- function(param){
              aa=param[1];  bb=param[2];  cc=param[3] # Derived via optimization
              PSIest=(((aa + bb *log(depths) - cc*depths^2)))
              div=sum((abs(SampledObs[jj,] - PSIest))^2, na.rm=TRUE)
              } 
        PSIoptim <- optim(c(200, 450, 250), PsiProfileFun) # optimize
        Profiles[,jj] <- ((((PSIoptim$par[1]+PSIoptim$par[2]*log(Z)-
                           PSIoptim$par[3]*Z^2))))*CTpsi 
                           # CTpsi to transfer data from MPa to m H2O
        # profile is considered unidirectional (note: SWIFT defines Z positive)
        ind = Z[which(Profiles[,jj] == max(Profiles[,jj]))]
        Profiles[which(Z >= ind[1]),jj] <- max(Profiles[,jj])
        }
      # Post correction - soil can't be more wet then saturated
      Profiles[Profiles >= PSIsat] <- PSIsat 
   }
 
  ### Situation B: for soil water isotope composition ### ----------------------
  if(DataType !="PSI"){
      depths <- SoilWaterIsotopes['Depth'][,1]
      SampledObs = matrix(NA,SWIFTitterations,nrow(SoilWaterIsotopes))
      colnames(SampledObs) <- depths
    
      # generate synthetic 'soil samples' given measured field data 
        ### for deuterium isotopes
        if(DataType == "D2H"){ 
          for (ii in 1:length(depths)){
              SampledObs[,ii]<-rnorm(SWIFTitterations, 
                                 SoilWaterIsotopes[ii,'AvgD2H'],
                               (SoilWaterIsotopes[ii, 'sdD2H']))
          }
          param0 <- c(-75, 0.15) # start value for curve fitting optimization
        }
        ### for oxygen isotopes
        if(DataType=="D18O"){ # for oxygen isotopes
          for (ii in 1:length(depths)){
              SampledObs[,ii]<-rnorm(SWIFTitterations, 
                                     SoilWaterIsotopes[ii,'AvgD18O'],
                                 (SoilWaterIsotopes[ii,'sdD18O']))
          }
          param0 <- c(-10, 0.15) # start value for curve fitting optimization
        }
    
      # Curve fit new soil profiles using the synthetic soil samples 
      for(jj in 1:SWIFTitterations){
          # Curve fitting function  
          isoProfileFun <- function(param){
              ll=param[1];  mm=param[2]  # Derived via optimization !!!
              ISOest=ll*(depths+dZ)^mm   # dZ is add to avoid infinity values
              div=sum(abs(SampledObs[jj,] - ISOest)^2, na.rm=TRUE)
          }
        ISOoptim <-optim(param0, isoProfileFun) # optimize
        Profiles[,jj] <- ISOoptim$par[1] * (Z+dZ)^ISOoptim$par[2]
        }
  }
  return(Profiles)
}

}