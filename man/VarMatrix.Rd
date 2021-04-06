\name{VarMatrix}
\alias{VarMatrix}
\title{ Creating a plant-traits variables matrix for SWIFT synthetic field data
generation, given user defined variable ranges}
\description{ 
The 'VarMatrix' function creates a matrix of SWIFT input variables related to 
plant-traits and sampling strategy. Considering user defined variable ranges 
(based on prior knowledge of the studied plants), values will be randomly extracted 
from these ranges and allocated to specific variable vector. The function returns 
a matrix of several SWIFT input varialbes or values related to obtain those 
(i.e., here defined for "DBH","Area_alg", "SFtot_alg","ARtot","LA", "Asapwood",
"Kr","SFdaily",'beta','tstud' and 'hos') where columns represent the 
individual variables and the rows the user defined number of syntehtic field 
sampples to generate.

IMPORTANT: Note that for this function to properly work, all variables of the 
SWIFT model should be declared locally or globally.} 
\usage{VarMatrix(beta.i=NULL)}
\arguments{
  \item{beta.i}{
 Description: Allocation parameter that describes the exponential decay of the 
              distribution of absorbing root surfaces for the focal plant;
 Functionality: used in the iSWIFT loglikelihood optimization  
 Structure:   Randomly assigned numeric value within defined ranges as provided 
              in the 'optim' function.}
}
\details{}
\value{
  The function returns 
  a matrix of several SWIFT input varialbes with a column for each considered  
  variables and the rows the user defined number of syntehtic field 
  sampples to generate.
}
\references{}
\author{ Hannes De Deurwaerder, Marco D. Visser, Matteo Detto, Pascal Boeckx, 
Felicien Meunier, Pedro Herve Fernandez and Hans Verbeeck}
\note{}
\seealso{}
\examples{
VarMatrix <- function(beta.i=NULL){
    #--------------------------------------------------------------------#
    # NOTE that for this function to properly work, all variables/ranges #
    # of the SWIFT model should be declared locally or globally.         #
    #--------------------------------------------------------------------#
    
  # create empty data matrix
  SWIFTvars=matrix(NA,SWIFTitterations,ncol=11)
  colnames(SWIFTvars)<-c("DBH","Area_alg", "SFtot_alg","ARtot",
                          "LA", "Asapwood","Kr","SFdaily",'beta','tstud','hos') 
  
  # DBH [in m]
    SWIFTvars[,'DBH']<-runif(SWIFTitterations, 
                           min = DBHrange[1], max = DBHrange[2])
  # Stem cross sectional area [in m2]
    SWIFTvars[,'Area_alg'] <- pi*(SWIFTvars[,'DBH']/2)^2   
  # Total absorbing root area [in m2] (growth form specific function) 
    SWIFTvars[,'ARtot'] <- ARcalc(SWIFTvars[,'DBH']) 
  # Lumen fraction [-]: 
    SWIFTvars[,'LA'] <- runif(SWIFTitterations, 
                              min = LArange[1], max = LArange[2])
  # Sapwood area [in m2] (logaritm of Meinzer et al(2001)):
    Asap <- 1.582*((100*SWIFTvars[,'DBH'])^1.764) *10^-4  
    SWIFTvars[,'Asapwood'] <- Asap*runif(SWIFTitterations, 
                                         min = Asaprange[1], max = Asaprange[2])
        # Correction: sapwood is never bigger then tree area
        if(any(SWIFTvars[,'Area_alg'] <= SWIFTvars[,'Asapwood'])){
          ind = which(SWIFTvars[,'Area_alg'] <= SWIFTvars[,'Asapwood'])
          SWIFTvars[ind,'Asapwood'] <- SWIFTvars[ind,'Area_alg']}
 
  # Daily total Sapflow volume derived from DBH and growth form specific rates
    SFday=mean(SapFlow) * (10^-3) * (3600*24)   # m3 of water per day
    SWIFTvars[,'SFtot_alg']=SFday*(SWIFTvars[,'Asapwood']) 
  # Total daily sapflow with some variation added  
    SWIFTvars[,'SFdaily'] <- SWIFTvars[,'SFtot_alg'] * 
          runif(SWIFTitterations, min = SFdailyrange[1], max =  SFdailyrange[2])
  # Kr, the effective root radial conductivity [in s-1]
    SWIFTvars[,'Kr'] <- runif(SWIFTitterations, 
                              min = Krrange[1], max = Krrange[2])
  # Beta
    SWIFTvars[,'beta'] <- beta.i          
  # Sampling time
    SWIFTvars[,'tstud'] <- round(runif(SWIFTitterations,
                  min=tstudrange[1]+TCOR, max=tstudrange[2]+TCOR),0)
                  # a time correction is add to account for the model spin-up
                  # in this specific example
  #Height of sampling
    SWIFTvars[,'hos'] <- runif(SWIFTitterations, min=Hosrange[1], max=Hosrange[2])  

  return(SWIFTvars)
}
}
