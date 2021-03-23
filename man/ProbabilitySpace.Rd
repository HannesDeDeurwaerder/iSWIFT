\name{ProbabilitySpace}
\alias{ProbabilitySpace}
\title{Creating synthetic field data by running the SWIFT model}
\description{ 
In this function, all required parameters of the SWIFT model are declared to 
perform a user-defined number of SWIFT runs, creating synthetic xylem water
isotope composition (D2H and D18O). 

NOTE that in this specific configuration, the function apeals to the 'VarMatrix'
and the 'SoilHeterogeneity' functions in order to retreive all SWIFT variables.
The user is encouraged to alter both functions according to the requirements of 
their study.

IMPORTANT: Note that for this function to properly work, all variables of the 
SWIFT model should be declared locally or globally. 
}
\usage{ ProbabilitySpace( beta.i=NULL, 
                          BiophysVariables=NULL,
                          PSIprofiles=NULL,
                          soilDeuterium=NULL, 
                          soilOxygen=NULL)}
\arguments{
  \item{beta.i}{
 Description: Allocation parameter that describes the exponential decay of the 
              distribution of absorbing root surfaces for the focal plant;
 Functionality: used in the iSWIFT loglikelihood optimization  
 Structure:   Randomly assigned numeric value within defined ranges as provided 
              in the 'optim' function.}
  \item{BiophysVariables}{
 Description: Matrix of plant-trait  and sampling strategy required variables 
              required to declare all SWIFT model variables;
 Structure:   Matrix of SWIFT input variables and/or variables required to derive 
              the SWIFT input variables (i.e. DBH, lumen fraction,...)
 Creation:    Via 'VarMatrix' function}              
  \item{PSIprofiles}{
 Description: Soil water potential values with soil depth over the discritized 
              soil layers
 Structure:   Matrix of distinct soil water potetnial profiles with depth. 
              Each column presents a soil water potential profile per planned 
              SWIFT run, rows correspond to the discritized soil layers
              each individual SWIFT run planned;
 Creation:    Via 'SoilHeterogeneity' function}              
  \item{soilDeuterium}{
 Description: Soil water deuterium composition profiles with soil depth over the 
              discritized soil layers
 Structure:   Matrix of distinct soil water deuterium composition with depth. 
              Each column presents a soil water deuterium composition profile 
              per planned SWIFT run, rows correspond to the discritized soil 
              layers each individual SWIFT run planned;
 Creation:    Via 'SoilHeterogeneity' function}  
  \item{soilOxygen}{
 Description: Soil water 18 Oxygen composition profiles with soil depth over the 
              discritized soil layers
 Structure:   Matrix of distinct soil water 18 Oxygen composition with depth. 
              Each column presents a soil water 18 Oxygen composition profile 
              per planned SWIFT run, rows correspond to the discritized soil 
              layers each individual SWIFT run planned;
 Creation:    Via 'SoilHeterogeneity' function}  

\details{
}
\value{
Returning values are SWIFT generated xylem water istopic composition, i.e., 
D2H and/or D18O. Hence providing synthetic field data corresponding the input
provided by the reader, but without consideration of potential extraction error}
\references{}
\author{ Hannes De Deurwaerder, Marco D. Visser, Matteo Detto, Pascal Boeckx, 
Felicien Meunier, Pedro Herve Fernandez and Hans Verbeeck}
\note{}
\seealso{}
\examples{
ProbabilitySpace <- function(beta.i, BiophysVariables, PSIprofiles, 
                             soilDeuterium, soilOxygen){
    #--------------------------------------------------------------------#
    # NOTE that for this function to properly work, all variables of the #
    # SWIFT model should be declared locally or globally.                #
    #--------------------------------------------------------------------#

  Iso_H <- rep(NA, SWIFTitterations )   # decalre empty vector of isotopic data

  # which isotope is considered
  if(ConsideredIsotopes=='D2H'){Isotopes <- soilDeuterium}  # deuterium
  if(ConsideredIsotopes!='D2H'){Isotopes <- soilOxygen} 
                                  # oxygen or dual using the LMWL approach
  for (ii in 1:SWIFTitterations){
    Ax <- BiophysVariables[ii,'LA']*BiophysVariables[ii,'Asapwood'] # lumen area
    SF <- relSF*BiophysVariables[ii,'SFdaily']* 3600/tF 
                                      # Sapflow in m3/time-frequency considered
    k  <- SoilRootCond(betaCom, BiophysVariables[ii,'Kr'], PSIprofiles[,ii], 
                       Z, 'Sandy Clay') # Total soil to root conductivity
    tmp <- beta.i^(100*Z)*(1-beta.i^(100*dZ))
    ARi <- BiophysVariables[ii,'ARtot']*tmp/sum(tmp) # root area distribution

    # run SWIFT for stem base and sampling at specific height    
    Iso_SB <- SWIFT_SB(ARi, Isotopes[,ii], k, PSIprofiles[,ii], SF, tt, Z)
    Iso_H[ii] <- SWIFT_H( Ax, Iso_SB, BiophysVariables[ii,'hos'], SF, 
                          BiophysVariables[ii,'tstud'], tF)
    # report progress of one run
    cat("\r","beta = ",beta.i, " | progress = ", 
        format((ii/SWIFTitterations)*100,nsmall=2),"%")
  }
  return(Iso_H)
}
}