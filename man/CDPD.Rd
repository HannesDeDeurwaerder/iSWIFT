\name{CDPD}
\alias{CDPD}
\title{Generating conditional density probability distributions using SWIFT model}
\description{ 
This function is used to feed the required data to the 'SWIFT model' in order to
generate synthetic field data (number of requested synthetic field data is 
specified by the user) to which random extraction errors are add.
From the resulting set of synthetic field data, the respective conditional 
density probability distribution can be derived in the 'LogLikOptim' Function.

IMPORTANT: Note that for this function to properly work, all variables of the 
SWIFT model should be declared locally or globally. 
}
\usage{ CDPD(beta.i=NULL)}
\arguments{
  \item{beta.i}{
 Description: Allocation parameter that describes the exponential decay of the 
              distribution of absorbing root surfaces for the focal plant;
 Functionality: used in the iSWIFT loglikelihood optimization  
 Structure:   Randomly assigned numeric value within defined ranges as provided 
              in the 'optim' function.
}
\details{
}
\value{
This function will create a matrix of synthetic field data of xylem water isotope 
compositions (D2H and/or D18O)(number of synthetic field data points are user 
defined = 'q') with account for extraction error, and by considering variable 
restrictions of both soil heterogeneity and plant traits which are input values 
of the 'SWIFT' model. Hence, a matrix of two colums (D2H and D18O) and 'q' rows
}
\references{
}
\author{ Hannes De Deurwaerder, Marco D. Visser, Matteo Detto, Pascal Boeckx, 
Felicien Meunier, Pedro Herve Fernandez and Hans Verbeeck}
\note{}
\seealso{}
\examples{
CDPD <- function(beta.i=NULL){
    #--------------------------------------------------------------------#
    # NOTE that for this function to properly work, all variables of the #
    # SWIFT model should be declared locally or globally.                #
    #--------------------------------------------------------------------#
    
  # create empty matrix
  PS <- matrix(NA, SWIFTitterations, 2)
  colnames(PS) <- c('D2H', 'D18O')
  
  # Create matrix of biophysical variables needed for the SWIFT model for which
  # values are random extracted from user defined variable ranges
  BiophysVariables <- VarMatrix(beta.i)
  
  # generate soil water potential heterogeneity
  PSIprofiles <- SoilHeterogeneity(DataType="PSI")
  
  # generate soil water isotope composition heterogeneity
        if(ConsideredIsotopes=='D2H'){
          soilDeuterium <- SoilHeterogeneity(DataType="D2H")
          PS[,'D2H'] <- ProbabilitySpace(beta.i, BiophysVariables, PSIprofiles, 
                                         soilDeuterium, soilOxygen=NULL)
        }else{
          soilOxygen <- SoilHeterogeneity(DataType="D18O")
          PS[,'D18O'] <- ProbabilitySpace(beta.i, BiophysVariables, PSIprofiles, 
                                         soilDeuterium=NULL, soilOxygen)}
        # If dual, deuterium can be obtained from oxygen via LMWL 
        # (as is considered in this example)
        if(ConsideredIsotopes == 'Dual'){
          PS[,'D2H'] <- LMWLslope*PS[,'D18O'] + LMWLintercept 
        }

  # Addition of extraction error to the produced isotope values of SWIFT
  if ( ConsideredIsotopes != 'D18O' ){
    PS[,'D2H']<- PS[,'D2H']+rsn(SWIFTitterations, xi=0, omega=1, alpha=-10^9)}
  if ( ConsideredIsotopes != 'D2H' ){
    PS[,'D18O']<-PS[,'D18O']+rsn(SWIFTitterations, xi=0, omega=.1, alpha=-10^9)}
  
  return(PS)  # Return as matrix
}
}