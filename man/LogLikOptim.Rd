\name{LogLikOptim}
\alias{LogLikOptim}
\title{The Log-likelihood optimizer of iSWIFT}
\description{ 
This function loops through beta.i provided to iSWIFT in order to minimize the 
log-likelyhood function, hence providing the best beta estimate of the studied 
plant group. Field data is matched to the conditional density probability 
distributions of generated synthetic field data using SWIFT model, which are
returned by the 'CDPD' function.

IMPORTANT: Note that for this function to properly work, all variables of the 
SWIFT model should be declared locally or globally. The 'optim' function, 
provided with type of optimizer (here 'Brent') and user provided ranges, will 
optimize this function where a minimum log-likelyhood is sought between field 
data and SWIFT generated synthetic field data.
}
\usage{ LogLikOptim(beta.i=NULL) }
\arguments{
  \item{beta.i}{
 Description: Allocation parameter that describes the exponential decay of the 
              distribution of absorbing root surfaces for the focal plant;
 Functionality: used in the iSWIFT loglikelihood optimization  
 Structure:   Randomly assigned numeric value within defined ranges as provided 
              in the 'optim' function.}
}
\details{
}
\value{
This function will result in an optimized 'beta' value, presenting the best 
estimate for the studied plant groups' beta (i.e., beta_hat)}
\references{
References for optim ranges and initial value can be found in: 
- Jackson et al (1996). A global analysis of root distributions for 
terrestrial biomes. Oecologia, 108(3), 389-411.
}
\author{ Hannes De Deurwaerder, Marco D. Visser, Matteo Detto, Pascal Boeckx, 
Felicien Meunier, Pedro Herve Fernandez and Hans Verbeeck}
\note{}
\seealso{}
\examples{
LogLikOptim <- function(beta.i=x){
    #--------------------------------------------------------------------#
    # NOTE that for this function to properly work, all variables of the #
    # SWIFT model should be declared locally or globally.                #
    #--------------------------------------------------------------------#
    
  # Generate conditional density probability distributions P(δxyl|β,χ) (CDPD)
      Csyn<-CDPD(beta.i)   
      
  # Calculate the loglikelihood
      # when modeled for D2H
      if(ConsideredIsotopes == "D2H"){
          Kernel<-sm.density(Csyn[,"D2H"], eval.points=FieldData[,"D2H"], 
                            display="none")} 
      # When modelled for D18O
      if(ConsideredIsotopes == "D18O"){
          Kernel<-sm.density(Csyn[,"D18O"], eval.points=FieldData[,"D18O"], 
                           display="none")}    
      # combined likelihood for when both are measured
      if(ConsideredIsotopes == "Dual"){
          Kernel<-kde(Csyn[,c("D2H","D18O")], 
                                  eval.points=FieldData[,c("D2H","D18O")])}
  
  PDest <- Kernel$estimate  # extract the best Kernel estimate
  PDest[which(PDest <= 10^-15)]<-10^-15 # remove extremely low values to avoid 
                                        # loglikelihood to produce 'NA' or 'inf'

  return(-sum(log(PDest), na.rm=T))     # return the log likelihood
}
}
