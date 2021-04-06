\name{iSWIFT-package}
\alias{iSWIFT-package}
\alias{iSWIFT}
\docType{package}
\title{
   \packageTitle{iSWIFT}
}
\description{
  \packageDescription{iSWIFT}
}
\details{	
  The DESCRIPTION file:
  \packageDESCRIPTION{iSWIFT}
  \packageIndices{iSWIFT}
This package provide an inverse approach to estimate the 
              distribution of the absorbing root surface from isotope data and 
              plant hydraulic traits considered. While this methodology can be 
              applied on all models combining a mixing model approach with plant
              hydraulic models, here it is applied on the SWIFT model 
              (De Deurwaerder et al. 2020, Biogeosciences).               
              'iSWIFT' provide a set of new function, with at its core the 
              'LogLikOptim' function which optimizes the log-likelihood function
              in order to obtain the best estimate of beta (beta.hat in 
              De Deurwaerder et al (submitted)). Random set of SWIFT variables
              are generated in the 'VarMatrix' module where values are randomly
              extracted from user defined ranges. 'SoilHeterogeneity' will
              generate synthetic varibility in the soil water potential and
              soil water isotope composition gradients with depth for soil
              profile fitting, also required by the SWIFT model. SWIFT model 
              outputs, with extraction error, are subsequently generated in the
              'ProbabilitySpace' function, providing synthetic field data for 
              which conditional desity probability distribution are calculated via 
              the 'CDPD'-function, required to evaluate the log likelihood of
              field data with the SWIFT model generated data for a specific beta
              within the 'LogLik' function.              
              Note that since the iSWIFT package builds further on the SWIFT 
              package, it contains the same functionalities and requires the 
              same variable and parameter declarations as decribe in more detail 
              in the 'SWIFT' package.The SWIFT model functionalities preserved
              are: Function 'SWIFT_SB' and 'SWIFT_H' can respectively be used to 
              calculate the isotopic signature over time at the stem base of the
              plant, or at a by the user defined height and time. Additional 
              functions are 'SoilRootCond', 'PSI0calc' and 'Bprep' which in this 
              order help the user to define (i) the soil to root conductivity 
              for every defined soil layer, (ii) the water potential at stem 
              base for every timestep, and (iii) the root length distribution 
              for every defined soil layer.
}
\author{	
  \packageAuthor{iSWIFT}
  Maintainer: 
  \packageMaintainer{iSWIFT}
}
\references{
   De Deurwaerder et al (in review)
}
