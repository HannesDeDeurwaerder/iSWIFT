#'##############################################################################
#'##############################################################################
#'
#' @title iSWIFT (inverse modeling methodology applied to the SWIFT model)
#' 
#' @param ARtot Description: The total absorbing root surface area of the focal 
#'                           plant [in m2];
#'              Functionality: Standard SWIFT parameter  
#'              Structure:   One value, representative of the studied plant. 
#' @param Ax    Description: The total lumen area [in m2];  
#'              Functionality: Standard SWIFT parameter  
#'              Structure:   One value, representative of the studied plant.
#' @param B     Description: The community plant Root length density at every 
#'                           discrete soil layer 
#'                           - Results from the function Bprep - [in m m-3]; 
#'              Functionality: Standard SWIFT parameter  
#'              Structure:   A discrete vector of n elements, with n the number 
#'                           of soil layers.
#' @param beta  Description: Allocation parameter that describes the exponential
#'                           decay of the distribution of absorbing root 
#'                           surfaces for the focal plant; 
#'              Functionality: Standard SWIFT parameter  
#'              Structure:   One value, derived from Jackson et al. (1996), 
#'                           table 1.
#' @param beta.i  Description: Allocation parameter that describes the exponential
#'                           decay of the distribution of absorbing root 
#'                           surfaces for the focal plant;
#'              Functionality: used in the iSWIFT loglikelihood optimization  
#'              Structure:   Randomly assigned integer value within defined 
#'                           ranges as provided in the 'optim' function.
#' @param betaCom  Description: Factor in the root length density distribution 
#'                           function of the plant community 
#'                           (required to distribute the root length over all 
#'                           soil layers) 
#'                 Functionality: Standard SWIFT parameter  
#'                 Structure:   One value, derived from Jackson et al. (1996), 
#'                              table 1.
#' @param D2Hxylem  Description: The deuterium signature at stem base - derived 
#'                               from the function SWIFT_SB [in permil, VSMOW];  
#'                  Functionality: Standard SWIFT parameter  
#'                  Structure:   A discrete vector of n elements, with n the 
#'                               number of timesteps sampled.
#' @param D2Hsoil   Description: Vector representing the stable water isotopic 
#'                               siganture profile with depth 
#'                               [in permil, VSMOW]; 
#'                  Functionality: Standard SWIFT parameter  
#'                  Structure:   A discrete vector of n elements, with n the 
#'                               number of soil layers.
#' @param dZ    Description: Thickness of the soil layer [in m];  
#'              Structure:   One value, Samples should be taken at discrete 
#'                           depths which have the same thickness.
#' @param hos   Description: A vector of studied measurement heights [in m]; 
#'              Functionality: Standard SWIFT parameter 'hom'  
#'              Structure:   The size of the vector is defined by the interest 
#'                           of the user. 
#' @param k     Description: Plant specific total soil-to-root conductance at
#'                           each specific specific soil layer [in s-1];  
#'              Functionality: Standard SWIFT parameter  
#'              Structure:   A vector - results from the function SoilRootCond
#'                           - and has n elements corresponding to the number of
#'                           soil layers.  
#' @param kr    Description: The root membrane permeability [s-1];  
#'              Functionality: Standard SWIFT parameter  
#'              Structure:   One value, representative of the studied plant.
#' @param PSIs  Description: Soil water potential at the each specific soil 
#'                           layer [in m]; 
#'              Functionality: Standard SWIFT parameter  
#'              Structure:   A vector of n element, where n corresponds to the 
#'                           number of soil layers;  
#'              Note:        The water potentials per soil layer is currently 
#'                           considered static in time.
#' @param R0    Description: Factor in the root lenght density distribution per 
#'                           unit of soil, for entire plant community;
#'              Functionality: Standard SWIFT parameter  
#'              NOTE:        Provide a Positive value!!  
#'              Structure:   One value,representative of the studied plant, here
#'                           derived from Huang et al, 2017.
#' @param SF    Description: Instantaneous sap flow over time [in m^3 s-1];  
#'              Functionality: Standard SWIFT parameter  
#'              Structure:   A vector of n elements, where n corresponds to the 
#'                           number of timesteps.
#' @param Soiltype  Description:  Type of soil layer (Sand; Loamy Sand; 
#'                                Sandy Loam; Silt Loam; Loam; Sandy Clay Loam;
#'                                Silty Clay Loam; Clay Loam;Sandy Clay;
#'                                Silty Clay; Clay); 
#'                  Functionality: Standard SWIFT parameter  
#'                  Structure:    A Character string. 
#' @param tt    Description:   Cumulative time steps [in s]; 
#'              Functionality: Standard SWIFT parameter 't' 
#'              Structure:     A vector.
#' @param tstud Description:   By the user defined timesteps to which an output 
#'                             is required [in second after initialisation]; 
#'              Functionality: Standard SWIFT parameter  
#'              Structure:     A vector with n elements defined by the user;  
#'              Note:          When an output is needed for all values, then 
#'                             tstud = t.
#' @param tF    Description:   Measurement frequence per hour [in measurement 
#'                             points per hour];  
#'              Functionality: Standard SWIFT parameter  
#'              Structure:     One value, defined by the user; 
#'              Note:          Lower numbers will increase the speed of the 
#'                             model, but for graphical display, higher values 
#'                             are prefered.
#' @param Z     Description:   A vector of Soil depth [in m]; 
#'              Functionality: Standard SWIFT parameter  
#'              Structure:     A discrete vector  of n elements, where n 
#'                             corresponds to the number of soil layers.
#'
#'   
#' @author Hannes P.T. De Deurwaerder, Marco D. Visser, Felicien Meunier, 
#' Matteo Detto, Pascal Boeckx, Pedro Herve Fernandez and  Hans Verbeeck.
#' 
#' @examples
#' \dontrun{
#' 
#'  # Initialize of Libraries ---------------------------------------------------
#'  require('SWIFT')
#'  require('lhs')
#'  require('GenSA') 
#'  require('sm')
#'  require('ks')
#'  require('sn')
#'  
#'  # Load associated Datasets --------------------------------------------------
#'  data(Laussat_SoilWaterIsotopes)  # soil water isotope composition 
#'                                   # [in permil, V-SMOW]
#'  data(Laussat_SoilWaterPotential) # soil water potential [in MPa]
#'  data(SapfluxData)  # sap flux density data [in kg m-2 s-1]
#'  data(DataLaussat) # Xylem water isotope composition [in perMil, V-SMOW]
#'  
#'  # Initialize global parameters ----------------------------------------------
#'  n   <<- 20            # Multiple number of days studied, needed for spin up 
#'                        # of the model.
#'  tF  <<- 60            # Time frequency of measurements per hour 
#'                        # [in measurments per h].
#'  tt  <<- seq(0,24*n,length.out = 24*tF*n)     # Discrete time vector [in h].
#'  dZ  <<- 0.001         # Thickness of sampled layer [in m].
#'  L   <<- 1             # maximum soil depth [in m].
#'  Z   <<- seq(dZ,L,dZ)  # Discrete depth vector centered [in m].
#'  PSIsat <<- -0.153     # Water potential at soil saturation for clay Sand 
#'                        # soil, [in m H2O](obtained from Clapp and Hornberger
#'                        # (1978))
#'  D2Hoffset <<- 6.855   # deuterium offset (calculated from 'DataLaussat')
#'  Bbeta <<- 0.962       # allocation parameter  that gives the distribution 
#'                        # of absorbing root surface area of the entire plant 
#'                        # community. Typical value for Tropical evergreen 
#'                        # forest obtained from Jackson et al (1996)
#'  Ltot <<- 10000        # Absorbing root area length Soethe et al (2006)
#'  BR0 <<- (Ltot*100)/(1-Bbeta^100)
#'  betaCom <<- Bprep(Bbeta, BR0, Z)  # Absorbing root length distribution of 
#'                                    # the entire plant community
#'  
#'  # Unit conversion factors ---------------------------------------------------
#'  CTpsi <<- 101.97      # Conversion factor between MPa and m H2O
#'  TCOR <<- 24*60*(n-2)  # time correction needed IN THIS EXAMPLE to account 
#'                        # for spin-up of the model
#'  cm2_to_m2 <<- 1/10000 # conversion from cm2 to m2!!
#'  
#'  # Growth form specific restriction schemes ----------------------------------
#'  # Liana restriction scheme - here in function format, declaring globally 
#'  LianaVariableRanges<- function(){
#'    DBHrange <<- c(0.63, 17.51)*10^-2 # DBH [in m] 
#'    # (range obtained from 'DataLaussat')
#'    LArange <<- c(0.14, 0.21)   # Lumen area fraction of sapwood [in m2 m-2]
#'    Asaprange <<- c(0.5,1.5)    # Sapwood area variability term 'a', see De 
#'                                # Deurwaerder et al (2021) Table S2
#'    Krrange <<- c(1,14)*10^-10  # kr, the effective root radial conductivity
#'                                # [in s-1]
#'    SFdailyrange <<- c(0.5,1.5) # Daily sapflow variability term 'a', see De 
#'                                # Deurwaerder et al (2021) Table S2
#'    Homrange <<- c(0,25)        # Sampling height [in m]
#'    tstudrange <<- c(9,14)*tF   # Timing of sampling
#'    relSF <<- rep(SapFlow/sum(SapFlow),n) # generic relative SF data 
#'                                          # repeated over 20 days
#'    # Growth form specific total absorbing root area calculation function
#'    ARcalc <<- function(DBH){
#'      #function parameters, see De Deurwaerder et al (2021)[Table S2]
#'      Mrl = 0.01; SRA = 40.94746;  cc = -7.094; dd = 1.690 #
#'      return(Mrl*SRA*exp(cc+dd*log(DBH*10^3)))}
#'  }
#'  # Tree restriction scheme - here in function format, declaring globally 
#'  TreeVariableRanges<- function(){
#'    DBHrange <<- c(9.87, 69.39)*10^-2 # DBH [in m] 
#'                                      # (range obtained from 'DataLaussat')
#'    LArange <<- c(0.19,0.41)    # Lumen area fraction of sapwood [in m2 m-2]
#'    Asaprange <<- c(0.5,1.5)    # Sapwood area variability term 'a', see De 
#'                                # Deurwaerder et al (2021) Table S2
#'    Krrange <<- c(1,14)*10^-10  # kr, the effective root radial conductivity
#'                                # [in s-1]
#'    SFdailyrange <<- c(0.5,1.5) # Daily sapflow variability term 'a', see De 
#'                                # Deurwaerder et al (2021) Table S2
#'    Hosrange <<- c(1.30, 1.30)  # Sampling height [in m]
#'    tstudrange <<- c(9,14)*tF   # Timing of sampling 
#'    relSF <<- rep(SapFlow/sum(SapFlow),n) # generic relative SF data 
#'                                          # repeated over 20 days
#'    # Growth form specific total absorbing root area calculation function
#'    ARcalc <<- function(DBH){
#'      return(exp(0.88*log(pi*(DBH*100/2)^2)-2))}
#'  }
#'  
#'  # Scenario details ----------------------------------------------------------
#'  ConsideredIsotopes <<- 'Dual'   
#'            # which isotope studied:'D2H','D18O','Dual'
#'  SoilWaterPotential <<- Laussat_SoilWaterPotential 
#'            # allocate soil water potential detail representative for Laussat
#'  SoilWaterIsotopes <<- Laussat_SoilWaterIsotopes
#'            # allocate soil water isotope compositions obtained in Laussat
#'  SWIFTitterations <<- 250 # Number of generated SWIFT datapoints
#'  LMWLslope <<- 7.748      # Local Meteoric water line slope in Laussat
#'  LMWLintercept <<- 4.930  # Local Meteoric water line intercepth in Laussat 
#'  
#'  # Field data to be evaluated ------------------------------------------------
#'  FD <<- DataLaussat # measured isotope data in Laussat
#'  FD$D2H <- FD$D2H + D2Hoffset   # remove D2Hoffset
#'  FD_L <- data.frame(D2H=FD$D2H[which(FD$GF=='L')],
#'                     D18O=FD$D18O[which(FD$GF=='L')])
#'  FD_T <- data.frame(D2H=FD$D2H[which(FD$GF=='T')],
#'                     D18O=FD$D18O[which(FD$GF=='T')])
#'  
#'  # Run iSWIFT per growth form-------------------------------------------------
#'  # run for 'LIANAS'
#'  FieldData <<- FD_L      # Allocate field data of lianas
#'  SapFlow <<- SapfluxData$Liana     # growth form specific sapflow 
#'                                    # rates normalized per sapwood area
#'                                    # Chen et al.(2017)[Fig 5b] 
#'                                    # [in kg m-2 s-1]
#'  LianaVariableRanges()             # load growht form specific info
#'  Beta.hat_L<-optim(0.80, fn=LogLikOptim, lower = 0.6, upper = 0.995, 
#'                    method = "Brent")   
#'  Csyn_L<-CDPD(beta.i=Beta.hat_L$par) # SWIFT generated output for best 
#'                                      # beta estimate
#'  # run for 'TREES'
#'  FieldData <<- FD_T      # Allocate field data of trees
#'  SapFlow <<- SapfluxData$Tree  # growth form specific sapflow 
#'                                # rates normalized per sapwood area
#'                                # Chen et al.(2017)[Fig 5b] 
#'                                # [in kg m-2 s-1]
#'  TreeVariableRanges()  # load growth form specific info
#'  Beta.hat_T<-optim(0.95, fn=LogLikOptim, lower = 0.8, upper = 0.995, 
#'                    method = "Brent")
#'  Csyn_T<-CDPD(beta.i=Beta.hat_T$par)  # SWIFT generated output for best 
#'                                       # beta estimate
#'  
#'  # Making a simple plot of the field data and iSWIFT output -----------------
#'  xlabel=expression(paste(delta,""^18,"O [","\211",", V-SMOW]"))
#'  ylabel=expression(paste(delta,""^2,"H [","\211",", V-SMOW]"))
#'  # Plot Laussat field data
#'  plot(FD_L[,"D18O"],FD_L[,"D2H"],  ylim=c(-45,10), xlim=c(-8,3),
#'      xlab=xlabel, ylab=ylabel, mgp = c(2, 0.8, 0), las=1, pch=20, 
#'       col='tan1')
#'  points( FD_T[,"D18O"],FD_T[,"D2H"], pch=20, col='yellowgreen')
#'  # plot SWIFT generated 'Tree' isotope data
#'  points(mean(Csyn_T[,"D18O"]),mean(Csyn_T[,"D2H"]),
#'         col='springgreen4',pch=15)
#'  text(mean(Csyn_T[,"D18O"]),mean(Csyn_T[,"D2H"])+2,
#'       c(round(Beta.hat_T$par,3)), col="springgreen4", pos=2, cex=0.7)
#'  arrows(mean(Csyn_T[,"D18O"]), mean(Csyn_T[,"D2H"])-2*sd(Csyn_T[,"D2H"]),
#'         mean(Csyn_T[,"D18O"]), mean(Csyn_T[,"D2H"])+2*sd(Csyn_T[,"D2H"]), 
#'         length = 0, col="springgreen4", lwd=1.5)
#'  arrows(mean(Csyn_T[,"D18O"])-2*sd(Csyn_T[,"D18O"]), mean(Csyn_T[,"D2H"]),
#'         mean(Csyn_T[,"D18O"])+2*sd(Csyn_T[,"D18O"]), mean(Csyn_T[,"D2H"]), 
#'         length = 0, col="springgreen4", lwd=1.5)
#'  # plot SWIFT generated 'Liana' isotope data
#'  points(mean(Csyn_L[,"D18O"]),mean(Csyn_L[,"D2H"]), 
#'         col='sienna3', pch=15)
#'  text(mean(Csyn_L[,"D18O"]),mean(Csyn_L[,"D2H"])+2,
#'       c(round(Beta.hat_L$par,3)), col="sienna3", pos=2, cex=0.7)
#'  arrows(mean(Csyn_L[,"D18O"]), mean(Csyn_L[,"D2H"])-2*sd(Csyn_L[,"D2H"]),
#'         mean(Csyn_L[,"D18O"]), mean(Csyn_L[,"D2H"])+2*sd(Csyn_L[,"D2H"]), 
#'         length = 0, col="sienna3", lwd=1.5)
#'  arrows(mean(Csyn_L[,"D18O"])-2*sd(Csyn_L[,"D18O"]), mean(Csyn_L[,"D2H"]),
#'         mean(Csyn_L[,"D18O"])+2*sd(Csyn_L[,"D18O"]), mean(Csyn_L[,"D2H"]), 
#'         length = 0, col="sienna3", lwd=1.5)
#'  # add legend
#'  legend('topleft', c(expression(paste('C'[liana])),
#'                      expression(paste(chi,''[liana])), expression(paste('C'[tree])),
#'                      expression(paste(chi,''[tree]))),  col=c('tan1','sienna3',
#'                      'yellowgreen','springgreen4'), pch=20, ncol=2, bty='n')
#'
#' load()
#'	}
#'	
#' @return Provided by function 
#'         - LogLikOptim:       The Log-likelihood function to be optimized
#'         - CDPD:              Generating conditional density probability 
#'                              distributions via the SWIFT model  
#'         - VarMatrix:         Create variable matrix of SWIFT variables using
#'                              user defined ranges
#'         - ProbabilitySpace:  Generating synthetic field data with SWIFT
#'         - SoilHeterogeneity: Generating some theoretical soil heterogeneity
#' @export
#' 
#'##############################################################################
#'##############################################################################


LogLikOptim <- function(beta.i=x){
#===============================================================================
#                         Log-likelihood function 
#===============================================================================
  # Function description: This is the log-likelihood function

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
  
  PDest <- Kernel$estimate #obtain best Kernel estimate
  PDest[which(PDest <= 10^-15)]<-10^-15 # remove extremely low values to avoid 
                                        # loglikelihood to produce 'NA' or 'inf'

  return(-sum(log(PDest), na.rm=T))  # return the log likelihood
}

################################################################################
################################################################################

CDPD <- function(beta.i=NULL){
#===============================================================================
#              Produce conditional density probability distribution 
#===============================================================================
  # Function description: Feeding data to the SWIFT model to generate synthetic
  # field data, representing the conditional density probability distribution
  
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

################################################################################
################################################################################

VarMatrix <- function(beta.i=NULL){
#-------------------------------------------------------------------------------
#                     Create a variable matrix
#===============================================================================
  # Function description: creating a matrix of SWIFT variables extracted 
  # at random from user defined variable ranges
  
  # create empty data matrix
  SWIFTvars=matrix(NA,SWIFTitterations,ncol=11)
  colnames(SWIFTvars)<-c("DBH","Area_alg", "SFtot_alg","ARtot",
                          "LA", "Asapwood","Kr","SFdaily",'beta', 'tstud','hos') 
  
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


################################################################################
################################################################################


ProbabilitySpace <- function(beta.i, BiophysVariables, PSIprofiles, 
                             soilDeuterium, soilOxygen){
#-------------------------------------------------------------------------------
#                     Create synthetic field data with SWIFT
#===============================================================================
  # Function description: feed data to SWIFT to create synthetic field data
  # hence, producing the probability spaces
  
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

################################################################################
################################################################################

SoilHeterogeneity <-function(DataType=NULL){
#-------------------------------------------------------------------------------
#             create some theoretical soil heterogeneity
#===============================================================================
  # Function description: creating some heterogeity in soil water potential and
  # soil isotopic curves to simulate soil local soil heterogeneity
  
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

################################################################################
################################################################################

#' 'DataLaussat' -> Trees and Liana xylem isotopic signature in Laussat
#' 
#' General Description:
#' This dataset contains isotopic data and sampling information of lianas
#' and trees sampled in Laussat, French Guiana, campaing 2017. Data consist of
#' timing of sampling, the plants' diameter [in cm], growth form type ('L' =
#' liana; 'T' = tree) and the isotopic deuterium and oxygen 18 composition
#' of sampled xylem water [in permil, V-SMOW]
#'
#' This dataset may not be used without previous written permission from the
#' owner.
#'
#' Ownership: Hannes De Deurwaerder  <Hannes_de_deurwaerder@hotmail.com>
#' 
#' Research Location:
#' Within the research forest of LAUSSAT, French Guyana.
#' Sampled species all are located on 'terra firme' subsoil.
#'
#' Dataset abstract: 
#' These data are being used to evaluate differences in xylem water isotopic
#' signatures between lianas and trees, and to inversely obtain the absorbing
#' root surface distribution of both growth forms
#' 
#' Methods:
#' The methods of sampling are in much detail described in 
#' De Deurwaerder et al. (submitted) 
#' 
#' The dataset contains the following labels (columns):
#' \itemize{
#' \item SoilType (Character) Indication of the soil type it was sampled
#' ('C' = 'terra firme clay soils).
#' \item GF (Character) indication of growth form ('L' = Liana; 'T'=Tree).
#' \item SampleTime (Numeric) decimal time of the day when sampled.
#' \item D2H (Numeric) Xylem deuterium isotopic composition [in permil, VSMOW]. 
#' \item D18O (Numeric) Xylem O18 isotopic composition [in permil, VSMOW]. 
#' \item DBH (Numeric) Xylem O18 isotopic composition [in cm]. 
#' 
#' \D1O
#' }
#'
#' @docType data
#' @keywords datasets
#' @name DataLaussat
#' @usage data(DataLaussat)
#' @format 
#' 
NULL


################################################################################
################################################################################

#' 'Laussat_SoilWaterPotential' -> Soil water potential measurements with depth
#' 
#' General Description:
#' This dataset provides soil water potential measurements averages and their 
#' standard deviations for multiple depths. This data is obtained in Paracou,
#' French Guiana, under similar edaphic and environmental conditions as the
#' Laussat sampling plots, hence considered representative for the Laussat setup
#'
#' This dataset may not be used without previous written permission from the
#' owner.
#'
#' Ownership: Hannes De Deurwaerder  <Hannes_de_deurwaerder@hotmail.com>
#' 
#' Research Location:
#' Within the research forest of PARACOU, French Guyana.
#' Compiled averages of 3 distinct setups of multiple psichrometers installed at
#' various depths.
#'
#' Dataset abstract: 
#' These data is being used to derive representative soil water profiles in the 
#' iSWIFT analysis
#' 
#' Methods:
#' The methods of sampling is described in much detail in 
#' De Deurwaerder et al. (submitted) 
#' 
#' The dataset contains the following labels (columns):
#' \itemize{
#' \item Depth (Numeric) Indication of soil depth [in m]
#' \item AvgPsi (Numeric) Average soil water potential at that depth, averaged 
#' over the monitoring period length and distinct setups [MPa]
#' \item sdPsi (Numeric) Standard deviation of soil water potential values at 
#' that depth, considering the monitoring period length and distinct setups [MPa]. 
#' 
#' \D1O
#' }
#'
#' @docType data
#' @keywords datasets
#' @name Laussat_SoilWaterPotential
#' @usage data(Laussat_SoilWaterPotential)
#' @format 
#' 
NULL
      
################################################################################
################################################################################


#' 'Laussat_SoilWaterIsotopes' -> Soil water isotopic composition obtained
#' in Laussat, French Guiana [in permil, VSMOW].
#' 
#' General Description:
#' This dataset contains isotopic data sampled in Laussat, French Guiana, in 
#' close vicinity with the sampled lianas and trees. Data is pooled from various
#' soil pits where soil is sampled at multiple depths.
#'
#' This dataset may not be used without previous written permission from the
#' owner.
#'
#' Ownership: Hannes De Deurwaerder  <Hannes_de_deurwaerder@hotmail.com>
#' 
#' Research Location:
#' Within the research forest of LAUSSAT, French Guyana.
#' Sampled soils are located on 'terra firme' subsoil and in close vicinity of 
#' the sampled lianas and trees.
#'
#' Dataset abstract: 
#' These data is being used to derive representative soil water isotopic 
#' composition profiles required for the iSWIFT analysis.
#' 
#' Methods:
#' The methods of sampling are in much detail described in 
#' De Deurwaerder et al. (submitted) 
#' 
#' The dataset contains the following labels (columns):
#' \itemize{
#' \item Depth (Numeric) Indication of soil depth [in m]
#' \item AvgD2H (Numeric) Averaged soil water deuterium composition 
#' [in permil, VSMOW]
#' #' \item AvgD18O (Numeric) Averaged soil water 18oxygen isotope composition 
#' [in permil, VSMOW]  
#' \item sdD2H (Numeric) Standard deviation of soil water deuterium composition 
#' [in permil, VSMOW]
#' \item sdD18O (Numeric) Standard deviation of soil water 18oxygen isotope 
#' composition [in permil, VSMOW]#'  
#' \D1O
#' }
#'
#' @docType data
#' @keywords datasets
#' @name Laussat_SoilWaterIsotopes
#' @usage data(Laussat_SoilWaterIsotopes)
#' @format 
#' 
NULL
################################################################################
################################################################################


#' 'SapfluxData' -> Sapflux densities for Trees and Liana normalized per sapwood
#' area. Data adapted from Chen et al. (2017) [Fig 5b, in kg m-2 s-1]
#' 
#' General Description:
#' This dataset shows the result of curve fitting optimization of the Chen et 
#' al. (2017) sapflow rates, extracted from fig 5b. Sapflow rates are provided
#' for each minute, but expressed in kg m-2 s-1.
#'
#' This dataset was digitized from the published paper by Chen et al. (2017). 
#' The original dataset can be obtained from the authors of the paper, and may 
#' not be used without previous written permission from the owner.
#'
#' Ownership: The orrginal belongs to the authors of the Chen et al. (2017) 
#' paper, and may be obtained on request to <zjl@xtbg.org.cn> or 
#' <caokf@xtbg.ac.cn>
#' 
#' #' Research Location:
#' This data is obtained in Xishuangbanna Tropical Botanical Garden 
#' (21°540N, 101°460E) in south of Yunnan Province, 
#' China.
#'
#' Dataset abstract: 
#' These data show the sap flux density per growth form, obtained from 4 
#' distinct tropical liana and 4 distinct tropical tree species (3 to 5 indiv.
#' per species). 
#' 
#' Methods:
#' Sap ﬂow was monitored during the end of the wet season (8–21 October 2012) 
#' using granier sensors. More details are provided in the respective paper 
#' by Chen et al. (2017). Functional Ecology. doi: 10.1111/1365-2435.12724 
#' 
#' The dataset contains the following labels (columns):
#' \itemize{
#' \item timeStep (Integer) Timestep in minutes during the day.
#' \item timeDec (Numeric) Decimal notation of the time of the day.
#' \item time (Character) Hour:Minute notation of the time of the day.
#' \item Liana (Numeric) Lianas' sap flux density per timestep [in kg m-2 s-1]. 
#' \item Tree (Numeric) Trees' sap flux density per timestep [in kg m-2 s-1]. 
#' 
#' \D1O
#' }
#'
#' @docType data
#' @keywords datasets
#' @name SapfluxData
#' @usage data(SapfluxData)
#' @format 
#' 
NULL


