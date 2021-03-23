iSWIFT (0.1) [Release notes](https://github.com/HannesDeDeurwaerder/iSWIFT/).
----------

The ***iSWIFT*** model presents an inverse method of estimating the allocation parameter 'beta' which defines he shape of the absorbing root surface distribution of a plant. Therefore, we apply an inverse modeling approach on the SWIFT model, which combines
a stable water source mixing model with a hydraulic plant model to characterize and track the isotopic composition within a plant. This approach can be pursued with prior knowledge of stable water isotope composition of xylem and soil water, sapflow dynamics and soil water potentials. Details on the iSWIFT model procedure are detailed in the paper by ***De Deurwaerder et al (In Review)***. Note that as this modeling approach build further on the SWIFT model, we recommend the user to familiarize themselves with the SWIFT model, its functionalities and variable requirements. We refer to the corresponding github [SWIFT](https://github.com/HannesDeDeurwaerder/SWIFT) package info or to the  article of [De Deurwaerder et al (2020)](https://doi.org/10.5194/bg-17-4853-2020]). 


## Functions

Functions ***VarMatrix*** and ***SoilHeterogeneity*** can respectively be used to generate variability in *(i)* plant trait variables, and in *(ii)* water potential and isotope composition soil profiles, both required by the SWIFT model to generate synthetic field data with consideration of extraction error via the ***ProbabilitySpace***. Note that variable ranges should be specified by the user. The corresponding conditional density probability distribution of model output for a 'generic' plant with specific 'beta' is obtained using the ***CDPD*** function. The ***LogLikOptim*** model matches the density probability of the real field data with the various generated synthetic density probabilities generated via the SWIFT model. Herefore, the LogLikOptim function provide an automated search for the 'beta' of the probability density that matches the field data best, and hence provides the best estimate of the 'true beta' of the measured plant individuals. 


### Quicklinks

-   [Quick start and tutorials](#quick-start-and-tutorials)
-   [Installation](#installation)
-   [Dependencies](#dependencies)
-   [Examples](#examples)
-   [Thanks](#thanks)
  

## Quick start and tutorials

An in depth description of the inverse methodology for estimating the relative distribution of absorbing root area, applied on the SWIFT model can be found in **De Deurwaerder et al (In Review)**. For more details on the SWIFT model and formula we refer the reader to [De Deurwaerder et al (2020)](https://doi.org/10.5194/bg-17-4853-2020]).  The example provided, in combination with the provided comments in the R-script should enable the user to successfully run the model. 


## Installation

Use the **devtools** package to install the current development version from R.

### Install from Github repository

	# PACKAGE DEVTOOLS REQUIRED
	#--------------------------
	require(devtools)

	# INSTALL PACKAGE FROM GITHUB REPOSITORY
	#---------------------------------------
	install_github("HannesDeDeurwaerder/iSWIFT")
	
	
### Dependencies

iSWIFT has no other dependencies other than the base R installation.


## Examples

The code below defines a simple example of the use of the iSWIFT model, based on
the case study example provided in **De Deurwaerder et al,(In Review)**. 


### Run the iSWIFT model

  # INITIALIZE LIBRARIES 
  #----------------------
  require(SWIFT)
  require(lhs)
  require(GenSA)
  require(sm)
  require(ks)
  require(sn)
  
  # Load associated Datasets 
  #-------------------------
  data(Laussat_SoilWaterIsotopes) # soil water isotope composition 
                                  #[in permil, V-SMOW]
  data(Laussat_SoilWaterPotential)  # soil water potential [in MPa]
  data(SapfluxData)   # sap flux density data [in kg m-2 s-1]
  data(DataLaussat)   # Xylem water isotope composition [in perMil, V-SMOW]
    
  # Initialize global parameters 
  #-----------------------------
  n   <<- 20            # Multiple number of days studied, needed for spin up
                        # of the model.
  tF  <<- 60            # Time frequency of measurements per hour                                    # [in measurments per h].
  tt  <<- seq(0,24*n,length.out = 24*tF*n)     # Discrete time vector [in h].
  dZ  <<- 0.001         # Thickness of sampled layer [in m].
  L   <<- 1             # maximum soil depth [in m].
  Z   <<- seq(dZ,L,dZ)  # Discrete depth vector centered [in m].
  PSIsat <<- -0.153     # Water potential at soil saturation for clay Sand 
                        # soil, [in m H2O](obtained from Clapp and Hornberger
                        # (1978))
  D2Hoffset <<- 6.855   # deuterium offset (calculated from 'DataLaussat')
  Bbeta <<- 0.962       # allocation parameter  that gives the distribution 
                        # of absorbing root surface area of the entire plant 
                        # community. Typical value for Tropical evergreen 
                        # forest obtained from Jackson et al (1996)
  Ltot <<- 10000        # Absorbing root area length Soethe et al (2006)
  BR0 <<- (Ltot*100)/(1-Bbeta^100)
  betaCom <<- Bprep(Bbeta, BR0, Z)  # Absorbing root length distribution of 
                                    # the entire plant community
  
  # Unit conversion factors 
  #------------------------
  CTpsi <<- 101.97      # Conversion factor between MPa and m H2O
  TCOR <<- 24*60*(n-2)  # time correction needed IN THIS EXAMPLE to account 
                        # for spin-up of the model
  cm2_to_m2 <<- 1/10000 # conversion from cm2 to m2!!
  
  # Growth form specific restriction schemes 
  #-----------------------------------------
    
    # Liana restriction scheme - here in function format, declaring globally 
    LianaVariableRanges<- function(){
    DBHrange <<- c(0.63, 17.51)*10^-2 # DBH [in m] 
                                      # (range obtained from 'DataLaussat')
    LArange <<- c(0.14, 0.21)   # Lumen area fraction of sapwood [in m2 m-2]
    Asaprange <<- c(0.5,1.5)    # Sapwood area variability term 'a', see De 
                                # Deurwaerder et al (2021) Table S2
    Krrange <<- c(1,14)*10^-10  # kr, the effective root radial conductivity
                                # [in s-1]
    SFdailyrange <<- c(0.5,1.5) # Daily sapflow variability term 'a', see De 
                                  # Deurwaerder et al (2021) Table S2
    Hosrange <<- c(0,25)        # Sampling height [in m]
    tstudrange <<- c(9,14)*tF   # Timing of sampling
    relSF <<- rep(SapFlow/sum(SapFlow),n) # generic relative SF data 
                                          # repeated over 20 days
      # Growth form specific total absorbing root area calculation function
      ARcalc <<- function(DBH){
        #function parameters, see De Deurwaerder et al (2021)[Table S2]
        Mrl = 0.01; SRA = 40.94746;  cc = -7.094; dd = 1.690 #
        return(Mrl*SRA*exp(cc+dd*log(DBH*10^3)))}
    }
    
    # Tree restriction scheme - here in function format, declaring globally 
    TreeVariableRanges<- function(){
    DBHrange <<- c(9.87, 69.39)*10^-2 # DBH [in m] 
                                            # (range obtained from 'DataLaussat')
          LArange <<- c(0.19,0.41)    # Lumen area fraction of sapwood [in m2 m-2]
          Asaprange <<- c(0.5,1.5)    # Sapwood area variability term 'a', see De 
                                      # Deurwaerder et al (2021) Table S2
          Krrange <<- c(1,14)*10^-10  # kr, the effective root radial conductivity
                                      # [in s-1]
          SFdailyrange <<- c(0.5,1.5) # Daily sapflow variability term 'a', see De 
                                      # Deurwaerder et al (2021) Table S2
          Hosrange <<- c(1.30, 1.30)  # Sampling height [in m]
          tstudrange <<- c(9,14)*tF   # Timing of sampling 
          relSF <<- rep(SapFlow/sum(SapFlow),n) # generic relative SF data 
                                                # repeated over 20 days
            # Growth form specific total absorbing root area calculation function
            ARcalc <<- function(DBH){
              return(exp(0.88*log(pi*(DBH*100/2)^2)-2))}
      }
    
    # Scenario details 
    #-----------------
        ConsideredIsotopes <<- 'Dual'   
                # which isotope studied:'D2H','D18O','Dual'
        SoilWaterPotential <<- Laussat_SoilWaterPotential 
                # allocate soil water potential detail representative for Laussat
        SoilWaterIsotopes <<- Laussat_SoilWaterIsotopes
                # allocate soil water isotope compositions obtained in Laussat
        SWIFTitterations <<- 250 # Number of generated SWIFT datapoints
        LMWLslope <<- 7.748      # Local Meteoric water line slope in Laussat
        LMWLintercept <<- 4.930  # Local Meteoric water line intercepth in Laussat 
        
    # Field data to be evaluated 
    #---------------------------
        FD <<- DataLaussat       # measured isotope data in Laussat
        FD$D2H <- FD$D2H + D2Hoffset   # remove D2Hoffset
        FD_L <- data.frame(D2H=FD$D2H[which(FD$GF=='L')],
                           D18O=FD$D18O[which(FD$GF=='L')])
        FD_T <- data.frame(D2H=FD$D2H[which(FD$GF=='T')],
                           D18O=FD$D18O[which(FD$GF=='T')])
    
    # Run iSWIFT per growth form 
    #---------------------------
        # run for 'LIANAS'
            FieldData <<- FD_L                # Allocate field data of lianas
            SapFlow <<- SapfluxData$Liana     # growth form specific sapflow 
                                              # rates normalized per sapwood area
                                              # Chen et al.(2017)[Fig 5b] 
                                              # [in kg m-2 s-1]
            LianaVariableRanges()             # load growht form specific info
            Beta.hat_L<-optim(0.80, fn=LogLikOptim, lower = 0.6, upper = 0.995, 
                          method = "Brent")   
            Csyn_L<-CDPD(beta.i=Beta.hat_L$par) # SWIFT generated output for best 
                                                # beta estimate
        # run for 'TREES'
            FieldData <<- FD_T            # Allocate field data of trees
            SapFlow <<- SapfluxData$Tree  # growth form specific sapflow 
                                          # rates normalized per sapwood area
                                          # Chen et al.(2017)[Fig 5b] 
                                          # [in kg m-2 s-1]
            TreeVariableRanges()          # load growth form specific info
            Beta.hat_T<-optim(0.95, fn=LogLikOptim, lower = 0.8, upper = 0.995, 
                      method = "Brent")
            Csyn_T<-CDPD(beta.i=Beta.hat_T$par)  # SWIFT generated output for best 
                                                 # beta estimate
    
    # Making a simple plot of the field data and iSWIFT output 
    #----------------------------------------------------------
            xlabel=expression(paste(delta,""^18,"O [","\211",", V-SMOW]"))
            ylabel=expression(paste(delta,""^2,"H [","\211",", V-SMOW]"))
      # Plot Laussat field data
            plot(FD_L[,"D18O"],FD_L[,"D2H"],  ylim=c(-45,10), xlim=c(-8,3),
                  xlab=xlabel, ylab=ylabel, mgp = c(2, 0.8, 0), las=1, pch=20, 
                  col='tan1')
            points( FD_T[,"D18O"],FD_T[,"D2H"], pch=20, col='yellowgreen')
      # plot SWIFT generated 'Tree' isotope data
            points(mean(Csyn_T[,"D18O"]),mean(Csyn_T[,"D2H"]),
                  col='springgreen4',pch=15)
            text(mean(Csyn_T[,"D18O"]),mean(Csyn_T[,"D2H"])+2,
                 c(round(Beta.hat_T$par,3)), col="springgreen4", pos=2, cex=0.7)
            arrows(mean(Csyn_T[,"D18O"]), mean(Csyn_T[,"D2H"])-2*sd(Csyn_T[,"D2H"]),
                   mean(Csyn_T[,"D18O"]), mean(Csyn_T[,"D2H"])+2*sd(Csyn_T[,"D2H"]), 
                   length = 0, col="springgreen4", lwd=1.5)
            arrows(mean(Csyn_T[,"D18O"])-2*sd(Csyn_T[,"D18O"]), mean(Csyn_T[,"D2H"]),
                   mean(Csyn_T[,"D18O"])+2*sd(Csyn_T[,"D18O"]), mean(Csyn_T[,"D2H"]), 
                   length = 0, col="springgreen4", lwd=1.5)
      # plot SWIFT generated 'Liana' isotope data
            points(mean(Csyn_L[,"D18O"]),mean(Csyn_L[,"D2H"]), 
                    col='sienna3', pch=15)
            text(mean(Csyn_L[,"D18O"]),mean(Csyn_L[,"D2H"])+2,
                 c(round(Beta.hat_L$par,3)), col="sienna3", pos=2, cex=0.7)
            arrows(mean(Csyn_L[,"D18O"]), mean(Csyn_L[,"D2H"])-2*sd(Csyn_L[,"D2H"]),
                   mean(Csyn_L[,"D18O"]), mean(Csyn_L[,"D2H"])+2*sd(Csyn_L[,"D2H"]), 
                   length = 0, col="sienna3", lwd=1.5)
            arrows(mean(Csyn_L[,"D18O"])-2*sd(Csyn_L[,"D18O"]), mean(Csyn_L[,"D2H"]),
                   mean(Csyn_L[,"D18O"])+2*sd(Csyn_L[,"D18O"]), mean(Csyn_L[,"D2H"]), 
                   length = 0, col="sienna3", lwd=1.5)
      # add legend
            legend('topleft', c(expression(paste('C'[liana])),
                   expression(paste(chi,''[liana])), expression(paste('C'[tree])),
                   expression(paste(chi,''[tree]))),  col=c('tan1','sienna3',
                   'yellowgreen','springgreen4'), pch=20, ncol=2, bty='n')
    
## Thanks
Special thanks all other co-authors for their help in developing the iSWIFT model.

