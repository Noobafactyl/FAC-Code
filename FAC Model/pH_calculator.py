'recovery save'
import parameter_input as pi
import constant_values as cv
import numpy as np

class BulkpHCalculator(object):
    def __init__(self, Section):
        #Equilibrium and Debye-Huckel constants - polynomials as a function of temperature 
        #Coeff1*x^4 + Coeff2*x^3 + Coeff3*x^2 + Coeff4*x + Coeff5, depends on # elements in coeff list
        self.DebyeHuckelConstant=(np.polyval(pi.DebyeHuckel,Section.Celsius)) 
        self.k_W=10**(np.polyval(pi.k_W,Section.Kelvin))
        self.k_Li=10**(np.polyval(pi.k_Li,Section.Kelvin)) 
        self.k_FeOH=10**(np.polyval(pi.k_FeOH,Section.Kelvin))
        self.k_FeOH2=10**(np.polyval(pi.k_FeOH2,Section.Kelvin))
        self.k_FeOH3=10**(np.polyval(pi.k_FeOH3,Section.Kelvin))
        self.k_NiOH=10**(np.polyval(pi.k_NiOH,Section.Kelvin))
        self.k_NiOH2=10**(np.polyval(pi.k_NiOH2,Section.Kelvin))
        self.k_NiOH3=10**(np.polyval(pi.k_NiOH3,Section.Kelvin))
        #Numerical method to calculate the pH based on known ConcentrationLi, eq'm constants(T), and an initial guess for ConcentrationH and activity coeff.    
    
    
    def BulkProtonConcentrationentration(self): #Bulk pH calculation 
        ConcentrationH = 0.000001  #initial guess [mol/kg]
        gamma_1 = 1 #initial guess 
        ConcentrationOH =  self.k_W / ((gamma_1**2) * (ConcentrationH))
        #At high temp, LiOH doesn't dissociate 100% - eq'm established: LiOH(aq) <-> Li+(aq) + OH-(aq)
        #KLi = ([Li+]gamma_1 *[OH-]gamma_1)/[LiOH]; ConcentrationLiTotal = ConcentrationLi  + LiConcentrationOH (sub in eq'm expression for LiConcentrationOH
        ConcentrationLi = self.k_Li * cv.ConcentrationLiTotal / (self.k_Li + ((gamma_1**2) * ConcentrationOH))
        #k_W = ConcentrationH*gamma_1*ConcentrationOH*gamma_1;   ConcentrationOH = k_W/ConcentrationH*gamma_1**2
        #H+ + Li+ = OH-;                    ConcentrationH = ConcentrationLi + k_W/ConcentrationH*gamma_1**2
        for i in range(20): #no more than 10 iterations should really be necessary for convergence at the provided error  
            FH = ConcentrationH**2 + (ConcentrationH * ConcentrationLi) - (self.k_W / (gamma_1**2)) #function of H
            DFH = 2 * ConcentrationH + ConcentrationLi       #derivative of function
            NewConcentrationH = ConcentrationH - (FH / DFH)
            RE = abs((NewConcentrationH - ConcentrationH) / NewConcentrationH)
            ConcentrationH = NewConcentrationH
            #print RE #optional value check of ConcentrationH at each iteration 
            ConcentrationOH = self.k_W / ((gamma_1**2) / ConcentrationH)
            ConcentrationLi = (self.k_Li * cv.ConcentrationLiTotal) / (self.k_Li + ((gamma_1**2) * ConcentrationOH))
            IonicStrength = ((1**2) * ConcentrationH + (1**2) * ConcentrationLi + (1**2) * ConcentrationOH) / 2
            #Davies equation loggamma_1 = -DebyeHuckConst*(z**2)*[(sqrt(I)/(1+sqrt(I)))-beta*I]
            gamma_1 = 10**(-self.DebyeHuckelConstant * (1**2) * (((IonicStrength ** 0.5) / (1 + (IonicStrength**0.5))) - 0.2 * IonicStrength))
            #All entries in ConcentrationH list must meet error convergence requirement (no matter if SG 22 element list or Inlet 7 element list)
            if RE.all() < 0.00000001: break 
        #print i #prints number of iterations before error minimized to desired level and exits loop
        #print -np.log10(ConcentrationH)
        #return ConcentrationH
        return ConcentrationH
