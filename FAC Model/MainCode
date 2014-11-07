import csv
import numpy as np
 
FileName1 = open('GeneralParameters.txt','r') #Don't need to specify file path if docs are in the same folder as this .py file     
Reader1 = list(csv.reader(FileName1, delimiter = ',')) #Assign file data to list, Reader[row][column], delimiter = comma 
FileName2 = open('ElectrochemicalParameters.txt','r'); Reader2 = list(csv.reader(FileName2, delimiter = ','))
#------------------------------------------------------------------------------------------------------------------------------- 
#class = data structure -holds data and methods to process data (i.e., functions)variables inside class=attributes; functions=methods
#Constants 
Beta= 0.5 #Symmetry coefficient
kH2=0.00078 #Henry#s law constant for H2 @ 298.15 K [mol/L*atm]
F=96485.3   #Faraday's constant [C/mol]
n=2 #number of electrons transferred 
R=8.314 #Ideal gas constant [J/K *mol]
kb=1.38E-23 #Boltzman constant
hp=6.626E-34 #Planck's constant [m^2*kg/s]
FeMolarMass=55.847 #[g/mol]
Fe3O4MolarMass=231.541 #[g/mol]
NiMolarMass=58.694 #[g/mol]
CrMolarMass=51.996 #[g/mol]
NiFe2O4MolarMass=234.388 #[g/mol]
FeCr2O4MoladMass=223.779 #[g/mol] 
H2MolarMass=2.016 #[g/mol]
Fe3O4Density=5.2 #[g/cm3]
NiFe2O4Density=5.368 #[g/cm3]
FeCr2O4Density=4.8 #[g/cm3]
H2Density=8.988e-5 #[g/cm3]
FeDensity=7.86 #[g/cm3]
NiDensity=8.908 #[g/cm3]
Alloy800Density=7.94 #[g/cm3]
FeDiffusivity=0.00041 #[g/cm3]
NiDiffusivity=0.00041 #[cm^2/s]
Fe3O4Porosity=0.25
FeCr2O4Porosity=0.15
Fe3O4Tortuosity=1.8
FeCr2O4Tortuosity=1.2 #both unitless                   
KpFe3O4=0.007 #[cm/s]
KdFe3O4=0.044 #[cm/s] Mohajery 
KpNi_NiFe2O4=0.0008 #[cm/s]
KdNi_NiFe2O4=0.002 #[cm/s]
KpFe_NiFe2O4=0.004 #[cm/s]
KdFe_NiFe2O4=0.0015 #[cm/s]
KpNi=0.001 #[cm/s]
KdNi=0.008 #[cm/s]                 
Kdep  = 0.07 #Burril [cm/s]
TimeIncrement=3600 #[s]
FracNi_NiFe2O4=0.25
FracFe_Fe3O4=0.723
InnerFeCr2O4Thickness=0.0000048  #[g/cm^2] (0.01 um or 10 nm) fixed passivation layer 
ConcLiTotal=0.00022585 #[mol/L]

#reads polynomial coefficients and sets up as lists                             
DebyeHuckPolynomial=[float(Reader1[i][0]) for i in range(64,69)]; KwPolynomial=[float(Reader1[i][1]) for i in range(64,69)] 
KLiPolynomial=[float(Reader1[i][2]) for i in range(64,67)]; KFeOHPolynomial=[float(Reader1[i][3]) for i in range(64,69)] 
KFeOH2Polynomial= [float(Reader1[i][4]) for i in range(64,69)]; KFeOH3Polynomial=[float(Reader1[i][5]) for i in range(64,69)] 
KNiOHPolynomial=[float(Reader1[i][6]) for i in range(64,69)]; KNiOH2Polynomial=[float(Reader1[i][7]) for i in range(64,69)]            
KNiOH3Polynomial=[float(Reader1[i][8]) for i in range(64,69)]                

class SectionParameters: #Defining each primary heat transport section as a class 
    def __init__(self,RowStart,RowEnd): 
        self.RowStart = RowStart; self.RowEnd = RowEnd
        #General section parameter input (interface specification not necessary)
        self.Diameter = [float(Reader1[i][1]) for i in range(self.RowStart,self.RowEnd)] #Reader([row][column]) #[cm]
        self.Velocity = [float(Reader1[i][2]) for i in range(self.RowStart,self.RowEnd)] #[cm/s]
        self.Length = [float(Reader1[i][3]) for i in range(self.RowStart,self.RowEnd)] #[cm]
        self.Celsius = [float(Reader1[i][4]) for i in range(self.RowStart, self.RowEnd)] #[oC]
        self.Kelvin = [x + 273.15 for x in self.Celsius] # [K]
        self.NernstConstant = [x*(2.303*R/(n*F)) for x in self.Kelvin] #2.303RT/nF
        
        self.Density = [float(Reader1[i][5]) for i in range(self.RowStart,self.RowEnd)] #[g/cm^3] = [kg/L]
        self.Viscosity = [float(Reader1[i][6]) for i in range(self.RowStart,self.RowEnd)] #[g/cm*s] 
        self.FeSolubility_Magnetite = [float(Reader1[i][7]) for i in range(self.RowStart,self.RowEnd)]#[mol/kg] (Tremaine & LeBlanc)
        self.FeSolubility_Chromite = [float(Reader1[i][8]) for i in range(self.RowStart,self.RowEnd)] #[mol/kg] chromite = FeCr2O4
        self.FeSolubility_Ferrite = [float(Reader1[i][9]) for i in range(self.RowStart,self.RowEnd)] #[mol/kg] (Sandler & Kunig)
        self.NiSolubility_Ferrite = [float(Reader1[i][10]) for i in range(self.RowStart,self.RowEnd)] #[mol/kg] (Sandler & Kunig)
        self.InnerOxideThickness = [float(Reader1[i][11]) for i in range(self.RowStart,self.RowEnd)] #[g/cm^2] 
        self.OuterMagnetiteThickness = [float(Reader1[i][12]) for i in range(self.RowStart,self.RowEnd)] #[g/cm^2]
        self.OuterFerriteThickness = [float(Reader1[i][13]) for i in range(self.RowStart,self.RowEnd)] #[g/cm^2] 
        self.OuterMetallicNickelThickness = [float(Reader1[i][14]) for i in range(self.RowStart,self.RowEnd)] #[g/cm^2] 
        self.OuterOxideThickness = [x+y+z for x,y,z in zip(self.OuterMagnetiteThickness,self.OuterFerriteThickness,\
                                    self.OuterMetallicNickelThickness)] #adds the lists together element by element
        
        #Electrochemical Parameters (Standard Potentials, interface specification not necessary)
        #Iron corrosion: Fe(s) <-> Fe2+(aq) + 2e-
        self.StandardEquilibriumPotentialFe = [float(Reader2[i+8][1]) for i in range(self.RowStart,self.RowEnd)] #[V] 
        #Magnetite precipitation: 3Fe(OH)2(s)<-> Fe3O4(s) + 2H2O(l) + 2H+(aq) + 2e-
        self.StandardEquilibriumPotentialFe3O4Precip =[float(Reader2[i+8][2]) for i in range(self.RowStart,self.RowEnd)] #[V] 
        #Magnetite dissolution: Fe3O4(s) + 2H2O(l) + 2H+(aq) + 2e- <-> 3Fe(OH)2(aq)
        self.StandardEquilibriumPotentialFe3O4Dissol = [float(Reader2[i+8][3]) for i in range(self.RowStart,self.RowEnd)] #[V] 
        #Nickel ferrite precipitation: 2.4Fe(OH)2(s) + 0.6Ni(OH)2(s) <-> Ni0.6Fe2.4O4(s) + 2H2O(l) + 2H+(aq) + 2e-
        self.StandardEquilibriumPotentialFerritePrecip = [float(Reader2[i+8][4]) for i in range(self.RowStart,self.RowEnd)]#[V] 
        #Nickel ferrite dissolution: Ni0.6Fe2.4O4(s) + 2H2O(l) + 2H+(aq) + 2e- <-> 2.4Fe(OH)2(aq) + 0.6Ni(OH)2(aq)
        self.StandardEquilibriumPotentialFerriteDissol = [float(Reader2[i+8][5]) for i in range(self.RowStart,self.RowEnd)] #[V] 
        #Hydrogen reaction: 2H+(aq) + 2e- <-> H2(g)
        self.StandardEquilibriumPotentialHydrogen = [float(Reader2[i+8][6]) for i in range(self.RowStart,self.RowEnd)] #[V] 
        #Nickel corrosion: Ni(s) <-> Ni2+(aq) + 2e-
        self.StandardEquilibriumPotentialNi = [float(Reader2[i+8][7]) for i in range(self.RowStart,self.RowEnd)]#[V] 
        #Nickel deposition: Ni(OH)2(aq)+ 2H+(aq) + 2e- <-> Ni(s) + 2H2O(l)
        self.StandardEquilibriumPotentialNiDeposition = [float(Reader2[i+8][8]) for i in range(self.RowStart,self.RowEnd)]#[V] 
#creates the 4 PHTS sections and their methods based on the Sections class template/blueprint
InletParameters=SectionParameters(4,11)
CoreParameters=SectionParameters(13,25)
OutletParameters=SectionParameters(27,36) 
SteamGeneratorParameters=SectionParameters(38,60)#ranges of rows for each sec. from txt file

#Equilibrium and Debye-Huckel constants - polynomials as a function of temperature 
def Constants(CelsiusTemp,KelvinTemp): #Coeff1*x^4 + Coeff2*x^3 + Coeff3*x^2 + Coeff4*x + Coeff5, depends on # elements in coeff list
    Constants.DebyeHuckConst=(np.polyval(DebyeHuckPolynomial,CelsiusTemp)) 
    Constants.Kw=10**(np.polyval(KwPolynomial,KelvinTemp))
    Constants.KLi=10**(np.polyval(KLiPolynomial,KelvinTemp)) 
    Constants.KFeOH=10**(np.polyval(KFeOHPolynomial,KelvinTemp))
    Constants.KFeOH2=10**(np.polyval(KFeOH2Polynomial,KelvinTemp))
    Constants.KFeOH3=10**(np.polyval(KFeOH3Polynomial,KelvinTemp))
    Constants.KNiOH=10**(np.polyval(KNiOHPolynomial,KelvinTemp))
    Constants.KNiOH2=10**(np.polyval(KNiOH2Polynomial,KelvinTemp))
    Constants.KNiOH3=10**(np.polyval(KNiOH3Polynomial,KelvinTemp))

#Numerical method to calculate the pH based on known ConcLi, eq'm constants(T), and an initial guess for ConcH and activity coeff.    
def NewtonRaphson(CelsiusTemp, KelvinTemp): #Bulk pH calculation
    Constants(CelsiusTemp, KelvinTemp)
    ConcH = 0.000001  #initial guess [mol/kg]
    gamma1 = 1 #initial guess 
    ConcOH =  Constants.Kw / ((gamma1**2) * (ConcH))
    #At high temp, LiOH doesn't dissociate 100% - eq'm established: LiOH(aq) <-> Li+(aq) + OH-(aq)
    #KLi = ([Li+]gamma1 *[OH-]gamma1)/[LiOH]; ConcLiTotal = ConcLi  + LiConcOH (sub in eq'm expression for LiConcOH
    ConcLi = Constants.KLi * ConcLiTotal / (Constants.KLi + ((gamma1**2) * ConcOH))
    #Kw = ConcH*gamma1*ConcOH*gamma1;   ConcOH = Kw/ConcH*gamma1**2
    #H+ + Li+ = OH-;                    ConcH = ConcLi + Kw/ConcH*gamma1**2
    for i in range(20): #no more than 10 iterations should really be necessary for convergence at the provided error  
        FH = ConcH**2 + (ConcH * ConcLi) - (Constants.Kw / (gamma1**2)) #function of H
        DFH = 2 * ConcH + ConcLi       #derivative of function
        NewConcH = ConcH - (FH / DFH)
        RE = abs((NewConcH - ConcH) / NewConcH)
        ConcH = NewConcH
        #print RE #optional value check of ConcH at each iteration 
        ConcOH = Constants.Kw / ((gamma1**2) / ConcH)
        ConcLi = (Constants.KLi * ConcLiTotal) / (Constants.KLi + ((gamma1**2) * ConcOH))
        IonicStrength = ((1**2) * ConcH + (1**2) * ConcLi + (1**2) * ConcOH) / 2
        #Davies equation loggamma1 = -DebyeHuckConst*(z**2)*[(sqrt(I)/(1+sqrt(I)))-beta*I]
        gamma1 = 10**(-Constants.DebyeHuckConst * (1**2) * (((IonicStrength ** 0.5) / (1 + (IonicStrength**0.5))) - 0.2 * IonicStrength))
        #All entries in ConcH list must meet error convergence requirement (no matter if SG 22 element list or Inlet 7 element list)
        if RE.all() < 0.000001: break #print i #prints number of iterations before error minimized to desired level and exits loop    
    NewtonRaphson.BulkConcH = ConcH                           
    NewtonRaphson.BulkConcOH = ConcOH
    
#---------------------------------- Can check pH at room temperature here 
# RTCelsius = 25 #oC
# RTKelvin = RTCelsius + 273 #Kelvin
# NewtonRaphson(RTCelsius,RTKelvin)
#----------------------------------
def Composition(CelsiusTemp, KelvinTemp, ConcH, ConcOH, gamma1, gamma2, FeTot, NiTot):
    Constants(CelsiusTemp, KelvinTemp)
    ConcOH = Constants.Kw / (ConcH * gamma1**2)
    Li = (ConcLiTotal * Constants.KLi) / (Constants.KLi + (ConcOH * gamma1**2))
    Composition.ConcFe2 = FeTot/(1 + (Constants.KFeOH * gamma2/(ConcH * gamma1**2)) + (Constants.KFeOH2 * gamma2 /((ConcH**2)*\
                        (gamma1**2))) + (Constants.KFeOH3 * gamma2/((ConcH**3) * (gamma1**4))))
    Composition.ConcFeOH = (Constants.KFeOH * gamma2 * Composition.ConcFe2) / (ConcH * gamma1**2)
    Composition.ConcFeOH2 = (Constants.KFeOH2 * gamma2 * Composition.ConcFe2) / ((ConcH**2) * gamma1**2)
    Composition.ConcFeOH3 = (Constants.KFeOH3 * gamma2 * Composition.ConcFe2) / ((ConcH**3) * gamma1**4)

    Composition.ConcNi2 = NiTot / (1 + (Constants.KNiOH * gamma2 / (ConcH * gamma1**2)) + (Constants.KNiOH2 * gamma2 / ((ConcH**2) * \
                        (gamma1**2))) + (Constants.KNiOH3 * gamma2 / ((ConcH**3) * (gamma1**4))))
    Composition.ConcNiOH = (Constants.KNiOH * gamma2 * Composition.ConcNi2) / (ConcH * gamma1**2)
    Composition.ConcNiOH2 = (Constants.KNiOH2 * gamma2 * Composition.ConcNi2) / ((ConcH**2) * gamma1**2)
    Composition.ConcNiOH3 = (Constants.KNiOH3 * gamma2 * Composition.ConcNi2) / ((ConcH**3) * gamma1**4)
    #FeOH2 and NiOH2 left out of ionic strength calc due to electroneutrality principle
    IonicStrength = ((1**2) * ConcH + ((-1)**2) * ConcOH + (2**2) * Composition.ConcFe2 + (1**2) * Composition.ConcFeOH + \
                     ((-1)**2) * Composition.ConcFeOH3 + (2**2) * Composition.ConcNi2 + (1**2) * Composition.ConcNiOH + \
                     ((-1)**2) * Composition.ConcNiOH3 + (1**2) * Li) / 2
    gamma1 = 10**(-Constants.DebyeHuckConst * (1**2) * (((IonicStrength**0.5) / (1 + (IonicStrength**0.5))) - 0.2 * IonicStrength))
    gamma2 = 10**(-Constants.DebyeHuckConst * ((-2)**2) * (((IonicStrength**0.5) / (1 + (IonicStrength**0.5))) - 0.2 * IonicStrength))    
    Composition.gamma1 = gamma1
    Composition.gamma2 = gamma2 
    Composition.ConcH = ConcH
    Composition.ConcOH = ConcOH
                
def SolutionOxide(SectionParameters):
    #Initial parameters based on bulk pH and oxide solubilities 
    #getattr(SectionParameters, "FeSolubility_Magnetite") self.NiTot = getattr(SectionParameters, "NiSolubility_Ferrite")
    gamma1 = 1 #initial estimate for activities (+/- 1 charged ions)
    gamma2 = 1 #(+/- 2 charged ions)   
    if SectionParameters == OutletParameters:
        #takes first element of inlet solubility list and converts to desired length (9 element for outlet; 12 element for core)
        SolutionOxide.FeTot = [InletParameters.FeSolubility_Magnetite[0]]*9 #Inlet total iron set up as < FeSat at outlet and core
        SolutionOxide.NiTot = [InletParameters.NiSolubility_Ferrite[0]]*9 #Inlet total nickel set up as < FeSat at outlet/core 
    elif SectionParameters == CoreParameters:
        SolutionOxide.FeTot = [InletParameters.FeSolubility_Magnetite[0]]*12 #Inlet total iron set up as < FeSat at outlet and core
        SolutionOxide.NiTot = [InletParameters.NiSolubility_Ferrite[0]]*12 #Inlet total nickel set up as < FeSat at outlet/core
    else: #SectionParameters in (InletParameters, SteamGeneratorParameters):  
        SolutionOxide.FeTot = SectionParameters.FeSolubility_Magnetite
        SolutionOxide.NiTot = SectionParameters.NiSolubility_Ferrite
    SolutionOxide.FeSat = SectionParameters.FeSolubility_Magnetite
    SolutionOxide.NiSat =SectionParameters.NiSolubility_Ferrite
    
    #values assigned to Fe and Ni Tot based on solutilities of dominant oxides at SO interface 
    Composition(SectionParameters.Celsius,SectionParameters.Kelvin,NewtonRaphson.BulkConcH,NewtonRaphson.BulkConcOH,gamma1, \
                gamma2, SolutionOxide.FeTot, SolutionOxide.NiTot)
    ConcFeOH2 = Composition.ConcFeOH2 #[mol/kg]
    ConcNiOH2 = Composition.ConcNiOH2 #[mol/kg]
    gamma1 = Composition.gamma1
    gamma2 = Composition.gamma2
    SolutionOxide.ConcH = Composition.ConcH #[mol/kg]    
    #Saturated feeders: 3Fe(OH)2(s) -> Fe3O4(s) + 2H+ +2H2O + 2e-
    if SolutionOxide.FeTot >= SolutionOxide.FeSat: # == comparison operator 
            #Define Q as (activity products/activity reactants)
            QFe3O4 = [(x*y*z)**2 for x,y,z in zip(SolutionOxide.ConcH,gamma1,SectionParameters.Density)]
            EquilibriumPotentialFe3O4 = [x+(y*z) for x,y,z in zip(SectionParameters.StandardEquilibriumPotentialFe3O4Precip, \
                                        SectionParameters.NernstConstant,np.log10(QFe3O4))]    
    #Undersaturated feeders        
    else:#if self.FeTot < self.FeSat:
            QFe3O4 = [((x*y)**3)/(z*y*q)**2 for x,y,z,q in zip(ConcFeOH2, SectionParameters.Density, SolutionOxide.ConcH, gamma1)]  
            EquilibriumPotentialFe3O4 = [x-(y*z) for x,y,z in zip(SectionParameters.StandardEquilibriumPotentialFe3O4Dissol, \
                                        SectionParameters.NernstConstant, np.log10(QFe3O4))]
    SolutionOxide.MixedPotential = EquilibriumPotentialFe3O4 #initial guess for concH @ MO interfacec adjustment
    
def MetalOxide(SectionParameters):
    gamma1 = 1 #initial estimate for activities (+/- 1 charged ions)
    gamma2 = 1 #(+/- 2 charged ions)
    MetalOxide.FeTot = SectionParameters.FeSolubility_Magnetite
    MetalOxide.NiTot = SectionParameters.NiSolubility_Ferrite
#     MetalOxide.FeSat = SectionParameters.FeSolubility_Magnetite
#     MetalOxide.NiSat =SectionParameters.NiSolubility_Ferrite 
    Composition(SectionParameters.Celsius,SectionParameters.Kelvin,NewtonRaphson.BulkConcH,NewtonRaphson.BulkConcOH,gamma1,\
                gamma2, MetalOxide.FeTot, MetalOxide.NiTot)
    ConcFe2 = Composition.ConcFe2 #[mol/kg]
    ConcNi2 = Composition.ConcNi2  #[mol/kg]
    gamma1 = Composition.gamma1
    gamma2 = Composition.gamma2
    MetalOxide.ConcH = Composition.ConcH#[mol/kg]
    
    QFe = [x*y*z for x,y,z in zip(ConcFe2,gamma2,SectionParameters.Density)] # activity products/activity reactants
    #Eeqm = Eo + (2.303RT/nF)*log(Q), [V]
    MetalOxide.EquilibriumPotentialFe =[x+(y*z) for x,y,z in zip(SectionParameters.StandardEquilibriumPotentialFe,\
                            SectionParameters.NernstConstant,np.log10(QFe))]
    MetalOxide.MixedPotential = MetalOxide.EquilibriumPotentialFe 

class ActivationEnergy: 
    def __init__(self, SectionParameters):
        NewtonRaphson(SectionParameters.Celsius,SectionParameters.Kelvin)
        SolutionOxide(SectionParameters)
        MetalOxide(SectionParameters) 
        #Effect of mixed potential on the ConcH at the MO interface 
        ElectrochemicalEffectConstant = [(1/x) * (F/R) for x in SectionParameters.Kelvin] #F/RT
        #no metal-oxide interface in the core 
        if SectionParameters in (InletParameters, OutletParameters, SteamGeneratorParameters): #if either of these sections is called
            DeltaPotential=[(x-y) for x,y in zip(MetalOxide.MixedPotential,SolutionOxide.MixedPotential)]
            #MetalOxide.ConcH =[x*y for x,y in zip(SolutionOxide.ConcH, np.exp([x*y for x,y in zip(DeltaPotential,ElectrochemicalEffectConstant)]))]
            MetalOxide.ConcH = [x*y for x,y in zip(SolutionOxide.ConcH,[1-x for x in [x*y for x,y in zip(DeltaPotential,ElectrochemicalEffectConstant)]])]
            #Initial estimates for activation energy calculation  
            ioFe = 4.467*10**-6 #A/cm^2 adjustable 
            MetalOxide.MixedPotential = -0.953 #V adjustable 
            ButlerVolmerConstant=[(n*Beta)*x for x in ElectrochemicalEffectConstant] #Beta*n*F/RT 
            print MetalOxide.EquilibriumPotentialFe
            ActivationEnergyFe 
            
InletEnergies=ActivationEnergy(InletParameters) 
OutletEnergies=ActivationEnergy(OutletParameters)
CoreEnergies = ActivationEnergy(CoreParameters)
SteamGeneratorEnergies=ActivationEnergy(SteamGeneratorParameters) 



# def Corrosion():
#     
# def MixedPotential():




# def main():
#     x = 3
#     f(x)
# def f(x):
#     print(x)  
# main()


