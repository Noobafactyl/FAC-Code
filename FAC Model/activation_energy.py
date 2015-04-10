
'Notes: oxidation and reduction reactions are left as is without re-writing all rxns in one convention. This is to prevent \
confusion about the reaction mechanism (e.g., reductive dissolution/oxidative precipitation). The Nernst equation is thus \
re-cast to accommodate the reduction or oxidation, i.e., the sign in the equation changes accordingly: \
E = Eo - (2.303*R*T/nF)*log(Q) versus E = Eo + (2.303*R*T/nF)*log(Q), for reduction and oxidation, respectively.\
Define Q as activity products/activity reactants.'

'Unit conversions: all starting concentrations are in mol/kg, which are converted to mol/L using the density of water, kg/L:\
[mol/kg]*[kg/L] = [mol/L] (need concentration in M to use in Nernst equation)\
For hydrogen gas (H2), the mol/kg concentration is converted from mol/L to atm: [mol/kg]*[kg/L] /[mol/L*atm] = [atm]'

'In ActivationEnergy function: A^z+ + ze- -> D,   where A = e- acceptor and D = e- donor (Bockris & Reddy - Modern Electrochemistry)\
If in terms of acceptor:  io = (Concentration_A*FkT/h)*exp(-ActivationEnergy/RT)*exp(-Eeqm*BnF/RT), \
If in terms of donor: io =(Concentration_D*FkT/h)*exp(-ActivationEnergy/RT)*exp(Eeqm*(1-B)nF/RT),  B term is positive in second exponential\
Concentration donor/acceptor is [mol/cm^2] (raise to the 2/3 bc of stoich)'

'recovery save'
import numpy as np
import constant_values as cv
import parameter_input as pi
import equilibrium_potential as ep
import composition as cmp
 
class SolutionOxide(cmp.Composition, ep.EquilibriumPotential):
    def __init__(self, Section):
         cmp.Composition.__init__(self,Section)
        
         #inlet/steam generators: ConcentrationFe set = SaturationFe to trigger precipitation 
         #outlet/core: ConcentrationFe set = InletSaturationFe, ensuring that ConcentrationFe < SaturationFe, to trigger dissolution 
         if Section == pi.OutletParameters: 
            self.ConcentrationFeTotal = [pi.InletParameters.SolubilityMagnetite_Fe[0]]*9 
            self.ConcentrationNiTotal = [pi.InletParameters.SolubilityFerrite_Ni[0]]*9 
         elif Section == pi.CoreParameters:
            self.ConcentrationFeTotal = [pi.InletParameters.SolubilityMagnetite_Fe[0]]*12
            self.ConcentrationNiTotal = [pi.InletParameters.SolubilityFerrite_Ni[0]]*12 
         else:
            self.ConcentrationFeTotal = Section.SolubilityMagnetite_Fe
            self.ConcentrationNiTotal = Section.SolubilityFerrite_Ni
         self.SaturationFe = Section.SolubilityMagnetite_Fe 
         self.SaturationNi = Section.SolubilityFerrite_Ni #initial estimate based on non-stoichiometric nickel ferrite solubility (S&K)
         #convert ConcentrationH2 from cm^3/kg to mol/kg, (cm^3/kg *  g/cm^3)/(g/mol) = [mol/kg]
         self.ConcentrationH2 = cv.ConcentrationH2Total*cv.H2Density/cv.H2MolarMass
         self.ConcentrationH = self.BulkProtonConcentrationentration()
         
         self.ConcentrationIteration(Section, self.ConcentrationFeTotal, self.ConcentrationNiTotal, self.SaturationFe, self.SaturationNi, self.ConcentrationH)
         
         ep.EquilibriumPotential.__init__(self,Section, self.ConcentrationFeTotal, self.SaturationFe, self.gamma_1, self.gamma_2,\
         self.ConcentrationNiTotal, self.SaturationNi, self.ConcentrationFeOH2, self.ConcentrationNiOH2, self.ConcentrationH, \
         self.ConcentrationH2, self.ConcentrationFe2, self.ConcentrationNi2)

class MetalOxide(cmp.Composition, ep.EquilibriumPotential):
    def __init__(self, Section):
         cmp.Composition.__init__(self,Section)
         
         #inlet/steam generators: ConcentrationFe set = SaturationFe to trigger precipitation 
         #outlet/core: ConcentrationFe set = InletSaturationFe, ensuring that ConcentrationFe < SaturationFe, to trigger dissolution 
         self.ConcentrationFeTotal = Section.SolubilityMagnetite_Fe
         self.ConcentrationNiTotal = Section.SolubilityFerrite_Ni
         self.SaturationFe = Section.SolubilityMagnetite_Fe 
         self.SaturationNi = Section.SolubilityFerrite_Ni #initial estimate based on non-stoichiometric nickel ferrite solubility (S&K)
         #convert ConcentrationH2 from cm^3/kg to mol/kg, (cm^3/kg *  g/cm^3)/(g/mol) = [mol/kg]
         self.ConcentrationH2 = cv.ConcentrationH2Total*cv.H2Density/cv.H2MolarMass
         self.ConcentrationH = self.BulkProtonConcentrationentration()
         
         self.ConcentrationIteration(Section, self.ConcentrationFeTotal, self.ConcentrationNiTotal, self.SaturationFe, self.SaturationNi, self.ConcentrationH)
         
         ep.EquilibriumPotential.__init__(self, Section, self.ConcentrationFeTotal, self.SaturationFe, self.gamma_1, self.gamma_2,\
         self.ConcentrationNiTotal, self.SaturationNi, self.ConcentrationFeOH2, self.ConcentrationNiOH2, self.ConcentrationH, \
         self.ConcentrationH2, self.ConcentrationFe2, self.ConcentrationNi2)
    
         
class ActivationEnergy(SolutionOxide, MetalOxide): #make a function of the interface (MO vs SO)
    def __init__(self, Section):
        SolutionOxide.__init__(self, Section)
        
        ButlerVolmerConstant = [(1/x) * (cv.Beta*cv.n*cv.F/cv.R) for x in Section.Kelvin] #BnF/RT       
         
        if self.ConcentrationFeTotal >= self.SaturationFe:
            MixedPotential_Fe3O4 = -0.927 #V initial estimate 
            ExchangeCurrent_Fe3O4 = 1E-7 #A/cm^2 initial estimate
        else:
            MixedPotential_Fe3O4 = -0.79 #V initial estimate  
            ExchangeCurrent_Fe3O4 = 1E-9 #A/cm^2 initial estimate
        
        if self.ConcentrationNiTotal >= self.SaturationNi:
            MixedPotential_Ferrite = -0.915 #V initial estimate 
            ExchangeCurrent_Ferrite = 1E-10 #A/cm^2 initial estimate
        else:
            MixedPotential_Ferrite = -0.59 #V initial estimate  
            ExchangeCurrent_Ferrite = 1E-9 #A/cm^2 initial estimate
         
        MixedPotential_Fe = -0.95 #V
        ExchangeCurrent_Fe = 4.467E-6 #A/cm^2
         
        MixedPotential_Ni = -0.88 #V
        ExchangeCurrent_Ni = 6.25E-6 #A/cm^2  
         
         
        def HydrogenExchangeCurrent(self, Section, EquilibriumPotentialAnode, EquilibriumPotentialCathode, ExchangeCurrent, MixedPotential):
            #Doesn't matter which is the anode and which is the cathode for the expression below - algebra works out the same
            DeltaAnodePotential = [MixedPotential - x for x in EquilibriumPotentialAnode]
            DeltaCathodePotential = [MixedPotential - x for x in EquilibriumPotentialCathode]
             
            I = [(x-y)/(z-q) for x,y,z,q in zip(np.exp([x*y for x,y in zip(ButlerVolmerConstant,DeltaAnodePotential)]),\
            np.exp([x*(-y) for x,y in zip(ButlerVolmerConstant,DeltaAnodePotential)]),\
            np.exp([x*(-y) for x,y in zip(ButlerVolmerConstant,DeltaCathodePotential)]), \
            np.exp([x*y for x,y in zip(ButlerVolmerConstant,DeltaCathodePotential)]))]
            return [ExchangeCurrent * x for x in I]
             
        ExchangeCurrent_H2onFe3O4 = HydrogenExchangeCurrent(self, Section, self.EquilibriumPotential_Fe3O4, self.EquilibriumPotential_H2, ExchangeCurrent_Fe3O4, MixedPotential_Fe3O4)
        
        ExchangeCurrent_H2onFerrite = HydrogenExchangeCurrent(self, Section, self.EquilibriumPotential_Ferrite, self.EquilibriumPotential_H2, ExchangeCurrent_Ferrite, MixedPotential_Ferrite)
        ExchangeCurrent_H2onFe = HydrogenExchangeCurrent(self, Section, self.EquilibriumPotential_Fe, self.EquilibriumPotential_H2, ExchangeCurrent_Fe, MixedPotential_Fe)
        ExchangeCurrent_H2onNi = HydrogenExchangeCurrent(self, Section, self.EquilibriumPotential_Ni, self.EquilibriumPotential_H2, ExchangeCurrent_Ni, MixedPotential_Ni)
        #print 'i_Ni:', ExchangeCurrent_Ni,'A/cm^2'
        
        
        def OxideExchangeCurrent(self, Section, Concentration, EquilibriumPotential, ExchangeCurrent):
            if Concentration == None:
                A =[(cv.F*cv.kb*x) for x in Section.Kelvin]
            else:
                A =[x*y for x,y in zip([(cv.F*cv.kb*x) for x in Section.Kelvin], \
                                     [(x/1000)**(2.0/3) for x in [x*y for x,y in zip(Concentration,Section.Density)]])]
            B = np.exp([x*y for x,y in zip(ButlerVolmerConstant,EquilibriumPotential)])  
            C = [x*y for x,y in zip(A,B)]
            
            if type(ExchangeCurrent)== list: 
                return [x*y for x,y in zip([cv.R*x for x in Section.Kelvin], \
                    -np.log([x*y for x,y in zip([cv.hp/x for x in C], ExchangeCurrent)]))]
            else:
                return [x*y for x,y in zip([cv.R*x for x in Section.Kelvin], -np.log([ExchangeCurrent*cv.hp/x for x in C]))]
         
        if self.ConcentrationFeTotal >= self.SaturationFe:
            self.ActivationEnergy_Fe3O4 = OxideExchangeCurrent(self,Section, None, self.EquilibriumPotential_Fe3O4,ExchangeCurrent_Fe3O4)
            self.ActivationEnergy_H2onFe3O4 = OxideExchangeCurrent(self,Section, self.ConcentrationH, [-1*x for x in self.EquilibriumPotential_H2],ExchangeCurrent_H2onFe3O4)
        else:
            self.ActivationEnergy_Fe3O4 = OxideExchangeCurrent(self,Section,self.ConcentrationFeOH2,self.EquilibriumPotential_Fe3O4,ExchangeCurrent_Fe3O4)
            self.ActivationEnergy_H2onFe3O4 = OxideExchangeCurrent(self,Section, self.ConcentrationH, [-1*x for x in self.EquilibriumPotential_H2] ,ExchangeCurrent_H2onFe3O4)
         
        if self.ConcentrationNiTotal >= self.SaturationNi:
            self.ActivationEnergy_Ferrite = OxideExchangeCurrent(self,Section, self.ConcentrationH, [-1*x for x in self.EquilibriumPotential_Ferrite],ExchangeCurrent_Ferrite)
            self.ActivationEnergy_H2onFeFerrite = OxideExchangeCurrent(self,Section, self.ConcentrationH, [-1*x for x in self.EquilibriumPotential_H2],ExchangeCurrent_H2onFerrite)
        else:
            self.ActivationEnergy_Ferrite = OxideExchangeCurrent(self,Section, self.ConcentrationH, [-1*x for x in self.EquilibriumPotential_Ferrite],ExchangeCurrent_Ferrite)
            self.ActivationEnergy_H2onFeFerrite = OxideExchangeCurrent(self,Section, self.ConcentrationH, [-1*x for x in self.EquilibriumPotential_H2],ExchangeCurrent_H2onFerrite)
        
        MetalOxide.__init__(self, Section)
        #print self.ConcentrationFe2[1]
        self.ActivationEnergy_Fe = OxideExchangeCurrent(self,Section,self.ConcentrationFe2, [-1*x for x in self.EquilibriumPotential_Fe], ExchangeCurrent_Fe)
        self.ActivationEnergy_H2onFe = OxideExchangeCurrent(self,Section, self.ConcentrationH, [-1*x for x in self.EquilibriumPotential_H2],ExchangeCurrent_H2onFe)
        
        if Section == pi.SteamGeneratorParameters:
            self.ActivationEnergy_Ni = OxideExchangeCurrent(self,Section,self.ConcentrationNi2, [-1*x for x in self.EquilibriumPotential_Ni], ExchangeCurrent_Ni)
            #self.ActivationEnergy_H2onNi = OxideExchangeCurrent(self,Section, self.ConcentrationH, [-1*x for x in self.EquilibriumPotential_H2],ExchangeCurrent_H2onNi)
        else:
            self.ActivationEnergy_Ni = 0
            self.ActivationEnergy_H2onNi = 0
        
        
InletAE = ActivationEnergy(pi.InletParameters)
OutletAE = ActivationEnergy(pi.OutletParameters)
CoreAE = ActivationEnergy(pi.CoreParameters)
SGAE = ActivationEnergy(pi.SteamGeneratorParameters)

print 'Inlet:', InletAE.ActivationEnergy_Ni, 'J/mol'
print 'SG:', SGAE.ActivationEnergy_Ni[1], 'J/mol'


#     MetalOxide.ConcentrationH = Composition.ConcentrationH#[mol/kg]
#     
#     QFe = [x*y*z for x,y,z in zip(ConcentrationFe2,gamma_2,SectionParameters.Density)] # activity products/activity reactants
#     #Eeqm = Eo + (2.303RT/nF)*log(Q), [V]
#     MetalOxide.EquilibriumPotentialFe =[x+(y*z) for x,y,z in zip(SectionParameters.StandardEquilibriumPotentialFe,\
#                             SectionParameters.NernstConstant,np.log10(QFe))]
#     MetalOxide.MixedPotential = MetalOxide.EquilibriumPotentialFe 
#


