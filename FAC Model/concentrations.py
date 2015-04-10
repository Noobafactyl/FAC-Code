#from scipy.optimize import fsolve
import numpy as np
import constant_values as cv
import parameter_input as pi
import equilibrium_potential as ep
import corrosion_rate as cr
'recovery save'
class Concentrations(cr.CorrosionRate): #if use MetalOxide, this will already include SolutionOxide (don't need to inherit both) 
    def __init__(self,Section, ConcentrationFeTotal, SaturationFe, gamma_1, gamma_2, ConcentrationNiTotal, SaturationNi, \
                 ConcentrationFeOH2,ConcentrationNiOH2, ConcentrationH, ConcentrationH2):
        ae.SolutionOxide.__init__(self,Section) #allows access to activation energies and equilibrium potentials, as calculated via
        #the initial concentration estimates of the concentrations (based on solubilities)
        #change the concentrations of aq species here 
        
        #if the self.concentrations change when called here from another class, do the ep.EquilibriumPotential values automatically
        #change or are these carried through from the original implementation of this class in SolutionOxide??
        
        
        ep.EquilibriumPotential.__init__(self, Section, self.ConcentrationFeTotal, self.SaturationFe, self.gamma_1, self.gamma_2,\
         self.ConcentrationNiTotal, self.SaturationNi, self.ConcentrationFeOH2, self.ConcentrationNiOH2, self.ConcentrationH,\
         self.ConcentrationH2)

Inlet = Concentrations(pi.InletParameters)
print Inlet.EquilibriumPotential_Ni