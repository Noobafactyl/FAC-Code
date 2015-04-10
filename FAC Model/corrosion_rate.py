import numpy as np
import parameter_input as pi
import constant_values as cv
import equilibrium_potential as ep
import activation_energy as ae

'recovery save'
class CorrosionRate(ae.ActivationEnergy):
    def __init__(self,Section, ConcentrationFeTotal, ConcentrationNiTotal, SaturationFe, SaturationNi, ConcentrationH):
        ae.ActivationEnergy.__init__(self,Section)
        
        #self.ConcentrationFeTotal = ConcentrationFeTotal #if want to update output of these concentrations from this function/class
        #then use self inside ConcentrationIteration function below
        self.ConcentrationIteration(Section, ConcentrationFeTotal, ConcentrationNiTotal, SaturationFe, SaturationNi, ConcentrationH)
#         print self.gamma_1[1]
#         print ConcentrationFeTotal
        if ConcentrationNiTotal == 0:
            None
        
        ep.EquilibriumPotential.__init__(self,Section, ConcentrationFeTotal, SaturationFe, self.gamma_1, self.gamma_2,\
        ConcentrationNiTotal, SaturationNi, self.ConcentrationFeOH2, self.ConcentrationNiOH2, self.ConcentrationH, \
        self.ConcentrationH2, self.ConcentrationFe2, self.ConcentrationNi2)
        
        print self.EquilibriumPotential_Fe[1]
       

InletCR = CorrosionRate(pi.InletParameters,1E-8,1E-9,1E-8,1E-9,1E-8)
#print InletCR.ConcentrationFeTotal