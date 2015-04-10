import pH_calculator as p
import constant_values as cv
'recovery save'
class Composition(p.BulkpHCalculator):
    def __init__(self, Section):
        p.BulkpHCalculator.__init__(self, Section)
        
    def SpeciesConcentrations(self, ConcentrationFeTotal, ConcentrationNiTotal, ConcentrationH, gamma_1, gamma_2):        
        #ConcentrationFeTotal, ConcentrationNiTotal, and ConcentrationH are not re-evaluated here, thus self. not used for them 
        self.ConcentrationOH = self.k_W / (ConcentrationH * self.gamma_1**2)
        self.ConcentrationLi = (cv.ConcentrationLiTotal * self.k_Li) / (self.k_Li + (self.ConcentrationOH * self.gamma_1**2))
        self.ConcentrationFe2 = ConcentrationFeTotal/(1 + (self.k_FeOH * self.gamma_2/(ConcentrationH * self.gamma_1**2)) + (self.k_FeOH2 * self.gamma_2 /
                        ((ConcentrationH**2)*(self.gamma_1**2)))+ (self.k_FeOH3 * self.gamma_2/((ConcentrationH**3) * (self.gamma_1**4))))
        self.ConcentrationFeOH = (self.k_FeOH * self.gamma_2 * self.ConcentrationFe2) / (ConcentrationH * self.gamma_1**2)
        self.ConcentrationFeOH2 = (self.k_FeOH2 * self.gamma_2 * self.ConcentrationFe2) / ((ConcentrationH**2) * self.gamma_1**2)
        self.ConcentrationFeOH3 = (self.k_FeOH3 * self.gamma_2 * self.ConcentrationFe2) / ((ConcentrationH**3) * self.gamma_1**4)
 
        self.ConcentrationNi2 = ConcentrationNiTotal / (1 + (self.k_NiOH * self.gamma_2 / (ConcentrationH * self.gamma_1**2)) + (self.k_NiOH2 * self.gamma_2 / 
                        ((ConcentrationH**2) *(self.gamma_1**2))) + (self.k_NiOH3 * self.gamma_2 / ((ConcentrationH**3) * (self.gamma_1**4))))
        self.ConcentrationNiOH = (self.k_NiOH * self.gamma_2 * self.ConcentrationNi2) / (ConcentrationH * self.gamma_1**2)
        self.ConcentrationNiOH2 = (self.k_NiOH2 * self.gamma_2 * self.ConcentrationNi2) / ((ConcentrationH**2) * self.gamma_1**2)
        self.ConcentrationNiOH3 = (self.k_NiOH3 * self.gamma_2 * self.ConcentrationNi2) / ((ConcentrationH**3) * self.gamma_1**4)
        #FeOH2 and NiOH2 left out of ionic strength calc due to electroneutrality principle
        IonicStrength = ((1**2) * ConcentrationH + ((-1)**2) * self.ConcentrationOH + (2**2) * self.ConcentrationFe2 + (1**2) * self.ConcentrationFeOH +((-1)**2) *
                          self.ConcentrationFeOH3 + (2**2) * self.ConcentrationNi2 + (1**2) * self.ConcentrationNiOH +((-1)**2) * self.ConcentrationNiOH3 + (1**2) * self.ConcentrationLi) / 2
        self.gamma_1 = 10**(-self.DebyeHuckelConstant * (1**2) * (((IonicStrength**0.5) / (1 + (IonicStrength**0.5))) - 0.2 * IonicStrength))
        self.gamma_2 = 10**(-self.DebyeHuckelConstant * ((-2)**2) * (((IonicStrength**0.5) / (1 + (IonicStrength**0.5))) - 0.2 * IonicStrength))    
            #print self.gamma_1, self.ConcentrationNi2
    
    
    def ConcentrationIteration(self, Section, ConcentrationFeTotal, ConcentrationNiTotal, SaturationFe, SaturationNi, ConcentrationH):
         self.gamma_1 = 1.0 #initial estimate for activities (+/- 1 charged ions)
         self.gamma_2 =1.0 #(+/- 2 charged ions)        
         for i in range(50):
           gamma_1itr = self.gamma_1
           gamma_2itr = self.gamma_2
           self.SpeciesConcentrations(ConcentrationFeTotal, ConcentrationNiTotal, ConcentrationH, self.gamma_1, self.gamma_2)
           RE = ((gamma_1itr - self.gamma_1) / gamma_1itr)
           RE2 = ((gamma_2itr - self.gamma_2) / gamma_2itr) 
           if RE.all() < 0.0000001 and RE2.all()< 0.0000001: break
           
