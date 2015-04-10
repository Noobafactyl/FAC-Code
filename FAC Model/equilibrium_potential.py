import numpy as np
import constant_values as cv
import parameter_input as pi
'test' 
class EquilibriumPotential():
        def __init__(self,Section, ConcentrationFeTotal, SaturationFe, gamma_1, gamma_2, ConcentrationNiTotal, SaturationNi, ConcentrationFeOH2,\
                          ConcentrationNiOH2, ConcentrationH, ConcentrationH2, ConcentrationFe2, ConcentrationNi2):
        
        #saturated feeders (inlet, part of steam generators)
            if ConcentrationFeTotal >= SaturationFe:
        #Fe(OH)2(s) ->Fe3O4(s) + 2e- + 2H+(aq) + 2H2O(l), oxidative precipitation
                Q_Fe3O4 = [(x*y*z)**2 for x,y,z in zip(ConcentrationH, gamma_1, Section.Density)]
                self.EquilibriumPotential_Fe3O4 = [x+(y*z) for x,y,z in zip(Section.StandardEquilibriumPotential_Fe3O4Precipitation, \
                                         Section.NernstConstant,np.log10(Q_Fe3O4))]
        #undersaturated feeders (outlet, part of core)        
            else:#if ConcentrationFeTotal < SaturationFe:
        #Fe3O4(s) + 2e- + 2H+ + 2H2O(l) -> 3Fe(OH)2(aq), reductive dissolution
                Q_Fe3O4 = [((x*y)**3)/(z*y*q)**2 for x,y,z,q in zip(ConcentrationFeOH2, Section.Density, ConcentrationH, gamma_1)]  
                self.EquilibriumPotential_Fe3O4 = [x-(y*z) for x,y,z in zip(Section.StandardEquilibriumPotential_Fe3O4Dissolution, \
                                         Section.NernstConstant, np.log10(Q_Fe3O4))]
        #saturated feeders (inlet, part of steam generators)
        #0.6Ni(OH)2(s) + 2.4Fe(OH)2(s) -> 2H+(aq) +2e- + 2H2O(l) + Ni0.6Fe2O4(s), oxidative precipitation    
            if ConcentrationNiTotal >= SaturationNi:
                Q_Ferrite = [(x*y*z)**2 for x,y,z in zip(ConcentrationH, gamma_1, Section.Density)]
                self.EquilibriumPotential_Ferrite = [x+(y*z) for x,y,z in zip(Section.StandardEquilibriumPotential_FerritePrecipitation, \
                                         Section.NernstConstant,np.log10(Q_Ferrite))]
        #undersaturated feeders
            else:#if ConcentrationNiTotal < SaturationNi:
        #Ni0.6Fe2O4(s) + 2H+(aq)+ 2e- +2H2O(l) -> 0.6Ni(OH)2(aq)+ 2.4Fe(OH)2(aq), reductive dissolution
                Q_Ferrite = [((x*y)**2.4)*((r*y)**0.6)/(z*y*q)**2 for x,y,r,z,q in zip(ConcentrationFeOH2, Section.Density, \
                        ConcentrationNiOH2, ConcentrationH, gamma_1)] 
                self.EquilibriumPotential_Ferrite = [x-(y*z) for x,y,z in zip(Section.StandardEquilibriumPotential_FerriteDissolution, \
                                         Section.NernstConstant, np.log10(Q_Ferrite))]
        #2H+(aq) + 2e- -> H2(g) 
            Q_H2 = [x/((y*z*q)**2) for x,y,z,q in zip([ConcentrationH2*i/cv.kH2 for i in Section.Density], \
                    ConcentrationH, gamma_1, Section.Density)]
            self.EquilibriumPotential_H2 = [x-(y*z) for x,y,z in zip(Section.StandardEquilibriumPotential_H2, \
                                         Section.NernstConstant,np.log10(Q_H2))]
        #Fe(s) -> Fe2+(aq) + 2e-
            Q_Fe = [1/(x*y*z) for x,y,z in zip(ConcentrationFe2, Section.Density, gamma_2)]
            self.EquilibriumPotential_Fe = [x-(y*z) for x,y,z in zip(Section.StandardEquilibriumPotential_Fe, \
                                         Section.NernstConstant,np.log10(Q_Fe))]
        #Ni(s) -> Ni2+(aq) + 2e-    
            Q_Ni = [1/(x*y*z) for x,y,z in zip(ConcentrationNi2, Section.Density, gamma_2)]
            self.EquilibriumPotential_Ni = [x-(y*z) for x,y,z in zip(Section.StandardEquilibriumPotential_Ni, \
                                         Section.NernstConstant,np.log10(Q_Ni))]
        #Ni(OH)02(aq)+ 2H+(aq) + 2e- -> Ni(s) + 2H2O(l)    
        #return {'E_Fe3O4':EquilibriumPotential_Fe3O4[1], 'E_Ferrite': EquilibriumPotential_Ferrite[1], \
         #       'E_H2':EquilibriumPotential_H2[1], 'E_Ni':EquilibriumPotential_Ni[1]}
    