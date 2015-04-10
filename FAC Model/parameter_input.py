'recovery save'
import csv
import constant_values as cv

GeneralParameters = open('GeneralParameters.txt','r') #Don't need to specify file path if docs are in the same folder as this .py file     
Reader1 = list(csv.reader(GeneralParameters, delimiter = ',')) #Assign file data to list, Reader[row][column], delimiter = comma 

ElectrochemicalParameters = open('ElectrochemicalParameters.txt','r')
Reader2 = list(csv.reader(ElectrochemicalParameters, delimiter = ','))

#class = data structure -holds data and methods to process data (i.e., functions)variables inside class=attributes; functions=methods
class SectionParameters(object): #Defining each primary heat transport section as a class 
    def __init__(self,RowStart=None,RowEnd=None): 
        self.RowEnd = RowEnd
        self.RowStart = RowStart
        
        self.Diameter = [float(Reader1[i][1]) for i in range(self.RowStart,self.RowEnd)] #Reader([row][column]) #[cm]
        self.Velocity = [float(Reader1[i][2]) for i in range(self.RowStart,self.RowEnd)] #[cm/s]
        self.Length = [float(Reader1[i][3]) for i in range(self.RowStart,self.RowEnd)] #[cm]
        self.Celsius = [float(Reader1[i][4]) for i in range(self.RowStart, self.RowEnd)] #[oC]
        self.Kelvin = [x + 273.15 for x in self.Celsius] # [K]
        self.NernstConstant = [x*(2.303*cv.R/(2*cv.F)) for x in self.Kelvin] #2.303RT/nF
        
        self.Density = [float(Reader1[i][5]) for i in range(self.RowStart,self.RowEnd)] #[g/cm^3] = [kg/L]
        self.Viscosity = [float(Reader1[i][6]) for i in range(self.RowStart,self.RowEnd)] #[g/cm*s] 
        self.SolubilityMagnetite_Fe = [float(Reader1[i][7]) for i in range(self.RowStart,self.RowEnd)]#[mol/kg] (Tremaine & LeBlanc)
        self.SolubilityChromite_Fe = [float(Reader1[i][8]) for i in range(self.RowStart,self.RowEnd)] #[mol/kg] chromite = FeCr2O4
        self.SolubilityFerrite_Fe = [float(Reader1[i][9]) for i in range(self.RowStart,self.RowEnd)] #[mol/kg] (Sandler & Kunig)
        self.SolubilityFerrite_Ni = [float(Reader1[i][10]) for i in range(self.RowStart,self.RowEnd)] #[mol/kg] (Sandler & Kunig)
        self.InnerOxideThickness = [float(Reader1[i][11]) for i in range(self.RowStart,self.RowEnd)] #[g/cm^2] 
        self.OuterMagnetiteThickness = [float(Reader1[i][12]) for i in range(self.RowStart,self.RowEnd)] #[g/cm^2]
        self.OuterFerriteThickness = [float(Reader1[i][13]) for i in range(self.RowStart,self.RowEnd)] #[g/cm^2] 
        self.OuterMetallicNickelThickness = [float(Reader1[i][14]) for i in range(self.RowStart,self.RowEnd)] #[g/cm^2] 
        self.OuterOxideThickness = [x+y+z for x,y,z in zip(self.OuterMagnetiteThickness,self.OuterFerriteThickness,\
                                    self.OuterMetallicNickelThickness)] #adds the lists together element by element
        
        #Electrochemical Parameters (Standard Potentials, interface specification not necessary)
        #Iron corrosion: Fe(s) <-> Fe2+(aq) + 2e-
        self.StandardEquilibriumPotential_Fe = [float(Reader2[i+8][1]) for i in range(self.RowStart,self.RowEnd)] #[V] 
        #Magnetite Precipitationitation: 3Fe(OH)2(s)<-> Fe3O4(s) + 2H2O(l) + 2H+(aq) + 2e-
        self.StandardEquilibriumPotential_Fe3O4Precipitation =[float(Reader2[i+8][2]) for i in range(self.RowStart,self.RowEnd)] #[V] 
        #Magnetite Dissolutionution: Fe3O4(s) + 2H2O(l) + 2H+(aq) + 2e- <-> 3Fe(OH)2(aq)
        self.StandardEquilibriumPotential_Fe3O4Dissolution = [float(Reader2[i+8][3]) for i in range(self.RowStart,self.RowEnd)] #[V] 
        #Nickel ferrite Precipitationitation: 2.4Fe(OH)2(s) + 0.6Ni(OH)2(s) <-> Ni0.6Fe2.4O4(s) + 2H2O(l) + 2H+(aq) + 2e-
        self.StandardEquilibriumPotential_FerritePrecipitation = [float(Reader2[i+8][4]) for i in range(self.RowStart,self.RowEnd)]#[V] 
        #Nickel ferrite Dissolutionution: Ni0.6Fe2.4O4(s) + 2H2O(l) + 2H+(aq) + 2e- <-> 2.4Fe(OH)2(aq) + 0.6Ni(OH)2(aq)
        self.StandardEquilibriumPotential_FerriteDissolution = [float(Reader2[i+8][5]) for i in range(self.RowStart,self.RowEnd)] #[V] 
        #Hydrogen reaction: 2H+(aq) + 2e- <-> H2(g)
        self.StandardEquilibriumPotential_H2 = [float(Reader2[i+8][6]) for i in range(self.RowStart,self.RowEnd)] #[V] 
        #Nickel corrosion: Ni(s) <-> Ni2+(aq) + 2e-
        self.StandardEquilibriumPotential_Ni = [float(Reader2[i+8][7]) for i in range(self.RowStart,self.RowEnd)]#[V] 
        #Nickel deposition: Ni(OH)2(aq)+ 2H+(aq) + 2e- <-> Ni(s) + 2H2O(l)
        self.StandardEquilibriumPotential_NiDeposition = [float(Reader2[i+8][8]) for i in range(self.RowStart,self.RowEnd)]#[V] 
        
        
DebyeHuckel=[float(Reader1[i][0]) for i in range(64,69)]
k_W=[float(Reader1[i][1]) for i in range(64,69)] 
k_Li=[float(Reader1[i][2]) for i in range(64,67)]
k_FeOH=[float(Reader1[i][3]) for i in range(64,69)] 
k_FeOH2= [float(Reader1[i][4]) for i in range(64,69)]
k_FeOH3=[float(Reader1[i][5]) for i in range(64,69)] 
k_NiOH=[float(Reader1[i][6]) for i in range(64,69)]
k_NiOH2=[float(Reader1[i][7]) for i in range(64,69)]            
k_NiOH3=[float(Reader1[i][8]) for i in range(64,69)]

#creates the 4 PHTS sections and their methods based on the Sections class template/blueprint
InletParameters=SectionParameters(4,11)
CoreParameters=SectionParameters(13,25)
OutletParameters=SectionParameters(27,36) 
SteamGeneratorParameters=SectionParameters(38,60)#ranges of rows for each sec. from txt file

