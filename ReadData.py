import os
import csv
import numpy as np
import matplotlib.pylab as pyl

# Helper functions
#can't put optional argument before required one 
def ourFit(x, y, deg, xi = None):
    if xi == None:
        xi = x
    FitCoeff = np.polyfit(x, y, deg)
    fitPoly = np.poly1d(FitCoeff)
    return fitPoly(xi)

dirName = 'C:\Users\Owner\Documents\Data'
filename = 'sampleData.dat'
#reading data from the input file 
with open(os.path.join(dirName,filename), 'r') as csvfile:
    dataReader = csv.reader(csvfile,delimiter = ',')
    data = list()
    for row in dataReader:
        data.append(row)
# data = data [1:] #entries from the first to the end
# print data
data.pop(0) #alternate way to achieve above

x = [int(data[i][0]) for i in range(len(data))]
y = [float(data[i][1]) for i in range(len(data))]
z = [float(data[i][2]) for i in range(len(data))]

#Fitting data to a polynomial
deg = 2
yFitValues = ourFit(x,y,deg)
xRandom = [4.5, 6.5, 8.5]
yFitRandom = ourFit(x,y,deg, xRandom)
zFitValues = ourFit(x,z,deg)

# yFitCoeff= np.polyfit(x, y, deg)
# yFitPoly = np.poly1d(yFitCoeff)
# yFitValues = fitPoly(x)
 
# zFitCoeff= np.polyfit(x, z, deg)
# zFitPoly = np.poly1d(zFitCoeff)
# zFitValues = zFitPoly(x)

# print yFitCoeff
# print yFitPoly
# print yFitValues

pyl.figure()
pyl.plot(x,y, 'o')
pyl.plot(x,yFitValues)
pyl.plot(x,z, 'o')

pyl.plot(xRandom,yFitRandom, 'x')

pyl.plot(x,zFitValues)
pyl.xlim(0,10)
pyl.show()


# with open(os.path.join(dirName,filename), 'w') as csvfile:
#     writer = csv.writer(csvfile, delimiter = ',')
#     writer.writerow(['x','y','z'])
#     writer.writerows(dataList)