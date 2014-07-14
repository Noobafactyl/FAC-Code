import os
import csv
import random as rand
import matplotlib.pylab as pyl
x = range(1,10)
print x 
y = [(i+rand.random())**2 for i in x]
z = [(i+rand.random()*2)**2 for i in x]
print y 
pyl.plot(x,y,'o') 
pyl.plot(x,z,'o') 
pyl.show()

dirName = 'C:\Users\Owner\Documents\Data'
filename = 'sampleData.dat'

if not os.path.exists(dirName):
    os.path.dirName

dataList = list()
[dataList.append([x[i],y[i],z[i]]) for i in range(len(x))]
print dataList 
with open(os.path.join(dirName,filename), 'w') as csvfile:
    writer = csv.writer(csvfile, delimiter = ',')
    writer.writerow(['x','y','z'])
    writer.writerows(dataList)