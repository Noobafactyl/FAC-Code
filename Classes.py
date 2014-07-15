



class simpleOperations(object):
    def __init__(self, x,y):
        self.x = x
        self.y = y 
    
    def productXY(self):
        prod = [self.x[i]*self.y[i] for i in range(len(self.x))]
        return prod
    
    def addOneX(self):
        print self.x
        self.x = [self.x[i]+1 for i in range(len(self.x))]
        return None
   
    def difference(self, x= None, y = None): # = means optional argument
        diff = [x[i]-y[i] for i in range(len(x))]
        return diff #always need to return something when using function
#         if (x==None and y!=None) or (x!=None or y==None):
#             raise NameError('ErrorNeeded')
#         if x!= None and y!=None:
#             if len(x)!=len(y):
#                 raise NameError('Length')
#         
#         if x== None:
#             x = self.x
#         if y==None:
#             y = self.y
        
    def difference(self):
        diff - [x[i]-y[i] for i in range(len(x))]
        return diff
    
    
    def printDefaults(self, x=1, y=1, z=1):
        print 'x=', x,'y=', y, 'z=', z
    
    def newFun(self):
        #increase values by one
        #calculate product 
        #addOneProd = [(x[i]+1)*y[i] for i in range(len(x))]
#       
        self.addOneX()
        prod = self.productXY()
        print 'x=', self.x, 'prod =', prod
        #product = [x[i]*y[i] for i in range(len(x))]
        return prod
#         return product
    
# x = [2,5,3]
# y = [4,5,7]

# ranObj = simpleOperations(x,y)
# ranObj.addOne()
# ranObj.product()

