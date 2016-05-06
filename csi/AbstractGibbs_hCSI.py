from numpy.random import choice as mysample;
import random;


import matplotlib.pyplot as plt;
class AbstractGibbs(object):
    #rvObjList is a list of random variables object, each object has
    #a function to compute its conditional probability given values for other
    #variable in the list
    #sampledValues is a list of lists where each 'inner' list is the value of a varialbe sampled by the sampler at each update step
    def __init__(self, rvObjList):
        self.rvl = rvObjList;
        self.indexList = range(len(self.rvl));
        self.sampledValues = [];
        for i in self.indexList:
            self.sampledValues.append([]);
    
    #update step
    def gibbsUpdate(self):
        #run through all variables
        random.seed(1);
        for i in self.indexList:
            #for each variable, build a distribution accross all possible values for that variable
            #based on current values of other variable 
            distribution = self.rvl[i].getConditionalDistribution(self.rvl[:i]+self.rvl[i+1:]);
            #sample a random value based on the distribution
            val = mysample(self.rvl[i].getValRange(),size=1,p=distribution);
            self.rvl[i].setCurrentValue(val[0]);

    #run Gibbs udpate with with "repeats" number of cycles        
    def sample(self, repeats=1):
        for i in range(repeats):
                if i % 100 == 0:
                        print("Done "+str(i));
                self.gibbsUpdate();
                for j in self.indexList:
                        self.sampledValues[j].append(self.rvl[j].getCurrentValue());
                


#class RandomVariable has 4 important methods:
#1) conditionalProbability
#2/3) get/setCurrentValue
#4) getConditionalDistribution
class RandomVariable(object):
        #currentValue is the value that the variable is taking on
        #valRange is a range of all possible values that this variable can take
        #normalize should probably be set to true, to make sure that distribution is always between [0,1]
        #distribution is a list of weight corresponding to each value in valRange
    def __init__(self,val = None, valRange = None, normalize = False):
        self.currentValue = val;
        self.valRange = valRange;
        self.normalize = normalize;
        self.distribution = list(range(len(self.valRange)));
    #getter for valRange
    def getValRange(self):
        return self.valRange;
    #getter for currentValue
    def getCurrentValue(self):
        return self.currentValue;
        
    #set currentValue and do any updates that might come afterwards
    def setCurrentValue(self, val):
        self.currentValue = val;
        self.updateAfterSettingValue(val);
        
    #to be implemented by anything extending these clases
    def updateAfterSettingValue(self, val):
        pass
    #compute the conditional probability of current RV given other RVs (and their values)        
    def conditionalProbability(self,value,rvList):
        pass
    
    #compute distribution
    def getConditionalDistribution(self,rvList):
        i = 0;
        #for each possible value of the RV
        for v in self.getValRange():
                #compute the conditional probability that this RV takes this value given other RVs (and their values)
            p =  self.conditionalProbability(v,rvList) ;
            #set the distribution
            self.distribution[i] = p;
            i = i + 1;
        #########IF NEED TO NORMALIZE DISTRIBUTION MATRIX (NOT ALWAYS NEEDED)########
        self.normalizeDistribution();
                
        return self.distribution;
    
    def normalizeDistribution(self):
        pass;


'''
#### CLASSes TO TEST

from matplotlib import mlab as mlab;
import math;
import random;
class IntegerRandomVariable(RandomVariable):
    def conditionalProbability(self,value,rvList):
        #return math.factorial(200)/math.factorial(value)/math.factorial(200-value)*math.pow(0.8,value)*math.pow(0.2,200-value)
        return mlab.normpdf(value,100 + 0.8*(rvList[0].getCurrentValue() - 100),math.sqrt(1-0.64));

class DiscreteRV1(RandomVariable):
    def conditionalProbability(self,value,rvList):
        return binomial(rvList[0].getCurrentValue(),16,value);
class DiscreteRV2(RandomVariable):
    def conditionalProbability(self,value,rvList):
        a = 2;
        b = 4;
        return beta(value,rvList[0].getCurrentValue() + a, 16 - rvList[0].getCurrentValue() + b);
    


from math import gamma as gamma;
def binomial(p,n,k):
    return math.factorial(n)/math.factorial(k)/math.factorial(n-k)*math.pow(p,k)*math.pow(1-p,n-k);
def beta(v, x,y):
    return gamma(x+y)/gamma(x)/gamma(y)*math.pow(v,x-1)*math.pow(1-v,y-1);
'''   
            
        
