#!/usr/bin/env python
# coding: utf-8

# # Notebook Version of Biological Age Code

# Edits 03/29/2020: V2.0
# * I made automatic groups detection (detects 1 or 2 groups)
# * Removed ggplot and changed Visualize.plots() to plot with matplotlib
# * The code is now compatible with pandas DataFrames. Just load in the DataFrame instead of converting to list.
# * Removed graph's ability to show SLR results of age, primarykey, samp_wt, sex against age (skip loop if idx belongs to important indices)
# 
# Edits 7/15/2019:  
# * I changed one line so that it's no longer reliable on column names being specifically 'age', 'seqn', 'group', etc.
# * Edited so it works as an importable notebook  
#     **Note: if you want to change the parameters GRAPHON or AGEON, you have to do it through the KDM method in the Methods Class in self.__init__()

# In[13]:


from numpy import mean, corrcoef, array
import math
from scipy import stats
import scipy as sp
from sklearn import datasets, linear_model
#import statsmodels.api as sm
import pandas as pd
import json
import matplotlib.pyplot as plt
import numpy as np
import datetime


# In[4]:


"""BIOLOGICAL AGE ALGORITHM
NOTE: RUN IN ANACONDA FOR PYTHON 2.7 ONLY

How to use commands:
- Create an instance of the data, then apply methods.
	- e.g. To find means, use Summary(dataset).mean()

Classes:
Summary:
	Methods: 
	.mean() [calculates mean]
	.corr(col1, col2) [calculates Pearson's Corr Coefficient & p-value] *col2 can be a list or one column. Both parameters are OPTIONAL
	.view [shows data as a dataframe]
Methods:
	Methods:
	.KDM() [gives Klemera & Doubal equation for biological age]
	.cleandata(col1,col2) [gets rid of all rows with missing values]
	.calcLR() [calculates linear regressions of all vars with age as predictor]
Visualize:
	.plot(X,Y,[stratifying_variable]) [view plot of biological age vs chronological age]
	.view [analyze data frame output]

"""


# ### Open Data

# In[15]:


################## CLASS SUMMARY #######################

#Double check means
#Make a dictionary of the means of each thing, turn into a function
class Summary(object):
    def __init__(self, dataframe, age, samp_wtIndex, primarykey):
        header_names = list(dataframe)
        temp = dataframe.values.tolist()
        temp.insert(0, header_names)
        
        self.primarykey = primarykey
        self.data = temp
        self.view = dataframe
        self.age = age
        self.samp_wt = samp_wtIndex

    def mean(self, lowerboundage=0, upperboundage=999):
        avgdict = {}
        floatdum = 3.0
        # cols 0 - 16
        for col in range(len(self.data[0])):
            avglist = []
            sampwtlist = []
            
            if col == self.primarykey:
                continue
            
            for line in self.data[1:]:

                if line[self.age] >= lowerboundage and line[self.age] <= upperboundage and type(line[col]) == type(floatdum):
                    avglist.append(line[col])
                    sampwtlist.append(line[self.samp_wt])
                
            avgdict[self.data[0][col]] = np.mean(avglist) #sum([i * j for i,j in zip(avglist, sampwtlist)])/len(avglist)

        self.make_pretty(avgdict,"MEANS:")
        return avgdict

    #returns a pretty row of results
    def make_pretty(self,dictionary,titlestring):
        print("\n" + titlestring)
        keylist = sorted(dictionary.keys())
        for key in keylist:
            if len(key) < 8 :
                print(key + '\t\t' + str(dictionary[key]))
            else:
                print(key + '\t' + str(dictionary[key]))

    def __str__(self):
        return "There are {} observations and {} variables in this dataset.".format(len(self.data),len(self.data[0]))

    #returns the Pearson Correlation coefficient of each row in the dataset to a variable
    #input parameters: column number of you want to compare to. Second one to compare to
    #first value can only be one value. If you put a #<0 for y, will be a list of all vals. default is all vals
    def corr(self, x=-2, y=-1):
        corrdict = {}
        repeat = False
        #establishes what we'll be correlating to x
        if y >= 0:
            collist = [y]
        elif y < 0:
            collist = range(len(self.data[0]))
        else:
            print("Please enter a valid input")
            quit()

        for col in collist:
            corrlistx = []
            corrlisty= []
            for line in self.data:
                if type(line[x]) == type(1.2):
                    corrlistx.append(line[x])
                elif type(line[x]) != type(1.2):
                    corrlistx.append('.')

                if type(line[col]) == type(1.2):
                    corrlisty.append(line[col])
                else:
                    corrlisty.append('.')

            try:
                corrdict[self.data[0][col]] = stats.pearsonr(*self.findmissing(corrlistx,corrlisty))
            except:
                print("Cannot perfom correlation on", self.data[0][col])

        self.make_pretty(corrdict,"CORRELATIONS: (r, p-value)")
        return corrdict

    #throws away values where at least one value is missing
    def findmissing(self,corrlist1,corrlist2):
        for i in range(len(corrlist1)):
            try:
                if corrlist1[i] != type(1.2):
                    del corrlist1[i]
                    del corrlist2[i]
            except: pass

        for i in range(len(corrlist2)):
            try:
                if corrlist2[i] != type(1.2):
                    del corrlist1[i]
                    del corrlist2[i]
            except: pass
        return (corrlist1, corrlist2)


# In[16]:



################## CLASS METHODS #######################
class Methods(object):
    def __init__(self, dataframe, cache_fname, age, genderindex, primarykey, samp_wt, output):
        # convert DataFrame to List Format
        header_names = list(dataframe)
        temp = dataframe.values.tolist()
        temp.insert(0, header_names)
        
        self.data = temp
        self.age = age
        self.genderindex = genderindex
        self.primarykey = primarykey
        self.samp_wt= samp_wt
        self.cache_fname = cache_fname
        #print("Looking for data in cache...")
        try:
            cache_fhnd = open(self.cache_fname,'r')
            self.CACHE_DICT = json.loads(cache_fhnd.read())
            cache_fhnd.close()
        except:
            self.CACHE_DICT = {}
        self.AGEON = True # age corrector toggle
        self.GRAPHON = True # do you want to see linear regression graphs for each variable on age?
        self.savepath = output
        
        # determine how groups there are
        self.GroupSet = pd.unique(dataframe[header_names[self.genderindex]].values).astype('int')

        if len(self.GroupSet) > 2:
            raise ValueError('There are more than 2 groups in this dataset')

###############################################################################################
########################## Klemera and Doubal Method ##########################################
###############################################################################################

#Note on CACHE_DICT vs cachedict: CACHE_DICT is the "global local" variable for the instance. cachedict is input as a parameter if the value is in CACHE_DICT
    def KDM(self):
        #have this one combine the two datasets, for men and for women
        #and append to a list
        finaldict = {}

        for sex in self.GroupSet:
            if sex == 0:
                sexname = "Males"
            else: sexname = "Females"

            if str(sex) in self.CACHE_DICT:
                #print("Collecting parameters from cache for %s..." % (sexname))
                self.calcLR(sex=sex,cachedict = self.CACHE_DICT[str(sex)]) #finds self.regressiondict and makes self.newdata
                finaldict[sex] = self.calcBA(correctiondict = self.correctionterm(self.calcBA()))
            else:
                #print("Could not find parameters in the cache for %s..." % (sexname))
                self.calcLR(sex=sex) #makes self.regressiondict
                finaldict[sex] = self.calcBA(correctiondict = self.correctionterm(self.calcBA()))

        #print("Merging data one last time...")
        #Double check this...
        if 'BA' not in self.data[0]:
            for row in self.data:
                for sexkey in finaldict: #going through each SEQN
                    for key in finaldict[sexkey]:
                        if key == row[self.primarykey]:          #if you found the right SEQN for that row, append the Corrected Biological Age to the end
                            row.append(finaldict[sexkey][key])
        if 'BA' not in self.data[0]:
            self.data[0].append('BA')
        if 'BAC' not in self.data[0]:
            self.data[0].append('BAC')

        #print("DONE!")

        path = self.savepath
        fopen = open(path,'w+')
        
        #print("Saving Data...")
        for row in self.data:
            for i in range(len(row)):
                row[i] = str(row[i])
            fopen.write(",".join(row) + "\n")

    #   fopen.write(json.dumps(finaldict))
        fopen.close()
        #print("DONE!")
        return self.data

    #REMEMBER TO CONSIDER SEX IN THIS...CALCLR CALCULATES FOR ONLY ONE SEX AT A TIME BUT YOU NEED A DICTIONARY THAT ACCOUNTS FOR BOTH
    #this uses the regression results to create all the baseline variables for KDM calculations
    def correctionterm(self, datatuple, agemax = 24, agemin = 36): # changed months to 2-3 years
        m = datatuple[1] #this tells you how many covariates we have
        BAdict = datatuple[0] #this was BasicBAdict from the calcBA function
        n = len(self.newdata)-1
        delcounter = 0

        CorrectedBAdict = {}
        
        #Merge BAdict with self.newdata (remember, this is the dataset with only one sex)
        #print("Merging data...")
        for row in self.newdata:
            try:
                row.append(BAdict[row[self.primarykey]][0]) #will append to the end of the row the BA value from the dictionary, using the primary key as the identifier
            except:
                del(row)
                delcounter += 1

        if delcounter > 0: 
            print("Warning: %d rows deleted. Check the integrity of the data or check the code." % (delcounter))

        #STANDARD DEVIATION CALCULATIONS
        #Calculating first term
        errcalc  = []
        #print("Calculating Correction Term...")
        for row in self.newdata:
            BA = row[-1] #the last value appended to each row
            errcalc.append(BA - row[self.age])
        errcalc = np.array(errcalc)
        #Calculating standard deviation

        for row in self.newdata:
            rchar = BAdict[row[self.primarykey]][1]
            stderr = (np.var(errcalc) - (((1-(rchar**2))/(rchar**2)) * (((agemax - agemin)**2)/(12*m))))
            
            # Extra step: Linearly transform so that SBA maintains same mean but now linearly increases with age, so difference is 5 between CAmax and CAmin
            if self.AGEON == True:
                agecorrector = 5/(agemax-agemin)
                std = np.sqrt(stderr) - 2.5 + (agecorrector * row[self.age]) # adjust based on age
                stderr = std ** 2
                #print(stderr)

            CorrectedBAdict[row[self.primarykey]] = stderr

        return CorrectedBAdict
    #############
    def cov(self, x, y, w):
        """Weighted Covariance"""
        return np.sum(w * (x - np.average(x, weights = w)) * (y - np.average(y, weights = w))) / np.sum(w)

    def corr(self, x, y, w):
        """Weighted Correlation"""
        return self.cov(x, y, w) / np.sqrt(self.cov(x, x, w) * self.cov(y, y, w))
    #############

    #This takes all the data from the Linear Regression function and calcuates the baseline predicted age and rchar.
    #It then passes all the variables to correctionterm method (if you use the .KDM() method), which will aggregate these variables into calculating the corrected BA.
    def calcBA(self, correctiondict={}):
        #append another value to the end of each row?
        BasicBAdict = {}
        CorrectedBAdict = {}

        if correctiondict == {}:
            #print("Calculating initial BA without Correction...")
            for row in self.newdata:
                numeratorlist = []
                denominatorlist = []
                rcharlistnumerator = []
                rcharlistdenominator = []
            #append BA data to each row
                for key in self.regressiondict:
                #append numerators to a list and denominators to a list, then sum and divide
                    covar = self.regressiondict[key] #references one column or variable in study
                    k = covar[0] #slope
                    s = covar[3] #MSE
                    q = covar[1] #intercept
                    r = covar[2] #r-value

                    colindex = self.heading.index(key) #find the column of the data

                    try:
                        numeratorlist.append((row[colindex]-q)*(k/(s**2)))
                        denominatorlist.append((k/s)**2)
                        rcharlistnumerator.append(((r**2)/(math.sqrt(1-(r**2)))))
                        rcharlistdenominator.append((r/(math.sqrt(1-(r**2)))))
                    except:
                        print("r:",r, "var:", key)
                        #print(self.regressiondict)
                        quit()
                #rcalculations
                rchar = sum(rcharlistnumerator)/sum(rcharlistdenominator)
                #create a dictionary with SEQN as key

                ####################print((sum(numeratorlist)/sum(denominatorlist), rchar))
                BasicBAdict[row[self.primarykey]] = (sum(numeratorlist)/sum(denominatorlist), rchar)
            #returns a tuple with dictionary first, m second
            return (BasicBAdict,len(self.regressiondict))
        else:
            #print("Calculating final Age with Corrected Values...")
            for row in self.newdata:
                numeratorlist = []
                denominatorlist = []

                #append BA data to each row
                for key in self.regressiondict:
                #append numerators to a list and denominators to a list, then sum and divide
                    covar = self.regressiondict[key] #references one column or variable in study
                    k = covar[0] #slope
                    s = covar[3] #MSE
                    q = covar[1] #intercept
                    r = covar[2] #r-value
                    colindex = self.heading.index(key) #find the column of the data

                    numeratorlist.append((row[colindex]-q)*(k/(s**2))+(row[self.age]/correctiondict[row[self.primarykey]]))
                    denominatorlist.append(((k/s)**2)+(1/correctiondict[row[self.primarykey]]))

                #create a dictionary with SEQN as key

                CorrectedBAdict[row[self.primarykey]] = sum(numeratorlist)/sum(denominatorlist)
            #returns a final dictionary of Corrected BA
            return CorrectedBAdict



    #remember female = 1, means female
    #this does the inital linear regressions: MODEL xx = age
    def calcLR(self, sex=0, cachedict = {}):
        #creates a new dataset for ONLY females or only males
        self.newdata = [line for line in self.data if line[self.genderindex] == "" or line[self.genderindex] == sex]
        self.df = pd.DataFrame(self.newdata[1:],columns=self.data[0])
        self.heading = list(self.df)
        self.regressiondict = {} #saves tuple: (slope,intercept,r_value,std_err)
        if sex == 0:
            sexname = "Males"
        else:
            sexname = "Females"
        #check if there is already a value in cachedict:
        if cachedict == {}:
            #print("Calculating regressions for %s..." % (sexname))
        

            # REGRESSION PART #
            for idx,col in enumerate(self.heading):
                #print(col)
                if idx in [self.age, self.primarykey, self.genderindex, self.samp_wt]:
                    continue
                
                regX = np.array(self.df[list(self.df)[self.age]]).reshape(-1, 1) #age
                regY = np.array(self.df[[col]]) #other var

                regr = linear_model.LinearRegression()

                try:
                # Train the model using the training sets
                    weight = np.array(self.get_column(self.samp_wt)).reshape(len(self.get_column(self.samp_wt)),1)

                    regr.fit(regX,regY,sample_weight=self.get_column(self.samp_wt))
                    #print regr.get_params(deep=False)

                    #The coefficients
                    #print 'Coefficients: \n', regr.coef_
                    #print 'Intercept: \n', regr.intercept_
        
                    # The mean squared error
                    #print "Mean squared error: %.2f" % mean((regr.predict(regX) - regY) ** 2) #MSE
                    # Explained variance score: 1 is perfect prediction
                    #print 'Variance score: %.2f' % (regr.score(regX, regY)) #variance
                    #need slope, intercept, r-value

                    # Plot outputs
                    if self.GRAPHON == True:
                        plt.scatter(regX, regY,  color='black')
                        plt.plot(regX, regr.predict(regX), color='blue', linewidth=3)
                        plt.xticks(())
                        plt.yticks(())
                        plt.title('%s vs Age' % (col))
                        plt.xlabel('age')
                        plt.ylabel(col)
                        plt.show()

                    #do linear regressions analysis
                    slope = regr.coef_[0][0]
                    intercept = regr.intercept_[0]
                    # weighted vs unweighted corr seems to make no difference
                    r_value = self.corr(regX, regY, weight) #stats.pearsonr(regX,regY)[0][0]
                    std_err = math.sqrt(mean((regr.predict(regX) - regY) ** 2))
                    #print(slope, intercept, r_value, std_err)

                except:
                    print("Cannot do LINEAR REGRESSION analysis for %s" % (col))
                    slope = 0
                    intercept = 0
                    r_value = 0
                    std_err = 0

                #CHECK THAT ONLY PROPER VALUES ARE ADDED
                if slope and intercept and r_value and std_err != 0:
                    #add to dictionaries
                    self.regressiondict[col] = (slope, intercept, r_value, std_err)



            #add to cache
            cache_fhnd = open(self.cache_fname,'w')
            self.CACHE_DICT[str(sex)] = self.regressiondict
            cache_fhnd.write(json.dumps(self.CACHE_DICT))
            cache_fhnd.close()
        #if there is, make it into self.regressiondict
        else:
            self.regressiondict = cachedict
        #self.make_pretty(self.regressiondict,"REGRESSION VALUES (slope, intercept, r_value, std_err)")
        return self.regressiondict


### Utility Functions ###

    #returns the mean of each row in the dataset
    def make_pretty(self,dictionary,titlestring):
        print("\n" + titlestring)
        keylist = sorted(dictionary.keys())
        for key in keylist:
            if len(key) < 8 :
                print(key + '\t\t' + str(dictionary[key]))
            else:
                print(key + '\t' + str(dictionary[key]))

    #extracts one column from the dataset
    def get_column(self,col):
        column = [line[col] for line in self.newdata]
        return column[1:]

    #throws away values where at least one value is missing
    def cleandata(self,col1,col2):
        for i in range(len(col1)):
            try:
                if col1[i] != type(1.2):
                    del col1[i]
                    del col2[i]
            except: pass

        for i in range(len(col2)):
            try:
                if col2[i] != type(1.2):
                    del col1[i]
                    del col2[i]
            except: pass
        if len(col1) == len(col2):
        #   print "Dimensions match for both columns"
            return (col1, col2)
        else:
            print("Cannot match data by length, check data dimensions.")
            return (len(col1), len(col2))


# In[17]:


class Visualize(object):
    def __init__(self,data):
        self.data = data

    def view(self):
        return self.data

    def plot(self, inp1, inp2, inp3): # (y, x, color by)
        new_df = pd.DataFrame()

        new_df[inp1] = pd.to_numeric(self.data[inp1])
        new_df[inp2] = pd.to_numeric(self.data[inp2])
        new_df[inp3] = self.data[inp3]

        # plot
        fig, ax = plt.subplots(figsize=(10,10))
        
        for group in pd.unique(new_df[inp3]):
            plt.scatter(new_df[new_df[inp3] == group][inp1], new_df[new_df[inp3] == group][inp2])
        
        ax.legend(pd.unique(new_df[inp3]))
        plt.ylabel(inp1)
        plt.xlabel(inp2)
        plt.show()


# ### Run Function

# In[18]:


def KDM_model(trainset, testset, cachename, output_filename, age_index, genderindex, primaryindex, samp_wt_index):

    # Build model
    output_filename_train = output_filename + '_train.csv'
    model = Methods(trainset, cachename, age = age_index, genderindex = genderindex, primarykey = primaryindex, samp_wt = samp_wt_index, output = output_filename_train)
    model.KDM() # train model

    # test using trained model
    output_filename_test = output_filename + '_test.csv'
    test_model = Methods(testset, cachename, age = age_index, genderindex = genderindex, primarykey = primaryindex, samp_wt = samp_wt_index, output = output_filename_test)
    results = test_model.KDM()
    results = pd.DataFrame(results[1:], columns = results[0])

    # calculate stats
    stats = return_stats(results['BAC'].astype('float64'), results[list(results)[age_index]].astype('float64'))
    
    return (stats)


# ### Return Stats

# In[19]:


#Return statistics for a correlation
def return_stats(x,y, dec=3):
    slope, intercept, r_value, p_value, std_err = sp.stats.linregress(x,y)
    
    Mr = abs(x-y).mean()
    Medr = np.median(abs(x-y))
    Mpr = np.median(100*abs((x-y)/x))
    rs = r_value**2
    print(r_value)
    
    return(Mr, Medr, Mpr, p_value, rs)


# In[20]:


print('You are running BA_NB_Final at', datetime.datetime.now())


# In[ ]:




