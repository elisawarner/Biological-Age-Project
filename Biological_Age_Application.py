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
from numpy import mean, corrcoef, array
import math
from scipy import stats
from sklearn import datasets, linear_model
#import statsmodels.api as sm
import pandas as ps
import ggplot as gg
import json
import matplotlib.pyplot as plt

#GLOBAL VARS
age = 8 #define the column 'age' for clarity
primarykey = 0  #SEQN or other primary identifier
genderindex = 12 #define the column for gender (should be binary)

#OPEN DATA
path = 'BioAgeSample_NHANESIII_without_fvc.txt'
#path = 'test2.csv'
fhnd = open(path,'r')
data = []
CACHE_FNAME = 'projectcache.txt' #create a cache for parameters dictionary to save time

#FORMAT DATA
for line in fhnd:
	linelist = line.strip().split('\t')
	data.append(linelist)

for line in data[1:]:
	for i in range(len(line)):
		try:
			line[i] = float(line[i])
		except:
			pass


#CHECK DATA
#for line in data[0:3]:
#	print line

fhnd.close()

################## CLASS SUMMARY #######################

#Double check means
#Make a dictionary of the means of each thing, turn into a function
class Summary(object):
	def __init__(self,data):
		self.data = data
		self.view = ps.DataFrame(self.data[1:],columns=self.data[0])

	def mean(self):
		avgdict = {}
		floatdum = 3.0
		# cols 0 - 16
		for col in range(len(self.data[0])):
			avglist = []
#			keyname = 'col_' + str(col)
			for line in self.data:
				if line[age] >= 30 and line[age] <= 75 and type(line[col]) == type(floatdum):
					avglist.append(line[col])
#			print "Mean age:", mean(avglist)
			avgdict[self.data[0][col]] = mean(avglist)

		self.make_pretty(avgdict,"MEANS:")
		return avgdict

	#returns a pretty row of results
	def make_pretty(self,dictionary,titlestring):
		print "\n" + titlestring
		keylist = sorted(dictionary.keys())
		for key in keylist:
			if len(key) < 8 :
				print key + '\t\t' + str(dictionary[key])
			else:
				print key + '\t' + str(dictionary[key])

	def __str__(self):
		return "There are {} observations and {} variables in this dataset.".format(len(self.data),len(self.data[0]))

	#returns the Pearson Correlation coefficient of each row in the dataset to a variable
	#input parameters: column number of you want to compare to. Second one to compare to
	#first value can only be one value. If you put a #<0 for y, will be a list of all vals. default is all vals
	def corr(self, x=age, y=-1):
		corrdict = {}
		repeat = False
		#establishes what we'll be correlating to x
		if y >= 0:
			collist = [y]
		elif y < 0:
			collist = range(len(self.data[0]))
		else:
			print "Please enter a valid input"
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
				print "Cannot perfom correlation on", self.data[0][col]

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


################## CLASS METHODS #######################
class Methods(object):
	def __init__(self,data):
		self.data = data
		print "Looking for data in cache..."
		try:
			cache_fhnd = open(CACHE_FNAME,'r')
			self.CACHE_DICT = json.loads(cache_fhnd.read())
			cache_fhnd.close()
		except:
			self.CACHE_DICT = {}


### Klemera and Doubal Method ###

#Note on CACHE_DICT vs cachedict: CACHE_DICT is the "global local" variable for the instance. cachedict is input as a parameter if the value is in CACHE_DICT
	def KDM(self):
		#have this one combine the two datasets, for men and for women
		#and append to a list
		finaldict = {}

		for sex in [1,0]:
			if sex == 0:
				sexname = "Males"
			else: sexname = "Females"

			if str(sex) in self.CACHE_DICT:
				print "Collecting parameters from cache for %s..." % (sexname)
				self.calcLR(sex=sex,cachedict = self.CACHE_DICT[str(sex)]) #finds self.regressiondict and makes self.newdata
				finaldict[sex] = self.calcBA(correctiondict = self.correctionterm(self.calcBA()))
			else:
				print "Could not find parameters in the cache for %s..." % (sexname)
				self.calcLR(sex=sex) #makes self.regressiondict
				finaldict[sex] = self.calcBA(correctiondict = self.correctionterm(self.calcBA()))

		print "Merging data one last time..."
		#Double check this...
		for row in self.data:
			for sexkey in finaldict: #going through each SEQN
				for key in finaldict[sexkey]:
					if key == row[primarykey]:          #if you found the right SEQN for that row, append the Corrected Biological Age to the end
						row.append(finaldict[sexkey][key])
		self.data[0].append('BA')
		self.data[0].append('BAC')

		print "DONE!"

		path = 'test2.csv'
		fopen = open(path,'w+')
		
		print "Saving Data..."
		for row in self.data:
			for i in range(len(row)):
				row[i] = str(row[i])
			fopen.write(",".join(row) + "\n")

	#	fopen.write(json.dumps(finaldict))
		fopen.close()
		print "DONE!"
		return self.data

	#REMEMBER TO CONSIDER SEX IN THIS...CALCLR CALCULATES FOR ONLY ONE SEX AT A TIME BUT YOU NEED A DICTIONARY THAT ACCOUNTS FOR BOTH
	#this uses the regression results to create all the baseline variables for KDM calculations
	def correctionterm(self, datatuple, agemax = 75, agemin = 30):
		m = datatuple[1] #this tells you how many covariates we have
		BAdict = datatuple[0] #this was BasicBAdict from the calcBA function
		n = len(self.newdata)-1
		delcounter = 0

		CorrectedBAdict = {}
		
		#Merge BAdict with self.newdata (remember, this is the dataset with only one sex)
		print "Merging data..."
		for row in self.newdata[1:]:
			try:
				row.append(BAdict[row[primarykey]][0]) #will append to the end of the row the BA value from the dictionary, using the primary key as the identifier
			except:
				del(row)
				delcounter += 1

		if delcounter > 0: 
			print "Warning: %d rows deleted. Check the integrity of the data or check the code." % (delcounter)

		#STANDARD DEVIATION CALCULATIONS
		#Calculating first term
		errcalc  = []
		print "Calculating Correction Term..."
		for row in self.newdata[1:]:
			BA = row[-1] #the last value appended to each row
			errcalc.append(BA - row[age])

		#Calculating standard deviation
		for row in self.newdata[1:]:
			rchar = BAdict[row[primarykey]][1]
			stderr = (((sum(errcalc)-(sum(errcalc)/n))**2)/n) - (1-(rchar**2))/(rchar**2) * (((agemax - agemin)**2)/(12*m))
			CorrectedBAdict[row[primarykey]] = stderr
		return CorrectedBAdict

	#This takes all the data from the Linear Regression function and calcuates the baseline predicted age and rchar.
	#It then passes all the variables to correctionterm method (if you use the .KDM() method), which will aggregate these variables into calculating the corrected BA.
	def calcBA(self, correctiondict={}):
		#append another value to the end of each row?
		BasicBAdict = {}
		CorrectedBAdict = {}

		if correctiondict == {}:
			print "Calculating initial BA without Correction..."
			for row in self.newdata[1:]:
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
					colindex = self.newdata[0].index(key) #find the column of the data

					try:
						numeratorlist.append((row[colindex]-q)*(k/(s**2)))
						denominatorlist.append((k/s)**2)
						rcharlistnumerator.append(((r**2)/(math.sqrt(1-(r**2)))))
						rcharlistdenominator.append((r/(math.sqrt(1-(r**2)))))
					except:
						print "r:",r, "var:", key
						print self.regressiondict
						quit()
				#rcalculations
				rchar = sum(rcharlistnumerator)/sum(rcharlistdenominator)
				#create a dictionary with SEQN as key

				BasicBAdict[row[primarykey]] = (sum(numeratorlist)/sum(denominatorlist), rchar)
			#returns a tuple with dictionary first, m second
			return (BasicBAdict,len(self.regressiondict))
		else:
			print "Calculating final Age with Corrected Values..."
			for row in self.newdata[1:]:
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
					colindex = self.newdata[0].index(key) #find the column of the data

					numeratorlist.append((row[colindex]-q)*(k/(s**2))+(row[age]/correctiondict[row[primarykey]]))
					denominatorlist.append(((k/s)**2)+(1/correctiondict[row[primarykey]]))

				#create a dictionary with SEQN as key
				CorrectedBAdict[row[primarykey]] = sum(numeratorlist)/sum(denominatorlist)
			#returns a final dictionary of Corrected BA
			return CorrectedBAdict



	#remember female = 1, means female
	#this does the inital linear regressions: MODEL xx = age
	def calcLR(self, sex=0, cachedict = {}):
		#creates a new dataset for ONLY females or only males
		self.newdata = [line for line in self.data if line[genderindex] == "female" or line[genderindex] == sex]
		self.df = ps.DataFrame(self.newdata[1:],columns=data[0])
		self.heading = self.newdata[0]
		self.regressiondict = {} #saves tuple: (slope,intercept,r_value,std_err)
		if sex == 0:
			sexname = "Males"
		else:
			sexname = "Females"
		#check if there is already a value in cachedict:
		if cachedict == {}:
			print "Calculating regressions for %s..." % (sexname)
		

			# REGRESSION PART #
			for col in self.heading:
				print col
				
				regX = self.df[['age']] #age
				regY = self.df[[col]] #other var


				regr = linear_model.LinearRegression()
				try:
					# Train the model using the training sets
					regr.fit(regX,regY,sample_weight=array(self.get_column(9)))
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
					r_value = regr.score(regX, regY)
					std_err = math.sqrt(mean((regr.predict(regX) - regY) ** 2))
					print slope, intercept, r_value, std_err

				except:
					print "Cannot do LINEAR REGRESSION analysis for %s" % (col)
					slope = 0
					intercept = 0
					r_value = 0
					std_err = 0

				#CHECK THAT ONLY PROPER VALUES ARE ADDED
				if slope and intercept and r_value and std_err != 0 and col != 'age' and col != 'seqn' and col != 'female' and col != 'samp_wt':
					#add to dictionaries
					self.regressiondict[col] = (slope, intercept, r_value, std_err)



			#add to cache
			cache_fhnd = open(CACHE_FNAME,'w')
			self.CACHE_DICT[str(sex)] = self.regressiondict
			cache_fhnd.write(json.dumps(self.CACHE_DICT))
			cache_fhnd.close()
		#if there is, make it into self.regressiondict
		else:
			self.regressiondict = cachedict
		self.make_pretty(self.regressiondict,"REGRESSION VALUES (slope, intercept, r_value, std_err)")
		return self.regressiondict


### Utility Functions ###

	#returns the mean of each row in the dataset
	def make_pretty(self,dictionary,titlestring):
		print "\n" + titlestring
		keylist = sorted(dictionary.keys())
		for key in keylist:
			if len(key) < 8 :
				print key + '\t\t' + str(dictionary[key])
			else:
				print key + '\t' + str(dictionary[key])

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
		#	print "Dimensions match for both columns"
			return (col1, col2)
		else:
			print "Cannot match data by length, check data dimensions."
			return (len(col1), len(col2))

class Visualize(object):
	def __init__(self,data):
		self.data = ps.DataFrame(data[1:],columns=data[0])

	def view(self):
		return self.data

	def plot(self, inp1, inp2, inp3=None):
		p = gg.ggplot(gg.aes(x=inp1, y=inp2, color=inp3), data=self.data) + \
		gg.geom_point()

		print p

#Test cases
"""
t1path = 'test2.csv'
t_fhnd = open(t1path,'r')
tdata = []

#FORMAT DATA
for line in t_fhnd:
	linelist = line.strip().split(',')
	tdata.append(linelist)

for line in tdata[1:]:
	for i in range(len(line)):
		try:
			line[i] = float(line[i])
		except:
			pass

for line in tdata:
	if len(line) != 17:
		print line
"""
b = Methods(data)
b.KDM()
#c = Visualize(tdata)
#c.plot('age','BAC','female')
#c = Visualize(b.KDM())
#c.plot('age','BAC','female')
#a = Summary(data)
#a.mean()
#print a
#a.corr()