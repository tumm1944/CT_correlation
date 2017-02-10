#!/usr/bin/env python
# Filename: distribution_plot.py
#
import os, sys
import numpy as np
import glob
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from sklearn import linear_model
import scipy.stats as scs
from scipy.stats import kurtosis
import statsmodels.sandbox.distributions.extras as extras
import numpy.random as npr
import pylab
from operator import itemgetter
import seaborn as seas
from math import cos, sin
from datetime import datetime

def bootstrap(data, num_samples, statistic, alpha):
    """Returns bootstrap estimate of 100.0*(1-alpha) CI for statistic."""
#   print data
    n = len(data)
    print np.amin(data), np.amax(data), n
    data1 = np.zeros((n,),dtype=np.float)
    for i in range(n):
	data1[i]=data[i]
    idx = npr.randint(0, n, (num_samples, n))
    samples = data1[idx]
    stat = np.sort(statistic(samples, 1))
    print stat, len(stat)
    return (stat[int((alpha/2.0)*num_samples)],
            stat[int((1-alpha/2.0)*num_samples)])


def ransac_ct(X, y, pr_Xmin, pr_Xmax, Label):
# Fit line using all data
#   print X[1], y[1]
#   print pr_Xmin, pr_Xmax, float( (pr_Xmax-pr_Xmin)/len(X) )
    X = X.reshape(454, -1)
#   y = y.reshape(1, -1)
    print len(X), len(y)
    model = linear_model.LinearRegression()
    model.fit(X, y)
    print "Linear fit done"

# Robustly fit linear model with RANSAC algorithm
    model_ransac = linear_model.RANSACRegressor(linear_model.LinearRegression(), residual_threshold=0.5)
    model_ransac.fit(X, y)
    inlier_mask = model_ransac.inlier_mask_
    print model_ransac.residual_threshold
    outlier_mask = np.logical_not(inlier_mask)
    print "Ransac fit done"

# Predict data of estimated models
    line_X = np.arange(pr_Xmin, pr_Xmax, float( (pr_Xmax-pr_Xmin)/len(X) ) )
    line_y = model.predict(line_X[:, np.newaxis])
    line_y_ransac = model_ransac.predict(line_X[:, np.newaxis])

# Compare estimated coefficients
    print("Estimated coefficients (normal, RANSAC):")
    print( model.coef_, model_ransac.estimator_.coef_)

    lw = 2
    plt.scatter(X[inlier_mask], y[inlier_mask], color='yellowgreen', marker='.',label='Inliers')
    plt.scatter(X[outlier_mask], y[outlier_mask], color='gold', marker='.',label='Outliers')
    plt.plot(line_X, line_y, color='navy', linestyle='-', linewidth=lw,label='Linear regressor')
    plt.plot(line_X, line_y_ransac, color='cornflowerblue', linestyle='-',linewidth=lw, label='RANSAC regressor')
    plt.legend(loc='lower right')
    plt.ylabel ('CT state Energy (eV)')
    plt.xlabel (Label)
    plt.show()


def call_dplot(argv):
	FileCounter = 0
	myPath = os.path.dirname(os.path.realpath(__file__))
	if len(argv) < 1:
        	print "\nUsage: CT_correlation.py  <filename without extension> " 
        	sys.exit(1)
	else:
        	inFile = sys.argv[1]
		FileCounter = len(glob.glob1(myPath,inFile+'.txt'))
		if FileCounter < 1:
			print ' File not present '
			sys.exit(1)
		else:
			print '# of FIles', FileCounter
			dplot(inFile)
			sys.exit(1)

#	if not os.path.isfile('%s.custom2' %inFile):
#            print ' %s.custom2 DOES NOT EXIST !!!' %inFile
#	if not os.path.isfile('%s.custom3' %inFile):
#            print ' %s.custom3 DOES NOT EXIST !!!' %inFile
	

def dplot(extname):
	Bond_P = []
	Dihedral_P = []
	Angle_PC = []
	Dist_PC = []
	Charge_P = []
	E_CT = []
	X1 = []
	X2 = []
	X3 = []
	X4 = []
	Y1 = []
	Apc = [1,2,3,4,5]
	Dpc = [6,7,8,9,10]
	Cp = [11,12,13,14,15]
	Bp = [16,17,18,19]
	Dip = [20,21,22,23]
	count = 0
	myPath = os.path.dirname(os.path.realpath(__file__))

#	Filenames = glob.glob('%s/*.%s',myPath,extname)

	Filename = extname +'.txt'

	with open(Filename, 'r') as inF:
	    for _ in xrange(2):
		next(inF)
	    for line in inF:
			line = line.strip().split()
			E_CT.append(line[0])
			Angle_PC.append(itemgetter(*Apc)(line))
#			print Angle_PC[len(Angle_PC)-1]
			Dist_PC.append(itemgetter(*Dpc)(line))
			Charge_P.append(itemgetter(*Cp)(line))
			Bond_P.append(itemgetter(*Bp)(line))
			Dihedral_P.append(itemgetter(*Dip)(line))

	Total_pts = len(Angle_PC)
	Training_pts = int(0.6*Total_pts)
	Test_pts = Total_pts - Training_pts
	Max_Index =  np.argmax(Charge_P,axis=1)
#	print Charge_P[0][Max_Index[0]], Charge_P[0]
	Max_Index = np.array(Max_Index)
#	print	Max_Index
#	print len(Max_Index)
	for i in range(len(Max_Index)):
#		print i, Charge_P[i][Max_Index[i]]	
#		print Angle_PC[i][Max_Index[i]], float(Angle_PC[i][Max_Index[i]])/57.296
#		X1.append(cos(float(Angle_PC[i][Max_Index[i]])/57.296))
		Bond_index = (float(Bond_P[i][0])/1.45)*(float(Bond_P[i][1])/1.45)*(float(Bond_P[i][2])/1.45)*(float(Bond_P[i][3])/1.45)
		Deloc_index = ( float( Charge_P[i][0])/0.2)**2  + ( float( Charge_P[i][1])/0.2)**2 + ( float( Charge_P[i][2])/0.2)**2
		Deloc_index = Deloc_index + ( float( Charge_P[i][3])/0.2)**2 + ( float( Charge_P[i][4])/0.2)**2 
		Deloc_index = Deloc_index/25
#		X2.append(float(Bond_P[i][Max_Index[i]-1])/1.45)
		X2.append(float(Bond_index))
		Elec_stat = float( Charge_P[i][0] ) / float (Dist_PC[i][0] ) + float( Charge_P[i][1] ) / float (Dist_PC[i][1] ) 
		Elec_stat = float( Charge_P[i][3] ) / float (Dist_PC[i][3] ) + float( Charge_P[i][2] ) / float (Dist_PC[i][2] )
		Elec_stat = float( Charge_P[i][4] ) / float (Dist_PC[i][4] )
		X3.append (float(Elec_stat))
#		print float(Bond_P[i][Max_Index[i]-1])/1.45
#		X4.append(sin(float(Dihedral_P[i][Max_Index[i]-1])/57.296))
#		X1.append( float(Bond_index) / sin ( float(Dihedral_P[i][Max_Index[i]-1] ) / 57.296 ) )
		X1.append( float(Deloc_index) / float(Bond_index) )

#	X1 = Cosine (Angle)
#	X2 = Bond length normalized with 1.45 Angstroms
# 	X3 = Charge/Distance
#	X4 = Sine (Dihedral)
#	Y (Ect) is proportional to X1, X2, X3, X4
	X1   =  np.array(X1).astype(np.float)
	X2   =  np.array(X2).astype(np.float)
	X3   =  np.array(X3).astype(np.float)
	X4   =  np.array(X4).astype(np.float)

	E_CT =  np.array(E_CT).astype(np.float)

	ransac_ct(X1[:Training_pts],E_CT[:Training_pts],min(X1), max(X1), 'Deloc/Bond')	
	ransac_ct(X2[:Training_pts],E_CT[:Training_pts],min(X2), max(X2), 'Norm. Bond Index')
	ransac_ct(X3[:Training_pts],E_CT[:Training_pts],min(X3), max(X3), 'Electrostatic Interactions (e*e/Angstroms)')
	return (Filename)

if __name__ == '__main__':
    print '{:%Y-%m-%d %H:%M}'.format(datetime.now())
    x = call_dplot(sys.argv[1:])
   # data of interest is bimodal and obviously not normal
#   x = np.concatenate([npr.normal(3, 1, 100), npr.normal(6, 2, 200)])

    # find mean 95% CI and 100,000 bootstrap samples
#    low, high = bootstrap(x, 100000, np.mean, 0.01)
#    print "mean = ", low, high
#    low, high = bootstrap(x, 100000, np.std, 0.01)
#    print "std = ", low, high


    # make plots
#    pylab.figure(figsize=(8,4))
#    pylab.subplot(121)
#    pylab.hist(x, 50, histtype='step')
#    pylab.title('Historgram of data')
#    pylab.subplot(122)
#    pylab.plot([-0.03,0.03], [np.mean(x), np.mean(x)], 'r', linewidth=2)
#    pylab.scatter(0.12+(npr.random(len(x))*1.8), x)
#    pylab.plot([0.19,0.21], [low, low], 'r', linewidth=2)
#    pylab.plot([0.19,0.21], [high, high], 'r', linewidth=2)
#    pylab.plot([0.2,0.2], [low, high], 'r', linewidth=2)
#    pylab.xlim([0, 2])
#    pylab.title('Bootstrap 95% CI for mean')
#    pylab.savefig('boostrap.png')

#####  End  #####
