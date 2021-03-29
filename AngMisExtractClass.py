#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 31 11:38:37 2018

@author: brentbelland
"""

from math import acos,sin,cos
#import sys
import numpy as np
#import matplotlib.pyplot as plt #Trying not to use the plt in the class function
#plt.rc(usetex = True)
import scipy.optimize as opt


class AngMisExtraction():
    
    def __init__(self,FRDlist,FRDerrlist,thetalist,philist):
        """Generates an Angular Misalignment Extraction Object"""
        
        #Generate a function to accept inputs other than self
        self.FRDlist = FRDlist
        self.FRDerrlist = FRDerrlist
        self.thetalist = thetalist
        self.philist = philist
        
        #Run checks to verify the format of FRDlist, thetalist, philist
        
        self.FRDcheck()
        self.thetacheck()
        self.phicheck()
        
        #Then test for angular misalignment.
        
    def FRDcheck(self):
        
        if not (isinstance(self.FRDlist,list) or isinstance(self.FRDlist,np.ndarray)):
            #Error
            raise Exception('FRDlist should be a np.ndarray or a list. FRDlist\'s current type is {}'.format(type(self.FRDlist))
            
        elif isinstance(self.FRDlist,list):
            self.FRDlist = np.array(self.FRDlist)
                            
    def FRDerrcheck(self):
        
        if not (isinstance(self.FRDerrlist,list) or isinstance(self.FRDerrlist,np.ndarray)):
            #Error
            raise Exception('FRDlist should be a np.ndarray or a list. FRDlist\'s current type is {}'.format(type(self.FRDlist))            
            
        elif isinstance(self.FRDlist,list):
            self.FRDlist = np.array(self.FRDlist)
                            
                            
        if len(FRDerrlist) != len(FRDlist):
            raise Exception('FRDerrlist should be the same length as the FRDlist')
        
    
    def thetacheck(self):
            
        if not (isinstance(self.thetalist,list) or isinstance(self.thetalist,np.ndarray) or isinstance(self.thetalist,np.ndarray)):
            #Error
            raise Exception('thetalist should be a np.ndarray, a list, or a float. thetalist\'s current type is {}'.format(type(self.thetalist))
            
        elif isinstance(self.thetalist,np.float):
            #Need to set thetalist as a constant list of length FRDlist if only one value input
            self.thetalist = self.thetalist*np.array(np.ones(np.shape(self.FRDlist)))
            
        elif isinstance(self.thetalist,list):
            self.thetalist = np.array(self.thetalist)
            
    def phicheck(self):
           
        if not (isinstance(self.philist,list) or isinstance(self.philist,np.ndarray) or isinstance(self.philist,np.ndarray)):
            #Error
            raise Exception('philist should be a np.ndarray, a list, or a float. philist\'s current type is {}'.format(type(self.philist))
            
        elif isinstance(self.philist,np.float):
            #Need to set philist as a constant list of length philist if only one value input
            self.philist = self.philist*np.array(np.ones(np.shape(self.philist)))
            
        elif isinstance(self.philist,list):
            self.philist = np.array(self.philist)
                            
    def netAngularMisalignment(Alpha,Beta,Gamma,Theta,Phi):
        #This function combines the known Theta, Phi with unknown initial angular offsets Alpha, Beta, Gamma in the Cobra to output a total angular misalignment. The purpose of this function is to utilize the additional information gained from FRD information to extract out the angular component. Assumes that theta and phi (and alpha, beta, gamma) are in units of radians.
        value = ((cos(c) - 1)*(sin(a)*sin(p)*cos(b) - sin(a)*sin(t)*cos(b)*cos(p - t) + sin(a)*sin(t)*cos(p - t) + sin(b)*sin(p - t)*cos(a))**2 - cos(c))*(sin(a)*sin(b)*cos(t) - cos(a)*cos(b)) - (sin(a)*cos(b) + sin(b)*cos(a)*cos(t))*((cos(c) - 1)*(sin(a)*sin(b)*sin(p - t) - sin(p)*cos(a)*cos(b) + sin(t)*cos(a)*cos(b)*cos(p - t) - sin(t)*cos(a)*cos(p - t))*(sin(a)*sin(p)*cos(b) - sin(a)*sin(t)*cos(b)*cos(p - t) + sin(a)*sin(t)*cos(p - t) + sin(b)*sin(p - t)*cos(a)) + (-sin(t)*sin(p - t)*cos(b) + sin(t)*sin(p - t) + cos(p))*sin(c)) - ((cos(c) - 1)*(-sin(t)*sin(p - t)*cos(b) + sin(t)*sin(p - t) + cos(p))*(sin(a)*sin(p)*cos(b) - sin(a)*sin(t)*cos(b)*cos(p - t) + sin(a)*sin(t)*cos(p - t) + sin(b)*sin(p - t)*cos(a)) - (sin(a)*sin(b)*sin(p - t) - sin(p)*cos(a)*cos(b) + sin(t)*cos(a)*cos(b)*cos(p - t) - sin(t)*cos(a)*cos(p - t))*sin(c))*sin(b)*sin(t)
        if value > 1:
            return 0
            #raise Exception('Error: Net angle found is not valid (Cos(angle) > 1)
        else:
            return acos(value)        
        
    def plottheta(a,b,c,tpoff,tphase,p):
        #This function's purpose is just to generate a nice plot of angular misalignments across theta inputs.
        thetalist = np.zeros([51,1])
        anglelist = np.zeros([51,1])
        for i in range(0,51):
            thetalist[i] = 2*np.pi*i/50
            anglelist[i] = netAngularMisalignment(a,b,c,2*np.pi*i/50+tphase,p+np.pi+tpoff+2*np.pi*i/50)
        return thetalist,anglelist

    def ResidualSum(abctotpo,args):
        (alpha, beta, gamma, thetaoffset, thetaphioffset) = abctotpo
        (angledata,UpperUncert,LowerUncert) = args
        """The sum of squares of residuals that scipy.optimize.minimize can iterate over"""
                            
        #angledata = np.array([[0,0,0,0,0,0,0,0,0,0,0],[0,6.3e-3,8.4e-3,9.0e-3,8.1e-3,6.9e-3,8.3e-3,9.4e-3,8.5e-3,7.3e-3,4.4e-3],[9.2e-3,10.4e-3,11.6e-3,11.4e-3,9.3e-3,9.2e-3,9.5e-3,11.7e-3,12.0e-3,10.8e-3,9.5e-3]])
        #UpperUncert = np.array([[2e-3,2e-3,2e-3,2e-3,2e-3,2e-3,2e-3,2e-3,2e-3,2e-3,2e-3],[10e-3,0.5679e-3,0.1635e-3,0.8141e-3,0.416e-3,0.3885e-3,0.6875e-3,0.5892e-3,0.1469e-3,0.7304e-3,10e-3],[0.1967e-3,0.4904e-3,0.3926e-3,0.1661e-3,0.7875e-3,0.2230e-3,0.1856e-3,0.4303e-3,0.9881e-3,0.2555e-3,0.7544e-3]])#For when the value is HIGHER than the estimated value (asymmetric error bars)
        #LowerUncert = np.array([[2e-3,2e-3,2e-3,2e-3,2e-3,2e-3,2e-3,2e-3,2e-3,2e-3,2e-3],[3e-3,0.1e-3,0.9e-3,0.1e-3,0.4e-3,1.1e-3,0.2e-3,1.5e-3,0.8e-3,1.5e-3,2e-3],[0.1967e-3,0.4904e-3,0.3926e-3,0.1661e-3,0.7875e-3,0.2230e-3,0.1856e-3,0.4303e-3,0.9881e-3,0.2555e-3,0.7544e-3]])#For when the value is LOWER than the estimated value
    
        squareresiduals = 0
        for theta_tot in range(0,11):#(0, 2*np.pi/5, 4*np.pi/5,6*np.pi/5,8*np.pi/5,2*np.pi):
            theta = 2*np.pi-np.abs(2*np.pi - 2*np.pi/5*theta_tot)
            for phi in range(0,3):#np.pi/2,np.pi):
                #Must split this into two cases: One where the the model overestimates the value, and one where the model underestimates the value
                residual = (value(alpha,beta,gamma,theta+thetaoffset,np.pi/2*phi+theta+np.pi+thetaphioffset)-angledata[phi][theta_tot])
            #index = theta+phi*11
                if residual > 0:
                    squareresiduals += residual**2/((UpperUncert[phi][theta_tot])**2)
                else:
                    squareresiduals += residual**2/((LowerUncert[phi][theta_tot])**2)


    def optimizeangles(guessvalues):
        #guessvalues should be in the form [alphaguess, betaguess, gammaguess, thetaoffset, netphioffset], all in radians
        #args should be a tuple with angledata (Which is really what I'm trying to solve for...), UpperUncert (on the angle) and LowerUncert (also on the angle).
        bestanglesResult = opt.minimize(self.ResidualSum,guessvalues,args=(),method='Nelder-Mead') #change to method = 'Nelder-Mead' instead of default BFGS?

        bestangles = bestanglesResult.x
        print('This truth value of this result being successful is {}'.format(bestanglesResult.success))
        print('Additional message(s): {}'.format(bestanglesResult.message))

                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            

#a = float(sys.argv[1])
#b = float(sys.argv[2])
#c = float(sys.argv[3])
#t = float(sys.argv[4])
#p = float(sys.argv[5])# + t

def value(a,b,c,t,p):
    value = ((cos(c) - 1)*(sin(a)*sin(p)*cos(b) - sin(a)*sin(t)*cos(b)*cos(p - t) + sin(a)*sin(t)*cos(p - t) + sin(b)*sin(p - t)*cos(a))**2 - cos(c))*(sin(a)*sin(b)*cos(t) - cos(a)*cos(b)) - (sin(a)*cos(b) + sin(b)*cos(a)*cos(t))*((cos(c) - 1)*(sin(a)*sin(b)*sin(p - t) - sin(p)*cos(a)*cos(b) + sin(t)*cos(a)*cos(b)*cos(p - t) - sin(t)*cos(a)*cos(p - t))*(sin(a)*sin(p)*cos(b) - sin(a)*sin(t)*cos(b)*cos(p - t) + sin(a)*sin(t)*cos(p - t) + sin(b)*sin(p - t)*cos(a)) + (-sin(t)*sin(p - t)*cos(b) + sin(t)*sin(p - t) + cos(p))*sin(c)) - ((cos(c) - 1)*(-sin(t)*sin(p - t)*cos(b) + sin(t)*sin(p - t) + cos(p))*(sin(a)*sin(p)*cos(b) - sin(a)*sin(t)*cos(b)*cos(p - t) + sin(a)*sin(t)*cos(p - t) + sin(b)*sin(p - t)*cos(a)) - (sin(a)*sin(b)*sin(p - t) - sin(p)*cos(a)*cos(b) + sin(t)*cos(a)*cos(b)*cos(p - t) - sin(t)*cos(a)*cos(p - t))*sin(c))*sin(b)*sin(t)
    if value > 1:
        #print(str.format('{0:.15f}',value))
        return 0
    else:
        return acos(value)
    
def value2(a,b,c,t,p): #This maps p->p+t, since rotation of the theta stage rotates the net angle that the phi stage rotates around.
    value = ((cos(c) - 1)*(sin(a)*sin(p+t)*cos(b) - sin(a)*sin(t)*cos(b)*cos(p) + sin(a)*sin(t)*cos(p) + sin(b)*sin(p)*cos(a))**2 - cos(c))*(sin(a)*sin(b)*cos(t) - cos(a)*cos(b)) - (sin(a)*cos(b) + sin(b)*cos(a)*cos(t))*((cos(c) - 1)*(sin(a)*sin(b)*sin(p) - sin(p+t)*cos(a)*cos(b) + sin(t)*cos(a)*cos(b)*cos(p) - sin(t)*cos(a)*cos(p))*(sin(a)*sin(p+t)*cos(b) - sin(a)*sin(t)*cos(b)*cos(p) + sin(a)*sin(t)*cos(p) + sin(b)*sin(p)*cos(a)) + (-sin(t)*sin(p)*cos(b) + sin(t)*sin(p) + cos(p+t))*sin(c)) - ((cos(c) - 1)*(-sin(t)*sin(p)*cos(b) + sin(t)*sin(p) + cos(p+t))*(sin(a)*sin(p+t)*cos(b) - sin(a)*sin(t)*cos(b)*cos(p) + sin(a)*sin(t)*cos(p) + sin(b)*sin(p)*cos(a)) - (sin(a)*sin(b)*sin(p) - sin(p+t)*cos(a)*cos(b) + sin(t)*cos(a)*cos(b)*cos(p) - sin(t)*cos(a)*cos(p))*sin(c))*sin(b)*sin(t)
    if value > 1:
        #print(str.format('{0:.15f}',value))
        return 0
    else:
        return acos(value)    
    
def nicetheta(a,b,c,tpoff,tphase,p):
    thetalist = np.zeros([51,1])
    anglelist = np.zeros([51,1])
    for i in range(0,51):
        thetalist[i] = 2*np.pi*i/50
        anglelist[i] = value(a,b,c,2*np.pi*i/50+tphase,p+np.pi+tpoff+2*np.pi*i/50)
    return thetalist,anglelist
                    

#value = acos(((cos(c) - 1)*(sin(a)*sin(p)*cos(b) - sin(a)*sin(t)*cos(b)*cos(p - t) + sin(a)*sin(t)*cos(p - t) + sin(b)*sin(p - t)*cos(a))**2 - cos(c))*(sin(a)*sin(b)*cos(t) - cos(a)*cos(b)) - (sin(a)*cos(b) + sin(b)*cos(a)*cos(t))*((cos(c) - 1)*(sin(a)*sin(b)*sin(p - t) - sin(p)*cos(a)*cos(b) + sin(t)*cos(a)*cos(b)*cos(p - t) - sin(t)*cos(a)*cos(p - t))*(sin(a)*sin(p)*cos(b) - sin(a)*sin(t)*cos(b)*cos(p - t) + sin(a)*sin(t)*cos(p - t) + sin(b)*sin(p - t)*cos(a)) + (-sin(t)*sin(p - t)*cos(b) + sin(t)*sin(p - t) + cos(p))*sin(c)) - ((cos(c) - 1)*(-sin(t)*sin(p - t)*cos(b) + sin(t)*sin(p - t) + cos(p))*(sin(a)*sin(p)*cos(b) - sin(a)*sin(t)*cos(b)*cos(p - t) + sin(a)*sin(t)*cos(p - t) + sin(b)*sin(p - t)*cos(a)) - (sin(a)*sin(b)*sin(p - t) - sin(p)*cos(a)*cos(b) + sin(t)*cos(a)*cos(b)*cos(p - t) - sin(t)*cos(a)*cos(p - t))*sin(c))*sin(b)*sin(t))


#print(value)

"""Phi = 0: Theta [0:700:3500, 2800:-700:0], angle [0,0,0,0,0,0,0,0,0,0,0]
Phi = pi/2: Theta [0:700:3500, 2800:-700:0], angle [0,6.3,8.4,9.0,8.1,6.9,8.3,9.4,8.5,7.3,4.4] *Not good to trust the 0 degree data*
Phi = pi: Theta [0:700:3500, 2800:-700:0], angle [9.2,10.4,11.6,9.3,9.2,9.5,11.7,12.0,10.8,9.5]"""

#angledata = np.array([[0,0,0,0,0,0],[4.4e-3,6.8e-3,8.45e-3,9.2e-3,8.2e-3,6.9e-3],[9.35e-3,10.6e-3,11.8e-3,11.7e-3,9.4e-3,9.35e-3]])
#NomUncert = np.array([[2e-3,2e-3,2e-3,2e-3,2e-3,2e-3],[3e-3,1e-3,1e-3,1e-3,1e-3,1e-3],[1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]])
#Uncertain = np.array([[2e-3,2e-3,2e-3,2e-3,2e-3,2e-3],[0.3e-3,0.7e-3,0.5e-3,0.8e-3,0.4e-3,0.7e-3],[0.3e-3,0.3e-3,0.6e-3,0.3e-3,0.465e-3,0.2e-3]])

#From top to bottom, phi = 0,pi/2,pi; from left to right, theta = 0,2pi/5,4pi/5,6pi/5,8pi/5,2pi

def ResidualSum(abctotpo):
    (alpha, beta, gamma, thetaoffset, thetaphioffset) = abctotpo
    """The sum of squares of residuals that scipy.optimize.minimize can iterate over"""
    angledata = np.array([[0,0,0,0,0,0,0,0,0,0,0],[0,6.3e-3,8.4e-3,9.0e-3,8.1e-3,6.9e-3,8.3e-3,9.4e-3,8.5e-3,7.3e-3,4.4e-3],[9.2e-3,10.4e-3,11.6e-3,11.4e-3,9.3e-3,9.2e-3,9.5e-3,11.7e-3,12.0e-3,10.8e-3,9.5e-3]])
    UpperUncert = np.array([[2e-3,2e-3,2e-3,2e-3,2e-3,2e-3,2e-3,2e-3,2e-3,2e-3,2e-3],[10e-3,0.5679e-3,0.1635e-3,0.8141e-3,0.416e-3,0.3885e-3,0.6875e-3,0.5892e-3,0.1469e-3,0.7304e-3,10e-3],[0.1967e-3,0.4904e-3,0.3926e-3,0.1661e-3,0.7875e-3,0.2230e-3,0.1856e-3,0.4303e-3,0.9881e-3,0.2555e-3,0.7544e-3]])#For when the value is HIGHER than the estimated value (asymmetric error bars)
    LowerUncert = np.array([[2e-3,2e-3,2e-3,2e-3,2e-3,2e-3,2e-3,2e-3,2e-3,2e-3,2e-3],[3e-3,0.1e-3,0.9e-3,0.1e-3,0.4e-3,1.1e-3,0.2e-3,1.5e-3,0.8e-3,1.5e-3,2e-3],[0.1967e-3,0.4904e-3,0.3926e-3,0.1661e-3,0.7875e-3,0.2230e-3,0.1856e-3,0.4303e-3,0.9881e-3,0.2555e-3,0.7544e-3]])#For when the value is LOWER than the estimated value
    
    squareresiduals = 0
    for theta_tot in range(0,11):#(0, 2*np.pi/5, 4*np.pi/5,6*np.pi/5,8*np.pi/5,2*np.pi):
        theta = 2*np.pi-np.abs(2*np.pi - 2*np.pi/5*theta_tot)
        for phi in range(0,3):#np.pi/2,np.pi):
            #Must split this into two cases: One where the the model overestimates the value, and one where the model underestimates the value
            residual = (value(alpha,beta,gamma,theta+thetaoffset,np.pi/2*phi+theta+np.pi+thetaphioffset)-angledata[phi][theta_tot])
            #index = theta+phi*11
            if residual > 0:
                squareresiduals += residual**2/((UpperUncert[phi][theta_tot])**2)
            else:
                squareresiduals += residual**2/((LowerUncert[phi][theta_tot])**2)


    return squareresiduals

bestanglesResult = opt.minimize(ResidualSum,[1e-3,5e-3,5e-3,0,0],method='Nelder-Mead') #change to method = 'Nelder-Mead' instead of default BFGS?

bestangles = bestanglesResult.x
print('This truth value of this result being successful is {}'.format(bestanglesResult.success))
print('Additional message(s): {}'.format(bestanglesResult.message))

#bestangles = (besta, bestb, bestc, (bestthetaoff, bestphioff))

#alphas = range(10,13,1)
#alphas = [i * 1e-4 for i in alphas]
#betas = range(-56,-53,1)
#betas = [i * 1e-4 for i in betas]
#gammas = range(-55,-44,1)
#gammas = [i * 1e-4 for i in gammas]
#thetaphioffsets = range(-20,21,1)
#thetaphioffsets = [i * 1e-3 for i in gammas]
#thetaphases = range(-20,21,1)
#thetaphases = [i * 1e-3 for i in gammas]
#
#bestresiduals = 999
#bestangles = (np.nan,np.nan,np.nan)
#for alpha in alphas:
#    for beta in betas:
#        for gamma in gammas:
#            for thetaphioffset in thetaphioffsets:
#                for thetaphase in thetaphases:
#                    squareresiduals = 0
#                    for theta in range(0,6):#(0, 2*np.pi/5, 4*np.pi/5,6*np.pi/5,8*np.pi/5,2*np.pi):
#                        for phi in range(0,3):#np.pi/2,np.pi):
#                            #print((alpha,beta,gamma))
#                            squareresiduals += ((value(alpha,beta,gamma,2*np.pi/5*theta+thetaphase,np.pi/2*phi+2*np.pi/5*theta+np.pi+thetaphioffset)-angledata[phi][theta])/(Uncertain[phi][theta]))**2
#                    if squareresiduals < bestresiduals:
#                        #print('Improved fit with a={}, b={}, c={}'.format(alpha,beta,gamma))
#                        bestresiduals = squareresiduals
#                        bestangles = (alpha,beta,gamma,thetaphase,thetaphioffset)

#opt.minimize(value,[])


print("The best group of (alpha, beta, gamma,offset_theta, offset_thetaphi) are {}.".format(bestangles))

angledata = np.array([[0,0,0,0,0,0,0,0,0,0,0],[0,6.3e-3,8.4e-3,9.0e-3,8.1e-3,6.9e-3,8.3e-3,9.4e-3,8.5e-3,7.3e-3,4.4e-3],[9.2e-3,10.4e-3,11.6e-3,11.4e-3,9.3e-3,9.2e-3,9.5e-3,11.7e-3,12.0e-3,10.8e-3,9.5e-3]])
UpperUncert = np.array([[2e-3,2e-3,2e-3,2e-3,2e-3,2e-3,2e-3,2e-3,2e-3,2e-3,2e-3],[10e-3,0.5679e-3,0.1635e-3,0.8141e-3,0.416e-3,0.3885e-3,0.6875e-3,0.5892e-3,0.1469e-3,0.7304e-3,10e-3],[0.1967e-3,0.4904e-3,0.3926e-3,0.1661e-3,0.7875e-3,0.2230e-3,0.1856e-3,0.4303e-3,0.9881e-3,0.2555e-3,0.7544e-3]])
LowerUncert = np.array([[2e-3,2e-3,2e-3,2e-3,2e-3,2e-3,2e-3,2e-3,2e-3,2e-3,2e-3],[10e-3,0.1e-3,0.9e-3,0.1e-3,0.4e-3,1.1e-3,0.2e-3,1.5e-3,0.8e-3,1.5e-3,10e-3],[0.1967e-3,0.4904e-3,0.3926e-3,0.1661e-3,0.7875e-3,0.2230e-3,0.1856e-3,0.4303e-3,0.9881e-3,0.2555e-3,0.7544e-3]])
simppos = np.array([0,2*np.pi/5,4*np.pi/5,6*np.pi/5,8*np.pi/5,2*np.pi,8*np.pi/5,6*np.pi/5,4*np.pi/5,2*np.pi/5,0])

plt.figure()
plt.errorbar(simppos,np.array(angledata[0]),yerr=[LowerUncert[0], UpperUncert[0]],xerr=np.array([np.pi/36,np.pi/36,np.pi/36,np.pi/36,np.pi/36,np.pi/36,np.pi/36,np.pi/36,np.pi/36,np.pi/36,np.pi/36]),fmt='b.')
plt.errorbar(simppos,np.array(angledata[1]),yerr=[LowerUncert[1], UpperUncert[1]],xerr=[np.pi/36,np.pi/36,np.pi/36,np.pi/36,np.pi/36,np.pi/36,np.pi/36,np.pi/36,np.pi/36,np.pi/36,np.pi/36],fmt='r.')
plt.errorbar(simppos,np.array(angledata[2]),yerr=[LowerUncert[2], UpperUncert[2]],xerr=[np.pi/36,np.pi/36,np.pi/36,np.pi/36,np.pi/36,np.pi/36,np.pi/36,np.pi/36,np.pi/36,np.pi/36,np.pi/36],fmt='k.')
thetas_philow, angles_philow = nicetheta(bestangles[0],bestangles[1],bestangles[2],bestangles[3],bestangles[4],0)
thetas_phimid, angles_phimid = nicetheta(bestangles[0],bestangles[1],bestangles[2],bestangles[3],bestangles[4],np.pi/2)
thetas_phihig, angles_phihig = nicetheta(bestangles[0],bestangles[1],bestangles[2],bestangles[3],bestangles[4],np.pi)
plt.plot(thetas_philow,angles_philow,'b-',lw=1)
plt.plot(thetas_phimid,angles_phimid,'r-',lw=1)
plt.plot(thetas_phihig,angles_phihig,'k-',lw=1)
plt.autoscale(enable=True, axis=u'both', tight=False)
plt.xlim(0,6.3)
plt.ylim(0,13e-3)
plt.title('Misalignment 1={0:.5f}mrad, Misalignment 2={1:.5f}mrad, Misalignment3={2:.5f}mrad, $\theta_0$ = {2:.5f}, $\phi_0$ = {3:.2f}'.format(bestangles[0],bestangles[1],bestangles[2],bestangles[3],bestangles[4]))
plt.savefig('Anglefitting.png')
plt.figure()
plt.title('Residuals')

print('Datapoints\n')
print(np.array(angledata[0]))
print(np.array(angledata[1]))
print(np.array(angledata[2]))

print('Errors (Low then high alternating) \n')
print(LowerUncert[0])
print(UpperUncert[0])
print(LowerUncert[1])
print(UpperUncert[1])
print(LowerUncert[2])
print(UpperUncert[2])

print('Fit\n')
print(angles_philow)
print(angles_phimid)
print(angles_phihig)


def valuegive(bestangledata,positions,phi):
    returnvector = np.zeros((1,len(positions)))[0]
    for position in range(len(positions)):
        returnvector[position] = value(bestangledata[0],bestangledata[1],bestangledata[2],positions[position]+bestangledata[3],positions[position]+phi+np.pi+bestangledata[4])
    return returnvector



plt.errorbar(simppos,np.array(angledata[0])-valuegive(bestangles,simppos,0),yerr=[LowerUncert[0], UpperUncert[0]],fmt='b.')
plt.errorbar(simppos,np.array(angledata[1])-valuegive(bestangles,simppos,np.pi/2),yerr=[LowerUncert[1], UpperUncert[1]],fmt='r.')
plt.errorbar(simppos,np.array(angledata[2])-valuegive(bestangles,simppos,np.pi),yerr=[LowerUncert[2], UpperUncert[2]],fmt='k.')
plt.ylim(-1e-3,1e-3)
plt.show()

