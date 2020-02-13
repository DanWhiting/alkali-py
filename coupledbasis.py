from __future__ import division
''' 
Author - Daniel J. Whiting 
Date modified - 10/08/2017

Contains functions for calculating Clebsh-Gordon coefficients and dipole
matrix elements of electric dipole transitions in the fully coupled basis
which is useful in the regime of weak magnetic fields.

WARNING: This code is currently untested for accuracy
'''

from uncoupledbasis import *

def CoupledEvals(I,L,S,J,B):
	out = np.zeros((int((2*J+1)*(2*I+1)),3))
	evals = EvalsEvecs(I,L,S,J,B)[0]
	dF = int((2*I+1)*(2*L+1)*(2*S+1))
	i = 0
	for F in np.arange(abs(J-I),abs(J+I)+1):
		for mF in np.arange(-F,F+1):
			index = int(1e-10 + (F*(F+1)-I*(I+1)-J*(J+1)) + (2*I+1)*(J*(J+1)-L*(L+1)-S*(S+1)) + mF + (dF+1)/2)
			out[i,:] = [F,mF,evals[index-1]/1e9] ###### [mL[mS[mI]]]
			i += 1
	return out

def CoupledEvecs(I,L,S,J,B):
	coupledevecs = []
	evecs = EvalsEvecs(I,L,S,J,B)[1]
	dF = int((2*I+1)*(2*L+1)*(2*S+1))
	for F in np.arange(abs(J-I),abs(J+I)+1):
		for mF in np.arange(-F,F+1):
			index = int(1e-10 + F*(F+1) - I*(I+1) - J*(J+1) + (2*I+1)*(J*(J+1)-L*(L+1)-S*(S+1)) + mF + (dF+1)/2)
			coupledevecs.append(evecs[:,index-1]) ###### [mL[mS[mI]]]
	return coupledevecs

def CoupledCGs(I,S,L1,J1,F1,mF1,L2,J2,F2,q):
	if abs(mF1+q)>F2:
		cg=0
	else:
		gvec = CoupledEvecs(I,L1,S,J1,F1,mF1,B)
		lg = len(gvec)
		evec = CoupledEvecs(I,L2,S,J2,F2,mF1+q,B)
		le = len(evec)
		lm = (le-lg)/2
		lp = (le+lg)/2
		cg = dot(gvec,evec[lm*(1-q):lp-lm*q])
	return cg

def CoupledDipole(I,S,L1,J1,F1,mF1,L2,J2,F2,q,B):
	'''
	returns relative dipole moment of entered transition
	i.e. the dipole moment in units of the (J -> J') dipole moment
	'''
	cg = get_cg(I,S,L1,J1,F1,mF1,L2,J2,F2,q,B)
	rdip = cg.real*((2*J1+1)/(2*J2+1))**.5 # state manifold reduction
	return rdip

if __name__ == "__main__":
	I = 3/2
	L = 0
	S = 1/2
	J = 1/2
	B = 0.62 # Tesla
	
	print CoupledEvals(I,L,S,J,B)
