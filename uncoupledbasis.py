'''
Author - Daniel J. Whiting 
Date modified - 10/08/2017

Contains functions for performing atomic structure calculations in the
fully uncoupled basis.
'''

from __future__ import division
from numpy import matrix,arange,zeros,kron,eye,pi,dot,linspace,exp,array
from scipy.linalg import eigh

import matplotlib
import matplotlib.pyplot as plt
import fractions
import atomdata
import numpy as np

from constants import *
h = 2*pi*hbar

def jmat(j):
	''' calculates the matrix J^2 from the angular momentum vector j '''
	d = int(2*j+1+1e-15)
	m = arange(j,-j-1,-1)
	jm = (j*(j+1)-m*(m+1))**.5
	Jp = matrix(zeros((d,d)))
	Jp[range(0,d-1),range(1,d)] = jm[1:]
	Jm = Jp.transpose()
	Jx = 1/2*(Jp+Jm)
	Jy = -1j/2*(Jp-Jm)
	Jz = 1/2*(Jp*Jm-Jm*Jp)
	JJ = Jx,Jy,Jz
	return JJ

def jsum(j1,j2):
	''' calculates the angular momentum vector J^2 from the sum
	of the angular momentum vectors j1 and j2 '''
	d1 = j1[0].shape[0]
	d2 = j2[0].shape[0]
	Jx = kron(j1[0],eye(d2))+kron(eye(d1),j2[0])
	Jy = kron(j1[1],eye(d2))+kron(eye(d1),j2[1])
	Jz = kron(j1[2],eye(d2))+kron(eye(d1),j2[2])
	JJ = Jx,Jy,Jz
	return JJ

def EvalsEvecs(I,L,S,J,B=0,IsotopeShift=0):
	''' Calculates the eigenvalues and eigenvectors of the atomic
	structure Hamiltonian. Including spin-orbit, nuclear, magnetic dipole, 
	electric quadrupole and Zeeman interactions.
	Calculations and results are in the fully uncoupled basis.
	I,L,S,J are the spin anglar momenta of the state of interest.
	B is the magnetic field strength in Tesla.
	IsotopeShift is a global shift to the resulting eigenvalues.''' 
	SS,LL,II = jmat(S),jmat(L),jmat(I)

	### The order in which ang mom's are combined here defines the order
	### of the decoupled angular momenta in the eigenvectors.
	### Current order gives [mL[mS[mI]]]
	JJ = jsum(LL,SS)
	FF = jsum(JJ,II)

	S2 = SS[0]**2+SS[1]**2+SS[2]**2
	L2 = LL[0]**2+LL[1]**2+LL[2]**2
	J2 = JJ[0]**2+JJ[1]**2+JJ[2]**2
	I2 = II[0]**2+II[1]**2+II[2]**2
	F2 = FF[0]**2+FF[1]**2+FF[2]**2

	dS = SS[0].shape[0]
	dL = LL[0].shape[0]
	dJ = JJ[0].shape[0]
	dI = II[0].shape[0]
	dF = FF[0].shape[0]

	#### Set up the hamiltonian in the uncoupled (S,L,I) basis ####
	H = 0
	
	A_fs, A_hfs, B_hfs, gI, gL = atomdata.atomic_structure_coefficients('Rb',I,L,J)
	
	# Fine structure
	LdotS = 1/2*kron(J2 - (kron(L2,eye(dS)) + kron(eye(dL),S2)),eye(dI))
	
	# Spin orbit interaction
	# recenter the appropriate J state on zero energy
	so_correction_factor = -(J*(J+1)-L*(L+1)-S*(S+1))/2 
	H += A_fs * (LdotS + so_correction_factor*eye(dF))

	# Hyperfine structure
	IdotJ = 1/2*(F2 - (kron(eye(dJ),I2) + kron(J2,eye(dI))))
	
	# Magnetic dipole interaction
	H += A_hfs * IdotJ
	
	# Electric quadrupole interaction
	if B_hfs != 0:
		H += B_hfs * (3*IdotJ**2 + 3/2*IdotJ - I*(I+1)*J*(J+1)*eye(dF))/(2*I*(2*I-1)*J*(2*J-1))

	# Zeeman interaction
	H -= muB*B/h * (gL*kron(kron(LL[2],eye(dS)),eye(dI)) + gS*kron(kron(eye(dL),SS[2]),eye(dI)) + gI*kron(kron(eye(dL),eye(dS)),II[2]))
	
	# Isotope shift
	H += IsotopeShift * eye(dF)

	return eigh(H) # return the eigenvalues and eigenvectors of the hamiltonian

def BreitRabi(I,L,S,J,Bmax=0.62,ylim=None):
	''' Calculates the Breit Rabi diagram for the evolution of the
	eigenvalues under an applied magnetic field.
	I,L,S,J are the spin anglar momenta of the state of interest.
	Bmax is the maximum magnetic field in Tesla.
	ylim (-ylim) is the maximum (minimum) energy/h displayed. '''
	Bs = linspace(0,Bmax,1e3)
	d = int((2*I+1)*(2*L+1)*(2*S+1))
	Es = zeros((len(Bs),d))
	for i in range(len(Bs)):
		Es[i,:] = EvalsEvecs(I,L,S,J,Bs[i])[0]
	fig = plt.figure(facecolor='white')
	ax = fig.add_subplot(111)
	ax.plot(Bs,Es/1e9,c='k',lw=1)
	ax.set_xlabel('Magnetic field (T)')
	ax.set_ylabel('Energy / h (GHz)')
	ax.set_xlim(0,Bmax)
	if ylim != None:
		ax.set_ylim(-ylim,ylim)
	plt.tight_layout()

def MomentumDecomp(I,L,S,J,B=0,ylim=None,cutoff=0.001):
	''' Calculates the uncoupled basis composition of a state defined
	by spin anglar momenta I,L,S,J.
	B is the magnetic field strength in Tesla.
	ylim (-ylim) is the maximum (minimum) energy/h displayed.
	cutoff is the minimum size of a component of a state that
	is considered significant enough to display. '''
	if ylim == None:
		ylim = np.inf
	evalsevecs = EvalsEvecs(I,L,S,J,B)
	evecs = abs(evalsevecs[1])
	evals = evalsevecs[0]/1e9
	evecs,evals=evecs[:,::-1],evals[::-1]						# reverse the array so highest energy first
	print 'Decoupled quantum numbers with fractions >', cutoff
	for i in range(evecs.shape[1]):
		if abs(evals[i])<ylim:
			print i
			evec = evecs[:,i]
			print '----------------------'
			print 'Energy:', evals[i],'GHz'
			print '       |mI,  mL,  mS >'
			for iI,mI in enumerate(arange(-I,I+1)):
				for iL,mL in enumerate(arange(-L,L+1)):
					for iS,mS in enumerate(arange(-S,S+1)):
						frac = evec[int(iL*(2*I+1)*(2*S+1)+iS*(2*I+1)+iI)]
						if frac > cutoff:
							print '%.4f' % frac, '|%2s,%4s,%4s >' % (str(fractions.Fraction(mI)),str(fractions.Fraction(mL)),str(fractions.Fraction(mS)))

if __name__=="__main__":
	S = 1/2
	I = 3/2
	L = 1
	J = 3/2
	B = 0.6 # Tesla

	g = 1 + (J*(J+1)+S*(S+1)-L*(L+1))/(2*J*(J+1)) # Lande g-factor
	ylim = 1.5*J*((g*muB*B/h)/1e9) # calculated ylim based on Lande g-factor
	BreitRabi(I,L,S,J,ylim=ylim,Bmax=B)
	MomentumDecomp(I,L,S,J,B,ylim=ylim)

	plt.show()
