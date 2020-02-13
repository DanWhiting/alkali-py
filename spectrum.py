'''
Author - Daniel J. Whiting 
Date modified - 10/08/2017

Contains functions for calculating weak probe transmission spectra.
'''

from __future__ import division
from uncoupledbasis import *

def StickSpectrum(transition='D2',B=0,pol_frac=0.5):
	''' A function for calculating stick spectra. '''
	L1 = 0; L2 = 1; J1 = 1/2
	Is = [3/2,5/2] ### Isotopes
	iso_fracs = [0.2783,0.7217] ### Relative isotopic abundances
	trans_strengths = []
	trans_freqs = []
	for Ii in range(len(Is)):
		I = Is[Ii]
		iso_frac = iso_fracs[Ii]
		if transition == 'D1':
			J2 = 1/2
			eibottom = 0
			eitop = int(2 * 2*I+1 + 1e-10)
			if Ii == 0:
				IsotopeShift = -56.077
			elif Ii == 1:
				IsotopeShift = 21.624
		elif transition == 'D2':
			J2 = 3/2
			eibottom = int(2 * (2*I+1) + 1e-10)
			eitop = int(3 * 2 * (2*I+1) -1 + 1e-10)
			if Ii == 0:
				IsotopeShift = -56.361
			elif Ii == 1:
				IsotopeShift = 21.734
		g_evals, g_evecs = EvalsEvecs(I,0,1/2,1/2,B,IsotopeShift)
		e_evals, e_evecs = EvalsEvecs(I,1,1/2,J2,B)
		lg = g_evecs.shape[0]
		le = e_evecs.shape[0]
		lm = (le-lg)/2
		lp = (le+lg)/2
		for q in [-1,+1]:
			for gi in range(len(g_evals)):
				for ei in range(eibottom,eitop+1):
					cleb = dot(g_evecs[:,gi],e_evecs[int(lm*(1-q)):int(lp-lm*q),ei])
					cleb2 = cleb*cleb
					if cleb2*pol_frac*iso_frac > 0.0005:
						strength = cleb2 * (2*J1+1)/(2*J2+1)  # Comes from ground state manifold reduction
						trans_strengths.append(abs(strength)*pol_frac*iso_frac)
						trans_freqs.append((e_evals[ei]-g_evals[gi])/1e9)
	cs = ['#AA2B4A','#006388','#7E317B'] # durham colours red,blue,purple
	fig, ax = plt.subplots(1,facecolor='white')
	ax.vlines(trans_freqs,0,trans_strengths,lw=4,color=cs[0])
	
	
if __name__=="__main__":
	B = 0 # Tesla
	
	StickSpectrum(B=B)

	plt.show()
