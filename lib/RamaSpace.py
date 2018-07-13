import numpy as np
from fractions import gcd
#from math import pi as pi

class RamaSP:
	def __init__(self, P_series = None, flag_normalized = False):
		self.S = dict()
		self.euler = dict()
		self.P2index = dict()
		if P_series is None:
			self.P_series = []
			return 
		self.P_series = P_series
		for p in P_series:
			self.put_in_period(p,flag_normalized)
	def put_in_period(self, p, flag_normalized):
		#sequential method
		epsnum = np.array([])
		eulernum = 0
		for k in range(1,int(p+1)):
		    if gcd(k,p) == 1:
		        eulernum = eulernum +1
		        epsnum = np.hstack( (epsnum,1j*2*np.pi*k/p) )
		total_element = np.matrix(np.arange(0,p)).T*np.matrix(epsnum)
		noshift = np.round(np.real(np.sum(np.exp(total_element), axis=1))).T;
		C = noshift
		for i in range(1,eulernum):
			C = np.vstack( (C,np.roll(noshift,i,axis=1)) )
		if flag_normalized:
			C = C/p
		self.S[int(p)] = C
		self.euler[int(p)] = eulernum
	def addSpace(self, P_series, flag_normalized = False):
		self.P_series = np.unique(np.hstack((self.P_series,P_series)))
		for p in P_series:
			if p not in self.S:
				self.put_in_period(p,flag_normalized)

	def GetSpaceByP(self, p):
		return self.S[p]
	def geteulerByP(self, p):
		return self.euler[int(p)]
	#Return a stantard row vector dictionary to 
	def DictTransform(self, Dicdim = None, ConvSpace = False):

		maxp = np.max(self.P_series)
		# euler = []
		# for p in self.P_series:
			# euler.append(self.euler[int(p)]) 
		if Dicdim is None or maxp > Dicdim:
			Dicdim = int(maxp*2)
		Dicdim = int(Dicdim)			
		if ConvSpace == False:
			D = np.zeros([ Dicdim,sum(self.euler.values())])
		else:
			D = np.zeros([Dicdim,len(self.P_series)])
		itr = 0
		for p in self.P_series:
			if ConvSpace == False:
				instan_size = self.euler[int(p)]
				padding = self.S[int(p)]
				padding = np.tile(padding,int(np.floor(Dicdim/p)))
				padding = np.hstack((padding,padding[:,0:int(Dicdim % p)]))
			else:
				instan_size = 1
				padding = (self.S[int(p)])[0,:]
				padding = np.tile(padding,int(np.floor(Dicdim/p)))
				padding = np.hstack((padding,padding[0,0:int(Dicdim % p)]))#np.pad(padding,((0,0),(0,int(Dicdim % p))),'constant')
			if ConvSpace == False:
				D[:,itr:itr+instan_size] = padding.T
				self.P2index[int(p)] = range(itr,itr+instan_size)
			else:
				D[:,itr] = padding
				self.P2index[int(p)] = itr
			itr += instan_size
		return D
#Unit test
if __name__ == '__main__':
	print("...Start to calculate the basises of RS subspace...")
	S = RamaSP(np.arange(1,11))
	#Test the generated p = 1-10 space which is [1] [1 -1] [2 -1 -1] ...
	print("...Output the basises...")
	for i in range(1,11):
		print (S.GetSpaceByP(i)[0])
	#Test the generated correspoding dictionary with shape (10,sum(euler(1:10))
	print("...Truncate the basises as an dictionary...(column vector)")
	D = S.DictTransform(ConvSpace = False)
	print (D)