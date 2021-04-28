import numpy as np
from scipy import optimize
import CostFunction as cf
import Heterogeneity as het

def splitCells(freq):
  C = 2
  N = freq.shape[1]
  tfreq = freq.T
  
  muVec = np.random.uniform(low = 0.45, high = 0.55, size = (N,1))
      
  bnd = optimize.Bounds(lb=0, ub=1)
  
  Out = optimize.minimize(cf.funcCost, 
                          x0 = muVec, 
                          args = (freq, tfreq), 
                          method = 'L-BFGS-B',
                          jac=cf.gradCost,
                          bounds=bnd)
                          
  mu = np.zeros(shape = (N, C))
  mu[:,range(C-1)] = Out.x.reshape(N, (C-1))
  mu[:,(C-1)] = 1 - mu.sum(1)
  
  return mu


def multiStartSplitCell(freq, multistart):
  
  freq = np.array(freq)
  multistart = int(multistart)
  
  maxScore = 0

  for i in range(multistart):
    tempMu = splitCells(freq)
    newIdent = tempMu.argmax(1)
    Score = het.intertype(freq, newIdent).sum()
    
    if Score > maxScore:
      maxmu = tempMu
      maxScore = Score
    
  return maxmu
    
