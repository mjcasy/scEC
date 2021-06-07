import numpy as np
from scipy import optimize

def ident2mu(ident):
  
  N = ident.size
  C = np.unique(ident).size
  mu = np.zeros((N,C), 'int')
  
  for i in range(C):
    mu[ident == i, i] = 1
    
  return mu
  

def pop(freq):
  
  freq = np.array(freq)
  
  N = freq.shape[1]
  
  popG = np.sum(freq * np.log(N * freq), 1)
  
  return popG
  
      
def intertype(freq, ident):
  
  ident = np.array(ident).astype(int)
  freq = np.array(freq)
  
  mu = ident2mu(ident)
  
  y = freq @ mu
  N = ident.size
  
  Nk = np.sum(mu, 0)
  logNk = np.log(Nk)
  
  intertypeG = np.sum(y * (np.log(N * y) - logNk), 1)
  
  return intertypeG


def funcCost(wVec, freq, tfreq):
  N = freq.shape[1]
  C = int((wVec.size / N))
  
  wVecar = np.array(wVec)
  W = wVecar.reshape(N, C)
  M = np.sum(np.exp(W), 1)
  mu = np.exp(W) / M[:,None]
  
  y = freq @ mu
  
  Nj = np.sum(mu, 0)
  logNj = np.log(Nj)
  
  Is = -1*np.sum(y * (np.log(N * y) - logNj))
  
  return Is
  
  
def gradCost(wVec, freq, tfreq):
  N = freq.shape[1]
  C = int(wVec.size / N)
  
  wVecar = np.array(wVec)
  W = wVecar.reshape(N, C)
  M = np.sum(np.exp(W), 1)
  mu = np.exp(W) / M[:,None]
  
  y = freq @ mu
  
  Nj = np.sum(mu, 0)
  yNj = y / Nj
  
  term1 = np.log(yNj) + np.log(N) + 1
  dIsdu = tfreq @ term1
  term2 = -1*np.sum(yNj, 0)
  dIsdu = dIsdu + term2
  dIsdw = dIsdu * mu * (1 - mu)
  
  iterC = np.arange(C)
  
  for j in iterC:
    ind = iterC[iterC!=j]
    term2 = np.transpose(-1 * np.transpose(mu[:,ind]) * mu[:,j])
    dIsdw[:,j] = dIsdw[:,j] + np.sum(dIsdu[:,ind] * term2, 1)
  
  dIsdwVec = -1 * dIsdw.reshape((dIsdw.size,))
  return dIsdwVec
  
  
def clusterCells(freq, numClusters):
  freq = np.array(freq)
  numClusters = int(numClusters)
  
  cN = freq.shape[1]
  tfreq = freq.T
  wVec = np.random.uniform(low = -0.5, high = 0.5, size = (numClusters*cN,))
      
  bounding = 3
  Bounds=optimize.Bounds(lb=-bounding, ub=bounding)
  
  Out = optimize.minimize(funcCost, 
                          x0 = wVec, 
                          args = (freq, tfreq), 
                          method = 'L-BFGS-B',
                          jac=gradCost,
                          bounds=Bounds)
                          
  W = Out.x.reshape(cN, numClusters)
  M = np.sum(np.exp(W), 1)
  mu = np.exp(W) / M[:,None]
  
  return mu


def multiStartClusterCells(freq, numClusters, multistart):
  
  freq = np.array(freq)
  numClusters = int(numClusters)
  multistart = int(multistart)
  
  maxScore = 0

  for i in range(multistart):
    tempMu = clusterCells(freq, numClusters)
    newIdent = tempMu.argmax(1)
    Score = intertype(freq, newIdent).sum()
    
    if Score > maxScore:
      maxmu = tempMu
      maxScore = Score
    
  return maxmu


def greedyFuncCost(muVec, freq, tfreq):
  N = freq.shape[1]
  C = int((muVec.size / N) + 1)
  
  mu = np.zeros(shape = (N, C))
  mu[:,range(C-1)] = muVec.reshape(N, (C-1))
  mu[:,(C-1)] = 1 - mu.sum(1)
  
  y = freq @ mu
  
  Nj = np.sum(mu, 0)

  logNj = np.zeros(shape = Nj.shape)
  logNj[Nj > 0] = np.log(Nj[Nj > 0])
  
  logNy = np.zeros(shape = y.shape)
  logNy[y > 0] = np.log(N * y[y > 0])
  
  Is = -1*np.sum(y * (logNy - logNj))
  
  return Is


def greedyGradCost(muVec, freq, tfreq):
  
  G = freq.shape[0]
  N = freq.shape[1]
  C = int((muVec.size / N) + 1)
  
  mu = np.zeros(shape = (N, C))
  mu[:,range(C-1)] = muVec.reshape(N, (C-1))
  mu[:,(C-1)] = 1 - mu.sum(1)
  
  y = freq @ mu
  Nj = np.sum(mu, 0)
  
  yNj = y / Nj
      
  logyNj = np.zeros(shape = yNj.shape)
  logyNj[y > 0] = np.log(yNj[y > 0])
  
  term1 = logyNj[:,range(C-1)] - logyNj[:,(C-1)].reshape((G,1))
  
  dIsdu = tfreq @ term1 - \
    np.sum(yNj[:,range(C-1)], 0) + \
    np.sum(yNj[:,(C-1)])
    
  dIsduVec = -1 * dIsdu.reshape((dIsdu.size,))
      
  return dIsduVec


def splitCells(freq):
  C = 2
  N = freq.shape[1]
  tfreq = freq.T
  
  muVec = np.random.uniform(low = 0.45, high = 0.55, size = (N,1))
      
  bnd = optimize.Bounds(lb=0.1, ub=0.9)
  
  Out = optimize.minimize(greedyFuncCost, 
                          x0 = muVec, 
                          args = (freq, tfreq), 
                          method = 'L-BFGS-B',
                          jac=greedyGradCost,
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
    Score = intertype(freq, newIdent).sum()
    
    if Score > maxScore:
      maxmu = tempMu
      maxScore = Score
    
  return maxmu
    
