
import numpy as np

def funcCost(muVec, freq, tfreq):
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


def gradCost(muVec, freq, tfreq):
  
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
  
