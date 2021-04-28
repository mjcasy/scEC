import numpy as np

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
