import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
class solver():

  def __init__(self,
              scheme_acronym, 
              fluid_velocity = 1,
              fluid_density = 0.5,
              Gamma = 0.5,
              domain_L = 1,
              inlet_bc = 100,
              outlet_bc = 20,
              gridpoints = 10
              ):
    
    # Global Attributes
    self.scheme = scheme_acronym
    
    self.u = fluid_velocity
    self.rho = fluid_density
    self.Gamma = Gamma
    self.L = domain_L
    self.dx = self.L/gridpoints
    self.global_Pe = self.rho * self.u * self.L / self.Gamma
    self.local_Pe = self.rho * self.u * self.dx / self.Gamma
    
    self.x = np.array([i*self.dx for i in range(gridpoints)])
    self.phi_0 = inlet_bc
    self.phi_L = outlet_bc
    
    self.analytical_phi = np.zeros(50*len(self.x))
    self.phi = np.zeros(gridpoints)
    self.P = np.zeros(gridpoints)
    self.Q = np.zeros(gridpoints)
    self.error = 0

    #Solve for phi and calc error
    if self.scheme == "CDS": self.CDS()
    elif self.scheme == "UDS": self.UDS()
    elif self.scheme == "PLDS": self.PLDS()
    self.analytical_soln()
    self.numerical_error()

  # Getter methods
  def get_local_Pe(self): return self.local_Pe
  def get_x(self): return self.x
  def get_P(self): return self.P
  def get_Q(self): return self.Q
  def get_numerical_phi(self): return self.phi
  def get_analytical_phi(self): return self.analytical_phi
  def get_rho(self): return self.rho
  def get_u(self): return self.u
  def get_Gamma(self): return self.Gamma
  def get_gridpoint_n(self): return int(self.L/self.dx)
  def get_error(self): return self.error


  # Discretization Schemes
  def CDS(self):
    D = self.Gamma / self.dx
    F = self.rho * self.u
    #Coefficients
    a = 2*D
    b = D - 0.5 * F
    c = D + 0.5 * F

    self.TDMA(a,b,c,0)
  
  def UDS(self):
    D = self.Gamma / self.dx
    F = self.rho * self.u
    #Coefficients
    a = 2*D + max(F,0) + max(-F,0)
    b = D + max(-F,0)
    c = D + max(F,0)
  
    self.TDMA(a,b,c,0)

  def PLDS(self):
    D = self.Gamma / self.dx
    F = self.rho * self.u
    #Coefficients
    b = D*max((1-0.1*self.local_Pe)**5,0) + max(-F,0)
    c = D*max((1-0.1*self.local_Pe)**5,0) + max(F,0)
    a = b + c
  
    self.TDMA(a,b,c,0)




  ## Use Linear Algebraic Equation Solution Algorithm TDMA
  def TDMA(self,a,b,c,d):
    # Set BCs
    self.P[0], self.P[-1] = 0,0
    self.Q[0], self.Q[-1] = self.phi_0, self.phi_L

    # Forward Pass to calculate P and Q
    for i in [x+1 for x in range(len(self.x)-2)]:
      self.P[i] = b / (a - c*self.P[i-1])
      self.Q[i] = c*self.Q[i-1] / (a - c*self.P[i-1])

    # Backward Pass to calculate phi
    self.phi[-1] = self.Q[-1]
    for i in [y for y in range(len(self.x)-1)][::-1]:
      self.phi[i] = self.P[i] * self.phi[i+1] + self.Q[i]
  



  def analytical_soln(self):
    x_domain = np.array([x/50/len(self.x) for x in range(50*len(self.x))])
    self.analytical_phi = self.phi_0 + (np.exp(self.global_Pe * x_domain/self.L)-1)/(np.exp(self.global_Pe)-1) * (self.phi_L-self.phi_0)

  def numerical_error(self):
    for i in range(len(self.phi)):
      self.error += 100*self.dx/self.L * (self.phi[i]-self.analytical_phi[10*i])/self.analytical_phi[10*i]
    

  