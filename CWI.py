import solver
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

import os
dir='Plots'
if not os.path.exists(dir):
    os.makedirs(dir)
# Run -defaults write org.python.python ApplePersistenceIgnoreState NO- before executing -python3 <filename>- if running from mac terminal



def comp_phi_plot(scheme):

  solve = solver(scheme)
  print(f"""Domain properties are defined as: 
          \t Gamma_phi = {solve.get_Gamma()}
          \t Fluid Velocity = {solve.get_u()} m/s
          \t Fluid Density = {solve.get_rho()} kg/m^3
          \t Number of Gridpoints = {solve.get_gridpoint_n()}""")
  print('Solving equation in the Domain ...')

  # Solve for phi
  analytical_phi = solve.get_analytical_phi()
  if scheme == "CDS": numerical_phi = solve.get_numerical_phi()
  elif scheme == "UDS": numerical_phi = solve.UDS()
  elif scheme == "PLDS": numerical_phi = solve.PLDS()
  else:
    print("Please enter a valid discrtization scheme acronym")
    return

  # Plot phi
  fig = plt.figure()
  plt.plot(solve.get_x(),analytical_phi,label='Analytical Solution')
  plt.plot(solve.get_x(),numerical_phi,label=f'Numerical Solution ({scheme})')

  plt.title(r'Distribution of $\phi$ along the domain')
  plt.xlabel('x')
  plt.ylabel(r'$\phi$')
  plt.legend()
  # plt.show()
  plt.savefig(f'{dir}/{scheme}',dpi=800)
  print(f"Plot of Distribution of phi along the domain was saved to {dir}/{scheme}\n\n")


if __name__ == '__main__':

  scheme_dict = {'CDS':'Central Differencing Scheme',
                'UDS':'Upwind Differencing Scheme',
                'PLDS':'Power Law Differencing Scheme'}

  end_condition = False
  while not end_condition:
    #Ask user what type of differencing scheme they would like to try
    while True:
      print("Please choose a discretization scheme. Type one of the following options into the command prompt: \n\tCDS \n\tUDS \n\tPLDS")
      scheme_text = input()
      if scheme_text in ['CDS','UDS','PLDS']: break
      print("Please ensure you have entered the acronym in all caps with no spaces!")
    print(f"You have chosen to discretize the convection component of the phi convection-diffusion equation with the {scheme_dict[scheme_text]}\n\n")

    #Plot the numerical solution against the analytical solution
    comp_phi_plot(scheme_text)

    while True:
      print("Do you want to terminate the program? \n\tyes \n\tno")
      break_cmd = input()
      if break_cmd == 'yes': end_condition = True
      if break_cmd in ['yes','no']: break
      print("Please type your response as presented above.")