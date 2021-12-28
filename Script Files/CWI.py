from Solver import solver
import time
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

import os
dir='Plots'
if not os.path.exists(dir):
    os.makedirs(dir)
# Run -defaults write org.python.python ApplePersistenceIgnoreState NO- before executing -python3 <filename>- if running from mac terminal



def phi_distribution_plot(scheme, velocity):

  solve = solver(scheme, velocity)
  print(f"""Domain properties are set to: 
        Gamma_phi = {solve.get_Gamma()}
        Fluid Velocity = {solve.get_u()} m/s
        Fluid Density = {solve.get_rho()} kg/m^3
        Number of Gridpoints = {solve.get_gridpoint_n()}""")
  print('Solving equation in the Domain ...')

  # Solve for phi
  analytical_phi = solve.get_analytical_phi()
  start_time = time.time()
  numerical_phi = solve.get_numerical_phi()
  end_time = time.time()
  print(f"It took {end_time-start_time}s to solve the equation.")
  print(f'% Error = {round(solve.get_error(),2)}')
  print(f'Local Pe = {solve.get_local_Pe()}')
  print(f"Saving plot of phi distribution along the domain to {dir}/{scheme}...\n\n")

  # Plot phi
  fig = plt.figure()
  plt.plot([x/10/len(solve.x) for x in range(10*len(solve.x))],analytical_phi,label='Analytical Solution')
  plt.plot(solve.get_x(),numerical_phi,label=f'Numerical Solution ({scheme})')

  plt.title(r'Distribution of $\phi$ along the domain')
  plt.xlabel('x')
  plt.ylabel(r'$\phi$')
  plt.legend()
  # plt.show()
  plt.savefig(f'{dir}/{scheme}',dpi=800)


if __name__ == '__main__':

  scheme_dict = {'CDS':'Central Differencing Scheme',
               'UDS':'Upwind Differencing Scheme',
               'PLDS':'Power Law Differencing Scheme'}

  #Ask user what fluid velocity they would like to try
  while True:
    print("What fluid velocity would you like to test? Please enter a float or an integer.")
    velocity = input()
    try: velocity = int(velocity)
    except ValueError:
      try: velocity = float(velocity)
      except ValueError: continue
    break
  print(type(velocity))
  print(f'The constant fluid velocity in the domain has been set to {velocity} m/s')

  #Ask user what type of differencing scheme they would like to try
  while True:
    print("Please choose a discretization scheme. Type one of the following options into the command prompt: \n\tCDS \n\tUDS \n\tPLDS")
    scheme_text = input()
    if scheme_text in ['CDS','UDS','PLDS']: break
    print("Please ensure you have entered the acronym in all caps with no spaces!")
  print(f"You have chosen to discretize the convection component of the phi convection-diffusion equation with the {scheme_dict[scheme_text]}\n\n")

  #Plot the numerical solution against the analytical solution
  phi_distribution_plot(scheme_text,velocity)
