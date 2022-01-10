import os
import time
import shutil
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

from Solver import solver
# Run -defaults write org.python.python ApplePersistenceIgnoreState NO- before executing -python3 <filename>- if running from mac terminal





def phi_distribution_plot(velocity, gridpoints, sweep=False):
  
  # Solve for phi
  start_time = time.time()
  solve_CDS = solver('CDS',fluid_velocity=velocity,gridpoints=gridpoints)
  end_time = time.time()
  CDS_time = end_time-start_time

  start_time = time.time()
  solve_UDS = solver('UDS',fluid_velocity=velocity,gridpoints=gridpoints)
  end_time = time.time()
  UDS_time = end_time-start_time

  start_time = time.time()
  solve_PLDS = solver('PLDS',fluid_velocity=velocity,gridpoints=gridpoints)
  end_time = time.time()
  PLDS_time = end_time-start_time

  status_str = f"""Domain properties are set to: 
  Gamma_phi = {solve_CDS.get_Gamma()}
  Fluid Density = {solve_CDS.get_rho()} kg/m^3
  Fluid Velocity = {solve_CDS.get_u()} m/s
  Number of Gridpoints = {solve_CDS.get_gridpoint_n()}
  |Local Pe| = {abs(solve_CDS.get_local_Pe())}
  
  It took {CDS_time}s to solve the equation w/ CDS
  It took {UDS_time}s to solve the equation w/ UDS
  It took {PLDS_time}s to solve the equation w/ PLDS
  CDS % Error = {round(solve_CDS.get_error(),2)}
  UDS % Error = {round(solve_UDS.get_error(),2)}
  PLDS % Error = {round(solve_PLDS.get_error(),2)}
  
  Saving plot of phi distribution along the domain to all_plots/phi_distribution_plots/Individual_test_plots/u{int(velocity)}_gridpoints{gridpoints}...\n\n
  """


  # Plot phi
  fig = plt.figure()
  plt.plot(solve_CDS.get_anx(),solve_CDS.get_analytical_phi(),label='Analytical Solution')
  plt.plot(solve_CDS.get_x(),solve_CDS.get_numerical_phi(),label=f'Numerical Solution (CDS)',linewidth=0.75)
  plt.plot(solve_UDS.get_x(),solve_UDS.get_numerical_phi(),label=f'Numerical Solution (UDS)',linewidth=0.75)
  plt.plot(solve_PLDS.get_x(),solve_PLDS.get_numerical_phi(),label=f'Numerical Solution (PLDS)',linewidth=0.75)
      

  plt.title(fr'Distribution of $\phi$ along the domain (u={int(velocity)}m/s, gridpoints={gridpoints})')
  plt.annotate(fr'|Local Pe| = {abs(round(velocity/gridpoints,1))}',xy=(0.5,50))
  plt.xlabel('x (m)')
  plt.ylabel(r'$\phi$')
  plt.legend()

  if sweep: plt.savefig(f'all_plots/phi_distribution_plots/sweep/{int(velocity)}m_per_s/u{int(velocity)}_gp{gridpoints}',dpi=800)
  else: plt.savefig(f'Individual_test_plots/u{int(velocity)}_gp{gridpoints}',dpi=800)
  return status_str




## Plot error (all schemes) for different u values
def scheme_errors():
  
  for velocity in np.linspace(-50,50,11):
    fig = plt.figure()
    for scheme in ['CDS','UDS','PLDS']:
      error_list = []
      gp_list = np.linspace(1,50,50-1)
      for gridpoints in gp_list:
        solve = solver(scheme,fluid_velocity=velocity,gridpoints=int(gridpoints))
        error_list.append(solve.get_error())
      
      plt.plot(gp_list,error_list,label=f'scheme = {scheme}')
      plt.title(f'Error of numerical Solution for u={velocity}m/s')
      plt.ylabel('% Error ')
      plt.xlabel('Number of Gridpoints (1/dx)')

    plt.legend()
    # plt.ylim(-20,20)
    plt.savefig(f'all_plots/error_plots/scheme_errors/error_(u{int(velocity)})',dpi=800)





## Plot error (range of u values) for different schemes
def scheme_errors_func_u():

  for scheme in ['CDS','UDS','PLDS']:
    fig = plt.figure()
    max_gp = 200
    for u in np.linspace(-10,10,11):
      error_list = []
      gp_list = np.linspace(1,max_gp,max_gp-1)
      for gridpoints in gp_list:
        solve = solver(scheme,gridpoints=int(gridpoints),fluid_velocity=u)
        error_list.append(solve.get_error())


      plt.plot(gp_list,error_list,linewidth=0.5,label=f'u = {u}m/s')
      plt.title(f'Error of numerical {scheme} Solution w.r.t the Analytical Solution')
      plt.ylabel('% Error ')
      plt.xlabel('Number of Gridpoints (1/dx)')

    plt.legend(loc='lower right',prop={'size': 6})
    plt.savefig(f'all_plots/error_plots/scheme_errors_func_u/{scheme}_error_f(u)',dpi=800)

 


if __name__ == '__main__':

  #Generate Plot Sweep if folder doesnt exist
  if not os.path.exists('Individual_test_plots'):
    os.makedirs('Individual_test_plots')
  while True:
    print("Do you want to generate a range of error plots and phi distribution plots for different parameter combinations?\n\t yes \n\t no")
    input_a = input()
    if input_a.lower() == 'yes':
      if not os.path.exists('all_plots'):
        os.makedirs('all_plots/error_plots/scheme_errors')
        os.makedirs('all_plots/error_plots/scheme_errors_func_u')
        
        for velocity in [-50,-20,-5,0,5,20,50]:
          os.makedirs(f'all_plots/phi_distribution_plots/{int(velocity)}m_per_s')
      
      print('Currently generating and saving in the all_plots directory. This will take around 2 mins ...')

      #Swept Distribution profile
      for count,velocity in enumerate([-50,-20,-5,0,5,20,50]):
        print(f"\n{round(100/8*(count))}% completed\n")
        for gridpoints in [5,10,20,50,500]:
          phi_distribution_plot(velocity,gridpoints,sweep=True)
      
      print("\n90% completed\n")
      #Error Plots
      scheme_errors_func_u()
      scheme_errors()
      print('Done!\n')

    elif input_a.lower() == 'no': break
    else:
      print('Please answer: yes or no')
      continue
    print('You can find plots of a range of parameters in the all_plots directory.')
    break

  print('This UI is for testing distribution of phi with specific (Velocity, Gridpoint Number) pairs.\n')

  #Ask user what fluid velocity and number of gridpoints they would like to try
  attribute_prompt = {'constant fluid velocity':('m/s',' or float'),
                      'number of gridpoints':('','')}
  for attribute, dec in attribute_prompt.items():
    while True:
      print(f"What {attribute} would you like to test? Please enter an integer{dec[1]}.")
      num_attribute = input()
      try: num_attribute = int(num_attribute)
      except ValueError:
        try: num_attribute = float(num_attribute)
        except ValueError: continue
      if attribute == 'constant fluid velocity': 
        velocity = num_attribute
        print(f'The {attribute} in the domain has been set to {velocity} {dec[0]}')
      if attribute == 'number of gridpoints': 
        gridpoints = int(num_attribute)
        print(f'The {attribute} in the domain has been set to {gridpoints} {dec[0]}')
      break
  

  print('Solving equation in the Domain ...')
  status_str = phi_distribution_plot(velocity,gridpoints)
  print(status_str)
    

