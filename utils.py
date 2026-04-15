# Third-party python imports
import shutil                       # for file copying
import subprocess                   # run bash processes through python
import numpy as np                  # vectorized operations
from scipy.optimize import brentq   # solve for decay scale for end of accretion

def calculate_tau(dM, mdot, t_f, fraction=0.8, mdot_min=1e-9):
    '''
        Function to solve for decay scale tau.
        This means the exponential decay at the end of accretion is fully defined by dM and mdot
        
        ! #### Important #### !
            We solve for tau, with a specified dM and t_f.
            dM influences result significantly. The great dM the longer t_f we can speficy.
            For example, if dM is small (~0.001 Msun) and t_f is large (~ 1 billion years),
            resulting decay with saturate way before 1 billion years.
            In the end we only accrete fraction * dM in order to fully capture an exponential decay.
            The remaining (1 - fraction) * dM is captured at mdot_min. 
        ! #### ######### #### !

        Inputs: 
            dM      :   amount of mass in solar masses to accrete in the final decay
            mdot    :   base accretion rate
            t_f     :   time-span for decay of accretion rate
        
        Outputs:
            tau     :   decay scale
    '''
    # mdot_min = 1e-9 # minimum accretion rate to aim for in decay (can be changed, but probably shouldn't)
    
    # Control what fraction of dM will be accreted in the exponential only,
    # i.e., to ensure that we reach mdot_min and get a smooth decay.
    dM *= fraction

    # ---- Safety checks ---- #
    if (dM - mdot*t_f) > 0:
        # No solution for this case
        raise ValueError("t_f too small to yield solution, increase t_f")

    if (dM - mdot*t_f) == 0:
        # Limiting case, solution for tau -> infinity
        raise ValueError("t_f limiting case, increase t_f")
    
    # We should be able to find a solution
    print("t_f sufficiently big, continue")

    # Check which regime we are in (calculate alpha)
    delta_mdot = mdot - mdot_min
    num = (dM - mdot_min * t_f)
    den = delta_mdot * t_f
    alpha = num / den
    print('alpha:', alpha)

    if alpha < 0.05:
        # Asymptotic regime
        print('Asymptotic regime')
        tau = (dM) / delta_mdot
        print(tau)
        return tau
    
    # We can't ignore mdot_min * t_f term
    print('Below asymptotic regime')
    def equation(beta):
        """Full transcendal equation we obtain after integration mdot profile"""
        delta_mdot = mdot - mdot_min
        num = (dM - mdot_min * t_f)
        den = delta_mdot * t_f
        alpha = num / den
        # beta = t_f/tau
        return beta * alpha - 1 + np.exp(-beta)

    # Solve for beta
    beta_solution = brentq(equation, 1e-8, 1e8)

    # Calculate and return decay scale (tau)
    tau = t_f / beta_solution

    return tau

def inplace_change(filename, old_string, new_string):
    # Safely read the input filename using 'with'
    with open(filename) as f:
        s = f.read()
        if old_string not in s:
            print('"{old_string}" not found in {filename}.'.format(**locals()))
            return

    # Safely write the changed content, if found in the file
    with open(filename, 'w') as f:
        print('Changing "{old_string}" to "{new_string}" in {filename}'.format(**locals()))
        s = s.replace(old_string, new_string)
        f.write(s)

def insert_inlist(main_inlist, placeholder, include_file):
    # Read main inlist
    with open(main_inlist, "r") as f:
        content = f.read()

    if placeholder not in content:
        print(f'Placeholder "{placeholder}" not found in {main_inlist}.')
        return

    # Read the file to insert
    with open(include_file, "r") as f:
        insert_text = f.read().rstrip()  # remove trailing newlines

    # Replace placeholder with content of other file
    new_content = content.replace(placeholder, insert_text)

    # Write back
    with open(main_inlist, "w") as f:
        f.write(new_content)

    print(f"Inserted contents of {include_file} into {main_inlist}.")


def check_mass(dM, mdot, tau, t_f, mdot_min=1e-9):
    delta_mdot = mdot - mdot_min
    M_calc = delta_mdot * tau * (1 - np.exp(-t_f / tau)) + mdot_min * t_f
    print(f"Target dM: {dM}, Calculated mass: {M_calc}, Difference: {dM - M_calc}")
    return M_calc
