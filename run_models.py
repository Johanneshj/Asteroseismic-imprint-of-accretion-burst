# Third-party python imports
import shutil                       # for file copying
import subprocess                   # run bash processes through python
import numpy as np                  # vectorized operations
from scipy.optimize import brentq   # solve for decay scale for end of accretion

# Copy-paste functions
from utils import calculate_tau, inplace_change, insert_inlist, check_mass

# Parameters we want to set for the run
# Initial model parameters
initial_model       = f'initial_model.mod'
seed_mass           = '0.01'
mesh_delta_coeff    = '0.8' # MESA's mesh_delta_coeff
initial_z           = '0.020'
initial_h2          = '20' # ppm
initial_he3         = '85' # ppm
seed_radius         = '1.5' # Relax to this radius

burst_mag           = '1d-3' # magnitude of burst in solar masses per year
time_delta_coeff    = '0.8' # MESA's time_delta_coeff

pre_burst_model     = f'pre_burst_model.mod' # name of pre_burst_model saved at end of run
post_burst_model    = f'post_burst_model.mod' # name of post_burst_model saved at end of run 
ZAMS_mass           = '4.0d0' # Final mass at the ZAMS in solar masses
dM                  = '0.05' # mass to accrete in decay of accretion rate (solar masses)
mdot                = '1d-6' # mass accretion rate outside of burst (must be given here as 1d-6 fx)
decay_time          = '1000000' # Time scale of exponential decay
log_dir             = f'LOGS'
png_dir             = f'png_pre_burst'
png_prefix          = f'pre_burst_'
photos_dir          = f'photos'        

# Should we do a burst?
do_burst = True
if do_burst:
    number_of_bursts = '1'
    int_ctrl_10 = '0'
else:
    number_of_bursts = '0'
    int_ctrl_10 = '2'

# Should we save models?
save_models = True
if save_models:
    save_models = '.true.'
else:
    save_models = '.false.'

# Copy clean template inlists into folders
shutil.copy('0_inlists_templates/inlist_common_base', '5_inlists_common/inlist_common')
shutil.copy('0_inlists_templates/inlist_common_single_burst_base', '5_inlists_common/inlist_common_single_burst')
shutil.copy('0_inlists_templates/inlist_initial_model_base', '1_inlists_initial_model/inlist_initial_model')
shutil.copy('0_inlists_templates/inlist_pre_burst_base', '2_inlists_pre_burst/inlist_pre_burst')
shutil.copy('0_inlists_templates/inlist_burst_base', '3_inlists_burst/inlist_burst')
shutil.copy('0_inlists_templates/inlist_post_burst_base', '4_inlists_post_burst/inlist_post_burst')
shutil.copy('0_inlists_templates/inlist_pgstar_pre_burst_base', '7_inlists_pgstar/inlist_pgstar_pre_burst')
shutil.copy('rn_template', 'rn')

# -------- CHANGE RN -------- # 
inplace_change('rn',
               '<<init_model>>',
               initial_model)
inplace_change('rn',
               '<<pre_burst_model>>',
               pre_burst_model)
inplace_change('rn',
               '<<post_burst_model>>',
               post_burst_model)

# -------- CHANGE INLISTS -------- #

# INITIAL MODEL INLIST
inlists_to_paste = ['inlist_common_mesh', 'inlist_common_mix']
strings_to_replace = ['<<mesh_controls>>', '<<mixing_controls>>']
for inlist_to_paste, string_to_replace in zip(inlists_to_paste, strings_to_replace):
    insert_inlist(
        '1_inlists_initial_model/inlist_initial_model',
        f'{string_to_replace}',
        f'5_inlists_common/{inlist_to_paste}'
    )

inplace_change(f'1_inlists_initial_model/inlist_initial_model', 
               '<<initial_model_name>>', 
               initial_model)

inplace_change(f'1_inlists_initial_model/inlist_initial_model', 
               '<<mesh_coeff>>', 
               mesh_delta_coeff)

inplace_change(f'1_inlists_initial_model/inlist_initial_model', 
               '<<seed_mass>>', 
               seed_mass)

inplace_change(f'1_inlists_initial_model/inlist_initial_model', 
               '<<mass_to_relax_from>>', 
               f'{float(seed_mass)+0.02}')

inplace_change(f'1_inlists_initial_model/inlist_initial_model', 
               '<<initial_z>>', 
               initial_z)

inplace_change(f'1_inlists_initial_model/inlist_initial_model', 
               '<<seed_radius>>', 
               seed_radius)

inplace_change(f'1_inlists_initial_model/inlist_initial_model', 
               '<<zbase>>', 
               initial_z)

inplace_change(f'1_inlists_initial_model/inlist_initial_model', 
               '<<initial_model_name>>', 
               initial_model)

inplace_change(f'1_inlists_initial_model/inlist_initial_model', 
               '<<h2>>', 
               initial_h2)

inplace_change(f'1_inlists_initial_model/inlist_initial_model', 
               '<<he3>>', 
               initial_he3)


               
# PRE BURS INLIST
inplace_change(f'2_inlists_pre_burst/inlist_pre_burst', 
               '<<save_models>>', 
               save_models)

inplace_change(f'2_inlists_pre_burst/inlist_pre_burst',
               '<<initial_model>>',
               initial_model)

inplace_change(f'2_inlists_pre_burst/inlist_pre_burst', 
               '<<pre_burst_model_name>>', 
               pre_burst_model)

inplace_change(f'2_inlists_pre_burst/inlist_pre_burst', 
               '<<int_ctrl_10>>', 
               int_ctrl_10)

# BURST INLIST
inplace_change(f'3_inlists_burst/inlist_burst', 
               '<<save_models>>', 
               save_models)

inplace_change(f'3_inlists_burst/inlist_burst', 
               '<<pre_burst_model_name>>', 
               pre_burst_model)

inplace_change(f'3_inlists_burst/inlist_burst', 
               '<<post_burst_model_name>>', 
               post_burst_model)

# POST BURST INLIST
inplace_change(f'4_inlists_post_burst/inlist_post_burst', 
               '<<post_burst_model_name>>', 
               post_burst_model)

# COMMON INLIST
inplace_change(f'5_inlists_common/inlist_common_single_burst', 
               '<<number_of_bursts>>', 
               number_of_bursts)

inplace_change(f'5_inlists_common/inlist_common_single_burst', 
               '<<burst_mag>>', 
               burst_mag)

inplace_change(f'5_inlists_common/inlist_common', 
               '<<ZAMS_mass>>', 
               ZAMS_mass)

inplace_change(f'5_inlists_common/inlist_common',
               '<<mdot>>',
               mdot)

inplace_change(f'5_inlists_common/inlist_common',
               '<<dM>>',
               dM)

inplace_change(f'5_inlists_common/inlist_common',
               '<<decay_time>>',
               decay_time)

inplace_change(f'5_inlists_common/inlist_common',
               '<<h2>>',
               initial_h2)

inplace_change(f'5_inlists_common/inlist_common',
               '<<he3>>',
               initial_he3)

inplace_change(f'5_inlists_common/inlist_common',
               '<<initial_z>>',
               initial_z)

inplace_change(f'5_inlists_common/inlist_common',
               '<<zbase>>',
               initial_z)

# Set mesh
inplace_change(f'5_inlists_common/inlist_common', 
               '<<mesh_delta_coeff>>', 
               mesh_delta_coeff)

# Set time delta coeff 
inplace_change(f'5_inlists_common/inlist_common', 
               '<<time_delta_coeff>>', 
               time_delta_coeff)
               
# Set decay scale for end of accretion
# We calculate this by assuming that some amount of user-specified mass will be accreted over 5e5 years
tau = calculate_tau(float(dM), float(mdot.replace('d', 'e')), float(decay_time))
inplace_change(f'5_inlists_common/inlist_common',
               '<<tau>>',
               str(tau))

# LOGS dir
inplace_change(f'5_inlists_common/inlist_common',
               '<<log_directory>>',
               log_dir)

# Photos dir
inplace_change(f'5_inlists_common/inlist_common',
               '<<photo_directory>>',
               photos_dir)

# INLIST PGSTAR
inplace_change(f'7_inlists_pgstar/inlist_pgstar_pre_burst', 
               '<<png_folder>>', 
               png_dir)
inplace_change(f'7_inlists_pgstar/inlist_pgstar_pre_burst', 
               '<<png_prefix>>', 
               png_prefix)
# -------- --------- #

# Copy inlist files into one big inlist common
inlist_files = ['inlist_common_mesh', 'inlist_common_mix',
                'inlist_common_single_burst', 'inlist_common_timestep'] 

replace_strings = ['<<mesh_controls>>', '<<mixing_controls>>',
                   '<<burst_parameters>>', '<<timestep_controls>>']

for inlist_file, replace_strings in zip(inlist_files, replace_strings):
    insert_inlist(
        '5_inlists_common/inlist_common',
        f'{replace_strings}',
        f'5_inlists_common/{inlist_file}'
    )

# Run model
subprocess.run(f"./clean; ./mk; ./rn", shell=True, check=True)