import os
import shutil
import sys
import subprocess
import pandas as pd

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

suffix = 'burst'
burst_mags = ['1d-3']#['5d-6', '1d-5', '5d-5', '1d-4', '5d-4', '1d-3']
for burstmag in burst_mags:
    directory = '.' #f'data/{burstmag}_data_new_profile_saves'
    logs_dir = f'{directory}/LOGS/'

    if not os.path.exists(f'{directory}/GYRE'):
        os.mkdir(f'{directory}/GYRE')
    gyre_dir = f'{directory}/GYRE'

    shutil.copy(f'gyre_template.in', f'{directory}/gyre_{burstmag}.in')

    inplace_change(f'{directory}/gyre_{burstmag}.in', 
                '<<directory>>', 
                directory)
    inplace_change(f'{directory}/gyre_{burstmag}.in', 
                '<<gyre_dir>>', 
                gyre_dir)

    index_file = pd.read_table(
                    os.path.join(logs_dir, f'profiles_{suffix}.index'), 
                    names=['model_number', 'priority', 'profile_number'],
                    skiprows=1, sep='\s+'
                )


    for profile_number in index_file['profile_number']:
        shutil.copy(f'{directory}/gyre_{burstmag}.in', 
                    f'{directory}/gyre_{burstmag}_now.in')
        file_name = f'profile_{suffix}{profile_number}.data.GYRE'
        inplace_change(f'{directory}/gyre_{burstmag}_now.in', 
                '<<file_name>>', 
                file_name)
        subprocess.run(f"$GYRE_DIR/bin/gyre {directory}/gyre_{burstmag}_now.in", 
                    shell=True, check=True)
        subprocess.run(f'rm {directory}/gyre_{burstmag}_now.in',
                    shell=True, check=True)
        
