# Python routines for reading MESA output files stolen from Earl Bellinger!
## https://github.com/earlbellinger/mesa-gyre-tutorial-2022

import os
import pandas as pd

class DataHandler:
    def __init__(self, data_dir, logs_dir='LOGS', identifier='_burst'):
        """
        Initialize the LogDataHandler with the directory containing log files.
        """
        self.logs_dir = f'{data_dir}/{logs_dir}'
        self.gyre_dir = f'{data_dir}/GYRE'
        self.identifier = identifier
        self.history_data = self.load_history_file()
        self.index_data = self.get_index()

    def load_history_file(self):
        """
        Load the history file from the logs directory.
        """
        return pd.read_table(
            os.path.join(self.logs_dir, f'history{self.identifier}.data'), 
            skiprows=5, sep='\s+'
        )

    def get_index(self):
        """
        Load the profiles index file from the logs directory.
        """
        return pd.read_table(
            os.path.join(self.logs_dir, f'profiles{self.identifier}.index'), 
            names=['model_number', 'priority', 'profile_number'],
            skiprows=1, sep='\s+'
        )

    def get_history(self, profile_number):
        """
        Get the history data for a specific profile number.
        """
        model_number = self.index_data[
            self.index_data.profile_number == profile_number
        ].model_number.values[0]
        return self.history_data[self.history_data.model_number == model_number]

    def load_profile(self, profile_number):
        """
        Load the profile data for a specific profile number.
        """
        return pd.read_table(
            os.path.join(self.logs_dir, f'profile{self.identifier}{profile_number}.data'), 
            skiprows=5, sep='\s+'
        )
    
    def load_nad_gyre_file(self, profile_number):
        return pd.read_table(
            os.path.join(self.gyre_dir, 
                            f'profile{profile_number}.data.GYRE.summary.nad.txt'), 
                            sep='\s+', skiprows=5
            )
    
    def load_ad_gyre_file(self, profile_number):
        return pd.read_table(
            os.path.join(self.gyre_dir, 
                            f'profile{profile_number}.data.GYRE.summary.ad.txt'), 
                            sep='\s+', skiprows=5
            )
