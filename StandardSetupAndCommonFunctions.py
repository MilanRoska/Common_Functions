# -*- coding: utf-8 -*-
# %% Packages
import matplotlib.pyplot as plt
import pandas as pd
import os

# %% Functions

# %% Figure setup


def standard_plot_parameters(ax, dark_mode=False):
    plt.tight_layout()
    # title and label size
    plt.rcParams['axes.titlesize'] = 20
    plt.rcParams['axes.labelsize'] = 18
    plt.rcParams['xtick.labelsize'] = 16
    plt.rcParams['ytick.labelsize'] = 16
    plt.rcParams['legend.fontsize'] = 12
    # legend positioning
    plt.legend(loc='best', handlelength=1, scatterpoints=1)
    # size of ticks
    ax.tick_params(axis='both', direction='in', length=5, width=1)
    ax.tick_params(which='major', size=5)  # Major ticks
    ax.tick_params(which='minor', size=5)   # Minor ticks
    # Show ticks on all four sides of the plot
    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')
    # switch to dark mode if selected
    if dark_mode:
        plt.style.use('dark_background')
    else:
        plt.style.use('classic')
    return


def standard_colors():
    orange = (244/255, 153/255, 62/255, 100/100)
    orange_light = (252/255, 230/255, 207/255, 100/100)
    orange_dark = (145/255, 77/255, 8/255, 100/100)
    orange_very_dark = (73/255, 39/255, 4/255, 100/100)
    purple = (215/255, 61/255, 245/255, 100/100)
    purple_light = (252/255, 240/255, 254/255, 100/100)
    purple_dark = (122/255, 8/255, 145/255, 100/100)
    purple_very_dark = (61/255, 4/255, 73/255, 100/100)
    return orange, orange_light, orange_dark, orange_very_dark, purple, purple_light, purple_dark, purple_very_dark


# %% Timestamp conversion
def convert_ldap_timestamp(ldap_timestamp):
    return pd.to_datetime(ldap_timestamp / 1e7 - 11644473600, unit='s')

def convert_unix_timestamp(nanoseconds):
    import pandas as pd
    # Convert nanoseconds to datetime
    return pd.to_datetime(nanoseconds, unit='ns')
# %% Load and write HDF files


# Read the HDF file and load DataFrames back into a dictionary
def open_hdf_to_dict(path):
    # generates empty dict
    opened_dict = {}
    with pd.HDFStore(path) as store:
        # iterates through HDF file keys and loads cotnent into dict
        for key in store.keys():
            opened_dict[key[1:]] = store[key]
    return opened_dict


# save dict as HDF file
def save_dict_to_hdf(dict_to_save, path):
    with pd.HDFStore(path) as store:
        # iterate through dict and save to hdf
        for key, df in dict_to_save.items():
            store[key] = df


# load multiple hdf files from directory and combine into one dict
def combine_hdf_dicts(folder_path):
    """
    Combines dictionaries of DataFrames from multiple HDF5 files in a folder.

    Parameters:
    folder_path (str): The path to the folder containing HDF5 files.

    Returns:
    dict: A dictionary where each key corresponds to a DataFrame that is the concatenation
          of all DataFrames with that key from the different HDF5 files.
    """
    combined_dict = {}

    # Iterate through all files in the folder
    for filename in os.listdir(folder_path):
        if filename.endswith('.h5'):  # Only consider HDF5 files
            file_path = os.path.join(folder_path, filename)
            print(f"Processing file: {file_path}")
            # Load the dictionary of DataFrames from the current HDF5 file
            current_dict = open_hdf_to_dict(file_path)
            # Combine the current dictionary with the combined dictionary
            for key, df in current_dict.items():
                if key in combined_dict:
                    # Append the DataFrame if the key already exists
                    combined_dict[key] = pd.concat([combined_dict[key], df], ignore_index=False)
                else:
                    # Otherwise, just add the DataFrame to the combined dictionary
                    combined_dict[key] = df

    return combined_dict


# convert dict to df and save df as hdf
def save_dict_via_df_to_hdf(dict_to_merge, path):
    """
    Saves the dictionary as a merged master DataFrame to an HDF5 file.
    Compression will save the former keys in a column labeled 'key'.

    Parameters:
    dict_to_merge (dict): A dictionary where keys are identifiers and values are DataFrames.
    path (str): The path where the HDF5 file will be saved.
    """
    master_df = pd.DataFrame()  # Initialize an empty DataFrame

    for key, df in dict_to_merge.items():
        df_copy = df.copy()  # Create a copy of the DataFrame to avoid modifying the original
        df_copy['key'] = key  # Add the key as a new column

        # Ensure all numeric columns are converted to float64
        for col in df_copy.select_dtypes(include=['int', 'float']).columns:
            df_copy[col] = df_copy[col].astype('float64')

        df_copy.reset_index(inplace=True)  # Save the index as a column
        master_df = pd.concat([master_df, df_copy], ignore_index=True)  # Append the DataFrame to the master

    if master_df.empty:
        raise ValueError("The provided dictionary is empty. No data to save.")

    master_df.to_hdf(path, key='master_df', mode='w', format='table')


def open_hdf_to_dict_via_df(path):
    """
    Opens the HDF5 file containing the master DataFrame and converts it back to a dictionary of DataFrames.

    Parameters:
    path (str): The path to the HDF5 file.

    Returns:
    dict: A dictionary where keys are the original keys and values are DataFrames.
    """
    # Load the master DataFrame from the HDF5 file
    master_df = pd.read_hdf(path, key='master_df')

    if master_df.empty:
        raise ValueError("The HDF5 file contains no data.")

    # Convert the master DataFrame back into a dictionary
    dict_from_master = {}
    for key in master_df['key'].unique():
        df = master_df[master_df['key'] == key].drop(columns=['key'])
        if 'index' in df.columns:
            df.set_index('index', inplace=True)  # Restore the index if it was saved as a column
        dict_from_master[key] = df

    return dict_from_master
