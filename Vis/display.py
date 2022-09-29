#!/usr/bin/env python

# -----------------------------------------------------------------------------
#  basic.py
#
#  Example how to load in a data file and use the data in a simple plot.
# -----------------------------------------------------------------------------

import argparse
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import uproot
import plotly.express as px
matplotlib.rcParams.update({'font.size': 8})

# parse arguments from command
parser = argparse.ArgumentParser(description="pixel event display")
parser.add_argument("file", type=str, default=None,
                    help="path to ROOT file")

args = parser.parse_args()
file_path = str(args.file)

with uproot.open(file_path) as f:

    #--------------------------------------------------------------------------
    # get event tree from ROOT file
    #--------------------------------------------------------------------------

    tree = f['event_tree']

    # list of branches to access
    branches = [

        # event number
        'event',

        # generator particle information [Q_PIX_GEANT4]
        'generator_initial_number_particles',
        'generator_initial_particle_pdg_code',
        'generator_initial_particle_energy',
        'generator_initial_particle_mass',
        'generator_final_number_particles',
        'generator_final_particle_pdg_code',
        'generator_final_particle_energy',
        'generator_final_particle_mass',

        # MC particle information [Q_PIX_GEANT4]
        'particle_track_id', 'particle_pdg_code',
        'particle_mass', 'particle_initial_energy',

        # MC hit information [Q_PIX_GEANT4]
        'hit_energy_deposit', 'hit_track_id', 'hit_process_key',
        'hit_start_x', 'hit_start_y', 'hit_start_z', 'hit_start_t',
        'hit_end_x', 'hit_end_y', 'hit_end_z', 'hit_end_t',


    ]

    #--------------------------------------------------------------------------
    # iterate through the event tree
    #--------------------------------------------------------------------------

    for arrays in tree.iterate(filter_name=branches):

        # get event number array
        event_array = arrays['event']

        # get generator particle arrays
        gen_ip_number_array   = arrays['generator_initial_number_particles']
        gen_ip_pdg_code_array = arrays['generator_initial_particle_pdg_code']
        gen_ip_energy_array   = arrays['generator_initial_particle_energy']
        gen_ip_mass_array     = arrays['generator_initial_particle_mass']
        gen_fp_number_array   = arrays['generator_final_number_particles']
        gen_fp_pdg_code_array = arrays['generator_final_particle_pdg_code']
        gen_fp_energy_array   = arrays['generator_final_particle_energy']
        gen_fp_mass_array     = arrays['generator_final_particle_mass']

        # get MC hit arrays
        hit_end_x_array = arrays['hit_end_x']
        hit_end_y_array = arrays['hit_end_y']
        hit_end_z_array = arrays['hit_end_z']
        hit_end_t_array = arrays['hit_end_t']
        hit_energy_deposit_array = arrays['hit_energy_deposit']

        # get number of events
        number_events = len(event_array)

        #----------------------------------------------------------------------
        # iterate over events
        #----------------------------------------------------------------------

        idx = 0
        while idx < 10:

            #------------------------------------------------------------------
            # fetch event information
            #------------------------------------------------------------------

            # get event number
            event = event_array[idx]
            print(idx)

            # get MC hit information for event
            hit_end_x = hit_end_x_array[idx]
            hit_end_y = hit_end_y_array[idx]
            hit_end_z = hit_end_z_array[idx]
            hit_end_t = hit_end_t_array[idx]
            hit_energy_deposit = hit_energy_deposit_array[idx]

            # get number of hits in event
            number_hits = len(hit_end_x)
            fig = plt.figure(figsize=(10, 10))
            ax = plt.axes(projection="3d")
            p=ax.scatter3D(hit_end_x, hit_end_y, hit_end_z, c=hit_end_t, cmap='viridis')
            fig.colorbar(p)
            # ax.set_xlim(0,600)
            # ax.set_ylim(0,600)
            # ax.set_zlim(0,600)

            plt.show()
            #df = pd.DataFrame({'x': hit_end_x, 'y': hit_end_y, 'z': hit_end_z, 't': hit_end_t, 'e': hit_energy_deposit})

            #fig = px.scatter_3d(df, x='x', y='y', z='z',color='t')

            #fig.show()


            idx +=1
            #------------------------------------------------------------------
        break
