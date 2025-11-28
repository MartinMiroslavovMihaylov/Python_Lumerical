
# Import some Python Libs 
import importlib
import numpy as np
import os
import os.path
import sys

# Import the Constructor 
from Constructor import Constructor, Charge
from Constructor import loadingBar
from Constructor import Logfile
from Constructor import Help

# =============================================================================
# Open File
# =============================================================================


Lumerical_Python_Path = "C:/Program Files/Lumerical/v242/api/python/lumapi.py"

Added_Material_Path = "...../Materials.mdf"


# ============================================================================
# Start Constructor 
# =============================================================================



# # If you need Help with hot to import the Constructor run the line below 
# Help()

# Call object Constructor
obj = Constructor(file, "FDE", MaterialPath )


# # Close Lumerical 
# obj.Close()



# # When Lumerical is started you can use the inbuild Help to construct and run your simulation 
obj.Help()
obj.Help("Objects")
obj.Help({"Objects": 1})

obj.Help('Solvers')
obj.Help({"Solvers": 1})
obj.Help('Start Simulation')
obj.Help({"Start Simulation": 2})






# Waveguide Parameters
Parameters = {}
# "LiNbO3_20deg_X cut" is an Custom material added from me 
# "SiO2 (Glass) - Palik" is an Lumerical Material
Parameters['Substrate Height'] = 1e-6
Parameters['WG Length'] = 5e-6
Parameters['WG Height'] = 0.3e-6
Parameters['WG Width'] = 1e-6
Parameters['angle'] = 30 
Parameters['Slab Height'] = 0.3e-6
Parameters['Material'] = ["SiO2 (Glass) - Palik", 'LiNbO3_20deg_X cut']

# Solver Parameters
Parameters['y res'] = 0.1e-6
Parameters['z res'] = 0.1e-6
Parameters['Wavelength'] = 1.550e-6
Parameters["Solver y span"] = 5e-6 
Parameters["Solver z span"] = 5e-6

# Create S-Bend
obj.Waveguide(Parameters)

# Create Solver and set boundaries
obj.Solver("Waveguide", "FDE",  Parameters)

# Strat Simulations
obj.StartFDESimulation() 