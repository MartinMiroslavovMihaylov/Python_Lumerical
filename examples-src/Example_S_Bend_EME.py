
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
obj = Constructor(file, "FDTD", MaterialPath )


# # Close Lumerical 
# obj.Close()



# # When Lumerical is started you can use the inbuild Help to construct and run your simulation 
# obj.Help()
# obj.Help("Objects")
# obj.Help({"Objects": 8})

# obj.Help('Solvers')
# obj.Help({"Solvers": 13})
# obj.Help('Start Simulation')
# obj.Help({"Start Simulation": 2})






# S-Bend Parameters
Parameters = {}
# "LiNbO3_20deg_X cut" is an Custom material added from me 
# "SiO2 (Glass) - Palik" is an Lumerical Material
Parameters['Substrate Height'] = 1e-6
Parameters['WG Height'] = 0.3e-6
Parameters['WG Width'] = 1e-6
Parameters['angle'] = 30
Parameters['Slab Height'] = 0.3e-6
Parameters['Material'] = ["SiO2 (Glass) - Palik", 'LiNbO3_20deg_X cut']
Parameters["Wavelength"] = 1.55e-6
Parameters["x span"] = 30e-6
Parameters["y span"] = 20e-6
Parameters["poles"] = True

# Solver Parameters
Parameters['x res'] = 0.1e-6
Parameters["Mode"] = "fundamental TE mode"
Parameters["Port Span"] = [2.5e-6, 2.5e-6, 2.5e-6]


# Create S-Bend
obj.BendWaveguide(Parameters) 

# Create Solver and set boundaries
obj.Solver("BendWaveguide", "FDTD", Parameters)

# Strat Simulations
obj.StartFDTDSimulation() 