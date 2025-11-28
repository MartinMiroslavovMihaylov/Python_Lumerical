
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


# =============================================================================
# Scripts Test
# =============================================================================



# # If you need Help with hot to import the Constructor run the line below 
# Help()

# Call object Constructor
obj = Constructor(file, "EME", MaterialPath )


# # Close Lumerical 
# obj.Close()



# # When Lumerical is started you can use the inbuild Help to construct and run your simulation 
# obj.Help()
# obj.Help("Objects")
# obj.Help({"Objects": 2})

# obj.Help('Solvers')
# obj.Help({"Solvers": 2})
# obj.Help('Start Simulation')
# obj.Help({"Start Simulation": 2})




# MMI Parameters
Parameters = {}
# "LiNbO3_20deg_X cut" is an Custom material added from me 
# "SiO2 (Glass) - Palik" is an Lumerical Material
Parameters['Material'] = ["SiO2 (Glass) - Palik", "LiNbO3_20deg_X cut"]
Parameters['Substrate Height'] = 1e-6
Parameters['MMI Width'] = 10e-6
Parameters['MMI Length'] = 20e-6
Parameters['angle'] = 30
Parameters['WG Height'] = 0.3e-6
Parameters['WG Width'] = 1e-6
Parameters['WG Length'] = 5e-6
Parameters['Position Offset'] = 3e-6
Parameters['Offset Input'] = 0
Parameters['Slab Height'] = 0.3e-6
Parameters['Taper'] = False
Parameters['Taper Length'] = 5e-6
Parameters['Taper Width'] = 3e-6



# Solver Parameters
Parameters['y res'] = 0.1e-6
Parameters['z res'] = 0.1e-6
Parameters['Wavelength'] = 1.550e-6
Parameters["Mode"] = "fundamental TE mode"
Parameters["Port Span"] = [2.5e-6, 2.5e-6, 2.5e-6]




# Create MMI
obj.MMI2x1(Parameters)  

# Create Solver and set boundaries
obj.Solver("MMI2x1", "EME", Parameters)

# Strat Simulations
obj.StartEMESimulation() 

        