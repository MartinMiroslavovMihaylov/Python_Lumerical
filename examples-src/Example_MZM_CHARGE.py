

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
obj = Charge(file, "CHARGE")


# # Close Lumerical 
# obj.Close()



# # When Lumerical is started you can use the inbuild Help to construct and run your simulation 
# obj.Help()
# obj.Help("Objects")
# obj.Help({"Objects": 12})

# obj.Help('Solvers')
# obj.Help({"Solvers": 17})
# obj.Help('Start Simulation')
# obj.Help({"Start Simulation": 4})






# MZM Parameters
Parameters = {}

Parameters['Substrate Height'] = 1e-6
Parameters["Optical"] = ["Au (Gold) - CRC", "SiO2 (Glass) - Palik"]
Parameters["Electrical"] = ["Air", "Au (Gold) - CRC", "LiNbO3 semiconductor - X/Y cut (Lithium Niobate)", "SiO2 (Glass) - Sze"]
Parameters['angle'] = 30
Parameters['Slab Height'] = 0.3e-6
Parameters['WG Height'] = 0.3e-6
Parameters['WG Width'] = 1e-6
Parameters['WG Length'] = 10e-6
Parameters["GND Electrodes Width"] = 5e-6
Parameters["Signal Electrodes Width"] = 3e-6
Parameters["Electrodes Height"] = 0.8e-6
Parameters["Gap"] = 1.4e-6
Parameters["Wavelength"] = 1.55e-6


# Create MZM
obj.MZM(Parameters) 

# Create Solver and set boundaries
obj.MZM_ChargeSolver(Parameters) 
obj.MZM_FEEMSolver(Parameters)


# # Strat Simulations
obj.StartCHARGESolver()
obj.StartFEEMSolver()

