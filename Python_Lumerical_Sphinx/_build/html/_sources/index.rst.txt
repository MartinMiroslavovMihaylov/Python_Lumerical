.. Python Lumercial Costructor Script documentation master file, created by
   sphinx-quickstart on Sun Jan 28 11:01:32 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Python Lumercial Costructor Script's documentation!
==============================================================

To isntall the library please clone the git repo https://github.com/MartinMiroslavovMihaylov/Python_Lumerical.git 
to a empty folder. The you can install the Costructor lib with the setup.py file. Type in open Terminal:

.. code-block:: python

	pip install -e .

How to use the script
===================================================================

First you need to locate your Lumerical installation folder and find the lumpi.py file. 
For example under Windows you can find the file under "C:\Program Files\Lumerical\v221\api\python\lumapi.py".

Once located yoy will need to add to Lumerical an Material File. The file will have an ending 
"Material_File_Example.mdf". If you dont have an Material file you can skip this part. 

Finaly go to your Python IDLE and import the Constructor lib that you just installed:

.. code-block:: python

	import os
	import os.path
	import sys
	# Import Constructor 
	import Constructor as CN
	# Import Help Menu
	from Constructor import Help
	#Import Loading Bar
	from Constructor import loadingBar
	# Import Logfile function extractor 
	from Constructor import Logfile
	
	# Give the Material and lumpi.py files paths
	Path = "C:/Program Files/Lumerical/v221/api/python/lumapi.py"
	MaterialPath = "C:/.../Materials.mdf"
	
	
	# Start the Lumerical FDTD API
	obj = CN.Constructor(file, MaterialPath,  "FDTD")
	
	
How to use the Help Menu
===================================================================	

If you need Help with the functions you can use the build in Help guide.
Just type:

.. code-block:: python

	Help()
	
	
How to define Materials and Parameters
===================================================================	

The Materials are taken from the Lumerical build in or imported Materials database. 
The user can give Materials by them names as an python list of str. For Example:

.. code-block:: python

	Material = ["SiO2 (Glass) - Palik", 'LiNbO3_20deg_X cut', "SU-8"]

The parameters are given as dictionary to the Costructor. For example to create an 1x2_MMI 
you can use this code snip:

.. code-block:: python	

	#Create empty dictionary 
	Parameters = {}
	
	Parameters['Substrate Height'] = 1.05e-6
	Parameters['MMI Width'] = 3.6e-6
	Parameters['MMI Length'] = 18e-6
	Parameters['angle'] = 20
	Parameters['WG Height'] = 0.3e-6
	Parameters ['WG Width'] = 0.5e-6
	Parameters['WG Length'] = 5e-6
	Parameters['Offset Input'] = 0
	Parameters['Position Offset'] = 1.49e-6
	Parameters['Slab Height'] = 0.6e-6 - Parameters['WG Height']
	Parameters['Material'] = Material
	Parameters['Taper Length'] = 10e-6
	Parameters['Taper Width'] =  4e-6
	Parameters['Taper'] = False
	
	# Create the MMI 
	obj.MMI2x1(Parameters)
	
	
	
How to use the Logfile function 
===================================================================	

The Logfile function can be used after each simulation to generate an .txt file with all the parameters used 
for this specific simulation and the time stemp of the simulation. To use the function simply type:

.. code-block:: python

	# Save solver info to an variable
	solver = obj.ReturnLogInfo()
	# create an list of solver info and parameters used for this simulation
	data = [solver, Parameters]
	# Write the information to an .txt file in the given location 
	path_for_logFile = "C:/..../Simulation_Data/"
	Logfile(data, path_for_logFile)



How to use the loadingBar function 
===================================================================	

The laoding bar function will create an loading bar in your python terminal whitch will 
give you more information about the current status of the simulation (if you are doing an 
sweep over couple of parameters). To show you how to use the function an for loop will be 
used :

.. code-block:: python

	import time as t
	Lenght = np.arange(20e-6,30e-6,2e-6)
	for i in range(len(Lenght)):
		#Set the loading bar function
	    loadingBar(i, len(Lenght))
	    Parameters['Taper Length'] = Lenght[i]
	    obj.StraightWaveguide(Parameters)
	    t.sleep(5)
		obj.removeObject()



.. include:: Constructor.rst


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
