Examples
========

Python code example how to start and call the Lumerical from Python. 
To call the following in your Python IDLEmß.

.. code-block:: python

	import os
	import os.path
	import sys
	# Import Constructor 
	import Constructor as Constructor
	# Import Help Menu
	from Constructor import Help
	#Import Loading Bar
	from Constructor import loadingBar
	# Import Logfile function extractor 
	from Constructor import Logfile
	
	# Give the lumpi.py files paths
	Path = "C:/Program Files/Lumerical/v221/api/python/lumapi.py"
	
	# Start the Lumerical FDTD API
	obj = Constructor.Constructor(file, "FDTD")
	
	
	## If you have an Material Library you wanna import you can use 
	#MaterialPath = "C:/.../Materials.mdf"
	#obj = Constructor.Constructor(file, "FDTD", MaterialPath )


Full Code Examples:
-------------------

.. toctree::
   :maxdepth: 1

   Example_Simple_Waveguide_FDE
   Example_MMI_EME
   Example_S_Bend_FDTD
   Example_MZM_CHARGE