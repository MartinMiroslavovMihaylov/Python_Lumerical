Examples
========

In this section the Help Menu, Loading Bar and Logger will be demonstrated.


Help Menu
==========

This library includes a Help system to guide you through available functions, structures, solvers, and simulations. 
It helps you navigate the library more easily and understand what options are available.

**To use the Help system:**

.. code-block:: python
		
	from Constructor import HelpSystem

	# Show the main help menu
	HelpSystem.menu_help()



Loading Bar
===========

The Loading bar is an option adder to track long sweeps. In this example the loading bar is 
used to track the progress of my Waveguide lenght sweep. 


.. code-block:: python
    
	Lenght = np.arange(1e-6, 5e-6, 0.5e-6)
	S_Param = {}
	
	for i in range(len(Lenght)):
		loadingBar(i, len(Lenght))
		Parameters['WG Length'] = Lenght[i]
		obj.StraightWaveguide(Parameters)
		obj.Solver("StraightWaveguide", "EME",  Parameters)
		obj.StartEMESimulation()
		dic1 , dic2 = obj.ExtrtactData('EME', Parameters)
		S11 = abs(dic2['user s matrix'][0][0])**2
		S12 = abs(dic2['user s matrix'][0][1])**2
		S21 = abs(dic2['user s matrix'][1][0])**2
		S22 = abs(dic2['user s matrix'][1][1])**2
		S_Param["S12"].append(S12)
		S_Param["S21"].append(S21)
		obj.removeObject()
		
		

Logger
======

The logger will create a .txt file with the Lumerical simulation object, Solver, and parameters given.
To use it, the user can call:

.. code-block:: python

    from Constructor import Logfile
    Path = "C:/Downloads/"
    solverInfo = obj.ReturnLogInfo()
    data = [solverInfo, Parameters]
    Logfile(data, Path)

If you need help with how to use it, you can check the function here or in the Help Menu.