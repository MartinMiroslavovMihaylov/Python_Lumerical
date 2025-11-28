Examples
========

In this section the Help Menu, Loading Bar and Logger will be demonstrated.



Help Menu
===========

When the obj.Constructur() is called the user can call the obj.Help(). The Help function can be called with the following 
arguments::

	obj.Help('Objects')
	obj.Help('Solvers')
	obj.Help('Start Simulation')
	obj.Help('Results')
	obj.Help('Loading Bar')
	obj.Help('Log File')

When for example obj.Help('Objects') is called an Menu with list of objects will be given in the python terminal. 
If the used for example want to design an 1x2 MMI, the following can be called next obj.Help({"Objects": 2}). 
This will geve more information how to call the 1x2 MMI Function. What parameters need to be passed to the 1x2 MMI 
so that Lumerical can create it afterwords. After creating the object the "Solver" Help Menu can be choosen so that 
the user can get the proper Solver function and also give the needed parameters to the solver. Follow the Help 
Menu or go to the online documentory of this repo if you want to know more about the functions and parameters needed 
to call them. 


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
===========

The logger will create an .txt file with the Lumerical simulation object, Solver and parameters given. 
To use it, the user can call ::

.. code-block:: python
	from Constructor import Logfile                               
	Path = "C:/Downloads/"                                        
	solverInfo = obj.ReturnLogInfo()                              
	data = [solverInfo, Parameters]                               
	Logfile(data, Path) 
