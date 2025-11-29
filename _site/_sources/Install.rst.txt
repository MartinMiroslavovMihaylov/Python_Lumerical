
Install
=======

This Python Libraries can be used with Python Version 3.10.0 or newer.

For the Moment you will need the following Python Libraries. ::

   pip install numpy
   pip install pandas
   pip install matplotlib

   
Lumerical Installation and file location 
==============================================================

To use a Python Lumerical script, you will need a licensed and installed copy of Ansys Lumerical.
It doesn’t matter whether you are on Windows or Linux—you need to locate where Lumerical is installed and search in::

   Windows: C:/Program Files/Lumerical/v.../api/python

   Linux: .../Lumerical/v.../api/python

Inside, you will find lumapi.py. You need to provide your Python script with the path to where lumapi.py is located on your system.


Lumerical Materials Dataset
==============================================================

The scripts can be used with an external material library file. If you have one, for example
Material_File_Example.mdf, you can import it at the beginning of the script.
To see how, please refer to the next section Examples.

