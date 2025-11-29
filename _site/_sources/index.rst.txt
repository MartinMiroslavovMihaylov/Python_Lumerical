.. Python Lumercial Costructor Script documentation master file, created by
   sphinx-quickstart on Sun Jan 28 11:01:32 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.
   
.. |date| date::

Last Updated on |date|

Welcome to Python Lumercial Costructor Script's documentation!
==============================================================

To install the library, please clone the git repo:
https://github.com/MartinMiroslavovMihaylov/Python_Lumerical.git
to an empty folder. Then you can install the Constructor lib with the setup.py file.

For installing Python_Lumerical use:

.. code-block:: python

    pip install git+https://github.com/MartinMiroslavovMihaylov/Python_Lumerical.git


Version numbers of this project might not always be incremented for minor changes.
To be sure, first uninstall:

.. code-block:: python

    pip uninstall Python_Lumerical
    pip install --no-cache-dir git+https://github.com/MartinMiroslavovMihaylov/Python_Lumerical.git


To install the library for development use:

.. code-block:: python

    git clone https://github.com/MartinMiroslavovMihaylov/Python_Lumerical.git
    cd Python_Lumerical
    pip install -e .

Option ``-e`` means **editable installation**.



How to use this Library
-----------------------

:doc:`Install`
   How to install.
   
:doc:`Examples/index`
   Code Examples.
   
:doc:`Extras/index`
	Code Help Menu, Loading Bar and Logger.

.. toctree::
   :maxdepth: 3
   :caption: Documentation

   Install
   Examples/index
   Extras/index
   Constructor