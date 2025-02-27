# -*- coding: utf-8 -*-
"""
Created on Fri Feb 10 08:31:42 2023

@author: Martin.Mihaylov
"""


import numpy as np
import sys
import matplotlib.pyplot as plt



def loadingBar(count, total, size=1):
    """
    

    Parameters
    ----------
    count : int
        count is the loop iteration variable. For example count = i in case "for i in range(...):" is used.
    total : int
        Total number of bars that must be loaded for the respective process. For example total = 10 in case "for i in range(10):" 
    size : int, optional
        Defines how many "=" signs are used when the loading bar is printed in the console. The default is 1.

    Returns
    -------
    None.

    """
    percent = float(count+1) / float(total) * 100
    sys.stdout.write("\r" + str(int(count+1)).rjust(3, '0') + "/" +
                     str(int(total)).rjust(3, '0') + ' [' + '=' * int(percent / 10) *
                     size + ' ' * (10 - int(percent / 10)) * size + ']')


# Array with arange including the last element
def arange(start, stop, step=1, endpoint=True):
    arr = np.arange(start, stop, step)

    if endpoint and arr[-1] != stop:
        arr = np.append(arr, stop)

    return arr


class Constructor:

    # Init Programm
    def __init__(self, file, Mode, MaterialLib = None):
        '''
        

        Parameters
        ----------
        file : str
            Path to the lumerical python file.
        Mode : str
            Solver Type, can be set to FDTD or EME.
        MaterialLib : str, optional
            Path to the lumerical material lib. The default is None.

        Raises
        ------
        ValueError
            Error when wrong solver is given.

        Returns
        -------
        None.

        '''

        self.file = file
        self.MaterialLib = MaterialLib
        if self.MaterialLib is None:
            pass
        else:
            self.MaterialLib = MaterialLib
            
         
        # Check Python Version !!! No import imp module after python 3.11
        PyVersion = sys.version
      
        try:
            if PyVersion.split(" ")[0] > "3.11.0":
                import importlib.util
                import importlib.machinery

                def load_source(modname, filename):
                    loader = importlib.machinery.SourceFileLoader(modname, filename)
                    spec = importlib.util.spec_from_file_location(modname, filename, loader=loader)
                    module = importlib.util.module_from_spec(spec)
                    # The module is always executed and not cached in sys.modules.
                    # Uncomment the following line to cache the module.
                    # sys.modules[module.__name__] = module
                    loader.exec_module(module)
                    return module
                
                self.lumpai = load_source('lumapi', self.file)
            else:
                import imp
                self.lumpi = imp.load_source('lumapi', self.file)
        except ValueError as e:
            print(f"ValueError encountered: {e}")
        except ImportError as e:
            print(f"ImportError encountered: {e}")
        except Exception as e:
            print(f"An unexpected error occurred: {e}")
            
        self.Mode = Mode
        self.Struct = None
        self.SolverInfo = {}

        if self.Mode == "FDTD":
            self.FDTD()
            self.SolverInfo["Solver Used"] = "FDTD"
        elif self.Mode == "EME":
            self.EME()
            self.SolverInfo["Solver Used"] = "EME"
        else:
            raise ValueError("Non Valid Solver was choosen. Please pass on one of the two supported solvers ['FDTD' or 'EME']")

    # Close Programm
    def Close(self):
        '''
        

        Returns
        -------
        Close Lumerical GUI

        '''
        self.lum.close()
        print('Lumerical API is closed')


    # Remove the object and solver region
    def removeObject(self):
        '''
        Switch from simulation to layout mode.
        Remove all objects.

        Returns
        -------
        None.

        '''
        self.lum.switchtolayout()
        self.lum.deleteall()


    def Help(self, Subject = None):
        if Subject is None:
            print("i am in None section")
            raise ValueError("Help can be called with Help(str(subject)). Subject can be choosen from 'Objects', 'Solvers', 'Start Simulation', 'Results', 'Loading Bar' and 'Log File' !! ")
            
        # Check Type
        elif type(Subject) == dict:
            print("i am in Dict section")
            StartNumber = [1,2,3]
            ResultNumber = [1,2,3,4]
            StrucNumbers = [1,2,3,4,5,6,7,8,9,10,11]
            SolverNumber = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]
            
            listSub = ['Objects', 'Solvers', 'Start Simulation', 'Results', 'Loading Bar', 'Log File']
            key = list(Subject.keys())
            val = list(Subject.values())
            _help = HelpSubject()
            if key[0] == "Objects" and val[0] in StrucNumbers:
                _help.Help_Objects()
                _help.Structures(int(val[0]))
            elif key[0] == "Solvers" and val[0] in SolverNumber:
                _help.Help_Solvers()
                _help.Solvers(int(val[0]))
            elif key[0] == "Start Simulation" and val[0] in StartNumber:
                _help.Help_StartSimulation()
                _help.StartSolver(int(val[0]))
            elif key[0] == "Results" and val[0] in ResultNumber:
                _help.Help_Results()
                _help.Result_extraction(int(val[0]))
            elif key[0] == "Loading Bar":
                _help.Help_LoadingBar()
            elif key[0] == "Remove Object":
                _help.Help_RemoveObject()
            elif key[0] == "Log File":
                _help.Help_LogFile()  
            elif key[0] not in listSub:
                raise ValueError("Help can be called with Help(str(subject)). Subject can be choosen from 'Objects', 'Solvers', 'Start Simulation', 'Results', 'Loading Bar' and 'Log File' !! ")
            #elif val[0] not in StartNumber or ResultNumber or StrucNumbers or StrucNumbers:
            
        elif type(Subject) == str :
            print("i am in STR section")
            _help = HelpSubject()
            if Subject == "Objects":
                _help.Help_Objects() 
            elif Subject == "Solvers":
                _help.Help_Solvers()
            elif Subject == "Start Simulation":
                _help.Help_StartSimulation()
            elif Subject == "Results":
                _help.Help_Results()
            elif Subject == "Loading Bar":
                _help.Help_LoadingBar()
            elif Subject == "Remove Object":
                _help.Help_RemoveObject()
            elif Subject == "Log File":
                _help.Help_LogFile()
            else:
                raise ValueError("Help can be called with Help(str(subject)). Subject can be choosen from 'Objects', 'Solvers', 'Start Simulation', 'Results', 'Loading Bar' and 'Log File' !! ")



                
            
         
        
     
      
            





    # Choose Solver! Only FDTD and EME supported
    def FDTD(self):
        '''
        
        Calls the FDTD Solver.

        Returns
        -------
        None.

        '''
        self.lum = self.lumpai.FDTD()
        if self.MaterialLib is None:
            pass
        else:
            self.lum.importmaterialdb(self.MaterialLib)
        print('Lumerical FDTD API is started')
        
        

    def EME(self):
        '''
        Calls the EME solver.

        Returns
        -------
        None.

        '''
        self.lum = self.lumpai.MODE()
        if self.MaterialLib is None:
            pass
        else:
            self.lum.importmaterialdb(self.MaterialLib)
        print('Lumerical EME API is started')
    


    #   Solvers Function
    def Solver(self, Structure, Type, Parameters):
        """
        

        Parameters
        ----------
        Structure : str
            Object that will be simulated. 
        Type : str
            Solver Type.
        Parameters : dict
            Dictionary of solver parameters.

        Raises
        ------
        ValueError
            Error when wrong solver type selected.
        Returns
        -------
        None.

        """
        if Type == "FDTD":
            self.FDTD_Solver(Structure, Parameters)
        elif Type == "EME":
            self.EME_Solver(Structure, Parameters)
        elif Type == "FDE":
            self.FDE_Solver(Structure, Parameters)
        else:
            raise ValueError(
                "Invalid Solver Type! Possible Options are FDE - For Waveguides, EME - for MMI , directional Couplers and FDTD - for MMI.")



    def StartSimFunction(self, Type):
        '''
        

        Parameters
        ----------
        Type : str
            Solver Type.

        Raises
        ------
        ValueError
            Error when wrong solver type selected.

        Returns
        -------
        None.

        '''
        if Type == "FDTD":
            self.StartFDTDSimulation()
            self.SolverInfo["Simulation Type"]  = "FDTD Simulation"
        elif Type == "EME":
            self.StartEMESimulation()
            self.SolverInfo["Simulation Type"]  = "EME Simulation"
        elif Type == "FDE":
            self.StartFDESimulation()
            self.SolverInfo["Simulation Type"]  = "FDE Simulation"
        else:
            raise ValueError("Invalid simulation type. Possibpe options for Type are: 'FDTD', 'EME', 'FDE' !")



    def ExtrtactData(self, Type, Parameters):
        '''
        

        Parameters
        ----------
        Type : str
            Solver Type.
        Parameters : boolen
            Need to be set True if S-Parameter analysis in FDTD solver is wanted. Otherwise False.

        Raises
        ------
        ValueError
            DESCRIPTION.

        Returns
        -------
        TYPE
            Error when wrong structure was selected.

        '''
        if Type == "FDTD":
            if self.Struct == "MMI2x1":
                S_Param3, Power3 = self.ExtractFDTDResults(3, Parameters)
                return S_Param3, Power3
            elif self.Struct == "MMI2x2":
                S_Param4 = self.ExtractFDTDResults(4, Parameters)
                return S_Param4
            elif self.Struct == "DirectionalCoupler":
                S_Param4, Power4 = self.ExtractFDTDResults(4, Parameters)
                return S_Param4, Power4
            elif self.Struct == "ArcWaveguide":
                S_Param2, Power2 = self.ExtractFDTDResults(2, Parameters)
                return S_Param2, Power2
            elif self.Struct == "BendWaveguide":
                S_Param22, Power22 = self.ExtractFDTDResults(2, Parameters)
                return S_Param22, Power22
            elif self.Struct == "StraightWaveguide":
                S_Param22, Power22 = self.ExtractFDTDResults(2, Parameters)
                return S_Param22, Power22
            else:
                raise ValueError("Incorect Structure was given. Possible structures are 'MMI2x1', 'MMI2x2', 'DirectionalCoupler', 'ArcWaveguide' and 'BendWaveguide' !")
        elif Type == "EME":
            OptionsMonitor, EME_Data = self.ExtractEMEResults('monitor', Type)
            return OptionsMonitor, EME_Data
        elif Type == "FDE":
            TEModes, TEData, TMModes, TMData = self.ExtractFDEResultsExtendet(Parameters['Effective Index'])
            return TEModes, TEData, TMModes, TMData
        else:
            raise ValueError(
                "Invalid solver Type was given. To extract data from the simulation please give a correct simulation Type! Possible Types are: 'FDTD', 'EME' and 'FDE'.")



    def OverlapFDEModes(self, Parameters):
        '''
        

        Parameters
        ----------
        Parameters : dict
            Dictionary with parameters needed for the FDE analysis.

        Returns
        -------
        None.

        '''
        self.Waveguide(Parameters)
        self.FDE_Solver("Waveguide", Parameters)
        self.StartFDESimulation()
        TEModes, TMModes = self.ExtractFDEModes(EffIndexValue=Parameters["Effective Index"])
        self.CoppyDcard(list(TEModes.keys())[0], list(TMModes.keys())[0])
        self.removeObject()
        print(f"TE Mode {list(TEModes.keys())[0]} and TM Mode {list(TMModes.keys())[0]} have been copyed to the Lumerical Dcard.")

    # =============================================================================
    # Selection of possible Functions
    # =============================================================================



    def FDTD_Solver(self, Structure, Parameters):
        '''
        

        Parameters
        ----------
        Structure : str
            Structure to be simulated.
        Parameters : dict
            Dictionary with solver parameters.

        Raises
        ------
        ValueError
            Error when wrong structure was selected.

        Returns
        -------
        None.

        '''
        if Structure == "MMI2x1":
            self.Struct = "MMI2x1"
            self.setMMI2x1FDTDSolver(Parameters)
            self.SolverInfo["Simulated Object"] = "MMI2x1"
        elif Structure == "MMI2x2":
            self.Struct = "MMI2x2"
            self.setMMI2x2FDTDSolver(Parameters)
            self.SolverInfo["Simulated Object"] = "MMI2x2"
        elif Structure == "DirectionalCoupler":
            self.Struct = "DirectionalCoupler"
            self.setDCFDTDSolver(Parameters)
            self.SolverInfo["Simulated Object"] = "Directional Coupler"
        elif Structure == "BendWaveguide":
            self.Struct = "BendWaveguide"
            self.setBendWaveguideFDTDSolver(Parameters)
            self.SolverInfo["Simulated Object"] = "Bend Waveguide"
        elif Structure == "ArcWaveguide":
            self.Struct = "ArcWaveguide"
            self.setArcWaveguideFDTDSolver(Parameters)
            self.SolverInfo["Simulated Object"] = "Arc Waveguide"
        elif Structure == "StraightWaveguide":
            self.Struct = "StraightWaveguide"
            self.setStraightWaveguideFDTDSolver(Parameters)
            self.SolverInfo["Simulated Object"] = "Straight Waveguide"
        # elif Structure == "CascadetMMI":
        #     self.Struct = "CascadetMMI"
        #     self.setCascadetMMIFDTDSolver(Parameters, SpaceX, SpaceY)
        elif Structure == "InverseTaper":
            self.Strcut = "InverseTaper"
            self.setInverseTaperFDTDSolver(Parameters)
            self.SolverInfo["Simulated Object"] = "Inverse Taper"
        elif Structure == "LenseFiber":
            self.Struct = "LenseFiber"
            self.setLenseFiberFDTD(Parameters)
            self.SolverInfo["Simulated Object"] = "Lense Fiber"
        elif Structure == "GratingCoupler":
            self.Struct = "GratingCoupler"
            self.setGratingCouplerFDTDSolver(Parameters)
            self.SolverInfo["Simulated Object"] = "GratingCoupler"
        elif Structure == "RingGratingCoupler":
            self.Struct = "RingGratingCoupler"
            self.setRingGratingCouplerFDTDSolver(Parameters)
            self.SolverInfo["Simulated Object"] = "RingGratingCoupler"   
        else:
            raise ValueError("Invalid Strucute for FDTD Solver is selected. Possible Strucures are MMI2x1, MMI2x2, DirectionalCoupler, BendWaveguide, ArcWaveguide, StraightWaveguide, InverseTaper, GratingCoupler and RingGratingCoupler")


    def EME_Solver(self, Structure, Parameters):
        '''
        

        Parameters
        ----------
        Structure : str
            Structure to be simulated.
        Parameters : dict
            Dictionary with solver parameters.

        Raises
        ------
        ValueError
            Error when wrong structure was selected.

        Returns
        -------
        None.

        '''
        if Structure == "MMI2x1":
            self.setMMI2x1EMESolver(Parameters)
            self.SolverInfo["Simulated Object"] = "MMI2x1"
        elif Structure == "MMI2x2":
            self.setMMI2x2EMESolver(Parameters)
            self.SolverInfo["Simulated Object"] = "MMI2x2"
        elif Structure == "DirectionalCoupler":
            self.setDCEMESolver(Parameters)
            self.SolverInfo["Simulated Object"] = "Directional Coupler"
        elif Structure == "WDM":
            self.setWDMEMESolver(Parameters)
            self.SolverInfo["Simulated Object"] = "Angled Wavelength Division Multiplexer"
        elif Structure == "StraightWaveguide":
            self.setStraightWaveguideEMESolver(Parameters)
            self.SolverInfo["Simulated Object"] = "Straight Waveguide"
        elif Structure == "InverseTaper":
            self.setInverseTaperEMESolver(Parameters)
            self.SolverInfo["Simulated Object"] = "Inverse Tapers"
        else:
            raise ValueError("Invalid Structure for EME Solver is selected. Possible Strucutres are MMI2x1, MMI2x2, DirectionalCoupler, WDM, StraightWaveguide andInverseTaper ")



    def FDE_Solver(self, Structure, Parameters):
        '''
        

        Parameters
        ----------
        Structure : str
            Structure to be simulated.
        Parameters : dict
            Dictionary with solver parameters.

        Raises
        ------
        ValueError
            Error when wrong structure was selected.

        Returns
        -------
        None

        '''
        if Structure == "Waveguide":
            self.setWaveguideFDESolver(Parameters)
            self.SolverInfo["Simulated Object"] = "Waveguide"
        else:
            raise ValueError("Invalid Structure for FDE Solver is selected. Possible Structure is Waveguide.")




    def ReturnLogInfo(self):
        '''
        

        Returns
        -------
        TYPE
            Solver log information from last simulation.

        '''
        return self.SolverInfo



    def Script(self):
        '''
        

        Returns
        -------
        myscript : str
            small script syntax for creating Tapers tooken from Lumerical . It is used only for check comperisons.

        '''
        myscript = 'delta_w = 2*thickness*tan((angle_side)*pi/180); \n'
        myscript = myscript + '?"width_l = " + num2str(width_l); \n'
        myscript = myscript + '?"width_r = " + num2str(width_r) + endl; \n'
        myscript = myscript + 'width_top_l = width_l - (1-hfrac_ref)*delta_w; \n'
        myscript = myscript + 'width_top_r = width_r - (1-hfrac_ref)*delta_w; \n'
        myscript = myscript + '?"width_top_l = " + num2str(width_top_l); \n'
        myscript = myscript + '?"width_top_r = " + num2str(width_top_r) + endl; \n'
        myscript = myscript + 'width_bot_l = width_l + hfrac_ref*delta_w; \n'
        myscript = myscript + 'width_bot_r = width_r + hfrac_ref*delta_w; \n'
        myscript = myscript + '?"width_bot_l = " + num2str(width_bot_l); \n'
        myscript = myscript + '?"width_bot_r = " + num2str(width_bot_r) + endl; \n'
        myscript = myscript + 'if ((hfrac_ref>1) or (hfrac_ref<0)){?"Error: hfrac_ref must be between 0 and 1.";break;} \n'
        myscript = myscript + 'if ((width_top_l<0) or (width_top_r<0) or (width_bot_l<0) or (width_bot_r<0)){?"Error: width and angle values are not correct.";break;} \n'
        myscript = myscript + 'zmin = -thickness/2; \n'
        myscript = myscript + 'zmax = thickness/2; \n'
        myscript = myscript + 'xmin = -len/2; \n'
        myscript = myscript + 'xmax = len/2; \n'
        myscript = myscript + 'ymin_bot_l = -width_bot_l/2; \n'
        myscript = myscript + 'ymax_bot_l = width_bot_l/2; \n'
        myscript = myscript + 'ymin_bot_r = -width_bot_r/2; \n'
        myscript = myscript + 'ymax_bot_r = width_bot_r/2; \n'
        myscript = myscript + 'ymin_top_l = -width_top_l/2; \n'
        myscript = myscript + 'ymax_top_l = width_top_l/2; \n'
        myscript = myscript + 'ymin_top_r = -width_top_r/2; \n'
        myscript = myscript + 'ymax_top_r = width_top_r/2; \n'
        myscript = myscript + 'vtx=    [xmin,ymin_bot_l,zmin; xmax,ymin_bot_r,zmin; xmax,ymax_bot_r,zmin; xmin,ymax_bot_l,zmin; xmin,ymin_top_l,zmax; xmax,ymin_top_r,zmax; xmax,ymax_top_r,zmax; xmin,ymax_top_l,zmax];   \n'
        myscript = myscript + 'a = cell(6); \n'
        myscript = myscript + 'for(i = 1:6){ a{i} = cell(1);}   \n'
        myscript = myscript + 'a{1}{1} = [1,4,3,2];   \n'
        myscript = myscript + 'a{2}{1} = [1,2,6,5];   \n'
        myscript = myscript + 'a{3}{1} = [2,3,7,6];   \n'
        myscript = myscript + 'a{4}{1} = [3,4,8,7];   \n'
        myscript = myscript + 'a{5}{1} = [1,5,8,4];   \n'
        myscript = myscript + 'a{6}{1} = [5,6,7,8];   \n'
        myscript = myscript + 'addplanarsolid(vtx,a);   \n'
        myscript = myscript + 'if (material=="<Object defined dielectric>"){setnamed("solid", "index",index);}   \n'
        myscript = myscript + 'else{setnamed("solid", "material",material);}   \n'

        return myscript




# =============================================================================
# Structures
# =============================================================================


    def Waveguide(self, Parameters):
        


        '''
        Parameters
        ----------
        Parameters : Dictionary
            Dictionary with all the data needet for the Wavaguide.
            
            Parameters['Substrate Height'] : int7float
                Height of the substrate.
            Parameters['WG Length']: int/float
                Length of the Waveguide.
            Parameters['WG Height'] : int/float
                Height of the Waveguide
            Parameters['WG Width'] : int/float
                Width of the Waveguide
            Parameters['angle'] : int
                Angle of the Waveguide Walls. it is calculated WG_angle = 90 - angle.
                For anfle = 90 we get a perfect rect!
            Parameters['Slab Height'] : int/float
                Slab height.
            Parameters['Material'] : list of str
                List of Materials. The list should be with names (str) of a valid Lumerical materials.
                Check the names in Lumerical Materials viewer.
                The List of materials must contain at least 2 materials! 
                Parameters['Material'] = ['Cladding/Substrat', 'Object Material'].
                For Example: Parameters['Material'] = ["SiO2 (Glass) - Palik", 'LiNbO3_20deg_X cut'].

        Raises
        ------
        ValueError
            Value Error

        Returns
        -------
        None.
        '''

        WG_Length = Parameters['WG Length']
        WG_Height = Parameters['WG Height']
        WG_Width = Parameters['WG Width']
        angle = Parameters['angle']
        Slab_Height = Parameters['Slab Height']
        Material = Parameters['Material']
        SubstrateHight = Parameters['Substrate Height']
        
        
        if "Cladding" in list(Parameters.keys()):
            Cladding = True
        else:
            Cladding = False


        # Material defginition
        if len(Material) < 2:
            raise ValueError("List of materials must contain at least 2 materials!, Parameters['Material'] = ['Cladding/Substrat', 'Object Material']")
        else:
            MaterialSub = Material[0]
            MaterialClad = Material[0]
            MaterialSlab = Material[1]
            MaterialWG = MaterialSlab


        maxWGL = WG_Length / 2
        minWGL = -WG_Length / 2
        min_subW = -15e-6
        max_subW = 15e-6
        # Procentage deviation of the height between slab and waveguide
        min_slabH = 0
        max_slabH = Slab_Height
        min_WGH = max_slabH
        max_WGH = WG_Height
        min_BoxH = -4e-6


        # Triangle EQ for waveguide Width
        x = abs(max_WGH / (np.cos((angle) * np.pi / 180)))  # in Radians
        extention = np.sqrt(x ** 2 - max_WGH ** 2)
        WG_W = WG_Width + 2 * extention

        if Slab_Height == 0:
            # Waveguide
            self.lum.addwaveguide()
            self.lum.set("name", 'Waveguide')
            self.lum.set("x", 0)
            self.lum.set("y", 0)
            self.lum.set("z", WG_Height/2)
            self.lum.set("base width", WG_W)
            self.lum.set("base height", WG_Height)
            self.lum.set("base angle", 90 - angle)
            pole = np.array([[maxWGL, 0], [minWGL, 0]])
            self.lum.set("poles", pole)
            self.lum.set("material", MaterialWG)

            # Mesh
            self.lum.addmesh()
            self.lum.set("name", 'Mesh_Waveguide')
            self.lum.set('based on a structure', 1)
            self.lum.set('structure', 'Waveguide')
            self.lum.set('override x mesh', 0)
            self.lum.set('dy', 0.01e-6)
            self.lum.set('dz', 0.005e-6)

            # Add SiO2 Slubstrate
            self.lum.addrect()
            self.lum.set("name", "Silicon")
            self.lum.set("y min", min_subW)
            self.lum.set("y max", max_subW)
            self.lum.set("z", -SubstrateHight/2)
            self.lum.set("z span",SubstrateHight)
            self.lum.set("x min", minWGL)
            self.lum.set("x max", maxWGL)
            self.lum.set("material", MaterialSub)

            if Cladding == True:
                # create_cover
                self.lum.addrect()
                self.lum.set("name", "cladding")
                self.lum.set("y min", min_subW)
                self.lum.set("y max", max_subW)
                self.lum.set('z', max_WGH / 2 + min_WGH)
                self.lum.set('z span', 2 * WG_Height)
                self.lum.set("x min", minWGL)
                self.lum.set("x max", maxWGL)
                self.lum.set("material", MaterialClad)
                self.lum.set("override mesh order from material database", True)
                self.lum.set("mesh order", 4)
                self.lum.set("alpha", 0.7)
            else: 
                pass


        else:
            # Slab
            self.lum.addrect()
            self.lum.set("name", "LN_slab")
            self.lum.set("y min", min_subW)
            self.lum.set("y max", max_subW)
            self.lum.set("z min", min_slabH)
            self.lum.set("z max", max_slabH)
            self.lum.set("x min", minWGL)
            self.lum.set("x max", maxWGL)
            self.lum.set("material", MaterialWG)

            # Waveguide
            self.lum.addwaveguide()
            self.lum.set("name", 'Waveguide')
            self.lum.set("x", 0)
            self.lum.set("y", 0)
            self.lum.set("z", max_WGH / 2 + min_WGH)
            self.lum.set("base width", WG_W)
            self.lum.set("base height", max_WGH)
            self.lum.set("base angle", 90 - angle)
            pole = np.array([[maxWGL, 0], [minWGL, 0]])
            self.lum.set("poles", pole)
            self.lum.set("material", MaterialWG)

            # Mesh
            self.lum.addmesh()
            self.lum.set("name", 'Mesh_Waveguide')
            self.lum.set('based on a structure', 1)
            self.lum.set('structure', 'Waveguide')
            # self.lum.set('dx', 0.01e-6)
            self.lum.set('dy', 0.01e-6)
            self.lum.set('dz', 0.005e-6)

            # Add SiO2 Box
            self.lum.addrect()
            self.lum.set("name", "Silicon")
            self.lum.set("y min", min_subW)
            self.lum.set("y max", max_subW)
            self.lum.set("z min", min_BoxH)
            self.lum.set("z max", min_slabH)
            self.lum.set("x min", minWGL)
            self.lum.set("x max", maxWGL)
            self.lum.set("material", MaterialSub)
            
            if Cladding == True:
                # Create Cladding
                self.lum.addrect()
                self.lum.set("name", "cladding")
                self.lum.set("y min", min_subW)
                self.lum.set("y max", max_subW)
                self.lum.set('z', max_WGH / 2 + min_WGH)
                self.lum.set('z span', 2 * WG_Height)
                self.lum.set("x min", minWGL)
                self.lum.set("x max", maxWGL)
                self.lum.set("material", MaterialClad)
                self.lum.set("override mesh order from material database", True)
                self.lum.set("mesh order", 4)
                self.lum.set("alpha", 0.7)
            else:
                pass
            
            





    def StraightWaveguide(self, Parameters):

        '''

        Parameters
        ----------
        Parameters : Dictionary
            Dictionary with all the data needet for the Straight Wavaguide.
            
            Parameters['Substrate Height'] : int float
                Substrate Height
            Parameters['WG Height'] : int/float
                Height of the Waveguide
            Parameters['WG Width'] : int/float
                Width of the Waveguide
            Parameters['WG Length'] : int/float
                Length of the Waveguide
            Parameters['Taper'] : boolen
                If Taper == False, only Waveguiedes will be constructed.
                If Taper == True, only an single Taper will be constructed
            Parameters['Taper Length'] : int/float
                Taper Length
            Parameters['Taper Width'] : int/float
                Taper backside width, the frontside width is the waveguide width
            Parameters['angle'] : int
                Angle of the Waveguide Walls. it is calculated WG_angle = 90 - angle.
                For anfle = 90 we get a perfect rect!
            Parameters['Slab Height'] : int/float
                Slab height.
            Parameters['Material'] : list
                List of Materials. The list should be with names (str) of a valid Lumerical materials.
                Check the names in Lumerical Materials viewer.
                The List of materials must contain at least 2 materials! 
                Parameters['Material'] = ['Cladding/Substrat', 'Object Material'].
                For Example: Parameters['Material'] = ["SiO2 (Glass) - Palik", 'LiNbO3_20deg_X cut'].
            Parameters["Wavelength"] : int/float
                Wavelength
            Parameters["Waveguide Angle"] : int/float
                Bending angle of the Waveguide. Set it to 0 to simulate the straight waveguide.
                If Waveguide Angle is different then 0, then the straight waveguide will be tilted
                at the choosen degrees.
            Parameters["Cladding"] : anything, optional
                This function will check if you have set Parameters["Cladding"] to anaything, for example "Parameters["Cladding"]=1" 
                and if so it will put cladding over your structure. If the user didnt give the "Cladding" as dictionary key no cladding 
                will be set.
            Parameters["Taper Type"] : anything, optional
                This function will check if you have set Parameters["Taper Type"] to anaything, for example "Parameters["Taper Type"]=1" 
                and if so it will design an Inverse Taper Structure with no Cladding. Here the option "Cladding" is not active and will be ignored.
                If the user didnt give the "Taper Type" as dictionary key, then an normal taper structure will be simulated.
                If Taper Type is selected the user need to provide additional information:
                    TaperWidthF = Parameters['PWB Taper Width Front']
                    TaperWidthB = Parameters['PWB Taper Width Back']
                    TaperHightB = Parameters['PWB Taper Hight Back']
                    TaperHightF = Parameters['PWB Taper Hight Front']
                    TaperLength_PWB = Parameters['PWB Taper Length']
                   
             
        Raises
        ------
        ValueError
            Value Error

        Returns
        -------
        None.

        '''


        Substrate_Height = Parameters['Substrate Height']
        WG_Height = Parameters['WG Height']
        WG_Width = Parameters['WG Width']
        WG_Length = Parameters["WG Length"]
        Taper = Parameters['Taper']
        TaperLength = Parameters['Taper Length']
        TaperWidth = Parameters['Taper Width']
        angle = Parameters['angle']
        Slab_Height = Parameters['Slab Height']
        Material = Parameters['Material']
        WaveLength = Parameters["Wavelength"]
        Side_Angle = Parameters["Waveguide Angle"]
        if "Cladding" in list(Parameters.keys()):
            Cladding = True
        else:
            Cladding = False
            
        if "Taper Type" in list(Parameters.keys()):
                TaperType = "Inverse"
                TaperWidthF = Parameters['PWB Taper Width Front']
                TaperWidthB = Parameters['PWB Taper Width Back']
                TaperHightB = Parameters['PWB Taper Hight Back']
                TaperHightF = Parameters['PWB Taper Hight Front']
                TaperLength_PWB = Parameters['PWB Taper Length']
        else:
            TaperType = "Normal"
        
            
        


        # Material definition
        if len(Material) < 2:
            raise ValueError(
                "List of materials must contain at least 2 materials!, Parameters['Material'] = ['Cladding/Substrat', 'Object Material']")
        else:
            MaterialSub = Material[0]
            MaterialClad = Material[0]
            MaterialSlab = Material[1]
            MaterialWG = MaterialSlab





        # Device specifications
        Device_Width = 2*WG_Length #+ 1e-6 * 2  

        # creating the substrate
        max_subH = Substrate_Height
        min_subH = -Substrate_Height



        if Taper == False:

            if Side_Angle == 0:

                self.lum.addrect()
                self.lum.set("name", "Substrate")
                self.lum.set("y", 0)
                self.lum.set("y span", Device_Width)
                self.lum.set("z min", min_subH)
                self.lum.set("z max", max_subH)
                self.lum.set("x", WG_Length / 2)
                self.lum.set("x span", WG_Length + WG_Width)
                self.lum.set("material", MaterialSub)
                
                
                if Slab_Height != 0:
                    
                    # creating the thin film
                    min_slabH = max_subH
                    max_slabH = max_subH + Slab_Height

                    self.lum.addrect()
                    self.lum.set("name", "Slab")
                    self.lum.set("y", 0)
                    self.lum.set("y span", Device_Width)
                    self.lum.set("z min", min_slabH)
                    self.lum.set("z max", max_slabH)
                    self.lum.set("x", WG_Length / 2)
                    self.lum.set("x span", WG_Length + WG_Width)
                    self.lum.set("material", MaterialSlab)

                    z_Offset = max_slabH + WG_Height / 2
                    
                    self.lum.select("Slab")
                    self.lum.addtogroup("Straight Waveguide")
                    
                    if Cladding == True:
                        # create_cover
                        self.lum.addrect()
                        self.lum.set("name", "cladding")
                        self.lum.set("material", MaterialClad)
                        self.lum.set("y", 0)
                        self.lum.set("y span", Device_Width)
                        self.lum.set("z min", max_slabH)
                        self.lum.set("z max", max_slabH*2)
                        self.lum.set("x", WG_Length / 2)
                        self.lum.set("x span", WG_Length + 1e-6)
                        self.lum.set("override mesh order from material database", True)
                        self.lum.set("mesh order", 4)
                        self.lum.set("alpha", 0.7)
                        self.lum.select("cladding")
                        self.lum.addtogroup("Straight Waveguide")
                    else:
                        pass

                    
                else:
                    z_Offset = max_subH + WG_Height / 2
                    
                    if Cladding == True:
                        # create_cover
                        self.lum.addrect()
                        self.lum.set("name", "cladding")
                        self.lum.set("material", MaterialClad)
                        self.lum.set("y", 0)
                        self.lum.set("y span", Device_Width)
                        self.lum.set("z min", max_subH)
                        self.lum.set("z max", max_subH*2)
                        self.lum.set("x", WG_Length / 2)
                        self.lum.set("x span", WG_Length + 1e-6)
                        self.lum.set("override mesh order from material database", True)
                        self.lum.set("mesh order", 4)
                        self.lum.set("alpha", 0.7)
                        self.lum.select("cladding")
                        self.lum.addtogroup("Straight Waveguide")
                    else:
                        pass
                        
                        
                    
                # Triangle EQ for waveguide Width
                x = abs(WG_Height / (np.cos((angle) * np.pi / 180)))  # in Radians
                extention = np.sqrt(x ** 2 - WG_Height ** 2)
                WG_W = WG_Width + 2 * extention
                



                names = ["Waveguide"]

                self.lum.addwaveguide()
                self.lum.set("name", names[0])
                self.lum.set("x", 0)
                self.lum.set("y", 0)
                self.lum.set("z", z_Offset)
                self.lum.set("base width", WG_W)
                self.lum.set("base height", WG_Height)
                self.lum.set("base angle", 90 - angle)
                pole = np.array([[-0.5e-6, 0], [WG_Length+0.5e-6, 0]])
                self.lum.set("poles", pole)
                self.lum.set("material", MaterialWG)
                self.lum.set("first axis", 'z')
                self.lum.set("rotation 1", Side_Angle)
                
                self.lum.select("Waveguide")
                self.lum.addtogroup("Straight Waveguide")
                
                

                
                

            else:
                self.lum.addrect()
                self.lum.set("name", "Substrate")
                self.lum.set("y", WG_Length / 4)
                self.lum.set("y span", Device_Width)
                self.lum.set("z min", min_subH)
                self.lum.set("z max", max_subH)
                self.lum.set("x", WG_Length / 2)
                self.lum.set("x span", WG_Length + WG_Width)
                self.lum.set("material", MaterialSub)


                if Slab_Height == 0:
                    z_Offset = max_subH + WG_Height / 2
                
                else:
                
                    # creating the thin film
                    min_slabH = max_subH
                    max_slabH = max_subH + Slab_Height

                    self.lum.addrect()
                    self.lum.set("name", "Slab")
                    self.lum.set("y", WG_Length / 4)
                    self.lum.set("y span", Device_Width)
                    self.lum.set("z min", min_slabH)
                    self.lum.set("z max", max_slabH)
                    self.lum.set("x", WG_Length / 2)
                    self.lum.set("x span", WG_Length + WG_Width)
                    self.lum.set("material", MaterialSlab)

                    z_Offset = max_slabH + WG_Height / 2
                    self.lum.select("Slab")
                    self.lum.addtogroup("Straight Waveguide")
                    
                    
                    
                    
                    
                # Triangle EQ for waveguide Width
                x = abs(WG_Height / (np.cos((angle) * np.pi / 180)))  # in Radians
                extention = np.sqrt(x ** 2 - WG_Height ** 2)
                WG_W = WG_Width + 2 * extention



                names = ["Waveguide"]

                self.lum.addwaveguide()
                self.lum.set("name", names[0])
                self.lum.set("x", 0)
                self.lum.set("y", 0)
                self.lum.set("z", z_Offset)
                self.lum.set("base width", WG_W)
                self.lum.set("base height", WG_Height)
                self.lum.set("base angle", 90 - angle)
                pole = np.array([[-1e-6, 0], [WG_Length+1e-6, 0]])
                self.lum.set("poles", pole)
                self.lum.set("material", MaterialWG)
                self.lum.set("first axis", 'z')
                self.lum.set("rotation 1", Side_Angle)
                
                self.lum.select("Waveguide")
                self.lum.addtogroup("Straight Waveguide")

                if Cladding == True:
                    # create_cover
                    self.lum.addrect()
                    self.lum.set("name", "cladding")
                    self.lum.set("material", MaterialClad)
                    self.lum.set("y", WG_Length / 4)
                    self.lum.set("y span", Device_Width)
                    self.lum.set("z min", min_slabH)
                    self.lum.set("z max", min_slabH*2)
                    self.lum.set("x", WG_Length / 2)
                    self.lum.set("x span", WG_Length + WG_Width)
                    self.lum.set("override mesh order from material database", True)
                    self.lum.set("mesh order", 4)
                    self.lum.set("alpha", 0.7)
                    self.lum.select("cladding")
                    self.lum.addtogroup("Straight Waveguide")
                else:
                    pass

        elif Taper == True:

            # Device specifications
            Device_Width = 2*TaperLength + WaveLength * 2  # MMI_Width

            # creating the substrate
            max_subH = Substrate_Height
            min_subH = -Substrate_Height

            self.lum.addrect()
            self.lum.set("name", "Substrate")
            self.lum.set("y", 0)
            self.lum.set("y span", Device_Width)
            self.lum.set("z min", min_subH)
            self.lum.set("z max", max_subH)
            self.lum.set("x", 0)
            self.lum.set("x span", TaperLength + WG_Width)
            self.lum.set("material", MaterialSub)
            
            if Slab_Height != 0:
                
                # creating the thin film
                min_slabH = max_subH
                max_slabH = max_subH + Slab_Height
                
                self.lum.addrect()
                self.lum.set("name", "Slab")
                self.lum.set("y", 0)
                self.lum.set("y span", Device_Width)
                self.lum.set("z min", min_slabH)
                self.lum.set("z max", max_slabH)
                self.lum.set("x", 0)
                self.lum.set("x span", TaperLength + WG_Width)
                self.lum.set("material", MaterialSlab)  

                
            
                z_Offset = max_slabH + WG_Height / 2

                # Taper Widths on Bott Cal
                x = abs(WG_Height / (np.cos((angle) * np.pi / 180)))  # in Radians
                extention = np.sqrt(x ** 2 - WG_Height ** 2)
                TaperSideWidth = TaperWidth + 2 * extention
                WG_W = WG_Width + 2 * extention
                
                z_Offset = max_slabH
                
                self.lum.select("Slab")
                self.lum.addtogroup("Straight Waveguide")
                
                
                if Cladding == True:
                    # create_cover
                    self.lum.addrect()
                    self.lum.set("name", "cladding")
                    self.lum.set("material", MaterialClad)
                    self.lum.set("y", 0)
                    self.lum.set("y span", Device_Width)
                    self.lum.set("z min", max_slabH)
                    self.lum.set("z max", max_slabH*2)
                    self.lum.set("x",0)
                    self.lum.set("x span", TaperLength + WG_Width)
                    self.lum.set("override mesh order from material database", True)
                    self.lum.set("mesh order", 4)
                    self.lum.set("alpha", 0.7)
                    self.lum.select("cladding")
                    self.lum.addtogroup("Straight Waveguide")
                else:
                    pass
                
            else:
            
                z_Offset = max_subH + WG_Height / 2

                # Taper Widths on Bott Cal
                x = abs(WG_Height / (np.cos((angle) * np.pi / 180)))  # in Radians
                extention = np.sqrt(x ** 2 - WG_Height ** 2)
                TaperSideWidth = TaperWidth + 2 * extention
                WG_W = WG_Width + 2 * extention
                
                z_Offset = max_subH
                
                if Cladding == True:
                    # create_cover
                    self.lum.addrect()
                    self.lum.set("name", "cladding")
                    self.lum.set("material", MaterialClad)
                    self.lum.set("y", 0)
                    self.lum.set("y span", Device_Width)
                    self.lum.set("z min", max_subH)
                    self.lum.set("z max", max_subH*2)
                    self.lum.set("x",0)
                    self.lum.set("x span", TaperLength + WG_Width)
                    self.lum.set("override mesh order from material database", True)
                    self.lum.set("mesh order", 4)
                    self.lum.set("alpha", 0.7)
                    self.lum.select("cladding")
                    self.lum.addtogroup("Straight Waveguide")
                else:
                    pass
            
            
            
            
            
            if TaperType == "Normal":
                # PWD Taper Hights
                TaperZmin = z_Offset
                TaperZmax = z_Offset + WG_Height

                # PWB Taper Length
                TaperXmax =  TaperLength/2
                TaperXmin =  -TaperLength/2

                # Create PWB Taper
                ymin_bot_l =  - WG_W / 2
                ymax_bot_l =    WG_W / 2

                ymin_bot_r =  - TaperSideWidth / 2
                ymax_bot_r =    TaperSideWidth / 2

                ymin_top_l =  - WG_Width / 2
                ymax_top_l =    WG_Width / 2

                ymin_top_r =  - TaperWidth / 2
                ymax_top_r =    TaperWidth / 2



                vtx = np.array([[TaperXmin, ymin_bot_l, TaperZmin],  # 1
                                [TaperXmax, ymin_bot_r, TaperZmin],  # 2
                                [TaperXmax, ymax_bot_r, TaperZmin],  # 3
                                [TaperXmin, ymax_bot_l, TaperZmin],  # 4
                                [TaperXmin, ymin_top_l, TaperZmax],  # 5
                                [TaperXmax, ymin_top_r, TaperZmax],  # 6
                                [TaperXmax, ymax_top_r, TaperZmax],  # 7
                                [TaperXmin, ymax_top_l, TaperZmax],  # 8
                                ])
                a = [[np.array([[1, 4, 3, 2]], dtype=object)], [np.array([[1, 5, 8, 4]], dtype=object)],
                     [np.array([[1, 2, 6, 5]], dtype=object)], [np.array([[2, 6, 7, 3]], dtype=object)],
                     [np.array([[3, 4, 8, 7]], dtype=object)], [np.array([[5, 6, 7, 8]], dtype=object)]]



                # Send Values to Lumerical and create solid
                self.lum.putv('vertices', vtx)
                self.lum.putv('facets', a)
                self.lum.addplanarsolid(vtx, a)
                self.lum.set('material', MaterialWG)
                self.lum.set('name', "Taper")
                
                self.lum.select("Taper")
                self.lum.addtogroup("Straight Waveguide")
                
                
            
            else:
            
                # PWB Taper Length
                PWB_TaperXmax =  TaperLength/2
                PWB_TaperXmin =  -TaperLength/2
                
                        
                # PWB Taper Hights
                TaperZmin = z_Offset
                TaperZmaxF =  TaperHightF + Substrate_Height
                TaperZmaxB =  TaperHightB + Substrate_Height
                
                        
                # PWD Taper Y-Parameters
                PWB_TaperPosYMax_BotR = [(TaperWidthF / 2)]
                PWB_TaperPosYMin_BotR = [(- TaperWidthF / 2)]
                PWB_TaperPosYMax_TopR = [(TaperWidthF / 2)]
                PWB_TaperPosYMin_TopR = [(- TaperWidthF / 2)]
                PWB_TaperPosYMax_BotL = [(TaperWidthB / 2)]
                PWB_TaperPosYMin_BotL = [(- TaperWidthB / 2)]
                PWB_TaperPosYMax_TopL = [(TaperWidthB / 2)]
                PWB_TaperPosYMin_TopL = [(- TaperWidthB / 2)]
                
                # Create PWB Taper
                ymin_bot_l = PWB_TaperPosYMin_BotL[0]
                ymax_bot_l = PWB_TaperPosYMax_BotL[0]

                ymin_bot_r = PWB_TaperPosYMin_BotR[0]
                ymax_bot_r = PWB_TaperPosYMax_BotR[0]

                ymin_top_l = PWB_TaperPosYMin_TopL[0]
                ymax_top_l = PWB_TaperPosYMax_TopL[0]

                ymin_top_r = PWB_TaperPosYMin_TopR[0]
                ymax_top_r = PWB_TaperPosYMax_TopR[0]

                vtx = np.array([[PWB_TaperXmin, ymin_bot_l, TaperZmin],  # 1
                                [PWB_TaperXmax, ymin_bot_r, TaperZmin],  # 2
                                [PWB_TaperXmax, ymax_bot_r, TaperZmin],  # 3
                                [PWB_TaperXmin, ymax_bot_l, TaperZmin],  # 4
                                [PWB_TaperXmin, ymin_top_l, TaperZmaxB],  # 5
                                [PWB_TaperXmax, ymin_top_r, TaperZmaxF],  # 6
                                [PWB_TaperXmax, ymax_top_r, TaperZmaxF],  # 7
                                [PWB_TaperXmin, ymax_top_l, TaperZmaxB],  # 8
                                ])
                a = [[np.array([[1, 4, 3, 2]], dtype=object)], [np.array([[1, 5, 8, 4]], dtype=object)],
                     [np.array([[1, 2, 6, 5]], dtype=object)], [np.array([[2, 6, 7, 3]], dtype=object)],
                     [np.array([[3, 4, 8, 7]], dtype=object)], [np.array([[5, 6, 7, 8]], dtype=object)]]

                # Send Values to Lumerical and create solid
                self.lum.putv('vertices', vtx)
                self.lum.putv('facets', a)
                self.lum.addplanarsolid(vtx, a)
                self.lum.set('material', MaterialWG)
                self.lum.set("override mesh order from material database",1)
                self.lum.set("mesh order",3)
                self.lum.set('name',  "Taper")

                self.lum.select("Taper")
                self.lum.addtogroup("Straight Waveguide")
                
                self.lum.select("Straight Waveguide::cladding")
                self.lum.delete()
                
                
               

    def BendWaveguide(self, Parameters):
        '''


        Parameters
        ----------
        Parameters : Dictionary
            Dictionary with all the data needet for the Bend Wavaguide. 
            
            Parameters['Substrate Height'] : int/float
                Substrate Height
            Parameters['WG Height'] : int/float
                Waveguide Height
            Parameters['WG Width'] : int/float
                Waveguide width
            Parameters['angle'] : int/float
                Angle of the Waveguide Walls. it is calculated WG_angle = 90 - angle.
                For anfle = 90 we get a perfect rect!
            Parameters['Slab Height'] : int/float
                Slab height
            Parameters['Material'] : list
                List of Materials. The list should be with names (str) of a valid Lumerical materials.
                Check the names in Lumerical Materials viewer.
                The List of materials must contain at least 2 materials! 
                Parameters['Material'] = ['Cladding/Substrat', 'Object Material'].
                For Example: Parameters['Material'] = ["SiO2 (Glass) - Palik", 'LiNbO3_20deg_X cut'].
            Parameters["x span"] : int/float
                Length of the S Curve. Span of the object.
            Parameters["y span"] : int/float
                Height of the curve. Difference between the input and output of the S-curve.
            Parameters["poles"] : boolen
                If Parameters["poles"] = True an Bezier Curbe will be made.
                if Parameters["poles"] = False an Cosinus Curve = y_span*(cos((pi/(2*x_span))*t)^2) will be made. Where
                t is in the range of 0 - y_span
            Parameters["Cladding"] : anything, optional
                This function will check if you have set Parameters["Cladding"] to anaything, for example "Parameters["Cladding"]=1" 
                and if so it will put cladding over your structure. If the user didnt give the "Cladding" as dictionary key no cladding 
                will be set.

        Raises
        ------
        ValueError
            Value Error Massage.

        Returns
        -------
        None.

        '''

        Substrate_Height = Parameters['Substrate Height']
        WG_Height = Parameters['WG Height']
        WG_Width = Parameters['WG Width']
        angle = Parameters['angle']
        Slab_Height = Parameters['Slab Height']
        Material = Parameters['Material']
        x_span = Parameters["x span"]
        y_span = Parameters["y span"]
        polesList = Parameters["poles"]
        if "Cladding" in list(Parameters.keys()):
            Cladding = True
        else:
            Cladding = False



        if polesList == False:
            pole = np.array([[0, 0], [x_span / 2, 0], [x_span / 2, y_span], [x_span, y_span]])
        elif polesList == True:
            K = ((np.pi -2)/np.pi)*100
            pole = np.array([[0, 0], [0+(Parameters["x span"]*(K/100)), 0], [(Parameters["x span"])-(Parameters["x span"]*(K/100)), Parameters["y span"]], [Parameters["x span"], Parameters["y span"]]])
        else:
            raise ValueError('Parameters["poles"] should be an boolen variable! Please set Parameters["poles"] to True if Bezier Curves needed. Set Parameters["poles"] to False for Cosinus Curve!')


        # Material definition
        if len(Material) < 2:
            raise ValueError(
                "List of materials must contain at least 2 materials!, Parameters['Material'] = ['Cladding/Substrat', 'Object Material']")
        else:
            MaterialSub = Material[0]
            MaterialClad = Material[0]
            MaterialSlab = Material[1]
            MaterialWG = MaterialSlab





        # Device specifications
        Device_Width = y_span + 3e-6  # MMI_Width

        # creating the substrate
        max_subH = Substrate_Height
        min_subH = -Substrate_Height
        min_subL = 0
        max_subL = x_span

        self.lum.addrect()
        self.lum.set("name", "Substrate")
        self.lum.set("y", y_span / 2)
        self.lum.set("y span", Device_Width)
        self.lum.set("z min", min_subH)
        self.lum.set("z max", max_subH)
        self.lum.set("x min", min_subL)
        self.lum.set("x max", max_subL)
        self.lum.set("material", MaterialSub)
        
        
        if Slab_Height == 0:
            z_Offset = max_subH + WG_Height / 2
        else:
        
            # creating the thin film
            min_slabH = max_subH
            max_slabH = max_subH + Slab_Height

            self.lum.addrect()
            self.lum.set("name", "Slab")
            self.lum.set("y", y_span / 2)
            self.lum.set("y span", Device_Width)
            self.lum.set("z min", min_slabH)
            self.lum.set("z max", max_slabH)
            self.lum.set("x min", min_subL)
            self.lum.set("x max", max_subL)
            self.lum.set("material", MaterialSlab)

            z_Offset = max_slabH + WG_Height / 2
            self.lum.select('Slab')
            self.lum.addtogroup("S Bend")
            
            
        # Triangle EQ for waveguide Width
        x = abs(WG_Height / (np.cos((angle) * np.pi / 180)))  # in Radians
        extention = np.sqrt(x ** 2 - WG_Height ** 2)
        WG_W = WG_Width + 2 * extention

        names = ["S-Bend"]

        self.lum.addwaveguide()
        self.lum.set("name", names[0])
        self.lum.set("x", 0)
        self.lum.set("y", 0)
        self.lum.set("z", z_Offset)
        self.lum.set("base width", WG_W)
        self.lum.set("base height", WG_Height)
        self.lum.set("base angle", 90 - angle)
        self.lum.set("poles", pole)
        self.lum.set("material", MaterialWG)


        # if polesList == None:
        #     pole = np.array([[0, 0], [x_span / 2, 0], [x_span / 2, y_span], [x_span, y_span]])
        #     self.lum.set("poles", pole)
        #     self.lum.set("material", MaterialWG)
        # else:
        #     pole = np.array(polesList)
        #     self.lum.set("poles", pole)
        #     self.lum.set("material", MaterialWG)
        
        if Cladding == True:
            # create_cover
            self.lum.addrect()
            self.lum.set("name", "Cladding")
            self.lum.set("material", MaterialClad)
            self.lum.set("y", y_span / 2)
            self.lum.set("y span", Device_Width)
            self.lum.set("z", z_Offset)
            self.lum.set("z span", WG_Height*2)
            self.lum.set("x", x_span / 2)
            self.lum.set("x span", max_subL)
            self.lum.set("override mesh order from material database", True)
            self.lum.set("mesh order", 4)
            self.lum.set("alpha", 0.7)
            self.lum.select("Cladding")
            self.lum.addtogroup("S Bend")
        else:
            pass

        
        
        self.lum.select("S-Bend")
        self.lum.addtogroup("S Bend")
        





    def ArcWaveguide(self, Parameters):
        '''


        Parameters
        ----------
        Parameters : Dictionary
            Dictionary with all the data needet for the Arc Wavaguide.
            
            Parameters['Substrate Height'] : int/float
                Substrate height
            Parameters['WG Height'] : int/float
                Waveguide height
            Parameters['WG Width'] : int/float
                Waveguide width
            Parameters['angle'] : int/float
                Angle of the Waveguide Walls. it is calculated WG_angle = 90 - angle.
                For anfle = 90 we get a perfect rect!
            Parameters['Slab Height'] : int/float
                Slab height
            Parameters['Material'] : list
                List of Materials. The list should be with names (str) of a valid Lumerical materials.
                Check the names in Lumerical Materials viewer.
                The List of materials must contain at least 2 materials! 
                Parameters['Material'] = ['Cladding/Substrat', 'Object Material'].
                For Example: Parameters['Material'] = ["SiO2 (Glass) - Palik", 'LiNbO3_20deg_X cut'].
            Parameters["Wavelength"] : int/float
                Wavelength
            Parameters["S_Band Radius"] : int/float
                Radius of the Circle in um
            Parameters['Arc deg'] : int
                Can be 90 or 180 for 1/4 of a cirle or 1/2 of a circle.
            Parameters["Cladding"] : anything, optional
                This function will check if you have set Parameters["Cladding"] to anaything, for example "Parameters["Cladding"]=1" 
                and if so it will put cladding over your structure. If the user didnt give the "Cladding" as dictionary key no cladding 
                will be set.

        Raises
        ------
        ValueError
            Value Error

        Returns
        -------
        None.

        '''

        Substrate_Height = Parameters['Substrate Height']
        WG_Height = Parameters['WG Height']
        WG_Width = Parameters['WG Width']
        angle = Parameters['angle']
        Slab_Height = Parameters['Slab Height']
        Material = Parameters['Material']
        radius = Parameters["S_Band Radius"]
        arc = Parameters['Arc deg']
        
        if "Cladding" in list(Parameters.keys()):
            Cladding = True
        else:
            Cladding = False



        # Material definition
        if len(Material) < 2:
            raise ValueError(
                "List of materials must contain at least 2 materials!, Parameters['Material'] = ['Cladding/Substrat', 'Object Material']")
        else:
            MaterialSub = Material[0]
            MaterialClad = Material[0]
            MaterialSlab = Material[1]
            MaterialWG = MaterialSlab




        if arc == 90:

           

            # magic number
            # The cubic Bezier curve using this magic number in the pole points approximates the semi-circile with least error
            m = 0.55191502449

            # creating the substrate
            max_subH = Substrate_Height
            min_subH = -Substrate_Height
            # max_subW = Device_Width / 2
            # min_subW = -Device_Width / 2
            # min_subL = 0
            # max_subL = radius

            # Device specifications
            # Device_Width = 2 * radius + WaveLength * 2 + WG_Width * 2  # MMI_Width

            # creating the substrate
            max_subH = Substrate_Height
            min_subH = -Substrate_Height
            # max_subW = Device_Width / 2
            # min_subW = -Device_Width / 2
            # min_subL = 0
            # max_subL = radius

            self.lum.addrect()
            self.lum.set("name", "Substrate")
            self.lum.set("y", (m * radius))
            self.lum.set("y span", m * radius * 2 + WG_Width)
            self.lum.set("z min", min_subH)
            self.lum.set("z max", max_subH)
            self.lum.set("x", m * radius)
            self.lum.set("x span", (m * radius * 2 + WG_Width))
            self.lum.set("material", MaterialSub)

            if Slab_Height == 0:
                # creating the thin film
                z_Offset = max_subH + WG_Height / 2
            else:
                # creating the thin film
                min_slabH = max_subH
                max_slabH = max_subH + Slab_Height

                self.lum.addrect()
                self.lum.set("name", "Slab")
                self.lum.set("y", radius * m)
                self.lum.set("y span", m * radius * 2 + WG_Width)
                self.lum.set("z min", min_slabH)
                self.lum.set("z max", max_slabH)
                self.lum.set("x", m * radius)
                self.lum.set("x span", (m * radius * 2 + WG_Width))
                self.lum.set("material", MaterialSlab)

                z_Offset = max_slabH + WG_Height / 2
                self.lum.select('Slab')
                self.lum.addtogroup("90 Grad Bend")
                
            # Triangle EQ for waveguide Width
            x = abs(WG_Height / (np.cos((angle) * np.pi / 180)))  # in Radians
            extention = np.sqrt(x ** 2 - WG_Height ** 2)
            WG_W = WG_Width + 2 * extention

            names = ['Bend_Waveguide']

            self.lum.addwaveguide()
            self.lum.set("name", names[0])
            self.lum.set("x", 0)
            self.lum.set("y", 0)
            self.lum.set("z", z_Offset)
            self.lum.set("base width", WG_W)
            self.lum.set("base height", WG_Height)
            self.lum.set("base angle", 90 - angle)
            pole = np.array([[radius * 0, radius * 1], [radius * m, radius * 1], [radius * 1, radius * m],
                             [radius * 1, radius * 0]])
            self.lum.set("poles", pole)
            self.lum.set("material", MaterialWG)
            

            self.lum.addwaveguide()
            self.lum.set("name", "In_STR")
            self.lum.set("x", 0)
            self.lum.set("y", 0)
            self.lum.set("z", z_Offset)
            self.lum.set("base width", WG_W)
            self.lum.set("base height", WG_Height)
            self.lum.set("base angle", 90 - angle)
            pole = np.array([[0, radius * 1], [-1e-6, radius * 1]])
            self.lum.set("poles", pole)
            self.lum.set("material", MaterialWG)

            self.lum.addwaveguide()
            self.lum.set("name", "Out_STR")
            self.lum.set("x", 0)
            self.lum.set("y", 0)
            self.lum.set("z", z_Offset)
            self.lum.set("base width", WG_W)
            self.lum.set("base height", WG_Height)
            self.lum.set("base angle", 90 - angle)
            pole = np.array([[radius * 1, 0], [radius * 1 , -1e-6]])
            self.lum.set("poles", pole)
            self.lum.set("material", MaterialWG)

            self.lum.select('Bend_Waveguide')
            self.lum.addtogroup("90 Grad Bend")
            self.lum.select('In_STR')
            self.lum.addtogroup("90 Grad Bend")
            self.lum.select('Out_STR')
            self.lum.addtogroup("90 Grad Bend")
            
            
            if Cladding == True:
               # create_cover
                self.lum.addrect()
                self.lum.set("name", "cladding")
                self.lum.set("material", MaterialClad)
                self.lum.set("y", radius * m)
                self.lum.set("y span", m * radius * 2 + WG_Width)
                self.lum.set("z min", max_slabH)
                self.lum.set("z max", max_slabH*2)
                self.lum.set("x", m * radius)
                self.lum.set("x span", (m * radius * 2 + WG_Width))
                self.lum.set("override mesh order from material database", True)
                self.lum.set("mesh order", 4)
                self.lum.set("alpha", 0.7)
                self.lum.select('cladding')
                self.lum.addtogroup("180 Grad Bend")
            else:
                pass



        elif arc == 180:

            # Device specifications
            # Device_Width = 2 * radius + WaveLength * 2 + WG_Width * 2  # MMI_Width

            # magic number
            # The cubic Bezier curve using this magic number in the pole points approximates the semi-circile with least error
            m = 0.55191502449

            # creating the substrate
            max_subH = Substrate_Height
            min_subH = -Substrate_Height


            self.lum.addrect()
            self.lum.set("name", "Substrate")
            self.lum.set("y", (m * radius))
            self.lum.set("y span", m * radius * 2 + WG_Width)
            self.lum.set("z min", min_subH)
            self.lum.set("z max", max_subH)
            self.lum.set("x", 0)
            self.lum.set("x span", (m * radius * 2 + WG_Width) * 2)
            self.lum.set("material", MaterialSub)
            
            
            if Slab_Height == 0:
                z_Offset = max_subH + WG_Height / 2
            else:
                # creating the thin film
                min_slabH = max_subH
                max_slabH = max_subH + Slab_Height

                self.lum.addrect()
                self.lum.set("name", "Slab")
                self.lum.set("y", radius * m)
                self.lum.set("y span", m * radius * 2 + WG_Width)
                self.lum.set("z min", min_slabH)
                self.lum.set("z max", max_slabH)
                self.lum.set("x", 0)
                self.lum.set("x span", (m * radius * 2 + WG_Width) * 2)
                self.lum.set("material", MaterialSlab)

                z_Offset = max_slabH + WG_Height / 2
                
                self.lum.select('Slab')
                self.lum.addtogroup("180 Grad Bend")
                
            # Triangle EQ for waveguide Width
            x = abs(WG_Height / (np.cos((angle) * np.pi / 180)))  # in Radians
            extention = np.sqrt(x ** 2 - WG_Height ** 2)
            WG_W = WG_Width + 2 * extention

            names = ['Arc Waveguide1', 'Arc Waveguie2']

            self.lum.addwaveguide()
            self.lum.set("name", names[0])
            self.lum.set("x", 0)
            self.lum.set("y", 0)
            self.lum.set("z", z_Offset)
            self.lum.set("base width", WG_W)
            self.lum.set("base height", WG_Height)
            self.lum.set("base angle", 90 - angle)
            pole = np.array(
                [[radius * 0, radius * 1], [radius * m, radius * 1], [radius * 1, radius * m],
                 [radius * 1, radius * 0]])
            self.lum.set("poles", pole)
            self.lum.set("material", MaterialWG)
            self.lum.set("first axis", 'z')
            self.lum.set("rotation 1", 0)

            self.lum.addwaveguide()
            self.lum.set("name", names[1])
            self.lum.set("x", 0)
            self.lum.set("y", 0)
            self.lum.set("z", z_Offset)
            self.lum.set("base width", WG_W)
            self.lum.set("base height", WG_Height)
            self.lum.set("base angle", 90 - angle)
            pole = np.array(
                [[radius * 0, radius * 1], [radius * m, radius * 1], [radius * 1, radius * m],
                 [radius * 1, radius * 0]])
            self.lum.set("poles", pole)
            self.lum.set("material", MaterialWG)
            self.lum.set("first axis", 'z')
            self.lum.set("rotation 1", 90)


            
            
            self.lum.select('Arc Waveguide1')
            self.lum.addtogroup("180 Grad Bend")
            self.lum.select('Arc Waveguie2')
            self.lum.addtogroup("180 Grad Bend")
            
            if Cladding == True:
               # create_cover
                self.lum.addrect()
                self.lum.set("name", "cladding")
                self.lum.set("material", MaterialClad)
                self.lum.set("y", radius * m)
                self.lum.set("y span", m * radius * 2 + WG_Width)
                self.lum.set("z min", max_slabH)
                self.lum.set("z max", max_slabH*2)
                self.lum.set("x", 0)
                self.lum.set("x span", (m * radius * 2 + WG_Width)*2)
                self.lum.set("override mesh order from material database", True)
                self.lum.set("mesh order", 4)
                self.lum.set("alpha", 0.7)
                self.lum.select('cladding')
                self.lum.addtogroup("180 Grad Bend")
            else:
                pass


        else:
            raise ValueError(
                "Incorrect Arc Value. The arc can only be an integer and can only be arc = 90 or arc = 180!")





    def MMI2x2(self, Parameters):
        '''


        Parameters
        ----------
        Parameters : Dictionary
            Dictionary with all the data needet for the 2x2 MMI.
            
            Parameters['Material'] : list of str
                List of Materials. The list should be with names (str) of a valid Lumerical materials.
                Check the names in Lumerical Materials viewer.
                The List of materials must contain at least 2 materials! 
                Parameters['Material'] = ['Cladding/Substrat', 'Object Material'].
                For Example: Parameters['Material'] = ["SiO2 (Glass) - Palik", 'LiNbO3_20deg_X cut'].
            Parameters['Substrate Height'] : int/float
                Substrate height.
            Parameters['MMI Width'] : int/float
                Width of the MMI.
            Parameters['MMI Length'] : int/float
                Length of the MMI.
            Parameters['angle'] : int/float
                Angle of the Waveguide Walls. it is calculated WG_angle = 90 - angle.
                For anfle = 90 we get a perfect rect!
            Parameters['WG Height'] : int/float
                Waveguide hight. Also the height of the MMI section
            Parameters['WG Width'] : int/float
                Waveguide width.
            Parameters['WG Length'] : int/float
                Waveguide length.
            Parameters['Position Offset'] : int/float
                Offset between the waveguides. If Taper == True then this become the offset
                betweent he tapers wider sides. Waveguide and Tapers cannot be placed ourside
                the MMI structure. The minimum distance between Taper and Waveguide is 1 um
                becouse of manufactering restrictions in the University.
            Parameters['Slab Height'] : int/float
                Height of the slab.
            Parameters['Taper'] : boolen
                Taper can be set to be True ot False.
                if Taper == False - No Taper used
                if Taper == True - Taper placed
            Parameters['Taper Length'] : int/float
                If Taper == True, then this will set the Tapers length. If Taper == False
                this will be ignored and some random value can be given.
            Parameters['Taper Width'] : int/float
                If Taper == True, then this will set the Tapers width. If Taper == False
                this will be ignored and some random value can be given.
            Parameters["Cladding"] : anything, optional
                    This function will check if you have set Parameters["Cladding"] to anaything, for example "Parameters["Cladding"]=1" 
                    and if so it will put cladding over your structure. If the user didnt give the "Cladding" as dictionary key no cladding 
                    will be set.

        Raises
        ------
        ValueError
            Value Error.

        Returns
        -------
        None.

        '''

        Substrate_Height = Parameters['Substrate Height']
        Material = Parameters['Material']
        MMI_Width = Parameters['MMI Width']
        MMI_Length = Parameters['MMI Length']
        angle = Parameters['angle']
        WG_Height = Parameters['WG Height']
        WG_Width = Parameters['WG Width']
        WG_Length = Parameters['WG Length']
        posOffset = Parameters['Position Offset']
        Slab_Height = Parameters['Slab Height']
        TaperLength = Parameters['Taper Length']
        TaperWidth = Parameters['Taper Width']
        Taper = Parameters['Taper']
        
        if "Cladding" in list(Parameters.keys()):
            Cladding = True
        else:
            Cladding = False



        # Material definition
        if len(Material) < 2:
            raise ValueError(
                "List of materials must contain at least 2 materials!, Parameters['Material'] = ['Cladding/Substrat', 'Object Material']")
        else:
            MaterialSub = Material[0]
            MaterialClad = Material[0]
            MaterialSlab = Material[1]
            MaterialWG = MaterialSlab


        # Device specifications
        Device_Length = MMI_Length + 2 * WG_Length
        Device_Width = MMI_Width + 3e-6  # MMI_Width

        # creating the substrate
        max_subH = Substrate_Height
        min_subH = -Substrate_Height
        max_subW = Device_Width / 2
        min_subW = -Device_Width / 2
        min_subL = -Device_Length / 2
        max_subL = Device_Length / 2


        self.lum.addrect()
        self.lum.set("name", "Substrate")
        self.lum.set("y min", min_subW)
        self.lum.set("y max", max_subW)
        self.lum.set("z min", min_subH)
        self.lum.set("z max", max_subH)
        self.lum.set("x min", min_subL)
        self.lum.set("x max", max_subL)
        self.lum.set("material", MaterialSub)
        
        # creating the MMI
        max_MMIH = WG_Height
        max_MMIL = MMI_Length / 2
        min_MMIL = -MMI_Length / 2

        if Slab_Height == 0:
            z_Offset = max_subH + max_MMIH / 2
        else:      
            # creating the thin film
            min_slabH = max_subH
            max_slabH = max_subH + Slab_Height

            self.lum.addrect()
            self.lum.set("name", "Slab")
            self.lum.set("y min", min_subW)
            self.lum.set("y max", max_subW)
            self.lum.set("z min", min_slabH)
            self.lum.set("z max", max_slabH)
            self.lum.set("x min", min_subL)
            self.lum.set("x max", max_subL)
            self.lum.set("material", MaterialSlab)

           
            z_Offset = max_slabH + max_MMIH / 2
            self.lum.select('Slab')
            self.lum.addtogroup("MMI Object")
            
            
        # Triangle EQ for MMI Width
        x = abs(max_MMIH / (np.cos((angle) * np.pi / 180)))  # in Radians
        extention = np.sqrt(x ** 2 - max_MMIH ** 2)
        MMI_Wid = MMI_Width + 2 * extention

        self.lum.addwaveguide()
        self.lum.set("name", "MMI")
        self.lum.set("base height", max_MMIH)
        self.lum.set("base angle", 90 - angle)
        self.lum.set("x", 0)
        self.lum.set("y", 0)
        self.lum.set("z", z_Offset)
        pole = np.array([[max_MMIL, 0], [min_MMIL, 0]])
        self.lum.set("poles", pole)
        self.lum.set("material", MaterialWG)
        self.lum.set("base width", MMI_Wid)

        # Positions of the Input and Output WGs
        # Triangle EQ for MMI Width
        x = abs(max_MMIH / (np.cos((angle) * np.pi / 180)))  # in Radians
        extention = np.sqrt(x ** 2 - max_MMIH ** 2)
        WG_W = WG_Width + 2 * extention
        WG_Width_top = WG_W
        OffMax = MMI_Width / 2

        if Taper == False:

            #Too Fara and Too close
            offset_WG = posOffset / 2 + WG_Width / 2 + WG_W / 2
            offset_WG2 = posOffset / 2



            if offset_WG > MMI_Wid / 2:
                self.lum.deleteall()
                raise ValueError('You are Trying to move the Waveguide outside the MMI. This is not possible!')

            # elif offset_WG2 <0.5e-6:
            #     self.lum.deleteall()
            #     raise ValueError('The distance between the Tapers is less then 1 um !')



            else:
                # Mirror the In and Out WG on both sides
                maxWGL = [WG_Length, WG_Length, 0, 0]
                minWGL = [0, 0, -WG_Length, -WG_Length]
                xPos = [max_MMIL, max_MMIL, min_MMIL, min_MMIL]
                yPos = [WG_Width / 2 + posOffset / 2, -(WG_Width / 2 + posOffset / 2), WG_Width / 2 + posOffset / 2,
                        -(WG_Width / 2 + posOffset / 2)]

                # Names of the WGs
                names = ['Input WG_L', 'Input WG_R', 'Output WG_L', 'Output WG_R']

                # create loop
                for i in range(len(xPos)):
                    self.lum.addwaveguide()
                    self.lum.set("name", names[i])
                    self.lum.set("x", xPos[i])
                    self.lum.set("y", yPos[i])
                    self.lum.set("z", z_Offset)
                    self.lum.set("base width", WG_Width_top)
                    self.lum.set("base height", max_MMIH)
                    self.lum.set("base angle", 90 - angle)
                    pole = np.array([[maxWGL[i], 0], [minWGL[i], 0]])
                    self.lum.set("poles", pole)
                    self.lum.set("material", MaterialSlab)
                    
            self.lum.select('MMI')
            self.lum.addtogroup("MMI Object")
            self.lum.select('Input WG_L')
            self.lum.addtogroup("MMI Object")
            self.lum.select('Input WG_R')
            self.lum.addtogroup("MMI Object")
            self.lum.select('Output WG_L')
            self.lum.addtogroup("MMI Object")
            self.lum.select('Output WG_R')
            self.lum.addtogroup("MMI Object")
                    
            if Cladding == True:
               # create_cover
                self.lum.addrect()
                self.lum.set("name", "cladding")
                self.lum.set("material", MaterialClad)
                self.lum.set("y min", min_subW)
                self.lum.set("y max", max_subW )
                self.lum.set("z", z_Offset)
                self.lum.set("z span", max_MMIH * 2)
                self.lum.set("x min", min_subL)
                self.lum.set("x max", max_subL)
                self.lum.set("override mesh order from material database", True)
                self.lum.set("mesh order", 4)
                self.lum.set("alpha", 0.7)
                self.lum.select('cladding')
                self.lum.addtogroup("MMI Object")
            else:
                pass

           
            


        elif Taper == True:

            # Delate the Structure to start new
            self.lum.deleteall()
            # Device specifications
            Device_Length = MMI_Length + 2 * WG_Length + 2 * TaperLength
            Device_Width = MMI_Width + 3e-6  # MMI_Width

            # creating the substrate
            max_subH = Substrate_Height
            min_subH = -Substrate_Height
            max_subW = Device_Width / 2
            min_subW = -Device_Width / 2
            min_subL = -Device_Length / 2
            max_subL = Device_Length / 2

            self.lum.addrect()
            self.lum.set("name", "Substrate")
            self.lum.set("y min", min_subW)
            self.lum.set("y max", max_subW)
            self.lum.set("z min", min_subH)
            self.lum.set("z max", max_subH)
            self.lum.set("x min", min_subL)
            self.lum.set("x max", max_subL)
            self.lum.set("material", MaterialSub)

            # creating the MMI
            max_MMIH = WG_Height
            max_MMIL = MMI_Length / 2
            min_MMIL = -MMI_Length / 2
            
            if Slab_Height == 0:
                z_Offset = max_subH + max_MMIH / 2
                
            else:
                # creating the thin film
                min_slabH = max_subH
                max_slabH = max_subH + Slab_Height

                self.lum.addrect()
                self.lum.set("name", "Slab")
                self.lum.set("y min", min_subW)
                self.lum.set("y max", max_subW)
                self.lum.set("z min", min_slabH)
                self.lum.set("z max", max_slabH)
                self.lum.set("x min", min_subL)
                self.lum.set("x max", max_subL)
                self.lum.set("material", MaterialSlab)

                
                z_Offset = max_slabH + max_MMIH / 2
                
                self.lum.select('Slab')
                self.lum.addtogroup("MMI Object")
            

            # Triangle EQ for MMI Width
            x = abs(max_MMIH / (np.cos((angle) * np.pi / 180)))  # in Radians
            extention = np.sqrt(x ** 2 - max_MMIH ** 2)
            MMI_Wid = MMI_Width + 2 * extention

            self.lum.addwaveguide()
            self.lum.set("name", "MMI")
            self.lum.set("x", 0)
            self.lum.set("y", 0)
            self.lum.set("z", z_Offset)
            self.lum.set("base height", max_MMIH)
            self.lum.set("base angle", 90 - angle)
            pole = np.array([[max_MMIL, 0], [min_MMIL, 0]])
            self.lum.set("poles", pole)
            self.lum.set("material", MaterialWG)
            self.lum.set("base width", MMI_Wid)

            # New x Length of the Tapers
            maxLength = max_MMIL + TaperLength
            minLength = min_MMIL - TaperLength

            # Mirror the In and Out WG on both sides
            maxWGL = [WG_Length, WG_Length, 0, 0]
            minWGL = [0, 0, -WG_Length, -WG_Length]
            xPos = [maxLength, maxLength, minLength, minLength]
            yPos = [WG_Width / 2 + posOffset / 2, -(WG_Width / 2 + posOffset / 2), WG_Width / 2 + posOffset / 2,
                    -(WG_Width / 2 + posOffset / 2)]

            # Names of the WGs
            names = ['Input WG_L', 'Input WG_R', 'Output WG_L', 'Output WG_R']
            TapersNames = ['Taper Input WG_L', 'Taper Input WG_R', 'Taper Output WG_L', 'Taper Output WG_R']

            # Taper Widths on Bott Cal
            x = abs(max_MMIH / (np.cos((angle) * np.pi / 180)))  # in Radians
            extention = np.sqrt(x ** 2 - max_MMIH ** 2)
            TaperSideWidth = TaperWidth + 2 * extention

            TaperPosXmin = [max_MMIL, max_MMIL, min_MMIL, min_MMIL]
            TaperPosXmax = [max_MMIL + TaperLength, max_MMIL + TaperLength, min_MMIL - TaperLength,
                            min_MMIL - TaperLength]

            PosOffset = [(posOffset / 2 + WG_Width / 2), - (posOffset / 2 + WG_Width / 2),
                         (posOffset / 2 + WG_Width / 2), -(posOffset / 2 + WG_Width / 2)]
            TaperPosYMax_BotR = [(WG_W / 2), (WG_W / 2), (- WG_W / 2), (- WG_W / 2)]
            TaperPosYMin_BotR = [(- WG_W / 2), (-WG_W / 2), (WG_W / 2), (WG_W / 2)]
            TaperPosYMax_TopR = [(WG_Width / 2), (WG_Width / 2), (-WG_Width / 2), (- WG_Width / 2)]
            TaperPosYMin_TopR = [(- WG_Width / 2), (-WG_Width / 2), (WG_Width / 2), (WG_Width / 2)]

            TaperPosYMax_BotL = [(TaperSideWidth / 2), (TaperSideWidth / 2), (- TaperSideWidth / 2),
                                 (- TaperSideWidth / 2)]
            TaperPosYMin_BotL = [(- TaperSideWidth / 2), (- TaperSideWidth / 2), (TaperSideWidth / 2),
                                 (TaperSideWidth / 2)]
            TaperPosYMax_TopL = [(TaperWidth / 2), (TaperWidth / 2), (- TaperWidth / 2), (- TaperWidth / 2)]
            TaperPosYMin_TopL = [(- TaperWidth / 2), (- TaperWidth / 2), (TaperWidth / 2), (TaperWidth / 2)]

            offset_Taper = posOffset / 2 + WG_Width / 2 + TaperWidth / 2  # + WG_W / 2
            BotCornerDistance = posOffset - TaperWidth / 2
            # offset_Set_R = posOffset/2+TaperWidth/2

            # if OffsetInpit <= OffMin or OffsetInpit >= OffMax or posOffset <= OffMin or posOffset >= OffMax:
            if offset_Taper > OffMax:
                self.lum.deleteall()
                raise ValueError('You are Trying to move the Taper outside the MMI. This is not possible!')
            # elif BotCornerDistance < 1e-6:
            #     self.lum.deleteall()
            #     raise ValueError('The distance between the Tapers is less then 1 um !')
            



            else:
                # Mirror the In and Out WG on both sides
                maxWGL = [WG_Length, WG_Length, 0, 0]
                minWGL = [0, 0, -WG_Length, -WG_Length]
                # xPos = [max_MMIL, max_MMIL, min_MMIL, min_MMIL]
                xPos = [maxLength, maxLength, minLength, minLength]

                for i in range(len(xPos)):
                    TaperZmin = max_slabH
                    TaperZmax = max_slabH + max_MMIH

                    TaperXmin = TaperPosXmin[i]
                    TaperXmax = TaperPosXmax[i]

                    ymin_bot_l = TaperPosYMin_BotL[i]
                    ymax_bot_l = TaperPosYMax_BotL[i]

                    ymin_bot_r = TaperPosYMin_BotR[i]
                    ymax_bot_r = TaperPosYMax_BotR[i]

                    ymin_top_l = TaperPosYMin_TopL[i]
                    ymax_top_l = TaperPosYMax_TopL[i]

                    ymin_top_r = TaperPosYMin_TopR[i]
                    ymax_top_r = TaperPosYMax_TopR[i]

                    vtx = np.array([[TaperXmin, ymin_bot_l, TaperZmin],  # 1
                                    [TaperXmax, ymin_bot_r, TaperZmin],  # 2
                                    [TaperXmax, ymax_bot_r, TaperZmin],  # 3
                                    [TaperXmin, ymax_bot_l, TaperZmin],  # 4
                                    [TaperXmin, ymin_top_l, TaperZmax],  # 5
                                    [TaperXmax, ymin_top_r, TaperZmax],  # 6
                                    [TaperXmax, ymax_top_r, TaperZmax],  # 7
                                    [TaperXmin, ymax_top_l, TaperZmax],  # 8
                                    ])
                    a = [[np.array([[1, 4, 3, 2]], dtype=object)], [np.array([[1, 5, 8, 4]], dtype=object)],
                         [np.array([[1, 2, 6, 5]], dtype=object)], [np.array([[2, 6, 7, 3]], dtype=object)],
                         [np.array([[3, 4, 8, 7]], dtype=object)], [np.array([[5, 6, 7, 8]], dtype=object)]]

                    # Send Values to Lumerical and create solid
                    self.lum.putv('vertices', vtx)
                    self.lum.putv('facets', a)
                    self.lum.addplanarsolid(vtx, a)
                    self.lum.set('material', MaterialWG)
                    self.lum.set('name', TapersNames[i])
                    self.lum.set('y', PosOffset[i])

                # create loop
                for i in range(len(xPos)):
                    self.lum.addwaveguide()
                    self.lum.set("name", names[i])
                    self.lum.set("x", xPos[i])
                    self.lum.set("y", yPos[i])
                    self.lum.set("z", z_Offset)
                    self.lum.set("base width", WG_Width_top)
                    self.lum.set("base height", max_MMIH)
                    self.lum.set("base angle", 90 - angle)
                    pole = np.array([[maxWGL[i], 0], [minWGL[i], 0]])
                    self.lum.set("poles", pole)
                    self.lum.set("material", MaterialSlab)
                

            self.lum.select('MMI')
            self.lum.addtogroup("MMI Object")
            self.lum.select('Input WG_L')
            self.lum.addtogroup("MMI Object")
            self.lum.select('Input WG_R')
            self.lum.addtogroup("MMI Object")
            self.lum.select('Output WG_L')
            self.lum.addtogroup("MMI Object")
            self.lum.select('Output WG_R')
            self.lum.addtogroup("MMI Object")
            self.lum.select('Taper Input WG_L')
            self.lum.addtogroup("MMI Object")
            self.lum.select('Taper Input WG_R')
            self.lum.addtogroup("MMI Object")
            self.lum.select('Taper Output WG_L')
            self.lum.addtogroup("MMI Object")
            self.lum.select('Taper Output WG_R')
            self.lum.addtogroup("MMI Object")
            
            if Cladding == True:
                # create_cover
                self.lum.addrect()
                self.lum.set("name", "cladding")
                self.lum.set("material", MaterialClad)
                self.lum.set("y min", min_subW)
                self.lum.set("y max", max_subW)
                self.lum.set("z", z_Offset)
                self.lum.set("z span", max_MMIH * 2)
                self.lum.set("x min", min_subL)
                self.lum.set("x max", max_subL)
                self.lum.set("override mesh order from material database", True)
                self.lum.set("mesh order", 4)
                self.lum.set("alpha", 0.7)
                
                self.lum.select('cladding')
                self.lum.addtogroup("MMI Object")
            else:
                pass


            
        else:
            raise ValueError(
                "Incorect Taper input. Taper must be an boolen. You can choose from Taper = True or Taper = False!")







    def MMI2x1(self, Parameters):
        '''

        Parameters
        ----------
        Parameters : Dictionary
            Dictionary with all the data needet for the 2x1 MMI. 
            
                Parameters['Material'] : list of str
                    List of Materials. The list should be with names (str) of a valid Lumerical materials.
                    Check the names in Lumerical Materials viewer.
                    The List of materials must contain at least 2 materials! 
                    Parameters['Material'] = ['Cladding/Substrat', 'Object Material'].
                    For Example: Parameters['Material'] = ["SiO2 (Glass) - Palik", 'LiNbO3_20deg_X cut'].
                Parameters['Substrate Height'] : int/float
                    Substrate height.
                Parameters['MMI Width'] : int/float
                    Width of the MMI.
                Parameters['MMI Length'] : int/float
                    Length of the MMI.
                Parameters['angle'] : int/float
                    Angle of the Waveguide Walls. it is calculated WG_angle = 90 - angle.
                    For anfle = 90 we get a perfect rect!
                Parameters['WG Height'] : int/float
                    Waveguide hight. Also the height of the MMI section
                Parameters['WG Width'] : int/float
                    Waveguide width.
                Parameters['WG Length'] : int/float
                    Waveguide length.
                Parameters['Position Offset'] : int/float
                    Offset between the waveguides. If Taper == True then this become the offset
                    betweent he tapers wider sides. Waveguide and Tapers cannot be placed ourside
                    the MMI structure. The minimum distance between Taper and Waveguide is 1 um
                    becouse of manufactering restrictions in the University.
                Parameters['Offset Input'] : int/float
                    Offset of the input waveguide.
                Parameters['Slab Height'] : int/float
                    Height of the slab.
                Parameters['Taper'] : boolen
                    Taper can be set to be True ot False.
                    if Taper == False - No Taper used
                    if Taper == True - Taper placed
                Parameters['Taper Length'] : int/float
                    If Taper == True, then this will set the Tapers length. If Taper == False
                    this will be ignored and some random value can be given.
                Parameters['Taper Width'] : int/float
                    If Taper == True, then this will set the Tapers width. If Taper == False
                    this will be ignored and some random value can be given.
                Parameters["Cladding"] : anything, optional
                    This function will check if you have set Parameters["Cladding"] to anaything, for example "Parameters["Cladding"]=1" 
                    and if so it will put cladding over your structure. If the user didnt give the "Cladding" as dictionary key no cladding 
                    will be set.
                Parameters["Offset Output"] : anything, optional
                    This function will allow the user to move the outputs in oposite direction. Please dont use it since is there only 
                    becouse the maschine of our physic departmant had some proiblems with the LNOI objects design. 

        Raises
        ------
        ValueError
            DESCRIPTION.

        Returns
        -------
        None.

        '''

        Material = Parameters['Material']
        Substrate_Height = Parameters['Substrate Height']
        MMI_Width = Parameters['MMI Width']
        MMI_Length = Parameters['MMI Length']
        angle = Parameters['angle']
        WG_Height = Parameters['WG Height']
        WG_Width = Parameters['WG Width']
        WG_Length = Parameters['WG Length']
        OffsetInput = Parameters['Offset Input']
        posOffset = Parameters['Position Offset']
        Slab_Height = Parameters['Slab Height']
        TaperLength = Parameters['Taper Length']
        TaperWidth = Parameters['Taper Width']
        Taper = Parameters['Taper']
        
        if "Cladding" in list(Parameters.keys()):
            Cladding = True
        else:
            Cladding = False
        
        if 'Offset Output' not in list(Parameters.keys()):
            OffsetOutput = None
        else:
            OffsetOutput = Parameters['Offset Output']




        # Material definition
        if len(Material) < 2:
            raise ValueError(
                "List of materials must contain at least 2 materials!, Parameters['Material'] = ['Cladding/Substrat', 'Object Material']")
        else:
            MaterialSub = Material[0]
            MaterialClad = Material[0]
            MaterialSlab = Material[1]
            MaterialWG = MaterialSlab



        # Device specifications
        Device_Length = MMI_Length + 2 * WG_Length
        Device_Width = MMI_Width + 3e-6  # MMI_Width

        # creating the substrate
        max_subH = Substrate_Height
        min_subH = -Substrate_Height
        max_subW = Device_Width / 2
        min_subW = -Device_Width / 2
        min_subL = -Device_Length / 2
        max_subL = Device_Length / 2

        self.lum.addrect()
        self.lum.set("name", "Substrate")
        self.lum.set("y min", min_subW)
        self.lum.set("y max", max_subW)
        self.lum.set("z min", min_subH)
        self.lum.set("z max", max_subH)
        self.lum.set("x min", min_subL)
        self.lum.set("x max", max_subL)
        self.lum.set("material", MaterialSub)
        
        # creating the MMI
        max_MMIH = WG_Height
        max_MMIL = MMI_Length / 2
        min_MMIL = -MMI_Length / 2
        
        if Slab_Height == 0:
            z_Offset = max_subH + max_MMIH / 2
            
        else:
            # creating the thin film
            min_slabH = max_subH
            max_slabH = max_subH + Slab_Height

            self.lum.addrect()
            self.lum.set("name", "Slab")
            self.lum.set("y min", min_subW)
            self.lum.set("y max", max_subW)
            self.lum.set("z min", min_slabH)
            self.lum.set("z max", max_slabH)
            self.lum.set("x min", min_subL)
            self.lum.set("x max", max_subL)
            self.lum.set("material", MaterialSlab)

            
            z_Offset = max_slabH + max_MMIH / 2
            
            self.lum.select("Slab")
            self.lum.addtogroup("MMI Object")
            
            

        # Triangle EQ for MMI Width
        x = abs(max_MMIH / (np.cos((angle) * np.pi / 180)))  # in Radians
        extention = np.sqrt(x ** 2 - max_MMIH ** 2)
        MMI_Wid = MMI_Width + 2 * extention

        self.lum.addwaveguide()
        self.lum.set("name", "MMI")
        self.lum.set("x", 0)
        self.lum.set("y", 0)
        self.lum.set("z", z_Offset)
        self.lum.set("base height", max_MMIH)
        self.lum.set("base angle", 90 - angle)
        pole = np.array([[max_MMIL, 0], [min_MMIL, 0]])
        self.lum.set("poles", pole)
        self.lum.set("material", MaterialWG)
        self.lum.set("base width", MMI_Wid)

        # Positions of the Input and Output WGs
        # Triangle EQ for MMI Width
        x = abs(max_MMIH / (np.cos((angle) * np.pi / 180)))  # in Radians
        extention = np.sqrt(x ** 2 - max_MMIH ** 2)
        WG_W = WG_Width + 2 * extention
        WG_Width_top = WG_W
        OffMax = MMI_Width / 2

        offset_Taper = posOffset / 2 + WG_Width / 2 + TaperWidth / 2  # + WG_W / 2
        BotCornerDistance = posOffset/2 - TaperWidth / 2
        offset_WG = posOffset / 2 + WG_Width / 2 + WG_W / 2
        offset_WG2 = posOffset / 2

        # if offset_WG2 < 0.5e-6:
            # self.lum.deleteall()
            # raise ValueError('The distance between the Tapers is less then 1 um !')
        # else:

        if Taper == False:

            # if offset_WG > OffMax:
                # self.lum.deleteall()
                # raise ValueError('You are Trying to move the Waveguide outside the MMI. This is not possible!')
            # else:
            # Mirror the In and Out WG on both sides
            maxWGL = [WG_Length, 0, 0]
            minWGL = [0, -WG_Length, -WG_Length]
            xPos = [max_MMIL, min_MMIL, min_MMIL]
            if OffsetOutput == None:
                yPos = [0 + OffsetInput, (WG_Width / 2 + posOffset / 2) , (- WG_Width / 2 - posOffset / 2)  ]
            else:
                yPos = [0 + OffsetInput, (WG_Width / 2 + posOffset / 2) + OffsetOutput, (- WG_Width / 2 - posOffset / 2) + OffsetOutput ]

            # Names of the WGs
            names = ['Input WG', 'Output WG_L', 'Output WG_R']

            # create loop
            for i in range(len(xPos)):
                self.lum.addwaveguide()
                self.lum.set("name", names[i])
                self.lum.set("x", xPos[i])
                self.lum.set("y", yPos[i])
                self.lum.set("z", z_Offset)
                self.lum.set("base width", WG_Width_top)
                self.lum.set("base height", max_MMIH)
                self.lum.set("base angle", 90 - angle)
                pole = np.array([[maxWGL[i], 0], [minWGL[i], 0]])
                self.lum.set("poles", pole)
                self.lum.set("material", MaterialSlab)


            self.lum.select("MMI")
            self.lum.addtogroup("MMI Object")
            self.lum.select("Input WG")
            self.lum.addtogroup("MMI Object")
            self.lum.select("Output WG_L")
            self.lum.addtogroup("MMI Object")
            self.lum.select("Output WG_R")
            self.lum.addtogroup("MMI Object")
            
            if Cladding == True:
                # create_cover
                self.lum.addrect()
                self.lum.set("name", "cladding")
                self.lum.set("material", MaterialClad)
                self.lum.set("y min", min_subW)
                self.lum.set("y max", max_subW)
                self.lum.set("z", z_Offset)
                self.lum.set("z span", max_MMIH + 0.4e-6)
                self.lum.set("x min", min_subL)
                self.lum.set("x max", max_subL)
                self.lum.set("override mesh order from material database", True)
                self.lum.set("mesh order", 4)
                self.lum.set("alpha", 0.7)
                
                self.lum.select("cladding")
                self.lum.addtogroup("MMI Object")
            else:
                pass
            
            

        elif Taper == True:
            if offset_Taper > OffMax:
                self.lum.deleteall()
                raise ValueError('You are Trying to move the Taper outside the MMI. This is not possible!')
            else:
                # Delate the Structure to start new
                self.lum.deleteall()
                # Device specifications
                Device_Length = MMI_Length + 2 * WG_Length + 2*TaperLength
                Device_Width = MMI_Width + 3e-6 # WaveLength * 2  

                # creating the substrate
                max_subH = Substrate_Height
                min_subH = -Substrate_Height
                max_subW = Device_Width / 2
                min_subW = -Device_Width / 2
                min_subL = -Device_Length / 2
                max_subL = Device_Length / 2

                self.lum.addrect()
                self.lum.set("name", "Substrate")
                self.lum.set("y min", min_subW)
                self.lum.set("y max", max_subW)
                self.lum.set("z min", min_subH)
                self.lum.set("z max", max_subH)
                self.lum.set("x min", min_subL)
                self.lum.set("x max", max_subL)
                self.lum.set("material", MaterialSub)
                
                # creating the MMI
                max_MMIH = WG_Height
                max_MMIL = MMI_Length / 2
                min_MMIL = -MMI_Length / 2
                
                if Slab_Height == 0:
                    z_Offset = max_subH + max_MMIH / 2
                else:
                    

                    # creating the thin film
                    min_slabH = max_subH
                    max_slabH = max_subH + Slab_Height

                    self.lum.addrect()
                    self.lum.set("name", "Slab")
                    self.lum.set("y min", min_subW)
                    self.lum.set("y max", max_subW)
                    self.lum.set("z min", min_slabH)
                    self.lum.set("z max", max_slabH)
                    self.lum.set("x min", min_subL)
                    self.lum.set("x max", max_subL)
                    self.lum.set("material", MaterialSlab)

                    
                    z_Offset = max_slabH + max_MMIH / 2
                    self.lum.select("Slab")
                    self.lum.addtogroup("MMI Object")
                    
                    

                # Triangle EQ for MMI Width
                x = abs(max_MMIH / (np.cos((angle) * np.pi / 180)))  # in Radians
                extention = np.sqrt(x ** 2 - max_MMIH ** 2)
                MMI_Wid = MMI_Width + 2 * extention

                self.lum.addwaveguide()
                self.lum.set("name", "MMI")
                self.lum.set("x", 0)
                self.lum.set("y", 0)
                self.lum.set("z", z_Offset)
                self.lum.set("base height", max_MMIH)
                self.lum.set("base angle", 90 - angle)
                pole = np.array([[max_MMIL, 0], [min_MMIL, 0]])
                self.lum.set("poles", pole)
                self.lum.set("material", MaterialWG)
                self.lum.set("base width", MMI_Wid)


                # New x Length of the Tapers
                maxLength = max_MMIL + TaperLength
                minLength = min_MMIL - TaperLength

                # Mirror the In and Out WG on both sides
                maxWGL = [WG_Length, 0, 0]
                minWGL = [0, -WG_Length, -WG_Length]
                xPos = [maxLength, minLength, minLength]
                yPos = [0 + OffsetInput, WG_Width / 2 + posOffset / 2, - WG_Width / 2 - posOffset / 2]

                # Names of the WGs
                names = ['Input WG', 'Output WG_L', 'Output WG_R']
                TapersNames = ['Taper Input WG', 'Taper Output WG_L', 'Taper Output WG_R']

                # Taper loop
                # Taper Widths on Bott Cal
                x = abs(max_MMIH / (np.cos((angle) * np.pi / 180)))  # in Radians
                extention = np.sqrt(x ** 2 - max_MMIH ** 2)
                TaperSideWidth = TaperWidth + 2 * extention

                TaperPosXmin = [max_MMIL, min_MMIL, min_MMIL]
                TaperPosXmax = [max_MMIL + TaperLength, min_MMIL - TaperLength, min_MMIL - TaperLength]

                PosOffset = [0, (posOffset / 2 + WG_Width / 2), -(posOffset / 2 + WG_Width / 2)]
                TaperPosYMax_BotR = [(0 + OffsetInput + WG_W / 2), (-WG_W / 2), (-WG_W / 2)]
                TaperPosYMin_BotR = [(0 + OffsetInput - WG_W / 2), (WG_W / 2), (WG_W / 2)]
                TaperPosYMax_TopR = [(0 + OffsetInput + WG_Width / 2), (-WG_Width / 2), (-WG_Width / 2)]
                TaperPosYMin_TopR = [(0 + OffsetInput - WG_Width / 2), (+WG_Width / 2), (+WG_Width / 2)]

                TaperPosYMax_BotL = [(0 + OffsetInput + TaperSideWidth / 2), (-TaperSideWidth / 2),
                                     (-TaperSideWidth / 2)]
                TaperPosYMin_BotL = [(0 + OffsetInput - TaperSideWidth / 2), (+TaperSideWidth / 2),
                                     (+TaperSideWidth / 2)]
                TaperPosYMax_TopL = [(0 + OffsetInput + TaperWidth / 2), (-TaperWidth / 2), (-TaperWidth / 2)]
                TaperPosYMin_TopL = [(0 + OffsetInput - TaperWidth / 2), (+TaperWidth / 2), (+TaperWidth / 2)]

                for i in range(len(xPos)):
                    TaperZmin = max_slabH
                    TaperZmax = max_slabH + max_MMIH

                    TaperXmin = TaperPosXmin[i]
                    TaperXmax = TaperPosXmax[i]

                    ymin_bot_l = TaperPosYMin_BotL[i]
                    ymax_bot_l = TaperPosYMax_BotL[i]

                    ymin_bot_r = TaperPosYMin_BotR[i]
                    ymax_bot_r = TaperPosYMax_BotR[i]

                    ymin_top_l = TaperPosYMin_TopL[i]
                    ymax_top_l = TaperPosYMax_TopL[i]

                    ymin_top_r = TaperPosYMin_TopR[i]
                    ymax_top_r = TaperPosYMax_TopR[i]

                    vtx = np.array([[TaperXmin, ymin_bot_l, TaperZmin],  # 1
                                    [TaperXmax, ymin_bot_r, TaperZmin],  # 2
                                    [TaperXmax, ymax_bot_r, TaperZmin],  # 3
                                    [TaperXmin, ymax_bot_l, TaperZmin],  # 4
                                    [TaperXmin, ymin_top_l, TaperZmax],  # 5
                                    [TaperXmax, ymin_top_r, TaperZmax],  # 6
                                    [TaperXmax, ymax_top_r, TaperZmax],  # 7
                                    [TaperXmin, ymax_top_l, TaperZmax],  # 8
                                    ])
                    a = [[np.array([[1, 4, 3, 2]], dtype=object)], [np.array([[1, 5, 8, 4]], dtype=object)],
                         [np.array([[1, 2, 6, 5]], dtype=object)], [np.array([[2, 6, 7, 3]], dtype=object)],
                         [np.array([[3, 4, 8, 7]], dtype=object)], [np.array([[5, 6, 7, 8]], dtype=object)]]



                    # Send Values to Lumerical and create solid
                    self.lum.putv('vertices', vtx)
                    self.lum.putv('facets', a)
                    self.lum.addplanarsolid(vtx, a)
                    self.lum.set('material', MaterialWG)
                    self.lum.set('name', TapersNames[i])
                    self.lum.set('y', PosOffset[i])

                # create loop
                for i in range(len(xPos)):
                    self.lum.addwaveguide()
                    self.lum.set("name", names[i])
                    self.lum.set("x", xPos[i])
                    self.lum.set("y", yPos[i])
                    self.lum.set("z", z_Offset)
                    self.lum.set("base width", WG_Width_top)
                    self.lum.set("base height", max_MMIH)
                    self.lum.set("base angle", 90 - angle)
                    pole = np.array([[maxWGL[i], 0], [minWGL[i], 0]])
                    self.lum.set("poles", pole)
                    self.lum.set("material", MaterialSlab)
                    
                    
            self.lum.select("MMI")
            self.lum.addtogroup("MMI Object")
            self.lum.select("Input WG")
            self.lum.addtogroup("MMI Object")
            self.lum.select("Output WG_L")
            self.lum.addtogroup("MMI Object")
            self.lum.select("Output WG_R")
            self.lum.addtogroup("MMI Object")
            self.lum.select("Taper Input WG")
            self.lum.addtogroup("MMI Object")
            self.lum.select("Taper Output WG_L")
            self.lum.addtogroup("MMI Object")
            self.lum.select("Taper Output WG_R")
            self.lum.addtogroup("MMI Object")
            
            if Cladding == True:
                # create_cover
                self.lum.addrect()
                self.lum.set("name", "cladding")
                self.lum.set("material", MaterialClad)
                self.lum.set("y min", min_subW)
                self.lum.set("y max", max_subW)
                self.lum.set("z", z_Offset)
                self.lum.set("z span", max_MMIH + 0.4e-6)
                self.lum.set("x min", min_subL)
                self.lum.set("x max", max_subL)
                self.lum.set("override mesh order from material database", True)
                self.lum.set("mesh order", 4)
                self.lum.set("alpha", 0.7)
                
                self.lum.select("cladding")
                self.lum.addtogroup("MMI Object")
            else:
                pass

        else:
            raise ValueError(
                "Incorect Taper input. Taper must be an boolen. You can choose from Taper = True or Taper = False!")






    def DirectionalCoupler(self, Parameters):
        '''


        Parameters
        ----------
        Parameters : Dictionary
            Dictionary with all the data needet for the Directional Coupler. 
            Parameters['Material'] : list of str
                List of Materials. The list should be with names (str) of a valid Lumerical materials.
                Check the names in Lumerical Materials viewer.
                The List of materials must contain at least 2 materials! 
                Parameters['Material'] = ['Cladding/Substrat', 'Object Material'].
                For Example: Parameters['Material'] = ["SiO2 (Glass) - Palik", 'LiNbO3_20deg_X cut'].
            Parameters['Substrate Height'] : int/float
                Substrate height.
            Parameters['Substrate Width'] : int/float
                Substrate Width.
            Parameters['DC Length'] : int/float
                Length of the directional coupler
            Parameters['angle'] : int/float
                Angle of the Waveguide Walls. it is calculated WG_angle = 90 - angle.
                For anfle = 90 we get a perfect rect!
            Parameters['WG Height'] : int/float
                Waveguide hight. Also the height of the MMI section
            Parameters['WG Width'] : int/float
                Waveguide width.
            Parameters['Position Offset'] : int/float
                Offset between the waveguides. The minimum distance between Waveguides is 1 um
                becouse of manufactering restrictions in the University.
            Parameters['Slab Height'] : int/float
                Height of the slab.
            Parameters['Wavelength'] : int/float
                Wavelength
            Parameters["Cladding"] : anything, optional
                This function will check if you have set Parameters["Cladding"] to anaything, for example "Parameters["Cladding"]=1" 
                and if so it will put cladding over your structure. If the user didnt give the "Cladding" as dictionary key no cladding 
                will be set.


        Raises
        ------
        ValueError
            DESCRIPTION.

        Returns
        -------
        None.

        '''

        Material = Parameters['Material']
        Substrate_Width = Parameters['Substrate Width']
        Substrate_Height = Parameters['Substrate Height']
        DC_Lenght = Parameters['DC Length']
        angle = Parameters['angle']
        WG_Height = Parameters['WG Height']
        WG_Width = Parameters['WG Width']
        posOffset = Parameters['Position Offset']
        Slab_Height = Parameters['Slab Height']

        if "Cladding" in list(Parameters.keys()):
            Cladding = True
        else:
            Cladding = False

        # Material definition
        if len(Material) < 2:
            raise ValueError(
                "List of materials must contain at least 2 materials!, Parameters['Material'] = ['Cladding/Substrat', 'Object Material']")
        else:
            MaterialSub = Material[0]
            MaterialClad = Material[0]
            MaterialSlab = Material[1]
            MaterialWG = MaterialSlab
            print(Material)




        # Device specifications
        Device_Length = DC_Lenght
        Device_Width = Substrate_Width

        # creating the substrate
        max_subH = Substrate_Height
        min_subH = -Substrate_Height
        max_subW = Device_Width / 2
        min_subW = -Device_Width / 2
        min_subL = -Device_Length / 2
        max_subL = Device_Length / 2

        self.lum.addrect()
        self.lum.set("name", "Substrate")
        self.lum.set("y min", min_subW)
        self.lum.set("y max", max_subW)
        self.lum.set("z min", min_subH)
        self.lum.set("z max", max_subH)
        self.lum.set("x min", min_subL)
        self.lum.set("x max", max_subL)
        self.lum.set("material", MaterialSub)
        
        
        if Slab_Height == 0:
            z_Offset = max_subH + WG_Height / 2
            
        else:
            # creating the thin film
            min_slabH = max_subH
            max_slabH = max_subH + Slab_Height

            self.lum.addrect()
            self.lum.set("name", "Slab")
            self.lum.set("y min", min_subW)
            self.lum.set("y max", max_subW)
            self.lum.set("z min", min_slabH)
            self.lum.set("z max", max_slabH)
            self.lum.set("x min", min_subL)
            self.lum.set("x max", max_subL)
            self.lum.set("material", MaterialSlab)
            z_Offset = max_slabH + WG_Height / 2
            
            self.lum.select("Slab")
            self.lum.addtogroup("Directional Coupler")

        # Positions of the Input and Output WGs
        # Triangle EQ for MMI Width
        x = abs(WG_Height / (np.cos((angle) * np.pi / 180)))  # in Radians
        extention = np.sqrt(x ** 2 - WG_Height ** 2)
        WG_W = WG_Width + 2 * extention
        WG_Width_bott = WG_W
        
        offset_WG = posOffset / 2 + WG_Width / 2 + WG_W / 2

        if offset_WG > Device_Width / 2:
            self.lum.deleteall()
            raise ValueError('You are Trying to move the Waveguide outside the MMI. This is not possible!')
        else:
            # Mirror the In and Out WG on both sides
            maxWGL = [DC_Lenght / 2, DC_Lenght / 2]
            minWGL = [-DC_Lenght / 2, -DC_Lenght / 2]
            xPos = [0, 0]
            yPos = [WG_Width / 2 + posOffset / 2, - WG_Width / 2 - posOffset / 2]

            # Names of the WGs
            names = ['Top WG', 'Bottom WG']

            # create loop
            for i in range(len(xPos)):
                self.lum.addwaveguide()
                self.lum.set("name", names[i])
                self.lum.set("x", xPos[i])
                self.lum.set("y", yPos[i])
                self.lum.set("z", z_Offset)
                self.lum.set("base width", WG_Width_bott)
                self.lum.set("base height", WG_Height)
                self.lum.set("base angle", 90 - angle)
                pole = np.array([[maxWGL[i], 0], [minWGL[i], 0]])
                self.lum.set("poles", pole)
                self.lum.set("material", MaterialWG)

            self.lum.select("Top WG")
            self.lum.addtogroup("Directional Coupler")
            self.lum.select("Bottom WG")
            self.lum.addtogroup("Directional Coupler")
            
            
            if Cladding == True:
                # create_cover
                self.lum.addrect()
                self.lum.set("name", "cladding")
                self.lum.set("material", MaterialClad)
                self.lum.set("y min", min_subW)
                self.lum.set("y max", max_subW)
                self.lum.set("z min", max_subH)
                self.lum.set("z max", 2e-6)
       
                self.lum.set("x min", min_subL)
                self.lum.set("x max", max_subL)
                self.lum.set("override mesh order from material database", True)
                self.lum.set("mesh order", 4)
                self.lum.set("alpha", 0.7)
                
                self.lum.select("cladding")
                self.lum.addtogroup("Directional Coupler")
            else:
                pass




    def WDM(self, Parameters):
        '''


        Parameters
        ----------
        Parameters : Dictionary
            Dictionary with all the data needet for the WDM. 
            Parameters['Material'] : list of str
                List of Materials. The list should be with names (str) of a valid Lumerical materials.
                Check the names in Lumerical Materials viewer.
                The List of materials must contain at least 2 materials! 
                Parameters['Material'] = ['Cladding/Substrat', 'Object Material'].
                For Example: Parameters['Material'] = ["SiO2 (Glass) - Palik", 'LiNbO3_20deg_X cut'].
            Parameters['Substrate Height'] : int/float
                Substrate height.
            Parameters['MMI Width'] : int/float
                Width of the MMI.
            Parameters['MMI Length'] : int/float
                Length of the MMI.
            Parameters['angle'] : int/float
                Angle of the Waveguide Walls. it is calculated WG_angle = 90 - angle.
                For anfle = 90 we get a perfect rect!
            Parameters['WG Height'] : int/float
                Waveguide hight. Also the height of the MMI section
            Parameters['WG Width'] : int/float
                Waveguide width.
            Parameters['Wavelength'] : int/float
                Wavelength.
            Parameters['WG Length'] : int/float
                Waveguide length.
            Parameters['Slab Height'] : int/float
                Height of the slab.
            Parameters['Angle Thetha'] :
                Input and output angle of the waveguide. This is only temporally
            Parameters['Taper Width'] : int/float
                Backside width of the Taper, frontside width is the waveguide width
            Parameters['Taper Length'] : int/float
                Length of the Taper.
            Parameters['Taper'] : boolen
                If Taper == False, no Taper will be placed with the Waveguides.
                If Taper == True, tapers will be placed with the waveguides
            Parameters["Cladding"] : anything, optional
                This function will check if you have set Parameters["Cladding"] to anaything, for example "Parameters["Cladding"]=1" 
                and if so it will put cladding over your structure. If the user didnt give the "Cladding" as dictionary key no cladding 
                will be set.

        Raises
        ------
        ValueError
            DESCRIPTION.

        Returns
        -------
        None.

        '''

        Material = Parameters['Material']
        Substrate_Height = Parameters['Substrate Height']
        MMI_Width = Parameters['MMI Width']
        MMI_Length = Parameters['MMI Length']
        WaveLength = Parameters['Wavelength']
        angle = Parameters['angle']
        WG_Height = Parameters['WG Height']
        WG_Width = Parameters['WG Width']
        WG_Length = Parameters['WG Length']
        Slab_Height = Parameters['Slab Height']
        angleTheta = Parameters['Angle Thetha']
        TaperLength = Parameters['Taper Length']
        TaperWidth = Parameters['Taper Width']
        Taper = Parameters['Taper']
        
        if "Cladding" in list(Parameters.keys()):
            Cladding = True
        else:
            Cladding = False

        # Material definition
        if len(Material) < 2:
            raise ValueError(
                "List of materials must contain at least 2 materials!, Parameters['Material'] = ['Cladding/Substrat', 'Object Material']")
        else:
            MaterialSub = Material[0]
            MaterialClad = Material[0]
            MaterialSlab = Material[1]
            MaterialWG = MaterialSlab



        if Taper == False:

            # Device specifications
            Device_Length = MMI_Length + 4 * WG_Length
            Device_Width = MMI_Width + 2*WG_Length + 3e-6  # MMI_Width

            # creating the substrate
            max_subH = Substrate_Height
            min_subH = -Substrate_Height


            self.lum.addrect()
            self.lum.set("name", "Substrate")
            self.lum.set("z min", min_subH)
            self.lum.set("z max", max_subH)
            self.lum.set("x", 0)
            self.lum.set("x span", Device_Length)
            self.lum.set("y", 0)
            self.lum.set("y span", Device_Width)
            self.lum.set("material", MaterialSub)
            self.lum.select("Substrate")
            self.lum.addtogroup("WDM")
            

            # creating the thin film
            min_slabH = max_subH
            max_slabH = max_subH + Slab_Height

            self.lum.addrect()
            self.lum.set("name", "Slab")
            self.lum.set("z min", min_slabH)
            self.lum.set("z max", max_slabH)
            self.lum.set("x", 0)
            self.lum.set("x span", Device_Length)
            self.lum.set("y", 0)
            self.lum.set("y span", Device_Width)
            self.lum.set("material", MaterialSlab)
            self.lum.select("Slab")
            self.lum.addtogroup("WDM")

            # creating the MMI
            # creating the MMI

            max_MMIH = WG_Height
            max_MMIL = MMI_Length / 2
            min_MMIL = -MMI_Length / 2
            z_Offset = max_slabH + max_MMIH / 2


            # Triangle EQ for MMI Width
            x = abs(max_MMIH / (np.cos((angle) * np.pi / 180)))  # in Radians
            extention = np.sqrt(x ** 2 - max_MMIH ** 2)
            MMI_Wid = MMI_Width + 2 * extention


            # self.lum.addwaveguide()
            # self.lum.set("name", "MMI")
            # self.lum.set("base height", max_MMIH)
            # self.lum.set("base angle", 90 - angle)
            # self.lum.set("x", 0)
            # self.lum.set("y", 0)
            # self.lum.set("z", z_Offset)
            # pole = np.array([[max_MMIL, 0], [min_MMIL, 0]])
            # self.lum.set("poles", pole)
            # self.lum.set("material", MaterialWG)
            # self.lum.set("base width", MMI_Wid)


            x = abs(max_MMIH / (np.cos((angle) * np.pi / 180)))  # in Radians
            extention = np.sqrt(x ** 2 - max_MMIH ** 2)
            WG_W = WG_Width + 2 * extention

            # Correction in Y-Axis
            difY = (TaperLength / 2) * np.cos(angleTheta * np.pi / 180)
            xLen = 2 * (TaperLength / 2) * np.sin((angleTheta / 2) * np.pi / 180)
            Diff = WG_W - WG_Width

            NewY = -MMI_Width / 2 + TaperWidth / 2 - np.sqrt((TaperLength / 2) ** 2 - difY ** 2)  # - xLen - Diff/2 #
            Li = (4 * 2.5 * MMI_Width ** 2) / WaveLength


            #Creaate MMI
            self.lum.addwaveguide()
            self.lum.set("name", "MMI")
            self.lum.set("base height", max_MMIH)
            self.lum.set("base angle", 90 - angle)
            self.lum.set("x", 0)
            self.lum.set("y", 0)
            self.lum.set("z", z_Offset)
            pole = np.array([[max_MMIL + WG_Length - xLen, 0], [min_MMIL - xLen / 2, 0]])
            self.lum.set("poles", pole)
            self.lum.set("material", MaterialWG)
            self.lum.set("base width", MMI_Wid)
            self.lum.select("MMI")
            self.lum.addtogroup("WDM")
            

            # Names of the WGs
            Waveguides = ['LN_Input', 'LN_Output']
            myscript = self.Script()

            self.lum.addstructuregroup()
            self.lum.set("name", Waveguides[0])
            self.lum.set("construction group", 1)
            self.lum.adduserprop("thickness", 2, WG_Height)
            self.lum.adduserprop("angle_side", 0, angle)
            self.lum.adduserprop("width_l", 2, WG_Width)
            self.lum.adduserprop("width_r", 2, WG_Width)
            self.lum.adduserprop("hfrac_ref", 0, 1)
            self.lum.adduserprop("len", 2, WG_Length)
            self.lum.adduserprop("material", 5, MaterialWG)
            self.lum.adduserprop("index", 0, 1)
            self.lum.set("script", myscript)
            self.lum.set("first axis", "z")
            self.lum.set("rotation 1", angleTheta)
            self.lum.set("x", -MMI_Length / 2 - TaperLength / 2)  # -MMI_Length / 2
            self.lum.set("z", z_Offset)
            self.lum.set("y", NewY)


            self.lum.addstructuregroup()
            self.lum.set("name", Waveguides[1])
            self.lum.set("construction group", 1)
            self.lum.adduserprop("thickness", 2, WG_Height)
            self.lum.adduserprop("angle_side", 0, angle)
            self.lum.adduserprop("width_l", 2, WG_Width)
            self.lum.adduserprop("width_r", 2, WG_Width)
            self.lum.adduserprop("hfrac_ref", 0, 1)
            self.lum.adduserprop("len", 2, WG_Length)
            self.lum.adduserprop("material", 5, MaterialWG)
            self.lum.adduserprop("index", 0, 1)
            self.lum.set("script", myscript)
            self.lum.set("first axis", "z")
            self.lum.set("rotation 1", angleTheta)
            self.lum.set("x", MMI_Length / 2 + TaperLength / 2)  # + MMI_Length / 2
            self.lum.set("z", z_Offset)
            self.lum.set("y", -NewY)
            
            
            
            if Cladding == True:
                # Create Cladding
                self.lum.addrect()
                self.lum.set("name", "cladding")
                self.lum.set("x", 0)
                self.lum.set("x span", Device_Length)
                self.lum.set("y", 0)
                self.lum.set("y span", Device_Width)
                self.lum.set('z min', max_slabH)
                self.lum.set('z max', max_slabH + 4 * WG_Height)
                self.lum.set("material", MaterialClad)
                self.lum.set("override mesh order from material database", True)
                self.lum.set("mesh order", 4)
                self.lum.set("alpha", 0.7)
                self.lum.select("cladding")
                self.lum.addtogroup("WDM")
                
                # self.lum.select("cladding")
                # self.lum.addtogroup("Directional Coupler")
            else:
                pass


     


        elif Taper == True:


            # Device specifications
            Device_Length = MMI_Length + 4 * WG_Length
            Device_Width = MMI_Width + 2*WG_Length + WaveLength * 2  

            # creating the substrate
            max_subH = Substrate_Height
            min_subH = -Substrate_Height

            # creating the thin film
            min_slabH = max_subH
            max_slabH = max_subH + Slab_Height

            max_MMIH = WG_Height
            max_MMIL = MMI_Length / 2
            min_MMIL = -MMI_Length / 2
            z_Offset = max_slabH + max_MMIH / 2


            # Triangle EQ for MMI Width
            x = abs(max_MMIH / (np.cos((angle) * np.pi / 180)))  # in Radians
            extention = np.sqrt(x ** 2 - max_MMIH ** 2)
            MMI_Wid = MMI_Width + 2 * extention

            TaperNames = ["Input_Taper", "Output_Taper"]
            Waveguides = ["Input_WG", "Output_WG"]

            x = abs(max_MMIH / (np.cos((angle) * np.pi / 180)))  # in Radians
            extention = np.sqrt(x ** 2 - max_MMIH ** 2)
            WG_W = WG_Width + 2 * extention

            # Correction in Y-Axis
            difY = (TaperLength/2) * np.cos(angleTheta*np.pi/180)
            xLen = 2*(TaperLength/2)* np.sin((angleTheta/2) * np.pi/180)
            Diff = WG_W - WG_Width

            NewY = -MMI_Width / 2 + TaperWidth/2- np.sqrt( (TaperLength/2)**2 - difY**2) # - xLen - Diff/2 #
            Li =( 4*2.5*MMI_Width**2)/WaveLength


            self.lum.addrect()
            self.lum.set("name", "Substrate")
            self.lum.set("z min", min_subH)
            self.lum.set("z max", max_subH)
            self.lum.set("x", 0)
            self.lum.set("x span", Device_Length)
            self.lum.set("y", 0)
            self.lum.set("y span", Device_Width)
            self.lum.set("material", MaterialSub)
            self.lum.select("Substrate")
            self.lum.addtogroup("WDM")
            



            self.lum.addrect()
            self.lum.set("name", "Slab")
            self.lum.set("z min", min_slabH)
            self.lum.set("z max", max_slabH)
            self.lum.set("x", 0)
            self.lum.set("x span", Device_Length)
            self.lum.set("y", 0)
            self.lum.set("y span", Device_Width)
            self.lum.set("material", MaterialSlab)
            self.lum.select("Slab")
            self.lum.addtogroup("WDM")
           

            # creating the MMI

            self.lum.addwaveguide()
            self.lum.set("name", "MMI")
            self.lum.set("base height", max_MMIH)
            self.lum.set("base angle", 90 - angle)
            self.lum.set("x", 0)
            self.lum.set("y", 0)
            self.lum.set("z", z_Offset)
            pole = np.array([[max_MMIL + WG_Length -  xLen , 0], [min_MMIL- xLen/2, 0]])
            self.lum.set("poles", pole)
            self.lum.set("material", MaterialWG)
            self.lum.set("base width", MMI_Wid)
            self.lum.select("MMI")
            self.lum.addtogroup("WDM")






            # Input and Output Tapers and WGs
            myscript = self.Script()


            self.lum.addstructuregroup()
            self.lum.set("name",TaperNames[0])
            self.lum.set("construction group",1)
            self.lum.adduserprop("thickness",2, WG_Height)
            self.lum.adduserprop("angle_side",0, angle)
            self.lum.adduserprop("width_l",2, WG_Width)
            self.lum.adduserprop("width_r",2, TaperWidth)
            self.lum.adduserprop("hfrac_ref",0,1)
            self.lum.adduserprop("len",2, TaperLength)
            self.lum.adduserprop("material",5,MaterialWG)
            self.lum.adduserprop("index",0,1)
            self.lum.set("script",myscript)
            self.lum.set("first axis", "z")
            self.lum.set("rotation 1",angleTheta)
            self.lum.set("x", -MMI_Length / 2 - TaperLength/2 ) #-MMI_Length / 2
            self.lum.set("z", z_Offset)
            self.lum.set("y", NewY)




            self.lum.addstructuregroup()
            self.lum.set("name",Waveguides[0])
            self.lum.set("construction group",1)
            self.lum.adduserprop("thickness",2, WG_Height)
            self.lum.adduserprop("angle_side",0, angle)
            self.lum.adduserprop("width_l",2, WG_Width)
            self.lum.adduserprop("width_r",2,WG_Width)
            self.lum.adduserprop("hfrac_ref",0,1)
            self.lum.adduserprop("len",2, WG_Length)
            self.lum.adduserprop("material",5,MaterialWG)
            self.lum.adduserprop("index",0,1)
            self.lum.set("script",myscript)
            self.lum.set("first axis", "z")
            self.lum.set("rotation 1",angleTheta)
            self.lum.set("x", -MMI_Length / 2 - TaperLength/2 ) #-MMI_Length / 2
            self.lum.set("z", z_Offset)
            self.lum.set("y", NewY )



            self.lum.addstructuregroup()
            self.lum.set("name",TaperNames[1])
            self.lum.set("construction group",1)
            self.lum.adduserprop("thickness",2, WG_Height)
            self.lum.adduserprop("angle_side",0, angle)
            self.lum.adduserprop("width_l",2, TaperWidth)
            self.lum.adduserprop("width_r",2,WG_Width)
            self.lum.adduserprop("hfrac_ref",0,1)
            self.lum.adduserprop("len",2, TaperLength)
            self.lum.adduserprop("material",5,MaterialWG)
            self.lum.adduserprop("index",0,1)
            self.lum.set("script",myscript)
            self.lum.set("first axis", "z")
            self.lum.set("rotation 1",angleTheta)
            self.lum.set("x", MMI_Length / 2 + TaperLength/2 )  #+ MMI_Length / 2
            self.lum.set("z", z_Offset)
            self.lum.set("y", -NewY )



            self.lum.addstructuregroup()
            self.lum.set("name",Waveguides[1])
            self.lum.set("construction group",1)
            self.lum.adduserprop("thickness",2, WG_Height)
            self.lum.adduserprop("angle_side",0, angle)
            self.lum.adduserprop("width_l",2, WG_Width)
            self.lum.adduserprop("width_r",2,WG_Width)
            self.lum.adduserprop("hfrac_ref",0,1)
            self.lum.adduserprop("len",2, WG_Length)
            self.lum.adduserprop("material",5,MaterialWG)
            self.lum.adduserprop("index",0,1)
            self.lum.set("script",myscript)
            self.lum.set("first axis", "z")
            self.lum.set("rotation 1",angleTheta)
            self.lum.set("x", MMI_Length / 2 + TaperLength/2 )  #+ MMI_Length / 2
            self.lum.set("z", z_Offset )
            self.lum.set("y", -NewY)
            
            
            
            if Cladding == True:
                # Create Cladding
                self.lum.addrect()
                self.lum.set("name", "cladding")
                self.lum.set("x", 0)
                self.lum.set("x span", Device_Length)
                self.lum.set("y", 0)
                self.lum.set("y span", Device_Width)
                self.lum.set('z min', max_slabH)
                self.lum.set('z max', max_slabH + 4 * WG_Height)
                self.lum.set("material", MaterialClad)
                self.lum.set("override mesh order from material database", True)
                self.lum.set("mesh order", 4)
                self.lum.set("alpha", 0.7)
                self.lum.select("cladding")
                self.lum.addtogroup("WDM")

            else:
                pass

            






    def InverseTaper(self, Parameters):
        '''


        Parameters
        ----------
        Parameters : Dictionary
            Dictionary with all the data needet for the InverseTaper.
            Parameters['Material'] : list of str
                List of Materials. The list should be with names (str) of a valid Lumerical materials.
                Check the names in Lumerical Materials viewer.
                The List of materials must contain at least 2 materials! 
                For this structure 3 Material will be needed
                Parameters['Material'] = ['Cladding/Substrat', 'Object Material'].
                For Example: Parameters['Material'] = ["SiO2 (Glass) - Palik", 'LiNbO3_20deg_X cut'].
            Parameters['Substrate Height'] : int/float
                Substrate height.
            Parameters['angle'] : int/float
                Angle of the Waveguide Walls. it is calculated WG_angle = 90 - angle.
                For anfle = 90 we get a perfect rect!
            Parameters['WG Height'] : int/float
                Waveguide hight. Also the height of the MMI section
            Parameters['WG Width'] : int/float
                Waveguide width. Also in this function and ONLY in this function this will be the
                ibverse Taper width!!!
            Parameters['Slab Height'] : int/float
                Slab height
            Parameters['Taper Length'] : int/float
                Inverse Taper Length
            Parameters['Taper Width'] : int/float
                Front Width of the inverse Taper!! In this Function and ONLY in this function, the
                Parameters['Taper Width'] is the frond width of the inverse Taper!
            Parameters['PWB Taper Width Back'] : int/float
                Photonic Wirebonding (PWB) Width back side (to the Photonic Wire Bonding)
            Parameters['PWB Taper Hight Back'] : int/float
                Photonic Wire Bonding Height back side (to the Photonic Wire Bonding)
            Parameters['PWB Taper Width Front'] : int/float
                Photonic Wirebonding (PWB) Width front side (to the photonic waveguide)
            Parameters['PWB Taper Hight Front'] : int/float
                Photonic Wire Bonding Height front side (to the photonic waveguide)
            Parameters['PWB Taper Length'] : int/float
                Length of the Photonic Wire Bonding Taper
            Parameters["SMF Core Diameter"] : int/float
                Single Mode Fiber core Diameter
            Parameters["SMF Cladding Diameter"] : int/float
                Single Mode Fiber Cladding Diameter
            Parameters["SMF Core Index"]
                Single Mode Fiber Core Index
            Parameters["SMF Cladding Index"]
                Single Mode Fiber Cladding Index
                
        Raises
        ------
        ValueError
            DESCRIPTION.

        Returns
        -------
        None.

        '''

        Material = Parameters['Material']
        Substrate_Height = Parameters['Substrate Height']
        angle = Parameters['angle']
        WG_Height = Parameters['WG Height']
        WG_Width = Parameters['WG Width']
        WG_Length = Parameters['WG Length']
        Slab_Height = Parameters['Slab Height']
        TaperLength = Parameters['Taper Length']
        TaperWidth = Parameters['Taper Width']
        TaperWidthF = Parameters['PWB Taper Width Front']
        TaperWidthB = Parameters['PWB Taper Width Back']
        TaperHightB = Parameters['PWB Taper Hight Back']
        TaperHightF = Parameters['PWB Taper Hight Front']
        TaperLength_PWB = Parameters['PWB Taper Length']
        
        # SMF Parameters
        CoreDiameter = Parameters["SMF Core Diameter"]
        CladdingDiameter = Parameters["SMF Cladding Diameter"]
        CoreIndex = Parameters["SMF Core Index"]
        CladdingIndex = Parameters["SMF Cladding Index"]

     


        # Material definition
        if len(Material) < 3:
            raise ValueError(
                "List of materials must contain at least 3 materials!, Parameters['Material'] = ['Cladding/Substrat', 'Object Material', 'Material Photonic Wire Bonding']")
        else:
            MaterialSub = Material[0]
            MaterialClad = Material[0]
            MaterialSlab = Material[1]
            MaterialWG = MaterialSlab
            MaterialPWB = Material[2]



        # creating the substrate
        max_subH = Substrate_Height/2
        min_subH = -Substrate_Height/2

        # make substrate
        self.lum.addrect()
        self.lum.set("name", "Substrate InverseTaper")
        self.lum.set("override mesh order from material database", 1)
        self.lum.set("mesh order", 4)
        self.lum.set("y", 0)
        self.lum.set("y span", TaperWidthB * 2)
        self.lum.set("z min", min_subH)
        self.lum.set("z max", max_subH)
        self.lum.set("x min", -TaperLength_PWB/2)
        self.lum.set("x max", TaperLength_PWB/2 + WG_Length )
        # self.lum.set("x", 10e-6)
        # self.lum.set("x span", TaperLength_PWB * 2)
        self.lum.set("material", MaterialSub)

        # creating the thin film
        min_slabH = max_subH
        max_slabH = max_subH + Slab_Height

        # Create Slab and Check if Slab is needed
        if Slab_Height == 0:
            pass
        else:
            self.lum.addrect()
            self.lum.set("name", "LN_slab")
            self.lum.set("override mesh order from material database", 1)
            self.lum.set("mesh order", 4)
            self.lum.set("y", 0e-12)
            self.lum.set("y span", TaperWidthB * 2)
            self.lum.set("z min", min_slabH)
            self.lum.set("z max", max_slabH)
            self.lum.set("x min", -TaperLength_PWB / 2 - 0.1e-6)
            self.lum.set("x max", TaperLength_PWB/2 + WG_Length)
            # self.lum.set("x", 10e-6)
            # self.lum.set("x span", TaperLength_PWB * 2)
            self.lum.set("material", MaterialSlab)

        # Names of the WGs
        TapersNames = ['Taper_PWB', 'Inverse_Taper']

        # Taper sideangle widths
        x = abs(WG_Height / (np.cos((angle) * np.pi / 180)))  # in Radians
        extention = np.sqrt(x ** 2 - WG_Height ** 2)
        TaperSideWidth = TaperWidth + 2 * extention
        WG_W = WG_Width + 2 * extention




        # PWD Taper Y-Parameters
        PWB_TaperPosYMax_BotR = [(TaperWidthF / 2)]
        PWB_TaperPosYMin_BotR = [(- TaperWidthF / 2)]
        PWB_TaperPosYMax_TopR = [(TaperWidthF / 2)]
        PWB_TaperPosYMin_TopR = [(- TaperWidthF / 2)]
        PWB_TaperPosYMax_BotL = [(TaperWidthB / 2)]
        PWB_TaperPosYMin_BotL = [(- TaperWidthB / 2)]
        PWB_TaperPosYMax_TopL = [(TaperWidthB / 2)]
        PWB_TaperPosYMin_TopL = [(- TaperWidthB / 2)]

        # Inverse   Taper Y-Parameters
        TaperPosYMax_BotR = [(WG_W / 2)]
        TaperPosYMin_BotR = [(- WG_W / 2)]
        TaperPosYMax_TopR = [(WG_Width / 2)]
        TaperPosYMin_TopR = [(- WG_Width / 2)]
        TaperPosYMax_BotL = [(TaperSideWidth / 2)]
        TaperPosYMin_BotL = [(- TaperSideWidth / 2)]
        TaperPosYMax_TopL = [(TaperWidth / 2)]
        TaperPosYMin_TopL = [(- TaperWidth / 2)]



        if Slab_Height == 0:
            z_Offset = max_subH
        else:
            z_Offset = max_slabH

        # PWB Taper Hights
        TaperZmin = z_Offset
        TaperZmaxF =  TaperHightF + TaperZmin
        TaperZmaxB =  TaperHightB + TaperZmin

        # Inverse Taper Hights
        TaperZmin = z_Offset
        TaperZmax = z_Offset + WG_Height

        # PWB Taper Length
        PWB_TaperXmin = -TaperLength_PWB / 2
        PWB_TaperXmax = TaperLength_PWB / 2

        # Inverse Taper Length
        TaperXmin = -TaperLength / 2
        TaperXmax = TaperLength / 2
        
        # Create Inverse Taper
        ymin_bot_l = TaperPosYMin_BotL[0]
        ymax_bot_l = TaperPosYMax_BotL[0]

        ymin_bot_r = TaperPosYMin_BotR[0]
        ymax_bot_r = TaperPosYMax_BotR[0]

        ymin_top_l = TaperPosYMin_TopL[0]
        ymax_top_l = TaperPosYMax_TopL[0]

        ymin_top_r = TaperPosYMin_TopR[0]
        ymax_top_r = TaperPosYMax_TopR[0]

        vtx = np.array([[TaperXmin, ymin_bot_l, TaperZmin],  # 1
                        [TaperXmax, ymin_bot_r, TaperZmin],  # 2
                        [TaperXmax, ymax_bot_r, TaperZmin],  # 3
                        [TaperXmin, ymax_bot_l, TaperZmin],  # 4
                        [TaperXmin, ymin_top_l, TaperZmax],  # 5
                        [TaperXmax, ymin_top_r, TaperZmax],  # 6
                        [TaperXmax, ymax_top_r, TaperZmax],  # 7
                        [TaperXmin, ymax_top_l, TaperZmax],  # 8
                        ])
        a = [[np.array([[1, 4, 3, 2]], dtype=object)], [np.array([[1, 5, 8, 4]], dtype=object)],
             [np.array([[1, 2, 6, 5]], dtype=object)], [np.array([[2, 6, 7, 3]], dtype=object)],
             [np.array([[3, 4, 8, 7]], dtype=object)], [np.array([[5, 6, 7, 8]], dtype=object)]]

        Offset_InvTaper = PWB_TaperXmin + TaperXmax

        if Offset_InvTaper - TaperXmax / 2 > PWB_TaperXmax:
            raise ValueError(
                "Inverse Taper is moved outside the PWB Taper! The maximal x-Offset of the Inverse Taper is Parameters['Offset Inverse Taper'] = " + str(
                    PWB_TaperXmax + TaperXmax))
        else:
            # Send Values to Lumerical and create solid
            self.lum.putv('vertices', vtx)
            self.lum.putv('facets', a)
            self.lum.addplanarsolid(vtx, a)
            self.lum.set('material', MaterialWG)
            self.lum.set('name', TapersNames[1])
            self.lum.set('second axis', 'z')
            self.lum.set('rotation 2', 0)
            self.lum.set('x', Offset_InvTaper)

        # Create PWB Taper
        ymin_bot_l = PWB_TaperPosYMin_BotL[0]
        ymax_bot_l = PWB_TaperPosYMax_BotL[0]

        ymin_bot_r = PWB_TaperPosYMin_BotR[0]
        ymax_bot_r = PWB_TaperPosYMax_BotR[0]

        ymin_top_l = PWB_TaperPosYMin_TopL[0]
        ymax_top_l = PWB_TaperPosYMax_TopL[0]

        ymin_top_r = PWB_TaperPosYMin_TopR[0]
        ymax_top_r = PWB_TaperPosYMax_TopR[0]

        vtx = np.array([[PWB_TaperXmin, ymin_bot_l, TaperZmin],  # 1
                        [PWB_TaperXmax, ymin_bot_r, TaperZmin],  # 2
                        [PWB_TaperXmax, ymax_bot_r, TaperZmin],  # 3
                        [PWB_TaperXmin, ymax_bot_l, TaperZmin],  # 4
                        [PWB_TaperXmin, ymin_top_l, TaperZmaxB],  # 5
                        [PWB_TaperXmax, ymin_top_r, TaperZmaxF],  # 6
                        [PWB_TaperXmax, ymax_top_r, TaperZmaxF],  # 7
                        [PWB_TaperXmin, ymax_top_l, TaperZmaxB],  # 8
                        ])
        a = [[np.array([[1, 4, 3, 2]], dtype=object)], [np.array([[1, 5, 8, 4]], dtype=object)],
             [np.array([[1, 2, 6, 5]], dtype=object)], [np.array([[2, 6, 7, 3]], dtype=object)],
             [np.array([[3, 4, 8, 7]], dtype=object)], [np.array([[5, 6, 7, 8]], dtype=object)]]

        # Send Values to Lumerical and create solid
        self.lum.putv('vertices', vtx)
        self.lum.putv('facets', a)
        self.lum.addplanarsolid(vtx, a)
        self.lum.set('material', MaterialPWB)
        self.lum.set("override mesh order from material database",1)
        self.lum.set("mesh order",3)
        self.lum.set('name', TapersNames[0])

        # Make Sqered WG-Extention of the PWB for mode Calculations

        # Extra Waveguide Lenght
        Ext_WGLength = 1e-6
        # PWD_x_Offset = PWB_TaperXmin

        # self.lum.addrect()
        # self.lum.set('x min', PWD_x_Offset - Ext_WGLength)
        # self.lum.set('x max',PWD_x_Offset)
        # self.lum.set('y', 0)
        # self.lum.set('y span', TaperWidthB)
        # self.lum.set('z min', TaperZmin)
        # self.lum.set('z max', TaperZmaxB)
        # self.lum.set("override mesh order from material database", 1)
        # self.lum.set("mesh order", 3)
        # self.lum.set('name', 'WG_Extention_PWB')
        # self.lum.set('material', MaterialPWB)

        # Make Sqred WG-Extention for the inverse Taper
        x_min = (Offset_InvTaper + TaperXmax)
        x_max =  PWB_TaperXmax + Ext_WGLength
        center = TaperXmax + Offset_InvTaper

        ZOffset_WG = TaperZmin+ WG_Height/2

        self.lum.addwaveguide()
        self.lum.set("name", 'WG_Extention_Inverse_Taper')
        self.lum.set("x", center )
        self.lum.set("y", 0)
        self.lum.set('z', ZOffset_WG)
        self.lum.set("base width", WG_W)
        self.lum.set("base height", WG_Height)
        self.lum.set("base angle", 90 - angle)
        pole = np.array([[0, 0], [WG_Length, 0]])
        self.lum.set("poles", pole)
        self.lum.set("material", MaterialSlab)
        
        # Create Taper abjects
        self.lum.select("WG_Extention_Inverse_Taper")
        self.lum.addtogroup('InverseTaper')
        self.lum.select(TapersNames[1])
        self.lum.addtogroup('InverseTaper')
        
        
 
        
        
        # Create the SMF 
        
        self.lum.addcircle()
        self.lum.set("name", "core")
        self.lum.set("override mesh order from material database", 1)
        self.lum.set("mesh order", 5)
        self.lum.set("first axis", "y")
        self.lum.set("rotation 1", 90)
        self.lum.set("alpha", 1)
        self.lum.set("radius", CoreDiameter/2)
        self.lum.set("index", CoreIndex)
        self.lum.set("x", PWB_TaperXmin)
        self.lum.set("y", 0)
        self.lum.set("z", CoreDiameter/2 + z_Offset)
        self.lum.set("z span", 2e-6)




        self.lum.addcircle()
        self.lum.set("name", "cladding")
        self.lum.set("override mesh order from material database", 1)
        self.lum.set("mesh order", 6)
        self.lum.set("first axis", "y")
        self.lum.set("rotation 1", 90)
        self.lum.set("alpha", 0.35)
        self.lum.set("radius", CladdingDiameter/2)
        self.lum.set("index", CladdingIndex)
        self.lum.set("x", PWB_TaperXmin )
        self.lum.set("y", 0)
        self.lum.set("z", CoreDiameter/2 + z_Offset)
        self.lum.set("z span", 2e-6)
        
        
        self.lum.select("core")
        self.lum.addtogroup('SMF')
        self.lum.select("cladding")
        self.lum.addtogroup('SMF')
        self.lum.select('SMF')
        self.lum.set("x", -1e-6)


    
        
        
        
    def CascadetMMI(self, Parameters, SpaceX, SpaceY):
        '''
        

        Parameters
        ----------
        Parameters : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        '''
        
        Material = Parameters['Material']
        Substrate_Height = Parameters['Substrate Height']
        MMI_Width = Parameters['MMI Width']
        MMI_Length = Parameters['MMI Length']
        angle = Parameters['angle']
        WG_Height = Parameters['WG Height']
        WG_Width = Parameters['WG Width']
        WG_Length = Parameters['WG Length']
        OffsetInput = Parameters['Offset Input']
        posOffset = Parameters['Position Offset']
        Slab_Height = Parameters['Slab Height']
        WaveLength = Parameters['Wavelength']
        TaperLength = Parameters['Taper Length']
        TaperWidth = Parameters['Taper Width']
        Taper = Parameters['Taper']
        x_span = Parameters["x span"]
        y_span = Parameters["y span"]
        polesList = Parameters["poles"]
        
        
        # Material definition
        if len(Material) < 2:
            raise ValueError(
                "List of materials must contain at least 2 materials!, Parameters['Material'] = ['Cladding/Substrat', 'Object Material']")
        else:
            MaterialSub = Material[0]
            MaterialClad = Material[0]
            MaterialSlab = Material[1]
            MaterialWG = MaterialSlab



        # Device specifications
        Device_Length = (2*MMI_Length + 2 * WG_Length + 2*SpaceX )
        Device_Width = (3*MMI_Width + 2 * WaveLength) + 3*SpaceY   # MMI_Width

        # creating the substrate
        max_subH = Substrate_Height
        min_subH = -Substrate_Height
        min_subL = -3*(Device_Length /2)
        max_subL = 3*(Device_Length / 2)
        

        self.lum.addrect()
        self.lum.set("name", "Substrate Main")
        self.lum.set("y", 0)
        self.lum.set("y span", Device_Width)
        self.lum.set("z min", min_subH)
        self.lum.set("z max", max_subH)
        self.lum.set("x min", min_subL)
        self.lum.set("x max", max_subL)
        self.lum.set("material", MaterialSub)
        
        # creating the MMI
        max_MMIH = WG_Height
        max_MMIL = MMI_Length / 2
        min_MMIL = -MMI_Length / 2
        
        if Slab_Height == 0:
            z_Offset = max_subH + max_MMIH / 2
        else:
            # creating the thin film
            min_slabH = max_subH
            max_slabH = max_subH + Slab_Height

            self.lum.addrect()
            self.lum.set("name", "Slab Main")
            self.lum.set("y", 0)
            self.lum.set("y span", Device_Width)
            self.lum.set("z min", min_slabH)
            self.lum.set("z max", max_slabH)
            self.lum.set("x min", min_subL)
            self.lum.set("x max", max_subL)
            self.lum.set("material", MaterialSlab)

            z_Offset = max_slabH + max_MMIH / 2
            self.lum.select("Slab")
            self.lum.addtogroup("Cascadet MMI")
            

        # Triangle EQ for MMI Width
        x = abs(max_MMIH / (np.cos((angle) * np.pi / 180)))  # in Radians
        extention = np.sqrt(x ** 2 - max_MMIH ** 2)
        MMI_Wid = MMI_Width + 2 * extention
        
        Parameters__Position_X = [0, -SpaceX -MMI_Length  - 2*WG_Length , -SpaceX -MMI_Length  - 2*WG_Length]
        Parameters__Position_Y = [0, (MMI_Wid + SpaceY)/2,  -(MMI_Wid + SpaceY)/2]
        

        names_MMI = ["MMI In", "MMI Out1", "MMI Out2"]
        for i in range(3):
            self.lum.addwaveguide()
            self.lum.set("name", names_MMI[i])
            self.lum.set("x", Parameters__Position_X[i])
            self.lum.set("y", Parameters__Position_Y[i])
            self.lum.set("z", z_Offset)
            self.lum.set("base height", max_MMIH)
            self.lum.set("base angle", 90 - angle)
            pole = np.array([[max_MMIL, 0], [min_MMIL, 0]])
            self.lum.set("poles", pole)
            self.lum.set("material", MaterialWG)
            self.lum.set("base width", MMI_Wid)

        # Positions of the Input and Output WGs
        # Triangle EQ for MMI Width
        x = abs(max_MMIH / (np.cos((angle) * np.pi / 180)))  # in Radians
        extention = np.sqrt(x ** 2 - max_MMIH ** 2)
        WG_W = WG_Width + 2 * extention
        WG_Width_top = WG_W
        OffMax = MMI_Width / 2

        offset_Taper = posOffset / 2 + WG_Width / 2 + TaperWidth / 2  # + WG_W / 2
        BotCornerDistance = posOffset/2 - TaperWidth / 2
        offset_WG = posOffset / 2 + WG_Width / 2 + WG_W / 2
        offset_WG2 = posOffset / 2

        if offset_WG2 < 0.5e-6:
            self.lum.deleteall()
            raise ValueError('The distance between the Tapers is less then 1 um !')
        else:

            if Taper == False:

                if offset_WG > OffMax:
                    self.lum.deleteall()
                    raise ValueError('You are Trying to move the Waveguide outside the MMI. This is not possible!')
                else:
                    # Mirror the In and Out WG on both sides
                    maxWGL = [WG_Length, 0, 0]
                    minWGL = [0, -WG_Length, -WG_Length]
                    xPos = [max_MMIL, min_MMIL, min_MMIL]
                    yPos = [0 + OffsetInput, WG_Width / 2 + posOffset / 2, - WG_Width / 2 - posOffset / 2]
                    
                    
                                
                    # Create two Bends Cos or Bezier depending on the poles option
                    x_span = SpaceX
                    y_span = (Parameters__Position_Y[1] + yPos[0]) - (Parameters__Position_Y[0] + yPos[1] ) 
                    
                    
                    if polesList == False:
                        pole = np.array([[0, 0], [x_span / 2, 0], [x_span / 2, y_span], [x_span, y_span]])
                        pole1 = np.array([[0, y_span], [x_span / 2, y_span], [x_span / 2, 0], [x_span,0]])
                    elif polesList == True:
                        K = ((np.pi -2)/np.pi)*100
                        pole = np.array([[0, 0], [0+(x_span*(K/100)), 0], [(x_span)-(x_span*(K/100)), y_span], [x_span, y_span]])
                        pole1 = np.array([[0, y_span], [0+(x_span*(K/100)), y_span], [(x_span)-(x_span*(K/100)), 0], [x_span, 0]])
                    else:
                        raise ValueError('Parameters["poles"] should be an boolen variable! Please set Parameters["poles"] to True if Bezier Curves needed. Set Parameters["poles"] to False for Cosinus Curve!')
        
        
                    names = ["S-Bend Top", "S-Bend Bot"]
                    S_Bend_X_Offset = [(-SpaceX - MMI_Length/2  - WG_Length), (SpaceX/2 - SpaceX/2 - MMI_Length/2 - WG_Length)]
                    S_Bend_Y_Offset = [( WG_Width / 2 + posOffset / 2) ,-( WG_Width / 2 + posOffset / 2)]
                    z_rotation = [0,180]
                    poles = [pole1, pole]
                    for i in range(2):
                        self.lum.addwaveguide()
                        self.lum.set("name", names[i])
                        self.lum.set("x", S_Bend_X_Offset[i])
                        self.lum.set("y", S_Bend_Y_Offset[i])
                        self.lum.set("z", z_Offset)
                        self.lum.set("base width", WG_W)
                        self.lum.set("base height", WG_Height)
                        self.lum.set("base angle", 90 - angle)
                        self.lum.set("poles", poles[i])
                        self.lum.set("material", MaterialWG)
                        self.lum.set("first axis","z")
                        self.lum.set("rotation 1",z_rotation[i])
                        
                        

                    # Names of the WGs
                    names1 = ['MMI In_WG', 'MMI Out_WG_Top', 'MMI Out_WG_Bot']
                    names2 = ['MMIOut1 In_WG', 'MMIOut1 Out_WG_Top', 'MMIOut1 Out_WG_Bot']
                    names3 = ['MMIOut2 In_WG', 'MMIOut2 Out_WG_Top', 'MMIOut2 Out_WG_Bot']
                    names = [names1, names2, names3]

                    # create loop
                    for j in range(3):
                        for i in range(len(xPos)):
                            self.lum.addwaveguide()
                            self.lum.set("name", names[j][i])
                            self.lum.set("x", Parameters__Position_X[j] + xPos[i])
                            self.lum.set("y", Parameters__Position_Y[j] + yPos[i])
                            self.lum.set("z", z_Offset)
                            self.lum.set("base width", WG_Width_top)
                            self.lum.set("base height", max_MMIH)
                            self.lum.set("base angle", 90 - angle)
                            pole = np.array([[maxWGL[i], 0], [minWGL[i], 0]])
                            self.lum.set("poles", pole)
                            self.lum.set("material", MaterialSlab)


                self.lum.select("MMI In")
                self.lum.addtogroup("Cascadet MMI")
                self.lum.select("MMI Out1")
                self.lum.addtogroup("Cascadet MMI")
                self.lum.select("MMI Out2")
                self.lum.addtogroup("Cascadet MMI")
                self.lum.select("MMI In_WG")
                self.lum.addtogroup("Cascadet MMI")
                self.lum.select('MMI Out_WG_Top')
                self.lum.addtogroup("Cascadet MMI")
                self.lum.select('MMI Out_WG_Bot')
                self.lum.addtogroup("Cascadet MMI")
                self.lum.select('MMIOut1 In_WG')
                self.lum.addtogroup("Cascadet MMI")
                self.lum.select('MMIOut1 Out_WG_Top')
                self.lum.addtogroup("Cascadet MMI")
                self.lum.select('MMIOut1 Out_WG_Bot')
                self.lum.addtogroup("Cascadet MMI")
                self.lum.select('MMIOut2 In_WG')
                self.lum.addtogroup("Cascadet MMI")
                self.lum.select('MMIOut2 Out_WG_Top')
                self.lum.addtogroup("Cascadet MMI")
                self.lum.select('MMIOut2 Out_WG_Bot')
                self.lum.addtogroup("Cascadet MMI")
                self.lum.select("S-Bend Top")
                self.lum.addtogroup("Cascadet MMI")
                self.lum.select("S-Bend Bot")
                self.lum.addtogroup("Cascadet MMI")
                
            elif Taper == True:
                if offset_Taper > OffMax:
                    self.lum.deleteall()
                    raise ValueError('You are Trying to move the Taper outside the MMI. This is not possible!')
                else:
                    # Delate the Structure to start new
                    self.lum.deleteall()
                    # Device specifications
              
                
                    Device_Length = 2*MMI_Length + 2 * WG_Length + 2*TaperLength + 2*SpaceX 
                    Device_Width = (2*MMI_Width + 2 * WaveLength) + 3*SpaceY  # MMI_Width

                    # creating the substrate
                    max_subH = Substrate_Height
                    min_subH = -Substrate_Height
                    min_subL = -3*(Device_Length / 2)
                    max_subL = 3*(Device_Length / 2)

                    self.lum.addrect()
                    self.lum.set("name", "Substrate_Main")
                    self.lum.set("y", 0)
                    self.lum.set("y span", Device_Width)
                    self.lum.set("z min", min_subH)
                    self.lum.set("z max", max_subH)
                    self.lum.set("x min", min_subL)
                    self.lum.set("x max", max_subL)
                    self.lum.set("material", MaterialSub)
                    
                    # creating the MMI
                    max_MMIH = WG_Height
                    max_MMIL = MMI_Length / 2
                    min_MMIL = -MMI_Length / 2
                    
                    if Slab_Height == 0:
                        z_Offset = max_subH + max_MMIH / 2
                    else:
                        # creating the thin film
                        min_slabH = max_subH
                        max_slabH = max_subH + Slab_Height

                        self.lum.addrect()
                        self.lum.set("name", "Slab Mian")
                        self.lum.set("y", 0)
                        self.lum.set("y span", Device_Width)
                        self.lum.set("z min", min_slabH)
                        self.lum.set("z max", max_slabH)
                        self.lum.set("x min", min_subL)
                        self.lum.set("x max", max_subL)
                        self.lum.set("material", MaterialSlab)

                        z_Offset = max_slabH + max_MMIH / 2
                        self.lum.select("MMI In")
                        self.lum.addtogroup("Cascadet MMI")
                        

                    # Triangle EQ for MMI Width
                    x = abs(max_MMIH / (np.cos((angle) * np.pi / 180)))  # in Radians
                    extention = np.sqrt(x ** 2 - max_MMIH ** 2)
                    MMI_Wid = MMI_Width + 2 * extention
                    
                    
                    Parameters__Position_X =  [0, -SpaceX -MMI_Length -2*TaperLength - 2*WG_Length , -SpaceX -MMI_Length -2*TaperLength - 2*WG_Length]
                    Parameters__Position_Y = [0, (MMI_Wid + SpaceY)/2,  -(MMI_Wid + SpaceY)/2]
                    
                    names_MMI = ["MMI In", "MMI Out1", "MMI Out2"]
                    for i in range(3):
                        self.lum.addwaveguide()
                        self.lum.set("name", names_MMI[i])
                        self.lum.set("x", Parameters__Position_X[i])
                        self.lum.set("y", Parameters__Position_Y[i])
                        self.lum.set("z", z_Offset)
                        self.lum.set("base height", max_MMIH)
                        self.lum.set("base angle", 90 - angle)
                        pole = np.array([[max_MMIL, 0], [min_MMIL, 0]])
                        self.lum.set("poles", pole)
                        self.lum.set("material", MaterialWG)
                        self.lum.set("base width", MMI_Wid)


                    # New x Length of the Tapers
                    maxLength = max_MMIL + TaperLength
                    minLength = min_MMIL - TaperLength

                    # Mirror the In and Out WG on both sides
                    maxWGL = [WG_Length, 0, 0]
                    minWGL = [0, -WG_Length, -WG_Length]
                    xPos = [maxLength, minLength, minLength]
                    yPos = [0 + OffsetInput, WG_Width / 2 + posOffset / 2, - WG_Width / 2 - posOffset / 2]


                    # Taper loop
                    # Taper Widths on Bott Cal
                    x = abs(max_MMIH / (np.cos((angle) * np.pi / 180)))  # in Radians
                    extention = np.sqrt(x ** 2 - max_MMIH ** 2)
                    TaperSideWidth = TaperWidth + 2 * extention

                    TaperPosXmin = [max_MMIL, min_MMIL, min_MMIL]
                    TaperPosXmax = [max_MMIL + TaperLength, min_MMIL - TaperLength, min_MMIL - TaperLength]

                    PosOffset = [0, (posOffset / 2 + WG_Width / 2), -(posOffset / 2 + WG_Width / 2)]
                    TaperPosYMax_BotR = [(0 + OffsetInput + WG_W / 2), (-WG_W / 2), (-WG_W / 2)]
                    TaperPosYMin_BotR = [(0 + OffsetInput - WG_W / 2), (WG_W / 2), (WG_W / 2)]
                    TaperPosYMax_TopR = [(0 + OffsetInput + WG_Width / 2), (-WG_Width / 2), (-WG_Width / 2)]
                    TaperPosYMin_TopR = [(0 + OffsetInput - WG_Width / 2), (+WG_Width / 2), (+WG_Width / 2)]

                    TaperPosYMax_BotL = [(0 + OffsetInput + TaperSideWidth / 2), (-TaperSideWidth / 2),
                                         (-TaperSideWidth / 2)]
                    TaperPosYMin_BotL = [(0 + OffsetInput - TaperSideWidth / 2), (+TaperSideWidth / 2),
                                         (+TaperSideWidth / 2)]
                    TaperPosYMax_TopL = [(0 + OffsetInput + TaperWidth / 2), (-TaperWidth / 2), (-TaperWidth / 2)]
                    TaperPosYMin_TopL = [(0 + OffsetInput - TaperWidth / 2), (+TaperWidth / 2), (+TaperWidth / 2)]
                    
                    
                    # Create two Bends Cos or Bezier depending on the poles option
                    x_span = SpaceX
                    y_span = (Parameters__Position_Y[1] + yPos[0]) - (Parameters__Position_Y[0] + yPos[1] ) 
                    
                    
                    if polesList == False:
                        pole = np.array([[0, 0], [x_span / 2, 0], [x_span / 2, y_span], [x_span, y_span]])
                        pole1 = np.array([[0, y_span], [x_span / 2, y_span], [x_span / 2, 0], [x_span,0]])
                    elif polesList == True:
                        K = ((np.pi -2)/np.pi)*100
                        pole = np.array([[0, 0], [0+(x_span*(K/100)), 0], [(x_span)-(x_span*(K/100)), y_span], [x_span, y_span]])
                        pole1 = np.array([[0, y_span], [0+(x_span*(K/100)), y_span], [(x_span)-(x_span*(K/100)), 0], [x_span, 0]])
                    else:
                        raise ValueError('Parameters["poles"] should be an boolen variable! Please set Parameters["poles"] to True if Bezier Curves needed. Set Parameters["poles"] to False for Cosinus Curve!')
        
        
                    names = ["S-Bend Top", "S-Bend Bot"]
                    S_Bend_X_Offset = [(-SpaceX - MMI_Length/2 - TaperLength - WG_Length), (SpaceX/2 - SpaceX/2 - MMI_Length/2 - TaperLength - WG_Length)]
                    S_Bend_Y_Offset = [( WG_Width / 2 + posOffset / 2) ,-( WG_Width / 2 + posOffset / 2)]
                    z_rotation = [0,180]
                    poles = [pole1, pole]
                    for i in range(2):
                        self.lum.addwaveguide()
                        self.lum.set("name", names[i])
                        self.lum.set("x", S_Bend_X_Offset[i])
                        self.lum.set("y", S_Bend_Y_Offset[i])
                        self.lum.set("z", z_Offset)
                        self.lum.set("base width", WG_W)
                        self.lum.set("base height", WG_Height)
                        self.lum.set("base angle", 90 - angle)
                        self.lum.set("poles", poles[i])
                        self.lum.set("material", MaterialWG)
                        self.lum.set("first axis","z")
                        self.lum.set("rotation 1",z_rotation[i])
                    
                    TapersNames1 = ['MMI Taper In_WG', 'MMI Taper Out WG_Top', 'MMI Taper Out WG_Bot']
                    TapersNames2 = ['MMIOut1 Taper In_WG', 'MMIOut1 Taper Out WG_Top', 'MMIOut1 Taper Out WG_Bot']
                    TapersNames3 = ['MMIOut2 Taper In_WG', 'MMIOut2 Taper Out WG_Top', 'MMIOut2 Taper ut WG_Bot']
                    TapersNames = [TapersNames1, TapersNames2, TapersNames3]
                    
                    
                    for j in range(3):
                        for i in range(len(xPos)):
                            TaperZmin = max_slabH
                            TaperZmax = max_slabH + max_MMIH
    
                            TaperXmin = TaperPosXmin[i]
                            TaperXmax = TaperPosXmax[i]
    
                            ymin_bot_l = TaperPosYMin_BotL[i]
                            ymax_bot_l = TaperPosYMax_BotL[i]
    
                            ymin_bot_r = TaperPosYMin_BotR[i]
                            ymax_bot_r = TaperPosYMax_BotR[i]
    
                            ymin_top_l = TaperPosYMin_TopL[i]
                            ymax_top_l = TaperPosYMax_TopL[i]
    
                            ymin_top_r = TaperPosYMin_TopR[i]
                            ymax_top_r = TaperPosYMax_TopR[i]
    
                            vtx = np.array([[TaperXmin, ymin_bot_l, TaperZmin],  # 1
                                            [TaperXmax, ymin_bot_r, TaperZmin],  # 2
                                            [TaperXmax, ymax_bot_r, TaperZmin],  # 3
                                            [TaperXmin, ymax_bot_l, TaperZmin],  # 4
                                            [TaperXmin, ymin_top_l, TaperZmax],  # 5
                                            [TaperXmax, ymin_top_r, TaperZmax],  # 6
                                            [TaperXmax, ymax_top_r, TaperZmax],  # 7
                                            [TaperXmin, ymax_top_l, TaperZmax],  # 8
                                            ])
                            a = [[np.array([[1, 4, 3, 2]], dtype=object)], [np.array([[1, 5, 8, 4]], dtype=object)],
                                 [np.array([[1, 2, 6, 5]], dtype=object)], [np.array([[2, 6, 7, 3]], dtype=object)],
                                 [np.array([[3, 4, 8, 7]], dtype=object)], [np.array([[5, 6, 7, 8]], dtype=object)]]
                            
                        
    
                            # Send Values to Lumerical and create solid
                            self.lum.putv('vertices', vtx)
                            self.lum.putv('facets', a)
                            self.lum.addplanarsolid(vtx, a)
                            self.lum.set('material', MaterialWG)
                            self.lum.set('name', TapersNames[j][i])
                            self.lum.set('y', Parameters__Position_Y[j] + PosOffset[i])
                            self.lum.set('x', Parameters__Position_X[j])
                            
                    # Names of the WGs
                    # Names of the WGs
                    names1 = ['MMI In_WG', 'MMI Out WG_Top', 'MMI Out WG_Bot']
                    names2 = ['MMIOut1 In_WG', 'MMIOut1 Out WG_Top', 'MMIOut1 Out WG_Bot']
                    names3 = ['MMIOut2 In_WG', 'MMIOut2 Out WG_Top', 'MMIOut2 In WG_Bot']
                    names = [names1, names2, names3]
                    


                    # create loop
                    for j in range(3):
                        for i in range(len(xPos)):
                            self.lum.addwaveguide()
                            self.lum.set("name", names[j][i])
                            self.lum.set("x", Parameters__Position_X[j] + xPos[i])
                            self.lum.set("y", Parameters__Position_Y[j] + yPos[i])
                            self.lum.set("z", z_Offset)
                            self.lum.set("base width", WG_Width_top)
                            self.lum.set("base height", max_MMIH)
                            self.lum.set("base angle", 90 - angle)
                            pole = np.array([[maxWGL[i], 0], [minWGL[i], 0]])
                            self.lum.set("poles", pole)
                            self.lum.set("material", MaterialSlab)
                

                self.lum.select("MMI In")
                self.lum.addtogroup("Cascadet MMI")
                self.lum.select("MMI Out1")
                self.lum.addtogroup("Cascadet MMI")
                self.lum.select("MMI Out2")
                self.lum.addtogroup("Cascadet MMI")
                
                self.lum.select("MMI In_WG")
                self.lum.addtogroup("Cascadet MMI")
                self.lum.select('MMI Out WG_Top')
                self.lum.addtogroup("Cascadet MMI")
                self.lum.select('MMI Out WG_Bot')
                self.lum.addtogroup("Cascadet MMI")
                
                self.lum.select('MMIOut1 In_WG')
                self.lum.addtogroup("Cascadet MMI")
                self.lum.select('MMIOut1 Out WG_Top')
                self.lum.addtogroup("Cascadet MMI")
                self.lum.select('MMIOut1 Out WG_Bot')
                self.lum.addtogroup("Cascadet MMI")
                
                self.lum.select('MMIOut2 In_WG')
                self.lum.addtogroup("Cascadet MMI")
                self.lum.select('MMIOut2 Out WG_Top')
                self.lum.addtogroup("Cascadet MMI")
                self.lum.select('MMIOut2 Out WG_Bot')
                self.lum.addtogroup("Cascadet MMI")
                
                self.lum.select("S-Bend Top")
                self.lum.addtogroup("Cascadet MMI")
                self.lum.select("S-Bend Bot")
                self.lum.addtogroup("Cascadet MMI")

                self.lum.select('MMI Taper In_WG')
                self.lum.addtogroup("Cascadet MMI")
                self.lum.select('MMI Taper Out WG_Top')
                self.lum.addtogroup("Cascadet MMI")
                self.lum.select('MMI Taper Out WG_Bot')
                self.lum.addtogroup("Cascadet MMI")
                
                self.lum.select('MMIOut1 Taper In_WG')
                self.lum.addtogroup("Cascadet MMI")
                self.lum.select('MMIOut1 Taper Out WG_Top')
                self.lum.addtogroup("Cascadet MMI")
                self.lum.select('MMIOut1 Taper Out WG_Bot')
                self.lum.addtogroup("Cascadet MMI")
                
                self.lum.select('MMIOut2 Taper In_WG')
                self.lum.addtogroup("Cascadet MMI")
                self.lum.select('MMIOut2 Taper Out WG_Top')
                self.lum.addtogroup("Cascadet MMI")
                self.lum.select('MMIOut2 Taper ut WG_Bot')
                self.lum.addtogroup("Cascadet MMI")
                    

            else:
                raise ValueError(
                    "Incorect Taper input. Taper must be an boolen. You can choose from Taper = True or Taper = False!")
        
            
            
            # create_cover
            self.lum.addrect()
            self.lum.set("name", "cladding")
            self.lum.set("material", MaterialClad)
            self.lum.set("y", 0)
            self.lum.set("y span", Device_Width)
            self.lum.set("z", z_Offset)
            self.lum.set("z span", max_MMIH * 2)
            self.lum.set("x min", min_subL)
            self.lum.set("x max", max_subL)
            self.lum.set("override mesh order from material database", True)
            self.lum.set("mesh order", 4)
            self.lum.set("alpha", 0.7)




    def GratingCoupler(self, Parameters):
        '''
       
        Parameters
        ----------
        Parameters : Dictionary
            Dictionary with all the data needet for the GratingCoupler.
            Parameters['Material GC'] : list of str
                List of Materials. the list should be with names (str) of a valid Lumerical materials.
                Check the names in Lumerical Materials viewer. For Example "Parameters['Material GC'] = ["SiO2 (Glass) - Palik", "Si (Silicon) - Palik"]"
            Parameters['Substrate Height'] : int/float
                Substrate height.
            Parameters["Length GC"]: int/float
                Lenght of the Grating Coupler Area
            Parameters["Width GC"]: int/float
                Widht of the Grating Coupler Area
            Parameters["Hight GC"]: int/float
                Hight of the Grating Coupler Material
            Parameters["Etch Depth GC"]: int/float
                How deep, taken from the Parameters["Hight GC"] will be the etchin depth of the gratings
            Parameters["Duty Cycle"]: int/float
                Duty cycle of the gratings. For example Parameters["Duty Cycle"] = 0.39 will result in 39% Duty Cycle 
            Parameters["Pitch GC"]: int/float
                Pitch of the Grating Coupler. For Example Parameters["Pitch GC"] = 0.6e-6 will result in 6um Etch Space + Rib Space = 0.6um.
            Parameters["Input Length GC"]: int/float
                An squere Waveguide with the same WG Height as the Grating coupler place before the Grating Coupler region will start. 
            Parameters["Output Length GC"]: int/float
                An squere Waveguide with the same WG Height as the Grating coupler place after the Grating Coupler region to finish the structure.
            Parameters['angle'] : int/float
                Angle of the Waveguide Walls. it is calculated WG_angle = 90 - angle.
                For anfle = 90 we get a perfect rect! In this case is used only when Parameters["Taper"] = True
            Parameters["Taper"] : boolen
                You can create an input Taper to your Grating Coupler structure
            Parameters['WG Height'] : int/float
                Waveguide hight. Also the height of the MMI section
            Parameters['WG Width'] : int/float
                Waveguide width. Also in this function and ONLY in this function this will be the
                ibverse Taper width!!!
            Parameters['Slab Height'] : int/float
                Slab height
            Parameters['Taper Length'] : int/float
                Taper Length
            Parameters["SMF Core Diameter"] : int/float
                Single Mode Fiber core Diameter
            Parameters["SMF Cladding Diameter"] : int/float
                Single Mode Fiber Cladding Diameter
            Parameters["SMF Core Index"]
                Single Mode Fiber Core Index
            Parameters["SMF Cladding Index"]
                Single Mode Fiber Cladding Index
            Parameters["SMF Theta"]: int/float
                Tilting Angle of the Single Mode Fiber to the Grating Coupler. Normaly we choose Parameters["SMF Theta"] = 15
            Parameters["SMF Z Span"]: int/float
                Lenght/Span of the Single Mode Fiber
            

        Returns
        -------
        None.

        '''

        # simplify variable names by removing spaces
        TargetLength = Parameters["Length GC"]
        WidthGC = Parameters["Width GC"]
        TaperLength = Parameters['Taper Length']
        Hight = Parameters["Hight GC"]
        EtchDepth = Parameters["Etch Depth GC"]
        DutyCycle = Parameters["Duty Cycle"]
        Pitch = Parameters["Pitch GC"]
        InputLlength = Parameters["Input Length GC"]
        OutputLength = Parameters["Output Length GC"]
        Material = Parameters["Material GC"]
        CoreDiameter = Parameters["SMF Core Diameter"]
        CladdingDiameter = Parameters["SMF Cladding Diameter"]
        ZSpan = Parameters["SMF Z Span"]
        Theta = Parameters["SMF Theta"]
        CoreIndex = Parameters["SMF Core Index"]
        CladdingIndex = Parameters["SMF Cladding Index"]
        Taper = Parameters["Taper"]
        SubstrateThickness = Parameters['Substrate Height']

        WG_Height = Parameters['WG Height']
        WG_Width = Parameters['WG Width']
        angle = Parameters['angle']



        # Create the MMI L2 Tapers for the Trapezoid
        GCNames = "Grating Coupler"
        FiberName = "SMF"
        

        # Make the Grating Coupler
        import math
        n_periods = math.floor(TargetLength / Pitch)
        fill_width = Pitch * DutyCycle
        etch_width = Pitch * (1 - DutyCycle)
        L = n_periods * Pitch + etch_width
        spanX = OutputLength + TargetLength + InputLlength

        if EtchDepth > Hight:
            EtchDepth = Hight
        elif EtchDepth < Hight:
            self.lum.addrect()
            self.lum.set("name", "lower layer")
            self.lum.set("x min", 0)
            self.lum.set("x max", L)
            self.lum.set("z min", 0)
            self.lum.set("z max", Hight - EtchDepth)
            
        

        self.lum.addrect()
        self.lum.set("name", "input waveguide")
        self.lum.set("x min", -InputLlength)
        self.lum.set("x max", 0)
        self.lum.set("z min", 0)
        self.lum.set("z max", Hight)

        self.lum.addrect()
        self.lum.set("name", "output waveguide")
        self.lum.set("x min", L)
        self.lum.set("x max", L + OutputLength)
        self.lum.set("z min", 0)
        self.lum.set("z max", Hight)


        for i in range(1, n_periods+1):
            self.lum.addrect()
            self.lum.set("name", "post")
            self.lum.set("x min", Pitch * (i - 1) + etch_width)
            self.lum.set("x max", Pitch * i)
            self.lum.set("z min", Hight - EtchDepth)
            self.lum.set("z max", Hight)

        self.lum.selectall()
        self.lum.set("material", Material[0])
        self.lum.set("y", 0)
        self.lum.set("y span", WidthGC)



        self.lum.selectall()
        self.lum.addtogroup(GCNames)
        self.lum.select(GCNames)
        self.lum.set("x", -TargetLength/2)





        # Build SMF

        core_index = CoreIndex
        cladding_index = CladdingIndex
        core_radius = CoreDiameter / 2
        cladding_radius = CladdingDiameter/2
        theta_rad = Theta / (180 / np.pi)
        L = ZSpan / np.cos(theta_rad)



        # Check if Material Cable is given
        # if len(Material) == 2:
        CoreIndex = CoreIndex
        CladdingIndex = CladdingIndex

        self.lum.addcircle()
        self.lum.set("name", "core")
        self.lum.set("override mesh order from material database", 1)
        self.lum.set("mesh order", 4)
        self.lum.set("first axis", "y")
        self.lum.set("rotation 1", 10)
        self.lum.set("alpha", 1)
        self.lum.set("radius", core_radius)
        self.lum.set("index", CoreIndex)
        self.lum.set("x", 0)
        self.lum.set("y", 0)
        self.lum.set("z", 0)
        self.lum.set("z span", L)
        self.lum.set("first axis", "y")
        self.lum.set("rotation 1", Theta)

        self.lum.addcircle()
        self.lum.set("name", "cladding")
        self.lum.set("override mesh order from material database", 1)
        self.lum.set("mesh order", 5)
        self.lum.set("first axis", "y")
        self.lum.set("rotation 1", 10)
        self.lum.set("alpha", 0.35)
        self.lum.set("radius", cladding_radius)
        self.lum.set("index", CladdingIndex)
        self.lum.set("x", 0)
        self.lum.set("y", 0)
        self.lum.set("z", 0)
        self.lum.set("z span", L)
        self.lum.set("first axis", "y")
        self.lum.set("rotation 1", Theta)

        # else:
            # CoreMaterial = Material[2]
            # CladdingMaterial = Material[3]
            # 
            # self.lum.addcircle()
            # self.lum.set("name", "core")
            # self.lum.set("override mesh order from material database", 1)
            # self.lum.set("mesh order", 4)
            # self.lum.set("first axis", "y")
            # self.lum.set("rotation 1", 10)
            # self.lum.set("alpha", 1)
            # self.lum.set("radius", core_radius)
            # self.lum.set("material", CoreMaterial)
            # self.lum.set("x", 0)
            # self.lum.set("y", 0)
            # self.lum.set("z", 0)
            # self.lum.set("z span", L)
            # self.lum.set("first axis", "y")
            # self.lum.set("rotation 1", Theta)
            # 
            # self.lum.addcircle()
            # self.lum.set("name", "cladding")
            # self.lum.set("override mesh order from material database", 1)
            # self.lum.set("mesh order", 5)
            # self.lum.set("first axis", "y")
            # self.lum.set("rotation 1", 10)
            # self.lum.set("alpha", 0.35)
            # self.lum.set("radius", cladding_radius)
            # self.lum.set("material", CladdingMaterial)
            # self.lum.set("x", 0)
            # self.lum.set("y", 0)
            # self.lum.set("z", 0)
            # self.lum.set("z span", L)
            # self.lum.set("first axis", "y")
            # self.lum.set("rotation 1", Theta)

        self.lum.select("core")
        self.lum.addtogroup('SMF')
        self.lum.select("cladding")
        self.lum.addtogroup('SMF')



        if Taper == True:
            # Add Taper
            TaperNames = "Taper"
            myscript = self.Script()
            spanX = OutputLength+TargetLength + InputLlength

            self.lum.addstructuregroup()
            self.lum.set("name", TaperNames)
            self.lum.set("construction group", 1)
            self.lum.adduserprop("thickness", 2, Hight)
            self.lum.adduserprop("angle_side", 0, angle)
            self.lum.adduserprop("width_l", 2, WG_Width)
            self.lum.adduserprop("width_r", 2, WidthGC)
            self.lum.adduserprop("hfrac_ref", 0, 1)
            self.lum.adduserprop("len", 2, TaperLength)
            self.lum.adduserprop("material", 5, Material[0])
            self.lum.adduserprop("index", 0, 1)
            self.lum.set("script", myscript)
            self.lum.set("x", -TargetLength/2 - InputLlength - TaperLength/2)
            self.lum.set("z", Hight/2)
            self.lum.set("y", 0)


            
        else:
            pass
            
        # global substrate and Si-Layer
        self.lum.addrect()
        self.lum.set("name", "SubstrateGlobal")
        self.lum.set("x", 0)
        self.lum.set("x span", 120e-6)
        self.lum.set("y", 0)
        self.lum.set("y span", 120e-6)
        self.lum.set("z", -SubstrateThickness/2)
        self.lum.set("z span", SubstrateThickness)
        self.lum.set("material", Material[1])
        
        self.lum.addrect()
        self.lum.set("name", "Si_Layer Global")
        self.lum.set("x", 0)
        self.lum.set("x span", 120e-6)
        self.lum.set("y", 0)
        self.lum.set("y span", 120e-6)
        self.lum.set("z", - SubstrateThickness -(2e-6) / 2)
        self.lum.set("z span", 2e-6)
        self.lum.set("material", Material[0])
        
        
        
        self.lum.addrect()
        self.lum.set("name", "Cladding Global")
        self.lum.set("x", 0)
        self.lum.set("x span", 120e-6)
        self.lum.set("y", 0)
        self.lum.set("y span", 120e-6)
        self.lum.set("z min", 0)
        self.lum.set("z max", Hight + 0.7e-6)
        self.lum.set("material", Material[1])
        self.lum.set("alpha", 0.7)
        self.lum.set("override mesh order from material database",1)
        self.lum.set("mesh order", 3)




    def RingGratingCoupler(self, Parameters):
    
        '''
       
        Parameters
        ----------
        Parameters : Dictionary
            Dictionary with all the data needet for the GratingCoupler.
            Parameters['Material GC'] : list of str
                List of Materials. the list should be with names (str) of a valid Lumerical materials.
                Check the names in Lumerical Materials viewer. For Example "Parameters['Material GC'] = ["SiO2 (Glass) - Palik", "Si (Silicon) - Palik"]"
            Parameters['Substrate Height'] : int/float
                Substrate height.
            Parameters["Length GC"]: int/float
                Lenght of the Grating Coupler Area
            Parameters["Width GC"]: int/float
                Widht of the Grating Coupler Area
            Parameters["Hight GC"]: int/float
                Hight of the Grating Coupler Material
            Parameters["GC Radius"]: int/float
                Radius of the Ring Grating Coupler in um. For Example "Parameters["GC Radius"] = 25e-6"
            Parameters["Etch Depth GC"]: int/float
                How deep, taken from the Parameters["Hight GC"] will be the etchin depth of the gratings
            Parameters["Duty Cycle"]: int/float
                Duty cycle of the gratings. For example Parameters["Duty Cycle"] = 0.39 will result in 39% Duty Cycle 
            Parameters["Pitch GC"]: int/float
                Pitch of the Grating Coupler. For Example Parameters["Pitch GC"] = 0.6e-6 will result in 6um Etch Space + Rib Space = 0.6um.
            Parameters["Input LengthGC"]: int/float
                An squere Waveguide with the same WG Height as the Grating coupler place before the Grating Coupler region will start. 
            Parameters["Output Length GC"]: int/float
                An squere Waveguide with the same WG Height as the Grating coupler place after the Grating Coupler region to finish the structure.
            Parameters['angle'] : int/float
                Angle of the Waveguide Walls. it is calculated WG_angle = 90 - angle.
                For anfle = 90 we get a perfect rect! 
            Parameters['WG Height'] : int/float
                Waveguide hight. Also the height of the MMI section
            Parameters['WG Width'] : int/float
                Waveguide width. 
            Parameters['Slab Height'] : int/float
                Slab height
            Parameters["SMF Core Diameter"] : int/float
                Single Mode Fiber core Diameter
            Parameters["SMF Cladding Diameter"] : int/float
                Single Mode Fiber Cladding Diameter
            Parameters["SMF Core Index"]
                Single Mode Fiber Core Index
            Parameters["SMF Cladding Index"]
                Single Mode Fiber Cladding Index
            Parameters["SMF Theta"]: int/float
                Tilting Angle of the Single Mode Fiber to the Grating Coupler. Normaly we choose Parameters["SMF Theta"] = 15
            Parameters["SMF Z Span"]: int/float
                Lenght/Span of the Single Mode Fiber
            

        Returns
        -------
        None.

        '''

        # simplify variable names by removing spaces
        TargetLength = Parameters["Length GC"]
        WidthGC = Parameters["Width GC"]
        TaperLength = Parameters['Taper Length']
        Hight = Parameters["Hight GC"]
        EtchDepth = Parameters["Etch Depth GC"]
        DutyCycle = Parameters["Duty Cycle"]
        Pitch = Parameters["Pitch GC"]
        InputLlength = Parameters["Input Length GC"]
        OutputLength = Parameters["Output Length GC"]
        Material = Parameters["Material GC"]
        CoreDiameter = Parameters["SMF Core Diameter"]
        CladdingDiameter = Parameters["SMF Cladding Diameter"]
        ZSpan = Parameters["SMF Z Span"]
        Theta = Parameters["SMF Theta"]
        CoreIndex = Parameters["SMF Core Index"]
        CladdingIndex = Parameters["SMF Cladding Index"]
        Taper = Parameters["Taper"]
        SubstrateThickness = Parameters['Substrate Height']
        GCRadius = Parameters["GC Radius"]

        
        WG_Height = Parameters['WG Height']
        WG_Width = Parameters['WG Width']
        angle = Parameters['angle']



        # Create the MMI L2 Tapers for the Trapezoid
        GCNames = "Grating Coupler"
        FiberName = "SMF"
    


        # Make the Grating Coupler
        import math
        n_periods = math.floor(TargetLength / Pitch)
        fill_width = Pitch * DutyCycle
        etch_width = Pitch * (1 - DutyCycle)
        L = n_periods * Pitch + etch_width
        spanX = OutputLength + TargetLength + InputLlength
        Theta = np.arcsin( 0.5*WidthGC/GCRadius ) * 180/np.pi
        

        if EtchDepth > Hight:
            EtchDepth = Hight
        elif EtchDepth < Hight:
            self.lum.addring()
            self.lum.set("name", "lower layer")
            self.lum.set("inner radius",GCRadius * np.cos(Theta*np.pi/180))
            self.lum.set("outer radius",GCRadius +L + OutputLength)
            self.lum.set("z min", 0)
            self.lum.set("z max", Hight - EtchDepth)
            self.lum.set("theta start",-Theta)
            self.lum.set("theta stop",Theta)
            


 
        # input section
        self.lum.addring()
        self.lum.set("name","input section")
        # self.lum.set("inner radius",GCRadius * np.cos(Theta*np.pi/180))
        self.lum.set("inner radius", 0)
        self.lum.set("outer radius",GCRadius)
        self.lum.set("x",0)
        self.lum.set("y",0)
        self.lum.set("z min",0)
        self.lum.set("z max",Hight)
        self.lum.set("theta start",-Theta)
        self.lum.set("theta stop",Theta)


        # output section
        self.lum.addring()
        self.lum.set("name","output section")
        self.lum.set("inner radius",GCRadius + L)
        self.lum.set("outer radius",GCRadius +L + OutputLength)
        self.lum.set("x",0)
        self.lum.set("y",0)
        self.lum.set("z min",0)
        self.lum.set("z max",Hight)
        self.lum.set("theta start",-Theta)
        self.lum.set("theta stop",Theta)

       
        # Add Grating
        for i in range(1, n_periods+1):
            self.lum.addring()
            self.lum.set("x",0)
            self.lum.set("y",0)
            self.lum.set("z min", Hight - EtchDepth)
            self.lum.set("z max",Hight)
            self.lum.set("theta start",-Theta)
            self.lum.set("theta stop",Theta)
            self.lum.set("inner radius",GCRadius +  Pitch * (i - 1) + etch_width)
            self.lum.set("outer radius",GCRadius + Pitch*i)
            self.lum.set("name", "post")
       

        self.lum.selectall()
        self.lum.set("material", Material[0])
 

        self.lum.selectall()
        self.lum.addtogroup(GCNames)
        self.lum.select(GCNames)
        self.lum.set("x", -TargetLength/2)





        # Build SMF

        core_index = CoreIndex
        cladding_index = CladdingIndex
        core_radius = CoreDiameter / 2
        cladding_radius = CladdingDiameter/2
        theta_rad = Theta / (180 / np.pi)
        L = ZSpan / np.cos(theta_rad)



        # Check if Material Cable is given
        # if len(Material) == 2:
        CoreIndex = CoreIndex
        CladdingIndex = CladdingIndex

        self.lum.addcircle()
        self.lum.set("name", "core")
        self.lum.set("override mesh order from material database", 1)
        self.lum.set("mesh order", 4)
        self.lum.set("first axis", "y")
        self.lum.set("rotation 1", 10)
        self.lum.set("alpha", 1)
        self.lum.set("radius", core_radius)
        self.lum.set("index", CoreIndex)
        self.lum.set("x", 0)
        self.lum.set("y", 0)
        self.lum.set("z", 0)
        self.lum.set("z span", L)



        self.lum.addcircle()
        self.lum.set("name", "cladding")
        self.lum.set("override mesh order from material database", 1)
        self.lum.set("mesh order", 5)
        self.lum.set("first axis", "y")
        self.lum.set("rotation 1", 10)
        self.lum.set("alpha", 0.35)
        self.lum.set("radius", cladding_radius)
        self.lum.set("index", CladdingIndex)
        self.lum.set("x", 0)
        self.lum.set("y", 0)
        self.lum.set("z", 0)
        self.lum.set("z span", L)


        # else:
        #     CoreMaterial = Material[2]
        #     CladdingMaterial = Material[3]
        # 
        #     self.lum.addcircle()
        #     self.lum.set("name", "core")
        #     self.lum.set("override mesh order from material database", 1)
        #     self.lum.set("mesh order", 4)
        #     self.lum.set("first axis", "y")
        #     self.lum.set("rotation 1", 10)
        #     self.lum.set("alpha", 1)
        #     self.lum.set("radius", core_radius)
        #     self.lum.set("material", CoreMaterial)
        #     self.lum.set("x", 0)
        #     self.lum.set("y", 0)
        #     self.lum.set("z", 0)
        #     self.lum.set("z span", L)
        # 
        # 
        #     self.lum.addcircle()
        #     self.lum.set("name", "cladding")
        #     self.lum.set("override mesh order from material database", 1)
        #     self.lum.set("mesh order", 5)
        #     self.lum.set("first axis", "y")
        #     self.lum.set("rotation 1", 10)
        #     self.lum.set("alpha", 0.35)
        #     self.lum.set("radius", cladding_radius)
        #     self.lum.set("material", CladdingMaterial)
        #     self.lum.set("x", 0)
        #     self.lum.set("y", 0)
        #     self.lum.set("z", 0)
        #     self.lum.set("z span", L)


        self.lum.select("core")
        self.lum.addtogroup('SMF')
        self.lum.select("cladding")
        self.lum.addtogroup('SMF')
        self.lum.select('SMF')
        self.lum.set("x", -TargetLength/2 +GCRadius + core_radius) #TargetLength

        
        
        
        # Triangle EQ for waveguide Width
        x = abs(Hight / (np.cos((angle) * np.pi / 180)))  # in Radians
        extention = np.sqrt(x ** 2 - Hight ** 2)
        WG_W = WG_Width + 2 * extention
    
      
        # Add Waveguide 
        names = ["Straight Waveguide"]
        self.lum.addwaveguide()
        self.lum.set("name", names[0])
        self.lum.set("x", -GCRadius/2)
        self.lum.set("y", 0)
        self.lum.set("z", Hight/2)
        self.lum.set("base width", WG_W)
        self.lum.set("base height", Hight)
        self.lum.set("base angle", 90 - angle)
        pole = np.array([[-20e-6, 0], [5e-6, 0]])
        self.lum.set("poles", pole)
        self.lum.set("material", Material[0])



        self.lum.select("Straight Waveguide")
        self.lum.addtogroup('Input Waveguide')
       

        # global substrate and Si-Layer
        self.lum.addrect()
        self.lum.set("name", "SubstrateGlobal")
        self.lum.set("x", 0)
        self.lum.set("x span", 120e-6)
        self.lum.set("y", 0)
        self.lum.set("y span", 120e-6)
        self.lum.set("z", -SubstrateThickness/2)
        self.lum.set("z span", SubstrateThickness)
        self.lum.set("material", Material[1])
        
        self.lum.addrect()
        self.lum.set("name", "Si_Layer Global")
        self.lum.set("x", 0)
        self.lum.set("x span", 120e-6)
        self.lum.set("y", 0)
        self.lum.set("y span", 120e-6)
        self.lum.set("z", - SubstrateThickness -(2e-6) / 2)
        self.lum.set("z span", 2e-6)
        self.lum.set("material", Material[0])
        
        
        
        self.lum.addrect()
        self.lum.set("name", "Cladding Global")
        self.lum.set("x", 0)
        self.lum.set("x span", 120e-6)
        self.lum.set("y", 0)
        self.lum.set("y span", 120e-6)
        self.lum.set("z min", 0)
        self.lum.set("z max", Hight + 0.7e-6)
        self.lum.set("material", Material[1])
        self.lum.set("alpha", 0.7)
        self.lum.set("override mesh order from material database",1)
        self.lum.set("mesh order", 3)
        




    
    
    
    def GratingCouplerNeff(self, Parameters):

        Hight = Parameters["Hight GC"]
        EtchDepth = Parameters["Etch Depth GC"]
        Pitch = Parameters["Pitch GC"]
        Material = Parameters["Material  GC"]



        # Make the Grating Coupler
        Gap = 1e-6 - Pitch


        # Create Ribs and Gaps
        # for i in range(1, 2 ):
        self.lum.addrect()
        self.lum.set("name", "Pitch")
        self.lum.set("x", 0)
        self.lum.set("x span", Pitch)
        self.lum.set("z min", Hight - EtchDepth)
        self.lum.set("z max", Hight)
        self.lum.set("material", Material[0])
        self.lum.set("y", 0)
        self.lum.set("y span", 5e-6)

        #Add Bottom Layer of Material
        self.lum.addrect()
        self.lum.set("name", "lower layer")
        self.lum.set("x", 0)
        self.lum.set("x span", Pitch + 2*Gap)
        self.lum.set("z min", 0)
        self.lum.set("z max", Hight - EtchDepth)
        self.lum.set("material", Material[0])
        self.lum.set("y", 0)
        self.lum.set("y span", 5e-6)

        # Add Substrate
        self.lum.addrect()
        self.lum.set("name", "Substrate")
        self.lum.set("x",  0)
        self.lum.set("x span", Pitch + 2*Gap)
        self.lum.set("z", -(1e-6) / 2)
        self.lum.set("z span", 1e-6)
        self.lum.set("y", 0)
        self.lum.set("y span", 5e-6)
        self.lum.set("material", Material[1])

        # Add Cladding
        self.lum.addrect()
        self.lum.set("name", "Cladding")
        self.lum.set("x", 0)
        self.lum.set("x span", Pitch + 2*Gap)
        self.lum.set("z min", Hight - EtchDepth)
        self.lum.set("z max", Hight + 0.7e-6)
        self.lum.set("y", 0)
        self.lum.set("y span", 5e-6)
        self.lum.set("material", Material[1])
        self.lum.set('override mesh order from material database', 1)
        self.lum.set("mesh order", 4)
        self.lum.set("alpha", 0.7)

        self.lum.selectall()
        self.lum.addtogroup("GC Rib")
        self.lum.select("GC Rib")
        self.lum.set("first axis", "z")
        self.lum.set("rotation 1", 90)



    def LenseModeSystem(self,Parameters):
        
        f_lense = Parameters["Focal Length"]
        Outter_R = Parameters["Ring Outter Radius"]
        Iner_R   = Parameters["Ring Iner Radius"]
        Lense_Thickness = Parameters["Lense Thickness"]
        x_res = Parameters['x res'] 
        Wavelength = Parameters['Wavelength'] 
        Port_Span = Parameters["Port Span"] 
        Materials = Parameters['Material']


        # SMF Parameters
        Core_Diameter = Parameters["SMF Core Diameter"] 
        Cladding_Diameter = Parameters["SMF Cladding Diameter"] 
        Core_Index = Parameters["SMF Core Index"]
        Cladding_Index = Parameters["SMF Cladding Index"] 
        
        
        self.lum.select("Straight Waveguide::Waveguide")
        z_Pos = self.lum.get("z")
        
        self.lum.addcircle()
        self.lum.set("name", "core")
        self.lum.set("override mesh order from material database", 1)
        self.lum.set("mesh order", 4)
        self.lum.set("alpha", 1)
        self.lum.set("radius", Core_Diameter / 2)
        self.lum.set("index", Core_Index)
        self.lum.set("x", 0)
        self.lum.set("y", 0)
        self.lum.set("z min", z_Pos + f_lense )
        self.lum.set("z max", z_Pos + f_lense + 2e-6)
        
        
        
        self.lum.addcircle()
        self.lum.set("name", "cladding")
        self.lum.set("override mesh order from material database", 1)
        self.lum.set("mesh order", 5)
        self.lum.set("alpha", 0.35)
        self.lum.set("radius", Cladding_Diameter / 2)
        self.lum.set("index", Cladding_Index)
        self.lum.set("x", 0)
        self.lum.set("y", 0)
        self.lum.set("z min", z_Pos + f_lense)
        self.lum.set("z max", z_Pos + f_lense + 2e-6)
        
        self.lum.select("core")
        self.lum.addtogroup('SMF')
        self.lum.select("cladding")
        self.lum.addtogroup('SMF')
        
        
        
        # Extract SMF Positions
        self.lum.select("SMF::core")
        x_Pos = self.lum.get("x")
        y_Pos = self.lum.get("y")
        self.lum.select("Straight Waveguide::Waveguide")
        z_Pos = self.lum.get("z")
        self.lum.select("Substrate")
        z_Sub = self.lum.get("z max")
        
        self.lum.addring()
        self.lum.set("name", "Support Ring Lense")
        self.lum.set("override mesh order from material database", 1)
        self.lum.set("mesh order", 4)
        self.lum.set("alpha", 0.7)
        self.lum.set("x", x_Pos)
        self.lum.set("y", y_Pos)
        self.lum.set("z min", z_Sub)
        self.lum.set("z max", z_Pos + f_lense)
        self.lum.set("outer radius", Outter_R)
        self.lum.set("inner radius", Iner_R)
        self.lum.set("theta start",0)
        self.lum.set("theta stop", 0)
        self.lum.set("material", Materials[2])
         
        self.lum.addsphere()
        self.lum.set("name", "Lense")
        self.lum.set("override mesh order from material database", 1)
        self.lum.set("mesh order", 4)
        self.lum.set("x", x_Pos)
        self.lum.set("y", y_Pos)
        self.lum.set("z", z_Pos + f_lense)
        self.lum.set("radius", Iner_R)
        self.lum.set("make ellipsoid", 1)
        self.lum.set("radius 2", Iner_R)
        self.lum.set("radius 3", Lense_Thickness)
        self.lum.set("material", Materials[2])
        
        self.lum.select("Support Ring Lense")
        self.lum.addtogroup('Funnel')
        self.lum.select("Lense")
        self.lum.addtogroup('Funnel')
        
        self.lum.select("Substrate")
        x_adj = self.lum.get("x min")
        self.lum.select("SMF")
        self.lum.set("x", (Parameters["SMF Core Diameter"]/2) -abs(x_adj))
        self.lum.select("Funnel")
        self.lum.set("x", (Parameters["SMF Core Diameter"]/2) - abs(x_adj))
        
        self.lum.select("SMF")
        x_Pos = self.lum.get("x")
        self.lum.select("SMF::core")
        y_Pos = self.lum.get("y")
        z_Pos = self.lum.get("z max")
        radius = self.lum.get("radius")
        
        self.lum.select("Substrate")
        Solver_X = self.lum.get("x")
        Solver_X_Span = self.lum.get("x span")
        Solver_X_max = self.lum.get("x max")
        Solver_X_min = self.lum.get("x min")
        Solver_Y = self.lum.get("y")
        Solver_Y_Span = self.lum.get("y span")
        Solver_Z = self.lum.get("z min")
        

        # Solver Object
        self.lum.addfdtd()
        # self.lum.set("x", Solver_X)
        # self.lum.set("x span", Solver_X_Span)
        self.lum.set("x max", Solver_X_max)
        self.lum.set("x min", Solver_X_min - radius)
        self.lum.set("y", Solver_Y)
        self.lum.set("y span", Solver_Y_Span)
        self.lum.set('simulation temperature', 273.15 + 20)
        self.lum.set("z min", Solver_Z)
        self.lum.set("z max", z_Pos + 0.5e-6)
        self.lum.set('z min bc', 'PML')
        self.lum.set('z max bc', 'PML')
        self.lum.set('mesh type', 'auto non-uniform')
        self.lum.set('min mesh step', x_res)
        self.lum.set('set simulation bandwidth', 0)
        self.lum.set('global source center wavelength', Wavelength)
        self.lum.set('global source wavelength span', 0)


        # Add Gausssian Source
        self.lum.addgaussian()
        self.lum.set("injection axis", "z-axis")
        self.lum.set("direction","Backward")
        self.lum.set("x", x_Pos)
        self.lum.set("x span", radius * 4)
        self.lum.set("y", y_Pos)
        self.lum.set("y span", radius * 4)
        self.lum.set("z", z_Pos)
        self.lum.set("waist radius w0", radius)
        
        
        self.lum.select("Straight Waveguide::Waveguide")
        WG_Z = self.lum.get("z")
        # Add Movie monitor
        self.lum.addmovie()
        self.lum.set("name", "Movie")
        self.lum.set("y", Solver_Y)
        self.lum.set("y span", Solver_Y_Span)
        self.lum.set("z", WG_Z)
        self.lum.set("x", Solver_X)
        self.lum.set("x span", Solver_X_Span)

        # Add Power and Freq Monitor
        self.lum.addpower()
        self.lum.set('monitor type', '2D Z-normal')
        self.lum.set("y", Solver_Y)
        self.lum.set("y span", Solver_Y_Span)
        self.lum.set("x", Solver_X)
        self.lum.set("x span", Solver_X_Span)
        self.lum.set("z", WG_Z)
        self.lum.set('output Px', 1)
        self.lum.set('output Py', 1)
        self.lum.set('output Pz', 1)
        self.lum.set('output power', 1)
        
        names = ["Input", "Output"]
        self.lum.select("Straight Waveguide::Waveguide")
        WG_X = self.lum.get("poles")
        x_Monitor = [WG_X[0][0] + 0.1e-6, WG_X[1][0] - 0.1e-6]
        
        self.lum.select("Substrate")
        z_sub = self.lum.get("z max")
        self.lum.select("SMF::core")
        z_Min = self.lum.get("z min")
        

        
        for i in range(len(names)):
            self.lum.addpower()
            self.lum.set('name', names[i])
            self.lum.set('monitor type', '2D X-normal')
            self.lum.set("x", x_Monitor[i] )
            self.lum.set("y", Solver_Y)
            self.lum.set("y span", Port_Span[1])
            self.lum.set("z", WG_Z)
            self.lum.set("z span", Port_Span[2])
            self.lum.set('output Px', 1)
            self.lum.set('output Py', 1)
            self.lum.set('output Pz', 1)
            
            
        for i in range(3):
            self.lum.addpower()
            self.lum.set('name', "Power Beam " + str(i))
            self.lum.set('monitor type', '2D Z-normal')
            self.lum.set("y", Solver_Y)
            self.lum.set("y span", Solver_Y_Span)
            self.lum.set("x", x_Pos)
            self.lum.set("x span", Solver_X_Span)
            self.lum.set("z", z_Min)
            self.lum.set('output Px', 1)
            self.lum.set('output Py', 1)
            self.lum.set('output Pz', 1)
            self.lum.set('output power', 1)
            z_Min = z_Min/2
        
        # Monitor Behind the Port 
        self.lum.addpower()
        self.lum.set('name', "Power Beam back")
        self.lum.set('monitor type', '2D Z-normal')
        self.lum.set("y", Solver_Y)
        self.lum.set("y span", Solver_Y_Span)
        self.lum.set("x", x_Pos)
        self.lum.set("x span", Solver_X_Span)
        self.lum.set("z", z_Min)
        self.lum.set('output Px', 1)
        self.lum.set('output Py', 1)
        self.lum.set('output Pz', 1)
        self.lum.set('output power', 1)
              

        # Monitor Mode  
        self.lum.addprofile()
        self.lum.set('name', "Beam Profile Monitor")
        self.lum.set('monitor type', '2D Y-normal')
        self.lum.set("y", 0)
        self.lum.set("x", x_Pos)
        self.lum.set("x span", Parameters["SMF Core Diameter"])
        self.lum.set("z min", 0)
        self.lum.set("z max", z_Pos + 0.5e-6)
        self.lum.set('output Px', 1)
        self.lum.set('output Py', 1)
        self.lum.set('output Pz', 1)
        self.lum.set('output power', 1)
        
        # Monitor Index  
        self.lum.addindex()
        self.lum.set('name', "index_monitor")
        self.lum.set('monitor type', '3D')
        self.lum.set("x max", Solver_X_max)
        self.lum.set("x min", Solver_X_min - radius)
        self.lum.set("y", Solver_Y)
        self.lum.set("y span", Solver_Y_Span)
        self.lum.set("z min", Solver_Z)
        self.lum.set("z max", z_Pos + 0.5e-6)
  
  
  
  
    def LenseFiber(self, Parameters):
        
        
        '''
       
        Parameters
        ----------
        Parameters : Dictionary
            Dictionary with all the data needet for the GratingCoupler.
            Take care since this object will create his own FDTD simulation object. No extra Solvers needed!
            Parameters['Material'] : list of str
                List of Materials. the list should be with names (str) of a valid Lumerical materials.
                Check the names in Lumerical Materials viewer. For Example "Parameters['Material'] = ["SiO2 (Glass) - Palik", "Si (Silicon) - Palik"]"
            Parameters["Lense Diameter"] : int/float
                Lense diameter.
            Parameters["Lense Thickness"]: int/float
                Lens thickness
            Parameters["SMF Core Diameter"] : int/float
                Single Mode Fiber core Diameter
            Parameters["SMF Cladding Diameter"] : int/float
                Single Mode Fiber Cladding Diameter
            Parameters["SMF Core Index"]
                Single Mode Fiber Core Index
            Parameters["SMF Cladding Index"]
                Single Mode Fiber Cladding Index
            Parameters['Support Cylunder Hight']: int/float
                Hight of the Cylinder holding the Lense. This is needed so that the Lense can have flat fine serface to be printed on.
            Parameters['Hexagon Hight']: optional int/float
                Hight of the Hexagon holding the Cylinder. If the Parameter is not set the Hexagon will be made with 4 um thickness. Where 2 um 
                will go into the fiber array and 2 um will be above to print the cylinder on it. This is how is done in NanoPrintX.
            Parameters['x res'] : int/float
                Mesh step size
            Parameters['Wavelength']: int/float
                Simulation Wavelenght
                

        Returns
        -------
        None.

        '''
    
    
    
        Lens_d   = Parameters["Lense Diameter"]
        Lense_Thickness = Parameters["Lense Thickness"]
        Materials = Parameters['Material']
        Cylinder_Hight = Parameters["Support Cylunder Hight"]
        
        if "Hexagon Hight" in list(Parameters.keys()):
            Hexagon_Hight = Parameters['Hexagon Hight']
        else:
            Hexagon_Hight = 4e-6
            
        


        # SMF Parameters
        Core_Diameter = Parameters["SMF Core Diameter"] 
        Cladding_Diameter = Parameters["SMF Cladding Diameter"] 
        Core_Index = Parameters["SMF Core Index"]
        Cladding_Index = Parameters["SMF Cladding Index"] 
        
       
        
        
        #Create SMF
        self.lum.addcircle()
        self.lum.set("name", "core")
        self.lum.set("override mesh order from material database", 0)
        self.lum.set("alpha", 1)
        self.lum.set("radius", Core_Diameter / 2)
        self.lum.set("index", Core_Index)
        self.lum.set("x", 0)
        self.lum.set("y", 0)
        self.lum.set("z", 0)
        self.lum.set("z span", 4e-6)
        
        
        self.lum.addcircle()
        self.lum.set("name", "cladding")
        self.lum.set("override mesh order from material database", 1)
        self.lum.set("mesh order", 3)
        self.lum.set("alpha", 0.35)
        self.lum.set("radius", Cladding_Diameter / 2)
        self.lum.set("index", Cladding_Index)
        self.lum.set("x", 0)
        self.lum.set("y", 0)
        self.lum.set("z", 0)
        self.lum.set("z span", 4e-6)
        
        self.lum.select("core")
        self.lum.addtogroup('SMF')
        self.lum.select("cladding")
        self.lum.addtogroup('SMF')
        
        
        
        # Extract SMF Positions
        self.lum.select("SMF::core")
        x_Pos = self.lum.get("x")
        y_Pos = self.lum.get("y")
        z_Pos = self.lum.get("z")
        
        
        
        # Create Hexagon support
        self.lum.addpoly()
        self.lum.set("name", "Support Hexagon")
        self.lum.set("override mesh order from material database", 1)
        self.lum.set("mesh order", 4)
        pole = np.array([[0, -Cladding_Diameter / 2], [np.sin(60*np.pi/180)*Cladding_Diameter / 2, -Cladding_Diameter / 4], [np.sin(60*np.pi/180)*Cladding_Diameter / 2, Cladding_Diameter / 4], [0,Cladding_Diameter / 2], [-np.sin(60*np.pi/180)*Cladding_Diameter / 2, Cladding_Diameter / 4], [-np.sin(60*np.pi/180)*Cladding_Diameter / 2, -Cladding_Diameter / 4]])
        self.lum.set("vertices", pole)
        self.lum.set("x", 0)
        self.lum.set("y", 0)
        self.lum.set("z", z_Pos + Hexagon_Hight/2)
        self.lum.set("z span", Hexagon_Hight)
        self.lum.set("material", Materials[1])
        
        
       
        # Set support ring with 0,5 um min thickness becouse of Voxel of Nanoscribe
        z_Pos = self.lum.get("z max")
        self.lum.addcircle()
        self.lum.set("name", "Support Ring Lense")
        self.lum.set("override mesh order from material database", 1)
        self.lum.set("mesh order", 4)
        self.lum.set("radius", Core_Diameter / 2)
        self.lum.set("material", Materials[1])
        self.lum.set("x", 0)
        self.lum.set("y", 0)
        self.lum.set("z", z_Pos + Cylinder_Hight/2)
        self.lum.set("z span", Cylinder_Hight)
        
        
        # Add lense to on top of the o,5um Voxel ring
        z_Pos = self.lum.get("z max")
        self.lum.addsphere()
        self.lum.set("name", "Lense")
        self.lum.set("override mesh order from material database", 1)
        self.lum.set("mesh order", 4)
        self.lum.set("x", 0)
        self.lum.set("y", 0)
        self.lum.set("z", z_Pos )
        self.lum.set("radius", Lens_d/2)
        self.lum.set("make ellipsoid", 1)
        self.lum.set("radius 2", Lens_d/2)
        self.lum.set("radius 3", Lense_Thickness)
        self.lum.set("material", Materials[1])
        

        
    


    # =============================================================================
    # Functions for the solvers
    # =============================================================================


    def setStraightWaveguideFDTDSolver(self, Parameters):
        '''


        Parameters
        ----------
        Parameters : Dictionary
            Dictionary with all the data needet for the StraightWaveguideFDTDSolver.
            Parameters : Dictionary
                Dictionary with all the data needet for the Bend Wavaguide. Data needet:
            Parameters['Substrate Height'] : int/float
               Substrate height.
            Parameters['WG Length'] : int/float
               Waveguide Length
            Parameters['WG Height'] : int/float
               Waveguide hight. Also the height of the MMI section
            Parameters['WG Width'] : int/float
               Waveguide width.
            Parameters["Taper"] : boolen
               If Taper == False, only straight Waveguide will be simulated,
               If Taper == True an Taper will be simulated
            Parameters['Taper Width'] : int/float
               Taper backside Width. Taper Fronside width is the width of the Waveguide
            Parameters['Taper Length'] : int/float
               Taper Length
            Parameters['x res'] : int/float
                 EME Mesh resolutio,
            Parameters['Slab Height'] : int/float
               Height of the slab.
            Parameters['Wavelength'] : int/float
               Wavelength
            Parameters["Waveguide Angle"] : int/float
               This Parameter will set the theta ratation angle of the port. It can be 90 or 180.
            Parameters["Port Span"] : list of floats/ints
               List of x,y and z span of the Ports. For this simulation only y and z parametes will be taken.
            Parameters["Mode"] : str
                Mode to choose from ("fundamental TE mode", "fundamental TM mode", "fundamental mode")
            Parameters["Taper Type"] : anything, optional
                    This function will check if you have set Parameters["Taper Type"] to anaything, for example "Parameters["Taper Type"]=1" 
                    and if so it will design an Inverse Taper Structure with no Cladding. Here the option "Cladding" is not active and will be ignored.
                    If the user didnt give the "Taper Type" as dictionary key, then an normal taper structure will be simulated.
                    
                    If Parameters["Taper Type"] is given, themn the user need to set couple more parameters:
                        Parameters['PWB Taper Width Back'] : int/float
                            Photonic Wirebonding (PWB) Width back side (to the Photonic Wire Bonding)
                        Parameters['PWB Taper Hight Back'] : int/float
                            Photonic Wire Bonding Height back side (to the Photonic Wire Bonding)
                        Parameters['PWB Taper Width Front'] : int/float
                            Photonic Wirebonding (PWB) Width front side (to the photonic waveguide)
                        Parameters['PWB Taper Hight Front'] : int/float
                            Photonic Wire Bonding Height front side (to the photonic waveguide)
                        Parameters['PWB Taper Length'] : int/float
                            Length of the Photonic Wire Bonding Taper
                    

        Returns
        -------
        None.

        '''


        Substrate_Height = Parameters['Substrate Height']
        WG_Length = Parameters['WG Length']
        WG_Height = Parameters['WG Height']
        WG_Width = Parameters['WG Width']
        angle = Parameters["Waveguide Angle"]
        x_res = Parameters['x res']
        Slab_Height = Parameters['Slab Height']
        WaveLength = Parameters['Wavelength']
        Mode = Parameters["Mode"]
        y_Port_Span = Parameters["Port Span"][1]
        z_Port_Span = Parameters["Port Span"][2]
        Taper = Parameters['Taper']
        TaperWidth = Parameters['Taper Width']
        TaperLength = Parameters['Taper Length']

        if "Taper Type" in list(Parameters.keys()):
                TaperType = "Inverse"
                TaperWidthF = Parameters['PWB Taper Width Front']
                TaperWidthB = Parameters['PWB Taper Width Back']
                TaperHightB = Parameters['PWB Taper Hight Back']
                TaperHightF = Parameters['PWB Taper Hight Front']
                TaperLength_PWB = Parameters['PWB Taper Length']
        else:
            TaperType = "Normal"
            
                
        # Device specifications
        Device_Width = 2*WG_Length + WaveLength * 2  # MMI_Width

        max_slabH = Slab_Height
        MonitorHeight = Substrate_Height + max_slabH + WG_Height / 2
        EME_WGLength = (WG_Length * np.cos(angle * np.pi / 180))




        if Taper == False:
            if angle == 0:

                # Adds a FDTD Solver
                self.lum.addfdtd()
                self.lum.set("x", WG_Length / 2)
                self.lum.set("x span", WG_Length)
                self.lum.set("y", 0)
                self.lum.set("y span", Device_Width)
                self.lum.set('simulation temperature', 273.15 + 20)
                self.lum.set("z", Substrate_Height)
                self.lum.set("z span", 4e-6)
                self.lum.set('x min bc', 'PML')
                self.lum.set('x max bc', 'PML')
                self.lum.set('y min bc', 'Anti-Symmetric')
                self.lum.set('y max bc', 'PML')
                self.lum.set('z min bc', 'PML')
                self.lum.set('z max bc', 'PML')

                self.lum.set('mesh type', 'auto non-uniform')
                self.lum.set('min mesh step', x_res)
                self.lum.set('set simulation bandwidth', 0)
                self.lum.set('global source center wavelength', WaveLength)
                self.lum.set('global source wavelength span', 0)
                


                # Define Ports
                x = [0,EME_WGLength]
                x_Monitor = [0+0.1e-6,EME_WGLength-0.1e-6]
                yPos = [0, 0]
                yPos_span = [y_Port_Span, y_Port_Span]
                theta = [0, 0]
                direction = ['Forward', 'Backward']
                name = ['Input', 'Output']

                for i in range(2):
                    self.lum.addport()
                    self.lum.set('name', name[i])
                    self.lum.set("injection axis", "x-axis")
                    self.lum.set("x", x[i])
                    self.lum.set("y", yPos[i])
                    self.lum.set("y span", yPos_span[i])
                    self.lum.set("z", MonitorHeight)
                    self.lum.set("z span", z_Port_Span)
                    self.lum.set('direction', direction[i])
                    self.lum.set('mode selection', Mode)



                    # Power Monitor Port 1
                    self.lum.addpower()
                    self.lum.set('monitor type', "2D X-normal")
                    self.lum.set('name', name[i])
                    self.lum.set("x", x_Monitor[i] )
                    self.lum.set("y", yPos[i])
                    self.lum.set("y span", yPos_span[i])
                    self.lum.set("z", MonitorHeight)
                    self.lum.set("z span", z_Port_Span)
                    self.lum.set('output Px', 1)
                    self.lum.set('output Py', 1)
                    self.lum.set('output Pz', 1)
                    self.lum.set('output power', 1)

                # Add Movie monitor
                self.lum.addmovie()
                self.lum.set("y", 0)
                self.lum.set("y span", Device_Width)
                self.lum.set("z", MonitorHeight)
                self.lum.set("x", WG_Length/2)
                self.lum.set("x span", WG_Length)

                # Add Power and Freq Monitor
                self.lum.addpower()
                self.lum.set('monitor type', '2D Z-normal')
                self.lum.set("y", 0)
                self.lum.set("y span", Device_Width)
                self.lum.set("x", WG_Length/2)
                self.lum.set("x span", WG_Length)
                self.lum.set("z", MonitorHeight)
                self.lum.set('output Px', 1)
                self.lum.set('output Py', 1)
                self.lum.set('output Pz', 1)
                self.lum.set('output power', 1)


            else:

                # Calc Output Loc
                sideLength = np.cos(angle*np.pi/180)* WG_Length
                sideHight = np.sqrt(WG_Length**2 - sideLength**2)

                # Adds a FDTD Solver
                self.lum.addfdtd()
                self.lum.set("x", WG_Length/2)
                self.lum.set("x span", WG_Length)
                self.lum.set("y", WG_Length / 4)
                self.lum.set("y span", Device_Width)
                self.lum.set('simulation temperature', 273.15 + 20)
                self.lum.set("z", Substrate_Height)
                self.lum.set("z span", 4e-6)
                self.lum.set('z min bc', 'PML')
                self.lum.set('z max bc', 'PML')
                self.lum.set('mesh type', 'auto non-uniform')
                self.lum.set('min mesh step', x_res)
                self.lum.set('set simulation bandwidth', 0)
                self.lum.set('global source center wavelength', WaveLength)
                self.lum.set('global source wavelength span', 0)


                # Define Ports
                x = [0, EME_WGLength]
                x_Monitor = [0+0.1e-6,EME_WGLength-0.1e-6]
                yPos = [0, 0 + sideHight]
                yPos_span = [y_Port_Span, y_Port_Span]
                direction = ['Forward', 'Backward']
                name = ['Input', 'Output']
                theta = [angle, angle]

                for i in range(2):
                    self.lum.addport()
                    self.lum.set('name', name[i])
                    self.lum.set("injection axis", "x-axis")
                    self.lum.set("x", x[i])
                    self.lum.set("y", yPos[i])
                    self.lum.set("y span", yPos_span[i])
                    self.lum.set("z", MonitorHeight)
                    self.lum.set("z span", z_Port_Span)
                    self.lum.set('direction', direction[i])
                    self.lum.set('mode selection', Mode)
                    self.lum.set('theta',theta[i])

                    # Power Monitor Port 1
                    self.lum.addtime()
                    self.lum.set('name', name[i])
                    self.lum.set("x", x_Monitor[i] )
                    self.lum.set("y", yPos[i])
                    self.lum.set("z", MonitorHeight)
                    self.lum.set('output Px', 1)
                    self.lum.set('output Py', 1)
                    self.lum.set('output Pz', 1)

                # Add Movie monitor
                self.lum.addmovie()
                self.lum.set("y", WG_Length/4)
                self.lum.set("y span", Device_Width)
                self.lum.set("z", MonitorHeight)
                self.lum.set("x", WG_Length/2)
                self.lum.set("x span", WG_Length)

                # Add Power and Freq Monitor
                self.lum.addpower()
                self.lum.set('monitor type', '2D Z-normal')
                self.lum.set("y", WG_Length/4)
                self.lum.set("y span", Device_Width)
                self.lum.set("x", WG_Length/2)
                self.lum.set("x span", WG_Length)
                self.lum.set("z", MonitorHeight)
                self.lum.set('output Px', 1)
                self.lum.set('output Py', 1)
                self.lum.set('output Pz', 1)
                self.lum.set('output power', 1)



        elif Taper == True:

            # Device specifications
            Device_Width = 2*TaperWidth + WaveLength * 2  # MMI_Width

            max_slabH = Slab_Height
            MonitorHeight = Substrate_Height + max_slabH + WG_Height/2
            EME_WGLength = TaperLength * np.cos(angle * np.pi / 180)


            
            
            
            if TaperType == "Normal":
            
                # Adds a FDTD Solver
                self.lum.addfdtd()
                self.lum.set("x", 0)
                self.lum.set("x span", TaperLength)
                self.lum.set("y", 0)
                self.lum.set("y span", Device_Width)
                self.lum.set('simulation temperature', 273.15 + 20)
                self.lum.set("z", Substrate_Height)
                self.lum.set("z span", 4e-6)
                self.lum.set('x min bc', 'PML')
                self.lum.set('x max bc', 'PML')
                self.lum.set('y min bc', 'Anti-Symmetric')
                self.lum.set('y max bc', 'PML')
                self.lum.set('z min bc', 'PML')
                self.lum.set('z max bc', 'PML')

                self.lum.set('mesh type', 'auto non-uniform')
                self.lum.set('min mesh step', x_res)
                self.lum.set('set simulation bandwidth', 0)
                self.lum.set('global source center wavelength', WaveLength)
                self.lum.set('global source wavelength span', 0)
            
                # Define Ports
                Diff_Span = y_Port_Span - WG_Width
                x = [-TaperLength/2+ 0.1e-6 , TaperLength/2 - 0.1e-6]
                x_Monitor = [-TaperLength/2+ 0.1e-6+0.1e-6,TaperLength/2 - 0.1e-6 - 0.1e-6]
                yPos = [0, 0]
                yPos_span = [y_Port_Span, TaperWidth + Diff_Span]
                theta = [0, 0]
                direction = ['Forward', 'Backward']
                name = ['Input', 'Output']
                
                
                for i in range(2):
                    self.lum.addport()
                    self.lum.set('name', name[i])
                    self.lum.set("injection axis", "x-axis")
                    self.lum.set("x", x[i])
                    self.lum.set("y", yPos[i])
                    self.lum.set("y span", yPos_span[i])
                    self.lum.set("z", MonitorHeight)
                    self.lum.set("z span", z_Port_Span)
                    self.lum.set('direction', direction[i])
                    self.lum.set('mode selection', Mode)



                    # Power Monitor Port 1
                    self.lum.addtime()
                    self.lum.set('name', name[i])
                    self.lum.set("x", x_Monitor[i] )
                    self.lum.set("y", yPos[i])
                    self.lum.set("z", MonitorHeight)
                    self.lum.set('output Px', 1)
                    self.lum.set('output Py', 1)
                    self.lum.set('output Pz', 1)
                    
                    
                    
                    # Add Movie monitor
                    self.lum.addmovie()
                    self.lum.set("y", 0)
                    self.lum.set("y span", 2*TaperWidth + WaveLength * 2)
                    self.lum.set("z", MonitorHeight)
                    self.lum.set("x", 0)
                    self.lum.set("x span", TaperLength)

                    # Add Power and Freq Monitor
                    self.lum.addpower()
                    self.lum.set('monitor type', '2D Z-normal')
                    self.lum.set("y",0)
                    self.lum.set("y span", 2*TaperWidth + WaveLength * 2)
                    self.lum.set("x", 0)
                    self.lum.set("x span", TaperLength)
                    self.lum.set("z", MonitorHeight)
                    self.lum.set('output Px', 1)
                    self.lum.set('output Py', 1)
                    self.lum.set('output Pz', 1)
                    self.lum.set('output power', 1)

                    
            else:
                # Adds a FDTD Solver
                self.lum.addfdtd()
                self.lum.set("x", 0)
                self.lum.set("x span", TaperLength)
                self.lum.set("y", 0)
                self.lum.set("y span", 2*TaperWidthB + 2e-6)
                self.lum.set('simulation temperature', 273.15 + 20)
                self.lum.set("z", Substrate_Height + max_slabH + TaperHightB/2)
                self.lum.set("z span", 2e-6 +  TaperHightB  )
                self.lum.set('x min bc', 'PML')
                self.lum.set('x max bc', 'PML')
                self.lum.set('y min bc', 'Anti-Symmetric')
                self.lum.set('y max bc', 'PML')
                self.lum.set('z min bc', 'PML')
                self.lum.set('z max bc', 'PML')

                self.lum.set('mesh type', 'auto non-uniform')
                self.lum.set('min mesh step', x_res)
                self.lum.set('set simulation bandwidth', 0)
                self.lum.set('global source center wavelength', WaveLength)
                self.lum.set('global source wavelength span', 0)
                
                
                # Define Ports
                Diff_Span = y_Port_Span - WG_Width
                x = [-TaperLength/2+ 0.1e-6 , TaperLength/2 - 0.1e-6]
                x_Monitor = [-TaperLength/2+ 0.1e-6+0.1e-6,TaperLength/2 - 0.1e-6 - 0.1e-6]
                yPos = [0, 0]
                yPos_span = [y_Port_Span, TaperWidthB + Diff_Span ]
                theta = [0, 0]
                direction = ['Forward', 'Backward']
                name = ['Input', 'Output']
                
                yPos_span = [y_Port_Span, TaperWidthF + Diff_Span ]
                z_Pos = [Substrate_Height + max_slabH + TaperHightB/2, Substrate_Height + max_slabH + TaperHightF/2 ]
                z_Span = [ TaperHightB + z_Port_Span , z_Port_Span]# TaperHightF/2  + z_Port_Span]
                for i in range(2):
                    self.lum.addport()
                    self.lum.set('name', name[i])
                    self.lum.set("injection axis", "x-axis")
                    self.lum.set("x", x[i])
                    self.lum.set("y", yPos[i])
                    self.lum.set("y span", yPos_span[i])
                    self.lum.set("z", z_Pos[i])
                    self.lum.set("z span", z_Span[i])
                    self.lum.set('direction', direction[i])
                    self.lum.set('mode selection', Mode)



                    # Power Monitor Port 1
                    self.lum.addtime()
                    self.lum.set('name', name[i])
                    self.lum.set("x", x_Monitor[i] )
                    self.lum.set("y", yPos[i])
                    self.lum.set("z", z_Pos[i])
                    self.lum.set('output Px', 1)
                    self.lum.set('output Py', 1)
                    self.lum.set('output Pz', 1)
                


                # Add Movie monitor
                self.lum.addmovie()
                self.lum.set("y", 0)
                self.lum.set("y span", TaperWidthB + Diff_Span)
                self.lum.set("z", z_Pos[1])
                self.lum.set("x", 0)
                self.lum.set("x span", TaperLength)

                # Add Power and Freq Monitor
                self.lum.addpower()
                self.lum.set('monitor type', '2D Z-normal')
                self.lum.set("y",0)
                self.lum.set("y span", TaperWidthB + Diff_Span)
                self.lum.set("x", 0)
                self.lum.set("x span", TaperLength)
                self.lum.set("z", z_Pos[1])
                self.lum.set('output Px', 1)
                self.lum.set('output Py', 1)
                self.lum.set('output Pz', 1)
                self.lum.set('output power', 1)





    def setArcWaveguideFDTDSolver(self, Parameters):
        '''


        Parameters
        ----------
        Parameters : Dictionary
            Dictionary with all the data needet for the ArcWaveguideFDTDSolver.
            Parameters['Substrate Height'] : int/float
                Substrate height.
            Parameters['WG Height'] : int/float
                Waveguide hight. Also the height of the MMI section
            Parameters['WG Width'] : int/float
                Waveguide width.
            Parameters['x res'] : int/float
                EME Mesh resolutio,
            Parameters['Slab Height'] : int/float
                Height of the slab.
            Parameters['Wavelength'] : int/float
                Wavelength
            Parameters["S_Band Radius"] : int/float
                S-Bend Radius in um.
            Parameters['Arc deg'] : int/float
                Arc define the Arc of the curve. It can be 90 or 180 degrees only.
                This two will define an 1/4 of a circle or 1/2 of a circle.
            Parameters["Mode"] : str
                Mode to choose from ("fundamental TE mode", "fundamental TM mode", "fundamental mode")
            Parameters["Port Span"]: list of int/floats
                Parameters["Port Span"] = [Span of Port in x direction, Span of Port in y direction, Span of Port in z direction]

        Returns
        -------
        None.

        '''

        Substrate_Height = Parameters['Substrate Height']
        WG_Height = Parameters['WG Height']
        WG_Width = Parameters['WG Width']
        x_res = Parameters['x res']
        Slab_Height = Parameters['Slab Height']
        WaveLength = Parameters['Wavelength']
        radius = Parameters["S_Band Radius"]
        arc = Parameters['Arc deg']
        Mode = Parameters["Mode"]
        y_Port_Span = Parameters["Port Span"][1]
        z_Port_Span = Parameters["Port Span"][2]



        if arc == 90:

            # Device dimentions
            # magic number
            # The cubic Bezier curve using this magic number in the pole points approximates the semi-circile with least error
            m = 0.55191502449
            MonitorHeight = Substrate_Height + Slab_Height + WG_Height / 2

            # Adds a Eigenmode Expansion (EME) solver region to the MODE simulation environment.
            self.lum.addfdtd()
            self.lum.set("x", m * radius)
            self.lum.set("x span", (m * radius * 2 + WG_Width)+2e-6)
            self.lum.set("y", radius * m)
            self.lum.set("y span", m * radius * 2 + WG_Width)
            self.lum.set('simulation temperature', 273.15 + 20)
            self.lum.set("z", Substrate_Height)
            self.lum.set("z span", Substrate_Height + WG_Height*2)
            self.lum.set('z min bc', 'PML')
            self.lum.set('z max bc', 'PML')
            self.lum.set('mesh type', 'auto non-uniform')
            self.lum.set('min mesh step', x_res)
            self.lum.set('set simulation bandwidth', 0)
            self.lum.set('global source center wavelength', WaveLength)
            self.lum.set('global source wavelength span', 0)

            x = [0-0.2e-6, radius]
            direction = ['Forward', 'Forward']
            name = ['Input', 'Output']
            y = [radius, 0-0.2e-6]

            # Port 1
            self.lum.addport()
            self.lum.set('name', name[0])
            self.lum.set("injection axis", "x-axis")
            self.lum.set("x", x[0])
            self.lum.set("y", y[0])
            self.lum.set("y span", y_Port_Span)
            self.lum.set("z", MonitorHeight)
            self.lum.set("z span", z_Port_Span)
            self.lum.set('direction', direction[0])
            self.lum.set('mode selection', Mode)
            # self.lum.set("bent waveguide", 1)
            # self.lum.set("bend radius", radius)



            # Time Monitor Port 1
            self.lum.addtime()
            self.lum.set('name', name[0])
            self.lum.set("x", x[0] + 0.2e-6)
            self.lum.set("y", y[0])
            self.lum.set("z", MonitorHeight)
            self.lum.set('output Px', 1)
            self.lum.set('output Py', 1)
            self.lum.set('output Pz', 1)



            # Power Monitor Port 1
            self.lum.addpower()
            self.lum.set('name', "2D X-mormal Input Power Monitor")
            self.lum.set("monitor type", "2D X-normal")
            self.lum.set("x", x[0] + 0.1e-6)
            self.lum.set("y", y[0])
            self.lum.set("y span", y_Port_Span)
            self.lum.set("z", MonitorHeight)
            self.lum.set("z span", z_Port_Span)
            self.lum.set('output Px', 1)
            self.lum.set('output Py', 1)
            self.lum.set('output Pz', 1)



            #  Port 2
            self.lum.addport()
            self.lum.set('name', name[1])
            self.lum.set("injection axis", "y-axis")
            self.lum.set("x", x[1])
            self.lum.set("x span", y_Port_Span) # Only becouse we just rotate the previus WG on 90 degrees!!!
            self.lum.set("y", y[1])
            self.lum.set("z", MonitorHeight)
            self.lum.set("z span", z_Port_Span)
            self.lum.set('direction', direction[1])
            self.lum.set('mode selection', Mode)
            # self.lum.set("bent waveguide", 1)
            # self.lum.set("bend radius", radius)
            # self.lum.set("bend orientation", 90)


            # Time Monitor Port 2
            self.lum.addtime()
            self.lum.set('name', name[1])
            self.lum.set("x", x[1])
            self.lum.set("y", y[1] + 0.2e-6)
            self.lum.set("z", MonitorHeight)
            self.lum.set('output Px', 1)
            self.lum.set('output Py', 1)
            self.lum.set('output Pz', 1)

            # Power Monitor Port 2
            self.lum.addpower()
            self.lum.set('name', "2D Y-mormal Input Power Monitor")
            self.lum.set("monitor type","2D Y-normal")
            self.lum.set("x", x[1])
            self.lum.set("x span", y_Port_Span)  # Only becouse we just rotate the previus WG on 90 degrees!!!
            self.lum.set("y", y[1]+0.1e-6)
            self.lum.set("z", MonitorHeight)
            self.lum.set("z span", z_Port_Span)
            self.lum.set('output Px', 1)
            self.lum.set('output Py', 1)
            self.lum.set('output Pz', 1)


            # Add Movie Monitor
            self.lum.addmovie()
            self.lum.set("y", (m * radius))
            self.lum.set("y span", m * radius * 2 + WG_Width)
            self.lum.set("z", MonitorHeight)
            self.lum.set("x", m * radius)
            self.lum.set("x span", (m * radius * 2 + WG_Width))

            # Add Power and Freq Monitor
            self.lum.addpower()
            self.lum.set('monitor type', '2D Z-normal')
            self.lum.set("y", (m * radius))
            self.lum.set("y span", m * radius * 2 + WG_Width)
            self.lum.set("x", m * radius)
            self.lum.set("x span", (m * radius * 2 + WG_Width))
            self.lum.set("z", MonitorHeight)
            self.lum.set('output Px', 1)
            self.lum.set('output Py', 1)
            self.lum.set('output Pz', 1)
            self.lum.set('output power', 1)

        elif arc == 180:

            # magic number
            # The cubic Bezier curve using this magic number in the pole points approximates the semi-circile with least error
            m = 0.55191502449
            MonitorHeight = Substrate_Height + Slab_Height + WG_Height / 2

            # Adds a Eigenmode Expansion (EME) solver region to the MODE simulation environment.
            self.lum.addfdtd()
            self.lum.set("x", 0)
            self.lum.set("x span", (m * radius * 2 + WG_Width) * 2)
            self.lum.set("y", radius * m)
            self.lum.set("y span", m * radius * 2 + WG_Width)
            self.lum.set('simulation temperature', 273.15 + 20)
            self.lum.set("z", Substrate_Height)
            self.lum.set("z span", Substrate_Height + WG_Height*2)
            self.lum.set('z min bc', 'PML')
            self.lum.set('z max bc', 'PML')
            self.lum.set('mesh type', 'auto non-uniform')
            self.lum.set('min mesh step', x_res)
            self.lum.set('set simulation bandwidth', 0)
            self.lum.set('global source center wavelength', WaveLength)
            self.lum.set('global source wavelength span', 0)

            x = [-radius, radius]
            direction = ['Forward', 'Forward']
            name = ['Input', 'Output']
            y = [0, 0]

            self.lum.addport()
            self.lum.set('name', name[0])
            self.lum.set("injection axis", "y-axis")
            self.lum.set("x", x[0])
            self.lum.set("x span", y_Port_Span) # Only becouse we just rotate the previus WG on 90 degrees!!!
            self.lum.set("y", y[0])
            self.lum.set("z", MonitorHeight)
            self.lum.set("z span", z_Port_Span)
            self.lum.set('direction', direction[0])
            self.lum.set('mode selection', Mode)
            self.lum.set("bent waveguide", 1)
            self.lum.set("bend radius", radius)
            self.lum.set("bend orientation", -90)

            # Power Monitor Port 1
            self.lum.addtime()
            self.lum.set('name', name[0])
            self.lum.set("x", x[0])
            self.lum.set("y", y[0] + 0.2e-6)
            self.lum.set("z", MonitorHeight)
            self.lum.set('output Px', 1)
            self.lum.set('output Py', 1)
            self.lum.set('output Pz', 1)

            self.lum.addport()
            self.lum.set('name', name[1])
            self.lum.set("injection axis", "y-axis")
            self.lum.set("x", x[1])
            self.lum.set("x span", y_Port_Span) # Only becouse we just rotate the previus WG on 90 degrees!!!
            self.lum.set("y", y[1])
            self.lum.set("z", MonitorHeight)
            self.lum.set("z span", z_Port_Span)
            self.lum.set('direction', direction[1])
            self.lum.set('mode selection', Mode)
            self.lum.set("bent waveguide", 1)
            self.lum.set("bend radius", radius)
            self.lum.set("bend orientation", 90)

            # Power Monitor Port 1
            self.lum.addtime()
            self.lum.set('name', name[1])
            self.lum.set("x", x[1])
            self.lum.set("y", y[1] + 0.2e-6)
            self.lum.set("z", MonitorHeight)
            self.lum.set('output Px', 1)
            self.lum.set('output Py', 1)
            self.lum.set('output Pz', 1)

            self.lum.addmovie()
            self.lum.set("x", 0)
            self.lum.set("x span", (m * radius * 2 + WG_Width) * 2)
            self.lum.set("y", radius * m)
            self.lum.set("y span", m * radius * 2 + WG_Width)
            self.lum.set("z", MonitorHeight)

            # Add Power and Freq Monitor
            self.lum.addpower()
            self.lum.set('monitor type', '2D Z-normal')
            self.lum.set("x", 0)
            self.lum.set("x span", (m * radius * 2 + WG_Width) * 2)
            self.lum.set("y", radius * m)
            self.lum.set("y span", m * radius * 2 + WG_Width)
            self.lum.set("z", MonitorHeight)
            self.lum.set('output Px', 1)
            self.lum.set('output Py', 1)
            self.lum.set('output Pz', 1)
            self.lum.set('output power', 1)





    def setBendWaveguideFDTDSolver(self, Parameters):

        '''


        Parameters
        ----------
        Parameters : Dictionary
            Dictionary with all the data needet for the BendWaveguideFDTDSolver.
            Parameters['Substrate Height'] : int/float
               Substrate height.
            Parameters['WG Height'] : int/float
               Waveguide hight. Also the height of the MMI section
            Parameters['x res'] : int/float
                 EME Mesh resolutio,
            Parameters['Slab Height'] : int/float
               Height of the slab.
            Parameters['Wavelength'] : int/float
               Wavelength
            Parameters["x span"] : int/float
               Length of the S-Bend Waveguide
            Parameters["y span"] : int/float
               Width of the S-Bend Waveguide
            Parameters["Mode"] : str
                Mode to choose from ("fundamental TE mode", "fundamental TM mode", "fundamental mode")
            Parameters["Port Span"]: list of int/floats
                Parameters["Port Span"] = [Span of Port in x direction, Span of Port in y direction, Span of Port in z direction]

        Returns
        -------
        None.

        '''

        Substrate_Height = Parameters['Substrate Height']
        WG_Height = Parameters['WG Height']
        x_res = Parameters['x res']
        Slab_Height = Parameters['Slab Height']
        WaveLength = Parameters['Wavelength']
        x_span = Parameters["x span"]
        y_span = Parameters["y span"]
        Mode = Parameters["Mode"]
        y_Port_Span = Parameters["Port Span"][1]
        z_Port_Span = Parameters["Port Span"][2]

        # Device specifications
        Device_Length = x_span
        Device_Width = y_span + WaveLength * 2

        max_slabH = Slab_Height
        MonitorHeight = Substrate_Height + max_slabH + WG_Height / 2
        
        # Check Object Geometry
        self.lum.select("S Bend::S-Bend")
        z_Structure = self.lum.get("z")
        z_Object = self.lum.get("Base Height")
        Z = z_Structure + z_Object

        # Adds a Eigenmode Expansion (EME) solver region to the MODE simulation environment.
        self.lum.addfdtd()
        self.lum.set("x", Device_Length / 2)
        self.lum.set("x span", Device_Length)
        self.lum.set("y", y_span / 2)
        self.lum.set("y span", Device_Width)
        self.lum.set('simulation temperature', 273.15 + 20)
        self.lum.set("y span", Device_Width)
        self.lum.set("z", Substrate_Height)
        self.lum.set("z span", Z + 4e-6)
        self.lum.set('z min bc', 'PML')
        self.lum.set('z max bc', 'PML')
        self.lum.set('mesh type', 'auto non-uniform')
        self.lum.set('min mesh step', x_res)
        self.lum.set('set simulation bandwidth', 0)
        self.lum.set('global source center wavelength', WaveLength)
        self.lum.set('global source wavelength span', 0)

        x = [0, x_span]
        direction = ['Forward', 'Backward']
        name = ['Input', 'Output']
        y = [0, y_span]

        self.lum.addport()
        self.lum.set('name', name[0])
        self.lum.set("x", x[0])
        self.lum.set("y", y[0])
        self.lum.set("y span", y_Port_Span)
        self.lum.set("z", MonitorHeight)
        self.lum.set("z span", z_Port_Span)
        self.lum.set('direction', direction[0])
        self.lum.set('mode selection', Mode)


        # Power Monitor Port 1
        self.lum.addpower()
        self.lum.set('name', name[0])
        self.lum.set('monitor type', '2D X-normal')
        self.lum.set("x", x[0]+0.1e-6)
        self.lum.set("y", y[0])
        self.lum.set("y span", y_Port_Span)
        self.lum.set("z", MonitorHeight)
        self.lum.set("z span", z_Port_Span)
        self.lum.set('output Px', 1)
        self.lum.set('output Py', 1)
        self.lum.set('output Pz', 1)
        self.lum.set('output power', 1)



        self.lum.addport()
        self.lum.set('name', name[1])
        self.lum.set("x", x[1])
        self.lum.set("y", y[1])
        self.lum.set("y span", y_Port_Span)
        self.lum.set("z", MonitorHeight)
        self.lum.set("z span", z_Port_Span)
        self.lum.set('direction', direction[1])
        self.lum.set('mode selection', Mode)

        # Power Monitor Port 1
        self.lum.addpower()
        self.lum.set('monitor type', '2D X-normal')
        self.lum.set('name', name[1])
        self.lum.set("x", x[1]-0.1e-6)
        self.lum.set("y", y[1])
        self.lum.set("y span", y_Port_Span)
        self.lum.set("z", MonitorHeight)
        self.lum.set("z span", z_Port_Span)
        self.lum.set('output Px', 1)
        self.lum.set('output Py', 1)
        self.lum.set('output Pz', 1)
        self.lum.set('output power', 1)



        self.lum.addmovie()
        self.lum.set('x', x_span / 2)
        self.lum.set("x span", x_span)
        self.lum.set("y", y_span / 2)
        self.lum.set("y span", Device_Width)
        self.lum.set("z", MonitorHeight)



        # Add Power and Freq Monitor
        self.lum.addpower()
        self.lum.set('name', 'FieldMonitor')
        self.lum.set('monitor type', '2D Z-normal')
        self.lum.set("x", x_span / 2)
        self.lum.set("x span", x_span)
        self.lum.set("y", y_span / 2)
        self.lum.set("y span", Device_Width)
        self.lum.set('z', MonitorHeight)
        self.lum.set('output Px', 1)
        self.lum.set('output Py', 1)
        self.lum.set('output Pz', 1)
        self.lum.set('output power', 1)





    def setMMI2x2FDTDSolver(self, Parameters):
        '''


        Parameters
        ----------
        Parameters : Dictionary
            Dictionary with all the data needet for the MMI2x2FDTDSolver.
            Parameters['Substrate Height'] : int/float
                Height of the slab.
            Parameters['MMI Width'] : int/float
                Width of MMI
            Parameters['MMI Length'] : int/float
                Length of MMI
            Parameters['WG Height'] : int/float
                Heigth of waveguide
            Parameters['WG Width'] : int/float
                Width og waveguide
            Parameters['WG Length'] : int/float
                Length of waveguide
            Parameters['Position Offset'] : int/float
                Offset between the waveguides. If Taper == True then this become the offset
                betweent he tapers wider sides. Waveguide and Tapers cannot be placed ourside
                the MMI structure. The minimum distance between Taper and Waveguide is 1 um
                becouse of manufactering restrictions in the University.
            Parameters["Taper"]: boolen 
                Add Taper to the structure on the input and output waveguids
            Parameters['Taper Length']: int/float
                Lenght of the Taper in Parameters["Taper"] = True
            Parameters['x res'] : int/float
                Mesh cell sizes.
            Parameters['Slab Height'] : int/float
                Height of the slab.
            Parameters['Wavelength'] : int/float
                Wavelength.
            Parameters["Mode"] : str
                Mode to choose from ("fundamental TE mode", "fundamental TM mode", "fundamental mode")
            Parameters["Port Span"]: list of int/floats
                Parameters["Port Span"] = [Span of Port in x direction, Span of Port in y direction, Span of Port in z direction]

        Returns
        -------
        None.

        '''

        Substrate_Height = Parameters['Substrate Height']
        MMI_Width = Parameters['MMI Width']
        MMI_Length = Parameters['MMI Length']
        WG_Height = Parameters['WG Height']
        WG_Width = Parameters['WG Width']
        WG_Length = Parameters['WG Length']
        posOffset = Parameters['Position Offset']
        x_res = Parameters['x res']
        Slab_Height = Parameters['Slab Height']
        WaveLength = Parameters['Wavelength']
        TaperLength = Parameters['Taper Length']
        Taper = Parameters['Taper']
        Mode = Parameters["Mode"]
        y_Port_Span = Parameters["Port Span"][1]
        z_Port_Span = Parameters["Port Span"][2]


        if Taper == False:

            # Device specifications
            Device_Length = MMI_Length + 2 * WG_Length
            Device_Width = MMI_Width + WaveLength * 2  # MMI_Width
            max_slabH = Slab_Height
            # Ports_mid = Substrate_Height + (max_slabH + WG_Height) / 2
            Ports_mid = Substrate_Height + max_slabH + WG_Height/2

            # Adds a Eigenmode Expansion (EME) solver region to the MODE simulation environment.
            self.lum.addfdtd()
            self.lum.set("x min", -(Device_Length / 2))
            self.lum.set("x max", (Device_Length / 2))
            self.lum.set("y", 0)
            self.lum.set('simulation temperature', 273.15 + 20)
            self.lum.set("y span", Device_Width)
            self.lum.set("z", Substrate_Height)
            self.lum.set("z span", 4e-6)
            self.lum.set('z min bc', 'PML')
            self.lum.set('z max bc', 'PML')
            self.lum.set('mesh type', 'auto non-uniform')
            self.lum.set('min mesh step', x_res)
            self.lum.set('set simulation bandwidth', 0)
            self.lum.set('global source center wavelength', WaveLength)
            self.lum.set('global source wavelength span', 0)

            x = [(Device_Length / 2), (Device_Length / 2), -(Device_Length / 2), -(Device_Length / 2)]
            direction = ['Backward', 'Backward', 'Forward', 'Forward']


            name = ['Input_L', 'Input_R', 'Output_L', 'Output_R']
            yPort_vec = [posOffset / 2 + WG_Width / 2, -(posOffset / 2 + WG_Width / 2), posOffset / 2 + WG_Width / 2, -(posOffset / 2 + WG_Width / 2)]

            overLapp = yPort_vec[0] - y_Port_Span/2
            if overLapp < 0:
                raise ValueError("!!! CAUTION !!! - The Ports are overlapping at the middle! Please change the Y Port Span or move the Waveguides away from each other!")
            else:
                pass

            self.lum.addmovie()
            self.lum.set("x min", -(0.5e-6 + Device_Length / 2))
            self.lum.set("x max", (0.5e-6 + Device_Length / 2))
            self.lum.set("y min", -(Device_Width + 1e-6))
            self.lum.set("y max", (Device_Width + 1e-6))
            self.lum.set("z", Ports_mid)


            PortCorrection = [-0.1e-6, -0.1e-6, 0.1e-6, 0.1e-6]
            for i in range(4):

                # Power Monitor Port 1
                self.lum.addpower()
                self.lum.set('name', "Power_"+ name[i])
                self.lum.set('monitor type', '2D X-normal')
                self.lum.set("x", x[i]+ PortCorrection[i])
                self.lum.set("y", yPort_vec[i])
                self.lum.set("y span", y_Port_Span)
                self.lum.set("z", Ports_mid)
                self.lum.set("z span", z_Port_Span)
                self.lum.set('output Px', 1)
                self.lum.set('output Py', 1)
                self.lum.set('output Pz', 1)
                self.lum.set('output power', 1)

                # Ports
                self.lum.addport()
                self.lum.set('name', name[i])
                self.lum.set("x", x[i])
                self.lum.set('y', yPort_vec[i])
                self.lum.set('y span', y_Port_Span)
                self.lum.set("z", Ports_mid)
                self.lum.set("z span", z_Port_Span)
                self.lum.set('direction', direction[i])
                self.lum.set('mode selection', Mode)

            # Add Power and Freq Monitor
            self.lum.addpower()
            self.lum.set('monitor type', '2D Z-normal')
            self.lum.set("x min", -(0.5e-6 + Device_Length / 2))
            self.lum.set("x max", (0.5e-6 + Device_Length / 2))
            self.lum.set("y min", -(Device_Width + 1e-6))
            self.lum.set("y max", (Device_Width + 1e-6))
            self.lum.set('z', Substrate_Height + Slab_Height + WG_Height / 2)
            self.lum.set('output Px', 1)
            self.lum.set('output Py', 1)
            self.lum.set('output Pz', 1)
            self.lum.set('output power', 1)


        elif Taper == True:
            # Device specifications
            Device_Length = MMI_Length + 2 * WG_Length + 2*TaperLength
            Device_Width = MMI_Width + WaveLength * 2  # MMI_Width
            max_slabH = Slab_Height
            # Ports_mid = Substrate_Height + (max_slabH + WG_Height) / 2
            Ports_mid = Substrate_Height + max_slabH + WG_Height/2

            # Adds a Eigenmode Expansion (EME) solver region to the MODE simulation environment.
            self.lum.addfdtd()
            self.lum.set("x min", -(Device_Length / 2))
            self.lum.set("x max", (Device_Length / 2))
            self.lum.set("y", 0)
            self.lum.set('simulation temperature', 273.15 + 20)
            self.lum.set("y span", Device_Width)
            self.lum.set("z", Substrate_Height)
            self.lum.set("z span", 4e-6)
            self.lum.set('z min bc', 'PML')
            self.lum.set('z max bc', 'PML')
            self.lum.set('mesh type', 'auto non-uniform')
            self.lum.set('min mesh step', x_res)
            self.lum.set('set simulation bandwidth', 0)
            self.lum.set('global source center wavelength', WaveLength)
            self.lum.set('global source wavelength span', 0)

            x = [(Device_Length / 2), (Device_Length / 2), -(Device_Length / 2), -(Device_Length / 2)]
            direction = ['Backward', 'Backward', 'Forward', 'Forward']


            name = ['Input_L', 'Input_R', 'Output_L', 'Output_R']
            yPort_vec = [posOffset / 2 + WG_Width / 2, -(posOffset / 2 + WG_Width / 2), posOffset / 2 + WG_Width / 2, -(posOffset / 2 + WG_Width / 2)]

            overLapp = yPort_vec[0] - y_Port_Span/2
            if overLapp < 0:
                raise ValueError("!!! CAUTION !!! - The Ports are overlapping at the middle! Please change the Y Port Span or move the Waveguides away from each other!")
            else:
                pass

            self.lum.addmovie()
            self.lum.set("x min", -(0.5e-6 + Device_Length / 2))
            self.lum.set("x max", (0.5e-6 + Device_Length / 2))
            self.lum.set("y min", -(Device_Width + 1e-6))
            self.lum.set("y max", (Device_Width + 1e-6))
            self.lum.set("z", Ports_mid)


            PortCorrection = [-0.1e-6, -0.1e-6, 0.1e-6, 0.1e-6]
            for i in range(4):

                # Power Monitor Port 1
                self.lum.addpower()
                self.lum.set('name', "Power_"+ name[i])
                self.lum.set('monitor type', '2D X-normal')
                self.lum.set("x", x[i]+ PortCorrection[i])
                self.lum.set("y", yPort_vec[i])
                self.lum.set("y span", y_Port_Span)
                self.lum.set("z", Ports_mid)
                self.lum.set("z span", z_Port_Span)
                self.lum.set('output Px', 1)
                self.lum.set('output Py', 1)
                self.lum.set('output Pz', 1)
                self.lum.set('output power', 1)

                # Ports
                self.lum.addport()
                self.lum.set('name', name[i])
                self.lum.set("x", x[i])
                self.lum.set('y', yPort_vec[i])
                self.lum.set('y span', y_Port_Span)
                self.lum.set("z", Ports_mid)
                self.lum.set("z span", z_Port_Span)
                self.lum.set('direction', direction[i])
                self.lum.set('mode selection', Mode)

            # Add Power and Freq Monitor
            self.lum.addpower()
            self.lum.set('monitor type', '2D Z-normal')
            self.lum.set("x min", -(0.5e-6 + Device_Length / 2))
            self.lum.set("x max", (0.5e-6 + Device_Length / 2))
            self.lum.set("y min", -(Device_Width + 1e-6))
            self.lum.set("y max", (Device_Width + 1e-6))
            self.lum.set('z', Substrate_Height + Slab_Height + WG_Height / 2)
            self.lum.set('output Px', 1)
            self.lum.set('output Py', 1)
            self.lum.set('output Pz', 1)
            self.lum.set('output power', 1)

        else:
            raise ValueError("Incorect Taper variable!. Possible Taper values are Taper = False or Taper = True.")





    def setMMI2x1FDTDSolver(self, Parameters):
        '''


        Parameters
        ----------
        Parameters : Dictionary
            Dictionary with all the data needet for the MMI2x1FDTDSolver.
            Parameters['Substrate Height'] : int/float
                Height of the slab.
            Parameters['MMI Width'] : int/float
                Width of MMI
            Parameters['MMI Length'] : int/float
                Length of MMI
            Parameters['WG Height'] : int/float
                Heigth of waveguide
            Parameters['WG Width'] : int/float
                Width og waveguide
            Parameters['WG Length'] : int/float
                Length of waveguide
            Parameters['Position Offset'] : int/float
                Offset between the waveguides. If Taper == True then this become the offset
                betweent he tapers wider sides. Waveguide and Tapers cannot be placed ourside
                the MMI structure. The minimum distance between Taper and Waveguide is 1 um
                becouse of manufactering restrictions in the University.
            Parameters['Offset Input'] : int/float
                Input waveguide/taper offset.
            Parameters["Taper"]: boolen 
                Add Taper to the structure on the input and output waveguids
            Parameters['Taper Length']: int/float
                Lenght of the Taper in Parameters["Taper"] = True
            Parameters['x res'] : int/float
                Mesh cell sizes.
            Parameters['Slab Height'] : int/float
                Height of the slab.
            Parameters['Wavelength'] : int/float
                Wavelength.
            Parameters["Mode"] : str
                Mode to choose from ("fundamental TE mode", "fundamental TM mode", "fundamental mode")
            Parameters["Port Span"]: list of int/floats
                Parameters["Port Span"] = [Span of Port in x direction, Span of Port in y direction, Span of Port in z direction]
            Parameters["Offset Output"] : anything, optional
                This function will allow the user to move the outputs in oposite direction. Please dont use it since is there only 
                becouse the maschine of our physic departmant had some proiblems with the LNOI objects design. 

        Returns
        -------
        None.

        '''

        Substrate_Height = Parameters['Substrate Height']
        MMI_Width = Parameters['MMI Width']
        MMI_Length = Parameters['MMI Length']

        WG_Height = Parameters['WG Height']
        WG_Width = Parameters['WG Width']
        WG_Length = Parameters['WG Length']
        OffsetInput = Parameters['Offset Input']
        posOffset = Parameters['Position Offset']
        x_res = Parameters['x res']
        Slab_Height = Parameters['Slab Height']
        WaveLength = Parameters['Wavelength']
        TaperLength = Parameters['Taper Length']
        Taper = Parameters['Taper']
        Mode = Parameters["Mode"]
        y_Port_Span = Parameters["Port Span"][1]
        z_Port_Span = Parameters["Port Span"][2]
        if 'Offset Output' not in list(Parameters.keys()):
            OffsetOutput = None
        else:
            OffsetOutput = Parameters['Offset Output']
        


        if Taper == False:



            # Device specifications
            Device_Length = MMI_Length + 2 * WG_Length
            Device_Width = MMI_Width + WaveLength * 2  # MMI_Width
            max_slabH = Slab_Height
            # Ports_mid = Substrate_Height + (max_slabH + WG_Height) / 2
            # Ports_mid = Substrate_Height + Slab_Height + WG_Height / 2
            # MonitorHeight = Substrate_Height + (max_slabH + WG_Height) / 2
            Ports_mid = Substrate_Height + max_slabH + WG_Height/2
            # WG_H = WG_Height

            # Adds a Eigenmode Expansion (EME) solver region to the MODE simulation environment.
            self.lum.addfdtd()
            self.lum.set("x min", -(Device_Length / 2))
            self.lum.set("x max", (Device_Length / 2))
            self.lum.set("y", 0)
            self.lum.set('simulation temperature', 273.15 + 20)
            self.lum.set("y span", Device_Width)
            self.lum.set("z", Substrate_Height)
            self.lum.set("z span", 4e-6)
            self.lum.set('z min bc', 'PML')
            self.lum.set('z max bc', 'PML')
            self.lum.set('mesh type', 'auto non-uniform')
            self.lum.set('min mesh step', x_res)
            self.lum.set('set simulation bandwidth', 0)
            self.lum.set('global source center wavelength', WaveLength)
            self.lum.set('global source wavelength span', 0)

            # Positions of the Input and Output WGs
            # Triangle EQ for MMI Width
            x = [(Device_Length / 2), -(Device_Length / 2), -(Device_Length / 2)]
            direction = ['Backward', 'Forward', 'Forward']

            name = ['Input', 'Output_L', 'Output_R']
            
            if OffsetOutput == None:
                yPort_vec = [OffsetInput, -(posOffset / 2 + WG_Width / 2) , (posOffset / 2 + WG_Width / 2)  ]
            else:
                yPort_vec = [OffsetInput, -(posOffset / 2 + WG_Width / 2) + OffsetOutput, (posOffset / 2 + WG_Width / 2) + OffsetOutput]

            overLapp = yPort_vec[2] - y_Port_Span/2
            if overLapp < 0:
                y_Port_Span_old = y_Port_Span
                y_Port_Span = yPort_vec[2]*2
                print((f"!!! CAUTION !!! - The Ports are overlapping at the middle! So the y Span of the ports will be reduze from y_span = {y_Port_Span_old} to y_span_new = {y_Port_Span}"))

                self.lum.addmovie()
                self.lum.set("x min", -(0.5e-6 + Device_Length / 2))
                self.lum.set("x max", (0.5e-6 + Device_Length / 2))
                self.lum.set("y min", -(Device_Width + 1e-6))
                self.lum.set("y max", (Device_Width + 1e-6))
                self.lum.set("z", Ports_mid)

                PortCorrection = [-0.1e-6, 0.1e-6, 0.1e-6]

                for i in range(3):
                    # Power Monitor Port 1
                    self.lum.addpower()
                    self.lum.set('name', "Power_" + name[i])
                    self.lum.set('monitor type', '2D X-normal')
                    self.lum.set("x", x[i] + PortCorrection[i])
                    self.lum.set("y", yPort_vec[i])
                    self.lum.set("y span", y_Port_Span_old)
                    self.lum.set("z", Ports_mid)
                    self.lum.set("z span", z_Port_Span)
                    self.lum.set('output Px', 1)
                    self.lum.set('output Py', 1)
                    self.lum.set('output Pz', 1)
                    self.lum.set('output power', 1)

                    self.lum.addport()
                    self.lum.set('name', name[i])
                    self.lum.set("x", x[i])
                    self.lum.set("y", yPort_vec[i])
                    self.lum.set("y span", y_Port_Span_old)
                    self.lum.set("z", Ports_mid)
                    self.lum.set("z span", z_Port_Span)
                    self.lum.set('direction', direction[i])
                    self.lum.set('mode selection', Mode)

                # Add Power and Freq Monitor
                self.lum.addpower()
                self.lum.set('monitor type', '2D Z-normal')
                self.lum.set("x min", -(0.5e-6 + Device_Length / 2))
                self.lum.set("x max", (0.5e-6 + Device_Length / 2))
                self.lum.set("y min", -(Device_Width + 1e-6))
                self.lum.set("y max", (Device_Width + 1e-6))
                self.lum.set('z', Substrate_Height + Slab_Height + WG_Height / 2)
                self.lum.set('output Px', 1)
                self.lum.set('output Py', 1)
                self.lum.set('output Pz', 1)
                self.lum.set('output power', 1)

                # raise ValueError("!!! CAUTION !!! - The Ports are overlapping at the middle! Please change the Y Port Span or move the Waveguides away from each other!")

            else:
                print("YOu are in the else part Now")
                self.lum.addmovie()
                self.lum.set("x min", -(0.5e-6 + Device_Length / 2))
                self.lum.set("x max", (0.5e-6 + Device_Length / 2))
                self.lum.set("y min", -(Device_Width + 1e-6))
                self.lum.set("y max", (Device_Width + 1e-6))
                self.lum.set("z", Ports_mid)

                PortCorrection = [-0.1e-6, 0.1e-6, 0.1e-6]

                for i in range(3):

                    # Power Monitor Port 1
                    self.lum.addpower()
                    self.lum.set('name', "Power_"+ name[i])
                    self.lum.set('monitor type', '2D X-normal')
                    self.lum.set("x", x[i] + PortCorrection[i])
                    self.lum.set("y", yPort_vec[i])
                    self.lum.set("y span", y_Port_Span)
                    self.lum.set("z", Ports_mid)
                    self.lum.set("z span", z_Port_Span)
                    self.lum.set('output Px', 1)
                    self.lum.set('output Py', 1)
                    self.lum.set('output Pz', 1)
                    self.lum.set('output power', 1)


                    self.lum.addport()
                    self.lum.set('name', name[i])
                    self.lum.set("x", x[i])
                    self.lum.set("y", yPort_vec[i])
                    self.lum.set("y span", y_Port_Span)
                    self.lum.set("z", Ports_mid)
                    self.lum.set("z span", z_Port_Span)
                    self.lum.set('direction', direction[i])
                    self.lum.set('mode selection', Mode)



                # Add Power and Freq Monitor
                self.lum.addpower()
                self.lum.set('monitor type', '2D Z-normal')
                self.lum.set("x min", -(0.5e-6 + Device_Length / 2))
                self.lum.set("x max", (0.5e-6 + Device_Length / 2))
                self.lum.set("y min", -(Device_Width + 1e-6))
                self.lum.set("y max", (Device_Width + 1e-6))
                self.lum.set('z', Substrate_Height + Slab_Height + WG_Height / 2)
                self.lum.set('output Px', 1)
                self.lum.set('output Py', 1)
                self.lum.set('output Pz', 1)
                self.lum.set('output power', 1)




        elif Taper == True:

            # Device specifications
            Device_Length = MMI_Length + 2 * WG_Length +2*TaperLength
            Device_Width = MMI_Width + WaveLength * 2  # MMI_Width
            max_slabH = Slab_Height
            # Ports_mid = Substrate_Height + (max_slabH + WG_Height) / 2
            Ports_mid = Substrate_Height + max_slabH + WG_Height/2
            # MonitorHeight = Substrate_Height + (max_slabH + WG_Height) / 2
            # WG_H = WG_Height

            # Adds a Eigenmode Expansion (EME) solver region to the MODE simulation environment.
            self.lum.addfdtd()
            self.lum.set("x min", -(Device_Length / 2))
            self.lum.set("x max", (Device_Length / 2))
            self.lum.set("y", 0)
            self.lum.set('simulation temperature', 273.15 + 20)
            self.lum.set("y span", Device_Width)
            self.lum.set("z", Substrate_Height)
            self.lum.set("z span", 4e-6)
            self.lum.set('z min bc', 'PML')
            self.lum.set('z max bc', 'PML')
            self.lum.set('mesh type', 'auto non-uniform')
            self.lum.set('min mesh step', x_res)
            self.lum.set('set simulation bandwidth', 0);
            self.lum.set('global source center wavelength', WaveLength)
            self.lum.set('global source wavelength span', 0)

            # Positions of the Input and Output WGs
            # Triangle EQ for MMI Width
            x = [(Device_Length / 2), -(Device_Length / 2), -(Device_Length / 2)]
            direction = ['Backward', 'Forward', 'Forward']

            name = ['Input', 'Output_L', 'Output_R']

            yPort_vec = [OffsetInput, -(posOffset / 2 + WG_Width / 2), posOffset / 2 + WG_Width / 2]

            overLapp = yPort_vec[2] - y_Port_Span/2
            if overLapp < 0:
                y_Port_Span_old = y_Port_Span
                y_Port_Span = yPort_vec[2] * 2
                print(f"!!! CAUTION !!! - The Ports are overlapping at the middle! So the y Span of the ports will be reduze from y_span = {y_Port_Span_old} to y_span_new = {y_Port_Span}")
                self.lum.addmovie()
                self.lum.set("x min", -(0.5e-6 + Device_Length / 2))
                self.lum.set("x max", (0.5e-6 + Device_Length / 2))
                self.lum.set("y min", -(Device_Width + 1e-6))
                self.lum.set("y max", (Device_Width + 1e-6))
                self.lum.set("z", Ports_mid)

                PortCorrection = [-0.1e-6, 0.1e-6, 0.1e-6]

                for i in range(3):
                    # Power Monitor Port 1
                    self.lum.addpower()
                    self.lum.set('name', "Power_" + name[i])
                    self.lum.set('monitor type', '2D X-normal')
                    self.lum.set("x", x[i] + PortCorrection[i])
                    self.lum.set("y", yPort_vec[i])
                    self.lum.set("y span", y_Port_Span)
                    self.lum.set("z", Ports_mid)
                    self.lum.set("z span", z_Port_Span)
                    self.lum.set('output Px', 1)
                    self.lum.set('output Py', 1)
                    self.lum.set('output Pz', 1)
                    self.lum.set('output power', 1)

                    self.lum.addport()
                    self.lum.set('name', name[i])
                    self.lum.set("x", x[i])
                    self.lum.set("y", yPort_vec[i])
                    self.lum.set("y span", y_Port_Span)
                    self.lum.set("z", Ports_mid)
                    self.lum.set("z span", z_Port_Span)
                    self.lum.set('direction', direction[i])
                    self.lum.set('mode selection', Mode)

                # Add Power and Freq Monitor
                self.lum.addpower()
                self.lum.set('monitor type', '2D Z-normal')
                self.lum.set("x min", -(0.5e-6 + Device_Length / 2))
                self.lum.set("x max", (0.5e-6 + Device_Length / 2))
                self.lum.set("y min", -(Device_Width + 1e-6))
                self.lum.set("y max", (Device_Width + 1e-6))
                self.lum.set('z', Substrate_Height + Slab_Height + WG_Height / 2)
                self.lum.set('output Px', 1)
                self.lum.set('output Py', 1)
                self.lum.set('output Pz', 1)
                self.lum.set('output power', 1)
                # raise ValueError("!!! CAUTION !!! - The Ports are overlapping at the middle! Please change the Y Port Span or move the Waveguides away from each other!")

            else:
                self.lum.addmovie()
                self.lum.set("x min", -(0.5e-6 + Device_Length / 2))
                self.lum.set("x max", (0.5e-6 + Device_Length / 2))
                self.lum.set("y min", -(Device_Width + 1e-6))
                self.lum.set("y max", (Device_Width + 1e-6))
                self.lum.set("z", Ports_mid)

                PortCorrection = [-0.1e-6, 0.1e-6, 0.1e-6]

                for i in range(3):

                    # Power Monitor Port 1
                    self.lum.addpower()
                    self.lum.set('name', "Power_"+ name[i])
                    self.lum.set('monitor type', '2D X-normal')
                    self.lum.set("x", x[i] + PortCorrection[i])
                    self.lum.set("y", yPort_vec[i])
                    self.lum.set("y span", y_Port_Span)
                    self.lum.set("z", Ports_mid)
                    self.lum.set("z span", z_Port_Span)
                    self.lum.set('output Px', 1)
                    self.lum.set('output Py', 1)
                    self.lum.set('output Pz', 1)
                    self.lum.set('output power', 1)


                    self.lum.addport()
                    self.lum.set('name', name[i])
                    self.lum.set("x", x[i])
                    self.lum.set("y", yPort_vec[i])
                    self.lum.set("y span", y_Port_Span)
                    self.lum.set("z", Ports_mid)
                    self.lum.set("z span", z_Port_Span)
                    self.lum.set('direction', direction[i])
                    self.lum.set('mode selection', Mode)



                # Add Power and Freq Monitor
                self.lum.addpower()
                self.lum.set('monitor type', '2D Z-normal')
                self.lum.set("x min", -(0.5e-6 + Device_Length / 2))
                self.lum.set("x max", (0.5e-6 + Device_Length / 2))
                self.lum.set("y min", -(Device_Width + 1e-6))
                self.lum.set("y max", (Device_Width + 1e-6))
                self.lum.set('z', Substrate_Height + Slab_Height + WG_Height / 2)
                self.lum.set('output Px', 1)
                self.lum.set('output Py', 1)
                self.lum.set('output Pz', 1)
                self.lum.set('output power', 1)


        else:
            raise ValueError("Incorect Taper variable!. Possible Taper values are Taper = False or Taper = True.")



    def setDCFDTDSolver(self, Parameters):
        '''

        Parameters
        ----------
        Parameters : Dictionary
            Dictionary with all the data needet for the DCFDTDSolver.
            Parameters['Substrate Height'] : float/int
                Height of the Substrate
            Parameters['Substrate Width'] : float/int
                Width of the MMI
            Parameters['DC Length'] : float/int
                Length of the Directional coupler
            Parameters['WG Height'] : float/int
                Height of the Waveguide
            Parameters['WG Width'] : float/int
                Waveguide Width
            Parameters['Position Offset'] : float/int
                Positional offser of the waveguides. If posOffset the two Waveguides
                will be offset of the middle position (y = 0) by the half of there
                Width. In this case they will not overlap if the Offset is 0.
            Parameters['x res'] : float/int
                Mesh resolution for the x-Axis
            Parameters['Slab Height'] : float/int
                Slab height.
            Parameters['Wavelength'] : float/int
                Wavelength
            Parameters["Port Span"]: list of int/floats
                Parameters["Port Span"] = [Span of Port in x direction, Span of Port in y direction, Span of Port in z direction]
            Parameters["Mode"] : str
                Mode to choose from ("fundamental TE mode", "fundamental TM mode", "fundamental mode")

        Returns
        -------
        None.

        '''

        Substrate_Height = Parameters['Substrate Height']
        Substrate_Width = Parameters['Substrate Width']
        DC_Lenght = Parameters['DC Length']
        WG_Height = Parameters['WG Height']
        WG_Width = Parameters['WG Width']
        posOffset = Parameters['Position Offset']
        x_res = Parameters['x res']
        Slab_Height = Parameters['Slab Height']
        WaveLength = Parameters['Wavelength']
        Mode = Parameters["Mode"]
        y_Port_Span = Parameters["Port Span"][1]
        z_Port_Span = Parameters["Port Span"][2]



        # Device specifications
        Device_Length = DC_Lenght
        Device_Width = Substrate_Width
        max_slabH = Slab_Height
        Ports_mid = Substrate_Height + max_slabH + WG_Height/2

        # Adds a Eigenmode Expansion (EME) solver region to the MODE simulation environment.
        self.lum.addfdtd()
        self.lum.set("x min", -(Device_Length / 2))
        self.lum.set("x max", (Device_Length / 2))
        self.lum.set("y", 0)
        self.lum.set('simulation temperature', 273.15 + 20)
        self.lum.set("y span", Device_Width)
        self.lum.set("z", Substrate_Height)
        self.lum.set("z span", 4e-6)
        self.lum.set('z min bc', 'PML')
        self.lum.set('z max bc', 'PML')
        self.lum.set('mesh type', 'auto non-uniform')
        self.lum.set('min mesh step', x_res)
        self.lum.set('set simulation bandwidth', 0)
        self.lum.set('global source center wavelength', WaveLength)
        self.lum.set('global source wavelength span', 0)

        x = [(Device_Length / 2), (Device_Length / 2), -(Device_Length / 2), -(Device_Length / 2)]
        direction = ['Backward', 'Backward', 'Forward', 'Forward']
        yPort_vec = [posOffset / 2 + WG_Width / 2, -(posOffset / 2 + WG_Width / 2), posOffset / 2 + WG_Width / 2, -(posOffset / 2 + WG_Width / 2)]
        name = ['Input_L', 'Input_R', 'Output_L', 'Output_R']

        overLapp = yPort_vec[0] - y_Port_Span/2
        if overLapp < 0:
            raise ValueError("!!! CAUTION !!! - The Ports are overlapping at the middle! Please change the Y Port Span or move the Waveguides away from each other!")
        else:
            pass

        self.lum.addmovie()
        self.lum.set("x min", -(0.5e-6 + Device_Length / 2))
        self.lum.set("x max", (0.5e-6 + Device_Length / 2))
        self.lum.set("y min", -(Device_Width + 1e-6))
        self.lum.set("y max", (Device_Width + 1e-6))
        self.lum.set("z", Ports_mid)

        for i in range(4):

            # Power Monitor Port 1
            self.lum.addpower()
            self.lum.set('name', "Power_"+ name[i])
            self.lum.set('monitor type', '2D X-normal')
            self.lum.set("x", x[i]+0.1e-6)
            self.lum.set("y", yPort_vec[i])
            self.lum.set("y span", y_Port_Span)
            self.lum.set("z", Ports_mid)
            self.lum.set("z span", z_Port_Span)
            self.lum.set('output Px', 1)
            self.lum.set('output Py', 1)
            self.lum.set('output Pz', 1)
            self.lum.set('output power', 1)


            self.lum.addport()
            self.lum.set('name', name[i])
            self.lum.set("x", x[i])
            self.lum.set('y', yPort_vec[i])
            self.lum.set('y span', y_Port_Span)
            self.lum.set("z", Ports_mid)
            self.lum.set("z span", z_Port_Span)
            self.lum.set('direction', direction[i])
            self.lum.set('mode selection', Mode)

        # Add Power and Freq Monitor
        self.lum.addpower()
        self.lum.set('monitor type', '2D Z-normal')
        self.lum.set("x min", -(0.5e-6 + Device_Length / 2))
        self.lum.set("x max", (0.5e-6 + Device_Length / 2))
        self.lum.set("y min", -(Device_Width + 1e-6))
        self.lum.set("y max", (Device_Width + 1e-6))
        self.lum.set('z', Substrate_Height + Slab_Height + WG_Height / 2)
        self.lum.set('output Px', 1)
        self.lum.set('output Py', 1)
        self.lum.set('output Pz', 1)
        self.lum.set('output power', 1)





    # def setWDMFDTDSolver(self, Parameters):
    #     '''


    #     Parameters
    #     ----------
    #     Parameters : Dictionary
    #         Dictionary with all the data needet for the Bend Wavaguide. Data needet:
    #             Parameters
    #             ----------
    #             Parameters['Substrate Height'] : int/float
    #                 Substrate height.
    #             Parameters['MMI Width'] : int/float
    #                 Width of the MMI.
    #             Parameters['MMI Length'] : int/float
    #                 Length of the MMI.
    #             Parameters['WG Height' : int/float
    #                 Waveguide hight. Also the height of the MMI section
    #             Parameters['WG Width'] : int/float
    #                 Waveguide width.
    #             Parameters['WG Length'] : int/float
    #                 Waveguide length.
    #             Parameters['y res'] : int/float
    #                 Mesh y-Axis
    #             Parameters['z res'] : int/float
    #                 Mesh z-Axis
    #             Parameters['Slab Height'] : int/float
    #                 Slab Height.
    #             Parameters['Wavelength'] : int/float
    #                 Wavelength
    #             Parameters['Angle Thetha'] : boolen
    #                 Angle for the input and output waveguides
    #             Parameters["Mode"] : str
    #                 Mode to choose from ("fundamental TE mode", "fundamental TM mode", "fundamental mode")

    #     Returns
    #     -------
    #     None.

    #     '''

    #     Substrate_Height = Parameters['Substrate Height']
    #     MMI_Width = Parameters['MMI Width']
    #     MMI_Length = Parameters['MMI Length']
    #     WG_Height = Parameters['WG Height']
    #     WG_Width = Parameters['WG Width']
    #     WG_Length = Parameters['WG Length']
    #     y_res = Parameters['y res']
    #     z_res = Parameters['z res']
    #     Slab_Height = Parameters['Slab Height']
    #     WaveLength = Parameters['Wavelength']
    #     angleTheta = Parameters['Angle Thetha']
    #     Mode = Parameters["Mode"]


    #     # Device specifications
    #     Device_Length = MMI_Length + 2 * WG_Length
    #     Device_Width = MMI_Width + 2*(WG_Length) + WaveLength * 2  # MMI_Width
    #     max_slabH = Slab_Height
    #     MonitorHeight = Substrate_Height + (max_slabH + WG_Height) / 2
    #     Ports_mid = (max_slabH + WG_Height) / 2

    #     # Adds a Eigenmode Expansion (EME) solver region to the MODE simulation environment.
    #     self.lum.addeme()
    #     self.lum.set('simulation temperature', 273.15 + 20)
    #     self.lum.set("x min", 1e-6)
    #     self.lum.set("y", 0)
    #     self.lum.set("y span", Device_Width)
    #     self.lum.set("z", Substrate_Height)
    #     self.lum.set("z span", 4e-6)
    #     self.lum.set("wavelength", WaveLength)
    #     self.lum.set("z min bc", "PML")
    #     self.lum.set("z max bc", "PML")
    #     self.lum.set("y min bc", "PML")
    #     self.lum.set("y max bc", "PML")
    #     # set cell properties
    #     self.lum.set("number of cell groups", 3)
    #     self.lum.set("group spans", np.array([[WG_Length - 1e-6], [MMI_Length], [WG_Length - 1e-6]]))
    #     self.lum.set("cells", np.array([[3], [3], [3]]))
    #     self.lum.set("subcell method", np.array([[0], [0], [0]]))

    #     # Modes to Calculate
    #     self.lum.set('number of modes for all cell groups', 20)

    #     # Mesh Cells
    #     self.lum.set("define y mesh by", "maximum mesh step")
    #     self.lum.set("dy", y_res)
    #     self.lum.set("define z mesh by", "maximum mesh step")
    #     self.lum.set("dz", z_res)
    #     self.lum.set('fit materials with multi-coefficient model', 1)
    #     self.lum.set('wavelength start', 0.4e-6)
    #     self.lum.set('wavelength stop', 2e-6)

    #     # max_yPos = [(WG_W / 2 + OffsetInput) * 2, 0, (WG_Width / 2 + posOffset / 2) * 2]
    #     # min_yPos = [-(WG_W / 2 + OffsetInput) * 2, -(WG_Width / 2 + posOffset / 2) * 2, 0]

    #     Input_yPos = -MMI_Width / 2 + WG_Width / 2
    #     y = WG_Length * np.tan(angleTheta * np.pi / 180)
    #     Input_Y = Input_yPos - y + WG_Width

    #     Output_yPos = MMI_Width / 2 - WG_Width / 2
    #     y2 = WG_Length * np.tan(angleTheta * np.pi / 180)
    #     Output_Y = Output_yPos + y2 - WG_Width / 2
    #     portLoc = ["left", "right"]

    #     self.lum.select("EME::Ports::port_" + str(1))
    #     self.lum.set("port location", portLoc[0])
    #     self.lum.set("use full simulation span", 0)
    #     # self.lum.set("y min", min_yPos[1])
    #     # self.lum.set("y max", max_yPos[1])
    #     self.lum.set("y", Input_Y)
    #     self.lum.set("y span", WG_Width + 2e-6)
    #     self.lum.set("z", Ports_mid)
    #     self.lum.set("z span", 2e-6)
    #     self.lum.set("mode selection", Mode)
    #     self.lum.set("theta", angleTheta)

    #     self.lum.select("EME::Ports::port_" + str(2))
    #     self.lum.set("port location", portLoc[1])
    #     self.lum.set("use full simulation span", 0)
    #     # self.lum.set("y min", min_yPos[2])
    #     # self.lum.set("y max", max_yPos[2])
    #     self.lum.set("y", Output_Y)
    #     self.lum.set("y span", WG_Width + 2e-6)
    #     self.lum.set("z", Ports_mid)
    #     self.lum.set("z span", 2e-6)
    #     self.lum.set("mode selection", Mode)
    #     self.lum.set("theta", angleTheta)

    #     # Add monitor
    #     # x_MMI = Device_Length / 2
    #     self.lum.addemeprofile()
    #     self.lum.set("x", Device_Length / 2)
    #     self.lum.set("x span", Device_Length)
    #     # self.lum.set("x min", -(Device_Length / 2))
    #     # self.lum.set("x max", (Device_Length / 2))
    #     self.lum.set("y", 0)
    #     self.lum.set("y span", Device_Width )
    #     # self.lum.set("y min", -Device_Width / 2)
    #     # self.lum.set("y max", Device_Width / 2)
    #     self.lum.set("z", MonitorHeight)
    
    
    
    def setInverseTaperFDTDSolver(self, Parameters):
        '''
          Parameters
          ----------
          Parameters : Dictionary
              Dictionary with all the data needet for the InverseTaperFDTDSolver.
              Parameters['Substrate Height'] : int/float
                  Substrate height.
              Parameters['WG Height' : int/float
                  Waveguide hight. Also the height of the MMI section
              Parameters['WG Width'] : int/float
                  Waveguide width.
              Parameters['Slab Height'] : int/float
                  Slab height
              Parameters['PWB Taper Width Back'] : int/float
                  Photonic Wirebonding (PWB) Width back side (to the Photonic Wire Bonding)
              Parameters['PWB Taper Hight Back'] : int/float
                  Photonic Wire Bonding Height back side (to the Photonic Wire Bonding)
              Parameters['PWB Taper Length'] : int/float
                  Length of the Photonic Wire Bonding Taper
              Parameters["SMF Core Diameter"] : int/float
                Single Mode Fiber core Diameter
              Parameters["SMF Cladding Diameter"] : int/float
                Single Mode Fiber Cladding Diameter
              Parameters['x res'] : int/float
                  Mesh x-Axis
              Parameters['Wavelength'] : int/float
                  Wavelength
              Parameters["Mode"] : str
                  Mode to choose from ("fundamental TE mode", "fundamental TM mode", "fundamental mode")
              Parameters["Port Span"] : list of floats/ints
                  List of x,y and z span of the Ports. For this simulation only y and z parametes will be taken.

        Returns
        -------
        None.

        '''

        Substrate_Height = Parameters['Substrate Height']
        WG_Height = Parameters['WG Height']
        WG_Width = Parameters['WG Width']
        WG_Length = Parameters['WG Length']
        Slab_Height = Parameters['Slab Height']
        TaperWidthB = Parameters['PWB Taper Width Back']
        TaperHightB = Parameters['PWB Taper Hight Back']
        TaperLength_PWB = Parameters['PWB Taper Length']
        WaveLength = Parameters['Wavelength']
        Mode = Parameters["Mode"]
        x_res = Parameters['x res']
        y_Port_Span = Parameters["Port Span"][1]
        z_Port_Span = Parameters["Port Span"][2]
        
        # SMF Parameters
        CoreDiameter = Parameters["SMF Core Diameter"]
        CladdingDiameter = Parameters["SMF Cladding Diameter"]


        if Slab_Height == 0:
            # # Device specifications

            Ports_mid = Substrate_Height/2 + WG_Height/2
            Ports_PWB_mid = Substrate_Height + CoreDiameter/2 # TaperHightB/2
            self.lum.select("SMF")
            xPos_SMF = self.lum.get("x")
            X_min = -TaperLength_PWB/2 - abs(xPos_SMF)
            
        else:
            # # Device specifications
            max_slabH = Slab_Height

            Ports_mid = max_slabH + Substrate_Height/2  + WG_Height/2
            Ports_PWB_mid = max_slabH + CoreDiameter/2 # TaperHightB/2
            self.lum.select("SMF")
            xPos_SMF = self.lum.get("x")
            X_min = -TaperLength_PWB/2 - abs(xPos_SMF)
            
            
        # Adds a Eigenmode Expansion (EME) solver region to the MODE simulation environment.
        self.lum.addfdtd()
        self.lum.set("x min", X_min)
        self.lum.set("x max", TaperLength_PWB/2 + WG_Length)
        self.lum.set("y", 0)
        self.lum.set("y span", 2*TaperWidthB )# + TaperWidthB/2
        self.lum.set('simulation temperature', 273.15 + 20)
        self.lum.set("z", Substrate_Height/2 + CoreDiameter/2 ) #Substrate_Height
        self.lum.set("z span", TaperHightB*2)
        self.lum.set('x min bc', 'PML')
        self.lum.set('x max bc', 'PML')
        self.lum.set('y min bc', 'Anti-Symmetric')
        self.lum.set('y max bc', 'PML')
        self.lum.set('z min bc', 'PML')
        self.lum.set('z max bc', 'PML')
        self.lum.set('mesh type', 'auto non-uniform')
        self.lum.set('min mesh step', x_res)
        self.lum.set('set simulation bandwidth', 0)
        self.lum.set('global source center wavelength', WaveLength)
        self.lum.set('global source wavelength span', 0)
        
        
        # Ports Positions
        xPort_Pos = [X_min, TaperLength_PWB]
        yPort_Pos = [CoreDiameter + 2e-6, y_Port_Span]
        zPort_Pos = [CoreDiameter+2e-6, z_Port_Span]
        direction = ['Backward', 'Forward', 'Forward', 'Forward', 'Forward']
        name = ['SMF Port', 'Waveguide Port']
        

   
        x =[X_min, TaperLength_PWB/2 + WG_Length]
        PortCorrection = [0.1e-6, -0.1e-6]
        zMid_Pos_Ports = [Ports_PWB_mid, Ports_mid]
        direction = ["Forward", "Backward"]
        
        # Gaussian Source
        self.lum.addgaussian()
        self.lum.set("injection axis", "x-axis")
        self.lum.set("direction", direction[0])
        self.lum.set("x", x[0])
        self.lum.set("y", 0)
        self.lum.set("y span",  2*TaperWidthB) # yPort_Pos[0]
        self.lum.set("z", zMid_Pos_Ports[0])
        self.lum.set("z span", TaperHightB*2) # zPort_Pos[0]
        self.lum.set("waist radius w0", CoreDiameter/2)
        self.lum. set("distance from waist",0)
    

        # Add Ports to Structure and power monitors
        for i in range(2):
            self.lum.addpower()
            self.lum.set('name', "Power_" + name[i])
            self.lum.set('monitor type', '2D X-normal')
            self.lum.set("x", x[i] + PortCorrection[i])
            self.lum.set("y", 0)
            self.lum.set("y span", yPort_Pos[i])
            self.lum.set("z", zMid_Pos_Ports[i])
            self.lum.set("z span", zPort_Pos[i])
            self.lum.set('output Px', 1)
            self.lum.set('output Py', 1)
            self.lum.set('output Pz', 1)
            self.lum.set('output power', 1)

            # self.lum.addport()
            # self.lum.set('name', name[i])
            # self.lum.set("x", x[i])
            # self.lum.set("y", 0)
            # self.lum.set("y span", yPort_Pos[i])
            # self.lum.set("z", zMid_Pos_Ports[i])
            # self.lum.set("z span", zPort_Pos[i])
            # self.lum.set('direction', direction[i])
            # self.lum.set('mode selection', Mode)
            

        
        # Add Z Monitor over stucute
        self.lum.addpower()
        self.lum.set('name', "Power 2D Z-Normal")
        self.lum.set('monitor type', '2D Z-normal')
        self.lum.set("x min", X_min)
        self.lum.set("x max", TaperLength_PWB/2 + WG_Length)
        self.lum.set("y", 0)
        self.lum.set("y span", TaperWidthB + TaperWidthB/2)
        self.lum.set("z", Ports_mid) #Substrate_Height
        self.lum.set('output Px', 1)
        self.lum.set('output Py', 1)
        self.lum.set('output Pz', 1)
        self.lum.set('output power', 1)
        
        # Add Movi Monitor Z-Normal over structure
        self.lum.addmovie()
        self.lum.set("x min", X_min)
        self.lum.set("x max", TaperLength_PWB/2 + WG_Length)
        self.lum.set("y", 0)
        self.lum.set("y span", TaperWidthB + TaperWidthB/2)
        self.lum.set("z", Ports_mid)
        
        
        # Add Y Monitor over stucute
        self.lum.addpower()
        self.lum.set('name', "Power 2D Y-Normal")
        self.lum.set('monitor type', '2D Y-normal')
        self.lum.set("x min", X_min)
        self.lum.set("x max", TaperLength_PWB/2 + WG_Length)
        self.lum.set("y", 0)
        self.lum.set("z", Substrate_Height/2 + CoreDiameter/2 ) #Substrate_Height
        self.lum.set("z span", TaperHightB*2)
        self.lum.set('output Px', 1)
        self.lum.set('output Py', 1)
        self.lum.set('output Pz', 1)
        self.lum.set('output power', 1)
        
        # #Add Extra Mesch
        # self.lum.addmesh()
        # self.lum.set("name", "Mesh PWB")
        # self.lum.set("based on a structure",1)
        # self.lum.set("structure", "Taper_PWB")
        # self.lum.set("set maximum mesh step",1)
        # self.lum.set("override x mesh",0)
        # # self.lum.set("dx", x_res)
        # self.lum.set("override y mesh",1)
        # self.lum.set("dy", x_res)
        # self.lum.set("override z mesh",1)
        # self.lum.set("dz", x_res)
        
        # self.lum.addmesh()
        # self.lum.set("name", "Mesh Inverse Taper")
        # self.lum.set("based on a structure",1)
        # self.lum.set("structure", "InverseTaper")
        # self.lum.set("set maximum mesh step",1)
        # self.lum.set("override x mesh",0)
        # # self.lum.set("dx", x_res)
        # self.lum.set("override y mesh",1)
        # self.lum.set("dy", x_res)
        # self.lum.set("override z mesh",1)
        # self.lum.set("dz", x_res)
        
        
    
    
    
    def setCascadetMMIFDTDSolver(self, Parameters, SpaceX, SpaceY):
        
        
            Substrate_Height = Parameters['Substrate Height']
            MMI_Width = Parameters['MMI Width']
            MMI_Length = Parameters['MMI Length']
    
            WG_Height = Parameters['WG Height']
            WG_Width = Parameters['WG Width']
            WG_Length = Parameters['WG Length']
            OffsetInput = Parameters['Offset Input']
            posOffset = Parameters['Position Offset']
            x_res = Parameters['x res']
            Slab_Height = Parameters['Slab Height']
            WaveLength = Parameters['Wavelength']
            TaperLength = Parameters['Taper Length']
            Taper = Parameters['Taper']
            Mode = Parameters["Mode"]
            y_Port_Span = Parameters["Port Span"][1]
            z_Port_Span = Parameters["Port Span"][2]
            
            
            
            if Taper == False:
                
        
        
                # Device specifications
                Device_Length = 2*MMI_Length + 4 * WG_Length + SpaceX
                Device_Width = 3*MMI_Width + WaveLength * 2  + 2*SpaceY# MMI_Width
                max_slabH = Slab_Height
                # Ports_mid = Substrate_Height + (max_slabH + WG_Height) / 2
                # Ports_mid = Substrate_Height + Slab_Height + WG_Height / 2
                # MonitorHeight = Substrate_Height + (max_slabH + WG_Height) / 2
                Ports_mid = Substrate_Height + max_slabH
                # WG_H = WG_Height
            
                x_Offset = MMI_Length + WG_Length
                # Adds a Eigenmode Expansion (EME) solver region to the MODE simulation environment.
                self.lum.addfdtd()
                self.lum.set("x", -x_Offset)
                self.lum.set("x span", Device_Length)
                self.lum.set("y", 0)
                self.lum.set('simulation temperature', 273.15 + 20)
                self.lum.set("y span", Device_Width)
                self.lum.set("z", Substrate_Height)
                self.lum.set("z span", 4e-6)
                self.lum.set('z min bc', 'PML')
                self.lum.set('z max bc', 'PML')
                self.lum.set('mesh type', 'auto non-uniform')
                self.lum.set('min mesh step', x_res)
                self.lum.set('set simulation bandwidth', 0)
                self.lum.set('global source center wavelength', WaveLength)
                self.lum.set('global source wavelength span', 0)
            
                # Positions of the Input and Output WGs
                # Triangle EQ for MMI Width

                yPos = [0 + OffsetInput, WG_Width / 2 + posOffset / 2, - WG_Width / 2 - posOffset / 2]
                MMI_Wid = MMI_Width 
                Parameters__Position_X = [0, -SpaceX -MMI_Length  - 2*WG_Length , -SpaceX -MMI_Length  - 2*WG_Length]
                Parameters__Position_Y = [0, (MMI_Wid + SpaceY),  -(MMI_Wid + SpaceY)]
                y_span = (Parameters__Position_Y[1] + yPos[0]) - (Parameters__Position_Y[0] + yPos[1] ) 
                x = [(MMI_Length/2 + WG_Length )  , -(MMI_Length/2 + MMI_Length+ 3*WG_Length + SpaceX ), -(MMI_Length/2 + MMI_Length+ 3*WG_Length + SpaceX),  -(MMI_Length/2 + MMI_Length+ 3*WG_Length + SpaceX ), -(MMI_Length/2 + MMI_Length+ 3*WG_Length + SpaceX )]
                direction = ['Backward', 'Forward', 'Forward', 'Forward', 'Forward']
            
                name = ['Input', 'Output1_Top', 'Output1_Bot', 'Output2_Top', 'Output2_Bot']
                
                yPort_vec = [OffsetInput,  y_span + WG_Width + posOffset +  WG_Width/2,  y_span + WG_Width / 2, -y_span - WG_Width / 2, - y_span - WG_Width- posOffset -  WG_Width/2]
                
                overLapp = yPort_vec[2] - y_Port_Span/2
                # if overLapp < 0:
                #     y_Port_Span_old = y_Port_Span
                #     y_Port_Span = yPort_vec[2]*2
                #     print((f"!!! CAUTION !!! - The Ports are overlapping at the middle! So the y Span of the ports will be reduze from y_span = {y_Port_Span_old} to y_span_new = {y_Port_Span}"))

                self.lum.addmovie()
                self.lum.set("x", -x_Offset)
                self.lum.set("x span", Device_Length)
                self.lum.set("y min", -(Device_Width + 1e-6))
                self.lum.set("y max", (Device_Width + 1e-6))
                self.lum.set("z", Ports_mid)

                PortCorrection = [-0.1e-6, 0.1e-6, 0.1e-6,0.1e-6, 0.1e-6]

                for i in range(5):
                    # Power Monitor Port 1
                    self.lum.addpower()
                    self.lum.set('name', "Power_" + name[i])
                    self.lum.set('monitor type', '2D X-normal')
                    self.lum.set("x", x[i] + PortCorrection[i])
                    self.lum.set("y", yPort_vec[i])
                    self.lum.set("y span", y_Port_Span)
                    self.lum.set("z", Ports_mid)
                    self.lum.set("z span", z_Port_Span)
                    self.lum.set('output Px', 1)
                    self.lum.set('output Py', 1)
                    self.lum.set('output Pz', 1)
                    self.lum.set('output power', 1)

                    self.lum.addport()
                    self.lum.set('name', name[i])
                    self.lum.set("x", x[i])
                    self.lum.set("y", yPort_vec[i])
                    self.lum.set("y span", y_Port_Span)
                    self.lum.set("z", Ports_mid)
                    self.lum.set("z span", z_Port_Span)
                    self.lum.set('direction', direction[i])
                    self.lum.set('mode selection', Mode)


                # Add Power and Freq Monitor
                self.lum.addpower()
                self.lum.set('monitor type', '2D Z-normal')
                self.lum.set("x", -(0.5e-6 + x_Offset))
                self.lum.set("x span", (0.5e-6 + Device_Length))
                self.lum.set("y min", -(Device_Width + 1e-6))
                self.lum.set("y max", (Device_Width + 1e-6))
                self.lum.set('z', Substrate_Height + Slab_Height + WG_Height / 2)
                self.lum.set('output Px', 1)
                self.lum.set('output Py', 1)
                self.lum.set('output Pz', 1)
                self.lum.set('output power', 1)




                # Add Power and Freq Monitor
                self.lum.addpower()
                self.lum.set('monitor type', '2D Z-normal')
                self.lum.set("x", -(0.5e-6 + x_Offset))
                self.lum.set("x span", (0.5e-6 + Device_Length ))
                self.lum.set("y min", -(Device_Width + 1e-6))
                self.lum.set("y max", (Device_Width + 1e-6))
                self.lum.set('z', Substrate_Height + Slab_Height + WG_Height / 2)
                self.lum.set('output Px', 1)
                self.lum.set('output Py', 1)
                self.lum.set('output Pz', 1)
                self.lum.set('output power', 1)



            
            elif Taper == True:
                
     
                # Device specifications
                Device_Length = 2*MMI_Length +4 * WG_Length + 4*TaperLength + SpaceX
                Device_Width = 3*MMI_Width + WaveLength * 2  + 2*SpaceY# MMI_Width
                max_slabH = Slab_Height
                # Ports_mid = Substrate_Height + (max_slabH + WG_Height) / 2
                Ports_mid = Substrate_Height + max_slabH
                # MonitorHeight = Substrate_Height + (max_slabH + WG_Height) / 2
                # WG_H = WG_Height
                x_Offset = MMI_Length + WG_Length + TaperLength
                # Adds a Eigenmode Expansion (EME) solver region to the MODE simulation environment.
                self.lum.addfdtd()
                self.lum.set("x", -x_Offset)
                self.lum.set("x span", Device_Length)
                self.lum.set("y", 0)
                self.lum.set('simulation temperature', 273.15 + 20)
                self.lum.set("y span", Device_Width)
                self.lum.set("z", Substrate_Height)
                self.lum.set("z span", 4e-6)
                self.lum.set('z min bc', 'PML')
                self.lum.set('z max bc', 'PML')
                self.lum.set('mesh type', 'auto non-uniform')
                self.lum.set('min mesh step', x_res)
                self.lum.set('set simulation bandwidth', 0);
                self.lum.set('global source center wavelength', WaveLength)
                self.lum.set('global source wavelength span', 0)
            
                # Positions of the Input and Output WGs
                # Triangle EQ for MMI Width
                yPos = [0 + OffsetInput, WG_Width / 2 + posOffset / 2, - WG_Width / 2 - posOffset / 2]
                MMI_Wid = MMI_Width 
                Parameters__Position_X = [0, -SpaceX -MMI_Length  - 2*WG_Length , -SpaceX -MMI_Length  - 2*WG_Length]
                Parameters__Position_Y = [0, (MMI_Wid + SpaceY),  -(MMI_Wid + SpaceY)]
                y_span = (Parameters__Position_Y[1] + yPos[0]) - (Parameters__Position_Y[0] + yPos[1] ) 
                x = [(MMI_Length/2 + WG_Length + TaperLength )  , -(MMI_Length/2 + MMI_Length+ 3*WG_Length + SpaceX+ 3*TaperLength ), -(MMI_Length/2 + MMI_Length+ 3*WG_Length + SpaceX+ 3*TaperLength),  -(MMI_Length/2 + MMI_Length+ 3*WG_Length + SpaceX+ 3*TaperLength ), -(MMI_Length/2 + MMI_Length+ 3*WG_Length + SpaceX+ 3*TaperLength )]
                direction = ['Backward', 'Forward', 'Forward', 'Forward', 'Forward']
            
                name = ['Input', 'Output1_Top', 'Output1_Bot', 'Output2_Top', 'Output2_Bot']
                
                yPort_vec = [OffsetInput,  y_span + WG_Width + posOffset +  WG_Width/2,  y_span + WG_Width / 2, -y_span - WG_Width / 2, - y_span - WG_Width- posOffset -  WG_Width/2]
                
                
                overLapp = yPort_vec[2] - y_Port_Span/2
                self.lum.addmovie()
                self.lum.set("x", -x_Offset)
                self.lum.set("x span", Device_Length)
                self.lum.set("y min", -(Device_Width + 1e-6))
                self.lum.set("y max", (Device_Width + 1e-6))
                self.lum.set("z", Ports_mid)

                PortCorrection = [-0.1e-6, 0.1e-6, 0.1e-6,0.1e-6, 0.1e-6]

                for i in range(5):
                    # Power Monitor Port 1
                    self.lum.addpower()
                    self.lum.set('name', "Power_" + name[i])
                    self.lum.set('monitor type', '2D X-normal')
                    self.lum.set("x", x[i] + PortCorrection[i])
                    self.lum.set("y", yPort_vec[i])
                    self.lum.set("y span", y_Port_Span)
                    self.lum.set("z", Ports_mid)
                    self.lum.set("z span", z_Port_Span)
                    self.lum.set('output Px', 1)
                    self.lum.set('output Py', 1)
                    self.lum.set('output Pz', 1)
                    self.lum.set('output power', 1)

                    self.lum.addport()
                    self.lum.set('name', name[i])
                    self.lum.set("x", x[i])
                    self.lum.set("y", yPort_vec[i])
                    self.lum.set("y span", y_Port_Span)
                    self.lum.set("z", Ports_mid)
                    self.lum.set("z span", z_Port_Span)
                    self.lum.set('direction', direction[i])
                    self.lum.set('mode selection', Mode)


                # Add Power and Freq Monitor
                self.lum.addpower()
                self.lum.set('monitor type', '2D Z-normal')
                self.lum.set("x", -(0.5e-6 + x_Offset))
                self.lum.set("x span", (0.5e-6 + Device_Length))
                self.lum.set("y min", -(Device_Width + 1e-6))
                self.lum.set("y max", (Device_Width + 1e-6))
                self.lum.set('z', Substrate_Height + Slab_Height + WG_Height / 2)
                self.lum.set('output Px', 1)
                self.lum.set('output Py', 1)
                self.lum.set('output Pz', 1)
                self.lum.set('output power', 1)




                # Add Power and Freq Monitor
                self.lum.addpower()
                self.lum.set('monitor type', '2D Z-normal')
                self.lum.set("x", -(0.5e-6 + x_Offset))
                self.lum.set("x span", (0.5e-6 + Device_Length ))
                self.lum.set("y min", -(Device_Width + 1e-6))
                self.lum.set("y max", (Device_Width + 1e-6))
                self.lum.set('z', Substrate_Height + Slab_Height + WG_Height / 2)
                self.lum.set('output Px', 1)
                self.lum.set('output Py', 1)
                self.lum.set('output Pz', 1)
                self.lum.set('output power', 1)


              

            else:
                raise ValueError("Incorect Taper variable!. Possible Taper values are Taper = False or Taper = True.")



    def setGratingCouplerFDTDSolver(self, Parameters):
        '''
       
        Parameters
        ----------
        Parameters : Dictionary
            Dictionary with all the data needet for the GratingCouplerFDTDSolver.
            Parameters['Substrate Height'] : int/float
                Substrate height.
            Parameters["Length GC"]: int/float
                Lenght of the Grating Coupler Area
            Parameters["Input Length GC"]: int/float
                An squere Waveguide with the same WG Height as the Grating coupler place before the Grating Coupler region will start. 
            Parameters["Output Length GC"]: int/float
                An squere Waveguide with the same WG Height as the Grating coupler place after the Grating Coupler region to finish the structure.
            Parameters["Width GC"]: int/float
                Widht of the Grating Coupler Area
            Parameters["Hight GC"]: int/float
                Hight of the Grating Coupler Material
            Parameters['Taper'] : boolen
                Add Taper to structure
            Parameters['Taper Length'] : int/float
                  Length of the Taper
            Parameters['Wavelength'] : int/float
                  Wavelength
            Parameters['x res'] : int/float
                  Mesh x-Axis
            Parameters["Mode"] : str
                Mode to choose from ("fundamental TE mode", "fundamental TM mode", "fundamental mode")
            Parameters["Port Span"] : list of floats/ints
                  List of x,y and z span of the Ports. For this simulation only y and z parametes will be taken.
            Parameters["GC Radius"]: int/float
                Radius of the Ring Grating Coupler in um. For Example "Parameters["GC Radius"] = 25e-6"
       
            

        Returns
        -------
        None.

        '''
        
       

        SubstrateThickness = Parameters['Substrate Height']
        GC_SectionLenght = Parameters["Length GC"]
        InputLenght = Parameters["Input Length GC"]
        OutputLenght = Parameters["Output Length GC"]
        WidthGC = Parameters["Width GC"]
        Hight = Parameters["Hight GC"]
        Taper = Parameters["Taper"]
        TaperLength = Parameters['Taper Length']
        ZSpan = Parameters["SMF Z Span"]
        Theta = Parameters["SMF Theta"]
        CoreDiameter = Parameters["SMF Core Diameter"]
        x_res = Parameters['x res']
        WaveLength = Parameters['Wavelength']
        Mode = Parameters["Mode"]
        y_Port_Span = Parameters["Port Span"][1]
        z_Port_Span = Parameters["Port Span"][2]
        CladdingThickness = 0.7e-6
        

        # Device specifications
        Device_Length = GC_SectionLenght + InputLenght + OutputLenght + TaperLength
        FDTD_ZSpan = Hight + CladdingThickness + SubstrateThickness


        # Define Ports
        Port_Names = ["Input_SMF_Port", "Output"]
       
        
        self.lum.select("SMF")
        fiber_xpos = self.lum.get("x")


        if Taper == True:

            # Adds a Finite-Difference Time-Domain  (FDTD) solver region to the MODE simulation environment.
            self.lum.addfdtd()
            self.lum.set("x min", -GC_SectionLenght/2 - InputLenght - TaperLength)
            self.lum.set("x max",  fiber_xpos + CoreDiameter )
            # self.lum.set("x", -TaperLength / 2)
            # self.lum.set("x span", Device_Length)
            self.lum.set("y", 0)
            self.lum.set("y span", WidthGC)
            self.lum.set("z", 0)
            self.lum.set("z span",  (ZSpan / 2))
            self.lum.set('simulation temperature', 273.15 + 20)
            self.lum.set('z min bc', 'PML')
            self.lum.set('z max bc', 'PML')
            self.lum.set('y min bc', 'Anti-Symmetric')
            self.lum.set('y max bc', 'PML')
            self.lum.set('mesh type', 'auto non-uniform')
            self.lum.set('min mesh step', x_res)
            self.lum.set('set simulation bandwidth', 0)
            self.lum.set('global source center wavelength', WaveLength)
            self.lum.set('global source wavelength span', 0)
            
            
            # Detect Fiber position for Port exact aligment
            self.lum.select("SMF")
            fiber_xpos = self.lum.get("x")
            fiber_ypos = self.lum.get("y")
            fiber_zpos = self.lum.get("z")
            self.lum.select("SMF::core")
            fiber_core_diameter = 2 * self.lum.get("radius")
            fiber_core_index = self.lum.get("index")
            fiber_theta = Theta

    

            # Faser Port
            self.lum.addport()
            self.lum.set('name', Port_Names[0])
            self.lum.set('injection axis', "z-axis")
            self.lum.set('direction', "Backward")
            self.lum.set('mode selection', Mode)
            self.lum.set('theta', Theta)
            self.lum.set("x", fiber_xpos)
            self.lum.set('x span', CoreDiameter + CoreDiameter / 2)
            # self.lum.set("x", 0)
            # self.lum.set('x span', CoreDiameter + CoreDiameter / 2)
            self.lum.set('y', 0)
            self.lum.set('y span', CoreDiameter + CoreDiameter / 2)
            self.lum.set("rotation offset", ZSpan / 4)
            self.lum.set("z", (SubstrateThickness + 2e-6) / 2 )


            # Output Port
            self.lum.addport()
            self.lum.set('name', Port_Names[1])
            self.lum.set("x", -GC_SectionLenght/2 - InputLenght - TaperLength)
            self.lum.set('x span', CoreDiameter)
            self.lum.set('y', 0)
            self.lum.set('y span', y_Port_Span)
            self.lum.set("z", Hight/2 )
            self.lum.set("z span", z_Port_Span )
            self.lum.set('injection axis', "x-axis")
            self.lum.set('direction', "Forward")
            self.lum.set('mode selection', Mode)


            # Power Monitor SMF Port
            self.lum.addpower()
            self.lum.set('name', "Power_"+ Port_Names[0])
            self.lum.set('monitor type', '2D Z-normal')
            # self.lum.set("x", 0)
            # self.lum.set('x span', CoreDiameter + CoreDiameter / 2)
            self.lum.set("x", fiber_xpos)
            self.lum.set('x span', CoreDiameter + CoreDiameter / 2)
            self.lum.set('y', 0)
            self.lum.set('y span', CoreDiameter + CoreDiameter / 2)
            self.lum.set("z", (ZSpan / 4) - 0.3e-6)
            self.lum.set('output Px', 1)
            self.lum.set('output Py', 1)
            self.lum.set('output Pz', 1)
            self.lum.set('output power', 1)



            # Power Monitor Output Port
            self.lum.addpower()
            self.lum.set('name', "Power_" + Port_Names[1])
            self.lum.set('monitor type', '2D X-normal')
            self.lum.set("x", -GC_SectionLenght/2 - InputLenght - TaperLength)
            self.lum.set('y', 0)
            self.lum.set('y span', CoreDiameter + CoreDiameter / 2)
            self.lum.set("z", Hight/2)
            self.lum.set("z span", z_Port_Span)
            self.lum.set('output Px', 1)
            self.lum.set('output Py', 1)
            self.lum.set('output Pz', 1)
            self.lum.set('output power', 1)
            
            
            # Add Global Power and Freq Monitor Z-Normal
            self.lum.addpower()
            self.lum.set('name', "Global_Power_Monitor Z-normal")
            self.lum.set('monitor type', '2D Z-normal')
            self.lum.set("x min", -GC_SectionLenght/2 - InputLenght - TaperLength)
            self.lum.set("x max",  fiber_xpos + CoreDiameter )
            self.lum.set("y", 0)
            self.lum.set("y span", WidthGC)
            self.lum.set('z', Hight/2)
            self.lum.set('output Px', 1)
            self.lum.set('output Py', 1)
            self.lum.set('output Pz', 1)
            self.lum.set('output power', 1)
        
        
            # Add Global Movie Monitor Z-Normal
            self.lum.addmovie()
            self.lum.set('name', "Global_Movie_Monitor Z-normal")
            self.lum.set('monitor type', '2D Z-normal')
            self.lum.set("x min", -GC_SectionLenght/2 - InputLenght - TaperLength)
            self.lum.set("x max",  fiber_xpos + CoreDiameter )
            self.lum.set("y", 0)
            self.lum.set("y span", WidthGC)
            self.lum.set('z', 0)
            
            # Add Global Power and Freq Monitor Y-Axis
            self.lum.addpower()
            self.lum.set('name', "Global_Power_Monitor Y-normal")
            self.lum.set('monitor type', '2D Y-normal')
            self.lum.set("x min", -GC_SectionLenght/2 - InputLenght - TaperLength)
            self.lum.set("x max",  fiber_xpos + CoreDiameter )
            self.lum.set("y", 0)
            self.lum.set('z', Hight / 2)
            self.lum.set("z span", z_Port_Span)
            self.lum.set('output Px', 1)
            self.lum.set('output Py', 1)
            self.lum.set('output Pz', 1)
            self.lum.set('output power', 1)
        
        
            # Add Global Movie Monitor Y-Axis
            self.lum.addmovie()
            self.lum.set('name', "Global_Movie_Monitor Y-normal")
            self.lum.set('monitor type', '2D Y-normal')
            self.lum.set("x min", -GC_SectionLenght/2 - InputLenght - TaperLength)
            self.lum.set("x max",  fiber_xpos + CoreDiameter )
            self.lum.set("y", 0)
            self.lum.set('z', Hight / 2)
            self.lum.set("z span", 2e-6)


            # Select Source
            self.lum.select('FDTD::ports')
            self.lum.set('source port', 'Input_SMF_Port')

        else:

            # Adds a Finite-Difference Time-Domain  (FDTD) solver region to the MODE simulation environment.
            self.lum.addfdtd()
            self.lum.set("x min", -GC_SectionLenght/2 - InputLenght )
            self.lum.set("x max",  fiber_xpos + CoreDiameter )
            self.lum.set("y", 0)
            self.lum.set("y span", WidthGC)
            self.lum.set("z", 0)
            self.lum.set("z span", (ZSpan / 2))
            self.lum.set('simulation temperature', 273.15 + 20)
            self.lum.set('z min bc', 'PML')
            self.lum.set('z max bc', 'PML')
            self.lum.set('y min bc', 'Anti-Symmetric')
            self.lum.set('y max bc', 'PML')
            self.lum.set('mesh type', 'auto non-uniform')
            self.lum.set('min mesh step', x_res)
            self.lum.set('set simulation bandwidth', 0)
            self.lum.set('global source center wavelength', WaveLength)
            self.lum.set('global source wavelength span', 0)
            
            
            
            # Detect Fiber position for Port exact aligment
            self.lum.select("SMF")
            fiber_xpos = self.lum.get("x")
            fiber_ypos = self.lum.get("y")
            fiber_zpos = self.lum.get("z")
            self.lum.select("SMF::core")
            fiber_core_diameter = 2 * self.lum.get("radius")
            fiber_core_index = self.lum.get("index")
            fiber_theta = Theta
            
            
            # Faser Port
            self.lum.addport()
            self.lum.set('name', Port_Names[0])
            self.lum.set('injection axis', "z-axis")
            self.lum.set('direction', "Backward")
            self.lum.set('mode selection', Mode)
            self.lum.set('theta', Theta)
            self.lum.set("x", fiber_xpos)
            self.lum.set('x span', CoreDiameter + CoreDiameter / 2)
            self.lum.set('y', 0)
            self.lum.set('y span', CoreDiameter + CoreDiameter / 2)
            self.lum.set("rotation offset", ZSpan / 4)
            self.lum.set("z", (ZSpan / 4) - 0.3e-6)


            # Output Port
            self.lum.addport()
            self.lum.set('name', Port_Names[1])
            self.lum.set("x", -GC_SectionLenght / 2 - InputLenght)
            self.lum.set('x span', CoreDiameter)
            self.lum.set('y', 0)
            self.lum.set('y span', WidthGC)
            self.lum.set("z", Hight / 2)
            self.lum.set("z span", z_Port_Span)
            self.lum.set('injection axis', "x-axis")
            self.lum.set('direction', "Forward")
            self.lum.set('mode selection', Mode)


            # Power Monitor SMF Port
            self.lum.addpower()
            self.lum.set('name', "Power_" + Port_Names[0])
            self.lum.set('monitor type', '2D Z-normal')
            self.lum.set("x", fiber_xpos)
            self.lum.set('x span', CoreDiameter + CoreDiameter / 2)
            self.lum.set('y', 0)
            self.lum.set('y span', CoreDiameter + CoreDiameter / 2)
            self.lum.set("z", (ZSpan / 4) - 0.3e-6)
            self.lum.set('output Px', 1)
            self.lum.set('output Py', 1)
            self.lum.set('output Pz', 1)
            self.lum.set('output power', 1)


            # Power Monitor Output Port
            self.lum.addpower()
            self.lum.set('name', "Power_" + Port_Names[1])
            self.lum.set('monitor type', '2D X-normal')
            self.lum.set("x", -GC_SectionLenght / 2 - InputLenght)
            self.lum.set('y', 0)
            self.lum.set('y span', WidthGC)
            self.lum.set("z", Hight / 2)
            self.lum.set("z span", z_Port_Span)
            self.lum.set('output Px', 1)
            self.lum.set('output Py', 1)
            self.lum.set('output Pz', 1)
            self.lum.set('output power', 1)
            
            
            
            # Add Global Power and Freq Monitor Z-Normal
            self.lum.addpower()
            self.lum.set('name', "Global_Power_Monitor Z-normal")
            self.lum.set('monitor type', '2D Z-normal')
            self.lum.set("x min", -GC_SectionLenght/2 - InputLenght)
            self.lum.set("x max",  fiber_xpos + CoreDiameter )
            self.lum.set("y", 0)
            self.lum.set("y span", WidthGC)
            self.lum.set('z', Hight/2)
            self.lum.set('output Px', 1)
            self.lum.set('output Py', 1)
            self.lum.set('output Pz', 1)
            self.lum.set('output power', 1)
        
        
            # Add Global Movie Monitor Z-Normal
            self.lum.addmovie()
            self.lum.set('name', "Global_Movie_Monitor Z-normal")
            self.lum.set('monitor type', '2D Z-normal')
            self.lum.set("x min", -GC_SectionLenght/2 - InputLenght)
            self.lum.set("x max",  fiber_xpos + CoreDiameter )
            self.lum.set("y", 0)
            self.lum.set("y span", WidthGC)
            self.lum.set('z', 0)
            
            
            # Add Global Power and Freq Monitor Y-Axis
            self.lum.addpower()
            self.lum.set('name', "Global_Power_Monitor Y-normal")
            self.lum.set('monitor type', '2D Y-normal')
            self.lum.set("x min", -GC_SectionLenght/2 - InputLenght)
            self.lum.set("x max",  fiber_xpos + CoreDiameter )
            self.lum.set("y", 0)
            self.lum.set('z', Hight / 2)
            self.lum.set("z span", z_Port_Span)
            self.lum.set('output Px', 1)
            self.lum.set('output Py', 1)
            self.lum.set('output Pz', 1)
            self.lum.set('output power', 1)
        
        
            # Add Global Movie Monitor Y-Axis
            self.lum.addmovie()
            self.lum.set('name', "Global_Movie_Monitor Y-normal")
            self.lum.set('monitor type', '2D Y-normal')
            self.lum.set("x min", -GC_SectionLenght/2 - InputLenght)
            self.lum.set("x max",  fiber_xpos + CoreDiameter )
            self.lum.set("y", 0)
            self.lum.set('z', Hight / 2)
            self.lum.set("z span", 2e-6)
            

 

            # Select Source
            self.lum.select('FDTD::ports')
            self.lum.set('source port', 'Input_SMF_Port')
            
         


         
    def setRingGratingCouplerFDTDSolver(self, Parameters):
        '''
       
        Parameters
        ----------
        Parameters : Dictionary
            Dictionary with all the data needet for the RingGratingCouplerFDTDSolver.
            Parameters['Substrate Height'] : int/float
                Substrate height.
            Parameters["Length GC"]: int/float
                Lenght of the Grating Coupler Area
            Parameters["Input Length GC"]: int/float
                An squere Waveguide with the same WG Height as the Grating coupler place before the Grating Coupler region will start. 
            Parameters["Output Length GC"]: int/float
                An squere Waveguide with the same WG Height as the Grating coupler place after the Grating Coupler region to finish the structure.
            Parameters["Width GC"]: int/float
                Widht of the Grating Coupler Area
            Parameters["Hight GC"]: int/float
                Hight of the Grating Coupler Material
            Parameters["GC Radius"]: int/float
                Radius of the Ring Grating Coupler in um. For Example "Parameters["GC Radius"] = 25e-6"
            Parameters['Taper Length'] : int/float
                  Length of the input Taper
            Parameters['Wavelength'] : int/float
                  Wavelength
            Parameters['x res'] : int/float
                  Mesh x-Axis
            Parameters["Mode"] : str
                Mode to choose from ("fundamental TE mode", "fundamental TM mode", "fundamental mode")
            Parameters["Port Span"] : list of floats/ints
                  List of x,y and z span of the Ports. For this simulation only y and z parametes will be taken.
        
       
            

        Returns
        -------
        None.

        '''
        
        SubstrateThickness = Parameters['Substrate Height']
        GC_SectionLenght = Parameters["Length GC"]
        InputLenght = Parameters["Input Length GC"]
        OutputLenght = Parameters["Output Length GC"]
        WidthGC = Parameters["Width GC"]
        Hight = Parameters["Hight GC"]
        GCRadius = Parameters["GC Radius"]
        TaperLength = Parameters['Taper Length']
        ZSpan = Parameters["SMF Z Span"]
        Theta = Parameters["SMF Theta"]
        CoreDiameter = Parameters["SMF Core Diameter"]
        x_res = Parameters['x res']
        WaveLength = Parameters['Wavelength']
        Mode = Parameters["Mode"]
        y_Port_Span = Parameters["Port Span"][1]
        z_Port_Span = Parameters["Port Span"][2]
        CladdingThickness = 0.7e-6
        
        
        
        
        

        # Device specifications
        Device_Length =  GC_SectionLenght + OutputLenght + TaperLength 
        Device_Width = WidthGC + ( GCRadius + GC_SectionLenght + OutputLenght  - GCRadius)
        FDTD_ZSpan = Hight + CladdingThickness + SubstrateThickness
       
        




        # Define Ports
        Port_Names = ["Input_SMF_Port", "Output"]

        self.lum.select("SMF")
        fiber_xpos = self.lum.get("x")

        # Adds a Finite-Difference Time-Domain  (FDTD) solver region to the MODE simulation environment.
        self.lum.addfdtd()
        self.lum.set("x min", -GC_SectionLenght/2 - 10e-6 /  2 - 1e-6)
        # self.lum.set("x max",  GCRadius + GC_SectionLenght/2 + OutputLenght +1e-6)
        self.lum.set("x max",  fiber_xpos + CoreDiameter )
        self.lum.set("y", 0)
        self.lum.set("y span", Device_Width/2)
        self.lum.set("z", 0)
        self.lum.set("z span",  SubstrateThickness + 2e-6)
        self.lum.set('simulation temperature', 273.15 + 20)
        self.lum.set('z min bc', 'PML')
        self.lum.set('z max bc', 'PML')
        self.lum.set('y min bc', 'Anti-Symmetric')
        self.lum.set('y max bc', 'PML')
        self.lum.set('mesh type', 'auto non-uniform')
        self.lum.set('min mesh step', x_res)
        self.lum.set('set simulation bandwidth', 0)
        self.lum.set('global source center wavelength', WaveLength)
        self.lum.set('global source wavelength span', 0)

        # Detect Fiber position for Port exact aligment
        self.lum.select("SMF")
        fiber_xpos = self.lum.get("x")
        fiber_ypos = self.lum.get("y")
        fiber_zpos = self.lum.get("z")
        self.lum.select("SMF::core")
        fiber_core_diameter = 2 * self.lum.get("radius")
        fiber_core_index = self.lum.get("index")
        fiber_theta = Theta

        # Faser Port
        self.lum.addport()
        self.lum.set('name', Port_Names[0])
        self.lum.set('injection axis', "z-axis")
        self.lum.set('direction', "Backward")
        self.lum.set('mode selection', Mode)
        self.lum.set('theta', Theta)
        self.lum.set("x", fiber_xpos)
        self.lum.set('x span', CoreDiameter + CoreDiameter / 2)
        self.lum.set('y', 0)
        self.lum.set('y span', CoreDiameter + CoreDiameter / 2)
        self.lum.set("z", (SubstrateThickness + 2e-6) / 2 )
        self.lum.set("rotation offset", ZSpan / 4)


        # Output Port
        self.lum.addport()
        self.lum.set('name', Port_Names[1])
        self.lum.set("x", -GC_SectionLenght/2 - 8e-6 /2)
        self.lum.set('x span', CoreDiameter)
        self.lum.set('y', 0)
        self.lum.set('y span', y_Port_Span)
        self.lum.set("z", Hight/2 )
        self.lum.set("z span", z_Port_Span )
        self.lum.set('injection axis', "x-axis")
        self.lum.set('direction', "Forward")
        self.lum.set('mode selection', Mode)




        # Power Monitor SMF Port
        self.lum.addpower()
        self.lum.set('name', "Power_"+ Port_Names[0])
        self.lum.set('monitor type', '2D Z-normal')
        self.lum.set("x", fiber_xpos)
        self.lum.set('x span', CoreDiameter + CoreDiameter / 2)
        self.lum.set('y', 0)
        self.lum.set('y span', CoreDiameter + CoreDiameter / 2)
        self.lum.set("z", (ZSpan / 4) - 0.3e-6)
        self.lum.set('output Px', 1)
        self.lum.set('output Py', 1)
        self.lum.set('output Pz', 1)
        self.lum.set('output power', 1)



        # Power Monitor Output Port
        self.lum.addpower()
        self.lum.set('name', "Power_" + Port_Names[1])
        self.lum.set('monitor type', '2D X-normal')
        self.lum.set("x", -GC_SectionLenght/2  - 8e-6 /2)
        self.lum.set('y', 0)
        self.lum.set('y span', y_Port_Span)
        self.lum.set("z", Hight/2)
        self.lum.set("z span", z_Port_Span)
        self.lum.set('output Px', 1)
        self.lum.set('output Py', 1)
        self.lum.set('output Pz', 1)
        self.lum.set('output power', 1)


        # Add Global Power and Freq Monitor
        self.lum.addpower()
        self.lum.set('name', "Global_Power_Monitor Z-normal")
        self.lum.set('monitor type', '2D Z-normal')
        self.lum.set("x min", -GC_SectionLenght/2  - 10e-6 /2 - 1e-6)
        # self.lum.set("x max",  GCRadius + GC_SectionLenght/2 + OutputLenght +1e-6)
        self.lum.set("x max",  fiber_xpos + CoreDiameter )
        self.lum.set("y", 0)
        self.lum.set("y span", Device_Width/2)
        self.lum.set('z', Hight/2)
        self.lum.set('output Px', 1)
        self.lum.set('output Py', 1)
        self.lum.set('output Pz', 1)
        self.lum.set('output power', 1)
        
        
        # Add Global Power and Freq Monitor
        self.lum.addmovie()
        self.lum.set('name', "Global_Movie_Monitor Z-normal")
        self.lum.set('monitor type', '2D Z-normal')
        self.lum.set("x min", -GC_SectionLenght/2  - 10e-6 /2 - 1e-6)
        # self.lum.set("x max",  GCRadius + GC_SectionLenght/2 + OutputLenght +1e-6)
        self.lum.set("x max",  fiber_xpos + CoreDiameter )
        self.lum.set("y", 0)
        self.lum.set("y span", Device_Width/2)
        self.lum.set('z', 0)
 


        # Add Global Power and Freq Monitor Y-Axis
        self.lum.addpower()
        self.lum.set('name', "Global_Power_Monitor Y-normal")
        self.lum.set('monitor type', '2D Y-normal')
        self.lum.set("x min", -GC_SectionLenght/2  - 10e-6 /2 - 1e-6)
        # self.lum.set("x max",  GCRadius + GC_SectionLenght/2 + OutputLenght + 1e-6)
        self.lum.set("x max",  fiber_xpos + CoreDiameter )
        self.lum.set("y", 0)
        self.lum.set('z', Hight / 2)
        self.lum.set("z span", z_Port_Span)
        self.lum.set('output Px', 1)
        self.lum.set('output Py', 1)
        self.lum.set('output Pz', 1)
        self.lum.set('output power', 1)
        
        
        # Add Global Movie Monitor Y-Axis
        self.lum.addmovie()
        self.lum.set('name', "Global_Movie_Monitor Y-normal")
        self.lum.set('monitor type', '2D Y-normal')
        self.lum.set("x min", -GC_SectionLenght/2  - 10e-6 /2 - 1e-6)
        # self.lum.set("x max",  GCRadius + GC_SectionLenght/2 + OutputLenght + 1e-6)
        self.lum.set("x max",  fiber_xpos + CoreDiameter )
        self.lum.set("y", 0)
        self.lum.set('z', Hight / 2)
        self.lum.set("z span", 2e-6)

        
        
        # Add Refractive index Monitor 
        self.lum.addindex()
        self.lum.set('name', "Refractive Index Monitor Y-normal")
        self.lum.set('monitor type', '2D Y-normal')
        self.lum.set("x min", -GC_SectionLenght/2  - 10e-6 /2 - 1e-6)
        # self.lum.set("x max",  GCRadius + GC_SectionLenght/2 + OutputLenght + 1e-6)
        self.lum.set("x max",  fiber_xpos + CoreDiameter )
        self.lum.set("y", 0)
        self.lum.set('z', Hight / 2)
        self.lum.set("z span", z_Port_Span)
        

        # Select Source
        self.lum.select('FDTD::ports')
        self.lum.set('source port', 'Input_SMF_Port')




    def setLenseFiberFDTD(self, Parameters):

        Lens_d = Parameters["Lense Diameter"]
        x_res = Parameters['x res']
        Wavelength = Parameters['Wavelength']

        # Set FDTD Solver directly
        self.lum.addfdtd()
        self.lum.set("x", 0)
        self.lum.set("x span", 15e-6)
        self.lum.set("y", 0)
        self.lum.set("y span", 15e-6)
        self.lum.set('simulation temperature', 273.15 + 20)
        self.lum.set("z min", -4e-6)
        self.lum.set("z max", 50e-6)
        self.lum.set('z min bc', 'PML')
        self.lum.set('z max bc', 'PML')
        self.lum.set('mesh type', 'auto non-uniform')
        self.lum.set('min mesh step', x_res)
        self.lum.set('set simulation bandwidth', 0)
        self.lum.set('global source center wavelength', Wavelength)
        self.lum.set('global source wavelength span', 0)

        # Add Gausssian Source
        self.lum.addgaussian()
        self.lum.set("injection axis", "z-axis")
        self.lum.set("direction", "Forward")
        self.lum.set("x", 0)
        self.lum.set("x span", Lens_d)
        self.lum.set("y", 0)
        self.lum.set("y span", Lens_d)
        self.lum.set("z", -2e-6)
        self.lum.set("waist radius w0", Lens_d / 2)

        # Add Index Monitor
        self.lum.addindex()
        self.lum.set("name", "index")
        self.lum.set("monitor type", "2D Y-normal")
        self.lum.set("x", 0)
        self.lum.set("x span", 10e-6)
        self.lum.set("y", 0)
        self.lum.set("z min", -3e-6)
        self.lum.set("z max", 50e-6)

        # Add Movie Monitor
        self.lum.addmovie()
        self.lum.set("name", "movie")
        self.lum.set("monitor type", "2D Y-normal")
        self.lum.set("x", 0)
        self.lum.set("x span", 10e-6)
        self.lum.set("y", 0)
        self.lum.set("z min", -3e-6)
        self.lum.set("z max", 50e-6)

        # Add Power and field 3D monitor
        self.lum.addpower()
        self.lum.set("name", "Power 3D monitor")
        self.lum.set("monitor type", "3D")
        self.lum.set("x", 0)
        self.lum.set("x span", 10e-6)
        self.lum.set("y", 0)
        self.lum.set("y span", 10e-6)
        self.lum.set("z min", -3e-6)
        self.lum.set("z max", 50e-6)
        self.lum.set("output Px", 1)
        self.lum.set("output Py", 1)
        self.lum.set("output Pz", 1)

        # Reflection 2D Monitor
        self.lum.addpower()
        self.lum.set("name", "2D Monitor")
        self.lum.set("monitor type", "2D X-normal")
        self.lum.set("x", 0)
        self.lum.set("y", 0)
        self.lum.set("y span", 10e-6)
        self.lum.set("z", -2e-6)
        self.lum.set("z max", 50e-6)
        self.lum.set("output Px", 1)
        self.lum.set("output Py", 1)
        self.lum.set("output Pz", 1)

        self.lum.addpower()
        self.lum.set("name", "R")
        self.lum.set("monitor type", "2D Z-normal")
        self.lum.set("x", 0)
        self.lum.set("x span", 10e-6)
        self.lum.set("y", 0)
        self.lum.set("y span", 10e-6)
        self.lum.set("z", -3e-6)
        self.lum.set("output Px", 1)
        self.lum.set("output Py", 1)
        self.lum.set("output Pz", 1)





    def setStraightWaveguideEMESolver(self, Parameters):
        '''


        Parameters
        ----------
        Parameters : Dictionary
            Dictionary with all the data needet for the StraightWaveguideEMESolver.
            Parameters['Substrate Height'] : int/float
               Substrate height.
            Parameters['WG Length'] : int/float
               Waveguide Length
            Parameters['WG Height'] : int/float
               Waveguide hight. Also the height of the MMI section
            Parameters['WG Width'] : int/float
               Waveguide width.
            Parameters["Taper"] : boolen
               If Taper == False, only straight Waveguide will be simulated,
               If Taper == True an Taper will be simulated
            Parameters['Taper Width'] : int/float
               Taper backside Width. Taper Fronside width is the width of the Waveguide
            Parameters['Taper Length'] : int/float
               Taper Length
            Parameters['y res']: int/float
                 EME Mesh resolutio,
            Parameters['z res']: int/float
                 EME Mesh resolutio,
            Parameters['Slab Height'] : int/float
               Height of the slab.
            Parameters['Wavelength'] : int/float
               Wavelength
            Parameters["Waveguide Angle"] : int/float
               This Parameter will set the theta ratation angle of the port. It can be 90 or 180.
            Parameters["Port Span"] : list of floats/ints
               List of x,y and z span of the Ports. For this simulation only y and z parametes will be taken.
            Parameters["Mode"] : str
                Mode to choose from ("fundamental TE mode", "fundamental TM mode", "fundamental mode")
            Parameters["Taper Type"] : anything, optional
                    This function will check if you have set Parameters["Taper Type"] to anaything, for example "Parameters["Taper Type"]=1" 
                    and if so it will design an Inverse Taper Structure with no Cladding. Here the option "Cladding" is not active and will be ignored.
                    If the user didnt give the "Taper Type" as dictionary key, then an normal taper structure will be simulated.
                    
                    If Parameters["Taper Type"] is given, themn the user need to set couple more parameters:
                        Parameters['PWB Taper Width Back'] : int/float
                            Photonic Wirebonding (PWB) Width back side (to the Photonic Wire Bonding)
                        Parameters['PWB Taper Hight Back'] : int/float
                            Photonic Wire Bonding Height back side (to the Photonic Wire Bonding)
                        Parameters['PWB Taper Width Front'] : int/float
                            Photonic Wirebonding (PWB) Width front side (to the photonic waveguide)
                        Parameters['PWB Taper Hight Front'] : int/float
                            Photonic Wire Bonding Height front side (to the photonic waveguide)
                        Parameters['PWB Taper Length'] : int/float
                            Length of the Photonic Wire Bonding Taper
                            
        Returns
        -------
        None.

        '''
    




        Substrate_Height = Parameters['Substrate Height']
        WG_Length = Parameters['WG Length']
        WG_Height = Parameters['WG Height']
        WG_Width = Parameters['WG Width']
        angle = Parameters["Waveguide Angle"]
        y_res = Parameters['y res']
        z_res = Parameters['z res']
        Slab_Height = Parameters['Slab Height']
        WaveLength = Parameters['Wavelength']
        Mode = Parameters["Mode"]
        y_Port_Span = Parameters["Port Span"][1]
        z_Port_Span = Parameters["Port Span"][2]
        Taper = Parameters['Taper']
        TaperWidth = Parameters['Taper Width']
        TaperLength = Parameters['Taper Length']
        


        # Device specifications
        Device_Length = WG_Length
        Device_Width = 2*WG_Length + WaveLength * 2  # MMI_Width

        max_slabH = Slab_Height
        # MonitorHeight = Substrate_Height + (max_slabH + WG_Height) / 2
        MonitorHeight = Substrate_Height + max_slabH + WG_Height/2
        # Ports_mid = (max_slabH + WG_Height) / 2
        Ports_mid = max_slabH + WG_Height/2
        EME_WGLength = WG_Length * np.cos(angle * np.pi / 180)


        if Taper == False:
            if angle == 0:

                # Adds a Eigenmode Expansion (EME) solver region to the MODE simulation environment.
                self.lum.addeme()
                self.lum.set('simulation temperature', 273.15 + 20)
                self.lum.set("x min", 0)
                self.lum.set("y", 0)
                self.lum.set("y span", Device_Width)
                self.lum.set("z", Substrate_Height)
                self.lum.set("z span", 4e-6)
                self.lum.set("wavelength", WaveLength)
                self.lum.set("z min bc", "PML")
                self.lum.set("z max bc", "PML")
                self.lum.set("y min bc", "PML")
                self.lum.set("y max bc", "PML")

                # set cell properties
                self.lum.set("number of cell groups", 1)
                self.lum.set("group spans", np.array([[EME_WGLength]]))
                self.lum.set("cells", np.array([[30]]))
                self.lum.set("subcell method", np.array([[1]]))

                # Modes to Calculate
                self.lum.set('number of modes for all cell groups', 20)

                # Mesh Cells
                self.lum.set("define y mesh by", "maximum mesh step")
                self.lum.set("dy", y_res)
                self.lum.set("define z mesh by", "maximum mesh step")
                self.lum.set("dz", z_res)
                self.lum.set('fit materials with multi-coefficient model', 1)
                self.lum.set('wavelength start', 0.4e-6)
                self.lum.set('wavelength stop', 2e-6)

                # Define Ports

                yPos = [0, 0]
                yPos_span = [y_Port_Span, y_Port_Span]
                portLoc = ["left", "right"]
                theta = [0, 0]

                for i in range(2):
                    self.lum.select("EME::Ports::port_" + str(i + 1))
                    self.lum.set("port location", portLoc[i])
                    self.lum.set("use full simulation span", 0)
                    self.lum.set("y", yPos[i])
                    self.lum.set("y span", yPos_span[i])
                    self.lum.set("z", Ports_mid)
                    self.lum.set("z span", z_Port_Span)
                    self.lum.set("mode selection", Mode)
                    self.lum.set("theta", theta[i])

                # Add monitor
                self.lum.addemeprofile()
                self.lum.set("x", (Device_Length / 2))
                self.lum.set("x span", Device_Length)
                self.lum.set("y", 0)
                self.lum.set("y span", Device_Width)
                self.lum.set("z", MonitorHeight)


            else:

                # Calc Output Loc
                sideLength = np.cos(angle*np.pi/180)* WG_Length
                sideHight = np.sqrt(WG_Length**2 - sideLength**2)

                # Adds a Eigenmode Expansion (EME) solver region to the MODE simulation environment.
                self.lum.addeme()
                self.lum.set("x min", 0)
                self.lum.set("y", WG_Length / 4)
                self.lum.set('simulation temperature', 273.15 + 20)
                self.lum.set("y span", Device_Width)
                self.lum.set("z", Substrate_Height)
                self.lum.set("z span", 4e-6)
                self.lum.set("wavelength", WaveLength)
                self.lum.set("z min bc", "PML")
                self.lum.set("z max bc", "PML")
                self.lum.set("y min bc", "PML")
                self.lum.set("y max bc", "PML")

                # set cell properties
                self.lum.set("number of cell groups", 1)
                self.lum.set("group spans", np.array([[EME_WGLength]]))
                self.lum.set("cells", np.array([[30]]))
                self.lum.set("subcell method", np.array([[1]]))

                # Modes to Calculate
                self.lum.set('number of modes for all cell groups', 20)

                # Mesh Cells
                self.lum.set("define y mesh by", "maximum mesh step")
                self.lum.set("dy", y_res)
                self.lum.set("define z mesh by", "maximum mesh step")
                self.lum.set("dz", z_res)
                self.lum.set('fit materials with multi-coefficient model', 1)
                self.lum.set('wavelength start', 0.4e-6)
                self.lum.set('wavelength stop', 2e-6)



                # Define Ports

                yPos = [-WG_Length / 4, -WG_Length / 4 + sideHight]
                yPos_span = [y_Port_Span, y_Port_Span]
                portLoc = ["left", "right"]
                theta = [angle, angle]


                for i in range(2):
                    self.lum.select("EME::Ports::port_" + str(i + 1))
                    self.lum.set("port location", portLoc[i])
                    self.lum.set("use full simulation span", 0)
                    self.lum.set("y", yPos[i])
                    self.lum.set("y span", yPos_span[i])
                    self.lum.set("z", Ports_mid)
                    self.lum.set("z span", z_Port_Span)
                    self.lum.set("mode selection", Mode)
                    self.lum.set("theta", theta[i])

                # Add monitor
                self.lum.addemeprofile()
                self.lum.set("x", (Device_Length / 2))
                self.lum.set("x span", Device_Length)
                self.lum.set("y", (Device_Length / 4))
                self.lum.set("y span", Device_Width)
                self.lum.set("z", MonitorHeight)



        elif Taper == True:
            # Device specifications
            Device_Length = TaperLength
            Device_Width = 2 * TaperLength + WaveLength * 2  # MMI_Width

            max_slabH = Slab_Height
            MonitorHeight = Substrate_Height + max_slabH + WG_Height/2
            Ports_mid = max_slabH + WG_Height/2
            EME_WGLength = TaperLength * np.cos(angle * np.pi / 180)
            
            if Slab_Height == 0:
                # # Device specifications
                MonitorHeight = Substrate_Height + WG_Height/2
                Ports_mid = Substrate_Height + WG_Height/2
              
            
            else:
                # # Device specifications
                max_slabH = Slab_Height
                MonitorHeight = Substrate_Height + max_slabH + WG_Height/2
                Ports_mid = max_slabH + Substrate_Height  + WG_Height/2
                
            if "Taper Type" in list(Parameters.keys()):
                TaperType = "Inverse"
                TaperWidthF = Parameters['PWB Taper Width Front']
                TaperWidthB = Parameters['PWB Taper Width Back']
                TaperHightB = Parameters['PWB Taper Hight Back']
                TaperHightF = Parameters['PWB Taper Hight Front']
                TaperLength_PWB = Parameters['PWB Taper Length']
            else:
                TaperType = "Normal"
                
                

            
            
            if TaperType == "Normal":
                # Adds a Eigenmode Expansion (EME) solver region to the MODE simulation environment.
                self.lum.addeme()
                self.lum.set('simulation temperature', 273.15 + 20)
                self.lum.set("x min", -TaperLength/2 + 0.1e-6)
                self.lum.set("y", 0)
                self.lum.set("y span", Device_Width)
                self.lum.set("z", Substrate_Height)
                self.lum.set("z span", 4e-6)
                self.lum.set("wavelength", WaveLength)
                self.lum.set("z min bc", "PML")
                self.lum.set("z max bc", "PML")
                self.lum.set("y min bc", "PML")
                self.lum.set("y max bc", "PML")

                # set cell properties
                self.lum.set("number of cell groups", 1)
                self.lum.set("group spans", np.array([[EME_WGLength - 0.2e-6]]))
                self.lum.set("cells", np.array([[30]]))
                self.lum.set("subcell method", np.array([[1]]))

                # Modes to Calculate
                self.lum.set('number of modes for all cell groups', 20)

                # Mesh Cells
                self.lum.set("define y mesh by", "maximum mesh step")
                self.lum.set("dy", y_res)
                self.lum.set("define z mesh by", "maximum mesh step")
                self.lum.set("dz", z_res)
                self.lum.set('fit materials with multi-coefficient model', 1)
                self.lum.set('wavelength start', 0.4e-6)
                self.lum.set('wavelength stop', 2e-6)

                # Define Ports
                Diff_Span = y_Port_Span - WG_Width
                yPos = [0, 0]
                yPos_span = [y_Port_Span, TaperWidth + Diff_Span]
                portLoc = ["left", "right"]
                theta = [0, 0]

                for i in range(2):
                    self.lum.select("EME::Ports::port_" + str(i + 1))
                    self.lum.set("port location", portLoc[i])
                    self.lum.set("use full simulation span", 0)
                    self.lum.set("y", yPos[i])
                    self.lum.set("y span", yPos_span[i])
                    self.lum.set("z", Ports_mid)
                    self.lum.set("z span", z_Port_Span)
                    self.lum.set("mode selection", Mode)
                    self.lum.set("theta", theta[i])

                # Add monitor
                self.lum.addemeprofile()
                self.lum.set("x", 0)
                self.lum.set("x span", Device_Length)
                self.lum.set("y", 0)
                self.lum.set("y span", Device_Width)
                self.lum.set("z", MonitorHeight)
           
                
           
            else:
                # Adds a Eigenmode Expansion (EME) solver region to the MODE simulation environment.
                self.lum.addeme()
                self.lum.set('simulation temperature', 273.15 + 20)
                self.lum.set("x min", -TaperLength/2 + 0.1e-6)
                self.lum.set("y", 0)
                self.lum.set("y span", Device_Width)
                # self.lum.set("z", Substrate_Height)
                # self.lum.set("z span", 4e-6)
                self.lum.set("z min", -Substrate_Height ) 
                self.lum.set("z max", TaperHightB*2)
                self.lum.set("wavelength", WaveLength)
                self.lum.set("z min bc", "PML")
                self.lum.set("z max bc", "PML")
                self.lum.set("y min bc", "PML")
                self.lum.set("y max bc", "PML")

                # set cell properties
                self.lum.set("number of cell groups", 1)
                self.lum.set("group spans", np.array([[EME_WGLength - 0.2e-6]]))
                self.lum.set("cells", np.array([[30]]))
                self.lum.set("subcell method", np.array([[1]]))

                # Modes to Calculate
                self.lum.set('number of modes for all cell groups', 20)

                # Mesh Cells
                self.lum.set("define y mesh by", "maximum mesh step")
                self.lum.set("dy", y_res)
                self.lum.set("define z mesh by", "maximum mesh step")
                self.lum.set("dz", z_res)
                self.lum.set('fit materials with multi-coefficient model', 1)
                self.lum.set('wavelength start', 0.4e-6)
                self.lum.set('wavelength stop', 2e-6)
                
                
                # Define Ports
                Diff_Span = y_Port_Span - WG_Width
                yPos = [0, 0]
                yPos_span = [y_Port_Span, TaperWidth + Diff_Span]
                theta = [0, 0]
                portLoc = ["left", "right"]
                name = ['Input', 'Output']
                
                self.lum.select("EME")
                Span_z = self.lum.get("z span")
                yPos_span = [ TaperWidthB + Diff_Span , y_Port_Span ]
                z_Pos = [ -Span_z/2 + Substrate_Height + max_slabH + TaperHightB/2,  -Span_z/2 + Substrate_Height + max_slabH + TaperHightF/2 ]
                z_Span = [ TaperHightB + z_Port_Span , TaperHightF/2  + z_Port_Span]
                for i in range(2):
                    self.lum.select("EME::Ports::port_" + str(i + 1))
                    self.lum.set("port location", portLoc[i])
                    self.lum.set("use full simulation span", 0)
                    self.lum.set("y", yPos[i])
                    self.lum.set("y span", yPos_span[i])
                    self.lum.set("z", z_Pos[i])
                    self.lum.set("z span", z_Span[i])
                    self.lum.set("mode selection", Mode)
                    self.lum.set("theta", theta[i])
                


                # Add monitor
                self.lum.addemeprofile()
                self.lum.set("x", 0)
                self.lum.set("x span", Device_Length)
                self.lum.set("y", 0)
                self.lum.set("y span", Device_Width)
                self.lum.set("z", MonitorHeight)
            
            


    def setDCEMESolver(self, Parameters):
        '''

        Parameters
        ----------
        Parameters : Dictionary
            Dictionary with all the data needet for the DCEMESolver.
            Parameters['Substrate Height'] : float/int
                Height of the Substrate
            Parameters['Substrate Width'] : float/int
                Width of the MMI
            Parameters['DC Length'] : float/int
                Length of the Directional coupler
            Parameters['WG Height'] : float/int
                Height of the Waveguide
            Parameters['WG Width'] : float/int
                Waveguide Width
            Parameters['Position Offset'] : float/int
                Positional offser of the waveguides. If posOffset the two Waveguides
                will be offset of the middle position (y = 0) by the half of there
                Width. In this case they will not overlap if the Offset is 0.
            Parameters['y res'] : float/int
                Mesh resolution for the y-Axis
            Parameters['z res'] : float/int
                Mesh resolution for the z Axis
            Parameters['Slab Height'] : float/int
                Slab height.
            Parameters['Wavelength'] : float/int
                Wavelength
            Parameters["Port Span"]: list of int/floats
                Parameters["Port Span"] = [Span of Port in x direction, Span of Port in y direction, Span of Port in z direction]
            Parameters["Mode"] : str
                Mode to choose from ("fundamental TE mode", "fundamental TM mode", "fundamental mode")

        Returns
        -------
        None.

        '''
        


        Substrate_Height = Parameters['Substrate Height']
        DC_Lenght = Parameters['DC Length']
        WG_Height = Parameters['WG Height']
        WG_Width = Parameters['WG Width']
        posOffset = Parameters['Position Offset']
        y_res = Parameters['y res']
        z_res = Parameters['z res']
        Slab_Height = Parameters['Slab Height']
        WaveLength = Parameters['Wavelength']
        Mode = Parameters["Mode"]
        y_Port_Span = Parameters["Port Span"][1]
        z_Port_Span = Parameters["Port Span"][2]


        # Device specifications
        Device_Length = DC_Lenght
        Device_Width = WG_Width * 10 + WaveLength * 2  # MMI_Width

        max_slabH = Slab_Height
        MonitorHeight = Substrate_Height + max_slabH + WG_Height/2
        Ports_mid = max_slabH + WG_Height/2

        # Adds a Eigenmode Expansion (EME) solver region to the MODE simulation environment.
        self.lum.addeme()
        self.lum.set("x min", -(Device_Length / 2))
        self.lum.set("y", 0)
        self.lum.set('simulation temperature', 273.15 + 20)
        self.lum.set("y span", Device_Width)
        self.lum.set("z", Substrate_Height)
        self.lum.set("z span", 4e-6)
        self.lum.set("wavelength", WaveLength)
        self.lum.set("z min bc", "PML")
        self.lum.set("z max bc", "PML")
        self.lum.set("y min bc", "PML")
        self.lum.set("y max bc", "PML")

        # set cell properties
        self.lum.set("number of cell groups", 1)
        self.lum.set("group spans", np.array([[DC_Lenght]]))
        self.lum.set("cells", np.array([[10]]))
        self.lum.set("subcell method", np.array([[0]]))

        # Modes to Calculate
        self.lum.set('number of modes for all cell groups', 20)

        # Mesh Cells
        self.lum.set("define y mesh by", "maximum mesh step")
        self.lum.set("dy", y_res)
        self.lum.set("define z mesh by", "maximum mesh step")
        self.lum.set("dz", z_res)
        self.lum.set('fit materials with multi-coefficient model', 1)
        self.lum.set('wavelength start', 0.4e-6)
        self.lum.set('wavelength stop', 2e-6)

        # Define Ports
        yPort_vec = [posOffset / 2 + WG_Width / 2, -(posOffset / 2 + WG_Width / 2), posOffset / 2 + WG_Width / 2, -(posOffset / 2 + WG_Width / 2)]
        portLoc = ["left", "left", "right", "right"]

        overLapp = yPort_vec[0] - y_Port_Span/2
        if overLapp < 0:
            raise ValueError("!!! CAUTION !!! - The Ports are overlapping at the middle! Please change the Y Port Span or move the Waveguides away from each other!")
        else:
            pass

        for i in range(2):
            self.lum.addemeport()

        for i in range(4):
            self.lum.select("EME::Ports::port_" + str(i + 1))
            self.lum.set("port location", portLoc[i])
            self.lum.set("use full simulation span", 0)
            self.lum.set('y', yPort_vec[i])
            self.lum.set('y span', y_Port_Span)
            self.lum.set("z", Ports_mid)
            self.lum.set("z span", z_Port_Span)
            self.lum.set("mode selection", Mode)

        # Add monitor
        self.lum.addemeprofile()
        self.lum.set("x min", -(Device_Length / 2))
        self.lum.set("x max", (Device_Length / 2))
        self.lum.set("y min", -Device_Width / 2)
        self.lum.set("y max", Device_Width / 2)
        self.lum.set("z", MonitorHeight)




    def setMMI2x1EMESolver(self, Parameters):
    
        '''


        Parameters
        ----------
        Parameters : Dictionary
            Dictionary with all the data needet for the MMI2x1EMESolver.
            Parameters['Substrate Height'] : int/float
                Height of the slab.
            Parameters['MMI Width'] : int/float
                Width of MMI
            Parameters['MMI Length'] : int/float
                Length of MMI
            Parameters['angle'] : int/float
                Angle of the Waveguide Walls. it is calculated WG_angle = 90 - angle.
                For anfle = 90 we get a perfect rect!
            Parameters['WG Height'] : int/float
                Heigth of waveguide
            Parameters['WG Width'] : int/float
                Width og waveguide
            Parameters['WG Length'] : int/float
                Length of waveguide
            Parameters['Position Offset'] : int/float
                Offset between the waveguides. If Taper == True then this become the offset
                betweent he tapers wider sides. Waveguide and Tapers cannot be placed ourside
                the MMI structure. The minimum distance between Taper and Waveguide is 1 um
                becouse of manufactering restrictions in the University.
            Parameters['Offset Input'] : int/float
                Input waveguide/taper offset.
            Parameters["Taper"]: boolen 
                Add Taper to the structure on the input and output waveguids
            Parameters['Taper Length']: int/float
                Lenght of the Taper in Parameters["Taper"] = True
            Parameters['y res'] : int/float
                Mesh cell sizes.
            Parameters['z res'] : int/float
                Mesh cell sizes.
            Parameters['Slab Height'] : int/float
                Height of the slab.
            Parameters['Wavelength'] : int/float
                Wavelength.
            Parameters["Mode"] : str
                Mode to choose from ("fundamental TE mode", "fundamental TM mode", "fundamental mode")
            Parameters["Port Span"]: list of int/floats
                Parameters["Port Span"] = [Span of Port in x direction, Span of Port in y direction, Span of Port in z direction]
            Parameters["Offset Output"] : anything, optional
                This function will allow the user to move the outputs in oposite direction. Please dont use it since is there only 
                becouse the maschine of our physic departmant had some proiblems with the LNOI objects design. 

        Returns
        -------
        None.

        '''


        Substrate_Height = Parameters['Substrate Height']
        MMI_Width = Parameters['MMI Width']
        MMI_Length = Parameters['MMI Length']
        angle = Parameters['angle']
        WG_Height = Parameters['WG Height']
        WG_Width = Parameters['WG Width']
        WG_Length = Parameters['WG Length']
        OffsetInput = Parameters['Offset Input']
        posOffset = Parameters['Position Offset']
        y_res = Parameters['y res']
        z_res = Parameters['z res']
        Slab_Height = Parameters['Slab Height']
        WaveLength = Parameters['Wavelength']
        TaperLength = Parameters['Taper Length']
        Taper = Parameters['Taper']
        Mode = Parameters["Mode"]
        y_Port_Span = Parameters["Port Span"][1]
        z_Port_Span = Parameters["Port Span"][2]



        if Taper == False:

            # Device specifications
            Device_Length = MMI_Length + 2 * WG_Length
            Device_Width = MMI_Width + WaveLength * 2  # MMI_Width

            # 1x2 MMI
            max_slabH = Slab_Height
            # MonitorHeight = Substrate_Height + (max_slabH + WG_Height) / 2
            MonitorHeight = Substrate_Height + max_slabH + WG_Height/2
            # Ports_mid = (max_slabH + WG_Height) / 2
            Ports_mid = max_slabH + WG_Height/2

            # Adds a Eigenmode Expansion (EME) solver region to the MODE simulation environment.
            self.lum.addeme()
            self.lum.set("x min", -(Device_Length / 2)+0.1e-6)
            self.lum.set("y", 0)
            self.lum.set('simulation temperature', 273.15 + 20)
            self.lum.set("y span", Device_Width)
            self.lum.set("z", Substrate_Height)
            self.lum.set("z span", 4e-6)
            self.lum.set("wavelength", WaveLength)
            self.lum.set("z min bc", "PML")
            self.lum.set("z max bc", "PML")
            self.lum.set("y min bc", "PML")
            self.lum.set("y max bc", "PML")

            # set cell properties no Taper
            self.lum.set("number of cell groups", 3)
            self.lum.set("group spans", np.array([[WG_Length-0.1e-6], [MMI_Length], [WG_Length-0.1e-6]]))
            self.lum.set("cells", np.array([[20], [3], [20]]))
            self.lum.set("subcell method", np.array([[0], [0], [0]]))

        elif Taper == True:

            # Device specifications
            Device_Length = MMI_Length + 2 * WG_Length + 2 * TaperLength
            Device_Width = MMI_Width + WaveLength * 2  # MMI_Width

            # 1x2 MMI
            max_slabH = Slab_Height
            # MonitorHeight = Substrate_Height + (max_slabH + WG_Height) / 2
            MonitorHeight = Substrate_Height + max_slabH + WG_Height/2
            # Ports_mid = (max_slabH + WG_Height) / 2
            Ports_mid = max_slabH + WG_Height/2

            # Adds a Eigenmode Expansion (EME) solver region to the MODE simulation environment.
            self.lum.addeme()
            self.lum.set("x min", -(Device_Length / 2) + 0.1e-6)
            self.lum.set("y", 0)
            self.lum.set('simulation temperature', 273.15 + 20)
            self.lum.set("y span", Device_Width)
            self.lum.set("z", Substrate_Height)
            self.lum.set("z span", 4e-6)
            self.lum.set("wavelength", WaveLength)
            self.lum.set("z min bc", "PML")
            self.lum.set("z max bc", "PML")
            self.lum.set("y min bc", "PML")
            self.lum.set("y max bc", "PML")

            Device_Length = MMI_Length + 2 * WG_Length + 2 * TaperLength
            # set cell properties
            self.lum.set("number of cell groups", 5)
            self.lum.set("group spans",
                         np.array([[WG_Length-0.1e-6], [TaperLength], [MMI_Length], [TaperLength], [WG_Length-0.1e-6]]))
            self.lum.set("cells", np.array([[3], [20, ], [3], [20], [3]]))
            self.lum.set("subcell method", np.array([[0], [1], [0], [1], [0]]))
        else:
            raise ValueError("Incorect Taper variable!. Possible Taper values are Taper = False or Taper = True.")

        # Modes to Calculate
        self.lum.set('number of modes for all cell groups', 20)

        # Mesh Cells
        self.lum.set("define y mesh by", "maximum mesh step")
        self.lum.set("dy", y_res)
        self.lum.set("define z mesh by", "maximum mesh step")
        self.lum.set("dz", z_res)
        self.lum.set('fit materials with multi-coefficient model', 1)
        self.lum.set('wavelength start', 0.4e-6)
        self.lum.set('wavelength stop', 2e-6)

        # creating the MMI
        max_MMIH = WG_Height

        # Triangle EQ for MMI Width
        x = abs(max_MMIH / (np.cos((angle) * np.pi / 180)))  # in Radians
        extention = np.sqrt(x ** 2 - max_MMIH ** 2)
        MMI_Wid = MMI_Width + 2 * extention

        # Positions of the Input and Output WGs
        # Triangle EQ for MMI Width
        x = abs(max_MMIH / (np.cos((angle) * np.pi / 180)))  # in Radians
        extention = np.sqrt(x ** 2 - max_MMIH ** 2)
        WG_W = WG_Width + 2 * extention
        OffMin = -MMI_Wid / 2
        OffMax = MMI_Wid / 2
        
        
        offset_Set_R = posOffset / 2 + WG_W / 2 + WG_Width / 2
        if OffsetInput - WG_W <= OffMin or OffsetInput + WG_W >= OffMax or offset_Set_R > OffMax:
            self.lum.deleteall()
            raise ValueError(
                'You are Trying to move the input or output Waveguide outside the MMI area. This is not possible! Max Input/Output Offset = ' + str(
                    OffMax) + ' and Min Input/output Offset = ' + str(OffMin))


        # Define Ports
        else:
            # if OffsetInput > 0:
            #     max_yPos = [(OffsetInput) * 2, 0, (WG_Width / 2 + posOffset / 2) * 2]
            #     min_yPos = [0, -(WG_Width / 2 + posOffset / 2) * 2, 0]
            #     portLoc = ["right", "left", "left"]
            # elif OffsetInput < 0:
            #     max_yPos = [0, 0, (WG_Width / 2 + posOffset / 2) * 2]
            #     min_yPos = [-(abs(OffsetInput)) * 2, -(WG_Width / 2 + posOffset / 2) * 2, 0]
            #     portLoc = ["right", "left", "left"]

            # else:
            #     max_yPos = [(WG_W / 2 + OffsetInput) * 2, 0, (WG_Width / 2 + posOffset / 2) * 2]
            #     min_yPos = [-(WG_W / 2 + OffsetInput) * 2, -(WG_Width / 2 + posOffset / 2) * 2, 0]
            #     portLoc = ["right", "left", "left"]

            yPort_vec = [OffsetInput, -(posOffset / 2 + WG_Width / 2), posOffset / 2 + WG_Width / 2]
            portLoc = ["right", "left", "left"]

            overLapp = yPort_vec[2] - y_Port_Span/2
            if overLapp < 0:
                raise ValueError("!!! CAUTION !!! - The Ports are overlapping at the middle! Please change the Y Port Span or move the Waveguides away from each other!")
            else:
                pass

            for i in range(1):
                self.lum.addemeport()

            for i in range(3):
                self.lum.select("EME::Ports::port_" + str(i + 1))
                self.lum.set("port location", portLoc[i])
                self.lum.set("use full simulation span", 0)
                self.lum.set('y', yPort_vec[i])
                self.lum.set('y span', y_Port_Span)
                # self.lum.set("y min", min_yPos[i])
                # self.lum.set("y max", max_yPos[i])
                self.lum.set("z", Ports_mid)
                self.lum.set("z span", z_Port_Span)
                self.lum.set("mode selection", Mode)

        # Define the Motinot
        self.lum.addemeprofile()
        self.lum.set("x min", -(Device_Length / 2))
        self.lum.set("x max", (Device_Length / 2))
        self.lum.set("y min", -Device_Width / 2)
        self.lum.set("y max", Device_Width / 2)
        self.lum.set("z", MonitorHeight)






    def setMMI2x2EMESolver(self, Parameters):
        '''


        Parameters
        ----------
        Parameters : Dictionary
            Dictionary with all the data needet for the MMI2x2EMESolver.
            Parameters['Substrate Height'] : int/float
                Height of the slab.
            Parameters['MMI Width'] : int/float
                Width of MMI
            Parameters['MMI Length'] : int/float
                Length of MMI
            Parameters['WG Height'] : int/float
                Heigth of waveguide
            Parameters['WG Width'] : int/float
                Width og waveguide
            Parameters['WG Length'] : int/float
                Length of waveguide
            Parameters['Position Offset'] : int/float
                Offset between the waveguides. If Taper == True then this become the offset
                betweent he tapers wider sides. Waveguide and Tapers cannot be placed ourside
                the MMI structure. The minimum distance between Taper and Waveguide is 1 um
                becouse of manufactering restrictions in the University.
            Parameters["Taper"]: boolen 
                Add Taper to the structure on the input and output waveguids
            Parameters['Taper Length']: int/float
                Lenght of the Taper in Parameters["Taper"] = True
            Parameters['Slab Height'] : int/float
                Height of the slab.
            Parameters['y res'] : int/float
                Mesh cell sizes.
            Parameters['z res'] : int/float
                Mesh cell size.
            Parameters['Wavelength'] : int/float
                Wavelength.
            Parameters["Mode"] : str
                Mode to choose from ("fundamental TE mode", "fundamental TM mode", "fundamental mode")
            Parameters["Port Span"]: list of int/floats
                Parameters["Port Span"] = [Span of Port in x direction, Span of Port in y direction, Span of Port in z direction]

        Returns
        -------
        None.

        '''
    
 

        Substrate_Height = Parameters['Substrate Height']
        MMI_Width = Parameters['MMI Width']
        MMI_Length = Parameters['MMI Length']
        WG_Height = Parameters['WG Height']
        WG_Width = Parameters['WG Width']
        WG_Length = Parameters['WG Length']
        posOffset = Parameters['Position Offset']
        y_res = Parameters['y res']
        z_res = Parameters['z res']
        Slab_Height = Parameters['Slab Height']
        WaveLength = Parameters['Wavelength']
        TaperLength = Parameters['Taper Length']
        Taper = Parameters['Taper']
        Mode = Parameters["Mode"]
        y_Port_Span = Parameters["Port Span"][1]
        z_Port_Span = Parameters["Port Span"][2]




        if Taper == False:

            # Device specifications
            Device_Length = MMI_Length + 2 * WG_Length
            Device_Width = MMI_Width + WaveLength * 2  # MMI_Width
            max_slabH = Slab_Height
            # MonitorHeight = Substrate_Height + (max_slabH + WG_Height) / 2
            MonitorHeight = Substrate_Height + max_slabH + WG_Height/2
            # Ports_mid = (max_slabH + WG_Height) / 2
            Ports_mid = max_slabH + WG_Height/2

            # Adds a Eigenmode Expansion (EME) solver region to the MODE simulation environment.
            self.lum.addeme()
            self.lum.set("x min", -(Device_Length / 2)+0.1e-6)
            self.lum.set("y", 0)
            self.lum.set('simulation temperature', 273.15 + 20)
            self.lum.set("y span", Device_Width)
            self.lum.set("z", Substrate_Height)
            self.lum.set("z span", 4e-6)
            self.lum.set("wavelength", WaveLength)
            self.lum.set("z min bc", "PML")
            self.lum.set("z max bc", "PML")
            self.lum.set("y min bc", "PML")
            self.lum.set("y max bc", "PML")


            # set cell properties
            self.lum.set("number of cell groups", 3)
            self.lum.set("group spans", np.array([[WG_Length-0.1e-6], [MMI_Length], [WG_Length-0.1e-6]]))
            self.lum.set("cells", np.array([[20], [3], [20]]))
            self.lum.set("subcell method", np.array([[0], [0], [0]]))

        elif Taper == True:

            # Device specifications
            Device_Length = MMI_Length + 2 * WG_Length + 2 * TaperLength
            Device_Width = MMI_Width + WaveLength * 2  # MMI_Width
            max_slabH = Slab_Height
            # MonitorHeight = Substrate_Height + (max_slabH + WG_Height) / 2
            MonitorHeight = Substrate_Height + max_slabH + WG_Height/2
            # Ports_mid = (max_slabH + WG_Height) / 2
            Ports_mid = max_slabH + WG_Height/2

            # Adds a Eigenmode Expansion (EME) solver region to the MODE simulation environment.
            self.lum.addeme()
            self.lum.set("x min", -(Device_Length / 2)+0.1e-6)
            self.lum.set("y", 0)
            self.lum.set('simulation temperature', 273.15 + 20)
            self.lum.set("y span", Device_Width)
            self.lum.set("z", Substrate_Height)
            self.lum.set("z span", 4e-6)
            self.lum.set("wavelength", WaveLength)
            self.lum.set("z min bc", "PML")
            self.lum.set("z max bc", "PML")
            self.lum.set("y min bc", "PML")
            self.lum.set("y max bc", "PML")

            Device_Length = MMI_Length + 2 * WG_Length + 2 * TaperLength
            # set cell properties
            self.lum.set("number of cell groups", 5)
            self.lum.set("group spans",
                         np.array([[WG_Length-0.1e-6], [TaperLength], [MMI_Length], [TaperLength], [WG_Length-0.1e-6]]))
            self.lum.set("cells", np.array([[3], [20], [3], [20], [3]]))
            self.lum.set("subcell method", np.array([[0], [1], [0], [1], [0]]))
        else:
            raise ValueError("Incorect Taper variable!. Possible Taper values are Taper = False or Taper = True.")

        # Modes to Calculate
        self.lum.set('number of modes for all cell groups', 20)

        # Mesh Cells
        self.lum.set("define y mesh by", "maximum mesh step")
        self.lum.set("dy", y_res)
        self.lum.set("define z mesh by", "maximum mesh step")
        self.lum.set("dz", z_res)
        self.lum.set('fit materials with multi-coefficient model', 1)
        self.lum.set('wavelength start', 0.4e-6)
        self.lum.set('wavelength stop', 2e-6)

        # define two values for upper and lower limit of the WG offsets
        OffMin = -MMI_Width / 2
        OffMax = MMI_Width / 2

        if posOffset <= OffMin or posOffset >= OffMax:
            self.lum.deleteall()
            raise ValueError('You are Trying to move the Waveguide outside the MMI. This is not possible!')
        else:
            # Define Ports
            yPort_vec = [posOffset / 2 + WG_Width / 2, -(posOffset / 2 + WG_Width / 2), posOffset / 2 + WG_Width / 2, -(posOffset / 2 + WG_Width / 2)]
            portLoc = ["left", "left", "right", "right"]


            overLapp = yPort_vec[0] - y_Port_Span/2
            if overLapp < 0:
                raise ValueError("!!! CAUTION !!! - The Ports are overlapping at the middle! Please change the Y Port Span or move the Waveguides away from each other!")
            else:
                pass

            for i in range(2):
                self.lum.addemeport()

            for i in range(4):
                self.lum.select("EME::Ports::port_" + str(i + 1))
                self.lum.set("port location", portLoc[i])
                self.lum.set("use full simulation span", 0)
                self.lum.set('y', yPort_vec[i])
                self.lum.set('y span', y_Port_Span)
                self.lum.set("z", Ports_mid)
                self.lum.set("z span", z_Port_Span)
                self.lum.set("mode selection", Mode)

            # Add monitor
            self.lum.addemeprofile()
            self.lum.set("x min", -(Device_Length / 2))
            self.lum.set("x max", (Device_Length / 2))
            self.lum.set("y min", -Device_Width / 2)
            self.lum.set("y max", Device_Width / 2)
            self.lum.set("z", MonitorHeight)






    def setWDMEMESolver(self, Parameters):
        '''

        Parameters
        ----------
        Parameters : Dictionary
            Dictionary with all the data needet for the WDMEMESolver.
            Parameters['Substrate Height'] : int/float
                Substrate height.
            Parameters['MMI Width'] : int/float
                Width of the MMI.
            Parameters['MMI Length'] : int/float
                Length of the MMI.
            Parameters['WG Height' : int/float
                Waveguide hight. Also the height of the MMI section
            Parameters['WG Length'] : int/float
                Waveguide length.
            Parameters['WG Width'] : int/float
                Waveguide width.
            Parameters['Taper Width'] : int/float
                Taper backside width, frontside width is the waveguide width.
            Parameters['Taper Length'] : int/float
                Taper Length
            Parameters['y res'] : int/float
                Mesh y-Axis
            Parameters['z res'] : int/float
                Mesh z-Axis
            Parameters['Slab Height'] : int/float
                Slab Height.
            Parameters['Wavelength'] : int/float
                Wavelength
            Parameters['Angle Thetha'] : boolen
                Angle for the input and output waveguides
            Parameters["Mode"] : str
                Mode to choose from ("fundamental TE mode", "fundamental TM mode", "fundamental mode")
            Parameters["Port Span"]: list of int/floats
                Parameters["Port Span"] = [Span of Port in x direction, Span of Port in y direction, Span of Port in z direction]

            Raises
            ------
            ValueError
                DESCRIPTION.

            Returns
            -------
            None.

                '''



        Substrate_Height = Parameters['Substrate Height']
        MMI_Width = Parameters['MMI Width']
        MMI_Length = Parameters['MMI Length']
        WG_Height = Parameters['WG Height']
        WG_Length = Parameters['WG Length']
        WG_Width = Parameters['WG Width']
        y_res = Parameters['y res']
        z_res = Parameters['z res']
        Slab_Height = Parameters['Slab Height']
        WaveLength = Parameters['Wavelength']
        Mode = Parameters["Mode"]
        angleTheta = Parameters['Angle Thetha']
        TaperWidth = Parameters['Taper Width']
        y_Port_Span = Parameters["Port Span"][1]
        z_Port_Span = Parameters["Port Span"][2]
        TaperLength = Parameters['Taper Length']


        # Device specifications
        Device_Length = MMI_Length + 2 * WG_Length
        # Device_Width = MMI_Width + 2*WG_Length + WaveLength * 2  # MMI_Width
        Device_Width = MMI_Width + 2*WG_Length + WaveLength * 2
        max_slabH = Slab_Height
        # MonitorHeight = Substrate_Height + (max_slabH + WG_Height) / 2
        MonitorHeight = Substrate_Height + max_slabH + WG_Height/2
        # Ports_mid = (max_slabH + WG_Height) / 2
        Ports_mid = max_slabH + WG_Height/2


        # EME Boundary Length
        BoardLen = np.cos(angleTheta * np.pi / 180) * TaperLength
        X_min = -BoardLen - MMI_Length/2

        # Adds a Eigenmode Expansion (EME) solver region to the MODE simulation environment.
        self.lum.addeme()
        self.lum.set('simulation temperature', 273.15 + 20)
        self.lum.set("y", 0)
        self.lum.set("y span", Device_Width)
        self.lum.set("x min", X_min + 0.1e-6)
        self.lum.set("z", Substrate_Height)
        self.lum.set("z span", 4e-6)
        self.lum.set("wavelength", WaveLength)
        self.lum.set("z min bc", "PML")
        self.lum.set("z max bc", "PML")
        self.lum.set("y min bc", "PML")
        self.lum.set("y max bc", "PML")
        # set cell properties
        self.lum.set("number of cell groups", 3)
        self.lum.set("group spans", np.array([[BoardLen ], [MMI_Length], [BoardLen ]]))
        self.lum.set("cells", np.array([[15], [10], [15]]))
        self.lum.set("subcell method", np.array([[1], [1], [1]]))

        # Modes to Calculate
        self.lum.set('number of modes for all cell groups', 20)

        # Mesh Cells
        self.lum.set("define y mesh by", "maximum mesh step")
        self.lum.set("dy", y_res)
        self.lum.set("define z mesh by", "maximum mesh step")
        self.lum.set("dz", z_res)
        self.lum.set('fit materials with multi-coefficient model', 1)
        self.lum.set('wavelength start', 0.4e-6)
        self.lum.set('wavelength stop', 2e-6)

        # max_yPos = [(WG_W / 2 + OffsetInput) * 2, 0, (WG_Width / 2 + posOffset / 2) * 2]
        # min_yPos = [-(WG_W / 2 + OffsetInput) * 2, -(WG_Width / 2 + posOffset / 2) * 2, 0]


        if angleTheta <= 20:
            Input_yPos = -MMI_Width / 2 +  TaperWidth/2   # Correction factor so that Ports are in the WG 0.1e-6
            Input_Y = Input_yPos - (np.sqrt((TaperLength)**2 - BoardLen**2))


            Output_yPos = MMI_Width / 2 -TaperWidth/2
            Output_Y = Output_yPos + (np.sqrt((TaperLength)**2 - BoardLen**2))
            portLoc = ["left", "right"]




            self.lum.select("EME::Ports::port_" + str(1))
            self.lum.set("port location", portLoc[0])
            self.lum.set("use full simulation span", 0)
            self.lum.set("y", Input_Y)
            self.lum.set("y span", y_Port_Span )
            self.lum.set("z", Ports_mid)
            self.lum.set("z span", z_Port_Span)
            self.lum.set("mode selection", Mode)
            self.lum.set("theta", angleTheta)

            self.lum.select("EME::Ports::port_" + str(2))
            self.lum.set("port location", portLoc[1])
            self.lum.set("use full simulation span", 0)
            self.lum.set("y", Output_Y)
            self.lum.set("y span", y_Port_Span )
            self.lum.set("z", Ports_mid)
            self.lum.set("z span", z_Port_Span)
            self.lum.set("mode selection", Mode)
            self.lum.set("theta", angleTheta)

            # Add monitor
            self.lum.addemeprofile()
            self.lum.set("x", 0)
            self.lum.set("x span", Device_Length)
            # self.lum.set("x min", -(Device_Length / 2))
            # self.lum.set("x max", (Device_Length / 2))

            self.lum.set("y min", -Device_Width / 2)
            self.lum.set("y max", Device_Width / 2)
            self.lum.set("z", MonitorHeight)


        else:

            # Correction in Y-Axis
            difY = (TaperLength / 2) * np.cos(angleTheta * np.pi / 180)
            NewY = -MMI_Width / 2 + TaperWidth / 2 - np.sqrt((TaperLength / 2) ** 2 - difY ** 2)  # - xLen - Diff/2 #

            Input_yPos = -MMI_Width / 2 +  TaperWidth/2
            Corr = TaperLength *np.sin((angleTheta*np.pi/180)/2)
            Input_Y = NewY - difY/2


            Output_yPos = MMI_Width / 2 - TaperWidth/2
            Output_Y = Output_yPos + (np.sqrt((TaperLength)**2 - BoardLen**2)) - WG_Width
            Output_Y = -NewY + difY / 2
            portLoc = ["left", "right"]




            self.lum.select("EME::Ports::port_" + str(1))
            self.lum.set("port location", portLoc[0])
            self.lum.set("use full simulation span", 0)
            self.lum.set("y", Input_Y)
            self.lum.set("y span", y_Port_Span*3 )
            self.lum.set("z", Ports_mid)
            self.lum.set("z span", z_Port_Span)
            self.lum.set("mode selection", Mode)
            self.lum.set("theta", angleTheta)

            self.lum.select("EME::Ports::port_" + str(2))
            self.lum.set("port location", portLoc[1])
            self.lum.set("use full simulation span", 0)
            self.lum.set("y", Output_Y)
            self.lum.set("y span", y_Port_Span*3 )
            self.lum.set("z", Ports_mid)
            self.lum.set("z span", z_Port_Span)
            self.lum.set("mode selection", Mode)
            self.lum.set("theta", angleTheta)

            # Add monitor
            self.lum.addemeprofile()
            self.lum.set("x", 0)
            self.lum.set("x span", Device_Length)
            # self.lum.set("x min", -(Device_Length / 2))
            # self.lum.set("x max", (Device_Length / 2))

            self.lum.set("y min", -Device_Width / 2)
            self.lum.set("y max", Device_Width / 2)
            self.lum.set("z", MonitorHeight)







    def setInverseTaperEMESolver(self, Parameters):
        '''
          Parameters
          ----------
          Parameters : Dictionary
              Dictionary with all the data needet for the InverseTaperEMESolver.
              Parameters['Substrate Height'] : int/float
                  Substrate height.
              Parameters['WG Height' : int/float
                  Waveguide hight. Also the height of the MMI section
              Parameters['WG Width'] : int/float
                  Waveguide width.
              Parameters['Slab Height'] : int/float
                  Slab height
              Parameters['PWB Taper Width Back'] : int/float
                  Photonic Wirebonding (PWB) Width back side (to the Photonic Wire Bonding)
              Parameters['PWB Taper Hight Back'] : int/float
                  Photonic Wire Bonding Height back side (to the Photonic Wire Bonding)
              Parameters['PWB Taper Length'] : int/float
                  Length of the Photonic Wire Bonding Taper
              Parameters["SMF Core Diameter"] : int/float
                Single Mode Fiber core Diameter
              Parameters["SMF Cladding Diameter"] : int/float
                Single Mode Fiber Cladding Diameter
              Parameters['y res'] : int/float
                  Mesh y-Axis
              Parameters['z res'] : int/float
                  Mesh z-Axis
              Parameters['Wavelength'] : int/float
                  Wavelength
              Parameters["Mode"] : str
                  Mode to choose from ("fundamental TE mode", "fundamental TM mode", "fundamental mode")
              Parameters["Port Span"] : list of floats/ints
                  List of x,y and z span of the Ports. For this simulation only y and z parametes will be taken.
          

        Returns
        -------
        None.

        '''
    
    
        Substrate_Height = Parameters['Substrate Height']
        WG_Height = Parameters['WG Height']
        WG_Width = Parameters['WG Width']
        Slab_Height = Parameters['Slab Height']
        TaperWidthB = Parameters['PWB Taper Width Back']
        TaperHightB = Parameters['PWB Taper Hight Back']
        TaperLength_PWB = Parameters['PWB Taper Length']
        WaveLength = Parameters['Wavelength']
        Mode = Parameters["Mode"]
        y_res = Parameters['y res']
        z_res = Parameters['z res']
        y_Port_Span = Parameters["Port Span"][1]
        z_Port_Span = Parameters["Port Span"][2]
        
        # SMF Parameters
        CoreDiameter = Parameters["SMF Core Diameter"]
        CladdingDiameter = Parameters["SMF Cladding Diameter"]


        if Slab_Height == 0:
            # # Device specifications
            MonitorHeight = Substrate_Height/2 + WG_Height/2
            Ports_mid = -(Substrate_Height/2 + CoreDiameter/2) + MonitorHeight
            Ports_PWB_mid = CoreDiameter/2 # TaperHightB/2
            self.lum.select("SMF")
            xPos_SMF = self.lum.get("x")
            X_min = -TaperLength_PWB/2 - abs(xPos_SMF)
            
        else:
            # # Device specifications
            max_slabH = Slab_Height
            MonitorHeight = Substrate_Height/2 + max_slabH + WG_Height/2
            Ports_mid = -(Substrate_Height/2 + CoreDiameter/2) + MonitorHeight
            #Ports_mid = max_slabH + Substrate_Height/2 - CoreDiameter/2 + WG_Width/2
            Ports_PWB_mid = max_slabH + CoreDiameter/2 # TaperHightB/2
            self.lum.select("SMF")
            xPos_SMF = self.lum.get("x")
            X_min = -TaperLength_PWB/2 - abs(xPos_SMF)
        



        # Adds a Eigenmode Expansion (EME) solver region to the MODE simulation environment.
        self.lum.addeme()
        self.lum.set("x min", X_min)
        self.lum.set("y", 0)
        self.lum.set("y span", TaperWidthB + TaperWidthB/2)
        self.lum.set('simulation temperature', 273.15 + 20)
        self.lum.set("z", Substrate_Height/2 + CoreDiameter/2 ) #Substrate_Height
        self.lum.set("z span", TaperHightB*2)
        self.lum.set("wavelength", WaveLength)
        self.lum.set("z min bc", "PML")
        self.lum.set("z max bc", "PML")
        self.lum.set("y min bc", "PML")
        self.lum.set("y max bc", "PML")
        # set cell properties
        self.lum.select("InverseTaper::WG_Extention_Inverse_Taper")
        ExtWG_Lenght = self.lum.get("poles")[1][0]
        self.lum.select("EME")
        self.lum.set("number of cell groups", 3)
        self.lum.set("group spans", np.array([[0.1e-6], [TaperLength_PWB], [ExtWG_Lenght]]))
        self.lum.set("cells", np.array([[3], [80], [1]]))
        self.lum.set("subcell method", np.array([[1], [1], [1]]))

        # Modes to Calculate
        self.lum.set('number of modes for all cell groups', 20)

        # Mesh Cells
        self.lum.set("define y mesh by", "maximum mesh step")
        self.lum.set("dy", y_res)
        self.lum.set("define z mesh by", "maximum mesh step")
        self.lum.set("dz", z_res)
        self.lum.set('fit materials with multi-coefficient model', 1)
        self.lum.set('wavelength start', 0.4e-6)
        self.lum.set('wavelength stop', 2e-6)


        portLoc = ["left", "right"]
        # z_span_PWB = CoreDiameter/2 + (z_Port_Span - WG_Height)/2

        self.lum.select("EME::Ports::port_" + str(1))
        self.lum.set("port location", portLoc[0])
        # self.lum.set("use full simulation span", 1)
        self.lum.set("use full simulation span", 0)
        self.lum.set("y", 0)
        self.lum.set("y span", CoreDiameter + 2e-6)
        self.lum.set("z", Substrate_Height/2)
        self.lum.set("z span", CoreDiameter + 2e-6 )
        self.lum.set("mode selection", Mode)


        self.lum.select("EME::Ports::port_" + str(2))
        self.lum.set("port location", portLoc[1])
        # self.lum.set("use full simulation span", 1)
        self.lum.set("use full simulation span", 0)
        self.lum.set("y", 0)
        self.lum.set("y span", y_Port_Span)
        self.lum.set("z", Ports_mid)
        self.lum.set("z span", z_Port_Span)
        self.lum.set("mode selection", Mode)



        # Add monitor Horizondal
        self.lum.addemeprofile()
        # self.lum.set("x", 0)
        # self.lum.set("x span", TaperLength_PWB+11e-6)
        # self.lum.set("x min", X_min - 5e-6)
        # self.lum.set("x max", TaperLength_PWB/2 + 15e-6)
        self.lum.set("x min", -TaperLength_PWB/2)
        self.lum.set("x max", TaperLength_PWB )
        self.lum.set("y", 0)
        self.lum.set("y span", TaperWidthB*2)
        self.lum.set("z", MonitorHeight)


        # Add monitor Vertical
        self.lum.addemeprofile()
        self.lum.set('monitor type', "2D Y-normal")
        # self.lum.set("x", 0)
        # self.lum.set("x span", TaperLength_PWB+11e-6)
        self.lum.set("x min", -TaperLength_PWB/2)
        self.lum.set("x max", TaperLength_PWB )
        # self.lum.set("x min", X_min - 5e-6)
        # self.lum.set("x max", TaperLength_PWB/2 + 15e-6)
        self.lum.set("y", 0)
        # self.lum.set("y span", TaperWidthB*2)
        self.lum.set("z", MonitorHeight)
        self.lum.set("z span", TaperHightB*2)





    def setWaveguideFDESolver(self, Parameters):
        '''


        Parameters
        ----------
        Parameters : Dictionary
            Dictionary with all the data needet for the WaveguideFDESolver.
            Parameters : Dictionary
                Dictionary with all the data needet for the Bend Wavaguide. Data needet:
            Parameters['Substrate Height']: int/float   
                Substrate height
            Parameters['WG Height'] : int/float
                Heigth of waveguide
            Parameters['y res'] : int/float
                Mesh cell sizes.
            Parameters['z res'] : int/float
                Mesh cell sizes.
            Parameters['Slab Height'] : int/float
                Height of the slab.
            Parameters['Wavelength'] : int/float
                Wavelength.

        Returns
        -------
        None.

        '''

        WG_Height = Parameters['WG Height']
        dy = Parameters['y res']
        dz = Parameters['z res']
        Slab_Height = Parameters['Slab Height']
        SubstrateHight = Parameters['Substrate Height']
        WaveLength = Parameters['Wavelength']
        

        # FDE Dimensions
        WG_Mid = Slab_Height + WG_Height / 2

        # Adds a Eigenmode Expansion (EME) solver region to the MODE simulation environment.
        self.lum.addfde()
        self.lum.set('solver type', "2D X normal")
        self.lum.set("x", 0)
        self.lum.set('define y mesh by', 'maximum mesh step')
        self.lum.set('dy', dy)
        self.lum.set('define z mesh by', 'maximum mesh step')
        self.lum.set('dz', dz)
        self.lum.set("z min bc", "PML")
        self.lum.set("z max bc", "PML")
        self.lum.set("y min bc", "PML")
        self.lum.set("y max bc", "PML")
        self.lum.set('fit materials with multi-coefficient model', 1)
        self.lum.set('wavelength start', 0.4e-6)
        self.lum.set('wavelength stop', 2e-6)
        self.lum.set("y max", 5e-6)
        self.lum.set("y min", -5e-6)
        self.lum.set('z', WG_Mid)
        self.lum.set('z span', WG_Height + Slab_Height + SubstrateHight/2)
        self.lum.set("wavelength", WaveLength)




    def setGratingCouplerNeffFDE(self, Parameters):

        GC_Height = Parameters["Hight GC"]
        Pitch = Parameters["Pitch GC"]
        dy = Parameters['y res']
        dz = Parameters['z res']
        WaveLength = Parameters['Wavelength']

        # FDE Dimensions
        Gap = 1e-6 - Pitch
        WG_Mid = GC_Height / 2


        # Adds a Eigenmode Expansion (EME) solver region to the MODE simulation environment.
        self.lum.addfde()
        self.lum.set('solver type', "2D X normal")
        self.lum.set("x", 0)
        self.lum.set('define y mesh by', 'maximum mesh step')
        self.lum.set('dy', dy)
        self.lum.set('define z mesh by', 'maximum mesh step')
        self.lum.set('dz', dz)
        self.lum.set("z min bc", "PML")
        self.lum.set("z max bc", "PML")
        self.lum.set("y min bc", "PML")
        self.lum.set("y max bc", "PML")
        self.lum.set('fit materials with multi-coefficient model', 1)
        self.lum.set('wavelength start', 0.4e-6)
        self.lum.set('wavelength stop', 2e-6)
        self.lum.set("y", 0)
        self.lum.set("y span", Pitch+Gap)
        self.lum.set('z', WG_Mid)
        self.lum.set('z span',1e-6)
        self.lum.set("wavelength", WaveLength)




    # =============================================================================
    # Functions
    # =============================================================================


    def StartEMESimulation(self):
        '''
        This Function will save a the construted object under the name "SimRun1"
        then will run the simulation. The simulation cannot be done without first
        saving the object. On the end the 'emepropagation' of the mode 1 port 1
        will be done.

        Returns
        -------
        None.

        '''
        self.lum.save('SimRun1')
        self.lum.run()
        self.lum.emepropagate()




    def StartFDTDSimulation(self):
        '''
        This Function will save a the construted object under the name "SimRun1"
        then will run the simulation. The simulation cannot be done without first
        saving the object. On the end the 'emepropagation' of the mode 1 port 1
        will be done.

        Returns
        -------
        None.

        '''
        self.lum.save('SimRun1')
        self.lum.run()
        
        
       
    def StartFDTDOptimizerRingGratingCoupler(self, Tranmission_on_Port, pos):
   
        self.lum.select("SMF")
        self.lum.set("x", pos)
        self.lum.select("FDTD::ports::Input_SMF_Port")
        self.lum.set("x", pos)
        self.lum.select("Power_Input_SMF_Port")
        self.lum.set("x", pos)
        self.lum.save('SimRun1')
        self.lum.run()
        data = np.abs(self.lum.getresult("FDTD::ports::" + str(Tranmission_on_Port), 'T')["T"])
        self.lum.switchtolayout()
        return data
  



    def StartFDESimulation(self):
        '''
        This Function will save a the construted object under the name "SimRun1"
        then will run the simulation. The simulation cannot be done without first
        saving the object. On the end the 'findmodes' will simulate all 20 modes
        that are set by defoult to be calculated.

        Returns
        -------
        None.

        '''
        self.lum.save('SimRun1')
        self.lum.run()
        self.lum.findmodes()
        




    def ExtractFDTDResults(self, Ports, Sparam):

        '''
        This fuction will add an S-Parameter Sweep to the FDTD simulation.
        This will allowed to sweep the S-parameter Matrix like it is done in
        EMe solver. Additionally will extract the Power on each port according to the time
        monitors that ware set in the solver functions.

        '''
        if Sparam == True:
            self.lum.addsweep(3)
            # Check "Excite all ports" option
            self.lum.setsweep("s-parameter sweep", "Excite all ports", 1)
            self.lum.save('SimRun1')
            # run s-parameter sweep
            self.lum.runsweep("s-parameter sweep")

            OptionsEME1 = {}
            Power = {}
            OptionsEME1 = self.lum.getsweepresult("s-parameter sweep", "S matrix")
            self.lum.deletesweep("s-parameter sweep")

            if Ports == 3:
                Power['Input'] = 0
                Power['Output_L'] = 0
                Power['Output_R'] = 0
                Power["Transmission Input"] = 0
                Power["Transmission Output_L"] = 0
                Power["Transmission Output_R"] = 0
                names = ['Power_Input', 'Power_Output_L', 'Power_Output_R']
                for i in range(len(names)):
                    Power[names[i]] = self.lum.getresult(names[i], 'P')
                for i in range(len(names)):
                    Power["Transmission " + names[i]] = self.lum.getresult(names[i], 'T')

            elif Ports == 2:
                Power["Input"] = 0
                Power["Output"] = 0
                Power["Transmission Input"] = self.lum.getresult('Input', 'T')
                Power["Transmission Output"] = self.lum.getresult('Output', 'T')
                names = ['Input', 'Output']
                for i in range(len(names)):
                    Power[names[i]] = self.lum.getresult(names[i], 'P')

            elif Ports == 4:
                Power['Input_L'] = 0
                Power['Input_R'] = 0
                Power['Output_L'] = 0
                Power['Output_R'] = 0
                Power["Transmission Input_L"] = 0
                Power["Transmission Input_R"] = 0
                Power["Transmission Output_L"] = 0
                Power["Transmission Output_R"] = 0
                names = ['Power_Input_L', 'Power_Input_R', 'Power_Output_L', 'Power_Output_R']
                for i in range(len(names)):
                    Power[names[i]] = self.lum.getresult(names[i], 'P')
                for i in range(len(names)):
                    Power["Transmission "+ names[i]] = self.lum.getresult(names[i], 'T')

            else:
                raise ValueError("Incorect Port number in ExtractFDTD() Function. Possible Port Numbers are int: 3 or 4! ")
            # S21 = abs(OptionsEME1['S'][0][2])**2
            # S31 = abs(OptionsEME1['S'][0][1])**2
            return OptionsEME1, Power

        else:
            Power = {}
            if Ports == 3:
                Power['Input'] = 0
                Power['Output_L'] = 0
                Power['Output_R'] = 0
                Power["Transmission Input"] = 0
                Power["Transmission Output_L"] = 0
                Power["Transmission Output_R"] = 0
                names = ['Power_Input', 'Power_Output_L', 'Power_Output_R']
                for i in range(len(names)):
                    Power[names[i]] = self.lum.getresult(names[i], 'P')
                for i in range(len(names)):
                    Power["Transmission " + names[i]] = self.lum.getresult(names[i], 'T')

            elif Ports == 2:
                Power["Input"] = 0
                Power["Output"] = 0
                Power["Transmission Input"] = self.lum.getresult('Input', 'T')
                Power["Transmission Output"] = self.lum.getresult('Output', 'T')
                names = ['Input', 'Output']
                for i in range(len(names)):
                    Power[names[i]] = self.lum.getresult(names[i], 'P')

            elif Ports == 4:
                Power['Input_L'] = 0
                Power['Input_R'] = 0
                Power['Output_L'] = 0
                Power['Output_R'] = 0
                Power["Transmission Input_L"] = 0
                Power["Transmission Input_R"] = 0
                Power["Transmission Output_L"] = 0
                Power["Transmission Output_R"] = 0
                names = ['Power_Input_L', 'Power_Input_R', 'Power_Output_L', 'Power_Output_R']
                for i in range(len(names)):
                    Power[names[i]] = self.lum.getresult(names[i], 'P')
                for i in range(len(names)):
                    Power["Transmission " + names[i]] = self.lum.getresult(names[i], 'T')

            else:
                raise ValueError(
                    "Incorect Port number in ExtractFDTD() Function. Possible Port Numbers are int: 3 or 4! ")
            return 0, Power





    def ExtractEMEResults(self, MonitorName, EMEName):
        '''


        Parameters
        ----------
        MonitorName : str
            The name of the monitor. By defoult the name will be 'monitor'.
            So MonitorName = 'monitor' can be used. In case of change in the name
            the new name of the monitor should be given.
        EMEName : str
            The name of the solver. In this case per defoult will be EMEName =
            'EME'. In case of change in the name, the new name should be given.

        Returns
        -------
        OptionsMonitor : dict
            Dictionary with:
                    'lambda'
                    'f'
                    'x'
                    'y'
                    'z'
                    'E'
                    'H'
                    'Lumerical_dataset'
        EME_Data : dict
            Dictionary with:
                    'user s matrix'
                    'internal s matrix'
                    'local diagnostics'
                    'coefficients'
                    'global diagnostics'

        '''

        # Define List with all the data from monitor and EME
        EME_Data = {}

        # Check all avaliable variables for the monitor and EME
        OptionsMonitor = self.lum.getresult(str(MonitorName), 'field profile')
        OptionsEME1 = self.lum.getresult(str(EMEName))
        OptionsEME = OptionsEME1.split('\n')

        for data in OptionsEME:
            EME_Data[data] = self.lum.getresult(str(EMEName), str(data))

        return OptionsMonitor, EME_Data




    def ExtractFDEModes(self, EffIndexValue):
        '''


        Parameters
        ----------
        EffIndexValue : float/int
            Effective index value. All modes with effective index smaller then
            the  Effective will be deleted. Effective is usualy the smalles effective
            index from the hole construction (usualy the cladding).

        Returns
        -------
        dictData : dict
                 Dictionary with the modes polarization numbers

        '''
        dictModes = {}
        dictModes2 = {}
        modeList = self.lum.nummodes()

        if modeList != 0:
            for i in range(1, int(modeList + 1)):
                if self.lum.getdata('FDE::data::mode' + str(i),
                                    'TE polarization fraction').real * 100 > 50 and self.lum.getdata(
                    'FDE::data::mode' + str(i), 'neff').real > EffIndexValue:
                    dictModes['mode' + str(i)] = self.lum.getdata('FDE::data::mode' + str(i),
                                                                  'TE polarization fraction') * 100
                elif self.lum.getdata('FDE::data::mode' + str(i),
                                      'TE polarization fraction').real * 100 < 50 and self.lum.getdata(
                    'FDE::data::mode' + str(i), 'neff').real > EffIndexValue:
                    dictModes2['mode' + str(i)] = self.lum.getdata('FDE::data::mode' + str(i),
                                                                   'TE polarization fraction') * 100
                else:
                    pass

        return dictModes, dictModes2



    def CoppyDcard(self, modeTE, modeTM):
        '''


        Parameters
        ----------
        modeTE : str
            The TE mode that will be copy to dCard to do the overlap analysis.

        modeTM : str
           The TM mode that will be copy to dCard to do the overlap analysis.

        Returns
        -------
        None.

        '''
        self.lum.copydcard(str(modeTE), 'TE')
        self.lum.copydcard(str(modeTM), 'TM')




    def Overlap(self):
        '''


        Returns
        -------
        TEModes : dict
            Dictionary of TE - Modes.
        dictDataTE : dict
            Dictionary of TE - Modes E field values
        SweepDictTE : dict
            Dictionary of TE - Modes:
                    1) overlap Procentage
                    2) TE polarization fraction number
                    3) TE polarization fraction
                    4) loss
                    5) Ex
                    6) Ey
                    7) Ez

        TMModes : dict
            Dictionary of TM - Modes..
        dictDataTM : dict
            Dictionary of TM - Modes E field values.
        SweepDictTM : dict
            Dictionary of TM - Modes:
                    1) overlap Procentage
                    2) TE polarization fraction number
                    3) TE polarization fraction
                    4) loss
                    5) Ex
                    6) Ey
                    7) Ez

        '''
        numModes = self.lum.nummodes()

        TEModes = {}
        dictDataTE = {}
        SweepDictTE = {}

        TMModes = {}
        dictDataTM = {}
        SweepDictTM = {}

        for i in range(1, int(numModes + 1)):
            if self.lum.overlap("TE", 'mode' + str(i))[0] > 0.80:
                SweepDictTE['mode ' + str(i) + ' overlap Procentage'] = self.lum.overlap("TE", 'mode' + str(i))[0] * 100
                SweepDictTE['mode ' + str(i) + ' TE polarization fraction num'] = self.lum.getdata(
                    'FDE::data::mode' + str(i), 'TE polarization fraction')
                TEModes['mode ' + str(i) + ' effective index num'] = self.lum.getdata('FDE::data::mode' + str(i),
                                                                                      'neff')
                SweepDictTE['mode ' + str(i) + ' effective index num'] = self.lum.getdata('FDE::data::mode' + str(i),
                                                                                          'neff')
                SweepDictTE['mode ' + str(i) + ' loss'] = self.lum.getdata('FDE::data::mode' + str(i), 'loss') / 100
                SweepDictTE['mode ' + str(i) + ' Ex'] = self.lum.getdata('FDE::data::mode' + str(i), 'Ex')
                SweepDictTE['mode ' + str(i) + ' Ey'] = self.lum.getdata('FDE::data::mode' + str(i), 'Ey')
                SweepDictTE['mode ' + str(i) + ' Ez'] = self.lum.getdata('FDE::data::mode' + str(i), 'Ez')
                dictDataTE['mode ' + str(i) + ' z'] = self.lum.getresult('FDE::data::mode' + str(i), 'E')
            else:
                pass

        for i in range(1, int(numModes + 1)):
            if self.lum.overlap("TM", 'mode' + str(i))[0] > 0.80:
                SweepDictTM['mode ' + str(i) + ' overlap Procentage'] = self.lum.overlap("TM", 'mode' + str(i))[0] * 100
                SweepDictTM['mode ' + str(i) + ' TE polarization fraction num'] = self.lum.getdata(
                    'FDE::data::mode' + str(i), 'TE polarization fraction')
                TMModes['mode ' + str(i) + ' effective index num'] = self.lum.getdata('FDE::data::mode' + str(i),
                                                                                      'neff')
                SweepDictTM['mode ' + str(i) + ' effective index num'] = self.lum.getdata('FDE::data::mode' + str(i),
                                                                                          'neff')
                SweepDictTM['mode ' + str(i) + ' loss'] = self.lum.getdata('FDE::data::mode' + str(i), 'loss') / 100
                SweepDictTM['mode ' + str(i) + ' Ex'] = self.lum.getdata('FDE::data::mode' + str(i), 'Ex')
                SweepDictTM['mode ' + str(i) + ' Ey'] = self.lum.getdata('FDE::data::mode' + str(i), 'Ey')
                SweepDictTM['mode ' + str(i) + ' Ez'] = self.lum.getdata('FDE::data::mode' + str(i), 'Ez')
                dictDataTM['mode ' + str(i) + ' z'] = self.lum.getresult('FDE::data::mode' + str(i), 'E')
            else:
                pass

        return TEModes, dictDataTE, SweepDictTE, TMModes, dictDataTM, SweepDictTM




    def ExtractFDEResultsExtendet(self, EffIndexValue):
        """


        Parameters
        ----------
        EffIndexValue : float/int
            Effective index value. All modes with effective index smaller then
            the  Effective will be deleted. Effective is usualy the smalles effective
            index from the hole construction (usualy the cladding).

        Returns
        -------
        dictModes : dict
                Dictionary with:
                        'TE polarization fraction num'
                        'effective area'
                        'effective index num'
                        'group index num'
                        'Ex'
                        'Ey'
                        'Ez'

        dictData : dict
                 Dictionary with:
                             'E'

        """

        dictModes = {}
        dictModes2 = {}
        dictData = {}
        dictData2 = {}
        modeList = self.lum.nummodes()

        if modeList != 0:
            for i in range(1, int(modeList + 1)):
                if self.lum.getdata('FDE::data::mode' + str(i),
                                    'TE polarization fraction').real * 100 > 50 and self.lum.getdata(
                    'FDE::data::mode' + str(i), 'neff').real > EffIndexValue:
                    dictModes['mode ' + str(i) + ' TE polarization fraction num'] = self.lum.getdata(
                        'FDE::data::mode' + str(i), 'TE polarization fraction') * 100

                elif self.lum.getdata('FDE::data::mode' + str(i),
                                      'TE polarization fraction').real * 100 < 50 and self.lum.getdata(
                    'FDE::data::mode' + str(i), 'neff').real > EffIndexValue:
                    dictModes2['mode ' + str(i) + ' TE polarization fraction num'] = self.lum.getdata(
                        'FDE::data::mode' + str(i),'TE polarization fraction') * 100
                else:
                    pass

        if modeList != 0:
            for i in range(1, int(modeList + 1)):
                if self.lum.getdata('FDE::data::mode' + str(i),
                                    'TE polarization fraction').real * 100 > 50 and self.lum.getdata(
                    'FDE::data::mode' + str(i), 'neff').real > EffIndexValue:  # 1.45625:
                    dictData['mode ' + str(i) + ' TE polarization fraction num'] = self.lum.getdata(
                        'FDE::data::mode' + str(i), 'TE polarization fraction')
                    dictData['mode ' + str(i) + ' effective index num'] = self.lum.getdata('FDE::data::mode' + str(i),
                                                                                           'neff')
                    dictData['mode ' + str(i) + ' group index'] = self.lum.getdata('FDE::data::mode' + str(i), 'ng')
                    dictData['mode ' + str(i) + ' loss'] = self.lum.getdata('FDE::data::mode' + str(i), 'loss') / 100
                    dictData['mode ' + str(i) + ' Ex'] = self.lum.getdata('FDE::data::mode' + str(i), 'Ex')
                    dictData['mode ' + str(i) + ' Ey'] = self.lum.getdata('FDE::data::mode' + str(i), 'Ey')
                    dictData['mode ' + str(i) + ' Ez'] = self.lum.getdata('FDE::data::mode' + str(i), 'Ez')
                elif self.lum.getdata('FDE::data::mode' + str(i),
                                      'TE polarization fraction').real * 100 < 50 and self.lum.getdata(
                    'FDE::data::mode' + str(i), 'neff').real > EffIndexValue:
                    dictData2['mode ' + str(i) + ' TE polarization fraction num'] = self.lum.getdata(
                        'FDE::data::mode' + str(i), 'TE polarization fraction')
                    dictData2['mode ' + str(i) + ' effective index num'] = self.lum.getdata('FDE::data::mode' + str(i),
                                                                                            'neff')
                    dictData2['mode ' + str(i) + ' group index'] = self.lum.getdata('FDE::data::mode' + str(i), 'ng')
                    dictData2['mode ' + str(i) + ' loss'] = self.lum.getdata('FDE::data::mode' + str(i), 'loss') / 100
                    dictData2['mode ' + str(i) + ' Ex'] = self.lum.getdata('FDE::data::mode' + str(i), 'Ex')
                    dictData2['mode ' + str(i) + ' Ey'] = self.lum.getdata('FDE::data::mode' + str(i), 'Ey')
                    dictData2['mode ' + str(i) + ' Ez'] = self.lum.getdata('FDE::data::mode' + str(i), 'Ez')
            else:
                pass
        return dictModes, dictData, dictModes2, dictData2




    def Plot3DFDTD(self, Parameters):

        MonitorName = Parameters["Monitor Name FDTD"]
        Field = Parameters["Field FDTD"]

        E = self.lum.getresult(MonitorName, Field)[Field]

        Shape = E.shape

        if Shape[0] == 1:
            Ex = self.lum.getresult(MonitorName, Field)["x"]
            Ey = []
            Ez = []

            for i in range(len(self.lum.getresult(MonitorName, Field)["y"])):
                Ey.append(self.lum.getresult(MonitorName, Field)["y"][i][0] * 1e6)
            for i in range(len(self.lum.getresult(MonitorName, Field)["z"])):
                Ez.append(self.lum.getresult(MonitorName, Field)["z"][i][0] * 1e6)
                
                
            E_a = abs(E) ** 2
            slice_2d = np.sqrt(E_a[0, :, :, 0, 0] ** 2 + E_a[0, :, :, 0, 1] ** 2 + E_a[0, :, :, 0, 2] ** 2)




        elif Shape[1] == 1:
            Ex = []
            Ey = self.lum.getresult(MonitorName, Field)["y"]
            Ez = []

            for i in range(len(self.lum.getresult(MonitorName, Field)["x"])):
                Ex.append(self.lum.getresult(MonitorName, Field)["x"][i][0] * 1e6)
            for i in range(len(self.lum.getresult(MonitorName, Field)["z"])):
                Ez.append(self.lum.getresult(MonitorName, Field)["z"][i][0] * 1e6)
            
            E_a = abs(E) ** 2
            slice_2d = np.sqrt(E_a[:, 0, :, 0, 0] ** 2 + E_a[:, 0, :, 0, 1] ** 2 + E_a[:, 0, :, 0, 2] ** 2)



        elif Shape[2] == 1:
            Ex = []
            Ey = []
            Ez = self.lum.getresult(MonitorName, Field)["z"]

            for i in range(len(self.lum.getresult(MonitorName, Field)["y"])):
                Ey.append(self.lum.getresult(MonitorName, Field)["y"][i][0] * 1e6)
            for i in range(len(self.lum.getresult(MonitorName, Field)["x"])):
                Ex.append(self.lum.getresult(MonitorName, Field)["x"][i][0] * 1e6)
                
                
            E_a = abs(E) ** 2
            slice_2d = np.sqrt(E_a[:, :, 0, 0, 0] ** 2 + E_a[:, :, 0, 0, 1] ** 2 + E_a[:, :, 0, 0, 2] ** 2)

        else:
            raise ValueError("This function can print only 2D fields and no 3D fields")

        # Calc Intensity and murge the 2D monitor data to one vector for printing
        

        # # Rotate the image by 90 degrees counterclockwise
        # slice_2d_rotated = np.rot90(slice_2d, k=3)


        # # Define the y and z axis ranges
        # z_range = np.linspace(-2.68e-6, 50e-6, slice_2d_rotated.shape[0]) * 1e6  # Convert to micrometers
        # y_range = np.linspace(-5e-6, 5e-6, slice_2d_rotated.shape[1]) * 1e6  # Convert to micrometers

        # Create the colormap plot
        fig = plt.figure(figsize=(10, 5))
        plt.pcolormesh( slice_2d, cmap='jet')
        # plt.pcolormesh(y_range, z_range, slice_2d_rotated, cmap='jet')
        plt.colorbar(label='$|'+ Field + '|^2$ Intensity')
        # plt.xlabel('Width / $\mu m$')
        # plt.ylabel('$|E|^2$')
        # plt.title('Lense Thickness = 2.5 $\mu m$')

        if "Save Image" in list(Parameters.keys()):
            # Save the plot as SVG
            plt.savefig(Parameters["Save Image"] + '.svg', format='svg', bbox_inches='tight')
            plt.savefig(Parameters["Save Image"] + '.png', format='png', bbox_inches='tight')
        else:
            pass
        plt.show()

        return fig



    def Plot3DFDE(self, Parameters):
        MonitorName = Parameters["Mode FDE"]
        Field = Parameters["Field FDE"]

        Mode = "TM_LNOI_1x0,5um"
        Ex = abs(self.lum.getdata("FDE::data::" + MonitorName, Field + "x")) ** 2
        Ey = abs(self.lum.getdata("FDE::data::" + MonitorName, Field + "y")) ** 2
        Ez = abs(self.lum.getdata("FDE::data::" + MonitorName, Field + "z")) ** 2

        y = self.lum.getdata("FDE::data::" + MonitorName, "y")[:, 0] * 1e6
        z = self.lum.getdata("FDE::data::" + MonitorName, "z")[:, 0] * 1e6

        slice_2d = np.sqrt(Ex[0, :, :, ] ** 2 + Ey[0, :, :] ** 2 + Ez[0, :, :] ** 2)
        # Rotate the image by 90 degrees counterclockwise
        slice_2d_rotated = np.rot90(slice_2d, k=3)

        # Create the colormap plot
        fig = plt.figure(figsize=(10, 5))
        # plt.pcolormesh(z_range, y_range, slice_2d, cmap='viridis')
        plt.pcolormesh(y, z, slice_2d_rotated[:, :, 0], cmap='jet')
        plt.colorbar(label='$|' + Field + '|^2$')


        if "Save Image" in list(Parameters.keys()):
            # Save the plot as SVG
            plt.savefig(Parameters["Save Image"] + '.svg', format='svg', bbox_inches='tight')
            plt.savefig(Parameters["Save Image"] + '.png', format='png', bbox_inches='tight')
        else:
            pass
        plt.show()

        return fig


class Charge(Constructor):
    def __init__(self, file, Mode, MaterialLib=None):
        '''


        Parameters
        ----------
        file : str
            Path to lumerical lumapi.py.
        Mode : str
            Can only be set to CHARGE. Later on it will be integrated to the rest of the solvers
        MaterialLib : str, optional
            If there is an material library providet it can be imported with passig MateriaLib to the path of where the Material lib is save in the system . The default is None.

        Raises
        ------
        ValueError
            DESCRIPTION.


        '''
        self.file = file
        self.Mode = Mode
        self.MaterialLib = MaterialLib
        self.Materials_Opt = []
        self.Materials_Electric = []

        # Check Python Version !!! No import imp module after python 3.11
        PyVersion = sys.version

        try:
            if PyVersion.split(" ")[0] > "3.11.0":
                import importlib.util
                import importlib.machinery

                def load_source(modname, filename):
                    loader = importlib.machinery.SourceFileLoader(modname, filename)
                    spec = importlib.util.spec_from_file_location(modname, filename, loader=loader)
                    module = importlib.util.module_from_spec(spec)
                    # The module is always executed and not cached in sys.modules.
                    # Uncomment the following line to cache the module.
                    # sys.modules[module.__name__] = module
                    loader.exec_module(module)
                    return module

                self.lumpai = load_source('lumapi', self.file)
            else:
                import imp
                self.lumpai = imp.load_source('lumapi', self.file)
        except ValueError as e:
            print(f"ValueError encountered: {e}")
        except ImportError as e:
            print(f"ImportError encountered: {e}")
        except Exception as e:
            print(f"An unexpected error occurred: {e}")

        self.Mode = Mode
        self.Struct = None
        self.SolverInfo = {}

        if self.Mode == "CHARGE":
            self.CHARGE()
            self.SolverInfo["Solver Used"] = "CHARGE"
            if self.MaterialLib is None:
                pass
            else:
                self.lum.importmaterialdb(self.MaterialLib)
        else:
            raise ValueError(
                "Non Valid Solver was choosen. Please pass on one of the two supported solvers ['FDTD' or 'EME']")

    def CHARGE(self):
        '''

        Calls the CHARGE Solver.

        Returns
        -------
        None.

        '''
        self.lum = self.lumpai.DEVICE()
        print('Lumerical CHARGE API is started')

    # Close Programm
    def Close(self):
        '''


        Returns
        -------
        Close Lumerical GUI

        '''
        self.lum.close()
        print('Lumerical API is closed')

    # Remove the object and solver region
    def removeObject(self):
        '''
        Switch from simulation to layout mode.
        Remove all objects.

        Returns
        -------
        None.

        '''
        self.lum.switchtolayout()
        self.lum.deleteall()

        file = self.file

    def Material(self, Material_Data):
        '''


        Parameters
        ----------
        Material_Data : dictionary
            Material Optical or Electrical data.
            Separate Materials into Electrical and Optical in List Material_Data["Electrical"] = ["mat1", "mat2"..]
            and optical Material_Data["Optical"] = ["omat1", "omat2", ...]

        Returns
        -------
        Materials_Added : list
            List of added materials to the simulation.

        '''

        # Separate Materials into Electrical and Optical in List Material_Data["Electrical"] = ["mat1", "mat2"..]
        # and optical Material_Data["Optical"] = ["omat1", "omat2", ...]
        Electrical_Materials = []
        Optical_Materials = []
        Materials_Added = []

        for i in Material_Data["Electrical"]:
            Electrical_Materials.append(i)
        for j in Material_Data["Optical"]:
            Optical_Materials.append(j)

        # Check for Items that are common for the optical and electrical Material
        Common_Materials = [item for item in Electrical_Materials if item in Optical_Materials]
        Electrical_Materials = [item for item in Electrical_Materials if item not in Common_Materials]
        Optical_Materials = [item for item in Optical_Materials if item not in Common_Materials]

        if Common_Materials:
            # Add all Electrical Materials to materials folder
            for i in range(len(Common_Materials)):
                self.lum.addmodelmaterial()
                self.lum.set("name", Common_Materials[i])
                self.lum.addmaterialproperties("CT", Common_Materials[i])  # importing from electrical material database
                self.lum.select("materials::" + Common_Materials[i])
                self.lum.addmaterialproperties("HT", Common_Materials[i])  # importing from thermal material database
                self.lum.select("materials::" + Common_Materials[i])
                self.lum.addmaterialproperties("EM", Common_Materials[i])
            # Add all Electrical Materials to materials folder
            Materials_Added = Materials_Added + Common_Materials
            for i in range(len(Electrical_Materials)):
                self.lum.addmodelmaterial()
                self.lum.set("name", Electrical_Materials[i])
                self.lum.addmaterialproperties("CT",
                                               Electrical_Materials[i])  # importing from electrical material database
                self.lum.select("materials::" + Electrical_Materials[i])
                self.lum.addmaterialproperties("HT",
                                               Electrical_Materials[i])  # importing from thermal material database
                self.lum.select("materials::" + Electrical_Materials[i])
                for _ in Optical_Materials:
                    if Electrical_Materials[i].split(" ")[0] == _.split(" ")[0]:
                        Optical_Materials.remove(_)
                        self.lum.select("materials::" + Electrical_Materials[i])
                        self.lum.addmaterialproperties("EM", _)
                    else:
                        pass
            Materials_Added = Materials_Added + Electrical_Materials
            # Add all Optical Materials to materials folder
            for i in range(len(Optical_Materials)):
                self.lum.addmodelmaterial()
                self.lum.set("name", Optical_Materials[i])
                self.lum.addmaterialproperties("EM", Optical_Materials[i])
            Materials_Added = Materials_Added + Optical_Materials

        else:
            # Add all Electrical Materials to materials folder
            for i in range(len(Electrical_Materials)):
                self.lum.addmodelmaterial()
                self.lum.set("name", Electrical_Materials[i])
                self.lum.addmaterialproperties("CT",
                                               Electrical_Materials[i])  # importing from electrical material database
                self.lum.select("materials::" + Electrical_Materials[i])
                self.lum.addmaterialproperties("HT",
                                               Electrical_Materials[i])  # importing from thermal material database
                self.lum.select("materials::" + Electrical_Materials[i])
                for _ in Optical_Materials:
                    if Electrical_Materials[i].split(" ")[0] == _.split(" ")[0]:
                        Optical_Materials.remove(_)
                        self.lum.select("materials::" + Electrical_Materials[i])
                        self.lum.addmaterialproperties("EM", _)
                    else:
                        pass
            Materials_Added = Materials_Added + Electrical_Materials
            # Add all Optical Materials to materials folder
            for i in range(len(Optical_Materials)):
                self.lum.addmodelmaterial()
                self.lum.set("name", Optical_Materials[i])
                self.lum.addmaterialproperties("EM", Optical_Materials[i])
            Materials_Added = Materials_Added + Optical_Materials

        return Materials_Added

    def MZM(self, Parameters):
        '''


        Parameters
        ----------
        Parameters : dictionary
            Dictionary with all the parameters needt for the MZM creation
            Parameters['Substrate Height'] : int/float
                Substrate Height
            Parameters["Optical"] : dictionary of str
                Optical Materials Dataset
            Parameters["Electrical"] : dictionary of str
                Electrical Materials Dataset
            Parameters['angle'] : int/float
                Side angle of the Waveguife
            Parameters['Slab Height'] : Slab Height
                Height of the Material slab. It can be set to 0 if no Slab is presented
            Parameters['WG Height'] : int/float
                Waveguide Height
            Parameters['WG Width'] : int/float
                Waveguide Width. Here the Top Waveguide width is considered
            Parameters['WG Length'] : int/float
                Waveguide lenght. This determin the structure length as well
            Parameters["GND Electrodes Width"] : int/float
                Ground Electrode width
            Parameters["Signal Electrodes Width"] : int/float
                Signal electrode Width
            Parameters["Electrodes Height"] : int/float
                Height of the Metal electrodes
            Parameters["Gap"] : int/float
                Gap between the waveguide and the electrodes. The Gab is set from bottom Wg corner to electrodes.

        Returns
        -------
        None.

        '''

        # Define Materials
        Substrate_Height = Parameters['Substrate Height']
        Optical_Material = Parameters["Optical"]
        Electrical_Material = Parameters["Electrical"]
        angle = Parameters['angle']
        Slab_Height = Parameters['Slab Height']
        WG_Height = Parameters['WG Height']
        WG_Width = Parameters['WG Width']
        WG_Length = Parameters['WG Length']
        Metal_GND_Width = Parameters["GND Electrodes Width"]
        Metal_Sig_Width = Parameters["Signal Electrodes Width"]
        Metal_Height = Parameters["Electrodes Height"]
        Gap = Parameters["Gap"]

        # Add materials to Simulation enviroment
        Materials_Dict = {}
        Materials_Dict["Electrical"] = Electrical_Material
        Materials_Dict["Optical"] = Optical_Material
        self.Material(Materials_Dict)

        # Material definition
        if "Air" in Electrical_Material:
            Electrical_Material.remove("Air")
            MaterialElectrodes = Electrical_Material[0]
        else:
            MaterialElectrodes = Electrical_Material[0]

        if "Si (Silicon)" in Electrical_Material:
            MaterialSub = Electrical_Material[2]
            MaterialClad = Electrical_Material[2]
            MaterialSlab = "Si (Silicon)"
            MaterialWG = MaterialSlab
        else:
            MaterialSub = Electrical_Material[2]
            MaterialClad = Electrical_Material[2]
            MaterialSlab = Electrical_Material[1]
            MaterialWG = MaterialSlab

        # Device Lenght
        MZM_Leght = WG_Length
        MZM_Width = WG_Width * 2 + 2 * Metal_GND_Width + Metal_Sig_Width + Gap * 2

        # creating the LN Handle
        self.lum.addrect()
        self.lum.set("name", "LN_Handle")
        self.lum.set("y", 0)
        self.lum.set("y span", MZM_Width + 4e-6)
        self.lum.set("z", -Substrate_Height / 2 - (5 / 2) * 1e-6)
        self.lum.set("z span", 5e-6)
        self.lum.set("x", 0)
        self.lum.set("x span", MZM_Leght)
        self.lum.set("material", MaterialSlab)

        # creating the substrate
        self.lum.addrect()
        self.lum.set("name", "Substrate")
        self.lum.set("y", 0)
        self.lum.set("y span", MZM_Width + 4e-6)
        self.lum.set("z", 0)
        self.lum.set("z span", Substrate_Height)
        self.lum.set("x", 0)
        self.lum.set("x span", MZM_Leght)
        self.lum.set("material", MaterialSub)
        # self.lum.set("color",[1; 1; 0; 0])

        # Position Thin Film and Waveguides
        if Slab_Height == 0:
            z_Offset = Substrate_Height / 2


        else:
            # creating the thin film
            self.lum.select("Substrate")
            zmax = self.lum.get("z max")
            z_Offset = zmax + Slab_Height

            self.lum.addrect()
            self.lum.set("name", "Slab")
            self.lum.set("y", 0)
            self.lum.set("y span", MZM_Width + 4e-6)
            self.lum.set("x", 0)
            self.lum.set("x span", MZM_Leght)
            self.lum.set("z min", Substrate_Height / 2)
            self.lum.set("z max", z_Offset)
            self.lum.set("material", MaterialSlab)
            self.lum.select("Slab")
            self.lum.addtogroup("LNOI")

        # Triangle EQ for MMI Width
        x = abs(WG_Height / (np.cos((angle) * np.pi / 180)))  # in Radians
        extention = np.sqrt(x ** 2 - WG_Height ** 2)
        WG_W = WG_Width + 2 * extention
        WG_Width_top = WG_W

        # Set offsets for the Optcal Waveguides
        WG_Y_Pos = Metal_Sig_Width / 2 + Gap + WG_Width / 2
        Metal_Y_Pos = WG_Y_Pos + Gap + WG_Width / 2 + Metal_GND_Width / 2  # Metal_GND_Width + Gap*2 + WG_Width

        # Add Waveguides
        # self.lum.addwaveguide()
        # self.lum.set("name", "Waveguide_Left")
        # self.lum.set("base height", WG_Height)
        # self.lum.set("base angle", 90 - angle)
        # self.lum.set("base width", WG_Width)
        # self.lum.set("x", 0)
        # self.lum.set("y", -WG_Y_Pos)
        # self.lum.set("z", z_Offset + WG_Height/2)
        # pole = np.array([[WG_Length/2, 0], [- WG_Length/2, 0]])
        # self.lum.set("poles", pole)
        # self.lum.set("material", MaterialWG)

        self.lum.addwaveguide()
        self.lum.set("name", "Waveguide_Right")
        self.lum.set("base height", WG_Height)
        self.lum.set("base angle", 90 - angle)
        self.lum.set("base width", WG_Width)
        self.lum.set("x", 0)
        self.lum.set("y", WG_Y_Pos)
        self.lum.set("z", z_Offset + WG_Height / 2)
        pole = np.array([[WG_Length / 2, 0], [- WG_Length / 2, 0]])
        self.lum.set("poles", pole)
        self.lum.set("material", MaterialWG)

        self.lum.select("Waveguide_Right")
        self.lum.addtogroup("LNOI")

        # Add Electrodes
        print("z min ", z_Offset)
        print("z min ", z_Offset + Metal_Height)

        self.lum.addrect()
        self.lum.set("name", "Sig")
        self.lum.set("y", 0)
        self.lum.set("y span", Metal_Sig_Width)
        self.lum.set("x", 0)
        self.lum.set("x span", MZM_Leght)
        self.lum.set("z min", z_Offset)
        self.lum.set("z max", z_Offset + Metal_Height)
        self.lum.set("material", MaterialElectrodes)

        # self.lum.addrect()
        # self.lum.set("name", "Ground_L")
        # self.lum.set("y", -Metal_Y_Pos)
        # self.lum.set("y span", Metal_GND_Width)
        # self.lum.set("x", 0)
        # self.lum.set("x span", MZM_Leght)
        # self.lum.set("z min", z_Offset)
        # self.lum.set("z max", z_Offset + Metal_Height)
        # self.lum.set("material", MaterialElectrodes)

        self.lum.addrect()
        self.lum.set("name", "Ground_R")
        self.lum.set("y", Metal_Y_Pos)
        self.lum.set("y span", Metal_GND_Width)
        self.lum.set("x", 0)
        self.lum.set("x span", MZM_Leght)
        self.lum.set("z min", z_Offset)
        self.lum.set("z max", z_Offset + Metal_Height)
        self.lum.set("material", MaterialElectrodes)

    def Set_SimulationRegion(self, Parameters):
        '''


        Parameters
        ----------
        Parameters : Dictionary of Parameters
            Dictionary with all the parameters needt for the MZM simulation region. The Simulation region is define
            only on the halft of the structure sice the MZM is symetrical component.
            Parameters['Substrate Height'] : int/float
                Substrate Heigh
            Parameters['Slab Height'] : Slab Height
                Height of the Material slab. It can be set to 0 if no Slab is presented
            Parameters['WG Height'] : int/float
                Waveguide Height
            Parameters['WG Width'] : int/float
                Waveguide Width. Here the Top Waveguide width is considered
            Parameters["GND Electrodes Width"] : int/float
                Ground Electrode width
            Parameters["Signal Electrodes Width"] : int/float
                Signal electrode Width
            Parameters["Electrodes Height"] : int/float
                Height of the Metal electrodes
            Parameters["Gap"] : int/float
                Gap between the waveguide and the electrodes. The Gab is set from bottom Wg corner to electrodes.


        Returns
        -------
        None.

        '''

        # Define Parameters
        Substrate_Height = Parameters['Substrate Height']
        Slab_Height = Parameters['Slab Height']
        WG_Height = Parameters['WG Height']
        WG_Width = Parameters['WG Width']
        Metal_GND_Width = Parameters["GND Electrodes Width"]
        Metal_Sig_Width = Parameters["Signal Electrodes Width"]
        Metal_Height = Parameters["Electrodes Height"]
        Gap = Parameters["Gap"]

        Metal_Y_Pos = Metal_GND_Width + Gap * 2 + WG_Width
        MZM_Width = WG_Width * 2 + 2 * Metal_GND_Width + Metal_Sig_Width + Gap * 4

        self.lum.select("simulation region")
        self.lum.set("dimension", "2D X-Normal")
        # self.lum.addmodelmaterial()
        # self.lum.set("name", "Air")
        # self.lum.addmaterialproperties("CT", "Air")  # importing from electrical material database
        # self.lum.addmaterialproperties("HT", "Air")  # importing from thermal material database)
        self.lum.set("background material", "Air")
        self.lum.set("x", 0)
        self.lum.set("y min", 0 - Metal_Sig_Width / 2)
        self.lum.set("y max", Metal_Y_Pos + Metal_GND_Width / 2)
        # self.lum.set("y",0)
        # self.lum.set("y span", MZM_Width + 1e-6)
        self.lum.set("z min", 0)
        self.lum.set("z max", Substrate_Height + Slab_Height + 2 * Metal_Height)

    def MZM_ChargeSolver(self, Parameters):

        # Parameters
        Substrate_Height = Parameters['Substrate Height']
        Slab_Height = Parameters['Slab Height']
        WG_Height = Parameters['WG Height']
        WG_Width = Parameters['WG Width']
        Metal_GND_Width = Parameters["GND Electrodes Width"]
        Metal_Sig_Width = Parameters["Signal Electrodes Width"]
        Metal_Height = Parameters["Electrodes Height"]
        Gap = Parameters["Gap"]

        MZM_Width = WG_Width * 2 + 2 * Metal_GND_Width + Metal_Sig_Width + Gap * 4
        Metal_Y_Pos = Metal_GND_Width + Gap * 2 + WG_Width

        # Set simulation region
        self.Set_SimulationRegion(Parameters)

        # Set charge Solver
        self.lum.addchargesolver()
        self.lum.set("min edge length", 0.04e-6)  # TODO
        # TODO Results Tab syntax

        # Set boundary
        self.lum.select("CHARGE::boundary conditions")
        # self.lum.addelectricalcontact()
        # self.lum.set("name", "Ground_L")
        # self.lum.set("surface type", "domain")
        # self.lum.set("domain", 1)
        # self.lum.set("name", "Ground_L")
        # self.lum.set("sweep type", "single")
        # self.lum.set("voltage", 0)
        # self.lum.set("surface type", "solid")
        # self.lum.set("solid", "Ground_R")
        # self.lum.set("outer surface only", 0)

        self.lum.addelectricalcontact()
        self.lum.set("name", "Signal")
        self.lum.set("sweep type", "range")
        self.lum.set("range start", 0)
        self.lum.set("range stop", 5)
        self.lum.set("range interval", 0.5)
        self.lum.set("range num points", 11)
        self.lum.set("surface type", "solid")
        self.lum.set("solid", "Sig")
        self.lum.set("outer surface only", 1)

        self.lum.addelectricalcontact()
        self.lum.set("name", "Ground_R")
        self.lum.set("sweep type", "single")
        self.lum.set("voltage", 0)
        self.lum.set("surface type", "solid")
        self.lum.set("solid", "Ground_R")
        self.lum.set("outer surface only", 0)

        # Set moditor
        self.lum.addefieldmonitor()
        self.lum.set("name", "CHARGE_Field_Monitor")
        self.lum.set("monitor type", "2D x-normal")
        self.lum.set("x", 0)
        # self.lum.set("y",0)
        # self.lum.set("y span", MZM_Width + 1e-6)
        self.lum.set("y min", 0 - Metal_Sig_Width / 2)
        self.lum.set("y max", Metal_Y_Pos + Metal_GND_Width / 2)
        self.lum.set("z min", -Substrate_Height / 2)
        self.lum.set("z max", Substrate_Height + Slab_Height + Metal_Height / 2)

        # Add fine mesh
        self.lum.addmesh()
        self.lum.set("name", "CHARGE_Volume_Mesh")
        self.lum.set("geometry type", "volume")
        self.lum.set("volume type", "solid")
        self.lum.set("volume solid", "Waveguide_Right")
        self.lum.set("max edge length", 0.01e-6)  # TODO need to be checked and made variable

    def MZM_FEEMSolver(self, Parameters):
        '''


        Parameters
        ----------
        Parameters : dict
            Dictionary with all the parameters needt for the MZM FEEM Solver region. The Simulation region is define
            only on the halft of the structure since the FEEM region is the simulation region too.
            Parameters['Substrate Height'] : int/float
                Substrate Heigh
            Parameters['Slab Height'] : Slab Height
                Height of the Material slab. It can be set to 0 if no Slab is presented
            Parameters['WG Height'] : int/float
                Waveguide Height
            Parameters['WG Width'] : int/float
                Waveguide Width. Here the Top Waveguide width is considered
            Parameters["GND Electrodes Width"] : int/float
                Ground Electrode width
            Parameters["Signal Electrodes Width"] : int/float
                Signal electrode Width
            Parameters["Electrodes Height"] : int/float
                Height of the Metal electrodes
            Parameters["Gap"] : int/float
                Gap between the waveguide and the electrodes. The Gab is set from bottom Wg corner to electrodes.
            Parameters["Wavelength"] : int/float
                Simulation wavelength

        Returns
        -------
        None.

        '''

        # Parameters
        Substrate_Height = Parameters['Substrate Height']
        Slab_Height = Parameters['Slab Height']
        WG_Height = Parameters['WG Height']
        WG_Width = Parameters['WG Width']
        Metal_GND_Width = Parameters["GND Electrodes Width"]
        Metal_Sig_Width = Parameters["Signal Electrodes Width"]
        Metal_Height = Parameters["Electrodes Height"]
        Gap = Parameters["Gap"]
        Wavelength = Parameters["Wavelength"]

        MZM_Width = WG_Width * 2 + 2 * Metal_GND_Width + Metal_Sig_Width + Gap * 4

        # Set Solver
        self.lum.addfeemsolver()
        self.lum.set("edges per wavelength", 2)  # TODO Check for Converchange need to be separate parameter
        self.lum.set("polynomial order", 2)  # TODO Check what is this parameter for
        self.lum.set("wavelength", Wavelength)
        self.lum.set("use max index", 0)
        self.lum.set("n", 2.02)  # TODO Check for what exactly supposable under limit of n_eff

        # Set Boundary
        self.lum.addpec()
        self.lum.set("surface type", "simulation region")
        self.lum.set("x min", 1)
        self.lum.set("x max", 1)
        self.lum.set("y min", 1)
        self.lum.set("y max", 1)
        self.lum.set("z min", 1)
        self.lum.set("z max", 1)

        self.lum.addpml()
        self.lum.set("sigma", 5)  # TODO Check for what exactly

        # Set monitor
        self.lum.addimportnk()
        self.lum.set("name", "nk import WG")
        self.lum.set("volume type", "solid")
        self.lum.set("volume solid", "Waveguide_Right")

        self.lum.addimportnk()
        self.lum.set("name", "nk import Slab")
        self.lum.set("volume type", "solid")
        self.lum.set("volume solid", "Slab")

    def StartCHARGESolver(self):
        '''
        This function will first save the simulation into the folder where the
        this script is saved. After saving the script the CHARGE solver will be started.

        Returns
        -------
        None.


        '''
        self.lum.save('SimRun1')
        self.lum.run("CHARGE")

    def StartFEEMSolver(self):
        '''
        This function will first save the simulation into the folder where the
        this script is saved. After saving the script the FEEM solver will be started.

        Returns
        -------
        None.

        '''
        self.lum.save("SimRun1")
        self.lum.run("FEEM")

    def ResultsCHARGE(self):
        '''

        Returns
        -------
        electro : dict
         Simulation results from the monitor of CHARGE. This function will create an Lumerical window into Charge solver.
         The user can use the window to see tha changes in EO index depending on the voltage apply on the electrodes.
         The function will return the monitor Data, so that the user can use it for some post precessing.

        '''

        # Lithium Niobate telecom permitivity
        eps_o = 2.21 ** 2
        eps_e = 2.14 ** 2

        # Lithium Niobate non;linear coefficents
        r_13 = 9.6e-12
        r_33 = 30.9e-12

        electro = self.lum.getresult("CHARGE::CHARGE_Field_Monitor", "electrostatics")

        # Get electrostatic results
        E = np.squeeze(electro["E"])
        Volt = electro["V_Signal"]

        # Intialize perturbation matrices same size as E
        dts = E.shape
        n_EO = np.zeros(dts)
        dn = np.zeros(dts)

        # Vector components work well with unstructured datasets, need to do this for voltage and wavelength
        for vv in range(len(Volt)):
            # Spatial index data (diagonal permittivity)
            eps_unperturbed = np.array([eps_e, eps_o, eps_o])[np.newaxis, :]  # Shape: (1, 3)

            # Pockels effect (perturbation term)
            deps_inv = np.array([r_33 * E[:, vv, 0], r_13 * E[:, vv, 0], r_13 * E[:, vv, 0]]).T  # Shape: (dts[0], 3)

            # Compute modified index using element-wise inversion (instead of np.linalg.inv)
            n_EO[:, vv, :] = np.sqrt(1 / (1 / eps_unperturbed + deps_inv))

        dn = np.copy(n_EO)
        dn[:, :, 0] = n_EO[:, :, 0] - np.sqrt(eps_e)
        dn[:, :, 1:3] = n_EO[:, :, 1:3] - np.sqrt(eps_o)

        # Add dn and n_EO to dataset and visualize
        electro['n_EO'] = n_EO  ## total index
        electro['dn'] = dn  ## total index
        self.lum.visualize(electro)

        return electro
        
    


        
        
        
        

class HelpSubject:
    
    @classmethod
    def Help_Objects(self):
        print("""
        # ============================================================================= #
        #                Welcome to the Help Menu for Lumerical Structure               #     
        #                You can ask for help for one of the following structures       #   
        #                                                                               #
        #                1) Waveguide - used for FDE simulation only!                   #  
        #                2) MMI2x1 - can be used in EME and FDTD solvers                #  
        #                3) MMI2x2 - can be used in EME and FDTD solvers                #
        #                4) InverseTaper - can be used in EME and FDTD solvers          #  
        #                5) Directional Coupler - can be used in EME and FDTD solvers   #  
        #                6) StraightWaveguide - can be used in EME and FDTD solvers     #     
        #                7) WDM (Wavelength Division Multiplexing, angeled MMI) -onyl   #
        #                   available for EME Solver                                    #  
        #                8) BendWaveguide - onyl available for FDTD Solver              # 
        #                9) ArcWaveguide - onyl available for FDTD Solver               #        
        #                10) GratingCoupler - onyl available for FDTD Solver            #         
        #                11) RingGratingCoupler - onyl available for FDTD Solver        #  
        #                                                                               #  
        #                To print information about the structure you desire please     #  
        #                use obj.Help({"Objects": Number}). For Example                 #  
        #                obj.Help({"Objects": 1}) will give you information about       #  
        #                Structure Waveguide                                            #

        # ============================================================================= #
        """)
        
    def Structures(self, Number):
        if Number == 1:
            print("""
                    Call Waveguide with -> obj.Waveguide(Parameters)
                        
                    Dictionary Parameters:
-----------------------------------------------------------------------------------------------------------
                    Parameters['Substrate Height'] : int7float
                        Height of the substrate.
                    Parameters['WG Length']: int/float
                        Length of the Waveguide.
                    Parameters['WG Height'] : int/float
                        Height of the Waveguide
                    Parameters['WG Width'] : int/float
                        Width of the Waveguide
                    Parameters['angle'] : int
                        Angle of the Waveguide Walls. it is calculated WG_angle = 90 - angle.
                        For anfle = 90 we get a perfect rect!
                    Parameters['Slab Height'] : int/float
                        Slab height.
                    Parameters['Material'] : list of str
                        List of Materials. The list should be with names (str) of a valid Lumerical materials.
                        Check the names in Lumerical Materials viewer.
                        The List of materials must contain at least 2 materials! 
                        Parameters['Material'] = ['Cladding/Substrat', 'Object Material'].
                        For Example: Parameters['Material'] = ["SiO2 (Glass) - Palik", 'LiNbO3_20deg_X cut'].
                        
-----------------------------------------------------------------------------------------------------------
                """)
        elif Number == 6:
            print("""
                    Call StraightWaveguide with -> obj.StraightWaveguide(Parameters)
                        
                    Dictionary Parameters:
-----------------------------------------------------------------------------------------------------------
                    Parameters['Substrate Height'] : int float
                        Substrate Height
                    Parameters['WG Height'] : int/float
                        Height of the Waveguide
                    Parameters['WG Width'] : int/float
                        Width of the Waveguide
                    Parameters['WG Length'] : int/float
                        Length of the Waveguide
                    Parameters['Taper'] : boolen
                        If Taper == False, only Waveguiedes will be constructed.
                        If Taper == True, only an single Taper will be constructed
                    Parameters['Taper Length'] : int/float
                        Taper Length
                    Parameters['Taper Width'] : int/float
                        Taper backside width, the frontside width is the waveguide width
                    Parameters['angle'] : int
                        Angle of the Waveguide Walls. it is calculated WG_angle = 90 - angle.
                        For anfle = 90 we get a perfect rect!
                    Parameters['Slab Height'] : int/float
                        Slab height.
                    Parameters['Material'] : list
                        List of Materials. The list should be with names (str) of a valid Lumerical materials.
                        Check the names in Lumerical Materials viewer.
                        The List of materials must contain at least 2 materials! 
                        Parameters['Material'] = ['Cladding/Substrat', 'Object Material'].
                        For Example: Parameters['Material'] = ["SiO2 (Glass) - Palik", 'LiNbO3_20deg_X cut'].
                    Parameters["Wavelength"] : int/float
                        Wavelength
                    Parameters["Waveguide Angle"] : int/float
                        Bending angle of the Waveguide. Set it to 0 to simulate the straight waveguide.
                        If Waveguide Angle is different then 0, then the straight waveguide will be tilted
                        at the choosen degrees.
                    Parameters["Cladding"] : anything, optional
                        This function will check if you have set Parameters["Cladding"] to anaything, for example "Parameters["Cladding"]=1" 
                        and if so it will put cladding over your structure. If the user didnt give the "Cladding" as dictionary key no cladding 
                        will be set.
                    Parameters["Taper Type"] : anything, optional
                        This function will check if you have set Parameters["Taper Type"] to anaything, for example "Parameters["Taper Type"]=1" 
                        and if so it will design an Inverse Taper Structure with no Cladding. Here the option "Cladding" is not active and will be ignored.
                        If the user didnt give the "Taper Type" as dictionary key, then an normal taper structure will be simulated.
                        If Taper Type is selected the user need to provide additional information:
                            TaperWidthF = Parameters['PWB Taper Width Front']
                            TaperWidthB = Parameters['PWB Taper Width Back']
                            TaperHightB = Parameters['PWB Taper Hight Back']
                            TaperHightF = Parameters['PWB Taper Hight Front']
                            TaperLength_PWB = Parameters['PWB Taper Length']
                        
-----------------------------------------------------------------------------------------------------------
                """)
        elif Number == 8:
            print("""
                    Call BendWaveguide with -> obj.BendWaveguide(Parameters) 
                        
                    Dictionary Parameters:
-----------------------------------------------------------------------------------------------------------
                    Parameters['Substrate Height'] : int/float
                        Substrate Height
                    Parameters['WG Height'] : int/float
                        Height of the Waveguide (Etching depth)
                    Parameters['WG Width'] : int/float
                        Top width of the Waveguide.
                    Parameters['angle'] : int/float
                        Angle of the Waveguide Walls. It is calculated WG_angle = 90 - angle.
                        For angle = 90 , a perfect rect is created!
                    Parameters['Slab Height'] : int/float
                        Slab height
                    Parameters['Material'] : list of str
                        List of Materials. The list should be with names (str) of a valid Lumerical materials.
                        Check the names in Lumerical Materials viewer.
                        The List of materials must contain at least 2 materials! 
                        Parameters['Material'] = ['Cladding/Substrat', 'Object Material'].
                        For Example: Parameters['Material'] = ["SiO2 (Glass) - Palik", 'LiNbO3_20deg_X cut'].
                    Parameters["Wavelength"] : int/float
                        Wavelength
                    Parameters["x span"] : int/float
                        Length of the S Curve. Span of the object in x direction.
                    Parameters["y span"] : int/float
                        Height of the curve. Span (difference between the input and output ports) of the S-curve in y direction.
                    Parameters["poles"] : boolen
                        If Parameters["poles"] = True an Bezier Curbe will be made.
                        if Parameters["poles"] = False an Cosinus Curve = y_span*(cos((pi/(2*x_span))*t)^2) will be made. Where
                        t is in the range of 0 - y_span
 -----------------------------------------------------------------------------------------------------------
                """)                       
        elif Number == 9:
            print("""
                    Call ArcWaveguide with -> obj.ArcWaveguide(Parameters) 
                        
                    Dictionary Parameters:
-----------------------------------------------------------------------------------------------------------
                    Parameters['Substrate Height'] : int/float
                        Substrate height
                    Parameters['WG Height'] : int/float
                        Waveguide height
                    Parameters['WG Width'] : int/float
                        Waveguide width
                    Parameters['angle'] : int/float
                        Angle of the Waveguide Walls. it is calculated WG_angle = 90 - angle.
                        For anfle = 90 we get a perfect rect!
                    Parameters['Slab Height'] : int/float
                        Slab height
                    Parameters['Material'] : list
                        List of Materials. The list should be with names (str) of a valid Lumerical materials.
                        Check the names in Lumerical Materials viewer.
                        The List of materials must contain at least 2 materials! 
                        Parameters['Material'] = ['Cladding/Substrat', 'Object Material'].
                        For Example: Parameters['Material'] = ["SiO2 (Glass) - Palik", 'LiNbO3_20deg_X cut'].
                    Parameters["Wavelength"] : int/float
                        Wavelength
                    Parameters["S_Band Radius"] : int/float
                        Radius of the Circle in um
                    Parameters['Arc deg'] : int
                        Can be 90 or 180 for 1/4 of a cirle or 1/2 of a circle.
                    Parameters["Cladding"] : anything, optional
                        This function will check if you have set Parameters["Cladding"] to anaything, for example "Parameters["Cladding"]=1" 
                        and if so it will put cladding over your structure. If the user didnt give the "Cladding" as dictionary key no cladding 
                        will be set.
 -----------------------------------------------------------------------------------------------------------
                """)                         
        elif Number == 2:
            print("""
                    Call MMI2x1 with -> obj.MMI2x1(Parameters)    
                        
                    Dictionary Parameters:
-----------------------------------------------------------------------------------------------------------
                     Parameters['Material'] : list of str
                        List of Materials. The list should be with names (str) of a valid Lumerical materials.
                        Check the names in Lumerical Materials viewer.
                        The List of materials must contain at least 2 materials! 
                        Parameters['Material'] = ['Cladding/Substrat', 'Object Material'].
                        For Example: Parameters['Material'] = ["SiO2 (Glass) - Palik", 'LiNbO3_20deg_X cut'].
                    Parameters['Substrate Height'] : int/float
                        Substrate height.
                    Parameters['MMI Width'] : int/float
                        Width of the MMI.
                    Parameters['MMI Length'] : int/float
                        Length of the MMI.
                    Parameters['angle'] : int/float
                        Angle of the Waveguide Walls. it is calculated WG_angle = 90 - angle.
                        For anfle = 90 we get a perfect rect!
                    Parameters['WG Height'] : int/float
                        Waveguide hight. Also the height of the MMI section
                     Parameters['WG Width'] : int/float
                        Waveguide width.
                    Parameters['WG Length'] : int/float
                        Waveguide length.
                    Parameters['Position Offset'] : int/float
                        Offset between the waveguides. If Taper == True then this become the offset
                        betweent he tapers wider sides. Waveguide and Tapers cannot be placed ourside
                        the MMI structure. The minimum distance between Taper and Waveguide is 1 um
                        becouse of manufactering restrictions in the University.
                    Parameters['Offset Input'] : int/float
                        Offset of the input waveguide.
                    Parameters['Slab Height'] : int/float
                        Height of the slab.
                    Parameters['Taper'] : boolen
                        Taper can be set to be True ot False.
                        if Taper == False - No Taper used
                        if Taper == True - Taper placed
                    Parameters['Taper Length'] : int/float
                        If Taper == True, then this will set the Tapers length. If Taper == False
                        this will be ignored and some random value can be given.
                    Parameters['Taper Width'] : int/float
                        If Taper == True, then this will set the Tapers width. If Taper == False
                        this will be ignored and some random value can be given.
                    Parameters["Cladding"] : anything, optional
                        This function will check if you have set Parameters["Cladding"] to anaything, for example "Parameters["Cladding"]=1" 
                        and if so it will put cladding over your structure. If the user didnt give the "Cladding" as dictionary key no cladding 
                        will be set.
                    Parameters["Offset Output"] : anything, optional
                        This function will allow the user to move the outputs in oposite direction. Please dont use it since is there only 
                        becouse the maschine of our physic departmant had some proiblems with the LNOI objects design. 
 -----------------------------------------------------------------------------------------------------------
                """)                          
        elif Number == 3:
            print("""
                    Call MMI2x2 with -> obj.MMI2x2(Parameters)   
                        
                    Dictionary Parameters:
-----------------------------------------------------------------------------------------------------------
                    Parameters['Material'] : list of str
                        List of Materials. The list should be with names (str) of a valid Lumerical materials.
                        Check the names in Lumerical Materials viewer.
                        The List of materials must contain at least 2 materials! 
                        Parameters['Material'] = ['Cladding/Substrat', 'Object Material'].
                        For Example: Parameters['Material'] = ["SiO2 (Glass) - Palik", 'LiNbO3_20deg_X cut'].
                     Parameters['Substrate Height'] : int/float
                        Substrate height.
                    Parameters['MMI Width'] : int/float
                        Width of the MMI.
                    Parameters['MMI Length'] : int/float
                        Length of the MMI.
                    Parameters['angle'] : int/float
                        Angle of the Waveguide Walls. it is calculated WG_angle = 90 - angle.
                        For anfle = 90 we get a perfect rect!
                    Parameters['WG Height'] : int/float
                        Waveguide hight. Also the height of the MMI section
                    Parameters['WG Width'] : int/float
                        Waveguide width.
                    Parameters['WG Length'] : int/float
                        Waveguide length.
                    Parameters['Position Offset'] : int/float
                        Offset between the waveguides. If Taper == True then this become the offset
                        betweent he tapers wider sides. Waveguide and Tapers cannot be placed ourside
                        the MMI structure. The minimum distance between Taper and Waveguide is 1 um
                        becouse of manufactering restrictions in the University.
                     Parameters['Slab Height'] : int/float
                        Height of the slab.
                    Parameters['Taper'] : boolen
                        Taper can be set to be True ot False.
                        if Taper == False - No Taper used
                        if Taper == True - Taper placed
                    Parameters['Taper Length'] : int/float
                        If Taper == True, then this will set the Tapers length. If Taper == False
                        this will be ignored and some random value can be given.
                    Parameters['Taper Width'] : int/float
                        If Taper == True, then this will set the Tapers width. If Taper == False
                        this will be ignored and some random value can be given.
                    Parameters["Cladding"] : anything, optional
                            This function will check if you have set Parameters["Cladding"] to anaything, for example "Parameters["Cladding"]=1" 
                            and if so it will put cladding over your structure. If the user didnt give the "Cladding" as dictionary key no cladding 
                            will be set.
 -----------------------------------------------------------------------------------------------------------
                """)                         
        elif Number == 7:
            print("""
                    Call WDM with -> obj.WDM(Parameters) 
                        
                    Dictionary Parameters:
-----------------------------------------------------------------------------------------------------------                      
                    Parameters['Material'] : list of str
                        List of Materials. The list should be with names (str) of a valid Lumerical materials.
                        Check the names in Lumerical Materials viewer.
                        The List of materials must contain at least 2 materials! 
                        Parameters['Material'] = ['Cladding/Substrat', 'Object Material'].
                        For Example: Parameters['Material'] = ["SiO2 (Glass) - Palik", 'LiNbO3_20deg_X cut'].
                    Parameters['Substrate Height'] : int/float
                        Substrate height.
                    Parameters['MMI Width'] : int/float
                        Width of the MMI.
                    Parameters['MMI Length'] : int/float
                        Length of the MMI.
                    Parameters['angle'] : int/float
                        Angle of the Waveguide Walls. it is calculated WG_angle = 90 - angle.
                        For anfle = 90 we get a perfect rect!
                    Parameters['WG Height'] : int/float
                        Waveguide hight. Also the height of the MMI section
                    Parameters['WG Width'] : int/float
                        Waveguide width.
                    Parameters['Wavelength'] : int/float
                        Wavelength.
                    Parameters['WG Length'] : int/float
                        Waveguide length.
                    Parameters['Slab Height'] : int/float
                        Height of the slab.
                    Parameters['Angle Thetha'] :
                        Input and output angle of the waveguide. This is only temporally
                    Parameters['Taper Width'] : int/float
                        Backside width of the Taper, frontside width is the waveguide width
                    Parameters['Taper Length'] : int/float
                        Length of the Taper.
                    Parameters['Taper'] : boolen
                        If Taper == False, no Taper will be placed with the Waveguides.
                        If Taper == True, tapers will be placed with the waveguides
                    Parameters["Cladding"] : anything, optional
                        This function will check if you have set Parameters["Cladding"] to anaything, for example "Parameters["Cladding"]=1" 
                        and if so it will put cladding over your structure. If the user didnt give the "Cladding" as dictionary key no cladding 
                        will be set.
 -----------------------------------------------------------------------------------------------------------
                """)   
        elif Number == 4:
            print("""
                    Call InverseTaper with -> obj.InverseTaper(Parameters) 
                        
                    Dictionary Parameters:
-----------------------------------------------------------------------------------------------------------  
                    Parameters['Material'] : list of str
                        List of Materials. The list should be with names (str) of a valid Lumerical materials.
                        Check the names in Lumerical Materials viewer.
                        The List of materials must contain at least 2 materials! 
                        For this structure 3 materials will be needed!
                        Parameters['Material'] = ['Cladding/Substrat', 'Object Material'].
                        For Example: Parameters['Material'] = ["SiO2 (Glass) - Palik", 'LiNbO3_20deg_X cut'].
                    Parameters['Substrate Height'] : int/float
                        Substrate height.
                    Parameters['angle'] : int/float
                        Angle of the Waveguide Walls. it is calculated WG_angle = 90 - angle.
                        For anfle = 90 we get a perfect rect!
                    Parameters['WG Height'] : int/float
                        Waveguide hight. Also the height of the MMI section
                    Parameters['WG Width'] : int/float
                        Waveguide width. Also in this function and ONLY in this function this will be the
                        ibverse Taper width!!!
                    Parameters['Slab Height'] : int/float
                        Slab height
                    Parameters['Taper Length'] : int/float
                        Inverse Taper Length
                    Parameters['Taper Width'] : int/float
                        Front Width of the inverse Taper!! In this Function and ONLY in this function, the
                        Parameters['Taper Width'] is the frond width of the inverse Taper!
                    Parameters['PWB Taper Width Back'] : int/float
                        Photonic Wirebonding (PWB) Width back side (to the Photonic Wire Bonding)
                    Parameters['PWB Taper Hight Back'] : int/float
                        Photonic Wire Bonding Height back side (to the Photonic Wire Bonding)
                    Parameters['PWB Taper Width Front'] : int/float
                        Photonic Wirebonding (PWB) Width front side (to the photonic waveguide)
                    Parameters['PWB Taper Hight Front'] : int/float
                        Photonic Wire Bonding Height front side (to the photonic waveguide)
                    Parameters['PWB Taper Length'] : int/float
                        Length of the Photonic Wire Bonding Taper
                    Parameters["SMF Core Diameter"] : int/float
                        Single Mode Fiber core Diameter
                    Parameters["SMF Cladding Diameter"] : int/float
                        Single Mode Fiber Cladding Diameter
                    Parameters["SMF Core Index"]
                        Single Mode Fiber Core Index
                    Parameters["SMF Cladding Index"]
                        Single Mode Fiber Cladding Index                  
 -----------------------------------------------------------------------------------------------------------
                """)  
        elif Number == 5:
            print("""
                    Call DirectionalCoupler with -> obj.DirectionalCoupler(Parameters) 
                        
                    Dictionary Parameters:
----------------------------------------------------------------------------------------------------------- 
                    Parameters['Material'] : list of str
                        List of Materials. The list should be with names (str) of a valid Lumerical materials.
                        Check the names in Lumerical Materials viewer.
                        The List of materials must contain at least 3 materials! 
                        Parameters['Material'] = ['Cladding/Substrat', 'Object Material', 'Photonic Wire Bonding'].
                        For Example: Parameters['Material'] = ["SiO2 (Glass) - Palik", 'LiNbO3_20deg_X cut', 'Si (Silicon) - Palik'].
                    Parameters['Substrate Height'] : int/float
                        Substrate height.
                    Parameters['Substrate Width'] : int/float
                        Substrate Width.
                    Parameters['DC Length'] : int/float
                        Length of the directional coupler
                    Parameters['angle'] : int/float
                        Angle of the Waveguide Walls. it is calculated WG_angle = 90 - angle.
                        For anfle = 90 we get a perfect rect!
                    Parameters['WG Height'] : int/float
                        Waveguide hight. Also the height of the MMI section
                    Parameters['WG Width'] : int/float
                        Waveguide width.
                    Parameters['Position Offset'] : int/float
                        Offset between the waveguides. The minimum distance between Waveguides is 1 um
                        becouse of manufactering restrictions in the University.
                    Parameters['Slab Height'] : int/float
                        Height of the slab.
                    Parameters['Wavelength'] : int/float
                        Wavelength
                    Parameters["Cladding"] : anything, optional
                        This function will check if you have set Parameters["Cladding"] to anaything, for example "Parameters["Cladding"]=1" 
                        and if so it will put cladding over your structure. If the user didnt give the "Cladding" as dictionary key no cladding 
                        will be set.
 -----------------------------------------------------------------------------------------------------------
                """) 
        elif Number == 10:
            print("""
                    Call GratingCoupler with -> obj.GratingCoupler(Parameters) 
                        
                    Dictionary Parameters:
----------------------------------------------------------------------------------------------------------- 
                    Parameters['Material GC'] : list of str
                        List of Materials. The list should be with names (str) of a valid Lumerical materials.
                        Check the names in Lumerical Materials viewer.
                        The List of materials must contain at least 3 materials! 
                        Parameters['Material'] = ['Cladding/Substrat', 'Object Material', 'Photonic Wire Bonding'].
                        For Example: Parameters['Material'] = ["SiO2 (Glass) - Palik", 'LiNbO3_20deg_X cut', 'Si (Silicon) - Palik'].
                    Parameters['Substrate Height'] : int/float
                        Substrate height.
                    Parameters["Length GC"]: int/float
                        Lenght of the Grating Coupler Area
                    Parameters["Width GC"]: int/float
                        Widht of the Grating Coupler Area
                    Parameters["Hight GC"]: int/float
                        Hight of the Grating Coupler Material
                    Parameters["Etch Depth GC"]: int/float
                        How deep, taken from the Parameters["Hight GC"] will be the etchin depth of the gratings
                    Parameters["Duty Cycle"]: int/float
                        Duty cycle of the gratings. For example Parameters["Duty Cycle"] = 0.39 will result in 39% Duty Cycle 
                    Parameters["Pitch GC"]: int/float
                        Pitch of the Grating Coupler. For Example Parameters["Pitch GC"] = 0.6e-6 will result in 6um Etch Space + Rib Space = 0.6um.
                    Parameters["Input LengthGC"]: int/float
                        An squere Waveguide with the same WG Height as the Grating coupler place before the Grating Coupler region will start. 
                    Parameters["Output Length GC"]: int/float
                        An squere Waveguide with the same WG Height as the Grating coupler place after the Grating Coupler region to finish the structure.
                    Parameters['angle'] : int/float
                        Angle of the Waveguide Walls. it is calculated WG_angle = 90 - angle.
                        For anfle = 90 we get a perfect rect! In this case is used only when Parameters["Taper"] = True
                    Parameters["Taper"] : boolen
                        You can create an input Taper to your Grating Coupler structure
                    Parameters['WG Height'] : int/float
                        Waveguide hight. Also the height of the MMI section
                    Parameters['WG Width'] : int/float
                        Waveguide width. Also in this function and ONLY in this function this will be the
                        ibverse Taper width!!!
                    Parameters['Slab Height'] : int/float
                        Slab height
                    Parameters['Taper Length'] : int/float
                        Taper Length
                    Parameters["SMF Core Diameter"] : int/float
                        Single Mode Fiber core Diameter
                    Parameters["SMF Cladding Diameter"] : int/float
                        Single Mode Fiber Cladding Diameter
                    Parameters["SMF Core Index"]
                        Single Mode Fiber Core Index
                    Parameters["SMF Cladding Index"]
                        Single Mode Fiber Cladding Index
                    Parameters["SMF Theta"]: int/float
                        Tilting Angle of the Single Mode Fiber to the Grating Coupler. Normaly we choose Parameters["SMF Theta"] = 15
                    Parameters["SMF Z Span"]: int/float
                        Lenght/Span of the Single Mode Fiber
 -----------------------------------------------------------------------------------------------------------
                """)  
        elif Number == 11:
            print("""
                    Call RingGratingCoupler with -> obj.RingGratingCoupler(Parameters) 
                        
                    Dictionary Parameters:
----------------------------------------------------------------------------------------------------------- 
                    Parameters['Material GC'] : list of str
                        List of Materials. The list should be with names (str) of a valid Lumerical materials.
                        Check the names in Lumerical Materials viewer.
                        The List of materials must contain at least 3 materials! 
                        Parameters['Material'] = ['Cladding/Substrat', 'Object Material', 'Photonic Wire Bonding'].
                        For Example: Parameters['Material'] = ["SiO2 (Glass) - Palik", 'LiNbO3_20deg_X cut', 'Si (Silicon) - Palik'].
                    Parameters['Substrate Height'] : int/float
                        Substrate height.
                    Parameters["Length GC"]: int/float
                        Lenght of the Grating Coupler Area
                    Parameters["Width GC"]: int/float
                        Widht of the Grating Coupler Area
                    Parameters["Hight GC"]: int/float
                        Hight of the Grating Coupler Material
                    Parameters["GC Radius"]: int/float
                        Radius of the Ring Grating Coupler in um. For Example "Parameters["GC Radius"] = 25e-6"
                    Parameters["Etch Depth GC"]: int/float
                        How deep, taken from the Parameters["Hight GC"] will be the etchin depth of the gratings
                    Parameters["Duty Cycle"]: int/float
                        Duty cycle of the gratings. For example Parameters["Duty Cycle"] = 0.39 will result in 39% Duty Cycle 
                    Parameters["Pitch GC"]: int/float
                        Pitch of the Grating Coupler. For Example Parameters["Pitch GC"] = 0.6e-6 will result in 6um Etch Space + Rib Space = 0.6um.
                    Parameters["Input LengthGC"]: int/float
                        An squere Waveguide with the same WG Height as the Grating coupler place before the Grating Coupler region will start. 
                    Parameters["Output Length GC"]: int/float
                        An squere Waveguide with the same WG Height as the Grating coupler place after the Grating Coupler region to finish the structure.
                    Parameters['angle'] : int/float
                        Angle of the Waveguide Walls. it is calculated WG_angle = 90 - angle.
                        For anfle = 90 we get a perfect rect! 
                    Parameters['WG Height'] : int/float
                        Waveguide hight. Also the height of the MMI section
                    Parameters['WG Width'] : int/float
                        Waveguide width. 
                    Parameters['Slab Height'] : int/float
                        Slab height
                    Parameters["SMF Core Diameter"] : int/float
                        Single Mode Fiber core Diameter
                    Parameters["SMF Cladding Diameter"] : int/float
                        Single Mode Fiber Cladding Diameter
                    Parameters["SMF Core Index"]
                        Single Mode Fiber Core Index
                    Parameters["SMF Cladding Index"]
                        Single Mode Fiber Cladding Index
                    Parameters["SMF Theta"]: int/float
                        Tilting Angle of the Single Mode Fiber to the Grating Coupler. Normaly we choose Parameters["SMF Theta"] = 15
                    Parameters["SMF Z Span"]: int/float
                        Lenght/Span of the Single Mode Fiber
 -----------------------------------------------------------------------------------------------------------
                """)                 
                  
                  
                  
                  
    def Help_Solvers(self):
        print("""
        # ============================================================================= #
        #                Welcome to the Help Menu for Lumerical Structure Solver. You   #     
        #                can ask for help for one of the following structural solvers   #   
        #                                                                               #
        #                1) Waveguide FDE Solver                                        #
        #                2) MMI2x1 EME Solver                                           # 
        #                3) MMI2x2 EME Solver                                           #
        #                4) InverseTaper EME Solver                                     #      
        #                5) Directional Coupler EME Solver                              # 
        #                6) StraightWaveguide EME Solver                                #       
        #                7) WDM (Wavelength Division Multiplexing, angeled MMI) -EME    #
        #                   Solver                                                      #          
        #                8) MMI2x1 FDTD Solver                                          # 
        #                9) MMI2x2 FDTD Solver                                          #
        #                10) InverseTaper FDTD Solver                                   #      
        #                11) Directional Coupler FDTD Solver                            # 
        #                12) StraightWaveguide FDTD Solver                              #     
        #                13) BendWaveguide FDTD Solver                                  # 
        #                14) ArcWaveguide FDTD Solver                                   #        
        #                15) GratingCoupler FDTD Solver                                 #         
        #                16) RingGratingCoupler FDTD Solver                             #  
        #                                                                               #  
        #                To print information about the solver you desire please        #  
        #                use obj.Help({"Solvers"}: Number). For Example                 #  
        #                obj.Help({"Solvers": 1}) will give you information about       #  
        #                FDE Solver Parameters for Waveguide Structure                  #
        # ============================================================================= #
        """)
    def Solvers(self, NumberSolver):
        if NumberSolver == 1:
            print("""
                    Call Waveguide FDE Solver with -> obj.Solver("Waveguide", "FDE",  Parameters)
                        
                    Dictionary Parameters:
-----------------------------------------------------------------------------------------------------------
                    Parameters['Substrate Height']: int/float   
                        Substrate height
                    Parameters['WG Height'] : int/float
                        Heigth of waveguide
                    Parameters['y res'] : int/float
                        Mesh cell sizes.
                    Parameters['z res'] : int/float
                        Mesh cell sizes.
                    Parameters['Slab Height'] : int/float
                        Height of the slab.
                    Parameters['Wavelength'] : int/float
                        Wavelength.
-----------------------------------------------------------------------------------------------------------
                """)                 
        elif NumberSolver == 2:
            print("""
                    Call MMI2x1 EME Solver with -> obj.Solver("MMI2x1", "EME", Parameters)
                        
                    Dictionary Parameters:
-----------------------------------------------------------------------------------------------------------       
                    Parameters['Substrate Height'] : int/float
                        Height of the slab.
                    Parameters['MMI Width'] : int/float
                        Width of MMI
                    Parameters['MMI Length'] : int/float
                        Length of MMI
                    Parameters['angle'] : int/float
                        Angle of the Waveguide Walls. it is calculated WG_angle = 90 - angle.
                        For anfle = 90 we get a perfect rect!
                    Parameters['WG Height'] : int/float
                        Heigth of waveguide
                    Parameters['WG Width'] : int/float
                        Width og waveguide
                    Parameters['WG Length'] : int/float
                        Length of waveguide
                    Parameters['Position Offset'] : int/float
                        Offset between the waveguides. If Taper == True then this become the offset
                        betweent he tapers wider sides. Waveguide and Tapers cannot be placed ourside
                        the MMI structure. The minimum distance between Taper and Waveguide is 1 um
                        becouse of manufactering restrictions in the University.
                    Parameters['Offset Input'] : int/float
                        Input waveguide/taper offset.
                    Parameters["Taper"]: boolen 
                        Add Taper to the structure on the input and output waveguids
                    Parameters['Taper Length']: int/float
                        Lenght of the Taper in Parameters["Taper"] = True
                    Parameters['y res'] : int/float
                        Mesh cell sizes.
                    Parameters['z res'] : int/float
                        Mesh cell sizes.
                    Parameters['Slab Height'] : int/float
                        Height of the slab.
                    Parameters['Wavelength'] : int/float
                        Wavelength.
                    Parameters["Mode"] : str
                        Mode to choose from ("fundamental TE mode", "fundamental TM mode", "fundamental mode")
                    Parameters["Port Span"]: list of int/floats
                        Parameters["Port Span"] = [Span of Port in x direction, Span of Port in y direction, Span of Port in z direction]
                    Parameters["Offset Output"] : anything, optional
                        This function will allow the user to move the outputs in oposite direction. Please dont use it since is there only 
                        becouse the maschine of our physic departmant had some proiblems with the LNOI objects design. 
-----------------------------------------------------------------------------------------------------------
                """)       
        elif NumberSolver == 3:
            print("""
                    Call MMI2x2 EME Solver with -> obj.Solver("MMI2x2", "EME", Parameters)
                        
                    Dictionary Parameters:
-----------------------------------------------------------------------------------------------------------      
                    Parameters['Substrate Height'] : int/float
                        Height of the slab.
                    Parameters['MMI Width'] : int/float
                        Width of MMI
                    Parameters['MMI Length'] : int/float
                        Length of MMI
                    Parameters['WG Height'] : int/float
                        Heigth of waveguide
                    Parameters['WG Width'] : int/float
                        Width og waveguide
                    Parameters['WG Length'] : int/float
                        Length of waveguide
                    Parameters['Position Offset'] : int/float
                        Offset between the waveguides. If Taper == True then this become the offset
                        betweent he tapers wider sides. Waveguide and Tapers cannot be placed ourside
                        the MMI structure. The minimum distance between Taper and Waveguide is 1 um
                        becouse of manufactering restrictions in the University.
                    Parameters["Taper"]: boolen 
                        Add Taper to the structure on the input and output waveguids
                    Parameters['Taper Length']: int/float
                        Lenght of the Taper in Parameters["Taper"] = True
                    Parameters['Slab Height'] : int/float
                        Height of the slab.
                    Parameters['y res'] : int/float
                        Mesh cell sizes.
                    Parameters['z res'] : int/float
                        Mesh cell size.
                    Parameters['Wavelength'] : int/float
                        Wavelength.
                    Parameters["Mode"] : str
                        Mode to choose from ("fundamental TE mode", "fundamental TM mode", "fundamental mode")
                    Parameters["Port Span"]: list of int/floats
                        Parameters["Port Span"] = [Span of Port in x direction, Span of Port in y direction, Span of Port in z direction]
-----------------------------------------------------------------------------------------------------------
                """)
        elif NumberSolver == 4:
            print("""
                    Call InverseTaper EME Solver with -> obj.Solver("InverseTaper", "EME", Parameters)
                        
                    Dictionary Parameters:
-----------------------------------------------------------------------------------------------------------     
                    Parameters['Substrate Height'] : int/float
                        Substrate height.
                    Parameters['WG Height' : int/float
                        Waveguide hight. Also the height of the MMI section
                    Parameters['WG Width'] : int/float
                        Waveguide width.
                    Parameters['Slab Height'] : int/float
                        Slab height
                    Parameters['PWB Taper Width Back'] : int/float
                        Photonic Wirebonding (PWB) Width back side (to the Photonic Wire Bonding)
                    Parameters['PWB Taper Hight Back'] : int/float
                        Photonic Wire Bonding Height back side (to the Photonic Wire Bonding)
                    Parameters['PWB Taper Length'] : int/float
                        Length of the Photonic Wire Bonding Taper
                    Parameters["SMF Core Diameter"] : int/float
                        Single Mode Fiber core Diameter
                    Parameters["SMF Cladding Diameter"] : int/float
                        Single Mode Fiber Cladding Diameter
                    Parameters['y res'] : int/float
                        Mesh y-Axis
                    Parameters['z res'] : int/float
                        Mesh z-Axis
                    Parameters['Wavelength'] : int/float
                        Wavelength
                    Parameters["Mode"] : str
                        Mode to choose from ("fundamental TE mode", "fundamental TM mode", "fundamental mode")
                    Parameters["Port Span"] : list of floats/ints
                        List of x,y and z span of the Ports. For this simulation only y and z parametes will be taken. 
-----------------------------------------------------------------------------------------------------------
                """)
        elif NumberSolver == 5:
            print("""
                    Call InverseTaper EME Solver with -> obj.Solver("DirectionalCoupler", "EME", Parameters)

                    Dictionary Parameters:
        -----------------------------------------------------------------------------------------------------------     
                    Parameters['Substrate Height'] : float/int
                        Height of the Substrate
                    Parameters['Substrate Width'] : float/int
                        Width of the MMI
                    Parameters['DC Length'] : float/int
                        Length of the Directional coupler
                    Parameters['WG Height'] : float/int
                        Height of the Waveguide
                    Parameters['WG Width'] : float/int
                        Waveguide Width
                    Parameters['Position Offset'] : float/int
                        Positional offser of the waveguides. If posOffset the two Waveguides
                        will be offset of the middle position (y = 0) by the half of there
                        Width. In this case they will not overlap if the Offset is 0.
                    Parameters['y res'] : float/int
                        Mesh resolution for the y-Axis
                    Parameters['z res'] : float/int
                        Mesh resolution for the z Axis
                    Parameters['Slab Height'] : float/int
                        Slab height.
                    Parameters['Wavelength'] : float/int
                        Wavelength
                    Parameters["Port Span"]: list of int/floats
                        Parameters["Port Span"] = [Span of Port in x direction, Span of Port in y direction, Span of Port in z direction]
                    Parameters["Mode"] : str
                        Mode to choose from ("fundamental TE mode", "fundamental TM mode", "fundamental mode") 
        -----------------------------------------------------------------------------------------------------------
                        """)
        elif NumberSolver == 6:
            print("""
                    Call StraightWaveguide EME Solver with -> obj.Solver("StraightWaveguide", "EME", Parameters)
                        
                    Dictionary Parameters:
-----------------------------------------------------------------------------------------------------------          
                    Parameters['Substrate Height'] : int/float
                        Substrate height.
                    Parameters['WG Length'] : int/float
                       Waveguide Length
                    Parameters['WG Height'] : int/float
                       Waveguide hight. Also the height of the MMI section
                    Parameters['WG Width'] : int/float
                       Waveguide width.
                    Parameters["Taper"] : boolen
                       If Taper == False, only straight Waveguide will be simulated,
                       If Taper == True an Taper will be simulated
                    Parameters['Taper Width'] : int/float
                       Taper backside Width. Taper Fronside width is the width of the Waveguide
                    Parameters['Taper Length'] : int/float
                       Taper Length
                    Parameters['y res']: int/float
                         EME Mesh resolutio,
                    Parameters['z res']: int/float
                         EME Mesh resolutio,
                    Parameters['Slab Height'] : int/float
                       Height of the slab.
                    Parameters['Wavelength'] : int/float
                       Wavelength
                    Parameters["Waveguide Angle"] : int/float
                       This Parameter will set the theta ratation angle of the port. It can be 90 or 180.
                    Parameters["Port Span"] : list of floats/ints
                       List of x,y and z span of the Ports. For this simulation only y and z parametes will be taken.
                    Parameters["Mode"] : str
                        Mode to choose from ("fundamental TE mode", "fundamental TM mode", "fundamental mode")
                    Parameters["Taper Type"] : anything, optional
                            This function will check if you have set Parameters["Taper Type"] to anaything, for example "Parameters["Taper Type"]=1" 
                            and if so it will design an Inverse Taper Structure with no Cladding. Here the option "Cladding" is not active and will be ignored.
                            If the user didnt give the "Taper Type" as dictionary key, then an normal taper structure will be simulated.
                            
                            If Parameters["Taper Type"] is given, themn the user need to set couple more parameters:
                                Parameters['PWB Taper Width Back'] : int/float
                                    Photonic Wirebonding (PWB) Width back side (to the Photonic Wire Bonding)
                                Parameters['PWB Taper Hight Back'] : int/float
                                    Photonic Wire Bonding Height back side (to the Photonic Wire Bonding)
                                Parameters['PWB Taper Width Front'] : int/float
                                    Photonic Wirebonding (PWB) Width front side (to the photonic waveguide)
                                Parameters['PWB Taper Hight Front'] : int/float
                                    Photonic Wire Bonding Height front side (to the photonic waveguide)
                                Parameters['PWB Taper Length'] : int/float
                                    Length of the Photonic Wire Bonding Taper
-----------------------------------------------------------------------------------------------------------
                """) 
        elif NumberSolver == 7:
            print("""
                    Call Wavelength Division Multiplexing EME Solver with -> obj.Solver("WDM", "EME", Parameters)
                        
                    Dictionary Parameters:
-----------------------------------------------------------------------------------------------------------       
                    Parameters['Substrate Height'] : int/float
                        Substrate height.
                    Parameters['MMI Width'] : int/float
                        Width of the MMI.
                    Parameters['MMI Length'] : int/float
                        Length of the MMI.
                    Parameters['WG Height' : int/float
                        Waveguide hight. Also the height of the MMI section
                    Parameters['WG Length'] : int/float
                        Waveguide length.
                    Parameters['WG Width'] : int/float
                        Waveguide width.
                    Parameters['Taper Width'] : int/float
                        Taper backside width, frontside width is the waveguide width.
                    Parameters['Taper Length'] : int/float
                        Taper Length
                    Parameters['y res'] : int/float
                        Mesh y-Axis
                    Parameters['z res'] : int/float
                        Mesh z-Axis
                    Parameters['Slab Height'] : int/float
                        Slab Height.
                    Parameters['Wavelength'] : int/float
                        Wavelength
                    Parameters['Angle Thetha'] : boolen
                        Angle for the input and output waveguides
                    Parameters["Mode"] : str
                        Mode to choose from ("fundamental TE mode", "fundamental TM mode", "fundamental mode")
                    Parameters["Port Span"]: list of int/floats
                        Parameters["Port Span"] = [Span of Port in x direction, Span of Port in y direction, Span of Port in z direction]
-----------------------------------------------------------------------------------------------------------
                """) 
        elif NumberSolver == 8:
            print("""
                    Call MMI2x1 FDTD Solver with -> obj.Solver("MMI2x1", "FDTD", Parameters)
                        
                    Dictionary Parameters:
-----------------------------------------------------------------------------------------------------------    
                    Parameters['Substrate Height'] : int/float
                        Height of the slab.
                    Parameters['MMI Width'] : int/float
                        Width of MMI
                    Parameters['MMI Length'] : int/float
                        Length of MMI
                    Parameters['WG Height'] : int/float
                        Heigth of waveguide
                    Parameters['WG Width'] : int/float
                        Width og waveguide
                    Parameters['WG Length'] : int/float
                        Length of waveguide
                    Parameters['Position Offset'] : int/float
                        Offset between the waveguides. If Taper == True then this become the offset
                        betweent he tapers wider sides. Waveguide and Tapers cannot be placed ourside
                        the MMI structure. The minimum distance between Taper and Waveguide is 1 um
                        becouse of manufactering restrictions in the University.
                    Parameters['Offset Input'] : int/float
                        Input waveguide/taper offset.
                    Parameters["Taper"]: boolen 
                        Add Taper to the structure on the input and output waveguids
                    Parameters['Taper Length']: int/float
                        Lenght of the Taper in Parameters["Taper"] = True
                    Parameters['x res'] : int/float
                        Mesh cell sizes.
                    Parameters['Slab Height'] : int/float
                        Height of the slab.
                    Parameters['Wavelength'] : int/float
                        Wavelength.
                    Parameters["Mode"] : str
                        Mode to choose from ("fundamental TE mode", "fundamental TM mode", "fundamental mode")
                    Parameters["Port Span"]: list of int/floats
                        Parameters["Port Span"] = [Span of Port in x direction, Span of Port in y direction, Span of Port in z direction]
                    Parameters["Offset Output"] : anything, optional
                        This function will allow the user to move the outputs in oposite direction. Please dont use it since is there only 
                        becouse the maschine of our physic departmant had some proiblems with the LNOI objects design.    
-----------------------------------------------------------------------------------------------------------
                """) 
        elif NumberSolver == 9:
            print("""
                    Call MMI2x2 FDTD Solver with -> obj.Solver("MMI2x2", "FDTD", Parameters)
                        
                    Dictionary Parameters:
-----------------------------------------------------------------------------------------------------------    
                    Parameters['Substrate Height'] : int/float
                        Height of the slab.
                    Parameters['MMI Width'] : int/float
                        Width of MMI
                    Parameters['MMI Length'] : int/float
                        Length of MMI
                    Parameters['WG Height'] : int/float
                        Heigth of waveguide
                    Parameters['WG Width'] : int/float
                        Width og waveguide
                    Parameters['WG Length'] : int/float
                        Length of waveguide
                    Parameters['Position Offset'] : int/float
                        Offset between the waveguides. If Taper == True then this become the offset
                        betweent he tapers wider sides. Waveguide and Tapers cannot be placed ourside
                        the MMI structure. The minimum distance between Taper and Waveguide is 1 um
                        becouse of manufactering restrictions in the University.
                    Parameters["Taper"]: boolen 
                        Add Taper to the structure on the input and output waveguids
                    Parameters['Taper Length']: int/float
                        Lenght of the Taper in Parameters["Taper"] = True
                    Parameters['x res'] : int/float
                        Mesh cell sizes.
                    Parameters['Slab Height'] : int/float
                        Height of the slab.
                    Parameters['Wavelength'] : int/float
                        Wavelength.
                    Parameters["Mode"] : str
                        Mode to choose from ("fundamental TE mode", "fundamental TM mode", "fundamental mode")
                    Parameters["Port Span"]: list of int/floats
                        Parameters["Port Span"] = [Span of Port in x direction, Span of Port in y direction, Span of Port in z direction]
-----------------------------------------------------------------------------------------------------------
                """) 
        elif NumberSolver == 10:
            print("""
                    Call InverseTaper FDTD Solver with -> obj.Solver("InverseTaper", "FDTD", Parameters)
                        
                    Dictionary Parameters:
-----------------------------------------------------------------------------------------------------------  
                    Parameters['Substrate Height'] : int/float
                      Substrate height.
                    Parameters['WG Height' : int/float
                      Waveguide hight. Also the height of the MMI section
                    Parameters['WG Width'] : int/float
                      Waveguide width.
                    Parameters['Slab Height'] : int/float
                      Slab height
                    Parameters['PWB Taper Width Back'] : int/float
                      Photonic Wirebonding (PWB) Width back side (to the Photonic Wire Bonding)
                    Parameters['PWB Taper Hight Back'] : int/float
                      Photonic Wire Bonding Height back side (to the Photonic Wire Bonding)
                    Parameters['PWB Taper Length'] : int/float
                      Length of the Photonic Wire Bonding Taper
                    Parameters["SMF Core Diameter"] : int/float
                    Single Mode Fiber core Diameter
                    Parameters["SMF Cladding Diameter"] : int/float
                    Single Mode Fiber Cladding Diameter
                    Parameters['x res'] : int/float
                      Mesh x-Axis
                    Parameters['Wavelength'] : int/float
                      Wavelength
                    Parameters["Mode"] : str
                      Mode to choose from ("fundamental TE mode", "fundamental TM mode", "fundamental mode")
                    Parameters["Port Span"] : list of floats/ints
                      List of x,y and z span of the Ports. For this simulation only y and z parametes will be taken.                 
-----------------------------------------------------------------------------------------------------------
                """)                         
        elif NumberSolver == 11:
            print("""
                    Call Directional FDTD Solver with -> obj.Solver("DirectionalCoupler", "FDTD", Parameters)
                        
                    Dictionary Parameters:
-----------------------------------------------------------------------------------------------------------  
                    Parameters['Substrate Height'] : float/int
                        Height of the Substrate
                    Parameters['Substrate Width'] : float/int
                        Width of the MMI
                    Parameters['DC Length'] : float/int
                        Length of the Directional coupler
                    Parameters['WG Height'] : float/int
                        Height of the Waveguide
                    Parameters['WG Width'] : float/int
                        Waveguide Width
                    Parameters['Position Offset'] : float/int
                        Positional offser of the waveguides. If posOffset the two Waveguides
                        will be offset of the middle position (y = 0) by the half of there
                        Width. In this case they will not overlap if the Offset is 0.
                    Parameters['x res'] : float/int
                        Mesh resolution for the x-Axis
                    Parameters['Slab Height'] : float/int
                        Slab height.
                    Parameters['Wavelength'] : float/int
                        Wavelength
                    Parameters["Port Span"]: list of int/floats
                        Parameters["Port Span"] = [Span of Port in x direction, Span of Port in y direction, Span of Port in z direction]
                    Parameters["Mode"] : str
                        Mode to choose from ("fundamental TE mode", "fundamental TM mode", "fundamental mode")
-----------------------------------------------------------------------------------------------------------
                """)   
        elif NumberSolver == 12:
            print("""
                    Call StraightWaveguide FDTD Solver with -> obj.Solver("StraightWaveguide", "FDTD", Parameters)
                        
                    Dictionary Parameters:
-----------------------------------------------------------------------------------------------------------         
                    Parameters['Substrate Height'] : int/float
                       Substrate height.
                    Parameters['WG Length'] : int/float
                       Waveguide Length
                    Parameters['WG Height'] : int/float
                       Waveguide hight. Also the height of the MMI section
                    Parameters['WG Width'] : int/float
                       Waveguide width.
                    Parameters["Taper"] : boolen
                       If Taper == False, only straight Waveguide will be simulated,
                       If Taper == True an Taper will be simulated
                    Parameters['Taper Width'] : int/float
                       Taper backside Width. Taper Fronside width is the width of the Waveguide
                    Parameters['Taper Length'] : int/float
                       Taper Length
                    Parameters['x res'] : int/float
                         EME Mesh resolutio,
                    Parameters['Slab Height'] : int/float
                       Height of the slab.
                    Parameters['Wavelength'] : int/float
                       Wavelength
                    Parameters["Waveguide Angle"] : int/float
                       This Parameter will set the theta ratation angle of the port. It can be 90 or 180.
                    Parameters["Port Span"] : list of floats/ints
                       List of x,y and z span of the Ports. For this simulation only y and z parametes will be taken.
                    Parameters["Mode"] : str
                        Mode to choose from ("fundamental TE mode", "fundamental TM mode", "fundamental mode")
                    Parameters["Taper Type"] : anything, optional
                            This function will check if you have set Parameters["Taper Type"] to anaything, for example "Parameters["Taper Type"]=1" 
                            and if so it will design an Inverse Taper Structure with no Cladding. Here the option "Cladding" is not active and will be ignored.
                            If the user didnt give the "Taper Type" as dictionary key, then an normal taper structure will be simulated.
                            
                            If Parameters["Taper Type"] is given, themn the user need to set couple more parameters:
                                Parameters['PWB Taper Width Back'] : int/float
                                    Photonic Wirebonding (PWB) Width back side (to the Photonic Wire Bonding)
                                Parameters['PWB Taper Hight Back'] : int/float
                                    Photonic Wire Bonding Height back side (to the Photonic Wire Bonding)
                                Parameters['PWB Taper Width Front'] : int/float
                                    Photonic Wirebonding (PWB) Width front side (to the photonic waveguide)
                                Parameters['PWB Taper Hight Front'] : int/float
                                    Photonic Wire Bonding Height front side (to the photonic waveguide)
                                Parameters['PWB Taper Length'] : int/float
                                    Length of the Photonic Wire Bonding Taper
-----------------------------------------------------------------------------------------------------------
                """) 
        elif NumberSolver == 13:
            print("""
                    Call BendWaveguide FDTD Solver with -> obj.Solver("BendWaveguide", "FDTD", Parameters)
                        
                    Dictionary Parameters:
-----------------------------------------------------------------------------------------------------------   
                    Parameters['Substrate Height'] : int/float
                       Substrate height.
                    Parameters['WG Height'] : int/float
                       Waveguide hight. Also the height of the MMI section
                    Parameters['x res'] : int/float
                         EME Mesh resolutio,
                    Parameters['Slab Height'] : int/float
                       Height of the slab.
                    Parameters['Wavelength'] : int/float
                       Wavelength
                    Parameters["x span"] : int/float
                       Length of the S-Bend Waveguide
                    Parameters["y span"] : int/float
                       Width of the S-Bend Waveguide
                    Parameters["Mode"] : str
                        Mode to choose from ("fundamental TE mode", "fundamental TM mode", "fundamental mode")
                    Parameters["Port Span"]: list of int/floats
                        Parameters["Port Span"] = [Span of Port in x direction, Span of Port in y direction, Span of Port in z direction]
-----------------------------------------------------------------------------------------------------------
                """) 
        elif NumberSolver == 14:
            print("""
                    Call ArcWaveguide FDTD Solver with -> obj.Solver("ArcWaveguide", "FDTD", Parameters)
                        
                    Dictionary Parameters:
-----------------------------------------------------------------------------------------------------------     
                    Parameters['Substrate Height'] : int/float
                        Substrate height.
                    Parameters['WG Height'] : int/float
                        Waveguide hight. Also the height of the MMI section
                    Parameters['WG Width'] : int/float
                        Waveguide width.
                    Parameters['x res'] : int/float
                        EME Mesh resolutio,
                    Parameters['Slab Height'] : int/float
                        Height of the slab.
                    Parameters['Wavelength'] : int/float
                        Wavelength
                    Parameters["S_Band Radius"] : int/float
                        S-Bend Radius in um.
                    Parameters['Arc deg'] : int/float
                        Arc define the Arc of the curve. It can be 90 or 180 degrees only.
                        This two will define an 1/4 of a circle or 1/2 of a circle.
                    Parameters["Mode"] : str
                        Mode to choose from ("fundamental TE mode", "fundamental TM mode", "fundamental mode")
                    Parameters["Port Span"]: list of int/floats
                        Parameters["Port Span"] = [Span of Port in x direction, Span of Port in y direction, Span of Port in z direction]
-----------------------------------------------------------------------------------------------------------
                """)
        elif NumberSolver == 15:
            print("""
                    Call GratingCoupler FDTD Solver with -> obj.Solver("GratingCoupler", "FDTD", Parameters)
                        
                    Dictionary Parameters:
-----------------------------------------------------------------------------------------------------------    
                    Parameters['Substrate Height'] : int/float
                        Substrate height.
                    Parameters["Length GC"]: int/float
                        Lenght of the Grating Coupler Area
                    Parameters["Input Length GC"]: int/float
                        An squere Waveguide with the same WG Height as the Grating coupler place before the Grating Coupler region will start. 
                    Parameters["Output Length GC"]: int/float
                        An squere Waveguide with the same WG Height as the Grating coupler place after the Grating Coupler region to finish the structure.
                    Parameters["Width GC"]: int/float
                        Widht of the Grating Coupler Area
                    Parameters["Hight GC"]: int/float
                        Hight of the Grating Coupler Material
                    Parameters['Taper'] : boolen
                        Add Taper to structure
                    Parameters['Taper Length'] : int/float
                          Length of the Taper
                    Parameters['Wavelength'] : int/float
                          Wavelength
                    Parameters['x res'] : int/float
                          Mesh x-Axis
                    Parameters["Mode"] : str
                        Mode to choose from ("fundamental TE mode", "fundamental TM mode", "fundamental mode")
                    Parameters["Port Span"] : list of floats/ints
                          List of x,y and z span of the Ports. For this simulation only y and z parametes will be taken.
                    Parameters["GC Radius"]: int/float
                        Radius of the Ring Grating Coupler in um. For Example "Parameters["GC Radius"] = 25e-6"
 -----------------------------------------------------------------------------------------------------------
                """)               
        elif NumberSolver == 16:
            print("""
                    Call RingGratingCoupler FDTD Solver with -> obj.Solver("RingGratingCoupler", "FDTD", Parameters)
                        
                    Dictionary Parameters:
-----------------------------------------------------------------------------------------------------------            
                    Parameters['Substrate Height'] : int/float
                        Substrate height.
                    Parameters["Length GC"]: int/float
                        Lenght of the Grating Coupler Area
                    Parameters["Input Length GC"]: int/float
                        An squere Waveguide with the same WG Height as the Grating coupler place before the Grating Coupler region will start. 
                    Parameters["Output Length GC"]: int/float
                        An squere Waveguide with the same WG Height as the Grating coupler place after the Grating Coupler region to finish the structure.
                    Parameters["Width GC"]: int/float
                        Widht of the Grating Coupler Area
                    Parameters["Hight GC"]: int/float
                        Hight of the Grating Coupler Material
                    Parameters["GC Radius"]: int/float
                        Radius of the Ring Grating Coupler in um. For Example "Parameters["GC Radius"] = 25e-6"
                    Parameters['Taper Length'] : int/float
                          Length of the input Taper
                    Parameters['Wavelength'] : int/float
                          Wavelength
                    Parameters['x res'] : int/float
                          Mesh x-Axis
                    Parameters["Mode"] : str
                        Mode to choose from ("fundamental TE mode", "fundamental TM mode", "fundamental mode")
                    Parameters["Port Span"] : list of floats/ints
                          List of x,y and z span of the Ports. For this simulation only y and z parametes will be taken. 
 -----------------------------------------------------------------------------------------------------------
                """)                           
     
              
    def Help_StartSimulation(self):
        print("""
            # ============================================================================= #
            #                Welcome to the Help Menu for Lumerical Simulation Run. You     #     
            #                can ask for help for one of the following simulations          #   
            #                                                                               #
            #                1) Save and start FDE Simulation                               #
            #                2) Save and start EME Simulation                               # 
            #                3) Save and start FDTD Simulation                              #
            #                                                                               #
            #                To print information about how to start an simulation in the   #  
            #                choosen Solver. With obj.Help({"Start Simulation"}, Number) an #
            #                save of your simulation will be made in your current working   #
            #                directory. Afterwards the Simulation will be started.For       #
            #                Example obj.Help({"Start Simulation": 1}) will start the FDE   #
            #                simulation.                                                    #
            # ============================================================================= #
            """)
    def StartSolver(self, NumberStart):
        if NumberStart == 1:
            print("""
                    Start FDE Simulation -> obj.StartFDESimulation() 
                        
                    Dictionary Parameters:
-----------------------------------------------------------------------------------------------------------      
                    None  
 -----------------------------------------------------------------------------------------------------------
                """) 
        elif NumberStart == 2:
            print("""
                    Start EME Simulation -> obj.StartEMESimulation() 
                        
                    Dictionary Parameters:
-----------------------------------------------------------------------------------------------------------      
                    None  
 -----------------------------------------------------------------------------------------------------------
                """) 
        elif NumberStart == 3:
            print("""
                    Start FDTD Simulation -> obj.StartFDTDSimulation() 
                        
                    Dictionary Parameters:
-----------------------------------------------------------------------------------------------------------      
                    None  
 -----------------------------------------------------------------------------------------------------------
                """) 

              
    def Help_Results(self):
        print("""
            # ============================================================================= #
            #                Welcome to the Help Menu for Lumerical Result Extraction.      #     
            #                Youcan ask for help for one of the following Result Extraction #   
            #                 options.                                                      #
            #                                                                               #
            #                1) Extract FDE Simulation Results                              #
            #                2) How to set Overlap Anaylsis with FDE Solver                 # 
            #                3) Extract EME Simulation Results                              #
            #                3) Extract FDTD Simulation Results                             #
            #                                                                               #
            #                To print information about how to extract results from the     #  
            #                simulation after Solver is done. With obj.Help({"Results"},    #
            #                Number) the simulated results will be extracted from the       #
            #                simulation workspace and transfered to python IDE.             #
            #                For Example obj.Help({"Results": 1}) will extract the FDE      #
            #                data into 4 dictionaries.                                      #
            # ============================================================================= #
            """)
    def Result_extraction(self, NumberResults):
        if NumberResults == 1:
            print("""
                    Extract FDE results -> obj.ExtrtactData('FDE', Parameters)
                        
                    Dictionary Parameters:
    -----------------------------------------------------------------------------------------------------------     
                    Parameters["Effective Index"] : float/int
                     Effective index value. All modes with effective index smaller then
                     the  Effective will be deleted. Effective is usualy the smalles effective
                     index from the hole construction (usualy the cladding). Since Lumerical will produce 
                     modes with very small Effektive index, it make no sance to extract them !
    -----------------------------------------------------------------------------------------------------------     


                    Returns
    -----------------------------------------------------------------------------------------------------------     
         dictModesTE : dict
                 Dictionary with:
                         'TE polarization fraction num'
                         'effective area'
                         'effective index num'
                         'group index num'
                         'Ex'
                         'Ey'
                         'Ez'
                         
         dictDataTE : dict
                 Dictionary with:
                             'E'
                             
                             
        dictModesTM : dict
                Dictionary with:
                        'TE polarization fraction num'
                        'effective area'
                        'effective index num'
                        'group index num'
                        'Ex'
                        'Ey'
                        'Ez'
                        
        dictDataTM : dict
                Dictionary with:
                            'E'  
    -----------------------------------------------------------------------------------------------------------
                """)
        elif NumberResults == 2:
            print("""
                    Performe an Mode Overlap analisys FDE 
                                    How to:
    -----------------------------------------------------------------------------------------------------------    
                After an FDE Simulation the user can:
                 
                 - Performe an Mode Overlap analisys. 
                 Here an Waveguide with an choosen Parameters[] will be placed in the FDE solver. 
                 Mode analysis is performed and the two fundamental modes TE and TM modes will be clipped to the Lumerical Board.
                 After the simulation the object will be removed. The two modes stay clipped to the board. 
                 The user can then place an new Waveguide and performe and FDE Analysis. After the analysis 
                 the modes can be compaired (overlapped) in Lumerical with the already clipped modes from the 
                 first simulation. 
    
                 With - OverlapFDEModes(Parameters)
                 Parameters["Effective Index"] = Reflective Effektive index for the Material, Effective index value. 
                 All modes with effective index smaller then the  Effective will be deleted. Effective is usualy the 
                 smalles effective index from the hole construction (usualy the cladding). Since Lumerical will produce 
                 modes with very small Effektive index, it make no sance to extract them ! 
    -----------------------------------------------------------------------------------------------------------
                """)
        elif NumberResults == 3:
            print("""
                    Extract EME results -> obj.ExtrtactData('EME', Parameters)  
                        
                    Dictionary Parameters:
    -----------------------------------------------------------------------------------------------------------   
                    None
    -----------------------------------------------------------------------------------------------------------     


                    Returns
    -----------------------------------------------------------------------------------------------------------   
                   MonitorData : dict
                       Dictionary with:
                               'lambda'
                               'f'
                               'x'
                               'y'
                               'z'
                               'E'
                               'H'
                               'Lumerical_dataset'
                   EME_Sim_Data : dict
                       Dictionary with:
                               'user s matrix'
                               'internal s matrix'
                               'local diagnostics'
                               'coefficients'
                               'global diagnostics'
    -----------------------------------------------------------------------------------------------------------
                """)
        elif NumberResults == 4:
            print("""
                    Extract FDTD results -> obj.ExtrtactData('FDTD', Parameters)
                        
                    Dictionary Parameters:
    -----------------------------------------------------------------------------------------------------------
                    Parameters: boolen
                        If Parameters = True, an extra S-Parameter sweep over the ports will be performed. This 
                        will add simulation time since the structute need to be sweeped one more time.
                        If Parameters = False, only FDTD simulation will be done and result will be extracted 
                        the "SParam_Dict" dictionary will return 0.
    -----------------------------------------------------------------------------------------------------------     


                    Returns
    -----------------------------------------------------------------------------------------------------------  
                    SParam_Dict: dict
                        Dictionary with:
                                "S Parameter"
                            
                    PowerData_Dict : dict
                        Dictionary with:
                                'Power Input Port/s'
                                'Power Output Port/s'
                                'Transmission'
                     
                
    -----------------------------------------------------------------------------------------------------------
                """)

    
     
              
    def Help_LoadingBar(self):
        print("""
            # ============================================================================= #
            #                Welcome to the Help Menu for the loading bar function          #     
            #                Here you can find a short discription on how to use it.        #   
            #                                                                               #
            #                If the user want to keep track of the simulation state an      #
            #                laoding bar is imported in this library. The loading bar is    #
            #                not part of the class, so you dont need to call it like an     #
            #                class.                                                         #
            #                                                                               #  
            #                loadingBar(count, total)                                       #
            #                count - is your loop variable                                  # 
            #                total - is the lenght of the variables that will be sweeped    #
            #                                                                               #
            #                Example python code of loading bar:                            # 
            #                                                                               #
            #      >from Constructor import loadingBar                                      #
            #      >                                                                        #
            #      >for i in range(10):                                                     #
            #      >    loadingBar(i, 10)                                                   #
            #                                                                               #
            # ============================================================================= #
            """)
        
        

    def Help_RemoveObject(self):
        print("""
            # ============================================================================= #
            #                Welcome to the Help Menu for how to remove simulation          #     
            #                from Lumrical                                                  #   
            #                                                                               #
            #                To Remove an Object simply use:                                #
            #                                                                               #
            #                >obj.removeObject()                                            #
            #                                                                               #  
            #                This will switch Lumerical from Simulation to                  #
            #                design and delete all the Structures, Monitors and Solver      #
            #                created by the user.                                           #
            #                                                                               #
            # ============================================================================= #
            """)
 
        
    def Help_LogFile(self):
        print("""
            # ============================================================================= #
            #                Welcome to the Help Menu for creating simulation log File      #
            #                from Lumrical                                                  #
            #                                                                               #
            #                When you load the Logfile function from the Constructor Library#
            #                the user can keep track with the simulation and parameters     #
            #                used. To use the log File you need to give an list with        #
            #                two dictionarays. The first can be the dictionary with         #
            #                Parameters given for the simulation like "Tapper Width",       #
            #                "WG Width"... The secound dictionary can be called direct      #
            #                after performing an analysis with obj.ReturnLogInfo(). This    #
            #                will extract the object that you simulate and the solver that  #
            #                you used. An text file will be created in the given Path. The  #
            #                name of the File will be "Lumerical_Simulation_Log_15_10_52"   #
            #                                                                               #
            #                An example can be seen here:                                   #
            #                                                                               #
            #                >from Constructor import Logfile                               #
            #                >Path = "C:/Downloads/"                                        #
            #                >solverInfo = obj.ReturnLogInfo()                              #
            #                >data = [solverInfo, Parameters]                               #
            #                >Logfile(data, Path)                                           #
            #                                                                               #
            # ============================================================================= #
            """)
    



def Help():
    print("""
    
        To start the Python Lumerical Constructor Library you need first to:
        from Constructor import Constructor
        from Constructor import loadingBar
        
        # Call the Constructor Object 
        obj = Constructor(file, Lumerical Solver Type)
        # or 
        obj = Constructor(file, MaterialPath, Lumerical Solver Type)
            # file = Path to lumapi.py for example -> "lumerical/2022-r2.4/api/python/lumapi.py"
            # MaterialPath = Path to Material Library for example -> "/Desktop/LNOI_Env/LiNbO3.mdf"-> THIS IS OPTIONAL
            # Lumerical Solver Type = Solver Type can be EME or FDTD for example -> "EME"
            
        # To Close Lumerical simply use
        obj.Close()
        
        # For Further Help after colling the Constructor Object please use:
        obj.Help()
            
    """)
    
# obj.Help('Objects') -> Shows how to build an photonic element
# obj.Help('Solvers') -> Shows hot to build an FDE, EME or FDTD Solver
# obj.Help('Start Simulation') -> Commands to start simulation
# obj.Help('Results') -> Shows how to extract resuslts form FDE, EME and FDTD Simulations
# obj.Help('Loading Bar') -> Shows how to use the Loading Bar Function 
# obj.Help('Remove Object') -> Romove object function
# obj.Help('Log File') -> Generate an log file from the simulation 



def Logfile(data, path):

    """
    

    Parameters
    ----------
    data : list
        List of dictionarys with all the data for the simulation (obj.ReturnLogInfo()) and the parameters (Parameters dictionary) given to the simulation. 
    Returns
    -------
    None.

    """
    
    
    from datetime import datetime
    # datetime object containing current date and time
    now = datetime.now()
    # dd/mm/YY H:M:S
    dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
    
    
    name = "Lumerical_Simulation_Log_"+(dt_string.split(" ")[1]).replace(":","_")
    
    
    if "EME" or "FDTD" in list(data[0].keys()):
        SimulationData = data[0]
        Parameters = data[1]
    else:
        Parameters = data[0]
        SimulationData = data[1]
        



    ParamKeys = list(Parameters.keys())
    SimulationDataKeys = list(SimulationData.keys())
    
    with open(path + name + '.txt','w') as data: 
        data.write("""
    ==========================================================================================================================================================
                   """ + "\n")
        data.write(f"This is an Simulation Log file from the Lumerical Simulation done on {dt_string}" +"\n")
        data.write("""
    ==========================================================================================================================================================
                   """ + "\n")
                   
        for keys in range(len(Parameters.keys())):
            data.write(ParamKeys[keys] + " : " + str(Parameters[ParamKeys[keys]]) + "\n")
            
            
        data.write("""
    ==========================================================================================================================================================
                   """ + "\n")
        data.write(f"Simulated Abject and used Solver" +"\n")
        data.write("""
    ==========================================================================================================================================================
                   """ + "\n")
        for keys in range(len(SimulationData.keys())):
            data.write(SimulationDataKeys[keys] + " : " + str(SimulationData[SimulationDataKeys[keys]]) + "\n")
    
    
