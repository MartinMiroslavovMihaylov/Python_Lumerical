# -*- coding: utf-8 -*-
"""
Created on Fri Feb 10 08:31:42 2023

@author: Martin.Mihaylov
"""

import imp
import numpy as np
import sys



def loadingBar(count, total, size=1):
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
    def __init__(self, file, MaterialLib, Mode):
        '''
        Path to the umerical python file
        '''
        self.file = file
        self.MaterialLib = MaterialLib
        self.lumpai = imp.load_source('lumapi', self.file)
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
            raise ValueError(
                "Non Valid Solver was choosen. Please pass on one of the two supported solvers ['FDTD' or 'EME']")

    # Close Programm
    def Close(self):
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
        
    
    def Help(self, Subject):

        if Subject is None:
            raise ValueError(
                "Help can be called with Help(str(subject)). Subject can be choosen from 'Objects', 'Solvers', 'Start Simulation', 'Results', 'Loading Bar' and 'Log File' !! ")
        else:
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
        self.lum = self.lumpai.FDTD()
        self.lum.importmaterialdb(self.MaterialLib)
        print('Lumerical FDTD API is started')

    def EME(self):
        self.lum = self.lumpai.MODE()
        self.lum.importmaterialdb(self.MaterialLib)
        print('Lumerical EME API is started')



    #   Solvers Function
    def Solver(self, Structure, Type, Parameters):
        if Type == "FDTD":
            self.FDTD_Solver(Structure, Parameters)
        elif Type == "EME":
            self.EME_Solver(Structure, Parameters)
        elif Type == "FDE":
            self.FDE_Solver(Structure, Parameters)
        elif Type == "varFDTD":
            self.varFDTD_Solver(Structure, Parameters)
        else:
            raise ValueError(
                "Invalid Solver Type! Possible Options are FDE - For Waveguides, EME - for MMI , directional Couplers, varFDTD - for MMI Structures and FDTD - for MMI.")



    def StartSimFunction(self, Type):
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
        if Type == "FDTD":
            if self.Struct == "MMI2x1" or self.Struct == "MMI2x1_Trapez":
                S_Param3, Power3 = self.ExtractFDTDResults(3, Parameters)
                return S_Param3, Power3
            elif self.Struct == "MMI2x2" or self.Struct == "MMI2x2_Trapez":
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
        self.Waveguide(Parameters)
        self.FDE_Solver("Waveguide", Parameters)
        self.StartFDESimulation()
        TEModes, TMModes = self.ExtractFDEModes(EffIndexValue=Parameters["Effective Index"])
        self.CoppyDcard(list(TEModes.keys())[0], list(TMModes.keys())[0])
        self.removeObject()
        print(
            f"TE Mode {list(TEModes.keys())[0]} and TM Mode {list(TMModes.keys())[0]} have been copyed to the Lumerical Dcard.")

    # =============================================================================
    # Selection of possible Functions
    # =============================================================================



    def FDTD_Solver(self, Structure, Parameters):
        if Structure == "MMI2x1" or Structure == "MMI2x1_Trapez":
            self.Struct = "MMI2x1"
            self.setMMI2x1FDTDSolver(Parameters)
            self.SolverInfo["Simulated Object"] = "MMI2x1"
        elif Structure == "MMI2x2" or Structure == "MMI2x2_Trapez":
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
        else:
            raise ValueError("Invalid Strucute for FDTD Solver is selected. Possible Strucures are MMI2x1 or MMI2x2")


    def EME_Solver(self, Structure, Parameters):
        if Structure == "MMI2x1" or Structure == "MMI2x1_Trapez":
            self.setMMI2x1EMESolver(Parameters)
            self.SolverInfo["Simulated Object"] = "MMI2x1"
        elif Structure == "MMI2x2" or Structure == "MMI2x2_Trapez":
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
            raise ValueError(
                "Invalid Structure for EME Solver is selected. Possible Strucutres are MMI2x1, MMI2x2 and DirectionalCoupler")



    def FDE_Solver(self, Structure, Parameters):
        if Structure == "Waveguide":
            self.setWaveguideFDESolver(Parameters)
            self.SolverInfo["Simulated Object"] = "Waveguide"
        else:
            raise ValueError("Invalid Structure for FDE Solver is selected. Possible Structure is Waveguide.")



    def varFDTD_Solver(self, Structure, Parameters):
        if Structure == "MMI2x1":
            self.setMMI1x2VarFDTDSolver(Parameters)
        elif Structure == "MMI2x2":
            self.setMMI2x2VarFDTDSolver(Parameters)
        else:
            raise ValueError("Invalid Strucute for varFDTD Solver is selected. Possible Strucures are MMI2x1 or MMI2x2")
            
            
    def ReturnLogInfo(self):
        return self.SolverInfo
    
    
    def Script(self):
        myscript =  'delta_w = 2*thickness*tan((angle_side)*pi/180); \n'
        myscript = myscript +  '?"width_l = " + num2str(width_l); \n'
        myscript = myscript +  '?"width_r = " + num2str(width_r) + endl; \n'
        myscript = myscript +  'width_top_l = width_l - (1-hfrac_ref)*delta_w; \n'
        myscript = myscript +  'width_top_r = width_r - (1-hfrac_ref)*delta_w; \n'
        myscript = myscript +  '?"width_top_l = " + num2str(width_top_l); \n'
        myscript = myscript +  '?"width_top_r = " + num2str(width_top_r) + endl; \n'
        myscript = myscript +  'width_bot_l = width_l + hfrac_ref*delta_w; \n'
        myscript = myscript +  'width_bot_r = width_r + hfrac_ref*delta_w; \n'
        myscript = myscript +  '?"width_bot_l = " + num2str(width_bot_l); \n'
        myscript = myscript +  '?"width_bot_r = " + num2str(width_bot_r) + endl; \n'
        myscript = myscript +  'if ((hfrac_ref>1) or (hfrac_ref<0)){?"Error: hfrac_ref must be between 0 and 1.";break;} \n'
        myscript = myscript +  'if ((width_top_l<0) or (width_top_r<0) or (width_bot_l<0) or (width_bot_r<0)){?"Error: width and angle values are not correct.";break;} \n'
        myscript = myscript +  'zmin = -thickness/2; \n'
        myscript = myscript +  'zmax = thickness/2; \n'
        myscript = myscript +  'xmin = -len/2; \n'
        myscript = myscript +  'xmax = len/2; \n'
        myscript = myscript +  'ymin_bot_l = -width_bot_l/2; \n'
        myscript = myscript +  'ymax_bot_l = width_bot_l/2; \n'
        myscript = myscript +  'ymin_bot_r = -width_bot_r/2; \n'
        myscript = myscript +  'ymax_bot_r = width_bot_r/2; \n'
        myscript = myscript +  'ymin_top_l = -width_top_l/2; \n'
        myscript = myscript +  'ymax_top_l = width_top_l/2; \n'
        myscript = myscript +  'ymin_top_r = -width_top_r/2; \n'
        myscript = myscript +  'ymax_top_r = width_top_r/2; \n'
        myscript = myscript +  'vtx=    [xmin,ymin_bot_l,zmin; xmax,ymin_bot_r,zmin; xmax,ymax_bot_r,zmin; xmin,ymax_bot_l,zmin; xmin,ymin_top_l,zmax; xmax,ymin_top_r,zmax; xmax,ymax_top_r,zmax; xmin,ymax_top_l,zmax];   \n'
        myscript = myscript +  'a = cell(6); \n'
        myscript = myscript +  'for(i = 1:6){ a{i} = cell(1);}   \n'
        myscript = myscript +  'a{1}{1} = [1,4,3,2];   \n'
        myscript = myscript +  'a{2}{1} = [1,2,6,5];   \n'
        myscript = myscript +  'a{3}{1} = [2,3,7,6];   \n'
        myscript = myscript +  'a{4}{1} = [3,4,8,7];   \n'
        myscript = myscript +  'a{5}{1} = [1,5,8,4];   \n'
        myscript = myscript +  'a{6}{1} = [5,6,7,8];   \n'
        myscript = myscript +  'addplanarsolid(vtx,a);   \n'
        myscript = myscript +  'if (material=="<Object defined dielectric>"){setnamed("solid", "index",index);}   \n'
        myscript = myscript +  'else{setnamed("solid", "material",material);}   \n'
        
        return myscript



# =============================================================================
# Structures
# =============================================================================


    def Waveguide(self, Parameters):
        
        
        '''
        Parameters
        ----------
        Parameters : Dictionary
            Dictionary with all the data needet for the Bend Wavaguide. Data needet:
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
        Parameters['Material'] : list
            List of Materials. the list should be with names (str) of a valid Lumerical materials.
            Check the names in Lumerical Materials viewer.
                   
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
            self.lum.set("z", max_slabH)
            self.lum.set("base width", WG_W)
            self.lum.set("base height", max_slabH + WG_Height / 2)
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
            self.lum.set("z min", min_BoxH)
            self.lum.set("z max", min_slabH)
            self.lum.set("x min", minWGL)
            self.lum.set("x max", maxWGL)
            self.lum.set("material", MaterialSub)

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
            




    def StraightWaveguide(self, Parameters):

        '''
        
        Parameters
        ----------
        Parameters : Dictionary
            Dictionary with all the data needet for the Bend Wavaguide. Data needet:
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
                List of Materials. the list should be with names (str) of a valid Lumerical materials.
                Check the names in Lumerical Materials viewer.
            Parameters["Wavelength"] : int/float
                Wavelength
            Parameters["Waveguide Angle"] : int/float
                Bending angle of the Waveguide. Set it to 0 to simulate the straight waveguide. 
                If Waveguide Angle is different then 0, then the straight waveguide will be tilted 
                at the choosen degrees. 
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
        Device_Width = 2*WG_Length + WaveLength * 2  # MMI_Width

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
    
                # creating the thin film
                min_slabH = max_subH
                max_slabH = max_subH + Slab_Height
    
                self.lum.addrect()
                self.lum.set("name", "LN_slab")
                self.lum.set("y", 0)
                self.lum.set("y span", Device_Width)
                self.lum.set("z min", min_slabH)
                self.lum.set("z max", max_slabH)
                self.lum.set("x", WG_Length / 2)
                self.lum.set("x span", WG_Length + WG_Width)
                self.lum.set("material", MaterialSlab)
    
                z_Offset = max_slabH + WG_Height / 2
                # Triangle EQ for waveguide Width
                x = abs(WG_Height / (np.cos((angle) * np.pi / 180)))  # in Radians
                extention = np.sqrt(x ** 2 - WG_Height ** 2)
                WG_W = WG_Width + 2 * extention
    
    
    
                names = ["Staight_WG"]
    
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
                
            
                # # create_cover
                # self.lum.addrect()
                # self.lum.set("name", "cladding")
                # self.lum.set("material", MaterialClad)
                # self.lum.set("y", 0)
                # self.lum.set("y span", Device_Width)
                # self.lum.set("z min", min_slabH)
                # self.lum.set("z max", min_slabH*2)
                # self.lum.set("x", WG_Length / 2)
                # self.lum.set("x span", WG_Length)
                # self.lum.set("override mesh order from material database", True)
                # self.lum.set("mesh order", 4)
                # self.lum.set("alpha", 0.7)
    
    
    
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
    
                # creating the thin film
                min_slabH = max_subH
                max_slabH = max_subH + Slab_Height
    
                self.lum.addrect()
                self.lum.set("name", "LN_slab")
                self.lum.set("y", WG_Length / 4)
                self.lum.set("y span", Device_Width)
                self.lum.set("z min", min_slabH)
                self.lum.set("z max", max_slabH)
                self.lum.set("x", WG_Length / 2)
                self.lum.set("x span", WG_Length + WG_Width)
                self.lum.set("material", MaterialSlab)
    
                z_Offset = max_slabH + WG_Height / 2
                # Triangle EQ for waveguide Width
                x = abs(WG_Height / (np.cos((angle) * np.pi / 180)))  # in Radians
                extention = np.sqrt(x ** 2 - WG_Height ** 2)
                WG_W = WG_Width + 2 * extention
    
    
    
                names = ["Staight_WG"]
    
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
                
                
                
                
                # # create_cover
                # self.lum.addrect()
                # self.lum.set("name", "cladding")
                # self.lum.set("material", MaterialClad)
                # self.lum.set("y", WG_Length / 4)
                # self.lum.set("y span", Device_Width)
                # self.lum.set("z min", min_slabH)
                # self.lum.set("z max", min_slabH*2)
                # self.lum.set("x", WG_Length / 2)
                # self.lum.set("x span", WG_Length + WG_Width)
                # self.lum.set("override mesh order from material database", True)
                # self.lum.set("mesh order", 4)
                # self.lum.set("alpha", 0.7)
                
                
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

            # creating the thin film
            min_slabH = max_subH
            max_slabH = max_subH + Slab_Height

            self.lum.addrect()
            self.lum.set("name", "LN_slab")
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

            
            





    def BendWaveguide(self, Parameters):
        '''


        Parameters
        ----------
        Parameters : Dictionary
            Dictionary with all the data needet for the Bend Wavaguide. Data needet:
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
                Bending angle of the Waveguide. Set it to 0 to simulate the straight waveguide. 
                If Waveguide Angle is different then 0, then the straight waveguide will be tilted 
                at the choosen degrees. 
            Parameters["Wavelength"] : int/float
                Wavelength
            Parameters["x span"] : int/float
                Length of the S Curve. Span of the object.
            Parameters["y span"] : int/float
                Height of the curve. Difference between the input and output of the S-curve.
            Parameters["poles"] : boolen
                If Parameters["poles"] = True an Bezier Curbe will be made.
                if Parameters["poles"] = False an Cosinus Curve = y_span*(cos((pi/(2*x_span))*t)^2) will be made. Where
                t is in the range of 0 - y_span
                
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
        WaveLength = Parameters["Wavelength"]
        x_span = Parameters["x span"]
        y_span = Parameters["y span"]
        polesList = Parameters["poles"]
        
        
        
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
        Device_Width = y_span + WaveLength * 2  # MMI_Width

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

        # creating the thin film
        min_slabH = max_subH
        max_slabH = max_subH + Slab_Height

        self.lum.addrect()
        self.lum.set("name", "LN_slab")
        self.lum.set("y", y_span / 2)
        self.lum.set("y span", Device_Width)
        self.lum.set("z min", min_slabH)
        self.lum.set("z max", max_slabH)
        self.lum.set("x min", min_subL)
        self.lum.set("x max", max_subL)
        self.lum.set("material", MaterialSlab)

        z_Offset = max_slabH + WG_Height / 2
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
            
        # create_cover
        self.lum.addrect()
        self.lum.set("name", "cladding")
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





    def ArcWaveguide(self, Parameters):
        '''


        Parameters
        ----------
        Parameters : Dictionary
            Dictionary with all the data needet for the Bend Wavaguide. Data needet:
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
                Bending angle of the Waveguide. Set it to 0 to simulate the straight waveguide. 
                If Waveguide Angle is different then 0, then the straight waveguide will be tilted 
                at the choosen degrees. 
            Parameters["Wavelength"] : int/float
                Wavelength
            Parameters["S_Band Radius"] : int/float
                Radius of the Circle in um
            Parameters['Arc deg'] : int
                Can be 90 or 180 for 1/4 of a cirle or 1/2 of a circle.

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

            # Device specifications
            # Device_Width = 2 * radius + WaveLength * 2 + WG_Width * 2  # MMI_Width

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

            # creating the thin film
            min_slabH = max_subH
            max_slabH = max_subH + Slab_Height

            self.lum.addrect()
            self.lum.set("name", "LN_slab")
            self.lum.set("y", radius * m)
            self.lum.set("y span", m * radius * 2 + WG_Width)
            self.lum.set("z min", min_slabH)
            self.lum.set("z max", max_slabH)
            self.lum.set("x", m * radius)
            self.lum.set("x span", (m * radius * 2 + WG_Width))
            self.lum.set("material", MaterialSlab)

            z_Offset = max_slabH + WG_Height / 2
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
            pole = np.array([[radius * 0, radius * 1], [radius * m, radius * 1], [radius * 1, radius * m], [radius * 1, radius * 0]])
            self.lum.set("poles", pole)
            self.lum.set("material", MaterialWG)

            
            # # create_cover
            # self.lum.addrect()
            # self.lum.set("name", "cladding")
            # self.lum.set("material", MaterialClad)
            # self.lum.set("y", radius * m)
            # self.lum.set("y span", m * radius * 2 + WG_Width)
            # self.lum.set("z min", max_slabH)
            # self.lum.set("z max", max_slabH*2)
            # self.lum.set("x", m * radius)
            # self.lum.set("x span", (m * radius * 2 + WG_Width))
            # self.lum.set("override mesh order from material database", True)
            # self.lum.set("mesh order", 4)
            # self.lum.set("alpha", 0.7)




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

            # creating the thin film
            min_slabH = max_subH
            max_slabH = max_subH + Slab_Height

            self.lum.addrect()
            self.lum.set("name", "LN_slab")
            self.lum.set("y", radius * m)
            self.lum.set("y span", m * radius * 2 + WG_Width)
            self.lum.set("z min", min_slabH)
            self.lum.set("z max", max_slabH)
            self.lum.set("x", 0)
            self.lum.set("x span", (m * radius * 2 + WG_Width) * 2)
            self.lum.set("material", MaterialSlab)

            z_Offset = max_slabH + WG_Height / 2
            # Triangle EQ for waveguide Width
            x = abs(WG_Height / (np.cos((angle) * np.pi / 180)))  # in Radians
            extention = np.sqrt(x ** 2 - WG_Height ** 2)
            WG_W = WG_Width + 2 * extention

            names = ['Arc_Waveguide1', 'Arc_Waveguie2']

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
            
            
            # # create_cover
            # self.lum.addrect()
            # self.lum.set("name", "cladding")
            # self.lum.set("material", MaterialClad)
            # self.lum.set("y", radius * m)
            # self.lum.set("y span", m * radius * 2 + WG_Width)
            # self.lum.set("z min", max_slabH)
            # self.lum.set("z max", max_slabH*2)
            # self.lum.set("x", 0)
            # self.lum.set("x span", (m * radius * 2 + WG_Width)*2)
            # self.lum.set("override mesh order from material database", True)
            # self.lum.set("mesh order", 4)
            # self.lum.set("alpha", 0.7)

        else:
            raise ValueError(
                "Incorrect Arc Value. The arc can only be an integer and can only be arc = 90 or arc = 180!")
            




    def MMI2x2(self, Parameters):
        '''


        Parameters
        ----------
        Parameters : Dictionary
            Dictionary with all the data needet for the Bend Wavaguide. Data needet:
                Parameters
                ----------
                Parameters['Material'] : list of str
                    List of Materials. the list should be with names (str) of a valid Lumerical materials.
                    Check the names in Lumerical Materials viewer.
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
                 Parameters['Wavelength'] : int/float
                    Wavelength
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
        WaveLength = Parameters['Wavelength']
        TaperLength = Parameters['Taper Length']
        TaperWidth = Parameters['Taper Width']
        Taper = Parameters['Taper']



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
        Device_Width = MMI_Width + WaveLength * 2  # MMI_Width

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

        # creating the thin film
        min_slabH = max_subH
        max_slabH = max_subH + Slab_Height

        self.lum.addrect()
        self.lum.set("name", "LN_slab")
        self.lum.set("y min", min_subW)
        self.lum.set("y max", max_subW)
        self.lum.set("z min", min_slabH)
        self.lum.set("z max", max_slabH)
        self.lum.set("x min", min_subL)
        self.lum.set("x max", max_subL)
        self.lum.set("material", MaterialSlab)

        # creating the MMI
        max_MMIH = WG_Height
        max_MMIL = MMI_Length / 2
        min_MMIL = -MMI_Length / 2
        z_Offset = max_slabH + max_MMIH / 2

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
            
            elif offset_WG2 <0.5e-6:
                self.lum.deleteall()
                raise ValueError('The distance between the Tapers is less then 1 um !')
                


            else:
                # Mirror the In and Out WG on both sides
                maxWGL = [WG_Length, WG_Length, 0, 0]
                minWGL = [0, 0, -WG_Length, -WG_Length]
                xPos = [max_MMIL, max_MMIL, min_MMIL, min_MMIL]
                yPos = [WG_Width / 2 + posOffset / 2, -(WG_Width / 2 + posOffset / 2), WG_Width / 2 + posOffset / 2,
                        -(WG_Width / 2 + posOffset / 2)]

                # Names of the WGs
                names = ['LN_input_WG_L', 'LN_input_WG_R', 'LN_output_WG_L', 'LN_output_WG_R']

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




        elif Taper == True:

            # Delate the Structure to start new
            self.lum.deleteall()
            # Device specifications
            Device_Length = MMI_Length + 2 * WG_Length + 2 * TaperLength
            Device_Width = MMI_Width + WaveLength * 2  # MMI_Width

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

            # creating the thin film
            min_slabH = max_subH
            max_slabH = max_subH + Slab_Height

            self.lum.addrect()
            self.lum.set("name", "LN_slab")
            self.lum.set("y min", min_subW)
            self.lum.set("y max", max_subW)
            self.lum.set("z min", min_slabH)
            self.lum.set("z max", max_slabH)
            self.lum.set("x min", min_subL)
            self.lum.set("x max", max_subL)
            self.lum.set("material", MaterialSlab)

            # creating the MMI
            max_MMIH = WG_Height
            max_MMIL = MMI_Length / 2
            min_MMIL = -MMI_Length / 2
            z_Offset = max_slabH + max_MMIH / 2

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
            names = ['LN_input_WG_L', 'LN_input_WG_R', 'LN_output_WG_L', 'LN_output_WG_R']
            TapersNames = ['Taper_input_WG_L', 'Taper_input_WG_R', 'Taper_output_WG_L', 'Taper_output_WG_R']

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
            if BotCornerDistance < 1e-6:
                self.lum.deleteall()
                raise ValueError('The distance between the Tapers is less then 1 um !')
            elif offset_Taper > OffMax:
                self.lum.deleteall()
                raise ValueError('You are Trying to move the Taper outside the MMI. This is not possible!')



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
        else:
            raise ValueError(
                "Incorect Taper input. Taper must be an boolen. You can choose from Taper = True or Taper = False!")






    def MMI2x2_Trapez(self, Parameters):
        '''


        Parameters
        ----------
        Parameters : Dictionary
            Dictionary with all the data needet for the Bend Wavaguide. Data needet:
                Parameters
                ----------
                Parameters['Material'] : list of str
                    List of Materials. the list should be with names (str) of a valid Lumerical materials.
                    Check the names in Lumerical Materials viewer.
                 Parameters['Substrate Height'] : int/float
                    Substrate height.
                Parameters['MMI Width'] : int/float
                    Width of the MMI.
                Parameters['MMI Length'] : int/float
                    Length of the MMI.
                Parameters['Middle MMI Length'] : int/float
                    Length of the MMI in the Middle Part due to inperfection by manufacturing
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
                 Parameters['Wavelength'] : int/float
                    Wavelength
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
        MMI_Length2 = Parameters['Middle MMI Length']
        angle = Parameters['angle']
        WG_Height = Parameters['WG Height']
        WG_Width = Parameters['WG Width']
        WG_Length = Parameters['WG Length']
        posOffset = Parameters['Position Offset']
        Slab_Height = Parameters['Slab Height']
        WaveLength = Parameters['Wavelength']
        TaperLength = Parameters['Taper Length']
        TaperWidth = Parameters['Taper Width']
        Taper = Parameters['Taper']



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
        Device_Width = MMI_Width + WaveLength * 2  # MMI_Width

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

        # creating the thin film
        min_slabH = max_subH
        max_slabH = max_subH + Slab_Height

        self.lum.addrect()
        self.lum.set("name", "LN_slab")
        self.lum.set("y min", min_subW)
        self.lum.set("y max", max_subW)
        self.lum.set("z min", min_slabH)
        self.lum.set("z max", max_slabH)
        self.lum.set("x min", min_subL)
        self.lum.set("x max", max_subL)
        self.lum.set("material", MaterialSlab)

        # creating the MMI
        max_MMIH = WG_Height
        max_MMIL = MMI_Length / 2
        min_MMIL = -MMI_Length / 2
        z_Offset = max_slabH + max_MMIH / 2

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
        
        
        #Create the MMI L2 Tapers for the Trapezoid
        TaperNames = ["Input_Trapez", "Output_Trapez"]
        myscript = self.Script()
        
        
        self.lum.addstructuregroup()
        self.lum.set("name",TaperNames[0])
        self.lum.set("construction group",1)
        self.lum.adduserprop("thickness",2, WG_Height)
        self.lum.adduserprop("angle_side",0, angle)
        self.lum.adduserprop("width_l",2, WG_Width)
        self.lum.adduserprop("width_r",2, MMI_Wid)
        self.lum.adduserprop("hfrac_ref",0,1)
        self.lum.adduserprop("len",2, MMI_Length2)
        self.lum.adduserprop("material",5,MaterialWG)
        self.lum.adduserprop("index",0,1)
        self.lum.set("script",myscript) 
        self.lum.set("x", -MMI_Length / 2 - MMI_Length2/2) 
        self.lum.set("z", z_Offset)
        self.lum.set("y", 0)
        
        
        self.lum.addstructuregroup()
        self.lum.set("name",TaperNames[1])
        self.lum.set("construction group",1)
        self.lum.adduserprop("thickness",2, WG_Height)
        self.lum.adduserprop("angle_side",0, angle)
        self.lum.adduserprop("width_l",2, WG_Width)
        self.lum.adduserprop("width_r",2, MMI_Wid)
        self.lum.adduserprop("hfrac_ref",0,1)
        self.lum.adduserprop("len",2, MMI_Length2)
        self.lum.adduserprop("material",5,MaterialWG)
        self.lum.adduserprop("index",0,1)
        self.lum.set("script",myscript) 
        self.lum.set("first axis", "z")
        self.lum.set("rotation 1",180)
        self.lum.set("x", MMI_Length / 2 + MMI_Length2/2)
        self.lum.set("z", z_Offset)
        self.lum.set("y", 0)



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
            
            elif offset_WG2 <0.5e-6:
                self.lum.deleteall()
                raise ValueError('The distance between the Tapers is less then 1 um !')
                


            else:
                # Mirror the In and Out WG on both sides
                maxWGL = [WG_Length, WG_Length, 0, 0]
                minWGL = [0, 0, -WG_Length, -WG_Length]
                xPos = [max_MMIL, max_MMIL, min_MMIL, min_MMIL]
                yPos = [WG_Width / 2 + posOffset / 2, -(WG_Width / 2 + posOffset / 2), WG_Width / 2 + posOffset / 2,
                        -(WG_Width / 2 + posOffset / 2)]

                # Names of the WGs
                names = ['LN_input_WG_L', 'LN_input_WG_R', 'LN_output_WG_L', 'LN_output_WG_R']

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




        elif Taper == True:

            # Delate the Structure to start new
            self.lum.deleteall()
            # Device specifications
            Device_Length = MMI_Length + 2 * WG_Length + 2 * TaperLength
            Device_Width = MMI_Width + WaveLength * 2  # MMI_Width

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

            # creating the thin film
            min_slabH = max_subH
            max_slabH = max_subH + Slab_Height

            self.lum.addrect()
            self.lum.set("name", "LN_slab")
            self.lum.set("y min", min_subW)
            self.lum.set("y max", max_subW)
            self.lum.set("z min", min_slabH)
            self.lum.set("z max", max_slabH)
            self.lum.set("x min", min_subL)
            self.lum.set("x max", max_subL)
            self.lum.set("material", MaterialSlab)

            # creating the MMI
            max_MMIH = WG_Height
            max_MMIL = MMI_Length / 2
            min_MMIL = -MMI_Length / 2
            z_Offset = max_slabH + max_MMIH / 2

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
            
            
            #Create the MMI L2 Tapers for the Trapezoid
            TaperNames = ["Input_Trapez", "Output_Trapez"]
            myscript = self.Script()
            
            
            self.lum.addstructuregroup()
            self.lum.set("name",TaperNames[0])
            self.lum.set("construction group",1)
            self.lum.adduserprop("thickness",2, WG_Height)
            self.lum.adduserprop("angle_side",0, angle)
            self.lum.adduserprop("width_l",2, WG_Width)
            self.lum.adduserprop("width_r",2, MMI_Wid)
            self.lum.adduserprop("hfrac_ref",0,1)
            self.lum.adduserprop("len",2, MMI_Length2)
            self.lum.adduserprop("material",5,MaterialWG)
            self.lum.adduserprop("index",0,1)
            self.lum.set("script",myscript) 
            self.lum.set("x", -MMI_Length / 2 - MMI_Length2/2) 
            self.lum.set("z", z_Offset)
            self.lum.set("y", 0)
            
            
            self.lum.addstructuregroup()
            self.lum.set("name",TaperNames[1])
            self.lum.set("construction group",1)
            self.lum.adduserprop("thickness",2, WG_Height)
            self.lum.adduserprop("angle_side",0, angle)
            self.lum.adduserprop("width_l",2, WG_Width)
            self.lum.adduserprop("width_r",2, MMI_Wid)
            self.lum.adduserprop("hfrac_ref",0,1)
            self.lum.adduserprop("len",2, MMI_Length2)
            self.lum.adduserprop("material",5,MaterialWG)
            self.lum.adduserprop("index",0,1)
            self.lum.set("script",myscript) 
            self.lum.set("first axis", "z")
            self.lum.set("rotation 1",180)
            self.lum.set("x", MMI_Length / 2 + MMI_Length2/2)
            self.lum.set("z", z_Offset)
            self.lum.set("y", 0)

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
            names = ['LN_input_WG_L', 'LN_input_WG_R', 'LN_output_WG_L', 'LN_output_WG_R']
            TapersNames = ['Taper_input_WG_L', 'Taper_input_WG_R', 'Taper_output_WG_L', 'Taper_output_WG_R']

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
            if BotCornerDistance < 1e-6:
                self.lum.deleteall()
                raise ValueError('The distance between the Tapers is less then 1 um !')
            elif offset_Taper > OffMax:
                self.lum.deleteall()
                raise ValueError('You are Trying to move the Taper outside the MMI. This is not possible!')



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
        else:
            raise ValueError(
                "Incorect Taper input. Taper must be an boolen. You can choose from Taper = True or Taper = False!")





    def MMI2x1(self, Parameters):
        '''


        Parameters
        ----------
        Parameters : Dictionary
            Dictionary with all the data needet for the Bend Wavaguide. Data needet:
                Parameters
                ----------
                Parameters['Material'] : list of str
                    List of Materials. the list should be with names (str) of a valid Lumerical materials.
                    Check the names in Lumerical Materials viewer.
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
                Parameters['Wavelength'] : int/float
                    Wavelength
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
        WaveLength = Parameters['Wavelength']
        TaperLength = Parameters['Taper Length']
        TaperWidth = Parameters['Taper Width']
        Taper = Parameters['Taper']



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
        Device_Width = MMI_Width + WaveLength * 2  # MMI_Width

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

        # creating the thin film
        min_slabH = max_subH
        max_slabH = max_subH + Slab_Height

        self.lum.addrect()
        self.lum.set("name", "LN_slab")
        self.lum.set("y min", min_subW)
        self.lum.set("y max", max_subW)
        self.lum.set("z min", min_slabH)
        self.lum.set("z max", max_slabH)
        self.lum.set("x min", min_subL)
        self.lum.set("x max", max_subL)
        self.lum.set("material", MaterialSlab)

        # creating the MMI
        max_MMIH = WG_Height
        max_MMIL = MMI_Length / 2
        min_MMIL = -MMI_Length / 2
        z_Offset = max_slabH + max_MMIH / 2

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
        #     self.lum.deleteall()
        #     raise ValueError('The distance between the Tapers is less then 1 um !')
        # else:

        if Taper == False:

            # if offset_WG > OffMax:
            #     self.lum.deleteall()
            #     raise ValueError('You are Trying to move the Waveguide outside the MMI. This is not possible!')
            # else:
            # Mirror the In and Out WG on both sides
            maxWGL = [WG_Length, 0, 0]
            minWGL = [0, -WG_Length, -WG_Length]
            xPos = [max_MMIL, min_MMIL, min_MMIL]
            yPos = [0 + OffsetInput, WG_Width / 2 + posOffset / 2, - WG_Width / 2 - posOffset / 2]

            # Names of the WGs
            names = ['LN_input_WG', 'LN_output_WG_L', 'LN_output_WG_R']

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




        elif Taper == True:
            # if offset_Taper > OffMax:
            #     self.lum.deleteall()
            #     raise ValueError('You are Trying to move the Taper outside the MMI. This is not possible!')
            # else:
            # Delate the Structure to start new
            self.lum.deleteall()
            # Device specifications
            Device_Length = MMI_Length + 2 * WG_Length + 2*TaperLength
            Device_Width = MMI_Width + WaveLength * 2  # MMI_Width

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

            # creating the thin film
            min_slabH = max_subH
            max_slabH = max_subH + Slab_Height

            self.lum.addrect()
            self.lum.set("name", "LN_slab")
            self.lum.set("y min", min_subW)
            self.lum.set("y max", max_subW)
            self.lum.set("z min", min_slabH)
            self.lum.set("z max", max_slabH)
            self.lum.set("x min", min_subL)
            self.lum.set("x max", max_subL)
            self.lum.set("material", MaterialSlab)

            # creating the MMI
            max_MMIH = WG_Height
            max_MMIL = MMI_Length / 2
            min_MMIL = -MMI_Length / 2
            z_Offset = max_slabH + max_MMIH / 2

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
            names = ['LN_input_WG', 'LN_output_WG_L', 'LN_output_WG_R']
            TapersNames = ['Taper_input_WG', 'Taper_output_WG_L', 'Taper_output_WG_R']

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

        else:
            raise ValueError(
                "Incorect Taper input. Taper must be an boolen. You can choose from Taper = True or Taper = False!")

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
    
    
    
    
    


    def MMI2x1_Trapez(self, Parameters):
        '''


        Parameters
        ----------
        Parameters : Dictionary
            Dictionary with all the data needet for the Bend Wavaguide. Data needet:
                Parameters
                ----------
                Parameters['Material'] : list of str
                    List of Materials. the list should be with names (str) of a valid Lumerical materials.
                    Check the names in Lumerical Materials viewer.
                Parameters['Substrate Height'] : int/float
                    Substrate height.
                Parameters['MMI Width'] : int/float
                    Width of the MMI.
                Parameters['MMI Length'] : int/float
                    Length of the MMI.
                Parameters['Middle MMI Length']: int/float
                    Length of the MMI in the Middle Part due to inperfection by manufacturing
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
                Parameters['Wavelength'] : int/float
                    Wavelength
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
        MMI_Length2 = Parameters['Middle MMI Length']
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
        Device_Width = MMI_Width + WaveLength * 2  # MMI_Width

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

        # creating the thin film
        min_slabH = max_subH
        max_slabH = max_subH + Slab_Height

        self.lum.addrect()
        self.lum.set("name", "LN_slab")
        self.lum.set("y min", min_subW)
        self.lum.set("y max", max_subW)
        self.lum.set("z min", min_slabH)
        self.lum.set("z max", max_slabH)
        self.lum.set("x min", min_subL)
        self.lum.set("x max", max_subL)
        self.lum.set("material", MaterialSlab)

        # creating the MMI
        max_MMIH = WG_Height
        max_MMIL = MMI_Length / 2
        min_MMIL = -MMI_Length / 2
        z_Offset = max_slabH + max_MMIH / 2

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
        
        
        
        #Create the MMI L2 Tapers for the Trapezoid
        TaperNames = ["Input_Trapez", "Output_Trapez"]
        myscript = self.Script()
        
        
        self.lum.addstructuregroup()
        self.lum.set("name",TaperNames[0])
        self.lum.set("construction group",1)
        self.lum.adduserprop("thickness",2, WG_Height)
        self.lum.adduserprop("angle_side",0, angle)
        self.lum.adduserprop("width_l",2, WG_Width)
        self.lum.adduserprop("width_r",2, MMI_Wid)
        self.lum.adduserprop("hfrac_ref",0,1)
        self.lum.adduserprop("len",2, MMI_Length2)
        self.lum.adduserprop("material",5,MaterialWG)
        self.lum.adduserprop("index",0,1)
        self.lum.set("script",myscript) 
        self.lum.set("x", -MMI_Length / 2 - MMI_Length2/2) 
        self.lum.set("z", z_Offset)
        self.lum.set("y", 0)
        
        
        self.lum.addstructuregroup()
        self.lum.set("name",TaperNames[1])
        self.lum.set("construction group",1)
        self.lum.adduserprop("thickness",2, WG_Height)
        self.lum.adduserprop("angle_side",0, angle)
        self.lum.adduserprop("width_l",2, WG_Width)
        self.lum.adduserprop("width_r",2, MMI_Wid)
        self.lum.adduserprop("hfrac_ref",0,1)
        self.lum.adduserprop("len",2, MMI_Length2)
        self.lum.adduserprop("material",5,MaterialWG)
        self.lum.adduserprop("index",0,1)
        self.lum.set("script",myscript) 
        self.lum.set("first axis", "z")
        self.lum.set("rotation 1",180)
        self.lum.set("x", MMI_Length / 2 + MMI_Length2/2)
        self.lum.set("z", z_Offset)
        self.lum.set("y", 0)
        
        
        
        

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

                    # Names of the WGs
                    names = ['LN_input_WG', 'LN_output_WG_L', 'LN_output_WG_R']

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




            elif Taper == True:
                if offset_Taper > OffMax:
                    self.lum.deleteall()
                    raise ValueError('You are Trying to move the Taper outside the MMI. This is not possible!')
                else:
                    # Delate the Structure to start new
                    self.lum.deleteall()
                    # Device specifications
                    Device_Length = MMI_Length + 2 * WG_Length + 2*TaperLength
                    Device_Width = MMI_Width + WaveLength * 2  # MMI_Width

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

                    # creating the thin film
                    min_slabH = max_subH
                    max_slabH = max_subH + Slab_Height

                    self.lum.addrect()
                    self.lum.set("name", "LN_slab")
                    self.lum.set("y min", min_subW)
                    self.lum.set("y max", max_subW)
                    self.lum.set("z min", min_slabH)
                    self.lum.set("z max", max_slabH)
                    self.lum.set("x min", min_subL)
                    self.lum.set("x max", max_subL)
                    self.lum.set("material", MaterialSlab)

                    # creating the MMI
                    max_MMIH = WG_Height
                    max_MMIL = MMI_Length / 2
                    min_MMIL = -MMI_Length / 2
                    z_Offset = max_slabH + max_MMIH / 2

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
                    
                    #Create the MMI L2 Tapers for the Trapezoid
                    TaperNames = ["Input_Trapez", "Output_Trapez"]
                    
                    self.lum.addstructuregroup()
                    self.lum.set("name",TaperNames[0])
                    self.lum.set("construction group",1)
                    self.lum.adduserprop("thickness",2, WG_Height)
                    self.lum.adduserprop("angle_side",0, angle)
                    self.lum.adduserprop("width_l",2, WG_Width)
                    self.lum.adduserprop("width_r",2, MMI_Wid)
                    self.lum.adduserprop("hfrac_ref",0,1)
                    self.lum.adduserprop("len",2, MMI_Length2)
                    self.lum.adduserprop("material",5,MaterialWG)
                    self.lum.adduserprop("index",0,1)
                    self.lum.set("script",myscript) 
                    # self.lum.set("first axis", "z")
                    # self.lum.set("rotation 1",angleTheta)
                    self.lum.set("x", -MMI_Length / 2 - MMI_Length2/2) #-MMI_Length / 2
                    self.lum.set("z", z_Offset)
                    self.lum.set("y", 0)
                    
                    
                    self.lum.addstructuregroup()
                    self.lum.set("name",TaperNames[1])
                    self.lum.set("construction group",1)
                    self.lum.adduserprop("thickness",2, WG_Height)
                    self.lum.adduserprop("angle_side",0, angle)
                    self.lum.adduserprop("width_l",2, WG_Width)
                    self.lum.adduserprop("width_r",2, MMI_Wid)
                    self.lum.adduserprop("hfrac_ref",0,1)
                    self.lum.adduserprop("len",2, MMI_Length2)
                    self.lum.adduserprop("material",5,MaterialWG)
                    self.lum.adduserprop("index",0,1)
                    self.lum.set("script",myscript) 
                    self.lum.set("first axis", "z")
                    self.lum.set("rotation 1",180)
                    self.lum.set("x", MMI_Length / 2 + MMI_Length2/2) #-MMI_Length / 2
                    self.lum.set("z", z_Offset)
                    self.lum.set("y", 0)


                    # New x Length of the Tapers
                    maxLength = max_MMIL + TaperLength
                    minLength = min_MMIL - TaperLength

                    # Mirror the In and Out WG on both sides
                    maxWGL = [WG_Length, 0, 0]
                    minWGL = [0, -WG_Length, -WG_Length]
                    xPos = [maxLength, minLength, minLength]
                    yPos = [0 + OffsetInput, WG_Width / 2 + posOffset / 2, - WG_Width / 2 - posOffset / 2]

                    # Names of the WGs
                    names = ['LN_input_WG', 'LN_output_WG_L', 'LN_output_WG_R']
                    TapersNames = ['Taper_input_WG', 'Taper_output_WG_L', 'Taper_output_WG_R']

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

            else:
                raise ValueError(
                    "Incorect Taper input. Taper must be an boolen. You can choose from Taper = True or Taper = False!")

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




    def DirectionalCoupler(self, Parameters):
        '''


        Parameters
        ----------
        Parameters : Dictionary
            Dictionary with all the data needet for the Bend Wavaguide. Data needet:
                Parameters
                ----------
                Parameters['Material'] : list of str
                    List of Materials. the list should be with names (str) of a valid Lumerical materials.
                    Check the names in Lumerical Materials viewer.
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

        # creating the thin film
        min_slabH = max_subH
        max_slabH = max_subH + Slab_Height

        self.lum.addrect()
        self.lum.set("name", "LN_slab")
        self.lum.set("y min", min_subW)
        self.lum.set("y max", max_subW)
        self.lum.set("z min", min_slabH)
        self.lum.set("z max", max_slabH)
        self.lum.set("x min", min_subL)
        self.lum.set("x max", max_subL)
        self.lum.set("material", MaterialSlab)

        # Positions of the Input and Output WGs
        # Triangle EQ for MMI Width
        x = abs(WG_Height / (np.cos((angle) * np.pi / 180)))  # in Radians
        extention = np.sqrt(x ** 2 - WG_Height ** 2)
        WG_W = WG_Width + 2 * extention
        WG_Width_bott = WG_W
        z_Offset = max_slabH + WG_Height / 2

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
            names = ['Top_WG', 'Bottom_WG']

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





    def WDM(self, Parameters):
        '''
        

        Parameters
        ----------
        Parameters : Dictionary
            Parameters['Material'] : list of str
                List of Materials. the list should be with names (str) of a valid Lumerical materials.
                Check the names in Lumerical Materials viewer.
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
            Parameters['Slab Height'] : int/float
                Height of the slab.
            Parameters['Wavelength'] : int/float
                Wavelength
            Parameters['Angle Thetha'] : 
                Input and output angle of the waveguide. This is only temporally 
            Parameters['Taper Width'] : int/float
                Backside width of the Taper, frontside width is the waveguide width
            Parameters['Taper Length'] : int/float
                Length of the Taper. 
            Parameters['Taper'] : boolen
                If Taper == False, no Taper will be placed with the Waveguides.
                If Taper == True, tapers will be placed with the waveguides
        
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
        Slab_Height = Parameters['Slab Height']
        WaveLength = Parameters['Wavelength']
        angleTheta = Parameters['Angle Thetha']
        TaperLength = Parameters['Taper Length']
        TaperWidth = Parameters['Taper Width']
        Taper = Parameters['Taper']

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
            Device_Width = MMI_Width + 2*WG_Length + WaveLength * 2  # MMI_Width

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
             
            # creating the thin film
            min_slabH = max_subH
            max_slabH = max_subH + Slab_Height
             
            self.lum.addrect()
            self.lum.set("name", "LN_slab")
            self.lum.set("z min", min_slabH)
            self.lum.set("z max", max_slabH)
            self.lum.set("x", 0)
            self.lum.set("x span", Device_Length)
            self.lum.set("y", 0)
            self.lum.set("y span", Device_Width)
            self.lum.set("material", MaterialSlab)
             
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



        elif Taper == True:
            
            
            # Device specifications
            Device_Length = MMI_Length + 4 * WG_Length
            Device_Width = MMI_Width + 2*WG_Length + WaveLength * 2  # MMI_Width

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
             
            
             
            self.lum.addrect()
            self.lum.set("name", "LN_slab")
            self.lum.set("z min", min_slabH)
            self.lum.set("z max", max_slabH)
            self.lum.set("x", 0)
            self.lum.set("x span", Device_Length)
            self.lum.set("y", 0)
            self.lum.set("y span", Device_Width)
            self.lum.set("material", MaterialSlab)
           
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
                
                
                
        
        
   
    def InverseTaper(self, Parameters):
        '''
        

        Parameters
        ----------
        Parameters['Material'] : list of str
            List of Materials. the list should be with names (str) of a valid Lumerical materials.
            Check the names in Lumerical Materials viewer.
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
        Slab_Height = Parameters['Slab Height']
        TaperLength = Parameters['Taper Length']
        TaperWidth = Parameters['Taper Width']
        TaperWidthF = Parameters['PWB Taper Width Front']
        TaperWidthB = Parameters['PWB Taper Width Back']
        TaperHightB = Parameters['PWB Taper Hight Back']
        TaperHightF = Parameters['PWB Taper Hight Front']
        TaperLength_PWB = Parameters['PWB Taper Length']


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
        max_subH = Substrate_Height
        min_subH = -Substrate_Height
        
        # make substrate
        self.lum.addrect()
        self.lum.set("name", "Substrate")
        self.lum.set("y", 0)
        self.lum.set("y span", TaperWidthB * 2)
        self.lum.set("z min", min_subH)
        self.lum.set("z max", max_subH)
        self.lum.set("x min", -TaperLength_PWB / 2 - 0.1e-6)
        self.lum.set("x max", TaperLength_PWB )
        # self.lum.set("x", 10e-6)
        # self.lum.set("x span", TaperLength_PWB * 2)
        self.lum.set("material", MaterialSub)
        
        # creating the thin film
        min_slabH = max_subH
        max_slabH = max_subH + Slab_Height
        
        # Slab
        self.lum.addrect()
        self.lum.set("name", "LN_slab")
        self.lum.set("y", 0e-12)
        self.lum.set("y span", TaperWidthB * 2)
        self.lum.set("z min", min_slabH)
        self.lum.set("z max", max_slabH)
        self.lum.set("x min", -TaperLength_PWB / 2 - 0.1e-6)
        self.lum.set("x max", TaperLength_PWB)
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
        

      
        
        z_Offset = max_slabH
        
        # PWD Taper Hights
        TaperZmin = z_Offset
        TaperZmaxF = z_Offset + TaperHightF
        TaperZmaxB = z_Offset + TaperHightB
        
        # Inverse Taper Hights
        TaperZmin = z_Offset
        TaperZmax = z_Offset + WG_Height
        
        # PWB Taper Length
        PWB_TaperXmin = -TaperLength_PWB / 2
        PWB_TaperXmax = TaperLength_PWB / 2
        
        # Inverse Taper Length
        TaperXmin = -TaperLength / 2
        TaperXmax = TaperLength / 2
        
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
        
        
        # Make Sqered WG-Extention of the PWB for mode Calculations
        
        # Extra Waveguide Lenght
        Ext_WGLength = 1e-6
        PWD_x_Offset = PWB_TaperXmin
        
        self.lum.addrect()
        self.lum.set('x min', PWD_x_Offset - Ext_WGLength)
        self.lum.set('x max',PWD_x_Offset)
        self.lum.set('y', 0)
        self.lum.set('y span', TaperWidthB)
        self.lum.set('z min', TaperZmin)
        self.lum.set('z max', TaperZmaxB)
        self.lum.set("override mesh order from material database", 1)
        self.lum.set("mesh order", 3)
        self.lum.set('name', 'WG_Extention_PWB')
        self.lum.set('material', MaterialPWB)
        
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
        pole = np.array([[0, 0], [x_max, 0]])
        self.lum.set("poles", pole)
        self.lum.set("material", MaterialSlab)







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
        self.lum.set("name", "Substrate_Main")
        self.lum.set("y", 0)
        self.lum.set("y span", Device_Width)
        self.lum.set("z min", min_subH)
        self.lum.set("z max", max_subH)
        self.lum.set("x min", min_subL)
        self.lum.set("x max", max_subL)
        self.lum.set("material", MaterialSub)

        # creating the thin film
        min_slabH = max_subH
        max_slabH = max_subH + Slab_Height

        self.lum.addrect()
        self.lum.set("name", "LN_slab_Main")
        self.lum.set("y", 0)
        self.lum.set("y span", Device_Width)
        self.lum.set("z min", min_slabH)
        self.lum.set("z max", max_slabH)
        self.lum.set("x min", min_subL)
        self.lum.set("x max", max_subL)
        self.lum.set("material", MaterialSlab)
        


        # creating the MMI
        max_MMIH = WG_Height
        max_MMIL = MMI_Length / 2
        min_MMIL = -MMI_Length / 2
        z_Offset = max_slabH + max_MMIH / 2

        # Triangle EQ for MMI Width
        x = abs(max_MMIH / (np.cos((angle) * np.pi / 180)))  # in Radians
        extention = np.sqrt(x ** 2 - max_MMIH ** 2)
        MMI_Wid = MMI_Width + 2 * extention
        
        Parameters__Position_X = [0, -SpaceX -MMI_Length  - 2*WG_Length , -SpaceX -MMI_Length  - 2*WG_Length]
        Parameters__Position_Y = [0, (MMI_Wid + SpaceY)/2,  -(MMI_Wid + SpaceY)/2]
        

        names_MMI = ["MMI_In", "MMI_Out1", "MMI_Out2"]
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
        
        
                    names = ["S-Bend_Top", "S-Bend_Bot"]
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
                    names1 = ['MMI_In_WG', 'MMI_Out_WG_Top', 'MMI_Out_WG_Bot']
                    names2 = ['MMIOut1_In_WG', 'MMIOut1_Out_WG_Top', 'MMIOut1_Out_WG_Bot']
                    names3 = ['MMIOut2_In_WG', 'MMIOut2_Out_WG_Top', 'MMIOut2_Out_WG_Bot']
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

                    # creating the thin film
                    min_slabH = max_subH
                    max_slabH = max_subH + Slab_Height

                    self.lum.addrect()
                    self.lum.set("name", "LN_slab_Mian")
                    self.lum.set("y", 0)
                    self.lum.set("y span", Device_Width)
                    self.lum.set("z min", min_slabH)
                    self.lum.set("z max", max_slabH)
                    self.lum.set("x min", min_subL)
                    self.lum.set("x max", max_subL)
                    self.lum.set("material", MaterialSlab)
                    

                    # creating the MMI
                    max_MMIH = WG_Height
                    max_MMIL = MMI_Length / 2
                    min_MMIL = -MMI_Length / 2
                    z_Offset = max_slabH + max_MMIH / 2

                    # Triangle EQ for MMI Width
                    x = abs(max_MMIH / (np.cos((angle) * np.pi / 180)))  # in Radians
                    extention = np.sqrt(x ** 2 - max_MMIH ** 2)
                    MMI_Wid = MMI_Width + 2 * extention
                    
                    
                    Parameters__Position_X =  [0, -SpaceX -MMI_Length -2*TaperLength - 2*WG_Length , -SpaceX -MMI_Length -2*TaperLength - 2*WG_Length]
                    Parameters__Position_Y = [0, (MMI_Wid + SpaceY)/2,  -(MMI_Wid + SpaceY)/2]
                    
                    names_MMI = ["MMI_In", "MMI_Out1", "MMI_Out2"]
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
        
        
                    names = ["S-Bend_Top", "S-Bend_Bot"]
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
                    
                    TapersNames1 = ['MMI_Taper_In_WG', 'MMI_Taper_Out_WG_Top', 'MMI_Taper_Out_WG_Bot']
                    TapersNames2 = ['MMIOut1_Taper_In_WG', 'MMIOut1_Taper_Out_WG_Top', 'MMIOut1_Taper_Out_WG_Bot']
                    TapersNames3 = ['MMIOut2_Taper_In_WG', 'MMIOut2_Taper_Out_WG_Top', 'MMIOut2_Taper_Out_WG_Bot']
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
                    names1 = ['MMI_In_WG', 'MMI_Out_WG_Top', 'MMI_Out_WG_Bot']
                    names2 = ['MMIOut1_In_WG', 'MMIOut1_Out_WG_Top', 'MMIOut1_Out_WG_Bot']
                    names3 = ['MMIOut2_In_WG', 'MMIOut2_Out_WG_Top', 'MMIOut2_In_WG_Bot']
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
        
        
        
        
       
        
        
     
            

 
                

    # =============================================================================
    # Functions for the solvers
    # =============================================================================


    def setStraightWaveguideFDTDSolver(self, Parameters):
        '''


        Parameters
        ----------
        Parameters : Dictionary
            Dictionary with all the data needet for the Bend Wavaguide. Data needet:
               Parameters['Substrate Height'] : int/float
                   Substrate height.
               Parameters['WG Lrngth'] : int/float
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




        # Device specifications
        Device_Width = 2*WG_Length + WaveLength * 2  # MMI_Width

        max_slabH = Slab_Height
        MonitorHeight = Substrate_Height + (max_slabH + WG_Height) / 2
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
                self.lum.set('z min bc', 'PML')
                self.lum.set('z max bc', 'PML')
                self.lum.set('mesh type', 'auto non-uniform')
                self.lum.set('min mesh step', x_res)
                self.lum.set('set simulation bandwidth', 0)
                self.lum.set('global source center wavelength', WaveLength)
                self.lum.set('global source wavelength span', 0)
    
    
    
                # Define Ports
                x = [0,EME_WGLength]
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
                    self.lum.addtime()
                    self.lum.set('name', name[i])
                    self.lum.set("x", x[i] )
                    self.lum.set("y", yPos[i])
                    self.lum.set("z", MonitorHeight)
                    self.lum.set('output Px', 1)
                    self.lum.set('output Py', 1)
                    self.lum.set('output Pz', 1)
    
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
                    self.lum.set("x", x[i] )
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
            Device_Width = 2*TaperLength + WaveLength * 2  # MMI_Width

            max_slabH = Slab_Height
            MonitorHeight = Substrate_Height + (max_slabH + WG_Height) / 2
            EME_WGLength = TaperLength * np.cos(angle * np.pi / 180)
            
            # Adds a FDTD Solver
            self.lum.addfdtd()
            self.lum.set("x", 0)
            self.lum.set("x span", TaperLength)
            self.lum.set("y", TaperLength/4)
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
            Diff_Span = y_Port_Span - WG_Width 
            x = [-TaperLength/2+ 0.1e-6 , TaperLength/2 - 0.1e-6]
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
                self.lum.set("x", x[i] )
                self.lum.set("y", yPos[i])
                self.lum.set("z", MonitorHeight)
                self.lum.set('output Px', 1)
                self.lum.set('output Py', 1)
                self.lum.set('output Pz', 1)
            
            
            # Add Movie monitor
            self.lum.addmovie()
            self.lum.set("y", TaperLength/4)
            self.lum.set("y span", Device_Width)
            self.lum.set("z", MonitorHeight)
            self.lum.set("x", 0)
            self.lum.set("x span", TaperLength)

            # Add Power and Freq Monitor
            self.lum.addpower()
            self.lum.set('monitor type', '2D Z-normal')
            self.lum.set("y", TaperLength/4)
            self.lum.set("y span", Device_Width)
            self.lum.set("x", 0)
            self.lum.set("x span", TaperLength)
            self.lum.set("z", MonitorHeight)
            self.lum.set('output Px', 1)
            self.lum.set('output Py', 1)
            self.lum.set('output Pz', 1)
            self.lum.set('output power', 1)





    def setArcWaveguideFDTDSolver(self, Parameters):
        '''


        Parameters
        ----------
        Parameters : Dictionary
            Dictionary with all the data needet for the Bend Wavaguide. Data needet:
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
            self.lum.set("x span", (m * radius * 2 + WG_Width))
            self.lum.set("y", radius * m)
            self.lum.set("y span", m * radius * 2 + WG_Width)
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

            x = [0, radius]
            direction = ['Forward', 'Forward']
            name = ['Input', 'Output']
            y = [radius, 0]

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
            self.lum.set("bent waveguide", 1)
            self.lum.set("bend radius", radius)
            


            # Power Monitor Port 1
            self.lum.addtime()
            self.lum.set('name', name[0])
            self.lum.set("x", x[0] + 0.2e-6)
            self.lum.set("y", y[0])
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
            self.lum.set("z span", 4e-6)
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
            Dictionary with all the data needet for the Bend Wavaguide. Data needet:
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
               Parameters["x span"] : int/float
                   Length of the S-Bend Waveguide
               Parameters["y span"] : int/float
                   Width of the S-Bend Waveguide
               Parameters["Mode"] : str
                    Mode to choose from ("fundamental TE mode", "fundamental TM mode", "fundamental mode")
               Parameters["Port Span"][0] : int/float
                   Span of Port in x direction
               Parameters["Port Span"][1] : int/float
                   Span of Port in y direction
               Parameters["Port Span"][2] ; int/float
                   Span of Port in z direction
                   

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

        # Adds a Eigenmode Expansion (EME) solver region to the MODE simulation environment.
        self.lum.addfdtd()
        self.lum.set("x", Device_Length / 2)
        self.lum.set("x span", Device_Length)
        self.lum.set("y", y_span / 2)
        self.lum.set("y span", Device_Width)
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
            Dictionary with all the data needet for the Bend Wavaguide. Data needet:
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
                Parameters['x res'] : int/float
                    Mesh cell sizes.
                Parameters['Slab Height'] : int/float
                    Height of the slab.
                Parameters['Wavelength'] : int/float
                    Wavelength.
                Parameters["Mode"] : str
                    Mode to choose from ("fundamental TE mode", "fundamental TM mode", "fundamental mode")



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
            Ports_mid = Substrate_Height + max_slabH
    
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
            Ports_mid = Substrate_Height + max_slabH
    
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
            Dictionary with all the data needet for the Bend Wavaguide. Data needet:
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
                Parameters['x res'] : int/float
                    Mesh cell sizes.
                Parameters['Slab Height'] : int/float
                    Height of the slab.
                Parameters['Wavelength'] : int/float
                    Wavelength.
                Parameters["Mode"] : str
                    Mode to choose from ("fundamental TE mode", "fundamental TM mode", "fundamental mode")
    
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
        
        
        if Taper == False:
            
    
    
            # Device specifications
            Device_Length = MMI_Length + 2 * WG_Length
            Device_Width = MMI_Width + WaveLength * 2  # MMI_Width
            max_slabH = Slab_Height
            # Ports_mid = Substrate_Height + (max_slabH + WG_Height) / 2
            # Ports_mid = Substrate_Height + Slab_Height + WG_Height / 2
            # MonitorHeight = Substrate_Height + (max_slabH + WG_Height) / 2
            Ports_mid = Substrate_Height + max_slabH
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
            
            yPort_vec = [OffsetInput, -(posOffset / 2 + WG_Width / 2), posOffset / 2 + WG_Width / 2]
            
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
            Ports_mid = Substrate_Height + max_slabH
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
            Dictionary with all the data needet for the Bend Wavaguide. Data needet:
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
        Ports_mid = Substrate_Height + (max_slabH + WG_Height) / 2

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





    def setStraightWaveguideEMESolver(self, Parameters):
        '''


        Parameters
        ----------
        Parameters : Dictionary
               Parameters['Substrate Height'] : int/float
                   Substrate height.
               Parameters['WG Height'] : int/float
                   Waveguide hight. Also the height of the MMI section
               Parameters['WG Length'] : int/float
                   Waveguide length
               Parameters['WG Length'] : int/float
                   Waveguide width
               Parameters["Taper"] : boolen
                   If Taper == False, only straight Waveguide will be simulated, 
                   If Taper == True an Taper will be simulated 
               Parameters['Taper Width'] : int/float
                   Taper backside Width. Taper Fronside width is the width of the Waveguide
               Parameters['Taper Length'] : int/float
                   Taper Length
               Parameters['y res']: int/float
                     EME Mesh resolutio,
               Parameters['y res']: int/float
                     EME Mesh resolutio,
               Parameters['Slab Height'] : int/float
                   Height of the slab.
               Parameters['Wavelength'] : int/float
                   Wavelength
               Parameters["Waveguide Angle"] : int/float
                   This Parameter will set the theta ratation angle of the port.
               Parameters["Port Span"] : list of floats/ints
                   List of x,y and z span of the Ports. For this simulation only y and z parametes will be taken.
               Parameters["Mode"] : str
                    Mode to choose from ("fundamental TE mode", "fundamental TM mode", "fundamental mode")




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
        MonitorHeight = Substrate_Height + max_slabH 
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
                    self.lum.set("z", max_slabH)
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
            MonitorHeight = Substrate_Height + (max_slabH + WG_Height) / 2
            Ports_mid = (max_slabH + WG_Height) / 2
            EME_WGLength = TaperLength * np.cos(angle * np.pi / 180)
            
            
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
            
            






    def setDCEMESolver(self, Parameters):

        '''

        Parameters
        ----------
        Parameters : Dictionary
            Dictionary with all the data needet for the Bend Wavaguide. Data needet:
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
                Parameters["Mode"] : str
                    Mode to choose from ("fundamental TE mode", "fundamental TM mode", "fundamental mode")
                Parameters["Port Span"] : list of floats/ints
                    List of x,y and z span of the Ports. For this simulation only y and z parametes will be taken.


        Raises
        ------
        ValueError
            DESCRIPTION.

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
        MonitorHeight = Substrate_Height + (max_slabH + WG_Height) / 2
        Ports_mid = (max_slabH + WG_Height) / 2

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
            Dictionary with all the data needet for the Bend Wavaguide. Data needet:
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
                Parameters["Port Span"] : list of floats/ints
                    List of x,y and z span of the Ports. For this simulation only y and z parametes will be taken.


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
            MonitorHeight = Substrate_Height + max_slabH
            # Ports_mid = (max_slabH + WG_Height) / 2
            Ports_mid = max_slabH

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
            MonitorHeight = Substrate_Height + max_slabH
            # Ports_mid = (max_slabH + WG_Height) / 2
            Ports_mid = max_slabH

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
            Dictionary with all the data needet for the Bend Wavaguide. Data needet:
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
                Parameters['y res'] : int/float
                    Mesh cell sizes.
                Parameters['z res'] : int/float
                    Mesh cell size.
                Parameters['Slab Height'] : int/float
                    Height of the slab.
                Parameters['Wavelength'] : int/float
                    Wavelength.
                Parameters["Mode"] : str
                    Mode to choose from ("fundamental TE mode", "fundamental TM mode", "fundamental mode")
                Parameters["Port Span"] : list of floats/ints
                    List of x,y and z span of the Ports. For this simulation only y and z parametes will be taken.



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
        WG_Width = Parameters['WG Width']
        WG_Length = Parameters['WG Length']
        posOffset = Parameters['Position Offset']
        # x_res            = Parameters['x res']
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
            MonitorHeight = Substrate_Height + max_slabH
            # Ports_mid = (max_slabH + WG_Height) / 2
            Ports_mid = max_slabH

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
            MonitorHeight = Substrate_Height + max_slabH
            # Ports_mid = (max_slabH + WG_Height) / 2
            Ports_mid = max_slabH

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
         Parameters["Port Span"] : list of floats/ints
            List of x,y and z span of the Ports. For this simulation only y and z parametes will be taken.



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
        MonitorHeight = Substrate_Height + max_slabH
        # Ports_mid = (max_slabH + WG_Height) / 2
        Ports_mid = max_slabH


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
    
    
        # # Device specifications
        max_slabH = Slab_Height
        MonitorHeight = Substrate_Height + (max_slabH + WG_Height) / 2
        Ports_mid = (max_slabH + WG_Width) / 2
        Ports_PWB_mid = (max_slabH + TaperHightB) / 2
        X_min = -TaperLength_PWB/2 - 0.5e-6
    
    
    
        # Adds a Eigenmode Expansion (EME) solver region to the MODE simulation environment.
        self.lum.addeme()
        self.lum.set("x min", X_min)
        self.lum.set("y", 0)
        self.lum.set("y span", TaperWidthB + TaperWidthB/2)
        self.lum.set('simulation temperature', 273.15 + 20)
        self.lum.set("z", TaperHightB/2) #Substrate_Height
        self.lum.set("z span", TaperHightB*2) 
        self.lum.set("wavelength", WaveLength)
        self.lum.set("z min bc", "PML")
        self.lum.set("z max bc", "PML")
        self.lum.set("y min bc", "PML")
        self.lum.set("y max bc", "PML")
        # set cell properties
        self.lum.set("number of cell groups", 3)
        self.lum.set("group spans", np.array([[0.1e-6], [TaperLength_PWB], [9.9e-6]]))
        self.lum.set("cells", np.array([[3], [80], [3]]))
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
        z_span_PWB = TaperHightB + (z_Port_Span - WG_Height)/2
    
        self.lum.select("EME::Ports::port_" + str(1))
        self.lum.set("port location", portLoc[0])
        self.lum.set("use full simulation span", 1)
        # self.lum.set("use full simulation span", 0)
        # self.lum.set("y", 0)
        # self.lum.set("y span", TaperWidthB + y_Port_Span/2)
        # self.lum.set("z", Ports_PWB_mid)
        # self.lum.set("z span", z_span_PWB)
        self.lum.set("mode selection", Mode)
    
        self.lum.select("EME::Ports::port_" + str(2))
        self.lum.set("port location", portLoc[1])
        self.lum.set("use full simulation span", 1)
        # self.lum.set("use full simulation span", 0)
        # self.lum.set("y", 0)
        # self.lum.set("y span", y_Port_Span)
        # self.lum.set("z", Ports_mid)
        # self.lum.set("z span", z_Port_Span)
        self.lum.set("mode selection", Mode)
    
    
    
        # Add monitor Horizondal
        self.lum.addemeprofile()
        # self.lum.set("x", 0)
        # self.lum.set("x span", TaperLength_PWB+11e-6)
        self.lum.set("x min", X_min - 5e-6)
        self.lum.set("x max", TaperLength_PWB/2 + 15e-6)
        self.lum.set("y", 0)
        self.lum.set("y span", TaperWidthB*2)
        self.lum.set("z", MonitorHeight)
        
        
        # Add monitor Vertical   
        self.lum.addemeprofile()
        self.lum.set('monitor type', "2D Y-normal")
        # self.lum.set("x", 0)
        # self.lum.set("x span", TaperLength_PWB+11e-6)
        self.lum.set("x min", X_min - 5e-6)
        self.lum.set("x max", TaperLength_PWB/2 + 15e-6)
        self.lum.set("y", 0)
        # self.lum.set("y span", TaperWidthB*2)
        self.lum.set("z", MonitorHeight)
        self.lum.set("z span", TaperHightB*2)







    def setMMI1x2VarFDTDSolver(self, Parameters):

        '''


        Slab_Height : int/float
            Height of the slab.
        MMI_Width : int/float
            Width of MMI
        MMI_Length : int/float
            Length of MMI
        angle : int/float
            Angle of the Waveguide Walls. it is calculated WG_angle = 90 - angle.
            For anfle = 90 we get a perfect rect!
        WG_Height : int/float
            Heigth of waveguide
        WG_Width : int/float
            Width og waveguide
        WG_Length : int/float
            Length of waveguide
        posOffset : int/float
            Offset between the waveguides. If Taper == True then this become the offset
            betweent he tapers wider sides. Waveguide and Tapers cannot be placed ourside
            the MMI structure. The minimum distance between Taper and Waveguide is 1 um
            becouse of manufactering restrictions in the University.
        OffsetInpit : int/float
            Input waveguide/taper offset.
        y_res : int/float
            Mesh cell sizes.
        z_res : int/float
            Mesh cell size.
        Slab_Height : int/float
            Height of the slab.
        WaveLength : int/float
            Wavelength.
        Taper : boolen, optional
            Defoult Taper == False
                if Taper == True Simulation region will ajust for simulation with tapers
        TaperLength : int/float
            Length of taper.
        Parameters["Mode"] : str
                    Mode to choose from ("fundamental TE mode", "fundamental TM mode", "fundamental mode")

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
        angle = Parameters['angle']
        WG_Height = Parameters['WG Height']
        WG_Width = Parameters['WG Width']
        WG_Length = Parameters['WG Length']
        OffsetInput = Parameters['Offset Input']
        posOffset = Parameters['Position Offset']
        x_res = Parameters['x res']
        y_res = Parameters['y res']
        z_res = Parameters['z res']
        Slab_Height = Parameters['Slab Height']
        WaveLength = Parameters['Wavelength']
        Mode = Parameters["Mode"]

        # Device specifications
        Device_Length = MMI_Length + 2 * WG_Length
        Device_Width = MMI_Width + WaveLength * 2  # MMI_Width
        max_subH = Substrate_Height / 2
        max_slabH = Slab_Height
        MonitorHeight = Substrate_Height + (max_slabH + WG_Height) / 2
        Lnoi_H = Slab_Height + WG_Height
        Ports_mid = Substrate_Height + (max_slabH + WG_Height) / 2
        WG_H = WG_Height

        # Adds a Eigenmode Expansion (EME) solver region to the MODE simulation environment.
        self.lum.addvarfdtd()
        self.lum.set("x min", -(Device_Length / 2))
        self.lum.set("x max", (Device_Length / 2))
        self.lum.set("y", 0)
        self.lum.set('simulation temperature', 273.15 + 20)
        self.lum.set("y span", Device_Width)
        self.lum.set("z", Substrate_Height)
        self.lum.set("z span", 4e-6)
        self.lum.set('number of test points', 1)
        self.lum.set('z min bc', 'PML')
        self.lum.set('z max bc', 'PML')
        self.lum.set('mesh type', 'auto non-uniform')
        self.lum.set('dx', x_res)
        self.lum.set('dy', y_res)
        self.lum.set('dz', z_res)
        self.lum.set('set simulation bandwidth', 0)
        self.lum.set('global source center wavelength', WaveLength)
        self.lum.set('global source wavelength span', 0)

        # define two values for upper and lower limit of the WG offsets
        OffMin = -MMI_Width / 2
        OffMax = MMI_Width / 2

        x = abs(WG_H / (np.cos((angle) * np.pi / 180)))  # in Radians
        extention = np.sqrt(x ** 2 - WG_H ** 2)
        WG_W = WG_Width + 2 * extention

        offset_Set_R = posOffset / 2 + WG_W / 2 + WG_Width / 2
        xPos = [(Device_Length / 2) - 2e-6, -(Device_Length / 2) + 2e-6, -(Device_Length / 2) + 2e-6]
        yPos_min = [0, -(WG_Width / 2 + posOffset / 2), WG_Width / 2 + posOffset / 2]
        name = ['Input', 'Output_L', 'Output_R']

        self.lum.addmodesource()
        self.lum.set('direction', 'Backward')
        self.lum.set("x", (Device_Length / 2) - 1e-7)
        self.lum.set("y min", -(WG_W / 2 + OffsetInput))
        self.lum.set("y max", (WG_W / 2 + OffsetInput))
        self.lum.set('override global source settings', 0)

        self.lum.addmovie()
        self.lum.set("x min", -(0.5e-6 + Device_Length / 2))
        self.lum.set("x max", (0.5e-6 + Device_Length / 2))
        self.lum.set("y min", -(Device_Width + 1e-6))
        self.lum.set("y max", (Device_Width + 1e-6))

        for i in range(3):
            self.lum.addtime()
            self.lum.set('name', name[i])
            self.lum.set("x", xPos[i])
            self.lum.set("y", yPos_min[i])
            self.lum.set('output Px', 1)
            self.lum.set('output Py', 1)
            self.lum.set('output Pz', 1)

        # # Add monitor
        # self.lum.addpower()
        # self.lum.set('monitor type', '2D X-normal')
        # self.lum.set("x min", -(Device_Length / 2))
        # self.lum.set("x max", (Device_Length / 2))
        # self.lum.set("y min", -Device_Width / 2)
        # self.lum.set("y max", Device_Width / 2)
        # self.lum.set("z", MonitorHeight)



    def setMMI2x2VarFDTDSolver(self, Parameters):

        '''


        Slab_Height : int/float
            Height of the slab.
        MMI_Width : int/float
            Width of MMI
        MMI_Length : int/float
            Length of MMI
        angle : int/float
            Angle of the Waveguide Walls. it is calculated WG_angle = 90 - angle.
            For anfle = 90 we get a perfect rect!
        WG_Height : int/float
            Heigth of waveguide
        WG_Width : int/float
            Width og waveguide
        WG_Length : int/float
            Length of waveguide
        posOffset : int/float
            Offset between the waveguides. If Taper == True then this become the offset
            betweent he tapers wider sides. Waveguide and Tapers cannot be placed ourside
            the MMI structure. The minimum distance between Taper and Waveguide is 1 um
            becouse of manufactering restrictions in the University.
        OffsetInpit : int/float
            Input waveguide/taper offset.
        y_res : int/float
            Mesh cell sizes.
        z_res : int/float
            Mesh cell size.
        Slab_Height : int/float
            Height of the slab.
        WaveLength : int/float
            Wavelength.
        Taper : boolen, optional
            Defoult Taper == False
                if Taper == True Simulation region will ajust for simulation with tapers
        TaperLength : int/float
            Length of taper.
        Parameters["Mode"] : str
                    Mode to choose from ("fundamental TE mode", "fundamental TM mode", "fundamental mode")

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
        angle = Parameters['angle']
        WG_Height = Parameters['WG Height']
        WG_Width = Parameters['WG Width']
        WG_Length = Parameters['WG Length']
        posOffset = Parameters['Position Offset']
        x_res = Parameters['x res']
        y_res = Parameters['y res']
        z_res = Parameters['z res']
        Slab_Height = Parameters['Slab Height']
        WaveLength = Parameters['Wavelength']
        Mode = Parameters["Mode"]


        # Device specifications
        Device_Length = MMI_Length + 2 * WG_Length
        Device_Width = MMI_Width + WaveLength * 2  # MMI_Width
        max_subH = Substrate_Height / 2
        max_slabH = Slab_Height
        MonitorHeight = Substrate_Height + (max_slabH + WG_Height) / 2
        Lnoi_H = Slab_Height + WG_Height
        Ports_mid = Substrate_Height + (max_slabH + WG_Height) / 2
        WG_H = WG_Height

        # Adds a Eigenmode Expansion (EME) solver region to the MODE simulation environment.
        self.lum.addvarfdtd()
        self.lum.set("x min", -(Device_Length / 2))
        self.lum.set("x max", (Device_Length / 2))
        self.lum.set("y", 0)
        self.lum.set('simulation temperature', 273.15 + 20)
        self.lum.set("y span", Device_Width)
        self.lum.set("z", Substrate_Height)
        self.lum.set("z span", 4e-6)
        self.lum.set('number of test points', 1)
        self.lum.set('z min bc', 'PML')
        self.lum.set('z max bc', 'PML')
        self.lum.set('mesh type', 'auto non-uniform')
        self.lum.set('dx', x_res)
        self.lum.set('dy', y_res)
        self.lum.set('dz', z_res)
        self.lum.set('set simulation bandwidth', 0)
        self.lum.set('global source center wavelength', WaveLength)
        self.lum.set('global source wavelength span', 0)

        # define two values for upper and lower limit of the WG offsets
        OffMin = -MMI_Width / 2
        OffMax = MMI_Width / 2

        x = abs(WG_H / (np.cos((angle) * np.pi / 180)))  # in Radians
        extention = np.sqrt(x ** 2 - WG_H ** 2)
        WG_W = WG_Width + 2 * extention

        offset_Set_R = posOffset / 2 + WG_W / 2 + WG_Width / 2
        xPos = [(Device_Length / 2) - 2e-6, (Device_Length / 2) - 2e-6, -(Device_Length / 2) + 2e-6,
                -(Device_Length / 2) + 2e-6]
        yPos_min = [-(WG_Width / 2 + posOffset / 2), WG_Width / 2 + posOffset / 2, -(WG_Width / 2 + posOffset / 2),
                    WG_Width / 2 + posOffset / 2]
        name = ['Input_L', 'Input_R', 'Output_L', 'Output_R']

        self.lum.addmodesource()
        self.lum.set('direction', 'Backward')
        self.lum.set("x", (Device_Length / 2) - 1e-7)
        self.lum.set("y", WG_Width / 2 + posOffset / 2)
        self.lum.set("y span", WG_Width * 2)
        self.lum.set('override global source settings', 0)

        self.lum.addmovie()
        self.lum.set("x min", -(0.5e-6 + Device_Length / 2))
        self.lum.set("x max", (0.5e-6 + Device_Length / 2))
        self.lum.set("y min", -(Device_Width + 1e-6))
        self.lum.set("y max", (Device_Width + 1e-6))

        for i in range(4):
            self.lum.addtime()
            self.lum.set('name', name[i])
            self.lum.set("x", xPos[i])
            self.lum.set("y", yPos_min[i])
            self.lum.set('output Px', 1)
            self.lum.set('output Py', 1)
            self.lum.set('output Pz', 1)

        # # Add monitor
        # self.lum.addpower()
        # self.lum.set('monitor type', '2D X-normal')
        # self.lum.set("x min", -(Device_Length / 2))
        # self.lum.set("x max", (Device_Length / 2))
        # self.lum.set("y min", -Device_Width / 2)
        # self.lum.set("y max", Device_Width / 2)
        # self.lum.set("z", MonitorHeight)




    def setWaveguideFDESolver(self, Parameters):
        '''


        Parameters
        ----------
        Parameters : Dictionary
            Dictionary with all the data needet for the Bend Wavaguide. Data needet:
                Parameters['angle'] : int/float
                    Angle of the Waveguide Walls. it is calculated WG_angle = 90 - angle.
                    For anfle = 90 we get a perfect rect!
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
        self.lum.set('z span', 2e-6)
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







class HelpSubject:
    
    @classmethod
    def Help_Objects(self):
            
            print("""
                  
                  Parameters for each object are given as dictionary!
                  
                  1) Waveguide(Parameters)
            
                    Parameters['WG Length']: int/float
                        Length of the Waveguide
                    Parameters['WG Height'] : int/float
                        Height of the Waveguide (Etching depth)
                    Parameters['WG Width'] : int/float
                        Top width of the Waveguide.
                    Parameters['angle'] : int
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
                        
                        
                        
                        
                        
                        
                  2) StraightWaveguide(Parameters)
                                                     
                    Parameters['Substrate Height'] : int/float
                        Substrate Height
                    Parameters['WG Height'] : int/float
                        Height of the Waveguide (Etching depth)
                    Parameters['WG Width'] : int/float
                        Top width of the Waveguide.
                    Parameters['WG Length'] : int/float
                        Length of the Waveguide
                    Parameters['Taper'] : boolen
                        Taper can be set to be True ot False.
                        if Taper == False - No Taper used, only waveguide will be used
                        if Taper == True - Taper placed, no waveguides only taper will be used
                    Parameters['Taper Length'] : int/float
                        Taper Length
                    Parameters['Taper Width'] : int/float
                        Taper backside width, the frontside width is the waveguide width
                    Parameters['angle'] : int
                        Angle of the Waveguide Walls. It is calculated WG_angle = 90 - angle.
                        For angle = 90 , a perfect rect is created!
                    Parameters['Slab Height'] : int/float
                        Slab height.
                    Parameters['Material'] : list of str
                        List of Materials. The list should be with names (str) of a valid Lumerical materials.
                        Check the names in Lumerical Materials viewer.
                        The List of materials must contain at least 2 materials! 
                        Parameters['Material'] = ['Cladding/Substrat', 'Object Material'].
                        For Example: Parameters['Material'] = ["SiO2 (Glass) - Palik", 'LiNbO3_20deg_X cut'].
                    Parameters["Wavelength"] : int/float
                        Wavelength
                    Parameters["Waveguide Angle"] : int/float
                        Tilting angle of the Waveguide. Set it to 0 to simulate the straight waveguide. 
                        If Waveguide Angle is different then 0, then the straight waveguide will be tilted 
                        at the choosen angle. 
                        
                        
                        
                  3) BendWaveguide(Parameters)  
        
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
                        
                        
                        
                  4) ArcWaveguide(Parameters) 
        
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
                    Parameters["S_Band Radius"] : int/float
                        Radius of the Circle. The Radius is given in um. 
                    Parameters['Arc deg'] : int
                        Can be 90 or 180 for 1/4 of a cirle or 1/2 of a circle.
                        
                        
                        
                  5) MMI2x1(Parameters)    
                    
                    Parameters['Material'] : list of str
                        List of Materials. The list should be with names (str) of a valid Lumerical materials.
                        Check the names in Lumerical Materials viewer.
                        The List of materials must contain at least 2 materials! 
                        Parameters['Material'] = ['Cladding/Substrat', 'Object Material'].
                        For Example: Parameters['Material'] = ["SiO2 (Glass) - Palik", 'LiNbO3_20deg_X cut'].
                    Parameters['Substrate Height'] : int/float
                        Substrate Height
                    Parameters['MMI Width'] : int/float
                        Top width of the MMI.
                    Parameters['MMI Length'] : int/float
                        Length of the MMI.
                    Parameters['Slab Height'] : int/float
                        Height of the LiNbO3 Slab
                    Parameters['angle'] : int/float
                        Angle of the Waveguide Walls. It is calculated WG_angle = 90 - angle.
                        For angle = 90 , a perfect rect is created!
                    Parameters['WG Height'] : int/float
                        Height of the Waveguide (Etching depth).
                    Parameters['WG Width'] : int/float
                        Top width of the Waveguide.
                    Parameters['WG Length'] : int/float
                        Waveguide length
                    Parameters['Taper'] : boolen
                        Taper can be set to be True ot False.
                        if Taper == False - No Taper used
                        if Taper == True - Taper placed
                    Parameters['Position Offset'] : int/float
                        Offset between the waveguides. If Taper == True then this become the offset
                        betweent he tapers wider sides. Waveguide and Tapers cannot be placed ourside
                        the MMI structure. The minimum distance between Taper and Waveguide is 1 um
                        becouse of manufactering restrictions in the University.
                    Parameters['Offset Input'] : int/float
                        Offset of the input waveguide.
                    Parameters['Wavelength'] : int/float
                        Wavelength
                    Parameters['Taper Length'] : int/float
                        If Taper == True, then this will set the Tapers length. If Taper == False
                        this will be ignored and some random value can be given.
                    Parameters['Taper Width'] : int/float
                        If Taper == True, then this will set the Tapers width. If Taper == False
                        this will be ignored and some random value can be given.
                        
                        
                        
                  6) MMI2x2(Parameters)   
    
                    Parameters['Material'] : list of str
                        List of Materials. The list should be with names (str) of a valid Lumerical materials.
                        Check the names in Lumerical Materials viewer.
                        The List of materials must contain at least 2 materials! 
                        Parameters['Material'] = ['Cladding/Substrat', 'Object Material'].
                        For Example: Parameters['Material'] = ["SiO2 (Glass) - Palik", 'LiNbO3_20deg_X cut'].
                    Parameters['Substrate Height'] : int/float
                        Substrate Height
                    Parameters['MMI Width'] : int/float
                        Top width of the MMI.
                    Parameters['MMI Length'] : int/float
                        Length of the MMI.
                    Parameters['Slab Height'] : int/float
                        Height of the LiNbO3 Slab
                    Parameters['angle'] : int/float
                        Angle of the Waveguide Walls. It is calculated WG_angle = 90 - angle.
                        For angle = 90 , a perfect rect is created!
                    Parameters['WG Height'] : int/float
                        Height of the Waveguide (Etching depth). This is also the
                    Parameters['WG Width'] : int/float
                        Top width of the Waveguide.
                    Parameters['WG Length'] : int/float
                        Waveguide length
                    Parameters['Position Offset'] : int/float
                        Offset between the waveguides. If Taper == True then this become the offset
                        betweent he tapers wider sides. Waveguide and Tapers cannot be placed ourside
                        the MMI structure. The minimum distance between Taper and Waveguide is 1 um
                        becouse of manufactering restrictions in the University.
                    Parameters['Wavelength'] : int/float
                        Wavelength
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
                        
                        
                        
                        
                  7) WDM(Parameters)
                                              
                    Parameters['Material'] : list of str
                        List of Materials. The list should be with names (str) of a valid Lumerical materials.
                        Check the names in Lumerical Materials viewer.
                        The List of materials must contain at least 2 materials! 
                        Parameters['Material'] = ['Cladding/Substrat', 'Object Material'].
                        For Example: Parameters['Material'] = ["SiO2 (Glass) - Palik", 'LiNbO3_20deg_X cut'].
                        
                    Parameters['Substrate Height'] : int/float
                       Substrate Height
                    Parameters['MMI Width'] : int/float
                        Top width of the MMI.
                    Parameters['MMI Length'] : int/float
                        Length of the MMI.
                    Parameters['angle'] : int/float
                        Angle of the Waveguide Walls. It is calculated WG_angle = 90 - angle.
                        For angle = 90 , a perfect rect is created!
                    Parameters['WG Height'] : int/float
                         Height of the Waveguide (Etching depth). This is also the
                    Parameters['WG Width'] : int/float
                         Top width of the Waveguide
                    Parameters['WG Length'] : int/float
                         Waveguide length. For the moment please choose the Parameters['WG Length'] to be 1e-6 bigger then Parameters['Taper Length']!!
                    Parameters['Slab Height'] : int/float
                        Height of the LiNbO3 Slab
                    Parameters['Wavelength'] : int/float
                        Wavelength
                    Parameters['Angle Thetha'] : 
                        Input and output angle of the waveguide. This is the rotation angle of the waveguide and taper. 
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
                        
                        
                        
                  8) InverseTaper(Paraemters)                                                                      

                    Parameters['Material'] : list of str
                        List of Materials. The list should be with names (str) of a valid Lumerical materials.
                        Check the names in Lumerical Materials viewer.
                        The List of materials must contain at least 3 materials! 
                        Parameters['Material'] = ['Cladding/Substrat', 'Object Material', 'Photonic Wire Bonding'].
                        For Example: Parameters['Material'] = ["SiO2 (Glass) - Palik", 'LiNbO3_20deg_X cut', 'Si (Silicon) - Palik'].
                    Parameters['Substrate Height'] : int/float
                       Substrate Height
                    Parameters['angle'] : int/float
                        Angle of the Waveguide Walls. It is calculated WG_angle = 90 - angle.
                        For angle = 90 , a perfect rect is created!
                        For anfle = 90 we get a perfect rect!
                    Parameters['WG Height'] : int/float
                        Height of the Waveguide (Etching depth). This is also the
                    Parameters['WG Width'] : int/float
                        Top width of the waveguide. Also in this function and ONLY in this function this will be the 
                        ibverse Taper backside width!!!
                    Parameters['Slab Height'] : int/float
                        Height of the LiNbO3 Slab
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
                        Photonic Wirebonding (PWB) Width front side (to the waveguide) 
                    Parameters['PWB Taper Hight Front'] : int/float
                        Photonic Wire Bonding Height front side (to the waveguide) 
                    Parameters['PWB Taper Length'] : int/float
                        Length of the Photonic Wire Bonding Taper       
                        
                        
                        
                9) MMI2x1_Trapez(Parameters)    
                  
                  Parameters['Material'] : list of str
                      List of Materials. The list should be with names (str) of a valid Lumerical materials.
                      Check the names in Lumerical Materials viewer.
                      The List of materials must contain at least 2 materials! 
                      Parameters['Material'] = ['Cladding/Substrat', 'Object Material'].
                      For Example: Parameters['Material'] = ["SiO2 (Glass) - Palik", 'LiNbO3_20deg_X cut'].
                  Parameters['Substrate Height'] : int/float
                      Substrate Height
                  Parameters['MMI Width'] : int/float
                      Top width of the MMI.
                  Parameters['MMI Length'] : int/float
                      Length of the MMI.
                  Parameters['Middle MMI Length']: int/float
                      Length of the MMI in the Middle Part due to inperfection by manufacturing
                  Parameters['Slab Height'] : int/float
                      Height of the LiNbO3 Slab
                  Parameters['angle'] : int/float
                      Angle of the Waveguide Walls. It is calculated WG_angle = 90 - angle.
                      For angle = 90 , a perfect rect is created!
                  Parameters['WG Height'] : int/float
                      Height of the Waveguide (Etching depth).
                  Parameters['WG Width'] : int/float
                      Top width of the Waveguide.
                  Parameters['WG Length'] : int/float
                      Waveguide length
                  Parameters['Taper'] : boolen
                      Taper can be set to be True ot False.
                      if Taper == False - No Taper used
                      if Taper == True - Taper placed
                  Parameters['Position Offset'] : int/float
                      Offset between the waveguides. If Taper == True then this become the offset
                      betweent he tapers wider sides. Waveguide and Tapers cannot be placed ourside
                      the MMI structure. The minimum distance between Taper and Waveguide is 1 um
                      becouse of manufactering restrictions in the University.
                  Parameters['Offset Input'] : int/float
                      Offset of the input waveguide.
                  Parameters['Wavelength'] : int/float
                      Wavelength
                  Parameters['Taper Length'] : int/float
                      If Taper == True, then this will set the Tapers length. If Taper == False
                      this will be ignored and some random value can be given.
                  Parameters['Taper Width'] : int/float
                      If Taper == True, then this will set the Tapers width. If Taper == False
                      this will be ignored and some random value can be given.
                      
                      
                      
                10) MMI2x2_Trapez(Parameters)   
  
                  Parameters['Material'] : list of str
                      List of Materials. The list should be with names (str) of a valid Lumerical materials.
                      Check the names in Lumerical Materials viewer.
                      The List of materials must contain at least 2 materials! 
                      Parameters['Material'] = ['Cladding/Substrat', 'Object Material'].
                      For Example: Parameters['Material'] = ["SiO2 (Glass) - Palik", 'LiNbO3_20deg_X cut'].
                  Parameters['Substrate Height'] : int/float
                      Substrate Height
                  Parameters['MMI Width'] : int/float
                      Top width of the MMI.
                  Parameters['MMI Length'] : int/float
                      Length of the MMI.
                  Parameters['Middle MMI Length']: int/float
                      Length of the MMI in the Middle Part due to inperfection by manufacturing
                  Parameters['Slab Height'] : int/float
                      Height of the LiNbO3 Slab
                  Parameters['angle'] : int/float
                      Angle of the Waveguide Walls. It is calculated WG_angle = 90 - angle.
                      For angle = 90 , a perfect rect is created!
                  Parameters['WG Height'] : int/float
                      Height of the Waveguide (Etching depth). This is also the
                  Parameters['WG Width'] : int/float
                      Top width of the Waveguide.
                  Parameters['WG Length'] : int/float
                      Waveguide length
                  Parameters['Position Offset'] : int/float
                      Offset between the waveguides. If Taper == True then this become the offset
                      betweent he tapers wider sides. Waveguide and Tapers cannot be placed ourside
                      the MMI structure. The minimum distance between Taper and Waveguide is 1 um
                      becouse of manufactering restrictions in the University.
                  Parameters['Wavelength'] : int/float
                      Wavelength
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
                        
                  
                  """)
                  
                  
    def Help_Solvers(self):
        print("""
         
                  1) Solver("Waveguide", "FDE",  Parameters)                                                        
                                                                                                                    
                    Parameters['angle'] : int/float                                                                 
                        Angle of the Waveguide Walls. It is calculated WG_angle = 90 - angle.                       
                        For angle = 90 , a perfect rect is created!                                                 
                    Parameters['WG Height'] : int/float                                                             
                        Height of the Waveguide (Etching depth).                                                    
                    Parameters['y res'] : int/float                                                                 
                        Maximum mesh step setting in y direction. (Size of the Mesh Box in y direction in um)       
                    Parameters['z res'] : int/float                                                                 
                        Maximum mesh step setting in z direction. (Size of the Mesh Box in y direction in um)       
                    Parameters['Slab Height'] : int/float                                                           
                        Height of the LiNbO3 Slab                                                                    
                    Parameters['Wavelength'] : int/float                                                            
                        Wavelength.                                                                                  
        
                    
        
                  2) Solver("MMI2x1", "EME", Parameters) and Solver("MMI2x2", "EME", Parameters)
                    
                    Parameters['Substrate Height'] : int/float
                        Substrate Height
                    Parameters['MMI Width'] : int/float
                        Top width of the MMI
                    Parameters['MMI Length'] : int/float
                        Length of MMI
                    Parameters['angle'] : int/float
                        Angle of the Waveguide Walls. It is calculated WG_angle = 90 - angle.
                        For angle = 90 , a perfect rect is created!
                    Parameters['WG Height'] : int/float
                        Height of the Waveguide (Etching depth)
                    Parameters['WG Width'] : int/float
                        Top width of the waveguide
                    Parameters['WG Length'] : int/float
                        Length of the Waveguide
                    Parameters['Position Offset'] : int/float
                        Offset between the waveguides. If Taper == True then this become the offset
                        betweent he tapers wider sides. Waveguide and Tapers cannot be placed ourside
                        the MMI structure. The minimum distance between Taper and Waveguide is 1 um
                        becouse of manufactering restrictions in the University.
                        
                    Only for 2x1MMI !!!!!
                        Parameters['Offset Input'] : int/float
                            Input waveguide/taper offset.
                    
                    Parameters['y res'] : int/float
                        Maximum mesh step setting in y direction. (Size of the Mesh Box in y direction in um)
                    Parameters['z res'] : int/float
                        Maximum mesh step setting in z direction. (Size of the Mesh Box in z direction in um)
                    Parameters['Slab Height'] : int/float
                        Height of the LiNbO3 Slab
                    Parameters['Wavelength'] : int/float
                        Wavelength.
                    Parameters["Mode"] : str
                        Mode to analyse. ("fundamental TE mode", "fundamental TM mode", "fundamental mode")
                    
                    
        
                  3) Solver('DirectionalCoupler', 'EME', Parameters)
        
                    Parameters['Substrate Height'] : float/int
                        Height of the Substrate
                    Parameters['Substrate Width'] : float/int
                        Width of the Substrate
                    Parameters['DC Length'] : float/int
                        Length of the Directional Coupler
                    Parameters['WG Height'] : float/int
                        Height of the Waveguide (Etching depth)
                    Parameters['WG Width'] : float/int
                        Top width of the Waveguide 
                    Parameters['Position Offset'] : float/int
                        Positional offser of the waveguides. If posOffset the two Waveguides
                        will be offset of the middle position (y = 0) by the half of there
                        Width. In this case they will not overlap if the Offset is 0.
                    Parameters['y res'] : float/int
                        Maximum mesh step setting in y direction. (Size of the Mesh Box in y direction in um)
                    Parameters['z res'] : float/int
                        Maximum mesh step setting in z direction. (Size of the Mesh Box in y direction in um)
                    Parameters['Slab Height'] : float/int
                        Height of the LiNbO3 Slab
                    Parameters['Wavelength'] : float/int
                        Wavelength
                    Parameters["Mode"] : str
                        Mode to analyse. ("fundamental TE mode", "fundamental TM mode", "fundamental mode")
                    
                    
        
                  4) Solver("StraightWaveguide", "EME",  Parameters)  
                     
                                                       
                    Parameters['Substrate Height'] : int/float
                        Substrate height.
                    Parameters['WG Height'] : int/float
                        Height of the Waveguide (Etching depth)
                    Parameters['WG Width'] : int/float
                        Top width of the Waveguide
                    Parameters['WG Length'] : int/float
                        Length of the Waveguide
                    Parameters['Taper'] : boolen
                        Taper can be set to be True ot False.
                        if Taper == False - No Taper used, only waveguide will be used
                        if Taper == True - Taper placed, no waveguides only taper will be used
                    Parameters['Taper Length'] : int/float
                        Taper Length
                    Parameters['Taper Width'] : int/float
                        Taper backside top waveguide width, the frontside top width of the waveguide is the waveguide top width
                    Parameters['y res']: int/float
                          Maximum mesh step setting in y direction. (Size of the Mesh Box in y direction in um)
                    Parameters['z res']: int/float
                          Maximum mesh step setting in z direction. (Size of the Mesh Box in z direction in um)
                    Parameters['Slab Height'] : int/float
                        Height of the LiNbO3 Slab
                    Parameters['Wavelength'] : int/float
                        Wavelength
                    Parameters["Waveguide Angle"] : int/float
                        This Parameter will set the theta ratation angle of the port. Please choose the same angle as you gave for the Object StraightWaveguide!
                    Parameters["Port Span"] : list of floats/ints
                        List of x,y and z span of the Ports. For this simulation only y and z parametes will be taken, since we use x cut material.
                    Parameters["Mode"] : str
                         Mode to analyse. ("fundamental TE mode", "fundamental TM mode", "fundamental mode")
                         
                    
                    
                  5) Solver("WDM", "EME",  Parameters)
                                                
                    Parameters['Substrate Height'] : int/float
                        Substrate height.
                    Parameters['MMI Width'] : int/float
                        Top width of the MMI.
                    Parameters['MMI Length'] : int/float
                        Length of the MMI.
                    Parameters['WG Height'] : int/float
                        Height of the Waveguide (Etching depth)
                    Parameters['WG Length'] : int/float
                        Length of the Waveguide
                    Parameters['WG Width'] : int/float
                        Top width of the Waveguide
                    Parameters['Taper Width'] : int/float
                        Taper backside top width of the waveguide, the frontside width is the top width of the waveguide
                    Parameters['Taper Length'] : int/float
                        Taper Length
                    Parameters['y res']: int/float
                          Maximum mesh step setting in y direction. (Size of the Mesh Box in y direction in um)
                    Parameters['z res']: int/float
                          Maximum mesh step setting in z direction. (Size of the Mesh Box in z direction in um)
                    Parameters['Slab Height'] : int/float
                        Slab Height.
                    Parameters['Wavelength'] : int/float
                        Wavelength
                    Parameters['Angle Thetha'] : boolen
                        Input and output angle of the waveguide. This is the rotation angle of the waveguide and taper. 
                    Parameters["Port Span"] : list of floats/ints
                        List of x,y and z span of the Ports. For this simulation only y and z parametes will be taken, since we use x cut material.
                    Parameters["Mode"] : str
                         Mode to analyse. ("fundamental TE mode", "fundamental TM mode", "fundamental mode") 
    
    
    
    
                  6) Solver("InverseTaper", "EME", Parameters)
                                                  
                    Parameters['Substrate Height'] : int/float
                        Substrate height.
                    Parameters['WG Height'] : int/float
                        Height of the Waveguide (Etching depth)
                    Parameters['WG Width'] : int/float
                        Top width of the Waveguide.
                    Parameters['Slab Height'] : int/float
                        Slab Height.
                    Parameters['PWB Taper Width Back'] : int/float
                        Photonic Wirebonding (PWB) Width back side (to the Photonic Wire Bonding) 
                    Parameters['PWB Taper Hight Back'] : int/float
                        Photonic Wire Bonding Height back side (to the Photonic Wire Bonding) 
                    Parameters['PWB Taper Length'] : int/float
                        Length of the Photonic Wire Bonding Taper      
                    Parameters['y res']: int/float
                          Maximum mesh step setting in y direction. (Size of the Mesh Box in y direction in um)
                    Parameters['z res']: int/float
                          Maximum mesh step setting in z direction. (Size of the Mesh Box in z direction in um)
                    Parameters['Wavelength'] : int/float
                        Wavelength
                    Parameters["Port Span"] : list of floats/ints
                        List of x,y and z span of the Ports. For this simulation only y and z parametes will be taken, since we use x cut material.
                    Parameters["Mode"] : str
                         Mode to analyse. ("fundamental TE mode", "fundamental TM mode", "fundamental mode") 
                         
      
                    
                    
        
                  7) Solver("MMI2x1", "FDTD", Parameters) and Solver("MMI2x2", "FDTD", Parameters) 
        
                    Parameters['Substrate Height'] : int/float
                        Substrate height.
                    Parameters['MMI Width'] : int/float
                        Top width of the MMI.
                    Parameters['MMI Length'] : int/float
                        Length of the MMI.
                    Parameters['WG Height' : int/float
                        Waveguide hight
                    Parameters['WG Width'] : int/float
                        Top width of the Waveguide.
                    Parameters['WG Length'] : int/float
                        Waveguide length.
                    Parameters['y res'] : int/float
                        Maximum mesh step setting in y direction. (Size of the Mesh Box in y direction in um)
                    Parameters['z res'] : int/float
                        Maximum mesh step setting in z direction. (Size of the Mesh Box in z direction in um)
                    Parameters['Slab Height'] : int/float
                        Height of the LiNbO3 Slab
                    Parameters['Wavelength'] : int/float
                        Wavelength
                    Parameters['Angle Thetha'] : int/float
                        Angle for the input and output waveguides 
                    Parameters["Mode"] : str
                        Mode to analyse. ("fundamental TE mode", "fundamental TM mode", "fundamental mode")
                   
                    
        
                  8) Solver("MMI2x1", "FDTD", Parameters) and Solver("MMI2x2", "FDTD", Parameters)
                    
                    Parameters['Substrate Height'] : int/float
                        Substrate height.
                    Parameters['MMI Width'] : int/float
                        Top width of the MMI
                    Parameters['MMI Length'] : int/float
                        Length of MMI
                    Parameters['WG Height'] : int/float
                        Heigth of the waveguide
                    Parameters['WG Width'] : int/float
                        Top width of the waveguide
                    Parameters['WG Length'] : int/float
                        Length of the waveguide
                    Parameters['Position Offset'] : int/float
                        Offset between the waveguides. If Taper == True then this become the offset
                        betweent he tapers wider sides. Waveguide and Tapers cannot be placed ourside
                        the MMI structure. The minimum distance between Taper and Waveguide is 1 um
                        becouse of manufactering restrictions in the University.
                    
                    Only for 2x1MMI !!!!!
                        Parameters['Offset Input'] : int/float
                            Input waveguide/taper offset.
         
                    Parameters['x res'] : int/float
                        Minimum mesh step setting in x direction. (Size of the Mesh Box in x direction in um)
                    Parameters['Slab Height'] : int/float
                        Height of the LiNbO3 Slab
                    Parameters['Wavelength'] : int/float
                        Wavelength.
                    Parameters["Mode"] : str
                        Mode to analyse. ("fundamental TE mode", "fundamental TM mode", "fundamental mode")
                    
                    
        
                  9) Solver("StraightWaveguide", "FDTD",  Parameters)    
                             
                    Parameters['Substrate Height'] : int/float
                        Substrate height.
                    Parameters['WG Length'] : int/float
                        Waveguide Length
                    Parameters['WG Height'] : int/float
                        Waveguide hight
                    Parameters['WG Width'] : int/float
                        Top width of the Waveguide
                    Parameters['WG Length'] : int/float
                        Length of the Waveguide
                    Parameters['Taper'] : boolen
                        Taper can be set to be True ot False.
                        if Taper == False - No Taper used, only waveguide will be used
                        if Taper == True - Taper placed, no waveguides only taper will be used
                    Parameters['Taper Length'] : int/float
                        Taper Length
                    Parameters['Taper Width'] : int/float
                        Taper backside top width of the waveguide, the frontside width is the top width of the waveguide
                    Parameters['x res'] : int/float
                          Minimum mesh step setting in x direction. (Size of the Mesh Box in x direction in um)
                    Parameters['Slab Height'] : int/float
                        Height of the LiNbO3 Slab
                    Parameters['Wavelength'] : int/float
                        Wavelength
                    Parameters["Waveguide Angle"] : int/float
                        This Parameter will set the theta ratation angle of the port
                    Parameters["Port Span"] : list of floats/ints
                        List of x,y and z span of the Ports. For this simulation only y and z parametes will be taken.
                    Parameters["Mode"] : str
                         Mode to analyse. ("fundamental TE mode", "fundamental TM mode", "fundamental mode")
                                         
                    
        
                  10) Solver("ArcWaveguide", "FDTD", Parameters)
                               
                     Parameters['Substrate Height'] : int/float
                         Substrate height.
                     Parameters['WG Height'] : int/float
                         Waveguide hight
                     Parameters['WG Width'] : int/float
                         Top width of the Waveguide.
                     Parameters['x res'] : int/float
                         Minimum mesh step setting in x direction. (Size of the Mesh Box in x direction in um)
                     Parameters['Slab Height'] : int/float
                         Height of the LiNbO3 Slab
                     Parameters['Wavelength'] : int/float
                         Wavelength
                     Parameters["S_Band Radius"] : int/float
                         S-Bend Radius in um.
                     Parameters['Arc deg'] : int/float
                         Arc define the Arc of the curve. It can be 90 or 180 degrees only.
                         This two will define an 1/4 of a circle or 1/2 of a circle.
                     Parameters["Mode"] : str
                          Mode to choose from ("fundamental TE mode", "fundamental TM mode", "fundamental mode")
                                                    
                    
        
                  11) Solver("BendWaveguide", "FDTD", Parameters)    
                               
                    Parameters['Substrate Height'] : int/float
                        Substrate height
                    Parameters['WG Height'] : int/float
                        Waveguide hight
                    Parameters['WG Width'] : int/float
                        Top width of the Waveguide
                    Parameters['x res'] : int/float
                        Minimum mesh step setting in x direction. (Size of the Mesh Box in x direction in um)
                    Parameters['Slab Height'] : int/float
                        Height of the LiNbO3 Slab
                    Parameters['Wavelength'] : int/float
                        Wavelength
                    Parameters["x span"] : int/float
                        Length of the S-Bend Waveguide
                    Parameters["y span"] : int/float
                        Width of the S-Bend Waveguide
                    Parameters["Mode"] : str
                         Mode to choose from ("fundamental TE mode", "fundamental TM mode", "fundamental mode")
                    Parameters["Port Span"] : list of floats/ints
                        List of x,y and z span of the Ports. For this simulation only y and z parametes will be taken, since we use x cut material.
                        
                   
                               
        
              """)
              
    def Help_StartSimulation(self):
        
        print("""
              
              Depending on the Solver Type ("EME" or "FDTD") simulation can be started with:
                  
                  - For "EME" Solver 
                  1) StartFDESimulation() 
                  2) StartEMESimulation()
                           
                  - For "FDTD" Solver
                  3) StartFDTDSimulation()
             
        
              """)
              
              
              
              
    def Help_Results(self):
        
        print("""
             Data can be Extracted from an FDE, EME or FDTD simulation. The Data will be writen in dictionaries!
             
             
    ############################################   FDE   ############################################
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
                 
                
                
                
                 
                 - Extract Results from Simulation.
                 
                 With - ExtrtactData('FDE', Parameters)
                 
                 
                 Parameters["Effective Index"] : float/int
                     Effective index value. All modes with effective index smaller then
                     the  Effective will be deleted. Effective is usualy the smalles effective
                     index from the hole construction (usualy the cladding). Since Lumerical will produce 
                     modes with very small Effektive index, it make no sance to extract them !
                     
                 Returns
                 -------
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
                    
    ############################################   FDE   ############################################
    
    
    
    
    
    ############################################   EME   ############################################
    
             After an EME Simulation the user can:
                 
                 
             - Extract all the data providet from the simulation.
             
             
               With - ExtrtactData('EME', Parameters)  
               
               Parameters here is an dummy variable, since you already give all the needed Parameter by the 
               object and solver. They will be automatically taken for the result extraction.
    
               Returns
               -------
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
                           
                           
    ############################################   EME   ############################################
    
    
    ############################################   FDTD   ############################################            
    
              After an FDTD Simulation the user can: 
            
            
                With this fuction will add an S-Parameter Sweep to the FDTD simulation.
                This will allowed to sweep the S-parameter Matrix like it is done in
                EMe solver. Additionally will extract the Power on each port according to the time
                monitors that ware set in the solver functions.
                
            With - ExtrtactData('FDTD', Parameters)
            
            Here Paraameters is an boolen operator!!
            Parameter = True - An S-Parameter extraction sweep will be added to the simulation. This will increase the simulation time, 
            since the S-Parameters need to be propagated thrue the strucutre. 
            Parameter = False - No S-Parameter extraction sweep. The dictionary will return 0. 
            
            Returns
            -------
             SParam_Dict: dict
                Dictionary with:
                        "S Parameter"
                    
             PowerData_Dict : dict
                Dictionary with:
                        'Power Input Port/s'
                        'Power Output Port/s'
                        'Transmission'
                        
                        
    ############################################   FDTD   ############################################ 
    
    
              """)
              
              
            
              
              
              
              
    def Help_LoadingBar(self):
        
        print("""
              
              If the user want to keep track of the simulation state an laoding bar is imported in this library.
              The loading bar is not part of the class, so you dont need to call it like an class. 
              
              loadingBar(count, total)
              count - is your loop variable
              total - is the lenght of the variables that will be sweeped
              
              
              Example python code of loading bar:
                  
                  from Constructor import loadingBar

                  for i in range(10):
                     loadingBar(i, 10)
                  

              """)

    def Help_RemoveObject(self):

        print("""
        
                To Remove an Object simply use:
                
                removeObject()
        
        """)
        
        
        
    def Help_LogFile(self):
        print(
            """
            When you load the Logfile function from the Constructor Library the user can keep track 
            with the simulation and parameters used. To use the log File you need to give an list with 
            two dictionarays. The first can be the dictionary with Parameters given for the simulation 
            like "Tapper Width", "WG Width"... The secound dictionary can be called direct after performing 
            an analysis with obj.ReturnLogInfo(). This will extract the object that you simulate and the 
            solver that you used. An text file will be created in the given Path. The name of the File 
            will be "Lumerical_Simulation_Log_15_10_52"
            
            An example can be seen here:
                
                from Constructor import Logfile
                Path = "C:/Downloads/"
                solverInfo = obj.ReturnLogInfo()
                data = [solverInfo, Parameters]
                Logfile(data, Path)
                
            
            """
            )


def Help():
    print("""
    
    To start the Python Lumerical Constructor Library you need first to:
        from Constructor import Constructor
        from Constructor import loadingBar
        
        # Call the Constructor Object 
        obj = Constructor(file, MaterialPath, Lumerical Solver Type)
            # file = Path to lumapi.py for example -> "lumerical/2022-r2.4/api/python/lumapi.py"
            # MaterialPath = Path to Material Library for example -> "/Desktop/LNOI_Env/LiNbO3.mdf"
            # Lumerical Solver Type = Solver Type can be EME or FDTD for example -> "EME"
            
        # To Close Lumerical simply use
        obj.Close()
        
        # For Further Help after colling the Constructor Object please use:
            # obj.Help('Objects') -> Shows how to build an photonic element
            # obj.Help('Solvers') -> Shows hot to build an FDE, EME or FDTD Solver
            # obj.Help('Start Simulation') -> Commands to start simulation
            # obj.Help('Results') -> Shows how to extract resuslts form FDE, EME and FDTD Simulations
            # obj.Help('Loading Bar') -> Shows how to use the Loading Bar Function 
            # obj.Help('Remove Object') -> Romove object function
            # obj.Help('Log File') -> Generate an log file from the simulation 
    """)
    




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
    
    