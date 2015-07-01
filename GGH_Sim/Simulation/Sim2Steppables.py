
from PySteppables import *
import CompuCell
import sys

import numpy as np

from PySteppablesExamples import MitosisSteppableBase


class ConstraintInitializerSteppable(SteppableBasePy):
    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
    def start(self):
        for cell in self.cellList:
            cell.targetVolume=25
            cell.lambdaVolume=2.0


growth_rate = 0.1

# include limited number of differentiation ( 4 cycles)
prolif_potential = 5
# include apoptosis (not quantified)
# track cells neighbourhood

class GrowthSteppable(SteppableBasePy):
    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
    def step(self,mcs):
        for cell in self.cellList:
            dice = np.random.random()
            if dice > growth_rate:
                continue
            if cell.type == 1:
                cell.targetVolume += 1
            elif cell.type == 2:
                if mcs % 100 == 0:
                    cell.targetVolume+=1
    # alternatively if you want to make growth a function of chemical concentration uncomment lines below and comment lines above
        # field=CompuCell.getConcentrationField(self.simulator,"PUT_NAME_OF_CHEMICAL_FIELD_HERE")
        # pt=CompuCell.Point3D()
        # for cell in self.cellList:
            # pt.x=int(cell.xCOM)
            # pt.y=int(cell.yCOM)
            # pt.z=int(cell.zCOM)
            # concentrationAtCOM=field.get(pt)
            # cell.targetVolume+=0.01*concentrationAtCOM  # you can use here any fcn of concentrationAtCOM



P_sr = 0.4 # symetric self renewing
P_ar = 0.4 # asymetric self renewing
P_sd = 1 - (P_sr + P_ar) # symetric differentiating

class MitosisSteppable(MitosisSteppableBase):
    def __init__(self,_simulator,_frequency=10):
        MitosisSteppableBase.__init__(self,_simulator, _frequency)

    def step(self,mcs):
        # print "INSIDE MITOSIS STEPPABLE"
        cells_to_divide=[]
        for cell in self.cellList:
            if cell.volume>50:
                cells_to_divide.append(cell)

        for cell in cells_to_divide:
            # to change mitosis mode leave one of the below lines uncommented
            self.divideCellRandomOrientation(cell)
            # self.divideCellOrientationVectorBased(cell,1,0,0)                 # this is a valid option
            # self.divideCellAlongMajorAxis(cell)                               # this is a valid option
            # self.divideCellAlongMinorAxis(cell)                               # this is a valid option

    def updateAttributes(self):
        parentCell = self.mitosisSteppable.parentCell
        childCell = self.mitosisSteppable.childCell
        parentCell.targetVolume = parentCell.targetVolume / 2
        childCell.targetVolume = parentCell.targetVolume
        childCell.lambdaVolume = parentCell.lambdaVolume
        if parentCell.type == 1:
            dice = np.random.random()
            if dice < P_sr:
                childCell.type = 1
            elif P_sr <= dice < P_ar + P_sr:
                childCell.type = 2
            elif dice >= P_ar + P_sr:
                childCell.type = 2
                parentCell.type = 2
        else:
            childCell.type=2
