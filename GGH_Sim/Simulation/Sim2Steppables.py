
from PySteppables import *
import CompuCell
import sys
import numpy as np
from PySteppablesExamples import MitosisSteppableBase


### TODO:
# include limited number of differentiation ( 4 cycles) for NCPs
# prolif_potential = 5
# include apoptosis (not quantified)
# track cells neighbourhood



### Limits the cell growth
growth_rate = 0.1


P_sr = 0.4 # symetric self renewing
P_ar = 0.4 # asymetric self renewing
P_sd = 1 - (P_sr + P_ar) # symetric differentiating

class ConstraintInitializerSteppable(SteppableBasePy):
    ''' Class used to initialize the cells

    '''
    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
    def start(self):
        for cell in self.cellList:
            cell.targetVolume=25
            cell.lambdaVolume=2.0

class GrowthSteppable(SteppableBasePy):
    ''' Class governing cell growth steps.

    '''
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


class MitosisSteppable(MitosisSteppableBase):
    ''' class defining the mitosis steps'''
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
            # self.divideCellOrientationVectorBased(cell,1,0,0)
            # self.divideCellAlongMajorAxis(cell)
            # self.divideCellAlongMinorAxis(cell)

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
