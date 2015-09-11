
from PySteppables import *
import CompuCell
import sys
import numpy as np
from PySteppablesExamples import MitosisSteppableBase

# TODO
# include apoptosis (not quantified)
# track cells neighbourhood

# *** All the global variables must be set bellow
# <parameter settings>

# Randomize cell growth
# This adds a probabilistic penalty to the growth process.
# if growth_rate is 1: deterministic growth
# if growth_rate is 1/n: cell will grow every n steps on average
growth_rate = 1.  # random process probability (uniform in [0, 1[)

# include limited number of differentiation ( 4 cycles) for NCPs
prolif_potential = 4  # maximum number of divisions (mother cell)

# Differentiation probabilities
P_sr = 0.4                  # symetric self renewing
P_ar = 0.4                  # asymetric self renewing
P_sd = 1 - (P_sr + P_ar)    # symetric differentiating

# Critical size to trigger mitosis
cell_critical_volume = 50

# </parameter setting>


class ConstraintInitializerSteppable(SteppableBasePy):
    ''' Class used to initialize the cells

    '''
    def __init__(self, _simulator, _frequency=1):
        SteppableBasePy.__init__(self, _simulator, _frequency)

    def start(self):

        for cell in self.cellList:
            cell.targetVolume = 25
            cell.lambdaVolume = 2.0
            cell.age = 0

class GrowthSteppable(SteppableBasePy):
    ''' Class governing cell growth steps.

    '''

    def __init__(self,_simulator,_frequency=1):

        SteppableBasePy.__init__(self,_simulator,_frequency)

    def step(self,mcs):
        for cell in self.cellList:
            # here is the random growth rate implementation
            dice = np.random.random()
            if dice > growth_rate:
                continue
            if cell.type == 1:
                cell.targetVolume += 1

            elif cell.type == 2:
                cell.targetVolume += 1


class CellGenerationSteppable(IdFieldVisualizationSteppable):
    ''' Cell level scalar field following the number of generations
        a cell has seen. It goes this way:

        cell0 cell1 cell2
        1
        2     1
        3     2     1
        4     3     2
    '''

    def __init__(self, _simulator, _frequency=1):

        SteppableBasePy.__init__(self, _simulator, _frequency)
        self.scalarCLField = self.createScalarFieldCellLevelPy("IdField")

    def step(self,mcs):
        self.scalarCLField.clear()

        for cell in self.cellList:
            self.scalarCLField[cell]=cell.id*random()


class MitosisSteppable(MitosisSteppableBase):
    ''' class defining the mitosis steps'''

    def __init__(self, _simulator, _frequency=10):

        MitosisSteppableBase.__init__(self, _simulator, _frequency)

    def step(self, mcs):
        # print "INSIDE MITOSIS STEPPABLE"
        cells_to_divide = []
        for cell in self.cellList:
            if cell.volume < cell_critical_volume:
                continue
            elif cell.type == 2 and age[cell] > prolif_potential:
                continue
            cells_to_divide.append(cell)

        for cell in cells_to_divide:
            # to change mitosis mode leave one of the below lines uncommented
            self.divideCellRandomOrientation(cell)
            # self.divideCellAlongMajorAxis(cell)
            # self.divideCellOrientationVectorBased(cell,1,0,0)
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
