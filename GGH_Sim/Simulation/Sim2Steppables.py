from PlayerPython import *
import CompuCellSetup

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
P_sr = 0.3                  # symetric self renewing
P_ar = 0.5                  # asymetric self renewing
P_sd = 1 - (P_sr + P_ar)    # symetric differentiating

# Critical size to trigger mitosis
cell_critical_volume = 50
targetVolume = 25
lambdaVolume = 10  # E_vol = lambdaVolume(V - V_cell)^2

# </parameter settings>

# Export parameters to a json file for interoperability

params = {'growth_rate': growth_rate,
          'P_ar': P_ar,
          'P_sr': P_sr,
          'P_sd': P_sd,
          'cell_critical_volume': cell_critical_volume,
          'targetVolume': targetVolume,
          'lambdaVolume': lambdaVolume}

json_fname = 'json_params.json'

with file(json_fname, 'w+') as json:
    json.dump(params)


class ConstraintInitializerSteppable(SteppableBasePy):
    ''' Class used to initialize the cells

    '''
    def __init__(self, _simulator, _frequency=1):
        SteppableBasePy.__init__(self, _simulator, _frequency)

    def start(self):

        for cell in self.cellList:
            cell.targetVolume = targetVolume
            cell.lambdaVolume = lambdaVolume
            cellDict = self.getDictionaryAttribute(cell)
            cellDict["age"] = 0


class GrowthSteppable(SteppableBasePy):
    ''' Class governing cell growth steps.

    '''

    def __init__(self, _simulator, _frequency=1):

        SteppableBasePy.__init__(self, _simulator, _frequency)

    def step(self, mcs):
        for cell in self.cellList:
            # here is the random growth rate implementation
            dice = np.random.random()
            if dice > growth_rate:
                continue
            if cell.type == 1:
                cell.targetVolume += 1

            elif cell.type == 2:
                cellDict = self.getDictionaryAttribute(cell)
                if cellDict['age'] < prolif_potential:
                    cell.targetVolume += 1


class MitosisSteppable(MitosisSteppableBase):
    ''' class defining the mitosis steps'''

    def __init__(self, _simulator, _frequency=10):

        MitosisSteppableBase.__init__(self, _simulator, _frequency)
        self.ages = self.createScalarFieldCellLevelPy("IdField")

    def step(self, mcs):
        # print "INSIDE MITOSIS STEPPABLE"
        cells_to_divide = []

        for cell in self.cellList:
            if cell.volume < cell_critical_volume:
                continue
            elif cell.type == 2 and self.ages[cell] > prolif_potential:
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

        self.ages[parentCell] += 1
        self.ages[childCell] = self.ages[parentCell]

        p_cellDict = self.getDictionaryAttribute(parentCell)
        p_cellDict['age'] = self.ages[parentCell]

        c_cellDict = self.getDictionaryAttribute(childCell)
        c_cellDict['age'] = self.ages[childCell]

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
            childCell.type = 2
