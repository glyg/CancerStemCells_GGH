from PlayerPython import *
import CompuCellSetup

from PySteppables import *
import CompuCell
import sys
import numpy as np
import json

from PySteppablesExamples import MitosisSteppableBase

# TODO
# include apoptosis (not quantified)
# track cells neighbourhood

# *** All the global variables must be set bellow
# <parameter settings>
params = {
    'growth_rate': 0.2,
    'P_sr': 0.5,
    'P_ar': 0.4,
    'cell_critical_volume': 50,
    'targetVolume': 25,
    'lambdaVolume': 10,
    'prolif_potential': 4,
    'neighbor_dep': False
    }
# </parameter settings>

# Randomize cell growth
# This adds a probabilistic penalty to the growth process.
# if growth_rate is 1: deterministic growth
# if growth_rate is 1/n: cell will grow every n steps on average
growth_rate = params['growth_rate']  # random  (uniform in [0, 1[)
# include limited number of differentiation ( 4 cycles) for NCPs
prolif_potential = params['prolif_potential']  # maximum number of divisions
# Differentiation probabilities
P_sr = params['P_sr']                 # symetric self renewing
P_ar = params['P_ar']                 # asymetric self renewing
P_sd = 1 - P_sr - P_ar    # symetric differentiating

# Critical size to trigger mitosis
cell_critical_volume = params['cell_critical_volume']

targetVolume = params['targetVolume']  # E_vol = V_\lambda(V - V_target)^2
lambdaVolume = params['lambdaVolume']  # (K in the Vertex model)

# TODO: Export parameters to a json file for interoperability
# This raises a IOError, investigate
# json_fname = 'json_params.json'


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
        self.ages = self.createScalarFieldCellLevelPy("CellAge")

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
            # TODO test this before ddiv / after div
            _P_sr, _P_ar = P_sr, P_ar
            if params['neighbor_dep']:
                n_types = [neighbor.type for neighbor, c_area
                           in self.getCellNeighborDataList(parentCell)
                           if neighbor is not None]
                n_types = np.array(n_types)
                if n_types.size > 0:
                    n_t1 = np.float((n_types == 1).sum())
                    _P_sr = (n_t1 / n_types.size)  # * P_sr / (1 + P_ar)
                    # _P_sr = max(_P_sr, P_sr)
                    # _P_ar = (n_t1 / n_types.size) * P_ar
                    print(_P_sr)
            dice = np.random.random()
            if dice < _P_sr:
                childCell.type = 1
            elif _P_sr <= dice < _P_ar + _P_sr:
                childCell.type = 2
            elif dice >= _P_ar + _P_sr:
                childCell.type = 2
                parentCell.type = 2
        else:
            childCell.type = 2
