from PlayerPython import *
import CompuCellSetup

from PySteppables import *
from PySteppables import SteppableBasePy
import CompuCell
import sys
import numpy as np
import json

from PySteppablesExamples import MitosisSteppableBase

# TODO
# include apoptosis (not quantified)
# track cells neighbourhood

# *** All the global variables must be set bellow
# Prameters set here prevail over

# <parameter settings>
params = {
    'growth_rate': 0.2,
    'P_sr': 0.8,
    'cell_critical_volume': 50,
    'targetVolume': 25,
    'lambdaVolume': 10,
    'prolif_potential': 4,
    'neighbor_dep_after': False,
    'neighbor_dep_before': False,
    }
# </parameter settings>

# Randomize cell growth
# This adds a probabilistic penalty to the growth process.
# if growth_rate is 1: deterministic growth
# if growth_rate is 1/n: cell will grow every n steps on average
growth_rate = params['growth_rate']  # random  (uniform in [0, 1[)
# include limited number of differentiation ( 4 cycles) for NCPs
prolif_potential = params['prolif_potential']  # maximum number of divisions

# Self renewing probability, if both 'neighbor_dep_after' and
#     'neighbor_dep_after' are False
P_sr0 = params['P_sr']

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
            cellDict["P_sr"] = params['P_sr']

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
            if cell.type == 1: # CancerStemCells
                cell.targetVolume += 1

            elif cell.type == 2: # Non prolif
                cellDict = self.getDictionaryAttribute(cell)
                if cellDict['age'] < prolif_potential:
                    cell.targetVolume += 1


class MitosisSteppable(MitosisSteppableBase):
    ''' class defining the mitosis steps'''

    def __init__(self, _simulator, _frequency=10):

        MitosisSteppableBase.__init__(self, _simulator, _frequency)

    def step(self, mcs):
        """Mitotic cells selection  performed at each step of the simulation
        """
        # print "INSIDE MITOSIS STEPPABLE"
        cells_to_divide = []

        for cell in self.cellList:
            cellDict = self.getDictionaryAttribute(cell)
            if cell.volume < cell_critical_volume:
                continue # jumps to next cell
            elif cell.type == 2 and cellDict['age'] > prolif_potential:
                continue
            cells_to_divide.append(cell)

            for cell in cells_to_divide:
                if params['neighbor_dep_before'] and cell.type == 1:
                    self.set_renewal_prob(cell)

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

        p_cellDict = self.getDictionaryAttribute(parentCell)
        p_cellDict['age'] = p_cellDict['age'] + 1

        c_cellDict = self.getDictionaryAttribute(childCell)
        c_cellDict['age'] = p_cellDict['age']
        c_cellDict['P_sr'] = p_cellDict['P_sr']

        if parentCell.type == 1:
            # TODO test this before ddiv / after div
            if params['neighbor_dep_after']:
                self.set_renewal_prob(parentCell)
                self.set_renewal_prob(childCell)
            self.choose_type(parentCell)
            self.choose_type(childCell)

        else: # Non proliferative cell don't differentiate
            c_cellDict['P_sr'] = p_cellDict['P_sr'] = 0
            parentCell.type = 2
            childCell.type = 2

    def get_neighbour_types(self, cell):
        n_types = []
        for neighbor, c_area in self.getCellNeighborDataList(cell):
            if neighbor is None:
                n_types.append(2)
            else:
                n_types.append(neighbor.type)
        return np.array(n_types)

    def set_renewal_prob(self, cell):
        """
        """
        cellDict = self.getDictionaryAttribute(cell)
        parent_types = self.get_neighbour_types(cell)
        if parent_types.size > 0:
            n_t1 = np.float((parent_types == 1).sum())
            P_sr = (n_t1 / parent_types.size)
            cellDict['P_sr'] = P_sr
        else:
            cellDict['P_sr'] = 0.5

    def choose_type(self, cell, P_sr=None):

        P_sr = self.getDictionaryAttribute(cell)['P_sr']

        dice = np.random.random() # random number between 0 and 1
        if dice < P_sr:
            cell.type = 1
        else:
            cell.type = 2

class RegisterSteppable(SteppableBasePy):
    ''' Class governing cell growth steps.

    '''

    def __init__(self, _simulator, _frequency=1):

        SteppableBasePy.__init__(self, _simulator, _frequency)

    def start(self):
        self.filename = 'sim_data.csv'
        self.data = []
        self.columns = ['mcs', 'id',
                        'type', 'age', 'P_sr',
                        'n_type1', 'n_type2', # Number of neighbors by types
                        'a_type1', 'a_type2'] # Total contact area by types



    def step(self, mcs):
        for cell in self.cellList:
            n_type1 = 0.0
            n_type2 = 0.0
            a_type1 = 0.0
            a_type2 = 0.0

            for neighbor, c_area in self.getCellNeighborDataList(cell):
                if neighbor is None:
                    continue
                if neighbor.type == 1:
                    n_type1 += 1
                    a_type1 += c_area
                elif neighbor.type == 2:
                    n_type2 += 1
                    a_type2 += c_area

            cellDict = self.getDictionaryAttribute(cell)

            step_data = [mcs, cell.id,
                         cell.type,
                         cellDict['age'],
                         cellDict['P_sr'],
                         n_type1, n_type2,
                         a_type1, a_type2]
            self.data.append(step_data)

    def finish(self):

        file_handle, full_fname = \
          self.openFileInSimulationOutputDirectory(self.filename, "w+")
        file_handle.close()

        metadata = full_fname
        col_names = ','.join(self.columns)
        header = '\n'.join([metadata, col_names])
        np.savetxt(full_fname, np.asarray(self.data), header=col_names, delimiter=',')
        print('saved sim data to {}'.format(full_fname))
