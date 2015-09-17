
import sys
from os import environ
from os import getcwd
import string

import CompuCellSetup


sys.path.append(environ["PYTHON_MODULE_PATH"])


sim, simthread = CompuCellSetup.getCoreSimulationObjects()

# add extra attributes here

CompuCellSetup.initializeSimulationObjects(sim, simthread)
# Definitions of additional Python-managed fields go here

# ### Python Steppables ####
steppableRegistry = CompuCellSetup.getSteppableRegistry()

# ## Constraint
from Sim2Steppables import ConstraintInitializerSteppable
ConstraintInitializerSteppableInstance = ConstraintInitializerSteppable(sim,_frequency=1)
steppableRegistry.registerSteppable(ConstraintInitializerSteppableInstance)

## Growth
from Sim2Steppables import GrowthSteppable
GrowthSteppableInstance = GrowthSteppable(sim,_frequency=1)
steppableRegistry.registerSteppable(GrowthSteppableInstance)

## Mitosis
from Sim2Steppables import MitosisSteppable
MitosisSteppableInstance = MitosisSteppable(sim,_frequency=1)
steppableRegistry.registerSteppable(MitosisSteppableInstance)

CompuCellSetup.mainLoop(sim, simthread, steppableRegistry)
