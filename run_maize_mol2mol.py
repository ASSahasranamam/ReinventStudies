import time
from pathlib import Path
from maize.core.workflow import Workflow
from maize.steps.io import LoadData, LogResult, Return, Void
from maize.steps.mai.docking.adv import AutoDockGPU
from maize.steps.mai.molecule import Gypsum, LoadMolecule
from maize.utilities.chem import IsomerCollection
from maize.steps.mai.misc import ReInvent
from maize.steps.mai.cheminformatics import RMSD, ExtractScores, TagIndex, SortByTag, TagSorter, LogTags
from maize.steps.plumbing import MergeLists

import numpy as np
from numpy.typing import NDArray

from maize.core.interface import Input
from maize.core.node import Node
import json


class ScoreLog(Node):
    """Logs scores in the form of NDArrays"""

    inp: Input[NDArray[np.float32]] = Input()

    def run(self) -> None:
        scores = self.inp.receive()
        # Format the output as expected by REINVENT
        output = {
            "version": 1,
            "payload": {
                "score": scores.tolist()  # Convert numpy array to list
            }
        }
        (json.dumps(output))  # Print JSON formatted output
        self.logger.info("Received scores: %s", scores)


flow = Workflow(name="dock", level="debug", cleanup_temp=False)
flow.config.update(Path("configs/Maize/maize-mol2mol-config.toml"))

rnve = flow.add(ReInvent)
embe = flow.add(Gypsum, loop=True)
indx = flow.add(TagIndex, loop=True)
dock = flow.add(AutoDockGPU, loop=True)
void = flow.add(Void)
load = flow.add(LoadMolecule)
rmsd = flow.add(RMSD, loop=True)
logt = flow.add(LogTags, loop=True)
# sort = flow.add(TagSorter, loop=True)
# dock_hp = flow.add(AutoDockGPU, name="dock-hp", loop=True)
# void_hp = flow.add(Void, name="void-hp", loop=True)
# merg = flow.add(MergeLists[IsomerCollection])
# sort_id = flow.add(SortByTag, loop=True)
scor = flow.add(ExtractScores, loop=True)
# save = flow.add(ScoreLog)
flow.connect_all(
    (rnve.out, embe.inp),
    (embe.out, indx.inp),
    (indx.out, dock.inp),
    (dock.out_scores, void.inp),
    (dock.out, rmsd.inp),
    (load.out, rmsd.inp_ref),
)
flow.connect_all(
    (rmsd.out, logt.inp),
    (logt.out, scor.inp),
    (scor.out, rnve.inp),
)

grid = Path("mols/Rad51/rad51_receptor.maps.fld")
ref = Path("mols/Rad51/Cam833HMdsRad51Docked_Cam833-Acid_3.sdf")
rnv_config = Path("configs/REINVENT/staged_learning_maize_mol2mol.toml")
prior = Path("priors/mol2mol_similarity.prior")

# Make sure these files exist and are readable
assert rnv_config.exists(), f"Config file not found: {rnv_config}"
assert prior.exists(), f"Prior model not found: {prior}"

# grid = Path("../maize/steps/mai/docking/data/1uyd.tar")


rnve.configuration.set(rnv_config)

rnve.prior.set(prior)
rnve.agent.set(prior)
# rnve.args.set("----loglevel DEBUG")

# The maximum number of RL epochs
rnve.max_epoch.set(3)
# rnve.max_epoch.set(3)

# The REINVENT configuration, excluding any entries for maize (these will be added automatically)
rnve.configuration.set(rnv_config)
rnve.prior.set(prior)
rnve.agent.set(prior)
# rnve.inp_smi.set(100)

rnve.distance_threshold.set(100)

# Settings to transform the docking score to a value between 0 and 1, with 1 being favourable, using a sigmoid
rnve.low.set(-10.0)
rnve.high.set(-5.0)
rnve.reverse.set(True)

# Number of molecules to generate each epoch
rnve.batch_size.set(128)

# Number of isomers to generate for each SMILES
embe.use_filters.set(False)
embe.n_variants.set(4)


# Docking grid for 1UYD
dock.grid_file.set(grid)
# dock_hp.grid_file.set(grid)

# Reference ligand for RMSD calculation
load.path.set(ref)

# Log the "rmsd" tag
logt.tag.set("rmsd")

# Send molecules with RMSD higher than 6 Ang to high-precision docking
# sort.sorter.set(["rmsd < 6.0", "rmsd >= 6.0"])

# More extensive search settings for "high-precision" docking
# dock_hp.nrun.set(50)
# dock_hp.population_size.set(300)
# dock_hp.lsit.set(500)

# Deactivate constraints from the grid
dock.constraints.set(False)
# dock_hp.constraints.set(False)

flow.check()
flow.execute()
# mols = retu.get()
# print(mols)
