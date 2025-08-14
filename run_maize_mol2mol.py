import time
from pathlib import Path
from maize.core.workflow import Workflow
from maize.steps.io import LoadData, LogResult, Return, Void
from maize.steps.mai.docking.adv import AutoDockGPU
from maize.steps.mai.molecule import Gypsum, LoadMolecule
from maize.utilities.chem import IsomerCollection, save_sdf_library
from maize.utilities.io import setup_workflow
from maize.steps.mai.misc import ReInvent
from maize.steps.mai.cheminformatics import RMSD, ExtractScores, TagIndex, SortByTag, TagSorter, LogTags
from maize.steps.plumbing import MergeLists

import numpy as np
from numpy.typing import NDArray

from maize.core.interface import Input, Output, Parameter
from maize.core.node import Node
import json

# NEW: imports for saving artifacts
import shutil
from datetime import datetime
from typing import Any

totalIsomerCollection: list[IsomerCollection] = []
# NEW: Pass-through node that writes docked structures to an SDF library

class SaveDockSDF(Node):
    """Saves the incoming docked payload to a timestamped SDF library and forwards it unchanged."""

    inp: Input[list[IsomerCollection]] = Input()
    out: Output[list[IsomerCollection]] = Output(optional=True)
    payload: Parameter[list[IsomerCollection]] = Parameter(optional=True)
    def run(self) -> None:
        self.payload = self.inp.receive()
        ts = time.strftime("%Y%m%d-%H%M%S")
        print(self.payload)

        out_loc_name = f"dock_outputs_{ts}.sdf"
        out_loc = Path(out_loc_name)

        totalIsomerCollection.extend(self.payload)
        print(totalIsomerCollection)
        ts = time.strftime("%Y%m%d-%H%M%S")

        try:
            save_sdf_library(mols=totalIsomerCollection, file=Path(f"dock_outputs_{ts}.sdf"))
            self.logger.info(Path(f"dock_outputs_{ts}.sdf"))
            if Path(f"dock_outputs_{ts}.sdf").is_file():
                self.logger.info(Path(f"dock_outputs_{ts}.sdf").absolute())
                self.logger.info(Path(f"dock_outputs_{ts}.sdf").resolve())
                shutil.copy(f"dock_outputs_{ts}.sdf", Path("/home/adi_sahasranamam/"))

                self.logger.info(f"dock_outputs_{ts}.sdf moved to {out_loc}")
        except Exception as e:
            self.logger.error(e)
        self.out.send(self.payload)

flow = Workflow(name="dock", level="debug", cleanup_temp=False)
flow.config.update(Path("configs/Maize/maize-mol2mol-config.toml"))

rnve = flow.add(ReInvent)
embe = flow.add(Gypsum, loop=True)
indx = flow.add(TagIndex, loop=True)
dock = flow.add(AutoDockGPU, loop=True, cleanup_temp=False)
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

# NEW: add SDF saver between docking and RMSD
saver = flow.add(SaveDockSDF, loop=True,cleanup_temp=False)

flow.connect_all(
    (rnve.out, embe.inp),
    (embe.out, indx.inp),
    (indx.out, dock.inp),
    (dock.out_scores, void.inp),
    # Modified: dock.out -> saver.inp -> rmsd.inp
    (dock.out, saver.inp),
    (saver.out, rmsd.inp),
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
rnve.distance_threshold.set(100)
rnve.sample_strategy.set("beamsearch")  # multinomial or beamsearch (deterministic)
rnve.input_smi.set("/home/adi_sahasranamam/workspace/reinventstudies/mols/Rad51/bdb_2.smi")  # multinomial or beamsearch (deterministic)

# Make sure these files exist and are readable
assert rnv_config.exists(), f"Config file not found: {rnv_config}"
assert prior.exists(), f"Prior model not found: {prior}"

# grid = Path("../maize/steps/mai/docking/data/1uyd.tar")


rnve.configuration.set(rnv_config)

rnve.prior.set(prior)
rnve.agent.set(prior)
# rnve.args.set("----loglevel DEBUG")

# The maximum number of RL epochs
# rnve.min_epoch.set(25)
rnve.max_epoch.set(1)
# rnve.max_score.set(0.85)



# Settings to transform the docking score to a value between 0 and 1, with 1 being favourable, using a sigmoid
rnve.low.set(-10.0)
rnve.high.set(-5.0)
rnve.reverse.set(True)

# Number of molecules to generate each epoch
rnve.batch_size.set(4)

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
# try:


flow.check()



setup_workflow(flow)


ts = time.strftime("%Y%m%d-%H%M%S")


print("HII, ", len(totalIsomerCollection))
print((totalIsomerCollection))

# save_sdf_library(mols=totalIsomerCollection,file=Path( f"dock_outputs_{ts}.sdf"))
# except Exception as e:
    # flow.logger.error(e)
    # flow._cleanup()

# mols = retu.get()
# print(mols)
