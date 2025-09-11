import time
import tomlkit
from pathlib import Path
from rdkit import Chem
from rdkit.Chem.ChemUtils.SDFToCSV import Convert
from rdkit.Chem import PandasTools
from rdkit import RDConfig
from exceptiongroup import catch
from maize.core.workflow import Workflow
from maize.steps.io import LoadData, LogResult, Return, Void
from maize.steps.mai.docking.adv import AutoDockGPU,Vina
from maize.steps.mai.molecule import Gypsum, LoadMolecule,SaveCSV,IsomerCollectionSaving
from maize.utilities.chem import IsomerCollection, save_sdf_library
from maize.utilities.io import setup_workflow
from maize.steps.mai.misc import ReInvent
from maize.steps.mai.cheminformatics import RMSD, ExtractScores, TagIndex, SortByTag, TagSorter, LogTags
from maize.steps.plumbing import MergeLists
import os
from maize.core.interface import Flag, FileParameter, Parameter, Suffix, Input, Output
import numpy as np
from numpy.typing import NDArray
from maize.steps.mai.misc import expose_reinvent
from maize.core.interface import Input, Output, Parameter
from maize.core.node import Node
import json
from  maize.graphs.mai.dock import Docking
from maize.core.graph import Graph
# NEW: imports for saving artifacts
import shutil
from datetime import datetime
from typing import Any
from maize.steps.mai.molecule import (
    Gypsum,
    File2Molecule,
    SaveMolecule,
    SaveScores,
    LoadSmiles,
    SaveLibrary,
    LoadMolecule,
    LoadLibrary,
    Ligprep,
    CombineMolecules,
    AggregateScores,
    SaveCSV,
)
totalIsomerCollection: list[IsomerCollection] = []
# NEW: Pass-through node that writes docked structures to an SDF library

import tomli as tomllib

class SaveDockSDF(Node):
    """Saves the incoming docked payload to a timestamped SDF library and forwards it unchanged."""

    inp: Input[list[IsomerCollection]] = Input()
    out: Output[list[IsomerCollection]] = Output()
    out2: Output[list[IsomerCollection]] = Output()

    payload: Parameter[list[IsomerCollection]] = Parameter(optional=True)
    def run(self) -> None:
        self.payload = self.inp.receive()
        ts = time.strftime("%m%d-%H%M%S")
        print(self.payload)

        # wd = f"/tmp/MC_R4_notebooks_output_{ts}"

        # ### Delete existing working directory and create a new one
        #
        # If the working directory already exists, it will be reused

        # +
        # if not os.path.isdir(wd):
        #     shutil.rmtree(wd, ignore_errors=True)
        #     os.mkdir(wd)

        out_loc_name = f"dock_outputs_{ts}.sdf"
        out_loc = Path(out_loc_name)

        totalIsomerCollection.extend(self.payload)
        print(totalIsomerCollection)

        try:
            save_sdf_library(mols=totalIsomerCollection, file=Path(f"dock_outputs_{ts}.sdf"))
            self.logger.info(Path(f"dock_outputs_{ts}.sdf"))
            if Path(f"dock_outputs_{ts}.sdf").is_file():
                self.logger.info(Path(f"dock_outputs_{ts}.sdf").absolute())
                self.logger.info(Path(f"dock_outputs_{ts}.sdf").resolve())
                shutil.copy(f"dock_outputs_{ts}.sdf", Path("/home/a/"))
                sdfFile = os.path.join(RDConfig.RDDataDir, '/home/a/dock_outputs_{ts}.sdf.sdf')

                frame = PandasTools.LoadSDF(sdfFile, smilesName='SMILES', molColName='Molecule',
                                            includeFingerprints=True)

                print(frame.info)
                self.logger.info(f"dock_outputs_{ts}.sdf moved to {out_loc}")
        except Exception as e:
            self.logger.error(e)
        self.out.send(self.payload)

flow = Workflow(name="dock", level="debug", cleanup_temp=False)
flow.config.update(Path("configs/Maize/maize-mol2mol-config_2.toml"))
flow.config.scratch = (Path("./test_scratch"))

# retu = flow.add(Return[list[IsomerCollection]])


rnve = flow.add(ReInvent )

# @expose_reinvent
class DockGPU(Graph):
    """Dock multiple molecules in the form of SMILES using AutoDockGPU and return their scores."""

    inp: Input[list[str]]
    out: Output[NDArray[np.float32]]

    # Gypsum
    n_variants: Parameter[int]
    ph_range: Parameter[tuple[float, float]]
    n_jobs: Parameter[int]
    thoroughness: Parameter[int]

    # ADGPU
    grid_file: FileParameter[Path]
    heurmax: Parameter[int]

    # Mol save location
    base_path: FileParameter[Path]

    def build(self) -> None:
        gyp = self.add(Gypsum)
        adg = self.add(AutoDockGPU)
        log = self.add(SaveLibrary)
        gyp.use_filters.set(False)
        adg.scores_only.set(False)
        adg.strict.set(False)
        gyp.ph_range.set((6.5, 7.5))

        self.connect_all((gyp.out, adg.inp), (adg.out, log.inp))
        self.map(
            gyp.n_variants,
            gyp.ph_range,
            gyp.n_jobs,
            gyp.thoroughness,
            adg.grid_file,
            adg.heurmax,
            log.base_path,
        )
        self.inp = self.map_port(gyp.inp)
        self.out = self.map_port(adg.out_scores)


dock = flow.add(DockGPU,name="DockGPU", loop=True)
flow.connect_all(
    (rnve.out, dock.inp),
    (dock.out, rnve.inp),
)

# grid = Path("mols/Rad51/rad51_receptor.maps.fld")
# ref = Path("mols/Rad51/Cam833HMdsRad51Docked_Cam833-Acid_3.sdf")


rnv_config = Path("configs/REINVENT/staged_learning_maize_reinvent.toml")
prior = Path("priors/reinvent.prior")


rnve.configuration.set(rnv_config)
rnve.prior.set(prior)
rnve.agent.set(prior)

# rnve.distance_threshold.set(100)
# rnve.sample_strategy.set("multinomial")  # multinomial or beamsearch (deterministic)
# rnve.input_smi.set("mols/Rad51/cam833.smi")  # multinomial or beamsearch (deterministic)
rnve.weight.set(1)

rnve.min_epoch.set(100)
rnve.max_epoch.set(101)
rnve.batch_size.set(4)

# Settings to transform the docking score to a value between 0 and 1, with 1 being favourable, using a sigmoid
rnve.low.set(-12.0)
rnve.high.set(-3.0)
rnve.k.set(0.2)
rnve.reverse.set(True)

assert rnv_config.exists(), f"Config file not found: {rnv_config}"
assert prior.exists(), f"Prior model not found: {prior}"

# dock.n_jobs.set(4)  # parallel jobs
# dock.n_poses.set(10)
dock.n_variants.set(4)
dock.grid_file.set(Path("/home/a/REINVENT4/ReinventStudies/mols/Rad51/maps/rad51_receptor.maps.fld"))
# dock.heurmax.set(10000)
# dock.thoroughness.set(1)
dock.ph_range.set((6.5, 7.5))
filename= "qed_tanimoto_docking_reinvent"
ts = time.strftime("%B_%d-%H%M%S")
dock.base_path.set(Path(f"/home/a/REINVENT4/ReinventStudies/results/{filename}_{ts}_run"))


flow.check()
setup_workflow(flow)


#
#
# #
#
# # dock.search_center.set((139.95,76.64,221.97))  # (center_x, center_y, center_z) in Å
# # dock.search_range.set((20.0, 20.0, 20.0))  # (size_x, size_y, size_z) in Å

#
# isoSave.isomer_output_location.set(Path("results"))
#
# flow.visualize()
#
# flow.check()
# try:
#     flow.execute()
# except Exception as e:
#     print("e")
#     print(e)
#
