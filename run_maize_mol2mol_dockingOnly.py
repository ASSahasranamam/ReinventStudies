import time
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


flow = Workflow(name="dock", level="debug", cleanup_temp=False)
flow.config.update(Path("configs/Maize/maize-mol2mol-config_2.toml"))
flow.config.scratch = (Path("./pleaseWork/justmol2mol_4_Compounds_ontop"))
# retu = flow.add(Return[list[IsomerCollection]])


rnve = flow.add(ReInvent)


class SaveLibraryEveryEpoch(Node):
    """Save a list of molecules to multiple SDF files."""

    tags = {"saving"}

    inp: Input[list[IsomerCollection]] = Input()
    """Molecule library input"""

    base_path: FileParameter[Path] = FileParameter(exist_required=False)
    """Base output file path name without a suffix, i.e. /path/to/output"""

    output_tags: Parameter[list[str]] = Parameter(optional=True)
    """Tags to write out"""
    epoch: Parameter[int] = Parameter(optional=True, default=0)
    def run(self) -> None:
        mols = self.inp.receive()
        base = self.base_path.value
        epoch = self.epoch.value
        for i, mol in enumerate(mols):
            print(epoch)
            file = base.with_name(f"{base.name}_{epoch}_{i}.sdf")
            if self.output_tags.is_set:
                # this if statement is better than using mol[0].tags
                # as a default since it would support ragged tags
                mol.to_sdf(file, tags=self.output_tags.value)
                self.epoch.set(epoch+1)
            else:
                mol.to_sdf(file)



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
        log = self.add(SaveLibraryEveryEpoch)
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

runName = "justmol2mol_4_Compounds_Staged2"
rnv_config = Path("configs/REINVENT/staged_learning_maize_mol2mol_2.toml")
prior = Path("priors/mol2mol_similarity.prior")
agent = Path("priors/justmol2mol_4_compounds_staged_learning_1.chkpt")


rnve.weight.set(0.5)
# Make sure these files exist and are readable
assert rnv_config.exists(), f"Config file not found: {rnv_config}"
assert prior.exists(), f"Prior model not found: {prior}"

# grid = Path("../maize/steps/mai/docking/data/1uyd.tar")


rnve.configuration.set(rnv_config)
rnve.prior.set(prior)
rnve.agent.set(agent)

#Mol2molOnly
inp_smi = "mols/Rad51/cam833.smi"
rnve.distance_threshold.set(100) 
rnve.sample_strategy.set("multinomial")  # multinomial or beamsearch (deterministic)
rnve.input_smi.set(inp_smi)


rnve.weight.set(0.5)

rnve.min_epoch.set(100)
rnve.max_epoch.set(102)
rnve.batch_size.set(32)

# Settings to transform the docking score to a value between 0 and 1, with 1 being favourable, using a sigmoid
rnve.low.set(-15.0)
rnve.high.set(-3.0)
rnve.k.set(0.2)
rnve.reverse.set(True)
rnve.maize_backend.set(True)
assert rnv_config.exists(), f"Config file not found: {rnv_config}"
assert prior.exists(), f"Prior model not found: {prior}"

# dock.n_jobs.set(4)  # parallel jobs
# dock.n_poses.set(10)

dock.n_variants.set(4)
dock.grid_file.set(Path("/home/a/REINVENT4/ReinventStudies/mols/Rad51/maps/rad51_receptor.maps.fld"))
# dock.heurmax.set(10000)
# dock.thoroughness.set(1)
dock.ph_range.set((6.5, 7.5))
ts = time.strftime("%B_%d-%H%M%S")
dock.base_path.set(Path(f"/home/a/REINVENT4/ReinventStudies/pleaseWork/{runName}_{ts}_run"))
# for smiles in ["mols/Rad51/cam833.smi", "mols/Rad51/cam833_compounds.smi"]:
#     # for searchStrategy in ["beamsearch","multinomial"]:
#


flow.check()
flow.visualize()
setup_workflow(flow)





