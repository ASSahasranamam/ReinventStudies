import subprocess

import time
from maize.utilities.execution import CommandRunner
from pathlib import Path
from rdkit import Chem
from rdkit.Chem.ChemUtils.SDFToCSV import Convert
from rdkit.Chem import PandasTools
from rdkit import RDConfig
from exceptiongroup import catch
from maize.core.workflow import Workflow
from maize.steps.io import LoadData, LogResult, Return, Void
from maize.steps.mai.docking.adv import AutoDockGPU
from maize.steps.mai.molecule import Gypsum, LoadMolecule, SaveCSV, IsomerCollectionSaving
from maize.utilities.chem import IsomerCollection, save_sdf_library
from maize.utilities.io import setup_workflow
from maize.steps.mai.misc import ReInvent
from maize.steps.mai.cheminformatics import RMSD, ExtractScores, TagIndex, SortByTag, TagSorter, LogTags
from maize.steps.plumbing import MergeLists
import os
import stat

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


class Boltz(Node):
    inp: Input[list[IsomerCollection]] = Input()
    # out: Output[list[str]] = Output()

    prot_seq: Parameter[str] = Parameter(optional=True)
    prot_ID: Parameter[str] = Parameter(optional=True)
    lig_ID: Parameter[str] = Parameter(optional=True)
    tags = {"chemistry", "ml", "generation"}

    required_callables = ["boltz"]

    scores = []

    def run(self) -> None:
        smiles: list[IsomerCollection] = self.inp.receive()

        print(smiles)
        for j in range(len(smiles)):
            smiles_code = smiles[j].to_smiles()[0]
            print(smiles_code)
            yaml_path = Path(f"./{j}.yaml")
            yaml_path.unlink(missing_ok=True)

            templateYaml = f"""
version: 1  # Optional, defaults to 1
sequences:
  - protein:
      id: {self.prot_ID.value}
      sequence: {self.prot_seq.value}
      msa: '/home/a/REINVENT4/ReinventStudies/mols/Rad51/rad51_receptor.a3m'
  - ligand:
      id: {self.lig_ID.value}
      smiles: '{smiles_code}'
properties:
  - affinity:
      binder: {self.lig_ID.value}
"""
            yaml_path.write_text(templateYaml, encoding="utf-8")
            os.chmod(yaml_path, mode=stat.S_IXUSR | stat.S_IRUSR | stat.S_IWUSR)

            command = (
                f"{self.runnable['boltz']} "
                f"predict {yaml_path.as_posix()} "
            )
            x = subprocess.run([ f"{self.runnable['boltz']} boltz --help "], capture_output=True)
            print(x.stdout)
            print(x.stderr)
            # cmd = CommandRunner(working_dir=self.work_dir, rm_config=self.config.batch_config)
            # worker = cmd.run_async(command)
            count = 0
            # while worker.is_alive() and not self.signal.is_set():
            #     self.logger.debug(f"Checking if ReInvent is still running, counter: {count}")
            #     target_output_path = self.work_dir / "summary_1.csv"
            #
            #
            #     time.sleep(60)
            #
            # self.logger.info("Loop complete, stopping worker")
            # worker.kill(timeout=5)
            # worker.wait()


class Readonly(Node):
    smiles: Parameter[list[str]] = Parameter()

    out: Output[list[IsomerCollection]] = Output()

    def run(self) -> None:
        isoColList: list[IsomerCollection] =[]

        for i in self.smiles.value:
            IsoCol = IsomerCollection.from_smiles(i)
            isoColList.append(IsoCol)
        self.out.send(isoColList)


flow = Workflow(name="dock", level="debug", cleanup_temp=False)
flow.config.update(Path("configs/Maize/maize-mol2mol-config_2.toml"))

# retu = flow.add(Return[list[IsomerCollection]])
# rnve = flow.add(ReInvent)
#

readonly = flow.add(Readonly)
blitz = flow.add(Boltz)
flow.connect(readonly.out, blitz.inp)

readonly.smiles.set(["C[C@@H](c1c(cc(cc1)OC)Cl)NC(=O)[C@@H]1C[C@H](CN1C(=O)CNC(=O)c1nc2c(cc1)cc(cc2)F)O"])
blitz.prot_seq.set("MAMQMQLEANADTSVEEESFGPQPISRLEQCGINANDVKKLEEAGFHTVEAVAYAPKKELINIKGISEAKADKILAEAAKLVPMGFTTATEFHQRRSEIIQITTGSKELDKLLQGGIETGSITEMFGEFRTGKTQICHTLAVTCQLPIDRGGGEGKAMYIDTEGTFRPERLLAVAERYGLSGSDVLDNVAYARAFNTDHQTQLLYQASAMMVESRYALLIVDSATALYRTDYSGRGELSARQMHLARFLRMLLRLADEFGVAVVITNQVVAQVDGAAMFAADPKKPIGGNIIAHASTTRLYLRKGRGETRICKIYDSPCLPEAEAMFAINADGVGDAKD")
blitz.prot_ID.set("A")
blitz.lig_ID.set("B")

flow.check()
try:
    flow.execute()
except Exception as e:
    print("e")
    print(e)

# flow.to_file("run_maize_mol2mol.yml")

ts = time.strftime("%Y%m%d-%H%M%S")

# save_sdf_library(mols=totalIsomerCollection,file=Path( f"dock_outputs_{ts}.sdf"))
# except Exception as e:
# flow.logger.error(e)
# flow._cleanup()
