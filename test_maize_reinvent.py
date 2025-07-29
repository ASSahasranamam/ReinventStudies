from pathlib import Path
from maize.core.workflow import Workflow
from maize.steps.io import LoadData, LogResult, Return, Void
from maize.steps.mai.docking.adv import AutoDockGPU
from maize.steps.mai.molecule import Gypsum
from maize.utilities.chem import IsomerCollection
from maize.steps.mai.misc import ReInvent
from maize.steps.mai.cheminformatics import RMSD, ExtractScores, TagIndex, SortByTag, TagSorter, LogTags


grid_FL = Path("/home/a/maize-contrib/maize/steps/mai/docking/data/1stp/1stp_protein.maps.fld")

flow = Workflow(name="dock", level="info", cleanup_temp=False)
flow.config.update(Path("test-config.toml"))

# load = flow.add(LoadData[list[str]])
rnve = flow.add(ReInvent)
embe = flow.add(Gypsum)
dock = flow.add(AutoDockGPU)
# void = flow.add(Void)
retu = flow.add(Return[list[IsomerCollection]])
# scor = flow.add(ExtractScores, loop=True)
# scor = flow.add(ExtractScores, loop=True)

# load.data.set(["Nc1ccc(ccc1N)C", "Nc1ccc(cc1N)C"])
embe.n_variants.set(2)



flow.connect(rnve.out, embe.inp)
# flow.connect(load.out, embe.inp)
flow.connect(embe.out, dock.inp)
flow.connect(dock.out, retu.inp)
flow.connect(dock.out_scores, rnve.inp)


grid = Path("/home/a/maize-contrib/maize/steps/mai/docking/data/1uyd.tar")
ref = Path("~/maize-contrib/maize/steps/mai/docking/data/1UYD_ligand.sdf")
rnv_config = Path("configs/staged_learning_maize.toml")
prior = Path("priors/reinvent.prior")


rnve.configuration.set(rnv_config)
rnve.prior.set(prior)
rnve.agent.set(prior)

# The maximum number of RL epochs
rnve.max_epoch.set(10)

# Settings to transform the docking score to a value between 0 and 1, with 1 being favourable, using a sigmoid
rnve.low.set(-10.0)
rnve.high.set(-5.0)
rnve.reverse.set(True)

# Number of molecules to generate each epoch
rnve.batch_size.set(32)

# Number of isomers to generate for each SMILES
embe.n_variants.set(4)

# Docking grid for
# dock.grid_file.set(grid_FL)
dock.inp_grid.set(grid)
dock.constraints.set(False)
# dock_hp.inp_grid.set(grid)

# Reference ligand for RMSD calculation
# load.path.set(ref)

flow.check()
flow.execute()