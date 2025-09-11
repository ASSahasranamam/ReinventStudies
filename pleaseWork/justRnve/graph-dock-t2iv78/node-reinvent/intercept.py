#!/usr/bin/env python
from pathlib import Path
import shutil
from time import sleep
import sys

with Path("inp.smi").open("w") as file:
    file.writelines(sys.stdin.readlines())
while not Path("out.json").exists():
    sleep(0.5)
with Path("out.json").open("r") as file:
    print(file.read())
Path("out.json").unlink()
