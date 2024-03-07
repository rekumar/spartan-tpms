from tqdm.notebook import tqdm
import pandas as pd
import os
from spartantpms.tpms.gyroid import sheet_gyroid_function
from spartantpms import generate_stl

from sdf import box
import pandas as pd
import os

FOLDER = "initial 25 sheet structures_medium"
samples = pd.read_csv("20240213 Sheet Samples with Lambda Driven.csv", index_col=0)

try:
    os.mkdir(FOLDER)
except:
    print("Folder already exists")

from pymeshlab import MeshSet
import numpy as np


def simplify_stl(fpath, outpath, ratio=0.5, num_steps=5):
    # create a new MeshSet
    ms = MeshSet()
    ms.load_new_mesh(fpath)
    ms.meshing_remove_connected_component_by_face_number(
        mincomponentsize=10000
    )  # value came from JiYoung single unit cells, not sure if 10k scales to larger models. I believe this removes the dangling bits of the repeat unit at the corners.

    current_facenum = ms.current_mesh().face_number()
    target_facenum = int(current_facenum * ratio)
    targets = np.linspace(current_facenum, target_facenum, num_steps)

    for target in targets:
        ms.meshing_decimation_quadric_edge_collapse(
            targetfacenum=int(target),
            preserveboundary=True,
            preservenormal=True,
            preservetopology=True,
            optimalplacement=True,
            planarquadric=True,
            autoclean=True,
        )
    ms.save_current_mesh(outpath)


for structure_number, row in samples.iterrows():
    f = sheet_gyroid_function(
        lambda_x=row["lambda"],
        lambda_y=row["lambda"],
        lambda_z=row["lambda"],
        theta_x=0,
        theta_y=0,
        theta_z=0,  # rotation about z is degenerate when stressing along z axis.
        porosity=row["porosity"],
    )

    f = f & box((45, 45, 135))  # total size is defined here

    fname = f"{FOLDER}/{structure_number}_sheet_45x45x90mm_{row['porosity']}porosity_{row['thickness']}sheetthickness"
    generate_stl(
        f=f,
        fpath=f"{fname}.stl",
        step=row["thickness"] / 6,
    )

    simplify_stl(
        fpath=f"{fname}.stl",
        outpath=f"{fname}_simplified.stl",
        ratio=0.05,  # reduce to 5% of original mesh count
        num_steps=5,
    )
