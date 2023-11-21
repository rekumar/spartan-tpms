from .marching_cubes import marching_cubes
import sdf
from typing import Optional, Tuple


def generate_stl(
    f: sdf.d3.SDF3,
    fpath: str,
    step: Optional[float] = 0.2,
    bounds: Optional[Tuple[Tuple[float]]] = None,
    verbose: Optional[bool] = True,
    clean: Optional[bool] = True,
):
    """Export an STL file of a given signed distance function

    Args:
        f (sdf.d3.SDF3): the structure to be exported as an STL
        fpath (str): filepath or filename to export to. Should end with ".stl"
        step (float, optional): The step size to render the STL to. Smaller values result in finer structures, but greatly increase the compute time. Defaults to 0.2.
        bounds (Tuple[Tuple[float]], optional): The bounds ((xmin, ymin, zmin), (xmax, ymax, zmax)) over which to generate the stl. If None, code will attempt to determine these bounds automatically. Defaults to None.
        verbose (bool, optional): Whether to print export status to the console. Defaults to True.
        clean (bool, optional): Whether to attempt cleaning dangling facets in the STL file. Defaults to True.
    """
    points = marching_cubes(f, step=step, bounds=bounds, verbose=verbose, clean=clean)
    sdf.write_binary_stl(fpath, points)
    if verbose:
        print(f"Wrote file to {fpath}")
