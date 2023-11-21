import numpy as np
from sdf.mesh import _estimate_bounds, _cartesian_product
import time
import PyScaffolder


# method to convert SDF to STL mesh.
def marching_cubes(f, step=0.01, bounds=None, verbose=True, clean=True):
    if not bounds:
        bounds = _estimate_bounds(f)

    (x0, y0, z0), (x1, y1, z1) = bounds

    try:
        dx, dy, dz = step
    except TypeError:
        dx = dy = dz = step

    if verbose:
        print("min %g, %g, %g" % (x0, y0, z0))
        print("max %g, %g, %g" % (x1, y1, z1))
        print("step %g, %g, %g" % (dx, dy, dz))

    X = np.arange(x0, x1, dx)
    Y = np.arange(y0, y1, dy)
    Z = np.arange(z0, z1, dz)
    P = _cartesian_product(X, Y, Z)

    try:
        # Since the PyScaffolder marching_cubes aceept FREP: F(x,y,z) > 0
        # Then the negative of implicit function is used
        Fxyz = -f(P)
        # Reshape to Fortran array (column-based) due to implementation of dualmc starting from z axis to x
        Fxyz = Fxyz.reshape((len(X), len(Y), len(Z))).reshape(-1, order="F")
        start = time.time()
        (v, f) = PyScaffolder.marching_cubes(
            Fxyz,
            grid_size=[len(X), len(Y), len(Z)],
            v_min=bounds[0],
            delta=step,
            clean=clean,
        )
        if verbose:
            seconds = time.time() - start
            # print('\n%d triangles in %g seconds' % (len(points) // 3, seconds))
        # merge vertices and faces into points
        return v[f].reshape((-1, 3))

    except Exception as e:
        print(e)
        return np.array([])
