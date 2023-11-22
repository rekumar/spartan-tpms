import numpy as np
from sdf import sdf3
import sdf
from .euler import rotate_about_xyz


def _gyroid_porosity_to_t_value(porosity: float) -> float:
    """Generate the `t` value to yield a gyroid with the desired porosity. Porosity is defined as the fraction of the volume that is empty space.

    Args:
        porosity (float): Desired porosity of the gyroid. Must be between 0 and 1.

    Returns:
        float: t-value to subtract from the signed distance function of the gyroid surface to yield a structure with the desired porosity.
    """
    if porosity < 0 or porosity > 1:
        raise ValueError("Porosity must be between 0 and 1.")

    polynomial_constants = [
        -9.86328589e02,
        4.43858792e03,
        -8.42057860e03,
        8.75829993e03,
        -5.43754382e03,
        2.05512365e03,
        -4.60683827e02,
        5.58905464e01,
        -5.79899099e00,
        1.51588868e00,
    ]
    return np.polyval(polynomial_constants, porosity)


@sdf3
def gyroid_function(
    lambda_x: float,
    lambda_y: float,
    lambda_z: float,
    porosity: float = 0.5,
    theta_x: float = 0,
    theta_y: float = 0,
    theta_z: float = 0,
) -> sdf.d3.SDF3:
    """Implicit surface function for the gyroid surface.

    Args:
        lambda_x (float): Gyroid wavelength in the x dimension.
        lambda_y (float): Gyroid wavelength in the y dimension.
        lambda_z (float): Gyroid wavelength in the z dimension.
        theta_x (float): Rotation about the x-axis in degrees.
        theta_y (float): Rotation about the y-axis in degrees.
        theta_z (float): Rotation about the z-axis in degrees.
        porosity (float, optional): Pore fraction of the structure. This controls the wall thickness of the gyroid lattice. Larger values = more porous = thinner walls.

    Returns:
        sdf.d3.SDF3: signed distance function representation of the structure.
    """

    t = _gyroid_porosity_to_t_value(porosity)

    def f(p):
        p = p * np.array(
            [
                2 * np.pi / lambda_per_dim
                for lambda_per_dim in [lambda_x, lambda_y, lambda_z]
            ]
        )  # put wavelength prefactors into each coordinate
        p = p @ rotate_about_xyz(
            theta_x, theta_y, theta_z, degrees=True
        )  # rotate the coordinate system
        x, y, z = p.T

        return np.cos(x) * np.sin(y) + np.cos(y) * np.sin(z) + np.cos(z) * np.sin(x) - t

    return f
