import numpy as np
from sdf import sdf3
from .euler import rotate_about_xyz


def _diamond_porosity_to_t_value(porosity: float) -> float:
    """Generate the `t` value to yield a diamond lattice with the desired porosity. Porosity is defined as the fraction of the volume that is empty space.

    Args:
        porosity (float): Desired porosity of the diamond lattice. Must be between 0 and 1.

    Returns:
        float: t-value to subtract from the signed distance function of the diamond surface to yield a structure with the desired porosity.
    """
    if porosity < 0 or porosity > 1:
        raise ValueError("Porosity must be between 0 and 1.")

    polynomial_constants = [
        -1.40181203e03,
        6.31442849e03,
        -1.20734024e04,
        1.27714436e04,
        -8.16588099e03,
        3.23757850e03,
        -7.84778926e02,
        1.09747587e02,
        -1.01136257e01,
        1.39501223e00,
    ]
    return np.polyval(polynomial_constants, porosity)


@sdf3
def diamond_function(
    lambda_x: float,
    lambda_y: float,
    lambda_z: float,
    theta_x: float,
    theta_y: float,
    theta_z: float,
    porosity: float = 0.5,
):
    """Implicit surface function for the diamond surface.

    Args:
        lambda_x (float): Diamond wavelength in the x dimension.
        lambda_y (float): Diamond wavelength in the y dimension.
        lambda_z (float): Diamond wavelength in the z dimension.
        theta_x (float): Rotation about the x-axis in degrees.
        theta_y (float): Rotation about the y-axis in degrees.
        theta_z (float): Rotation about the z-axis in degrees.
        porosity (float, optional): Pore fraction of the structure. This controls the wall thickness of the diamond lattice. Larger values = more porous = thinner walls.
    """

    t = _diamond_porosity_to_t_value(porosity)

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

        return (
            np.sin(x) * np.sin(y) * np.sin(z)
            + np.sin(x) * np.cos(y) * np.cos(z)
            + np.cos(x) * np.sin(y) * np.cos(z)
            + np.cos(x) * np.cos(y) * np.sin(z)
            - t
        )

    return f
