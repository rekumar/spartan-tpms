import numpy as np
import sdf
from spartantpms.tpms import gyroid_function, diamond_function


def gyroid_box(
    lambda_x: float,
    lambda_y: float,
    lambda_z: float,
    theta_x: float,
    theta_y: float,
    theta_z: float,
    porosity: float,
    n_periods: int = 4,
) -> sdf.d3.SDF3:
    """Generates an implicit representation of the test structure, which is a gyroid lattice inside of a rectangular volume.

    Args:
        lambda_x (float): Gyroid wavelength in the x dimension.
        lambda_y (float): Gyroid wavelength in the y dimension.
        lambda_z (float): Gyroid wavelength in the z dimension.
        theta_x (float): Rotation about the x-axis in degrees.
        theta_y (float): Rotation about the y-axis in degrees.
        theta_z (float): Rotation about the z-axis in degrees.
        porosity (float, optional): Pore fraction of the structure. This controls the wall thickness of the gyroid lattice. Larger values = more porous = thinner walls.
        n_periods (int, optional): Number of periods of the gyroid to include in the structure. This defines the total size of the structure as multiples of lambda_x, lambda_y, and lambda_z. Defaults to 4.

    Returns:
        sdf.d3.SDF3: signed distance function representation of the structure.
    """
    dx = lambda_x * n_periods
    dy = lambda_y * n_periods
    dz = lambda_z * n_periods

    infill = gyroid_function(
        lambda_x=lambda_x,
        lambda_y=lambda_y,
        lambda_z=lambda_z,
        theta_x=theta_x,
        theta_y=theta_y,
        theta_z=theta_z,
        porosity=porosity,
    )

    box = sdf.box(size=(dx, dy, dz))

    return box & infill


def diamond_box(
    lambda_x: float,
    lambda_y: float,
    lambda_z: float,
    theta_x: float,
    theta_y: float,
    theta_z: float,
    porosity: float,
    n_periods: int = 4,
) -> sdf.d3.SDF3:
    """Generates an implicit representation of the test structure, which is a diamond lattice inside of a rectangular volume.

    Args:
        lambda_x (float): Gyroid wavelength in the x dimension.
        lambda_y (float): Gyroid wavelength in the y dimension.
        lambda_z (float): Gyroid wavelength in the z dimension.
        theta_x (float): Rotation about the x-axis in degrees.
        theta_y (float): Rotation about the y-axis in degrees.
        theta_z (float): Rotation about the z-axis in degrees.
        porosity (float, optional): Pore fraction of the structure. This controls the wall thickness of the diamond lattice. Larger values = more porous = thinner walls.
        n_periods (int, optional): Number of periods of the diamond lattice to include in the structure. This defines the total size of the structure as multiples of lambda_x, lambda_y, and lambda_z. Defaults to 4.

    Returns:
        sdf.d3.SDF3: signed distance function representation of the structure.
    """
    dx = lambda_x * n_periods
    dy = lambda_y * n_periods
    dz = lambda_z * n_periods

    infill = diamond_function(
        lambda_x=lambda_x,
        lambda_y=lambda_y,
        lambda_z=lambda_z,
        theta_x=theta_x,
        theta_y=theta_y,
        theta_z=theta_z,
        porosity=porosity,
    )

    box = sdf.box(size=(dx, dy, dz))

    return box & infill
