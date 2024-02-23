import numpy as np
from sdf import sdf3
import sdf
from .euler import rotate_about_xyz
from typing import Optional
from warnings import warn


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
    theta_x: float = 0,
    theta_y: float = 0,
    theta_z: float = 0,
    porosity: float = 0.5,
    isovalue: Optional[float] = None,
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
        isovalue (float, optional): Isovalue to use for the implicit surface function. If None, the isovalue is set to the t-value that yields the desired porosity. Default is None.

    Returns:
        sdf.d3.SDF3: signed distance function representation of the structure.
    """

    if isovalue is None:
        t = _gyroid_porosity_to_t_value(porosity)
    else:
        if porosity != 0.5:
            raise ValueError(
                "Porosity can only be set if isovalue is set to None. Do not provide both values!"
            )
        t = isovalue

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


def _sheet_thickness_to_t_value(thickness: float, lambda_xyz: float) -> float:
    """Generate the `t` value to yield a sheet gyroid with the desired thickness. Thickness is defined as the fraction of the volume that is empty space.

    Args:
        thickness (float): Desired thickness of the sheet gyroid (mm). Must be between 0 and 0.88 times lambda_xyz.
        lambda_xyz (float): Wavelength (mm) of the gyroid in the x, y, and z dimensions.

    Returns:
        float: t-value to input to `sheet_gyroid_function` to yield a structure with the desired thickness.
    """

    t_norm = thickness / lambda_xyz
    if t_norm <= 0 or t_norm >= 0.88:
        raise ValueError("Thickness must be between 0 and 0.88 times lambda_xyz.")

    thickness_normalized = [
        0.004,
        0.012,
        0.02,
        0.028,
        0.04,
        0.048,
        0.056,
        0.064,
        0.076,
        0.084,
        0.092,
        0.104,
        0.112,
        0.12,
        0.132,
        0.14,
        0.152,
        0.16,
        0.172,
        0.184,
        0.192,
        0.204,
        0.216,
        0.228,
        0.24,
        0.252,
        0.264,
        0.276,
        0.292,
        0.308,
        0.324,
        0.34,
        0.36,
        0.384,
        0.412,
        0.452,
        0.52878729,
        0.55713553,
        0.58060658,
        0.60113226,
        0.62096699,
        0.64089937,
        0.66020603,
        0.67976467,
        0.69945121,
        0.71929966,
        0.74095074,
        0.76459663,
        0.79095891,
        0.82341727,
        0.8830855,
    ]

    isovalues = [
        1.000000e-04,
        5.659800e-02,
        1.130960e-01,
        1.695940e-01,
        2.260920e-01,
        2.825900e-01,
        3.390880e-01,
        3.955860e-01,
        4.520840e-01,
        5.085820e-01,
        5.650800e-01,
        6.215780e-01,
        6.780760e-01,
        7.345740e-01,
        7.910720e-01,
        8.475700e-01,
        9.040680e-01,
        9.605660e-01,
        1.017064e00,
        1.073562e00,
        1.130060e00,
        1.186558e00,
        1.243056e00,
        1.299554e00,
        1.356052e00,
        1.412550e00,
        1.469048e00,
        1.525546e00,
        1.582044e00,
        1.638542e00,
        1.695040e00,
        1.751538e00,
        1.808036e00,
        1.864534e00,
        1.921032e00,
        1.977530e00,
        2.034028e00,
        2.090526e00,
        2.147024e00,
        2.203522e00,
        2.260020e00,
        2.316518e00,
        2.373016e00,
        2.429514e00,
        2.486012e00,
        2.542510e00,
        2.599008e00,
        2.655506e00,
        2.712004e00,
        2.768502e00,
        2.825000e00,
    ]

    return np.interp(t_norm, thickness_normalized, isovalues)


@sdf3
def sheet_gyroid_function(
    lambda_x: float,
    lambda_y: float,
    lambda_z: float,
    theta_x: float = 0,
    theta_y: float = 0,
    theta_z: float = 0,
    porosity: Optional[float] = 0.5,
    isovalue: Optional[float] = None,
    sheet_thickness: Optional[float] = None,
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
        isovalue (float, optional): Isovalue to use for the implicit surface function. If None, the isovalue is set to the t-value that yields the desired porosity. Default is None.
        sheet_thickness (float, optional): Thickness of the sheet gyroid. Must be between 0 and 0.88 times the wavelength in the x, y, or z dimension. Default is None.

    Returns:
        sdf.d3.SDF3: signed distance function representation of the structure.
    """

    # porosity argument
    if isovalue is None and sheet_thickness is None:
        d_porosity = (1 - porosity) / 2

        f_inner = gyroid_function(
            lambda_x,
            lambda_y,
            lambda_z,
            theta_x,
            theta_y,
            theta_z,
            0.5 + d_porosity,
        )
        f_outer = gyroid_function(
            lambda_x,
            lambda_y,
            lambda_z,
            theta_x,
            theta_y,
            theta_z,
            0.5 - d_porosity,
        )
        return f_outer - f_inner

    # sheet_thickness argument
    if sheet_thickness is not None:
        if isovalue is not None:
            raise ValueError(
                "Sheet thickness can only be set if isovalue is set to None. Do not provide both values!"
            )

        if lambda_x != lambda_y or lambda_x != lambda_z or lambda_y != lambda_z:
            warn(
                "The `sheet_thickness` argument uses a method that was derived for isotropic gyroids (lambda_x = lambda_y = lambda_z). You provided an anisotropic structure, so the sheet thickness may not be accurate. Using the average wavelength for all dimensions to calculate sheet thickness."
            )
        isovalue = _sheet_thickness_to_t_value(
            sheet_thickness, np.mean([lambda_x, lambda_y, lambda_z])
        )

    if isovalue <= 0:
        raise ValueError(
            "Isovalue for a sheet gyroid must be positive. This value is split in half and mirrored about 0 to yield the sheet gyroid. For example, if isovalue=0.5, the sheet gyroid will be generated as the difference of gyroids with isovalue=0.25 and isovalue=-0.25."
        )
    if porosity != 0.5:
        raise ValueError(
            "Porosity can only be set if isovalue is set to None. Do not provide both values!"
        )
    d_t = isovalue / 2

    f_inner = gyroid_function(
        lambda_x, lambda_y, lambda_z, theta_x, theta_y, theta_z, isovalue=-d_t
    )
    f_outer = gyroid_function(
        lambda_x, lambda_y, lambda_z, theta_x, theta_y, theta_z, isovalue=d_t
    )

    return f_outer - f_inner
