import numpy as np
from sdf import sdf3

POROSITY = [
    1.0,
    0.99726682,
    0.99099292,
    0.98260702,
    0.97316512,
    0.96080369,
    0.94682719,
    0.92794339,
    0.90347675,
    0.879938,
    0.85341759,
    0.82876072,
    0.80354479,
    0.77733497,
    0.75664976,
    0.72975664,
    0.70963049,
    0.68292373,
    0.66391569,
    0.63714681,
    0.61627524,
    0.59105931,
    0.57000139,
    0.54348098,
    0.52335483,
    0.5022804,
    0.47664517,
    0.45651902,
    0.42999861,
    0.40894069,
    0.38372476,
    0.36285319,
    0.33608431,
    0.31707627,
    0.29036951,
    0.27024336,
    0.24335024,
    0.22266503,
    0.19645521,
    0.17123928,
    0.14658241,
    0.120062,
    0.09652325,
    0.07205661,
    0.05317281,
    0.03919631,
    0.02683488,
    0.01739298,
    0.00900708,
    0.00273318,
    0.0,
][::-1]

DIAMOND_T_VALUE = [
    -1.41421356,
    -1.35764502,
    -1.30107648,
    -1.24450793,
    -1.18793939,
    -1.13137085,
    -1.07480231,
    -1.01823376,
    -0.96166522,
    -0.90509668,
    -0.84852814,
    -0.79195959,
    -0.73539105,
    -0.67882251,
    -0.62225397,
    -0.56568542,
    -0.50911688,
    -0.45254834,
    -0.3959798,
    -0.33941125,
    -0.28284271,
    -0.22627417,
    -0.16970563,
    -0.11313708,
    -0.05656854,
    0.0,
    0.05656854,
    0.11313708,
    0.16970563,
    0.22627417,
    0.28284271,
    0.33941125,
    0.3959798,
    0.45254834,
    0.50911688,
    0.56568542,
    0.62225397,
    0.67882251,
    0.73539105,
    0.79195959,
    0.84852814,
    0.90509668,
    0.96166522,
    1.01823376,
    1.07480231,
    1.13137085,
    1.18793939,
    1.24450793,
    1.30107648,
    1.35764502,
    1.41421356,
][::-1]


def _diamond_porosity_to_t_value(porosity: float) -> float:
    return np.interp(porosity, POROSITY, DIAMOND_T_VALUE)


@sdf3
def diamond_function(
    lambda_x: float, lambda_y: float, lambda_z: float, porosity: float = 0.5
):
    """Implicit surface function for the diamond surface.

    Args:
        lambda_x (float): Diamond wavelength in the x dimension.
        lambda_y (float): Diamond wavelength in the y dimension.
        lambda_z (float): Diamond wavelength in the z dimension.
        porosity (float, optional): Pore fraction of the structure. This controls the wall thickness of the diamond lattice. Larger values = more porous = thinner walls.
    """

    t = _diamond_porosity_to_t_value(porosity)

    def f(p):
        x, y, z = (
            (2 * np.pi / lambda_per_dim) * p[:, i]
            for i, lambda_per_dim in enumerate([lambda_x, lambda_y, lambda_z])
        )

        return (
            np.sin(x) * np.sin(y) * np.sin(z)
            + np.sin(x) * np.cos(y) * np.cos(z)
            + np.cos(x) * np.sin(y) * np.cos(z)
            + np.cos(x) * np.cos(y) * np.sin(z)
            - t
        )

    return f
