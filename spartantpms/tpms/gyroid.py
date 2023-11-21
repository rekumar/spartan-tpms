import numpy as np
from sdf import sdf3

GYROID_T_VALUE = [
    -1.55563492,
    -1.49340952,
    -1.43118413,
    -1.36895873,
    -1.30673333,
    -1.24450793,
    -1.18228254,
    -1.12005714,
    -1.05783174,
    -0.99560635,
    -0.93338095,
    -0.87115555,
    -0.80893016,
    -0.74670476,
    -0.68447936,
    -0.62225397,
    -0.56002857,
    -0.49780317,
    -0.43557778,
    -0.37335238,
    -0.31112698,
    -0.24890159,
    -0.18667619,
    -0.12445079,
    -0.0622254,
    0.0,
    0.0622254,
    0.12445079,
    0.18667619,
    0.24890159,
    0.31112698,
    0.37335238,
    0.43557778,
    0.49780317,
    0.56002857,
    0.62225397,
    0.68447936,
    0.74670476,
    0.80893016,
    0.87115555,
    0.93338095,
    0.99560635,
    1.05783174,
    1.12005714,
    1.18228254,
    1.24450793,
    1.30673333,
    1.36895873,
    1.43118413,
    1.49340952,
    1.55563492,
][::-1]

POROSITY = [
    1.00000000e00,
    9.99751529e-01,
    9.89253626e-01,
    9.66176875e-01,
    9.41842238e-01,
    9.17787132e-01,
    8.95191793e-01,
    8.72394572e-01,
    8.50637823e-01,
    8.27875543e-01,
    8.05070557e-01,
    7.85709225e-01,
    7.63568122e-01,
    7.41427020e-01,
    7.21040744e-01,
    7.01054352e-01,
    6.80730194e-01,
    6.59889683e-01,
    6.41196116e-01,
    6.20091604e-01,
    5.99938270e-01,
    5.79392818e-01,
    5.59472426e-01,
    5.40452741e-01,
    5.20699291e-01,
    5.00032515e-01,
    4.79300709e-01,
    4.59547259e-01,
    4.40527574e-01,
    4.20607182e-01,
    4.00061730e-01,
    3.79908396e-01,
    3.58803884e-01,
    3.40110317e-01,
    3.19269806e-01,
    2.98945648e-01,
    2.78959256e-01,
    2.58572980e-01,
    2.36431878e-01,
    2.14290775e-01,
    1.94929443e-01,
    1.72124457e-01,
    1.49362177e-01,
    1.27605428e-01,
    1.04808207e-01,
    8.22128679e-02,
    5.81577617e-02,
    3.38231255e-02,
    1.07463741e-02,
    2.48471078e-04,
    0.00000000e00,
][::-1]


def _gyroid_porosity_to_t_value(porosity: float) -> float:
    return np.interp(porosity, POROSITY, GYROID_T_VALUE)


@sdf3
def gyroid_function(
    lambda_x: float, lambda_y: float, lambda_z: float, porosity: float = 0.5
):
    """Implicit surface function for the gyroid surface.

    Args:
        lambda_x (float): Gyroid wavelength in the x dimension.
        lambda_y (float): Gyroid wavelength in the y dimension.
        lambda_z (float): Gyroid wavelength in the z dimension.
        porosity (float, optional): Pore fraction of the structure. This controls the wall thickness of the diamond lattice. Larger values = more porous = thinner walls.
    """

    t = _gyroid_porosity_to_t_value(porosity)

    def f(p):
        x, y, z = (
            (2 * np.pi / lambda_per_dim) * p[:, i]
            for i, lambda_per_dim in enumerate([lambda_x, lambda_y, lambda_z])
        )

        return np.cos(x) * np.sin(y) + np.cos(y) * np.sin(z) + np.cos(z) * np.sin(x) - t

    return f
