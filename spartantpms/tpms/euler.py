"""
Use Euler angles and rotation matrices to rotate vectors in R3.

We use this to rotate coordinate systems that define the TPMS lattice -> rotate the structure.

https://en.wikipedia.org/wiki/Rotation_matrix#In_three_dimensions
"""
import numpy as np


def rotate_about_x(theta: float, degrees: bool = True) -> np.array:
    """Returns a rotation matrix about the x-axis.

    Args:
        theta (float): Angle of rotation in radians.
        degrees (bool, optional): Whether theta is in degrees (True) or radians (False). Defaults to True.

    Returns:
        np.array: 3x3 rotation matrix. Multiply an R3 vector by this matrix to rotate it about the x-axis.
    """
    if degrees:
        theta = np.deg2rad(theta)

    return np.array(
        [
            [1, 0, 0],
            [0, np.cos(theta), -np.sin(theta)],
            [0, np.sin(theta), np.cos(theta)],
        ]
    )


def rotate_about_y(theta: float, degrees: bool = True) -> np.array:
    """Returns a rotation matrix about the y-axis.

    Args:
        theta (float): Angle of rotation in radians.
        degrees (bool, optional): Whether theta is in degrees (True) or radians (False). Defaults to True.

    Returns:
        np.array: 3x3 rotation matrix. Multiply an R3 vector by this matrix to rotate it about the y-axis.
    """
    if degrees:
        theta = np.deg2rad(theta)

    return np.array(
        [
            [np.cos(theta), 0, np.sin(theta)],
            [0, 1, 0],
            [-np.sin(theta), 0, np.cos(theta)],
        ]
    )


def rotate_about_z(theta: float, degrees: bool = True) -> np.array:
    """Returns a rotation matrix about the z-axis.

    Args:
        theta (float): Angle of rotation in radians.
        degrees (bool, optional): Whether theta is in degrees (True) or radians (False). Defaults to True.

    Returns:
        np.array: 3x3 rotation matrix. Multiply an R3 vector by this matrix to rotate it about the z-axis.
    """
    if degrees:
        theta = np.deg2rad(theta)

    return np.array(
        [
            [np.cos(theta), -np.sin(theta), 0],
            [np.sin(theta), np.cos(theta), 0],
            [0, 0, 1],
        ]
    )


def rotate_about_xyz(
    theta_x: float, theta_y: float, theta_z: float, degrees: bool = True
) -> np.array:
    """Returns a rotation matrix about the x, y, and z axes.

    Args:
        theta_x (float): Angle of rotation about the x-axis in radians.
        theta_y (float): Angle of rotation about the y-axis in radians.
        theta_z (float): Angle of rotation about the z-axis in radians.
        degrees (bool, optional): Whether the angles are in degrees (True) or radians (False). Defaults to True.

    Returns:
        np.array: 3x3 rotation matrix. Multiply an R3 vector by this matrix to rotate it about the x, y, and z axes.
    """
    if degrees:
        theta_x = np.deg2rad(theta_x)
        theta_y = np.deg2rad(theta_y)
        theta_z = np.deg2rad(theta_z)

    return (
        rotate_about_x(theta_x, degrees=False)
        @ rotate_about_y(theta_y, degrees=False)
        @ rotate_about_z(theta_z, degrees=False)
    )
