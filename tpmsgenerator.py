import argparse
from spartantpms.stl import generate_stl
from spartantpms.optimization_functions import gyroid_box, diamond_box, sheet_gyroid_box

parser = argparse.ArgumentParser(
    prog="SPARTAN TPMS Generator",
    description="Generate STL files for gyroid or diamond TPMS structures for the SPARTAN project.",
)
parser.add_argument("filename", type=str, help="Name of the stl file to export")
parser.add_argument(
    "tpms",
    type=str,
    help="Type of TPMS to generate (`gyroid`, `sheet_gyroid`, or `diamond`)",
)
parser.add_argument("lambda_x", type=float, help="TPMS wavelength in the x dimension")
parser.add_argument("lambda_y", type=float, help="TPMS wavelength in the y dimension")
parser.add_argument("lambda_z", type=float, help="TPMS wavelength in the z dimension")
parser.add_argument(
    "theta_x",
    type=float,
    help="TPMS lattice rotation (degrees) about the x axis",
)
parser.add_argument(
    "theta_y",
    type=float,
    help="TPMS lattice rotation (degrees) about the y axis",
)
parser.add_argument(
    "theta_z",
    type=float,
    help="TPMS lattice rotation (degrees) about the z axis",
)
parser.add_argument(
    "porosity",
    type=float,
    help="Porosity of the lattice (0-1). This is the fraction of free space in the structure.",
)
parser.add_argument(
    "-s",
    "--step_size",
    dest="step_size",
    type=float,
    help="Step size at which to render the stl. This is in the same units as the `lambda_x,y,z` values. Smaller values yield higher resolution structures but increase computational cost. Default value is 0.2.",
    required=False,
    default=0.2,
)
parser.add_argument(
    "-n",
    "--num_periods",
    dest="num_periods",
    type=int,
    help="This determines the overall size of the test structure as a multiple of `lambda_x,y,z`. This is the number of periods/unit cells to generate. Default value is 4.",
    required=False,
    default=4,
)

parser.add_argument(
    "-nx",
    "--num_x_periods",
    dest="num_x_periods",
    type=int,
    help="This determines the overall size of the test structure in the x dimension as a multiple of `lambda_x`. If you provide this value, you must also provide `ny` and `nz`.",
    required=False,
    default=None,
)

parser.add_argument(
    "-ny",
    "--num_y_periods",
    dest="num_y_periods",
    type=int,
    help="This determines the overall size of the test structure in the y dimension as a multiple of `lambda_y`. If you provide this value, you must also provide `nx` and `nz`.",
    required=False,
    default=None,
)

parser.add_argument(
    "-nz",
    "--num_z_periods",
    dest="num_z_periods",
    type=int,
    help="This determines the overall size of the test structure in the z dimension as a multiple of `lambda_z`. If you provide this value, you must also provide `nx` and `ny`.",
    required=False,
    default=None,
)

if __name__ == "__main__":
    args = parser.parse_args()

    if args.tpms == "gyroid":
        tpms_function = gyroid_box
    elif args.tpms == "diamond":
        tpms_function = diamond_box
    elif args.tpms == "sheet_gyroid":
        tpms_function = sheet_gyroid_box
    else:
        raise ValueError(
            f"The `tpms` argument must be either `gyroid` or `diamond` -- you provided {args.tpms}"
        )

    for arg in ["lambda_x", "lambda_y", "lambda_z"]:
        val = getattr(args, arg)
        if val < 0:
            raise ValueError(
                f"The `{arg}` argument must be greater than 0 -- you provided {val}"
            )

    for arg in ["theta_x", "theta_y", "theta_z"]:
        val = getattr(args, arg)
        if val < 0 or val > 90:
            raise ValueError(
                f"The `{arg}` argument must be between 0 and 90 -- you provided {val}"
            )

    if args.porosity < 0 or args.porosity > 1:
        raise ValueError(
            f"The `porosity` argument must be between 0 and 1 -- you provided {args.porosity}"
        )

    if any([args.num_x_periods, args.num_y_periods, args.num_z_periods]):
        if not all([args.num_x_periods, args.num_y_periods, args.num_z_periods]):
            raise ValueError(
                "If you provide any of `num_x_periods`, `num_y_periods`, or `num_z_periods`, you must provide all three."
            )
        if args.num_periods != 4:
            raise ValueError(
                "If you provide any of `num_x_periods`, `num_y_periods`, or `num_z_periods`, you must not provide `num_periods`."
            )
        num_periods = (args.num_x_periods, args.num_y_periods, args.num_z_periods)
    else:
        num_periods = (args.num_periods, args.num_periods, args.num_periods)

    sdf = tpms_function(
        args.lambda_x,
        args.lambda_y,
        args.lambda_z,
        args.theta_x,
        args.theta_y,
        args.theta_z,
        args.porosity,
        n_periods=num_periods,
    )
    generate_stl(
        f=sdf,
        fpath=args.filename,
        step=args.step_size,
        bounds=None,
        verbose=False,
        clean=True,
    )
