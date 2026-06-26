import argparse
import os
import sys
import pyJMAK


def main():
    parser = argparse.ArgumentParser(
        description="Run phase transformation simulation using pyJMAK"
    )
    parser.add_argument(
        "-t", "--tempfile", required=True,
        help="Temperature file (time vs temperature)"
    )
    parser.add_argument(
        "-i", "--inputfile", required=True,
        help="Abaqus-formatted phase transformation input file"
    )
    parser.add_argument(
        "-e", "--exp", nargs="+", metavar=("PHASE_NAME", "FILE"),
        help="Pairs of phase name and experimental data file. Example: -e Austenite austenite.txt Martensite martensite.txt"
    )
    parser.add_argument(
        "-o", "--output", default=None,
        help="Output figure filename (PNG)"
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true",
        help="Enable verbose output"
    )

    args = parser.parse_args()

    # Change working dir to current script location
    abspath = os.path.abspath(__file__)
    path = os.path.dirname(abspath)
    os.chdir(path)

    # Create model
    model = pyJMAK.jmak(args.verbose)
    model.loadInputFiles(args.tempfile, args.inputfile)
    model.run_phase_trans_sim()

    # Experimental comparison if provided
    if args.exp and len(args.exp) % 2 == 0:
        exp_data_files = {
            args.exp[i]: args.exp[i + 1] for i in range(0, len(args.exp), 2)
        }
        model.setExperimentalComparison(True, exp_data_files)

    # Plot results
    model.plot_results(filename=args.output if args.output else None)

    if args.verbose:
        print("Phase names:", model.ph_names)


if __name__ == "__main__":
    main()
