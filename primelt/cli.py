import argparse
from .core import run_primelt3, run_primarsmelt

def main():
    parser = argparse.ArgumentParser(description="PRIMELT3-P / PRIMARSMELT CLI Tool")
    parser.add_argument("mode", choices=["primelt3", "primarsmelt"], help="Calculation mode")
    parser.add_argument("input_file", help="Input CSV file")
    parser.add_argument("--mgo", type=float, help="MgO value for source")
    parser.add_argument("--feo", type=float, help="FeO value for source")
    parser.add_argument("--output", default=None, help="Output directory")
    args = parser.parse_args()

    if args.mode == "primelt3":
        run_primelt3(args.input_file, args.mgo, args.feo, args.output)
    else:
        run_primarsmelt(args.input_file, args.mgo, args.feo, args.output)

if __name__ == "__main__":
    main()