import argparse

def main_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="AmberTI",
        # formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="\n"
    )
    subparsers = parser.add_subparsers(title="Valid subcommands", dest="command")
    parser_maketop = subparsers.add_parser(
        "maketop",
        help="generate topology from sdf",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser_maketop.add_argument("MOL")
    parser_maketop.add_argument("--forcefield", "-f", default="gaff2")
    parser_maketop.add_argument("--charge", "-c", default="bcc")
    parser_maketop.add_argument("--name", "-n", default="MOL")
    parser_maketop.add_argument("--charge-number", "-nc", default=0)
    
    parser_run = subparsers.add_parser(
        "run",
        help="run FEP!",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser_run.add_argument("protein")
    parser_run.add_argument("lig1")
    parser_run.add_argument("lig2")

    return parser

def parse_args(args=None):
    parser = main_parser()

    parsed_args = parser.parse_args(args=args)
    if parsed_args.command is None:
        parser.print_help()

    return parsed_args



if __name__ == "__main__":

    args = parse_args()
