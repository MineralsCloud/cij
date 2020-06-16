import argparse
import logging

def run(config_fname: str):
    import cij.core.calculator
    calculator = cij.core.calculator.Calculator(config_fname)
    calculator.write_output()

def parse_args():
    import cij

    parser = argparse.ArgumentParser()

    parser.add_argument(dest="config_filename", metavar="CONF.yml", action="store", help="name of the configuration file")
    parser.add_argument('-V', '--version', action='version', version=f'Cij v.{cij.__version__}')

    return parser.parse_args()

def main():
    parsed = parse_args()
    run(parsed.config_filename)