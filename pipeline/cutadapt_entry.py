"""Console entry point frozen as the Cutadapt executable shipped with the app."""

from importlib import import_module
import sys


if __name__ == "__main__":
    import_module("cutadapt.cli").main(sys.argv[1:])
