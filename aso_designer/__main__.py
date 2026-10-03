"""Entry point for `python -m aso_designer`.

Allows running ASODesigner as a module:
    python -m aso_designer design --target BRCA1 --chemistry gapmer
"""

import sys

from aso_designer.cli import main

if __name__ == "__main__":
    sys.exit(main())
