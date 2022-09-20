#!/bin/bash
set -e

fc-cache --force
echo
echo "matplotlib cache dir:"
echo
python3 -c "import matplotlib; print(matplotlib.get_cachedir())"
echo
echo "matplotlib fonts:"
echo
python3 -c "import matplotlib.font_manager ; [print(f) for f in matplotlib.font_manager.findSystemFonts(fontpaths=None)]"
echo
exec "$@"
