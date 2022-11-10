#!/bin/bash
set -e

# delete and rebuild matplotlib font cache
# rm -rf /home/${NB_USER}/.cache/matplotlib
fc-cache --force
PATH=$(echo "$PATH" | sed -e 's|/opt/conda[^:]*:||g') python3 -c 'import matplotlib.font_manager; matplotlib.font_manager._load_fontmanager(try_read_cache=False)'

echo
echo "matplotlib cache dir:"
echo
PATH=$(echo "$PATH" | sed -e 's|/opt/conda[^:]*:||g') python3 -c "import matplotlib; print(matplotlib.get_cachedir())"
echo
echo "matplotlib fonts:"
echo
PATH=$(echo "$PATH" | sed -e 's|/opt/conda[^:]*:||g') python3 -c "import matplotlib.font_manager ; [print(f) for f in matplotlib.font_manager.findSystemFonts(fontpaths=None)]"
echo
exec "$@"
