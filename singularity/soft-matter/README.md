Stacked images for building Jupyter notebook allrounder for soft matter stuff.

Refer to Makefile for dependencies.

Build all images with

    make all

Definition file `jupyter_with_gromacs_imteksim.def` requires an archive named

    gmx_top.tar.gz

within this folder. The archive is to containe GROMACS' *top* folder and can be
provided via 

    git clone https://github.com/gromacs/gromacs.git
    cd gromacs/share \
      && tar czvf gmx_top.tar.gz top \
      && mv gmx_top.tar.gz ../.. \
      && cd ../..

The purpose of this mechanism is to allow injecting custom modification on the
force fields into the container image at build time.
