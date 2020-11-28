.RECIPEPREFIX := $(.RECIPEPREFIX) 
SHELL := /bin/bash

def = $(wildcard *.def)
sif = $(def:.def=.sif)

all: $(sif)

openmpi.sif:
    sudo -E singularity build openmpi.sif openmpi.def

jupyter_base.sif: openmpi.sif
    sudo -E singularity build jupyter_base.sif jupyter_base.def

jupyter_nbextensions.sif: jupyter_base.sif
    sudo -E singularity build jupyter_nbextensions.sif jupyter_nbextensions.def

jupyter_imteksim_jotelha.sif: jupyter_nbextensions.sif
    sudo -E singularity build jupyter_imteksim_jotelha.sif jupyter_imteksim_jotelha.def

jupyter_with_gromacs_imteksim.sif: jupyter_imteksim_jotelha.sif
    sudo -E singularity build jupyter_with_gromacs_imteksim.sif jupyter_with_gromacs_imteksim.def

.PHONY: clean
clean:
    rm -f $(sif)
