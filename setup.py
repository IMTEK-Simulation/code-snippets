#!/usr/bin/env python

from setuptools import setup, find_packages
import os
import versioneer

__author__     = "Lars Pastewka, Johannes Hörmann"
__copyright__  = "Copyright 2020, IMTEK Simulation, University of Freiburg"
__maintainer__ = "IMTEK Simulation"
__email__      = "johannes.hoermann@imtek.uni-freiburg.de"
__date__       = "Mar 17, 2020"

module_dir = os.path.dirname(os.path.abspath(__file__))

if __name__ == "__main__":
    setup(
        name='imteksimcs',
        version=versioneer.get_version(),
        cmdclass=versioneer.get_cmdclass(),
        description='This repository contains a random collection of code snippets for pre- and postprocessing simulation runs with Ovito, ASE and other tools.',
        long_description=open(os.path.join(module_dir, 'README.md')).read(),
        url='https://github.com/IMTEK-Simulation/code-snippets',
        author='Lars Pastewka, Johannes Hörmann',
        author_email='johannes.hoermann@imtek.uni-freiburg.de',
        license='MIT',
        packages=find_packages(),
        # package_data={'': ['ChangeLog.md']},
        include_package_data=True,
        python_requires='>=3.6.8',
        zip_safe=False,
        install_requires=[
            'numpy>=1.15.2',
            'scipy>=1.2.1',
            'pandas>=0.24.2',
        ],
        extras_require={
            "GMX":  [
                "ParmEd>=3.1.0",
                "GromacsWrapper>=0.8.0",
            ],
            "MPI": [
                "mpi4py>=3.0.1",
            ],
            "NetCDF": [
                "netCDF4>=1.5.1.2",
            ],
            "Ovito": [
                "ovito>=3.1.2",
            ],
            # "Pizza.py": [
            #     "Pizza.py", # not sure how to specify dependency for Pizza.py
            # ],
        },
        scripts=[
            'imteksimcs/cli/LAMMPS_data/merge.py', # requires Pizza.py (Python 2.7)
            'LAMMPS_thermo/extract_thermo.sh',
        ],
        entry_points={
            'console_scripts': [
                'gmx_tools = imteksimcs.cli.GROMACS.gmx_tools:main [GMX]',
                'join_thermo = imteksimcs.cli.LAMMPS_thermo.join_thermo:main',
                'lmp_extract_property = imteksimcs.cli.LAMMPS_data.lmp_extract_property:main [Ovito]',
                'make_self_affine = imteksimcs.cli.Random_Fields.make_self_affine:main',
                'netcdf2data = imteksimcs.cli.NetCDF.netcdf2data:main [Ovito]',
                'ncfilter = imteksimcs.cli.NetCDF.ncfilter:main [GMX,NetCDF,MPI]',
                'ncjoin = imteksimcs.cli.NetCDF.ncjoin:main [NetCDF]', # requires NetCDF4
                'strip_comments = imteksimcs.cli.LAMMPS_data.strip_comments:main',
                'to_hybrid = imteksimcs.cli.LAMMPS_data.to_hybrid:main',
            ],
        },
    )
