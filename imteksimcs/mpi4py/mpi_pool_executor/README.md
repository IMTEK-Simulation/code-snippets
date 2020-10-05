# Run mpi4py-parallel Python functions from command line without overhead.

## Examples

Serial function call with

```console with
$ python -m imteksimcs.mpi4py.mpi_pool_executor --debug print "Hello world!"
INFO:root:This is '/mnt/home/git/code-snippets/imteksimcs/mpi4py/mpi_pool_executor/run.py' : 'imteksimcs.mpi4py.mpi_pool_executor.run'.
INFO:root:Argument args: '['Hello world!']' (<class 'list'>)
INFO:root:Argument callable: 'print' (<class 'str'>)
INFO:root:Argument debug: 'True' (<class 'bool'>)
INFO:root:Argument log: 'None' (<class 'NoneType'>)
INFO:root:Argument verbose: 'False' (<class 'bool'>)
DEBUG:imteksimcs.mpi4py.mpi_pool_executor.mpi_pool_executor:MPI version: (3, 1)
DEBUG:imteksimcs.mpi4py.mpi_pool_executor.mpi_pool_executor:Current MPI size is: 1
DEBUG:imteksimcs.mpi4py.mpi_pool_executor.mpi_pool_executor:MPI universe size is: 1
Hello world!
```

or statically parallel with

```
$ mpirun -n 4 python -m mpi4py.futures -m imteksimcs.mpi4py.mpi_pool_executor --debug print "Hello world!"
INFO:root:This is '/mnt/home/git/code-snippets/imteksimcs/mpi4py/mpi_pool_executor/run.py' : 'imteksimcs.mpi4py.mpi_pool_executor.run'.
INFO:root:Argument args: '['Hello world!']' (<class 'list'>)
INFO:root:Argument callable: 'print' (<class 'str'>)
INFO:root:Argument debug: 'True' (<class 'bool'>)
INFO:root:Argument log: 'None' (<class 'NoneType'>)
INFO:root:Argument verbose: 'False' (<class 'bool'>)
DEBUG:imteksimcs.mpi4py.mpi_pool_executor.mpi_pool_executor:MPI version: (3, 1)
DEBUG:imteksimcs.mpi4py.mpi_pool_executor.mpi_pool_executor:Current MPI size is: 4
DEBUG:imteksimcs.mpi4py.mpi_pool_executor.mpi_pool_executor:MPI universe size is: 1
Hello world!
Hello world!
Hello world!
Hello world!
```


Instead of wrapping a call like

```python
from imteksimcs.mpi4py.mpi_pool_executor import call
from imteksimcs.GROMACS.gmx_mpi_rdf import atom_atom_rdf

call(atom_atom_rdf, gro='default.gro', trr='default.trr',
     out='OW_NA_rdf_parallel.txt', atom_name_a='OW', atom_name_b='NA')
```

into some script `run.py` for execution via

```bash
mpirun -np 4 python -m mpi4py.futures test.py
```

run the desired `atom_atom_rdf` callable directly MPI-parallel via

```bash
mpirun -n 4 python -m mpi4py.futures -m imteksimcs.mpi4py.mpi_pool_executor --debug \
    imteksimcs.GROMACS.gmx_mpi_rdf.atom_atom_rdf default.gro default.trr OW_NA_rdf_parallel.txt OW NA
```