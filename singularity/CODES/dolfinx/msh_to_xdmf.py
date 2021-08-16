import meshio
import sys
import numpy

mesh3D_from_msh = meshio.read(sys.argv[1]+".msh")

def create_mesh(mesh, cell_type, prune_z=False):
    cells = mesh.get_cells_type(cell_type)
    cell_data = mesh.get_cell_data("gmsh:physical", cell_type)
    out_mesh = meshio.Mesh(points=mesh.points, cells={cell_type: cells}, cell_data={"name_to_read":[cell_data]})
    if prune_z:
        out_mesh.prune_z_0()
    return out_mesh


tetra_mesh = create_mesh(mesh3D_from_msh, "tetra")
meshio.write(sys.argv[1]+".xdmf", tetra_mesh)
