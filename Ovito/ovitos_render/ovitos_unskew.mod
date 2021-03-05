# unskew the cell and wrap atoms back at PBC for clarity -- require creation of a new modifier

def unskew_xz(frame, input, output):
    # Original simulation cell is passed through by default.
    # Output simulation cell is just a reference to the input cell.
    assert(output.cell is input.cell)

    # Make a copy of the simulation cell:
    cell = output.copy_if_needed(output.cell)

    # copy_if_needed() made a deep copy of the simulation cell object.
    # Now the the input and output each point to different objects.
    assert(cell is output.cell)
    assert(cell is not input.cell)

    # Now it's safe to modify the object copy:
    mat = cell.matrix.copy()
    mat[0,2] = 0.0
    cell.matrix = mat


# Inserting it into the modification pipeline of the node:
node.modifiers.append(PythonScriptModifier(function = unskew_xz))

modifier = WrapPeriodicImagesModifier()
node.modifiers.append(modifier)
