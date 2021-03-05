# unskew the cell, shift cell origine half-z up and wrap atoms back at PBC for clarity and get velocity = 0 at the median plane -- require creation of a new modifier

def unskew_xz(frame, input, output):
    global home, ID
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
    
    # Using the marray to access a modifiable array containing the initial input information
    pos = input.particle_properties.position.marray
    
    fn = home+"/flip_{0}.tmp".format(ID)
    try:
        flip = np.memmap(fn, dtype='float32', mode='r+',shape=(2)) 
    except:
        flip = np.memmap(fn, dtype='float32', mode='w+',shape=(2))
        
    if frame == 0:
        flip[0] = 0    
    elif mat[0,2] < 0 and flip[1] > 0:
        flip[0] += 1 
    
    pos[:,0] -= mat[0,2]/2 + flip[0]*mat[0,0]/2
    ## Inform pipeline that input has changed.
    ## Failing to do so would lead to incorrect results below. OVITO would assume the 
    ## cached pipeline output is  still valid and wouldn't re-evaluate the modifiers.
    #input.particle_properties.position.changed()
    
    flip[1] = mat[0,2]
    mat[0,2] = 0.0
    cell.matrix = mat
