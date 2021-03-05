# compute atomic strain with reference on frame 0 
modstrain = AtomicStrainModifier()
modstrain.reference.load(config.get('files','infile'))
modstrain.eliminate_cell_deformation = True
modstrain.cutoff = 3.5
modstrain.reference_frame = 0
node.modifiers.append(modstrain)

# adding color coding at-strain 0-1 magma coloring
modcolor = ColorCodingModifier(
    particle_property = "Shear Strain",
    gradient = ColorCodingModifier.Magma(), 
	start_value = 0,
	end_value = 1
)
node.modifiers.append(modcolor)

# Select and remove indenter type (5,6)
modifier = SelectParticleTypeModifier(property = "Particle Type")
modifier.types = {5, 6}
node.modifiers.append(modifier)
delmod = DeleteSelectedParticlesModifier()
node.modifiers.append(delmod)
