# Select and remove indenter type (5,6)
modifier = SelectParticleTypeModifier(property = "Particle Type")
modifier.types = {5, 6}
node.modifiers.append(modifier)
delmod = DeleteSelectedParticlesModifier()
node.modifiers.append(delmod)

# Apply DXA
moddxa = DislocationAnalysisModifier()
moddxa.input_crystal_structure = DislocationAnalysisModifier.Lattice.FCC
node.modifiers.append(moddxa)

# reduce defect mesh opacity
node.compute()
node.output.surface.display.enabled = True
node.output.surface.display.surface_color = (1.0, 0.70588, 0.63137)
node.output.surface.display.surface_transparency = 0.9

# remove particles rendering
node.source.particle_properties.position.display.enabled = False
