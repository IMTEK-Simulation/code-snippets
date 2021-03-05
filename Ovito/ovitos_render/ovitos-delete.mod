# Select and remove indenter type (5,6)
modifier = SelectParticleTypeModifier(property = "Particle Type")
modifier.types = {5, 6}
node.modifiers.append(modifier)
delmod = DeleteSelectedParticlesModifier()
node.modifiers.append(delmod)
