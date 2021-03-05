# Ovitos modifiers script
# Run polyhedral template matching modifer to find local structure and alloy types.
# Remove all atoms that are not L12 / L10 alloy type
# Color code remaining atoms with blue-white-red of alloy types 2-3-4
# alloy type 2 : L1_0 : Blue
# alloy type 3 : L1_2 (Cu-type) : White
# alloy type 4 : L1_2 (Au-type) : Red

 Run PTM analysis
modptm = PolyhedralTemplateMatchingModifier()
modptm.output_alloy_types = True
node.modifiers.append(modptm)

# Select particles different from L10, L12_Cu and L12_Au
modselection = SelectExpressionModifier(expression = 'AlloyType == 0 || AlloyType == 1 || AlloyType == 5')
node.modifiers.append(modselection)

# Delete selected particles
delmod = DeleteSelectedParticlesModifier()
node.modifiers.append(delmod)

# Color code remaining particles
modcolor = ColorCodingModifier(
    particle_property = "Alloy Type",
    gradient = ColorCodingModifier.BlueWhiteRed(), 
	start_value = 2,
	end_value = 4
)
node.modifiers.append(modcolor)
