import numpy as np
from osp.core import cobramm
from osp.wrappers.cobrammwrapper import CobrammSession
from osp.core.utils import pretty_print

material = cobramm.material()

# the following is the geometry of a formaldehyde molecule
atoms_to_add = [["C", [0.60720,  0.00000, -0.00040]],
                ["O", [-0.60040, 0.00000,  0.00010]],
                ["H", [1.14720,  0.93530,  0.00160]],
                ["H", [1.14720, -0.93530,  0.00160]]]

# define the material by specifying atom symbols and cartesian coordinates, atom-by-atom
for atom_symbol, atom_coords in atoms_to_add:
    atom = material.add(cobramm.atom())
    atom.add(cobramm.position(vector_value=atom_coords, unit="Ang"))
    atom.add(cobramm.element_id(atom_symbol=atom_symbol))

# add other physical properties: temperature, pressure, solvent (with commercial name)
material.add(cobramm.temperature(scalar_value=300., unit="K"))
material.add(cobramm.pressure(scalar_value=1., unit="bar"))
material.add(cobramm.solvent(commercial_name="water"))

# pretty-print the material so far
pretty_print(material)

# now open a new CobrammSession to execute the simulation
with CobrammSession() as session:
    wrapper = cobramm.wrapper(session=session)
    material_wrapper = wrapper.add(material, rel=cobramm.has_part)

    material_wrapper.session.run()
#
#     print("-" * 100)
#
#     atom = material_wrapper.get(rel=cobramm.has_part)[0]
#     print(atom)
#     atom.get(oclass=cobramm.position)[0].vector_value = [0, 0, 0]
#     wrapper.session.run()
#
#     print("-" * 100)
#
#     pretty_print(material_wrapper)