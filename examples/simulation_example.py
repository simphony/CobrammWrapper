
from osp.core import cobramm
from osp.wrappers.cobrammwrapper import CobrammSession
from osp.core.utils import pretty_print


# the following is the geometry of a formaldehyde molecule
atoms_to_add = [["C", [0.60720,  0.00000, -0.00040]],
                ["O", [-0.60040, 0.00000,  0.00010]],
                ["H", [1.14720,  0.93530,  0.00160]],
                ["H", [1.14720, -0.93530,  0.00160]]]


# now open a new CobrammSession to execute the simulation
with CobrammSession() as session:
    wrapper = cobramm.wrapper(session=session)

    # Create input ontological object
    system = cobramm.system()
    material = cobramm.material()

    # define the material by specifying atom symbols and cartesian coordinates, atom-by-atom
    for atom_symbol, atom_coords in atoms_to_add:
        atom = material.add(cobramm.atom())
        atom.add(cobramm.position(vector_value=atom_coords, unit="Ang"))
        atom.add(cobramm.element_id(atom_symbol=atom_symbol))

    material.add(cobramm.charge(integer_value=0))
    system.add(material)

    molecule = cobramm.molecule(commercial_name="water")
    solvent = cobramm.solvent()
    solvent.add(molecule)
    system.add(solvent)

    temperature = cobramm.temperature(scalar_value=300, unit="K")
    system.add(temperature)

    pressure = cobramm.pressure(scalar_value=1, unit="bar")
    system.add(pressure)

    # Set accuracy
    accuracy = cobramm.accuracy(accuracy_value=1)

    # Set case
    case = cobramm.case(case_name="absorption")

    system_wrapper = wrapper.add(system, rel=cobramm.has_part)
    accuracy_wrapper = wrapper.add(accuracy, rel=cobramm.has_part)
    case_wrapper = wrapper.add(case, rel=cobramm.has_part)

    # pretty-print the material so far
    print("\n")
    pretty_print(system)
    print("\n")

    # run the simulation
    wrapper.session.run()

    # extract the spectrum
    spectrum = (wrapper.get(oclass=cobramm.system)[0]).get(oclass=cobramm.spectrum)[0]

    # pretty-print the final spectrum
    print("\n")
    pretty_print(spectrum)
    print("\n")
