import unittest
import json
import math
from collections import OrderedDict

from uuid import UUID
from osp.core import cobramm
from osp.core.utils import pretty_print, export_cuds, import_cuds
from osp.wrappers.cobrammwrapper import CobrammSession
from osp.wrappers.cobrammwrapper import CobrammSimulationEngine

from osp.core.utils import delete_cuds_object_recursively

"""
This test suite is intended to test the session class in two ways:
 * integrations tests related to the interaction with the simulation engine
 * unit test of single methods (part of the wrapper API)
"""
class Session_TestCase(unittest.TestCase):

    def setUp(self):
        """
        Set up input ontological structure with default values.
        Actual values will be setup in the tests.
        """
        
        # Set up material
        material = cobramm.material()

        # Read dummy molecule
        filename = "./data/molecule_dummy.txt"
        with open(filename, "rt") as f:
            content = f.read().splitlines()
        atoms_to_add=[]
        for line in content:
            atoms_to_add.append([line.split(", ")[0],line.split(", ")[1:]])

        # define the material by specifying atom symbols and cartesian coordinates, atom-by-atom
        for atom_symbol, atom_coords in atoms_to_add:
            print("\nAdding\n ", atom_coords)
            atom = material.add(cobramm.atom())
            atom.add(cobramm.position(vector_value=atom_coords, unit="Ang"))
            atom.add(cobramm.element_id(atom_symbol=atom_symbol))

        ##################
        # Set system
        self.system = cobramm.system()

        # system.add(cobramm.spectrum())

        charge = cobramm.charge(integer_value=0)
        material.add(charge)
        self.system.add(material)

        molecule = cobramm.molecule(commercial_name="water")
        solvent = cobramm.solvent()
        solvent.add(molecule)
        self.system.add(solvent)

        temperature = cobramm.temperature(scalar_value=300, unit="K")
        self.system.add(temperature)

        pressure = cobramm.pressure(scalar_value=1, unit="bar")
        self.system.add(pressure)


        # Set accuracy
        self.accuracy = cobramm.accuracy(accuracy_value=1)

        # Set case
        self.case = cobramm.case(case_name="absorption")

    #
    # TEST
    #

    def test_init(self):
        """
        Test init method of the session
        """
        with CobrammSession() as session:
            self.assertTrue(isinstance(session._engine, CobrammSimulationEngine)) 

    def test_str(self):
        """
        Test str method of the session
        """
        with CobrammSession() as session:
            self.assertEqual(str(session), "COBRAMM Wrapper Session")

    def test_load_from_backend(self):
        """
        Test _load_from_backend method of the session
        """
        
        with CobrammSession() as session:
            wrapper = cobramm.wrapper(session=session)

            temperature = cobramm.temperature(scalar_value=300, unit="K")
            fake_uid = UUID(int=8)

            print(list(session._load_from_backend([temperature.uid])))

            self.assertListEqual([None], list(session._load_from_backend([fake_uid])))
            self.assertListEqual([temperature.uid], list(x.uid for x in session._load_from_backend([temperature.uid])))

    def test_run(self):
        """
        Test run method of the session
        """
        
        with CobrammSession() as session:
            wrapper = cobramm.wrapper(session=session)

            # Setup actual values
            self.accuracy.accuracy_value = 1
            self.case.case_name = "absorption"

            self.system.get(oclass=cobramm.temperature)[0].scalar_value = 330
            self.system.get(oclass=cobramm.pressure)[0].scalar_value = 1
            (self.system.get(oclass=cobramm.solvent)[0]).get(oclass=cobramm.molecule)[0].commercial_name = "water"
            (self.system.get(oclass=cobramm.material)[0]).get(oclass=cobramm.charge)[0].integer_value = 0

            wrapper.add(self.system, rel=cobramm.has_part)
            wrapper.add(self.accuracy, rel=cobramm.has_part)
            wrapper.add(self.case, rel=cobramm.has_part)

            wrapper.session.run()

            self.assertTrue(session._ran)

            spectrum = (wrapper.get(oclass=cobramm.system)[0]).get(oclass=cobramm.spectrum)
            print("SPECTRUM ", spectrum)
            self.assertEqual(len(spectrum), 1, "Spectrum not created")
            self.assertEqual(spectrum[0].get(oclass=cobramm.accuracy, rel=cobramm.has_property)[0].accuracy_value, 1, "Wrong accuracy")           



    def test_applyadded(self):
        """
        Test apply_added method of the session
        """

        with CobrammSession() as session:
            wrapper = cobramm.wrapper(session=session)

            # Setup actual values
            self.accuracy.accuracy_value = 1
            self.case.case_name = "absorption"

            self.system.get(oclass=cobramm.temperature)[0].scalar_value = 330
            self.system.get(oclass=cobramm.pressure)[0].scalar_value = 1
            (self.system.get(oclass=cobramm.solvent)[0]).get(oclass=cobramm.molecule)[0].commercial_name = "water"
            (self.system.get(oclass=cobramm.material)[0]).get(oclass=cobramm.charge)[0].integer_value = 0

            wrapper.add(self.system, rel=cobramm.has_part)
            wrapper.add(self.accuracy, rel=cobramm.has_part)
            wrapper.add(self.case, rel=cobramm.has_part)

            wrapper.session.run()

            # Check if apply_added has populated correctly the engine's input
            self.assertEqual(session._engine.input_setup["accuracy"], 1, "Accuracy is wrong")
            self.assertEqual(session._engine.input_setup["case"], "absorption", "Case is wrong")
            self.assertEqual(session._engine.input_setup["charge"], 0, "Charge is wrong")
            self.assertEqual(session._engine.input_setup["temperature"], 330, "Temperature is wrong")
            self.assertEqual(session._engine.input_setup["pressure"], 1, "Pressure is wrong")
            self.assertEqual(session._engine.input_setup["solvent"], "water", "Solvent is wrong")

            atoms = []
            positions = []
            input_atoms = ["C", "O", "H", "H"]
            for uuid in list(session._engine.atoms.keys()):
                atoms.append(session._engine.atoms[uuid][0])
                positions.append(session._engine.atoms[uuid][1])

            self.assertEqual(len(atoms), len(input_atoms), "Wrogn number of atoms")
            for atom in atoms:
                self.assertTrue(atom in input_atoms)
                input_atoms.remove(atom)


    def test_applyupdated(self):
        """
        Test apply_updated method of the session
        """

        with CobrammSession() as session:
            wrapper = cobramm.wrapper(session=session)

            # Setup actual values
            self.accuracy.accuracy_value = 1
            self.case.case_name = "absorption"

            self.system.get(oclass=cobramm.temperature)[0].scalar_value = 330
            self.system.get(oclass=cobramm.pressure)[0].scalar_value = 1
            (self.system.get(oclass=cobramm.solvent)[0]).get(oclass=cobramm.molecule)[0].commercial_name = "water"
            (self.system.get(oclass=cobramm.material)[0]).get(oclass=cobramm.charge)[0].integer_value = 0

            wrapper.add(self.system, rel=cobramm.has_part)
            wrapper.add(self.accuracy, rel=cobramm.has_part)
            wrapper.add(self.case, rel=cobramm.has_part)

            wrapper.session.run()

            (wrapper.get(oclass=cobramm.system)[0]).get(oclass=cobramm.temperature)[0].scalar_value = 420
            ((wrapper.get(oclass=cobramm.system)[0]).get(oclass=cobramm.solvent)[0]).get(oclass=cobramm.molecule)[0].commercial_name = "methanol"

            wrapper.session.run()

            self.assertEqual(session._engine.input_setup["temperature"], 420, "Temperature is wrong")
            self.assertEqual(session._engine.input_setup["solvent"], "methanol", "Solvent is wrong")


    @unittest.expectedFailure
    def test_applydeleted(self):
        """
        Test apply_deleted method of the session
        """ 

        with CobrammSession() as session:
            wrapper = cobramm.wrapper(session=session)

            # Setup actual values
            self.accuracy.accuracy_value = 1
            self.case.case_name = "absorption"

            self.system.get(oclass=cobramm.temperature)[0].scalar_value = 330
            self.system.get(oclass=cobramm.pressure)[0].scalar_value = 1
            (self.system.get(oclass=cobramm.solvent)[0]).get(oclass=cobramm.molecule)[0].commercial_name = "water"
            (self.system.get(oclass=cobramm.material)[0]).get(oclass=cobramm.charge)[0].integer_value = 0

            wrapper.add(self.system, rel=cobramm.has_part)
            wrapper.add(self.accuracy, rel=cobramm.has_part)
            wrapper.add(self.case, rel=cobramm.has_part)

            wrapper.session.run()

            delete_cuds_object_recursively((wrapper.get(oclass=cobramm.system)[0]).get(oclass=cobramm.temperature)[0])
            delete_cuds_object_recursively((wrapper.get(oclass=cobramm.system)[0]).get(oclass=cobramm.pressure)[0])

            wrapper.session.run()


