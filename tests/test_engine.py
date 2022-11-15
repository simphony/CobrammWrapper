import unittest
import json
import math
from collections import OrderedDict

from osp.core import cobramm
from osp.core.utils import pretty_print, export_cuds, import_cuds
from osp.wrappers.cobrammwrapper import CobrammSimulationEngine

"""
This test suite is intended to unit test the simulation engine:
"""
class Session_TestCase(unittest.TestCase):

    def setUp(self):
        """
        Set up input ontological structure with default values.
        Actual values will be setup in the tests.
        """

        self.simu_engine = CobrammSimulationEngine()
        
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
        Test the init method of the simulation engine
        """

        self.assertFalse(self.simu_engine.executed, "Executed is true")
        self.assertFalse(self.simu_engine.atoms, "Atoms are not empty")
        self.assertFalse(self.simu_engine.input_setup, "Input is not empty")
        self.assertEqual(self.simu_engine.accuracy, 2, "Accuracy is not 2")
        self.assertEqual(self.simu_engine.case, "transient", "Case is not transient")
        self.assertFalse(self.simu_engine.grid1, "Grid1 is not empty") 
        self.assertFalse(self.simu_engine.spectrum1, "Spectrum1 is not empty") 
        self.assertFalse(self.simu_engine.grid2, "Grid2 is not empty") 
        self.assertFalse(self.simu_engine.spectrum2, "Spectrum2 is not empty") 
        self.assertFalse(self.simu_engine.grid3, "Grid3 is not empty") 
        self.assertFalse(self.simu_engine.spectrum3, "Spectrum3 is not empty") 

    def test_str(self):
        """
        Test the str method of the simulation engine
        """
        self.assertEqual(str(self.simu_engine), "COBRAMM Engine Connection")

    def test_addatom(self):
        """
        Test the add_atom method of the simulation engine
        """

        uuid = "6ae7082c-2fbc-4871-951d-a39f501a8a47"
        self.simu_engine.add_atom(uuid, "H", ["1.14720",  "0.93530",  "0.00160"])

        self.assertTrue(uuid in self.simu_engine.atoms, "Missing atom")
        self.assertEqual(self.simu_engine.atoms[uuid][0], "H", "Wrong symbol")
        self.assertEqual(self.simu_engine.atoms[uuid][1], ["1.14720",  "0.93530",  "0.00160"], "Wrong positions")


    def test_deleteatom(self):
        """
        Test the delete_atom method of the simulation engine
        """

        uuid = "6ae7082c-2fbc-4871-951d-a39f501a8a47"
        self.simu_engine.add_atom(uuid, "H", ["1.14720",  "0.93530",  "0.00160"])
        self.simu_engine.delete_atom(uuid)

        self.assertFalse(uuid in self.simu_engine.atoms, "Atom not deleted")

    def test_addproperty(self):
        """
        Test the add_property method of the simulation engine
        """

        self.simu_engine.add_property("temperature", 400)

        self.assertTrue("temperature" in self.simu_engine.input_setup, "Missing temperature")
        self.assertEqual(self.simu_engine.input_setup["temperature"], 400, "Wrong temperature")

    def test_deleteproperty(self):
        """
        Test the delete_property method of the simulation engine
        """

        self.simu_engine.add_property("temperature", 400)
        self.simu_engine.delete_property("temperature")

        self.assertFalse("temperature" in self.simu_engine.input_setup, "Temperature not deleted")

    def test_run_getspectrum(self):
        """
        Test both the run and get_spectrum methods of the simulation engine
        """

        self.assertIsNone(self.simu_engine.get_spectrum(), "Spectrum is not None")

        self.simu_engine.add_property("case", "absorption")
        self.simu_engine.add_property("accuracy", 1)
        self.simu_engine.add_property("temperature", 330)
        self.simu_engine.add_property("pressure", 1)
        self.simu_engine.add_property("solvent", "water")
        self.simu_engine.add_property("charge", 0)
        self.simu_engine.add_atom("5f7afbe1-d112-4443-8726-78b5fb50fa0a", "C", [0.60720,  0.00000, -0.00040])
        self.simu_engine.add_atom("abf35cde-374e-4477-96e1-54adbcff68ec", "O", [-0.60040, 0.00000,  0.00010])
        self.simu_engine.add_atom("62750a1c-7049-4783-a315-68e3e8a2fb28", "H", [1.14720,  0.93530,  0.00160])
        self.simu_engine.add_atom("8c38c263-93d5-4ac0-aa6f-8eb77ce37ffb", "H", [1.14720, -0.93530,  0.00160])

        self.simu_engine.run()

        spectrum = self.simu_engine.get_spectrum()

        self.assertIsNotNone(spectrum)
        self.assertTrue(spectrum[0])
        self.assertTrue(spectrum[1])
        self.assertEqual(len(spectrum[0]),len(spectrum[1]))
        self.assertFalse(spectrum[2]) 
        self.assertFalse(spectrum[3]) 
        self.assertFalse(spectrum[4]) 
        self.assertFalse(spectrum[5])



        