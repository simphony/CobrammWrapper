import unittest
import json
import math
from collections import OrderedDict

from osp.core import cobramm
from osp.core.utils import pretty_print, export_cuds, import_cuds
from osp.wrappers.cobrammwrapper import CobrammSession

"""
This test suite is intended to test the end to end interaction
application <--> wrapper <--> simulation engine.
"""
class EndToEnd_TestCase(unittest.TestCase):

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
    # Utility method for testing
    #

    def read_spectrum_file(self, filename):
        """
        Read a file containing the spectrum and returns a dictionary
        """

        lines = []
        spectrum = {}
        with open(filename) as file:
            lines = [line.rstrip() for line in file]

        for line in lines:
            wave, intensity = line.split(" ")
            spectrum[float(wave)] = float(intensity)

        # with open('./data/expected_spectrum.txt', 'w') as file:
        #     file.write(json.dumps(spectrum))

        return spectrum

    def read_spectrum_ontobj(self, ontobj):
        """
        Read a specturm ontological object and returns a dictionary
        """

        spectrum = {}

        bins = ontobj.get(oclass=cobramm.bin) 
        for bin in bins:
            intensity = bin.get(oclass=cobramm.intensity)[0].scalar_value
            wavelength = bin.get(oclass=cobramm.wavelength)[0].scalar_value
            spectrum[float(wavelength)] = float(intensity)

        # with open('./data/actual_spectrum.txt', 'w') as file:
        #     file.write(json.dumps(spectrum))

        return spectrum

    def get_spectrum_byaccuracy(self, level, ontobj):
        """
        Filter a list of spectrum by accuracy level
        """
        spectra = ontobj.get(oclass=cobramm.spectrum)
        for spectrum in spectra:
            if spectrum.get(oclass=cobramm.accuracy, rel=cobramm.has_property)[0].accuracy_value == level:
                print("\nSELECTED SPECTRUM\n")
                pretty_print(spectrum.get(oclass=cobramm.accuracy, rel=cobramm.has_property)[0])
                print("\n")
                return spectrum

        return cobramm.spectrum()

    #
    # TEST
    #

    def test_calcA(self):
        """
        First simulation - test level 1 accuracy
        absorption 1 330 1 water 0 molecule_dummy.txt
        """

        spectrum_file = "./data/CALC_A/spectrum_nm_level_1.dat"

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

            
            print("\n")
            pretty_print(wrapper)
            print("\n")

            wrapper.session.run()

            #(wrapper.get(oclass=cobramm.system)[0]).get(oclass=cobramm.spectrum)[0]
            actual_spectrum = self.read_spectrum_ontobj(self.get_spectrum_byaccuracy(1, (wrapper.get(oclass=cobramm.system)[0])))
            # actual_spectrum = self.read_spectrum_file(spectrum_file)
            expected_spectrum = self.read_spectrum_file(spectrum_file)

            print("Comparing...")
            self.assertEqual(len(list(actual_spectrum.keys())), len(list(expected_spectrum.keys())))
            for (k1,v1) in expected_spectrum.items():
                self.assertTrue(k1 in actual_spectrum)
                actual_value = actual_spectrum[k1]
                self.assertTrue(math.isclose(abs(v1 - actual_value), 0, abs_tol=1e-2))


    def test_calcB(self):
        """
        Second simulation - test level 2 accuracy
        absorption 2 330 1 methanol 0 molecule_dummy.txt
        """
        
        spectrum_file = "./data/CALC_B/spectrum_nm_level_2.dat"

        with CobrammSession() as session:
            wrapper = cobramm.wrapper(session=session)

            # Setup actual values
            self.accuracy.accuracy_value = 2
            self.case.case_name = "absorption"

            self.system.get(oclass=cobramm.temperature)[0].scalar_value = 330
            self.system.get(oclass=cobramm.pressure)[0].scalar_value = 1
            (self.system.get(oclass=cobramm.solvent)[0]).get(oclass=cobramm.molecule)[0].commercial_name = "methanol"
            (self.system.get(oclass=cobramm.material)[0]).get(oclass=cobramm.charge)[0].integer_value = 0

            wrapper.add(self.system, rel=cobramm.has_part)
            wrapper.add(self.accuracy, rel=cobramm.has_part)
            wrapper.add(self.case, rel=cobramm.has_part)

            
            print("\n")
            pretty_print(wrapper)
            print("\n")

            wrapper.session.run()

            actual_spectrum = self.read_spectrum_ontobj(self.get_spectrum_byaccuracy(2, (wrapper.get(oclass=cobramm.system)[0])))
            expected_spectrum = self.read_spectrum_file(spectrum_file)

            print("Comparing ...")
            self.assertEqual(len(list(actual_spectrum.keys())), len(list(expected_spectrum.keys())))
            for (k1,v1) in expected_spectrum.items():
                self.assertTrue(k1 in actual_spectrum)
                actual_value = actual_spectrum[k1]
                self.assertTrue(math.isclose(abs(v1 - actual_value), 0, abs_tol=1e-2))

    def test_calcC(self):
        """
        Third simulation - test both level 1 and level 2 accuracy
        emission 2 250 1 water 0 molecule_dummy.txt
        """
        
        spectrum_file_1 = "./data/CALC_C/spectrum_nm_level_1.dat"
        spectrum_file_2 = "./data/CALC_C/spectrum_nm_level_2.dat"

        with CobrammSession() as session:
            wrapper = cobramm.wrapper(session=session)

            # Setup actual values
            self.accuracy.accuracy_value = 2
            self.case.case_name = "emission"

            self.system.get(oclass=cobramm.temperature)[0].scalar_value = 250
            self.system.get(oclass=cobramm.pressure)[0].scalar_value = 1
            (self.system.get(oclass=cobramm.solvent)[0]).get(oclass=cobramm.molecule)[0].commercial_name = "water"
            (self.system.get(oclass=cobramm.material)[0]).get(oclass=cobramm.charge)[0].integer_value = 0

            wrapper.add(self.system, rel=cobramm.has_part)
            wrapper.add(self.accuracy, rel=cobramm.has_part)
            wrapper.add(self.case, rel=cobramm.has_part)

            
            print("\n")
            pretty_print(wrapper)
            print("\n")

            wrapper.session.run()

            actual_spectrum = self.read_spectrum_ontobj(self.get_spectrum_byaccuracy(1, (wrapper.get(oclass=cobramm.system)[0])))
            expected_spectrum = self.read_spectrum_file(spectrum_file_1)

            print("Comparing lev 1 ...")
            self.assertEqual(len(list(actual_spectrum.keys())), len(list(expected_spectrum.keys())))
            for (k1,v1) in expected_spectrum.items():
                self.assertTrue(k1 in actual_spectrum)
                actual_value = actual_spectrum[k1]
                self.assertTrue(math.isclose(abs(v1 - actual_value), 0, abs_tol=1e-2))

            actual_spectrum = self.read_spectrum_ontobj(self.get_spectrum_byaccuracy(2, (wrapper.get(oclass=cobramm.system)[0])))
            expected_spectrum = self.read_spectrum_file(spectrum_file_2)

            print("Comparing lev 2...")
            self.assertEqual(len(list(actual_spectrum.keys())), len(list(expected_spectrum.keys())))
            for (k1,v1) in expected_spectrum.items():
                self.assertTrue(k1 in actual_spectrum)
                actual_value = actual_spectrum[k1]
                self.assertTrue(math.isclose(abs(v1 - actual_value), 0, abs_tol=1e-2))