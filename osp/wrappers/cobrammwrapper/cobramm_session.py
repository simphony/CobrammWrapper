# noinspection PyUnresolvedReferences
from osp.core import cobramm
from osp.core.session import SimWrapperSession
from osp.wrappers.cobrammwrapper.simulation_engine import CobrammSimulationEngine
from osp.core.utils import pretty_print
from pprint import pprint
import csv


class CobrammSession(SimWrapperSession):
    """
    Session class for some engine.
    """

    def __init__(self, engine=None, **kwargs):
        super().__init__(engine or CobrammSimulationEngine(), **kwargs)

    def __str__(self):
        return "COBRAMM Wrapper Session"

    # OVERRIDE
    def _run(self, root_cuds_object):
        """Call the run command of the engine."""
        self._engine.run()

        case = (root_cuds_object.get(oclass=cobramm.case)[0]).case_name
        accuracy = (root_cuds_object.get(oclass=cobramm.accuracy)[0]).accuracy_value
        print("Case is " + case)

       
        if(case == "transient"):
            for accuracy_level in range(1, accuracy + 1):
                res = {}
                idx = -1
                interval = 19             
                with open("./results/PPheat_2D_data_level_{}/data_PPheatmap_t2.txt".format(accuracy_level)) as csv_file:
                    csv_reader = csv.reader(csv_file, delimiter='\t')
                    for row in csv_reader:
                        if idx % interval == 0 and len(row) == 3:
                            print("Analyzing {}".format(row))
                            if row[0] not in res.keys():
                                res[row[0]] = []

                            res[row[0]].append((float(row[1]), float(row[2])))

                        idx = idx + 1

                print("KEYS {}".format(res.keys()))

                # res = [el for el in res if len(el) == 3]
                # res = [[[float(el[0]), float(el[1]), float(el[2])] for el in entry] for entry in res]

                if(len(res.keys()) > 0):
                    for time in res:
                        i = 0
                        value = res[time]
                        spectrum = (root_cuds_object.get(oclass=cobramm.system)[0]).add(cobramm.spectrum())
                        spectrum.add(cobramm.accuracy(accuracy_value=accuracy_level), rel=cobramm.has_property)
                        spectrum.add(cobramm.time(scalar_value=time, unit="ns"))
                        for el in value:
                            print(el, i)
                            bin = spectrum.add(cobramm.bin())
                            # spectrum.add(cobramm.time(scalar_value=el[0], unit="ns"))
                            bin.add(cobramm.intensity(scalar_value=el[0], unit="arbitrary"))
                            bin.add(cobramm.wavelength(scalar_value=el[1], unit="nm"))
                            i = i + 1
        else:
            n = 2
            spc = list(self._engine.get_spectrum())

            if spc is not None:
                acc = 0
                for i in range(0, len(spc), n):
                    chunk = spc[i:i+n]
                    zipped = list(zip(*chunk))
                    acc = acc + 1

                    if len(zipped) > 0:
                        spectrum = (root_cuds_object.get(oclass=cobramm.system)[0]).add(cobramm.spectrum())
                        spectrum.add(cobramm.accuracy(accuracy_value=acc), rel=cobramm.has_property)
                        for wav, ints in zipped:
                            bin = spectrum.add(cobramm.bin())
                            bin.add(cobramm.intensity(scalar_value=ints, unit="arbitrary"))
                            bin.add(cobramm.wavelength(scalar_value=wav, unit="nm"))

    # OVERRIDE
    def _load_from_backend(self, uids, expired=None):
        """Loads the cuds object from the simulation engine"""
        for uid in uids:
            if uid in self._registry:
                yield self._registry.get(uid)
            else:
                yield None

        # print("ROOT {}".format(self._registry.get(self.root)))
        # print("REGISTRY {}".format(self._registry))
        # print("LEN {}".format(len(self._registry)))
        # print("UIDS {}".format(uids))
        # for uid in uids:
        #     if uid in self._registry:
        #         obj = self._registry.get(uid)
        #         root = self._registry.get(self.root)

        #         print(type(obj))

        #         if obj.is_a(cobramm.system):
        #             print("Found system")
        #             pretty_print(self._registry.get(self.root))
        #             print("\n")
        #             pretty_print(obj)
        #             spc = self._engine.get_spectrum()
        #             print("Check OBJ {}".format(obj.get(oclass=cobramm.spectrum)))
        #             print("Check ROOT {}".format(root.get(oclass=cobramm.spectrum)))
        #             if obj.get(oclass=cobramm.spectrum):
        #                 print("Found existing spectrum...removing")
        #                 obj.remove(oclass=cobramm.spectrum)
        #             if spc is not None:
        #                 print("Adding new spectrum")
        #                 obj.add(cobramm.spectrum())

        #         # if obj.is_a(cobramm.spectrum):
        #         #     spc = self._engine.get_spectrum()
        #         #     # print("Found spectrum {} ".format(spc))
        #         #     # V5
        #         #     if spc is not None:
        #         #         # print("Entered {}".format(list(zip(*spc))))
        #         #         if obj.get(oclass=cobramm.intensity):
        #         #             obj.remove(oclass=cobramm.intensity)
        #         #         if obj.get(oclass=cobramm.wavelength):
        #         #             obj.remove(oclass=cobramm.wavelength)
        #         #         for wav, ints in zip(*spc):
        #         #             # print("Adding {} {}".format(wav, ints))
        #         #             obj.add(cobramm.intensity(scalar_value=ints, unit="arbitrary"))
        #         #             obj.add(cobramm.wavelength(scalar_value=wav, unit="nm"))
                       
        #             #V6
        #             # if spc is not None:
        #             #     if obj.get(oclass=cobramm.spectrum_element):
        #             #         obj.remove(oclass=cobramm.spectrum_element)
        #             #     for wav, ints in zip(*spc):
        #             #         value = obj.add(cobramm.spectrum_value())
        #             #         value.add(cobramm.intensity(scalar_value=ints, unit="arbitrary"))
        #             #         value.add(cobramm.wavelength(scalar_value=wav, unit="nm"))

        #         yield obj
        #     else:
        #         yield None

    # OVERRIDE
    def _apply_added(self, root_obj, buffer):
        """Adds the added cuds to the engine."""
        # print("\nbuffer\n", buffer.values())
        for obj in buffer.values():
            # add a new atom
            if obj.is_a(cobramm.atom):
                if len(obj.get(oclass=cobramm.position)) != 0:
                    # print("\n\n-------------------------", obj.get(oclass=cobramm.position))
                    pos = obj.get(oclass=cobramm.position)[0].vector_value
                    sym = obj.get(oclass=cobramm.element_id)[0].atom_symbol
                    self._engine.add_atom(obj.uid, sym, pos)
            # specify temperature of the simulation
            elif obj.is_a(cobramm.temperature):
                self._engine.add_property("temperature", obj.scalar_value)
            # specify pressure of the simulation
            elif obj.is_a(cobramm.pressure):
                self._engine.add_property("pressure", obj.scalar_value)
            # give solvent commercial name
            elif obj.is_a(cobramm.molecule):
                self._engine.add_property("solvent", obj.commercial_name)
            elif obj.is_a(cobramm.accuracy):
                self._engine.add_property("accuracy", obj.accuracy_value)
            elif obj.is_a(cobramm.charge):
                self._engine.add_property("charge", obj.integer_value)
            elif obj.is_a(cobramm.case) and obj.case_name in ["transient", "emission", "absorption"]:
                self._engine.add_property("case", obj.case_name)

    # OVERRIDE
    def _apply_updated(self, root_obj, buffer):
        """Updates the updated cuds in the engine."""
        for obj in buffer.values():
            # check if position has been updated
            if obj.is_a(cobramm.position):
                atom = obj.get(rel=cobramm.is_part_of)[0]
                sym = atom.get(oclass=cobramm.element_id)[0].atom_symbol
                self._engine.add_atom(atom.uid, sym, obj.vector_value)
            # check if atomic symbol has been updated
            elif obj.is_a(cobramm.element_id):
                atom = obj.get(rel=cobramm.is_part_of)[0]
                pos = atom.get(oclass=cobramm.position)[0].vector_value
                self._engine.add_atom(atom.uid, obj.atom_symbol, pos)
            # check if temperature has changed
            elif obj.is_a(cobramm.temperature):
                self._engine.add_property("temperature", obj.scalar_value)
            # check if pressure has changed
            elif obj.is_a(cobramm.pressure):
                self._engine.add_property("pressure", obj.scalar_value)
            # check if solvent has changed
            elif obj.is_a(cobramm.solvent):
                self._engine.add_property("solvent", obj.commercial_name)

    # OVERRIDE
    def _apply_deleted(self, root_obj, buffer):
        """Deletes the deleted cuds from the engine."""
        for obj in buffer.values():
            if obj.is_a(cobramm.atom):
                self._engine.delete_atom(obj.uid)
            elif obj.is_a(cobramm.temperature):
                self._engine.delete_property("temperature")
            elif obj.is_a(cobramm.pressure):
                self._engine.delete_property("pressure")
            elif obj.is_a(cobramm.solvent):
                self._engine.delete_property("solvent")

    # OVERRIDE
    def _initialise(self, root_obj):
        """Initialise the session.
        This method is executed before the first run.

        Args:
            root_obj (Cuds): The wrapper cuds object.
        """
        pass
