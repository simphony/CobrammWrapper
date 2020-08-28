# noinspection PyUnresolvedReferences
from osp.core import cobramm
from osp.core.session import SimWrapperSession
from osp.wrappers.cobrammwrapper.simulation_engine import CobrammSimulationEngine


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

    # OVERRIDE
    def _load_from_backend(self, uids, expired=None):
        """Loads the cuds object from the simulation engine"""
        for uid in uids:
            if uid in self._registry:
                obj = self._registry.get(uid)

                if obj.is_a(cobramm.material):
                    spc = self._engine.get_spectrum()
                    if spc is not None:
                        obj.add(cobramm.spectrum())

                if obj.is_a(cobramm.spectrum):
                    spc = self._engine.get_spectrum()
                    if spc is not None:
                        # first delete previous spectrum possibly stored
                        if obj.get(oclass=cobramm.absorption_value):
                            obj.remove(oclass=cobramm.absorption_value)
                        # and add the new values
                        for wav, ints in zip(*spc):
                            value = obj.add(cobramm.absorption_value())
                            value.add(cobramm.intensity(scalar_value=ints, unit="arbitrary"))
                            value.add(cobramm.wavelength(scalar_value=wav, unit="nm"))

                yield obj
            else:
                yield None

    # OVERRIDE
    def _apply_added(self, root_obj, buffer):
        """Adds the added cuds to the engine."""
        for obj in buffer.values():
            # add a new atom
            if obj.is_a(cobramm.atom):
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
            elif obj.is_a(cobramm.solvent):
                self._engine.add_property("solvent", obj.commercial_name)

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
