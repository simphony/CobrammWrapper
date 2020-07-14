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

                # check if position has been updated
                if obj.is_a(cobramm.position):
                    atom = obj.get(rel=cobramm.is_part_of)[0]
                    pos = self._engine.get_position(atom.uid)
                if pos is not None:
                    obj.value = pos

            yield obj
        else:
            yield None

    # OVERRIDE
    def _apply_added(self, root_obj, buffer):
        """Adds the added cuds to the engine."""
        for obj in buffer.values():
            if obj.is_a(cobramm.atom):
                pos = obj.get(oclass=cobramm.position)[0].value
                self._engine.add_atom(obj.uid, pos)

    # OVERRIDE
    def _apply_updated(self, root_obj, buffer):
        """Updates the updated cuds in the engine."""
        for obj in buffer.values():
            # check if position has been updated
            if obj.is_a(cobramm.position):
                atom = obj.get(rel=cobramm.is_part_of)[0]
                self._engine.update_position(atom.uid, obj.value)

    # OVERRIDE
    def _apply_deleted(self, root_obj, buffer):
        """Deletes the deleted cuds from the engine."""
        for cuds_object in buffer.values():
            if cuds_object.is_a(cobramm.atom):
                self._engine.delete_atom(cuds_object.uid)

    # OVERRIDE
    def _initialise(self, root_obj):
        """Initialise the session.
        This method is executed before the first run.

        Args:
            root_obj (Cuds): The wrapper cuds object.
        """
        # TODO: Is there any initialisation required in the engine
        #  before the run method is called?
        # This method is optional, it is not always necessary
