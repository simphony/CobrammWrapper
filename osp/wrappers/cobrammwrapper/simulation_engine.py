import numpy as np


class CobrammSimulationEngine:
    """
    Simple engine sample code.
    """

    def __init__(self):
        self.executed = False
        print("Engine instantiated!")

    def __str__(self):
        return "Some Engine Connection"

    def run(self):
        """Call the run command of the engine."""
        print("Now the engine is running")
        self.executed = True

    def add_atom(self, uid, position, velocity):
        """"""
        print("Added atom %s with position %s and velocity %s"
              % (uid, position, velocity))

    def update_position(self, uid, position):
        """"""
        print("Update atom %s. Setting position to %s"
              % (uid, position))

    def update_velocity(self, uid, velocity):
        """"""
        print("Update atom %s. Setting velocity to %s"
              % (uid, velocity))

    def delete_atom(self, uid):
        """"""
        print("Deleted atom %s")

    def get_velocity(self, uid):
        if self.executed:
            return np.array([42, 42, 42])

    def get_position(self, uid):
        if self.executed:
            return np.array([-42, -42, -42])
