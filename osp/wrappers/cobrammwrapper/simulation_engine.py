import numpy as np
import sys
import os

# add $COBRAM_PATH/cobramm directory to import COBRAMM modules
try:
    sys.path.append(os.path.join(os.environ['COBRAM_PATH'], 'cobramm'))
except KeyError:
    raise RuntimeError('cannot import cobramm module.\n'
                       'Please check if $COBRAM_PATH is defined in your environment, if not run:\n'
                       'export COBRAM_PATH=/path/to/COBRAMM')

from amberCalculator import AmberCalculator, AmberInput
# from cobrammCalculator import CobrammCalculator, CobrammInput


class CobrammSimulationEngine:
    """
    Simple engine sample code.
    """

    def __init__(self):

        # variables that tells whether the engine has been executed
        self.executed = False
        # dictionary to store the definition of the atoms
        self.atoms = {}
        # dictionary to store the other options of the simulation
        self.input_setup = {}
        # define integer values to set the requested accuracy of the calculation (1, 2 or 3, increasing accuracy)
        self.accuracy = 1

        # initialize Amber wrapper class
        AmberCalculator()

        print("COBRAMM Engine instantiated!")

    def __str__(self):
        return "COBRAMM Engine Connection"

    def add_atom(self, uid, sym, position):
        """Add or update atom with a given unique identifier,
        specifying atomic symbol and position (cartesian coordinates)"""
        self.atoms[uid] = [sym, position]

    def delete_atom(self, uid):
        """Delete an atom given its unique identifier"""
        del self.atoms[uid]

    def add_property(self, label, value):
        """Add or update an entry to the input_setup dictionary"""
        self.input_setup[label] = value

    def delete_property(self, label):
        """Delete an entry of input_setup using its unique label"""
        del self.input_setup[label]

    def run(self):
        """Call the run command of the engine."""
        print("Now COBRAMM is running")

        # TODO: add here the calls to COBRAMM to execute the
        #  spectroscopy simulation workflow

        # extract lists of atomic labels and coordinates
        atomlabels, atomcoords = [], []
        for lab, crd in self.atoms.values():
            atomlabels.append(lab)
            atomcoords.append(crd)
        # extract solvent name, temperature and pressure
        solv = self.input_setup["solvent"]
        temp = self.input_setup["temperature"]
        prs = self.input_setup["pressure"]

        # define common parameters that are required to run the simulation
        # gs_functional = "b3lyp"
        # td_functional = "cam-b3lyp"

        # define additional parameters based on the required accuracy
        if self.accuracy == 1:
            boxsize = 9.   # size of the MD simulation box, in Ang
            cutoff = 5.  # cutoff for potential evaluation in the MD simulations
            optsteps = 100  # number of initial optimization step in the MD equilibration
            timestep = 0.002  # time step of the MD run (in ps)
            heatingtime = 5.  # time (in ps) of the thermalization step
            equiltime = 25.  # time (in ps) of the final equilibration step
        elif self.accuracy == 2:
            boxsize = 11.   # size of the MD simulation box, in Ang
            cutoff = 7.  # cutoff for potential evaluation in the MD simulations
            optsteps = 500  # number of initial optimization step in the MD equilibration
            timestep = 0.002  # time step of the MD run (in ps)
            heatingtime = 10.  # time (in ps) of the thermalization step
            equiltime = 50.  # time (in ps) of the final equilibration step
        else:
            boxsize = 13.   # size of the MD simulation box, in Ang
            cutoff = 9.  # cutoff for potential evaluation in the MD simulations
            optsteps = 1000  # number of initial optimization step in the MD equilibration
            timestep = 0.002  # time step of the MD run (in ps)
            heatingtime = 20.  # time (in ps) of the thermalization step
            equiltime = 100.  # time (in ps) of the final equilibration step

        # directories where to run AMBER calculations
        _CALC01_DIR = "01_opt"
        _CALC02_DIR = "02_heat"
        _CALC03_DIR = "03_equil"
        # output files with plot of the equilibration steps
        _CALC01_PLOT = "optimization.pdf"
        _CALC02_PLOT = "thermalization.pdf"
        _CALC03_PLOT = "equilibration.pdf"
        
        print("\nCreate AMBER topology and coordinates for a molecule in a box of solvent\n")

        # create a random box of solvent with AMBER topology and parameters
        init_snap, topology = AmberCalculator.createSolvatedMolecule(atomlabels, atomcoords, solvsize=boxsize,
                                                                     solvent=solv)

        print("\nGeometry optimization at the MM level, {} steps\n".format(optsteps))

        # define whether PBC are defined
        coord_has_pbc = init_snap.unitcell is not None

        # prepare the input for a minimization run
        input01 = AmberInput(minimize=True, nrSteps=optsteps, nprntsteps=1, cutoff=cutoff, usePBC=coord_has_pbc)

        # run the optimization, and do not remove directory at the end
        output01 = AmberCalculator.run(input01, topology, init_snap, calcDir=_CALC01_DIR, nCores=1,  store=True)

        # analysis of the optimization
        output01.optimizationPlot(_CALC01_PLOT)
        print("Plot of the energy optimization written to file {0}\n\n".format(_CALC01_PLOT))

        # extract optimized geometry
        opt_snap = output01.snapshot(-1)

        print("Heating the system at {0} K for {1} ps".format(temp, heatingtime))

        # prepare the input for the heating
        nr_steps = int(heatingtime / timestep)
        prt_steps = int(0.1 / timestep)  # print snapshot every 100 fs
        input02 = AmberInput(minimize=False, nrSteps=nr_steps, deltaT=timestep, temperature=temp,
                             pressure=None, cutoff=cutoff, freezeH=True, nprntsteps=prt_steps, usePBC=coord_has_pbc)

        # run the optimization, and do not remove directory at the end
        output02 = AmberCalculator.run(input02, topology, opt_snap, calcDir=_CALC02_DIR, nCores=1, store=True)

        # analysis of the optimization
        output02.dynamicsPlot(_CALC02_PLOT)
        print("Plot of the thermalization dynamics written to file {0}\n\n".format(_CALC02_PLOT))

        # extract thermalized geometry
        thermal_snap = output02.snapshot(-1)

        print("Equilibrating the system at {0} K and {1} bar for {2} ps\n".format(temp, prs, equiltime))
        print("files are stored in the directory {0}".format(_CALC03_DIR))

        # prepare the input for the heating
        nr_steps = int(equiltime / timestep)
        prt_steps = int(0.1 / timestep)  # print snapshot every 100 fs
        input03 = AmberInput(minimize=False, nrSteps=nr_steps, deltaT=timestep, temperature=temp, pressure=prs,
                             cutoff=cutoff, freezeH=True, nprntsteps=prt_steps, usevelocity=True, usePBC=coord_has_pbc)

        # run the equilibration MD with AMBER
        output03 = AmberCalculator.run(input03, topology, thermal_snap, calcDir=_CALC03_DIR, nCores=1, store=True)

        # analysis of the optimization
        output03.dynamicsPlot(_CALC03_PLOT)
        print("Plot of the equilibration dynamics written to file {0}\n\n".format(_CALC03_PLOT))

        # extract final equilibration geometry and write it to file
        equil_snap = output02.snapshot(-1)

        # # use AmberCalculator to create the Amber snapshot and topology objects that represent
        # # a molecule surrounded by a solvent droplet of radius give by the argument spherRad
        # droplet_snap, droplet_topo = AmberCalculator.cutdroplet(topology, equil_snap, spherRad=solvradius)
        #
        # # now we can create input for COBRAMM with the qmmmLayers method of CobrammCalculator
        # geometry, modelh_topo, real_topo = CobrammCalculator.qmmmLayers(
        #     droplet_snap, droplet_topo, mobileThreshold=mradius, reorderRes=True)
        #
        # # define cobram.command for a ground state optimization with DFT
        # cobcomm = CobrammInput(calc_type="optxg", qm_basis=qmbasis, qm_functional=qmfunctional, nr_steps=7)
        #
        # # run the cobramm optimization and extract the optimized geometry
        # optresult = CobrammCalculator.run(cobcomm, geometry, modelh_topo, real_topo,
        #                                   calc_dir="optimization", store=True)
        # opt_geometry = optresult.snapshot(-1)
        #
        # # # define cobram.command for a single point calculation, with TD-DFT and 5 electronic states
        # cobcomm = CobrammInput(n_el_states=5, qm_basis=qmbasis, qm_functional=qmfunctional)
        #
        # # run the cobramm optimization and extract the optimized geometry
        # tddftresult = CobrammCalculator.run(cobcomm, opt_geometry, modelh_topo, real_topo,
        #                                     calc_dir="tddft", store=True)
        #
        # # now process the electronic state results to print the spectrum
        # grid = np.linspace(0., 0.5, 200)
        # spect = tddftresult.eletronicspectrum(grid)

        # COBRAMM has been executed
        self.executed = True

    def get_position(self, uid):
        if self.executed:
            return np.array([-42, -42, -42])
        print("Get position of atom {0}".format(uid))
