import numpy as np
import sys
import os
import copy

# add $COBRAM_PATH/cobramm directory to import COBRAMM modules
try:
    sys.path.append(os.path.join(os.environ['COBRAM_PATH'], 'cobramm'))
except KeyError:
    raise RuntimeError('cannot import cobramm module.\n'
                       'Please check if $COBRAM_PATH is defined in your environment, if not run:\n'
                       'export COBRAM_PATH=/path/to/COBRAMM')

from amberCalculator import AmberCalculator, AmberInput
from cobrammCalculator import CobrammCalculator, CobrammInput
import constants
from layers import Layers
from harmonicSampling import HarmonicSampling

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

        # initialize attributes to store the spectral results of the simulation
        self.grid, self.spectrum = [], []

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
        qm_basis = "sto-3g"
        gs_functional = "b3lyp"
        td_functional = "cam-b3lyp"
        qmmm_opt_maxsteps = 50
        sampling_type = "wigner"

        # define additional parameters based on the required accuracy
        if self.accuracy == 1:
            boxsize = 9.   # size of the MD simulation box, in Ang
            cutoff = 5.  # cutoff for potential evaluation in the MD simulations, in Ang
            solvradius = 7.  # size of the solvent droplet cut from the MD simulation box, in Ang
            mradius = 5.  # radius of the mobile layer in the QM/MM scheme, in Ang
            optsteps = 100  # number of initial optimization step in the MD equilibration
            timestep = 0.002  # time step of the MD run (in ps)
            heatingtime = 5.  # time (in ps) of the thermalization step
            equiltime = 25.  # time (in ps) of the final equilibration step
            nsamples = 10  # number of samples of the wigner sampling
        elif self.accuracy == 2:
            boxsize = 11.   # size of the MD simulation box, in Ang
            cutoff = 7.  # cutoff for potential evaluation in the MD simulations, in Ang
            solvradius = 9.  # size of the solvent droplet cut from the MD simulation box, in Ang
            mradius = 7.  # radius of the mobile layer in the QM/MM scheme, in Ang
            optsteps = 500  # number of initial optimization step in the MD equilibration
            timestep = 0.002  # time step of the MD run (in ps)
            heatingtime = 10.  # time (in ps) of the thermalization step
            equiltime = 50.  # time (in ps) of the final equilibration step
            nsamples = 20  # number of samples of the wigner sampling
        else:
            boxsize = 13.   # size of the MD simulation box, in Ang
            cutoff = 9.  # cutoff for potential evaluation in the MD simulations, in Ang
            solvradius = 11.  # size of the solvent droplet cut from the MD simulation box, in Ang
            mradius = 9.  # radius of the mobile layer in the QM/MM scheme, in Ang
            optsteps = 1000  # number of initial optimization step in the MD equilibration
            timestep = 0.002  # time step of the MD run (in ps)
            heatingtime = 20.  # time (in ps) of the thermalization step
            equiltime = 100.  # time (in ps) of the final equilibration step
            nsamples = 100  # number of samples of the wigner sampling

        # directories where to run AMBER calculations
        _CALC01_DIR = "01_opt"
        _CALC02_DIR = "02_heat"
        _CALC03_DIR = "03_equil"
        # output files with plot of the equilibration steps
        _CALC01_PLOT = "optimization.pdf"
        _CALC02_PLOT = "thermalization.pdf"
        _CALC03_PLOT = "equilibration.pdf"

        heatingtime = 1.  # time (in ps) of the thermalization step
        equiltime = 1.  # time (in ps) of the final equilibration step
        qmmm_opt_maxsteps = 5
        nsamples = 5

        print("\n * Create AMBER topology and coordinates for a molecule in a box of solvent")

        # create a random box of solvent with AMBER topology and parameters
        init_snap, topology = AmberCalculator.createSolvatedMolecule(atomlabels, atomcoords, solvsize=boxsize,
                                                                     solvent=solv)

        print("\n * Geometry optimization at the MM level, {} steps".format(optsteps))

        # define whether PBC are defined
        coord_has_pbc = init_snap.unitcell is not None

        # prepare the input for a minimization run
        input01 = AmberInput(minimize=True, nrSteps=optsteps, nprntsteps=1, cutoff=cutoff, usePBC=coord_has_pbc)

        # run the optimization, and do not remove directory at the end
        output01 = AmberCalculator.run(input01, topology, init_snap, calcDir=_CALC01_DIR, nCores=1,  store=True)

        # analysis of the optimization
        output01.optimizationPlot(_CALC01_PLOT)
        print("   Plot of the energy optimization written to file {0}".format(_CALC01_PLOT))

        # extract optimized geometry
        opt_snap = output01.snapshot(-1)

        print("\n * Heating the system at {0} K for {1} ps".format(temp, heatingtime))

        # prepare the input for the heating
        nr_steps = int(heatingtime / timestep)
        prt_steps = int(0.1 / timestep)  # print snapshot every 100 fs
        input02 = AmberInput(minimize=False, nrSteps=nr_steps, deltaT=timestep, temperature=temp,
                             pressure=None, cutoff=cutoff, freezeH=True, nprntsteps=prt_steps, usePBC=coord_has_pbc)

        # run the optimization, and do not remove directory at the end
        output02 = AmberCalculator.run(input02, topology, opt_snap, calcDir=_CALC02_DIR, nCores=1, store=True)

        # analysis of the optimization
        output02.dynamicsPlot(_CALC02_PLOT)
        print("   Plot of the thermalization dynamics written to file {0}".format(_CALC02_PLOT))

        # extract thermalized geometry
        thermal_snap = output02.snapshot(-1)

        print("\n * Equilibrating the system at {0} K and {1} bar for {2} ps".format(temp, prs, equiltime))

        # prepare the input for the heating
        nr_steps = int(equiltime / timestep)
        prt_steps = int(0.1 / timestep)  # print snapshot every 100 fs
        input03 = AmberInput(minimize=False, nrSteps=nr_steps, deltaT=timestep, temperature=temp, pressure=prs,
                             cutoff=cutoff, freezeH=True, nprntsteps=prt_steps, usevelocity=True, usePBC=coord_has_pbc)

        # run the equilibration MD with AMBER
        output03 = AmberCalculator.run(input03, topology, thermal_snap, calcDir=_CALC03_DIR, nCores=1, store=True)

        # analysis of the optimization
        output03.dynamicsPlot(_CALC03_PLOT)
        print("   Plot of the equilibration dynamics written to file {0}".format(_CALC03_PLOT))

        # extract final equilibration geometry and write it to file
        equil_snap = output02.snapshot(-1)

        print("\n * Extracting a droplet of solvent of radius {0} Ang\n"
              "   with mobile molecules within {1} Ang from the chromophore".format(solvradius, mradius))

        # use AmberCalculator to create the Amber snapshot and topology objects that represent
        # a molecule surrounded by a solvent droplet of radius give by the argument spherRad
        droplet_snap, droplet_topo = AmberCalculator.cutdroplet(topology, equil_snap, spherRad=solvradius)

        # now we can create input for COBRAMM with the qmmmLayers method of CobrammCalculator
        inp_geometry, modelh_topo, real_topo = CobrammCalculator.qmmmLayers(
            droplet_snap, droplet_topo, mobileThreshold=mradius, reorderRes=True)

        print("\n * Geometry optimization at the QM/MM level, max {} steps".format(qmmm_opt_maxsteps))

        # define cobram.command for a ground state optimization with DFT
        cobcomm = CobrammInput(calc_type="optxg", qm_basis=qm_basis, qm_functional=gs_functional,
                               nr_steps=qmmm_opt_maxsteps)

        # run the cobramm optimization and extract the optimized geometry
        optresult = CobrammCalculator.run(cobcomm, inp_geometry, modelh_topo, real_topo,
                                          calc_dir="optimization", store=True)
        opt_geometry = optresult.snapshot(-1)

        grid_e, spectrum_e = None, None
        if sampling_type == "vertical":

            print("\n * Computing electronic transitions for the optimized geometry\n"
                  "   and computing the electronic spectrum by gaussian convolution ")

            # # define cobram.command for a single point calculation, with TD-DFT and 5 electronic states
            cobcomm = CobrammInput(n_el_states=5, qm_basis=qm_basis, qm_functional=td_functional)

            # run the cobramm optimization and extract the optimized geometry
            tddftresult = CobrammCalculator.run(cobcomm, opt_geometry, modelh_topo, real_topo,
                                                calc_dir="tddft", store=True)

            # now process the electronic state results to get the spectrum as a function of the energy
            grid_e = np.linspace(0.1 / constants.Hartree2eV, 20. / constants.Hartree2eV, 1000)
            spectrum_e = tddftresult.eletronicspectrum(grid_e, width=0.1 / constants.Hartree2eV)

        elif sampling_type == "wigner":

            # normal mode computation for the wigner sampling
            print("\n * Computing the normal modes of the molecule at the optimized geometry")

            # redefine M layer as L layer, to freeze atoms in freq calculation
            standard_text = opt_geometry.reallayertext
            modified_text = ""
            for line in standard_text.splitlines():
                if line.split()[5] == "M":
                    modified_text += line.replace(" M", " L") + "\n"
                else:
                    modified_text += line + "\n"
            frozen_mlayer = Layers.from_real_layers_xyz(modified_text)

            # define cobram.command for a ground state frequency calculation
            cobcomm = CobrammInput(calc_type="freqxg", qm_basis=qm_basis, qm_functional=gs_functional)

            # run COBRAMM
            freqresults = CobrammCalculator.run(cobcomm, frozen_mlayer, modelh_topo, real_topo,
                                                calc_dir="frequencies", store=True)

            print("\n * Computing electronic transitions for the a set of {} sample geometries\n"
                  "   according to a Wigner sampling of the harmonic vibrational wavefunction".format(nsamples))

            # extract equilibrium geometry and format as a simple 1D vector of 3N elements
            geomvector = []
            for x, y, z in zip(*frozen_mlayer.model):
                geomvector.append(x), geomvector.append(y), geomvector.append(z)
            # define harmonic sampling of the molecular oscillations
            mol_oscillator = HarmonicSampling(geomvector, freqresults.coord_masses, freqresults.force_matrix, temp)

            # define directory where to store displaced geometries
            displ_dir = "sampling"
            # check the existence of the directory where to store samples, in case create the directory
            if not os.path.isdir(displ_dir):
                os.mkdir(displ_dir)

            # now randomly generate snapshots for subsequent calculations
            displ_geometries = []
            for i_geom in range(nsamples):
                # get sample of the wigner distribution
                newgeomvector = mol_oscillator.get_sample()
                # create a new Layers object to store the new snapshot, move the H layer and store the snapshot
                new_geometry = Layers.from_real_layers_xyz(opt_geometry.reallayertext)
                new_geometry.updateHlayer([newgeomvector[0::3], newgeomvector[1::3], newgeomvector[2::3]])
                displ_geometries.append(new_geometry)

            # define cobram.command for a single point calculation, with TD-DFT and 5 electronic states
            cobcomm = CobrammInput(n_el_states=5, qm_basis=qm_basis, qm_functional=td_functional)

            # now run the calculations
            displ_results = []
            for i, new_geometry in enumerate(displ_geometries):

                # define name of the directory where to store the sample point
                path_snap = os.path.join(displ_dir, "sample_{0:04d}".format(i))

                # run the cobramm optimization and extract the optimized geometry
                tddftresult = CobrammCalculator.run(cobcomm, new_geometry, modelh_topo, real_topo,
                                                    calc_dir=path_snap, store=True)
                displ_results.append(tddftresult)

            # compute grid for absorption spectrum
            grid_e = np.linspace(0.1 / constants.Hartree2eV, 20. / constants.Hartree2eV, 1000)
            # compute and accumulate spectrum
            spectrum_e = np.array([0.0] * len(grid_e))
            for tddftresult in displ_results:
                spectrum_e += tddftresult.eletronicspectrum(grid_e, width=0.1 / constants.Hartree2eV)
            # normalize spectrum by the number of sample TD-DFT calculations
            spectrum_e /= float(nsamples)

        # convert the spectrum to wavelength and store the results
        for en, ints_e in zip(grid_e, spectrum_e):
            wavelength = 10000000. / (en / constants.wavnr2au)
            self.grid.append(wavelength)
            self.spectrum.append(ints_e / wavelength**2)

        file_name = "spectrum_nm.dat"
        with open(file_name, "w") as f:
            for wav, ints in zip(self.grid, self.spectrum):
                f.write("{} {}\n".format(wav, ints))

        print("\nCOBRAMM run has been completed!")

        # COBRAMM has been executed
        self.executed = True

    def get_spectrum(self):
        if self.executed:
            return self.grid, self.spectrum
        else:
            return None
