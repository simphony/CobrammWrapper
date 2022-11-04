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
from cobrammCalculator import CobrammCalculator, CobrammInput, CobrammOutput
# from spectronCalculator import readCobrammOut, SpectronCalculator
import constants
from layers import Layers
from harmonicSampling import HarmonicSampling
# from transientCalculator import PumpProbeCalculator, MultiwfnCalculator
from transientCalculator import PumpProbeCalculator, MultiwfnCalculator, readCobrammOut, SpectronCalculator



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
        self.accuracy = 2
        # define user case: absorption or emission or transient
        self.case = "transient"
        # initialize attributes to store the spectral results of the simulation
        self.grid1, self.spectrum1 = [], []
        self.grid2, self.spectrum2 = [], []
        self.grid3, self.spectrum3 = [], []

        # initialize Amber wrapper class
        AmberCalculator()
        # initialize COBRAMM wrapper class
        CobrammCalculator()

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

        # set accuracy and case
        self.accuracy = self.input_setup["accuracy"]
        if "case" in self.input_setup:
            self.case = self.input_setup["case"]

        print("[DEBUG]: Accuracy = {}".format(self.accuracy))
        print("[DEBUG]: Charge = {}".format(self.input_setup["charge"]))
        print("[DEBUG]: Case = {}".format(self.case))

        # extract lists of atomic labels and coordinates
        atomlabels, atomcoords = [], []
        for lab, crd in self.atoms.values():
            atomlabels.append(lab)
            atomcoords.append(crd)
        # extract solvent name, temperature and pressure
        solv = self.input_setup["solvent"]
        temp = self.input_setup["temperature"]
        prs = self.input_setup["pressure"]
        chrg = self.input_setup["charge"]

        ##### PARAMETERS FOR USER CASE TRANSIENT #####

        ## DEFAULT PARAMETERS
        default_initial_state = "bright" #bright or int
        default_final_state = 0 #int
        default_decay = "slow"
        dafault_tau = 0.00001  # customized tau in case of decay = custom

        ## SIMULATION PARAMETERS
        initial_state = default_initial_state
        final_state = default_final_state
        decay = default_decay # slow,medium,fast custom.
        tau = dafault_tau
        adia_sim_time = 250
        adia_time_step = 10
        PP_sim_time = 250
        PP_time_step = 10
        tsh_simulation_time = 100
        ###############################################

        print("\nNow COBRAMM is running, computing the emission spectrum of the input molecule")
        print("Solvent: {}".format(solv))
        print("Temperature: {} K and Pressure: {} bar".format(temp, prs))
        
        # define common parameters that are required to run the simulation
        qm_basis_low = "3-21g"
        qm_basis_high = "6-31g*"
        gs_functional = "b3lyp"
        td_functional = "cam-b3lyp"
        pre_opt_maxsteps = 25     # number of steps of pre-optimisation in S1
        qmmm_opt1_maxsteps = 25   # number of steps optimisation steps for all the levels of accuracy
        qmmm_opt2_maxsteps = 25  # number of steps of higher level optimisation for level 2 and 3
        ## Davide: added new fixed parameters and removed the flexible parameters
        boxsize = 16
        cutoff = 8
        solvradius = 12
        mradius = 8
        optsteps = 1000
        timestep = 0.002
        heatingtime = 20
        if solv == "acetonitrile":
            equiltime = 200
        else:
            equiltime = 100
        nsamples = 50
        ntrj = 5

        #state to compute the gradient for in the 3rd user case in Gaussian notation (S1=1, S2=2 ecc)
        ##grad_to_compute = [1, 2]

        # directories where to run AMBER calculations
        _CALC01_DIR = "01_opt"
        _CALC02_DIR = "02_heat"
        _CALC03_DIR = "03_equil"
        # output files with plot of the equilibration steps
        _CALC01_PLOT = "optimization.pdf"
        _CALC02_PLOT = "thermalization.pdf"
        _CALC03_PLOT = "equilibration.pdf"


        ## run a pre-optimisation of S1 in gas phase with COBRAMM
        if self.case == "emission":
            print("\ncoordinate prima di pre-equilibration sono: {}".format(atomcoords)) ##da cancellare
            ## S1 pre optimisation ##
            print("\n * Create AMBER topology and coordinates for a molecule in gas phase and input files for COBRAMM")

            # create a random empty box of solvent with AMBER topology and parameters
            gas_snap, gas_topology = AmberCalculator.createSolvatedMolecule(atomlabels, atomcoords, solvsize=boxsize,
                                                                     dry=True, net_charge=chrg)

            # now we can create input for COBRAMM with the qmmmLayers method of CobrammCalculator
            gas_inp_geometry, gas_modelh_topo, gas_real_topo = CobrammCalculator.qmmmLayers(
                gas_snap, gas_topology, mobileThreshold=mradius, reorderRes=True, dry=True)

            # define cobram.command for an S1 state optimization with DFT
            cobcomm = CobrammInput(calc_type="optxg", qm_charge=chrg, qm_basis=qm_basis_low, qm_functional=td_functional,
                               nr_steps=pre_opt_maxsteps, n_root=1)

            # run the cobramm pre-optimization and extract the pre-optimized geometry
            preoptresult = CobrammCalculator.run(cobcomm, gas_inp_geometry, gas_modelh_topo, gas_real_topo,
                                          calc_dir="pre_optimization", store=True, n_cores=4)
            #opt_geometry_pre = preoptresult.snapshot(-1)

            ## reads the new charges and put in a list will be read later.
            # TODO: improving it...
            allQMcharges = preoptresult.qmcharges
            reverse=[]
            for natoms in range(len(atomlabels)):
                lastcharge=allQMcharges.pop()
                reverse.append(lastcharge)
            QMcharges=(reverse[::-1])

            # reads and stores the new coordinates
            pre_optcoords = preoptresult.coords

        print("\n * Create AMBER topology and coordinates for a molecule in a box of solvent")

        # create a random box of solvent with AMBER topology and parameters
        if self.case == "absorption" or self.case == "transient":
            init_snap, topology = AmberCalculator.createSolvatedMolecule(atomlabels, atomcoords, solvsize=boxsize,
                                                                     solvent=solv, net_charge=chrg)
        elif self.case == "emission":
            init_snap, topology = AmberCalculator.createSolvatedMolecule(atomlabels, pre_optcoords, solvsize=boxsize,
                                                                         solvent=solv, charges=True, newcharges=QMcharges, net_charge=chrg)


        print("\n * Geometry optimization at the MM level, {} steps".format(optsteps))

        # define whether PBC are defined
        coord_has_pbc = init_snap.unitcell is not None

        # prepare the input for a minimization run

        if self.case == "absorption" or self.case == "transient":
            input01 = AmberInput(minimize=True, nrSteps=optsteps, nprntsteps=1, cutoff=cutoff, usePBC=coord_has_pbc)
        elif self.case == "emission":
            input01 = AmberInput(minimize=True, nrSteps=optsteps, nprntsteps=1, cutoff=cutoff, freezesolute=True,
                                     usePBC=coord_has_pbc)

        # run the optimization, and do not remove directory at the end
        output01 = AmberCalculator.run(input01, topology, init_snap, calcDir=_CALC01_DIR, nCores=4,  store=True, ref=True)

        # analysis of the optimization
        output01.optimizationPlot(_CALC01_PLOT)
        print("   Plot of the energy optimization written to file {0}".format(_CALC01_PLOT))

        # extract optimized geometry
        opt_snap = output01.snapshot(-1)

        print("\n * Heating the system at {0} K for {1} ps".format(temp, heatingtime))

        # prepare the input for the heating
        nr_steps = int(heatingtime / timestep)
        prt_steps = int(0.1 / timestep)  # print snapshot every 100 fs

        if self.case == "absorption" or self.case == "transient":
            input02 = AmberInput(minimize=False, nrSteps=nr_steps, deltaT=timestep, temperature=temp,
                             pressure=None, cutoff=cutoff, freezeH=True, nprntsteps=prt_steps, usePBC=coord_has_pbc)
        elif self.case == "emission":
            input02 = AmberInput(minimize=False, nrSteps=nr_steps, deltaT=timestep, temperature=temp,
                             pressure=None, cutoff=cutoff, freezeH=True, freezesolute=True, nprntsteps=prt_steps, usePBC=coord_has_pbc)

        # run the heating, and do not remove directory at the end
        output02 = AmberCalculator.run(input02, topology, opt_snap, calcDir=_CALC02_DIR, nCores=4, store=True, ref=True)

        # analysis of the heating
        output02.dynamicsPlot(_CALC02_PLOT)
        print("   Plot of the thermalization dynamics written to file {0}".format(_CALC02_PLOT))

        # extract thermalized geometry
        thermal_snap = output02.snapshot(-1)

        print("\n * Equilibrating the system at {0} K and {1} bar for {2} ps".format(temp, prs, equiltime))

        # prepare the input for the equilibration
        nr_steps = int(equiltime / timestep)
        prt_steps = int(0.1 / timestep)  # print snapshot every 100 fs

        if self.case == "absorption" or self.case == "transient":
            input03 = AmberInput(minimize=False, nrSteps=nr_steps, deltaT=timestep, temperature=temp, pressure=prs,
                             cutoff=cutoff, freezeH=True, nprntsteps=prt_steps, usevelocity=True, usePBC=coord_has_pbc)
        elif self.case == "emission":
            input03 = AmberInput(minimize=False, nrSteps=nr_steps, deltaT=timestep, temperature=temp, pressure=prs,
                             cutoff=cutoff, freezeH=True, freezesolute=True, nprntsteps=prt_steps, usevelocity=True, usePBC=coord_has_pbc)

        # run the equilibration MD with AMBER
        output03 = AmberCalculator.run(input03, topology, thermal_snap, calcDir=_CALC03_DIR, nCores=4, store=True, ref=True)

        # analysis of the equilibration
        output03.dynamicsPlot(_CALC03_PLOT)
        print("   Plot of the equilibration dynamics written to file {0}".format(_CALC03_PLOT))

        # extract final equilibration geometry and write it to file
        equil_snap = output03.snapshot(-1)

        print("\n * Extracting a droplet of solvent of radius {0} Ang\n"
              "   with mobile molecules within {1} Ang from the chromophore".format(solvradius, mradius))

        # use AmberCalculator to create the Amber snapshot and topology objects that represent
        # a molecule surrounded by a solvent droplet of radius give by the argument spherRad
        droplet_snap, droplet_topo = AmberCalculator.cutdroplet(topology, equil_snap, spherRad=solvradius)

        # now we can create input for COBRAMM with the qmmmLayers method of CobrammCalculator
        inp_geometry, modelh_topo, real_topo = CobrammCalculator.qmmmLayers(
            droplet_snap, droplet_topo, mobileThreshold=mradius, reorderRes=True)

        print("\n * Geometry optimization at the QM/MM level, max {} steps".format(qmmm_opt1_maxsteps))

        # define cobram.command for a ground or excited state optimization with DFT
        if self.case == "absorption" or self.case == "transient":
            cobcomm = CobrammInput(calc_type="optxg", qm_charge=chrg, qm_basis=qm_basis_low, qm_functional=gs_functional,
                               nr_steps=qmmm_opt1_maxsteps)
        if self.case == "emission":
            cobcomm = CobrammInput(calc_type="optxg", qm_charge=chrg, qm_basis=qm_basis_low, qm_functional=td_functional,
                               nr_steps=qmmm_opt1_maxsteps, n_root=1)

        # run the cobramm optimization and extract the optimized geometry
        optresult = CobrammCalculator.run(cobcomm, inp_geometry, modelh_topo, real_topo,
                                          calc_dir="optimization_low", store=True, n_cores=4)
        opt_geometry_low = optresult.snapshot(-1)

        if self.case == "absorption" or self.case == "emission":

            print("\n * Computing electronic transitions for the optimized geometry\n"
                  "   and computing the electronic spectrum by gaussian convolution at level of accuracy = 1")

            # define cobram.command for a single point calculation, with TD-DFT and 5 electronic states
            if self.case == "absorption":
                cobcomm = CobrammInput(n_el_states=5, qm_charge=chrg, qm_basis=qm_basis_low, qm_functional=td_functional)
            elif self.case == "emission":
                cobcomm = CobrammInput(n_el_states=1, qm_charge=chrg, qm_basis=qm_basis_low, qm_functional=td_functional, n_root=1)

            # run the cobramm optimization and extract the optimized geometry
            tddftresult = CobrammCalculator.run(cobcomm, opt_geometry_low, modelh_topo, real_topo,
                                            calc_dir="tddft_level_1", store=True, n_cores=4)
            grid_e, spectrum_e = None, None
            # now process the electronic state results to get the spectrum as a function of the energy
            grid_e = np.linspace(0.1 / constants.Hartree2eV, 20. / constants.Hartree2eV, 1000)
            spectrum_e = tddftresult.eletronicspectrum(grid_e, width=0.1 / constants.Hartree2eV)

            for en, ints_e in zip(grid_e, spectrum_e):
                wavelength = 10000000. / (en / constants.wavnr2au)
                self.grid1.append(wavelength)
                self.spectrum1.append(ints_e / wavelength**2)

            file_name = "spectrum_nm_level_1.dat"
            with open(file_name, "w") as f:
                for wav, ints in zip(self.grid1, self.spectrum1):
                    f.write("{} {}\n".format(wav, ints))

        ## Davide: run an additional optimisation for level 2 or 3

            if self.accuracy == 2 or self.accuracy == 3:
                print("\n * Additional geometry QM/MM optimization for level of accuracy {0}, max {1} steps".format(
                    self.accuracy,qmmm_opt2_maxsteps))

                if self.case == "absorption":
                # define cobram.command for an additional ground state optimization with DFT
                    cobcomm = CobrammInput(calc_type="optxg", qm_charge=chrg, qm_basis=qm_basis_high, qm_functional=gs_functional,
                                   nr_steps=qmmm_opt2_maxsteps)
                elif self.case == "emission":
                    cobcomm = CobrammInput(calc_type="optxg", qm_charge=chrg, qm_basis=qm_basis_high, qm_functional=td_functional,
                               nr_steps=qmmm_opt2_maxsteps, n_root=1)

                # run the cobramm optimization and extract the optimized geometry
                optresult = CobrammCalculator.run(cobcomm, opt_geometry_low, modelh_topo, real_topo,
                                              calc_dir="optimization_high", store=True, n_cores=4)  # ## check directory,
                # now it should simply overwrite the dir "optimization"
                opt_geometry_high = optresult.snapshot(-1)
                grid_e, spectrum_e = None, None

                print("\n * Computing electronic transitions for the optimized geometry\n"
                      "   and computing the electronic spectrum by gaussian convolution at level of accuracy = 2")

                # define cobram.command for a single point calculation, with TD-DFT and 5 or 1 electronic states
                if self.case == "absorption":
                    cobcomm = CobrammInput(n_el_states=5, qm_charge=chrg, qm_basis=qm_basis_high, qm_functional=td_functional)
                elif self.case == "emission":
                    cobcomm = CobrammInput(n_el_states=1, qm_charge=chrg, qm_basis=qm_basis_high, qm_functional=td_functional, n_root=1)
                    # run the cobramm optimization and extract the optimized geometry

                tddftresult = CobrammCalculator.run(cobcomm, opt_geometry_high, modelh_topo, real_topo,
                                                calc_dir="tddft_level_2", store=True, n_cores=4)

                # now process the electronic state results to get the spectrum as a function of the energy
                grid_e = np.linspace(0.1 / constants.Hartree2eV, 20. / constants.Hartree2eV, 1000)
                spectrum_e = tddftresult.eletronicspectrum(grid_e, width=0.1 / constants.Hartree2eV)

                for en, ints_e in zip(grid_e, spectrum_e):
                    wavelength = 10000000. / (en / constants.wavnr2au)
                    self.grid2.append(wavelength)
                    self.spectrum2.append(ints_e / wavelength**2)

                file_name = "spectrum_nm_level_2.dat"
                with open(file_name, "w") as f:
                    for wav, ints in zip(self.grid2, self.spectrum2):
                        f.write("{} {}\n".format(wav, ints))

                if self.accuracy == 3:

                # normal mode computation for the wigner sampling
                    print("\n * Computing the normal modes of the molecule at the optimized geometry")

                # redefine M layer as L layer, to freeze atoms in freq calculation
                    standard_text = opt_geometry_high.reallayertext
                    modified_text = ""
                    for line in standard_text.splitlines():
                        if line.split()[5] == "M":
                            modified_text += line.replace(" M", " L") + "\n"
                        else:
                            modified_text += line + "\n"
                    frozen_mlayer = Layers.from_real_layers_xyz(modified_text)

                    if self.case == "absorption":
                    # define cobram.command for a ground state frequency calculation
                        cobcomm = CobrammInput(calc_type="freqxg", qm_charge=chrg, qm_basis=qm_basis_high, qm_functional=gs_functional)
                    elif self.case == "emission":
                        cobcomm = CobrammInput(calc_type="freqxg", qm_charge=chrg, qm_basis=qm_basis_high, qm_functional=td_functional,n_root=1)

                # run COBRAMM
                    freqresults = CobrammCalculator.run(cobcomm, frozen_mlayer, modelh_topo, real_topo,
                                                    calc_dir="frequencies", store=True, n_cores=4)

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
                        new_geometry = Layers.from_real_layers_xyz(opt_geometry_high.reallayertext)
                        new_geometry.updateHlayer([newgeomvector[0::3], newgeomvector[1::3], newgeomvector[2::3]])
                        displ_geometries.append(new_geometry)

                # define cobram.command for a single point calculation, with TD-DFT and 5 electronic states
                    if self.case == "absorption":
                        cobcomm = CobrammInput(n_el_states=5, qm_charge=chrg, qm_basis=qm_basis_high, qm_functional=td_functional)
                    elif self.case == "emission":
                        cobcomm = CobrammInput(n_el_states=1, qm_charge=chrg, qm_basis=qm_basis_high, qm_functional=td_functional, n_root=1)
                # now run the calculations
                    displ_results = []
                    for i, new_geometry in enumerate(displ_geometries):

                    # define name of the directory where to store the sample point
                        path_snap = os.path.join(displ_dir, "sample_{0:04d}".format(i))

                    # run the cobramm optimization and extract the optimized geometry
                        tddftresult = CobrammCalculator.run(cobcomm, new_geometry, modelh_topo, real_topo,
                                                    calc_dir=path_snap, store=True, n_cores=4)
                        displ_results.append(tddftresult)

                    grid_e, spectrum_e = None, None
                # compute grid for absorption spectrum
                    grid_e = np.linspace(0.1 / constants.Hartree2eV, 20. / constants.Hartree2eV, 1000)
                # compute and accumulate spectrum
                    spectrum_e = np.zeros(len(grid_e))
                    for tddftresult in displ_results:
                        spectrum_e += tddftresult.eletronicspectrum(grid_e, width=0.2 / constants.Hartree2eV)
                # normalize spectrum by the number of sample TD-DFT calculations
                    spectrum_e /= float(nsamples)

            # convert the spectrum to wavelength and store the results
                    for en, ints_e in zip(grid_e, spectrum_e):
                        wavelength = 10000000. / (en / constants.wavnr2au)
                        self.grid3.append(wavelength)
                        self.spectrum3.append(ints_e / wavelength**2)

                    file_name = "spectrum_nm_level_3.dat"
                    with open(file_name, "w") as f:
                        for wav, ints in zip(self.grid3, self.spectrum3):
                            f.write("{} {}\n".format(wav, ints))

        elif self.case == "transient":

                #run additional optimisation
            print("\n * Additional geometry QM/MM optimization, max {1} steps".format(
                    self.accuracy,qmmm_opt2_maxsteps))

            # define cobram.command for an additional ground state optimization with DFT
            cobcomm = CobrammInput(calc_type="optxg", qm_charge=chrg, qm_basis=qm_basis_high, qm_functional=gs_functional,
                                nr_steps=qmmm_opt2_maxsteps)

            # run the cobramm optimization and extract the optimized geometry
            optresult = CobrammCalculator.run(cobcomm, opt_geometry_low, modelh_topo, real_topo,
                                           calc_dir="optimization_high", store=True, n_cores=4)  # ## check directory,

            opt_geometry_high = optresult.snapshot(-1)

            if self.accuracy == 1 or self.accuracy == 2:

                #compute excited states calculation to get energies and dipoles
                print("\n * Computing electronic transitions for the optimized geometry\n")

            # define cobram.command for a single point calculation, with TD-DFT and 5 states
                cobcomm = CobrammInput(n_el_states=5, qm_charge=chrg, qm_basis=qm_basis_high, qm_functional=td_functional, tda=True)

                CobrammCalculator.run(cobcomm, opt_geometry_high, modelh_topo, real_topo,
                                                            calc_dir="tddft", store=True, n_cores=4)

                readCobrammOut(energies=True, en_dir="tddft", es_dipoles=True, tda_dir="tddft").write_files()

                #compute the ground state frequencies
                print("\n * Computing the normal modes of the molecule at the optimized geometry")

                # redefine M layer as L layer, to freeze atoms in freq calculation
                standard_text = opt_geometry_high.reallayertext
                modified_text = ""
                for line in standard_text.splitlines():
                    if line.split()[5] == "M":
                        modified_text += line.replace(" M", " L") + "\n"
                    else:
                        modified_text += line + "\n"
                frozen_mlayer = Layers.from_real_layers_xyz(modified_text)

                # define cobram.command for a ground state frequency calculation
                cobcomm = CobrammInput(calc_type="freqxg", qm_charge=chrg, qm_basis=qm_basis_high, qm_functional=gs_functional)

           # run COBRAMM
                CobrammCalculator.run(cobcomm, frozen_mlayer, modelh_topo, real_topo,
                                                   calc_dir="frequencies", store=True, n_cores=4)
                freqfile = (os.path.join("frequencies", "geometry.log"))


                #compute the gradients for n states
                print("\n * Computing the the gradient for the excited states of interest")

                gradient_files_string = ""
            #gradient_files_string = "gradient_S1.txt gradient_S2.txt gradient_S3.txt gradient_S4.txt gradient_S5.txt "
                for singlet in range(1, 6):
                    cobcomm = CobrammInput(calc_type="optxg", qm_charge=chrg, qm_basis=qm_basis_high, qm_functional=td_functional,
                               nr_steps=1, n_root=singlet)
                    CobrammCalculator.run(cobcomm, opt_geometry_high, modelh_topo, real_topo,
                                              calc_dir="grad_singlet_{}".format(singlet), store=True, n_cores=4)

                    readCobrammOut(gradient=True, grad_dir="grad_singlet_{}".format(singlet)).write_files()

                    gradient_files_string += "gradient_S{}.txt ".format(singlet)

            
            # LINEAR ABSORPTION

            #run iSPECTRON
                print("\n * Preparing input files for SPECTRON through iSPECTRON")

                ispectronresult_LA=SpectronCalculator.runiSpectronA(grad_files=gradient_files_string, modes=freqfile, free=6, signal = "LA", tag="iSp")

                # run spectron
                print("\n * Computing Linear Absorption")
                SpectronCalculator.runSpectron(ispectronresult_LA)

            # ADIABATIC PP

            #run iSPECTRON
                print("\n * Computing Time resolved signals")

                ispectronresult_PP1 = SpectronCalculator.runiSpectronA(grad_files=gradient_files_string, modes=freqfile, free=6, signal="PP", tag="iSp", diagrams="ESA GSB SE", t2=0)

                spectronOut = SpectronCalculator.run_pp(inpDir=ispectronresult_PP1, outDir="2D_data_level_1", simulation_time=adia_sim_time, time_step=adia_time_step)

                SpectronCalculator.runiSpectronB(spectronOut)
            
                if self.accuracy == 2:
                
                    if initial_state == "bright":
                        initial_state = readCobrammOut.get_bright_state("tddft/QM_data/qmALL.log")

                #ispectronresult_PP2=SpectronCalculator.runiSpectronA(grad_files=gradient_files_string, modes=freqfile, free=6, signal="PP", tag="_level_2", diagrams="SE", t2=0)

                    SpectronCalculator.insert_monoexp_k(ispectronresult_PP1, decay_time=decay, tau=tau, nstates=6, active_state=initial_state, final_state=final_state)
                
                    spectronOut = SpectronCalculator.run_pp(inpDir=ispectronresult_PP1, outDir="2D_data_level_2", simulation_time=PP_sim_time, time_step=PP_time_step)

                    SpectronCalculator.runiSpectronB(spectronOut)


            if self.accuracy == 3:

            # normal mode computation for the wigner sampling
                print("\n * Computing the normal modes of the molecule at the optimized geometry")

            # redefine M layer as L layer, to freeze atoms in freq calculation
                standard_text = opt_geometry_high.reallayertext
                modified_text = ""
                for line in standard_text.splitlines():
                    if line.split()[5] == "M":
                        modified_text += line.replace(" M", " L") + "\n"
                    else:
                        modified_text += line + "\n"
                frozen_mlayer = Layers.from_real_layers_xyz(modified_text)

            # run COBRAMM
                cobcomm = CobrammInput(calc_type="freqxg", qm_charge=chrg, qm_basis=qm_basis_high, qm_functional=gs_functional)
                freqresults = CobrammCalculator.run(cobcomm, frozen_mlayer, modelh_topo, real_topo,
                                                    calc_dir="frequencies", store=True, n_cores=4)

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
                for i_geom in range(ntrj):
                # get sample of the wigner distribution
                    newgeomvector = mol_oscillator.get_sample()
                    # create a new Layers object to store the new snapshot, move the H layer and store the snapshot
                    new_geometry = Layers.from_real_layers_xyz(opt_geometry_high.reallayertext)
                    new_geometry.updateHlayer([newgeomvector[0::3], newgeomvector[1::3], newgeomvector[2::3]])
                    displ_geometries.append(new_geometry)

                # now run the calculations
                displ_results = []
                for i, new_geometry in enumerate(displ_geometries):

                    # define name of the directory where to store the sample point
                    path_snap = os.path.join(displ_dir, "sample_{0:04d}".format(i))
                    top_dir=os.getcwd()
                    if not os.path.isdir(path_snap):
                        os.mkdir(path_snap)
                    os.chdir(path_snap)

                    # vertical excitation and tsh

                    cobcommVE = CobrammInput(n_el_states=5, qm_charge=chrg, qm_basis=qm_basis_high,
                                             qm_functional=td_functional, tda=True)

                    CobrammCalculator.run(cobcommVE, opt_geometry_high, modelh_topo, real_topo,
                                          calc_dir="VE", store=True, n_cores=4)

                    readCobrammOut(energies=True, en_dir="VE", es_dipoles=True, tda_dir="VE").write_files()

                    initial_state = readCobrammOut.get_bright_state("VE/QM_data/qmALL.log")

                    # define cobram.command for a tsh
                    cobcommTSH = CobrammInput(calc_type="tsh", n_el_states=5, qm_charge=chrg, qm_basis=qm_basis_high,
                                              qm_functional=td_functional, nr_steps=200, n_root=initial_state)

                    CobrammCalculator.run(cobcommTSH, opt_geometry_high, modelh_topo, real_topo,
                                          calc_dir="TSH", store=True, n_cores=4)

                    os.chdir(top_dir)

                trj_folders = []
                for trj in range(ntrj):
                    trj_folders.append("sample_{0:04d}/TSH".format(trj))

                pp = PumpProbeCalculator(trjdir_list=trj_folders, simulation_time=100, delta_t=10, nstates=20)
                upper_dir = os.getcwd()
                print("\nSetting up vertical excitations with a time step of {} fs\n".format(pp.delta_t))
                os.chdir("sampling")
                sampling_dir = os.getcwd()
                # setup and run single point calculations
                for folder in trj_folders:

                    print("\nEntering in the trajectory folder : {}".format(folder)) #sample_{0:04d}\n".format(folder))

                    os.chdir("{}".format(folder))
                    pp.setup_vertical_excitations(ncores=4, basis_set=qm_basis_high)

                    os.chdir(sampling_dir)

                os.chdir(sampling_dir)

                pp = PumpProbeCalculator(simulation_time=100, delta_t=10, nstates=20,
                                 trjdir_list=trj_folders, en_min=0.01, en_max=10.0,
                                 en_width=0.3, t_width=0.3)

                print("Convolution of the spectra for each time\n")

                all_init_tdm = pp.extract_init_tdm()

                for step in range(0, pp.nsteps + 1, pp.delta):
                    t = step * pp.time_step

                    print("Convolution of the spectrum for time {}\n".format(t))

                    collected_values = pp.collect_values_single_time(nstep=step, all_init_tdm=all_init_tdm)
                    pp.get_spectrum_single_time_all_traj(values=collected_values, currtime=t)

                os.chdir(upper_dir)

                print("Convulationg the final pump-probe spectrum")

                pp.time_convolution(output="spectrumPP.txt")


        print("\nCOBRAMM run has been completed!")

        # COBRAMM has been executed
        self.executed = True

    def get_spectrum(self):
        if self.executed:
            return self.grid1, self.spectrum1, self.grid2, self.spectrum2, self.grid3, self.spectrum3
        else:
            return None
