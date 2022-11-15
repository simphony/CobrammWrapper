# CobrammWrapper
SimPhoNy wrapper for the COBRAMM engine developed at UNIBO.

*Contact:* marco.garavelli@unibo.it (UNIBO)

## Requirements
- Install [COBRAMM version 2.3](https://gitlab.com/cobrammgroup/cobramm/-/wikis/Installation-Guide) and required software - Amber and Gaussian.  
- Add the directories of the required tools in the PATH environmental variable


## Compatibility
The following table describes the version compatibility between the [OSP core](https://github.com/simphony/osp-core) package and documentation presented in this project.

| __OSP core__ |
|:------------:|
|     3.4.0    |

The releases of OSP core are available [here](https://github.com/simphony/osp-core/releases).


## Description
COBRAMM wrapper for modelling the electronic (UV-VIS) linear absorption (LA) spectrum of a solvated neutral molecule in its singlet ground state. 
The [QM/MM COBRAMM engine](https://site.unibo.it/cobramm/en) developed at UNIBO is employed. COBRAMM performs hybrid QM/MM calculations by linking the Gaussian (QM) and Amber (MM) software packages. 

<img src="https://site.unibo.it/cobramm/en/@@images/97dcd2e2-e470-4cf4-a0e7-4d15b2e44018.png" width="200
" >

An explicit atomistic MM model (Amber force field) is used for the solvent, while an explicit atomistic QM model (DFT/TDDFT) is used for the solvated chromophore. A classical equilibration of the solvent-solute system followed by a QM/MM optimization with Wigner sampling and QM/MM excited state calculations is employed to generate the LA spectrum and its bandshape. Vibrational progressions are not considered at this stage. The wrapper makes available to the user three different schemes of increasing complexity and computational time (Low, Medium, High) that do set parameters for MD simulations (boxsize, cutoff, initial optimization steps, time step of the MD run, heating time, time of the final equilibration steps) and QM/MM simulations (solvent droplet radius, mobile layer radius, number of samples produced in the Wigner sampling) of increasing complexity, accuracy and computational cost. See the wrapper readme file for further info and the workflow details.

## Workflow
Molecular structure, solvent, temperature, pressure are provided as input by the user through the wrapper GUI (while the accuracy level is currently hard-coded in the wrapper). This info is employed to run the following workflow:

1) Solvent-solute MM optimization, thermalization and equilibration dynamics (Amber force field for the solvent and GAFF for the solute) with the charges estimated by the AM1-BCC method from antechamber.

2) Set-up of the high-medium-low layers for QM/MM COBRAMM driven computations from the last MD snapshot produced in 1) (which is supposed to return the better solvent organization): the solvent MM shell is selected as a spherical droplet around the chromophore,

3) QM(DFT: b3lyp, STO-3G)/MM optimization for locating the representative ground state equilibrium structure of the solvated chromophore,

4) QM(DFT: b3lyp, STO-3G)/MM frequency calculation on the QM/MM minimum optimized in 3),

5) Wigner ensemble generation,

6) Single point QM(TDDFT: cam-b3lyp, STO-3G, nroot=5)/MM calculations on top of each equilibrated Wigner snapshots found in 7),

7) Sum over gaussian-functions centered at the excitation energy values found in 6) (with the gaussian maximum set at the oscillator strength value and a fixed 0.2 ev standard deviation) divided by the number of Wigner samples.

The Low/Medium/High accuracy level do set parameters for MD simulations (boxsize, cutoff, initial optimization steps, time step of the MD run, heating time, time of the final equilibration steps) and QM/MM simulations (solvent droplet radius, mobile layer radius, number of samples produced in the Wigner sampling) of increasing complexity, accuracy and computational cost.

## The actual wrapper implementation
The code that is necessary to carry out the linear absorption simulation workflow is included in the main COBRAMM code (COBRAMM engine), that is also available independently of the wrapper using some utility scripts.
The wrapper has been developed according to the indications outlined by the SymPhoNy developers.
The information that is supposed to come from the SimDOME GUI (molecule, solvent, temperature and pressure, level of accuracy of the calculation) is currently hard-coded in the wrapper that run the simulation workflow for the spectral simulations. We have developed the code assuming that the molecular structure (atoms and coordinates) is known to the wrapper, although in the SimDOME GUI the molecule and solvent will be specified by their common names (e.g., "water"). For this reason, additional code will be necessary to handle the conversion from the common name to the molecular structure. Since this mapping might be relevant for multiple wrappers, this could be part of the tool that processes the input from the user, that connects to a DB with these structures and fetches the correct one. 
Regarding the solvent, the actual COBRAMM implementation of the simulation workflow treat this information as a parameter to be chosen among a predefined list of options, specified by their common names. This choice was done because of two reasons. First, the simulation workflow requires the MM parameters for the solvent molecules, and these cannot be generated automatically but require a substantial amount of non-standardizable work. Second, the definition of the solvent molecules is a complex procedure that is done using external tools from the AMBER package, and thus belongs in a substantial measure to the third-party software. In this light, at the level of the wrapper it is not needed to know the molecular structure of the solvent.
The wrapper operates with a minimal ontology (see  file cobramm.ontology.yml) that defines the I/O information of the wrapper: the molecular structure (a list of atoms each specified by cartesian coordinates and atomic symbol), the name of the solvent, temperature and pressure, the absorption spectrum (a list of absorption values, each considered as a couple wavelength-intensity).
The wrapper is now composed of three python code files:
* osp/wrapper/cobrammwrapper/simulation_engine.py which is the syntactic layer of the wrapper, i.e. a class that specifies the operations needed to run the simulation workflow in the COBRAMM language.
* osp/wrapper/cobrammwrapper/cobramm_session.py which is the interoperability layer, that connects the syntactic layer to the osp core.
* examples/simulation_example.py, which is a simple example in which we compute semantically the linear absorption spectrum of a small molecule.
In the following we will describe in detail the syntactic layer and the example. The interoperability layer has been built following the standard rules of the SimPhoNy wrapper implementation and thus requires no particular description. 

## The syntactic layer
The concept behind the CobrammSimulationEngine class that is contained in the syntactic layer is quite simple. It contains a few attributes to store the information coming from the semantic layer: 
* CobrammSimulationEngine.atoms: a dictionary to store the definition of the molecular structure.
* CobrammSimulationEngine.input_setup: a dictionary to store the other input properties (temperature, pressure, solvent).
* CobrammSimulationEngine.accuracy: an integer that specify the accuracy level of the calculation (supposed to be 1, 2 or 3).
* CobrammSimulationEngine.grid and CobrammSimulationEngine.spectrum: list to store the linear absorption spectrum.
These attributes are set and deleted with appropriate methods (lines 49-64). 
The CobrammSimulationEngine.run method runs the simulation workflow, using the input parameters stored in the attributes of the class. In lines 91-120, the accuracy level is mapped to a set of parameters that specify the overall accuracy and computation cost of the calculation: by increasing the value from 1 to 3, we get larger simulation boxes, longer simulation times, smaller time steps and more configuration samples.
Once the required input information is retrieved and the relevant parameters are set, we can start with the simulation workflow. From line 131 on, the function contains the whole simulation workflow that has been already described in the section 2). This CobrammSimulationEngine class does not directly contains the operations of the simulation workflow. Instead, it links the code for the various simulation steps that has been included directly in COBRAMM. As an example, let us briefly consider the MM optimization and equilibration step, lines 137-189. This part of the syntactic layer makes use of the classes, AmberInput and AmberCalculator, that are included in COBRAMM and constitute a smart limited-purpose wrapper to Amber MD.
Inspecting the import statements at the beginning of simulation_engine.py we have a list of the parts of COBRAMM that are used by this layer: 
* The classes AmberCalculator and AmberInput constitute the COBRAMM internal wrapper to Amber.
* The classes CobrammCalculator and CobrammInput constitute a wrapper to COBRAMM, available in the same COBRAMM code.
* The class HarmonicSampling contains the Wigner sampling code.
After the CobrammSimulationEngine.run method has been run, an attribute flag CobrammSimulationEngine.executed is set to True and the spectrum computed with the simulation workflow becomes available for the interoperability layer.

## The example
In the example (simulation_example.py) contained in the actual wrapper implementation, we use the COBRAMM wrapper to compute the linear absorption spectrum of Formaldehyde in water at 1 atm and 300 K. This is obviously done semantically using SimPhoNy.
Lines 15-18: we include in the material the four atoms of the formaldehyde specifying their cartesian coordinates in Angstrom.
Lines 21-23: we specify the physical conditions, i.e., temperature, pressure, and solvent.  
Lines 31-33: we add the material specified above to the COBRAMM wrapper.
Line 36: we run the simulation.
Line 39: we extract the simulated spectrum.

## How to install
1. Clone the wrapper repository.
2. Install the ontology.
``` 
~/comsolwrapper$ pico install cobramm.ontology.yml
```
3. In order to install the wrapper and all its requirements please run the following commands. 
``` 
~/comsolwrapper$ python setup.py install 
```

## How to test
After the installation
1. Go inside the tests folder
2. Execute the command

```python
python -m unittest
```

The suite is composed of three different tests, it will takes approximately 1 hour to run them all. If you prefer you can just run one of them by specifying the name of the file.

```python
python -m unittest test_engine.py
```

### Structure of the tests
The suite is composed of three different tests
* **test_end_to_end**: It is intended to test the end to end interaction between application, wrapper and the simulation engine.
* **test_session**: It comprises both integrations tests related to the interaction with the simulation engine and unit tests of single methods (part of the wrapper API).
* **test_engine**: This test suite is intended to unit test the simulation engine.
