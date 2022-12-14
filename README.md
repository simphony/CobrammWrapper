# COBRAMM SimPhoNy Wrapper

<!---introduction-start-4a4f2d70-->

[SimPhoNy](https://github.com/simphony/simphony-osp) wrapper for the 
[COBRAMM engine](https://site.unibo.it/cobramm/en) developed at 
[UNIBO](https://www.unibo.it). 

<!---introduction-end-4a4f2d70-->

<img src="https://site.unibo.it/cobramm/en/@@images/97dcd2e2-e470-4cf4-a0e7-4d15b2e44018.png" width="200">

<!---introduction-start-f8a2dfb4-->

The COBRAMM SimPhoNy wrapper enables modelling the electronic (UV-VIS) linear
absorption (LA) spectrum of a solvated neutral molecule in its singlet ground
state. The [QM/MM COBRAMM engine](https://site.unibo.it/cobramm/en) developed 
at UNIBO is employed. COBRAMM performs hybrid QM/MM calculations by linking 
the [Gaussian](https://gaussian.com/) (QM) and [Amber](https://ambermd.org/)
(MM) software packages. 

An explicit atomistic MM model (Amber force field) is used for the solvent, 
while an explicit atomistic QM model (DFT/TDDFT) is used for the solvated 
chromophore. A classical equilibration of the solvent-solute system followed 
by a QM/MM optimization with Wigner sampling and QM/MM excited state 
calculations is employed to generate the LA spectrum and its bandshape. 
Vibrational progressions are not considered at this stage. The wrapper makes 
available to the user three different schemes of increasing complexity and 
computational time (Low, Medium, High) that do set parameters for MD 
simulations (boxsize, cutoff, initial optimization steps, time step of the 
MD run, heating time, time of the final equilibration steps) and QM/MM 
simulations (solvent droplet radius, mobile layer radius, number of samples 
produced in the Wigner sampling) of increasing complexity, accuracy and 
computational cost. Have a look at the
[documentation](https://cobrammwrapper.readthedocs.io) for further
information and the workflow details.

*Contact:* [Marco Garavelli](mailto:marco.garavelli@unibo.it) ([UNIBO](https://www.unibo.it))

<!---introduction-end-f8a2dfb4-->

## Installation

<!---installation-start-d47a8de6-->

This wrapper requires a working installation of COBRAMM 2.3. Please follow the
[installation guide](https://gitlab.com/cobrammgroup/cobramm/-/wikis/Installation-Guide)
and make sure that you have added the required directories to the `PATH` 
environment variable as described on it.

After that, clone the wrapper repository

```shell
git clone https://gitlab.cc-asp.fraunhofer.de/simphony/wrappers/cobrammwrapper.git
```

install the wrapper

```shell
pip install cobrammwrapper
```

and install the ontology.

```shell
pico install cobrammwrapper/cobramm.ontology.yml
```

<!---installation-end-d47a8de6-->

## Documentation

Visit the [COBRAMM Wrapper's documentation](https://cobrammwrapper.readthedocs.io)
to learn how to use the wrapper. You may also have a look at the
[`examples` folder](https://github.com/simphony/cobrammwrapper/tree/master/examples)
for additional examples.

The documentation may also be built and displayed locally running the commands
below.

```shell
pip install -r docs/requirements.txt
sphinx-autobuild --host 127.0.0.1 docs build
```

## Testing

<!---tests-start-b16e3913-->

If you wish to verify that the wrapper is working as expected, you can run its
test suite. Install the wrapper first and, after the installation:

1. Go inside the tests folder `cd tests`
2. Execute the command `python -m unittest`

This will run all tests in the test suite, and it will take approximately 
1 hour to run them all. If you prefer you can just run one of them by 
specifying the name of the file.

```python
python -m unittest test_engine.py
```

The suite is composed of three different tests
- **test_end_to_end**: It is intended to test the end to end interaction 
  between application, wrapper and the simulation engine.
- **test_session**: It comprises both integration tests related to the 
  interaction with the simulation engine and unit tests of the methods that are
  part of the SimPhoNy wrapper API.
- **test_engine**: This test suite contains unit tests for the simulation 
  engine.

<!---tests-end-b16e3913-->


