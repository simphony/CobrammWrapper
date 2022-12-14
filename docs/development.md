# Development

This section is aimed at developers and advanced users. It briefly explains the
wrapper implementation and how to run the tests.

## Wrapper implementation
The necessary code to carry out the linear absorption simulation workflow is 
included in the [COBRAMM engine](https://site.unibo.it/cobramm/en),
which is available independently of the wrapper in the form of some utility 
scripts.

The wrapper has been developed according to the indications outlined by the 
[SimPhoNy](https://github.com/simphony/simphony-osp) developers. Regarding the
solvent, the COBRAMM implementation of the simulation workflow treats this 
information as a parameter to be chosen among a predefined list of options, 
specified by their common names. This choice was done because of two reasons.
First, the simulation workflow requires the MM parameters for the solvent 
molecules, and these cannot be generated automatically but require a 
substantial amount of non-standardizable work. Second, the definition of the
solvent molecules is a complex procedure that is done using external tools from
the AMBER package, and thus belongs in a substantial measure to the third-party
software. In this light, at the level of the wrapper it is not needed to know
the molecular structure of the solvent.

The wrapper operates with a minimal ontology (see file `cobramm.ontology.yml`)
that is used to define the inputs and outputs of the wrapper: the molecular
structure (a list of atoms each specified by cartesian coordinates and atomic
symbol), the name of the solvent, temperature and pressure, the absorption 
spectrum (a list of absorption values, each considered as a couple 
wavelength-intensity).

The wrapper is composed of three Python code files:
* `osp/wrapper/cobrammwrapper/simulation_engine.py` which is the 
  [syntactic layer](https://simphony.readthedocs.io/en/v3.9.0/detailed_design.html#syntactic-layer)
  of the wrapper, i.e. a class that specifies the operations needed to run the
  simulation workflow in the COBRAMM language.
* `osp/wrapper/cobrammwrapper/cobramm_session.py` which is the 
  [interoperability layer](https://simphony.readthedocs.io/en/v3.9.0/detailed_design.html#interoperability-layer),
  that connects the syntactic layer to the OSP-core.
* `examples/simulation_example.py`, which is a simple example in which the 
  linear absorption spectrum of a small molecule is computed.

The syntactic layer is further described in detail next. The interoperability
layer has been built following the 
[SimPhoNy Wrapper API specification](https://simphony.readthedocs.io/en/v3.9.0/wrapper_development.html)
and thus requires no specific description. 

### Syntactic layer
The concept behind the `CobrammSimulationEngine` class that is contained in the
syntactic layer is quite simple. It contains a few attributes to store the
information coming from the semantic layer: 

* `CobrammSimulationEngine.atoms`: a dictionary to store the definition of the 
  molecular structure.
* `CobrammSimulationEngine.input_setup`: a dictionary to store the other input
  parameters (temperature, pressure, solvent).
* `CobrammSimulationEngine.accuracy`: an integer that specify the accuracy
  level of the calculation (supposed to be 1, 2 or 3).
* `CobrammSimulationEngine.grid` and `CobrammSimulationEngine.spectrum`: list
   to store the linear absorption spectrum.

These attributes are set and deleted with appropriate methods (lines 56-71). 

```{include} ../osp/wrappers/cobrammwrapper/simulation_engine.py
   :start-line: 55
   :end-line: 72
   :code: python
   :number-lines: 56
```

The `CobrammSimulationEngine.run` method runs the simulation workflow, using 
the input parameters stored in the attributes of the class. The accuracy level
is mapped to a set of parameters that specify the overall accuracy and 
computation cost of the calculation: by increasing the value from 1 to 3, we 
get larger simulation boxes, longer simulation times, smaller time steps and 
more configuration samples.

Once the required input information is retrieved and the relevant parameters 
are set, we can start with the simulation workflow. For that, the function also
contains the whole simulation workflow that has been already described in 
section 2). This `CobrammSimulationEngine` class does not contain the code that
executes the individual operations of the simulation workflow. Instead, it 
calls the code for the various simulation steps that has been included in the
[COBRAMM engine](https://site.unibo.it/cobramm/en).

As an example, let us briefly consider the MM optimization and equilibration
step. This part of the syntactic layer makes use of the classes `AmberInput`
and `AmberCalculator` that are included in COBRAMM and constitute a smart 
limited-purpose wrapper to Amber MD. Inspecting the import statements at the 
beginning of `simulation_engine.py` we have a list of the parts of COBRAMM that
are used by this layer: 

* The classes `AmberCalculator` and `AmberInput` constitute the COBRAMM
  internal wrapper to Amber.
* The classes `CobrammCalculator` and `CobrammInput` constitute a wrapper to 
  COBRAMM, available in the same COBRAMM code.
* The class `HarmonicSampling` contains the Wigner sampling code.

After the `CobrammSimulationEngine.run` method has been run, an attribute flag 
`CobrammSimulationEngine.executed` is set to `True` and the spectrum computed 
with the simulation workflow becomes available for the interoperability layer.

## Tests

```{include} ../README.md
   :start-after: <!---tests-start-b16e3913-->
   :end-before: <!---tests-end-b16e3913-->
```
