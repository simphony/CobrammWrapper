# Usage

The simulation inputs that the wrapper needs to be provided with are the 
molecular structure, solvent, temperature and pressure. Such simulation inputs
must be instantiated as CUDS objects in the wrapper session using the
namespace from the wrapper ontology. The `cobramm` namespace and the wrapper 
session can be imported as follows:

```python
from osp.core.namespaces import cobramm
from osp.wrappers.cobrammwrapper import CobrammSession
```

The pattern of CUDS objects that the COBRAMM wrapper expects as input is 
depicted in the figure at the end of this section.


The input information is employed to run the following workflow:

1) Solvent-solute MM optimization, thermalization and equilibration dynamics 
   (Amber force field for the solvent and GAFF for the solute) with the charges
   estimated by the AM1-BCC method from antechamber.

2) Set-up of the high-medium-low layers for QM/MM COBRAMM driven computations 
   from the last MD snapshot produced in 1) (which is supposed to return the 
   better solvent organization): the solvent MM shell is selected as a
   spherical droplet around the chromophore,

3) QM(DFT: b3lyp, STO-3G)/MM optimization for locating the representative 
   ground state equilibrium structure of the solvated chromophore,

4) QM(DFT: b3lyp, STO-3G)/MM frequency calculation on the QM/MM minimum 
   optimized in 3),

5) Wigner ensemble generation,

6) Single point QM(TDDFT: cam-b3lyp, STO-3G, nroot=5)/MM calculations on top of
   each equilibrated Wigner snapshots found in 7),

7) Sum over gaussian-functions centered at the excitation energy values found 
   in 6) (with the gaussian maximum set at the oscillator strength value and a
   fixed 0.2 ev standard deviation) divided by the number of Wigner samples.

The Low/Medium/High accuracy level do set parameters for MD simulations 
(boxsize, cutoff, initial optimization steps, time step of the MD run, 
heating time, time of the final equilibration steps) and QM/MM simulations 
(solvent droplet radius, mobile layer radius, number of samples produced in the
Wigner sampling) of increasing complexity, accuracy and computational cost.

After the simulation is run, the simulated spectrum is available as a 
collection CUDS objects attached to the `SYSTEM` object. Their structure is
depicted in the figure below.

<figure style="display: table; text-align:center; margin-left: auto; margin-right:auto">

![COBRAMM Wrapper input](./static/graph_pattern.drawio.svg)

<figcaption style="display: table-caption; caption-side: bottom; text-align:center">

_Diagram showing the pattern of input CUDS objects that the COBRAMM wrapper
expects to find in the session's knowledge graph and the CUDS objects that are
produced as simulation outputs._

</figcaption>
    
</figure>

## Example

An example of how to use the wrapper is available in 
`examples/simulation_example.py`. In the example, the COBRAMM wrapper is used
to compute the linear absorption spectrum of Formaldehyde in water at 1 atm and
300 K.

```{include} ../examples/simulation_example.py
   :code: python
   :number-lines: 1
```

Lines 22-25: we include in the material the four atoms of the formaldehyde 
specifying their cartesian coordinates in Angstroms.

Lines 30-39: we specify the physical conditions, i.e., temperature, pressure,
and solvent.  

Lines 47-49: we add the material specified above to the COBRAMM wrapper.

Line 57: we run the simulation.

Line 60: we extract the simulated spectrum.
