# A two-cell circuit around a tumour in 3D Matrigel

Code for the simulation of sensor and actuator T cells searching a 96-well gel, with and without chemokine recruitment.

A tumour hemisphere of radius 1 mm sits on the floor of a 96-well plate under 80 µL of Matrigel. Sensor cells and actuator cells enter from the top of the gel and crawl. A sensor cell is active while it is inside the tumour. An actuator cell switches on permanently the first time it touches an active sensor. The readout is the fraction of actuators switched on. Two conditions are compared, one where the cells cannot read the tumour's CXCL12 and one where they can.

## Running it

Julia 1.10 or later.

```
julia --project=. -e 'import Pkg; Pkg.instantiate()'
julia --project=. figure09_matrigel_circuit.jl
```

This writes three files into `figures/`.

| file | holds |
|:--|:--|
| `figure09_matrigel_circuit.png` | the figure |
| `figure09_matrigel_circuit.pdf` | the same figure as vector art, with editable text |
| `figure09_matrigel_circuit_data.csv` | every record the figure is built from |

The run takes about 20 minutes on one core. Panel D is the slow part, since it sweeps the sensor number to 40,000 over 5 independent gels per point.

## The files

| file | holds |
|:--|:--|
| `figure09_matrigel_circuit.jl` | the entry point |
| `core_gel.jl` | the gel, the tumour, the chemokine field, the motility and the circuit |
| `makie_common.jl` | the CairoMakie theme and the save helpers |
| `common.jl` | the parula palette and the data-file schema |

## The data file

`figure09_matrigel_circuit_data.csv` is one long table. A `panel` column says which panel a row belongs to and a `record` column says what the row is.

| `record` | one row per | carries |
|:--|:--|:--|
| `trajectory` | simulation step of each recorded cell | its position and the time |
| `trace_switch` | recorded actuator that switched on | the time it switched |
| `switch_time` | actuator in a population gel | the time it switched on, blank if it never did |
| `replicate` | individual gel behind panel D | its seed, its count switched on and its fraction |

Nothing in the file is a summary. Panel C's curve is the empirical distribution of the switch times the file holds, and panel D's means and error bars are computed from the individual gels, which the `replicate` column separates. The full parameter list sits above the data as comment lines, which `pandas.read_csv(comment='#')` and `read.csv(comment.char='#')` both skip.

## Parameters

Every cell and chemokine parameter is in the `CHEMO` and `GEL` tables at the top of `core_gel.jl`, and every one is written into the header of the data file. The chemokine decay rate is the measured CXCL12 half-life of 26.5 min under DPPIV digestion (Takekoshi et al. 2012). The diffusion coefficient is 500 µm²/min from single-molecule tracking in collagen (Miller et al. 2018). The speed is 10 µm/min from two-photon tracking in lymph node (Miller et al. 2003). The concentration at the tumour surface, 22.37 nM, comes from the tumour's own CXCL12 ELISA.
