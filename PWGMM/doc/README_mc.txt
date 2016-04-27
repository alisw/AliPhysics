/*! \page pwgmm_mc Monte Carlo event generators

# available generators

generator | in package | interface | versions
--- | --- | --- | ---
AMPT | AliRoot | n/a | v1.25t3/v2.25t3
DIPSY | ThePEG | HepMC | v2015-08-11 (private comm.)
DPMJET | AliRoot | n/a | 3.0.5
EPOS | EPOS | EPOS reader | 3.111
EPOS LHC | CRMC | HepMC | 1.5.4
Herwig | AliRoot | n/a | 6.507, 6.510
Herwig++ | ... | ... | ...
HIJING | AliRoot | n/a | 1.35, 1.36
JEWEL | JEWEL | HepMC | 2.0.2
POWHEG | POWHEG | ... | ...
Pythia 6 | AliRoot | n/a | 6.4.21, 6.4.24, 6.4.28
Pythia 8 | AliRoot, pythia | n/a, HepMC | 8.175, 8.205, 8.210
Sherpa | SHERPA | HepMC | ...
Starlight | AliRoot | n/a | r193

# on-the-fly LEGO trains

...

## internal generators

Controlled by dedicated Add... macros

## external generators

Use generic macro to add generator:
<pre>
ANALYSIS/macros/train/AddMCGenExtExec.C
</pre>
then define the script for a specific generator:
```
((AliGenExtExec*) generator)->SetPathScript("$ALICE_PHYSICS/PWG/MCLEGO/ThePEG/gen_dipsy_ropes.sh");
```

*/
