/*! \page pwgmm_mc Monte Carlo event generators

# available generators

generator | in package | versions
--- | --- | ---
AMPT | AliRoot | v1.25t3/v2.25t3
DIPSY | ThePEG | v2015-08-11 (private comm.)
DPMJET | AliRoot | 3.0.5
EPOS | EPOS | 3.111
EPOS LHC | CRMC | 1.5.4
Herwig | AliRoot | 6.507, 6.510
HIJING | AliRoot | 1.35, 1.36
JEWEL | JEWEL | 2.0.2
Pythia 6 | AliRoot | 6.4.21, 6.4.24, 6.4.28
Pythia 8 | AliRoot/pythia | 8.175, 8.205, 8.210
Starlight | AliRoot | r193

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
