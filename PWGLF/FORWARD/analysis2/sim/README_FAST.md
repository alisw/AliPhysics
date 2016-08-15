# Fast simulations anchored in Runs {#pwglf_fwd_fastsim}

Christian Holm Christensen
<cholm@nbi.dk>

## Introduction

The files in this archive represents a generic approach to running
FAST simulations - that is - only generator level data (no swimming
through ALICE).  The idea is that the user/committer only specifies
the run number for which an anchor run is needed, and possibly an
event generator.

## Usage:

    RunFast.C(URL)


Here `URL` is a string specifying the job. The `URL` argument has the
format 

    PROTOCOL://[HOST]/?[OPTIONS]

where `PROTOCOL` can be one of

* `local`: Run a single-thread simulation on the local machine.
* `lite:`: Run a multi-thread simulation on the local machine, using
  the ProofLite functionality
* `proof`: Run a multi-thread simulation on a PROOF farm specified by
  `HOST`.

The `HOST` part only makes sense for `PROTOCOL=proof`.

`OPTIONS` is a list of options separated by ampersands &amp;, and can contain

* `run=RUN`: Specify which run to anchor in.  Beam conditions (beam
  species, collision energy, and such is taken from GRP for that
  run).
* `events=NUMBER`: Number of events to simulate
* `eg=NAME`: The event generator to use. Note, if this is not
  specified, then a default EG based on the run parameters is
  chosen. Names are typically of the form

> _name_[_sub-type_[_sub-sub-type_]]

  where _name_ could be `Hijing`, `Pythia`, `DpmJet`, and so on, while
  _sub-type_ and _sub-sub-type_ specifies tunes or variants of the
  main model.
* `monitor=SECONDS`: Specifies that monitor histograms should be
  updated every `SECONDS`.  If 0 or smaller or not specified, then no
  monitoring is done.
* `b=RANGE`: Specifies the impact parameter range.
* `override=LIST`: Specifies a comma-separated (,) list of settings
  that should be used rather than the ones being read from GRP. The
  can contain.  Each item in the list is a key, value pair separated
  by a colon (:)
  * `beamEnergy:BEAMENERGY`: Total beam energy (not collision energy)
    in GeV
  * `energy=ENERGY`: Collision energy in GeV
  * `period=IDENTIFIER`: ALICE running period (e.g., `LHC15a`)
  * `run=NUMBER`: Run number
  * `beam1.a=ATOMIC_NUMBER`: Atomic number of particles in beam1
  * `beam1.z=CHARGE`: Charge of particles in beam1
  * `beam2.a=ATOMIC_NUMBER`: Atomic number of particles in beam2
  * `beam2.z=CHARGE`: Charge of particles in beam2
* `save=MODE`: Only relevant for Proof(Lite). Values can be
  * `none`: Do not retrieve the final `galice.root` and
    `Kinematics.root` files.
  * `split`: Return the final `galice.root` and `Kinematics.root`
    files - one for each worker.  The files are moved to a
    sub-directory on the client machine, and an collection of `TUrl`
    objects is written to `index.root`.  One can easily define a chain
    using this information. 
  * `merge`: Does not work
  
`OPTIONS` can also contain options for the execution environment, such
as `workers=N` for ProofLite.

### Example:

    RunFast.C("lite:///?events=10000&eg=default&run=138190")

will make 10000 events using the default generator (Hijing), anchored
in run 138190 (LHC10h, PbPb @ 2.76TeV).  The output file will be

    Hijing_000138190_AA_02760_10k.root

## Analysis:

    ProcessFast.C(URL,OUTPUT)

where, again the URL argument specifies the input and execution
environment.  It has the format

    PROTOCOL://[HOST]/[INPUT]?[OPTIONS]


where `PROTOCOL` and `HOST` is as above. The `INPUT` part specifies
the input file to read. `OPTIONS` is again a list of options separated
by ampersands.  Valid options are

* `events=NUMBER`: Number of events to analyze
* `type=NAME` Type of analysis to perform.  Defined types are
   * `INEL`
   * `NSD`
   * `CENTV0M`
   * `CENTV0A`
   * `CENTV0C`
   * `MULTRefMult00d80`
   * `MULTRefMult00d50`

`OPTIONS` can also contain options for the execution environment, such
as `workers=N` for ProofLite.

The `OUTPUT` argument specifies which ROOT file to write the results
to.

### Example

Using the output from the above example

    ProcessFast.C("lite:///${PWD}/Hijing_000138190_AA_02760_10k.root?events=1000&type=CENTV0M","out.root")

will analyze 1000 events from the above run, and build dN/deta per
centrality bin, using a simulated V0M centrality estimator. 

### Centrality selection

The centrality estimation is simulated by counting the number of
charged particles in a given region - e.g., the eta region
corresponding to the V0 acceptance.  The distribution of that signal
is then integrated from the top to the bottom and the incremental
integral is written to a histogram.  Furthermore, the estimator
quantity is saved as a separate branch in the output tree.

All this happens during `RunFast.C`.  When looping over the data with
`ProcessFast.C` we retrieve the integral histogram from the input
file, and for each event we use the stored estimator quantity to look
up the centrality of that event.

In this way, we do not need an additional pass to extract the
centrality information, and we are sure that the data is
self-consistent because of we replay the estimator quantity using the
stored value. 

## Technicalities

The fast simulation uses the same strategy as for full simulations.
For more details see README.md in this directory.

<!-- Local Variables: -->
<!--   mode: markdown -->
<!--   ispell-dictionary: "british" -->
<!-- End: -->
<!--  LocalWords:  RunFast multi ProofLite GRP eg Hijing Pythia TeV
 -->
<!--  LocalWords:  DpmJet ProcessFast PbPb LHC PWD INEL NSD CENTV dN
 -->
<!--  LocalWords:  MULTRefMult deta LocalWords
 -->
