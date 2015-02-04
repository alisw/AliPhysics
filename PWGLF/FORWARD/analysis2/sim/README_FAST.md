# Fast simulations anchored in Runs

Christian Holm Christensen
<cholm@nbi.dk>

## Introduction

The files in this archive represents a generic approach to running
FAST simulations - that is - only generator level data (no swimming
through ALICE).  The idea is that the user/committer only specifies
the run number for which an anchor run is needed, and possibly an
event generator.

## Usage:

    RunFast`.C(URL,OPTIONS)


Here `URL` is a string specifying the job, and `OPTIONS` is a string
of optional options.

The `URL` argument has the format

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

`OPTIONS` can also contain options for the execution environment, such
as `workers=N` for ProofLite.

## Technicalities

The fast simulation uses the same strategy as for full simulations.
For more details see README.md in this directory.

Local Variables:
  mode: markdown
End:







 
