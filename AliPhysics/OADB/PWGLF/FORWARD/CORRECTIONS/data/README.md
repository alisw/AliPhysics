Corrections for PWGLF/FORWARD Analysis
======================================

This directory contain 2 files:

* fmd_corrections.root: Database of corrections used by the FMD part
  of PWGLF/FORWARD
* spd_corrections.root: Database of corrections used by the SPD part
  of PWGLF/FORWARD

General Remarks
---------------

The code to store and retrieve lives in AliCorrectionManagerBase with
specific sub-classes AliFMDCorrectionManager and
AliSPDCorrectionManager. The base class - and it's contained classes -
has all the functionality for retrieving and storing the corrections
in a TTree (one tree for each correction) while the sub-classes merely
registers which corrections exists and provide convenience functions.

The data is stored in TTrees - one for each correction. Each entry in
the trees correspond to a single instance of a correction. The
correction object itself is wrapped in a AliOADBForward::Entry
object, which provide meta data such as

* run number
* collision system (1: pp, 2: PbPb, 3: pPb)
* collision energy in GeV
* L3 magnetic field in kG
* whether the correction can be used for satellite analysis
* whether the correction is for simulated or real data
* date of submission
* author of the submission.

The code to handle these trees and the entries lives in
AliOADBForward.

Entries in the trees are indexed by one or more of the above mentioned
meta fields.  Which fields are used for a given correction is defined
by the sub-classes of AliCorrectionManagerBase.

In general, a query will only succeed if all fields of the query
matches an entry.  If this returns more than one result, the most
recent one is used.  However, there are some exceptions to this.

Run Number Selection
--------------------
In general, one cannot expect to get an entry with the exact run
number queried.  AliOADBForward therefor defines a set of policies for
how to match a given run number

* EXACT - run number must match exactly
* NEAR - the run number of the selected entries must not be too far
  from the run number given in the query.
* OLDER - entries with a run number less than or equal to the queried
  run number are selected
* NEWER - entries with a run number greater than or equal to the
  queried run number are selected

Which mode to use for a given correction is defined on creation of the
database, but can be overridden later on.

If no entry is found with the specified run number and mode, then the
correction may allow for a fall-back option.  These are entries with
run number equal to 0.  In this case, all other active meta fields
must match.

Fall-back Mode
--------------

If a general fall-back mode is enabled, then if no entry was selected
using the full set of fields, then we can ignore the collision
energy.  If this still does not return anything, we can ignore all
fields, which means we will select the most recent entry -
irrespective of other meta fields.



<!--  Local Variables: -->
<!--   mode: markdown -->
<!--  End --> 
<!--  LocalWords:  AliCorrectionManagerBase AliFMDCorrectionManager -->
<!--  LocalWords:  AliSPDCorrectionManager TTree TTrees PbPb pPb GeV
<!--  LocalWords:  AliOADBForward kG
 -->
