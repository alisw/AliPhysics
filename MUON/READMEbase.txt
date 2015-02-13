// $Id$

/*! \page README_base MUON Data definition and access

Both the simulation and reconstruction use containers (called stores in the MUON jargon) 
to hold the data we're dealing with: hits, (s)digits, trigger, clusters, tracks and 
trigger tracks. All those stores share some commonalities, in particular with respect
to how they are read/written from/to TTree. @ref AliMUONVStore "More..."


\section base_s1 How to dump the content of Root data files 

To check the content of Root data files, the AliMUON*DataInterface classes
provides the functions to produce an ASCII output on the screen
which can be redirected on the file:

for MC information, use AliMUONMCDataInterface :

<pre>
> aliroot (or root with just the loading of MUON libs, see loadlibs.C)
root [0] AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
root [1] AliMUONMCDataInterface mcdi("galice.root");
root [2] mcdi.DumpKine(5);       > dump.kine
root [3] mcdi.DumpHits(5);       > dump.hits
root [4] mcdi.DumpTrackRefs(5);  > dump.trackrefs
</pre>

for all other information, use AliMUONDataInterface :

<pre>
> aliroot
root [0] AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
root [1] AliMUONDataInterface di("galice.root");
root [2] di.DumpDigits(5);     > dump.digits
root [3] di.DumpSDigits(5);    > dump.sdigits
root [4] di.DumpRecPoints(5);  > dump.recpoints
root [5] di.DumpTrigger(5); > dump.rectrigger
</pre>

Remind that during simulation and reconstruction two 
differents galice.root are generated: one for the generation 
(simulation) and other during the reconstruction.

If you open the wrong galice.root file you could get:
<pre>
root [0] AliMUONMCDataInterface mcdi("galice.root");
root [1] mcdi.DumpKine(5);
W-AliRunLoader::GetEvent: Stack not found in header
E-TFile::TFile: file ./Kinematics.root does not exist
</pre>

\section basee_s2 Macro MUONCheckDI.C

MUONCheckDI.C performs a consistency check on the methods of the 
AliMUONMCDataInterface and AliMUONDataInterface classes. There are several 
helper methods in these classes which make it easier to fetch data, which 
means there are at least two ways of fetching the data within the same class 
interface. The macro checks to see that the results given by these different 
methods are identical, as they should be.

The macro also inherently exercises the AliMUONMCDataInterface and 
AliMUONDataInterface classes and should be run after any modifications to 
these classes to see if things still work. Putting it another way: 
MUONCheckDI.C is a testing facility for developers of these two classes.

This chapter is defined in the READMEbase.txt file (although it describes
also the code in the evaluation library.)

*/
