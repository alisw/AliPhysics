// $Id$

/*! 

\page README_raw Raw data
 
\section raw_s1 How to read & decode raw data 

These macros can read & decode DDL files and root files.
DATE files are no more supported.

For tracker raw data, full output
<pre>
.includepath $ALICE_ROOT/STEER
.includepath $ALICE_ROOT/MUON
.includepath $ALICE_ROOT/RAW
.L $ALICE_ROOT/MUON/MUONRawStreamTracker.C+
MUONRawStreamTrackerExpert(rawFileName, maxEvent, firstDDL, lastDDL)
</pre>


For tracker raw data, just digit output
<pre>
.includepath $ALICE_ROOT/STEER
.includepath $ALICE_ROOT/MUON
.includepath $ALICE_ROOT/RAW
.L $ALICE_ROOT/MUON/MUONRawStreamTracker.C+
MUONRawStreamTrackerSimple(rawFileName, maxEvent, firstDDL, lastDDL)
</pre>

For trigger raw data, full output
<pre>
.includepath $ALICE_ROOT/STEER
.includepath $ALICE_ROOT/MUON
.includepath $ALICE_ROOT/RAW
.L $ALICE_ROOT/MUON/MUONRawStreamTrigger.C+ 
MUONRawStreamTrigger(maxEvent, firstDDL, lastDDL, rawFileName)
</pre>

For trigger raw data, local response output
<pre>
.includepath $ALICE_ROOT/STEER
.includepath $ALICE_ROOT/MUON
.includepath $ALICE_ROOT/RAW
.L $ALICE_ROOT/MUON/MUONRawStreamTrigger.C+ 
MUONRawStreamTriggerSimple(maxEvent, firstDDL, lastDDL, rawFileName)
</pre>


Default wise the macro read all DDL files from the current directory for 1000 events.
For root file rawFileName must end with .root. For DDL files you have to specified the directory 
where the raw0...n subdirectories are located:
<pre>
MUONRawStreamTracker(..)(maxEvent, firstDDL, lastDDL, "$YOUR_WORKING_DIRECTORY/"); //Do not forget the slash at the end!
</pre>

This chapter is defined in the READMEraw.txt file.

*/

