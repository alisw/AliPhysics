// $Id$

/*! 

\page README_raw Raw data
 
\section raw_s1 How to read & decode raw data 

The macros MUONRawStreamTracker.C and MUONRawStreamTrigger.C can 
be used to read & decode DDL files and root files.
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

\section raw_s2 Timing of the raw data decoders

The macro MUONTimeRawStreamTracker.C can used to check the timing (speed) performance 
of the existing offline decoder for the tracker DDLs (AliMUONRawStreamTracker) and also 
for the new high performance decoder (AliMUONRawStreamTrackerHP). It can be invoked 
as follows:

<pre>
 $ aliroot
.L $ALICE_ROOT/MUON/MUONTimeRawStreamTracker.C+
 MUONTimeRawStreamTracker(filename, maxEvent);
</pre>

where \em filename is the name of a file containing the raw data, or alternatively
the directory containing rawX (X being an integer) paths with the raw DDL
data. The \em maxEvent value is the maximum event to process (default set to
1000). Thus the macro will time the algorithm for all events in the range
[0 .. maxEvent-1].


This chapter is defined in the READMEraw.txt file.

*/

