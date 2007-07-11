// $Id$

/*! 

\page README_raw README raw
 
\section raw_s1 How to read & decode raw data 

These macros can read & decode DDL files, root and DATE files.
Nevertheless for the two latter, aliroot has to be compiled with DATE.

For tracker raw data
<pre>
.includepath $ALICE_ROOT/STEER
.includepath $ALICE_ROOT/MUON
.includepath $ALICE_ROOT/RAW
.L $ALICE_ROOT/MUON/MUONRawStreamTracker.C+
MUONRawStreamTracker(maxEvent, firstDDL, lastDDL, rawFileName)
</pre>

For trigger raw data
<pre>
.includepath $ALICE_ROOT/STEER
.includepath $ALICE_ROOT/MUON
.includepath $ALICE_ROOT/RAW
.L $ALICE_ROOT/MUON/MUONRawStreamTrigger.C+ 
MUONRawStreamTrigger(maxEvent, firstDDL, lastDDL, rawFileName)
</pre>

Default wise the macro read all DDL files from the current directory for 1000 events.
For root file rawFileName must end with .root, for date file rawFileName 
must be no extention. For DDL files you have to specified the directory 
where the raw0...n subdirectories are located:
<pre>
MUONRawStreamTracker(maxEvent, "$YOUR_WORKING_DIRECTORY/"); //Do not forget the slash at the end!
</pre>

*/

