// $Id$

/*! 

\page README_raw MUON Raw data
 
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

There is also a high performance decoder available, which can be run as follows
for full output:
<pre>
.includepath $ALICE_ROOT/STEER
.includepath $ALICE_ROOT/MUON
.includepath $ALICE_ROOT/RAW
.L $ALICE_ROOT/MUON/MUONRawStreamTracker.C+
MUONRawStreamTrackerHPExpert(rawFileName, maxEvent, firstDDL, lastDDL)
</pre>

And just for digit output like so:
<pre>
.includepath $ALICE_ROOT/STEER
.includepath $ALICE_ROOT/MUON
.includepath $ALICE_ROOT/RAW
.L $ALICE_ROOT/MUON/MUONRawStreamTracker.C+
MUONRawStreamTrackerHPSimple(rawFileName, maxEvent, firstDDL, lastDDL)
</pre>

The MUONRawStreamTracker.C macro also provides other alternative implementations
for fetching the decoded data from the high performance decoder's interface.
These generate the same output, but show how to write code to fetch the data in
various ways from the interface. Developers should consult the macro as an example
of how to use the interface. The alternate methods are called:
<pre>
MUONRawStreamTrackerHPExpert2()
MUONRawStreamTrackerHPExpert3()
MUONRawStreamTrackerHPSimple2()
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
MUONRawStreamTriggerSimple(maxEvent, rawFileName)
</pre>

Similarly there is a high performance decoder available for the trigger DDLs,
which can be run as follows for full output:
<pre>
.includepath $ALICE_ROOT/STEER
.includepath $ALICE_ROOT/MUON
.includepath $ALICE_ROOT/RAW
.L $ALICE_ROOT/MUON/MUONRawStreamTrigger.C+
MUONRawStreamTriggerHPExpert(maxEvent, firstDDL, lastDDL, rawFileName)
</pre>

And just for local response output like so:
<pre>
.includepath $ALICE_ROOT/STEER
.includepath $ALICE_ROOT/MUON
.includepath $ALICE_ROOT/RAW
.L $ALICE_ROOT/MUON/MUONRawStreamTrigger.C+
MUONRawStreamTriggerHPSimple(maxEvent, rawFileName)
</pre>

The MUONRawStreamTrigger.C macro also provides other alternative implementations
for fetching the decoded data from the trigger high performance decoder's interface.
These generate the same output, but show how to write code to fetch the data in
various ways from the interface. Developers should consult the macro as an example
of how to use the interface. The alternate methods are called:
<pre>
MUONRawStreamTriggerHPExpert2()
MUONRawStreamTriggerHPExpert3()
MUONRawStreamTriggerHPSimple2()
</pre>


Default wise the macros read all DDL files from the current directory for 1000 events.
For root file rawFileName must end with .root. For DDL files you have to specified the directory 
where the raw0...n subdirectories are located:
<pre>
MUONRawStreamTracker(..)(maxEvent, firstDDL, lastDDL, "$YOUR_WORKING_DIRECTORY/"); //Do not forget the slash at the end!
</pre>

\section raw_s2 Timing of the raw data decoders

The macro MUONTimeRawStreamTracker.C and MUONTimeRawStreamTrigger.C can used to
check the timing (speed) performance of the existing offline decoders compared to
the new high performance decoders.
For the tracker DDLs the MUONTimeRawStreamTracker.C macro compares the timing of
AliMUONRawStreamTracker against the high performance AliMUONRawStreamTrackerHP
decoder.
Similarly the MUONTimeRawStreamTrigger.C macro compares the timing of the
existing AliMUONRawStreamTrigger decoder against the high performance
AliMUONRawStreamTrackerHP decoder.
The macros can be invoked as follows:

<pre>
 $ aliroot
.L $ALICE_ROOT/MUON/MUONTimeRawStreamTracker.C+
 MUONTimeRawStreamTracker(filename, maxEvent);
.L $ALICE_ROOT/MUON/MUONTimeRawStreamTrigger.C+
 MUONTimeRawStreamTrigger(filename, maxEvent);
</pre>

where \em filename is the name of a file containing the raw data, or alternatively
the directory containing rawX (X being an integer) paths with the raw DDL
data. The \em maxEvent value is the maximum event to process (default set to
1000). Thus the macro will time the algorithm for all events in the range
[0 .. maxEvent-1].


\section raw_s3 Special flags for high performance tracker and trigger DDL decoders.

There are three flags that are available through the AliMUONRawStreamTrackerHP
interface, which effect the decoding of raw data in special ways.
\li TryRecover (default = false)
\li AutoDetectTrailer (default = true)
\li CheckForTrailer  (default = true)

The TryRecover flag is used to enable special logic in the decoder that tries
to recover from a partially corrupt raw data structure header, or a corrupt/missing
end of DDL trailer. Normally if the header is found to be corrupt (by seeing
that the block lengths do not correspond),
then the whole corresponding structure is skipped and we move on to the
next structure that we are sure of, or skip the rest of the DDL.
With the option to recover turned on, we attempt to figure out what the correct
values for the structure length should have been. This is possible because there
is redundant information in the headers giving the structure size. This recovery
procedure will work for single bit flips found inside the headers, but obviously
with more severe corruption, the header will not be recoverable.
The TryRecover flag is set to kFALSE by default to disable this feature, meaning
any structures in the raw data that cause problems in any way are skipped by
default. More details about the corrupt header recovery procedure can be found
in the method documentation AliMUONTrackerDDLDecoder::TryRecoverStruct().

The AutoDetectTrailer and CheckForTrailer flags are used to deal with the
difference in raw data format from the real tracker chambers and the simulated
data generated with AliRoot up to March 2008. The end of DDL trailer words are
missing in the simulated raw data in the older versions of AliRoot. To tell the
decoder that it should expect the trailer words as generated by the real readout
electronics then the CheckForTrailer flag should be set to kTRUE (this is the
default). If the trailer words are missing in the raw data then CheckForTrailer
should be kFALSE. Alternatively one can set the flag AutoDetectTrailer to kTRUE
(which is the default) so that the decoder will try to autodetect if the trailer
words are there or not. When AutoDetectTrailer is set to kTRUE then the
CheckForTrailer flag is ignored. Only when AutoDetectTrailer is set to kFALSE
will the CheckForTrailer flag have an effect.

Each of these decoder flags can be set to kTRUE or kFALSE, or tested by the
corresponding setter and getter methods respectively. For example, to set the
TryRecover flag to kTRUE and fetch its value afterwards, use the following code:
\code
AliMUONRawStreamTrackerHP decoder;
// to set the flag:
decoder.TryRecover(kTRUE);
// or to get the value of the flag:
Bool_t value = decoder.TryRecover();
\endcode
Similarly, the other flags are manipulated with corresponding methods having
the same name as the flag.

The trigger DDL decoder AliMUONRawStreamTriggerHP has the following flag available:
\li TryRecover (default = false)

For the AliMUONRawStreamTriggerHP decoder the TryRecover flag can be set in the
same way as for AliMUONRawStreamTrackerHP, with a call to the TryRecover() method.
For trigger DDLs this option will enable logic, which attempts to find the next
correct header / structure marker key in the DDL stream, whenever such a marker
has been found corrupt or missing. Decoding then continues from the new location
found or stops if no good key was found. The default setting is to disable this
logic, since it is only useful to try recover corrupted data.

\note Raw data containing software scalars (Start-of-Data events for example)
from the trigger detector taken during the Feb-March 2008 cosmics run is corrupt,
but can be successfully decoded by enabling this TryRecover flag for the
trigger DDL decoder.


This chapter is defined in the READMEraw.txt file.

*/

