// $Id$

/*! \page README_cosmics MUON The Software testing on Cosmics Test Data


Seeveral macros has been developed during the Offline shifts
during the Feb-March 2008 cosmics runs. On this page we summarize
how to use these macros.

\section cosmics_s1 Software installation    

- For installing AliRoot see ALICE Offline Installation Web page at
http://aliceinfo.cern.ch/Offline/AliRoot/Installation.html

- and for installing alien see ALICE Offline Tutorial, slides 131 - 145 at
http://aliceinfo.cern.ch/Offline/Analysis/Tutorial

If we need a fast fix, it may be provided via a patch.txt file,
which should be then applied in this way: 
<pre> $&gt; cd $ALICE_ROOT
 $&gt; patch -p0 --posix &lt; patch.txt
 $&gt; make
</pre>

Eventually, we may get in situation when an important update of the MUON code 
is not yet included in the release, and we may then need to replace the MUON
from the release with the MUON from the trunk:
<pre> $&gt; cd $ALICE_ROOT
 $&gt; mv MUON MUON.release
 $&gt; svn co http://alisoft.cern.ch/AliRoot/trunk/MUON MUON
 $&gt; rm -fr $ALICE_ROOT/lib/tgt_$ALICE_TARGET/libMUON*
 $&gt; make
</pre>

\section cosmics_s2 Running reconstruction   

During the offline shifts, there has been added a a new macro 
runDataReconstruction.C. You have first to edit the macro and change
the value of <em>input</em> with the path to raw data file
which you want to reconstruct. Then you prepare an output directory
and run reconstruction in this way:

<pre> $&gt; alien-token-int
 $&gt; . /tmp/gclient_env_$UID
 $&gt; cd $ALICE_ROOT/MUON
 $&gt; mkdir runXYZ_c1_out
 $&gt; cp rootlogon.C .rootrc runDataReconstruction.C runXYZ_c1_out
 $&gt; cd runXYZ_c1_out
 $&gt; root
 root[0] .x runDataReconstruction.C(1) &gt;&amp; runReco.out
</pre>

This will run reconstruction with the first calibration option
selected ("NOGAIN"); to run the same with the second calibration option
("GAINCONSTANTCAPA"):
<pre> $&gt; cd $ALICE_ROOT/MUON
 $&gt; mkdir runXYZ_c2_out
 $&gt; cp rootlogon.C .rootrc runDataReconstruction.C runXYZ_c2_out
 $&gt; cd runXYZ_c2_out
 $&gt; root
 root[0] .x runDataReconstruction.C(2) &gt;&amp; runReco.out
 </pre>

\section cosmics_s3 Inspecting data with mchview

The new macro MUONOfflineShift.C will process the data
and generate the Root output file which can be then open with
the \em mchview program:

<pre> $&gt; alien-token-int
 $&gt; . /tmp/gclient_env_$UID
 $&gt; cd $ALICE_ROOT/MUON
 root[0] .L MUONOfflineShift.C+
 root[1] MUONOfflineShift("path_to_raw_file","basename of output file");
 .q
 
 $&gt; mchview --use basename.root
</pre>
    

\section cosmics_s4 Inspecting rec points

There has been added a new macro TestRecPoints.C.

The analysis of the trigger part needs only RecPoints: the digits are then re-created 
out of the local trigger information. The analysis of the tracker part needs the 
MUON.Digits, which can be created during reconstruction with the "SAVEDIGITS" option. 
(This option is switched on by default in the runDataReconstruction.C macro.)
The macro, then, performs the clusterization on the fly and analyses data.
This was necessary since, normally, the clusters are saved in the AliESDs.root 
but only for reconstructed tracks (which I guess will be very few in the cosmic run). 
Re-performing clusterization is surely more time expensive, but (I guess) is the 
only way to get information on all clusters even when no track is created.

To use the macro: 
<pre> $&gt; cd $ALICE_ROOT/MUON
 $&gt; root
 root[0] .L TestRecPoints.C+
 root[1] TestRecPoints(runXYZ_c1);
</pre>

It is also possible to check only tracker:
<pre> root[1] TestRecPoints("pathToData","outputDirectory",kOnlyTracker);
</pre>

or only the trigger
<pre> root[0] TestRecPoints("pathToData","outputDirectory",kOnlyTrigger);
</pre>
    
\section cosmics_s5 Event display

See EVE/README_MUON. 

This chapter is defined in the READMEcosmics.txt file.

*/
