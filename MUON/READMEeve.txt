// $Id$

/*!

\page README_eve MUON Event display

This README file shows how to use both the old and new EVE event display for
the MUON detector. The explained procedure is a complement to the visualisation
of the full ALICE detector discussed on the ALICE Offline page, "Visualisation": 
http://aliceinfo.cern.ch/Offline/Activities/Visualisation/

\section eve_s1 The old EVE event display

\subsection eve_s1_sub1 Macros

Found in $ALICE_ROOT/EVE/macros and $ALICE_ROOT/EVE/macros.
- MUON_displaySimu.C , to be used with simulations (shows also reference 
tracks, Monte-Carlo tracks and hits)
- MUON_displayData.C , to be used with reconstructed raw data (with the option
SAVEDIGITS in the reconstruction, "simulation" digits are also produced and 
can be visualized)
- alieve_init.C , initialisation of the EVE environment
- event_goto.C , navigation in the tree of events
- geom_gentle_muon.C , draws a "light" geometry of the MUON detector obtained 
from the sensitive volumes of all MUON chambers; the default geometry of the
event display is made of contours of the chambers extracted from the mapping

\subsection eve_s1_sub2 Usage

-# Launch the "alieve" executable
<pre>
alieve
</pre>
-# Load the following macros (.L macro.C)
<pre>
alieve_init.C
event_goto.C
MUON_displayData.C
or
MUON_displaySimu.C
geom_gentle_muon.C
</pre>
-# Initialize the EVE environment with your data
<pre>
alieve_init("local://$ALICE_ROOT/OCDB","directory_to_data",event_number)
</pre>
-# Display the current event
<pre>
MUON_displaySimu(0,0)   - do not show tracks
MUON_displaySimu(0,1)   - show tracks
MUON_displaySimu(1,0)   - do not tracks, read digits from produced raw data
MUON_displaySimu(1,1)   - show tracks, read digits from produced raw data
MUON_displayData(1,0)   - do not show tracks
MUON_displayData(1,1)   - show tracks
MUON_displayData(0,0)   - do not show tracks, read the saved digits
MUON_displayData(0,1)   - show tracks, read the saved digits
</pre>
-# Draw the "gentle" geometry (sensitive volumes)
<pre>
geom_gentle_muon()
</pre>
-# Navigate in the event list
<pre>
event_goto(n)
</pre>
-# Shift + right mouse button one a muon track opens the context menu (track 
and trigger information is available)
<pre></pre>
-# A new flag argument is added to the MUON_display macros for displaying the 
the clusters taken from the ESD (default), only those attached to the tracks, 
or from MUONRecPoints, all reconstructed clusters
<pre>
MUON_displaySimu(x,x,1) - clusters from ESD
MUON_displaySimu(x,x,0) - clusters from rec points
MUON_displayData(x,x,1) - clusters from ESD
MUON_displayData(x,x,0) - clusters from rec points
</pre>

\section eve_s2 The new EVE event display

\subsection eve_s1_sub3 Macros

Found in $ALICE_ROOT/EVE/macros and $ALICE_ROOT/EVE/macros.
- muon_init.C , to launch the new display (based on visscan_init.C)
- muon_raw.C , display digits from raw data
- muon_digits.C , display digits from MUON.Digits.root file
- muon_clusters.C , display clusters from MUON.RecPoints.root file
- esd_muon_tracks.C , display tracks, clusters attached to tracks and digits
attached to clusters (if any) from ESD
- muon_trackRefs.C , display simulated tracks and hits in the MUON chambers
- kine_tracks.C , display all the simulated tracks
- geom_gentle_muon.C , draws a "light" geometry of the MUON detector obtained 
from the sensitive volumes of all MUON chambers; the default geometry of the
event display is made of contours of the chambers extracted from the mapping

\subsection eve_s1_sub4 Usage

-# Launch the "alieve" executable
<pre>
alieve
</pre>
-# Initialize the EVE environment with the path to the OCDB and your data.
   -# Reading data locally:
<pre>
.x muon_init.C("local://$ALICE_ROOT/OCDB","directory_to_data")
</pre>
   -# Reading data on the grid (link to corresponding raw data is made automatically):
<pre>
.x muon_init.C("raw://","alien:///alice/data/2009/LHC09c/000084039/ESDs/pass1/09000084039008.10")
</pre>
-# Change what to draw: \n 
Just enable/disable the corresponding macros in the tab "DataSelection"
by clicking on the macro name and checking the box "Active"
<pre></pre>
-# Navigate in the event list: \n
Just use the buttons or choose the event number in the bottom panel. You can also scan the events by checking the box "Autoload" and change the time between 2 drawing.
<pre></pre>
-# Leave the cursor on a track to pop-up its characteristics

\section eve_s3 The simplified "gentle" geometry for the event display

Execute from aliroot the macro MUONGenerateGentleGeometry.C and place the resulting file in EVE/alice-data. This file is used by the macro geom_gentle_muon.C
from EVE/macros.

This chapter is defined in the READMEeve.txt file.

*/

