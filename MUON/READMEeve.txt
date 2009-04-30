// $Id$

/*!

\page README_eve Event display

\section eve_s1 The EVE event display

This README file shows how to use the EVE event display for the MUON detector. 
The explained procedure is a complement to the visualisation of the full ALICE
detector discussed on the ALICE Offline page, "Visualisation": 
http://aliceinfo.cern.ch/Offline/Activities/Visualisation/

\subsection eve_s1_sub1 Macros for the event display

Found in $ALICE_ROOT/EVE/alice-macros and $ALICE_ROOT/EVE/macros.
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

- launch the "alieve" executable

<pre>
   alieve
</pre>

- load the following macros (.L macro.C)

<pre>
   alieve_init.C
   event_goto.C
   MUON_displayData.C
   or
   MUON_displaySimu.C
   geom_gentle_muon.C
</pre>

- initialize the EVE environment with your data

<pre>
   alieve_init("directory_to_data",event_number)
</pre>

- display the current event

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

- draw the "gentle" geometry (sensitive volumes)

   geom_gentle_muon()

- navigate in the event list

   event_goto(n)

- Shift + right mouse button one a muon track opens the context menu (track 
and trigger information is available)

- a new flag argument is added to the MUON_display macros for displaying the 
the clusters taken from the ESD (default), only those attached to the tracks, 
or from MUONRecPoints, all reconstructed clusters

<pre>
MUON_displaySimu(x,x,1) - clusters from ESD
MUON_displaySimu(x,x,0) - clusters from rec points
MUON_displayData(x,x,1) - clusters from ESD
MUON_displayData(x,x,0) - clusters from rec points
</pre>

This chapter is defined in the READMEeve.txt file.

*/

