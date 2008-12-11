// $Id$

/*! 

\page README_trigger Trigger


\section trigger_s1  How to reprocess trigger decision from already produced digits

The MUONTrigger.C macro can be used to check the trigger algorithm w/o 
having to (re-)perform simulation and digitalization. 
It loads the digits, erase TreeR and store the current trigger output in 
TreeR.

The different trigger outputs can be compared by looking at the GLT branch 
of TreeD (filled during simulation) and the TC branch of TreeR (filled from 
a copy of TreeD during reconstruction or with this macro).
Note: rec points from tracking chamber will be lost.

Usage:
<pre>
root [0] .L $ALICE_ROOT/MUON/MUONTrigger.C+
root [1] MUONTrigger("galice.root")
</pre>


\section trigger_s2 OFFLINE trigger Graphical User Interface (GUI) data quality and debugging tool

- read digits and local trigger decision from simulated/real data
- display the strips in local boards
- reprocess trigger decision inside AliRoot
- set x/y strips interactively on boards and transfer them to the AliRoot
  TriggerElectronics, execute trigger algorithm and recover the local trigger
  decision

Usage (aliroot):
<pre>
root [0] new AliMUONTriggerGUI
</pre>

Main window shows the map of the local boards as seen from the I.P.
The main window is position sensitive: by focusing a board, the "tip text"
shows the board name, the crate name, the board ID and the board
internal number in the GUI.

Menus:


\subsection trigger_s2_sub1 File

- "Run input" - open a file and start with a given event number:
              - "your_path/galice.root" to use simulated (or re-created) aliroot digits 
              - confirm with "Apply (galice)"              
              - "your_path/rawfilename.root" to use raw data in root format
              - confirm with "Apply (raw)"
- "Control"   - navigate in the tree with events
- "Exit"      - exit the main application


\subsection trigger_s2_sub2 Maps

- "Digits map"   - graphical map showing digits in the four chambers, MT11 ... MT22
- "Reset digits" - clean the digits map

\subsection trigger_s2_sub3 Chambers digit maps

- "Update" - update the map after:
             - loading of another event
             - changing interactively the strip signals in boards GUI


\subsection trigger_s2_sub4 Circuit

- "Open"   - open a board GUI by circuit (board) number


\subsection trigger_s2_sub5 TriggerDSET, run the trigger algorithm with interactively set strips 

- this is an interface to AliMUONTriggerElectronics
- "Digit store"    - create a digit store (object) with the current board digits (input from the circuit GUI)
- "Trigger store"  - create a trigger store (object) from the digit store, using AliMUONTriggerElectronics; 
each type of store can be "cleared" and "printed"
- "Front End Test" - simulate a FET (all strips fire) on all boards or on a regional board
- "Write raw data" - save the trigger store in a raw data file (DATE and ROOT) format

\subsection trigger_s2_sub6 Circuit GUI ("Circuit/Open" or click on boards map)

- the window title shows the name of the board, the number and the status
- "Draw" visualize x/y strips
- "Set/unset" x (or) y strips, on mouse click
- "Digits" create board digits from the actual configuration created in the GUI
- "Reset" reset modification on strips done interactively

The sequence to test the trigger algorithm is:
- open a board GUI
- set some x/y strips
- press "Digits" to create board digits
- from the TrigerDSET menu create a digit store, then a trigger store and print it to see the trigger
decision
- write raw data for further investigation

\section trigger_s3 How to check integrated trigger efficiency

The MUONTriggerEfficiency.C macro (included in the check scripts) calculates
the trigger efficiency for the 2 pt cuts. 
The output is stored in MUONTriggerEfficiency.out file.

Usage:
<pre>
root [0] .L $ALICE_ROOT/MUON/MUONTriggerEfficiency.C+
root [1] MUONTriggerEfficiency()
</pre>

For the CVS default version of the trigger LUT (i.e. lutAptLpt1Hpt1p7.root),
Two LUTs are stored in the CDB (/Calib/TriggerLut/)
Run0_999999999_v0_s0.root   with Lpt 1.0 GeV and Hpt 1.7 GeV
and
Run0_999999999_v1_s0.root   with Lpt 0.0 GeV and Hpt 1.0 GeV (default)

The reference for J/psi and Upsilon is as below
 For 1000 Jpsi events with:
<pre>
    AliGenParam *gener = new AliGenParam(1, AliGenMUONlib::kJpsi);
    gener->SetMomentumRange(0,999);
    gener->SetPtRange(0,100.);
    gener->SetPhiRange(0., 360.);
    gener->SetCutOnChild(1);
    gener->SetChildPhiRange(0.,360.);
    gener->SetChildThetaRange(171.0,178.0);
    gener->SetOrigin(0,0,0);          
    gener->SetForceDecay(kDiMuon);
    gener->SetTrackingFlag(1);
</pre>

 the output should be 
<pre>
  Efficiency Lpt cut = 0.9061 +/- 0.0456
  Efficiency Hpt cut = 0.6943 +/- 0.0376
</pre>

 Similarly, for 1000 Upsilon events, the output should be

<pre>
  Efficiency Lpt cut = 0.9872 +/- 0.0458
  Efficiency Hpt cut = 0.9851 +/- 0.0457
</pre>


\section trigger_s4 How to check single muon trigger efficiency versus pt

The MUONTriggerEfficiencyPt.C macro produces trigger single muon efficiency 
versus pt plots for the 2 pt cuts. 
Results are compared to the reference (red curves).   
To be used with (at least) 10000 events as follows
<pre>
   AliGenBox * gener = new AliGenBox(1);
   gener->SetPtRange(0.,10.);
   gener->SetPhiRange(0., 360.);         
   gener->SetThetaRange(171.000,178.001);
   gener->SetPart(13);           // or -13
   gener->SetOrigin(0.,0., 0.);  
   gener->SetSigma(0.0, 0.0, 0.0);     
</pre>

Outputs are stored in MUONTriggerEfficiencyPt.gif/eps/out files
Important note: this macro works with one (real) muon track per event only

Usage:
<pre>
root [0] .L $ALICE_ROOT/MUON/MUONTriggerEfficiencyPt.C+
root [1] MUONTriggerEfficiencyPt()
</pre>


\section trigger_s5 How to get trigger chamber efficiency from data

Trigger chamber efficiency map is calculated during reconstruction and saved in AliESDs.root
In order to view and save the map, use macro MUONTriggerChamberEfficiency.C

To compile MUONTriggerChamberEfficiency.C
<pre>
.includepath $ALICE_ROOT/MUON
.L $ALICE_ROOT/MUON/MUONTriggerChamberEfficiency.C+
</pre>

To run MUONTriggerChamberEfficiency.C
<pre>
MUONTriggerChamberEfficiency();
</pre>

If you want to make the calculated map available for next simulation use option kTRUE, i.e.
<pre>
MUONTriggerChamberEfficiency(kTRUE);
</pre>

When running next simulation, please remember to activate trigger efficiency
by adding in Config.C:
<pre>
MUON->SetTriggerEffCells(1);
</pre>

This chapter is defined in the READMEtrigger.txt file.

*/

