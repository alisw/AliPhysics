// $Id$

/*! 

\page README_trigger README trigger


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


\section trigger_s2 OFFLINE trigger GUI data quality and debugging tool

- read digits and local trigger decision from simulated/real data
- display
- reprocess trigger decision inside AliRoot
- set x/y strips interactively on boards and transfer them to the AliRoot
  TriggerElectronics, execute trigger algorithm and recover the local trigger
  decision

To run:

<pre>
aliroot
new AliMUONTriggerGUI()
</pre>

Main window shows the map of the local boards as seen from the I.P.
The main window is position sensitive (after file initialization) and it is 
possible to open a GUI for a circuit.

By menus:


\subsection trigger_s2_sub1 File

- Run     - open a file and start with a given event number
          takes the full path <path>/galice.root
- Control - navigate in the tree with events
- Exit    - exit the main application


\subsection trigger_s2_sub2 Maps

- Digits map   - graphical map with digits in the four chambers, MT11 ... MT22

\subsection trigger_s2_sub3 Chambers digit maps window

- Update - update the map after:
             - loading of another event
             - changing interactively the strip signals in boards GUI


\subsection trigger_s2_sub4 Circuit

- Open   - open a board GUI by circuit number


\subsection trigger_s2_sub5 Trigger

- Trigger DSET  - (re)run the trigger algorithm


\subsection trigger_s2_sub6 Circuit GUI window

- visualize x/y strips
- "Set/unset" x (or) y strips
- "Digits" create board digits from the actual configuration created in the GUI
- "Reset" reset modification on strips done interactively


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
  Efficiency Lpt cut = 0.7362 +/- 0.0391
  Efficiency Hpt cut = 0.2662 +/- 0.0201
</pre>

 Similarly, for 1000 Upsilon events, the output should be

<pre>
  Efficiency Lpt cut = 0.9806 +/- 0.0457
  Efficiency Hpt cut = 0.9537 +/- 0.0448
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

*/

