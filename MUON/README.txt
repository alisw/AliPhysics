// $Id$

/*! 

\page README_main README main
 
Please add  to this README file all information concerning 
general items about the MUON code or the libraries 
\ref base, \ref sim and \ref rec.

Other README pages of the muon code are
- \ref README_raw 
- \ref README_mapping 
- \ref README_calib 
- \ref README_geometry 
- \ref README_trigger 
- \ref README_shuttle 
- \ref README_evaluation 
- \ref README_fast

\section s1 How to check that your aliroot is working well

There is a script file AlirootRun_MUONtest.sh which 
allows for simulating, reconstructing and making the
invariant analysis of the generated Upsilon (1S).
The used configuration file is Config.C in MUON 
directory.

You have to type :
<pre>
$ALICE_ROOT/MUON/AlirootRun_MUONtest.sh [option]
</pre>

The complete list of the option is printed when you call
the script with whatever non valid option, .eg. h:

<pre>
./AlirootRun_MUONtest.sh h
ERROR : extra option not recognized
Usage: AlirootRun_MUONtest.sh options (-SRXsrxn:p:d:c:)
       -S (-s) perform (or not) simulation (default is do it, i.e -S)
       -R (-r) perform (or not) reconstruction (default is do it, i.e. -R)
       -X event (-x) perform (or not) checks and dumps (default is do it for event 5, i.e. -X 5)
       -n nevents (int) number of events to simulate (default 100)
       -p recoptions (quotified string) reconstruction options to use (default "SAVEDIGITS")
       -d full path to output directory (default /work/projects/alice/dev/AliRoot/MUON/test_out.100)
       -c full path to configuration file for simulation (default /work/projects/alice/dev/AliRoot/MUON/Config.C)
</pre>

The results of this test are saved in test_out.nevent directory.
Please note that the CDB (Condition DataBase) is now always *required* 
to perform either simulation or reconstruction. For the moment, a version
of that CDB is stored in CVS, so you should have one already in MUON/Calib
subdirectories.


\section s2 How to check that your aliroot is working VERY well

There is a script file AlirootRun_MUONlongtest.sh which
allows for simulating, reconstructing and making the
-+invariant analysis of the generated Upsilon (1S).
This script generates a large number of Upsilon (20k) 
in order to access differential quantities. 
The used configuration file is Config.C in MUON
directory.

One should really run this script to check if the MUON 
code can process a large number of events WITHOUT errors,
in particular before making important commits !!

You have to type :
<pre>
$ALICE_ROOT/MUON/AlirootRun_MUONtestlong.sh
</pre>
The results of this test are saved in testlong_out/ directory
and will be kept in CVS

(NOTE: the macros performing the calculations/plots MUONefficiency.C 
and MUONplotefficiency.C are also able to handle J/Psi if 
Config.C is modified accordingly )

\section s3  How to run a MUON generation

You only need to run the simulation part of the test script
AlirootRun_MUONtest.sh

\section s4 How to dump the content of Root data files 

To check the content of Root data files, the AliMUON*DataInterface classes
provides the functions to produce an ASCII output on the screen
which can be redirected on the file:

for MC information, use AliMUONMCDataInterface :

<pre>
> aliroot (or root with just the loading of MUON libs, see loadlibs.C)
root [0] AliMUONMCDataInterface mcdi("galice.root");
root [1] mcdi.DumpKine(5);       > dump.kine
root [2] mcdi.DumpHits(5);       > dump.hits
root [3] mcdi.DumpTrackRefs(5);  > dump.trackrefs
</pre>

for all other information, use AliMUONDataInterface :

<pre>
> aliroot
root [0] AliMUONDataInterface di("galice.root");
root [1] di.DumpDigits(5);     > dump.digits
root [2] di.DumpSDigits(5);    > dump.sdigits
root [3] di.DumpRecPoints(5);  > dump.recpoints
root [4] di.DumpTrigger(5); > dump.rectrigger
</pre>

Remind that during simulation and reconstruction two 
differents galice.root are generated: one for the generation 
(simulation) and other during the reconstruction.

If you open the wrong galice.root file you could get:
root [0] AliMUONMCDataInterface mcdi("galice.root");
root [1] mcdi.DumpKine(5);
W-AliRunLoader::GetEvent: Stack not found in header
E-TFile::TFile: file ./Kinematics.root does not exist

\section s5 Tracking parameters, cuts, energy loss and physics processes

Tracking parameters in MUON are automatically defined by GEANT
MUON takes the default values of CUTs  and physics processes
defined by the Config files, except for the gas mixture medium 
of the tracking chambers. The CUT's and physics processes of
the gas mixture medium  is then defined in the galice.cuts file
in the data directory. In particular ILOSS parameter MUST be
equal unity (1) in order simulate a realistic energy loss
distribution (mean value and fluctuations) in the active gas.

\section s6  Tracking of particle in the magnetic field

GEANT has two ways for tracking charged particles in the 
magnetic field: HELIX et RKUTA.
HELIX is faster and works well if the gradient of magnetic 
field is small. 
For MUON, HELIX is a not a good approximation and we must 
use RKUTA to get the optimal mass resolution of the 
spectrometer. The choice of HELIX or RKUTA is done in the
config file when the magnetic field is defined:
<pre>
  AliMagFMaps* field = new AliMagFMaps("Maps","Maps", TRACKING, FACTOR, MAXB, AliMagFMaps::k5kG);
  gAlice->SetField(field);
</pre>  
TRACKING must be 1 for RKUTA and 2 for HELIX (the default value for aliroot is 2 (HELIX))
FACTOR allows you to set the magnetic field to 0, just putting FACTOR=0. Default value is 1.
MAXB is the maximum magnetic field which is 10.T

\section s7 How to tune muon track reconstruction

Several options and adjustable parameters are available for both Kalman and Original
tracking algorithms (hard coded for the moment in AliMUONVTrackReconstructor.cxx):
- *fgkSigmaToCutForTracking* : quality cut used to select new clusters to be
  attached to the track candidate and to select good tracks.
- *fgkMakeTrackCandidatesFast* : if this flag is set to 'true', the track candidates
  are made assuming linear propagation between stations 4 and 5.
- *fgkTrackAllTracks* : according to the value of this flag, in case that several
  new clusters pass the quality cut, either we consider all the possibilities
  (duplicating tracks) or we attach only the best cluster.
- *fgkRecoverTracks* : if this flag is set to 'true', we try to recover the tracks
  lost during the tracking by removing the worst of the 2 clusters attached in the
  previous station.
- *fgkImproveTracks* : if this flag is set to 'true', we try to improve the quality
  of the tracks at the end of the tracking by removing clusters that do not pass
  new quality cut (within the limit that we must keep at least one cluster per
  the station).
- *fgkSigmaToCutForImprovement* : quality cut used when we try to improve the
  quality of the tracks.

\section s8 MUON cocktail for physics ..............

There is a MUON cocktail generator of the muon sources in the
EVGEN directory. This class derives from AliGenCocktail.
In the init of this class I have filled the cocktail with 
the muon sources: J/Psi, Upsilon, Open Charm, Open Beauty, 
Pion, Kaons. The code needs only the production cross section 
at 4pi (for the moment this values are in the code since I 
prefere them do not be modified), and the code calculates the  
rate of particles in the acceptance, making the scaling based 
on the number of collisions for the hard probes and on the  
number of participants for soft sources: Pions and Kaons.

In the Genereate of this class all entries in the cocktail 
are called and we define a "primordial trigger" with requires 
a minimum number of muons above a Pt cut in the required acceptance.
In order to normalized to the real number of simulated events, 
there are 2 data members in the class fNsuceeded adn fNGenerate 
which tell us what is the biais source.

Enclose an example to use this generator:   
<pre>
AliGenMUONCocktail * gener = new AliGenMUONCocktail();
gener->SetPtRange(1.,100.);       // Transverse momentum range  
gener->SetPhiRange(0.,360.);    // Azimuthal angle range 
gener->SetYRange(-4.0,-2.4);
gener->SetMuonPtCut(1.);
gener->SetMuonThetaCut(171.,178.);
gener->SetMuonMultiplicity(2);
gener->SetImpactParameterRange(0.,5.); // 10% most centra PbPb collisions
gener->SetVertexSmear(kPerTrack);  
gener->SetOrigin(0,0,0);        // Vertex position
gener->SetSigma(0,0,0.0);       // Sigma in (X,Y,Z) (cm) on IP position
gener->Init();
</pre>
 

\section s9 How to simulate events with misaligned geometry in local CDB

If you want to use a misaligned geometry to simulate some
events you can use a local CDB. For this need to follow
the next steps:

- Generate misaligned data in local CDB.
You can use MUONGenerateGeometryData.C as described above in
the corresponding section. Let's assume you used the default
residual misalignment settings, then you have a local CDB in
your working directory called ResMisAlignCDB containing
misalignement data (ResMisAlignCDB/MUON/Align).

- Tell AliSimulation you want to use your local CDB for 
MUON/Align/Data
To do this you need to instantiate the AliCDBManager, set the
default storage and set the specific storage for MUON/Align/Data,
before instantiating AliSimulation (see for example the commented
lines AlirootRun_MUONtest.sh).

<pre>
aliroot -b  >& testSim.out << EOF
AliCDBManager* man = AliCDBManager::Instance();
man->SetDefaultStorage("local://$ALICE_ROOT");
man->SetSpecificStorage("MUON/align/Data","local://ResMisAlignCDB");
AliSimulation MuonSim("$ALICE_ROOT/MUON/Config.C");
MuonSim.SetWriteRawData("MUON");
MuonSim.Run(10);
.q
EOF
</pre>

\section s10 How to Merge events

You can merge 2 types of simulated events. For example, 
you can simulate Hijing events, and then simulate muons
merging both.

Merging is done at the sdigits level, so Kinematics files 
of the merged events will just correspond to the 
Config.C simulated file).

You must, first, do the Hijing simulation and store it 
in directory $HIJING_SIM. Note that for merging you 
won't need Kinematics files of the Hijing simulation...

Hijing simulation

<pre>
aliroot -b << EOF
AliSimulation HijingSim("$HIJING_SIM/YourConfigForHIJING.C")
HijingSim.Run(5)
.q
EOF
</pre>

You cand build YourConfigFroHIJING.C File from the 
ConfigPPR file in AliRoot/macros module.

Then you can do muon simulation and reconstruction
merging both simulated events. In next example, we are
merging 20 times each Hijing event in order to simulate 
100 muons merged with 5 Hijing events.

<pre>
aliroot -b << EOF
AliSimulation MuonSim("$ALICE_ROOT/MUON/Config.C")
MuonSim.MergeWith("$HIJING_SIM/galice.root",20) //parameters are the Hijing simulation file and the number of times we use each Hijing event
MuonSim.Run(100) // number of muon (Config.C) events
.q
EOF


aliroot -b << EOF
TPluginManager * pluginmanager = gROOT->GetPluginManager()
pluginmanager->AddHandler("AliReconstructor","MUON","AliMUONReconstructor","MUON","AliMUONReconstructor()")
AliReconstruction  MuonRec("galice.root")
MuonRec.SetRunTracking("")
MuonRec.SetRunVertexFinder(kFALSE)
MuonRec.SetRunLocalReconstruction("MUON")
MuonRec.SetFillESD("MUON")
MuonRec.Run()
.q
EOF
</pre>

\section s11 ...On track numbering 

All generated particles, including primary and secondary
particles are put on the stack. The secondary particles are kept
in the stack only if they gave a hit in *any* of the ALICE detectors
The number of all particles placed on the stack for a given event 
can be obtained with
Int_t nPart = AliStack::GetNtrack();
Looping from 0 to nPart via AliStack::Particle(ipart)
gives the particle listing as obtained from the particle generator (primaries) 
and Monte Carlo (secondaries).

The particle response in the detector, a hit, is registered
in the hits tree and the hits are filled with each primary track.
The total number of "tracks" (fills of the tree) can be obtained
with ntracks = AliMUONMCDataInterface::NumberOfTracks(event) and is usually smaller than "nPart".
Since particles can also deposit hits in other detectors than 
the MUON spectrometer, there will be many "tracks" (fills) in the hit-tree
without a hit in MUON.

The correspondence between "track ID" in the hits-tree ("itr") and the
particle ID for particles on the stack (i.e. generated particles) can be
obtained via:
<pre>
for (Int_t itr = 0; itr < ntracks; itr++) {
    AliMUONVHitStore* hitStore = mcDataInterface.HitStore(event,itr);
    //track "itr" of the hits-tree
    Int_t nhitstot = hitStore->GetSize();
    AliMUONHit* mHit; 
    TIter next(hitStore->CreateIterator());
    while ( ( mHit = static_cast<AliMUONHit*>(next()) ) )
    {   
       Int_t id = mHit->Track(); //gives particle ID on stack
       TParticle* particle = mcDataInterface.Stack(event)->Particle(id);
    }  
}
</pre>

where mcDataInterface has been obtained by
AliMUONMCDataInterface mcDataInterface("galice.root");

During the procedure to go from hits to digits, the hits 
are summed up such that more than one track can contribute
to a given digit. As a consequence the method
Int_t AliMUONDigit::Track(Int_t trackID)
takes an argument, where "trackID" runs from 0 to 
AliMUONDigit::Ntracks() to provide the reference to *all*
tracks that contributed to it. The returned track ID is the one 
referred to in the hit-tree. To know which is the generated particle
that deposited a given digit one has to follow the sequence of the kind:
(shown here using the simple, but not fast, DataInterface interfaces) :

<pre>
AliMUONMCDataInterface mcdi("galice.root");
AliMUONDataInterface di("galice.root");

AliMUONVDigitStore* digitStore = di.DigitStore(event);
AliMUONVDigit* mDigit = ... get some digit from the digitStore

for (int tr = 0; tr < mDigit->Ntracks(); tr++)
{
   Int_t hitTrackID = mDigit->Track(tr);
   // get the hits corresponding to this trackID
   AliMUONHitStore* hitStore = mcdi.HitStore(event,hitTrackID);
   // loop over the hits
   TIter hNext(hitStore->CreateIterator());
   AliMUONHit* mHit;
   while ( ( mHit = static_cast<AliMUONHit*>(hNext()) ) )
   {
    Int_t numPart = mHit->Track(); //gives ID of particle on the stack
    Int_t idTrack = mHit->Particle(); //gives flavour code of the particle
   }
}
</pre>


\section s12 How to process invariant mass spectra for J/psi or Upsilon

The macro MUONmassPlot_ESD.C reads back the MUON ESD informations and compute 
the invariant mass spectra and corresponding uncorelated background. 
Moreover gives the number of event in the resonance peak and the number of triggers.
<pre>
Usage:
root [0] .L $ALICE_ROOT/MUON/MUONmassPlot_ESD.C+
root [1] MUONmassPlot_ESD(ExtrapToVertex, 
			geoFilenam, filename
                        FirstEvent, LastEvent, 
	                esdFileName,
			ResType, Chi2Cut,
			PtCutMin, PtCutMax,
        		massMin, massMax)

with:
ExtrapToVertex (default -1)
      <0: no extrapolation;
      =0: extrapolation to (0,0,0);
      >0: extrapolation to ESDVertex if available, else to (0,0,0)
geoFilename (default "geometry.root") geometry file name needed to extrap to vertex
filename    (default "galice.root") galice root file name
FirstEvent  (default 0)
LastEvent   (default 10000)
esdFileName (default "AliESDs.root") esd root file name
ResType     (default 553):   553 for Upsilon, anything else for J/Psi
Chi2Cut     (default 100):   keep only tracks with chi2 per d.o.f. < Chi2Cut
PtCutMin    (default 1):     keep only tracks with transverse momentum > PtCutMin
PtCutMax    (default 10000): keep only tracks with transverse momentum < PtCutMax
massMin     (default 9.17 for Upsilon) keep only invariant masses with 
massMax     (default 9.77 for Upsilon) massMin < mass < massMax
</pre>

\section s13  Still working ..............

*/
