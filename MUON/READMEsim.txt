// $Id$

/*! \page README_sim MUON Simulation

The simulation encompasses the following tasks :

- Generation of MC particles (the kinematics of the event ends up in the TreeK
 of Kinematics#.root)                              
  
- Tracking particles through the detector using 
the Virtual Monte Carlo, producing AliMUONHit objects, that end up in 
 the TreeH of MUON.Hits#.root file(s). This part is steered by AliMUON and its child
AliMUONv1 classes.

- Converting MC hits into AliMUONVDigit, called SDigits, that end up in the TreeS
 of the MUON.SDigits#.root file(s). A S(ummable)Digit is a pad with its associated
charge, but no noise or electronics response function applied. Steered by AliMUONSDigitizerV2 class.

- Converting SDigits into Digits, by applying electronics calibrations. Actually, we de-calibrate
 the digits at this stage, by adding a pedestal and dividing by a gain, more or less. Steered
 by AliMUONDigitizerV3 class. Digits end up in TreeD of MUON.Digits#.root file(s). In addition,
 for the trigger, we create AliMUONLocalTrigger, AliMUONRegionalTrigger and AliMUONGlobalTrigger objects 
at this stage, that ends up in TreeD as well.

- Convert the Digits into RAW data, in a format that should be exactly the same as real data from the
DAQ. Performed by AliMUONRawWriter.

From there on, the reconstruction can proceed, in the very same way for real or simulated data,
 as long as they are in RAW format.

\section sim_s1  How to run a MUON generation

You only need to run the simulation part of the test script
AlirootRun_MUONtest.sh


\section sim_s2 Tracking parameters, cuts, energy loss and physics processes

Tracking parameters in MUON are automatically defined by GEANT
MUON takes the default values of CUTs  and physics processes
defined by the Config files, except for the gas mixture medium 
of the tracking chambers. The CUT's and physics processes of
the gas mixture medium  is then defined in the galice.cuts file
in the data directory. In particular ILOSS parameter MUST be
equal unity (1) in order simulate a realistic energy loss
distribution (mean value and fluctuations) in the active gas.

\section sim_s3  Tracking of particle in the magnetic field

GEANT has two integration methods for tracking charged particles in the 
magnetic field: HELIX et RKUTA.
HELIX is faster and works well if the gradient of magnetic 
field is small. 
For MUON, HELIX is a not a good approximation and we must 
use RKUTA to get the optimal mass resolution of the 
spectrometer. The choice of HELIX or RKUTA is done in the
config file when the magnetic field is defined:
<pre>
  TGeoGlobalMagField::Instance()
    ->SetField(new AliMagF("Maps","Maps", INTEG, FACTOR_SOL, FACTOR_DIP, MAXB, AliMagF::k5kG));  
</pre>  
INTEG must be 1 for RKUTA and 2 for HELIX (the default value for aliroot is 2 (HELIX)).
FACTOR_SOL, FACTOR_DIP allow you to set the multiplicative factor for solenoid
and dipole, respectively; just putting FACTOR_SOL=0 or FACTOR_DIP=0 will set the
magnetic field for solenoid or dipole to 0. Default values are 1.0, 1.0. 
MAXB is the maximum magnetic field which default value is 10.T

\section sim_s4 Tailing effect

The control to turn on/off the parametrized tailing effect: 
<pre>
AliMUON::SetTailEffect(Bool_t),
</pre>

The parameter to tune increase/decrease the tailing effect is kept inside,
AliMUONResponseV0::DisIntegrate(). This parameter is an integer number 
(excluding zero and four), the higher the value is the less is the tailing 
effect:
<pre>
Int_t para = 5; 
</pre>
Zero is excluded because it gives straight line transformation, and four
is excluded because the AliRoot simulation chain spends VERY VERY long
time in AliMUONResponseV0::DisIntegrate method, which reason was not yet
understood. The parameter for 1, 2, 3, 5, 6, 7, 8, 9, 10 were checked with 
no slowing down problem, however parameters greater than 6 give almost no 
tailing effect since they basically correspond to higher order polynomial 
transform.


\section sim_s5 MUON cocktail generator

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
 
\section sim_s6 How to simulate events with misaligned geometry in local CDB

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
man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
man->SetSpecificStorage("MUON/align/Data","local://ResMisAlignCDB");
AliSimulation MuonSim("$ALICE_ROOT/MUON/Config.C");
MuonSim.SetWriteRawData("MUON");
MuonSim.Run(10);
.q
EOF
</pre>

\section sim_s7 How to Merge events

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

\section sim_s8 On track numbering 

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

This chapter is defined in the READMEsim.txt file.

*/
