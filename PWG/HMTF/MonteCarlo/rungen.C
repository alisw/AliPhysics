#if !defined (__CINT__) || (defined(__MAKECINT__))
#include <iostream>
#include "TStopwatch.h"
#include "TSystem.h"
#include "TString.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "TROOT.h"
#include "AliPDG.h"
#include "TRandom.h"
#include "AliRunLoader.h"
#include "AliSimulation.h"
#include "AliGenerator.h"
//#include "AliGenPythiaPlus.h"
#include "AliGenPythia.h"
//#include "AliPythia8.h"
#include "AliGenDPMjet.h"
#include "AliHeader.h"
#include "AliGenEventHeader.h"
#include "AliStack.h"
#include "AliRun.h"
#include "TTree.h"
#include "AliMC.h"
#include "AliGenReaderHepMC.h"
#include "AliGenExtFile.h"
#include "AliGenEposReader.h"
#include "AliGenEpos3EventHeader.h"

#endif

// This macro produces generator-level Monte Carlo for different
// generator/tunes.  It can downscale high multiplicity events, based
// on the thresholds and scaling factor defined in the multthresholds
// and scaling strings (these are comma separated list of values).
// Events above the highest threshold are not downscaled.

TString defaultMultThresholds = "10,20,30";
TString defaultScaling        = "1,1,1";
TString hepMCFile             = "";

const char * fEposFilename1 = "epos1.root";
const char * fEposFilename2 = "epos2.root";
const char * flagFile1 = "flag1";
const char * flagFile2 = "flag2";
const char * finishEpos = "finishepos";

typedef enum {kPhojet = -1, kHepMC = -2, kEpos3111 = -3111, kPyTuneMonash2013=14, kPyTuneCDFA=100,kPyTuneAtlasCSC=306, kPyTuneCMS6D6T=109, kPyTunePerugia0=320,  kPyTunePerugia2011 = 350} Tune_t;
void fastGen(Tune_t tune = kPyTuneCDFA , Float_t energy=7000, Int_t nev = 1,const Int_t nMultClasses=2, const Int_t* multThresholds=0, const Int_t* scaling=0, Int_t maxEposevents=100);

AliGenerator*  CreateGenerator(Tune_t tune, Float_t energy);

void rungen(Tune_t tune = kPyTuneCDFA,
            Float_t energy=7000,
            Int_t nev=10,
            const TString sMultThresholds=defaultMultThresholds.Data(), //array of thresholds in multiplicity
            const TString sScaling=defaultScaling.Data(),		//array of downscaling factors
            const char * hepmcfileloc = "dummy.hepmc",
	    Int_t maxEposevents=100 					//number of events per epos file
           )
{
  // Simulation and reconstruction
  TStopwatch timer;
  timer.Start();
  hepMCFile = hepmcfileloc;
  // Process the thresholds and scaling parameters, parsing the strings into int arrays
  Int_t* locMultThresholds;
  Int_t* locScaling;
  // Tokenize strings
  TObjArray* sobjMultThresholds=sMultThresholds.Tokenize(",");
  TObjArray* sobjScaling=sScaling.Tokenize(",");
  // Make sure scalings and thresholds have the same size
  if (sobjMultThresholds->GetEntries()!=sobjScaling->GetEntries()) {
    std::cout<<"Error: Number of scaling factors does not correspond to the number of thresholds."<<std::endl;
    return;
  }
  // Allocte arrays
  const Int_t nMultClasses=sobjMultThresholds->GetEntries(); //number of downscaled classes in multiplicity
  locMultThresholds = new Int_t[nMultClasses]; //array of thresholds in multiplicity
  locScaling        = new Int_t[nMultClasses]; //array of downscaling factors
  for (Int_t i=0; i<sobjMultThresholds->GetEntries(); i++) {
   locMultThresholds[i] =(((TObjString *)(sobjMultThresholds ->At(i)))->String()).Atoi();
   locScaling[i]        =(((TObjString *)(sobjScaling        ->At(i)))->String()).Atoi();
  }
  // Delete arrays of tokens
  delete sobjMultThresholds;
  delete sobjScaling;
  // Generate!
  fastGen(tune, energy, nev, nMultClasses, locMultThresholds, locScaling, maxEposevents);

  timer.Stop();
  timer.Print();
}

void fastGen(Tune_t tune  , Float_t energy, Int_t nev , const Int_t nMultClasses, const Int_t* multThresholds, const Int_t* scaling, Int_t maxEposevents)
{
  // Add all particles to the PDG database
  AliPDG::AddParticlesToPdgDataBase();

  // set the random seed
  TDatime date;

  // Create alirun object. Needed to run without aliroot, with plain
  // root. This was done to avoid any conflict with the libraries
  // loaded authomatically by aliroot (e.g. pythia).
  new AliRun("gAlice","The ALICE Off-line Simulation Framework");

  //  Runloader  
  AliRunLoader* rl = AliRunLoader::Open("galice.root", "FASTRUN","recreate");

  TString tmp(gSystem->Getenv("ALIEN_PROC_ID"));
  UInt_t seed =  tmp.Atoi();
  if ( seed == 0){
    seed = time(0);
  }

  gRandom->SetSeed(seed);
  std::cout<<"Seed for random number generation= "<<seed << "(ALIEN_PROC_ID="<<tmp.Data()<<")"<<std::endl; 
    
  rl->SetCompressionLevel(2);
  // The number of events per file has a huge impact on the speed of
  // the simulation. If you want to reuse the Kinematics files for a
  // full simulation, it is advisable to limit this to 200. In any
  // case, numbers > 1000 have a large impact on the processing speed.
  rl->SetNumberOfEventsPerFile(500);
  //  rl->SetNumberOfEventsPerFile(200);
  rl->LoadKinematics("RECREATE");
  rl->MakeTree("E");
  gAlice->SetRunLoader(rl);

  //  Create stack
  rl->MakeStack();
  AliStack* stack      = rl->Stack();
 
  //  Header
  AliHeader* header = rl->GetHeader();
  //
  Int_t iev;

  Int_t multArray[nMultClasses]; //array of counters of rejected events for different multiplicity classes
  for (Int_t i=0; i<nMultClasses; i++) {
    multArray[i]=0;
    std::cout<<"Thresholds and scaling factors: nMultClasses="<<nMultClasses<<"; multArray["<<i<<"]="<<multArray[i]<<"; multThresholds["<<i<<"]="<<multThresholds[i]<<"; scaling["<<i<<"]="<<scaling[i]<<std::endl;
  }

  Int_t flagGenerate=1; //flag for Epos3 processing
  Int_t Eposevents=0; //counter of processed Epos events (non-downscaled)
  //  Create and Initialize Generator
  if (tune==kEpos3111) {
    flagGenerate=0;
    while (!flagGenerate) {
      if (!(gSystem->AccessPathName(flagFile1))) flagGenerate=1;
      else sleep(1);
    }
  }

  AliGenerator *gener = CreateGenerator(tune,energy);

  gener->Init();
  // if nsd switch off single diffraction
  gener->SetStack(stack);
    
  //
  //                        Event Loop
  //

  for (iev = 0; iev < nev; iev++) { // Warning: In case an event is skipped because of the
                                    // downscaling thresholds below the iev counter is decremented 
	
    //  Initialize event
    header->Reset(0,iev);
    rl->SetEventNumber(iev);
    stack->Reset();
    rl->MakeTree("K");
    //	stack->ConnectTree();
    while (!flagGenerate) {     //  this loop is activated only for tune==kEpos3111
      if ( !(gSystem->AccessPathName(flagFile1)) ) {
        flagGenerate=1;
        ( (AliGenEposReader*) ((AliGenExtFile*)gener)->Reader() )->ChangeFile(fEposFilename1);
      }
      else if ( !(gSystem->AccessPathName(flagFile2)) ) {
        flagGenerate=2;
        ( (AliGenEposReader*) ((AliGenExtFile*)gener)->Reader() )->ChangeFile(fEposFilename2);
      }
      else sleep(1);
    }
     
    //  Generate event
    gener->Generate();

    if (tune==kEpos3111) {
      Eposevents++;
    //  cout<<flagGenerate<<" "<<Eposevents<<endl;
      if (Eposevents==maxEposevents) {
        Eposevents=0;
        if (flagGenerate==1) { gSystem->Unlink(flagFile1); gSystem->Unlink(fEposFilename1); flagGenerate=0; }
        else if (flagGenerate==2) { gSystem->Unlink(flagFile2); gSystem->Unlink(fEposFilename2); flagGenerate=0; }
      }
    }
    //  Analysis
    //    Int_t npart = stack->GetNprimary();
    Int_t npart = stack->GetNtransported();
    Int_t thisClass = -1; //the the multiplicity class number for this event
    Float_t evScaling = 1; // default scaling

    for (Int_t i=0; i<nMultClasses; i++) {
      if (npart<multThresholds[i]) {
        thisClass=i;
        break;
      }
    }
    
    if (thisClass>=0) {
      multArray[thisClass]++; //counter of rejected events
      // we have to decrement the event counter here, so that at the end of the day we get the nev events which was requested
      if (multArray[thisClass]<scaling[thisClass]) {iev--; continue;} 
      // we keep this one! Let's reset the counters
      if (multArray[thisClass]==scaling[thisClass]) multArray[thisClass]=0;
      // Set the event scaling to a non-default value
      evScaling = scaling[thisClass];
    }
    //    std::cout << "Ev " << iev << ", Mult:"<< npart <<", Class: " << thisClass << ", Scaling: " << evScaling << std::endl;

    if(!(iev%500)) printf("\n \n Event number %d \n \n", iev);
	
    //  Finish event
    header->SetNprimary(stack->GetNprimary());

    // Set weight of event
    // (taking into account what might already come from the generator!)
    AliGenEventHeader* eventheader = (AliGenEventHeader*) header->GenEventHeader();
    eventheader->SetEventWeight(evScaling * eventheader->EventWeight());
//    eventheader->SetNtransported(stack->GetNtransported());      
    header->SetNtrack(stack->GetNtrack());  
    //      I/O
    //	
    stack->FinishEvent();
    header->SetStack(stack);
    rl->TreeE()->Fill();
    rl->WriteKinematics("OVERWRITE"); // Can I move this (or not call it for every event?)

  } // event loop

  if (tune==kEpos3111) gSystem->Exec(Form("touch %s",finishEpos));
    //
    //                         Termination
    //  Generator
  gener->FinishRun();
  //  Write file
  rl->WriteHeader("OVERWRITE");
  gener->Write();
  rl->Write();

}


AliGenerator*  CreateGenerator(Tune_t tune, Float_t energy)
{

  

  if (tune == -1) {

    // phojet
    AliGenDPMjet* gener = new AliGenDPMjet(1);
    gener->SetProcess(kDpmMb);
    gener->SetProjectile("P", 1, 1);
    gener->SetTarget("P", 1, 1);

    gener->SetEnergyCMS(energy);
    return gener;
  }

  else if (tune == kHepMC) {

    AliGenReaderHepMC *reader = new AliGenReaderHepMC();
    //reader->SetFileName("crmc_eposlhc_447719733_Pb_Pb_1577.hepmc");    
    //    reader->SetFileName("/Users/mfloris/Work/ALICE/ANALYSIS/current/HMTF/MonteCarlo/DebugHepMC/pythiaPerugia2011.hepmc");// FIXME: make file name settable
    reader->SetFileName(hepMCFile.Data());
    AliGenExtFile *gener = new AliGenExtFile(-1);
    gener->SetReader(reader);
    
    return gener;
    
  }

  else if (tune == kEpos3111) {

    AliGenEposReader *reader = new AliGenEposReader();
    reader->SetFileName(hepMCFile.Data());
    AliGenExtFile *gener = new AliGenExtFile(-1);
    gener->SetReader(reader);

    return gener;    
  }

  if (tune != kPyTuneMonash2013 && tune != kPyTuneAtlasCSC && tune != kPyTuneCDFA && tune != kPyTuneCMS6D6T && tune != kPyTunePerugia0 && tune != kPyTunePerugia2011) {
    
    Printf("Unknown pythia tune, quitting");
    exit(1);

  }
  else {

    AliGenPythia * gener =  new AliGenPythia(1);
    //
    //  set pythia tune

    gener->SetTune(tune);

    //   structure function  
    if(tune == kPyTuneAtlasCSC) {
      std::cout << "Setting structure function" << std::endl;      
      gener->SetStrucFunc(kCTEQ61);
    }
    if(tune == kPyTuneCMS6D6T) {
      std::cout << "Setting structure function" << std::endl;      
      gener->SetStrucFunc(kCTEQ6l);
    }
    if(tune == kPyTunePerugia0) {
      std::cout << "Setting new parton shower" << std::endl;      
      gener->UseNewMultipleInteractionsScenario();
    }
    //   charm, beauty, charm_unforced, beauty_unforced, jpsi, jpsi_chi, mb
    gener->SetProcess(kPyMb);
    //   Centre of mass energy 
    gener->SetEnergyCMS(energy);

    // Set target/projectile // Is this needded?
    gener->SetProjectile("P", 1, 1);
    gener->SetTarget("P", 1, 1);

    return gener;
  }
}
