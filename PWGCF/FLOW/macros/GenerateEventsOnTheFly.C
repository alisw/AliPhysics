#include "TStopwatch.h"
#include "Riostream.h"
#include "TFile.h"
#include "TTimer.h" 
#include "TPythia6.h"
#include "TMath.h"

/*
 * macro to generate flow events on the fly 
 * if the macro is called with kFALSE as argument,
 * an example function of how to do analysis on the output is called
 *
 * Redmer Alexander Bertens (rbertens@cern.ch)
 * Carlos Eugenio Perez Lara
 * Andrea Dubla
 */

class AliFlowOnTheFlyGenerator;

// main function
//
// you can call this macro in different modes. 
// a) default mode will launch the event genetor and write the output to file
// b) called with kFALSE the macro will read events from file and do an example analysis 
// this could be expanded to support more modes (e.g. not writing to file but immediately
// doing the desired v2 analysis)
void GenerateEventsOnTheFly(Bool_t generate = kTRUE) {
  if(!generate) {
      AnalyseEventsFromFile();
      return;
  }
  TStopwatch aa;
  aa.Start();
  // load libraries
  if(!LoadLibraries()) {
      printf(" > problem loading libraries, maybe you need to fix your includes or pythia path <");
      exit(0);
  }
   ///========= GENERATOR INTERFACE ======///
  // this is the generator interface. here you can specify what your event should look like in terms of 
  // mother particles (define the available species and their multiplicities and define a background)
  // see the decayer interface for decayer specifics
  // in this example macro, pythia is used as a decayer. 
  // however, any decayer that inherits from TVirualMCDecayer can be used
  // 
  // by default, all tracks are made rp's, so 'offline' tagging is necessary. 
  // decayed particles are not saved in the output events. 
  // output is built up from files containing AliFlowEventSimple objects, which in turn consist of 
  // AliFlowSimpleTracks, in which a few things should be noted:
  // 1) Charge() tells whether a particle is a primary (-1) or a secondary (1)
  // 2) Weight() stores the longitudinal momentum (Pz)
  // 3) GetID() gives the pdg value of the genreated track
  // ===================================///
  //
  // a) select which species you want to have in your initial event (pdg code)
  Int_t nSpecies        = 7; // number of poi species
  Int_t species[]       = {     22, 333, 221, 311, 3122, 3334, 313}; // gamma, phi eta, k0, lambda, omega, kstar
  // b) specify their multiplicities per event
  Int_t multiplicities[] = {    10, 10, 10, 10, 10, 10, 10}; // so we get 10 particles of each species here
  // c) set the multiplicity for charged background (mixture of pi, k, pi)
  Int_t setChargedBackground = 300;
  // d) poissonian fluctuations for above multiplicities
  Bool_t fluctuate      = kFALSE;
  // e) which tracks should get a flow boost ?
  Bool_t addMotherV2    = kTRUE;
  Bool_t addDaughterV2  = kFALSE;
  // f) for each species (poi or rp, or species that will be created by the decayer) you can specify 'custom' 
  //    functions for pt, v2 or v3. this will override defaults
  //    exapmle: modify the v2 for species 310
  TString       myFunction   = "0.05*TMath::Sin(x)+.07";
  Int_t         modifyV2    = 310;
  // and call - once an instance of the class is created (see below)
  // AliFlowOnTheFlyEventGenerator::SetPtDependentV2(myFunction.Data(), modifyV2)

  // == done with the flow stuff, proceed to the decayer, qa and writing to file == //

  // g) decayer modes. select which decays needs to be forced (in the example by pythia)
  // note that by default gamma's are 'decayed' to e+ e- pairs
  Int_t decayModes[]    = { TPythia6Decayer::kHadronicD };     // some decay modes
  // h) number of events you want to generate.
  const int events      = 1000;
  // i) write output to file
  Bool_t writeOutput    = kTRUE;
  // j) do qa. qa histograms will be created at the end of generation adn written to a rootfile in your working directory
  Bool_t qa             = kTRUE;
  // k) write analysis to an output file
  TString trunkName             =       "OnTheFlyEvent"; // trunk of output file name
  // l) specify the max number of events per file (new files will be generated automatically)
  const int maxEventsPerFile    =    500;          // specify the maximum number of events that is stored per file
                                                     // note that events are stored temporarily in RAM memory before writing to file
                                                     // so setting this number to a very large value could lead to an unresponsive system
  ///====== END OF GENERATOR INTERFACE===///
  // from here there's nothing that needs to be changed for the event generation 
  
  int nFileCounter(0);
  TClonesArray *event = new TClonesArray("TParticle",5000);
  TObjArray *dataContainer = new TObjArray(maxEventsPerFile);
  dataContainer->SetOwner(kTRUE);
  TFile* _tmpfile = 0x0;
  // setup the decayer
  TPythia6Decayer* decayer = new TPythia6Decayer();
  Int_t nDecayModes(sizeof(decayModes)/sizeof(Int_t));
  // get an estimate of no of tracks per event
  Int_t nMult(setChargedBackground);
  for(Int_t i(0); i < nSpecies; i++) nMult += species[i];
  for(Int_t i(0); i < nDecayModes; i++) decayer->SetForceDecay(decayModes[i]);
  // get an instance of the class that we'll use to generate events
  AliFlowOnTheFlyEventGenerator* eventGenerator = new AliFlowOnTheFlyEventGenerator(    qa,             // make some QA plots
                                                                                        nMult,          // set initial value for evnet multiplcity
                                                                                        decayer,        // specify the decayer (NULL is no decay)
                                                                                        addMotherV2,    // add v2 for mothers
                                                                                        NULL,           // add v3 for mothers  (not implemented yet) 
                                                                                        addDaughterV2,  // add v2 for daughters
                                                                                        NULL);          // add v3 for daughters (not implemented yet)
  // example of how to pass a custom function for v2 for a certain species to the generator
  //  eventGenerator->SetPtDependentV2(myFunction.Data(), modifyV2);
  for(int i(0); i < events; i++) {      // event generator loop
     // prepare I/O necessities 
     if(!_tmpfile) _tmpfile = new TFile(Form("%s_%i.root", trunkName.Data(), nFileCounter), "RECREATE");       // new file, overwrite if already exists
     // for each event, one can embed a TClonesArray of type TParticle*, e.g. do
     //       eventGenerator->EmbedEvent(myEmbeddedEvent);
     // now that the generator is prepared, generate the actual event
     AliFlowEventSimple* _tmpevent = eventGenerator->GenerateOnTheFlyEvent(event, nSpecies, species, multiplicities, setChargedBackground, fluctuate);
     if(i==0) printf(" \n\n > Generating events on the fly, please be patient ... < \n");
     _tmpevent->Write();
     // check if we need to open a new file, or if we've reached the end of the generator cycle
     if((i-nFileCounter*maxEventsPerFile)==(maxEventsPerFile-1) || i == (events-1)) {
         printf("   - writing memory buffer %i to output file ... \n", nFileCounter);
         dataContainer->Write();
         dataContainer->Clear();
         _tmpfile->Close();
	 _tmpfile->Delete();
         _tmpfile = 0x0;
         nFileCounter++;
     }
  }
  if(qa) eventGenerator->DoGeneratorQA(kTRUE, kFALSE);
  delete eventGenerator;
  delete decayer;
  delete event;
  aa.Stop();
  printf("\n > events have been generated and written to file ! < \n");
  aa.Print();
}
//_____________________________________________________________________________
void AnalyseEventsFromFile() {
    // example function: read events from file and do some analysis.
    gSystem->Load("libPWGflowBase");
    Int_t nFiles(0);                                    // file counter
    // setup analysis methods
    AliFlowAnalysisWithQCumulants *qc(0x0);             // for cumulants
    AliFlowAnalysisWithScalarProduct *sp(0x0);          // for scalar product
    sp = new AliFlowAnalysisWithScalarProduct();
    sp->SetHarmonic(2);
    sp->SetApplyCorrectionForNUA(kFALSE);
    sp->Init();
    qc = new AliFlowAnalysisWithQCumulants();
    qc->SetHarmonic(2);
    qc->SetCalculateDiffFlow(kTRUE);
    qc->SetCalculate2DDiffFlow(kFALSE); // vs (pt,eta)
    qc->SetApplyCorrectionForNUA(kFALSE);
    qc->SetFillMultipleControlHistograms(kFALSE);     
    qc->SetMultiplicityWeight("combinations"); // default (other supported options are "unit" and "multiplicity")
    qc->SetCalculateCumulantsVsM(kFALSE);
    qc->SetCalculateAllCorrelationsVsM(kFALSE); // calculate all correlations in mixed harmonics "vs M"
    qc->SetnBinsMult(10000);
    qc->SetMinMult(0);
    qc->SetMaxMult(10000);      
    qc->SetBookOnlyBasicCCH(kFALSE); // book only basic common control histograms
    qc->SetCalculateDiffFlowVsEta(kTRUE); // if you set kFALSE only differential flow vs pt is calculated
    qc->SetCalculateMixedHarmonics(kFALSE); // calculate all multi-partice mixed-harmonics correlators
    qc->Init();  

    while(kTRUE) {
        TFile file(Form("OnTheFlyEvent_%i.root",nFiles));         // open the root file
        if(file.IsZombie()) break;
        TIter iter(file.GetListOfKeys());           // get a list of keys

        // loop over all the events in the file
        while(kTRUE) {                              // infinite loop ...
            TKey* key = iter();                     // get the next key from the file
            if(!key) break;                         // ... and exit the loop if the key is empty
            AliFlowEventSimple* event = (AliFlowEventSimple*)key->ReadObj();        // cast to a flow event
            printf("   - read event with %i tracks\n", event->NumberOfTracks());
            // track loop
            for(int i(0); i < event->NumberOfTracks(); i++) {
                // do some offline fun stuff with the events
                AliFlowTrackSimple* track = event->GetTrack(i);
                if(!track) continue;                // break the loop when we don't find a track
                                                    // do some next stuff with the tracks
                // example: retag the primary kaons as poi's and the rest as reference
                if(track->Charge()!=-1) continue;           // reject secondary tracks UGLY !!!
                if(TMath::Abs(track->GetID()==321)) {       // find the kaons
                    if(track->InRPSelection()) {            // if the kaon is an rp (it should be by default)
                        track->SetForPOISelection(kTRUE);   // set it as poi ...
                        track->SetForRPSelection(kFALSE);   // ... and 'detag' it as rp
                        // make sure to correct the reference multiplicity
                        event->SetReferenceMultiplicity(-1+event->GetReferenceMultiplicity());
                    }
                    else track->SetForPOISelection(kTRUE);  // tag as poi
                }
                else {                                      // if track is not a kaon
                    if(!track->InRPSelection()) {           // check if it's an rp
                        track->SetForRPSelection(kTRUE);    // if not, tag as rp (it should be there, but just to be safe ... )
                        // and again correct the reference multiplicity
                        event->SetReferenceMultiplicity(1+event->GetReferenceMultiplicity());
                    }
                }
            }
            // insert the poor guy into the flow anaysis method
            qc->Make(event);
            sp->Make(event);
            delete event;
        }
        nFiles++;
    }
    // prepare output for the flow package analyses
    TString outputFileName = "AnalysisResults.root";  
    TFile *outputFile = new TFile(outputFileName.Data(),"RECREATE");
    const Int_t nMethods(2);
    TString method[] = {"SP", "QC"};
    TDirectoryFile *dirFileFinal[nMethods] = {NULL};
    TString fileName[nMethods]; 
    for(Int_t i(0); i < nMethods; i++) {
         fileName[i]+="output";
         fileName[i]+=method[i].Data();
         fileName[i]+="analysis";
         dirFileFinal[i] = new TDirectoryFile(fileName[i].Data(),fileName[i].Data());
    } 
    sp->Finish();sp->WriteHistograms(dirFileFinal[0]);
    qc->Finish();qc->WriteHistograms(dirFileFinal[1]);
    outputFile->Close();
    delete outputFile;
}
//_____________________________________________________________________________
Bool_t LoadLibraries() {
    // load libraries
    gSystem->SetIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TPC -I$ALICE_ROOT/CONTAINERS -I$ALICE_ROOT/STEER -I$ALICE_ROOT/TRD -I$ALICE_ROOT/macros -I$ALICE_ROOT/ANALYSIS  -I$ALICE_ROOT/OADB -I$ALICE_ROOT/PWGHF -I$ALICE_ROOT/PWGHF/base -I$ALICE_ROOT/PWGHF/vertexingHF -I$ALICE_ROOT/PWG/FLOW/Base -I$ALICE_ROOT/PWG/FLOW/Tasks -g");
    gSystem->Load("libGeom");
    gSystem->Load("libVMC");
    gSystem->Load("libXMLIO");
    gSystem->Load("libPhysics");
    gSystem->Load("libEG");
    gSystem->Load("libPWGflowBase");
    if(gSystem->Load("$ALICE_ROOT/lib/tgt_linuxx8664gcc/libpythia6")!=0) {
        printf(" \n\n\n *** fatal error, couldn't load pythia, exiting ! *** \n\n\n");
        return kFALSE;
    }
    return kTRUE;
} 
//_____________________________________________________________________________
