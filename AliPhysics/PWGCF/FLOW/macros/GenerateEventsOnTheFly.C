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
 * at the end of the macro, an example function is given that can be used 
 * to connect the on the fly events to the AliAnalysistwoParticleResonanceFlowTask
 *
 * Redmer Alexander Bertens (rbertens@cern.ch)
 * Carlos Eugenio Perez Lara
 * Andrea Dubla
 */

class AliFlowOnTheFlyGenerator;
class AliAnalysisTwoParticleResonanceFlowTask;
class AliFlowAnalysisWithQCumulants;
// main function
//
// you can call this macro in different modes. 
// a) default mode will launch the event genetor and write the output to file
// b) called with kFALSE the macro will read events from file and do an example analysis 
// this could be expanded to support more modes (e.g. not writing to file but immediately
// doing the desired v2 analysis)
void GenerateEventsOnTheFly(Bool_t generate = kTRUE) {
  if(!generate) {
      LoadLibraries();
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
  const int events      = 50000;
  // i) write output to file
  Bool_t writeOutput    = kTRUE;
  // j) do qa. qa histograms will be created at the end of generation adn written to a rootfile in your working directory
  Bool_t qa             = kTRUE;
  // k) write analysis to an output file
  TString trunkName             =       "OnTheFlyEvent"; // trunk of output file name
  // l) specify the max number of events per file (new files will be generated automatically)
  const int maxEventsPerFile    =    1000;          // specify the maximum number of events that is stored per file
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
                                                                                        3,              // introduce flow fluctuations
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
    Int_t nFiles(0);                                    // file counter
    // setup analysis methods
    AliFlowAnalysisWithQCumulants* qc = new AliFlowAnalysisWithQCumulants();
    PrepareQCumulants(qc);

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
            delete event;
        }
        nFiles++;
    }
    // prepare output for the flow package analyses
    TString outputFileName = "AnalysisResults.root";  
    TFile *outputFile = new TFile(outputFileName.Data(),"RECREATE");
    const Int_t nMethods(1);
    TString method[] = {"QC"};
    TDirectoryFile *dirFileFinal[nMethods] = {NULL};
    TString fileName[nMethods]; 
    for(Int_t i(0); i < nMethods; i++) {
         fileName[i]+="output";
         fileName[i]+=method[i].Data();
         fileName[i]+="analysis";
         dirFileFinal[i] = new TDirectoryFile(fileName[i].Data(),fileName[i].Data());
    } 
    qc->Finish();qc->WriteHistograms(dirFileFinal[0]);
    outputFile->Close();
    delete outputFile;
}

/* ====================================================
 *
 * the following two functions serve as an exmaple
 * of how to run the phi and kstar analysis on on the 
 * fly events 
 * 
 * ====================================================
 */

//_____________________________________________________________________________
void TwoParticleResonanceFlowOnTheFly() {
    
    // load additional libraries
    gSystem->Load("libGeom");
    gSystem->Load("libVMC");
    gSystem->Load("libXMLIO");
    gSystem->Load("libPhysics");
    gSystem->Load("libANALYSIS");
    gSystem->Load("libANALYSISalice");
    gSystem->Load("libPWGflowBase");
    gSystem->Load("libPWGflowTasks");

    TString dirName = "reconstruction";
    Int_t nFiles(0);
    // define the poi cuts, see AddTwoParticelResonanceFlowTask.C for more info
    const Int_t mb(30); // no of massbands available for analysis
    Float_t minMass(.99), maxMass(1.092);       // upper and lower bound
    AliFlowTrackSimpleCuts* POIfilterQC[mb];    // pointers to poi filters
    Double_t flowBands[2][mb];                  // define the invm regins
    Double_t _inc = (maxMass-minMass)/(float)mb;
    for (Int_t _mb = 0; _mb < mb; _mb++) {      // create the poi filters
       flowBands[0][_mb] = minMass + _mb * _inc;
       flowBands[1][_mb] = minMass + (_mb + 1) * _inc;
       POIfilterQC[_mb] = new AliFlowTrackSimpleCuts(Form("FilterPOIQC_MB%d", _mb));
       POIfilterQC[_mb]->SetEtaMin(-0.8);       // eta range
       POIfilterQC[_mb]->SetEtaMax(+0.8);       
       POIfilterQC[_mb]->SetMassMin(flowBands[0][_mb]); // invm range lower bound
       POIfilterQC[_mb]->SetMassMax(flowBands[1][_mb]); // invm rnage upper bound
    }
    // do the flow analysis
    AliFlowAnalysisWithQCumulants* qc[mb];
    for(int i(0); i < mb; i++) {        // init the q-cumulant tasks, one for each invm bin
        qc[i] = new AliFlowAnalysisWithQCumulants();
        PrepareQCumulants(qc[i]);
        printf("  > init qc task %i < \n", i);
    }
    AliAnalysisTwoParticleResonanceFlowTask* task = new AliAnalysisTwoParticleResonanceFlowTask("onthefly");
    SetUpTask(task, minMass,maxMass);// setup the task
    // open files
    while(kTRUE) {      // infinite loop which we will break when we've looked through all the files
        TFile file(Form("OnTheFlyEvent_%i.root", nFiles));   // open the root file
        if(nFiles==0&&file.IsZombie()) {                     // something went wrong ...
            printf(" > cannot get input files ... exiting < \n");
            exit(0);
        }
        if(file.IsZombie()) break;                  // break the loop on the first empty file
        TIter iter(file.GetListOfKeys());           // get a list of keys
        while(kTRUE) {                              // infinite loop over events ...
            TKey* key = iter();                     // get the next key from the file
            if(!key) break;                         // ... and exit the loop if the key is empty
            AliFlowEventSimple* event = (AliFlowEventSimple*)key->ReadObj();        // cast to a flow event
            printf("   - read event with %i tracks from file %i \n       > task ", event->NumberOfTracks(), nFiles);
            task->DoAnalysisOnTheFly(event);                       // do the on the fly analysis
            AliFlowEventSimple* flowEvent = task->GetFlowEvent();  // retrieve the flow event
            for(int j(0); j < mb; j++) {                           // loop over all invm bands
                flowEvent->TagPOI(POIfilterQC[j]);                 // 'offline' tagging of poi's in certain mass range
                flowEvent->TagSubeventsInEta(-.8, 0, 0, .8);       // setup subevents
                qc[j]->Make(flowEvent);                            // do qc analysis
                printf(" %i", j);
            }
        }
        nFiles++;
    }
    // prepare output for the flow package analyses
    TFile *outputFile = new TFile("AnalysisResults.root","RECREATE");   // common outputfile    
    TDirectoryFile* dirFileFinal[mb];                                   // tdirectory files
    TString fileName[mb];                                               // dir names
    for(Int_t i(0); i < mb; i++) {                                      // loop over all bands ...
         fileName[i]+=Form("QC_minv_%i", i);
         dirFileFinal[i] = new TDirectoryFile(fileName[i].Data(),fileName[i].Data());
         qc[i]->Finish();                               // finalize the method
         qc[i]->WriteHistograms(dirFileFinal[i]);       // and write it to file
    }
    // write the output of the phi reconstruction to the same file
    TDirectoryFile* dir = new TDirectoryFile(dirName.Data(), dirName.Data());
    task->DoAnalysisOnTheFly(dir);
    // end of analysis
    outputFile->Close();
    delete outputFile;
}
//_____________________________________________________________________________
void SetUpTask(AliAnalysisTwoParticleResonanceFlowTask* task, Float_t minMass, Float_t maxMass) 
{
    // some magic which is necessary to 'trick' the analysis task into thinking
    // that we're actually doing an analysis on aod events.
    // some of these configurations are used (e.g. setupSpecies, common constants
    // and the binning in pt)
    // others have no meaning (PID, DCA, RP and POI cuts)
    // note that UserCreateOutputObjects is necessary since it initializes
    // the output of the analysis
    gSystem->Load("libPWGflowTasks");
    Float_t PIDconfig[] = {0,0,0,0,0,0,0.3};                    // not used
    task->SetPIDConfiguration(PIDconfig);                       // not used
    task->SetupSpeciesA(321, 1, 4.93676999e-01, 0.15, 5.);      // pid and charge 
    task->SetupSpeciesB(-321, -1, 4.93676999e-01, 0.15, 5.);    // pid and charge
    task->SetRPCuts(new AliFlowTrackCuts("unused_rp_cuts"));    // not used
    task->SetPOICuts(new AliFlowTrackCuts("unused_poi_cuts"));  // not used
    Float_t POIDCA[] = {0,0,0,0,0};                             // not used
    task->SetPOIDCAXYZ(POIDCA);                                 // not used
    task->SetCommonConstants(30, minMass, maxMass);             // necesssary
    Float_t _pt[] = {0., 0.6, 1.2, 1.8, 2.4, 3., 4., 5., 6., 7.};       // same
    task->SetPtBins(_pt, (Int_t)(sizeof(_pt)/sizeof(_pt[1]))-1);        // same
    task->UserCreateOutputObjects();                            // necessary
}
//_____________________________________________________________________________
void PrepareQCumulants(AliFlowAnalysisWithQCumulants* qc) 
{
    // prepare the cumulants task
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
    qc->Init();
}
//_____________________________________________________________________________
Bool_t LoadLibraries() {
    // load libraries
    gSystem->SetIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TPC -I$ALICE_ROOT/CONTAINERS -I$ALICE_ROOT/STEER -I$ALICE_ROOT/TRD -I$ALICE_PHYSICS/macros -I$ALICE_ROOT/ANALYSIS  -I$ALICE_PHYSICS/OADB -I$ALICE_PHYSICS/PWGHF -I$ALICE_PHYSICS/PWGHF/base -I$ALICE_PHYSICS/PWGHF/vertexingHF -I$ALICE_PHYSICS/PWG/FLOW/Base -I$ALICE_PHYSICS/PWG/FLOW/Tasks -g");
    gSystem->AddIncludePath("-I$ALICE_PHYSICS/include");
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
