#include "TStopwatch.h"
#include "Riostream.h"
#include "TFile.h"


int runFlowOnTheFlyExample(Int_t nEvts=2000, Int_t mult=1000, Float_t v2=0.05, Int_t iseed=7669)
{
  TStopwatch timer;
  timer.Start();
  
  // Load the needed libraries for root (in AliRoot already loaded)
  gSystem->Load("libGeom");
  gSystem->Load("libVMC");
  gSystem->Load("libXMLIO");
  gSystem->Load("libPhysics");
  gSystem->Load("libPWGflowBase");

  fMyTRandom3 = new TRandom3(iseed);   
  gRandom->SetSeed(fMyTRandom3->Integer(65539));

  // Initialize the flow methods for default analysis:
  AliFlowAnalysisWithMCEventPlane *mcep = new AliFlowAnalysisWithMCEventPlane();
  //mcep->SetHarmonic(2); // default is v2
  mcep->Init();

  AliFlowAnalysisWithQCumulants* qc = new AliFlowAnalysisWithQCumulants();
  // qc->SetHarmonic(2); // default is v2
  qc->Init();
  

  // set cuts for the Reference Particles and Particles Of Interest:
  AliFlowTrackSimpleCuts *cutsRP = new AliFlowTrackSimpleCuts();
  //  cutsRP->SetPtMax(ptMaxRP);
  AliFlowTrackSimpleCuts *cutsPOI = new AliFlowTrackSimpleCuts();
  cutsPOI->SetPtMin(0.2);
  cutsPOI->SetPtMax(2.0);

  Printf("starting the main event loop..");
  // create and analyze events 'on the fly':
  for(Int_t i=0; i<nEvts; i++)
    {
      // creating the event with above settings:
      AliFlowEventSimple* event = new AliFlowEventSimple(mult,AliFlowEventSimple::kGenerate);
       event->AddV2(v2);
      //event->TagTracks(cutsRP, cutsPOI);
      event->TagPOI(cutsPOI);
      event->TagRP(cutsRP);
      // event->Print();
      
      // do flow analysis for various methods:
      mcep->Make(event);
      qc->Make(event);
      cout <<"Event: " << i+1 << "\r"; cout.flush();
      delete event;
    } // end of for(Int_t i=0;i<nEvts;i++)
  
 
  // calculate the final results
  mcep->Finish();
  qc->Finish();
 
  // open a new file which will hold the final results of all methods:
  TString outputFileName = "AnalysisResults.root";
  TFile *outputFile = new TFile(outputFileName.Data(),"RECREATE");
  const Int_t nMethods = 2;
  TString method[nMethods] = {"MCEP","QC"};
  TDirectoryFile *dirFileFinal[nMethods] = {NULL};
  TString fileName[nMethods];
  for(Int_t i=0; i<nMethods; i++)
    {
      // form a file name for each method:
      fileName[i]+="output";
      fileName[i]+=method[i].Data();
      fileName[i]+="analysis";
      dirFileFinal[i] = new TDirectoryFile(fileName[i].Data(),fileName[i].Data());
    }
  
  // store the final results
  mcep->WriteHistograms(dirFileFinal[0]);
  qc->WriteHistograms(dirFileFinal[1]);
  
  outputFile->Close();
  
  delete outputFile;
  cout<<endl; cout<<" ---- Fini ---- "<<endl; cout<<endl;
  
  timer.Stop(); timer.Print();
}




