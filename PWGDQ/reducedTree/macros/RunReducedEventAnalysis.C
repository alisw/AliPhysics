const Char_t* gkFriendFilename="";

#include "AliReducedAnalysisTaskSE.h"
#include "AliHistogramManager.h"
#include "AliReducedVarManager.h"
#include "AliReducedEventInfo.h"
#include "AliReducedEventPlaneInfo.h"
#include "TSystem.h"
#include "TChain.h"
#include "TFile.h"
#include "TTimeStamp.h"
#include <iostream>

//_________________________________________________________________
void RunReducedEventAnalysis(AliReducedAnalysisTaskSE* analysis, const Char_t* outputdir, const Char_t* inputfilename, Int_t run=-1, Int_t howMany=10000000, Int_t offset=0) {
  //
  // Configure and run an analysis locally on local trees of AliReducedEvent's
  //
  TTimeStamp start;

  if(run==-1) return;

  TFile* saveFile = new TFile("histograms.root", "RECREATE");
  
  //cout << "Histograms defined:" << endl;
  //analysis.GetHistogramManager()->Print();
  
  // create the tree chain(s)
  Long64_t entries=0;
  TChain* friendChain=0x0;
  //if(gkFriendFilename[0]!='\0') friendChain = new TChain("DstFriendTree");
  TChain* chain = AliReducedVarManager::GetChain(inputfilename, howMany, offset, entries, friendChain, gkFriendFilename);
  if(!chain) return;
  //entries=3;
  //if(kTRUE){  // For merging dstTrees
  //  saveFile->cd();
  //  chain->CloneTree(-1,"fast");
  //  saveFile->Write();
  //  saveFile->Close();
  //  return;
  //}
  AliReducedEventInfo* event = new AliReducedEventInfo();
  AliReducedEventPlaneInfo* eventPlane = 0x0;
  if(friendChain) eventPlane = new AliReducedEventPlaneInfo();
  chain->SetBranchAddress("Event",&event);

  //if(friendChain)
    //friendChain->SetBranchAddress("Event",&eventF);
  
  TTimeStamp startEventLoop;
  
  cout << "Looping over " << entries << " events" << endl;
  
  for(Int_t ie=0; ie<entries; ++ie) {
    chain->GetEntry(ie); 
    if(friendChain) friendChain->GetEntry(ie);
    if(ie%1000==0) 
      cout << "event " << ie << endl;
    
      analysis->SetEvent(event);
      analysis->Process();
  }
  
  TTimeStamp stopEventLoop;
  
  analysis->Finish();
  analysis->GetHistogramManager()->WriteOutput(saveFile);
  
  cout << "Initialization time: " << startEventLoop.GetSec() - start.GetSec() << " seconds" << endl;
  cout << "Looping time: " << stopEventLoop.GetSec() - startEventLoop.GetSec() << " seconds" << endl;
  cout << "Speed       : " << Double_t(stopEventLoop.GetSec() - startEventLoop.GetSec())/Double_t(entries+1.0e-5) << " sec./event" << endl;
}

