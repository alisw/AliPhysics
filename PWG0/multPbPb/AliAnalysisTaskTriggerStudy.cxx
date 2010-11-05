// AliAnalysisTaskTriggerStudy

// Author: Michele Floris, CERN
// TODO:
// - Add chi2/cluster plot for primary, secondaries and fakes


#include "AliAnalysisTaskTriggerStudy.h"
#include "AliESDInputHandler.h"
#include "AliHistoListWrapper.h"
#include "AliAnalysisManager.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "TH1I.h"
#include "TH3D.h"
#include "AliMCParticle.h"
#include "AliGenEventHeader.h"
#include "AliESDCentrality.h"

#include <iostream>
#include "AliTriggerAnalysis.h"
#include "AliMultiplicity.h"
#include "TFile.h"
#include "AliLog.h"

using namespace std;

ClassImp(AliAnalysisTaskTriggerStudy)

const char * AliAnalysisTaskTriggerStudy::kVDNames[] = {"C0MBS2","C0VBA","C0VBC","C0OM2"};       

AliAnalysisTaskTriggerStudy::AliAnalysisTaskTriggerStudy()
: AliAnalysisTaskSE("TaskTriggerStudy"),
  fESD(0),fHistoList(0),fIsMC(0),fTriggerAnalysis(0),fHistoSuffix(""),fNTrackletsCut(1000000)
{
  // constructor

  DefineOutput(1, AliHistoListWrapper::Class());

}
AliAnalysisTaskTriggerStudy::AliAnalysisTaskTriggerStudy(const char * name)
  : AliAnalysisTaskSE(name),
    fESD(0),fHistoList(0),fIsMC(0),fTriggerAnalysis(0),fHistoSuffix(""),fNTrackletsCut(1000000)
{
  //
  // Standard constructur which should be used
  //

  DefineOutput(1, AliHistoListWrapper::Class());

}

AliAnalysisTaskTriggerStudy::AliAnalysisTaskTriggerStudy(const AliAnalysisTaskTriggerStudy& obj) : 
  AliAnalysisTaskSE(obj) ,fESD (0), fIsMC(0), fTriggerAnalysis(0),fHistoSuffix("")
{
  //copy ctor
  fESD = obj.fESD ;
  fHistoList = obj.fHistoList;
  fTriggerAnalysis = obj.fTriggerAnalysis;
  fHistoSuffix = obj.fHistoSuffix;
}

AliAnalysisTaskTriggerStudy::~AliAnalysisTaskTriggerStudy(){
  // destructor

  if(!AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
    if(fHistoList) {
      delete fHistoList;
      fHistoList = 0;
    }
    if(fTriggerAnalysis) {
      delete fTriggerAnalysis;
      fHistoList = 0;
    }
  }
  // Histo list should not be destroyed: fListWrapper is owner!

}
void AliAnalysisTaskTriggerStudy::UserCreateOutputObjects()
{
  // Called once
  fHistoList = new AliHistoListWrapper("histoList","histogram list for trigger studies");
  fTriggerAnalysis = new AliTriggerAnalysis();
}


void AliAnalysisTaskTriggerStudy::UserExec(Option_t *)
{
  // User code

  /* PostData(0) is taken care of by AliAnalysisTaskSE */
  PostData(1,fHistoList);

  fESD = dynamic_cast<AliESDEvent*>(fInputEvent);
  if (strcmp(fESD->ClassName(),"AliESDEvent")) {
    AliFatal("Not processing ESDs");
  }

  
  // get the multiplicity object
  const AliMultiplicity* mult = fESD->GetMultiplicity();
  Int_t ntracklets = mult->GetNumberOfTracklets();

  if(ntracklets > fNTrackletsCut) return;

  // Reset histo suffix and fill reference histograms without any suffix
  fHistoSuffix = "";
  GetHistoTracklets("all","All events")->Fill(ntracklets);

  // Fast or in the outer layer  
  Int_t nFastOrOnline  = fTriggerAnalysis->SPDFiredChips(fESD, 1, 0, 2); // offline
  Int_t nFastOrOffline = fTriggerAnalysis->SPDFiredChips(fESD, 0, 0, 2); // online
  
  Bool_t c0sm1 = nFastOrOffline >= 1;
  Bool_t c0sm2 = nFastOrOffline >= 2;
  Bool_t c0sm3 = nFastOrOffline >= 3;
  Bool_t c0sm4 = nFastOrOffline >= 4;
  Bool_t c0sm5 = nFastOrOffline >= 5;
  
  // V0 triggers
  Bool_t c0v0A       = fTriggerAnalysis->IsOfflineTriggerFired(fESD, AliTriggerAnalysis::kV0A);
  Bool_t c0v0C       = fTriggerAnalysis->IsOfflineTriggerFired(fESD, AliTriggerAnalysis::kV0C);
  Bool_t v0AHW     = (fTriggerAnalysis->V0Trigger(fESD, AliTriggerAnalysis::kASide, kTRUE) == AliTriggerAnalysis::kV0BB);// should replay hw trigger
  Bool_t v0CHW     = (fTriggerAnalysis->V0Trigger(fESD, AliTriggerAnalysis::kCSide, kTRUE) == AliTriggerAnalysis::kV0BB);// should replay hw trigger

  // TOF triggers 
  // FIXME: move to triggeranalysis?
  AliESDHeader*h = fESD->GetHeader(); // taken the header from AliESDEvent 
  Bool_t c0OM2 = h->IsTriggerInputFired("0OM2"); // thr >= 2 (input 19)
  Bool_t c0OM3 = h->IsTriggerInputFired("0OM3"); // thr >= 3 (input 20)

  // Some macros for the online triggers
  Bool_t cMBS2A = c0sm2 && c0v0A;
  Bool_t cMBS2C = c0sm2 && c0v0C;
  Bool_t cMBAC  = c0v0A && c0v0C;
  

  Bool_t vdArray[kNVDEntries];
  vdArray[kVDC0MBS2] = c0sm2;
  vdArray[kVDC0VBA]  = c0v0A;
  vdArray[kVDC0VBC]  = c0v0C;
  vdArray[kVDC0OM2]  = c0OM2;

  FillTriggerOverlaps("All", "All Events",vdArray);
  

  // loop over trigger classes in the event
  TObjArray * tokens = 0;
  if(fIsMC) {
    // in case of montecarlo I override the trigger class
    tokens = new TObjArray;
    tokens->SetOwner();
    //    tokens->Add(new TObjString("CINT1B-ABCE-NOPF-ALL")); 
    tokens->Add(new TObjString("MC")); 
  }
  else {  
    TString trgClasses = fESD->GetFiredTriggerClasses();
    tokens = trgClasses.Tokenize(" ");
  }
  TIter iter(tokens);
    
  while(TObjString * tok = (TObjString*) iter.Next()){
    // clean up trigger name
    TString trg = tok->GetString();
    trg.Strip(TString::kTrailing, ' ');
    trg.Strip(TString::kLeading, ' ');

    fHistoSuffix = "_";
    fHistoSuffix += trg;

    // Fill histograms mismatchs
    // TODO: check mismatch trigger class 
    if(nFastOrOffline != nFastOrOnline) {
      GetHistoTracklets("mismatchingFastOr", "Events where fast or offline differs from fast-or online")->Fill(ntracklets);
    }
    
    if (c0v0A != v0AHW){
      GetHistoTracklets("mismatchingV0A", "Events where V0A offline differs from V0A online")->Fill(ntracklets);
    }
    
    if (c0v0C != v0CHW){
      GetHistoTracklets("mismatchingV0C", "Events where V0C offline differs from V0C online")->Fill(ntracklets);
    }    
    
    // Fill a tracklet histo for each trigger type
    if(c0sm1)  GetHistoTracklets("c0sm1" ,"Events were trigger c0sm1 fired" )->Fill(ntracklets);
    if(c0sm2)  GetHistoTracklets("c0sm2" ,"Events were trigger c0sm2 fired" )->Fill(ntracklets);
    if(c0sm3)  GetHistoTracklets("c0sm3" ,"Events were trigger c0sm3 fired" )->Fill(ntracklets);
    if(c0sm4)  GetHistoTracklets("c0sm4" ,"Events were trigger c0sm4 fired" )->Fill(ntracklets);
    if(c0sm5)  GetHistoTracklets("c0sm5" ,"Events were trigger c0sm5 fired" )->Fill(ntracklets);
    if(c0OM2)  GetHistoTracklets("c0OM2" ,"Events were trigger c0OM2 fired" )->Fill(ntracklets);
    if(c0OM3)  GetHistoTracklets("c0OM3" ,"Events were trigger c0OM3 fired" )->Fill(ntracklets);
    if(c0v0A)  GetHistoTracklets("c0v0A" ,"Events were trigger c0v0A fired" )->Fill(ntracklets);
    if(c0v0C)  GetHistoTracklets("c0v0C" ,"Events were trigger c0v0C fired" )->Fill(ntracklets);
    if(cMBS2A) GetHistoTracklets("cMBS2A","Events were trigger cMBS2A fired")->Fill(ntracklets);
    if(cMBS2C) GetHistoTracklets("cMBS2C","Events were trigger cMBS2C fired")->Fill(ntracklets);
    if(cMBAC ) GetHistoTracklets("cMBAC ","Events were trigger cMBAC  fired")->Fill(ntracklets);
    //  if() GetHistoTracklets("","Events were trigger  fired");
    
    // Fill trigger overlaps
    FillTriggerOverlaps("All", "All Events in trigger class",vdArray);

    delete tokens;
  }
    
    // if (fIsMC) {
    

  //   if (!fMCEvent) {
  //     AliError("No MC info found");
  //   } else {
      
  //     //loop on the MC event
  //     //      Int_t nMCTracks = fMCEvent->GetNumberOfTracks();
  //     Int_t offset    = fMCEvent->GetPrimaryOffset();
  //     Int_t nMCTracks = fMCEvent->GetNumberOfPrimaries()+offset;
  //     for (Int_t ipart=offset; ipart<nMCTracks; ipart++) { 
	
  // 	AliMCParticle *mcPart  = (AliMCParticle*)fMCEvent->GetTrack(ipart);
	
  // 	// We don't care about neutrals and non-physical primaries
  // 	if(mcPart->Charge() == 0) continue;

  // PHYSICAL PRIMARY
  // 	// Get MC vertex
    // 	TArrayF   vertex;
  // 	fMCEvent->GenEventHeader()->PrimaryVertex(vertex);
  // 	Float_t zv = vertex[2];
  // 	//	Float_t zv = vtxESD->GetZ();
  // 	// Fill generated histo
  // 	hTracks[AliAnalysisMultPbTrackHistoManager::kHistoGen]->Fill(mcPart->Pt(),mcPart->Eta(),zv);
	
  //     }
  //   }
  // }
  



}

void   AliAnalysisTaskTriggerStudy::Terminate(Option_t *){
  // terminate
  // Save output in a more friendly format
  fHistoList = dynamic_cast<AliHistoListWrapper*> (GetOutputData(1));
  if (!fHistoList){
    Printf("ERROR: fHistoList not available");
    return;
  }
  TFile * f = new TFile("trigger_study.root", "recreate");
  fHistoList->GetList()->Write();
  f->Close();

}

TH1 *   AliAnalysisTaskTriggerStudy::GetHistoTracklets(const char * name, const char * title){
  // Book histo of events vs ntracklets, if needed

  TString hname = "hTracklets_";
  hname+=name;  
  hname+=fHistoSuffix;
  TH1 * h = (TH1*) fHistoList->GetList()->FindObject(hname.Data());
  
  if(!h) {
    AliInfo(Form("Booking histo %s",hname.Data()));
    Bool_t oldStatus = TH1::AddDirectoryStatus();
    TH1::AddDirectory(kFALSE);
    h = new TH1F (hname.Data(), title, 50, 0.5, 200);
    h->Sumw2();
    h->SetXTitle("ntracklets");
    fHistoList->GetList()->Add(h);
    TH1::AddDirectory(oldStatus);
  }
  return h;
}

void AliAnalysisTaskTriggerStudy::FillTriggerOverlaps (const char * name, const char * title, Bool_t * vdArray){
  //Fills a histo with the different trigger statistics in a venn like diagramm. Books it if needed.

  // Get or book histo
  TString hname = "hTrigStat_";
  hname+=name;  
  hname+=fHistoSuffix;
  TH1 * h = (TH1*) fHistoList->GetList()->FindObject(hname.Data());
  
  if(!h) {
    AliInfo(Form("Booking histo %s",hname.Data()));
    Bool_t oldStatus = TH1::AddDirectoryStatus();
    TH1::AddDirectory(kFALSE);
    Int_t nbins = 0;
    for(Int_t ientry = 0; ientry < kNVDEntries; ientry++){
      nbins = nbins | (1<<ientry);
    }
    
    h = new TH1I (hname, title, nbins, -0.5, nbins-0.5);
    fHistoList->GetList()->Add(h);
    TH1::AddDirectory(oldStatus);
  
    // we look at the combinations of n triggers
    // We set a bit for each trigger to fill the diagram
    // This is much simpler and faster than any recursive function
    h->GetXaxis()->SetBinLabel(1,"NONE"); 
    for(Int_t ibin = 1; ibin < nbins; ibin++){
      TString binname = "";
      Bool_t first = kTRUE;
      for(Int_t ivdentry = 0; ivdentry < kNVDEntries; ivdentry++){
	if (ibin & (1<<ivdentry)) {
	  if(!first) binname += " & ";
	  binname += kVDNames[ivdentry];
	  first=kFALSE;
	}
      }
      h->GetXaxis()->SetBinLabel(ibin+1,binname.Data());
    }
    
  }

  UInt_t mask = 0;
  for(Int_t ivdentry = 0; ivdentry < kNVDEntries; ivdentry++){
    if(vdArray[ivdentry]) {
      mask  = mask | (1<<ivdentry);
      //      cout << " 1 "   ;
    } //else cout << " 0 ";
  }
  //  cout << hex << " = " << mask << endl;
  
  h->Fill(mask);

}
