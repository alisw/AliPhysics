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

AliAnalysisTaskTriggerStudy::AliAnalysisTaskTriggerStudy()
: AliAnalysisTaskSE("TaskTriggerStudy"),
  fESD(0),fHistoList(0),fIsMC(0),fTriggerAnalysis(0)
{
  // constructor

  DefineOutput(1, AliHistoListWrapper::Class());

}
AliAnalysisTaskTriggerStudy::AliAnalysisTaskTriggerStudy(const char * name)
  : AliAnalysisTaskSE(name),
    fESD(0),fHistoList(0),fIsMC(0),fTriggerAnalysis(0)
{
  //
  // Standard constructur which should be used
  //

  DefineOutput(1, AliHistoListWrapper::Class());

}

AliAnalysisTaskTriggerStudy::AliAnalysisTaskTriggerStudy(const AliAnalysisTaskTriggerStudy& obj) : 
  AliAnalysisTaskSE(obj) ,fESD (0), fIsMC(0), fTriggerAnalysis(0)
{
  //copy ctor
  fESD = obj.fESD ;
  fHistoList = obj.fHistoList;
  fTriggerAnalysis = obj.fTriggerAnalysis;
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

  // FIXME: two options here: either we add a loop setting the name of
  // the histos with the trigger class in them (there may be a global
  // SetHistoNamePrefix) or I put a cut to select only collision
  // classes
  
  // get the multiplicity object
  const AliMultiplicity* mult = fESD->GetMultiplicity();
  Int_t ntracklets = mult->GetNumberOfTracklets();

  GetHistoTracklets("all","All events")->Fill(ntracklets);

  // Fast or in the outer layer  
  Int_t nFastOrOnline  = fTriggerAnalysis->SPDFiredChips(fESD, 1, 0, 2); // offline
  Int_t nFastOrOffline = fTriggerAnalysis->SPDFiredChips(fESD, 0, 0, 2); // online

  if(nFastOrOffline != nFastOrOnline) {
    GetHistoTracklets("mismatchingFastOr", "Events where fast or offline differs from fast-or online")->Fill(ntracklets);
  }
  
  Bool_t c0sm1 = nFastOrOffline >= 1;
  Bool_t c0sm2 = nFastOrOffline >= 2;
  Bool_t c0sm3 = nFastOrOffline >= 3;
  Bool_t c0sm4 = nFastOrOffline >= 4;
  Bool_t c0sm5 = nFastOrOffline >= 5;
  
  // V0 triggers
  Bool_t c0v0A       = fTriggerAnalysis->IsOfflineTriggerFired(fESD, AliTriggerAnalysis::kV0A);
  Bool_t c0v0C       = fTriggerAnalysis->IsOfflineTriggerFired(fESD, AliTriggerAnalysis::kV0C);
  Bool_t v0AHW     = (fTriggerAnalysis->V0Trigger(fESD, AliTriggerAnalysis::kASide, kTRUE) == AliTriggerAnalysis::kV0BB);// should replay hw trigger
  Bool_t v0CHW     = (fTriggerAnalysis->V0Trigger(fESD, AliTriggerAnalysis::kCSide, kTRUE) == AliTriggerAnalysis::kV0BB);// should replay hw tr

  if (c0v0A != v0AHW){
    GetHistoTracklets("mismatchingV0A", "Events where V0A offline differs from V0A online")->Fill(ntracklets);
  }

  if (c0v0C != v0CHW){
    GetHistoTracklets("mismatchingV0C", "Events where V0C offline differs from V0C online")->Fill(ntracklets);
  }
  
  // TOF triggers 
  // FIXME: implement real triggers
  Bool_t c0OM2 = kFALSE;
  Bool_t c0OM3 = kFALSE;

  // Some macros for the online triggers
  Bool_t cMBS2A = c0sm2 && c0v0A;
  Bool_t cMBS2C = c0sm2 && c0v0C;
  Bool_t cMBAC  = c0v0A && c0v0C;

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

  // 	// FIXME: add kTransportBit (uncomment below)
  // 	if(!fMCEvent->Stack()->IsPhysicalPrimary(ipart)) continue;

  // 	//check if current particle is a physical primary
  // 	// Bool_t physprim=fMCEvent->IsPhysicalPrimary(label);
  // 	// if (!physprim) continue;
  // 	// if (!track) return kFALSE;
  // 	// Bool_t transported = mcPart->Particle()->TestBit(kTransportBit);
  // 	// if(!transported) return kFALSE;
 
  // 	// Get MC vertex
  // 	//FIXME: which vertex do I take for MC?
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
  // terminate

  TString hname = "h";
  hname+=name;
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


