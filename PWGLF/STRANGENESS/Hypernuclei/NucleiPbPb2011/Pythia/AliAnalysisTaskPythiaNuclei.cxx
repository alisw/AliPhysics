#include <TChain.h>
#include <TH1I.h>
#include <TH2F.h>
#include <TList.h>

#include "AliAnalysisManager.h"
#include "AliMCEvent.h"
#include "AliVEvent.h"
#include "AliVHeader.h"
#include "AliHeader.h"
#include "AliAODMCHeader.h"
#include "AliAnalysisTaskPythiaNuclei.h"
#include "AliGenEventHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliVParticle.h"
#include "AliStack.h"
#include "AliInputEventHandler.h"
#include "AliMCEvent.h"

ClassImp(AliAnalysisTaskPythiaNuclei)

AliAnalysisTaskPythiaNuclei::AliAnalysisTaskPythiaNuclei()
:fMcEvent(0x0) 
,fMcHandler(0x0)
,fOutputList(0)
,fYPtA(0)
,fYPtM(0)
,fYPA(0)
,fYPM(0) {
  
}

AliAnalysisTaskPythiaNuclei::AliAnalysisTaskPythiaNuclei(const char* name)
:AliAnalysisTaskSE(name)
,fMcEvent(0x0) 
,fMcHandler(0x0)
,fOutputList(0)
,fYPtA(0)
,fYPtM(0)
,fYPA(0)
,fYPM(0) {
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());  
}

AliAnalysisTaskPythiaNuclei::~AliAnalysisTaskPythiaNuclei() {
  if(fOutputList) 
    delete fOutputList;
}

void AliAnalysisTaskPythiaNuclei::UserCreateOutputObjects() {  
  fOutputList = new TList();
  fOutputList->SetOwner(kTRUE);

  fCounter = new TH1I("fCounter",";;Number of Events",1,-0.5,0.5);
  fYPtM = new TH2F("fYPtM",";y;#it{p}_{T} (GeV/#it{c})",4,-1,1,400,0,10);
  fYPM = new TH2F("fYPM",";y;#it{p} (GeV/#it{c})",4,-1,1,400,0,10);
  fYPtA = new TH2F("fYPtA",";y;#it{p}_{T} (GeV/#it{c})",4,-1,1,400,0,10);
  fYPA = new TH2F("fYPA",";y;#it{p} (GeV/#it{c})",4,-1,1,400,0,10);

  fOutputList->Add(fCounter);
  fOutputList->Add(fYPtM);
  fOutputList->Add(fYPM);
  fOutputList->Add(fYPtA);
  fOutputList->Add(fYPA);

  PostData(1,fOutputList);
}

void AliAnalysisTaskPythiaNuclei::Init() {
  fMcHandler = dynamic_cast<AliInputEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()); 
}

void AliAnalysisTaskPythiaNuclei::LocalInit() {
  Init();
}

void AliAnalysisTaskPythiaNuclei::UserExec(Option_t*) {
  Init();

  if (!fMcHandler) {
    PostData(1,fOutputList);
    return;
  }
  
  fMcEvent = fMcHandler->MCEvent(); 
  if (!fMcEvent){
    PostData(1,fOutputList);
    return;
  }

  fCounter->Fill(0);

  for (Int_t iTracks = 0; iTracks < fMCEvent->GetNumberOfTracks(); ++iTracks) {
    AliVParticle* track = fMCEvent->GetTrack(iTracks);
    if (!track) {
      continue;
    }
     
    if (fMCEvent->IsPhysicalPrimary(iTracks) && 
        TMath::Abs(track->Eta()) < 0.8 &&
        TMath::Abs(track->PdgCode()) == 1000010020) {
      if (track->PdgCode() > 0) {
        fYPtM->Fill(track->Y(),track->Pt());
        fYPM->Fill(track->Y(),track->P()); 
      } else {
        fYPtA->Fill(track->Y(),track->Pt());
        fYPA->Fill(track->Y(),track->P()); 
      }
    }
  }

  PostData(1, fOutputList);
}

void AliAnalysisTaskPythiaNuclei::Terminate(Option_t*) {
  AliAnalysisTaskSE::Terminate();
}
