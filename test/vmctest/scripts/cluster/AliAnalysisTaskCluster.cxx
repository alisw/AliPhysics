#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TCanvas.h"
#include "TList.h"
#include "TParticle.h"
#include "TParticlePDG.h"
#include "TProfile.h"
#include "TNtuple.h"
#include "TFile.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliESDEvent.h"
#include "AliLog.h"
#include "AliESDVertex.h"
#include "AliESDInputHandler.h"
#include "AliESDtrackCuts.h"
#include "AliMultiplicity.h"

#include "AliAnalysisTaskCluster.h"
#include "AliExternalTrackParam.h"
#include "AliTrackReference.h"
#include "AliHeader.h"
#include "AliGenEventHeader.h"
#include "AliGenDPMjetEventHeader.h"

// Analysis Task for basic QA on cluster numbers

// Authors: Eva Sicking

ClassImp(AliAnalysisTaskCluster)

  //________________________________________________________________________
  AliAnalysisTaskCluster::AliAnalysisTaskCluster(const char *name) 
    : AliAnalysisTaskSE(name) 
    ,fTrackType(0)
    ,fFieldOn(kTRUE)
    ,fHists(0)
  
    ,fHistRECpt(0)
    ,fITSncl(0)
    ,fTPCncl(0)
    ,fTPCnclSA(0)
    ,fTRDncl(0)

    ,pEtaITSncl(0)
    ,pEtaTPCncl(0)
    ,pEtaTPCnclR(0)
    ,pEtaTRDncl(0)

    ,pPhiITSncl(0)
    ,pPhiTPCncl(0)
    ,pPhiTPCnclR(0)
    ,pPhiTRDncl(0)

    ,pPtITSncl(0)
    ,pPtTPCncl(0)
    ,pPtTPCnclR(0)
    ,pPtTRDncl(0)

    ,fITSlayer(0)
    ,fITSlayerPhi(0)
    ,fITSlayerEta(0)

    ,fCuts(0)

{
 

  DefineOutput(1,  TList::Class()); 

  
  
}


//________________________________________________________________________
void AliAnalysisTaskCluster::UserCreateOutputObjects()
{
  // Create histograms
  // Called once

  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);


  Double_t pt = 20.;

  fHists = new TList();
  fHistRECpt   = new TH1F("fHistRECpt", " p_{T}", 100, 0., pt);

  fITSncl = new TH1F("fITSncl", "fITSncl", 11, -1.5, 9.5);
  fTPCncl = new TH1F("fTPCncl", "fTPCncl", 182, -1.5, 180.5);
  fTPCnclSA = new TH1F("fTPCnclSA", "fTPCnclSA", 182, -1.5, 180.5);
  fTRDncl = new TH1F("fTRDncl", "fTRDncl", 182, -1.5, 180.5);

  //Create Profiles
  pEtaITSncl = new TProfile("pEtaITSncl", "pEtaITSncl", 20, -2, 2);
  pEtaTPCncl = new TProfile("pEtaTPCncl", "pEtaTPCncl", 20, -2, 2);
  pEtaTPCnclR = new TProfile("pEtaTPCnclR", "pEtaTPCnclR", 20, -2, 2, 90, 180);
  pEtaTRDncl = new TProfile("pEtaTRDncl", "pEtaTRDncl", 20, -2, 2);

  pPhiITSncl = new TProfile("pPhiITSncl", "pPhiITSncl", 90, 0, 2*TMath::Pi());
  pPhiTPCncl = new TProfile("pPhiTPCncl", "pPhiTPCncl", 90, 0, 2*TMath::Pi());
  pPhiTPCnclR = new TProfile("pPhiTPCnclR", "pPhiTPCnclR", 90, 0, 2*TMath::Pi(), 90, 180);
  pPhiTRDncl = new TProfile("pPhiTRDncl", "pPhiTRDncl", 90, 0, 2*TMath::Pi());

  pPtITSncl = new TProfile("pPtITSncl", "pPtITSncl", 40, 0, 20);
  pPtTPCncl = new TProfile("pPtTPCncl", "pPtTPCncl", 40, 0, 20);
  pPtTPCnclR = new TProfile("pPtTPCnclR", "pPtTPCnclR", 40, 0, 20, 90, 180);
  pPtTRDncl = new TProfile("pPtTRDncl", "pPtTRDncl", 40, 0, 20);

  fITSlayer = new TH1F("fITSlayer", "fITSlayer", 8, -1.5, 6.5);
  fITSlayerPhi = new TH2F("fITSlayerPhi", "fITSlayerPhi", 8, -1.5, 6.5,90, 0, 2*TMath::Pi());
  fITSlayerEta = new TH2F("fITSlayerEta", "fITSlayerEta", 8, -1.5, 6.5,100, -2, 2);

  fHists->SetOwner();
  fHists->Add(fHistRECpt);
  fHists->Add(fITSncl);
  fHists->Add(fTPCncl);
  fHists->Add(fTPCnclSA);
  fHists->Add(fTRDncl);

  fHists->Add(pEtaITSncl);
  fHists->Add(pEtaTPCncl);
  fHists->Add(pEtaTPCnclR);
  fHists->Add(pEtaTRDncl);

  fHists->Add(pPhiITSncl);
  fHists->Add(pPhiTPCncl);
  fHists->Add(pPhiTPCnclR);
  fHists->Add(pPhiTRDncl);

  fHists->Add(pPtITSncl);
  fHists->Add(pPtTPCncl);
  fHists->Add(pPtTPCnclR);
  fHists->Add(pPtTRDncl);

  fHists->Add(fITSlayer);
  fHists->Add(fITSlayerPhi);
  fHists->Add(fITSlayerEta);

    
  //   for (Int_t i=0; i<fHists->GetEntries(); ++i) {
  //     TH1 *h1 = dynamic_cast<TH1*>(fHists->At(i));
  //     if (h1){
  //       // Printf("%s ",h1->GetName());
  //       h1->Sumw2();
  //     }
  //   }


  TH1::AddDirectory(oldStatus);
}

//__________________________________________________________

void AliAnalysisTaskCluster::UserExec(Option_t *) 
{
  AliVEvent *event = InputEvent();
  if (!event) {
    Printf("ERROR: Could not retrieve event");
    return;
  }


  if(Entry()==0){
    AliESDEvent* esd = dynamic_cast<AliESDEvent*>(event);
    if(esd){
      Printf("We are reading from ESD");
    }
   
  }


  if(fDebug>1)Printf("There are %d tracks in this event", event->GetNumberOfTracks());

  
  const AliVVertex* vertex = event->GetPrimaryVertex();
  Float_t vz = vertex->GetZ();
  if (TMath::Abs(vz) > 10.) return;

  for (Int_t iTrack = 0; iTrack < event->GetNumberOfTracks(); iTrack++) {
   
    AliVParticle *track = event->GetTrack(iTrack);
    AliESDtrack *esdtrack =  dynamic_cast<AliESDtrack*>(track);


    if (!track) {
      Printf("ERROR: Could not receive track %d", iTrack);
      continue;
    }
    
    if (!fCuts->AcceptTrack(esdtrack)) continue;
  

    fHistRECpt->Fill(esdtrack->Pt());
    fITSncl->Fill(esdtrack->GetNcls(0)); //0=ITS
    fTPCncl->Fill(esdtrack->GetNcls(1)); //1=TPC
    fTRDncl->Fill(esdtrack->GetNcls(2)); //2=TRD
    
    const AliExternalTrackParam *tpcP = 0x0;
    tpcP = esdtrack->GetTPCInnerParam();
    if(tpcP){
      fTPCnclSA->Fill(esdtrack->GetNcls(1));
    }

    
    if(esdtrack->GetNcls(0)!=0){
      for(Int_t layer=0;layer<6;layer++){
	if(esdtrack->HasPointOnITSLayer(layer)){
	  fITSlayer->Fill(layer);
	  fITSlayerPhi->Fill(layer, esdtrack->Phi());
	  fITSlayerEta->Fill(layer, esdtrack->Eta());
	}
      }

      pEtaITSncl->Fill(esdtrack->Eta(),esdtrack->GetNcls(0));
      pPhiITSncl->Fill(esdtrack->Phi(),esdtrack->GetNcls(0));
      pPtITSncl->Fill(esdtrack->Pt(),esdtrack->GetNcls(0));
    }
    if(esdtrack->GetNcls(1)!=0){
      pEtaTPCncl->Fill(esdtrack->Eta(),esdtrack->GetNcls(1));
      pPhiTPCncl->Fill(esdtrack->Phi(),esdtrack->GetNcls(1));
      pPtTPCncl->Fill(esdtrack->Pt(),esdtrack->GetNcls(1));

      pEtaTPCnclR->Fill(esdtrack->Eta(),esdtrack->GetNcls(1));
      pPhiTPCnclR->Fill(esdtrack->Phi(),esdtrack->GetNcls(1));
      pPtTPCnclR->Fill(esdtrack->Pt(),esdtrack->GetNcls(1));
    }
    if(esdtrack->GetNcls(2)!=0){
      pPhiTRDncl->Fill(esdtrack->Phi(),esdtrack->GetNcls(2));	
      pEtaTRDncl->Fill(esdtrack->Eta(),esdtrack->GetNcls(2));
      pPtTRDncl->Fill(esdtrack->Pt(),esdtrack->GetNcls(2));
    }


  
  }//first track loop

  
  // Post output data.
  // PostData(1, fHistPt);
  PostData(1, fHists);
}      



//________________________________________________________________________
void AliAnalysisTaskCluster::Terminate(Option_t *) 
{


}  





