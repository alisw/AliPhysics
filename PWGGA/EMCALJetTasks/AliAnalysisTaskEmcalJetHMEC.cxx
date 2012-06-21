// To make jet hadron correlations
// Author: Megan Connors

#include "TChain.h"
#include "TTree.h"
#include "TList.h"
#include "TH1F.h"
#include "TH2F.h"
#include "THnSparse.h"
#include "TCanvas.h"
#include <TClonesArray.h>
#include <TParticle.h>
#include "AliVTrack.h"
#include "TParameter.h"

#include "AliAODEvent.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDCaloCluster.h"
#include "AliESDVertex.h"
#include "AliCentrality.h"
#include "AliAODJet.h"
#include "AliEmcalJet.h"
#include "AliESDtrackCuts.h"
#include "AliAnalysisTaskEmcalJetHMEC.h"
#include "TVector3.h"

ClassImp(AliAnalysisTaskEmcalJetHMEC)

//________________________________________________________________________
AliAnalysisTaskEmcalJetHMEC::AliAnalysisTaskEmcalJetHMEC() : 
  AliAnalysisTaskSE(),
  fTracksName("tracks"),
  fJetsName("jets"),
  fPhimin(-10), 
  fPhimax(10),
  fEtamin(-0.9), 
  fEtamax(0.9),
  fAreacut(0.0),
  fESD(0), 
  fOutputList(0), 
  fHistTrackPt(0),
  fHistCentrality(0), 
  fHistJetEtaPhi(0), 
  fHistTrackEtaPhi(0), 
  fHistJetHEtaPhi(0) 
{
  // Default Constructor

  for(Int_t icent = 0; icent<6; ++icent){
    fHistJetPt[icent]=0;
    fHistJetPt[icent]=0;
    for(Int_t iptjet = 0; iptjet<3; ++iptjet){
      for(Int_t ieta = 0; ieta<3; ++ieta){	
	fHistJetH[icent][iptjet][ieta]=0;
	fHistJetHBias[icent][iptjet][ieta]=0;
      }
    }
  }

}
//________________________________________________________________________
AliAnalysisTaskEmcalJetHMEC::AliAnalysisTaskEmcalJetHMEC(const char *name) : 
  AliAnalysisTaskSE(name),
  fTracksName("tracks"),
  fJetsName("jets"),
  fPhimin(-10), 
  fPhimax(10),
  fEtamin(-0.9), 
  fEtamax(0.9),
  fAreacut(0.0),
  fESD(0), 
  fOutputList(0), 
  fHistTrackPt(0),
  fHistCentrality(0), 
  fHistJetEtaPhi(0), 
  fHistTrackEtaPhi(0), 
  fHistJetHEtaPhi(0) 
{
  // Constructor
  for(Int_t icent = 0; icent<6; ++icent){
    fHistJetPt[icent]=0;
    fHistJetPtBias[icent]=0;
    for(Int_t iptjet = 0; iptjet<3; ++iptjet){
      for(Int_t ieta = 0; ieta<3; ++ieta){	
	fHistJetH[icent][iptjet][ieta]=0;
	fHistJetHBias[icent][iptjet][ieta]=0;
      }
    }
  }

  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());

}

//________________________________________________________________________
void AliAnalysisTaskEmcalJetHMEC::UserCreateOutputObjects()
{
  // Called once

 
  AliVEventHandler* handler = AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler();
  if (!handler) {
    AliError("Input handler not available!");
    return;
  }

  OpenFile(1);
  fOutputList = new TList();
  fOutputList->SetOwner();

  // Create histograms

  fHistTrackPt = new TH1F("fHistTrackPt", "P_{T} distribution", 1000, 0.0, 100.0);



  fHistCentrality = new TH1F("fHistCentrality","centrality",100,0,100);

  fHistJetEtaPhi = new TH2F("fHistJetEtaPhi","Jet eta-phi",900,-1.8,1.8,640,-3.2,3.2);
  fHistTrackEtaPhi = new TH2F("fHistTrackEtaPhi","Track eta-phi",900,-1.8,1.8,640,-3.2,3.2);
  fHistJetHEtaPhi = new TH2F("fHistJetHEtaPhi","Jet-Hadron deta-dphi",900,-1.8,1.8,640,-1.6,4.8);

  char name[200];

 
  for(Int_t icent = 0; icent<6; ++icent){
    sprintf(name,"fHistJetPt_%i",icent);   
    fHistJetPt[icent] = new TH1F(name,name,200,0,200);
    fOutputList->Add(fHistJetPt[icent]);

    sprintf(name,"fHistJetPtBias_%i",icent);   
    fHistJetPtBias[icent] = new TH1F(name,name,200,0,200);
    fOutputList->Add(fHistJetPtBias[icent]);

    for(Int_t iptjet = 0; iptjet<3; ++iptjet){
      for(Int_t ieta = 0; ieta<3; ++ieta){	
	sprintf(name,"fHistJetH_%i_%i_%i",icent,iptjet,ieta);   
	fHistJetH[icent][iptjet][ieta]=new TH2F(name,name,64,-0.5*TMath::Pi(),1.5*TMath::Pi(),300,0,30);
	fOutputList->Add(fHistJetH[icent][iptjet][ieta]);

	sprintf(name,"fHistJetHBias_%i_%i_%i",icent,iptjet,ieta);   
	fHistJetHBias[icent][iptjet][ieta]=new TH2F(name,name,64,-0.5*TMath::Pi(),1.5*TMath::Pi(),300,0,30);
	fOutputList->Add(fHistJetHBias[icent][iptjet][ieta]);

      }
    }
  }



  fOutputList->Add(fHistTrackPt);
  fOutputList->Add(fHistCentrality);
  fOutputList->Add(fHistJetEtaPhi);
  fOutputList->Add(fHistTrackEtaPhi);
  fOutputList->Add(fHistJetHEtaPhi);


  PostData(1, fOutputList);

}

//________________________________________________________________________

Double_t AliAnalysisTaskEmcalJetHMEC:: RelativePhi(Double_t mphi,Double_t vphi) {


  if (vphi < -1*TMath::Pi()) vphi += (2*TMath::Pi());
  else if (vphi > TMath::Pi()) vphi -= (2*TMath::Pi());
  if (mphi < -1*TMath::Pi()) mphi += (2*TMath::Pi());
  else if (mphi > TMath::Pi()) mphi -= (2*TMath::Pi());
  Double_t dphi = vphi-mphi;
  if (dphi < -1*TMath::Pi()) dphi += (2*TMath::Pi());
  else if (dphi > TMath::Pi()) dphi -= (2*TMath::Pi());

  return dphi;//dphi in [-Pi, Pi]                                                                                                    
}


//________________________________________________________________________
Int_t AliAnalysisTaskEmcalJetHMEC::GetCentBin(Double_t cent) const 
{
  // Get centrality bin.

  Int_t centbin = -1;
  if (cent>=0 && cent<10)
    centbin = 0;
  else if (cent>=10 && cent<20)
    centbin = 1;
  else if (cent>=20 && cent<30)
    centbin = 2;
  else if (cent>=30 && cent<40)
    centbin = 3;
  else if (cent>=40 && cent<50)
    centbin = 4;
  else if (cent>=50 && cent<90)
    centbin = 5;
  return centbin;
}


//________________________________________________________________________
Int_t AliAnalysisTaskEmcalJetHMEC::GetEtaBin(Double_t eta) const 
{
  // Get eta bin for histos.

  Int_t etabin = -1;
  if (TMath::Abs(eta)<=0.4)
    etabin = 0;
  else if (TMath::Abs(eta)>0.4 && TMath::Abs(eta)<0.8)
    etabin = 1;
  else if (TMath::Abs(eta)>=0.8)
    etabin = 2;
  return etabin;
}
//________________________________________________________________________
Int_t AliAnalysisTaskEmcalJetHMEC::GetpTjetBin(Double_t jetpt) const 
{
  // Get eta bin for histos.

  Int_t jetptbin = -1;
  if (jetpt>20)
    jetptbin = 0;
  else if (jetpt>30)
    jetptbin = 1;
  else if (jetpt>45)
    jetptbin = 2;

  return jetptbin;
}


//________________________________________________________________________
void AliAnalysisTaskEmcalJetHMEC::UserExec(Option_t *) 
{


  // Main loop called for each event
 // esd or aod mode
  Bool_t esdMode = kTRUE;
  if (dynamic_cast<AliAODEvent*>(InputEvent()))
    esdMode = kFALSE;


  if (esdMode) {
    // optimization in case autobranch loading is off
    AliAnalysisManager *am = AliAnalysisManager::GetAnalysisManager();
    if (fTracksName == "Tracks")
      am->LoadBranch("Tracks");
  }
 

  //get centrality
  TList *list = InputEvent()->GetList();
  AliCentrality *centrality = InputEvent()->GetCentrality() ;
  Double_t fcent=-1; 
  if(centrality)
    fcent = centrality->GetCentralityPercentile("V0M");
  else
    fcent=99;//probably pp data

  if (fcent<0) {
    AliError(Form("Centrality negative: %f", fcent));
    return;
  }


  fHistCentrality->Fill(fcent);
  Int_t centbin = GetCentBin(fcent);

  TClonesArray *jets = 0;
  TClonesArray *tracks = 0;

  tracks = dynamic_cast<TClonesArray*>(list->FindObject(fTracksName));
  if (!tracks) {
    AliError(Form("Pointer to tracks %s == 0", fTracksName.Data() ));
    return;
  }
  const Int_t Ntracks=tracks->GetEntries();

  jets= dynamic_cast<TClonesArray*>(list->FindObject(fJetsName));
  const Int_t Njets = jets->GetEntries();
 

  

  for (Int_t ijet = 0; ijet < Njets; ijet++)
    {
      AliEmcalJet *jet = static_cast<AliEmcalJet*>(jets->At(ijet));
      
      if (!jet)
	continue;
      
      //pt,eta,phi,centrality
      float jetphi = jet->Phi();
      if (jetphi>TMath::Pi())
	jetphi = jetphi-2*TMath::Pi();
      
      if ((jet->Phi()<fPhimin)||(jet->Phi()>fPhimax))
	continue;
      if ((jet->Eta()<fEtamin)||(jet->Eta()>fEtamax))
	continue;
      //fHistAreavsRawPt[centbin]->Fill(jet->Pt(),jet->Area());
      if (jet->Area()<fAreacut)
	continue;
      //prevents 0 area jets from sneaking by when area cut == 0
      if (jet->Area()==0)
	continue;


      Double_t jetPt = jet->Pt();


      fHistJetPt[centbin]->Fill(jet->Pt());

      fHistJetEtaPhi->Fill(jet->Eta(),jetphi);

      //fHistDeltaPtvsArea->Fill(jetPt,jet->Area());

      if (jetPt<20) 
	continue;

      for (Int_t iTracks = 0; iTracks < Ntracks; iTracks++) 
	{
	  AliVTrack* track = static_cast<AliVTrack*>(tracks->At(iTracks));
	  if (!track) {
	    printf("ERROR: Could not receive track %d\n", iTracks);
	    continue;
	  }

	  if(TMath::Abs(track->Eta())>0.9) continue;

	  fHistTrackPt->Fill(track->Pt());
	  	  
	  if (track->Pt()<0.5)
	    continue;
	  
	  Double_t trackphi = track->Phi();
	  if (trackphi > TMath::Pi())
	    trackphi = trackphi-2*TMath::Pi();

	  Double_t tracketa=track->Eta();
	  Double_t jeteta=jet->Eta();

	  Double_t deta=tracketa-jeteta;
	  Int_t ieta=GetEtaBin(deta);

	  //Jet pt, track pt, dPhi,deta,fcent
	  Double_t dphijh = RelativePhi(jetphi,trackphi);
	  if (dphijh < -0.5*TMath::Pi())
	    dphijh= dphijh+ 2*TMath::Pi();

	  Int_t iptjet=-1;
	  iptjet=GetpTjetBin(jetPt);

	  fHistJetH[centbin][iptjet][ieta]->Fill(dphijh,track->Pt());
	  fHistJetHEtaPhi->Fill(deta,dphijh);
	  fHistTrackEtaPhi->Fill(tracketa,trackphi);
	  if ((jet->MaxTrackPt()>6) || (jet->MaxClusterPt()>6)){
	    fHistJetHBias[centbin][iptjet][ieta]->Fill(dphijh,track->Pt());
	    fHistJetPtBias[centbin]->Fill(jet->Pt());

	  }

	} //track loop
      
    }//jet loop
	  
  
  PostData(1, fOutputList);
}      

//________________________________________________________________________
void AliAnalysisTaskEmcalJetHMEC::Terminate(Option_t *) 
{
  //just terminate

}


