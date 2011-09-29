#include <exception>
#include "TChain.h"
#include "TTree.h"
#include "TVector3.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TH2.h"
#include "TH1.h"
#include "TH3.h"

#include "AliLog.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskPi0v2.h"

#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliEventplane.h"
#include "AliCentrality.h"
#include <iostream>


// Author Daniel Lohner (Daniel.Lohner@cern.ch)

using namespace std;

ClassImp(AliAnalysisTaskPi0v2)


//________________________________________________________________________
    AliAnalysisTaskPi0v2::AliAnalysisTaskPi0v2(const char *name) : AliAnalysisTaskPi0Reconstruction(name),
    fNBinsPhi(1)
{

    // Define input and output slots here
    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());
}

//________________________________________________________________________
AliAnalysisTaskPi0v2::~AliAnalysisTaskPi0v2(){

}

//________________________________________________________________________
void AliAnalysisTaskPi0v2::UserCreateOutputObjects()
{
  // Create histograms

  AliAnalysisTaskPi0Reconstruction::UserCreateOutputObjects();

  TList *fESDList=new TList();
  fESDList->SetName("ESD histograms");
  fESDList->SetOwner(kTRUE);
  fOutputList->Add(fESDList);
  TList *fBGList=new TList();
  fBGList->SetName("Background histograms");
  fBGList->SetOwner(kTRUE);
  fOutputList->Add(fBGList);

  //Adding the histograms to the output container
  Int_t kGCnXBinsSpectra = Int_t((Pi0MassRange[1]-Pi0MassRange[0])*500);  //500 for range 0 - 1
  Double_t kGCfirstXBinSpectra = Pi0MassRange[0];
  Double_t kGClastXBinSpectra = Pi0MassRange[1];

  Int_t kGCnYBinsSpectra = 250;
  Double_t kGCfirstYBinSpectra = 0.;
  Double_t kGClastYBinSpectra = 25.;

  hInvMassPt=new TH2F("ESD_Mother_InvMass_vs_Pt","ESD Invariant Mass vs Pt",kGCnXBinsSpectra, kGCfirstXBinSpectra, kGClastXBinSpectra,kGCnYBinsSpectra, kGCfirstYBinSpectra, kGClastYBinSpectra);
  hInvMassPt->Sumw2();
  fESDList->Add(hInvMassPt);
  hInvMass=new TH1F("ESD_Mother_InvMass","ESD Invariant Mass",kGCnXBinsSpectra, kGCfirstXBinSpectra, kGClastXBinSpectra);
  hInvMass->Sumw2();
  fESDList->Add(hInvMass);
  hBGPt=new TH2F("ESD_Background_InvMass_vs_Pt","ESD Invariant Mass vs Pt",kGCnXBinsSpectra, kGCfirstXBinSpectra, kGClastXBinSpectra,kGCnYBinsSpectra, kGCfirstYBinSpectra, kGClastYBinSpectra);
  hBGPt->Sumw2();
  fBGList->Add(hBGPt);
  hBG=new TH1F("ESD_Background_InvMass","ESD Invariant Mass",kGCnXBinsSpectra, kGCfirstXBinSpectra, kGClastXBinSpectra);
  hBG->Sumw2();
  fBGList->Add(hBG);

  hInvMassmphi=new TH2F**[fNBinsPhi];
  hBGmphi=new TH2F**[fNBinsPhi];
  for(Int_t phi=0;phi<fNBinsPhi;phi++){
      hInvMassmphi[phi]=new TH2F*[fNCentralityBins];
      hBGmphi[phi]=new TH2F*[fNCentralityBins];
      for(Int_t m=0;m<fNCentralityBins;m++){
	  hInvMassmphi[phi][m]=new TH2F(Form("%d%dESD_Mother_InvMass_vs_Pt",phi,m),"ESD Invariant Mass vs Pt",kGCnXBinsSpectra, kGCfirstXBinSpectra, kGClastXBinSpectra,kGCnYBinsSpectra, kGCfirstYBinSpectra, kGClastYBinSpectra);
          hInvMassmphi[phi][m]->Sumw2();
	  fESDList->Add(hInvMassmphi[phi][m]);
	  hBGmphi[phi][m]=new TH2F(Form("%d%dESD_Background_InvMass_vs_Pt",phi,m),"ESD Invariant Mass vs Pt",kGCnXBinsSpectra, kGCfirstXBinSpectra, kGClastXBinSpectra,kGCnYBinsSpectra, kGCfirstYBinSpectra, kGClastYBinSpectra);
          hBGmphi[phi][m]->Sumw2();
	  fBGList->Add(hBGmphi[phi][m]);
      }
  }

  hInclv2PtInvMassCentrality=new TH2F*[fNCentralityBins];
  
  for(Int_t m=0;m<fNCentralityBins;m++){
      hInclv2PtInvMassCentrality[m]=new TH2F(Form("%dInclv2_vs_MInv_vs_Pt",m),"" , kGCnXBinsSpectra, kGCfirstXBinSpectra, kGClastXBinSpectra,kGCnYBinsSpectra, kGCfirstYBinSpectra, kGClastYBinSpectra);
      hInclv2PtInvMassCentrality[m]->Sumw2();
      fESDList->Add(hInclv2PtInvMassCentrality[m]);
  }

  // RP Calculation
  TList *fRPList=new TList();
  fRPList->SetName("Reaction Plane");
  fRPList->SetOwner(kTRUE);
  fOutputList->Add(fRPList);

  hRP=new TH1F("RP" ,"Reaction Plane Angle" , 360, 0, TMath::Pi());
  fRPList->Add(hRP);
  hRPCentrality=new TH2F("RP_Centrality" ,"Reaction Plane Angle" , 20, 0.001,100.001, 180, 0, TMath::Pi());
  fRPList->Add(hRPCentrality);
  hRPSubevents=new TH2F("RP_Subevents" ,"Reaction Plane Angle" , 180, 0, TMath::Pi(), 180, 0, TMath::Pi());
  fRPList->Add(hRPSubevents);
  hRPDeltaRP=new TH2F("DeltaRP_Centrality" ,"Delta Reaction Plane Angle" , 100, -TMath::Pi()/2, TMath::Pi()/2,20,0.001,100.001);
  fRPList->Add(hRPDeltaRP);
  hRPCosDeltaRP=new TH2F("Cos(DeltaRP)" ,"" , 20, 0.001,100.001,110,-1.1,1.1); // range -1.1 - 1.1 else cos(0) will go to the overflow bin!!!!!
  fRPList->Add(hRPCosDeltaRP);

   // Other
  TList *fOtherList=new TList();
  fOtherList->SetName("Charged");
  fOtherList->SetOwner(kTRUE);
  fOutputList->Add(fOtherList);

  hChargedPt=new TH1F*[fNCentralityBins];

  for(Int_t m=0;m<fNCentralityBins;m++){
      hChargedPt[m]=new TH1F(Form("ChargedSpectrum_Pt_%d",m),"Charged Pt",kGCnYBinsSpectra, kGCfirstYBinSpectra, kGClastYBinSpectra);
      hChargedPt[m]->Sumw2();
      fOtherList->Add(hChargedPt[m]);

  }

  // Gamma

  TList *fGammaList=new TList();
  fGammaList->SetName("Gammav2");
  fGammaList->SetOwner(kTRUE);
  fOutputList->Add(fGammaList);

  hGammav2=new TH2F*[fNCentralityBins];

  for(int i=0;i<fNCentralityBins;i++){
      hGammav2[i]=new TH2F(Form("%d_Gamma_v2_vs_Pt",i) ,Form("%d_Gamma_v2_vs_Pt",i),kGCnYBinsSpectra, kGCfirstYBinSpectra, kGClastYBinSpectra,6,0,TMath::Pi()/2);
      hGammav2[i]->Sumw2();
      fGammaList->Add(hGammav2[i]);
  }

  PostData(1, fOutputList);
}

//________________________________________________________________________
void AliAnalysisTaskPi0v2::UserExec(Option_t *) 
{

    AliAnalysisTaskPi0Reconstruction::UserExec("");


    // Event Selection
    if(fEventIsSelected){

	ProcessChargedParticles();

	ProcessGammas();

	ProcessPi0s();

        ProcessEventPlaneResolution();

    }

 
    PostData(1, fOutputList);
}

//________________________________________________________________________
void AliAnalysisTaskPi0v2::ProcessChargedParticles()
{
    AliVParticle *track=0x0;

	for(Int_t ii=0;ii<fInputEvent->GetNumberOfTracks();ii++){
	    track=fInputEvent->GetTrack(ii);
	    hChargedPt[fCentralityBin]->Fill(track->Pt());
	}
}

//________________________________________________________________________
void AliAnalysisTaskPi0v2::ProcessGammas(){
       
    // Process Reconstructed Gammas
   
    for(Int_t firstGammaIndex=0;firstGammaIndex<fConversionGammas->GetEntriesFast();firstGammaIndex++){

	AliAODConversionPhoton *gamma=dynamic_cast<AliAODConversionPhoton*>(fConversionGammas->At(firstGammaIndex));

	// Gamma Phi wrt EP

	Double_t phiwrt=TMath::ACos(TMath::Abs(TMath::Cos(gamma->Phi()-fEPAngle)));
	hGammav2[fCentralityBin]->Fill(gamma->Pt(),phiwrt);
    }

}


//________________________________________________________________________
void AliAnalysisTaskPi0v2::ProcessPi0s(){


    for(Int_t ii=0;ii<fPi0Candidates->GetEntriesFast();ii++){

	AliAODConversionMother *pi0cand=dynamic_cast<AliAODConversionMother*>(fPi0Candidates->At(ii));

	hInvMassPt->Fill(pi0cand->M() ,pi0cand->Pt());
	hInvMass->Fill(pi0cand->M());

	hInclv2PtInvMassCentrality[fCentralityBin]->Fill(pi0cand->M(),pi0cand->Pt(),TMath::Cos(2*(pi0cand->Phi()-fEPAngle)));

	Int_t phibin=GetPhiBin(pi0cand->Phi()-fEPAngle);
	if(!(phibin<0||phibin>=fNBinsPhi)){
	    hInvMassmphi[phibin][fCentralityBin]->Fill(pi0cand->M(),pi0cand->Pt());}

   
    }

    ProcessBGPi0s();
}

//________________________________________________________________________
void AliAnalysisTaskPi0v2::Terminate(Option_t *) 
{
 
    // Draw result to the screen
    fOutputList->Print();
    // Called once at the end of the query
}


//________________________________________________________________________
void AliAnalysisTaskPi0v2::ProcessBGPi0s(){

    for(Int_t ii=0;ii<fBGPi0s->GetEntriesFast();ii++){

	AliAODConversionMother *pi0candidate=dynamic_cast<AliAODConversionMother*>(fBGPi0s->At(ii));

	hBGPt->Fill(pi0candidate->M(),pi0candidate->Pt(),pi0candidate->GetWeight());
	hBG->Fill(pi0candidate->M(),pi0candidate->GetWeight());

	Int_t phibin=GetPhiBin(pi0candidate->Phi()-fEPAngle);
	if(!(phibin<0||phibin>=fNBinsPhi)){
	    hBGmphi[phibin][fCentralityBin]->Fill(pi0candidate->M(),pi0candidate->Pt(),pi0candidate->GetWeight()); }
    }
}


//________________________________________________________________________
void AliAnalysisTaskPi0v2::ProcessEventPlaneResolution()
{
    AliEventplane *fEP=GetEventPlane();

    if(!fEP||!fEP->GetQsub1()||!fEP->GetQsub2()) return;

    // Subevents for Resolution

    Double_t fPsiRP1=fEP->GetQsub1()->Phi()/2;
    Double_t fPsiRP2=fEP->GetQsub2()->Phi()/2;

    // Calculations for Resolution
    Double_t fDeltaPsiRP=fPsiRP1-fPsiRP2;

    // reactionplaneangle + Pi() is the same angle
    if(TMath::Abs(fDeltaPsiRP)>TMath::Pi()/2){
	if(fDeltaPsiRP>0)fDeltaPsiRP-=TMath::Pi();
	else fDeltaPsiRP+=TMath::Pi();
    }

    Double_t cos2deltaRP=TMath::Cos(2*fDeltaPsiRP);

    // FillHistograms
    hRPSubevents->Fill(fPsiRP1,fPsiRP2);
    hRP->Fill(fEPAngle);
    hRPCentrality->Fill(fCentrality,fEPAngle);
    hRPDeltaRP->Fill(fDeltaPsiRP,fCentrality);
    //hRPQxQy->Fill(mQx,mQy);
    hRPCosDeltaRP->Fill(fCentrality,cos2deltaRP);

}


 //________________________________________________________________________
Int_t AliAnalysisTaskPi0v2::GetPhiBin(Double_t phi)
{
    Int_t phibin=-1;


    if(!(TMath::Abs(phi)>=0&&TMath::Abs(phi)<=2*TMath::Pi())){AliError("phi w.r.t. RP out of defined range");return -1;}

    Double_t phiwrtrp=TMath::ACos(TMath::Abs(TMath::Cos(phi)));

    phibin=Int_t(fNBinsPhi*phiwrtrp/(0.5*TMath::Pi()));

    if(phibin==-1){AliError("Phi Bin not defined");}
    return phibin;
}
