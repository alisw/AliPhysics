#include "TH1D.h"
#include "AliAODEvent.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
#include "AliAnalysisTaskMFTExample.h"
#include "TDatabasePDG.h"
#include "AliMUONTrackParam.h"
#include "AliMUONTrackExtrap.h"
#include "AliAODDimuon.h"
#include "AliAODTrack.h"
#include "AliAODHeader.h"
#include "TClonesArray.h"
#include "AliMFTConstants.h"
#include "AliMFTAnalysisTools.h"
#include "TRandom.h"
#include "TList.h"

ClassImp(AliAnalysisTaskMFTExample)

//====================================================================================================================================================

AliAnalysisTaskMFTExample::AliAnalysisTaskMFTExample() : 
  AliAnalysisTaskSE(),
  fVertexMode(0),
  fHistogramList(0),
  fHistPtSingleMuons(0),
  fHistPtSingleMuonsFromJpsi(0),
  fHistPtDimuonsOS(0),
  fHistMassDimuonsOS(0),
  fHistPtDimuonsJpsi(0),
  fHistMassDimuonsJpsi(0),
  fHistResidualXVtxJpsi(0),
  fHistResidualYVtxJpsi(0),
  fHistResidualZVtxJpsi(0) {
  
  // Default constructor

  for (Int_t i=0; i<3; i++) fPrimaryVertex[i] = 0;
  fVtxResolutionITS[0] = 5.e-4;
  fVtxResolutionITS[1] = 5.e-4;
  fVtxResolutionITS[2] = 4.e-4;

}

//====================================================================================================================================================

AliAnalysisTaskMFTExample::AliAnalysisTaskMFTExample(const Char_t *name) : 
  AliAnalysisTaskSE(name),
  fVertexMode(0),
  fHistogramList(0),
  fHistPtSingleMuons(0),
  fHistPtSingleMuonsFromJpsi(0),
  fHistPtDimuonsOS(0),
  fHistMassDimuonsOS(0),
  fHistPtDimuonsJpsi(0),
  fHistMassDimuonsJpsi(0),
  fHistResidualXVtxJpsi(0),
  fHistResidualYVtxJpsi(0),
  fHistResidualZVtxJpsi(0) {

  // Constructor

  for (Int_t i=0; i<3; i++) fPrimaryVertex[i] = 0;
  fVtxResolutionITS[0] = 5.e-4;
  fVtxResolutionITS[1] = 5.e-4;
  fVtxResolutionITS[2] = 4.e-4;

  // Define input and output slots here
  DefineOutput(1, TList::Class());

}

//====================================================================================================================================================

void AliAnalysisTaskMFTExample::UserCreateOutputObjects() {

  // Called once

  fHistogramList = new TList();
  fHistogramList->SetOwner(kTRUE);
  
  fHistPtSingleMuons = new TH1D("fHistPtSingleMuons","p_{T} of single muons (All)", 100, 0, 10);
  fHistPtSingleMuons -> SetXTitle("p_{T}  [GeV/c]");
  fHistPtSingleMuons -> Sumw2();

  fHistPtSingleMuonsFromJpsi = new TH1D("fHistPtSingleMuonsFromJpsi", "p_{T} of single muons (from J/#psi)", 100, 0, 10);
  fHistPtSingleMuonsFromJpsi -> SetXTitle("p_{T}  [GeV/c]");
  fHistPtSingleMuonsFromJpsi -> Sumw2();

  fHistMassDimuonsOS = new TH1D("fHistMassDimuonsOS", "Mass of OS dimuons", 500, 0, 10);
  fHistMassDimuonsOS -> SetXTitle("Mass  [GeV/c^{2}]");
  fHistMassDimuonsOS -> Sumw2();

  fHistMassDimuonsJpsi = new TH1D("fHistMassDimuonsJpsi", "Mass of J/#psi dimuons", 500, 0, 10);
  fHistMassDimuonsJpsi -> SetXTitle("Mass  [GeV/c^{2}]");
  fHistMassDimuonsJpsi -> Sumw2();

  fHistPtDimuonsOS = new TH1D("fHistPtDimuonsOS", "p_{T} of OS dimuons", 100, 0, 10);
  fHistPtDimuonsOS -> SetXTitle("p_{T}  [GeV/c]");
  fHistPtDimuonsOS -> Sumw2();

  fHistPtDimuonsJpsi = new TH1D("fHistPtDimuonsJpsi", "p_{T} of J/#psi dimuons", 100, 0, 10);
  fHistPtDimuonsJpsi -> SetXTitle("p_{T}  [GeV/c]");
  fHistPtDimuonsJpsi -> Sumw2();

  fHistResidualXVtxJpsi = new TH1D("fHistResidualXVtxJpsi", "J/#psi vertex residuals along x", 100, -100, 100);
  fHistResidualXVtxJpsi -> SetXTitle("Residual  [#mum]");
  fHistResidualXVtxJpsi -> Sumw2();

  fHistResidualYVtxJpsi = new TH1D("fHistResidualYVtxJpsi", "J/#psi vertex residuals along y", 100, -100, 100);
  fHistResidualYVtxJpsi -> SetXTitle("Residual  [#mum]");
  fHistResidualYVtxJpsi -> Sumw2();

  fHistResidualZVtxJpsi = new TH1D("fHistResidualZVtxJpsi", "J/#psi vertex residuals along z", 100, -1000, 1000);
  fHistResidualZVtxJpsi -> SetXTitle("Residual  [#mum]");
  fHistResidualZVtxJpsi -> Sumw2();

  fHistogramList->Add(fHistPtSingleMuons);
  fHistogramList->Add(fHistPtSingleMuonsFromJpsi);
  fHistogramList->Add(fHistMassDimuonsOS);
  fHistogramList->Add(fHistMassDimuonsJpsi);
  fHistogramList->Add(fHistPtDimuonsOS);
  fHistogramList->Add(fHistPtDimuonsJpsi);
  fHistogramList->Add(fHistResidualXVtxJpsi);
  fHistogramList->Add(fHistResidualYVtxJpsi);
  fHistogramList->Add(fHistResidualZVtxJpsi);

  PostData(1, fHistogramList);

}

//====================================================================================================================================================

void AliAnalysisTaskMFTExample::UserExec(Option_t *) {

  // Main loop
  // Called for each event	

  AliAODEvent *aodEv = dynamic_cast<AliAODEvent*>(InputEvent());  
  if (!aodEv) return;

  aodEv->GetHeader()->InitMagneticField();
  AliMUONTrackExtrap::SetField();

  TClonesArray *stackMC = (TClonesArray*) (aodEv->GetList()->FindObject(AliAODMCParticle::StdBranchName()));

  // Getting primary vertex, either from the generation or from the reconstruction -------------------

  if (fVertexMode == kGenerated) {
    AliAODMCHeader *mcHeader = (AliAODMCHeader*) (aodEv->GetList()->FindObject(AliAODMCHeader::StdBranchName()));
    mcHeader->GetVertex(fPrimaryVertex);
    for (Int_t i=0; i<3; i++) fPrimaryVertex[i] = gRandom->Gaus(fPrimaryVertex[i], fVtxResolutionITS[i]);
  }
  else if (fVertexMode == kReconstructed) {
    aodEv->GetPrimaryVertex()->GetXYZ(fPrimaryVertex);
  }

  // -------------------------------------------------------------------------------------------------

  AliAODTrack *recMuon1=0, *recMuon2=0;
  AliAODMCParticle *mcMuon1=0, *mcMuon2=0;
  AliAODMCParticle *mcMother1=0;
  Bool_t isMuon1FromJpsi=0;

  // Loop over MUON+MFT muons

  //--------------------------------------------------------------------------------------------------

  for (Int_t iTrack=0; iTrack<aodEv->GetNumberOfTracks(); iTrack++) { 

    if (!(aodEv->GetTrack(iTrack)->IsMuonGlobalTrack())) continue;

    recMuon1 = aodEv->GetTrack(iTrack);
    isMuon1FromJpsi = kFALSE;
    
    // pt all muons
    fHistPtSingleMuons -> Fill(recMuon1->Pt());

    // pt muons from J/psi
    if (recMuon1->GetLabel() >= 0) {
      mcMuon1 = (AliAODMCParticle*) stackMC->At(recMuon1->GetLabel());
      if (mcMuon1) {
	if (mcMuon1->GetMother()>=0) {
	  mcMother1 = (AliAODMCParticle*) stackMC->At(mcMuon1->GetMother());
	  if (mcMother1->PdgCode()==443) {
	    isMuon1FromJpsi = kTRUE;
	    fHistPtSingleMuonsFromJpsi -> Fill(recMuon1->Pt());
	  }
	}
      }
    }

    for (Int_t jTrack=0; jTrack<iTrack; jTrack++) { 

      if (!(aodEv->GetTrack(jTrack)->IsMuonGlobalTrack())) continue;
      
      recMuon2 = aodEv->GetTrack(jTrack);
      
      AliAODDimuon *dimuon = new AliAODDimuon(recMuon1, recMuon2);
      if (dimuon->Charge()) {
	delete dimuon;
	continue;
      }

      // pt and mass all OS dimuons
      fHistPtDimuonsOS   -> Fill(dimuon->Pt());
      fHistMassDimuonsOS -> Fill(dimuon->Mass());

      delete dimuon;

      // pt and mass J/psi dimuons
      if (!isMuon1FromJpsi) continue;
      if (recMuon2->GetLabel() >= 0) {
	mcMuon2 = (AliAODMCParticle*) stackMC->At(recMuon2->GetLabel());
	if (mcMuon2) {
	  if (mcMuon2->GetMother() == mcMuon1->GetMother()) {
	    AliAODDimuon *dimuonJpsi = new AliAODDimuon;
	    dimuonJpsi->SetMuons(recMuon1,recMuon2);
	    Double_t pca[3]={0};
	    Double_t pcaQuality=0;
	    TLorentzVector kinem(0,0,0,0);
	    if (!AliMFTAnalysisTools::CalculatePCA(dimuonJpsi, pca, pcaQuality, kinem)) {
	      delete dimuonJpsi;
	      continue;
	    }
	    fHistPtDimuonsJpsi    -> Fill(kinem.Pt());
	    fHistMassDimuonsJpsi  -> Fill(kinem.M());
	    fHistResidualXVtxJpsi -> Fill(1.e4*(pca[0] - fPrimaryVertex[0]));
	    fHistResidualYVtxJpsi -> Fill(1.e4*(pca[1] - fPrimaryVertex[1]));
	    fHistResidualZVtxJpsi -> Fill(1.e4*(pca[2] - fPrimaryVertex[2]));
	    delete dimuonJpsi;
	  }
	}
      }
   
    }   // end of loop on 2nd muon

  }   // end of loop on 1st muon
 
  //--------------------------------------------------------------------------------------------------
   
  PostData(1, fHistogramList);

}

//====================================================================================================================================================

void AliAnalysisTaskMFTExample::Terminate(Option_t *) {

  // Draw result to the screen
  // Called once at the end of the query
  
}

//====================================================================================================================================================
