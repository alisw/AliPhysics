/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Coiibutors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id: AliAnalysisTaskQuarkoniumTreeEmbedding.cxx $ */

//-----------------------------------------------------------------------------
// Analysis task to create a tree containing MC embedding information
// for both generated and reconstructed particles
// R. Arnaldi
//-----------------------------------------------------------------------------

// ROOT includes
#include "TROOT.h"
#include "TMath.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "TLorentzVector.h"

#include "AliInputEventHandler.h"
#include "AliAODHeader.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliAODMCParticle.h"
#include "AliMultSelection.h"
#include "AliMuonTrackCuts.h"

#include "AliAnalysisTaskQuarkoniumTreeEmbedding.h"

ClassImp(AliAnalysisTaskQuarkoniumTreeEmbedding)
//__________________________________________________________________________
AliAnalysisTaskQuarkoniumTreeEmbedding::AliAnalysisTaskQuarkoniumTreeEmbedding() :
  AliAnalysisTaskSE(),
  fBeamEnergy(0.),
  fkAnalysisType(0x0),
  fPeriod(0x0),
  fAODEvent(0x0),
  fOutputTree(0x0),
  fMuonTrackCuts(0x0),
  fResonance(0x0),
  fNMuons_gen(0x0),
  fNDimu_gen(0x0),
  fNMuons_rec(0x0),
  fPercentV0M(0x0)

{
  /// Default ctor.

  for(Int_t i=0; i<100;i++){
    fPt_rec[i]=999.;
    fE_rec[i]=999.;
    fPx_rec[i]=999;
    fPy_rec[i]=999;
    fPz_rec[i]=999;
    fY_rec[i]=999.;
    fEta_rec[i]=999.;
    fMatchTrig_rec[i]=999.;
    fTrackChi2_rec[i]=999.;
    fMatchTrigChi2_rec[i]=999.;
    fCharge_rec[i]=999;
    fRAtAbsEnd_rec[i]=999;
    fpDCA_rec[i] = 999.;
  }
  for(Int_t i=0; i<1000;i++){
    fDimuPt_gen[i]=999.;
    fDimuPx_gen[i]=999.;
    fDimuPy_gen[i]=999.;
    fDimuPz_gen[i]=999.;
    fDimuY_gen[i]=999.;
    fDimuMass_gen[i]=999.;
    fDimuCharge_gen[i]=999;

    fDimuPt_rec[i]=999.;
    fDimuPx_rec[i]=999.;
    fDimuPy_rec[i]=999.;
    fDimuPz_rec[i]=999.;
    fDimuY_rec[i]=999.;
    fDimuMass_rec[i]=999.;
    fDimuCharge_rec[i]=999;
    fDimuMatch_rec[i]=999;
    for(Int_t k=0;k<2;k++) fDimuMu_rec[i][k]=999;
  }

  fMuonTrackCuts = new AliMuonTrackCuts("StandardMuonTrackCuts", "StandardMuonTrackCuts");
  fMuonTrackCuts->SetFilterMask(AliMuonTrackCuts::kMuPdca);
  fMuonTrackCuts->SetAllowDefaultParams(kTRUE);
}

//__________________________________________________________________________
AliAnalysisTaskQuarkoniumTreeEmbedding::AliAnalysisTaskQuarkoniumTreeEmbedding(const char *name) :
  AliAnalysisTaskSE(name),
  fBeamEnergy(0.),
  fkAnalysisType(0x0),
  fPeriod(0x0),
  fAODEvent(0x0),
  fOutputTree(0x0),
  fMuonTrackCuts(0x0),
  fResonance(0x0),
  fNMuons_gen(0x0),
  fNDimu_gen(0x0),
  fNMuons_rec(0x0),
  fPercentV0M(0x0)
{
  //
  // Constructor. Initialization of Inputs and Outputs
  //
  Info("AliAnalysisTaskQuarkoniumTreeEmbedding","Calling Constructor");

   for(Int_t i=0; i<100;i++){
     fPt_rec[i]=999.;
     fE_rec[i]=999.;
     fPx_rec[i]=999;
     fPy_rec[i]=999;
     fPz_rec[i]=999;
     fY_rec[i]=999.;
     fEta_rec[i]=999.;
     fMatchTrig_rec[i]=999.;
     fTrackChi2_rec[i]=999.;
     fMatchTrigChi2_rec[i]=999.;
     fCharge_rec[i]=999;
     fRAtAbsEnd_rec[i]=999;
     fpDCA_rec[i] = 999.;
   }
   for(Int_t i=0; i<1000;i++){
     fDimuPt_gen[i]=999.;
     fDimuPx_gen[i]=999.;
     fDimuPy_gen[i]=999.;
     fDimuPz_gen[i]=999.;
     fDimuY_gen[i]=999.;
     fDimuMass_gen[i]=999.;
     fDimuCharge_gen[i]=999;

     fDimuPt_rec[i]=999.;
     fDimuPx_rec[i]=999.;
     fDimuPy_rec[i]=999.;
     fDimuPz_rec[i]=999.;
     fDimuY_rec[i]=999.;
     fDimuMass_rec[i]=999.;
     fDimuCharge_rec[i]=999;
     fDimuMatch_rec[i]=999;
     for(Int_t k=0;k<2;k++) fDimuMu_rec[i][k]=999;
   }

  fMuonTrackCuts = new AliMuonTrackCuts("StandardMuonTrackCuts", "StandardMuonTrackCuts");
  fMuonTrackCuts->SetFilterMask(AliMuonTrackCuts::kMuPdca);
  fMuonTrackCuts->SetAllowDefaultParams(kTRUE);

  DefineOutput(1,TTree::Class());
}

//___________________________________________________________________________
AliAnalysisTaskQuarkoniumTreeEmbedding& AliAnalysisTaskQuarkoniumTreeEmbedding::operator=(const AliAnalysisTaskQuarkoniumTreeEmbedding& c)
{
  //
  // Assignment operator
  //
  if (this!=&c) {
    AliAnalysisTaskSE::operator=(c) ;
  }
  return *this;
}

//___________________________________________________________________________
AliAnalysisTaskQuarkoniumTreeEmbedding::AliAnalysisTaskQuarkoniumTreeEmbedding(const AliAnalysisTaskQuarkoniumTreeEmbedding& c) :
  AliAnalysisTaskSE(c),
  fBeamEnergy(c.fBeamEnergy),
  fkAnalysisType(c.fkAnalysisType),
  fPeriod(c.fPeriod),
  fAODEvent(c.fAODEvent),
  fOutputTree(c.fOutputTree),
  fMuonTrackCuts(c.fMuonTrackCuts),
  fResonance(c.fResonance),
  fNMuons_gen(c.fNMuons_gen),
  fNMuons_rec(c.fNMuons_rec),
  fNDimu_rec(c.fNDimu_rec),
  fPercentV0M(c.fPercentV0M)
 {
  //
  // Copy Constructor
  //
 }

//___________________________________________________________________________
AliAnalysisTaskQuarkoniumTreeEmbedding::~AliAnalysisTaskQuarkoniumTreeEmbedding() {
  //
  //destructor
  //
  Info("~AliAnalysisTaskQuarkoniumTreeEmbedding","Calling Destructor");
  if (AliAnalysisManager::GetAnalysisManager()->GetAnalysisType() != AliAnalysisManager::kProofAnalysis) delete fOutputTree;
}

//___________________________________________________________________________
void AliAnalysisTaskQuarkoniumTreeEmbedding::NotifyRun()
{
  fMuonTrackCuts->SetRun(fInputHandler);
}

//___________________________________________________________________________
void AliAnalysisTaskQuarkoniumTreeEmbedding::UserCreateOutputObjects(){

 if (fOutputTree) return;

 OpenFile(1,"RECREATE");
 fOutputTree = new TTree("MCTree","Data Tree");

 fOutputTree->Branch("PercentV0M",&fPercentV0M,"PercentV0M/D");
 fOutputTree->Branch("NMuons_gen",&fNMuons_gen,"NMuons_gen/I");
 fOutputTree->Branch("NDimu_gen",&fNDimu_gen,"NDimu_gen/I");
 fOutputTree->Branch("DimuPt_gen",fDimuPt_gen,"DimuPt_gen[NDimu_gen]/D");
 fOutputTree->Branch("DimuPx_gen",fDimuPx_gen,"DimuPx_gen[NDimu_gen]/D");
 fOutputTree->Branch("DimuPy_gen",fDimuPy_gen,"DimuPy_gen[NDimu_gen]/D");
 fOutputTree->Branch("DimuPz_gen",fDimuPz_gen,"DimuPz_gen[NDimu_gen]/D");
 fOutputTree->Branch("DimuY_gen",fDimuY_gen,"DimuY_gen[NDimu_gen]/D");
 fOutputTree->Branch("DimuMass_gen",fDimuMass_gen,"DimuMass_gen[NDimu_gen]/D");
 fOutputTree->Branch("DimuCharge_gen",fDimuCharge_gen,"DimuCharge_gen[NDimu_gen]/I");

 fOutputTree->Branch("NMuons_rec",&fNMuons_rec,"NMuons_rec/I");
 fOutputTree->Branch("Pt_rec",fPt_rec,"Pt_rec[NMuons_rec]/D");
 fOutputTree->Branch("E_rec",fE_rec,"E_rec[NMuons_rec]/D");
 fOutputTree->Branch("Px_rec",fPx_rec,"Px_rec[NMuons_rec]/D");
 fOutputTree->Branch("Py_rec",fPy_rec,"Py_rec[NMuons_rec]/D");
 fOutputTree->Branch("Pz_rec",fPz_rec,"Pz_rec[NMuons_rec]/D");
 fOutputTree->Branch("Y_rec",fY_rec,"Y_rec[NMuons_rec]/D");
 fOutputTree->Branch("Eta_rec",fEta_rec,"Eta_rec[NMuons_rec]/D");
 fOutputTree->Branch("MatchTrig_rec",fMatchTrig_rec,"MatchTrig_rec[NMuons_rec]/I");
 fOutputTree->Branch("TrackChi2_rec",fTrackChi2_rec,"TrackChi2_rec[NMuons_rec]/D");
 fOutputTree->Branch("MatchTrigChi2_rec",fMatchTrigChi2_rec,"MatchTrigChi2_rec[NMuons_rec]/D");
 fOutputTree->Branch("Charge_rec",fCharge_rec,"Charge_rec[NMuons_rec]/I");
 fOutputTree->Branch("RAtAbsEnd_rec",fRAtAbsEnd_rec,"RAtAbsEnd_rec[NMuons_rec]/D");
 fOutputTree->Branch("pDCA_rec",fpDCA_rec,"pDCA[NMuons_rec]/I");

 fOutputTree->Branch("NDimu_rec",&fNDimu_rec,"NDimu_rec/I");
 fOutputTree->Branch("DimuMu_rec",fDimuMu_rec,"DimuMu_rec[NDimu_rec][2]/I");
 fOutputTree->Branch("DimuPt_rec",fDimuPt_rec,"DimuPt_rec[NDimu_rec]/D");
 fOutputTree->Branch("DimuPx_rec",fDimuPx_rec,"DimuPx_rec[NDimu_rec]/D");
 fOutputTree->Branch("DimuPy_rec",fDimuPy_rec,"DimuPy_rec[NDimu_rec]/D");
 fOutputTree->Branch("DimuPz_rec",fDimuPz_rec,"DimuPz_rec[NDimu_rec]/D");
 fOutputTree->Branch("DimuY_rec",fDimuY_rec,"DimuY_rec[NDimu_rec]/D");
 fOutputTree->Branch("DimuMass_rec",fDimuMass_rec,"DimuMass_rec[NDimu_rec]/D");
 fOutputTree->Branch("DimuCharge_rec",fDimuCharge_rec,"DimuCharge_rec[NDimu_rec]/I");
 fOutputTree->Branch("DimuMatch_rec",fDimuMatch_rec,"DimuMatch_rec[NDimu_rec]/I");

 fOutputTree->ls();

 PostData(1,fOutputTree);
}

//_________________________________________________
void AliAnalysisTaskQuarkoniumTreeEmbedding::UserExec(Option_t *)
{
  fNMuons_gen=0;
  fNDimu_gen=0;
  fNMuons_rec=0;
  fNDimu_rec=0;
  for(Int_t i=0; i<100;i++){
    fPt_rec[i]=999.;
    fE_rec[i]=999.;
    fPx_rec[i]=999;
    fPy_rec[i]=999;
    fPz_rec[i]=999;
    fY_rec[i]=999.;
    fEta_rec[i]=999.;
    fMatchTrig_rec[i]=999.;
    fTrackChi2_rec[i]=999.;
    fMatchTrigChi2_rec[i]=999.;
    fCharge_rec[i]=999;
    fRAtAbsEnd_rec[i]=999;
    fpDCA_rec[i]=999;
  }
  for(Int_t i=0; i<1000;i++){
    fDimuPt_gen[i]=999.;
    fDimuPx_gen[i]=999.;
    fDimuPy_gen[i]=999.;
    fDimuPz_gen[i]=999.;
    fDimuY_gen[i]=999.;
    fDimuMass_gen[i]=999.;
    fDimuCharge_gen[i]=999.;

    fDimuPt_rec[i]=999.;
    fDimuPx_rec[i]=999.;
    fDimuPy_rec[i]=999.;
    fDimuPz_rec[i]=999.;
    fDimuY_rec[i]=999.;
    fDimuMass_rec[i]=999.;
    fDimuCharge_rec[i]=999.;
    for(Int_t k=0;k<2;k++) fDimuMu_rec[i][k]=999;
  }

//---------------------------------------------------
// Execute analysis for current event
//---------------------------------------------------
  fAODEvent = dynamic_cast<AliAODEvent*> (InputEvent());
  if ( ! fAODEvent ) {
    AliError ("AOD event not found. Nothing done!");
    return;
  }

  //-----------------------------------------------
  // loop on MC generated events
  //-----------------------------------------------
  TClonesArray *mcarray = dynamic_cast<TClonesArray*>(fAODEvent->FindListObject(AliAODMCParticle::StdBranchName()));

  AliAODHeader *aodheader=dynamic_cast<AliAODHeader*>(fAODEvent->GetHeader());
  AliMultSelection *MultSelection = (AliMultSelection * ) fAODEvent->FindListObject("MultSelection");
  fPercentV0M = (Double_t) MultSelection->GetMultiplicityPercentile("V0M");

  //check resonance type
  Int_t PDGCode;
  if(fResonance == "JPsi") PDGCode = 443;
  if(fResonance == "Psi2S") PDGCode = 100443;

  Int_t numdimu_gen=0;
  for (Int_t i=0;i<mcarray->GetEntries();i++){
    AliAODMCParticle *mcp = (AliAODMCParticle *)mcarray->At(i);
    if(mcp->GetPdgCode()==PDGCode){
      fDimuPt_gen[i]=mcp->Pt();
      fDimuPx_gen[i]=mcp->Px();
      fDimuPy_gen[i]=mcp->Py();
      fDimuPz_gen[i]=mcp->Pz();
      fDimuY_gen[i]=mcp->Y();
      fDimuMass_gen[i]=mcp->M();
      fDimuCharge_gen[i]=mcp->Charge();
      numdimu_gen++;
    }
  }
  fNDimu_gen=numdimu_gen;

//-----------------------------------------------
// loop on reconstructed particles
//-----------------------------------------------
Int_t numtracks = fAODEvent->GetNumberOfTracks();
Int_t ndimu=0;
Int_t nmu=0;
Int_t Label0[numtracks];
Int_t Label1[numtracks];
Bool_t GoodTrack[numtracks];

for(int i=0;i<numtracks;i++){
  Label0[i]=999;
  Label1[i]=999;
  GoodTrack[i]=kFALSE;
}


for(int i=0;i<numtracks;i++) {
  AliAODTrack *mu0=(AliAODTrack*)fAODEvent->GetTrack(i);
  if(mu0->GetLabel()== -1) continue;
  AliAODMCParticle *mctrack0 = (AliAODMCParticle*) mcarray->At(mu0->GetLabel());
  if(mu0->GetLabel() != mctrack0->GetLabel()) continue;
  if(TMath::Abs(mctrack0->GetPdgCode())!=13) continue;
  Int_t mum_num0 = mctrack0->GetMother();
  if(mum_num0 == -1) continue;
  AliAODMCParticle *mcpart0 = (AliAODMCParticle *)mcarray->At(mum_num0);
  if(mcpart0->GetPdgCode()!=PDGCode) continue;
  GoodTrack[i]=kTRUE;

  for(int j=i+1;j<numtracks;j++){
    AliAODTrack *mu1=(AliAODTrack*)fAODEvent->GetTrack(j);
    if(mu1->GetLabel()== -1) continue;
    AliAODMCParticle *mctrack1 = (AliAODMCParticle*) mcarray->At(mu1->GetLabel());
    if(mu1->GetLabel() != mctrack1->GetLabel()) continue;
    if(TMath::Abs(mctrack1->GetPdgCode())!=13) continue;
    Int_t mum_num1 = mctrack1->GetMother();
    if(mum_num1 == -1) continue;
    AliAODMCParticle *mcpart1 = (AliAODMCParticle *)mcarray->At(mum_num1);
    if(mcpart1->GetPdgCode()!=PDGCode) continue;
    GoodTrack[j]=kTRUE;

    AliAODDimuon *dimu = new AliAODDimuon(mu0,mu1);
    fDimuMass_rec[ndimu] = dimu->Mass();
    fDimuPt_rec[ndimu] = dimu->Pt();
    fDimuPx_rec[ndimu] = dimu->Px();
    fDimuPy_rec[ndimu] = dimu->Py();
    fDimuPz_rec[ndimu] = dimu->Pz();
    fDimuY_rec[ndimu] = dimu->Y();
    fDimuCharge_rec[ndimu]= dimu->Charge();

    if(mu0->GetMatchTrigger()>1 || mu1->GetMatchTrigger()>1) fDimuMatch_rec[ndimu]=1;
    if(mu0->GetMatchTrigger()>1 && mu1->GetMatchTrigger()>1) fDimuMatch_rec[ndimu]=2;

    Label0[ndimu]=i;
    Label1[ndimu]=j;

    fDimuMu_rec[ndimu][0]=i;
    fDimuMu_rec[ndimu][1]=j;

    delete dimu;
    ndimu++;
    }
  }
// save muon tracks surviving cuts

 for(int i=0;i<numtracks;i++) {
    AliAODTrack *mu0=(AliAODTrack*)fAODEvent->GetTrack(i);
    if(GoodTrack[i]){
        fPt_rec[nmu] = mu0->Pt();
        fPx_rec[nmu] = mu0->Px();
        fPy_rec[nmu] = mu0->Py();
        fPz_rec[nmu] = mu0->Pz();
        fY_rec[nmu]  = mu0->Y();
        fEta_rec[nmu]= mu0->Eta();
        fE_rec[nmu] = mu0->E();
        fMatchTrig_rec[nmu]   = mu0->GetMatchTrigger();
        fMatchTrigChi2_rec[nmu]= mu0->GetChi2MatchTrigger();
        fRAtAbsEnd_rec[nmu] = mu0->GetRAtAbsorberEnd();
        if(fMuonTrackCuts->IsSelected(mu0)) fpDCA_rec[nmu] = 1;
      }
      nmu++;
    }

     fNMuons_rec = nmu;
     fNDimu_rec = ndimu;

    fOutputTree->Fill();
    PostData(1,fOutputTree);
}

//________________________________________________________________________
void AliAnalysisTaskQuarkoniumTreeEmbedding::Terminate(Option_t *)
{
}
