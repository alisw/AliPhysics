/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *									  *
 * Authors: Svein Lindal, Daniel Lohner					  *
 * Version 1.0								  *
 *									  *
 * Permission to use, copy, modify and distribute this software and its	  *
 * documentation strictly for non-commercial purposes is hereby granted	  *
 * without fee, provided that the above copyright notice appears in all	  *
 * copies and that both the copyright notice and this permission notice	  *
 * appear in the supporting documentation. The authors make no claims	  *
 * about the suitability of this software for any purpose. It is	  *
 * provided "as is" without express or implied warranty.		  *
 **************************************************************************/

////////////////////////////////////////////////
//---------------------------------------------
// QA Task for V0 Reader V1
//---------------------------------------------
////////////////////////////////////////////////

#include "AliAnalysisTaskConversionTree.h"
#include "TChain.h"
#include "TRandom.h"
#include "AliAnalysisManager.h"
#include "TParticle.h"
#include "TVectorF.h"
#include "AliPIDResponse.h"
#include "TFile.h"
#include "AliESDtrackCuts.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
#include "AliAODEvent.h"

class iostream;

using namespace std;

ClassImp(AliAnalysisTaskConversionTree)

//________________________________________________________________________
AliAnalysisTaskConversionTree::AliAnalysisTaskConversionTree() : AliAnalysisTaskSE(),
  fV0Reader(NULL),
  fV0ReaderName("V0ReaderV1"),
  fConversionGammas(NULL),
  fConversionCuts(NULL),
  fEventCuts(NULL),
  fInputEvent(NULL),
  fMCEvent(NULL),
  fAnalysisTree(NULL),
  fIsHeavyIon(kFALSE),
  ffillTree(-100),
  fOutputList(NULL),

  fBuffer_NConversionCandidates(0),
  fBuffer_ConversionCandidate_E(0),
  fBuffer_ConversionCandidate_Px(0),
  fBuffer_ConversionCandidate_Py(0),
  fBuffer_ConversionCandidate_Pz(0),
  fBuffer_ConversionCandidate_Qt(0),
  fBuffer_ConversionCandidate_Alpha(0),
  fBuffer_ConversionCandidate_PsiPair(0),
  fBuffer_ConversionCandidate_Chi2(0),
  fBuffer_ConversionCandidate_CosPA(0),
  fBuffer_ConversionCandidate_Eta(0),
  fBuffer_ConversionCandidate_Phi(0),
  fBuffer_ConversionCandidate_ConvPointX(0),
  fBuffer_ConversionCandidate_ConvPointY(0),
  fBuffer_ConversionCandidate_ConvPointZ(0),
  fBuffer_ConversionCandidate_MC_Type(0),
  fBuffer_ConversionCandidate_MC_Mother_ID(0),
  fBuffer_ConversionCandidate_MC_Mother_PDG(0),
  fBuffer_ConversionCandidate_MC_Mother_Type(0),
  fBuffer_ConversionCandidate_MC_Mother_TruePt(0),
  fBuffer_ConversionCandidate_MC_GrandMother_PDG(0),
  fIsMC(kFALSE),
  fnGammaCandidates(1),
  fMCStackPos(NULL),
  fMCStackNeg(NULL)
{
  fBuffer_ConversionCandidate_E               = new Float_t[kMaxConvCandidates];
  fBuffer_ConversionCandidate_Px              = new Float_t[kMaxConvCandidates];
  fBuffer_ConversionCandidate_Py              = new Float_t[kMaxConvCandidates];
  fBuffer_ConversionCandidate_Pz              = new Float_t[kMaxConvCandidates];
  fBuffer_ConversionCandidate_Qt              = new Float_t[kMaxConvCandidates];
  fBuffer_ConversionCandidate_Alpha           = new Float_t[kMaxConvCandidates];
  fBuffer_ConversionCandidate_PsiPair         = new Float_t[kMaxConvCandidates];
  fBuffer_ConversionCandidate_Chi2            = new Float_t[kMaxConvCandidates];
  fBuffer_ConversionCandidate_CosPA           = new Float_t[kMaxConvCandidates];
  fBuffer_ConversionCandidate_Eta             = new Float_t[kMaxConvCandidates];
  fBuffer_ConversionCandidate_Phi             = new Float_t[kMaxConvCandidates];
  fBuffer_ConversionCandidate_ConvPointX      = new Float_t[kMaxConvCandidates];
  fBuffer_ConversionCandidate_ConvPointY      = new Float_t[kMaxConvCandidates];
  fBuffer_ConversionCandidate_ConvPointZ      = new Float_t[kMaxConvCandidates];
  fBuffer_ConversionCandidate_MC_Type         = new UInt_t[kMaxConvCandidates];
  fBuffer_ConversionCandidate_MC_Mother_ID    = new Int_t[kMaxConvCandidates];
  fBuffer_ConversionCandidate_MC_Mother_PDG   = new Int_t[kMaxConvCandidates];
  fBuffer_ConversionCandidate_MC_Mother_Type  = new UInt_t[kMaxConvCandidates];
  fBuffer_ConversionCandidate_MC_Mother_TruePt= new Float_t[kMaxConvCandidates];
  fBuffer_ConversionCandidate_MC_GrandMother_PDG= new Int_t[kMaxConvCandidates];
}

AliAnalysisTaskConversionTree::AliAnalysisTaskConversionTree(const char *name) : AliAnalysisTaskSE(name),
  fV0Reader(NULL),
  fV0ReaderName("V0ReaderV1"),
  fConversionGammas(NULL),
  fConversionCuts(NULL),
  fEventCuts(NULL),
  fInputEvent(NULL),
  fMCEvent(NULL),
  fAnalysisTree(NULL),
  fIsHeavyIon(kFALSE),
  ffillTree(-100),
  fOutputList(NULL),

  fBuffer_NConversionCandidates(0),
  fBuffer_ConversionCandidate_E(0),
  fBuffer_ConversionCandidate_Px(0),
  fBuffer_ConversionCandidate_Py(0),
  fBuffer_ConversionCandidate_Pz(0),
  fBuffer_ConversionCandidate_Qt(0),
  fBuffer_ConversionCandidate_Alpha(0),
  fBuffer_ConversionCandidate_PsiPair(0),
  fBuffer_ConversionCandidate_Chi2(0),
  fBuffer_ConversionCandidate_CosPA(0),
  fBuffer_ConversionCandidate_Eta(0),
  fBuffer_ConversionCandidate_Phi(0),
  fBuffer_ConversionCandidate_ConvPointX(0),
  fBuffer_ConversionCandidate_ConvPointY(0),
  fBuffer_ConversionCandidate_ConvPointZ(0),
  fBuffer_ConversionCandidate_MC_Type(0),
  fBuffer_ConversionCandidate_MC_Mother_ID(0),
  fBuffer_ConversionCandidate_MC_Mother_PDG(0),
  fBuffer_ConversionCandidate_MC_Mother_Type(0),
  fBuffer_ConversionCandidate_MC_Mother_TruePt(0),
  fBuffer_ConversionCandidate_MC_GrandMother_PDG(0),
  fIsMC(kFALSE),
  fnGammaCandidates(1),
  fMCStackPos(NULL),
  fMCStackNeg(NULL)
{
  // Default constructor
  fBuffer_ConversionCandidate_E               = new Float_t[kMaxConvCandidates];
  fBuffer_ConversionCandidate_Px              = new Float_t[kMaxConvCandidates];
  fBuffer_ConversionCandidate_Py              = new Float_t[kMaxConvCandidates];
  fBuffer_ConversionCandidate_Pz              = new Float_t[kMaxConvCandidates];
  fBuffer_ConversionCandidate_Qt              = new Float_t[kMaxConvCandidates];
  fBuffer_ConversionCandidate_Alpha           = new Float_t[kMaxConvCandidates];
  fBuffer_ConversionCandidate_PsiPair         = new Float_t[kMaxConvCandidates];
  fBuffer_ConversionCandidate_Chi2            = new Float_t[kMaxConvCandidates];
  fBuffer_ConversionCandidate_CosPA           = new Float_t[kMaxConvCandidates];
  fBuffer_ConversionCandidate_Eta             = new Float_t[kMaxConvCandidates];
  fBuffer_ConversionCandidate_Phi             = new Float_t[kMaxConvCandidates];
  fBuffer_ConversionCandidate_ConvPointX      = new Float_t[kMaxConvCandidates];
  fBuffer_ConversionCandidate_ConvPointY      = new Float_t[kMaxConvCandidates];
  fBuffer_ConversionCandidate_ConvPointZ      = new Float_t[kMaxConvCandidates];
  fBuffer_ConversionCandidate_MC_Type         = new UInt_t[kMaxConvCandidates];
  fBuffer_ConversionCandidate_MC_Mother_ID    = new Int_t[kMaxConvCandidates];
  fBuffer_ConversionCandidate_MC_Mother_PDG   = new Int_t[kMaxConvCandidates];
  fBuffer_ConversionCandidate_MC_Mother_Type  = new UInt_t[kMaxConvCandidates];
  fBuffer_ConversionCandidate_MC_Mother_TruePt= new Float_t[kMaxConvCandidates];
  fBuffer_ConversionCandidate_MC_GrandMother_PDG= new Int_t[kMaxConvCandidates];

  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());
}

//________________________________________________________________________
AliAnalysisTaskConversionTree::~AliAnalysisTaskConversionTree()
{
  // default deconstructor
  
}
//________________________________________________________________________
void AliAnalysisTaskConversionTree::UserCreateOutputObjects()
{
  // Create User Output Objects

  if(fOutputList != NULL){
    delete fOutputList;
    fOutputList = NULL;
  }
  if(fOutputList == NULL){
    fOutputList = new TList();
    fOutputList->SetOwner(kTRUE);
  }


  if(ffillTree>=1.0){
    fAnalysisTree = new TTree(Form("PhotonTree_%s_%s",(fEventCuts->GetCutNumber()).Data(),(fConversionCuts->GetCutNumber()).Data()),Form("PhotonTree_%s_%s",(fEventCuts->GetCutNumber()).Data(),(fConversionCuts->GetCutNumber()).Data()));

    fAnalysisTree->Branch("NConversionCandidates",           &fBuffer_NConversionCandidates,             "NConversionCandidates/I"); // max 200 for now
    fAnalysisTree->Branch("ConversionCandidate_E",            fBuffer_ConversionCandidate_E,             "ConversionCandidate_E[NConversionCandidates]/F");
    fAnalysisTree->Branch("ConversionCandidate_Px",           fBuffer_ConversionCandidate_Px,            "ConversionCandidate_Px[NConversionCandidates]/F");
    fAnalysisTree->Branch("ConversionCandidate_Py",           fBuffer_ConversionCandidate_Py,            "ConversionCandidate_Py[NConversionCandidates]/F");
    fAnalysisTree->Branch("ConversionCandidate_Pz",           fBuffer_ConversionCandidate_Pz,            "ConversionCandidate_Pz[NConversionCandidates]/F");
    fAnalysisTree->Branch("ConversionCandidate_Qt",           fBuffer_ConversionCandidate_Qt,            "ConversionCandidate_Qt[NConversionCandidates]/F");
    fAnalysisTree->Branch("ConversionCandidate_Alpha",        fBuffer_ConversionCandidate_Alpha,         "ConversionCandidate_Alpha[NConversionCandidates]/F");
    fAnalysisTree->Branch("ConversionCandidate_PsiPair",      fBuffer_ConversionCandidate_PsiPair,       "ConversionCandidate_PsiPair[NConversionCandidates]/F");
    fAnalysisTree->Branch("ConversionCandidate_Chi2",         fBuffer_ConversionCandidate_Chi2,          "ConversionCandidate_Chi2[NConversionCandidates]/F");
    fAnalysisTree->Branch("ConversionCandidate_CosPA",        fBuffer_ConversionCandidate_CosPA,         "ConversionCandidate_CosPA[NConversionCandidates]/F");
    fAnalysisTree->Branch("ConversionCandidate_Eta",          fBuffer_ConversionCandidate_Eta,           "ConversionCandidate_Eta[NConversionCandidates]/F");
    fAnalysisTree->Branch("ConversionCandidate_Phi",          fBuffer_ConversionCandidate_Phi,           "ConversionCandidate_Phi[NConversionCandidates]/F");
    fAnalysisTree->Branch("ConversionCandidate_ConvPointX",   fBuffer_ConversionCandidate_ConvPointX,    "ConversionCandidate_ConvPointX[NConversionCandidates]/F");
    fAnalysisTree->Branch("ConversionCandidate_ConvPointY",   fBuffer_ConversionCandidate_ConvPointY,    "ConversionCandidate_ConvPointY[NConversionCandidates]/F");
    fAnalysisTree->Branch("ConversionCandidate_ConvPointZ",   fBuffer_ConversionCandidate_ConvPointZ,    "ConversionCandidate_ConvPointZ[NConversionCandidates]/F");

    if (fIsMC) {
      fAnalysisTree->Branch("ConversionCandidate_MC_Type",          fBuffer_ConversionCandidate_MC_Type,          "ConversionCandidate_MC_Type[NConversionCandidates]/b");
      fAnalysisTree->Branch("ConversionCandidate_MC_Mother_ID",     fBuffer_ConversionCandidate_MC_Mother_ID,     "ConversionCandidate_MC_Mother_ID[NConversionCandidates]/I");
      fAnalysisTree->Branch("ConversionCandidate_MC_Mother_PDG",    fBuffer_ConversionCandidate_MC_Mother_PDG,    "ConversionCandidate_MC_Mother_PDG[NConversionCandidates]/I");
      fAnalysisTree->Branch("ConversionCandidate_MC_Mother_Type",   fBuffer_ConversionCandidate_MC_Mother_Type,   "ConversionCandidate_MC_Mother_Type[NConversionCandidates]/b");
      fAnalysisTree->Branch("ConversionCandidate_MC_Mother_TruePt", fBuffer_ConversionCandidate_MC_Mother_TruePt, "ConversionCandidate_MC_Mother_TruePt[NConversionCandidates]/F");
      fAnalysisTree->Branch("ConversionCandidate_MC_GrandMother_PDG",    fBuffer_ConversionCandidate_MC_GrandMother_PDG,    "ConversionCandidate_MC_GrandMother_PDG[NConversionCandidates]/I");
    }
  }


  fV0Reader=(AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()->GetTask(fV0ReaderName.Data());
  
  if(fV0Reader && fV0Reader->GetProduceV0FindingEfficiency())
    if (fV0Reader->GetV0FindingEfficiencyHistograms())
      fOutputList->Add(fV0Reader->GetV0FindingEfficiencyHistograms());


  PostData(1, fOutputList);
  if(ffillTree>=1.0){
      OpenFile(2);
      PostData(2, fAnalysisTree);
  }
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskConversionTree::Notify()
{
    if (fEventCuts->GetPeriodEnum() == AliConvEventCuts::kNoPeriod && ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetPeriodEnum() != AliConvEventCuts::kNoPeriod){        
        fEventCuts->SetPeriodEnumExplicit(((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetPeriodEnum());
    } else if (fEventCuts->GetPeriodEnum() == AliConvEventCuts::kNoPeriod ){
      fEventCuts->SetPeriodEnum(fV0Reader->GetPeriodName());
    }  
 
  
  if(!fEventCuts->GetDoEtaShift()) return kTRUE; // No Eta Shift requested, continue
    
  if(fEventCuts->GetEtaShift() == 0.0){ // Eta Shift requested but not set, get shift automatically
    fEventCuts->GetCorrectEtaShiftFromPeriod();
    fEventCuts->DoEtaShift(kFALSE); // Eta Shift Set, make sure that it is called only once
    return kTRUE;
  }
  else{
    printf(" Gamma Conversion QA Task %s :: Eta Shift Manually Set to %f \n\n",
        (fEventCuts->GetCutNumber()).Data(),fEventCuts->GetEtaShift());
    fEventCuts->DoEtaShift(kFALSE); // Eta Shift Set, make sure that it is called only once
  }
  
  return kTRUE;
}
//________________________________________________________________________
void AliAnalysisTaskConversionTree::UserExec(Option_t *){

  Int_t eventQuality = ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetEventQuality();
  if(eventQuality != 0){// Event Not Accepted
    return;
  }
  fInputEvent = InputEvent();
  if(fIsMC) fMCEvent = MCEvent();

  Int_t eventNotAccepted =
    fEventCuts->IsEventAcceptedByCut(fV0Reader->GetEventCuts(),fInputEvent,fMCEvent,fIsHeavyIon,kFALSE);
  if(eventNotAccepted) return; // Check Centrality, PileUp, SDD and V0AND --> Not Accepted => eventQuality = 1

  fConversionGammas=fV0Reader->GetReconstructedGammas();

  // if a JetJet MC is used, check conversion products for their MC header origin
  if(fMCEvent){
    if(fEventCuts->GetSignalRejection() != 0){
      if(fInputEvent->IsA()==AliESDEvent::Class()){
        fEventCuts->GetNotRejectedParticles(fEventCuts->GetSignalRejection(),
                          fEventCuts->GetAcceptedHeader(),
                          fMCEvent);
      }
      else if(fInputEvent->IsA()==AliAODEvent::Class()){
        fEventCuts->GetNotRejectedParticles(fEventCuts->GetSignalRejection(),
                          fEventCuts->GetAcceptedHeader(),
                          fInputEvent);
      }
    }
  }

  if(fMCEvent && fInputEvent->IsA()==AliAODEvent::Class() && !(fV0Reader->AreAODsRelabeled())){
    RelabelAODPhotonCandidates(kTRUE);  // In case of AODMC relabeling MC
    fV0Reader->RelabelAODs(kTRUE);
  }

  // reduce event statistics in the tree by a factor ffilltree
  Bool_t ffillTreeNew = kFALSE;
  if(ffillTree>=1.0) {
    ffillTreeNew = kTRUE;
    if (ffillTree>1.0) {
      gRandom->SetSeed(0);
      if(gRandom->Uniform(ffillTree)>1.0) {
	      ffillTreeNew = kFALSE;
      }
    }
  }

  // reset tree buffer for each event
  ResetBuffer();

  // loop over conversion candidates
  for(Int_t firstGammaIndex=0;firstGammaIndex<fConversionGammas->GetEntriesFast();firstGammaIndex++){
    AliAODConversionPhoton *gamma=dynamic_cast<AliAODConversionPhoton*>(fConversionGammas->At(firstGammaIndex));
    if (gamma==NULL) continue;
    if(fMCEvent && fEventCuts->GetSignalRejection() != 0){
      if(!fEventCuts->IsParticleFromBGEvent(gamma->GetMCLabelPositive(), fMCEvent, fInputEvent))
        continue;
      if(!fEventCuts->IsParticleFromBGEvent(gamma->GetMCLabelNegative(), fMCEvent, fInputEvent))
        continue;
    }
    if(!fConversionCuts->PhotonIsSelected(gamma,fInputEvent)){
      continue;
    }

    if(ffillTreeNew) ProcessQATree(gamma);
  }
  // fill tree for each accepted event
  if (fAnalysisTree){
    fAnalysisTree->Fill();
  }

  if(fMCEvent && fInputEvent->IsA()==AliAODEvent::Class() && !(fV0Reader->AreAODsRelabeled())){
    RelabelAODPhotonCandidates(kFALSE); // Back to ESDMC Label
    fV0Reader->RelabelAODs(kFALSE);
  }

  PostData(1, fOutputList);
}


///________________________________________________________________________
void AliAnalysisTaskConversionTree::ProcessQATree(AliAODConversionPhoton *gamma){

  // Fill Histograms for QA and MC
  AliVEvent* event = (AliVEvent*) InputEvent();
    
  fBuffer_ConversionCandidate_E[fBuffer_NConversionCandidates] = gamma->GetPhotonP();
  fBuffer_ConversionCandidate_Px[fBuffer_NConversionCandidates] = gamma->GetPx();
  fBuffer_ConversionCandidate_Py[fBuffer_NConversionCandidates] = gamma->GetPy();
  fBuffer_ConversionCandidate_Pz[fBuffer_NConversionCandidates] = gamma->GetPz();
  fBuffer_ConversionCandidate_Qt[fBuffer_NConversionCandidates] = gamma->GetArmenterosQt();
  fBuffer_ConversionCandidate_Alpha[fBuffer_NConversionCandidates] = gamma->GetArmenterosAlpha();
  fBuffer_ConversionCandidate_PsiPair[fBuffer_NConversionCandidates] = gamma->GetPsiPair();
  fBuffer_ConversionCandidate_Chi2[fBuffer_NConversionCandidates] = gamma->GetChi2perNDF();
  fBuffer_ConversionCandidate_CosPA[fBuffer_NConversionCandidates] = fConversionCuts->GetCosineOfPointingAngle(gamma,event);
  fBuffer_ConversionCandidate_Eta[fBuffer_NConversionCandidates] = gamma->GetPhotonEta();
  fBuffer_ConversionCandidate_Phi[fBuffer_NConversionCandidates] = gamma->GetPhotonPhi();
  fBuffer_ConversionCandidate_ConvPointX[fBuffer_NConversionCandidates] = gamma->GetConversionX();
  fBuffer_ConversionCandidate_ConvPointY[fBuffer_NConversionCandidates] = gamma->GetConversionY();
  fBuffer_ConversionCandidate_ConvPointZ[fBuffer_NConversionCandidates] = gamma->GetConversionZ();

  fBuffer_ConversionCandidate_MC_Type[fBuffer_NConversionCandidates] = 9;
  fBuffer_ConversionCandidate_MC_Mother_Type[fBuffer_NConversionCandidates] = 9;
  if(fMCEvent && fInputEvent->IsA()==AliESDEvent::Class()){
    fBuffer_ConversionCandidate_MC_Type[fBuffer_NConversionCandidates] = IsTruePhotonESD(gamma);
    fBuffer_ConversionCandidate_MC_Mother_Type[fBuffer_NConversionCandidates] = GetTrueMotherInfoESD(gamma);
  } else if (fMCEvent && fInputEvent->IsA()==AliAODEvent::Class()){
    fBuffer_ConversionCandidate_MC_Type[fBuffer_NConversionCandidates] = IsTruePhotonAOD(gamma);
    fBuffer_ConversionCandidate_MC_Mother_Type[fBuffer_NConversionCandidates] = GetTrueMotherInfoAOD(gamma);
  }
  fBuffer_NConversionCandidates++;
}


UInt_t AliAnalysisTaskConversionTree::IsTruePhotonESD(AliAODConversionPhoton *TruePhotonCandidate)
{
  UInt_t kind = 9;
  TParticle *posDaughter = TruePhotonCandidate->GetPositiveMCDaughter(fMCEvent);
  TParticle *negDaughter = TruePhotonCandidate->GetNegativeMCDaughter(fMCEvent);
  Int_t pdgCodePos = 0; 
  Int_t pdgCodeNeg = 0; 
  Int_t pdgCode = 0; 

  const AliVVertex* primVtxMC 	= fMCEvent->GetPrimaryVertex();
  Double_t mcProdVtxX 	= primVtxMC->GetX();
  Double_t mcProdVtxY 	= primVtxMC->GetY();
  Double_t mcProdVtxZ 	= primVtxMC->GetZ();

  
  if(posDaughter == NULL || negDaughter == NULL) {
    kind = 9;
    //		return kFALSE; // One particle does not exist
  
  } else if( posDaughter->GetMother(0) != negDaughter->GetMother(0)  || (posDaughter->GetMother(0) == negDaughter->GetMother(0) && posDaughter->GetMother(0) ==-1)) {
    kind = 1;
    // 	  	return 1;
    pdgCodePos=TMath::Abs(posDaughter->GetPdgCode());
    pdgCodeNeg=TMath::Abs(negDaughter->GetPdgCode());
    if(pdgCodePos==11 && pdgCodeNeg==11) return 10; //Electron Combinatorial
    if(pdgCodePos==11 && pdgCodeNeg==11 && 
      (posDaughter->GetMother(0) == negDaughter->GetMother(0) && posDaughter->GetMother(0) ==-1)) return 15; //direct Electron Combinatorial
        
    if(pdgCodePos==211 && pdgCodeNeg==211) kind = 11; //Pion Combinatorial
    if((pdgCodePos==211 && pdgCodeNeg==2212) ||(pdgCodePos==2212 && pdgCodeNeg==211))	kind = 12; //Pion, Proton Combinatorics
    if((pdgCodePos==11 && pdgCodeNeg==2212) ||(pdgCodePos==2212 && pdgCodeNeg==11))	kind = 16; //electron, Proton Combinatorics
    if((pdgCodePos==11 && pdgCodeNeg==321) ||(pdgCodePos==321 && pdgCodeNeg==11))	kind = 17; //electron, kaon
    if((pdgCodePos==211 && pdgCodeNeg==321) ||(pdgCodePos==321 && pdgCodeNeg==211))	kind = 18; //pion, kaon
    if((pdgCodePos==211 && pdgCodeNeg==11) ||(pdgCodePos==11 && pdgCodeNeg==211)) kind = 13; //Pion, Electron Combinatorics
    if(pdgCodePos==321 && pdgCodeNeg==321) kind = 14; //Kaon,Kaon combinatorics
  }else{		
    pdgCodePos=posDaughter->GetPdgCode();
    pdgCodeNeg=negDaughter->GetPdgCode();
    Bool_t gammaIsPrimary = fEventCuts->IsConversionPrimaryESD( fMCEvent, posDaughter->GetMother(0), mcProdVtxX, mcProdVtxY, mcProdVtxZ);
    if ( TruePhotonCandidate->GetMCParticle(fMCEvent)->GetPdgCode()) pdgCode = TruePhotonCandidate->GetMCParticle(fMCEvent)->GetPdgCode();

    if(TMath::Abs(pdgCodePos)!=11 || TMath::Abs(pdgCodeNeg)!=11) return 2; // true from hadronic decays
    else if ( !(pdgCodeNeg==pdgCodePos)){
      if(pdgCode == 111) return 3; // pi0 Dalitz
      else if (pdgCode == 221) return 4; // eta Dalitz
      else if (!(negDaughter->GetUniqueID() != 5 || posDaughter->GetUniqueID() !=5)){
        if(pdgCode == 22 && gammaIsPrimary){
          return 0; // primary photons
        } else if (pdgCode == 22){
          return 5; //secondary photons
        }
      }
    }
  }

  return kind;
}

//________________________________________________________________________
UInt_t AliAnalysisTaskConversionTree::IsTruePhotonAOD(AliAODConversionPhoton *TruePhotonCandidate)
{   

  UInt_t kind = 9;
  TClonesArray *AODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
  if (AODMCTrackArray!=NULL && TruePhotonCandidate!=NULL){
    AliAODMCParticle *posDaughter = (AliAODMCParticle*) AODMCTrackArray->At(TruePhotonCandidate->GetMCLabelPositive());
    AliAODMCParticle *negDaughter = (AliAODMCParticle*) AODMCTrackArray->At(TruePhotonCandidate->GetMCLabelNegative());
    Int_t pdgCodePos = 0; 
    Int_t pdgCodeNeg = 0; 
    Int_t pdgCode = 0; 
    if(posDaughter == NULL || negDaughter == NULL) {
      kind = 9;
    } else if( posDaughter->GetMother() != negDaughter->GetMother()  || (posDaughter->GetMother() == negDaughter->GetMother() && posDaughter->GetMother() ==-1)) {
      kind = 1;
      pdgCodePos=TMath::Abs(posDaughter->GetPdgCode());
      pdgCodeNeg=TMath::Abs(negDaughter->GetPdgCode());
      if(pdgCodePos==11 && pdgCodeNeg==11)	kind = 10; //Electron Combinatorial
      if(pdgCodePos==11 && pdgCodeNeg==11 && 
        (posDaughter->GetMother() == negDaughter->GetMother() && posDaughter->GetMother() ==-1))kind = 15; //direct Electron Combinatorial

      if(pdgCodePos==211 && pdgCodeNeg==211) kind = 11; //Pion Combinatorial
      if((pdgCodePos==211 && pdgCodeNeg==2212) ||(pdgCodePos==2212 && pdgCodeNeg==211))	kind = 12; //Pion, Proton Combinatorics
      if((pdgCodePos==11 && pdgCodeNeg==2212) ||(pdgCodePos==2212 && pdgCodeNeg==11))	kind = 16; //electron, Proton Combinatorics
      if((pdgCodePos==11 && pdgCodeNeg==321) ||(pdgCodePos==321 && pdgCodeNeg==11))	kind = 17; //electron, kaon
      if((pdgCodePos==211 && pdgCodeNeg==321) ||(pdgCodePos==321 && pdgCodeNeg==211))	kind = 18; //pion, kaon
      if((pdgCodePos==211 && pdgCodeNeg==11) ||(pdgCodePos==11 && pdgCodeNeg==211)) kind = 13; //Pion, Electron Combinatorics
      if(pdgCodePos==321 && pdgCodeNeg==321) kind = 14; //Kaon,Kaon combinatorics
    }else{		
      AliAODMCParticle *Photon = (AliAODMCParticle*) AODMCTrackArray->At(posDaughter->GetMother());
      pdgCodePos=posDaughter->GetPdgCode();
      pdgCodeNeg=negDaughter->GetPdgCode();

      if ( Photon->GetPdgCode()) 
        pdgCode = Photon->GetPdgCode(); 
      if(TMath::Abs(pdgCodePos)!=11 || TMath::Abs(pdgCodeNeg)!=11) kind = 2; // true from hadronic decays
      else if ( !(pdgCodeNeg==pdgCodePos)){
        if(pdgCode == 111) kind = 3; // pi0 Dalitz
        else if (pdgCode == 221) kind = 4; // eta Dalitz
        else if (!(negDaughter->GetMCProcessCode() != 5 || posDaughter->GetMCProcessCode() !=5)){
          const AliVVertex* primVtxMC 	= fMCEvent->GetPrimaryVertex();
          Double_t mcProdVtxX 	= primVtxMC->GetX();
          Double_t mcProdVtxY 	= primVtxMC->GetY();
          Double_t mcProdVtxZ 	= primVtxMC->GetZ();
          Bool_t isPrimary = fEventCuts->IsConversionPrimaryAOD(fInputEvent, Photon, mcProdVtxX, mcProdVtxY, mcProdVtxZ);

          if(pdgCode == 22 && isPrimary){
            kind = 0; // primary photons
          } else if (pdgCode == 22){
            kind = 5; //secondary photons
          }
        }
      }
    }

    return kind;
  }	
  return kind;
}

UInt_t AliAnalysisTaskConversionTree::GetTrueMotherInfoESD(AliAODConversionPhoton *TruePhotonCandidate)
{
  UInt_t kind = 9;
  if (TruePhotonCandidate!=NULL){
    TParticle * negativeMC = (TParticle*)TruePhotonCandidate->GetNegativeMCDaughter(fMCEvent);
    TParticle * positiveMC = (TParticle*)TruePhotonCandidate->GetPositiveMCDaughter(fMCEvent);


    Int_t gammaMCLabel = -1;
    if(!positiveMC||!negativeMC)
      return kind;

    if(positiveMC->GetFirstMother()>-1&&(negativeMC->GetFirstMother() == positiveMC->GetFirstMother())){
      gammaMCLabel = positiveMC->GetFirstMother();
    }

    TParticle *photonCandMCParticle = fMCEvent->Particle(gammaMCLabel);
    Int_t pdgCodeMother = 0;

    if(photonCandMCParticle == NULL) {
      // particle does not exist
      kind = 9;
    } else if( photonCandMCParticle->GetFirstMother() ==-1) {
      // particle mother is no particle
      kind = 9;
    }else{
      // get mother and save ID of mother
      TParticle *photonMotherParticle = fMCEvent->Particle(photonCandMCParticle->GetFirstMother());
      fBuffer_ConversionCandidate_MC_Mother_ID[fBuffer_NConversionCandidates] = photonCandMCParticle->GetFirstMother();

      // get true Pt of mother particle
      fBuffer_ConversionCandidate_MC_Mother_TruePt[fBuffer_NConversionCandidates] = photonMotherParticle->Pt();

      // get PDG code of mother and save it
      pdgCodeMother = photonMotherParticle->GetPdgCode();
      fBuffer_ConversionCandidate_MC_Mother_PDG[fBuffer_NConversionCandidates] = pdgCodeMother;

      if(TMath::Abs(pdgCodeMother)==111){
        if(photonMotherParticle->GetFirstMother() != -1){
          TParticle *photonGrandMotherParticle = fMCEvent->Particle(photonMotherParticle->GetFirstMother());
          fBuffer_ConversionCandidate_MC_GrandMother_PDG[fBuffer_NConversionCandidates] = photonGrandMotherParticle->GetPdgCode();
          return 5;
        } else {
          return 0;
        }
      }
    }
  }
  return kind;
}

//________________________________________________________________________
UInt_t AliAnalysisTaskConversionTree::GetTrueMotherInfoAOD(AliAODConversionPhoton *TruePhotonCandidate)
{
  UInt_t kind = 9;
  TClonesArray *AODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
  if (AODMCTrackArray!=NULL && TruePhotonCandidate!=NULL){
    if(TruePhotonCandidate->GetMCLabelPositive() == -1 || TruePhotonCandidate->GetMCLabelNegative() == -1){
      return kind;
    }
    AliAODMCParticle *positiveMC = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TruePhotonCandidate->GetMCLabelPositive()));
    AliAODMCParticle *negativeMC = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TruePhotonCandidate->GetMCLabelNegative()));

    Int_t gammaMCLabel = -1;
    if(!positiveMC||!negativeMC)
      return kind;

    if(positiveMC->GetMother()>-1&&(negativeMC->GetMother() == positiveMC->GetMother())){
      gammaMCLabel = positiveMC->GetMother();
    }
    AliAODMCParticle *photonCandMCParticle = (AliAODMCParticle*) AODMCTrackArray->At(gammaMCLabel);
    Int_t pdgCodeMother = 0;

    if(photonCandMCParticle == NULL) {
      // particle does not exist
      kind = 9;
    } else if( photonCandMCParticle->GetMother() ==-1) {
      // particle mother is no particle
      kind = 9;
    }else{
      // get mother and save ID of mother
      AliAODMCParticle *photonMotherParticle = (AliAODMCParticle*) AODMCTrackArray->At(photonCandMCParticle->GetMother());
      if(photonMotherParticle != NULL){
        fBuffer_ConversionCandidate_MC_Mother_ID[fBuffer_NConversionCandidates] = photonCandMCParticle->GetMother();
        // get true Pt of mother particle
        fBuffer_ConversionCandidate_MC_Mother_TruePt[fBuffer_NConversionCandidates] = photonMotherParticle->Pt();

        // get PDG code of mother and save it
        pdgCodeMother = photonMotherParticle->GetPdgCode();
        fBuffer_ConversionCandidate_MC_Mother_PDG[fBuffer_NConversionCandidates] = pdgCodeMother;

        if(TMath::Abs(pdgCodeMother)==111){
          if(photonMotherParticle->GetMother() != -1){
            AliAODMCParticle *photonGrandMotherParticle = (AliAODMCParticle*) AODMCTrackArray->At(photonMotherParticle->GetMother());
            fBuffer_ConversionCandidate_MC_GrandMother_PDG[fBuffer_NConversionCandidates] = photonGrandMotherParticle->GetPdgCode();
            return 5;
          } else {
            return 0;
          }
        }
      }
    }
  }
  return kind;
}

//________________________________________________________________________
void AliAnalysisTaskConversionTree::RelabelAODPhotonCandidates(Bool_t mode){

  // Relabeling For AOD Event
  // ESDiD -> AODiD
  // MCLabel -> AODMCLabel
  
  if(mode){
    fMCStackPos = new Int_t[fConversionGammas->GetEntries()];
    fMCStackNeg = new Int_t[fConversionGammas->GetEntries()];
  }
  
  for(Int_t iGamma = 0;iGamma<fConversionGammas->GetEntries();iGamma++){
    AliAODConversionPhoton* PhotonCandidate = (AliAODConversionPhoton*) fConversionGammas->At(iGamma);
    if(!PhotonCandidate) continue;
    if(!mode){// Back to ESD Labels
      PhotonCandidate->SetMCLabelPositive(fMCStackPos[iGamma]);
      PhotonCandidate->SetMCLabelNegative(fMCStackNeg[iGamma]);
      //PhotonCandidate->IsAODMCLabel(kFALSE);
      continue;
    }
    fMCStackPos[iGamma] =  PhotonCandidate->GetMCLabelPositive();
    fMCStackNeg[iGamma] =  PhotonCandidate->GetMCLabelNegative();

    Bool_t AODLabelPos = kFALSE;
    Bool_t AODLabelNeg = kFALSE;

    for(Int_t i = 0; i<fInputEvent->GetNumberOfTracks();i++){
      AliAODTrack *tempDaughter = static_cast<AliAODTrack*>(fInputEvent->GetTrack(i));
      if(!AODLabelPos){
        if( tempDaughter->GetID() == PhotonCandidate->GetTrackLabelPositive() ){
        PhotonCandidate->SetMCLabelPositive(TMath::Abs(tempDaughter->GetLabel()));
        AODLabelPos = kTRUE;
        }
      }
      if(!AODLabelNeg){
        if( tempDaughter->GetID() == PhotonCandidate->GetTrackLabelNegative()){
        PhotonCandidate->SetMCLabelNegative(TMath::Abs(tempDaughter->GetLabel()));
        AODLabelNeg = kTRUE;
        }
      }
      if(AODLabelNeg && AODLabelPos){
        break;
      }
    } // Both ESD Tracks have AOD Tracks with Positive IDs
    if(!AODLabelPos || !AODLabelNeg){
      for(Int_t i = 0; i<fInputEvent->GetNumberOfTracks();i++){
        AliAODTrack *tempDaughter = static_cast<AliAODTrack*>(fInputEvent->GetTrack(i));
        if(tempDaughter->GetID()<0){
        if(!AODLabelPos){
          if( (TMath::Abs(tempDaughter->GetID())-1) == PhotonCandidate->GetTrackLabelPositive()){
            PhotonCandidate->SetMCLabelPositive(TMath::Abs(tempDaughter->GetLabel()));
            AODLabelPos = kTRUE;
          }
        }
        if(!AODLabelNeg){
          if( (TMath::Abs(tempDaughter->GetID())-1) == PhotonCandidate->GetTrackLabelNegative()){
            PhotonCandidate->SetMCLabelNegative(TMath::Abs(tempDaughter->GetLabel()));
            AODLabelNeg = kTRUE;
          }
        }
        }
        if(AODLabelNeg && AODLabelPos){
        break;
        }
      }
      if(!AODLabelPos || !AODLabelNeg){
        cout<<"WARNING!!! AOD TRACKS NOT FOUND FOR"<<endl;
        if(!AODLabelNeg){
          PhotonCandidate->SetMCLabelNegative(-999999);
          PhotonCandidate->SetLabelNegative(-999999);
        }
        if(!AODLabelPos){
          PhotonCandidate->SetMCLabelPositive(-999999);
          PhotonCandidate->SetLabelPositive(-999999);
        }
      }
    }
  }
  
  if(!mode){
    delete[] fMCStackPos;
    delete[] fMCStackNeg;
  }
}

//________________________________________________________________________
void AliAnalysisTaskConversionTree::Terminate(Option_t *)
{

}

void AliAnalysisTaskConversionTree::ResetBuffer(){

  fBuffer_NConversionCandidates = 0;

  for(Int_t ccand = 0; ccand < kMaxConvCandidates; ccand++){
    fBuffer_ConversionCandidate_E[ccand] = 0;
    fBuffer_ConversionCandidate_Px[ccand] = 0;
    fBuffer_ConversionCandidate_Py[ccand] = 0;
    fBuffer_ConversionCandidate_Pz[ccand] = 0;
    fBuffer_ConversionCandidate_Qt[ccand] = 0;
    fBuffer_ConversionCandidate_Alpha[ccand] = 0;
    fBuffer_ConversionCandidate_PsiPair[ccand] = 0;
    fBuffer_ConversionCandidate_Chi2[ccand] = 0;
    fBuffer_ConversionCandidate_CosPA[ccand] = 0;
    fBuffer_ConversionCandidate_Eta[ccand] = 0;
    fBuffer_ConversionCandidate_Phi[ccand] = 0;
    fBuffer_ConversionCandidate_ConvPointX[ccand] = 0;
    fBuffer_ConversionCandidate_ConvPointY[ccand] = 0;
    fBuffer_ConversionCandidate_ConvPointZ[ccand] = 0;
    if(fIsMC){
      fBuffer_ConversionCandidate_MC_Type[ccand] = 0;
      fBuffer_ConversionCandidate_MC_Mother_ID[ccand] = -1;
      fBuffer_ConversionCandidate_MC_Mother_PDG[ccand] = -1;
      fBuffer_ConversionCandidate_MC_Mother_Type[ccand] = 0;
      fBuffer_ConversionCandidate_MC_Mother_TruePt[ccand] = -1;
      fBuffer_ConversionCandidate_MC_GrandMother_PDG[ccand] = -1;
    }
  }
}