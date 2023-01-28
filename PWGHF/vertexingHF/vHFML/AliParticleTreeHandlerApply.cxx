/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//*************************************************************************
// \class AliParticleTreeHandlerApply
// \brief helper class to handle a tree for tracks with pion/kaon/proton hypotheses
// \authors:
// L. Vermunt, luuk.vermunt@cern.ch
/////////////////////////////////////////////////////////////

#include <TString.h>
#include <TDatabasePDG.h>
#include "AliPIDResponse.h"
#include "AliParticleTreeHandlerApply.h"

/// \cond CLASSIMP
ClassImp(AliParticleTreeHandlerApply);
/// \endcond

//________________________________________________________________
AliParticleTreeHandlerApply::AliParticleTreeHandlerApply():
  TObject(),
  fTreeVar(nullptr),
  fTrackSel(0),
  fPt(-9999.),
  fY(-9999.),
  fEta(-9999.),
  fPhi(-9999.),
  fPtGen(-9999.),
  fYGen(-9999.),
  fEtaGen(-9999.),
  fPhiGen(-9999.),
  fCharge(-9999),
  fID(-9999),
  fPDG(-9999),
  fMCLabel(-9999),
  fEvID(-9999),
  fEvIDExt(-9999),
  fEvIDLong(-9999),
  fRunNumber(-9999),
  fOnlyDedicatedBranches(false),
  fIsMC(false),
  fDebugMode(false),
  fTPCCls(-9999),
  fDCAXY(-9999),
  fDCAZ(-9999),
  fDCAXYProp(-9999),
  fDCAZProp(-9999),
  fTPCCrRows(-9999),
  fTPCClsToFnd(-9999),
  fChi2(-9999),
  fPTPC(-9999),
  fNSigTPC(-9999),
  fNSigTOF(-9999),
  fNSigCombTPCTOF(-9999)
{
  //
  // Default constructor
  //

}

//________________________________________________________________
AliParticleTreeHandlerApply::~AliParticleTreeHandlerApply()
{
  //
  // Destructor
  //

  delete fTreeVar;
}

//________________________________________________________________
TTree* AliParticleTreeHandlerApply::BuildTree(TString name, TString title)
{
  if(fTreeVar) {
    delete fTreeVar;
    fTreeVar=nullptr;
  }
  fTreeVar = new TTree(name.Data(),title.Data());

  fTreeVar->Branch("run_number", &fRunNumber);
  fTreeVar->Branch("ev_id",&fEvID);
  fTreeVar->Branch("ev_id_ext",&fEvIDExt);
  if(!fOnlyDedicatedBranches) fTreeVar->Branch("ev_id_long",&fEvIDLong);
  fTreeVar->Branch("track_sel",&fTrackSel);
  fTreeVar->Branch("pt",&fPt);
  fTreeVar->Branch("y",&fY);
  fTreeVar->Branch("eta",&fEta);
  fTreeVar->Branch("phi",&fPhi);
  fTreeVar->Branch("charge",&fCharge);
  fTreeVar->Branch("id",&fID);

  if(fIsMC){
    fTreeVar->Branch("pt_gen",&fPtGen);
    fTreeVar->Branch("y_gen",&fYGen);
    fTreeVar->Branch("eta_gen",&fEtaGen);
    fTreeVar->Branch("phi_gen",&fPhiGen);
    fTreeVar->Branch("pdg",&fPDG);
    fTreeVar->Branch("mc_label",&fMCLabel);
  }
  
  if(fDebugMode){
    fTreeVar->Branch("tpc_cls",&fTPCCls);
    fTreeVar->Branch("dca_xy",&fDCAXY);
    fTreeVar->Branch("dca_z",&fDCAZ);
    fTreeVar->Branch("dca_xy_prop",&fDCAXYProp);
    fTreeVar->Branch("dca_z_prop",&fDCAZProp);
    fTreeVar->Branch("tpc_cr_rows",&fTPCCrRows);
    fTreeVar->Branch("tpc_cls_to_fnd",&fTPCClsToFnd);
    fTreeVar->Branch("chi2",&fChi2);
    fTreeVar->Branch("p_tpc",&fPTPC);
    fTreeVar->Branch("nsig_TPC",&fNSigTPC);
    fTreeVar->Branch("nsig_TOF",&fNSigTOF);
    fTreeVar->Branch("nsig_comb_TPCTOF",&fNSigCombTPCTOF);
  }

  return fTreeVar;
}

//________________________________________________________________
bool AliParticleTreeHandlerApply::SetVariables(int runnumber, int eventID, int eventID_Ext, Long64_t eventID_Long, AliAODTrack* tr, AliAODTrack* trGlobal)
{
  if(!tr) return false;
  
  fRunNumber=runnumber;
  fEvID=eventID;
  fEvIDExt=eventID_Ext;
  fEvIDLong=eventID_Long;

  fPt=tr->Pt();
  fY=tr->Y();
  fEta=tr->Eta();
  fPhi=tr->Phi();
  fCharge=tr->Charge();
  fID=tr->GetID();

  fTPCCls=tr->GetTPCNcls();
  Float_t trDCAXYProp, trDCAZProp;
  trGlobal->GetImpactParameters(trDCAXYProp,trDCAZProp);
  fDCAXYProp=trDCAXYProp;
  fDCAZProp=trDCAZProp;
  fDCAXY=tr->DCA();
  fDCAZ=tr->ZAtDCA();
  fTPCCrRows=tr->GetTPCClusterInfo(2, 1);
  if (!(tr->GetTPCNclsF() > 0)) {
    fTPCClsToFnd = 0.;
  } else {
    fTPCClsToFnd = tr->GetTPCClusterInfo(2, 1) / float(tr->GetTPCNclsF());
  }
  fChi2=tr->Chi2perNDF();
  fPTPC=trGlobal->GetTPCmomentum();
    
  return true;
}

//________________________________________________________________
void AliParticleTreeHandlerApply::SetSelectionType(int part, bool isstd, bool ispidloose, bool ispidtight) {

  if(part == 0){
    if(isstd) fTrackSel |= kSelectedProton;
    else      fTrackSel &= ~kSelectedProton;
    if(ispidloose) fTrackSel |= kSelectedProtonPIDLoose;
    else           fTrackSel &= ~kSelectedProtonPIDLoose;
    if(ispidtight) fTrackSel |= kSelectedProtonPIDTight;
    else           fTrackSel &= ~kSelectedProtonPIDTight;
  } else if(part == 1){
    if(isstd) fTrackSel |= kSelectedKaon;
    else      fTrackSel &= ~kSelectedKaon;
    if(ispidloose) fTrackSel |= kSelectedKaonPIDLoose;
    else           fTrackSel &= ~kSelectedKaonPIDLoose;
    if(ispidtight) fTrackSel |= kSelectedKaonPIDTight;
    else           fTrackSel &= ~kSelectedKaonPIDTight;
  } else if(part == 2){
    if(isstd) fTrackSel |= kSelectedPion;
    else      fTrackSel &= ~kSelectedPion;
    if(ispidloose) fTrackSel |= kSelectedPionPIDLoose;
    else           fTrackSel &= ~kSelectedPionPIDLoose;
    if(ispidtight) fTrackSel |= kSelectedPionPIDTight;
    else           fTrackSel &= ~kSelectedPionPIDTight;
  } else if(part == 3){
    if(isstd) fTrackSel |= kSelectedNClsTPCStd;
    else      fTrackSel &= ~kSelectedNClsTPCStd;
    if(ispidtight) fTrackSel |= kSelectedNClsTPCTight;
    else           fTrackSel &= ~kSelectedNClsTPCTight;
  }
}

//________________________________________________________________
bool AliParticleTreeHandlerApply::SetMCGenVariables(AliAODMCParticle* mcpart, int mc_label) {

  fPtGen = -9999.;
  fYGen = -9999.;
  fEtaGen = -9999.;
  fPhiGen = -9999.;
  fPDG = -9999;
  fMCLabel = -9999;

  if(!mcpart) return false;
  
  fPtGen = mcpart->Pt();
  fYGen = mcpart->Y();
  fEtaGen = mcpart->Eta();
  fPhiGen = mcpart->Phi();
  fPDG = mcpart->GetPdgCode();
  fMCLabel = mc_label;

  return true;
}

//________________________________________________________________
bool AliParticleTreeHandlerApply::SetPIDVariables(float nsig_TPC, float nsig_TOF, float nsig_comb_TPCTOF) {

  fNSigTPC = nsig_TPC;
  fNSigTOF = nsig_TOF;
  fNSigCombTPCTOF = nsig_comb_TPCTOF;

  return true;
}
