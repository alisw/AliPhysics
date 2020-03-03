/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//*************************************************************************
// \class AliHFTreeHandlerBstoDspi
// \brief helper class to handle a tree for Bs cut optimisation and MVA analyses
// \authors:
// L. Vermunt, luuk.vermunt@cern.ch
// G. Innocenti, gian.michele.innocenti@cern.ch
/////////////////////////////////////////////////////////////

#include <TString.h>
#include <TDatabasePDG.h>
#include "AliHFTreeHandlerBstoDspi.h"
#include "AliAODRecoDecayHF3Prong.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliESDVertex.h"
#include "AliESDtrackCuts.h"
#include "AliRDHFCutsDstoKKpi.h"

/// \cond CLASSIMP
ClassImp(AliHFTreeHandlerBstoDspi);
/// \endcond

//________________________________________________________________
AliHFTreeHandlerBstoDspi::AliHFTreeHandlerBstoDspi():
  AliHFTreeHandler(),
  fCosThetaStar(-9999.),
  fImpParProd(-9999.),
  fNormd0MeasMinusExp(-9999.),
  fInvMass_Ds(-9999.),
  fPt_Ds(-9999.),
  fY_Ds(-9999.),
  fEta_Ds(-9999.),
  fPhi_Ds(-9999.),
  fDecayLength_Ds(-9999.),
  fDecayLengthXY_Ds(-9999.),
  fNormDecayLengthXY_Ds(-9999.),
  fCosP_Ds(-9999.),
  fCosPXY_Ds(-9999.),
  fImpParXY_Ds(-9999.),
  fDCA_Ds(-9999.),
  fSigmaVertex_Ds(-9999.),
  fMassKK_Ds(-9999.),
  fCosPiDs_Ds(-9999.),
  fCosPiKPhi_Ds(-9999.),
  fNormd0MeasMinusExp_Ds(-9999.),
  fInvMassBsCut(0.3),
  fPtBsCut(-1),
  fImpParProdBsCut(9999.),
  fCosPBsCut(-9999.),
  fCosPXYBsCut(-9999.)
{
  //
  // Default constructor
  //

  //Only used in for-loops with a self-made AliAODtrack vector, so "Bs pion + 3 Ds-prongs" is fine
  fNProngs=4; // --> cannot be changed
  for(unsigned int iProng=0; iProng<fNProngs; iProng++) 
    fImpParProng[iProng] = -9999.;
}

//________________________________________________________________
AliHFTreeHandlerBstoDspi::AliHFTreeHandlerBstoDspi(int PIDopt):
  AliHFTreeHandler(PIDopt),
  fCosThetaStar(-9999.),
  fImpParProd(-9999.),
  fNormd0MeasMinusExp(-9999.),
  fInvMass_Ds(-9999.),
  fPt_Ds(-9999.),
  fY_Ds(-9999.),
  fEta_Ds(-9999.),
  fPhi_Ds(-9999.),
  fDecayLength_Ds(-9999.),
  fDecayLengthXY_Ds(-9999.),
  fNormDecayLengthXY_Ds(-9999.),
  fCosP_Ds(-9999.),
  fCosPXY_Ds(-9999.),
  fImpParXY_Ds(-9999.),
  fDCA_Ds(-9999.),
  fSigmaVertex_Ds(-9999.),
  fMassKK_Ds(-9999.),
  fCosPiDs_Ds(-9999.),
  fCosPiKPhi_Ds(-9999.),
  fNormd0MeasMinusExp_Ds(-9999.),
  fInvMassBsCut(0.3),
  fPtBsCut(-1),
  fImpParProdBsCut(9999.),
  fCosPBsCut(-9999.),
  fCosPXYBsCut(-9999.)
{
  //
  // Standard constructor
  //
    
  //Only used in for-loops with a self-made AliAODtrack vector, so "Bs pion + 3 Ds-prongs" is fine
  fNProngs=4; // --> cannot be changed
  for(unsigned int iProng=0; iProng<fNProngs; iProng++) 
    fImpParProng[iProng] = -9999.;
}

//________________________________________________________________
AliHFTreeHandlerBstoDspi::~AliHFTreeHandlerBstoDspi()
{
  //
  // Default Destructor
  //
}

//________________________________________________________________
TTree* AliHFTreeHandlerBstoDspi::BuildTree(TString name, TString title)
{
  fIsMCGenTree=false;
  
  if(fTreeVar) {
    delete fTreeVar;
    fTreeVar=nullptr;
  }
  fTreeVar = new TTree(name.Data(),title.Data());

  //set common variables
  AddCommonDmesonVarBranches();

  //set Bs variables
  fTreeVar->Branch("cos_t_star",&fCosThetaStar);
  fTreeVar->Branch("imp_par_prod",&fImpParProd);
  fTreeVar->Branch("max_norm_d0d0exp",&fNormd0MeasMinusExp);
  for(unsigned int iProng=0; iProng<fNProngs; iProng++){
    fTreeVar->Branch(Form("imp_par_prong%d",iProng),&fImpParProng[iProng]);
  }

  //set Ds variables
  fTreeVar->Branch("inv_mass_Ds",&fInvMass_Ds);
  fTreeVar->Branch("pt_Ds",&fPt_Ds);
  fTreeVar->Branch("y_Ds",&fY_Ds);
  fTreeVar->Branch("eta_Ds",&fEta_Ds);
  fTreeVar->Branch("phi_Ds",&fPhi_Ds);
  fTreeVar->Branch("d_len_Ds",&fDecayLength_Ds);
  fTreeVar->Branch("d_len_xy_Ds",&fDecayLengthXY_Ds);
  fTreeVar->Branch("norm_dl_xy_Ds",&fNormDecayLengthXY_Ds);
  fTreeVar->Branch("cos_p_Ds",&fCosP_Ds);
  fTreeVar->Branch("cos_p_xy_Ds",&fCosPXY_Ds);
  fTreeVar->Branch("imp_par_xy_Ds",&fImpParXY_Ds);
  fTreeVar->Branch("dca_Ds",&fDCA_Ds);
  fTreeVar->Branch("sig_vert_Ds",&fSigmaVertex_Ds);
  fTreeVar->Branch("delta_mass_KK_Ds",&fMassKK_Ds);
  fTreeVar->Branch("cos_PiDs_Ds",&fCosPiDs_Ds);
  fTreeVar->Branch("cos_PiKPhi_3_Ds",&fCosPiKPhi_Ds);
  fTreeVar->Branch("max_norm_d0d0exp_Ds",&fNormd0MeasMinusExp_Ds);

  //set single-track variables
  AddSingleTrackBranches();
  if (fFillJets) AddJetBranches();

  //set PID variables
  if(fPidOpt!=kNoPID) AddPidBranches(true,true,false,true,true);

  return fTreeVar;
}

//________________________________________________________________
bool AliHFTreeHandlerBstoDspi::SetVariables(int runnumber, int eventID, int eventID_Ext, Long64_t eventID_Long, float ptgen, AliAODRecoDecayHF* cand, float bfield, int masshypo/*used for Ds*/, AliPIDResponse* pidrespo)
{
  fIsMCGenTree=false;

  if(!cand) return false;
  if(fFillOnlySignal) { //if fill only signal and not signal candidate, do not store
    if(!(fCandType&kSignal || fCandType&kRefl)) return true;
  }
  fRunNumber=runnumber;
  fEvID=eventID;
  fEvIDExt=eventID_Ext;
  fEvIDLong=eventID_Long;
  fPtGen=ptgen;
  
  AliAODRecoDecayHF3Prong* candDs = (AliAODRecoDecayHF3Prong*)cand->GetDaughter(0); //Ds
  
  //topological variables
  //common (Bs -> Ds pi)
  fPt=((AliAODRecoDecayHF2Prong*)cand)->Pt();
  fY=((AliAODRecoDecayHF2Prong*)cand)->Y(531);
  fEta=((AliAODRecoDecayHF2Prong*)cand)->Eta();
  fPhi=((AliAODRecoDecayHF2Prong*)cand)->Phi();
  fDecayLength=((AliAODRecoDecayHF2Prong*)cand)->DecayLength();
  fDecayLengthXY=((AliAODRecoDecayHF2Prong*)cand)->DecayLengthXY();
  fNormDecayLengthXY=((AliAODRecoDecayHF2Prong*)cand)->NormalizedDecayLengthXY();
  fCosP=((AliAODRecoDecayHF2Prong*)cand)->CosPointingAngle();
  fCosPXY=((AliAODRecoDecayHF2Prong*)cand)->CosPointingAngleXY();
  fImpParXY=((AliAODRecoDecayHF2Prong*)cand)->ImpParXY();
  fDCA=((AliAODRecoDecayHF2Prong*)cand)->GetDCA();

  UInt_t prongs[2];
  prongs[0] = 431; prongs[1] = 211;
  fInvMass=((AliAODRecoDecayHF2Prong*)cand)->InvMass(2,prongs);
    
  fCosThetaStar=((AliAODRecoDecayHF2Prong*)cand)->CosThetaStar(0,531,421,211);
  fImpParProd=((AliAODRecoDecayHF2Prong*)cand)->Prodd0d0();
  fNormd0MeasMinusExp=ComputeMaxd0MeasMinusExp(cand,bfield);

  //Ds -> K K pi variables
  fPt_Ds=candDs->Pt();
  fY_Ds=candDs->Y(431);
  fEta_Ds=candDs->Eta();
  fPhi_Ds=candDs->Phi();
  fDecayLength_Ds=candDs->DecayLength();
  fDecayLengthXY_Ds=candDs->DecayLengthXY();
  fNormDecayLengthXY_Ds=candDs->NormalizedDecayLengthXY();
  fCosP_Ds=candDs->CosPointingAngle();
  fCosPXY_Ds=candDs->CosPointingAngleXY();
  fImpParXY_Ds=candDs->ImpParXY();
  fDCA_Ds=candDs->GetDCA();
  fNormd0MeasMinusExp_Ds=ComputeMaxd0MeasMinusExp(candDs,bfield);
  fSigmaVertex_Ds=((AliAODRecoDecayHF3Prong*)candDs)->GetSigmaVert();

  float cospikphi=-2;
  float massPhi = TDatabasePDG::Instance()->GetParticle(333)->Mass();
  if(masshypo==0){ //phiKKpi
    fInvMass_Ds=((AliAODRecoDecayHF3Prong*)candDs)->InvMassDsKKpi();
    fMassKK_Ds=TMath::Abs(((AliAODRecoDecayHF3Prong*)candDs)->InvMass2Prongs(0,1,321,321)-massPhi);
    fCosPiDs_Ds=((AliAODRecoDecayHF3Prong*)candDs)->CosPiDsLabFrameKKpi();
    cospikphi = ((AliAODRecoDecayHF3Prong*)candDs)->CosPiKPhiRFrameKKpi();
  }
  else if(masshypo==1){ //phipiKK
    fInvMass_Ds=((AliAODRecoDecayHF3Prong*)candDs)->InvMassDspiKK();
    fMassKK_Ds=TMath::Abs(((AliAODRecoDecayHF3Prong*)candDs)->InvMass2Prongs(1,2,321,321)-massPhi);
    fCosPiDs_Ds=((AliAODRecoDecayHF3Prong*)candDs)->CosPiDsLabFramepiKK();
    cospikphi = ((AliAODRecoDecayHF3Prong*)candDs)->CosPiKPhiRFramepiKK();
  }
  fCosPiKPhi_Ds=cospikphi*cospikphi*cospikphi;
  
  for(unsigned int iProng=0; iProng<3; iProng++) {
    fImpParProng[iProng]=candDs->Getd0Prong(iProng);
  }
  fImpParProng[3]=cand->Getd0Prong(1);

  AliAODTrack* prongtracks[4];
  for(unsigned int iProng=0; iProng<3; iProng++) prongtracks[iProng] = (AliAODTrack*)candDs->GetDaughter(iProng);
  prongtracks[3] = (AliAODTrack*)cand->GetDaughter(1);

  //single track variables
  bool setsingletrack = SetSingleTrackVars(prongtracks);
  if(!setsingletrack) return false;

  //pid variables
  if(fPidOpt==kNoPID) return true;

  bool setpid = SetPidVars(prongtracks,pidrespo,true,true,false,true,true);
  if(!setpid) return false;

  return true;
}

//________________________________________________________________
Int_t AliHFTreeHandlerBstoDspi::IsBsPionSelected(TObject* obj, AliRDHFCutsDstoKKpi* cutsDs, AliAODPidHF* fPidHFDs, AliAODEvent* aod, AliAODVertex *vtx) {

  AliAODTrack* candidatePion = (AliAODTrack*)obj;
  if (!candidatePion){ AliWarning("No pion object. Track rejected."); return 0; }
  
  AliESDtrackCuts* fTrackCuts = cutsDs->GetTrackCuts();
  if (!fTrackCuts){ AliWarning("No fTrackCuts object. Track rejected.");  return 0; }

  Double_t pos[3],cov[6];
  vtx->GetXYZ(pos);
  vtx->GetCovarianceMatrix(cov);
  const AliESDVertex vESD(pos,cov,100.,100);
  
  if(!cutsDs->IsDaughterSelected(candidatePion,&vESD,fTrackCuts,aod)) return 0;

  if(!cutsDs->GetIsUsePID()) return 1;
  else {
    
    if(!fPidHFDs){ AliWarning("AliAODPidHF not created. Track accepted"); return 1; }
    
    Int_t isPion=fPidHFDs->MakeRawPid(candidatePion,AliPID::kPion);
    if(isPion) return 1;
    else       return 0;
  }
  
  return 1;
}

//________________________________________________________________
Int_t AliHFTreeHandlerBstoDspi::IsBsSelected(AliAODRecoDecayHF2Prong* bs) {

  if (!bs){ AliWarning("No Bs AliAODRecoDecayHF2Prong object. Candidate rejected."); return 0; }

  UInt_t pdgDgBstoDspiUInt[2] = {431,211};
  Double_t invmassBs = bs->InvMass(2, pdgDgBstoDspiUInt);
  Double_t massBsPDG = TDatabasePDG::Instance()->GetParticle(531)->Mass();
  if(TMath::Abs(invmassBs-massBsPDG) > fInvMassBsCut) return 0;

  Double_t ptBs = bs->Pt();
  if(ptBs < fPtBsCut) return 0;

  Double_t impparprodBs = bs->Prodd0d0();
  if(impparprodBs > fImpParProdBsCut) return 0;

  Double_t cospBs = bs->CosPointingAngle();
  if(cospBs < fCosPBsCut) return 0;

  Double_t cospxyBs = bs->CosPointingAngleXY();
  if(cospxyBs < fCosPXYBsCut) return 0;

  return 1;
}

//________________________________________________________________
void AliHFTreeHandlerBstoDspi::SetDsBackgroundShapeType(bool isPr, bool isFDBplus, bool isFDB0, bool isFDLb0, bool isFDBs0) {

  if(isPr)      fCandType |= kDsPrompt;
  else          fCandType &= ~kDsPrompt;
  if(isFDBplus) fCandType |= kDsFDBplus;
  else          fCandType &= ~kDsFDBplus;
  if(isFDB0)    fCandType |= kDsFDB0;
  else          fCandType &= ~kDsFDB0;
  if(isFDLb0)   fCandType |= kDsFDLb0;
  else          fCandType &= ~kDsFDLb0;
  if(isFDBs0)   fCandType |= kDsFDBs0;
  else          fCandType &= ~kDsFDBs0;

}
