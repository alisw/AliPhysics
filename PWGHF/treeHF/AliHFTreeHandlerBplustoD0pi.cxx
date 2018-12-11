/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//*************************************************************************
// \class AliHFTreeHandlerBplustoD0pi
// \brief helper class to handle a tree for B+ cut optimisation and MVA analyses
// \authors:
// F. Catalano, fabio.catalano@cern.ch
// A. Festanti, andrea.festanti@cern.ch
// F. Grosa, fabrizio.grosa@cern.ch
// G. Innocenti, gian.michele.innocenti@cern.ch
// F. Prino, prino@to.infn.it
// L. Vermunt, luuk.vermunt@cern.ch
// L. van Doremalen, lennart.van.doremalen@cern.ch
/////////////////////////////////////////////////////////////

#include <TString.h>
#include "AliHFTreeHandlerBplustoD0pi.h"
#include "AliAODRecoDecayHF2Prong.h"

/// \cond CLASSIMP
ClassImp(AliHFTreeHandlerBplustoD0pi);
/// \endcond

//________________________________________________________________
AliHFTreeHandlerBplustoD0pi::AliHFTreeHandlerBplustoD0pi():
  AliHFTreeHandler(),
  fImpParProng(),
  fCosThetaStar(),
  fImpParProd(),
  fNormd0MeasMinusExp(),
  fDCA(),
  fAngleProngs(),
  fInvMass_D0(),
  fPt_D0(),
  fY_D0(),
  fEta_D0(),
  fPhi_D0(),
  fDecayLength_D0(),
  fDecayLengthXY_D0(),
  fNormDecayLengthXY_D0(),
  fCosP_D0(),
  fCosPXY_D0(),
  fImpParXY_D0(),
  fCosThetaStar_D0(),
  fImpParProd_D0(),
  fNormd0MeasMinusExp_D0(),
  fDCA_D0(),
  fAngleProngs_D0()
{
  //
  // Default constructor
  //

  //Only used in for-loops with a self-made AliAODtrack vector, so "Bplus pion + 2 D0-prongs" is fine
  fNProngs=3; // --> cannot be changed
}

//________________________________________________________________
AliHFTreeHandlerBplustoD0pi::AliHFTreeHandlerBplustoD0pi(int PIDopt):
  AliHFTreeHandler(PIDopt),
  fImpParProng(),
  fCosThetaStar(),
  fImpParProd(),
  fNormd0MeasMinusExp(),
  fDCA(),
  fAngleProngs(),
  fInvMass_D0(),
  fPt_D0(),
  fY_D0(),
  fEta_D0(),
  fPhi_D0(),
  fDecayLength_D0(),
  fDecayLengthXY_D0(),
  fNormDecayLengthXY_D0(),
  fCosP_D0(),
  fCosPXY_D0(),
  fImpParXY_D0(),
  fCosThetaStar_D0(),
  fImpParProd_D0(),
  fNormd0MeasMinusExp_D0(),
  fDCA_D0(),
  fAngleProngs_D0()
{
  //
  // Standard constructor
  //
    
  //Only used in for-loops with a self-made AliAODtrack vector, so "Bplus pion + 2 D0-prongs" is fine
  fNProngs=3; // --> cannot be changed
}

//________________________________________________________________
AliHFTreeHandlerBplustoD0pi::~AliHFTreeHandlerBplustoD0pi()
{
  //
  // Default Destructor
  //
}

//________________________________________________________________
TTree* AliHFTreeHandlerBplustoD0pi::BuildTree(TString name, TString title) 
{
  fIsMCGenTree=false;
  
  if(fTreeVar) {
    delete fTreeVar;
    fTreeVar=0x0;
  }
  fTreeVar = new TTree(name.Data(),title.Data());

  //set common variables
  AddCommonDmesonVarBranches();

  //set Bplus variables
  fTreeVar->Branch("cos_t_star",&fCosThetaStar);
  fTreeVar->Branch("imp_par_prod",&fImpParProd);
  fTreeVar->Branch("max_norm_d0d0exp",&fNormd0MeasMinusExp);
  for(unsigned int iProng=0; iProng<fNProngs; iProng++){
    fTreeVar->Branch(Form("imp_par_prong%d",iProng),&fImpParProng[iProng]);
  }
  fTreeVar->Branch("dca",&fDCA);
  fTreeVar->Branch("angle_prongs",&fAngleProngs);

  //set D0 variables
  fTreeVar->Branch("inv_mass_D0",&fInvMass_D0);
  fTreeVar->Branch("pt_D0",&fPt_D0);
  fTreeVar->Branch("y_D0",&fY_D0);
  fTreeVar->Branch("eta_D0",&fEta_D0);
  fTreeVar->Branch("phi_D0",&fPhi_D0);
  fTreeVar->Branch("d_len_D0",&fDecayLength_D0);
  fTreeVar->Branch("d_len_xy_D0",&fDecayLengthXY_D0);
  fTreeVar->Branch("norm_dl_xy_D0",&fNormDecayLengthXY_D0);
  fTreeVar->Branch("cos_p_D0",&fCosP_D0);
  fTreeVar->Branch("cos_p_xy_D0",&fCosPXY_D0);
  fTreeVar->Branch("imp_par_xy_D0",&fImpParXY_D0);
  fTreeVar->Branch("cos_t_star_D0",&fCosThetaStar_D0);
  fTreeVar->Branch("imp_par_prod_D0",&fImpParProd_D0);
  fTreeVar->Branch("max_norm_d0d0exp_D0",&fNormd0MeasMinusExp_D0);
  fTreeVar->Branch("dca_D0",&fDCA_D0);
  fTreeVar->Branch("angle_prongs_D0",&fAngleProngs_D0);
    
  //set single-track variables
  AddSingleTrackBranches();

  //set PID variables
  if(fPidOpt!=kNoPID) AddPidBranches(true,true,false,true,true);

  return fTreeVar;
}

//________________________________________________________________
bool AliHFTreeHandlerBplustoD0pi::SetVariables(AliAODRecoDecayHF* cand, float bfield, int /*masshypo*/, AliAODPidHF* pidHF)
{
  fIsMCGenTree=false;

  if(!cand) return false;
  if(fFillOnlySignal) { //if fill only signal and not signal candidate, do not store
    if(!(fCandTypeMap&kSignal || fCandTypeMap&kRefl)) return true;
  }
  fNCandidates++;

  AliAODTrack* cand_pr0 = (AliAODTrack*)cand->GetDaughter(0); //Bplus pion
  AliAODRecoDecayHF2Prong* candD0 = (AliAODRecoDecayHF2Prong*)cand->GetDaughter(1); //D0
  AliAODTrack *candD0_pr0 = (AliAODTrack*)candD0->GetDaughter(0); //Daughter 0
  AliAODTrack *candD0_pr1 = (AliAODTrack*)candD0->GetDaughter(1); //Daughter 1

  Double_t angleProngs = (candD0->Px() * cand_pr0->Px() + candD0->Py() * cand_pr0->Py() + candD0->Pz() * cand_pr0->Pz()) /(candD0->P() * cand_pr0->P());
  Double_t angleProngs_D0  = (candD0_pr1->Px() * candD0_pr0->Px() + candD0_pr1->Py() * candD0_pr0->Py() + candD0_pr1->Pz() * candD0_pr0->Pz()) /(candD0_pr1->P() * candD0_pr0->P());
    
  //topological variables
  //common (B+ -> D0 pi)
  fCandType.push_back(fCandTypeMap);
  fPt.push_back(((AliAODRecoDecayHF2Prong*)cand)->Pt());
  fY.push_back(((AliAODRecoDecayHF2Prong*)cand)->Y(521));
  fEta.push_back(((AliAODRecoDecayHF2Prong*)cand)->Eta());
  fPhi.push_back(((AliAODRecoDecayHF2Prong*)cand)->Phi());
  fDecayLength.push_back(((AliAODRecoDecayHF2Prong*)cand)->DecayLength());
  fDecayLengthXY.push_back(((AliAODRecoDecayHF2Prong*)cand)->DecayLengthXY());
  fNormDecayLengthXY.push_back(((AliAODRecoDecayHF2Prong*)cand)->NormalizedDecayLengthXY());
  fCosP.push_back(((AliAODRecoDecayHF2Prong*)cand)->CosPointingAngle());
  fCosPXY.push_back(((AliAODRecoDecayHF2Prong*)cand)->CosPointingAngleXY());
  fImpParXY.push_back(((AliAODRecoDecayHF2Prong*)cand)->ImpParXY());

  UInt_t prongs[2];
  prongs[0] = 211; prongs[1] = 421;
  fInvMass.push_back(((AliAODRecoDecayHF2Prong*)cand)->InvMass(2,prongs));
    
  fCosThetaStar.push_back(((AliAODRecoDecayHF2Prong*)cand)->CosThetaStar(0,521,211,421));
  fImpParProd.push_back(((AliAODRecoDecayHF2Prong*)cand)->Prodd0d0());
  fNormd0MeasMinusExp.push_back(ComputeMaxd0MeasMinusExp(cand,bfield));
  fDCA.push_back(((AliAODRecoDecayHF2Prong*)cand)->GetDCA());
  fAngleProngs.push_back(angleProngs);

  AliAODTrack* prongtracks[3];
  prongtracks[0] = (AliAODTrack*)cand->GetDaughter(0);
  fImpParProng[0].push_back(cand->Getd0Prong(0));

  //D0 -> K pi variables
  fPt_D0.push_back(candD0->Pt());
  fY_D0.push_back(candD0->Y(421));
  fEta_D0.push_back(candD0->Eta());
  fPhi_D0.push_back(candD0->Phi());
  fDecayLength_D0.push_back(candD0->DecayLength());
  fDecayLengthXY_D0.push_back(candD0->DecayLengthXY());
  fNormDecayLengthXY_D0.push_back(candD0->NormalizedDecayLengthXY());
  fCosP_D0.push_back(candD0->CosPointingAngle());
  fCosPXY_D0.push_back(candD0->CosPointingAngleXY());
  fImpParXY_D0.push_back(candD0->ImpParXY());
  fImpParProd_D0.push_back(candD0->Prodd0d0());
  fNormd0MeasMinusExp_D0.push_back(ComputeMaxd0MeasMinusExp(candD0,bfield));
  fDCA_D0.push_back(candD0->GetDCA());
  fAngleProngs_D0.push_back(angleProngs_D0);
    
  if(((AliAODRecoDecayHF2Prong*)cand)->Charge()==-1) {
    fInvMass_D0.push_back(candD0->InvMassD0());
    fCosThetaStar_D0.push_back(candD0->CosThetaStarD0());
      
    fImpParProng[1].push_back(candD0->Getd0Prong(0));
    fImpParProng[2].push_back(candD0->Getd0Prong(1));
    prongtracks[1] = (AliAODTrack*)candD0->GetDaughter(0);
    prongtracks[2] = (AliAODTrack*)candD0->GetDaughter(1);
  } else {
    fInvMass_D0.push_back(candD0->InvMassD0bar());
    fCosThetaStar_D0.push_back(candD0->CosThetaStarD0bar());

    fImpParProng[1].push_back(candD0->Getd0Prong(1));
    fImpParProng[2].push_back(candD0->Getd0Prong(0));
    prongtracks[1] = (AliAODTrack*)candD0->GetDaughter(1);
    prongtracks[2] = (AliAODTrack*)candD0->GetDaughter(0);
  }

  //single track variables
  bool setsingletrack = SetSingleTrackVars(prongtracks);
  if(!setsingletrack) return false;

  //pid variables
  if(fPidOpt==kNoPID) return true;

  bool setpid = SetPidVars(prongtracks,pidHF,true,true,false,true,true);
  if(!setpid) return false;

  return true;
}

//________________________________________________________________
void AliHFTreeHandlerBplustoD0pi::FillTree() {
  fTreeVar->Fill();

  //VERY IMPORTANT: CLEAR ALL VECTORS
  if(!fIsMCGenTree) {
    ResetDmesonCommonVarVectors();
    fCosThetaStar.clear();
    fImpParProd.clear();
    fNormd0MeasMinusExp.clear();
    fDCA.clear();
    fAngleProngs.clear();
    fInvMass_D0.clear();
    fPt_D0.clear();
    fY_D0.clear();
    fEta_D0.clear();
    fPhi_D0.clear();
    fDecayLength_D0.clear();
    fDecayLengthXY_D0.clear();
    fNormDecayLengthXY_D0.clear();
    fCosP_D0.clear();
    fCosPXY_D0.clear();
    fImpParXY_D0.clear();
    fCosThetaStar_D0.clear();
    fImpParProd_D0.clear();
    fNormd0MeasMinusExp_D0.clear();
    fDCA_D0.clear();
    fAngleProngs_D0.clear();
    for(unsigned int iProng=0; iProng<fNProngs; iProng++) fImpParProng[iProng].clear();
    ResetSingleTrackVarVectors();
    if(fPidOpt!=kNoPID) ResetPidVarVectors();
  }
  else {
    ResetMCGenVectors();
  }
  fCandTypeMap=0;
  fNCandidates=0;
}