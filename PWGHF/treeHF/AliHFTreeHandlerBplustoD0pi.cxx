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
// J. Norman, jaime.norman@cern.ch
// G. Luparello, grazia.luparello@cern.ch
/////////////////////////////////////////////////////////////

#include <TString.h>
#include "AliHFTreeHandlerBplustoD0pi.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliESDVertex.h"
#include "AliESDtrackCuts.h"
#include "AliRDHFCutsD0toKpi.h"

/// \cond CLASSIMP
ClassImp(AliHFTreeHandlerBplustoD0pi);
/// \endcond

//________________________________________________________________
AliHFTreeHandlerBplustoD0pi::AliHFTreeHandlerBplustoD0pi():
  AliHFTreeHandler(),
  fCosThetaStar(-9999.),
  fImpParProd(-9999.),
  fNormd0MeasMinusExp(-9999.),
  fAngleProngs(-9999.),
  fInvMass_D0(-9999.),
  fPt_D0(-9999.),
  fY_D0(-9999.),
  fEta_D0(-9999.),
  fPhi_D0(-9999.),
  fDecayLength_D0(-9999.),
  fDecayLengthXY_D0(-9999.),
  fNormDecayLengthXY_D0(-9999.),
  fCosP_D0(-9999.),
  fCosPXY_D0(-9999.),
  fImpParXY_D0(-9999.),
  fCosThetaStar_D0(-9999.),
  fImpParProd_D0(-9999.),
  fNormd0MeasMinusExp_D0(-9999.),
  fDCA_D0(-9999.),
  fAngleProngs_D0(-9999.)
{
  //
  // Default constructor
  //

  //Only used in for-loops with a self-made AliAODtrack vector, so "Bplus pion + 2 D0-prongs" is fine
  fNProngs=3; // --> cannot be changed
  for(unsigned int iProng=0; iProng<fNProngs; iProng++) 
    fImpParProng[iProng] = -9999.;
}

//________________________________________________________________
AliHFTreeHandlerBplustoD0pi::AliHFTreeHandlerBplustoD0pi(int PIDopt):
  AliHFTreeHandler(PIDopt),
  fCosThetaStar(-9999.),
  fImpParProd(-9999.),
  fNormd0MeasMinusExp(-9999.),
  fAngleProngs(-9999.),
  fInvMass_D0(-9999.),
  fPt_D0(-9999.),
  fY_D0(-9999.),
  fEta_D0(-9999.),
  fPhi_D0(-9999.),
  fDecayLength_D0(-9999.),
  fDecayLengthXY_D0(-9999.),
  fNormDecayLengthXY_D0(-9999.),
  fCosP_D0(-9999.),
  fCosPXY_D0(-9999.),
  fImpParXY_D0(-9999.),
  fCosThetaStar_D0(-9999.),
  fImpParProd_D0(-9999.),
  fNormd0MeasMinusExp_D0(-9999.),
  fDCA_D0(-9999.),
  fAngleProngs_D0(-9999.)
{
  //
  // Standard constructor
  //
    
  //Only used in for-loops with a self-made AliAODtrack vector, so "Bplus pion + 2 D0-prongs" is fine
  fNProngs=3; // --> cannot be changed
  for(unsigned int iProng=0; iProng<fNProngs; iProng++) 
    fImpParProng[iProng] = -9999.;
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
    fTreeVar=nullptr;
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
  if (fFillJets) AddJetBranches();

  //set PID variables
  if(fPidOpt!=kNoPID) AddPidBranches(true,true,false,true,true);

  return fTreeVar;
}

//________________________________________________________________
bool AliHFTreeHandlerBplustoD0pi::SetVariables(int runnumber, int eventID, int eventID_Ext, Long64_t eventID_Long, float ptgen, AliAODRecoDecayHF* cand, float bfield, int /*masshypo*/, AliPIDResponse* pidrespo)
{
  fRunNumber=runnumber;
  fEvID=eventID;
  fEvIDExt=eventID_Ext;
  fEvIDLong=eventID_Long;
  fIsMCGenTree=false;

  if(!cand) return false;
  if(fFillOnlySignal) { //if fill only signal and not signal candidate, do not store
    if(!(fCandType&kSignal)) return true;
  }

  fPtGen=ptgen;
  
  fCandType &= ~kRefl; //protection --> Bplus -> D0pi cannot be reflected

  AliAODTrack* cand_pr0 = (AliAODTrack*)cand->GetDaughter(0); //Bplus pion
  AliAODRecoDecayHF2Prong* candD0 = (AliAODRecoDecayHF2Prong*)cand->GetDaughter(1); //D0
  AliAODTrack *candD0_pr0 = (AliAODTrack*)candD0->GetDaughter(0); //Daughter 0
  AliAODTrack *candD0_pr1 = (AliAODTrack*)candD0->GetDaughter(1); //Daughter 1

  Double_t angleProngs = (candD0->Px() * cand_pr0->Px() + candD0->Py() * cand_pr0->Py() + candD0->Pz() * cand_pr0->Pz()) /(candD0->P() * cand_pr0->P());
  Double_t angleProngs_D0  = (candD0_pr1->Px() * candD0_pr0->Px() + candD0_pr1->Py() * candD0_pr0->Py() + candD0_pr1->Pz() * candD0_pr0->Pz()) /(candD0_pr1->P() * candD0_pr0->P());
    
  //topological variables
  //common (B+ -> D0 pi)
  fPt=((AliAODRecoDecayHF2Prong*)cand)->Pt();
  fY=((AliAODRecoDecayHF2Prong*)cand)->Y(521);
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
  prongs[0] = 211; prongs[1] = 421;
  fInvMass=((AliAODRecoDecayHF2Prong*)cand)->InvMass(2,prongs);
    
  fCosThetaStar=((AliAODRecoDecayHF2Prong*)cand)->CosThetaStar(0,521,211,421);
  fImpParProd=((AliAODRecoDecayHF2Prong*)cand)->Prodd0d0();
  fNormd0MeasMinusExp=ComputeMaxd0MeasMinusExp(cand,bfield);
  fAngleProngs=angleProngs;

  AliAODTrack* prongtracks[3];
  prongtracks[0] = (AliAODTrack*)cand->GetDaughter(0);
  fImpParProng[0]=cand->Getd0Prong(0);

  //D0 -> K pi variables
  fPt_D0=candD0->Pt();
  fY_D0=candD0->Y(421);
  fEta_D0=candD0->Eta();
  fPhi_D0=candD0->Phi();
  fDecayLength_D0=candD0->DecayLength();
  fDecayLengthXY_D0=candD0->DecayLengthXY();
  fNormDecayLengthXY_D0=candD0->NormalizedDecayLengthXY();
  fCosP_D0=candD0->CosPointingAngle();
  fCosPXY_D0=candD0->CosPointingAngleXY();
  fImpParXY_D0=candD0->ImpParXY();
  fImpParProd_D0=candD0->Prodd0d0();
  fNormd0MeasMinusExp_D0=ComputeMaxd0MeasMinusExp(candD0,bfield);
  fDCA_D0=candD0->GetDCA();
  fAngleProngs_D0=angleProngs_D0;
    
  if(((AliAODRecoDecayHF2Prong*)cand)->Charge()==-1) {
    fInvMass_D0=candD0->InvMassD0();
    fCosThetaStar_D0=candD0->CosThetaStarD0();
      
    fImpParProng[1]=candD0->Getd0Prong(0);
    fImpParProng[2]=candD0->Getd0Prong(1);
    prongtracks[1] = (AliAODTrack*)candD0->GetDaughter(0);
    prongtracks[2] = (AliAODTrack*)candD0->GetDaughter(1);
  } else {
    fInvMass_D0=candD0->InvMassD0bar();
    fCosThetaStar_D0=candD0->CosThetaStarD0bar();

    fImpParProng[1]=candD0->Getd0Prong(1);
    fImpParProng[2]=candD0->Getd0Prong(0);
    prongtracks[1] = (AliAODTrack*)candD0->GetDaughter(1);
    prongtracks[2] = (AliAODTrack*)candD0->GetDaughter(0);
  }

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
Int_t AliHFTreeHandlerBplustoD0pi::IsBplusPionSelected(TObject* obj, AliRDHFCutsD0toKpi* cutsD0, AliAODPidHF* fPidHFD0, AliAODEvent* aod, AliAODVertex *vtx) {
  
  AliAODTrack* candidatePion = (AliAODTrack*)obj;
  if (!candidatePion){ AliWarning("No pion object. Track rejected."); return 0; }
  
  AliESDtrackCuts* fTrackCuts = cutsD0->GetTrackCuts();
  if (!fTrackCuts){ AliWarning("No fTrackCuts object. Track rejected.");  return 0; }
  
  Double_t pos[3],cov[6];
  vtx->GetXYZ(pos);
  vtx->GetCovarianceMatrix(cov);
  const AliESDVertex vESD(pos,cov,100.,100);
  
  if(!cutsD0->IsDaughterSelected(candidatePion,&vESD,fTrackCuts,aod)) return 0;
  
  if(!cutsD0->GetIsUsePID()) return 1;
  else {
    
    if(!fPidHFD0){ AliWarning("AliAODPidHF not created. Track accepted"); return 1; }
    
    Int_t isPion=fPidHFD0->MakeRawPid(candidatePion,AliPID::kPion);
    if(isPion) return 1;
    else       return 0;
  }
  
  return 1;
}
