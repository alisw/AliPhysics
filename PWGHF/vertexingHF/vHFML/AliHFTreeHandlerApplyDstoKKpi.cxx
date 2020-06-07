/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//*************************************************************************
// \class AliHFTreeHandlerApplyDstoKKpi
// \brief helper class to handle a tree for Ds cut optimisation and MVA analyses
// \authors:
// L. Vermunt, luuk.vermunt@cern.ch
/////////////////////////////////////////////////////////////

#include <TString.h>
#include <TDatabasePDG.h>
#include "AliHFTreeHandlerApplyDstoKKpi.h"
#include "AliAODRecoDecayHF3Prong.h"

/// \cond CLASSIMP
ClassImp(AliHFTreeHandlerApplyDstoKKpi);
/// \endcond

//________________________________________________________________
AliHFTreeHandlerApplyDstoKKpi::AliHFTreeHandlerApplyDstoKKpi():
AliHFTreeHandlerApplyDstoKKpi(kNsigmaPID)
{
  //
  // Default constructor
  //
}

//________________________________________________________________
AliHFTreeHandlerApplyDstoKKpi::AliHFTreeHandlerApplyDstoKKpi(int PIDopt):
AliHFTreeHandlerApply(PIDopt),
fImpParProng{},
fSigmaVertex(-9999.),
fMassKK(-9999.),
fCosPiDs(-9999.),
fCosPiKPhi(-9999.),
fNormd0MeasMinusExp(-9999.),
fMassKKOpt(kMassKK)
{  
  //
  // Standard constructor
  //
  
  fNProngs=3; // --> cannot be changed
  for(unsigned int iProng=0; iProng<fNProngs; iProng++)
    fImpParProng[iProng] = -9999.;
}

//________________________________________________________________
AliHFTreeHandlerApplyDstoKKpi::~AliHFTreeHandlerApplyDstoKKpi()
{
  //
  // Default Destructor
  //
}

//________________________________________________________________
TTree* AliHFTreeHandlerApplyDstoKKpi::BuildTree(TString name, TString title)
{
  fIsMCGenTree=false;
  
  if(fTreeVar) {
    delete fTreeVar;
    fTreeVar=nullptr;
  }
  fTreeVar = new TTree(name.Data(),title.Data());
  
  //set common variables
  AddCommonDmesonVarBranches();
  
  //set Ds variables
  TString massKKname="";
  if(fMassKKOpt==kMassKK) massKKname = "mass_KK";
  else if(fMassKKOpt==kDeltaMassKKPhi) massKKname = "delta_mass_KK";
  fTreeVar->Branch("sig_vert",&fSigmaVertex);
  fTreeVar->Branch(massKKname.Data(),&fMassKK);
  fTreeVar->Branch("cos_PiDs",&fCosPiDs);
  fTreeVar->Branch("cos_PiKPhi_3",&fCosPiKPhi);
  fTreeVar->Branch("max_norm_d0d0exp",&fNormd0MeasMinusExp);
  for(unsigned int iProng=0; iProng<fNProngs; iProng++){
    fTreeVar->Branch(Form("imp_par_prong%d",iProng),&fImpParProng[iProng]);
  }
  
  //set single-track variables
  AddSingleTrackBranches();
  
  //sed pid variables
  if(fPidOpt!=kNoPID) AddPidBranches(true,true,false,true,true);
  
  return fTreeVar;
}

//________________________________________________________________
bool AliHFTreeHandlerApplyDstoKKpi::SetVariables(int runnumber, int eventID, int eventID_Ext, Long64_t eventID_Long, float ptgen, float mlprob, AliAODRecoDecayHF* cand, float bfield, int masshypo, AliPIDResponse *pidrespo)
{
  if(!cand) return false;
  if(fFillOnlySignal) { //if fill only signal and not signal candidate, do not store
    if(!(fCandType&kSignal || fCandType&kRefl)) return true;
  }
  fRunNumber=runnumber;
  fEvID=eventID;
  fEvIDExt=eventID_Ext;
  fEvIDLong=eventID_Long;
  fPtGen=ptgen;
  fMLProb=mlprob;
  
  //topological variables
  //common
  fPt=cand->Pt();
  fY=cand->Y(431);
  fEta=cand->Eta();
  fPhi=cand->Phi();
  fDecayLength=cand->DecayLength();
  fDecayLengthXY=cand->DecayLengthXY();
  fNormDecayLengthXY=cand->NormalizedDecayLengthXY();
  fCosP=cand->CosPointingAngle();
  fCosPXY=cand->CosPointingAngleXY();
  fImpParXY=cand->ImpParXY();
  fDCA=cand->GetDCA();
  fNormd0MeasMinusExp=ComputeMaxd0MeasMinusExp(cand,bfield);
  
  //Ds+ -> KKpi variables
  fSigmaVertex=((AliAODRecoDecayHF3Prong*)cand)->GetSigmaVert();
  float massPhi = 0;
  float cospikphi=-2;
  if(fMassKKOpt==kDeltaMassKKPhi) massPhi = TDatabasePDG::Instance()->GetParticle(333)->Mass();
  if(masshypo==0){ //phiKKpi
    fInvMass=((AliAODRecoDecayHF3Prong*)cand)->InvMassDsKKpi();
    fMassKK=TMath::Abs(((AliAODRecoDecayHF3Prong*)cand)->InvMass2Prongs(0,1,321,321)-massPhi);
    fCosPiDs=((AliAODRecoDecayHF3Prong*)cand)->CosPiDsLabFrameKKpi();
    cospikphi = ((AliAODRecoDecayHF3Prong*)cand)->CosPiKPhiRFrameKKpi();
  }
  else if(masshypo==1){ //phipiKK
    fInvMass=((AliAODRecoDecayHF3Prong*)cand)->InvMassDspiKK();
    fMassKK=TMath::Abs(((AliAODRecoDecayHF3Prong*)cand)->InvMass2Prongs(1,2,321,321)-massPhi);
    fCosPiDs=((AliAODRecoDecayHF3Prong*)cand)->CosPiDsLabFramepiKK();
    cospikphi = ((AliAODRecoDecayHF3Prong*)cand)->CosPiKPhiRFramepiKK();
  }
  fCosPiKPhi=cospikphi*cospikphi*cospikphi;
  for(unsigned int iProng=0; iProng<fNProngs; iProng++) {
    fImpParProng[iProng]=cand->Getd0Prong(iProng);
  }
  
  //single-track variables
  AliAODTrack* prongtracks[3];
  for(unsigned int iProng=0; iProng<fNProngs; iProng++) prongtracks[iProng] = (AliAODTrack*)cand->GetDaughter(iProng);
  bool setsingletrack = SetSingleTrackVars(prongtracks);  
  if(!setsingletrack) return false;
  
  //pid variables
  if(fPidOpt==kNoPID) return true;
  
  bool setpid = SetPidVars(prongtracks,pidrespo,true,true,false,true,true);
  if(!setpid) return false;
  
  return true;
}
