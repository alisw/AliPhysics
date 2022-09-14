/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//*************************************************************************
// \class AliHFTreeHandlerApplyDstartoKpipi
// \brief helper class to handle a tree for Dstar cut optimisation and MVA analyses
// \authors:
// L. Vermunt, luuk.vermunt@cern.ch
/////////////////////////////////////////////////////////////

#include <TString.h>
#include <TDatabasePDG.h>
#include "AliHFTreeHandlerApplyDstartoKpipi.h"
#include "AliAODRecoCascadeHF.h"

/// \cond CLASSIMP
ClassImp(AliHFTreeHandlerApplyDstartoKpipi);
/// \endcond

//________________________________________________________________
AliHFTreeHandlerApplyDstartoKpipi::AliHFTreeHandlerApplyDstartoKpipi():
AliHFTreeHandlerApplyDstartoKpipi(kNsigmaPID)
{
  //
  // Default constructor
  //
}

//________________________________________________________________
AliHFTreeHandlerApplyDstartoKpipi::AliHFTreeHandlerApplyDstartoKpipi(int PIDopt):
  AliHFTreeHandlerApply(PIDopt),
  fCharge(-9999),
  fCosThetaStar(-9999.),
  fImpParProd(-9999.),
  fNormd0MeasMinusExp(-9999.),
  fAngleD0dkpPisoft(-9999.),
  fDeltaInvMassD0(-9999.)
{
  //
  // Standard constructor
  //

  //Only used in for-loops with a self-made AliAODtrack vector, so "Dstar pion + 2 D0-prongs" is fine
  fNProngs=3; // --> cannot be changed
  for(unsigned int iProng=0; iProng<fNProngs; iProng++){ 
    fIDProng[iProng] = -9999;
    fChargeProng[iProng] = -9999;
    fImpParProng[iProng] = -9999.;
  }
}

//________________________________________________________________
AliHFTreeHandlerApplyDstartoKpipi::~AliHFTreeHandlerApplyDstartoKpipi()
{
  //
  // Default Destructor
  //
}

//________________________________________________________________
TTree* AliHFTreeHandlerApplyDstartoKpipi::BuildTree(TString name, TString title)
{
  fIsMCGenTree=false;

  if(fTreeVar) {
    delete fTreeVar;
    fTreeVar=nullptr;
  }
  fTreeVar = new TTree(name.Data(),title.Data());

  //set common variables (please pay attention, filled with D0 and Dstar variables)
  AddCommonDmesonVarBranches();

  //set Dstar variables
  if(!fOnlyDedicatedBranches){
    fTreeVar->Branch("cos_t_star",&fCosThetaStar);
    fTreeVar->Branch("angle_D0dkpPisoft",&fAngleD0dkpPisoft);
    fTreeVar->Branch("delta_mass_D0",&fDeltaInvMassD0);
    for(unsigned int iProng=0; iProng<fNProngs; iProng++){
      fTreeVar->Branch(Form("imp_par_prong%d",iProng),&fImpParProng[iProng]);
    }
  }
  fTreeVar->Branch("imp_par_prod",&fImpParProd);
  fTreeVar->Branch("max_norm_d0d0exp",&fNormd0MeasMinusExp);
  fTreeVar->Branch("charge_cand",&fCharge);
  if(fFillOnlySignal) fTreeVar->Branch("mc_label", &fMCLabel);
  for(unsigned int iProng=0; iProng<fNProngs; iProng++){
    fTreeVar->Branch(Form("charge_prong%d",iProng),&fChargeProng[iProng]);
    fTreeVar->Branch(Form("id_prong%d",iProng),&fIDProng[iProng]);
  }

  //set single-track variables
  AddSingleTrackBranches();

  //set PID variables
  bool prongusepid[3] = {false, true, true};
  if(!fOnlyDedicatedBranches) prongusepid[0] = true;
  if(fPidOpt!=kNoPID) AddPidBranches(prongusepid,true,true,false,true,true);

  return fTreeVar;
}

//________________________________________________________________
bool AliHFTreeHandlerApplyDstartoKpipi::SetVariables(int runnumber, int eventID, int eventID_Ext, Long64_t eventID_Long, float ptgen, float mlprob, AliAODRecoDecayHF* cand, float bfield, int /*masshypo*/, AliPIDResponse *pidrespo, AliAODPidHF* pidhf)
{
  if(!cand) return false;
  if(fFillOnlySignal) { //if fill only signal and not signal candidate, do not store
    if(!(fCandType&kSignal)) return true;
  }
  fRunNumber=runnumber;
  fEvID=eventID;
  fEvIDExt=eventID_Ext;
  fEvIDLong=eventID_Long;
  fPtGen=ptgen;
  fMLProb=mlprob;
  
  fCandType &= ~kRefl; //protection --> Dstar ->Kpipi cannot be reflected

  //topological variables (Dstar and D0 variables combined)
  //common (Dstar)
  fPt=((AliAODRecoCascadeHF*)cand)->Pt();
  fY=((AliAODRecoCascadeHF*)cand)->YDstar();
  fEta=((AliAODRecoCascadeHF*)cand)->Eta();
  fPhi=((AliAODRecoCascadeHF*)cand)->Phi();
  fCharge=((AliAODRecoCascadeHF*)cand)->Charge();
  //common (D0)
  AliAODRecoDecayHF2Prong *d0 = ((AliAODRecoCascadeHF*)cand)->Get2Prong();
  fDecayLength=d0->DecayLength();
  fDecayLengthXY=d0->DecayLengthXY();
  fNormDecayLengthXY=d0->NormalizedDecayLengthXY()*(d0->P()/d0->Pt());
  fCosP=d0->CosPointingAngle();
  fCosPXY=d0->CosPointingAngleXY();
  fImpParXY=d0->ImpParXY();
  fDCA=d0->GetDCA();
  fNormd0MeasMinusExp=ComputeMaxd0MeasMinusExp(d0,bfield);
  fAngleD0dkpPisoft=((AliAODRecoCascadeHF*)cand)->AngleD0dkpPisoft();

  AliAODTrack* prongtracks[3];
  prongtracks[0] = (AliAODTrack*)((AliAODRecoCascadeHF*)cand)->GetBachelor();
  fImpParProng[0]=cand->Getd0Prong(0);
    
  //D* -> D0 pi variables
  fInvMass=((AliAODRecoCascadeHF*)cand)->DeltaInvMass();
  fImpParProd=d0->Prodd0d0();

  double massD0 = TDatabasePDG::Instance()->GetParticle(421)->Mass();
  if( (((AliAODRecoCascadeHF*)cand)->Charge()) > 0) {
    fCosThetaStar=d0->CosThetaStarD0();
    fDeltaInvMassD0=TMath::Abs(d0->InvMassD0() - massD0);
      
    fImpParProng[1]=d0->Getd0Prong(0);
    fImpParProng[2]=d0->Getd0Prong(1);
    prongtracks[1] = (AliAODTrack*)d0->GetDaughter(0);
    prongtracks[2] = (AliAODTrack*)d0->GetDaughter(1);
  } else {
    fCosThetaStar=d0->CosThetaStarD0bar();
    fDeltaInvMassD0=TMath::Abs(d0->InvMassD0bar() - massD0);

    fImpParProng[1]=d0->Getd0Prong(1);
    fImpParProng[2]=d0->Getd0Prong(0);
    prongtracks[1] = (AliAODTrack*)d0->GetDaughter(1);
    prongtracks[2] = (AliAODTrack*)d0->GetDaughter(0);
  }

  for(int ipr = 0; ipr < fNProngs; ipr++){
    fIDProng[ipr] = prongtracks[ipr]->GetID();
    fChargeProng[ipr] = prongtracks[ipr]->Charge();
  }
    
  //single track variables
  bool setsingletrack = SetSingleTrackVars(prongtracks);
  if(!setsingletrack) return false;

  //pid variables
  if(fPidOpt==kNoPID) return true;

  bool setpid = SetPidVars(prongtracks,pidrespo,true,true,false,true,true,pidhf);
  if(!setpid) return false;

  return true;
}
