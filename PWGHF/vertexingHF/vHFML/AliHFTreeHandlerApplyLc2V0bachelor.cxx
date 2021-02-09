/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//*************************************************************************
// \class AliHFTreeHandlerApplyLc2V0bachelor
// \brief helper class to handle a tree for Lc->K0sp cut optimisation and MVA analyses
// \authors:
// L. Vermunt, luuk.vermunt@cern.ch
/////////////////////////////////////////////////////////////

#include <TString.h>
#include <TDatabasePDG.h>
#include "AliHFTreeHandlerApplyLc2V0bachelor.h"
#include "AliAODRecoCascadeHF.h"
#include "AliESDtrack.h"

/// \cond CLASSIMP
ClassImp(AliHFTreeHandlerApplyLc2V0bachelor);
/// \endcond

//________________________________________________________________
AliHFTreeHandlerApplyLc2V0bachelor::AliHFTreeHandlerApplyLc2V0bachelor():
AliHFTreeHandlerApplyLc2V0bachelor(kNsigmaPID)
{
  //
  // Default constructor
  //
}

//________________________________________________________________
AliHFTreeHandlerApplyLc2V0bachelor::AliHFTreeHandlerApplyLc2V0bachelor(int PIDopt):
AliHFTreeHandlerApply(PIDopt),
fImpParProng{},
fITSRefitProng{},
fImpParK0s(-9999.),
fDecayLengthK0s(-9999.),
fInvMassK0s(-9999.),
fDCAK0s(-9999.),
fPtK0s(-9999.),
fEtaK0s(-9999.),
fPhiK0s(-9999.),
fcTauK0s(-9999.),
fV0PointingAngle(-9999.),
fCosThetaStar(-9999.),
fsignd0(-9999.),
fArmqTOverAlpha(-9999.),
fCalcSecoVtx(0),
fReducePbPbBranches(false)
{
  //
  // Standard constructor
  //
  
  fNProngs = 3; // --> cannot be changed (prong 0 is the bachelor, 1,2 are prongs of the K0s)
  for(unsigned int iProng=0; iProng<fNProngs; iProng++){
    fImpParProng[iProng] = -9999.;
    fITSRefitProng[iProng] = -9999;
  }
}

//________________________________________________________________
AliHFTreeHandlerApplyLc2V0bachelor::~AliHFTreeHandlerApplyLc2V0bachelor()
{
  //
  // Default Destructor
  //
}

//________________________________________________________________
TTree* AliHFTreeHandlerApplyLc2V0bachelor::BuildTree(TString name, TString title)
{
  fIsMCGenTree=false;
  
  if (fTreeVar) {
    delete fTreeVar;
    fTreeVar = nullptr;
  }
  fTreeVar = new TTree(name.Data(), title.Data());
  
  //set common variables
  AddCommonDmesonVarBranches(fCalcSecoVtx);
  
  //set Lc variables
  fTreeVar->Branch("cos_t_star", &fCosThetaStar);
  fTreeVar->Branch("signd0", &fsignd0);
  fTreeVar->Branch("inv_mass_K0s", &fInvMassK0s);
  fTreeVar->Branch("dca_K0s", &fDCAK0s);
  fTreeVar->Branch("imp_par_K0s", &fImpParK0s);
  fTreeVar->Branch("d_len_K0s", &fDecayLengthK0s);
  fTreeVar->Branch("armenteros_K0s", &fArmqTOverAlpha);
  fTreeVar->Branch("ctau_K0s", &fcTauK0s);
  fTreeVar->Branch("cos_p_K0s", &fV0PointingAngle);
  fTreeVar->Branch("pt_K0s", &fPtK0s);
  if(!fReducePbPbBranches){
    fTreeVar->Branch("eta_K0s", &fEtaK0s);
    fTreeVar->Branch("phi_K0s", &fPhiK0s);
  }
  for(unsigned int iProng=0; iProng<fNProngs; iProng++){
    fTreeVar->Branch(Form("imp_par_prong%d", iProng), &fImpParProng[iProng]);
  }
  fTreeVar->Branch("its_refit_prong1",&fITSRefitProng[1]);
  fTreeVar->Branch("its_refit_prong2",&fITSRefitProng[2]);
  
  //set single-track variables
  AddSingleTrackBranches();
  
  //set PID variables
  if(fPidOpt != kNoPID){
    if(!fReducePbPbBranches){
      bool prongusepid[3] = {true, true, true};
      AddPidBranches(prongusepid, true, true, true, true, true);
    } else {
      bool prongusepid[3] = {true, false, false};
      AddPidBranches(prongusepid, false, false, true, true, true);
    }
  }
  
  return fTreeVar;
}

//________________________________________________________________
bool AliHFTreeHandlerApplyLc2V0bachelor::SetVariables(int runnumber, int eventID, int eventID_Ext, Long64_t eventID_Long, float ptgen, float mlprob, AliAODRecoDecayHF* cand, float bfield, int masshypo, AliPIDResponse* pidrespo, AliAODPidHF* pidhf)
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

  //topological variables
  //common
  fPt=cand->Pt();
  fY=cand->Y(4122);
  fEta=cand->Eta();
  fPhi=cand->Phi();
  if(fCalcSecoVtx){
    fDecayLength=cand->DecayLength();
    fDecayLengthXY=cand->DecayLengthXY();
    fNormDecayLengthXY=cand->NormalizedDecayLengthXY();
    fCosP=cand->CosPointingAngle();
    fCosPXY=cand->CosPointingAngleXY();
    fImpParXY=cand->ImpParXY();
    fDCA=cand->GetDCA();
  }
  
  //Lc -> K0sp variables
  fInvMass=((AliAODRecoCascadeHF*)cand)->InvMassLctoK0sP();
  
  AliAODv0 * v0part = ((AliAODRecoCascadeHF*)cand)->Getv0();
  for(unsigned int iProng = 0; iProng < fNProngs; iProng++) {
    if(iProng==0) fImpParProng[iProng]=cand->Getd0Prong(iProng);
    else          fImpParProng[iProng]=v0part->Getd0Prong(iProng-1);
  }
  fImpParK0s=cand->Getd0Prong(1);
  fDecayLengthK0s=((AliAODRecoCascadeHF*)cand)->DecayLengthV0();
  fInvMassK0s=v0part->MassK0Short();
  fDCAK0s=v0part->GetDCA();
  fPtK0s=v0part->Pt();
  fEtaK0s=v0part->Eta();
  fPhiK0s=v0part->Phi();
  
  fcTauK0s=((AliAODRecoCascadeHF*)cand)->DecayLengthV0()*0.497/(v0part->P());
  fV0PointingAngle=((AliAODRecoCascadeHF*)cand)->CosV0PointingAngle();
  
  // Cosine of proton emission angle (theta*) in the rest frame of the mother particle
  // (from AliRDHFCutsLctoV0)
  TLorentzVector vpr, vk0s,vlc;
  Double_t massK0SPDG = TDatabasePDG::Instance()->GetParticle(310)->Mass();    // mass K0S PDG
  Double_t massPrPDG = TDatabasePDG::Instance()->GetParticle(2212)->Mass();    // mass Proton PDG
  vpr.SetXYZM(cand->PxProng(0), cand->PyProng(0), cand->PzProng(0), massPrPDG);
  vk0s.SetXYZM(cand->PxProng(1), cand->PyProng(1), cand->PzProng(1), massK0SPDG);
  vlc = vpr + vk0s;
  TVector3 vboost = vlc.BoostVector();
  vpr.Boost(-vboost);
  Double_t cts = TMath::Cos(vpr.Angle(vlc.Vect()));
  fCosThetaStar=cts;
  
  // Sign of d0 proton (different from regular d0)
  // (from AliRDHFCutsLctoV0)
  AliAODTrack *bachelor = (AliAODTrack*)((AliAODRecoCascadeHF*)cand)->GetBachelor();
  AliAODVertex *primvert = dynamic_cast<AliAODVertex*>(cand->GetPrimaryVtx());
  Double_t d0z0bach[2], covd0z0bach[3];
  bachelor->PropagateToDCA(primvert, bfield, kVeryBig, d0z0bach, covd0z0bach); // HOW DO WE SET THE B FIELD?; kVeryBig should come from AliExternalTrackParam
  Double_t tx[3];
  bachelor->GetXYZ(tx);
  tx[0] -= primvert->GetX();
  tx[1] -= primvert->GetY();
  tx[2] -= primvert->GetZ();
  Double_t innerpro = tx[0]*cand->Px()+tx[1]*cand->Py();
  Double_t signd0 = 1.;
  if(innerpro<0.) signd0 = -1.;
  signd0 = signd0*TMath::Abs(d0z0bach[0]);
  fsignd0=signd0;
  
  //Armenteros qT/|alpha|
  fArmqTOverAlpha= v0part->PtArmV0()/TMath::Abs(v0part->AlphaV0());
  
  //single track variables
  AliAODTrack* prongtracks[3];
  for(unsigned int iProng = 0; iProng < fNProngs; iProng++){
    if(iProng==0) prongtracks[iProng] = (AliAODTrack*)cand->GetDaughter(iProng);
    else          prongtracks[iProng] = (AliAODTrack*)v0part->GetDaughter(iProng-1);
    fITSRefitProng[iProng] = (prongtracks[iProng]->GetStatus()&AliESDtrack::kITSrefit);
  }
  bool setsingletrack = SetSingleTrackVars(prongtracks);
  if(!setsingletrack) return false;
  
  //pid variables
  if(fPidOpt == kNoPID) return true;
  
  bool setpid;
  if(!fReducePbPbBranches) setpid = SetPidVars(prongtracks, pidrespo, true, true, true, true, true, pidhf);
  else                     setpid = SetPidVars(prongtracks, pidrespo, false, false, true, true, true, pidhf);

  if(!setpid) return false;
  
  return true;
}
