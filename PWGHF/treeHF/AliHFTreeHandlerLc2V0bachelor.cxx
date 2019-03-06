/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//*************************************************************************
// \class AliHFTreeHandlerLc2V0bachelor
// \brief helper class to handle a tree for Lc->K0sp cut optimisation and MVA analyses
// \authors:
// C. Zampolli (Chiara.Zampolli@cern.ch)
/////////////////////////////////////////////////////////////

#include <TString.h>
#include <TDatabasePDG.h>
#include "AliHFTreeHandlerLc2V0bachelor.h"
#include "AliAODRecoCascadeHF.h"

/// \cond CLASSIMP
ClassImp(AliHFTreeHandlerLc2V0bachelor);
/// \endcond

//________________________________________________________________
AliHFTreeHandlerLc2V0bachelor::AliHFTreeHandlerLc2V0bachelor():
  AliHFTreeHandler(),
  fImpParProng(),
  fImpParK0s(),
  fDecayLengthK0s(),
  fInvMassK0s(),
  fDCAK0s(),
  fPtK0s(),
  fEtaK0s(),
  fPhiK0s(),
  fcTauK0s(),
  fV0PointingAngle(),
  fCosThetaStar(),
  fsignd0(),
  fArmqTOverAlpha(),
  fCalcSecoVtx(0)
{
  //
  // Default constructor
  //

  fNProngs = 3; // --> cannot be changed (prong 0 is the bachelor, 1,2 are prongs of the K0s)
}

//________________________________________________________________
AliHFTreeHandlerLc2V0bachelor::AliHFTreeHandlerLc2V0bachelor(int PIDopt):
  AliHFTreeHandler(PIDopt),
  fImpParProng(),
  fImpParK0s(),
  fDecayLengthK0s(),
  fInvMassK0s(),
  fDCAK0s(),
  fPtK0s(),
  fEtaK0s(),
  fPhiK0s(),
  fcTauK0s(),
  fV0PointingAngle(),
  fCosThetaStar(),
  fsignd0(),
  fArmqTOverAlpha(),
  fCalcSecoVtx(0)
{
  //
  // Standard constructor
  //

  fNProngs = 3; // --> cannot be changed (prong 0 is the bachelor, 1,2 are prongs of the K0s)
}

//________________________________________________________________
AliHFTreeHandlerLc2V0bachelor::~AliHFTreeHandlerLc2V0bachelor()
{
  //
  // Default Destructor
  //
}

//________________________________________________________________
TTree* AliHFTreeHandlerLc2V0bachelor::BuildTree(TString name, TString title) 
{
  fIsMCGenTree=false;

  if (fTreeVar) {
    delete fTreeVar;
    fTreeVar = 0x0;
  }
  fTreeVar = new TTree(name.Data(), title.Data());

  //set common variables
  AddCommonDmesonVarBranches();

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
  fTreeVar->Branch("eta_K0s", &fEtaK0s);
  fTreeVar->Branch("phi_K0s", &fPhiK0s);
  for(unsigned int iProng=0; iProng<fNProngs; iProng++){
    fTreeVar->Branch(Form("imp_par_prong%d", iProng), &fImpParProng[iProng]);
  }

  //set single-track variables
  AddSingleTrackBranches();

  //set PID variables
  if(fPidOpt != kNoPID) AddPidBranches(true, false, true, true, true);

  return fTreeVar;
}

//________________________________________________________________
bool AliHFTreeHandlerLc2V0bachelor::SetVariables(int runnumber, unsigned int eventID, AliAODRecoDecayHF* cand, float bfield, int masshypo, AliPIDResponse* pidrespo) 
{
  if(!cand) return false;
  if(fFillOnlySignal) { //if fill only signal and not signal candidate, do not store
    if(!(fCandTypeMap&kSignal)) return true;
  }
  fNCandidates++;
  fRunNumber.push_back(runnumber);
  fEvID.push_back(eventID);
    
  //topological variables
  //common
  fCandType.push_back(fCandTypeMap);
  fCandTypeMap=0; //reset candtype
  fPt.push_back(cand->Pt());
  fY.push_back(cand->Y(411));
  fEta.push_back(cand->Eta());
  fPhi.push_back(cand->Phi());
  if(fCalcSecoVtx){
    fDecayLength.push_back(cand->DecayLength());
    fDecayLengthXY.push_back(cand->DecayLengthXY());
    fNormDecayLengthXY.push_back(cand->NormalizedDecayLengthXY());
    fCosP.push_back(cand->CosPointingAngle());
    fCosPXY.push_back(cand->CosPointingAngleXY());
    fImpParXY.push_back(cand->ImpParXY());
    fDCA.push_back(cand->GetDCA());
  }
  
  //Lc -> K0sp variables  
  fInvMass.push_back(((AliAODRecoCascadeHF*)cand)->InvMassLctoK0sP());
  
  AliAODv0 * v0part = ((AliAODRecoCascadeHF*)cand)->Getv0();
  for(unsigned int iProng = 0; iProng < fNProngs; iProng++) {
    if(iProng==0) fImpParProng[iProng].push_back(cand->Getd0Prong(iProng));
    else          fImpParProng[iProng].push_back(v0part->Getd0Prong(iProng-1));
  }
  fImpParK0s.push_back(cand->Getd0Prong(1));
  fDecayLengthK0s.push_back(((AliAODRecoCascadeHF*)cand)->DecayLengthV0());
  fInvMassK0s.push_back(v0part->MassK0Short());
  fDCAK0s.push_back(v0part->GetDCA());
  fPtK0s.push_back(v0part->Pt());
  fEtaK0s.push_back(v0part->Eta());
  fPhiK0s.push_back(v0part->Phi());

  fcTauK0s.push_back(((AliAODRecoCascadeHF*)cand)->DecayLengthV0()*0.497/(v0part->P()));
  fV0PointingAngle.push_back(((AliAODRecoCascadeHF*)cand)->CosV0PointingAngle());

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
  fCosThetaStar.push_back(cts);

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
  fsignd0.push_back(signd0);
  
  //Armenteros qT/|alpha|
  fArmqTOverAlpha.push_back( v0part->PtArmV0()/TMath::Abs(v0part->AlphaV0()) );
  
  //single track variables
  AliAODTrack* prongtracks[3];
  for(unsigned int iProng = 0; iProng < fNProngs; iProng++){
    if(iProng==0) prongtracks[iProng] = (AliAODTrack*)cand->GetDaughter(iProng);
    else          prongtracks[iProng] = (AliAODTrack*)v0part->GetDaughter(iProng-1);
  }
  bool setsingletrack = SetSingleTrackVars(prongtracks);  
  if(!setsingletrack) return false;

  //pid variables
  if(fPidOpt == kNoPID) return true;

  bool setpid = SetPidVars(prongtracks, pidrespo, true, false, true, true, true);
  if(!setpid) return false;

  return true;
}

//________________________________________________________________
void AliHFTreeHandlerLc2V0bachelor::FillTree() {
  fTreeVar->Fill();
  
  //VERY IMPORTANT: CLEAR ALL VECTORS
  if(!fIsMCGenTree) {
    ResetDmesonCommonVarVectors();
    fImpParK0s.clear();
    fDecayLengthK0s.clear();
    fInvMassK0s.clear();
    fDCAK0s.clear();
    fPtK0s.clear();
    fEtaK0s.clear();
    fPhiK0s.clear();
    fcTauK0s.clear();
    fV0PointingAngle.clear();
    fCosThetaStar.clear();
    fsignd0.clear();
    fArmqTOverAlpha.clear();
    for(unsigned int iProng = 0; iProng<fNProngs; iProng++) fImpParProng[iProng].clear();
    ResetSingleTrackVarVectors();
    if(fPidOpt != kNoPID) ResetPidVarVectors();
  }
  else {
    ResetMCGenVectors();
  }
  fCandTypeMap = 0;
  fNCandidates = 0;
}
