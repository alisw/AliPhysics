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
  fInvMassK0s(),
  fcTauK0s(),
  fV0PointingAngle(),
  fCosThetaStar(),
  fsignd0()
{
  //
  // Default constructor
  //

  fNProngs = 2; // --> cannot be changed (prong 0 is the bachelor, prong 1 is the K0s)
}

//________________________________________________________________
AliHFTreeHandlerLc2V0bachelor::AliHFTreeHandlerLc2V0bachelor(int PIDopt):
  AliHFTreeHandler(PIDopt),
  fImpParProng(),
  fInvMassK0s(),
  fcTauK0s(),
  fV0PointingAngle(),
  fCosThetaStar(),
  fsignd0()
{
  //
  // Standard constructor
  //

  fNProngs = 2; // --> cannot be changed
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
  fTreeVar->Branch("invmassK0s", &fInvMassK0s);
  for(unsigned int iProng=0; iProng<fNProngs; iProng++){
    fTreeVar->Branch(Form("imp_par_prong%d", iProng), &fImpParProng[iProng]);
  }
  fTreeVar->Branch("cTauK0s", &fcTauK0s);
  fTreeVar->Branch("V0PointingAngle", &fV0PointingAngle);
  fTreeVar->Branch("cosThetaStar", &fCosThetaStar);
  fTreeVar->Branch("signd0", &fsignd0);

  //set single-track variables
  AddSingleTrackBranches();

  //set PID variables
  if(fPidOpt != kNoPID) AddPidBranches(false, false, true, true, true);

  return fTreeVar;
}

//________________________________________________________________
bool AliHFTreeHandlerLc2V0bachelor::SetVariables(AliAODRecoCascadeHF* cand, float bfield, int masshypo, AliPIDResponse* pidrespo) 
{
  if(!cand) return false;
  if(fFillOnlySignal) { //if fill only signal and not signal candidate, do not store
    if(!(fCandTypeMap&kSignal)) return true;
  }
  fNCandidates++;

  //topological variables
  //common
  fCandType.push_back(fCandTypeMap);
  fPt.push_back(cand->Pt());
  fY.push_back(cand->Y(411));
  fEta.push_back(cand->Eta());
  fPhi.push_back(cand->Phi());
  fDecayLength.push_back(cand->DecayLength());
  fDecayLengthXY.push_back(cand->DecayLengthXY());
  fNormDecayLengthXY.push_back(cand->NormalizedDecayLengthXY());
  fCosP.push_back(cand->CosPointingAngle());
  fCosPXY.push_back(cand->CosPointingAngleXY());
  fImpParXY.push_back(cand->ImpParXY());

  //Lc -> K0sp variables  
  fInvMass.push_back(((AliAODRecoCascadeHF*)cand)->InvMassLctoK0sP());

  for(unsigned int iProng = 0; iProng < fNProngs; iProng++) {
    fImpParProng[iProng].push_back(cand->Getd0Prong(iProng));
  }
  
  AliAODv0 * v0part = ((AliAODRecoCascadeHF*)cand)->Getv0();
  fInvMassK0s.push_back(v0part->MassK0Short());
  
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

  AliAODTrack *bachelor = (AliAODTrack*)cand->GetBachelor();
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
			     
  //single track variables
  AliAODTrack* prongtracks[2];
  for(unsigned int iProng = 0; iProng < fNProngs; iProng++) prongtracks[iProng] = (AliAODTrack*)cand->GetDaughter(iProng);
  bool setsingletrack = SetSingleTrackVars(prongtracks);  
  if(!setsingletrack) return false;

  //pid variables
  if(fPidOpt == kNoPID) return true;

  bool setpid = SetPidVars(prongtracks, pidrespo, true, true, true, true, true);
  if(!setpid) return false;

  return true;
}

//________________________________________________________________
void AliHFTreeHandlerLc2V0bachelor::FillTree() {
  fTreeVar->Fill();
  
  //VERY IMPORTANT: CLEAR ALL VECTORS
  if(!fIsMCGenTree) {
    ResetDmesonCommonVarVectors();
    fInvMassK0s.clear();
    fcTauK0s.clear();
    fV0PointingAngle.clear();
    fCosThetaStar.clear();
    fsignd0.clear();
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
