/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//*************************************************************************
// \class AliHFTreeHandlerD0toKpi
// \brief helper class to handle a tree for D+ cut optimisation and MVA analyses
// \authors:
// F. Catalano, fabio.catalano@cern.ch
// A. Festanti, andrea.festanti@cern.ch
// F. Grosa, fabrizio.grosa@cern.ch
// G. Innocenti, gian.michele.innocenti@cern.ch
// F. Prino, prino@to.infn.it
// L. Vermunt, luuk.vermunt@cern.ch
/////////////////////////////////////////////////////////////

#include <TString.h>
#include "AliHFTreeHandlerD0toKpi.h"
#include "AliAODRecoDecayHF2Prong.h"

/// \cond CLASSIMP
ClassImp(AliHFTreeHandlerD0toKpi);
/// \endcond

//________________________________________________________________
AliHFTreeHandlerD0toKpi::AliHFTreeHandlerD0toKpi():
  AliHFTreeHandler(),
  fImpParProng(),
  fCosThetaStar(),
  fImpParProd(),
  fNormd0MeasMinusExp()
{
  //
  // Default constructor
  //

  fNProngs=2; // --> cannot be changed
}

//________________________________________________________________
AliHFTreeHandlerD0toKpi::AliHFTreeHandlerD0toKpi(int PIDopt):
  AliHFTreeHandler(PIDopt),
  fImpParProng(),
  fCosThetaStar(),
  fImpParProd(),
  fNormd0MeasMinusExp()
{
  //
  // Standard constructor
  //

  fNProngs=2; // --> cannot be changed
}

//________________________________________________________________
AliHFTreeHandlerD0toKpi::~AliHFTreeHandlerD0toKpi()
{
  //
  // Default Destructor
  //
}

//________________________________________________________________
TTree* AliHFTreeHandlerD0toKpi::BuildTree(TString name, TString title) 
{
  fIsMCGenTree=false;
  
  if(fTreeVar) {
    delete fTreeVar;
    fTreeVar=0x0;
  }
  fTreeVar = new TTree(name.Data(),title.Data());

  //set common variables
  AddCommonDmesonVarBranches();

  //set D0 variables
  fTreeVar->Branch("cos_t_star",&fCosThetaStar);
  fTreeVar->Branch("imp_par_prod",&fImpParProd);
  fTreeVar->Branch("max_norm_d0d0exp",&fNormd0MeasMinusExp);
  for(unsigned int iProng=0; iProng<fNProngs; iProng++){
    fTreeVar->Branch(Form("imp_par_prong%d",iProng),&fImpParProng[iProng]);
  }

  //set single-track variables
  AddSingleTrackBranches();

  //set PID variables
  if(fPidOpt!=kNoPID) AddPidBranches(true,true,false,true,true);

  return fTreeVar;
}

//________________________________________________________________
bool AliHFTreeHandlerD0toKpi::SetVariables(AliAODRecoDecayHF* cand, float bfield, int masshypo, AliAODPidHF* pidHF) 
{
  fIsMCGenTree=false;

  if(!cand) return false;
  if(fFillOnlySignal) { //if fill only signal and not signal candidate, do not store
    if(!(fCandTypeMap&kSignal || fCandTypeMap&kRefl)) return true;
  }
  fNCandidates++;

  //topological variables
  //common
  fCandType.push_back(fCandTypeMap);
  fPt.push_back(cand->Pt());
  fY.push_back(cand->Y(421));
  fEta.push_back(cand->Eta());
  fPhi.push_back(cand->Phi());
  fDecayLength.push_back(cand->DecayLength());
  fDecayLengthXY.push_back(cand->DecayLengthXY());
  fNormDecayLengthXY.push_back(cand->NormalizedDecayLengthXY());
  fCosP.push_back(cand->CosPointingAngle());
  fCosPXY.push_back(cand->CosPointingAngleXY());
  fImpParXY.push_back(cand->ImpParXY());
  fNormd0MeasMinusExp.push_back(ComputeMaxd0MeasMinusExp(cand,bfield));
  
  //D0 -> Kpi variables
  fImpParProd.push_back(((AliAODRecoDecayHF2Prong*)cand)->Prodd0d0());
  if(masshypo==0) { //D0 -> Kpi
    fInvMass.push_back(((AliAODRecoDecayHF2Prong*)cand)->InvMassD0());
    fCosThetaStar.push_back(((AliAODRecoDecayHF2Prong*)cand)->CosThetaStarD0());
  }
  else if(masshypo==1) { //D0 -> piK
    fInvMass.push_back(((AliAODRecoDecayHF2Prong*)cand)->InvMassD0bar());
    fCosThetaStar.push_back(((AliAODRecoDecayHF2Prong*)cand)->CosThetaStarD0bar());
  }
  for(unsigned int iProng=0; iProng<fNProngs; iProng++) {
    fImpParProng[iProng].push_back(cand->Getd0Prong(iProng));
  }
    
  //single track variables
  AliAODTrack* prongtracks[2];
  for(unsigned int iProng=0; iProng<fNProngs; iProng++) prongtracks[iProng] = (AliAODTrack*)cand->GetDaughter(iProng);
  bool setsingletrack = SetSingleTrackVars(prongtracks);  
  if(!setsingletrack) return false;

  //pid variables
  if(fPidOpt==kNoPID) return true;

  bool setpid = SetPidVars(prongtracks,pidHF,true,true,false,true,true);
  if(!setpid) return false;

  return true;
}

//________________________________________________________________
void AliHFTreeHandlerD0toKpi::FillTree() {
  fTreeVar->Fill();
  
  //VERY IMPORTANT: CLEAR ALL VECTORS
  if(!fIsMCGenTree) {
    ResetDmesonCommonVarVectors();
    fCosThetaStar.clear();
    fImpParProd.clear();
    fNormd0MeasMinusExp.clear();
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
