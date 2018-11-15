/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//*************************************************************************
// \class AliHFTreeHandlerDstoKKpi
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
#include <TDatabasePDG.h>
#include "AliHFTreeHandlerDstoKKpi.h"
#include "AliAODRecoDecayHF3Prong.h"

/// \cond CLASSIMP
ClassImp(AliHFTreeHandlerDstoKKpi);
/// \endcond

//________________________________________________________________
AliHFTreeHandlerDstoKKpi::AliHFTreeHandlerDstoKKpi():
  AliHFTreeHandler(),
  fImpParProng(),
  fSigmaVertex(),
  fMassKK(),
  fCosPiDs(),
  fCosPiKPhi(),
  fNormd0MeasMinusExp(),
  fMassKKOpt(kMassKK)
{
  //
  // Default constructor
  //

  fNProngs=3; // --> cannot be changed
}

//________________________________________________________________
AliHFTreeHandlerDstoKKpi::AliHFTreeHandlerDstoKKpi(int PIDopt):
  AliHFTreeHandler(PIDopt),
  fImpParProng(),
  fSigmaVertex(),
  fMassKK(),
  fCosPiDs(),
  fCosPiKPhi(),
  fNormd0MeasMinusExp(),
  fMassKKOpt(kMassKK)
{
  //
  // Standard constructor
  //

  fNProngs=3; // --> cannot be changed
}

//________________________________________________________________
AliHFTreeHandlerDstoKKpi::~AliHFTreeHandlerDstoKKpi()
{
  //
  // Default Destructor
  //
}

//________________________________________________________________
TTree* AliHFTreeHandlerDstoKKpi::BuildTree(TString name, TString title) 
{
  fIsMCGenTree=false;

  if(fTreeVar) {
    delete fTreeVar;
    fTreeVar=0x0;
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
bool AliHFTreeHandlerDstoKKpi::SetVariables(AliAODRecoDecayHF* cand, float bfield, int masshypo, AliAODPidHF* pidHF) 
{
  if(!cand) return false;
  if(fFillOnlySignal) { //if fill only signal and not signal candidate, do not store
    if(!(fCandTypeMap&kSignal || fCandTypeMap&kRefl)) return true;
  }
  fNCandidates++;

  //topological variables
  //common
  fCandType.push_back(fCandTypeMap);
  fPt.push_back(cand->Pt());
  fY.push_back(cand->Y(431));
  fEta.push_back(cand->Eta());
  fPhi.push_back(cand->Phi());
  fDecayLength.push_back(cand->DecayLength());
  fDecayLengthXY.push_back(cand->DecayLengthXY());
  fNormDecayLengthXY.push_back(cand->NormalizedDecayLengthXY());
  fCosP.push_back(cand->CosPointingAngle());
  fCosPXY.push_back(cand->CosPointingAngleXY());
  fImpParXY.push_back(cand->ImpParXY());
  fNormd0MeasMinusExp.push_back(ComputeMaxd0MeasMinusExp(cand,bfield));

  //Ds+ -> KKpi variables
  fSigmaVertex.push_back(((AliAODRecoDecayHF3Prong*)cand)->GetSigmaVert());
  float massPhi = 0;
  float cospikphi=-2;
  if(fMassKKOpt==kDeltaMassKKPhi) massPhi = TDatabasePDG::Instance()->GetParticle(333)->Mass();
  if(masshypo==0){ //phiKKpi
    fInvMass.push_back(((AliAODRecoDecayHF3Prong*)cand)->InvMassDsKKpi());
    fMassKK.push_back(TMath::Abs(((AliAODRecoDecayHF3Prong*)cand)->InvMass2Prongs(0,1,321,321)-massPhi));
    fCosPiDs.push_back(((AliAODRecoDecayHF3Prong*)cand)->CosPiDsLabFrameKKpi());
    cospikphi = ((AliAODRecoDecayHF3Prong*)cand)->CosPiKPhiRFrameKKpi();
  }
  else if(masshypo==1){ //phipiKK
    fInvMass.push_back(((AliAODRecoDecayHF3Prong*)cand)->InvMassDspiKK());
    fMassKK.push_back(TMath::Abs(((AliAODRecoDecayHF3Prong*)cand)->InvMass2Prongs(1,2,321,321)-massPhi));
    fCosPiDs.push_back(((AliAODRecoDecayHF3Prong*)cand)->CosPiDsLabFramepiKK());
    cospikphi = ((AliAODRecoDecayHF3Prong*)cand)->CosPiKPhiRFramepiKK();
  }
  fCosPiKPhi.push_back(cospikphi*cospikphi*cospikphi);
  for(unsigned int iProng=0; iProng<fNProngs; iProng++) {
    fImpParProng[iProng].push_back(cand->Getd0Prong(iProng));
  }
    
  //single-track variables
  AliAODTrack* prongtracks[3];
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
void AliHFTreeHandlerDstoKKpi::FillTree() {
  fTreeVar->Fill();
  
  //VERY IMPORTANT: CLEAR ALL VECTORS
  if(!fIsMCGenTree) {
    ResetDmesonCommonVarVectors();
    fSigmaVertex.clear();
    fMassKK.clear();
    fCosPiDs.clear();
    fCosPiKPhi.clear();
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
