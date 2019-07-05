/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//*************************************************************************
// \class AliHFTreeHandlerLctopKpi
// \brief helper class to handle a tree for Lc->pKpi cut optimisation and MVA analyses
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
#include "AliHFTreeHandlerLctopKpi.h"
#include "AliAODRecoDecayHF3Prong.h"

/// \cond CLASSIMP
ClassImp(AliHFTreeHandlerLctopKpi);
/// \endcond

//________________________________________________________________
AliHFTreeHandlerLctopKpi::AliHFTreeHandlerLctopKpi():
  AliHFTreeHandler(),
  fSigmaVertex(-9999.),
  fDist12toPrim(-9999.),
  fDist23toPrim(-9999.),
  fNormd0MeasMinusExp(-9999.),
  fSumImpParProngs(-9999.)
{
  //
  // Default constructor
  //

  fNProngs=3; // --> cannot be changed
  for(unsigned int iProng=0; iProng<fNProngs; iProng++){
    fImpParProng[iProng] = -9999.;
    fDCAProng[iProng] = -9999.;
  }
}

//________________________________________________________________
AliHFTreeHandlerLctopKpi::AliHFTreeHandlerLctopKpi(int PIDopt):
  AliHFTreeHandler(PIDopt),
  fSigmaVertex(-9999.),
  fDist12toPrim(-9999.),
  fDist23toPrim(-9999.),
  fNormd0MeasMinusExp(-9999.),
  fSumImpParProngs(-9999.)
{
  //
  // Standard constructor
  //

  fNProngs=3; // --> cannot be changed
  for(unsigned int iProng=0; iProng<fNProngs; iProng++){
    fImpParProng[iProng] = -9999.;
    fDCAProng[iProng] = -9999.;
  }
}

//________________________________________________________________
AliHFTreeHandlerLctopKpi::~AliHFTreeHandlerLctopKpi()
{
  //
  // Default Destructor
  //
}

//________________________________________________________________
TTree* AliHFTreeHandlerLctopKpi::BuildTree(TString name, TString title) 
{
  fIsMCGenTree=false;

  if(fTreeVar) {
    delete fTreeVar;
    fTreeVar=nullptr;
  }
  fTreeVar = new TTree(name.Data(),title.Data());

  //set common variables
  AddCommonDmesonVarBranches();

  //set Lc variables
  fTreeVar->Branch("sig_vert",&fSigmaVertex);
  fTreeVar->Branch("dist_12",&fDist12toPrim);
  fTreeVar->Branch("dist_23",&fDist23toPrim);
  fTreeVar->Branch("max_norm_d0d0exp",&fNormd0MeasMinusExp);
  fTreeVar->Branch("sum_d0d0_prongs",&fSumImpParProngs);
  for(unsigned int iProng=0; iProng<fNProngs; iProng++){
    fTreeVar->Branch(Form("imp_par_prong%d",iProng),&fImpParProng[iProng]);
    fTreeVar->Branch(Form("dca_prong%d",iProng),&fDCAProng[iProng]);
  }

  //set single-track variables
  AddSingleTrackBranches();

  //set PID variables
  if(fPidOpt!=kNoPID) AddPidBranches(true,true,true,true,true);

  return fTreeVar;
}

//________________________________________________________________
bool AliHFTreeHandlerLctopKpi::SetVariables(int runnumber, unsigned int eventID, float ptgen, AliAODRecoDecayHF* cand, float bfield, int masshypo, AliPIDResponse *pidrespo)
{
  if(!cand) return false;
  if(fFillOnlySignal) { //if fill only signal and not signal candidate, do not store
    if(!(fCandType&kSignal)) return true;
  }
  fRunNumber=runnumber;
  fEvID=eventID;
  fPtGen=ptgen;
  
  //topological variables
  //common
  fPt=cand->Pt();
  fY=cand->Y(4122);
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
  fSumImpParProngs=cand->Getd0Prong(0)*cand->Getd0Prong(0)+cand->Getd0Prong(1)*cand->Getd0Prong(1)+cand->Getd0Prong(2)*cand->Getd0Prong(2);
  
  //Lc -> pKpi variables
  if(masshypo==1){ //pKpi
    fInvMass=((AliAODRecoDecayHF3Prong*)cand)->InvMassLcpKpi();
  }
  else if(masshypo==2){ //piKp
    fInvMass=((AliAODRecoDecayHF3Prong*)cand)->InvMassLcpiKp();
  }
  else return false;
  fSigmaVertex=((AliAODRecoDecayHF3Prong*)cand)->GetSigmaVert();
  fDist12toPrim=((AliAODRecoDecayHF3Prong*)cand)->GetDist12toPrim();
  fDist23toPrim=((AliAODRecoDecayHF3Prong*)cand)->GetDist23toPrim();
  for(unsigned int iProng=0; iProng<fNProngs; iProng++) {
    fImpParProng[iProng]=cand->Getd0Prong(iProng);
    fDCAProng[iProng]=cand->GetDCA(iProng);
  }
    
  //single track variables
  AliAODTrack* prongtracks[3];
  for(unsigned int iProng=0; iProng<fNProngs; iProng++) prongtracks[iProng] = (AliAODTrack*)cand->GetDaughter(iProng);
  bool setsingletrack = SetSingleTrackVars(prongtracks);  
  if(!setsingletrack) return false;

  //pid variables
  if(fPidOpt==kNoPID) return true;

  bool setpid = SetPidVars(prongtracks,pidrespo,true,true,true,true,true);
  if(!setpid) return false;

  return true;
}
