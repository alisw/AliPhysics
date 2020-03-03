/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//*************************************************************************
// \class AliHFTreeHandlerDplustoKpipi
// \brief helper class to handle a tree for D+ cut optimisation and MVA analyses
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
#include "AliHFTreeHandlerDplustoKpipi.h"
#include "AliAODRecoDecayHF3Prong.h"

/// \cond CLASSIMP
ClassImp(AliHFTreeHandlerDplustoKpipi);
/// \endcond

//________________________________________________________________
AliHFTreeHandlerDplustoKpipi::AliHFTreeHandlerDplustoKpipi():
  AliHFTreeHandler(),
  fSigmaVertex(-9999.),
  fNormd0MeasMinusExp(-9999.)
{
  //
  // Default constructor
  //

  fNProngs=3; // --> cannot be changed
  for(unsigned int iProng=0; iProng<fNProngs; iProng++) 
    fImpParProng[iProng] = -9999.;
}

//________________________________________________________________
AliHFTreeHandlerDplustoKpipi::AliHFTreeHandlerDplustoKpipi(int PIDopt):
  AliHFTreeHandler(PIDopt),
  fSigmaVertex(-9999.),
  fNormd0MeasMinusExp(-9999.)
{
  //
  // Standard constructor
  //

  fNProngs=3; // --> cannot be changed
  for(unsigned int iProng=0; iProng<fNProngs; iProng++) 
    fImpParProng[iProng] = -9999.;
}

//________________________________________________________________
AliHFTreeHandlerDplustoKpipi::~AliHFTreeHandlerDplustoKpipi()
{
  //
  // Default Destructor
  //
}

//________________________________________________________________
TTree* AliHFTreeHandlerDplustoKpipi::BuildTree(TString name, TString title) 
{
  fIsMCGenTree=false;

  if(fTreeVar) {
    delete fTreeVar;
    fTreeVar=nullptr;
  }
  fTreeVar = new TTree(name.Data(),title.Data());

  //set common variables
  AddCommonDmesonVarBranches();

  //set D+ variables
  fTreeVar->Branch("sig_vert",&fSigmaVertex);
  fTreeVar->Branch("max_norm_d0d0exp",&fNormd0MeasMinusExp);
  for(unsigned int iProng=0; iProng<fNProngs; iProng++){
    fTreeVar->Branch(Form("imp_par_prong%d",iProng),&fImpParProng[iProng]);
  }

  //set single-track variables
  AddSingleTrackBranches();
  if (fFillJets) AddJetBranches();

  //set PID variables
  if(fPidOpt!=kNoPID) AddPidBranches(true,true,false,true,true);

  return fTreeVar;
}

//________________________________________________________________
bool AliHFTreeHandlerDplustoKpipi::SetVariables(int runnumber, int eventID, int eventID_Ext, Long64_t eventID_Long, float ptgen, AliAODRecoDecayHF* cand, float bfield, int /*masshypo*/, AliPIDResponse *pidrespo) 
{
  fRunNumber=runnumber;
  fEvID=eventID;
  fEvIDExt=eventID_Ext;
  fEvIDLong=eventID_Long;

  if(!cand) return false;
  if(fFillOnlySignal) { //if fill only signal and not signal candidate, do not store
    if(!(fCandType&kSignal)) return true;
  }

  fPtGen=ptgen;
  
  fCandType &= ~kRefl; //protection --> D+ ->Kpipi cannot be reflected

  //topological variables
  //common
  fPt=cand->Pt();
  fY=cand->Y(411);
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

  //D+ -> Kpipi variables
  fInvMass=((AliAODRecoDecayHF3Prong*)cand)->InvMassDplus();
  fSigmaVertex=((AliAODRecoDecayHF3Prong*)cand)->GetSigmaVert();
  for(unsigned int iProng=0; iProng<fNProngs; iProng++) {
    fImpParProng[iProng]=cand->Getd0Prong(iProng);
  }
    
  //single track variables
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
