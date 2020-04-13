/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//*************************************************************************
// \class AliHFTreeHandlerD0toKpi
// \brief helper class to handle a tree for D0 cut optimisation and MVA analyses
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
#include "AliHFTreeHandlerD0toKpi.h"
#include "AliAODRecoDecayHF2Prong.h"

/// \cond CLASSIMP
ClassImp(AliHFTreeHandlerD0toKpi);
/// \endcond

//________________________________________________________________
AliHFTreeHandlerD0toKpi::AliHFTreeHandlerD0toKpi():
  AliHFTreeHandler(),
  fCosThetaStar(-9999.),
  fImpParProd(-9999.),
  fNormd0MeasMinusExp(-9999.),
  fNormDecayLength(-9999.)
{
  //
  // Default constructor
  //

  fNProngs=2; // --> cannot be changed
  for(unsigned int iProng=0; iProng<fNProngs; iProng++) {
    fImpParProng[iProng] = -9999.;
    fImpParErrProng[iProng] = -9999.;
  }
}

//________________________________________________________________
AliHFTreeHandlerD0toKpi::AliHFTreeHandlerD0toKpi(int PIDopt):
  AliHFTreeHandler(PIDopt),
  fCosThetaStar(-9999.),
  fImpParProd(-9999.),
  fNormd0MeasMinusExp(-9999.),
  fNormDecayLength(-9999.)
{
  //
  // Standard constructor
  //

  fNProngs=2; // --> cannot be changed
  for(unsigned int iProng=0; iProng<fNProngs; iProng++) { 
    fImpParProng[iProng] = -9999.;
    fImpParErrProng[iProng] = -9999.;
  }
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
    fTreeVar=nullptr;
  }
  fTreeVar = new TTree(name.Data(),title.Data());

  //set common variables
  AddCommonDmesonVarBranches();

  //set D0 variables
  fTreeVar->Branch("cos_t_star",&fCosThetaStar);
  fTreeVar->Branch("imp_par_prod",&fImpParProd);
  fTreeVar->Branch("max_norm_d0d0exp",&fNormd0MeasMinusExp);
  fTreeVar->Branch("norm_dl",&fNormDecayLength);
  for(unsigned int iProng=0; iProng<fNProngs; iProng++){
    fTreeVar->Branch(Form("imp_par_prong%d",iProng),&fImpParProng[iProng]);
    fTreeVar->Branch(Form("imp_par_err_prong%d",iProng),&fImpParErrProng[iProng]);
  }

  //set single-track variables
  AddSingleTrackBranches();
  if (fFillJets) AddJetBranches();
  
  //set PID variables
  if(fPidOpt!=kNoPID) AddPidBranches(true,true,false,true,true);

  return fTreeVar;
}

//________________________________________________________________
bool AliHFTreeHandlerD0toKpi::SetVariables(int runnumber, int eventID, int eventID_Ext, Long64_t eventID_Long, float ptgen, AliAODRecoDecayHF* cand, float bfield, int masshypo, AliPIDResponse *pidrespo) 
{
  fRunNumber=runnumber;
  fEvID=eventID;
  fEvIDExt=eventID_Ext;
  fEvIDLong=eventID_Long;
  fIsMCGenTree=false;

  if(!cand) return false;
  if(fFillOnlySignal) { //if fill only signal and not signal candidate, do not store
    if(!(fCandType&kSignal || fCandType&kRefl)) return true;
  }
  fNCandidates++;
  fPtGen=ptgen;
  
  //topological variables
  //common
  fPt=cand->Pt();
  fY=cand->Y(421);
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
  fNormDecayLength=cand->NormalizedDecayLength();
  
  //D0 -> Kpi variables
  fImpParProd=((AliAODRecoDecayHF2Prong*)cand)->Prodd0d0();
  if(masshypo==0) { //D0 -> Kpi
    fInvMass=((AliAODRecoDecayHF2Prong*)cand)->InvMassD0();
    fCosThetaStar=((AliAODRecoDecayHF2Prong*)cand)->CosThetaStarD0();
  }
  else if(masshypo==1) { //D0 -> piK
    fInvMass=((AliAODRecoDecayHF2Prong*)cand)->InvMassD0bar();
    fCosThetaStar=((AliAODRecoDecayHF2Prong*)cand)->CosThetaStarD0bar();
  }
  for(unsigned int iProng=0; iProng<fNProngs; iProng++) {
    fImpParProng[iProng]=cand->Getd0Prong(iProng);
    fImpParErrProng[iProng]=cand->Getd0errProng(iProng);
  }
    
  //single track variables
  AliAODTrack* prongtracks[2];
  for(unsigned int iProng=0; iProng<fNProngs; iProng++) prongtracks[iProng] = (AliAODTrack*)cand->GetDaughter(iProng);
  bool setsingletrack = SetSingleTrackVars(prongtracks);  
  if(!setsingletrack) return false;

  //pid variables
  if(fPidOpt==kNoPID) return true;

  bool setpid = SetPidVars(prongtracks,pidrespo,true,true,false,true,true);
  if(!setpid) return false;

  return true;
}

//________________________________________________________________
void AliHFTreeHandlerD0toKpi::SetIsDzeroDzeroBar(int isSel, int isSelTopo, int isSelPID, int isSelFilt, int isSelTopoFilt, int isSelPIDFilt) {
    
    //analysis cuts
    //combined selection PID topo
    if(isSel==0){
        fCandType &= ~kDzeroComb;
        fCandType &= ~kDzeroBarComb;
    }
    else if(isSel==1){
        fCandType |= kDzeroComb;
        fCandType &= ~kDzeroBarComb;
    }
    else if(isSel==2){
        fCandType |= kDzeroBarComb;
        fCandType &= ~kDzeroComb;
    }
    else if(isSel==3){
        fCandType |= kDzeroComb;
        fCandType |= kDzeroBarComb;
    }
    //topol selection
    if(isSelTopo==0){
        fCandType &= ~kDzeroTopo;
        fCandType &= ~kDzeroBarTopo;
    }
    else if(isSelTopo==1){
        fCandType |= kDzeroTopo;
        fCandType &= ~kDzeroBarTopo;
    }
    else if(isSelTopo==2){
        fCandType |= kDzeroBarTopo;
        fCandType &= ~kDzeroTopo;
    }
    else if(isSelTopo==3){
        fCandType |= kDzeroTopo;
        fCandType |= kDzeroBarTopo;
    }
    //PID selection
    if(isSelPID==0){
        fCandType &= ~kDzeroPID;
        fCandType &= ~kDzeroBarPID;
    }
    else if(isSelPID==1){
        fCandType |= kDzeroPID;
        fCandType &= ~kDzeroBarPID;
    }
    else if(isSelPID==2){
        fCandType |= kDzeroBarPID;
        fCandType &= ~kDzeroPID;
    }
    else if(isSelPID==3){
        fCandType |= kDzeroPID;
        fCandType |= kDzeroBarPID;
    }
    
    //filtering cuts
    //combined selection PID topo
    if(isSelFilt==0){
        fCandType &= ~kDzeroCombFilt;
        fCandType &= ~kDzeroBarCombFilt;
    }
    else if(isSelFilt==1){
        fCandType |= kDzeroCombFilt;
        fCandType &= ~kDzeroBarCombFilt;
    }
    else if(isSelFilt==2){
        fCandType |= kDzeroBarCombFilt;
        fCandType &= ~kDzeroCombFilt;
    }
    else if(isSelFilt==3){
        fCandType |= kDzeroCombFilt;
        fCandType |= kDzeroBarCombFilt;
    }
    //topol selection
    if(isSelTopoFilt==0){
        fCandType &= ~kDzeroTopoFilt;
        fCandType &= ~kDzeroBarTopoFilt;
    }
    else if(isSelTopoFilt==1){
        fCandType |= kDzeroTopoFilt;
        fCandType &= ~kDzeroBarTopoFilt;
    }
    else if(isSelTopoFilt==2){
        fCandType |= kDzeroBarTopoFilt;
        fCandType &= ~kDzeroTopoFilt;
    }
    else if(isSelTopoFilt==3){
        fCandType |= kDzeroTopoFilt;
        fCandType |= kDzeroBarTopoFilt;
    }
    //PID selection
    if(isSelPIDFilt==0){
        fCandType &= ~kDzeroPIDFilt;
        fCandType &= ~kDzeroBarPIDFilt;
    }
    else if(isSelPIDFilt==1){
        fCandType |= kDzeroPIDFilt;
        fCandType &= ~kDzeroBarPIDFilt;
    }
    else if(isSelPIDFilt==2){
        fCandType |= kDzeroBarPIDFilt;
        fCandType &= ~kDzeroPIDFilt;
    }
    else if(isSelPIDFilt==3){
        fCandType |= kDzeroPIDFilt;
        fCandType |= kDzeroBarPIDFilt;
    }    
}
