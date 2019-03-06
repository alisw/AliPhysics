/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//*************************************************************************
// \class AliHFTreeHandlerDstartoKpipi
// \brief helper class to handle a tree for Dstar cut optimisation and MVA analyses
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
#include "AliHFTreeHandlerDstartoKpipi.h"
#include "AliAODRecoCascadeHF.h"

/// \cond CLASSIMP
ClassImp(AliHFTreeHandlerDstartoKpipi);
/// \endcond

//________________________________________________________________
AliHFTreeHandlerDstartoKpipi::AliHFTreeHandlerDstartoKpipi():
  AliHFTreeHandler(),
  fImpParProng(),
  fCosThetaStar(),
  fImpParProd(),
  fNormd0MeasMinusExp(),
  fAngleD0dkpPisoft(),
  fInvMass_D0(),
  fPt_D0(),
  fY_D0(),
  fEta_D0(),
  fPhi_D0()
{
  //
  // Default constructor
  //

  //Only used in for-loops with a self-made AliAODtrack vector, so "Dstar pion + 2 D0-prongs" is fine
  fNProngs=3; // --> cannot be changed
}

//________________________________________________________________
AliHFTreeHandlerDstartoKpipi::AliHFTreeHandlerDstartoKpipi(int PIDopt):
  AliHFTreeHandler(PIDopt),
  fImpParProng(),
  fCosThetaStar(),
  fImpParProd(),
  fNormd0MeasMinusExp(),
  fAngleD0dkpPisoft(),
  fInvMass_D0(),
  fPt_D0(),
  fY_D0(),
  fEta_D0(),
  fPhi_D0()
{
  //
  // Standard constructor
  //

  //Only used in for-loops with a self-made AliAODtrack vector, so "Dstar pion + 2 D0-prongs" is fine
  fNProngs=3; // --> cannot be changed
}

//________________________________________________________________
AliHFTreeHandlerDstartoKpipi::~AliHFTreeHandlerDstartoKpipi()
{
  //
  // Default Destructor
  //
}

//________________________________________________________________
TTree* AliHFTreeHandlerDstartoKpipi::BuildTree(TString name, TString title)
{
  fIsMCGenTree=false;

  if(fTreeVar) {
    delete fTreeVar;
    fTreeVar=0x0;
  }
  fTreeVar = new TTree(name.Data(),title.Data());

  //set common variables (please pay attention, filled with D0 and Dstar variables)
  AddCommonDmesonVarBranches();

  //set Dstar variables
  fTreeVar->Branch("cos_t_star",&fCosThetaStar);
  fTreeVar->Branch("imp_par_prod",&fImpParProd);
  fTreeVar->Branch("max_norm_d0d0exp",&fNormd0MeasMinusExp);
  fTreeVar->Branch("angle_D0dkpPisoft",&fAngleD0dkpPisoft);
  for(unsigned int iProng=0; iProng<fNProngs; iProng++){
    fTreeVar->Branch(Form("imp_par_prong%d",iProng),&fImpParProng[iProng]);
  }

  //set D0 variables
  fTreeVar->Branch("inv_mass_D0",&fInvMass_D0);
  fTreeVar->Branch("pt_D0",&fPt_D0);
  fTreeVar->Branch("y_D0",&fY_D0);
  fTreeVar->Branch("eta_D0",&fEta_D0);
  fTreeVar->Branch("phi_D0",&fPhi_D0);

  //set single-track variables
  AddSingleTrackBranches();

  //set PID variables
  if(fPidOpt!=kNoPID) AddPidBranches(true,true,false,true,true);

  return fTreeVar;
}

//________________________________________________________________
bool AliHFTreeHandlerDstartoKpipi::SetVariables(int runnumber, unsigned int eventID, AliAODRecoDecayHF* cand, float bfield, int /*masshypo*/, AliPIDResponse *pidrespo)
{
  fIsMCGenTree=false;

  if(!cand) return false;
  if(fFillOnlySignal) { //if fill only signal and not signal candidate, do not store
    if(!(fCandTypeMap&kSignal)) return true;
  }
  fNCandidates++;
  fRunNumber.push_back(runnumber);
  fEvID.push_back(eventID);
  
  fCandTypeMap &= ~kRefl; //protection --> Dstar ->Kpipi cannot be reflected

  AliAODRecoDecayHF2Prong *d0 = ((AliAODRecoCascadeHF*)cand)->Get2Prong();

  //topological variables (Dstar and D0 variables combined)
  //common (Dstar)
  fCandType.push_back(fCandTypeMap);
  fCandTypeMap=0; //reset candtype
  fPt.push_back(((AliAODRecoCascadeHF*)cand)->Pt());
  fY.push_back(((AliAODRecoCascadeHF*)cand)->YDstar());
  fEta.push_back(((AliAODRecoCascadeHF*)cand)->Eta());
  fPhi.push_back(((AliAODRecoCascadeHF*)cand)->Phi());
  //common (D0)
  fDecayLength.push_back(d0->DecayLength());
  fDecayLengthXY.push_back(d0->DecayLengthXY());
  fNormDecayLengthXY.push_back(d0->NormalizedDecayLengthXY()*(d0->P()/d0->Pt()));
  fCosP.push_back(d0->CosPointingAngle());
  fCosPXY.push_back(d0->CosPointingAngleXY());
  fImpParXY.push_back(d0->ImpParXY());
  fDCA.push_back(d0->GetDCA());
  fNormd0MeasMinusExp.push_back(ComputeMaxd0MeasMinusExp(d0,bfield));
  fAngleD0dkpPisoft.push_back(((AliAODRecoCascadeHF*)cand)->AngleD0dkpPisoft());

  AliAODTrack* prongtracks[3];
  prongtracks[0] = (AliAODTrack*)((AliAODRecoCascadeHF*)cand)->GetBachelor();
  fImpParProng[0].push_back(cand->Getd0Prong(0));
    
  //D* -> D0 pi variables
  fInvMass.push_back(((AliAODRecoCascadeHF*)cand)->DeltaInvMass());
  fImpParProd.push_back(d0->Prodd0d0());
  if( (((AliAODRecoCascadeHF*)cand)->Charge()) >0)   fCosThetaStar.push_back(d0->CosThetaStarD0());
  else   fCosThetaStar.push_back(d0->CosThetaStarD0bar());

  //D0 -> K pi variables
  fPt_D0.push_back(d0->Pt());
  fY_D0.push_back(d0->Y(421));
  fEta_D0.push_back(d0->Eta());
  fPhi_D0.push_back(d0->Phi());
  if( (((AliAODRecoCascadeHF*)cand)->Charge()) > 0) {
    fInvMass_D0.push_back(d0->InvMassD0());
      
    fImpParProng[1].push_back(d0->Getd0Prong(0));
    fImpParProng[2].push_back(d0->Getd0Prong(1));
    prongtracks[1] = (AliAODTrack*)d0->GetDaughter(0);
    prongtracks[2] = (AliAODTrack*)d0->GetDaughter(1);
  } else {
    fInvMass_D0.push_back(d0->InvMassD0bar());

    fImpParProng[1].push_back(d0->Getd0Prong(1));
    fImpParProng[2].push_back(d0->Getd0Prong(0));
    prongtracks[1] = (AliAODTrack*)d0->GetDaughter(1);
    prongtracks[2] = (AliAODTrack*)d0->GetDaughter(0);
  }
    
  //single track variables
  bool setsingletrack = SetSingleTrackVars(prongtracks);
  if(!setsingletrack) return false;

  //pid variables
  if(fPidOpt==kNoPID) return true;

  bool setpid = SetPidVars(prongtracks,pidrespo,true,true,false,true,true);
  if(!setpid) return false;

  return true;
}

//________________________________________________________________
void AliHFTreeHandlerDstartoKpipi::FillTree() {
  fTreeVar->Fill();
  
  //VERY IMPORTANT: CLEAR ALL VECTORS
  if(!fIsMCGenTree) {
    ResetDmesonCommonVarVectors();
    fCosThetaStar.clear();
    fImpParProd.clear();
    fNormd0MeasMinusExp.clear();
    fAngleD0dkpPisoft.clear();
    fInvMass_D0.clear();
    fPt_D0.clear();
    fY_D0.clear();
    fEta_D0.clear();
    fPhi_D0.clear();
    for(unsigned int iProng=0; iProng<fNProngs; iProng++) fImpParProng[iProng].clear();
    ResetSingleTrackVarVectors();
    if(fPidOpt!=kNoPID) ResetPidVarVectors();
  } else {
    ResetMCGenVectors();
  }
  fCandTypeMap=0;
  fNCandidates=0;
}
