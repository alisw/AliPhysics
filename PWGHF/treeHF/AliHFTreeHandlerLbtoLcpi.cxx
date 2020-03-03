/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//*************************************************************************
// \class AliHFTreeHandlerLbtoLcpi
// \brief helper class to handle a tree for Lb cut optimisation and MVA analyses
// \authors:
// D. Andreou, dimitra.andreou@cern.ch
// L. Vermunt, luuk.vermunt@cern.ch
// G. Innocenti, gian.michele.innocenti@cern.ch
/////////////////////////////////////////////////////////////

#include <TString.h>
#include <TDatabasePDG.h>
#include "AliHFTreeHandlerLbtoLcpi.h"
#include "AliAODRecoDecayHF3Prong.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliESDVertex.h"
#include "AliESDtrackCuts.h"
#include "AliRDHFCutsLctopKpi.h"

/// \cond CLASSIMP
ClassImp(AliHFTreeHandlerLbtoLcpi);
/// \endcond

//________________________________________________________________
AliHFTreeHandlerLbtoLcpi::AliHFTreeHandlerLbtoLcpi():
  AliHFTreeHandler(),
  fCosThetaStar(-9999.),
  fImpParProd(-9999.),
  fcTau(-9999.),
  fInvMass_Lc(-9999.),
  fImpPar_Lc(-9999.),
  fPt_Lc(-9999.),
  fY_Lc(-9999.),
  fEta_Lc(-9999.),
  fPhi_Lc(-9999.),
  fDecayLength_Lc(-9999.),
  fDecayLengthXY_Lc(-9999.),
  fNormDecayLengthXY_Lc(-9999.),
  fCosP_Lc(-9999.),
  fCosPXY_Lc(-9999.),
  fImpParXY_Lc(-9999.),
  fDCA_Lc(-9999.),
  fSigmaVertex_Lc(-9999.),
  fDist12toPrim_Lc(-9999.),
  fDist23toPrim_Lc(-9999.),
  fNormd0MeasMinusExp_Lc(-9999.),
  fSumImpParProngs_Lc(-9999.),
  fChi2OverNDF(-9999.)
{
  //
  // Default constructor
  //

  //Only used in for-loops with a self-made AliAODtrack vector, so "Lb pion + 3 Lc-prongs" is fine
  fNProngs=4; // --> cannot be changed
  for(unsigned int iProng=0; iProng<fNProngs; iProng++) 
    fImpParProng[iProng] = -9999.;
  
  //Lc specific variable, so only three prongs
  for(unsigned int iProng=0; iProng<3; iProng++)
    fDCAProng_Lc[iProng] = -9999.;
}

//________________________________________________________________
AliHFTreeHandlerLbtoLcpi::AliHFTreeHandlerLbtoLcpi(int PIDopt):
  AliHFTreeHandler(PIDopt),
  fCosThetaStar(-9999.),
  fImpParProd(-9999.),
  fcTau(-9999.),
  fInvMass_Lc(-9999.),
  fImpPar_Lc(-9999.),
  fPt_Lc(-9999.),
  fY_Lc(-9999.),
  fEta_Lc(-9999.),
  fPhi_Lc(-9999.),
  fDecayLength_Lc(-9999.),
  fDecayLengthXY_Lc(-9999.),
  fNormDecayLengthXY_Lc(-9999.),
  fCosP_Lc(-9999.),
  fCosPXY_Lc(-9999.),
  fImpParXY_Lc(-9999.),
  fDCA_Lc(-9999.),
  fSigmaVertex_Lc(-9999.),
  fDist12toPrim_Lc(-9999.),
  fDist23toPrim_Lc(-9999.),
  fNormd0MeasMinusExp_Lc(-9999.),
  fSumImpParProngs_Lc(-9999.),
  fChi2OverNDF(-9999.)
{
  //
  // Standard constructor
  //
    
  //Only used in for-loops with a self-made AliAODtrack vector, so "Lb pion + 3 Lc-prongs" is fine
  fNProngs=4; // --> cannot be changed
  for(unsigned int iProng=0; iProng<fNProngs; iProng++)
    fImpParProng[iProng] = -9999.;
  
  //Lc specific variable, so only three prongs
  for(unsigned int iProng=0; iProng<3; iProng++)
    fDCAProng_Lc[iProng] = -9999.;
}

//________________________________________________________________
AliHFTreeHandlerLbtoLcpi::~AliHFTreeHandlerLbtoLcpi()
{
  //
  // Default Destructor
  //
}

//________________________________________________________________
TTree* AliHFTreeHandlerLbtoLcpi::BuildTree(TString name, TString title)
{
  fIsMCGenTree=false;
  
  if(fTreeVar) {
    delete fTreeVar;
    fTreeVar=nullptr;
  }
  fTreeVar = new TTree(name.Data(),title.Data());

  //set common variables
  AddCommonDmesonVarBranches();

  //set Lb variables
  fTreeVar->Branch("cos_t_star",&fCosThetaStar);
  fTreeVar->Branch("imp_par_prod",&fImpParProd);
  fTreeVar->Branch("ctau",&fcTau);
  fTreeVar->Branch("chi2_over_ndf",&fChi2OverNDF);
  for(unsigned int iProng=0; iProng<fNProngs; iProng++){
    fTreeVar->Branch(Form("imp_par_prong%d",iProng),&fImpParProng[iProng]);
  }

  //set Lc variables
  fTreeVar->Branch("inv_mass_Lc",&fInvMass_Lc);
  fTreeVar->Branch("imp_par_Lc", &fImpPar_Lc);
  fTreeVar->Branch("pt_Lc",&fPt_Lc);
  fTreeVar->Branch("y_Lc",&fY_Lc);
  fTreeVar->Branch("eta_Lc",&fEta_Lc);
  fTreeVar->Branch("phi_Lc",&fPhi_Lc);
  fTreeVar->Branch("d_len_Lc",&fDecayLength_Lc);
  fTreeVar->Branch("d_len_xy_Lc",&fDecayLengthXY_Lc);
  fTreeVar->Branch("norm_dl_xy_Lc",&fNormDecayLengthXY_Lc);
  fTreeVar->Branch("cos_p_Lc",&fCosP_Lc);
  fTreeVar->Branch("cos_p_xy_Lc",&fCosPXY_Lc);
  fTreeVar->Branch("imp_par_xy_Lc",&fImpParXY_Lc);
  fTreeVar->Branch("dca_Lc",&fDCA_Lc);
  fTreeVar->Branch("sig_vert_Lc",&fSigmaVertex_Lc);
  fTreeVar->Branch("dist_12_Lc",&fDist12toPrim_Lc);
  fTreeVar->Branch("dist_23_Lc",&fDist23toPrim_Lc);
  fTreeVar->Branch("max_norm_d0d0exp_Lc",&fNormd0MeasMinusExp_Lc);
  fTreeVar->Branch("sum_d0d0_prongs_Lc",&fSumImpParProngs_Lc);
  for(unsigned int iProng=0; iProng<3; iProng++){
    fTreeVar->Branch(Form("dca_prong%d_Lc",iProng),&fDCAProng_Lc[iProng]);
  }

  //set single-track variables
  AddSingleTrackBranches();
  if (fFillJets) AddJetBranches();

  //set PID variables
  if(fPidOpt!=kNoPID) AddPidBranches(true,true,true,true,true);

  return fTreeVar;
}

//________________________________________________________________
bool AliHFTreeHandlerLbtoLcpi::SetVariables(int runnumber, int eventID, int eventID_Ext, Long64_t eventID_Long, float ptgen, AliAODRecoDecayHF* cand, float bfield, int masshypo/*used for Lc*/, AliPIDResponse* pidrespo)
{
  fIsMCGenTree=false;

  if(!cand) return false;
  if(fFillOnlySignal) { //if fill only signal and not signal candidate, do not store
    if(!(fCandType&kSignal)) return true;
  }
  fRunNumber=runnumber;
  fEvID=eventID;
  fEvIDExt=eventID_Ext;
  fEvIDLong=eventID_Long;
  fPtGen=ptgen;
 
  AliAODRecoDecayHF3Prong* candLc = (AliAODRecoDecayHF3Prong*)cand->GetDaughter(0); //Lc
  AliAODVertex* vtxLb = candLc->GetSecondaryVtx();
  
  //topological variables
  //common (Lb -> Lc pi)
  fPt=((AliAODRecoDecayHF2Prong*)cand)->Pt();
  fY=((AliAODRecoDecayHF2Prong*)cand)->Y(5122);
  fEta=((AliAODRecoDecayHF2Prong*)cand)->Eta();
  fPhi=((AliAODRecoDecayHF2Prong*)cand)->Phi();
  fDecayLength=((AliAODRecoDecayHF2Prong*)cand)->DecayLength();
  fDecayLengthXY=((AliAODRecoDecayHF2Prong*)cand)->DecayLengthXY();
  fNormDecayLengthXY=((AliAODRecoDecayHF2Prong*)cand)->NormalizedDecayLengthXY();
  fCosP=((AliAODRecoDecayHF2Prong*)cand)->CosPointingAngle();
  fCosPXY=((AliAODRecoDecayHF2Prong*)cand)->CosPointingAngleXY();
  fImpParXY=((AliAODRecoDecayHF2Prong*)cand)->ImpParXY();
  fDCA=((AliAODRecoDecayHF2Prong*)cand)->GetDCA();
  fCosThetaStar=cand->CosThetaStar(0,5122,4122,211);
  fImpParProd=cand->Prodd0d0();
  fcTau=cand->Ct(5122);
  fChi2OverNDF = vtxLb->GetChi2perNDF();
  
  UInt_t prongs[2];
  prongs[0] = 4122; prongs[1] = 211;
  fInvMass=((AliAODRecoDecayHF2Prong*)cand)->InvMass(2,prongs);
  
  for(unsigned int iProng=0; iProng<3; iProng++) {
    fImpParProng[iProng]=candLc->Getd0Prong(iProng);
  }
  fImpParProng[3]=cand->Getd0Prong(1);
  fImpPar_Lc=cand->Getd0Prong(0);

  //Lc -> p K pi variables
  fPt_Lc=candLc->Pt();
  fY_Lc=candLc->Y(4122);
  fEta_Lc=candLc->Eta();
  fPhi_Lc=candLc->Phi();
  fDecayLength_Lc=candLc->DecayLength();
  fDecayLengthXY_Lc=candLc->DecayLengthXY();
  fNormDecayLengthXY_Lc=candLc->NormalizedDecayLengthXY();
  fCosP_Lc=candLc->CosPointingAngle();
  fCosPXY_Lc=candLc->CosPointingAngleXY();
  fImpParXY_Lc=candLc->ImpParXY();
  fDCA_Lc=candLc->GetDCA();
  fSigmaVertex_Lc=candLc->GetSigmaVert();
  fDist12toPrim_Lc=candLc->GetDist12toPrim();
  fDist23toPrim_Lc=candLc->GetDist23toPrim();
  fNormd0MeasMinusExp_Lc=ComputeMaxd0MeasMinusExp(candLc,bfield);
  fSumImpParProngs_Lc=candLc->Getd0Prong(0)*candLc->Getd0Prong(0)+candLc->Getd0Prong(1)*candLc->Getd0Prong(1)+candLc->Getd0Prong(2)*candLc->Getd0Prong(2);
  
  if(masshypo==1){ //pKpi
    fInvMass_Lc=candLc->InvMassLcpKpi();
  } else{ //piKp
    fInvMass_Lc=candLc->InvMassLcpiKp();
  }
  
  for(unsigned int iProng=0; iProng<3; iProng++) {
    fDCAProng_Lc[iProng]=candLc->GetDCA(iProng);
  }

  //single track variables
  AliAODTrack* prongtracks[4];
  for(unsigned int iProng=0; iProng<3; iProng++) prongtracks[iProng] = (AliAODTrack*)candLc->GetDaughter(iProng);
  prongtracks[3] = (AliAODTrack*)cand->GetDaughter(1);

  bool setsingletrack = SetSingleTrackVars(prongtracks);
  if(!setsingletrack) return false;

  //pid variables
  if(fPidOpt==kNoPID) return true;

  bool setpid = SetPidVars(prongtracks,pidrespo,true,true,true,true,true);
  if(!setpid) return false;

  return true;
}

//________________________________________________________________
Int_t AliHFTreeHandlerLbtoLcpi::IsLbPionSelected(TObject* obj, AliRDHFCutsLctopKpi* cutsLc, AliAODPidHF* fPidHFLc, AliAODEvent* aod, AliAODVertex *vtx) {

  AliAODTrack* candidatePion = (AliAODTrack*)obj;
  if (!candidatePion){ AliWarning("No pion object. Track rejected."); return 0; }
  
  AliESDtrackCuts* fTrackCuts = cutsLc->GetTrackCuts();
  if (!fTrackCuts){ AliWarning("No fTrackCuts object. Track rejected.");  return 0; }

  Double_t pos[3],cov[6];
  vtx->GetXYZ(pos);
  vtx->GetCovarianceMatrix(cov);
  const AliESDVertex vESD(pos,cov,100.,100);
  
  if(!cutsLc->IsDaughterSelected(candidatePion,&vESD,fTrackCuts,aod)) return 0;

  if(!cutsLc->GetIsUsePID()) return 1;
  else {
    
    if(!fPidHFLc){ AliWarning("AliAODPidHF not created. Track accepted"); return 1; }
    
    Int_t isPion=fPidHFLc->MakeRawPid(candidatePion,AliPID::kPion);
    if(isPion) return 1;
    else       return 0;
  }
  
  return 1;
}
