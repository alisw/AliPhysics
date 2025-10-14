// std
#include <iostream>
#include <stdio.h>

// ROOT includes
#include <THashList.h>
#include <TChain.h>
#include <TMath.h>
#include <TH1.h>
#include <TFile.h>
#include <TDatabasePDG.h>

// AliRoot includes
#include <AliAnalysisTask.h>
#include <AliAnalysisManager.h>
#include <AliAnalysisUtils.h>
#include <AliAODEvent.h>
#include <AliESDEvent.h>
#include <AliVEvent.h>
#include <AliAODVZERO.h>
#include <AliAODVertex.h>
#include <AliAODInputHandler.h> 
#include <AliAODHandler.h> 
#include <AliAODHeader.h>
#include <AliAODTrack.h> 
#include <AliAODPid.h> 
#include <AliPIDResponse.h>
#include <AliAODv0.h>
#include <AliAODcascade.h>
#include <AliAODTracklets.h>
#include <AliEventPoolManager.h>
//#include <AliCentrality>

//20220524 Flow study
#include <AliMultSelection.h>
#include <AliEventplane.h>
#include <AliQnCorrectionsManager.h>
#include <AliAnalysisTaskFlowVectorCorrections.h>
#include <AliQnCorrectionsQnVector.h>

#include "AliAnalysisTaskpOmegaDibaryon.h"

//class AliAnalysisTaskpOmegaDibaryon;

using namespace std;

ClassImp(AliAnalysisTaskpOmegaDibaryon)

//_____________________________________________________________________________
AliAnalysisTaskpOmegaDibaryon::AliAnalysisTaskpOmegaDibaryon() : 
AliAnalysisTaskSE(),
  fEvent(0x0),
  fESD(0x0),
  fAOD(0),
  fHeader(0),
  fPIDResponse(0),
  ftrigBit(0),
  fEstimator("V0M"),
  fMultSelection(0x0),
  fCentralityMain(0.),
  fCentralityMin(0.),
  fCentralityMax(0.),
  fIsFlowTask(kFALSE),
  fFlowQnVectorMgr(0x0),
//==================== variable
  nevt(0),
//========== Event selection
  isEvtMB(kFALSE),
  isPileupEvt(kFALSE),
  isGoodEvt(kFALSE),
  centralityV0M(-99),
  centralityCL0(-99),
  centralityCL1(-99),
  centralityV0A(-99),
centralityV0C(-99),
centralityZNA(-99),
centralityZNC(-99),
centralityMain(-99),
//========== Event
  nTracks(0),
  nV0s(0),
  nCascades(0),
  vecTarget(),
  kMagF(0),
  PrimVertexKF(),
  vertex(),
  primVertex(),
  primKFVertex(),
//========== Track 
  lorentzsum(),
  particleProton(),
  particlePion(),
  h(),
  np(0),
  nap(0),
  npi(0),
  nppi(0),
  track(),
  dcaxy(-99),
  dcaz(-99),
  dca(-99),
  dcaxyn(-99),
  dcazn(-99),
  dcan(-99),
  dDCA(),
  cDCA(),
  v0Vtx(),
  track_pT(),
  Daugproton(kFALSE),
  Daugpion(kFALSE),
  isProton(kFALSE),
  isAntiProton(kFALSE),
  isNPion(kFALSE),
  isPPion(kFALSE),
  isBachpion(kFALSE),
  isBachkaon(kFALSE),
  isLambdadecayProton(kFALSE),
  isLambdadecayAntiProton(kFALSE),
  isLambdadecayPPion(kFALSE),
  isLambdadecayNPion(kFALSE),
  isBachPPion(kFALSE),
  isBachNPion(kFALSE),
  isBachPKaon(kFALSE),
  isBachNKaon(kFALSE),
  pTCut_Proton(kFALSE),
  pTCut_Pion(kFALSE),
  pTCut_Kaon(kFALSE),
//========== Array 
  ProtonArray(),
  AntiProtonArray(),
  NPionArray(),
  PPionArray(),
  LambdadecayProtonArray(),
  LambdadecayAntiProtonArray(),
  LambdadecayPPionArray(),
  LambdadecayNPionArray(),
  XidecayBachPPionArray(),
  XidecayBachNPionArray(),
  XidecayBachPKaonArray(),
  XidecayBachNKaonArray(),
//========== TrackParam
  exProtonTrack(),
  exPionTrack(),
  exBachPionTrack(),
  exBachKaonTrack(),
//========== KF particle 
  ProtonTrack(),
  PionTrack(),
  KFProtonTrack(),
  KFPPionTrack(),
  KFPionTrack(),
  KFMother(),
  KFMother_K0s(),
  SecVertexKF(),
  SecVertexErrKF(),
  KFsecVertex(),
  fDCADaugter(-99),
  fDcaSecDaughterProton(-99),
  fDcaSecDaughterPion(-99),
  impar(),
  dd(),
  fPA(-99),
  DecLength(-99),
//========== V0
  v0(),
  ptrack(),
  ntrack(),
  dcaxyp_v0(-99),
  dcazp_v0(-99),
  dcap_v0(-99),
  dcaxyn_v0(-99),
  dcazn_v0(-99),
  dcan_v0(-99),
  dDCAp_v0(),
  dDCAn_v0(),
  cDCAp_v0(),
  cDCAn_v0(),
  v0Vtx_v0(),
  energy_v0(-99),
  fcpa_v0(-99),
  ftransradius_v0(-99),
  Daugproton_v0(kFALSE),
  Daugantiproton_v0(kFALSE),
  Daugpion_v0(kFALSE),
  Daugnpion_v0(kFALSE),
  lambda_v0(kFALSE),
  antilambda_v0(kFALSE),
//========== KF particle
  ndp(-99),
  ndap(-99),
  ndpp(-99),
  ndnp(-99),
  nbpp(-99),
  nbnp(-99),
  nbpk(-99),
  nbnk(-99),
//DCAtoPrimVtxXY(-99),
//DCAtoPrimVtxZ(-99),
  DCAtoPrimVtx(-99),
  BachNPionTrack(),
  KFSubMother(),
  KFBachPionTrack(),
  KFBachKaonTrack(),
  KFMother_Omega(),
  SecVertexKF_Xi(),
  SecVertexErrKF_Xi(),
  KFsecVertex_Xi(),
  SecVertexKF_Omega(),
  SecVertexErrKF_Omega(),
  KFsecVertex_Omega(),
  lorentzsum_Xi(),
  particleBachPion(),
  h_Xi(),
  dd_Xi(),
  lorentzsum_Omega(),
  particleBachKaon(),
  fPA_Xi(-99),
  DecLength_Xi(-99),
  DCAdauXY_Xi(-99),
  DCAdauZ_Xi(-99),
  DCAdau_Xi(-99),
  ftrkID_daughter1(-99),
  ftrkID_daughter2(-99),
  ftrkID_daughter3(-99),
//========== Cascade
  casc(),
  pTrackXi(),
  nTrackXi(),
  bachTrackXi(),
  dcaxyp(-99),
  dcazp(-99),
  dcap(-99),
  dcaxyb(-99),
  dcazb(-99),
  dcab(-99),
  dDCAp(),
  cDCAp(),
  dDCAn(),
  cDCAn(),
  dDCAb(),
  cDCAb(),
  v0Vtx_casc(),
  xivertex(),
  lambdacpa(-99),
  lambdatransradius(-99),
  xicpa(-99),
  xitransradius(-99),
  v0vtxx(-99),
  v0vtxy(-99),
  v0vtxz(-99),
  xivtxx(-99),
  xivtxy(-99),
  xivtxz(-99),
//TTree_omegam
  fTree_omegam(0),
  fevtID_omegam(-99),
  fzvertex_omegam(-99),
  fcentrality_omegam(-99),
//Daughter track                                                                                                     
  fdca_daughter1_omegam(-99),
  fdca_daughter2_omegam(-99),
  fdca_daughter3_omegam(-99),
  ftrkID_daughter1_omegam(-99),
  ftrkID_daughter2_omegam(-99),
  ftrkID_daughter3_omegam(-99),
//Lambda                                                                                                             
  fcpa_lambda_omegam(-99),
  fpcaxy_lambda_omegam(-99),
  fpca_lambda_omegam(-99),
  fdeclength_lambda_omegam(-99),
  fchi2_lambda_omegam(-99),
  fmass_lamdba_omegam(-99),
  fdcatoPVxy_lambda_omegam(-99),
  fdcatoPV_lambda_omegam(-99),
  fpt_lambda_omegam(-99),
//Omega                                                                                                                 
  fcpa_omega_omegam(-99),
  fpcaxy_omega_omegam(-99),
  fpca_omega_omegam(-99),
  fdeclength_omega_omegam(-99),
  fchi2_omega_omegam(-99),
  fmass_omega_omegam(-99),
  fpx_omega_omegam(-99),
  fpy_omega_omegam(-99),
  fpz_omega_omegam(-99),
  fpt_omega_omegam(-99),
  feta_omega_omegam(-99),
  fphi_omega_omegam(-99),
  fdcatoPVxy_omega_omegam(-99),
  fdcatoPV_omega_omegam(-99),
  fct_omega_omegam(-99),
  fpx_omega_omegam_mc(-99),
  fpy_omega_omegam_mc(-99),
  fpz_omega_omegam_mc(-99),
  fpt_omega_omegam_mc(-99),
//TTree_omegap
  fTree_omegap(0),
  fevtID_omegap(-99),
  fzvertex_omegap(-99),
  fcentrality_omegap(-99),
//Daughter track                                                                                                     
  fdca_daughter1_omegap(-99),
  fdca_daughter2_omegap(-99),
  fdca_daughter3_omegap(-99),
  ftrkID_daughter1_omegap(-99),
  ftrkID_daughter2_omegap(-99),
  ftrkID_daughter3_omegap(-99),
//Lambda                                                                                                             
  fcpa_lambda_omegap(-99),
  fpcaxy_lambda_omegap(-99),
  fpca_lambda_omegap(-99),
  fdeclength_lambda_omegap(-99),
  fchi2_lambda_omegap(-99),
  fmass_lamdba_omegap(-99),
  fdcatoPVxy_lambda_omegap(-99),
  fdcatoPV_lambda_omegap(-99),
  fpt_lambda_omegap(-99),
//Omega                                                                                                                 
  fcpa_omega_omegap(-99),
  fpcaxy_omega_omegap(-99),
  fpca_omega_omegap(-99),
  fdeclength_omega_omegap(-99),
  fchi2_omega_omegap(-99),
  fmass_omega_omegap(-99),
  fpx_omega_omegap(-99),
  fpy_omega_omegap(-99),
  fpz_omega_omegap(-99),
  fpt_omega_omegap(-99),
  feta_omega_omegap(-99),
  fphi_omega_omegap(-99),
  fdcatoPVxy_omega_omegap(-99),
  fdcatoPV_omega_omegap(-99),
  fct_omega_omegap(-99),
  fpx_omega_omegap_mc(-99),
  fpy_omega_omegap_mc(-99),
  fpz_omega_omegap_mc(-99),
  fpt_omega_omegap_mc(-99),
//TTree sidebandm   
  fTree_sidebandm(0),
  fevtID_sidebandm(-99),
  fzvertex_sidebandm(-99),
  fcentrality_sidebandm(-99),
//Daughter track                                                                                          
  fdca_daughter1_sidebandm(-99),
  fdca_daughter2_sidebandm(-99),
  fdca_daughter3_sidebandm(-99),
  ftrkID_daughter1_sidebandm(-99),
  ftrkID_daughter2_sidebandm(-99),
  ftrkID_daughter3_sidebandm(-99),
//Lambda                                                                                                   
  fcpa_lambda_sidebandm(-99),
  fpcaxy_lambda_sidebandm(-99),
  fpca_lambda_sidebandm(-99),
  fdeclength_lambda_sidebandm(-99),
  fchi2_lambda_sidebandm(-99),
  fmass_lamdba_sidebandm(-99),
  fdcatoPVxy_lambda_sidebandm(-99),
  fdcatoPV_lambda_sidebandm(-99),
  fpt_lambda_sidebandm(-99),
//Omega
  fcpa_omega_sidebandm(-99),
  fpcaxy_omega_sidebandm(-99),
  fpca_omega_sidebandm(-99),
  fdeclength_omega_sidebandm(-99),
  fchi2_omega_sidebandm(-99),
  fmass_omega_sidebandm(-99),
  fpx_omega_sidebandm(-99),
  fpy_omega_sidebandm(-99),
  fpz_omega_sidebandm(-99),
  fpt_omega_sidebandm(-99),
  feta_omega_sidebandm(-99),
  fphi_omega_sidebandm(-99),
  fdcatoPVxy_omega_sidebandm(-99),
  fdcatoPV_omega_sidebandm(-99),
  fct_omega_sidebandm(-99),
  fpx_omega_sidebandm_mc(-99),
  fpy_omega_sidebandm_mc(-99),
  fpz_omega_sidebandm_mc(-99),
  fpt_omega_sidebandm_mc(-99),
//TTree    
  fTree_sidebandp(0),
  fevtID_sidebandp(-99),
  fzvertex_sidebandp(-99),
  fcentrality_sidebandp(-99),
//Daughter track  
  fdca_daughter1_sidebandp(-99),
  fdca_daughter2_sidebandp(-99),
  fdca_daughter3_sidebandp(-99),
  ftrkID_daughter1_sidebandp(-99),
  ftrkID_daughter2_sidebandp(-99),
  ftrkID_daughter3_sidebandp(-99),
//Lambda
  fcpa_lambda_sidebandp(-99),
  fpcaxy_lambda_sidebandp(-99),
  fpca_lambda_sidebandp(-99),
  fdeclength_lambda_sidebandp(-99),
  fchi2_lambda_sidebandp(-99),
  fmass_lamdba_sidebandp(-99),
  fdcatoPVxy_lambda_sidebandp(-99),
  fdcatoPV_lambda_sidebandp(-99),
  fpt_lambda_sidebandp(-99),
//Omega 
  fcpa_omega_sidebandp(-99),
  fpcaxy_omega_sidebandp(-99),
  fpca_omega_sidebandp(-99),
  fdeclength_omega_sidebandp(-99),
  fchi2_omega_sidebandp(-99),
  fmass_omega_sidebandp(-99),
  fpx_omega_sidebandp(-99),
  fpy_omega_sidebandp(-99),
  fpz_omega_sidebandp(-99),
  fpt_omega_sidebandp(-99),
  feta_omega_sidebandp(-99),
  fphi_omega_sidebandp(-99),
  fdcatoPVxy_omega_sidebandp(-99),
  fdcatoPV_omega_sidebandp(-99),
  fct_omega_sidebandp(-99),
  fpx_omega_sidebandp_mc(-99),
  fpy_omega_sidebandp_mc(-99),
  fpz_omega_sidebandp_mc(-99),
  fpt_omega_sidebandp_mc(-99),
//TTree cascade    
  fTree_cascade(0),
  fpt_daughter1(0),
  fpt_daughter2(0),
  fpt_daughter3(0),
  fpt_daughter1_cascade(0),
  fpt_daughter2_cascade(0),
  fpt_daughter3_cascade(0),
  fcpa_lambda_cascade(0),
  fpca_lambda_cascade(0),
  fpt_lambda_cascade(0),
  fcpa_omega_cascade(0),
  fpca_omega_cascade(0),
  fmass_omega_cascade(0),
  fpt_omega_cascade(0),
//TTree_proton
  fTree_proton(0),
  fevtID_proton(-99),
  fzvertex_proton(-99),
  fcentrality_proton(-99),
  ftrkID_proton(-99),
  fdcaxy_proton(-99),
  fdcaz_proton(-99),
  fpx_proton(-99),
  fpy_proton(-99),
  fpz_proton(-99),
  fpt_proton(-99),
  fpTPC_proton(-99),
  fnSigmaTPC_proton(-99),
  fnSigmaTOF_proton(-99),
  fTrkMassTOF_proton(-99),
  fnSigmaITS_proton(-99),
  feta_proton(-99),
  fphi_proton(-99),
//TTree_antiproton
  fTree_antiproton(0),
  fevtID_antiproton(-99),
  fzvertex_antiproton(-99),
  fcentrality_antiproton(-99),
  ftrkID_antiproton(-99),
  fdcaxy_antiproton(-99),
  fdcaz_antiproton(-99),
  fpx_antiproton(-99),
  fpy_antiproton(-99),
  fpz_antiproton(-99),
  fpt_antiproton(-99),
  fpTPC_antiproton(-99),
  fnSigmaTPC_antiproton(-99),
  fnSigmaTOF_antiproton(-99),
  fTrkMassTOF_antiproton(-99),
  fnSigmaITS_antiproton(-99),
  feta_antiproton(-99),
  fphi_antiproton(-99),
  fOutputList(0)
{
 
}
 
//_____________________________________________________________________________
AliAnalysisTaskpOmegaDibaryon::AliAnalysisTaskpOmegaDibaryon(const char *name) : 
  AliAnalysisTaskSE(name),
  fEvent(0x0),
  fESD(0x0),
  fAOD(0),
  fHeader(0),
  fPIDResponse(0),
  ftrigBit(0),
  fEstimator("V0M"),
  fMultSelection(0x0),
  fCentralityMain(0.),
  fCentralityMin(0.),
  fCentralityMax(0.),
  fIsFlowTask(kFALSE),
  fFlowQnVectorMgr(0x0),
  //==================== variable
  nevt(0),
//========== Event selection
  isEvtMB(kFALSE),
  isPileupEvt(kFALSE),
  isGoodEvt(kFALSE),
  centralityV0M(-99),
  centralityCL0(-99),
  centralityCL1(-99),
  centralityV0A(-99),
  centralityV0C(-99),
  centralityZNA(-99),
  centralityZNC(-99),
  centralityMain(-99),
  //========== Event
  nTracks(0),
  nV0s(0),
  nCascades(0),
  vecTarget(),
  kMagF(0),
  PrimVertexKF(),
  vertex(),
  primVertex(),
  primKFVertex(),
  //========== Track 
  lorentzsum(),
  particleProton(),
  particlePion(),
  h(),
  np(0),
  nap(0),
  npi(0),
  nppi(0),
  track(),
  dcaxy(-99),
  dcaz(-99),
  dca(-99),
  dcaxyn(-99),
  dcazn(-99),
  dcan(-99),
  dDCA(),
  cDCA(),
  v0Vtx(),
  track_pT(),
  Daugproton(kFALSE),
  Daugpion(kFALSE),
  isProton(kFALSE),
  isAntiProton(kFALSE),
  isNPion(kFALSE),
  isPPion(kFALSE),
  isBachpion(kFALSE),
  isBachkaon(kFALSE),
  isLambdadecayProton(kFALSE),
  isLambdadecayAntiProton(kFALSE),
  isLambdadecayPPion(kFALSE),
  isLambdadecayNPion(kFALSE),
  isBachPPion(kFALSE),
  isBachNPion(kFALSE),
  isBachPKaon(kFALSE),
  isBachNKaon(kFALSE),
  pTCut_Proton(kFALSE),
  pTCut_Pion(kFALSE),
  pTCut_Kaon(kFALSE),
//========== Array 
  ProtonArray(),
  AntiProtonArray(),
  NPionArray(),
  PPionArray(),
  LambdadecayProtonArray(),
  LambdadecayAntiProtonArray(),
  LambdadecayPPionArray(),
  LambdadecayNPionArray(),
  XidecayBachPPionArray(),
  XidecayBachNPionArray(),
  XidecayBachPKaonArray(),
  XidecayBachNKaonArray(),
//========== TrackParam
  exProtonTrack(),
  exPionTrack(),
  exBachPionTrack(),
  exBachKaonTrack(),
//========== KF particle 
  ProtonTrack(),
  PionTrack(),
  KFProtonTrack(),
  KFPPionTrack(),
  KFPionTrack(),
  KFMother(),
  KFMother_K0s(),
  SecVertexKF(),
  SecVertexErrKF(),
  KFsecVertex(),
  fDCADaugter(-99),
  fDcaSecDaughterProton(-99),
  fDcaSecDaughterPion(-99),
  impar(),
  dd(),
  fPA(-99),
  DecLength(-99),
//========== V0
  v0(),
  ptrack(),
  ntrack(),
  dcaxyp_v0(-99),
  dcazp_v0(-99),
  dcap_v0(-99),
  dcaxyn_v0(-99),
  dcazn_v0(-99),
  dcan_v0(-99),
  dDCAp_v0(),
  dDCAn_v0(),
  cDCAp_v0(),
  cDCAn_v0(),
  v0Vtx_v0(),
  energy_v0(-99),
  fcpa_v0(-99),
  ftransradius_v0(-99),
  Daugproton_v0(kFALSE),
  Daugantiproton_v0(kFALSE),
  Daugpion_v0(kFALSE),
  Daugnpion_v0(kFALSE),
  lambda_v0(kFALSE),
  antilambda_v0(kFALSE),
//========== KF particle
  ndp(-99),
  ndap(-99),
  ndpp(-99),
  ndnp(-99),
  nbpp(-99),
  nbnp(-99),
  nbpk(-99),
  nbnk(-99),
//DCAtoPrimVtxXY(-99),
//DCAtoPrimVtxZ(-99),
  DCAtoPrimVtx(-99),
  BachNPionTrack(),
  KFSubMother(),
  KFBachPionTrack(),
  KFBachKaonTrack(),
  KFMother_Omega(),
  SecVertexKF_Xi(),
  SecVertexErrKF_Xi(),
  KFsecVertex_Xi(),
  SecVertexKF_Omega(),
  SecVertexErrKF_Omega(),
  KFsecVertex_Omega(),
  lorentzsum_Xi(),
  particleBachPion(),
  h_Xi(),
  dd_Xi(),
  lorentzsum_Omega(),
  particleBachKaon(),
  fPA_Xi(-99),
  DecLength_Xi(-99),
  DCAdauXY_Xi(-99),
  DCAdauZ_Xi(-99),
  DCAdau_Xi(-99),
  ftrkID_daughter1(-99),
  ftrkID_daughter2(-99),
  ftrkID_daughter3(-99),
//========== Cascade
  casc(),
  pTrackXi(),
  nTrackXi(),
  bachTrackXi(),
  dcaxyp(-99),
  dcazp(-99),
  dcap(-99),
  dcaxyb(-99),
  dcazb(-99),
  dcab(-99),
  dDCAp(),
  cDCAp(),
  dDCAn(),
  cDCAn(),
  dDCAb(),
  cDCAb(),
  v0Vtx_casc(),
  xivertex(),
  lambdacpa(-99),
  lambdatransradius(-99),
  xicpa(-99),
  xitransradius(-99),
  v0vtxx(-99),
  v0vtxy(-99),
  v0vtxz(-99),
  xivtxx(-99),
  xivtxy(-99),
  xivtxz(-99),
//TTree_omegam
  fTree_omegam(0),
  fevtID_omegam(-99),
  fzvertex_omegam(-99),
  fcentrality_omegam(-99),
//Daughter track                                                                                                     
  fdca_daughter1_omegam(-99),
  fdca_daughter2_omegam(-99),
  fdca_daughter3_omegam(-99),
  ftrkID_daughter1_omegam(-99),
  ftrkID_daughter2_omegam(-99),
  ftrkID_daughter3_omegam(-99),
//Lambda                                                                                                             
  fcpa_lambda_omegam(-99),
  fpcaxy_lambda_omegam(-99),
  fpca_lambda_omegam(-99),
  fdeclength_lambda_omegam(-99),
  fchi2_lambda_omegam(-99),
  fmass_lamdba_omegam(-99),
  fdcatoPVxy_lambda_omegam(-99),
  fdcatoPV_lambda_omegam(-99),
  fpt_lambda_omegam(-99),
//Omega                                                                                                                 
  fcpa_omega_omegam(-99),
  fpcaxy_omega_omegam(-99),
  fpca_omega_omegam(-99),
  fdeclength_omega_omegam(-99),
  fchi2_omega_omegam(-99),
  fmass_omega_omegam(-99),
  fpx_omega_omegam(-99),
  fpy_omega_omegam(-99),
  fpz_omega_omegam(-99),
  fpt_omega_omegam(-99),
  feta_omega_omegam(-99),
  fphi_omega_omegam(-99),
  fdcatoPVxy_omega_omegam(-99),
  fdcatoPV_omega_omegam(-99),
  fct_omega_omegam(-99),
  fpx_omega_omegam_mc(-99),
  fpy_omega_omegam_mc(-99),
  fpz_omega_omegam_mc(-99),
  fpt_omega_omegam_mc(-99),
//TTree_omegap
  fTree_omegap(0),
  fevtID_omegap(-99),
  fzvertex_omegap(-99),
  fcentrality_omegap(-99),
//Daughter track                                                                                                     
  fdca_daughter1_omegap(-99),
  fdca_daughter2_omegap(-99),
  fdca_daughter3_omegap(-99),
  ftrkID_daughter1_omegap(-99),
  ftrkID_daughter2_omegap(-99),
  ftrkID_daughter3_omegap(-99),
//Lambda                                                                                                             
  fcpa_lambda_omegap(-99),
  fpcaxy_lambda_omegap(-99),
  fpca_lambda_omegap(-99),
  fdeclength_lambda_omegap(-99),
  fchi2_lambda_omegap(-99),
  fmass_lamdba_omegap(-99),
  fdcatoPVxy_lambda_omegap(-99),
  fdcatoPV_lambda_omegap(-99),
  fpt_lambda_omegap(-99),
//Omega                                                                                                                 
  fcpa_omega_omegap(-99),
  fpcaxy_omega_omegap(-99),
  fpca_omega_omegap(-99),
  fdeclength_omega_omegap(-99),
  fchi2_omega_omegap(-99),
  fmass_omega_omegap(-99),
  fpx_omega_omegap(-99),
  fpy_omega_omegap(-99),
  fpz_omega_omegap(-99),
  fpt_omega_omegap(-99),
  feta_omega_omegap(-99),
  fphi_omega_omegap(-99),
  fdcatoPVxy_omega_omegap(-99),
  fdcatoPV_omega_omegap(-99),
  fct_omega_omegap(-99),
  fpx_omega_omegap_mc(-99),
  fpy_omega_omegap_mc(-99),
  fpz_omega_omegap_mc(-99),
  fpt_omega_omegap_mc(-99),
//TTree    
  fTree_sidebandm(0),
  fevtID_sidebandm(-99),
  fzvertex_sidebandm(-99),
  fcentrality_sidebandm(-99),
//Daughter track  
  fdca_daughter1_sidebandm(-99),
  fdca_daughter2_sidebandm(-99),
  fdca_daughter3_sidebandm(-99),
  ftrkID_daughter1_sidebandm(-99),
  ftrkID_daughter2_sidebandm(-99),
  ftrkID_daughter3_sidebandm(-99),
//Lambda
  fcpa_lambda_sidebandm(-99),
  fpcaxy_lambda_sidebandm(-99),
  fpca_lambda_sidebandm(-99),
  fdeclength_lambda_sidebandm(-99),
  fchi2_lambda_sidebandm(-99),
  fmass_lamdba_sidebandm(-99),
  fdcatoPVxy_lambda_sidebandm(-99),
  fdcatoPV_lambda_sidebandm(-99),
  fpt_lambda_sidebandm(-99),
//Omega 
  fcpa_omega_sidebandm(-99),
  fpcaxy_omega_sidebandm(-99),
  fpca_omega_sidebandm(-99),
  fdeclength_omega_sidebandm(-99),
  fchi2_omega_sidebandm(-99),
  fmass_omega_sidebandm(-99),
  fpx_omega_sidebandm(-99),
  fpy_omega_sidebandm(-99),
  fpz_omega_sidebandm(-99),
  fpt_omega_sidebandm(-99),
  feta_omega_sidebandm(-99),
  fphi_omega_sidebandm(-99),
  fdcatoPVxy_omega_sidebandm(-99),
  fdcatoPV_omega_sidebandm(-99),
  fct_omega_sidebandm(-99),
  fpx_omega_sidebandm_mc(-99),
  fpy_omega_sidebandm_mc(-99),
  fpz_omega_sidebandm_mc(-99),
  fpt_omega_sidebandm_mc(-99),
//TTree    
  fTree_sidebandp(0),
  fevtID_sidebandp(-99),
  fzvertex_sidebandp(-99),
  fcentrality_sidebandp(-99),
//Daughter track  
  fdca_daughter1_sidebandp(-99),
  fdca_daughter2_sidebandp(-99),
  fdca_daughter3_sidebandp(-99),
  ftrkID_daughter1_sidebandp(-99),
  ftrkID_daughter2_sidebandp(-99),
  ftrkID_daughter3_sidebandp(-99),
//Lambda
  fcpa_lambda_sidebandp(-99),
  fpcaxy_lambda_sidebandp(-99),
  fpca_lambda_sidebandp(-99),
  fdeclength_lambda_sidebandp(-99),
  fchi2_lambda_sidebandp(-99),
  fmass_lamdba_sidebandp(-99),
  fdcatoPVxy_lambda_sidebandp(-99),
  fdcatoPV_lambda_sidebandp(-99),
  fpt_lambda_sidebandp(-99),
//Omega 
  fcpa_omega_sidebandp(-99),
  fpcaxy_omega_sidebandp(-99),
  fpca_omega_sidebandp(-99),
  fdeclength_omega_sidebandp(-99),
  fchi2_omega_sidebandp(-99),
  fmass_omega_sidebandp(-99),
  fpx_omega_sidebandp(-99),
  fpy_omega_sidebandp(-99),
  fpz_omega_sidebandp(-99),
  fpt_omega_sidebandp(-99),
  feta_omega_sidebandp(-99),
  fphi_omega_sidebandp(-99),
  fdcatoPVxy_omega_sidebandp(-99),
  fdcatoPV_omega_sidebandp(-99),
  fct_omega_sidebandp(-99),
  fpx_omega_sidebandp_mc(-99),
  fpy_omega_sidebandp_mc(-99),
  fpz_omega_sidebandp_mc(-99),
  fpt_omega_sidebandp_mc(-99),
//TTree cascade    
  fTree_cascade(0),
  fpt_daughter1(0),
  fpt_daughter2(0),
  fpt_daughter3(0),
  fpt_daughter1_cascade(0),
  fpt_daughter2_cascade(0),
  fpt_daughter3_cascade(0),
  fcpa_lambda_cascade(0),
  fpca_lambda_cascade(0),
  fpt_lambda_cascade(0),
  fcpa_omega_cascade(0),
  fpca_omega_cascade(0),
  fmass_omega_cascade(0),
  fpt_omega_cascade(0),
//TTree_proton
  fTree_proton(0),
  fevtID_proton(-99),
  fzvertex_proton(-99),
  fcentrality_proton(-99),
  ftrkID_proton(-99),
  fdcaxy_proton(-99),
  fdcaz_proton(-99),
  fpx_proton(-99),
  fpy_proton(-99),
  fpz_proton(-99),
  fpt_proton(-99),
  fpTPC_proton(-99),
  fnSigmaTPC_proton(-99),
  fnSigmaTOF_proton(-99),
  fTrkMassTOF_proton(-99),
  fnSigmaITS_proton(-99),
  feta_proton(-99),
  fphi_proton(-99),
//TTree_antiproton
  fTree_antiproton(0),
  fevtID_antiproton(-99),
  fzvertex_antiproton(-99),
  fcentrality_antiproton(-99),
  ftrkID_antiproton(-99),
  fdcaxy_antiproton(-99),
  fdcaz_antiproton(-99),
  fpx_antiproton(-99),
  fpy_antiproton(-99),
  fpz_antiproton(-99),
  fpt_antiproton(-99),
  fpTPC_antiproton(-99),
  fnSigmaTPC_antiproton(-99),
  fnSigmaTOF_antiproton(-99),
  fTrkMassTOF_antiproton(-99),
  fnSigmaITS_antiproton(-99),
  feta_antiproton(-99),
  fphi_antiproton(-99),
  fOutputList(0)
{
  // Constructor
  DefineInput(0,TChain::Class());
  DefineOutput(1,THashList::Class());
  DefineOutput(2,TTree::Class());
  DefineOutput(3,TTree::Class());
  DefineOutput(4,TTree::Class());
  DefineOutput(5,TTree::Class());
  DefineOutput(6,TTree::Class());
  DefineOutput(7,TTree::Class());
  DefineOutput(8,TTree::Class());

}

//________________________________________________________________________
AliAnalysisTaskpOmegaDibaryon::~AliAnalysisTaskpOmegaDibaryon()
{
  // destructor
  if(fOutputList){
    delete fOutputList;
  }
  
}

//________________________________________________________________________
const Int_t AliAnalysisTaskpOmegaDibaryon::fgkPdgCode[] = {
  211,
  2212,
  3122,
  310,
  3312,
  3334
};

//________________________________________________________________________
void AliAnalysisTaskpOmegaDibaryon::UserCreateOutputObjects()
{

  //===============================
  fOutputList=new THashList();
  fOutputList->SetOwner(kTRUE);

  //========== Event cut ==========
  fOutputList->Add(new TH1F("fEvtCounter","",1,0,1));
  fOutputList->Add(new TH1F("fEvtPassCut","",1,0,1));
  fOutputList->Add(new TH1F("fEvtVtxX","",400,-20.,20.));
  fOutputList->Add(new TH1F("fEvtVtxY","",400,-20.,20.));
  fOutputList->Add(new TH1F("fEvtVtxZ","",400,-20.,20.));
  fOutputList->Add(new TH1F("fEvtVtxTrk","",1000,0.,1000.));
  fOutputList->Add(new TH1F("fSPDVtxZ","",400,-20.,20.));
  fOutputList->Add(new TH1F("fSPDVtxTrk","",1000,0.,1000.));
  fOutputList->Add(new TH2F("fSPDVtxCor","",400,-20.,20.,400,-20.,20.));
  fOutputList->Add(new TH1F("fSPDVtxDisp","",300,0.,1.5));

  fOutputList->Add(new TH1F("hCentralityV0M","centrality V0M",101,0,101));
  fOutputList->Add(new TH1F("hCentralityCL0","centrality CL0",101,0,101));
  fOutputList->Add(new TH1F("hCentralityCL1","centrality CL1",101,0,101));
  fOutputList->Add(new TH1F("hCentralityZNA","centrality ZNA",101,0,101));
  fOutputList->Add(new TH1F("hCentralityZNC","centrality ZNC",101,0,101));
  fOutputList->Add(new TH1F("hCentralityV0A","centrality V0A",101,0,101));
  fOutputList->Add(new TH1F("hCentralityV0C","centrality V0C",101,0,101));
  fOutputList->Add(new TH1F("hCentralityMain","centrality Main",101,0,101));

  fOutputList->Add(new TH1F("hCentralityV0M_cent","centrality V0M MostCent",101,0,101));
  fOutputList->Add(new TH1F("hCentralityV0M_semicent","centrality V0M SemiCent",101,0,101));
  fOutputList->Add(new TH1F("hCentralityV0M_mb","centrality V0M MB",101,0,101));

  fOutputList->Add(new TH1F("hCentralityV0M_cut","centrality V0M cut",101,0,101));
  //========== KF vertex ==========
  fOutputList->Add(new TH1F("fZvertex_aod","",100,-20,20));
  //using KF
  fOutputList->Add(new TH1F("fTrkDCA_KF","",1000,0,10));
  fOutputList->Add(new TH1F("fSecvertex_KF","",3000,0,150));
  fOutputList->Add(new TH1F("fDCADaug_KF","",1000,0,50));
  fOutputList->Add(new TH1F("fCPA_KF","",1000,0,1));
  fOutputList->Add(new TH1F("fTransradius_KF","",2000,0,200));
  fOutputList->Add(new TH1F("fChi2_KF","",1000,0,100));
  fOutputList->Add(new TH1F("fInvMassLambda_KF","",1000,1,2));
  //mass window
  fOutputList->Add(new TH2F("fpTvsMassLambda_KF_mass","",1000,0,10,1000,1,2));
  fOutputList->Add(new TH1F("fSecvertex_KF_mass","",3000,0,150));
  fOutputList->Add(new TH1F("fCPA_KF_mass","",10000,0,1));
  fOutputList->Add(new TH1F("fTransradius_KF_mass","",2000,0,200));
  fOutputList->Add(new TH1F("fDCADaug_KF_mass","",1000,0,5));
  fOutputList->Add(new TH1F("fChi2_KF_mass","",1000,0,100));
  //side
  fOutputList->Add(new TH2F("fpTvsMassLambda_KF_side","",1000,0,10,1000,1,2));
  fOutputList->Add(new TH1F("fSecvertex_KF_side","",3000,0,150));
  fOutputList->Add(new TH1F("fCPA_KF_side","",10000,0,1));
  fOutputList->Add(new TH1F("fTransradius_KF_side","",2000,0,200));
  fOutputList->Add(new TH1F("fDCADaug_KF_side","",1000,0,5));
  fOutputList->Add(new TH1F("fChi2_KF_side","",1000,0,100));
  //using V0
  fOutputList->Add(new TH1F("fTrkDCA_V0","",1000,0,10));
  fOutputList->Add(new TH1F("fSecvertex_V0","",100,-20,20));
  fOutputList->Add(new TH1F("fDCADaug_V0","",500,0,50));
  fOutputList->Add(new TH1F("fCPA_V0","",1000,0,1));
  fOutputList->Add(new TH1F("fTransradius_V0","",2000,0,200));
  fOutputList->Add(new TH1F("fInvMassLambda_V0","",1000,1,2));
  fOutputList->Add(new TH2F("fpTvsMassLambda_V0","",1000,0,10,1000,1,2));
  fOutputList->Add(new TH1F("fpTLambda_V0","",1000,0,10));
  //KF particle cascade V0
  fOutputList->Add(new TH1F("fpT_CascDecayProton_KF","",1000,0,10));
  fOutputList->Add(new TH1F("fpT_CascDecayPion_KF","",1000,0,10));
  fOutputList->Add(new TH1F("fSecvertex_CascDecay_KF","",3000,0,300));
  fOutputList->Add(new TH1F("fDCADaug_CascDecay_KF","",1000,0,10));
  fOutputList->Add(new TH1F("fDCAtoPrimVtx_CascDecay_KF","",1000,0,10));
  fOutputList->Add(new TH1F("fCPA_CascDecay_KF","",1000,0,1));
  fOutputList->Add(new TH1F("fChi2_CascDecay_KF","",1000,0,100));
  fOutputList->Add(new TH1F("fDecaylength_CascDecay_KF","",2000,0,200));
  fOutputList->Add(new TH1F("fInvMassLambda_CascDecay_KF","",1000,1,2));
  //KF particle cascade Xi 
  fOutputList->Add(new TH1F("fSecvertex_Casc_KF","",3000,0,300));
  fOutputList->Add(new TH1F("fCPA_Casc_KF","",1000,0,1));
  fOutputList->Add(new TH1F("fDCAdau_Casc_KF","",10000,0,10));
  fOutputList->Add(new TH1F("fTransradius_Casc_KF","",2000,0,200));
  fOutputList->Add(new TH1F("fInvMassXi_Casc_KF","",1000,1,2));
  fOutputList->Add(new TH1F("fpT_CascBach_Pion_KF","",1000,0,10));
  fOutputList->Add(new TH1F("fpT_CascDecay_Lambda_KF","",1000,0,10)); 
  fOutputList->Add(new TH1F("fpT_Casc_Xi_KF","",1000,0,10));
  fOutputList->Add(new TH1F("fpT_Casc_AntiXi_KF","",1000,0,10));
  fOutputList->Add(new TH1F("fpT_CascBach_PPion_KF","",1000,0,10));
  fOutputList->Add(new TH1F("fpT_CascDecay_AntiLambda_KF","",1000,0,10)); 
  //using Cascade
  //cascade daughter select contents 
  fOutputList->Add(new TH1F("fcascdauEta","",100,-2,2));
  fOutputList->Add(new TH1F("fcascdauEtab","",100,-2,2));
  fOutputList->Add(new TH1F("fcascdauTPC","",200,0,200));
  fOutputList->Add(new TH1F("fcascdauTPCb","",200,0,200));
  fOutputList->Add(new TH1F("fcascdauPt","",1000,0,10));
  fOutputList->Add(new TH1F("fcascdauPtb","",1000,0,10));
  fOutputList->Add(new TH1F("fcascdaupT_Proton","",1000,0,10));
  fOutputList->Add(new TH1F("fcascdaupT_Pion","",1000,0,10));
  fOutputList->Add(new TH1F("fcascdaudcab","",100,0,1));
  fOutputList->Add(new TH1F("fcascdauPIDproton","",60,0,6));
  fOutputList->Add(new TH1F("fcascdauPIDpion","",60,0,6));
  fOutputList->Add(new TH1F("fcascbachPIDpion","",60,0,6));
  //cascade V0 select contents 
  fOutputList->Add(new TH1F("fcascV0Transradius","",3000,0,300));
  fOutputList->Add(new TH1F("fcascV0daugdca","",1000,0,10));
  fOutputList->Add(new TH1F("fcascV0dca","",1000,0,10));
  fOutputList->Add(new TH1F("fcascV0cpa","",1000,0,1));
  fOutputList->Add(new TH1F("hInvMassXidecayLambda","",10000,0,10));
  //Xi select contents  
  fOutputList->Add(new TH1F("XiTranverseradius","",2000,0,200));
  fOutputList->Add(new TH1F("dcaXiDaughters","",1000,0,10));
  fOutputList->Add(new TH1F("XiCPA","",1000,0,1));
  //fOutputList->Add(new TH1F("hInvMassOmega","",10000,0,10));
  //invariant mass xi 
  fOutputList->Add(new TH1F("hInvMassXi","",1000,1,2));
  fOutputList->Add(new TH1F("fPtcascade_xi","",1000,0,10));
  fOutputList->Add(new TH1F("fPtcascade_antixi","",1000,0,10));
  fOutputList->Add(new TH1F("fPtcascdecay_lambda","",1000,0,10));
  fOutputList->Add(new TH1F("fPtcascdecay_antilambda","",1000,0,10));
  fOutputList->Add(new TH1F("fPtBachelor","",1000,0,10));
  fOutputList->Add(new TH1F("fPhiXidecayLambda","",1000,0,10));
  fOutputList->Add(new TH1F("fPhi-Psi_xi","",2000,-10,10));
  fOutputList->Add(new TH2F("fTrkMasspT_Proton","",1000,0,10,10000,0,10));
  fOutputList->Add(new TH1F("fSigma_Proton","",100,0,100));
  fOutputList->Add(new TH1F("fSigma_AntiProton","",100,0,100));
  fOutputList->Add(new TH2F("fPIDITS_proton","",1000,0,10,5000,0,500));
  fOutputList->Add(new TH1F("fSigmaITS_Proton","",200,-100,100));
  fOutputList->Add(new TH1F("fSigmaITS_Proton_cut","",200,-100,100));

  //================= Eta vs Phi
  fOutputList->Add(new TH1F("fEta_omega","",200,-10,10));
  fOutputList->Add(new TH1F("fPhi_omega","",200,-10,10));
  fOutputList->Add(new TH1F("fEta_proton","",200,-10,10));
  fOutputList->Add(new TH1F("fPhi_proton","",200,-10,10));

  //================= TTree_omegam
  fTree_omegam = new TTree("fTree_omegam", "fTree_omegam");
  fTree_omegam->Branch("fevtID", &fevtID_omegam, "fevtID/L");
  fTree_omegam->Branch("fzvertex", &fzvertex_omegam, "fzvertex/D");
  fTree_omegam->Branch("fcentrality", &fcentrality_omegam, "fcentrality/D");

  //Daughter track  
  fTree_omegam->Branch("fdca_daughter1", &fdca_daughter1_omegam, "fdca_daughter1/F");
  fTree_omegam->Branch("fdca_daughter2", &fdca_daughter2_omegam, "fdca_daughter2/F");
  fTree_omegam->Branch("fdca_daughter3", &fdca_daughter3_omegam, "fdca_daughter3/F");
  fTree_omegam->Branch("ftrkID_daughter1", &ftrkID_daughter1_omegam, "ftrkID_daughter1/I");
  fTree_omegam->Branch("ftrkID_daughter2", &ftrkID_daughter2_omegam, "ftrkID_daughter2/I");
  fTree_omegam->Branch("ftrkID_daughter3", &ftrkID_daughter3_omegam, "ftrkID_daughter3/I");

  //Lambda  
  fTree_omegam->Branch("fcpa_lambda", &fcpa_lambda_omegam, "fcpa_lambda/D");
  fTree_omegam->Branch("fpcaxy_lambda", &fpcaxy_lambda_omegam, "fpcaxy_lambda/D");
  fTree_omegam->Branch("fpca_lambda", &fpca_lambda_omegam, "fpca_lambda/D");
  fTree_omegam->Branch("fdeclength_lambda", &fdeclength_lambda_omegam, "fdeclength_lambda/D");
  fTree_omegam->Branch("fchi2_lambda", &fchi2_lambda_omegam, "fchi2_lambda/F");
  fTree_omegam->Branch("fmass_lamdba", &fmass_lamdba_omegam, "fmass_lamdba/D");
  fTree_omegam->Branch("fdcatoPVxy_lambda", &fdcatoPVxy_lambda_omegam, "fdcatoPVxy_lambda/D");
  fTree_omegam->Branch("fdcatoPV_lambda", &fdcatoPV_lambda_omegam, "fdcatoPV_lambda/D");
  fTree_omegam->Branch("fpt_lambda", &fpt_lambda_omegam, "fpt_lambda/D");

  //Omega  
  fTree_omegam->Branch("fcpa_omega", &fcpa_omega_omegam, "fcpa_omega/D");
  fTree_omegam->Branch("fpcaxy_omega", &fpcaxy_omega_omegam, "fpcaxy_omega/D");
  fTree_omegam->Branch("fpca_omega", &fpca_omega_omegam, "fpca_omega/D");
  fTree_omegam->Branch("fdeclength_omega", &fdeclength_omega_omegam, "fdeclength_omega/D");
  fTree_omegam->Branch("fchi2_omega", &fchi2_omega_omegam, "fchi2_omega/F");
  fTree_omegam->Branch("fmass_omega", &fmass_omega_omegam, "fmass_omega/D");
  fTree_omegam->Branch("fpx_omega", &fpx_omega_omegam, "fpx_omega/D");
  fTree_omegam->Branch("fpy_omega", &fpy_omega_omegam, "fpy_omega/D");
  fTree_omegam->Branch("fpz_omega", &fpz_omega_omegam, "fpz_omega/D");
  fTree_omegam->Branch("fpt_omega", &fpt_omega_omegam, "fpt_omega/D");
  fTree_omegam->Branch("feta_omega", &feta_omega_omegam, "feta_omega/D");
  fTree_omegam->Branch("fphi_omega", &fphi_omega_omegam, "fphi_omega/D");
  fTree_omegam->Branch("fdcatoPVxy_omega", &fdcatoPVxy_omega_omegam, "fdcatoPVxy_omega/D");
  fTree_omegam->Branch("fdcatoPV_omega", &fdcatoPV_omega_omegam, "fdcatoPV_omega/D");
  fTree_omegam->Branch("fct_omega", &fct_omega_omegam, "fct_omega/D");
  fTree_omegam->Branch("fpx_omega_mc", &fpx_omega_omegam_mc, "fpx_omega_mc/D");
  fTree_omegam->Branch("fpy_omega_mc", &fpy_omega_omegam_mc, "fpy_omega_mc/D"); 
  fTree_omegam->Branch("fpz_omega_mc", &fpz_omega_omegam_mc, "fpz_omega_mc/D");
  fTree_omegam->Branch("fpt_omega_mc", &fpt_omega_omegam_mc, "fpt_omega_mc/D");

  //================= TTree_omegap
  fTree_omegap = new TTree("fTree_omegap", "fTree_omegap");
  fTree_omegap->Branch("fevtID", &fevtID_omegap, "fevtID/L");
  fTree_omegap->Branch("fzvertex", &fzvertex_omegap, "fzvertex/D");
  fTree_omegap->Branch("fcentrality", &fcentrality_omegap, "fcentrality/D");

  //Daughter track  
  fTree_omegap->Branch("fdca_daughter1", &fdca_daughter1_omegap, "fdca_daughter1/F");
  fTree_omegap->Branch("fdca_daughter2", &fdca_daughter2_omegap, "fdca_daughter2/F");
  fTree_omegap->Branch("fdca_daughter3", &fdca_daughter3_omegap, "fdca_daughter3/F");
  fTree_omegap->Branch("ftrkID_daughter1", &ftrkID_daughter1_omegap, "ftrkID_daughter1/I");
  fTree_omegap->Branch("ftrkID_daughter2", &ftrkID_daughter2_omegap, "ftrkID_daughter2/I");
  fTree_omegap->Branch("ftrkID_daughter3", &ftrkID_daughter3_omegap, "ftrkID_daughter3/I");

  //Lambda  
  fTree_omegap->Branch("fcpa_lambda", &fcpa_lambda_omegap, "fcpa_lambda/D");
  fTree_omegap->Branch("fpcaxy_lambda", &fpcaxy_lambda_omegap, "fpcaxy_lambda/D");
  fTree_omegap->Branch("fpca_lambda", &fpca_lambda_omegap, "fpca_lambda/D");
  fTree_omegap->Branch("fdeclength_lambda", &fdeclength_lambda_omegap, "fdeclength_lambda/D");
  fTree_omegap->Branch("fchi2_lambda", &fchi2_lambda_omegap, "fchi2_lambda/F");
  fTree_omegap->Branch("fmass_lamdba", &fmass_lamdba_omegap, "fmass_lamdba/D");
  fTree_omegap->Branch("fdcatoPVxy_lambda", &fdcatoPVxy_lambda_omegap, "fdcatoPVxy_lambda/D");
  fTree_omegap->Branch("fdcatoPV_lambda", &fdcatoPV_lambda_omegap, "fdcatoPV_lambda/D");
  fTree_omegap->Branch("fpt_lambda", &fpt_lambda_omegap, "fpt_lambda/D");

  //Omega  
  fTree_omegap->Branch("fcpa_omega", &fcpa_omega_omegap, "fcpa_omega/D");
  fTree_omegap->Branch("fpcaxy_omega", &fpcaxy_omega_omegap, "fpcaxy_omega/D");
  fTree_omegap->Branch("fpca_omega", &fpca_omega_omegap, "fpca_omega/D");
  fTree_omegap->Branch("fdeclength_omega", &fdeclength_omega_omegap, "fdeclength_omega/D");
  fTree_omegap->Branch("fchi2_omega", &fchi2_omega_omegap, "fchi2_omega/F");
  fTree_omegap->Branch("fmass_omega", &fmass_omega_omegap, "fmass_omega/D");
  fTree_omegap->Branch("fpx_omega", &fpx_omega_omegap, "fpx_omega/D");
  fTree_omegap->Branch("fpy_omega", &fpy_omega_omegap, "fpy_omega/D");
  fTree_omegap->Branch("fpz_omega", &fpz_omega_omegap, "fpz_omega/D");
  fTree_omegap->Branch("fpt_omega", &fpt_omega_omegap, "fpt_omega/D");
  fTree_omegap->Branch("feta_omega", &feta_omega_omegap, "feta_omega/D");
  fTree_omegap->Branch("fphi_omega", &fphi_omega_omegap, "fphi_omega/D");
  fTree_omegap->Branch("fdcatoPVxy_omega", &fdcatoPVxy_omega_omegap, "fdcatoPVxy_omega/D");
  fTree_omegap->Branch("fdcatoPV_omega", &fdcatoPV_omega_omegap, "fdcatoPV_omega/D");
  fTree_omegap->Branch("fct_omega", &fct_omega_omegap, "fct_omega/D");
  fTree_omegap->Branch("fpx_omega_mc", &fpx_omega_omegap_mc, "fpx_omega_mc/D");
  fTree_omegap->Branch("fpy_omega_mc", &fpy_omega_omegap_mc, "fpy_omega_mc/D");
  fTree_omegap->Branch("fpz_omega_mc", &fpz_omega_omegap_mc, "fpz_omega_mc/D"); 
  fTree_omegap->Branch("fpt_omega_mc", &fpt_omega_omegap_mc, "fpt_omega_mc/D");

  //================= TTree_sidebandm
  fTree_sidebandm = new TTree("fTree_sidebandm", "fTree_sidebandm");
  fTree_sidebandm->Branch("fevtID", &fevtID_sidebandm, "fevtID/L");
  fTree_sidebandm->Branch("fzvertex", &fzvertex_sidebandm, "fzvertex/D");
  fTree_sidebandm->Branch("fcentrality", &fcentrality_sidebandm, "fcentrality/D");

  //Daughter track 
  fTree_sidebandm->Branch("fdca_daughter1", &fdca_daughter1_sidebandm, "fdca_daughter1/F");
  fTree_sidebandm->Branch("fdca_daughter2", &fdca_daughter2_sidebandm, "fdca_daughter2/F");
  fTree_sidebandm->Branch("fdca_daughter3", &fdca_daughter3_sidebandm, "fdca_daughter3/F");
  fTree_sidebandm->Branch("ftrkID_daughter1", &ftrkID_daughter1_sidebandm, "ftrkID_daughter1/I");
  fTree_sidebandm->Branch("ftrkID_daughter2", &ftrkID_daughter2_sidebandm, "ftrkID_daughter2/I");
  fTree_sidebandm->Branch("ftrkID_daughter3", &ftrkID_daughter3_sidebandm, "ftrkID_daughter3/I");

  //Lambda  
  fTree_sidebandm->Branch("fcpa_lambda", &fcpa_lambda_sidebandm, "fcpa_lambda/D");
  fTree_sidebandm->Branch("fpcaxy_lambda", &fpcaxy_lambda_sidebandm, "fpcaxy_lambda/D");
  fTree_sidebandm->Branch("fpca_lambda", &fpca_lambda_sidebandm, "fpca_lambda/D");
  fTree_sidebandm->Branch("fdeclength_lambda", &fdeclength_lambda_sidebandm, "fdeclength_lambda/D");
  fTree_sidebandm->Branch("fchi2_lambda", &fchi2_lambda_sidebandm, "fchi2_lambda/F");
  fTree_sidebandm->Branch("fmass_lamdba", &fmass_lamdba_sidebandm, "fmass_lamdba/D");
  fTree_sidebandm->Branch("fdcatoPVxy_lambda", &fdcatoPVxy_lambda_sidebandm, "fdcatoPVxy_lambda/D");
  fTree_sidebandm->Branch("fdcatoPV_lambda", &fdcatoPV_lambda_sidebandm, "fdcatoPV_lambda/D");
  fTree_sidebandm->Branch("fpt_lambda", &fpt_lambda_sidebandm, "fpt_lambda/D");

  //Omega  
  fTree_sidebandm->Branch("fcpa_omega", &fcpa_omega_sidebandm, "fcpa_omega/D");
  fTree_sidebandm->Branch("fpcaxy_omega", &fpcaxy_omega_sidebandm, "fpcaxy_omega/D");
  fTree_sidebandm->Branch("fpca_omega", &fpca_omega_sidebandm, "fpca_omega/D");
  fTree_sidebandm->Branch("fdeclength_omega", &fdeclength_omega_sidebandm, "fdeclength_omega/D");
  fTree_sidebandm->Branch("fchi2_omega", &fchi2_omega_sidebandm, "fchi2_omega/F");
  fTree_sidebandm->Branch("fmass_omega", &fmass_omega_sidebandm, "fmass_omega/D");  
  fTree_sidebandm->Branch("fpx_omega", &fpx_omega_sidebandm, "fpx_omega/D");  
  fTree_sidebandm->Branch("fpy_omega", &fpy_omega_sidebandm, "fpy_omega/D");  
  fTree_sidebandm->Branch("fpz_omega", &fpz_omega_sidebandm, "fpz_omega/D");  
  fTree_sidebandm->Branch("fpt_omega", &fpt_omega_sidebandm, "fpt_omega/D");  
  fTree_sidebandm->Branch("feta_omega", &feta_omega_sidebandm, "feta_omega/D");  
  fTree_sidebandm->Branch("fphi_omega", &fphi_omega_sidebandm, "fphi_omega/D");  
  fTree_sidebandm->Branch("fdcatoPVxy_omega", &fdcatoPVxy_omega_sidebandm, "fdcatoPVxy_omega/D");
  fTree_sidebandm->Branch("fdcatoPV_omega", &fdcatoPV_omega_sidebandm, "fdcatoPV_omega/D");
  fTree_sidebandm->Branch("fct_omega", &fct_omega_sidebandm, "fct_omega/D");
  fTree_sidebandm->Branch("fpx_omega_mc", &fpx_omega_sidebandm_mc, "fpx_omega_mc/D");
  fTree_sidebandm->Branch("fpy_omega_mc", &fpy_omega_sidebandm_mc, "fpy_omega_mc/D");
  fTree_sidebandm->Branch("fpz_omega_mc", &fpz_omega_sidebandm_mc, "fpz_omega_mc/D");
  fTree_sidebandm->Branch("fpt_omega_mc", &fpt_omega_sidebandm_mc, "fpt_omega_mc/D");

  //================= TTree_sidebandp
  fTree_sidebandp = new TTree("fTree_sidebandp", "fTree_sidebandp");
  fTree_sidebandp->Branch("fevtID", &fevtID_sidebandp, "fevtID/L");
  fTree_sidebandp->Branch("fzvertex", &fzvertex_sidebandp, "fzvertex/D");
  fTree_sidebandp->Branch("fcentrality", &fcentrality_sidebandp, "fcentrality/D");

  //Daughter track 
  fTree_sidebandp->Branch("fdca_daughter1", &fdca_daughter1_sidebandp, "fdca_daughter1/F");
  fTree_sidebandp->Branch("fdca_daughter2", &fdca_daughter2_sidebandp, "fdca_daughter2/F");
  fTree_sidebandp->Branch("fdca_daughter3", &fdca_daughter3_sidebandp, "fdca_daughter3/F");
  fTree_sidebandp->Branch("ftrkID_daughter1", &ftrkID_daughter1_sidebandp, "ftrkID_daughter1/I");
  fTree_sidebandp->Branch("ftrkID_daughter2", &ftrkID_daughter2_sidebandp, "ftrkID_daughter2/I");
  fTree_sidebandp->Branch("ftrkID_daughter3", &ftrkID_daughter3_sidebandp, "ftrkID_daughter3/I");

  //Lambda  
  fTree_sidebandp->Branch("fcpa_lambda", &fcpa_lambda_sidebandp, "fcpa_lambda/D");
  fTree_sidebandp->Branch("fpcaxy_lambda", &fpcaxy_lambda_sidebandp, "fpcaxy_lambda/D");
  fTree_sidebandp->Branch("fpca_lambda", &fpca_lambda_sidebandp, "fpca_lambda/D");
  fTree_sidebandp->Branch("fdeclength_lambda", &fdeclength_lambda_sidebandp, "fdeclength_lambda/D");
  fTree_sidebandp->Branch("fchi2_lambda", &fchi2_lambda_sidebandp, "fchi2_lambda/F");
  fTree_sidebandp->Branch("fmass_lamdba", &fmass_lamdba_sidebandp, "fmass_lamdba/D");
  fTree_sidebandp->Branch("fdcatoPVxy_lambda", &fdcatoPVxy_lambda_sidebandp, "fdcatoPVxy_lambda/D");
  fTree_sidebandp->Branch("fdcatoPV_lambda", &fdcatoPV_lambda_sidebandp, "fdcatoPV_lambda/D");
  fTree_sidebandp->Branch("fpt_lambda", &fpt_lambda_sidebandp, "fpt_lambda/D");

  //Omega  
  fTree_sidebandp->Branch("fcpa_omega", &fcpa_omega_sidebandp, "fcpa_omega/D");
  fTree_sidebandp->Branch("fpcaxy_omega", &fpcaxy_omega_sidebandp, "fpcaxy_omega/D");
  fTree_sidebandp->Branch("fpca_omega", &fpca_omega_sidebandp, "fpca_omega/D");
  fTree_sidebandp->Branch("fdeclength_omega", &fdeclength_omega_sidebandp, "fdeclength_omega/D");
  fTree_sidebandp->Branch("fchi2_omega", &fchi2_omega_sidebandp, "fchi2_omega/F");
  fTree_sidebandp->Branch("fmass_omega", &fmass_omega_sidebandp, "fmass_omega/D");  
  fTree_sidebandp->Branch("fpx_omega", &fpx_omega_sidebandp, "fpx_omega/D");  
  fTree_sidebandp->Branch("fpy_omega", &fpy_omega_sidebandp, "fpy_omega/D");  
  fTree_sidebandp->Branch("fpz_omega", &fpz_omega_sidebandp, "fpz_omega/D");  
  fTree_sidebandp->Branch("fpt_omega", &fpt_omega_sidebandp, "fpt_omega/D");  
  fTree_sidebandp->Branch("feta_omega", &feta_omega_sidebandp, "feta_omega/D");  
  fTree_sidebandp->Branch("fphi_omega", &fphi_omega_sidebandp, "fphi_omega/D");
  fTree_sidebandp->Branch("fdcatoPVxy_omega", &fdcatoPVxy_omega_sidebandp, "fdcatoPVxy_omega/D");    
  fTree_sidebandp->Branch("fdcatoPV_omega", &fdcatoPV_omega_sidebandp, "fdcatoPV_omega/D");    
  fTree_sidebandp->Branch("fct_omega", &fct_omega_sidebandp, "fct_omega/D");
  fTree_sidebandp->Branch("fpx_omega_mc", &fpx_omega_sidebandp_mc, "fpx_omega_mc/D");
  fTree_sidebandp->Branch("fpy_omega_mc", &fpy_omega_sidebandp_mc, "fpy_omega_mc/D");
  fTree_sidebandp->Branch("fpz_omega_mc", &fpz_omega_sidebandp_mc, "fpz_omega_mc/D"); 
  fTree_sidebandp->Branch("fpt_omega_mc", &fpt_omega_sidebandp_mc, "fpt_omega_mc/D");

  //================= TTree_cascade
  fTree_cascade = new TTree("fTree_cascade", "fTree_cascade");

  fTree_cascade->Branch("fpt_daughter1", &fpt_daughter1_cascade, "fpt_daughter1/D");
  fTree_cascade->Branch("fpt_daughter2", &fpt_daughter2_cascade, "fpt_daughter2/D");
  fTree_cascade->Branch("fpt_daughter3", &fpt_daughter3_cascade, "fpt_daughter3/D");

  //Lambda
  fTree_cascade->Branch("fcpa_lambda", &fcpa_lambda_cascade, "fcpa_lambda/D");
  fTree_cascade->Branch("fpca_lambda", &fpca_lambda_cascade, "fpca_lambda/D");
  //Omega
  fTree_cascade->Branch("fcpa_omega", &fcpa_omega_cascade, "fcpa_omega/D");
  fTree_cascade->Branch("fpca_omega", &fpca_omega_cascade, "fpca_omega/D");
  fTree_cascade->Branch("fmass_omega", &fmass_omega_cascade, "fmass_omega/D");
  fTree_cascade->Branch("fpt_omega", &fpt_omega_cascade, "fpt_omega/D");

  //================= TTree_proton
  fTree_proton = new TTree("fTree_proton", "fTree_proton");
  fTree_proton->Branch("fevtID", &fevtID_proton, "fevtID/L");
  fTree_proton->Branch("fzvertex", &fzvertex_proton, "fzvertex/D");
  fTree_proton->Branch("fcentrality", &fcentrality_proton, "fcentrality/D");

  fTree_proton->Branch("ftrkID_proton", &ftrkID_proton, "ftrkID_proton/I");
  fTree_proton->Branch("fdcaxy_proton", &fdcaxy_proton, "fdcaxy_proton/D");
  fTree_proton->Branch("fdcaz_proton", &fdcaz_proton, "fdcaz_proton/D");
  fTree_proton->Branch("fpx_proton", &fpx_proton, "fpx_proton/D");
  fTree_proton->Branch("fpy_proton", &fpy_proton, "fpy_proton/D");
  fTree_proton->Branch("fpz_proton", &fpz_proton, "fpz_proton/D");
  fTree_proton->Branch("fpt_proton", &fpt_proton, "fpt_proton/D");
  fTree_proton->Branch("fpTPC_proton", &fpTPC_proton, "fpTPC_proton/D");
  fTree_proton->Branch("fnSigmaTPC_proton", &fnSigmaTPC_proton, "fnSigmaTPC_proton/D");
  fTree_proton->Branch("fnSigmaTOF_proton", &fnSigmaTOF_proton, "fnSigmaTOF_proton/D");
  fTree_proton->Branch("fTrkMassTOF_proton", &fTrkMassTOF_proton, "fTrkMassTOF_proton/D");
  fTree_proton->Branch("fnSigmaITS_proton", &fnSigmaITS_proton, "fnSigmaITS_proton/D");
  fTree_proton->Branch("feta_proton", &feta_proton, "feta_proton/D");
  fTree_proton->Branch("fphi_proton", &fphi_proton, "fphi_proton/D");

  //================= TTree_antiproton
  fTree_antiproton = new TTree("fTree_antiproton", "fTree_antiproton");
  fTree_antiproton->Branch("fevtID", &fevtID_antiproton, "fevtID/L");
  fTree_antiproton->Branch("fzvertex", &fzvertex_antiproton, "fzvertex/D");
  fTree_antiproton->Branch("fcentrality", &fcentrality_antiproton, "fcentrality/D");

  fTree_antiproton->Branch("ftrkID_antiproton", &ftrkID_antiproton, "ftrkID_antiproton/I");
  fTree_antiproton->Branch("fdcaxy_antiproton", &fdcaxy_antiproton, "fdcaxy_antiproton/D");
  fTree_antiproton->Branch("fdcaz_antiproton", &fdcaz_antiproton, "fdcaz_antiproton/D");
  fTree_antiproton->Branch("fpx_antiproton", &fpx_antiproton, "fpx_antiproton/D");
  fTree_antiproton->Branch("fpy_antiproton", &fpy_antiproton, "fpy_antiproton/D");
  fTree_antiproton->Branch("fpz_antiproton", &fpz_antiproton, "fpz_antiproton/D");
  fTree_antiproton->Branch("fpt_antiproton", &fpt_antiproton, "fpt_antiproton/D");
  fTree_antiproton->Branch("fpTPC_antiproton", &fpTPC_antiproton, "fpTPC_antiproton/D");
  fTree_antiproton->Branch("fnSigmaTPC_antiproton", &fnSigmaTPC_antiproton, "fnSigmaTPC_antiproton/D");
  fTree_antiproton->Branch("fnSigmaTOF_antiproton", &fnSigmaTOF_antiproton, "fnSigmaTOF_antiproton/D");
  fTree_antiproton->Branch("fTrkMassTOF_antiproton", &fTrkMassTOF_antiproton, "fTrkMassTOF_antiproton/D");
  fTree_antiproton->Branch("fnSigmaITS_antiproton", &fnSigmaITS_antiproton, "fnSigmaITS_antiproton/D");
  fTree_antiproton->Branch("feta_antiproton", &feta_antiproton, "feta_antiproton/D");
  fTree_antiproton->Branch("fphi_antiproton", &fphi_antiproton, "fphi_antiproton/D");

  PostData(1,fOutputList);
  PostData(2,fTree_omegam);
  PostData(3,fTree_omegap);
  PostData(4,fTree_sidebandm);
  PostData(5,fTree_sidebandp);
  PostData(6,fTree_cascade);
  PostData(7,fTree_proton);
  PostData(8,fTree_antiproton);

}

//________________________________________________________________________
void AliAnalysisTaskpOmegaDibaryon::UserExec(Option_t *)
{
  //cout<<"++++++++++++++++++++"<<endl;
  //cout<<"++ Analysis Start ++"<<endl;
  //cout<<"++++++++++++++++++++"<<endl;

  // Input event
  //AliAnalysisManager   *man         =AliAnalysisManager::GetAnalysisManager();
  //AliInputEventHandler *inputHandler=(AliInputEventHandler*)(man->GetInputEventHandler());

  fEvent = dynamic_cast<AliVEvent*>(InputEvent());
  if(!fEvent){
    AliError("event is not available.");
    return;
  }
  fESD = dynamic_cast<AliESDEvent*>(fEvent);
  fAOD = dynamic_cast<AliAODEvent*>(fEvent);
 
  if(!fAOD){
    printf("ERROR: fAOD not available\n");
    return;
  }
  fHeader=(AliAODHeader*)fAOD->GetHeader();
  if(!fHeader){
    printf("ERROR: fHeader not available\n");
    return;
  }
  if(fInputHandler){
    fPIDResponse=fInputHandler->GetPIDResponse();
    if(!fPIDResponse){
      printf("ERROR: fPIDResponse not available\n");
      return;
    }
  }
  
  ////=============== Event selection
  isEvtMB=((((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected()) & ftrigBit);
  if(!isEvtMB) return;

  //printf("before Event selection\n");
  dynamic_cast<TH1F*>(fOutputList->FindObject("fEvtCounter"))->Fill(0.5);

  // Event pile-up
  isPileupEvt = fAOD->IsPileupFromSPD();
  isGoodEvt=EventSelection(fAOD);
  if(!isGoodEvt || isPileupEvt){
    PostData(1,fOutputList);
    PostData(2,fTree_omegam);
    PostData(3,fTree_omegap);
    PostData(4,fTree_sidebandm);
    PostData(5,fTree_sidebandp);
    PostData(6,fTree_cascade);
    PostData(7,fTree_proton);
    PostData(8,fTree_antiproton);
    return; 
  }
  dynamic_cast<TH1F*>(fOutputList->FindObject("fEvtPassCut"))->Fill(0.5);

  nevt ++;
  //cout<<"nevt = "<<nevt<<endl;

  centralityV0M = -1.;
  centralityCL0 = -1.;
  centralityCL1 = -1.;
  centralityV0A = -1.;
  centralityV0C = -1.;
  centralityZNA = -1.;
  centralityZNC = -1.;
  centralityMain = -1.;

  //Get Centrality
  fMultSelection = (AliMultSelection*)fEvent->FindListObject("MultSelection");
  if(!fMultSelection){
    //If you get this warning (and fCentralityV0M 300) please check that the AliMultSelectionTask actually ran (before your task)
    AliWarning("AliMultSelection object not found!");
    return;
  }
  else{
    centralityV0M  = fMultSelection->GetMultiplicityPercentile("V0M");
    centralityCL0  = fMultSelection->GetMultiplicityPercentile("CL0");
    centralityCL1  = fMultSelection->GetMultiplicityPercentile("CL1");
    centralityV0A  = fMultSelection->GetMultiplicityPercentile("V0A");
    centralityV0C  = fMultSelection->GetMultiplicityPercentile("V0C");
    centralityZNA  = fMultSelection->GetMultiplicityPercentile("ZNA");
    centralityZNC  = fMultSelection->GetMultiplicityPercentile("ZNC");
    centralityMain = fMultSelection->GetMultiplicityPercentile(fEstimator);
  }

  dynamic_cast<TH1F*>(fOutputList->FindObject("hCentralityV0M")) ->Fill(centralityV0M);
  dynamic_cast<TH1F*>(fOutputList->FindObject("hCentralityV0A")) ->Fill(centralityV0A);
  dynamic_cast<TH1F*>(fOutputList->FindObject("hCentralityV0C")) ->Fill(centralityV0C);
  dynamic_cast<TH1F*>(fOutputList->FindObject("hCentralityZNA")) ->Fill(centralityZNA);
  dynamic_cast<TH1F*>(fOutputList->FindObject("hCentralityZNC")) ->Fill(centralityZNC);
  dynamic_cast<TH1F*>(fOutputList->FindObject("hCentralityCL0")) ->Fill(centralityCL0);
  dynamic_cast<TH1F*>(fOutputList->FindObject("hCentralityCL1")) ->Fill(centralityCL1);
  dynamic_cast<TH1F*>(fOutputList->FindObject("hCentralityMain"))->Fill(centralityMain);

  /*
  if(centralityV0M >= 90){
    AliWarning("Centtality is high!");
    return;
  }
  */

  //Central
  if(centralityV0M > 10){
    AliWarning("Centrality is high!");
    return;
  }

  //Semi-central
  /*
  if(centralityV0M < 30){
    AliWarning("Centrality is high!");
    return;
  }

  if(centralityV0M > 50){
    AliWarning("Centrality is high!");
    return;
  }
  */
  dynamic_cast<TH1F*>(fOutputList->FindObject("hCentralityV0M_cut")) ->Fill(centralityV0M);

  //============= variable definition
  nTracks   =fAOD->GetNumberOfTracks();
  nV0s      =fAOD->GetNumberOfV0s();
  nCascades =fAOD->GetNumberOfCascades();

  vecTarget[3]={0.};
  vecTarget[0]=fAOD->GetPrimaryVertex()->GetX(); // primary vertex
  vecTarget[1]=fAOD->GetPrimaryVertex()->GetY();
  vecTarget[2]=fAOD->GetPrimaryVertex()->GetZ();

  kMagF =fAOD->GetMagneticField();
  KFParticle::SetField(kMagF);

  //=============== Vertex definition & initialization  
  //get initial primary vertex for secondary vertex using KF 
  vertex       = fAOD->GetPrimaryVertexSPD();
  primVertex   = AliAnalysisTaskpOmegaDibaryon::AODToESDVertex(*vertex);
  primKFVertex  = CreateKFVertex(*primVertex);
  PrimVertexKF[0]   = primKFVertex.GetX();
  PrimVertexKF[1]   = primKFVertex.GetY();
  PrimVertexKF[2]   = primKFVertex.GetZ();

  // __ initialize some utils __ //
  lorentzsum     = new TLorentzVector(0., 0., 0., 0.);
  particleProton = new TLorentzVector(0., 0., 0., 0.);
  particlePion   = new TLorentzVector(0., 0., 0., 0.);
  h              = new TVector3(0., 0., 0.);  

  lorentzsum_Xi    = new TLorentzVector(0., 0., 0., 0.);
  particleBachPion = new TLorentzVector(0., 0., 0., 0.);
  h_Xi             = new TVector3(0., 0., 0.);

  lorentzsum_Omega    = new TLorentzVector(0., 0., 0., 0.);
  particleBachKaon    = new TLorentzVector(0., 0., 0., 0.);

  // __ initialize __ //
  ProtonArray.clear(); 
  AntiProtonArray.clear();
  NPionArray.clear(); 
  PPionArray.clear();

  LambdadecayProtonArray.clear();
  LambdadecayAntiProtonArray.clear();
  LambdadecayPPionArray.clear();
  LambdadecayNPionArray.clear();
  XidecayBachPPionArray.clear();
  XidecayBachNPionArray.clear();
  XidecayBachPKaonArray.clear();
  XidecayBachNKaonArray.clear();

  np   = 0;
  nap  = 0;
  npi  = 0;
  nppi = 0;

  ndp  = 0;
  ndap = 0;
  ndpp = 0;
  ndnp = 0;
  nbpp = 0;
  nbnp = 0;
  nbpk = 0;
  nbnk = 0;

  exProtonTrack   = new AliExternalTrackParam();
  exPionTrack     = new AliExternalTrackParam();
  exBachPionTrack = new AliExternalTrackParam();
  exBachKaonTrack = new AliExternalTrackParam();

  //============== KF particle
  //========== daughter selection
  for(Int_t i=0; i<nTracks; i++){

    AliAODTrack* track=dynamic_cast<AliAODTrack*>(fAOD->GetTrack(i));
    if(!track) continue;
    dcaxy=0.,dcaz=0.,dca=0.; // proton 
    dcaxyn=0.,dcazn=0.,dcan=0.; // pion
    dDCA[2]={0.};              // DCA to the vertex xy(d) and z
    cDCA[3]={0.};              // convariance of impact parameters 
    v0Vtx[3]={0.};
    track->GetImpactParameters(dDCA,cDCA);
    dcaxy          =dDCA[0];
    dcaz           =dDCA[1];
    dca            =sqrt(dcaxy*dcaxy+dcaz*dcaz);

    track_pT     = track->Pt();
    pTCut_Proton = track_pT < 3.0;
    pTCut_Pion   = track_pT < 1.5;
    pTCut_Kaon   = track_pT < 2.0;
    
    isLambdadecayProton     = (fabs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton)) < 3 && track->Charge() > 0 && pTCut_Proton); 
    isLambdadecayAntiProton = (fabs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton)) < 3 && track->Charge() < 0 && pTCut_Proton); 
    isLambdadecayPPion      = (fabs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion)) < 3 && track->Charge() > 0 && pTCut_Pion);
    isLambdadecayNPion      = (fabs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion)) < 3 && track->Charge() < 0 && pTCut_Pion);  
    isBachPPion             = (fabs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion)) < 3 && track->Charge() > 0 && pTCut_Pion);
    isBachNPion             = (fabs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion)) < 3 && track->Charge() < 0 && pTCut_Pion);
    isBachPKaon             = (fabs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon)) < 3 && track->Charge() > 0 && pTCut_Kaon);
    isBachNKaon             = (fabs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon)) < 3 && track->Charge() < 0 && pTCut_Kaon);
    
    if(isLambdadecayProton){
      dynamic_cast<TH1F*>(fOutputList->FindObject("fpT_CascDecayProton_KF")) ->Fill(track->Pt());
    }

    if(isLambdadecayNPion){
      dynamic_cast<TH1F*>(fOutputList->FindObject("fpT_CascDecayPion_KF")) ->Fill(track->Pt());
    }

    if(!(IsQualityTrack(track))) continue;
    if(dca < 0.05)               continue;
    if(dca > 3)                  continue;
    
    // Xi
    if(isLambdadecayProton){ 
      LambdadecayProtonArray.push_back(i); 
      ndp++;
    }  
    if(isLambdadecayNPion){ 
      LambdadecayNPionArray.push_back(i); 
      ndnp++;
    }  
    if(isBachNKaon){ 
      XidecayBachNPionArray.push_back(i); 
      nbnp++;
    }  
    if(isLambdadecayAntiProton){ 
      LambdadecayAntiProtonArray.push_back(i); 
      ndap++;
    }  
    if(isLambdadecayPPion){ 
      LambdadecayPPionArray.push_back(i); 
      ndpp++;
    }  
   if(isBachPKaon){ 
      XidecayBachPPionArray.push_back(i); 
      nbpp++;
    }   

  }

  //============== KF Anti Lambda 
  for(Int_t i=0; i<ndap; i++){

    exProtonTrack->Reset();
    ProtonTrack=dynamic_cast<AliAODTrack*>(fAOD->GetTrack(LambdadecayAntiProtonArray[i])); 
    if(!ProtonTrack) continue;
    ftrkID_daughter1 = 0;
    ftrkID_daughter1 = ProtonTrack->GetID();
    if(ftrkID_daughter1 < 0) ftrkID_daughter1 =-ftrkID_daughter1-1;
    exProtonTrack->CopyFromVTrack(ProtonTrack);
    KFProtonTrack = AliAnalysisTaskpOmegaDibaryon::CreateKFParticle(*exProtonTrack, AliPID::ParticleMass(AliPID::kProton),-1);  
    if(TMath::Abs(KFProtonTrack.GetDistanceFromVertexXY(primKFVertex)) < 0.05) continue;
    if(TMath::Abs(KFProtonTrack.GetDistanceFromVertexXY(primKFVertex)) > 3.)   continue;
    
    for(Int_t j=0; j<ndpp; j++){

      exPionTrack->Reset(); 
      PionTrack=dynamic_cast<AliAODTrack*>(fAOD->GetTrack(LambdadecayPPionArray[j]));
      if(!PionTrack) continue;
      ftrkID_daughter2 = 0;
      ftrkID_daughter2 = PionTrack->GetID();
      if(ftrkID_daughter2 < 0) ftrkID_daughter2 =-ftrkID_daughter2-1;
      exPionTrack->CopyFromVTrack(PionTrack);
      KFPionTrack = AliAnalysisTaskpOmegaDibaryon::CreateKFParticle(*exPionTrack, AliPID::ParticleMass(AliPID::kPion),1);  
      
      if(TMath::Abs(KFPionTrack.GetDistanceFromVertexXY(primKFVertex)) < 0.05) continue;
      if(TMath::Abs(KFPionTrack.GetDistanceFromVertexXY(primKFVertex)) > 3.)   continue;

      KFParticleDibaryon KFSubMother;
      KFSubMother.ActivateWarnings();
      
      KFSubMother.AddDaughter(KFProtonTrack);
      if(!KFSubMother.CheckDaughter(KFPionTrack)) continue;
      KFSubMother.AddDaughter(KFPionTrack);
      
      //Transport the particle to its decay vertex
      KFSubMother.TransportToDecayVertex(); 
      
      SecVertexKF[0]     = KFSubMother.GetX();
      SecVertexKF[1]     = KFSubMother.GetY();
      SecVertexKF[2]     = KFSubMother.GetZ();
      SecVertexErrKF[0]  = KFSubMother.GetErrX();
      SecVertexErrKF[1]  = KFSubMother.GetErrY();
      SecVertexErrKF[2]  = KFSubMother.GetErrZ();      
      KFsecVertex   = new AliESDVertex(SecVertexKF, SecVertexErrKF);
      
      exProtonTrack->PropagateToDCA(KFsecVertex, kMagF, 10);
      exPionTrack->PropagateToDCA(KFsecVertex, kMagF, 10);
      if (KFsecVertex) delete KFsecVertex;
      
      fDcaSecDaughterProton  = TMath::Abs(exProtonTrack->GetD(SecVertexKF[0], SecVertexKF[1], kMagF));
      fDcaSecDaughterPion    = TMath::Abs(exPionTrack->GetD(SecVertexKF[0], SecVertexKF[1], kMagF));
      //fDCADaugter = fDcaSecDaughterProton + fDcaSecDaughterPion;
      fDCADaugter = KFProtonTrack.GetDistanceFromParticleXY(KFPionTrack);
      //if(fDCADaugter > 1.5)       continue;
      DCAtoPrimVtx  = (Float_t)KFSubMother.GetDistanceFromVertexXY(primKFVertex);
      //if(KFSubMother.Chi2() > 10) continue;
      //if(DCAtoPrimVtx < 0.2) continue;
      
      //cos(PA) 
      lorentzsum    ->SetXYZM(0., 0., 0., 0.);
      particleProton->SetXYZM(0., 0., 0., 0.);
      particlePion  ->SetXYZM(0., 0., 0., 0.);
      
      particleProton->SetXYZM(exProtonTrack->Px(), exProtonTrack->Py(), exProtonTrack->Pz(), AliPID::ParticleMass(AliPID::kProton));
      particlePion  ->SetXYZM(exPionTrack->Px(), exPionTrack->Py(), exPionTrack->Pz(), AliPID::ParticleMass(AliPID::kPion));
      *lorentzsum = *particleProton + *particlePion;

      dynamic_cast<TH1F*>(fOutputList->FindObject("fSecvertex_CascDecay_KF")) ->Fill(fabs(SecVertexKF[0]));
      dynamic_cast<TH1F*>(fOutputList->FindObject("fDCADaug_CascDecay_KF")) ->Fill(fDCADaugter);
      dynamic_cast<TH1F*>(fOutputList->FindObject("fDCAtoPrimVtx_CascDecay_KF")) ->Fill(DCAtoPrimVtx);
      dynamic_cast<TH1F*>(fOutputList->FindObject("fChi2_CascDecay_KF")) ->Fill(KFSubMother.Chi2());
      
      //if(fabs(SecVertexKF[0]) > 100) continue;
      //if(fabs(SecVertexKF[1]) > 100) continue;
      //if(fabs(SecVertexKF[2]) > 100) continue;
      //if(DCAtoPrimVtx < 0.07)  continue;

      dynamic_cast<TH1F*>(fOutputList->FindObject("fInvMassLambda_CascDecay_KF")) ->Fill(lorentzsum->M());
      if(TMath::Abs(lorentzsum->M() - TDatabasePDG::Instance()->GetParticle(3122)->Mass()) > 0.00624465) continue; 

      // KF Anti Xi
      for(Int_t k=0; k<nbpp; k++){

	//if(j==k) continue;
	BachNPionTrack=dynamic_cast<AliAODTrack*>(fAOD->GetTrack(XidecayBachPPionArray[k]));
	if(!BachNPionTrack) continue;
	ftrkID_daughter3 = 0;
	ftrkID_daughter3 = BachNPionTrack->GetID();
	if(ftrkID_daughter3 < 0) ftrkID_daughter3 =-ftrkID_daughter3-1;
	exBachPionTrack->CopyFromVTrack(BachNPionTrack);
	exBachKaonTrack->CopyFromVTrack(BachNPionTrack);
	KFBachPionTrack = AliAnalysisTaskpOmegaDibaryon::CreateKFParticle(*exBachPionTrack, AliPID::ParticleMass(AliPID::kPion),1); 
	KFBachKaonTrack = AliAnalysisTaskpOmegaDibaryon::CreateKFParticle(*exBachPionTrack, AliPID::ParticleMass(AliPID::kKaon),1); 
	
	if(TMath::Abs(KFBachPionTrack.GetDistanceFromVertexXY(primKFVertex)) < 0.05) continue;
	if(TMath::Abs(KFBachPionTrack.GetDistanceFromVertexXY(primKFVertex)) > 3.)   continue;
	if(TMath::Abs(KFBachKaonTrack.GetDistanceFromVertexXY(primKFVertex)) < 0.05) continue;
	if(TMath::Abs(KFBachKaonTrack.GetDistanceFromVertexXY(primKFVertex)) > 3.)   continue;
	
	KFParticleDibaryon KFMother;
	KFParticleDibaryon KFMother_Omega;

	KFMother.AddDaughter(KFSubMother);
	if(!KFMother.CheckDaughter(KFBachKaonTrack)) continue;
	KFMother.AddDaughter(KFBachKaonTrack);
	KFMother.TransportToDecayVertex();

	//KFMother_Omega.AddDaughter(KFSubMother);
	//if(!KFMother.CheckDaughter(KFBachKaonTrack)) continue;
	//KFMother_Omega.AddDaughter(KFBachKaonTrack);
	//KFMother_Omega.TransportToDecayVertex();
	
	SecVertexKF_Xi[0]     = KFMother.GetX();
	SecVertexKF_Xi[1]     = KFMother.GetY();
	SecVertexKF_Xi[2]     = KFMother.GetZ();
	SecVertexErrKF_Xi[0]  = KFMother.GetErrX();
	SecVertexErrKF_Xi[1]  = KFMother.GetErrY();
	SecVertexErrKF_Xi[2]  = KFMother.GetErrZ();	
	KFsecVertex_Xi   = new AliESDVertex(SecVertexKF_Xi, SecVertexErrKF_Xi);
	
	exBachPionTrack->PropagateToDCA(KFsecVertex_Xi, kMagF, 10);
	if (KFsecVertex_Xi) delete KFsecVertex_Xi;
	
	SecVertexKF_Omega[0]     = KFMother_Omega.GetX();
	SecVertexKF_Omega[1]     = KFMother_Omega.GetY();
	SecVertexKF_Omega[2]     = KFMother_Omega.GetZ();
	SecVertexErrKF_Omega[0]  = KFMother_Omega.GetErrX();
	SecVertexErrKF_Omega[1]  = KFMother_Omega.GetErrY();
	SecVertexErrKF_Omega[2]  = KFMother_Omega.GetErrZ();	
	KFsecVertex_Omega   = new AliESDVertex(SecVertexKF_Omega, SecVertexErrKF_Omega);
	
	exBachKaonTrack->PropagateToDCA(KFsecVertex_Omega, kMagF, 10);
	if (KFsecVertex_Omega) delete KFsecVertex_Omega;
	
	//cos(PA) 
	lorentzsum_Xi     ->SetXYZM(0., 0., 0., 0.);
	particleBachPion  ->SetXYZM(0., 0., 0., 0.);
	h_Xi              ->SetXYZ(0., 0., 0.);
	h                 ->SetXYZ(0., 0., 0.);	

	lorentzsum_Omega     ->SetXYZM(0., 0., 0., 0.);
	particleBachKaon     ->SetXYZM(0., 0., 0., 0.);
	
	particleBachPion ->SetXYZM(exBachPionTrack->Px(), exBachPionTrack->Py(), exBachPionTrack->Pz(), AliPID::ParticleMass(AliPID::kKaon));
	*lorentzsum_Xi = *lorentzsum + *particleBachPion;
	
	// new lambda CPA 
	dd[0] = SecVertexKF_Xi[0] - SecVertexKF[0];
	dd[1] = SecVertexKF_Xi[1] - SecVertexKF[1];
	dd[2] = SecVertexKF_Xi[2] - SecVertexKF[2];
	h->SetXYZ(-dd[0], -dd[1], -dd[2]);
	fPA          = TMath::Cos(lorentzsum->Angle(*h));
	DecLength    = TMath::Sqrt(TMath::Power(dd[0], 2) + TMath::Power(dd[1], 2) + TMath::Power(dd[2], 2));
	dynamic_cast<TH1F*>(fOutputList->FindObject("fCPA_CascDecay_KF")) ->Fill(fPA);
	dynamic_cast<TH1F*>(fOutputList->FindObject("fDecaylength_CascDecay_KF")) ->Fill(DecLength);
	if(fPA < 0.97)           continue;	
	//if(DecLength < 1.4)      continue;
	//if(DecLength > 200)      continue;

	dd_Xi[0] = PrimVertexKF[0] - SecVertexKF_Xi[0];
	dd_Xi[1] = PrimVertexKF[1] - SecVertexKF_Xi[1];
	dd_Xi[2] = PrimVertexKF[2] - SecVertexKF_Xi[2];
	h_Xi->SetXYZ(-dd_Xi[0], -dd_Xi[1], -dd_Xi[2]);
	
	particleBachKaon ->SetXYZM(exBachKaonTrack->Px(), exBachKaonTrack->Py(), exBachKaonTrack->Pz(), AliPID::ParticleMass(AliPID::kKaon));
	//*lorentzsum_Omega = *lorentzsum + *particleBachKaon;
	//if((lorentzsum_Omega->M() > 1.667) && (lorentzsum_Omega->M() < 1.677)) continue;   
	
	fPA_Xi = TMath::Cos(lorentzsum_Xi->Angle(*h_Xi));
	//Tranceverse decay length
	DecLength_Xi = TMath::Sqrt(TMath::Power(dd_Xi[0], 2) + TMath::Power(dd_Xi[1], 2));	
	DCAdauXY_Xi    = (Float_t)KFBachPionTrack.GetDistanceFromParticleXY(KFSubMother);
	DCAdauZ_Xi     = (Float_t)KFBachPionTrack.GetDistanceFromParticle(KFSubMother);
	DCAdau_Xi      = TMath::Sqrt(TMath::Power(DCAdauXY_Xi, 2) + TMath::Power(DCAdauZ_Xi, 2));
	
	dynamic_cast<TH1F*>(fOutputList->FindObject("fSecvertex_Casc_KF")) ->Fill(fabs(SecVertexKF_Xi[0]));
	dynamic_cast<TH1F*>(fOutputList->FindObject("fCPA_Casc_KF")) ->Fill(fPA_Xi);
	dynamic_cast<TH1F*>(fOutputList->FindObject("fTransradius_Casc_KF")) ->Fill(DecLength_Xi);
	dynamic_cast<TH1F*>(fOutputList->FindObject("fDCAdau_Casc_KF")) ->Fill(DCAdau_Xi);
	//dynamic_cast<TH1F*>(fOutputList->FindObject("fDCADaug_Casc_KF")) ->Fill(fDCADaugter);
	
	if(fDCADaugter > 1.5)      continue;
	//if(DecLength_Xi < 0.8)     continue;
	//if(DecLength_Xi > 200)     continue;
	if(fPA_Xi < 0.9)           continue;
	if(DCAdauXY_Xi > 1.5)      continue;
	
	dynamic_cast<TH1F*>(fOutputList->FindObject("fInvMassXi_Casc_KF")) ->Fill(lorentzsum_Xi->M());
	if(TMath::Abs(lorentzsum_Xi->M() - TDatabasePDG::Instance()->GetParticle(3334)->Mass()) < 0.00907){

	  dynamic_cast<TH1F*>(fOutputList->FindObject("fpT_Casc_AntiXi_KF")) ->Fill(lorentzsum_Xi->Pt());
	  dynamic_cast<TH1F*>(fOutputList->FindObject("fpT_CascBach_PPion_KF")) ->Fill(particleBachPion->Pt());
	  dynamic_cast<TH1F*>(fOutputList->FindObject("fpT_CascDecay_AntiLambda_KF")) ->Fill(lorentzsum->Pt());	
	  dynamic_cast<TH1F*>(fOutputList->FindObject("fEta_omega")) ->Fill(lorentzsum->Eta());	
	  dynamic_cast<TH1F*>(fOutputList->FindObject("fPhi_omega")) ->Fill(lorentzsum->Phi());	

	  //== TTree_omegap
	  fevtID_omegap         = fHeader->GetEventIdAsLong();
	  fzvertex_omegap       = vertex->GetZ();
	  fcentrality_omegap    = centralityV0M;
	  //=== Daughter track
	  fdca_daughter1_omegap   = TMath::Abs(KFProtonTrack.GetDistanceFromVertexXY(primKFVertex));
	  fdca_daughter2_omegap   = TMath::Abs(KFPionTrack.GetDistanceFromVertexXY(primKFVertex));
	  fdca_daughter3_omegap   = TMath::Abs(KFBachPionTrack.GetDistanceFromVertexXY(primKFVertex));
	  ftrkID_daughter1_omegap = ftrkID_daughter1;	
	  ftrkID_daughter2_omegap = ftrkID_daughter2;
	  ftrkID_daughter3_omegap = ftrkID_daughter3;

	  //=== Lambda
	  fcpa_lambda_omegap       = fPA;
	  fpcaxy_lambda_omegap     = KFProtonTrack.GetDistanceFromParticleXY(KFPionTrack);
	  fpca_lambda_omegap       = KFProtonTrack.GetDistanceFromParticle(KFPionTrack);
	  fdeclength_lambda_omegap = DecLength;
	  fchi2_lambda_omegap      = KFSubMother.Chi2();
	  fmass_lamdba_omegap      = lorentzsum->M();
	  fdcatoPVxy_lambda_omegap = fabs((Float_t)KFSubMother.GetDistanceFromVertexXY(primKFVertex));
	  fdcatoPV_lambda_omegap   = fabs((Float_t)KFSubMother.GetDistanceFromVertex(primKFVertex));
	  fpt_lambda_omegap        = lorentzsum->Pt();
	  
	  //=== Omega
	  fcpa_omega_omegap       = fPA_Xi;
	  fpcaxy_omega_omegap     = KFBachPionTrack.GetDistanceFromParticleXY(KFSubMother);
	  fpca_omega_omegap       = KFBachPionTrack.GetDistanceFromParticle(KFSubMother);
	  fdeclength_omega_omegap = DecLength_Xi;
	  fchi2_omega_omegap      = KFMother.Chi2();
	  fmass_omega_omegap      = lorentzsum_Xi->M();
	  fpx_omega_omegap        = lorentzsum_Xi->Px();
	  fpy_omega_omegap        = lorentzsum_Xi->Py();
	  fpz_omega_omegap        = lorentzsum_Xi->Pz();
	  fpt_omega_omegap        = lorentzsum_Xi->Pt();
	  feta_omega_omegap       = lorentzsum_Xi->Eta();
	  fphi_omega_omegap       = lorentzsum_Xi->Phi();
	  fdcatoPVxy_omega_omegap = fabs((Float_t)KFMother.GetDistanceFromVertexXY(primKFVertex));
	  fdcatoPV_omega_omegap   = fabs((Float_t)KFMother.GetDistanceFromVertex(primKFVertex));
	  fct_omega_omegap        = (lorentzsum_Xi->M() * TMath::Sqrt(TMath::Power(dd_Xi[0], 2) + TMath::Power(dd_Xi[1], 2) + TMath::Power(dd_Xi[2], 2))) / lorentzsum_Xi->P();
	  fpx_omega_omegap_mc     = -99;
	  fpy_omega_omegap_mc     = -99;
	  fpz_omega_omegap_mc     = -99;
	  fpt_omega_omegap_mc     = -99;

	  fTree_omegap->Fill();

	}
	
	else if(TMath::Abs(lorentzsum_Xi->M() - TDatabasePDG::Instance()->GetParticle(3334)->Mass()) < 0.010884){

	  //== TTree sidebandpm
	  fevtID_sidebandp         = fHeader->GetEventIdAsLong();
	  fzvertex_sidebandp       = vertex->GetZ();
	  fcentrality_sidebandp    = centralityV0M;
          //=== Daughter track
          fdca_daughter1_sidebandp   = TMath::Abs(KFProtonTrack.GetDistanceFromVertexXY(primKFVertex));
          fdca_daughter2_sidebandp   = TMath::Abs(KFPionTrack.GetDistanceFromVertexXY(primKFVertex));
          fdca_daughter3_sidebandp   = TMath::Abs(KFBachPionTrack.GetDistanceFromVertexXY(primKFVertex));
	  ftrkID_daughter1_sidebandp = ftrkID_daughter1;	
	  ftrkID_daughter2_sidebandp = ftrkID_daughter2;
	  ftrkID_daughter3_sidebandp = ftrkID_daughter3;

          //=== Lambda
          fcpa_lambda_sidebandp       = fPA;
          fpcaxy_lambda_sidebandp     = KFProtonTrack.GetDistanceFromParticleXY(KFPionTrack);
          fpca_lambda_sidebandp       = KFProtonTrack.GetDistanceFromParticle(KFPionTrack);
          fdeclength_lambda_sidebandp = DecLength;
          fchi2_lambda_sidebandp      = KFSubMother.Chi2();
          fmass_lamdba_sidebandp      = lorentzsum->M();
          fdcatoPVxy_lambda_sidebandp = fabs((Float_t)KFSubMother.GetDistanceFromVertexXY(primKFVertex));
          fdcatoPV_lambda_sidebandp   = fabs((Float_t)KFSubMother.GetDistanceFromVertex(primKFVertex));
	  fpt_lambda_sidebandp        = lorentzsum->Pt();

          //=== Omega
          fcpa_omega_sidebandp       = fPA_Xi;
          fpcaxy_omega_sidebandp     = KFBachPionTrack.GetDistanceFromParticleXY(KFSubMother);
          fpca_omega_sidebandp       = KFBachPionTrack.GetDistanceFromParticle(KFSubMother);
          fdeclength_omega_sidebandp = DecLength_Xi;
          fchi2_omega_sidebandp      = KFMother.Chi2();
          fmass_omega_sidebandp      = lorentzsum_Xi->M();
	  fpx_omega_sidebandp        = lorentzsum_Xi->Px();
	  fpy_omega_sidebandp        = lorentzsum_Xi->Py();
	  fpz_omega_sidebandp        = lorentzsum_Xi->Pz();
	  fpt_omega_sidebandp        = lorentzsum_Xi->Pt();
	  feta_omega_sidebandp       = lorentzsum_Xi->Eta();
	  fphi_omega_sidebandp       = lorentzsum_Xi->Phi();
	  fdcatoPVxy_omega_sidebandp = fabs((Float_t)KFMother.GetDistanceFromVertexXY(primKFVertex));
	  fdcatoPV_omega_sidebandp   = fabs((Float_t)KFMother.GetDistanceFromVertex(primKFVertex));
	  fct_omega_sidebandp        = (lorentzsum_Xi->M() * TMath::Sqrt(TMath::Power(dd_Xi[0], 2) + TMath::Power(dd_Xi[1], 2) + TMath::Power(dd_Xi[2], 2))) / lorentzsum_Xi->P();
	  fpx_omega_sidebandp_mc     = -99;
	  fpy_omega_sidebandp_mc     = -99;
	  fpz_omega_sidebandp_mc     = -99;
	  fpt_omega_sidebandp_mc     = -99;

          fTree_sidebandp->Fill();
        }

	KFMother.Delete();
	KFMother_Omega.Delete();
      }
      KFSubMother.Delete();
    }
    
  }
  
  //============== KF Lambda
  for(Int_t i=0; i<ndp; i++){

    exProtonTrack->Reset();
    ProtonTrack=dynamic_cast<AliAODTrack*>(fAOD->GetTrack(LambdadecayProtonArray[i])); 
    if(!ProtonTrack) continue;
    ftrkID_daughter1 = 0;
    ftrkID_daughter1 = ProtonTrack->GetID();	
    if(ftrkID_daughter1 < 0) ftrkID_daughter1 =-ftrkID_daughter1-1;
    exProtonTrack->CopyFromVTrack(ProtonTrack);
    KFProtonTrack = AliAnalysisTaskpOmegaDibaryon::CreateKFParticle(*exProtonTrack, AliPID::ParticleMass(AliPID::kProton),1);  
    if(TMath::Abs(KFProtonTrack.GetDistanceFromVertexXY(primKFVertex)) < 0.05) continue;
    if(TMath::Abs(KFProtonTrack.GetDistanceFromVertexXY(primKFVertex)) > 3.)   continue;    

    for(Int_t j=0; j<ndnp; j++){

      exPionTrack->Reset(); 
      PionTrack=dynamic_cast<AliAODTrack*>(fAOD->GetTrack(LambdadecayNPionArray[j]));
      if(!PionTrack) continue;
      ftrkID_daughter2 = 0;
      ftrkID_daughter2 = PionTrack->GetID();	
      if(ftrkID_daughter2 < 0) ftrkID_daughter2 =-ftrkID_daughter2-1;
      exPionTrack->CopyFromVTrack(PionTrack);
      KFPionTrack = AliAnalysisTaskpOmegaDibaryon::CreateKFParticle(*exPionTrack, AliPID::ParticleMass(AliPID::kPion),-1);  
      if(TMath::Abs(KFPionTrack.GetDistanceFromVertexXY(primKFVertex)) < 0.05) continue;
      if(TMath::Abs(KFPionTrack.GetDistanceFromVertexXY(primKFVertex)) > 3.)   continue;      

      KFParticleDibaryon KFSubMother;
      KFSubMother.ActivateWarnings();

      KFSubMother.AddDaughter(KFProtonTrack);
      if(!KFSubMother.CheckDaughter(KFPionTrack)) continue;
      KFSubMother.AddDaughter(KFPionTrack);
      
      //Transport the particle to its decay vertex
      KFSubMother.TransportToDecayVertex(); 
      
      SecVertexKF[0]     = KFSubMother.GetX();
      SecVertexKF[1]     = KFSubMother.GetY();
      SecVertexKF[2]     = KFSubMother.GetZ();
      SecVertexErrKF[0]  = KFSubMother.GetErrX();
      SecVertexErrKF[1]  = KFSubMother.GetErrY();
      SecVertexErrKF[2]  = KFSubMother.GetErrZ();
      KFsecVertex   = new AliESDVertex(SecVertexKF, SecVertexErrKF);
      
      exProtonTrack->PropagateToDCA(KFsecVertex, kMagF, 10);
      exPionTrack->PropagateToDCA(KFsecVertex, kMagF, 10);
      if (KFsecVertex) delete KFsecVertex;
      
      fDcaSecDaughterProton  = TMath::Abs(exProtonTrack->GetD(SecVertexKF[0], SecVertexKF[1], kMagF));
      fDcaSecDaughterPion    = TMath::Abs(exPionTrack->GetD(SecVertexKF[0], SecVertexKF[1], kMagF));
      //fDCADaugter = fDcaSecDaughterProton + fDcaSecDaughterPion;
      fDCADaugter = KFProtonTrack.GetDistanceFromParticleXY(KFPionTrack);
      //if(fDCADaugter > 1.5)       continue;
      DCAtoPrimVtx  = (Float_t)KFSubMother.GetDistanceFromVertexXY(primKFVertex);      
      //if(KFSubMother.Chi2() > 10) continue;
      //if(DCAtoPrimVtx < 0.2) continue;
      
      //cos(PA) 
      lorentzsum    ->SetXYZM(0., 0., 0., 0.);
      particleProton->SetXYZM(0., 0., 0., 0.);
      particlePion  ->SetXYZM(0., 0., 0., 0.);
      
      particleProton->SetXYZM(exProtonTrack->Px(), exProtonTrack->Py(), exProtonTrack->Pz(), AliPID::ParticleMass(AliPID::kProton));
      particlePion  ->SetXYZM(exPionTrack->Px(), exPionTrack->Py(), exPionTrack->Pz(), AliPID::ParticleMass(AliPID::kPion));
      *lorentzsum = *particleProton + *particlePion;

      dynamic_cast<TH1F*>(fOutputList->FindObject("fSecvertex_CascDecay_KF")) ->Fill(fabs(SecVertexKF[0]));
      dynamic_cast<TH1F*>(fOutputList->FindObject("fDCADaug_CascDecay_KF")) ->Fill(fDCADaugter);
      dynamic_cast<TH1F*>(fOutputList->FindObject("fDCAtoPrimVtx_CascDecay_KF")) ->Fill(DCAtoPrimVtx);
      dynamic_cast<TH1F*>(fOutputList->FindObject("fChi2_CascDecay_KF")) ->Fill(KFSubMother.Chi2());
      
      //if(fabs(SecVertexKF[0]) > 100) continue;
      //if(fabs(SecVertexKF[1]) > 100) continue;
      //if(fabs(SecVertexKF[2]) > 100) continue;
      //if(DCAtoPrimVtx < 0.07)  continue;

      dynamic_cast<TH1F*>(fOutputList->FindObject("fInvMassLambda_CascDecay_KF")) ->Fill(lorentzsum->M());
      if(TMath::Abs(lorentzsum->M() - TDatabasePDG::Instance()->GetParticle(3122)->Mass()) > 0.00624465) continue; 

      // KF Xi
      for(Int_t k=0; k<nbnp; k++){
	//if(j==k) continue;
	BachNPionTrack=dynamic_cast<AliAODTrack*>(fAOD->GetTrack(XidecayBachNPionArray[k]));
	if(!BachNPionTrack) continue;
	ftrkID_daughter3 = 0;
	ftrkID_daughter3 = BachNPionTrack->GetID();	
	if(ftrkID_daughter3 < 0) ftrkID_daughter3 =-ftrkID_daughter3-1;
	//cout<<"+++++++++++ Bach pion Px = "<<XidecayBachNPionArray[k]<<":"<<BachNPionTrack->Px()<<"+++++++++++"<<endl;   
	//cout<<"Xi!"<<endl;
	exBachPionTrack->CopyFromVTrack(BachNPionTrack);
	exBachKaonTrack->CopyFromVTrack(BachNPionTrack);
	KFBachPionTrack = AliAnalysisTaskpOmegaDibaryon::CreateKFParticle(*exBachPionTrack, AliPID::ParticleMass(AliPID::kPion),-1); 
	KFBachKaonTrack = AliAnalysisTaskpOmegaDibaryon::CreateKFParticle(*exBachPionTrack, AliPID::ParticleMass(AliPID::kKaon),-1); 
	if(TMath::Abs(KFBachPionTrack.GetDistanceFromVertexXY(primKFVertex)) < 0.05) continue;
	if(TMath::Abs(KFBachPionTrack.GetDistanceFromVertexXY(primKFVertex)) > 3.)   continue;      
	if(TMath::Abs(KFBachKaonTrack.GetDistanceFromVertexXY(primKFVertex)) < 0.05) continue;
	if(TMath::Abs(KFBachKaonTrack.GetDistanceFromVertexXY(primKFVertex)) > 3.)   continue;      	

	KFParticleDibaryon KFMother;
	KFParticleDibaryon KFMother_Omega;
	//KFSubMother.SetConstructMethod(2);

	KFMother.AddDaughter(KFSubMother);
	if(!KFMother.CheckDaughter(KFBachKaonTrack)) continue;
	KFMother.AddDaughter(KFBachKaonTrack);
	KFMother.TransportToDecayVertex();

	//KFMother_Omega.AddDaughter(KFSubMother);
	//if(!KFMother_Omega.CheckDaughter(KFBachKaonTrack)) continue;
	//KFMother_Omega.AddDaughter(KFBachKaonTrack);
	//KFMother_Omega.TransportToDecayVertex();
	
	SecVertexKF_Xi[0]     = KFMother.GetX();
	SecVertexKF_Xi[1]     = KFMother.GetY();
	SecVertexKF_Xi[2]     = KFMother.GetZ();
	SecVertexErrKF_Xi[0]  = KFMother.GetErrX();
	SecVertexErrKF_Xi[1]  = KFMother.GetErrY();
	SecVertexErrKF_Xi[2]  = KFMother.GetErrZ();	
	KFsecVertex_Xi   = new AliESDVertex(SecVertexKF_Xi, SecVertexErrKF_Xi);
	
	exBachPionTrack->PropagateToDCA(KFsecVertex_Xi, kMagF, 10);
	if (KFsecVertex_Xi) delete KFsecVertex_Xi;
	
	SecVertexKF_Omega[0]     = KFMother_Omega.GetX();
	SecVertexKF_Omega[1]     = KFMother_Omega.GetY();
	SecVertexKF_Omega[2]     = KFMother_Omega.GetZ();
	SecVertexErrKF_Omega[0]  = KFMother_Omega.GetErrX();
	SecVertexErrKF_Omega[1]  = KFMother_Omega.GetErrY();
	SecVertexErrKF_Omega[2]  = KFMother_Omega.GetErrZ();	
	KFsecVertex_Omega   = new AliESDVertex(SecVertexKF_Omega, SecVertexErrKF_Omega);
	
	exBachKaonTrack->PropagateToDCA(KFsecVertex_Omega, kMagF, 10);
	if (KFsecVertex_Omega) delete KFsecVertex_Omega;
	
	//cos(PA) 
	lorentzsum_Xi     ->SetXYZM(0., 0., 0., 0.);
	particleBachPion  ->SetXYZM(0., 0., 0., 0.);
	h_Xi              ->SetXYZ(0., 0., 0.);
	h                 ->SetXYZ(0., 0., 0.);	

	lorentzsum_Omega     ->SetXYZM(0., 0., 0., 0.);
	particleBachKaon     ->SetXYZM(0., 0., 0., 0.);
	
	particleBachPion ->SetXYZM(exBachPionTrack->Px(), exBachPionTrack->Py(), exBachPionTrack->Pz(), AliPID::ParticleMass(AliPID::kKaon));
	*lorentzsum_Xi = *lorentzsum + *particleBachPion;
	
	// new lambda CPA 
	dd[0] = SecVertexKF_Xi[0] - SecVertexKF[0];
	dd[1] = SecVertexKF_Xi[1] - SecVertexKF[1];
	dd[2] = SecVertexKF_Xi[2] - SecVertexKF[2];
	h->SetXYZ(-dd[0], -dd[1], -dd[2]);
	fPA          = TMath::Cos(lorentzsum->Angle(*h));
	DecLength    = TMath::Sqrt(TMath::Power(dd[0], 2) + TMath::Power(dd[1], 2) + TMath::Power(dd[2], 2));
	dynamic_cast<TH1F*>(fOutputList->FindObject("fCPA_CascDecay_KF")) ->Fill(fPA);
	dynamic_cast<TH1F*>(fOutputList->FindObject("fDecaylength_CascDecay_KF")) ->Fill(DecLength);
	//if(i==62 && j==56 && k==447) cout<<"fPA = "<<fPA<<", nbnp = "<<nbnp<<endl;
	if(fPA < 0.97)           continue;	
	//if(DecLength < 1.4)      continue;
	//if(DecLength > 200)      continue;

	dd_Xi[0] = PrimVertexKF[0] - SecVertexKF_Xi[0];
	dd_Xi[1] = PrimVertexKF[1] - SecVertexKF_Xi[1];
	dd_Xi[2] = PrimVertexKF[2] - SecVertexKF_Xi[2];
	h_Xi->SetXYZ(-dd_Xi[0], -dd_Xi[1], -dd_Xi[2]);
	
	particleBachKaon ->SetXYZM(exBachKaonTrack->Px(), exBachKaonTrack->Py(), exBachKaonTrack->Pz(), AliPID::ParticleMass(AliPID::kKaon));
	//*lorentzsum_Omega = *lorentzsum + *particleBachKaon;
	//if((lorentzsum_Omega->M() > 1.667) && (lorentzsum_Omega->M() < 1.677)) continue;   
	
	fPA_Xi = TMath::Cos(lorentzsum_Xi->Angle(*h_Xi));
	//Tranceverse decay length
	DecLength_Xi = TMath::Sqrt(TMath::Power(dd_Xi[0], 2) + TMath::Power(dd_Xi[1], 2));	
	DCAdauXY_Xi    = (Float_t)KFBachPionTrack.GetDistanceFromParticleXY(KFSubMother);
	DCAdauZ_Xi     = (Float_t)KFBachPionTrack.GetDistanceFromParticle(KFSubMother);
	DCAdau_Xi      = TMath::Sqrt(TMath::Power(DCAdauXY_Xi, 2) + TMath::Power(DCAdauZ_Xi, 2));
	
	dynamic_cast<TH1F*>(fOutputList->FindObject("fSecvertex_Casc_KF")) ->Fill(fabs(SecVertexKF_Xi[0]));
	dynamic_cast<TH1F*>(fOutputList->FindObject("fCPA_Casc_KF")) ->Fill(fPA_Xi);
	dynamic_cast<TH1F*>(fOutputList->FindObject("fTransradius_Casc_KF")) ->Fill(DecLength_Xi);
	dynamic_cast<TH1F*>(fOutputList->FindObject("fDCAdau_Casc_KF")) ->Fill(DCAdau_Xi);
	//dynamic_cast<TH1F*>(fOutputList->FindObject("fDCADaug_Casc_KF")) ->Fill(fDCADaugter);
	
	if(fDCADaugter > 1.5)      continue;
	//if(DecLength_Xi < 0.8)     continue;
	//if(DecLength_Xi > 200)     continue;
	if(fPA_Xi < 0.9)           continue;
	if(DCAdauXY_Xi > 1.5)      continue;
	
	dynamic_cast<TH1F*>(fOutputList->FindObject("fInvMassXi_Casc_KF")) ->Fill(lorentzsum_Xi->M());
	if(TMath::Abs(lorentzsum_Xi->M() - TDatabasePDG::Instance()->GetParticle(3334)->Mass()) < 0.00907){
	  dynamic_cast<TH1F*>(fOutputList->FindObject("fpT_Casc_AntiXi_KF")) ->Fill(lorentzsum_Xi->Pt());
	  dynamic_cast<TH1F*>(fOutputList->FindObject("fpT_CascBach_PPion_KF")) ->Fill(particleBachPion->Pt());
	  dynamic_cast<TH1F*>(fOutputList->FindObject("fpT_CascDecay_AntiLambda_KF")) ->Fill(lorentzsum->Pt());	
	  dynamic_cast<TH1F*>(fOutputList->FindObject("fEta_omega")) ->Fill(lorentzsum->Eta());	
	  dynamic_cast<TH1F*>(fOutputList->FindObject("fPhi_omega")) ->Fill(lorentzsum->Phi());		  

	  //== TTree_omegam
	  fevtID_omegam         = fHeader->GetEventIdAsLong();
	  fzvertex_omegam       = vertex->GetZ();
	  fcentrality_omegam    = centralityV0M;
	  //=== Daughter track
	  fdca_daughter1_omegam   = TMath::Abs(KFProtonTrack.GetDistanceFromVertexXY(primKFVertex));
	  fdca_daughter2_omegam   = TMath::Abs(KFPionTrack.GetDistanceFromVertexXY(primKFVertex));
	  fdca_daughter3_omegam   = TMath::Abs(KFBachPionTrack.GetDistanceFromVertexXY(primKFVertex));
	  ftrkID_daughter1_omegam = ftrkID_daughter1;	
	  ftrkID_daughter2_omegam = ftrkID_daughter2;
	  ftrkID_daughter3_omegam = ftrkID_daughter3;

	  //=== Lambda
	  fcpa_lambda_omegam       = fPA;
	  fpcaxy_lambda_omegam     = KFProtonTrack.GetDistanceFromParticleXY(KFPionTrack);
	  fpca_lambda_omegam       = KFProtonTrack.GetDistanceFromParticle(KFPionTrack);
	  fdeclength_lambda_omegam = DecLength;
	  fchi2_lambda_omegam      = KFSubMother.Chi2();
	  fmass_lamdba_omegam      = lorentzsum->M();
	  fdcatoPVxy_lambda_omegam = fabs((Float_t)KFSubMother.GetDistanceFromVertexXY(primKFVertex));
	  fdcatoPV_lambda_omegam   = fabs((Float_t)KFSubMother.GetDistanceFromVertex(primKFVertex));
	  fpt_lambda_omegam        = lorentzsum->Pt();
	  
	  //=== Omega
	  fcpa_omega_omegam       = fPA_Xi;
	  fpcaxy_omega_omegam     = KFBachPionTrack.GetDistanceFromParticleXY(KFSubMother);
	  fpca_omega_omegam       = KFBachPionTrack.GetDistanceFromParticle(KFSubMother);
	  fdeclength_omega_omegam = DecLength_Xi;
	  fchi2_omega_omegam      = KFMother.Chi2();
	  fmass_omega_omegam      = lorentzsum_Xi->M();
	  fpx_omega_omegam        = lorentzsum_Xi->Px();
	  fpy_omega_omegam        = lorentzsum_Xi->Py();
	  fpz_omega_omegam        = lorentzsum_Xi->Pz();
	  fpt_omega_omegam        = lorentzsum_Xi->Pt();
	  feta_omega_omegam       = lorentzsum_Xi->Eta();
	  fphi_omega_omegam       = lorentzsum_Xi->Phi();
	  fdcatoPVxy_omega_omegam = fabs((Float_t)KFMother.GetDistanceFromVertexXY(primKFVertex));
	  fdcatoPV_omega_omegam   = fabs((Float_t)KFMother.GetDistanceFromVertex(primKFVertex));
	  fct_omega_omegam        = (lorentzsum_Xi->M() * TMath::Sqrt(TMath::Power(dd_Xi[0], 2) + TMath::Power(dd_Xi[1], 2) + TMath::Power(dd_Xi[2], 2))) / lorentzsum_Xi->P();
	  fpx_omega_omegam_mc     = -99;
	  fpy_omega_omegam_mc     = -99;
	  fpz_omega_omegam_mc     = -99;
	  fpt_omega_omegam_mc     = -99;

	  fTree_omegam->Fill();

	}
	
	else if(TMath::Abs(lorentzsum_Xi->M() - TDatabasePDG::Instance()->GetParticle(3334)->Mass()) < 0.010884){

	  //== TTree sidebandm
	  fevtID_sidebandm         = fHeader->GetEventIdAsLong();
	  fzvertex_sidebandm       = vertex->GetZ();
	  fcentrality_sidebandm    = centralityV0M;
          //=== Daughter track
          fdca_daughter1_sidebandm   = TMath::Abs(KFProtonTrack.GetDistanceFromVertexXY(primKFVertex));
          fdca_daughter2_sidebandm   = TMath::Abs(KFPionTrack.GetDistanceFromVertexXY(primKFVertex));
          fdca_daughter3_sidebandm   = TMath::Abs(KFBachPionTrack.GetDistanceFromVertexXY(primKFVertex));
	  ftrkID_daughter1_sidebandm = ftrkID_daughter1;	
	  ftrkID_daughter2_sidebandm = ftrkID_daughter2;
	  ftrkID_daughter3_sidebandm = ftrkID_daughter3;

          //=== Lambda
          fcpa_lambda_sidebandm       = fPA;
          fpcaxy_lambda_sidebandm     = KFProtonTrack.GetDistanceFromParticleXY(KFPionTrack);
          fpca_lambda_sidebandm       = KFProtonTrack.GetDistanceFromParticle(KFPionTrack);
          fdeclength_lambda_sidebandm = DecLength;
          fchi2_lambda_sidebandm      = KFSubMother.Chi2();
          fmass_lamdba_sidebandm      = lorentzsum->M();
          fdcatoPVxy_lambda_sidebandm = fabs((Float_t)KFSubMother.GetDistanceFromVertexXY(primKFVertex));
          fdcatoPV_lambda_sidebandm   = fabs((Float_t)KFSubMother.GetDistanceFromVertex(primKFVertex));
	  fpt_lambda_sidebandm        = lorentzsum->Pt();

          //=== Omega
          fcpa_omega_sidebandm       = fPA_Xi;
          fpcaxy_omega_sidebandm     = KFBachPionTrack.GetDistanceFromParticleXY(KFSubMother);
          fpca_omega_sidebandm       = KFBachPionTrack.GetDistanceFromParticle(KFSubMother);
          fdeclength_omega_sidebandm = DecLength_Xi;
          fchi2_omega_sidebandm      = KFMother.Chi2();
          fmass_omega_sidebandm      = lorentzsum_Xi->M();
	  fpx_omega_sidebandm        = lorentzsum_Xi->Px();
	  fpy_omega_sidebandm        = lorentzsum_Xi->Py();
	  fpz_omega_sidebandm        = lorentzsum_Xi->Pz();
	  fpt_omega_sidebandm        = lorentzsum_Xi->Pt();
	  feta_omega_sidebandm       = lorentzsum_Xi->Eta();
	  fphi_omega_sidebandm       = lorentzsum_Xi->Phi();
	  fdcatoPVxy_omega_sidebandm = fabs((Float_t)KFMother.GetDistanceFromVertexXY(primKFVertex));
	  fdcatoPV_omega_sidebandm   = fabs((Float_t)KFMother.GetDistanceFromVertex(primKFVertex));
	  fct_omega_sidebandm        = (lorentzsum_Xi->M() * TMath::Sqrt(TMath::Power(dd_Xi[0], 2) + TMath::Power(dd_Xi[1], 2) + TMath::Power(dd_Xi[2], 2))) / lorentzsum_Xi->P();
	  fpx_omega_sidebandm_mc     = -99;
	  fpy_omega_sidebandm_mc     = -99;
	  fpz_omega_sidebandm_mc     = -99;
	  fpt_omega_sidebandm_mc     = -99;

          fTree_sidebandm->Fill();
        }

	KFMother.Delete();
	KFMother_Omega.Delete();
      }
      KFSubMother.Delete();
    }
    
  }

  //============== proton loop
  for(Int_t i=0; i<nTracks; i++){
    AliAODTrack* track=dynamic_cast<AliAODTrack*>(fAOD->GetTrack(i));
    if(!track) continue;
    //if(!track->TestFilterBit(128)) continue;

    dcaxy=0.,dcaz=0.,dca=0.;
    dDCA[2]={0.};              // DCA to the vertex xy(d) and z
    cDCA[3]={0.};              // convariance of impact parameters 
    v0Vtx[3]={0.};
    track->GetImpactParameters(dDCA,cDCA);
    dcaxy          =dDCA[0];
    dcaz           =dDCA[1];
    dca            =sqrt(dcaxy*dcaxy+dcaz*dcaz);
    track_pT       =track->Pt();

    //proton loose cuts
    if(!(IsQualityProton(track)))       continue;
    if(fabs(dcaxy) > 0.1)               continue;
    if(fabs(dcaz)  > 0.2)               continue;
    dynamic_cast<TH2F*>(fOutputList->FindObject("fPIDITS_proton"))
      ->Fill(track->GetTPCmomentum(),track->GetITSsignal());

    //Layer hit requirements
    if(!(track->HasPointOnITSLayer(0)) && !(track->HasPointOnITSLayer(1))) continue;

    //PID
    if(fabs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton)) > 3.0) continue; //TPC
    if((track->GetTPCmomentum() < 0.7) && (track->GetTOFsignal()!=99999)){
      if(fabs(CalculateProtonSigma(track->Pt(),CalculateMassSquareTOF(*track),track->Charge())) > 3.0) continue; //TOF
    }
    if(track->GetTPCmomentum() > 0.7){
      if(fabs(CalculateProtonSigma(track->Pt(),CalculateMassSquareTOF(*track),track->Charge())) > 3.0) continue; //TOF
    }

    //=== ITS PID
    AliPIDResponse::EDetPidStatus statusITS = fPIDResponse->CheckPIDStatus(AliPIDResponse::kITS,track);
    Bool_t ITSisOK = false;
    if(statusITS == AliPIDResponse::kDetPidOk) ITSisOK = true;
    if(ITSisOK){
      Double_t ITS_nCluster    = track->GetITSNcls();
      if(ITS_nCluster < 0) ITS_nCluster = 0;
      if(ITS_nCluster < 2) continue;
      Double_t ITS_dEdx        = track->GetITSsignal();
      Double_t ITS_dEdx_nSigma = CalculateSigmadEdxITS(*track);
      dynamic_cast<TH1F*>(fOutputList->FindObject("fSigmaITS_Proton"))
	->Fill(ITS_dEdx_nSigma);
      if(fabs(ITS_dEdx_nSigma) > 3) continue;
      dynamic_cast<TH1F*>(fOutputList->FindObject("fSigmaITS_Proton_cut"))
	->Fill(ITS_dEdx_nSigma);
    }
    
    if(track->Charge()>0) {
      Double_t phiAOD = track->Phi();
      Double_t phiTLV = phiAOD;
      if (phiAOD > TMath::Pi()) {
	phiTLV -= 2 * TMath::Pi();
      }
      if (phiTLV < -TMath::Pi()) {
	phiTLV += 2 * TMath::Pi();
      }
      dynamic_cast<TH1F*>(fOutputList->FindObject("fEta_proton"))->Fill(track->Eta());
      dynamic_cast<TH1F*>(fOutputList->FindObject("fPhi_proton"))->Fill(phiTLV);

      //== TTree_proton
      fevtID_proton         = fHeader->GetEventIdAsLong();
      fzvertex_proton       = vertex->GetZ();
      fcentrality_proton    = centralityV0M;
      //=== Daughter track
      ftrkID_proton     = track->GetID();
      if(ftrkID_proton < 0) ftrkID_proton =-ftrkID_proton-1;
      fdcaxy_proton     = dcaxy;
      fdcaz_proton      = dcaz;
      fpx_proton        = track->Px();
      fpy_proton        = track->Py();
      fpz_proton        = track->Pz();
      fpt_proton        = track_pT;
      fpTPC_proton      = track->GetTPCmomentum();
      fnSigmaTPC_proton = fabs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton));
      feta_proton       = track->Eta();
      fphi_proton       = phiTLV;
      if(track->GetTOFsignal()!=99999){
	fnSigmaTOF_proton  = fabs(CalculateProtonSigma(track->Pt(),CalculateMassSquareTOF(*track),track->Charge()));
	fTrkMassTOF_proton = CalculateMassSquareTOF(*track)/TMath::Power(track->Charge(),2);
	dynamic_cast<TH2F*>(fOutputList->FindObject("fTrkMasspT_Proton")) 
	  ->Fill(track->Pt(),CalculateMassSquareTOF(*track)/TMath::Power(track->Charge(),2));
	dynamic_cast<TH1F*>(fOutputList->FindObject("fSigma_Proton"))
	  ->Fill(CalculateProtonSigma(track->Pt(),CalculateMassSquareTOF(*track),track->Charge()));
      }
      else {
	fnSigmaTOF_proton  = -999;
	fTrkMassTOF_proton = -999;
      }
      if(ITSisOK){
	fnSigmaITS_proton = CalculateSigmadEdxITS(*track);
      }
      else{
	fnSigmaITS_proton = -999;
      }

      fTree_proton->Fill();
    }
    
    if(track->Charge()<0) {
      Double_t phiAOD = track->Phi();
      Double_t phiTLV = phiAOD;
      if (phiAOD > TMath::Pi()) {
	phiTLV -= 2 * TMath::Pi();
      }
      if (phiTLV < -TMath::Pi()) {
	phiTLV += 2 * TMath::Pi();
      }
      dynamic_cast<TH1F*>(fOutputList->FindObject("fEta_proton"))->Fill(track->Eta());
      dynamic_cast<TH1F*>(fOutputList->FindObject("fPhi_proton"))->Fill(phiTLV);
      //== TTree_proton
      fevtID_antiproton         = fHeader->GetEventIdAsLong();
      fzvertex_antiproton       = vertex->GetZ();
      fcentrality_antiproton    = centralityV0M;
      //=== Daughter track
      ftrkID_antiproton     = track->GetID();
      if(ftrkID_antiproton < 0) ftrkID_antiproton =-ftrkID_antiproton-1;
      fdcaxy_antiproton     = dcaxy;
      fdcaz_antiproton      = dcaz;
      fpx_antiproton        = track->Px();
      fpy_antiproton        = track->Py();
      fpz_antiproton        = track->Pz();
      fpt_antiproton        = track_pT;
      fpTPC_antiproton      = track->GetTPCmomentum();
      fnSigmaTPC_antiproton = fabs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton));
      feta_antiproton       = track->Eta();
      fphi_antiproton       = phiTLV;
      if(track->GetTOFsignal()!=99999){
	fnSigmaTOF_antiproton = fabs(CalculateProtonSigma(track->Pt(),CalculateMassSquareTOF(*track),track->Charge()));
	fTrkMassTOF_antiproton = CalculateMassSquareTOF(*track)/TMath::Power(track->Charge(),2);
	dynamic_cast<TH2F*>(fOutputList->FindObject("fTrkMasspT_Proton")) 
	  ->Fill(track->Pt(),CalculateMassSquareTOF(*track)/TMath::Power(track->Charge(),2));
	dynamic_cast<TH1F*>(fOutputList->FindObject("fSigma_AntiProton"))
	  ->Fill(CalculateProtonSigma(track->Pt(),CalculateMassSquareTOF(*track),track->Charge()));
      }
      else {
	fnSigmaTOF_antiproton  = -999;
	fTrkMassTOF_antiproton = -999;
      }
      if(ITSisOK){
	fnSigmaITS_antiproton = CalculateSigmadEdxITS(*track);
      }
      else{
	fnSigmaITS_antiproton = -999;
      }

      fTree_antiproton->Fill();
    }

  }

  //============== cascade loop  
  for (int iCasc=0; iCasc< nCascades; iCasc++){
    AliAODcascade *casc= (AliAODcascade*)fAOD->GetCascade(iCasc);
    if(!casc) continue;
    
    AliAODTrack* pTrackXi   = dynamic_cast<AliAODTrack*>(casc->GetDaughter(0));                     // daughter track + 
    AliAODTrack* nTrackXi   = dynamic_cast<AliAODTrack*>(casc->GetDaughter(1));                     // daughter track -
    AliAODTrack* bachTrackXi= dynamic_cast<AliAODTrack*>(casc->GetDecayVertexXi()->GetDaughter(0)); // bachlor track
    if(!(IsQualityTrack(pTrackXi)))    continue;
    if(!(IsQualityTrack(nTrackXi)))    continue;
    if(!(IsQualityTrack(bachTrackXi))) continue;
    dcaxyp=0.,dcazp=0.,dcap=0.;  // daughter positive track
    dcaxyn=0.,dcazn=0.,dcan=0.;  // daughter negative track
    dcaxyb=0.,dcazb=0.,dcab=0.;  // bachelor track
    dDCAp[2]    ={0.};
    cDCAp[3]    ={0.};
    dDCAn[2]    ={0.};
    cDCAn[3]    ={0.};
    dDCAb[2]    ={0.};
    cDCAb[3]    ={0.};
    v0Vtx_casc[3]    ={0.};
    xivertex[3] ={0.};
    lambdacpa=0.,lambdatransradius=0.;
    xicpa=0.,xitransradius=0.;
    v0vtxx=0.,v0vtxy=0.,v0vtxz=0.,xivtxx=0.,xivtxy=0.,xivtxz=0.;
    
    if((pTrackXi->Charge()) < (nTrackXi->Charge())) continue;

    //========== Daughter track selection                                                                             
    pTrackXi->GetImpactParameters(dDCAp,cDCAp);
    dcaxyp          =dDCAp[0];
    dcazp           =dDCAp[1];
    dcap            =sqrt(dcaxyp*dcaxyp+dcazp*dcazp);
    nTrackXi->GetImpactParameters(dDCAn,cDCAn);
    dcaxyn          =dDCAn[0];
    dcazn           =dDCAn[1];
    dcan            =sqrt(dcaxyn*dcaxyn+dcazn*dcazn);
    bachTrackXi->GetImpactParameters(dDCAb,cDCAb);
    dcaxyb          =dDCAb[0];
    dcazb           =dDCAb[1];
    dcab            =sqrt(dcaxyb*dcaxyb+dcazb*dcazb);

    //cascade daughter select contents                                                                                                                          
    dynamic_cast<TH1F*>(fOutputList->FindObject("fcascdaupT_Proton"))->Fill(pTrackXi->Pt());
    dynamic_cast<TH1F*>(fOutputList->FindObject("fcascdaupT_Pion"))->Fill(nTrackXi->Pt());
    dynamic_cast<TH1F*>(fOutputList->FindObject("fcascdaudcab"))->Fill(bachTrackXi->Pt());

    if(fabs(pTrackXi   ->Eta()) > 0.8)       continue;
    if(fabs(nTrackXi   ->Eta()) > 0.8)       continue;
    if(fabs(bachTrackXi->Eta()) > 0.8)       continue;
    if(pTrackXi   ->GetTPCNcls() < 70) continue;
    if(nTrackXi   ->GetTPCNcls() < 70) continue;
    if(bachTrackXi->GetTPCNcls() < 70) continue;
    //if(pTrackXi   ->Pt() > 0.3)        continue;
    //if(nTrackXi   ->Pt() > 0.3)        continue;
    //if(bachTrackXi->Pt() > 0.3)        continue;
    if(dcap < 0.05) continue;
    if(dcan < 0.05) continue;
    if(dcab < 0.05) continue;
    if(dcap > 3) continue;
    if(dcan > 3) continue;
    if(dcab > 3) continue;

    dynamic_cast<TH1F*>(fOutputList->FindObject("fcascdauEta"))->Fill(pTrackXi->Eta());
    dynamic_cast<TH1F*>(fOutputList->FindObject("fcascdauEta"))->Fill(nTrackXi->Eta());
    dynamic_cast<TH1F*>(fOutputList->FindObject("fcascdauEtab"))->Fill(bachTrackXi->Eta());
    dynamic_cast<TH1F*>(fOutputList->FindObject("fcascdauTPC"))->Fill(pTrackXi->GetTPCNcls());
    dynamic_cast<TH1F*>(fOutputList->FindObject("fcascdauTPC"))->Fill(nTrackXi->GetTPCNcls());
    dynamic_cast<TH1F*>(fOutputList->FindObject("fcascdauTPCb"))->Fill(bachTrackXi->GetTPCNcls());
    dynamic_cast<TH1F*>(fOutputList->FindObject("fcascdauPt"))->Fill(pTrackXi->Pt());
    dynamic_cast<TH1F*>(fOutputList->FindObject("fcascdauPt"))->Fill(nTrackXi->Pt());
    dynamic_cast<TH1F*>(fOutputList->FindObject("fcascdauPtb"))->Fill(bachTrackXi->Pt());
    dynamic_cast<TH1F*>(fOutputList->FindObject("fcascdauPIDproton"))->Fill(fabs(fPIDResponse->NumberOfSigmasTPC(pTrackXi, AliPID::kProton)));
    dynamic_cast<TH1F*>(fOutputList->FindObject("fcascdauPIDproton"))->Fill(fabs(fPIDResponse->NumberOfSigmasTPC(nTrackXi, AliPID::kProton)));
    dynamic_cast<TH1F*>(fOutputList->FindObject("fcascdauPIDpion"))->Fill(fabs(fPIDResponse->NumberOfSigmasTPC(pTrackXi, AliPID::kPion)));
    dynamic_cast<TH1F*>(fOutputList->FindObject("fcascdauPIDpion"))->Fill(fabs(fPIDResponse->NumberOfSigmasTPC(nTrackXi, AliPID::kPion)));
    dynamic_cast<TH1F*>(fOutputList->FindObject("fcascbachPIDpion"))->Fill(fabs(fPIDResponse->NumberOfSigmasTPC(bachTrackXi, AliPID::kPion)));

    //========== PID
    isProton     =(fabs(fPIDResponse->NumberOfSigmasTPC(pTrackXi, AliPID::kProton)) < 3);  // proton
    isAntiProton =(fabs(fPIDResponse->NumberOfSigmasTPC(nTrackXi, AliPID::kProton)) < 3);  // proton-
    isPPion      =(fabs(fPIDResponse->NumberOfSigmasTPC(pTrackXi, AliPID::kPion)) < 3);    // pion+
    isNPion      =(fabs(fPIDResponse->NumberOfSigmasTPC(nTrackXi, AliPID::kPion)) < 3);    // pion-
    isBachpion   =(fabs(fPIDResponse->NumberOfSigmasTPC(bachTrackXi, AliPID::kPion)) < 3); // bacelor (pion)
    isBachkaon   =(fabs(fPIDResponse->NumberOfSigmasTPC(bachTrackXi, AliPID::kKaon)) < 3); // bacelor (kaon) 
    
    Bool_t isXi         =(isProton     && isNPion && bachTrackXi->Charge() < 0 && isBachpion);
    Bool_t isXibar      =(isAntiProton && isPPion && bachTrackXi->Charge() > 0 && isBachpion);
    Bool_t isOmega      =(isProton     && isNPion && bachTrackXi->Charge() < 0 && isBachkaon);
    Bool_t isOmegabar   =(isAntiProton && isPPion && bachTrackXi->Charge() > 0 && isBachkaon);

    //========== lambda selection
    v0vtxx               =casc->DecayVertexV0X();
    v0vtxy               =casc->DecayVertexV0Y();
    v0vtxz               =casc->DecayVertexV0Z();
    v0Vtx_casc[0]        =v0vtxx;
    v0Vtx_casc[1]        =v0vtxy;
    v0Vtx_casc[2]        =v0vtxz;
    lambdacpa            =LambdaCosPointingAngle(casc,v0Vtx_casc,vecTarget);
    lambdatransradius    =DecayLengthXY(v0Vtx_casc,vecTarget);

    //cascade V0 select contents
    dynamic_cast<TH1F*>(fOutputList->FindObject("fcascV0Transradius"))->Fill(lambdatransradius);
    dynamic_cast<TH1F*>(fOutputList->FindObject("fcascV0daugdca"))->Fill(casc->DcaV0Daughters());
    dynamic_cast<TH1F*>(fOutputList->FindObject("fcascV0dca"))->Fill(casc->DcaV0ToPrimVertex());
    dynamic_cast<TH1F*>(fOutputList->FindObject("fcascV0cpa"))->Fill(lambdacpa);

    if(casc->GetOnFlyStatus())           continue; // select offline
    if(lambdacpa < 0.97)                 continue;

    Bool_t CascV0isXi    = true;
    Bool_t CascV0isOmega = true;
    
    if(lambdatransradius < 1.4)          CascV0isXi = false;
    if(lambdatransradius > 200)          CascV0isXi = false;
    if(casc->DcaV0Daughters() > 1.5)     CascV0isXi = false;
    //if(casc->DcaV0ToPrimVertex() < 0.07) CascV0isXi = false;

    //if(lambdatransradius < 1.0)          CascV0isOmega = false;
    //if(lambdatransradius > 200)          CascV0isOmega = false;
    //if(casc->DcaV0Daughters() > 1.2)     CascV0isOmega = false;
    //if(casc->DcaV0ToPrimVertex() < 0.06) CascV0isOmega = false;

    if(CascV0isXi && isXi){
      dynamic_cast<TH1F*>(fOutputList->FindObject("hInvMassXidecayLambda"))->Fill(casc->MassLambda());
    }
    if(CascV0isXi && isXibar){
      dynamic_cast<TH1F*>(fOutputList->FindObject("hInvMassXidecayLambda"))->Fill(casc->MassAntiLambda());
    }

    //===== xi selection =====
    xivtxx        =casc->DecayVertexXiX();
    xivtxy        =casc->DecayVertexXiY();
    xivtxz        =casc->DecayVertexXiZ();
    xivertex[0]   =xivtxx;
    xivertex[1]   =xivtxy;
    xivertex[2]   =xivtxz;
    xicpa         =casc->CosPointingAngleXi(vecTarget[0],vecTarget[1],vecTarget[2]);
    xitransradius =xiDecayLengthXY(xivertex,vecTarget);

    //Xi select contents
    dynamic_cast<TH1F*>(fOutputList->FindObject("XiTranverseradius"))->Fill(xitransradius);
    dynamic_cast<TH1F*>(fOutputList->FindObject("dcaXiDaughters"))->Fill(casc->DcaXiDaughters());
    dynamic_cast<TH1F*>(fOutputList->FindObject("XiCPA"))->Fill(xicpa);
    //dynamic_cast<TH1F*>(fOutputList->FindObject("hInvMassOmega"))->Fill(InvMassOmega(casc));

    Bool_t CascisXi    = true;
    Bool_t CascisOmega = true;

    if(xicpa < 0.999)                 CascisXi = false;
    if(xitransradius < 0.8)           CascisXi = false;
    if(xitransradius > 200)           CascisXi = false;
    if(casc->DcaXiDaughters() > 0.1)  CascisXi = false;
    if((casc->MassOmega() > 1.667) && (casc->MassOmega() < 1.677)) CascisXi = false;

    if(xicpa < 0.995)                CascisOmega = false;
    if(xitransradius < 0.2)          CascisOmega = false;
    if(xitransradius > 200)          CascisOmega = false;
    if(casc->DcaXiDaughters() > 0.8) CascisOmega = false;
    if((casc->MassXi() > 1.317) && (casc->MassXi() < 1.327)) CascisOmega = false;
    //=========== lambda + pion- -> xi-
    //if(isXi && CascV0isXi && CascisXi){
    if(isOmega && CascV0isXi && CascisOmega){
      if(TMath::Abs(casc->MassLambda() - TDatabasePDG::Instance()->GetParticle(3122)->Mass()) > 0.00624465) continue;
      dynamic_cast<TH1F*>(fOutputList->FindObject("hInvMassXi"))->Fill(casc->MassOmega());
      if(TMath::Abs(casc->MassOmega() - TDatabasePDG::Instance()->GetParticle(3334)->Mass()) > 0.00907) continue;
      //Side band
      //if(TMath::Abs(casc->MassXi() - TDatabasePDG::Instance()->GetParticle(3312)->Mass()) < 0.005) continue;      
      //if(TMath::Abs(casc->MassXi() - TDatabasePDG::Instance()->GetParticle(3312)->Mass()) > 0.01)  continue;
      dynamic_cast<TH1F*>(fOutputList->FindObject("fPtcascade_xi"))->Fill(sqrt(casc->Pt2Xi()));
      dynamic_cast<TH1F*>(fOutputList->FindObject("fPtcascdecay_lambda"))->Fill(sqrt(casc->Pt2V0()));
      dynamic_cast<TH1F*>(fOutputList->FindObject("fPtBachelor"))->Fill(sqrt(pow(casc->MomBachX(),2)+pow(casc->MomBachY(),2)));
      dynamic_cast<TH1F*>(fOutputList->FindObject("fPhiXidecayLambda"))->Fill(casc->OpenAngleV0());

      //== TTree cascade 
      fpt_daughter1_cascade = pTrackXi->Pt();
      fpt_daughter2_cascade = nTrackXi->Pt();
      fpt_daughter3_cascade = bachTrackXi->Pt();
      //=== Lambda                                                                                               
      fcpa_lambda_cascade = lambdacpa;
      fpca_lambda_cascade = casc->DcaV0Daughters();
      fpt_lambda_cascade  = sqrt(casc->Pt2V0());
      //=== Omega                                                                                                   
      fcpa_omega_cascade     = xicpa;
      fpca_omega_cascade     = casc->DcaXiDaughters();
      fmass_omega_cascade    = casc->MassOmega(); 
      fpt_omega_cascade      = sqrt(casc->Pt2Xi()); 

      fTree_cascade->Fill();
    }
    //========== antilambda + pion+ -> xi+
    //if(isXibar && CascV0isXi && CascisXi){
    if(isOmegabar && CascV0isXi && CascisOmega){
      if(TMath::Abs(casc->MassAntiLambda() - TDatabasePDG::Instance()->GetParticle(3122)->Mass()) > 0.00624465) continue;
      dynamic_cast<TH1F*>(fOutputList->FindObject("hInvMassXi"))->Fill(casc->MassOmega());

      if(TMath::Abs(casc->MassOmega() - TDatabasePDG::Instance()->GetParticle(3334)->Mass()) > 0.00907) continue;
      //Side band
      //if(TMath::Abs(casc->MassXi() - TDatabasePDG::Instance()->GetParticle(3312)->Mass()) < 0.005) continue;      
      //if(TMath::Abs(casc->MassXi() - TDatabasePDG::Instance()->GetParticle(3312)->Mass()) > 0.01) continue; 
      dynamic_cast<TH1F*>(fOutputList->FindObject("fPtcascade_antixi"))->Fill(sqrt(casc->Pt2Xi()));
      dynamic_cast<TH1F*>(fOutputList->FindObject("fPtcascdecay_antilambda"))->Fill(sqrt(casc->Pt2V0()));
      dynamic_cast<TH1F*>(fOutputList->FindObject("fPtBachelor"))->Fill(sqrt(pow(casc->MomBachX(),2)+pow(casc->MomBachY(),2)));
      dynamic_cast<TH1F*>(fOutputList->FindObject("fPhiXidecayLambda"))->Fill(casc->OpenAngleV0());

      //== TTree cascade                                                                                        
      fpt_daughter1_cascade = nTrackXi->Pt();
      fpt_daughter2_cascade = pTrackXi->Pt();
      fpt_daughter3_cascade = bachTrackXi->Pt();
      //=== Lambda                                                                                               
      fcpa_lambda_cascade = lambdacpa;
      fpca_lambda_cascade = casc->DcaV0Daughters();
      fpt_lambda_cascade  = sqrt(casc->Pt2V0());
      //=== Omega                                                                                                   
      fcpa_omega_cascade     = xicpa;
      fpca_omega_cascade     = casc->DcaXiDaughters();
      fmass_omega_cascade    = casc->MassOmega(); 
      fpt_omega_cascade      = sqrt(casc->Pt2Xi()); 

      fTree_cascade->Fill();

    }
  }
  //============== cascade loop end


  // Post output data
  PostData(1,fOutputList);
  PostData(2,fTree_omegam);
  PostData(3,fTree_omegap);
  PostData(4,fTree_sidebandm);
  PostData(5,fTree_sidebandp);
  PostData(6,fTree_cascade);
  PostData(7,fTree_proton);
  PostData(8,fTree_antiproton);
  //cout<<"+++++++++++++++++++++++++++++++++++++++++++++"<<endl;
  //cout<<"+++++++++++++++ Analysis finish +++++++++++++"<<endl;
  //cout<<"+++++++++++++++++++++++++++++++++++++++++++++"<<endl;

}


//________________________________________________________________________
Bool_t AliAnalysisTaskpOmegaDibaryon::EventSelection(AliAODEvent* data)
{
  const AliVVertex *vtx   =data->GetPrimaryVertex();
  const AliVVertex *vtxSPD=data->GetPrimaryVertexSPD();  
  
  Double_t xvtx=0.,yvtx=0.,zvtx=0.;
  Double_t xvtxSPD=0.,yvtxSPD=0.,zvtxSPD=0.;
  Int_t    ncont=0,ncontSPD=0;
  Double_t vdisp=0.;
  // Event vertex
  xvtx =vtx->GetX();
  yvtx =vtx->GetY();
  zvtx =vtx->GetZ();
  ncont=vtx->GetNContributors();

  dynamic_cast<TH1F*>(fOutputList->FindObject("fEvtVtxX"))->Fill(xvtx);
  dynamic_cast<TH1F*>(fOutputList->FindObject("fEvtVtxY"))->Fill(yvtx);
  dynamic_cast<TH1F*>(fOutputList->FindObject("fEvtVtxZ"))->Fill(zvtx);
  dynamic_cast<TH1F*>(fOutputList->FindObject("fEvtVtxTrk"))->Fill(ncont);
  // SPD vertex
  xvtxSPD =vtxSPD->GetX();
  yvtxSPD =vtxSPD->GetY();
  zvtxSPD =vtxSPD->GetZ();
  ncontSPD=vtxSPD->GetNContributors();
  vdisp=fabs(zvtx-zvtxSPD);

  dynamic_cast<TH1F*>(fOutputList->FindObject("fSPDVtxZ"))->Fill(zvtxSPD);
  dynamic_cast<TH1F*>(fOutputList->FindObject("fSPDVtxTrk"))->Fill(ncontSPD);
  dynamic_cast<TH2F*>(fOutputList->FindObject("fSPDVtxCor"))->Fill(zvtx,zvtxSPD);
  dynamic_cast<TH1F*>(fOutputList->FindObject("fSPDVtxDisp"))->Fill(vdisp);

  Bool_t zvtx_cut=kFALSE;
  //zvtx_cut=(fabs(zvtx)<10. && ncont>0 && ncontSPD>0 && vdisp<0.5);
  zvtx_cut=(abs(zvtx) < 10. && ncont > 0);
  return zvtx_cut;
}
//________________________________________________________________________

//invariant mass Lambda cascade class
Double_t AliAnalysisTaskpOmegaDibaryon::InvMassLambda(AliAODcascade *casc)
{  
  Double_t pPx=0.,pPy=0.,pPz=0.,pE=0.; //proton
  Double_t nPx=0.,nPy=0.,nPz=0.,nE=0.; //pion-
  Float_t  energysum=0.,psum2=0.;
  Double_t invMass =0.;  
  pPx=casc->MomPosX();
  pPy=casc->MomPosY();
  pPz=casc->MomPosZ();
  nPx=casc->MomNegX();
  nPy=casc->MomNegY();
  nPz=casc->MomNegZ();
  pE =sqrt(0.938*0.938+pPx*pPx+pPy*pPy+pPz*pPz); //proton
  nE =sqrt(0.14*0.14+nPx*nPx+nPy*nPy+nPz*nPz);   //pion-
  energysum =pE+nE;
  psum2     =(pPx+nPx)*(pPx+nPx)+(pPy+nPy)*(pPy+nPy)+(pPz+nPz)*(pPz+nPz);
  invMass   =((energysum*energysum-psum2)>0)?sqrt((energysum*energysum)-psum2):0; 
  return invMass;
}
//invariant mass anti lambda
Double_t AliAnalysisTaskpOmegaDibaryon::InvMassAntiLambda(AliAODcascade *casc)
{  
  Double_t nPx=0.,nPy=0.,nPz=0.,nE=0.; //anti proton
  Double_t pPx=0.,pPy=0.,pPz=0.,pE=0.; //pion+
  Float_t  energysum=0.,psum2=0.;
  Double_t invMass =0.;  
  nPx=casc->MomNegX();
  nPy=casc->MomNegY();
  nPz=casc->MomNegZ();
  pPx=casc->MomPosX();
  pPy=casc->MomPosY();
  pPz=casc->MomPosZ();
  nE =sqrt(0.938*0.938+nPx*nPx+nPy*nPy+nPz*nPz); //anti proton
  pE =sqrt(0.14*0.14+pPx*pPx+pPy*pPy+pPz*pPz);   //pion+
  energysum =nE+pE;
  psum2     =(nPx+pPx)*(nPx+pPx)+(nPy+pPy)*(nPy+pPy)+(nPz+pPz)*(nPz+pPz);
  invMass   =((energysum*energysum-psum2)>0)?sqrt((energysum*energysum)-psum2):0; 
  return invMass;
}
//invariant mass Xi-
Double_t AliAnalysisTaskpOmegaDibaryon::InvMassXi(AliAODcascade *casc)
{  
  Double_t lPx=0.,lPy=0.,lPz=0.,lE=0.; //lambda
  Double_t bPx=0.,bPy=0.,bPz=0.,bE=0.; //bachelor pion-
  Float_t  energysum=0.,psum2=0.;
  Double_t invMass =0.;  
  lPx=casc->MomV0X();
  lPy=casc->MomV0Y();
  lPz=casc->MomV0Z();
  bPx=casc->MomBachX();
  bPy=casc->MomBachY();
  bPz=casc->MomBachZ();
  lE =sqrt(1.1157*1.1157+lPx*lPx+lPy*lPy+lPz*lPz); //lambda 1.115
  bE =sqrt(0.14*0.14+bPx*bPx+bPy*bPy+bPz*bPz);   //pion-
  energysum =lE+bE;
  psum2     =(lPx+bPx)*(lPx+bPx)+(lPy+bPy)*(lPy+bPy)+(lPz+bPz)*(lPz+bPz);
  invMass   =((energysum*energysum-psum2)>0)?sqrt((energysum*energysum)-psum2):0; 
  return invMass;
}
//omega rejection
Double_t AliAnalysisTaskpOmegaDibaryon::InvMassOmega(AliAODcascade *casc)
{  
  Double_t lPx=0.,lPy=0.,lPz=0.,lE=0.; //lambda
  Double_t bPx=0.,bPy=0.,bPz=0.,bE=0.; //bachelor kaon
  Float_t  energysum=0.,psum2=0.;
  Double_t invMass =0.;  
  lPx=casc->MomV0X();
  lPy=casc->MomV0Y();
  lPz=casc->MomV0Z();
  bPx=casc->MomBachX();
  bPy=casc->MomBachY();
  bPz=casc->MomBachZ();
  lE =sqrt(1.1157*1.1157+lPx*lPx+lPy*lPy+lPz*lPz); //lambda 1.116
  bE =sqrt(0.494*0.494+bPx*bPx+bPy*bPy+bPz*bPz); //kaon-
  energysum =lE+bE;
  psum2     =(lPx+bPx)*(lPx+bPx)+(lPy+bPy)*(lPy+bPy)+(lPz+bPz)*(lPz+bPz);
  invMass   =((energysum*energysum-psum2)>0)?sqrt((energysum*energysum)-psum2):0; 
  return invMass;
}
//lambda cosine pointing angle
Double_t AliAnalysisTaskpOmegaDibaryon::LambdaCosPointingAngle(AliAODcascade *casc,const Double_t *DecayVtx,
							    const Float_t *point) const 
{  
  TVector3 v0Mom(casc->MomV0X(),casc->MomV0Y(),casc->MomV0Z());
  TVector3 fline(DecayVtx[0] - point[0], DecayVtx[1] - point[1],
                 DecayVtx[2] - point[2]);
  Double_t ptot2 = v0Mom.Mag2() * fline.Mag2();
  if (ptot2 <= 0) {
    return 0.0;
  }
  else {
    Double_t cos = v0Mom.Dot(fline) / TMath::Sqrt(ptot2);
    if (cos > 1.0)
      cos = 1.0;
    if (cos < -1.0)
      cos = -1.0;
    return cos;
  }
}
//transverse radius of the lambda decay vertex
Double_t AliAnalysisTaskpOmegaDibaryon::DecayLengthXY(const Double_t *DecayVtx,const Float_t *point) const 
{ return TMath::Sqrt( (point[0] - DecayVtx[0]) * (point[0] - DecayVtx[0])
		      + (point[1] - DecayVtx[1]) * (point[1] - DecayVtx[1]));
}
//transverse radius of the xi decay vertex
Double_t AliAnalysisTaskpOmegaDibaryon::xiDecayLengthXY(const Double_t *xiDecayVtx,const Float_t *point) const 
{ return TMath::Sqrt( (point[0] - xiDecayVtx[0]) * (point[0] - xiDecayVtx[0])
		      + (point[1] - xiDecayVtx[1]) * (point[1] - xiDecayVtx[1]));
}
//invariant mass Lambda v0 class
Double_t AliAnalysisTaskpOmegaDibaryon::InvMasslambda(AliAODv0 *v0)
{  
  Double_t pPx=0.,pPy=0.,pPz=0.,pE=0.; //proton
  Double_t nPx=0.,nPy=0.,nPz=0.,nE=0.; //pion
  Float_t  energysum=0.,psum2=0.;
  Double_t invMass =0.;
  pPx=v0->MomPosX();
  pPy=v0->MomPosY();
  pPz=v0->MomPosZ();
  nPx=v0->MomNegX();
  nPy=v0->MomNegY();
  nPz=v0->MomNegZ();
  pE =sqrt(0.938*0.938+pPx*pPx+pPy*pPy+pPz*pPz); //proton
  nE =sqrt(0.14*0.14+nPx*nPx+nPy*nPy+nPz*nPz);   //pion-
  energysum =pE+nE;
  psum2     =(pPx+nPx)*(pPx+nPx)+(pPy+nPy)*(pPy+nPy)+(pPz+nPz)*(pPz+nPz);
  invMass   =((energysum*energysum-psum2)>0)?sqrt((energysum*energysum)-psum2):0; 
  return invMass;
}

//K0 rejection                                                                                         
Double_t AliAnalysisTaskpOmegaDibaryon::InvMassK0(AliAODv0 *v0)
{
  Double_t pPx=0.,pPy=0.,pPz=0.,pE=0.; //pion+                                                         
  Double_t nPx=0.,nPy=0.,nPz=0.,nE=0.; //pion                                                          
  Float_t  energysum=0.,psum2=0.;
  Double_t invMass =0.;
  pPx=v0->MomPosX();
  pPy=v0->MomPosY();
  pPz=v0->MomPosZ();
  nPx=v0->MomNegX();
  nPy=v0->MomNegY();
  nPz=v0->MomNegZ();
  pE =sqrt(0.14*0.14+pPx*pPx+pPy*pPy+pPz*pPz); //pion+                                                 
  nE =sqrt(0.14*0.14+nPx*nPx+nPy*nPy+nPz*nPz); //pion-                                                 
  energysum =pE+nE;
  psum2     =(pPx+nPx)*(pPx+nPx)+(pPy+nPy)*(pPy+nPy)+(pPz+nPz)*(pPz+nPz);
  invMass   =((energysum*energysum-psum2)>0)?sqrt((energysum*energysum)-psum2):0;
  return invMass;
}

//K0 rejection                                                                                         
Double_t AliAnalysisTaskpOmegaDibaryon::InvMassSelf(AliExternalTrackParam *track1, AliExternalTrackParam *track2, Double_t daughtermass1, Double_t daughtermass2)
{
  Double_t pPx=0.,pPy=0.,pPz=0.,pE=0.; //pion+                                                         
  Double_t nPx=0.,nPy=0.,nPz=0.,nE=0.; //pion                                                          
  Float_t  energysum=0.,psum2=0.;
  Double_t invMass =0.;
  pPx=track1->Px();
  pPy=track1->Py();
  pPz=track1->Pz();
  nPx=track2->Px();
  nPy=track2->Py();
  nPz=track2->Pz();
  pE =sqrt(TMath::Power(daughtermass1,2)+pPx*pPx+pPy*pPy+pPz*pPz); //pion+                                                 
  nE =sqrt(TMath::Power(daughtermass2,2)+nPx*nPx+nPy*nPy+nPz*nPz); //pion-                                                 
  energysum =pE+nE;
  psum2     =(pPx+nPx)*(pPx+nPx)+(pPy+nPy)*(pPy+nPy)+(pPz+nPz)*(pPz+nPz);
  invMass   =((energysum*energysum-psum2)>0)?sqrt((energysum*energysum)-psum2):0;
  return invMass;
}

//invariant mass antiLambda v0 class
Double_t AliAnalysisTaskpOmegaDibaryon::InvMassAntilambda(AliAODv0 *v0)
{
  Double_t pPx=0.,pPy=0.,pPz=0.,pE=0.; //pion+                                                                                               
  Double_t nPx=0.,nPy=0.,nPz=0.,nE=0.; //antiproton                                                                                            
  Float_t  energysum=0.,psum2=0.;
  Double_t invMass =0.;
  pPx=v0->MomPosX();
  pPy=v0->MomPosY();
  pPz=v0->MomPosZ();
  nPx=v0->MomNegX();
  nPy=v0->MomNegY();
  nPz=v0->MomNegZ();
  pE =sqrt(0.14*0.14+pPx*pPx+pPy*pPy+pPz*pPz);     //pion+                                                                                    
  nE =sqrt(0.938*0.938+nPx*nPx+nPy*nPy+nPz*nPz);   //antiproton                                                                              
  energysum =pE+nE;
  psum2     =(pPx+nPx)*(pPx+nPx)+(pPy+nPy)*(pPy+nPy)+(pPz+nPz)*(pPz+nPz);
  invMass   =((energysum*energysum-psum2)>0)?sqrt((energysum*energysum)-psum2):0;
  return invMass;
}
Double_t AliAnalysisTaskpOmegaDibaryon::CosPointingAngle(AliAODv0 *v0,const Double_t *DecayVtx,
						      const Float_t *point) const {  
  /// Cosine of pointing angle in space assuming it is produced at "point"
  TVector3 v0Mom(v0->MomV0X(),v0->MomV0Y(),v0->MomV0Z());
  TVector3 fline(DecayVtx[0] - point[0], DecayVtx[1] - point[1],
                 DecayVtx[2] - point[2]);
  Double_t ptot2 = v0Mom.Mag2() * fline.Mag2();
  if (ptot2 <= 0) {
    return 0.0;
  }
  else {
    Double_t cos = v0Mom.Dot(fline) / TMath::Sqrt(ptot2);
    if (cos > 1.0)
      cos = 1.0;
    if (cos < -1.0)
      cos = -1.0;
    return cos;
  }
}

Double_t AliAnalysisTaskpOmegaDibaryon::OpenAngle(Double_t px1,Double_t py1,Double_t pz1,
					       Double_t px2,Double_t py2,Double_t pz2){
  Double_t lScalPtot1Ptot2=0.,lPtot1xPtot2=0.;

  lScalPtot1Ptot2 = (px1*px2)+(py1*py2)+(pz1*pz2);
  lPtot1xPtot2 = (sqrt(pow(px1,2)+pow(py1,2)+pow(pz1,2)))*(sqrt(pow(px2,2)+pow(py2,2)+pow(pz2,2)));
  
  return acos(lScalPtot1Ptot2/lPtot1xPtot2);
}

Double_t AliAnalysisTaskpOmegaDibaryon::InvariantMass(Double_t px1,Double_t py1,Double_t pz1,
						   Double_t px2,Double_t py2,Double_t pz2,Double_t energysum){

  Double_t psum2=0.,pt=0.,invMass=0.;

  psum2     =pow((px1+px2),2)+pow((py1+py2),2)+pow((pz1+pz2),2);
  pt        =sqrt(pow(px1+px2,2)+pow(py1+py2,2));
  invMass   =sqrt((energysum*energysum)-psum2);

  return invMass;

}

//Bool_t AliAnalysisTaskpOmegaDibaryon::ExtractQnVector()
Double_t AliAnalysisTaskpOmegaDibaryon::ExtractQnVector()
{
  Int_t fHarmonics =2;
  Int_t fQnDetectorMain =0;

  //TString fTPCEPName[0] = "TPC";
  //TString fTPCEPName[1] = "TPCNegEta";
  //TString fTPCEPName[2] = "TPCPosEta";
  TString fTPCEPName[3] = {"TPC","TPCNegEta","TPCPosEta"};

  //fV0EPName[0] = "VZERO";
  //fV0EPName[1] = "VZEROA";
  //fV0EPName[2] = "VZEROC";

  TString fV0EPName[3] = {"VZERO","VZEROA","VZEROC"};
  TString fQNormalization = "QoverM";

  if(fHarmonics < 0){
    AliError(Form("Qn Flow vector correction flag is ON, but fHarmonics is not set. (it is %d now).",fHarmonics));
    return kFALSE;
  }

  TList* qnlist = fFlowQnVectorMgr->GetQnVectorList();

  const AliQnCorrectionsQnVector *QnVectorTPCDet[3];
  Double_t TPCEP[3] = {};
  for(Int_t i=0;i<3;i++){
    QnVectorTPCDet[i] = GetQnVectorFromList(qnlist,Form("%s%s",fTPCEPName[i].Data(),fQNormalization.Data()),"latest","plain");
    if(!QnVectorTPCDet[i]){
      AliInfo("Event is rejected because event plane is not found or bad event plane quality in TPC.");
      Printf("Event is rejected because event plane is not found or bad event plane quality in TPC.");
      return kFALSE;//Qn vector correction does not exist or bad quality.
    }
    TPCEP[i] = QnVectorTPCDet[i]->EventPlane(fHarmonics);
    if(TPCEP[i] < 0) TPCEP[i] += 2./(Double_t) fHarmonics * TMath::Pi();
    AliInfo(Form("harmonics %d | TPC sub detector name %s%s : event plane = %f (rad).",fHarmonics,fTPCEPName[i].Data(),fQNormalization.Data(),TPCEP[i]));
  }

  const AliQnCorrectionsQnVector *QnVectorV0Det[3];
  Double_t V0EP[3]  = {};
  for(Int_t i=0;i<3;i++){
    QnVectorV0Det[i]  = GetQnVectorFromList(qnlist,Form("%s%s",fV0EPName[i].Data(),fQNormalization.Data()),"latest","raw");
    if(!QnVectorV0Det[i]){
      AliInfo("Event is rejected because event plane is not found or bad event plane quality in VZERO.");
      Printf("Event is rejected because event plane is not found or bad event plane quality in VZERO.");
      return kFALSE;//Qn vector correction does not exist or bad quality.
    }
    V0EP[i] = QnVectorV0Det[i]->EventPlane(fHarmonics);
    if(V0EP[i] < 0)  V0EP[i]  += 2./(Double_t) fHarmonics * TMath::Pi();
    AliInfo(Form("harmonics %d | V0  sub detector name %s%s : event plane = %f (rad).",fHarmonics,fV0EPName[i].Data(),fQNormalization.Data(),V0EP[i]));
    //Printf("harmonics %d | V0  sub detector name %s%s : event plane = %f (rad).",fHarmonics,fV0EPName[i].Data(),fQNormalization.Data(),V0EP[i]);
  }

  //0 < event plane < 2*pi/fHarmonics.
  Double_t EP2 = -999; 
  Double_t EP3 = -999; 

  Double_t Q1[2] = {};//for Main
  Double_t Q2[2] = {};//for Sub1
  Double_t Q3[2] = {};//for Sub2

  Float_t fEventPlane = 999;

  if(fQnDetectorMain == AliAnalysisTaskpOmegaDibaryon::kFullTPC){
    Q1[0] = QnVectorTPCDet[0]->Qx(fHarmonics);//FullTPC
    Q1[1] = QnVectorTPCDet[0]->Qy(fHarmonics);//FullTPC
    Q2[0] = QnVectorV0Det[1]->Qx(fHarmonics);//V0A
    Q2[1] = QnVectorV0Det[1]->Qy(fHarmonics);//V0A
    Q3[0] = QnVectorV0Det[2]->Qx(fHarmonics);//V0C
    Q3[1] = QnVectorV0Det[2]->Qy(fHarmonics);//V0C
    fEventPlane = TPCEP[0];
    EP2 = V0EP[1];
    EP3 = V0EP[2];
  }
  else if(fQnDetectorMain == AliAnalysisTaskpOmegaDibaryon::kTPCNegEta){
    Q1[0] = QnVectorTPCDet[1]->Qx(fHarmonics);//TPCNegEta
    Q1[1] = QnVectorTPCDet[1]->Qy(fHarmonics);//TPCNegEta
    Q2[0] = QnVectorV0Det[1]->Qx(fHarmonics);//V0A
    Q2[1] = QnVectorV0Det[1]->Qy(fHarmonics);//V0A
    Q3[0] = QnVectorV0Det[2]->Qx(fHarmonics);//V0C
    Q3[1] = QnVectorV0Det[2]->Qy(fHarmonics);//V0C
    fEventPlane = TPCEP[2];
    EP2 = V0EP[1];
    EP3 = V0EP[2];
  }
  else if(fQnDetectorMain == AliAnalysisTaskpOmegaDibaryon::kTPCPosEta){
    Q1[0] = QnVectorTPCDet[2]->Qx(fHarmonics);//TPCPosEta
    Q1[1] = QnVectorTPCDet[2]->Qy(fHarmonics);//TPCPosEta
    Q2[0] = QnVectorV0Det[1]->Qx(fHarmonics);//V0A
    Q2[1] = QnVectorV0Det[1]->Qy(fHarmonics);//V0A
    Q3[0] = QnVectorV0Det[2]->Qx(fHarmonics);//V0C
    Q3[1] = QnVectorV0Det[2]->Qy(fHarmonics);//V0C
    fEventPlane = TPCEP[1];
    EP2 = V0EP[1];
    EP3 = V0EP[2];
  }
  else if(fQnDetectorMain == AliAnalysisTaskpOmegaDibaryon::kFullV0){
    Q1[0] = QnVectorV0Det[0]->Qx(fHarmonics);//FullV0
    Q1[1] = QnVectorV0Det[0]->Qy(fHarmonics);//FullV0
    Q2[0] = QnVectorTPCDet[1]->Qx(fHarmonics);//TPCNegEta
    Q2[1] = QnVectorTPCDet[1]->Qy(fHarmonics);//TPCNegEta
    Q3[0] = QnVectorTPCDet[2]->Qx(fHarmonics);//TPCPosEta
    Q3[1] = QnVectorTPCDet[2]->Qy(fHarmonics);//TPCPosEta
    fEventPlane = V0EP[0];
    EP2 = TPCEP[1];
    EP3 = TPCEP[2];
  }
  else if(fQnDetectorMain == AliAnalysisTaskpOmegaDibaryon::kV0A){
    Q1[0] = QnVectorV0Det[1]->Qx(fHarmonics);//V0A
    Q1[1] = QnVectorV0Det[1]->Qy(fHarmonics);//V0A
    Q2[0] = QnVectorV0Det[2]->Qx(fHarmonics);//V0C
    Q2[1] = QnVectorV0Det[2]->Qy(fHarmonics);//V0C
    Q3[0] = QnVectorTPCDet[0]->Qx(fHarmonics);//full acceptance of TPC
    Q3[1] = QnVectorTPCDet[0]->Qy(fHarmonics);//full acceptance of TPC
    fEventPlane = V0EP[1];
    EP2 = V0EP[2];
    EP3 = TPCEP[0];
  }
  else if(fQnDetectorMain == AliAnalysisTaskpOmegaDibaryon::kV0C){
    Q1[0] = QnVectorV0Det[2]->Qx(fHarmonics);//V0C
    Q1[1] = QnVectorV0Det[2]->Qy(fHarmonics);//V0C
    Q2[0] = QnVectorV0Det[1]->Qx(fHarmonics);//V0A
    Q2[1] = QnVectorV0Det[1]->Qy(fHarmonics);//V0A
    Q3[0] = QnVectorTPCDet[0]->Qx(fHarmonics);//full acceptance of TPC
    Q3[1] = QnVectorTPCDet[0]->Qy(fHarmonics);//full acceptance of TPC
    fEventPlane = V0EP[2];
    EP2 = V0EP[1];
    EP3 = TPCEP[0];
  }

  //fQVector1.Set(Q1[0],Q1[1]);
  TVector2 fQVector1(Q1[0],Q1[1]);
  TVector2 QVector2(Q2[0],Q2[1]);
  TVector2 QVector3(Q3[0],Q3[1]);
  Double_t sp12 = fQVector1 *  QVector2;//scalar product between Q1 vector and Q2 vector
  Double_t sp23 =  QVector2 *  QVector3;//scalar product between Q2 vector and Q3 vector
  Double_t sp31 =  QVector3 * fQVector1;//scalar product between Q3 vector and Q1 vector

  //AliInfo(Form("Q1x = %e , Q1y = %e , Q2x = %e , Q2y = %e , Q3x = %e , Q3y = %e ,  SP12 = %e ,  SP23 = %e ,  SP31 = %e",Q1[0],Q1[1],Q2[0],Q2[1],Q3[0],Q3[1],sp12,sp23,sp31));
  //Printf("Q1x = %e , Q1y = %e , Q2x = %e , Q2y = %e , Q3x = %e , Q3y = %e ,  SP12 = %e ,  SP23 = %e ,  SP31 = %e",Q1[0],Q1[1],Q2[0],Q2[1],Q3[0],Q3[1],sp12,sp23,sp31);
  
  //dynamic_cast<TH1F*>(fOutputList->FindObject("hEventPlane"))->Fill(EP2);
  dynamic_cast<TH1F*>(fOutputList->FindObject("hEventPlane"))->Fill(fEventPlane);
  //cout<<"========== EventPlane!! =============="<<endl;
  //cout<<fEventPlane<<endl;

  //  const Double_t delta = 2. * TMath::Pi() / Double_t(fHarmonics) / 12.;
  //  fEPBin = (Int_t)((fEventPlane) / delta);//it should be 0-11.
  //  if(fEPBin < 0)  fEPBin =  0;//protection to avoid fEPBin = -1.
  //  if(fEPBin > 11) fEPBin = 11;//protection to avoid fEPBin = 12.

  //return kTRUE;
  return fEventPlane;
}
const AliQnCorrectionsQnVector *AliAnalysisTaskpOmegaDibaryon::GetQnVectorFromList(const TList *list, const char* subdetector, const char *expcorr, const char *altcorr)
{
  AliQnCorrectionsQnVector *theQnVector = NULL;

  TList *pQvecList = dynamic_cast<TList*> (list->FindObject(subdetector));
  if(pQvecList != NULL){//sub detector is found
    if(TString(expcorr).EqualTo("latest")) theQnVector = (AliQnCorrectionsQnVector*) pQvecList->First();
    else                                   theQnVector = (AliQnCorrectionsQnVector*) pQvecList->FindObject(expcorr);
    //theQnVector = (AliQnCorrectionsQnVector*) pQvecList->FindObject(expcorr);

    if(theQnVector == NULL || !(theQnVector->IsGoodQuality()) || !(theQnVector->GetN() != 0)){ //the Qn vector for the expected correction is not found
      AliInfo(Form("expected correction (%s) is not found. use %s as an alternative step in %s.",expcorr,altcorr,subdetector));
      if(TString(altcorr).EqualTo("latest")) theQnVector = (AliQnCorrectionsQnVector*) pQvecList->First();
      else                                   theQnVector = (AliQnCorrectionsQnVector*) pQvecList->FindObject(altcorr);
      //theQnVector = (AliQnCorrectionsQnVector*) pQvecList->FindObject(altcorr);
    }

    //check the Qn vector quality
    if(!(theQnVector->IsGoodQuality()) || !(theQnVector->GetN() != 0)) theQnVector = NULL; //bad quality, discarded

  }
  return theQnVector;

}

AliESDVertex *AliAnalysisTaskpOmegaDibaryon::AODToESDVertex(const AliAODVertex &aodVert){ 
  Double_t covmatrix[6];  
  Double_t position[3]; 
  Double_t chi2; 
  Int_t nContributors;
  aodVert.GetCovMatrix(covmatrix); 
  aodVert.GetXYZ(position);
  chi2 = aodVert.GetChi2(); 
  nContributors = aodVert.GetNDF(); 
  nContributors = (nContributors + 3) / 2;
  AliESDVertex *esdVert = new AliESDVertex(position, covmatrix, chi2, nContributors, "primaryVertex");
  return esdVert;
}

KFVertex AliAnalysisTaskpOmegaDibaryon::CreateKFVertex(const AliVVertex& vertex) {
 
  /// GetTrack parameters
  Double_t param[6];
  Double_t cov[6];
 
  vertex.GetXYZ(param);
  vertex.GetCovarianceMatrix(cov);
 
  KFPVertex kfpVtx;
  /// Set the values
  Float_t paramF[3] = { (Float_t)param[0],(Float_t)param[1],(Float_t)param[2] };
  kfpVtx.SetXYZ(paramF);
  Float_t covF[6] = { (Float_t)cov[0],(Float_t)cov[1],(Float_t)cov[2],
		      (Float_t)cov[3],(Float_t)cov[4],(Float_t)cov[5] };
  kfpVtx.SetCovarianceMatrix(covF);
  KFVertex KFVtx(kfpVtx);
  return KFVtx;
}

KFParticle AliAnalysisTaskpOmegaDibaryon::CreateKFParticle(AliExternalTrackParam& track, Double_t Mass, Int_t Charge) {

  Double_t fP[6];
  track.GetXYZ(fP);
  track.PxPyPz(fP + 3);
  Int_t fQ = track.Charge() * TMath::Abs(Charge);
  fP[3] *= TMath::Abs(Charge);
  fP[4] *= TMath::Abs(Charge);
  fP[5] *= TMath::Abs(Charge);

  Double_t pt = 1. / TMath::Abs(track.GetParameter()[4]) * TMath::Abs(Charge);
  Double_t cs = TMath::Cos(track.GetAlpha()), sn = TMath::Sin(track.GetAlpha());
  Double_t r = TMath::Sqrt((1. - track.GetParameter()[2]) * (1. + track.GetParameter()[2]));

  Double_t m00 = -sn, m10 = cs;
  Double_t m23 = -pt * (sn + track.GetParameter()[2] * cs / r), m43 = -pt * pt * (r * cs - track.GetParameter()[2] * sn);
  Double_t m24 = pt * (cs - track.GetParameter()[2] * sn / r), m44 = -pt * pt * (r * sn + track.GetParameter()[2] * cs);
  Double_t m35 = pt, m45 = -pt * pt * track.GetParameter()[3];

  m43 *= track.GetSign();
  m44 *= track.GetSign();
  m45 *= track.GetSign();

  const Double_t* cTr = track.GetCovariance();
  Double_t fC[21];
  fC[0] = cTr[0] * m00 * m00;
  fC[1] = cTr[0] * m00 * m10;
  fC[2] = cTr[0] * m10 * m10;
  fC[3] = cTr[1] * m00;
  fC[4] = cTr[1] * m10;
  fC[5] = cTr[2];
  fC[6] = m00 * (cTr[3] * m23 + cTr[10] * m43);
  fC[7] = m10 * (cTr[3] * m23 + cTr[10] * m43);
  fC[8] = cTr[4] * m23 + cTr[11] * m43;
  fC[9] = m23 * (cTr[5] * m23 + cTr[12] * m43) + m43 * (cTr[12] * m23 + cTr[14] * m43);
  fC[10] = m00 * (cTr[3] * m24 + cTr[10] * m44);
  fC[11] = m10 * (cTr[3] * m24 + cTr[10] * m44);
  fC[12] = cTr[4] * m24 + cTr[11] * m44;
  fC[13] = m23 * (cTr[5] * m24 + cTr[12] * m44) + m43 * (cTr[12] * m24 + cTr[14] * m44);
  fC[14] = m24 * (cTr[5] * m24 + cTr[12] * m44) + m44 * (cTr[12] * m24 + cTr[14] * m44);
  fC[15] = m00 * (cTr[6] * m35 + cTr[10] * m45);
  fC[16] = m10 * (cTr[6] * m35 + cTr[10] * m45);
  fC[17] = cTr[7] * m35 + cTr[11] * m45;
  fC[18] = m23 * (cTr[8] * m35 + cTr[12] * m45) + m43 * (cTr[13] * m35 + cTr[14] * m45);
  fC[19] = m24 * (cTr[8] * m35 + cTr[12] * m45) + m44 * (cTr[13] * m35 + cTr[14] * m45);
  fC[20] = m35 * (cTr[9] * m35 + cTr[13] * m45) + m45 * (cTr[13] * m35 + cTr[14] * m45);

  KFParticle part;
  part.Create(fP, fC, fQ, Mass);
  return part;
}

//Track quality cut                                                                                                     
Bool_t AliAnalysisTaskpOmegaDibaryon::IsQualityTrack(AliAODTrack *track)
{
  if (!track->IsOn(AliAODTrack::kTPCrefit)) return kFALSE;
  if (track->GetTPCncls() < 70) return kFALSE;
  if (track->GetTPCchi2perCluster() > 2.5) return kFALSE;
  Float_t nCrossedRowsTPC = track->GetTPCClusterInfo(2,1);
  if (nCrossedRowsTPC < 70) return kFALSE;
  Int_t findable=track->GetTPCNclsF();
  if (findable <= 0) return kFALSE;
  if (nCrossedRowsTPC/findable < 0.8) return kFALSE;
  if (TMath::Abs(track->Eta()) > 0.8) return kFALSE;

  return kTRUE;
}

Double_t AliAnalysisTaskpOmegaDibaryon::GetTOFSignal(AliAODTrack& trackHe, Double_t tP) {
  // returns the mass calculated from TOF Signal                                                            
  Double_t mass = 0, time = -1, beta = 0, gamma = 0, length = 0, time0 = 0;
  length = trackHe.GetIntegratedLength();
  time0 = fPIDResponse->GetTOFResponse().GetStartTime(trackHe.P());
  time = trackHe.GetTOFsignal() - time0;
  //cout<<"time = "<<time<<endl;
  if (time > 0 && length > 0) {
    beta = length / (2.99792457999999984e-02 * time);
    if(beta < 1){
      gamma = 1 / TMath::Sqrt(1 - beta*beta);
      mass = (tP / TMath::Sqrt(gamma*gamma - 1));
      return mass*mass;
    }
  }
  return -1;
}
Double_t AliAnalysisTaskpOmegaDibaryon::CalculateBetaTOF(AliAODTrack &track)
{
  Double_t length = track.GetIntegratedLength(); // cm                                                                  
  if(TMath::IsNaN(length)) return -999.0;
  if(length <= 350.0) return -999.0;

  Double_t c = TMath::C(); // m/s                                                                                       
  Double_t end_time = track.GetTOFsignal(); // ps                                                                       
  Double_t start_time = fPIDResponse->GetTOFResponse().GetStartTime(track.GetTPCmomentum()); // ps                      

  if(TMath::IsNaN(end_time)) return -999.0;
  if(TMath::IsNaN(start_time)) return -999.0;

  Double_t time = (end_time - start_time) * 1e-12; // ps -> s                                                           
  Double_t velocity = (length*0.01) / time; // m/s                                                                      
  Double_t beta = velocity / c;

  return beta;
}
Double_t AliAnalysisTaskpOmegaDibaryon::CalculateMassSquareTOF(AliAODTrack &track)
{
  Double_t p = track.P();
  Double_t beta = CalculateBetaTOF(track);

  if(TMath::IsNaN(p)) return -999.0;
  if(TMath::IsNaN(beta)) return -999.0;

  Double_t mass2 = -999.0;

  if(beta > 0.0){
    mass2 = (1/(beta*beta)-1) * (p*p);
  }
  return mass2;
}
Double_t AliAnalysisTaskpOmegaDibaryon::CalculateProtonSigma(Double_t pT, Double_t massSq, Double_t charge)
{
  double SigmaParticle = -999.0;
  if(massSq < -990.0) return SigmaParticle;
  
  TF1 *Mean = new TF1("Mean","[0] + ([1] * (x)) + [2] * pow((1 -([3] / (x))),[4])",0.0,6.0);
  TF1 *Sigma = new TF1("Sigma","[0] + ([1] * (x)) + [2] * pow((1 -([3] / (x))),[4])",0.0,6.0);

  // MetaLHC18q PbPb data 
  if(charge>0){
    Mean->FixParameter(0,0.886972);
    Mean->FixParameter(1,0.00847106);
    Mean->FixParameter(2,-1.94201e-08);
    Mean->FixParameter(3,93.4605);
    Mean->FixParameter(4,3);
    
    Sigma->FixParameter(0,-0.717964);
    Sigma->FixParameter(1,0.0455849);
    Sigma->FixParameter(2,0.654092);
    Sigma->FixParameter(3,-0.213367);
    Sigma->FixParameter(4,0.378298);
  }

  if(charge<0){
    Mean->FixParameter(0,0.883118);
    Mean->FixParameter(1,0.00949616);
    Mean->FixParameter(2,-1.47418e-06);
    Mean->FixParameter(3,22.7463);
    Mean->FixParameter(4,3);

    Sigma->FixParameter(0,0.103316);
    Sigma->FixParameter(1,0.0422522);
    Sigma->FixParameter(2,-0.156961);
    Sigma->FixParameter(3,0.0952412);
    Sigma->FixParameter(4,3.17989);
  }

  double mean  = Mean->Eval(pT);
  double sigma = Sigma->Eval(pT);

  Mean->Delete();
  Sigma->Delete();

  SigmaParticle = (massSq - mean)/(sigma);
  return SigmaParticle;
}

//Proton quality cut                                                                                                     
Bool_t AliAnalysisTaskpOmegaDibaryon::IsQualityProton(AliAODTrack *track)
{
  if (TMath::Abs(track->Eta()) > 0.8) return kFALSE;
  if (track->Pt() > 4.0) return kFALSE;
  if (track->Pt() < 0.5) return kFALSE;
  if (track->GetTPCncls() < 80) return kFALSE;
  if (track->GetTPCnclsS() > 0) return kFALSE;
  if (track->GetTPCCrossedRows() < 70)     return kFALSE;
  if (track->GetTPCchi2perCluster() > 4.0) return kFALSE;
  if (track->GetTPCchi2perNDF() > 4.0)     return kFALSE;
  Float_t nCrossedRowsTPC = track->GetTPCClusterInfo(2,1);
  if (nCrossedRowsTPC < 70) return kFALSE;
  Int_t findable=track->GetTPCNclsF();
  if (findable <= 0) return kFALSE;
  if (nCrossedRowsTPC/findable < 0.83) return kFALSE;
  //if (track->GetITSNcls() < 2) return kFALSE;
  return kTRUE;
}

Double_t AliAnalysisTaskpOmegaDibaryon::CalculateSigmadEdxITS(AliAODTrack &Track)
{
  double SigmaParticle = -999.0;
  double SignalITS = Track.GetITSsignal();
  if(TMath::IsNaN(SignalITS)) return SigmaParticle;

  double p = Track.P();
  double Mass = AliPID::ParticleMass(AliPID::kProton);

  TF1 *Mean = new TF1("Mean","[5]*[5]*AliExternalTrackParam::BetheBlochGeant([5]*x/([6]),[0],[1],[2],[3],[4])",0.0,6.0);
  Mean->FixParameter(0,2.36861e-07);
  Mean->FixParameter(1,-55831.1);
  Mean->FixParameter(2,-238672);
  Mean->FixParameter(3,9.55834);
  Mean->FixParameter(4,17081);
  Mean->FixParameter(5,1);
  Mean->FixParameter(6,Mass);
  
  double mean = Mean->Eval(p);
  Mean->Delete();

  double Resolution = 0.10;//LHC18r,q
  
  double ScaleFactor = 1.0-(Resolution);
  double sigma = (mean*ScaleFactor) - mean;
  if(TMath::Abs(sigma) < 0.0001) return -999.0;
  SigmaParticle = (mean - SignalITS) / (sigma);
  return SigmaParticle;
}
