/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include "TChain.h"
#include "TList.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "THnSparse.h"
#include "TCanvas.h"
#include "TRandom.h"
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TList.h>
#include <TMath.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TProfile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include "AliAnalysisManager.h"
#include "AliAODHandler.h"
#include "AliAODInputHandler.h"
#include "AliInputEventHandler.h"
#include "AliLog.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliMultiplicity.h"
#include "AliCentrality.h"
#include "AliAnalysisTaskPIDBFDptDpt.h"
#include "AliESDVertex.h"
#include "AliESDInputHandler.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAODInputHandler.h"
#include "AliESD.h"
#include "AliStack.h"
#include "AliESDtrackCuts.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliGenHijingEventHeader.h"
#include "AliGenEventHeader.h"
#include "AliPID.h"
#include "AliAODPid.h"
#include "AliPIDResponse.h"
#include "AliAODpidUtil.h"
#include "AliPIDCombined.h"
#include "AliVTrack.h"
#include "AliVParticle.h"
#include "AliMultSelection.h"
#include "AliOADBMultSelection.h"
#include "AliMultEstimator.h"
#include "AliMultVariable.h"
#include "AliMultInput.h"
#include "AliAnalysisUtils.h"


using namespace std;

ClassImp(AliAnalysisTaskPIDBFDptDpt)

AliAnalysisTaskPIDBFDptDpt::AliAnalysisTaskPIDBFDptDpt()
: AliAnalysisTaskSE(),
  fNSigmaPID( 2 ),
  ptUpperLimit( 2.0 ),
  electronNSigmaVeto( 1.0 ),
  fRemoveTracksT0Fill( 0 ),
  fHasTOFPID( 0 ),
  fAODEvent(0),
  fESDEvent(0),
  fInputHandler(0),
  fPIDResponse(0x0),
  _outputHistoList(0),
  _twoPi         ( 2.0 * 3.1415927),
  _eventCount    ( 0),
  _debugLevel    ( 0),
  _singlesOnly   ( 0),
  PIDparticle   ( 0),
  NoContamination   ( 0),
  _useWeights    ( 0),
  _useRapidity   ( 0),
  _sameFilter    ( false),
  _rejectPileup  ( 1),
  _rejectPairConversion ( 0),
  _vertexZMin           ( -6),
  _vertexZMax           (  6),
  _vertexZWidth         (0.5),
  _etaWidth             (0.1),
  _vertexXYMin          ( -10),
  _vertexXYMax          (  10),
  _centralityMethod     (  4),
  _centralityMin        (  0.),
  _centralityMax        (  0.),
  _requestedCharge_1    (   1),
  _requestedCharge_2    (  -1),
  _dcaZMin              ( -3),
  _dcaZMax              (  3.),
  _dcaXYMin             ( -2.4),
  _dcaXYMax             (  2.4),
  _dedxMin              ( 0),
  _dedxMax              ( 100000),
  _nClusterMin          ( 80),
  _trackFilterBit       (0),
  fAnalysisType         ( "RealData" ),
  fSystemType           ( "PbPb" ),
  fExcludeResonancesInMC ( kFALSE ),
  fExcludeElectronsInMC ( kFALSE ),
  particleSpecies       ( 0 ),
  _tpcnclus             ( 50),
  _chi2ndf              (5.),
  _field    ( 1.),
  _nTracks  ( 0 ),
  _mult0    ( 0 ),
  _mult1    ( 0 ),
  _mult2    ( 0 ),
  _mult3    ( 0 ),
  _mult4    ( 0 ),
  _mult5    ( 0 ),
  _mult6    ( 0 ),
  _mult7    ( 0 ),
  _mult8    ( 0 ),
  arraySize ( 5000),
  _id_1(0),
  _charge_1(0),
  _iEtaPhi_1(0),
  _iPt_1(0),
  _pt_1(0),
  _px_1(0),
  _py_1(0),
  _pz_1(0),
  _correction_1(0),
  _dedx_1(0),
  _id_2(0),
  _charge_2(0),
  _iEtaPhi_2(0),
  _iPt_2(0),
  _pt_2(0),
  _px_2(0),
  _py_2(0),
  _pz_2(0),
  _correction_2(0),
  _dedx_2(0),
  _correctionWeight_1(0),
  _correctionWeight_2(0),
  _nBins_M0(500),       _min_M0(0),        _max_M0(10000),          _width_M0(20),
  _nBins_M1(500),       _min_M1(0),        _max_M1(10000),          _width_M1(20),
  _nBins_M2(500),       _min_M2(0),        _max_M2(10000),          _width_M2(20),
  _nBins_M3(500),       _min_M3(0),        _max_M3(10000),          _width_M3(20),
  _nBins_M4(100),       _min_M4(0),        _max_M4(1),              _width_M4(0.01),
  _nBins_M5(100),       _min_M5(0),        _max_M5(1),              _width_M5(0.01),
  _nBins_M6(100),       _min_M6(0),        _max_M6(1),              _width_M6(0.01),
  _nBins_M7(100),       _min_M7(0),        _max_M7(1),              _width_M7(0.01),
  _nBins_M8(100),       _min_M8(0),        _max_M8(1),              _width_M8(0.01),
  _nBins_vertexZ(120),   _min_vertexZ(-6), _max_vertexZ(6),        _width_vertexZ(0.1),

  _nBins_pt_1(18),      _min_pt_1(0.2),    _max_pt_1(2.0),          _width_pt_1(0.1),
  _nBins_phi_1(36),     _min_phi_1(0),     _max_phi_1(2.*3.1415927),_width_phi_1(2.*3.1415927/36.),
  _nBins_eta_1(0),      _min_eta_1(0),  _max_eta_1(0),           _width_eta_1(0.1),

  _nBins_etaPhi_1(0),
  _nBins_etaPhiPt_1(0),
  _nBins_zEtaPhiPt_1(0),

  _nBins_pt_2(18),     _min_pt_2(0.2),     _max_pt_2(2.0),          _width_pt_2(0.1),
  _nBins_phi_2(36),    _min_phi_2(0),      _max_phi_2(2.*3.1415927),_width_phi_2(2.*3.1415927/36.),
  _nBins_eta_2(0),     _min_eta_2(0),     _max_eta_2(0),           _width_eta_2(0.1),

  _nBins_etaPhi_2(0),
  _nBins_etaPhiPt_2(0),
  _nBins_zEtaPhiPt_2(0),
  _nBins_etaPhi_12(0),
  __n1_1(0),
  __n1_2(0),
  __n2_12(0),
  __s1pt_1(0),
  __s1pt_2(0),
  __s2ptpt_12(0),
  __s2NPt_12(0),
  __s2PtN_12(0),
  __n1Nw_1(0),
  __n1Nw_2(0),
  __n2Nw_12(0),
  __s1ptNw_1(0),
  __s1ptNw_2(0),
  __s2ptptNw_12(0),
  __s2NPtNw_12(0),
  __s2PtNNw_12(0),
  __n1_1_vsPt(0),
  __n1_1_vsEtaPhi(0),
  __s1pt_1_vsEtaPhi(0),
  __n1_1_vsZEtaPhiPt(0),
  __n1_2_vsPt(0),
  __n1_2_vsEtaPhi(0),
  __s1pt_2_vsEtaPhi(0),
  __n1_2_vsZEtaPhiPt(0),
  __n2_12_vsPtPt(0),
  __n2_12_vsEtaPhi(0),
  __s2ptpt_12_vsEtaPhi(0),
  __s2PtN_12_vsEtaPhi(0),
  __s2NPt_12_vsEtaPhi(0),
  _weight_1      ( 0    ),
  _weight_2      ( 0    ),
  _eventAccounting ( 0),
  _m0 ( 0),
  _m1 ( 0),
  _m2 ( 0),
  _m3 ( 0),
  _m4 ( 0),
  _m5 ( 0),
  _m6 ( 0),
  _m7 ( 0),
  _m8 ( 0),
  _vertexZ ( 0),

  _Ncluster1  ( 0),
  _Ncluster2  ( 0),

  _t0_1d (0),
  _timeTOF_1d (0),
  _realTOF_1d (0),
  _trackLength (0),
  _trackLength_GetIntegratedLength(0),
  _t0_1d_POI (0),
  _timeTOF_1d_POI (0),
  _realTOF_1d_POI (0),
  _trackLength_POI (0),
  _trackLength_GetIntegratedLength_POI(0),

  _nsigmakaon_1d (0),
  _nsigmaTOFkaon_1d (0), 
  
  _dedx_p (0),
  _beta_p (0), 
  _inverse_beta_p (0),
  _msquare_p (0),

  _dedx_p_POI_AliHelperPID (0),
  _beta_p_POI_AliHelperPID (0),
  _inverse_beta_p_POI_AliHelperPID (0),
  _msquare_p_POI_AliHelperPID (0),
  _nSigmaTOF_p_POI (0),
  _nSigmaTOF_p (0),

  _dedx_p_AliHelperPID_no_Undefined (0),
  _beta_p_AliHelperPID_no_Undefined (0),
  _inverse_beta_p_AliHelperPID_no_Undefined (0),
  _msquare_p_AliHelperPID_no_Undefined (0),
  
  _etadis_POI_AliHelperPID ( 0),
  _ydis_POI_AliHelperPID ( 0),
  _etadis_before_any_cuts ( 0),

//_vZ_y_Pt_POI_AliHelperPID ( 0),
//_vZ_y_eta_POI_AliHelperPID ( 0),

  _y_Pt_AllCh_MCAODTruth ( 0 ),
  _y_Pt_POI_MCAODTruth ( 0 ),
  
  _phidis_POI_AliHelperPID ( 0),
  _phidis_before_any_cuts ( 0),

  _dcaz   ( 0),
  _dcaxy  ( 0),
  _n1_1_vsPt         ( 0),
  _n1_1_vsEtaVsPhi   ( 0),
  _s1pt_1_vsEtaVsPhi ( 0),
  _n1_1_vsZVsEtaVsPhiVsPt ( 0),
  _n1_1_vsM          ( 0),
  _s1pt_1_vsM        ( 0),
  _n1Nw_1_vsM        ( 0),
  _s1ptNw_1_vsM      ( 0),
  _dedxVsP_1         ( 0),
  _corrDedxVsP_1     ( 0),
  _betaVsP_1         ( 0),
  _n1_2_vsPt         ( 0),
  _n1_2_vsEtaVsPhi   ( 0),
  _s1pt_2_vsEtaVsPhi ( 0),
  _n1_2_vsZVsEtaVsPhiVsPt ( 0),
  _n1_2_vsM          ( 0),
  _s1pt_2_vsM        ( 0),
  _n1Nw_2_vsM        ( 0),
  _s1ptNw_2_vsM      ( 0),
  _dedxVsP_2         ( 0),
  _corrDedxVsP_2     ( 0),
  _betaVsP_2         ( 0),
  _n2_12_vsEtaPhi    ( 0),
  _n2_12_vsPtVsPt    ( 0),
  _s2PtPt_12_vsEtaPhi( 0),
  _s2PtN_12_vsEtaPhi ( 0),
  _s2NPt_12_vsEtaPhi ( 0),
  _n2_12_vsM         ( 0),
  _s2PtPt_12_vsM     ( 0),
  _s2PtN_12_vsM      ( 0),
  _s2NPt_12_vsM      ( 0),
  _n2Nw_12_vsM       ( 0),
  _s2PtPtNw_12_vsM   ( 0),
  _s2PtNNw_12_vsM    ( 0),
  _s2NPtNw_12_vsM    ( 0),
  _invMass           ( 0),
  _invMassElec       ( 0),
  n1Name("NA"),
  n1NwName("NA"),
  n2Name("NA"),
  n2NwName("NA"),
  n3Name("NA"),
  n1n1Name("NA"),
  n1n1n1Name("NA"),
  n2n1Name("NA"),
  r1Name("NA"),
  r2Name("NA"),
  r3Name("NA"),
  r2r1Name("NA"),
  c2Name("NA"),
  c3Name("NA"),
  d3Name("NA"),
  p3Name("NA"),
  cName("NA"),

  intR2Name("NA"),
  binCorrName("NA"),
  intBinCorrName("NA"),

  countsName("NA"),
  part_1_Name("NA"),
  part_2_Name("NA"),
  part_3_Name("NA"),
  pair_12_Name("NA"),
  pair_13_Name("NA"),
  pair_23_Name("NA"),
  tripletName("NA"),

  avg("NA"),
  avgName("NA"),
  sumName("NA"),
  s1ptName("NA"),
  s1ptNwName("NA"),
  s1DptName("NA"),

  s2PtPtName("NA"),
  s2NPtName("NA"),
  s2PtNName("NA"),
  s2DptDptName("NA"),

  s2PtPtNwName("NA"),
  s2NPtNwName("NA"),
  s2PtNNwName("NA"),

  ptName("NA"),
  ptptName("NA"),
  pt1pt1Name("NA"),
  DptName("NA"),
  DptDptName("NA"),
  RDptDptName("NA"),
  nPtName("NA"),
  ptNName("NA"),
  seanName("NA"),

  _title_counts("NA"),

  _title_m0("NA"),
  _title_m1("NA"),
  _title_m2("NA"),
  _title_m3("NA"),
  _title_m4("NA"),
  _title_m5("NA"),
  _title_m6("NA"),
  _title_m7("NA"),
  _title_m8("NA"),

  _title_eta_1("NA"),
  _title_phi_1("NA"),
  _title_pt_1("NA"),
  _title_etaPhi_1("NA"),
  _title_n_1("NA"),
  _title_SumPt_1("NA"),
  _title_AvgPt_1("NA"),
  _title_AvgN_1("NA"),
  _title_AvgSumPt_1("NA"),

  _title_eta_2("NA"),
  _title_phi_2("NA"),
  _title_pt_2("NA"),
  _title_etaPhi_2("NA"),
  _title_n_2("NA"),
  _title_SumPt_2("NA"),
  _title_AvgPt_2("NA"),
  _title_AvgN_2("NA"),
  _title_AvgSumPt_2("NA"),

  _title_etaPhi_12("NA"),

  _title_AvgN2_12("NA"),
  _title_AvgSumPtPt_12("NA"),
  _title_AvgSumPtN_12("NA"),
  _title_AvgNSumPt_12("NA"),

  vsZ("NA"),
  vsM("NA"),
  vsPt("NA"),
  vsPhi("NA"),
  vsEta("NA"),
  vsEtaPhi("NA"),
  vsPtVsPt("NA"),
  fUtils(0)
{
  printf("Default constructor called \n");  
  printf("passed \n ");
}

AliAnalysisTaskPIDBFDptDpt::AliAnalysisTaskPIDBFDptDpt(const TString & name)
: AliAnalysisTaskSE(name),
  fNSigmaPID( 2 ),
  ptUpperLimit( 2.0 ),
  electronNSigmaVeto( 1.0 ),
  fRemoveTracksT0Fill( 0 ),
  fHasTOFPID( 0 ),
  fAODEvent(0),
  fESDEvent(0),
  fInputHandler(0),
  fPIDResponse(0),
  _outputHistoList(0),
  _twoPi         ( 2.0 * 3.1415927),
  _eventCount    ( 0),
  _debugLevel    ( 0),
  _singlesOnly   ( 0),
  PIDparticle    ( 0),
  NoContamination   ( 0),
  _useWeights    ( 0),
  _useRapidity   ( 0),
  _sameFilter    ( false),
  _rejectPileup  ( 1),
  _rejectPairConversion ( 0),
  _vertexZMin           ( -6.),
  _vertexZMax           (  6.),
  _vertexZWidth         (0.5 ),
  _etaWidth             (0.1),
  _vertexXYMin          ( -10.),
  _vertexXYMax          (  10.),
  _centralityMethod     (  4),
  _centralityMin        (  0.),
  _centralityMax        (  1.),
  _requestedCharge_1    (   1),
  _requestedCharge_2    (  -1),
  _dcaZMin              ( -3),
  _dcaZMax              (  3.),
  _dcaXYMin             ( -2.4),
  _dcaXYMax             (  2.4),
  _dedxMin              ( 0),
  _dedxMax              ( 100000),
  _nClusterMin          ( 80),
  _trackFilterBit       ( 0),
  fAnalysisType         ( "RealData" ),
  fSystemType           ( "PbPb" ),
  fExcludeResonancesInMC ( kFALSE ),
  fExcludeElectronsInMC ( kFALSE ),
  particleSpecies       ( 0 ),
  _tpcnclus             ( 50),
  _chi2ndf              (5.),
  _field    ( 1.),
  _nTracks  ( 0 ),
  _mult0    ( 0 ),
  _mult1    ( 0 ),
  _mult2    ( 0 ),
  _mult3    ( 0 ),
  _mult4    ( 0 ),
  _mult5    ( 0 ),
  _mult6    ( 0 ),
  _mult7    ( 0 ),
  _mult8    ( 0 ),
  arraySize ( 5000),
  _id_1(0),
  _charge_1(0),
  _iEtaPhi_1(0),
  _iPt_1(0),
  _pt_1(0),
  _px_1(0),
  _py_1(0),
  _pz_1(0),
  _correction_1(0),
  _dedx_1(0),
  _id_2(0),
  _charge_2(0),
  _iEtaPhi_2(0),
  _iPt_2(0),
  _pt_2(0),
  _px_2(0),
  _py_2(0),
  _pz_2(0),
  _correction_2(0),
  _dedx_2(0),
  _correctionWeight_1(0),
  _correctionWeight_2(0),
  _nBins_M0(500),       _min_M0(0),        _max_M0(10000),          _width_M0(20),
  _nBins_M1(500),       _min_M1(0),        _max_M1(10000),          _width_M1(20),
  _nBins_M2(500),       _min_M2(0),        _max_M2(10000),          _width_M2(20),
  _nBins_M3(500),       _min_M3(0),        _max_M3(10000),          _width_M3(20),
  _nBins_M4(100),       _min_M4(0),        _max_M4(1),              _width_M4(0.01),
  _nBins_M5(100),       _min_M5(0),        _max_M5(1),              _width_M5(0.01),
  _nBins_M6(100),       _min_M6(0),        _max_M6(1),              _width_M6(0.01),
  _nBins_M7(100),       _min_M7(0),        _max_M7(1),              _width_M7(0.01),
  _nBins_M8(100),       _min_M8(0),        _max_M8(1),              _width_M8(0.01),
  _nBins_vertexZ(120),   _min_vertexZ(-6), _max_vertexZ(6),        _width_vertexZ(0.1),

  _nBins_pt_1(18),      _min_pt_1(0.2),    _max_pt_1(2.0),          _width_pt_1(0.1),
  _nBins_phi_1(36),     _min_phi_1(0),     _max_phi_1(2.*3.1415927),_width_phi_1(2.*3.1415927/36.),
  _nBins_eta_1(0),      _min_eta_1(0),    _max_eta_1(0),           _width_eta_1(0.1),

  _nBins_etaPhi_1(0),
  _nBins_etaPhiPt_1(0),
  _nBins_zEtaPhiPt_1(0),

  _nBins_pt_2(18),     _min_pt_2(0.2),     _max_pt_2(2.0),          _width_pt_2(0.1),
  _nBins_phi_2(36),    _min_phi_2(0),      _max_phi_2(2.*3.1415927),_width_phi_2(2.*3.1415927/36.),
  _nBins_eta_2(0),    _min_eta_2(0),     _max_eta_2(0),           _width_eta_2(0.1),

  _nBins_etaPhi_2(0),
  _nBins_etaPhiPt_2(0),
  _nBins_zEtaPhiPt_2(0),
  _nBins_etaPhi_12(0),
  __n1_1(0),
  __n1_2(0),
  __n2_12(0),
  __s1pt_1(0),
  __s1pt_2(0),
  __s2ptpt_12(0),
  __s2NPt_12(0),
  __s2PtN_12(0),
  __n1Nw_1(0),
  __n1Nw_2(0),
  __n2Nw_12(0),
  __s1ptNw_1(0),
  __s1ptNw_2(0),
  __s2ptptNw_12(0),
  __s2NPtNw_12(0),
  __s2PtNNw_12(0),
  __n1_1_vsPt(0),
  __n1_1_vsEtaPhi(0),
  __s1pt_1_vsEtaPhi(0),
  __n1_1_vsZEtaPhiPt(0),
  __n1_2_vsPt(0),
  __n1_2_vsEtaPhi(0),
  __s1pt_2_vsEtaPhi(0),
  __n1_2_vsZEtaPhiPt(0),
  __n2_12_vsPtPt(0),
  __n2_12_vsEtaPhi(0),
  __s2ptpt_12_vsEtaPhi(0),
  __s2PtN_12_vsEtaPhi(0),
  __s2NPt_12_vsEtaPhi(0),
  _weight_1        ( 0    ),
  _weight_2        ( 0    ),
  _eventAccounting ( 0),
  _m0 ( 0),
  _m1 ( 0),
  _m2 ( 0),
  _m3 ( 0),
  _m4 ( 0),
  _m5 ( 0),
  _m6 ( 0),
  _m7 ( 0),
  _m8 ( 0),
  _vertexZ ( 0),
  _Ncluster1  ( 0),
  _Ncluster2  ( 0),

  _t0_1d (0),
  _timeTOF_1d (0),
  _realTOF_1d (0),
  _trackLength (0),
  _trackLength_GetIntegratedLength(0),
  _t0_1d_POI (0),
  _timeTOF_1d_POI (0),
  _realTOF_1d_POI (0),
  _trackLength_POI (0),
  _trackLength_GetIntegratedLength_POI(0),
  
  _nsigmakaon_1d (0),
  _nsigmaTOFkaon_1d (0), 
  
  _dedx_p (0),
  _beta_p (0), 
  _inverse_beta_p (0),
  _msquare_p (0),

  _dedx_p_POI_AliHelperPID (0),
  _beta_p_POI_AliHelperPID (0),
  _inverse_beta_p_POI_AliHelperPID (0),
  _msquare_p_POI_AliHelperPID (0),
  _nSigmaTOF_p_POI (0),
  _nSigmaTOF_p (0),

  _dedx_p_AliHelperPID_no_Undefined (0),
  _beta_p_AliHelperPID_no_Undefined (0),
  _inverse_beta_p_AliHelperPID_no_Undefined (0),
  _msquare_p_AliHelperPID_no_Undefined (0),
  
  _etadis_POI_AliHelperPID ( 0),
  _ydis_POI_AliHelperPID ( 0),
  _etadis_before_any_cuts ( 0),

  //_vZ_y_Pt_POI_AliHelperPID ( 0),
  //_vZ_y_eta_POI_AliHelperPID ( 0),

  _y_Pt_AllCh_MCAODTruth ( 0 ),
  _y_Pt_POI_MCAODTruth ( 0 ),
  
  _phidis_POI_AliHelperPID ( 0),
  _phidis_before_any_cuts ( 0),

  _dcaz ( 0),
  _dcaxy ( 0),
  _n1_1_vsPt         ( 0),
  _n1_1_vsEtaVsPhi   ( 0),
  _s1pt_1_vsEtaVsPhi ( 0),
  _n1_1_vsZVsEtaVsPhiVsPt ( 0),
  _n1_1_vsM          ( 0),
  _s1pt_1_vsM        ( 0),
  _n1Nw_1_vsM        ( 0),
  _s1ptNw_1_vsM      ( 0),
  _dedxVsP_1         ( 0),
  _corrDedxVsP_1     ( 0),
  _betaVsP_1         ( 0),
  _n1_2_vsPt         ( 0),
  _n1_2_vsEtaVsPhi   ( 0),
  _s1pt_2_vsEtaVsPhi ( 0),
  _n1_2_vsZVsEtaVsPhiVsPt ( 0),
  _n1_2_vsM          ( 0),
  _s1pt_2_vsM        ( 0),
  _n1Nw_2_vsM        ( 0),
  _s1ptNw_2_vsM      ( 0),
  _dedxVsP_2         ( 0),
  _corrDedxVsP_2     ( 0),
  _betaVsP_2         ( 0),
  _n2_12_vsEtaPhi    ( 0),
  _n2_12_vsPtVsPt    ( 0),
  _s2PtPt_12_vsEtaPhi( 0),
  _s2PtN_12_vsEtaPhi ( 0),
  _s2NPt_12_vsEtaPhi ( 0),
  _n2_12_vsM         ( 0),
  _s2PtPt_12_vsM     ( 0),
  _s2PtN_12_vsM      ( 0),
  _s2NPt_12_vsM      ( 0),
  _n2Nw_12_vsM       ( 0),
  _s2PtPtNw_12_vsM   ( 0),
  _s2PtNNw_12_vsM    ( 0),
  _s2NPtNw_12_vsM    ( 0),
  _invMass           ( 0),
  _invMassElec       ( 0),
  n1Name("NA"),
  n1NwName("NA"),
  n2Name("NA"),
  n2NwName("NA"),
  n3Name("NA"),
  n1n1Name("NA"),
  n1n1n1Name("NA"),
  n2n1Name("NA"),
  r1Name("NA"),
  r2Name("NA"),
  r3Name("NA"),
  r2r1Name("NA"),
  c2Name("NA"),
  c3Name("NA"),
  d3Name("NA"),
  p3Name("NA"),
  cName("NA"),

  intR2Name("NA"),
  binCorrName("NA"),
  intBinCorrName("NA"),

  countsName("NA"),
  part_1_Name("NA"),
  part_2_Name("NA"),
  part_3_Name("NA"),
  pair_12_Name("NA"),
  pair_13_Name("NA"),
  pair_23_Name("NA"),
  tripletName("NA"),

  avg("NA"),
  avgName("NA"),
  sumName("NA"),
  s1ptName("NA"),
  s1ptNwName("NA"),
  s1DptName("NA"),

  s2PtPtName("NA"),
  s2NPtName("NA"),
  s2PtNName("NA"),
  s2DptDptName("NA"),

  s2PtPtNwName("NA"),
  s2NPtNwName("NA"),
  s2PtNNwName("NA"),

  ptName("NA"),
  ptptName("NA"),
  pt1pt1Name("NA"),
  DptName("NA"),
  DptDptName("NA"),
  RDptDptName("NA"),
  nPtName("NA"),
  ptNName("NA"),
  seanName("NA"),

  _title_counts("NA"),

  _title_m0("NA"),
  _title_m1("NA"),
  _title_m2("NA"),
  _title_m3("NA"),
  _title_m4("NA"),
  _title_m5("NA"),
  _title_m6("NA"),
  _title_m7("NA"),
  _title_m8("NA"),

  _title_eta_1("NA"),
  _title_phi_1("NA"),
  _title_pt_1("NA"),
  _title_etaPhi_1("NA"),
  _title_n_1("NA"),
  _title_SumPt_1("NA"),
  _title_AvgPt_1("NA"),
  _title_AvgN_1("NA"),
  _title_AvgSumPt_1("NA"),

  _title_eta_2("NA"),
  _title_phi_2("NA"),
  _title_pt_2("NA"),
  _title_etaPhi_2("NA"),
  _title_n_2("NA"),
  _title_SumPt_2("NA"),
  _title_AvgPt_2("NA"),
  _title_AvgN_2("NA"),
  _title_AvgSumPt_2("NA"),

  _title_etaPhi_12("NA"),

  _title_AvgN2_12("NA"),
  _title_AvgSumPtPt_12("NA"),
  _title_AvgSumPtN_12("NA"),
  _title_AvgNSumPt_12("NA"),

  vsZ("NA"),
  vsM("NA"),
  vsPt("NA"),
  vsPhi("NA"),
  vsEta("NA"),
  vsEtaPhi("NA"),
  vsPtVsPt("NA"),
  fUtils(0)
{
  // Au-Au added this block of code to use his own PID functions
  for( Int_t ipart = 0; ipart < 4; ipart++ )
    for( Int_t ipid = 0; ipid < 2; ipid++ )
      fnsigmas[ipart][ipid] = 999.;
  
  printf("2nd constructor called ");    
  DefineOutput(0, TList::Class());   
  printf("passed  ");   
}

AliAnalysisTaskPIDBFDptDpt::~AliAnalysisTaskPIDBFDptDpt()
{
    
}

void AliAnalysisTaskPIDBFDptDpt::UserCreateOutputObjects()
{
  OpenFile(0);
  _outputHistoList = new TList();
  _outputHistoList->SetOwner();
    
  //if ( _singlesOnly )   _outputHistoList -> Add( fHelperPID -> GetOutputList() ); // add AliHelperPIDBFDptDpt object output list to task output list only for singles
    
  _nBins_M0 = 500; _min_M0   = 0.;    _max_M0    = 5000.;  _width_M0 = (_max_M0-_min_M0)/_nBins_M0;
  _nBins_M1 = 500; _min_M1   = 0.;    _max_M1    = 5000.;  _width_M1 = (_max_M1-_min_M1)/_nBins_M1;
  _nBins_M2 = 500; _min_M2   = 0.;    _max_M2    = 5000.;  _width_M2 = (_max_M2-_min_M2)/_nBins_M2;
  _nBins_M3 = 500; _min_M3   = 0.;    _max_M3    = 5000.;  _width_M3 = (_max_M3-_min_M3)/_nBins_M3;
  _nBins_M4 = 100; _min_M4   = 0.;    _max_M4    = 100.;   _width_M4 = (_max_M4-_min_M4)/_nBins_M4;
  _nBins_M5 = 100; _min_M5   = 0.;    _max_M5    = 100.;   _width_M5 = (_max_M5-_min_M5)/_nBins_M5;
  _nBins_M6 = 100; _min_M6   = 0.;    _max_M6    = 100.;   _width_M6 = (_max_M6-_min_M6)/_nBins_M6;
  _nBins_M7 = 100; _min_M7   = 0.;    _max_M7    = 100.;   _width_M7 = (_max_M7-_min_M7)/_nBins_M7;
  _nBins_M8 = 100; _min_M8   = 0.;    _max_M8    = 100.;   _width_M8 = (_max_M8-_min_M8)/_nBins_M8;
    
  _min_vertexZ       = _vertexZMin;
  _max_vertexZ       = _vertexZMax;
  _width_vertexZ     = _vertexZWidth;
  _width_eta_1       = _etaWidth;
  _width_eta_2       = _etaWidth;
    
  _nBins_vertexZ     = int( 0.5 + ( _max_vertexZ - _min_vertexZ) / _width_vertexZ );
  _nBins_pt_1        = int( 0.5 + ( _max_pt_1 -_min_pt_1 ) / _width_pt_1 );
  _nBins_eta_1       = int( 0.5 + ( _max_eta_1-_min_eta_1 ) / _width_eta_1 );
  _width_phi_1       = ( _max_phi_1  - _min_phi_1 )  / _nBins_phi_1;
  _nBins_etaPhi_1    = _nBins_phi_1    * _nBins_eta_1;
  _nBins_etaPhiPt_1  = _nBins_etaPhi_1 * _nBins_pt_1;
  _nBins_zEtaPhiPt_1 = _nBins_vertexZ  * _nBins_etaPhiPt_1;
    
  _nBins_pt_2        =  int(0.5+ (_max_pt_2 -_min_pt_2 )/_width_pt_2);
  _nBins_eta_2       = int(0.5+ (_max_eta_2-_min_eta_2)/_width_eta_2);
  _width_phi_2       = (_max_phi_2  - _min_phi_2)  /_nBins_phi_2;
  _nBins_etaPhi_2    = _nBins_phi_2    * _nBins_eta_2;
  _nBins_etaPhiPt_2  = _nBins_etaPhi_2 * _nBins_pt_2;
  _nBins_zEtaPhiPt_2 = _nBins_vertexZ  * _nBins_etaPhiPt_2;
  _nBins_etaPhi_12   = _nBins_etaPhi_1 * _nBins_etaPhi_2;
    
  _id_1       = new int[arraySize];
  _charge_1   = new int[arraySize];
  _iEtaPhi_1  = new int[arraySize];
  _iPt_1      = new int[arraySize];
  _pt_1       = new float[arraySize];
  _px_1       = new float[arraySize];
  _py_1       = new float[arraySize];
  _pz_1       = new float[arraySize];
  _correction_1 = new float[arraySize];
  _dedx_1     = new float[arraySize];
  __n1_1_vsPt              = getDoubleArray(_nBins_pt_1,        0.);
  __n1_1_vsEtaPhi          = getDoubleArray(_nBins_etaPhi_1,    0.);
  __s1pt_1_vsEtaPhi        = getDoubleArray(_nBins_etaPhi_1,    0.);
  __n1_1_vsZEtaPhiPt       = getFloatArray(_nBins_zEtaPhiPt_1,  0.);
    
    
  if (_requestedCharge_2!=_requestedCharge_1)
    {
      _sameFilter = 0;
      //particle 2
      _id_2       = new int[arraySize];
      _charge_2   = new int[arraySize];
      _iEtaPhi_2  = new int[arraySize];
      _iPt_2      = new int[arraySize];
      _pt_2       = new float[arraySize];
      _px_2       = new float[arraySize];
      _py_2       = new float[arraySize];
      _pz_2       = new float[arraySize];
      _correction_2 = new float[arraySize];
      _dedx_2       = new float[arraySize];
      __n1_2_vsPt              = getDoubleArray(_nBins_pt_2,        0.);
      __n1_2_vsEtaPhi          = getDoubleArray(_nBins_etaPhi_2,    0.);
      __s1pt_2_vsEtaPhi        = getDoubleArray(_nBins_etaPhi_2,    0.);
      __n1_2_vsZEtaPhiPt       = getFloatArray(_nBins_zEtaPhiPt_2, 0.);
        
    }
    
  __n2_12_vsPtPt           = getDoubleArray(_nBins_pt_1*_nBins_pt_2,0.);
  __n2_12_vsEtaPhi         = getFloatArray(_nBins_etaPhi_12,       0.);
  __s2ptpt_12_vsEtaPhi     = getFloatArray(_nBins_etaPhi_12,       0.);
  __s2PtN_12_vsEtaPhi      = getFloatArray(_nBins_etaPhi_12,       0.);
  __s2NPt_12_vsEtaPhi      = getFloatArray(_nBins_etaPhi_12,       0.);
    
  // Setup all the labels needed.
    
  part_1_Name   = "_1";
  part_2_Name   = "_2";
  pair_12_Name  = "_12";
    
  n1Name     = "n1";
  n2Name     = "n2";
  n1NwName   = "n1Nw";
  n2NwName   = "n2Nw";
  r1Name     = "r1";
  r2Name     = "r2";
  r3Name     = "r3";
  r2r1Name   = "r2r1";
  c2Name     = "c2";
  c3Name     = "c3";
  d3Name     = "d3";
  p3Name     = "p3";
  cName      = "sean";
    
  intR2Name       = "intR2";
  binCorrName     = "binCorr";
  intBinCorrName  = "intBinCorr";
    
  avgName      = "avg";
  sumName      = "sum";
  s1ptName     = "sumPt";
  s1ptNwName   = "sumPtNw";
  s1DptName    = "sumDpt";
  s2PtPtName   = "sumPtPt";
  s2PtPtNwName = "sumPtPtNw";
  s2DptDptName = "sumDptDpt";
  s2NPtName    = "sumNPt";
  s2NPtNwName  = "sumNPtNw";
  s2PtNName    = "sumPtN";
  s2NPtNwName  = "sumNPtNw";
  s2PtNNwName  = "sumPtNNw";
  ptName       = "avgPt";
  ptptName     = "avgPtPt";
  pt1pt1Name   = "avgPtavgPt";
  DptName      = "avgDpt";
  DptDptName   = "avgDptDpt";
  RDptDptName  = "relDptDpt"; // ratio of avgDptDpt by avgPt*avgPt
  nPtName      = "avgNpt";
  ptNName      = "avgPtN";
  seanName     = "seanC";
    
  _title_counts = "yield";
    
  _title_m0     = "M_{0}";
  _title_m1     = "M_{1}";
  _title_m2     = "M_{2}";
  _title_m3     = "M_{3}";
  _title_m4     = "V0Centrality";
  _title_m5     = "TrkCentrality";
  _title_m6     = "SpdCentrality";
  _title_m7     = "V0ACentrality";
  _title_m8     = "V0CCentrality";
    
  _title_eta_1       = "#eta_{1}";
  _title_phi_1       = "#varphi_{1} (radian)";
  _title_etaPhi_1    = "#eta_{1}#times#varphi_{1}";
  _title_pt_1        = "p_{t,1} (GeV/c)";
  _title_n_1         = "n_{1}";
  _title_SumPt_1     = "#Sigma p_{t,1} (GeV/c)";
  _title_AvgPt_1     = "#LT p_{t,1} #GT (GeV/c)";
  _title_AvgN_1      = "#LT n_{1} #GT";
  _title_AvgSumPt_1  = "#LT #Sigma p_{t,1} #GT (GeV/c)";
    
  _title_eta_2       = "#eta_{2}";
  _title_phi_2       = "#varphi_{2} (radian)";
  _title_etaPhi_2    = "#eta_{2}#times#varphi_{2}";
  _title_pt_2        = "p_{t,2} (GeV/c)";
  _title_n_2         = "n_{2}";
  _title_SumPt_2     = "#Sigma p_{t,1} (GeV/c)";
  _title_AvgPt_2     = "#LT p_{t,2} #GT (GeV/c)";
  _title_AvgN_2      = "#LT n_{2} #GT";
  _title_AvgSumPt_2  = "#LT #Sigma p_{t,2} #GT (GeV/c)";
    
  _title_etaPhi_12   = "#eta_{1}#times#varphi_{1}#times#eta_{2}#times#varphi_{2}";
    
  _title_AvgN2_12       = "#LT n_{2} #GT";;
  _title_AvgSumPtPt_12  = "#LT #Sigma p_{t,1}p_{t,2} #GT";;
  _title_AvgSumPtN_12   = "#LT #Sigma p_{t,1}N #GT";;
  _title_AvgNSumPt_12   = "#LT N#Sigma p_{t,2} #GT";;
    
    
  vsZ         = "_vsZ";
  vsM         = "_vsM";
  vsPt         = "_vsPt";
  vsPhi        = "_vsPhi";
  vsEta        = "_vsEta";
  vsEtaPhi     = "_vsEtaPhi";
  vsPtVsPt     = "_vsPtVsPt";
    
    
  if (_useWeights)
    {
      int iZ, iEtaPhi, iPt;
      int iZ1,iEtaPhi1,iPt1;
      int a, b;
      if (_weight_1)
        {
	  _correctionWeight_1 = new float[_nBins_vertexZ*_nBins_etaPhi_1*_nBins_pt_1];
	  a = _nBins_pt_1;
	  b = _nBins_etaPhi_1*_nBins_pt_1;
	  for (iZ=0,iZ1=1; iZ<_nBins_vertexZ; iZ++, iZ1++)
            {
	      for (iEtaPhi=0,iEtaPhi1=1; iEtaPhi<_nBins_etaPhi_1; iEtaPhi++, iEtaPhi1++)
                {
		  for (iPt=0,iPt1=1; iPt<_nBins_pt_1; iPt++, iPt1++)
                    {
		      _correctionWeight_1[iZ*b+iEtaPhi*a+iPt] = _weight_1->GetBinContent(iZ1,iEtaPhi1,iPt1);
                    }
                }
            }
        } // _weight_1
      else
        {
	  AliError("AliAnalysisTaskPIDBFDptDpt:: _weight_1 is a null pointer.");
	  return;
        }
      if (!_sameFilter)
        {
	  if (_weight_2)
            {
	      _correctionWeight_2 = new float[_nBins_vertexZ*_nBins_etaPhi_2*_nBins_pt_2];
	      a = _nBins_pt_2;
	      b = _nBins_etaPhi_2*_nBins_pt_2;
	      for (iZ=0,iZ1=1; iZ<_nBins_vertexZ; iZ++, iZ1++)
                {
		  for (iEtaPhi=0,iEtaPhi1=1; iEtaPhi<_nBins_etaPhi_2; iEtaPhi++, iEtaPhi1++)
                    {
		      for (iPt=0,iPt1=1; iPt<_nBins_pt_2; iPt++, iPt1++)
                        {
			  _correctionWeight_2[iZ*b+iEtaPhi*a+iPt] = _weight_2->GetBinContent(iZ1,iEtaPhi1,iPt1);
                        }
                    }
                }
            } // _weight_2
	  else
            {
	      AliError("AliAnalysisTaskPIDBFDptDpt:: _weight_1 is a null pointer.");
	      return;
            }
        }
    }
    
  fUtils = new AliAnalysisUtils();
	
  createHistograms();

  PostData(0,_outputHistoList);
  
  //cout<< "AliAnalysisTaskPIDBFDptDpt::CreateOutputObjects() DONE " << endl;
    
}

void  AliAnalysisTaskPIDBFDptDpt::createHistograms()
{
  AliInfo(" AliAnalysisTaskPIDBFDptDpt::createHistoHistograms() Creating Event Histos");
  TString name;

  // histos for events:
  name = "eventAccounting";  _eventAccounting      = createHisto1D(name,name,10, -0.5, 9.5, "event Code", _title_counts);    
  name = "m0"; _m0      = createHisto1D(name,name,_nBins_M1, _min_M1, _max_M1, _title_m0, _title_counts);
  name = "m1"; _m1      = createHisto1D(name,name,_nBins_M1, _min_M1, _max_M1, _title_m1, _title_counts);
  name = "m2"; _m2      = createHisto1D(name,name,_nBins_M2, _min_M2, _max_M2, _title_m2, _title_counts);
  name = "m3"; _m3      = createHisto1D(name,name,_nBins_M3, _min_M3, _max_M3, _title_m3, _title_counts);
  name = "m4"; _m4      = createHisto1D(name,name,_nBins_M4, _min_M4, _max_M4, _title_m4, _title_counts);
  name = "m5"; _m5      = createHisto1D(name,name,_nBins_M5, _min_M5, _max_M5, _title_m5, _title_counts);
  name = "m6"; _m6      = createHisto1D(name,name,_nBins_M6, _min_M6, _max_M6, _title_m6, _title_counts);
  name = "m7"; _m7      = createHisto1D(name,name,_nBins_M7, _min_M7, _max_M7, _title_m7, _title_counts);
  name = "m8"; _m8      = createHisto1D(name,name,_nBins_M8, _min_M8, _max_M8, _title_m8, _title_counts);
  name = "zV"; _vertexZ = createHisto1D(name,name,_nBins_vertexZ, _min_vertexZ, _max_vertexZ, "z-Vertex (cm)", _title_counts);

  // histos for tracks:
  if ( _singlesOnly )
    {
      name = "etadis_POI_AliHelperPID";          _etadis_POI_AliHelperPID   = createHisto1F(name,name, 200, -1.0, 1.0, "#eta","counts");
      name = "ydis_POI_AliHelperPID";            _ydis_POI_AliHelperPID   = createHisto1F(name,name, 200, -1.0, 1.0, "y","counts");
      name = "etadis_before_any_cuts";            _etadis_before_any_cuts   = createHisto1F(name,name, 200, -1.0, 1.0, "#eta","counts");       
      name = "phidis_POI_AliHelperPID";          _phidis_POI_AliHelperPID   = createHisto1F(name,name, 360, 0.0, 6.4, "#phi","counts");
      name = "phidis_before_any_cuts";            _phidis_before_any_cuts   = createHisto1F(name,name, 360, 0.0, 6.4, "#phi","counts");   
      name = "DCAz";    _dcaz     = createHisto1F(name,name, 500, -5.0, 5.0, "dcaZ","counts");
      name = "DCAxy";   _dcaxy    = createHisto1F(name,name, 500, -5.0, 5.0, "dcaXY","counts");    
      name = "Nclus1";   _Ncluster1    = createHisto1F(name,name, 200, 0, 200, "Ncluster1","counts");
      name = "Nclus2";   _Ncluster2    = createHisto1F(name,name, 200, 0, 200, "Ncluster2","counts");
      name = "T0";       _t0_1d    = createHisto1F(name,name, 20000, -10000, 10000, "T0","counts");
      name = "timeTOF";  _timeTOF_1d    = createHisto1F(name,name, 32000, -2000, 30000, "timeTOF","counts");
      name = "realTOF";   _realTOF_1d    = createHisto1F(name,name, 32000, -2000, 30000, "realTOF","counts");
      name = "trackLength";   _trackLength    = createHisto1F(name,name, 20000, -10, 10, "track Length","counts");
      name = "trackLength_GetIntegratedLength";   _trackLength_GetIntegratedLength = createHisto1F(name,name, 20000, -1000, 1000, "track Length by GetIntegratedLength()","counts");
      name = "T0_POI";       _t0_1d_POI    = createHisto1F(name,name, 20000, -10000, 10000, "T0","counts");
      name = "timeTOF_POI";  _timeTOF_1d_POI    = createHisto1F(name,name, 32000, -2000, 30000, "timeTOF","counts");
      name = "realTOF_POI";   _realTOF_1d_POI    = createHisto1F(name,name, 32000, -2000, 30000, "realTOF","counts");
      name = "trackLength_POI";   _trackLength_POI    = createHisto1F(name,name, 20000, -10, 10, "track Length","counts");
      name = "trackLength_GetIntegratedLength_POI";   _trackLength_GetIntegratedLength_POI = createHisto1F(name,name, 20000, -1000, 1000, "track Length by GetIntegratedLength()","counts");
      name = "dedx_p";   _dedx_p = createHisto2F(name,name, 1980, 0.2, 20,  200, 0, 200,  "p", "dedx","counts");
      name = "beta_p";   _beta_p = createHisto2F(name,name, 500, 0, 5,  100, 0.1, 1.1,  "p", "beta","counts");
      name = "inverse_beta_p";   _inverse_beta_p = createHisto2F(name,name, 500, 0, 5, 200, 0.6, 2.6,  "p", "1/#beta","counts");
      name = "msquare_p";   _msquare_p = createHisto2F(name,name, 500, 0, 5,  200, -0.5, 1.5,  "p", "mass square",  "counts");   
      name = "nsigmakaon_1d";   _nsigmakaon_1d    = createHisto1F(name,name, 2000, 0, 200, "nsigmakaon","counts");
      name = "nsigmaTOFkaon_1d";   _nsigmaTOFkaon_1d    = createHisto1F(name,name, 2000, 0, 200, "nsigmaTOFkaon","counts");  
      name = "msquare_p_POI_AliHelperPID";   _msquare_p_POI_AliHelperPID = createHisto2F(name,name, 500, 0, 5,  200, -0.5, 1.5,  "p", "mass square","counts");
      name = "msquare_p_AliHelperPID_no_Undefined";   _msquare_p_AliHelperPID_no_Undefined = createHisto2F(name,name, 500, 0, 5,  200, -0.5, 1.5,  "p", "mass square","counts");
      name = "nSigmaTOF_p";       _nSigmaTOF_p     = createHisto2F(name,name, 500, 0, 5, 120, -20, 100.0, "p", "nSigmaTOF","counts");
      name = "nSigmaTOF_p_POI";   _nSigmaTOF_p_POI = createHisto2F(name,name, 500, 0, 5, 80, -4.0, 4.0,  "p", "nSigmaTOF_POI","counts");
      name = "dedx_p_POI_AliHelperPID";   _dedx_p_POI_AliHelperPID = createHisto2F(name,name, 1980, 0.2, 20,  200, 0, 200,  "p", "dedx","counts");
      name = "dedx_p_AliHelperPID_no_Undefined";   _dedx_p_AliHelperPID_no_Undefined = createHisto2F(name,name, 1980, 0.2, 20,  200, 0, 200,  "p", "dedx","counts");    
      name = "beta_p_POI_AliHelperPID";   _beta_p_POI_AliHelperPID = createHisto2F(name,name, 500, 0, 5,  100, 0.1, 1.1,  "p", "beta","counts");
      name = "beta_p_AliHelperPID_no_Undefined";   _beta_p_AliHelperPID_no_Undefined = createHisto2F(name,name, 500, 0, 5,  100, 0.1, 1.1,  "p", "beta","counts");
      name = "inverse_beta_p_POI_AliHelperPID";   _inverse_beta_p_POI_AliHelperPID = createHisto2F(name,name, 500, 0, 5, 200, 0.6, 2.6,  "p", "1/#beta","counts");
      name = "inverse_beta_p_AliHelperPID_no_Undefined";   _inverse_beta_p_AliHelperPID_no_Undefined = createHisto2F(name,name, 500, 0, 5, 200, 0.6, 2.6,  "p", "1/#beta","counts");
      //name = "vZ_y_Pt_POI_AliHelperPID";        _vZ_y_Pt_POI_AliHelperPID = createHisto3F( name, name, _nBins_vertexZ, _min_vertexZ, _max_vertexZ, _nBins_eta_1, _min_eta_1, _max_eta_1, _nBins_pt_1, _min_pt_1, _max_pt_1, "zVertex", "y", "p_{T}" );
      //name = "vZ_y_eta_POI_AliHelperPID";        _vZ_y_eta_POI_AliHelperPID = createHisto3F( name, name, _nBins_vertexZ, _min_vertexZ, _max_vertexZ, _nBins_eta_1, _min_eta_1, _max_eta_1, _nBins_eta_1, _min_eta_1, _max_eta_1, "zVertex", "y", "#eta" );
      name = "y_Pt_AllCh_MCAODTruth";       _y_Pt_AllCh_MCAODTruth = createHisto2F( name, name, _nBins_eta_1, _min_eta_1, _max_eta_1, _nBins_pt_1, _min_pt_1, _max_pt_1, "y", "p_{T}", "counts" );
      name = "y_Pt_POI_MCAODTruth";        _y_Pt_POI_MCAODTruth = createHisto2F( name, name, _nBins_eta_1, _min_eta_1, _max_eta_1, _nBins_pt_1, _min_pt_1, _max_pt_1, "y", "p_{T}", "counts" );      
      name = n1Name + part_1_Name + vsPt;                         _n1_1_vsPt = createHisto1F( name, name, _nBins_pt_1, _min_pt_1, _max_pt_1, _title_pt_1, _title_AvgN_1 );
      name = n1Name + part_1_Name + vsZ + vsEtaPhi + vsPt;        _n1_1_vsZVsEtaVsPhiVsPt = createHisto3F( name, name, _nBins_vertexZ, _min_vertexZ, _max_vertexZ, _nBins_etaPhi_1, 0., double(_nBins_etaPhi_1), _nBins_pt_1, _min_pt_1, _max_pt_1, "zVertex", _title_etaPhi_1,  _title_pt_1);        
      name = n1Name + part_2_Name + vsPt;                         _n1_2_vsPt = createHisto1F( name, name, _nBins_pt_2, _min_pt_2, _max_pt_2, _title_pt_2, _title_AvgN_2 );
      name = n1Name + part_2_Name + vsZ + vsEtaPhi + vsPt;        _n1_2_vsZVsEtaVsPhiVsPt = createHisto3F( name, name, _nBins_vertexZ, _min_vertexZ, _max_vertexZ, _nBins_etaPhi_2, 0., double(_nBins_etaPhi_2), _nBins_pt_2, _min_pt_2, _max_pt_2, "zVertex", _title_etaPhi_2,  _title_pt_2);	
    }
  else
    {
      name = n1Name+part_1_Name+vsEtaPhi;       _n1_1_vsEtaVsPhi      = createHisto2F(name,name, _nBins_eta_1, _min_eta_1, _max_eta_1,  _nBins_phi_1, _min_phi_1, _max_phi_1,  _title_eta_1,  _title_phi_1,  _title_AvgN_1);
      name = s1ptName+part_1_Name+vsEtaPhi;     _s1pt_1_vsEtaVsPhi    = createHisto2F(name,name, _nBins_eta_1, _min_eta_1, _max_eta_1,  _nBins_phi_1, _min_phi_1, _max_phi_1,  _title_eta_1,  _title_phi_1,  _title_AvgSumPt_1);
      name = n1Name+part_1_Name+vsM;            _n1_1_vsM             = createProfile(name,name, _nBins_M4, _min_M4, _max_M4, _title_m4, _title_AvgN_1);
      name = s1ptName+part_1_Name+vsM;          _s1pt_1_vsM           = createProfile(name,name, _nBins_M4, _min_M4, _max_M4, _title_m4, _title_AvgSumPt_1);
      name = n1NwName+part_1_Name+vsM;          _n1Nw_1_vsM           = createProfile(name,name, _nBins_M4, _min_M4, _max_M4, _title_m4, _title_AvgN_1);
      name = s1ptNwName+part_1_Name+vsM;        _s1ptNw_1_vsM         = createProfile(name,name, _nBins_M4, _min_M4, _max_M4, _title_m4, _title_AvgSumPt_1);        
      name = n1Name+part_2_Name+vsEtaPhi;       _n1_2_vsEtaVsPhi      = createHisto2F(name,name, _nBins_eta_2, _min_eta_2, _max_eta_2,  _nBins_phi_2, _min_phi_2, _max_phi_2,  _title_eta_2,  _title_phi_2,  _title_AvgN_2);
      name = s1ptName+part_2_Name+vsEtaPhi;     _s1pt_2_vsEtaVsPhi    = createHisto2F(name,name, _nBins_eta_2, _min_eta_2, _max_eta_2,  _nBins_phi_2, _min_phi_2, _max_phi_2,  _title_eta_2,  _title_phi_2,  _title_AvgSumPt_2);
      name = n1Name+part_2_Name + vsM;          _n1_2_vsM             = createProfile(name,name, _nBins_M4, _min_M4, _max_M4, _title_m4, _title_AvgN_2);
      name = s1ptName+part_2_Name + vsM;        _s1pt_2_vsM           = createProfile(name,name, _nBins_M4, _min_M4, _max_M4, _title_m4, _title_AvgSumPt_2);
      name = n1NwName+part_2_Name+vsM;          _n1Nw_2_vsM           = createProfile(name,name, _nBins_M4, _min_M4, _max_M4, _title_m4, _title_AvgN_1);
      name = s1ptNwName+part_2_Name+vsM;        _s1ptNw_2_vsM         = createProfile(name,name, _nBins_M4, _min_M4, _max_M4, _title_m4, _title_AvgSumPt_1);        
      name = n2Name+pair_12_Name+vsEtaPhi;      _n2_12_vsEtaPhi       = createHisto1F(name,name, _nBins_etaPhi_12, 0.,        double(_nBins_etaPhi_12), _title_etaPhi_12, _title_AvgN2_12);
      name = s2PtPtName+pair_12_Name + vsEtaPhi;_s2PtPt_12_vsEtaPhi   = createHisto1F(name,name, _nBins_etaPhi_12, 0.,        double(_nBins_etaPhi_12), _title_etaPhi_12,  _title_AvgSumPtPt_12);
      name = s2PtNName+pair_12_Name + vsEtaPhi; _s2PtN_12_vsEtaPhi    = createHisto1F(name,name, _nBins_etaPhi_12, 0.,        double(_nBins_etaPhi_12), _title_etaPhi_12,  _title_AvgSumPtN_12);
      name = s2NPtName+pair_12_Name + vsEtaPhi; _s2NPt_12_vsEtaPhi    = createHisto1F(name,name, _nBins_etaPhi_12, 0.,        double(_nBins_etaPhi_12), _title_etaPhi_12,  _title_AvgNSumPt_12);
      name = n2Name+pair_12_Name+vsPtVsPt;      _n2_12_vsPtVsPt       = createHisto2F(name,name, _nBins_pt_1, _min_pt_1, _max_pt_1, _nBins_pt_2, _min_pt_2, _max_pt_2, _title_pt_1, _title_pt_2, _title_AvgN2_12);        
      name = n2Name+pair_12_Name + vsM;         _n2_12_vsM            = createProfile(name,name, _nBins_M4, _min_M4, _max_M4, _title_m4, _title_AvgN2_12);
      name = s2PtPtName+pair_12_Name + vsM;     _s2PtPt_12_vsM        = createProfile(name,name, _nBins_M4, _min_M4, _max_M4, _title_m4, _title_AvgSumPtPt_12);
      name = s2PtNName+pair_12_Name + vsM;      _s2PtN_12_vsM         = createProfile(name,name, _nBins_M4, _min_M4, _max_M4, _title_m4, _title_AvgSumPtN_12);
      name = s2NPtName+pair_12_Name + vsM;      _s2NPt_12_vsM         = createProfile(name,name, _nBins_M4, _min_M4, _max_M4, _title_m4, _title_AvgNSumPt_12);        
      name = n2NwName+pair_12_Name + vsM;       _n2Nw_12_vsM          = createProfile(name,name, _nBins_M4, _min_M4, _max_M4, _title_m4, _title_AvgN2_12);
      name = s2PtPtNwName+pair_12_Name + vsM;   _s2PtPtNw_12_vsM      = createProfile(name,name, _nBins_M4, _min_M4, _max_M4, _title_m4, _title_AvgSumPtPt_12);
      name = s2PtNNwName+pair_12_Name + vsM;    _s2PtNNw_12_vsM       = createProfile(name,name, _nBins_M4, _min_M4, _max_M4, _title_m4, _title_AvgSumPtN_12);
      name = s2NPtNwName+pair_12_Name + vsM;    _s2NPtNw_12_vsM       = createProfile(name,name, _nBins_M4, _min_M4, _max_M4, _title_m4, _title_AvgNSumPt_12);        
      name = "mInv";     _invMass     = createHisto1F(name,name, 50, 0.41, 0.55, "M_{inv}","counts");
      name = "mInvElec"; _invMassElec = createHisto1F(name,name, 500, 0., 1.000, "M_{inv}","counts");
    }
    
  AliInfo(" AliAnalysisTaskPIDBFDptDpt::createHistoHistograms() All Done");
}
//-----------------------//

void  AliAnalysisTaskPIDBFDptDpt::finalizeHistograms()
{
    
  AliInfo("AliAnalysisTaskPIDBFDptDpt::finalizeHistograms() starting");
  AliInfo(Form("CorrelationAnalyzers::finalizeHistograms()   _eventCount : %d",int(_eventCount)));
  if (_singlesOnly)
    {
      if (_sameFilter)
        {
	  fillHistoWithArray(_n1_1_vsPt,              __n1_1_vsPt,        _nBins_pt_1);
	  fillHistoWithArray(_n1_1_vsZVsEtaVsPhiVsPt, __n1_1_vsZEtaPhiPt, _nBins_vertexZ, _nBins_etaPhi_1, _nBins_pt_1);
	  fillHistoWithArray(_n1_2_vsPt,              __n1_1_vsPt,        _nBins_pt_1);
	  fillHistoWithArray(_n1_2_vsZVsEtaVsPhiVsPt, __n1_1_vsZEtaPhiPt, _nBins_vertexZ, _nBins_etaPhi_1, _nBins_pt_1);
        }
      else
        {
	  fillHistoWithArray(_n1_1_vsPt,              __n1_1_vsPt,        _nBins_pt_1);
	  fillHistoWithArray(_n1_1_vsZVsEtaVsPhiVsPt, __n1_1_vsZEtaPhiPt, _nBins_vertexZ, _nBins_etaPhi_1, _nBins_pt_1);
	  fillHistoWithArray(_n1_2_vsPt,              __n1_2_vsPt,        _nBins_pt_2);
	  fillHistoWithArray(_n1_2_vsZVsEtaVsPhiVsPt, __n1_2_vsZEtaPhiPt, _nBins_vertexZ, _nBins_etaPhi_2, _nBins_pt_2);
        }
    }
  else
    {
      if (_sameFilter)
        {
	  fillHistoWithArray(_n1_1_vsEtaVsPhi,        __n1_1_vsEtaPhi,    _nBins_eta_1,   _nBins_phi_1);
	  fillHistoWithArray(_s1pt_1_vsEtaVsPhi,      __s1pt_1_vsEtaPhi,  _nBins_eta_1,   _nBins_phi_1);
	  fillHistoWithArray(_n1_2_vsEtaVsPhi,        __n1_1_vsEtaPhi,    _nBins_eta_1,   _nBins_phi_1);
	  fillHistoWithArray(_s1pt_2_vsEtaVsPhi,      __s1pt_1_vsEtaPhi,  _nBins_eta_1,   _nBins_phi_1);
        }
      else
        {
	  fillHistoWithArray(_n1_1_vsEtaVsPhi,        __n1_1_vsEtaPhi,    _nBins_eta_1,   _nBins_phi_1);
	  fillHistoWithArray(_s1pt_1_vsEtaVsPhi,      __s1pt_1_vsEtaPhi,  _nBins_eta_1,   _nBins_phi_1);
	  fillHistoWithArray(_n1_2_vsEtaVsPhi,        __n1_2_vsEtaPhi,    _nBins_eta_2,   _nBins_phi_2);
	  fillHistoWithArray(_s1pt_2_vsEtaVsPhi,      __s1pt_2_vsEtaPhi,  _nBins_eta_2,   _nBins_phi_2);
        }
      fillHistoWithArray(_n2_12_vsEtaPhi,     __n2_12_vsEtaPhi,     _nBins_etaPhi_12);
      fillHistoWithArray(_s2PtPt_12_vsEtaPhi, __s2ptpt_12_vsEtaPhi, _nBins_etaPhi_12);
      fillHistoWithArray(_s2PtN_12_vsEtaPhi,  __s2PtN_12_vsEtaPhi,  _nBins_etaPhi_12);
      fillHistoWithArray(_s2NPt_12_vsEtaPhi,  __s2NPt_12_vsEtaPhi,  _nBins_etaPhi_12);
      fillHistoWithArray(_n2_12_vsPtVsPt,     __n2_12_vsPtPt,       _nBins_pt_1,    _nBins_pt_2);
        
    }
  AliInfo("AliAnalysisTaskPIDBFDptDpt::finalizeHistograms()  Done ");
}
//--------------//


void  AliAnalysisTaskPIDBFDptDpt::UserExec(Option_t */*option*/)
{    
  int    k1,k2;
  int    iPhi, iEta, iEtaPhi, iPt, charge;
  float  q, phi, pt, eta, y, y_direct, mass, corr, corrPt, px, py, pz, dedx,p,l, timeTOF, beta, t0, msquare; // Au-Au put p here to make _dedx_p and _beta_p plots. 
  int    ij;
  int    id_1, q_1, iEtaPhi_1, iPt_1;
  float  pt_1, px_1, py_1, pz_1, corr_1, dedx_1;
  int    id_2, q_2, iEtaPhi_2, iPt_2;
  float  pt_2, px_2, py_2, pz_2, corr_2, dedx_2;
  float  ptpt;
  int    iVertex, iVertexP1, iVertexP2;
  int    iZEtaPhiPt;
  float  massElecSq = 1.94797849000000016e-02;
  //double b[2];
  //double bCov[3];
  const  AliAODVertex*	vertex;
  int    nClus;
  bool   bitOK;
  const float mpion   = 0.139570; // GeV/c2
  const float mkaon   = 0.493677; // GeV/c2
  const float mproton = 0.938272; // GeV/c2
  Double_t c = TMath::C() * 1.E-9;// m/ns
    
  float slope_vZ = 0;
  float vZ_bin_center = 0;
  int iVertexZ_plus1 = 0;
  int iVertexZ_minus1 = 0;
  int iVertexP1_vZplus1 = 0;
  int iVertexP1_vZminus1 = 0;
  int iVertexP2_vZplus1 = 0;
  int iVertexP2_vZminus1 = 0;
  int iZEtaPhiPt_vZplus1 = 0;
  int iZEtaPhiPt_vZminus1 = 0;

  float slope_y = 0;
  float y_bin_center = 0;
  int iEta_plus1 = 0;
  int iEta_minus1 = 0;
  int iEtaPhi_Etaplus1 = 0;
  int iEtaPhi_Etaminus1 = 0;
  int iZEtaPhiPt_Etaplus1 = 0;
  int iZEtaPhiPt_Etaminus1 = 0;

  
  AliAnalysisManager* manager = AliAnalysisManager::GetAnalysisManager();
  if ( !manager ) { return; }

  
  AliAODInputHandler* inputHandler = dynamic_cast<AliAODInputHandler*> (manager->GetInputEventHandler());
  if ( !inputHandler ) { return; }

  
  fAODEvent = dynamic_cast<AliAODEvent*>(InputEvent()); // create pointer to event
    
  if ( !fAODEvent ) { return; }

  
  fPIDResponse = inputHandler -> GetPIDResponse();
  if (!fPIDResponse){
    AliFatal("This Task needs the PID response attached to the inputHandler");
    return;
  }

  
  _eventCount++;
    
  if ( _eventAccounting )  _eventAccounting -> Fill( 0 );
  else   return;

  //reset single particle counters
  k1 = k2 = 0;
  __n1_1 = __n1_2 = __s1pt_1 = __s1pt_2 = __n1Nw_1 = __n1Nw_2 = __s1ptNw_1 = __s1ptNw_2 = 0;
    
  float v0Centr  = -999.;
  float v0ACentr = -999.;
  float v0CCentr = -999.;
  float trkCentr = -999.;
  float spdCentr = -999.;
    
  float vertexX  = -999;
  float vertexY  = -999;
  float vertexZ  = -999;
  //float vertexXY = -999;
  float dcaZ     = -999;
  float dcaXY    = -999;
  float centrality = -999;
    
  if( fAODEvent )
    {
      
      if( fSystemType == "PbPb" || fSystemType == "pPb" )
	{
	  AliCentrality * centralityObject =  ( ( AliVAODHeader * ) fAODEvent -> GetHeader() ) -> GetCentralityP();
	  	
	  if ( centralityObject )
	    {
	      v0Centr  = centralityObject->GetCentralityPercentile("V0M");
	      v0ACentr = centralityObject->GetCentralityPercentile("V0A");
	      v0CCentr = centralityObject->GetCentralityPercentile("V0C");
	      trkCentr = centralityObject->GetCentralityPercentile("TRK");
	      spdCentr = centralityObject->GetCentralityPercentile("CL1");   
	    }
	}

      else if ( fSystemType == "pp_V0A_kMB" || fSystemType == "pp_V0C_kMB" || fSystemType == "pp_V0_kMB" )
	{ 
	  AliMultSelection *multSelection = (AliMultSelection*) fAODEvent->FindListObject("MultSelection");
	  //If get this warning, please check that the AliMultSelectionTask actually ran (before your task)
	  if (!multSelection)    AliWarning("MultSelection not found in input event");
	  else{
	    v0Centr  = multSelection->GetMultiplicityPercentile("V0M");
	    v0ACentr = multSelection->GetMultiplicityPercentile("V0A");
	    v0CCentr = multSelection->GetMultiplicityPercentile("V0C");
	  }
        }
      else if ( fSystemType == "pp_V0A_kMB_Utils" || fSystemType == "pp_V0C_kMB_Utils" || fSystemType == "pp_V0_kMB_Utils" )
	{
	  //remove Pile-Up events
	  fUtils->SetUseMVPlpSelection(kTRUE);
	  fUtils->SetUseOutOfBunchPileUp(kTRUE);
	  if(fUtils->IsPileUpEvent(fAODEvent))   return;
	  _eventAccounting -> Fill( 1 ); // number of events ( no pile-up ) 
	  
	  v0Centr  = fUtils->GetMultiplicityPercentile(fAODEvent,"V0MEq");
	  v0ACentr = fUtils->GetMultiplicityPercentile(fAODEvent,"V0AEq");
	  v0CCentr = fUtils->GetMultiplicityPercentile(fAODEvent,"V0CEq");
        }

      
      _nTracks  = fAODEvent -> GetNumberOfTracks();
        
      _mult3    = _nTracks;
      _mult4    = v0Centr;
      _mult5    = trkCentr;
      _mult6    = spdCentr;
      _mult7    = v0ACentr;
      _mult8    = v0CCentr;
      _field    = fAODEvent -> GetMagneticField();


      switch ( _centralityMethod )
        {
	case 0: centrality = _mult0;  break;
	case 1: centrality = _mult1;  break;
	case 2: centrality = _mult2;  break;
	case 3: centrality = _mult3;  break;
	case 4: centrality = _mult4;  break;
	case 5: centrality = _mult5;  break;
	case 6: centrality = _mult6;  break;
	case 7: centrality = _mult7;  break;
	case 8: centrality = _mult8;  break;
        }

      
      if      ( fSystemType == "PbPb" )
	{ if  ( centrality < _centralityMin || centrality > _centralityMax || fabs( v0Centr - trkCentr ) > 5.0 )  return; }
      else if ( fSystemType == "pPb" || fSystemType == "pp" || fSystemType == "pp_V0A_kMB" || fSystemType == "pp_V0_kMB" || fSystemType == "pp_V0C_kMB" || fSystemType == "pp_V0A_kMB_Utils" || fSystemType == "pp_V0_kMB_Utils" || fSystemType == "pp_V0C_kMB_Utils" )
	{ if  ( centrality < _centralityMin || centrality > _centralityMax )  return; }
      else    return;
	
      _eventAccounting -> Fill( 2 ); // count all events with right centrality
        
      // filter on z and xy vertex
      vertex = ( AliAODVertex * ) fAODEvent -> GetPrimaryVertex();
        
      if( vertex )
        {
	  Double32_t fCov[6];
	  vertex->GetCovarianceMatrix(fCov);
	  if(vertex->GetNContributors() > 0)
            {
	      if(fCov[5] != 0)
                {
		  vertexX = vertex->GetX();
		  vertexY = vertex->GetY();
		  vertexZ = vertex->GetZ();
                    
		  if( TMath::Abs( vertexZ ) > _max_vertexZ )	return;  // Z-Vertex Cut
                }
            }
        }
        
      iVertex = int( ( vertexZ - _min_vertexZ ) / _width_vertexZ );
      iVertexP1 = iVertex*_nBins_etaPhiPt_1;
      iVertexP2 = iVertex*_nBins_etaPhiPt_2;
      if ( iVertex<0 || iVertex >= _nBins_vertexZ )
        {
	  AliError("AliAnalysisTaskPIDBFDptDpt::Exec(Option_t *option) iVertex<0 || iVertex>=_nBins_vertexZ ");
	  return;
        }
	
      _eventAccounting -> Fill( 3 ); // count all events with right centrality & vertexZ

      //=========================================================================================================
      //*********************************************************************************************************            
      // "RealData" & "MCAODreco" share this same piece of code; if needed, seperate these 2 parts later 
      // "RealData" -- for a real dataset, e.g. LHC10h
      // "MCAODreco" -- for a MC reconstructed production dataset, e.g. Hijing_PbPb_LHC10h
      if ( fAnalysisType == "RealData" || fAnalysisType == "MCAODreco" )
	{
	  //Track Loop starts here
	  for (int iTrack = 0; iTrack < _nTracks; iTrack++ )
	    {
	      AliAODTrack * t = dynamic_cast<AliAODTrack *>( fAODEvent -> GetTrack( iTrack ) );
	      if ( !t ) {
                AliError( Form( "Could not receive track %d", iTrack ) );
                continue;
	      }

	      bitOK  = t -> TestFilterBit( _trackFilterBit );
	      if ( !bitOK ) continue;
            
	      q      = t -> Charge();
	      charge = int( q );
	      phi    = t -> Phi();
	      p      = t -> P(); // momentum magnitude
	      pt     = t -> Pt();
	      px     = t -> Px();
	      py     = t -> Py();
	      pz     = t -> Pz();
	      eta    = t -> Eta();
	      if ( fAnalysisType == "RealData" )    dedx = t -> GetTPCsignal();
	      else if ( fAnalysisType == "MCAODreco" )  dedx = fPIDResponse -> GetTPCsignalTunedOnData( t );
	    
	      // QA for all the particles in the event
	      if ( _singlesOnly )
		{
		  _etadis_before_any_cuts -> Fill( eta );
		  _phidis_before_any_cuts -> Fill( phi );

		  CheckTOF( t );
		  if ( fHasTOFPID )
		    {
		      _t0_1d -> Fill( fPIDResponse->GetTOFResponse().GetStartTime(t->P()) );
		      if ( fAnalysisType == "RealData" )
			{
			  _timeTOF_1d -> Fill( t->GetTOFsignal() );
			  _realTOF_1d -> Fill( t->GetTOFsignal() - fPIDResponse->GetTOFResponse().GetStartTime(t->P()) );
			}
		      else if ( fAnalysisType == "MCAODreco" )
			{
			  _timeTOF_1d -> Fill( t->GetTOFsignalTunedOnData() );
			  _realTOF_1d -> Fill( t->GetTOFsignalTunedOnData() - fPIDResponse->GetTOFResponse().GetStartTime(t->P()) );
			}
		      _trackLength -> Fill( fPIDResponse->GetTOFResponse().GetExpectedSignal(t, AliPID::kElectron)*1E-3*c );
		      _trackLength_GetIntegratedLength -> Fill( t -> GetIntegratedLength() );
		      _msquare_p -> Fill( p, massSquareCalculation(t) );
		      _beta_p -> Fill( p, TOFBetaCalculation(t) );
		      _nSigmaTOF_p -> Fill( p, fPIDResponse->NumberOfSigmasTOF(t,AliPID::kPion) );
		      _inverse_beta_p -> Fill( p, 1./TOFBetaCalculation(t)  );
		    }
		  else  _dedx_p -> Fill( p, dedx );
		}
	    
	      // Kinematics cuts begins:
	      if( charge == 0 ) continue;
	    
	      Double_t pos[3];
	      t -> GetXYZ(pos);
	      Double_t DCAX = pos[0] - vertexX;
	      Double_t DCAY = pos[1] - vertexY;
	      Double_t DCAZ = pos[2] - vertexZ;             
	      Double_t DCAXY = TMath::Sqrt((DCAX*DCAX) + (DCAY*DCAY));             
	      if (DCAZ     <  _dcaZMin ||
		  DCAZ     >  _dcaZMax ||
		  DCAXY    >  _dcaXYMax ) continue;       

	      nClus = t -> GetTPCNcls();
	      if ( nClus < _nClusterMin ) continue; // Kinematics cuts ends.
	    
	      if ( _singlesOnly )    _Ncluster1 -> Fill( nClus );
	      //_Ncluster2->Fill(nClus);   //_Ncluster2 same with _Ncluster1

	      //******************************************************************************************************************************************************************
	      if ( PIDparticle )
		{
		  Int_t IDrec = TellParticleSpecies( t ); //returns 0, 1, 2 for Pion, Kaon, Proton, respectively.


		  // QA for all identified hadrons
		  if ( IDrec != 999 )
		    {
		      if ( _singlesOnly )
			{
			  CheckTOF( t );
			  if ( fHasTOFPID  && ( pt > ptTOFlowerBoundary ) && ( pt <= ptUpperLimit ) )
			    {
			      _beta_p_AliHelperPID_no_Undefined -> Fill( p, TOFBetaCalculation(t) );
			      _inverse_beta_p_AliHelperPID_no_Undefined -> Fill( p, 1./TOFBetaCalculation(t) );
			      _msquare_p_AliHelperPID_no_Undefined -> Fill( p, massSquareCalculation(t) );
			    }
			  else if ( ( pt >= _min_pt_1 ) && ( pt <= ptTOFlowerBoundary ) )   _dedx_p_AliHelperPID_no_Undefined -> Fill( p, dedx );
			}
		    }

		
		  if ( IDrec != particleSpecies ) continue; // select POI
		
		  // This block of code used to get PdgCode for MC reco tracks
		  if ( fAnalysisType == "MCAODreco" && NoContamination == 1 )
		    {
		      Int_t lab =  t -> GetLabel();
		      //cout << "step 1 Au-Au: particle Label: " << lab << endl;
		      TClonesArray * arr = dynamic_cast<TClonesArray*>( fAODEvent -> FindListObject( "mcparticles" ) );
		      if( !arr ) continue;
		      //cout << "step 2 Au-Au: arr: " << arr->GetEntriesFast() << endl;
		      AliAODMCParticle * particle = ( AliAODMCParticle * ) arr -> At( TMath::Abs(lab) );
		      //cout << "step 3 Au-Au: particle ID: " << particle -> GetPdgCode() << endl;
		      if ( TMath::Abs(particle->GetPdgCode()) != 211 ) continue;
		      //cout << "step 4 Au-Au: pion ID: " << particle -> GetPdgCode() << endl;
		    }
		
		  if ( particleSpecies == 0 )  mass = mpion;
		  if ( particleSpecies == 1 )  mass = mkaon;
		  if ( particleSpecies == 2 )  mass = mproton;
		  y = log( ( sqrt(mass*mass + pt*pt*cosh(eta)*cosh(eta)) + pt*sinh(eta) ) / sqrt(mass*mass + pt*pt) ); // convert eta to y // CAVEAT: y is not right for non-POI @ this step
		  if( y < _min_eta_1 || y > _max_eta_1 ) continue;

		  //Filling QA plots ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		  // QA for POI
		  if ( _singlesOnly )
		    {
		      _dcaz                      -> Fill( DCAZ );
		      _dcaxy                     -> Fill( DCAXY );
		      _etadis_POI_AliHelperPID   -> Fill( eta );    //Eta dist. for POI distribution after AliHelperPID cuts
		      _ydis_POI_AliHelperPID     -> Fill( y );
		      _phidis_POI_AliHelperPID   -> Fill( phi );
		      //_vZ_y_Pt_POI_AliHelperPID  -> Fill( vertexZ, y, pt );
		      //_vZ_y_eta_POI_AliHelperPID -> Fill( vertexZ, y, eta );
		    
		      CheckTOF( t );
		      if ( fHasTOFPID  && ( pt > ptTOFlowerBoundary ) && ( pt <= ptUpperLimit ) )
			{
			  _t0_1d_POI -> Fill( fPIDResponse->GetTOFResponse().GetStartTime(t->P()) );
			  if ( fAnalysisType == "RealData" )
			    {
			      _timeTOF_1d_POI -> Fill( t->GetTOFsignal() );
			      _realTOF_1d_POI -> Fill( t->GetTOFsignal() - fPIDResponse->GetTOFResponse().GetStartTime(t->P()) );
			    }
			  else if ( fAnalysisType == "MCAODreco" )
			    {
			      _timeTOF_1d_POI -> Fill( t->GetTOFsignalTunedOnData() );
			      _realTOF_1d_POI -> Fill( t->GetTOFsignalTunedOnData() - fPIDResponse->GetTOFResponse().GetStartTime(t->P()) );
			    }
			  _trackLength_POI -> Fill( fPIDResponse->GetTOFResponse().GetExpectedSignal(t, AliPID::kElectron)*1E-3*c );
			  _trackLength_GetIntegratedLength_POI -> Fill( t -> GetIntegratedLength() );
			  _nSigmaTOF_p_POI -> Fill( p, fPIDResponse->NumberOfSigmasTOF(t,AliPID::kPion) );
			  _msquare_p_POI_AliHelperPID -> Fill( p, massSquareCalculation(t) );
			  _beta_p_POI_AliHelperPID -> Fill( p, TOFBetaCalculation(t) );
			  _inverse_beta_p_POI_AliHelperPID -> Fill( p, 1./TOFBetaCalculation(t) );
			}
		      else if ( ( pt >= _min_pt_1 ) && ( pt <= ptTOFlowerBoundary ) )   _dedx_p_POI_AliHelperPID -> Fill( p, dedx );
		    }
		}
	      else //all charged particles 
		{
		  if( pt < _min_pt_1 || pt > ptUpperLimit ) continue;
		  CalculateTPCNSigmasElectron( t );
		  Double_t nsigmaElectron = 999.;
		  nsigmaElectron =  TMath::Abs( fnsigmas[3][0] ); //Electron_TPC
		  if( ( pt >= _min_pt_1 ) && ( pt <= ptTOFlowerBoundary ) && ( nsigmaElectron < electronNSigmaVeto ) )    continue;  // reject TPC region electrons
		  if( eta < _min_eta_1 || eta > _max_eta_1 ) continue;

		  // QA for POI
		  if ( _singlesOnly )
		    {
		      _dcaz                      -> Fill( DCAZ );
		      _dcaxy                     -> Fill( DCAXY );
		      _etadis_POI_AliHelperPID   -> Fill( eta );
		      _phidis_POI_AliHelperPID   -> Fill( phi );
		      _dedx_p_POI_AliHelperPID   -> Fill( p, dedx ); 
		    }
		}   
	      //Filling QA plots ends ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	    
	      if ( _useRapidity )  eta = y;  //switch from eta to y
	      //*************************************************************************************************************************************************************
	    
	      //Particle 1
	      if ( _requestedCharge_1 == charge )
		//&& dedx >=  _dedxMin && dedx < _dedxMax)
		{
		  iPhi   = int( phi/_width_phi_1);
                
		  if (iPhi<0 || iPhi>=_nBins_phi_1 )
		    {
		      AliWarning("AliAnalysisTaskPIDBFDptDpt::analyze() iPhi<0 || iPhi>=_nBins_phi_1");
		      return;
		    }
                
		  iEta    = int((eta-_min_eta_1)/_width_eta_1);
		  if (iEta<0 || iEta>=_nBins_eta_1)
		    {
		      AliWarning(Form("AliAnalysisTaskPIDBFDptDpt::analyze(AliceEvent * event) Mismatched iEta: %d", iEta));
		      continue;
		    }
		  iPt     = int((pt -_min_pt_1 )/_width_pt_1 );
		  if (iPt<0  || iPt >=_nBins_pt_1)
		    {
		      AliWarning(Form("AliAnalysisTaskPIDBFDptDpt::analyze(AliceEvent * event) Mismatched iPt: %d",iPt));
		      continue;
		    }
		  iEtaPhi = iEta * _nBins_phi_1 + iPhi;
		  iZEtaPhiPt = iVertexP1 + iEtaPhi * _nBins_pt_1 + iPt;

		  if ( _correctionWeight_1 )   corr = _correctionWeight_1[ iZEtaPhiPt ];
		  else   corr = 1;				

		  if ( iZEtaPhiPt < 0 || iZEtaPhiPt >= _nBins_zEtaPhiPt_1 )
		    {
		      AliWarning("AliAnalysisTaskPIDBFDptDpt::analyze(AliceEvent * event) iZEtaPhiPt<0 || iZEtaPhiPt>=_nBins_zEtaPhiPt_1");
		      continue;
		    }                
                
		  if ( _singlesOnly )
		    {
		      __n1_1_vsPt[iPt]               += corr;          //cout << "step 15" << endl;
		      __n1_1_vsZEtaPhiPt[iZEtaPhiPt] += corr;       //cout << "step 12" << endl;       
		    }
		  else
		    {
		      corrPt                      = corr*pt;
		      _id_1[k1]                   = iTrack;
		      _charge_1[k1]               = charge;
		      _iEtaPhi_1[k1]              = iEtaPhi;
		      _iPt_1[k1]                  = iPt;
		      _pt_1[k1]                   = pt;
		      _px_1[k1]                   = px;
		      _py_1[k1]                   = py;
		      _pz_1[k1]                   = pz;
		      _correction_1[k1]           = corr;
		      __n1_1                      += corr;
		      __n1_1_vsEtaPhi[iEtaPhi]    += corr;
		      __s1pt_1                    += corrPt;
		      __s1pt_1_vsEtaPhi[iEtaPhi]  += corrPt;
		      __n1Nw_1                    += 1;
		      __s1ptNw_1                  += pt;
		      ++k1;
		      if (k1>=arraySize)
			{
			  AliError(Form("AliAnalysisTaskPIDBFDptDpt::analyze(AliceEvent * event) k1 >=arraySize; arraySize: %d",arraySize));
			  return;
			}
		    }
		}//if ( _requestedCharge_1 == charge )
            
	      if ( !_sameFilter && _requestedCharge_2 == charge )
		//&& dedx >=  _dedxMin && dedx < _dedxMax)                
		{
		  iPhi   = int( phi/_width_phi_2);
                
		  if (iPhi<0 || iPhi>=_nBins_phi_2 )
		    {
		      AliWarning("AliAnalysisTaskPIDBFDptDpt::analyze() iPhi<0 || iPhi>=_nBins_phi_1");
		      return;
		    }
                
		  iEta    = int((eta-_min_eta_2)/_width_eta_2);
		  if (iEta<0 || iEta>=_nBins_eta_2)
		    {
		      AliWarning(Form("AliAnalysisTaskPIDBFDptDpt::analyze(AliceEvent * event) Mismatched iEta: %d", iEta));
		      continue;
		    }
		  iPt     = int((pt -_min_pt_2 )/_width_pt_2 );
		  if (iPt<0  || iPt >=_nBins_pt_2)
		    {
		      AliWarning(Form("AliAnalysisTaskPIDBFDptDpt::analyze(AliceEvent * event) Mismatched iPt: %d",iPt));
		      continue;
		    }
                
		  iEtaPhi = iEta*_nBins_phi_2+iPhi;
		  iZEtaPhiPt = iVertexP2 + iEtaPhi*_nBins_pt_2 + iPt;
		  if (iZEtaPhiPt<0 || iZEtaPhiPt>=_nBins_zEtaPhiPt_2)
		    {
		      AliWarning("AliAnalysisTaskPIDBFDptDpt::analyze(AliceEvent * event) iZEtaPhiPt<0 || iZEtaPhiPt>=_nBins_zEtaPhiPt_2");
		      continue;
		    }

		  if ( _correctionWeight_2 )   corr = _correctionWeight_2[ iZEtaPhiPt ];
		  else   corr = 1;		

		  if (_singlesOnly)
		    {
		      __n1_2_vsPt[iPt]               += corr;          //cout << "step 15" << endl;
		      __n1_2_vsZEtaPhiPt[iZEtaPhiPt] += corr;       //cout << "step 12" << endl;
		    }
		  else
		    {
		      corrPt                      = corr*pt;
		      _id_2[k2]                   = iTrack;         //cout << "step 1" << endl;
		      _charge_2[k2]               = charge;         //cout << "step 2" << endl;
		      _iEtaPhi_2[k2]              = iEtaPhi;        //cout << "step 3" << endl;
		      _iPt_2[k2]                  = iPt;            //cout << "step 4" << endl;
		      _pt_2[k2]                   = pt;             //cout << "step 5" << endl;
		      _px_2[k2]                   = px;             //cout << "step 6" << endl;
		      _py_2[k2]                   = py;             //cout << "step 7" << endl;
		      _pz_2[k2]                   = pz;             //cout << "step 8" << endl;
		      _correction_2[k2]           = corr;           //cout << "step 9" << endl;
		      __n1_2                      += corr;          //cout << "step 10" << endl;
		      __s1pt_2                    += corrPt;        //cout << "step 13" << endl;
		      __n1Nw_2                    += 1;
		      __n1_2_vsEtaPhi[iEtaPhi]    += corr;          //cout << "step 11" << endl;
		      __s1pt_2_vsEtaPhi[iEtaPhi]  += corrPt;        //cout << "step 14" << endl;
		      __s1ptNw_2                  += pt;
		      ++k2;
		      if (k2>=arraySize)
			{
			  AliWarning(Form("-W- k2 >=arraySize; arraySize: %d",arraySize));
			  return;
			}
		    }
		} //if ( !_sameFilter && _requestedCharge_2 == charge )
	    } //Track Loop ends here //for (int iTrack=0; iTrack< _nTracks; iTrack++)
	} //end of "if ( fAnalysisType == "RealData" || fAnalysisType == "MCAODreco" )"

      //=========================================================================================================
      //*********************************************************************************************************           
      // MC AOD Truth
      else if ( fAnalysisType == "MCAOD" )
	{
	  AliMCEvent * mcEvent = MCEvent();
	  _nTracks = mcEvent -> GetNumberOfTracks();
	
	  //Track Loop starts here
	  for ( int iTrack = 0; iTrack < _nTracks; iTrack++ )
	    {
	      AliAODMCParticle * t = ( AliAODMCParticle * ) mcEvent -> GetTrack( iTrack );
	  
	      if ( !t ) {
                AliError(Form("Could not receive track %d", iTrack));
                continue;
	      }

	      if( !t -> IsPhysicalPrimary() ) continue;
	    
	      q      = t -> Charge();

	      //cout << "step 1 Au-Au: q = " << q << endl; 
	    
	      charge = int( q/3. ); // particle charges are 3s in HIJING truth

	      //cout << "step 2 Au-Au: charge = " << charge << endl;
	    
	      phi    = t -> Phi();
	      pt     = t -> Pt();
	      eta    = t -> Eta();
	      //y_direct = t -> Y(); // rapidity
	      y      = t -> Y();

	      if ( _singlesOnly )
		{
		  _etadis_before_any_cuts -> Fill( eta );
		  _phidis_before_any_cuts -> Fill( phi );
		}
	    
	      /*
		const float mpion = 0.1395701835; // GeV/c2
		const float mkaon = 0.493667; // GeV/c2
		const float mproton = 0.93827204621; // GeV/c2
		if ( particleSpecies == 0 )  mass = mpion;
		if ( particleSpecies == 1 )  mass = mkaon;
		if ( particleSpecies == 2 )  mass = mproton;
		y = log( ( sqrt(mass*mass + pt*pt*cosh(eta)*cosh(eta)) + pt*sinh(eta) ) / sqrt(mass*mass + pt*pt) ); // convert eta to y
	    
		//check if 2 ways of getting rapidity give same results
		if( TMath::Abs( t -> GetPdgCode() ) == 321  )
		{
		cout << " y = " << y << endl;
		cout << " y_direct = " << y_direct << endl;
		}
	      */	    

	      // Kinematics cuts:
	      if( charge == 0 ) continue;
	      if( pt < _min_pt_1 || pt > ptUpperLimit ) continue;         
	      if( y < _min_eta_1 || y > _max_eta_1 ) continue;

	      /*
		int pdg = t -> GetPdgCode();
		// Compare to http://pdg.lbl.gov/2007/reviews/montecarlorpp.pdf
		// Important ones:
		// proton: +/- 2212
		// neutron 2212 (maybe +/- for anti-neutrons)
		// charged kaon: +/- 321
		// charged pion: +/- 211
		// electron: +/- 11
		// muon: +/- 13
		// If you get any other number, double-check!
		double mass = TParticlePDG::Mass( pdg );
	      */    
	    
	      // fill track QA histograms
	      if ( _singlesOnly )  _y_Pt_AllCh_MCAODTruth -> Fill( y, pt ); // All Charged particles	    

	      if( particleSpecies == 0 )
		{ if( TMath::Abs( t -> GetPdgCode() ) != 211  )  continue; }
	      else if( particleSpecies == 1 )
		{ if( TMath::Abs( t -> GetPdgCode() ) != 321  )  continue; }
	      else if( particleSpecies == 2 )
		{ if( TMath::Abs( t -> GetPdgCode() ) != 2212 )  continue; }
	      else   return;

	      if ( _singlesOnly )
		{
		  _etadis_POI_AliHelperPID   -> Fill( eta );    //Eta dist. for POI distribution after AliHelperPID cuts
		  _ydis_POI_AliHelperPID     -> Fill( y );
		  _phidis_POI_AliHelperPID   -> Fill( phi );
		  _y_Pt_POI_MCAODTruth -> Fill( y, pt ); //POI
		}
	      //cout << "step 3 Au-Au: particle ID: " << t -> GetPdgCode() << endl;
	    
	      //Exclude resonances
	      if( fExcludeResonancesInMC )
		{
		  Int_t gMotherIndex = t -> GetMother();
		  if( gMotherIndex != -1 ) {
		    AliAODMCParticle * motherTrack = dynamic_cast<AliAODMCParticle *>( mcEvent -> GetTrack( gMotherIndex ) );
		    if( motherTrack ) {
		      Int_t pdgCodeOfMother = motherTrack -> GetPdgCode();
                       
		      if( pdgCodeOfMother == 311 || pdgCodeOfMother == -311 // K0
			  || pdgCodeOfMother == 310
			  || pdgCodeOfMother == 3122 || pdgCodeOfMother == -3122 // Lambda
			  || pdgCodeOfMother == 111 // pi0
			  || pdgCodeOfMother == 22 // photon
			  ) continue;
		    }
		  }
		}

	      //cout << "step 4 Au-Au" << endl;

	      /*
	      //Exclude electrons with PDG                                                              
	      if( fExcludeElectronsInMC ) {
	      if( TMath::Abs( t -> GetPdgCode() ) == 11 ) continue;
	      }
	      */	    

	      if ( _useRapidity )  eta = y;  //switch from eta to y	    

	      //cout << "step 5 Au-Au" << endl;
	    
	      //Particle 1
	      if ( _requestedCharge_1 == charge )
		//&& dedx >=  _dedxMin && dedx < _dedxMax)
		{
		  iPhi   = int( phi/_width_phi_1);
                
		  if (iPhi<0 || iPhi>=_nBins_phi_1 )
		    {
		      AliWarning("AliAnalysisTaskPIDBFDptDpt::analyze() iPhi<0 || iPhi>=_nBins_phi_1");
		      return;
		    }

		  //cout << "step 6 Au-Au" << endl;
		
		  iEta    = int((eta-_min_eta_1)/_width_eta_1);
		  if (iEta<0 || iEta>=_nBins_eta_1)
		    {
		      AliWarning(Form("AliAnalysisTaskPIDBFDptDpt::analyze(AliceEvent * event) Mismatched iEta: %d", iEta));
		      continue;
		    }

		  //cout << "step 7 Au-Au" << endl;
		
		  iPt     = int((pt -_min_pt_1 )/_width_pt_1 );
		  if (iPt<0  || iPt >=_nBins_pt_1)
		    {
		      AliWarning(Form("AliAnalysisTaskPIDBFDptDpt::analyze(AliceEvent * event) Mismatched iPt: %d",iPt));
		      continue;
		    }
		  iEtaPhi = iEta*_nBins_phi_1+iPhi;
		  iZEtaPhiPt = iVertexP1 + iEtaPhi*_nBins_pt_1 + iPt;
                
		  if (_correctionWeight_1)  corr = _correctionWeight_1[iZEtaPhiPt];
		  else  corr = 1;

		  //cout << "step 8 Au-Au" << endl;
		
		  if (iZEtaPhiPt<0 || iZEtaPhiPt>=_nBins_zEtaPhiPt_1)
		    {
		      AliWarning("AliAnalysisTaskPIDBFDptDpt::analyze(AliceEvent * event) iZEtaPhiPt<0 || iZEtaPhiPt>=_nBins_zEtaPhiPt_1");
		      continue;
		    }
                                
		  if (_singlesOnly)
		    {
		      __n1_1_vsPt[iPt]               += corr;       //cout << "step 9 Au-Au: __n1_1_vsPt[iPt] = " << __n1_1_vsPt[iPt] << endl;
		      __n1_1_vsZEtaPhiPt[iZEtaPhiPt] += corr;       //cout << "step 12" << endl;   
		    }
		  else
		    {
		      corrPt                      = corr*pt;
		      _id_1[k1]                   = iTrack;
		      _charge_1[k1]               = charge;
		      _iEtaPhi_1[k1]              = iEtaPhi;
		      _iPt_1[k1]                  = iPt;
		      _pt_1[k1]                   = pt;
		      _px_1[k1]                   = px;
		      _py_1[k1]                   = py;
		      _pz_1[k1]                   = pz;
		      _correction_1[k1]           = corr;
		      __n1_1                      += corr;
		      __n1_1_vsEtaPhi[iEtaPhi]    += corr;
		      __s1pt_1                    += corrPt;
		      __s1pt_1_vsEtaPhi[iEtaPhi]  += corrPt;
		      __n1Nw_1                    += 1;
		      __s1ptNw_1                  += pt;
		      ++k1;
		      if (k1>=arraySize)
			{
			  AliError(Form("AliAnalysisTaskPIDBFDptDpt::analyze(AliceEvent * event) k1 >=arraySize; arraySize: %d",arraySize));
			  return;
			}
		    }
		} // if ( _requestedCharge_1 == charge )
            
	      if ( !_sameFilter && _requestedCharge_2 == charge )
		//&& dedx >=  _dedxMin && dedx < _dedxMax)  
		{
		  iPhi   = int( phi/_width_phi_2);
                
		  if (iPhi<0 || iPhi>=_nBins_phi_2 )
		    {
		      AliWarning("AliAnalysisTaskPIDBFDptDpt::analyze() iPhi<0 || iPhi>=_nBins_phi_1");
		      return;
		    }
                
		  iEta    = int((eta-_min_eta_2)/_width_eta_2);
		  if (iEta<0 || iEta>=_nBins_eta_2)
		    {
		      AliWarning(Form("AliAnalysisTaskPIDBFDptDpt::analyze(AliceEvent * event) Mismatched iEta: %d", iEta));
		      continue;
		    }
		  iPt     = int((pt -_min_pt_2 )/_width_pt_2 );
		  if (iPt<0  || iPt >=_nBins_pt_2)
		    {
		      AliWarning(Form("AliAnalysisTaskPIDBFDptDpt::analyze(AliceEvent * event) Mismatched iPt: %d",iPt));
		      continue;
		    }
                
		  iEtaPhi = iEta*_nBins_phi_2+iPhi;
		  iZEtaPhiPt = iVertexP2 + iEtaPhi*_nBins_pt_2 + iPt;
		  if (iZEtaPhiPt<0 || iZEtaPhiPt>=_nBins_zEtaPhiPt_2)
		    {
		      AliWarning("AliAnalysisTaskPIDBFDptDpt::analyze(AliceEvent * event) iZEtaPhiPt<0 || iZEtaPhiPt>=_nBins_zEtaPhiPt_2");
		      continue;
		    }   
                
		  if (_correctionWeight_2)  corr = _correctionWeight_2[iZEtaPhiPt];
		  else corr = 1;
                
		  if (_singlesOnly)
		    {
		      __n1_2_vsPt[iPt]               += corr;       //cout << "step 15" << endl;
		      __n1_2_vsZEtaPhiPt[iZEtaPhiPt] += corr;       //cout << "step 12" << endl;
		    }
		  else
		    {
		      corrPt                      = corr*pt;
		      _id_2[k2]                   = iTrack;         //cout << "step 1" << endl;
		      _charge_2[k2]               = charge;         //cout << "step 2" << endl;
		      _iEtaPhi_2[k2]              = iEtaPhi;        //cout << "step 3" << endl;
		      _iPt_2[k2]                  = iPt;            //cout << "step 4" << endl;
		      _pt_2[k2]                   = pt;             //cout << "step 5" << endl;
		      _px_2[k2]                   = px;             //cout << "step 6" << endl;
		      _py_2[k2]                   = py;             //cout << "step 7" << endl;
		      _pz_2[k2]                   = pz;             //cout << "step 8" << endl;
		      _correction_2[k2]           = corr;           //cout << "step 9" << endl;
		      __n1_2                      += corr;          //cout << "step 10" << endl;
		      __s1pt_2                    += corrPt;        //cout << "step 13" << endl;
		      __n1Nw_2                    += 1;
		      __n1_2_vsEtaPhi[iEtaPhi]    += corr;          //cout << "step 11" << endl;
		      __s1pt_2_vsEtaPhi[iEtaPhi]  += corrPt;        //cout << "step 14" << endl;
		      __s1ptNw_2                  += pt;
		      ++k2;
		      if (k2>=arraySize)
			{
			  AliWarning(Form("-W- k2 >=arraySize; arraySize: %d",arraySize));
			  return;
			}
		    }
		} //if ( !_sameFilter && _requestedCharge_2 == charge )
	    } //Track Loop ends here //for (int iTrack=0; iTrack< _nTracks; iTrack++)
	} //end of "if ( fAnalysisType == "MCAOD" )" 
    } //if( fAODEvent )

  // Fill event QA histos
  _m0->Fill(_mult0);
  _m1->Fill(_mult1);
  _m2->Fill(_mult2);
  _m3->Fill(_mult3);
  _m4->Fill(_mult4);
  _m5->Fill(_mult5);
  _m6->Fill(_mult6);
  _m7->Fill(_mult7);
  _m8->Fill(_mult8);
  _vertexZ->Fill(vertexZ);
    
  if ( _singlesOnly )
    {
      // nothing to do here.
    }
  else
    {
      if (_sameFilter)
        {
	  _n1_1_vsM->Fill(centrality,      __n1_1);
	  _s1pt_1_vsM->Fill(centrality,    __s1pt_1);
	  _n1Nw_1_vsM->Fill(centrality,    __n1Nw_1);
	  _s1ptNw_1_vsM->Fill(centrality,  __s1ptNw_1);
	  _n1_2_vsM->Fill(centrality,      __n1_1);
	  _s1pt_2_vsM->Fill(centrality,    __s1pt_1);
	  _n1Nw_2_vsM->Fill(centrality,    __n1Nw_1);
	  _s1ptNw_2_vsM->Fill(centrality,  __s1ptNw_1);
	  // reset pair counters
	  __n2_12   = __s2ptpt_12   = __s2NPt_12    = __s2PtN_12    = 0;
	  __n2Nw_12 = __s2ptptNw_12 = __s2NPtNw_12  = __s2PtNNw_12  = 0;
	  if (_field>0)
            {
	      for (int i1=0; i1<k1; i1++)
                {
		  ////cout << "         i1:" << i1 << endl;
		  id_1      = _id_1[i1];           ////cout << "       id_1:" << id_1 << endl;
		  q_1       = _charge_1[i1];       ////cout << "        q_1:" << q_1 << endl;
		  iEtaPhi_1 = _iEtaPhi_1[i1];      ////cout << "  iEtaPhi_1:" << iEtaPhi_1 << endl;
		  iPt_1     = _iPt_1[i1];          ////cout << "      iPt_1:" << iPt_1 << endl;
		  corr_1    = _correction_1[i1];   ////cout << "     corr_1:" << corr_1 << endl;
		  pt_1      = _pt_1[i1];           ////cout << "       pt_1:" << pt_1 << endl;
		  dedx_1    = _dedx_1[i1];         ////cout << "     dedx_1:" << dedx_1 << endl;
		  //1 and 2
		  for (int i2=i1+1; i2<k1; i2++)
                    {
		      ////cout << "         i2:" << i2 << endl;
		      id_2      = _id_1[i2];              ////cout << "       id_2:" << id_2 << endl;
		      if (id_1!=id_2)
                        {
			  q_2       = _charge_1[i2];     ////cout << "        q_1:" << q_1 << endl;
			  iEtaPhi_2 = _iEtaPhi_1[i2];    ////cout << "  iEtaPhi_1:" << iEtaPhi_1 << endl;
			  iPt_2     = _iPt_1[i2];        ////cout << "      iPt_1:" << iPt_1 << endl;
			  corr_2    = _correction_1[i2]; ////cout << "     corr_1:" << corr_1 << endl;
			  pt_2      = _pt_1[i2];         ////cout << "       pt_1:" << pt_1 << endl;
			  dedx_2    = _dedx_1[i2];       ////cout << "     dedx_2:" << dedx_2 << endl;
			  corr      = corr_1*corr_2;
			  if (q_2>q_1 || (q_1>0 && q_2>0 && pt_2<=pt_1) || (q_1<0 && q_2<0 && pt_2>=pt_1))
                            {
			      ij = iEtaPhi_1*_nBins_etaPhi_1 + iEtaPhi_2;   ////cout << " ij:" << ij<< endl;
                            }
			  else // swap particles
                            {
			      ij = iEtaPhi_2*_nBins_etaPhi_1 + iEtaPhi_1;   ////cout << " ij:" << ij<< endl;
                            }
                            
			  __n2_12                  += corr;
			  __n2_12_vsEtaPhi[ij]     += corr;
			  ptpt                     = pt_1*pt_2;
			  __s2ptpt_12              += corr*ptpt;
			  __s2PtN_12               += corr*pt_1;
			  __s2NPt_12               += corr*pt_2;
			  __s2ptpt_12_vsEtaPhi[ij] += corr*ptpt;
			  __s2PtN_12_vsEtaPhi[ij]  += corr*pt_1;
			  __s2NPt_12_vsEtaPhi[ij]  += corr*pt_2;
			  __n2_12_vsPtPt[iPt_1*_nBins_pt_2 + iPt_2] += corr;
                            
			  __n2Nw_12                  += 1;
			  __s2ptptNw_12              += ptpt;
			  __s2PtNNw_12               += pt_1;
			  __s2NPtNw_12               += pt_2;
                            
                        }
                    } //i2
                } //i1
            }
	  else // field<0
            {
	      for (int i1=0; i1<k1; i1++)
                {
		  ////cout << "         i1:" << i1 << endl;
		  id_1      = _id_1[i1];           ////cout << "       id_1:" << id_1 << endl;
		  q_1       = _charge_1[i1];       ////cout << "        q_1:" << q_1 << endl;
		  iEtaPhi_1 = _iEtaPhi_1[i1];      ////cout << "  iEtaPhi_1:" << iEtaPhi_1 << endl;
		  iPt_1     = _iPt_1[i1];          ////cout << "      iPt_1:" << iPt_1 << endl;
		  corr_1    = _correction_1[i1];   ////cout << "     corr_1:" << corr_1 << endl;
		  pt_1      = _pt_1[i1];           ////cout << "       pt_1:" << pt_1 << endl;
		  dedx_1    = _dedx_1[i1];         ////cout << "     dedx_1:" << dedx_1 << endl;
		  //1 and 2
		  for (int i2=i1+1; i2<k1; i2++)
                    {
		      ////cout << "         i2:" << i2 << endl;
		      id_2      = _id_1[i2];              ////cout << "       id_2:" << id_2 << endl;
		      if (id_1!=id_2)
                        {
			  q_2       = _charge_1[i2];     ////cout << "        q_2:" << q_2 << endl;
			  iEtaPhi_2 = _iEtaPhi_1[i2];    ////cout << "  iEtaPhi_2:" << iEtaPhi_2 << endl;
			  iPt_2     = _iPt_1[i2];        ////cout << "      iPt_2:" << iPt_2 << endl;
			  corr_2    = _correction_1[i2]; ////cout << "     corr_2:" << corr_2 << endl;
			  pt_2      = _pt_1[i2];         ////cout << "       pt_2:" << pt_2 << endl;
			  dedx_2    = _dedx_1[i2];       ////cout << "     dedx_2:" << dedx_2 << endl;
			  corr      = corr_1*corr_2;
			  if ( q_2<q_1 || (q_1>0 && q_2>0 && pt_2>=pt_1) || (q_1<0 && q_2<0 && pt_2<=pt_1))
                            {
			      ij = iEtaPhi_1*_nBins_etaPhi_1 + iEtaPhi_2;   ////cout << " ij:" << ij<< endl;
                            }
			  else // swap particles
                            {
			      ij = iEtaPhi_2*_nBins_etaPhi_1 + iEtaPhi_1;   ////cout << " ij:" << ij<< endl;
                            }
                            
			  __n2_12                  += corr;
			  __n2_12_vsEtaPhi[ij]     += corr;
			  ptpt                     = pt_1*pt_2;
			  __s2ptpt_12              += corr*ptpt;
			  __s2PtN_12               += corr*pt_1;
			  __s2NPt_12               += corr*pt_2;
			  __s2ptpt_12_vsEtaPhi[ij] += corr*ptpt;
			  __s2PtN_12_vsEtaPhi[ij]  += corr*pt_1;
			  __s2NPt_12_vsEtaPhi[ij]  += corr*pt_2;
			  __n2_12_vsPtPt[iPt_1*_nBins_pt_2 + iPt_2] += corr;
                            
			  __n2Nw_12                  += 1;
			  __s2ptptNw_12              += ptpt;
			  __s2PtNNw_12               += pt_1;
			  __s2NPtNw_12               += pt_2;
                            
                        }
                    } //i2
                } //i1
            }
        }
      else  // filter 1 and 2 are different -- must do all particle pairs...
        {
	  _n1_1_vsM->Fill(centrality,      __n1_1);
	  _s1pt_1_vsM->Fill(centrality,    __s1pt_1);
	  _n1Nw_1_vsM->Fill(centrality,    __n1Nw_1);
	  _s1ptNw_1_vsM->Fill(centrality,  __s1ptNw_1);
	  _n1_2_vsM->Fill(centrality,      __n1_2);
	  _s1pt_2_vsM->Fill(centrality,    __s1pt_2);
	  _n1Nw_2_vsM->Fill(centrality,    __n1Nw_2);
	  _s1ptNw_2_vsM->Fill(centrality,  __s1ptNw_2);
	  // reset pair counters
	  __n2_12   = __s2ptpt_12   = __s2NPt_12    = __s2PtN_12    = 0;
	  __n2Nw_12 = __s2ptptNw_12 = __s2NPtNw_12  = __s2PtNNw_12  = 0;
	  for (int i1=0; i1<k1; i1++)
            {
	      ////cout << "         i1:" << i1 << endl;
	      id_1      = _id_1[i1];           ////cout << "       id_1:" << id_1 << endl;
	      q_1       = _charge_1[i1];       ////cout << "        q_1:" << q_1 << endl;
	      iEtaPhi_1 = _iEtaPhi_1[i1];      ////cout << "  iEtaPhi_1:" << iEtaPhi_1 << endl;
	      iPt_1     = _iPt_1[i1];          ////cout << "      iPt_1:" << iPt_1 << endl;
	      corr_1    = _correction_1[i1];   ////cout << "     corr_1:" << corr_1 << endl;
	      pt_1      = _pt_1[i1];           ////cout << "       pt_1:" << pt_1 << endl;
	      px_1      = _px_1[i1];          ////cout << "      px_1:" << px_1 << endl;
	      py_1      = _py_1[i1];          ////cout << "      py_1:" << py_1 << endl;
	      pz_1      = _pz_1[i1];          ////cout << "      pz_1:" << pz_1 << endl;
	      dedx_1    = _dedx_1[i1];        ////cout << "     dedx_1:" << dedx_1 << endl;
                
	      //1 and 2
	      for (int i2=0; i2<k2; i2++)
                {
		  ////cout << "         i2:" << i2 << endl;
		  id_2   = _id_2[i2];              ////cout << "       id_2:" << id_2 << endl;
		  if (id_1!=id_2)  // exclude auto correlation
                    {
		      q_2       = _charge_2[i2];     ////cout << "        q_2:" << q_2 << endl;
		      iEtaPhi_2 = _iEtaPhi_2[i2];    ////cout << "  iEtaPhi_2:" << iEtaPhi_2 << endl;
		      iPt_2     = _iPt_2[i2];        ////cout << "      iPt_2:" << iPt_2 << endl;
		      corr_2    = _correction_2[i2]; ////cout << "     corr_2:" << corr_2 << endl;
		      pt_2      = _pt_2[i2];         ////cout << "       pt_2:" << pt_2 << endl;
		      px_2      = _px_2[i2];          ////cout << "      px_2:" << px_2 << endl;
		      py_2      = _py_2[i2];          ////cout << "      py_2:" << py_2 << endl;
		      pz_2      = _pz_2[i2];          ////cout << "      pz_2:" << pz_2 << endl;
		      dedx_2    = _dedx_2[i2];        ////cout << "     dedx_2:" << dedx_2 << endl;                   

		      corr      = corr_1*corr_2;
		      ij        = iEtaPhi_1*_nBins_etaPhi_1 + iEtaPhi_2;   ////cout << " ij:" << ij<< endl;
		      __n2_12                  += corr;
		      __n2_12_vsEtaPhi[ij]     += corr;
		      ptpt                     = pt_1*pt_2;
		      __s2ptpt_12              += corr*ptpt;
		      __s2PtN_12               += corr*pt_1;
		      __s2NPt_12               += corr*pt_2;
		      __s2ptpt_12_vsEtaPhi[ij] += corr*ptpt;
		      __s2PtN_12_vsEtaPhi[ij]  += corr*pt_1;
		      __s2NPt_12_vsEtaPhi[ij]  += corr*pt_2;
		      __n2_12_vsPtPt[iPt_1*_nBins_pt_2 + iPt_2] += corr;
		      __n2Nw_12                  += 1;
		      __s2ptptNw_12              += ptpt;
		      __s2PtNNw_12               += pt_1;
		      __s2NPtNw_12               += pt_2;
                    }
                } //i2
            } //i1
        }
        
      _n2_12_vsM->Fill(centrality,     __n2_12);
      _s2PtPt_12_vsM->Fill(centrality, __s2ptpt_12);
      _s2PtN_12_vsM->Fill(centrality,  __s2NPt_12);
      _s2NPt_12_vsM->Fill(centrality,  __s2PtN_12);
        
      _n2Nw_12_vsM->Fill(centrality,     __n2Nw_12);
      _s2PtPtNw_12_vsM->Fill(centrality, __s2ptptNw_12);
      _s2PtNNw_12_vsM->Fill(centrality,  __s2NPtNw_12);
      _s2NPtNw_12_vsM->Fill(centrality,  __s2PtNNw_12);   
    }    
    
  AliInfo("AliAnalysisTaskPIDBFDptDpt::UserExec()   -----------------Event Done ");
  PostData(0,_outputHistoList);    
} // End of UserExec

void   AliAnalysisTaskPIDBFDptDpt::FinishTaskOutput()
{
  AliInfo("AliAnalysisTaskPIDBFDptDpt::FinishTaskOutput() Starting.");
  Printf("= 0 ====================================================================");
  finalizeHistograms();
  AliInfo("= 1 ====================================================================");
  PostData(0,_outputHistoList);
  AliInfo("= 2 ====================================================================");
  AliInfo("AliAnalysisTaskPIDBFDptDpt::FinishTaskOutput() Done.");
}

void   AliAnalysisTaskPIDBFDptDpt::Terminate(Option_t* /*option*/)
{
  AliInfo("AliAnalysisTaskPIDBFDptDpt::Terminate() Starting/Done.");
}


//Tools
//===================================================================================================
void  AliAnalysisTaskPIDBFDptDpt::fillHistoWithArray(TH1 * h, float * array, int size)
{
  int i, i1;
  float v1, ev1, v2, ev2, sum, esum;
  for (i=0, i1=1; i<size; ++i,++i1)
    {
      v1  = array[i]; ev1 = sqrt(v1);
      v2  = h->GetBinContent(i1);
      ev2 = h->GetBinError(i1);
      sum = v1 + v2;
      esum = sqrt(ev1*ev1+ev2*ev2);
      h->SetBinContent(i1,sum);
      h->SetBinError(i1,esum);
    }
}

void  AliAnalysisTaskPIDBFDptDpt::fillHistoWithArray(TH2 * h, float * array, int size1, int size2)
{
  int i, i1;
  int j, j1;
  float v1, ev1, v2, ev2, sum, esum;
  for (i=0, i1=1; i<size1; ++i,++i1)
    {
      for (j=0, j1=1; j<size2; ++j,++j1)
        {
	  v1  = array[i*size2+j]; ev1 = sqrt(v1);
	  v2  = h->GetBinContent(i1,j1);
	  ev2 = h->GetBinError(i1,j1);
	  sum = v1 + v2;
	  esum = sqrt(ev1*ev1+ev2*ev2);
	  h->SetBinContent(i1,j1,sum);
	  h->SetBinError(i1,j1,esum);
        }
    }
}

void  AliAnalysisTaskPIDBFDptDpt::fillHistoWithArray(TH3 * h, float * array, int size1, int size2, int size3)
{
  int i, i1;
  int j, j1;
  int k, k1;
  float v1, ev1, v2, ev2, sum, esum;
  int size23 = size2*size3;
  for (i=0, i1=1; i<size1; ++i,++i1)
    {
      for (j=0, j1=1; j<size2; ++j,++j1)
        {
	  for (k=0, k1=1; k<size3; ++k,++k1)
            {
	      v1  = array[i*size23+j*size3+k]; ev1 = sqrt(v1);
	      v2  = h->GetBinContent(i1,j1,k1);
	      ev2 = h->GetBinError(i1,j1,k1);
	      sum = v1 + v2;
	      esum = sqrt(ev1*ev1+ev2*ev2);
	      h->SetBinContent(i1,j1,k1,sum);
	      h->SetBinError(i1,j1,k1,esum);
            }
        }
    }
}

void  AliAnalysisTaskPIDBFDptDpt::fillHistoWithArray(TH1 * h, double * array, int size)
{
  int i, i1;
  double v1, ev1, v2, ev2, sum, esum;
  for (i=0, i1=1; i<size; ++i,++i1)
    {
      v1  = array[i]; ev1 = sqrt(v1);
      v2  = h->GetBinContent(i1);
      ev2 = h->GetBinError(i1);
      sum = v1 + v2;
      esum = sqrt(ev1*ev1+ev2*ev2);
      h->SetBinContent(i1,sum);
      h->SetBinError(i1,esum);
    }
}

void  AliAnalysisTaskPIDBFDptDpt::fillHistoWithArray(TH2 * h, double * array, int size1, int size2)
{
  int i, i1;
  int j, j1;
  double v1, ev1, v2, ev2, sum, esum;
  for (i=0, i1=1; i<size1; ++i,++i1)
    {
      for (j=0, j1=1; j<size2; ++j,++j1)
        {
	  v1  = array[i*size2+j]; ev1 = sqrt(v1);
	  v2  = h->GetBinContent(i1,j1);
	  ev2 = h->GetBinError(i1,j1);
	  sum = v1 + v2;
	  esum = sqrt(ev1*ev1+ev2*ev2);
	  h->SetBinContent(i1,j1,sum);
	  h->SetBinError(i1,j1,esum);
        }
    }
}

void  AliAnalysisTaskPIDBFDptDpt::fillHistoWithArray(TH3 * h, double * array, int size1, int size2, int size3)
{
  int i, i1;
  int j, j1;
  int k, k1;
  double v1, ev1, v2, ev2, sum, esum;
  int size23 = size2*size3;
  for (i=0, i1=1; i<size1; ++i,++i1)
    {
      for (j=0, j1=1; j<size2; ++j,++j1)
        {
	  for (k=0, k1=1; k<size3; ++k,++k1)
            {
	      v1  = array[i*size23+j*size3+k]; ev1 = sqrt(v1);
	      v2  = h->GetBinContent(i1,j1,k1);
	      ev2 = h->GetBinError(i1,j1,k1);
	      sum = v1 + v2;
	      esum = sqrt(ev1*ev1+ev2*ev2);
	      h->SetBinContent(i1,j1,k1,sum);
	      h->SetBinError(i1,j1,k1,esum);
            }
        }
    }
}

//________________________________________________________________________
double *  AliAnalysisTaskPIDBFDptDpt::getDoubleArray(int size, double v)
{
  /// Allocate an array of type double with n values
  /// Initialize the array to the given value
  double * array = new double [size];
  for (int i=0;i<size;++i) array[i]=v;
  return array;
}

//________________________________________________________________________
float *  AliAnalysisTaskPIDBFDptDpt::getFloatArray(int size, float v)
{
  /// Allocate an array of type float with n values
  /// Initialize the array to the given value
  float * array = new float [size];
  for (int i=0;i<size;++i) array[i]=v;
  return array;
}


//________________________________________________________________________
TH1D * AliAnalysisTaskPIDBFDptDpt::createHisto1D(const TString &  name, const TString &  title, 
						 int n, double xMin, double xMax, 
						 const TString &  xTitle, const TString &  yTitle)
{
  //CreateHisto new 1D historgram
  AliInfo(Form("createHisto 1D histo %s   nBins: %d  xMin: %f   xMax: %f",name.Data(),n,xMin,xMax));
  TH1D * h = new TH1D(name,title,n,xMin,xMax);
  h->GetXaxis()->SetTitle(xTitle);
  h->GetYaxis()->SetTitle(yTitle);
  addToList(h);
  return h;
}


//________________________________________________________________________
TH1D * AliAnalysisTaskPIDBFDptDpt::createHisto1D(const TString &  name, const TString &  title, 
						 int n, double * bins, 
						 const TString &  xTitle, const TString &  yTitle)
{
  AliInfo(Form("createHisto 1D histo %s   with %d non uniform nBins",name.Data(),n));
  TH1D * h = new TH1D(name,title,n,bins);
  h->GetXaxis()->SetTitle(xTitle);
  h->GetYaxis()->SetTitle(yTitle);
  addToList(h);
  return h;
}


//________________________________________________________________________
TH2D * AliAnalysisTaskPIDBFDptDpt::createHisto2D(const TString &  name, const TString &  title, 
						 int nx, double xMin, double xMax, int ny, double yMin, double yMax, 
						 const TString &  xTitle, const TString &  yTitle, const TString &  zTitle)
{
  AliInfo(Form("createHisto 2D histo %s  nx: %d  xMin: %f10.4 xMax: %f10.4  ny: %d   yMin: %f10.4 yMax: %f10.4",name.Data(),nx,xMin,xMax,ny,yMin,yMax));
  TH2D * h = new TH2D(name,title,nx,xMin,xMax,ny,yMin,yMax);
  h->GetXaxis()->SetTitle(xTitle);
  h->GetYaxis()->SetTitle(yTitle);
  h->GetZaxis()->SetTitle(zTitle);
  addToList(h);
  return h;
}

//________________________________________________________________________
TH2D * AliAnalysisTaskPIDBFDptDpt::createHisto2D(const TString &  name, const TString &  title, 
						 int nx, double* xbins, int ny, double yMin, double yMax, 
						 const TString &  xTitle, const TString &  yTitle, const TString &  zTitle)
{
  AliInfo(Form("createHisto 2D histo %s   with %d non uniform nBins",name.Data(),nx));
  TH2D * h;
  h = new TH2D(name,title,nx,xbins,ny,yMin,yMax);
  h->GetXaxis()->SetTitle(xTitle);
  h->GetYaxis()->SetTitle(yTitle);
  h->GetZaxis()->SetTitle(zTitle);
  addToList(h);
  return h;
}

//// F /////
//________________________________________________________________________
TH1F * AliAnalysisTaskPIDBFDptDpt::createHisto1F(const TString &  name, const TString &  title, 
						 int n, double xMin, double xMax, 
						 const TString &  xTitle, const TString &  yTitle)
{
  //CreateHisto new 1D historgram
  AliInfo(Form("createHisto 1D histo %s   nBins: %d  xMin: %f   xMax: %f",name.Data(),n,xMin,xMax));
  TH1F * h = new TH1F(name,title,n,xMin,xMax);
  h->GetXaxis()->SetTitle(xTitle);
  h->GetYaxis()->SetTitle(yTitle);
  addToList(h);
  return h;
}


//________________________________________________________________________
TH1F * AliAnalysisTaskPIDBFDptDpt::createHisto1F(const TString &  name, const TString &  title, 
						 int n, double * bins, 
						 const TString &  xTitle, const TString &  yTitle)
{
  AliInfo(Form("createHisto 1D histo %s   with %d non uniform nBins",name.Data(),n));
  TH1F * h = new TH1F(name,title,n,bins);
  h->GetXaxis()->SetTitle(xTitle);
  h->GetYaxis()->SetTitle(yTitle);
  addToList(h);
  return h;
}


//________________________________________________________________________
TH2F * AliAnalysisTaskPIDBFDptDpt::createHisto2F(const TString &  name, const TString &  title, 
						 int nx, double xMin, double xMax, int ny, double yMin, double yMax, 
						 const TString &  xTitle, const TString &  yTitle, const TString &  zTitle)
{
  AliInfo(Form("createHisto 2D histo %s  nx: %d  xMin: %f10.4 xMax: %f10.4  ny: %d   yMin: %f10.4 yMax: %f10.4",name.Data(),nx,xMin,xMax,ny,yMin,yMax));
  TH2F * h = new TH2F(name,title,nx,xMin,xMax,ny,yMin,yMax);
  h->GetXaxis()->SetTitle(xTitle);
  h->GetYaxis()->SetTitle(yTitle);
  h->GetZaxis()->SetTitle(zTitle);
  addToList(h);
  return h;
}

//________________________________________________________________________
TH2F * AliAnalysisTaskPIDBFDptDpt::createHisto2F(const TString &  name, const TString &  title, 
						 int nx, double* xbins, int ny, double yMin, double yMax, 
						 const TString &  xTitle, const TString &  yTitle, const TString &  zTitle)
{
  AliInfo(Form("createHisto 2D histo %s   with %d non uniform nBins",name.Data(),nx));
  TH2F * h;
  h = new TH2F(name,title,nx,xbins,ny,yMin,yMax);
  h->GetXaxis()->SetTitle(xTitle);
  h->GetYaxis()->SetTitle(yTitle);
  h->GetZaxis()->SetTitle(zTitle);
  addToList(h);
  return h;
}


//________________________________________________________________________
TH3F * AliAnalysisTaskPIDBFDptDpt::createHisto3F(const TString &  name, const TString &  title, 
						 int nx, double xMin, double xMax, 
						 int ny, double yMin, double yMax, 
						 int nz, double zMin, double zMax, 
						 const TString &  xTitle, const TString &  yTitle, const TString &  zTitle)
{
  AliInfo(Form("createHisto 2D histo %s  nx: %d  xMin: %f10.4 xMax: %f10.4  ny: %d   yMin: %f10.4 yMax: %f10.4 nz: %d   zMin: %f10.4 zMax: %f10.4",name.Data(),nx,xMin,xMax,ny,yMin,yMax,nz,zMin,zMax));
  TH3F * h = new TH3F(name,title,nx,xMin,xMax,ny,yMin,yMax,nz,zMin,zMax);
  h->GetXaxis()->SetTitle(xTitle);
  h->GetYaxis()->SetTitle(yTitle);
  h->GetZaxis()->SetTitle(zTitle);
  addToList(h);
  return h;
}


//________________________________________________________________________
TProfile * AliAnalysisTaskPIDBFDptDpt::createProfile(const TString & name, const TString & description, 
						     int nx,double xMin,double xMax,
						     const TString &  xTitle, const TString &  yTitle)
{
  AliInfo(Form("createHisto 1D profile %s   nBins: %d  xMin: %f10.4 xMax: %f10.4",name.Data(),nx,xMin,xMax));
  TProfile * h = new TProfile(name,description,nx,xMin,xMax);
  h->GetXaxis()->SetTitle(xTitle);
  h->GetYaxis()->SetTitle(yTitle);
  addToList(h);
  return h;  
}

//________________________________________________________________________
TProfile * AliAnalysisTaskPIDBFDptDpt::createProfile(const TString &  name,const TString &  description, 
						     int nx,  double* bins,
						     const TString &  xTitle, const TString &  yTitle)
{
  AliInfo(Form("createHisto 1D profile %s  with %d non uniform bins",name.Data(),nx));
  TProfile * h = new TProfile(name,description,nx,bins);
  h->GetXaxis()->SetTitle(xTitle);
  h->GetYaxis()->SetTitle(yTitle);
  addToList(h);
  return h;
}

//________________________________________________________________________
void   AliAnalysisTaskPIDBFDptDpt::addToList(TH1 *h)
{
  if (_outputHistoList)
    {
      _outputHistoList->Add(h);
    }
  else
    AliInfo("addToList(TH1 *h) _outputHistoList is null!!!!! Should abort ship");
    
}

//____________________________________________________________________________________________________
void AliAnalysisTaskPIDBFDptDpt::CalculateNSigmas( AliVTrack * trk )   //defines data member fnsigmas
{  
  // Compute nsigma for each hypthesis
  AliVParticle * inEvHMain = dynamic_cast < AliVParticle * > ( trk );
  // --- TPC
  Double_t nsigmaTPCkProton   = fPIDResponse -> NumberOfSigmasTPC( inEvHMain, AliPID::kProton );
  Double_t nsigmaTPCkKaon     = fPIDResponse -> NumberOfSigmasTPC( inEvHMain, AliPID::kKaon ); 
  Double_t nsigmaTPCkPion     = fPIDResponse -> NumberOfSigmasTPC( inEvHMain, AliPID::kPion );
  Double_t nsigmaTPCkElectron = fPIDResponse -> NumberOfSigmasTPC( inEvHMain, AliPID::kElectron );
  //set data member fnsigmas
  fnsigmas[0][0] = nsigmaTPCkPion;
  fnsigmas[1][0] = nsigmaTPCkKaon;
  fnsigmas[2][0] = nsigmaTPCkProton;
  fnsigmas[3][0] = nsigmaTPCkElectron;  
  
  // --- TOF
  Double_t nsigmaTOFkProton = 999., nsigmaTOFkKaon = 999., nsigmaTOFkPion = 999., nsigmaTOFkElectron = 999.;

  CheckTOF( trk ); // check TOF matching

  if( fHasTOFPID )
    {  // use TOF information when there is TOF matching for a track
      nsigmaTOFkProton   = fPIDResponse -> NumberOfSigmasTOF( inEvHMain, AliPID::kProton   );
      nsigmaTOFkKaon     = fPIDResponse -> NumberOfSigmasTOF( inEvHMain, AliPID::kKaon     ); 
      nsigmaTOFkPion     = fPIDResponse -> NumberOfSigmasTOF( inEvHMain, AliPID::kPion     );
      nsigmaTOFkElectron = fPIDResponse -> NumberOfSigmasTOF( inEvHMain, AliPID::kElectron );

      //set data member fnsigmas
      fnsigmas[0][1] = nsigmaTOFkPion;
      fnsigmas[1][1] = nsigmaTOFkKaon;
      fnsigmas[2][1] = nsigmaTOFkProton;
      fnsigmas[3][1] = nsigmaTOFkElectron;
    }
}

//____________________________________________________________________________________________________
void AliAnalysisTaskPIDBFDptDpt::CalculateTPCNSigmasElectron( AliVTrack * trk )   //defines data member fnsigmas
{  
  // Compute nsigma for each hypthesis
  AliVParticle * inEvHMain = dynamic_cast < AliVParticle * > ( trk );
  // --- TPC
  Double_t nsigmaTPCkElectron = fPIDResponse -> NumberOfSigmasTPC( inEvHMain, AliPID::kElectron );
  //set data member fnsigmas
  fnsigmas[3][0] = nsigmaTPCkElectron;
}

//________________________________________________________________________________________________________
void AliAnalysisTaskPIDBFDptDpt::CheckTOF( AliVTrack * trk )     //check if the particle has TOF Matching
{ 
  //get the PIDResponse
  if( fPIDResponse -> CheckPIDStatus( AliPIDResponse::kTOF, trk ) == 0 )   fHasTOFPID = kFALSE;
  else fHasTOFPID = kTRUE;
  
  if( fRemoveTracksT0Fill ) //what does this block of code do???
    {
      Int_t startTimeMask = fPIDResponse -> GetTOFResponse().GetStartTimeMask( trk->P() );
      if ( startTimeMask < 0 )   fHasTOFPID = kFALSE; 
    }
}

//________________________________________________________________________________________________________________
Int_t AliAnalysisTaskPIDBFDptDpt::TellParticleSpecies( AliVTrack * trk )  //function to determine the particle ID
{  
  CalculateNSigmas( trk );

  Double_t nsigmaPion = 999., nsigmaKaon = 999., nsigmaProton = 999., nsigmaElectron = 999.;

  CheckTOF( trk );
  
  if( fHasTOFPID && ( trk->Pt() > ptTOFlowerBoundary ) )
    {
      nsigmaPion     =  TMath::Abs( fnsigmas[0][1] ); //Pion_TOF
      nsigmaKaon     =  TMath::Abs( fnsigmas[1][1] ); //Kaon_TOF
      nsigmaProton   =  TMath::Abs( fnsigmas[2][1] ); //Proton_TOF
      nsigmaElectron =  TMath::Abs( fnsigmas[3][1] ); //Electron_TOF

      if( ( nsigmaPion < fNSigmaPID ) && ( nsigmaKaon > fNSigmaPID ) && ( nsigmaProton > fNSigmaPID ) && ( trk->Pt() <= ptUpperLimit ) )  return 0; //Pion
      if( ( nsigmaKaon < fNSigmaPID ) && ( nsigmaPion > fNSigmaPID ) && ( nsigmaProton > fNSigmaPID ) && ( trk->Pt() <= ptUpperLimit ) )  return 1; //Kaon
      if( ( nsigmaProton < fNSigmaPID ) && ( nsigmaPion > fNSigmaPID ) && ( nsigmaKaon > fNSigmaPID ) && ( trk->Pt() <= ptUpperLimit ) && ( massSquareCalculation(trk) > 0.6 ) && ( massSquareCalculation(trk) < 1.1 ) )   return 2; //Proton // need to add mass square cut in the source code as well!!!
    }
  else
    {
      nsigmaPion     =  TMath::Abs( fnsigmas[0][0] ); //Pion_TPC
      nsigmaKaon     =  TMath::Abs( fnsigmas[1][0] ); //Kaon_TPC
      nsigmaProton   =  TMath::Abs( fnsigmas[2][0] ); //Proton_TPC
      nsigmaElectron =  TMath::Abs( fnsigmas[3][0] ); //Electron_TPC

      if( ( nsigmaPion < fNSigmaPID ) && ( nsigmaKaon > fNSigmaPID ) && ( nsigmaProton > fNSigmaPID ) && ( nsigmaElectron > electronNSigmaVeto ) && ( trk->Pt() >= _min_pt_1 ) && ( trk->Pt() <= ptTOFlowerBoundary ) )   return 0; //Pion
      if( ( nsigmaKaon < fNSigmaPID ) && ( nsigmaPion > fNSigmaPID ) && ( nsigmaProton > fNSigmaPID ) && ( nsigmaElectron > electronNSigmaVeto ) && ( trk->Pt() >= _min_pt_1 ) && ( trk->Pt() <= ptTOFlowerBoundary ) )   return 1; //Kaon
      if( ( nsigmaProton < fNSigmaPID ) && ( nsigmaPion > fNSigmaPID ) && ( nsigmaKaon > fNSigmaPID ) && ( nsigmaElectron > electronNSigmaVeto ) && ( trk->Pt() >= _min_pt_1 ) && ( trk->Pt() <= ptTOFlowerBoundary ) )   return 2; //Proton
    }

  // else, return undefined
  return 999;
}

//________________________________________________________________________________________________________________
Double_t AliAnalysisTaskPIDBFDptDpt::TOFBetaCalculation( AliVTrack * track ) const
{
  Double_t tofTime = 0;
  if ( fAnalysisType == "MCAODreco" )       tofTime = track -> GetTOFsignalTunedOnData();
  else if ( fAnalysisType == "RealData" )   tofTime = track -> GetTOFsignal();  

  Double_t c = TMath::C() * 1.E-9;// m/ns
  Float_t startTime = fPIDResponse -> GetTOFResponse().GetStartTime( track->P() );//in ps
  Double_t length = fPIDResponse -> GetTOFResponse().GetExpectedSignal( track, AliPID::kElectron ) * 1E-3 * c;
  tofTime -= startTime;      // subtract startTime to the signal
  Double_t tof = tofTime * 1E-3; // ns, average T0 fill subtracted, no info from T0detector 	 
  tof = tof * c;
  return length / tof;
}

//________________________________________________________________________________________________________________
Double_t AliAnalysisTaskPIDBFDptDpt::massSquareCalculation( AliVTrack * track ) const
{
  Double_t massSquare = track->P() * track->P() * ( 1 / ( TOFBetaCalculation(track) * TOFBetaCalculation(track) ) - 1 ); 
  return massSquare;
}
