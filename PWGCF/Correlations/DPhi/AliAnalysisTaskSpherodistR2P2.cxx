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
#include <fstream>
#include <iostream>
#include "AliAnalysisManager.h"
#include "AliAODHandler.h"
#include "AliAODInputHandler.h"
#include "AliInputEventHandler.h"
#include "AliLog.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliMultiplicity.h"
#include "AliCentrality.h"
#include "AliAnalysisTaskSpherodistR2P2.h"
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
#include "TVector2.h"
#include "AliEventCuts.h"
#include <AliVVZERO.h>
#include "AliHelperPID.h"
#include "TBits.h"
#include "AliPPVsMultUtils.h"
#include "AliVEvent.h"


using namespace AliHelperPIDNameSpace;
using namespace std;

ClassImp(AliAnalysisTaskSpherodistR2P2)

AliAnalysisTaskSpherodistR2P2::AliAnalysisTaskSpherodistR2P2()
: AliAnalysisTaskSE(),
  fNSigmaPID( 2 ),
  fNSigmaPID_veto( 2 ),
  ptUpperLimit( 2.0 ),
  ptTOFlowerBoundary( 0.5),
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
  _pileUpEventPP13 ( 0),
  _pileUpTrackPP13 ( 0),
  PIDparticle   ( 0),
  use_pT_cut   ( 0),
  useAliHelperPID( 0),
  useCircularCutPID( 0),
  fHelperPID(0),
  NoContamination   ( 0),
  NoContaminationWeak   ( 0),
  NoContaminationWeakMaterial   ( 0),
  Closure_NoMisIDWeakMaterial   ( 0),
  _usePtEff    ( 0),
  _useRapidity   ( 0),
  _useEventPlane   ( 0),
  EP_min( -3.1415927/6 ),
  EP_max( 3.1415927/6 ),
  _psi_EventPlane ( 0),
  _sameFilter    ( false),
  _rejectPileup  ( 1),
  _rejectPairConversion ( 0),
  _vertexZMin           ( -6),
  _vertexZMax           (  6),
  _NchMin               (5  ),
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
  _mult4a   ( 0 ),
  _mult5    ( 0 ),
  _mult6    ( 0 ),
  _mult7    ( 0 ),
  _mult8    ( 0 ),
  arraySize ( 5000),
  _id_1(0),
  _charge_1(0),
  _iEtaPhi_1(0),
  _iPt_1(0),
  _TrackArray(0),
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
  _correctionPtEff_1(0),
  _correctionPtEff_2(0),
  _nBins_M0(500),       _min_M0(0),        _max_M0(10000),          _width_M0(20),
  _nBins_M1(500),       _min_M1(0),        _max_M1(10000),          _width_M1(20),
  _nBins_M2(500),       _min_M2(0),        _max_M2(10000),          _width_M2(20),
  _nBins_M3(500),       _min_M3(0),        _max_M3(10000),          _width_M3(20),
  _nBins_M4(100),       _min_M4(0),        _max_M4(1),              _width_M4(0.01),
  _nBins_M5(100),       _min_M5(0),        _max_M5(1),              _width_M5(0.01),
  _nBins_M6(100),       _min_M6(0),        _max_M6(1),              _width_M6(0.01),
  _nBins_M7(100),       _min_M7(0),        _max_M7(1),              _width_M7(0.01),
  _nBins_M8(100),       _min_M8(0),        _max_M8(1),              _width_M8(0.01),
  _nBins_vertexZ(0),   _min_vertexZ(0), _max_vertexZ(0),        _width_vertexZ(0.5),
  
  _nBins_pt_1(18),      _min_pt_1(0.2),    _max_pt_1(2.0),          _width_pt_1(0.1),
  _nBins_phi_1(72),     _min_phi_1(0),     _max_phi_1(2.*3.1415927),_width_phi_1(2.*3.1415927/72.),
  _nBins_eta_1(0),      _min_eta_1(0),  _max_eta_1(0),           _width_eta_1(0.1),
  
  _nBins_etaPhi_1(0),
  _nBins_etaPhiPt_1(0),
  _nBins_zEtaPhiPt_1(0),
  
  _nBins_pt_2(18),     _min_pt_2(0.2),     _max_pt_2(2.0),          _width_pt_2(0.1),
  _nBins_phi_2(72),    _min_phi_2(0),      _max_phi_2(2.*3.1415927),_width_phi_2(2.*3.1415927/72.),
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
  __n1_1_vsPt_pdg(0),
  __n1_1_vsPt_pdg_Weak(0),
  __n1_1_vsPt_pdg_Weak_Material(0),
  __n1_1_vsPt_Weak(0),
  __n1_1_vsPt_Material(0),
  __n1_1_vsEtaPhi(0),
  __s1pt_1_vsEtaPhi(0),
  // __n1_1_vsZEtaPhiPt(0),
  __n1_2_vsPt(0),
  __n1_2_vsPt_pdg(0),
  __n1_2_vsPt_pdg_Weak(0),
  __n1_2_vsPt_pdg_Weak_Material(0),
  __n1_2_vsPt_Weak(0),
  __n1_2_vsPt_Material(0),
  __n1_2_vsEtaPhi(0),
  __s1pt_2_vsEtaPhi(0),
  // __n1_2_vsZEtaPhiPt(0),
  __n2_12_vsPtPt(0),
  __n2_12_vsEtaPhi(0),
  __s2ptpt_12_vsEtaPhi(0),
  __s2PtN_12_vsEtaPhi(0),
  __s2NPt_12_vsEtaPhi(0),
  _hPtEff_1      ( 0    ),
  _hPtEff_2      ( 0    ),
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
  _m9 ( 0),
  _m10 ( 0),
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
  
  _dedx_p (0),
  _beta_p (0), 
  _inverse_beta_p (0),
  _msquare_p (0),

  _dedx_p_POI_AliHelperPID (0),
  _beta_p_POI_AliHelperPID (0),
  _inverse_beta_p_POI_AliHelperPID (0),
  _msquare_p_POI_AliHelperPID (0),
  _nSigmaTOF_p_pion_before (0),
  _nSigmaTOF_p_pion_after (0),
  _nSigmaTOF_p_kaon_before (0),
  _nSigmaTOF_p_kaon_after (0),
  _nSigmaTOF_p_proton_before (0),
  _nSigmaTOF_p_proton_after (0),
  _nSigmaTPC_nSigmaTOF_Pion_before (0),
  _nSigmaTPC_nSigmaTOF_Pion_after (0),
  _nSigmaTPC_nSigmaTOF_Kaon_before (0),
  _nSigmaTPC_nSigmaTOF_Kaon_after (0),
  _nSigmaTPC_nSigmaTOF_Proton_before (0),
  _nSigmaTPC_nSigmaTOF_Proton_after (0),
  
  _dedx_p_AliHelperPID_no_Undefined (0),
  _beta_p_AliHelperPID_no_Undefined (0),
  _inverse_beta_p_AliHelperPID_no_Undefined (0),
  _msquare_p_AliHelperPID_no_Undefined (0),
  _fhV0MvsTracksTPCout_after(0),

  _etadis_POI_AliHelperPID ( 0),
  _ydis_POI_AliHelperPID ( 0),
  _etadis_before_any_cuts ( 0),
  
  _y_Pt_AllCh_MCAODTruth ( 0 ),
  _y_Pt_POI_MCAODTruth ( 0 ),
  
  _phidis_POI_AliHelperPID ( 0),
  _phidis_before_any_cuts ( 0),
  
  _dcaz   ( 0),
  _dcaxy  ( 0),
  _n1_1_vsPt         ( 0),
  _n1_1_vsPt_pdg(0),
  _n1_1_vsPt_pdg_Weak(0),
  _n1_1_vsPt_pdg_Weak_Material(0),
  _n1_1_vsPt_Weak(0),
  _n1_1_vsPt_Material(0),
  _n1_1_vsEtaVsPhi   ( 0),
  _s1pt_1_vsEtaVsPhi ( 0),
  _n1_1_vsM          ( 0),
  _s1pt_1_vsM        ( 0),
  _n1Nw_1_vsM        ( 0),
  _s1ptNw_1_vsM      ( 0),
  _dedxVsP_1         ( 0),
  _corrDedxVsP_1     ( 0),
  _betaVsP_1         ( 0),
  _n1_2_vsPt         ( 0),
  _n1_2_vsPt_pdg(0),
  _n1_2_vsPt_pdg_Weak(0),
  _n1_2_vsPt_pdg_Weak_Material(0),
  _n1_2_vsPt_Weak(0),
  _n1_2_vsPt_Material(0),
  _n1_2_vsEtaVsPhi   ( 0),
  _s1pt_2_vsEtaVsPhi ( 0),
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
  _invMassKaon       ( 0),
  _invMassKaonSq     ( 0),
  _invMassElec       ( 0),
  fSharedfraction_Pair_cut( 0 ),
  _ClusterSharedFraction_beforeCut (0),
  _ClusterSharedFraction_afterCut (0),
  _ClusterSharedFraction_3by3Bins_beforeCut (0),
  _ClusterSharedFraction_3by3Bins_afterCut (0),
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
  _title_m9("NA"),
  _title_m10("NA"),

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
  pdg("NA"),
  Weak("NA"),
  Material("NA"),
  fUtils(0),
  fPPVsMultUtils(0),
  f2015V0MtoTrkTPCout(NULL),
  f2015V0MtoTrkTPCout_Upper(NULL),
  fV0Multiplicity(0),
  fV0Multiplicity_Victor(0),
  fNoOfTPCoutTracks(0),
  fEventCut(0),
  hSpherocity(0),
  h1i_n1_multSphero(0)
    {
      // Au-Au added this block of code to use his own PID functions
      for( Int_t ipart = 0; ipart < 4; ipart++ )
      for( Int_t ipid = 0; ipid < 2; ipid++ )
      fnsigmas[ipart][ipid] = 999.;

      printf("Default constructor called \n");  
      printf("passed \n ");
    }

AliAnalysisTaskSpherodistR2P2::AliAnalysisTaskSpherodistR2P2(const TString & name)
: AliAnalysisTaskSE(name),
fNSigmaPID( 2 ),
fNSigmaPID_veto( 2 ),
ptUpperLimit( 2.0 ),
ptTOFlowerBoundary( 0.5),
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
_pileUpEventPP13 ( 0),
_pileUpTrackPP13 ( 0),
PIDparticle    ( 0),
use_pT_cut     ( 0),
useAliHelperPID( 0),
useCircularCutPID( 0),
fHelperPID(0),
NoContamination   ( 0),
NoContaminationWeak   ( 0),
NoContaminationWeakMaterial   ( 0),
Closure_NoMisIDWeakMaterial   ( 0),
_usePtEff    ( 0),
_useRapidity   ( 0),
_useEventPlane   ( 0),
EP_min( -3.1415927/6 ),
EP_max( 3.1415927/6 ),
_psi_EventPlane ( 0),
_sameFilter    ( false),
_rejectPileup  ( 1),
_rejectPairConversion ( 0),
_vertexZMin           ( -8.),
_vertexZMax           (  8.),
_NchMin               (5  ),
_vertexZWidth         (0.5 ),
_etaWidth             (0.1),
_vertexXYMin          ( -10.),
_vertexXYMax          (  10.),
_centralityMethod     (  4),
_centralityMin        (  0.),
_centralityMax        (  1.),
_requestedCharge_1    (   1),
_requestedCharge_2    (  -1),
_dcaZMin              ( -0.2),
_dcaZMax              (  0.2),
_dcaXYMin             ( -0.2),
_dcaXYMax             (  0.2),
_dedxMin              ( 0),
_dedxMax              ( 100000),
_nClusterMin          ( 70),
_trackFilterBit       ( 0),
fAnalysisType         ( "RealData" ),
fSystemType           ( "pp18_V0_kMB_kFALSE" ),
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
_mult4a   ( 0 ),
_mult5    ( 0 ),
_mult6    ( 0 ),
_mult7    ( 0 ),
_mult8    ( 0 ),
arraySize ( 5000),
_id_1(0),
_charge_1(0),
_iEtaPhi_1(0),
_iPt_1(0),
_TrackArray(0),
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
_correctionPtEff_1(0),
_correctionPtEff_2(0),
_nBins_M0(10000),       _min_M0(0),        _max_M0(10000),          _width_M0(1),
_nBins_M1(10000),       _min_M1(0),        _max_M1(10000),          _width_M1(1),
_nBins_M2(10000),       _min_M2(0),        _max_M2(10000),          _width_M2(1),
_nBins_M3(10000),       _min_M3(0),        _max_M3(10000),          _width_M3(1),
_nBins_M4(100),       _min_M4(0),        _max_M4(100),              _width_M4(1),
_nBins_M5(100),       _min_M5(0),        _max_M5(1),              _width_M5(0.01),
_nBins_M6(100),       _min_M6(0),        _max_M6(1),              _width_M6(0.01),
_nBins_M7(100),       _min_M7(0),        _max_M7(1),              _width_M7(0.01),
_nBins_M8(100),       _min_M8(0),        _max_M8(1),              _width_M8(0.01),
_nBins_vertexZ(0),   _min_vertexZ(0), _max_vertexZ(0),        _width_vertexZ(0.5),

_nBins_pt_1(18),      _min_pt_1(0.2),    _max_pt_1(2.0),          _width_pt_1(0.1),
_nBins_phi_1(72),     _min_phi_1(0),     _max_phi_1(2.*3.1415927),_width_phi_1(2.*3.1415927/72.),
_nBins_eta_1(0),      _min_eta_1(0),    _max_eta_1(0),           _width_eta_1(0.1),

_nBins_etaPhi_1(0),
_nBins_etaPhiPt_1(0),
_nBins_zEtaPhiPt_1(0),

_nBins_pt_2(18),     _min_pt_2(0.2),     _max_pt_2(2.0),          _width_pt_2(0.1),
_nBins_phi_2(72),    _min_phi_2(0),      _max_phi_2(2.*3.1415927),_width_phi_2(2.*3.1415927/72.),
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
__n1_1_vsPt_pdg(0),
__n1_1_vsPt_pdg_Weak(0),
__n1_1_vsPt_pdg_Weak_Material(0),
__n1_1_vsPt_Weak(0),
__n1_1_vsPt_Material(0),
__n1_1_vsEtaPhi(0),
__s1pt_1_vsEtaPhi(0),
// __n1_1_vsZEtaPhiPt(0),
__n1_2_vsPt(0),
__n1_2_vsPt_pdg(0),
__n1_2_vsPt_pdg_Weak(0),
__n1_2_vsPt_pdg_Weak_Material(0),
__n1_2_vsPt_Weak(0),
__n1_2_vsPt_Material(0),
__n1_2_vsEtaPhi(0),
__s1pt_2_vsEtaPhi(0),
//  __n1_2_vsZEtaPhiPt(0),
__n2_12_vsPtPt(0),
__n2_12_vsEtaPhi(0),
__s2ptpt_12_vsEtaPhi(0),
__s2PtN_12_vsEtaPhi(0),
__s2NPt_12_vsEtaPhi(0),
_hPtEff_1      ( 0    ),
_hPtEff_2      ( 0    ),
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
_m9 ( 0),
_m10 ( 0),

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
  
_dedx_p (0),
_beta_p (0), 
_inverse_beta_p (0),
_msquare_p (0),

_dedx_p_POI_AliHelperPID (0),
_beta_p_POI_AliHelperPID (0),
_inverse_beta_p_POI_AliHelperPID (0),
_msquare_p_POI_AliHelperPID (0),
_nSigmaTOF_p_pion_before (0),
_nSigmaTOF_p_pion_after (0),
_nSigmaTOF_p_kaon_before (0),
_nSigmaTOF_p_kaon_after (0),
_nSigmaTOF_p_proton_before (0),
_nSigmaTOF_p_proton_after (0),
_nSigmaTPC_nSigmaTOF_Pion_before (0),
_nSigmaTPC_nSigmaTOF_Pion_after (0),
_nSigmaTPC_nSigmaTOF_Kaon_before (0),
_nSigmaTPC_nSigmaTOF_Kaon_after (0),
_nSigmaTPC_nSigmaTOF_Proton_before (0),
_nSigmaTPC_nSigmaTOF_Proton_after (0),

_dedx_p_AliHelperPID_no_Undefined (0),
_beta_p_AliHelperPID_no_Undefined (0),
_inverse_beta_p_AliHelperPID_no_Undefined (0),
_msquare_p_AliHelperPID_no_Undefined (0),
_fhV0MvsTracksTPCout_after(0),

_etadis_POI_AliHelperPID ( 0),
_ydis_POI_AliHelperPID ( 0),
_etadis_before_any_cuts ( 0),

_y_Pt_AllCh_MCAODTruth ( 0 ),
_y_Pt_POI_MCAODTruth ( 0 ),
  
_phidis_POI_AliHelperPID ( 0),
_phidis_before_any_cuts ( 0),

_dcaz ( 0),
_dcaxy ( 0),
_n1_1_vsPt         ( 0),
_n1_1_vsPt_pdg(0),
_n1_1_vsPt_pdg_Weak(0),
_n1_1_vsPt_pdg_Weak_Material(0),
_n1_1_vsPt_Weak(0),
_n1_1_vsPt_Material(0),
_n1_1_vsEtaVsPhi   ( 0),
_s1pt_1_vsEtaVsPhi ( 0),
_n1_1_vsM          ( 0),
_s1pt_1_vsM        ( 0),
_n1Nw_1_vsM        ( 0),
_s1ptNw_1_vsM      ( 0),
_dedxVsP_1         ( 0),
_corrDedxVsP_1     ( 0),
_betaVsP_1         ( 0),
_n1_2_vsPt         ( 0),
_n1_2_vsPt_pdg(0),
_n1_2_vsPt_pdg_Weak(0),
_n1_2_vsPt_pdg_Weak_Material(0),
_n1_2_vsPt_Weak(0),
_n1_2_vsPt_Material(0),
_n1_2_vsEtaVsPhi   ( 0),
_s1pt_2_vsEtaVsPhi ( 0),
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
_invMassKaon       ( 0),
_invMassKaonSq     ( 0),
_invMassElec       ( 0),
fSharedfraction_Pair_cut( 0 ),
_ClusterSharedFraction_beforeCut (0),
_ClusterSharedFraction_afterCut (0),
_ClusterSharedFraction_3by3Bins_beforeCut (0),
_ClusterSharedFraction_3by3Bins_afterCut (0),
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
_title_m9("NA"),
_title_m10("NA"),

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
pdg("NA"),
Weak("NA"),
Material("NA"),
fUtils(0),
fPPVsMultUtils(0),
f2015V0MtoTrkTPCout(NULL),
f2015V0MtoTrkTPCout_Upper(NULL),
fV0Multiplicity(0),
fV0Multiplicity_Victor(0),
fNoOfTPCoutTracks(0),
fEventCut(0),
hSpherocity(0),
h1i_n1_multSphero(0)
{
    // Au-Au added this block of code to use his own PID functions
    for( Int_t ipart = 0; ipart < 4; ipart++ )
    for( Int_t ipid = 0; ipid < 2; ipid++ )
    fnsigmas[ipart][ipid] = 999.;
  
    printf("2nd constructor called ");
    //DefineOutput(0, TList::Class()); 
    DefineOutput(1, TList::Class());   
    printf("passed  "); 
}

AliAnalysisTaskSpherodistR2P2::~AliAnalysisTaskSpherodistR2P2()
{

}

void AliAnalysisTaskSpherodistR2P2::UserCreateOutputObjects()
{
    _outputHistoList = new TList();
    _outputHistoList->SetOwner();

    Int_t nBin = 1000; 
    Double_t minVal = -1., maxVal = 999.;
    hSpherocity  = new TH1D("hSpherocity", "Spherocity", 100, 0, 1.0);
    h1i_n1_multSphero    = new TH1D("h1i_n1_multSphero", "N_{ch}", nBin, minVal, maxVal );//500,0,500);

    _outputHistoList->Add(hSpherocity);
    _outputHistoList->Add(h1i_n1_multSphero);
    


    
   
   if ( fSystemType == "PbPb_2015_kTRUE" || fSystemType == "PbPb_2015_kFALSE" )
    {
        fEventCut = new AliEventCuts();
        fEventCut->AddQAplotsToList(_outputHistoList, kTRUE);
        //fEventCut->SetManualMode();
        fEventCut->fUseVariablesCorrelationCuts = true;
        fEventCut->fUseStrongVarCorrelationCut = true;
    }

    if ( _singlesOnly )   _outputHistoList -> Add( fHelperPID -> GetOutputList() ); // add AliHelperPIDBFDptDpt object output list to task output list only for singles
    
    _nBins_M0 = 2000; _min_M0   = 0.;    _max_M0    = 2000.;  _width_M0 = (_max_M0-_min_M0)/_nBins_M0;
    _nBins_M1 = 2000; _min_M1   = 0.;    _max_M1    = 2000.;  _width_M1 = (_max_M1-_min_M1)/_nBins_M1;
    _nBins_M2 = 2000; _min_M2   = 0.;    _max_M2    = 2000.;  _width_M2 = (_max_M2-_min_M2)/_nBins_M2;
    _nBins_M3 = 2000; _min_M3   = 0.;    _max_M3    = 2000.;  _width_M3 = (_max_M3-_min_M3)/_nBins_M3; 
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
    _TrackArray = new AliAODTrack * [arraySize];
    __n1_1_vsPt                   = getDoubleArray(_nBins_pt_1, 0.);
    __n1_1_vsPt_pdg               = getDoubleArray(_nBins_pt_1, 0.);
    __n1_1_vsPt_pdg_Weak          = getDoubleArray(_nBins_pt_1, 0.);
    __n1_1_vsPt_pdg_Weak_Material = getDoubleArray(_nBins_pt_1, 0.);
    __n1_1_vsPt_Weak              = getDoubleArray(_nBins_pt_1, 0.);
    __n1_1_vsPt_Material          = getDoubleArray(_nBins_pt_1, 0.);
    __n1_1_vsEtaPhi          = getDoubleArray(_nBins_etaPhi_1,    0.);
    __s1pt_1_vsEtaPhi        = getDoubleArray(_nBins_etaPhi_1,    0.);
    // __n1_1_vsZEtaPhiPt       = getFloatArray(_nBins_zEtaPhiPt_1,  0.);


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
        __n1_2_vsPt_pdg               = getDoubleArray(_nBins_pt_2, 0.);
        __n1_2_vsPt_pdg_Weak          = getDoubleArray(_nBins_pt_2, 0.);
        __n1_2_vsPt_pdg_Weak_Material = getDoubleArray(_nBins_pt_2, 0.);
        __n1_2_vsPt_Weak              = getDoubleArray(_nBins_pt_2, 0.);
        __n1_2_vsPt_Material          = getDoubleArray(_nBins_pt_2, 0.);
        __n1_2_vsEtaPhi          = getDoubleArray(_nBins_etaPhi_2,    0.);
        __s1pt_2_vsEtaPhi        = getDoubleArray(_nBins_etaPhi_2,    0.);
       
    }//if (_requestedCharge_2!=_requestedCharge_1)

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
    _title_m9     = "M_{9}";
    _title_m10     = "M_{10}";


  
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

    _title_AvgN2_12       = "#LT n_{2} #GT";
    _title_AvgSumPtPt_12  = "#LT #Sigma p_{t,1}p_{t,2} #GT";
    _title_AvgSumPtN_12   = "#LT #Sigma p_{t,1}N #GT";
    _title_AvgNSumPt_12   = "#LT N#Sigma p_{t,2} #GT";
    
  
    vsZ         = "_vsZ";
    vsM         = "_vsM";
    vsPt        = "_vsPt";
    vsPhi       = "_vsPhi";
    vsEta       = "_vsEta";
    vsEtaPhi    = "_vsEtaPhi";
    vsPtVsPt    = "_vsPtVsPt";
    pdg         = "_pdg";
    Weak        = "_Weak";
    Material    = "_Material";
    
    //---------pt-efficiency starts ---------------
    if (_usePtEff)
    {
        Int_t  iPt;
        Int_t iPt1;
      
        if (_hPtEff_1)
        {
            _correctionPtEff_1 = new Float_t[_nBins_pt_1];

            for (iPt=0,iPt1=1; iPt<_nBins_pt_1; iPt++, iPt1++)
            {
                 _correctionPtEff_1[iPt] = _hPtEff_1->GetBinContent(iPt1);
            }

        }//if (_hPtEff_1)

        else
        {
            AliError("AliAnalysisTaskSpherodistR2P2:: _hPtEff_1 is a null pointer.");
	        return;
        }

        if (!_sameFilter)
        {
            if (_hPtEff_2)
            {   
                _correctionPtEff_2 = new Float_t[_nBins_pt_2];

                for (iPt=0,iPt1=1; iPt<_nBins_pt_2; iPt++, iPt1++)
                {
                    _correctionPtEff_2[iPt] = _hPtEff_2->GetBinContent(iPt1);
                }

            }//if (_hPtEff_2)

            else
            {
                AliError("AliAnalysisTaskSpherodistR2P2:: _hPtEff_2 is a null pointer.");
	            return;
            }

        }//if (!_sameFilter)

    }//if (_usePtEff)

    //-----------pt-efficiency ends------------
   
    fUtils = new AliAnalysisUtils();
  
    createHistograms();
  
    PostData(1,_outputHistoList);
  
    cout<< "AliAnalysisTaskSpherodistR2P2::CreateOutputObjects() DONE " << endl;
  

} //void AliAnalysisTaskSpherodistR2P2::UserCreateOutputObjects()


void  AliAnalysisTaskSpherodistR2P2::createHistograms()
{
    AliInfo("AliAnalysisTaskSpherodistR2P2::createHistoHistograms() Creating Event Histos");
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
    name = "m9"; _m9      = createHisto1D(name,name,_nBins_M2, _min_M2, _max_M2, _title_m9, _title_counts);
    name = "m10"; _m10      = createHisto1D(name,name,_nBins_M2, _min_M2, _max_M2, _title_m10, _title_counts);

    name = "m8"; _m8      = createHisto1D(name,name,_nBins_M8, _min_M8, _max_M8, _title_m8, _title_counts);
    name = "zV"; _vertexZ = createHisto1D(name,name,_nBins_vertexZ, _min_vertexZ, _max_vertexZ, "z-Vertex (cm)", _title_counts);
    name = "psi_EventPlane"; _psi_EventPlane = createHisto1F(name,name, 360, 0.0, 6.4, "#psi","counts");
    name = "V0MvsTracksTPCout_after";   _fhV0MvsTracksTPCout_after = createHisto2F(name,name,200,0,20000,400,0,40000, "# tracks with kTPCout on", "V0 multiplicity","counts"); //V0 multiplicity vs tracks with kTPCout on after cut

  //--------------------------------------------------------------------

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
        name = "msquare_p_POI_AliHelperPID";   _msquare_p_POI_AliHelperPID = createHisto2F(name,name, 500, 0, 5,  200, -0.5, 1.5,  "p", "mass square","counts");
        name = "msquare_p_AliHelperPID_no_Undefined";   _msquare_p_AliHelperPID_no_Undefined = createHisto2F(name,name, 500, 0, 5,  200, -0.5, 1.5,  "p", "mass square","counts");
        name = "nSigmaTOF_p_pion_before";   _nSigmaTOF_p_pion_before   = createHisto2F(name,name, 500, 0, 5, 120, -20, 100.0, "p", "n#sigma_{TOF}^{pion}","counts");
        name = "nSigmaTOF_p_pion_after";    _nSigmaTOF_p_pion_after    = createHisto2F(name,name, 500, 0, 5, 120, -20, 100.0, "p", "n#sigma_{TOF}^{pion}","counts");
        name = "nSigmaTOF_p_kaon_before";   _nSigmaTOF_p_kaon_before   = createHisto2F(name,name, 500, 0, 5, 120, -20, 100.0, "p", "n#sigma_{TOF}^{kaon}","counts");
        name = "nSigmaTOF_p_kaon_after";    _nSigmaTOF_p_kaon_after    = createHisto2F(name,name, 500, 0, 5, 120, -20, 100.0, "p", "n#sigma_{TOF}^{kaon}","counts");
        name = "nSigmaTOF_p_proton_before"; _nSigmaTOF_p_proton_before = createHisto2F(name,name, 500, 0, 5, 120, -20, 100.0, "p", "n#sigma_{TOF}^{proton}","counts");
        name = "nSigmaTOF_p_proton_after";  _nSigmaTOF_p_proton_after  = createHisto2F(name,name, 500, 0, 5, 120, -20, 100.0, "p", "n#sigma_{TOF}^{proton}","counts");
        name = "nSigmaTPC_nSigmaTOF_Pion_before";   _nSigmaTPC_nSigmaTOF_Pion_before   = createHisto2F(name,name,1000,-5,5,1000,-5,5,"n#sigma_{TPC}^{pion}","n#sigma_{TOF}^{pion}","counts");
        name = "nSigmaTPC_nSigmaTOF_Pion_after";    _nSigmaTPC_nSigmaTOF_Pion_after    = createHisto2F(name,name,1000,-5,5,1000,-5,5,"n#sigma_{TPC}^{pion}","n#sigma_{TOF}^{pion}","counts");
        name = "nSigmaTPC_nSigmaTOF_Kaon_before";   _nSigmaTPC_nSigmaTOF_Kaon_before   = createHisto2F(name,name,1000,-5,5,1000,-5,5,"n#sigma_{TPC}^{kaon}","n#sigma_{TOF}^{kaon}","counts");
        name = "nSigmaTPC_nSigmaTOF_Kaon_after";    _nSigmaTPC_nSigmaTOF_Kaon_after    = createHisto2F(name,name,1000,-5,5,1000,-5,5,"n#sigma_{TPC}^{kaon}","n#sigma_{TOF}^{kaon}","counts");
        name = "nSigmaTPC_nSigmaTOF_Proton_before"; _nSigmaTPC_nSigmaTOF_Proton_before = createHisto2F(name,name,1000,-5,5,1000,-5,5,"n#sigma_{TPC}^{proton}","n#sigma_{TOF}^{proton}","counts");
        name = "nSigmaTPC_nSigmaTOF_Proton_after";  _nSigmaTPC_nSigmaTOF_Proton_after  = createHisto2F(name,name,1000,-5,5,1000,-5,5,"n#sigma_{TPC}^{proton}","n#sigma_{TOF}^{proton}","counts");

        name = "dedx_p_POI_AliHelperPID";   _dedx_p_POI_AliHelperPID = createHisto2F(name,name, 1980, 0.2, 20,  200, 0, 200,  "p", "dedx","counts");
        name = "dedx_p_AliHelperPID_no_Undefined";   _dedx_p_AliHelperPID_no_Undefined = createHisto2F(name,name, 1980, 0.2, 20,  200, 0, 200,  "p", "dedx","counts");    
        name = "beta_p_POI_AliHelperPID";   _beta_p_POI_AliHelperPID = createHisto2F(name,name, 500, 0, 5,  100, 0.1, 1.1,  "p", "beta","counts");
        name = "beta_p_AliHelperPID_no_Undefined";   _beta_p_AliHelperPID_no_Undefined = createHisto2F(name,name, 500, 0, 5,  100, 0.1, 1.1,  "p", "beta","counts");
        name = "inverse_beta_p_POI_AliHelperPID";   _inverse_beta_p_POI_AliHelperPID = createHisto2F(name,name, 500, 0, 5, 200, 0.6, 2.6,  "p", "1/#beta","counts");
        name = "inverse_beta_p_AliHelperPID_no_Undefined";   _inverse_beta_p_AliHelperPID_no_Undefined = createHisto2F(name,name, 500, 0, 5, 200, 0.6, 2.6,  "p", "1/#beta","counts");
      
        name = "y_Pt_AllCh_MCAODTruth";       _y_Pt_AllCh_MCAODTruth = createHisto2F( name, name, _nBins_eta_1, _min_eta_1, _max_eta_1, _nBins_pt_1, _min_pt_1, _max_pt_1, "y", "p_{T}", "counts" );
        name = "y_Pt_POI_MCAODTruth";        _y_Pt_POI_MCAODTruth = createHisto2F( name, name, _nBins_eta_1, _min_eta_1, _max_eta_1, _nBins_pt_1, _min_pt_1, _max_pt_1, "y", "p_{T}", "counts" );      
      
        name = n1Name + part_1_Name + vsPt;                         _n1_1_vsPt = createHisto1F( name, name, _nBins_pt_1, _min_pt_1, _max_pt_1, _title_pt_1, _title_AvgN_1 );
        name = n1Name + part_1_Name + vsPt + pdg;                   _n1_1_vsPt_pdg = createHisto1F( name, name, _nBins_pt_1, _min_pt_1, _max_pt_1, _title_pt_1, _title_AvgN_1 );
        name = n1Name + part_1_Name + vsPt + pdg + Weak;            _n1_1_vsPt_pdg_Weak = createHisto1F( name, name, _nBins_pt_1, _min_pt_1, _max_pt_1, _title_pt_1, _title_AvgN_1 );
        name = n1Name + part_1_Name + vsPt + pdg + Weak + Material; _n1_1_vsPt_pdg_Weak_Material = createHisto1F( name, name, _nBins_pt_1, _min_pt_1, _max_pt_1, _title_pt_1, _title_AvgN_1 );
        name = n1Name + part_1_Name + vsPt + Weak;                  _n1_1_vsPt_Weak = createHisto1F( name, name, _nBins_pt_1, _min_pt_1, _max_pt_1, _title_pt_1, _title_AvgN_1 );
        name = n1Name + part_1_Name + vsPt + Material;              _n1_1_vsPt_Material = createHisto1F( name, name, _nBins_pt_1, _min_pt_1, _max_pt_1, _title_pt_1, _title_AvgN_1 );
      
        name = n1Name + part_2_Name + vsPt;                         _n1_2_vsPt = createHisto1F( name, name, _nBins_pt_2, _min_pt_2, _max_pt_2, _title_pt_2, _title_AvgN_2 );
        name = n1Name + part_2_Name + vsPt + pdg;                   _n1_2_vsPt_pdg = createHisto1F( name, name, _nBins_pt_2, _min_pt_2, _max_pt_2, _title_pt_2, _title_AvgN_2 );
        name = n1Name + part_2_Name + vsPt + pdg + Weak;            _n1_2_vsPt_pdg_Weak = createHisto1F( name, name, _nBins_pt_2, _min_pt_2, _max_pt_2, _title_pt_2, _title_AvgN_2 );
        name = n1Name + part_2_Name + vsPt + pdg + Weak + Material; _n1_2_vsPt_pdg_Weak_Material = createHisto1F( name, name, _nBins_pt_2, _min_pt_2, _max_pt_2, _title_pt_2, _title_AvgN_2 );
        name = n1Name + part_2_Name + vsPt + Weak;                  _n1_2_vsPt_Weak = createHisto1F( name, name, _nBins_pt_2, _min_pt_2, _max_pt_2, _title_pt_2, _title_AvgN_2 );
        name = n1Name + part_2_Name + vsPt + Material;              _n1_2_vsPt_Material = createHisto1F( name, name, _nBins_pt_2, _min_pt_2, _max_pt_2, _title_pt_2, _title_AvgN_2 );
      
        
    }//if ( _singlesOnly )

    else
    {
        name = n1Name + part_1_Name + vsPt;                         _n1_1_vsPt = createHisto1F( name, name, _nBins_pt_1, _min_pt_1, _max_pt_1, _title_pt_1, _title_AvgN_1 );
        name = n1Name + part_2_Name + vsPt;                         _n1_2_vsPt = createHisto1F( name, name, _nBins_pt_2, _min_pt_2, _max_pt_2, _title_pt_2, _title_AvgN_2 );
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
        name = "mInvKaon";   _invMassKaon   = createHisto1F(name,name, 80, 0.98, 1.06, "M_{KK}","counts");
        name = "mInvKaonSq"; _invMassKaonSq = createHisto1F(name,name, 120, 0.98, 1.10, "M_{KK}^2","counts");
        name = "mInvElec"; _invMassElec = createHisto1F(name,name, 500, 0., 1.000, "M_{inv}","counts");
        name = "ClusterSharedFraction_beforeCut"; _ClusterSharedFraction_beforeCut = createHisto1F(name,name, 400, 0., 1.0, "Cluster Shared Fraction","Pair Counts");
        name = "ClusterSharedFraction_afterCut"; _ClusterSharedFraction_afterCut = createHisto1F(name,name, 400, 0., 1.0, "Cluster Shared Fraction after Cut","Pair Counts");
        name = "ClusterSharedFraction_3by3Bins_beforeCut"; _ClusterSharedFraction_3by3Bins_beforeCut = createHisto1F(name,name, 400, 0., 1.0, "Cluster Shared Fraction","Pair Counts");
        name = "ClusterSharedFraction_3by3Bins_afterCut"; _ClusterSharedFraction_3by3Bins_afterCut = createHisto1F(name,name, 400, 0., 1.0, "Cluster Shared Fraction after Cut","Pair Counts");
    }//ELSE

    AliInfo(" AliAnalysisTaskSpherodistR2P2::createHistoHistograms() All Done");
}//void  AliAnalysisTaskSpherodistR2P2::createHistograms()


void  AliAnalysisTaskSpherodistR2P2::finalizeHistograms()
{
    AliInfo("\n\n\n\n\n\n AliAnalysisTaskSpherodistR2P2::finalizeHistograms() starting");

    AliInfo(Form("CorrelationAnalyzers::finalizeHistograms()   _eventCount : %d",int(_eventCount)));
    if (_singlesOnly)
    {
        if (_sameFilter)
        {
            fillHistoWithArray(_n1_1_vsPt,              __n1_1_vsPt,        _nBins_pt_1);
	        fillHistoWithArray(_n1_1_vsPt_pdg,          __n1_1_vsPt_pdg,    _nBins_pt_1);
            fillHistoWithArray(_n1_1_vsPt_pdg_Weak,            __n1_1_vsPt_pdg_Weak,            _nBins_pt_1);
            fillHistoWithArray(_n1_1_vsPt_pdg_Weak_Material,   __n1_1_vsPt_pdg_Weak_Material,   _nBins_pt_1);
            fillHistoWithArray(_n1_1_vsPt_Weak,                __n1_1_vsPt_Weak,                _nBins_pt_1);
            fillHistoWithArray(_n1_1_vsPt_Material,            __n1_1_vsPt_Material,            _nBins_pt_1);
	        
            fillHistoWithArray(_n1_2_vsPt,              __n1_1_vsPt,        _nBins_pt_1);
	        fillHistoWithArray(_n1_2_vsPt_pdg,          __n1_1_vsPt_pdg,    _nBins_pt_1);
            fillHistoWithArray(_n1_2_vsPt_pdg_Weak,            __n1_1_vsPt_pdg_Weak,            _nBins_pt_1);
            fillHistoWithArray(_n1_2_vsPt_pdg_Weak_Material,   __n1_1_vsPt_pdg_Weak_Material,   _nBins_pt_1);
            fillHistoWithArray(_n1_2_vsPt_Weak,                __n1_1_vsPt_Weak,                _nBins_pt_1);
            fillHistoWithArray(_n1_2_vsPt_Material,            __n1_1_vsPt_Material,            _nBins_pt_1);
	        
        }//if (_sameFilter)

        else
        {
            fillHistoWithArray(_n1_1_vsPt,              __n1_1_vsPt,        _nBins_pt_1);
	        fillHistoWithArray(_n1_1_vsPt_pdg,          __n1_1_vsPt_pdg,    _nBins_pt_1);
            fillHistoWithArray(_n1_1_vsPt_pdg_Weak,            __n1_1_vsPt_pdg_Weak,            _nBins_pt_1);
            fillHistoWithArray(_n1_1_vsPt_pdg_Weak_Material,   __n1_1_vsPt_pdg_Weak_Material,   _nBins_pt_1);
            fillHistoWithArray(_n1_1_vsPt_Weak,                __n1_1_vsPt_Weak,                _nBins_pt_1);
            fillHistoWithArray(_n1_1_vsPt_Material,            __n1_1_vsPt_Material,            _nBins_pt_1);
	        
            fillHistoWithArray(_n1_2_vsPt,              __n1_2_vsPt,        _nBins_pt_2);
	        fillHistoWithArray(_n1_2_vsPt_pdg,          __n1_2_vsPt_pdg,    _nBins_pt_2);
            fillHistoWithArray(_n1_2_vsPt_pdg_Weak,            __n1_2_vsPt_pdg_Weak,            _nBins_pt_2);
            fillHistoWithArray(_n1_2_vsPt_pdg_Weak_Material,   __n1_2_vsPt_pdg_Weak_Material,   _nBins_pt_2);
            fillHistoWithArray(_n1_2_vsPt_Weak,                __n1_2_vsPt_Weak,                _nBins_pt_2);
            fillHistoWithArray(_n1_2_vsPt_Material,            __n1_2_vsPt_Material,            _nBins_pt_2);
	        
        }//ELSE


    }//if (_singlesOnly)

    else
    {
        if (_sameFilter)
        {
            fillHistoWithArray(_n1_1_vsPt,              __n1_1_vsPt,        _nBins_pt_1);
	        fillHistoWithArray(_n1_2_vsPt,              __n1_1_vsPt,        _nBins_pt_1);
	        fillHistoWithArray(_n1_1_vsEtaVsPhi,        __n1_1_vsEtaPhi,    _nBins_eta_1,   _nBins_phi_1);
	        fillHistoWithArray(_s1pt_1_vsEtaVsPhi,      __s1pt_1_vsEtaPhi,  _nBins_eta_1,   _nBins_phi_1);
	        fillHistoWithArray(_n1_2_vsEtaVsPhi,        __n1_1_vsEtaPhi,    _nBins_eta_1,   _nBins_phi_1);
	        fillHistoWithArray(_s1pt_2_vsEtaVsPhi,      __s1pt_1_vsEtaPhi,  _nBins_eta_1,   _nBins_phi_1);

        }//if (_sameFilter)

        else
        {

            fillHistoWithArray(_n1_1_vsPt,              __n1_1_vsPt,        _nBins_pt_1);
	        fillHistoWithArray(_n1_2_vsPt,              __n1_2_vsPt,        _nBins_pt_2);
	        fillHistoWithArray(_n1_1_vsEtaVsPhi,        __n1_1_vsEtaPhi,    _nBins_eta_1,   _nBins_phi_1);
	        fillHistoWithArray(_s1pt_1_vsEtaVsPhi,      __s1pt_1_vsEtaPhi,  _nBins_eta_1,   _nBins_phi_1);
	        fillHistoWithArray(_n1_2_vsEtaVsPhi,        __n1_2_vsEtaPhi,    _nBins_eta_2,   _nBins_phi_2);
	        fillHistoWithArray(_s1pt_2_vsEtaVsPhi,      __s1pt_2_vsEtaPhi,  _nBins_eta_2,   _nBins_phi_2);

        }//else

        fillHistoWithArray(_n2_12_vsEtaPhi,     __n2_12_vsEtaPhi,     _nBins_etaPhi_12);
        fillHistoWithArray(_s2PtPt_12_vsEtaPhi, __s2ptpt_12_vsEtaPhi, _nBins_etaPhi_12);
        fillHistoWithArray(_s2PtN_12_vsEtaPhi,  __s2PtN_12_vsEtaPhi,  _nBins_etaPhi_12);
        fillHistoWithArray(_s2NPt_12_vsEtaPhi,  __s2NPt_12_vsEtaPhi,  _nBins_etaPhi_12);
        fillHistoWithArray(_n2_12_vsPtVsPt,     __n2_12_vsPtPt,       _nBins_pt_1,    _nBins_pt_2);

    }//ELSE

    AliInfo("AliAnalysisTaskSpherodistR2P2::finalizeHistograms()  Done ");

}//void  AliAnalysisTaskSpherodistR2P2::finalizeHistograms()



void  AliAnalysisTaskSpherodistR2P2::UserExec(Option_t */*option*/)
{

    vector<Double_t> nx;
    vector<Double_t> ny;
    vector<Double_t> pex;
    vector<Double_t> pey;
    Double_t spherocity = -999.;
    const Double_t pi2by4 = (TMath::Pi()*TMath::Pi())/4.0;
    Int_t nParts;

    int    k1 = 0, k2 = 0;
    int    iPhi=0, iEta=0, iEtaPhi=0, iPt=0, charge=0;
    int    IDrec=0;
    float  q=0., phi=0., pt=0., eta=0., y=0., y_direct=0., mass=0., corr=0., corrPt=0., px=0., py=0., pz=0., dedx=0.,p=0.,l=0., timeTOF=0., beta=0., t0=0., msquare=0.; // Au-Au put p here to make _dedx_p and _beta_p plots.
    Double_t  chi2ndf = 0.;
    Int_t countMult = 0;
   
   



    int    ij=0;
    int    id_1=0, q_1=0, iEtaPhi_1=0, iPt_1=0;
    float  pt_1=0., px_1=0., py_1=0., pz_1=0., corr_1=0., dedx_1=0.;
    int    id_2=0, q_2=0, iEtaPhi_2=0, iPt_2=0;
    float  pt_2=0., px_2=0., py_2=0., pz_2=0., corr_2=0., dedx_2=0.;
    float  ptpt=0.;
    int    iVertex=0, iVertexP1=0, iVertexP2=0;
    int    iZEtaPhiPt=0;
    float  massElecSq = 1.94797849000000016e-02;
    const  AliAODVertex*	vertex;
    int    nClus=0;
    bool   bitOK;
    const float mpion   = 0.139570; // GeV/c2
    const float mkaon   = 0.493677; // GeV/c2
    const float massKaonSq = 0.2437169803; // GeV/c2
    const float mproton = 0.938272; // GeV/c2
    Double_t c = TMath::C() * 1.E-9;// m/ns
    double EP = 0;

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
    Double_t nsigmaElectron = 999.;

    AliAnalysisManager* manager = AliAnalysisManager::GetAnalysisManager();
    if ( !manager ) 
    { 
        return; 
    }

    AliAODInputHandler* inputHandler = dynamic_cast<AliAODInputHandler*> (manager->GetInputEventHandler());

    if ( !inputHandler )
    {
        return; 
    }

    AliVEvent* event = (AliVEvent*)InputEvent();
    if (!event)
    {
        Printf("ERROR: Could not retrieve event");
        return;
    }

    if ( fSystemType == "pp18_V0_kMB_kFALSE")
    {
        
        fAODEvent = dynamic_cast<AliAODEvent*>(event); // create pointer to event

    }
    
    if ( !fAODEvent )
    {
        return;
    }

    fPIDResponse = inputHandler -> GetPIDResponse();
    if (!fPIDResponse)
    {
        AliFatal("This Task needs the PID response attached to the inputHandler");
        return;
    }

    if ( fSystemType == "pp18_V0_kMB_kFALSE")
    {
        //Discard incomplete events
        if(event->IsIncompleteDAQ()) return;

        if(! fPPVsMultUtils )
        {
            fPPVsMultUtils = new AliPPVsMultUtils();
        }

        
        Double_t vtx0[3] = {0., 0., 0.}; // don't rely on AOD vertex, assume (0,0,0)                        
      
        Double_t vtxBest[3];

        vtxBest[0] = event->GetPrimaryVertex()->GetX();
        vtxBest[1] = event->GetPrimaryVertex()->GetY();
        vtxBest[2] = event->GetPrimaryVertex()->GetZ();

        if (event->GetPrimaryVertex() && event->GetPrimaryVertex()->GetNContributors()>0.5)
        {
	        if (!(TMath::Abs(event->GetPrimaryVertex()->GetZ()) < _max_vertexZ) ) return;
            if ( _pileUpEventPP13 )
            {
                
                if (event->IsPileupFromSPD(3, 0.8, 3., 2., 5.)) return;
                if (event->IsPileupFromSPDInMultBins()) return;
                if(fUtils->IsPileUpMV(event)) return;
            }
            if(!(fPPVsMultUtils->IsINELgtZERO(event))) return;
            
            if(!(fPPVsMultUtils->IsAcceptedVertexPosition(event))) return;
            


        }//if (event->GetPrimaryVertex()

    }//if ( fSystemType == "pp18_V0_kMB_kFALSE")

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
    Int_t id = -999;
    Int_t label = -999;
    Double_t DCAX = -999;
    Double_t DCAY = -999;
    Double_t DCAZ = -999;
    Double_t DCAXY= -999;

    if( fAODEvent )
    {
        if ( fSystemType == "PbPb_2015_kTRUE" || fSystemType == "PbPb_2015_kFALSE")
        {
            if (!fEventCut->AcceptEvent(fAODEvent))
            {
                PostData(1, _outputHistoList);
	            return;
            }

            StoreEventMultiplicities( fAODEvent );
	        f2015V0MtoTrkTPCout = new TFormula(Form("f2015V0MtoTrkTPCout_%s",""),"-1250+3.125*x");
	        f2015V0MtoTrkTPCout_Upper = new TFormula(Form("f2015V0MtoTrkTPCout_Upper_%s",""),"2000+3.1666667*x");
	        //cout << "fV0Multiplicity = " << fV0Multiplicity <<endl;
	        //cout << "fV0Multiplicity_Victor = " << fV0Multiplicity_Victor <<endl;
	        //cout << "fNoOfTPCoutTracks = " << fNoOfTPCoutTracks <<endl;
	        if(Is2015PileUpEvent()) return;
            
            _fhV0MvsTracksTPCout_after->Fill(fNoOfTPCoutTracks, fV0Multiplicity);
        }//if ( fSystemType == "PbPb_2015_kTRUE" || fSystemType == "PbPb_2015_kFALSE")

        _eventAccounting -> Fill( 1 ); // count all events afer "official" pile-up cut

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

        }//if( fSystemType == "PbPb" || fSystemType == "pPb" )


        else if ( fSystemType == "PbPb_2015_kTRUE" || fSystemType == "pp_V0A_kMB_kTRUE" || fSystemType == "pp_V0C_kMB_kTRUE" || fSystemType == "pp_V0_kMB_kTRUE" )
        {
            AliMultSelection *multSelection = (AliMultSelection*) fAODEvent->FindListObject("MultSelection");
	        //If get this warning, please check that the AliMultSelectionTask actually ran (before your task)
	        if (!multSelection)    AliWarning("MultSelection not found in input event");
            else
            {
                v0Centr  = multSelection->GetMultiplicityPercentile("V0M", kTRUE);
	            v0ACentr = multSelection->GetMultiplicityPercentile("V0A", kTRUE);
	            v0CCentr = multSelection->GetMultiplicityPercentile("V0C", kTRUE);
            }

        }//else if ( fSystemType == "PbPb_2015_kTRUE"...)

        else if ( fSystemType == "PbPb_2015_kFALSE" || fSystemType == "pp_V0A_kMB_kFALSE" || fSystemType == "pp_V0C_kMB_kFALSE" || fSystemType == "pp_V0_kMB_kFALSE" || fSystemType == "pp18_V0_kMB_kFALSE")

        {
            
            AliMultSelection *multSelection = (AliMultSelection*) fAODEvent->FindListObject("MultSelection");
	        //If get this warning, please check that the AliMultSelectionTask actually ran (before your task)
	        if (!multSelection)    AliWarning("MultSelection not found in input event");

            else
            {
                v0Centr  = multSelection->GetMultiplicityPercentile("V0M", kFALSE);
	            v0ACentr = multSelection->GetMultiplicityPercentile("V0A", kFALSE);
	            v0CCentr = multSelection->GetMultiplicityPercentile("V0C", kFALSE);
            }

        }//else if ( fSystemType == "PbPb_2015_kFALSE"...)


        else if ( fSystemType == "pp_V0A_kMB_Utils" || fSystemType == "pp_V0C_kMB_Utils" || fSystemType == "pp_V0_kMB_Utils" )
        {
            //remove Pile-Up events
	        fUtils->SetUseMVPlpSelection(kTRUE);
	        fUtils->SetUseOutOfBunchPileUp(kTRUE);
	        if(fUtils->IsPileUpEvent(fAODEvent))   return;
            _eventAccounting -> Fill( 4 ); // number of events ( no pile-up ) 

	        v0Centr  = fUtils->GetMultiplicityPercentile(fAODEvent,"V0MEq");
	        v0ACentr = fUtils->GetMultiplicityPercentile(fAODEvent,"V0AEq");
	        v0CCentr = fUtils->GetMultiplicityPercentile(fAODEvent,"V0CEq");

        }//else if ( fSystemType == "pp_V0A_kMB_Utils" ...)

        else return;

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
        }//switch ( _centralityMethod )

        if( fSystemType == "PbPb" )
        {
            if  ( centrality < _centralityMin || centrality > _centralityMax || fabs( v0Centr - trkCentr ) > 5.0 )  return;
        }

        else if ( fSystemType == "PbPb_2015_kTRUE" || fSystemType == "PbPb_2015_kFALSE" || fSystemType == "pPb" || fSystemType == "pp" || fSystemType == "pp_V0A_kMB_kTRUE" || fSystemType == "pp_V0A_kMB_kFALSE" || fSystemType == "pp_V0_kMB_kTRUE" || fSystemType == "pp_V0_kMB_kFALSE" || fSystemType == "pp18_V0_kMB_kFALSE" || fSystemType == "pp_V0C_kMB_kTRUE" || fSystemType == "pp_V0C_kMB_kFALSE" || fSystemType == "pp_V0A_kMB_Utils" || fSystemType == "pp_V0_kMB_Utils" || fSystemType == "pp_V0C_kMB_Utils" )
        {
            if  ( centrality < _centralityMin || centrality > _centralityMax )  return;
        }

        else return;
        _eventAccounting -> Fill( 2 ); // count all events with right centrality


        
        // filter on z and xy vertex
        vertex = ( AliAODVertex * ) fAODEvent -> GetPrimaryVertex();

        if( vertex )
        {
            Double32_t fCov[6];
	        vertex->GetCovarianceMatrix(fCov);
	        if(vertex->GetNContributors() > 0.5)
            {
                if(fCov[5] != 0)
                {
                    vertexX = vertex->GetX();
		            vertexY = vertex->GetY();
		            vertexZ = vertex->GetZ();
                    if( TMath::Abs( vertexZ ) > _max_vertexZ )	return;  // Z-Vertex Cut
                }
            }
        }//if( vertex )


        iVertex = int(( vertexZ - _min_vertexZ ) / _width_vertexZ );
        iVertexP1 = iVertex*_nBins_etaPhiPt_1;
        iVertexP2 = iVertex*_nBins_etaPhiPt_2;
        if ( iVertex<0 || iVertex >= _nBins_vertexZ )
        {
            AliError("AliAnalysisTaskSpherodistR2P2::Exec(Option_t *option) iVertex<0 || iVertex>=_nBins_vertexZ ");

	        return;
        }
        _eventAccounting -> Fill( 3 ); // count all events with right centrality & vertexZ



        //*************************************************************************************************************//
        //                                                                                                            //
        ///////////////////////                 To include the spherocity here                //////////////////////////
        //                                                                                                            //
        //*************************************************************************************************************//
    
        if ( fSystemType == "pp18_V0_kMB_kFALSE")
        {
            if ( fAnalysisType == "RealData" || fAnalysisType == "MCAODreco" )
            {   

                //cout<<"Number of tracks for AOD reco ="<<_nTracks<<endl;


                nParts = 0;       

                Int_t countMult1 = 0;
                Double_t sumpt=0;

                //Track loop starts
                for (int iTrack = 0; iTrack < _nTracks; iTrack++ )
                {   

                    AliAODTrack * t = dynamic_cast<AliAODTrack *>( fAODEvent -> GetTrack( iTrack ) );

                    if ( !t )
                    {
                        AliError( Form( "Could not receive track %d", iTrack ) );
                        continue;
                    }

                    bitOK  = t -> TestFilterBit( _trackFilterBit );
                    if ( !bitOK ) continue;

                    if ( fAnalysisType == "RealData")
                    {
                        q      = t -> Charge();
                        charge = int( q );
                        phi    = t -> Phi();
                        p      = t -> P(); // momentum magnitude
                        pt     = t -> Pt();
                        px     = t -> Px();
                        py     = t -> Py();
                        pz     = t -> Pz();
                        eta    = t -> Eta();

                    }//if ( fAnalysisType == "RealData")

                    else if(fAnalysisType == "MCAODreco" )
                    {
                        
                        Int_t lab =  t -> GetLabel();
                        TClonesArray * arr = dynamic_cast<TClonesArray*>( fAODEvent -> FindListObject( "mcparticles" ) );
                        if( !arr ) continue;
                        AliAODMCParticle * particle = ( AliAODMCParticle * ) arr -> At( TMath::Abs(lab) );

                        if(!particle-> IsPhysicalPrimary()) continue;
                        //if(!particle->IsPrimary()) continue;
                        

                        q      =particle  -> Charge();
                        charge = int( q );
                        phi    = particle -> Phi();
                        p      = particle -> P(); // momentum magnitude
                        pt     = particle -> Pt();
                        px     = particle -> Px();
                        py     = particle -> Py();
                        pz     = particle -> Pz();
                        eta    = particle -> Eta();

                    }//else if (fAnalysisType == "MCAODreco" )

                        

                    if( q == 0 ) continue;

                        
                    if((pt >= 0.15) && (fabs(eta) <= 0.8))
                    {
                        

                        if(phi < 0) phi += 2.0 * TMath::Pi();  

                        pt=1;
                
                        sumpt += pt;

                        nx.push_back(TMath::Cos(phi));                                                                                               
                        ny.push_back(TMath::Sin(phi));                                                                                               
                        pex.push_back(pt*nx[countMult1]);                                                                                              
                        pey.push_back(pt*ny[countMult1]); 

                        ++countMult1;
                        
                    }

                }////Track loop ends

                //cout<<"Number of tracks for AOD reco ="<<nParts<< " "<<nParts1<< endl;

                if(countMult1 >= _NchMin)
                {
                    spherocity=1;
                    for(Int_t iMult = 0; iMult < countMult1; iMult++)
                    {
                        Double_t num = 0;

                        for(Int_t jMult = 0; jMult < countMult1; jMult++)
                        {
                        num += TMath::Abs( ny[iMult]*pex[jMult] - nx[iMult]*pey[jMult] );
                        }

                        if (sumpt==0)continue;
                        Double_t sFull = TMath::Power(( num / sumpt ),2 );
                        if(sFull <= spherocity) spherocity = sFull;

                    }//for(Int_t iMult = 0; iMult < countMult1; iMult++)

                    if( spherocity >= 0.0 && spherocity <= 1.0)
                    {
                        spherocity = spherocity * pi2by4; //normalization                                                                                        
                        hSpherocity->Fill(spherocity);
                    }

                }//if(countMult1 >= 5)

                    
                // h1i_n1_multSphero->Fill(countMult1);
                
            }//if ( fAnalysisType == "RealData" || fAnalysisType == "MCAODreco" )



           
            else if ( fAnalysisType == "MCAOD" )
            { 
                nParts = 0;

                AliMCEvent * mcEvent = MCEvent();
                _nTracks = mcEvent -> GetNumberOfTracks();                      

                Int_t countMult1 = 0;
                Double_t sumpt=0;


                for ( int iTrack = 0; iTrack < _nTracks; iTrack++ )

                {

                    AliAODMCParticle * t = ( AliAODMCParticle * ) mcEvent -> GetTrack( iTrack );
                    if ( !t )
                    {
                        AliError(Form("Could not receive track %d", iTrack));
                        continue;
                    }

                    if(!t-> IsPhysicalPrimary())  continue;

                    q      = t -> Charge();
	                charge = int( q/3. ); // particle charges are 3s in HIJING truth
                    phi    = t -> Phi();
		            p      = t -> P(); // momentum magnitude
		            pt     = t -> Pt();
		            px     = t -> Px();
		            py     = t -> Py();
		            pz     = t -> Pz();
		            eta    = t -> Eta();

                    if( q == 0 ) continue;

                    if((pt >= 0.15) && (fabs(eta) <= 0.8))
                    {
                    

                        if(phi < 0) phi += 2.0 * TMath::Pi();  

                        pt=1;
                    
                        sumpt += pt;

                        nx.push_back(TMath::Cos(phi));                                                                                               
                        ny.push_back(TMath::Sin(phi));                                                                                               
                        pex.push_back(pt*nx[countMult1]);                                                                                              
                        pey.push_back(pt*ny[countMult1]); 

                        ++countMult1;
                        
                    }

                }//Track loop 1

                //cout<<"Number of tracks for AOD gen ="<<nParts<<endl;

                
                if(countMult1 >= _NchMin)
                {
                    spherocity=1;

                    for(Int_t iMult = 0; iMult < countMult1; iMult++)
                    {
                        Double_t num = 0;

                        for(Int_t jMult = 0; jMult < countMult1; jMult++)
                        {
                            num += TMath::Abs( ny[iMult]*pex[jMult] - nx[iMult]*pey[jMult] );
                        }

                        if (sumpt==0)continue;
                        Double_t sFull = TMath::Power(( num / sumpt ),2 );
                        if(sFull <= spherocity) spherocity = sFull;

                    }//for(Int_t iMult = 0; iMult < countMult1; iMult++)

                    if (spherocity >= 0.0 ||spherocity <= 1.0)
                    {
                        spherocity = spherocity * pi2by4; //normalization                                                                                        
                        hSpherocity->Fill(spherocity);
                    }

                }//if(countMult1 >= 5)

            }//else if ( fAnalysisType == "MCAOD" )

        } //( fSystemType == "pp18_V0_kMB_kFALSE") spherocity
            

    }//if( fAODEvent )
    
    
    AliInfo("AliAnalysisTaskSpherodistR2P2::UserExec()   -----------------Event Done ");
    PostData(1,_outputHistoList);

}//void  AliAnalysisTaskSpherodistR2P2::UserExec(Option_t */*option*/)




void   AliAnalysisTaskSpherodistR2P2::FinishTaskOutput()
{
    AliInfo("AliAnalysisTaskSpherodistR2P2::FinishTaskOutput() Starting.");
    Printf("= 0 ====================================================================");
    finalizeHistograms();
    AliInfo("= 1 ====================================================================");
    PostData(1,_outputHistoList);
    AliInfo("= 2 ====================================================================");
    AliInfo("AliAnalysisTaskSpherodistR2P2::FinishTaskOutput() Done.");

}//void   AliAnalysisTaskSpherodistR2P2::FinishTaskOutput()


void   AliAnalysisTaskSpherodistR2P2::Terminate(Option_t* /*option*/)
{
  AliInfo("AliAnalysisTaskSpherodistR2P2::Terminate() Starting/Done.");
}


//===================================================================================================//
//                                             Tools                                                 //
//===================================================================================================//


void  AliAnalysisTaskSpherodistR2P2::fillHistoWithArray(TH1 * h, float * array, int size)
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


void  AliAnalysisTaskSpherodistR2P2::fillHistoWithArray(TH2 * h, float * array, int size1, int size2)
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


void  AliAnalysisTaskSpherodistR2P2::fillHistoWithArray(TH3 * h, float * array, int size1, int size2, int size3)
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

        }//for (j=0, j1=1; j<size2; ++j,++j1)
    }//for (i=0, i1=1; i<size1; ++i,++i1)


}//void  AliAnalysisTaskSpherodistR2P2::fillHistoWithArray(TH3 * h, float * array, int size1, int size2, int size3)



void  AliAnalysisTaskSpherodistR2P2::fillHistoWithArray(TH1 * h, double * array, int size)
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



void  AliAnalysisTaskSpherodistR2P2::fillHistoWithArray(TH2 * h, double * array, int size1, int size2)
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


void  AliAnalysisTaskSpherodistR2P2::fillHistoWithArray(TH3 * h, double * array, int size1, int size2, int size3)
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
double *  AliAnalysisTaskSpherodistR2P2::getDoubleArray(int size, double v)
{
  /// Allocate an array of type double with n values
  /// Initialize the array to the given value
  double * array = new double [size];
  for (int i=0;i<size;++i) array[i]=v;
  return array;
}

//________________________________________________________________________
float *  AliAnalysisTaskSpherodistR2P2::getFloatArray(int size, float v)
{
  /// Allocate an array of type float with n values
  /// Initialize the array to the given value
  float * array = new float [size];
  for (int i=0;i<size;++i) array[i]=v;
  return array;
}


//________________________________________________________________________
TH1D * AliAnalysisTaskSpherodistR2P2::createHisto1D(const TString &  name, const TString &  title, int n, double xMin, double xMax, const TString &  xTitle, const TString &  yTitle)
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
TH1D * AliAnalysisTaskSpherodistR2P2::createHisto1D(const TString &  name, const TString &  title, int n, double * bins, const TString &  xTitle, const TString &  yTitle)
{
  AliInfo(Form("createHisto 1D histo %s   with %d non uniform nBins",name.Data(),n));
  TH1D * h = new TH1D(name,title,n,bins);
  h->GetXaxis()->SetTitle(xTitle);
  h->GetYaxis()->SetTitle(yTitle);
  addToList(h);
  return h;
}


//________________________________________________________________________
TH2D * AliAnalysisTaskSpherodistR2P2::createHisto2D(const TString &  name, const TString &  title, int nx, double xMin, double xMax, int ny, double yMin, double yMax, const TString &  xTitle, const TString &  yTitle, const TString &  zTitle)
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
TH2D * AliAnalysisTaskSpherodistR2P2::createHisto2D(const TString &  name, const TString &  title, int nx, double* xbins, int ny, double yMin, double yMax, const TString &  xTitle, const TString &  yTitle, const TString &  zTitle)
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
TH1F * AliAnalysisTaskSpherodistR2P2::createHisto1F(const TString &  name, const TString &  title, int n, double xMin, double xMax, const TString &  xTitle, const TString &  yTitle)
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
TH1F * AliAnalysisTaskSpherodistR2P2::createHisto1F(const TString &  name, const TString &  title, int n, double * bins, const TString &  xTitle, const TString &  yTitle)
{
  AliInfo(Form("createHisto 1D histo %s   with %d non uniform nBins",name.Data(),n));
  TH1F * h = new TH1F(name,title,n,bins);
  h->GetXaxis()->SetTitle(xTitle);
  h->GetYaxis()->SetTitle(yTitle);
  addToList(h);
  return h;
}


//________________________________________________________________________
TH2F * AliAnalysisTaskSpherodistR2P2::createHisto2F(const TString &  name, const TString &  title, int nx, double xMin, double xMax, int ny, double yMin, double yMax,const TString &  xTitle, const TString &  yTitle, const TString &  zTitle)
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
TH2F * AliAnalysisTaskSpherodistR2P2::createHisto2F(const TString &  name, const TString &  title, int nx, double* xbins, int ny, double yMin, double yMax, const TString &  xTitle, const TString &  yTitle, const TString &  zTitle)
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
TH3F * AliAnalysisTaskSpherodistR2P2::createHisto3F(const TString &  name, const TString &  title, int nx, double xMin, double xMax, int ny, double yMin, double yMax, int nz, double zMin, double zMax, const TString &  xTitle, const TString &  yTitle, const TString &  zTitle)
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
TProfile * AliAnalysisTaskSpherodistR2P2::createProfile(const TString & name, const TString & description, int nx,double xMin,double xMax,const TString &  xTitle, const TString &  yTitle)
{
  AliInfo(Form("createHisto 1D profile %s   nBins: %d  xMin: %f10.4 xMax: %f10.4",name.Data(),nx,xMin,xMax));
  TProfile * h = new TProfile(name,description,nx,xMin,xMax);
  h->GetXaxis()->SetTitle(xTitle);
  h->GetYaxis()->SetTitle(yTitle);
  addToList(h);
  return h;  
}

//________________________________________________________________________
TProfile * AliAnalysisTaskSpherodistR2P2::createProfile(const TString &  name,const TString &  description, int nx,  double* bins,const TString &  xTitle, const TString &  yTitle)
{
  AliInfo(Form("createHisto 1D profile %s  with %d non uniform bins",name.Data(),nx));
  TProfile * h = new TProfile(name,description,nx,bins);
  h->GetXaxis()->SetTitle(xTitle);
  h->GetYaxis()->SetTitle(yTitle);
  addToList(h);
  return h;
}

void   AliAnalysisTaskSpherodistR2P2::addToList(TH1 *h)
{

    if (_outputHistoList)
    {
        _outputHistoList->Add(h);
    }
    else
    AliInfo("addToList(TH1 *h) _outputHistoList is null!!!!! Should abort ship");
}

void AliAnalysisTaskSpherodistR2P2::CalculateNSigmas( AliVTrack * trk )   //defines data member fnsigmas
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
void AliAnalysisTaskSpherodistR2P2::CalculateTPCNSigmasElectron( AliVTrack * trk )   //defines data member fnsigmas
{  
  // Compute nsigma for each hypthesis
  AliVParticle * inEvHMain = dynamic_cast < AliVParticle * > ( trk );
  // --- TPC
  Double_t nsigmaTPCkElectron = fPIDResponse -> NumberOfSigmasTPC( inEvHMain, AliPID::kElectron );
  //set data member fnsigmas
  fnsigmas[3][0] = nsigmaTPCkElectron;
}


//________________________________________________________________________________________________________
void AliAnalysisTaskSpherodistR2P2::CheckTOF( AliVTrack * trk )     //check if the particle has TOF Matching
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



Int_t AliAnalysisTaskSpherodistR2P2::TellParticleSpecies( AliVTrack * trk )  //function to determine the particle ID
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
      
      if( ( nsigmaPion < fNSigmaPID ) && ( nsigmaKaon > fNSigmaPID_veto ) && ( nsigmaProton > fNSigmaPID_veto ) && ( trk->Pt() <= ptUpperLimit ) )  return 0; //Pion
      if( ( nsigmaKaon < fNSigmaPID ) && ( nsigmaPion > fNSigmaPID_veto ) && ( nsigmaProton > fNSigmaPID_veto ) && ( trk->Pt() <= ptUpperLimit ) )  return 1; //Kaon
      if( ( nsigmaProton < fNSigmaPID ) && ( nsigmaPion > fNSigmaPID_veto ) && ( nsigmaKaon > fNSigmaPID_veto ) && ( trk->Pt() <= ptUpperLimit ) && ( massSquareCalculation(trk) > 0.6 ) && ( massSquareCalculation(trk) < 1.1 ) )   return 2; //Proton // need to add mass square cut in the source code as well!!!
    }
  else
    {
      nsigmaPion     =  TMath::Abs( fnsigmas[0][0] ); //Pion_TPC
      nsigmaKaon     =  TMath::Abs( fnsigmas[1][0] ); //Kaon_TPC
      nsigmaProton   =  TMath::Abs( fnsigmas[2][0] ); //Proton_TPC
      nsigmaElectron =  TMath::Abs( fnsigmas[3][0] ); //Electron_TPC
      
      if( ( nsigmaPion < fNSigmaPID ) && ( nsigmaKaon > fNSigmaPID_veto ) && ( nsigmaProton > fNSigmaPID_veto ) && ( nsigmaElectron > electronNSigmaVeto ) && ( trk->Pt() >= _min_pt_1 ) && ( trk->Pt() <= ptTOFlowerBoundary ) )   return 0; //Pion
      if( ( nsigmaKaon < fNSigmaPID ) && ( nsigmaPion > fNSigmaPID_veto ) && ( nsigmaProton > fNSigmaPID_veto ) && ( nsigmaElectron > electronNSigmaVeto ) && ( trk->Pt() >= _min_pt_1 ) && ( trk->Pt() <= ptTOFlowerBoundary ) )   return 1; //Kaon
      if( ( nsigmaProton < fNSigmaPID ) && ( nsigmaPion > fNSigmaPID_veto ) && ( nsigmaKaon > fNSigmaPID_veto ) && ( nsigmaElectron > electronNSigmaVeto ) && ( trk->Pt() >= _min_pt_1 ) && ( trk->Pt() <= ptTOFlowerBoundary ) )   return 2; //Proton
    }
  
  // else, return undefined
  return 999;
}




Int_t AliAnalysisTaskSpherodistR2P2::TellParticleSpecies_CircularCut( AliVTrack * trk )  //function to determine the particle ID
{  
  CalculateNSigmas( trk );

  Double_t nsigmaPion = 999., nsigmaKaon = 999., nsigmaProton = 999., nsigmaElectron = 999.;
  
  CheckTOF( trk );
  
  if( fHasTOFPID && ( trk->Pt() > ptTOFlowerBoundary ) )
    {
      Double_t d2Pion   = fnsigmas[0][0]*fnsigmas[0][0] + fnsigmas[0][1]*fnsigmas[0][1];
      Double_t d2Kaon   = fnsigmas[1][0]*fnsigmas[1][0] + fnsigmas[1][1]*fnsigmas[1][1];
      Double_t d2Proton = fnsigmas[2][0]*fnsigmas[2][0] + fnsigmas[2][1]*fnsigmas[2][1];
      
      nsigmaPion    =  TMath::Sqrt(d2Pion); //Pion_TPC+TOF_circular
      nsigmaKaon    =  TMath::Sqrt(d2Kaon); //Kaon_TPC+TOF_circular
      nsigmaProton  =  TMath::Sqrt(d2Proton); //Proton_TPC+TOF_circular
      
      if( ( nsigmaPion < fNSigmaPID ) && ( nsigmaKaon > fNSigmaPID_veto ) && ( nsigmaProton > fNSigmaPID_veto ) && ( trk->Pt() <= ptUpperLimit ) )  return 0; //Pion
      if( ( nsigmaKaon < fNSigmaPID ) && ( nsigmaPion > fNSigmaPID_veto ) && ( nsigmaProton > fNSigmaPID_veto ) && ( trk->Pt() <= ptUpperLimit ) )  return 1; //Kaon
      if( ( nsigmaProton < fNSigmaPID ) && ( nsigmaPion > fNSigmaPID_veto ) && ( nsigmaKaon > fNSigmaPID_veto ) && ( trk->Pt() <= ptUpperLimit ) && ( massSquareCalculation(trk) > 0.6 ) && ( massSquareCalculation(trk) < 1.1 ) )   return 2; //Proton // need to add mass square cut in the source code as well!!!
    }
  else
    {
      nsigmaPion     =  TMath::Abs( fnsigmas[0][0] ); //Pion_TPC
      nsigmaKaon     =  TMath::Abs( fnsigmas[1][0] ); //Kaon_TPC
      nsigmaProton   =  TMath::Abs( fnsigmas[2][0] ); //Proton_TPC
      nsigmaElectron =  TMath::Abs( fnsigmas[3][0] ); //Electron_TPC
      
      if( ( nsigmaPion < fNSigmaPID ) && ( nsigmaKaon > fNSigmaPID_veto ) && ( nsigmaProton > fNSigmaPID_veto ) && ( nsigmaElectron > electronNSigmaVeto ) && ( trk->Pt() >= _min_pt_1 ) && ( trk->Pt() <= ptTOFlowerBoundary ) )   return 0; //Pion
      if( ( nsigmaKaon < fNSigmaPID ) && ( nsigmaPion > fNSigmaPID_veto ) && ( nsigmaProton > fNSigmaPID_veto ) && ( nsigmaElectron > electronNSigmaVeto ) && ( trk->Pt() >= _min_pt_1 ) && ( trk->Pt() <= ptTOFlowerBoundary ) )   return 1; //Kaon
      if( ( nsigmaProton < fNSigmaPID ) && ( nsigmaPion > fNSigmaPID_veto ) && ( nsigmaKaon > fNSigmaPID_veto ) && ( nsigmaElectron > electronNSigmaVeto ) && ( trk->Pt() >= _min_pt_1 ) && ( trk->Pt() <= ptTOFlowerBoundary ) )   return 2; //Proton
    }
  
  // else, return undefined
  return 999;
}

//________________________________________________________________________________________________________________
Int_t AliAnalysisTaskSpherodistR2P2::TellParticleSpecies_by_P( AliVTrack * trk )  //function to determine the particle ID
{  
  CalculateNSigmas( trk );
  
  Double_t nsigmaPion = 999., nsigmaKaon = 999., nsigmaProton = 999., nsigmaElectron = 999.;
  
  CheckTOF( trk );
  
  if( fHasTOFPID && ( trk->P() > ptTOFlowerBoundary ) ) // !!! diff from TellParticleSpecies
    {
      nsigmaPion     =  TMath::Abs( fnsigmas[0][1] ); //Pion_TOF
      nsigmaKaon     =  TMath::Abs( fnsigmas[1][1] ); //Kaon_TOF
      nsigmaProton   =  TMath::Abs( fnsigmas[2][1] ); //Proton_TOF
      nsigmaElectron =  TMath::Abs( fnsigmas[3][1] ); //Electron_TOF
      
      if( ( nsigmaPion < fNSigmaPID ) && ( nsigmaKaon > fNSigmaPID_veto ) && ( nsigmaProton > fNSigmaPID_veto ) && ( trk->Pt() <= ptUpperLimit ) )  return 0; //Pion
      if( ( nsigmaKaon < fNSigmaPID ) && ( nsigmaPion > fNSigmaPID_veto ) && ( nsigmaProton > fNSigmaPID_veto ) && ( trk->Pt() <= ptUpperLimit ) )  return 1; //Kaon
      if( ( nsigmaProton < fNSigmaPID ) && ( nsigmaPion > fNSigmaPID_veto ) && ( nsigmaKaon > fNSigmaPID_veto ) && ( trk->Pt() <= ptUpperLimit ) && ( massSquareCalculation(trk) > 0.6 ) && ( massSquareCalculation(trk) < 1.1 ) )   return 2; //Proton // need to add mass square cut in the source code as well!!!
    }
  else
    {
      nsigmaPion     =  TMath::Abs( fnsigmas[0][0] ); //Pion_TPC
      nsigmaKaon     =  TMath::Abs( fnsigmas[1][0] ); //Kaon_TPC
      nsigmaProton   =  TMath::Abs( fnsigmas[2][0] ); //Proton_TPC
      nsigmaElectron =  TMath::Abs( fnsigmas[3][0] ); //Electron_TPC
      
      if( ( nsigmaPion < fNSigmaPID ) && ( nsigmaKaon > fNSigmaPID_veto ) && ( nsigmaProton > fNSigmaPID_veto ) && ( nsigmaElectron > electronNSigmaVeto ) && ( trk->Pt() >= _min_pt_1 ) && ( trk->P() <= ptTOFlowerBoundary ) )   return 0; //Pion 
      if( ( nsigmaKaon < fNSigmaPID ) && ( nsigmaPion > fNSigmaPID_veto ) && ( nsigmaProton > fNSigmaPID_veto ) && ( nsigmaElectron > electronNSigmaVeto ) && ( trk->Pt() >= _min_pt_1 ) && ( trk->P() <= ptTOFlowerBoundary ) )   return 1; //Kaon
      if( ( nsigmaProton < fNSigmaPID ) && ( nsigmaPion > fNSigmaPID_veto ) && ( nsigmaKaon > fNSigmaPID_veto ) && ( nsigmaElectron > electronNSigmaVeto ) && ( trk->Pt() >= _min_pt_1 ) && ( trk->P() <= ptTOFlowerBoundary ) )   return 2; //Proton
    }
  
  // else, return undefined
  return 999;
}



//________________________________________________________________________________________________________________
Int_t AliAnalysisTaskSpherodistR2P2::TellParticleSpecies_by_P_CircularCut( AliVTrack * trk )  //function to determine the particle ID
{  
  CalculateNSigmas( trk );
  
  Double_t nsigmaPion = 999., nsigmaKaon = 999., nsigmaProton = 999., nsigmaElectron = 999.;
  
  CheckTOF( trk );
  
  if( fHasTOFPID && ( trk->P() > ptTOFlowerBoundary ) ) // !!! diff from TellParticleSpecies
    {
      Double_t d2Pion   = fnsigmas[0][0]*fnsigmas[0][0] + fnsigmas[0][1]*fnsigmas[0][1];
      Double_t d2Kaon   = fnsigmas[1][0]*fnsigmas[1][0] + fnsigmas[1][1]*fnsigmas[1][1];
      Double_t d2Proton = fnsigmas[2][0]*fnsigmas[2][0] + fnsigmas[2][1]*fnsigmas[2][1];
    
      nsigmaPion    =  TMath::Sqrt(d2Pion); //Pion_TPC+TOF_circular
      nsigmaKaon    =  TMath::Sqrt(d2Kaon); //Kaon_TPC+TOF_circular
      nsigmaProton  =  TMath::Sqrt(d2Proton); //Proton_TPC+TOF_circular

      if( ( nsigmaPion < fNSigmaPID ) && ( nsigmaKaon > fNSigmaPID_veto ) && ( nsigmaProton > fNSigmaPID_veto ) && ( trk->Pt() <= ptUpperLimit ) )  return 0; //Pion
      if( ( nsigmaKaon < fNSigmaPID ) && ( nsigmaPion > fNSigmaPID_veto ) && ( nsigmaProton > fNSigmaPID_veto ) && ( trk->Pt() <= ptUpperLimit ) )  return 1; //Kaon
      if( ( nsigmaProton < fNSigmaPID ) && ( nsigmaPion > fNSigmaPID_veto ) && ( nsigmaKaon > fNSigmaPID_veto ) && ( trk->Pt() <= ptUpperLimit ) && ( massSquareCalculation(trk) > 0.6 ) && ( massSquareCalculation(trk) < 1.1 ) )   return 2; //Proton // need to add mass square cut in the source code as well!!!
    }
  else
    {
      nsigmaPion     =  TMath::Abs( fnsigmas[0][0] ); //Pion_TPC
      nsigmaKaon     =  TMath::Abs( fnsigmas[1][0] ); //Kaon_TPC
      nsigmaProton   =  TMath::Abs( fnsigmas[2][0] ); //Proton_TPC
      nsigmaElectron =  TMath::Abs( fnsigmas[3][0] ); //Electron_TPC

      if( ( nsigmaPion < fNSigmaPID ) && ( nsigmaKaon > fNSigmaPID_veto ) && ( nsigmaProton > fNSigmaPID_veto ) && ( nsigmaElectron > electronNSigmaVeto ) && ( trk->Pt() >= _min_pt_1 ) && ( trk->P() <= ptTOFlowerBoundary ) )   return 0; //Pion 
      if( ( nsigmaKaon < fNSigmaPID ) && ( nsigmaPion > fNSigmaPID_veto ) && ( nsigmaProton > fNSigmaPID_veto ) && ( nsigmaElectron > electronNSigmaVeto ) && ( trk->Pt() >= _min_pt_1 ) && ( trk->P() <= ptTOFlowerBoundary ) )   return 1; //Kaon
      if( ( nsigmaProton < fNSigmaPID ) && ( nsigmaPion > fNSigmaPID_veto ) && ( nsigmaKaon > fNSigmaPID_veto ) && ( nsigmaElectron > electronNSigmaVeto ) && ( trk->Pt() >= _min_pt_1 ) && ( trk->P() <= ptTOFlowerBoundary ) )   return 2; //Proton
    }
  
  // else, return undefined
  return 999;
}



//________________________________________________________________________________________________________________
Double_t AliAnalysisTaskSpherodistR2P2::TOFBetaCalculation( AliVTrack * track ) const
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
Double_t AliAnalysisTaskSpherodistR2P2::massSquareCalculation( AliVTrack * track ) const
{
  Double_t massSquare = track->P() * track->P() * ( 1 / ( TOFBetaCalculation(track) * TOFBetaCalculation(track) ) - 1 ); 
  return massSquare;
}

//________________________________________________________________________________________________________________
Float_t AliAnalysisTaskSpherodistR2P2::TPC_EventPlane(AliAODEvent *fAOD)
{
  Int_t ntracks = fAOD->GetNumberOfTracks();
  
  TVector2 mQ;
  float mypsi=0;
  float mQx=0, mQy=0;
  
  for(Int_t i=0; i<ntracks; i++)
    {
      AliAODTrack *track =(AliAODTrack *)fAOD->GetTrack(i); // pointer to reconstructed to track
      if(!track)
        {
	        AliError(Form("ERROR: Could not retrieve aod track %d\n",i));
	        continue;
        }                        
      mQx += (cos(2*track->Phi()))*(track->Pt());
      mQy += (sin(2*track->Phi()))*(track->Pt());
    }
  
  mQ.Set(mQx,mQy);
  
  mypsi=(TMath::ATan2(-mQy,-mQx))/2.;
  //cout<<"mypsi_fct :"<<mypsi<<endl;
  return mypsi+TMath::Pi();
}


//________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskSpherodistR2P2::Is2015PileUpEvent()
{
  Bool_t IsPileUpEvent2015 = kTRUE;
  if ( (fV0Multiplicity<f2015V0MtoTrkTPCout->Eval(fNoOfTPCoutTracks)) || (fV0Multiplicity>f2015V0MtoTrkTPCout_Upper->Eval(fNoOfTPCoutTracks)) ) IsPileUpEvent2015 = kTRUE;
  else IsPileUpEvent2015 = kFALSE;
  return IsPileUpEvent2015;
}


Bool_t AliAnalysisTaskSpherodistR2P2::StoreEventMultiplicities(AliVEvent *event)
{  
  AliAODEvent *aodEvent=dynamic_cast<AliAODEvent*>(event);
  
  fNoOfTPCoutTracks = 0;
  //fV0Multiplicity_Victor = aodEvent->GetVZEROData()->GetMTotV0A()+event->GetVZEROData()->GetMTotV0C(); //Victor's way
  
  AliVVZERO *vzero = (AliVVZERO*)event->GetVZEROData();
  if(vzero)
    {
      fV0Multiplicity = 0;
      for(int ich=0; ich < 64; ich++)
      {
        fV0Multiplicity += vzero->GetMultiplicity(ich);
      } //ich loop   
      
    } // AliEventCuts
  
  Int_t nTracks = 0;
  nTracks = aodEvent->GetNumberOfTracks();
  
  for (Int_t itrk = 0; itrk < nTracks; itrk++)
    {
      AliAODTrack *aodt = dynamic_cast<AliAODTrack*>(aodEvent->GetTrack(itrk));
      if ((aodt->GetStatus() & AliVTrack::kTPCout) && aodt->GetID() > 0 )   fNoOfTPCoutTracks++;

      
      
    } 
  return kTRUE;
}



//________________________________________________________________________________________________________________
Double_t AliAnalysisTaskSpherodistR2P2::CalculateSharedFraction (const TBits *triggerClusterMap,const TBits *assocClusterMap,const TBits *triggerShareMap,const TBits *assocShareMap)
{
  Double_t nofhits=0;
  Double_t nofsharedhits=0;
  
  for(UInt_t imap=0;imap< (triggerClusterMap->GetNbits() );imap++)
    {
        //if they are in same pad
        //cout<<triggerClusterMap->TestBitNumber(imap)<<"    "<< assocClusterMap->TestBitNumber(imap)<<endl;
        if (triggerClusterMap->TestBitNumber(imap) && assocClusterMap->TestBitNumber(imap))
	    {
	        //if they share
	        //cout<<triggerShareMap->TestBitNumber(imap)<<"   "<<assocShareMap->TestBitNumber(imap)<<endl;
	        if (triggerShareMap->TestBitNumber(imap) && assocShareMap->TestBitNumber(imap))
	        { //cout<<triggerShareMap->TestBitNumber(imap)<<"   "<<assocShareMap->TestBitNumber(imap)<<endl;
	            nofhits+=2;
	            nofsharedhits+=2;
	        }
	        //not shared
	        else { nofhits+=2; }
	    }
      
        //different pad
        //cout<< (triggerClusterMap->TestBitNumber(imap) || assocClusterMap->TestBitNumber(imap))<<endl;
        else if (triggerClusterMap->TestBitNumber(imap) || assocClusterMap->TestBitNumber(imap))
	    {// One track has a hit, the other does not
	    nofhits++;
	    //cout<<"No hits :"<<nofhits<<endl;
	    }
    }
  
  Double_t SharedFraction=0.0;
  if(nofhits>0) SharedFraction=(nofsharedhits/nofhits);
  
  //cout<<"Fraction shared hits :"<<SharedFraction<<endl;
  
  return SharedFraction;
  
}
