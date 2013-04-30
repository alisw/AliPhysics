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

/* $Id:$ */

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
#include <TRandom.h>
#include "AliAnalysisManager.h"

#include "AliAODHandler.h"
#include "AliAODInputHandler.h"
//#include "AliAODMCParticle.h"
#include "AliInputEventHandler.h"
#include "AliLog.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliMultiplicity.h"
#include "AliCentrality.h"
#include "AliAnalysisTaskDptDptCorrelations.h"
using std::cout;
using std::endl;
using std::cerr;
ClassImp(AliAnalysisTaskDptDptCorrelations)

AliAnalysisTaskDptDptCorrelations::AliAnalysisTaskDptDptCorrelations()
: AliAnalysisTaskSE(),
fAODEvent(0), 
fESDEvent(0),             //! ESD Event 
fInputHandler(0),
_outputHistoList(0),
_twoPi         ( 2.0 * 3.1415927),
_eventCount    ( 0), 
_debugLevel    ( 0),
_singlesOnly   ( 0), 
_useWeights    ( 0), 
_sameFilter    ( false),
_rejectPileup  ( 1), 
_rejectPairConversion ( 0), 
_vertexZMin           ( -10), 
_vertexZMax           (  10), 
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
_nClusterMin          ( 70), 
_trackFilterBit       ( 128),
_field    ( 1.),
_nTracks  ( 0 ),
_mult0    ( 0 ),
_mult1    ( 0 ),
_mult2    ( 0 ),
_mult3    ( 0 ),
_mult4    ( 0 ),
_mult5    ( 0 ),
_mult6    ( 0 ),
arraySize ( 2000),
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
_nBins_vertexZ(40),   _min_vertexZ(-10), _max_vertexZ(10),        _width_vertexZ(0.5),

_nBins_pt_1(18),      _min_pt_1(0.2),    _max_pt_1(2.0),          _width_pt_1(0.1),
_nBins_phi_1(72),     _min_phi_1(0),     _max_phi_1(2.*3.1415927),_width_phi_1(2.*3.1415927/72.),
_nBins_eta_1(20),     _min_eta_1(-1),    _max_eta_1(1),           _width_eta_1(0.1),

_nBins_etaPhi_1(0), 
_nBins_etaPhiPt_1(0),
_nBins_zEtaPhiPt_1(0),

_nBins_pt_2(18),     _min_pt_2(0.2),     _max_pt_2(2.0),          _width_pt_2(0.1),
_nBins_phi_2(72),    _min_phi_2(0),      _max_phi_2(2.*3.1415927),_width_phi_2(2.*3.1415927/72),
_nBins_eta_2(20),    _min_eta_2(-1),     _max_eta_2(1),           _width_eta_2(0.1),

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

_etadis ( 0),
_phidis ( 0),
_dcaz   ( 0),
_dcaxy  ( 0),

_vertexZ           ( 0),
_n1_1_vsPt         ( 0),         
_n1_1_vsEtaVsPhi   ( 0),
_s1pt_1_vsEtaVsPhi ( 0), 
_n1_1_vsZVsEtaVsPhiVsPt ( 0),
_n1_1_vsM          ( 0),  // w/ weight
_s1pt_1_vsM        ( 0),
_n1Nw_1_vsM        ( 0), // w/o weight
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
vsPtVsPt("NA")
{
  printf("Default constructor called \n");
  
  printf("passed \n ");
  
}

AliAnalysisTaskDptDptCorrelations::AliAnalysisTaskDptDptCorrelations(const TString & name)
: AliAnalysisTaskSE(name),
fAODEvent(0), 
fESDEvent(0),  
fInputHandler(0),
_outputHistoList(0),
_twoPi         ( 2.0 * 3.1415927),
_eventCount    ( 0), 
_debugLevel    ( 0),
_singlesOnly   ( 0), 
_useWeights    ( 0), 
_sameFilter    ( false),
_rejectPileup  ( 1), 
_rejectPairConversion ( 0), 
_vertexZMin           ( -10.), 
_vertexZMax           (  10.), 
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
_nClusterMin          ( 70), 
_trackFilterBit       ( 128),
_field    ( 1.),
_nTracks  ( 0 ),
_mult0    ( 0 ),
_mult1    ( 0 ),
_mult2    ( 0 ),
_mult3    ( 0 ),
_mult4    ( 0 ),
_mult5    ( 0 ),
_mult6    ( 0 ),
arraySize ( 2000),
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
_nBins_vertexZ(40),   _min_vertexZ(-10), _max_vertexZ(10),        _width_vertexZ(0.5),

_nBins_pt_1(18),      _min_pt_1(0.2),    _max_pt_1(2.0),          _width_pt_1(0.1),
_nBins_phi_1(72),     _min_phi_1(0),     _max_phi_1(2.*3.1415927),_width_phi_1(2.*3.1415927/72.),
_nBins_eta_1(20),     _min_eta_1(-1),    _max_eta_1(1),           _width_eta_1(0.1),

_nBins_etaPhi_1(0), 
_nBins_etaPhiPt_1(0),
_nBins_zEtaPhiPt_1(0),

_nBins_pt_2(18),     _min_pt_2(0.2),     _max_pt_2(2.0),          _width_pt_2(0.1),
_nBins_phi_2(72),    _min_phi_2(0),      _max_phi_2(2.*3.1415927),_width_phi_2(2.*3.1415927/72),
_nBins_eta_2(20),    _min_eta_2(-1),     _max_eta_2(1),           _width_eta_2(0.1),

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
_etadis ( 0),
_phidis ( 0),
_dcaz ( 0),
_dcaxy ( 0),
_vertexZ           ( 0),
_n1_1_vsPt         ( 0),         
_n1_1_vsEtaVsPhi   ( 0),
_s1pt_1_vsEtaVsPhi ( 0), 
_n1_1_vsZVsEtaVsPhiVsPt ( 0),
_n1_1_vsM          ( 0),  // w/ weight
_s1pt_1_vsM        ( 0),
_n1Nw_1_vsM        ( 0), // w/o weight
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
vsPtVsPt("NA")
{
  printf("2nd constructor called ");
  
  DefineOutput(0, TList::Class());
  
  printf("passed  ");
  
}

AliAnalysisTaskDptDptCorrelations::~AliAnalysisTaskDptDptCorrelations()
{
  /*
   delete _id_1;
   delete _charge_1; 
   delete _iEtaPhi_1; 
   delete _iPt_1;  
   delete _pt_1;    
   delete _px_1;   
   delete _py_1;   
   delete _pz_1;   
   delete _correction_1;
   delete _dedx_1;
   delete __n1_1_vsPt;
   delete __n1_1_vsEtaPhi;
   delete __s1pt_1_vsEtaPhi;
   delete __n1_1_vsZEtaPhiPt; 
   if (_correctionWeight_1) delete _correctionWeight_1;
   
   if (!_sameFilter)
   {
   delete _id_2;
   delete _charge_2; 
   delete _iEtaPhi_2; 
   delete _iPt_2;  
   delete _pt_2;    
   delete _px_2;   
   delete _py_2;   
   delete _pz_2;   
   delete _correction_2;
   delete _dedx_2;
   delete __n1_2_vsPt;
   delete __n1_2_vsEtaPhi;
   delete __s1pt_2_vsEtaPhi;
   delete __n1_2_vsZEtaPhiPt; 
   if (_correctionWeight_2) delete _correctionWeight_2;
   }
   
   if (!_singlesOnly)
   {
   delete __n2_12_vsPtPt;
   delete __n2_12_vsEtaPhi;
   delete __s2ptpt_12_vsEtaPhi;
   delete __s2PtN_12_vsEtaPhi;
   delete __s2NPt_12_vsEtaPhi;
   }
   */
}

void AliAnalysisTaskDptDptCorrelations::UserCreateOutputObjects()
{
  cout<< "AliAnalysisTaskDptDptCorrelations::CreateOutputObjects() Starting " << endl;
  OpenFile(0);
  _outputHistoList = new TList();
  _outputHistoList->SetOwner();
  
  //if (_useWeights) DefineInput(2, TList::Class());
  
  //Setup the parameters of the histograms
  _nBins_M0 = 500; _min_M0   = 0.;    _max_M0    = 5000.;  _width_M0 = (_max_M0-_min_M0)/_nBins_M0;
  _nBins_M1 = 500; _min_M1   = 0.;    _max_M1    = 5000.;  _width_M1 = (_max_M1-_min_M1)/_nBins_M1;
  _nBins_M2 = 500; _min_M2   = 0.;    _max_M2    = 5000.;  _width_M2 = (_max_M2-_min_M2)/_nBins_M2;
  _nBins_M3 = 500; _min_M3   = 0.;    _max_M3    = 5000.;  _width_M3 = (_max_M3-_min_M3)/_nBins_M3;
  _nBins_M4 = 100; _min_M4   = 0.;    _max_M4    = 100.;   _width_M4 = (_max_M4-_min_M4)/_nBins_M4;
  _nBins_M5 = 100; _min_M5   = 0.;    _max_M5    = 100.;   _width_M5 = (_max_M5-_min_M5)/_nBins_M5;
  _nBins_M6 = 100; _min_M6   = 0.;    _max_M6    = 100.;   _width_M6 = (_max_M6-_min_M6)/_nBins_M6;
  
  _min_vertexZ       = _vertexZMin;  
  _max_vertexZ       = _vertexZMax;  
  _width_vertexZ     = 0.5;
  _nBins_vertexZ     = int(0.5+ (_max_vertexZ - _min_vertexZ)/_width_vertexZ); 
  _nBins_pt_1        = int(0.5+ (_max_pt_1 -_min_pt_1 )/_width_pt_1); 
  _nBins_eta_1       = int(0.5+ (_max_eta_1-_min_eta_1)/_width_eta_1);  
  _width_phi_1       = (_max_phi_1  - _min_phi_1)  /_nBins_phi_1;
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
  
  //setup the work arrays
  
  _id_1       = new int[arraySize];   
  _charge_1   = new int[arraySize]; 
  //_iPhi_1 = new int[arraySize]; 
  //_iEta_1 = new int[arraySize]; 
  _iEtaPhi_1  = new int[arraySize]; 
  _iPt_1      = new int[arraySize];  
  _pt_1       = new float[arraySize];    
  _px_1       = new float[arraySize];   
  _py_1       = new float[arraySize];   
  _pz_1       = new float[arraySize];   
  //_phi_1 = new float[arraySize];  
  //_eta_1 = new float[arraySize];  
  _correction_1 = new float[arraySize];
  _dedx_1     = new float[arraySize];
  
  __n1_1_vsPt              = getDoubleArray(_nBins_pt_1,        0.);
  __n1_1_vsEtaPhi          = getDoubleArray(_nBins_etaPhi_1,    0.);
  __s1pt_1_vsEtaPhi        = getDoubleArray(_nBins_etaPhi_1,    0.);
  __n1_1_vsZEtaPhiPt       = getFloatArray(_nBins_zEtaPhiPt_1,  0.);
  
  cout << "==========================================================================================" << endl;
  cout << "=============== Booking for particle 1 done." << endl;
  cout << "_requestedCharge_1: " << _requestedCharge_1 << endl;
  cout << "_requestedCharge_2: " << _requestedCharge_2 << endl;
  
  if (_requestedCharge_2!=_requestedCharge_1)
    {
    cout << " creating arrays for particle 2 with size: " << arraySize << endl;
    _sameFilter = 0;
    //particle 2
    _id_2       = new int[arraySize];   
    _charge_2   = new int[arraySize]; 
    //_iPhi_2   = new int[arraySize]; 
    //_iEta_2   = new int[arraySize]; 
    _iEtaPhi_2  = new int[arraySize]; 
    _iPt_2      = new int[arraySize];  
    _pt_2       = new float[arraySize];   
    _px_2       = new float[arraySize];   
    _py_2       = new float[arraySize];   
    _pz_2       = new float[arraySize];   
    //_phi_2    = new float[arraySize];  
    //_eta_2    = new float[arraySize];  
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
      AliError("AliAnalysisTaskDptDptCorrelations:: _weight_1 is a null pointer.");
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
        AliError("AliAnalysisTaskDptDptCorrelations:: _weight_1 is a null pointer.");
        return;
        }
      }
    }
  
  createHistograms();
  PostData(0,_outputHistoList);
  
  cout<< "AliAnalysisTaskDptDptCorrelations::CreateOutputObjects() DONE " << endl;
  
}

void  AliAnalysisTaskDptDptCorrelations::createHistograms()
{
  AliInfo(" AliAnalysisTaskDptDptCorrelations::createHistoHistograms() Creating Event Histos");
  TString name;
  
  name = "eventAccounting";
  
  // bin index : what it is...
  //         0 :  number of event submitted
  //         1 :  number accepted by centrality cut
  //         2 :  number accepted by centrality cut and z cut
  //         3 :  total number of particles that satisfied filter 1
  //         4 :  total number of particles that satisfied filter 2
  _eventAccounting      = createHisto1D(name,name,10, -0.5, 9.5, "event Code", _title_counts);
  
  name = "m0"; _m0      = createHisto1D(name,name,_nBins_M1, _min_M1, _max_M1, _title_m0, _title_counts);
  name = "m1"; _m1      = createHisto1D(name,name,_nBins_M1, _min_M1, _max_M1, _title_m1, _title_counts);
  name = "m2"; _m2      = createHisto1D(name,name,_nBins_M2, _min_M2, _max_M2, _title_m2, _title_counts);
  name = "m3"; _m3      = createHisto1D(name,name,_nBins_M3, _min_M3, _max_M3, _title_m3, _title_counts);
  name = "m4"; _m4      = createHisto1D(name,name,_nBins_M4, _min_M4, _max_M4, _title_m4, _title_counts);
  name = "m5"; _m5      = createHisto1D(name,name,_nBins_M5, _min_M5, _max_M5, _title_m5, _title_counts);
  name = "m6"; _m6      = createHisto1D(name,name,_nBins_M6, _min_M6, _max_M6, _title_m6, _title_counts);
  name = "zV"; _vertexZ = createHisto1D(name,name,_nBins_vertexZ, _min_vertexZ, _max_vertexZ, "z-Vertex (cm)", _title_counts);
  

  name = "Eta";     _etadis   = createHisto1F(name,name, 200, -1.1, 1.1, "#eta","counts");
  name = "Phi";     _phidis   = createHisto1F(name,name, 360, 0.0, 6.4, "#phi","counts");
  name = "DCAz";    _dcaz     = createHisto1F(name,name, 500, -5.0, 5.0, "dcaZ","counts");
  name = "DCAxy";   _dcaxy    = createHisto1F(name,name, 500, -5.0, 5.0, "dcaXY","counts");


  if (_singlesOnly)
    {
    name = n1Name+part_1_Name+vsPt;              _n1_1_vsPt              = createHisto1F(name,name, _nBins_pt_1,  _min_pt_1,  _max_pt_1,   _title_pt_1,  _title_AvgN_1);
    name = n1Name+part_1_Name+vsZ+vsEtaPhi+vsPt; _n1_1_vsZVsEtaVsPhiVsPt = createHisto3F(name,name, _nBins_vertexZ,_min_vertexZ,_max_vertexZ, _nBins_etaPhi_1, 0., double(_nBins_etaPhi_1), _nBins_pt_1, _min_pt_1, _max_pt_1, "zVertex", _title_etaPhi_1,  _title_pt_1);
    //name = "dedxVsP_1";                          _dedxVsP_1              = createHisto2F(name,name,400,-2.,2.,120,0.,120.,"p (GeV/c)", "dedx", "counts");
    //name = "corrDedxVsP_1";                      _corrDedxVsP_1          = createHisto2F(name,name,400,-2.,2.,120,0.,120.,"p (GeV/c)", "dedx", "counts");
    //name = "betaVsP_1";                          _betaVsP_1              = createHisto2F(name,name,400,-2.,2.,120,0.5,1.1,"p (GeV/c)", "beta", "counts");

    name = n1Name+part_2_Name+vsPt;              _n1_2_vsPt              = createHisto1F(name,name, _nBins_pt_2,  _min_pt_2,  _max_pt_2,   _title_pt_2,  _title_AvgN_2);
    name = n1Name+part_2_Name+vsZ+vsEtaPhi+vsPt; _n1_2_vsZVsEtaVsPhiVsPt = createHisto3F(name,name, _nBins_vertexZ,_min_vertexZ,_max_vertexZ, _nBins_etaPhi_2, 0., double(_nBins_etaPhi_2), _nBins_pt_2, _min_pt_2, _max_pt_2, "zVertex", _title_etaPhi_2,  _title_pt_2);
    //name = "dedxVsP_2";                          _dedxVsP_2              = createHisto2F(name,name,400,-2.,2.,120,0.,120.,"p (GeV/c)", "dedx", "counts");
    //name = "corrDedxVsP_2";                      _corrDedxVsP_2          = createHisto2F(name,name,400,-2.,2.,120,0.,120.,"p (GeV/c)", "dedx", "counts");
    //name = "betaVsP_2";                          _betaVsP_2              = createHisto2F(name,name,400,-2.,2.,120,0.5,1.1,"p (GeV/c)", "beta", "counts");

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
    
    name = "mInv";     _invMass     = createHisto1F(name,name, 500, 0., 1.000, "M_{inv}","counts"); 
    name = "mInvElec"; _invMassElec = createHisto1F(name,name, 500, 0., 1.000, "M_{inv}","counts"); 
    }
  
  AliInfo(" AliAnalysisTaskDptDptCorrelations::createHistoHistograms() All Done"); 
}
//-----------------------//

void  AliAnalysisTaskDptDptCorrelations::finalizeHistograms()
{
  
  AliInfo("AliAnalysisTaskDptDptCorrelations::finalizeHistograms() starting");
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
  AliInfo("AliAnalysisTaskDptDptCorrelations::finalizeHistograms()  Done ");
}
//--------------//


void  AliAnalysisTaskDptDptCorrelations::UserExec(Option_t */*option*/)
{
  
  int    k1,k2;
  int    iPhi, iEta, iEtaPhi, iPt, charge;
  float  q, p, phi, pt, eta, corr, corrPt, dedx, px, py, pz;
  int    ij;
  int    id_1, q_1, iEtaPhi_1, iPt_1;
  float  pt_1, px_1, py_1, pz_1, corr_1, dedx_1;
  int    id_2, q_2, iEtaPhi_2, iPt_2;
  float  pt_2, px_2, py_2, pz_2, corr_2, dedx_2;
  float  ptpt;
  int    iVertex, iVertexP1, iVertexP2;
  int    iZEtaPhiPt;
  float  massElecSq = 2.5e-7;
  double b[2];
  double bCov[3];
  const  AliAODVertex*	vertex;
  int    nClus;
  bool   bitOK;
    
  // count all events looked at here
  _eventCount++;

  // A. Adare - Fix from above to avoid coverity complaint about null dereference 
  if (_eventAccounting)
    {
      _eventAccounting->Fill(0);// count all calls to this function
    }
  else 
    {
    cout << "AliAnalysisTaskDptDptCorrelations::UserExec(Option_t *option) - !_eventAccounting" << endl;
    cout << "AliAnalysisTaskDptDptCorrelations::UserExec(Option_t *option) - Skip this task" << endl;
    return;
    }
    
  //cout << "AliAnalysisTaskDptDptCorrelations::UserExec(Option_t *option) - 1" << endl;
  
  //Get pointer to current event
  //fESDEvent = dynamic_cast<AliESDEvent*> (InputEvent());
  fAODEvent = dynamic_cast<AliAODEvent*> (InputEvent());
  
  //  if(!fAODEvent && !fESDEvent)
  if(!fAODEvent)
    {
    AliError("AliAnalysisTaskDptDptCorrelations::Exec(Option_t *option) !fAODEvent");
    return;
    }
  
  _eventAccounting->Fill(1);// count all calls to this function with a valid pointer
  //reset single particle counters
  k1 = k2 = 0;
  __n1_1 = __n1_2 = __s1pt_1 = __s1pt_2 = __n1Nw_1 = __n1Nw_2 = __s1ptNw_1 = __s1ptNw_2 = 0;
  
  float v0Centr  = -999.;
  float trkCentr = -999.;
  float spdCentr = -999.;
  
  float vertexX  = -999;
  float vertexY  = -999;
  float vertexZ  = -999;
  float vertexXY = -999;
  //float dcaX     = -999;
  //float dcaY     = -999;
  float dcaZ     = -999;
  float dcaXY    = -999;
  float centrality = -999;
  
  if(fAODEvent)
    {
    //cout << "AliAnalysisTaskDptDptCorrelations::UserExec(Option_t *option) - 5" << endl;
    
    //Centrality
    AliCentrality* centralityObject =  fAODEvent->GetHeader()->GetCentralityP();
    if (centralityObject)
      {
      //cout << "AliAnalysisTaskDptDptCorrelations::UserExec(Option_t *option) - 6" << endl;
      
      v0Centr  = centralityObject->GetCentralityPercentile("V0M");
      trkCentr = centralityObject->GetCentralityPercentile("TRK"); 
      spdCentr = centralityObject->GetCentralityPercentile("CL1");
      //cout << "AliAnalysisTaskDptDptCorrelations::UserExec(Option_t *option) - 7" << endl;
      
      }
    //cout << "AliAnalysisTaskDptDptCorrelations::UserExec(Option_t *option) - 8" << endl;
    
    _nTracks  = fAODEvent->GetNTracks();
    //_mult0    = 0;
    //_mult1    = 0;
    //_mult2    = 0;
    _mult3    = _nTracks; 
    _mult4    = v0Centr;
    _mult5    = trkCentr;
    _mult6    = spdCentr;
    _field    = fAODEvent->GetMagneticField(); 
    
    //cout << "AliAnalysisTaskDptDptCorrelations::UserExec(Option_t *option) - _field:" << _field << endl;
    
    //_centralityMethod
    switch (_centralityMethod)
      {
        case 0: centrality = _mult0; break;
        case 1: centrality = _mult1; break;
        case 2: centrality = _mult2; break;
        case 3: centrality = _mult3; break;
        case 4: centrality = _mult4; break;
        case 5: centrality = _mult5; break;
        case 6: centrality = _mult6; break;
      }
    
    //cout << "AliAnalysisTaskDptDptCorrelations::UserExec(Option_t *option) - 10" << endl;
    
    //filter on centrality
    // require the v0 and trk centralities to agree within 5%.
    
    //if ((Float_t(GetTPCMult(fAOD)) > (-40.3+1.22*GetGlobalMult(fAOD))) && (Float_t(GetTPCMult(fAOD)) < (32.1+1.59*GetGlobalMult(fAOD)))) {//3 sigma

    
    if ( centrality < _centralityMin ||  
         centrality > _centralityMax ||
         fabs(v0Centr-trkCentr)>5.0)
      {
      //cout << "AliAnalysisTaskDptDptCorrelations::UserExec(Option_t *option) - 11" << endl;
      
      // we dont analyze this event here... 
      return;
      }
    _eventAccounting->Fill(2);// count all events with right centrality
    //cout << "AliAnalysisTaskDptDptCorrelations::UserExec(Option_t *option) - 12" << endl;
    
    // filter on z and xy vertex
    vertex = (AliAODVertex*) fAODEvent->GetPrimaryVertexSPD();
    if (!vertex || vertex->GetNContributors()<1)
      {
      vertexZ  = -999;
      vertexXY = -999;
      cout << "AliAnalysisTaskDptDptCorrelations::UserExec(Option_t *option) - No valid vertex object or poor vertex" << endl;
      }
    else
      { 
        vertexX = vertex->GetX();
        vertexY = vertex->GetY();
        vertexZ = vertex->GetZ();
        vertexXY = sqrt(vertexX*vertexX+vertexY*vertexY);
      }
    if (!vertex ||
        vertexZ  < _vertexZMin  || 
        vertexZ  > _vertexZMax  ||
        vertexXY < _vertexXYMin || 
        vertexXY > _vertexXYMax)
      return;
    
    iVertex = int((vertexZ-_min_vertexZ)/_width_vertexZ);
    iVertexP1 = iVertex*_nBins_etaPhiPt_1;
    iVertexP2 = iVertex*_nBins_etaPhiPt_2;
    if (iVertex<0 || iVertex>=_nBins_vertexZ)
      {
      AliError("AliAnalysisTaskDptDptCorrelations::Exec(Option_t *option) iVertex<0 || iVertex>=_nBins_vertexZ ");
      return;
      }
    _eventAccounting->Fill(3);// count all calls to this function with a valid pointer
    
    // loop over all particles for singles
    //debug() << "_nTracks:"<< _nTracks << endl;
    for (int iTrack=0; iTrack< _nTracks; iTrack++)
      {
      //cout << "AliAnalysisTaskDptDptCorrelations::UserExec(Option_t *option) - 16 iTrack:" << iTrack << endl;
      
      AliAODTrack * t = (AliAODTrack *) fAODEvent->GetTrack(iTrack);
      if (!t) 
        {
        AliError(Form("AliAnalysisTaskDptDptCorrelations::Exec(Option_t *option) No track ofr iTrack=%d", iTrack));
        continue;
        }
      
      q      = t->Charge();
      charge = int(q);
      phi    = t->Phi();
      p      = t->P();
      pt     = t->Pt();
      px     = t->Px();
      py     = t->Py();
      pz     = t->Pz();
      eta    = t->Eta();
      dedx   = t->GetTPCsignal();
      nClus  = t->GetTPCNcls();
      bitOK  = t->TestFilterBit(_trackFilterBit);
      
      //cout << "_trackFilterBit:" << _trackFilterBit << " Track returns:" << bitOK << endl;
      //cout << "    q:" << q    << " _requestedCharge_1:" << _requestedCharge_1 << endl;
      //cout << "   pt:" << pt   << " _min_pt_1:" << _min_pt_1 << " _max_pt_1:" << _max_pt_1<< endl;
      //cout << "  phi:" << phi  << endl;
      //cout << "  eta:" << eta  << " _min_eta_1:" << _min_eta_1 << " _max_eta_1:" << _max_eta_1<< endl;
      //cout << " dedx:" << dedx << " _dedxMin:" << _dedxMin << " _dedxMax:" << _dedxMax << endl;
      //cout << "nclus:" << nClus<< " _nClusterMin:" << _nClusterMin << endl;
      
      //common cuts first
      
      if (!bitOK || dedx<_dedxMin || dedx>_dedxMax || nClus<_nClusterMin) continue;
      //Particle 1
      if (_requestedCharge_1 == charge &&
          pt       >=  _min_pt_1 && 
          pt       <   _max_pt_1 && 
          eta      >=  _min_eta_1 && 
          eta      <   _max_eta_1) 
        {
        
	  /*
        if (t->PropagateToDCA(vertex, _field, 100., b, bCov) )
          {
          dcaXY = b[0];
          dcaZ  = b[1];
          } 
        else
          {
          dcaXY = -999999;
          dcaZ  = -999999;
          }
	  */
	  
	  dcaXY = t->DCA(); //new change Prabhat                                        
	  dcaZ  = t->ZAtDCA(); //new change Prabhat  
	  
	   // skip track if DCA too large
	     if (dcaZ  <  _dcaZMin || 
	     dcaZ     >   _dcaZMax || 
	     dcaXY    <  _dcaXYMin || 
	     dcaXY    >   _dcaXYMax)
	     continue; //track does not have a valid DCA
	  

	  //if (dcaZ > _dcaZMax) continue;
	     //cout <<"Prabhat=============="<<"  "<<_dcaZMax<<"    "<<dcaZ<<endl;
	_etadis->Fill(eta); //prabhat QA                                  
	_phidis->Fill(phi);
	_dcaz->Fill(dcaZ); //Prabhat QA             
        _dcaxy->Fill(dcaXY);


        iPhi   = int( phi/_width_phi_1);
        
        if (iPhi<0 || iPhi>=_nBins_phi_1 ) 
          {
          AliWarning("AliAnalysisTaskDptDptCorrelations::analyze() iPhi<0 || iPhi>=_nBins_phi_1");
          return;
          }
        
        iEta    = int((eta-_min_eta_1)/_width_eta_1); 
        if (iEta<0 || iEta>=_nBins_eta_1) 
          {
          AliWarning(Form("AliAnalysisTaskDptDptCorrelations::analyze(AliceEvent * event) Mismatched iEta: %d", iEta));
          continue;
          }
        iPt     = int((pt -_min_pt_1 )/_width_pt_1 ); 
        if (iPt<0  || iPt >=_nBins_pt_1)
          {
          AliWarning(Form("AliAnalysisTaskDptDptCorrelations::analyze(AliceEvent * event) Mismatched iPt: %d",iPt));
          continue;
          }
        iEtaPhi = iEta*_nBins_phi_1+iPhi;
        iZEtaPhiPt = iVertexP1 + iEtaPhi*_nBins_pt_1 + iPt;
        
        if (_correctionWeight_1)
          corr = _correctionWeight_1[iZEtaPhiPt];
        else
          corr = 1;
        if (iZEtaPhiPt<0 || iZEtaPhiPt>=_nBins_zEtaPhiPt_1)
          {
          AliWarning("AliAnalysisTaskDptDptCorrelations::analyze(AliceEvent * event) iZEtaPhiPt<0 || iZEtaPhiPt>=_nBins_zEtaPhiPt_1");
          continue;
          }
        

        if (_singlesOnly)
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
            AliError(Form("AliAnalysisTaskDptDptCorrelations::analyze(AliceEvent * event) k1 >=arraySize; arraySize: %d",arraySize));
            return;
            }
          }
        }
      
      if (!_sameFilter &&
          _requestedCharge_2 == charge &&
          pt       >=  _min_pt_2 && 
          pt       <   _max_pt_2 && 
          eta      >=  _min_eta_2 && 
          eta      <   _max_eta_2)  
        {
	  /* 
	     if (t->PropagateToDCA(vertex, _field, 100., b, bCov) )
	     {
	     dcaXY = b[0];
	     dcaZ  = b[1];
	     } 
	     else
	     {
	     dcaXY = -999999;
	     dcaZ  = -999999;
	     }
	  */
	  
	  dcaXY = t->DCA(); //new change Prabhat                                        
	  dcaZ  = t->ZAtDCA(); //new change Prabhat  
	  //if (dcaZ > _dcaZMax) continue;
	  //cout <<"Prabhat=============="<<"  "<<_dcaZMax<<"    "<<dcaZ<<endl;

	  
	  if (dcaZ     <  _dcaZMin || 
	      dcaZ     >   _dcaZMax || 
	      dcaXY    <  _dcaXYMin || 
	      dcaXY    >   _dcaXYMax)
	    continue; //track does not have a valid DCA
	  

	  iPhi   = int( phi/_width_phi_2);
	  
        if (iPhi<0 || iPhi>=_nBins_phi_2 ) 
          {
          AliWarning("AliAnalysisTaskDptDptCorrelations::analyze() iPhi<0 || iPhi>=_nBins_phi_1");
          return;
          }
        
        iEta    = int((eta-_min_eta_2)/_width_eta_2);
        if (iEta<0 || iEta>=_nBins_eta_2) 
          {
          AliWarning(Form("AliAnalysisTaskDptDptCorrelations::analyze(AliceEvent * event) Mismatched iEta: %d", iEta));
          continue;
          }
        iPt     = int((pt -_min_pt_2 )/_width_pt_2 ); 
        if (iPt<0  || iPt >=_nBins_pt_2)
          {
          AliWarning(Form("AliAnalysisTaskDptDptCorrelations::analyze(AliceEvent * event) Mismatched iPt: %d",iPt));
          continue;
          }
        
        iEtaPhi = iEta*_nBins_phi_2+iPhi;
        iZEtaPhiPt = iVertexP2 + iEtaPhi*_nBins_pt_2 + iPt;
        if (iZEtaPhiPt<0 || iZEtaPhiPt>=_nBins_zEtaPhiPt_2)
          {
          AliWarning("AliAnalysisTaskDptDptCorrelations::analyze(AliceEvent * event) iZEtaPhiPt<0 || iZEtaPhiPt>=_nBins_zEtaPhiPt_2");
          continue;
          }
        

        if (_correctionWeight_2)
          corr = _correctionWeight_2[iZEtaPhiPt];
        else
          corr = 1;
        //dpt = pt - (charge>0) ? _avgPt_vsEtaPhi_2p[iEtaPhi] : _avgPt_vsEtaPhi_2m[iEtaPhi];       

        if (_singlesOnly)
          {
          //_dedxVsP_2->Fill(p*q,dedx);
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
      
      //cout << "done with track" << endl;
      } //iTrack
    } //aod 
    }

  
  //cout << "Filling histograms now" << endl;
  _m0->Fill(_mult0);
  _m1->Fill(_mult1);
  _m2->Fill(_mult2);
  _m3->Fill(_mult3);
  _m4->Fill(_mult4);
  _m5->Fill(_mult5);
  _m6->Fill(_mult6);
  _vertexZ->Fill(vertexZ);
      
  if (_singlesOnly)
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
              if (_rejectPairConversion)
                {
                float e1Sq = massElecSq + pt_1*pt_1 + pz_1*pz_1;
                float e2Sq = massElecSq + pt_2*pt_2 + pz_2*pz_2;
                float mInvSq = 2*(massElecSq + sqrt(e1Sq*e2Sq) - px_1*px_2 - py_1*py_2 - pz_1*pz_2 );
                float mInv = sqrt(mInvSq);
                _invMass->Fill(mInv);
                if (mInv<0.5)
                  {
                  if (dedx_1>75. && dedx_2>75.)
                    {
                    _invMassElec->Fill(mInv);
                    if (mInv<0.05) continue;
                    }
                  }
                }
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
  
  
  cout << "AliAnalysisTaskDptDptCorrelations::UserExec()   -----------------Event Done " << endl;
  PostData(0,_outputHistoList);
  
}

void   AliAnalysisTaskDptDptCorrelations::FinishTaskOutput()
{
  cout << "AliAnalysisTaskDptDptCorrelations::FinishTaskOutput() Starting." << endl;
  cout << "= 0 ====================================================================" << endl;
  finalizeHistograms();
  cout << "= 1 ====================================================================" << endl;
  PostData(0,_outputHistoList);
  cout << "= 2 ====================================================================" << endl;
  cout << "AliAnalysisTaskDptDptCorrelations::FinishTaskOutput() Done." << endl;
}

void   AliAnalysisTaskDptDptCorrelations::Terminate(Option_t* /*option*/)
{
  cout << "AliAnalysisTaskDptDptCorrelations::Terminate() Starting/Done." << endl;
}


//Tools
//===================================================================================================
void  AliAnalysisTaskDptDptCorrelations::fillHistoWithArray(TH1 * h, float * array, int size)
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

void  AliAnalysisTaskDptDptCorrelations::fillHistoWithArray(TH2 * h, float * array, int size1, int size2)
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

void  AliAnalysisTaskDptDptCorrelations::fillHistoWithArray(TH3 * h, float * array, int size1, int size2, int size3)
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

void  AliAnalysisTaskDptDptCorrelations::fillHistoWithArray(TH1 * h, double * array, int size)
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

void  AliAnalysisTaskDptDptCorrelations::fillHistoWithArray(TH2 * h, double * array, int size1, int size2)
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

void  AliAnalysisTaskDptDptCorrelations::fillHistoWithArray(TH3 * h, double * array, int size1, int size2, int size3)
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
double *  AliAnalysisTaskDptDptCorrelations::getDoubleArray(int size, double v)
{
  /// Allocate an array of type double with n values
  /// Initialize the array to the given value
  double * array = new double [size];
  for (int i=0;i<size;++i) array[i]=v;
  return array;
}

//________________________________________________________________________
float *  AliAnalysisTaskDptDptCorrelations::getFloatArray(int size, float v)
{
  /// Allocate an array of type float with n values
  /// Initialize the array to the given value
  float * array = new float [size];
  for (int i=0;i<size;++i) array[i]=v;
  return array;
}


//________________________________________________________________________
TH1D * AliAnalysisTaskDptDptCorrelations::createHisto1D(const TString &  name, const TString &  title, 
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
TH1D * AliAnalysisTaskDptDptCorrelations::createHisto1D(const TString &  name, const TString &  title, 
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
TH2D * AliAnalysisTaskDptDptCorrelations::createHisto2D(const TString &  name, const TString &  title, 
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
TH2D * AliAnalysisTaskDptDptCorrelations::createHisto2D(const TString &  name, const TString &  title, 
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
TH1F * AliAnalysisTaskDptDptCorrelations::createHisto1F(const TString &  name, const TString &  title, 
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
TH1F * AliAnalysisTaskDptDptCorrelations::createHisto1F(const TString &  name, const TString &  title, 
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
TH2F * AliAnalysisTaskDptDptCorrelations::createHisto2F(const TString &  name, const TString &  title, 
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
TH2F * AliAnalysisTaskDptDptCorrelations::createHisto2F(const TString &  name, const TString &  title, 
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
TH3F * AliAnalysisTaskDptDptCorrelations::createHisto3F(const TString &  name, const TString &  title, 
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
TProfile * AliAnalysisTaskDptDptCorrelations::createProfile(const TString & name, const TString & description, 
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
TProfile * AliAnalysisTaskDptDptCorrelations::createProfile(const TString &  name,const TString &  description, 
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


void   AliAnalysisTaskDptDptCorrelations::addToList(TH1 *h)
{
  if (_outputHistoList)
    {
    _outputHistoList->Add(h);
    }
  else
    cout << "addToList(TH1 *h) _outputHistoList is null!!!!! Shoudl abort ship" << endl;

}



