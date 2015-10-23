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
#include "AliAnalysisTask3PCorrelations.h"
#include <iostream>

using std::cout;
using std::endl;

ClassImp(AliAnalysisTask3PCorrelations)

AliAnalysisTask3PCorrelations::AliAnalysisTask3PCorrelations()
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
_rejectPileup  ( 1), 
_rejectPairConversion ( 0), 
_vertexZMin           ( -10), 
_vertexZMax           (  10), 
_vertexXYMin          ( -10),
_vertexXYMax          (  10),
_centralityMethod     (  4),
_centralityMin        (  0.),
_centralityMax        (  0.),
_dcaZMin              ( -3), 
_dcaZMax              (  3.), 
_dcaXYMin             ( -3.), 
_dcaXYMax             (  3.),
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
_iPhi_1(0),    
_pt_1(0),       
_px_1(0),      
_py_1(0),      
_pz_1(0),      
_correction_1(0),   
_dedx_1(0),   
_id_2(0),      
_charge_2(0),    
_iPhi_2(0),    
_pt_2(0),      
_px_2(0),      
_py_2(0),      
_pz_2(0),      
_correction_2(0),   
_dedx_2(0),   
_id_3(0),      
_charge_3(0),    
_iPhi_3(0),    
_pt_3(0),      
_px_3(0),      
_py_3(0),      
_pz_3(0),      
_correction_3(0),   
_dedx_3(0),   

_correctionWeight_1P(0),   
_correctionWeight_1M(0),   
_correctionWeight_2P(0),   
_correctionWeight_2M(0),   
_correctionWeight_3P(0),   
_correctionWeight_3M(0),   
_nBins_M0(500),       _min_M0(0),        _max_M0(10000),          _width_M0(20),
_nBins_M1(500),       _min_M1(0),        _max_M1(10000),          _width_M1(20),
_nBins_M2(500),       _min_M2(0),        _max_M2(10000),          _width_M2(20),
_nBins_M3(500),       _min_M3(0),        _max_M3(10000),          _width_M3(20),
_nBins_M4(100),       _min_M4(0),        _max_M4(1),              _width_M4(0.01),
_nBins_M5(100),       _min_M5(0),        _max_M5(1),              _width_M5(0.01),
_nBins_M6(100),       _min_M6(0),        _max_M6(1),              _width_M6(0.01),
_nBins_vsM(0),        _min_vsM(0),       _max_vsM(0),             _width_vsM(0),
_nBins_vertexZ(40),   _min_vertexZ(-10), _max_vertexZ(10),        _width_vertexZ(0.5),

_nBins_pt_1(8),       _min_pt_1(2.),    _max_pt_1(10.0),          _width_pt_1(1.0),
_nBins_phi_1(72),     _min_phi_1(0),     _max_phi_1(2.*3.1415927),_width_phi_1(2.*3.1415927/72.),
_nBins_eta_1(1),      _min_eta_1(-0.25), _max_eta_1(0.25),        _width_eta_1(0.5),
_nBins_etaPhi_1(0), 
_nBins_etaPhiPt_1(0),
_nBins_zEtaPhiPt_1(0),

_nBins_pt_2(5),       _min_pt_2(0.5),     _max_pt_2(1.0),          _width_pt_2(0.1),
_nBins_phi_2(72),     _min_phi_2(0),      _max_phi_2(2.*3.1415927),_width_phi_2(2.*3.1415927/72),
_nBins_eta_2(1),      _min_eta_2(-1),     _max_eta_2(-0.7),        _width_eta_2(0.3),
_nBins_etaPhi_2(0), 
_nBins_etaPhiPt_2(0),
_nBins_zEtaPhiPt_2(0),

_nBins_pt_3(5),       _min_pt_3(0.5),     _max_pt_3(1.0),          _width_pt_3(0.1),
_nBins_phi_3(72),     _min_phi_3(0),      _max_phi_3(2.*3.1415927),_width_phi_3(2.*3.1415927/72),
_nBins_eta_3(1),      _min_eta_3(0.7),    _max_eta_3(1.),          _width_eta_3(0.3),
_nBins_etaPhi_3(0), 
_nBins_etaPhiPt_3(0),
_nBins_zEtaPhiPt_3(0),

_nBins_phi_12(0),
_nBins_phi_13(0),
_nBins_phi_23(0),
_nBins_phi_123(0),


_nBins_etaPhi_12(0),
_nBins_etaPhi_13(0),
_nBins_etaPhi_23(0),
_nBins_etaPhi_123(0),

__n1_1(0),
__n1_2(0),
__n1_3(0),
__n2_12(0),   
__n2_13(0),   
__n2_23(0),   
__n3_123(0),   
__n1Nw_1(0),
__n1Nw_2(0),
__n1Nw_3(0),
__n2Nw_12(0),   
__n2Nw_13(0),   
__n2Nw_23(0),   
__n3Nw_123(0),  

__n1_1_vsPt(0),
__n1_1_vsPhi(0), 
__n1_1P_vsZEtaPhiPt(0),
__n1_1M_vsZEtaPhiPt(0),

__n1_2_vsPt(0),
__n1_2_vsPhi(0), 
__n1_2P_vsZEtaPhiPt(0),
__n1_2M_vsZEtaPhiPt(0),

__n1_3_vsPt(0),
__n1_3_vsPhi(0), 
__n1_3P_vsZEtaPhiPt(0),
__n1_3M_vsZEtaPhiPt(0),

__n2_12_vsPhi(0),
__n2_13_vsPhi(0),
__n2_23_vsPhi(0),
__n3_123_vsPhi(0),

_weight_1P      ( 0    ),
_weight_1M      ( 0    ),
_weight_2P      ( 0    ),
_weight_2M      ( 0    ),
_weight_3P      ( 0    ),
_weight_3M      ( 0    ),
_eventAccounting ( 0),
_m0 ( 0),
_m1 ( 0),
_m2 ( 0),
_m3 ( 0),
_m4 ( 0),
_m5 ( 0),
_m6 ( 0),
_vertexZ        ( 0),

_n1_1_vsPhi     ( 0),
_n1_1_vsM          ( 0), // w/ weight
_n1Nw_1_vsM        ( 0), // w/o weight

_n1_1_vsPt      ( 0),         
_n1_1P_vsZVsEtaVsPhiVsPt ( 0),
_n1_1M_vsZVsEtaVsPhiVsPt ( 0),
_dedxVsP_1         ( 0),
_corrDedxVsP_1     ( 0),
_betaVsP_1         ( 0),

_n1_2_vsPhi        ( 0),
_n1_2_vsM          ( 0),
_n1Nw_2_vsM        ( 0),
_n1_2_vsPt         ( 0),       
_n1_2P_vsZVsEtaVsPhiVsPt ( 0), 
_n1_2M_vsZVsEtaVsPhiVsPt ( 0), 
_dedxVsP_2         ( 0),
_corrDedxVsP_2     ( 0),
_betaVsP_2         ( 0),

_n1_3_vsPhi        ( 0),
_n1_3_vsM          ( 0),
_n1Nw_3_vsM        ( 0),
_n1_3_vsPt         ( 0),       
_n1_3P_vsZVsEtaVsPhiVsPt ( 0), 
_n1_3M_vsZVsEtaVsPhiVsPt ( 0), 
_dedxVsP_3         ( 0),
_corrDedxVsP_3     ( 0),
_betaVsP_3         ( 0),

_n2_12_vsPhi    ( 0),  
_n2_13_vsPhi    ( 0),  
_n2_23_vsPhi    ( 0),  
_n3_123_vsPhi   ( 0),  

_n2_12_vsM         ( 0),        
_n2_13_vsM         ( 0),        
_n2_23_vsM         ( 0),        
_n3_123_vsM        ( 0),        

_n2Nw_12_vsM       ( 0),        
_n2Nw_13_vsM       ( 0),        
_n2Nw_23_vsM       ( 0),        
_n3Nw_123_vsM      ( 0),        

n1Name("NA"),
n2Name("NA"),
n3Name("NA"),
n1NwName("NA"),
n2NwName("NA"),
n3NwName("NA"),
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

intR2Name("NA"),
binCorrName("NA"),
intBinCorrName("NA"),

countsName("NA"),
part_1_Name("NA"),
part_1P_Name("NA"),
part_1M_Name("NA"),
part_2_Name("NA"),
part_2P_Name("NA"),
part_2M_Name("NA"),
part_3_Name("NA"),
part_3P_Name("NA"),
part_3M_Name("NA"),
pair_12_Name("NA"),
pair_13_Name("NA"),
pair_23_Name("NA"),
triplet_Name("NA"),

_title_counts("NA"),
_title_m0("NA"),
_title_m1("NA"),
_title_m2("NA"),
_title_m3("NA"),
_title_m4("NA"),
_title_m5("NA"),
_title_m6("NA"),

_title_phi_1("NA"),
_title_phi_2("NA"),
_title_phi_3("NA"),
_title_phi_12("NA"),
_title_phi_13("NA"),
_title_phi_23("NA"),
_title_phi_123("NA"),

vsZ("NA"),
vsM("NA"),
vsPt("NA"),
vsPhi("NA"), 
vsEta("NA"), 
vsEtaPhi("NA")
{
  printf("Default constructor called \n");
  
  printf("passed \n ");
  
}

AliAnalysisTask3PCorrelations::AliAnalysisTask3PCorrelations(const TString & name)
: AliAnalysisTaskSE(name),
fAODEvent(0), 
fESDEvent(0),             //! ESD Event 
fInputHandler(0),
_outputHistoList(0),
_twoPi         ( 2.0 * 3.1415927),
_eventCount    ( 0), 
_debugLevel    ( 0),
_singlesOnly   ( 0), 
_useWeights    ( 0), 
_rejectPileup  ( 1), 
_rejectPairConversion ( 0), 
_vertexZMin           ( -10), 
_vertexZMax           (  10), 
_vertexXYMin          ( -10),
_vertexXYMax          (  10),
_centralityMethod     (  4),
_centralityMin        (  0.),
_centralityMax        (  0.),
_dcaZMin              ( -3), 
_dcaZMax              (  3.), 
_dcaXYMin             ( -3.), 
_dcaXYMax             (  3.),
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
_iPhi_1(0),    
_pt_1(0),       
_px_1(0),      
_py_1(0),      
_pz_1(0),      
_correction_1(0),   
_dedx_1(0),   
_id_2(0),      
_charge_2(0),    
_iPhi_2(0),    
_pt_2(0),      
_px_2(0),      
_py_2(0),      
_pz_2(0),      
_correction_2(0),   
_dedx_2(0),   
_id_3(0),      
_charge_3(0),    
_iPhi_3(0),    
_pt_3(0),      
_px_3(0),      
_py_3(0),      
_pz_3(0),      
_correction_3(0),   
_dedx_3(0),   

_correctionWeight_1P(0),   
_correctionWeight_1M(0),   
_correctionWeight_2P(0),   
_correctionWeight_2M(0),   
_correctionWeight_3P(0),   
_correctionWeight_3M(0),   
_nBins_M0(500),       _min_M0(0),        _max_M0(10000),          _width_M0(20),
_nBins_M1(500),       _min_M1(0),        _max_M1(10000),          _width_M1(20),
_nBins_M2(500),       _min_M2(0),        _max_M2(10000),          _width_M2(20),
_nBins_M3(500),       _min_M3(0),        _max_M3(10000),          _width_M3(20),
_nBins_M4(100),       _min_M4(0),        _max_M4(1),              _width_M4(0.01),
_nBins_M5(100),       _min_M5(0),        _max_M5(1),              _width_M5(0.01),
_nBins_M6(100),       _min_M6(0),        _max_M6(1),              _width_M6(0.01),
_nBins_vsM(0),        _min_vsM(0),       _max_vsM(0),             _width_vsM(0),

_nBins_vertexZ(40),   _min_vertexZ(-10), _max_vertexZ(10),        _width_vertexZ(0.5),

_nBins_pt_1(8),       _min_pt_1(2.),    _max_pt_1(10.0),          _width_pt_1(1.0),
_nBins_phi_1(72),     _min_phi_1(0),     _max_phi_1(2.*3.1415927),_width_phi_1(2.*3.1415927/72.),
_nBins_eta_1(1),      _min_eta_1(-0.25), _max_eta_1(0.25),        _width_eta_1(0.5),
_nBins_etaPhi_1(0), 
_nBins_etaPhiPt_1(0),
_nBins_zEtaPhiPt_1(0),

_nBins_pt_2(5),       _min_pt_2(0.5),     _max_pt_2(1.0),          _width_pt_2(0.1),
_nBins_phi_2(72),     _min_phi_2(0),      _max_phi_2(2.*3.1415927),_width_phi_2(2.*3.1415927/72),
_nBins_eta_2(1),      _min_eta_2(-1),     _max_eta_2(-0.7),        _width_eta_2(0.3),
_nBins_etaPhi_2(0), 
_nBins_etaPhiPt_2(0),
_nBins_zEtaPhiPt_2(0),

_nBins_pt_3(5),       _min_pt_3(0.5),     _max_pt_3(1.0),          _width_pt_3(0.1),
_nBins_phi_3(72),     _min_phi_3(0),      _max_phi_3(2.*3.1415927),_width_phi_3(2.*3.1415927/72),
_nBins_eta_3(1),      _min_eta_3(0.7),    _max_eta_3(1.),          _width_eta_3(0.3),
_nBins_etaPhi_3(0), 
_nBins_etaPhiPt_3(0),
_nBins_zEtaPhiPt_3(0),

_nBins_phi_12(0),
_nBins_phi_13(0),
_nBins_phi_23(0),
_nBins_phi_123(0),

_nBins_etaPhi_12(0),
_nBins_etaPhi_13(0),
_nBins_etaPhi_23(0),
_nBins_etaPhi_123(0),

__n1_1(0),
__n1_2(0),
__n1_3(0),
__n2_12(0),   
__n2_13(0),   
__n2_23(0),   
__n3_123(0),   
__n1Nw_1(0),
__n1Nw_2(0),
__n1Nw_3(0),
__n2Nw_12(0),   
__n2Nw_13(0),   
__n2Nw_23(0),   
__n3Nw_123(0),  

__n1_1_vsPt(0),
__n1_1_vsPhi(0), 
__n1_1P_vsZEtaPhiPt(0),
__n1_1M_vsZEtaPhiPt(0),

__n1_2_vsPt(0),
__n1_2_vsPhi(0), 
__n1_2P_vsZEtaPhiPt(0),
__n1_2M_vsZEtaPhiPt(0),

__n1_3_vsPt(0),
__n1_3_vsPhi(0), 
__n1_3P_vsZEtaPhiPt(0),
__n1_3M_vsZEtaPhiPt(0),

__n2_12_vsPhi(0),
__n2_13_vsPhi(0),
__n2_23_vsPhi(0),
__n3_123_vsPhi(0),

_weight_1P      ( 0    ),
_weight_1M      ( 0    ),
_weight_2P      ( 0    ),
_weight_2M      ( 0    ),
_weight_3P      ( 0    ),
_weight_3M      ( 0    ),
_eventAccounting ( 0),
_m0 ( 0),
_m1 ( 0),
_m2 ( 0),
_m3 ( 0),
_m4 ( 0),
_m5 ( 0),
_m6 ( 0),
_vertexZ        ( 0),

_n1_1_vsPhi     ( 0),
_n1_1_vsM          ( 0), // w/ weight
_n1Nw_1_vsM        ( 0), // w/o weight

_n1_1_vsPt      ( 0),         
_n1_1P_vsZVsEtaVsPhiVsPt ( 0),
_n1_1M_vsZVsEtaVsPhiVsPt ( 0),
_dedxVsP_1         ( 0),
_corrDedxVsP_1     ( 0),
_betaVsP_1         ( 0),

_n1_2_vsPhi        ( 0),
_n1_2_vsM          ( 0),
_n1Nw_2_vsM        ( 0),
_n1_2_vsPt         ( 0),       
_n1_2P_vsZVsEtaVsPhiVsPt ( 0), 
_n1_2M_vsZVsEtaVsPhiVsPt ( 0), 
_dedxVsP_2         ( 0),
_corrDedxVsP_2     ( 0),
_betaVsP_2         ( 0),

_n1_3_vsPhi        ( 0),
_n1_3_vsM          ( 0),
_n1Nw_3_vsM        ( 0),
_n1_3_vsPt         ( 0),       
_n1_3P_vsZVsEtaVsPhiVsPt ( 0), 
_n1_3M_vsZVsEtaVsPhiVsPt ( 0), 
_dedxVsP_3         ( 0),
_corrDedxVsP_3     ( 0),
_betaVsP_3         ( 0),

_n2_12_vsPhi    ( 0),  
_n2_13_vsPhi    ( 0),  
_n2_23_vsPhi    ( 0),  
_n3_123_vsPhi   ( 0),  

_n2_12_vsM         ( 0),        
_n2_13_vsM         ( 0),        
_n2_23_vsM         ( 0),        
_n3_123_vsM        ( 0),        

_n2Nw_12_vsM       ( 0),        
_n2Nw_13_vsM       ( 0),        
_n2Nw_23_vsM       ( 0),        
_n3Nw_123_vsM      ( 0),        

n1Name("NA"),
n2Name("NA"),
n3Name("NA"),
n1NwName("NA"),
n2NwName("NA"),
n3NwName("NA"),
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

intR2Name("NA"),
binCorrName("NA"),
intBinCorrName("NA"),

countsName("NA"),
part_1_Name("NA"),
part_1P_Name("NA"),
part_1M_Name("NA"),
part_2_Name("NA"),
part_2P_Name("NA"),
part_2M_Name("NA"),
part_3_Name("NA"),
part_3P_Name("NA"),
part_3M_Name("NA"),
pair_12_Name("NA"),
pair_13_Name("NA"),
pair_23_Name("NA"),
triplet_Name("NA"),

_title_counts("NA"),
_title_m0("NA"),
_title_m1("NA"),
_title_m2("NA"),
_title_m3("NA"),
_title_m4("NA"),
_title_m5("NA"),
_title_m6("NA"),

_title_phi_1("NA"),
_title_phi_2("NA"),
_title_phi_3("NA"),

_title_phi_12("NA"),
_title_phi_13("NA"),
_title_phi_23("NA"),
_title_phi_123("NA"),

vsZ("NA"),
vsM("NA"),
vsPt("NA"),
vsPhi("NA"), 
vsEta("NA"), 
vsEtaPhi("NA")
{
  printf("2nd constructor called ");
  
  DefineOutput(0, TList::Class());
  
  printf("passed  ");
  
}

AliAnalysisTask3PCorrelations::~AliAnalysisTask3PCorrelations()
{
  // no ops --- at least for now....
}

void AliAnalysisTask3PCorrelations::UserCreateOutputObjects()
{
  cout<< "AliAnalysisTask3PCorrelations::CreateOutputObjects() Starting " << endl;
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
  _nBins_zEtaPhiPt_1 = _nBins_vertexZ  * _nBins_phi_1 * _nBins_eta_1 * _nBins_pt_1;
  
  _nBins_pt_2        =  int(0.5+ (_max_pt_2 -_min_pt_2 )/_width_pt_2);  
  _nBins_eta_2       = int(0.5+ (_max_eta_2-_min_eta_2)/_width_eta_2); 
  _width_phi_2       = (_max_phi_2  - _min_phi_2)  /_nBins_phi_2;
  _nBins_etaPhi_2    = _nBins_phi_2    * _nBins_eta_2;
  _nBins_zEtaPhiPt_2 = _nBins_vertexZ   * _nBins_phi_2 * _nBins_eta_2 * _nBins_pt_2;
  
  _nBins_pt_3        =  int(0.5+ (_max_pt_3 -_min_pt_3 )/_width_pt_3);  
  _nBins_eta_3       = int(0.5+ (_max_eta_3-_min_eta_3)/_width_eta_3); 
  _width_phi_3       = (_max_phi_3  - _min_phi_3)  /_nBins_phi_3;
  _nBins_etaPhi_3    = _nBins_phi_3    * _nBins_eta_3;
  _nBins_zEtaPhiPt_3 = _nBins_vertexZ   * _nBins_phi_3 * _nBins_eta_3 * _nBins_pt_3;
  
  
  _nBins_etaPhi_12   = _nBins_etaPhi_1 * _nBins_etaPhi_2;
  _nBins_etaPhi_13   = _nBins_etaPhi_1 * _nBins_etaPhi_3;
  _nBins_etaPhi_23   = _nBins_etaPhi_2 * _nBins_etaPhi_3;
  _nBins_etaPhi_123  = _nBins_etaPhi_1 * _nBins_etaPhi_2 * _nBins_etaPhi_3;
  
  _nBins_phi_12   = _nBins_phi_1 * _nBins_phi_2;
  _nBins_phi_13   = _nBins_phi_1 * _nBins_phi_3;
  _nBins_phi_23   = _nBins_phi_2 * _nBins_phi_3;
  _nBins_phi_123  = _nBins_phi_1 * _nBins_phi_2 * _nBins_phi_3;
  
  
  //setup the work arrays
  if (_singlesOnly)
    {
    __n1_1_vsPt         = getDoubleArray(_nBins_pt_1,         0.);
    __n1_1P_vsZEtaPhiPt = getFloatArray (_nBins_zEtaPhiPt_1,  0.);
    __n1_1M_vsZEtaPhiPt = getFloatArray (_nBins_zEtaPhiPt_1,  0.);
    __n1_2_vsPt         = getDoubleArray(_nBins_pt_2,         0.);
    __n1_2P_vsZEtaPhiPt = getFloatArray (_nBins_zEtaPhiPt_2,  0.);
    __n1_2M_vsZEtaPhiPt = getFloatArray (_nBins_zEtaPhiPt_2,  0.);
    __n1_3_vsPt         = getDoubleArray(_nBins_pt_3,         0.);
    __n1_3P_vsZEtaPhiPt = getFloatArray (_nBins_zEtaPhiPt_3,  0.);
    __n1_3M_vsZEtaPhiPt = getFloatArray (_nBins_zEtaPhiPt_3,  0.);
    }
  else
    {
    _id_1          = new int[arraySize];   
    //_charge_1      = new int[arraySize]; 
    _iPhi_1        = new int[arraySize]; 
    //_pt_1          = new float[arraySize];    
    //_px_1          = new float[arraySize];   
    //_py_1          = new float[arraySize];   
    //_pz_1          = new float[arraySize];   
    _correction_1  = new float[arraySize];
    //_dedx_1        = new float[arraySize];
    __n1_1_vsPhi   = getDoubleArray(_nBins_phi_1,    0.);
    
    _id_2          = new int[arraySize];   
    //_charge_2      = new int[arraySize]; 
    _iPhi_2        = new int[arraySize]; 
    //_pt_2          = new float[arraySize];   
    //_px_2          = new float[arraySize];   
    //_py_2          = new float[arraySize];   
    //_pz_2          = new float[arraySize];   
    _correction_2  = new float[arraySize];
    //_dedx_2        = new float[arraySize];
    __n1_2_vsPhi   = getDoubleArray(_nBins_phi_2,       0.);
    
    _id_3          = new int[arraySize];   
    //_charge_3      = new int[arraySize]; 
    _iPhi_3        = new int[arraySize]; 
    //_pt_3          = new float[arraySize];   
    //_px_3          = new float[arraySize];   
    //_py_3          = new float[arraySize];   
    //_pz_3          = new float[arraySize];   
    _correction_3  = new float[arraySize];
    //_dedx_3        = new float[arraySize];
    __n1_3_vsPhi        = getDoubleArray(_nBins_phi_3,       0.);
    
    __n2_12_vsPhi    = getFloatArray(_nBins_phi_12,    0.);
    __n2_13_vsPhi    = getFloatArray(_nBins_phi_13,    0.);
    __n2_23_vsPhi    = getFloatArray(_nBins_phi_23,    0.);
    __n3_123_vsPhi   = getFloatArray(_nBins_phi_123,   0.);
    }
  
  
  // Setup all the labels needed.
  part_1_Name   = "_1";
  part_1P_Name  = "_1P";
  part_1M_Name  = "_1M";
  part_2_Name   = "_2";
  part_2P_Name  = "_2P";
  part_2M_Name  = "_2M";
  part_3_Name   = "_3";
  part_3P_Name  = "_3P";
  part_3M_Name  = "_3M";
  pair_12_Name  = "_12";
  pair_13_Name  = "_13";
  pair_23_Name  = "_23";
  triplet_Name   = "_123";
  
  n1Name     = "n1";
  n2Name     = "n2";
  n3Name     = "n3";
  n1NwName   = "n1Nw";
  n2NwName   = "n2Nw";
  n3NwName   = "n3Nw";
  
  r1Name     = "r1";
  r2Name     = "r2";
  r3Name     = "r3";
  r2r1Name   = "r2r1";
  c2Name     = "c2";
  c3Name     = "c3";
  d3Name     = "d3";
  p3Name     = "p3";
  
  intR2Name       = "intR2";
  binCorrName     = "binCorr";
  intBinCorrName  = "intBinCorr";
  
  
  _title_counts = "yield";
  _title_m0     = "M_{0}";
  _title_m1     = "M_{1}";
  _title_m2     = "M_{2}";
  _title_m3     = "M_{3}";
  _title_m4     = "V0Centrality";
  _title_m5     = "TrkCentrality";
  _title_m6     = "SpdCentrality";
  
  _title_phi_1       = "#varphi_{1} (radian)";
  _title_phi_2       = "#varphi_{2} (radian)";
  _title_phi_3       = "#varphi_{3} (radian)";
  _title_phi_12   = "#eta_{1}#times#varphi_{1}#times#eta_{2}#times#varphi_{2}";
  _title_phi_13   = "#eta_{1}#times#varphi_{1}#times#eta_{3}#times#varphi_{3}";
  _title_phi_23   = "#eta_{2}#times#varphi_{2}#times#eta_{3}#times#varphi_{3}";
  
  vsZ         = "_vsZ";
  vsM         = "_vsM";
  vsPt         = "_vsPt";
  vsPhi        = "_vsPhi"; 
  vsEta        = "_vsEta"; 
  vsEtaPhi     = "vsPhi"; 
  
  
  if (_useWeights)
    {
    int iZ, iPhi, iPt;
    int iZ1,iPhi1,iPt1;
    int a, b;
    if (_weight_1P)
      {
      _correctionWeight_1P = new float[_nBins_vertexZ*_nBins_etaPhi_1*_nBins_pt_1];
      _correctionWeight_1M = new float[_nBins_vertexZ*_nBins_etaPhi_1*_nBins_pt_1];
      a = _nBins_pt_1;
      b = _nBins_etaPhi_1*_nBins_pt_1;
      for (iZ=0,iZ1=1; iZ<_nBins_vertexZ; iZ++, iZ1++)
        {
        for (iPhi=0,iPhi1=1; iPhi<_nBins_etaPhi_1; iPhi++, iPhi1++)
          {
          for (iPt=0,iPt1=1; iPt<_nBins_pt_1; iPt++, iPt1++)
            {
            _correctionWeight_1P[iZ*b+iPhi*a+iPt] = _weight_1P->GetBinContent(iZ1,iPhi1,iPt1);
            _correctionWeight_1M[iZ*b+iPhi*a+iPt] = _weight_1M->GetBinContent(iZ1,iPhi1,iPt1);
            }      
          }
        }
      } // _weight_1
    else
      {
      AliError("AliAnalysisTask3PCorrelations:: _weight_1 is a null pointer.");
      return;
      }
    if (_weight_2P)
      {
      _correctionWeight_2P = new float[_nBins_vertexZ*_nBins_etaPhi_2*_nBins_pt_2];
      _correctionWeight_2M = new float[_nBins_vertexZ*_nBins_etaPhi_2*_nBins_pt_2];
      a = _nBins_pt_2;
      b = _nBins_etaPhi_2*_nBins_pt_2;
      for (iZ=0,iZ1=1; iZ<_nBins_vertexZ; iZ++, iZ1++)
        {
        for (iPhi=0,iPhi1=1; iPhi<_nBins_etaPhi_2; iPhi++, iPhi1++)
          {
          for (iPt=0,iPt1=1; iPt<_nBins_pt_2; iPt++, iPt1++)
            {
            _correctionWeight_2P[iZ*b+iPhi*a+iPt] = _weight_2P->GetBinContent(iZ1,iPhi1,iPt1);
            _correctionWeight_2M[iZ*b+iPhi*a+iPt] = _weight_2M->GetBinContent(iZ1,iPhi1,iPt1);
            }      
          }
        }
      } // _weight_2
    else
      {
      AliError("AliAnalysisTask3PCorrelations:: _weight_2 is a null pointer.");
      return;
      }
    if (_weight_3P)
      {
      _correctionWeight_3P = new float[_nBins_vertexZ*_nBins_etaPhi_3*_nBins_pt_3];
      _correctionWeight_3M = new float[_nBins_vertexZ*_nBins_etaPhi_3*_nBins_pt_3];
      a = _nBins_pt_3;
      b = _nBins_etaPhi_3*_nBins_pt_3;
      for (iZ=0,iZ1=1; iZ<_nBins_vertexZ; iZ++, iZ1++)
        {
        for (iPhi=0,iPhi1=1; iPhi<_nBins_etaPhi_3; iPhi++, iPhi1++)
          {
          for (iPt=0,iPt1=1; iPt<_nBins_pt_3; iPt++, iPt1++)
            {
            _correctionWeight_3P[iZ*b+iPhi*a+iPt] = _weight_3P->GetBinContent(iZ1,iPhi1,iPt1);
            _correctionWeight_3M[iZ*b+iPhi*a+iPt] = _weight_3M->GetBinContent(iZ1,iPhi1,iPt1);
            }      
          }
        }
      } // _weight_2
    else
      {
      AliError("AliAnalysisTask3PCorrelations:: _weight_2 is a null pointer.");
      return;
      }
    }
  
  createHistograms();
  PostData(0,_outputHistoList);
  
  cout<< "AliAnalysisTask3PCorrelations::CreateOutputObjects() DONE " << endl;
  
}

void  AliAnalysisTask3PCorrelations::createHistograms()
{
  AliInfo(" AliAnalysisTask3PCorrelations::createHistoHistograms() Creating Event Histos");
  TString name;
  
  // bin index : what it is...
  //         0 :  number of event submitted
  //         1 :  number accepted by centrality cut
  //         2 :  number accepted by centrality cut and z cut
  //         3 :  total number of particles that satisfied filter 1
  //         4 :  total number of particles that satisfied filter 2
  //         5 :  total number of particles that satisfied filter 3
  name = "eventAccounting";
  _eventAccounting      = createHisto1D(name,name,10, -0.5, 9.5, "event Code", _title_counts);
  
  name = "m0"; _m0      = createHisto1D(name,name,_nBins_M1, _min_M1, _max_M1, _title_m0, _title_counts);
  name = "m1"; _m1      = createHisto1D(name,name,_nBins_M1, _min_M1, _max_M1, _title_m1, _title_counts);
  name = "m2"; _m2      = createHisto1D(name,name,_nBins_M2, _min_M2, _max_M2, _title_m2, _title_counts);
  name = "m3"; _m3      = createHisto1D(name,name,_nBins_M3, _min_M3, _max_M3, _title_m3, _title_counts);
  name = "m4"; _m4      = createHisto1D(name,name,_nBins_M4, _min_M4, _max_M4, _title_m4, _title_counts);
  name = "m5"; _m5      = createHisto1D(name,name,_nBins_M5, _min_M5, _max_M5, _title_m5, _title_counts);
  name = "m6"; _m6      = createHisto1D(name,name,_nBins_M6, _min_M6, _max_M6, _title_m6, _title_counts);
  name = "zV"; _vertexZ = createHisto1D(name,name,_nBins_vertexZ, _min_vertexZ, _max_vertexZ, "z-Vertex (cm)", _title_counts);
  
  TString title_vsM;
  //_centralityMethod
  switch (_centralityMethod)
  {
    case 0: _nBins_vsM = _nBins_M0; _min_vsM = _min_M0; _max_vsM = _max_M0; _width_vsM = _width_M0; title_vsM = _title_m0; break;
    case 1: _nBins_vsM = _nBins_M1; _min_vsM = _min_M1; _max_vsM = _max_M1; _width_vsM = _width_M1; title_vsM = _title_m1; break;
    case 2: _nBins_vsM = _nBins_M2; _min_vsM = _min_M2; _max_vsM = _max_M2; _width_vsM = _width_M2; title_vsM = _title_m2; break;
    case 3: _nBins_vsM = _nBins_M3; _min_vsM = _min_M3; _max_vsM = _max_M3; _width_vsM = _width_M3; title_vsM = _title_m3; break;
    case 4: _nBins_vsM = _nBins_M4; _min_vsM = _min_M4; _max_vsM = _max_M4; _width_vsM = _width_M4; title_vsM = _title_m4; break;
    case 5: _nBins_vsM = _nBins_M5; _min_vsM = _min_M5; _max_vsM = _max_M5; _width_vsM = _width_M5; title_vsM = _title_m5; break;
    case 6: _nBins_vsM = _nBins_M6; _min_vsM = _min_M6; _max_vsM = _max_M6; _width_vsM = _width_M6; title_vsM = _title_m6; break;
  }
  
  AliInfo(" AliAnalysisTask3PCorrelations::createHistoHistograms() Creating Part 1 Histos");
  
  if (_singlesOnly) 
    {
    name = n1Name+part_1_Name+vsPt;               _n1_1_vsPt               = createHisto1F(name,name, _nBins_pt_1,  _min_pt_1,  _max_pt_1,   "pt1", "counts");
    name = "dedxVsP_1";                           _dedxVsP_1               = createHisto2F(name,name,400,-2.,2.,120,0.,120.,"p (GeV/c)", "dedx", "counts");
    name = "betaVsP_1";                           _betaVsP_1               = createHisto2F(name,name,400,-2.,2.,120,0.5,1.1,"p (GeV/c)", "beta", "counts");
    name = n1Name+part_1P_Name+vsZ+vsEtaPhi+vsPt; _n1_1P_vsZVsEtaVsPhiVsPt = createHisto3F(name,name, _nBins_vertexZ,_min_vertexZ,_max_vertexZ, _nBins_etaPhi_1, 0., double(_nBins_etaPhi_1), _nBins_pt_1, _min_pt_1, _max_pt_1, "zVertex", _title_phi_1,  "pt1");
    name = n1Name+part_1M_Name+vsZ+vsEtaPhi+vsPt; _n1_1M_vsZVsEtaVsPhiVsPt = createHisto3F(name,name, _nBins_vertexZ,_min_vertexZ,_max_vertexZ, _nBins_etaPhi_1, 0., double(_nBins_etaPhi_1), _nBins_pt_1, _min_pt_1, _max_pt_1, "zVertex", _title_phi_1,  "pt1");
    
    name = n1Name+part_2_Name+vsPt;               _n1_2_vsPt               = createHisto1F(name,name, _nBins_pt_2,  _min_pt_2,  _max_pt_2,   "pt2",  "counts");
    name = "dedxVsP_2";                           _dedxVsP_2               = createHisto2F(name,name,400,-2.,2.,120,0.,120.,"p (GeV/c)", "dedx", "counts");
    name = "betaVsP_2";                           _betaVsP_2               = createHisto2F(name,name,400,-2.,2.,120,0.5,1.1,"p (GeV/c)", "beta", "counts");
    name = n1Name+part_2P_Name+vsZ+vsEtaPhi+vsPt; _n1_2P_vsZVsEtaVsPhiVsPt = createHisto3F(name,name, _nBins_vertexZ,_min_vertexZ,_max_vertexZ, _nBins_etaPhi_2, 0., double(_nBins_etaPhi_2), _nBins_pt_2, _min_pt_2, _max_pt_2, "zVertex", _title_phi_2,  "pt2");
    name = n1Name+part_2M_Name+vsZ+vsEtaPhi+vsPt; _n1_2M_vsZVsEtaVsPhiVsPt = createHisto3F(name,name, _nBins_vertexZ,_min_vertexZ,_max_vertexZ, _nBins_etaPhi_2, 0., double(_nBins_etaPhi_2), _nBins_pt_2, _min_pt_2, _max_pt_2, "zVertex", _title_phi_2,  "pt2");
    
    name = n1Name+part_3_Name+vsPt;               _n1_3_vsPt               = createHisto1F(name,name, _nBins_pt_3,  _min_pt_3,  _max_pt_3,   "pt3",  "counts");
    name = "dedxVsP_3";                           _dedxVsP_3               = createHisto2F(name,name,400,-2.,2.,120,0.,120.,"p (GeV/c)", "dedx", "counts");
    name = "betaVsP_3";                           _betaVsP_3               = createHisto2F(name,name,400,-2.,2.,120,0.5,1.1,"p (GeV/c)", "beta", "counts");
    name = n1Name+part_3P_Name+vsZ+vsEtaPhi+vsPt; _n1_3P_vsZVsEtaVsPhiVsPt = createHisto3F(name,name, _nBins_vertexZ,_min_vertexZ,_max_vertexZ, _nBins_etaPhi_3, 0., double(_nBins_etaPhi_3), _nBins_pt_3, _min_pt_3, _max_pt_3, "zVertex", _title_phi_3,  "pt3");
    name = n1Name+part_3M_Name+vsZ+vsEtaPhi+vsPt; _n1_3M_vsZVsEtaVsPhiVsPt = createHisto3F(name,name, _nBins_vertexZ,_min_vertexZ,_max_vertexZ, _nBins_etaPhi_3, 0., double(_nBins_etaPhi_3), _nBins_pt_3, _min_pt_3, _max_pt_3, "zVertex", _title_phi_3,  "pt3");
    }
  else
    {
    name = n1Name+part_1_Name+vsPhi;  _n1_1_vsPhi   = createHisto1F(name,name, _nBins_phi_1, _min_phi_1, _max_phi_1, _title_phi_1,  "counts");
    name = n1Name+part_2_Name+vsPhi;  _n1_2_vsPhi   = createHisto1F(name,name, _nBins_phi_2, _min_phi_2, _max_phi_2, _title_phi_2,  "counts");
    name = n1Name+part_3_Name+vsPhi;  _n1_3_vsPhi   = createHisto1F(name,name, _nBins_phi_3, _min_phi_3, _max_phi_3, _title_phi_3,  "counts");
    
    name = n1Name+part_1_Name + vsM;  _n1_1_vsM     = createProfile(name,name, _nBins_vsM,   _min_vsM,   _max_vsM,  title_vsM,      "<N_{1}>");
    name = n1Name+part_2_Name + vsM;  _n1_2_vsM     = createProfile(name,name, _nBins_vsM,   _min_vsM,   _max_vsM,  title_vsM,      "<N_{2}>");
    name = n1Name+part_3_Name + vsM;  _n1_3_vsM     = createProfile(name,name, _nBins_vsM,   _min_vsM,   _max_vsM,  title_vsM,      "<N_{3}>");
    
    name = n1NwName+part_1_Name+vsM;  _n1Nw_1_vsM   = createProfile(name,name, _nBins_vsM,   _min_vsM,   _max_vsM,  title_vsM,      "<N_{1}>");
    name = n1NwName+part_2_Name+vsM;  _n1Nw_2_vsM   = createProfile(name,name, _nBins_vsM,   _min_vsM,   _max_vsM,  title_vsM,      "<N_{2}>");
    name = n1NwName+part_3_Name+vsM;  _n1Nw_3_vsM   = createProfile(name,name, _nBins_vsM,   _min_vsM,   _max_vsM,  title_vsM,      "<N_{3}>");
    
    name = n2Name+pair_12_Name+vsPhi; _n2_12_vsPhi  = createHisto2F(name,name, _nBins_phi_1, _min_phi_1, _max_phi_1, _nBins_phi_2, _min_phi_2, _max_phi_2, _title_phi_1, _title_phi_2, "<N_{12}>");
    name = n2Name+pair_13_Name+vsPhi; _n2_13_vsPhi  = createHisto2F(name,name, _nBins_phi_1, _min_phi_1, _max_phi_1, _nBins_phi_3, _min_phi_3, _max_phi_3, _title_phi_2, _title_phi_3, "<N_{13}>");        
    name = n2Name+pair_23_Name+vsPhi; _n2_23_vsPhi  = createHisto2F(name,name, _nBins_phi_2, _min_phi_2, _max_phi_2, _nBins_phi_3, _min_phi_3, _max_phi_3, _title_phi_1, _title_phi_3, "<N_{23}>");
    name = n3Name+triplet_Name+vsPhi; _n3_123_vsPhi = createHisto3F(name,name, _nBins_phi_1, _min_phi_1, _max_phi_1, _nBins_phi_2, _min_phi_2, _max_phi_2, _nBins_phi_3, _min_phi_3, _max_phi_3, _title_phi_1, _title_phi_2, _title_phi_3);        
    
    name = n2Name+pair_12_Name+vsM;  _n2_12_vsM     = createProfile(name,name, _nBins_vsM,   _min_vsM,   _max_vsM,  title_vsM,      "<N_{12}>");
    name = n2Name+pair_13_Name+vsM;  _n2_13_vsM     = createProfile(name,name, _nBins_vsM,   _min_vsM,   _max_vsM,  title_vsM,      "<N_{13}>");
    name = n2Name+pair_23_Name+vsM;  _n2_23_vsM     = createProfile(name,name, _nBins_vsM,   _min_vsM,   _max_vsM,  title_vsM,      "<N_{23}>");
    name = n3Name+pair_23_Name+vsM;  _n3_123_vsM    = createProfile(name,name, _nBins_vsM,   _min_vsM,   _max_vsM,  title_vsM,      "<N_{123}>");
    name = n2NwName+pair_12_Name+vsM; _n2Nw_12_vsM  = createProfile(name,name, _nBins_vsM,   _min_vsM,   _max_vsM,  title_vsM,      "<N_{12}>");
    name = n2NwName+pair_13_Name+vsM; _n2Nw_13_vsM  = createProfile(name,name, _nBins_vsM,   _min_vsM,   _max_vsM,  title_vsM,      "<N_{13}>");
    name = n2NwName+pair_23_Name+vsM; _n2Nw_23_vsM  = createProfile(name,name, _nBins_vsM,   _min_vsM,   _max_vsM,  title_vsM,      "<N_{23}>");
    name = n3NwName+triplet_Name+vsM; _n3Nw_123_vsM = createProfile(name,name, _nBins_vsM,   _min_vsM,   _max_vsM,  title_vsM,      "<N_{123}>");
    
    }
  
  AliInfo(" AliAnalysisTask3PCorrelations::createHistoHistograms() All Done"); 
}
//-----------------------//

void  AliAnalysisTask3PCorrelations::finalizeHistograms()
{
  
  AliInfo("AliAnalysisTask3PCorrelations::finalizeHistograms() starting");
  AliInfo(Form("CorrelationAnalyzers::finalizeHistograms()   _eventCount : %d",int(_eventCount)));
  
  if (_singlesOnly) 
    {
    fillHistoWithArray(_n1_1_vsPt,               __n1_1_vsPt,         _nBins_pt_1);
    fillHistoWithArray(_n1_1P_vsZVsEtaVsPhiVsPt, __n1_1P_vsZEtaPhiPt, _nBins_vertexZ, _nBins_etaPhi_1, _nBins_pt_1);
    fillHistoWithArray(_n1_1M_vsZVsEtaVsPhiVsPt, __n1_1M_vsZEtaPhiPt, _nBins_vertexZ, _nBins_etaPhi_1, _nBins_pt_1);
    
    fillHistoWithArray(_n1_2_vsPt,               __n1_2_vsPt,         _nBins_pt_2);
    fillHistoWithArray(_n1_2P_vsZVsEtaVsPhiVsPt, __n1_2P_vsZEtaPhiPt, _nBins_vertexZ, _nBins_etaPhi_2, _nBins_pt_2);
    fillHistoWithArray(_n1_2M_vsZVsEtaVsPhiVsPt, __n1_2M_vsZEtaPhiPt, _nBins_vertexZ, _nBins_etaPhi_2, _nBins_pt_2);
    
    fillHistoWithArray(_n1_3_vsPt,               __n1_3_vsPt,         _nBins_pt_3);
    fillHistoWithArray(_n1_3P_vsZVsEtaVsPhiVsPt, __n1_3P_vsZEtaPhiPt, _nBins_vertexZ, _nBins_etaPhi_3, _nBins_pt_3);
    fillHistoWithArray(_n1_3M_vsZVsEtaVsPhiVsPt, __n1_3M_vsZEtaPhiPt, _nBins_vertexZ, _nBins_etaPhi_3, _nBins_pt_3);
    
    }
  else
    {
    fillHistoWithArray(_n1_1_vsPhi,    __n1_1_vsPhi,   _nBins_phi_1);
    fillHistoWithArray(_n1_2_vsPhi,    __n1_2_vsPhi,   _nBins_phi_2);
    fillHistoWithArray(_n1_3_vsPhi,    __n1_3_vsPhi,   _nBins_phi_3);
    fillHistoWithArray(_n2_12_vsPhi,   __n2_12_vsPhi,  _nBins_phi_1, _nBins_phi_2);
    fillHistoWithArray(_n2_13_vsPhi,   __n2_13_vsPhi,  _nBins_phi_1, _nBins_phi_3);
    fillHistoWithArray(_n2_23_vsPhi,   __n2_23_vsPhi,  _nBins_phi_2, _nBins_phi_3);
    fillHistoWithArray(_n3_123_vsPhi,  __n3_123_vsPhi, _nBins_phi_1, _nBins_phi_2, _nBins_phi_3);
    }
  AliInfo("AliAnalysisTask3PCorrelations::finalizeHistograms()  Done ");
}
//--------------//


void  AliAnalysisTask3PCorrelations::UserExec(Option_t */*option*/)
{
  // cout << "AliAnalysisTask3PCorrelations::UserExec(Option_t *option) - Starting!!!!" << endl;
  
  int    k1=0,k2=0,k3=0;
  int    iPhi=0, iEta=0, iEtaPhi=0, iPt=0, charge=0;
  float  q=0, p=0, phi=0, pt=0, eta=0, corr=0, dedx=0, px=0, py=0, pz=0;
  int    ij_12=0, ij_13=0, ij_23=0, ijk_123=0;
  int    id_1=0, id_2=0, id_3=0, iPhi_1=0, iPhi_2=0, iPhi_3=0;
  float  corr_1=0, corr_2=0, corr_3=0;
  float  corr_12=0, corr_13=0, corr_23=0, corr_123=0;
  int    iVertex=0, iVertexP1=0, iVertexP2=0, iVertexP3=0;
  int    iZEtaPhiPt=0;
  //float  massElecSq = 2.5e-7=0;
  double b[2];
  double bCov[3];
  const  AliAODVertex*	vertex=0;
  int    nClus=0;
  bool   bitOK=0;
  
  // count all events looked at here
  _eventCount++;
  
  if (_eventAccounting)
    {
    _eventAccounting->Fill(0);// count all calls to this function
    }
  else
    {
    cout << "AliAnalysisTask3PCorrelations::UserExec(Option_t *option) - !_eventAccounting" << endl;
    cout << "AliAnalysisTask3PCorrelations::UserExec(Option_t *option) - Skip this task" << endl;
    return;
    }
  
  //Get pointer to current event
  //fESDEvent = dynamic_cast<AliESDEvent*> (InputEvent());
  fAODEvent = dynamic_cast<AliAODEvent*> (InputEvent());
  
  // cout << "AliAnalysisTask3PCorrelations::UserExec(Option_t *option) - 2" << endl;
  
  if(!fAODEvent)
    {
    AliError("AliAnalysisTask3PCorrelations::Exec(Option_t *option) !fAODEvent"); 
    return;
    }
  
  
  float v0Centr  = -999.;
  float trkCentr = -999.;
  float spdCentr = -999.;
  float vertexX  = -999;
  float vertexY  = -999;
  float vertexZ  = -999;
  float vertexXY = -999;
  float dcaZ     = -999;
  float dcaXY    = -999;
  float centrality = -999;
  k1 = k2 = k3 = 0;
  __n1_1 = __n1_2 = __n1_3 = __n1Nw_1 = __n1Nw_2 = __n1Nw_3 = 0;
  
  _eventAccounting->Fill(1);// count all calls to this function with a valid pointer
  
  //Centrality
  AliCentrality* centralityObject =  ((AliVAODHeader*)fAODEvent->GetHeader())->GetCentralityP();
  if (centralityObject)
    {
    // cout << "AliAnalysisTask3PCorrelations::UserExec(Option_t *option) - 6" << endl;
    
    v0Centr  = centralityObject->GetCentralityPercentile("V0M");
    trkCentr = centralityObject->GetCentralityPercentile("TRK"); 
    spdCentr = centralityObject->GetCentralityPercentile("CL1");
    // cout << "AliAnalysisTask3PCorrelations::UserExec(Option_t *option) - 7" << endl;
    
    }
  // cout << "AliAnalysisTask3PCorrelations::UserExec(Option_t *option) - 8" << endl;
  
  _nTracks  = fAODEvent->GetNumberOfTracks();
  _mult0    = 0;
  _mult1    = 0;
  _mult2    = 0;
  _mult3    = _nTracks; 
  _mult4    = v0Centr;
  _mult5    = trkCentr;
  _mult6    = spdCentr;
  _field    = fAODEvent->GetMagneticField(); 
  
  // cout << "AliAnalysisTask3PCorrelations::UserExec(Option_t *option) - _field:" << _field << endl;
  
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
  
  // cout << "AliAnalysisTask3PCorrelations::UserExec(Option_t *option) - 10" << endl;
  
  //filter on centrality
  if ( centrality < _centralityMin ||
      centrality > _centralityMax)
    {
    // cout << "AliAnalysisTask3PCorrelations::UserExec(Option_t *option) - 11" << endl;
    
    // we dont analyze this event here... 
    return;
    }
  _eventAccounting->Fill(2);// count all events with right centrality
                            // cout << "AliAnalysisTask3PCorrelations::UserExec(Option_t *option) - 12" << endl;
  
  // filter on z and xy vertex
  vertex = (AliAODVertex*) fAODEvent->GetPrimaryVertexSPD();
  if (!vertex || vertex->GetNContributors()<1)
    {
    vertexZ  = -999;
    vertexXY = -999;
    // cout << "AliAnalysisTask3PCorrelations::UserExec(Option_t *option) - No valid vertex object or poor vertex" << endl;
    }
  else
    { 
      vertexX = vertex->GetX();
      vertexY = vertex->GetY();
      vertexZ = vertex->GetZ();
      vertexXY = sqrt(vertexX*vertexX+vertexY*vertexY);
    }
  // cout << "AliAnalysisTask3PCorrelations::UserExec(Option_t *option) - 13" << endl;
  
  // cout << "vertexZ : " << vertexZ  << endl;
  // cout << "vertexXY: " << vertexXY << endl;
  // cout << "_vertexZMin: " << _vertexZMin << endl;
  // cout << "_vertexZMax: " << _vertexZMax << endl;
  // cout << "_vertexXYMin: " << _vertexXYMin << endl;
  // cout << "_vertexXYMax: " << _vertexXYMax << endl;
  
  if (!vertex ||
      vertexZ  < _vertexZMin  || 
      vertexZ  > _vertexZMax  ||
      vertexXY < _vertexXYMin || 
      vertexXY > _vertexXYMax)
    return;
  
  
  // cout << "AliAnalysisTask3PCorrelations::UserExec(Option_t *option) - 14" << endl;
  
  iVertex = int((vertexZ-_min_vertexZ)/_width_vertexZ);
  iVertexP1 = iVertex*_nBins_etaPhiPt_1;
  iVertexP2 = iVertex*_nBins_etaPhiPt_2;
  iVertexP3 = iVertex*_nBins_etaPhiPt_3;
  if (iVertex<0 || iVertex>=_nBins_vertexZ)
    {
    AliError("AliAnalysisTask3PCorrelations::Exec(Option_t *option) iVertex<0 || iVertex>=_nBins_vertexZ ");
    return;
    }
  _eventAccounting->Fill(3);// count all calls to this function with a valid pointer
                            // cout << "AliAnalysisTask3PCorrelations::UserExec(Option_t *option) - 15" << endl;
  
  
  // loop over all particles for singles
  //debug() << "_nTracks:"<< _nTracks << endl;
  //reset single particle counters
  for (int iTrack=0; iTrack< _nTracks; iTrack++)
    {
    // cout << "AliAnalysisTask3PCorrelations::UserExec(Option_t *option) - 16 iTrack:" << iTrack << endl;
    
    AliAODTrack * t = (AliAODTrack *) fAODEvent->GetTrack(iTrack);
    if (!t) 
      {
      AliError(Form("AliAnalysisTask3PCorrelations::Exec(Option_t *option) No track ofr iTrack=%d", iTrack));
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
    if (!bitOK ||
        dedx  < _dedxMin   ||
        dedx  > _dedxMax   ||
        nClus < _nClusterMin) continue;
    
    // cout << "_trackFilterBit:" << _trackFilterBit << " Track returns:" << bitOK << endl;
    // cout << "   pt:" << pt   << " _min_pt_1:" << _min_pt_1 << " _max_pt_1:" << _max_pt_1<< endl;
    // cout << "  phi:" << phi  << endl;
    // cout << "  eta:" << eta  << " _min_eta_1:" << _min_eta_1 << " _max_eta_1:" << _max_eta_1<< endl;
    // cout << " dedx:" << dedx << " _dedxMin:" << _dedxMin << " _dedxMax:" << _dedxMax << endl;
    // cout << "nclus:" << nClus<< " _nClusterMin:" << _nClusterMin << endl;
    if (pt       >=  _min_pt_1  && 
        pt       <   _max_pt_1  && 
        eta      >=  _min_eta_1 && 
        eta      <   _max_eta_1) 
      {
      // cout << "AliAnalysisTask3PCorrelations::UserExec(Option_t *option) - check vertex  for 1:" << endl;
      // Get the dca information
      AliAODTrack* clone = (AliAODTrack*) t->Clone("trk_clone"); //need clone, in order not to change track parameters
      if (clone->PropagateToDCA(vertex, _field, 100., b, bCov) )
        {
        dcaXY = b[0];
        dcaZ  = b[1];
        } 
      else
        {
        dcaXY = -999999;
        dcaZ  = -999999;
        }
      delete clone;
      // cout << "1 dcaZ:" << dcaZ << " _dcaZMin:" << _dcaZMin << " _dcaZMax:" << _dcaZMax << endl;
      // cout << "1 dcaXY:" << dcaXY << " _dcaXYMin:" << _dcaXYMin << " _dcaXYMax:" << _dcaXYMax << endl;
      
      // skip track if DCA too large
      if (dcaZ     >=  _dcaZMin && 
          dcaZ     <   _dcaZMax && 
          dcaXY    >=  _dcaXYMin && 
          dcaXY    <   _dcaXYMax)
        continue; //track does not have a valid DCA
                  // cout << "keep track:" << endl;
      
      iPhi   = int( phi/_width_phi_1);
      // cout << " AliAnalysisTask3PCorrelations::analyze(Event * event) -1- iTrack:" << iTrack<< endl<< "pt:" << pt << " phi:" <<  phi << " eta:" << eta << endl;
      if (iPhi<0 || iPhi>=_nBins_phi_1 ) 
        {
        AliWarning("AliAnalysisTask3PCorrelations::analyze() iPhi<0 || iPhi>=_nBins_phi_1");
        return;
        }
      
      iEta    = int((eta-_min_eta_1)/_width_eta_1); 
      if (iEta<0 || iEta>=_nBins_eta_1) 
        {
        AliWarning(Form("AliAnalysisTask3PCorrelations::analyze(AliceEvent * event) Mismatched iEta: %d", iEta));
        continue;
        }
      iPt     = int((pt -_min_pt_1 )/_width_pt_1 ); 
      if (iPt<0  || iPt >=_nBins_pt_1)
        {
        AliWarning(Form("AliAnalysisTask3PCorrelations::analyze(AliceEvent * event) Mismatched iPt: %d",iPt));
        continue;
        }
      iEtaPhi    = iEta*_nBins_phi_1+iPhi;
      iZEtaPhiPt = iVertexP1 + iEtaPhi*_nBins_pt_1 + iPt;
      
      if (charge>0 && _correctionWeight_1P)
        corr = _correctionWeight_1P[iZEtaPhiPt];
      else if (charge<0 && _correctionWeight_1M)
        corr = _correctionWeight_1M[iZEtaPhiPt];
      else
        corr = 1;
      // cout << "all good; process track:" << endl;
      if (_singlesOnly)
        {
        __n1_1_vsPt[iPt]                += corr;
        if (charge>0) 
          __n1_1P_vsZEtaPhiPt[iZEtaPhiPt] += corr;
        else if (charge>0) 
          __n1_1M_vsZEtaPhiPt[iZEtaPhiPt] += corr;
        }
      else
        {
        _id_1[k1]            = iTrack;     
        //_charge_1[k1]        = charge;
        _iPhi_1[k1]          = iPhi; 
        _correction_1[k1]    = corr; 
        __n1_1               += corr;
        __n1_1_vsPhi[iPhi]   += corr; 
        __n1Nw_1             += 1;
        ++k1;
        if (k1>=arraySize)
          {
          AliError(Form("AliAnalysisTask3PCorrelations::analyze(AliceEvent * event) k1 >=arraySize; arraySize: %d",arraySize));
          return;
          }
        }
      }
    
    else if (pt  >= _min_pt_2  && 
             pt  <  _max_pt_2  && 
             eta >= _min_eta_2 && 
             eta <  _max_eta_2)  
      {
      // Get the dca information
      // cout << "AliAnalysisTask3PCorrelations::UserExec(Option_t *option) - check vertex for 2:" << endl;
      AliAODTrack* clone = (AliAODTrack*) t->Clone("trk_clone"); //need clone, in order not to change track parameters
      if (clone->PropagateToDCA(vertex, _field, 100., b, bCov) )
        {
        dcaXY = b[0];
        dcaZ  = b[1];
        } 
      else
        {
        dcaXY = -999999;
        dcaZ  = -999999;
        }
      delete clone;
      // cout << "2 dcaZ:" << dcaZ << " _dcaZMin:" << _dcaZMin << " _dcaZMax:" << _dcaZMax << endl;
      // cout << "2 dcaXY:" << dcaXY << " _dcaXYMin:" << _dcaXYMin << " _dcaXYMax:" << _dcaXYMax << endl;
      // skip track if DCA too large
      if (dcaZ     >=  _dcaZMin && 
          dcaZ     <   _dcaZMax && 
          dcaXY    >=  _dcaXYMin && 
          dcaXY    <   _dcaXYMax)
        continue; //track does not have a valid DCA
      
      iPhi   = int( phi/_width_phi_2);
      //cout << " AliAnalysisTask3PCorrelations::analyze(Event * event) -1- iTrack:" << iTrack  << endl
      //<< "pt:" << pt << " phi:" <<  phi << " eta:" << eta << endl;
      if (iPhi<0 || iPhi>=_nBins_phi_2 ) 
        {
        AliWarning("AliAnalysisTask3PCorrelations::analyze() iPhi<0 || iPhi>=_nBins_phi_1");
        return;
        }
      
      iEta    = int((eta-_min_eta_2)/_width_eta_2);
      if (iEta<0 || iEta>=_nBins_eta_2) 
        {
        AliWarning(Form("AliAnalysisTask3PCorrelations::analyze(AliceEvent * event) Mismatched iEta: %d", iEta));
        continue;
        }
      iPt     = int((pt -_min_pt_2 )/_width_pt_2 ); 
      if (iPt<0  || iPt >=_nBins_pt_2)
        {
        AliWarning(Form("AliAnalysisTask3PCorrelations::analyze(AliceEvent * event) Mismatched iPt: %d",iPt));
        continue;
        }
      
      iEtaPhi    = iEta*_nBins_phi_2+iPhi;
      iZEtaPhiPt = iVertexP2 + iEtaPhi*_nBins_pt_2 + iPt;
      if (iZEtaPhiPt<0 || iZEtaPhiPt>=_nBins_zEtaPhiPt_2)
        {
        AliWarning("AliAnalysisTask3PCorrelations::analyze(AliceEvent * event) iZEtaPhiPt<0 || iZEtaPhiPt>=_nBins_zEtaPhiPt_2");
        continue;
        }
      
      if (charge>0 && _correctionWeight_2P)
        corr = _correctionWeight_2P[iZEtaPhiPt];
      else if (charge<0 && _correctionWeight_2M)
        corr = _correctionWeight_2M[iZEtaPhiPt];
      else
        corr = 1;
      
      if (_singlesOnly)
        {
        //_dedxVsP_2->Fill(p*q,dedx);
        __n1_2_vsPt[iPt]               += corr;          // cout << "step 15" << endl;
        if (charge>0) 
          __n1_2P_vsZEtaPhiPt[iZEtaPhiPt] += corr;
        else if (charge>0) 
          __n1_2M_vsZEtaPhiPt[iZEtaPhiPt] += corr;          }
      else
        {
        _id_2[k2]           = iTrack;         // cout << "step 1" << endl;
                                              //_charge_2[k2]       = charge;         // cout << "step 2" << endl;
        _iPhi_2[k2]         = iPhi;        // cout << "step 3" << endl;
        _correction_2[k2]   = corr;           // cout << "step 9" << endl;
        __n1_2              += corr;          // cout << "step 10" << endl;
        __n1Nw_2            += 1;
        __n1_2_vsPhi[iPhi]  += corr;          // cout << "step 11" << endl;
        ++k2;
        if (k2>=arraySize)
          {
          AliWarning(Form("-W- k2 >=arraySize; arraySize: %d",arraySize)); 
          return;
          }
        }
      }
    
    else if (pt  >= _min_pt_3  && 
             pt  <  _max_pt_3  && 
             eta >= _min_eta_3 && 
             eta <  _max_eta_3)  
      {
      // Get the dca information
      // cout << "AliAnalysisTask3PCorrelations::UserExec(Option_t *option) - check vertex for 2:" << endl;
      AliAODTrack* clone = (AliAODTrack*) t->Clone("trk_clone"); //need clone, in order not to change track parameters
      if (clone->PropagateToDCA(vertex, _field, 100., b, bCov) )
        {
        dcaXY = b[0];
        dcaZ  = b[1];
        } 
      else
        {
        dcaXY = -999999;
        dcaZ  = -999999;
        }
      delete clone;
      // skip track if DCA too large
      if (dcaZ     >=  _dcaZMin && 
          dcaZ     <   _dcaZMax && 
          dcaXY    >=  _dcaXYMin && 
          dcaXY    <   _dcaXYMax)
        continue; //track does not have a valid DCA
      
      iPhi   = int( phi/_width_phi_3);
      if (iPhi<0 || iPhi>=_nBins_phi_3 ) 
        {
        AliWarning("AliAnalysisTask3PCorrelations::analyze() iPhi<0 || iPhi>=_nBins_phi_3");
        return;
        }
      
      iEta    = int((eta-_min_eta_3)/_width_eta_3);
      if (iEta<0 || iEta>=_nBins_eta_3) 
        {
        AliWarning(Form("AliAnalysisTask3PCorrelations::analyze(AliceEvent * event) Mismatched iEta: %d", iEta));
        continue;
        }
      iEtaPhi    = iEta*_nBins_phi_3+iPhi;
      iZEtaPhiPt = iVertexP3 + iEtaPhi*_nBins_pt_3 + iPt;
      if (iZEtaPhiPt<0 || iZEtaPhiPt>=_nBins_zEtaPhiPt_3)
        {
        AliWarning("AliAnalysisTask3PCorrelations::analyze(AliceEvent * event) iZEtaPhiPt<0 || iZEtaPhiPt>=_nBins_zEtaPhiPt_3");
        continue;
        }
      if (charge>0 && _correctionWeight_3P)
        corr = _correctionWeight_3P[iZEtaPhiPt];
      else if (charge<0 && _correctionWeight_3M)
        corr = _correctionWeight_3M[iZEtaPhiPt];
      else
        corr = 1;
      
      if (_singlesOnly)
        {
        __n1_3_vsPt[iPt]               += corr;          // cout << "step 15" << endl;
        if (charge>0) 
          __n1_3P_vsZEtaPhiPt[iZEtaPhiPt] += corr;
        else if (charge>0) 
          __n1_3M_vsZEtaPhiPt[iZEtaPhiPt] += corr;          }
      else
        {
        _id_3[k3]            = iTrack;         // cout << "step 1" << endl;
                                               //_charge_3[k3]        = charge;         // cout << "step 2" << endl;
        _iPhi_3[k3]          = iPhi;        // cout << "step 3" << endl;
        _correction_3[k3]    = corr;           // cout << "step 9" << endl;
        __n1_3               += corr;          // cout << "step 10" << endl;
        __n1Nw_3             += 1;
        __n1_3_vsPhi[iPhi]   += corr;          // cout << "step 11" << endl;
        ++k3;
        if (k3>=arraySize)
          {
          AliWarning(Form("-W- k3 >=arraySize; arraySize: %d",arraySize)); 
          return;
          }
        }
      }
    
    // cout << "done with track" << endl;
    } //iTrack
  
  // cout << "Filling histograms now" << endl;
  _m0->Fill(_mult0);
  _m1->Fill(_mult1);
  _m2->Fill(_mult2);
  _m3->Fill(_mult3);
  _m4->Fill(_mult4);
  _m5->Fill(_mult5);
  _m6->Fill(_mult6);
  _vertexZ->Fill(vertexZ);
  
  
  //if singels only selected, do not fill pair histograms.
  if (_singlesOnly) 
    {
    
    }
  else
    {
    // reset pair counters
    _n1_1_vsM->Fill(    centrality, __n1_1);
    _n1Nw_1_vsM->Fill(  centrality, __n1Nw_1);
    _n1_2_vsM->Fill(    centrality, __n1_2);
    _n1Nw_2_vsM->Fill(  centrality, __n1Nw_2);
    _n1_3_vsM->Fill(    centrality, __n1_3);
    _n1Nw_3_vsM->Fill(  centrality, __n1Nw_3);

    __n2_12 = __n2Nw_12 = __n2_13 = __n2Nw_13 = __n2_23 = __n2Nw_23 = __n3_123 = __n3Nw_123 = 0;
    
    for (int i1=0; i1<k1; i1++)
      {
      // cout << "         i1:" << i1 << endl;
      id_1    = _id_1[i1]; 
      iPhi_1  = _iPhi_1[i1];
      corr_1  = _correction_1[i1];
      //1 and 2
      for (int i2=0;i2<k2;i2++)
        {        
          id_2 = _id_2[i2];
          if (id_1!=id_2)
            {
            iPhi_2  = _iPhi_2[i2];  
            corr_2  = _correction_2[i2];
            corr_12 = corr_1*corr_2;
            ij_12 = iPhi_1*_nBins_phi_1 + iPhi_2;
            __n2_12_vsPhi[ij_12] += corr_12;
            __n2_12              += corr_12;
            __n2Nw_12            += 1;
            for (int i3=0;i3<k3;i3++)
              {        
                id_3 = _id_3[i3];
                if (id_1!=id_3 && id_2!=id_3)
                  {
                  iPhi_3   = _iPhi_3[i3];  
                  corr_3   = _correction_3[i3];
                  corr_123 = corr_12*corr_3;
                  ijk_123  = ij_12 *_nBins_phi_2 + iPhi_3;
                  __n3_123_vsPhi[ijk_123] += corr_123;
                  __n3_123                += corr_123;
                  __n3Nw_123              += 1;
                  }
              } //i3       
            }
        } //i2       
          //1 and 3
      for (int i3=0;i3<k3;i3++)
        {        
          id_3 = _id_3[i3];
          if (id_1!=id_3)
            {
            iPhi_3  = _iPhi_3[i3];  
            corr_3  = _correction_3[i3];
            corr_13 = corr_1*corr_3;
            ij_13   = iPhi_1*_nBins_phi_1 + iPhi_3;
            __n2_13_vsPhi[ij_13] += corr_13;
            __n2_13              += corr_13;
            __n2Nw_13            += 1;
            }
        } //i3       
      } //i1  
    
    //2 & 3
    for (int i2=0; i2<k2; i2++)
      {
      id_2    = _id_2[i2]; 
      iPhi_2  = _iPhi_2[i2];
      corr_2  = _correction_2[i2];
      for (int i3=0;i3<k3;i3++)
        {        
          id_3 = _id_3[i3];
          if (id_2!=id_3)
            {
            iPhi_3  = _iPhi_3[i3];  
            corr_3  = _correction_3[i3];
            corr_23 = corr_2*corr_3;
            ij_23   = iPhi_2*_nBins_phi_2 + iPhi_3;
            __n2_23_vsPhi[ij_23] += corr_23;
            __n2_23              += corr_23;
            __n2Nw_23            += 1;
            }
        } //i2
      }
    _n2_12_vsM->Fill(   centrality, __n2_12);
    _n2_13_vsM->Fill(   centrality, __n2_13);
    _n2_23_vsM->Fill(   centrality, __n2_23);
    _n3_123_vsM->Fill(  centrality, __n3_123);
    _n2Nw_12_vsM->Fill( centrality, __n2Nw_12);
    _n2Nw_13_vsM->Fill( centrality, __n2Nw_13);
    _n2Nw_23_vsM->Fill( centrality, __n2Nw_23);
    _n3Nw_123_vsM->Fill(centrality, __n3_123);
    }
  
  // cout << "Event Done " << endl;
  PostData(0,_outputHistoList);
  
}

void   AliAnalysisTask3PCorrelations::FinishTaskOutput()
{
  cout << "AliAnalysisTask3PCorrelations::FinishTaskOutput() Starting." << endl;
  finalizeHistograms();
  PostData(0,_outputHistoList);
  cout << "AliAnalysisTask3PCorrelations::FinishTaskOutput() Done." << endl;
}

void   AliAnalysisTask3PCorrelations::Terminate(Option_t* /*option*/)
{
  // no ops
}


//Tools
//===================================================================================================
void  AliAnalysisTask3PCorrelations::fillHistoWithArray(TH1 * h, float * array, int size)
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

void  AliAnalysisTask3PCorrelations::fillHistoWithArray(TH2 * h, float * array, int size1, int size2)
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

void  AliAnalysisTask3PCorrelations::fillHistoWithArray(TH3 * h, float * array, int size1, int size2, int size3)
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

void  AliAnalysisTask3PCorrelations::fillHistoWithArray(TH1 * h, double * array, int size)
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

void  AliAnalysisTask3PCorrelations::fillHistoWithArray(TH2 * h, double * array, int size1, int size2)
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

void  AliAnalysisTask3PCorrelations::fillHistoWithArray(TH3 * h, double * array, int size1, int size2, int size3)
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
double *  AliAnalysisTask3PCorrelations::getDoubleArray(int size, double v)
{
  /// Allocate an array of type double with n values
  /// Initialize the array to the given value
  double * array = new double [size];
  for (int i=0;i<size;++i) array[i]=v;
  return array;
}

//________________________________________________________________________
float *  AliAnalysisTask3PCorrelations::getFloatArray(int size, float v)
{
  /// Allocate an array of type float with n values
  /// Initialize the array to the given value
  float * array = new float [size];
  for (int i=0;i<size;++i) array[i]=v;
  return array;
}


//________________________________________________________________________
TH1D * AliAnalysisTask3PCorrelations::createHisto1D(const TString &  name, const TString &  title, 
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
TH1D * AliAnalysisTask3PCorrelations::createHisto1D(const TString &  name, const TString &  title, 
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
TH2D * AliAnalysisTask3PCorrelations::createHisto2D(const TString &  name, const TString &  title, 
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
TH2D * AliAnalysisTask3PCorrelations::createHisto2D(const TString &  name, const TString &  title, 
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
TH1F * AliAnalysisTask3PCorrelations::createHisto1F(const TString &  name, const TString &  title, 
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
TH1F * AliAnalysisTask3PCorrelations::createHisto1F(const TString &  name, const TString &  title, 
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
TH2F * AliAnalysisTask3PCorrelations::createHisto2F(const TString &  name, const TString &  title, 
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
TH2F * AliAnalysisTask3PCorrelations::createHisto2F(const TString &  name, const TString &  title, 
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
TH3F * AliAnalysisTask3PCorrelations::createHisto3F(const TString &  name, const TString &  title, 
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
TProfile * AliAnalysisTask3PCorrelations::createProfile(const TString & name, const TString & description, 
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
TProfile * AliAnalysisTask3PCorrelations::createProfile(const TString &  name,const TString &  description, 
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


void   AliAnalysisTask3PCorrelations::addToList(TH1 *h)
{
  if (_outputHistoList)
    {
    _outputHistoList->Add(h);
    }
  else
    cout << "addToList(TH1 *h) _outputHistoList is null!!!!! Shoudl abort ship" << endl;
  
}



