//////////////////////////////////////////////////////////////////////
//
// $Id$
//
// Author: Emanuele Simili
//
//////////////////////////////////////////////////////////////////////
//_____________________________________________________________
//
// Description: 
//         the AliFlowAnalysisMaker class loops over AliFlowEvents, 
// calculates/extracts interesting observables, produces and fills the 
// related histograms for flow analysis. 
// Set AliFlowEvents file via SetInputFileName("file.root"), and the 
// histograms file via SetHistFileName("hist.root"), defalts are given. 
// It also calculates the phi weights to make the event plane isotropic 
// in the lab, and the bayesian weight for p.id. calculation, and saves 
// them in a file (set it via SetPhiWgtFileName("flowPhiWgt.hist.root")). 
// Those weights, if enabled, will be used only in the next run of the 
// analysis. In the first run or whenever the wgt file is missing, all 
// weights (phi & bayesian) are set to 1 .
// The most important features of the Flow Analysis are:
//  - Reaction Plane calculation: for different harmonics and selections
//   the R.P. is simply extracted via the appropriate AliFlowEvent methods.
//  - Correlation Analysis: particles are related to the R.P. and Vn 
//   coefficient are extracted. The AliFlowAnalysisMaker also removes 
//   auto-correlations of each particle with respect to the event plane, 
//   by subtracting track by track from the Q vector before the correlation 
//   study. Observable Flow is plotted.
//  - Differential Flow: for each harmonic and selection, flow values 
//   (Vn) are plotted versus rapidity/pseudorapidity and pT. 
//  - Resulution Estimate: event plane resolution is estimated via the 
//   sub-event method (res ~ sqrt(cos(psi1-psi2))), and flow values are 
//   then corrected for it. "True" Flow is plotted. 
//  - Integrated Flow: for the doubly-integrated (pT & eta) Vn values 
//   the error on the event plane resolution is folded in to the error. 
//  - Phi and Baiesian weights calculation (from event observables). 
//
// The AliFlowAnalysisMaker Class is adapted from the original 
// StFlowAnalysisMaker, succesfully used to study flow in the 
// STAR experiment at RICH.
// Original Authors:                 Raimond Snellings & Art Poskanzer
//

#include "AliFlowAnalysisMaker.h"
#include "AliFlowConstants.h"
#include "TVector2.h"
#include <iostream>
#include <cmath>     //PK: <cmath> needed for fabs() and fmod() on some systems.
using namespace std; //required for resolving the 'cout' symbol

ClassImp(AliFlowAnalysisMaker)
//-----------------------------------------------------------------------
AliFlowAnalysisMaker::AliFlowAnalysisMaker(const Char_t* name)
{
 // default constructor 

 pFlowSelect = new AliFlowSelection();
 mName = name ;
 InitDefault() ;
 cout << FlowAnalysis << " !!! No flowSelection defined .  Will use default cuts  (class AliFlowSelection) !!! " << endl ;
}
//-----------------------------------------------------------------------
AliFlowAnalysisMaker::AliFlowAnalysisMaker(const Char_t* name, const AliFlowSelection& flowSelect)
{
 //copy constructor

 pFlowSelect = new AliFlowSelection(flowSelect);  
 mName = name ;
 InitDefault() ;
}
//-----------------------------------------------------------------------
AliFlowAnalysisMaker::~AliFlowAnalysisMaker() 
{
 // default distructor (no actions) 
}
//-----------------------------------------------------------------------
void AliFlowAnalysisMaker::InitDefault()
{
 // sets some defaults for the analysis
 
 SetDebugg(1) ;  // Debug options (1,2,3)
 if(Debug0) { cout << FlowAnalysis << "Constructor" << endl ; }
 
 // input-output file (default output is given via the name given to the AliFlowAnalysis Object)
 SetInputFileName("AliFlowEVT.root") ;
 TString histfilename = mName ; histfilename += "Plot"  ; histfilename += ".root" ; 
 SetHistFileName(histfilename.Data());
 fFileNames.Clear() ;

 // internal settings (for event loop & printing purpose)
 FlowAnalysis = " [" ; FlowAnalysis += mName ; FlowAnalysis += "] - " ;
 spaces = " " ; for(Int_t s=0;s<(FlowAnalysis.Length()-1);s++) { spaces += " " ; }
 evtN = 0 ;
 pFlowEventsFile = 0 ;
 pFlowEvent      = 0 ;
 for(Int_t ii=0;ii<3;ii++){ vertex[ii] = 0 ; }
 // for labels histogram
 mLabelling = kTRUE ; 
 mLabellingFirst = kTRUE ;
 maxLabel = 0 ;

 // analysis settings 
 SetPhiWgtFileName("flowPhiWgt.hist.root");
 
 mV0           = kTRUE ;  // correlation analysis is done also for neutral secundary vertex
 mShuffle      = kFALSE ; // randomly reshuffles tracks
 mV1Ep1Ep2     = kFALSE ; // disabled by default
 mEtaSub       = kFALSE ; // eta subevents
 mPhiWgt       = kTRUE ;  // if flowPhiWgt file is there -> Phi Weights are used
 mBayWgt       = kFALSE ; // Bayesian Weights for P.Id. not used  
 mRePid        = kFALSE ; // recalculates the P.Id. (becomes kTRUE if new bayesian weights are plugged in)  

 mPtWgt        = kFALSE ; // pT as a weight
 mEtaWgt       = kFALSE ; // eta as a weight
 mOnePhiWgt    = kTRUE  ; // one/three phi-wgt histogram(s)

 mRedoWgt      = kFALSE ; // to recalculate Wgt histograms (even if flowPhi.hist.root exist)
 mMakeAll      = kTRUE ;  // claculates all events vars in one shoot (should run faster) 
}
//------------------------------------------------------------------------
// ###
//------------------------------------------------------------------------
Int_t AliFlowAnalysisMaker::Init()                    
{
 // Opens input, output and wgt files, books histograms, sets limits & bin 
 
 if(Debug0) { cout << FlowAnalysis << "Init() . " << endl ; }

 // Check for Read/Write weights ... if weight file is there, wgt arrays are filled (phiwgt & bayesian)
 WgtChk() ;  // -> FillWgtArrays(wgtFile)
 
 // Open output files (->plots)
 histFile = new TFile(GetHistFileName().Data(), "RECREATE");
 histFile->cd() ;

 // Variable setting
 float ptMaxPart = 5. ;      if(pFlowSelect->PtMaxPart())  { ptMaxPart = pFlowSelect->PtMaxPart() ; }
 int nPtBinsPart = 40  ;     if(pFlowSelect->PtBinsPart()) { nPtBinsPart = pFlowSelect->PtBinsPart() ; }
 xLabel = "Pseudorapidity" ; if(strlen(pFlowSelect->PidPart()) != 0) { xLabel = "Rapidity"; }
 // -			     
 const float etaMin  = Flow::fEtaMin ; 
 const float etaMax  = Flow::fEtaMax ; 
 const float ptMin   = Flow::fPtMin ;	     
 const float ptMax   = Flow::fPtMax ;	     
 const int   EtaBins  = Flow::nEtaBins ;      
 const int   nPtBins  = Flow::nPtBins ;      
 const int   nPhiBins = Flow::nPhiBins ;     
 // -
 SetHistoRanges(etaMin, etaMax, EtaBins);
 SetPtRange_for_vEta(1.,0.);   //(ptMin, ptMax);
 SetEtaRange_for_vPt(1.,0.);   //(etaMin, etaMax);
 // -
 const float triggerMin      =  -1.5;
 const float triggerMax      =  10.5;
 const float chargeMin       =  -2.5;
 const float chargeMax       =   2.5; 
 const float dcaMin	     =   0.0;
 const float dcaMax	     =  10.0; //  0.5; 
 const float glDcaMax	     =  50.0; 
 const float chi2Min	     =    0.;
 const float chi2Max	     =  500.; 
 const float chi2normMax     =   50.; 
 const float chi2MaxC	     =   30.; 
 const float chi2MaxI	     =   60.; 
 const float fitPtsMin       =  -0.5;
 const float fitPtsMaxTpc    = 160.5;
 const float fitPtsMaxIts    =   6.5;
 const float fitPtsMaxTrd    = 130.5;
 const float fitPtsMaxTof    =   1.5;
 const float fitOverMaxMin   =    0.;
 const float fitOverMaxMax   =   1.1;
 const float multOverOrigMin =    0.;
 const float multOverOrigMax =    1.;
 const float vertexZMin      =  -10.;
 const float vertexZMax      =   10.; 
 const float vertexXYMin     =  -0.1;
 const float vertexXYMax     =   0.1; 
 const float etaSymMinPart   =   -1.;
 const float etaSymMaxPart   =    1.;
 const float etaSymMin       =   -2.;
 const float etaSymMax       =    2.;
 const float phiMin	     =    0.;
 const float phiMax	     = 2*TMath::Pi() ; 
 const float psiMin	     =    0.;
 const float psiMax	     = 2*TMath::Pi() ; 
 const float multMin	     =    0. ;
 const float multMax	     =25000. ;
 const float multV0	     = 5000. ;
 const float qMin	     =    0. ;
 const float qMax	     =   20. ;
 const float pidMin	     =    0. ;
 const float pidMax	     =    1. ;
 const float centMin	     =  -0.5 ;
 const float centMax	     =   9.5 ;
 const float logpMin	     =   -2.5 ;
 const float logpMax	     =    2.5 ;
 const float pMin	     =     0. ;
 const float pMax	     =   100. ;
 const float dEdxMax	     =  1000. ;
 const float dEdxMaxTPC	     =  1500. ;
 const float TOFmin	     =  10000. ;
 const float TOFmax	     =  50000. ;
 const float TRDmax	     =  2000. ;
 const float lgMin	     =  -1000. ;
 const float lgMax	     =  1000. ;
 const float lgMinV0	     =  0. ;
 const float lgMaxV0	     =  50. ;
 const float massMin	     =  -0.1 ;
 const float massMax	     =  2.1 ;
 const float ZDCpartMin      = -0.5 ;
 const float ZDCpartMax	     =  1000.5 ;
 const float ZDCeMin	     =  0. ;
 const float ZDCeMax         =  10000. ;
 // - 
 enum { nTriggerBins	  = 12,
	nChargeBins	  = 5,
	nDcaBins	  = 1000,
	nChi2Bins	  = 150,
	nFitPtsBinsTpc    = 161,
	nFitPtsBinsIts    = 7,
	nFitPtsBinsTrd    = 131,
	nFitPtsBinsTof    = 2,
	nFitOverMaxBins   = 55,
	nMultOverOrigBins =100,
	nVertexZBins	  = 200,
	nVertexXYBins	  = 1000,
	nEtaSymBins	  = 200,
	nPhi3DBins	  = nPhiBins/6 ,
	nPsiBins	  = 36,
	nMultBins	  = 250,
	nPidBins	  = 100,
        nCentBins	  = 10,
	nDedxBins	 = 1000,
	nTofBins	 = 1000,
	nTrdBins	 = 100,
	nMomenBins	 = 500,
	n_qBins 	 =  50,
	nLgBins 	 = 500,
	nMassBins	 = 2200,
 	nZDCpartBins	 = 1001,
 	nZDCeBins	 = 200
      } ;

// *temp*   
  int maxMaxLabel = maxLabel ;
  if(Debug1) { cout << spaces << " Creating Labels hist at the beginning ( max bin by hand = " << maxMaxLabel << " ) . " << endl ; }
  TString* labTitle = new TString("TracksLabel_full") ; 
  mLabHist = new TH2F(labTitle->Data(),labTitle->Data(),Flow::nSels+4,0.,Flow::nSels+4.,maxMaxLabel+1,0.,maxMaxLabel+1.) ;
  if(Flow::nSels==2) { mLabHist->SetXTitle("tracks , sel0 , sel1 , selPart , v0s , selV0s ") ; }
  else  { mLabHist->SetXTitle("tracks , sel... , selPart , v0s , selV0s ") ; }
  mLabHist->SetYTitle("MC label") ;
  delete labTitle ;
  mLabellingFirst = kFALSE ; 
// *temp* 

 // Histograms Booking ...
 // Trigger
 mHistTrigger = new TH1F("Flow_Trigger", "Flow_Trigger", nTriggerBins, triggerMin, triggerMax);
 mHistTrigger->SetXTitle("Trigger: 1 minbias, 2 central, 3 laser, 10 other");
 mHistTrigger->SetYTitle("Counts");

 // Charge
 mHistCharge = new TH1F("Flow_Charge", "Flow_Charge", nChargeBins, chargeMin, chargeMax);
 mHistCharge->SetXTitle("Charge");
 mHistCharge->SetYTitle("Counts");

 // Reconstructed Momentum (constrained, if so)
 mHistPtot = new TH1F("Flow_totalP","Flow_totalP", nPtBins, ptMin, ptMax);
 mHistPtot->Sumw2();
 mHistPtot->SetXTitle("P (GeV/c)");
 mHistPtot->SetYTitle("Counts");

 // Reconstructed pT (constrained, if so)
 mHistPt = new TH1F("Flow_Pt","Flow_Pt", nMomenBins, pMin, pMax);
 mHistPt->Sumw2();
 mHistPt->SetXTitle("Pt (GeV/c)");
 mHistPt->SetYTitle("Counts");
   
 // Reconstructed P.id vs pT 
 mHistPidPt = new TH2F("Flow_PidPt","Flow_PidPt", 12, 0., 12., nPtBins, ptMin, ptMax);
 mHistPidPt->Sumw2();
 mHistPidPt->SetXTitle("e-, e+, mu-, mu+, pi-, pi+, K-, K+, pr-, pr+, d-, d+");
 mHistPidPt->SetYTitle("Pt (GeV/c)");
   
 // Distance of closest approach - Transverse, Unsigned
 mHistDca = new TH1F("Flow_TransDCA", "Flow_TransverseDCA", nDcaBins, dcaMin, dcaMax);
 mHistDca->Sumw2();
 mHistDca->SetXTitle("|Track's Signed dca (cm)|");
 mHistDca->SetYTitle("Counts");

 // Distance of closest approach - Transverse
 mHistTransDca = new TH1F("Flow_TransDCA_Signed", "Flow_TransverseDCA_Signed", nDcaBins, -dcaMax, dcaMax);
 mHistTransDca->Sumw2();
 mHistTransDca->SetXTitle("Track's Transverse dca (cm)");
 mHistTransDca->SetYTitle("Counts");

 // Distance of closest approach - 3d
 mHistDcaGlobal = new TH1F("Flow_3dDcaGlobal", "Flow_3dDcaGlobal", nDcaBins, dcaMin, glDcaMax);
 mHistDcaGlobal->Sumw2();
 mHistDcaGlobal->SetXTitle("|3d Global Track's dca to Vertex (cm)|");
 mHistDcaGlobal->SetYTitle("Counts");

 // Distance of closest approach for particles correlated with the event plane - 3d
 mHistDcaGlobalPart = new TH1F("Flow_3dDcaGlobalPart", "Flow_3dDcaGlobalPart", nDcaBins, dcaMin, glDcaMax);
 mHistDcaGlobalPart->Sumw2();
 mHistDcaGlobalPart->SetXTitle("|3d Global Track's dca to Vertex (cm)|");
 mHistDcaGlobalPart->SetYTitle("Counts");

 // Distance of closest approach for particles NOT SELECTED for correlation with the event plane - 3d
 mHistDcaGlobalOut = new TH1F("Flow_3dDcaGlobalOut", "Flow_3dDcaGlobalOut", nDcaBins, dcaMin, glDcaMax);
 mHistDcaGlobalOut->Sumw2();
 mHistDcaGlobalOut->SetXTitle("|3d Global Track's dca to Vertex (cm)|");
 mHistDcaGlobalOut->SetYTitle("Counts");

 // Yield Pt-Phi for constrained tracks
 mHistPhiPtCon = new TH2D("Flow_PtPhi_Con", "Flow_PtPhi_Con", nPhiBins, phiMin, phiMax, nPtBins, ptMin, ptMax);
 mHistPhiPtCon->Sumw2();
 mHistPhiPtCon->SetXTitle("Phi");
 mHistPhiPtCon->SetYTitle("Pt (GeV/c)");
 
 // Yield Pt-Phi for UNconstrained tracks
 mHistPhiPtUnc = new TH2D("Flow_PtPhi_Unc", "Flow_PtPhi_Unc", nPhiBins, phiMin, phiMax, nPtBins, ptMin, ptMax);
 mHistPhiPtUnc->Sumw2();
 mHistPhiPtUnc->SetXTitle("Phi");
 mHistPhiPtUnc->SetYTitle("Pt (GeV/c)");
 
 // Chi2 
 // at the main vertex
 mHistChi2 = new TH1F("Flow_Chi2", "Flow_Chi2", nChi2Bins, chi2Min, chi2MaxC);
 mHistChi2->SetXTitle("Chi square at Main Vertex");
 mHistChi2->SetYTitle("Counts");
 // TPC
 mHistChi2TPC = new TH1F("Flow_Chi2_TPC", "Flow_Chi2_TPC", nChi2Bins, chi2Min, chi2Max);
 mHistChi2TPC->SetXTitle("Chi square for TPC");
 mHistChi2TPC->SetYTitle("Counts");
 // ITS
 mHistChi2ITS = new TH1F("Flow_Chi2_ITS", "Flow_Chi2_ITS", nChi2Bins, chi2Min, chi2MaxI);
 mHistChi2ITS->SetXTitle("Chi square for ITS");
 mHistChi2ITS->SetYTitle("Counts");
 // TRD
 mHistChi2TRD = new TH1F("Flow_Chi2_TRD", "Flow_Chi2_TRD", nChi2Bins, chi2Min, chi2Max);
 mHistChi2TRD->SetXTitle("Chi square for TRD");
 mHistChi2TRD->SetYTitle("Counts");
 // TOF
 mHistChi2TOF = new TH1F("Flow_Chi2_TOF", "Flow_Chi2_TOF", nChi2Bins, chi2Min, chi2Max);
 mHistChi2TOF->SetXTitle("Chi square for TOF");
 mHistChi2TOF->SetYTitle("Counts");

 // Chi2 normalized
 // TPC
 mHistChi2normTPC = new TH1F("Flow_Chi2norm_TPC", "Flow_Chi2norm_TPC", nChi2Bins, chi2Min, chi2normMax);
 mHistChi2normTPC->SetXTitle("Normalized #Chi^{2} for TPC");
 mHistChi2normTPC->SetYTitle("Counts");
 // ITS
 mHistChi2normITS = new TH1F("Flow_Chi2norm_ITS", "Flow_Chi2norm_ITS", nChi2Bins, chi2Min, chi2normMax);
 mHistChi2normITS->SetXTitle("Normalized #Chi^{2} for ITS");
 mHistChi2normITS->SetYTitle("Counts");
 // TRD
 mHistChi2normTRD = new TH1F("Flow_Chi2norm_TRD", "Flow_Chi2norm_TRD", nChi2Bins, chi2Min, chi2normMax);
 mHistChi2normTRD->SetXTitle("Normalized #Chi^{2} for TRD");
 mHistChi2normTRD->SetYTitle("Counts");
 // TOF
 mHistChi2normTOF = new TH1F("Flow_Chi2norm_TOF", "Flow_Chi2norm_TOF", nChi2Bins, chi2Min, chi2normMax);
 mHistChi2normTOF->SetXTitle("Normalized #Chi^{2} for TOF");
 mHistChi2normTOF->SetYTitle("Counts");
 
 // FitPts
 // TPC
 mHistFitPtsTPC = new TH1F("Flow_FitPts_TPC", "Flow_FitPts_TPC", nFitPtsBinsTpc, fitPtsMin, fitPtsMaxTpc);
 mHistFitPtsTPC->SetXTitle("Fit Points");
 mHistFitPtsTPC->SetYTitle("Counts");
 // ITS
 mHistFitPtsITS = new TH1F("Flow_HitPts_ITS", "Flow_HitPts_ITS", nFitPtsBinsIts, fitPtsMin, fitPtsMaxIts);
 mHistFitPtsITS->SetXTitle("Fit Points");
 mHistFitPtsITS->SetYTitle("Counts");
 // TRD
 mHistFitPtsTRD = new TH1F("Flow_FitPts_TRD", "Flow_FitPts_TRD", nFitPtsBinsTrd, fitPtsMin, fitPtsMaxTrd);
 mHistFitPtsTRD->SetXTitle("Fit Points");
 mHistFitPtsTRD->SetYTitle("Counts");
 // TOF
 mHistFitPtsTOF = new TH1F("Flow_HitPts_TOF", "Flow_HitPts_TOF", nFitPtsBinsTof, fitPtsMin, fitPtsMaxTof);
 mHistFitPtsTOF->SetXTitle("Fit Points");
 mHistFitPtsTOF->SetYTitle("Counts");

 // MaxPts 
 // TPC
 mHistMaxPtsTPC = new TH1F("Flow_MaxPts_TPC", "Flow_MaxPts_TPC", nFitPtsBinsTpc, fitPtsMin, fitPtsMaxTpc);
 mHistMaxPtsTPC->SetXTitle("Max Points");
 mHistMaxPtsTPC->SetYTitle("Counts");
 // ITS
 mHistMaxPtsITS = new TH1F("Flow_MaxPts_ITS", "Flow_MaxPts_ITS", nFitPtsBinsIts, fitPtsMin, fitPtsMaxIts);
 mHistMaxPtsITS->SetXTitle("Max Points");
 mHistMaxPtsITS->SetYTitle("Counts");
 // TRD
 mHistMaxPtsTRD = new TH1F("Flow_MaxPts_TRD", "Flow_MaxPts_TRD", nFitPtsBinsTrd, fitPtsMin, fitPtsMaxTrd);
 mHistMaxPtsTRD->SetXTitle("Max Points");
 mHistMaxPtsTRD->SetYTitle("Counts");
 // TOF
 mHistMaxPtsTOF = new TH1F("Flow_MaxPts_TOF", "Flow_MaxPts_TOF", nFitPtsBinsTof, fitPtsMin, fitPtsMaxTof);
 mHistMaxPtsTOF->SetXTitle("Max Points");
 mHistMaxPtsTOF->SetYTitle("Counts");

 // FitOverMax 
 // Tpc
 mHistFitOverMaxTPC = new TH1F("Flow_FitOverMax_TPC", "Flow_FitOverMax_TPC", nFitOverMaxBins, fitOverMaxMin, fitOverMaxMax);
 mHistFitOverMaxTPC->SetXTitle("(Fit Points - 1) / Max Points");
 mHistFitOverMaxTPC->SetYTitle("Counts");  
 // All
 mHistFitOverMax = new TH1F("Flow_FitOverMax", "Flow_FitOverMax", nFitOverMaxBins, fitOverMaxMin, fitOverMaxMax);
 mHistFitOverMax->SetXTitle("(Fit Points - 1) / Max Points");
 mHistFitOverMax->SetYTitle("Counts");  

 // lenght
 mHistLenght = new TH1F("Flow_TrackLenght", "Flow_TrackLenght", nLgBins, lgMin, lgMax);
 mHistLenght->SetXTitle("Lenght of the Track (cm)");
 mHistLenght->SetYTitle("Counts");

  // OrigMult
 mHistOrigMult = new TH1F("Flow_OrigMult", "Flow_OrigMult", nMultBins, multMin, multMax);
 mHistOrigMult->SetXTitle("Original Mult");
 mHistOrigMult->SetYTitle("Counts");

 // MultEta
 mHistMultEta = new TH1F("Flow_MultEta", "Flow_MultEta", nMultBins, multMin, multMax);
 mHistMultEta->SetXTitle("Mult for Centrality");
 mHistMultEta->SetYTitle("Counts");
   
 // Mult
 mHistMult = new TH1F("Flow_Mult", "Flow_Mult", nMultBins, multMin, multMax);
 mHistMult->SetXTitle("Mult");
 mHistMult->SetYTitle("Counts");

 // V0s multiplicity
 mHistV0Mult = new TH1F("FlowV0_Mult","FlowV0_Mult", nMultBins, multMin, multV0);
 mHistV0Mult->SetXTitle("V0s Multiplicity");
 mHistV0Mult->SetYTitle("Counts");
   
 // MultOverOrig
 mHistMultOverOrig = new TH1F("Flow_MultOverOrig", "Flow_MultOverOrig", nMultOverOrigBins, multOverOrigMin, multOverOrigMax);
 mHistMultOverOrig->SetXTitle("Mult / Orig. Mult");
 mHistMultOverOrig->SetYTitle("Counts");
   
 // Mult correlated with the event planes
 mHistMultPart = new TH1F("Flow_MultPart", "Flow_MultPart", 2*nMultBins, multMin, multMax);
 mHistMultPart->SetXTitle("Mult of Correlated Particles");
 mHistMultPart->SetYTitle("Counts");
   
 // Mult correlated with the event planes in 1 unit rapidity
 mHistMultPartUnit = new TH1F("Flow_MultPartUnit", "Flow_MultPartUnit", 2*nMultBins, multMin, multMax/2);
 mHistMultPartUnit->SetXTitle("Mult of Correlated Particles (-0.5 < eta < 0.5)");
 mHistMultPartUnit->SetYTitle("Counts");
   
 // Mult of V0s correlated with the event planes
 mHistV0MultPart = new TH1F("FlowV0_MultPart", "FlowV0_MultPart", nMultBins, multMin, multV0);
 mHistV0MultPart->SetXTitle("Mult of Correlated V0s");
 mHistV0MultPart->SetYTitle("Counts");
   
 // VertexZ
 mHistVertexZ = new TH1F("Flow_VertexZ", "Flow_VertexZ", nVertexZBins, vertexZMin, vertexZMax);
 mHistVertexZ->SetXTitle("Vertex Z (cm)");
 mHistVertexZ->SetYTitle("Counts");
   
 // VertexXY
 mHistVertexXY2D = new TH2F("Flow_VertexXY2D", "Flow_VertexXY2D", nVertexXYBins, vertexXYMin, vertexXYMax, nVertexXYBins, vertexXYMin, vertexXYMax);
 mHistVertexXY2D->SetXTitle("Vertex X (cm)");
 mHistVertexXY2D->SetYTitle("Vertex Y (cm)");
   
 // EtaSym vs. Vertex Z Tpc
 mHistEtaSymVerZ2D = new TH2F("Flow_EtaSymVerZ2D", "Flow_EtaSymVerZ2D", nVertexZBins, vertexZMin, vertexZMax, nEtaSymBins, etaSymMin, etaSymMax);
 mHistEtaSymVerZ2D->SetXTitle("Vertex Z (cm)");
 mHistEtaSymVerZ2D->SetYTitle("Eta Symmetry TPC");

 // EtaSym Tpc
 mHistEtaSym = new TH1F("Flow_EtaSym_TPC", "Flow_EtaSym_TPC", nEtaSymBins, etaSymMin, etaSymMax);
 mHistEtaSym->SetXTitle("Eta Symmetry Ratio TPC");
 mHistEtaSym->SetYTitle("Counts");
   
 // EtaSym vs. Vertex Z Tpc - correlation analysis
 mHistEtaSymVerZ2DPart = new TH2F("Flow_EtaSymVerZ2D_part", "Flow_EtaSymVerZ2D_part", nVertexZBins, vertexZMin, vertexZMax, nEtaSymBins, etaSymMinPart, etaSymMaxPart);
 mHistEtaSymVerZ2DPart->SetXTitle("Vertex Z (cm)");
 mHistEtaSymVerZ2DPart->SetYTitle("Eta Symmetry TPC");

 // EtaSym Tpc - correlation analysis
 mHistEtaSymPart = new TH1F("Flow_EtaSym_TPC_part", "Flow_EtaSym_TPC_part", nEtaSymBins, etaSymMinPart, etaSymMaxPart);
 mHistEtaSymPart->SetXTitle("Eta Symmetry Ratio TPC");
 mHistEtaSymPart->SetYTitle("Counts");
   
 // phi , whatever
 mHistPhi = new TH1F("Flow_xPhi", "Flow_xPhi (all)", nPhiBins, phiMin, phiMax);
 mHistPhi->SetXTitle("Phi (rad)");
 mHistPhi->SetYTitle("Counts");

 // phi constrained
 mHistPhiCons = new TH1F("Flow_cPhi", "Flow_cPhi", nPhiBins, phiMin, phiMax);
 mHistPhiCons->SetXTitle("cPhi (rad)");
 mHistPhiCons->SetYTitle("Counts");
   
 // EtaPtPhi , whatever
 mHistAllEtaPtPhi3D = new TH3F("Flow_EtaPtPhi3Dall", "Flow_EtaPtPhi3Dall (whatever)", mNEtaBins, mEtaMin, mEtaMax, nPtBins, ptMin, ptMax, nPhi3DBins, phiMin, phiMax);
 mHistAllEtaPtPhi3D->SetXTitle("Eta");
 mHistAllEtaPtPhi3D->SetYTitle("Pt (GeV/c)");
 mHistAllEtaPtPhi3D->SetZTitle("Phi (rad)");
   
 // Constrained EtaPtPhi
 mHistConsEtaPtPhi3D = new TH3F("Flow_consEtaPtPhi3D", "Flow_consEtaPtPhi3D (constrainable)", mNEtaBins, mEtaMin, mEtaMax, nPtBins, ptMin, ptMax, nPhi3DBins, phiMin, phiMax);
 mHistConsEtaPtPhi3D->SetXTitle("cEta");
 mHistConsEtaPtPhi3D->SetYTitle("cPt (GeV/c)");
 mHistConsEtaPtPhi3D->SetZTitle("cPhi (rad)");
   
 // Global EtaPtPhi
 mHistGlobEtaPtPhi3D = new TH3F("Flow_globEtaPtPhi3D", "Flow_globEtaPtPhi3D (constrainable)", mNEtaBins, mEtaMin, mEtaMax, nPtBins, ptMin, ptMax, nPhi3DBins, phiMin, phiMax);
 mHistGlobEtaPtPhi3D->SetXTitle("gEta");
 mHistGlobEtaPtPhi3D->SetYTitle("gPt (GeV/c)");
 mHistGlobEtaPtPhi3D->SetZTitle("gPhi (rad)");

 // UnConstrained EtaPtPhi
 mHistUncEtaPtPhi3D = new TH3F("Flow_uncEtaPtPhi3D", "Flow_uncEtaPtPhi3D (un-constrainable)", mNEtaBins, mEtaMin, mEtaMax, nPtBins, ptMin, ptMax, nPhi3DBins, phiMin, phiMax);
 mHistUncEtaPtPhi3D->SetXTitle("gEta");
 mHistUncEtaPtPhi3D->SetYTitle("gPt (GeV/c)");
 mHistUncEtaPtPhi3D->SetZTitle("gPhi (rad)");
   
 // EtaPtPhi for particles correlated with the event plane
 mHistEtaPtPhi3DPart = new TH3F("Flow_EtaPtPhi3Dpart", "Flow_EtaPtPhi3Dpart (selected part)", mNEtaBins, mEtaMin, mEtaMax, nPtBins, ptMin, ptMax, nPhi3DBins, phiMin, phiMax);
 mHistEtaPtPhi3DPart->SetXTitle("Eta");
 mHistEtaPtPhi3DPart->SetYTitle("Pt (GeV/c)");
 mHistEtaPtPhi3DPart->SetZTitle("Phi (rad)");
   
 // EtaPtPhi for particles NOT SELECTED for correlation with the event plane
 mHistEtaPtPhi3DOut = new TH3F("Flow_EtaPtPhi3Dout", "Flow_EtaPtPhi3Dout (NOT selected part)", mNEtaBins, mEtaMin, mEtaMax, nPtBins, ptMin, ptMax, nPhi3DBins, phiMin, phiMax);
 mHistEtaPtPhi3DOut->SetXTitle("Eta");
 mHistEtaPtPhi3DOut->SetYTitle("Pt (GeV/c)");
 mHistEtaPtPhi3DOut->SetZTitle("Phi (rad)");
   
 // Yield Pt-Phi for all positive
 mHistPtPhiPos = new TH2D("Flow_PtPhi_Plus", "Flow_PtPhi_Plus", nPhiBins, phiMin, phiMax, 50, ptMin, ptMax);
 mHistPtPhiPos->Sumw2();
 mHistPtPhiPos->SetXTitle("Phi");
 mHistPtPhiPos->SetYTitle("Pt (GeV/c)");

 // Yield Pt-Phi for all negative
 mHistPtPhiNeg = new TH2D("Flow_PtPhi_Minus", "Flow_PtPhi_Minus", nPhiBins, phiMin, phiMax, 50, ptMin, ptMax);
 mHistPtPhiNeg->Sumw2();
 mHistPtPhiNeg->SetXTitle("Phi");
 mHistPtPhiNeg->SetYTitle("Pt (GeV/c)");
   
 // Yield for all particles
 mHistYieldAll2D = new TH2D("Flow_YieldAll2D", "Flow_YieldAll2D", mNEtaBins, mEtaMin, mEtaMax, nPtBins, ptMin, ptMax);
 mHistYieldAll2D->Sumw2();
 mHistYieldAll2D->SetXTitle("Pseudorapidty");
 mHistYieldAll2D->SetYTitle("Pt (GeV/c)");

 // Yield for constrainable tracks
 mHistYieldCon2D = new TH2D("Flow_YieldCons2D", "Flow_YieldCons2D", mNEtaBins, mEtaMin, mEtaMax, nPtBins, ptMin, ptMax);
 mHistYieldCon2D->Sumw2();
 mHistYieldCon2D->SetXTitle("Pseudorapidty");
 mHistYieldCon2D->SetYTitle("Pt (GeV/c)");

 // Yield for un-constrainable tracks
 mHistYieldUnc2D = new TH2D("Flow_YieldUnc2D", "Flow_YieldUnc2D", mNEtaBins, mEtaMin, mEtaMax, nPtBins, ptMin, ptMax);
 mHistYieldUnc2D->Sumw2();
 mHistYieldUnc2D->SetXTitle("Pseudorapidty");
 mHistYieldUnc2D->SetYTitle("Pt (GeV/c)");

 // Yield for particles correlated with the event plane
 mHistYieldPart2D = new TH2D("Flow_YieldPart2D", "Flow_YieldPart2D (selected part)", mNEtaBins, mEtaMin, mEtaMax, nPtBinsPart, ptMin, ptMaxPart);
 mHistYieldPart2D->Sumw2();
 mHistYieldPart2D->SetXTitle((char*)xLabel.Data());
 mHistYieldPart2D->SetYTitle("Pt (GeV/c)");

 // Yield for particles NOT SELECTED for correlation with the event plane
 mHistYieldOut2D = new TH2D("Flow_YieldOut2D", "Flow_YieldOut2D (NOT selected part)", mNEtaBins, mEtaMin, mEtaMax, nPtBins, ptMin, ptMax);
 mHistYieldOut2D->Sumw2();
 mHistYieldOut2D->SetXTitle("Pseudorapidty");
 mHistYieldOut2D->SetYTitle("Pt (GeV/c)");
 
 // invariant Mass for all particles (from TOF)
 mHistInvMass = new TH1F("Flow_InvMass", "Flow_InvMass (tof)", nMassBins, massMin, massMax);
 mHistInvMass->SetXTitle("Invariant Mass (GeV)");
 mHistInvMass->SetYTitle("Counts");
 
 // invariant Mass for particles correlated with the event plane (from TOF)
 mHistInvMassPart = new TH1F("Flow_InvMassPart", "Flow_InvMassPart (tof)", nMassBins, massMin, massMax);
 mHistInvMassPart->SetXTitle("Invariant Mass (GeV)");
 mHistInvMassPart->SetYTitle("Counts");
 
 // invariant Mass for particles NOT SELECTED for correlation with the event plane (from TOF)
 mHistInvMassOut = new TH1F("Flow_InvMassOut", "Flow_InvMassOut (tof)", nMassBins, massMin, massMax);
 mHistInvMassOut->SetXTitle("Invariant Mass (GeV)");
 mHistInvMassOut->SetYTitle("Counts");

 // Mean Eta in each bin for particles correlated with the event plane
 mHistBinEta = new TProfile("Flow_Bin_Eta", "Flow_Bin_Eta_part (selected part)", mNEtaBins, mEtaMin, mEtaMax, mEtaMin, mEtaMax, "");
 mHistBinEta->SetXTitle((char*)xLabel.Data());
 mHistBinEta->SetYTitle((char*)xLabel.Data());
 
 // Mean Pt in each bin for particles correlated with the event plane
 mHistBinPt = new TProfile("Flow_Bin_Pt", "Flow_Bin_Pt_part (selected part)", nPtBinsPart, ptMin, ptMaxPart, ptMin, ptMaxPart, "");
 mHistBinPt->SetXTitle("Pt (GeV/c)");
 mHistBinPt->SetYTitle("<Pt> (GeV/c)");
 
 // cos(n*phiLab)
 mHistCosPhi = new TProfile("Flow_CosPhiLab", "Flow_CosPhiLab", Flow::nHars, 0.5, (float)(Flow::nHars) + 0.5, -100., 100., "");
 mHistCosPhi->SetXTitle("Harmonic");
 mHistCosPhi->SetYTitle("<cos(n*PhiLab)> (%)");
   
 // PID pi+
 mHistPidPiPlus = new TH1F("Flow_PidPiPlus", "Flow_PidPiPlus", nPidBins, pidMin, pidMax);
 mHistPidPiPlus->SetXTitle("ALICE P.Id.");
 mHistPidPiPlus->SetYTitle("Counts");
   
 // PID pi-
 mHistPidPiMinus = new TH1F("Flow_PidPiMinus", "Flow_PidPiMinus", nPidBins, pidMin, pidMax);
 mHistPidPiMinus->SetXTitle("ALICE P.Id.");
 mHistPidPiMinus->SetYTitle("Counts");
   
 // PID proton
 mHistPidProton = new TH1F("Flow_PidProton", "Flow_PidProton", nPidBins, pidMin, pidMax);
 mHistPidProton->SetXTitle("ALICE P.Id.");
 mHistPidProton->SetYTitle("Counts");

 // PID anti proton
 mHistPidAntiProton = new TH1F("Flow_PidAntiProton", "Flow_PidAntiProton", nPidBins, pidMin, pidMax);
 mHistPidAntiProton->SetXTitle("ALICE P.Id.");
 mHistPidAntiProton->SetYTitle("Counts");

 // PID Kplus
 mHistPidKplus = new TH1F("Flow_PidKplus", "Flow_PidKplus", nPidBins, pidMin, pidMax);
 mHistPidKplus->SetXTitle("ALICE P.Id.");
 mHistPidKplus->SetYTitle("Counts");

 // PID Kminus
 mHistPidKminus = new TH1F("Flow_PidKminus", "Flow_PidKminus", nPidBins, pidMin, pidMax);
 mHistPidKminus->SetXTitle("ALICE P.Id.");
 mHistPidKminus->SetYTitle("Counts");

 // PID deuteron
 mHistPidDeuteron = new TH1F("Flow_PidDeuteron", "Flow_PidDeuteron", nPidBins, pidMin, pidMax);
 mHistPidDeuteron->SetXTitle("ALICE P.Id.");
 mHistPidDeuteron->SetYTitle("Counts");

 // PID anti deuteron
 mHistPidAntiDeuteron = new TH1F("Flow_PidAntiDeuteron", "Flow_PidAntiDeuteron", nPidBins, pidMin, pidMax);
 mHistPidAntiDeuteron->SetXTitle("ALICE P.Id.");
 mHistPidAntiDeuteron->SetYTitle("Counts");

 // PID electron
 mHistPidElectron = new TH1F("Flow_PidElectron", "Flow_PidElectron", nPidBins, pidMin, pidMax);
 mHistPidElectron->SetXTitle("ALICE P.Id.");
 mHistPidElectron->SetYTitle("Counts");

 // PID positron
 mHistPidPositron = new TH1F("Flow_PidPositron", "Flow_PidPositron", nPidBins, pidMin, pidMax);
 mHistPidPositron->SetXTitle("ALICE P.Id.");
 mHistPidPositron->SetYTitle("Counts");

 // PID Muon+
 mHistPidMuonPlus = new TH1F("Flow_PidMuonPlus", "Flow_PidMuonPlus", nPidBins, pidMin, pidMax);
 mHistPidMuonPlus->SetXTitle("ALICE P.Id.");
 mHistPidMuonPlus->SetYTitle("Counts");

 // PID Muon-
 mHistPidMuonMinus = new TH1F("Flow_PidMuonMinus", "Flow_PidMuonMinus", nPidBins, pidMin, pidMax);
 mHistPidMuonMinus->SetXTitle("ALICE P.Id.");
 mHistPidMuonMinus->SetYTitle("Counts");

 // PID pi+ selected
 mHistPidPiPlusPart = new TH1F("Flow_PidPiPlusPart", "Flow_PidPiPlusPart", nPidBins, pidMin, pidMax);
 mHistPidPiPlusPart->SetXTitle("ALICE P.Id. (pid = pi+)");
 mHistPidPiPlusPart->SetYTitle("Counts");
   
 // PID pi- selected
 mHistPidPiMinusPart = new TH1F("Flow_PidPiMinusPart", "Flow_PidPiMinusPart", nPidBins, pidMin, pidMax);
 mHistPidPiMinusPart->SetXTitle("ALICE P.Id. (pid = pi-)");
 mHistPidPiMinusPart->SetYTitle("Counts");
   
 // PID proton selected
 mHistPidProtonPart = new TH1F("Flow_PidProtonPart", "Flow_PidProtonPart", nPidBins, pidMin, pidMax);
 mHistPidProtonPart->SetXTitle("ALICE P.Id. (pid = p)");
 mHistPidProtonPart->SetYTitle("Counts");

 // PID anti proton selected
 mHistPidAntiProtonPart = new TH1F("Flow_PidAntiProtonPart", "Flow_PidAntiProtonPart", nPidBins, pidMin, pidMax);
 mHistPidAntiProtonPart->SetXTitle("ALICE P.Id. (pid = p-)");
 mHistPidAntiProtonPart->SetYTitle("Counts");

 // PID Kplus selected
 mHistPidKplusPart = new TH1F("Flow_PidKplusPart", "Flow_PidKplusPart", nPidBins, pidMin, pidMax);
 mHistPidKplusPart->SetXTitle("ALICE P.Id. (pid = K+)");
 mHistPidKplusPart->SetYTitle("Counts");

 // PID Kminus selected
 mHistPidKminusPart = new TH1F("Flow_PidKminusPart", "Flow_PidKminusPart", nPidBins, pidMin, pidMax);
 mHistPidKminusPart->SetXTitle("ALICE P.Id. (pid = K-)");
 mHistPidKminusPart->SetYTitle("Counts");

 // PID deuteron selected
 mHistPidDeuteronPart = new TH1F("Flow_PidDeuteronPart", "Flow_PidDeuteronPart", nPidBins, pidMin, pidMax);
 mHistPidDeuteronPart->SetXTitle("ALICE P.Id. (pid = d)");
 mHistPidDeuteronPart->SetYTitle("Counts");

 // PID anti deuteron selected
 mHistPidAntiDeuteronPart = new TH1F("Flow_PidAntiDeuteronPart", "Flow_PidAntiDeuteronPart", nPidBins, pidMin, pidMax);
 mHistPidAntiDeuteronPart->SetXTitle("ALICE P.Id. (pid = d--)");
 mHistPidAntiDeuteronPart->SetYTitle("Counts");

 // PID electron selected
 mHistPidElectronPart = new TH1F("Flow_PidElectronPart", "Flow_PidElectronPart", nPidBins, pidMin, pidMax);
 mHistPidElectronPart->SetXTitle("ALICE P.Id. (pid = e-)");
 mHistPidElectronPart->SetYTitle("Counts");

 // PID positron selected
 mHistPidPositronPart = new TH1F("Flow_PidPositronPart", "Flow_PidPositronPart", nPidBins, pidMin, pidMax);
 mHistPidPositronPart->SetXTitle("ALICE P.Id. (pid = e+)");
 mHistPidPositronPart->SetYTitle("Counts");

 // PID Muon+ selected
 mHistPidMuonPlusPart = new TH1F("Flow_PidMuonPlusPart", "Flow_PidMuonPlusPart", nPidBins, pidMin, pidMax);
 mHistPidMuonPlusPart->SetXTitle("ALICE P.Id. (pid = mu+)");
 mHistPidMuonPlusPart->SetYTitle("Counts");

 // PID Muon- selected
 mHistPidMuonMinusPart = new TH1F("Flow_PidMuonMinusPart", "Flow_PidMuonMinusPart", nPidBins, pidMin, pidMax);
 mHistPidMuonMinusPart->SetXTitle("ALICE P.Id. (pid = mu-)");
 mHistPidMuonMinusPart->SetYTitle("Counts");

 // PID multiplicities (all)
 mHistPidMult = new TProfile("Flow_PidMult", "Flow_PidMult", 15, 0.5, 15.5, "");
 mHistPidMult->SetXTitle("all, h+, h-, pi+, pi-, pr+, pr-, K+, K-, d+, d-, e-, e+, mu-, mu+");
 mHistPidMult->SetYTitle("Multiplicity");

 // PID for all tracks
 mHistBayPidMult = new TH1F("Flow_BayPidMult","Flow_BayPidMult",Flow::nPid,-0.5,((float)Flow::nPid-0.5));
 mHistBayPidMult->Sumw2() ;
 mHistBayPidMult->SetXTitle("e+/- , mu+/- , pi+/- , K+/- , p+/- , d+/- ");
 mHistBayPidMult->SetYTitle("Counts");
   
 // PID for particles correlated with the event plane
 mHistBayPidMultPart = new TH1F("Flow_BayPidMult_Part","Flow_BayPidMult_Part",Flow::nPid,-0.5,((float)Flow::nPid-0.5));
 mHistBayPidMultPart->Sumw2() ;
 mHistBayPidMultPart->SetXTitle("e+/- , mu+/- , pi+/- , K+/- , p+/- , d+/- ");
 mHistBayPidMultPart->SetYTitle("Counts");

 // Centrality
 mHistCent = new TH1F("Flow_Cent", "Flow_Cent", nCentBins, centMin, centMax);
 mHistCent->SetXTitle("Centrality Bin");
 mHistCent->SetYTitle("Counts");
   
 // Deposited Energy in ZDC
 mHistEnergyZDC = new TH2F("Flow_ZDC_E", "Flow_ZDC_E", 3, 0., 3., nZDCeBins, ZDCeMin, ZDCeMax);
 mHistEnergyZDC->SetXTitle("neutron , proton , e.m.");
 mHistEnergyZDC->SetYTitle("ZDC energy");

 // Part. versus Energy in ZDC
 mHistPartZDC = new TH1F("Flow_ZDC_Participants", "Flow_ZDC_Participants", nZDCpartBins, ZDCpartMin, ZDCpartMax);
 mHistPartZDC->SetXTitle("ZDC part");
 mHistPartZDC->SetYTitle("Counts");

 // MeanDedxPos TPC
 mHistMeanDedxPos2D = new TH2F("Flow_MeanDedxPos2D_TPC","Flow_MeanDedxPos2D_TPC", nMomenBins, logpMin, logpMax, nDedxBins, 0, dEdxMaxTPC);
 mHistMeanDedxPos2D->SetXTitle("log(momentum) (GeV/c)");
 mHistMeanDedxPos2D->SetYTitle("mean dEdx (+)");
 
 // MeanDedxNeg TPC
 mHistMeanDedxNeg2D = new TH2F("Flow_MeanDedxNeg2D_TPC","Flow_MeanDedxNeg2D_TPC", nMomenBins, logpMin, logpMax, nDedxBins, 0, dEdxMaxTPC);
 mHistMeanDedxNeg2D->SetXTitle("log(momentum) (GeV/c)");
 mHistMeanDedxNeg2D->SetYTitle("mean dEdx (-)");

 // MeanDedxPos ITS
 mHistMeanDedxPos2DITS = new TH2F("Flow_MeanDedxPos2D_ITS","Flow_MeanDedxPos2D_ITS", nMomenBins, logpMin, logpMax, nDedxBins, 0, dEdxMax);
 mHistMeanDedxPos2DITS->SetXTitle("log(momentum) (GeV/c)");
 mHistMeanDedxPos2DITS->SetYTitle("mean dEdx (+)");
 
 // MeanDedxNeg ITS
 mHistMeanDedxNeg2DITS = new TH2F("Flow_MeanDedxNeg2D_ITS","Flow_MeanDedxNeg2D_ITS", nMomenBins, logpMin, logpMax, nDedxBins, 0, dEdxMax);
 mHistMeanDedxNeg2DITS->SetXTitle("log(momentum) (GeV/c)");
 mHistMeanDedxNeg2DITS->SetYTitle("mean dEdx (-)");

 // MeanDedxPos TPC 3D selected Part
 mHistMeanDedxPos3DPart = new TH3F("Flow_MeanDedxPos3DPart_TPC","Flow_MeanDedxPos3DPart_TPC", nMomenBins, logpMin, logpMax, nDedxBins, 0., dEdxMaxTPC, Flow::nPid, 0., Flow::nPid);
 mHistMeanDedxPos3DPart->SetXTitle("log(momentum) (GeV/c)");
 mHistMeanDedxPos3DPart->SetYTitle("mean dEdx (+)");
 mHistMeanDedxPos3DPart->SetZTitle("e , mu , pi , k , p , d");
 
 // MeanDedxNeg TPC 3D selected Part
 mHistMeanDedxNeg3DPart = new TH3F("Flow_MeanDedxNeg3DPart_TPC","Flow_MeanDedxNeg3DPart_TPC", nMomenBins, logpMin, logpMax, nDedxBins, 0., dEdxMaxTPC, Flow::nPid, 0., Flow::nPid);
 mHistMeanDedxNeg3DPart->SetXTitle("log(momentum) (GeV/c)");
 mHistMeanDedxNeg3DPart->SetYTitle("mean dEdx (-)");
 mHistMeanDedxNeg3DPart->SetZTitle("e , mu , pi , k , p , d");

 // MeanDedxPos ITS 3D selected Part
 mHistMeanDedxPos3DPartITS = new TH3F("Flow_MeanDedxPos3DPart_ITS","Flow_MeanDedxPos3DPart_ITS", nMomenBins, logpMin, logpMax, nDedxBins, 0., dEdxMax, Flow::nPid, 0., Flow::nPid);
 mHistMeanDedxPos3DPartITS->SetXTitle("log(momentum) (GeV/c)");
 mHistMeanDedxPos3DPartITS->SetYTitle("mean dEdx (+)");
 mHistMeanDedxPos3DPartITS->SetZTitle("e , mu , pi , k , p , d");
 
 // MeanDedxNeg ITS 3D selected Part
 mHistMeanDedxNeg3DPartITS = new TH3F("Flow_MeanDedxNeg3DPart_ITS","Flow_MeanDedxNeg3DPart_ITS", nMomenBins, logpMin, logpMax, nDedxBins, 0., dEdxMax, Flow::nPid, 0., Flow::nPid);
 mHistMeanDedxNeg3DPartITS->SetXTitle("log(momentum) (GeV/c)");
 mHistMeanDedxNeg3DPartITS->SetYTitle("mean dEdx (-)");
 mHistMeanDedxNeg3DPartITS->SetZTitle("e , mu , pi , k , p , d");

 // MeanSignalPos TRD
 mHistMeanDedxPos2DTRD = new TH2F("Flow_MeanSignalPos2D_TRD","Flow_MeanSignalPos2D_TRD", nMomenBins, logpMin, logpMax, nTrdBins, 0, TRDmax);
 mHistMeanDedxPos2DTRD->SetXTitle("log(momentum) (GeV/c)");
 mHistMeanDedxPos2DTRD->SetYTitle("mean TRD (+)");
 
 // MeanSignalNeg TRD
 mHistMeanDedxNeg2DTRD = new TH2F("Flow_MeanSignalNeg2D_TRD","Flow_MeanSignalNeg2D_TRD", nMomenBins, logpMin, logpMax, nTrdBins, 0, TRDmax);
 mHistMeanDedxNeg2DTRD->SetXTitle("log(momentum) (GeV/c)");
 mHistMeanDedxNeg2DTRD->SetYTitle("mean TRD (-)");

 // MeanTimePos TOF
 mHistMeanDedxPos2DTOF = new TH2F("Flow_MeanTimePos2D_TOF","Flow_MeanTimePos2D_TOF", nMomenBins, logpMin, logpMax, nTofBins, TOFmin, TOFmax);
 mHistMeanDedxPos2DTOF->SetXTitle("log(momentum) (GeV/c)");
 mHistMeanDedxPos2DTOF->SetYTitle("mean Time of Flight (+)");
 
 // MeanTimeNeg TOF
 mHistMeanDedxNeg2DTOF = new TH2F("Flow_MeanTimeNeg2D_TOF","Flow_MeanTimeNeg2D_TOF", nMomenBins, logpMin, logpMax, nTofBins, TOFmin, TOFmax);
 mHistMeanDedxNeg2DTOF->SetXTitle("log(momentum) (GeV/c)");
 mHistMeanDedxNeg2DTOF->SetYTitle("mean Time of Flight (-)");
 
 // Detector response for each particle { ... 
 // MeanDedx PiPlus - TPC
 mHistMeanTPCPiPlus = new TH2F("Flow_MeanDedxPiPlusTPC", "Flow_MeanDedxPiPlusTPC", nMomenBins, logpMin, logpMax, nDedxBins, 0, dEdxMaxTPC);
 mHistMeanTPCPiPlus->SetXTitle("Log10(p) (GeV/c)");
 mHistMeanTPCPiPlus->SetYTitle("mean dEdx (pi+)");
 // MeanDedx PiMinus
 mHistMeanTPCPiMinus = new TH2F("Flow_MeanDedxPiMinusTPC", "Flow_MeanDedxPiMinusTPC", nMomenBins, logpMin, logpMax, nDedxBins, 0, dEdxMaxTPC);
 mHistMeanTPCPiMinus->SetXTitle("Log10(p) (GeV/c)");
 mHistMeanTPCPiMinus->SetYTitle("mean dEdx (pi-)");
 // MeanDedx Proton
 mHistMeanTPCProton = new TH2F("Flow_MeanDedxProtonTPC", "Flow_MeanDedxProtonTPC", nMomenBins, logpMin, logpMax, nDedxBins, 0, dEdxMaxTPC);
 mHistMeanTPCProton->SetXTitle("Log10(p) (GeV/c)");
 mHistMeanTPCProton->SetYTitle("mean dEdx (pr+)");
 // MeanDedx Pbar
 mHistMeanTPCPbar = new TH2F("Flow_MeanDedxPbarTPC", "Flow_MeanDedxPbarTPC", nMomenBins, logpMin, logpMax, nDedxBins, 0, dEdxMaxTPC);
 mHistMeanTPCPbar->SetXTitle("Log10(p) (GeV/c)");
 mHistMeanTPCPbar->SetYTitle("mean dEdx (pr-)");
 // MeanDedx Kplus
 mHistMeanTPCKplus = new TH2F("Flow_MeanDedxKplusTPC", "Flow_MeanDedxKplusTPC", nMomenBins, logpMin, logpMax, nDedxBins, 0, dEdxMaxTPC);
 mHistMeanTPCKplus->SetXTitle("Log10(p) (GeV/c)");
 mHistMeanTPCKplus->SetYTitle("mean dEdx (K+)");
 // MeanDedx Kminus
 mHistMeanTPCKminus = new TH2F("Flow_MeanDedxKminusTPC", "Flow_MeanDedxKminusTPC", nMomenBins, logpMin, logpMax, nDedxBins, 0, dEdxMaxTPC);
 mHistMeanTPCKminus->SetXTitle("Log10(p) (GeV/c)");
 mHistMeanTPCKminus->SetYTitle("mean dEdx (K-)");
 // MeanDedx Deuteron
 mHistMeanTPCDeuteron = new TH2F("Flow_MeanDedxDeuteronTPC", "Flow_MeanDedxDeuteronTPC", nMomenBins, logpMin, logpMax, nDedxBins, 0, dEdxMaxTPC);
 mHistMeanTPCDeuteron->SetXTitle("Log10(p) (GeV/c)");
 mHistMeanTPCDeuteron->SetYTitle("mean dEdx (d+)");
 // MeanDedx AntiDeuteron
 mHistMeanTPCAntiDeuteron = new TH2F("Flow_MeanDedxAntiDeuteronTPC", "Flow_MeanDedxAntiDeuteronTPC", nMomenBins, logpMin, logpMax, nDedxBins, 0, dEdxMaxTPC);
 mHistMeanTPCAntiDeuteron->SetXTitle("Log10(p) (GeV/c)");
 mHistMeanTPCAntiDeuteron->SetYTitle("mean dEdx (d-)");
 // MeanDedxElectron
 mHistMeanTPCElectron = new TH2F("Flow_MeanDedxElectronTPC", "Flow_MeanDedxElectronTPC", nMomenBins, logpMin, logpMax, nDedxBins, 0, dEdxMaxTPC);
 mHistMeanTPCElectron->SetXTitle("Log10(p) (GeV/c)");
 mHistMeanTPCElectron->SetYTitle("mean dEdx (e-)");
 // MeanDedx Positron
 mHistMeanTPCPositron = new TH2F("Flow_MeanDedxPositronTPC", "Flow_MeanDedxPositronTPC", nMomenBins, logpMin, logpMax, nDedxBins, 0, dEdxMaxTPC);
 mHistMeanTPCPositron->SetXTitle("Log10(p) (GeV/c)");
 mHistMeanTPCPositron->SetYTitle("mean dEdx (e+)");
 // MeanDedx MuonPlus
 mHistMeanTPCMuonPlus = new TH2F("Flow_MeanDedxMuonPlusTPC", "Flow_MeanDedxMuonPlusTPC", nMomenBins, logpMin, logpMax, nDedxBins, 0, dEdxMaxTPC);
 mHistMeanTPCMuonPlus->SetXTitle("Log10(p) (GeV/c)");
 mHistMeanTPCMuonPlus->SetYTitle("mean dEdx (mu+)");
 // MeanDedx MuonMinus
 mHistMeanTPCMuonMinus = new TH2F("Flow_MeanDedxMuonMinusTPC", "Flow_MeanDedxMuonMinusTPC", nMomenBins, logpMin, logpMax, nDedxBins, 0, dEdxMaxTPC);
 mHistMeanTPCMuonMinus->SetXTitle("Log10(p) (GeV/c)");
 mHistMeanTPCMuonMinus->SetYTitle("mean dEdx (mu-)");
 
 // MeanDedx PiPlus - ITS
 mHistMeanITSPiPlus = new TH2F("Flow_MeanDedxPiPlusITS", "Flow_MeanDedxPiPlusITS", nMomenBins, logpMin, logpMax, nDedxBins, 0, dEdxMax);
 mHistMeanITSPiPlus->SetXTitle("Log10(p) (GeV/c)");
 mHistMeanITSPiPlus->SetYTitle("mean dEdx (pi+)");
 // MeanDedx PiMinus
 mHistMeanITSPiMinus = new TH2F("Flow_MeanDedxPiMinusITS", "Flow_MeanDedxPiMinusITS", nMomenBins, logpMin, logpMax, nDedxBins, 0, dEdxMax);
 mHistMeanITSPiMinus->SetXTitle("Log10(p) (GeV/c)");
 mHistMeanITSPiMinus->SetYTitle("mean dEdx (pi-)");
 // MeanDedx Proton
 mHistMeanITSProton = new TH2F("Flow_MeanDedxProtonITS", "Flow_MeanDedxProtonITS", nMomenBins, logpMin, logpMax, nDedxBins, 0, dEdxMax);
 mHistMeanITSProton->SetXTitle("Log10(p) (GeV/c)");
 mHistMeanITSProton->SetYTitle("mean dEdx (pr+)");
 // MeanDedx Pbar
 mHistMeanITSPbar = new TH2F("Flow_MeanDedxPbarITS", "Flow_MeanDedxPbarITS", nMomenBins, logpMin, logpMax, nDedxBins, 0, dEdxMax);
 mHistMeanITSPbar->SetXTitle("Log10(p) (GeV/c)");
 mHistMeanITSPbar->SetYTitle("mean dEdx (pr-)");
 // MeanDedx Kplus
 mHistMeanITSKplus = new TH2F("Flow_MeanDedxKplusITS", "Flow_MeanDedxKplusITS", nMomenBins, logpMin, logpMax, nDedxBins, 0, dEdxMax);
 mHistMeanITSKplus->SetXTitle("Log10(p) (GeV/c)");
 mHistMeanITSKplus->SetYTitle("mean dEdx (K+)");
 // MeanDedx Kminus
 mHistMeanITSKminus = new TH2F("Flow_MeanDedxKminusITS", "Flow_MeanDedxKminusITS", nMomenBins, logpMin, logpMax, nDedxBins, 0, dEdxMax);
 mHistMeanITSKminus->SetXTitle("Log10(p) (GeV/c)");
 mHistMeanITSKminus->SetYTitle("mean dEdx (K-)");
 // MeanDedx Deuteron
 mHistMeanITSDeuteron = new TH2F("Flow_MeanDedxDeuteronITS", "Flow_MeanDedxDeuteronITS", nMomenBins, logpMin, logpMax, nDedxBins, 0, dEdxMax);
 mHistMeanITSDeuteron->SetXTitle("Log10(p) (GeV/c)");
 mHistMeanITSDeuteron->SetYTitle("mean dEdx (d+)");
 // MeanDedx AntiDeuteron
 mHistMeanITSAntiDeuteron = new TH2F("Flow_MeanDedxAntiDeuteronITS", "Flow_MeanDedxAntiDeuteronITS", nMomenBins, logpMin, logpMax, nDedxBins, 0, dEdxMax);
 mHistMeanITSAntiDeuteron->SetXTitle("Log10(p) (GeV/c)");
 mHistMeanITSAntiDeuteron->SetYTitle("mean dEdx (d-)");
 // MeanDedx Electron
 mHistMeanITSElectron = new TH2F("Flow_MeanDedxElectronITS", "Flow_MeanDedxElectronITS", nMomenBins, logpMin, logpMax, nDedxBins, 0, dEdxMax);
 mHistMeanITSElectron->SetXTitle("Log10(p) (GeV/c)");
 mHistMeanITSElectron->SetYTitle("mean dEdx (e-)");
 // MeanDedx Positron
 mHistMeanITSPositron = new TH2F("Flow_MeanDedxPositronITS", "Flow_MeanDedxPositronITS", nMomenBins, logpMin, logpMax, nDedxBins, 0, dEdxMax);
 mHistMeanITSPositron->SetXTitle("Log10(p) (GeV/c)");
 mHistMeanITSPositron->SetYTitle("mean dEdx (e+)");
 // MeanDedx MuonPlus
 mHistMeanITSMuonPlus = new TH2F("Flow_MeanDedxMuonPlusITS", "Flow_MeanDedxMuonPlusITS", nMomenBins, logpMin, logpMax, nDedxBins, 0, dEdxMax);
 mHistMeanITSMuonPlus->SetXTitle("Log10(p) (GeV/c)");
 mHistMeanITSMuonPlus->SetYTitle("mean dEdx (mu+)");
 // MeanDedx MuonMinus
 mHistMeanITSMuonMinus = new TH2F("Flow_MeanDedxMuonMinusITS", "Flow_MeanDedxMuonMinusITS", nMomenBins, logpMin, logpMax, nDedxBins, 0, dEdxMax);
 mHistMeanITSMuonMinus->SetXTitle("Log10(p) (GeV/c)");
 mHistMeanITSMuonMinus->SetYTitle("mean dEdx (mu-)");
 
 // MeanTrd PiPlus - TRD
 mHistMeanTRDPiPlus = new TH2F("Flow_MeanTrdPiPlusTRD", "Flow_MeanTrdPiPlusTRD", nMomenBins, logpMin, logpMax, nTrdBins, 0, TRDmax);
 mHistMeanTRDPiPlus->SetXTitle("Log10(p) (GeV/c)");
 mHistMeanTRDPiPlus->SetYTitle("mean t.r.[] (pi+)");
 // MeanTrd PiMinus
 mHistMeanTRDPiMinus = new TH2F("Flow_MeanTrdPiMinusTRD", "Flow_MeanTrdPiMinusTRD", nMomenBins, logpMin, logpMax, nTrdBins, 0, TRDmax);
 mHistMeanTRDPiMinus->SetXTitle("Log10(p) (GeV/c)");
 mHistMeanTRDPiMinus->SetYTitle("mean t.r.[] (pi-)");
 // MeanTrd Proton
 mHistMeanTRDProton = new TH2F("Flow_MeanTrdProtonTRD", "Flow_MeanTrdProtonTRD", nMomenBins, logpMin, logpMax, nTrdBins, 0, TRDmax);
 mHistMeanTRDProton->SetXTitle("Log10(p) (GeV/c)");
 mHistMeanTRDProton->SetYTitle("mean t.r.[] (pr+)");
 // MeanTrd Pbar
 mHistMeanTRDPbar = new TH2F("Flow_MeanTrdPbarTRD", "Flow_MeanTrdPbarTRD", nMomenBins, logpMin, logpMax, nTrdBins, 0, TRDmax);
 mHistMeanTRDPbar->SetXTitle("Log10(p) (GeV/c)");
 mHistMeanTRDPbar->SetYTitle("mean t.r.[] (pr-)");
 // MeanTrd Kplus
 mHistMeanTRDKplus = new TH2F("Flow_MeanTrdKplusTRD", "Flow_MeanTrdKplusTRD", nMomenBins, logpMin, logpMax, nTrdBins, 0, TRDmax);
 mHistMeanTRDKplus->SetXTitle("Log10(p) (GeV/c)");
 mHistMeanTRDKplus->SetYTitle("mean t.r.[] (K+)");
 // MeanTrd Kminus
 mHistMeanTRDKminus = new TH2F("Flow_MeanTrdKminusTRD", "Flow_MeanTrdKminusTRD", nMomenBins, logpMin, logpMax, nTrdBins, 0, TRDmax);
 mHistMeanTRDKminus->SetXTitle("Log10(p) (GeV/c)");
 mHistMeanTRDKminus->SetYTitle("mean t.r.[] (K-)");
 // MeanTrd Deuteron
 mHistMeanTRDDeuteron = new TH2F("Flow_MeanTrdDeuteronTRD", "Flow_MeanTrdDeuteronTRD", nMomenBins, logpMin, logpMax, nTrdBins, 0, TRDmax);
 mHistMeanTRDDeuteron->SetXTitle("Log10(p) (GeV/c)");
 mHistMeanTRDDeuteron->SetYTitle("mean t.r.[] (d+)");
 // MeanTrd AntiDeuteron
 mHistMeanTRDAntiDeuteron = new TH2F("Flow_MeanTrdAntiDeuteronTRD", "Flow_MeanTrdAntiDeuteronTRD", nMomenBins, logpMin, logpMax, nTrdBins, 0, TRDmax);
 mHistMeanTRDAntiDeuteron->SetXTitle("Log10(p) (GeV/c)");
 mHistMeanTRDAntiDeuteron->SetYTitle("mean t.r.[] (d-)");
 // MeanTrd Electron
 mHistMeanTRDElectron = new TH2F("Flow_MeanTrdElectronTRD", "Flow_MeanTrdElectronTRD", nMomenBins, logpMin, logpMax, nTrdBins, 0, TRDmax);
 mHistMeanTRDElectron->SetXTitle("Log10(p) (GeV/c)");
 mHistMeanTRDElectron->SetYTitle("mean t.r.[] (e-)");
 // MeanTrd Positron
 mHistMeanTRDPositron = new TH2F("Flow_MeanTrdPositronTRD", "Flow_MeanTrdPositronTRD", nMomenBins, logpMin, logpMax, nTrdBins, 0, TRDmax);
 mHistMeanTRDPositron->SetXTitle("Log10(p) (GeV/c)");
 mHistMeanTRDPositron->SetYTitle("mean t.r.[] (e+)");
 // MeanTrd MuonPlus
 mHistMeanTRDMuonPlus = new TH2F("Flow_MeanTrdMuonPlusTRD", "Flow_MeanTrdMuonPlusTRD", nMomenBins, logpMin, logpMax, nTrdBins, 0, TRDmax);
 mHistMeanTRDMuonPlus->SetXTitle("Log10(p) (GeV/c)");
 mHistMeanTRDMuonPlus->SetYTitle("mean t.r.[] (mu+)");
 // MeanTrd MuonMinus
 mHistMeanTRDMuonMinus = new TH2F("Flow_MeanTrdMuonMinusTRD", "Flow_MeanTrdMuonMinusTRD", nMomenBins, logpMin, logpMax, nTrdBins, 0, TRDmax);
 mHistMeanTRDMuonMinus->SetXTitle("Log10(p) (GeV/c)");
 mHistMeanTRDMuonMinus->SetYTitle("mean t.r.[] (mu-)");
 
 // T.O.F. PiPlus - TOF
 mHistMeanTOFPiPlus = new TH2F("Flow_MeanTofPiPlusTOF", "Flow_MeanTofPiPlusTOF", nMassBins, massMin, massMax, nTofBins, TOFmin, TOFmax);
 mHistMeanTOFPiPlus->SetXTitle("invariant mass (GeV)");
 mHistMeanTOFPiPlus->SetYTitle("mean t.o.f.[psec] (pi+)");
 // MeanTof PiMinus
 mHistMeanTOFPiMinus = new TH2F("Flow_MeanTofPiMinusTOF", "Flow_MeanTofPiMinusTOF", nMassBins, massMin, massMax, nTofBins, TOFmin, TOFmax);
 mHistMeanTOFPiMinus->SetXTitle("invariant mass (GeV)");
 mHistMeanTOFPiMinus->SetYTitle("mean t.o.f.[psec] (pi-)");
 // MeanTof Proton
 mHistMeanTOFProton = new TH2F("Flow_MeanTofProtonTOF", "Flow_MeanTofProtonTOF", nMassBins, massMin, massMax, nTofBins, TOFmin, TOFmax);
 mHistMeanTOFProton->SetXTitle("invariant mass (GeV)");
 mHistMeanTOFProton->SetYTitle("mean t.o.f.[psec] (pr+)");
 // Mean TofPbar
 mHistMeanTOFPbar = new TH2F("Flow_MeanTofPbarTOF", "Flow_MeanTofPbarTOF", nMassBins, massMin, massMax, nTofBins, TOFmin, TOFmax);
 mHistMeanTOFPbar->SetXTitle("invariant mass (GeV)");
 mHistMeanTOFPbar->SetYTitle("mean t.o.f.[psec] (pr-)");
 // mean t.o.f.[psec]Kplus
 mHistMeanTOFKplus = new TH2F("Flow_MeanTofKplusTOF", "Flow_MeanTofKplusTOF", nMassBins, massMin, massMax, nTofBins, TOFmin, TOFmax);
 mHistMeanTOFKplus->SetXTitle("invariant mass (GeV)");
 mHistMeanTOFKplus->SetYTitle("mean t.o.f.[psec] (K+)");
 // mean t.o.f.[psec]Kminus
 mHistMeanTOFKminus = new TH2F("Flow_MeanTofKminusTOF", "Flow_MeanTofKminusTOF", nMassBins, massMin, massMax, nTofBins, TOFmin, TOFmax);
 mHistMeanTOFKminus->SetXTitle("invariant mass (GeV)");
 mHistMeanTOFKminus->SetYTitle("mean t.o.f.[psec] (K-)");
 // MeanTof Deuteron
 mHistMeanTOFDeuteron = new TH2F("Flow_MeanTofDeuteronTOF", "Flow_MeanTofDeuteronTOF", nMassBins, massMin, massMax, nTofBins, TOFmin, TOFmax);
 mHistMeanTOFDeuteron->SetXTitle("invariant mass (GeV)");
 mHistMeanTOFDeuteron->SetYTitle("mean t.o.f.[psec] (d+)");
 // MeanTof AntiDeuteron
 mHistMeanTOFAntiDeuteron = new TH2F("Flow_MeanTofAntiDeuteronTOF", "Flow_MeanTofAntiDeuteronTOF", nMassBins, massMin, massMax, nTofBins, TOFmin, TOFmax);
 mHistMeanTOFAntiDeuteron->SetXTitle("invariant mass (GeV)");
 mHistMeanTOFAntiDeuteron->SetYTitle("mean t.o.f.[psec] (d-)");
 // MeanTof Electron
 mHistMeanTOFElectron = new TH2F("Flow_MeanTofElectronTOF", "Flow_MeanTofElectronTOF", nMassBins, massMin, massMax, nTofBins, TOFmin, TOFmax);
 mHistMeanTOFElectron->SetXTitle("invariant mass (GeV)");
 mHistMeanTOFElectron->SetYTitle("mean t.o.f.[psec] (e-)");
 // MeanTof Positron
 mHistMeanTOFPositron = new TH2F("Flow_MeanTofPositronTOF", "Flow_MeanTofPositronTOF", nMassBins, massMin, massMax, nTofBins, TOFmin, TOFmax);
 mHistMeanTOFPositron->SetXTitle("invariant mass (GeV)");
 mHistMeanTOFPositron->SetYTitle("mean t.o.f.[psec] (e+)");
 // MeanTof MuonPlus
 mHistMeanTOFMuonPlus = new TH2F("Flow_MeanTofMuonPlusTOF", "Flow_MeanTofMuonPlusTOF", nMassBins, massMin, massMax, nTofBins, TOFmin, TOFmax);
 mHistMeanTOFMuonPlus->SetXTitle("invariant mass (GeV)");
 mHistMeanTOFMuonPlus->SetYTitle("mean t.o.f.[psec] (mu+)");
 // MeanTof MuonMinus
 mHistMeanTOFMuonMinus = new TH2F("Flow_MeanTofMuonMinusTOF", "Flow_MeanTofMuonMinusTOF", nMassBins, massMin, massMax, nTofBins, TOFmin, TOFmax);
 mHistMeanTOFMuonMinus->SetXTitle("invariant mass (GeV)");
 mHistMeanTOFMuonMinus->SetYTitle("mean t.o.f.[psec] (mu-)");
 // ... }

 TString* histTitle;
 for (int n = 0; n < Flow::nSubs; n++) // for sub-events
 {
  for (int k = 0; k < Flow::nSels; k++) 
  {
   for (int j = 0; j < Flow::nHars; j++) 
   {
    float order = (float)(j + 1);
    int i = Flow::nSubs * k + n ;

    // event planes
    histTitle = new TString("Flow_Psi_Sub");
    *histTitle += n+1;
    histTitle->Append("_Sel");
    *histTitle += k+1;
    histTitle->Append("_Har");
    *histTitle += j+1;
    histSub[i].histSubHar[j].mHistPsiSubs = new TH1F(histTitle->Data(),histTitle->Data(), nPsiBins, psiMin, (psiMax / order));
    histSub[i].histSubHar[j].mHistPsiSubs->SetXTitle("Event Plane Angle (rad)");
    histSub[i].histSubHar[j].mHistPsiSubs->SetYTitle("Counts");
    delete histTitle;
   }
  }
 }
 
 if(mV0)            // All V0s (if there, if flag on)
 {
 // Mass
  mHistV0Mass = new TH1F("FlowV0_InvMass", "FlowV0_InvMass", nMassBins, massMin, massMax);
  mHistV0Mass->SetXTitle("Invariant Mass (GeV)");
  mHistV0Mass->SetYTitle("Counts");
 // Distance of closest approach
  mHistV0Dca = new TH1F("FlowV0_Dca", "FlowV0_Dca", nDcaBins, dcaMin, glDcaMax);
  mHistV0Dca->SetXTitle("dca between tracks (cm)");
  mHistV0Dca->SetYTitle("Counts");
 // lenght
  mHistV0Lenght = new TH1F("FlowV0_Lenght", "FlowV0_Lenght", nLgBins, lgMinV0, lgMaxV0);
  mHistV0Lenght->SetXTitle("Distance of V0s (cm)");
  mHistV0Lenght->SetYTitle("Counts");
 // Sigma for all particles
  mHistV0Sigma = new TH1F("FlowV0_Sigma", "FlowV0_Sigma", nLgBins, lgMinV0, lgMaxV0 );
  mHistV0Sigma->SetXTitle("Sigma");
  mHistV0Sigma->SetYTitle("Counts");
 // Chi2 
  mHistV0Chi2 = new TH1F("FlowV0_Chi2", "FlowV0_Chi2", nChi2Bins, chi2Min, chi2MaxC);
  mHistV0Chi2->SetXTitle("Chi square at Main Vertex");
  mHistV0Chi2->SetYTitle("Counts");
 // EtaPtPhi
  mHistV0EtaPtPhi3D = new TH3F("FlowV0_EtaPtPhi3D", "FlowV0_EtaPtPhi3D", mNEtaBins, mEtaMin, mEtaMax, nPtBins, ptMin, ptMax, nPhi3DBins, phiMin, phiMax);
  mHistV0EtaPtPhi3D->SetXTitle("Eta");
  mHistV0EtaPtPhi3D->SetYTitle("Pt (GeV/c)");
  mHistV0EtaPtPhi3D->SetZTitle("Phi (rad)");
 // Yield for all v0s
  mHistV0YieldAll2D = new TH2D("FlowV0_YieldAll2D", "FlowV0_YieldAll2D", mNEtaBins, mEtaMin, mEtaMax, nPtBins, ptMin, ptMax);
  mHistV0YieldAll2D->Sumw2();
  mHistV0YieldAll2D->SetXTitle("Pseudorapidty");
  mHistV0YieldAll2D->SetYTitle("Pt (GeV/c)");
 // Mass slices on pT
  mHistV0MassPtSlices = new TH2D("FlowV0_MassPtSlices", "FlowV0_MassPtSlices", nMassBins, massMin, massMax, nPtBins, ptMin, ptMax);
  mHistV0MassPtSlices->Sumw2();
  mHistV0MassPtSlices->SetXTitle("Invariant Mass (GeV)");
  mHistV0MassPtSlices->SetYTitle("Pt (GeV/c)");

 // Selected V0s ...  
 // Yield 
  mHistV0YieldPart2D = new TH2D("FlowV0_Yield2Dsel", "FlowV0_Yield2Dsel", mNEtaBins, mEtaMin, mEtaMax, nPtBins, ptMin, ptMax);
  mHistV0YieldPart2D->Sumw2();
  mHistV0YieldPart2D->SetXTitle("Pseudorapidty");
  mHistV0YieldPart2D->SetYTitle("Pt (GeV/c)");
 // Mass Window
  mHistV0MassWin = new TH1F("FlowV0_MassWinPart", "FlowV0_MassWinPart", nMassBins, massMin, massMax);
  mHistV0MassWin->SetXTitle("Invariant Mass (GeV)");
  mHistV0MassWin->SetYTitle("Counts");
 // EtaPtPhi
  mHistV0EtaPtPhi3DPart = new TH3F("FlowV0_EtaPtPhi3Dpart", "FlowV0_EtaPtPhi3Dpart", mNEtaBins, mEtaMin, mEtaMax, nPtBins, ptMin, ptMax, nPhi3DBins, phiMin, phiMax);
  mHistV0EtaPtPhi3DPart->SetXTitle("Eta");
  mHistV0EtaPtPhi3DPart->SetYTitle("Pt (GeV/c)");
  mHistV0EtaPtPhi3DPart->SetZTitle("Phi (rad)");
 // Distance of closest approach
  mHistV0DcaPart = new TH1F("FlowV0_DcaPart", "FlowV0_DcaPart", nDcaBins, dcaMin, dcaMax);
  mHistV0DcaPart->Sumw2();
  mHistV0DcaPart->SetXTitle("dca between tracks (cm)");
  mHistV0DcaPart->SetYTitle("Counts");
 // lenght
  mHistV0LenghtPart = new TH1F("FlowV0_LenghtPart", "FlowV0_LenghtPart", nLgBins, lgMinV0, lgMaxV0);
  mHistV0LenghtPart->SetXTitle("Distance of V0s (cm)");
  mHistV0LenghtPart->SetYTitle("Counts");
 // SideBand Mass (sidebands)
  mHistV0sbMassSide = new TH1F("FlowV0sb_MassWinSideBands", "FlowV0sb_MassWinSideBands", nMassBins, massMin, massMax);
  mHistV0sbMassSide->SetXTitle("Invariant Mass (GeV)");
  mHistV0sbMassSide->SetYTitle("Counts");
 // EtaPtPhi (sidebands)
  mHistV0sbEtaPtPhi3DPart = new TH3F("FlowV0sb_EtaPtPhi3D", "FlowV0sb_EtaPtPhi3D", mNEtaBins, mEtaMin, mEtaMax, nPtBins, ptMin, ptMax, nPhi3DBins, phiMin, phiMax);
  mHistV0sbEtaPtPhi3DPart->SetXTitle("Eta");
  mHistV0sbEtaPtPhi3DPart->SetYTitle("Pt (GeV/c)");
  mHistV0sbEtaPtPhi3DPart->SetZTitle("Phi (rad)");
 // Distance of closest approach (sidebands)
  mHistV0sbDcaPart = new TH1F("FlowV0sb_Dca", "FlowV0sb_Dca", nDcaBins, dcaMin, dcaMax);
  mHistV0sbDcaPart->Sumw2();
  mHistV0sbDcaPart->SetXTitle("dca between tracks (cm)");
  mHistV0sbDcaPart->SetYTitle("Counts");
 // lenght (sidebands)
  mHistV0sbLenghtPart = new TH1F("FlowV0sb_Lenght", "FlowV0sb_Lenght", nLgBins, lgMinV0, lgMaxV0);
  mHistV0sbLenghtPart->SetXTitle("Distance of V0s (cm)");
  mHistV0sbLenghtPart->SetYTitle("Counts");

 // Mean Eta in each bin
  mHistV0BinEta = new TProfile("FlowV0_Bin_Eta", "FlowV0_Bin_Eta", mNEtaBins, mEtaMin, mEtaMax, mEtaMin, mEtaMax, "");
  mHistV0BinEta->SetXTitle((char*)xLabel.Data());
  mHistV0BinEta->SetYTitle("<Eta>");
 // Mean Pt in each bin
  mHistV0BinPt = new TProfile("FlowV0_Bin_Pt", "FlowV0_Bin_Pt", nPtBinsPart, ptMin, ptMaxPart, ptMin, ptMaxPart, "");
  mHistV0BinPt->SetXTitle("Pt (GeV/c)");
  mHistV0BinPt->SetYTitle("<Pt> (GeV/c)");
 // Mean Eta in each bin (sidebands)
  mHistV0sbBinEta = new TProfile("FlowV0sb_Bin_Eta", "FlowV0sb_Bin_Eta", mNEtaBins, mEtaMin, mEtaMax, mEtaMin, mEtaMax, "");
  mHistV0sbBinEta->SetXTitle((char*)xLabel.Data());
  mHistV0sbBinEta->SetYTitle("<Eta>");
 // Mean Pt in each bin (sidebands)
  mHistV0sbBinPt = new TProfile("FlowV0sb_Bin_Pt", "FlowV0sb_Bin_Pt", nPtBinsPart, ptMin, ptMaxPart, ptMin, ptMaxPart, "");
  mHistV0sbBinPt->SetXTitle("Pt (GeV/c)");
  mHistV0sbBinPt->SetYTitle("<Pt> (GeV/c)");
 }

 for (int k = 0; k < Flow::nSels; k++) // for each selection
 {
  // cos(n*delta_Psi)
  histTitle = new TString("Flow_Cos_Sel");
  *histTitle += k+1;
  histFull[k].mHistCos = new TProfile(histTitle->Data(), histTitle->Data(), Flow::nHars, 0.5, (float)(Flow::nHars) + 0.5, -1., 1., "");
  histFull[k].mHistCos->SetXTitle("Harmonic");
  histFull[k].mHistCos->SetYTitle("<cos(n*delta_Psi)>");
  delete histTitle;
   
  // resolution
  histTitle = new TString("Flow_Res_Sel");
  *histTitle += k+1;
  histFull[k].mHistRes = new TH1F(histTitle->Data(), histTitle->Data(), Flow::nHars, 0.5, (float)(Flow::nHars) + 0.5);
  histFull[k].mHistRes->SetXTitle("Harmonic");
  histFull[k].mHistRes->SetYTitle("Resolution");
  delete histTitle;

  // vObs
  histTitle = new TString("Flow_vObs_Sel");
  *histTitle += k+1;
  histFull[k].mHist_vObs = new TProfile(histTitle->Data(), histTitle->Data(), Flow::nHars, 0.5, (float)(Flow::nHars) + 0.5, -100., 100., "");
  histFull[k].mHist_vObs->SetXTitle("Harmonic");
  histFull[k].mHist_vObs->SetYTitle("vObs (%)");
  delete histTitle;

  // vObs V0
  histTitle = new TString("FlowV0_vObs_Sel");
  *histTitle += k+1;
  histFull[k].mHistV0_vObs = new TProfile(histTitle->Data(), histTitle->Data(), Flow::nHars, 0.5, (float)(Flow::nHars) + 0.5, -100., 100., "");
  histFull[k].mHistV0_vObs->SetXTitle("Harmonic");
  histFull[k].mHistV0_vObs->SetYTitle("vObs (%)");
  delete histTitle;
   
  // vObs V0 sideband SX
  histTitle = new TString("FlowV0sb_vObs_sx_Sel");
  *histTitle += k+1;
  histFull[k].mHistV0sb_vObs_sx = new TProfile(histTitle->Data(), histTitle->Data(), Flow::nHars, 0.5, (float)(Flow::nHars) + 0.5, -100., 100., "");
  histFull[k].mHistV0sb_vObs_sx->SetXTitle("Harmonic");
  histFull[k].mHistV0sb_vObs_sx->SetYTitle("vObs (%)");
  delete histTitle;
   
  // vObs V0 sideband DX
  histTitle = new TString("FlowV0sb_vObs_dx_Sel");
  *histTitle += k+1;
  histFull[k].mHistV0sb_vObs_dx = new TProfile(histTitle->Data(), histTitle->Data(), Flow::nHars, 0.5, (float)(Flow::nHars) + 0.5, -100., 100., "");
  histFull[k].mHistV0sb_vObs_dx->SetXTitle("Harmonic");
  histFull[k].mHistV0sb_vObs_dx->SetYTitle("vObs (%)");
  delete histTitle;
   
  // PID for tracks used in R.P.
  histTitle = new TString("Flow_BayPidMult_Sel");
  *histTitle += k+1;
  histFull[k].mHistBayPidMult = new TH1F(histTitle->Data(), histTitle->Data(),Flow::nPid,-0.5,((float)Flow::nPid-0.5));
  histFull[k].mHistBayPidMult->Sumw2() ;
  histFull[k].mHistBayPidMult->SetXTitle("e+/-  ,  mu+/-  ,  pi+/-  ,  K+/-  ,  p+/-  ,  d+/- ");
  histFull[k].mHistBayPidMult->SetYTitle("Counts");
  delete histTitle;

  for (int j = 0; j < Flow::nHars; j++)   // for each harmonic
  {
   float order  = (float)(j+1);

   // multiplicity
   histTitle = new TString("Flow_Mul_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   histFull[k].histFullHar[j].mHistMult = new TH1F(histTitle->Data(),histTitle->Data(), nMultBins, multMin, multMax);
   histFull[k].histFullHar[j].mHistMult->SetXTitle("Multiplicity");
   histFull[k].histFullHar[j].mHistMult->SetYTitle("Counts");
   delete histTitle;

   // event plane
   histTitle = new TString("Flow_Psi_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   histFull[k].histFullHar[j].mHistPsi = new TH1F(histTitle->Data(), histTitle->Data(), nPsiBins, psiMin, psiMax / order);
   histFull[k].histFullHar[j].mHistPsi->SetXTitle("Event Plane Angle (rad)");
   histFull[k].histFullHar[j].mHistPsi->SetYTitle("Counts");
   delete histTitle;
     
   // event plane difference of two selections
   histTitle = new TString("Flow_Psi_Diff_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   if (k == 0 ) 
   {
    Int_t my_order = 1;
    if (j == 1) { my_order = 2 ; }
    histFull[k].histFullHar[j].mHistPsi_Diff = new TH1F(histTitle->Data(), histTitle->Data(), nPsiBins, -psiMax/my_order/2., psiMax/my_order/2.);
   } 
   else 
   {
    histFull[k].histFullHar[j].mHistPsi_Diff = new TH1F(histTitle->Data(), histTitle->Data(), nPsiBins, -psiMax/2., psiMax/2.);
   }
   if (k == 0) 
   {
    if (j == 0) 
    {
     histFull[k].histFullHar[j].mHistPsi_Diff->SetXTitle("#Psi_{1,Sel1} - #Psi_{1,Sel2}(rad)");
    } 
    else if (j == 1) 
    {
     histFull[k].histFullHar[j].mHistPsi_Diff->SetXTitle("#Psi_{2,Sel1} - #Psi_{2,Sel2}(rad)");
    }
   } 
   else if (k == 1) 
   {
    if (j == 0) 
    {  
     histFull[k].histFullHar[j].mHistPsi_Diff->SetXTitle("#Psi_{1,Sel1} - #Psi_{2,Sel2}(rad)");
    } 
    else if (j == 1) 
    {
     histFull[k].histFullHar[j].mHistPsi_Diff->SetXTitle("#Psi_{1,Sel1} - #Psi_{2,Sel1}(rad)");
    }
   }
   histFull[k].histFullHar[j].mHistPsi_Diff->SetYTitle("Counts");
   delete histTitle;

   // correlation of sub-event planes
   histTitle = new TString("Flow_Psi_Sub_Corr_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   histFull[k].histFullHar[j].mHistPsiSubCorr = new TH1F(histTitle->Data(), histTitle->Data(), nPsiBins, psiMin, psiMax / order);
   histFull[k].histFullHar[j].mHistPsiSubCorr->Sumw2();
   histFull[k].histFullHar[j].mHistPsiSubCorr->SetXTitle("Sub-Event Correlation (rad)");
   histFull[k].histFullHar[j].mHistPsiSubCorr->SetYTitle("Counts");
   delete histTitle;
     
   // correlation of sub-event planes of different order
   histTitle = new TString("Flow_Psi_Sub_Corr_Diff_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   histFull[k].histFullHar[j].mHistPsiSubCorrDiff = new TH1F(histTitle->Data(), histTitle->Data(), nPsiBins, psiMin, psiMax / (order+1.));
   histFull[k].histFullHar[j].mHistPsiSubCorrDiff->Sumw2();
   histFull[k].histFullHar[j].mHistPsiSubCorrDiff->SetXTitle("Sub-Event Correlation (rad)");
   histFull[k].histFullHar[j].mHistPsiSubCorrDiff->SetYTitle("Counts");
   delete histTitle;
     
   // q
   histTitle = new TString("Flow_NormQ_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   histFull[k].histFullHar[j].mHist_q = new TH1F(histTitle->Data(), histTitle->Data(), n_qBins, qMin, qMax);
   histFull[k].histFullHar[j].mHist_q->Sumw2();
   histFull[k].histFullHar[j].mHist_q->SetXTitle("q = |Q|/sqrt(Mult)");
   histFull[k].histFullHar[j].mHist_q->SetYTitle("Counts");
   delete histTitle;

    // particle-plane azimuthal correlation
   histTitle = new TString("Flow_Phi_Corr_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   histFull[k].histFullHar[j].mHistPhiCorr = new TH1F(histTitle->Data(), histTitle->Data(), nPhiBins, phiMin, phiMax / order);
   histFull[k].histFullHar[j].mHistPhiCorr->Sumw2();
   histFull[k].histFullHar[j].mHistPhiCorr->SetXTitle("Particle-Plane Correlation (rad)");
   histFull[k].histFullHar[j].mHistPhiCorr->SetYTitle("Counts");
   delete histTitle;
     
   // neutral particle-plane azimuthal correlation
   histTitle = new TString("FlowV0_Phi_Corr_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   histFull[k].histFullHar[j].mHistV0PhiCorr = new TH1F(histTitle->Data(), histTitle->Data(), nPhiBins, phiMin, phiMax / order);
   histFull[k].histFullHar[j].mHistV0PhiCorr->Sumw2();
   histFull[k].histFullHar[j].mHistV0PhiCorr->SetXTitle("V0-Plane Correlation (rad)");
   histFull[k].histFullHar[j].mHistV0PhiCorr->SetYTitle("Counts");
   delete histTitle;

   // neutral sidebands-plane azimuthal correlation
   histTitle = new TString("FlowV0sb_Phi_Corr_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   histFull[k].histFullHar[j].mHistV0sbPhiCorr = new TH1F(histTitle->Data(), histTitle->Data(), nPhiBins, phiMin, phiMax / order);
   histFull[k].histFullHar[j].mHistV0sbPhiCorr->Sumw2();
   histFull[k].histFullHar[j].mHistV0sbPhiCorr->SetXTitle("V0sideBands-Plane Correlation (rad)");
   histFull[k].histFullHar[j].mHistV0sbPhiCorr->SetYTitle("Counts");
   delete histTitle;

   // Yield(pt)
   histTitle = new TString("Flow_YieldPt_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   histFull[k].histFullHar[j].mHistYieldPt = new TH1F(histTitle->Data(), histTitle->Data(), nPtBins, ptMin, ptMax);
   histFull[k].histFullHar[j].mHistYieldPt->Sumw2();
   histFull[k].histFullHar[j].mHistYieldPt->SetXTitle("Pt (GeV/c)");
   histFull[k].histFullHar[j].mHistYieldPt->SetYTitle("Yield");
   delete histTitle;
     
   // EtaPtPhi 
   histTitle = new TString("Flow_EtaPtPhi3D_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   histFull[k].histFullHar[j].mHistEtaPtPhi3D = new TH3F(histTitle->Data(), histTitle->Data(), mNEtaBins, mEtaMin, mEtaMax, nPtBins, ptMin, ptMax, nPhi3DBins, phiMin, phiMax);
   histFull[k].histFullHar[j].mHistEtaPtPhi3D->SetXTitle("Eta");
   histFull[k].histFullHar[j].mHistEtaPtPhi3D->SetYTitle("Pt (GeV/c)");
   histFull[k].histFullHar[j].mHistEtaPtPhi3D->SetZTitle("Phi (rad)");
   delete histTitle;

   // Yield(eta,pt)
   histTitle = new TString("Flow_Yield2D_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   histFull[k].histFullHar[j].mHistYield2D = new TH2D(histTitle->Data(), histTitle->Data(), mNEtaBins, mEtaMin, mEtaMax, nPtBins, ptMin, ptMax);
   histFull[k].histFullHar[j].mHistYield2D->Sumw2();
   histFull[k].histFullHar[j].mHistYield2D->SetXTitle("Pseudorapidty");
   histFull[k].histFullHar[j].mHistYield2D->SetYTitle("Pt (GeV/c)");
   delete histTitle;

   // Dca - 3D
   histTitle = new TString("Flow_3dDca_Sel") ;
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   histFull[k].histFullHar[j].mHistDcaGlob = new TH1F(histTitle->Data(), histTitle->Data(), nDcaBins, dcaMin, glDcaMax);
   histFull[k].histFullHar[j].mHistDcaGlob->Sumw2();
   histFull[k].histFullHar[j].mHistDcaGlob->SetXTitle("|3d Global Track's dca to Vertex (cm)|");
   delete histTitle;

   // Yield(pt) - excluded from R.P.
   histTitle = new TString("Flow_YieldPt_outSel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   histFull[k].histFullHar[j].mHistYieldPtout = new TH1F(histTitle->Data(), histTitle->Data(), nPtBins, ptMin, ptMax);
   histFull[k].histFullHar[j].mHistYieldPtout->Sumw2();
   histFull[k].histFullHar[j].mHistYieldPtout->SetXTitle("Pt (GeV/c)");
   histFull[k].histFullHar[j].mHistYieldPtout->SetYTitle("Yield");
   delete histTitle;
     
   // EtaPtPhi - excluded from R.P. 
   histTitle = new TString("Flow_EtaPtPhi3D_outSel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   histFull[k].histFullHar[j].mHistEtaPtPhi3Dout = new TH3F(histTitle->Data(), histTitle->Data(), mNEtaBins, mEtaMin, mEtaMax, nPtBins, ptMin, ptMax, nPhi3DBins, phiMin, phiMax);
   histFull[k].histFullHar[j].mHistEtaPtPhi3Dout->SetXTitle("Eta");
   histFull[k].histFullHar[j].mHistEtaPtPhi3Dout->SetYTitle("Pt (GeV/c)");
   histFull[k].histFullHar[j].mHistEtaPtPhi3Dout->SetZTitle("Phi (rad)");
   delete histTitle;

   // Yield(eta,pt) - excluded from R.P.
   histTitle = new TString("Flow_Yield2D_outSel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   histFull[k].histFullHar[j].mHistYield2Dout = new TH2D(histTitle->Data(), histTitle->Data(), mNEtaBins, mEtaMin, mEtaMax, nPtBins, ptMin, ptMax);
   histFull[k].histFullHar[j].mHistYield2Dout->Sumw2();
   histFull[k].histFullHar[j].mHistYield2Dout->SetXTitle("Pseudorapidty");
   histFull[k].histFullHar[j].mHistYield2Dout->SetYTitle("Pt (GeV/c)");
   delete histTitle;

   // Dca - 3D - excluded from R.P.
   histTitle = new TString("Flow_3dDca_outSel") ;
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   histFull[k].histFullHar[j].mHistDcaGlobout = new TH1F(histTitle->Data(), histTitle->Data(), nDcaBins, dcaMin, glDcaMax);
   histFull[k].histFullHar[j].mHistDcaGlobout->Sumw2();
   histFull[k].histFullHar[j].mHistDcaGlobout->SetXTitle("|3d Global Track's dca to Vertex (cm)|");
   delete histTitle;

   // Flow observed - v_obs,Pt,Eta
   histTitle = new TString("Flow_vObs2D_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   histFull[k].histFullHar[j].mHist_vObs2D = new TProfile2D(histTitle->Data(), histTitle->Data(), mNEtaBins, mEtaMin, mEtaMax, nPtBinsPart, ptMin, ptMaxPart, -100., 100., "");
   histFull[k].histFullHar[j].mHist_vObs2D->SetXTitle((char*)xLabel.Data());
   histFull[k].histFullHar[j].mHist_vObs2D->SetYTitle("Pt (GeV/c)");
   delete histTitle;

   // v_obs,Eta
   histTitle = new TString("Flow_vObsEta_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   histFull[k].histFullHar[j].mHist_vObsEta = new TProfile(histTitle->Data(), histTitle->Data(), mNEtaBins, mEtaMin, mEtaMax, -100., 100., "");
   histFull[k].histFullHar[j].mHist_vObsEta->SetXTitle((char*)xLabel.Data());
   histFull[k].histFullHar[j].mHist_vObsEta->SetYTitle("v (%)");
   delete histTitle;

   // v_obs,Pt
   histTitle = new TString("Flow_vObsPt_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   histFull[k].histFullHar[j].mHist_vObsPt = new TProfile(histTitle->Data(), histTitle->Data(), nPtBinsPart, ptMin, ptMaxPart, -100., 100., "");
   histFull[k].histFullHar[j].mHist_vObsPt->SetXTitle("Pt (GeV/c)");
   histFull[k].histFullHar[j].mHist_vObsPt->SetYTitle("v (%)");
   delete histTitle;

   // neutral Flow observed - Pt,Eta
   histTitle = new TString("FlowV0_vObs2D_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   histFull[k].histFullHar[j].mHistV0_vObs2D = new TProfile2D(histTitle->Data(), histTitle->Data(), mNEtaBins, mEtaMin, mEtaMax, nPtBinsPart, ptMin, ptMaxPart, -100., 100., "");
   histFull[k].histFullHar[j].mHistV0_vObs2D->SetXTitle((char*)xLabel.Data());
   histFull[k].histFullHar[j].mHistV0_vObs2D->SetYTitle("Pt (GeV/c)");
   delete histTitle;

   // neutral Flow observed - Eta
   histTitle = new TString("FlowV0_vObsEta_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   histFull[k].histFullHar[j].mHistV0_vObsEta = new TProfile(histTitle->Data(), histTitle->Data(), mNEtaBins, mEtaMin, mEtaMax, -100., 100., "");
   histFull[k].histFullHar[j].mHistV0_vObsEta->SetXTitle((char*)xLabel.Data());
   histFull[k].histFullHar[j].mHistV0_vObsEta->SetYTitle("v (%)");
   delete histTitle;

   // neutral Flow observed - Pt
   histTitle = new TString("FlowV0_vObsPt_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   histFull[k].histFullHar[j].mHistV0_vObsPt = new TProfile(histTitle->Data(), histTitle->Data(), nPtBinsPart, ptMin, ptMaxPart, -100., 100., "");
   histFull[k].histFullHar[j].mHistV0_vObsPt->SetXTitle("Pt (GeV/c)");
   histFull[k].histFullHar[j].mHistV0_vObsPt->SetYTitle("v (%)");
   delete histTitle;

   // neutral sidebands Flow observed - Pt,Eta
   histTitle = new TString("FlowV0sb_vObs2D_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   histFull[k].histFullHar[j].mHistV0sb_vObs2D = new TProfile2D(histTitle->Data(), histTitle->Data(), mNEtaBins, mEtaMin, mEtaMax, nPtBinsPart, ptMin, ptMaxPart, -100., 100., "");
   histFull[k].histFullHar[j].mHistV0sb_vObs2D->SetXTitle((char*)xLabel.Data());
   histFull[k].histFullHar[j].mHistV0sb_vObs2D->SetYTitle("Pt (GeV/c)");
   delete histTitle;

   // neutral sidebands Flow observed - Eta
   histTitle = new TString("FlowV0sb_vObsEta_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   histFull[k].histFullHar[j].mHistV0sb_vObsEta = new TProfile(histTitle->Data(), histTitle->Data(), mNEtaBins, mEtaMin, mEtaMax, -100., 100., "");
   histFull[k].histFullHar[j].mHistV0sb_vObsEta->SetXTitle((char*)xLabel.Data());
   histFull[k].histFullHar[j].mHistV0sb_vObsEta->SetYTitle("v (%)");
   delete histTitle;

   // neutral sidebands Flow observed - Pt
   histTitle = new TString("FlowV0sb_vObsPt_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   histFull[k].histFullHar[j].mHistV0sb_vObsPt = new TProfile(histTitle->Data(), histTitle->Data(), nPtBinsPart, ptMin, ptMaxPart, -100., 100., "");
   histFull[k].histFullHar[j].mHistV0sb_vObsPt->SetXTitle("Pt (GeV/c)");
   histFull[k].histFullHar[j].mHistV0sb_vObsPt->SetYTitle("v (%)");
   delete histTitle;

   // SX neutral sidebands Flow observed - Eta
   histTitle = new TString("FlowV0sb_vObsEta_sx_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   histFull[k].histFullHar[j].mHistV0sb_vObsEta_sx = new TProfile(histTitle->Data(), histTitle->Data(), mNEtaBins, mEtaMin, mEtaMax, -100., 100., "");
   histFull[k].histFullHar[j].mHistV0sb_vObsEta_sx->SetXTitle((char*)xLabel.Data());
   histFull[k].histFullHar[j].mHistV0sb_vObsEta_sx->SetYTitle("v (%)");
   delete histTitle;

   // SX neutral sidebands Flow observed - Pt
   histTitle = new TString("FlowV0sb_vObsPt_sx_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   histFull[k].histFullHar[j].mHistV0sb_vObsPt_sx = new TProfile(histTitle->Data(), histTitle->Data(), nPtBinsPart, ptMin, ptMaxPart, -100., 100., "");
   histFull[k].histFullHar[j].mHistV0sb_vObsPt_sx->SetXTitle("Pt (GeV/c)");
   histFull[k].histFullHar[j].mHistV0sb_vObsPt_sx->SetYTitle("v (%)");
   delete histTitle;

   // DX neutral sidebands Flow observed - Eta
   histTitle = new TString("FlowV0sb_vObsEta_dx_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   histFull[k].histFullHar[j].mHistV0sb_vObsEta_dx = new TProfile(histTitle->Data(), histTitle->Data(), mNEtaBins, mEtaMin, mEtaMax, -100., 100., "");
   histFull[k].histFullHar[j].mHistV0sb_vObsEta_dx->SetXTitle((char*)xLabel.Data());
   histFull[k].histFullHar[j].mHistV0sb_vObsEta_dx->SetYTitle("v (%)");
   delete histTitle;

   // DX neutral sidebands Flow observed - Pt
   histTitle = new TString("FlowV0sb_vObsPt_dx_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   histFull[k].histFullHar[j].mHistV0sb_vObsPt_dx = new TProfile(histTitle->Data(), histTitle->Data(), nPtBinsPart, ptMin, ptMaxPart, -100., 100., "");
   histFull[k].histFullHar[j].mHistV0sb_vObsPt_dx->SetXTitle("Pt (GeV/c)");
   histFull[k].histFullHar[j].mHistV0sb_vObsPt_dx->SetYTitle("v (%)");
   delete histTitle;

   // Phi lab
   // Tpc (plus)
   histTitle = new TString("Flow_Phi_TPCplus_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   histFull[k].histFullHar[j].mHistPhiPlus = new TH1D(histTitle->Data(), histTitle->Data(), nPhiBins, phiMin, phiMax);
   histFull[k].histFullHar[j].mHistPhiPlus->SetXTitle("Azimuthal Angles (rad)");
   histFull[k].histFullHar[j].mHistPhiPlus->SetYTitle("Counts");
   delete histTitle;

   // Tpc (minus)
   histTitle = new TString("Flow_Phi_TPCminus_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   histFull[k].histFullHar[j].mHistPhiMinus = new TH1D(histTitle->Data(), histTitle->Data(), nPhiBins, phiMin, phiMax);
   histFull[k].histFullHar[j].mHistPhiMinus->SetXTitle("Azimuthal Angles (rad)");
   histFull[k].histFullHar[j].mHistPhiMinus->SetYTitle("Counts");
   delete histTitle;
  
   // Tpc (cross)
   histTitle = new TString("Flow_Phi_TPCcross_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   histFull[k].histFullHar[j].mHistPhiAll = new TH1D(histTitle->Data(), histTitle->Data(), nPhiBins, phiMin, phiMax);
   histFull[k].histFullHar[j].mHistPhiAll->SetXTitle("Azimuthal Angles (rad)");
   histFull[k].histFullHar[j].mHistPhiAll->SetYTitle("Counts");
   delete histTitle;
  
   // Tpc
   histTitle = new TString("Flow_Phi_TPC_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   histFull[k].histFullHar[j].mHistPhi = new TH1D(histTitle->Data(), histTitle->Data(), nPhiBins, phiMin, phiMax);
   histFull[k].histFullHar[j].mHistPhi->SetXTitle("Azimuthal Angles (rad)");
   histFull[k].histFullHar[j].mHistPhi->SetYTitle("Counts");
   delete histTitle;
  
   // Phi lab flattened
   // Tpc (Plus)
   histTitle = new TString("Flow_Phi_Flat_TPCplus_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   histFull[k].histFullHar[j].mHistPhiFlatPlus = new TH1D(histTitle->Data(), histTitle->Data(), nPhiBins, phiMin, phiMax);
   histFull[k].histFullHar[j].mHistPhiFlatPlus->SetXTitle("Azimuthal Angles (rad)");
   histFull[k].histFullHar[j].mHistPhiFlatPlus->SetYTitle("Counts");
   delete histTitle;

   // Tpc (Minus)
   histTitle = new TString("Flow_Phi_Flat_TPCminus_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   histFull[k].histFullHar[j].mHistPhiFlatMinus = new TH1D(histTitle->Data(), histTitle->Data(), nPhiBins, phiMin, phiMax);
   histFull[k].histFullHar[j].mHistPhiFlatMinus->SetXTitle("Azimuthal Angles (rad)");
   histFull[k].histFullHar[j].mHistPhiFlatMinus->SetYTitle("Counts");
   delete histTitle;

   // Tpc (cross)
   histTitle = new TString("Flow_Phi_Flat_TPCcross_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   histFull[k].histFullHar[j].mHistPhiFlatAll = new TH1D(histTitle->Data(), histTitle->Data(), nPhiBins, phiMin, phiMax);
   histFull[k].histFullHar[j].mHistPhiFlatAll->SetXTitle("Azimuthal Angles (rad)");
   histFull[k].histFullHar[j].mHistPhiFlatAll->SetYTitle("Counts");
   delete histTitle;

   // Tpc
   histTitle = new TString("Flow_Phi_Flat_TPC_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   histFull[k].histFullHar[j].mHistPhiFlat = new TH1D(histTitle->Data(), histTitle->Data(), nPhiBins, phiMin, phiMax);
   histFull[k].histFullHar[j].mHistPhiFlat->SetXTitle("Azimuthal Angles (rad)");
   histFull[k].histFullHar[j].mHistPhiFlat->SetYTitle("Counts");
   delete histTitle;
  }
 }
 if(Debug1)  { cout << spaces << "Histograms booked" << endl ; }
 
 if(fOneInputFile)
 {
  if(Open(GetInputFileName().Data()))  { if(Debug1)  { cout << spaces << "Events File Opened" << endl ; } }
  else 	     			       { cout << FlowAnalysis << "ERROR OPENING THE FILE!" << endl ; return 0 ; }
 }
 else
 {
  if(Debug1)  
  { 
   cout << spaces << "Events Files : " << endl ; 
   Int_t tot = fFileNames.GetEntries() ;
   for(Int_t jj=0;jj<tot;jj++) 
   { 
    cout << spaces << "          " << jj << " : " << ((TObjString*)fFileNames.At(jj))->GetString().Data() << endl ;  
   }
   cout << spaces << " tot : " << tot << " . " << endl ; 
  }
 }

 if(Debug1)  { cout << FlowAnalysis << "Initialized" << endl ; }
 
 return 1 ;
}
//-----------------------------------------------------------------------
Int_t AliFlowAnalysisMaker::Make()                   
{
 // Performs the event loop (0..nEvents) and makes call to : 
 //  - get the event and checks it with AliFlowSelection (dummy) ,
 //  - plugs in the Phi and the Bayesian weights from file (if there) ,
 //  - re-sets the r.p. calculation cuts and re-does the tracks' selection loop (*) ,
 //  - re-shuffles the tracks' array and re-devide in sub-events (*) ,
 //  - fills event histograms (FillEventHistograms()) ,
 //  - fills particle histograms (FillParticleHistograms()) ,
 //  - fills v0 histograms (FillV0Histograms()) (*) .
 // Then makes call to:
 //  - calculate the global resolution and correct flow coefficients (Resolution()) ,
 //  - calculate Phi & Bayesian weights and fill the weight histograms (Weightening()) .
 
 if(Debug0) { cout << FlowAnalysis << "Make (starting events loop) . " << endl ; cout << endl ; }

 Int_t succEv = 0 ;
 Int_t tot = 1 ;
 if(!fOneInputFile) { tot = fFileNames.GetEntries() ; }
 for(Int_t ll=0;ll<tot;ll++)
 {
  if(!fOneInputFile) 
  { 
   TString nextFile = ((TObjString*)fFileNames.At(ll))->GetString().Data() ;
   if(Open(nextFile.Data()))  { if(Debug0)  { cout << spaces << "Just Opened : " << nextFile.Data() << endl ; } }
   else 	     	      { cout << FlowAnalysis << "ERROR OPENING THE FILE: " << nextFile.Data() << endl ; continue ; }
  }
  for(int ie=0;ie<nEvents;ie++)
  {
   evtN = ie ;
   pFlowEvent = GetEvt(evtN) ;
   if(pFlowEvent) 
   {
    if(Analyze(pFlowEvent)) { cout << FlowAnalysis << " ev. " << evtN << " taken" << endl ; succEv++ ; } 
    else                    { cout << FlowAnalysis << " ev. " << evtN << " discarded" << endl ; }
   }
   else 
   {
    if(Debug1) { cout << FlowAnalysis << "Problem opening  ev" << evtN << " . " << endl ; }
    delete pFlowEvent  ; pFlowEvent = 0 ; 
    continue ;
   }
   delete pFlowEvent  ; pFlowEvent = 0 ; 
   pFlowTracks = 0 ;
   pFlowV0s = 0 ;   
  }
  if(!fOneInputFile) 
  {
   pFlowEventsFile->Close() ; 
   if(Debug1) { cout << spaces << "...going for next    " << ll << endl ; } 
  }
 }
 if(!succEv) 
 { 
  if(Debug0) { cout << FlowAnalysis << "NO EVENT HAS BEEN PROCESSED !!!" << endl ; } 
  return kFALSE ; 
 }
 Resolution() ;				                 // calculates the resolution (and corrected v's)
 if(mWritePhiWgt || mRedoWgt) { Weightening() ; }	 // calculates the phi weights 
 
 if(Debug1) { cout << FlowAnalysis << "Made" << endl ; }
  
 return 1 ;
}
//-----------------------------------------------------------------------
Int_t AliFlowAnalysisMaker::Finish()                   
{
 // Close the analysis and saves the histograms on the histFile .

 if(Debug0) { cout << FlowAnalysis << "Finish . " << endl ; } 
 cout << FlowAnalysis << "Output File :  " << GetHistFileName() << endl ;
 cout << endl ;

 // Write all histograms
 histFile->cd() ; histFile->Write() ; 
 histFile->cd() ; mVnResHistList->Write();

// *temp*  --  !uncomment!
 // if(mLabelling) { histFile->cd() ; mLabHist->Write(); delete mLabHist ; }
// *temp*  --  !uncomment!

 histFile->Close() ;
  
 delete mVnResHistList ;
 delete pFlowSelect ;
 
 if(Debug1) { cout << FlowAnalysis << "Finished" << endl ; } 
 return 1 ;
}
//-----------------------------------------------------------------------
// ###
//------------------------------------------------------------------------
void AliFlowAnalysisMaker::WgtChk()
{
 // Check for Read/Write weights ...
 
 TFile* wgtFile = new TFile(GetPhiWgtFileName().Data(),"READ");
 // In case the file does not exist or is not a valid ROOT file, it is made a Zombie. 
 // One can detect this situation with a code like:  f.IsZombie()
 if(wgtFile->IsZombie()) 
 {
  if(Debug0) { cout << spaces << "Weights File ( " << GetPhiWgtFileName().Data() << " ) not present . A new one will be created ." << endl ; }
  mReadPhiWgt = kFALSE ; mWritePhiWgt = kTRUE  ; mRedoWgt = kFALSE ; 
 }
 else 
 { 
  if(Debug0) { cout << spaces << "Reading Weights from File : " << GetPhiWgtFileName().Data() << " . " << endl ; }
  mReadPhiWgt = kTRUE  ; mWritePhiWgt = kFALSE ; 
  FillWgtArrays(wgtFile) ; 			
  if(Debug2) { cout << spaces << "PhiWgt Arrays filled . " << endl ; } 
 }
 wgtFile->Close() ;
}
//-----------------------------------------------------------------------
void AliFlowAnalysisMaker::FillWgtArrays(TFile* wgtFile)
{
 // Loads PhiWeights & Bayesian particles' abundance from file (default: 
 // flowPhiWgt.hist.root). Weights are stored in some AliFlowAnalysisMaker 
 // static data member, ready to be plugged into the AliFlowEvent .
 // This is called at the beginning of the analysis (if wgt file is there).
 
 TString* histTitle ;
 TH1D* TPC_all ; TH1D* TPC_plus ; TH1D* TPC_minus ; TH1D* TPC_cross ;
 TH1D* PID_bay ;
 const int nPhiBins = Flow::nPhiBins ;

 for(int k=0;k<Flow::nSels;k++)
 {
  for(int j=0;j<Flow::nHars;j++) 
  {
   // Tpc (plus)
   histTitle = new TString("Flow_Phi_Weight_TPCplus_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   TPC_plus = (TH1D*)wgtFile->Get(histTitle->Data());
   delete histTitle;
   // Tpc (minus)
   histTitle = new TString("Flow_Phi_Weight_TPCminus_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   TPC_minus = (TH1D*)wgtFile->Get(histTitle->Data());
   delete histTitle;
   // Tpc (cross)
   histTitle = new TString("Flow_Phi_Weight_TPCcross_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   TPC_cross = (TH1D*)wgtFile->Get(histTitle->Data());
   delete histTitle;

   // Tpc
   histTitle = new TString("Flow_Phi_Weight_TPC_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   TPC_all = (TH1D*)wgtFile->Get(histTitle->Data());
   delete histTitle;
  
   for (int n=0;n<nPhiBins;n++) 
   { 
    nPhiWgtPlus[k][j][n]  = TPC_plus->GetBinContent(n+1) ;
    nPhiWgtMinus[k][j][n] = TPC_minus->GetBinContent(n+1) ;
    nPhiWgtCross[k][j][n] = TPC_cross->GetBinContent(n+1) ;
    nPhiWgt[k][j][n]	  = TPC_all->GetBinContent(n+1) ;
    if(Debug2) { cout << "Weights: " << nPhiWgt[k][j][n] << " ; " << nPhiWgtPlus[k][j][n] << " | " << nPhiWgtMinus[k][j][n] << " | " << nPhiWgtCross[k][j][n] << endl ; }
   } 
  }
  // Bayesian weights
  histTitle = new TString("Flow_BayPidMult_Sel");
  *histTitle += k+1;
  PID_bay = (TH1D*)wgtFile->Get(histTitle->Data());
  delete histTitle;
  Double_t totCount = PID_bay->GetSumOfWeights() ;
  for (int n=0;n<Flow::nPid;n++) 
  { 
   if(totCount) { nBayWgt[k][n] = PID_bay->GetBinContent(n+1) / totCount ; }
   else 	{ nBayWgt[k][n] = 1. ; }
   if(Debug2)	{ cout << "Bayesian Weights (" << n << ") : " << nBayWgt[k][n] << endl ; }
  }
 }
 delete TPC_all ; delete TPC_plus ; delete TPC_minus ; delete TPC_cross ;
 delete PID_bay ;
}
//-----------------------------------------------------------------------
void AliFlowAnalysisMaker::FillEvtPhiWgt()
{
 // Plugs Wgt arrays into the current AliFlowEvent.
 
 pFlowEvent->SetPhiWeight(nPhiWgt);
 pFlowEvent->SetPhiWeightPlus(nPhiWgtPlus);
 pFlowEvent->SetPhiWeightMinus(nPhiWgtMinus);
 pFlowEvent->SetPhiWeightCross(nPhiWgtCross); 
}
//-----------------------------------------------------------------------
void AliFlowAnalysisMaker::FillBayesianWgt(int sel)
{
 // Plugs Bayesian particle abundance into the current AliFlowEvent.
 // In principle a bayesian vector is stored for each selection and each
 // harmonic. In practice just the 1st one is then plugged into the event
 // (AliFlowEvent::mBayesianCs[6] is a 1-dimensional array).

 Double_t bayes[Flow::nPid] ; 
 Double_t bayCheck = 0. ;
 for (int n=0;n<Flow::nPid;n++) 
 {
  bayes[n] = nBayWgt[sel][n] ;
  bayCheck += bayes[n] ;   
  if(Debug2) { cout << "Bayesian V[" << n << "]  =  " << nBayWgt[sel][n] << endl ; } 
 }
 if(bayCheck)   { pFlowEvent->SetBayesian(bayes) ; }
 else		{ cout << "An empty bayesian vector is stored !!! - Bayesian weights = {1,1,1,1,1,1} " << endl ; }
 mRePid = kTRUE ;
}
//-----------------------------------------------------------------------
void AliFlowAnalysisMaker::Weightening()
{
 // Calculates weights, and fills PhiWgt histograms and saves them into
 // the PhiWgt file (by default: flowPhiWgt.hist.root) .
 // This is called at the end of the analysis.
 
 if(Debug0) { cout << FlowAnalysis << "Calculating/Filling weights . " << endl ; } 
 cout << endl ;
 
 // PhiWgt histogram collection
 mPhiWgtHistList = new TOrdCollection(Flow::nSels*Flow::nHars) ;
 
 const float phiMin   = 0. ;
 const float phiMax   = 2*TMath::Pi() ; 
 const int   nPhiBins = Flow::nPhiBins ;  // 120 
 TString* histTitle ;

 for(int k = 0; k < Flow::nSels; k++)
 {
  for(int j = 0; j < Flow::nHars; j++) 
  {
   // Creates PhiWgt Histograms

   // Tpc (plus)
   histTitle = new TString("Flow_Phi_Weight_TPCplus_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   histFull[k].histFullHar[j].mHistPhiWgtPlus = new TH1D(histTitle->Data(),histTitle->Data(), nPhiBins, phiMin, phiMax);
   histFull[k].histFullHar[j].mHistPhiWgtPlus->Sumw2();
   histFull[k].histFullHar[j].mHistPhiWgtPlus->SetXTitle("Azimuthal Angles (rad)");
   histFull[k].histFullHar[j].mHistPhiWgtPlus->SetYTitle("PhiWgt");
   delete histTitle;
   // Tpc (minus)
   histTitle = new TString("Flow_Phi_Weight_TPCminus_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   histFull[k].histFullHar[j].mHistPhiWgtMinus = new TH1D(histTitle->Data(),histTitle->Data(), nPhiBins, phiMin, phiMax);
   histFull[k].histFullHar[j].mHistPhiWgtMinus->Sumw2();
   histFull[k].histFullHar[j].mHistPhiWgtMinus->SetXTitle("Azimuthal Angles (rad)");
   histFull[k].histFullHar[j].mHistPhiWgtMinus->SetYTitle("PhiWgt");
   delete histTitle;
   // Tpc (cross)
   histTitle = new TString("Flow_Phi_Weight_TPCcross_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   histFull[k].histFullHar[j].mHistPhiWgtAll = new TH1D(histTitle->Data(),histTitle->Data(), nPhiBins, phiMin, phiMax);
   histFull[k].histFullHar[j].mHistPhiWgtAll->Sumw2();
   histFull[k].histFullHar[j].mHistPhiWgtAll->SetXTitle("Azimuthal Angles (rad)");
   histFull[k].histFullHar[j].mHistPhiWgtAll->SetYTitle("PhiWgt");
   delete histTitle;

   // Tpc
   histTitle = new TString("Flow_Phi_Weight_TPC_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   histFull[k].histFullHar[j].mHistPhiWgt = new TH1D(histTitle->Data(),histTitle->Data(), nPhiBins, phiMin, phiMax);
   histFull[k].histFullHar[j].mHistPhiWgt->Sumw2();
   histFull[k].histFullHar[j].mHistPhiWgt->SetXTitle("Azimuthal Angles (rad)");
   histFull[k].histFullHar[j].mHistPhiWgt->SetYTitle("PhiWgt");
   delete histTitle;

   // Calculate PhiWgt
   double meanPlus  = histFull[k].histFullHar[j].mHistPhiPlus->Integral() / (double)nPhiBins ;
   double meanMinus = histFull[k].histFullHar[j].mHistPhiMinus->Integral() / (double)nPhiBins ;
   double meanCross = histFull[k].histFullHar[j].mHistPhiAll->Integral() / (double)nPhiBins ;

   double meanTPC = histFull[k].histFullHar[j].mHistPhi->Integral() / (double)nPhiBins ;

   // Tpc
   for (int i=0;i<nPhiBins;i++) 
   {
    histFull[k].histFullHar[j].mHistPhiWgtPlus->SetBinContent(i+1,meanPlus);
    histFull[k].histFullHar[j].mHistPhiWgtPlus->SetBinError(i+1, 0.);
    histFull[k].histFullHar[j].mHistPhiWgtMinus->SetBinContent(i+1,meanMinus);
    histFull[k].histFullHar[j].mHistPhiWgtMinus->SetBinError(i+1, 0.);
    histFull[k].histFullHar[j].mHistPhiWgtAll->SetBinContent(i+1,meanCross);
    histFull[k].histFullHar[j].mHistPhiWgtAll->SetBinError(i+1, 0.);

    histFull[k].histFullHar[j].mHistPhiWgt->SetBinContent(i+1,meanTPC);
    histFull[k].histFullHar[j].mHistPhiWgt->SetBinError(i+1, 0.);
   }
   histFull[k].histFullHar[j].mHistPhiWgtPlus->Divide(histFull[k].histFullHar[j].mHistPhiPlus);
   mPhiWgtHistList->AddLast(histFull[k].histFullHar[j].mHistPhiWgtPlus);
   histFull[k].histFullHar[j].mHistPhiWgtMinus->Divide(histFull[k].histFullHar[j].mHistPhiMinus);
   mPhiWgtHistList->AddLast(histFull[k].histFullHar[j].mHistPhiWgtMinus);
   histFull[k].histFullHar[j].mHistPhiWgtAll->Divide(histFull[k].histFullHar[j].mHistPhiAll);
   mPhiWgtHistList->AddLast(histFull[k].histFullHar[j].mHistPhiWgtAll);

   histFull[k].histFullHar[j].mHistPhiWgt->Divide(histFull[k].histFullHar[j].mHistPhi);
   mPhiWgtHistList->AddLast(histFull[k].histFullHar[j].mHistPhiWgt);
  }
 }
 if(Debug2) { mPhiWgtHistList->ls() ; }

 // Write PhiWgt histograms
 TString wgtName ;
 if(mRedoWgt) { wgtName = "new." ; }
 else 	      { wgtName = "" ; }
 wgtName += GetPhiWgtFileName() ;
 TFile* phiWgtFile = new TFile(wgtName.Data(), "RECREATE");
 phiWgtFile->cd() ; mPhiWgtHistList->Write();
 delete mPhiWgtHistList ; delete wgtName ;
 
 // Write Bayesian Weights for P.Id.
 phiWgtFile->cd() ; 
 for(int k=0;k<Flow::nSels;k++) { histFull[k].mHistBayPidMult->Write() ; }
 phiWgtFile->Close();

 if(Debug1) { cout << FlowAnalysis << "Done . " << endl ; cout << endl ; }

 return ;
}
//----------------------------------------------------------------------
void AliFlowAnalysisMaker::GetRunBayesian(Double_t bayes[Flow::nPid],int sel)  
{
 // Returns the normalized particle abundance of all the events.
 // This is called at the end of the analysis.

 if(sel>Flow::nSels) 
 { 
  if(Debug1) { cout << "Wrong Selection! " << endl ; } 
  return ; 
 }
 Double_t totCount = histFull[sel].mHistBayPidMult->GetSumOfWeights() ;
 for(int i=0;i<Flow::nPid;i++)
 {
  if(totCount) { bayes[i] = histFull[sel].mHistBayPidMult->GetBinContent(i+1) / totCount ; }
  else         { bayes[i] = 1. ; }
 }
 return ;
}
//----------------------------------------------------------------------
void AliFlowAnalysisMaker::PrintRunBayesian(int sel)  
{
 // Prints the normalized particle abundance of all the events (at this step).

 if(sel>Flow::nSels) 
 { 
  if(Debug1) { cout << "Wrong Selection! " << endl ; } 
  return ; 
 }
 Double_t totCount = histFull[sel].mHistBayPidMult->GetSumOfWeights() ;
 Char_t* names[Flow::nPid] = {"e","mu","pi","k","p","d"} ;
 Double_t bayes = 0. ;
 cout << "  " ;
 for(int i=0;i<Flow::nPid;i++)
 {
  if(totCount) { bayes = histFull[sel].mHistBayPidMult->GetBinContent(i+1) / totCount ; }
  else         { bayes = 1. ; }
  cout << bayes << "_" << names[i] << " ; " ;
 }
 cout << endl ;
 return ;
}
//-----------------------------------------------------------------------
// ###
//-----------------------------------------------------------------------
Bool_t AliFlowAnalysisMaker::Open(const Char_t* filename)              
{
 // Checks/Opens the FlowEvents file for analysis. Sets the n. of events 
 // (for the loop) and sets the internal counter evtN to 0.

 if(Debug1) { cout << spaces << "...opening file . " << endl ; }
 //if(pFlowEventsFile) { if(Debug1) { cout << spaces << "...file already open (I'll go on) . " << endl ; } return kTRUE ; } 
 
 if(Debug0)  { cout << FlowAnalysis << "Opening input-File :  " << filename << "  . " << endl ; }
 pFlowEventsFile = new TFile(filename ,"READ") ; 
 if(!pFlowEventsFile || pFlowEventsFile->IsZombie())
 {
  if(Debug0) { cout << FlowAnalysis << "flow events file NOT FOUND ( " << filename << " ) . " << endl ; }
  return kFALSE ;
 }
 else
 {
  nEvents = pFlowEventsFile->GetNkeys() ; 
  pFlowEventsList = (TList*)pFlowEventsFile->GetListOfKeys() ; 
  evtN = 0 ;  
  if(Debug0) { cout << spaces << "file: " << GetInputFileName().Data() << " , found " << nEvents << " Flow Events . " << endl ; }
  if(Debug2) { pFlowEventsFile->ls() ; } //if(Debug2) { cout << " Lista : " << endl ; pFlowEventsList->Dump() ; }
 }
 if(Debug0) { cout << endl ; }

 return kTRUE ;
}
//-----------------------------------------------------------------------
AliFlowEvent* AliFlowAnalysisMaker::GetEvt(Int_t evt)         
{
 // Gets the FlowEvent n.evt from file, loads Tracks and V0s arrays into memory
 // if evt is not specified, gets the next event basing on the internal counter evtN (0...N)
 
 if(Debug1) { cout << spaces << "Getting Event " << evt << " . " << endl ; } 
 
 if(!pFlowEventsFile) { return 0 ; }  	 // if no file open -> quit!
 else
 {
  if(evt>=0) { evtN = evt ; }	 // specified event "evt" is selected, the iterator "evtN" changes accordingly
  //evt_name = "ev" ; evt_name += evtN ;
  if(evtN<nEvents) //&& (pFlowEventsFile->GetKey(evt_name)))
  {
   // gets the name of the key
   evt_name = pFlowEventsList->At(evtN)->GetName() ;
   if(Debug1) { cout << "####" << endl ; cout << evt_name << endl ; cout << "####" << endl ; }

   // gets the event from file
   pFlowEventsFile->GetObject(evt_name.Data(),pFlowEvent) ; 			// new way
   //pFlowEvent = (AliFlowEvent*)pFlowEventsFile->Get(evt_name.Data()) ; 	// old way

   if(Debug2) { cout << "     dumping: " << endl ; pFlowEvent->Dump() ; cout << endl ; }
    
   // loading tracks and v0s collections
   pFlowTracks = pFlowEvent->TrackCollection() ; 
   nTracks = pFlowTracks->GetEntries() ;
   if(Debug1) { cout << FlowAnalysis << "event: " << evt_name.Data() << " , found " << nTracks << " Flow Tracks . " << endl ; }
   pFlowV0s = pFlowEvent->V0Collection() ; 
   nV0s = pFlowV0s->GetEntries() ;
   if(Debug1) { cout << FlowAnalysis << "event: " << evt_name.Data() << " , found " << nV0s << " Flow V0s . " << endl ; }
  }
  else 
  { 
   cout << FlowAnalysis << "! no event !" << endl ; 
   return 0 ;
  }
 }
 return pFlowEvent ;
}
//-----------------------------------------------------------------------
Bool_t AliFlowAnalysisMaker::Analyze(AliFlowEvent* pFlowEvent)         
{
 // Runs the analysis on a single AliFlowEvent (pointer given). 
 // In principle this method can be used externally by the 
 // AliSelectorFlow for on-fly analysis.
 
//
cout << "* 1 . " << endl ;
//
 if(pFlowSelect->Select(pFlowEvent))		   // event selected - here below the ANALYSIS FLAGS are setted -
 {
//
cout << "* 2 . " << endl ;
//
  if(mReadPhiWgt)					 // if weights file is there (previously checked by WgtChk()
  { 
   if(mPhiWgt)    { FillEvtPhiWgt() ;  }		 // pFlowEvent->SetPhiWgt(...)
   if(mOnePhiWgt) { pFlowEvent->SetOnePhiWgt() ; }	 // one phi-wgt histogram
   else 	  { pFlowEvent->SetFirstLastPhiWgt() ; } // three phi-wgt histogram
  }
  if(mBayWgt)	  { FillBayesianWgt() ; }		 // pFlowEvent->SetBayesian(...)
  if(mRePid)	  { pFlowEvent->SetPids() ; }		 // re-calculate all p.id. hypotesis with the (new) bayesian array

  if(mPtWgt)	  { pFlowEvent->SetPtWgt(); ; } 	 // pT as a weight
  if(mEtaWgt)	  { pFlowEvent->SetEtaWgt() ; } 	 // eta as a weight

  pFlowEvent->SetSelections(pFlowSelect) ;		 // does the selection of tracks for r.p. calculation (sets flags in AliFlowTrack)
  if(mShuffle)    { pFlowEvent->RandomShuffle() ; }	 // tracks re-shuffling
  pFlowEvent->SetEtaSubs(mEtaSub) ;			 // setting for the subevents (eta or random)
  pFlowEvent->MakeSubEvents() ; 			 // makes the subevent, eta or random basing on the previous flag

  if(mMakeAll)
  {
   if(Debug1) { cout << FlowAnalysis << "All calculation in one shoot for ev." << evtN << " . " << endl ; }
   pFlowEvent->MakeAll() ; 
  }
  else 
  {
   if(Debug1) { cout << FlowAnalysis << "Single calls for each Q,Psi,... for ev." << evtN << " . " << endl ; }
  }
  if(FillFromFlowEvent())			   // calculates event quantities
  {
//  
cout << "* 3 . " << endl ;
//
   FillEventHistograms();			   // fill histograms from AliFlowEvents
   FillParticleHistograms();			   // fill histograms from AliFlowTracks
   if(mV0)	  { FillV0Histograms() ; }	   // fill histograms from AliFlowV0s 
   if(mLabelling) { FillLabels() ; }		   // fill the histogram of MC labels (from the simulation)
  }
  else 
  {
   cout << FlowAnalysis << "Event psi = 0 . " << endl ; 
   cout << FlowAnalysis << "				Skipping! " << endl ; 
   return kFALSE ;
  }
  
 }
 else 
 {
  if(Debug1) { cout << FlowAnalysis << "Event  ev" << evtN << " discarded . " << endl ; }
  delete pFlowEvent  ; pFlowEvent = 0 ; 
  return kFALSE ;
 }
 
 return kTRUE ;
}
//-----------------------------------------------------------------------
Bool_t AliFlowAnalysisMaker::FillFromFlowEvent()
{
 // gets event quantities ;

 if(Debug1) { cout << spaces << "Getting Event quantities . " << endl ; } 
 
 Int_t selCheck = 0 ;
 for(int k = 0; k < Flow::nSels; k++) 
 {
  pFlowSelect->SetSelection(k) ;
  for(int j = 0; j < Flow::nHars; j++) 
  {
   pFlowSelect->SetHarmonic(j) ;
   for(int n = 0; n < Flow::nSubs; n++) 
   {
    pFlowSelect->SetSubevent(n) ;
    mPsiSub[n][k][j] = pFlowEvent->Psi(pFlowSelect) ;  	// sub-event quantities
    mMultSub[n][k][j] = pFlowEvent->Mult(pFlowSelect) ;
   }
   pFlowSelect->SetSubevent(-1);
   mQ[k][j]    = pFlowEvent->Q(pFlowSelect) ;	   	// full event quantities
   mPsi[k][j]  = pFlowEvent->Psi(pFlowSelect) ;
   m_q[k][j]   = pFlowEvent->NormQ(pFlowSelect).Mod() ; // was: pFlowEvent->q(pFlowSelect) ; // but the normalization was bad (no pT,eta weight)
   mMult[k][j] = pFlowEvent->Mult(pFlowSelect) ;
   selCheck += mMult[k][j] ; 
  }
 }
 if(!selCheck) { return kFALSE ; }           		// if there are no particles in the selection -> skip the event
 
 if(Debug1) 						// prints event by event calculated quantities
 {
  cout << FlowAnalysis << "Read from FlowEvent : " << endl ;
  cout << endl ; cout << spaces << "mQ[k][j] = " ; 
  for(int k=0;k<Flow::nSels;k++) 
  {
   for(int j=0;j<Flow::nHars;j++) 
   { 
    cout << (Float_t)mQ[k][j].X() << "," << (Float_t)mQ[k][j].Y() << " ; " ; 
   }
   cout << endl ;
   cout << spaces << "            " ; 
  }
  cout << endl ; cout << spaces << "mPsi[k][j] = " ; 
  for(int k=0;k<Flow::nSels;k++) 
  {
   for(int j=0;j<Flow::nHars;j++) { Float_t aaa = (Float_t)mPsi[k][j] ; cout << aaa << " , " ; }
   cout << endl ;
   cout << spaces << "              " ; 
  }
  cout << endl ; cout << spaces << "m_q[k][j] = " ; 
  for(int k=0;k<Flow::nSels;k++) 
  {
   for(int j=0;j<Flow::nHars;j++) { Float_t aaa = (Float_t)m_q[k][j] ; cout << aaa << " , " ; }
   cout << endl ;
   cout << spaces << "             " ; 
  }
  cout << endl ; cout << spaces << "mMult[k][j] = " ; 
  for(int k=0;k<Flow::nSels;k++) 
  {
   for(int j=0;j<Flow::nHars;j++) { Float_t aaa = (Float_t)mMult[k][j] ; cout << aaa << " , " ; } 
   cout << endl ;
   cout << spaces << "               " ; 
  }
  if(Debug1)
  {
   cout << endl ; cout << spaces << "mMultSub[n][k][j] = " ; 
   for(int n=0;n<Flow::nSubs;n++) 
   {
    for(int k=0;k<Flow::nSels;k++) 
    {
     for(int j=0;j<Flow::nHars;j++) { Float_t aaa = mMultSub[n][k][j] ; cout << aaa << " , " ; } 
     cout << endl ;
     cout << spaces << "                 " ; 
    }
   }
   cout << endl ; cout << spaces << "mPsiSub[n][k][j] = " ; 
   for(int n=0;n<Flow::nSubs;n++) 
   {
    for(int k=0;k<Flow::nSels;k++) 
    {
     for(int j=0;j<Flow::nHars;j++) { Float_t aaa = mPsiSub[n][k][j] ; cout << aaa << " , " ; } 
     cout << endl ;
     cout << spaces << "                 " ; 
    }
   }
   cout << endl ; cout << spaces << "Delta_PsiSub[k][j] = " ; 
   for(int k=0;k<Flow::nSels;k++) 
   {
    for(int j=0;j<Flow::nHars;j++) { Float_t aaa = mPsiSub[0][k][j]-mPsiSub[1][k][j] ; cout << aaa << " , " ; } 
    cout << endl ;
    cout << spaces << " 		  " ; 
   }
  }
  cout << endl ;
 }
 return kTRUE ;
}
//-------------------------------------------------------------
void AliFlowAnalysisMaker::FillEventHistograms()    
{
 // event histograms

 if(Debug1) { cout << spaces << "Event Histograms . " << endl ; } 

 Float_t trigger = (Float_t)pFlowEvent->L0TriggerWord() ;
 mHistTrigger->Fill(trigger);

 // no selections: OrigMult, Centrality, Mult, MultOverOrig, VertexZ, VertexXY
 int origMult = pFlowEvent->OrigMult();
 mHistOrigMult->Fill((float)origMult);
 mHistMultEta->Fill((float)pFlowEvent->MultEta());

 int cent = pFlowEvent->Centrality();
 mHistCent->Fill((float)cent);

 mHistMult->Fill((float)nTracks) ;
 mHistV0Mult->Fill((float)nV0s) ;
 if(origMult) { mHistMultOverOrig->Fill((float)nTracks/(float)origMult) ; }

 pFlowEvent->VertexPos(vertex) ;
 mHistVertexZ->Fill(vertex[2]) ;
 mHistVertexXY2D->Fill(vertex[0],vertex[1]) ;

 // ZDC info
 mHistPartZDC->Fill(pFlowEvent->ZDCpart()) ;
 for(int ii=0;ii<3;ii++) { mHistEnergyZDC->Fill(ii,pFlowEvent->ZDCenergy(ii)) ; }

 // sub-event Psi_Subs
 for(int k = 0; k < Flow::nSels; k++) 
 {
  for(int j = 0; j < Flow::nHars; j++) 
  {
   for(int n = 0; n < Flow::nSubs; n++) 
   {
    int iii = Flow::nSubs * k + n ;    //cout << "  " << k << j << n << " , " << iii << endl ;
    histSub[iii].histSubHar[j].mHistPsiSubs->Fill(mPsiSub[n][k][j]) ;
   }
  }
 }

 // full event Psi, PsiSubCorr, PsiSubCorrDiff, cos, mult, q
 for(int k = 0; k < Flow::nSels; k++) 
 {
  for(int j = 0; j < Flow::nHars; j++) 
  {
   float order = (float)(j+1);
   histFull[k].histFullHar[j].mHistPsi->Fill(mPsi[k][j]);
   if(k<2 && j<2)
   {
    if(k==0) 
    { 
     float psi1 = mPsi[0][j] ; 
     float psi2 = mPsi[1][j] ;
     float diff = psi1 - psi2 ;
     if(diff < -TMath::Pi()/(j+1))      { diff += 2*TMath::Pi()/(j+1) ; } 
     else if(diff > +TMath::Pi()/(j+1)) { diff -= 2*TMath::Pi()/(j+1) ; }
     histFull[k].histFullHar[j].mHistPsi_Diff->Fill(diff) ; // k=0
    } 
    else if(k==1) 
    {
     float psi1 ; float psi2 ;
     if (j==0)     { psi1 = mPsi[0][0] ; psi2 = mPsi[1][1] ; }  
     else if(j==1) { psi1 = mPsi[0][0] ; psi2 = mPsi[0][1] ; }
     float diff = psi1 - psi2 ;
     diff = (TMath::Abs(diff) > TMath::Pi()) ? ((diff > 0.) ? -(2*TMath::Pi()-diff) : -(diff+2*TMath::Pi())) : diff ;
     histFull[k].histFullHar[j].mHistPsi_Diff->Fill(diff) ; // k=1
    }	   
   }

   if(mPsiSub[0][k][j] != 0. && mPsiSub[1][k][j] != 0.)
   {
    float psiSubCorr;    			// this is:  delta_Psi
    if(mV1Ep1Ep2 == kFALSE || order != 1) 
    {
     psiSubCorr = mPsiSub[0][k][j] - mPsiSub[1][k][j];
    }
    else // i.e. (mV1Ep1Ep2 == kTRUE && order == 1)
    { 
     psiSubCorr = mPsiSub[0][k][0] + mPsiSub[1][k][0] - 2*mPsi[k][1];
    }
    histFull[k].mHistCos->Fill(order, (float)cos(order * psiSubCorr)) ;
    if(psiSubCorr < 0.) 		 { psiSubCorr += 2*TMath::Pi()/order ; }
    if(psiSubCorr > 2*TMath::Pi()/order) { psiSubCorr -= 2*TMath::Pi()/order ; } // for v1Ep1Ep2 which gives -2*TMath::Pi() < psiSubCorr < 2*2*TMath::Pi()
    histFull[k].histFullHar[j].mHistPsiSubCorr->Fill(psiSubCorr);
   }

   if(j < Flow::nHars - 1) // subevents of different harmonics
   {
    int j1, j2;
    float psiSubCorrDiff;
    if(j==0)	  { j1 = 1, j2 = 2 ; } 
    else if(j==1) { j1 = 1, j2 = 3 ; } 
    else if(j==2) { j1 = 2, j2 = 4 ; }
    psiSubCorrDiff = fmod((double)mPsiSub[0][k][j1-1],2*TMath::Pi()/(double)j2)-fmod((double)mPsiSub[1][k][j2-1],2*TMath::Pi()/(double)j2) ;
    if(psiSubCorrDiff < 0.) { psiSubCorrDiff += 2*TMath::Pi()/(float)j2 ; }
    histFull[k].histFullHar[j].mHistPsiSubCorrDiff->Fill(psiSubCorrDiff) ;
    psiSubCorrDiff = fmod((double)mPsiSub[0][k][j2-1],2*TMath::Pi()/(double)j2)-fmod((double)mPsiSub[1][k][j1-1],2*TMath::Pi()/(double)j2) ;
    if(psiSubCorrDiff < 0.) { psiSubCorrDiff += 2*TMath::Pi()/(float)j2 ; }
    histFull[k].histFullHar[j].mHistPsiSubCorrDiff->Fill(psiSubCorrDiff) ;
   }
   
   histFull[k].histFullHar[j].mHistMult->Fill((float)mMult[k][j]) ;
   histFull[k].histFullHar[j].mHist_q->Fill(m_q[k][j]) ;
  }
 }
 if(Debug1) { cout << FlowAnalysis << "Event  " << (pFlowEvent->GetName()) << endl ; }
}
//-----------------------------------------------------------------------
void AliFlowAnalysisMaker::FillParticleHistograms() 
{
 // tracks histograms

 if(Debug0) { cout << spaces << "Tracks Loop . " << endl ; } 

 float corrMultUnit   = 0. ;
 float corrMultN      = 0. ;
 float etaSymPosTpcN  = 0. ;
 float etaSymNegTpcN  = 0. ;
 float etaSymPosTpcNpart = 0. ;
 float etaSymNegTpcNpart = 0. ;
 float hPlusN	      = 0. ;
 float hMinusN        = 0. ;
 float piPlusN        = 0. ;
 float piMinusN       = 0. ;
 float protonN        = 0. ;
 float pbarN	      = 0. ;
 float kMinusN        = 0. ;
 float kPlusN	      = 0. ;
 float deuteronN      = 0. ;
 float dbarN	      = 0. ;
 float electronN      = 0. ;
 float positronN      = 0. ;
 float muonMinusN     = 0. ;
 float muonPlusN      = 0. ;
 int numPid = -1 ;

 mOnlyConstrainable = pFlowSelect->JustLoopConstrainable() ;
 if(Debug1) { cout << spaces << "Looping over " << nTracks << "  flow tracks  " ; if(mOnlyConstrainable) { cout << " (just constrainable) " ; } cout << ". " << endl ; } 
 int itr ;
 for(itr=0;itr<nTracks;itr++) 
 {
  pFlowTrack = (AliFlowTrack*)pFlowTracks->At(itr) ; 	if(Debug2) { cout << "Track n. " << itr << endl ; pFlowTrack->Dump() ; }
  
  bool constrainable = pFlowTrack->IsConstrainable() ;  if(Debug2) {  
  							 if(constrainable) { cout << "Constrainable" << endl ; }
  							 else		   { cout << "UnConstrainable" << endl ; }
  							}
  if(mOnlyConstrainable && !constrainable) { continue ; } 
  
  int label = pFlowTrack->Label() ; 			if(Debug2) { cout << itr << " label: " << label << endl ; }
  if(label>maxLabel) { maxLabel = label ; }             // this has just the purpose to fill a labels histogram
  // -
  float dcaGlobal = TMath::Abs(pFlowTrack->Dca()) ; 	if(Debug2) { cout << "dcaGlob "  << dcaGlobal << endl ; }	        
  float dcaSigned = pFlowTrack->TransDcaSigned() ;	if(Debug2) { cout << "dcaSigned " << dcaSigned  << endl ; }	        
  float dcaTrans = TMath::Abs(dcaSigned) ; 		if(Debug2) { cout << "dcaTrans "  << dcaTrans << endl ; }
  float eta = pFlowTrack->Eta() ; 			if(Debug2) { cout << "eta "  << eta << endl ; }
  							if(Debug2) { cout << "rapidity " << pFlowTrack->Y() << endl ; }
  float phi = pFlowTrack->Phi() ; 			if(Debug2) { cout << "Phi "  << phi<< endl ; } 
  float pt = pFlowTrack->Pt() ; 			if(Debug2) { cout << "Pt  "   << pt << endl ; }
  float totalp = pFlowTrack->P() ; 			if(Debug2) { cout << "P "    << totalp << endl ; }
  // float logp = TMath::Log10(totalp) ;
  float etaGlob = pFlowTrack->EtaGlobal() ; 		if(Debug2) { cout << "etaGlobal "  << eta << endl ; }
  							if(Debug2) { cout << "rapidityGlobal " << pFlowTrack->YGlobal() << endl ; }
  float phiGlob = pFlowTrack->PhiGlobal() ; 		if(Debug2) { cout << "PhiGlobal "  << phi<< endl ; } 
  float ptGlob = pFlowTrack->PtGlobal() ; 		if(Debug2) { cout << "PtGlobal  "   << pt << endl ; }
  float zFirstPoint = pFlowTrack->ZFirstPoint() ; 	if(Debug2) { cout << "Zfirst " << zFirstPoint << endl ; }
  float zLastPoint = pFlowTrack->ZLastPoint() ; 	if(Debug2) { cout << "Zlast " << zLastPoint << endl ; }
  float lenght = pFlowTrack->TrackLength() ;     	if(Debug2) { cout << "lenght " << lenght << endl ; }
  int	charge = pFlowTrack->Charge() ; 		if(Debug2) { cout << "charge " << charge << endl ; }
  float chi2 = pFlowTrack->Chi2() ; 			if(Debug2) { cout << "chi2 " << chi2 << endl ; } 
  int	fitPtsTPC = (int)((float)pFlowTrack->FitPtsTPC()) ;	if(Debug2) { cout << "fitPtsTPC " << fitPtsTPC << endl ; }
  int	maxPtsTPC = pFlowTrack->MaxPtsTPC() ;		if(Debug2) { cout << "maxPtsTPC " << maxPtsTPC << endl ; }	    
  float chi2TPC = pFlowTrack->Chi2TPC() ; 		if(Debug2) { cout << "chi2TPC " << chi2TPC << endl ; }
  int	fitPtsITS = pFlowTrack->FitPtsITS() ;		if(Debug2) { cout << "fitPtsITS " << fitPtsITS << endl ; } 	
  int	maxPtsITS = pFlowTrack->MaxPtsITS() ;	        if(Debug2) { cout << "maxPtsITS " << maxPtsITS << endl ; }  
  float chi2ITS = pFlowTrack->Chi2ITS() ;	        if(Debug2) { cout << "chi2ITS " << chi2ITS << endl ; }
  int	fitPtsTRD = pFlowTrack->NhitsTRD() ;		if(Debug2) { cout << "fitPtsTRD " << fitPtsTRD << endl ; }
  int	maxPtsTRD = pFlowTrack->MaxPtsTRD() ;		if(Debug2) { cout << "maxPtsTRD " << maxPtsTRD << endl ; }	    
  float chi2TRD = pFlowTrack->Chi2TRD() ;	        if(Debug2) { cout << "chi2TRD " << chi2TRD << endl ; }
  int	fitPtsTOF = pFlowTrack->NhitsTOF() ;		if(Debug2) { cout << "fitPtsTOF " << fitPtsTOF << endl ; }
  int	maxPtsTOF = pFlowTrack->MaxPtsTOF() ;		if(Debug2) { cout << "maxPtsTOF " << maxPtsTOF << endl ; }	    
  float chi2TOF = pFlowTrack->Chi2TOF() ;	        if(Debug2) { cout << "chi2TOF " << chi2TOF << endl ; }
  float dEdx = pFlowTrack->DedxTPC() ; 			if(Debug2) { cout << "dEdxTPC " << dEdx << endl ; }
  float its = pFlowTrack->DedxITS() ;		        if(Debug2) { cout << "dEdxITS " << its << endl ; }
  float trd = pFlowTrack->SigTRD() ;		        if(Debug2) { cout << "SigTRD " << trd << endl ; }
  float tof = pFlowTrack->TofTOF() ;		        if(Debug2) { cout << "TofTOF " << tof << endl ; }

  float lpTPC = 0 ; if(pFlowTrack->PatTPC()>0) { lpTPC = TMath::Log10(pFlowTrack->PatTPC()) ; }
  							if(Debug2) { cout << "log(pTPC) " << lpTPC << endl ; }							
  float lpITS = 0 ; if(pFlowTrack->PatITS()>0) { lpITS = TMath::Log10(pFlowTrack->PatITS()) ; }
  							if(Debug2) { cout << "log(pTPC) " << lpITS << endl ; }
  float lpTRD = 0 ; if(pFlowTrack->PatTRD()>0) { lpTRD = TMath::Log10(pFlowTrack->PatTRD()) ; }
  							if(Debug2) { cout << "log(pTPC) " << lpTRD << endl ; }
  float lpTOF = 0 ; if(pFlowTrack->PatTOF()>0) { lpTOF = TMath::Log10(pFlowTrack->PatTOF()) ; }
  							if(Debug2) { cout << "log(pTPC) " << lpTOF << endl ; }  
  
  float invMass = pFlowTrack->InvMass() ;	        if(Debug2) { cout << "InvMass " << invMass << endl ; }
  Char_t pid[10]="0" ; strcpy(pid,pFlowTrack->Pid()) ; 	if(Debug2) { cout << "pid "  << pid << endl ; }   

  // no selections: Charge, Dca, Chi2, FitPts, MaxPts, FitOverMax, PID
  mHistCharge->Fill((float)charge);
  mHistDcaGlobal->Fill(dcaGlobal);
  mHistDca->Fill(dcaTrans) ;
  mHistTransDca->Fill(dcaSigned);
  mHistChi2->Fill(chi2);
  
  // - here ITS   (chi2 & nHits)
  mHistChi2ITS->Fill(chi2ITS);
  if (fitPtsITS>0)
    mHistChi2normITS->Fill(chi2ITS/((float)fitPtsITS));  // 5 is the # of parameters in the track fit
  mHistFitPtsITS->Fill((float)fitPtsITS);
  mHistMaxPtsITS->Fill((float)maxPtsITS);

  // - here TPC   (chi2 & nHits)
  mHistChi2TPC->Fill(chi2TPC);
  if (fitPtsTPC>0)
    mHistChi2normTPC->Fill(chi2TPC/((float)fitPtsTPC));
  mHistFitPtsTPC->Fill((float)fitPtsTPC);
  mHistMaxPtsTPC->Fill((float)maxPtsTPC);
  if(maxPtsTPC>0)
    mHistFitOverMaxTPC->Fill((float)(fitPtsTPC)/(float)maxPtsTPC);

  // - here TRD  (chi2 & nHits)
  mHistChi2TRD->Fill(chi2TRD);
  if (fitPtsTRD>0)
    mHistChi2normTRD->Fill(chi2TRD/((float)fitPtsTRD));
  mHistFitPtsTRD->Fill((float)fitPtsTRD);
  mHistMaxPtsTRD->Fill((float)maxPtsTRD);

  // - here TOF   (chi2 & nHits)
  mHistChi2TOF->Fill(chi2TOF);
  if (fitPtsTOF>0)
    mHistChi2normTOF->Fill(chi2TOF/((float)fitPtsTOF));
  mHistFitPtsTOF->Fill((float)fitPtsTOF);
  mHistMaxPtsTOF->Fill((float)maxPtsTOF);
  
  // fit over max (all)
  int maxPts = maxPtsITS + maxPtsTPC + maxPtsTRD + maxPtsTOF ;
  int fitPts = fitPtsITS + fitPtsTPC + fitPtsTRD + fitPtsTOF ;
  if(maxPts>0)
    mHistFitOverMax->Fill((float)(fitPts)/(float)maxPts) ;
  
  // lenght
  mHistLenght->Fill(lenght) ;
  mHistInvMass->Fill(invMass) ;

  // PID histograms & multiplicity count (for bayesian histogram)
  if(charge == 1) 
  {
   hPlusN++ ;
   mHistMeanDedxPos2D->Fill(lpTPC, dEdx) ;
   mHistMeanDedxPos2DITS->Fill(lpITS, its) ;
   mHistMeanDedxPos2DTRD->Fill(lpTRD, trd) ;
   mHistMeanDedxPos2DTOF->Fill(lpTOF, tof) ;
   //
   float positron  = pFlowTrack->ElectronPositronProb() ;
   mHistPidPositron->Fill(positron) ;
   if(strcmp(pid, "e+") == 0) 
   {
    numPid = 0 ; positronN++ ; mHistPidPt->Fill(1.5,pt) ;
    mHistMeanTPCPositron->Fill(lpTPC, dEdx) ;
    mHistMeanITSPositron->Fill(lpITS, its);
    mHistMeanTRDPositron->Fill(lpTRD, trd);
    mHistMeanTOFPositron->Fill(invMass, tof);
    mHistPidPositronPart->Fill(positron) ;
   }
   float muonPlus  = pFlowTrack->MuonPlusMinusProb() ;
   mHistPidMuonPlus->Fill(muonPlus) ;
   if(strcmp(pid, "mu+") == 0) 
   {
    numPid = 1 ; muonPlusN++ ; mHistPidPt->Fill(3.5,pt) ; 
    mHistMeanTPCMuonPlus->Fill(lpTPC, dEdx) ;
    mHistMeanITSMuonPlus->Fill(lpITS, its);
    mHistMeanTRDMuonPlus->Fill(lpTRD, trd);
    mHistMeanTOFMuonPlus->Fill(invMass, tof);
    mHistPidMuonPlusPart->Fill(muonPlus) ;
   }
   float piPlus = pFlowTrack->PionPlusMinusProb() ;
   mHistPidPiPlus->Fill(piPlus) ;
   if(strcmp(pid, "pi+") == 0) 
   { 
    numPid = 2 ; piPlusN++ ; mHistPidPt->Fill(5.5,pt) ;
    mHistMeanTPCPiPlus->Fill(lpTPC, dEdx) ;
    mHistMeanITSPiPlus->Fill(lpITS, its);
    mHistMeanTRDPiPlus->Fill(lpTRD, trd);
    mHistMeanTOFPiPlus->Fill(invMass, tof);
    mHistPidPiPlusPart->Fill(piPlus) ;
   }
   float kplus  = pFlowTrack->KaonPlusMinusProb() ;
   mHistPidKplus->Fill(kplus) ;
   if(strcmp(pid, "k+") == 0) 
   {
    numPid = 3 ; kPlusN++ ; mHistPidPt->Fill(7.5,pt) ; 
    mHistMeanTPCKplus->Fill(lpTPC, dEdx) ;
    mHistMeanITSKplus->Fill(lpITS, its);
    mHistMeanTRDKplus->Fill(lpTRD, trd);
    mHistMeanTOFKplus->Fill(invMass, tof);
    mHistPidKplusPart->Fill(kplus) ;
   }
   float proton  = pFlowTrack->ProtonPbarProb() ;
   mHistPidProton->Fill(proton) ;
   if(strcmp(pid, "pr+") == 0) 
   {
    numPid = 4 ; protonN++ ; mHistPidPt->Fill(9.5,pt) ;  
    mHistMeanTPCProton->Fill(lpTPC, dEdx) ;
    mHistMeanITSProton->Fill(lpITS, its);
    mHistMeanTRDProton->Fill(lpTRD, trd);
    mHistMeanTOFProton->Fill(invMass, tof);
    mHistPidProtonPart->Fill(proton) ;
   }
   float deuteron  = pFlowTrack->DeuteriumAntiDeuteriumProb() ;
   mHistPidDeuteron->Fill(deuteron) ;
   if(strcmp(pid, "d+") == 0) 
   {
    numPid = 6 ; deuteronN++ ; mHistPidPt->Fill(11.5,pt) ; 
    mHistMeanTPCDeuteron->Fill(lpTPC, dEdx) ;
    mHistMeanITSDeuteron->Fill(lpITS, its);
    mHistMeanTRDDeuteron->Fill(lpTRD, trd);
    mHistMeanTOFDeuteron->Fill(invMass, tof);
    mHistPidDeuteronPart->Fill(deuteron) ;
   }
  } 
  else if(charge == -1) 
  {
   hMinusN++ ;
   mHistMeanDedxNeg2D->Fill(lpTPC, dEdx) ;
   mHistMeanDedxNeg2DITS->Fill(lpITS, its) ;
   mHistMeanDedxNeg2DTRD->Fill(lpTRD, trd) ;
   mHistMeanDedxNeg2DTOF->Fill(lpTOF, tof) ;
   //
   float electron  = pFlowTrack->ElectronPositronProb() ;
   mHistPidElectron->Fill(electron);
   if(strcmp(pid, "e-") == 0) 
   {
    numPid = 0 ; electronN++ ; mHistPidPt->Fill(0.5,pt) ; 
    mHistMeanTPCElectron->Fill(lpTPC, dEdx);
    mHistMeanITSElectron->Fill(lpITS, its);
    mHistMeanTRDElectron->Fill(lpTRD, trd);
    mHistMeanTOFElectron->Fill(invMass, tof);
    mHistPidElectronPart->Fill(electron);
   }
   float muonMinus  = pFlowTrack->MuonPlusMinusProb() ;
   mHistPidMuonMinus->Fill(muonMinus) ;
   if(strcmp(pid, "mu-") == 0) 
   {
    numPid = 1 ; muonMinusN++ ; mHistPidPt->Fill(2.5,pt) ;
    mHistMeanTPCMuonMinus->Fill(lpTPC, dEdx) ;
    mHistMeanITSMuonMinus->Fill(lpITS, its);
    mHistMeanTRDMuonMinus->Fill(lpTRD, trd);
    mHistMeanTOFMuonMinus->Fill(invMass, tof);
    mHistPidMuonMinusPart->Fill(muonMinus) ;
   }
   float piMinus = pFlowTrack->PionPlusMinusProb() ;
   mHistPidPiMinus->Fill(piMinus) ;
   if(strcmp(pid, "pi-") == 0) 
   {
    numPid = 2 ; piMinusN++ ; mHistPidPt->Fill(4.5,pt) ;
    mHistMeanTPCPiMinus->Fill(lpTPC, dEdx);
    mHistMeanITSPiMinus->Fill(lpITS, its);
    mHistMeanTRDPiMinus->Fill(lpTRD, trd);
    mHistMeanTOFPiMinus->Fill(invMass, tof);
    mHistPidPiMinusPart->Fill(piMinus);
   }
   float kminus  = pFlowTrack->KaonPlusMinusProb() ;
   mHistPidKminus->Fill(kminus);
   if(strcmp(pid, "k-") == 0) 
   {
    numPid = 3 ; kMinusN++ ; mHistPidPt->Fill(6.5,pt) ;
    mHistMeanTPCKminus->Fill(lpTPC, dEdx);
    mHistMeanITSKminus->Fill(lpITS, its);
    mHistMeanTRDKminus->Fill(lpTRD, trd);
    mHistMeanTOFKminus->Fill(invMass, tof);
    mHistPidKminusPart->Fill(kminus);
   }
   float antiproton  = pFlowTrack->ProtonPbarProb() ;
   mHistPidAntiProton->Fill(antiproton);
   if(strcmp(pid, "pr-") == 0) 
   {
    numPid = 4 ; pbarN++ ; mHistPidPt->Fill(8.5,pt) ;
    mHistMeanTPCPbar->Fill(lpTPC, dEdx);
    mHistMeanITSPbar->Fill(lpITS, its);
    mHistMeanTRDPbar->Fill(lpTRD, trd);
    mHistMeanTOFPbar->Fill(invMass, tof);
    mHistPidAntiProtonPart->Fill(antiproton);
   }
   float antideuteron  = pFlowTrack->DeuteriumAntiDeuteriumProb() ;
   mHistPidAntiDeuteron->Fill(antideuteron);
   if(strcmp(pid, "d-") == 0) 
   {
    numPid = 6 ; dbarN++ ; mHistPidPt->Fill(10.5,pt) ;
    mHistMeanTPCAntiDeuteron->Fill(lpTPC, dEdx);
    mHistMeanITSAntiDeuteron->Fill(lpITS, its);
    mHistMeanTRDAntiDeuteron->Fill(lpTRD, trd);
    mHistMeanTOFAntiDeuteron->Fill(invMass, tof);
    mHistPidAntiDeuteronPart->Fill(antideuteron);
   }
  }
  
  // Yield3D, Yield2D, Eta, Pt, Phi, bayP.Id.
  mHistPtot->Fill(totalp) ;
  mHistPt->Fill(pt) ;
  mHistPhi->Fill(phi);  
  mHistAllEtaPtPhi3D->Fill(eta, pt, phi) ;
  mHistYieldAll2D->Fill(eta, pt) ;
  mHistBayPidMult->Fill(numPid) ;
  if(constrainable) 
  { 
   mHistPhiCons->Fill(phi);  
   mHistPhiPtCon->Fill(phi, pt);  
   mHistYieldCon2D->Fill(eta, pt) ;
   mHistConsEtaPtPhi3D->Fill(eta, pt, phi) ;
   mHistGlobEtaPtPhi3D->Fill(etaGlob, ptGlob, phiGlob) ;
  }
  else 
  { 
   mHistYieldUnc2D->Fill(etaGlob, ptGlob) ;
   mHistUncEtaPtPhi3D->Fill(etaGlob, ptGlob, phiGlob) ;
   mHistPhiPtUnc->Fill(phiGlob, ptGlob) ; 
  }
  if(pFlowTrack->Charge()>0)       { mHistPtPhiPos->Fill(phi, pt); }
  else if(pFlowTrack->Charge()<0)  { mHistPtPhiNeg->Fill(phi, pt); }

  // fills selected part histograms
  if(pFlowSelect->SelectPart(pFlowTrack)) 
  {
   if(strlen(pFlowSelect->PidPart()) != 0) 
   {
    float rapidity = pFlowTrack->Y();
    mHistBinEta->Fill(rapidity, rapidity);
    mHistYieldPart2D->Fill(rapidity, pt);
   }
   else 
   {
    mHistBinEta->Fill(eta, eta) ;
    mHistYieldPart2D->Fill(eta, pt) ;
   }
   mHistBayPidMultPart->Fill(numPid) ;
   mHistBinPt->Fill(pt, pt) ;
   mHistEtaPtPhi3DPart->Fill(eta,pt,phi) ;
   mHistDcaGlobalPart->Fill(dcaGlobal) ;
   mHistInvMassPart->Fill(invMass) ;
   if(charge == 1) 
   {
    mHistMeanDedxPos3DPart->Fill(lpITS, dEdx, numPid) ;
    mHistMeanDedxPos3DPartITS->Fill(lpITS, its, numPid) ;
   }
   else if(charge == -1) 
   {
    mHistMeanDedxNeg3DPart->Fill(lpITS, dEdx, numPid) ;
    mHistMeanDedxNeg3DPartITS->Fill(lpITS, its, numPid) ;
   }
   if(eta > 0.) { etaSymPosTpcNpart++ ; }  // for mHistEtaSymPart 
   else         { etaSymNegTpcNpart++ ; }
  }
  else         
  {
   mHistEtaPtPhi3DOut->Fill(eta,pt,phi) ;
   mHistYieldOut2D->Fill(eta, pt) ;
   mHistDcaGlobalOut->Fill(dcaGlobal) ;
   mHistInvMassOut->Fill(invMass) ;
  }

  // cos(n*phiLab)
  for(int j = 0; j < Flow::nHars; j++) 
  {
   bool oddHar = (j+1) % 2 ;
   float order = (float)(j+1) ;
   float vIn   = 100 * cos((double)order * phi) ;
   if(eta < 0 && oddHar) { vIn *= -1 ; }
   mHistCosPhi->Fill(order, vIn);
  }

  //For Eta symmetry TPC
  if(pFlowTrack->FitPtsTPC())
  {
   if(eta > 0.) { etaSymPosTpcN++ ; } // for mHistEtaSym 
   else         { etaSymNegTpcN++ ; }
  }
  
  // not to call it twice ...
  int nEtaS = HarmonicsLoop(eta,phi,pt,numPid) ;

  //For Correlated Multiplicity 
  corrMultN += (float)nEtaS ;

  //For Correlated Multiplicity in 1 unit rapidity
  if(TMath::Abs(eta) <= 0.5) { corrMultUnit += (float)nEtaS ; }

  //delete pointer
  pFlowTrack = 0 ;
 }      											 // end of tracks loop
						     
 // EtaSym
 float etaSymTpc = 0 ; 
 if(etaSymPosTpcN || etaSymNegTpcN) 	    { etaSymTpc  = (etaSymPosTpcN  - etaSymNegTpcN)  / (etaSymPosTpcN  + etaSymNegTpcN); }
 float etaSymTpcPart = 0 ; 
 if(etaSymPosTpcNpart || etaSymNegTpcNpart) { etaSymTpcPart  = (etaSymPosTpcNpart  - etaSymNegTpcNpart)  / (etaSymPosTpcNpart  + etaSymNegTpcNpart) ; }
 Float_t vertexZ = vertex[2] ;
 // all
 mHistEtaSym->Fill(etaSymTpc);
 mHistEtaSymVerZ2D->Fill(vertexZ , etaSymTpc);
 // selected
 mHistEtaSymPart->Fill(etaSymTpc);
 mHistEtaSymVerZ2DPart->Fill(vertexZ , etaSymTpcPart);

 // PID multiplicities
 float totalMult = (float)pFlowTracks->GetEntries() ;
 mHistPidMult->Fill(1., totalMult);
 mHistPidMult->Fill(2., hPlusN);
 mHistPidMult->Fill(3., hMinusN);
 mHistPidMult->Fill(4., piPlusN);
 mHistPidMult->Fill(5., piMinusN);
 mHistPidMult->Fill(6., protonN);
 mHistPidMult->Fill(7., pbarN);
 mHistPidMult->Fill(8., kPlusN);
 mHistPidMult->Fill(9., kMinusN);
 mHistPidMult->Fill(10., deuteronN);
 mHistPidMult->Fill(11., dbarN);
 mHistPidMult->Fill(12., electronN);
 mHistPidMult->Fill(13., positronN);
 mHistPidMult->Fill(14., muonMinusN);
 mHistPidMult->Fill(15., muonPlusN);

 // Multiplicity of particles correlated with the event planes
 corrMultN /= (float)(Flow::nHars * Flow::nSels) ; 
 mHistMultPart->Fill(corrMultN) ;
 // ...in one unit rapidity
 corrMultUnit /= (float)(Flow::nHars * Flow::nSels) ; 
 mHistMultPartUnit->Fill(corrMultUnit) ;

 if(Debug1) { cout << FlowAnalysis << "Tracks Loop...    " << totalMult << " tracks done . (" << corrMultN << ") . " << endl ; }
}
//-----------------------------------------------------------------------
int AliFlowAnalysisMaker::HarmonicsLoop(float eta, float phi, float pt, int numPid)
{
 // HarmonicsLoop ...

 if(Debug2) { cout << FlowAnalysis << "Harmonic Loop . " << endl ; }

 int corrMultN = 0 ;
 float zFirstPoint = pFlowTrack->ZFirstPoint() ; 
 // float zLastPoint = pFlowTrack->ZLastPoint() ;

 // Looping over Selections and Harmonics
 for (int k = 0; k < Flow::nSels; k++) 
 {
  pFlowSelect->SetSelection(k) ;
  for (int j = 0; j < Flow::nHars; j++) 
  {
   bool oddHar = (j+1) % 2;
   pFlowSelect->SetHarmonic(j);
   double order  = (double)(j+1);
   float psi_i, psi_2;
   if(pFlowEvent->EtaSubs())	       // particles with the opposite subevent
   {
    if(eta > 0) { psi_i = mPsiSub[1][k][j] ; }     //check
    else	{ psi_i = mPsiSub[0][k][j] ; }
   } 
   else if(order > 3. && !oddHar) 
   {
    psi_i = mPsi[k][1];  // 2nd harmomic event plane
    if(psi_i > 2*TMath::Pi()/order) { psi_i -= 2*TMath::Pi()/order ; } 
    if(psi_i > 2*TMath::Pi()/order) { psi_i -= 2*TMath::Pi()/order ; }
   } 
   else  // random subevents
   {
    psi_i = mPsi[k][j] ;
   }

   if(pFlowSelect->Select(pFlowTrack)) // Get detID
   {
    Bool_t kTpcPlus  = kFALSE ;
    Bool_t kTpcMinus = kFALSE ;
    Bool_t kTpcAll   = kFALSE ;

    histFull[k].histFullHar[j].mHistYieldPt->Fill(pt);
    histFull[k].histFullHar[j].mHistEtaPtPhi3D->Fill(eta, pt, phi);
    histFull[k].histFullHar[j].mHistYield2D->Fill(eta, pt);
    histFull[k].histFullHar[j].mHistDcaGlob->Fill(TMath::Abs(pFlowTrack->Dca()));

    // Set Tpc (+ and -)
    if(pFlowTrack->FitPtsTPC())    //OR*: AliESDtrack:: "TBits& GetTPCClusterMap()" or "Int_t GetTPCclusters(Int_t* idx)" ...
    {
     if(zFirstPoint >= 0. && eta > 0.)	    { kTpcPlus  = kTRUE ; } 
     else if(zFirstPoint <= 0. && eta < 0.) { kTpcMinus = kTRUE ; }
     else				    { kTpcAll   = kTRUE ; }	  
    }
    else 
    {
     if(Debug2) { cout << FlowAnalysis << "!!! no TPC hits ( track ??? ) ." << endl ; }
    }

    // PID Multiplicities (particle for R.P.) - done just one time for each selection
    if(j==0) { histFull[k].mHistBayPidMult->Fill(numPid) ; }

    // Calculate weights for filling histograms   
    float wt = 1. ; // TMath::Abs(pFlowEvent->Weight(k, j, pFlowTrack)) ; 

    // Fill histograms with selections
    if(kTpcPlus)       { histFull[k].histFullHar[j].mHistPhiPlus->Fill(phi,wt) ;  } 
    else if(kTpcMinus) { histFull[k].histFullHar[j].mHistPhiMinus->Fill(phi,wt) ; } 
    else if(kTpcAll)   { histFull[k].histFullHar[j].mHistPhiAll->Fill(phi,wt) ;   } 
    histFull[k].histFullHar[j].mHistPhi->Fill(phi,wt) ;

    // Get phiWgt from file
    double phiWgt;
    if(order > 3. && !oddHar) { phiWgt = pFlowEvent->PhiWeight(k, 1, pFlowTrack) ; } 
    else		      { phiWgt = pFlowEvent->PhiWeight(k, j, pFlowTrack) ; }
    if(oddHar && eta<0.)      { phiWgt /= -1. ; } // only for flat hists
    // chk
    if(Debug2) { cout << " Weight [" << k << "][" << j << "] (" << pFlowTrack->GetName() << ") = " << phiWgt << " or " << wt << " for phi = " << phi << endl ; }

    // Fill Flat histograms
    if(kTpcPlus)       { histFull[k].histFullHar[j].mHistPhiFlatPlus->Fill(phi, phiWgt) ; } 
    else if(kTpcMinus) { histFull[k].histFullHar[j].mHistPhiFlatMinus->Fill(phi, phiWgt) ; } 
    else if(kTpcAll)   { histFull[k].histFullHar[j].mHistPhiFlatAll->Fill(phi, phiWgt) ; } 
    histFull[k].histFullHar[j].mHistPhiFlat->Fill(phi,phiWgt) ;
    
    if(oddHar && eta<0.)  { phiWgt *= -1. ; } // restore value

    // Remove autocorrelations
    TVector2 Q_i;
    if(!pFlowEvent->EtaSubs())   // random subevents
    {
     if(order > 3. && !oddHar) // 2nd harmonic event plane
     { 
      Q_i.Set(phiWgt * cos(phi * 2), phiWgt * sin(phi * 2));
      TVector2 mQ_i = mQ[k][1] - Q_i;
      psi_i = mQ_i.Phi() / 2;
      if(psi_i < 0.) { psi_i += TMath::Pi() ; } 
     } 
     else 
     {
      Q_i.Set(phiWgt * cos(phi * order), phiWgt * sin(phi * order));
      TVector2 mQ_i = mQ[k][j] - Q_i;
      psi_i = mQ_i.Phi() / order;
      if(psi_i < 0.) { psi_i += 2*TMath::Pi()/order ; }
     }
    }
          
    // Remove autocorrelations of the second order 'particles' which are used for v1{EP1,EP2}.
    if (mV1Ep1Ep2 == kTRUE && order == 1) 
    {
     AliFlowSelection usedForPsi2 = *pFlowSelect ;
     usedForPsi2.SetHarmonic(1);
     if(usedForPsi2.Select(pFlowTrack))  // particle was used for Psi2
     {
      Q_i.Set(phiWgt * cos(phi * 2), phiWgt * sin(phi * 2));
      TVector2 mQ_i = mQ[k][1] - Q_i;
      psi_2 = mQ_i.Phi() / 2;
      if(psi_2 < 0.) { psi_2 += TMath::Pi() ; }
     }
     else				 // particle was not used for Psi2
     { 
      psi_2 = mPsi[k][1];
     }
    }
   }
   else
   {
    histFull[k].histFullHar[j].mHistYieldPtout->Fill(pt);
    histFull[k].histFullHar[j].mHistEtaPtPhi3Dout->Fill(eta, pt, phi);
    histFull[k].histFullHar[j].mHistYield2Dout->Fill(eta, pt);
    histFull[k].histFullHar[j].mHistDcaGlobout->Fill(TMath::Abs(pFlowTrack->Dca()));
   }

   // Caculate v for all particles selected for correlation analysis
   if(pFlowSelect->SelectPart(pFlowTrack)) 
   {
    corrMultN++;
    
    float v;
    if (mV1Ep1Ep2 == kFALSE || order != 1) 
    {
     v = 100 * cos(order * (phi - psi_i)) ;
    }
    else // i.e. (mV1Ep1Ep2 == kTRUE && order == 1)
    {
     v = 100 * cos(phi + psi_i - 2*psi_2) ;
    }
    
    float vFlip = v;
    if(eta < 0 && oddHar) { vFlip *= -1 ; }
    if(strlen(pFlowSelect->PidPart()) != 0) // pid, fill rapidity 
    {
     float rapidity = pFlowTrack->Y();
     histFull[k].histFullHar[j].mHist_vObs2D->Fill(rapidity, pt, v);
   
     if(mPtRange_for_vEta[1] > mPtRange_for_vEta[0])  // cut is used
     {
      if(pt < mPtRange_for_vEta[1] && pt >= mPtRange_for_vEta[0]) // check cut range, fill if in range
      {
       histFull[k].histFullHar[j].mHist_vObsEta->Fill(rapidity, v);
      }
     }
     else // cut is not used, fill in any case
     { 
      histFull[k].histFullHar[j].mHist_vObsEta->Fill(rapidity, v);
     }
    } 
    else  // no pid, fill eta
    {
     histFull[k].histFullHar[j].mHist_vObs2D->Fill(eta, pt, v);

     if(mPtRange_for_vEta[1] > mPtRange_for_vEta[0]) // cut is used
     {
      if(pt < mPtRange_for_vEta[1] && pt >= mPtRange_for_vEta[0]) // check cut range, fill if in range
      {
       histFull[k].histFullHar[j].mHist_vObsEta->Fill(eta, v);
      }
     }
     else // cut is not used, fill in any case
     { 
      histFull[k].histFullHar[j].mHist_vObsEta->Fill(eta, v);
     }
    }

    if(mEtaRange_for_vPt[1] > mEtaRange_for_vPt[0]) // cut is used
    { 
     if(TMath::Abs(eta) < mEtaRange_for_vPt[1] && TMath::Abs(eta) >= mEtaRange_for_vPt[0]) // check cut range, fill if in range
     {
      histFull[k].histFullHar[j].mHist_vObsPt->Fill(pt, vFlip);  // for odd harmonis /-/
     }
    }
    else  // cut is not used, fill in any case
    {
     histFull[k].histFullHar[j].mHist_vObsPt->Fill(pt, vFlip);
    }

    // v_
    Bool_t etaPtNoCut = kTRUE;
    if(mPtRange_for_vEta[1] > mPtRange_for_vEta[0] && (pt < mPtRange_for_vEta[0] || pt >= mPtRange_for_vEta[1])) 
    {
     etaPtNoCut = kFALSE;
    }
    if(mEtaRange_for_vPt[1] > mEtaRange_for_vPt[0] && (TMath::Abs(eta) < mEtaRange_for_vPt[0] || TMath::Abs(eta) >= mEtaRange_for_vPt[1]))
    {
     etaPtNoCut = kFALSE;
    }
    if(etaPtNoCut) { histFull[k].mHist_vObs->Fill(order, vFlip) ; }
 
    // Correlation of Phi of selected particles with Psi
    float phi_i = phi;
    if(eta < 0 && oddHar) 
    {
     phi_i += TMath::Pi() ; // backward particle and odd harmonic
     if(phi_i > 2*TMath::Pi()) { phi_i -= 2*TMath::Pi() ; }
    }
    float dPhi = phi_i - psi_i;
    if(dPhi < 0.)	       { dPhi += 2*TMath::Pi() ; }
    histFull[k].histFullHar[j].mHistPhiCorr->Fill(fmod((double)dPhi, 2*TMath::Pi() / order));
   }
  }
 }
 return corrMultN ;
}
//-----------------------------------------------------------------------
void AliFlowAnalysisMaker::FillV0Histograms()
{
 // v0s histograms

 if(Debug0) { cout << spaces << "V0s Loop . " << endl ; } 

 int corrMultV0 = 0 ;

 if(Debug1) { cout << spaces << "Looping over " << nV0s << "  flow V0s  . " << endl ; } 
 for(Int_t vloop=0;vloop<nV0s;vloop++) 
 {
  pFlowV0 = (AliFlowV0*)pFlowV0s->At(vloop) ;	    if(Debug3) { cout << "v0 n. " << vloop << endl ; pFlowV0->Dump() ; }

  // // LABELS FOR V0s ARE NOT DEFINED !!!
  // int label = pFlowV0->Label() ; 		    if(Debug3) { cout << vloop << " label: " << label << endl ; }

  float mass = pFlowV0->Mass() ;   	   	    if(Debug3) { cout << "mass " << mass << endl ; } 
  float eta = pFlowV0->Eta() ;   	   	    if(Debug3) { cout << "eta "  << eta << endl ; } 
  float rapidity = pFlowV0->Y() ;   	   	    if(Debug3) { cout << "Y   "  << rapidity << endl ; } 
  float phi = pFlowV0->Phi() ;     	   	    if(Debug3) { cout << "Phi "  << phi<< endl ; } 
  float pt = pFlowV0->Pt() ;      	   	    if(Debug3) { cout << "Pt "	<< pt << endl ; }  
  float totalp = pFlowV0->P() ;  	   	    if(Debug3) { cout << "P "	<< totalp << endl ; }	
  int	charge = pFlowV0->Charge() ;	    	    if(Debug3) { cout << "charge " << charge << endl ; }
  float dca = pFlowV0->Dca() ;		   	    if(Debug3) { cout << "dca "  << dca << endl ; }  
  float lenght = pFlowV0->V0Lenght() ;     	    if(Debug3) { cout << "lenght " << lenght << endl ; }
  float sigma = pFlowV0->Sigma() ; 	    	    if(Debug3) { cout << "sigma " << sigma << endl ; }  
  float chi2 = pFlowV0->Chi2() ; 	    	    if(Debug3) { cout << "chi2 " << chi2 << endl ; }  
  Char_t pid[10] ; strcpy(pid, pFlowV0->Pid()) ;    if(Debug3) { cout << "pid "  << pid << endl ; }   
  AliFlowTrack* daughterPlus = pFlowV0->DaughterP() ;   if(Debug3) { cout << "daughter + :  " << daughterPlus << endl ; }
  AliFlowTrack* daughterMinus = pFlowV0->DaughterN() ;  if(Debug3) { cout << "daughter - :  " << daughterMinus << endl ; }
  // -
  mHistV0Mass->Fill(mass) ;
  mHistV0EtaPtPhi3D->Fill(eta, pt, phi) ;
  mHistV0YieldAll2D->Fill(eta, pt) ;
  mHistV0Dca->Fill(dca);
  mHistV0Chi2->Fill(chi2);
  mHistV0Lenght->Fill(lenght);
  mHistV0Sigma->Fill(sigma);
  mHistV0MassPtSlices->Fill(mass,pt);
  if(Debug3) { cout << spaces <<  "global histograms filled . " << endl ; }

  if(pFlowSelect->SelectPart(pFlowV0))
  {
   bool inWin = pFlowSelect->SelectV0Part(pFlowV0) ;
   bool sx = pFlowSelect->SelectV0sxSide(pFlowV0) ;
   bool dx = pFlowSelect->SelectV0dxSide(pFlowV0) ;
   corrMultV0++ ;

   mHistV0YieldPart2D->Fill(eta, pt) ;
   if(inWin) 
   { 
    mHistV0EtaPtPhi3DPart->Fill(eta, pt, phi) ;
    mHistV0LenghtPart->Fill(lenght);
    mHistV0DcaPart->Fill(dca);
    mHistV0MassWin->Fill(mass) ; 
    mHistV0BinEta->Fill(eta, eta);
    mHistV0BinPt->Fill(pt, pt);   
   }
   else      
   { 
    mHistV0sbEtaPtPhi3DPart->Fill(eta, pt, phi) ;
    mHistV0sbLenghtPart->Fill(lenght);
    mHistV0sbDcaPart->Fill(dca);
    mHistV0sbMassSide->Fill(mass) ; 
    mHistV0sbBinEta->Fill(eta, eta);
    mHistV0sbBinPt->Fill(pt, pt);
   }

   for(int k = 0; k < Flow::nSels; k++)  // sort of HarmonicsLoop - selection number used
   {
    pFlowSelect->SetSelection(k) ;
    for(Int_t j=0;j<Flow::nHars;j++) 
    {
     Bool_t oddHar = (j+1) % 2 ;
     Float_t order = (Float_t)(j+1) ;
     pFlowSelect->SetHarmonic(j);			  if(Debug3) { cout << spaces <<  "start correlation computation .   " << j << endl ; }

     // Remove autocorrelations
     Float_t psi_i, psi_2 ;
     TVector2 Q_1, Q_2 ;
     Float_t phiDaughter1 = 0. ; Float_t phiDaughter2 = 0. ;
     Double_t phiWgt1 = 0. ; Double_t phiWgt2 = 0. ;
     // -
     if(daughterPlus)
     { // 1
      pFlowTrack = daughterPlus ;
      phiDaughter1 = pFlowTrack->Phi() ;		  if(Debug3) { cout << spaces <<  "track 1 . " << daughterPlus << " . phi1 = " << phiDaughter1 << " . " << endl ; } 
      if(pFlowSelect->Select(pFlowTrack)) // Get phiWgt from file
      {
       if(order > 3. && !oddHar) { phiWgt1 = pFlowEvent->PhiWeight(k, 1, pFlowTrack) ; } 
       else			 { phiWgt1 = pFlowEvent->PhiWeight(k, j, pFlowTrack) ; }
      }
      if(Debug3) { cout << spaces <<  "track 1 . " << daughterPlus << " . wgt1 = " << phiWgt1 << " . " << endl ; } 
     }
     if(daughterMinus)
     { // 2
      pFlowTrack = daughterMinus ;
      phiDaughter2 = pFlowTrack->Phi() ;		  if(Debug3) { cout << spaces <<  "track 2 . " << daughterMinus << " . phi2 = " << phiDaughter2 << " . " << endl ; }
      if(pFlowSelect->Select(pFlowTrack)) // Get phiWgt from file
      {
       if(order > 3. && !oddHar) { phiWgt2 = pFlowEvent->PhiWeight(k, 1, pFlowTrack) ; } 
       else			 { phiWgt2 = pFlowEvent->PhiWeight(k, j, pFlowTrack) ; }
      }
      if(Debug3) { cout << spaces <<  "track 2 . " << daughterMinus << " . wgt2 = " << phiWgt2 << " . " << endl ; } 
     }

     // psi_2
     Q_1.Set(phiWgt1 * cos(phiDaughter1 * 2), phiWgt1 * sin(phiDaughter1 * 2));
     Q_2.Set(phiWgt2 * cos(phiDaughter2 * 2), phiWgt2 * sin(phiDaughter2 * 2));
     TVector2 mQ_i = mQ[k][1] ; mQ_i -= Q_1 ; mQ_i -= Q_2 ; 
     psi_2 = mQ_i.Phi() / 2 ;	  
     if(psi_2 < 0.)	     { psi_2 += TMath::Pi() ; }
     if(psi_2 > TMath::Pi()) { psi_2 -= TMath::Pi() ; } 
     if(Debug3) { cout << spaces <<  "psi2 = " << psi_2 << " . " << endl ; }

     // psi_i
     if(order > 3. && !oddHar) { psi_i = psi_2 ; } 
     else 
     {
      Q_1.Set(phiWgt1 * cos(phiDaughter1 * order), phiWgt1 * sin(phiDaughter1 * order));
      Q_2.Set(phiWgt2 * cos(phiDaughter2 * order), phiWgt2 * sin(phiDaughter2 * order));
      TVector2 mQ_i = mQ[k][j] ; mQ_i -= Q_1 ; mQ_i -= Q_2 ;
      psi_i = mQ_i.Phi()/order ; 
      if(psi_i < 0.)		      { psi_i += 2*TMath::Pi()/order ; }	
      if(psi_i > 2*TMath::Pi()/order) { psi_i -= 2*TMath::Pi()/order ; } 
     }
     if(Debug3) { cout << spaces <<  "psi_i = " << psi_i << " . " << endl ; }

    // Caculate v for all V0s selected for correlation analysis
     float v ;
     if(mV1Ep1Ep2 == kFALSE || order != 1) { v = 100 * cos(order * (phi - psi_i)) ; }
     else { v = 100 * cos(phi + psi_i - 2*psi_2) ; }
     float vFlip = v ; if(eta < 0 && oddHar) { vFlip *= -1 ; }

     // invariant mass windows & sidebands
     if(inWin) { histFull[k].histFullHar[j].mHistV0_vObs2D->Fill(eta, pt, v) ; }
     else      { histFull[k].histFullHar[j].mHistV0sb_vObs2D->Fill(eta, pt, v) ; }

     if(mPtRange_for_vEta[1] > mPtRange_for_vEta[0]) // cut is used
     {
      if(pt < mPtRange_for_vEta[1] && pt >= mPtRange_for_vEta[0]) // check cut range, fill if in range
      {
       if(inWin) { histFull[k].histFullHar[j].mHistV0_vObsEta->Fill(eta, v) ; }
       else      
       { 
        histFull[k].histFullHar[j].mHistV0sb_vObsEta->Fill(eta, v) ; 
	if(sx)      { histFull[k].histFullHar[j].mHistV0sb_vObsEta_sx->Fill(eta, v) ; } 
	else if(dx) { histFull[k].histFullHar[j].mHistV0sb_vObsEta_dx->Fill(eta, v) ; } 
       }
      }
     }
     else // cut is not used, fill in any case
     { 
      if(inWin) { histFull[k].histFullHar[j].mHistV0_vObsEta->Fill(eta, v) ; }
      else	
      { 
       histFull[k].histFullHar[j].mHistV0sb_vObsEta->Fill(eta, v) ; 
       if(sx)	   { histFull[k].histFullHar[j].mHistV0sb_vObsEta_sx->Fill(eta, v) ; } 
       else if(dx) { histFull[k].histFullHar[j].mHistV0sb_vObsEta_dx->Fill(eta, v) ; } 
      }
     } 
     if(mEtaRange_for_vPt[1] > mEtaRange_for_vPt[0]) // cut is used
     {
      if(TMath::Abs(eta) < mEtaRange_for_vPt[1] && TMath::Abs(eta) >= mEtaRange_for_vPt[0]) // check cut range, fill if in range
      {          
       if(inWin) { histFull[k].histFullHar[j].mHistV0_vObsPt->Fill(pt, vFlip) ; }      // for odd harmonis /-/
       else      
       { 
        histFull[k].histFullHar[j].mHistV0sb_vObsPt->Fill(pt, vFlip) ; 
        if(sx)	    { histFull[k].histFullHar[j].mHistV0sb_vObsPt_sx->Fill(eta, v) ; } 
        else if(dx) { histFull[k].histFullHar[j].mHistV0sb_vObsPt_dx->Fill(eta, v) ; } 	
       }
      }
     }
     else  // cut is not used, fill in any case
     {
      if(inWin) { histFull[k].histFullHar[j].mHistV0_vObsPt->Fill(pt, vFlip) ; } 
      else	
      { 
       histFull[k].histFullHar[j].mHistV0sb_vObsPt->Fill(pt, vFlip) ; 
       if(sx)	   { histFull[k].histFullHar[j].mHistV0sb_vObsPt_sx->Fill(eta, v) ; } 
       else if(dx) { histFull[k].histFullHar[j].mHistV0sb_vObsPt_dx->Fill(eta, v) ; }  
      }
     }
     // v_
     Bool_t etaPtNoCut = kTRUE;
     if(mPtRange_for_vEta[1] > mPtRange_for_vEta[0] && (pt < mPtRange_for_vEta[0] || pt >= mPtRange_for_vEta[1])) 
     {
      etaPtNoCut = kFALSE;
     }
     if(mEtaRange_for_vPt[1] > mEtaRange_for_vPt[0] && (TMath::Abs(eta) < mEtaRange_for_vPt[0] || TMath::Abs(eta) >= mEtaRange_for_vPt[1]))
     {
      etaPtNoCut = kFALSE;
     }
     if(etaPtNoCut) 
     { 
      if(inWin) { histFull[k].mHistV0_vObs->Fill(order, vFlip) ; }
      else	
      { 
       if(sx)	   { histFull[k].mHistV0sb_vObs_sx->Fill(order, vFlip) ; } 
       else if(dx) { histFull[k].mHistV0sb_vObs_dx->Fill(order, vFlip) ; }  
      }
     }
   
     // Correlation of Phi of selected v0s with Psi
     float phi_i = phi;
     if(eta < 0 && oddHar)
     {
      phi_i += TMath::Pi() ; // backward particle and odd harmonic
      if(phi_i > 2*TMath::Pi()) { phi_i -= 2*TMath::Pi() ; }
     }
     float dPhi = phi_i - psi_i;
     if(dPhi < 0.)  { dPhi += 2*TMath::Pi() ; }

     if(inWin) { histFull[k].histFullHar[j].mHistV0PhiCorr->Fill(fmod((double)dPhi, 2*TMath::Pi() / order)) ; } 
     else      { histFull[k].histFullHar[j].mHistV0sbPhiCorr->Fill(fmod((double)dPhi, 2*TMath::Pi() / order)) ; }
    }
   }
  }
  //delete pFlowV0 ; pFlowV0 = 0 ;
  //delete pFlowTrack ; pFlowTrack = 0 ;
 }
 mHistV0MultPart->Fill(corrMultV0) ;

 if(Debug1) { cout << FlowAnalysis << "V0s Loop...    " << nV0s << " v0s done . (" << corrMultV0 << ") . " << endl ; }
}
//-----------------------------------------------------------------------
void AliFlowAnalysisMaker::FillLabels()  
{
 // Reads the tracks label (i.e. the link from ESD tracking to KineTree)
 // and puts them in a TH2F (nEvt x nLabels) for later use. The histogram 
 // is then  saved in the standar output.
 
 if(Debug1) { cout << FlowAnalysis << "Filling Labels . " << endl ; } 

// *temp*  --  !uncomment!
//  if(mLabellingFirst) 	
//  {
//   if(Debug1) { cout << spaces << " Creating Labels hist ( max label = " << maxLabel << " ) . " << endl ; }
//   mLabellingFirst = kFALSE ; 
//   TString* labTitle = new TString("TracksLabel_full") ; 
//   int maxMaxLabel = maxLabel + (int)((float)maxLabel / 10.) ;
//   mLabHist = new TH2F(labTitle->Data(),labTitle->Data(),Flow::nSels+4,0.,Flow::nSels+4.,maxMaxLabel+1,0.,maxMaxLabel+1.) ;
//   if(Flow::nSels==2) { mLabHist->SetXTitle("tracks , sel0 , sel1 , selPart , v0s , selV0s ") ; }
//   else  { mLabHist->SetXTitle("tracks , sel... , selPart , v0s , selV0s ") ; }
//   mLabHist->SetYTitle("MC label") ;
//   delete labTitle ;
//   if(Debug1) { cout << spaces << " Labels histogram created ( max bin = " << maxMaxLabel << " ) . " << endl ; }
//  }
// *temp*  --  !uncomment! 

 // Labels histogram
 float lBin = 0. ; float label = 0. ; 
 for(int tk=0;tk<nTracks;tk++) 
 {
  pFlowTrack = (AliFlowTrack*)pFlowTracks->At(tk) ;
  label = (float)pFlowTrack->Label() ; 				if(Debug4) { cout << " trak. " << tk  << " label. " << label  << endl ; }
  lBin = 0.1 ; 
  mLabHist->Fill(lBin,label) ;                           	if(Debug4) { cout << "    allTr	-	bin " << lBin << endl ; } 
  for(int ss=0;ss<Flow::nSels;ss++) 
  {
   lBin = (float)ss + 1.1 ;
   if(pFlowTrack->Select(0,ss,-1)) 	  
   { 
    mLabHist->Fill(lBin,label) ;                               	if(Debug4) { cout << "    sel " << ss << "	-	bin " << lBin << endl ; } 
   }
  }
  lBin = (float)Flow::nSels + 1.1 ;
  if(pFlowSelect->SelectPart(pFlowTrack)) 
  { 
   mLabHist->Fill(lBin,label) ;                              	if(Debug4) { cout << "    part	-	bin " << lBin << endl ; }  
  }
 }
 // LABELS FOR V0s ARE NOT DEFINED (...but for each daughter)
 if(mV0)
 {
  float label1 = 0. ; float label2 = 0. ; 
  for(int vk=0;vk<nV0s;vk++) 
  {
   pFlowV0 = (AliFlowV0*)pFlowV0s->At(vk) ;    		      if(Debug4) { cout << " v0  . " << vk  << " label from daughters .... " << endl ; }
   //label = (float)pFlowV0->Label() ;			      if(Debug4) { cout << " v0  . " << vk  << " label. " << label  << endl ; }
   pFlowTrack = (AliFlowTrack*)pFlowV0->DaughterP() ; label1 = (float)pFlowTrack->Label() ;
   pFlowTrack = (AliFlowTrack*)pFlowV0->DaughterN() ; label2 = (float)pFlowTrack->Label() ;
   					   		      if(Debug4) { cout << " daughter1 = " << label1  << " , daughter2 = " << label2  << endl ; }
   lBin = (float)Flow::nSels + 2.1 ; 
   mLabHist->Fill(lBin,label1) ; 			      if(Debug4) { cout << "	allV0 -       bin " << lBin << endl ; } 
   mLabHist->Fill(lBin,label2) ; 			      if(Debug4) { cout << "	allV0 -       bin " << lBin << endl ; } 
   lBin = (float)Flow::nSels + 3.1 ;   
   if(pFlowSelect->SelectPart(pFlowV0)) 
   {
    mLabHist->Fill(lBin,label1) ;		      	      if(Debug4) { cout << "	V0 part        -       bin " << lBin << endl ; }  
    mLabHist->Fill(lBin,label2) ;		      	      if(Debug4) { cout << "	V0 part        -       bin " << lBin << endl ; }  
   }
  }
 }
}
//----------------------------------------------------------------------
// ###
//----------------------------------------------------------------------
void AliFlowAnalysisMaker::Resolution()
{
 // Calculates the resolution, corrects the flow values, and stores all 
 // of them in few new histograms (see mVnResHistList). Histogram are 
 // then saved in the same otput file.
 
 if(Debug0) { cout << FlowAnalysis << "Calculating Resolution . " << endl ; cout << endl ; } 
 
 // VnRes histogram collection
 mVnResHistList = new TOrdCollection(Flow::nSels*Flow::nHars);
 
 // Calculate resolution from sqrt(mHistCos)
 double cosPair[Flow::nSels][Flow::nHars];
 double cosPairErr[Flow::nSels][Flow::nHars];
 double content, contentV0;
 double error, errorV0;
 double totalError;
 TString* histTitle;

 for(int k = 0; k < Flow::nSels; k++)
 {
  // v for tracks (corrected for resolution)
  histTitle = new TString("Flow_v_Sel");
  *histTitle += k+1;
  histFull[k].mHist_v = histFull[k].mHist_vObs->ProjectionX(histTitle->Data());
  histFull[k].mHist_v->SetTitle(histTitle->Data());
  histFull[k].mHist_v->SetXTitle("Harmonic");
  histFull[k].mHist_v->SetYTitle("v (%)");
  delete histTitle;
  mVnResHistList->AddLast(histFull[k].mHist_v);
   
  // v for v0s (corrected for resolution)
  histTitle = new TString("FlowV0_v_Sel");
  *histTitle += k+1;
  histFull[k].mHistV0_v = histFull[k].mHistV0_vObs->ProjectionX(histTitle->Data());
  histFull[k].mHistV0_v->SetTitle(histTitle->Data());
  histFull[k].mHistV0_v->SetXTitle("Harmonic");
  histFull[k].mHistV0_v->SetYTitle("v (%)");
  delete histTitle;
  mVnResHistList->AddLast(histFull[k].mHistV0_v);
   
  for(int j = 0; j < Flow::nHars; j++) 
  {
   double order = (double)(j+1);
   cosPair[k][j]    = histFull[k].mHistCos->GetBinContent(j+1);
   cosPairErr[k][j] = histFull[k].mHistCos->GetBinError(j+1);
   if(cosPair[k][j] > 0.) 
   {
    double resSub = 0. ;
    double resSubErr = 0. ;
    if(mV1Ep1Ep2 == kTRUE && order == 1) // calculate resolution of second order event plane first
    {
     double res2 = 0. ;
     double res2error = 0. ;
     if(histFull[k].mHistCos->GetBinContent(2) > 0.) 
     {
      if(mEtaSub)  // sub res only
      {
       res2 = TMath::Sqrt(histFull[k].mHistCos->GetBinContent(2));
       res2error = histFull[k].mHistCos->GetBinError(2) / (2. * res2);
      }
      else 
      {
       if (histFull[k].mHistCos->GetBinContent(2) > 0.92)  // resolution saturates
       {
        res2 = 0.99;
        res2error = 0.007;
       } 
       else 
       {
        double deltaRes2Sub = 0.005;  // differential for the error propagation
        double res2Sub = TMath::Sqrt(histFull[k].mHistCos->GetBinContent(2));
        double res2SubErr = histFull[k].mHistCos->GetBinError(2) / (2. * res2Sub);
        double chiSub2 = chi(res2Sub);
        double chiSub2Delta = chi(res2Sub + deltaRes2Sub);
        res2 = resEventPlane(TMath::Sqrt(2.) * chiSub2); // full event plane res.
        double mRes2Delta = resEventPlane(TMath::Sqrt(2.) * chiSub2Delta);
        res2error = res2SubErr * fabs((double)res2 - mRes2Delta) / deltaRes2Sub;
       }
      }
     }
     else 
     {
      res2 = 0.;
      res2error = 0.;
     }
     // now put everything together with first order event plane
     mRes[k][j]    = TMath::Sqrt(cosPair[k][0]*res2);
     mResErr[k][j] = 1./(2.*mRes[k][j]) * TMath::Sqrt(cosPairErr[k][0]*cosPairErr[k][0] + res2error*res2error); // Gaussian error propagation
    }
    else if(mEtaSub)  // sub res only
    {
     resSub = TMath::Sqrt(cosPair[k][j]);                
     resSubErr = cosPairErr[k][j] / (2. * resSub);
     mRes[k][j]    = resSub;
     mResErr[k][j] = resSubErr;
    } 
    else if(order==4. || order==6.|| order==8.)  // 2nd harmonic event plane
    {
     double deltaResSub = 0.005;  // differential for the error propagation
     double resSub = TMath::Sqrt(cosPair[k][1]);
     double resSubErr = cosPairErr[k][1] / (2. * resSub);
     double chiSub = chi(resSub);
     double chiSubDelta = chi(resSub + deltaResSub);
     double mResDelta;
     if (order==4.) 
     {
      mRes[k][j] = resEventPlaneK2(TMath::Sqrt(2.) * chiSub); // full event plane res.
      mResDelta = resEventPlaneK2(TMath::Sqrt(2.) * chiSubDelta);
     } 
     else if(order==6.) 
     {
      mRes[k][j] = resEventPlaneK3(TMath::Sqrt(2.) * chiSub); // full event plane res.
      mResDelta = resEventPlaneK3(TMath::Sqrt(2.) * chiSubDelta);
     } 
     else 
     {
      mRes[k][j] = resEventPlaneK4(TMath::Sqrt(2.) * chiSub); // full event plane res.
      mResDelta = resEventPlaneK4(TMath::Sqrt(2.) * chiSubDelta);
     }
     mResErr[k][j] = resSubErr * fabs((double)mRes[k][j] - mResDelta) / deltaResSub;
    }
    else 
    {
     if(cosPair[k][j] > 0.92) // resolution saturates
     {
      mRes[k][j]    = 0.99;
      mResErr[k][j] = 0.007;
     } 
     else 
     {
      double deltaResSub = 0.005;  // differential for the error propagation
      double resSub = TMath::Sqrt(cosPair[k][j]);
      double resSubErr = cosPairErr[k][j] / (2. * resSub);
      double chiSub = chi(resSub);
      double chiSubDelta = chi(resSub + deltaResSub);
      mRes[k][j] = resEventPlane(TMath::Sqrt(2.) * chiSub); // full event plane res.
      double mResDelta = resEventPlane(TMath::Sqrt(2.) * chiSubDelta);
      mResErr[k][j] = resSubErr * TMath::Abs((double)mRes[k][j] - mResDelta) / deltaResSub;
     }
    }
   }
   else  // subevent correlation must be positive
   {
    mRes[k][j]    = 0.;
    mResErr[k][j] = 0.;
   }
   histFull[k].mHistRes->SetBinContent(j+1, mRes[k][j]);
   histFull[k].mHistRes->SetBinError(j+1, mResErr[k][j]);

   // Create the v 2D histogram (Flow corrected for res.) - v,Pt,Eta
   histTitle = new TString("Flow_v2D_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   histFull[k].histFullHar[j].mHist_v2D = histFull[k].histFullHar[j].mHist_vObs2D->ProjectionXY(histTitle->Data());
   histFull[k].histFullHar[j].mHist_v2D->SetTitle(histTitle->Data());
   histFull[k].histFullHar[j].mHist_v2D->SetXTitle((char*)xLabel.Data());
   histFull[k].histFullHar[j].mHist_v2D->SetYTitle("Pt (GeV/c)");
   histFull[k].histFullHar[j].mHist_v2D->SetZTitle("v (%)");
   delete histTitle;
   mVnResHistList->AddLast(histFull[k].histFullHar[j].mHist_v2D);

   // Create the 1D v histograms - v,Eta
   histTitle = new TString("Flow_vEta_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   histFull[k].histFullHar[j].mHist_vEta = histFull[k].histFullHar[j].mHist_vObsEta->ProjectionX(histTitle->Data());
   histFull[k].histFullHar[j].mHist_vEta->SetTitle(histTitle->Data());
   histFull[k].histFullHar[j].mHist_vEta->SetXTitle((char*)xLabel.Data());
   histFull[k].histFullHar[j].mHist_vEta->SetYTitle("v (%)");
   delete histTitle;
   mVnResHistList->AddLast(histFull[k].histFullHar[j].mHist_vEta);
   
   // v,Pt
   histTitle = new TString("Flow_vPt_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   histFull[k].histFullHar[j].mHist_vPt = histFull[k].histFullHar[j].mHist_vObsPt->ProjectionX(histTitle->Data());
   histFull[k].histFullHar[j].mHist_vPt->SetTitle(histTitle->Data());
   histFull[k].histFullHar[j].mHist_vPt->SetXTitle("Pt (GeV/c)");
   histFull[k].histFullHar[j].mHist_vPt->SetYTitle("v (%)");
   delete histTitle;
   mVnResHistList->AddLast(histFull[k].histFullHar[j].mHist_vPt);

   // Create the v 2D histogram (V0s Flow corrected for res.) - v,Pt,Eta
   histTitle = new TString("FlowV0_v2D_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   histFull[k].histFullHar[j].mHistV0_v2D = histFull[k].histFullHar[j].mHistV0_vObs2D->ProjectionXY(histTitle->Data());
   histFull[k].histFullHar[j].mHistV0_v2D->SetTitle(histTitle->Data());
   histFull[k].histFullHar[j].mHistV0_v2D->SetXTitle((char*)xLabel.Data());
   histFull[k].histFullHar[j].mHistV0_v2D->SetYTitle("Pt (GeV/c)");
   histFull[k].histFullHar[j].mHistV0_v2D->SetZTitle("v (%)");
   delete histTitle;
   mVnResHistList->AddLast(histFull[k].histFullHar[j].mHistV0_v2D);

   // Create the 1D v histograms (V0s) - v,Eta
   histTitle = new TString("FlowV0_vEta_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   histFull[k].histFullHar[j].mHistV0_vEta = histFull[k].histFullHar[j].mHistV0_vObsEta->ProjectionX(histTitle->Data());
   histFull[k].histFullHar[j].mHistV0_vEta->SetTitle(histTitle->Data());
   histFull[k].histFullHar[j].mHistV0_vEta->SetXTitle((char*)xLabel.Data());
   histFull[k].histFullHar[j].mHistV0_vEta->SetYTitle("v (%)");
   delete histTitle;
   mVnResHistList->AddLast(histFull[k].histFullHar[j].mHistV0_vEta);
   
   // (V0s) v,Pt
   histTitle = new TString("FlowV0_vPt_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   histFull[k].histFullHar[j].mHistV0_vPt = histFull[k].histFullHar[j].mHistV0_vObsPt->ProjectionX(histTitle->Data());
   histFull[k].histFullHar[j].mHistV0_vPt->SetTitle(histTitle->Data());
   histFull[k].histFullHar[j].mHistV0_vPt->SetXTitle("Pt (GeV/c)");
   histFull[k].histFullHar[j].mHistV0_vPt->SetYTitle("v (%)");
   delete histTitle;
   mVnResHistList->AddLast(histFull[k].histFullHar[j].mHistV0_vPt);

   // Create the v 2D histogram (V0s sidebands Flow corrected for res.) - v,Pt,Eta
   histTitle = new TString("FlowV0sb_v2D_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   histFull[k].histFullHar[j].mHistV0sb_v2D = histFull[k].histFullHar[j].mHistV0sb_vObs2D->ProjectionXY(histTitle->Data());
   histFull[k].histFullHar[j].mHistV0sb_v2D->SetTitle(histTitle->Data());
   histFull[k].histFullHar[j].mHistV0sb_v2D->SetXTitle((char*)xLabel.Data());
   histFull[k].histFullHar[j].mHistV0sb_v2D->SetYTitle("Pt (GeV/c)");
   histFull[k].histFullHar[j].mHistV0sb_v2D->SetZTitle("v (%)");
   delete histTitle;
   mVnResHistList->AddLast(histFull[k].histFullHar[j].mHistV0sb_v2D);

   // Create the 1D v histograms (V0s sidebands) - v,Eta
   histTitle = new TString("FlowV0sb_vEta_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   histFull[k].histFullHar[j].mHistV0sb_vEta = histFull[k].histFullHar[j].mHistV0sb_vObsEta->ProjectionX(histTitle->Data());
   histFull[k].histFullHar[j].mHistV0sb_vEta->SetTitle(histTitle->Data());
   histFull[k].histFullHar[j].mHistV0sb_vEta->SetXTitle((char*)xLabel.Data());
   histFull[k].histFullHar[j].mHistV0sb_vEta->SetYTitle("v (%)");
   delete histTitle;
   mVnResHistList->AddLast(histFull[k].histFullHar[j].mHistV0sb_vEta);
   
   // (V0s sidebands) v,Pt
   histTitle = new TString("FlowV0sb_vPt_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   histFull[k].histFullHar[j].mHistV0sb_vPt = histFull[k].histFullHar[j].mHistV0sb_vObsPt->ProjectionX(histTitle->Data());
   histFull[k].histFullHar[j].mHistV0sb_vPt->SetTitle(histTitle->Data());
   histFull[k].histFullHar[j].mHistV0sb_vPt->SetXTitle("Pt (GeV/c)");
   histFull[k].histFullHar[j].mHistV0sb_vPt->SetYTitle("v (%)");
   delete histTitle;
   mVnResHistList->AddLast(histFull[k].histFullHar[j].mHistV0sb_vPt);

  // Calulate v = vObs / Resolution
   if(mRes[k][j]) 
   {
    cout << FlowAnalysis << "Resolution of the " << j+1 << "th harmonic = " << mRes[k][j] << " +/- " << mResErr[k][j] << endl;
    // selected tracks
    histFull[k].histFullHar[j].mHist_v2D->Scale(1. / mRes[k][j]);
    histFull[k].histFullHar[j].mHist_vEta->Scale(1. / mRes[k][j]);
    histFull[k].histFullHar[j].mHist_vPt->Scale(1. / mRes[k][j]);
    content = histFull[k].mHist_v->GetBinContent(j+1);        
    content /=  mRes[k][j]; 
    histFull[k].mHist_v->SetBinContent(j+1, content);            
    // selected V0s    
    histFull[k].histFullHar[j].mHistV0_v2D->Scale(1. / mRes[k][j]);
    histFull[k].histFullHar[j].mHistV0_vEta->Scale(1. / mRes[k][j]);
    histFull[k].histFullHar[j].mHistV0_vPt->Scale(1. / mRes[k][j]);
    contentV0 = histFull[k].mHistV0_v->GetBinContent(j+1);        
    contentV0 /=  mRes[k][j];
    histFull[k].mHistV0_v->SetBinContent(j+1, contentV0);            
    // V0s sidebands    
    histFull[k].histFullHar[j].mHistV0sb_v2D->Scale(1. / mRes[k][j]);
    histFull[k].histFullHar[j].mHistV0sb_vEta->Scale(1. / mRes[k][j]);
    histFull[k].histFullHar[j].mHistV0sb_vPt->Scale(1. / mRes[k][j]);
    // The systematic error of the resolution is folded in.
    // tracks
    error = histFull[k].mHist_v->GetBinError(j+1) ; 
    error /= mRes[k][j];                                         
    totalError = fabs(content) * TMath::Sqrt((error/content)*(error/content) + (mResErr[k][j]/mRes[k][j])*(mResErr[k][j]/mRes[k][j]));
    histFull[k].mHist_v->SetBinError(j+1, totalError);
    cout << FlowAnalysis << "v" << j+1 << "= (" << content << " +/- " << error << " +/- " << totalError << "(with syst.)) %" << endl;
    // V0s
    errorV0 = histFull[k].mHistV0_v->GetBinError(j+1) ; 
    errorV0 /= mRes[k][j]; 
    if (contentV0>0)
      totalError = fabs(contentV0) * TMath::Sqrt((errorV0/contentV0)*(errorV0/contentV0) + (mResErr[k][j]/mRes[k][j])*(mResErr[k][j]/mRes[k][j]));
    else
      totalError = 0;
    histFull[k].mHistV0_v->SetBinError(j+1, totalError);
   } 
   else 
   {
    cout << "##### Resolution of the " << j+1 << "th harmonic was zero." << endl;
    // tracks hist
    histFull[k].histFullHar[j].mHist_v2D->Reset();
    histFull[k].histFullHar[j].mHist_vEta->Reset();
    histFull[k].histFullHar[j].mHist_vPt->Reset();
    histFull[k].mHist_v->SetBinContent(j+1, 0.);
    histFull[k].mHist_v->SetBinError(j+1, 0.);
    // v0s hist
    histFull[k].histFullHar[j].mHistV0_v2D->Reset();
    histFull[k].histFullHar[j].mHistV0_vEta->Reset();
    histFull[k].histFullHar[j].mHistV0_vPt->Reset();
    histFull[k].histFullHar[j].mHistV0sb_v2D->Reset();
    histFull[k].histFullHar[j].mHistV0sb_vEta->Reset();
    histFull[k].histFullHar[j].mHistV0sb_vPt->Reset();
    histFull[k].mHistV0_v->SetBinContent(j+1, 0.);
    histFull[k].mHistV0_v->SetBinError(j+1, 0.);
   }
  }
 }
 if(Debug1) { cout << FlowAnalysis << "Done . " << endl ; cout << endl ; }

 return ;
}
//-----------------------------------------------------------------------
Double_t AliFlowAnalysisMaker::chi(double res)          // chi from the event plane resolution
{
 double chi   = 2.0;
 double delta = 1.0;
 for(int i = 0; i < 15; i++) 
 {
  if(resEventPlane(chi) < res) { chi = chi + delta ; }
  else                         { chi = chi - delta ; }
  delta = delta / 2.;
 }

 return chi;
}
//-----------------------------------------------------------------------
Double_t AliFlowAnalysisMaker::resEventPlane(double chi)    // plane resolution as function of chi
{
 double con = TMath::Sqrt(TMath::Pi()/2)/2 ;   // ~ 0.626657
 double arg = chi * chi / 4.;
 Double_t res = con * chi * exp(-arg) * (TMath::BesselI0(arg) + TMath::BesselI1(arg)); 

 return res ;
}
//-----------------------------------------------------------------------
Double_t AliFlowAnalysisMaker::resEventPlaneK2(double chi)  // ... for the case k=2.
{
 double con = TMath::Sqrt(TMath::Pi()/2)/2 ;   // ~ 0.626657
 double arg = chi * chi / 4.;

 double besselOneHalf = TMath::Sqrt(arg/(TMath::Pi()/2)) * sinh(arg)/arg;
 double besselThreeHalfs = TMath::Sqrt(arg/(TMath::Pi()/2)) * (cosh(arg)/arg - sinh(arg)/(arg*arg));
 Double_t res = con * chi * exp(-arg) * (besselOneHalf + besselThreeHalfs); 

 return res;
}
//-----------------------------------------------------------------------
Double_t AliFlowAnalysisMaker::resEventPlaneK3(double chi) // ...for the case k=3.
{
 double con = TMath::Sqrt(TMath::Pi()/2)/2 ;   // ~ 0.626657
 double arg = chi * chi / 4.;
 Double_t res = con * chi * exp(-arg) * (TMath::BesselI1(arg) + TMath::BesselI(2, arg)); 

 return res;
}
//-----------------------------------------------------------------------
Double_t AliFlowAnalysisMaker::resEventPlaneK4(double chi) // ...for the case k=4.
{
 double con = TMath::Sqrt(TMath::Pi()/2)/2 ;   // ~ 0.626657
 double arg = chi * chi / 4.;

 double besselOneHalf = TMath::Sqrt(arg/(TMath::Pi()/2)) * sinh(arg)/arg;
 double besselThreeHalfs = TMath::Sqrt(arg/(TMath::Pi()/2)) * (cosh(arg)/arg - sinh(arg)/(arg*arg));
 double besselFiveHalfs = besselOneHalf - 3*besselThreeHalfs/arg;
 Double_t res = con * chi * exp(-arg) * (besselThreeHalfs + besselFiveHalfs); 

 return res;
}
//-----------------------------------------------------------------------
Float_t AliFlowAnalysisMaker::Res(Int_t eventN, Int_t harN) 	const 	{ return mRes[eventN][harN]; }
Float_t AliFlowAnalysisMaker::ResErr(Int_t eventN, Int_t harN) 	const 	{ return mResErr[eventN][harN]; }
//-----------------------------------------------------------------------
// ###
//-----------------------------------------------------------------------
TString AliFlowAnalysisMaker::GetHistFileName() 		const	{ return mHistFileName ; }
TString AliFlowAnalysisMaker::GetPhiWgtFileName() 		const 	{ return mPhiWgtFileName ; }
TString AliFlowAnalysisMaker::GetInputFileName() 		const	{ return mFlowEvtFileName ; }
//-----------------------------------------------------------------------
void AliFlowAnalysisMaker::SetV1Ep1Ep2(Bool_t v1Ep1Ep2) 		{ mV1Ep1Ep2 = v1Ep1Ep2 ; }
void AliFlowAnalysisMaker::SetFlowForV0(Bool_t fV0) 			{ mV0 = fV0 ; }
void AliFlowAnalysisMaker::SetEtaSub() 					{ mEtaSub = kTRUE ; }
void AliFlowAnalysisMaker::SetHistFileName(TString name) 		{ mHistFileName = name ; }
void AliFlowAnalysisMaker::SetPhiWgtFileName(TString name) 		{ mPhiWgtFileName = name ; }
void AliFlowAnalysisMaker::SetInputFileName(TString name) 		{ mFlowEvtFileName = name ; }
void AliFlowAnalysisMaker::SetShuffle(Bool_t sh) 			{ mShuffle = sh ; }
void AliFlowAnalysisMaker::SetRedoWgt(Bool_t rd) 			{ mRedoWgt = rd ; }

void AliFlowAnalysisMaker::SetUsePhiWgt(Bool_t pw) 			{ mPhiWgt = pw ; }
void AliFlowAnalysisMaker::SetUseBayWgt(Bool_t bw) 			{ mBayWgt = bw ; }
void AliFlowAnalysisMaker::SetUsePtWgt(Bool_t ptw) 			{ mPtWgt = ptw ; }
void AliFlowAnalysisMaker::SetUseEtaWgt(Bool_t etw) 			{ mEtaWgt = etw ; }
void AliFlowAnalysisMaker::SetUseOnePhiWgt() 				{ mOnePhiWgt = kTRUE ; }
void AliFlowAnalysisMaker::SetUseFirstLastPhiWgt() 			{ mOnePhiWgt = kFALSE ; }

void AliFlowAnalysisMaker::SetFillLabels(Bool_t lb) 			{ mLabelling = lb ; }
void AliFlowAnalysisMaker::SetMaxLabel(Int_t bin) 			{ maxLabel = bin ; }
void AliFlowAnalysisMaker::SetOneInputFile(Bool_t in) 			{ fOneInputFile = in ; }
void AliFlowAnalysisMaker::SetMakeAll(Bool_t mka) 			{ mMakeAll = mka ; }
//-----------------------------------------------------------------------
void AliFlowAnalysisMaker::SetInputFileNames(TString name) 		
{ 
 TObjString* objName = new TObjString(name.Data()) ; 
 fFileNames.Add(objName) ;
}
//-----------------------------------------------------------------------
void AliFlowAnalysisMaker::SetHistoRanges(float etaMin, float etaMax, int EtaBins)    
{
 // sets eta limits for histogram 

 mEtaMin = etaMin ;
 mEtaMax = etaMax ;
 mNEtaBins = EtaBins ;
}
//------------------------------------------------------------------------
void AliFlowAnalysisMaker::SetPtRange_for_vEta(Float_t lo, Float_t hi)  
{
 // pt range for the v(eta) histograms

 mPtRange_for_vEta[0] = lo ;
 mPtRange_for_vEta[1] = hi ;
}
//------------------------------------------------------------------------
void AliFlowAnalysisMaker::SetEtaRange_for_vPt(Float_t lo, Float_t hi)
{
 // |eta| range for the v(pt) histograms 
 
 mEtaRange_for_vPt[0] = lo ;
 mEtaRange_for_vPt[1] = hi ;
}
//------------------------------------------------------------------------
void AliFlowAnalysisMaker::SetDebugg(Int_t db)  
{
 // debug options (0..3) 
 
 mDebug = db ; 
 Debug0 = kTRUE ; Debug1 = kFALSE ; Debug2 = kFALSE ; Debug3 = kFALSE ; 
 if(mDebug == 1) { Debug1 = kTRUE ; }
 if(mDebug == 2) { Debug1 = kTRUE ; Debug2 = kTRUE ; } 
 if(mDebug == 3) { Debug1 = kTRUE ; Debug3 = kTRUE ; } 
 if(mDebug == 4) { Debug1 = kTRUE ; Debug4 = kTRUE ; } 
 if(mDebug > 4)  { Debug1 = kTRUE ; Debug2 = kTRUE ; Debug3 = kTRUE ; Debug4 = kTRUE ; } 
}
//-----------------------------------------------------------------------

