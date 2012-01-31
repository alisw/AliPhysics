//////////////////////////////////////////////////////////////////////
//
// $Id: AliFlowAnalyser.cxx 18618 2007-05-16 15:38:22Z snelling $
//
// Author: Emanuele Simili
//
//////////////////////////////////////////////////////////////////////
//_____________________________________________________________
//
// Description: 
//         the AliFlowAnalyser provides the methods to perform an Event 
// plane flow analysis over AliFlowEvents. 
// - The method Init() generates the histograms.
// - The method Analyze(AliFlowEvent*) calculates/extracts flow observables,  
//  fills some histograms and performs the E.P. analysis.
// - The method Resolution() calculates the resolution correction and 
//  applies it to observed v_n values.
// - The method Finish() saves the histograms and closes the file.
// Weights calculation has been moved to the class AliFlowWeighter.
//
// The most important features of the Flow Analysis are:
//  - Reaction Plane calculation: for different harmonics and selections
//   with user defined cuts. (The calculation of the R.P. is a method in
//   the AliFlowEvent class).
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
//////////////////////////////////////////////////////////////////////

#ifndef ALIFLOWANALYSER_CXX
#define ALIFLOWANALYSER_CXX


// ROOT things
#include <TROOT.h>
#include <TMath.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TString.h>
#include <TObject.h>
#include <TClonesArray.h>
#include <TOrdCollection.h>
#include <TParticle.h>
#include <TParticlePDG.h>
#include <TDatabasePDG.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TVector2.h>

// Flow things
#include "AliFlowEvent.h"
#include "AliFlowTrack.h"
#include "AliFlowV0.h"
#include "AliFlowConstants.h"
#include "AliFlowSelection.h"
#include "AliFlowAnalyser.h"

// ANSI things
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
using namespace std; //required for resolving the 'cout' symbol

ClassImp(AliFlowAnalyser) 
//-----------------------------------------------------------------------
AliFlowAnalyser::AliFlowAnalyser(const AliFlowSelection* flowSelect):
  fTrackLoop(kTRUE), 
  fV0loop(kTRUE), 
  fShuffle(kFALSE), 
  fV1Ep1Ep2(kFALSE), 
  fEtaSub(kFALSE),
  fReadPhiWgt(kFALSE), 
  fBayWgt(kFALSE), 
  fRePid(kFALSE), 
  fCustomRespFunc(kFALSE), 
  fSub(0),
  fPtWgt(kFALSE),
  fEtaWgt(kFALSE), 
  fOnePhiWgt(kTRUE), 
  fHistFile(0x0), fPhiWgtFile(0x0),
  fHistFileName("flowAnalysPlot.root"), 
  fEventNumber(0), 
  fTrackNumber(0), 
  fV0Number(0), 
  fNumberOfTracks(0), 
  fNumberOfV0s(0),
  fPidId(0), 
  fSelParts(0), 
  fSelV0s(0),
  fTotalNumberOfEvents(0),
  fTotalNumberOfTracks(0),
  fTotalNumberOfV0s(0),
  fFlowEvent(0x0), 
  fFlowTrack(0x0), 
  fFlowV0(0x0), 
  fFlowSelect(0x0), 
  fFlowTracks(0x0), 
  fFlowV0s(0x0),  
  fPhiBins(AliFlowConstants::kPhiBins), 
  fPhiMin(0.), 
  fPhiMax(2*TMath::Pi()),
  fLabel(0), 
  fEtaMin(0.), 
  fEtaMax(0.), 
  fPtMin(0.), 
  fPtMax(0.), 
  fPtMaxPart(0),
  fEtaBins(0), 
  fPtBins(0), 
  fPtBinsPart(0), 
  fMaxLabel(0),
  fPhiWgtHistList(0x0), 
  fVnResHistList(0x0), 
  fHistTrigger(0),
  fHistMult(0), 
  fHistV0Mult(0), 
  fHistOrigMult(0), 
  fHistMultOverOrig(0), 
  fHistMultEta(0), 
  fHistCent(0), 
  fHistVertexZ(0), 
  fHistVertexXY2D(0), 
  fHistEnergyZDC(0), 
  fHistPartZDC(0), 
  fHistPidMult(0), 
  fHistBayPidMult(0), 
  fHistEtaSym(0), 
  fHistEtaSymPart(0), 
  fHistEtaSymVerZ2D(0), 
  fHistEtaSymVerZ2DPart(0),
  fHistMultPart(0),
  fHistV0MultPart(0),
  fHistBayPidMultPart(0),
  fHistMultPartUnit(0),
  fHistPtot(0),
  fHistPt(0),
  fHistCharge(0),
  fHistDcaGlobal(0),
  fHistDca(0),
  fHistTransDca(0),
  fHistChi2(0),
  fHistLenght(0),
  fHistInvMass(0),
  fHistFitOverMax(0),
  fHistPhiPtCon(0),
  fHistPhiPtUnc(0),
  fHistPtPhiPos(0),
  fHistPtPhiNeg(0),
  fHistAllEtaPtPhi3D(0),
  fHistCosPhi(0),
  fHistPidPt(0),
  fHistPhi(0),
  fHistPhiCons(0),
  fHistYieldAll2D(0),
  fHistYieldCon2D(0),
  fHistYieldUnc2D(0),
  fHistConsEtaPtPhi3D(0),
  fHistGlobEtaPtPhi3D(0),
  fHistUncEtaPtPhi3D(0),
  fHistChi2ITS(0),
  fHistChi2normITS(0),
  fHistFitPtsITS(0),
  fHistMaxPtsITS(0),
  fHistMeanDedxPos2DITS(0),
  fHistMeanDedxNeg2DITS(0),
  fHistChi2TPC(0),
  fHistChi2normTPC(0),
  fHistFitPtsTPC(0),
  fHistMaxPtsTPC(0),
  fHistFitOverMaxTPC(0),
  fHistMeanDedxPos2D(0),
  fHistMeanDedxNeg2D(0),
  fHistChi2TRD(0),
  fHistChi2normTRD(0),
  fHistFitPtsTRD(0),
  fHistMaxPtsTRD(0),
  fHistMeanDedxPos2DTRD(0),
  fHistMeanDedxNeg2DTRD(0),
  fHistChi2TOF(0),
  fHistChi2normTOF(0),
  fHistFitPtsTOF(0),
  fHistMaxPtsTOF(0),
  fHistMeanDedxPos2DTOF(0),
  fHistMeanDedxNeg2DTOF(0),
  fHistMeanTPCPiPlus(0),
  fHistMeanTPCPiMinus(0),
  fHistMeanTPCProton(0),
  fHistMeanTPCPbar(0),
  fHistMeanTPCKplus(0),
  fHistMeanTPCKminus(0),
  fHistMeanTPCDeuteron(0),
  fHistMeanTPCAntiDeuteron(0),
  fHistMeanTPCPositron(0),
  fHistMeanTPCElectron(0),
  fHistMeanTPCMuonPlus(0),
  fHistMeanTPCMuonMinus(0),
  fHistMeanITSPiPlus(0),
  fHistMeanITSPiMinus(0),
  fHistMeanITSProton(0),
  fHistMeanITSPbar(0),
  fHistMeanITSKplus(0),
  fHistMeanITSKminus(0),
  fHistMeanITSDeuteron(0),
  fHistMeanITSAntiDeuteron(0),
  fHistMeanITSPositron(0),
  fHistMeanITSElectron(0),
  fHistMeanITSMuonPlus(0),
  fHistMeanITSMuonMinus(0),
  fHistMeanTOFPiPlus(0),
  fHistMeanTOFPiMinus(0),
  fHistMeanTOFProton(0),
  fHistMeanTOFPbar(0),
  fHistMeanTOFKplus(0),
  fHistMeanTOFKminus(0),
  fHistMeanTOFDeuteron(0),
  fHistMeanTOFAntiDeuteron(0),
  fHistMeanTOFPositron(0),
  fHistMeanTOFElectron(0),
  fHistMeanTOFMuonPlus(0),
  fHistMeanTOFMuonMinus(0),
  fHistMeanTRDPiPlus(0),
  fHistMeanTRDPiMinus(0),
  fHistMeanTRDProton(0),
  fHistMeanTRDPbar(0),
  fHistMeanTRDKplus(0),
  fHistMeanTRDKminus(0),
  fHistMeanTRDDeuteron(0),
  fHistMeanTRDAntiDeuteron(0),
  fHistMeanTRDPositron(0),
  fHistMeanTRDElectron(0),
  fHistMeanTRDMuonPlus(0),
  fHistMeanTRDMuonMinus(0),
  fHistPidPiPlus(0),
  fHistPidPiMinus(0),
  fHistPidProton(0),
  fHistPidAntiProton(0),
  fHistPidKplus(0),
  fHistPidKminus(0),
  fHistPidDeuteron(0),
  fHistPidAntiDeuteron(0),
  fHistPidElectron(0),
  fHistPidPositron(0),
  fHistPidMuonMinus(0),
  fHistPidMuonPlus(0),
  fHistPidPiPlusPart(0),
  fHistPidPiMinusPart(0),
  fHistPidProtonPart(0),
  fHistPidAntiProtonPart(0),
  fHistPidKplusPart(0),
  fHistPidKminusPart(0),
  fHistPidDeuteronPart(0),
  fHistPidAntiDeuteronPart(0),
  fHistPidElectronPart(0),
  fHistPidPositronPart(0),
  fHistPidMuonMinusPart(0),
  fHistPidMuonPlusPart(0),
  fHistBinEta(0),
  fHistBinPt(0),
  fHistEtaPtPhi3DPart(0),
  fHistYieldPart2D(0),
  fHistDcaGlobalPart(0),
  fHistInvMassPart(0),
  fHistEtaPtPhi3DOut(0),
  fHistYieldOut2D(0),
  fHistDcaGlobalOut(0),
  fHistInvMassOut(0),
  fHistMeanDedxPos3DPart(0),
  fHistMeanDedxNeg3DPart(0),
  fHistMeanDedxPos3DPartITS(0),
  fHistMeanDedxNeg3DPartITS(0),
  fHistV0Mass(0),
  fHistV0EtaPtPhi3D(0),
  fHistV0YieldAll2D(0),
  fHistV0PYall2D(0),
  fHistV0Dca(0),
  fHistV0Chi2(0),
  fHistV0Lenght(0),
  fHistV0Sigma(0),
  fHistV0CosPhi(0),
  fHistV0MassPtSlices(0),
  fHistV0BinEta(0),
  fHistV0BinPt(0),
  fHistV0sbBinEta(0),
  fHistV0sbBinPt(0),
  fHistV0MassWin(0),
  fHistV0EtaPtPhi3DPart(0),
  fHistV0YieldPart2D(0),
  fHistV0DcaPart(0),
  fHistV0LenghtPart(0),
  fHistV0sbMassSide(0),
  fHistV0sbEtaPtPhi3DPart(0),
  fHistV0sbYieldPart2D(0),
  fHistV0sbDcaPart(0),
  fHistV0sbLenghtPart(0)
{
 // default constructor (selection given or default selection)

 if(flowSelect) { fFlowSelect = new AliFlowSelection(*flowSelect) ; }
 else 		{ fFlowSelect = new AliFlowSelection() ; }
}
//-----------------------------------------------------------------------
AliFlowAnalyser::~AliFlowAnalyser()
{
 // default distructor (no actions) 
}
//-----------------------------------------------------------------------
Bool_t AliFlowAnalyser::Init() 
{
 // sets some defaults for the analysis
 
 cout << "* FlowAnalysis *  -  Init()" << endl ; cout << endl ; 
 
 // Open output files (->plots)
 fHistFile = new TFile(fHistFileName.Data(), "RECREATE");
 fHistFile->cd() ;

 // counters and pointers to 0
 fEventNumber = 0 ;
 fTrackNumber = 0 ; 	
 fV0Number    = 0 ;
 fPidId       = -1 ;
 //fNumberOfEvents = 0 ;	
 fNumberOfTracks = 0 ;	
 fNumberOfV0s    = 0 ; 	
 fFlowEvent  = 0 ;
 fFlowTrack  = 0 ;
 fFlowV0     = 0 ;
 fFlowTracks = 0 ;
 fFlowV0s    = 0 ;
 fSelParts = 0 ;
 fSelV0s   = 0 ;
 for(Int_t ii=0;ii<3;ii++) { fVertex[ii] = 0 ; }  
 
 fTotalNumberOfEvents = 0 ;
 fTotalNumberOfTracks = 0 ;
 fTotalNumberOfV0s = 0 ;

 // Check for Reading weights: if weight file is not there, wgt arrays are not used (phiwgt & bayesian)
 if(!fPhiWgtFile) 
 { 
  fReadPhiWgt = kFALSE ; 
  fBayWgt     = kFALSE ;
  cout << "  No wgt file. Phi & Bayesian weights will not be used ! " << endl ;
 }
 
 // Histogram settings
 fLabel = "Pseudorapidity" ;  			if(strlen(fFlowSelect->PidPart()) != 0) 	   { fLabel = "Rapidity"; }
 fPtMaxPart = AliFlowConstants::fgPtMaxPart ;  	if(fFlowSelect->PtMaxPart())    { fPtMaxPart = fFlowSelect->PtMaxPart() ; }
 fPtBinsPart = AliFlowConstants::kPtBinsPart ;  if(fFlowSelect->PtBinsPart()) { fPtBinsPart = fFlowSelect->PtBinsPart() ; }
 fPtMin = AliFlowConstants::fgPtMin ;       
 fPtMax = AliFlowConstants::fgPtMax ;       
 fEtaMin = AliFlowConstants::fgEtaMin ;
 fEtaMax = AliFlowConstants::fgEtaMax ;
 fPtBins  = AliFlowConstants::kPtBins ;      
 fEtaBins = AliFlowConstants::kEtaBins ;      
 // -
 SetPtRangevEta(1.,0.);   // (AliFlowConstants::fgPtMin , AliFlowConstants::fgPtMax) ;
 SetEtaRangevPt(1.,0.);   // (AliFlowConstants::fgEtaMin , AliFlowConstants::fgEtaMax) ;
 // -
 float triggerMin      =  -1.5;
 float triggerMax      =  10.5;
 float chargeMin       =  -2.5;
 float chargeMax       =   2.5; 
 float dcaMin	       =   0.0;
 float dcaMax	       =  10.0; //  0.5; 
 float glDcaMax        =  50.0; 
 float chi2Min         =    0.;
 float chi2Max         =  500.; 
 float chi2normMax     =   50.; 
 float chi2MaxC        =   30.; 
 float chi2MaxI        =   60.; 
 float fitPtsMin       =  -0.5;
 float fitPtsMaxTpc    = 160.5;
 float fitPtsMaxIts    =   6.5;
 float fitPtsMaxTrd    = 130.5;
 float fitPtsMaxTof    =   1.5;
 float fitOverMaxMin   =    0.;
 float fitOverMaxMax   =   1.1;
 float multOverOrigMin =    0.;
 float multOverOrigMax =    1.;
 float vertexZMin      =  -10.;
 float vertexZMax      =   10.; 
 float vertexXYMin     =  -0.1;
 float vertexXYMax     =   0.1; 
 float etaSymMinPart   =   -1.;
 float etaSymMaxPart   =    1.;
 float etaSymMin       =   -2.;
 float etaSymMax       =    2.;
 float psiMin	       =    0.;
 float psiMax	       = 2*TMath::Pi() ; 
 float multMin         =    0. ;
 float multMax         =25000. ;
 float multV0	       = 5000. ;
 float qMin	       =    0. ;
 float qMax	       =   20. ;
 float pidMin	       =    0. ;
 float pidMax	       =    1. ;
 float centMin         =  -0.5 ;
 float centMax         =   9.5 ;
 float logpMin         =  -2.5 ;
 float logpMax         =   2.5 ;
 float pMin	       =    0. ;
 float pMax	       =   10. ;
 float dEdxMax         = 1000. ;
 float dEdxMaxTPC      = 1500. ;
 float tofmin	       = 10000. ;
 float tofmax	       = 50000. ;
 float trdmax	       = 2000. ;
 float lgMin	       = -1000. ;
 float lgMax	       = 1000. ;
 float lgMinV0         =    0. ;
 float lgMaxV0         =   50. ;
 float massMin         =  -0.1 ;
 float massMax         =   2.1 ;
 float zdcpartMin      =  -0.5 ;
 float zdcpartMax      = 1000.5 ;
 float zdceMin         =    0. ;
 float zdceMax         = 10000. ;
 // - 
 enum { kTriggerBins	  = 12,
	kChargeBins	  = 5,
	kDcaBins	  = 1000,
	kChi2Bins	  = 150,
	kFitPtsBinsTpc    = 161,
	kFitPtsBinsIts    = 7,
	kFitPtsBinsTrd    = 131,
	kFitPtsBinsTof    = 2,
	kFitOverMaxBins   = 55,
	kMultOverOrigBins =100,
	kVertexZBins	  = 200,
	kVertexXYBins	  = 1000,
	kEtaSymBins	  = 200,
	kPhi3DBins	  = 60,
	kPsiBins	  = 36,
	kMultBins	  = 250,
	kPidBins	  = 100,
        kCentBins	  = 10,
	kDedxBins	  = 1000,
	kTofBins	  = 1000,
	kTrdBins	  = 100,
	kMomenBins	  = 100,
	kQbins 	 	  =  50,
	kLgBins 	  = 500,
	kMassBins	  = 2200,
 	kZDCpartBins	  = 1001,
 	kZDCeBins	  = 200
      } ;

 // Histograms Booking ...
 // Trigger
 fHistTrigger = new TH1F("Flow_Trigger", "Flow_Trigger", kTriggerBins, triggerMin, triggerMax);
 fHistTrigger->SetXTitle("Trigger: 1 minbias, 2 central, 3 laser, 10 other");
 fHistTrigger->SetYTitle("Counts");

 // Charge
 fHistCharge = new TH1F("Flow_Charge", "Flow_Charge", kChargeBins, chargeMin, chargeMax);
 fHistCharge->SetXTitle("Charge");
 fHistCharge->SetYTitle("Counts");

 // Reconstructed Momentum (constrained, if so)
 fHistPtot = new TH1F("Flow_totalP","Flow_totalP", fPtBins, fPtMin, fPtMax);
 fHistPtot->Sumw2();
 fHistPtot->SetXTitle("P (GeV/c)");
 fHistPtot->SetYTitle("Counts");

 // Reconstructed pT (constrained, if so)
 fHistPt = new TH1F("Flow_Pt","Flow_Pt", kMomenBins, pMin, pMax);
 fHistPt->Sumw2();
 fHistPt->SetXTitle("Pt (GeV/c)");
 fHistPt->SetYTitle("Counts");
   
 // Reconstructed P.id vs pT 
 fHistPidPt = new TH2F("Flow_PidPt","Flow_PidPt", 12, 0., 12., fPtBins, fPtMin, fPtMax);
 fHistPidPt->Sumw2();
 fHistPidPt->SetXTitle("e-, e+, mu-, mu+, pi-, pi+, K-, K+, pr-, pr+, d-, d+");
 fHistPidPt->SetYTitle("Pt (GeV/c)");
   
 // Distance of closest approach - Transverse, Unsigned
 fHistDca = new TH1F("Flow_TransDCA", "Flow_TransverseDCA", kDcaBins, dcaMin, dcaMax);
 fHistDca->Sumw2();
 fHistDca->SetXTitle("|Track's Signed dca (cm)|");
 fHistDca->SetYTitle("Counts");

 // Distance of closest approach - Transverse
 fHistTransDca = new TH1F("Flow_TransDCA_Signed", "Flow_TransverseDCA_Signed", kDcaBins, -dcaMax, dcaMax);
 fHistTransDca->Sumw2();
 fHistTransDca->SetXTitle("Track's Transverse dca (cm)");
 fHistTransDca->SetYTitle("Counts");

 // Distance of closest approach - 3d
 fHistDcaGlobal = new TH1F("Flow_3dDcaGlobal", "Flow_3dDcaGlobal", kDcaBins, dcaMin, glDcaMax);
 fHistDcaGlobal->Sumw2();
 fHistDcaGlobal->SetXTitle("|3d Global Track's dca to Vertex (cm)|");
 fHistDcaGlobal->SetYTitle("Counts");

 // Distance of closest approach for particles correlated with the event plane - 3d
 fHistDcaGlobalPart = new TH1F("Flow_3dDcaGlobalPart", "Flow_3dDcaGlobalPart", kDcaBins, dcaMin, glDcaMax);
 fHistDcaGlobalPart->Sumw2();
 fHistDcaGlobalPart->SetXTitle("|3d Global Track's dca to Vertex (cm)|");
 fHistDcaGlobalPart->SetYTitle("Counts");

 // Distance of closest approach for particles NOT SELECTED for correlation with the event plane - 3d
 fHistDcaGlobalOut = new TH1F("Flow_3dDcaGlobalOut", "Flow_3dDcaGlobalOut", kDcaBins, dcaMin, glDcaMax);
 fHistDcaGlobalOut->Sumw2();
 fHistDcaGlobalOut->SetXTitle("|3d Global Track's dca to Vertex (cm)|");
 fHistDcaGlobalOut->SetYTitle("Counts");

 // Yield Pt-Phi for constrained tracks
 fHistPhiPtCon = new TH2D("Flow_PtPhi_Con", "Flow_PtPhi_Con", fPhiBins, fPhiMin, fPhiMax, fPtBins, fPtMin, fPtMax);
 fHistPhiPtCon->Sumw2();
 fHistPhiPtCon->SetXTitle("Phi");
 fHistPhiPtCon->SetYTitle("Pt (GeV/c)");
 
 // Yield Pt-Phi for UNconstrained tracks
 fHistPhiPtUnc = new TH2D("Flow_PtPhi_Unc", "Flow_PtPhi_Unc", fPhiBins, fPhiMin, fPhiMax, fPtBins, fPtMin, fPtMax);
 fHistPhiPtUnc->Sumw2();
 fHistPhiPtUnc->SetXTitle("Phi");
 fHistPhiPtUnc->SetYTitle("Pt (GeV/c)");
 
 // Chi2 
 // at the main vertex
 fHistChi2 = new TH1F("Flow_Chi2", "Flow_Chi2", kChi2Bins, chi2Min, chi2MaxC);
 fHistChi2->SetXTitle("Chi square at Main Vertex");
 fHistChi2->SetYTitle("Counts");
 // TPC
 fHistChi2TPC = new TH1F("Flow_Chi2_TPC", "Flow_Chi2_TPC", kChi2Bins, chi2Min, chi2Max);
 fHistChi2TPC->SetXTitle("Chi square for TPC");
 fHistChi2TPC->SetYTitle("Counts");
 // ITS
 fHistChi2ITS = new TH1F("Flow_Chi2_ITS", "Flow_Chi2_ITS", kChi2Bins, chi2Min, chi2MaxI);
 fHistChi2ITS->SetXTitle("Chi square for ITS");
 fHistChi2ITS->SetYTitle("Counts");
 // TRD
 fHistChi2TRD = new TH1F("Flow_Chi2_TRD", "Flow_Chi2_TRD", kChi2Bins, chi2Min, chi2Max);
 fHistChi2TRD->SetXTitle("Chi square for TRD");
 fHistChi2TRD->SetYTitle("Counts");
 // TOF
 fHistChi2TOF = new TH1F("Flow_Chi2_TOF", "Flow_Chi2_TOF", kChi2Bins, chi2Min, chi2Max);
 fHistChi2TOF->SetXTitle("Chi square for TOF");
 fHistChi2TOF->SetYTitle("Counts");

 // Chi2 normalized
 // TPC
 fHistChi2normTPC = new TH1F("Flow_Chi2norm_TPC", "Flow_Chi2norm_TPC", kChi2Bins, chi2Min, chi2normMax);
 fHistChi2normTPC->SetXTitle("Normalized #Chi^{2} for TPC");
 fHistChi2normTPC->SetYTitle("Counts");
 // ITS
 fHistChi2normITS = new TH1F("Flow_Chi2norm_ITS", "Flow_Chi2norm_ITS", kChi2Bins, chi2Min, chi2normMax);
 fHistChi2normITS->SetXTitle("Normalized #Chi^{2} for ITS");
 fHistChi2normITS->SetYTitle("Counts");
 // TRD
 fHistChi2normTRD = new TH1F("Flow_Chi2norm_TRD", "Flow_Chi2norm_TRD", kChi2Bins, chi2Min, chi2normMax);
 fHistChi2normTRD->SetXTitle("Normalized #Chi^{2} for TRD");
 fHistChi2normTRD->SetYTitle("Counts");
 // TOF
 fHistChi2normTOF = new TH1F("Flow_Chi2norm_TOF", "Flow_Chi2norm_TOF", kChi2Bins, chi2Min, chi2normMax);
 fHistChi2normTOF->SetXTitle("Normalized #Chi^{2} for TOF");
 fHistChi2normTOF->SetYTitle("Counts");
 
 // FitPts
 // TPC
 fHistFitPtsTPC = new TH1F("Flow_FitPts_TPC", "Flow_FitPts_TPC", kFitPtsBinsTpc, fitPtsMin, fitPtsMaxTpc);
 fHistFitPtsTPC->SetXTitle("Fit Points");
 fHistFitPtsTPC->SetYTitle("Counts");
 // ITS
 fHistFitPtsITS = new TH1F("Flow_HitPts_ITS", "Flow_HitPts_ITS", kFitPtsBinsIts, fitPtsMin, fitPtsMaxIts);
 fHistFitPtsITS->SetXTitle("Fit Points");
 fHistFitPtsITS->SetYTitle("Counts");
 // TRD
 fHistFitPtsTRD = new TH1F("Flow_FitPts_TRD", "Flow_FitPts_TRD", kFitPtsBinsTrd, fitPtsMin, fitPtsMaxTrd);
 fHistFitPtsTRD->SetXTitle("Fit Points");
 fHistFitPtsTRD->SetYTitle("Counts");
 // TOF
 fHistFitPtsTOF = new TH1F("Flow_HitPts_TOF", "Flow_HitPts_TOF", kFitPtsBinsTof, fitPtsMin, fitPtsMaxTof);
 fHistFitPtsTOF->SetXTitle("Fit Points");
 fHistFitPtsTOF->SetYTitle("Counts");

 // MaxPts 
 // TPC
 fHistMaxPtsTPC = new TH1F("Flow_MaxPts_TPC", "Flow_MaxPts_TPC", kFitPtsBinsTpc, fitPtsMin, fitPtsMaxTpc);
 fHistMaxPtsTPC->SetXTitle("Max Points");
 fHistMaxPtsTPC->SetYTitle("Counts");
 // ITS
 fHistMaxPtsITS = new TH1F("Flow_MaxPts_ITS", "Flow_MaxPts_ITS", kFitPtsBinsIts, fitPtsMin, fitPtsMaxIts);
 fHistMaxPtsITS->SetXTitle("Max Points");
 fHistMaxPtsITS->SetYTitle("Counts");
 // TRD
 fHistMaxPtsTRD = new TH1F("Flow_MaxPts_TRD", "Flow_MaxPts_TRD", kFitPtsBinsTrd, fitPtsMin, fitPtsMaxTrd);
 fHistMaxPtsTRD->SetXTitle("Max Points");
 fHistMaxPtsTRD->SetYTitle("Counts");
 // TOF
 fHistMaxPtsTOF = new TH1F("Flow_MaxPts_TOF", "Flow_MaxPts_TOF", kFitPtsBinsTof, fitPtsMin, fitPtsMaxTof);
 fHistMaxPtsTOF->SetXTitle("Max Points");
 fHistMaxPtsTOF->SetYTitle("Counts");

 // FitOverMax 
 // Tpc
 fHistFitOverMaxTPC = new TH1F("Flow_FitOverMax_TPC", "Flow_FitOverMax_TPC", kFitOverMaxBins, fitOverMaxMin, fitOverMaxMax);
 fHistFitOverMaxTPC->SetXTitle("(Fit Points - 1) / Max Points");
 fHistFitOverMaxTPC->SetYTitle("Counts");  
 // All
 fHistFitOverMax = new TH1F("Flow_FitOverMax", "Flow_FitOverMax", kFitOverMaxBins, fitOverMaxMin, fitOverMaxMax);
 fHistFitOverMax->SetXTitle("(Fit Points - 1) / Max Points");
 fHistFitOverMax->SetYTitle("Counts");  

 // lenght
 fHistLenght = new TH1F("Flow_TrackLenght", "Flow_TrackLenght", kLgBins, lgMin, lgMax);
 fHistLenght->SetXTitle("Lenght of the Track (cm)");
 fHistLenght->SetYTitle("Counts");

  // OrigMult
 fHistOrigMult = new TH1F("Flow_OrigMult", "Flow_OrigMult", kMultBins, multMin, multMax);
 fHistOrigMult->SetXTitle("Original Mult");
 fHistOrigMult->SetYTitle("Counts");

 // MultEta
 fHistMultEta = new TH1F("Flow_MultEta", "Flow_MultEta", kMultBins, multMin, multMax);
 fHistMultEta->SetXTitle("Mult for Centrality");
 fHistMultEta->SetYTitle("Counts");
   
 // Mult
 fHistMult = new TH1F("Flow_Mult", "Flow_Mult", kMultBins, multMin, multMax);
 fHistMult->SetXTitle("Mult");
 fHistMult->SetYTitle("Counts");

 // V0s multiplicity
 fHistV0Mult = new TH1F("FlowV0_Mult","FlowV0_Mult", kMultBins, multMin, multV0);
 fHistV0Mult->SetXTitle("V0s Multiplicity");
 fHistV0Mult->SetYTitle("Counts");
   
 // MultOverOrig
 fHistMultOverOrig = new TH1F("Flow_MultOverOrig", "Flow_MultOverOrig", kMultOverOrigBins, multOverOrigMin, multOverOrigMax);
 fHistMultOverOrig->SetXTitle("Mult / Orig. Mult");
 fHistMultOverOrig->SetYTitle("Counts");
   
 // Mult correlated with the event planes
 fHistMultPart = new TH1F("Flow_MultPart", "Flow_MultPart", 2*kMultBins, multMin, multMax);
 fHistMultPart->SetXTitle("Mult of Correlated Particles");
 fHistMultPart->SetYTitle("Counts");
   
 // Mult correlated with the event planes in 1 unit rapidity
 fHistMultPartUnit = new TH1F("Flow_MultPartUnit", "Flow_MultPartUnit", 2*kMultBins, multMin, multMax/2);
 fHistMultPartUnit->SetXTitle("Mult of Correlated Particles (-0.5 < eta < 0.5)");
 fHistMultPartUnit->SetYTitle("Counts");
   
 // Mult of V0s correlated with the event planes
 fHistV0MultPart = new TH1F("FlowV0_MultPart", "FlowV0_MultPart", kMultBins, multMin, multV0);
 fHistV0MultPart->SetXTitle("Mult of Correlated V0s");
 fHistV0MultPart->SetYTitle("Counts");
   
 // VertexZ
 fHistVertexZ = new TH1F("Flow_VertexZ", "Flow_VertexZ", kVertexZBins, vertexZMin, vertexZMax);
 fHistVertexZ->SetXTitle("Vertex Z (cm)");
 fHistVertexZ->SetYTitle("Counts");
   
 // VertexXY
 fHistVertexXY2D = new TH2F("Flow_VertexXY2D", "Flow_VertexXY2D", kVertexXYBins, vertexXYMin, vertexXYMax, kVertexXYBins, vertexXYMin, vertexXYMax);
 fHistVertexXY2D->SetXTitle("Vertex X (cm)");
 fHistVertexXY2D->SetYTitle("Vertex Y (cm)");
   
 // EtaSym vs. Vertex Z Tpc
 fHistEtaSymVerZ2D = new TH2F("Flow_EtaSymVerZ2D", "Flow_EtaSymVerZ2D", kVertexZBins, vertexZMin, vertexZMax, kEtaSymBins, etaSymMin, etaSymMax);
 fHistEtaSymVerZ2D->SetXTitle("Vertex Z (cm)");
 fHistEtaSymVerZ2D->SetYTitle("Eta Symmetry TPC");

 // EtaSym Tpc
 fHistEtaSym = new TH1F("Flow_EtaSym_TPC", "Flow_EtaSym_TPC", kEtaSymBins, etaSymMin, etaSymMax);
 fHistEtaSym->SetXTitle("Eta Symmetry Ratio TPC");
 fHistEtaSym->SetYTitle("Counts");
   
 // EtaSym vs. Vertex Z Tpc - correlation analysis
 fHistEtaSymVerZ2DPart = new TH2F("Flow_EtaSymVerZ2D_part", "Flow_EtaSymVerZ2D_part", kVertexZBins, vertexZMin, vertexZMax, kEtaSymBins, etaSymMinPart, etaSymMaxPart);
 fHistEtaSymVerZ2DPart->SetXTitle("Vertex Z (cm)");
 fHistEtaSymVerZ2DPart->SetYTitle("Eta Symmetry TPC");

 // EtaSym Tpc - correlation analysis
 fHistEtaSymPart = new TH1F("Flow_EtaSym_TPC_part", "Flow_EtaSym_TPC_part", kEtaSymBins, etaSymMinPart, etaSymMaxPart);
 fHistEtaSymPart->SetXTitle("Eta Symmetry Ratio TPC");
 fHistEtaSymPart->SetYTitle("Counts");
   
 // phi , whatever
 fHistPhi = new TH1F("Flow_xPhi", "Flow_xPhi (all)", fPhiBins, fPhiMin, fPhiMax);
 fHistPhi->SetXTitle("Phi (rad)");
 fHistPhi->SetYTitle("Counts");

 // phi constrained
 fHistPhiCons = new TH1F("Flow_cPhi", "Flow_cPhi", fPhiBins, fPhiMin, fPhiMax);
 fHistPhiCons->SetXTitle("cPhi (rad)");
 fHistPhiCons->SetYTitle("Counts");
   
 // EtaPtPhi , whatever
 fHistAllEtaPtPhi3D = new TH3F("Flow_EtaPtPhi3Dall", "Flow_EtaPtPhi3Dall (whatever)", fEtaBins, fEtaMin, fEtaMax, fPtBins, fPtMin, fPtMax, kPhi3DBins, fPhiMin, fPhiMax);
 fHistAllEtaPtPhi3D->SetXTitle("Eta");
 fHistAllEtaPtPhi3D->SetYTitle("Pt (GeV/c)");
 fHistAllEtaPtPhi3D->SetZTitle("Phi (rad)");
   
 // Constrained EtaPtPhi
 fHistConsEtaPtPhi3D = new TH3F("Flow_consEtaPtPhi3D", "Flow_consEtaPtPhi3D (constrainable)", fEtaBins, fEtaMin, fEtaMax, fPtBins, fPtMin, fPtMax, kPhi3DBins, fPhiMin, fPhiMax);
 fHistConsEtaPtPhi3D->SetXTitle("cEta");
 fHistConsEtaPtPhi3D->SetYTitle("cPt (GeV/c)");
 fHistConsEtaPtPhi3D->SetZTitle("cPhi (rad)");
   
 // Global EtaPtPhi
 fHistGlobEtaPtPhi3D = new TH3F("Flow_globEtaPtPhi3D", "Flow_globEtaPtPhi3D (constrainable)", fEtaBins, fEtaMin, fEtaMax, fPtBins, fPtMin, fPtMax, kPhi3DBins, fPhiMin, fPhiMax);
 fHistGlobEtaPtPhi3D->SetXTitle("gEta");
 fHistGlobEtaPtPhi3D->SetYTitle("gPt (GeV/c)");
 fHistGlobEtaPtPhi3D->SetZTitle("gPhi (rad)");

 // UnConstrained EtaPtPhi
 fHistUncEtaPtPhi3D = new TH3F("Flow_uncEtaPtPhi3D", "Flow_uncEtaPtPhi3D (un-constrainable)", fEtaBins, fEtaMin, fEtaMax, fPtBins, fPtMin, fPtMax, kPhi3DBins, fPhiMin, fPhiMax);
 fHistUncEtaPtPhi3D->SetXTitle("gEta");
 fHistUncEtaPtPhi3D->SetYTitle("gPt (GeV/c)");
 fHistUncEtaPtPhi3D->SetZTitle("gPhi (rad)");
   
 // EtaPtPhi for particles correlated with the event plane
 fHistEtaPtPhi3DPart = new TH3F("Flow_EtaPtPhi3Dpart", "Flow_EtaPtPhi3Dpart (selected part)", fEtaBins, fEtaMin, fEtaMax, fPtBins, fPtMin, fPtMax, kPhi3DBins, fPhiMin, fPhiMax);
 fHistEtaPtPhi3DPart->SetXTitle("Eta");
 fHistEtaPtPhi3DPart->SetYTitle("Pt (GeV/c)");
 fHistEtaPtPhi3DPart->SetZTitle("Phi (rad)");
   
 // EtaPtPhi for particles NOT SELECTED for correlation with the event plane
 fHistEtaPtPhi3DOut = new TH3F("Flow_EtaPtPhi3Dout", "Flow_EtaPtPhi3Dout (NOT selected part)", fEtaBins, fEtaMin, fEtaMax, fPtBins, fPtMin, fPtMax, kPhi3DBins, fPhiMin, fPhiMax);
 fHistEtaPtPhi3DOut->SetXTitle("Eta");
 fHistEtaPtPhi3DOut->SetYTitle("Pt (GeV/c)");
 fHistEtaPtPhi3DOut->SetZTitle("Phi (rad)");
   
 // Yield Pt-Phi for all positive
 fHistPtPhiPos = new TH2D("Flow_PtPhi_Plus", "Flow_PtPhi_Plus", fPhiBins, fPhiMin, fPhiMax, 50, fPtMin, fPtMax);
 fHistPtPhiPos->Sumw2();
 fHistPtPhiPos->SetXTitle("Phi");
 fHistPtPhiPos->SetYTitle("Pt (GeV/c)");

 // Yield Pt-Phi for all negative
 fHistPtPhiNeg = new TH2D("Flow_PtPhi_Minus", "Flow_PtPhi_Minus", fPhiBins, fPhiMin, fPhiMax, 50, fPtMin, fPtMax);
 fHistPtPhiNeg->Sumw2();
 fHistPtPhiNeg->SetXTitle("Phi");
 fHistPtPhiNeg->SetYTitle("Pt (GeV/c)");
   
 // Yield for all particles
 fHistYieldAll2D = new TH2D("Flow_YieldAll2D", "Flow_YieldAll2D", fEtaBins, fEtaMin, fEtaMax, fPtBins, fPtMin, fPtMax);
 fHistYieldAll2D->Sumw2();
 fHistYieldAll2D->SetXTitle("Pseudorapidty");
 fHistYieldAll2D->SetYTitle("Pt (GeV/c)");

 // Yield for constrainable tracks
 fHistYieldCon2D = new TH2D("Flow_YieldCons2D", "Flow_YieldCons2D", fEtaBins, fEtaMin, fEtaMax, fPtBins, fPtMin, fPtMax);
 fHistYieldCon2D->Sumw2();
 fHistYieldCon2D->SetXTitle("Pseudorapidty");
 fHistYieldCon2D->SetYTitle("Pt (GeV/c)");

 // Yield for un-constrainable tracks
 fHistYieldUnc2D = new TH2D("Flow_YieldUnc2D", "Flow_YieldUnc2D", fEtaBins, fEtaMin, fEtaMax, fPtBins, fPtMin, fPtMax);
 fHistYieldUnc2D->Sumw2();
 fHistYieldUnc2D->SetXTitle("Pseudorapidty");
 fHistYieldUnc2D->SetYTitle("Pt (GeV/c)");

 // Yield for particles correlated with the event plane
 fHistYieldPart2D = new TH2D("Flow_YieldPart2D", "Flow_YieldPart2D (selected part)", fEtaBins, fEtaMin, fEtaMax, fPtBinsPart, fPtMin, fPtMaxPart);
 fHistYieldPart2D->Sumw2();
 fHistYieldPart2D->SetXTitle((char*)fLabel.Data());
 fHistYieldPart2D->SetYTitle("Pt (GeV/c)");

 // Yield for particles NOT SELECTED for correlation with the event plane
 fHistYieldOut2D = new TH2D("Flow_YieldOut2D", "Flow_YieldOut2D (NOT selected part)", fEtaBins, fEtaMin, fEtaMax, fPtBins, fPtMin, fPtMax);
 fHistYieldOut2D->Sumw2();
 fHistYieldOut2D->SetXTitle("Pseudorapidty");
 fHistYieldOut2D->SetYTitle("Pt (GeV/c)");
 
 // invariant Mass for all particles (from TOF)
 fHistInvMass = new TH1F("Flow_InvMass", "Flow_InvMass (tof)", kMassBins, massMin, massMax);
 fHistInvMass->SetXTitle("Invariant Mass (GeV)");
 fHistInvMass->SetYTitle("Counts");
 
 // invariant Mass for particles correlated with the event plane (from TOF)
 fHistInvMassPart = new TH1F("Flow_InvMassPart", "Flow_InvMassPart (tof)", kMassBins, massMin, massMax);
 fHistInvMassPart->SetXTitle("Invariant Mass (GeV)");
 fHistInvMassPart->SetYTitle("Counts");
 
 // invariant Mass for particles NOT SELECTED for correlation with the event plane (from TOF)
 fHistInvMassOut = new TH1F("Flow_InvMassOut", "Flow_InvMassOut (tof)", kMassBins, massMin, massMax);
 fHistInvMassOut->SetXTitle("Invariant Mass (GeV)");
 fHistInvMassOut->SetYTitle("Counts");

 // Mean Eta in each bin for particles correlated with the event plane
 fHistBinEta = new TProfile("Flow_Bin_Eta", "Flow_Bin_Eta_part (selected part)", fEtaBins, fEtaMin, fEtaMax, fEtaMin, fEtaMax, "");
 fHistBinEta->SetXTitle((char*)fLabel.Data());
 fHistBinEta->SetYTitle((char*)fLabel.Data());
 
 // Mean Pt in each bin for particles correlated with the event plane
 fHistBinPt = new TProfile("Flow_Bin_Pt", "Flow_Bin_Pt_part (selected part)", fPtBinsPart, fPtMin, fPtMaxPart, fPtMin, fPtMaxPart, "");
 fHistBinPt->SetXTitle("Pt (GeV/c)");
 fHistBinPt->SetYTitle("<Pt> (GeV/c)");
 
 // cos(n*phiLab)
 fHistCosPhi = new TProfile("Flow_CosPhiLab", "Flow_CosPhiLab", AliFlowConstants::kHars, 0.5, (float)(AliFlowConstants::kHars) + 0.5, -100., 100., "");
 fHistCosPhi->SetXTitle("Harmonic");
 fHistCosPhi->SetYTitle("<cos(n*PhiLab)> (%)");
   
 // PID pi+
 fHistPidPiPlus = new TH1F("Flow_PidPiPlus", "Flow_PidPiPlus", kPidBins, pidMin, pidMax);
 fHistPidPiPlus->SetXTitle("ALICE P.Id.");
 fHistPidPiPlus->SetYTitle("Counts");
   
 // PID pi-
 fHistPidPiMinus = new TH1F("Flow_PidPiMinus", "Flow_PidPiMinus", kPidBins, pidMin, pidMax);
 fHistPidPiMinus->SetXTitle("ALICE P.Id.");
 fHistPidPiMinus->SetYTitle("Counts");
   
 // PID proton
 fHistPidProton = new TH1F("Flow_PidProton", "Flow_PidProton", kPidBins, pidMin, pidMax);
 fHistPidProton->SetXTitle("ALICE P.Id.");
 fHistPidProton->SetYTitle("Counts");

 // PID anti proton
 fHistPidAntiProton = new TH1F("Flow_PidAntiProton", "Flow_PidAntiProton", kPidBins, pidMin, pidMax);
 fHistPidAntiProton->SetXTitle("ALICE P.Id.");
 fHistPidAntiProton->SetYTitle("Counts");

 // PID Kplus
 fHistPidKplus = new TH1F("Flow_PidKplus", "Flow_PidKplus", kPidBins, pidMin, pidMax);
 fHistPidKplus->SetXTitle("ALICE P.Id.");
 fHistPidKplus->SetYTitle("Counts");

 // PID Kminus
 fHistPidKminus = new TH1F("Flow_PidKminus", "Flow_PidKminus", kPidBins, pidMin, pidMax);
 fHistPidKminus->SetXTitle("ALICE P.Id.");
 fHistPidKminus->SetYTitle("Counts");

 // PID deuteron
 fHistPidDeuteron = new TH1F("Flow_PidDeuteron", "Flow_PidDeuteron", kPidBins, pidMin, pidMax);
 fHistPidDeuteron->SetXTitle("ALICE P.Id.");
 fHistPidDeuteron->SetYTitle("Counts");

 // PID anti deuteron
 fHistPidAntiDeuteron = new TH1F("Flow_PidAntiDeuteron", "Flow_PidAntiDeuteron", kPidBins, pidMin, pidMax);
 fHistPidAntiDeuteron->SetXTitle("ALICE P.Id.");
 fHistPidAntiDeuteron->SetYTitle("Counts");

 // PID electron
 fHistPidElectron = new TH1F("Flow_PidElectron", "Flow_PidElectron", kPidBins, pidMin, pidMax);
 fHistPidElectron->SetXTitle("ALICE P.Id.");
 fHistPidElectron->SetYTitle("Counts");

 // PID positron
 fHistPidPositron = new TH1F("Flow_PidPositron", "Flow_PidPositron", kPidBins, pidMin, pidMax);
 fHistPidPositron->SetXTitle("ALICE P.Id.");
 fHistPidPositron->SetYTitle("Counts");

 // PID Muon+
 fHistPidMuonPlus = new TH1F("Flow_PidMuonPlus", "Flow_PidMuonPlus", kPidBins, pidMin, pidMax);
 fHistPidMuonPlus->SetXTitle("ALICE P.Id.");
 fHistPidMuonPlus->SetYTitle("Counts");

 // PID Muon-
 fHistPidMuonMinus = new TH1F("Flow_PidMuonMinus", "Flow_PidMuonMinus", kPidBins, pidMin, pidMax);
 fHistPidMuonMinus->SetXTitle("ALICE P.Id.");
 fHistPidMuonMinus->SetYTitle("Counts");

 // PID pi+ selected
 fHistPidPiPlusPart = new TH1F("Flow_PidPiPlusPart", "Flow_PidPiPlusPart", kPidBins, pidMin, pidMax);
 fHistPidPiPlusPart->SetXTitle("ALICE P.Id. (pid = pi+)");
 fHistPidPiPlusPart->SetYTitle("Counts");
   
 // PID pi- selected
 fHistPidPiMinusPart = new TH1F("Flow_PidPiMinusPart", "Flow_PidPiMinusPart", kPidBins, pidMin, pidMax);
 fHistPidPiMinusPart->SetXTitle("ALICE P.Id. (pid = pi-)");
 fHistPidPiMinusPart->SetYTitle("Counts");
   
 // PID proton selected
 fHistPidProtonPart = new TH1F("Flow_PidProtonPart", "Flow_PidProtonPart", kPidBins, pidMin, pidMax);
 fHistPidProtonPart->SetXTitle("ALICE P.Id. (pid = p)");
 fHistPidProtonPart->SetYTitle("Counts");

 // PID anti proton selected
 fHistPidAntiProtonPart = new TH1F("Flow_PidAntiProtonPart", "Flow_PidAntiProtonPart", kPidBins, pidMin, pidMax);
 fHistPidAntiProtonPart->SetXTitle("ALICE P.Id. (pid = p-)");
 fHistPidAntiProtonPart->SetYTitle("Counts");

 // PID Kplus selected
 fHistPidKplusPart = new TH1F("Flow_PidKplusPart", "Flow_PidKplusPart", kPidBins, pidMin, pidMax);
 fHistPidKplusPart->SetXTitle("ALICE P.Id. (pid = K+)");
 fHistPidKplusPart->SetYTitle("Counts");

 // PID Kminus selected
 fHistPidKminusPart = new TH1F("Flow_PidKminusPart", "Flow_PidKminusPart", kPidBins, pidMin, pidMax);
 fHistPidKminusPart->SetXTitle("ALICE P.Id. (pid = K-)");
 fHistPidKminusPart->SetYTitle("Counts");

 // PID deuteron selected
 fHistPidDeuteronPart = new TH1F("Flow_PidDeuteronPart", "Flow_PidDeuteronPart", kPidBins, pidMin, pidMax);
 fHistPidDeuteronPart->SetXTitle("ALICE P.Id. (pid = d)");
 fHistPidDeuteronPart->SetYTitle("Counts");

 // PID anti deuteron selected
 fHistPidAntiDeuteronPart = new TH1F("Flow_PidAntiDeuteronPart", "Flow_PidAntiDeuteronPart", kPidBins, pidMin, pidMax);
 fHistPidAntiDeuteronPart->SetXTitle("ALICE P.Id. (pid = d--)");
 fHistPidAntiDeuteronPart->SetYTitle("Counts");

 // PID electron selected
 fHistPidElectronPart = new TH1F("Flow_PidElectronPart", "Flow_PidElectronPart", kPidBins, pidMin, pidMax);
 fHistPidElectronPart->SetXTitle("ALICE P.Id. (pid = e-)");
 fHistPidElectronPart->SetYTitle("Counts");

 // PID positron selected
 fHistPidPositronPart = new TH1F("Flow_PidPositronPart", "Flow_PidPositronPart", kPidBins, pidMin, pidMax);
 fHistPidPositronPart->SetXTitle("ALICE P.Id. (pid = e+)");
 fHistPidPositronPart->SetYTitle("Counts");

 // PID Muon+ selected
 fHistPidMuonPlusPart = new TH1F("Flow_PidMuonPlusPart", "Flow_PidMuonPlusPart", kPidBins, pidMin, pidMax);
 fHistPidMuonPlusPart->SetXTitle("ALICE P.Id. (pid = mu+)");
 fHistPidMuonPlusPart->SetYTitle("Counts");

 // PID Muon- selected
 fHistPidMuonMinusPart = new TH1F("Flow_PidMuonMinusPart", "Flow_PidMuonMinusPart", kPidBins, pidMin, pidMax);
 fHistPidMuonMinusPart->SetXTitle("ALICE P.Id. (pid = mu-)");
 fHistPidMuonMinusPart->SetYTitle("Counts");

 // PID multiplicities (all)
 fHistPidMult = new TProfile("Flow_PidMult", "Flow_PidMult", 15, 0.5, 15.5, "");
 fHistPidMult->SetXTitle("all, h+, h-, pi+, pi-, pr+, pr-, K+, K-, d+, d-, e-, e+, mu-, mu+");
 fHistPidMult->SetYTitle("Multiplicity");

 // PID for all tracks
 fHistBayPidMult = new TH1F("Flow_BayPidMult","Flow_BayPidMult",AliFlowConstants::kPid,-0.5,((float)AliFlowConstants::kPid-0.5));
 fHistBayPidMult->Sumw2() ;
 fHistBayPidMult->SetXTitle("e+/- , mu+/- , pi+/- , K+/- , p+/- , d+/- ");
 fHistBayPidMult->SetYTitle("Counts");
   
 // PID for particles correlated with the event plane
 fHistBayPidMultPart = new TH1F("Flow_BayPidMult_Part","Flow_BayPidMult_Part",AliFlowConstants::kPid,-0.5,((float)AliFlowConstants::kPid-0.5));
 fHistBayPidMultPart->Sumw2() ;
 fHistBayPidMultPart->SetXTitle("e+/- , mu+/- , pi+/- , K+/- , p+/- , d+/- ");
 fHistBayPidMultPart->SetYTitle("Counts");

 // Centrality
 fHistCent = new TH1F("Flow_Cent", "Flow_Cent", kCentBins, centMin, centMax);
 fHistCent->SetXTitle("Centrality Bin");
 fHistCent->SetYTitle("Counts");
   
 // Deposited Energy in ZDC
 fHistEnergyZDC = new TH2F("Flow_ZDC_E", "Flow_ZDC_E", 3, 0., 3., kZDCeBins, zdceMin, zdceMax);
 fHistEnergyZDC->SetXTitle("neutron , proton , e.m.");
 fHistEnergyZDC->SetYTitle("ZDC energy");

 // Part. versus Energy in ZDC
 fHistPartZDC = new TH1F("Flow_ZDC_Participants", "Flow_ZDC_Participants", kZDCpartBins, zdcpartMin, zdcpartMax);
 fHistPartZDC->SetXTitle("ZDC part");
 fHistPartZDC->SetYTitle("Counts");

 // MeanDedxPos TPC
 fHistMeanDedxPos2D = new TH2F("Flow_MeanDedxPos2D_TPC","Flow_MeanDedxPos2D_TPC", kMomenBins, logpMin, logpMax, kDedxBins, 0, dEdxMaxTPC);
 fHistMeanDedxPos2D->SetXTitle("log(momentum) (GeV/c)");
 fHistMeanDedxPos2D->SetYTitle("mean dEdx (+)");
 
 // MeanDedxNeg TPC
 fHistMeanDedxNeg2D = new TH2F("Flow_MeanDedxNeg2D_TPC","Flow_MeanDedxNeg2D_TPC", kMomenBins, logpMin, logpMax, kDedxBins, 0, dEdxMaxTPC);
 fHistMeanDedxNeg2D->SetXTitle("log(momentum) (GeV/c)");
 fHistMeanDedxNeg2D->SetYTitle("mean dEdx (-)");

 // MeanDedxPos ITS
 fHistMeanDedxPos2DITS = new TH2F("Flow_MeanDedxPos2D_ITS","Flow_MeanDedxPos2D_ITS", kMomenBins, logpMin, logpMax, kDedxBins, 0, dEdxMax);
 fHistMeanDedxPos2DITS->SetXTitle("log(momentum) (GeV/c)");
 fHistMeanDedxPos2DITS->SetYTitle("mean dEdx (+)");
 
 // MeanDedxNeg ITS
 fHistMeanDedxNeg2DITS = new TH2F("Flow_MeanDedxNeg2D_ITS","Flow_MeanDedxNeg2D_ITS", kMomenBins, logpMin, logpMax, kDedxBins, 0, dEdxMax);
 fHistMeanDedxNeg2DITS->SetXTitle("log(momentum) (GeV/c)");
 fHistMeanDedxNeg2DITS->SetYTitle("mean dEdx (-)");

 // MeanDedxPos TPC 3D selected Part
 fHistMeanDedxPos3DPart = new TH3F("Flow_MeanDedxPos3DPart_TPC","Flow_MeanDedxPos3DPart_TPC", kMomenBins, logpMin, logpMax, kDedxBins, 0., dEdxMaxTPC, AliFlowConstants::kPid, 0., AliFlowConstants::kPid);
 fHistMeanDedxPos3DPart->SetXTitle("log(momentum) (GeV/c)");
 fHistMeanDedxPos3DPart->SetYTitle("mean dEdx (+)");
 fHistMeanDedxPos3DPart->SetZTitle("e , mu , pi , k , p , d");
 
 // MeanDedxNeg TPC 3D selected Part
 fHistMeanDedxNeg3DPart = new TH3F("Flow_MeanDedxNeg3DPart_TPC","Flow_MeanDedxNeg3DPart_TPC", kMomenBins, logpMin, logpMax, kDedxBins, 0., dEdxMaxTPC, AliFlowConstants::kPid, 0., AliFlowConstants::kPid);
 fHistMeanDedxNeg3DPart->SetXTitle("log(momentum) (GeV/c)");
 fHistMeanDedxNeg3DPart->SetYTitle("mean dEdx (-)");
 fHistMeanDedxNeg3DPart->SetZTitle("e , mu , pi , k , p , d");

 // MeanDedxPos ITS 3D selected Part
 fHistMeanDedxPos3DPartITS = new TH3F("Flow_MeanDedxPos3DPart_ITS","Flow_MeanDedxPos3DPart_ITS", kMomenBins, logpMin, logpMax, kDedxBins, 0., dEdxMax, AliFlowConstants::kPid, 0., AliFlowConstants::kPid);
 fHistMeanDedxPos3DPartITS->SetXTitle("log(momentum) (GeV/c)");
 fHistMeanDedxPos3DPartITS->SetYTitle("mean dEdx (+)");
 fHistMeanDedxPos3DPartITS->SetZTitle("e , mu , pi , k , p , d");
 
 // MeanDedxNeg ITS 3D selected Part
 fHistMeanDedxNeg3DPartITS = new TH3F("Flow_MeanDedxNeg3DPart_ITS","Flow_MeanDedxNeg3DPart_ITS", kMomenBins, logpMin, logpMax, kDedxBins, 0., dEdxMax, AliFlowConstants::kPid, 0., AliFlowConstants::kPid);
 fHistMeanDedxNeg3DPartITS->SetXTitle("log(momentum) (GeV/c)");
 fHistMeanDedxNeg3DPartITS->SetYTitle("mean dEdx (-)");
 fHistMeanDedxNeg3DPartITS->SetZTitle("e , mu , pi , k , p , d");

 // MeanSignalPos TRD
 fHistMeanDedxPos2DTRD = new TH2F("Flow_MeanSignalPos2D_TRD","Flow_MeanSignalPos2D_TRD", kMomenBins, logpMin, logpMax, kTrdBins, 0, trdmax);
 fHistMeanDedxPos2DTRD->SetXTitle("log(momentum) (GeV/c)");
 fHistMeanDedxPos2DTRD->SetYTitle("mean TRD (+)");
 
 // MeanSignalNeg TRD
 fHistMeanDedxNeg2DTRD = new TH2F("Flow_MeanSignalNeg2D_TRD","Flow_MeanSignalNeg2D_TRD", kMomenBins, logpMin, logpMax, kTrdBins, 0, trdmax);
 fHistMeanDedxNeg2DTRD->SetXTitle("log(momentum) (GeV/c)");
 fHistMeanDedxNeg2DTRD->SetYTitle("mean TRD (-)");

 // MeanTimePos TOF
 fHistMeanDedxPos2DTOF = new TH2F("Flow_MeanTimePos2D_TOF","Flow_MeanTimePos2D_TOF", kMomenBins, logpMin, logpMax, kTofBins, tofmin, tofmax);
 fHistMeanDedxPos2DTOF->SetXTitle("log(momentum) (GeV/c)");
 fHistMeanDedxPos2DTOF->SetYTitle("mean Time of Flight (+)");
 
 // MeanTimeNeg TOF
 fHistMeanDedxNeg2DTOF = new TH2F("Flow_MeanTimeNeg2D_TOF","Flow_MeanTimeNeg2D_TOF", kMomenBins, logpMin, logpMax, kTofBins, tofmin, tofmax);
 fHistMeanDedxNeg2DTOF->SetXTitle("log(momentum) (GeV/c)");
 fHistMeanDedxNeg2DTOF->SetYTitle("mean Time of Flight (-)");
 
 // Detector response for each particle { ... 
 // MeanDedx PiPlus - TPC
 fHistMeanTPCPiPlus = new TH2F("Flow_MeanDedxPiPlusTPC", "Flow_MeanDedxPiPlusTPC", kMomenBins, logpMin, logpMax, kDedxBins, 0, dEdxMaxTPC);
 fHistMeanTPCPiPlus->SetXTitle("Log10(p) (GeV/c)");
 fHistMeanTPCPiPlus->SetYTitle("mean dEdx (pi+)");
 // MeanDedx PiMinus
 fHistMeanTPCPiMinus = new TH2F("Flow_MeanDedxPiMinusTPC", "Flow_MeanDedxPiMinusTPC", kMomenBins, logpMin, logpMax, kDedxBins, 0, dEdxMaxTPC);
 fHistMeanTPCPiMinus->SetXTitle("Log10(p) (GeV/c)");
 fHistMeanTPCPiMinus->SetYTitle("mean dEdx (pi-)");
 // MeanDedx Proton
 fHistMeanTPCProton = new TH2F("Flow_MeanDedxProtonTPC", "Flow_MeanDedxProtonTPC", kMomenBins, logpMin, logpMax, kDedxBins, 0, dEdxMaxTPC);
 fHistMeanTPCProton->SetXTitle("Log10(p) (GeV/c)");
 fHistMeanTPCProton->SetYTitle("mean dEdx (pr+)");
 // MeanDedx Pbar
 fHistMeanTPCPbar = new TH2F("Flow_MeanDedxPbarTPC", "Flow_MeanDedxPbarTPC", kMomenBins, logpMin, logpMax, kDedxBins, 0, dEdxMaxTPC);
 fHistMeanTPCPbar->SetXTitle("Log10(p) (GeV/c)");
 fHistMeanTPCPbar->SetYTitle("mean dEdx (pr-)");
 // MeanDedx Kplus
 fHistMeanTPCKplus = new TH2F("Flow_MeanDedxKplusTPC", "Flow_MeanDedxKplusTPC", kMomenBins, logpMin, logpMax, kDedxBins, 0, dEdxMaxTPC);
 fHistMeanTPCKplus->SetXTitle("Log10(p) (GeV/c)");
 fHistMeanTPCKplus->SetYTitle("mean dEdx (K+)");
 // MeanDedx Kminus
 fHistMeanTPCKminus = new TH2F("Flow_MeanDedxKminusTPC", "Flow_MeanDedxKminusTPC", kMomenBins, logpMin, logpMax, kDedxBins, 0, dEdxMaxTPC);
 fHistMeanTPCKminus->SetXTitle("Log10(p) (GeV/c)");
 fHistMeanTPCKminus->SetYTitle("mean dEdx (K-)");
 // MeanDedx Deuteron
 fHistMeanTPCDeuteron = new TH2F("Flow_MeanDedxDeuteronTPC", "Flow_MeanDedxDeuteronTPC", kMomenBins, logpMin, logpMax, kDedxBins, 0, dEdxMaxTPC);
 fHistMeanTPCDeuteron->SetXTitle("Log10(p) (GeV/c)");
 fHistMeanTPCDeuteron->SetYTitle("mean dEdx (d+)");
 // MeanDedx AntiDeuteron
 fHistMeanTPCAntiDeuteron = new TH2F("Flow_MeanDedxAntiDeuteronTPC", "Flow_MeanDedxAntiDeuteronTPC", kMomenBins, logpMin, logpMax, kDedxBins, 0, dEdxMaxTPC);
 fHistMeanTPCAntiDeuteron->SetXTitle("Log10(p) (GeV/c)");
 fHistMeanTPCAntiDeuteron->SetYTitle("mean dEdx (d-)");
 // MeanDedxElectron
 fHistMeanTPCElectron = new TH2F("Flow_MeanDedxElectronTPC", "Flow_MeanDedxElectronTPC", kMomenBins, logpMin, logpMax, kDedxBins, 0, dEdxMaxTPC);
 fHistMeanTPCElectron->SetXTitle("Log10(p) (GeV/c)");
 fHistMeanTPCElectron->SetYTitle("mean dEdx (e-)");
 // MeanDedx Positron
 fHistMeanTPCPositron = new TH2F("Flow_MeanDedxPositronTPC", "Flow_MeanDedxPositronTPC", kMomenBins, logpMin, logpMax, kDedxBins, 0, dEdxMaxTPC);
 fHistMeanTPCPositron->SetXTitle("Log10(p) (GeV/c)");
 fHistMeanTPCPositron->SetYTitle("mean dEdx (e+)");
 // MeanDedx MuonPlus
 fHistMeanTPCMuonPlus = new TH2F("Flow_MeanDedxMuonPlusTPC", "Flow_MeanDedxMuonPlusTPC", kMomenBins, logpMin, logpMax, kDedxBins, 0, dEdxMaxTPC);
 fHistMeanTPCMuonPlus->SetXTitle("Log10(p) (GeV/c)");
 fHistMeanTPCMuonPlus->SetYTitle("mean dEdx (mu+)");
 // MeanDedx MuonMinus
 fHistMeanTPCMuonMinus = new TH2F("Flow_MeanDedxMuonMinusTPC", "Flow_MeanDedxMuonMinusTPC", kMomenBins, logpMin, logpMax, kDedxBins, 0, dEdxMaxTPC);
 fHistMeanTPCMuonMinus->SetXTitle("Log10(p) (GeV/c)");
 fHistMeanTPCMuonMinus->SetYTitle("mean dEdx (mu-)");
 
 // MeanDedx PiPlus - ITS
 fHistMeanITSPiPlus = new TH2F("Flow_MeanDedxPiPlusITS", "Flow_MeanDedxPiPlusITS", kMomenBins, logpMin, logpMax, kDedxBins, 0, dEdxMax);
 fHistMeanITSPiPlus->SetXTitle("Log10(p) (GeV/c)");
 fHistMeanITSPiPlus->SetYTitle("mean dEdx (pi+)");
 // MeanDedx PiMinus
 fHistMeanITSPiMinus = new TH2F("Flow_MeanDedxPiMinusITS", "Flow_MeanDedxPiMinusITS", kMomenBins, logpMin, logpMax, kDedxBins, 0, dEdxMax);
 fHistMeanITSPiMinus->SetXTitle("Log10(p) (GeV/c)");
 fHistMeanITSPiMinus->SetYTitle("mean dEdx (pi-)");
 // MeanDedx Proton
 fHistMeanITSProton = new TH2F("Flow_MeanDedxProtonITS", "Flow_MeanDedxProtonITS", kMomenBins, logpMin, logpMax, kDedxBins, 0, dEdxMax);
 fHistMeanITSProton->SetXTitle("Log10(p) (GeV/c)");
 fHistMeanITSProton->SetYTitle("mean dEdx (pr+)");
 // MeanDedx Pbar
 fHistMeanITSPbar = new TH2F("Flow_MeanDedxPbarITS", "Flow_MeanDedxPbarITS", kMomenBins, logpMin, logpMax, kDedxBins, 0, dEdxMax);
 fHistMeanITSPbar->SetXTitle("Log10(p) (GeV/c)");
 fHistMeanITSPbar->SetYTitle("mean dEdx (pr-)");
 // MeanDedx Kplus
 fHistMeanITSKplus = new TH2F("Flow_MeanDedxKplusITS", "Flow_MeanDedxKplusITS", kMomenBins, logpMin, logpMax, kDedxBins, 0, dEdxMax);
 fHistMeanITSKplus->SetXTitle("Log10(p) (GeV/c)");
 fHistMeanITSKplus->SetYTitle("mean dEdx (K+)");
 // MeanDedx Kminus
 fHistMeanITSKminus = new TH2F("Flow_MeanDedxKminusITS", "Flow_MeanDedxKminusITS", kMomenBins, logpMin, logpMax, kDedxBins, 0, dEdxMax);
 fHistMeanITSKminus->SetXTitle("Log10(p) (GeV/c)");
 fHistMeanITSKminus->SetYTitle("mean dEdx (K-)");
 // MeanDedx Deuteron
 fHistMeanITSDeuteron = new TH2F("Flow_MeanDedxDeuteronITS", "Flow_MeanDedxDeuteronITS", kMomenBins, logpMin, logpMax, kDedxBins, 0, dEdxMax);
 fHistMeanITSDeuteron->SetXTitle("Log10(p) (GeV/c)");
 fHistMeanITSDeuteron->SetYTitle("mean dEdx (d+)");
 // MeanDedx AntiDeuteron
 fHistMeanITSAntiDeuteron = new TH2F("Flow_MeanDedxAntiDeuteronITS", "Flow_MeanDedxAntiDeuteronITS", kMomenBins, logpMin, logpMax, kDedxBins, 0, dEdxMax);
 fHistMeanITSAntiDeuteron->SetXTitle("Log10(p) (GeV/c)");
 fHistMeanITSAntiDeuteron->SetYTitle("mean dEdx (d-)");
 // MeanDedx Electron
 fHistMeanITSElectron = new TH2F("Flow_MeanDedxElectronITS", "Flow_MeanDedxElectronITS", kMomenBins, logpMin, logpMax, kDedxBins, 0, dEdxMax);
 fHistMeanITSElectron->SetXTitle("Log10(p) (GeV/c)");
 fHistMeanITSElectron->SetYTitle("mean dEdx (e-)");
 // MeanDedx Positron
 fHistMeanITSPositron = new TH2F("Flow_MeanDedxPositronITS", "Flow_MeanDedxPositronITS", kMomenBins, logpMin, logpMax, kDedxBins, 0, dEdxMax);
 fHistMeanITSPositron->SetXTitle("Log10(p) (GeV/c)");
 fHistMeanITSPositron->SetYTitle("mean dEdx (e+)");
 // MeanDedx MuonPlus
 fHistMeanITSMuonPlus = new TH2F("Flow_MeanDedxMuonPlusITS", "Flow_MeanDedxMuonPlusITS", kMomenBins, logpMin, logpMax, kDedxBins, 0, dEdxMax);
 fHistMeanITSMuonPlus->SetXTitle("Log10(p) (GeV/c)");
 fHistMeanITSMuonPlus->SetYTitle("mean dEdx (mu+)");
 // MeanDedx MuonMinus
 fHistMeanITSMuonMinus = new TH2F("Flow_MeanDedxMuonMinusITS", "Flow_MeanDedxMuonMinusITS", kMomenBins, logpMin, logpMax, kDedxBins, 0, dEdxMax);
 fHistMeanITSMuonMinus->SetXTitle("Log10(p) (GeV/c)");
 fHistMeanITSMuonMinus->SetYTitle("mean dEdx (mu-)");
 
 // MeanTrd PiPlus - TRD
 fHistMeanTRDPiPlus = new TH2F("Flow_MeanTrdPiPlusTRD", "Flow_MeanTrdPiPlusTRD", kMomenBins, logpMin, logpMax, kTrdBins, 0, trdmax);
 fHistMeanTRDPiPlus->SetXTitle("Log10(p) (GeV/c)");
 fHistMeanTRDPiPlus->SetYTitle("mean t.r.[] (pi+)");
 // MeanTrd PiMinus
 fHistMeanTRDPiMinus = new TH2F("Flow_MeanTrdPiMinusTRD", "Flow_MeanTrdPiMinusTRD", kMomenBins, logpMin, logpMax, kTrdBins, 0, trdmax);
 fHistMeanTRDPiMinus->SetXTitle("Log10(p) (GeV/c)");
 fHistMeanTRDPiMinus->SetYTitle("mean t.r.[] (pi-)");
 // MeanTrd Proton
 fHistMeanTRDProton = new TH2F("Flow_MeanTrdProtonTRD", "Flow_MeanTrdProtonTRD", kMomenBins, logpMin, logpMax, kTrdBins, 0, trdmax);
 fHistMeanTRDProton->SetXTitle("Log10(p) (GeV/c)");
 fHistMeanTRDProton->SetYTitle("mean t.r.[] (pr+)");
 // MeanTrd Pbar
 fHistMeanTRDPbar = new TH2F("Flow_MeanTrdPbarTRD", "Flow_MeanTrdPbarTRD", kMomenBins, logpMin, logpMax, kTrdBins, 0, trdmax);
 fHistMeanTRDPbar->SetXTitle("Log10(p) (GeV/c)");
 fHistMeanTRDPbar->SetYTitle("mean t.r.[] (pr-)");
 // MeanTrd Kplus
 fHistMeanTRDKplus = new TH2F("Flow_MeanTrdKplusTRD", "Flow_MeanTrdKplusTRD", kMomenBins, logpMin, logpMax, kTrdBins, 0, trdmax);
 fHistMeanTRDKplus->SetXTitle("Log10(p) (GeV/c)");
 fHistMeanTRDKplus->SetYTitle("mean t.r.[] (K+)");
 // MeanTrd Kminus
 fHistMeanTRDKminus = new TH2F("Flow_MeanTrdKminusTRD", "Flow_MeanTrdKminusTRD", kMomenBins, logpMin, logpMax, kTrdBins, 0, trdmax);
 fHistMeanTRDKminus->SetXTitle("Log10(p) (GeV/c)");
 fHistMeanTRDKminus->SetYTitle("mean t.r.[] (K-)");
 // MeanTrd Deuteron
 fHistMeanTRDDeuteron = new TH2F("Flow_MeanTrdDeuteronTRD", "Flow_MeanTrdDeuteronTRD", kMomenBins, logpMin, logpMax, kTrdBins, 0, trdmax);
 fHistMeanTRDDeuteron->SetXTitle("Log10(p) (GeV/c)");
 fHistMeanTRDDeuteron->SetYTitle("mean t.r.[] (d+)");
 // MeanTrd AntiDeuteron
 fHistMeanTRDAntiDeuteron = new TH2F("Flow_MeanTrdAntiDeuteronTRD", "Flow_MeanTrdAntiDeuteronTRD", kMomenBins, logpMin, logpMax, kTrdBins, 0, trdmax);
 fHistMeanTRDAntiDeuteron->SetXTitle("Log10(p) (GeV/c)");
 fHistMeanTRDAntiDeuteron->SetYTitle("mean t.r.[] (d-)");
 // MeanTrd Electron
 fHistMeanTRDElectron = new TH2F("Flow_MeanTrdElectronTRD", "Flow_MeanTrdElectronTRD", kMomenBins, logpMin, logpMax, kTrdBins, 0, trdmax);
 fHistMeanTRDElectron->SetXTitle("Log10(p) (GeV/c)");
 fHistMeanTRDElectron->SetYTitle("mean t.r.[] (e-)");
 // MeanTrd Positron
 fHistMeanTRDPositron = new TH2F("Flow_MeanTrdPositronTRD", "Flow_MeanTrdPositronTRD", kMomenBins, logpMin, logpMax, kTrdBins, 0, trdmax);
 fHistMeanTRDPositron->SetXTitle("Log10(p) (GeV/c)");
 fHistMeanTRDPositron->SetYTitle("mean t.r.[] (e+)");
 // MeanTrd MuonPlus
 fHistMeanTRDMuonPlus = new TH2F("Flow_MeanTrdMuonPlusTRD", "Flow_MeanTrdMuonPlusTRD", kMomenBins, logpMin, logpMax, kTrdBins, 0, trdmax);
 fHistMeanTRDMuonPlus->SetXTitle("Log10(p) (GeV/c)");
 fHistMeanTRDMuonPlus->SetYTitle("mean t.r.[] (mu+)");
 // MeanTrd MuonMinus
 fHistMeanTRDMuonMinus = new TH2F("Flow_MeanTrdMuonMinusTRD", "Flow_MeanTrdMuonMinusTRD", kMomenBins, logpMin, logpMax, kTrdBins, 0, trdmax);
 fHistMeanTRDMuonMinus->SetXTitle("Log10(p) (GeV/c)");
 fHistMeanTRDMuonMinus->SetYTitle("mean t.r.[] (mu-)");
 
 // T.O.F. PiPlus - TOF
 fHistMeanTOFPiPlus = new TH2F("Flow_MeanTofPiPlusTOF", "Flow_MeanTofPiPlusTOF", kMassBins, massMin, massMax, kTofBins, tofmin, tofmax);
 fHistMeanTOFPiPlus->SetXTitle("invariant mass (GeV)");
 fHistMeanTOFPiPlus->SetYTitle("mean t.o.f.[psec] (pi+)");
 // MeanTof PiMinus
 fHistMeanTOFPiMinus = new TH2F("Flow_MeanTofPiMinusTOF", "Flow_MeanTofPiMinusTOF", kMassBins, massMin, massMax, kTofBins, tofmin, tofmax);
 fHistMeanTOFPiMinus->SetXTitle("invariant mass (GeV)");
 fHistMeanTOFPiMinus->SetYTitle("mean t.o.f.[psec] (pi-)");
 // MeanTof Proton
 fHistMeanTOFProton = new TH2F("Flow_MeanTofProtonTOF", "Flow_MeanTofProtonTOF", kMassBins, massMin, massMax, kTofBins, tofmin, tofmax);
 fHistMeanTOFProton->SetXTitle("invariant mass (GeV)");
 fHistMeanTOFProton->SetYTitle("mean t.o.f.[psec] (pr+)");
 // Mean TofPbar
 fHistMeanTOFPbar = new TH2F("Flow_MeanTofPbarTOF", "Flow_MeanTofPbarTOF", kMassBins, massMin, massMax, kTofBins, tofmin, tofmax);
 fHistMeanTOFPbar->SetXTitle("invariant mass (GeV)");
 fHistMeanTOFPbar->SetYTitle("mean t.o.f.[psec] (pr-)");
 // mean t.o.f.[psec]Kplus
 fHistMeanTOFKplus = new TH2F("Flow_MeanTofKplusTOF", "Flow_MeanTofKplusTOF", kMassBins, massMin, massMax, kTofBins, tofmin, tofmax);
 fHistMeanTOFKplus->SetXTitle("invariant mass (GeV)");
 fHistMeanTOFKplus->SetYTitle("mean t.o.f.[psec] (K+)");
 // mean t.o.f.[psec]Kminus
 fHistMeanTOFKminus = new TH2F("Flow_MeanTofKminusTOF", "Flow_MeanTofKminusTOF", kMassBins, massMin, massMax, kTofBins, tofmin, tofmax);
 fHistMeanTOFKminus->SetXTitle("invariant mass (GeV)");
 fHistMeanTOFKminus->SetYTitle("mean t.o.f.[psec] (K-)");
 // MeanTof Deuteron
 fHistMeanTOFDeuteron = new TH2F("Flow_MeanTofDeuteronTOF", "Flow_MeanTofDeuteronTOF", kMassBins, massMin, massMax, kTofBins, tofmin, tofmax);
 fHistMeanTOFDeuteron->SetXTitle("invariant mass (GeV)");
 fHistMeanTOFDeuteron->SetYTitle("mean t.o.f.[psec] (d+)");
 // MeanTof AntiDeuteron
 fHistMeanTOFAntiDeuteron = new TH2F("Flow_MeanTofAntiDeuteronTOF", "Flow_MeanTofAntiDeuteronTOF", kMassBins, massMin, massMax, kTofBins, tofmin, tofmax);
 fHistMeanTOFAntiDeuteron->SetXTitle("invariant mass (GeV)");
 fHistMeanTOFAntiDeuteron->SetYTitle("mean t.o.f.[psec] (d-)");
 // MeanTof Electron
 fHistMeanTOFElectron = new TH2F("Flow_MeanTofElectronTOF", "Flow_MeanTofElectronTOF", kMassBins, massMin, massMax, kTofBins, tofmin, tofmax);
 fHistMeanTOFElectron->SetXTitle("invariant mass (GeV)");
 fHistMeanTOFElectron->SetYTitle("mean t.o.f.[psec] (e-)");
 // MeanTof Positron
 fHistMeanTOFPositron = new TH2F("Flow_MeanTofPositronTOF", "Flow_MeanTofPositronTOF", kMassBins, massMin, massMax, kTofBins, tofmin, tofmax);
 fHistMeanTOFPositron->SetXTitle("invariant mass (GeV)");
 fHistMeanTOFPositron->SetYTitle("mean t.o.f.[psec] (e+)");
 // MeanTof MuonPlus
 fHistMeanTOFMuonPlus = new TH2F("Flow_MeanTofMuonPlusTOF", "Flow_MeanTofMuonPlusTOF", kMassBins, massMin, massMax, kTofBins, tofmin, tofmax);
 fHistMeanTOFMuonPlus->SetXTitle("invariant mass (GeV)");
 fHistMeanTOFMuonPlus->SetYTitle("mean t.o.f.[psec] (mu+)");
 // MeanTof MuonMinus
 fHistMeanTOFMuonMinus = new TH2F("Flow_MeanTofMuonMinusTOF", "Flow_MeanTofMuonMinusTOF", kMassBins, massMin, massMax, kTofBins, tofmin, tofmax);
 fHistMeanTOFMuonMinus->SetXTitle("invariant mass (GeV)");
 fHistMeanTOFMuonMinus->SetYTitle("mean t.o.f.[psec] (mu-)");
 // ... }

 TString* histTitle;
 for (int n = 0; n < AliFlowConstants::kSubs; n++) // for sub-events
 {
  for (int k = 0; k < AliFlowConstants::kSels; k++) 
  {
   for (int j = 0; j < AliFlowConstants::kHars; j++) 
   {
    float order = (float)(j + 1);
    int i = AliFlowConstants::kSubs * k + n ;

    // event planes
    histTitle = new TString("Flow_Psi_Sub");
    *histTitle += n+1;
    histTitle->Append("_Sel");
    *histTitle += k+1;
    histTitle->Append("_Har");
    *histTitle += j+1;
    fHistSub[i].fHistSubHar[j].fHistPsiSubs = new TH1F(histTitle->Data(),histTitle->Data(), kPsiBins, psiMin, (psiMax / order));
    fHistSub[i].fHistSubHar[j].fHistPsiSubs->SetXTitle("Event Plane Angle (rad)");
    fHistSub[i].fHistSubHar[j].fHistPsiSubs->SetYTitle("Counts");
    delete histTitle;
   }
  }
 }
 
 if(fV0loop)            // All V0s (if there, if flag on)
 {
 // Mass
  fHistV0Mass = new TH1F("FlowV0_InvMass", "FlowV0_InvMass", kMassBins, massMin, massMax);
  fHistV0Mass->SetXTitle("Invariant Mass (GeV)");
  fHistV0Mass->SetYTitle("Counts");
 // Distance of closest approach
  fHistV0Dca = new TH1F("FlowV0_Dca", "FlowV0_Dca", kDcaBins, dcaMin, glDcaMax);
  fHistV0Dca->SetXTitle("dca between tracks (cm)");
  fHistV0Dca->SetYTitle("Counts");
 // lenght
  fHistV0Lenght = new TH1F("FlowV0_Lenght", "FlowV0_Lenght", kLgBins, lgMinV0, lgMaxV0);
  fHistV0Lenght->SetXTitle("Distance of V0s (cm)");
  fHistV0Lenght->SetYTitle("Counts");
 // Sigma for all particles
  fHistV0Sigma = new TH1F("FlowV0_Sigma", "FlowV0_Sigma", kLgBins, lgMinV0, lgMaxV0 );
  fHistV0Sigma->SetXTitle("Sigma");
  fHistV0Sigma->SetYTitle("Counts");
 // Chi2 
  fHistV0Chi2 = new TH1F("FlowV0_Chi2", "FlowV0_Chi2", kChi2Bins, chi2Min, chi2MaxC);
  fHistV0Chi2->SetXTitle("Chi square at Main Vertex");
  fHistV0Chi2->SetYTitle("Counts");
 // EtaPtPhi
  fHistV0EtaPtPhi3D = new TH3F("FlowV0_EtaPtPhi3D", "FlowV0_EtaPtPhi3D", fEtaBins, fEtaMin, fEtaMax, fPtBins, fPtMin, fPtMax, kPhi3DBins, fPhiMin, fPhiMax);
  fHistV0EtaPtPhi3D->SetXTitle("Eta");
  fHistV0EtaPtPhi3D->SetYTitle("Pt (GeV/c)");
  fHistV0EtaPtPhi3D->SetZTitle("Phi (rad)");
 // Yield for all v0s
  fHistV0YieldAll2D = new TH2D("FlowV0_PtEta2Dall", "FlowV0_PtEta2Dall", fEtaBins, fEtaMin, fEtaMax, fPtBins, fPtMin, fPtMax);
  fHistV0YieldAll2D->Sumw2();
  fHistV0YieldAll2D->SetXTitle("Pseudorapidty");
  fHistV0YieldAll2D->SetYTitle("Pt (GeV/c)");
 // Mass slices on pT
  fHistV0MassPtSlices = new TH2D("FlowV0_MassPtSlices", "FlowV0_MassPtSlices", kMassBins, massMin, massMax, fPtBins, fPtMin, fPtMax);
  fHistV0MassPtSlices->Sumw2();
  fHistV0MassPtSlices->SetXTitle("Invariant Mass (GeV)");
  fHistV0MassPtSlices->SetYTitle("Pt (GeV/c)");
 // Yield v0s (total P and Rapidity)
  fHistV0PYall2D = new TH2D("FlowV0_PYall2D", "FlowV0_PYall2D", fEtaBins, fEtaMin, fEtaMax, fPtBins, fPtMin, fPtMax);
  fHistV0PYall2D->Sumw2();
  fHistV0PYall2D->SetXTitle("Rapidty");
  fHistV0PYall2D->SetYTitle("P (GeV/c)");

 // Selected V0s ...  
 // Yield 
  fHistV0YieldPart2D = new TH2D("FlowV0_PtEta2Dpart", "FlowV0_PtEta2Dpart", fEtaBins, fEtaMin, fEtaMax, fPtBins, fPtMin, fPtMax);
  fHistV0YieldPart2D->Sumw2();
  fHistV0YieldPart2D->SetXTitle("Pseudorapidty");
  fHistV0YieldPart2D->SetYTitle("Pt (GeV/c)");
 // Mass Window
  fHistV0MassWin = new TH1F("FlowV0_MassWinPart", "FlowV0_MassWinPart", kMassBins, massMin, massMax);
  fHistV0MassWin->SetXTitle("Invariant Mass (GeV)");
  fHistV0MassWin->SetYTitle("Counts");
 // EtaPtPhi
  fHistV0EtaPtPhi3DPart = new TH3F("FlowV0_EtaPtPhi3Dpart", "FlowV0_EtaPtPhi3Dpart", fEtaBins, fEtaMin, fEtaMax, fPtBins, fPtMin, fPtMax, kPhi3DBins, fPhiMin, fPhiMax);
  fHistV0EtaPtPhi3DPart->SetXTitle("Eta");
  fHistV0EtaPtPhi3DPart->SetYTitle("Pt (GeV/c)");
  fHistV0EtaPtPhi3DPart->SetZTitle("Phi (rad)");
 // Distance of closest approach
  fHistV0DcaPart = new TH1F("FlowV0_DcaPart", "FlowV0_DcaPart", kDcaBins, dcaMin, dcaMax);
  fHistV0DcaPart->Sumw2();
  fHistV0DcaPart->SetXTitle("dca between tracks (cm)");
  fHistV0DcaPart->SetYTitle("Counts");
 // lenght
  fHistV0LenghtPart = new TH1F("FlowV0_LenghtPart", "FlowV0_LenghtPart", kLgBins, lgMinV0, lgMaxV0);
  fHistV0LenghtPart->SetXTitle("Distance of V0s (cm)");
  fHistV0LenghtPart->SetYTitle("Counts");
 // SideBand Mass (sidebands)
  fHistV0sbMassSide = new TH1F("FlowV0sb_MassWinSideBands", "FlowV0sb_MassWinSideBands", kMassBins, massMin, massMax);
  fHistV0sbMassSide->SetXTitle("Invariant Mass (GeV)");
  fHistV0sbMassSide->SetYTitle("Counts");
 // EtaPtPhi (sidebands)
  fHistV0sbEtaPtPhi3DPart = new TH3F("FlowV0sb_EtaPtPhi3D", "FlowV0sb_EtaPtPhi3D", fEtaBins, fEtaMin, fEtaMax, fPtBins, fPtMin, fPtMax, kPhi3DBins, fPhiMin, fPhiMax);
  fHistV0sbEtaPtPhi3DPart->SetXTitle("Eta");
  fHistV0sbEtaPtPhi3DPart->SetYTitle("Pt (GeV/c)");
  fHistV0sbEtaPtPhi3DPart->SetZTitle("Phi (rad)");
 // Distance of closest approach (sidebands)
  fHistV0sbDcaPart = new TH1F("FlowV0sb_Dca", "FlowV0sb_Dca", kDcaBins, dcaMin, dcaMax);
  fHistV0sbDcaPart->Sumw2();
  fHistV0sbDcaPart->SetXTitle("dca between tracks (cm)");
  fHistV0sbDcaPart->SetYTitle("Counts");
 // lenght (sidebands)
  fHistV0sbLenghtPart = new TH1F("FlowV0sb_Lenght", "FlowV0sb_Lenght", kLgBins, lgMinV0, lgMaxV0);
  fHistV0sbLenghtPart->SetXTitle("Distance of V0s (cm)");
  fHistV0sbLenghtPart->SetYTitle("Counts");

 // Mean Eta in each bin
  fHistV0BinEta = new TProfile("FlowV0_Bin_Eta", "FlowV0_Bin_Eta", fEtaBins, fEtaMin, fEtaMax, fEtaMin, fEtaMax, "");
  fHistV0BinEta->SetXTitle((char*)fLabel.Data());
  fHistV0BinEta->SetYTitle("<Eta>");
 // Mean Pt in each bin
  fHistV0BinPt = new TProfile("FlowV0_Bin_Pt", "FlowV0_Bin_Pt", fPtBinsPart, fPtMin, fPtMaxPart, fPtMin, fPtMaxPart, "");
  fHistV0BinPt->SetXTitle("Pt (GeV/c)");
  fHistV0BinPt->SetYTitle("<Pt> (GeV/c)");
 // Mean Eta in each bin (sidebands)
  fHistV0sbBinEta = new TProfile("FlowV0sb_Bin_Eta", "FlowV0sb_Bin_Eta", fEtaBins, fEtaMin, fEtaMax, fEtaMin, fEtaMax, "");
  fHistV0sbBinEta->SetXTitle((char*)fLabel.Data());
  fHistV0sbBinEta->SetYTitle("<Eta>");
 // Mean Pt in each bin (sidebands)
  fHistV0sbBinPt = new TProfile("FlowV0sb_Bin_Pt", "FlowV0sb_Bin_Pt", fPtBinsPart, fPtMin, fPtMaxPart, fPtMin, fPtMaxPart, "");
  fHistV0sbBinPt->SetXTitle("Pt (GeV/c)");
  fHistV0sbBinPt->SetYTitle("<Pt> (GeV/c)");
 }

 for (int k = 0; k < AliFlowConstants::kSels; k++) // for each selection
 {
  // cos(n*delta_Psi)
  histTitle = new TString("Flow_Cos_Sel");
  *histTitle += k+1;
  fHistFull[k].fHistCos = new TProfile(histTitle->Data(), histTitle->Data(), AliFlowConstants::kHars, 0.5, (float)(AliFlowConstants::kHars) + 0.5, -1., 1., "");
  fHistFull[k].fHistCos->SetXTitle("Harmonic");
  fHistFull[k].fHistCos->SetYTitle("<cos(n*delta_Psi)>");
  delete histTitle;
   
  // resolution
  histTitle = new TString("Flow_Res_Sel");
  *histTitle += k+1;
  fHistFull[k].fHistRes = new TH1F(histTitle->Data(), histTitle->Data(), AliFlowConstants::kHars, 0.5, (float)(AliFlowConstants::kHars) + 0.5);
  fHistFull[k].fHistRes->SetXTitle("Harmonic");
  fHistFull[k].fHistRes->SetYTitle("Resolution");
  delete histTitle;

  // vObs
  histTitle = new TString("Flow_vObs_Sel");
  *histTitle += k+1;
  fHistFull[k].fHistvObs = new TProfile(histTitle->Data(), histTitle->Data(), AliFlowConstants::kHars, 0.5, (float)(AliFlowConstants::kHars) + 0.5, -100., 100., "");
  fHistFull[k].fHistvObs->SetXTitle("Harmonic");
  fHistFull[k].fHistvObs->SetYTitle("vObs (%)");
  delete histTitle;

  // vObs V0
  histTitle = new TString("FlowV0_vObs_Sel");
  *histTitle += k+1;
  fHistFull[k].fHistV0vObs = new TProfile(histTitle->Data(), histTitle->Data(), AliFlowConstants::kHars, 0.5, (float)(AliFlowConstants::kHars) + 0.5, -100., 100., "");
  fHistFull[k].fHistV0vObs->SetXTitle("Harmonic");
  fHistFull[k].fHistV0vObs->SetYTitle("vObs (%)");
  delete histTitle;
   
  // vObs V0 sideband SX
  histTitle = new TString("FlowV0sb_vObs_sx_Sel");
  *histTitle += k+1;
  fHistFull[k].fHistV0sbvObsSx = new TProfile(histTitle->Data(), histTitle->Data(), AliFlowConstants::kHars, 0.5, (float)(AliFlowConstants::kHars) + 0.5, -100., 100., "");
  fHistFull[k].fHistV0sbvObsSx->SetXTitle("Harmonic");
  fHistFull[k].fHistV0sbvObsSx->SetYTitle("vObs (%)");
  delete histTitle;
   
  // vObs V0 sideband DX
  histTitle = new TString("FlowV0sb_vObs_dx_Sel");
  *histTitle += k+1;
  fHistFull[k].fHistV0sbvObsDx = new TProfile(histTitle->Data(), histTitle->Data(), AliFlowConstants::kHars, 0.5, (float)(AliFlowConstants::kHars) + 0.5, -100., 100., "");
  fHistFull[k].fHistV0sbvObsDx->SetXTitle("Harmonic");
  fHistFull[k].fHistV0sbvObsDx->SetYTitle("vObs (%)");
  delete histTitle;
   
  // PID for tracks used in R.P.
  histTitle = new TString("Flow_BayPidMult_Sel");
  *histTitle += k+1;
  fHistFull[k].fHistBayPidMult = new TH1F(histTitle->Data(), histTitle->Data(),AliFlowConstants::kPid,-0.5,((float)AliFlowConstants::kPid-0.5));
  fHistFull[k].fHistBayPidMult->Sumw2() ;
  fHistFull[k].fHistBayPidMult->SetXTitle("e+/-  ,  mu+/-  ,  pi+/-  ,  K+/-  ,  p+/-  ,  d+/- ");
  fHistFull[k].fHistBayPidMult->SetYTitle("Counts");
  delete histTitle;

  for (int j = 0; j < AliFlowConstants::kHars; j++)   // for each harmonic
  {
   float order  = (float)(j+1);

   // multiplicity
   histTitle = new TString("Flow_Mul_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   fHistFull[k].fHistFullHar[j].fHistMult = new TH1F(histTitle->Data(),histTitle->Data(), kMultBins, multMin, multMax);
   fHistFull[k].fHistFullHar[j].fHistMult->SetXTitle("Multiplicity");
   fHistFull[k].fHistFullHar[j].fHistMult->SetYTitle("Counts");
   delete histTitle;

   // event plane
   histTitle = new TString("Flow_Psi_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   fHistFull[k].fHistFullHar[j].fHistPsi = new TH1F(histTitle->Data(), histTitle->Data(), kPsiBins, psiMin, psiMax / order);
   fHistFull[k].fHistFullHar[j].fHistPsi->SetXTitle("Event Plane Angle (rad)");
   fHistFull[k].fHistFullHar[j].fHistPsi->SetYTitle("Counts");
   delete histTitle;
     
   // event plane difference of two selections
   histTitle = new TString("Flow_Psi_Diff_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   if (k == 0 ) 
   {
    Int_t myOrder = 1;
    if (j == 1) { myOrder = 2 ; }
    fHistFull[k].fHistFullHar[j].fHistPsiDiff = new TH1F(histTitle->Data(), histTitle->Data(), kPsiBins, -psiMax/myOrder/2., psiMax/myOrder/2.);
   } 
   else 
   {
    fHistFull[k].fHistFullHar[j].fHistPsiDiff = new TH1F(histTitle->Data(), histTitle->Data(), kPsiBins, -psiMax/2., psiMax/2.);
   }
   if (k == 0) 
   {
    if (j == 0) 
    {
     fHistFull[k].fHistFullHar[j].fHistPsiDiff->SetXTitle("#Psi_{1,Sel1} - #Psi_{1,Sel2}(rad)");
    } 
    else if (j == 1) 
    {
     fHistFull[k].fHistFullHar[j].fHistPsiDiff->SetXTitle("#Psi_{2,Sel1} - #Psi_{2,Sel2}(rad)");
    }
   } 
   else if (k == 1) 
   {
    if (j == 0) 
    {  
     fHistFull[k].fHistFullHar[j].fHistPsiDiff->SetXTitle("#Psi_{1,Sel1} - #Psi_{2,Sel2}(rad)");
    } 
    else if (j == 1) 
    {
     fHistFull[k].fHistFullHar[j].fHistPsiDiff->SetXTitle("#Psi_{1,Sel1} - #Psi_{2,Sel1}(rad)");
    }
   }
   fHistFull[k].fHistFullHar[j].fHistPsiDiff->SetYTitle("Counts");
   delete histTitle;

   // correlation of sub-event planes
   histTitle = new TString("Flow_Psi_Sub_Corr_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   fHistFull[k].fHistFullHar[j].fHistPsiSubCorr = new TH1F(histTitle->Data(), histTitle->Data(), kPsiBins, psiMin, psiMax / order);
   fHistFull[k].fHistFullHar[j].fHistPsiSubCorr->Sumw2();
   fHistFull[k].fHistFullHar[j].fHistPsiSubCorr->SetXTitle("Sub-Event Correlation (rad)");
   fHistFull[k].fHistFullHar[j].fHistPsiSubCorr->SetYTitle("Counts");
   delete histTitle;
     
   // correlation of sub-event planes of different order
   histTitle = new TString("Flow_Psi_Sub_Corr_Diff_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   fHistFull[k].fHistFullHar[j].fHistPsiSubCorrDiff = new TH1F(histTitle->Data(), histTitle->Data(), kPsiBins, psiMin, psiMax / (order+1.));
   fHistFull[k].fHistFullHar[j].fHistPsiSubCorrDiff->Sumw2();
   fHistFull[k].fHistFullHar[j].fHistPsiSubCorrDiff->SetXTitle("Sub-Event Correlation (rad)");
   fHistFull[k].fHistFullHar[j].fHistPsiSubCorrDiff->SetYTitle("Counts");
   delete histTitle;
     
   // q
   histTitle = new TString("Flow_NormQ_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   fHistFull[k].fHistFullHar[j].fHistQnorm = new TH1F(histTitle->Data(), histTitle->Data(), kQbins, qMin, qMax);
   fHistFull[k].fHistFullHar[j].fHistQnorm->Sumw2();
   fHistFull[k].fHistFullHar[j].fHistQnorm->SetXTitle("q = |Q|/sqrt(Mult)");
   fHistFull[k].fHistFullHar[j].fHistQnorm->SetYTitle("Counts");
   delete histTitle;

    // particle-plane azimuthal correlation
   histTitle = new TString("Flow_Phi_Corr_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   fHistFull[k].fHistFullHar[j].fHistPhiCorr = new TH1F(histTitle->Data(), histTitle->Data(), fPhiBins, fPhiMin, fPhiMax / order);
   fHistFull[k].fHistFullHar[j].fHistPhiCorr->Sumw2();
   fHistFull[k].fHistFullHar[j].fHistPhiCorr->SetXTitle("Particle-Plane Correlation (rad)");
   fHistFull[k].fHistFullHar[j].fHistPhiCorr->SetYTitle("Counts");
   delete histTitle;
     
   // neutral particle-plane azimuthal correlation
   histTitle = new TString("FlowV0_Phi_Corr_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   fHistFull[k].fHistFullHar[j].fHistV0PhiCorr = new TH1F(histTitle->Data(), histTitle->Data(), fPhiBins, fPhiMin, fPhiMax / order);
   fHistFull[k].fHistFullHar[j].fHistV0PhiCorr->Sumw2();
   fHistFull[k].fHistFullHar[j].fHistV0PhiCorr->SetXTitle("V0-Plane Correlation (rad)");
   fHistFull[k].fHistFullHar[j].fHistV0PhiCorr->SetYTitle("Counts");
   delete histTitle;

   // neutral sidebands-plane azimuthal correlation
   histTitle = new TString("FlowV0sb_Phi_Corr_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   fHistFull[k].fHistFullHar[j].fHistV0sbPhiCorr = new TH1F(histTitle->Data(), histTitle->Data(), fPhiBins, fPhiMin, fPhiMax / order);
   fHistFull[k].fHistFullHar[j].fHistV0sbPhiCorr->Sumw2();
   fHistFull[k].fHistFullHar[j].fHistV0sbPhiCorr->SetXTitle("V0sideBands-Plane Correlation (rad)");
   fHistFull[k].fHistFullHar[j].fHistV0sbPhiCorr->SetYTitle("Counts");
   delete histTitle;

   // Yield(pt)
   histTitle = new TString("Flow_YieldPt_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   fHistFull[k].fHistFullHar[j].fHistYieldPt = new TH1F(histTitle->Data(), histTitle->Data(), fPtBins, fPtMin, fPtMax);
   fHistFull[k].fHistFullHar[j].fHistYieldPt->Sumw2();
   fHistFull[k].fHistFullHar[j].fHistYieldPt->SetXTitle("Pt (GeV/c)");
   fHistFull[k].fHistFullHar[j].fHistYieldPt->SetYTitle("Yield");
   delete histTitle;
     
   // EtaPtPhi 
   histTitle = new TString("Flow_EtaPtPhi3D_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   fHistFull[k].fHistFullHar[j].fHistEtaPtPhi3D = new TH3F(histTitle->Data(), histTitle->Data(), fEtaBins, fEtaMin, fEtaMax, fPtBins, fPtMin, fPtMax, kPhi3DBins, fPhiMin, fPhiMax);
   fHistFull[k].fHistFullHar[j].fHistEtaPtPhi3D->SetXTitle("Eta");
   fHistFull[k].fHistFullHar[j].fHistEtaPtPhi3D->SetYTitle("Pt (GeV/c)");
   fHistFull[k].fHistFullHar[j].fHistEtaPtPhi3D->SetZTitle("Phi (rad)");
   delete histTitle;

   // Yield(eta,pt)
   histTitle = new TString("Flow_Yield2D_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   fHistFull[k].fHistFullHar[j].fHistYield2D = new TH2D(histTitle->Data(), histTitle->Data(), fEtaBins, fEtaMin, fEtaMax, fPtBins, fPtMin, fPtMax);
   fHistFull[k].fHistFullHar[j].fHistYield2D->Sumw2();
   fHistFull[k].fHistFullHar[j].fHistYield2D->SetXTitle("Pseudorapidty");
   fHistFull[k].fHistFullHar[j].fHistYield2D->SetYTitle("Pt (GeV/c)");
   delete histTitle;

   // Dca - 3D
   histTitle = new TString("Flow_3dDca_Sel") ;
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   fHistFull[k].fHistFullHar[j].fHistDcaGlob = new TH1F(histTitle->Data(), histTitle->Data(), kDcaBins, dcaMin, glDcaMax);
   fHistFull[k].fHistFullHar[j].fHistDcaGlob->Sumw2();
   fHistFull[k].fHistFullHar[j].fHistDcaGlob->SetXTitle("|3d Global Track's dca to Vertex (cm)|");
   delete histTitle;

   // Yield(pt) - excluded from R.P.
   histTitle = new TString("Flow_YieldPt_outSel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   fHistFull[k].fHistFullHar[j].fHistYieldPtout = new TH1F(histTitle->Data(), histTitle->Data(), fPtBins, fPtMin, fPtMax);
   fHistFull[k].fHistFullHar[j].fHistYieldPtout->Sumw2();
   fHistFull[k].fHistFullHar[j].fHistYieldPtout->SetXTitle("Pt (GeV/c)");
   fHistFull[k].fHistFullHar[j].fHistYieldPtout->SetYTitle("Yield");
   delete histTitle;
     
   // EtaPtPhi - excluded from R.P. 
   histTitle = new TString("Flow_EtaPtPhi3D_outSel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   fHistFull[k].fHistFullHar[j].fHistEtaPtPhi3Dout = new TH3F(histTitle->Data(), histTitle->Data(), fEtaBins, fEtaMin, fEtaMax, fPtBins, fPtMin, fPtMax, kPhi3DBins, fPhiMin, fPhiMax);
   fHistFull[k].fHistFullHar[j].fHistEtaPtPhi3Dout->SetXTitle("Eta");
   fHistFull[k].fHistFullHar[j].fHistEtaPtPhi3Dout->SetYTitle("Pt (GeV/c)");
   fHistFull[k].fHistFullHar[j].fHistEtaPtPhi3Dout->SetZTitle("Phi (rad)");
   delete histTitle;

   // Yield(eta,pt) - excluded from R.P.
   histTitle = new TString("Flow_Yield2D_outSel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   fHistFull[k].fHistFullHar[j].fHistYield2Dout = new TH2D(histTitle->Data(), histTitle->Data(), fEtaBins, fEtaMin, fEtaMax, fPtBins, fPtMin, fPtMax);
   fHistFull[k].fHistFullHar[j].fHistYield2Dout->Sumw2();
   fHistFull[k].fHistFullHar[j].fHistYield2Dout->SetXTitle("Pseudorapidty");
   fHistFull[k].fHistFullHar[j].fHistYield2Dout->SetYTitle("Pt (GeV/c)");
   delete histTitle;

   // Dca - 3D - excluded from R.P.
   histTitle = new TString("Flow_3dDca_outSel") ;
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   fHistFull[k].fHistFullHar[j].fHistDcaGlobout = new TH1F(histTitle->Data(), histTitle->Data(), kDcaBins, dcaMin, glDcaMax);
   fHistFull[k].fHistFullHar[j].fHistDcaGlobout->Sumw2();
   fHistFull[k].fHistFullHar[j].fHistDcaGlobout->SetXTitle("|3d Global Track's dca to Vertex (cm)|");
   delete histTitle;

   // Flow observed - v_obs,Pt,Eta
   histTitle = new TString("Flow_vObs2D_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   fHistFull[k].fHistFullHar[j].fHistvObs2D = new TProfile2D(histTitle->Data(), histTitle->Data(), fEtaBins, fEtaMin, fEtaMax, fPtBinsPart, fPtMin, fPtMaxPart, -100., 100., "");
   fHistFull[k].fHistFullHar[j].fHistvObs2D->SetXTitle((char*)fLabel.Data());
   fHistFull[k].fHistFullHar[j].fHistvObs2D->SetYTitle("Pt (GeV/c)");
   delete histTitle;

   // v_obs,Eta
   histTitle = new TString("Flow_vObsEta_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   fHistFull[k].fHistFullHar[j].fHistvObsEta = new TProfile(histTitle->Data(), histTitle->Data(), fEtaBins, fEtaMin, fEtaMax, -100., 100., "");
   fHistFull[k].fHistFullHar[j].fHistvObsEta->SetXTitle((char*)fLabel.Data());
   fHistFull[k].fHistFullHar[j].fHistvObsEta->SetYTitle("v (%)");
   delete histTitle;

   // v_obs,Pt
   histTitle = new TString("Flow_vObsPt_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   fHistFull[k].fHistFullHar[j].fHistvObsPt = new TProfile(histTitle->Data(), histTitle->Data(), fPtBinsPart, fPtMin, fPtMaxPart, -100., 100., "");
   fHistFull[k].fHistFullHar[j].fHistvObsPt->SetXTitle("Pt (GeV/c)");
   fHistFull[k].fHistFullHar[j].fHistvObsPt->SetYTitle("v (%)");
   delete histTitle;

   // neutral Flow observed - Pt,Eta
   histTitle = new TString("FlowV0_vObs2D_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   fHistFull[k].fHistFullHar[j].fHistV0vObs2D = new TProfile2D(histTitle->Data(), histTitle->Data(), fEtaBins, fEtaMin, fEtaMax, fPtBinsPart, fPtMin, fPtMaxPart, -100., 100., "");
   fHistFull[k].fHistFullHar[j].fHistV0vObs2D->SetXTitle((char*)fLabel.Data());
   fHistFull[k].fHistFullHar[j].fHistV0vObs2D->SetYTitle("Pt (GeV/c)");
   delete histTitle;

   // neutral Flow observed - Eta
   histTitle = new TString("FlowV0_vObsEta_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   fHistFull[k].fHistFullHar[j].fHistV0vObsEta = new TProfile(histTitle->Data(), histTitle->Data(), fEtaBins, fEtaMin, fEtaMax, -100., 100., "");
   fHistFull[k].fHistFullHar[j].fHistV0vObsEta->SetXTitle((char*)fLabel.Data());
   fHistFull[k].fHistFullHar[j].fHistV0vObsEta->SetYTitle("v (%)");
   delete histTitle;

   // neutral Flow observed - Pt
   histTitle = new TString("FlowV0_vObsPt_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   fHistFull[k].fHistFullHar[j].fHistV0vObsPt = new TProfile(histTitle->Data(), histTitle->Data(), fPtBinsPart, fPtMin, fPtMaxPart, -100., 100., "");
   fHistFull[k].fHistFullHar[j].fHistV0vObsPt->SetXTitle("Pt (GeV/c)");
   fHistFull[k].fHistFullHar[j].fHistV0vObsPt->SetYTitle("v (%)");
   delete histTitle;

   // neutral sidebands Flow observed - Pt,Eta
   histTitle = new TString("FlowV0sb_vObs2D_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   fHistFull[k].fHistFullHar[j].fHistV0sbvObs2D = new TProfile2D(histTitle->Data(), histTitle->Data(), fEtaBins, fEtaMin, fEtaMax, fPtBinsPart, fPtMin, fPtMaxPart, -100., 100., "");
   fHistFull[k].fHistFullHar[j].fHistV0sbvObs2D->SetXTitle((char*)fLabel.Data());
   fHistFull[k].fHistFullHar[j].fHistV0sbvObs2D->SetYTitle("Pt (GeV/c)");
   delete histTitle;

   // neutral sidebands Flow observed - Eta
   histTitle = new TString("FlowV0sb_vObsEta_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   fHistFull[k].fHistFullHar[j].fHistV0sbvObsEta = new TProfile(histTitle->Data(), histTitle->Data(), fEtaBins, fEtaMin, fEtaMax, -100., 100., "");
   fHistFull[k].fHistFullHar[j].fHistV0sbvObsEta->SetXTitle((char*)fLabel.Data());
   fHistFull[k].fHistFullHar[j].fHistV0sbvObsEta->SetYTitle("v (%)");
   delete histTitle;

   // neutral sidebands Flow observed - Pt
   histTitle = new TString("FlowV0sb_vObsPt_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   fHistFull[k].fHistFullHar[j].fHistV0sbvObsPt = new TProfile(histTitle->Data(), histTitle->Data(), fPtBinsPart, fPtMin, fPtMaxPart, -100., 100., "");
   fHistFull[k].fHistFullHar[j].fHistV0sbvObsPt->SetXTitle("Pt (GeV/c)");
   fHistFull[k].fHistFullHar[j].fHistV0sbvObsPt->SetYTitle("v (%)");
   delete histTitle;

   // SX neutral sidebands Flow observed - Eta
   histTitle = new TString("FlowV0sb_vObsEta_sx_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   fHistFull[k].fHistFullHar[j].fHistV0sbvObsEtaSx = new TProfile(histTitle->Data(), histTitle->Data(), fEtaBins, fEtaMin, fEtaMax, -100., 100., "");
   fHistFull[k].fHistFullHar[j].fHistV0sbvObsEtaSx->SetXTitle((char*)fLabel.Data());
   fHistFull[k].fHistFullHar[j].fHistV0sbvObsEtaSx->SetYTitle("v (%)");
   delete histTitle;

   // SX neutral sidebands Flow observed - Pt
   histTitle = new TString("FlowV0sb_vObsPt_sx_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   fHistFull[k].fHistFullHar[j].fHistV0sbvObsPtSx = new TProfile(histTitle->Data(), histTitle->Data(), fPtBinsPart, fPtMin, fPtMaxPart, -100., 100., "");
   fHistFull[k].fHistFullHar[j].fHistV0sbvObsPtSx->SetXTitle("Pt (GeV/c)");
   fHistFull[k].fHistFullHar[j].fHistV0sbvObsPtSx->SetYTitle("v (%)");
   delete histTitle;

   // DX neutral sidebands Flow observed - Eta
   histTitle = new TString("FlowV0sb_vObsEta_dx_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   fHistFull[k].fHistFullHar[j].fHistV0sbvObsEtaDx = new TProfile(histTitle->Data(), histTitle->Data(), fEtaBins, fEtaMin, fEtaMax, -100., 100., "");
   fHistFull[k].fHistFullHar[j].fHistV0sbvObsEtaDx->SetXTitle((char*)fLabel.Data());
   fHistFull[k].fHistFullHar[j].fHistV0sbvObsEtaDx->SetYTitle("v (%)");
   delete histTitle;

   // DX neutral sidebands Flow observed - Pt
   histTitle = new TString("FlowV0sb_vObsPt_dx_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   fHistFull[k].fHistFullHar[j].fHistV0sbvObsPtDx = new TProfile(histTitle->Data(), histTitle->Data(), fPtBinsPart, fPtMin, fPtMaxPart, -100., 100., "");
   fHistFull[k].fHistFullHar[j].fHistV0sbvObsPtDx->SetXTitle("Pt (GeV/c)");
   fHistFull[k].fHistFullHar[j].fHistV0sbvObsPtDx->SetYTitle("v (%)");
   delete histTitle;

   // Phi lab
   // Tpc (plus)
   histTitle = new TString("Flow_Phi_TPCplus_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   fHistFull[k].fHistFullHar[j].fHistPhiPlus = new TH1D(histTitle->Data(), histTitle->Data(), fPhiBins, fPhiMin, fPhiMax);
   fHistFull[k].fHistFullHar[j].fHistPhiPlus->SetXTitle("Azimuthal Angles (rad)");
   fHistFull[k].fHistFullHar[j].fHistPhiPlus->SetYTitle("Counts");
   delete histTitle;

   // Tpc (minus)
   histTitle = new TString("Flow_Phi_TPCminus_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   fHistFull[k].fHistFullHar[j].fHistPhiMinus = new TH1D(histTitle->Data(), histTitle->Data(), fPhiBins, fPhiMin, fPhiMax);
   fHistFull[k].fHistFullHar[j].fHistPhiMinus->SetXTitle("Azimuthal Angles (rad)");
   fHistFull[k].fHistFullHar[j].fHistPhiMinus->SetYTitle("Counts");
   delete histTitle;
  
   // Tpc (cross)
   histTitle = new TString("Flow_Phi_TPCcross_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   fHistFull[k].fHistFullHar[j].fHistPhiAll = new TH1D(histTitle->Data(), histTitle->Data(), fPhiBins, fPhiMin, fPhiMax);
   fHistFull[k].fHistFullHar[j].fHistPhiAll->SetXTitle("Azimuthal Angles (rad)");
   fHistFull[k].fHistFullHar[j].fHistPhiAll->SetYTitle("Counts");
   delete histTitle;
  
   // Tpc
   histTitle = new TString("Flow_Phi_TPC_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   fHistFull[k].fHistFullHar[j].fHistPhi = new TH1D(histTitle->Data(), histTitle->Data(), fPhiBins, fPhiMin, fPhiMax);
   fHistFull[k].fHistFullHar[j].fHistPhi->SetXTitle("Azimuthal Angles (rad)");
   fHistFull[k].fHistFullHar[j].fHistPhi->SetYTitle("Counts");
   delete histTitle;
  
   // Phi lab flattened
   // Tpc (Plus)
   histTitle = new TString("Flow_Phi_Flat_TPCplus_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   fHistFull[k].fHistFullHar[j].fHistPhiFlatPlus = new TH1D(histTitle->Data(), histTitle->Data(), fPhiBins, fPhiMin, fPhiMax);
   fHistFull[k].fHistFullHar[j].fHistPhiFlatPlus->SetXTitle("Azimuthal Angles (rad)");
   fHistFull[k].fHistFullHar[j].fHistPhiFlatPlus->SetYTitle("Counts");
   delete histTitle;

   // Tpc (Minus)
   histTitle = new TString("Flow_Phi_Flat_TPCminus_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   fHistFull[k].fHistFullHar[j].fHistPhiFlatMinus = new TH1D(histTitle->Data(), histTitle->Data(), fPhiBins, fPhiMin, fPhiMax);
   fHistFull[k].fHistFullHar[j].fHistPhiFlatMinus->SetXTitle("Azimuthal Angles (rad)");
   fHistFull[k].fHistFullHar[j].fHistPhiFlatMinus->SetYTitle("Counts");
   delete histTitle;

   // Tpc (cross)
   histTitle = new TString("Flow_Phi_Flat_TPCcross_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   fHistFull[k].fHistFullHar[j].fHistPhiFlatAll = new TH1D(histTitle->Data(), histTitle->Data(), fPhiBins, fPhiMin, fPhiMax);
   fHistFull[k].fHistFullHar[j].fHistPhiFlatAll->SetXTitle("Azimuthal Angles (rad)");
   fHistFull[k].fHistFullHar[j].fHistPhiFlatAll->SetYTitle("Counts");
   delete histTitle;

   // Tpc
   histTitle = new TString("Flow_Phi_Flat_TPC_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   fHistFull[k].fHistFullHar[j].fHistPhiFlat = new TH1D(histTitle->Data(), histTitle->Data(), fPhiBins, fPhiMin, fPhiMax);
   fHistFull[k].fHistFullHar[j].fHistPhiFlat->SetXTitle("Azimuthal Angles (rad)");
   fHistFull[k].fHistFullHar[j].fHistPhiFlat->SetYTitle("Counts");
   delete histTitle;
  }
 }
 cout << "    Init()  -  Histograms booked" << endl ; 

 return kTRUE ;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
Bool_t AliFlowAnalyser::Finish() 
{
 // Close the analysis and saves the histograms on the histFile .
 
 cout << "* FlowAnalysis *  -  Finish()" << endl ; cout << endl ;

 // Write all histograms
 fHistFile->cd() ; 
 fHistFile->Write() ; 
 
 // Write Resolution corrected histograms
 if(fVnResHistList)
 { 
  fVnResHistList->Write(); 
  delete fVnResHistList ;
 }
 else { cout << "  E.P. resolution has not been calculated. No v_n histograms!" << endl ; } 

 // Write PhiWgt histograms
 if(fPhiWgtHistList)
 { 
  fPhiWgtHistList->Write();
  delete fPhiWgtHistList ; 
 }
 
 fFlowSelect->Write(); 
 // delete fFlowSelect ;

 fHistFile->Close() ;

 cout << endl ;
 cout << "    Finish()  -  tot. " << fTotalNumberOfEvents << " have been analyzed (with " << fTotalNumberOfTracks << " tracks and " << fTotalNumberOfV0s << " v0s) ."  << endl ; 
 cout << endl ; 
 cout << "    Finish()  -  Histograms saved : " << fHistFileName.Data() << endl ; cout << endl ; 

 return kTRUE ;
}
//-----------------------------------------------------------------------
// ###
//----------------------------------------------------------------------
Float_t AliFlowAnalyser::GetRunBayesian(Int_t nPid, Int_t selN) const 
{
 // Returns the normalized particle abundance of "e","mu","pi","k","p","d"
 // in all the analysed events (in selection selN).
 // Call at the end of the analysis.

 if(selN>AliFlowConstants::kSels) { selN = 0 ; }
 Double_t totCount = (fHistFull[selN].fHistBayPidMult)->GetSumOfWeights() ;
 if(totCount) { return (fHistFull[selN].fHistBayPidMult->GetBinContent(nPid+1) / totCount) ; }
 else         { return 1. ; }
}
//-----------------------------------------------------------------------
void AliFlowAnalyser::PrintRunBayesian(Int_t selN) const
{
 // Prints the normalized particle abundance of all the analysed events 
 // (in selection selN).

 if(selN>AliFlowConstants::kSels) { selN = 0 ; } 
 Char_t* names[AliFlowConstants::kPid] = {"e","mu","pi","k","p","d"} ;
 Double_t bayes = 0. ;
 cout << " selN = " << selN << " particles normalized abundance : " ;
 for(int i=0;i<AliFlowConstants::kPid;i++)
 {
  bayes = GetRunBayesian(i, selN) ;
  cout << bayes << "_" << names[i] << " ; " ;
 }
 cout << endl ;
 
 return ;
}
//-----------------------------------------------------------------------
void AliFlowAnalyser::FillWgtArrays(TFile* wgtFile)
{
 // Loads PhiWeights & Bayesian particles' abundance from file (default: flowPhiWgt.hist.root). 
 // Weights are stored in a static TH1D* data member, ready to be plugged into the AliFlowEvent.
 // The plugging is done by the method ::FillEvtPhiWgt() (if wgt file is there).

 fPhiWgtFile = wgtFile ;

 TString* histTitle ;
 TH1D* hTPCall ; TH1D* hTPCplus ; TH1D* hTPCminus ; TH1D* hTPCcross ;
 TH1D* hPIDbay ;
 for(int k=0;k<AliFlowConstants::kSels;k++)
 {
  for(int j=0;j<AliFlowConstants::kHars;j++) 
  {
  // Tpc (plus)
   histTitle = new TString("Flow_Phi_Weight_TPCplus_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   hTPCplus = (TH1D*)fPhiWgtFile->Get(histTitle->Data());
   delete histTitle;
   // Tpc (minus)
   histTitle = new TString("Flow_Phi_Weight_TPCminus_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   hTPCminus = (TH1D*)fPhiWgtFile->Get(histTitle->Data());
   delete histTitle;
   // Tpc (cross)
   histTitle = new TString("Flow_Phi_Weight_TPCcross_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   hTPCcross = (TH1D*)fPhiWgtFile->Get(histTitle->Data());
   delete histTitle;

   // Tpc
   histTitle = new TString("Flow_Phi_Weight_TPC_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   hTPCall = (TH1D*)fPhiWgtFile->Get(histTitle->Data());
   delete histTitle;

   for(int n=0;n<fPhiBins;n++) 
   {
    fPhiWgtPlus[k][j][n]  = hTPCplus->GetBinContent(n+1) ;
    fPhiWgtMinus[k][j][n] = hTPCminus->GetBinContent(n+1) ;
    fPhiWgtCross[k][j][n] = hTPCcross->GetBinContent(n+1) ;
    fPhiWgt[k][j][n]	  = hTPCall->GetBinContent(n+1) ;
    // cout << "  Weights: " << fPhiWgt[k][j][n] << " ; " << fPhiWgtPlus[k][j][n] << " | " << fPhiWgtMinus[k][j][n] << " | " << fPhiWgtCross[k][j][n] << endl ; 
   }
  }
  
  // Bayesian weights
  histTitle = new TString("Flow_BayPidMult_Sel");
  *histTitle += k+1;
  hPIDbay = (TH1D*)fPhiWgtFile->Get(histTitle->Data());
  delete histTitle;
  Double_t totCount = hPIDbay->GetSumOfWeights() ;
  for (int n=0;n<AliFlowConstants::kPid;n++) 
  { 
   if(totCount) { fBayesianWgt[k][n] = hPIDbay->GetBinContent(n+1) / totCount ; }
   else 	{ fBayesianWgt[k][n] = 1. ; }
   // cout << "  Bayesian Weights (" << n << ") : " << fBayesianWgt[k][n] << endl ; 
  }
 }

 delete hTPCall ; delete hTPCplus ; delete hTPCminus ; delete hTPCcross ;
 delete hPIDbay ;

 return ;
}
//-----------------------------------------------------------------------
void AliFlowAnalyser::FillEvtPhiWgt(AliFlowEvent* fFlowEvent)
{
 // Plugs phi weights into the static dwgt data member of the AliFlowEvent class.
 // Weights are given in special AliFlowConstants::PhiWgt_t arrays (see AliFlowConstants), 
 // which are read from the wgt histograms by the method FillWgtArrays(...).
 
 fFlowEvent->SetPhiWeight(fPhiWgt);
 fFlowEvent->SetPhiWeightPlus(fPhiWgtPlus);
 fFlowEvent->SetPhiWeightMinus(fPhiWgtMinus);
 fFlowEvent->SetPhiWeightCross(fPhiWgtCross); 
 
 return ;
}
//-----------------------------------------------------------------------
void AliFlowAnalyser::FillBayesianWgt(AliFlowEvent* fFlowEvent)
{
 // Plugs Bayesian particle abundance into the current AliFlowEvent.
 // A bayesian vector should be used for the PId of any different selection 
 // (different sets of cuts -> different particle abundance), but for now  
 // just the Selection n.0 (with no cuts) is used .
 // (AliFlowEvent::mBayesianCs[6] is a 1-dimensional array, change that first!).

 Double_t bayes[AliFlowConstants::kPid] ; 
 Double_t bayCheck = 0. ;
 for (int n=0;n<AliFlowConstants::kPid;n++) 
 {
  bayes[n] = fBayesianWgt[0][n] ;
  bayCheck += bayes[n] ;   
  // cout << "Bayesian V[" << n << "]  =  " << fBayesianWgt[0][n] << endl ; 
 }
 if(bayCheck)   { fFlowEvent->SetBayesian(bayes) ; fRePid = kTRUE ; }
 else		{ cout << "An empty bayesian vector is stored !!! - Bayesian weights = {1,1,1,1,1,1} " << endl ; }

 return ;
}
//-----------------------------------------------------------------------
void AliFlowAnalyser::Weightening()
{
 // Calculates weights, and fills PhiWgt histograms .
 // This is called at the end of the event analysis.
 
 cout << " AliFlowAnalyser::Weightening() " << endl ; cout << endl ;
 
 // PhiWgt histogram collection
 fPhiWgtHistList = new TOrdCollection(4*AliFlowConstants::kSels*AliFlowConstants::kHars) ;
 
 // Creates PhiWgt Histograms
 TString* histTitle ;
 for(int k = 0; k < AliFlowConstants::kSels; k++)
 {
  for(int j = 0; j < AliFlowConstants::kHars; j++) 
  {
   // Tpc (plus)
   histTitle = new TString("Flow_Phi_Weight_TPCplus_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   fHistFull[k].fHistFullHar[j].fHistPhiWgtPlus = new TH1D(histTitle->Data(),histTitle->Data(), fPhiBins, fPhiMin, fPhiMax);
   fHistFull[k].fHistFullHar[j].fHistPhiWgtPlus->Sumw2();
   fHistFull[k].fHistFullHar[j].fHistPhiWgtPlus->SetXTitle("Azimuthal Angles (rad)");
   fHistFull[k].fHistFullHar[j].fHistPhiWgtPlus->SetYTitle("PhiWgt");
   delete histTitle;
   // Tpc (minus)
   histTitle = new TString("Flow_Phi_Weight_TPCminus_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   fHistFull[k].fHistFullHar[j].fHistPhiWgtMinus = new TH1D(histTitle->Data(),histTitle->Data(), fPhiBins, fPhiMin, fPhiMax);
   fHistFull[k].fHistFullHar[j].fHistPhiWgtMinus->Sumw2();
   fHistFull[k].fHistFullHar[j].fHistPhiWgtMinus->SetXTitle("Azimuthal Angles (rad)");
   fHistFull[k].fHistFullHar[j].fHistPhiWgtMinus->SetYTitle("PhiWgt");
   delete histTitle;
   // Tpc (cross)
   histTitle = new TString("Flow_Phi_Weight_TPCcross_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   fHistFull[k].fHistFullHar[j].fHistPhiWgtAll = new TH1D(histTitle->Data(),histTitle->Data(), fPhiBins, fPhiMin, fPhiMax);
   fHistFull[k].fHistFullHar[j].fHistPhiWgtAll->Sumw2();
   fHistFull[k].fHistFullHar[j].fHistPhiWgtAll->SetXTitle("Azimuthal Angles (rad)");
   fHistFull[k].fHistFullHar[j].fHistPhiWgtAll->SetYTitle("PhiWgt");
   delete histTitle;
   // Tpc
   histTitle = new TString("Flow_Phi_Weight_TPC_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   fHistFull[k].fHistFullHar[j].fHistPhiWgt = new TH1D(histTitle->Data(),histTitle->Data(), fPhiBins, fPhiMin, fPhiMax);
   fHistFull[k].fHistFullHar[j].fHistPhiWgt->Sumw2();
   fHistFull[k].fHistFullHar[j].fHistPhiWgt->SetXTitle("Azimuthal Angles (rad)");
   fHistFull[k].fHistFullHar[j].fHistPhiWgt->SetYTitle("PhiWgt");
   delete histTitle;

   // Calculate PhiWgt
   double meanPlus  = fHistFull[k].fHistFullHar[j].fHistPhiPlus->Integral() / (double)fPhiBins ;
   double meanMinus = fHistFull[k].fHistFullHar[j].fHistPhiMinus->Integral() / (double)fPhiBins ;
   double meanCross = fHistFull[k].fHistFullHar[j].fHistPhiAll->Integral() / (double)fPhiBins ;
   double meanTPC = fHistFull[k].fHistFullHar[j].fHistPhi->Integral() / (double)fPhiBins ;

   // Tpc
   for (int i=0;i<fPhiBins;i++) 
   {
    fHistFull[k].fHistFullHar[j].fHistPhiWgtPlus->SetBinContent(i+1,meanPlus);
    fHistFull[k].fHistFullHar[j].fHistPhiWgtPlus->SetBinError(i+1, 0.);
    fHistFull[k].fHistFullHar[j].fHistPhiWgtMinus->SetBinContent(i+1,meanMinus);
    fHistFull[k].fHistFullHar[j].fHistPhiWgtMinus->SetBinError(i+1, 0.);
    fHistFull[k].fHistFullHar[j].fHistPhiWgtAll->SetBinContent(i+1,meanCross);
    fHistFull[k].fHistFullHar[j].fHistPhiWgtAll->SetBinError(i+1, 0.);
    fHistFull[k].fHistFullHar[j].fHistPhiWgt->SetBinContent(i+1,meanTPC);
    fHistFull[k].fHistFullHar[j].fHistPhiWgt->SetBinError(i+1, 0.);
   }
   
   if(meanTPC==0) { cout << " Sel." << k << " , Har." << j << " :  empty phi histogram ! " << endl ; }
   else 
   {
    fHistFull[k].fHistFullHar[j].fHistPhiWgtPlus->Divide(fHistFull[k].fHistFullHar[j].fHistPhiPlus);
    fHistFull[k].fHistFullHar[j].fHistPhiWgtMinus->Divide(fHistFull[k].fHistFullHar[j].fHistPhiMinus);
    fHistFull[k].fHistFullHar[j].fHistPhiWgtAll->Divide(fHistFull[k].fHistFullHar[j].fHistPhiAll);
    fHistFull[k].fHistFullHar[j].fHistPhiWgt->Divide(fHistFull[k].fHistFullHar[j].fHistPhi);
   }
   
   fPhiWgtHistList->AddLast(fHistFull[k].fHistFullHar[j].fHistPhiWgtPlus);
   fPhiWgtHistList->AddLast(fHistFull[k].fHistFullHar[j].fHistPhiWgtMinus);
   fPhiWgtHistList->AddLast(fHistFull[k].fHistFullHar[j].fHistPhiWgtAll);
   fPhiWgtHistList->AddLast(fHistFull[k].fHistFullHar[j].fHistPhiWgt);
  }
 }

 return ;
}
//-----------------------------------------------------------------------
// ###
//-----------------------------------------------------------------------
Bool_t AliFlowAnalyser::Analyze(AliFlowEvent* flowEvent)         
{
 // Runs the analysis on the AliFlowEvent (* fFlowEvent). 
 // This method can be inserted in a loop over a collection of 
 // AliFlowEvents or for on-fly analysis on AliESDs if used toghether 
 // with the AliFlowMaker.
 
 cout << " AliFlowAnalyser::Analyze(" << fFlowEvent << " )   -   " << fEventNumber << endl ;
 if(!flowEvent) { return kFALSE ; }
 else  { fFlowEvent = flowEvent ; }

 if(fFlowSelect->Select(fFlowEvent))	 // event selected - here below the ANALYSIS FLAGS are setted -
 {
  cout << " * 1 . Load event (track & v0s) and set flags . " << endl ;
  fFlowTracks = fFlowEvent->TrackCollection() ; 
  fNumberOfTracks = fFlowTracks->GetEntries() ;
  fFlowV0s = fFlowEvent->V0Collection() ; 
  fNumberOfV0s = fFlowV0s->GetEntries() ;
  cout << "       event ID = " << fFlowEvent->EventID() << " :  found " << fNumberOfTracks << " AliFlowTracks, and " << fNumberOfV0s << " AliFlowV0s . " << endl ;  
  	    						    
  if(fReadPhiWgt) { FillEvtPhiWgt(fFlowEvent) ; }	    // phi and bayesian weights are filled previous to the loop (FillWgtArrays(TFile*))
  else 	    	  { fFlowEvent->SetNoWgt() ; }		    // phi weights can be used or not , this plugs them into the event

  fFlowEvent->SetCustomRespFunc(fCustomRespFunc) ;	    // a custom "detector response function" is used for assigning P.Id

  if(fBayWgt)	  { FillBayesianWgt(fFlowEvent) ; }   	    // bayesian weights can be used or not , this plugs them into the event
  if(fRePid)	  { fFlowEvent->SetPids() ; }		    // re-calculate all p.id. hypotesis with the (new) bayesian array

  if(fOnePhiWgt)  { fFlowEvent->SetOnePhiWgt() ; }	    // one phi-wgt histogram
  else 	  	  { fFlowEvent->SetFirstLastPhiWgt() ; }    // three phi-wgt histogram
  if(fPtWgt)	  { fFlowEvent->SetPtWgt(); ; } 	    // pT as a weight
  if(fEtaWgt)	  { fFlowEvent->SetEtaWgt() ; } 	    // eta as a weight

  if(fShuffle)    { fFlowEvent->RandomShuffle() ; }	    // tracks re-shuffling

  fFlowEvent->SetSelections(fFlowSelect) ;		    // does the selection of tracks for r.p. calculation (sets flags in AliFlowTrack)
  
  if(fSub == 1)       { fFlowEvent->SetEtaSubs() ; }	    // setting for the subevents (eta, random, charged)
  else if(fSub == -1) { fFlowEvent->SetChrSubs() ; }	    // 
  else 		      { fFlowEvent->SetRndSubs() ; }	    // ...
  
  fEtaSub = fFlowEvent->EtaSubs() ;

  fFlowEvent->MakeSubEvents() ; 			    // makes the subevent, eta or random basing on the previous flag

  cout << " * 2 . Calculating event quantities all in one shoot . " << endl ; 
  fFlowEvent->MakeAll() ; 
  
  if(FillFromFlowEvent(fFlowEvent))			    // calculates event quantities
  {
   cout << " * 3 . Event Histograms and Particles loop . " << endl ;
   FillEventHistograms(fFlowEvent);			    // fill histograms from AliFlowEvents
   if(fTrackLoop) { FillParticleHistograms(fFlowTracks) ; } // fill histograms from AliFlowTracks
   if(fV0loop)    { FillV0Histograms(fFlowV0s) ; }	    // fill histograms from AliFlowV0s 
   //FillLabels() ; 		  	 		    // fill the histogram of MC labels (from the simulation)
  }
  else 
  {
   cout << " * 3 . Event psi = 0   -   Skipping! " << endl ; 
   return kFALSE ;
  }
 }
 else 
 {
  cout << " * 0 . Event " << fEventNumber << " (event ID = " << fFlowEvent->EventID() << ") discarded . " << endl ; 
  // delete fFlowEvent  ; fFlowEvent = 0 ; 
  return kFALSE ;
 }
 fEventNumber++ ;
 
 return kTRUE ;
}
//-----------------------------------------------------------------------
Bool_t AliFlowAnalyser::FillFromFlowEvent(AliFlowEvent* fFlowEvent)
{
 // gets event quantities, returns kFALSE if Q vector is always 0 

 Int_t selCheck = 0 ;
 for(Int_t k = 0; k < AliFlowConstants::kSels; k++) 
 {
  fFlowSelect->SetSelection(k) ;
  for(Int_t j = 0; j < AliFlowConstants::kHars; j++) 
  {
   fFlowSelect->SetHarmonic(j) ;
   for(Int_t n = 0; n < AliFlowConstants::kSubs; n++) 
   {
    fFlowSelect->SetSubevent(n) ;
    fPsiSub[n][k][j] = fFlowEvent->Psi(fFlowSelect) ;  	// sub-event quantities
    fMultSub[n][k][j] = fFlowEvent->Mult(fFlowSelect) ;
   }
   fFlowSelect->SetSubevent(-1);
   fQ[k][j]    = fFlowEvent->Q(fFlowSelect) ;	   	// full event quantities
   fPsi[k][j]  = fFlowEvent->Psi(fFlowSelect) ;
   fQnorm[k][j]   = fFlowEvent->NormQ(fFlowSelect).Mod() ; // was: fFlowEvent->q(fFlowSelect) ; // but the normalization was bad (no pT,eta weight)
   fMult[k][j] = fFlowEvent->Mult(fFlowSelect) ;
   selCheck += fMult[k][j] ; 
  }
 }
 if(!selCheck) { return kFALSE ; }           		// if there are no particles in the selection -> skip the event

 return kTRUE ;
}
//-----------------------------------------------------------------------
void AliFlowAnalyser::FillEventHistograms(AliFlowEvent* fFlowEvent)    
{
 // fill event histograms

 cout << " Fill Event Histograms ... " << endl ; 

 float trigger = (float)fFlowEvent->L0TriggerWord() ;
 fHistTrigger->Fill(trigger);

 // no selections: OrigMult, Centrality, Mult, MultOverOrig, VertexZ, VertexXY
 int origMult = fFlowEvent->OrigMult();
 fHistOrigMult->Fill((float)origMult);
 fHistMultEta->Fill((float)fFlowEvent->MultEta());

 int cent = fFlowEvent->Centrality();
 fHistCent->Fill((float)cent);

 fHistMult->Fill((float)fNumberOfTracks) ;
 fHistV0Mult->Fill((float)fNumberOfV0s) ;
 if(origMult) { fHistMultOverOrig->Fill((float)fNumberOfTracks/(float)origMult) ; }

 fFlowEvent->VertexPos(fVertex) ;
 fHistVertexZ->Fill(fVertex[2]) ;
 fHistVertexXY2D->Fill(fVertex[0],fVertex[1]) ;

 // ZDC info
 fHistPartZDC->Fill(fFlowEvent->ZDCpart()) ;
 for(int ii=0;ii<3;ii++) { fHistEnergyZDC->Fill(ii,fFlowEvent->ZDCenergy(ii)) ; }

 // sub-event Psi_Subs
 for(int k = 0; k < AliFlowConstants::kSels; k++) 
 {
  for(int j = 0; j < AliFlowConstants::kHars; j++) 
  {
   for(int n = 0; n < AliFlowConstants::kSubs; n++) 
   {
    int iii = AliFlowConstants::kSubs * k + n ;    //cout << "  " << k << j << n << " , " << iii << endl ;
    fHistSub[iii].fHistSubHar[j].fHistPsiSubs->Fill(fPsiSub[n][k][j]) ;
   }
  }
 }

 // full event Psi, PsiSubCorr, PsiSubCorrDiff, cos, mult, q
 for(int k = 0; k < AliFlowConstants::kSels; k++) 
 {
  for(int j = 0; j < AliFlowConstants::kHars; j++) 
  {
   float order = (float)(j+1);
   fHistFull[k].fHistFullHar[j].fHistPsi->Fill(fPsi[k][j]);
   if(k<2 && j<2)
   {
    if(k==0) 
    { 
     float psi1 = fPsi[0][j] ; 
     float psi2 = fPsi[1][j] ;
     float diff = psi1 - psi2 ;
     if(diff < -TMath::Pi()/(j+1))      { diff += 2*TMath::Pi()/(j+1) ; } 
     else if(diff > +TMath::Pi()/(j+1)) { diff -= 2*TMath::Pi()/(j+1) ; }
     fHistFull[k].fHistFullHar[j].fHistPsiDiff->Fill(diff) ; // k=0
    } 
    else if(k==1) 
    {
     float psi1 = 0. ; float psi2 = 0. ;
     if (j==0)     { psi1 = fPsi[0][0] ; psi2 = fPsi[1][1] ; }  
     else if(j==1) { psi1 = fPsi[0][0] ; psi2 = fPsi[0][1] ; }
     float diff = psi1 - psi2 ;
     diff = (TMath::Abs(diff) > TMath::Pi()) ? ((diff > 0.) ? -(2*TMath::Pi()-diff) : -(diff+2*TMath::Pi())) : diff ;
     fHistFull[k].fHistFullHar[j].fHistPsiDiff->Fill(diff) ; // k=1
    }	   
   }

   if(fPsiSub[0][k][j] != 0. && fPsiSub[1][k][j] != 0.)
   {
    float psiSubCorr;    			// this is:  delta_Psi
    if(fV1Ep1Ep2 == kFALSE || order != 1) 
    {
     psiSubCorr = fPsiSub[0][k][j] - fPsiSub[1][k][j];
    }
    else // i.e. (fV1Ep1Ep2 == kTRUE && order == 1)
    { 
     psiSubCorr = fPsiSub[0][k][0] + fPsiSub[1][k][0] - 2*fPsi[k][1];
    }
    fHistFull[k].fHistCos->Fill(order, (float)cos(order * psiSubCorr)) ;
    if(psiSubCorr < 0.) 		 { psiSubCorr += 2*TMath::Pi()/order ; }
    if(psiSubCorr > 2*TMath::Pi()/order) { psiSubCorr -= 2*TMath::Pi()/order ; } // for v1Ep1Ep2 which gives -2*TMath::Pi() < psiSubCorr < 2*2*TMath::Pi()
    fHistFull[k].fHistFullHar[j].fHistPsiSubCorr->Fill(psiSubCorr);
   }

   if(j < AliFlowConstants::kHars - 1) // subevents of different harmonics
   {
    int j1 = 0 ; int j2 = 0 ;
    float psiSubCorrDiff;
    if(j==0)	  { j1 = 1, j2 = 2 ; } 
    else if(j==1) { j1 = 1, j2 = 3 ; } 
    else if(j==2) { j1 = 2, j2 = 4 ; }
    psiSubCorrDiff = fmod((double)fPsiSub[0][k][j1-1],2*TMath::Pi()/(double)j2)-fmod((double)fPsiSub[1][k][j2-1],2*TMath::Pi()/(double)j2) ;
    if(psiSubCorrDiff < 0.) { psiSubCorrDiff += 2*TMath::Pi()/(float)j2 ; }
    fHistFull[k].fHistFullHar[j].fHistPsiSubCorrDiff->Fill(psiSubCorrDiff) ;
    psiSubCorrDiff = fmod((double)fPsiSub[0][k][j2-1],2*TMath::Pi()/(double)j2)-fmod((double)fPsiSub[1][k][j1-1],2*TMath::Pi()/(double)j2) ;
    if(psiSubCorrDiff < 0.) { psiSubCorrDiff += 2*TMath::Pi()/(float)j2 ; }
    fHistFull[k].fHistFullHar[j].fHistPsiSubCorrDiff->Fill(psiSubCorrDiff) ;
   }
   
   fHistFull[k].fHistFullHar[j].fHistMult->Fill((float)fMult[k][j]) ;
   fHistFull[k].fHistFullHar[j].fHistQnorm->Fill(fQnorm[k][j]) ;
  }
 } 
 fTotalNumberOfEvents++ ;
 
 return ;
}
//-----------------------------------------------------------------------
void AliFlowAnalyser::FillParticleHistograms(TClonesArray* fFlowTracks) 
{
 // fills tracks histograms

 cout << " Tracks Loop . " << endl ; 

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

 fSelParts = 0 ;
 fSelV0s   = 0 ;
 for(fTrackNumber=0;fTrackNumber<fNumberOfTracks;fTrackNumber++) 
 {
  //fFlowTrack = (AliFlowTrack*)fFlowTracks->At(fTrackNumber) ;
  
  fFlowTrack = (AliFlowTrack*)(fFlowTracks->AddrAt(fTrackNumber)) ;
  //cout << "Track n. " << fTrackNumber << endl ; // fFlowTrack->Dump() ; 
  
  bool constrainable = fFlowTrack->IsConstrainable() ;
  // int label = fFlowTrack->Label() ; 			
  float dcaGlobal = TMath::Abs(fFlowTrack->Dca()) ; 	
  float dcaSigned = fFlowTrack->TransDcaSigned() ;	
  float dcaTrans = TMath::Abs(dcaSigned) ; 		
  float eta = fFlowTrack->Eta() ; 			
  float phi = fFlowTrack->Phi() ; 			
  float pt = fFlowTrack->Pt() ; 			
  float etaGlob = fFlowTrack->EtaGlobal() ; 		 							
  float phiGlob = fFlowTrack->PhiGlobal() ; 		
  float ptGlob = fFlowTrack->PtGlobal() ; 		
  float totalp = fFlowTrack->P() ; 			
  // float logp = TMath::Log10(totalp) ;
  // float zFirstPoint = fFlowTrack->ZFirstPoint() ; 	
  // float zLastPoint = fFlowTrack->ZLastPoint() ; 	
  float lenght = fFlowTrack->TrackLength() ;     	
  int	charge = fFlowTrack->Charge() ; 		
  float chi2 = fFlowTrack->Chi2() ; 			
  int	fitPtsTPC = (int)((float)fFlowTrack->FitPtsTPC()) ;
  int	maxPtsTPC = fFlowTrack->MaxPtsTPC() ;		
  float chi2TPC = fFlowTrack->Chi2TPC() ; 		
  int	fitPtsITS = fFlowTrack->FitPtsITS() ;		
  int	maxPtsITS = fFlowTrack->MaxPtsITS() ;	        
  float chi2ITS = fFlowTrack->Chi2ITS() ;	        
  int	fitPtsTRD = fFlowTrack->NhitsTRD() ;		
  int	maxPtsTRD = fFlowTrack->MaxPtsTRD() ;		
  float chi2TRD = fFlowTrack->Chi2TRD() ;	        
  int	fitPtsTOF = fFlowTrack->NhitsTOF() ;		
  int	maxPtsTOF = fFlowTrack->MaxPtsTOF() ;		
  float chi2TOF = fFlowTrack->Chi2TOF() ;	        
  float dEdx = fFlowTrack->DedxTPC() ; 			
  float its = fFlowTrack->DedxITS() ;		        
  float trd = fFlowTrack->SigTRD() ;		        
  float tof = fFlowTrack->TofTOF() ;
  float lpTPC = 0 ; if(fFlowTrack->PatTPC()>0) { lpTPC = TMath::Log10(fFlowTrack->PatTPC()) ; }
  float lpITS = 0 ; if(fFlowTrack->PatITS()>0) { lpITS = TMath::Log10(fFlowTrack->PatITS()) ; }
  float lpTRD = 0 ; if(fFlowTrack->PatTRD()>0) { lpTRD = TMath::Log10(fFlowTrack->PatTRD()) ; }
  float lpTOF = 0 ; if(fFlowTrack->PatTOF()>0) { lpTOF = TMath::Log10(fFlowTrack->PatTOF()) ; }  
  float invMass = fFlowTrack->InvMass() ;	       
  Char_t pid[10]="0" ; strcpy(pid,fFlowTrack->Pid()) ; 
  fPidId        = -1 ; // assigned later

  // no selections: Charge, Dca, Chi2, FitPts, MaxPts, FitOverMax, PID
  fHistCharge->Fill((float)charge);
  fHistDcaGlobal->Fill(dcaGlobal);
  fHistDca->Fill(dcaTrans) ;
  fHistTransDca->Fill(dcaSigned);
  fHistChi2->Fill(chi2);
  
  // - here TPC   (chi2 & nHits)
  fHistChi2TPC->Fill(chi2TPC);
  if(fitPtsTPC>0) { fHistChi2normTPC->Fill(chi2TPC/((float)fitPtsTPC)) ; }
  fHistFitPtsTPC->Fill((float)fitPtsTPC);
  fHistMaxPtsTPC->Fill((float)maxPtsTPC);
  if(maxPtsTPC>0) { fHistFitOverMaxTPC->Fill((float)(fitPtsTPC)/(float)maxPtsTPC) ; }

  // - here ITS   (chi2 & nHits)
  fHistChi2ITS->Fill(chi2ITS);
  if(fitPtsITS>0) { fHistChi2normITS->Fill(chi2ITS/((float)fitPtsITS)) ; }
  fHistFitPtsITS->Fill((float)fitPtsITS);
  fHistMaxPtsITS->Fill((float)maxPtsITS);

  // - here TRD   (chi2 & nHits)
  fHistChi2TRD->Fill(chi2TRD);
  if(fitPtsTRD>0) { fHistChi2normTRD->Fill(chi2TRD/((float)fitPtsTRD)) ; }
  fHistFitPtsTRD->Fill((float)fitPtsTRD);
  fHistMaxPtsTRD->Fill((float)maxPtsTRD);

  // - here TOF   (chi2 & nHits)
  fHistChi2TOF->Fill(chi2TOF);
  if(fitPtsTOF>0) { fHistChi2normTOF->Fill(chi2TOF/((float)fitPtsTOF)) ; }
  fHistFitPtsTOF->Fill((float)fitPtsTOF);
  fHistMaxPtsTOF->Fill((float)maxPtsTOF);
  
  // fit over max (all)
  int maxPts = maxPtsITS + maxPtsTPC + maxPtsTRD + maxPtsTOF ;
  int fitPts = fitPtsITS + fitPtsTPC + fitPtsTRD + fitPtsTOF ;
  if(maxPts>0)    { fHistFitOverMax->Fill((float)(fitPts)/(float)maxPts) ; }
  
  // lenght
  fHistLenght->Fill(lenght) ;
  fHistInvMass->Fill(invMass) ;

  // PID histograms & multiplicity count (for bayesian histogram)
  if(charge == 1) 
  {
   hPlusN++ ;
   fHistMeanDedxPos2D->Fill(lpTPC, dEdx) ;
   fHistMeanDedxPos2DITS->Fill(lpITS, its) ;
   fHistMeanDedxPos2DTRD->Fill(lpTRD, trd) ;
   fHistMeanDedxPos2DTOF->Fill(lpTOF, tof) ;
   //
   float positron  = fFlowTrack->ElectronPositronProb() ;
   fHistPidPositron->Fill(positron) ;
   if(strcmp(pid, "e+") == 0) 
   {
    fPidId = 0 ; positronN++ ; fHistPidPt->Fill(1.5,pt) ;
    fHistMeanTPCPositron->Fill(lpTPC, dEdx) ;
    fHistMeanITSPositron->Fill(lpITS, its);
    fHistMeanTRDPositron->Fill(lpTRD, trd);
    fHistMeanTOFPositron->Fill(invMass, tof);
    fHistPidPositronPart->Fill(positron) ;
   }
   float muonPlus  = fFlowTrack->MuonPlusMinusProb() ;
   fHistPidMuonPlus->Fill(muonPlus) ;
   if(strcmp(pid, "mu+") == 0) 
   {
    fPidId = 1 ; muonPlusN++ ; fHistPidPt->Fill(3.5,pt) ; 
    fHistMeanTPCMuonPlus->Fill(lpTPC, dEdx) ;
    fHistMeanITSMuonPlus->Fill(lpITS, its);
    fHistMeanTRDMuonPlus->Fill(lpTRD, trd);
    fHistMeanTOFMuonPlus->Fill(invMass, tof);
    fHistPidMuonPlusPart->Fill(muonPlus) ;
   }
   float piPlus = fFlowTrack->PionPlusMinusProb() ;
   fHistPidPiPlus->Fill(piPlus) ;
   if(strcmp(pid, "pi+") == 0) 
   { 
    fPidId = 2 ; piPlusN++ ; fHistPidPt->Fill(5.5,pt) ;
    fHistMeanTPCPiPlus->Fill(lpTPC, dEdx) ;
    fHistMeanITSPiPlus->Fill(lpITS, its);
    fHistMeanTRDPiPlus->Fill(lpTRD, trd);
    fHistMeanTOFPiPlus->Fill(invMass, tof);
    fHistPidPiPlusPart->Fill(piPlus) ;
   }
   float kplus  = fFlowTrack->KaonPlusMinusProb() ;
   fHistPidKplus->Fill(kplus) ;
   if(strcmp(pid, "k+") == 0) 
   {
    fPidId = 3 ; kPlusN++ ; fHistPidPt->Fill(7.5,pt) ; 
    fHistMeanTPCKplus->Fill(lpTPC, dEdx) ;
    fHistMeanITSKplus->Fill(lpITS, its);
    fHistMeanTRDKplus->Fill(lpTRD, trd);
    fHistMeanTOFKplus->Fill(invMass, tof);
    fHistPidKplusPart->Fill(kplus) ;
   }
   float proton  = fFlowTrack->ProtonPbarProb() ;
   fHistPidProton->Fill(proton) ;
   if(strcmp(pid, "pr+") == 0) 
   {
    fPidId = 4 ; protonN++ ; fHistPidPt->Fill(9.5,pt) ;  
    fHistMeanTPCProton->Fill(lpTPC, dEdx) ;
    fHistMeanITSProton->Fill(lpITS, its);
    fHistMeanTRDProton->Fill(lpTRD, trd);
    fHistMeanTOFProton->Fill(invMass, tof);
    fHistPidProtonPart->Fill(proton) ;
   }
   float deuteron  = fFlowTrack->DeuteriumAntiDeuteriumProb() ;
   fHistPidDeuteron->Fill(deuteron) ;
   if(strcmp(pid, "d+") == 0) 
   {
    fPidId = 6 ; deuteronN++ ; fHistPidPt->Fill(11.5,pt) ; 
    fHistMeanTPCDeuteron->Fill(lpTPC, dEdx) ;
    fHistMeanITSDeuteron->Fill(lpITS, its);
    fHistMeanTRDDeuteron->Fill(lpTRD, trd);
    fHistMeanTOFDeuteron->Fill(invMass, tof);
    fHistPidDeuteronPart->Fill(deuteron) ;
   }
  }
  else if(charge == -1) 
  {
   hMinusN++ ;
   fHistMeanDedxNeg2D->Fill(lpTPC, dEdx) ;
   fHistMeanDedxNeg2DITS->Fill(lpITS, its) ;
   fHistMeanDedxNeg2DTRD->Fill(lpTRD, trd) ;
   fHistMeanDedxNeg2DTOF->Fill(lpTOF, tof) ;
   //
   float electron  = fFlowTrack->ElectronPositronProb() ;
   fHistPidElectron->Fill(electron);
   if(strcmp(pid, "e-") == 0) 
   {
    fPidId = 0 ; electronN++ ; fHistPidPt->Fill(0.5,pt) ; 
    fHistMeanTPCElectron->Fill(lpTPC, dEdx);
    fHistMeanITSElectron->Fill(lpITS, its);
    fHistMeanTRDElectron->Fill(lpTRD, trd);
    fHistMeanTOFElectron->Fill(invMass, tof);
    fHistPidElectronPart->Fill(electron);
   }
   float muonMinus  = fFlowTrack->MuonPlusMinusProb() ;
   fHistPidMuonMinus->Fill(muonMinus) ;
   if(strcmp(pid, "mu-") == 0) 
   {
    fPidId = 1 ; muonMinusN++ ; fHistPidPt->Fill(2.5,pt) ;
    fHistMeanTPCMuonMinus->Fill(lpTPC, dEdx) ;
    fHistMeanITSMuonMinus->Fill(lpITS, its);
    fHistMeanTRDMuonMinus->Fill(lpTRD, trd);
    fHistMeanTOFMuonMinus->Fill(invMass, tof);
    fHistPidMuonMinusPart->Fill(muonMinus) ;
   }
   float piMinus = fFlowTrack->PionPlusMinusProb() ;
   fHistPidPiMinus->Fill(piMinus) ;
   if(strcmp(pid, "pi-") == 0) 
   {
    fPidId = 2 ; piMinusN++ ; fHistPidPt->Fill(4.5,pt) ;
    fHistMeanTPCPiMinus->Fill(lpTPC, dEdx);
    fHistMeanITSPiMinus->Fill(lpITS, its);
    fHistMeanTRDPiMinus->Fill(lpTRD, trd);
    fHistMeanTOFPiMinus->Fill(invMass, tof);
    fHistPidPiMinusPart->Fill(piMinus);
   }
   float kminus  = fFlowTrack->KaonPlusMinusProb() ;
   fHistPidKminus->Fill(kminus);
   if(strcmp(pid, "k-") == 0) 
   {
    fPidId = 3 ; kMinusN++ ; fHistPidPt->Fill(6.5,pt) ;
    fHistMeanTPCKminus->Fill(lpTPC, dEdx);
    fHistMeanITSKminus->Fill(lpITS, its);
    fHistMeanTRDKminus->Fill(lpTRD, trd);
    fHistMeanTOFKminus->Fill(invMass, tof);
    fHistPidKminusPart->Fill(kminus);
   }
   float antiproton  = fFlowTrack->ProtonPbarProb() ;
   fHistPidAntiProton->Fill(antiproton);
   if(strcmp(pid, "pr-") == 0) 
   {
    fPidId = 4 ; pbarN++ ; fHistPidPt->Fill(8.5,pt) ;
    fHistMeanTPCPbar->Fill(lpTPC, dEdx);
    fHistMeanITSPbar->Fill(lpITS, its);
    fHistMeanTRDPbar->Fill(lpTRD, trd);
    fHistMeanTOFPbar->Fill(invMass, tof);
    fHistPidAntiProtonPart->Fill(antiproton);
   }
   float antideuteron  = fFlowTrack->DeuteriumAntiDeuteriumProb() ;
   fHistPidAntiDeuteron->Fill(antideuteron);
   if(strcmp(pid, "d-") == 0) 
   {
    fPidId = 6 ; dbarN++ ; fHistPidPt->Fill(10.5,pt) ;
    fHistMeanTPCAntiDeuteron->Fill(lpTPC, dEdx);
    fHistMeanITSAntiDeuteron->Fill(lpITS, its);
    fHistMeanTRDAntiDeuteron->Fill(lpTRD, trd);
    fHistMeanTOFAntiDeuteron->Fill(invMass, tof);
    fHistPidAntiDeuteronPart->Fill(antideuteron);
   }
  }
  
  // Yield3D, Yield2D, Eta, Pt, Phi, bayP.Id.
  fHistPtot->Fill(totalp) ;
  fHistPt->Fill(pt) ;
  fHistPhi->Fill(phi);  
  fHistAllEtaPtPhi3D->Fill(eta, pt, phi) ;
  fHistYieldAll2D->Fill(eta, pt) ;
  fHistBayPidMult->Fill(fPidId) ;
  if(constrainable) 
  { 
   fHistPhiCons->Fill(phi);  
   fHistPhiPtCon->Fill(phi, pt);  
   fHistYieldCon2D->Fill(eta, pt) ;
   fHistConsEtaPtPhi3D->Fill(eta, pt, phi) ;
   fHistGlobEtaPtPhi3D->Fill(etaGlob, ptGlob, phiGlob) ;
  }
  else 
  { 
   fHistYieldUnc2D->Fill(etaGlob, ptGlob) ;
   fHistUncEtaPtPhi3D->Fill(etaGlob, ptGlob, phiGlob) ;
   fHistPhiPtUnc->Fill(phiGlob, ptGlob) ; 
  }
  if(fFlowTrack->Charge()>0)       { fHistPtPhiPos->Fill(phi, pt); }
  else if(fFlowTrack->Charge()<0)  { fHistPtPhiNeg->Fill(phi, pt); }

  // fills selected part histograms
  if(fFlowSelect->SelectPart(fFlowTrack)) 
  {
   fSelParts++ ;
  
   if(strlen(fFlowSelect->PidPart()) != 0) 
   {
    float rapidity = fFlowTrack->Y();
    fHistBinEta->Fill(rapidity, rapidity);
    fHistYieldPart2D->Fill(rapidity, pt);
   }
   else 
   {
    fHistBinEta->Fill(eta, eta) ;
    fHistYieldPart2D->Fill(eta, pt) ;
   }
   fHistBayPidMultPart->Fill(fPidId) ;
   fHistBinPt->Fill(pt, pt) ;
   fHistEtaPtPhi3DPart->Fill(eta,pt,phi) ;
   fHistDcaGlobalPart->Fill(dcaGlobal) ;
   fHistInvMassPart->Fill(invMass) ;
   if(charge == 1) 
   {
    fHistMeanDedxPos3DPart->Fill(lpITS, dEdx, fPidId) ;
    fHistMeanDedxPos3DPartITS->Fill(lpITS, its, fPidId) ;
   }
   else if(charge == -1) 
   {
    fHistMeanDedxNeg3DPart->Fill(lpITS, dEdx, fPidId) ;
    fHistMeanDedxNeg3DPartITS->Fill(lpITS, its, fPidId) ;
   }
   if(eta > 0.) { etaSymPosTpcNpart++ ; }  // for fHistEtaSymPart 
   else         { etaSymNegTpcNpart++ ; }
  }
  else         
  {
   fHistEtaPtPhi3DOut->Fill(eta,pt,phi) ;
   fHistYieldOut2D->Fill(eta, pt) ;
   fHistDcaGlobalOut->Fill(dcaGlobal) ;
   fHistInvMassOut->Fill(invMass) ;
  }

  // cos(n*phiLab)
  for(int j = 0; j < AliFlowConstants::kHars; j++) 
  {
   bool oddHar = (j+1) % 2 ;
   float order = (float)(j+1) ;
   float vIn   = 100 * cos((double)order * phi) ;
   if(eta < 0 && oddHar) { vIn *= -1 ; }
   fHistCosPhi->Fill(order, vIn);
  }

  //For Eta symmetry TPC
  if(fFlowTrack->FitPtsTPC())
  {
   if(eta > 0.) { etaSymPosTpcN++ ; } // for fHistEtaSym 
   else         { etaSymNegTpcN++ ; }
  }
  
  // not to call it twice ...
  int nEtaS = HarmonicsLoop(fFlowTrack) ;

  //For Correlated Multiplicity 
  corrMultN += (float)nEtaS ;

  //For Correlated Multiplicity in 1 unit rapidity
  if(TMath::Abs(eta) <= 0.5) { corrMultUnit += (float)nEtaS ; }

  //delete pointer
  fFlowTrack = 0 ;
 }      											 // end of tracks loop
						     
 // EtaSym
 float etaSymTpc = 0 ; 
 if(etaSymPosTpcN || etaSymNegTpcN) 	    { etaSymTpc  = (etaSymPosTpcN  - etaSymNegTpcN)  / (etaSymPosTpcN  + etaSymNegTpcN); }
 float etaSymTpcPart = 0 ; 
 if(etaSymPosTpcNpart || etaSymNegTpcNpart) { etaSymTpcPart  = (etaSymPosTpcNpart  - etaSymNegTpcNpart)  / (etaSymPosTpcNpart  + etaSymNegTpcNpart) ; }
 Float_t vertexZ = fVertex[2] ;
 // all
 fHistEtaSym->Fill(etaSymTpc);
 fHistEtaSymVerZ2D->Fill(vertexZ , etaSymTpc);
 // selected
 fHistEtaSymPart->Fill(etaSymTpc);
 fHistEtaSymVerZ2DPart->Fill(vertexZ , etaSymTpcPart);

 // PID multiplicities
 float totalMult = (float)fFlowTracks->GetEntries() ;
 fHistPidMult->Fill(1., totalMult);
 fHistPidMult->Fill(2., hPlusN);
 fHistPidMult->Fill(3., hMinusN);
 fHistPidMult->Fill(4., piPlusN);
 fHistPidMult->Fill(5., piMinusN);
 fHistPidMult->Fill(6., protonN);
 fHistPidMult->Fill(7., pbarN);
 fHistPidMult->Fill(8., kPlusN);
 fHistPidMult->Fill(9., kMinusN);
 fHistPidMult->Fill(10., deuteronN);
 fHistPidMult->Fill(11., dbarN);
 fHistPidMult->Fill(12., electronN);
 fHistPidMult->Fill(13., positronN);
 fHistPidMult->Fill(14., muonMinusN);
 fHistPidMult->Fill(15., muonPlusN);

 // Multiplicity of particles correlated with the event planes
 corrMultN /= (float)(AliFlowConstants::kHars * AliFlowConstants::kSels) ; 
 if(fSelParts == corrMultN) { fHistMultPart->Fill(corrMultN) ; } // crosscheck
 else  			    { fHistMultPart->Fill(-1) ; }
 // ...in one unit rapidity
 corrMultUnit /= (float)(AliFlowConstants::kHars * AliFlowConstants::kSels) ; 
 fHistMultPartUnit->Fill(corrMultUnit) ;

 fTotalNumberOfTracks += fNumberOfTracks ;

 return ;
}
//-----------------------------------------------------------------------
Int_t AliFlowAnalyser::HarmonicsLoop(AliFlowTrack* fFlowTrack)
{
 // HarmonicsLoop, does the correlation analysis, fills histograms
 // (internally called by ::FillParticleHistograms(..) loop) .
 
 float eta = fFlowTrack->Eta() ; 			
 float phi = fFlowTrack->Phi() ; 			
 float pt = fFlowTrack->Pt() ; 			
 float zFirstPoint = fFlowTrack->ZFirstPoint() ; 
 // float zLastPoint = fFlowTrack->ZLastPoint() ;

 // Looping over Selections and Harmonics
 int corrMultN = 0 ;
 for (int k = 0; k < AliFlowConstants::kSels; k++) 
 {
  fFlowSelect->SetSelection(k) ;
  for (int j = 0; j < AliFlowConstants::kHars; j++) 
  {
   bool oddHar = (j+1) % 2;
   fFlowSelect->SetHarmonic(j);
   double order  = (double)(j+1);
   float psii = 0. ; float psi2 = 0. ;
   if(fEtaSub)	       // particles with the opposite subevent
   {
    if(eta > 0) { psii = fPsiSub[1][k][j] ; }     //check
    else	{ psii = fPsiSub[0][k][j] ; }
   } 
   else if(order > 3. && !oddHar) 
   {
    psii = fPsi[k][1];  // 2nd harmomic event plane
    if(psii > 2*TMath::Pi()/order) { psii -= 2*TMath::Pi()/order ; } 
    if(psii > 2*TMath::Pi()/order) { psii -= 2*TMath::Pi()/order ; }
   } 
   else  // random subevents
   {
    psii = fPsi[k][j] ;
   }

   if(fFlowSelect->Select(fFlowTrack)) // Get detID
   {
    Bool_t kTpcPlus  = kFALSE ;
    Bool_t kTpcMinus = kFALSE ;
    Bool_t kTpcAll   = kFALSE ;

    fHistFull[k].fHistFullHar[j].fHistYieldPt->Fill(pt);
    fHistFull[k].fHistFullHar[j].fHistEtaPtPhi3D->Fill(eta, pt, phi);
    fHistFull[k].fHistFullHar[j].fHistYield2D->Fill(eta, pt);
    fHistFull[k].fHistFullHar[j].fHistDcaGlob->Fill(TMath::Abs(fFlowTrack->Dca()));

    // Set Tpc (+ and -)
    if(fFlowTrack->FitPtsTPC())    //OR*: AliESDtrack:: "TBits& GetTPCClusterMap()" or "Int_t GetTPCclusters(Int_t* idx)" ...
    {
     if(zFirstPoint >= 0. && eta > 0.)	    { kTpcPlus  = kTRUE ; } 
     else if(zFirstPoint <= 0. && eta < 0.) { kTpcMinus = kTRUE ; }
     else				    { kTpcAll   = kTRUE ; }	  
    }

    // PID Multiplicities (particle for R.P.) - done just one time for each selection
    if(j==0) { fHistFull[k].fHistBayPidMult->Fill(fPidId) ; }

    // Calculate weights for filling histograms   
    float wt = 1. ; // TMath::Abs(fFlowEvent->Weight(k, j, fFlowTrack)) ; 

    // Fill histograms with selections
    if(kTpcPlus)       { fHistFull[k].fHistFullHar[j].fHistPhiPlus->Fill(phi,wt) ;  } 
    else if(kTpcMinus) { fHistFull[k].fHistFullHar[j].fHistPhiMinus->Fill(phi,wt) ; } 
    else if(kTpcAll)   { fHistFull[k].fHistFullHar[j].fHistPhiAll->Fill(phi,wt) ;   } 
    fHistFull[k].fHistFullHar[j].fHistPhi->Fill(phi,wt) ;

    // Get phiWgt from file
    double phiWgt;
    if(order > 3. && !oddHar) { phiWgt = fFlowEvent->PhiWeight(k, 1, fFlowTrack) ; } 
    else		      { phiWgt = fFlowEvent->PhiWeight(k, j, fFlowTrack) ; }
    if(oddHar && eta<0.)      { phiWgt /= -1. ; } // only for flat hists
    // cout << " Weight [" << k << "][" << j << "] (" << fFlowTrack->GetName() << ") = " << phiWgt << " or " << wt << " for phi = " << phi << endl ; // chk

    // Fill Flat histograms
    if(kTpcPlus)       { fHistFull[k].fHistFullHar[j].fHistPhiFlatPlus->Fill(phi, phiWgt) ; } 
    else if(kTpcMinus) { fHistFull[k].fHistFullHar[j].fHistPhiFlatMinus->Fill(phi, phiWgt) ; } 
    else if(kTpcAll)   { fHistFull[k].fHistFullHar[j].fHistPhiFlatAll->Fill(phi, phiWgt) ; } 
    fHistFull[k].fHistFullHar[j].fHistPhiFlat->Fill(phi,phiWgt) ;
    
    if(oddHar && eta<0.)  { phiWgt *= -1. ; } // restore value

    // Remove autocorrelations
    TVector2 qi;
    if(!fEtaSub)   // random subevents (or other)
    {
     if(order > 3. && !oddHar) // 2nd harmonic event plane
     { 
      qi.Set(phiWgt * cos(phi * 2), phiWgt * sin(phi * 2));
      TVector2 mQi = fQ[k][1] - qi;
      psii = mQi.Phi() / 2;
      if(psii < 0.) { psii += TMath::Pi() ; } 
     } 
     else 
     {
      qi.Set(phiWgt * cos(phi * order), phiWgt * sin(phi * order));
      TVector2 mQi = fQ[k][j] - qi;
      psii = mQi.Phi() / order;
      if(psii < 0.) { psii += 2*TMath::Pi()/order ; }
     }
    }
          
    // Remove autocorrelations of the second order 'particles' which are used for v1{EP1,EP2}.
    if (fV1Ep1Ep2 == kTRUE && order == 1) 
    {
     AliFlowSelection usedForPsi2 = *fFlowSelect ;
     usedForPsi2.SetHarmonic(1);
     if(usedForPsi2.Select(fFlowTrack))  // particle was used for Psi2
     {
      qi.Set(phiWgt * cos(phi * 2), phiWgt * sin(phi * 2));
      TVector2 mQi = fQ[k][1] - qi;
      psi2 = mQi.Phi() / 2;
      if(psi2 < 0.) { psi2 += TMath::Pi() ; }
     }
     else				 // particle was not used for Psi2
     { 
      psi2 = fPsi[k][1];
     }
    }
   }
   else
   {
    fHistFull[k].fHistFullHar[j].fHistYieldPtout->Fill(pt);
    fHistFull[k].fHistFullHar[j].fHistEtaPtPhi3Dout->Fill(eta, pt, phi);
    fHistFull[k].fHistFullHar[j].fHistYield2Dout->Fill(eta, pt);
    fHistFull[k].fHistFullHar[j].fHistDcaGlobout->Fill(TMath::Abs(fFlowTrack->Dca()));
   }

   // Caculate v for all particles selected for correlation analysis
   if(fFlowSelect->SelectPart(fFlowTrack)) 
   {
    corrMultN++;
    
    float v;
    if (fV1Ep1Ep2 == kFALSE || order != 1) 
    {
     v = 100 * cos(order * (phi - psii)) ;
    }
    else // i.e. (fV1Ep1Ep2 == kTRUE && order == 1)
    {
     v = 100 * cos(phi + psii - 2*psi2) ;
    }
    
    float vFlip = v;
    if(eta < 0 && oddHar) { vFlip *= -1 ; }
    if(strlen(fFlowSelect->PidPart()) != 0) // pid, fill rapidity 
    {
     float rapidity = fFlowTrack->Y();
     fHistFull[k].fHistFullHar[j].fHistvObs2D->Fill(rapidity, pt, v);
   
     if(fPtRangevEta[1] > fPtRangevEta[0])  // cut is used
     {
      if(pt < fPtRangevEta[1] && pt >= fPtRangevEta[0]) // check cut range, fill if in range
      {
       fHistFull[k].fHistFullHar[j].fHistvObsEta->Fill(rapidity, v);
      }
     }
     else // cut is not used, fill in any case
     { 
      fHistFull[k].fHistFullHar[j].fHistvObsEta->Fill(rapidity, v);
     }
    } 
    else  // no pid, fill eta
    {
     fHistFull[k].fHistFullHar[j].fHistvObs2D->Fill(eta, pt, v);

     if(fPtRangevEta[1] > fPtRangevEta[0]) // cut is used
     {
      if(pt < fPtRangevEta[1] && pt >= fPtRangevEta[0]) // check cut range, fill if in range
      {
       fHistFull[k].fHistFullHar[j].fHistvObsEta->Fill(eta, v);
      }
     }
     else // cut is not used, fill in any case
     { 
      fHistFull[k].fHistFullHar[j].fHistvObsEta->Fill(eta, v);
     }
    }

    if(fEtaRangevPt[1] > fEtaRangevPt[0]) // cut is used
    { 
     if(TMath::Abs(eta) < fEtaRangevPt[1] && TMath::Abs(eta) >= fEtaRangevPt[0]) // check cut range, fill if in range
     {
      fHistFull[k].fHistFullHar[j].fHistvObsPt->Fill(pt, vFlip);  // for odd harmonis /-/
     }
    }
    else  // cut is not used, fill in any case
    {
     fHistFull[k].fHistFullHar[j].fHistvObsPt->Fill(pt, vFlip);
    }

    // v_
    Bool_t etaPtNoCut = kTRUE;
    if(fPtRangevEta[1] > fPtRangevEta[0] && (pt < fPtRangevEta[0] || pt >= fPtRangevEta[1])) 
    {
     etaPtNoCut = kFALSE;
    }
    if(fEtaRangevPt[1] > fEtaRangevPt[0] && (TMath::Abs(eta) < fEtaRangevPt[0] || TMath::Abs(eta) >= fEtaRangevPt[1]))
    {
     etaPtNoCut = kFALSE;
    }
    if(etaPtNoCut) { fHistFull[k].fHistvObs->Fill(order, vFlip) ; }
 
    // Correlation of Phi of selected particles with Psi
    float phii = phi;
    if(eta < 0 && oddHar) 
    {
     phii += TMath::Pi() ; // backward particle and odd harmonic
     if(phii > 2*TMath::Pi()) { phii -= 2*TMath::Pi() ; }
    }
    float dPhi = phii - psii;
    if(dPhi < 0.)	       { dPhi += 2*TMath::Pi() ; }
    fHistFull[k].fHistFullHar[j].fHistPhiCorr->Fill(fmod((double)dPhi, 2*TMath::Pi() / order));
   }
  }
 }
 
 return corrMultN ;
}
//-----------------------------------------------------------------------
void AliFlowAnalyser::FillV0Histograms(TClonesArray* fFlowV0s)
{
 // v0s histograms

 cout << " V0s Loop . " << endl ;

 fSelV0s = 0 ;
 for(fV0Number=0;fV0Number<fNumberOfV0s;fV0Number++) 
 {
  //fFlowV0 = (AliFlowV0*)fFlowV0s->At(fV0Number) ;	  
  
  fFlowV0 = (AliFlowV0*)(fFlowV0s->AddrAt(fV0Number)) ;
  // cout << " looping v0 n. " << fV0Number << " * " << fFlowV0 << endl ;  

  float mass = fFlowV0->Mass() ;   	   	    
  float eta = fFlowV0->Eta() ;   	   	    
  float rapidity = fFlowV0->Y() ;   	   	    
  float phi = fFlowV0->Phi() ;     	   	    
  float pt = fFlowV0->Pt() ;      	   	    
  float totalp = fFlowV0->P() ;  	   	    
  //int	charge = fFlowV0->Charge() ;	    	    
  float dca = fFlowV0->Dca() ;		   	    
  float lenght = fFlowV0->V0Lenght() ;     	    
  float sigma = fFlowV0->Sigma() ; 	    	    
  float chi2 = fFlowV0->Chi2() ; 	
  Char_t pid[10] ; strcpy(pid, fFlowV0->Pid()) ;    
  int labelPlus =  -1 ; 
  int labelMinus = -1 ;   	    
  AliFlowTrack* daughterPlus = (AliFlowTrack*)(fFlowTracks->AddrAt(fFlowV0->DaughterP())) ;
  AliFlowTrack* daughterMinus = (AliFlowTrack*)(fFlowTracks->AddrAt(fFlowV0->DaughterN())) ; ; 
  if(daughterPlus)  { labelPlus = daughterPlus->Label() ; }		    
  if(daughterMinus) { labelMinus = daughterMinus->Label() ; }			    

  fHistV0Mass->Fill(mass) ;
  fHistV0EtaPtPhi3D->Fill(eta, pt, phi) ;
  fHistV0YieldAll2D->Fill(eta, pt) ;
  fHistV0Dca->Fill(dca);
  fHistV0Chi2->Fill(chi2);
  fHistV0Lenght->Fill(lenght);
  fHistV0Sigma->Fill(sigma);
  fHistV0MassPtSlices->Fill(mass,pt);
  fHistV0PYall2D->Fill(rapidity, totalp) ;

  if(fFlowSelect->SelectPart(fFlowV0))
  {
   bool inWin = fFlowSelect->SelectV0Part(fFlowV0) ;
   bool sx = fFlowSelect->SelectV0sxSide(fFlowV0) ;
   bool dx = fFlowSelect->SelectV0dxSide(fFlowV0) ;
   fSelV0s++ ;

   fHistV0YieldPart2D->Fill(eta, pt) ;
   if(inWin) 
   { 
    fHistV0EtaPtPhi3DPart->Fill(eta, pt, phi) ;
    fHistV0LenghtPart->Fill(lenght);
    fHistV0DcaPart->Fill(dca);
    fHistV0MassWin->Fill(mass) ; 
    fHistV0BinEta->Fill(eta, eta);
    fHistV0BinPt->Fill(pt, pt);   
   }
   else      
   { 
    fHistV0sbEtaPtPhi3DPart->Fill(eta, pt, phi) ;
    fHistV0sbLenghtPart->Fill(lenght);
    fHistV0sbDcaPart->Fill(dca);
    fHistV0sbMassSide->Fill(mass) ; 
    fHistV0sbBinEta->Fill(eta, eta);
    fHistV0sbBinPt->Fill(pt, pt);
   }

   for(int k = 0; k < AliFlowConstants::kSels; k++)  // sort of HarmonicsLoop - selection number used
   {
    fFlowSelect->SetSelection(k) ;
    for(Int_t j=0;j<AliFlowConstants::kHars;j++) 
    {
     Bool_t oddHar = (j+1) % 2 ;
     Float_t order = (Float_t)(j+1) ;
     fFlowSelect->SetHarmonic(j);	

     // Remove autocorrelations
     Float_t psii, psi2 ;
     TVector2 q1, q2 ;
     Float_t phiDaughter1 = 0. ; Float_t phiDaughter2 = 0. ;
     Double_t phiWgt1 = 0. ; Double_t phiWgt2 = 0. ;
     // -
     if(daughterPlus)
     { // 1
      fFlowTrack = daughterPlus ;
      phiDaughter1 = fFlowTrack->Phi() ;		  
      if(fFlowSelect->Select(fFlowTrack)) // Get phiWgt from file
      {
       if(order > 3. && !oddHar) { phiWgt1 = fFlowEvent->PhiWeight(k, 1, fFlowTrack) ; } 
       else			 { phiWgt1 = fFlowEvent->PhiWeight(k, j, fFlowTrack) ; }
      }
     }
     if(daughterMinus)
     { // 2
      fFlowTrack = daughterMinus ;
      phiDaughter2 = fFlowTrack->Phi() ;
      if(fFlowSelect->Select(fFlowTrack)) // Get phiWgt from file
      {
       if(order > 3. && !oddHar) { phiWgt2 = fFlowEvent->PhiWeight(k, 1, fFlowTrack) ; } 
       else			 { phiWgt2 = fFlowEvent->PhiWeight(k, j, fFlowTrack) ; }
      } 
     }

     // psi2
     q1.Set(phiWgt1 * cos(phiDaughter1 * 2), phiWgt1 * sin(phiDaughter1 * 2));
     q2.Set(phiWgt2 * cos(phiDaughter2 * 2), phiWgt2 * sin(phiDaughter2 * 2));
     TVector2 mQi = fQ[k][1] ; mQi -= q1 ; mQi -= q2 ; 
     psi2 = mQi.Phi() / 2 ;	  
     if(psi2 < 0.)	     { psi2 += TMath::Pi() ; }
     if(psi2 > TMath::Pi())  { psi2 -= TMath::Pi() ; } 

     // psii
     if(order > 3. && !oddHar) { psii = psi2 ; } 
     else 
     {
      q1.Set(phiWgt1 * cos(phiDaughter1 * order), phiWgt1 * sin(phiDaughter1 * order));
      q2.Set(phiWgt2 * cos(phiDaughter2 * order), phiWgt2 * sin(phiDaughter2 * order));
      TVector2 mQi = fQ[k][j] ; mQi -= q1 ; mQi -= q2 ;
      psii = mQi.Phi()/order ; 
      if(psii < 0.)		     { psii += 2*TMath::Pi()/order ; }	
      if(psii > 2*TMath::Pi()/order) { psii -= 2*TMath::Pi()/order ; } 
     }

    // Caculate v for all V0s selected for correlation analysis
     float v ;
     if(fV1Ep1Ep2 == kFALSE || order != 1) { v = 100 * cos(order * (phi - psii)) ; }
     else { v = 100 * cos(phi + psii - 2*psi2) ; }
     float vFlip = v ; if(eta < 0 && oddHar) { vFlip *= -1 ; }

     // invariant mass windows & sidebands
     if(inWin) { fHistFull[k].fHistFullHar[j].fHistV0vObs2D->Fill(eta, pt, v) ; }
     else      { fHistFull[k].fHistFullHar[j].fHistV0sbvObs2D->Fill(eta, pt, v) ; }

     if(fPtRangevEta[1] > fPtRangevEta[0]) // cut is used
     {
      if(pt < fPtRangevEta[1] && pt >= fPtRangevEta[0]) // check cut range, fill if in range
      {
       if(inWin) { fHistFull[k].fHistFullHar[j].fHistV0vObsEta->Fill(eta, v) ; }
       else      
       { 
        fHistFull[k].fHistFullHar[j].fHistV0sbvObsEta->Fill(eta, v) ; 
	if(sx)      { fHistFull[k].fHistFullHar[j].fHistV0sbvObsEtaSx->Fill(eta, v) ; } 
	else if(dx) { fHistFull[k].fHistFullHar[j].fHistV0sbvObsEtaDx->Fill(eta, v) ; } 
       }
      }
     }
     else // cut is not used, fill in any case
     { 
      if(inWin) { fHistFull[k].fHistFullHar[j].fHistV0vObsEta->Fill(eta, v) ; }
      else	
      { 
       fHistFull[k].fHistFullHar[j].fHistV0sbvObsEta->Fill(eta, v) ; 
       if(sx)	   { fHistFull[k].fHistFullHar[j].fHistV0sbvObsEtaSx->Fill(eta, v) ; } 
       else if(dx) { fHistFull[k].fHistFullHar[j].fHistV0sbvObsEtaDx->Fill(eta, v) ; } 
      }
     } 
     if(fEtaRangevPt[1] > fEtaRangevPt[0]) // cut is used
     {
      if(TMath::Abs(eta) < fEtaRangevPt[1] && TMath::Abs(eta) >= fEtaRangevPt[0]) // check cut range, fill if in range
      {          
       if(inWin) { fHistFull[k].fHistFullHar[j].fHistV0vObsPt->Fill(pt, vFlip) ; }      // for odd harmonis /-/
       else      
       { 
        fHistFull[k].fHistFullHar[j].fHistV0sbvObsPt->Fill(pt, vFlip) ; 
        if(sx)	    { fHistFull[k].fHistFullHar[j].fHistV0sbvObsPtSx->Fill(eta, v) ; } 
        else if(dx) { fHistFull[k].fHistFullHar[j].fHistV0sbvObsPtDx->Fill(eta, v) ; } 	
       }
      }
     }
     else  // cut is not used, fill in any case
     {
      if(inWin) { fHistFull[k].fHistFullHar[j].fHistV0vObsPt->Fill(pt, vFlip) ; } 
      else	
      { 
       fHistFull[k].fHistFullHar[j].fHistV0sbvObsPt->Fill(pt, vFlip) ; 
       if(sx)	   { fHistFull[k].fHistFullHar[j].fHistV0sbvObsPtSx->Fill(eta, v) ; } 
       else if(dx) { fHistFull[k].fHistFullHar[j].fHistV0sbvObsPtDx->Fill(eta, v) ; }  
      }
     }
     // v_
     Bool_t etaPtNoCut = kTRUE;
     if(fPtRangevEta[1] > fPtRangevEta[0] && (pt < fPtRangevEta[0] || pt >= fPtRangevEta[1])) 
     {
      etaPtNoCut = kFALSE;
     }
     if(fEtaRangevPt[1] > fEtaRangevPt[0] && (TMath::Abs(eta) < fEtaRangevPt[0] || TMath::Abs(eta) >= fEtaRangevPt[1]))
     {
      etaPtNoCut = kFALSE;
     }
     if(etaPtNoCut) 
     { 
      if(inWin) { fHistFull[k].fHistV0vObs->Fill(order, vFlip) ; }
      else	
      { 
       if(sx)	   { fHistFull[k].fHistV0sbvObsSx->Fill(order, vFlip) ; } 
       else if(dx) { fHistFull[k].fHistV0sbvObsDx->Fill(order, vFlip) ; }  
      }
     }
   
     // Correlation of Phi of selected v0s with Psi
     float phii = phi;
     if(eta < 0 && oddHar)
     {
      phii += TMath::Pi() ; // backward particle and odd harmonic
      if(phii > 2*TMath::Pi()) { phii -= 2*TMath::Pi() ; }
     }
     float dPhi = phii - psii;
     if(dPhi < 0.)  { dPhi += 2*TMath::Pi() ; }

     if(inWin) { fHistFull[k].fHistFullHar[j].fHistV0PhiCorr->Fill(fmod((double)dPhi, 2*TMath::Pi() / order)) ; } 
     else      { fHistFull[k].fHistFullHar[j].fHistV0sbPhiCorr->Fill(fmod((double)dPhi, 2*TMath::Pi() / order)) ; }
    }
   }
  }
  //delete fFlowV0 ; fFlowV0 = 0 ;
  //delete daughterPlus ; daughterPlus = 0 ;  
  //delete daughterMinus ; daughterMinus = 0 ;  
 }
 fHistV0MultPart->Fill(fSelV0s) ;
 fTotalNumberOfV0s += fNumberOfV0s ;

 return ;
}
//-----------------------------------------------------------------------
// void AliFlowAnalyser::FillLabels()  
// {
//  // Reads the tracks label (i.e. the link from ESD tracking to KineTree)
//  // and ...
//  
//  cout << " Fill Labels . " << endl ;  
// 
//  fLabHist ;
//  return ;
// }
//-----------------------------------------------------------------------
void AliFlowAnalyser::PrintEventQuantities() const
{
 // prints event by event calculated quantities

 cout << endl ; cout << " # " << fEventNumber << " - Event quantities : " << endl ; cout << endl ; 
 
 cout << "fQ[k][j] = " ; 
 for(int k=0;k<AliFlowConstants::kSels;k++) 
 {
  for(int j=0;j<AliFlowConstants::kHars;j++) 
  { 
   cout << (Float_t)fQ[k][j].X() << "," << (Float_t)fQ[k][j].Y() << " ; " ; 
  }
  cout << endl ;
  cout << "	       " ; 
 }
 cout << endl ; cout << "fPsi[k][j] = " ; 
 for(int k=0;k<AliFlowConstants::kSels;k++) 
 {
  for(int j=0;j<AliFlowConstants::kHars;j++) { Float_t aaa = (Float_t)fPsi[k][j] ; cout << aaa << " , " ; }
  cout << endl ;
  cout << "		 " ; 
 }
 cout << endl ; cout << "fQnorm[k][j] = " ; 
 for(int k=0;k<AliFlowConstants::kSels;k++) 
 {
  for(int j=0;j<AliFlowConstants::kHars;j++) { Float_t aaa = (Float_t)fQnorm[k][j] ; cout << aaa << " , " ; }
  cout << endl ;
  cout << "		" ; 
 }
 cout << endl ; cout << "fMult[k][j] = " ; 
 for(int k=0;k<AliFlowConstants::kSels;k++) 
 {
  for(int j=0;j<AliFlowConstants::kHars;j++) { Float_t aaa = (Float_t)fMult[k][j] ; cout << aaa << " , " ; } 
  cout << endl ;
  cout << "		  " ; 
 }
 cout << endl ; cout << "fMultSub[n][k][j] = " ; 
 for(int n=0;n<AliFlowConstants::kSubs;n++) 
 {
  for(int k=0;k<AliFlowConstants::kSels;k++) 
  {
   for(int j=0;j<AliFlowConstants::kHars;j++) { Float_t aaa = fMultSub[n][k][j] ; cout << aaa << " , " ; } 
   cout << endl ;
   cout << "		     " ; 
  }
 }
 cout << endl ; cout << "fPsiSub[n][k][j] = " ; 
 for(int n=0;n<AliFlowConstants::kSubs;n++) 
 {
  for(int k=0;k<AliFlowConstants::kSels;k++) 
  {
   for(int j=0;j<AliFlowConstants::kHars;j++) { Float_t aaa = fPsiSub[n][k][j] ; cout << aaa << " , " ; } 
   cout << endl ;
   cout << "		     " ; 
  }
 }
 cout << endl ; cout << "Delta_PsiSub[k][j] = " ; 
 for(int k=0;k<AliFlowConstants::kSels;k++) 
 {
  for(int j=0;j<AliFlowConstants::kHars;j++) { Float_t aaa = fPsiSub[0][k][j]-fPsiSub[1][k][j] ; cout << aaa << " , " ; } 
  cout << endl ;
  cout << "		" ; 
 }
 cout << endl ;

 cout << " N. of selected tracks - SelPart(AliFlowTrack*) = " << fSelParts << endl ;	  		    		//! n. of tracks selected for correlation analysis
 cout << " N. of selected V0s    - SelPart(AliFlowV0*)    = " << fSelV0s << endl ;	  		    		//! n. of v0s selected for correlation analysis

 cout << endl ;
}
//----------------------------------------------------------------------
// ###
//----------------------------------------------------------------------
Bool_t AliFlowAnalyser::Resolution()
{
 // Calculates the resolution from cos(dPsi) and corrects the observed flow.
 // New histograms are then created and filled with the corrected vn 
 // (see fVnResHistList), and saved in the otput file.
 
 cout << " AliFlowAnalyser::Resolution()" << endl ;
 
 // VnRes histogram collection
 fVnResHistList = new TOrdCollection(AliFlowConstants::kSels*AliFlowConstants::kHars);
 
 // Calculate resolution from sqrt(fHistCos)
 double cosPair[AliFlowConstants::kSels][AliFlowConstants::kHars];
 double cosPairErr[AliFlowConstants::kSels][AliFlowConstants::kHars];
 double content, contentV0;
 double error, errorV0;
 double totalError;
 TString* histTitle;

 for(int k = 0; k < AliFlowConstants::kSels; k++)
 {
  // v for tracks (corrected for resolution)
  histTitle = new TString("Flow_v_Sel");
  *histTitle += k+1;
  fHistFull[k].fHistv = fHistFull[k].fHistvObs->ProjectionX(histTitle->Data());
  fHistFull[k].fHistv->SetTitle(histTitle->Data());
  fHistFull[k].fHistv->SetXTitle("Harmonic");
  fHistFull[k].fHistv->SetYTitle("v (%)");
  delete histTitle;
  fVnResHistList->AddLast(fHistFull[k].fHistv);
   
  // v for v0s (corrected for resolution)
  histTitle = new TString("FlowV0_v_Sel");
  *histTitle += k+1;
  fHistFull[k].fHistV0v = fHistFull[k].fHistV0vObs->ProjectionX(histTitle->Data());
  fHistFull[k].fHistV0v->SetTitle(histTitle->Data());
  fHistFull[k].fHistV0v->SetXTitle("Harmonic");
  fHistFull[k].fHistV0v->SetYTitle("v (%)");
  delete histTitle;
  fVnResHistList->AddLast(fHistFull[k].fHistV0v);
   
  for(int j = 0; j < AliFlowConstants::kHars; j++) 
  {
   double order = (double)(j+1);
   cosPair[k][j]    = fHistFull[k].fHistCos->GetBinContent(j+1);
   cosPairErr[k][j] = fHistFull[k].fHistCos->GetBinError(j+1);
   if(cosPair[k][j] > 0.) 
   {
    double resSub = 0. ;
    double resSubErr = 0. ;
    if(fV1Ep1Ep2 == kTRUE && order == 1) // calculate resolution of second order event plane first
    {
     double res2 = 0. ;
     double res2error = 0. ;
     if(fHistFull[k].fHistCos->GetBinContent(2) > 0.) 
     {
      if(fEtaSub)  // sub res only
      {
       res2 = TMath::Sqrt(fHistFull[k].fHistCos->GetBinContent(2));
       res2error = fHistFull[k].fHistCos->GetBinError(2) / (2. * res2);
      }
      else 
      {
       if (fHistFull[k].fHistCos->GetBinContent(2) > 0.92)  // resolution saturates
       {
        res2 = 0.99;
        res2error = 0.007;
       } 
       else 
       {
        double deltaRes2Sub = 0.005;  // differential for the error propagation
        double res2Sub = TMath::Sqrt(fHistFull[k].fHistCos->GetBinContent(2));
        double res2SubErr = fHistFull[k].fHistCos->GetBinError(2) / (2. * res2Sub);
        double chiSub2 = Chi(res2Sub);
        double chiSub2Delta = Chi(res2Sub + deltaRes2Sub);
        res2 = ResEventPlane(TMath::Sqrt(2.) * chiSub2); // full event plane res.
        double mRes2Delta = ResEventPlane(TMath::Sqrt(2.) * chiSub2Delta);
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
     fRes[k][j]    = TMath::Sqrt(cosPair[k][0]*res2);
     fResErr[k][j] = 1./(2.*fRes[k][j]) * TMath::Sqrt(cosPairErr[k][0]*cosPairErr[k][0] + res2error*res2error); // Gaussian error propagation
    }
    else if(fEtaSub)  // sub res only
    {
     resSub = TMath::Sqrt(cosPair[k][j]);                
     resSubErr = cosPairErr[k][j] / (2. * resSub);
     fRes[k][j]    = resSub;
     fResErr[k][j] = resSubErr;
    } 
    else if(order==4. || order==6.|| order==8.)  // 2nd harmonic event plane
    {
     double deltaResSub = 0.005;  // differential for the error propagation
     double resSub = TMath::Sqrt(cosPair[k][1]);
     double resSubErr = cosPairErr[k][1] / (2. * resSub);
     double chiSub = Chi(resSub);
     double chiSubDelta = Chi(resSub + deltaResSub);
     double mResDelta;
     if (order==4.) 
     {
      fRes[k][j] = ResEventPlaneK2(TMath::Sqrt(2.) * chiSub); // full event plane res.
      mResDelta = ResEventPlaneK2(TMath::Sqrt(2.) * chiSubDelta);
     } 
     else if(order==6.) 
     {
      fRes[k][j] = ResEventPlaneK3(TMath::Sqrt(2.) * chiSub); // full event plane res.
      mResDelta = ResEventPlaneK3(TMath::Sqrt(2.) * chiSubDelta);
     } 
     else 
     {
      fRes[k][j] = ResEventPlaneK4(TMath::Sqrt(2.) * chiSub); // full event plane res.
      mResDelta = ResEventPlaneK4(TMath::Sqrt(2.) * chiSubDelta);
     }
     fResErr[k][j] = resSubErr * fabs((double)fRes[k][j] - mResDelta) / deltaResSub;
    }
    else 
    {
     if(cosPair[k][j] > 0.92) // resolution saturates
     {
      fRes[k][j]    = 0.99;
      fResErr[k][j] = 0.007;
     } 
     else 
     {
      double deltaResSub = 0.005;  // differential for the error propagation
      double resSub = TMath::Sqrt(cosPair[k][j]);
      double resSubErr = cosPairErr[k][j] / (2. * resSub);
      double chiSub = Chi(resSub);
      double chiSubDelta = Chi(resSub + deltaResSub);
      fRes[k][j] = ResEventPlane(TMath::Sqrt(2.) * chiSub); // full event plane res.
      double mResDelta = ResEventPlane(TMath::Sqrt(2.) * chiSubDelta);
      fResErr[k][j] = resSubErr * TMath::Abs((double)fRes[k][j] - mResDelta) / deltaResSub;
     }
    }
   }
   else  // subevent correlation must be positive
   {
    fRes[k][j]    = 0.;
    fResErr[k][j] = 0.;
   }
   fHistFull[k].fHistRes->SetBinContent(j+1, fRes[k][j]);
   fHistFull[k].fHistRes->SetBinError(j+1, fResErr[k][j]);

   // Create the v 2D histogram (Flow corrected for res.) - v,Pt,Eta
   histTitle = new TString("Flow_v2D_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   fHistFull[k].fHistFullHar[j].fHistv2D = fHistFull[k].fHistFullHar[j].fHistvObs2D->ProjectionXY(histTitle->Data());
   fHistFull[k].fHistFullHar[j].fHistv2D->SetTitle(histTitle->Data());
   fHistFull[k].fHistFullHar[j].fHistv2D->SetXTitle((char*)fLabel.Data());
   fHistFull[k].fHistFullHar[j].fHistv2D->SetYTitle("Pt (GeV/c)");
   fHistFull[k].fHistFullHar[j].fHistv2D->SetZTitle("v (%)");
   delete histTitle;
   fVnResHistList->AddLast(fHistFull[k].fHistFullHar[j].fHistv2D);

   // Create the 1D v histograms - v,Eta
   histTitle = new TString("Flow_vEta_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   fHistFull[k].fHistFullHar[j].fHistvEta = fHistFull[k].fHistFullHar[j].fHistvObsEta->ProjectionX(histTitle->Data());
   fHistFull[k].fHistFullHar[j].fHistvEta->SetTitle(histTitle->Data());
   fHistFull[k].fHistFullHar[j].fHistvEta->SetXTitle((char*)fLabel.Data());
   fHistFull[k].fHistFullHar[j].fHistvEta->SetYTitle("v (%)");
   delete histTitle;
   fVnResHistList->AddLast(fHistFull[k].fHistFullHar[j].fHistvEta);
   
   // v,Pt
   histTitle = new TString("Flow_vPt_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   fHistFull[k].fHistFullHar[j].fHistvPt = fHistFull[k].fHistFullHar[j].fHistvObsPt->ProjectionX(histTitle->Data());
   fHistFull[k].fHistFullHar[j].fHistvPt->SetTitle(histTitle->Data());
   fHistFull[k].fHistFullHar[j].fHistvPt->SetXTitle("Pt (GeV/c)");
   fHistFull[k].fHistFullHar[j].fHistvPt->SetYTitle("v (%)");
   delete histTitle;
   fVnResHistList->AddLast(fHistFull[k].fHistFullHar[j].fHistvPt);

   // Create the v 2D histogram (V0s Flow corrected for res.) - v,Pt,Eta
   histTitle = new TString("FlowV0_v2D_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   fHistFull[k].fHistFullHar[j].fHistV0v2D = fHistFull[k].fHistFullHar[j].fHistV0vObs2D->ProjectionXY(histTitle->Data());
   fHistFull[k].fHistFullHar[j].fHistV0v2D->SetTitle(histTitle->Data());
   fHistFull[k].fHistFullHar[j].fHistV0v2D->SetXTitle((char*)fLabel.Data());
   fHistFull[k].fHistFullHar[j].fHistV0v2D->SetYTitle("Pt (GeV/c)");
   fHistFull[k].fHistFullHar[j].fHistV0v2D->SetZTitle("v (%)");
   delete histTitle;
   fVnResHistList->AddLast(fHistFull[k].fHistFullHar[j].fHistV0v2D);

   // Create the 1D v histograms (V0s) - v,Eta
   histTitle = new TString("FlowV0_vEta_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   fHistFull[k].fHistFullHar[j].fHistV0vEta = fHistFull[k].fHistFullHar[j].fHistV0vObsEta->ProjectionX(histTitle->Data());
   fHistFull[k].fHistFullHar[j].fHistV0vEta->SetTitle(histTitle->Data());
   fHistFull[k].fHistFullHar[j].fHistV0vEta->SetXTitle((char*)fLabel.Data());
   fHistFull[k].fHistFullHar[j].fHistV0vEta->SetYTitle("v (%)");
   delete histTitle;
   fVnResHistList->AddLast(fHistFull[k].fHistFullHar[j].fHistV0vEta);
   
   // (V0s) v,Pt
   histTitle = new TString("FlowV0_vPt_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   fHistFull[k].fHistFullHar[j].fHistV0vPt = fHistFull[k].fHistFullHar[j].fHistV0vObsPt->ProjectionX(histTitle->Data());
   fHistFull[k].fHistFullHar[j].fHistV0vPt->SetTitle(histTitle->Data());
   fHistFull[k].fHistFullHar[j].fHistV0vPt->SetXTitle("Pt (GeV/c)");
   fHistFull[k].fHistFullHar[j].fHistV0vPt->SetYTitle("v (%)");
   delete histTitle;
   fVnResHistList->AddLast(fHistFull[k].fHistFullHar[j].fHistV0vPt);

   // Create the v 2D histogram (V0s sidebands Flow corrected for res.) - v,Pt,Eta
   histTitle = new TString("FlowV0sb_v2D_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   fHistFull[k].fHistFullHar[j].fHistV0sbv2D = fHistFull[k].fHistFullHar[j].fHistV0sbvObs2D->ProjectionXY(histTitle->Data());
   fHistFull[k].fHistFullHar[j].fHistV0sbv2D->SetTitle(histTitle->Data());
   fHistFull[k].fHistFullHar[j].fHistV0sbv2D->SetXTitle((char*)fLabel.Data());
   fHistFull[k].fHistFullHar[j].fHistV0sbv2D->SetYTitle("Pt (GeV/c)");
   fHistFull[k].fHistFullHar[j].fHistV0sbv2D->SetZTitle("v (%)");
   delete histTitle;
   fVnResHistList->AddLast(fHistFull[k].fHistFullHar[j].fHistV0sbv2D);

   // Create the 1D v histograms (V0s sidebands) - v,Eta
   histTitle = new TString("FlowV0sb_vEta_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   fHistFull[k].fHistFullHar[j].fHistV0sbvEta = fHistFull[k].fHistFullHar[j].fHistV0sbvObsEta->ProjectionX(histTitle->Data());
   fHistFull[k].fHistFullHar[j].fHistV0sbvEta->SetTitle(histTitle->Data());
   fHistFull[k].fHistFullHar[j].fHistV0sbvEta->SetXTitle((char*)fLabel.Data());
   fHistFull[k].fHistFullHar[j].fHistV0sbvEta->SetYTitle("v (%)");
   delete histTitle;
   fVnResHistList->AddLast(fHistFull[k].fHistFullHar[j].fHistV0sbvEta);
   
   // (V0s sidebands) v,Pt
   histTitle = new TString("FlowV0sb_vPt_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   fHistFull[k].fHistFullHar[j].fHistV0sbvPt = fHistFull[k].fHistFullHar[j].fHistV0sbvObsPt->ProjectionX(histTitle->Data());
   fHistFull[k].fHistFullHar[j].fHistV0sbvPt->SetTitle(histTitle->Data());
   fHistFull[k].fHistFullHar[j].fHistV0sbvPt->SetXTitle("Pt (GeV/c)");
   fHistFull[k].fHistFullHar[j].fHistV0sbvPt->SetYTitle("v (%)");
   delete histTitle;
   fVnResHistList->AddLast(fHistFull[k].fHistFullHar[j].fHistV0sbvPt);

  // Calulate v = vObs / Resolution
   if(fRes[k][j]) 
   {
    cout << "***** Resolution of the " << j+1 << "th harmonic = " << fRes[k][j] << " +/- " << fResErr[k][j] << endl;
    // selected tracks
    fHistFull[k].fHistFullHar[j].fHistv2D->Scale(1. / fRes[k][j]);
    fHistFull[k].fHistFullHar[j].fHistvEta->Scale(1. / fRes[k][j]);
    fHistFull[k].fHistFullHar[j].fHistvPt->Scale(1. / fRes[k][j]);
    content = fHistFull[k].fHistv->GetBinContent(j+1);        
    content /=  fRes[k][j]; 
    fHistFull[k].fHistv->SetBinContent(j+1, content);            
    // selected V0s    
    fHistFull[k].fHistFullHar[j].fHistV0v2D->Scale(1. / fRes[k][j]);
    fHistFull[k].fHistFullHar[j].fHistV0vEta->Scale(1. / fRes[k][j]);
    fHistFull[k].fHistFullHar[j].fHistV0vPt->Scale(1. / fRes[k][j]);
    contentV0 = fHistFull[k].fHistV0v->GetBinContent(j+1);        
    contentV0 /=  fRes[k][j];
    fHistFull[k].fHistV0v->SetBinContent(j+1, contentV0);            
    // V0s sidebands    
    fHistFull[k].fHistFullHar[j].fHistV0sbv2D->Scale(1. / fRes[k][j]);
    fHistFull[k].fHistFullHar[j].fHistV0sbvEta->Scale(1. / fRes[k][j]);
    fHistFull[k].fHistFullHar[j].fHistV0sbvPt->Scale(1. / fRes[k][j]);
    // The systematic error of the resolution is folded in.
    // tracks
    error = fHistFull[k].fHistv->GetBinError(j+1) ; 
    error /= fRes[k][j];                                         
    totalError = fabs(content) * TMath::Sqrt((error/content)*(error/content) + (fResErr[k][j]/fRes[k][j])*(fResErr[k][j]/fRes[k][j]));
    fHistFull[k].fHistv->SetBinError(j+1, totalError);
    cout << "***** v" << j+1 << "= (" << content << " +/- " << error << " +/- " << totalError << "(with syst.)) %" << endl;
    // V0s
    errorV0 = fHistFull[k].fHistV0v->GetBinError(j+1) ; 
    errorV0 /= fRes[k][j]; 
    if (contentV0>0)
      totalError = fabs(contentV0) * TMath::Sqrt((errorV0/contentV0)*(errorV0/contentV0) + (fResErr[k][j]/fRes[k][j])*(fResErr[k][j]/fRes[k][j]));
    else
      totalError = 0;
    fHistFull[k].fHistV0v->SetBinError(j+1, totalError);
   } 
   else 
   {
    cout << "##### Resolution of the " << j+1 << "th harmonic was zero." << endl;
    // tracks hist
    fHistFull[k].fHistFullHar[j].fHistv2D->Reset();
    fHistFull[k].fHistFullHar[j].fHistvEta->Reset();
    fHistFull[k].fHistFullHar[j].fHistvPt->Reset();
    fHistFull[k].fHistv->SetBinContent(j+1, 0.);
    fHistFull[k].fHistv->SetBinError(j+1, 0.);
    // v0s hist
    fHistFull[k].fHistFullHar[j].fHistV0v2D->Reset();
    fHistFull[k].fHistFullHar[j].fHistV0vEta->Reset();
    fHistFull[k].fHistFullHar[j].fHistV0vPt->Reset();
    fHistFull[k].fHistFullHar[j].fHistV0sbv2D->Reset();
    fHistFull[k].fHistFullHar[j].fHistV0sbvEta->Reset();
    fHistFull[k].fHistFullHar[j].fHistV0sbvPt->Reset();
    fHistFull[k].fHistV0v->SetBinContent(j+1, 0.);
    fHistFull[k].fHistV0v->SetBinError(j+1, 0.);
   }
  }
 }

 return kTRUE ;
}
//-----------------------------------------------------------------------
Double_t AliFlowAnalyser::Chi(Double_t res)
{          
 // chi from the event plane resolution

 Double_t chi   = 2.0;
 Double_t delta = 1.0;
 for(int i = 0; i < 15; i++) 
 {
  if(ResEventPlane(chi) < res) { chi = chi + delta ; }
  else                         { chi = chi - delta ; }
  delta = delta / 2.;
 }

 return chi ;
}
//-----------------------------------------------------------------------
Double_t AliFlowAnalyser::ResEventPlane(Double_t chi)    
{
 // plane resolution as function of chi

 Double_t con = TMath::Sqrt(TMath::Pi()/2)/2 ;   // ~ 0.626657
 Double_t arg = chi * chi / 4.;
 Double_t res = con * chi * exp(-arg) * (TMath::BesselI0(arg) + TMath::BesselI1(arg)); 

 return res ;
}
//-----------------------------------------------------------------------
Double_t AliFlowAnalyser::ResEventPlaneK2(Double_t chi)  
{
 // ... for the case k=2.

 Double_t con = TMath::Sqrt(TMath::Pi()/2)/2 ;   // ~ 0.626657
 Double_t arg = chi * chi / 4.;

 Double_t besselOneHalf = TMath::Sqrt(arg/(TMath::Pi()/2)) * sinh(arg)/arg;
 Double_t besselThreeHalfs = TMath::Sqrt(arg/(TMath::Pi()/2)) * (cosh(arg)/arg - sinh(arg)/(arg*arg));
 Double_t res = con * chi * exp(-arg) * (besselOneHalf + besselThreeHalfs); 

 return res ;
}
//-----------------------------------------------------------------------
Double_t AliFlowAnalyser::ResEventPlaneK3(Double_t chi) 
{
 // ...for the case k=3.

 Double_t con = TMath::Sqrt(TMath::Pi()/2)/2 ;   // ~ 0.626657
 Double_t arg = chi * chi / 4.;
 Double_t res = con * chi * exp(-arg) * (TMath::BesselI1(arg) + TMath::BesselI(2, arg)); 

 return res ;
}
//-----------------------------------------------------------------------
Double_t AliFlowAnalyser::ResEventPlaneK4(Double_t chi) 
{
 // ...for the case k=4.

 Double_t con = TMath::Sqrt(TMath::Pi()/2)/2 ;   // ~ 0.626657
 Double_t arg = chi * chi / 4.;

 Double_t besselOneHalf = TMath::Sqrt(arg/(TMath::Pi()/2)) * sinh(arg)/arg;
 Double_t besselThreeHalfs = TMath::Sqrt(arg/(TMath::Pi()/2)) * (cosh(arg)/arg - sinh(arg)/(arg*arg));
 Double_t besselFiveHalfs = besselOneHalf - 3*besselThreeHalfs/arg;
 Double_t res = con * chi * exp(-arg) * (besselThreeHalfs + besselFiveHalfs); 

 return res ;
}
//-----------------------------------------------------------------------
// ###
//-----------------------------------------------------------------------


#endif
