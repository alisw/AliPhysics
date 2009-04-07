/*************************************************************************
* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
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

/********************************** 
 * flow analysis with Q-cumulants * 
 *                                * 
 * author:  Ante Bilandzic        * 
 *           (anteb@nikhef.nl)    *
 *********************************/ 

#define AliFlowAnalysisWithQCumulants_cxx

#include "Riostream.h"
#include "AliFlowCommonConstants.h"
#include "AliFlowCommonHist.h"
#include "AliFlowCommonHistResults.h"
#include "TChain.h"
#include "TFile.h"
#include "TList.h"
#include "TGraph.h"
#include "TParticle.h"
#include "TRandom3.h"
#include "TStyle.h"
#include "TProfile.h"
#include "TProfile2D.h" 
#include "TProfile3D.h"
#include "TMath.h"
#include "TArrow.h"
#include "TPaveLabel.h"
#include "TCanvas.h"
#include "AliFlowEventSimple.h"
#include "AliFlowTrackSimple.h"
#include "AliFlowAnalysisWithQCumulants.h"
#include "TArrayD.h"
#include "TRandom.h"

class TH1;
class TH2;
class TGraph;
class TPave;
class TLatex;
class TMarker;
class TRandom3;
class TObjArray;
class TList;
class TCanvas;
class TSystem;
class TROOT;
class AliFlowVector;
class TVector;

//================================================================================================================

ClassImp(AliFlowAnalysisWithQCumulants)

AliFlowAnalysisWithQCumulants::AliFlowAnalysisWithQCumulants():  
 fTrack(NULL),
 fHistList(NULL),
 fDiffFlowList(NULL),
 fWeightsList(NULL),
 fResultsList(NULL),
 fAvMultIntFlowQC(NULL),
 fQvectorComponents(NULL),
 fDiffFlowResults2ndOrderQC(NULL),
 fDiffFlowResults4thOrderQC(NULL),
 fCovariances(NULL),
 fQvectorForEachEventX(NULL),//to be removed
 fQvectorForEachEventY(NULL),//to be removed
 fQCorrelations(NULL),
 fQCorrelationsW(NULL),
 fQProduct(NULL),
 fDirectCorrelations(NULL),
 fDirectCorrelationsW(NULL),
 fDirectCorrelationsDiffFlow(NULL),
 fDirectCorrelationsDiffFlowW(NULL),
 f2PerPtBin1n1nPOI(NULL),
 f4PerPtBin1n1n1n1nPOI(NULL),
 f2PerEtaBin1n1nPOI(NULL),
 f4PerEtaBin1n1n1n1nPOI(NULL), 
 f2WPerPtBin1n1nPOI(NULL),
 f4WPerPtBin1n1n1n1nPOI(NULL),
 f2WPerEtaBin1n1nPOI(NULL),
 f4WPerEtaBin1n1n1n1nPOI(NULL),
 f2PerPtBin1n1nRP(NULL),
 f4PerPtBin1n1n1n1nRP(NULL),
 f2PerEtaBin1n1nRP(NULL),
 f4PerEtaBin1n1n1n1nRP(NULL),
 f2WPerPtBin1n1nRP(NULL),
 f4WPerPtBin1n1n1n1nRP(NULL),
 f2WPerEtaBin1n1nRP(NULL),
 f4WPerEtaBin1n1n1n1nRP(NULL),
 fCommonHists2nd(NULL),
 fCommonHists4th(NULL),
 fCommonHists6th(NULL),
 fCommonHists8th(NULL),
 fCommonHistsResults2nd(NULL),
 fCommonHistsResults4th(NULL),
 fCommonHistsResults6th(NULL),
 fCommonHistsResults8th(NULL),
 f2pDistribution(NULL),
 f4pDistribution(NULL),
 f6pDistribution(NULL),
 f8pDistribution(NULL),
 fnBinsPt(0),
 fPtMin(0),
 fPtMax(0),
 fnBinsEta(0),
 fEtaMin(0),
 fEtaMax(0),
 fEventCounter(0),
 fUsePhiWeights(kFALSE),
 fUsePtWeights(kFALSE),
 fUseEtaWeights(kFALSE),
 fUseWeights(kFALSE), 
 fUseWeightsBits(NULL), 

 // ...................................................................................................................
 // Q_{n,k} and S^M_{n,k}:    
 fReQ(NULL),
 fImQ(NULL),
 fSMpk(NULL),
 
 // q_n and m: 
 fReqnPtEta(NULL),  
 fImqnPtEta(NULL),
 fmPtEta(NULL),       
 
 // non-weighted q''_{n} and q''_{2n}:
 fReqPrimePrime1nPtEta(NULL),   
 fImqPrimePrime1nPtEta(NULL),
 fReqPrimePrime2nPtEta(NULL), 
 fImqPrimePrime2nPtEta(NULL), 
 
 // weighted q''_{n,2k} and q''_{2n,k}:
 fReqPrimePrime1n2kPtEta(NULL),   
 fImqPrimePrime1n2kPtEta(NULL),
 fReqPrimePrime2n1kPtEta(NULL), 
 fImqPrimePrime2n1kPtEta(NULL), 

 // m''   
 fmPrimePrimePtEta(NULL), 
 
 // S^{m''}_{n,k}
 fSmPrimePrime1p1kPtEta(NULL),
 fSmPrimePrime1p2kPtEta(NULL),
 fSmPrimePrime1p3kPtEta(NULL),
 
 // non-weighted q_RP{n} and q_RP{2n}:
 fReqRP1nPtEta(NULL), 
 fImqRP1nPtEta(NULL), 
 fReqRP2nPtEta(NULL), 
 fImqRP2nPtEta(NULL), 
  
 // weighted q_RP{n,2k} and q_RP{2n,k} (for each (pt,eta) bin for RPs)
 fReqRP1n2kPtEta(NULL), 
 fImqRP1n2kPtEta(NULL), 
 fReqRP2n1kPtEta(NULL),
 fImqRP2n1kPtEta(NULL), 
  
 // m_RP:
 fmRPPtEta(NULL), // # of particles which are RPs for each (pt,eta) bin
  
 // S^{m_RP}_{p,k} (for each (pt,eta) bin for RPs):
 fSmRP1p1kPtEta(NULL), 
 fSmRP1p2kPtEta(NULL), 
 fSmRP1p3kPtEta(NULL),
 
 // ----- RESULTS ----
 
 // non-weighted integrated flow:
 fIntFlowResultsQC(NULL),
 fIntFlowResultsPOIQC(NULL),
 fIntFlowResultsRPQC(NULL),
 
 // weighted integrated flow:
 fIntFlowResultsQCW(NULL),
 fIntFlowResultsPOIQCW(NULL),
 fIntFlowResultsRPQCW(NULL),

 // non-weighted correlations for each (pt,eta) bin for POIs:
 f2pPtEtaPOI(NULL),
 f4pPtEtaPOI(NULL),
 f6pPtEtaPOI(NULL),
 f8pPtEtaPOI(NULL),
 
 // non-weighted final results for differential flow for POIs:
 // 3D (pt,eta)
 fvn2ndPtEtaPOI(NULL),
 fvn4thPtEtaPOI(NULL),
 fvn6thPtEtaPOI(NULL),
 fvn8thPtEtaPOI(NULL),
 // 2D (pt)
 fvn2ndPtPOI(NULL),
 fvn4thPtPOI(NULL),
 fvn6thPtPOI(NULL),
 fvn8thPtPOI(NULL),
 // 2D (eta)
 fvn2ndEtaPOI(NULL),
 fvn4thEtaPOI(NULL),
 fvn6thEtaPOI(NULL),
 fvn8thEtaPOI(NULL),

 // weighted correlations for each (pt,eta) bin for POIs:
 f2pPtEtaPOIW(NULL),
 f4pPtEtaPOIW(NULL),
 f6pPtEtaPOIW(NULL),
 f8pPtEtaPOIW(NULL),
 
 // weighted final results for differential flow for POIs:
 // 3D (pt,eta)
 fvn2ndPtEtaPOIW(NULL),
 fvn4thPtEtaPOIW(NULL),
 fvn6thPtEtaPOIW(NULL),
 fvn8thPtEtaPOIW(NULL),
 // 2D (pt)
 fvn2ndPtPOIW(NULL),
 fvn4thPtPOIW(NULL),
 fvn6thPtPOIW(NULL),
 fvn8thPtPOIW(NULL),
 // 2D (eta)
 fvn2ndEtaPOIW(NULL),
 fvn4thEtaPOIW(NULL),
 fvn6thEtaPOIW(NULL),
 fvn8thEtaPOIW(NULL),
 
 // non-weighted correlations for each (pt,eta) bin for RPs:
 f2pPtEtaRP(NULL),
 f4pPtEtaRP(NULL),
 f6pPtEtaRP(NULL),
 f8pPtEtaRP(NULL),
 
 // non-weighted final results for differential flow for RPs:
 // 3D (pt,eta)
 fvn2ndPtEtaRP(NULL),
 fvn4thPtEtaRP(NULL),
 fvn6thPtEtaRP(NULL),
 fvn8thPtEtaRP(NULL),
 // 2D (pt)
 fvn2ndPtRP(NULL),
 fvn4thPtRP(NULL),
 fvn6thPtRP(NULL),
 fvn8thPtRP(NULL),
 // 2D (eta)
 fvn2ndEtaRP(NULL),
 fvn4thEtaRP(NULL),
 fvn6thEtaRP(NULL),
 fvn8thEtaRP(NULL),

 // weighted correlations for each (pt,eta) bin for RPs:
 f2pPtEtaRPW(NULL),
 f4pPtEtaRPW(NULL),
 f6pPtEtaRPW(NULL),
 f8pPtEtaRPW(NULL),
 
 // weighted final results for differential flow for RPs:
 // 3D (pt,eta)
 fvn2ndPtEtaRPW(NULL),
 fvn4thPtEtaRPW(NULL),
 fvn6thPtEtaRPW(NULL),
 fvn8thPtEtaRPW(NULL),
 // 2D (pt)
 fvn2ndPtRPW(NULL),
 fvn4thPtRPW(NULL),
 fvn6thPtRPW(NULL),
 fvn8thPtRPW(NULL),
 // 2D (eta)
 fvn2ndEtaRPW(NULL),
 fvn4thEtaRPW(NULL),
 fvn6thEtaRPW(NULL),
 fvn8thEtaRPW(NULL)
 // ...................................................................................................................
 
{
 // constructor 
 fHistList = new TList();
 fDiffFlowList = new TList(); 
 fDiffFlowList->SetName("DifferentialFlow"); 
 fWeightsList = new TList();
 fWeightsList->SetName("Weights");
 fResultsList = new TList();
 fResultsList->SetName("Results");
  
 fnBinsPt = AliFlowCommonConstants::GetNbinsPt();
 fPtMin   = AliFlowCommonConstants::GetPtMin();	     
 fPtMax   = AliFlowCommonConstants::GetPtMax();
 
 fnBinsEta = AliFlowCommonConstants::GetNbinsEta();
 fEtaMin   = AliFlowCommonConstants::GetEtaMin();	     
 fEtaMax   = AliFlowCommonConstants::GetEtaMax();
}

AliFlowAnalysisWithQCumulants::~AliFlowAnalysisWithQCumulants()
{
 //destructor
 delete fHistList; 
 delete fDiffFlowList;
 delete fWeightsList; 
 delete fResultsList; 
}

//================================================================================================================

void AliFlowAnalysisWithQCumulants::Init()
{
 //various output histograms
 //avarage multiplicity 
 fAvMultIntFlowQC = new TProfile("fAvMultIntFlowQC","Average Multiplicity",1,0,1,"s");
 fAvMultIntFlowQC->SetXTitle("");
 fAvMultIntFlowQC->SetYTitle("");
 fAvMultIntFlowQC->SetLabelSize(0.06);
 fAvMultIntFlowQC->SetMarkerStyle(25);
 fAvMultIntFlowQC->SetLabelOffset(0.01);
 (fAvMultIntFlowQC->GetXaxis())->SetBinLabel(1,"Average Multiplicity");
 fHistList->Add(fAvMultIntFlowQC);
 
 //Q-vector stuff
 fQvectorComponents = new TProfile("fQvectorComponents","Avarage of Q-vector components",44,0.,44.,"s");
 fQvectorComponents->SetXTitle("");
 fQvectorComponents->SetYTitle("");
 //fHistList->Add(fQvectorComponents);
 
 //final results for differential flow from 2nd order Q-cumulant
 fDiffFlowResults2ndOrderQC = new TH1D("fDiffFlowResults2ndOrderQC","Differential Flow from 2nd Order Q-cumulant",fnBinsPt,fPtMin,fPtMax);
 fDiffFlowResults2ndOrderQC->SetXTitle("p_{t} [GeV]");
 //fDiffFlowResults2ndOrderQC->SetYTitle("Differential Flow");
 fHistList->Add(fDiffFlowResults2ndOrderQC);
 
 //final results for differential flow from 4th order Q-cumulant
 fDiffFlowResults4thOrderQC = new TH1D("fDiffFlowResults4thOrderQC","Differential Flow from 4th Order Q-cumulant",fnBinsPt,fPtMin,fPtMax);
 fDiffFlowResults4thOrderQC->SetXTitle("p_{t} [GeV]");
 //fDiffFlowResults4thOrderQC->SetYTitle("Differential Flow");
 fHistList->Add(fDiffFlowResults4thOrderQC);
 
 //final results for covariances (1st bin: <2*4>-<2>*<4>, 2nd bin: <2*6>-<2>*<6>, ...)
 fCovariances = new TH1D("fCovariances","Covariances",6,0,6);
 //fCovariances->SetXTitle("");
 //fCovariances->SetYTitle("<covariance>");
 fCovariances->SetLabelSize(0.04);
 fCovariances->SetTickLength(1);
 fCovariances->SetMarkerStyle(25);
 (fCovariances->GetXaxis())->SetBinLabel(1,"Cov(2,4)");
 (fCovariances->GetXaxis())->SetBinLabel(2,"Cov(2,6)");
 (fCovariances->GetXaxis())->SetBinLabel(3,"Cov(2,8)");
 (fCovariances->GetXaxis())->SetBinLabel(4,"Cov(4,6)");
 (fCovariances->GetXaxis())->SetBinLabel(5,"Cov(4,8)");
 (fCovariances->GetXaxis())->SetBinLabel(6,"Cov(6,8)");
 fHistList->Add(fCovariances);
  
 //xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 //                   !!!! to be removed !!!! 
 //profile containing the x-components of Q-vectors from all events 
 fQvectorForEachEventX = new TProfile("fQvectorForEachEventX","x-components of Q-vectors",44000,1,44000,"s");
 fHistList->Add(fQvectorForEachEventX);
 
 //profile containing the y-components of Q-vectors from all events 
 fQvectorForEachEventY = new TProfile("fQvectorForEachEventY","y-components of Q-vectors",44000,1,44000,"s");
 fHistList->Add(fQvectorForEachEventY);
 //xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    
 // multi-particle correlations calculated from Q-vectors
 fQCorrelations = new TProfile("fQCorrelations","multi-particle correlations from Q-vectors",32,0,32,"s");
 fQCorrelations->SetTickLength(-0.01,"Y");
 fQCorrelations->SetMarkerStyle(25);
 fQCorrelations->SetLabelSize(0.03);
 fQCorrelations->SetLabelOffset(0.01,"Y");
 // 2-p:
 (fQCorrelations->GetXaxis())->SetBinLabel(1,"<<2>>_{n|n}");
 (fQCorrelations->GetXaxis())->SetBinLabel(2,"<<2>>_{2n|2n}");
 (fQCorrelations->GetXaxis())->SetBinLabel(3,"<<2>>_{3n|3n}");
 (fQCorrelations->GetXaxis())->SetBinLabel(4,"<<2>>_{4n|4n}");
 // 3-p:
 (fQCorrelations->GetXaxis())->SetBinLabel(6,"<<3>>_{2n|n,n}");
 (fQCorrelations->GetXaxis())->SetBinLabel(7,"<<3>>_{3n|2n,n}");
 (fQCorrelations->GetXaxis())->SetBinLabel(8,"<<3>>_{4n|2n,2n}");
 (fQCorrelations->GetXaxis())->SetBinLabel(9,"<<3>>_{4n|3n,n}");
 // 4-p:
 (fQCorrelations->GetXaxis())->SetBinLabel(11,"<<4>>_{n,n|n,n}"); 
 (fQCorrelations->GetXaxis())->SetBinLabel(12,"<<4>>_{2n,n|2n,n}");
 (fQCorrelations->GetXaxis())->SetBinLabel(13,"<<4>>_{2n,2n|2n,2n}");
 (fQCorrelations->GetXaxis())->SetBinLabel(14,"<<4>>_{3n|n,n,n}");
 (fQCorrelations->GetXaxis())->SetBinLabel(15,"<<4>>_{3n,n|3n,n}");
 (fQCorrelations->GetXaxis())->SetBinLabel(16,"<<4>>_{3n,n|2n,2n}"); 
 (fQCorrelations->GetXaxis())->SetBinLabel(17,"<<4>>_{4n|2n,n,n}");
 // 5-p:
 (fQCorrelations->GetXaxis())->SetBinLabel(19,"<<5>>_{2n|n,n,n,n}"); 
 (fQCorrelations->GetXaxis())->SetBinLabel(20,"<<5>>_{2n,2n|2n,n,n}");
 (fQCorrelations->GetXaxis())->SetBinLabel(21,"<<5>>_{3n,n|2n,n,n}");
 (fQCorrelations->GetXaxis())->SetBinLabel(22,"<<5>>_{4n|n,n,n,n}");
 // 6-p:
 (fQCorrelations->GetXaxis())->SetBinLabel(24,"<<6>>_{n,n,n|n,n,n}");
 (fQCorrelations->GetXaxis())->SetBinLabel(25,"<<6>>_{2n,n,n|2n,n,n}");
 (fQCorrelations->GetXaxis())->SetBinLabel(26,"<<6>>_{2n,2n|n,n,n,n}");
 (fQCorrelations->GetXaxis())->SetBinLabel(27,"<<6>>_{3n,n|n,n,n,n}");
 // 7-p:
 (fQCorrelations->GetXaxis())->SetBinLabel(29,"<<7>>_{2n,n,n|n,n,n,n}");
 // 8-p:
 (fQCorrelations->GetXaxis())->SetBinLabel(31,"<<8>>_{n,n,n,n|n,n,n,n}");
 // add fQCorrelations to the main list:
 fHistList->Add(fQCorrelations);
 
 //.........................................................................
 //weighted multi-particle correlations calculated from Q-vectors
 fQCorrelationsW = new TProfile("fQCorrelationsW","weighted multi-particle correlations from Q-vectors",200,0,200,"s");
 fQCorrelationsW->SetTickLength(-0.01,"Y");
 fQCorrelationsW->SetMarkerStyle(25);
 fQCorrelationsW->SetLabelSize(0.03);
 fQCorrelationsW->SetLabelOffset(0.01,"Y");
 // 2-p:
 (fQCorrelationsW->GetXaxis())->SetBinLabel(1,"<w_{1}w_{2}cos(n(#phi_{1}-#phi_{2}))>");
 (fQCorrelationsW->GetXaxis())->SetBinLabel(2,"<w_{1}^{2}w_{2}^{2}cos(2n(#phi_{1}-#phi_{2}))>");
 (fQCorrelationsW->GetXaxis())->SetBinLabel(3,"<w_{1}^{3}w_{2}^{3}cos(3n(#phi_{1}-#phi_{2}))>");
 (fQCorrelationsW->GetXaxis())->SetBinLabel(4,"<w_{1}^{4}w_{2}^{4}cos(4n(#phi_{1}-#phi_{2}))>");
 (fQCorrelationsW->GetXaxis())->SetBinLabel(5,"<w_{1}^{3}w_{2}cos(n(#phi_{1}-#phi_{2}))>");
 (fQCorrelationsW->GetXaxis())->SetBinLabel(6,"<w_{1}^{2}w_{2}w_{3}cos(n(#phi_{1}-#phi_{2}))>");
 // 3-p:
 (fQCorrelationsW->GetXaxis())->SetBinLabel(21,"<w_{1}w_{2}w_{3}^{2}cos(n(2#phi_{1}-#phi_{2}-#phi_{3}))>");
 // 4-p:
 (fQCorrelationsW->GetXaxis())->SetBinLabel(41,"<w_{1}w_{2}w_{3}w_{4}cos(n(#phi_{1}+#phi_{2}-#phi_{3}-#phi_{4}))>");
 // add fQCorrelationsW to the main list:
 fHistList->Add(fQCorrelationsW);
 //.........................................................................
 
 
 //average products
 fQProduct = new TProfile("fQProduct","average of products",6,0,6,"s");
 fQProduct->SetTickLength(-0.01,"Y");
 fQProduct->SetMarkerStyle(25);
 fQProduct->SetLabelSize(0.03);
 fQProduct->SetLabelOffset(0.01,"Y");
 (fQProduct->GetXaxis())->SetBinLabel(1,"<<2*4>>");
 (fQProduct->GetXaxis())->SetBinLabel(2,"<<2*6>>");
 (fQProduct->GetXaxis())->SetBinLabel(3,"<<2*8>>");
 (fQProduct->GetXaxis())->SetBinLabel(4,"<<4*6>>");
 (fQProduct->GetXaxis())->SetBinLabel(5,"<<4*8>>");
 (fQProduct->GetXaxis())->SetBinLabel(6,"<<6*8>>");
 fQProduct->SetXTitle("");
 fQProduct->SetYTitle("");
 fHistList->Add(fQProduct);
 
 // multi-particle correlations calculated with nested loops (needed for int. flow)
 fDirectCorrelations = new TProfile("fDirectCorrelations","multi-particle correlations with nested loops",100,0,100,"s");
 fDirectCorrelations->SetXTitle("");
 fDirectCorrelations->SetYTitle("correlations");
 fHistList->Add(fDirectCorrelations);
 
 // multi-particle correlations calculated with nested loops (needed for weighted int. flow)
 fDirectCorrelationsW = new TProfile("fDirectCorrelationsW","multi-particle correlations with nested loops",200,0,200,"s");
 fDirectCorrelationsW->SetXTitle("");
 fDirectCorrelationsW->SetYTitle("correlations");
 fHistList->Add(fDirectCorrelationsW);
 
 // multi-particle correlations calculated with nested loops (needed for diff. flow)
 fDirectCorrelationsDiffFlow = new TProfile("fDirectCorrelationsDiffFlow","multi-particle correlations with nested loops",200,0,200,"s");
 fDirectCorrelationsDiffFlow->SetXTitle("");
 fDirectCorrelationsDiffFlow->SetYTitle("correlations");
 fHistList->Add(fDirectCorrelationsDiffFlow);
 
 // multi-particle correlations calculated with nested loops (needed for weighted diff. flow)
 fDirectCorrelationsDiffFlowW = new TProfile("fDirectCorrelationsDiffFlowW","multi-particle correlations with nested loops",200,0,200,"s");
 fDirectCorrelationsDiffFlowW->SetXTitle("");
 fDirectCorrelationsDiffFlowW->SetYTitle("correlations");
 fHistList->Add(fDirectCorrelationsDiffFlowW);
 
 //f2PerPtBin1n1nRP
 f2PerPtBin1n1nRP = new TProfile("f2PerPtBin1n1nRP","<2'>_{n|n}",fnBinsPt,fPtMin,fPtMax,"s");
 f2PerPtBin1n1nRP->SetXTitle("p_{t} [GeV]");
 fDiffFlowList->Add(f2PerPtBin1n1nRP);
 
 //f4PerPtBin1n1n1n1nRP
 f4PerPtBin1n1n1n1nRP = new TProfile("f4PerPtBin1n1n1n1nRP","<4'>_{n,n|n,n}",fnBinsPt,fPtMin,fPtMax,"s");
 f4PerPtBin1n1n1n1nRP->SetXTitle("p_{t} [GeV]");
 fDiffFlowList->Add(f4PerPtBin1n1n1n1nRP);
 
 //f2PerEtaBin1n1nRP
 f2PerEtaBin1n1nRP = new TProfile("f2PerEtaBin1n1nRP","<2'>_{n|n}",fnBinsEta,fEtaMin,fEtaMax,"s");
 f2PerEtaBin1n1nRP->SetXTitle("#eta");
 fDiffFlowList->Add(f2PerEtaBin1n1nRP);
 
 //f4PerEtaBin1n1n1n1nRP
 f4PerEtaBin1n1n1n1nRP = new TProfile("f4PerEtaBin1n1n1n1nRP","<4'>_{n,n|n,n}",fnBinsEta,fEtaMin,fEtaMax,"s");
 f4PerEtaBin1n1n1n1nRP->SetXTitle("#eta");
 fDiffFlowList->Add(f4PerEtaBin1n1n1n1nRP);
 
 //f2PerPtBin1n1nPOI
 f2PerPtBin1n1nPOI = new TProfile("f2PerPtBin1n1nPOI","<2'>_{n|n}",fnBinsPt,fPtMin,fPtMax,"s");
 f2PerPtBin1n1nPOI->SetXTitle("#eta");
 fDiffFlowList->Add(f2PerPtBin1n1nPOI);
 
 //f4PerPtBin1n1n1n1nPOI
 f4PerPtBin1n1n1n1nPOI = new TProfile("f4PerPtBin1n1n1n1nPOI","<4'>_{n,n|n,n}",fnBinsPt,fPtMin,fPtMax,"s");
 f4PerPtBin1n1n1n1nPOI->SetXTitle("p_{t} [GeV]");
 fDiffFlowList->Add(f4PerPtBin1n1n1n1nPOI);
 
 //f2PerEtaBin1n1nPOI
 f2PerEtaBin1n1nPOI = new TProfile("f2PerEtaBin1n1nPOI","<2'>_{n|n}",fnBinsEta,fEtaMin,fEtaMax,"s");
 f2PerEtaBin1n1nPOI->SetXTitle("#eta");
 fDiffFlowList->Add(f2PerEtaBin1n1nPOI);
 
 //f4PerEtaBin1n1n1n1nPOI
 f4PerEtaBin1n1n1n1nPOI = new TProfile("f4PerEtaBin1n1n1n1nPOI","<4'>_{n,n|n,n}",fnBinsEta,fEtaMin,fEtaMax,"s");
 f4PerEtaBin1n1n1n1nPOI->SetXTitle("#eta");
 fDiffFlowList->Add(f4PerEtaBin1n1n1n1nPOI);
 
 //f2WPerPtBin1n1nPOI
 f2WPerPtBin1n1nPOI = new TProfile("f2WPerPtBin1n1nPOI","<2'>_{n|n}",fnBinsPt,fPtMin,fPtMax,"s");
 f2WPerPtBin1n1nPOI->SetXTitle("#pt");
 fDiffFlowList->Add(f2WPerPtBin1n1nPOI);
 
 //f4WPerPtBin1n1n1n1nPOI
 f4WPerPtBin1n1n1n1nPOI = new TProfile("f4WPerPtBin1n1n1n1nPOI","<4'>_{n,n|n,n}",fnBinsPt,fPtMin,fPtMax,"s");
 f4WPerPtBin1n1n1n1nPOI->SetXTitle("#Pt");
 fDiffFlowList->Add(f4WPerPtBin1n1n1n1nPOI);
 
 //f2WPerEtaBin1n1nPOI
 f2WPerEtaBin1n1nPOI = new TProfile("f2WPerEtaBin1n1nPOI","<2'>_{n|n}",fnBinsEta,fEtaMin,fEtaMax,"s");
 f2WPerEtaBin1n1nPOI->SetXTitle("#eta");
 fDiffFlowList->Add(f2WPerEtaBin1n1nPOI);
 
 //f4WPerEtaBin1n1n1n1nPOI
 f4WPerEtaBin1n1n1n1nPOI = new TProfile("f4WPerEtaBin1n1n1n1nPOI","<4'>_{n,n|n,n}",fnBinsEta,fEtaMin,fEtaMax,"s");
 f4WPerEtaBin1n1n1n1nPOI->SetXTitle("#eta");
 fDiffFlowList->Add(f4WPerEtaBin1n1n1n1nPOI);
 
 //f2WPerPtBin1n1nRP
 f2WPerPtBin1n1nRP = new TProfile("f2WPerPtBin1n1nRP","<2'>_{n|n}",fnBinsPt,fPtMin,fPtMax,"s");
 f2WPerPtBin1n1nRP->SetXTitle("#pt");
 fDiffFlowList->Add(f2WPerPtBin1n1nRP);
 
 //f4WPerPtBin1n1n1n1nRP
 f4WPerPtBin1n1n1n1nRP = new TProfile("f4WPerPtBin1n1n1n1nRP","<4'>_{n,n|n,n}",fnBinsPt,fPtMin,fPtMax,"s");
 f4WPerPtBin1n1n1n1nRP->SetXTitle("#Pt");
 fDiffFlowList->Add(f4WPerPtBin1n1n1n1nRP);
 
 //f2WPerEtaBin1n1nRP
 f2WPerEtaBin1n1nRP = new TProfile("f2WPerEtaBin1n1nRP","<2'>_{n|n}",fnBinsEta,fEtaMin,fEtaMax,"s");
 f2WPerEtaBin1n1nRP->SetXTitle("#eta");
 fDiffFlowList->Add(f2WPerEtaBin1n1nRP);
 
 //f4WPerEtaBin1n1n1n1nRP
 f4WPerEtaBin1n1n1n1nRP = new TProfile("f4WPerEtaBin1n1n1n1nRP","<4'>_{n,n|n,n}",fnBinsEta,fEtaMin,fEtaMax,"s");
 f4WPerEtaBin1n1n1n1nRP->SetXTitle("#eta");
 fDiffFlowList->Add(f4WPerEtaBin1n1n1n1nRP);
 
 //common control histogram (2nd order)
 fCommonHists2nd = new AliFlowCommonHist("AliFlowCommonHist2ndOrderQC");
 fHistList->Add(fCommonHists2nd);  
 
 //common control histogram (4th order)
 fCommonHists4th = new AliFlowCommonHist("AliFlowCommonHist4thOrderQC");
 fHistList->Add(fCommonHists4th);  
 
 //common control histogram (6th order)
 fCommonHists6th = new AliFlowCommonHist("AliFlowCommonHist6thOrderQC");
 fHistList->Add(fCommonHists6th);  
 
 //common control histogram (8th order)
 fCommonHists8th = new AliFlowCommonHist("AliFlowCommonHist8thOrderQC");
 fHistList->Add(fCommonHists8th);  
  
 //common histograms for final results (2nd order)
 fCommonHistsResults2nd = new AliFlowCommonHistResults("AliFlowCommonHistResults2ndOrderQC");
 fHistList->Add(fCommonHistsResults2nd); 
 
 //common histograms for final results (4th order)
 fCommonHistsResults4th = new AliFlowCommonHistResults("AliFlowCommonHistResults4thOrderQC");
 fHistList->Add(fCommonHistsResults4th);
 
 //common histograms for final results (6th order)
 fCommonHistsResults6th = new AliFlowCommonHistResults("AliFlowCommonHistResults6thOrderQC");
 fHistList->Add(fCommonHistsResults6th); 
 
 //common histograms for final results (8th order)
 fCommonHistsResults8th = new AliFlowCommonHistResults("AliFlowCommonHistResults8thOrderQC");
 fHistList->Add(fCommonHistsResults8th); 
 
 //weighted <2>_{n|n} distribution
 f2pDistribution = new TH1D("f2pDistribution","<2>_{n|n} distribution",100000,-0.02,0.1);
 f2pDistribution->SetXTitle("<2>_{n|n}");
 f2pDistribution->SetYTitle("Counts");
 fHistList->Add(f2pDistribution);

 //weighted <4>_{n,n|n,n} distribution
 f4pDistribution = new TH1D("f4pDistribution","<4>_{n,n|n,n} distribution",100000,-0.00025,0.002);
 f4pDistribution->SetXTitle("<4>_{n,n|n,n}");
 f4pDistribution->SetYTitle("Counts");
 fHistList->Add(f4pDistribution); 
 
 //weighted <6>_{n,n,n|n,n,n} distribution
 f6pDistribution = new TH1D("f6pDistribution","<6>_{n,n,n|n,n,n} distribution",100000,-0.000005,0.000025);
 f6pDistribution->SetXTitle("<6>_{n,n,n|n,n,n}");
 f6pDistribution->SetYTitle("Counts");
 fHistList->Add(f6pDistribution);
 
 //weighted <8>_{n,n,n,n|n,n,n,n} distribution
 f8pDistribution = new TH1D("f8pDistribution","<8>_{n,n,n,n|n,n,n,n} distribution",100000,-0.000000001,0.00000001);
 f8pDistribution->SetXTitle("<8>_{n,n,n,n|n,n,n,n}");
 f8pDistribution->SetYTitle("Counts");
 fHistList->Add(f8pDistribution);
 
 
 
 
 
 
 
 
 // .......................................................................................................................................
 // Q_{n,k} and S^M_{n,k}:    
 fReQ  = new TMatrixD(4,9);
 fImQ  = new TMatrixD(4,9);
 fSMpk = new TMatrixD(8,9);
 
 // q'_{n}:
 fReqnPtEta = new TH2D("fReqnPtEta","Re[q_{n}(p_{t},#eta)]",fnBinsPt,fPtMin,fPtMax,fnBinsEta,fEtaMin,fEtaMax);  
 fImqnPtEta = new TH2D("fImqnPtEta","Im[q_{n}(p_{t},#eta)]",fnBinsPt,fPtMin,fPtMax,fnBinsEta,fEtaMin,fEtaMax); 
 fmPtEta    = new TH2D("fmPtEta","m(p_{t},#eta)",fnBinsPt,fPtMin,fPtMax,fnBinsEta,fEtaMin,fEtaMax);      
 
 // non-weighted q''_{n} and q''_{2n}:
 fReqPrimePrime1nPtEta = new TH2D("fReqPrimePrime1nPtEta","Re[q''_{n}(p_{t},#eta)]",fnBinsPt,fPtMin,fPtMax,fnBinsEta,fEtaMin,fEtaMax);   
 fImqPrimePrime1nPtEta = new TH2D("fImqPrimePrime1nPtEta","Im[q''_{n}(p_{t},#eta)]",fnBinsPt,fPtMin,fPtMax,fnBinsEta,fEtaMin,fEtaMax);   
 fReqPrimePrime2nPtEta = new TH2D("fReqPrimePrime2nPtEta","Re[q''_{2n}(p_{t},#eta)]",fnBinsPt,fPtMin,fPtMax,fnBinsEta,fEtaMin,fEtaMax); 
 fImqPrimePrime2nPtEta = new TH2D("fImqPrimePrime2nPtEta","Im[q''_{2n}(p_{t},#eta)]",fnBinsPt,fPtMin,fPtMax,fnBinsEta,fEtaMin,fEtaMax); 

 // weighted q''_{n,2k} and q''_{2n,k}:
 fReqPrimePrime1n2kPtEta = new TH2D("fReqPrimePrime1n2kPtEta","Re[q''_{n,2}(p_{t},#eta)]",fnBinsPt,fPtMin,fPtMax,fnBinsEta,fEtaMin,fEtaMax);   
 fImqPrimePrime1n2kPtEta = new TH2D("fImqPrimePrime1n2kPtEta","Im[q''_{n,2}(p_{t},#eta)]",fnBinsPt,fPtMin,fPtMax,fnBinsEta,fEtaMin,fEtaMax);   
 fReqPrimePrime2n1kPtEta = new TH2D("fReqPrimePrime2n1kPtEta","Re[q''_{2n,1(p_{t},#eta)}]",fnBinsPt,fPtMin,fPtMax,fnBinsEta,fEtaMin,fEtaMax); 
 fImqPrimePrime2n1kPtEta = new TH2D("fImqPrimePrime2n1kPtEta","Im[q''_{2n,1}(p_{t},#eta)]",fnBinsPt,fPtMin,fPtMax,fnBinsEta,fEtaMin,fEtaMax);  
 
 // m'':
 fmPrimePrimePtEta = new TH2D("fmPrimePrimePtEta","m''(p_{t},#eta)",fnBinsPt,fPtMin,fPtMax,fnBinsEta,fEtaMin,fEtaMax); 
 
 // S^{m''}_{p,k}:
 fSmPrimePrime1p1kPtEta = new TH2D("fSmPrimePrime1p1kPtEta","S^{m''}_{1,1}(p_{t},#eta)",fnBinsPt,fPtMin,fPtMax,fnBinsEta,fEtaMin,fEtaMax); 
 fSmPrimePrime1p2kPtEta = new TH2D("fSmPrimePrime1p2kPtEta","S^{m''}_{1,2}(p_{t},#eta)",fnBinsPt,fPtMin,fPtMax,fnBinsEta,fEtaMin,fEtaMax);
 fSmPrimePrime1p3kPtEta = new TH2D("fSmPrimePrime1p3kPtEta","S^{m''}_{1,3}(p_{t},#eta)",fnBinsPt,fPtMin,fPtMax,fnBinsEta,fEtaMin,fEtaMax);
 
 // non-weighted q_RP{n} and q_RP{2n}:
 fReqRP1nPtEta = new TH2D("fReqRP1nPtEta","Re[q_{n}(p_{t},#eta)] for RPs",fnBinsPt,fPtMin,fPtMax,fnBinsEta,fEtaMin,fEtaMax);   
 fImqRP1nPtEta = new TH2D("fImqRP1nPtEta","Im[q_{n}(p_{t},#eta)] for RPs",fnBinsPt,fPtMin,fPtMax,fnBinsEta,fEtaMin,fEtaMax);   
 fReqRP2nPtEta = new TH2D("fReqRP2nPtEta","Re[q_{2n}(p_{t},#eta)] for RPs",fnBinsPt,fPtMin,fPtMax,fnBinsEta,fEtaMin,fEtaMax); 
 fImqRP2nPtEta = new TH2D("fImqRP2nPtEta","Im[q_{2n}(p_{t},#eta)] for RPs",fnBinsPt,fPtMin,fPtMax,fnBinsEta,fEtaMin,fEtaMax); 

 // weighted q_RP{n,2k} and q_RP{2n,k}:
 fReqRP1n2kPtEta = new TH2D("fReqRP1n2kPtEta","Re[q_{n,2}(p_{t},#eta)] for RPs",fnBinsPt,fPtMin,fPtMax,fnBinsEta,fEtaMin,fEtaMax);   
 fImqRP1n2kPtEta = new TH2D("fImqRP1n2kPtEta","Im[q_{n,2}(p_{t},#eta)] for RPs",fnBinsPt,fPtMin,fPtMax,fnBinsEta,fEtaMin,fEtaMax);   
 fReqRP2n1kPtEta = new TH2D("fReqRP2n1kPtEta","Re[q_{2n,1(p_{t},#eta)}] for RPs",fnBinsPt,fPtMin,fPtMax,fnBinsEta,fEtaMin,fEtaMax); 
 fImqRP2n1kPtEta = new TH2D("fImqRP2n1kPtEta","Im[q_{2n,1}(p_{t},#eta)] for RPs",fnBinsPt,fPtMin,fPtMax,fnBinsEta,fEtaMin,fEtaMax);  
 
 // mRP:
 fmRPPtEta = new TH2D("fmRPPtEta","m(p_{t},#eta) for RPs",fnBinsPt,fPtMin,fPtMax,fnBinsEta,fEtaMin,fEtaMax); 
 
 // S^{mRP}_{p,k}:
 fSmRP1p1kPtEta = new TH2D("fSmRP1p1kPtEta","S^{m}_{1,1}(p_{t},#eta) for RPs",fnBinsPt,fPtMin,fPtMax,fnBinsEta,fEtaMin,fEtaMax); 
 fSmRP1p2kPtEta = new TH2D("fSmRP1p2kPtEta","S^{m}_{1,2}(p_{t},#eta) for RPs",fnBinsPt,fPtMin,fPtMax,fnBinsEta,fEtaMin,fEtaMax);
 fSmRP1p3kPtEta = new TH2D("fSmRP1p3kPtEta","S^{m}_{1,3}(p_{t},#eta) for RPs",fnBinsPt,fPtMin,fPtMax,fnBinsEta,fEtaMin,fEtaMax);
 
 // ----- RESULTS ----
 
 // final results for non-weighted no-name integrated flow:
 fIntFlowResultsQC = new TH1D("fIntFlowResultsQC","Integrated Flow from Q-cumulants",4,0,4);
 fIntFlowResultsQC->SetLabelSize(0.06);
 fIntFlowResultsQC->SetMarkerStyle(25);
 (fIntFlowResultsQC->GetXaxis())->SetBinLabel(1,"v_{n}{2}");
 (fIntFlowResultsQC->GetXaxis())->SetBinLabel(2,"v_{n}{4}");
 (fIntFlowResultsQC->GetXaxis())->SetBinLabel(3,"v_{n}{6}");
 (fIntFlowResultsQC->GetXaxis())->SetBinLabel(4,"v_{n}{8}");
 fResultsList->Add(fIntFlowResultsQC);
 
 // final results for non-weighted POIs integrated flow:
 fIntFlowResultsPOIQC = new TH1D("fIntFlowResultsPOIQC","Integrated Flow (POI) from Q-cumulants",4,0,4);
 fIntFlowResultsPOIQC->SetLabelSize(0.06);
 fIntFlowResultsPOIQC->SetMarkerStyle(25);
 (fIntFlowResultsPOIQC->GetXaxis())->SetBinLabel(1,"v_{n}{2}");
 (fIntFlowResultsPOIQC->GetXaxis())->SetBinLabel(2,"v_{n}{4}");
 (fIntFlowResultsPOIQC->GetXaxis())->SetBinLabel(3,"v_{n}{6}");
 (fIntFlowResultsPOIQC->GetXaxis())->SetBinLabel(4,"v_{n}{8}");
 fResultsList->Add(fIntFlowResultsPOIQC);
 
 // final results for non-weighted RPs integrated flow:
 fIntFlowResultsRPQC = new TH1D("fIntFlowResultsRPQC","Integrated Flow (RP) from Q-cumulants",4,0,4);
 fIntFlowResultsRPQC->SetLabelSize(0.06);
 fIntFlowResultsRPQC->SetMarkerStyle(25);
 (fIntFlowResultsRPQC->GetXaxis())->SetBinLabel(1,"v_{n}{2}");
 (fIntFlowResultsRPQC->GetXaxis())->SetBinLabel(2,"v_{n}{4}");
 (fIntFlowResultsRPQC->GetXaxis())->SetBinLabel(3,"v_{n}{6}");
 (fIntFlowResultsRPQC->GetXaxis())->SetBinLabel(4,"v_{n}{8}");
 fResultsList->Add(fIntFlowResultsRPQC);
 
 // final results for weighted no-name integrated flow:
 fIntFlowResultsQCW = new TH1D("fIntFlowResultsQCW","Integrated Flow from Q-cumulants with Weights",4,0,4);
 fIntFlowResultsQCW->SetLabelSize(0.06);
 fIntFlowResultsQCW->SetMarkerStyle(25);
 (fIntFlowResultsQCW->GetXaxis())->SetBinLabel(1,"v_{n}{2}");
 (fIntFlowResultsQCW->GetXaxis())->SetBinLabel(2,"v_{n}{4}");
 (fIntFlowResultsQCW->GetXaxis())->SetBinLabel(3,"v_{n}{6}");
 (fIntFlowResultsQCW->GetXaxis())->SetBinLabel(4,"v_{n}{8}");
 fResultsList->Add(fIntFlowResultsQCW);
 
 // final results for weighted POIs integrated flow:
 fIntFlowResultsPOIQCW = new TH1D("fIntFlowResultsPOIQCW","Integrated Flow (POI) from Q-cumulants with Weights",4,0,4);
 fIntFlowResultsPOIQCW->SetLabelSize(0.06);
 fIntFlowResultsPOIQCW->SetMarkerStyle(25);
 (fIntFlowResultsPOIQCW->GetXaxis())->SetBinLabel(1,"v_{n}{2}");
 (fIntFlowResultsPOIQCW->GetXaxis())->SetBinLabel(2,"v_{n}{4}");
 (fIntFlowResultsPOIQCW->GetXaxis())->SetBinLabel(3,"v_{n}{6}");
 (fIntFlowResultsPOIQCW->GetXaxis())->SetBinLabel(4,"v_{n}{8}");
 fResultsList->Add(fIntFlowResultsPOIQCW);
 
 // final results for weighted RPs integrated flow:
 fIntFlowResultsRPQCW = new TH1D("fIntFlowResultsRPQCW","Integrated Flow (RP) from Q-cumulants with Weights",4,0,4);
 fIntFlowResultsRPQCW->SetLabelSize(0.06);
 fIntFlowResultsRPQCW->SetMarkerStyle(25);
 (fIntFlowResultsRPQCW->GetXaxis())->SetBinLabel(1,"v_{n}{2}");
 (fIntFlowResultsRPQCW->GetXaxis())->SetBinLabel(2,"v_{n}{4}");
 (fIntFlowResultsRPQCW->GetXaxis())->SetBinLabel(3,"v_{n}{6}");
 (fIntFlowResultsRPQCW->GetXaxis())->SetBinLabel(4,"v_{n}{8}");
 fResultsList->Add(fIntFlowResultsRPQCW);
 
 // <cos n(psi1-phi2)> for POIs:
 f2pPtEtaPOI = new TProfile2D("f2pPtEtaPOI","<cos n(#psi_{1}-#phi_{2})> (p_{t},#eta) for POIs",
                               fnBinsPt,fPtMin,fPtMax,fnBinsEta,fEtaMin,fEtaMax,"s");
 f2pPtEtaPOI->SetXTitle("p_{t}");
 f2pPtEtaPOI->SetYTitle("#eta");
 fDiffFlowList->Add(f2pPtEtaPOI);
 
 // <cos n(psi1+phi2-phi3-phi4)> for POIs:
 f4pPtEtaPOI = new TProfile2D("f4pPtEtaPOI","<cos n(#psi_{1}+#phi_{2}-#phi_{3}-#phi_{4})> (p_{t},#eta) for POIs",
                               fnBinsPt,fPtMin,fPtMax,fnBinsEta,fEtaMin,fEtaMax,"s");
 f4pPtEtaPOI->SetXTitle("p_{t}");
 f4pPtEtaPOI->SetYTitle("#eta");
 fDiffFlowList->Add(f4pPtEtaPOI);
 
 // <cos n(psi1+phi2+phi3-phi4-phi5-phi6)> for POIs:
 f6pPtEtaPOI = new TProfile2D("f6pPtEtaPOI","<cos n(#psi_{1}+#phi_{2}+#phi_{3}-#phi_{4}-#phi_{5}-#phi_{6})> (p_{t},#eta) for POIs",
                               fnBinsPt,fPtMin,fPtMax,fnBinsEta,fEtaMin,fEtaMax,"s");
 f6pPtEtaPOI->SetXTitle("p_{t}");
 f6pPtEtaPOI->SetYTitle("#eta");
 fDiffFlowList->Add(f6pPtEtaPOI);
 
 // <cos n(psi1+phi2+phi3+phi4-phi5-phi6-phi7-phi8)> for POIs:
 f8pPtEtaPOI = new TProfile2D("f8pPtEtaPOI","<cos n(#psi_{1}+#phi_{2}+#phi_{3}+#phi_{4}-#phi_{5}-#phi_{6}-#phi_{7}-#phi_{8})> (p_{t},#eta) for POIs",
                               fnBinsPt,fPtMin,fPtMax,fnBinsEta,fEtaMin,fEtaMax,"s");
 f8pPtEtaPOI->SetXTitle("p_{t}");
 f8pPtEtaPOI->SetYTitle("#eta");
 fDiffFlowList->Add(f8pPtEtaPOI);
 
 // non-weighted v'_{n}{2,QC} (pt,eta) for POIs
 fvn2ndPtEtaPOI = new TH2D("fvn2ndPtEtaPOI","v'_{n}{2,QC} (p_{t},#eta) for POIs",fnBinsPt,fPtMin,fPtMax,fnBinsEta,fEtaMin,fEtaMax);
 fvn2ndPtEtaPOI->SetXTitle("p_{t}");
 fvn2ndPtEtaPOI->SetYTitle("#eta");
 fResultsList->Add(fvn2ndPtEtaPOI);
 
 // non-weighted v'_{n}{4,QC} (pt,eta) for POIs
 fvn4thPtEtaPOI = new TH2D("fvn4thPtEtaPOI","v'_{n}{4,QC} (p_{t},#eta) for POIs",fnBinsPt,fPtMin,fPtMax,fnBinsEta,fEtaMin,fEtaMax);
 fvn4thPtEtaPOI->SetXTitle("p_{t}");
 fvn4thPtEtaPOI->SetYTitle("#eta");
 fResultsList->Add(fvn4thPtEtaPOI);

 // non-weighted v'_{n}{6,QC} (pt,eta) for POIs
 fvn6thPtEtaPOI = new TH2D("fvn6thPtEtaPOI","v'_{n}{6,QC} (p_{t},#eta) for POIs",fnBinsPt,fPtMin,fPtMax,fnBinsEta,fEtaMin,fEtaMax);
 fvn6thPtEtaPOI->SetXTitle("p_{t}");
 fvn6thPtEtaPOI->SetYTitle("#eta");
 fResultsList->Add(fvn6thPtEtaPOI);

 // non-weighted v'_{n}{8,QC} (pt,eta) for POIs
 fvn8thPtEtaPOI = new TH2D("fvn8thPtEtaPOI","v'_{n}{8,QC} (p_{t},#eta) for POIs",fnBinsPt,fPtMin,fPtMax,fnBinsEta,fEtaMin,fEtaMax);
 fvn8thPtEtaPOI->SetXTitle("p_{t}");
 fvn8thPtEtaPOI->SetYTitle("#eta");
 fResultsList->Add(fvn8thPtEtaPOI);
 
 // non-weighted v'_{n}{2,QC} (pt) for POIs
 fvn2ndPtPOI = new TH1D("fvn2ndPtPOI","v'_{n}{2,QC} (p_{t}) for POIs",fnBinsPt,fPtMin,fPtMax);
 fvn2ndPtPOI->SetXTitle("p_{t}");
 fResultsList->Add(fvn2ndPtPOI);
 
 // non-weighted v'_{n}{4,QC} (pt) for POIs
 fvn4thPtPOI = new TH1D("fvn4thPtPOI","v'_{n}{4,QC} (p_{t}) for POIs",fnBinsPt,fPtMin,fPtMax);
 fvn4thPtPOI->SetXTitle("p_{t}");
 fvn4thPtPOI->SetYTitle("#eta");
 fResultsList->Add(fvn4thPtPOI);
 
 // non-weighted v'_{n}{6,QC} (pt) for POIs
 fvn6thPtPOI = new TH1D("fvn6thPtPOI","v'_{n}{6,QC} (p_{t}) for POIs",fnBinsPt,fPtMin,fPtMax);
 fvn6thPtPOI->SetXTitle("p_{t}");
 fResultsList->Add(fvn6thPtPOI);
 
 // non-weighted v'_{n}{8,QC} (pt) for POIs
 fvn8thPtPOI = new TH1D("fvn8thPtPOI","v'_{n}{8,QC} (p_{t}) for POIs",fnBinsPt,fPtMin,fPtMax);
 fvn8thPtPOI->SetXTitle("p_{t}");
 fResultsList->Add(fvn8thPtPOI);
 
 // non-weighted v'_{n}{2,QC} (eta) for POIs
 fvn2ndEtaPOI = new TH1D("fvn2ndEtaPOI","v'_{n}{2,QC} (#eta) for POIs",fnBinsEta,fEtaMin,fEtaMax);
 fvn2ndEtaPOI->SetXTitle("#eta");
 fResultsList->Add(fvn2ndEtaPOI);
 
 // non-weighted v'_{n}{4,QC} (eta) for POIs
 fvn4thEtaPOI = new TH1D("fvn4thEtaPOI","v'_{n}{4,QC} (#eta) for POIs",fnBinsEta,fEtaMin,fEtaMax);
 fvn4thEtaPOI->SetXTitle("#eta");
 fResultsList->Add(fvn4thEtaPOI);
 
 // non-weighted v'_{n}{6,QC} (eta) for POIs
 fvn6thEtaPOI = new TH1D("fvn6thEtaPOI","v'_{n}{6,QC} (#eta) for POIs",fnBinsEta,fEtaMin,fEtaMax);
 fvn6thEtaPOI->SetXTitle("#eta");
 fResultsList->Add(fvn6thEtaPOI);
 
 // non-weighted v'_{n}{8,QC} (eta) for POIs
 fvn8thEtaPOI = new TH1D("fvn8thEtaPOI","v'_{n}{8,QC} (#eta) for POIs",fnBinsEta,fEtaMin,fEtaMax);
 fvn8thEtaPOI->SetXTitle("p_{t}");
 fResultsList->Add(fvn8thEtaPOI);
 
 // <w2 cos n(psi1-phi2)> for POIs:
 f2pPtEtaPOIW = new TProfile2D("f2pPtEtaPOIW","<w_{2} cos n(#psi_{1}-#phi_{2})> (p_{t},#eta) for POIs",
                               fnBinsPt,fPtMin,fPtMax,fnBinsEta,fEtaMin,fEtaMax,"s");
 f2pPtEtaPOIW->SetXTitle("p_{t}");
 fDiffFlowList->Add(f2pPtEtaPOIW);
 
 // <w2 w3 w4 cos n(psi1+phi2-phi3-phi4)> for POIs:
 f4pPtEtaPOIW = new TProfile2D("f4pPtEtaPOIW","<w_{2}w_{3}w_{4} cos n(#psi_{1}+#phi_{2}-#phi_{3}-#phi_{4})> (p_{t},#eta) for POIs",
                               fnBinsPt,fPtMin,fPtMax,fnBinsEta,fEtaMin,fEtaMax,"s");
 f4pPtEtaPOIW->SetXTitle("p_{t}");
 fDiffFlowList->Add(f4pPtEtaPOIW);
 
 // <w2 w3 w4 w5 w6 cos n(psi1+phi2+phi3-phi4-phi5-phi6)> for POIs:
 f6pPtEtaPOIW = new TProfile2D("f6pPtEtaPOIW","<w_{2}w_{3}w_{4}w_{5}w_{6} cos n(#psi_{1}+#phi_{2}+#phi_{3}-#phi_{4}-#phi_{5}-#phi_{6})> (p_{t},#eta) for POIs",fnBinsPt,fPtMin,fPtMax,fnBinsEta,fEtaMin,fEtaMax,"s");
 f6pPtEtaPOIW->SetXTitle("p_{t}");
 fDiffFlowList->Add(f6pPtEtaPOIW);
 
 // <w2 w3 w4 w5 w6 w7 w8 cos n(psi1+phi2+phi3+phi4-phi5-phi6-phi7-phi8)> for POIs:
 f8pPtEtaPOIW = new TProfile2D("f8pPtEtaPOIW","<w_{2}w_{3}w_{4}w_{5}w_{6}w_{7}w_{8} cos n(#psi_{1}+#phi_{2}+#phi_{3}+#phi_{4}-#phi_{5}-#phi_{6}-#phi_{7}-#phi_{8})> (p_{t},#eta) for POIs",fnBinsPt,fPtMin,fPtMax,fnBinsEta,fEtaMin,fEtaMax,"s");
 f8pPtEtaPOIW->SetXTitle("p_{t}");
 f8pPtEtaPOIW->SetYTitle("#eta");
 fDiffFlowList->Add(f8pPtEtaPOIW);

 // weighted v'_{n}{2,QC} (pt,eta) for POIs
 fvn2ndPtEtaPOIW = new TH2D("fvn2ndPtEtaPOIW","weighted v'_{n}{2,QC} (p_{t},#eta) for POIs",fnBinsPt,fPtMin,fPtMax,fnBinsEta,fEtaMin,fEtaMax);
 fvn2ndPtEtaPOIW->SetXTitle("p_{t}");
 fvn2ndPtEtaPOIW->SetYTitle("#eta");
 fResultsList->Add(fvn2ndPtEtaPOIW);
 
 // weighted v'_{n}{4,QC} (pt,eta) for POIs
 fvn4thPtEtaPOIW = new TH2D("fvn4thPtEtaPOIW","weighted v'_{n}{4,QC} (p_{t},#eta) for POIs",fnBinsPt,fPtMin,fPtMax,fnBinsEta,fEtaMin,fEtaMax);
 fvn4thPtEtaPOIW->SetXTitle("p_{t}");
 fvn4thPtEtaPOIW->SetYTitle("#eta");
 fResultsList->Add(fvn4thPtEtaPOIW);

 // weighted v'_{n}{6,QC} (pt,eta) for POIs
 fvn6thPtEtaPOIW = new TH2D("fvn6thPtEtaPOIW","weighted v'_{n}{6,QC} (p_{t},#eta) for POIs",fnBinsPt,fPtMin,fPtMax,fnBinsEta,fEtaMin,fEtaMax);
 fvn6thPtEtaPOIW->SetXTitle("p_{t}");
 fvn6thPtEtaPOIW->SetYTitle("#eta");
 fResultsList->Add(fvn6thPtEtaPOIW);

 // weighted v'_{n}{8,QC} (pt,eta) for POIs
 fvn8thPtEtaPOIW = new TH2D("fvn8thPtEtaPOIW","weighted v'_{n}{8,QC} (p_{t},#eta) for POIs",fnBinsPt,fPtMin,fPtMax,fnBinsEta,fEtaMin,fEtaMax);
 fvn8thPtEtaPOIW->SetXTitle("p_{t}");
 fvn8thPtEtaPOIW->SetYTitle("#eta");
 fResultsList->Add(fvn8thPtEtaPOIW);
 
  // weighted v'_{n}{2,QC} (pt) for POIs
 fvn2ndPtPOIW = new TH1D("fvn2ndPtPOIW","weighted v'_{n}{2,QC} (p_{t}) for POIs",fnBinsPt,fPtMin,fPtMax);
 fvn2ndPtPOIW->SetXTitle("p_{t}");
 fResultsList->Add(fvn2ndPtPOIW);
 
 // weighted v'_{n}{4,QC} (pt) for POIs
 fvn4thPtPOIW = new TH1D("fvn4thPtPOIW","weighted v'_{n}{4,QC} (p_{t}) for POIs",fnBinsPt,fPtMin,fPtMax);
 fvn4thPtPOIW->SetXTitle("p_{t}");
 fResultsList->Add(fvn4thPtPOIW);
 
 // weighted v'_{n}{6,QC} (pt) for POIs
 fvn6thPtPOIW = new TH1D("fvn6thPtPOIW","weighted v'_{n}{6,QC} (p_{t}) for POIs",fnBinsPt,fPtMin,fPtMax);
 fvn6thPtPOIW->SetXTitle("p_{t}");
 fResultsList->Add(fvn6thPtPOIW);
 
 // weighted v'_{n}{8,QC} (pt) for POIs
 fvn8thPtPOIW = new TH1D("fvn8thPtPOIW","weighted v'_{n}{8,QC} (p_{t}) for POIs",fnBinsPt,fPtMin,fPtMax);
 fvn8thPtPOIW->SetXTitle("p_{t}");
 fResultsList->Add(fvn8thPtPOIW);
 
 // weighted v'_{n}{2,QC} (eta) for POIs
 fvn2ndEtaPOIW = new TH1D("fvn2ndEtaPOIW","weighted v'_{n}{2,QC} (#eta) for POIs",fnBinsEta,fEtaMin,fEtaMax);
 fvn2ndEtaPOIW->SetXTitle("#eta");
 fResultsList->Add(fvn2ndEtaPOIW);
 
 // weighted v'_{n}{4,QC} (eta) for POIs
 fvn4thEtaPOIW = new TH1D("fvn4thEtaPOIW","weighted v'_{n}{4,QC} (#eta) for POIs",fnBinsEta,fEtaMin,fEtaMax);
 fvn4thEtaPOIW->SetXTitle("#eta");
 fResultsList->Add(fvn4thEtaPOIW);
 
 // weighted v'_{n}{6,QC} (eta) for POIs
 fvn6thEtaPOIW = new TH1D("fvn6thEtaPOIW","weighted v'_{n}{6,QC} (#eta) for POIs",fnBinsEta,fEtaMin,fEtaMax);
 fvn6thEtaPOIW->SetXTitle("#eta");
 fResultsList->Add(fvn6thEtaPOIW);
 
 // weighted v'_{n}{8,QC} (eta) for POIs
 fvn8thEtaPOIW = new TH1D("fvn8thEtaPOIW","weighted v'_{n}{8,QC} (#eta) for POIs",fnBinsEta,fEtaMin,fEtaMax);
 fvn8thEtaPOIW->SetXTitle("#eta");
 fResultsList->Add(fvn8thEtaPOIW);
 
 // <cos n(psi1-phi2)> for RPs:
 f2pPtEtaRP = new TProfile2D("f2pPtEtaRP","<cos n(#psi_{1}-#phi_{2})> (p_{t},#eta) for RPs",
                               fnBinsPt,fPtMin,fPtMax,fnBinsEta,fEtaMin,fEtaMax,"s");
 f2pPtEtaRP->SetXTitle("p_{t}");
 f2pPtEtaRP->SetYTitle("#eta");
 fDiffFlowList->Add(f2pPtEtaRP);
 
 // <cos n(psi1+phi2-phi3-phi4)> for RPs:
 f4pPtEtaRP = new TProfile2D("f4pPtEtaRP","<cos n(#psi_{1}+#phi_{2}-#phi_{3}-#phi_{4})> (p_{t},#eta) for RPs",
                               fnBinsPt,fPtMin,fPtMax,fnBinsEta,fEtaMin,fEtaMax,"s");
 f4pPtEtaRP->SetXTitle("p_{t}");
 f4pPtEtaRP->SetYTitle("#eta");
 fDiffFlowList->Add(f4pPtEtaRP);
 
 // <cos n(psi1+phi2+phi3-phi4-phi5-phi6)> for RPs:
 f6pPtEtaRP = new TProfile2D("f6pPtEtaRP","<cos n(#psi_{1}+#phi_{2}+#phi_{3}-#phi_{4}-#phi_{5}-#phi_{6})> (p_{t},#eta) for RPs",
                               fnBinsPt,fPtMin,fPtMax,fnBinsEta,fEtaMin,fEtaMax,"s");
 f6pPtEtaRP->SetXTitle("p_{t}");
 f6pPtEtaRP->SetYTitle("#eta");
 fDiffFlowList->Add(f6pPtEtaRP);
 
 // <cos n(psi1+phi2+phi3+phi4-phi5-phi6-phi7-phi8)> for RPs:
 f8pPtEtaRP = new TProfile2D("f8pPtEtaRP","<cos n(#psi_{1}+#phi_{2}+#phi_{3}+#phi_{4}-#phi_{5}-#phi_{6}-#phi_{7}-#phi_{8})> (p_{t},#eta) for RPs",
                               fnBinsPt,fPtMin,fPtMax,fnBinsEta,fEtaMin,fEtaMax,"s");
 f8pPtEtaRP->SetXTitle("p_{t}");
 f8pPtEtaRP->SetYTitle("#eta");
 fDiffFlowList->Add(f8pPtEtaRP);
 
 // non-weighted v'_{n}{2,QC} (pt,eta) for RPs
 fvn2ndPtEtaRP = new TH2D("fvn2ndPtEtaRP","v'_{n}{2,QC} (p_{t},#eta) for RPs",fnBinsPt,fPtMin,fPtMax,fnBinsEta,fEtaMin,fEtaMax);
 fvn2ndPtEtaRP->SetXTitle("p_{t}");
 fvn2ndPtEtaRP->SetYTitle("#eta");
 fResultsList->Add(fvn2ndPtEtaRP);
 
 // non-weighted v'_{n}{4,QC} (pt,eta) for RPs
 fvn4thPtEtaRP = new TH2D("fvn4thPtEtaRP","v'_{n}{4,QC} (p_{t},#eta) for RPs",fnBinsPt,fPtMin,fPtMax,fnBinsEta,fEtaMin,fEtaMax);
 fvn4thPtEtaRP->SetXTitle("p_{t}");
 fvn4thPtEtaRP->SetYTitle("#eta");
 fResultsList->Add(fvn4thPtEtaRP);

 // non-weighted v'_{n}{6,QC} (pt,eta) for RPs
 fvn6thPtEtaRP = new TH2D("fvn6thPtEtaRP","v'_{n}{6,QC} (p_{t},#eta) for RPs",fnBinsPt,fPtMin,fPtMax,fnBinsEta,fEtaMin,fEtaMax);
 fvn6thPtEtaRP->SetXTitle("p_{t}");
 fvn6thPtEtaRP->SetYTitle("#eta");
 fResultsList->Add(fvn6thPtEtaRP);

 // non-weighted v'_{n}{8,QC} (pt,eta) for RPs
 fvn8thPtEtaRP = new TH2D("fvn8thPtEtaRP","v'_{n}{8,QC} (p_{t},#eta) for RPs",fnBinsPt,fPtMin,fPtMax,fnBinsEta,fEtaMin,fEtaMax);
 fvn8thPtEtaRP->SetXTitle("p_{t}");
 fvn8thPtEtaRP->SetYTitle("#eta");
 fResultsList->Add(fvn8thPtEtaRP);
 
 // non-weighted v'_{n}{2,QC} (pt) for RPs
 fvn2ndPtRP = new TH1D("fvn2ndPtRP","v'_{n}{2,QC} (p_{t}) for RPs",fnBinsPt,fPtMin,fPtMax);
 fvn2ndPtRP->SetXTitle("p_{t}");
 fResultsList->Add(fvn2ndPtRP);
 
 // non-weighted v'_{n}{4,QC} (pt) for RPs
 fvn4thPtRP = new TH1D("fvn4thPtRP","v'_{n}{4,QC} (p_{t}) for RPs",fnBinsPt,fPtMin,fPtMax);
 fvn4thPtRP->SetXTitle("p_{t}");
 fResultsList->Add(fvn4thPtRP);

 // non-weighted v'_{n}{6,QC} (pt) for RPs
 fvn6thPtRP = new TH1D("fvn6thPtRP","v'_{n}{6,QC} (p_{t}) for RPs",fnBinsPt,fPtMin,fPtMax);
 fvn6thPtRP->SetXTitle("p_{t}");
 fResultsList->Add(fvn6thPtRP);

 // non-weighted v'_{n}{8,QC} (pt) for RPs
 fvn8thPtRP = new TH1D("fvn8thPtRP","v'_{n}{8,QC} (p_{t}) for RPs",fnBinsPt,fPtMin,fPtMax);
 fvn8thPtRP->SetXTitle("p_{t}");
 fResultsList->Add(fvn8thPtRP);

 // non-weighted v'_{n}{2,QC} (eta) for RPs
 fvn2ndEtaRP = new TH1D("fvn2ndEtaRP","v'_{n}{2,QC} (#eta) for RPs",fnBinsEta,fEtaMin,fEtaMax);
 fvn2ndEtaRP->SetXTitle("#eta");
 fResultsList->Add(fvn2ndEtaRP);
 
 // non-weighted v'_{n}{4,QC} (eta) for RPs
 fvn4thEtaRP = new TH1D("fvn4thEtaRP","v'_{n}{4,QC} (#eta) for RPs",fnBinsEta,fEtaMin,fEtaMax);
 fvn4thEtaRP->SetXTitle("#eta");
 fResultsList->Add(fvn4thEtaRP);

 // non-weighted v'_{n}{6,QC} (eta) for RPs
 fvn6thEtaRP = new TH1D("fvn6thEtaRP","v'_{n}{6,QC} (#eta) for RPs",fnBinsEta,fEtaMin,fEtaMax);
 fvn6thEtaRP->SetXTitle("#eta");
 fResultsList->Add(fvn6thEtaRP);

 // non-weighted v'_{n}{8,QC} (eta) for RPs
 fvn8thEtaRP = new TH1D("fvn8thEtaRP","v'_{n}{8,QC} (#eta) for RPs",fnBinsEta,fEtaMin,fEtaMax);
 fvn8thEtaRP->SetXTitle("#eta");
 fResultsList->Add(fvn8thEtaRP);
 
 // <w2 cos n(psi1-phi2)> for RPs:
 f2pPtEtaRPW = new TProfile2D("f2pPtEtaRPW","<w_{2} cos n(#psi_{1}-#phi_{2})> (p_{t},#eta) for RPs",
                               fnBinsPt,fPtMin,fPtMax,fnBinsEta,fEtaMin,fEtaMax,"s");
 f2pPtEtaRPW->SetXTitle("p_{t}");
 f2pPtEtaRPW->SetYTitle("#eta");
 fDiffFlowList->Add(f2pPtEtaRPW);
 
 // <w2 w3 w4 cos n(psi1+phi2-phi3-phi4)> for RPs:
 f4pPtEtaRPW = new TProfile2D("f4pPtEtaRPW","<w_{2}w_{3}w_{4} cos n(#psi_{1}+#phi_{2}-#phi_{3}-#phi_{4})> (p_{t},#eta) for RPs",
                               fnBinsPt,fPtMin,fPtMax,fnBinsEta,fEtaMin,fEtaMax,"s");
 f4pPtEtaRPW->SetXTitle("p_{t}");
 f4pPtEtaRPW->SetYTitle("#eta");
 fDiffFlowList->Add(f4pPtEtaRPW);
 
 // <w2 w3 w4 w5 w6 cos n(psi1+phi2+phi3-phi4-phi5-phi6)> for RPs:
 f6pPtEtaRPW = new TProfile2D("f6pPtEtaRPW","<w_{2}w_{3}w_{4}w_{5}w_{6} cos n(#psi_{1}+#phi_{2}+#phi_{3}-#phi_{4}-#phi_{5}-#phi_{6})> (p_{t},#eta) for RPs",fnBinsPt,fPtMin,fPtMax,fnBinsEta,fEtaMin,fEtaMax,"s");
 f6pPtEtaRPW->SetXTitle("p_{t}");
 f6pPtEtaRPW->SetYTitle("#eta");
 fDiffFlowList->Add(f6pPtEtaRPW);
 
 // <w2 w3 w4 w5 w6 w7 w8 cos n(psi1+phi2+phi3+phi4-phi5-phi6-phi7-phi8)> for RPs:
 f8pPtEtaRPW = new TProfile2D("f8pPtEtaRPW","<w_{2}w_{3}w_{4}w_{5}w_{6}w_{7}w_{8} cos n(#psi_{1}+#phi_{2}+#phi_{3}+#phi_{4}-#phi_{5}-#phi_{6}-#phi_{7}-#phi_{8})> (p_{t},#eta) for RPs",fnBinsPt,fPtMin,fPtMax,fnBinsEta,fEtaMin,fEtaMax,"s");
 f8pPtEtaRPW->SetXTitle("p_{t}");
 f8pPtEtaRPW->SetYTitle("#eta");
 fDiffFlowList->Add(f8pPtEtaRPW);
 
 // weighted v'_{n}{2,QC} (pt,eta) for RPs
 fvn2ndPtEtaRPW = new TH2D("fvn2ndPtEtaRPW","weighted v'_{n}{2,QC} (p_{t},#eta) for RPs",fnBinsPt,fPtMin,fPtMax,fnBinsEta,fEtaMin,fEtaMax);
 fvn2ndPtEtaRPW->SetXTitle("p_{t}");
 fvn2ndPtEtaRPW->SetYTitle("#eta");
 fResultsList->Add(fvn2ndPtEtaRPW);
 
 // weighted v'_{n}{4,QC} (pt,eta) for RPs
 fvn4thPtEtaRPW = new TH2D("fvn4thPtEtaRPW","weighted v'_{n}{4,QC} (p_{t},#eta) for RPs",fnBinsPt,fPtMin,fPtMax,fnBinsEta,fEtaMin,fEtaMax);
 fvn4thPtEtaRPW->SetXTitle("p_{t}");
 fvn4thPtEtaRPW->SetYTitle("#eta");
 fResultsList->Add(fvn4thPtEtaRPW);

 // weighted v'_{n}{6,QC} (pt,eta) for RPs
 fvn6thPtEtaRPW = new TH2D("fvn6thPtEtaRPW","weighted v'_{n}{6,QC} (p_{t},#eta) for RPs",fnBinsPt,fPtMin,fPtMax,fnBinsEta,fEtaMin,fEtaMax);
 fvn6thPtEtaRPW->SetXTitle("p_{t}");
 fvn6thPtEtaRPW->SetYTitle("#eta");
 fResultsList->Add(fvn6thPtEtaRPW);

 // weighted v'_{n}{8,QC} (pt,eta) for RPs
 fvn8thPtEtaRPW = new TH2D("fvn8thPtEtaRPW","weighted v'_{n}{8,QC} (p_{t},#eta) for RPs",fnBinsPt,fPtMin,fPtMax,fnBinsEta,fEtaMin,fEtaMax);
 fvn8thPtEtaRPW->SetXTitle("p_{t}");
 fvn8thPtEtaRPW->SetYTitle("#eta");
 fResultsList->Add(fvn8thPtEtaRPW);
 
 // weighted v'_{n}{2,QC} (pt) for RPs
 fvn2ndPtRPW = new TH1D("fvn2ndPtRPW","weighted v'_{n}{2,QC} (p_{t}) for RPs",fnBinsPt,fPtMin,fPtMax);
 fvn2ndPtRPW->SetXTitle("p_{t}");
 fResultsList->Add(fvn2ndPtRPW);
 
 // weighted v'_{n}{4,QC} (pt) for RPs
 fvn4thPtRPW = new TH1D("fvn4thPtRPW","weighted v'_{n}{4,QC} (p_{t}) for RPs",fnBinsPt,fPtMin,fPtMax);
 fvn4thPtRPW->SetXTitle("p_{t}");
 fResultsList->Add(fvn4thPtRPW);

 // weighted v'_{n}{6,QC} (pt) for RPs
 fvn6thPtRPW = new TH1D("fvn6thPtRPW","weighted v'_{n}{6,QC} (p_{t}) for RPs",fnBinsPt,fPtMin,fPtMax);
 fvn6thPtRPW->SetXTitle("p_{t}");
 fResultsList->Add(fvn6thPtRPW);

 // weighted v'_{n}{8,QC} (pt) for RPs
 fvn8thPtRPW = new TH1D("fvn8thPtRPW","weighted v'_{n}{8,QC} (p_{t}) for RPs",fnBinsPt,fPtMin,fPtMax);
 fvn8thPtRPW->SetXTitle("p_{t}");
 fResultsList->Add(fvn8thPtRPW);

 // weighted v'_{n}{2,QC} (eta) for RPs
 fvn2ndEtaRPW = new TH1D("fvn2ndEtaRPW","weighted v'_{n}{2,QC} (#eta) for RPs",fnBinsEta,fEtaMin,fEtaMax);
 fvn2ndEtaRPW->SetXTitle("#eta");
 fResultsList->Add(fvn2ndEtaRPW);
 
 // weighted v'_{n}{4,QC} (eta) for RPs
 fvn4thEtaRPW = new TH1D("fvn4thEtaRPW","weighted v'_{n}{4,QC} (#eta) for RPs",fnBinsEta,fEtaMin,fEtaMax);
 fvn4thEtaRPW->SetXTitle("#eta");
 fResultsList->Add(fvn4thEtaRPW);

 // weighted v'_{n}{6,QC} (eta) for RPs
 fvn6thEtaRPW = new TH1D("fvn6thEtaRPW","weighted v'_{n}{6,QC} (#eta) for RPs",fnBinsEta,fEtaMin,fEtaMax);
 fvn6thEtaRPW->SetXTitle("#eta");
 fResultsList->Add(fvn6thEtaRPW);

 // weighted v'_{n}{8,QC} (eta) for RPs
 fvn8thEtaRP = new TH1D("fvn8thEtaEtaRP","weighted v'_{n}{8,QC} (#eta) for RPs",fnBinsEta,fEtaMin,fEtaMax);
 fvn8thEtaRP->SetXTitle("#eta");
 fResultsList->Add(fvn8thEtaRP);
 // .....................................................................................................................................
 
 
 
 
 // add fUseWeightsBits to the main list (to be improved)
 fUseWeightsBits = new TBits(1);
 fHistList->Add(fUseWeightsBits);
 
 // add list fWeightsList with weights to the main list
 fHistList->Add(fWeightsList);
  
 // add list fDiffFlowList with histograms and profiles needed for differential flow to the main list 
 fHistList->Add(fDiffFlowList); 
 
 // add list fResultsList with final results to the main list 
 fHistList->Add(fResultsList); 

 
}//end of Init()


//================================================================================================================


void AliFlowAnalysisWithQCumulants::Make(AliFlowEventSimple* anEvent)
{
 // running over data only in this method
 
 
 
 
 //                                     *********************************************
 //                                     **** ACCESS THE OUTPUT FILE WITH WEIGHTS ****
 //                                     ********************************************* 
 
 fUseWeights = fUsePhiWeights||fUsePtWeights||fUseEtaWeights;
 fUseWeightsBits->SetBitNumber(1,fUseWeights); // to be improved (how to pass boolean to Finish()?)
 
 TH1F *phiWeights = NULL; // histogram with phi weights
 TH1D *ptWeights  = NULL; // histogram with pt weights
 TH1D *etaWeights = NULL; // histogram with eta weights
   
 if(fUseWeights)
 {
  if(!fWeightsList)
  {
   cout<<" WARNING: fWeightsList is NULL pointer in AFAWQC::Make(). "<<endl;
   exit(0);
  }
  if(fUsePhiWeights) 
  {
   phiWeights = dynamic_cast<TH1F *>(fWeightsList->FindObject("phi_weights"));
   if(!phiWeights)
   {
    cout<<" WARNING: couldn't access the histogram with phi weights in AFAWQC::Make(). "<<endl;
    exit(0);
   } 
  } 
  if(fUsePtWeights) 
  { 
   ptWeights = dynamic_cast<TH1D *>(fWeightsList->FindObject("pt_weights"));
   if(!ptWeights) 
   {
    cout<<" WARNING: couldn't access the histogram with pt weights in AFAWQC::Make(). "<<endl;
    exit(0);
   } 
  } 
  if(fUseEtaWeights) 
  {
   etaWeights = dynamic_cast<TH1D *>(fWeightsList->FindObject("eta_weights"));
   if(!etaWeights) 
   {
    cout<<" WARNING: couldn't access the histogram with eta weights in AFAWQC::Make(). "<<endl;
    exit(0);
   }
  } 
 } 
  
 Int_t nBinsPhi = 0; 
 Double_t dBinWidthPt = 0.;
 Double_t dBinWidthEta = 0.;
 
 if(fnBinsPt)
 {
  dBinWidthPt=(fPtMax-fPtMin)/fnBinsPt;  
 } 
 
 if(fnBinsEta)
 {
  dBinWidthEta=(fEtaMax-fEtaMin)/fnBinsEta;  
 } 
 
 if(fWeightsList)
 {
  if(fUsePhiWeights)
  {
   if(phiWeights) nBinsPhi = phiWeights->GetNbinsX();
  }          
  if(fUsePtWeights)
  {
   if(ptWeights)
   {
    Double_t dBinWidthPtW = ptWeights->GetBinWidth(1); // assuming that all bins have the same width
    if(dBinWidthPtW != dBinWidthPt)
    {
     cout<<" WARNING: dBinWidthPtW != dBinWidthPt in AFAWQC::Make()."<<endl;
     exit(0);
    }
    Double_t dPtMinW = (ptWeights->GetXaxis())->GetXmin();
    if(dPtMinW != fPtMin)
    {
     cout<<" WARNING: dPtMinW != fPtMin in AFAWQC::Make()."<<endl;
     exit(0);
    }
   } 
  }       
  if(fUseEtaWeights)
  {
   if(etaWeights)
   {
    Double_t dBinWidthEtaW = etaWeights->GetBinWidth(1); // assuming that all bins have the same width
    if(dBinWidthEtaW != dBinWidthEta)
    {
     cout<<" WARNING: dBinWidthEtaW != dBinWidthEta in AFAWQC::Make()."<<endl;
     exit(0);
    }
    Double_t dEtaMinW = (etaWeights->GetXaxis())->GetXmin();
    if(dEtaMinW != fEtaMin)
    {
     cout<<" WARNING: dEtaMinW != fEtaMin in AFAWQC::Make()."<<endl;
     exit(0);
    }
   } 
  }          
 } // end of if(weightsList)
 
 Double_t dPhi = 0.; // azumithal angle in the laboratory frame
 Double_t dPt  = 0.; // transverse momentum
 Double_t dEta = 0.; // pseudorapidity

 Double_t wPhi = 1.; // phi weight
 Double_t wPt  = 1.; // pt weight
 Double_t wEta = 1.; // eta weight
 
                                 
                                                                 
                                                                                                                                 
 //                                     ********************************************
 //                                     **** FILL THE COMMON CONTROL HISTOGRAMS ****
 //                                     ********************************************
                                         
 Int_t nRP = anEvent->GetEventNSelTracksRP(); 
 if(nRP>1)
 {
  fCommonHists2nd->FillControlHistograms(anEvent);                                        
  if(nRP>3)
  {
   fCommonHists4th->FillControlHistograms(anEvent);                                        
   if(nRP>5)
   {
    fCommonHists6th->FillControlHistograms(anEvent);                                        
    if(nRP>7)
    {
     fCommonHists8th->FillControlHistograms(anEvent);                                        
    } // end of if(nRP>7)  
   } // end of if(nRP>5) 
  } // end of if(nRP>3)                                                                                                                      
 } // end of if(nRP>1) 
 
 
 
                                                                 
 //                                     ***************************
 //                                     **** LOOPING OVER DATA ****
 //                                     ***************************
 
 Int_t nPrim = anEvent->NumberOfTracks(); 
 
 // nPrim = total number of primary tracks, i.e. nPrim = nRP + nPOI + rest, where:
 // nRP   = # of particles used to determine the reaction plane;
 // nPOI  = # of particles of interest for a detailed flow analysis;
 // rest  = # of particles which are niether RPs not POIs.  
 
 for(Int_t i=0;i<nPrim;i++) 
 { 
  fTrack=anEvent->GetTrack(i);  
  if(fTrack)
  {
   if(!(fTrack->InRPSelection() || fTrack->InPOISelection())) continue;
 
   // checking the RP condition:
   if(fTrack->InRPSelection())
   {    
    dPhi = fTrack->Phi();
    dPt  = fTrack->Pt();
    dEta = fTrack->Eta();
  
    // determine phi weight for this particle: 
    if(phiWeights && nBinsPhi)
    {
     wPhi = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(dPhi*nBinsPhi/TMath::TwoPi())));
    }
    // determine pt weight for this particle:    
    if(ptWeights && dBinWidthPt)
    {
     wPt = ptWeights->GetBinContent(1+(Int_t)(TMath::Floor((dPt-fPtMin)/dBinWidthPt))); 
    }            
    // determine eta weight for this particle:    
    if(etaWeights && dBinWidthEta)
    {
     wEta = etaWeights->GetBinContent(1+(Int_t)(TMath::Floor((dEta-fEtaMin)/dBinWidthEta))); 
    } 

    // fill Re[Q_{n,k}] and Im[Q_{n,k}]:
    for(Int_t n=0;n<4;n++)
    {
     for(Int_t k=0;k<9;k++)
     {
      (*fReQ)(n,k)+=pow(wPhi*wPt*wEta,k)*TMath::Cos(2*(n+1)*dPhi);
      (*fImQ)(n,k)+=pow(wPhi*wPt*wEta,k)*TMath::Sin(2*(n+1)*dPhi);
     } 
    }
     
    // fill S^{M}_{p,k}:
    for(Int_t p=0;p<8;p++)
    {
     for(Int_t k=0;k<9;k++)
     {     
      (*fSMpk)(p,k)+=pow(wPhi*wPt*wEta,k);
     }
    } 
     
    Int_t n = 2; // to be improved (add setter for harmonic) 
     
    // fill non-weighted q_RPs 
    fReqRP1nPtEta->Fill(dPt,dEta,TMath::Cos(1.*n*dPhi));   
    fImqRP1nPtEta->Fill(dPt,dEta,TMath::Sin(1.*n*dPhi));
    fReqRP2nPtEta->Fill(dPt,dEta,TMath::Cos(2.*n*dPhi)); 
    fImqRP2nPtEta->Fill(dPt,dEta,TMath::Sin(2.*n*dPhi));
     
    // mRP:
    fmRPPtEta->Fill(dPt,dEta,1); 
    
    // fill weighted q_RPs 
    if(fUseWeights)
    {
     n = 2; // to be improved (add setter for harmonic) 
     
     // qRP_{n,k} (weighted qRP):
     fReqRP1n2kPtEta->Fill(dPt,dEta,pow(wPhi*wPt*wEta,2.)*TMath::Cos(1.*n*dPhi));   
     fImqRP1n2kPtEta->Fill(dPt,dEta,pow(wPhi*wPt*wEta,2.)*TMath::Sin(1.*n*dPhi));
     fReqRP2n1kPtEta->Fill(dPt,dEta,pow(wPhi*wPt*wEta,1.)*TMath::Cos(2.*n*dPhi)); 
     fImqRP2n1kPtEta->Fill(dPt,dEta,pow(wPhi*wPt*wEta,1.)*TMath::Sin(2.*n*dPhi)); 
  
     // S^{mRP}_{p,k}: 
     fSmRP1p1kPtEta->Fill(dPt,dEta,pow(wPhi*wPt*wEta,1.));
     fSmRP1p2kPtEta->Fill(dPt,dEta,pow(wPhi*wPt*wEta,2.));
     fSmRP1p3kPtEta->Fill(dPt,dEta,pow(wPhi*wPt*wEta,3.));     
    } 
    
    // checking if RP particle is also POI particle:
    if(fTrack->InPOISelection())
    { 
     n = 2; // to be improved (add setter for harmonic)  
     
     // q''_{n} (non-weighted q''):
     fReqPrimePrime1nPtEta->Fill(dPt,dEta,TMath::Cos(1.*n*dPhi));   
     fImqPrimePrime1nPtEta->Fill(dPt,dEta,TMath::Sin(1.*n*dPhi));
     fReqPrimePrime2nPtEta->Fill(dPt,dEta,TMath::Cos(2.*n*dPhi)); 
     fImqPrimePrime2nPtEta->Fill(dPt,dEta,TMath::Sin(2.*n*dPhi));
     
     // m'':
     fmPrimePrimePtEta->Fill(dPt,dEta,1); 

     if(fUseWeights)
     {
      // q''_{n,k} (weighted q''):
      fReqPrimePrime1n2kPtEta->Fill(dPt,dEta,pow(wPhi*wPt*wEta,2.)*TMath::Cos(1.*n*dPhi));   
      fImqPrimePrime1n2kPtEta->Fill(dPt,dEta,pow(wPhi*wPt*wEta,2.)*TMath::Sin(1.*n*dPhi));
      fReqPrimePrime2n1kPtEta->Fill(dPt,dEta,pow(wPhi*wPt*wEta,1.)*TMath::Cos(2.*n*dPhi)); 
      fImqPrimePrime2n1kPtEta->Fill(dPt,dEta,pow(wPhi*wPt*wEta,1.)*TMath::Sin(2.*n*dPhi)); 
  
      // S^{m''}_{p,k}: 
      fSmPrimePrime1p1kPtEta->Fill(dPt,dEta,pow(wPhi*wPt*wEta,1.));
      fSmPrimePrime1p2kPtEta->Fill(dPt,dEta,pow(wPhi*wPt*wEta,2.));
      fSmPrimePrime1p3kPtEta->Fill(dPt,dEta,pow(wPhi*wPt*wEta,3.));     
     }          
    } // end of if(fTrack->InPOISelection())
   } // end of if(pTrack->InRPSelection())
   
   // checking the POI condition:
   if(fTrack->InPOISelection())
   {
    Int_t n = 2; // to be improved (add setter for harmonic)  
    
    dPhi = fTrack->Phi();
    dPt  = fTrack->Pt();
    dEta = fTrack->Eta();
   
    // q_n:   
    fReqnPtEta->Fill(dPt,dEta,TMath::Cos(1.*n*dPhi));
    fImqnPtEta->Fill(dPt,dEta,TMath::Sin(1.*n*dPhi));
    
    // m: 
    fmPtEta->Fill(dPt,dEta,1);
      
   } // end of if(pTrack->InPOISelection() )  
  } // end of if(fTrack)
  else{
       cout<<endl;
       cout<<" WARNING: no particle! (i.e. fTrack is a NULL pointer in AFAWQC::Make().)"<<endl;
       cout<<endl;       
      }
 } // end of for(Int_t i=0;i<nPrim;i++) 
  
 // calculate the final expressions for S^{M}_{p,k} = (sum_{i=1}^{M} w_{i}^{k})^{p}:
 for(Int_t p=0;p<8;p++)
 {
  for(Int_t k=0;k<9;k++)
  {
   (*fSMpk)(p,k)=pow((*fSMpk)(p,k),p+1);
  }  
 } 
 
 
 
 
 //                                     *****************************
 //                                     **** CALLING THE METHODS ****
 //                                     *****************************

 // nested loops (needed for cross-checking the results):
 Bool_t evaluateNestedLoopsForIntegratedFlow = kFALSE;    // to be improved / removed
 Bool_t evaluateNestedLoopsForDifferentialFlow = kFALSE; // to be improved / removed
 
 if(evaluateNestedLoopsForIntegratedFlow && nPrim>0 && nPrim<14) // to be improved / removed (eventually I would not need this if())
 {
  // calculate all correlations needed for 'no-name' integrated flow WITHOUT weights 
  // (the results are stored in 1D profile fQCorrelations) 
  if(!(fUseWeights)) this->CalculateCorrelationsForIntegratedFlow();
 
  // calculate all correlations needed for 'no-name' integrated flow WITH weights 
  // (the results are stored in 1D profile fQCorrelationsW) 
  if(fUseWeights) this->CalculateWeightedCorrelationsForIntegratedFlow();
 }
 else if (!evaluateNestedLoopsForIntegratedFlow)
 {
  this->CalculateCorrelationsForIntegratedFlow();
  if(fUseWeights) this->CalculateWeightedCorrelationsForIntegratedFlow();
 }

 if(evaluateNestedLoopsForDifferentialFlow && nPrim>0 && nPrim<14 ) // to be improved / removed (eventually I would not need this if())
 {
  // calculate all correlations needed for differential flow WITHOUT weights 
  // and store the results in 2D profiles (pt,eta): 
  // a) POIs: f2pPtEtaPOI, f4pPtEtaPOI, f6pPtEtaPOI and f8pPtEtaPOI; 
  // b) RPs: f2pPtEtaRP, f4pPtEtaRP, f6pPtEtaRP and f8pPtEtaRP. 
  if(!(fUseWeights))
  {
   this->CalculateCorrelationsForDifferentialFlow("POI");
   this->CalculateCorrelationsForDifferentialFlow("RP");
  }
  // calculate all correlations needed for differential flow WITH weights
  // and store the results in 2D profiles (pt,eta): 
  // a) POIs: f2pPtEtaPOIW, f4pPtEtaPOIW, f6pPtEtaPOIW and f8pPtEtaPOIW; 
  // b) RPs: f2pPtEtaRPW, f4pPtEtaRPW, f6pPtEtaRPW and f8pPtEtaRPW.
  if(fUseWeights)
  { 
   this->CalculateWeightedCorrelationsForDifferentialFlow("POI");
   this->CalculateWeightedCorrelationsForDifferentialFlow("RP");
  } 
 }
 else if (!evaluateNestedLoopsForDifferentialFlow)
 {
  this->CalculateCorrelationsForDifferentialFlow("POI");
  this->CalculateCorrelationsForDifferentialFlow("RP");
  
  if(fUseWeights) 
  {
   this->CalculateWeightedCorrelationsForDifferentialFlow("POI");
   this->CalculateWeightedCorrelationsForDifferentialFlow("RP");
  }
   
 } 
 
 if(evaluateNestedLoopsForIntegratedFlow && nPrim>0 && nPrim<14) // to be improved / removed (eventually I would not need this if())
 {
  this->EvaluateNestedLoopsForIntegratedFlow(anEvent);  
 }
 
 if(evaluateNestedLoopsForDifferentialFlow && nPrim>0 && nPrim<14) // to be improved / removed (eventually I would not need this if())
 {
  this->EvaluateNestedLoopsForDifferentialFlow(anEvent);  
 }
 
 
 
 
 //                                     ********************************
 //                                     **** RESET E-B-E QUANTITIES ****
 //                                     ********************************
  
 fReQ->Zero();
 fImQ->Zero();
 fSMpk->Zero();
 fReqnPtEta->Reset();  
 fImqnPtEta->Reset(); 
 fmPtEta->Reset();  
 fReqPrimePrime1nPtEta->Reset();   
 fImqPrimePrime1nPtEta->Reset(); 
 fReqPrimePrime2nPtEta->Reset(); 
 fImqPrimePrime2nPtEta->Reset();     
 fmPrimePrimePtEta->Reset();  
 fReqPrimePrime1n2kPtEta->Reset();   
 fImqPrimePrime1n2kPtEta->Reset(); 
 fReqPrimePrime2n1kPtEta->Reset(); 
 fImqPrimePrime2n1kPtEta->Reset();    
 fSmPrimePrime1p1kPtEta->Reset(); 
 fSmPrimePrime1p2kPtEta->Reset(); 
 fSmPrimePrime1p3kPtEta->Reset(); 
 // qRPs (to be improved - notation)
 fReqRP1nPtEta->Reset();  
 fImqRP1nPtEta->Reset(); 
 fReqRP2nPtEta->Reset(); 
 fImqRP2nPtEta->Reset(); 
 fmRPPtEta->Reset(); 
 fReqRP1n2kPtEta->Reset();  
 fImqRP1n2kPtEta->Reset(); 
 fReqRP2n1kPtEta->Reset(); 
 fImqRP2n1kPtEta->Reset(); 
 fSmRP1p1kPtEta->Reset(); 
 fSmRP1p2kPtEta->Reset(); 
 fSmRP1p3kPtEta->Reset();  
 
} // end of AliFlowAnalysisWithQCumulants::Make(AliFlowEventSimple* anEvent)


//================================================================================================================================


void AliFlowAnalysisWithQCumulants::CalculateCorrelationsForIntegratedFlow()
{
 // calculate all correlations needed for 'no-name' integrated flow // to be improved (name)
 
 // multiplicity:
 Double_t dMult = (*fSMpk)(0,0);
 
 // real and imaginary parts of non-weighted Q-vectors evaluated in harmonics n, 2n, 3n and 4n: 
 Double_t dReQ1n = (*fReQ)(0,0);
 Double_t dReQ2n = (*fReQ)(1,0);
 Double_t dReQ3n = (*fReQ)(2,0);
 Double_t dReQ4n = (*fReQ)(3,0);
 Double_t dImQ1n = (*fImQ)(0,0);
 Double_t dImQ2n = (*fImQ)(1,0);
 Double_t dImQ3n = (*fImQ)(2,0);
 Double_t dImQ4n = (*fImQ)(3,0);
  
 // real and imaginary parts of some expressions involving various combinations of Q-vectors evaluated in harmonics n, 2n, 3n and 4n:
 // (these expression appear in the Eqs. for the multi-particle correlations bellow)
 
 // Re[Q_{2n} Q_{n}^* Q_{n}^*]
 Double_t reQ2nQ1nstarQ1nstar = pow(dReQ1n,2.)*dReQ2n + 2.*dReQ1n*dImQ1n*dImQ2n - pow(dImQ1n,2.)*dReQ2n; 
 
 // Im[Q_{2n} Q_{n}^* Q_{n}^*]
 //Double_t imQ2nQ1nstarQ1nstar = pow(dReQ1n,2.)*dImQ2n-2.*dReQ1n*dImQ1n*dReQ2n-pow(dImQ1n,2.)*dImQ2n; 
 
 // Re[Q_{n} Q_{n} Q_{2n}^*] = Re[Q_{2n} Q_{n}^* Q_{n}^*]
 Double_t reQ1nQ1nQ2nstar = reQ2nQ1nstarQ1nstar; 
 
 // Re[Q_{3n} Q_{n} Q_{2n}^* Q_{2n}^*]
 Double_t reQ3nQ1nQ2nstarQ2nstar = (pow(dReQ2n,2.)-pow(dImQ2n,2.))*(dReQ3n*dReQ1n-dImQ3n*dImQ1n) 
                                 + 2.*dReQ2n*dImQ2n*(dReQ3n*dImQ1n+dImQ3n*dReQ1n);

 // Im[Q_{3n} Q_{n} Q_{2n}^* Q_{2n}^*]                                                                  
 //Double_t imQ3nQ1nQ2nstarQ2nstar = calculate and implement this (deleteMe)
  
 // Re[Q_{2n} Q_{2n} Q_{3n}^* Q_{1n}^*] = Re[Q_{3n} Q_{n} Q_{2n}^* Q_{2n}^*]
 Double_t reQ2nQ2nQ3nstarQ1nstar = reQ3nQ1nQ2nstarQ2nstar;
  
 // Re[Q_{4n} Q_{2n}^* Q_{2n}^*]
 Double_t reQ4nQ2nstarQ2nstar = pow(dReQ2n,2.)*dReQ4n+2.*dReQ2n*dImQ2n*dImQ4n-pow(dImQ2n,2.)*dReQ4n;

 // Im[Q_{4n} Q_{2n}^* Q_{2n}^*]
 //Double_t imQ4nQ2nstarQ2nstar = calculate and implement this (deleteMe)
 
 // Re[Q_{2n} Q_{2n} Q_{4n}^*] =  Re[Q_{4n} Q_{2n}^* Q_{2n}^*]
 Double_t reQ2nQ2nQ4nstar = reQ4nQ2nstarQ2nstar;
 
 // Re[Q_{4n} Q_{3n}^* Q_{n}^*]
 Double_t reQ4nQ3nstarQ1nstar = dReQ4n*(dReQ3n*dReQ1n-dImQ3n*dImQ1n)+dImQ4n*(dReQ3n*dImQ1n+dImQ3n*dReQ1n);
 
 // Re[Q_{3n} Q_{n} Q_{4n}^*] = Re[Q_{4n} Q_{3n}^* Q_{n}^*]
 Double_t reQ3nQ1nQ4nstar = reQ4nQ3nstarQ1nstar;
 
 // Im[Q_{4n} Q_{3n}^* Q_{n}^*]
 //Double_t imQ4nQ3nstarQ1nstar = calculate and implement this (deleteMe)

 // Re[Q_{3n} Q_{2n}^* Q_{n}^*]
 Double_t reQ3nQ2nstarQ1nstar = dReQ3n*dReQ2n*dReQ1n-dReQ3n*dImQ2n*dImQ1n+dImQ3n*dReQ2n*dImQ1n
                              + dImQ3n*dImQ2n*dReQ1n;
                              
 // Re[Q_{2n} Q_{n} Q_{3n}^*] = Re[Q_{3n} Q_{2n}^* Q_{n}^*]
 Double_t reQ2nQ1nQ3nstar = reQ3nQ2nstarQ1nstar;
 
 // Im[Q_{3n} Q_{2n}^* Q_{n}^*]
 //Double_t imQ3nQ2nstarQ1nstar; //calculate and implement this (deleteMe)
 
 // Re[Q_{3n} Q_{n}^* Q_{n}^* Q_{n}^*]
 Double_t reQ3nQ1nstarQ1nstarQ1nstar = dReQ3n*pow(dReQ1n,3)-3.*dReQ1n*dReQ3n*pow(dImQ1n,2)
                                     + 3.*dImQ1n*dImQ3n*pow(dReQ1n,2)-dImQ3n*pow(dImQ1n,3);

 // Im[Q_{3n} Q_{n}^* Q_{n}^* Q_{n}^*]
 //Double_t imQ3nQ1nstarQ1nstarQ1nstar; //calculate and implement this (deleteMe)
 
 // |Q_{2n}|^2 |Q_{n}|^2
 Double_t dQ2nQ1nQ2nstarQ1nstar = (pow(dReQ2n,2.)+pow(dImQ2n,2.))*(pow(dReQ1n,2.)+pow(dImQ1n,2.));
 
 // Re[Q_{4n} Q_{2n}^* Q_{n}^* Q_{n}^*]
 Double_t reQ4nQ2nstarQ1nstarQ1nstar = (dReQ4n*dReQ2n+dImQ4n*dImQ2n)*(pow(dReQ1n,2)-pow(dImQ1n,2))
                                     + 2.*dReQ1n*dImQ1n*(dImQ4n*dReQ2n-dReQ4n*dImQ2n); 
 
 // Im[Q_{4n} Q_{2n}^* Q_{n}^* Q_{n}^*]
 //Double_t imQ4nQ2nstarQ1nstarQ1nstar; //calculate and implement this (deleteMe)
 
 // Re[Q_{2n} Q_{n} Q_{n}^* Q_{n}^* Q_{n}^*]
 Double_t reQ2nQ1nQ1nstarQ1nstarQ1nstar = (dReQ2n*dReQ1n-dImQ2n*dImQ1n)*(pow(dReQ1n,3)-3.*dReQ1n*pow(dImQ1n,2))
                                        + (dReQ2n*dImQ1n+dReQ1n*dImQ2n)*(3.*dImQ1n*pow(dReQ1n,2)-pow(dImQ1n,3));

 // Im[Q_{2n} Q_{n} Q_{n}^* Q_{n}^* Q_{n}^*] 
 //Double_t imQ2nQ1nQ1nstarQ1nstarQ1nstar; //calculate and implement this (deleteMe)
 
 // Re[Q_{2n} Q_{2n} Q_{2n}^* Q_{n}^* Q_{n}^*]
 Double_t reQ2nQ2nQ2nstarQ1nstarQ1nstar = (pow(dReQ2n,2.)+pow(dImQ2n,2.))
                                        * (dReQ2n*(pow(dReQ1n,2.)-pow(dImQ1n,2.)) + 2.*dImQ2n*dReQ1n*dImQ1n);

 // Im[Q_{2n} Q_{2n} Q_{2n}^* Q_{n}^* Q_{n}^*]
 //Double_t imQ2nQ2nQ2nstarQ1nstarQ1nstar = (pow(dReQ2n,2.)+pow(dImQ2n,2.))
 //                                       * (dImQ2n*(pow(dReQ1n,2.)-pow(dImQ1n,2.)) - 2.*dReQ2n*dReQ1n*dImQ1n);
 
 // Re[Q_{4n} Q_{n}^* Q_{n}^* Q_{n}^* Q_{n}^*]
 Double_t reQ4nQ1nstarQ1nstarQ1nstarQ1nstar = pow(dReQ1n,4.)*dReQ4n-6.*pow(dReQ1n,2.)*dReQ4n*pow(dImQ1n,2.)
                                            + pow(dImQ1n,4.)*dReQ4n+4.*pow(dReQ1n,3.)*dImQ1n*dImQ4n
                                            - 4.*pow(dImQ1n,3.)*dReQ1n*dImQ4n;
                                            
 // Im[Q_{4n} Q_{n}^* Q_{n}^* Q_{n}^* Q_{n}^*]
 //Double_t imQ4nQ1nstarQ1nstarQ1nstarQ1nstar = pow(dReQ1n,4.)*dImQ4n-6.*pow(dReQ1n,2.)*dImQ4n*pow(dImQ1n,2.)
 //                                           + pow(dImQ1n,4.)*dImQ4n+4.*pow(dImQ1n,3.)*dReQ1n*dReQ4n
 //                                           - 4.*pow(dReQ1n,3.)*dImQ1n*dReQ4n;
 
 // Re[Q_{3n} Q_{n} Q_{2n}^* Q_{n}^* Q_{n}^*]
 Double_t reQ3nQ1nQ2nstarQ1nstarQ1nstar = (pow(dReQ1n,2.)+pow(dImQ1n,2.))
                                        * (dReQ1n*dReQ2n*dReQ3n-dReQ3n*dImQ1n*dImQ2n+dReQ2n*dImQ1n*dImQ3n+dReQ1n*dImQ2n*dImQ3n);
 
 // Im[Q_{3n} Q_{n} Q_{2n}^* Q_{n}^* Q_{n}^*]
 //Double_t imQ3nQ1nQ2nstarQ1nstarQ1nstar = (pow(dReQ1n,2.)+pow(dImQ1n,2.))
 //                                       * (-dReQ2n*dReQ3n*dImQ1n-dReQ1n*dReQ3n*dImQ2n+dReQ1n*dReQ2n*dImQ3n-dImQ1n*dImQ2n*dImQ3n);
 
 
 // Re[Q_{2n} Q_{2n} Q_{n}^* Q_{n}^* Q_{n}^* Q_{n}^*]
 Double_t reQ2nQ2nQ1nstarQ1nstarQ1nstarQ1nstar = (pow(dReQ1n,2.)*dReQ2n-2.*dReQ1n*dReQ2n*dImQ1n-dReQ2n*pow(dImQ1n,2.)
                                               + dImQ2n*pow(dReQ1n,2.)+2.*dReQ1n*dImQ1n*dImQ2n-pow(dImQ1n,2.)*dImQ2n)
                                               * (pow(dReQ1n,2.)*dReQ2n+2.*dReQ1n*dReQ2n*dImQ1n-dReQ2n*pow(dImQ1n,2.)
                                               - dImQ2n*pow(dReQ1n,2.)+2.*dReQ1n*dImQ1n*dImQ2n+pow(dImQ1n,2.)*dImQ2n);
 
 // Im[Q_{2n} Q_{2n} Q_{n}^* Q_{n}^* Q_{n}^* Q_{n}^*]
 //Double_t imQ2nQ2nQ1nstarQ1nstarQ1nstarQ1nstar = 2.*(pow(dReQ1n,2.)*dReQ2n-dReQ2n*pow(dImQ1n,2.)
 //                                              + 2.*dReQ1n*dImQ1n*dImQ2n)*(pow(dReQ1n,2.)*dImQ2n
 //                                              - 2.*dReQ1n*dImQ1n*dReQ2n-pow(dImQ1n,2.)*dImQ2n);
 
 // Re[Q_{3n} Q_{n} Q_{n}^* Q_{n}^* Q_{n}^* Q_{n}^*]
 Double_t reQ3nQ1nQ1nstarQ1nstarQ1nstarQ1nstar = (pow(dReQ1n,2.)+pow(dImQ1n,2.))
                                               * (pow(dReQ1n,3.)*dReQ3n-3.*dReQ1n*dReQ3n*pow(dImQ1n,2.)
                                               + 3.*pow(dReQ1n,2.)*dImQ1n*dImQ3n-pow(dImQ1n,3.)*dImQ3n);
  
 // Im[Q_{3n} Q_{n} Q_{n}^* Q_{n}^* Q_{n}^* Q_{n}^*]                                                                                           
 //Double_t imQ3nQ1nQ1nstarQ1nstarQ1nstarQ1nstar = (pow(dReQ1n,2.)+pow(dImQ1n,2.))
 //                                              * (pow(dImQ1n,3.)*dReQ3n-3.*dImQ1n*dReQ3n*pow(dReQ1n,2.)
 //                                              - 3.*pow(dImQ1n,2.)*dReQ1n*dImQ3n+pow(dReQ1n,3.)*dImQ3n);
 
 // |Q_{2n}|^2 |Q_{n}|^4
 Double_t dQ2nQ1nQ1nQ2nstarQ1nstarQ1nstar = (pow(dReQ2n,2.)+pow(dImQ2n,2.))*pow((pow(dReQ1n,2.)+pow(dImQ1n,2.)),2.);
 
 // Re[Q_{2n} Q_{n} Q_{n} Q_{n}^* Q_{n}^* Q_{n}^* Q_{n}^*]
 Double_t reQ2nQ1nQ1nQ1nstarQ1nstarQ1nstarQ1nstar = pow((pow(dReQ1n,2.)+pow(dImQ1n,2.)),2.)
                                                  * (pow(dReQ1n,2.)*dReQ2n-dReQ2n*pow(dImQ1n,2.)
                                                  + 2.*dReQ1n*dImQ1n*dImQ2n);
                                                  
 // Im[Q_{2n} Q_{n} Q_{n} Q_{n}^* Q_{n}^* Q_{n}^* Q_{n}^*]                                                  
 //Double_t imQ2nQ1nQ1nQ1nstarQ1nstarQ1nstarQ1nstar = pow((pow(dReQ1n,2.)+pow(dImQ1n,2.)),2.)
 //                                                 * (pow(dReQ1n,2.)*dImQ2n-dImQ2n*pow(dImQ1n,2.)
 //                                                 - 2.*dReQ1n*dReQ2n*dImQ1n);
 
  
 
       
 //                                        **************************************
 //                                        **** multi-particle correlations: ****
 //                                        **************************************
 //
 // Remark 1: multi-particle correlations calculated with non-weighted Q-vectors are stored in 1D profile fQCorrelations.
 // Remark 2: binning of fQCorrelations is organized as follows:
 // --------------------------------------------------------------------------------------------------------------------
 //  1st bin: <2>_{1n|1n} = two1n1n = cos(n*(phi1-phi2))>
 //  2nd bin: <2>_{2n|2n} = two2n2n = cos(2n*(phi1-phi2))>
 //  3rd bin: <2>_{3n|3n} = two3n3n = cos(3n*(phi1-phi2))> 
 //  4th bin: <2>_{4n|4n} = two4n4n = cos(4n*(phi1-phi2))>
 //  5th bin:           ----  EMPTY ----
 //  6th bin: <3>_{2n|1n,1n} = three2n1n1n = <cos(n*(2.*phi1-phi2-phi3))>
 //  7th bin: <3>_{3n|2n,1n} = three3n2n1n = <cos(n*(3.*phi1-2.*phi2-phi3))>
 //  8th bin: <3>_{4n|2n,2n} = three4n2n2n = <cos(n*(4.*phi1-2.*phi2-2.*phi3))>
 //  9th bin: <3>_{4n|3n,1n} = three4n3n1n = <cos(n*(4.*phi1-3.*phi2-phi3))>
 // 10th bin:           ----  EMPTY ----
 // 11th bin: <4>_{1n,1n|1n,1n} = four1n1n1n1n = <cos(n*(phi1+phi2-phi3-phi4))>
 // 12th bin: <4>_{2n,1n|2n,1n} = four2n1n2n1n = <cos(2.*n*(phi1+phi2-phi3-phi4))>
 // 13th bin: <4>_{2n,2n|2n,2n} = four2n2n2n2n = <cos(n*(2.*phi1+phi2-2.*phi3-phi4))>
 // 14th bin: <4>_{3n|1n,1n,1n} = four3n1n1n1n = <cos(n*(3.*phi1-phi2-phi3-phi4))> 
 // 15th bin: <4>_{3n,1n|3n,1n} = four3n1n3n1n = <cos(n*(4.*phi1-2.*phi2-phi3-phi4))>
 // 16th bin: <4>_{3n,1n|2n,2n} = four3n1n2n2n = <cos(n*(3.*phi1+phi2-2.*phi3-2.*phi4))>
 // 17th bin: <4>_{4n|2n,1n,1n} = four4n2n1n1n = <cos(n*(3.*phi1+phi2-3.*phi3-phi4))> 
 // 18th bin:           ----  EMPTY ----
 // 19th bin: <5>_{2n|1n,1n,1n,1n} = five2n1n1n1n1n = <cos(n*(2.*phi1+phi2-phi3-phi4-phi5))>
 // 20th bin: <5>_{2n,2n|2n,1n,1n} = five2n2n2n1n1n = <cos(n*(2.*phi1+2.*phi2-2.*phi3-phi4-phi5))>
 // 21st bin: <5>_{3n,1n|2n,1n,1n} = five3n1n2n1n1n = <cos(n*(3.*phi1+phi2-2.*phi3-phi4-phi5))>
 // 22nd bin: <5>_{4n|1n,1n,1n,1n} = five4n1n1n1n1n = <cos(n*(4.*phi1-phi2-phi3-phi4-phi5))>
 // 23rd bin:           ----  EMPTY ----
 // 24th bin: <6>_{1n,1n,1n|1n,1n,1n} = six1n1n1n1n1n1n = <cos(n*(phi1+phi2+phi3-phi4-phi5-phi6))>
 // 25th bin: <6>_{2n,1n,1n|2n,1n,1n} = six2n1n1n2n1n1n = <cos(n*(2.*phi1+2.*phi2-phi3-phi4-phi5-phi6))>
 // 26th bin: <6>_{2n,2n|1n,1n,1n,1n} = six2n2n1n1n1n1n = <cos(n*(3.*phi1+phi2-phi3-phi4-phi5-phi6))>
 // 27th bin: <6>_{3n,1n|1n,1n,1n,1n} = six3n1n1n1n1n1n = <cos(n*(2.*phi1+phi2+phi3-2.*phi4-phi5-phi6))>
 // 28th bin:           ----  EMPTY ----
 // 29th bin: <7>_{2n,1n,1n|1n,1n,1n,1n} = seven2n1n1n1n1n1n1n =  <cos(n*(2.*phi1+phi2+phi3-phi4-phi5-phi6-phi7))>
 // 30th bin:           ----  EMPTY ----
 // 31st bin: <8>_{1n,1n,1n,1n|1n,1n,1n,1n} = eight1n1n1n1n1n1n1n1n = <cos(n*(phi1+phi2+phi3+phi4-phi5-phi6-phi7-phi8))>
 // --------------------------------------------------------------------------------------------------------------------
    
 // 2-particle:
 Double_t two1n1n = 0.; // <cos(n*(phi1-phi2))>
 Double_t two2n2n = 0.; // <cos(2n*(phi1-phi2))>
 Double_t two3n3n = 0.; // <cos(3n*(phi1-phi2))>
 Double_t two4n4n = 0.; // <cos(4n*(phi1-phi2))>
 
 if(dMult>1)
 {
  two1n1n = (pow(dReQ1n,2.)+pow(dImQ1n,2.)-dMult)/(dMult*(dMult-1.)); 
  two2n2n = (pow(dReQ2n,2.)+pow(dImQ2n,2.)-dMult)/(dMult*(dMult-1.)); 
  two3n3n = (pow(dReQ3n,2.)+pow(dImQ3n,2.)-dMult)/(dMult*(dMult-1.)); 
  two4n4n = (pow(dReQ4n,2.)+pow(dImQ4n,2.)-dMult)/(dMult*(dMult-1.)); 
    
  fQCorrelations->Fill(0.,two1n1n,dMult*(dMult-1.));  
  fQCorrelations->Fill(1.,two2n2n,dMult*(dMult-1.)); 
  fQCorrelations->Fill(2.,two3n3n,dMult*(dMult-1.)); 
  fQCorrelations->Fill(3.,two4n4n,dMult*(dMult-1.)); 
  
  // distribution of <cos(n*(phi1-phi2))>:
  f2pDistribution->Fill(two1n1n,dMult*(dMult-1.)); 
 } // end of if(dMult>1)
 
 // 3-particle:
 Double_t three2n1n1n = 0.; // <cos(n*(2.*phi1-phi2-phi3))>
 Double_t three3n2n1n = 0.; // <cos(n*(3.*phi1-2.*phi2-phi3))>
 Double_t three4n2n2n = 0.; // <cos(n*(4.*phi1-2.*phi2-2.*phi3))>
 Double_t three4n3n1n = 0.; // <cos(n*(4.*phi1-3.*phi2-phi3))>
 
 if(dMult>2)
 {
  three2n1n1n = (reQ2nQ1nstarQ1nstar-2.*(pow(dReQ1n,2.)+pow(dImQ1n,2.))
              - (pow(dReQ2n,2.)+pow(dImQ2n,2.))+2.*dMult)
              / (dMult*(dMult-1.)*(dMult-2.));              
  three3n2n1n = (reQ3nQ2nstarQ1nstar-(pow(dReQ3n,2.)+pow(dImQ3n,2.))
              - (pow(dReQ2n,2.)+pow(dImQ2n,2.))
              - (pow(dReQ1n,2.)+pow(dImQ1n,2.))+2.*dMult)
              / (dMult*(dMult-1.)*(dMult-2.));
  three4n2n2n = (reQ4nQ2nstarQ2nstar-2.*(pow(dReQ2n,2.)+pow(dImQ2n,2.))
              - (pow(dReQ4n,2.)+pow(dImQ4n,2.))+2.*dMult)
              / (dMult*(dMult-1.)*(dMult-2.)); 
  three4n3n1n = (reQ4nQ3nstarQ1nstar-(pow(dReQ4n,2.)+pow(dImQ4n,2.))
              - (pow(dReQ3n,2.)+pow(dImQ3n,2.))
              - (pow(dReQ1n,2.)+pow(dImQ1n,2.))+2.*dMult)
              / (dMult*(dMult-1.)*(dMult-2.)); 
              
  fQCorrelations->Fill(5.,three2n1n1n,dMult*(dMult-1.)*(dMult-2.)); 
  fQCorrelations->Fill(6.,three3n2n1n,dMult*(dMult-1.)*(dMult-2.));
  fQCorrelations->Fill(7.,three4n2n2n,dMult*(dMult-1.)*(dMult-2.)); 
  fQCorrelations->Fill(8.,three4n3n1n,dMult*(dMult-1.)*(dMult-2.));    
 } // end of if(dMult>2)
 
 // 4-particle:
 Double_t four1n1n1n1n = 0.; // <cos(n*(phi1+phi2-phi3-phi4))>
 Double_t four2n2n2n2n = 0.; // <cos(2.*n*(phi1+phi2-phi3-phi4))>
 Double_t four2n1n2n1n = 0.; // <cos(n*(2.*phi1+phi2-2.*phi3-phi4))> 
 Double_t four3n1n1n1n = 0.; // <cos(n*(3.*phi1-phi2-phi3-phi4))> 
 Double_t four4n2n1n1n = 0.; // <cos(n*(4.*phi1-2.*phi2-phi3-phi4))> 
 Double_t four3n1n2n2n = 0.; // <cos(n*(3.*phi1+phi2-2.*phi3-2.*phi4))> 
 Double_t four3n1n3n1n = 0.; // <cos(n*(3.*phi1+phi2-3.*phi3-phi4))>   
 
 if(dMult>3)
 {
  four1n1n1n1n = (2.*dMult*(dMult-3.)+pow((pow(dReQ1n,2.)+pow(dImQ1n,2.)),2.)-4.*(dMult-2.)*(pow(dReQ1n,2.)
               + pow(dImQ1n,2.))-2.*reQ2nQ1nstarQ1nstar+(pow(dReQ2n,2.)+pow(dImQ2n,2.)))
               / (dMult*(dMult-1)*(dMult-2.)*(dMult-3.));     
  four2n2n2n2n = (2.*dMult*(dMult-3.)+pow((pow(dReQ2n,2.)+pow(dImQ2n,2.)),2.)-4.*(dMult-2.)*(pow(dReQ2n,2.)
               + pow(dImQ2n,2.))-2.*reQ4nQ2nstarQ2nstar+(pow(dReQ4n,2.)+pow(dImQ4n,2.)))
               / (dMult*(dMult-1)*(dMult-2.)*(dMult-3.));
  four2n1n2n1n = (dQ2nQ1nQ2nstarQ1nstar-2.*reQ3nQ2nstarQ1nstar-2.*reQ2nQ1nstarQ1nstar)
               / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.))
               - ((dMult-5.)*(pow(dReQ1n,2.)+pow(dImQ1n,2.))
               + (dMult-4.)*(pow(dReQ2n,2.)+pow(dImQ2n,2.))-(pow(dReQ3n,2.)+pow(dImQ3n,2.)))
               / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.))
               + (dMult-6.)/((dMult-1.)*(dMult-2.)*(dMult-3.));
  four3n1n1n1n = (reQ3nQ1nstarQ1nstarQ1nstar-3.*reQ3nQ2nstarQ1nstar-3.*reQ2nQ1nstarQ1nstar)
               / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.))
               + (2.*(pow(dReQ3n,2.)+pow(dImQ3n,2.))+3.*(pow(dReQ2n,2.)+pow(dImQ2n,2.))
               + 6.*(pow(dReQ1n,2.)+pow(dImQ1n,2.))-6.*dMult)
               / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.));
  four4n2n1n1n = (reQ4nQ2nstarQ1nstarQ1nstar-2.*reQ4nQ3nstarQ1nstar-reQ4nQ2nstarQ2nstar-2.*reQ3nQ2nstarQ1nstar)
               / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.))
               - (reQ2nQ1nstarQ1nstar-2.*(pow(dReQ4n,2.)+pow(dImQ4n,2.))-2.*(pow(dReQ3n,2.)+pow(dImQ3n,2.))
               - 3.*(pow(dReQ2n,2.)+pow(dImQ2n,2.))-4.*(pow(dReQ1n,2.)+pow(dImQ1n,2.)))
               / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.))
               - 6./((dMult-1.)*(dMult-2.)*(dMult-3.));
  four3n1n2n2n = (reQ3nQ1nQ2nstarQ2nstar-reQ4nQ2nstarQ2nstar-reQ3nQ1nQ4nstar-2.*reQ3nQ2nstarQ1nstar)
               / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.))
               - (2.*reQ1nQ1nQ2nstar-(pow(dReQ4n,2.)+pow(dImQ4n,2.))-2.*(pow(dReQ3n,2.)+pow(dImQ3n,2.))
               - 4.*(pow(dReQ2n,2.)+pow(dImQ2n,2.))-4.*(pow(dReQ1n,2.)+pow(dImQ1n,2.)))
               / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.))
               - 6./((dMult-1.)*(dMult-2.)*(dMult-3.)); 
  four3n1n3n1n = ((pow(dReQ3n,2.)+pow(dImQ3n,2.))*(pow(dReQ1n,2.)+pow(dImQ1n,2.))
               - 2.*reQ4nQ3nstarQ1nstar-2.*reQ3nQ2nstarQ1nstar)
               / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.))
               + ((pow(dReQ4n,2.)+pow(dImQ4n,2.))-(dMult-4.)*(pow(dReQ3n,2.)+pow(dImQ3n,2.))
               + (pow(dReQ2n,2.)+pow(dImQ2n,2.))-(dMult-4.)*(pow(dReQ1n,2.)+pow(dImQ1n,2.)))
               / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.))
               + (dMult-6.)/((dMult-1.)*(dMult-2.)*(dMult-3.));
               
  fQCorrelations->Fill(10.,four1n1n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.));
  fQCorrelations->Fill(11.,four2n1n2n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.));
  fQCorrelations->Fill(12.,four2n2n2n2n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.));
  fQCorrelations->Fill(13.,four3n1n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.));
  fQCorrelations->Fill(14.,four3n1n3n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.));
  fQCorrelations->Fill(15.,four3n1n2n2n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.));  
  fQCorrelations->Fill(16.,four4n2n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)); 
  
  // distribution of <cos(n*(phi1+phi2-phi3-phi4))>
  f4pDistribution->Fill(four1n1n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.));
  
  // fQProduct->Fill(0.,two1n1n*four1n1n1n1n,dMult*(dMult-1.)*dMult*(dMult-1.)*(dMult-2.)*(dMult-3.));
 } // end of if(dMult>3)

 // 5-particle:
 Double_t five2n1n1n1n1n = 0.; // <cos(n*(2.*phi1+phi2-phi3-phi4-phi5))>
 Double_t five2n2n2n1n1n = 0.; // <cos(n*(2.*phi1+2.*phi2-2.*phi3-phi4-phi5))>
 Double_t five3n1n2n1n1n = 0.; // <cos(n*(3.*phi1+phi2-2.*phi3-phi4-phi5))>
 Double_t five4n1n1n1n1n = 0.; // <cos(n*(4.*phi1-phi2-phi3-phi4-phi5))>
 
 if(dMult>4)
 {
  five2n1n1n1n1n = (reQ2nQ1nQ1nstarQ1nstarQ1nstar-reQ3nQ1nstarQ1nstarQ1nstar+6.*reQ3nQ2nstarQ1nstar)
                 / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.))
                 - (reQ2nQ1nQ3nstar+3.*(dMult-6.)*reQ2nQ1nstarQ1nstar+3.*reQ1nQ1nQ2nstar)
                 / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.))
                 - (2.*(pow(dReQ3n,2.)+pow(dImQ3n,2.))
                 + 3.*(pow(dReQ1n,2.)+pow(dImQ1n,2.))*(pow(dReQ2n,2.)+pow(dImQ2n,2.))     
                 - 3.*(dMult-4.)*(pow(dReQ2n,2.)+pow(dImQ2n,2.)))
                 / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.))
                 - 3.*(pow((pow(dReQ1n,2.)+pow(dImQ1n,2.)),2.)
                 - 2.*(2*dMult-5.)*(pow(dReQ1n,2.)+pow(dImQ1n,2.))+2.*dMult*(dMult-4.))
                 / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.));
                 
  five2n2n2n1n1n = (reQ2nQ2nQ2nstarQ1nstarQ1nstar-reQ4nQ2nstarQ1nstarQ1nstar-2.*reQ2nQ2nQ3nstarQ1nstar)
                 / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.))
                 + 2.*(reQ4nQ2nstarQ2nstar+4.*reQ3nQ2nstarQ1nstar+reQ3nQ1nQ4nstar)
                 / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.))
                 + (reQ2nQ2nQ4nstar-2.*(dMult-5.)*reQ2nQ1nstarQ1nstar+2.*reQ1nQ1nQ2nstar)
                 / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.))
                 - (2.*(pow(dReQ4n,2.)+pow(dImQ4n,2.))+4.*(pow(dReQ3n,2.)+pow(dImQ3n,2.))
                 + 1.*pow((pow(dReQ2n,2.)+pow(dImQ2n,2.)),2.)
                 - 2.*(3.*dMult-10.)*(pow(dReQ2n,2.)+pow(dImQ2n,2.)))
                 / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.))
                 - (4.*(pow(dReQ1n,2.)+pow(dImQ1n,2.))*(pow(dReQ2n,2.)+pow(dImQ2n,2.))
                 - 4.*(dMult-5.)*(pow(dReQ1n,2.)+pow(dImQ1n,2.))+4.*dMult*(dMult-6.))
                 / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)); 

  five4n1n1n1n1n = (reQ4nQ1nstarQ1nstarQ1nstarQ1nstar-6.*reQ4nQ2nstarQ1nstarQ1nstar-4.*reQ3nQ1nstarQ1nstarQ1nstar)
                 / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.))
                 + (8.*reQ4nQ3nstarQ1nstar+3.*reQ4nQ2nstarQ2nstar+12.*reQ3nQ2nstarQ1nstar+12.*reQ2nQ1nstarQ1nstar)
                 / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.))
                 - (6.*(pow(dReQ4n,2.)+pow(dImQ4n,2.))+8.*(pow(dReQ3n,2.)+pow(dImQ3n,2.))
                 + 12.*(pow(dReQ2n,2.)+pow(dImQ2n,2.))+24.*(pow(dReQ1n,2.)+pow(dImQ1n,2.))-24.*dMult)
                 / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.));
  
  five3n1n2n1n1n = (reQ3nQ1nQ2nstarQ1nstarQ1nstar-reQ4nQ2nstarQ1nstarQ1nstar-reQ3nQ1nstarQ1nstarQ1nstar)
                 / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.))
                 - (reQ3nQ1nQ2nstarQ2nstar-3.*reQ4nQ3nstarQ1nstar-reQ4nQ2nstarQ2nstar)
                 / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.))
                 - ((2.*dMult-13.)*reQ3nQ2nstarQ1nstar-reQ3nQ1nQ4nstar-9.*reQ2nQ1nstarQ1nstar)
                 / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.))
                 - (2.*reQ1nQ1nQ2nstar+2.*(pow(dReQ4n,2.)+pow(dImQ4n,2.))
                 - 2.*(dMult-5.)*(pow(dReQ3n,2.)+pow(dImQ3n,2.))+2.*(pow(dReQ3n,2.)
                 + pow(dImQ3n,2.))*(pow(dReQ1n,2.)+pow(dImQ1n,2.)))
                 / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.))
                 + (2.*(dMult-6.)*(pow(dReQ2n,2.)+pow(dImQ2n,2.))
                 - 2.*(pow(dReQ2n,2.)+pow(dImQ2n,2.))*(pow(dReQ1n,2.)+pow(dImQ1n,2.))
                 - pow((pow(dReQ1n,2.)+pow(dImQ1n,2.)),2.)
                 + 2.*(3.*dMult-11.)*(pow(dReQ1n,2.)+pow(dImQ1n,2.)))
                 / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.))
                 - 4.*(dMult-6.)/((dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.));
                 
  fQCorrelations->Fill(18.,five2n1n1n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)); 
  fQCorrelations->Fill(19.,five2n2n2n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.));
  fQCorrelations->Fill(20.,five3n1n2n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.));
  fQCorrelations->Fill(21.,five4n1n1n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.));
 } // end of if(dMult>4)
    
 // 6-particle:
 Double_t six1n1n1n1n1n1n = 0.; // <cos(n*(phi1+phi2+phi3-phi4-phi5-phi6))>
 Double_t six2n2n1n1n1n1n = 0.; // <cos(n*(2.*phi1+2.*phi2-phi3-phi4-phi5-phi6))>
 Double_t six3n1n1n1n1n1n = 0.; // <cos(n*(3.*phi1+phi2-phi3-phi4-phi5-phi6))>
 Double_t six2n1n1n2n1n1n = 0.; // <cos(n*(2.*phi1+phi2+phi3-2.*phi4-phi5-phi6))>
 
 if(dMult>5)
 {
  six1n1n1n1n1n1n = (pow(pow(dReQ1n,2.)+pow(dImQ1n,2.),3.)+9.*dQ2nQ1nQ2nstarQ1nstar-6.*reQ2nQ1nQ1nstarQ1nstarQ1nstar)
                  / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)*(dMult-5.))
                  + 4.*(reQ3nQ1nstarQ1nstarQ1nstar-3.*reQ3nQ2nstarQ1nstar)
                  / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)*(dMult-5.))
                  + 2.*(9.*(dMult-4.)*reQ2nQ1nstarQ1nstar+2.*(pow(dReQ3n,2.)+pow(dImQ3n,2.)))
                  / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)*(dMult-5.))
                  - 9.*(pow((pow(dReQ1n,2.)+pow(dImQ1n,2.)),2.)+(pow(dReQ2n,2.)+pow(dImQ2n,2.)))
                  / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-5.))
                  + (18.*(pow(dReQ1n,2.)+pow(dImQ1n,2.)))
                  / (dMult*(dMult-1)*(dMult-3)*(dMult-4))
                  - 6./((dMult-1.)*(dMult-2.)*(dMult-3.));
                  
  six2n1n1n2n1n1n = (dQ2nQ1nQ1nQ2nstarQ1nstarQ1nstar-dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)
                  * (2.*five2n2n2n1n1n+4.*five2n1n1n1n1n+4.*five3n1n2n1n1n+4.*four2n1n2n1n+1.*four1n1n1n1n)
                  - dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(4.*four1n1n1n1n+4.*two1n1n
                  + 2.*three2n1n1n+2.*three2n1n1n+4.*four3n1n1n1n+8.*three2n1n1n+2.*four4n2n1n1n
                  + 4.*four2n1n2n1n+2.*two2n2n+8.*four2n1n2n1n+4.*four3n1n3n1n+8.*three3n2n1n
                  + 4.*four3n1n2n2n+4.*four1n1n1n1n+4.*four2n1n2n1n+1.*four2n2n2n2n)
                  - dMult*(dMult-1.)*(dMult-2.)*(2.*three2n1n1n+8.*two1n1n+4.*two1n1n+2.
                  + 4.*two1n1n+4.*three2n1n1n+2.*two2n2n+4.*three2n1n1n+8.*three3n2n1n
                  + 8.*two2n2n+4.*three4n3n1n+4.*two3n3n+4.*three3n2n1n+4.*two1n1n
                  + 8.*three2n1n1n+4.*two1n1n+4.*three3n2n1n+4.*three2n1n1n+2.*two2n2n
                  + 4.*three3n2n1n+2.*three4n2n2n)-dMult*(dMult-1.)
                  * (4.*two1n1n+4.+4.*two1n1n+2.*two2n2n+1.+4.*two1n1n+4.*two2n2n+4.*two3n3n
                  + 1.+2.*two2n2n+1.*two4n4n)-dMult)
                  / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)*(dMult-5.)); // to be improved (direct formula needed)
 
  six2n2n1n1n1n1n = (reQ2nQ2nQ1nstarQ1nstarQ1nstarQ1nstar-dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)
                  * (five4n1n1n1n1n+8.*five2n1n1n1n1n+6.*five2n2n2n1n1n)-dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)
                  * (4.*four3n1n1n1n+6.*four4n2n1n1n+12.*three2n1n1n+12.*four1n1n1n1n+24.*four2n1n2n1n
                  + 4.*four3n1n2n2n+3.*four2n2n2n2n)-dMult*(dMult-1.)*(dMult-2.)*(6.*three2n1n1n+12.*three3n2n1n
                  + 4.*three4n3n1n+3.*three4n2n2n+8.*three2n1n1n+24.*two1n1n+12.*two2n2n+12.*three2n1n1n+8.*three3n2n1n
                  + 1.*three4n2n2n)-dMult*(dMult-1.)*(4.*two1n1n+6.*two2n2n+4.*two3n3n+1.*two4n4n+2.*two2n2n+8.*two1n1n+6.)-dMult)
                  / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)*(dMult-5.)); // to be improved (direct formula needed)
   
  six3n1n1n1n1n1n = (reQ3nQ1nQ1nstarQ1nstarQ1nstarQ1nstar-dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)
                  * (five4n1n1n1n1n+4.*five2n1n1n1n1n+6.*five3n1n2n1n1n+4.*four3n1n1n1n)
                  - dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(4.*four3n1n1n1n+6.*four4n2n1n1n+6.*four1n1n1n1n
                  + 12.*three2n1n1n+12.*four2n1n2n1n+6.*four3n1n1n1n+12.*three3n2n1n+4.*four3n1n3n1n+3.*four3n1n2n2n)
                  - dMult*(dMult-1.)*(dMult-2.)*(6.*three2n1n1n+12.*three3n2n1n+4.*three4n3n1n+3.*three4n2n2n+4.*two1n1n
                  + 12.*two1n1n+6.*three2n1n1n+12.*three2n1n1n+4.*three3n2n1n+12.*two2n2n+4.*three3n2n1n+4.*two3n3n+1.*three4n3n1n
                  + 6.*three3n2n1n)-dMult*(dMult-1.)*(4.*two1n1n+6.*two2n2n+4.*two3n3n+1.*two4n4n+1.*two1n1n+4.+6.*two1n1n+4.*two2n2n
                  + 1.*two3n3n)-dMult)/(dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)*(dMult-5.)); // to be improved (direct formula needed)
   
  fQCorrelations->Fill(23.,six1n1n1n1n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)*(dMult-5.)); 
  fQCorrelations->Fill(24.,six2n1n1n2n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)*(dMult-5.)); 
  fQCorrelations->Fill(25.,six2n2n1n1n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)*(dMult-5.));
  fQCorrelations->Fill(26.,six3n1n1n1n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)*(dMult-5.)); 

  // distribution of <cos(n*(phi1+phi2+phi3-phi4-phi5-phi6))>
  f6pDistribution->Fill(six1n1n1n1n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)*(dMult-5.)); 
  
  //fQProduct->Fill(1.,two1n1n*six1n1n1n1n1n1n,dMult*(dMult-1.)*dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)*(dMult-5.));
  //fQProduct->Fill(3.,four1n1n1n1n*six1n1n1n1n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)*(dMult-5.));
 } // end of if(dMult>5)
 
 // 7-particle:
 Double_t seven2n1n1n1n1n1n1n = 0.; // <cos(n*(2.*phi1+phi2+phi3-phi4-phi5-phi6-phi7))>
 
 if(dMult>6)
 {
  seven2n1n1n1n1n1n1n = (reQ2nQ1nQ1nQ1nstarQ1nstarQ1nstarQ1nstar-dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)*(dMult-5.)
                      * (2.*six3n1n1n1n1n1n+4.*six1n1n1n1n1n1n+1.*six2n2n1n1n1n1n+6.*six2n1n1n2n1n1n+8.*five2n1n1n1n1n)
                      - dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)*(1.*five4n1n1n1n1n +8.*five2n1n1n1n1n+8.*four3n1n1n1n
                      + 12.*five3n1n2n1n1n+4.*five2n1n1n1n1n+3.*five2n2n2n1n1n+6.*five2n2n2n1n1n+6.*four1n1n1n1n+24.*four1n1n1n1n
                      + 12.*five2n1n1n1n1n+12.*five2n1n1n1n1n+12.*three2n1n1n+24.*four2n1n2n1n+4.*five3n1n2n1n1n+4.*five2n1n1n1n1n)
                      - dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(4.*four3n1n1n1n+6.*four4n2n1n1n+12.*four1n1n1n1n+24.*three2n1n1n
                      + 24.*four2n1n2n1n+12.*four3n1n1n1n+24.*three3n2n1n+8.*four3n1n3n1n+6.*four3n1n2n2n+6.*three2n1n1n+12.*four1n1n1n1n
                      + 12.*four2n1n2n1n+6.*three2n1n1n+12.*four2n1n2n1n+4.*four3n1n2n2n+3.*four2n2n2n2n+4.*four1n1n1n1n+6.*three2n1n1n
                      + 24.*two1n1n+24.*four1n1n1n1n+4.*four3n1n1n1n+24.*two1n1n+24.*three2n1n1n+12.*two2n2n+24.*three2n1n1n+12.*four2n1n2n1n
                      + 8.*three3n2n1n+8.*four2n1n2n1n+1.*four4n2n1n1n)-dMult*(dMult-1.)*(dMult-2.)*(6.*three2n1n1n+1.*three2n1n1n+8.*two1n1n
                      + 12.*three3n2n1n+24.*two1n1n+12.*three2n1n1n+4.*three2n1n1n+8.*two1n1n+4.*three4n3n1n+24.*three2n1n1n+8.*three3n2n1n
                      + 12.*two1n1n+12.*two1n1n+3.*three4n2n2n+24.*two2n2n+6.*two2n2n+12.+12.*three3n2n1n+8.*two3n3n+12.*three2n1n1n+24.*two1n1n
                      + 4.*three3n2n1n+8.*three3n2n1n+2.*three4n3n1n+12.*two1n1n+8.*three2n1n1n+4.*three2n1n1n+2.*three3n2n1n+6.*two2n2n+8.*two2n2n
                      + 1.*three4n2n2n+4.*three3n2n1n+6.*three2n1n1n)-dMult*(dMult-1.)*(4.*two1n1n+2.*two1n1n+6.*two2n2n+8.+1.*two2n2n+4.*two3n3n
                      + 12.*two1n1n+4.*two1n1n+1.*two4n4n+8.*two2n2n+6.+2.*two3n3n+4.*two1n1n+1.*two2n2n)-dMult)
                      / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)*(dMult-5.)*(dMult-6.)); // to be improved (direct formula needed)
        
  fQCorrelations->Fill(28.,seven2n1n1n1n1n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)*(dMult-5.)*(dMult-6.));
 } // end of if(dMult>6)
 
 // 8-particle:
 Double_t eight1n1n1n1n1n1n1n1n = 0.; // <cos(n*(phi1+phi2+phi3+phi4-phi5-phi6-phi7-phi8))>
 if(dMult>7)
 {
  eight1n1n1n1n1n1n1n1n = (pow(pow(dReQ1n,2.)+pow(dImQ1n,2.),4.)-dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)*(dMult-5.)*(dMult-6.)
                        * (12.*seven2n1n1n1n1n1n1n+16.*six1n1n1n1n1n1n)-dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)*(dMult-5.)
                        * (8.*six3n1n1n1n1n1n+48.*six1n1n1n1n1n1n+6.*six2n2n1n1n1n1n+96.*five2n1n1n1n1n+72.*four1n1n1n1n+36.*six2n1n1n2n1n1n)
                        - dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)*(2.*five4n1n1n1n1n+32.*five2n1n1n1n1n+36.*four1n1n1n1n
                        + 32.*four3n1n1n1n+48.*five2n1n1n1n1n+48.*five3n1n2n1n1n+144.*five2n1n1n1n1n+288.*four1n1n1n1n+36.*five2n2n2n1n1n
                        + 144.*three2n1n1n+96.*two1n1n+144.*four2n1n2n1n)-dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)
                        * (8.*four3n1n1n1n+48.*four1n1n1n1n+12.*four4n2n1n1n+96.*four2n1n2n1n+96.*three2n1n1n+72.*three2n1n1n+144.*two1n1n
                        + 16.*four3n1n3n1n+48.*four3n1n1n1n+144.*four1n1n1n1n+72.*four1n1n1n1n+96.*three3n2n1n+24.*four3n1n2n2n+144.*four2n1n2n1n
                        + 288.*two1n1n+288.*three2n1n1n+9.*four2n2n2n2n+72.*two2n2n+24.)-dMult*(dMult-1.)*(dMult-2.)*(12.*three2n1n1n+16.*two1n1n
                        + 24.*three3n2n1n+48.*three2n1n1n+96.*two1n1n+8.*three4n3n1n+32.*three3n2n1n+96.*three2n1n1n+144.*two1n1n+6.*three4n2n2n
                        + 96.*two2n2n+36.*two2n2n+72.+48.*three3n2n1n+16.*two3n3n+72.*three2n1n1n+144.*two1n1n)-dMult*(dMult-1.)*(8.*two1n1n
                        + 12.*two2n2n+16.+8.*two3n3n+48.*two1n1n+1.*two4n4n+16.*two2n2n+18.)-dMult)
                        / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)*(dMult-5.)*(dMult-6.)*(dMult-7.)); // to be improved (direct formula needed)
  
  fQCorrelations->Fill(30.,eight1n1n1n1n1n1n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)*(dMult-5.)*(dMult-6.)*(dMult-7.));
 
  // distribution of <cos(n*(phi1+phi2+phi3+phi4-phi5-phi6-phi7-phi8))>
  f8pDistribution->Fill(eight1n1n1n1n1n1n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)*(dMult-5.)*(dMult-6.)*(dMult-7.));
  
 } // end of if(dMult>7) 
 
} // end of AliFlowAnalysisWithQCumulants::CalculateCorrelationsForIntegratedFlow()


//================================================================================================================================


void AliFlowAnalysisWithQCumulants::CalculateWeightedCorrelationsForIntegratedFlow()
{
 // calculate all weighted correlations needed for 'no-name' integrated flow and store them in 1D profile fQCorrelationsW
 
 // Remark 1: binning of fQCorrelationsW is organized as follows:
 //..............................................................................................
 //       ---- bins 1-20: 2-particle correlations ----
 // 1st bin: two1n1nW1W1 = <w1 w2 cos(n*(phi1-phi2))>
 // 2nd bin: two2n2nW2W2 = <w1^2 w2^2 cos(2n*(phi1-phi2))>
 // 3rd bin: two3n3nW3W3 = <w1^3 w2^3 cos(3n*(phi1-phi2))>
 // 4th bin: two4n4nW4W4 = <w1^4 w2^4 cos(4n*(phi1-phi2))>
 // 5th bin: two1n1nW3W1 = <w1^3 w2 cos(n*(phi1-phi2))>
 // 6th bin: two1n1nW1W1W2 = <w1 w2 w3^2 cos(n*(phi1-phi2))>  
 //       ---- bins 21-40: 3-particle correlations ----
 // 21st bin: three2n1n1nW2W1W1 = <w1^2 w2 w3 cos(n*(2phi1-phi2-phi3))> 
 //       ---- bins 41-60: 4-particle correlations ----
 // 41st bin: four1n1n1n1nW1W1W1W1 = <w1 w2 w3 w4 cos(n*(phi1+phi2-phi3-phi4))>
 //       ---- bins 61-80: 5-particle correlations ---- 
 //       ---- bins 81-100: 6-particle correlations ----
 //       ---- bins 101-120: 7-particle correlations ----
 //       ---- bins 121-140: 8-particle correlations ----
 //..............................................................................................
 
 // multiplicity (number of particles used to determine the reaction plane)
 Double_t dMult = (*fSMpk)(0,0);
 
 // real and imaginary parts of weighted Q-vectors evaluated in harmonics n, 2n, 3n and 4n: 
 Double_t dReQ1n1k = (*fReQ)(0,1);
 Double_t dReQ2n2k = (*fReQ)(1,2);
 Double_t dReQ3n3k = (*fReQ)(2,3);
 Double_t dReQ4n4k = (*fReQ)(3,4);
 Double_t dReQ1n3k = (*fReQ)(0,3);
 Double_t dImQ1n1k = (*fImQ)(0,1);
 Double_t dImQ2n2k = (*fImQ)(1,2);
 Double_t dImQ3n3k = (*fImQ)(2,3);
 Double_t dImQ4n4k = (*fImQ)(3,4);
 Double_t dImQ1n3k = (*fImQ)(0,3);

 // dMs are variables introduced in order to simplify some Eqs. bellow:
 //..............................................................................................
 Double_t dM11 = (*fSMpk)(1,1)-(*fSMpk)(0,2); // dM11 = sum_{i,j=1,i!=j}^M w_i w_j
 Double_t dM22 = (*fSMpk)(1,2)-(*fSMpk)(0,4); // dM22 = sum_{i,j=1,i!=j}^M w_i^2 w_j^2
 Double_t dM33 = (*fSMpk)(1,3)-(*fSMpk)(0,6); // dM33 = sum_{i,j=1,i!=j}^M w_i^3 w_j^3
 Double_t dM44 = (*fSMpk)(1,4)-(*fSMpk)(0,8); // dM44 = sum_{i,j=1,i!=j}^M w_i^4 w_j^4
 Double_t dM31 = (*fSMpk)(0,3)*(*fSMpk)(0,1)-(*fSMpk)(0,4); // dM31 = sum_{i,j=1,i!=j}^M w_i^3 w_j
 Double_t dM211 = (*fSMpk)(0,2)*(*fSMpk)(1,1)-2.*(*fSMpk)(0,3)*(*fSMpk)(0,1)
                - (*fSMpk)(1,2)+2.*(*fSMpk)(0,4); // dM211 = sum_{i,j,k=1,i!=j!=k}^M w_i^2 w_j w_k
 Double_t dM1111 = (*fSMpk)(3,1)-6.*(*fSMpk)(0,2)*(*fSMpk)(1,1)  
                  + 8.*(*fSMpk)(0,3)*(*fSMpk)(0,1)
                  + 3.*(*fSMpk)(1,2)-6.*(*fSMpk)(0,4); // dM1111 = sum_{i,j,k,l=1,i!=j!=k!=l}^M w_i w_j w_k w_l
 //..............................................................................................


 

 //                                        ***********************************************
 //                                        **** weighted multi-particle correlations: ****
 //                                        ***********************************************
 //.............................................................................................. 
 // weighted 2-particle correlations:
 Double_t two1n1nW1W1 = 0.; // <w1 w2 cos(n*(phi1-phi2))>
 Double_t two2n2nW2W2 = 0.; // <w1^2 w2^2 cos(2n*(phi1-phi2))>
 Double_t two3n3nW3W3 = 0.; // <w1^3 w2^3 cos(3n*(phi1-phi2))>
 Double_t two4n4nW4W4 = 0.; // <w1^4 w2^4 cos(4n*(phi1-phi2))>
 Double_t two1n1nW3W1 = 0.; // <w1^3 w2 cos(n*(phi1-phi2))>
 Double_t two1n1nW1W1W2 = 0.; // <w1 w2 w3^2 cos(n*(phi1-phi2))> 
 
 if(dMult>1) 
 { 
  if(dM11)
  {
   two1n1nW1W1 = (pow(dReQ1n1k,2)+pow(dImQ1n1k,2)-(*fSMpk)(0,2))/dM11; 
   fQCorrelationsW->Fill(0.,two1n1nW1W1,dM11);
  }
  if(dM22)
  {
   two2n2nW2W2 = (pow(dReQ2n2k,2)+pow(dImQ2n2k,2)-(*fSMpk)(0,4))/dM22; 
   fQCorrelationsW->Fill(1.,two2n2nW2W2,dM22); 
  }
  if(dM33)
  {
   two3n3nW3W3 = (pow(dReQ3n3k,2)+pow(dImQ3n3k,2)-(*fSMpk)(0,6))/dM33;
   fQCorrelationsW->Fill(2.,two3n3nW3W3,dM33);
  }
  if(dM44)
  {
   two4n4nW4W4 = (pow(dReQ4n4k,2)+pow(dImQ4n4k,2)-(*fSMpk)(0,8))/dM44; 
   fQCorrelationsW->Fill(3.,two4n4nW4W4,dM44);  
  } 
  if(dM31)
  {
   two1n1nW3W1 = (dReQ1n3k*dReQ1n1k+dImQ1n3k*dImQ1n1k-(*fSMpk)(0,4))/dM31; 
   fQCorrelationsW->Fill(4.,two1n1nW3W1,dM31);  
  } 
  if(dM211)
  {
   two1n1nW1W1W2 = ((*fSMpk)(0,2)*(pow(dReQ1n1k,2)+pow(dImQ1n1k,2)-(*fSMpk)(0,2))
                 - 2.*(dReQ1n3k*dReQ1n1k+dImQ1n3k*dImQ1n1k
                 - (*fSMpk)(0,4)))/dM211;
   fQCorrelationsW->Fill(5.,two1n1nW1W1W2,dM211);  
  }  
 } // end of if(dMult>1)
 //..............................................................................................
 
 //..............................................................................................
 // weighted 3-particle correlations:
 Double_t three2n1n1nW2W1W1 = 0.; // <w1^2 w2 w3 cos(n*(2phi1-phi2-phi3))>
 
 if(dMult>2) 
 { 
  if(dM211)
  {                                                       
   three2n1n1nW2W1W1 = (pow(dReQ1n1k,2.)*dReQ2n2k+2.*dReQ1n1k*dImQ1n1k*dImQ2n2k-pow(dImQ1n1k,2.)*dReQ2n2k
                     - 2.*(dReQ1n3k*dReQ1n1k+dImQ1n3k*dImQ1n1k)
                     - pow(dReQ2n2k,2)-pow(dImQ2n2k,2)
                     + 2.*(*fSMpk)(0,4))/dM211;                                                                               
   fQCorrelationsW->Fill(20.,three2n1n1nW2W1W1,dM211);
  } 
 } // end of if(dMult>2) 
 //..............................................................................................
 
 //..............................................................................................
 // weighted 4-particle correlations:
 Double_t four1n1n1n1nW1W1W1W1 = 0.; // <w1 w2 w3 w4 cos(n*(phi1+phi2-phi3-phi4))>
 if(dMult>3) 
 { 
  if(dM1111)
  {      
   four1n1n1n1nW1W1W1W1 = (pow(pow(dReQ1n1k,2.)+pow(dImQ1n1k,2.),2)
                        - 2.*(pow(dReQ1n1k,2.)*dReQ2n2k+2.*dReQ1n1k*dImQ1n1k*dImQ2n2k-pow(dImQ1n1k,2.)*dReQ2n2k)
                        + 8.*(dReQ1n3k*dReQ1n1k+dImQ1n3k*dImQ1n1k)
                        + (pow(dReQ2n2k,2)+pow(dImQ2n2k,2))
                        - 4.*(*fSMpk)(0,2)*(pow(dReQ1n1k,2)+pow(dImQ1n1k,2))
                        - 6.*(*fSMpk)(0,4)+2.*(*fSMpk)(1,2))/dM1111;                                       
   fQCorrelationsW->Fill(40.,four1n1n1n1nW1W1W1W1,dM1111);
  } 
 } // end of if(dMult>3) 
 //..............................................................................................
 
} // end of AliFlowAnalysisWithQCumulants::CalculateWeightedCorrelationsForIntegratedFlow()


//================================================================================================================================


void AliFlowAnalysisWithQCumulants::CalculateCorrelationsForDifferentialFlow(TString type)
{
 // calculate all correlations needed for differential flow for each (pt,eta) bin: 
 
 // pt and eta bin width:
 Double_t dBinWidthPt = 0.; // to be improved (should I promote this variable to data members?)
 Double_t dBinWidthEta = 0.; // to be improved (should I promote this variable to data members?)
 
 if(fnBinsPt) dBinWidthPt=(fPtMax-fPtMin)/fnBinsPt;  
 if(fnBinsEta) dBinWidthEta=(fEtaMax-fEtaMin)/fnBinsEta;  
 
 // multiplicity:
 Double_t dMult = (*fSMpk)(0,0);
 
 // real and imaginary parts of non-weighted Q-vectors evaluated in harmonics n, 2n, 3n and 4n: 
 Double_t dReQ1n = (*fReQ)(0,0);
 Double_t dReQ2n = (*fReQ)(1,0);
 //Double_t dReQ3n = (*fReQ)(2,0);
 //Double_t dReQ4n = (*fReQ)(3,0);
 Double_t dImQ1n = (*fImQ)(0,0);
 Double_t dImQ2n = (*fImQ)(1,0);
 //Double_t dImQ3n = (*fImQ)(2,0);
 //Double_t dImQ4n = (*fImQ)(3,0);

 // looping over all (pt,eta) bins and calculating correlations needed for differential flow: 
 for(Int_t p=1;p<=fnBinsPt;p++)
 {
  for(Int_t e=1;e<=fnBinsEta;e++)
  {
   // real and imaginary parts of q_n (non-weighted Q-vector evaluated only for POIs in harmonic n for each (pt,eta) bin): 
   Double_t dReqnPtEta = 0.;
   Double_t dImqnPtEta = 0.;

   // number of POIs in each (pt,eta) bin:
   Double_t dmPtEta = 0.;

   // real and imaginary parts of q''_{n}, q''_{2n}, ... 
   // (non-weighted Q-vectors evaluated only for particles which are both RPs and POIs in harmonic n, 2n, ... for each (pt,eta) bin): 
   Double_t dReqPrimePrime1nPtEta = 0.;
   Double_t dImqPrimePrime1nPtEta = 0.;
   Double_t dReqPrimePrime2nPtEta = 0.;
   Double_t dImqPrimePrime2nPtEta = 0.;

   // number of particles which are both RPs and POIs in each (pt,eta) bin:
   Double_t dmPrimePrimePtEta = 0.;
   
   if(type == "POI")
   {
    // q''_{n}, q''_{2n}:
    //...............................................................................................
    dReqPrimePrime1nPtEta = fReqPrimePrime1nPtEta->GetBinContent(fReqPrimePrime1nPtEta->GetBin(p,e));
    dImqPrimePrime1nPtEta = fImqPrimePrime1nPtEta->GetBinContent(fImqPrimePrime1nPtEta->GetBin(p,e));
    dReqPrimePrime2nPtEta = fReqPrimePrime2nPtEta->GetBinContent(fReqPrimePrime2nPtEta->GetBin(p,e));
    dImqPrimePrime2nPtEta = fImqPrimePrime2nPtEta->GetBinContent(fImqPrimePrime2nPtEta->GetBin(p,e));
    //...............................................................................................
   
    // m'':
    dmPrimePrimePtEta = fmPrimePrimePtEta->GetBinContent(fmPrimePrimePtEta->GetBin(p,e));
   
    // q'_{n}: 
    dReqnPtEta = fReqnPtEta->GetBinContent(fReqnPtEta->GetBin(p,e));
    dImqnPtEta = fImqnPtEta->GetBinContent(fImqnPtEta->GetBin(p,e));
    dmPtEta    = fmPtEta->GetBinContent(fmPtEta->GetBin(p,e));
   }
   else if(type == "RP")
   {
    // q_RP{n}, q_RP{2n}:
    //...............................................................................................
    dReqPrimePrime1nPtEta = fReqRP1nPtEta->GetBinContent(fReqRP1nPtEta->GetBin(p,e));
    dImqPrimePrime1nPtEta = fImqRP1nPtEta->GetBinContent(fImqRP1nPtEta->GetBin(p,e));
    dReqPrimePrime2nPtEta = fReqRP2nPtEta->GetBinContent(fReqRP2nPtEta->GetBin(p,e));
    dImqPrimePrime2nPtEta = fImqRP2nPtEta->GetBinContent(fImqRP2nPtEta->GetBin(p,e));
    //...............................................................................................
   
    // m'':
    dmPrimePrimePtEta = fmRPPtEta->GetBinContent(fmRPPtEta->GetBin(p,e));
   
    dReqnPtEta = fReqRP1nPtEta->GetBinContent(fReqRP1nPtEta->GetBin(p,e)); // not a bug ;-)
    dImqnPtEta = fImqRP1nPtEta->GetBinContent(fImqRP1nPtEta->GetBin(p,e)); // not a bug ;-)
    dmPtEta    = fmRPPtEta->GetBinContent(fmRPPtEta->GetBin(p,e));         // not a bug ;-) 
   }
   
   // 2'-particle correlation:
   Double_t two1n1nPtEta = 0.;
   if(dmPtEta*dMult-dmPrimePrimePtEta)
   {
    two1n1nPtEta = (dReqnPtEta*dReQ1n+dImqnPtEta*dImQ1n-dmPrimePrimePtEta)
                 / (dmPtEta*dMult-dmPrimePrimePtEta);
   
    // fill the 2D profile to get the average correlation for each (pt, eta) bin:
    if(type == "POI")
    {
     f2pPtEtaPOI->Fill(fPtMin+(p-1)*dBinWidthPt,fEtaMin+(e-1)*dBinWidthEta,two1n1nPtEta,dmPtEta*dMult-dmPrimePrimePtEta);
    }
    else if(type == "RP")
    {
     f2pPtEtaRP->Fill(fPtMin+(p-1)*dBinWidthPt,fEtaMin+(e-1)*dBinWidthEta,two1n1nPtEta,dmPtEta*dMult-dmPrimePrimePtEta);   
    }
   } // end of if(dmPtEta*dMult-dmPrimePrimePtEta)
   
   // 4'-particle correlation:
   Double_t four1n1n1n1nPtEta = 0.;
   if((dmPtEta-dmPrimePrimePtEta)*dMult*(dMult-1.)*(dMult-2.)
       + dmPrimePrimePtEta*(dMult-1.)*(dMult-2.)*(dMult-3.)) // to be improved (introduce a new variable for this expression)
   {
    four1n1n1n1nPtEta = ((pow(dReQ1n,2.)+pow(dImQ1n,2.))*(dReqnPtEta*dReQ1n+dImqnPtEta*dImQ1n)
                      - dReqPrimePrime2nPtEta*(pow(dReQ1n,2.)-pow(dImQ1n,2.))
                      - 2.*dImqPrimePrime2nPtEta*dReQ1n*dImQ1n
                      - dReqnPtEta*(dReQ1n*dReQ2n+dImQ1n*dImQ2n)
                      + dImqnPtEta*(dImQ1n*dReQ2n-dReQ1n*dImQ2n)
                      - 2.*dMult*(dReqnPtEta*dReQ1n+dImqnPtEta*dImQ1n)
                      - 2.*(pow(dReQ1n,2.)+pow(dImQ1n,2.))*dmPrimePrimePtEta                      
                      + 6.*(dReqPrimePrime1nPtEta*dReQ1n+dImqPrimePrime1nPtEta*dImQ1n)                                            
                      + 1.*(dReqPrimePrime2nPtEta*dReQ2n+dImqPrimePrime2nPtEta*dImQ2n)                      
                      + 2.*(dReqnPtEta*dReQ1n+dImqnPtEta*dImQ1n)                       
                      + 2.*dmPrimePrimePtEta*dMult                      
                      - 6.*dmPrimePrimePtEta)        
                      / ((dmPtEta-dmPrimePrimePtEta)*dMult*(dMult-1.)*(dMult-2.)
                          + dmPrimePrimePtEta*(dMult-1.)*(dMult-2.)*(dMult-3.)); 
    
    // fill the 2D profile to get the average correlation for each (pt, eta) bin:
    if(type == "POI")
    {
     f4pPtEtaPOI->Fill(fPtMin+(p-1)*dBinWidthPt,fEtaMin+(e-1)*dBinWidthEta,four1n1n1n1nPtEta,
                       (dmPtEta-dmPrimePrimePtEta)*dMult*(dMult-1.)*(dMult-2.)
                        + dmPrimePrimePtEta*(dMult-1.)*(dMult-2.)*(dMult-3.));
    }
    else if(type == "RP")
    {
     f4pPtEtaRP->Fill(fPtMin+(p-1)*dBinWidthPt,fEtaMin+(e-1)*dBinWidthEta,four1n1n1n1nPtEta,
                      (dmPtEta-dmPrimePrimePtEta)*dMult*(dMult-1.)*(dMult-2.)
                       + dmPrimePrimePtEta*(dMult-1.)*(dMult-2.)*(dMult-3.));   
    }
   } // end of if((dmPtEta-dmPrimePrimePtEta)*dMult*(dMult-1.)*(dMult-2.)
     //            +dmPrimePrimePtEta*(dMult-1.)*(dMult-2.)*(dMult-3.))
   
  } // end of for(Int_t e=1;e<=fnBinsEta;e++)
 } // end of for(Int_t p=1;p<=fnBinsPt;p++)
 
} // end of AliFlowAnalysisWithQCumulants::CalculateCorrelationsForDifferentialFlow()


//================================================================================================================================


void AliFlowAnalysisWithQCumulants::CalculateWeightedCorrelationsForDifferentialFlow(TString type)
{
 // calculate all weighted correlations needed for differential flow 
 
 // pt and eta bin width:
 Double_t dBinWidthPt = 0.; // to be improved (should I promote this variable to data members?)
 Double_t dBinWidthEta = 0.; // to be improved (should I promote this variable to data members?)
 
 if(fnBinsPt) dBinWidthPt=(fPtMax-fPtMin)/fnBinsPt;  
 if(fnBinsEta) dBinWidthEta=(fEtaMax-fEtaMin)/fnBinsEta; 

 // real and imaginary parts of weighted Q-vectors evaluated in harmonics n, 2n, 3n and 4n: 
 Double_t dReQ1n1k = (*fReQ)(0,1);
 Double_t dReQ2n2k = (*fReQ)(1,2);
 Double_t dReQ1n3k = (*fReQ)(0,3);
 //Double_t dReQ4n4k = (*fReQ)(3,4);
 Double_t dImQ1n1k = (*fImQ)(0,1);
 Double_t dImQ2n2k = (*fImQ)(1,2);
 Double_t dImQ1n3k = (*fImQ)(0,3);
 //Double_t dImQ4n4k = (*fImQ)(3,4);
 
 // S^M_{p,k} (see .h file for the definition of fSMpk):
 Double_t dSM1p1k = (*fSMpk)(0,1);
 Double_t dSM1p2k = (*fSMpk)(0,2);
 Double_t dSM1p3k = (*fSMpk)(0,3);
 Double_t dSM2p1k = (*fSMpk)(1,1);
 Double_t dSM3p1k = (*fSMpk)(2,1);
 
 // looping over all (pt,eta) bins and calculating weighted correlations needed for differential flow: 
 for(Int_t p=1;p<=fnBinsPt;p++)
 {
  for(Int_t e=1;e<=fnBinsEta;e++)
  {
   // real and imaginary parts of q_n (non-weighted Q-vector evaluated only for POIs in harmonic n for each (pt,eta) bin): 
   Double_t dReqnPtEta = 0.;
   Double_t dImqnPtEta = 0.;

   // number of POIs in each (pt,eta) bin:
   Double_t dmPtEta = 0.;

   // real and imaginary parts of q''_{n,2k}, q''_{2n,1k}, ... 
   // (weighted Q-vectors evaluated only for particles which are both RPs and POIs in harmonic n, 2n, ... for each (pt,eta) bin): 
   Double_t dReqPrimePrime1n2kPtEta = 0.;
   Double_t dImqPrimePrime1n2kPtEta = 0.;
   Double_t dReqPrimePrime2n1kPtEta = 0.;
   Double_t dImqPrimePrime2n1kPtEta = 0.;

   // S^{m''}_{1,1}, S^{m''}_{1,2}, S^{m''}_{1,3}... (see .h file for the definition): 
   Double_t dSmPrimePrime1p1kPtEta = 0.; 
   Double_t dSmPrimePrime1p2kPtEta = 0.; 
   Double_t dSmPrimePrime1p3kPtEta = 0.; 
   
   // M0111 from Eq. (118) in QC2c (to be improved (notation))
   Double_t dM0111 = 0.;
 
   // qPOI_{n}: // to be improved (notation)
   if(type == "POI")
   {
    dReqnPtEta = fReqnPtEta->GetBinContent(fReqnPtEta->GetBin(p,e));
    dImqnPtEta = fImqnPtEta->GetBinContent(fImqnPtEta->GetBin(p,e));
    dmPtEta    = fmPtEta->GetBinContent(fmPtEta->GetBin(p,e));
    
    //...............................................................................................
    // q''_{n,2k}, q''_{2n,1k}:
    dReqPrimePrime1n2kPtEta = fReqPrimePrime1n2kPtEta->GetBinContent(fReqPrimePrime1n2kPtEta->GetBin(p,e));
    dImqPrimePrime1n2kPtEta = fImqPrimePrime1n2kPtEta->GetBinContent(fImqPrimePrime1n2kPtEta->GetBin(p,e));
    dReqPrimePrime2n1kPtEta = fReqPrimePrime2n1kPtEta->GetBinContent(fReqPrimePrime2n1kPtEta->GetBin(p,e));
    dImqPrimePrime2n1kPtEta = fImqPrimePrime2n1kPtEta->GetBinContent(fImqPrimePrime2n1kPtEta->GetBin(p,e));
   
    // S^{m''}_{1,1}, S^{m''}_{1,2}, S^{m''}_{1,3}...: 
    dSmPrimePrime1p1kPtEta = fSmPrimePrime1p1kPtEta->GetBinContent(fSmPrimePrime1p1kPtEta->GetBin(p,e)); 
    dSmPrimePrime1p2kPtEta = fSmPrimePrime1p2kPtEta->GetBinContent(fSmPrimePrime1p2kPtEta->GetBin(p,e)); 
    dSmPrimePrime1p3kPtEta = fSmPrimePrime1p3kPtEta->GetBinContent(fSmPrimePrime1p3kPtEta->GetBin(p,e));
   
    // M0111 from Eq. (118) in QC2c (to be improved (notation)):
    dM0111 = dmPtEta*(dSM3p1k-3.*dSM1p1k*dSM1p2k+2.*dSM1p3k)
           - 3.*(dSmPrimePrime1p1kPtEta*(dSM2p1k-dSM1p2k)
           + 2.*(dSmPrimePrime1p3kPtEta-dSmPrimePrime1p2kPtEta*dSM1p1k));
    //...............................................................................................   
   }
   else if(type == "RP")
   {
    dReqnPtEta = fReqRP1nPtEta->GetBinContent(fReqRP1nPtEta->GetBin(p,e)); // not a bug ;-)
    dImqnPtEta = fImqRP1nPtEta->GetBinContent(fImqRP1nPtEta->GetBin(p,e)); // not a bug ;-)
    dmPtEta    = fmRPPtEta->GetBinContent(fmRPPtEta->GetBin(p,e));         // not a bug ;-) 
    
    //...............................................................................................
    // q''_{n,2k}, q''_{2n,1k}: (to be improved (notation)):
    dReqPrimePrime1n2kPtEta = fReqRP1n2kPtEta->GetBinContent(fReqRP1n2kPtEta->GetBin(p,e));
    dImqPrimePrime1n2kPtEta = fImqRP1n2kPtEta->GetBinContent(fImqRP1n2kPtEta->GetBin(p,e));
    dReqPrimePrime2n1kPtEta = fReqRP2n1kPtEta->GetBinContent(fReqRP2n1kPtEta->GetBin(p,e));
    dImqPrimePrime2n1kPtEta = fImqRP2n1kPtEta->GetBinContent(fImqRP2n1kPtEta->GetBin(p,e));
   
    // S^{m''}_{1,1}, S^{m''}_{1,2}, S^{m''}_{1,3}...:  (to be improved (notation)):
    dSmPrimePrime1p1kPtEta = fSmRP1p1kPtEta->GetBinContent(fSmRP1p1kPtEta->GetBin(p,e)); 
    dSmPrimePrime1p2kPtEta = fSmRP1p2kPtEta->GetBinContent(fSmRP1p2kPtEta->GetBin(p,e)); 
    dSmPrimePrime1p3kPtEta = fSmRP1p3kPtEta->GetBinContent(fSmRP1p3kPtEta->GetBin(p,e));
   
    // M0111 from Eq. (118) in QC2c (to be improved (notation)):
    dM0111 = dmPtEta*(dSM3p1k-3.*dSM1p1k*dSM1p2k+2.*dSM1p3k)
           - 3.*(dSmPrimePrime1p1kPtEta*(dSM2p1k-dSM1p2k)
           + 2.*(dSmPrimePrime1p3kPtEta-dSmPrimePrime1p2kPtEta*dSM1p1k));
    //...............................................................................................   
   }
   
   // 2'-particle correlation:
   Double_t two1n1nW0W1PtEta = 0.;
   if(dmPtEta*dSM1p1k-dSmPrimePrime1p1kPtEta)
   {
    two1n1nW0W1PtEta = (dReqnPtEta*dReQ1n1k+dImqnPtEta*dImQ1n1k-dSmPrimePrime1p1kPtEta)
                 / (dmPtEta*dSM1p1k-dSmPrimePrime1p1kPtEta);
   
    // fill the 2D profile to get the average correlation for each (pt, eta) bin:
    if(type == "POI")
    {
     f2pPtEtaPOIW->Fill(fPtMin+(p-1)*dBinWidthPt,fEtaMin+(e-1)*dBinWidthEta,two1n1nW0W1PtEta,
                        dmPtEta*dSM1p1k-dSmPrimePrime1p1kPtEta);
    }
    else if(type == "RP")
    {
     f2pPtEtaRPW->Fill(fPtMin+(p-1)*dBinWidthPt,fEtaMin+(e-1)*dBinWidthEta,two1n1nW0W1PtEta,
                       dmPtEta*dSM1p1k-dSmPrimePrime1p1kPtEta);   
    }
   } // end of if(dmPtEta*dMult-dmPrimePrimePtEta)
   
   // 4'-particle correlation:
   Double_t four1n1n1n1nW0W1W1W1PtEta = 0.;
   if(dM0111)
   {
    four1n1n1n1nW0W1W1W1PtEta = ((pow(dReQ1n1k,2.)+pow(dImQ1n1k,2.))*(dReqnPtEta*dReQ1n1k+dImqnPtEta*dImQ1n1k)
                      - dReqPrimePrime2n1kPtEta*(pow(dReQ1n1k,2.)-pow(dImQ1n1k,2.))
                      - 2.*dImqPrimePrime2n1kPtEta*dReQ1n1k*dImQ1n1k
                      - dReqnPtEta*(dReQ1n1k*dReQ2n2k+dImQ1n1k*dImQ2n2k)
                      + dImqnPtEta*(dImQ1n1k*dReQ2n2k-dReQ1n1k*dImQ2n2k)
                      - 2.*dSM1p2k*(dReqnPtEta*dReQ1n1k+dImqnPtEta*dImQ1n1k)
                      - 2.*(pow(dReQ1n1k,2.)+pow(dImQ1n1k,2.))*dSmPrimePrime1p1kPtEta                                            
                      + 6.*(dReqPrimePrime1n2kPtEta*dReQ1n1k+dImqPrimePrime1n2kPtEta*dImQ1n1k)                                           
                      + 1.*(dReqPrimePrime2n1kPtEta*dReQ2n2k+dImqPrimePrime2n1kPtEta*dImQ2n2k)                         
                      + 2.*(dReqnPtEta*dReQ1n3k+dImqnPtEta*dImQ1n3k)                      
                      + 2.*dSmPrimePrime1p1kPtEta*dSM1p2k                                      
                      - 6.*dSmPrimePrime1p3kPtEta)        
                      / dM0111; // to be imropoved (notation of dM0111)
   
    // fill the 2D profile to get the average correlation for each (pt, eta) bin:
    if(type == "POI")
    {
     f4pPtEtaPOIW->Fill(fPtMin+(p-1)*dBinWidthPt,fEtaMin+(e-1)*dBinWidthEta,four1n1n1n1nW0W1W1W1PtEta,dM0111);
    }
    else if(type == "RP")
    {
     f4pPtEtaRPW->Fill(fPtMin+(p-1)*dBinWidthPt,fEtaMin+(e-1)*dBinWidthEta,four1n1n1n1nW0W1W1W1PtEta,dM0111);   
    }
   } // end of if(dM0111)
  
  } // end of for(Int_t e=1;e<=fnBinsEta;e++)
 } // end of for(Int_t p=1;p<=fnBinsPt;p++)
  
} // end of AliFlowAnalysisWithQCumulants::CalculateWeightedCorrelationsForDifferentialFlow(TString type)


//================================================================================================================================


void AliFlowAnalysisWithQCumulants::EvaluateNestedLoopsForIntegratedFlow(AliFlowEventSimple* anEvent)
{
 // evaluate the nested loops relevant for integrated flow (needed for cross-checking the results)
 
 Int_t nPrim = anEvent->NumberOfTracks(); 
 
 TH1F *phiWeights = NULL; // histogram with phi weights
 Int_t nBinsPhi = 0; 
 
 if(fUsePhiWeights)
 {
  if(!fWeightsList)
  {
   cout<<" WARNING: fWeightsList is NULL pointer in AFAWQC::ENLFIF(). "<<endl;
   exit(0);
  }
  phiWeights = dynamic_cast<TH1F *>(fWeightsList->FindObject("phi_weights"));
  if(!phiWeights)
  {
   cout<<" WARNING: couldn't access the histogram with phi weights in AFAWQC::ENLFIF(). "<<endl;
   exit(0);
  }
  nBinsPhi = phiWeights->GetNbinsX();
 } 
 
 Double_t phi1=0., phi2=0., phi3=0., phi4=0., phi5=0., phi6=0., phi7=0., phi8=0.;
 Double_t wPhi1=1., wPhi2=1., wPhi3=1., wPhi4=1., wPhi5=1., wPhi6=1., wPhi7=1., wPhi8=1.;
 
 Int_t n=2; // to be improved
 
 //                                          ******************************************
 //                                          **** NESTED LOOPS FOR INTEGRATED FLOW ****
 //                                          ****************************************** 
 //
 // Remark 1: multi-particle correlations calculated with nested loops without weights are stored in 1D profile fDirectCorrelations;
 // Remark 2: multi-particle correlations calculated with nested loops with weights are stored in 1D profile fDirectCorrelationsW;
 
 // Remark 3: binning of fDirectCorrelations is organized as follows:
 // --------------------------------------------------------------------------------------------------------------------
 //  1st bin: <2>_{1n|1n} = two1n1n = cos(n*(phi1-phi2))>
 //  2nd bin: <2>_{2n|2n} = two2n2n = cos(2n*(phi1-phi2))>
 //  3rd bin: <2>_{3n|3n} = two3n3n = cos(3n*(phi1-phi2))> 
 //  4th bin: <2>_{4n|4n} = two4n4n = cos(4n*(phi1-phi2))>
 //  5th bin:           ----  EMPTY ----
 //  6th bin: <3>_{2n|1n,1n} = three2n1n1n = <cos(n*(2.*phi1-phi2-phi3))>
 //  7th bin: <3>_{3n|2n,1n} = three3n2n1n = <cos(n*(3.*phi1-2.*phi2-phi3))>
 //  8th bin: <3>_{4n|2n,2n} = three4n2n2n = <cos(n*(4.*phi1-2.*phi2-2.*phi3))>
 //  9th bin: <3>_{4n|3n,1n} = three4n3n1n = <cos(n*(4.*phi1-3.*phi2-phi3))>
 // 10th bin:           ----  EMPTY ----
 // 11th bin: <4>_{1n,1n|1n,1n} = four1n1n1n1n = <cos(n*(phi1+phi2-phi3-phi4))>
 // 12th bin: <4>_{2n,1n|2n,1n} = four2n1n2n1n = <cos(2.*n*(phi1+phi2-phi3-phi4))>
 // 13th bin: <4>_{2n,2n|2n,2n} = four2n2n2n2n = <cos(n*(2.*phi1+phi2-2.*phi3-phi4))>
 // 14th bin: <4>_{3n|1n,1n,1n} = four3n1n1n1n = <cos(n*(3.*phi1-phi2-phi3-phi4))> 
 // 15th bin: <4>_{3n,1n|3n,1n} = four3n1n3n1n = <cos(n*(4.*phi1-2.*phi2-phi3-phi4))>
 // 16th bin: <4>_{3n,1n|2n,2n} = four3n1n2n2n = <cos(n*(3.*phi1+phi2-2.*phi3-2.*phi4))>
 // 17th bin: <4>_{4n|2n,1n,1n} = four4n2n1n1n = <cos(n*(3.*phi1+phi2-3.*phi3-phi4))> 
 // 18th bin:           ----  EMPTY ----
 // 19th bin: <5>_{2n|1n,1n,1n,1n} = five2n1n1n1n1n = <cos(n*(2.*phi1+phi2-phi3-phi4-phi5))>
 // 20th bin: <5>_{2n,2n|2n,1n,1n} = five2n2n2n1n1n = <cos(n*(2.*phi1+2.*phi2-2.*phi3-phi4-phi5))>
 // 21st bin: <5>_{3n,1n|2n,1n,1n} = five3n1n2n1n1n = <cos(n*(3.*phi1+phi2-2.*phi3-phi4-phi5))>
 // 22nd bin: <5>_{4n|1n,1n,1n,1n} = five4n1n1n1n1n = <cos(n*(4.*phi1-phi2-phi3-phi4-phi5))>
 // 23rd bin:           ----  EMPTY ----
 // 24th bin: <6>_{1n,1n,1n|1n,1n,1n} = six1n1n1n1n1n1n = <cos(n*(phi1+phi2+phi3-phi4-phi5-phi6))>
 // 25th bin: <6>_{2n,1n,1n|2n,1n,1n} = six2n1n1n2n1n1n = <cos(n*(2.*phi1+2.*phi2-phi3-phi4-phi5-phi6))>
 // 26th bin: <6>_{2n,2n|1n,1n,1n,1n} = six2n2n1n1n1n1n = <cos(n*(3.*phi1+phi2-phi3-phi4-phi5-phi6))>
 // 27th bin: <6>_{3n,1n|1n,1n,1n,1n} = six3n1n1n1n1n1n = <cos(n*(2.*phi1+phi2+phi3-2.*phi4-phi5-phi6))>
 // 28th bin:           ----  EMPTY ----
 // 29th bin: <7>_{2n,1n,1n|1n,1n,1n,1n} = seven2n1n1n1n1n1n1n =  <cos(n*(2.*phi1+phi2+phi3-phi4-phi5-phi6-phi7))>
 // 30th bin:           ----  EMPTY ----
 // 31st bin: <8>_{1n,1n,1n,1n|1n,1n,1n,1n} = eight1n1n1n1n1n1n1n1n = <cos(n*(phi1+phi2+phi3+phi4-phi5-phi6-phi7-phi8))>
 // --------------------------------------------------------------------------------------------------------------------
 
 // Remark 4: binning of fDirectCorrelationsW is organized as follows:
 //..............................................................................................
 //       ---- bins 1-20: 2-particle correlations ----
 // 1st bin: two1n1nW1W1 = <w1 w2 cos(n*(phi1-phi2))>
 // 2nd bin: two2n2nW2W2 = <w1^2 w2^2 cos(2n*(phi1-phi2))>
 // 3rd bin: two3n3nW3W3 = <w1^3 w2^3 cos(3n*(phi1-phi2))>
 // 4th bin: two4n4nW4W4 = <w1^4 w2^4 cos(4n*(phi1-phi2))>
 // 5th bin: two1n1nW3W1 = <w1^3 w2 cos(n*(phi1-phi2))>
 // 6th bin: two1n1nW1W1W2 = <w1 w2 w3^2 cos(n*(phi1-phi2))>  
 //       ---- bins 21-40: 3-particle correlations ----
 // 21st bin: three2n1n1nW2W1W1 = <w1^2 w2 w3 cos(n*(2phi1-phi2-phi3))> 
 //       ---- bins 41-60: 4-particle correlations ----
 // 41st bin: four1n1n1n1nW1W1W1W1 = <w1 w2 w3 w4 cos(n*(phi1+phi2-phi3-phi4))>
 //       ---- bins 61-80: 5-particle correlations ---- 
 //       ---- bins 81-100: 6-particle correlations ----
 //       ---- bins 101-120: 7-particle correlations ----
 //       ---- bins 121-140: 8-particle correlations ----
 //..............................................................................................
 
 
 // 2-particle correlations:       
 for(Int_t i1=0;i1<nPrim;i1++)
 {
  fTrack=anEvent->GetTrack(i1);
  if(!(fTrack->InRPSelection())) continue;
  phi1=fTrack->Phi();
  if(phiWeights) wPhi1 = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi1*nBinsPhi/TMath::TwoPi())));
  for(Int_t i2=0;i2<nPrim;i2++)
  {
   if(i2==i1)continue;
   fTrack=anEvent->GetTrack(i2);
   if(!(fTrack->InRPSelection())) continue;
   phi2=fTrack->Phi();
   if(phiWeights) wPhi2 = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi2*nBinsPhi/TMath::TwoPi())));
    
   // non-weighted: 
   //------------------------------------------------------------------------------
   fDirectCorrelations->Fill(0.,cos(n*(phi1-phi2)),1.);    // <cos(n*(phi1-phi2))>
   fDirectCorrelations->Fill(1.,cos(2.*n*(phi1-phi2)),1.); // <cos(2n*(phi1-phi2))>
   fDirectCorrelations->Fill(2.,cos(3.*n*(phi1-phi2)),1.); // <cos(3n*(phi1-phi2))>
   fDirectCorrelations->Fill(3.,cos(4.*n*(phi1-phi2)),1.); // <cos(4n*(phi1-phi2))> 
   //------------------------------------------------------------------------------
   
   // weighted:
   //................................................................................................................
   fDirectCorrelationsW->Fill(0.,cos(n*(phi1-phi2)),wPhi1*wPhi2);                  // <w1   w2   cos( n*(phi1-phi2))>
   fDirectCorrelationsW->Fill(1.,cos(2.*n*(phi1-phi2)),pow(wPhi1,2)*pow(wPhi2,2)); // <w1^2 w2^2 cos(2n*(phi1-phi2))>
   fDirectCorrelationsW->Fill(2.,cos(3.*n*(phi1-phi2)),pow(wPhi1,3)*pow(wPhi2,3)); // <w1^3 w2^3 cos(3n*(phi1-phi2))>
   fDirectCorrelationsW->Fill(3.,cos(4.*n*(phi1-phi2)),pow(wPhi1,4)*pow(wPhi2,4)); // <w1^4 w2^4 cos(4n*(phi1-phi2))> 
   fDirectCorrelationsW->Fill(4.,cos(n*(phi1-phi2)),pow(wPhi1,3)*wPhi2);           // <w1^3 w2 cos(n*(phi1-phi2))>
   //................................................................................................................
   
  }
 }  

 // 3-particle correlations:         
 for(Int_t i1=0;i1<nPrim;i1++)
 {
  fTrack=anEvent->GetTrack(i1);
  if(!(fTrack->InRPSelection())) continue;
  phi1=fTrack->Phi();
  if(phiWeights) wPhi1 = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi1*nBinsPhi/TMath::TwoPi())));
  for(Int_t i2=0;i2<nPrim;i2++)
  {
   if(i2==i1)continue;
   fTrack=anEvent->GetTrack(i2);
   if(!(fTrack->InRPSelection())) continue;
   phi2=fTrack->Phi();
   if(phiWeights) wPhi2 = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi2*nBinsPhi/TMath::TwoPi())));
   for(Int_t i3=0;i3<nPrim;i3++)
   {
    if(i3==i1||i3==i2)continue;
    fTrack=anEvent->GetTrack(i3);
    if(!(fTrack->InRPSelection())) continue;
    phi3=fTrack->Phi();
    if(phiWeights) wPhi3 = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi3*nBinsPhi/TMath::TwoPi())));
    
    // non-weighted:
    //-----------------------------------------------------------------------------------
    fDirectCorrelations->Fill(5.,cos(2.*n*phi1-n*(phi2+phi3)),1.);       //<3>_{2n|nn,n}
    fDirectCorrelations->Fill(6.,cos(3.*n*phi1-2.*n*phi2-n*phi3),1.);    //<3>_{3n|2n,n}
    fDirectCorrelations->Fill(7.,cos(4.*n*phi1-2.*n*phi2-2.*n*phi3),1.); //<3>_{4n|2n,2n}
    fDirectCorrelations->Fill(8.,cos(4.*n*phi1-3.*n*phi2-n*phi3),1.);    //<3>_{4n|3n,n}
    //-----------------------------------------------------------------------------------
    
    // weighted:
    //..............................................................................................................................
    // 2-p:
    fDirectCorrelationsW->Fill(5.,cos(n*(phi1-phi2)),wPhi1*wPhi2*pow(wPhi3,2)); // <w1 w2 w3^2 cos(n*(phi1-phi2))>
    
    // 3-p:
    fDirectCorrelationsW->Fill(20.,cos(2.*n*phi1-n*(phi2+phi3)),pow(wPhi1,2)*wPhi2*wPhi3); // <w1^2 w2 w3 cos(n*(2phi1-phi2-phi3))>
    //..............................................................................................................................
    
   }
  }
 }

 // 4-particle correlations:       
 for(Int_t i1=0;i1<nPrim;i1++)
 {
  fTrack=anEvent->GetTrack(i1);
  if(!(fTrack->InRPSelection())) continue;
  phi1=fTrack->Phi();
  if(phiWeights) wPhi1 = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi1*nBinsPhi/TMath::TwoPi())));
  for(Int_t i2=0;i2<nPrim;i2++)
  {
   if(i2==i1)continue;
   fTrack=anEvent->GetTrack(i2);
   if(!(fTrack->InRPSelection())) continue;
   phi2=fTrack->Phi();
   if(phiWeights) wPhi2 = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi2*nBinsPhi/TMath::TwoPi())));
   for(Int_t i3=0;i3<nPrim;i3++)
   {
    if(i3==i1||i3==i2)continue;
    fTrack=anEvent->GetTrack(i3);
    if(!(fTrack->InRPSelection())) continue;
    phi3=fTrack->Phi();
    if(phiWeights) wPhi3 = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi3*nBinsPhi/TMath::TwoPi())));
    for(Int_t i4=0;i4<nPrim;i4++)
    {
     if(i4==i1||i4==i2||i4==i3)continue;
     fTrack=anEvent->GetTrack(i4);
     if(!(fTrack->InRPSelection())) continue;
     phi4=fTrack->Phi();
     if(phiWeights) wPhi4 = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi4*nBinsPhi/TMath::TwoPi())));
     
     // non-weighted:
     //-------------------------------------------------------------------------------------------------
     fDirectCorrelations->Fill(10.,cos(n*phi1+n*phi2-n*phi3-n*phi4),1.);            // <4>_{n,n|n,n} 
     fDirectCorrelations->Fill(11.,cos(2.*n*phi1+n*phi2-2.*n*phi3-n*phi4),1.);      // <4>_{2n,n|2n,n}
     fDirectCorrelations->Fill(12.,cos(2.*n*phi1+2*n*phi2-2.*n*phi3-2.*n*phi4),1.); // <4>_{2n,2n|2n,2n}
     fDirectCorrelations->Fill(13.,cos(3.*n*phi1-n*phi2-n*phi3-n*phi4),1.);         // <4>_{3n|n,n,n}
     fDirectCorrelations->Fill(14.,cos(3.*n*phi1+n*phi2-3.*n*phi3-n*phi4),1.);      // <4>_{3n,n|3n,n}   
     fDirectCorrelations->Fill(15.,cos(3.*n*phi1+n*phi2-2.*n*phi3-2.*n*phi4),1.);   // <4>_{3n,n|2n,2n}
     fDirectCorrelations->Fill(16.,cos(4.*n*phi1-2.*n*phi2-n*phi3-n*phi4),1.);      // <4>_{4n|2n,n,n}
     //-------------------------------------------------------------------------------------------------
     
     // weighted:
     //.......................................................................................
     // 4-p:
     fDirectCorrelationsW->Fill(40.,cos(n*phi1+n*phi2-n*phi3-n*phi4),wPhi1*wPhi2*wPhi3*wPhi4);              
     //.......................................................................................
     
    }  
   }
  }
 }

 // 5-particle correlations:      
 for(Int_t i1=0;i1<nPrim;i1++)
 {
  //cout<<"i1 = "<<i1<<endl;
  fTrack=anEvent->GetTrack(i1);
  if(!(fTrack->InRPSelection())) continue;  
  phi1=fTrack->Phi();
  if(phiWeights) wPhi1 = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi1*nBinsPhi/TMath::TwoPi())));
  for(Int_t i2=0;i2<nPrim;i2++)
  {
   if(i2==i1)continue;
   fTrack=anEvent->GetTrack(i2);
   if(!(fTrack->InRPSelection())) continue;
   phi2=fTrack->Phi();
   if(phiWeights) wPhi2 = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi2*nBinsPhi/TMath::TwoPi())));
   for(Int_t i3=0;i3<nPrim;i3++)
   {
    if(i3==i1||i3==i2)continue;
    fTrack=anEvent->GetTrack(i3);
    if(!(fTrack->InRPSelection())) continue;
    phi3=fTrack->Phi();
    if(phiWeights) wPhi3 = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi3*nBinsPhi/TMath::TwoPi())));
    for(Int_t i4=0;i4<nPrim;i4++)
    {
     if(i4==i1||i4==i2||i4==i3)continue;
     fTrack=anEvent->GetTrack(i4);
     if(!(fTrack->InRPSelection())) continue;
     phi4=fTrack->Phi();
     if(phiWeights) wPhi4 = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi4*nBinsPhi/TMath::TwoPi())));
     for(Int_t i5=0;i5<nPrim;i5++)
     {
      if(i5==i1||i5==i2||i5==i3||i5==i4)continue;
      fTrack=anEvent->GetTrack(i5);
      if(!(fTrack->InRPSelection())) continue;
      phi5=fTrack->Phi();
      if(phiWeights) wPhi5 = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi5*nBinsPhi/TMath::TwoPi())));
      
      // non-weighted:
      //------------------------------------------------------------------------------------------------------
      fDirectCorrelations->Fill(18.,cos(2.*n*phi1+n*phi2-n*phi3-n*phi4-n*phi5),1.);       //<5>_{2n,n|n,n,n}
      fDirectCorrelations->Fill(19.,cos(2.*n*phi1+2.*n*phi2-2.*n*phi3-n*phi4-n*phi5),1.); //<5>_{2n,2n|2n,n,n}
      fDirectCorrelations->Fill(20.,cos(3.*n*phi1+n*phi2-2.*n*phi3-n*phi4-n*phi5),1.);    //<5>_{3n,n|2n,n,n}
      fDirectCorrelations->Fill(21.,cos(4.*n*phi1-n*phi2-n*phi3-n*phi4-n*phi5),1.);       //<5>_{4n|n,n,n,n}
      //------------------------------------------------------------------------------------------------------
      
      // weighted:
      //..............................................................................................................
      // 5-p:
      fDirectCorrelationsW->Fill(60.,cos(2.*n*phi1+n*phi2-n*phi3-n*phi4-n*phi5),pow(wPhi1,2)*wPhi2*wPhi3*wPhi4*wPhi5);     
      //..............................................................................................................
      
     }
    }  
   }
  }
 }
 
 // 6-particle correlations:
 for(Int_t i1=0;i1<nPrim;i1++)
 {
  //cout<<"i1 = "<<i1<<endl;
  fTrack=anEvent->GetTrack(i1);
  if(!(fTrack->InRPSelection())) continue;
  phi1=fTrack->Phi();
  if(phiWeights) wPhi1 = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi1*nBinsPhi/TMath::TwoPi())));
  for(Int_t i2=0;i2<nPrim;i2++)
  {
   if(i2==i1)continue;
   fTrack=anEvent->GetTrack(i2);
   if(!(fTrack->InRPSelection())) continue;
   phi2=fTrack->Phi();
   if(phiWeights) wPhi2 = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi2*nBinsPhi/TMath::TwoPi())));
   for(Int_t i3=0;i3<nPrim;i3++)
   {
    if(i3==i1||i3==i2)continue;
    fTrack=anEvent->GetTrack(i3);
    if(!(fTrack->InRPSelection())) continue;
    phi3=fTrack->Phi();
    if(phiWeights) wPhi3 = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi3*nBinsPhi/TMath::TwoPi())));
    for(Int_t i4=0;i4<nPrim;i4++)
    {
     if(i4==i1||i4==i2||i4==i3)continue;
     fTrack=anEvent->GetTrack(i4);
     if(!(fTrack->InRPSelection())) continue;
     phi4=fTrack->Phi();
     if(phiWeights) wPhi4 = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi4*nBinsPhi/TMath::TwoPi())));
     for(Int_t i5=0;i5<nPrim;i5++)
     {
      if(i5==i1||i5==i2||i5==i3||i5==i4)continue;
      fTrack=anEvent->GetTrack(i5);
      if(!(fTrack->InRPSelection())) continue;
      phi5=fTrack->Phi();
      if(phiWeights) wPhi5 = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi5*nBinsPhi/TMath::TwoPi())));
      for(Int_t i6=0;i6<nPrim;i6++)
      {
       if(i6==i1||i6==i2||i6==i3||i6==i4||i6==i5)continue;
       fTrack=anEvent->GetTrack(i6);
       if(!(fTrack->InRPSelection())) continue;
       phi6=fTrack->Phi(); 
       if(phiWeights) wPhi6 = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi6*nBinsPhi/TMath::TwoPi())));
       
       // non-weighted:
       //-----------------------------------------------------------------------------------------------------------
       fDirectCorrelations->Fill(23.,cos(n*phi1+n*phi2+n*phi3-n*phi4-n*phi5-n*phi6),1.);       //<6>_{n,n,n|n,n,n}
       fDirectCorrelations->Fill(24.,cos(2.*n*phi1+n*phi2+n*phi3-2.*n*phi4-n*phi5-n*phi6),1.); //<6>_{2n,n,n|2n,n,n}
       fDirectCorrelations->Fill(25.,cos(2.*n*phi1+2.*n*phi2-n*phi3-n*phi4-n*phi5-n*phi6),1.); //<6>_{2n,2n|n,n,n,n}
       fDirectCorrelations->Fill(26.,cos(3.*n*phi1+n*phi2-n*phi3-n*phi4-n*phi5-n*phi6),1.);    //<6>_{3n,n|n,n,n,n}  
       //-----------------------------------------------------------------------------------------------------------

       // weighted:
       //.................................................................................................................
       // 6-p:
       fDirectCorrelationsW->Fill(80.,cos(n*phi1+n*phi2+n*phi3-n*phi4-n*phi5-n*phi6),wPhi1*wPhi2*wPhi3*wPhi4*wPhi5*wPhi6);
       //.................................................................................................................       
          
      } 
     }
    }  
   }
  }
 }
 
 // 7-particle correlations:
 for(Int_t i1=0;i1<nPrim;i1++)
 {
  //cout<<"i1 = "<<i1<<endl;
  fTrack=anEvent->GetTrack(i1);
  if(!(fTrack->InRPSelection())) continue;
  phi1=fTrack->Phi();
  if(phiWeights) wPhi1 = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi1*nBinsPhi/TMath::TwoPi())));
  for(Int_t i2=0;i2<nPrim;i2++)
  {
   if(i2==i1)continue;
   fTrack=anEvent->GetTrack(i2);
   if(!(fTrack->InRPSelection())) continue;
   phi2=fTrack->Phi();
   if(phiWeights) wPhi2 = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi2*nBinsPhi/TMath::TwoPi())));
   for(Int_t i3=0;i3<nPrim;i3++)
   {
    if(i3==i1||i3==i2)continue;
    fTrack=anEvent->GetTrack(i3);
    if(!(fTrack->InRPSelection())) continue;
    phi3=fTrack->Phi();
    if(phiWeights) wPhi3 = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi3*nBinsPhi/TMath::TwoPi())));
    for(Int_t i4=0;i4<nPrim;i4++)
    {
     if(i4==i1||i4==i2||i4==i3)continue;
     fTrack=anEvent->GetTrack(i4);
     if(!(fTrack->InRPSelection())) continue;
     phi4=fTrack->Phi();
     if(phiWeights) wPhi4 = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi4*nBinsPhi/TMath::TwoPi())));
     for(Int_t i5=0;i5<nPrim;i5++)
     {
      if(i5==i1||i5==i2||i5==i3||i5==i4)continue;
      fTrack=anEvent->GetTrack(i5);
      if(!(fTrack->InRPSelection())) continue;
      phi5=fTrack->Phi();
      if(phiWeights) wPhi5 = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi5*nBinsPhi/TMath::TwoPi())));
      for(Int_t i6=0;i6<nPrim;i6++)
      {
       if(i6==i1||i6==i2||i6==i3||i6==i4||i6==i5)continue;
       fTrack=anEvent->GetTrack(i6);
       if(!(fTrack->InRPSelection())) continue;
       phi6=fTrack->Phi(); 
       if(phiWeights) wPhi6 = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi6*nBinsPhi/TMath::TwoPi())));
       for(Int_t i7=0;i7<nPrim;i7++)
       {
        if(i7==i1||i7==i2||i7==i3||i7==i4||i7==i5||i7==i6)continue;
        fTrack=anEvent->GetTrack(i7);
        if(!(fTrack->InRPSelection())) continue;
        phi7=fTrack->Phi(); 
        if(phiWeights) wPhi7 = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi7*nBinsPhi/TMath::TwoPi())));
        
        // non-weighted:
        //---------------------------------------------------------------------------------------------------------------
        fDirectCorrelations->Fill(28.,cos(2.*n*phi1+n*phi2+n*phi3-n*phi4-n*phi5-n*phi6-n*phi7),1.);//<7>_{2n,n,n|n,n,n,n}
        //---------------------------------------------------------------------------------------------------------------
        
        // weighted:
        //..........................................................................................................................................
        fDirectCorrelationsW->Fill(100.,cos(2.*n*phi1+n*phi2+n*phi3-n*phi4-n*phi5-n*phi6-n*phi7),pow(wPhi1,2.)*wPhi2*wPhi3*wPhi4*wPhi5*wPhi6*wPhi7);
        //..........................................................................................................................................
        
       } 
      } 
     }
    }  
   }
  }
 }
 
 // 8-particle correlations:
 for(Int_t i1=0;i1<nPrim;i1++)
 {
  cout<<"i1 = "<<i1<<endl;
  fTrack=anEvent->GetTrack(i1);
  if(!(fTrack->InRPSelection())) continue;
  phi1=fTrack->Phi();
  if(phiWeights) wPhi1 = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi1*nBinsPhi/TMath::TwoPi())));
  for(Int_t i2=0;i2<nPrim;i2++)
  {
   if(i2==i1)continue;
   fTrack=anEvent->GetTrack(i2);
   if(!(fTrack->InRPSelection())) continue;
   phi2=fTrack->Phi();
   if(phiWeights) wPhi2 = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi2*nBinsPhi/TMath::TwoPi())));
   for(Int_t i3=0;i3<nPrim;i3++)
   {
    if(i3==i1||i3==i2)continue;
    fTrack=anEvent->GetTrack(i3);
    if(!(fTrack->InRPSelection())) continue;
    phi3=fTrack->Phi();
    if(phiWeights) wPhi3 = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi3*nBinsPhi/TMath::TwoPi())));
    for(Int_t i4=0;i4<nPrim;i4++)
    {
     if(i4==i1||i4==i2||i4==i3)continue;
     fTrack=anEvent->GetTrack(i4);
     if(!(fTrack->InRPSelection())) continue;
     phi4=fTrack->Phi();
     if(phiWeights) wPhi4 = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi4*nBinsPhi/TMath::TwoPi())));
     for(Int_t i5=0;i5<nPrim;i5++)
     {
      if(i5==i1||i5==i2||i5==i3||i5==i4)continue;
      fTrack=anEvent->GetTrack(i5);
      if(!(fTrack->InRPSelection())) continue;
      phi5=fTrack->Phi();
      if(phiWeights) wPhi5 = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi5*nBinsPhi/TMath::TwoPi())));
      for(Int_t i6=0;i6<nPrim;i6++)
      {
       if(i6==i1||i6==i2||i6==i3||i6==i4||i6==i5)continue;
       fTrack=anEvent->GetTrack(i6);
       if(!(fTrack->InRPSelection())) continue;
       phi6=fTrack->Phi();
       if(phiWeights) wPhi6 = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi6*nBinsPhi/TMath::TwoPi()))); 
       for(Int_t i7=0;i7<nPrim;i7++)
       {
        if(i7==i1||i7==i2||i7==i3||i7==i4||i7==i5||i7==i6)continue;
        fTrack=anEvent->GetTrack(i7);
        if(!(fTrack->InRPSelection())) continue;
        phi7=fTrack->Phi();
        if(phiWeights) wPhi7 = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi7*nBinsPhi/TMath::TwoPi()))); 
        for(Int_t i8=0;i8<nPrim;i8++)
        {
         if(i8==i1||i8==i2||i8==i3||i8==i4||i8==i5||i8==i6||i8==i7)continue;
         fTrack=anEvent->GetTrack(i8);
         if(!(fTrack->InRPSelection())) continue;
         phi8=fTrack->Phi();
         if(phiWeights) wPhi8 = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi8*nBinsPhi/TMath::TwoPi()))); 
          
         // non-weighted: 
         //--------------------------------------------------------------------------------------------------------------------
         fDirectCorrelations->Fill(30.,cos(n*phi1+n*phi2+n*phi3+n*phi4-n*phi5-n*phi6-n*phi7-n*phi8),1.);//<8>_{n,n,n,n|n,n,n,n}
         //--------------------------------------------------------------------------------------------------------------------
         
         // weighted: 
         //...........................................................................................................................................
         fDirectCorrelations->Fill(120.,cos(n*phi1+n*phi2+n*phi3+n*phi4-n*phi5-n*phi6-n*phi7-n*phi8),wPhi1*wPhi2*wPhi3*wPhi4*wPhi5*wPhi6*wPhi7*wPhi8);
         //...........................................................................................................................................
     
        } 
       } 
      } 
     }
    }  
   }
  }
 }
 
} // end of AliFlowAnalysisWithQCumulants::EvaluateNestedLoopsForIntegratedFlow(AliFlowEventSimple* anEvent)


//================================================================================================================================


void AliFlowAnalysisWithQCumulants::EvaluateNestedLoopsForDifferentialFlow(AliFlowEventSimple* anEvent)
{
 // evaluate the nested loops relevant for differential flow (needed for cross-checking the results)
 
 Int_t nPrim = anEvent->NumberOfTracks(); 

 TH1F *phiWeights = NULL; // histogram with phi weights
 Int_t nBinsPhi = 0; 
 
 if(fUsePhiWeights)
 {
  if(!fWeightsList)
  {
   cout<<" WARNING: fWeightsList is NULL pointer in AFAWQC::ENLFDF(). "<<endl;
   exit(0);
  }
  phiWeights = dynamic_cast<TH1F *>(fWeightsList->FindObject("phi_weights"));
  if(!phiWeights)
  {
   cout<<" WARNING: couldn't access the histogram with phi weights in AFAWQC::ENLFDF(). "<<endl;
   exit(0);
  }
  nBinsPhi = phiWeights->GetNbinsX();
 } 
 
 Double_t psi1=0., phi2=0., phi3=0., phi4=0.;// phi5=0., phi6=0., phi7=0., phi8=0.;
 Double_t wPhi1=1., wPhi2=1., wPhi3=1., wPhi4=1.;// wPhi5=1., wPhi6=1., wPhi7=1., wPhi8=1.;
 
 Int_t n=2; // to be improved
 
 //                                          ********************************************
 //                                          **** NESTED LOOPS FOR DIFFERENTIAL FLOW ****
 //                                          ******************************************** 
 
 // Remark 1: (pt,eta) bin in which the cross-checking will be performed is given by 1.1 < pt < 1.2 GeV and -0.55 < eta < -0.525 
 
 // Remark 2: multi-particle correlations needed for diff. flow calculated with nested loops without weights are stored in 1D profile  
 //           fDirectCorrelationsDiffFlow
 
 // Remark 3: multi-particle correlations needed for diff. flow calculated with nested loops with weights are stored in 1D profile  
 //           fDirectCorrelationsDiffFlowW;
 
 // Remark 4: binning of fDirectCorrelationsDiffFlow is organized as follows:
 //......................................................................................
 //       ---- bins 1-20: 2-particle correlations ----
 //  1st bin: <2'>_{1n|1n} = twoPrime1n1n = <cos(n*(psi1-phi2))>
 //       ---- bins 21-40: 3-particle correlations ----
 //       ---- bins 41-60: 4-particle correlations ----
 // 41st bin: <4'>_{1n,1n|1n,1n} = fourPrime1n1n1n1n  = <cos(n*(psi1+phi2-phi3-phi4))>
 //......................................................................................
 
 // Remark 5: binning of fDirectCorrelationsDiffFlow is organized as follows:
 //......................................................................................
 //       ---- bins 1-20: 2-particle correlations ----
 //  1st bin: twoPrime1n1nW0W1 = <w2 cos(n*(psi1-phi2))>
 //       ---- bins 21-40: 3-particle correlations ----
 //       ---- bins 41-60: 4-particle correlations ----
 // 41st bin: fourPrime1n1n1n1nW0W1W1W1 = <w2 w3 w4 cos(n*(psi1+phi2-phi3-phi4))>
 //......................................................................................
 
 // 2'-particle:
 for(Int_t i1=0;i1<nPrim;i1++)
 {
  fTrack=anEvent->GetTrack(i1);
  // POI condition (first particle in the correlator must be POI): 
  if(!((fTrack->Pt()>=1.1 && fTrack->Pt()<1.2) && (fTrack->Eta()>=-0.55 && fTrack->Eta()<-0.525) && (fTrack->InPOISelection())))continue; 
  psi1=fTrack->Phi();
  if(phiWeights) wPhi1 = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(psi1*nBinsPhi/TMath::TwoPi())));
  for(Int_t i2=0;i2<nPrim;i2++)
  {
   if(i2==i1)continue;
   fTrack=anEvent->GetTrack(i2);
   // RP condition (!(first) particle in the correlator must be RP):
   if(!(fTrack->InRPSelection()))continue;
   phi2=fTrack->Phi();
   if(phiWeights) wPhi2 = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi2*nBinsPhi/TMath::TwoPi())));
    
   // non-weighted: 
   //.....................................................................................  
   fDirectCorrelationsDiffFlow->Fill(0.,cos(1.*n*(psi1-phi2)),1.); // <cos(n*(psi1-phi2))>
   //.....................................................................................  
   
   // weighted:
   //.....................................................................................   
   fDirectCorrelationsDiffFlowW->Fill(0.,cos(1.*n*(psi1-phi2)),wPhi2); // <w2 cos(n*(psi1-phi2))>
   //.....................................................................................  
   
   //fDirectCorrelations->Fill(103.,cos(1.*n*(phi1-phi2)),pow(wPhi1,2)*wPhi2);//<2'>_{n,n}
   //fDirectCorrelations->Fill(104.,cos(2.*n*(phi1-phi2)),wPhi1*pow(wPhi2,2));//<2'>_{n,n}
   //fDirectCorrelations->Fill(105.,cos(1.*n*(phi1-phi2)),pow(wPhi2,3));//<2'>_{n,n}  
   //fDirectCorrelations->Fill(41.,cos(2.*n*(phi1-phi2)),1);//<2'>_{2n,2n}
   //fDirectCorrelations->Fill(42.,cos(3.*n*(phi1-phi2)),1);//<2'>_{3n,3n}
   //fDirectCorrelations->Fill(43.,cos(4.*n*(phi1-phi2)),1);//<2'>_{4n,4n}   
    
  }//end of for(Int_t i2=0;i2<nPrim;i2++)
 }//end of for(Int_t i1=0;i1<nPrim;i1++)
 
 
 
 /*
 
 //<3'>_{2n|n,n}
 for(Int_t i1=0;i1<nPrim;i1++)
 {
  fTrack=anEvent->GetTrack(i1);
  if(!((fTrack->Pt()>=0.5&&fTrack->Pt()<0.6)&&(fTrack->InPOISelection())))continue;//POI condition
  psi1=fTrack->Phi();
  if(phiWeights) wPhi1 = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(psi1*nBinsPhi/TMath::TwoPi())));
  for(Int_t i2=0;i2<nPrim;i2++)
  {
   if(i2==i1)continue;
   fTrack=anEvent->GetTrack(i2);
   if(!(fTrack->InRPSelection()))continue;//RP condition
   phi2=fTrack->Phi();
   if(phiWeights) wPhi2 = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi2*nBinsPhi/TMath::TwoPi())));
   for(Int_t i3=0;i3<nPrim;i3++)
   {
    if(i3==i1||i3==i2)continue;
    fTrack=anEvent->GetTrack(i3);
    if(!(fTrack->InRPSelection()))continue;//RP condition
    phi3=fTrack->Phi();
    if(phiWeights) wPhi3 = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi3*nBinsPhi/TMath::TwoPi())));
    //fill the fDirectCorrelations:     
    
    // 2-p
    //fDirectCorrelations->Fill(101.,cos(n*(phi2-phi3)),wPhi1*wPhi2*wPhi3); // <w1 w2 w3 cos(n(phi2-phi3))>
    //fDirectCorrelations->Fill(102.,cos(n*(phi1-phi3)),pow(wPhi2,2.)*wPhi3); // <w2^2 w3 cos(n(psi1-phi2))>
    
    // 3-p            
    //fDirectCorrelations->Fill(110.,cos(n*(2.*phi1-phi2-phi3)),wPhi1*wPhi2*wPhi3); // <w1 w2 w3 cos(n(2psi1-phi2-phi3))>
    //fDirectCorrelations->Fill(111.,cos(n*(phi1+phi2-2.*phi3)),wPhi2*pow(wPhi3,2.)); // <w2 w3^2 cos(n(psi1+phi2-2.*phi3))>
    
    
    //fDirectCorrelations->Fill(46.,cos(n*(phi1+phi2-2.*phi3)),1);//<3'>_{n,n|2n}    
   }//end of for(Int_t i3=0;i3<nPrim;i3++)  
  }//end of for(Int_t i2=0;i2<nPrim;i2++)  
 }//end of for(Int_t i1=0;i1<nPrim;i1++)
 */
 
 // 4'-particle:
 for(Int_t i1=0;i1<nPrim;i1++)
 {
  fTrack=anEvent->GetTrack(i1);
  // POI condition (first particle in the correlator must be POI): 
  if(!((fTrack->Pt()>=1.1 && fTrack->Pt()<1.2) && (fTrack->Eta()>=-0.55 && fTrack->Eta()<-0.525) && (fTrack->InPOISelection())))continue; 
  psi1=fTrack->Phi();
  if(phiWeights) wPhi1 = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(psi1*nBinsPhi/TMath::TwoPi())));
  for(Int_t i2=0;i2<nPrim;i2++)
  {
   if(i2==i1)continue;
   fTrack=anEvent->GetTrack(i2);
   // RP condition (!(first) particle in the correlator must be RP): 
   if(!(fTrack->InRPSelection()))continue;
   phi2=fTrack->Phi();
   if(phiWeights) wPhi2 = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi2*nBinsPhi/TMath::TwoPi())));
   for(Int_t i3=0;i3<nPrim;i3++)
   { 
    if(i3==i1||i3==i2)continue;
    fTrack=anEvent->GetTrack(i3);
    // RP condition (!(first) particle in the correlator must be RP):
    if(!(fTrack->InRPSelection()))continue;
    phi3=fTrack->Phi();
    if(phiWeights) wPhi3 = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi3*nBinsPhi/TMath::TwoPi())));
    for(Int_t i4=0;i4<nPrim;i4++)
    {
     if(i4==i1||i4==i2||i4==i3)continue;
     fTrack=anEvent->GetTrack(i4);
     // RP condition (!(first) particle in the correlator must be RP):
     if(!(fTrack->InRPSelection()))continue;  
     phi4=fTrack->Phi();
     if(phiWeights) wPhi4 = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi4*nBinsPhi/TMath::TwoPi())));
     
     // non-weighted:
     //.........................................................................................................................
     fDirectCorrelationsDiffFlow->Fill(40.,cos(n*(psi1+phi2-phi3-phi4)),1.); // <cos(n(psi1+phi1-phi2-phi3))> 
     //.........................................................................................................................     

     // weighted:
     //...............................................................................................................................
     fDirectCorrelationsDiffFlowW->Fill(40.,cos(n*(psi1+phi2-phi3-phi4)),wPhi2*wPhi3*wPhi4); // <w2 w3 w4 cos(n(psi1+phi2-phi3-phi4))> 
     //...............................................................................................................................     
          
    }//end of for(Int_t i4=0;i4<nPrim;i4++)
   }//end of for(Int_t i3=0;i3<nPrim;i3++)
  }//end of for(Int_t i2=0;i2<nPrim;i2++) 
 }//end of for(Int_t i1=0;i1<nPrim;i1++)
 
 /*                
 //<5'>_{2n,n|n,n,n}
 for(Int_t i1=0;i1<nPrim;i1++)
 {
  fTrack=anEvent->GetTrack(i1);
  if(!((fTrack->Pt()>=0.5&&fTrack->Pt()<0.6)&&(fTrack->InPOISelection())))continue;//POI condition
  phi1=fTrack->Phi();
  for(Int_t i2=0;i2<nPrim;i2++)
  {
   if(i2==i1)continue;
   fTrack=anEvent->GetTrack(i2);
   if(!(fTrack->InRPSelection()))continue;//RP condition   
   phi2=fTrack->Phi();
   for(Int_t i3=0;i3<nPrim;i3++)
   { 
    if(i3==i1||i3==i2)continue;
    fTrack=anEvent->GetTrack(i3);
    if(!(fTrack->InRPSelection()))continue;//RP condition   
    phi3=fTrack->Phi();
    for(Int_t i4=0;i4<nPrim;i4++)
    {
     if(i4==i1||i4==i2||i4==i3)continue;
     fTrack=anEvent->GetTrack(i4);
     if(!(fTrack->InRPSelection()))continue;//RP condition  
     phi4=fTrack->Phi();//
     for(Int_t i5=0;i5<nPrim;i5++)
     {
      if(i5==i1||i5==i2||i5==i3||i5==i4)continue;
      fTrack=anEvent->GetTrack(i5);
      if(!(fTrack->InRPSelection()))continue;//RP condition  
      phi5=fTrack->Phi();    
      //fill the fDirectCorrelations:if(bNestedLoops)
      //fDirectCorrelations->Fill(55.,cos(2.*n*phi1+n*phi2-n*phi3-n*phi4-n*phi5),1);//<5'>_{2n,n|n,n,n}
     }//end of for(Int_t i5=0;i5<nPrim;i5++)  
    }//end of for(Int_t i4=0;i4<nPrim;i4++)
   }//end of for(Int_t i3=0;i3<nPrim;i3++)
  }//end of for(Int_t i2=0;i2<nPrim;i2++) 
 }//end of for(Int_t i1=0;i1<nPrim;i1++)
 

  
 */
 /*
 
 
 
 //<6'>_{n,n,n|n,n,n}
 for(Int_t i1=0;i1<nPrim;i1++)
 {
  fTrack=anEvent->GetTrack(i1);
  if(!((fTrack->Pt()>=0.5&&fTrack->Pt()<0.6)&&(fTrack->InPOISelection())))continue;//POI condition
  phi1=fTrack->Phi();
  for(Int_t i2=0;i2<nPrim;i2++)
  {
   if(i2==i1)continue;
   fTrack=anEvent->GetTrack(i2);
   if(!(fTrack->InRPSelection()))continue;//RP condition   
   phi2=fTrack->Phi();
   for(Int_t i3=0;i3<nPrim;i3++)
   { 
    if(i3==i1||i3==i2)continue;
    fTrack=anEvent->GetTrack(i3);
    if(!(fTrack->InRPSelection()))continue;//RP condition   
    phi3=fTrack->Phi();
    for(Int_t i4=0;i4<nPrim;i4++)
    {
     if(i4==i1||i4==i2||i4==i3)continue;
     fTrack=anEvent->GetTrack(i4);
     if(!(fTrack->InRPSelection()))continue;//RP condition  
     phi4=fTrack->Phi();
     for(Int_t i5=0;i5<nPrim;i5++)
     {
      if(i5==i1||i5==i2||i5==i3||i5==i4)continue;
      fTrack=anEvent->GetTrack(i5);
      if(!(fTrack->InRPSelection()))continue;//RP condition  
      phi5=fTrack->Phi();    
      for(Int_t i6=0;i6<nPrim;i6++)
      {
       if(i6==i1||i6==i2||i6==i3||i6==i4||i6==i5)continue;
       fTrack=anEvent->GetTrack(i6);
       if(!(fTrack->InRPSelection()))continue;//RP condition  
       phi6=fTrack->Phi();  
       //fill the fDirectCorrelations:
       //fDirectCorrelations->Fill(60.,cos(n*(phi1+phi2+phi3-phi4-phi5-phi6)),1);//<6'>_{n,n,n|n,n,n}
      }//end of for(Int_t i6=0;i6<nPrim;i6++)   
     }//end of for(Int_t i5=0;i5<nPrim;i5++)  
    }//end of for(Int_t i4=0;i4<nPrim;i4++)
   }//end of for(Int_t i3=0;i3<nPrim;i3++)
  }//end of for(Int_t i2=0;i2<nPrim;i2++) 
 }//end of for(Int_t i1=0;i1<nPrim;i1++)

 
 */
 /* 
   
     
 //<7'>_{2n,n,n|n,n,n,n}
 for(Int_t i1=0;i1<nPrim;i1++)
 {
  fTrack=anEvent->GetTrack(i1);
  if(!((fTrack->Pt()>=0.5&&fTrack->Pt()<0.6)&&(fTrack->InPOISelection())))continue;//POI condition
  phi1=fTrack->Phi();
  for(Int_t i2=0;i2<nPrim;i2++)
  {
   if(i2==i1)continue;
   fTrack=anEvent->GetTrack(i2);
   if(!(fTrack->InRPSelection()))continue;//RP condition   
   phi2=fTrack->Phi();
   for(Int_t i3=0;i3<nPrim;i3++)
   { 
    if(i3==i1||i3==i2)continue;
    fTrack=anEvent->GetTrack(i3);
    if(!(fTrack->InRPSelection()))continue;//RP condition   
    phi3=fTrack->Phi();
    for(Int_t i4=0;i4<nPrim;i4++)
    {
     if(i4==i1||i4==i2||i4==i3)continue;
     fTrack=anEvent->GetTrack(i4);
     if(!(fTrack->InRPSelection()))continue;//RP condition  
     phi4=fTrack->Phi();
     for(Int_t i5=0;i5<nPrim;i5++)
     {
      if(i5==i1||i5==i2||i5==i3||i5==i4)continue;
      fTrack=anEvent->GetTrack(i5);
      if(!(fTrack->InRPSelection()))continue;//RP condition  
      phi5=fTrack->Phi();    
      for(Int_t i6=0;i6<nPrim;i6++)
      {
       if(i6==i1||i6==i2||i6==i3||i6==i4||i6==i5)continue;
       fTrack=anEvent->GetTrack(i6);
       if(!(fTrack->InRPSelection()))continue;//RP condition  
       phi6=fTrack->Phi();
       for(Int_t i7=0;i7<nPrim;i7++)
       {
        if(i7==i1||i7==i2||i7==i3||i7==i4||i7==i5||i7==i6)continue;
        fTrack=anEvent->GetTrack(i7);
        if(!(fTrack->InRPSelection()))continue;//RP condition  
        phi7=fTrack->Phi();   
        //fill the fDirectCorrelations:
        //fDirectCorrelations->Fill(65.,cos(2.*n*phi1+n*phi2+n*phi3-n*phi4-n*phi5-n*phi6-n*phi7),1);//<7'>_{2n,n,n|n,n,n,n}
       }//end of for(Int_t i7=0;i7<nPrim;i7++)  
      }//end of for(Int_t i6=0;i6<nPrim;i6++)   
     }//end of for(Int_t i5=0;i5<nPrim;i5++)  
    }//end of for(Int_t i4=0;i4<nPrim;i4++)
   }//end of for(Int_t i3=0;i3<nPrim;i3++)
  }//end of for(Int_t i2=0;i2<nPrim;i2++) 
 }//end of for(Int_t i1=0;i1<nPrim;i1++)

 
  
 */
 /*  
    
     
       
 //<8'>_{n,n,n,n|n,n,n,n}
 for(Int_t i1=0;i1<nPrim;i1++)
 {
  fTrack=anEvent->GetTrack(i1);
  if(!((fTrack->Pt()>=0.5&&fTrack->Pt()<0.6)&&(fTrack->InPOISelection())))continue;//POI condition
  phi1=fTrack->Phi();
  for(Int_t i2=0;i2<nPrim;i2++)
  {
   if(i2==i1)continue;
   fTrack=anEvent->GetTrack(i2);
   if(!(fTrack->InRPSelection()))continue;//RP condition   
   phi2=fTrack->Phi();
   for(Int_t i3=0;i3<nPrim;i3++)
   { 
    if(i3==i1||i3==i2)continue;
    fTrack=anEvent->GetTrack(i3);
    if(!(fTrack->InRPSelection()))continue;//RP condition   
    phi3=fTrack->Phi();
    for(Int_t i4=0;i4<nPrim;i4++)
    {
     if(i4==i1||i4==i2||i4==i3)continue;
     fTrack=anEvent->GetTrack(i4);
     if(!(fTrack->InRPSelection()))continue;//RP condition  
     phi4=fTrack->Phi();
     for(Int_t i5=0;i5<nPrim;i5++)
     {
      if(i5==i1||i5==i2||i5==i3||i5==i4)continue;
      fTrack=anEvent->GetTrack(i5);
      if(!(fTrack->InRPSelection()))continue;//RP condition  
      phi5=fTrack->Phi();    
      for(Int_t i6=0;i6<nPrim;i6++)
      {
       if(i6==i1||i6==i2||i6==i3||i6==i4||i6==i5)continue;
       fTrack=anEvent->GetTrack(i6);
       if(!(fTrack->InRPSelection()))continue;//RP condition  
       phi6=fTrack->Phi();
       for(Int_t i7=0;i7<nPrim;i7++)
       {
        if(i7==i1||i7==i2||i7==i3||i7==i4||i7==i5||i7==i6)continue;
        fTrack=anEvent->GetTrack(i7);
        if(!(fTrack->InRPSelection()))continue;//RP condition  
        phi7=fTrack->Phi();
        for(Int_t i8=0;i8<nPrim;i8++)
        {
         if(i8==i1||i8==i2||i8==i3||i8==i4||i8==i5||i8==i6||i8==i7)continue;
         fTrack=anEvent->GetTrack(i8);
         if(!(fTrack->InRPSelection()))continue;//RP condition  
         phi8=fTrack->Phi();           
         //fill the fDirectCorrelations:
         //fDirectCorrelations->Fill(70.,cos(n*(phi1+phi2+phi3+phi4-phi5-phi6-phi7-phi8)),1);//<8'>_{n,n,n,n|n,n,n,n}
        }//end of for(Int_t i8=0;i8<nPrim;i8++) 
       }//end of for(Int_t i7=0;i7<nPrim;i7++)  
      }//end of for(Int_t i6=0;i6<nPrim;i6++)   
     }//end of for(Int_t i5=0;i5<nPrim;i5++)  
    }//end of for(Int_t i4=0;i4<nPrim;i4++)
   }//end of for(Int_t i3=0;i3<nPrim;i3++)
  }//end of for(Int_t i2=0;i2<nPrim;i2++) 
 }//end of for(Int_t i1=0;i1<nPrim;i1++)
 
 
 
 */ 
 
 
 
 
} // end of AliFlowAnalysisWithQCumulants::EvaluateNestedLoopsForDifferentialFlow(AliFlowEventSimple* anEvent)


//================================================================================================================================


void AliFlowAnalysisWithQCumulants::Finish()
{
 // calculate the final results
 
 fUseWeights = fUseWeightsBits->TestBitNumber(1); // to be improved 
 
 // compare correlations needed for integrated flow:
 if(fDirectCorrelations->GetBinContent(1) != 0 || fDirectCorrelationsW->GetBinContent(1) != 0) 
 {
  this->CompareDirectAndQCorrelationsForIntegratedFlow(fUseWeights);
 } 
 // compare correlations needed for integrated flow:
 if(fDirectCorrelationsDiffFlow->GetBinContent(1) != 0 || fDirectCorrelationsDiffFlowW->GetBinContent(1) != 0) 
 {
  this->CompareDirectAndQCorrelationsForDifferentialFlow(fUseWeights);
 } 
 
  
 
     
 //                      *************************************
 //                      **** CALCULATE THE FINAL RESULTS ****
 //                      *************************************    
  
 // integrated flow ('no-name') without weights:
 // calculate final results for no-name integrated flow without weights:
 this->CalculateFinalResultsForNoNameIntegratedFlow(kFALSE);
 
 // integrated flow ('no-name') with weights:
 // calculate final results for no-name integrated flow with weights:
 if(fUseWeights) this->CalculateFinalResultsForNoNameIntegratedFlow(fUseWeights);
 
 
 //            **** POI ****
 
 // differential flow (POI) without weights:
 // calculate final results for 2nd order differential flow of POIs without weights:
 this->CalculateFinalResultsForDifferentialFlow(fvn2ndPtEtaPOI,fvn2ndPtPOI,fvn2ndEtaPOI,f2pPtEtaPOI);
 // calculate final results for 4th order differential flow of POIs without weights:
 this->CalculateFinalResultsForDifferentialFlow(fvn4thPtEtaPOI,fvn4thPtPOI,fvn4thEtaPOI,f2pPtEtaPOI,f4pPtEtaPOI);
 // calculate final results for 6th order differential flow of POIs without weights:
 // this->CalculateFinalResultsForDifferentialFlow(fvn6thPtEtaPOI,fvn6thPtPOI,fvn6thEtaPOI,f2pPtEtaPOI,f4pPtEtaPOI,f6pPtEtaPOI);
 // calculate final results for 8th order differential flow of POIs without weights:
 // this->CalculateFinalResultsForDifferentialFlow(fvn8thPtEtaPOI,fvn8thPtPOI,fvn8thEtaPOI,f2pPtEtaPOI,f4pPtEtaPOI,f6pPtEtaPOI,f8pPtEtaPOI);
 
 // differential flow (POI) with weights:
 // calculate final results for 2nd order differential flow of POIs with weights:
 if(fUseWeights) this->CalculateFinalResultsForDifferentialFlow(fvn2ndPtEtaPOIW,fvn2ndPtPOIW,fvn2ndEtaPOIW,f2pPtEtaPOIW);
 // calculate final results for 4th order differential flow of POIs without weights:
 if(fUseWeights) this->CalculateFinalResultsForDifferentialFlow(fvn4thPtEtaPOIW,fvn4thPtPOIW,fvn4thEtaPOIW,f2pPtEtaPOIW,f4pPtEtaPOIW);
 // calculate final results for 6th order differential flow of POIs with weights:
 // if(fUseWeights) this->CalculateFinalResultsForDifferentialFlow(fvn6thPtEtaPOIW,fvn6thPtPOIW,fvn6thEtaPOIW,f2pPtEtaPOIW,f4pPtEtaPOIW,f6pPtEtaPOIW);
 // calculate final results for 8th order differential flow of POIs without weights:
 // if(fUseWeights) this->CalculateFinalResultsForDifferentialFlow(fvn8thPtEtaPOIW,fvn8thPtPOIW,fvn8thEtaPOIW,f2pPtEtaPOIW,f4pPtEtaPOIW,f6pPtEtaPOIW,f8pPtEtaPOIW);
  
 // integrated flow (POI) without weights:
 // calculate final results for integrated flow of POIs without weights:
 this->CalculateFinalResultsForRPandPOIIntegratedFlow(kFALSE,"POI");

 // integrated flow (POI) with weights:
 // calculate final results for integrated flow of POIs with weights:
 this->CalculateFinalResultsForRPandPOIIntegratedFlow(kTRUE,"POI");
 
 
 //            **** RP ****
 
 // differential flow (RP) without weights:
 // calculate final results for 2nd order differential flow of RPs without weights:
 this->CalculateFinalResultsForDifferentialFlow(fvn2ndPtEtaRP,fvn2ndPtRP,fvn2ndEtaRP,f2pPtEtaRP);
 // calculate final results for 4th order differential flow of RPs without weights:
 this->CalculateFinalResultsForDifferentialFlow(fvn4thPtEtaRP,fvn4thPtRP,fvn4thEtaRP,f2pPtEtaRP,f4pPtEtaRP);
 // calculate final results for 6th order differential flow of RPs without weights:
 // this->CalculateFinalResultsForDifferentialFlow(fvn6thPtEtaRP,fvn6thPtRP,fvn6thEtaRP,f2pPtEtaRP,f4pPtEtaRP,f6pPtEtaRP);
 // calculate final results for 8th order differential flow of RPs without weights:
 // this->CalculateFinalResultsForDifferentialFlow(fvn8thPtEtaRP,fvn8thPtRP,fvn8thEtaRP,f2pPtEtaRP,f4pPtEtaRP,f6pPtEtaRP,f8pPtEtaRP);
 
 // differential flow (RP) with weights:
 // calculate final results for 2nd order differential flow of RPs with weights:
 if(fUseWeights) this->CalculateFinalResultsForDifferentialFlow(fvn2ndPtEtaRPW,fvn2ndPtRPW,fvn2ndEtaRPW,f2pPtEtaRPW);
 // calculate final results for 4th order differential flow of RPs without weights:
 if(fUseWeights) this->CalculateFinalResultsForDifferentialFlow(fvn4thPtEtaRPW,fvn4thPtRPW,fvn4thEtaRPW,f2pPtEtaRPW,f4pPtEtaRPW);
 // calculate final results for 6th order differential flow of RPs with weights:
 // if(fUseWeights) this->CalculateFinalResultsForDifferentialFlow(fvn6thPtEtaRPW,fvn6thPtRPW,fvn6thEtaRPW,f2pPtEtaRPW,f4pPtEtaRPW,f6pPtEtaRPW);
 // calculate final results for 8th order differential flow of RPs without weights:
 // if(fUseWeights) this->CalculateFinalResultsForDifferentialFlow(fvn8thPtEtaRPW,fvn8thPtRPW,fvn8thEtaRPW,f2pPtEtaRPW,f4pPtEtaRPW,f6pPtEtaRPW,f8pPtEtaRPW);
 
 // integrated flow (RP) without weights: 
 // calculate final results for integrated flow of RPs without weights:
 this->CalculateFinalResultsForRPandPOIIntegratedFlow(kFALSE,"RP");

 // integrated flow (RP) with weights: 
 // calculate final results for integrated flow of POIs with weights:
 this->CalculateFinalResultsForRPandPOIIntegratedFlow(kTRUE,"RP");
 
   
 
 //             *****************************************************
 //             **** PRINT THE FINAL RESULTS FOR INTEGRATED FLOW ****
 //             *****************************************************       
  
 // print the final results for 'no-name' integrated flow without weights:
 this->PrintFinalResultsForIntegratedFlow(kFALSE,"NONAME"); // OK tested (just still nEvts and AvM)

 // print the final results for 'no-name' integrated flow with weights: 
 if(fUseWeights) this->PrintFinalResultsForIntegratedFlow(fUseWeights,"NONAME"); // OK tested (just still nEvts and AvM)

 // print the final results for RPs integrated flow without weights:
 this->PrintFinalResultsForIntegratedFlow(kFALSE,"RP");
 
 // print the final results for RPs integrated flow with weights:
 if(fUseWeights) this->PrintFinalResultsForIntegratedFlow(kTRUE,"RP");
 
 // print the final results for POIs integrated flow without weights:
 this->PrintFinalResultsForIntegratedFlow(kFALSE,"POI");
 
 // print the final results for POIs integrated flow with weights:
 if(fUseWeights) this->PrintFinalResultsForIntegratedFlow(kTRUE,"POI"); 
 
 //this->TempDeleteMe();
                                                                
} // end of AliFlowAnalysisWithQCumulants::Finish()


//================================================================================================================================


TProfile* AliFlowAnalysisWithQCumulants::MakePtProjection(TProfile2D *profilePtEta) const
{
 // project 2D profile onto pt axis to get 1D profile
 
 Int_t nBinsPt   = profilePtEta->GetNbinsX();
 Double_t dPtMin = (profilePtEta->GetXaxis())->GetXmin();
 Double_t dPtMax = (profilePtEta->GetXaxis())->GetXmax();
 
 Int_t nBinsEta   = profilePtEta->GetNbinsY();
 
 TProfile *profilePt = new TProfile("","",nBinsPt,dPtMin,dPtMax); 
 
 for(Int_t p=1;p<=nBinsPt;p++)
 {
  Double_t contentPt = 0.;
  Double_t entryPt = 0.;
  for(Int_t e=1;e<=nBinsEta;e++)
  {
   contentPt += (profilePtEta->GetBinContent(profilePtEta->GetBin(p,e)))
              * (profilePtEta->GetBinEntries(profilePtEta->GetBin(p,e)));
   entryPt   += (profilePtEta->GetBinEntries(profilePtEta->GetBin(p,e)));
  }
  profilePt->SetBinContent(p,contentPt);
  profilePt->SetBinEntries(p,entryPt);
 }
 
 return profilePt;
 
} // end of TProfile* AliFlowAnalysisWithQCumulants::MakePtProjection(TProfile2D *profilePtEta)


//================================================================================================================================


TProfile* AliFlowAnalysisWithQCumulants::MakeEtaProjection(TProfile2D *profilePtEta) const
{
 // project 2D profile onto eta axis to get 1D profile
 
 Int_t nBinsEta   = profilePtEta->GetNbinsY();
 Double_t dEtaMin = (profilePtEta->GetYaxis())->GetXmin();
 Double_t dEtaMax = (profilePtEta->GetYaxis())->GetXmax();
 
 Int_t nBinsPt = profilePtEta->GetNbinsX();
 
 TProfile *profileEta = new TProfile("","",nBinsEta,dEtaMin,dEtaMax); 
 
 for(Int_t e=1;e<=nBinsEta;e++)
 {
  Double_t contentEta = 0.;
  Double_t entryEta = 0.;
  for(Int_t p=1;p<=nBinsPt;p++)
  {
   contentEta += (profilePtEta->GetBinContent(profilePtEta->GetBin(p,e)))
              * (profilePtEta->GetBinEntries(profilePtEta->GetBin(p,e)));
   entryEta   += (profilePtEta->GetBinEntries(profilePtEta->GetBin(p,e)));
  }
  profileEta->SetBinContent(e,contentEta);
  profileEta->SetBinEntries(e,entryEta);
 }
 
 return profileEta;
 
} // end of TProfile* AliFlowAnalysisWithQCumulants::MakeEtaProjection(TProfile2D *profilePtEta)


//================================================================================================================================


void AliFlowAnalysisWithQCumulants::CalculateFinalResultsForNoNameIntegratedFlow(Bool_t useWeights)
{
 // calculate final results for 'no-name' integrated flow
 
 // 2-, 4-, 6- and 8-particle azimuthal correlation:
 Double_t two   = 0.; // <<2>>_{n|n}
 Double_t four  = 0.; // <<4>>_{n,n|n,n}
 Double_t six   = 0.; // <<6>>_{n,n,n|n,n,n}
 Double_t eight = 0.; // <<8>>_{n,n,n,n|n,n,n,n}
 
 if(!(useWeights))
 {
  two   = fQCorrelations->GetBinContent(1);  
  four  = fQCorrelations->GetBinContent(11); 
  six   = fQCorrelations->GetBinContent(24); 
  eight = fQCorrelations->GetBinContent(31); 
 }

 if(useWeights)
 {
  two   = fQCorrelationsW->GetBinContent(1);  
  four  = fQCorrelationsW->GetBinContent(41); 
  six   = fQCorrelationsW->GetBinContent(81); 
  eight = fQCorrelationsW->GetBinContent(121);  
 }
 
 // 2nd, 4th, 6th and 8th order Q-cumulant:
 Double_t secondOrderQCumulant = two; // c_n{2} 
 Double_t fourthOrderQCumulant = four-2.*pow(two,2.); // c_n{4}
 Double_t sixthOrderQCumulant  = six-9.*two*four+12.*pow(two,3.); // c_n{6}
 Double_t eightOrderQCumulant  = eight-16.*two*six-18.*pow(four,2.)+144.*pow(two,2.)*four-144.*pow(two,4.); // c_n{8} 
 
 if(useWeights) sixthOrderQCumulant = 0.; // to be removed (once 6th order with weights is calculated)
 if(useWeights) eightOrderQCumulant = 0.; // to be removed (once 8th order with weights is calculated)
 
 // "no-name" integrated flow estimates from Q-cumulants:
 Double_t dVn2 = 0.,dVn4 = 0.,dVn6 = 0.,dVn8 = 0.;
 // Double_t sd2=0.,sd4=0.,sd6=0.,sd8=0.; // to be improved (errors needed)
 if(secondOrderQCumulant>0.)
 {
  // v_n{2}
  dVn2 = pow(secondOrderQCumulant,0.5); 
  if(!(useWeights)) 
  { 
   fIntFlowResultsQC->SetBinContent(1,dVn2);
   fIntFlowResultsQC->SetBinError(1,0.); // to be improved
  }
  if(useWeights)
  { 
   fIntFlowResultsQCW->SetBinContent(1,dVn2);
   fIntFlowResultsQCW->SetBinError(1,0.); // to be improved 
  }
  
  // fill common histogram:
  fCommonHistsResults2nd->FillIntegratedFlow(dVn2, 0.); // to be improved 
  
 } 
 if(fourthOrderQCumulant<0.)
 {
  // v_n{4}
  dVn4 = pow(-fourthOrderQCumulant,1./4.); 
  if(!(useWeights)) 
  {
   fIntFlowResultsQC->SetBinContent(2,dVn4);
   fIntFlowResultsQC->SetBinError(2,0.); // to be improved 
  } 
  if(useWeights) 
  {
   fIntFlowResultsQCW->SetBinContent(2,dVn4);
   fIntFlowResultsQCW->SetBinError(2,0.); // to be improved 
  }
   
  // fill common histogram:
  fCommonHistsResults4th->FillIntegratedFlow(dVn4, 0.); // to be improved 
  
 } 
 if(sixthOrderQCumulant>0.)
 {
  // v_n{6}
  dVn6 = pow((1./4.)*sixthOrderQCumulant,1./6.); 
  if(!(useWeights)) 
  {
   fIntFlowResultsQC->SetBinContent(3,dVn6);
   fIntFlowResultsQC->SetBinError(3,0.); // to be improved
  } 
  if(useWeights)
  {
   fIntFlowResultsQCW->SetBinContent(3,dVn6);
   fIntFlowResultsQCW->SetBinError(3,0.); // to be improved
  }
  
  // fill common histogram:
  fCommonHistsResults6th->FillIntegratedFlow(dVn6, 0.); // to be improved 
  
 } 
 if(eightOrderQCumulant<0.)
 {
  // v_n{8}
  dVn8 = pow((-1./33.)*eightOrderQCumulant,1./8.); 
  if(!(useWeights))
  {
   fIntFlowResultsQC->SetBinContent(4,dVn8);
   fIntFlowResultsQC->SetBinError(4,0.); // to be improved
  } 
  if(useWeights)
  {
   fIntFlowResultsQCW->SetBinContent(4,dVn8);
   fIntFlowResultsQCW->SetBinError(4,0.); // to be improved
  } 
  
  // fill common histogram:
  fCommonHistsResults8th->FillIntegratedFlow(dVn8, 0.); // to be improved 

 }
 
} // end of AliFlowAnalysisWithQCumulants::CalculateFinalResultsForNoNameIntegratedFlow(Bool_t useWeights)


//================================================================================================================================


void AliFlowAnalysisWithQCumulants::CalculateFinalResultsForRPandPOIIntegratedFlow(Bool_t useWeights, TString type)
{
 // calculate final results for integrated flow of RPs and POIs 
     
 TH1F *yield2ndPt = NULL;
 TH1F *yield4thPt = NULL;
 TH1F *yield6thPt = NULL;
 TH1F *yield8thPt = NULL;
 
 if(type == "POI")
 {
  yield2ndPt = new TH1F(*(fCommonHists2nd->GetHistPtPOI()));
  yield4thPt = new TH1F(*(fCommonHists4th->GetHistPtPOI()));
  yield6thPt = new TH1F(*(fCommonHists6th->GetHistPtPOI()));
  yield8thPt = new TH1F(*(fCommonHists8th->GetHistPtPOI()));
 } 
 else if (type == "RP")
 {
  yield2ndPt = new TH1F(*(fCommonHists2nd->GetHistPtRP()));
  yield4thPt = new TH1F(*(fCommonHists4th->GetHistPtRP()));
  yield6thPt = new TH1F(*(fCommonHists6th->GetHistPtRP()));
  yield8thPt = new TH1F(*(fCommonHists8th->GetHistPtRP()));
 } 
 
 Int_t nBinsPt = yield2ndPt->GetNbinsX();
 
 TH1D *flow2ndPt = NULL;
 TH1D *flow4thPt = NULL;
 TH1D *flow6thPt = NULL;
 TH1D *flow8thPt = NULL;
 
 if(!(useWeights))
 { 
  if(type == "POI")
  {
   flow2ndPt = new TH1D(*fvn2ndPtPOI);
   flow4thPt = new TH1D(*fvn4thPtPOI);
   flow6thPt = new TH1D(*fvn6thPtPOI);
   flow8thPt = new TH1D(*fvn8thPtPOI);
  }
  else if (type == "RP")
  {
   flow2ndPt = new TH1D(*fvn2ndPtRP);
   flow4thPt = new TH1D(*fvn4thPtRP);
   flow6thPt = new TH1D(*fvn6thPtRP);
   flow8thPt = new TH1D(*fvn8thPtRP);
  }
 }  
 else if (useWeights)
 {
  if(type == "POI")
  {
   flow2ndPt = new TH1D(*fvn2ndPtPOIW);
   flow4thPt = new TH1D(*fvn4thPtPOIW);
   flow6thPt = new TH1D(*fvn6thPtPOIW);
   flow8thPt = new TH1D(*fvn8thPtPOIW);
  }
  else if (type == "RP")
  {
   flow2ndPt = new TH1D(*fvn2ndPtRPW);
   flow4thPt = new TH1D(*fvn4thPtRPW);
   flow6thPt = new TH1D(*fvn6thPtRPW);
   flow8thPt = new TH1D(*fvn8thPtRPW);
  } 
 } 
 
 Double_t dvn2nd = 0., dvn4th = 0., dvn6th = 0., dvn8th = 0.; // differential flow
 Double_t dVn2nd = 0., dVn4th = 0., dVn6th = 0., dVn8th = 0.; // integrated flow 
 Double_t dSd2nd = 0., dSd4th = 0., dSd6th = 0., dSd8th = 0.; // error on integrated flow (to be improved - calculation needed)

 Double_t dYield2nd = 0., dYield4th = 0., dYield6th = 0., dYield8th = 0.; // pt yield 
 Double_t dSum2nd = 0., dSum4th = 0., dSum6th = 0., dSum8th = 0.; // needed for normalizing integrated flow
 
 // looping over pt bins:
 for(Int_t p=1;p<nBinsPt+1;p++)
 {
 
  dvn2nd = flow2ndPt->GetBinContent(p);
  dvn4th = flow4thPt->GetBinContent(p);
  dvn6th = flow6thPt->GetBinContent(p);
  dvn8th = flow8thPt->GetBinContent(p);

  dYield2nd = yield2ndPt->GetBinContent(p);
  dYield4th = yield4thPt->GetBinContent(p);
  dYield6th = yield6thPt->GetBinContent(p);
  dYield8th = yield8thPt->GetBinContent(p);
  
  dVn2nd += dvn2nd*dYield2nd;
  dVn4th += dvn4th*dYield4th;
  dVn6th += dvn6th*dYield6th;
  dVn8th += dvn8th*dYield8th;
  
  dSum2nd += dYield2nd;
  dSum4th += dYield4th;
  dSum6th += dYield6th;
  dSum8th += dYield8th;
  
  // ... to be improved - errors needed to be calculated   
  
 } // end of for(Int_t p=1;p<nBinsPt+1;p++)

 // normalizing the results for integrated flow:
 if(dSum2nd) dVn2nd/=dSum2nd;
 if(dSum4th) dVn4th/=dSum4th;
 if(dSum6th) dVn6th/=dSum6th;
 if(dSum8th) dVn8th/=dSum8th;
 
 // storing the results for integrated flow:
 if(!(useWeights))
 {
  if(type == "POI")
  {
   // 2nd:
   fIntFlowResultsPOIQC->SetBinContent(1,dVn2nd);
   fIntFlowResultsPOIQC->SetBinError(1,dSd2nd);
   // 4th:
   fIntFlowResultsPOIQC->SetBinContent(2,dVn4th);
   fIntFlowResultsPOIQC->SetBinError(2,dSd4th);
   // 6th:
   fIntFlowResultsPOIQC->SetBinContent(3,dVn6th);
   fIntFlowResultsPOIQC->SetBinError(3,dSd6th);
   // 8th:
   fIntFlowResultsPOIQC->SetBinContent(4,dVn8th);
   fIntFlowResultsPOIQC->SetBinError(4,dSd8th);
  }
  else if (type == "RP")
  {
   // 2nd:
   fIntFlowResultsRPQC->SetBinContent(1,dVn2nd);
   fIntFlowResultsRPQC->SetBinError(1,dSd2nd);
   // 4th:
   fIntFlowResultsRPQC->SetBinContent(2,dVn4th);
   fIntFlowResultsRPQC->SetBinError(2,dSd4th);
   // 6th:
   fIntFlowResultsRPQC->SetBinContent(3,dVn6th);
   fIntFlowResultsRPQC->SetBinError(3,dSd6th);
   // 8th:
   fIntFlowResultsRPQC->SetBinContent(4,dVn8th);
   fIntFlowResultsRPQC->SetBinError(4,dSd8th);
  }
 } 
 else if (useWeights)
 {
  if(type == "POI")
  {
   // 2nd:
   fIntFlowResultsPOIQCW->SetBinContent(1,dVn2nd);
   fIntFlowResultsPOIQCW->SetBinError(1,dSd2nd);
   // 4th:
   fIntFlowResultsPOIQCW->SetBinContent(2,dVn4th);
   fIntFlowResultsPOIQCW->SetBinError(2,dSd4th);
   // 6th:
   fIntFlowResultsPOIQCW->SetBinContent(3,dVn6th);
   fIntFlowResultsPOIQCW->SetBinError(3,dSd6th);
   // 8th:
   fIntFlowResultsPOIQCW->SetBinContent(4,dVn8th);
   fIntFlowResultsPOIQCW->SetBinError(4,dSd8th);
  }
  else if (type == "RP")
  {
   // 2nd:
   fIntFlowResultsRPQCW->SetBinContent(1,dVn2nd);
   fIntFlowResultsRPQCW->SetBinError(1,dSd2nd);
   // 4th:
   fIntFlowResultsRPQCW->SetBinContent(2,dVn4th);
   fIntFlowResultsRPQCW->SetBinError(2,dSd4th);
   // 6th:
   fIntFlowResultsRPQCW->SetBinContent(3,dVn6th);
   fIntFlowResultsRPQCW->SetBinError(3,dSd6th);
   // 8th:
   fIntFlowResultsRPQCW->SetBinContent(4,dVn8th);
   fIntFlowResultsRPQCW->SetBinError(4,dSd8th);
  }
 }
 
 // storing the results for integrated flow in common histos:
 // to be improved - now they are being filled twice ...
 if(type == "POI")
 {
  fCommonHistsResults2nd->FillIntegratedFlowPOI(dVn2nd,0.); // to be improved (errors)
  fCommonHistsResults4th->FillIntegratedFlowPOI(dVn4th,0.); // to be improved (errors)
  fCommonHistsResults6th->FillIntegratedFlowPOI(dVn6th,0.); // to be improved (errors)
  fCommonHistsResults8th->FillIntegratedFlowPOI(dVn8th,0.); // to be improved (errors)
 }
 else if (type == "RP")
 {
  fCommonHistsResults2nd->FillIntegratedFlowRP(dVn2nd,0.); // to be improved (errors)
  fCommonHistsResults4th->FillIntegratedFlowRP(dVn4th,0.); // to be improved (errors)
  fCommonHistsResults6th->FillIntegratedFlowRP(dVn6th,0.); // to be improved (errors)
  fCommonHistsResults8th->FillIntegratedFlowRP(dVn8th,0.); // to be improved (errors)
 }
 
 delete flow2ndPt;
 delete flow4thPt;
 delete flow6thPt;
 delete flow8thPt;
 
 delete yield2ndPt;
 delete yield4thPt;
 delete yield6thPt;
 delete yield8thPt;
  
} // end of AliFlowAnalysisWithQCumulants::CalculateFinalResultsForRPandPOIIntegratedFlow(Bool_t useWeights, TString type)


//==================================================================================================================================


void AliFlowAnalysisWithQCumulants::CalculateFinalResultsForDifferentialFlow(
                                                        TH2D *flowPtEta, TH1D *flowPt, TH1D *flowEta, 
                                                        TProfile2D *profile2ndPtEta, TProfile2D *profile4thPtEta, 
                                                        TProfile2D *profile6thPtEta, TProfile2D *profile8thPtEta)
{
 // calculate and store the final results for integrated flow
 
 TString *namePtEta = new TString();
 TString *type = new TString();
 TString *order2nd = new TString();
 TString *order4th = new TString();
 TString *order6th = new TString();
 TString *order8th = new TString(); 
 TString *w = new TString();

 if(profile2ndPtEta) *namePtEta = profile2ndPtEta->GetName();
 if(namePtEta->Contains("POI")) *type = "POI";
 if(namePtEta->Contains("RP")) *type  = "RP";
 if(namePtEta->Contains("W")) *w      = "W";
 if(namePtEta->Contains("2")) *order2nd  = "2";
 
 if(profile4thPtEta) *namePtEta = profile4thPtEta->GetName();
 if(namePtEta->Contains("4")) *order4th = "4";

 if(profile6thPtEta) *namePtEta = profile6thPtEta->GetName();
 if(namePtEta->Contains("6")) *order6th = "6";
 
 if(profile8thPtEta) *namePtEta = profile8thPtEta->GetName();
 if(namePtEta->Contains("8")) *order8th = "8";
 
 TProfile *profile2ndPt = NULL;
 TProfile *profile4thPt = NULL;
 TProfile *profile6thPt = NULL;
 TProfile *profile8thPt = NULL;

 TProfile *profile2ndEta = NULL;
 TProfile *profile4thEta = NULL;
 TProfile *profile6thEta = NULL;
 TProfile *profile8thEta = NULL;
  
 if(*order2nd == "2")
 {
  profile2ndPt  = new TProfile(*(this->MakePtProjection(profile2ndPtEta))); 
  profile2ndEta = new TProfile(*(this->MakeEtaProjection(profile2ndPtEta))); 
  if(*order4th == "4")
  {
   profile4thPt  = new TProfile(*(this->MakePtProjection(profile4thPtEta))); 
   profile4thEta = new TProfile(*(this->MakeEtaProjection(profile4thPtEta))); 
   if(*order6th == "6")
   {
    profile6thPt  = new TProfile(*(this->MakePtProjection(profile6thPtEta))); 
    profile6thEta = new TProfile(*(this->MakeEtaProjection(profile6thPtEta))); 
    if(*order8th == "8")
    {
     profile8thPt  = new TProfile(*(this->MakePtProjection(profile8thPtEta))); 
     profile8thEta = new TProfile(*(this->MakeEtaProjection(profile8thPtEta))); 
    }     
   }    
  } 
 }
 
 Int_t nBinsPt  = profile2ndPt->GetNbinsX();
 Int_t nBinsEta = profile2ndEta->GetNbinsX();
 
 Double_t dV2 = 0.;
 Double_t dV4 = 0.;
 Double_t dV6 = 0.;
 Double_t dV8 = 0.; 
 
 if(!(*w == "W"))
 {
  dV2 = fIntFlowResultsQC->GetBinContent(1);  
  dV4 = fIntFlowResultsQC->GetBinContent(2); 
  dV6 = fIntFlowResultsQC->GetBinContent(3); 
  dV8 = fIntFlowResultsQC->GetBinContent(4); 
 } 
 else if(*w == "W")
 {
  dV2 = fIntFlowResultsQCW->GetBinContent(1);  
  dV4 = fIntFlowResultsQCW->GetBinContent(2); 
  dV6 = fIntFlowResultsQCW->GetBinContent(3); 
  dV8 = fIntFlowResultsQCW->GetBinContent(4); 
 }    
 
 // 3D (pt,eta): 
 Double_t twoPrimePtEta   = 0.; // <<2'>> (pt,eta) 
 Double_t fourPrimePtEta  = 0.; // <<4'>> (pt,eta)  
 //Double_t sixPrimePtEta   = 0.; // <<6'>> (pt,eta) 
 //Double_t eightPrimePtEta = 0.; // <<8'>> (pt,eta) 
 Double_t secondOrderDiffFlowCumulantPtEta = 0.; // d_n{2,Q} (pt,eta)
 Double_t fourthOrderDiffFlowCumulantPtEta = 0.; // d_n{4,Q} (pt,eta) 
 //Double_t sixthOrderDiffFlowCumulantPtEta = 0.; // d_n{6,Q} (pt,eta)
 //Double_t eightOrderDiffFlowCumulantPtEta = 0.; // d_n{8,Q} (pt,eta)2nd
 Double_t dv2PtEta = 0.; // v'_n{2} (pt,eta) 
 Double_t dv4PtEta = 0.; // v'_n{4} (pt,eta) 
 //Double_t dv6PtEta = 0.; // v'_n{6} (pt,eta) 
 //Double_t dv8PtEta = 0.; // v'_n{8} (pt,eta)  

 // 2D (pt):   
 Double_t twoPrimePt   = 0.; // <<2'>> (pt)  
 Double_t fourPrimePt  = 0.; // <<4'>> (pt) 
 //Double_t sixPrimePt   = 0.; // <<6'>> (pt) 
 //Double_t eightPrimePt = 0.; // <<8'>> (pt)          
 Double_t secondOrderDiffFlowCumulantPt = 0.; // d_n{2,Q} (pt) 
 Double_t fourthOrderDiffFlowCumulantPt = 0.; // d_n{4,Q} (pt)  
 //Double_t sixthOrderDiffFlowCumulantPt = 0.; // d_n{6,Q} (pt)
 //Double_t eightOrderDiffFlowCumulantPt = 0.; // d_n{8,Q} (pt)
 Double_t dv2Pt = 0.; // v'_n{2} (pt)
 Double_t dv4Pt = 0.; // v'_n{4} (pt)
 //Double_t dv6Pt = 0.; // v'_n{6} (pt) 
 //Double_t dv8Pt = 0.; // v'_n{8} (pt)  

 // 2D (eta):           
 Double_t twoPrimeEta   = 0.; // <<2'>> (eta)  
 Double_t fourPrimeEta  = 0.; // <<4>> (eta) 
 //Double_t sixPrimeEta   = 0.; // <<6>> (eta) 
 //Double_t eightPrimeEta = 0.; // <<8'>> (eta)  
 Double_t secondOrderDiffFlowCumulantEta = 0.; // d_n{2,Q} (eta)
 Double_t fourthOrderDiffFlowCumulantEta = 0.; // d_n{4,Q} (eta) 
 //Double_t sixthOrderDiffFlowCumulantEta = 0.; // d_n{6,Q} (eta) 
 //Double_t eightOrderDiffFlowCumulantEta = 0.; // d_n{8,Q} (eta) 
 Double_t dv2Eta = 0.; // v'_n{2} (eta)
 Double_t dv4Eta = 0.; // v'_n{4} (eta)
 //Double_t dv6Eta = 0.; // v'_n{6} (eta) 
 //Double_t dv8Eta = 0.; // v'_n{8} (eta)
 

 // looping over (pt,eta) bins to calculate v'(pt,eta) 
 for(Int_t p=1;p<nBinsPt+1;p++)
 {
  for(Int_t e=1;e<nBinsEta+1;e++)
  {
  
   // 2nd order: 
   twoPrimePtEta = profile2ndPtEta->GetBinContent(profile2ndPtEta->GetBin(p,e));
   secondOrderDiffFlowCumulantPtEta = twoPrimePtEta;
   if(dV2)
   {
    dv2PtEta = secondOrderDiffFlowCumulantPtEta/dV2;
    if(*order2nd == "2") 
    {
     flowPtEta->SetBinContent(p,e,dv2PtEta);   
    } 
   }
   
   // 4th order: 
   if(*order4th == "4" || *order6th == "6" || *order8th == "8")
   {
    fourPrimePtEta = profile4thPtEta->GetBinContent(profile4thPtEta->GetBin(p,e));
    fourthOrderDiffFlowCumulantPtEta = fourPrimePtEta - 2.*twoPrimePtEta*pow(dV2,2.); // to be improved (correlations instead of pow(dV2,2.))
    if(dV4)
    {
     dv4PtEta = -fourthOrderDiffFlowCumulantPtEta/pow(dV4,3);
     if(*order4th == "4")
     {
      flowPtEta->SetBinContent(p,e,dv4PtEta);
     } 
    }
   }    
   
  } // end of for(Int_t e=1;e<nBinsEta+1;e++)
 } // end of for(Int_t p=1;p<nBinsPt+1;p++) 
   
   
 // looping over (pt) bins to calcualate v'(pt)
 for(Int_t p=1;p<nBinsPt+1;p++)
 {
 
  // 2nd order: 
  twoPrimePt = profile2ndPt->GetBinContent(p);
  secondOrderDiffFlowCumulantPt = twoPrimePt;
  if(dV2)
  {
   dv2Pt = secondOrderDiffFlowCumulantPt/dV2;
   if(*order2nd == "2") 
   {
    flowPt->SetBinContent(p,dv2Pt);
   }
   
   // common control histos: (to be improved fill only once. now they are filled first without weights and then with weights):
   if(namePtEta->Contains("POI") && *order2nd == "2")
   {
    fCommonHistsResults2nd->FillDifferentialFlowPtPOI(p,dv2Pt,0.); //to be improved (errors && bb or bb+1 ?)
   } 
   else if(namePtEta->Contains("RP") && *order2nd == "2")
   {
    fCommonHistsResults2nd->FillDifferentialFlowPtRP(p,dv2Pt,0.); //to be improved (errors && bb or bb+1 ?)
   }
   
  }
  
  // 4th order: 
  if(*order4th == "4" || *order6th == "6" || *order8th == "8")
  {
   fourPrimePt = profile4thPt->GetBinContent(profile4thPt->GetBin(p));
   fourthOrderDiffFlowCumulantPt = fourPrimePt - 2.*twoPrimePt*pow(dV2,2.); // to be improved (correlations instead of pow(dV2,2.))
   if(dV4)
   {
    dv4Pt = -fourthOrderDiffFlowCumulantPt/pow(dV4,3);
    if(*order4th == "4") 
    {
     flowPt->SetBinContent(p,dv4Pt);
    }
    
    // common control histos: (to be improved):
    if(namePtEta->Contains("POI") && *order4th == "4")
    {
     fCommonHistsResults4th->FillDifferentialFlowPtPOI(p,dv4Pt,0.); //to be improved (errors && bb or bb+1 ?)
    } 
    else if(namePtEta->Contains("RP") && *order4th == "4" )
    {
     fCommonHistsResults4th->FillDifferentialFlowPtRP(p,dv4Pt,0.); //to be improved (errors && bb or bb+1 ?)
    }
        
   }
  }    
  
 } // end of for(Int_t p=1;p<nBinsPt+1;p++)  
 
 
 // looping over (eta) bins to calcualate v'(eta)
 for(Int_t e=1;e<nBinsEta+1;e++)
 {
 
  // 2nd order: 
  twoPrimeEta = profile2ndEta->GetBinContent(e);
  secondOrderDiffFlowCumulantEta = twoPrimeEta;
  if(dV2)
  {
   dv2Eta = secondOrderDiffFlowCumulantEta/dV2;
   if(*order2nd == "2") 
   {
    flowEta->SetBinContent(e,dv2Eta);
   }
   
   // common control histos: (to be improved):
   if(namePtEta->Contains("POI") && *order2nd == "2")
   {
    fCommonHistsResults2nd->FillDifferentialFlowEtaPOI(e,dv2Eta,0.); //to be improved (errors && bb or bb+1 ?)
   } 
   else if(namePtEta->Contains("RP") && *order2nd == "2")
   {
    fCommonHistsResults2nd->FillDifferentialFlowEtaRP(e,dv2Eta,0.); //to be improved (errors && bb or bb+1 ?)
   }
     
  }
  
  // 4th order: 
  if(*order4th == "4" || *order6th == "6" || *order8th == "8")
  {
   fourPrimeEta = profile4thEta->GetBinContent(profile4thEta->GetBin(e));
   fourthOrderDiffFlowCumulantEta = fourPrimeEta - 2.*twoPrimeEta*pow(dV2,2.); // to be improved (correlations instead of pow(dV2,2.))
   if(dV4)
   {
    dv4Eta = -fourthOrderDiffFlowCumulantEta/pow(dV4,3);
    if(*order4th == "4")
    {
     flowEta->SetBinContent(e,dv4Eta);
    }
    
    // common control histos: (to be improved):
    if(namePtEta->Contains("POI") && *order4th == "4")
    {
     fCommonHistsResults4th->FillDifferentialFlowEtaPOI(e,dv4Eta,0.); //to be improved (errors && bb or bb+1 ?)
    } 
    else if(namePtEta->Contains("RP") && *order4th == "4")
    {
     fCommonHistsResults4th->FillDifferentialFlowEtaRP(e,dv4Eta,0.); //to be improved (errors && bb or bb+1 ?)
    }
   
   }
  }    
  
 } // end of for(Int_t e=1;e<nBinsEta+1;e++)    
    
 delete namePtEta;
 delete type;
 delete order2nd;
 delete order4th;
 delete order6th;
 delete order8th;
 delete w;
 delete profile2ndPt;
 delete profile4thPt;
 delete profile6thPt;
 delete profile8thPt;
 delete profile2ndEta;
 delete profile4thEta;
 delete profile6thEta;
 delete profile8thEta;

} // end of AliFlowAnalysisWithQCumulants::CalculateFinalResultsForDifferentialFlow(Bool_t useWeights, TString type)


//================================================================================================================================


void AliFlowAnalysisWithQCumulants::PrintFinalResultsForIntegratedFlow(Bool_t useWeights, TString type)
{
 // printing on the screen the final results for integrated flow ('no-name', POI and RP, without/with weights)
 
 Int_t n = 2; // to be improved / removed
 
 Double_t nEvtsNoName = (fCommonHists2nd->GetHistMultRP())->GetEntries(); // to be improved 
 Double_t dMultNoName = (fCommonHists2nd->GetHistMultRP())->GetMean(); // to be improved 
 Double_t nEvtsPOI = (fCommonHists2nd->GetHistMultPOI())->GetEntries(); // to be improved 
 Double_t dMultPOI = (fCommonHists2nd->GetHistMultPOI())->GetMean(); // to be improved 
 Double_t nEvtsRP = (fCommonHists2nd->GetHistMultRP())->GetEntries(); // to be improved 
 Double_t dMultRP = (fCommonHists2nd->GetHistMultRP())->GetMean(); // to be improved 
 
 TH1D *finalResultsIntFlow = NULL;
 
 if(!(useWeights))
 {
  if(type == "NONAME") finalResultsIntFlow = new TH1D(*fIntFlowResultsQC);
  if(type == "POI") finalResultsIntFlow = new TH1D(*fIntFlowResultsPOIQC);
  if(type == "RP") finalResultsIntFlow = new TH1D(*fIntFlowResultsRPQC);
 }
 
 if(useWeights)
 {
  if(type == "NONAME") finalResultsIntFlow = new TH1D(*fIntFlowResultsQCW);
  if(type == "POI") finalResultsIntFlow = new TH1D(*fIntFlowResultsPOIQCW);
  if(type == "RP") finalResultsIntFlow = new TH1D(*fIntFlowResultsRPQCW);
 }
 
 Double_t dVn[4] = {0.}; // array to hold Vn{2}, Vn{4}, Vn{6} and Vn{8}   
 Double_t dVnErr[4] = {0.}; // array to hold errors of Vn{2}, Vn{4}, Vn{6} and Vn{8}   
  
 if(finalResultsIntFlow)
 {
  for(Int_t i=0;i<4;i++)
  {
   dVn[i] = finalResultsIntFlow->GetBinContent(i+1);  
   dVnErr[i] = finalResultsIntFlow->GetBinError(i+1);
  }
 }
  
 TString title = " flow estimates from Q-cumulants"; 
 TString subtitle = "    ("; 
 
 if(!(useWeights))
 {
  subtitle.Append(type);
  subtitle.Append(", without weights)");
 }
 
 if(useWeights)
 {
  subtitle.Append(type);
  subtitle.Append(", with weights)");
 }
  
 cout<<endl;
 cout<<"**********************************"<<endl;
 cout<<"**********************************"<<endl;
 cout<<title.Data()<<endl; 
 cout<<subtitle.Data()<<endl; 
 cout<<endl;
  
 for(Int_t i=0;i<4;i++)
 {
  if(dVn[i]>=0.)
  {
   cout<<"  v_"<<n<<"{"<<2*(i+1)<<"} = "<<dVn[i]<<" +/- "<<dVnErr[i]<<endl;
  }
  else
  {
   cout<<"  v_"<<n<<"{"<<2*(i+1)<<"} = Im"<<endl;
  }  
 }

 cout<<endl;
 if(type == "NONAME")
 {
  cout<<"     nEvts = "<<nEvtsNoName<<", AvM = "<<dMultNoName<<endl; // to be improved
 }
 else if (type == "RP")
 {
  cout<<"     nEvts = "<<nEvtsRP<<", AvM = "<<dMultRP<<endl; // to be improved  
 } 
 else if (type == "POI")
 {
  cout<<"     nEvts = "<<nEvtsPOI<<", AvM = "<<dMultPOI<<endl; // to be improved  
 } 
 cout<<"**********************************"<<endl;
 cout<<"**********************************"<<endl;
 cout<<endl; 
  
}// end of AliFlowAnalysisWithQCumulants::PrintFinalResultsForIntegratedFlow(Bool_t useWeights=kTRUE, TString type="NONAME");


//================================================================================================================================


void AliFlowAnalysisWithQCumulants::CompareDirectAndQCorrelationsForIntegratedFlow(Bool_t useWeights)
{
 // compare correlations needed for int. flow calculated with nested loops and those calculated from Q-vectors

 cout<<endl;
 cout<<"   *************************************"<<endl;
 cout<<"   **** cross-checking the formulas ****"<<endl;
 cout<<"   ****     for integrated flow     ****"<<endl;
 cout<<"   *************************************"<<endl;
 cout<<endl;
 
 if(!(useWeights))
 {
  cout<<"<2>_{1n,1n} from Q-vectors    = "<<fQCorrelations->GetBinContent(1)<<endl;
  cout<<"<2>_{1n,1n} from nested loops = "<<fDirectCorrelations->GetBinContent(1)<<endl;
  cout<<endl;
  cout<<"<2>_{2n,2n} from Q-vectors    = "<<fQCorrelations->GetBinContent(2)<<endl;
  cout<<"<2>_{2n,2n} from nested loops = "<<fDirectCorrelations->GetBinContent(2)<<endl;
  cout<<endl;
  cout<<"<2>_{3n,3n} from Q-vectors    = "<<fQCorrelations->GetBinContent(3)<<endl;
  cout<<"<2>_{3n,3n} from nested loops = "<<fDirectCorrelations->GetBinContent(3)<<endl;
  cout<<endl;
  cout<<"<2>_{4n,4n} from Q-vectors    = "<<fQCorrelations->GetBinContent(4)<<endl;
  cout<<"<2>_{4n,4n} from nested loops = "<<fDirectCorrelations->GetBinContent(4)<<endl;
  cout<<endl; 
  cout<<"<3>_{2n|1n,1n} from Q-vectors    = "<<fQCorrelations->GetBinContent(6)<<endl;
  cout<<"<3>_{2n|1n,1n} from nested loops = "<<fDirectCorrelations->GetBinContent(6)<<endl;
  cout<<endl;
  cout<<"<3>_{3n|2n,1n} from Q-vectors    = "<<fQCorrelations->GetBinContent(7)<<endl;
  cout<<"<3>_{3n|2n,1n} from nested loops = "<<fDirectCorrelations->GetBinContent(7)<<endl;
  cout<<endl;
  cout<<"<3>_{4n,2n,2n} from Q-vectors    = "<<fQCorrelations->GetBinContent(8)<<endl;
  cout<<"<3>_{4n,2n,2n} from nested loops = "<<fDirectCorrelations->GetBinContent(8)<<endl;
  cout<<endl;
  cout<<"<3>_{4n,3n,1n} from Q-vectors    = "<<fQCorrelations->GetBinContent(9)<<endl;
  cout<<"<3>_{4n,3n,1n} from nested loops = "<<fDirectCorrelations->GetBinContent(9)<<endl;
  cout<<endl; 
  cout<<"<4>_{1n,1n|1n,1n} from Q-vectors    = "<<fQCorrelations->GetBinContent(11)<<endl;
  cout<<"<4>_{1n,1n|1n,1n} from nested loops = "<<fDirectCorrelations->GetBinContent(11)<<endl;
  cout<<endl;
  cout<<"<4>_{2n,1n|2n,1n} from Q-vectors    = "<<fQCorrelations->GetBinContent(12)<<endl;
  cout<<"<4>_{2n,1n|2n,1n} from nested loops = "<<fDirectCorrelations->GetBinContent(12)<<endl;
  cout<<endl;
  cout<<"<4>_{2n,2n|2n,2n} from Q-vectors    = "<<fQCorrelations->GetBinContent(13)<<endl;
  cout<<"<4>_{2n,2n|2n,2n} from nested loops = "<<fDirectCorrelations->GetBinContent(13)<<endl;
  cout<<endl;
  cout<<"<4>_{3n|1n,1n,1n} from Q-vectors    = "<<fQCorrelations->GetBinContent(14)<<endl;
  cout<<"<4>_{3n|1n,1n,1n} from nested loops = "<<fDirectCorrelations->GetBinContent(14)<<endl;
  cout<<endl;
  cout<<"<4>_{3n,1n|3n,1n} from Q-vectors    = "<<fQCorrelations->GetBinContent(15)<<endl;
  cout<<"<4>_{3n,1n|3n,1n} from nested loops = "<<fDirectCorrelations->GetBinContent(15)<<endl;
  cout<<endl;
  cout<<"<4>_{3n,1n|2n,2n} from Q-vectors    = "<<fQCorrelations->GetBinContent(16)<<endl;
  cout<<"<4>_{3n,1n|2n,2n} from nested loops = "<<fDirectCorrelations->GetBinContent(16)<<endl;
  cout<<endl; 
  cout<<"<4>_{4n|2n,1n,1n} from Q-vectors    = "<<fQCorrelations->GetBinContent(17)<<endl;
  cout<<"<4>_{4n|2n,1n,1n} from nested loops = "<<fDirectCorrelations->GetBinContent(17)<<endl;
  cout<<endl;
  cout<<"<5>_{2n,1n|1n,1n,1n} from Q-vectors    = "<<fQCorrelations->GetBinContent(19)<<endl;
  cout<<"<5>_{2n,1n|1n,1n,1n} from nested loops = "<<fDirectCorrelations->GetBinContent(19)<<endl;
  cout<<endl;
  cout<<"<5>_{2n,2n|2n,1n,1n} from Q-vectors    = "<<fQCorrelations->GetBinContent(20)<<endl;
  cout<<"<5>_{2n,2n|2n,1n,1n} from nested loops = "<<fDirectCorrelations->GetBinContent(20)<<endl;
  cout<<endl;
  cout<<"<5>_{3n,1n|2n,1n,1n} from Q-vectors    = "<<fQCorrelations->GetBinContent(21)<<endl;
  cout<<"<5>_{3n,1n|2n,1n,1n} from nested loops = "<<fDirectCorrelations->GetBinContent(21)<<endl;
  cout<<endl;
  cout<<"<5>_{4n|1n,1n,1n,1n} from Q-vectors    = "<<fQCorrelations->GetBinContent(22)<<endl;
  cout<<"<5>_{4n|1n,1n,1n,1n} from nested loops = "<<fDirectCorrelations->GetBinContent(22)<<endl;
  cout<<endl;
  cout<<"<6>_{1n,1n,1n|1n,1n,1n} from Q-vectors    = "<<fQCorrelations->GetBinContent(24)<<endl;
  cout<<"<6>_{1n,1n,1n|1n,1n,1n} from nested loops = "<<fDirectCorrelations->GetBinContent(24)<<endl;
  cout<<endl; 
  cout<<"<6>_{2n,1n,1n|2n,1n,1n} from Q-vectors    = "<<fQCorrelations->GetBinContent(25)<<endl;
  cout<<"<6>_{2n,1n,1n|2n,1n,1n} from nested loops = "<<fDirectCorrelations->GetBinContent(25)<<endl;
  cout<<endl;
  cout<<"<6>_{2n,2n|1n,1n,1n,1n} from Q-vectors    = "<<fQCorrelations->GetBinContent(26)<<endl;
  cout<<"<6>_{2n,2n|1n,1n,1n,1n} from nested loops = "<<fDirectCorrelations->GetBinContent(26)<<endl;
  cout<<endl; 
  cout<<"<6>_{3n,1n|1n,1n,1n,1n} from Q-vectors    = "<<fQCorrelations->GetBinContent(27)<<endl;
  cout<<"<6>_{3n,1n|1n,1n,1n,1n} from nested loops = "<<fDirectCorrelations->GetBinContent(27)<<endl;
  cout<<endl; 
  cout<<"<7>_{2n,1n,1n|1n,1n,1n,1n} from Q-vectors    = "<<fQCorrelations->GetBinContent(29)<<endl;
  cout<<"<7>_{2n,1n,1n|1n,1n,1n,1n} from nested loops = "<<fDirectCorrelations->GetBinContent(29)<<endl;
  cout<<endl; 
  cout<<"<8>_{1n,1n,1n,1n|1n,1n,1n,1n} from Q-vectors    = "<<fQCorrelations->GetBinContent(31)<<endl;
  cout<<"<8>_{1n,1n,1n,1n|1n,1n,1n,1n} from nested loops = "<<fDirectCorrelations->GetBinContent(31)<<endl;
  cout<<endl; 
 }
 
 if(useWeights)
 {
  //.........................................................................................
  cout<<"<w1 w2 cos(n*(phi1-phi2))> from Q-vectors         = "<<fQCorrelationsW->GetBinContent(1)<<endl;
  cout<<"<<w1 w2 cos(n*(phi1-phi2))> from nested loops     = "<<fDirectCorrelationsW->GetBinContent(1)<<endl;
  cout<<endl;
  cout<<"<w1^2 w2^2 cos(2n*(phi1-phi2))> from Q-vectors    = "<<fQCorrelationsW->GetBinContent(2)<<endl;
  cout<<"<w1^2 w2^2 cos(2n*(phi1-phi2))> from nested loops = "<<fDirectCorrelationsW->GetBinContent(2)<<endl;
  cout<<endl;
  cout<<"<w1^3 w2^3 cos(3n*(phi1-phi2))> from Q-vectors    = "<<fQCorrelationsW->GetBinContent(3)<<endl;
  cout<<"<w1^3 w2^3 cos(3n*(phi1-phi2))> from nested loops = "<<fDirectCorrelationsW->GetBinContent(3)<<endl;
  cout<<endl;
  cout<<"<w1^4 w2^4 cos(4n*(phi1-phi2))> from Q-vectors    = "<<fQCorrelationsW->GetBinContent(4)<<endl;
  cout<<"<w1^4 w2^4 cos(4n*(phi1-phi2))> from nested loops = "<<fDirectCorrelationsW->GetBinContent(4)<<endl;
  cout<<endl;  
  cout<<"<w1^3 w2 cos(n*(phi1-phi2))> from Q-vectors       = "<<fQCorrelationsW->GetBinContent(5)<<endl;
  cout<<"<w1^3 w2 cos(n*(phi1-phi2))> from nested loops    = "<<fDirectCorrelationsW->GetBinContent(5)<<endl;
  cout<<endl;
  cout<<"<w1 w2 w3^2 cos(n*(phi1-phi2))> from Q-vectors    = "<<fQCorrelationsW->GetBinContent(6)<<endl;
  cout<<"<w1 w2 w3^2 cos(n*(phi1-phi2))> from nested loops = "<<fDirectCorrelationsW->GetBinContent(6)<<endl;
  cout<<endl;
  cout<<"<w1^2 w2 w3 cos(n*(2phi1-phi2-phi3))> from Q-vectors    = "<<fQCorrelationsW->GetBinContent(21)<<endl;
  cout<<"<w1^2 w2 w3 cos(n*(2phi1-phi2-phi3))> from nested loops = "<<fDirectCorrelationsW->GetBinContent(21)<<endl;
  cout<<endl;
  cout<<"<w1 w2 w3 w4 cos(n*(phi1+phi2-phi3-phi4))> from Q-vectors    = "<<fQCorrelationsW->GetBinContent(41)<<endl;
  cout<<"<w1 w2 w3 w4 cos(n*(phi1+phi2-phi3-phi4))> from nested loops = "<<fDirectCorrelationsW->GetBinContent(41)<<endl;
  cout<<endl;
  //.........................................................................................
 }
 
} // end of AliFlowAnalysisWithQCumulants::CompareDirectAndQCorrelationsForIntegratedFlow(Bool_t useWeights)


//================================================================================================================================


void AliFlowAnalysisWithQCumulants::CompareDirectAndQCorrelationsForDifferentialFlow(Bool_t useWeights)
{
 // compare correlations needed for diff. flow calculated with nested loops and those calculated from Q-vectors

 cout<<endl;
 cout<<"   *************************************"<<endl;
 cout<<"   **** cross-checking the formulas ****"<<endl;
 cout<<"   ****    for differential flow    ****"<<endl;
 cout<<"   ****                             ****"<<endl;
 cout<<"   ****        (pt,eta) bin:        ****"<<endl; 
 cout<<"   ****    1.1  < pt  <  1.2  GeV   ****"<<endl;  
 cout<<"   ****   -0.55 < eta < -0.525      ****"<<endl; 
 cout<<"   *************************************"<<endl;                             
 cout<<endl;
 
 if(!useWeights)
 {                                       
  cout<<"<cos(n(psi1-phi2))> from Q-vectors    = "<<f2pPtEtaPOI->GetBinContent(f2pPtEtaPOI->GetBin(12,19))<<endl;
  cout<<"<cos(n(psi1-phi2))> from nested loops = "<<fDirectCorrelationsDiffFlow->GetBinContent(1)<<endl;
  cout<<endl;  
  cout<<"<cos(n(psi1+phi2-phi3-phi4))> from Q-vectors    = "<<f4pPtEtaPOI->GetBinContent(f4pPtEtaPOI->GetBin(12,19))<<endl;
  cout<<"<cos(n(psi1+phi2-phi3-phi4))> from nested loops = "<<fDirectCorrelationsDiffFlow->GetBinContent(41)<<endl;
  cout<<endl;  
 }
 
 if(useWeights)
 {
  cout<<"<w2 cos(n(psi1-phi2))> from Q-vectors    = "<<f2pPtEtaPOIW->GetBinContent(f2pPtEtaPOIW->GetBin(12,19))<<endl;
  cout<<"<w2 cos(n(psi1-phi2))> from nested loops = "<<fDirectCorrelationsDiffFlowW->GetBinContent(1)<<endl;
  cout<<endl;  
  cout<<"<w2 w3 w4 cos(n(psi1+phi2-phi3-phi4))> from Q-vectors    = "<<f4pPtEtaPOIW->GetBinContent(f4pPtEtaPOIW->GetBin(12,19))<<endl;
  cout<<"<w2 w3 w4 cos(n(psi1+phi2-phi3-phi4))> from nested loops = "<<fDirectCorrelationsDiffFlowW->GetBinContent(41)<<endl;
  cout<<endl;   
 }
 
} // end of void AliFlowAnalysisWithQCumulants::CompareDirectAndQCorrelationsForDifferentialFlow()


//================================================================================================================================


void AliFlowAnalysisWithQCumulants::WriteHistograms(TString outputFileName)
{
 //store the final results in output .root file
 TFile *output = new TFile(outputFileName.Data(),"RECREATE");
 output->WriteObject(fHistList, "cobjQC","SingleKey");
 delete output;
}


//================================================================================================================================


void AliFlowAnalysisWithQCumulants::TempDeleteMe()
{
 /*
 //xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx 
 TCanvas* qvectorPlot = new TCanvas("qvectorPlot","Q-vector Plot",1000,1000);
 
 qvectorPlot->cd(1);
 
 TH1D* style = new TH1D("style","Q-vectors",100,-244,244); 
 (style->GetYaxis())->SetRangeUser(-244,244);
 
 style->Draw();
  
 Int_t nBins=fQvectorForEachEventX->GetNbinsX();
 Double_t qxxx=0.,qyyy=0.;
 //cout<<"nBins = "<<nBins<<endl;   
 //cout<<fQvectorForEachEventX->GetBinEntries(4)<<endl;    
 //cout<<fQvectorForEachEventY->GetBinEntries(4)<<endl;     
 
 for(Int_t b=1;b<nBins+1;b++)
 {
  if(fQvectorForEachEventX->GetBinEntries(b)==1 && fQvectorForEachEventY->GetBinEntries(b)==1)
  {
   qxxx=fQvectorForEachEventX->GetBinContent(b);
   qyyy=fQvectorForEachEventY->GetBinContent(b);
   //cout<<qxxx<<" "<<qyyy<<endl;
   TArrow *qvector = new TArrow(0.0,0.0,qxxx,qyyy,0.0144,"|>");
   qvector->SetAngle(40);
   qvector->SetLineWidth(2);
   qvector->Draw("");
  }
 }  
 //xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx       
 */                             
                                                                                                                                                              
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    

 //xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx                          
 Int_t nBinsPt=3, nBinsEta=2;
 Double_t ptMin=0., ptMax=3.;
 Double_t etaMin=0., etaMax=2.;
                         
                           
 //avarage of the generating function for integrated flow <G[p][q]>
 TProfile2D *tempPtEta = new TProfile2D("tempPtEta","<2'>(pt,eta)",nBinsPt,ptMin,ptMax,nBinsEta,etaMin,etaMax);
 tempPtEta->SetXTitle("pt");
 tempPtEta->SetYTitle("eta");
 
 // (1,1):
 tempPtEta->Fill(0.5,0.67,0.4,2);
 tempPtEta->Fill(0.1,0.44,0.6,3); 
 
 // (3,1):
 tempPtEta->Fill(2.5,0.01,2.2,4);
 tempPtEta->Fill(2.1,0.74,2.6,3.7); 
  
 //tempPtEta->Fill(2.5,0.5,1,2);
 //tempPtEta->Fill(2.5,1.5,3,1);
 //tempPtEta->Fill(2.6,0.6,7,3);
 //tempPtEta->Fill(2.5,0.5,1,1);

  
 
 TCanvas* tempCanvas = new TCanvas("tempCanvas","tempCanvas",1000,600);
 
  tempCanvas->Divide(1,2);
 
  tempCanvas->cd(1); 
  tempPtEta->Draw("SURF1");
  
  (tempCanvas->cd(2))->Divide(1,2);
  
  tempCanvas->cd(1);  
  //tempPt->Draw();
 
  tempCanvas->cd(2);  
  //tempEta->Draw();
 
 /*
 cout<<tempPtEta->GetBinContent(tempPtEta->GetBin(1,1))<<endl;
 cout<<tempPtEta->GetBinContent(tempPtEta->GetBin(3,1))<<endl;
 cout<<tempEta->GetBinContent(1)<<endl;
 cout<<tempEta->GetBinEntries(1)<<endl;
 cout<<tempEta->GetBinContent(2)<<endl;
 cout<<tempEta->GetBinEntries(2)<<endl; 
 cout<<endl; 
 */
 /*
 cout<<tempPtEta->GetBinContent(3,1)<<endl;
 cout<<tempPtEta->GetBinEntries(tempPtEta->GetBin(3,1))<<endl;
 cout<<endl;
 
 cout<<tempPtEta->GetBinContent(1,2)<<endl;
 cout<<tempPtEta->GetBinEntries(tempPtEta->GetBin(1,2))<<endl;
 cout<<endl;
 
 cout<<"xy"<<endl;
 cout<<tempPt->GetBinContent(1)<<endl;
 //cout<<tempPt->GetBinEntries(1)<<endl;
 cout<<tempPt->GetBinContent(3)<<endl;
 //cout<<tempPt->GetBinEntries(3)<<endl; 
 cout<<endl;                           
                  
 cout<<tempEta->GetBinContent(1)<<endl;    
 //cout<<tempEta->GetBinEntries(1)<<endl;
 cout<<tempEta->GetBinContent(2)<<endl;    
 //cout<<tempEta->GetBinEntries(2)<<endl;
 cout<<endl;                                          
     
 //tempPtEta->Draw("LEGO2");    
 */
      
 //xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx   


}


//================================================================================================================================


