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
 fAvMultIntFlowQC(NULL),
 fQvectorComponents(NULL),
 fIntFlowResultsQC(NULL),
 fDiffFlowResults2ndOrderQC(NULL),
 fDiffFlowResults4thOrderQC(NULL),
 fCovariances(NULL),
 fQvectorForEachEventX(NULL),//to be removed
 fQvectorForEachEventY(NULL),//to be removed
 fQCorrelations(NULL),
 fWeightedQCorrelations(NULL),
 fQProduct(NULL),
 fDirectCorrelations(NULL),
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
 fUseEtaWeights(kFALSE)
{
 //constructor 
 fHistList = new TList();
 fDiffFlowList = new TList(); 
 fDiffFlowList->SetName("DifferentialFlow"); 
 fWeightsList = new TList();
 fWeightsList->SetName("Weights");
  
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
 
 //final results for integrated flow from Q-cumulants
 fIntFlowResultsQC = new TH1D("fIntFlowResultsQC","Integrated Flow from Q-cumulants",4,0,4);
 //fIntFlowResults->SetXTitle("");
 //fIntFlowResultsQC->SetYTitle("IntegFALSrated Flow");
 fIntFlowResultsQC->SetLabelSize(0.06);
 //fIntFlowResultsQC->SetTickLength(1);
 fIntFlowResultsQC->SetMarkerStyle(25);
 (fIntFlowResultsQC->GetXaxis())->SetBinLabel(1,"v_{n}{2}");
 (fIntFlowResultsQC->GetXaxis())->SetBinLabel(2,"v_{n}{4}");
 (fIntFlowResultsQC->GetXaxis())->SetBinLabel(3,"v_{n}{6}");
 (fIntFlowResultsQC->GetXaxis())->SetBinLabel(4,"v_{n}{8}");
 fHistList->Add(fIntFlowResultsQC);

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
    
 //multi-particle correlations calculated from Q-vectors
 fQCorrelations = new TProfile("fQCorrelations","multi-particle correlations from Q-vectors",32,0,32,"s");
 //fQCorrelations->SetXTitle("correlations");
 //fQCorrelations->SetYTitle("");
 fQCorrelations->SetTickLength(-0.01,"Y");
 fQCorrelations->SetMarkerStyle(25);
 fQCorrelations->SetLabelSize(0.03);
 fQCorrelations->SetLabelOffset(0.01,"Y");
 
 (fQCorrelations->GetXaxis())->SetBinLabel(1,"<<2>>_{n|n}");
 (fQCorrelations->GetXaxis())->SetBinLabel(2,"<<2>>_{2n|2n}");
 (fQCorrelations->GetXaxis())->SetBinLabel(3,"<<2>>_{3n|3n}");
 (fQCorrelations->GetXaxis())->SetBinLabel(4,"<<2>>_{4n|4n}");
 
 (fQCorrelations->GetXaxis())->SetBinLabel(6,"<<3>>_{2n|n,n}");
 (fQCorrelations->GetXaxis())->SetBinLabel(7,"<<3>>_{3n|2n,n}");
 (fQCorrelations->GetXaxis())->SetBinLabel(8,"<<3>>_{4n|2n,2n}");
 (fQCorrelations->GetXaxis())->SetBinLabel(9,"<<3>>_{4n|3n,n}");
 
 (fQCorrelations->GetXaxis())->SetBinLabel(11,"<<4>>_{n,n|n,n}"); 
 (fQCorrelations->GetXaxis())->SetBinLabel(12,"<<4>>_{2n,n|2n,n}");
 (fQCorrelations->GetXaxis())->SetBinLabel(13,"<<4>>_{2n,2n|2n,2n}");
 (fQCorrelations->GetXaxis())->SetBinLabel(14,"<<4>>_{3n|n,n,n}");
 (fQCorrelations->GetXaxis())->SetBinLabel(15,"<<4>>_{3n,n|3n,n}");
 (fQCorrelations->GetXaxis())->SetBinLabel(16,"<<4>>_{3n,n|2n,2n}"); 
 (fQCorrelations->GetXaxis())->SetBinLabel(17,"<<4>>_{4n|2n,n,n}");
 
 (fQCorrelations->GetXaxis())->SetBinLabel(19,"<<5>>_{2n|n,n,n,n}"); 
 (fQCorrelations->GetXaxis())->SetBinLabel(20,"<<5>>_{2n,2n|2n,n,n}");
 (fQCorrelations->GetXaxis())->SetBinLabel(21,"<<5>>_{3n,n|2n,n,n}");
 (fQCorrelations->GetXaxis())->SetBinLabel(22,"<<5>>_{4n|n,n,n,n}");
 
 (fQCorrelations->GetXaxis())->SetBinLabel(24,"<<6>>_{n,n,n|n,n,n}");
 (fQCorrelations->GetXaxis())->SetBinLabel(25,"<<6>>_{2n,n,n|2n,n,n}");
 (fQCorrelations->GetXaxis())->SetBinLabel(26,"<<6>>_{2n,2n|n,n,n,n}");
 (fQCorrelations->GetXaxis())->SetBinLabel(27,"<<6>>_{3n,n|n,n,n,n}");
 
 (fQCorrelations->GetXaxis())->SetBinLabel(29,"<<7>>_{2n,n,n|n,n,n,n}");

 (fQCorrelations->GetXaxis())->SetBinLabel(31,"<<8>>_{n,n,n,n|n,n,n,n}");
 
 fHistList->Add(fQCorrelations);
 
 
 
 
 //weighted multi-particle correlations calculated from Q-vectors
 fWeightedQCorrelations = new TProfile("fWeightedQCorrelations","weighted multi-particle correlations from Q-vectors",100,0,100,"s");
 //fWeightedQCorrelations->SetXTitle("correlations");
 //fWeightedQCorrelations->SetYTitle("");
 fWeightedQCorrelations->SetTickLength(-0.01,"Y");
 fWeightedQCorrelations->SetMarkerStyle(25);
 fWeightedQCorrelations->SetLabelSize(0.03);
 fWeightedQCorrelations->SetLabelOffset(0.01,"Y");
 
 (fWeightedQCorrelations->GetXaxis())->SetBinLabel(1,"<w_{1}w_{2}cos(n(#phi_{1}-#phi_{2}))>");
 (fWeightedQCorrelations->GetXaxis())->SetBinLabel(2,"<w_{1}^{2}w_{2}^{2}cos(2n(#phi_{1}-#phi_{2}))>");
 (fWeightedQCorrelations->GetXaxis())->SetBinLabel(3,"<w_{1}^{3}w_{2}^{3}cos(3n(#phi_{1}-#phi_{2}))>");
 (fWeightedQCorrelations->GetXaxis())->SetBinLabel(4,"<w_{1}^{4}w_{2}^{4}cos(4n(#phi_{1}-#phi_{2}))>");
 (fWeightedQCorrelations->GetXaxis())->SetBinLabel(5,"<w_{1}^{3}w_{2}cos(n(#phi_{1}-#phi_{2}))>");
 (fWeightedQCorrelations->GetXaxis())->SetBinLabel(6,"<w_{1}^{2}w_{2}w_{3}cos(n(#phi_{1}-#phi_{2}))>");
 
 (fWeightedQCorrelations->GetXaxis())->SetBinLabel(11,"<w_{1}w_{2}w_{3}^{2}cos(n(2#phi_{1}-#phi_{2}-#phi_{3}))>");
 
 (fWeightedQCorrelations->GetXaxis())->SetBinLabel(21,"<w_{1}w_{2}w_{3}w_{4}cos(n(#phi_{1}+#phi_{2}-#phi_{3}-#phi_{4}))>");
 
 /*
 (fWeightedQCorrelations->GetXaxis())->SetBinLabel(7,"<<3>>_{3n|2n,n}");
 (fWeightedQCorrelations->GetXaxis())->SetBinLabel(8,"<<3>>_{4n|2n,2n}");
 (fWeightedQCorrelations->GetXaxis())->SetBinLabel(9,"<<3>>_{4n|3n,n}");
 
 (fWeightedQCorrelations->GetXaxis())->SetBinLabel(11,"<<4>>_{n,n|n,n}"); 
 (fWeightedQCorrelations->GetXaxis())->SetBinLabel(12,"<<4>>_{2n,n|2n,n}");
 (fWeightedQCorrelations->GetXaxis())->SetBinLabel(13,"<<4>>_{2n,2n|2n,2n}");
 (fWeightedQCorrelations->GetXaxis())->SetBinLabel(14,"<<4>>_{3n|n,n,n}");
 (fWeightedQCorrelations->GetXaxis())->SetBinLabel(15,"<<4>>_{3n,n|3n,n}");
 (fWeightedQCorrelations->GetXaxis())->SetBinLabel(16,"<<4>>_{3n,n|2n,2n}"); 
 (fWeightedQCorrelations->GetXaxis())->SetBinLabel(17,"<<4>>_{4n|2n,n,n}");
 
 (fWeightedQCorrelations->GetXaxis())->SetBinLabel(19,"<<5>>_{2n|n,n,n,n}"); 
 (fWeightedQCorrelations->GetXaxis())->SetBinLabel(20,"<<5>>_{2n,2n|2n,n,n}");
 (fWeightedQCorrelations->GetXaxis())->SetBinLabel(21,"<<5>>_{3n,n|2n,n,n}");
 (fWeightedQCorrelations->GetXaxis())->SetBinLabel(22,"<<5>>_{4n|n,n,n,n}");
 
 (fWeightedQCorrelations->GetXaxis())->SetBinLabel(24,"<<6>>_{n,n,n|n,n,n}");
 (fWeightedQCorrelations->GetXaxis())->SetBinLabel(25,"<<6>>_{2n,n,n|2n,n,n}");
 (fWeightedQCorrelations->GetXaxis())->SetBinLabel(26,"<<6>>_{2n,2n|n,n,n,n}");
 (fWeightedQCorrelations->GetXaxis())->SetBinLabel(27,"<<6>>_{3n,n|n,n,n,n}");
 
 (fWeightedQCorrelations->GetXaxis())->SetBinLabel(29,"<<7>>_{2n,n,n|n,n,n,n}");

 (fWeightedQCorrelations->GetXaxis())->SetBinLabel(31,"<<8>>_{n,n,n,n|n,n,n,n}");
 */
 
 fHistList->Add(fWeightedQCorrelations);
 
 
 
 
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
 
 //weighted multi-particle correlations calculated with nested loops (0..100 integrated flow; 100..200 differential flow)
 fDirectCorrelations = new TProfile("fDirectCorrelations","multi-particle correlations with nested loops",200,0,200,"s");
 fDirectCorrelations->SetXTitle("");
 fDirectCorrelations->SetYTitle("correlations");
 fHistList->Add(fDirectCorrelations);
 
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
 
 // add list fWeightsList with weights to the main list
 fHistList->Add(fWeightsList);
  
 // add list fDiffFlowList with histograms and profiles needed for differential flow to the main list 
 fHistList->Add(fDiffFlowList); 
}//end of Init()

//================================================================================================================

void AliFlowAnalysisWithQCumulants::Make(AliFlowEventSimple* anEvent)
{
 // running over data 
 
 Int_t nPrim = anEvent->NumberOfTracks(); // nPrim = nRP + nPOI + rest  
 Int_t nRP = anEvent->GetEventNSelTracksIntFlow(); // nRP = number of particles used to determine the reaction plane
 
 Int_t n = 2; // int flow harmonic (to be improved)
 
 //needed for debugging: (to be improved - add explanation here) 
 //Bool_t bNestedLoops=kTRUE;
 //if(!(bNestedLoops)||(nPrim>0&&nPrim<12))
 //{
 //if(nPrim>0&&nPrim<10)
 //{
 
 
 
 //---------------------------------------------------------------------------------------------------------
 // non-weighted and weighted Q-vectors of an event built-up from RP particles evaluated in harmonics n, 2n, 3n and 4n:
 AliFlowVector afvQvector1n, afvQvector2n, afvQvector3n, afvQvector4n;

 // non-weighted Q-vector in harmonic n: 
 afvQvector1n.Set(0.,0.);
 afvQvector1n.SetMult(0);
 afvQvector1n = anEvent->GetQ(1*n); 
 
 // non-weighted Q-vector in harmonic 2n: 
 afvQvector2n.Set(0.,0.);
 afvQvector2n.SetMult(0);
 afvQvector2n = anEvent->GetQ(2*n); // to be improved: weights   
          
 // non-weighted Q-vector in harmonic 3n:                                                                 
 afvQvector3n.Set(0.,0.);
 afvQvector3n.SetMult(0);
 afvQvector3n = anEvent->GetQ(3*n); // to be improved: weights
 
 // non-weighted Q-vector in harmonic 4n:
 afvQvector4n.Set(0.,0.);
 afvQvector4n.SetMult(0);
 afvQvector4n = anEvent->GetQ(4*n); // to be improved: weights
            
            
 //xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 //                        !!!! to be removed !!!!
 fQvectorForEachEventX->Fill(1.*(++fEventCounter),afvQvector1n.X());
 fQvectorForEachEventY->Fill(1.*(fEventCounter),afvQvector1n.Y()); 
 //xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
                       
                                  
                                                        
 //---------------------------------------------------------------------------------------------------------
 
 //multiplicity of RP particles:
 Double_t dMult = afvQvector1n.GetMult(); // to be improved (name, this is actually weighted multiplicity)
 
 fAvMultIntFlowQC->Fill(0.,dMult,1.); // to be removed (this info is also stored in one of control histograms)
 
 //---------------------------------------------------------------------------------------------------------
 //
 //                                          *******************
 //                                          **** Q-vectors ****
 //                                          *******************
 //
 Double_t reQ2nQ1nstarQ1nstar = pow(afvQvector1n.X(),2.)*afvQvector2n.X()+2.*afvQvector1n.X()*afvQvector1n.Y()*afvQvector2n.Y()-pow(afvQvector1n.Y(),2.)*afvQvector2n.X();//Re[Q_{2n} Q_{n}^* Q_{n}^*]
 //Double_t imQ2nQ1nstarQ1nstar = pow(Qvector1n.X(),2.)*Qvector2n.Y()-2.*Qvector1n.X()*Qvector1n.Y()*Qvector2n.X()-pow(Qvector1n.Y(),2.)*Qvector2n.Y();//Im[Q_{2n} Q_{n}^* Q_{n}^*]
 Double_t reQ1nQ1nQ2nstar = reQ2nQ1nstarQ1nstar;//Re[Q_{n} Q_{n} Q_{2n}^*] = Re[Q_{2n} Q_{n}^* Q_{n}^*]
 Double_t reQ3nQ1nQ2nstarQ2nstar = (pow(afvQvector2n.X(),2.)-pow(afvQvector2n.Y(),2.))*(afvQvector3n.X()*afvQvector1n.X()-afvQvector3n.Y()*afvQvector1n.Y())+2.*afvQvector2n.X()*afvQvector2n.Y()*(afvQvector3n.X()*afvQvector1n.Y()+afvQvector3n.Y()*afvQvector1n.X());
 //Double_t imQ3nQ1nQ2nstarQ2nstar = calculate and implement this (deleteMe) 
 Double_t reQ2nQ2nQ3nstarQ1nstar = reQ3nQ1nQ2nstarQ2nstar;
 Double_t reQ4nQ2nstarQ2nstar = pow(afvQvector2n.X(),2.)*afvQvector4n.X()+2.*afvQvector2n.X()*afvQvector2n.Y()*afvQvector4n.Y()-pow(afvQvector2n.Y(),2.)*afvQvector4n.X();//Re[Q_{4n} Q_{2n}^* Q_{2n}^*]
 //Double_t imQ4nQ2nstarQ2nstar = calculate and implement this (deleteMe)
 Double_t reQ2nQ2nQ4nstar = reQ4nQ2nstarQ2nstar;
 Double_t reQ4nQ3nstarQ1nstar = afvQvector4n.X()*(afvQvector3n.X()*afvQvector1n.X()-afvQvector3n.Y()*afvQvector1n.Y())+afvQvector4n.Y()*(afvQvector3n.X()*afvQvector1n.Y()+afvQvector3n.Y()*afvQvector1n.X());//Re[Q_{4n} Q_{3n}^* Q_{n}^*]
 Double_t reQ3nQ1nQ4nstar = reQ4nQ3nstarQ1nstar;//Re[Q_{3n} Q_{n} Q_{4n}^*] = Re[Q_{4n} Q_{3n}^* Q_{n}^*]
 //Double_t imQ4nQ3nstarQ1nstar = calculate and implement this (deleteMe)
 Double_t reQ3nQ2nstarQ1nstar = afvQvector3n.X()*afvQvector2n.X()*afvQvector1n.X()-afvQvector3n.X()*afvQvector2n.Y()*afvQvector1n.Y()+afvQvector3n.Y()*afvQvector2n.X()*afvQvector1n.Y()+afvQvector3n.Y()*afvQvector2n.Y()*afvQvector1n.X();//Re[Q_{3n} Q_{2n}^* Q_{n}^*]
 Double_t reQ2nQ1nQ3nstar = reQ3nQ2nstarQ1nstar;//Re[Q_{2n} Q_{n} Q_{3n}^*] = Re[Q_{3n} Q_{2n}^* Q_{n}^*]
 //Double_t imQ3nQ2nstarQ1nstar; //calculate and implement this (deleteMe)
 Double_t reQ3nQ1nstarQ1nstarQ1nstar = afvQvector3n.X()*pow(afvQvector1n.X(),3)-3.*afvQvector1n.X()*afvQvector3n.X()*pow(afvQvector1n.Y(),2)+3.*afvQvector1n.Y()*afvQvector3n.Y()*pow(afvQvector1n.X(),2)-afvQvector3n.Y()*pow(afvQvector1n.Y(),3);//Re[Q_{3n} Q_{n}^* Q_{n}^* Q_{n}^*]
 //Double_t imQ3nQ1nstarQ1nstarQ1nstar; //calculate and implement this (deleteMe)
 Double_t xQ2nQ1nQ2nstarQ1nstar = pow(afvQvector2n.Mod()*afvQvector1n.Mod(),2);//|Q_{2n}|^2 |Q_{n}|^2
 Double_t reQ4nQ2nstarQ1nstarQ1nstar = (afvQvector4n.X()*afvQvector2n.X()+afvQvector4n.Y()*afvQvector2n.Y())*(pow(afvQvector1n.X(),2)-pow(afvQvector1n.Y(),2))+2.*afvQvector1n.X()*afvQvector1n.Y()*(afvQvector4n.Y()*afvQvector2n.X()-afvQvector4n.X()*afvQvector2n.Y());//Re[Q_{4n} Q_{2n}^* Q_{n}^* Q_{n}^*] 
 //Double_t imQ4nQ2nstarQ1nstarQ1nstar; //calculate and implement this (deleteMe)
 Double_t reQ2nQ1nQ1nstarQ1nstarQ1nstar = (afvQvector2n.X()*afvQvector1n.X()-afvQvector2n.Y()*afvQvector1n.Y())*(pow(afvQvector1n.X(),3)-3.*afvQvector1n.X()*pow(afvQvector1n.Y(),2))+(afvQvector2n.X()*afvQvector1n.Y()+afvQvector1n.X()*afvQvector2n.Y())*(3.*afvQvector1n.Y()*pow(afvQvector1n.X(),2)-pow(afvQvector1n.Y(),3));//Re[Q_{2n} Q_{n} Q_{n}^* Q_{n}^* Q_{n}^*]
 //Double_t imQ2nQ1nQ1nstarQ1nstarQ1nstar; //calculate and implement this (deleteMe)
 Double_t reQ2nQ2nQ2nstarQ1nstarQ1nstar = pow(afvQvector2n.Mod(),2.)*(afvQvector2n.X()*(pow(afvQvector1n.X(),2.)-pow(afvQvector1n.Y(),2.))+2.*afvQvector2n.Y()*afvQvector1n.X()*afvQvector1n.Y());//Re[Q_{2n} Q_{2n} Q_{2n}^* Q_{n}^* Q_{n}^*]
 //Double_t imQ2nQ2nQ2nstarQ1nstarQ1nstar = pow(Qvector2n.Mod(),2.)*(Qvector2n.Y()*(pow(Qvector1n.X(),2.)-pow(Qvector1n.Y(),2.))-2.*Qvector2n.X()*Qvector1n.X()*Qvector1n.Y());//Im[Q_{2n} Q_{2n} Q_{2n}^* Q_{n}^* Q_{n}^*]
 Double_t reQ4nQ1nstarQ1nstarQ1nstarQ1nstar = pow(afvQvector1n.X(),4.)*afvQvector4n.X()-6.*pow(afvQvector1n.X(),2.)*afvQvector4n.X()*pow(afvQvector1n.Y(),2.)+pow(afvQvector1n.Y(),4.)*afvQvector4n.X()+4.*pow(afvQvector1n.X(),3.)*afvQvector1n.Y()*afvQvector4n.Y()-4.*pow(afvQvector1n.Y(),3.)*afvQvector1n.X()*afvQvector4n.Y();//Re[Q_{4n} Q_{n}^* Q_{n}^* Q_{n}^* Q_{n}^*]
 //Double_t imQ4nQ1nstarQ1nstarQ1nstarQ1nstar = pow(Qvector1n.X(),4.)*Qvector4n.Y()-6.*pow(Qvector1n.X(),2.)*Qvector4n.Y()*pow(Qvector1n.Y(),2.)+pow(Qvector1n.Y(),4.)*Qvector4n.Y()+4.*pow(Qvector1n.Y(),3.)*Qvector1n.X()*Qvector4n.X()-4.*pow(Qvector1n.X(),3.)*Qvector1n.Y()*Qvector4n.X();//Im[Q_{4n} Q_{n}^* Q_{n}^* Q_{n}^* Q_{n}^*]
 Double_t reQ3nQ1nQ2nstarQ1nstarQ1nstar = pow(afvQvector1n.Mod(),2.)*(afvQvector1n.X()*afvQvector2n.X()*afvQvector3n.X()-afvQvector3n.X()*afvQvector1n.Y()*afvQvector2n.Y()+afvQvector2n.X()*afvQvector1n.Y()*afvQvector3n.Y()+afvQvector1n.X()*afvQvector2n.Y()*afvQvector3n.Y());//Re[Q_{3n} Q_{n} Q_{2n}^* Q_{n}^* Q_{n}^*]
 //Double_t imQ3nQ1nQ2nstarQ1nstarQ1nstar = pow(afvQvector1n.Mod(),2.)*(-afvQvector2n.X()*afvQvector3n.X()*afvQvector1n.Y()-afvQvector1n.X()*afvQvector3n.X()*afvQvector2n.Y()+afvQvector1n.X()*afvQvector2n.X()*afvQvector3n.Y()-afvQvector1n.Y()*afvQvector2n.Y()*afvQvector3n.Y());//Im[Q_{3n} Q_{n} Q_{2n}^* Q_{n}^* Q_{n}^*]
 Double_t reQ2nQ2nQ1nstarQ1nstarQ1nstarQ1nstar = (pow(afvQvector1n.X(),2.)*afvQvector2n.X()-2.*afvQvector1n.X()*afvQvector2n.X()*afvQvector1n.Y()-afvQvector2n.X()*pow(afvQvector1n.Y(),2.)+afvQvector2n.Y()*pow(afvQvector1n.X(),2.)+2.*afvQvector1n.X()*afvQvector1n.Y()*afvQvector2n.Y()-pow(afvQvector1n.Y(),2.)*afvQvector2n.Y())*(pow(afvQvector1n.X(),2.)*afvQvector2n.X()+2.*afvQvector1n.X()*afvQvector2n.X()*afvQvector1n.Y()-afvQvector2n.X()*pow(afvQvector1n.Y(),2.)-afvQvector2n.Y()*pow(afvQvector1n.X(),2.)+2.*afvQvector1n.X()*afvQvector1n.Y()*afvQvector2n.Y()+pow(afvQvector1n.Y(),2.)*afvQvector2n.Y());//Re[Q_{2n} Q_{2n} Q_{n}^* Q_{n}^* Q_{n}^* Q_{n}^*]
 //Double_t imQ2nQ2nQ1nstarQ1nstarQ1nstarQ1nstar = 2.*(pow(afvQvector1n.X(),2.)*afvQvector2n.X()-afvQvector2n.X()*pow(afvQvector1n.Y(),2.)+2.*afvQvector1n.X()*afvQvector1n.Y()*afvQvector2n.Y())*(pow(afvQvector1n.X(),2.)*afvQvector2n.Y()-2.*afvQvector1n.X()*afvQvector1n.Y()*afvQvector2n.X()-pow(afvQvector1n.Y(),2.)*afvQvector2n.Y());//Im[Q_{2n} Q_{2n} Q_{n}^* Q_{n}^* Q_{n}^* Q_{n}^*]
 Double_t reQ3nQ1nQ1nstarQ1nstarQ1nstarQ1nstar = pow(afvQvector1n.Mod(),2.)*(pow(afvQvector1n.X(),3.)*afvQvector3n.X()-3.*afvQvector1n.X()*afvQvector3n.X()*pow(afvQvector1n.Y(),2.)+3.*pow(afvQvector1n.X(),2.)*afvQvector1n.Y()*afvQvector3n.Y()-pow(afvQvector1n.Y(),3.)*afvQvector3n.Y());//Re[Q_{3n} Q_{n} Q_{n}^* Q_{n}^* Q_{n}^* Q_{n}^*]
 //Double_t imQ3nQ1nQ1nstarQ1nstarQ1nstarQ1nstar = pow(afvQvector1n.Mod(),2.)*(pow(afvQvector1n.Y(),3.)*afvQvector3n.X()-3.*afvQvector1n.Y()*afvQvector3n.X()*pow(afvQvector1n.X(),2.)-3.*pow(afvQvector1n.Y(),2.)*afvQvector1n.X()*afvQvector3n.Y()+pow(afvQvector1n.X(),3.)*afvQvector3n.Y());//Im[Q_{3n} Q_{n} Q_{n}^* Q_{n}^* Q_{n}^* Q_{n}^*]
 Double_t xQ2nQ1nQ1nQ2nstarQ1nstarQ1nstar = pow(afvQvector2n.Mod(),2.)*pow(afvQvector1n.Mod(),4.);//|Q_{2n}|^2 |Q_{n}|^4
 Double_t reQ2nQ1nQ1nQ1nstarQ1nstarQ1nstarQ1nstar = pow(afvQvector1n.Mod(),4.)*(pow(afvQvector1n.X(),2.)*afvQvector2n.X()-afvQvector2n.X()*pow(afvQvector1n.Y(),2.)+2.*afvQvector1n.X()*afvQvector1n.Y()*afvQvector2n.Y());//Re[Q_{2n} Q_{n} Q_{n} Q_{n}^* Q_{n}^* Q_{n}^* Q_{n}^*]
 //Double_t imQ2nQ1nQ1nQ1nstarQ1nstarQ1nstarQ1nstar = pow(afvQvector1n.Mod(),4.)*(pow(afvQvector1n.X(),2.)*afvQvector2n.Y()-afvQvector2n.Y()*pow(afvQvector1n.Y(),2.)-2.*afvQvector1n.X()*afvQvector2n.X()*afvQvector1n.Y());//Re[Q_{2n} Q_{n} Q_{n} Q_{n}^* Q_{n}^* Q_{n}^* Q_{n}^*]
 //---------------------------------------------------------------------------------------------------------
 
 //---------------------------------------------------------------------------------------------------------
 //
 //                                        **************************************
 //                                        **** multi-particle correlations: ****
 //                                        **************************************
 //
 // Remark 1: multi-particle correlations calculated with Q-vectors are stored in fQCorrelations.
 // Remark 2: binning of fQCorrelations is organized as follows: 
 //
 // 1st bin: <2>_{n|n} = two1n1n
 // 2nd bin: <2>_{2n|2n} = two2n2n
 // 3rd bin: <2>_{3n|3n} = two3n3n
 // 4th bin: <2>_{4n|4n} = two4n4n
 // 5th bin: --  EMPTY --
 // 6th bin: <3>_{2n|n,n} = three2n1n1n
 // 7th bin: <3>_{3n|2n,n} = three3n2n1n
 // 8th bin: <3>_{4n|2n,2n} = three4n2n2n
 // 9th bin: <3>_{4n|3n,n} = three4n3n1n
 //10th bin: --  EMPTY --
 //11th bin: <4>_{n,n|n,n} = four1n1n1n1n
 //12th bin: <4>_{2n,n|2n,n} = four2n1n2n1n
 //13th bin: <4>_{2n,2n|2n,2n} = four2n2n2n2n
 //14th bin: <4>_{3n|n,n,n} = four3n1n1n1n
 //15th bin: <4>_{3n,n|3n,n} = four3n1n3n1n 
 //16th bin: <4>_{3n,n|2n,2n} = four3n1n2n2n
 //17th bin: <4>_{4n|2n,n,n} = four4n2n1n1n
 //18th bin: --  EMPTY --
 //19th bin: <5>_{2n|n,n,n,n} = five2n1n1n1n1n
 //20th bin: <5>_{2n,2n|2n,n,n} = five2n2n2n1n1n
 //21st bin: <5>_{3n,n|2n,n,n} = five3n1n2n1n1n
 //22nd bin: <5>_{4n|n,n,n,n} = five4n1n1n1n1n  
 //23rd bin: --  EMPTY --
 //24th bin: <6>_{n,n,n|n,n,n} = six1n1n1n1n1n1n
 //25th bin: <6>_{2n,n,n|2n,n,n} = six2n1n1n2n1n1n
 //26th bin: <6>_{2n,2n|n,n,n,n} = six2n2n1n1n1n1n
 //27th bin: <6>_{3n,n|n,n,n,n} = six3n1n1n1n1n1n
 //28th bin: --  EMPTY --
 //29th bin: <7>_{2n,n,n|n,n,n,n} = seven2n1n1n1n1n1n1n
 //30th bin: --  EMPTY --
 //31st bin: <8>_{n,n,n,n|n,n,n,n} = eight1n1n1n1n1n1n1n1n

 
 // binning of fQProduct (all correlations are evaluated in harmonic n): 
 // 1st bin: <2>*<4>
 // 2nd bin: <2>*<6>
 // 3rd bin: <2>*<8> 
 // 4th bin: <4>*<6>
 // 5th bin: <4>*<8>
 // 6th bin: <6>*<8>
         
 // 2-particle
 Double_t two1n1n = 0., two2n2n = 0., two3n3n = 0., two4n4n = 0.; 
 
 if(dMult>1)
 {
  //fill the common control histogram (2nd order): 
  fCommonHists2nd->FillControlHistograms(anEvent); 
  
  two1n1n = (pow(afvQvector1n.Mod(),2.)-dMult)/(dMult*(dMult-1.)); // <2>_{n|n} = <cos(n*(phi1-phi2))>
  two2n2n = (pow(afvQvector2n.Mod(),2.)-dMult)/(dMult*(dMult-1.)); // <2>_{2n|2n} = <cos(2n*(phi1-phi2))>
  two3n3n = (pow(afvQvector3n.Mod(),2.)-dMult)/(dMult*(dMult-1.)); // <2>_{3n|3n} = <cos(3n*(phi1-phi2))>
  two4n4n = (pow(afvQvector4n.Mod(),2.)-dMult)/(dMult*(dMult-1.)); // <2>_{4n|4n} = <cos(4n*(phi1-phi2))>
    
  fQCorrelations->Fill(0.,two1n1n,dMult*(dMult-1.));  
  fQCorrelations->Fill(1.,two2n2n,dMult*(dMult-1.)); 
  fQCorrelations->Fill(2.,two3n3n,dMult*(dMult-1.)); 
  fQCorrelations->Fill(3.,two4n4n,dMult*(dMult-1.)); 
  
  f2pDistribution->Fill(two1n1n,dMult*(dMult-1.)); 
 }
 
 // 3-particle
 Double_t three2n1n1n=0., three3n2n1n=0., three4n2n2n=0., three4n3n1n=0.;
 if(dMult>2)
 {
  three2n1n1n = (reQ2nQ1nstarQ1nstar-2.*pow(afvQvector1n.Mod(),2.)-pow(afvQvector2n.Mod(),2.)+2.*dMult)/(dMult*(dMult-1.)*(dMult-2.)); //Re[<3>_{2n|n,n}] = Re[<3>_{n,n|2n}] = <cos(n*(2.*phi1-phi2-phi3))>
  three3n2n1n = (reQ3nQ2nstarQ1nstar-pow(afvQvector3n.Mod(),2.)-pow(afvQvector2n.Mod(),2.)-pow(afvQvector1n.Mod(),2.)+2.*dMult)/(dMult*(dMult-1.)*(dMult-2.)); //Re[<3>_{3n|2n,n}] = Re[<3>_{2n,n|3n}] = <cos(n*(3.*phi1-2.*phi2-phi3))>
  three4n2n2n = (reQ4nQ2nstarQ2nstar-2.*pow(afvQvector2n.Mod(),2.)-pow(afvQvector4n.Mod(),2.)+2.*dMult)/(dMult*(dMult-1.)*(dMult-2.)); //Re[<3>_{4n|2n,2n}] = Re[<3>_{2n,2n|4n}] = <cos(n*(4.*phi1-2.*phi2-2.*phi3))>
  three4n3n1n = (reQ4nQ3nstarQ1nstar-pow(afvQvector4n.Mod(),2.)-pow(afvQvector3n.Mod(),2.)-pow(afvQvector1n.Mod(),2.)+2.*dMult)/(dMult*(dMult-1.)*(dMult-2.)); //Re[<3>_{4n|3n,n}] = Re[<3>_{3n,n|4n}] = <cos(n*(4.*phi1-3.*phi2-phi3))>
 
  fQCorrelations->Fill(5.,three2n1n1n,dMult*(dMult-1.)*(dMult-2.)); 
  fQCorrelations->Fill(6.,three3n2n1n,dMult*(dMult-1.)*(dMult-2.));
  fQCorrelations->Fill(7.,three4n2n2n,dMult*(dMult-1.)*(dMult-2.)); 
  fQCorrelations->Fill(8.,three4n3n1n,dMult*(dMult-1.)*(dMult-2.));    
 }
 
 //4-particle
 Double_t four1n1n1n1n=0., four2n2n2n2n=0., four2n1n2n1n=0., four3n1n1n1n=0., four4n2n1n1n=0., four3n1n2n2n=0., four3n1n3n1n=0.;  
 if(dMult>3)
 {
  //fill the common control histogram (4th order):
  fCommonHists4th->FillControlHistograms(anEvent); 

  four1n1n1n1n = (2.*dMult*(dMult-3.)+pow(afvQvector1n.Mod(),4.)-4.*(dMult-2.)*pow(afvQvector1n.Mod(),2.)-2.*reQ2nQ1nstarQ1nstar+pow(afvQvector2n.Mod(),2.))/(dMult*(dMult-1)*(dMult-2.)*(dMult-3.));//<4>_{n,n|n,n}
  four2n2n2n2n = (2.*dMult*(dMult-3.)+pow(afvQvector2n.Mod(),4.)-4.*(dMult-2.)*pow(afvQvector2n.Mod(),2.)-2.*reQ4nQ2nstarQ2nstar+pow(afvQvector4n.Mod(),2.))/(dMult*(dMult-1)*(dMult-2.)*(dMult-3.));//<4>_{2n,2n|2n,2n}
  four2n1n2n1n = (xQ2nQ1nQ2nstarQ1nstar-2.*reQ3nQ2nstarQ1nstar-2.*reQ2nQ1nstarQ1nstar)/(dMult*(dMult-1.)*(dMult-2.)*(dMult-3.))-((dMult-5.)*pow(afvQvector1n.Mod(),2.)+(dMult-4.)*pow(afvQvector2n.Mod(),2.)-pow(afvQvector3n.Mod(),2.))/(dMult*(dMult-1.)*(dMult-2.)*(dMult-3.))+(dMult-6.)/((dMult-1.)*(dMult-2.)*(dMult-3.));//Re[<4>_{2n,n|2n,n}]
  four3n1n1n1n = (reQ3nQ1nstarQ1nstarQ1nstar-3.*reQ3nQ2nstarQ1nstar-3.*reQ2nQ1nstarQ1nstar)/(dMult*(dMult-1.)*(dMult-2.)*(dMult-3.))+(2.*pow(afvQvector3n.Mod(),2.)+3.*pow(afvQvector2n.Mod(),2.)+6.*pow(afvQvector1n.Mod(),2.)-6.*dMult)/(dMult*(dMult-1.)*(dMult-2.)*(dMult-3.));//Re[<4>_{3n|n,n,n}]
  four4n2n1n1n = (reQ4nQ2nstarQ1nstarQ1nstar-2.*reQ4nQ3nstarQ1nstar-reQ4nQ2nstarQ2nstar-2.*reQ3nQ2nstarQ1nstar)/(dMult*(dMult-1.)*(dMult-2.)*(dMult-3.))-(reQ2nQ1nstarQ1nstar-2.*pow(afvQvector4n.Mod(),2.)-2.*pow(afvQvector3n.Mod(),2.)-3.*pow(afvQvector2n.Mod(),2.)-4.*pow(afvQvector1n.Mod(),2.))/(dMult*(dMult-1.)*(dMult-2.)*(dMult-3.))-(6.)/((dMult-1.)*(dMult-2.)*(dMult-3.));//Re[<4>_{4n|2n,n,n}]
  four3n1n2n2n = (reQ3nQ1nQ2nstarQ2nstar-reQ4nQ2nstarQ2nstar-reQ3nQ1nQ4nstar-2.*reQ3nQ2nstarQ1nstar)/(dMult*(dMult-1.)*(dMult-2.)*(dMult-3.))-(2.*reQ1nQ1nQ2nstar-pow(afvQvector4n.Mod(),2.)-2.*pow(afvQvector3n.Mod(),2.)-4.*pow(afvQvector2n.Mod(),2.)-4.*pow(afvQvector1n.Mod(),2.))/(dMult*(dMult-1.)*(dMult-2.)*(dMult-3.))-(6.)/((dMult-1.)*(dMult-2.)*(dMult-3.));//Re[<4>_{3n,n|2n,2n}] 
  four3n1n3n1n = (pow(afvQvector3n.Mod(),2.)*pow(afvQvector1n.Mod(),2.)-2.*reQ4nQ3nstarQ1nstar-2.*reQ3nQ2nstarQ1nstar)/(dMult*(dMult-1.)*(dMult-2.)*(dMult-3.))+(pow(afvQvector4n.Mod(),2.)-(dMult-4.)*pow(afvQvector3n.Mod(),2.)+pow(afvQvector2n.Mod(),2.)-(dMult-4.)*pow(afvQvector1n.Mod(),2.))/(dMult*(dMult-1.)*(dMult-2.)*(dMult-3.))+(dMult-6.)/((dMult-1.)*(dMult-2.)*(dMult-3.));//<4>_{3n,n|3n,n}
  //four_3n1n3n1n = Q3nQ1nQ3nstarQ1nstar/(M*(M-1.)*(M-2.)*(M-3.))-(2.*three_3n2n1n+2.*three_4n3n1n)/(M-3.)-(two_4n4n+M*two_3n3n+two_2n2n+M*two_1n1n)/((M-2.)*(M-3.))-M/((M-1.)*(M-2.)*(M-3.));//<4>_{3n,n|3n,n}
  
  fQCorrelations->Fill(10.,four1n1n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.));
  fQCorrelations->Fill(11.,four2n1n2n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.));
  fQCorrelations->Fill(12.,four2n2n2n2n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.));
  fQCorrelations->Fill(13.,four3n1n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.));
  fQCorrelations->Fill(14.,four3n1n3n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.));
  fQCorrelations->Fill(15.,four3n1n2n2n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.));  
  fQCorrelations->Fill(16.,four4n2n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)); 
  
  f4pDistribution->Fill(four1n1n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.));
  
  fQProduct->Fill(0.,two1n1n*four1n1n1n1n,dMult*(dMult-1.)*dMult*(dMult-1.)*(dMult-2.)*(dMult-3.));
 }

 //5-particle
 Double_t five2n1n1n1n1n=0., five2n2n2n1n1n=0., five3n1n2n1n1n=0., five4n1n1n1n1n=0.;
 if(dMult>4)
 {
  five2n1n1n1n1n = (reQ2nQ1nQ1nstarQ1nstarQ1nstar-reQ3nQ1nstarQ1nstarQ1nstar+6.*reQ3nQ2nstarQ1nstar)/(dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.))-(reQ2nQ1nQ3nstar+3.*(dMult-6.)*reQ2nQ1nstarQ1nstar+3.*reQ1nQ1nQ2nstar)/(dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.))-(2.*pow(afvQvector3n.Mod(),2.)+3.*pow(afvQvector2n.Mod()*afvQvector1n.Mod(),2.)-3.*(dMult-4.)*pow(afvQvector2n.Mod(),2.))/(dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.))-3.*(pow(afvQvector1n.Mod(),4.)-2.*(2*dMult-5.)*pow(afvQvector1n.Mod(),2.)+2.*dMult*(dMult-4.))/(dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.));//Re[<5>_{2n,n|n,n,n}]
  
  five2n2n2n1n1n = (reQ2nQ2nQ2nstarQ1nstarQ1nstar-reQ4nQ2nstarQ1nstarQ1nstar-2.*reQ2nQ2nQ3nstarQ1nstar)/(dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.))+2.*(reQ4nQ2nstarQ2nstar+4.*reQ3nQ2nstarQ1nstar+reQ3nQ1nQ4nstar)/(dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.))+(reQ2nQ2nQ4nstar-2.*(dMult-5.)*reQ2nQ1nstarQ1nstar+2.*reQ1nQ1nQ2nstar)/(dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.))-(2.*pow(afvQvector4n.Mod(),2.)+4.*pow(afvQvector3n.Mod(),2.)+1.*pow(afvQvector2n.Mod(),4.)-2.*(3.*dMult-10.)*pow(afvQvector2n.Mod(),2.))/(dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.))-(4.*pow(afvQvector1n.Mod(),2.)*pow(afvQvector2n.Mod(),2.)-4.*(dMult-5.)*pow(afvQvector1n.Mod(),2.)+4.*dMult*(dMult-6.))/(dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.));//Re[<5>_{2n,2n|2n,n,n}]  

  //five_2n2n2n1n1n = reQ2nQ2nQ2nstarQ1nstarQ1nstar/(M*(M-1.)*(M-2.)*(M-3.)*(M-4.))-(4.*four_2n1n2n1n+2.*four_3n1n2n2n+1.*four_2n2n2n2n+four_4n2n1n1n)/(M-4.)-(2.*three_4n3n1n+three_4n2n2n+three_4n2n2n+2.*three_3n2n1n)/((M-3.)*(M-4.))-(4.*three_3n2n1n+(2.*M-1.)*three_2n1n1n+2.*three_2n1n1n)/((M-3.)*(M-4.))-(two_4n4n+2.*two_3n3n+4.*(M-1.)*two_2n2n+2.*(2.*M-1.)*two_1n1n)/((M-2.)*(M-3.)*(M-4.))-(2.*M-1.)/((M-1.)*(M-2.)*(M-3.)*(M-4.)); //OK! 
   
  five4n1n1n1n1n = (reQ4nQ1nstarQ1nstarQ1nstarQ1nstar-6.*reQ4nQ2nstarQ1nstarQ1nstar-4.*reQ3nQ1nstarQ1nstarQ1nstar)/(dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.))+(8.*reQ4nQ3nstarQ1nstar+3.*reQ4nQ2nstarQ2nstar+12.*reQ3nQ2nstarQ1nstar+12.*reQ2nQ1nstarQ1nstar)/(dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.))-(6.*pow(afvQvector4n.Mod(),2.)+8.*pow(afvQvector3n.Mod(),2.)+12.*pow(afvQvector2n.Mod(),2.)+24.*pow(afvQvector1n.Mod(),2.)-24.*dMult)/(dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.));//Re[<5>_{4n|n,n,n,n}] 
  
  //five_4n1n1n1n1n = reQ4nQ1nstarQ1nstarQ1nstarQ1nstar/(M*(M-1.)*(M-2.)*(M-3.)*(M-4.)) -  (4.*four_3n1n1n1n+6.*four_4n2n1n1n)/(M-4.)  -  (6.*three_2n1n1n  + 12.*three_3n2n1n + 4.*three_4n3n1n + 3.*three_4n2n2n)/((M-3.)*(M-4.))  -  (4.*two_1n1n + 6.*two_2n2n + 4.*two_3n3n + 1.*two_4n4n)/((M-2.)*(M-3.)*(M-4.)) - 1./((M-1.)*(M-2.)*(M-3.)*(M-4.)); //OK!
  
  five3n1n2n1n1n = (reQ3nQ1nQ2nstarQ1nstarQ1nstar-reQ4nQ2nstarQ1nstarQ1nstar-reQ3nQ1nstarQ1nstarQ1nstar)/(dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.))-(reQ3nQ1nQ2nstarQ2nstar-3.*reQ4nQ3nstarQ1nstar-reQ4nQ2nstarQ2nstar)/(dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.))-((2.*dMult-13.)*reQ3nQ2nstarQ1nstar-reQ3nQ1nQ4nstar-9.*reQ2nQ1nstarQ1nstar)/(dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.))-(2.*reQ1nQ1nQ2nstar+2.*pow(afvQvector4n.Mod(),2.)-2.*(dMult-5.)*pow(afvQvector3n.Mod(),2.)+2.*pow(afvQvector3n.Mod(),2.)*pow(afvQvector1n.Mod(),2.))/(dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.))+(2.*(dMult-6.)*pow(afvQvector2n.Mod(),2.)-2.*pow(afvQvector2n.Mod(),2.)*pow(afvQvector1n.Mod(),2.)-pow(afvQvector1n.Mod(),4.)+2.*(3.*dMult-11.)*pow(afvQvector1n.Mod(),2.))/(dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.))-4.*(dMult-6.)/((dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.));//Re[<5>_{3n,n|2n,n,n}] 
  
  //five3n1n2n1n1n = reQ3nQ1nQ2nstarQ1nstarQ1nstar/(dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)) -  (four4n2n1n1n+four1n1n1n1n+four3n1n1n1n+2.*four2n1n2n1n+2.*three3n2n1n+2.*four3n1n3n1n+four3n1n2n2n)/(dMult-4.)  -   (2.*three4n3n1n+three4n2n2n+6.*three3n2n1n+three4n3n1n+2.*three3n2n1n+3.*three2n1n1n+2.*three2n1n1n+4.*two1n1n+2.*two2n2n+2.*two3n3n)/((dMult-3.)*(dMult-4.))  -  (5.*two1n1n + 4.*two2n2n + 3.*two3n3n + 1.*two4n4n + 2.)/((dMult-2.)*(dMult-3.)*(dMult-4.))  - 1./((dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)) ;//Re[<5>_{3n,n|2n,n,n}] //OK!
  
  //five3n1n2n1n1n = reQ3nQ1nQ2nstarQ1nstarQ1nstar/(dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)) -  (four4n2n1n1n+four1n1n1n1n+four3n1n1n1n+2.*four2n1n2n1n+2.*four3n1n3n1n+four3n1n2n2n)/(dMult-4.)  -      (2.*three4n3n1n+three4n2n2n+2.*dMult*three3n2n1n+three4n3n1n+2.*three3n2n1n+3.*three2n1n1n+2.*three2n1n1n)/((dMult-3.)*(dMult-4.))  -  ((4.*dMult-3.)*two1n1n + 2.*dMult*two2n2n + (2.*dMult-1.)*two3n3n + two4n4n)/((dMult-2.)*(dMult-3.)*(dMult-4.))  - (2.*dMult-1.)/((dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)) ;//Re[<5>_{3n,n|2n,n,n}] //OK!
   
  fQCorrelations->Fill(18.,five2n1n1n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)); 
  fQCorrelations->Fill(19.,five2n2n2n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.));
  fQCorrelations->Fill(20.,five3n1n2n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.));
  fQCorrelations->Fill(21.,five4n1n1n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.));
 }

 //xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 //      !!!! to be removed: temporary fix !!!!
 if(dMult>1)
 {
  two1n1n = (pow(afvQvector1n.Mod(),2.)-dMult)/(dMult*(dMult-1.));
 }  
 //xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    
 //6-particle
 Double_t six1n1n1n1n1n1n=0., six2n2n1n1n1n1n=0., six3n1n1n1n1n1n=0., six2n1n1n2n1n1n=0.;
 if(dMult>5)
 {
  //fill the common control histogram (6th order):
  fCommonHists6th->FillControlHistograms(anEvent); 
 
  six1n1n1n1n1n1n = (pow(afvQvector1n.Mod(),6.)+9.*xQ2nQ1nQ2nstarQ1nstar-6.*reQ2nQ1nQ1nstarQ1nstarQ1nstar)/(dMult*(dMult-1)*(dMult-2)*(dMult-3)*(dMult-4)*(dMult-5))+4.*(reQ3nQ1nstarQ1nstarQ1nstar-3.*reQ3nQ2nstarQ1nstar)/(dMult*(dMult-1)*(dMult-2)*(dMult-3)*(dMult-4)*(dMult-5))+2.*(9.*(dMult-4.)*reQ2nQ1nstarQ1nstar+2.*pow(afvQvector3n.Mod(),2.))/(dMult*(dMult-1)*(dMult-2)*(dMult-3)*(dMult-4)*(dMult-5))-9.*(pow(afvQvector1n.Mod(),4.)+pow(afvQvector2n.Mod(),2.))/(dMult*(dMult-1)*(dMult-2)*(dMult-3)*(dMult-5))+(18.*pow(afvQvector1n.Mod(),2.))/(dMult*(dMult-1)*(dMult-3)*(dMult-4))-(6.)/((dMult-1)*(dMult-2)*(dMult-3));//<6>_{n,n,n|n,n,n}
  
  six2n1n1n2n1n1n = (xQ2nQ1nQ1nQ2nstarQ1nstarQ1nstar-dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)*(2.*five2n2n2n1n1n+4.*five2n1n1n1n1n+4.*five3n1n2n1n1n+4.*four2n1n2n1n+1.*four1n1n1n1n)-dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(4.*four1n1n1n1n+4.*two1n1n+2.*three2n1n1n+2.*three2n1n1n+4.*four3n1n1n1n+8.*three2n1n1n+2.*four4n2n1n1n+4.*four2n1n2n1n+2.*two2n2n+8.*four2n1n2n1n+4.*four3n1n3n1n+8.*three3n2n1n+4.*four3n1n2n2n+4.*four1n1n1n1n+4.*four2n1n2n1n+1.*four2n2n2n2n)-dMult*(dMult-1.)*(dMult-2.)*(2.*three2n1n1n+8.*two1n1n+4.*two1n1n+2.+4.*two1n1n+4.*three2n1n1n+2.*two2n2n+4.*three2n1n1n+8.*three3n2n1n+8.*two2n2n+4.*three4n3n1n+4.*two3n3n+4.*three3n2n1n+4.*two1n1n+8.*three2n1n1n+4.*two1n1n+4.*three3n2n1n+4.*three2n1n1n+2.*two2n2n+4.*three3n2n1n+2.*three4n2n2n)-dMult*(dMult-1.)*(4.*two1n1n+4.+4.*two1n1n+2.*two2n2n+1.+4.*two1n1n+4.*two2n2n+4.*two3n3n+   1.+2.*two2n2n+1.*two4n4n)-dMult)/(dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)*(dMult-5.));//<6>_{2n,n,n|2n,n,n}
 
  six2n2n1n1n1n1n = (reQ2nQ2nQ1nstarQ1nstarQ1nstarQ1nstar-dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)*(five4n1n1n1n1n+8.*five2n1n1n1n1n+6.*five2n2n2n1n1n)-dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(4.*four3n1n1n1n+6.*four4n2n1n1n+12.*three2n1n1n+12.*four1n1n1n1n+24.*four2n1n2n1n+4.*four3n1n2n2n+3.*four2n2n2n2n)-dMult*(dMult-1.)*(dMult-2.)*(6.*three2n1n1n+12.*three3n2n1n+4.*three4n3n1n+3.*three4n2n2n+8.*three2n1n1n+24.*two1n1n+12.*two2n2n+12.*three2n1n1n+8.*three3n2n1n+1.*three4n2n2n)-dMult*(dMult-1.)*(4.*two1n1n+6.*two2n2n+4.*two3n3n+1.*two4n4n+2.*two2n2n+8.*two1n1n+6.)-dMult)/(dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)*(dMult-5.));//<6>_{2n,2n,n|n,n,n}
   
  six3n1n1n1n1n1n = (reQ3nQ1nQ1nstarQ1nstarQ1nstarQ1nstar-dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)*(five4n1n1n1n1n+4.*five2n1n1n1n1n+6.*five3n1n2n1n1n+4.*four3n1n1n1n)-dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(4.*four3n1n1n1n+6.*four4n2n1n1n+6.*four1n1n1n1n+12.*three2n1n1n+12.*four2n1n2n1n+6.*four3n1n1n1n+12.*three3n2n1n+4.*four3n1n3n1n+3.*four3n1n2n2n)-dMult*(dMult-1.)*(dMult-2.)*(6.*three2n1n1n+12.*three3n2n1n+4.*three4n3n1n+3.*three4n2n2n+4.*two1n1n+12.*two1n1n+6.*three2n1n1n+12.*three2n1n1n+4.*three3n2n1n+12.*two2n2n+4.*three3n2n1n+4.*two3n3n+1.*three4n3n1n+6.*three3n2n1n)-dMult*(dMult-1.)*(4.*two1n1n+6.*two2n2n+4.*two3n3n+1.*two4n4n+1.*two1n1n+4.+6.*two1n1n+4.*two2n2n+1.*two3n3n)-dMult)/(dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)*(dMult-5.));//<6>_{3n,n|n,n,n,n}
   
  fQCorrelations->Fill(23.,six1n1n1n1n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)*(dMult-5.)); 
  fQCorrelations->Fill(24.,six2n1n1n2n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)*(dMult-5.)); 
  fQCorrelations->Fill(25.,six2n2n1n1n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)*(dMult-5.));
  fQCorrelations->Fill(26.,six3n1n1n1n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)*(dMult-5.)); 

  f6pDistribution->Fill(six1n1n1n1n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)*(dMult-5.)); 
  
  fQProduct->Fill(1.,two1n1n*six1n1n1n1n1n1n,dMult*(dMult-1.)*dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)*(dMult-5.));
  fQProduct->Fill(3.,four1n1n1n1n*six1n1n1n1n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)*(dMult-5.));
 }
 
 //7-particle
 Double_t seven2n1n1n1n1n1n1n=0.;
 if(dMult>6)
 {
  seven2n1n1n1n1n1n1n = (reQ2nQ1nQ1nQ1nstarQ1nstarQ1nstarQ1nstar-dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)*(dMult-5.)*(2.*six3n1n1n1n1n1n+4.*six1n1n1n1n1n1n+1.*six2n2n1n1n1n1n+6.*six2n1n1n2n1n1n+8.*five2n1n1n1n1n)-dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)*(1.*five4n1n1n1n1n +8.*five2n1n1n1n1n+8.*four3n1n1n1n+12.*five3n1n2n1n1n+4.*five2n1n1n1n1n+3.*five2n2n2n1n1n+6.*five2n2n2n1n1n+6.*four1n1n1n1n+24.*four1n1n1n1n+12.*five2n1n1n1n1n+12.*five2n1n1n1n1n+12.*three2n1n1n+24.*four2n1n2n1n+4.*five3n1n2n1n1n+4.*five2n1n1n1n1n)-dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(4.*four3n1n1n1n+6.*four4n2n1n1n+12.*four1n1n1n1n+24.*three2n1n1n+24.*four2n1n2n1n+12.*four3n1n1n1n+24.*three3n2n1n+8.*four3n1n3n1n+6.*four3n1n2n2n+6.*three2n1n1n+12.*four1n1n1n1n+12.*four2n1n2n1n+6.*three2n1n1n+12.*four2n1n2n1n+4.*four3n1n2n2n+3.*four2n2n2n2n+4.*four1n1n1n1n+6.*three2n1n1n+24.*two1n1n+24.*four1n1n1n1n+4.*four3n1n1n1n+24.*two1n1n+24.*three2n1n1n+12.*two2n2n+24.*three2n1n1n+12.*four2n1n2n1n+8.*three3n2n1n+8.*four2n1n2n1n+1.*four4n2n1n1n)-dMult*(dMult-1.)*(dMult-2.)*(6.*three2n1n1n+1.*three2n1n1n+8.*two1n1n+12.*three3n2n1n+24.*two1n1n+12.*three2n1n1n+4.*three2n1n1n+8.*two1n1n+4.*three4n3n1n+24.*three2n1n1n+8.*three3n2n1n+12.*two1n1n+12.*two1n1n+3.*three4n2n2n+24.*two2n2n+6.*two2n2n+12.+12.*three3n2n1n+8.*two3n3n+12.*three2n1n1n+24.*two1n1n+4.*three3n2n1n+8.*three3n2n1n+2.*three4n3n1n+12.*two1n1n+8.*three2n1n1n+4.*three2n1n1n+2.*three3n2n1n+6.*two2n2n+8.*two2n2n+1.*three4n2n2n+4.*three3n2n1n+6.*three2n1n1n)-dMult*(dMult-1.)*(4.*two1n1n+2.*two1n1n+6.*two2n2n+8.+1.*two2n2n+4.*two3n3n+12.*two1n1n+4.*two1n1n+1.*two4n4n+8.*two2n2n+6.+2.*two3n3n+4.*two1n1n+1.*two2n2n)-dMult)/(dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)*(dMult-5.)*(dMult-6.));
        
  fQCorrelations->Fill(28.,seven2n1n1n1n1n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)*(dMult-5.)*(dMult-6.));
 }
 
 //8-particle
 Double_t eight1n1n1n1n1n1n1n1n=0.;
 if(dMult>7)
 {
  //fill the common control histogram (8th order):
  fCommonHists8th->FillControlHistograms(anEvent); 
  
  eight1n1n1n1n1n1n1n1n = (pow(afvQvector1n.Mod(),8)-dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)*(dMult-5.)*(dMult-6.)*(12.*seven2n1n1n1n1n1n1n+16.*six1n1n1n1n1n1n)-dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)*(dMult-5.)*(8.*six3n1n1n1n1n1n+48.*six1n1n1n1n1n1n+6.*six2n2n1n1n1n1n+96.*five2n1n1n1n1n+72.*four1n1n1n1n+36.*six2n1n1n2n1n1n)-dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)*(2.*five4n1n1n1n1n+32.*five2n1n1n1n1n+36.*four1n1n1n1n+32.*four3n1n1n1n+48.*five2n1n1n1n1n+48.*five3n1n2n1n1n+144.*five2n1n1n1n1n+288.*four1n1n1n1n+36.*five2n2n2n1n1n+144.*three2n1n1n+96.*two1n1n+144.*four2n1n2n1n)-dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(8.*four3n1n1n1n+48.*four1n1n1n1n+12.*four4n2n1n1n+96.*four2n1n2n1n+96.*three2n1n1n+72.*three2n1n1n+144.*two1n1n+16.*four3n1n3n1n+48.*four3n1n1n1n+144.*four1n1n1n1n+72.*four1n1n1n1n+96.*three3n2n1n+24.*four3n1n2n2n+144.*four2n1n2n1n+288.*two1n1n+288.*three2n1n1n+9.*four2n2n2n2n+72.*two2n2n+24.)-dMult*(dMult-1.)*(dMult-2.)*(12.*three2n1n1n+16.*two1n1n+24.*three3n2n1n+48.*three2n1n1n+96.*two1n1n+8.*three4n3n1n+32.*three3n2n1n+96.*three2n1n1n+144.*two1n1n+6.*three4n2n2n+96.*two2n2n+36.*two2n2n+72.+48.*three3n2n1n+16.*two3n3n+72.*three2n1n1n+144.*two1n1n)-dMult*(dMult-1.)*(8.*two1n1n+12.*two2n2n+16.+8.*two3n3n+48.*two1n1n+1.*two4n4n+16.*two2n2n+18.)-dMult)/(dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)*(dMult-5.)*(dMult-6.)*(dMult-7.));
  
  fQCorrelations->Fill(30.,eight1n1n1n1n1n1n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)*(dMult-5.)*(dMult-6.)*(dMult-7.));
 
  f8pDistribution->Fill(eight1n1n1n1n1n1n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)*(dMult-5.)*(dMult-6.)*(dMult-7.));
 } 
 //---------------------------------------------------------------------------------------------------------
 
 
 
 
 //---------------------------------------------------------------------------------------------------------
 // weights:
 Bool_t useWeights = fUsePhiWeights||fUsePtWeights||fUseEtaWeights;

 TH1F *phiWeights = NULL; // histogram with phi weights
 TH1D *ptWeights  = NULL; // histogram with pt weights
 TH1D *etaWeights = NULL; // histogram with eta weights
   
 if(useWeights)
 {
  if(!fWeightsList)
  {
   cout<<" WARNING: fWeightsList is NULL pointer. "<<endl;
   exit(0);
  }
  if(fUsePhiWeights) 
  {
   phiWeights = dynamic_cast<TH1F *>(fWeightsList->FindObject("phi_weights"));
   if(!phiWeights)
   {
    cout<<" WARNING: couldn't access the histogram with phi weights. "<<endl;
    exit(0);
   } 
  } 
  if(fUsePtWeights) 
  { 
   ptWeights = dynamic_cast<TH1D *>(fWeightsList->FindObject("pt_weights"));
   if(!ptWeights) 
   {
    cout<<" WARNING: couldn't access the histogram with pt weights. "<<endl;
    exit(0);
   } 
  } 
  if(fUseEtaWeights) 
  {
   etaWeights = dynamic_cast<TH1D *>(fWeightsList->FindObject("eta_weights"));
   if(!etaWeights) 
   {
    cout<<" WARNING: couldn't access the histogram with eta weights. "<<endl;
    exit(0);
   }
  } 
 } 
  
 Int_t nBinsPhi = 0; 
 Double_t dBinWidthPt=0.;
 Double_t dBinWidthEta=0.;
  
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
      cout<<" WARNING: dBinWidthPtW != dBinWidthPt in AFAWQC::Make()"<<endl;
      exit(0);
     }
     Double_t dPtMinW = (ptWeights->GetXaxis())->GetXmin();
     if(dPtMinW != fPtMin)
     {
      cout<<" WARNING: dPtMinW != fPtMin in AFAWQC::Make()"<<endl;
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
      cout<<" WARNING: dBinWidthEtaW != dBinWidthEta in AFAWQC::Make()"<<endl;
      exit(0);
     }
     Double_t dEtaMinW = (etaWeights->GetXaxis())->GetXmin();
     if(dEtaMinW != fEtaMin)
     {
      cout<<" WARNING: dEtaMinW != fEtaMin in AFAWQC::Make()"<<endl;
      exit(0);
     }
    } 
   }          
  } // end of if(weightsList)
 
 Double_t wPhi = 1.; // phi weight
 Double_t wPt  = 1.; // pt weight
 Double_t wEta = 1.; // eta weight
 
 Double_t dPhi = 0.;
 Double_t dPt  = 0.;
 Double_t dEta = 0.;

 Double_t dQnkX[4][8] = {{0.}}; // sum_{i=1}^{M} w_i^k cos(n phi_i)
 Double_t dQnkY[4][8] = {{0.}}; // sum_{i=1}^{M} w_i^k sin(n phi_i)
 Double_t dSnk[4][8] = {{0.}};  //(sum_{i=1}^{M} w_i^k)^n
 
  for(Int_t i=0;i<nPrim;i++) // loop over all particles
  { 
   fTrack=anEvent->GetTrack(i);   
   if(fTrack)
   {
    if(fTrack->UseForIntegratedFlow()) 
    {
     dPhi = fTrack->Phi();
     dPt  = fTrack->Pt();
     dEta = fTrack->Eta();
     
     // determine Phi weight: 
     if(phiWeights && nBinsPhi)
     {
      wPhi = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(dPhi*nBinsPhi/TMath::TwoPi())));
     }
     // determine pt weight:    
     if(ptWeights && dBinWidthPt)
     {
      wPt = ptWeights->GetBinContent(1+(Int_t)(TMath::Floor((dPt-fPtMin)/dBinWidthPt))); 
     }            
     // determine eta weight:    
     if(etaWeights && dBinWidthEta)
     {
      wEta = etaWeights->GetBinContent(1+(Int_t)(TMath::Floor((dEta-fEtaMin)/dBinWidthEta))); 
     } 

     // (Q_{n,k})_x, (Q_{n,k})_y and S_{n,k}
     for(Int_t nn=0;nn<4;nn++)
     {
      for(Int_t k=0;k<8;k++)
      {
       dQnkX[nn][k]+=pow(wPhi*wPt*wEta,k+1)*TMath::Cos(2*(nn+1)*dPhi);
       dQnkY[nn][k]+=pow(wPhi*wPt*wEta,k+1)*TMath::Sin(2*(nn+1)*dPhi);
       dSnk[nn][k]+=pow(wPhi*wPt*wEta,k+1);
      } 
     }
     
    } // end of if (pTrack->UseForIntegratedFlow())
   } // end of if (pTrack)
   else {cerr << "no particle!!!"<<endl;}
  } // loop over particles
  
  for(Int_t nn=0;nn<4;nn++)
  {
   for(Int_t k=0;k<8;k++)
   {
    dSnk[nn][k]=pow(dSnk[nn][k],nn+1);
   }  
  }
  
 //.............................................................................................. 
 // Ms (introduced in order to simplify some Eqs. bellow)
 Double_t dM11 = dSnk[1][0]-dSnk[0][1]; // dM11 = sum_{i,j=1,i!=j}^M w_i w_j
 Double_t dM22 = dSnk[1][1]-dSnk[0][3]; // dM22 = sum_{i,j=1,i!=j}^M w_i^2 w_j^2
 Double_t dM33 = dSnk[1][2]-dSnk[0][5]; // dM33 = sum_{i,j=1,i!=j}^M w_i^3 w_j^3
 Double_t dM44 = dSnk[1][3]-dSnk[0][7]; // dM44 = sum_{i,j=1,i!=j}^M w_i^4 w_j^4
 Double_t dM31 = dSnk[0][2]*dSnk[0][0]-dSnk[0][3]; // dM31 = sum_{i,j=1,i!=j}^M w_i^3 w_j
 Double_t dM211 = dSnk[0][1]*dSnk[1][0]-2.*dSnk[0][2]*dSnk[0][0]-dSnk[1][1]+2.*dSnk[0][3]; // dM211 = sum_{i,j,k=1,i!=j!=k}^M w_i^2 w_j w_k
 Double_t dM1111 = dSnk[3][0]-6.*dM211-4.*dM31-3.*dM22-dSnk[0][3]; // dM1111 = sum_{i,j,k,l=1,i!=j!=k!=l}^M w_i w_j w_k w_l
 //..............................................................................................

 //---------------------------------------------------------------------------------------------------------
 //
 //                                        ***********************************************
 //                                        **** weighted multi-particle correlations: ****
 //                                        ***********************************************
 //
 // Remark 1: weighted multi-particle correlations calculated with Q-vectors are stored in fWeightedQCorrelations.
 // Remark 2: binning of fWeightedQCorrelations is organized as follows: 
 
 // binning
 //..............................................................................................
 // 1st bin: weighted <2>_{n|n}   = <w1   w2   cos( n*(phi1-phi2))>
 // 2nd bin: weighted <2>_{2n|2n} = <w1^2 w2^2 cos(2n*(phi1-phi2))>
 // 3rd bin: weighted <2>_{3n|3n} = <w1^3 w2^3 cos(3n*(phi1-phi2))>
 // 4th bin: weighted <2>_{4n|4n} = <w1^4 w2^4 cos(4n*(phi1-phi2))>
 // 5th bin: weighted <2>_{n|n} = <w1^3 w2 cos(n*(phi1-phi2))>
 // 6th bin: weighted <2>_{n|n} = <w1 w2 w3^2 cos(n*(phi1-phi2))>
 
 // 11th bin: weighted <3>_{2n|n,n} = <w1^2 w2 w3 cos(n*(2phi1-phi2-phi3))> 
 
 // 21st bin: weighted <4>_{n,n|n,n} = <w1 w2 w3 w4 cos(n*(phi1+phi2-phi3-phi4))>
 //..............................................................................................
 
 //.............................................................................................. 
 // weighted 2-particle correlations:
 Double_t two1n1nW1W1=0., two2n2nW2W2=0., two3n3nW3W3=0., two4n4nW4W4=0., two1n1nW3W1=0., two1n1nW2W1W1=0.;
 
 if(nRP>1) // nRP = number of particles used to determine the reaction plane
 { 
  if(dM11 != 0)
  {
   two1n1nW1W1 = (pow(dQnkX[0][0],2)+pow(dQnkY[0][0],2)-dSnk[0][1])/dM11; // <2>_{n|n}=<w1 w2 cos(n*(phi1-phi2))>
   fWeightedQCorrelations->Fill(0.,two1n1nW1W1,dM11);
  }
  if(dM22 != 0)
  {
   two2n2nW2W2 = (pow(dQnkX[1][1],2)+pow(dQnkY[1][1],2)-dSnk[0][3])/dM22; // <2>_{2n|2n}=<w1^2 w2^2 cos(2n*(phi1-phi2))>
   fWeightedQCorrelations->Fill(1.,two2n2nW2W2,dM22); 
  }
  if(dM33 != 0)
  {
   two3n3nW3W3 = (pow(dQnkX[2][2],2)+pow(dQnkY[2][2],2)-dSnk[0][5])/dM33; // <2>_{2n|2n}=<w1^3 w2^3 cos(3n*(phi1-phi2))>
   fWeightedQCorrelations->Fill(2.,two3n3nW3W3,dM33);
  }
  if(dM44 != 0)
  {
   two4n4nW4W4 = (pow(dQnkX[3][3],2)+pow(dQnkY[3][3],2)-dSnk[0][7])/dM44; // <2>_{2n|2n}=<w1^4 w2^4 cos(4n*(phi1-phi2))>
   fWeightedQCorrelations->Fill(3.,two4n4nW4W4,dM44);  
  } 
  if(dM31 != 0)
  {
   two1n1nW3W1 = (dQnkX[0][2]*dQnkX[0][0]+dQnkY[0][2]*dQnkY[0][0]-dSnk[0][3])/dM31; // <2>_{n|n}=<w1^3 w2 cos(n*(phi1-phi2))>
   fWeightedQCorrelations->Fill(4.,two1n1nW3W1,dM31);  
  } 
  if(dM211 != 0)
  {
   two1n1nW2W1W1 = (dSnk[0][1]*dM11*two1n1nW1W1-2.*dM31*two1n1nW3W1)/dM211; // <2>_{n|n}=<w1^2 w2 w3 cos(n*(phi1-phi2))>
   fWeightedQCorrelations->Fill(5.,two1n1nW2W1W1,dM211);  
  } 
 } // end of if(nRP>1)
 //..............................................................................................

 //..............................................................................................
 // weighted 3-particle correlations:
 Double_t three2n1n1nW2W1W1=0.;
 
 if(nRP>2) // nRP = number of particles used to determine the reaction plane
 { 
  if(dM211 != 0)
  {
   three2n1n1nW2W1W1 = (pow(dQnkX[0][0],2.)*dQnkX[1][1]+2.*dQnkX[0][0]*dQnkY[0][0]*dQnkY[1][1]-pow(dQnkY[0][0],2.)*dQnkX[1][1]-2.*dM31*two1n1nW3W1-dM22*two2n2nW2W2-dSnk[0][3])/dM211;
   fWeightedQCorrelations->Fill(10.,three2n1n1nW2W1W1,dM211);
  } 
 } // end of if(nRP>2) 
 //..............................................................................................
 
 //..............................................................................................
 // weighted 4-particle correlations:
 Double_t four1n1n1n1nW1W1W1W1=0.;
 
 if(nRP>3) // nRP = number of particles used to determine the reaction plane
 { 
  if(dM1111 != 0)
  {
   four1n1n1n1nW1W1W1W1 = (pow(pow(dQnkX[0][0],2.)+pow(dQnkY[0][0],2.),2)-2.*dM211*three2n1n1nW2W1W1-4.*dM211*two1n1nW2W1W1-4.*dM31*two1n1nW3W1-dM22*two2n2nW2W2-2.*dM22-dSnk[0][3])/dM1111;
   fWeightedQCorrelations->Fill(20.,four1n1n1n1nW1W1W1W1,dM1111);
  } 
 } // end of if(nRP>3) 
 //..............................................................................................
     
 
 
 
 //---------------------------------------------------------------------------------------------------------
 //
 //                                        *******************************************************
 //                                        **** weighted reduced multi-particle correlations: ****
 //                                        *******************************************************
 // 
 // pt POI
 TProfile *ptReq1nPrime = new TProfile("ptReq1nPrime","Re[q_{n}^{''}]",fnBinsPt,fPtMin,fPtMax,"s");
 TProfile *ptImq1nPrime = new TProfile("ptImq1nPrime","Im[q_{n}^{''}]",fnBinsPt,fPtMin,fPtMax,"s");
 TProfile *ptReq2nPrime = new TProfile("ptReq2nPrime","Re[q_{2n}^{''}]",fnBinsPt,fPtMin,fPtMax,"s");
 TProfile *ptImq2nPrime = new TProfile("ptImq2nPrime","Im[q_{2n}^{''}]",fnBinsPt,fPtMin,fPtMax,"s");
 
 TProfile *ptReq1nPrimePrime = new TProfile("ptReq1nPrimePrime","Re[q_{n}^{''}]",fnBinsPt,fPtMin,fPtMax,"s");
 TProfile *ptImq1nPrimePrime = new TProfile("ptImq1nPrimePrime","Im[q_{n}^{''}]",fnBinsPt,fPtMin,fPtMax,"s");
 TProfile *ptReq2nPrimePrime = new TProfile("ptReq2nPrimePrime","Re[q_{2n}^{''}]",fnBinsPt,fPtMin,fPtMax,"s");
 TProfile *ptImq2nPrimePrime = new TProfile("ptImq2nPrimePrime","Im[q_{2n}^{''}]",fnBinsPt,fPtMin,fPtMax,"s");
 
 TProfile *req1nW2PrimePrimePt = new TProfile("req1nW2PrimePrimePt","#sum_{i=1}^{m''} w_{i}^{2} cos(n(#phi_{i}))''",fnBinsPt,fPtMin,fPtMax,"s");
 TProfile *imq1nW2PrimePrimePt = new TProfile("imq1nW2PrimePrimePt","#sum_{i=1}^{m''} w_{i}^{2} sin(n(#phi_{i}))''",fnBinsPt,fPtMin,fPtMax,"s");
 TProfile *req2nW1PrimePrimePt = new TProfile("req2nW1PrimePrimePt","#sum_{i=1}^{m''} w_{i} cos(2n(#phi_{i}))''",fnBinsPt,fPtMin,fPtMax,"s");
 TProfile *imq2nW1PrimePrimePt = new TProfile("imq2nW1PrimePrimePt","#sum_{i=1}^{m''} w_{i} sin(2n(#phi_{i}))''",fnBinsPt,fPtMin,fPtMax,"s");
 
 TProfile *sumOfW1upTomPrimePrimePt = new TProfile("sumOfW1upTomPrimePrimePt","#sum_{i=1}^{m''} w_{i}''",fnBinsPt,fPtMin,fPtMax,"s");
 TProfile *sumOfW2upTomPrimePrimePt = new TProfile("sumOfW2upTomPrimePrimePt","#sum_{i=1}^{m''} w_{i}^{2}''",fnBinsPt,fPtMin,fPtMax,"s");
 TProfile *sumOfW3upTomPrimePrimePt = new TProfile("sumOfW3upTomPrimePrimePt","#sum_{i=1}^{m''} w_{i}^{3}''",fnBinsPt,fPtMin,fPtMax,"s");
 
 // eta POI 
 TProfile *etaReq1nPrime = new TProfile("etaReq1nPrime","Re[q_{n}^{''}]",fnBinsEta,fEtaMin,fEtaMax,"s");
 TProfile *etaImq1nPrime = new TProfile("etaImq1nPrime","Im[q_{n}^{''}]",fnBinsEta,fEtaMin,fEtaMax,"s");
 TProfile *etaReq2nPrime = new TProfile("etaReq2nPrime","Re[q_{2n}^{''}]",fnBinsEta,fEtaMin,fEtaMax,"s");
 TProfile *etaImq2nPrime = new TProfile("etaImq2nPrime","Im[q_{2n}^{''}]",fnBinsEta,fEtaMin,fEtaMax,"s");
 
 TProfile *etaReq1nPrimePrime = new TProfile("etaReq1nPrimePrime","Re[q_{n}^{''}]",fnBinsEta,fEtaMin,fEtaMax,"s");
 TProfile *etaImq1nPrimePrime = new TProfile("etaImq1nPrimePrime","Im[q_{n}^{''}]",fnBinsEta,fEtaMin,fEtaMax,"s");
 TProfile *etaReq2nPrimePrime = new TProfile("etaReq2nPrimePrime","Re[q_{2n}^{''}]",fnBinsEta,fEtaMin,fEtaMax,"s");
 TProfile *etaImq2nPrimePrime = new TProfile("etaImq2nPrimePrime","Im[q_{2n}^{''}]",fnBinsEta,fEtaMin,fEtaMax,"s");
 
 TProfile *req1nW2PrimePrimeEta = new TProfile("req1nW2PrimePrimeEta","#sum_{i=1}^{m''} w_{i}^{2} cos(n(#phi_{i}))''",fnBinsEta,fEtaMin,fEtaMax,"s");
 TProfile *imq1nW2PrimePrimeEta = new TProfile("imq1nW2PrimePrimeEta","#sum_{i=1}^{m''} w_{i}^{2} sin(n(#phi_{i}))''",fnBinsEta,fEtaMin,fEtaMax,"s");
 TProfile *req2nW1PrimePrimeEta = new TProfile("req2nW1PrimePrimeEta","#sum_{i=1}^{m''} w_{i} cos(2n(#phi_{i}))''",fnBinsEta,fEtaMin,fEtaMax,"s");
 TProfile *imq2nW1PrimePrimeEta = new TProfile("imq2nW1PrimePrimeEta","#sum_{i=1}^{m''} w_{i} sin(2n(#phi_{i}))''",fnBinsEta,fEtaMin,fEtaMax,"s");
 
 TProfile *sumOfW1upTomPrimePrimeEta = new TProfile("sumOfW1upTomPrimePrimeEta","#sum_{i=1}^{m''} w_{i}''",fnBinsEta,fEtaMin,fEtaMax,"s");
 TProfile *sumOfW2upTomPrimePrimeEta = new TProfile("sumOfW2upTomPrimePrimeEta","#sum_{i=1}^{m''} w_{i}^{2}''",fnBinsEta,fEtaMin,fEtaMax,"s");
 TProfile *sumOfW3upTomPrimePrimeEta = new TProfile("sumOfW3upTomPrimePrimeEta","#sum_{i=1}^{m''} w_{i}^{3}''",fnBinsEta,fEtaMin,fEtaMax,"s");
 
 // pt RP
 TProfile *ptReq1n = new TProfile("ptReq1n","Re[q_{n}]",fnBinsPt,fPtMin,fPtMax,"s");
 TProfile *ptImq1n = new TProfile("ptImq1n","Im[q_{n}]",fnBinsPt,fPtMin,fPtMax,"s");
 TProfile *ptReq2n = new TProfile("ptReq2n","Re[q_{2n}]",fnBinsPt,fPtMin,fPtMax,"s");
 TProfile *ptImq2n = new TProfile("ptImq2n","Im[q_{2n}]",fnBinsPt,fPtMin,fPtMax,"s");
 
 TProfile *req1nW2Pt = new TProfile("req1nW2Pt","#sum_{i=1}^{m} w_{i}^{2} cos(n(#phi_{i}))''",fnBinsPt,fPtMin,fPtMax,"s");
 TProfile *imq1nW2Pt = new TProfile("imq1nW2Pt","#sum_{i=1}^{m} w_{i}^{2} sin(n(#phi_{i}))''",fnBinsPt,fPtMin,fPtMax,"s");
 TProfile *req2nW1Pt = new TProfile("req2nW1Pt","#sum_{i=1}^{m} w_{i} cos(2n(#phi_{i}))''",fnBinsPt,fPtMin,fPtMax,"s");
 TProfile *imq2nW1Pt = new TProfile("imq2nW1Pt","#sum_{i=1}^{m} w_{i} sin(2n(#phi_{i}))''",fnBinsPt,fPtMin,fPtMax,"s");
 
 TProfile *sumOfW1upTomPt = new TProfile("sumOfW1upTomPt","#sum_{i=1}^{m} w_{i}''",fnBinsPt,fPtMin,fPtMax,"s");
 TProfile *sumOfW2upTomPt = new TProfile("sumOfW2upTomPt","#sum_{i=1}^{m} w_{i}^{2}''",fnBinsPt,fPtMin,fPtMax,"s");
 TProfile *sumOfW3upTomPt = new TProfile("sumOfW3upTomPt","#sum_{i=1}^{m} w_{i}^{3}''",fnBinsPt,fPtMin,fPtMax,"s");
 
 // eta RP
 TProfile *etaReq1n = new TProfile("etaReq1n","Re[q_{n}]",fnBinsEta,fEtaMin,fEtaMax,"s");
 TProfile *etaImq1n = new TProfile("etaImq1n","Im[q_{n}]",fnBinsEta,fEtaMin,fEtaMax,"s");
 TProfile *etaReq2n = new TProfile("etaReq2n","Re[q_{2n}]",fnBinsEta,fEtaMin,fEtaMax,"s");
 TProfile *etaImq2n = new TProfile("etaImq2n","Im[q_{2n}]",fnBinsEta,fEtaMin,fEtaMax,"s");
 
 TProfile *req1nW2Eta = new TProfile("req1nW2Eta","#sum_{i=1}^{m} w_{i}^{2} cos(n(#phi_{i}))''",fnBinsEta,fEtaMin,fEtaMax,"s");
 TProfile *imq1nW2Eta = new TProfile("imq1nW2Eta","#sum_{i=1}^{m} w_{i}^{2} sin(n(#phi_{i}))''",fnBinsEta,fEtaMin,fEtaMax,"s");
 TProfile *req2nW1Eta = new TProfile("req2nW1Eta","#sum_{i=1}^{m} w_{i} cos(2n(#phi_{i}))''",fnBinsEta,fEtaMin,fEtaMax,"s");
 TProfile *imq2nW1Eta = new TProfile("imq2nW1Eta","#sum_{i=1}^{m} w_{i} sin(2n(#phi_{i}))''",fnBinsEta,fEtaMin,fEtaMax,"s");
 
 TProfile *sumOfW1upTomEta = new TProfile("sumOfW1upTomEta","#sum_{i=1}^{m} w_{i}''",fnBinsEta,fEtaMin,fEtaMax,"s");
 TProfile *sumOfW2upTomEta = new TProfile("sumOfW2upTomEta","#sum_{i=1}^{m} w_{i}^{2}''",fnBinsEta,fEtaMin,fEtaMax,"s");
 TProfile *sumOfW3upTomEta = new TProfile("sumOfW3upTomEta","#sum_{i=1}^{m} w_{i}^{3}''",fnBinsEta,fEtaMin,fEtaMax,"s");

 for(Int_t i=0;i<nPrim;i++) // loop over all particles  
 { 
  fTrack=anEvent->GetTrack(i);
  if(fTrack)
  {
   if(fTrack->UseForDifferentialFlow()) // checking if particle is POI 
   {
    if(fTrack->UseForIntegratedFlow()) // checking if particle is both POI and RP 
    {
     // get azimuthal angle, momentum and pseudorapidity of a particle:
     dPhi = fTrack->Phi();
     dPt  = fTrack->Pt();
     dEta = fTrack->Eta();
     // phi weights:
     if(fUsePhiWeights) 
     {
      nBinsPhi = phiWeights->GetNbinsX();
      if(nBinsPhi) 
      {
       wPhi = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(dPhi*nBinsPhi/TMath::TwoPi())));
      }
     } 
     // pt weights:
     if(fUsePtWeights)
     {          
      if(dBinWidthPt) 
      {
       wPt = ptWeights->GetBinContent(1+(Int_t)(TMath::Floor((dPt-fPtMin)/dBinWidthPt)));
      }
     }             
     // eta weights:
     if(fUseEtaWeights)
     {    
      if(dBinWidthEta)
      {
       wEta = etaWeights->GetBinContent(1+(Int_t)(TMath::Floor((dEta-fEtaMin)/dBinWidthEta))); 
      }
     }
     // pt:
     ptReq1nPrimePrime->Fill(dPt,cos(n*dPhi),1.); 
     ptImq1nPrimePrime->Fill(dPt,sin(n*dPhi),1.);
     ptReq2nPrimePrime->Fill(dPt,cos(2.*n*dPhi),1.);
     ptImq2nPrimePrime->Fill(dPt,sin(2.*n*dPhi),1.);
     // weighted pt:
     req1nW2PrimePrimePt->Fill(dPt,cos(n*dPhi),pow(wPhi*wPt*wEta,2.));
     imq1nW2PrimePrimePt->Fill(dPt,sin(n*dPhi),pow(wPhi*wPt*wEta,2.));
     req2nW1PrimePrimePt->Fill(dPt,cos(2.*n*dPhi),wPhi*wPt*wEta);
     imq2nW1PrimePrimePt->Fill(dPt,sin(2.*n*dPhi),wPhi*wPt*wEta);
     sumOfW1upTomPrimePrimePt->Fill(dPt,wPhi*wPt*wEta,1.);
     sumOfW2upTomPrimePrimePt->Fill(dPt,pow(wPhi*wPt*wEta,2),1.);
     sumOfW3upTomPrimePrimePt->Fill(dPt,pow(wPhi*wPt*wEta,3),1.);
     
     // eta:
     etaReq1nPrimePrime->Fill(dEta,cos(n*dPhi),1.); 
     etaImq1nPrimePrime->Fill(dEta,sin(n*dPhi),1.);
     etaReq2nPrimePrime->Fill(dEta,cos(2.*n*dPhi),1.);
     etaImq2nPrimePrime->Fill(dEta,sin(2.*n*dPhi),1.);
     // weighted eta:
     req1nW2PrimePrimeEta->Fill(dEta,cos(n*dPhi),pow(wPhi*wPt*wEta,2.));
     imq1nW2PrimePrimeEta->Fill(dEta,sin(n*dPhi),pow(wPhi*wPt*wEta,2.));
     req2nW1PrimePrimeEta->Fill(dEta,cos(2.*n*dPhi),wPhi*wPt*wEta);
     imq2nW1PrimePrimeEta->Fill(dEta,sin(2.*n*dPhi),wPhi*wPt*wEta);
     sumOfW1upTomPrimePrimeEta->Fill(dEta,wPhi*wPt*wEta,1.);
     sumOfW2upTomPrimePrimeEta->Fill(dEta,pow(wPhi*wPt*wEta,2),1.);
     sumOfW3upTomPrimePrimeEta->Fill(dEta,pow(wPhi*wPt*wEta,3),1.);

    }else if(!(fTrack->UseForIntegratedFlow())) // checking if particles is POI and not RP  
     {
      // get azimuthal angle, momentum and pseudorapidity of a particle:
      dPhi = fTrack->Phi();
      dPt  = fTrack->Pt();
      dEta = fTrack->Eta();
      // pt:
      ptReq1nPrime->Fill(dPt,cos(n*dPhi),1.);   
      ptImq1nPrime->Fill(dPt,sin(n*dPhi),1.);
      ptReq2nPrime->Fill(dPt,cos(2.*n*dPhi),1.);
      ptImq2nPrime->Fill(dPt,sin(2.*n*dPhi),1.);

      // eta:
      etaReq1nPrime->Fill(dEta,cos(n*dPhi),1.); 
      etaImq1nPrime->Fill(dEta,sin(n*dPhi),1.);
      etaReq2nPrime->Fill(dEta,cos(2.*n*dPhi),1.);
      etaImq2nPrime->Fill(dEta,sin(2.*n*dPhi),1.);

     } // end of else if(!(fTrack->UseForIntegratedFlow())) // checking if particles is POI and not RP 
   } // end of if(fTrack->UseForDifferentialFlow()) // checking if particle is POI 
     
   if(fTrack->UseForIntegratedFlow()) // checking if particles is only RP:
   {
    dPhi = fTrack->Phi();
    dPt  = fTrack->Pt();
    dEta = fTrack->Eta();
     
    // phi weights:
    if(fUsePhiWeights) 
    {
     nBinsPhi = phiWeights->GetNbinsX();
     if(nBinsPhi) 
     {
      wPhi = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(dPhi*nBinsPhi/TMath::TwoPi())));
     }
    } 
    // pt weights:
    if(fUsePtWeights)
    {          
     if(dBinWidthPt) 
     {
      wPt = ptWeights->GetBinContent(1+(Int_t)(TMath::Floor((dPt-fPtMin)/dBinWidthPt)));
     }
    }             
    // eta weights:
    if(fUseEtaWeights)
    {    
     if(dBinWidthEta)
     {
      wEta = etaWeights->GetBinContent(1+(Int_t)(TMath::Floor((dEta-fEtaMin)/dBinWidthEta))); 
     }
    }
    // pt:
    ptReq1n->Fill(dPt,cos(n*dPhi),1.); 
    ptImq1n->Fill(dPt,sin(n*dPhi),1.);
    ptReq2n->Fill(dPt,cos(2.*n*dPhi),1.);
    ptImq2n->Fill(dPt,sin(2.*n*dPhi),1.);
    // weighted pt:
    req1nW2Pt->Fill(dPt,cos(n*dPhi),pow(wPhi*wPt*wEta,2.));
    imq1nW2Pt->Fill(dPt,sin(n*dPhi),pow(wPhi*wPt*wEta,2.));
    req2nW1Pt->Fill(dPt,cos(2.*n*dPhi),wPhi*wPt*wEta);
    imq2nW1Pt->Fill(dPt,sin(2.*n*dPhi),wPhi*wPt*wEta);
    sumOfW1upTomPt->Fill(dPt,wPhi*wPt*wEta,1.);
    sumOfW2upTomPt->Fill(dPt,pow(wPhi*wPt*wEta,2),1.);
    sumOfW3upTomPt->Fill(dPt,pow(wPhi*wPt*wEta,3),1.);
     
    // eta:
    etaReq1n->Fill(dEta,cos(n*dPhi),1.); 
    etaImq1n->Fill(dEta,sin(n*dPhi),1.);
    etaReq2n->Fill(dEta,cos(2.*n*dPhi),1.);
    etaImq2n->Fill(dEta,sin(2.*n*dPhi),1.);
    // weighted eta:
    req1nW2Eta->Fill(dEta,cos(n*dPhi),pow(wPhi*wPt*wEta,2.));
    imq1nW2Eta->Fill(dEta,sin(n*dPhi),pow(wPhi*wPt*wEta,2.));
    req2nW1Eta->Fill(dEta,cos(2.*n*dPhi),wPhi*wPt*wEta);
    imq2nW1Eta->Fill(dEta,sin(2.*n*dPhi),wPhi*wPt*wEta);
    sumOfW1upTomEta->Fill(dEta,wPhi*wPt*wEta,1.);
    sumOfW2upTomEta->Fill(dEta,pow(wPhi*wPt*wEta,2),1.);
    sumOfW3upTomEta->Fill(dEta,pow(wPhi*wPt*wEta,3),1.);
   } // end of if(fTrack->UseForIntegratedFlow()) // checking if particles is only RP:
  } // end of if(fTrack}      
 } // end of for(Int_t i=0;i<nPrim;i++) 
 
 //...........................................................................................................
 // PrimePrime Pt POI
 Double_t qxPrimePrimePtPOI=0.,qyPrimePrimePtPOI=0.,q2xPrimePrimePtPOIHere=0.,q2yPrimePrimePtPOIHere=0.;//add comments for these variable
 Double_t qxW2PrimePrimePtPOI=0.,qyW2PrimePrimePtPOI=0.,q2xW1PrimePrimePtPOI=0.,q2yW1PrimePrimePtPOI=0.;//add comments for these variable
 Double_t dS11mPrimePrimePtPOI=0.; // to be improved (name)
 Double_t dS12mPrimePrimePtPOI=0.; // to be improved (name)
 Double_t dS13mPrimePrimePtPOI=0.; // to be improved (name)
 Double_t mPrimePrimePtPOI=0.; // to be improved (name)

 Double_t dM1pp11PtPOI=0.; // to be improved (name)
 Double_t dM0pp111PtPOI=0.; // to be improved (name)
 Double_t dM0pp12PtPOI=0.;
 Double_t dM2pp1PtPOI=0.;
 Double_t dM1pp2PtPOI=0.;
 Double_t dM0pp3PtPOI=0.;
 
 // Prime Pt POI
 Double_t qxPrimePtPOI=0.,qyPrimePtPOI=0.;
 Double_t mPrimePtPOI=0.;
 Double_t dM0p111PtPOI=0.; // to be improved (name)
  
 for(Int_t bin=1;bin<(fnBinsPt+1);bin++) // loop over pt-bins 
 {     
  // q'':                    
  qxPrimePrimePtPOI = (ptReq1nPrimePrime->GetBinContent(bin))*(ptReq1nPrimePrime->GetBinEntries(bin));
  qyPrimePrimePtPOI = (ptImq1nPrimePrime->GetBinContent(bin))*(ptImq1nPrimePrime->GetBinEntries(bin)); 
  q2xPrimePrimePtPOIHere = (ptReq2nPrimePrime->GetBinContent(bin))*(ptReq2nPrimePrime->GetBinEntries(bin));  
  q2yPrimePrimePtPOIHere = (ptImq2nPrimePrime->GetBinContent(bin))*(ptImq2nPrimePrime->GetBinEntries(bin)); 
  
  qxW2PrimePrimePtPOI = (req1nW2PrimePrimePt->GetBinContent(bin))*(req1nW2PrimePrimePt->GetBinEntries(bin));
  qyW2PrimePrimePtPOI = (imq1nW2PrimePrimePt->GetBinContent(bin))*(imq1nW2PrimePrimePt->GetBinEntries(bin));
  
  q2xW1PrimePrimePtPOI = (req2nW1PrimePrimePt->GetBinContent(bin))*(req2nW1PrimePrimePt->GetBinEntries(bin));
  q2yW1PrimePrimePtPOI = (imq2nW1PrimePrimePt->GetBinContent(bin))*(imq2nW1PrimePrimePt->GetBinEntries(bin));
  
  dS11mPrimePrimePtPOI = (sumOfW1upTomPrimePrimePt->GetBinContent(bin))*(sumOfW1upTomPrimePrimePt->GetBinEntries(bin));
  dS12mPrimePrimePtPOI = (sumOfW2upTomPrimePrimePt->GetBinContent(bin))*(sumOfW2upTomPrimePrimePt->GetBinEntries(bin));
  dS13mPrimePrimePtPOI = (sumOfW3upTomPrimePrimePt->GetBinContent(bin))*(sumOfW3upTomPrimePrimePt->GetBinEntries(bin));

  mPrimePrimePtPOI = sumOfW1upTomPrimePrimePt->GetBinEntries(bin); // to be improved
      
  dM1pp11PtPOI=dS11mPrimePrimePtPOI*(dSnk[1][0]-dSnk[0][1])-2.*dS12mPrimePrimePtPOI*dSnk[0][0]+2.*dS13mPrimePrimePtPOI;
  dM1pp2PtPOI=dS11mPrimePrimePtPOI*dSnk[0][1]-dS13mPrimePrimePtPOI;
  dM2pp1PtPOI=dS12mPrimePrimePtPOI*dSnk[0][0]-dS13mPrimePrimePtPOI;
  dM0pp3PtPOI=mPrimePrimePtPOI*dSnk[0][2]-dS13mPrimePrimePtPOI;
  dM0pp12PtPOI=mPrimePrimePtPOI*dSnk[0][0]*dSnk[0][1]-dM2pp1PtPOI-dM1pp2PtPOI-dM0pp3PtPOI-dS13mPrimePrimePtPOI;
  dM0pp111PtPOI=mPrimePrimePtPOI*dSnk[2][0]-3.*dM1pp11PtPOI-3.*dM0pp12PtPOI-3.*dM2pp1PtPOI-3.*dM1pp2PtPOI-dM0pp3PtPOI-dS13mPrimePrimePtPOI;
  
  // q':
  qxPrimePtPOI = (ptReq1nPrime->GetBinContent(bin))*(ptReq1nPrime->GetBinEntries(bin));
  qyPrimePtPOI = (ptImq1nPrime->GetBinContent(bin))*(ptImq1nPrime->GetBinEntries(bin));
  
  mPrimePtPOI = ptReq1nPrime->GetBinEntries(bin); // to be improved
  dM0p111PtPOI=mPrimePtPOI*(dSnk[2][0]-3.*dSnk[0][1]*dSnk[0][0]+2.*dSnk[0][2]);
 
  // 2-p the needed one
  Double_t two1n1nWPerPtBinPOI=0.; 
  if((mPrimePrimePtPOI+mPrimePtPOI)*dSnk[0][0]-dS11mPrimePrimePtPOI>0)
  {
   two1n1nWPerPtBinPOI = (qxPrimePrimePtPOI*dQnkX[0][0]+qyPrimePrimePtPOI*dQnkY[0][0]+qxPrimePtPOI*dQnkX[0][0]+qyPrimePtPOI*dQnkY[0][0]-dS11mPrimePrimePtPOI)/((mPrimePrimePtPOI+mPrimePtPOI)*dSnk[0][0]-dS11mPrimePrimePtPOI); 
  
   f2WPerPtBin1n1nPOI->Fill(fPtMin+(bin-1)*dBinWidthPt,two1n1nWPerPtBinPOI,(mPrimePrimePtPOI+mPrimePtPOI)*dSnk[0][0]-dS11mPrimePrimePtPOI);     
  }
 
  // 2-p temporary one
  Double_t two1n1nW1ppW1W1PtPOI=0.; // <w1 w2 w3 cos(n(phi2-phi3))> // OK!!!
  if(dM1pp11PtPOI)
  {
   two1n1nW1ppW1W1PtPOI = ((pow(dQnkX[0][0],2.)+pow(dQnkY[0][0],2.))*dS11mPrimePrimePtPOI-2.*(qxW2PrimePrimePtPOI*dQnkX[0][0]+qyW2PrimePrimePtPOI*dQnkY[0][0]) - dS11mPrimePrimePtPOI*dSnk[0][1]+2.*dS13mPrimePrimePtPOI)/dM1pp11PtPOI; // CORRECT !!! <w1 w2 w3 cos(n(phi2-phi3))>
  }
  
  // 2-p temporary one
  Double_t two1npp1nW1W2PtPOI=0.; // <w2 w3^2 cos(n(psi1-phi2))> // OK !!!
  if(dM0pp12PtPOI)
  {
   two1npp1nW1W2PtPOI = (dSnk[0][1]*(qxPrimePrimePtPOI*dQnkX[0][0]+qyPrimePrimePtPOI*dQnkY[0][0])-(qxW2PrimePrimePtPOI*dQnkX[0][0]+qyW2PrimePrimePtPOI*dQnkY[0][0])-dM1pp2PtPOI-(qxPrimePrimePtPOI*dQnkX[0][2]+qyPrimePrimePtPOI*dQnkY[0][2])+dS13mPrimePrimePtPOI)/dM0pp12PtPOI; // CORRECT !!! <w2 w3^2 cos(n(psi1-phi2))>
  }
  
  // 2-p temporary one
  Double_t two1npp1nW2ppW1PtPOI=0.; // <w1^2 w2 cos(n(psi1-phi2))> // OK !!!
  if(dM2pp1PtPOI)
  {
   two1npp1nW2ppW1PtPOI = ((qxW2PrimePrimePtPOI*dQnkX[0][0]+qyW2PrimePrimePtPOI*dQnkY[0][0])-dS13mPrimePrimePtPOI)/dM2pp1PtPOI; // CORRECT !!! <w1^2 w2 cos(n(psi1-phi2))> 
  }
  
  // 2-p temporary one
  Double_t two2npp2nW1ppW2PtPOI=0.; // <w1 w2^2 cos(2n(psi1-phi2))> OK !!!
  if(dM1pp2PtPOI)
  {
   two2npp2nW1ppW2PtPOI = ((q2xW1PrimePrimePtPOI*dQnkX[1][1]+q2yW1PrimePrimePtPOI*dQnkY[1][1])-dS13mPrimePrimePtPOI)/dM1pp2PtPOI; // CORRECT !!! <w1 w2^2 cos(2n(psi1-phi2))> 
  }
  
  // 2-p temporary one
  Double_t two1npp1nW3PtPOI=0.; // <w2^3 cos(n(psi1-phi2))> // OK !!!
  if(dM0pp3PtPOI)
  {
   two1npp1nW3PtPOI = (qxPrimePrimePtPOI*dQnkX[0][2]+qyPrimePrimePtPOI*dQnkY[0][2]-dS13mPrimePrimePtPOI)/dM0pp3PtPOI; // CORRECT !!! <w2^3 cos(n(psi1-phi2))>
  }
   
  // 3-p temporary one
  Double_t three2npp1n1nW1ppW1W1PtPOI=0.; // <w1 w2 w3 cos(n(2psi1-phi2-phi3))> // OK!!!
  if(dM1pp11PtPOI)
  {
   three2npp1n1nW1ppW1W1PtPOI = (q2xW1PrimePrimePtPOI*(dQnkX[0][0]*dQnkX[0][0]-dQnkY[0][0]*dQnkY[0][0])+2.*q2yW1PrimePrimePtPOI*dQnkX[0][0]*dQnkY[0][0]-2.*(qxW2PrimePrimePtPOI*dQnkX[0][0]+qyW2PrimePrimePtPOI*dQnkY[0][0])-(q2xW1PrimePrimePtPOI*dQnkX[1][1]+q2yW1PrimePrimePtPOI*dQnkY[1][1])+2.*dS13mPrimePrimePtPOI)/dM1pp11PtPOI; // CORRECT !!! <w1 w2 w3 cos(n(2psi1-phi2-phi3))> 
  }
  
  // 3-p temporary one
  Double_t three1npp1n2nW0ppW1W2PtPOI=0.; // <w2 w3^2 cos(n(psi1+phi2-2*phi3))> // OK!!!
  if(dM0pp12PtPOI)
  {
   three1npp1n2nW0ppW1W2PtPOI = (qxPrimePrimePtPOI*(dQnkX[0][0]*dQnkX[1][1]+dQnkY[0][0]*dQnkY[1][1])-qyPrimePrimePtPOI*(dQnkY[0][0]*dQnkX[1][1]-dQnkX[0][0]*dQnkY[1][1])-(qxW2PrimePrimePtPOI*dQnkX[0][0]+qyW2PrimePrimePtPOI*dQnkY[0][0])-(q2xW1PrimePrimePtPOI*dQnkX[1][1]+q2yW1PrimePrimePtPOI*dQnkY[1][1])-(qxPrimePrimePtPOI*dQnkX[0][2]+qyPrimePrimePtPOI*dQnkY[0][2])+2.*dS13mPrimePrimePtPOI)/dM0pp12PtPOI; // CORRECT !!! <w2 w3^2 cos(n(psi1+phi2-2.*phi3))>
  }
 
  /*
  // 4-p RP part
  Double_t four1npp1n1n1nW1W1W1=0.; // <w1 w2 w3 cos(n(psi1+phi1-phi2-phi3))>
  if(dM0pp111PtPOI)
  {   
   four1npp1n1n1nW1W1W1 = ((pow(dQnkX[0][0],2.)+pow(dQnkY[0][0],2.))*(qxPrimePrimePtPOI*dQnkX[0][0]+qyPrimePrimePtPOI*dQnkY[0][0])-2.*dM1pp11PtPOI*two1n1nW1ppW1W1PtPOI-dM1pp11PtPOI*three2npp1n1nW1ppW1W1PtPOI-dM0pp12PtPOI*three1npp1n2nW0ppW1W2PtPOI-2.*dM0pp12PtPOI*two1npp1nW1W2PtPOI-3.*dM2pp1PtPOI*two1npp1nW2ppW1PtPOI-2.*dM1pp2PtPOI-dM1pp2PtPOI*two2npp2nW1ppW2PtPOI-dM0pp3PtPOI*two1npp1nW3PtPOI-dS13mPrimePrimePtPOI)/(dM0pp111PtPOI); 
  }
  */
  
  /*
  // 4-p POI part
  Double_t four1npp1n1n1nW1W1W1POI=0.;
  if(dM0p111PtPOI>0&&mPrimePtPOI>0&&nRP>0)
  {
   four1npp1n1n1nW1W1W1POI = ((pow(dQnkX[0][0],2.)+pow(dQnkY[0][0],2.))*(qxPrimePtPOI*dQnkX[0][0]+qyPrimePtPOI*dQnkY[0][0])-2.*dSnk[0][1]* (qxPrimePtPOI*dQnkX[0][0]+qyPrimePtPOI*dQnkY[0][0])+2.*(qxPrimePtPOI*dQnkX[0][2]+qyPrimePtPOI*dQnkY[0][2])-qxPrimePtPOI*(dQnkX[0][0]*dQnkX[1][1]+dQnkY[0][0]*dQnkY[1][1])+qyPrimePtPOI*(dQnkY[0][0]*dQnkX[1][1]-dQnkX[0][0]*dQnkY[1][1]))/dM0p111PtPOI;
  }
  */
 
  // 4-p RP and POI in all combinations (full, partial and no overlap)
  Double_t four1npp1n1n1nW1W1W1PtPOI=0.;
  if(dM0pp111PtPOI+dM0p111PtPOI)
  {
  four1npp1n1n1nW1W1W1PtPOI = ((pow(dQnkX[0][0],2.)+pow(dQnkY[0][0],2.))*(qxPrimePrimePtPOI*dQnkX[0][0]+qyPrimePrimePtPOI*dQnkY[0][0])-2.*dM1pp11PtPOI*two1n1nW1ppW1W1PtPOI-dM1pp11PtPOI*three2npp1n1nW1ppW1W1PtPOI-dM0pp12PtPOI*three1npp1n2nW0ppW1W2PtPOI-2.*dM0pp12PtPOI*two1npp1nW1W2PtPOI-3.*dM2pp1PtPOI*two1npp1nW2ppW1PtPOI-2.*dM1pp2PtPOI-dM1pp2PtPOI*two2npp2nW1ppW2PtPOI-dM0pp3PtPOI*two1npp1nW3PtPOI-dS13mPrimePrimePtPOI+(pow(dQnkX[0][0],2.)+pow(dQnkY[0][0],2.))*(qxPrimePtPOI*dQnkX[0][0]+qyPrimePtPOI*dQnkY[0][0])-2.*dSnk[0][1]* (qxPrimePtPOI*dQnkX[0][0]+qyPrimePtPOI*dQnkY[0][0])+2.*(qxPrimePtPOI*dQnkX[0][2]+qyPrimePtPOI*dQnkY[0][2])-qxPrimePtPOI*(dQnkX[0][0]*dQnkX[1][1]+dQnkY[0][0]*dQnkY[1][1])+qyPrimePtPOI*(dQnkY[0][0]*dQnkX[1][1]-dQnkX[0][0]*dQnkY[1][1]))/(dM0pp111PtPOI+dM0p111PtPOI);
 
   f4WPerPtBin1n1n1n1nPOI->Fill(fPtMin+(bin-1)*dBinWidthPt,four1npp1n1n1nW1W1W1PtPOI,dM0pp111PtPOI+dM0p111PtPOI);
  } // end of if(dM0pp111PtPOI+dM0p111PtPOI)
 } // for(Int_t bin=1;bin<(fnBinsPt+1);bin++) // loop over pt-bins
 //...........................................................................................................

 //...........................................................................................................  
 // PrimePrime Eta POI
 Double_t qxPrimePrimeEtaPOI=0.,qyPrimePrimeEtaPOI=0.,q2xPrimePrimeEtaPOIHere=0.,q2yPrimePrimeEtaPOIHere=0.;//add comments for these variable
 Double_t qxW2PrimePrimeEtaPOI=0.,qyW2PrimePrimeEtaPOI=0.,q2xW1PrimePrimeEtaPOI=0.,q2yW1PrimePrimeEtaPOI=0.;//add comments for these variable
 Double_t dS11mPrimePrimeEtaPOI=0.; // to be improved (name)
 Double_t dS12mPrimePrimeEtaPOI=0.; // to be improved (name)
 Double_t dS13mPrimePrimeEtaPOI=0.; // to be improved (name)
 Double_t mPrimePrimeEtaPOIHere=0.; // to be improved (name)

 Double_t dM1pp11EtaPOI=0.; // to be improved (name)
 Double_t dM0pp111EtaPOI=0.; // to be improved (name)
 Double_t dM0pp12EtaPOI=0.;
 Double_t dM2pp1EtaPOI=0.;
 Double_t dM1pp2EtaPOI=0.;
 Double_t dM0pp3EtaPOI=0.;
 
 // Prime Eta POI
 Double_t qxPrimeEtaPOI=0.,qyPrimeEtaPOI=0.;
 Double_t mPrimeEtaPOI=0.;
 Double_t dM0p111EtaPOI=0.; // to be improved (name)
 
 for(Int_t bin=1;bin<(fnBinsEta+1);bin++) // loop over eta-bins 
 {     
  // q'':                    
  qxPrimePrimeEtaPOI = (etaReq1nPrimePrime->GetBinContent(bin))*(etaReq1nPrimePrime->GetBinEntries(bin));
  qyPrimePrimeEtaPOI = (etaImq1nPrimePrime->GetBinContent(bin))*(etaImq1nPrimePrime->GetBinEntries(bin)); 
  q2xPrimePrimeEtaPOIHere = (etaReq2nPrimePrime->GetBinContent(bin))*(etaReq2nPrimePrime->GetBinEntries(bin));  
  q2yPrimePrimeEtaPOIHere = (etaImq2nPrimePrime->GetBinContent(bin))*(etaImq2nPrimePrime->GetBinEntries(bin)); 
  
  qxW2PrimePrimeEtaPOI = (req1nW2PrimePrimeEta->GetBinContent(bin))*(req1nW2PrimePrimeEta->GetBinEntries(bin));
  qyW2PrimePrimeEtaPOI = (imq1nW2PrimePrimeEta->GetBinContent(bin))*(imq1nW2PrimePrimeEta->GetBinEntries(bin));
  
  q2xW1PrimePrimeEtaPOI = (req2nW1PrimePrimeEta->GetBinContent(bin))*(req2nW1PrimePrimeEta->GetBinEntries(bin));
  q2yW1PrimePrimeEtaPOI = (imq2nW1PrimePrimeEta->GetBinContent(bin))*(imq2nW1PrimePrimeEta->GetBinEntries(bin));
  
  dS11mPrimePrimeEtaPOI = (sumOfW1upTomPrimePrimeEta->GetBinContent(bin))*(sumOfW1upTomPrimePrimeEta->GetBinEntries(bin));
  dS12mPrimePrimeEtaPOI = (sumOfW2upTomPrimePrimeEta->GetBinContent(bin))*(sumOfW2upTomPrimePrimeEta->GetBinEntries(bin));
  dS13mPrimePrimeEtaPOI = (sumOfW3upTomPrimePrimeEta->GetBinContent(bin))*(sumOfW3upTomPrimePrimeEta->GetBinEntries(bin));

  mPrimePrimeEtaPOIHere = sumOfW1upTomPrimePrimeEta->GetBinEntries(bin); // to be improved
      
  dM1pp11EtaPOI=dS11mPrimePrimeEtaPOI*(dSnk[1][0]-dSnk[0][1])-2.*dS12mPrimePrimeEtaPOI*dSnk[0][0]+2.*dS13mPrimePrimeEtaPOI;
  dM1pp2EtaPOI=dS11mPrimePrimeEtaPOI*dSnk[0][1]-dS13mPrimePrimeEtaPOI;
  dM2pp1EtaPOI=dS12mPrimePrimeEtaPOI*dSnk[0][0]-dS13mPrimePrimeEtaPOI;
  dM0pp3EtaPOI=mPrimePrimeEtaPOIHere*dSnk[0][2]-dS13mPrimePrimeEtaPOI;
  dM0pp12EtaPOI=mPrimePrimeEtaPOIHere*dSnk[0][0]*dSnk[0][1]-dM2pp1EtaPOI-dM1pp2EtaPOI-dM0pp3EtaPOI-dS13mPrimePrimeEtaPOI;
  dM0pp111EtaPOI=mPrimePrimeEtaPOIHere*dSnk[2][0]-3.*dM1pp11EtaPOI-3.*dM0pp12EtaPOI-3.*dM2pp1EtaPOI-3.*dM1pp2EtaPOI-dM0pp3EtaPOI-dS13mPrimePrimeEtaPOI;
  
  // q':
  qxPrimeEtaPOI = (etaReq1nPrime->GetBinContent(bin))*(etaReq1nPrime->GetBinEntries(bin));
  qyPrimeEtaPOI = (etaImq1nPrime->GetBinContent(bin))*(etaImq1nPrime->GetBinEntries(bin));
  
  mPrimeEtaPOI = etaReq1nPrime->GetBinEntries(bin); // to be improved
  dM0p111EtaPOI=mPrimeEtaPOI*(dSnk[2][0]-3.*dSnk[0][1]*dSnk[0][0]+2.*dSnk[0][2]);
 
  // 2-p the needed one
  Double_t two1n1nWPerEtaBinPOI=0.; 
  if((mPrimePrimeEtaPOIHere+mPrimeEtaPOI)*dSnk[0][0]-dS11mPrimePrimeEtaPOI>0)
  {
   two1n1nWPerEtaBinPOI = (qxPrimePrimeEtaPOI*dQnkX[0][0]+qyPrimePrimeEtaPOI*dQnkY[0][0]+qxPrimeEtaPOI*dQnkX[0][0]+qyPrimeEtaPOI*dQnkY[0][0]-dS11mPrimePrimeEtaPOI)/((mPrimePrimeEtaPOIHere+mPrimeEtaPOI)*dSnk[0][0]-dS11mPrimePrimeEtaPOI); 
  
   f2WPerEtaBin1n1nPOI->Fill(fEtaMin+(bin-1)*dBinWidthEta,two1n1nWPerEtaBinPOI,(mPrimePrimeEtaPOIHere+mPrimeEtaPOI)*dSnk[0][0]-dS11mPrimePrimeEtaPOI); // <2'>_{n|n} 
  }
 
  // 2-p the temporary one
  Double_t two1n1nW1ppW1W1EtaPOI=0.; // <w1 w2 w3 cos(n(phi2-phi3))> // OK!!!
  if(dM1pp11EtaPOI)
  {
   two1n1nW1ppW1W1EtaPOI = ((pow(dQnkX[0][0],2.)+pow(dQnkY[0][0],2.))*dS11mPrimePrimeEtaPOI-2.*(qxW2PrimePrimeEtaPOI*dQnkX[0][0]+qyW2PrimePrimeEtaPOI*dQnkY[0][0])-dS11mPrimePrimeEtaPOI*dSnk[0][1]+2.*dS13mPrimePrimeEtaPOI)/dM1pp11EtaPOI; // CORRECT !!! <w1 w2 w3 cos(n(phi2-phi3))>
  }
  
  // 2-p the temporary one
  Double_t two1npp1nW1W2EtaPOI=0.; // <w2 w3^2 cos(n(psi1-phi2))> // OK !!!
  if(dM0pp12EtaPOI)
  {
   two1npp1nW1W2EtaPOI = (dSnk[0][1]*(qxPrimePrimeEtaPOI*dQnkX[0][0]+qyPrimePrimeEtaPOI*dQnkY[0][0])-(qxW2PrimePrimeEtaPOI*dQnkX[0][0]+qyW2PrimePrimeEtaPOI*dQnkY[0][0])-dM1pp2EtaPOI-(qxPrimePrimeEtaPOI*dQnkX[0][2]+qyPrimePrimeEtaPOI*dQnkY[0][2])+dS13mPrimePrimeEtaPOI)/dM0pp12EtaPOI; // CORRECT !!! <w2 w3^2 cos(n(psi1-phi2))>
  }
  
  // 2-p the temporary one 
  Double_t two1npp1nW2ppW1EtaPOI=0.; // <w1^2 w2 cos(n(psi1-phi2))> // OK !!!
  if(dM2pp1EtaPOI)
  {
   two1npp1nW2ppW1EtaPOI = ((qxW2PrimePrimeEtaPOI*dQnkX[0][0]+qyW2PrimePrimeEtaPOI*dQnkY[0][0])-dS13mPrimePrimeEtaPOI)/dM2pp1EtaPOI; // CORRECT !!! <w1^2 w2 cos(n(psi1-phi2))> 
  }
  
  // 2-p the temporary one
  Double_t two2npp2nW1ppW2EtaPOI=0.; // <w1 w2^2 cos(2n(psi1-phi2))> OK !!!
  if(dM1pp2EtaPOI)
  {
   two2npp2nW1ppW2EtaPOI = ((q2xW1PrimePrimeEtaPOI*dQnkX[1][1]+q2yW1PrimePrimeEtaPOI*dQnkY[1][1])-dS13mPrimePrimeEtaPOI)/dM1pp2EtaPOI; // CORRECT !!! <w1 w2^2 cos(2n(psi1-phi2))> 
  }
  
  // 2-p the temporary one 
  Double_t two1npp1nW3EtaPOI=0.; // <w2^3 cos(n(psi1-phi2))> // OK !!!
  if(dM0pp3EtaPOI)
  {
   two1npp1nW3EtaPOI = (qxPrimePrimeEtaPOI*dQnkX[0][2]+qyPrimePrimeEtaPOI*dQnkY[0][2]-dS13mPrimePrimeEtaPOI)/dM0pp3EtaPOI; // CORRECT !!! <w2^3 cos(n(psi1-phi2))>
  }
   
  // 3-p the temporary one 
  Double_t three2npp1n1nW1ppW1W1EtaPOI=0.; // <w1 w2 w3 cos(n(2psi1-phi2-phi3))> // OK!!!
  if(dM1pp11EtaPOI)
  {
   three2npp1n1nW1ppW1W1EtaPOI = (q2xW1PrimePrimeEtaPOI*(dQnkX[0][0]*dQnkX[0][0]-dQnkY[0][0]*dQnkY[0][0])+2.*q2yW1PrimePrimeEtaPOI*dQnkX[0][0]*dQnkY[0][0]-2.*(qxW2PrimePrimeEtaPOI*dQnkX[0][0]+qyW2PrimePrimeEtaPOI*dQnkY[0][0])-(q2xW1PrimePrimeEtaPOI*dQnkX[1][1]+q2yW1PrimePrimeEtaPOI*dQnkY[1][1])+2.*dS13mPrimePrimeEtaPOI)/dM1pp11EtaPOI; // CORRECT !!! <w1 w2 w3 cos(n(2psi1-phi2-phi3))> 
  }
  
  // 3-p the temporary one 
  Double_t three1npp1n2nW0ppW1W2EtaPOI=0.; // <w2 w3^2 cos(n(psi1+phi2-2*phi3))> // OK!!!
  if(dM0pp12EtaPOI)
  {
   three1npp1n2nW0ppW1W2EtaPOI = (qxPrimePrimeEtaPOI*(dQnkX[0][0]*dQnkX[1][1]+dQnkY[0][0]*dQnkY[1][1])-qyPrimePrimeEtaPOI*(dQnkY[0][0]*dQnkX[1][1]-dQnkX[0][0]*dQnkY[1][1])-(qxW2PrimePrimeEtaPOI*dQnkX[0][0]+qyW2PrimePrimeEtaPOI*dQnkY[0][0])-(q2xW1PrimePrimeEtaPOI*dQnkX[1][1]+q2yW1PrimePrimeEtaPOI*dQnkY[1][1])-(qxPrimePrimeEtaPOI*dQnkX[0][2]+qyPrimePrimeEtaPOI*dQnkY[0][2])+2.*dS13mPrimePrimeEtaPOI)/dM0pp12EtaPOI; // CORRECT !!! <w2 w3^2 cos(n(psi1+phi2-2.*phi3))>
  }
 
  /*
  // 4-p RP part
  Double_t four1npp1n1n1nW1W1W1=0.; // <w1 w2 w3 cos(n(psi1+phi1-phi2-phi3))>
  if(dM0pp111EtaPOI)
  {   
   four1npp1n1n1nW1W1W1 = ((pow(dQnkX[0][0],2.)+pow(dQnkY[0][0],2.))*(qxPrimePrimeEtaPOI*dQnkX[0][0]+qyPrimePrimeEtaPOI*dQnkY[0][0])-2.*dM1pp11EtaPOI*two1n1nW1ppW1W1EtaPOI-dM1pp11EtaPOI*three2npp1n1nW1ppW1W1EtaPOI-dM0pp12EtaPOI*three1npp1n2nW0ppW1W2EtaPOI-2.*dM0pp12EtaPOI*two1npp1nW1W2EtaPOI-3.*dM2pp1EtaPOI*two1npp1nW2ppW1EtaPOI-2.*dM1pp2EtaPOI-dM1pp2EtaPOI*two2npp2nW1ppW2EtaPOI-dM0pp3EtaPOI*two1npp1nW3EtaPOI-dS13mPrimePrimeEtaPOI)/(dM0pp111EtaPOI); 
  }
  */
  /*
  // 4-p POI part
  Double_t four1npp1n1n1nW1W1W1POI=0.;
  if(dM0p111EtaPOI>0&&mPrimeEtaPOI>0&&nRP>0)
  {
   four1npp1n1n1nW1W1W1POI = ((pow(dQnkX[0][0],2.)+pow(dQnkY[0][0],2.))*(qxPrimeEtaPOI*dQnkX[0][0]+qyPrimeEtaPOI*dQnkY[0][0])-2.*dSnk[0][1]* (qxPrimeEtaPOI*dQnkX[0][0]+qyPrimeEtaPOI*dQnkY[0][0])+2.*(qxPrimeEtaPOI*dQnkX[0][2]+qyPrimeEtaPOI*dQnkY[0][2])-qxPrimeEtaPOI*(dQnkX[0][0]*dQnkX[1][1]+dQnkY[0][0]*dQnkY[1][1])+qyPrimeEtaPOI*(dQnkY[0][0]*dQnkX[1][1]-dQnkX[0][0]*dQnkY[1][1]))/dM0p111EtaPOI;
  }
  */
 
  // 4-p RP and POI in all combinations (full, partial and no overlap)
  Double_t four1npp1n1n1nW1W1W1EtaPOI=0.;
 
  if(dM0pp111EtaPOI+dM0p111EtaPOI)
  {
   four1npp1n1n1nW1W1W1EtaPOI = ((pow(dQnkX[0][0],2.)+pow(dQnkY[0][0],2.))*(qxPrimePrimeEtaPOI*dQnkX[0][0]+qyPrimePrimeEtaPOI*dQnkY[0][0])-2.*dM1pp11EtaPOI*two1n1nW1ppW1W1EtaPOI-dM1pp11EtaPOI*three2npp1n1nW1ppW1W1EtaPOI-dM0pp12EtaPOI*three1npp1n2nW0ppW1W2EtaPOI-2.*dM0pp12EtaPOI*two1npp1nW1W2EtaPOI-3.*dM2pp1EtaPOI*two1npp1nW2ppW1EtaPOI-2.*dM1pp2EtaPOI-dM1pp2EtaPOI*two2npp2nW1ppW2EtaPOI-dM0pp3EtaPOI*two1npp1nW3EtaPOI-dS13mPrimePrimeEtaPOI+(pow(dQnkX[0][0],2.)+pow(dQnkY[0][0],2.))*(qxPrimeEtaPOI*dQnkX[0][0]+qyPrimeEtaPOI*dQnkY[0][0])-2.*dSnk[0][1]*(qxPrimeEtaPOI*dQnkX[0][0]+qyPrimeEtaPOI*dQnkY[0][0])+2.*(qxPrimeEtaPOI*dQnkX[0][2]+qyPrimeEtaPOI*dQnkY[0][2])-qxPrimeEtaPOI*(dQnkX[0][0]*dQnkX[1][1]+dQnkY[0][0]*dQnkY[1][1])+qyPrimeEtaPOI*(dQnkY[0][0]*dQnkX[1][1]-dQnkX[0][0]*dQnkY[1][1]))/(dM0pp111EtaPOI+dM0p111EtaPOI);
   
   f4WPerEtaBin1n1n1n1nPOI->Fill(fEtaMin+(bin-1)*dBinWidthEta,four1npp1n1n1nW1W1W1EtaPOI,dM0pp111EtaPOI+dM0p111EtaPOI);
  }
 }
 //...........................................................................................................    
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 // RPs
 ptReq1nPrime->Reset(); // to be improved 
 ptImq1nPrime->Reset(); // to be improved
 ptReq2nPrime->Reset(); // to be improved 
 ptImq2nPrime->Reset(); // to be improved
 
 etaReq1nPrime->Reset(); // to be improved 
 etaImq1nPrime->Reset(); // to be improved
 etaReq2nPrime->Reset(); // to be improved 
 etaImq2nPrime->Reset(); // to be improved

 
 //...........................................................................................................
 // PrimePrime Pt RP
 Double_t qxPtRP=0.,qyPtRP=0.,q2xPtRP=0.,q2yPtRP=0.;//add comments for these variable
 Double_t qxW2PtRP=0.,qyW2PtRP=0.,q2xW1PtRP=0.,q2yW1PtRP=0.;//add comments for these variable
 Double_t dS11mPtRP=0.; // to be improved (name)
 Double_t dS12mPtRP=0.; // to be improved (name)
 Double_t dS13mPtRP=0.; // to be improved (name)
 Double_t mPtRP=0.; // to be improved (name)

 Double_t dM1pp11PtRP=0.; // to be improved (name)
 Double_t dM0pp111PtRP=0.; // to be improved (name)
 Double_t dM0pp12PtRP=0.;
 Double_t dM2pp1PtRP=0.;
 Double_t dM1pp2PtRP=0.;
 Double_t dM0pp3PtRP=0.;
 
 // Prime Pt RP
 Double_t qxPrimePtRP=0.,qyPrimePtRP=0.;
 Double_t mPrimePtRP=0.;
 Double_t dM0p111PtRP=0.; // to be improved (name)
  
 for(Int_t bin=1;bin<(fnBinsPt+1);bin++) // loop over pt-bins 
 {     
  // q'':                    
  qxPtRP = (ptReq1n->GetBinContent(bin))*(ptReq1n->GetBinEntries(bin));
  qyPtRP = (ptImq1n->GetBinContent(bin))*(ptImq1n->GetBinEntries(bin)); 
  q2xPtRP = (ptReq2n->GetBinContent(bin))*(ptReq2n->GetBinEntries(bin));  
  q2yPtRP = (ptImq2n->GetBinContent(bin))*(ptImq2n->GetBinEntries(bin)); 
  
  qxW2PtRP = (req1nW2Pt->GetBinContent(bin))*(req1nW2Pt->GetBinEntries(bin));
  qyW2PtRP = (imq1nW2Pt->GetBinContent(bin))*(imq1nW2Pt->GetBinEntries(bin));
  
  q2xW1PtRP = (req2nW1Pt->GetBinContent(bin))*(req2nW1Pt->GetBinEntries(bin));
  q2yW1PtRP = (imq2nW1Pt->GetBinContent(bin))*(imq2nW1Pt->GetBinEntries(bin));
  
  dS11mPtRP = (sumOfW1upTomPt->GetBinContent(bin))*(sumOfW1upTomPt->GetBinEntries(bin));
  dS12mPtRP = (sumOfW2upTomPt->GetBinContent(bin))*(sumOfW2upTomPt->GetBinEntries(bin));
  dS13mPtRP = (sumOfW3upTomPt->GetBinContent(bin))*(sumOfW3upTomPt->GetBinEntries(bin));

  mPtRP = sumOfW1upTomPt->GetBinEntries(bin); // to be improved
      
  dM1pp11PtRP=dS11mPtRP*(dSnk[1][0]-dSnk[0][1])-2.*dS12mPtRP*dSnk[0][0]+2.*dS13mPtRP;
  dM1pp2PtRP=dS11mPtRP*dSnk[0][1]-dS13mPtRP;
  dM2pp1PtRP=dS12mPtRP*dSnk[0][0]-dS13mPtRP;
  dM0pp3PtRP=mPtRP*dSnk[0][2]-dS13mPtRP;
  dM0pp12PtRP=mPtRP*dSnk[0][0]*dSnk[0][1]-dM2pp1PtRP-dM1pp2PtRP-dM0pp3PtRP-dS13mPtRP;
  dM0pp111PtRP=mPtRP*dSnk[2][0]-3.*dM1pp11PtRP-3.*dM0pp12PtRP-3.*dM2pp1PtRP-3.*dM1pp2PtRP-dM0pp3PtRP-dS13mPtRP;
  
  // q':
  qxPrimePtRP = (ptReq1nPrime->GetBinContent(bin))*(ptReq1nPrime->GetBinEntries(bin));
  qyPrimePtRP = (ptImq1nPrime->GetBinContent(bin))*(ptImq1nPrime->GetBinEntries(bin));
  
  mPrimePtRP = ptReq1nPrime->GetBinEntries(bin); // to be improved
  dM0p111PtRP=mPrimePtRP*(dSnk[2][0]-3.*dSnk[0][1]*dSnk[0][0]+2.*dSnk[0][2]);
  
  // 2-p the needed one
  Double_t two1n1nWPerPtBinRP=0.; 
  if((mPtRP+mPrimePtRP)*dSnk[0][0]-dS11mPtRP>0)
  {
   two1n1nWPerPtBinRP = (qxPtRP*dQnkX[0][0]+qyPtRP*dQnkY[0][0]+qxPrimePtRP*dQnkX[0][0]+qyPrimePtRP*dQnkY[0][0]-dS11mPtRP)/((mPtRP+mPrimePtRP)*dSnk[0][0]-dS11mPtRP); 
   f2WPerPtBin1n1nRP->Fill(fPtMin+(bin-1)*dBinWidthPt,two1n1nWPerPtBinRP,(mPtRP+mPrimePtRP)*dSnk[0][0]-dS11mPtRP);  
  }
 
  // 2-p temporary one
  Double_t two1n1nW1ppW1W1PtRP=0.; // <w1 w2 w3 cos(n(phi2-phi3))> // OK!!!
  if(dM1pp11PtRP)
  {
   two1n1nW1ppW1W1PtRP = ((pow(dQnkX[0][0],2.)+pow(dQnkY[0][0],2.))*dS11mPtRP-2.*(qxW2PtRP*dQnkX[0][0]+qyW2PtRP*dQnkY[0][0])-dS11mPtRP*dSnk[0][1]+2.*dS13mPtRP)/dM1pp11PtRP; // CORRECT !!! <w1 w2 w3 cos(n(phi2-phi3))>
  }
  
  // 2-p temporary one
  Double_t two1npp1nW1W2PtRP=0.; // <w2 w3^2 cos(n(psi1-phi2))> // OK !!!
  if(dM0pp12PtRP)
  {
   two1npp1nW1W2PtRP = (dSnk[0][1]*(qxPtRP*dQnkX[0][0]+qyPtRP*dQnkY[0][0])-(qxW2PtRP*dQnkX[0][0]+qyW2PtRP*dQnkY[0][0])-dM1pp2PtRP-(qxPtRP*dQnkX[0][2]+qyPtRP*dQnkY[0][2])+dS13mPtRP)/dM0pp12PtRP; // CORRECT !!! <w2 w3^2 cos(n(psi1-phi2))>
  }
  
  // 2-p temporary one
  Double_t two1npp1nW2ppW1PtRP=0.; // <w1^2 w2 cos(n(psi1-phi2))> // OK !!!
  if(dM2pp1PtRP)
  {
   two1npp1nW2ppW1PtRP = ((qxW2PtRP*dQnkX[0][0]+qyW2PtRP*dQnkY[0][0])-dS13mPtRP)/dM2pp1PtRP; // CORRECT !!! <w1^2 w2 cos(n(psi1-phi2))> 
  }
  
  // 2-p temporary one
  Double_t two2npp2nW1ppW2PtRP=0.; // <w1 w2^2 cos(2n(psi1-phi2))> OK !!!
  if(dM1pp2PtRP)
  {
   two2npp2nW1ppW2PtRP = ((q2xW1PtRP*dQnkX[1][1]+q2yW1PtRP*dQnkY[1][1])-dS13mPtRP)/dM1pp2PtRP; // CORRECT !!! <w1 w2^2 cos(2n(psi1-phi2))> 
  }
  
  // 2-p temporary one
  Double_t two1npp1nW3PtRP=0.; // <w2^3 cos(n(psi1-phi2))> // OK !!!
  if(dM0pp3PtRP)
  {
   two1npp1nW3PtRP = (qxPtRP*dQnkX[0][2]+qyPtRP*dQnkY[0][2]-dS13mPtRP)/dM0pp3PtRP; // CORRECT !!! <w2^3 cos(n(psi1-phi2))>
  }
   
  // 3-p temporary one
  Double_t three2npp1n1nW1ppW1W1PtRP=0.; // <w1 w2 w3 cos(n(2psi1-phi2-phi3))> // OK!!!
  if(dM1pp11PtRP)
  {
   three2npp1n1nW1ppW1W1PtRP = (q2xW1PtRP*(dQnkX[0][0]*dQnkX[0][0]-dQnkY[0][0]*dQnkY[0][0])+2.*q2yW1PtRP*dQnkX[0][0]*dQnkY[0][0]-2.*(qxW2PtRP*dQnkX[0][0]+qyW2PtRP*dQnkY[0][0])-(q2xW1PtRP*dQnkX[1][1]+q2yW1PtRP*dQnkY[1][1])+2.*dS13mPtRP)/dM1pp11PtRP; // CORRECT !!! <w1 w2 w3 cos(n(2psi1-phi2-phi3))> 
  }
  
  // 3-p temporary one
  Double_t three1npp1n2nW0ppW1W2PtRP=0.; // <w2 w3^2 cos(n(psi1+phi2-2*phi3))> // OK!!!
  if(dM0pp12PtRP)
  {
   three1npp1n2nW0ppW1W2PtRP = (qxPtRP*(dQnkX[0][0]*dQnkX[1][1]+dQnkY[0][0]*dQnkY[1][1])-qyPtRP*(dQnkY[0][0]*dQnkX[1][1]-dQnkX[0][0]*dQnkY[1][1])-(qxW2PtRP*dQnkX[0][0]+qyW2PtRP*dQnkY[0][0])-(q2xW1PtRP*dQnkX[1][1]+q2yW1PtRP*dQnkY[1][1])-(qxPtRP*dQnkX[0][2]+qyPtRP*dQnkY[0][2])+2.*dS13mPtRP)/dM0pp12PtRP; // CORRECT !!! <w2 w3^2 cos(n(psi1+phi2-2.*phi3))>
  }
 
  /*
  // 4-p RP part
  Double_t four1npp1n1n1nW1W1W1=0.; // <w1 w2 w3 cos(n(psi1+phi1-phi2-phi3))>
  if(dM0pp111PtRP)
  {   
   four1npp1n1n1nW1W1W1 = ((pow(dQnkX[0][0],2.)+pow(dQnkY[0][0],2.))*(qxPtRP*dQnkX[0][0]+qyPtRP*dQnkY[0][0])-2.*dM1pp11PtRP*two1n1nW1ppW1W1PtRP-dM1pp11PtRP*three2npp1n1nW1ppW1W1PtRP-dM0pp12PtRP*three1npp1n2nW0ppW1W2PtRP-2.*dM0pp12PtRP*two1npp1nW1W2PtRP-3.*dM2pp1PtRP*two1npp1nW2ppW1PtRP-2.*dM1pp2PtRP-dM1pp2PtRP*two2npp2nW1ppW2PtRP-dM0pp3PtRP*two1npp1nW3PtRP-dS13mPtRP)/(dM0pp111PtRP); 
  }
  */
  
  /*
  // 4-p POI part
  Double_t four1npp1n1n1nW1W1W1RP=0.;
  if(dM0p111PtRP>0&&mPrimePtRP>0&&nRP>0)
  {
   four1npp1n1n1nW1W1W1RP = ((pow(dQnkX[0][0],2.)+pow(dQnkY[0][0],2.))*(qxPrimePtRP*dQnkX[0][0]+qyPrimePtRP*dQnkY[0][0])-2.*dSnk[0][1]* (qxPrimePtRP*dQnkX[0][0]+qyPrimePtRP*dQnkY[0][0])+2.*(qxPrimePtRP*dQnkX[0][2]+qyPrimePtRP*dQnkY[0][2])-qxPrimePtRP*(dQnkX[0][0]*dQnkX[1][1]+dQnkY[0][0]*dQnkY[1][1])+qyPrimePtRP*(dQnkY[0][0]*dQnkX[1][1]-dQnkX[0][0]*dQnkY[1][1]))/dM0p111PtRP;
  }
  */
 
  // 4-p RP and POI in all combinations (full, partial and no overlap)
  Double_t four1npp1n1n1nW1W1W1PtRP=0.;
  if(dM0pp111PtRP+dM0p111PtRP)
  {
  four1npp1n1n1nW1W1W1PtRP = ((pow(dQnkX[0][0],2.)+pow(dQnkY[0][0],2.))*(qxPtRP*dQnkX[0][0]+qyPtRP*dQnkY[0][0])-2.*dM1pp11PtRP*two1n1nW1ppW1W1PtRP-dM1pp11PtRP*three2npp1n1nW1ppW1W1PtRP-dM0pp12PtRP*three1npp1n2nW0ppW1W2PtRP-2.*dM0pp12PtRP*two1npp1nW1W2PtRP-3.*dM2pp1PtRP*two1npp1nW2ppW1PtRP-2.*dM1pp2PtRP-dM1pp2PtRP*two2npp2nW1ppW2PtRP-dM0pp3PtRP*two1npp1nW3PtRP-dS13mPtRP+(pow(dQnkX[0][0],2.)+pow(dQnkY[0][0],2.))*(qxPrimePtRP*dQnkX[0][0]+qyPrimePtRP*dQnkY[0][0])-2.*dSnk[0][1]* (qxPrimePtRP*dQnkX[0][0]+qyPrimePtRP*dQnkY[0][0])+2.*(qxPrimePtRP*dQnkX[0][2]+qyPrimePtRP*dQnkY[0][2])-qxPrimePtRP*(dQnkX[0][0]*dQnkX[1][1]+dQnkY[0][0]*dQnkY[1][1])+qyPrimePtRP*(dQnkY[0][0]*dQnkX[1][1]-dQnkX[0][0]*dQnkY[1][1]))/(dM0pp111PtRP+dM0p111PtRP);
 
   f4WPerPtBin1n1n1n1nRP->Fill(fPtMin+(bin-1)*dBinWidthPt,four1npp1n1n1nW1W1W1PtRP,dM0pp111PtRP+dM0p111PtRP);
  } // end of if(dM0pp111PtRP+dM0p111PtRP)
 } // for(Int_t bin=1;bin<(fnBinsPt+1);bin++) // loop over pt-bins
 
 delete ptReq1nPrime;
 delete ptImq1nPrime;
 delete ptReq2nPrime;
 delete ptImq2nPrime;
 delete ptReq1n;
 delete ptImq1n;
 delete ptReq2n;
 delete ptImq2n;
 
 delete req1nW2Pt;
 delete imq1nW2Pt;
 delete req2nW1Pt;
 delete imq2nW1Pt;
 delete sumOfW1upTomPt;
 delete sumOfW2upTomPt;
 delete sumOfW3upTomPt;
 //...........................................................................................................

 //...........................................................................................................  
 // PrimePrime Eta RP
 Double_t qxEtaRP=0.,qyEtaRP=0.,q2xEtaRP=0.,q2yEtaRP=0.;//add comments for these variable
 Double_t qxW2EtaRP=0.,qyW2EtaRP=0.,q2xW1EtaRP=0.,q2yW1EtaRP=0.;//add comments for these variable
 Double_t dS11mEtaRP=0.; // to be improved (name)
 Double_t dS12mEtaRP=0.; // to be improved (name)
 Double_t dS13mEtaRP=0.; // to be improved (name)
 Double_t mEtaRPHere=0.; // to be improved (name)

 Double_t dM1pp11EtaRP=0.; // to be improved (name)
 Double_t dM0pp111EtaRP=0.; // to be improved (name)
 Double_t dM0pp12EtaRP=0.;
 Double_t dM2pp1EtaRP=0.;
 Double_t dM1pp2EtaRP=0.;
 Double_t dM0pp3EtaRP=0.;
 
 // Prime Eta RP
 Double_t qxPrimeEtaRPHere=0.,qyPrimeEtaRPHere=0.;
 Double_t mPrimeEtaRPHere=0.;
 Double_t dM0p111EtaRP=0.; // to be improved (name)
 
 for(Int_t bin=1;bin<(fnBinsEta+1);bin++) // loop over eta-bins 
 {     
  // q'':                    
  qxEtaRP = (etaReq1n->GetBinContent(bin))*(etaReq1n->GetBinEntries(bin));
  qyEtaRP = (etaImq1n->GetBinContent(bin))*(etaImq1n->GetBinEntries(bin)); 
  q2xEtaRP = (etaReq2n->GetBinContent(bin))*(etaReq2n->GetBinEntries(bin));  
  q2yEtaRP = (etaImq2n->GetBinContent(bin))*(etaImq2n->GetBinEntries(bin)); 
  
  qxW2EtaRP = (req1nW2Eta->GetBinContent(bin))*(req1nW2Eta->GetBinEntries(bin));
  qyW2EtaRP = (imq1nW2Eta->GetBinContent(bin))*(imq1nW2Eta->GetBinEntries(bin));
  
  q2xW1EtaRP = (req2nW1Eta->GetBinContent(bin))*(req2nW1Eta->GetBinEntries(bin));
  q2yW1EtaRP = (imq2nW1Eta->GetBinContent(bin))*(imq2nW1Eta->GetBinEntries(bin));
  
  dS11mEtaRP = (sumOfW1upTomEta->GetBinContent(bin))*(sumOfW1upTomEta->GetBinEntries(bin));
  dS12mEtaRP = (sumOfW2upTomEta->GetBinContent(bin))*(sumOfW2upTomEta->GetBinEntries(bin));
  dS13mEtaRP = (sumOfW3upTomEta->GetBinContent(bin))*(sumOfW3upTomEta->GetBinEntries(bin));

  mEtaRPHere = sumOfW1upTomEta->GetBinEntries(bin); // to be improved
      
  dM1pp11EtaRP=dS11mEtaRP*(dSnk[1][0]-dSnk[0][1])-2.*dS12mEtaRP*dSnk[0][0]+2.*dS13mEtaRP;
  dM1pp2EtaRP=dS11mEtaRP*dSnk[0][1]-dS13mEtaRP;
  dM2pp1EtaRP=dS12mEtaRP*dSnk[0][0]-dS13mEtaRP;
  dM0pp3EtaRP=mEtaRPHere*dSnk[0][2]-dS13mEtaRP;
  dM0pp12EtaRP=mEtaRPHere*dSnk[0][0]*dSnk[0][1]-dM2pp1EtaRP-dM1pp2EtaRP-dM0pp3EtaRP-dS13mEtaRP;
  dM0pp111EtaRP=mEtaRPHere*dSnk[2][0]-3.*dM1pp11EtaRP-3.*dM0pp12EtaRP-3.*dM2pp1EtaRP-3.*dM1pp2EtaRP-dM0pp3EtaRP-dS13mEtaRP;
  
  // q':
  qxPrimeEtaRPHere = (etaReq1nPrime->GetBinContent(bin))*(etaReq1nPrime->GetBinEntries(bin));
  qyPrimeEtaRPHere = (etaImq1nPrime->GetBinContent(bin))*(etaImq1nPrime->GetBinEntries(bin));
  
  mPrimeEtaRPHere = etaReq1nPrime->GetBinEntries(bin); // to be improved
  dM0p111EtaRP=mPrimeEtaRPHere*(dSnk[2][0]-3.*dSnk[0][1]*dSnk[0][0]+2.*dSnk[0][2]);
 
  // 2-p the needed one
  Double_t two1n1nWPerEtaBinRP=0.; 
  if((mEtaRPHere+mPrimeEtaRPHere)*dSnk[0][0]-dS11mEtaRP>0)
  {
   two1n1nWPerEtaBinRP = (qxEtaRP*dQnkX[0][0]+qyEtaRP*dQnkY[0][0]+qxPrimeEtaRPHere*dQnkX[0][0]+qyPrimeEtaRPHere*dQnkY[0][0]-dS11mEtaRP)/((mEtaRPHere+mPrimeEtaRPHere)*dSnk[0][0]-dS11mEtaRP); 
  
   f2WPerEtaBin1n1nRP->Fill(fEtaMin+(bin-1)*dBinWidthEta,two1n1nWPerEtaBinRP,(mEtaRPHere+mPrimeEtaRPHere)*dSnk[0][0]-dS11mEtaRP); // <2'>_{n|n} 
  }
 
  // 2-p the temporary one
  Double_t two1n1nW1ppW1W1EtaRP=0.; // <w1 w2 w3 cos(n(phi2-phi3))> // OK!!!
  if(dM1pp11EtaRP)
  {
   two1n1nW1ppW1W1EtaRP = ((pow(dQnkX[0][0],2.)+pow(dQnkY[0][0],2.))*dS11mEtaRP-2.*(qxW2EtaRP*dQnkX[0][0]+qyW2EtaRP*dQnkY[0][0])-dS11mEtaRP*dSnk[0][1]+2.*dS13mEtaRP)/dM1pp11EtaRP; // CORRECT !!! <w1 w2 w3 cos(n(phi2-phi3))>
  }
  
  // 2-p the temporary one
  Double_t two1npp1nW1W2EtaRP=0.; // <w2 w3^2 cos(n(psi1-phi2))> // OK !!!
  if(dM0pp12EtaRP)
  {
   two1npp1nW1W2EtaRP = (dSnk[0][1]*(qxEtaRP*dQnkX[0][0]+qyEtaRP*dQnkY[0][0])-(qxW2EtaRP*dQnkX[0][0]+qyW2EtaRP*dQnkY[0][0])-dM1pp2EtaRP-(qxEtaRP*dQnkX[0][2]+qyEtaRP*dQnkY[0][2])+dS13mEtaRP)/dM0pp12EtaRP; // CORRECT !!! <w2 w3^2 cos(n(psi1-phi2))>
  }
  
  // 2-p the temporary one 
  Double_t two1npp1nW2ppW1EtaRP=0.; // <w1^2 w2 cos(n(psi1-phi2))> // OK !!!
  if(dM2pp1EtaRP)
  {
   two1npp1nW2ppW1EtaRP = ((qxW2EtaRP*dQnkX[0][0]+qyW2EtaRP*dQnkY[0][0])-dS13mEtaRP)/dM2pp1EtaRP; // CORRECT !!! <w1^2 w2 cos(n(psi1-phi2))> 
  }
  
  // 2-p the temporary one
  Double_t two2npp2nW1ppW2EtaRP=0.; // <w1 w2^2 cos(2n(psi1-phi2))> OK !!!
  if(dM1pp2EtaRP)
  {
   two2npp2nW1ppW2EtaRP = ((q2xW1EtaRP*dQnkX[1][1]+q2yW1EtaRP*dQnkY[1][1])-dS13mEtaRP)/dM1pp2EtaRP; // CORRECT !!! <w1 w2^2 cos(2n(psi1-phi2))> 
  }
  
  // 2-p the temporary one 
  Double_t two1npp1nW3EtaRP=0.; // <w2^3 cos(n(psi1-phi2))> // OK !!!
  if(dM0pp3EtaRP)
  {
   two1npp1nW3EtaRP = (qxEtaRP*dQnkX[0][2]+qyEtaRP*dQnkY[0][2]-dS13mEtaRP)/dM0pp3EtaRP; // CORRECT !!! <w2^3 cos(n(psi1-phi2))>
  }
   
  // 3-p the temporary one 
  Double_t three2npp1n1nW1ppW1W1EtaRP=0.; // <w1 w2 w3 cos(n(2psi1-phi2-phi3))> // OK!!!
  if(dM1pp11EtaRP)
  {
   three2npp1n1nW1ppW1W1EtaRP = (q2xW1EtaRP*(dQnkX[0][0]*dQnkX[0][0]-dQnkY[0][0]*dQnkY[0][0])+2.*q2yW1EtaRP*dQnkX[0][0]*dQnkY[0][0]-2.*(qxW2EtaRP*dQnkX[0][0]+qyW2EtaRP*dQnkY[0][0])-(q2xW1EtaRP*dQnkX[1][1]+q2yW1EtaRP*dQnkY[1][1])+2.*dS13mEtaRP)/dM1pp11EtaRP; // CORRECT !!! <w1 w2 w3 cos(n(2psi1-phi2-phi3))> 
  }
  
  // 3-p the temporary one 
  Double_t three1npp1n2nW0ppW1W2EtaRP=0.; // <w2 w3^2 cos(n(psi1+phi2-2*phi3))> // OK!!!
  if(dM0pp12EtaRP)
  {
   three1npp1n2nW0ppW1W2EtaRP = (qxEtaRP*(dQnkX[0][0]*dQnkX[1][1]+dQnkY[0][0]*dQnkY[1][1])-qyEtaRP*(dQnkY[0][0]*dQnkX[1][1]-dQnkX[0][0]*dQnkY[1][1])-(qxW2EtaRP*dQnkX[0][0]+qyW2EtaRP*dQnkY[0][0])-(q2xW1EtaRP*dQnkX[1][1]+q2yW1EtaRP*dQnkY[1][1])-(qxEtaRP*dQnkX[0][2]+qyEtaRP*dQnkY[0][2])+2.*dS13mEtaRP)/dM0pp12EtaRP; // CORRECT !!! <w2 w3^2 cos(n(psi1+phi2-2.*phi3))>
  }
 
  /*
  // 4-p RP part
  Double_t four1npp1n1n1nW1W1W1=0.; // <w1 w2 w3 cos(n(psi1+phi1-phi2-phi3))>
  if(dM0pp111EtaRP)
  {   
   four1npp1n1n1nW1W1W1 = ((pow(dQnkX[0][0],2.)+pow(dQnkY[0][0],2.))*(qxEtaRP*dQnkX[0][0]+qyEtaRP*dQnkY[0][0])-2.*dM1pp11EtaRP*two1n1nW1ppW1W1EtaRP-dM1pp11EtaRP*three2npp1n1nW1ppW1W1EtaRP-dM0pp12EtaRP*three1npp1n2nW0ppW1W2EtaRP-2.*dM0pp12EtaRP*two1npp1nW1W2EtaRP-3.*dM2pp1EtaRP*two1npp1nW2ppW1EtaRP-2.*dM1pp2EtaRP-dM1pp2EtaRP*two2npp2nW1ppW2EtaRP-dM0pp3EtaRP*two1npp1nW3EtaRP-dS13mEtaRP)/(dM0pp111EtaRP); 
  }
  */
  /*
  // 4-p POI part
  Double_t four1npp1n1n1nW1W1W1RP=0.;
  if(dM0p111EtaRP>0&&mPrimeEtaRPHere>0&&nRP>0)
  {
   four1npp1n1n1nW1W1W1RP = ((pow(dQnkX[0][0],2.)+pow(dQnkY[0][0],2.))*(qxPrimeEtaRPHere*dQnkX[0][0]+qyPrimeEtaRPHere*dQnkY[0][0])-2.*dSnk[0][1]* (qxPrimeEtaRPHere*dQnkX[0][0]+qyPrimeEtaRPHere*dQnkY[0][0])+2.*(qxPrimeEtaRPHere*dQnkX[0][2]+qyPrimeEtaRPHere*dQnkY[0][2])-qxPrimeEtaRPHere*(dQnkX[0][0]*dQnkX[1][1]+dQnkY[0][0]*dQnkY[1][1])+qyPrimeEtaRPHere*(dQnkY[0][0]*dQnkX[1][1]-dQnkX[0][0]*dQnkY[1][1]))/dM0p111EtaRP;
  }
  */
 
  // 4-p RP and POI in all combinations (full, partial and no overlap)
  Double_t four1npp1n1n1nW1W1W1EtaRP=0.;
 
  if(dM0pp111EtaRP+dM0p111EtaRP)
  {
   four1npp1n1n1nW1W1W1EtaRP = ((pow(dQnkX[0][0],2.)+pow(dQnkY[0][0],2.))*(qxEtaRP*dQnkX[0][0]+qyEtaRP*dQnkY[0][0])-2.*dM1pp11EtaRP*two1n1nW1ppW1W1EtaRP-dM1pp11EtaRP*three2npp1n1nW1ppW1W1EtaRP-dM0pp12EtaRP*three1npp1n2nW0ppW1W2EtaRP-2.*dM0pp12EtaRP*two1npp1nW1W2EtaRP-3.*dM2pp1EtaRP*two1npp1nW2ppW1EtaRP-2.*dM1pp2EtaRP-dM1pp2EtaRP*two2npp2nW1ppW2EtaRP-dM0pp3EtaRP*two1npp1nW3EtaRP-dS13mEtaRP+(pow(dQnkX[0][0],2.)+pow(dQnkY[0][0],2.))*(qxPrimeEtaRPHere*dQnkX[0][0]+qyPrimeEtaRPHere*dQnkY[0][0])-2.*dSnk[0][1]*(qxPrimeEtaRPHere*dQnkX[0][0]+qyPrimeEtaRPHere*dQnkY[0][0])+2.*(qxPrimeEtaRPHere*dQnkX[0][2]+qyPrimeEtaRPHere*dQnkY[0][2])-qxPrimeEtaRPHere*(dQnkX[0][0]*dQnkX[1][1]+dQnkY[0][0]*dQnkY[1][1])+qyPrimeEtaRPHere*(dQnkY[0][0]*dQnkX[1][1]-dQnkX[0][0]*dQnkY[1][1]))/(dM0pp111EtaRP+dM0p111EtaRP);
   
   f4WPerEtaBin1n1n1n1nRP->Fill(fEtaMin+(bin-1)*dBinWidthEta,four1npp1n1n1nW1W1W1EtaRP,dM0pp111EtaRP+dM0p111EtaRP);
  }
 }

 delete etaReq1nPrime;
 delete etaImq1nPrime;
 delete etaReq2nPrime;
 delete etaImq2nPrime;
 delete etaReq1n;
 delete etaImq1n;
 delete etaReq2n;
 delete etaImq2n;

 delete req1nW2Eta;
 delete imq1nW2Eta;
 delete req2nW1Eta;
 delete imq2nW1Eta;
 delete sumOfW1upTomEta;
 delete sumOfW2upTomEta;
 delete sumOfW3upTomEta;
 //........................................................................................................... 
 

















 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 

























 /*

 //if(bNestedLoops)to be improved
 //{ to be improved               
 //-------------------------------------------------------------------------------------------------------------------------------- 
 //
 //                                          **********************
 //                                          **** NESTED LOOPS ****
 //                                          ********************** 
 //
 // Remark 1: multi-particle correlations calculated with nested loops are stored in fDirectCorrelations.
 // Remark 2: binning of fDirectCorrelations: bins 0..100 - correlations needed for integrated flow; bins 100..200 - correlations needed for differential flow (taking as an example bin 0.5 < pt < 0.6)
 //
 // binning details of fDirectCorrelations (integrated flow):
 //..........................................................................
 // 1st bin: weighted <2>_{n|n}   = <w1   w2   cos( n*(phi1-phi2))>
 // 2nd bin: weighted <2>_{2n|2n} = <w1^2 w2^2 cos(2n*(phi1-phi2))>
 // 3rd bin: weighted <2>_{3n|3n} = <w1^3 w2^3 cos(3n*(phi1-phi2))>
 // 4th bin: weighted <2>_{4n|4n} = <w1^4 w2^4 cos(4n*(phi1-phi2))>
 // 5th bin: weighted <2>_{n|n} = <w1^3 w2 cos(n*(phi1-phi2))>
 
 // 11th bin: weighted <3>_{2n|n,n} = <w1^2 w2 w3 cos(n*(2phi1-phi2-phi3))>
 //..........................................................................
 

 Double_t phi1=0., phi2=0., phi3=0., phi4=0., phi5=0., phi6=0., phi7=0., phi8=0.;
 Double_t wPhi1=1., wPhi2=1., wPhi3=1., wPhi4=1., wPhi5=1., wPhi6=1., wPhi7=1., wPhi8=1.;



 


 Double_t tempLoop = 0.;
 Int_t tempCounter = 0;
 
 
 
 
 
 for(Int_t i1=0;i1<dMult;i1++)
 {
  fTrack=anEvent->GetTrack(i1);
  phi1=fTrack->Phi();
  if(phiWeights) wPhi1 = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi1*nBinsPhi/TMath::TwoPi())));
  for(Int_t i2=0;i2<dMult;i2++)
  {
   if(i2==i1)continue;
   fTrack=anEvent->GetTrack(i2);
   phi2=fTrack->Phi();
   if(phiWeights) wPhi2 = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi2*nBinsPhi/TMath::TwoPi())));
   // 2-p
   fDirectCorrelations->Fill(0.,cos(n*(phi1-phi2)),wPhi1*wPhi2);                  // <w1   w2   cos( n*(phi1-phi2))>
   fDirectCorrelations->Fill(1.,cos(2.*n*(phi1-phi2)),pow(wPhi1,2)*pow(wPhi2,2)); // <w1^2 w2^2 cos(2n*(phi1-phi2))>
   fDirectCorrelations->Fill(2.,cos(3.*n*(phi1-phi2)),pow(wPhi1,3)*pow(wPhi2,3)); // <w1^3 w2^3 cos(3n*(phi1-phi2))>
   fDirectCorrelations->Fill(3.,cos(4.*n*(phi1-phi2)),pow(wPhi1,4)*pow(wPhi2,4)); // <w1^4 w2^4 cos(4n*(phi1-phi2))> 
   fDirectCorrelations->Fill(4.,cos(n*(phi1-phi2)),pow(wPhi1,3)*wPhi2); // <w1^3 w2 cos(n*(phi1-phi2))>
  }
 }  
 





 for(Int_t i1=0;i1<dMult;i1++)
 {
  fTrack=anEvent->GetTrack(i1);
  phi1=fTrack->Phi();
  if(phiWeights) wPhi1 = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi1*nBinsPhi/TMath::TwoPi())));
  for(Int_t i2=0;i2<dMult;i2++)
  {
   if(i2==i1)continue;
   fTrack=anEvent->GetTrack(i2);
   phi2=fTrack->Phi();
   if(phiWeights) wPhi2 = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi2*nBinsPhi/TMath::TwoPi())));
   for(Int_t i3=0;i3<dMult;i3++)
   {
    if(i3==i1||i3==i2)continue;
    fTrack=anEvent->GetTrack(i3);
    phi3=fTrack->Phi();
    if(phiWeights) wPhi3 = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi3*nBinsPhi/TMath::TwoPi())));
    // 2-p
    fDirectCorrelations->Fill(5.,cos(n*(phi1-phi2)),wPhi1*wPhi2*pow(wPhi3,2)); // <w1 w2 w3^2 cos(n*(phi1-phi2))>
    
    // 3-p
    fDirectCorrelations->Fill(10.,cos(2.*n*phi1-n*(phi2+phi3)),pow(wPhi1,2)*wPhi2*wPhi3); // <w1^2 w2 w3 cos(n*(2phi1-phi2-phi3))>
    
    
    fDirectCorrelations->Fill(6.,cos(3.*n*phi1-2.*n*phi2-n*phi3),pow(wPhi1,3)*pow(wPhi2,2)*wPhi3);           //<3>_{3n|2n,n}
    fDirectCorrelations->Fill(7.,cos(4.*n*phi1-2.*n*phi2-2.*n*phi3),pow(wPhi1,4)*pow(wPhi2,2)*pow(wPhi3,2)); //<3>_{4n|2n,2n}
    fDirectCorrelations->Fill(8.,cos(4.*n*phi1-3.*n*phi2-n*phi3),pow(wPhi1,4)*pow(wPhi2,3)*wPhi3);           //<3>_{4n|3n,n}
    
    
   }
  }
 }
 

 
 


  
 //<4>_{n,n|n,n}, <4>_{2n,n|2n,n}, <4>_{2n,2n|2n,2n}, <4>_{3n|n,n,n}, <4>_{3n,n|3n,n}, <4>_{3n,n|2n,2n} and <4>_{4n|2n,n,n} 
 for(Int_t i1=0;i1<dMult;i1++)
 {
  fTrack=anEvent->GetTrack(i1);
  phi1=fTrack->Phi();
  if(phiWeights) wPhi1 = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi1*nBinsPhi/TMath::TwoPi())));
  for(Int_t i2=0;i2<dMult;i2++)
  {
   if(i2==i1)continue;
   fTrack=anEvent->GetTrack(i2);
   phi2=fTrack->Phi();
   if(phiWeights) wPhi2 = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi2*nBinsPhi/TMath::TwoPi())));
   for(Int_t i3=0;i3<dMult;i3++)
   {
    if(i3==i1||i3==i2)continue;
    fTrack=anEvent->GetTrack(i3);
    phi3=fTrack->Phi();
    if(phiWeights) wPhi3 = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi3*nBinsPhi/TMath::TwoPi())));
    for(Int_t i4=0;i4<dMult;i4++)
    {
     if(i4==i1||i4==i2||i4==i3)continue;
     fTrack=anEvent->GetTrack(i4);
     phi4=fTrack->Phi();
     if(phiWeights) wPhi4 = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi4*nBinsPhi/TMath::TwoPi())));
     fDirectCorrelations->Fill(20.,cos(n*phi1+n*phi2-n*phi3-n*phi4),wPhi1*wPhi2*wPhi3*wPhi4); // <4>_{n,n|n,n} = <w1 w2 w3 w4 cos(n*(phi1+phi2-phi3-phi4))>
     
    
     fDirectCorrelations->Fill(11.,cos(2.*n*phi1+n*phi2-2.*n*phi3-n*phi4),wPhi1*wPhi2*wPhi3*wPhi4);      //<4>_{2n,n|2n,n}
     fDirectCorrelations->Fill(12.,cos(2.*n*phi1+2*n*phi2-2.*n*phi3-2.*n*phi4),wPhi1*wPhi2*wPhi3*wPhi4); //<4>_{2n,2n|2n,2n}
     fDirectCorrelations->Fill(13.,cos(3.*n*phi1-n*phi2-n*phi3-n*phi4),wPhi1*wPhi2*wPhi3*wPhi4);         //<4>_{3n|n,n,n}
     fDirectCorrelations->Fill(14.,cos(3.*n*phi1+n*phi2-3.*n*phi3-n*phi4),wPhi1*wPhi2*wPhi3*wPhi4);      //<4>_{3n,n|3n,n}   
     fDirectCorrelations->Fill(15.,cos(3.*n*phi1+n*phi2-2.*n*phi3-2.*n*phi4),wPhi1*wPhi2*wPhi3*wPhi4);   //<4>_{3n,n|2n,2n}
     fDirectCorrelations->Fill(16.,cos(4.*n*phi1-2.*n*phi2-n*phi3-n*phi4),wPhi1*wPhi2*wPhi3*wPhi4);      //<4>_{4n|2n,n,n}
    
    }  
   }
  }
 }
 
 
 

 
 //<5>_{2n,n,n,n,n}, //<5>_{2n,2n|2n,n,n}, <5>_{3n,n|2n,n,n} and <5>_{4n|n,n,n,n}
 for(Int_t i1=0;i1<dMult;i1++)
 {
  //cout<<"i1 = "<<i1<<endl;
  fTrack=anEvent->GetTrack(i1);
  phi1=fTrack->Phi();
  if(phiWeights) wPhi1 = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi1*nBinsPhi/TMath::TwoPi())));
  for(Int_t i2=0;i2<dMult;i2++)
  {
   if(i2==i1)continue;
   fTrack=anEvent->GetTrack(i2);
   phi2=fTrack->Phi();
   if(phiWeights) wPhi2 = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi2*nBinsPhi/TMath::TwoPi())));
   for(Int_t i3=0;i3<dMult;i3++)
   {
    if(i3==i1||i3==i2)continue;
    fTrack=anEvent->GetTrack(i3);
    phi3=fTrack->Phi();
    if(phiWeights) wPhi3 = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi3*nBinsPhi/TMath::TwoPi())));
    for(Int_t i4=0;i4<dMult;i4++)
    {
     if(i4==i1||i4==i2||i4==i3)continue;
     fTrack=anEvent->GetTrack(i4);
     phi4=fTrack->Phi();
     if(phiWeights) wPhi4 = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi4*nBinsPhi/TMath::TwoPi())));
     for(Int_t i5=0;i5<dMult;i5++)
     {
      if(i5==i1||i5==i2||i5==i3||i5==i4)continue;
      fTrack=anEvent->GetTrack(i5);
      phi5=fTrack->Phi();
      if(phiWeights) wPhi5 = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi5*nBinsPhi/TMath::TwoPi())));
      fDirectCorrelations->Fill(18.,cos(2.*n*phi1+n*phi2-n*phi3-n*phi4-n*phi5),wPhi1*wPhi2*wPhi3*wPhi4*wPhi5);       //<5>_{2n,n|n,n,n}
      fDirectCorrelations->Fill(19.,cos(2.*n*phi1+2.*n*phi2-2.*n*phi3-n*phi4-n*phi5),wPhi1*wPhi2*wPhi3*wPhi4*wPhi5); //<5>_{2n,2n|2n,n,n}
      fDirectCorrelations->Fill(20.,cos(3.*n*phi1+n*phi2-2.*n*phi3-n*phi4-n*phi5),wPhi1*wPhi2*wPhi3*wPhi4*wPhi5);    //<5>_{3n,n|2n,n,n}
      fDirectCorrelations->Fill(21.,cos(4.*n*phi1-n*phi2-n*phi3-n*phi4-n*phi5),wPhi1*wPhi2*wPhi3*wPhi4*wPhi5);       //<5>_{4n|n,n,n,n}
     }
    }  
   }
  }
 }
 
 
 //<6>_{n,n,n,n,n,n}, <6>_{2n,n,n|2n,n,n}, <6>_{2n,2n|n,n,n,n} and <6>_{3n,n|n,n,n,n}
 for(Int_t i1=0;i1<dMult;i1++)
 {
  //cout<<"i1 = "<<i1<<endl;
  fTrack=anEvent->GetTrack(i1);
  phi1=fTrack->Phi();
  if(phiWeights) wPhi1 = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi1*nBinsPhi/TMath::TwoPi())));
  for(Int_t i2=0;i2<dMult;i2++)
  {
   if(i2==i1)continue;
   fTrack=anEvent->GetTrack(i2);
   phi2=fTrack->Phi();
   if(phiWeights) wPhi2 = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi2*nBinsPhi/TMath::TwoPi())));
   for(Int_t i3=0;i3<dMult;i3++)
   {
    if(i3==i1||i3==i2)continue;
    fTrack=anEvent->GetTrack(i3);
    phi3=fTrack->Phi();
    if(phiWeights) wPhi3 = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi3*nBinsPhi/TMath::TwoPi())));
    for(Int_t i4=0;i4<dMult;i4++)
    {
     if(i4==i1||i4==i2||i4==i3)continue;
     fTrack=anEvent->GetTrack(i4);
     phi4=fTrack->Phi();
     if(phiWeights) wPhi4 = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi4*nBinsPhi/TMath::TwoPi())));
     for(Int_t i5=0;i5<dMult;i5++)
     {
      if(i5==i1||i5==i2||i5==i3||i5==i4)continue;
      fTrack=anEvent->GetTrack(i5);
      phi5=fTrack->Phi();
      if(phiWeights) wPhi5 = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi5*nBinsPhi/TMath::TwoPi())));
      for(Int_t i6=0;i6<dMult;i6++)
      {
       if(i6==i1||i6==i2||i6==i3||i6==i4||i6==i5)continue;
       fTrack=anEvent->GetTrack(i6);
       phi6=fTrack->Phi(); 
       if(phiWeights) wPhi6 = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi6*nBinsPhi/TMath::TwoPi())));
       fDirectCorrelations->Fill(23.,cos(n*phi1+n*phi2+n*phi3-n*phi4-n*phi5-n*phi6),wPhi1*wPhi2*wPhi3*wPhi4*wPhi5*wPhi6);       //<6>_{n,n,n|n,n,n}
       fDirectCorrelations->Fill(24.,cos(2.*n*phi1+n*phi2+n*phi3-2.*n*phi4-n*phi5-n*phi6),wPhi1*wPhi2*wPhi3*wPhi4*wPhi5*wPhi6); //<6>_{2n,n,n|2n,n,n}
       fDirectCorrelations->Fill(25.,cos(2.*n*phi1+2.*n*phi2-n*phi3-n*phi4-n*phi5-n*phi6),wPhi1*wPhi2*wPhi3*wPhi4*wPhi5*wPhi6); //<6>_{2n,2n|n,n,n,n}
       fDirectCorrelations->Fill(26.,cos(3.*n*phi1+n*phi2-n*phi3-n*phi4-n*phi5-n*phi6),wPhi1*wPhi2*wPhi3*wPhi4*wPhi5*wPhi6);    //<6>_{3n,n|n,n,n,n}     
      } 
     }
    }  
   }
  }
 }


 
 //<7>_{2n,n,n|n,n,n,n}
 for(Int_t i1=0;i1<dMult;i1++)
 {
  //cout<<"i1 = "<<i1<<endl;
  fTrack=anEvent->GetTrack(i1);
  phi1=fTrack->Phi();
  if(phiWeights) wPhi1 = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi1*nBinsPhi/TMath::TwoPi())));
  for(Int_t i2=0;i2<dMult;i2++)
  {
   if(i2==i1)continue;
   fTrack=anEvent->GetTrack(i2);
   phi2=fTrack->Phi();
   if(phiWeights) wPhi2 = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi2*nBinsPhi/TMath::TwoPi())));
   for(Int_t i3=0;i3<dMult;i3++)
   {
    if(i3==i1||i3==i2)continue;
    fTrack=anEvent->GetTrack(i3);
    phi3=fTrack->Phi();
    if(phiWeights) wPhi3 = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi3*nBinsPhi/TMath::TwoPi())));
    for(Int_t i4=0;i4<dMult;i4++)
    {
     if(i4==i1||i4==i2||i4==i3)continue;
     fTrack=anEvent->GetTrack(i4);
     phi4=fTrack->Phi();
     if(phiWeights) wPhi4 = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi4*nBinsPhi/TMath::TwoPi())));
     for(Int_t i5=0;i5<dMult;i5++)
     {
      if(i5==i1||i5==i2||i5==i3||i5==i4)continue;
      fTrack=anEvent->GetTrack(i5);
      phi5=fTrack->Phi();
      if(phiWeights) wPhi5 = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi5*nBinsPhi/TMath::TwoPi())));
      for(Int_t i6=0;i6<dMult;i6++)
      {
       if(i6==i1||i6==i2||i6==i3||i6==i4||i6==i5)continue;
       fTrack=anEvent->GetTrack(i6);
       phi6=fTrack->Phi(); 
       if(phiWeights) wPhi6 = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi6*nBinsPhi/TMath::TwoPi())));
       for(Int_t i7=0;i7<dMult;i7++)
       {
        if(i7==i1||i7==i2||i7==i3||i7==i4||i7==i5||i7==i6)continue;
        fTrack=anEvent->GetTrack(i7);
        phi7=fTrack->Phi(); 
        if(phiWeights) wPhi7 = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi7*nBinsPhi/TMath::TwoPi())));
        fDirectCorrelations->Fill(28.,cos(2.*n*phi1+n*phi2+n*phi3-n*phi4-n*phi5-n*phi6-n*phi7),wPhi1*wPhi2*wPhi3*wPhi4*wPhi5*wPhi6*wPhi7);//<7>_{2n,n,n|n,n,n,n}
       } 
      } 
     }
    }  
   }
  }
 }
 
 
 
 //<8>_{n,n,n,n|n,n,n,n}
 for(Int_t i1=0;i1<dMult;i1++)
 {
  cout<<"i1 = "<<i1<<endl;
  fTrack=anEvent->GetTrack(i1);
  phi1=fTrack->Phi();
  if(phiWeights) wPhi1 = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi1*nBinsPhi/TMath::TwoPi())));
  for(Int_t i2=0;i2<dMult;i2++)
  {
   if(i2==i1)continue;
   fTrack=anEvent->GetTrack(i2);
   phi2=fTrack->Phi();
   if(phiWeights) wPhi2 = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi2*nBinsPhi/TMath::TwoPi())));
   for(Int_t i3=0;i3<dMult;i3++)
   {
    if(i3==i1||i3==i2)continue;
    fTrack=anEvent->GetTrack(i3);
    phi3=fTrack->Phi();
    if(phiWeights) wPhi3 = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi3*nBinsPhi/TMath::TwoPi())));
    for(Int_t i4=0;i4<dMult;i4++)
    {
     if(i4==i1||i4==i2||i4==i3)continue;
     fTrack=anEvent->GetTrack(i4);
     phi4=fTrack->Phi();
     if(phiWeights) wPhi4 = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi4*nBinsPhi/TMath::TwoPi())));
     for(Int_t i5=0;i5<dMult;i5++)
     {
      if(i5==i1||i5==i2||i5==i3||i5==i4)continue;
      fTrack=anEvent->GetTrack(i5);
      phi5=fTrack->Phi();
      if(phiWeights) wPhi5 = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi5*nBinsPhi/TMath::TwoPi())));
      for(Int_t i6=0;i6<dMult;i6++)
      {
       if(i6==i1||i6==i2||i6==i3||i6==i4||i6==i5)continue;
       fTrack=anEvent->GetTrack(i6);
       phi6=fTrack->Phi();
       if(phiWeights) wPhi6 = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi6*nBinsPhi/TMath::TwoPi()))); 
       for(Int_t i7=0;i7<dMult;i7++)
       {
        if(i7==i1||i7==i2||i7==i3||i7==i4||i7==i5||i7==i6)continue;
        fTrack=anEvent->GetTrack(i7);
        phi7=fTrack->Phi();
        if(phiWeights) wPhi7 = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi7*nBinsPhi/TMath::TwoPi()))); 
        for(Int_t i8=0;i8<dMult;i8++)
        {
         if(i8==i1||i8==i2||i8==i3||i8==i4||i8==i5||i8==i6||i8==i7)continue;
         fTrack=anEvent->GetTrack(i8);
         phi8=fTrack->Phi();
         if(phiWeights) wPhi8 = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi8*nBinsPhi/TMath::TwoPi())));  
         fDirectCorrelations->Fill(30.,cos(n*phi1+n*phi2+n*phi3+n*phi4-n*phi5-n*phi6-n*phi7-n*phi8),wPhi1*wPhi2*wPhi3*wPhi4*wPhi5*wPhi6*wPhi7*wPhi8);//<8>_{n,n,n,n|n,n,n,n}
        } 
       } 
      } 
     }
    }  
   }
  }
 }
 
 */

 

 // binning details of fDirectCorrelations (differential flow):
 //
 //101st bin: <2'>_{n|n}
 
 
 /*
 
 //<2'>_{n|n}
 for(Int_t i1=0;i1<nPrim;i1++)
 {
  fTrack=anEvent->GetTrack(i1);
  if(!((fTrack->Pt()>=0.5&&fTrack->Pt()<0.6)&&(fTrack->UseForDifferentialFlow())))continue;//POI condition
  phi1=fTrack->Phi();
  if(phiWeights) wPhi1 = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi1*nBinsPhi/TMath::TwoPi())));
  for(Int_t i2=0;i2<nPrim;i2++)
  {
   if(i2==i1)continue;
   fTrack=anEvent->GetTrack(i2);
   if(!(fTrack->UseForIntegratedFlow()))continue;//RP condition
   phi2=fTrack->Phi();
   if(phiWeights) wPhi2 = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi2*nBinsPhi/TMath::TwoPi())));
    //cout<<"1st = "<<i1<<"     "<< (anEvent->GetTrack(i1))->Eta() << " " << (anEvent->GetTrack(i1))->Pt()<<endl;
    //cout<<"2nd = "<<i2<<"     "<< (anEvent->GetTrack(i2))->Eta() << " " << (anEvent->GetTrack(i2))->Pt()<<endl; 
   //fill the fDirectCorrelations:    
   fDirectCorrelations->Fill(100.,cos(1.*n*(phi1-phi2)),wPhi2);//<2'>_{n,n}
   fDirectCorrelations->Fill(103.,cos(1.*n*(phi1-phi2)),pow(wPhi1,2)*wPhi2);//<2'>_{n,n}
   fDirectCorrelations->Fill(104.,cos(2.*n*(phi1-phi2)),wPhi1*pow(wPhi2,2));//<2'>_{n,n}
   fDirectCorrelations->Fill(105.,cos(1.*n*(phi1-phi2)),pow(wPhi2,3));//<2'>_{n,n}
   
   
   
   fDirectCorrelations->Fill(41.,cos(2.*n*(phi1-phi2)),1);//<2'>_{2n,2n}
   fDirectCorrelations->Fill(42.,cos(3.*n*(phi1-phi2)),1);//<2'>_{3n,3n}
   fDirectCorrelations->Fill(43.,cos(4.*n*(phi1-phi2)),1);//<2'>_{4n,4n}   
    
  }//end of for(Int_t i2=0;i2<nPrim;i2++)
 }//end of for(Int_t i1=0;i1<nPrim;i1++)

 */ 
  
 
 
 /*
 
 //<3'>_{2n|n,n}
 for(Int_t i1=0;i1<nPrim;i1++)
 {
  fTrack=anEvent->GetTrack(i1);
  if(!((fTrack->Pt()>=0.5&&fTrack->Pt()<0.6)&&(fTrack->UseForDifferentialFlow())))continue;//POI condition
  phi1=fTrack->Phi();
  if(phiWeights) wPhi1 = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi1*nBinsPhi/TMath::TwoPi())));
  for(Int_t i2=0;i2<nPrim;i2++)
  {
   if(i2==i1)continue;
   fTrack=anEvent->GetTrack(i2);
   if(!(fTrack->UseForIntegratedFlow()))continue;//RP condition
   phi2=fTrack->Phi();
   if(phiWeights) wPhi2 = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi2*nBinsPhi/TMath::TwoPi())));
   for(Int_t i3=0;i3<nPrim;i3++)
   {
    if(i3==i1||i3==i2)continue;
    fTrack=anEvent->GetTrack(i3);
    if(!(fTrack->UseForIntegratedFlow()))continue;//RP condition
    phi3=fTrack->Phi();
    if(phiWeights) wPhi3 = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi3*nBinsPhi/TMath::TwoPi())));
    //fill the fDirectCorrelations:     
    
    // 2-p
    fDirectCorrelations->Fill(101.,cos(n*(phi2-phi3)),wPhi1*wPhi2*wPhi3); // <w1 w2 w3 cos(n(phi2-phi3))>
    fDirectCorrelations->Fill(102.,cos(n*(phi1-phi3)),pow(wPhi2,2.)*wPhi3); // <w2^2 w3 cos(n(psi1-phi2))>
    
    // 3-p            
    fDirectCorrelations->Fill(110.,cos(n*(2.*phi1-phi2-phi3)),wPhi1*wPhi2*wPhi3); // <w1 w2 w3 cos(n(2psi1-phi2-phi3))>
    fDirectCorrelations->Fill(111.,cos(n*(phi1+phi2-2.*phi3)),wPhi2*pow(wPhi3,2.)); // <w2 w3^2 cos(n(psi1+phi2-2.*phi3))>
    
    
    //fDirectCorrelations->Fill(46.,cos(n*(phi1+phi2-2.*phi3)),1);//<3'>_{n,n|2n}    
   }//end of for(Int_t i3=0;i3<nPrim;i3++)  
  }//end of for(Int_t i2=0;i2<nPrim;i2++)  
 }//end of for(Int_t i1=0;i1<nPrim;i1++) 
 
 
 */
 
 
 /*
 
 //<4'>_{n,n|n,n}
 for(Int_t i1=0;i1<nPrim;i1++)
 {
  fTrack=anEvent->GetTrack(i1);
  if(!((fTrack->Pt()>=0.5&&fTrack->Pt()<0.6)&&(fTrack->UseForDifferentialFlow())))continue;//POI condition
  tempCounter++;
  phi1=fTrack->Phi();
  if(phiWeights) wPhi1 = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi1*nBinsPhi/TMath::TwoPi())));
  for(Int_t i2=0;i2<nPrim;i2++)
  {
   if(i2==i1)continue;
   fTrack=anEvent->GetTrack(i2);
   if(!(fTrack->UseForIntegratedFlow()))continue;//RP condition   
   phi2=fTrack->Phi();
   if(phiWeights) wPhi2 = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi2*nBinsPhi/TMath::TwoPi())));
   for(Int_t i3=0;i3<nPrim;i3++)
   { 
    if(i3==i1||i3==i2)continue;
    fTrack=anEvent->GetTrack(i3);
    if(!(fTrack->UseForIntegratedFlow()))continue;//RP condition   
    phi3=fTrack->Phi();
    if(phiWeights) wPhi3 = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi3*nBinsPhi/TMath::TwoPi())));
    for(Int_t i4=0;i4<nPrim;i4++)
    {
     if(i4==i1||i4==i2||i4==i3)continue;
     fTrack=anEvent->GetTrack(i4);
     if(!(fTrack->UseForIntegratedFlow()))continue;//RP condition  
     phi4=fTrack->Phi();
     if(phiWeights) wPhi4 = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi4*nBinsPhi/TMath::TwoPi())));
     //fill the fDirectCorrelations:
     // 4-p            
     fDirectCorrelations->Fill(120.,cos(n*(phi1+phi2-phi3-phi4)),wPhi2*wPhi3*wPhi4); // <w1 w2 w3 cos(n(psi1+phi1-phi2-phi3))> 
     
       tempLoop+=wPhi2*wPhi3*wPhi4;
     
    }//end of for(Int_t i4=0;i4<nPrim;i4++)
   }//end of for(Int_t i3=0;i3<nPrim;i3++)
  }//end of for(Int_t i2=0;i2<nPrim;i2++) 
 }//end of for(Int_t i1=0;i1<nPrim;i1++)
  
 */
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 /*
 
 
 
 
 
 
 
 //<5'>_{2n,n|n,n,n}
 for(Int_t i1=0;i1<nPrim;i1++)
 {
  fTrack=anEvent->GetTrack(i1);
  if(!((fTrack->Pt()>=0.5&&fTrack->Pt()<0.6)&&(fTrack->UseForDifferentialFlow())))continue;//POI condition
  phi1=fTrack->Phi();
  for(Int_t i2=0;i2<nPrim;i2++)
  {
   if(i2==i1)continue;
   fTrack=anEvent->GetTrack(i2);
   if(!(fTrack->UseForIntegratedFlow()))continue;//RP condition   
   phi2=fTrack->Phi();
   for(Int_t i3=0;i3<nPrim;i3++)
   { 
    if(i3==i1||i3==i2)continue;
    fTrack=anEvent->GetTrack(i3);
    if(!(fTrack->UseForIntegratedFlow()))continue;//RP condition   
    phi3=fTrack->Phi();
    for(Int_t i4=0;i4<nPrim;i4++)
    {
     if(i4==i1||i4==i2||i4==i3)continue;
     fTrack=anEvent->GetTrack(i4);
     if(!(fTrack->UseForIntegratedFlow()))continue;//RP condition  
     phi4=fTrack->Phi();//
     for(Int_t i5=0;i5<nPrim;i5++)
     {
      if(i5==i1||i5==i2||i5==i3||i5==i4)continue;
      fTrack=anEvent->GetTrack(i5);
      if(!(fTrack->UseForIntegratedFlow()))continue;//RP condition  
      phi5=fTrack->Phi();    
      //fill the fDirectCorrelations:if(bNestedLoops)
      fDirectCorrelations->Fill(55.,cos(2.*n*phi1+n*phi2-n*phi3-n*phi4-n*phi5),1);//<5'>_{2n,n|n,n,n}
     }//end of for(Int_t i5=0;i5<nPrim;i5++)  
    }//end of for(Int_t i4=0;i4<nPrim;i4++)
   }//end of for(Int_t i3=0;i3<nPrim;i3++)
  }//end of for(Int_t i2=0;i2<nPrim;i2++) 
 }//end of for(Int_t i1=0;i1<nPrim;i1++)
 
 //<6'>_{n,n,n|n,n,n}
 for(Int_t i1=0;i1<nPrim;i1++)
 {
  fTrack=anEvent->GetTrack(i1);
  if(!((fTrack->Pt()>=0.5&&fTrack->Pt()<0.6)&&(fTrack->UseForDifferentialFlow())))continue;//POI condition
  phi1=fTrack->Phi();
  for(Int_t i2=0;i2<nPrim;i2++)
  {
   if(i2==i1)continue;
   fTrack=anEvent->GetTrack(i2);
   if(!(fTrack->UseForIntegratedFlow()))continue;//RP condition   
   phi2=fTrack->Phi();
   for(Int_t i3=0;i3<nPrim;i3++)
   { 
    if(i3==i1||i3==i2)continue;
    fTrack=anEvent->GetTrack(i3);
    if(!(fTrack->UseForIntegratedFlow()))continue;//RP condition   
    phi3=fTrack->Phi();
    for(Int_t i4=0;i4<nPrim;i4++)
    {
     if(i4==i1||i4==i2||i4==i3)continue;
     fTrack=anEvent->GetTrack(i4);
     if(!(fTrack->UseForIntegratedFlow()))continue;//RP condition  
     phi4=fTrack->Phi();
     for(Int_t i5=0;i5<nPrim;i5++)
     {
      if(i5==i1||i5==i2||i5==i3||i5==i4)continue;
      fTrack=anEvent->GetTrack(i5);
      if(!(fTrack->UseForIntegratedFlow()))continue;//RP condition  
      phi5=fTrack->Phi();    
      for(Int_t i6=0;i6<nPrim;i6++)
      {
       if(i6==i1||i6==i2||i6==i3||i6==i4||i6==i5)continue;
       fTrack=anEvent->GetTrack(i6);
       if(!(fTrack->UseForIntegratedFlow()))continue;//RP condition  
       phi6=fTrack->Phi();  
       //fill the fDirectCorrelations:
       fDirectCorrelations->Fill(60.,cos(n*(phi1+phi2+phi3-phi4-phi5-phi6)),1);//<6'>_{n,n,n|n,n,n}
      }//end of for(Int_t i6=0;i6<nPrim;i6++)   
     }//end of for(Int_t i5=0;i5<nPrim;i5++)  
    }//end of for(Int_t i4=0;i4<nPrim;i4++)
   }//end of for(Int_t i3=0;i3<nPrim;i3++)
  }//end of for(Int_t i2=0;i2<nPrim;i2++) 
 }//end of for(Int_t i1=0;i1<nPrim;i1++)
 
 //<7'>_{2n,n,n|n,n,n,n}
 for(Int_t i1=0;i1<nPrim;i1++)
 {
  fTrack=anEvent->GetTrack(i1);
  if(!((fTrack->Pt()>=0.5&&fTrack->Pt()<0.6)&&(fTrack->UseForDifferentialFlow())))continue;//POI condition
  phi1=fTrack->Phi();
  for(Int_t i2=0;i2<nPrim;i2++)
  {
   if(i2==i1)continue;
   fTrack=anEvent->GetTrack(i2);
   if(!(fTrack->UseForIntegratedFlow()))continue;//RP condition   
   phi2=fTrack->Phi();
   for(Int_t i3=0;i3<nPrim;i3++)
   { 
    if(i3==i1||i3==i2)continue;
    fTrack=anEvent->GetTrack(i3);
    if(!(fTrack->UseForIntegratedFlow()))continue;//RP condition   
    phi3=fTrack->Phi();
    for(Int_t i4=0;i4<nPrim;i4++)
    {
     if(i4==i1||i4==i2||i4==i3)continue;
     fTrack=anEvent->GetTrack(i4);
     if(!(fTrack->UseForIntegratedFlow()))continue;//RP condition  
     phi4=fTrack->Phi();
     for(Int_t i5=0;i5<nPrim;i5++)
     {
      if(i5==i1||i5==i2||i5==i3||i5==i4)continue;
      fTrack=anEvent->GetTrack(i5);
      if(!(fTrack->UseForIntegratedFlow()))continue;//RP condition  
      phi5=fTrack->Phi();    
      for(Int_t i6=0;i6<nPrim;i6++)
      {
       if(i6==i1||i6==i2||i6==i3||i6==i4||i6==i5)continue;
       fTrack=anEvent->GetTrack(i6);
       if(!(fTrack->UseForIntegratedFlow()))continue;//RP condition  
       phi6=fTrack->Phi();
       for(Int_t i7=0;i7<nPrim;i7++)
       {
        if(i7==i1||i7==i2||i7==i3||i7==i4||i7==i5||i7==i6)continue;
        fTrack=anEvent->GetTrack(i7);
        if(!(fTrack->UseForIntegratedFlow()))continue;//RP condition  
        phi7=fTrack->Phi();   
        //fill the fDirectCorrelations:
        fDirectCorrelations->Fill(65.,cos(2.*n*phi1+n*phi2+n*phi3-n*phi4-n*phi5-n*phi6-n*phi7),1);//<7'>_{2n,n,n|n,n,n,n}
       }//end of for(Int_t i7=0;i7<nPrim;i7++)  
      }//end of for(Int_t i6=0;i6<nPrim;i6++)   
     }//end of for(Int_t i5=0;i5<nPrim;i5++)  
    }//end of for(Int_t i4=0;i4<nPrim;i4++)
   }//end of for(Int_t i3=0;i3<nPrim;i3++)
  }//end of for(Int_t i2=0;i2<nPrim;i2++) 
 }//end of for(Int_t i1=0;i1<nPrim;i1++)
 
 //<8'>_{n,n,n,n|n,n,n,n}
 for(Int_t i1=0;i1<nPrim;i1++)
 {
  fTrack=anEvent->GetTrack(i1);
  if(!((fTrack->Pt()>=0.5&&fTrack->Pt()<0.6)&&(fTrack->UseForDifferentialFlow())))continue;//POI condition
  phi1=fTrack->Phi();
  for(Int_t i2=0;i2<nPrim;i2++)
  {
   if(i2==i1)continue;
   fTrack=anEvent->GetTrack(i2);
   if(!(fTrack->UseForIntegratedFlow()))continue;//RP condition   
   phi2=fTrack->Phi();
   for(Int_t i3=0;i3<nPrim;i3++)
   { 
    if(i3==i1||i3==i2)continue;
    fTrack=anEvent->GetTrack(i3);
    if(!(fTrack->UseForIntegratedFlow()))continue;//RP condition   
    phi3=fTrack->Phi();
    for(Int_t i4=0;i4<nPrim;i4++)
    {
     if(i4==i1||i4==i2||i4==i3)continue;
     fTrack=anEvent->GetTrack(i4);
     if(!(fTrack->UseForIntegratedFlow()))continue;//RP condition  
     phi4=fTrack->Phi();
     for(Int_t i5=0;i5<nPrim;i5++)
     {
      if(i5==i1||i5==i2||i5==i3||i5==i4)continue;
      fTrack=anEvent->GetTrack(i5);
      if(!(fTrack->UseForIntegratedFlow()))continue;//RP condition  
      phi5=fTrack->Phi();    
      for(Int_t i6=0;i6<nPrim;i6++)
      {
       if(i6==i1||i6==i2||i6==i3||i6==i4||i6==i5)continue;
       fTrack=anEvent->GetTrack(i6);
       if(!(fTrack->UseForIntegratedFlow()))continue;//RP condition  
       phi6=fTrack->Phi();
       for(Int_t i7=0;i7<nPrim;i7++)
       {
        if(i7==i1||i7==i2||i7==i3||i7==i4||i7==i5||i7==i6)continue;
        fTrack=anEvent->GetTrack(i7);
        if(!(fTrack->UseForIntegratedFlow()))continue;//RP condition  
        phi7=fTrack->Phi();
        for(Int_t i8=0;i8<nPrim;i8++)
        {
         if(i8==i1||i8==i2||i8==i3||i8==i4||i8==i5||i8==i6||i8==i7)continue;
         fTrack=anEvent->GetTrack(i8);
         if(!(fTrack->UseForIntegratedFlow()))continue;//RP condition  
         phi8=fTrack->Phi();           
         //fill the fDirectCorrelations:
         fDirectCorrelations->Fill(70.,cos(n*(phi1+phi2+phi3+phi4-phi5-phi6-phi7-phi8)),1);//<8'>_{n,n,n,n|n,n,n,n}
        }//end of for(Int_t i8=0;i8<nPrim;i8++) 
       }//end of for(Int_t i7=0;i7<nPrim;i7++)  
      }//end of for(Int_t i6=0;i6<nPrim;i6++)   
     }//end of for(Int_t i5=0;i5<nPrim;i5++)  
    }//end of for(Int_t i4=0;i4<nPrim;i4++)
   }//end of for(Int_t i3=0;i3<nPrim;i3++)
  }//end of for(Int_t i2=0;i2<nPrim;i2++) 
 }//end of for(Int_t i1=0;i1<nPrim;i1++)
 
 
 
 
 
 
 
 
 
 
 */
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 

 
 
 
 
 
 
 
 
 
//--------------------------------------------------------------------------------------------------------------------------------  

 //}//end of if(nPrim>0&&nPrim<12)
}//end of Make()

//================================================================================================================

void AliFlowAnalysisWithQCumulants::Finish()
{
 //calculate the final results
 
 // harmonics
 Int_t n = 2; // to be improved 
 
 //--------------------------------------------------------------------------------------------------------- 
 // avarage multiplicity
 Double_t AvMPOI = (fCommonHists2nd->GetHistMultDiff())->GetMean(); // to be improved 
 Double_t AvMRP  = (fCommonHists2nd->GetHistMultInt())->GetMean(); // to be improved 

 // number of events
 Double_t nEvtsPOI = (fCommonHists2nd->GetHistMultDiff())->GetEntries(); // to be improved
 Double_t nEvtsRP  = (fCommonHists2nd->GetHistMultInt())->GetEntries(); // to be improved
 //---------------------------------------------------------------------------------------------------------
 
 //---------------------------------------------------------------------------------------------------------
 // 2-, 4-, 6- and 8-particle azimuthal correlation:
 Double_t two   = fQCorrelations->GetBinContent(1);  //<<2>>_{n|n}
 Double_t four  = fQCorrelations->GetBinContent(11); //<<4>>_{n,n|n,n}
 Double_t six   = fQCorrelations->GetBinContent(24); //<<6>>_{n,n,n|n,n,n}
 Double_t eight = fQCorrelations->GetBinContent(31); //<<8>>_{n,n,n,n|n,n,n,n}
 
 // 2nd, 4th, 6th and 8th order Q-cumulant:
 Double_t secondOrderQCumulant = two; //c_n{2} 
 Double_t fourthOrderQCumulant = four-2.*pow(two,2.); //c_n{4}
 Double_t sixthOrderQCumulant  = six-9.*two*four+12.*pow(two,3.); //c_n{6}
 Double_t eightOrderQCumulant  = eight-16.*two*six-18.*pow(four,2.)+144.*pow(two,2.)*four-144.*pow(two,4.); //c_n{8} 
 
 // "no-name" integrated flow estimates from Q-cumulants:
 Double_t vn2=0.,vn4=0.,vn6=0.,vn8=0.;
 // Double_t sd2=0.,sd4=0.,sd6=0.,sd8=0.; 
 Double_t sd6=0.,sd8=0.; // to be improved/removed
 if(secondOrderQCumulant>0.)
 {
  vn2 = pow(secondOrderQCumulant,0.5); //v_n{2}
 } 
 if(fourthOrderQCumulant<0.)
 {
  vn4 = pow(-fourthOrderQCumulant,1./4.); //v_n{4}
 } 
 if(sixthOrderQCumulant>0.)
 {
  vn6 = pow((1./4.)*sixthOrderQCumulant,1./6.); //v_n{6}
 } 
 if(eightOrderQCumulant<0.)
 {
  vn8 = pow((-1./33.)*eightOrderQCumulant,1./8.); //v_n{8}
 }
 //---------------------------------------------------------------------------------------------------------

 //--------------------------------------------------------------------------------------------------------- 
 // weighted 2-, 4-, 6- and 8-particle azimuthal correlation:
 Double_t twoW   = fWeightedQCorrelations->GetBinContent(1);  //<<2>>_{n|n}
 Double_t fourW  = fWeightedQCorrelations->GetBinContent(21); //<<4>>_{n,n|n,n}
 // Double_t sixW   = fWeightedQCorrelations->GetBinContent(24); //<<6>>_{n,n,n|n,n,n}
 // Double_t eightW = fWeightedQCorrelations->GetBinContent(31); //<<8>>_{n,n,n,n|n,n,n,n}
 
 // 2nd, 4th, 6th and 8th order weighted Q-cumulant:
 Double_t secondOrderQCumulantW = twoW; //c_n{2} 
 Double_t fourthOrderQCumulantW = fourW-2.*pow(twoW,2.); //c_n{4}
 // Double_t sixthOrderQCumulantW  = sixW-9.*twoW*fourW+12.*pow(twoW,3.); //c_n{6}
 // Double_t eightOrderQCumulantW  = eightW-16.*twoW*sixW-18.*pow(fourW,2.)+144.*pow(twoW,2.)*fourW-144.*pow(twoW,4.); //c_n{8} 
 
 // "no-name" integrated flow estimates from weighted Q-cumulants:
 cout<<endl;
 cout<<"**************************************"<<endl;
 cout<<"**************************************"<<endl;
 cout<<"flow estimates from Q-cumulants :"<<endl;
 cout<<endl;
 
 Double_t vn2W=0.,vn4W=0.;
 Double_t sd2W=0.,sd4W=0.; 
 if(secondOrderQCumulantW>0.){
  vn2W = pow(secondOrderQCumulantW,0.5); // weighted v_n{2}
  //sd2W = 0.5*pow(secondOrderQCumulantW,-0.5)*secondOrderQCumulantErrorW; // to be improved (correct treatment of errors needed)
  cout<<" v_"<<n<<"{2} = "<<vn2W<<" +/- "<<sd2W<<endl;
  fIntFlowResultsQC->SetBinContent(1,vn2W);
  //fIntFlowResultsQC->SetBinError(1,sd2W);
  //common histograms:
  fCommonHistsResults2nd->FillIntegratedFlow(vn2W,sd2W);
  //fCommonHistsResults2nd->FillChi(vn2W*pow(AvM,0.5)); // to be removed
 }else{
  cout<<" v_"<<n<<"{2} = Im"<<endl;
 }          
 if(fourW!=0. && fourthOrderQCumulantW<0.){
  vn4W = pow(-fourthOrderQCumulantW,1./4.); //v_n{4}
  //sd4W = 0.25*pow(-fourthOrderQCumulantW,-3./4.)*fourthOrderQCumulantErrorW; // to be improved (correct treatment of errors needed)
  cout<<" v_"<<n<<"{4} = "<<vn4W<<" +/- "<<sd4W<<endl;
  fIntFlowResultsQC->SetBinContent(2,vn4W);
  //fIntFlowResultsQC->SetBinError(2,sd4W);
  //common histograms:
  fCommonHistsResults4th->FillIntegratedFlow(vn4W,sd4W);
  //fCommonHistsResults4th->FillChi(vn4W*pow(AvM,0.5)); // to be removed
 }else{
  cout<<" v_"<<n<<"{4} = Im"<<endl;
 }
 // !!! to be improved (6th and 8th order are without weights) !!!
 if(six!=0. && sixthOrderQCumulant>0.){
  vn6 = pow((1./4.)*sixthOrderQCumulant,1./6.); //v_n{6}
  //sd6 = (1./6.)*pow(2.,-1./3.)*pow(sixthOrderQCumulant,-5./6.)*sixthOrderQCumulantError;
  cout<<" v_"<<n<<"{6} = "<<vn6<<" +/- "<<sd6<<endl;
  fIntFlowResultsQC->SetBinContent(3,vn6);
  //fIntFlowResultsQC->SetBinError(3,sd6);
  //common histograms:
  fCommonHistsResults6th->FillIntegratedFlow(vn6,sd6);
  //fCommonHistsResults6th->FillChi(vn6*pow(AvM,0.5));//to be removed
 }else{
  cout<<" v_"<<n<<"{6} = Im"<<endl;
 }
 if(eight!=0. && eightOrderQCumulant<0.){
  vn8 = pow((-1./33.)*eightOrderQCumulant,1./8.); //v_n{8}
  cout<<" v_"<<n<<"{8} = "<<vn8<<" +/- "<<sd8<<endl;
  fIntFlowResultsQC->SetBinContent(4,vn8);
  //fIntFlowResultsQC->SetBinError(4,sd8);
  //common histograms:
  fCommonHistsResults8th->FillIntegratedFlow(vn8,sd8);
  //fCommonHistsResults8th->FillChi(vn8*pow(AvM,0.5));//to be removed
 }else{
  cout<<" v_"<<n<<"{8} = Im"<<endl;
 }
 cout<<endl;
 cout<<"   nEvts = "<<nEvtsRP<<", AvM = "<<AvMRP<<endl; // to be improved
 cout<<"**************************************"<<endl;
 cout<<"**************************************"<<endl;
 cout<<endl; 
//--------------------------------------------------------------------------------------------------------- 
 
//---------------------------------------------------------------------------------------------------------
// differential flow (POI)
Int_t nBinsPtPOI = f2WPerPtBin1n1nPOI->GetNbinsX();
Int_t nBinsEtaPOI = f4WPerPtBin1n1n1n1nPOI->GetNbinsX();

// Pt:
Double_t secondOrderQCumulantDiffFlowPtPOIW = 0.;
Double_t fourthOrderQCumulantDiffFlowPtPOIW = 0.;

Double_t dVn2ndPOIW=0.,dSd2ndPOIW=0.,dDiffvn2ndPOIW=0.,dYield2ndPOIW=0.,dSum2ndPOIW=0.;
Double_t dVn4thPOIW=0.,dSd4thPOIW=0.,dDiffvn4thPOIW=0.,dYield4thPOIW=0.,dSum4thPOIW=0.;

for(Int_t bb=1;bb<nBinsPtPOI+1;bb++)
{
 // QC{2}
 if(f2WPerPtBin1n1nPOI->GetBinEntries(bb)>0.&&vn2W!=0)
 { 
  secondOrderQCumulantDiffFlowPtPOIW = f2WPerPtBin1n1nPOI->GetBinContent(bb); // with weights
  fDiffFlowResults2ndOrderQC->SetBinContent(bb,secondOrderQCumulantDiffFlowPtPOIW/vn2W);
  // common histogram:
  fCommonHistsResults2nd->FillDifferentialFlowPtPOI(bb,secondOrderQCumulantDiffFlowPtPOIW/vn2W, 0.); //to be improved (errors && bb or bb+1 ?)
  // -------------------------------------------------------------------
  // integrated flow (weighted, POI, Pt, 2nd order):
  dDiffvn2ndPOIW=(fCommonHistsResults2nd->GetHistDiffFlowPtPOI())->GetBinContent(bb);
  dYield2ndPOIW=(fCommonHists2nd->GetHistPtDiff())->GetBinContent(bb);
  dVn2ndPOIW+=dDiffvn2ndPOIW*dYield2ndPOIW;
  dSum2ndPOIW+=dYield2ndPOIW;
  // -------------------------------------------------------------------
 }
 // QC{4]
 if(f4WPerPtBin1n1n1n1nPOI->GetBinEntries(bb)>0.&&vn4W!=0.)
 {
  fourthOrderQCumulantDiffFlowPtPOIW = f4WPerPtBin1n1n1n1nPOI->GetBinContent(bb)-2.*f2WPerPtBin1n1nPOI->GetBinContent(bb)*pow(vn2W,2.); // with weights
  fDiffFlowResults4thOrderQC->SetBinContent(bb,-1.*fourthOrderQCumulantDiffFlowPtPOIW/pow(vn4W,3.));
  //common histogram:
  fCommonHistsResults4th->FillDifferentialFlowPtPOI(bb,-1.*fourthOrderQCumulantDiffFlowPtPOIW/pow(vn4W,3.), 0.); //to be improved (errors)
  // -------------------------------------------------------------------
  //integrated flow (POI, Pt, 4th order):
  dDiffvn4thPOIW=(fCommonHistsResults4th->GetHistDiffFlowPtPOI())->GetBinContent(bb);
  dYield4thPOIW=(fCommonHists4th->GetHistPtDiff())->GetBinContent(bb);
  dVn4thPOIW+=dDiffvn4thPOIW*dYield4thPOIW;
  dSum4thPOIW+=dYield4thPOIW;
  // -------------------------------------------------------------------
 }
}      

cout<<endl;
cout<<"**************************************"<<endl;
cout<<"**************************************"<<endl;
cout<<"flow estimates from Q-cumulants (POI):"<<endl;
cout<<endl;
//storing the final results for integrated flow (POI):
// QC{2}
if(dSum2ndPOIW && fCommonHistsResults2nd)
{
 dVn2ndPOIW/=dSum2ndPOIW;
 fCommonHistsResults2nd->FillIntegratedFlowPOI(dVn2ndPOIW,0.); // to be improved (errors)
 cout<<" v_"<<n<<"{2} = "<<dVn2ndPOIW<<" +/- "<<dSd2ndPOIW<<endl;
}else 
 {
  cout<<" v_"<<n<<"{2} = Im"<<endl;
 }

// QC{4}
if(dSum4thPOIW && fCommonHistsResults4th)
{
 dVn4thPOIW/=dSum4thPOIW;
 fCommonHistsResults4th->FillIntegratedFlowPOI(dVn4thPOIW,0.); // to be improved (errors)
 cout<<" v_"<<n<<"{4} = "<<dVn4thPOIW<<" +/- "<<dSd4thPOIW<<endl;
}else
 {
  cout<<" v_"<<n<<"{4} = Im"<<endl;
 }

cout<<endl;
cout<<"   nEvts = "<<nEvtsPOI<<", AvM = "<<AvMPOI<<endl;
cout<<"**************************************"<<endl;
cout<<"**************************************"<<endl;
cout<<endl;  

//Eta:
Double_t secondOrderQCumulantDiffFlowEtaPOIW = 0.;
Double_t fourthOrderQCumulantDiffFlowEtaPOIW = 0.;

for(Int_t bb=1;bb<nBinsEtaPOI+1;bb++)
{
 if(f2WPerEtaBin1n1nPOI->GetBinEntries(bb)>0.&&vn2W!=0)
 {
  secondOrderQCumulantDiffFlowEtaPOIW = f2WPerEtaBin1n1nPOI->GetBinContent(bb); // with weights
  fDiffFlowResults2ndOrderQC->SetBinContent(bb,secondOrderQCumulantDiffFlowEtaPOIW/vn2W);
  //common histogram:
  fCommonHistsResults2nd->FillDifferentialFlowEtaPOI(bb,secondOrderQCumulantDiffFlowEtaPOIW/vn2W, 0.);//to be improved (errors)
 }
 if(f4WPerEtaBin1n1n1n1nPOI->GetBinEntries(bb)>0.&&vn4W!=0.)
 {
  fourthOrderQCumulantDiffFlowEtaPOIW = f4WPerEtaBin1n1n1n1nPOI->GetBinContent(bb)-2.*f2WPerEtaBin1n1nPOI->GetBinContent(bb)*pow(vn2W,2.); // with weights
  fDiffFlowResults4thOrderQC->SetBinContent(bb,-1.*fourthOrderQCumulantDiffFlowEtaPOIW/pow(vn4W,3.));
  //common histogram:
  fCommonHistsResults4th->FillDifferentialFlowEtaPOI(bb,-1.*fourthOrderQCumulantDiffFlowEtaPOIW/pow(vn4W,3.), 0.);//to be improved (errors)
 }
}      
//------------------------------------------------------------
 
 
 
 
 
 //---------------------------------------------------------------------------------------------------------
// differential flow (RP)
Int_t nBinsPtRP = f2WPerPtBin1n1nRP->GetNbinsX();
Int_t nBinsEtaRP = f4WPerPtBin1n1n1n1nRP->GetNbinsX();

// Pt:
Double_t secondOrderQCumulantDiffFlowPtRPW = 0.;
Double_t fourthOrderQCumulantDiffFlowPtRPW = 0.;

Double_t dVn2ndRPW=0.,dSd2ndRPW=0.,dDiffvn2ndRPW=0.,dYield2ndRPW=0.,dSum2ndRPW=0.;
Double_t dVn4thRPW=0.,dSd4thRPW=0.,dDiffvn4thRPW=0.,dYield4thRPW=0.,dSum4thRPW=0.;

for(Int_t bb=1;bb<nBinsPtRP+1;bb++)
{
 // QC{2}
 if(f2WPerPtBin1n1nRP->GetBinEntries(bb)>0.&&vn2W!=0)
 { 
  secondOrderQCumulantDiffFlowPtRPW = f2WPerPtBin1n1nRP->GetBinContent(bb); // with weights
  fDiffFlowResults2ndOrderQC->SetBinContent(bb,secondOrderQCumulantDiffFlowPtRPW/vn2W);
  // common histogram:
  fCommonHistsResults2nd->FillDifferentialFlowPtRP(bb,secondOrderQCumulantDiffFlowPtRPW/vn2W, 0.); //to be improved (errors && bb or bb+1 ?)
  // -------------------------------------------------------------------
  // integrated flow (weighted, RP, Pt, 2nd order):
  dDiffvn2ndRPW=(fCommonHistsResults2nd->GetHistDiffFlowPtRP())->GetBinContent(bb);
  dYield2ndRPW=(fCommonHists2nd->GetHistPtInt())->GetBinContent(bb);
  dVn2ndRPW+=dDiffvn2ndRPW*dYield2ndRPW;
  dSum2ndRPW+=dYield2ndRPW;
  // -------------------------------------------------------------------
 }
 // QC{4]
 if(f4WPerPtBin1n1n1n1nRP->GetBinEntries(bb)>0.&&vn4W!=0.)
 {
  fourthOrderQCumulantDiffFlowPtRPW = f4WPerPtBin1n1n1n1nRP->GetBinContent(bb)-2.*f2WPerPtBin1n1nRP->GetBinContent(bb)*pow(vn2W,2.); // with weights
  fDiffFlowResults4thOrderQC->SetBinContent(bb,-1.*fourthOrderQCumulantDiffFlowPtRPW/pow(vn4W,3.));
  //common histogram:
  fCommonHistsResults4th->FillDifferentialFlowPtRP(bb,-1.*fourthOrderQCumulantDiffFlowPtRPW/pow(vn4W,3.), 0.); //to be improved (errors)
  // -------------------------------------------------------------------
  //integrated flow (RP, Pt, 4th order):
  dDiffvn4thRPW=(fCommonHistsResults4th->GetHistDiffFlowPtRP())->GetBinContent(bb);
  dYield4thRPW=(fCommonHists4th->GetHistPtInt())->GetBinContent(bb);
  dVn4thRPW+=dDiffvn4thRPW*dYield4thRPW;
  dSum4thRPW+=dYield4thRPW;
  // -------------------------------------------------------------------
 }
}      

cout<<endl;
cout<<"**************************************"<<endl;
cout<<"**************************************"<<endl;
cout<<"flow estimates from Q-cumulants (RP):"<<endl;
cout<<endl;
//storing the final results for integrated flow (RP):
// QC{2}
if(dSum2ndRPW && fCommonHistsResults2nd)
{
 dVn2ndRPW/=dSum2ndRPW;
 fCommonHistsResults2nd->FillIntegratedFlowRP(dVn2ndRPW,0.); // to be improved (errors)
 cout<<" v_"<<n<<"{2} = "<<dVn2ndRPW<<" +/- "<<dSd2ndRPW<<endl;
}else 
 {
  cout<<" v_"<<n<<"{2} = Im"<<endl;
 }

// QC{4}
if(dSum4thRPW && fCommonHistsResults4th)
{
 dVn4thRPW/=dSum4thRPW;
 fCommonHistsResults4th->FillIntegratedFlowRP(dVn4thRPW,0.); // to be improved (errors)
 cout<<" v_"<<n<<"{4} = "<<dVn4thRPW<<" +/- "<<dSd4thRPW<<endl;
}else
 {
  cout<<" v_"<<n<<"{4} = Im"<<endl;
 }

cout<<endl;
cout<<"   nEvts = "<<nEvtsRP<<", AvM = "<<AvMRP<<endl;
cout<<"**************************************"<<endl;
cout<<"**************************************"<<endl;
cout<<endl;  

//Eta:
Double_t secondOrderQCumulantDiffFlowEtaRPW = 0.;
Double_t fourthOrderQCumulantDiffFlowEtaRPW = 0.;

for(Int_t bb=1;bb<nBinsEtaRP+1;bb++)
{
 if(f2WPerEtaBin1n1nRP->GetBinEntries(bb)>0.&&vn2W!=0)
 {
  secondOrderQCumulantDiffFlowEtaRPW = f2WPerEtaBin1n1nRP->GetBinContent(bb); // with weights
  fDiffFlowResults2ndOrderQC->SetBinContent(bb,secondOrderQCumulantDiffFlowEtaRPW/vn2W);
  //common histogram:
  fCommonHistsResults2nd->FillDifferentialFlowEtaRP(bb,secondOrderQCumulantDiffFlowEtaRPW/vn2W, 0.);//to be improved (errors)
 }
 if(f4WPerEtaBin1n1n1n1nRP->GetBinEntries(bb)>0.&&vn4W!=0.)
 {
  fourthOrderQCumulantDiffFlowEtaRPW = f4WPerEtaBin1n1n1n1nRP->GetBinContent(bb)-2.*f2WPerEtaBin1n1nRP->GetBinContent(bb)*pow(vn2W,2.); // with weights
  fDiffFlowResults4thOrderQC->SetBinContent(bb,-1.*fourthOrderQCumulantDiffFlowEtaRPW/pow(vn4W,3.));
  //common histogram:
  fCommonHistsResults4th->FillDifferentialFlowEtaRP(bb,-1.*fourthOrderQCumulantDiffFlowEtaRPW/pow(vn4W,3.), 0.);//to be improved (errors)
 }
}      
//------------------------------------------------------------
 
 
 
 /*
 cout<<"0.5 < Pt < 0.6 GeV"<<endl;                                
 cout<<endl;                                       
 cout<<"<w2 cos(n(psi1-phi2))> from Q-vectors                    = "<<f2WPerPtBin1n1nPOI->GetBinContent(6)<<endl;
 //cout<<"<w2 cos(n(psi1-phi2))> from Q-vectors                    = "<<f2WPerPtBin1n1nRP->GetName()<<endl;
 cout<<"<w2 cos(n(psi1-phi2))> from Q-vectors                    = "<<fDirectCorrelations->GetBinContent(101)<<endl;
 cout<<endl;  
 cout<<"<w2 w3 w4 cos(n(psi1+phi2-phi3-phi4))> from Q-vectors    = "<<f4WPerPtBin1n1n1n1nPOI->GetBinContent(6)<<endl;
 //cout<<"<w2 w3 w4 cos(n(psi1+phi2-phi3-phi4))> from Q-vectors    = "<<f4WPerPtBin1n1n1n1nRP->GetName()<<endl;
 cout<<"<w2 w3 w4 cos(n(psi1+phi2-phi3-phi4))> from nested loops = "<<fDirectCorrelations->GetBinContent(121)<<endl;   
 cout<<endl; 
 
 
 
 
 
  */
 
 
 
 /*
   
    
 //xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 //                   !!!! to be removed !!!!  
 
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
                           
}

//================================================================================================================

void AliFlowAnalysisWithQCumulants::WriteHistograms(TString outputFileName)
{
 //store the final results in output .root file
 TFile *output = new TFile(outputFileName.Data(),"RECREATE");
 output->WriteObject(fHistList, "cobjQC","SingleKey");
 delete output;
}

//================================================================================================================



