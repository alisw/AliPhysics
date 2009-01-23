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
#include "TParticle.h"
#include "TRandom3.h"
#include "TStyle.h"
#include "TProfile.h"
#include "TProfile2D.h" 
#include "TProfile3D.h"
#include "TMath.h"
#include "AliFlowEventSimple.h"
#include "AliFlowTrackSimple.h"
#include "AliFlowAnalysisWithQCumulants.h"
#include "AliQCumulantsFunctions.h"

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
 fAvMultIntFlowQC(NULL),
 fQvectorComponents(NULL),
 fIntFlowResultsQC(NULL),
 fDiffFlowResults2ndOrderQC(NULL),
 fDiffFlowResults4thOrderQC(NULL),
 fCovariances(NULL),
 fQCorrelations(NULL),
 fQProduct(NULL),
 fDirectCorrelations(NULL),
 fPtReq1nRP(NULL),
 fPtImq1nRP(NULL),
 fPtReq2nRP(NULL),
 fPtImq2nRP(NULL),
 f2PerPtBin1n1nRP(NULL),
 f2PerPtBin2n2nRP(NULL),
 f3PerPtBin2n1n1nRP(NULL),
 f3PerPtBin1n1n2nRP(NULL),
 f4PerPtBin1n1n1n1nRP(NULL),
 fEtaReq1nRP(NULL),
 fEtaImq1nRP(NULL),
 fEtaReq2nRP(NULL),
 fEtaImq2nRP(NULL),
 f2PerEtaBin1n1nRP(NULL),
 f2PerEtaBin2n2nRP(NULL),
 f3PerEtaBin2n1n1nRP(NULL),
 f3PerEtaBin1n1n2nRP(NULL),
 f4PerEtaBin1n1n1n1nRP(NULL),
 fPtReq1nPOI(NULL),
 fPtImq1nPOI(NULL),
 fPtReq2nPOI(NULL),
 fPtImq2nPOI(NULL),
 f2PerPtBin1n1nPOI(NULL),
 f2PerPtBin2n2nPOI(NULL),
 f3PerPtBin2n1n1nPOI(NULL),
 f3PerPtBin1n1n2nPOI(NULL),
 f4PerPtBin1n1n1n1nPOI(NULL),
 fEtaReq1nPOI(NULL),
 fEtaImq1nPOI(NULL),
 fEtaReq2nPOI(NULL),
 fEtaImq2nPOI(NULL),
 f2PerEtaBin1n1nPOI(NULL),
 f2PerEtaBin2n2nPOI(NULL),
 f3PerEtaBin2n1n1nPOI(NULL),
 f3PerEtaBin1n1n2nPOI(NULL),
 f4PerEtaBin1n1n1n1nPOI(NULL), 
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
 fEtaMax(0)
{
 //constructor 
 fHistList = new TList(); 
 
 fnBinsPt = AliFlowCommonConstants::GetNbinsPt();
 fPtMin   = AliFlowCommonConstants::GetPtMin();	     
 fPtMax   = AliFlowCommonConstants::GetPtMax();
 
 fnBinsEta = AliFlowCommonConstants::GetNbinsEta();
 fEtaMin   = AliFlowCommonConstants::GetEtaMin();	     
 fEtaMax   = AliFlowCommonConstants::GetEtaMax();
}

AliFlowAnalysisWithQCumulants::~AliFlowAnalysisWithQCumulants()
{
 //desctructor
 delete fHistList; 
}

//================================================================================================================

void AliFlowAnalysisWithQCumulants::CreateOutputObjects()
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
 //fIntFlowResultsQC->SetYTitle("Integrated Flow");
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
 
 //multi-particle correlations calculated with nested loops (0..40 integrated flow; 40..80 differential flow)
 fDirectCorrelations = new TProfile("fDirectCorrelations","multi-particle correlations with nested loops",80,0,80,"s");
 fDirectCorrelations->SetXTitle("");
 fDirectCorrelations->SetYTitle("correlations");
 fHistList->Add(fDirectCorrelations);
 
 //fPtReq1nRP
 fPtReq1nRP = new TProfile("fPtReq1nRP","Re[q_n]",fnBinsPt,fPtMin,fPtMax,"s");
 fPtReq1nRP->SetXTitle("p_{t} [GeV]");
 fPtReq1nRP->SetYTitle("Re[q_n]");
 //fHistList->Add(fPtReq1nRP);
 
 //fPtImq1nRP
 fPtImq1nRP = new TProfile("fPtImq1nRP","Im[q_n]",fnBinsPt,fPtMin,fPtMax,"s");
 fPtImq1nRP->SetXTitle("p_{t} [GeV]");
 fPtImq1nRP->SetYTitle("Im[q_n]");
 //fHistList->Add(fPtImq1nRP);
 
 //fPtReq2nRP
 fPtReq2nRP = new TProfile("fPtReq2nRP","Re[q_2n]",fnBinsPt,fPtMin,fPtMax,"s");
 fPtReq2nRP->SetXTitle("p_{t} [GeV]");
 fPtReq2nRP->SetYTitle("Im[D]");
 //fHistList->Add(fPtReq2nRP);
 
 //fPtImq2nRP
 fPtImq2nRP = new TProfile("fPtImq2nRP","Im[q_2n]",fnBinsPt,fPtMin,fPtMax,"s");
 fPtImq2nRP->SetXTitle("p_{t} [GeV]");
 fPtImq2nRP->SetYTitle("Im[q_2n]");
 //fHistList->Add(fPtImq2nRP);
 
 //f2PerPtBin1n1nRP
 f2PerPtBin1n1nRP = new TProfile("f2PerPtBin1n1nRP","<2'>_{n|n}",fnBinsPt,fPtMin,fPtMax,"s");
 f2PerPtBin1n1nRP->SetXTitle("p_{t} [GeV]");
 //f2PerPtBin1n1n->SetYTitle("<2'>_{n|n}");
 fHistList->Add(f2PerPtBin1n1nRP);
 
 //f2PerPtBin2n2nRP
 f2PerPtBin2n2nRP = new TProfile("f2PerPtBin2n2nRP","<2'>_{2n|2n}",fnBinsPt,fPtMin,fPtMax,"s");
 f2PerPtBin2n2nRP->SetXTitle("p_{t} [GeV]");
 //f2PerPtBin2n2nRP->SetYTitle("<2'>_{2n|2n}");
 fHistList->Add(f2PerPtBin2n2nRP);
 
 //f3PerPtBin2n1n1nRP
 f3PerPtBin2n1n1nRP = new TProfile("f3PerPtBin2n1n1nRP","<3'>_{2n|n,n}",fnBinsPt,fPtMin,fPtMax,"s");
 f3PerPtBin2n1n1nRP->SetXTitle("p_{t} [GeV]");
 //f3PerPtBin2n1n1nRP->SetYTitle("<3'>_{2n|n,n}");
 fHistList->Add(f3PerPtBin2n1n1nRP);
 
 //f3PerPtBin1n1n2nRP
 f3PerPtBin1n1n2nRP = new TProfile("f3PerPtBin1n1n2nRP","<3'>_{n,n|2n}",fnBinsPt,fPtMin,fPtMax,"s");
 f3PerPtBin1n1n2nRP->SetXTitle("p_{t} [GeV]");
 //f3PerPtBin1n1n2nRP->SetYTitle("<3'>_{n,n|2n}");
 fHistList->Add(f3PerPtBin1n1n2nRP);
 
 //f4PerPtBin1n1n1n1nRP
 f4PerPtBin1n1n1n1nRP = new TProfile("f4PerPtBin1n1n1n1nRP","<4'>_{n,n|n,n}",fnBinsPt,fPtMin,fPtMax,"s");
 f4PerPtBin1n1n1n1nRP->SetXTitle("p_{t} [GeV]");
 //f4PerPtBin1n1n1n1nRP->SetYTitle("<4'>_{n,n|n,n}");
 fHistList->Add(f4PerPtBin1n1n1n1nRP);
 
 //fEtaReq1nRP
 fEtaReq1nRP = new TProfile("fEtaReq1nRP","Re[q_n]",fnBinsEta,fEtaMin,fEtaMax,"s");
 fEtaReq1nRP->SetXTitle("#eta");
 fEtaReq1nRP->SetYTitle("Re[q_n]");
 fHistList->Add(fEtaReq1nRP);
 
 //fEtaImq1nRP
 fEtaImq1nRP = new TProfile("fEtaImq1nRP","Im[q_n]",fnBinsEta,fEtaMin,fEtaMax,"s");
 fEtaImq1nRP->SetXTitle("#eta");
 fEtaImq1nRP->SetYTitle("Im[q_n]");
 //fHistList->Add(fEtaImq1nRP);
 
 //fEtaReq2nRP
 fEtaReq2nRP = new TProfile("fEtaReq2nRP","Re[q_2n]",fnBinsEta,fEtaMin,fEtaMax,"s");
 fEtaReq2nRP->SetXTitle("#eta");
 fEtaReq2nRP->SetYTitle("Im[D]");
 //fHistList->Add(fEtaReq2nRP);
 
 //fEtaImq2nRP
 fEtaImq2nRP = new TProfile("fEtaImq2nRP","Im[q_2n]",fnBinsEta,fEtaMin,fEtaMax,"s");
 fEtaImq2nRP->SetXTitle("#eta");
 fEtaImq2nRP->SetYTitle("Im[q_2n]");
 //fHistList->Add(fEtaImq2nRP);
 
 //f2PerEtaBin1n1nRP
 f2PerEtaBin1n1nRP = new TProfile("f2PerEtaBin1n1nRP","<2'>_{n|n}",fnBinsEta,fEtaMin,fEtaMax,"s");
 f2PerEtaBin1n1nRP->SetXTitle("#eta");
 //f2PerEtaBin1n1nRP->SetYTitle("<2'>_{n|n}");
 fHistList->Add(f2PerEtaBin1n1nRP);
 
 //f2PerEtaBin2n2nRP
 f2PerEtaBin2n2nRP = new TProfile("f2PerEtaBin2n2nRP","<2'>_{2n|2n}",fnBinsEta,fEtaMin,fEtaMax,"s");
 f2PerEtaBin2n2nRP->SetXTitle("#eta");
 //f2PerEtaBin2n2nRP->SetYTitle("<2'>_{2n|2n}");
 fHistList->Add(f2PerEtaBin2n2nRP);
 
 //f3PerEtaBin2n1n1nRP
 f3PerEtaBin2n1n1nRP = new TProfile("f3PerEtaBin2n1n1nRP","<3'>_{2n|n,n}",fnBinsEta,fEtaMin,fEtaMax,"s");
 f3PerEtaBin2n1n1nRP->SetXTitle("#eta");
 //f3PerEtaBin2n1n1nRP->SetYTitle("<3'>_{2n|n,n}");
 fHistList->Add(f3PerEtaBin2n1n1nRP);
 
 //f3PerEtaBin1n1n2nRP
 f3PerEtaBin1n1n2nRP = new TProfile("f3PerEtaBin1n1n2RP","<3'>_{n,n|2n}",fnBinsEta,fEtaMin,fEtaMax,"s");
 f3PerEtaBin1n1n2nRP->SetXTitle("#eta");
 //f3PerEtaBin1n1n2n->SetYTitle("<3'>_{n,n|2n}");
 fHistList->Add(f3PerEtaBin1n1n2nRP);
 
 //f4PerEtaBin1n1n1n1nRP
 f4PerEtaBin1n1n1n1nRP = new TProfile("f4PerEtaBin1n1n1n1nRP","<4'>_{n,n|n,n}",fnBinsEta,fEtaMin,fEtaMax,"s");
 f4PerEtaBin1n1n1n1nRP->SetXTitle("#eta");
 //f4PerEtaBin1n1n1n1nRP->SetYTitle("<4'>_{n,n|n,n}");
 fHistList->Add(f4PerEtaBin1n1n1n1nRP);
 
 //fPtReq1nPOI
 fPtReq1nPOI = new TProfile("fPtReq1nPOI","Re[q_n]",fnBinsPt,fPtMin,fPtMax,"s");
 fPtReq1nPOI->SetXTitle("p_{t} [GeV]");
 fPtReq1nPOI->SetYTitle("Re[q_n]");
 //fHistList->Add(fPtReq1nPOI);
 
 //fPtImq1nPOI
 fPtImq1nPOI = new TProfile("fPtImq1nPOI","Im[q_n]",fnBinsPt,fPtMin,fPtMax,"s");
 fPtImq1nPOI->SetXTitle("p_{t} [GeV]");
 fPtImq1nPOI->SetYTitle("Im[q_n]");
 //fHistList->Add(fPtImq1nPOI);
 
 //fPtReq2nPOI
 fPtReq2nPOI = new TProfile("fPtReq2nPOI","Re[q_2n]",fnBinsPt,fPtMin,fPtMax,"s");
 fPtReq2nPOI->SetXTitle("p_{t} [GeV]");
 fPtReq2nPOI->SetYTitle("Im[D]");
 //fHistList->Add(fPtReq2nPOI);
 
 //fPtImq2nPOI
 fPtImq2nPOI = new TProfile("fPtImq2nPOI","Im[q_2n]",fnBinsPt,fPtMin,fPtMax,"s");
 fPtImq2nPOI->SetXTitle("p_{t} [GeV]");
 fPtImq2nPOI->SetYTitle("Im[q_2n]");
 //fHistList->Add(fPtImq2nPOI);
 
 //f2PerPtBin1n1nPOI
 f2PerPtBin1n1nPOI = new TProfile("f2PerPtBin1n1nPOI","<2'>_{n|n}",fnBinsPt,fPtMin,fPtMax,"s");
 f2PerPtBin1n1nPOI->SetXTitle("p_{t} [GeV]");
 //f2PerPtBin1n1n->SetYTitle("<2'>_{n|n}");
 fHistList->Add(f2PerPtBin1n1nPOI);
 
 //f2PerPtBin2n2nPOI
 f2PerPtBin2n2nPOI = new TProfile("f2PerPtBin2n2nPOI","<2'>_{2n|2n}",fnBinsPt,fPtMin,fPtMax,"s");
 f2PerPtBin2n2nPOI->SetXTitle("p_{t} [GeV]");
 //f2PerPtBin2n2nPOI->SetYTitle("<2'>_{2n|2n}");
 fHistList->Add(f2PerPtBin2n2nPOI);
 
 //f3PerPtBin2n1n1nPOI
 f3PerPtBin2n1n1nPOI = new TProfile("f3PerPtBin2n1n1nPOI","<3'>_{2n|n,n}",fnBinsPt,fPtMin,fPtMax,"s");
 f3PerPtBin2n1n1nPOI->SetXTitle("p_{t} [GeV]");
 //f3PerPtBin2n1n1nPOI->SetYTitle("<3'>_{2n|n,n}");
 fHistList->Add(f3PerPtBin2n1n1nPOI);
 
 //f3PerPtBin1n1n2nPOI
 f3PerPtBin1n1n2nPOI = new TProfile("f3PerPtBin1n1n2nPOI","<3'>_{n,n|2n}",fnBinsPt,fPtMin,fPtMax,"s");
 f3PerPtBin1n1n2nPOI->SetXTitle("p_{t} [GeV]");
 //f3PerPtBin1n1n2nPOI->SetYTitle("<3'>_{n,n|2n}");
 fHistList->Add(f3PerPtBin1n1n2nPOI);
 
 //f4PerPtBin1n1n1n1nPOI
 f4PerPtBin1n1n1n1nPOI = new TProfile("f4PerPtBin1n1n1n1nPOI","<4'>_{n,n|n,n}",fnBinsPt,fPtMin,fPtMax,"s");
 f4PerPtBin1n1n1n1nPOI->SetXTitle("p_{t} [GeV]");
 //f4PerPtBin1n1n1n1nPOI->SetYTitle("<4'>_{n,n|n,n}");
 fHistList->Add(f4PerPtBin1n1n1n1nPOI);
 
 //fEtaReq1nPOI
 fEtaReq1nPOI = new TProfile("fEtaReq1nPOI","Re[q_n]",fnBinsEta,fEtaMin,fEtaMax,"s");
 fEtaReq1nPOI->SetXTitle("#eta");
 fEtaReq1nPOI->SetYTitle("Re[q_n]");
 fHistList->Add(fEtaReq1nPOI);
 
 //fEtaImq1nPOI
 fEtaImq1nPOI = new TProfile("fEtaImq1nPOI","Im[q_n]",fnBinsEta,fEtaMin,fEtaMax,"s");
 fEtaImq1nPOI->SetXTitle("#eta");
 fEtaImq1nPOI->SetYTitle("Im[q_n]");
 //fHistList->Add(fEtaImq1nPOI);
 
 //fEtaReq2nPOI
 fEtaReq2nPOI = new TProfile("fEtaReq2nPOI","Re[q_2n]",fnBinsEta,fEtaMin,fEtaMax,"s");
 fEtaReq2nPOI->SetXTitle("#eta");
 fEtaReq2nPOI->SetYTitle("Im[D]");
 //fHistList->Add(fEtaReq2nPOI);
 
 //fEtaImq2nPOI
 fEtaImq2nPOI = new TProfile("fEtaImq2nPOI","Im[q_2n]",fnBinsEta,fEtaMin,fEtaMax,"s");
 fEtaImq2nPOI->SetXTitle("#eta");
 fEtaImq2nPOI->SetYTitle("Im[q_2n]");
 //fHistList->Add(fEtaImq2nPOI);
 
 //f2PerEtaBin1n1nPOI
 f2PerEtaBin1n1nPOI = new TProfile("f2PerEtaBin1n1nPOI","<2'>_{n|n}",fnBinsEta,fEtaMin,fEtaMax,"s");
 f2PerEtaBin1n1nPOI->SetXTitle("#eta");
 //f2PerEtaBin1n1nPOI->SetYTitle("<2'>_{n|n}");
 fHistList->Add(f2PerEtaBin1n1nPOI);
 
 //f2PerEtaBin2n2nPOI
 f2PerEtaBin2n2nPOI = new TProfile("f2PerEtaBin2n2nPOI","<2'>_{2n|2n}",fnBinsEta,fEtaMin,fEtaMax,"s");
 f2PerEtaBin2n2nPOI->SetXTitle("#eta");
 //f2PerEtaBin2n2nPOI->SetYTitle("<2'>_{2n|2n}");
 fHistList->Add(f2PerEtaBin2n2nPOI);
 
 //f3PerEtaBin2n1n1nPOI
 f3PerEtaBin2n1n1nPOI = new TProfile("f3PerEtaBin2n1n1nPOI","<3'>_{2n|n,n}",fnBinsEta,fEtaMin,fEtaMax,"s");
 f3PerEtaBin2n1n1nPOI->SetXTitle("#eta");
 //f3PerEtaBin2n1n1nPOI->SetYTitle("<3'>_{2n|n,n}");
 fHistList->Add(f3PerEtaBin2n1n1nPOI);
 
 //f3PerEtaBin1n1n2nPOI
 f3PerEtaBin1n1n2nPOI = new TProfile("f3PerEtaBin1n1n2POI","<3'>_{n,n|2n}",fnBinsEta,fEtaMin,fEtaMax,"s");
 f3PerEtaBin1n1n2nPOI->SetXTitle("#eta");
 //f3PerEtaBin1n1n2n->SetYTitle("<3'>_{n,n|2n}");
 fHistList->Add(f3PerEtaBin1n1n2nPOI);
 
 //f4PerEtaBin1n1n1n1nPOI
 f4PerEtaBin1n1n1n1nPOI = new TProfile("f4PerEtaBin1n1n1n1nPOI","<4'>_{n,n|n,n}",fnBinsEta,fEtaMin,fEtaMax,"s");
 f4PerEtaBin1n1n1n1nPOI->SetXTitle("#eta");
 //f4PerEtaBin1n1n1n1nPOI->SetYTitle("<4'>_{n,n|n,n}");
 fHistList->Add(f4PerEtaBin1n1n1n1nPOI);
 
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
 
}//end of CreateOutputObjects()

//================================================================================================================

void AliFlowAnalysisWithQCumulants::Make(AliFlowEventSimple* anEvent)
{
 //running over data         
          
 //get the total multiplicity of event:
 //Int_t nPrim = anEvent->NumberOfTracks();//line needed only for nested loops
 //cout<<"nPrim = "<<nPrim<<endl;

 //if(nPrim>8&&nPrim<12)//line needed only for nested loops  
 //{                    //line needed only for nested loops
 
 //get the selected multiplicity (i.e. number of particles used for int. flow):
 //Int_t nEventNSelTracksIntFlow = anEvent->GetEventNSelTracksIntFlow();
 
 Int_t n=2; //int flow harmonic (to be improved)
 
 //---------------------------------------------------------------------------------------------------------
 //Q-vectors of an event evaluated in harmonics n, 2n, 3n and 4n:
 AliFlowVector afvQvector1n, afvQvector2n, afvQvector3n, afvQvector4n;
 
 afvQvector1n.Set(0.,0.);
 afvQvector1n.SetMult(0);
 afvQvector1n=anEvent->GetQ(1*n); 
 
 afvQvector2n.Set(0.,0.);
 afvQvector2n.SetMult(0);
 afvQvector2n=anEvent->GetQ(2*n);                          
 
 afvQvector3n.Set(0.,0.);
 afvQvector3n.SetMult(0);
 afvQvector3n=anEvent->GetQ(3*n);       
 
 afvQvector4n.Set(0.,0.);
 afvQvector4n.SetMult(0);
 afvQvector4n=anEvent->GetQ(4*n);       
 //---------------------------------------------------------------------------------------------------------
 
 //multiplicity (to be improved, because I already have nEventNSelTracksIntFlow and nPrim)
 Double_t dMult = afvQvector1n.GetMult();
 
 fAvMultIntFlowQC->Fill(0.,dMult,1.);
 
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
         
 //2-particle
 Double_t two1n1n=0., two2n2n=0., two3n3n=0., two4n4n=0.; 
 if(dMult>1)
 {
  //fill the common control histogram (2nd order):
  fCommonHists2nd->FillControlHistograms(anEvent); 
 
  two1n1n = (pow(afvQvector1n.Mod(),2.)-dMult)/(dMult*(dMult-1.)); //<2>_{n|n}   = <cos(n*(phi1-phi2))>
  two2n2n = (pow(afvQvector2n.Mod(),2.)-dMult)/(dMult*(dMult-1.)); //<2>_{2n|2n} = <cos(2n*(phi1-phi2))>
  two3n3n = (pow(afvQvector3n.Mod(),2.)-dMult)/(dMult*(dMult-1.)); //<2>_{3n|3n} = <cos(3n*(phi1-phi2))>
  two4n4n = (pow(afvQvector4n.Mod(),2.)-dMult)/(dMult*(dMult-1.)); //<2>_{4n|4n} = <cos(4n*(phi1-phi2))>
    
  fQCorrelations->Fill(0.,two1n1n,dMult*(dMult-1.)); 
  fQCorrelations->Fill(1.,two2n2n,dMult*(dMult-1.)); 
  fQCorrelations->Fill(2.,two3n3n,dMult*(dMult-1.)); 
  fQCorrelations->Fill(3.,two4n4n,dMult*(dMult-1.)); 
  
  f2pDistribution->Fill(two1n1n,dMult*(dMult-1.)); 
 }
 
 //3-particle
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
 // DIFFERENTIAL FLOW
 
 Double_t dQ1nx = afvQvector1n.X();
 Double_t dQ1ny = afvQvector1n.Y();
 Double_t dQ2nx = afvQvector2n.X();
 Double_t dQ2ny = afvQvector2n.Y();
 
 Double_t dBinWidthPt=1.*(fPtMax-fPtMin)/fnBinsPt;
 Double_t dBinWidthEta=1.*(fEtaMax-fEtaMin)/fnBinsEta;

 //RP:
 Double_t qxPtRP=0.,qyPtRP=0.,q2xPtRP=0.,q2yPtRP=0.,mPtRP=0.;//add comments for these variables (deleteMe)
 Double_t qxEtaRP=0.,qyEtaRP=0.,q2xEtaRP=0.,q2yEtaRP=0.,mEtaRP=0.;//add comments for these variables (deleteMe)
  
 for(Int_t i=0;i<dMult;i++) //check if nPrim == M
 { 
  fTrack=anEvent->GetTrack(i);
  if(fTrack && fTrack->UseForIntegratedFlow())//checking RP condition 
  {
   //Pt:
   fPtReq1nRP->Fill(fTrack->Pt(),cos(n*(fTrack->Phi())),1.);
   fPtImq1nRP->Fill(fTrack->Pt(),sin(n*(fTrack->Phi())),1.);
   fPtReq2nRP->Fill(fTrack->Pt(),cos(2.*n*(fTrack->Phi())),1.);
   fPtImq2nRP->Fill(fTrack->Pt(),sin(2.*n*(fTrack->Phi())),1.);
   //Eta:
   fEtaReq1nRP->Fill(fTrack->Eta(),cos(n*(fTrack->Phi())),1.);
   fEtaImq1nRP->Fill(fTrack->Eta(),sin(n*(fTrack->Phi())),1.);
   fEtaReq2nRP->Fill(fTrack->Eta(),cos(2.*n*(fTrack->Phi())),1.);
   fEtaImq2nRP->Fill(fTrack->Eta(),sin(2.*n*(fTrack->Phi())),1.);
  }
 } 
  
 //Pt:
 Double_t twoDiffPt1n1nRP=0.,twoDiffPt2n2nRP=0.,threeDiffPt2n1n1nRP=0.,threeDiffPt1n1n2nRP=0.,fourDiffPt1n1n1n1nRP=0.;
 
 for(Int_t bin=1;bin<(fnBinsPt+1);bin++)//loop over pt-bins 
 { 
  qxPtRP = (fPtReq1nRP->GetBinContent(bin))*(fPtReq1nRP->GetBinEntries(bin));
  qyPtRP = (fPtImq1nRP->GetBinContent(bin))*(fPtImq1nRP->GetBinEntries(bin)); 
  q2xPtRP = (fPtReq2nRP->GetBinContent(bin))*(fPtReq2nRP->GetBinEntries(bin));  
  q2yPtRP = (fPtImq2nRP->GetBinContent(bin))*(fPtImq2nRP->GetBinEntries(bin)); 
  mPtRP = fPtReq1nRP->GetBinEntries(bin);          
 
  if(mPtRP>0&&dMult>1)
  {
   twoDiffPt1n1nRP = (qxPtRP*dQ1nx+qyPtRP*dQ1ny-mPtRP)/(mPtRP*(dMult-1.));
   f2PerPtBin1n1nRP->Fill(fPtMin+(bin-1)*dBinWidthPt,twoDiffPt1n1nRP,mPtRP*(dMult-1.));//<2'>_{n|n}
   
   twoDiffPt2n2nRP = (q2xPtRP*dQ2nx+q2yPtRP*dQ2ny-mPtRP)/(mPtRP*(dMult-1.));
   f2PerPtBin2n2nRP->Fill(fPtMin+(bin-1)*dBinWidthPt,twoDiffPt2n2nRP,mPtRP*(dMult-1.));//<2'>_{2n|2n} 
  }
  
  if(mPtRP>0&&dMult>2)
  {
   threeDiffPt2n1n1nRP = (q2xPtRP*(dQ1nx*dQ1nx-dQ1ny*dQ1ny)+2.*q2yPtRP*dQ1nx*dQ1ny-2.*(qxPtRP*dQ1nx+qyPtRP*dQ1ny)-(q2xPtRP*dQ2nx+q2yPtRP*dQ2ny)+2.*mPtRP)/(mPtRP*(dMult-1.)*(dMult-2.));
   f3PerPtBin2n1n1nRP->Fill(fPtMin+(bin-1)*dBinWidthPt,threeDiffPt2n1n1nRP,mPtRP*(dMult-1.)*(dMult-2.));//Re[<3'>_{2n|n,n}]
   
   threeDiffPt1n1n2nRP = (dQ2nx*(qxPtRP*dQ1nx-qyPtRP*dQ1ny)+dQ2ny*(qxPtRP*dQ1ny+qyPtRP*dQ1nx)-2.*(qxPtRP*dQ1nx+qyPtRP*dQ1ny)-(q2xPtRP*dQ2nx+q2yPtRP*dQ2ny)+2.*mPtRP)/(mPtRP*(dMult-1.)*(dMult-2.));
   f3PerPtBin1n1n2nRP->Fill(fPtMin+(bin-1)*dBinWidthPt,threeDiffPt1n1n2nRP,mPtRP*(dMult-1.)*(dMult-2.));//Re[<3'>_{n,n|2n}]
  }
  
  if(mPtRP>0&&dMult>3)
  {
   fourDiffPt1n1n1n1nRP = ((dQ1nx*dQ1nx+dQ1ny*dQ1ny)*(qxPtRP*dQ1nx+qyPtRP*dQ1ny)-(q2xPtRP*(dQ1nx*dQ1nx-dQ1ny*dQ1ny)+2.*q2yPtRP*dQ1nx*dQ1ny)-(dQ2nx*(qxPtRP*dQ1nx-qyPtRP*dQ1ny)+dQ2ny*(qxPtRP*dQ1ny+qyPtRP*dQ1nx))+(q2xPtRP*dQ2nx+q2yPtRP*dQ2ny)-2.*(dMult-3.)*(qxPtRP*dQ1nx+qyPtRP*dQ1ny)-2.*mPtRP*(dQ1nx*dQ1nx+dQ1ny*dQ1ny)+2.*(dQ1nx*qxPtRP+dQ1ny*qyPtRP)+2.*mPtRP*(dMult-3.))/(mPtRP*(dMult-1.)*(dMult-2.)*(dMult-3.));
   f4PerPtBin1n1n1n1nRP->Fill(fPtMin+(bin-1)*dBinWidthPt,fourDiffPt1n1n1n1nRP,mPtRP*(dMult-1.)*(dMult-2.)*(dMult-3.));//Re[<4'>_{n,n|n,n}]
  }
   
 } 
  
 fPtReq1nRP->Reset();
 fPtImq1nRP->Reset();
 fPtReq2nRP->Reset();
 fPtImq2nRP->Reset();
 
 //Eta:
 Double_t twoDiffEta1n1nRP=0.,twoDiffEta2n2nRP=0.,threeDiffEta2n1n1nRP=0.,threeDiffEta1n1n2nRP=0.,fourDiffEta1n1n1n1nRP=0.;
 
 for(Int_t bin=1;bin<(fnBinsEta+1);bin++)//loop over eta-bins 
 { 
  qxEtaRP = (fEtaReq1nRP->GetBinContent(bin))*(fEtaReq1nRP->GetBinEntries(bin));
  qyEtaRP = (fEtaImq1nRP->GetBinContent(bin))*(fEtaImq1nRP->GetBinEntries(bin)); 
  q2xEtaRP = (fEtaReq2nRP->GetBinContent(bin))*(fEtaReq2nRP->GetBinEntries(bin));  
  q2yEtaRP = (fEtaImq2nRP->GetBinContent(bin))*(fEtaImq2nRP->GetBinEntries(bin)); 
  mEtaRP = fEtaReq1nRP->GetBinEntries(bin); 
  
  if(mEtaRP>0&&dMult>1)
  {
   twoDiffEta1n1nRP = (qxEtaRP*dQ1nx+qyEtaRP*dQ1ny-mEtaRP)/(mEtaRP*(dMult-1.));
   f2PerEtaBin1n1nRP->Fill(fEtaMin+(bin-1)*dBinWidthEta,twoDiffEta1n1nRP,mEtaRP*(dMult-1.));//<2'>_{n|n}
   
   twoDiffEta2n2nRP = (q2xEtaRP*dQ2nx+q2yEtaRP*dQ2ny-mEtaRP)/(mEtaRP*(dMult-1.));
   f2PerEtaBin2n2nRP->Fill(fEtaMin+(bin-1)*dBinWidthEta,twoDiffEta2n2nRP,mEtaRP*(dMult-1.));//<2'>_{2n|2n} 
  }
  
  if(mEtaRP>0&&dMult>2)
  {
   threeDiffEta2n1n1nRP = (q2xEtaRP*(dQ1nx*dQ1nx-dQ1ny*dQ1ny)+2.*q2yEtaRP*dQ1nx*dQ1ny-2.*(qxEtaRP*dQ1nx+qyEtaRP*dQ1ny)-(q2xEtaRP*dQ2nx+q2yEtaRP*dQ2ny)+2.*mEtaRP)/(mEtaRP*(dMult-1.)*(dMult-2.));
   f3PerEtaBin2n1n1nRP->Fill(fEtaMin+(bin-1)*dBinWidthEta,threeDiffEta2n1n1nRP,mEtaRP*(dMult-1.)*(dMult-2.));//Re[<3'>_{2n|n,n}]
   
   threeDiffEta1n1n2nRP = (dQ2nx*(qxEtaRP*dQ1nx-qyEtaRP*dQ1ny)+dQ2ny*(qxEtaRP*dQ1ny+qyEtaRP*dQ1nx)-2.*(qxEtaRP*dQ1nx+qyEtaRP*dQ1ny)-(q2xEtaRP*dQ2nx+q2yEtaRP*dQ2ny)+2.*mEtaRP)/(mEtaRP*(dMult-1.)*(dMult-2.));
   f3PerEtaBin1n1n2nRP->Fill(fEtaMin+(bin-1)*dBinWidthEta,threeDiffEta1n1n2nRP,mEtaRP*(dMult-1.)*(dMult-2.));//Re[<3'>_{n,n|2n}]
  }
  
  if(mEtaRP>0&&dMult>3)
  {
   fourDiffEta1n1n1n1nRP = ((dQ1nx*dQ1nx+dQ1ny*dQ1ny)*(qxEtaRP*dQ1nx+qyEtaRP*dQ1ny)-(q2xEtaRP*(dQ1nx*dQ1nx-dQ1ny*dQ1ny)+2.*q2yEtaRP*dQ1nx*dQ1ny)-(dQ2nx*(qxEtaRP*dQ1nx-qyEtaRP*dQ1ny)+dQ2ny*(qxEtaRP*dQ1ny+qyEtaRP*dQ1nx))+(q2xEtaRP*dQ2nx+q2yEtaRP*dQ2ny)-2.*(dMult-3.)*(qxEtaRP*dQ1nx+qyEtaRP*dQ1ny)-2.*mEtaRP*(dQ1nx*dQ1nx+dQ1ny*dQ1ny)+2.*(dQ1nx*qxEtaRP+dQ1ny*qyEtaRP)+2.*mEtaRP*(dMult-3.))/(mEtaRP*(dMult-1.)*(dMult-2.)*(dMult-3.));
   f4PerEtaBin1n1n1n1nRP->Fill(fEtaMin+(bin-1)*dBinWidthEta,fourDiffEta1n1n1n1nRP,mEtaRP*(dMult-1.)*(dMult-2.)*(dMult-3.));//Re[<4'>_{n,n|n,n}]
  }
   
 } 
  
 fEtaReq1nRP->Reset();
 fEtaImq1nRP->Reset();
 fEtaReq2nRP->Reset();
 fEtaImq2nRP->Reset();
 
 //POI:
 Double_t qxPtPOI=0.,qyPtPOI=0.,q2xPtPOI=0.,q2yPtPOI=0.,mPtPOI=0.;//add comments for these variables (deleteMe)
 Double_t qxEtaPOI=0.,qyEtaPOI=0.,q2xEtaPOI=0.,q2yEtaPOI=0.,mEtaPOI=0.;//add comments for these variables (deleteMe)
 Int_t iOverlap=0; 
  
 for(Int_t i=0;i<dMult;i++) //check if nPrim == M
 { 
  fTrack=anEvent->GetTrack(i);
  if(fTrack && fTrack->UseForDifferentialFlow())//checking POI condition 
  {
   //Pt:
   fPtReq1nPOI->Fill(fTrack->Pt(),cos(n*(fTrack->Phi())),1.);
   fPtImq1nPOI->Fill(fTrack->Pt(),sin(n*(fTrack->Phi())),1.);
   fPtReq2nPOI->Fill(fTrack->Pt(),cos(2.*n*(fTrack->Phi())),1.);
   fPtImq2nPOI->Fill(fTrack->Pt(),sin(2.*n*(fTrack->Phi())),1.);
   //Eta:
   fEtaReq1nPOI->Fill(fTrack->Eta(),cos(n*(fTrack->Phi())),1.);
   fEtaImq1nPOI->Fill(fTrack->Eta(),sin(n*(fTrack->Phi())),1.);
   fEtaReq2nPOI->Fill(fTrack->Eta(),cos(2.*n*(fTrack->Phi())),1.);
   fEtaImq2nPOI->Fill(fTrack->Eta(),sin(2.*n*(fTrack->Phi())),1.);       
   if(fTrack->UseForDifferentialFlow() && fTrack->UseForIntegratedFlow())//counting the overlap between RP and POI
   {
    iOverlap++;
   } 
  }
 } 
 
 //Pt:
 Double_t twoDiffPt1n1nPOI=0.,twoDiffPt2n2nPOI=0.,threeDiffPt2n1n1nPOI=0.,threeDiffPt1n1n2nPOI=0.,fourDiffPt1n1n1n1nPOI=0.;
 
 for(Int_t bin=1;bin<(fnBinsPt+1);bin++)//loop over pt-bins 
 { 
  qxPtPOI = (fPtReq1nPOI->GetBinContent(bin))*(fPtReq1nPOI->GetBinEntries(bin));
  qyPtPOI = (fPtImq1nPOI->GetBinContent(bin))*(fPtImq1nPOI->GetBinEntries(bin)); 
  q2xPtPOI = (fPtReq2nPOI->GetBinContent(bin))*(fPtReq2nPOI->GetBinEntries(bin));  
  q2yPtPOI = (fPtImq2nPOI->GetBinContent(bin))*(fPtImq2nPOI->GetBinEntries(bin)); 
  mPtPOI = fPtReq1nPOI->GetBinEntries(bin);          
 
  if(mPtPOI>0&&dMult>1)
  {
   //twoDiffPt1n1nPOI = (qxPtPOI*dQ1nx+qyPtPOI*dQ1ny-mPtPOI)/(mPtPOI*(dMult-1.));
   twoDiffPt1n1nPOI = (qxPtPOI*dQ1nx+qyPtPOI*dQ1ny-iOverlap)/(mPtPOI*(dMult-1.));   
   f2PerPtBin1n1nPOI->Fill(fPtMin+(bin-1)*dBinWidthPt,twoDiffPt1n1nPOI,mPtPOI*(dMult-1.));//<2'>_{n|n}
   
   //twoDiffPt2n2nPOI = (q2xPtPOI*dQ2nx+q2yPtPOI*dQ2ny-mPtPOI)/(mPtPOI*(dMult-1.));
   twoDiffPt2n2nPOI = (q2xPtPOI*dQ2nx+q2yPtPOI*dQ2ny-iOverlap)/(mPtPOI*(dMult-1.));   
   f2PerPtBin2n2nPOI->Fill(fPtMin+(bin-1)*dBinWidthPt,twoDiffPt2n2nPOI,mPtPOI*(dMult-1.));//<2'>_{2n|2n} 
  }
  
  if(mPtPOI>0&&dMult>2)
  {
   threeDiffPt2n1n1nPOI = (q2xPtPOI*(dQ1nx*dQ1nx-dQ1ny*dQ1ny)+2.*q2yPtPOI*dQ1nx*dQ1ny-2.*(qxPtPOI*dQ1nx+qyPtPOI*dQ1ny)-(q2xPtPOI*dQ2nx+q2yPtPOI*dQ2ny)+2.*mPtPOI)/(mPtPOI*(dMult-1.)*(dMult-2.));//to be improved (correct formula)
   //f3PePOItBin2n1n1nPOI->Fill(fPtMin+(bin-1)*dBinWidthPt,threeDiffPt2n1n1nPOI,mPtPOI*(dMult-1.)*(dMult-2.));//Re[<3'>_{2n|n,n}]
   
   threeDiffPt1n1n2nPOI = (dQ2nx*(qxPtPOI*dQ1nx-qyPtPOI*dQ1ny)+dQ2ny*(qxPtPOI*dQ1ny+qyPtPOI*dQ1nx)-2.*(qxPtPOI*dQ1nx+qyPtPOI*dQ1ny)-(q2xPtPOI*dQ2nx+q2yPtPOI*dQ2ny)+2.*mPtPOI)/(mPtPOI*(dMult-1.)*(dMult-2.));//to be improved (correct formula)
   //f3PePOItBin1n1n2nPOI->Fill(fPtMin+(bin-1)*dBinWidthPt,threeDiffPt1n1n2nPOI,mPtPOI*(dMult-1.)*(dMult-2.));//Re[<3'>_{n,n|2n}]
  }
  
  if(mPtPOI>0&&dMult>3)
  {
   fourDiffPt1n1n1n1nPOI = ((dQ1nx*dQ1nx+dQ1ny*dQ1ny)*(qxPtPOI*dQ1nx+qyPtPOI*dQ1ny)-(q2xPtPOI*(dQ1nx*dQ1nx-dQ1ny*dQ1ny)+2.*q2yPtPOI*dQ1nx*dQ1ny)-(dQ2nx*(qxPtPOI*dQ1nx-qyPtPOI*dQ1ny)+dQ2ny*(qxPtPOI*dQ1ny+qyPtPOI*dQ1nx))+(q2xPtPOI*dQ2nx+q2yPtPOI*dQ2ny)-2.*(dMult-3.)*(qxPtPOI*dQ1nx+qyPtPOI*dQ1ny)-2.*mPtPOI*(dQ1nx*dQ1nx+dQ1ny*dQ1ny)+2.*(dQ1nx*qxPtPOI+dQ1ny*qyPtPOI)+2.*mPtPOI*(dMult-3.))/(mPtPOI*(dMult-1.)*(dMult-2.)*(dMult-3.));//to be improved (correct formula)
   //f4PePOItBin1n1n1n1nPOI->Fill(fPtMin+(bin-1)*dBinWidthPt,fourDiffPt1n1n1n1nPOI,mPtPOI*(dMult-1.)*(dMult-2.)*(dMult-3.));//Re[<4'>_{n,n|n,n}]
  }
   
 } 
  
 fPtReq1nPOI->Reset();
 fPtImq1nPOI->Reset();
 fPtReq2nPOI->Reset();
 fPtImq2nPOI->Reset();
 
 //Eta:
 Double_t twoDiffEta1n1nPOI=0.,twoDiffEta2n2nPOI=0.,threeDiffEta2n1n1nPOI=0.,threeDiffEta1n1n2nPOI=0.,fourDiffEta1n1n1n1nPOI=0.;
 
 for(Int_t bin=1;bin<(fnBinsEta+1);bin++)//loop over eta-bins 
 { 
  qxEtaPOI = (fEtaReq1nPOI->GetBinContent(bin))*(fEtaReq1nPOI->GetBinEntries(bin));
  qyEtaPOI = (fEtaImq1nPOI->GetBinContent(bin))*(fEtaImq1nPOI->GetBinEntries(bin)); 
  q2xEtaPOI = (fEtaReq2nPOI->GetBinContent(bin))*(fEtaReq2nPOI->GetBinEntries(bin));  
  q2yEtaPOI = (fEtaImq2nPOI->GetBinContent(bin))*(fEtaImq2nPOI->GetBinEntries(bin)); 
  mEtaPOI = fEtaReq1nPOI->GetBinEntries(bin); 
  
  if(mEtaPOI>0&&dMult>1)
  {
   //twoDiffEta1n1nPOI = (qxEtaPOI*dQ1nx+qyEtaPOI*dQ1ny-mEtaPOI)/(mEtaPOI*(dMult-1.));
   twoDiffEta1n1nPOI = (qxEtaPOI*dQ1nx+qyEtaPOI*dQ1ny-iOverlap)/(mEtaPOI*(dMult-1.));
   f2PerEtaBin1n1nPOI->Fill(fEtaMin+(bin-1)*dBinWidthEta,twoDiffEta1n1nPOI,mEtaPOI*(dMult-1.));//<2'>_{n|n}
   
   //twoDiffEta2n2nPOI = (q2xEtaPOI*dQ2nx+q2yEtaPOI*dQ2ny-mEtaPOI)/(mEtaPOI*(dMult-1.));
   twoDiffEta2n2nPOI = (q2xEtaPOI*dQ2nx+q2yEtaPOI*dQ2ny-iOverlap)/(mEtaPOI*(dMult-1.));
   f2PerEtaBin2n2nPOI->Fill(fEtaMin+(bin-1)*dBinWidthEta,twoDiffEta2n2nPOI,mEtaPOI*(dMult-1.));//<2'>_{2n|2n} 
  }
  
  if(mEtaPOI>0&&dMult>2)
  {
   threeDiffEta2n1n1nPOI = (q2xEtaPOI*(dQ1nx*dQ1nx-dQ1ny*dQ1ny)+2.*q2yEtaPOI*dQ1nx*dQ1ny-2.*(qxEtaPOI*dQ1nx+qyEtaPOI*dQ1ny)-(q2xEtaPOI*dQ2nx+q2yEtaPOI*dQ2ny)+2.*mEtaPOI)/(mEtaPOI*(dMult-1.)*(dMult-2.));//to be improved (correct formula)
   //f3PerEtaBin2n1n1nPOI->Fill(fEtaMin+(bin-1)*dBinWidthEta,threeDiffEta2n1n1nPOI,mEtaPOI*(dMult-1.)*(dMult-2.));//Re[<3'>_{2n|n,n}]
   
   threeDiffEta1n1n2nPOI = (dQ2nx*(qxEtaPOI*dQ1nx-qyEtaPOI*dQ1ny)+dQ2ny*(qxEtaPOI*dQ1ny+qyEtaPOI*dQ1nx)-2.*(qxEtaPOI*dQ1nx+qyEtaPOI*dQ1ny)-(q2xEtaPOI*dQ2nx+q2yEtaPOI*dQ2ny)+2.*mEtaPOI)/(mEtaPOI*(dMult-1.)*(dMult-2.));//to be improved (correct formula)
   //f3PerEtaBin1n1n2nPOI->Fill(fEtaMin+(bin-1)*dBinWidthEta,threeDiffEta1n1n2nPOI,mEtaPOI*(dMult-1.)*(dMult-2.));//Re[<3'>_{n,n|2n}]
  }
  
  if(mEtaPOI>0&&dMult>3)
  {
   fourDiffEta1n1n1n1nPOI = ((dQ1nx*dQ1nx+dQ1ny*dQ1ny)*(qxEtaPOI*dQ1nx+qyEtaPOI*dQ1ny)-(q2xEtaPOI*(dQ1nx*dQ1nx-dQ1ny*dQ1ny)+2.*q2yEtaPOI*dQ1nx*dQ1ny)-(dQ2nx*(qxEtaPOI*dQ1nx-qyEtaPOI*dQ1ny)+dQ2ny*(qxEtaPOI*dQ1ny+qyEtaPOI*dQ1nx))+(q2xEtaPOI*dQ2nx+q2yEtaPOI*dQ2ny)-2.*(dMult-3.)*(qxEtaPOI*dQ1nx+qyEtaPOI*dQ1ny)-2.*mEtaPOI*(dQ1nx*dQ1nx+dQ1ny*dQ1ny)+2.*(dQ1nx*qxEtaPOI+dQ1ny*qyEtaPOI)+2.*mEtaPOI*(dMult-3.))/(mEtaPOI*(dMult-1.)*(dMult-2.)*(dMult-3.));//to be improved (correct formula)
   //f4PerEtaBin1n1n1n1nPOI->Fill(fEtaMin+(bin-1)*dBinWidthEta,fourDiffEta1n1n1n1nPOI,mEtaPOI*(dMult-1.)*(dMult-2.)*(dMult-3.));//Re[<4'>_{n,n|n,n}]
  }
   
 } 
  
 fEtaReq1nPOI->Reset();
 fEtaImq1nPOI->Reset();
 fEtaReq2nPOI->Reset();
 fEtaImq2nPOI->Reset();
 
//---------------------------------------------------------------------------------------------------------
























/*








 //-------------------------------------------------------------------------------------------------------------------------------- 
 //
 //                                          **********************
 //                                          **** NESTED LOOPS ****
 //                                          ********************** 
 //
 // Remark 1: multi-particle correlations calculated with nested loops are stored in fDirectCorrelations.
 // Remark 2: binning of fDirectCorrelations: bins 0..40 - correlations needed for integrated flow; bins 40..80 - correlations needed for differential flow (taking as an example bin 0.5 < pt < 0.6)
 //
 // binning details of fDirectCorrelations (integrated flow):
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
 
 Double_t phi1=0., phi2=0., phi3=0., phi4=0., phi5=0., phi6=0., phi7=0., phi8=0.;

 //<2>_{k*n|k*n} (k=1,2,3 and 4)
 for(Int_t i1=0;i1<dMult;i1++)
 {
  fTrack=anEvent->GetTrack(i1);
  phi1=fTrack->Phi();
  for(Int_t i2=0;i2<dMult;i2++)
  {
   if(i2==i1)continue;
   fTrack=anEvent->GetTrack(i2);
   phi2=fTrack->Phi();
   fDirectCorrelations->Fill(0.,cos(n*(phi1-phi2)),1);    //<2>_{n|n}
   fDirectCorrelations->Fill(1.,cos(2.*n*(phi1-phi2)),1); //<2>_{2n|2n}
   fDirectCorrelations->Fill(2.,cos(3.*n*(phi1-phi2)),1); //<2>_{3n|3n}
   fDirectCorrelations->Fill(3.,cos(4.*n*(phi1-phi2)),1); //<2>_{4n|4n} 
  }
 }  
     
 //<3>_{2n|n,n}, <3>_{3n|2n,n}, <3>_{4n|2n,2n} and <3>_{4n|3n,n}
 for(Int_t i1=0;i1<dMult;i1++)
 {
  fTrack=anEvent->GetTrack(i1);
  phi1=fTrack->Phi();
  for(Int_t i2=0;i2<dMult;i2++)
  {
   if(i2==i1)continue;
   fTrack=anEvent->GetTrack(i2);
   phi2=fTrack->Phi();
   for(Int_t i3=0;i3<dMult;i3++)
   {
    if(i3==i1||i3==i2)continue;
    fTrack=anEvent->GetTrack(i3);
    phi3=fTrack->Phi();
    fDirectCorrelations->Fill(5.,cos(2*n*phi1-n*(phi2+phi3)),1);        //<3>_{2n|n,n}
    fDirectCorrelations->Fill(6.,cos(3.*n*phi1-2.*n*phi2-n*phi3),1);    //<3>_{3n|2n,n}
    fDirectCorrelations->Fill(7.,cos(4.*n*phi1-2.*n*phi2-2.*n*phi3),1); //<3>_{4n|2n,2n}
    fDirectCorrelations->Fill(8.,cos(4.*n*phi1-3.*n*phi2-n*phi3),1);    //<3>_{4n|3n,n}
   }
  }
 }
  
 //<4>_{n,n|n,n}, <4>_{2n,n|2n,n}, <4>_{2n,2n|2n,2n}, <4>_{3n|n,n,n}, <4>_{3n,n|3n,n}, <4>_{3n,n|2n,2n} and <4>_{4n|2n,n,n} 
 for(Int_t i1=0;i1<dMult;i1++)
 {
  fTrack=anEvent->GetTrack(i1);
  phi1=fTrack->Phi();
  for(Int_t i2=0;i2<dMult;i2++)
  {
   if(i2==i1)continue;
   fTrack=anEvent->GetTrack(i2);
   phi2=fTrack->Phi();
   for(Int_t i3=0;i3<dMult;i3++)
   {
    if(i3==i1||i3==i2)continue;
    fTrack=anEvent->GetTrack(i3);
    phi3=fTrack->Phi();
    for(Int_t i4=0;i4<dMult;i4++)
    {
     if(i4==i1||i4==i2||i4==i3)continue;
     fTrack=anEvent->GetTrack(i4);
     phi4=fTrack->Phi();
     fDirectCorrelations->Fill(10.,cos(n*phi1+n*phi2-n*phi3-n*phi4),1);            //<4>_{n,n|n,n}
     fDirectCorrelations->Fill(11.,cos(2.*n*phi1+n*phi2-2.*n*phi3-n*phi4),1);      //<4>_{2n,n|2n,n}
     fDirectCorrelations->Fill(12.,cos(2.*n*phi1+2*n*phi2-2.*n*phi3-2.*n*phi4),1); //<4>_{2n,2n|2n,2n}
     fDirectCorrelations->Fill(13.,cos(3.*n*phi1-n*phi2-n*phi3-n*phi4),1);         //<4>_{3n|n,n,n}
     fDirectCorrelations->Fill(14.,cos(3.*n*phi1+n*phi2-3.*n*phi3-n*phi4),1);      //<4>_{3n,n|3n,n}   
     fDirectCorrelations->Fill(15.,cos(3.*n*phi1+n*phi2-2.*n*phi3-2.*n*phi4),1);   //<4>_{3n,n|2n,2n}
     fDirectCorrelations->Fill(16.,cos(4.*n*phi1-2.*n*phi2-n*phi3-n*phi4),1);      //<4>_{4n|2n,n,n}
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
  for(Int_t i2=0;i2<dMult;i2++)
  {
   if(i2==i1)continue;
   fTrack=anEvent->GetTrack(i2);
   phi2=fTrack->Phi();
   for(Int_t i3=0;i3<dMult;i3++)
   {
    if(i3==i1||i3==i2)continue;
    fTrack=anEvent->GetTrack(i3);
    phi3=fTrack->Phi();
    for(Int_t i4=0;i4<dMult;i4++)
    {
     if(i4==i1||i4==i2||i4==i3)continue;
     fTrack=anEvent->GetTrack(i4);
     phi4=fTrack->Phi();
     for(Int_t i5=0;i5<dMult;i5++)
     {
      if(i5==i1||i5==i2||i5==i3||i5==i4)continue;
      fTrack=anEvent->GetTrack(i5);
      phi5=fTrack->Phi();
      fDirectCorrelations->Fill(18.,cos(2.*n*phi1+n*phi2-n*phi3-n*phi4-n*phi5),1);       //<5>_{2n,n|n,n,n}
      fDirectCorrelations->Fill(19.,cos(2.*n*phi1+2.*n*phi2-2.*n*phi3-n*phi4-n*phi5),1); //<5>_{2n,2n|2n,n,n}
      fDirectCorrelations->Fill(20.,cos(3.*n*phi1+n*phi2-2.*n*phi3-n*phi4-n*phi5),1);    //<5>_{3n,n|2n,n,n}
      fDirectCorrelations->Fill(21.,cos(4.*n*phi1-n*phi2-n*phi3-n*phi4-n*phi5),1);       //<5>_{4n|n,n,n,n}
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
  for(Int_t i2=0;i2<dMult;i2++)
  {
   if(i2==i1)continue;
   fTrack=anEvent->GetTrack(i2);
   phi2=fTrack->Phi();
   for(Int_t i3=0;i3<dMult;i3++)
   {
    if(i3==i1||i3==i2)continue;
    fTrack=anEvent->GetTrack(i3);
    phi3=fTrack->Phi();
    for(Int_t i4=0;i4<dMult;i4++)
    {
     if(i4==i1||i4==i2||i4==i3)continue;
     fTrack=anEvent->GetTrack(i4);
     phi4=fTrack->Phi();
     for(Int_t i5=0;i5<dMult;i5++)
     {
      if(i5==i1||i5==i2||i5==i3||i5==i4)continue;
      fTrack=anEvent->GetTrack(i5);
      phi5=fTrack->Phi();
      for(Int_t i6=0;i6<dMult;i6++)
      {
       if(i6==i1||i6==i2||i6==i3||i6==i4||i6==i5)continue;
       fTrack=anEvent->GetTrack(i6);
       phi6=fTrack->Phi(); 
       fDirectCorrelations->Fill(23.,cos(n*phi1+n*phi2+n*phi3-n*phi4-n*phi5-n*phi6),1);       //<6>_{n,n,n|n,n,n}
       fDirectCorrelations->Fill(24.,cos(2.*n*phi1+n*phi2+n*phi3-2.*n*phi4-n*phi5-n*phi6),1); //<6>_{2n,n,n|2n,n,n}
       fDirectCorrelations->Fill(25.,cos(2.*n*phi1+2.*n*phi2-n*phi3-n*phi4-n*phi5-n*phi6),1); //<6>_{2n,2n|n,n,n,n}
       fDirectCorrelations->Fill(26.,cos(3.*n*phi1+n*phi2-n*phi3-n*phi4-n*phi5-n*phi6),1);    //<6>_{3n,n|n,n,n,n}     
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
  for(Int_t i2=0;i2<dMult;i2++)
  {
   if(i2==i1)continue;
   fTrack=anEvent->GetTrack(i2);
   phi2=fTrack->Phi();
   for(Int_t i3=0;i3<dMult;i3++)
   {
    if(i3==i1||i3==i2)continue;
    fTrack=anEvent->GetTrack(i3);
    phi3=fTrack->Phi();
    for(Int_t i4=0;i4<dMult;i4++)
    {
     if(i4==i1||i4==i2||i4==i3)continue;
     fTrack=anEvent->GetTrack(i4);
     phi4=fTrack->Phi();
     for(Int_t i5=0;i5<dMult;i5++)
     {
      if(i5==i1||i5==i2||i5==i3||i5==i4)continue;
      fTrack=anEvent->GetTrack(i5);
      phi5=fTrack->Phi();
      for(Int_t i6=0;i6<dMult;i6++)
      {
       if(i6==i1||i6==i2||i6==i3||i6==i4||i6==i5)continue;
       fTrack=anEvent->GetTrack(i6);
       phi6=fTrack->Phi(); 
       for(Int_t i7=0;i7<dMult;i7++)
       {
        if(i7==i1||i7==i2||i7==i3||i7==i4||i7==i5||i7==i6)continue;
        fTrack=anEvent->GetTrack(i7);
        phi7=fTrack->Phi(); 
        fDirectCorrelations->Fill(28.,cos(2.*n*phi1+n*phi2+n*phi3-n*phi4-n*phi5-n*phi6-n*phi7),1);//<7>_{2n,n,n|n,n,n,n}
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
  for(Int_t i2=0;i2<dMult;i2++)
  {
   if(i2==i1)continue;
   fTrack=anEvent->GetTrack(i2);
   phi2=fTrack->Phi();
   for(Int_t i3=0;i3<dMult;i3++)
   {
    if(i3==i1||i3==i2)continue;
    fTrack=anEvent->GetTrack(i3);
    phi3=fTrack->Phi();
    for(Int_t i4=0;i4<dMult;i4++)
    {
     if(i4==i1||i4==i2||i4==i3)continue;
     fTrack=anEvent->GetTrack(i4);
     phi4=fTrack->Phi();
     for(Int_t i5=0;i5<dMult;i5++)
     {
      if(i5==i1||i5==i2||i5==i3||i5==i4)continue;
      fTrack=anEvent->GetTrack(i5);
      phi5=fTrack->Phi();
      for(Int_t i6=0;i6<dMult;i6++)
      {
       if(i6==i1||i6==i2||i6==i3||i6==i4||i6==i5)continue;
       fTrack=anEvent->GetTrack(i6);
       phi6=fTrack->Phi(); 
       for(Int_t i7=0;i7<dMult;i7++)
       {
        if(i7==i1||i7==i2||i7==i3||i7==i4||i7==i5||i7==i6)continue;
        fTrack=anEvent->GetTrack(i7);
        phi7=fTrack->Phi(); 
        for(Int_t i8=0;i8<dMult;i8++)
        {
         if(i8==i1||i8==i2||i8==i3||i8==i4||i8==i5||i8==i6||i8==i7)continue;
         fTrack=anEvent->GetTrack(i8);
         phi8=fTrack->Phi();  
         fDirectCorrelations->Fill(30.,cos(n*phi1+n*phi2+n*phi3+n*phi4-n*phi5-n*phi6-n*phi7-n*phi8),1);//<8>_{n,n,n,n|n,n,n,n}
        } 
       } 
      } 
     }
    }  
   }
  }
 }
 
 // binning details of fDirectCorrelations (differential flow):
 //
 //41st bin: <2'>_{n|n}
 //42nd bin: <2'>_{2n|2n}
 //46th bin: <3'>_{2n|n,n}
 //47th bin: <3'>_{n,n|2n}
 //51st bin: <4'>_{n,n|n,n}
 
 //<2'>_{n|n}
 for(Int_t i1=0;i1<dMult;i1++)
 {
  fTrack=anEvent->GetTrack(i1);
  if(fTrack->Pt()>=0.5&&fTrack->Pt()<0.6)
  {
   phi1=fTrack->Phi();
   for(Int_t i2=0;i2<dMult;i2++)
   {
    if(i2==i1)continue;
    fTrack=anEvent->GetTrack(i2);
    phi2=fTrack->Phi(); 
    fDirectCorrelations->Fill(40.,cos(1.*n*(phi1-phi2)),1); //<2'>_{n,n}
    fDirectCorrelations->Fill(41.,cos(2.*n*(phi1-phi2)),1); //<2'>_{2n,2n}  
   }
  }
 }  

 //<3'>_{2n|n,n}
 for(Int_t i1=0;i1<dMult;i1++)
 {
  fTrack=anEvent->GetTrack(i1);
  if(fTrack->Pt()>=0.5&&fTrack->Pt()<0.6)
  {
   phi1=fTrack->Phi();
   for(Int_t i2=0;i2<dMult;i2++)
   {
    if(i2==i1)continue;
    fTrack=anEvent->GetTrack(i2);
    phi2=fTrack->Phi();
    for(Int_t i3=0;i3<dMult;i3++)
    {
     if(i3==i1||i3==i2)continue;
     fTrack=anEvent->GetTrack(i3);
     phi3=fTrack->Phi();           
     fDirectCorrelations->Fill(45.,cos(n*(2.*phi1-phi2-phi3)),1); //<3'>_{2n|n,n}
     fDirectCorrelations->Fill(46.,cos(n*(phi1+phi2-2.*phi3)),1); //<3'>_{n,n|2n}    
    }
   }  
  }  
 }
  
 //<4'>_{n,n|n,n}
 for(Int_t i1=0;i1<dMult;i1++)
 {
  fTrack=anEvent->GetTrack(i1);
  if(fTrack->Pt()>=0.5&&fTrack->Pt()<0.6)
  {
   phi1=fTrack->Phi();
   for(Int_t i2=0;i2<dMult;i2++)
   {
    if(i2==i1)continue;
    fTrack=anEvent->GetTrack(i2);
    phi2=fTrack->Phi();
    for(Int_t i3=0;i3<dMult;i3++)
    {
     if(i3==i1||i3==i2)continue;
     fTrack=anEvent->GetTrack(i3);
     phi3=fTrack->Phi();
     for(Int_t i4=0;i4<dMult;i4++)
     {
      if(i4==i1||i4==i2||i4==i3)continue;
      fTrack=anEvent->GetTrack(i4);
      phi4=fTrack->Phi();
      fDirectCorrelations->Fill(50.,cos(n*(phi1+phi2-phi3-phi4)),1); //<4'>_{n,n|n,n}   
     } 
    }
   }  
  }  
 }
 //--------------------------------------------------------------------------------------------------------------------------------  




*/




//}//line needed only for nested loops - end of if(nPrim>8&&nPrim<14) 

}//end of Make()

//================================================================================================================

void AliFlowAnalysisWithQCumulants::Finish()
{
 //calculate the final results
 AliQCumulantsFunctions finalResults(fIntFlowResultsQC,fDiffFlowResults2ndOrderQC,fDiffFlowResults4thOrderQC,fCovariances,fAvMultIntFlowQC,fQvectorComponents,fQCorrelations, fQProduct,fDirectCorrelations,f2PerPtBin1n1nRP,f2PerPtBin2n2nRP,f3PerPtBin2n1n1nRP,f3PerPtBin1n1n2nRP,f4PerPtBin1n1n1n1nRP, f2PerEtaBin1n1nRP,f2PerEtaBin2n2nRP,f3PerEtaBin2n1n1nRP,f3PerEtaBin1n1n2nRP,f4PerEtaBin1n1n1n1nRP, f2PerPtBin1n1nPOI,f2PerPtBin2n2nPOI,f3PerPtBin2n1n1nPOI,f3PerPtBin1n1n2nPOI,f4PerPtBin1n1n1n1nPOI, f2PerEtaBin1n1nPOI,f2PerEtaBin2n2nPOI,f3PerEtaBin2n1n1nPOI,f3PerEtaBin1n1n2nPOI,f4PerEtaBin1n1n1n1nPOI,fCommonHistsResults2nd, fCommonHistsResults4th,fCommonHistsResults6th,fCommonHistsResults8th);
         
 finalResults.Calculate();  
}

//================================================================================================================

void AliFlowAnalysisWithQCumulants::WriteHistograms(TString* outputFileName)
{
 //store the final results in output .root file
 TFile *output = new TFile(outputFileName->Data(),"RECREATE");
 output->WriteObject(fHistList, "cobjQC","SingleKey");
 delete output;
}

//================================================================================================================

