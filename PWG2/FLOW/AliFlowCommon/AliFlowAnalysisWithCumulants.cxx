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

/******************************** 
 * flow analysis with cumulants * 
 *                              *
 * author: Ante Bilandzic       * 
 *          (anteb@nikhef.nl)   *
 *******************************/ 

#define AliFlowAnalysisWithCumulants_cxx

#include "Riostream.h"
#include "AliFlowCommonConstants.h"
#include "AliFlowCommonHist.h"
#include "AliFlowCommonHistResults.h"
#include "TChain.h"
#include "TFile.h"
#include "TList.h" //NEW
#include "TParticle.h"
#include "TRandom3.h"
#include "TProfile.h"
#include "TProfile2D.h" 
#include "TProfile3D.h"
#include "AliFlowEventSimple.h"
#include "AliFlowTrackSimple.h"
#include "AliFlowAnalysisWithCumulants.h"
#include "AliFlowCumuConstants.h"
#include "AliCumulantsFunctions.h"

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

ClassImp(AliFlowAnalysisWithCumulants)

AliFlowAnalysisWithCumulants::AliFlowAnalysisWithCumulants():  
 fTrack(NULL),
 fHistList(NULL),
 fWeightsList(NULL),
 fR0(0),
 fPtMax(0),
 fPtMin(0),
 fBinWidthPt(0),
 fgknBinsPt(0),
 fEtaMax(0),
 fEtaMin(0),
 fBinWidthEta(0),
 fgknBinsEta(0),
 fAvQx(0),
 fAvQy(0),
 fAvQ2x(0),
 fAvQ2y(0),
 fAvMultIntFlowGFC(NULL),
 fQVectorComponentsGFC(NULL),
 fIntFlowResultsGFC(NULL),
 fDiffFlowResults2ndOrderGFC(NULL),
 fDiffFlowResults4thOrderGFC(NULL),
 fDiffFlowResults6thOrderGFC(NULL),
 fDiffFlowResults8thOrderGFC(NULL),
 fCommonHistsResults2nd(NULL),
 fCommonHistsResults4th(NULL),
 fCommonHistsResults6th(NULL),
 fCommonHistsResults8th(NULL),
 fIntFlowGenFun(NULL),
 fIntFlowGenFun4(NULL),//(only for other system of Eq.)
 fIntFlowGenFun6(NULL),//(only for other system of Eq.)
 fIntFlowGenFun8(NULL),//(only for other system of Eq.)
 fIntFlowGenFun16(NULL),//(only for other system of Eq.)
 fAvMultIntFlow4GFC(NULL),//(only for other system of Eq.)
 fAvMultIntFlow6GFC(NULL),//(only for other system of Eq.)
 fAvMultIntFlow8GFC(NULL),//(only for other system of Eq.)
 fAvMultIntFlow16GFC(NULL),//(only for other system of Eq.)
 fDiffFlowPtRPGenFunRe(NULL),
 fDiffFlowPtRPGenFunIm(NULL),
 fPtBinRPNoOfParticles(NULL),
 fDiffFlowEtaRPGenFunRe(NULL),
 fDiffFlowEtaRPGenFunIm(NULL),
 fEtaBinRPNoOfParticles(NULL),
 fDiffFlowPtPOIGenFunRe(NULL),
 fDiffFlowPtPOIGenFunIm(NULL),
 fPtBinPOINoOfParticles(NULL),
 fDiffFlowEtaPOIGenFunRe(NULL),
 fDiffFlowEtaPOIGenFunIm(NULL),
 fEtaBinPOINoOfParticles(NULL),
 /*
 fDiffFlowGenFunRe0(NULL),
 fDiffFlowGenFunRe1(NULL),
 fDiffFlowGenFunRe2(NULL),
 fDiffFlowGenFunRe3(NULL),
 fDiffFlowGenFunRe4(NULL),
 fDiffFlowGenFunRe5(NULL),
 fDiffFlowGenFunRe6(NULL),
 fDiffFlowGenFunRe7(NULL),
 fDiffFlowGenFunIm0(NULL),
 fDiffFlowGenFunIm1(NULL),
 fDiffFlowGenFunIm2(NULL),
 fDiffFlowGenFunIm3(NULL),
 fDiffFlowGenFunIm4(NULL),
 fDiffFlowGenFunIm5(NULL),
 fDiffFlowGenFunIm6(NULL),
 fDiffFlowGenFunIm7(NULL),
 */
 fCommonHists(NULL),
 fOtherEquations(kFALSE),
 fUsePhiWeights(kFALSE),
 fUsePtWeights(kFALSE),
 fUseEtaWeights(kFALSE),
 fAverageOfSquaredWeight(NULL)
{
 //constructor 
 fHistList = new TList();
 fWeightsList = new TList();
 fWeightsList->SetName("Weights");
 fR0=AliFlowCumuConstants::fgR0;
 //Pt:
 fPtMax=AliFlowCommonConstants::GetPtMax(); 
 fPtMin=AliFlowCommonConstants::GetPtMin();
 fgknBinsPt=AliFlowCommonConstants::GetNbinsPt();
 if(fgknBinsPt)
 {
  fBinWidthPt=(fPtMax-fPtMin)/fgknBinsPt;  
 } 
 //Eta: 
 fEtaMax=AliFlowCommonConstants::GetEtaMax(); 
 fEtaMin=AliFlowCommonConstants::GetEtaMin();
 fgknBinsEta=AliFlowCommonConstants::GetNbinsEta();
 if(fgknBinsEta)
 {
  fBinWidthEta=(fEtaMax-fEtaMin)/fgknBinsEta;  
 } 
 
 fOtherEquations=AliFlowCumuConstants::fgOtherEquations;
}

AliFlowAnalysisWithCumulants::~AliFlowAnalysisWithCumulants()
{
 //desctructor
 delete fHistList; 
 delete fWeightsList;
}

//================================================================================================================

void AliFlowAnalysisWithCumulants::Init()
{
 //various output histograms
 
 //average multiplicity 
 fAvMultIntFlowGFC = new TProfile("fAvMultIntFlowGFC","Average Weighted Multiplicity",1,0,1,"s");
 fAvMultIntFlowGFC->SetXTitle("");
 fAvMultIntFlowGFC->SetYTitle("");
 fAvMultIntFlowGFC->SetLabelSize(0.06);
 fAvMultIntFlowGFC->SetMarkerStyle(25);
 fAvMultIntFlowGFC->SetLabelOffset(0.01);
 (fAvMultIntFlowGFC->GetXaxis())->SetBinLabel(1,"Average Weighted Multiplicity");
 fHistList->Add(fAvMultIntFlowGFC);
 
 //averages of Q-vector components (1st bin: <Q_x>, 2nd bin: <Q_y>, 3rd bin: <(Q_x)^2>, 4th bin: <(Q_y)^2>)
 fQVectorComponentsGFC = new TProfile("fQVectorComponentsGFC","Average of Q-vector components",4,0.,4.);
 fQVectorComponentsGFC->SetXTitle("");
 fQVectorComponentsGFC->SetYTitle("");
 fQVectorComponentsGFC->SetLabelSize(0.06);
 fQVectorComponentsGFC->SetMarkerStyle(25);
 (fQVectorComponentsGFC->GetXaxis())->SetBinLabel(1,"<Q_{x}>");
 (fQVectorComponentsGFC->GetXaxis())->SetBinLabel(2,"<Q_{y}>");
 (fQVectorComponentsGFC->GetXaxis())->SetBinLabel(3,"<Q_{x}^{2}>");
 (fQVectorComponentsGFC->GetXaxis())->SetBinLabel(4,"<Q_{y}^{2}>");
 fHistList->Add(fQVectorComponentsGFC);
 
 //final results for integrated flow (v_n{2}, v_n{4},..., v_n{16}) from cumulants (by default n=2) 
 fIntFlowResultsGFC = new TH1D("fIntFlowResultsGFC","Integrated Flow From Cumulants (Generating Function)",4,0,4);
 fIntFlowResultsGFC->SetXTitle("");
 fIntFlowResultsGFC->SetYTitle("");
 fIntFlowResultsGFC->SetLabelSize(0.06);
 //fIntFlowResultsGFC->SetTickLength(1);
 fIntFlowResultsGFC->SetMarkerStyle(25);
 (fIntFlowResultsGFC->GetXaxis())->SetBinLabel(1,"v_{n}{2}");
 (fIntFlowResultsGFC->GetXaxis())->SetBinLabel(2,"v_{n}{4}");
 (fIntFlowResultsGFC->GetXaxis())->SetBinLabel(3,"v_{n}{6}");
 (fIntFlowResultsGFC->GetXaxis())->SetBinLabel(4,"v_{n}{8}");
 fHistList->Add(fIntFlowResultsGFC);
  
 //final results for differential flow v_p/n{2} (by default p=n=2)
 fDiffFlowResults2ndOrderGFC = new TH1D("fDiffFlowResults2ndOrderGFC","v'_2/2{2}",fgknBinsPt,fPtMin,fPtMax);
 fDiffFlowResults2ndOrderGFC->SetXTitle("pt [GeV]");
 fDiffFlowResults2ndOrderGFC->SetYTitle("");
 fHistList->Add(fDiffFlowResults2ndOrderGFC);

 //final results for differential flow v_p/n{4} (by default p=n=2) 
 fDiffFlowResults4thOrderGFC = new TH1D("fDiffFlowResults4thOrderGFC","v'_2/2{4}",fgknBinsPt,fPtMin,fPtMax);
 fDiffFlowResults4thOrderGFC->SetXTitle("pt [GeV]");
 fDiffFlowResults4thOrderGFC->SetYTitle("");
 fHistList->Add(fDiffFlowResults4thOrderGFC);
 
 //final results for differential flow v_p/n{6} (by default p=n=2)  
 fDiffFlowResults6thOrderGFC = new TH1D("fDiffFlowResults6thOrderGFC","v'_2/2{6}",fgknBinsPt,fPtMin,fPtMax);
 fDiffFlowResults6thOrderGFC->SetXTitle("pt [GeV]");
 fDiffFlowResults6thOrderGFC->SetYTitle("");
 fHistList->Add(fDiffFlowResults6thOrderGFC);
 
 //final results for differential flow v_p/n{8} (by default p=n=2)
 fDiffFlowResults8thOrderGFC = new TH1D("fDiffFlowResults8thOrderGFC","v'_2/2{8}",fgknBinsPt,fPtMin,fPtMax);
 fDiffFlowResults8thOrderGFC->SetXTitle("pt [GeV]");
 fDiffFlowResults8thOrderGFC->SetYTitle("");
 fHistList->Add(fDiffFlowResults8thOrderGFC);
  
 //avarage of the generating function for integrated flow <G[p][q]>
 fIntFlowGenFun = new TProfile2D("fIntFlowGenFun","<G[p][q]>",fgkPmax,0.,(Double_t)fgkPmax,fgkQmax,0.,(Double_t)fgkQmax);
 fIntFlowGenFun->SetXTitle("p");
 fIntFlowGenFun->SetYTitle("q");
 fHistList->Add(fIntFlowGenFun);

 if(fOtherEquations)
 {
  //avarage of the generating function for integrated flow <G[p][q]> (only for other system of Eq. - up to 4th order)
  fIntFlowGenFun4 = new TProfile2D("fIntFlowGenFun4","<G4[p4][q4]>",fgkPmax4,0.,(Double_t)fgkPmax4,fgkQmax4,0.,(Double_t)fgkQmax4);
  fIntFlowGenFun4->SetXTitle("p4");
  fIntFlowGenFun4->SetYTitle("q4");
  fHistList->Add(fIntFlowGenFun4);

  //avarage of the generating function for integrated flow <G[p][q]> (only for other system of Eq. - up to 6th order) 
  fIntFlowGenFun6 = new TProfile2D("fIntFlowGenFun6","<G6[p6][q6]>",fgkPmax6,0.,(Double_t)fgkPmax6,fgkQmax6,0.,(Double_t)fgkQmax6);
  fIntFlowGenFun6->SetXTitle("p6");
  fIntFlowGenFun6->SetYTitle("q6");
  fHistList->Add(fIntFlowGenFun6);

  //avarage of the generating function for integrated flow <G[p][q]> (only for other system of Eq. - up to 8th order)
  fIntFlowGenFun8 = new TProfile2D("fIntFlowGenFun8","<G8[p8][q8]>",fgkPmax8,0.,(Double_t)fgkPmax8,fgkQmax8,0.,(Double_t)fgkQmax8);
  fIntFlowGenFun8->SetXTitle("p8");
  fIntFlowGenFun8->SetYTitle("q8");
  fHistList->Add(fIntFlowGenFun8);

  //avarage of the generating function for integrated flow <G[p][q]> (only for other system of Eq. - up to 16th order)
  fIntFlowGenFun16 = new TProfile2D("fIntFlowGenFun16","<G16[p16][q16]>",fgkPmax16,0.,(Double_t)fgkPmax16,fgkQmax16,0.,(Double_t)fgkQmax16);
  fIntFlowGenFun16->SetXTitle("p16");
  fIntFlowGenFun16->SetYTitle("q16");
  fHistList->Add(fIntFlowGenFun16);
 
  //average multiplicity (only for other system of Eq. - up to 4th order)
  fAvMultIntFlow4GFC = new TProfile("fAvMultIntFlow4GFC","Average Multiplicity",1,0,1,"s");
  fAvMultIntFlow4GFC->SetXTitle("");
  fAvMultIntFlow4GFC->SetYTitle("");
  fAvMultIntFlow4GFC->SetLabelSize(0.06);
  fAvMultIntFlow4GFC->SetMarkerStyle(25);
  fAvMultIntFlow4GFC->SetLabelOffset(0.01);
  (fAvMultIntFlow4GFC->GetXaxis())->SetBinLabel(1,"Average Multiplicity");
  fHistList->Add(fAvMultIntFlow4GFC);
 
  //average multiplicity (only for other system of Eq. - up to 6th order)
  fAvMultIntFlow6GFC = new TProfile("fAvMultIntFlow6GFC","Average Multiplicity",1,0,1,"s");
  fAvMultIntFlow6GFC->SetXTitle("");
  fAvMultIntFlow6GFC->SetYTitle("");
  fAvMultIntFlow6GFC->SetLabelSize(0.06);
  fAvMultIntFlow6GFC->SetMarkerStyle(25);
  fAvMultIntFlow6GFC->SetLabelOffset(0.01);
  (fAvMultIntFlow6GFC->GetXaxis())->SetBinLabel(1,"Average Multiplicity");
  fHistList->Add(fAvMultIntFlow6GFC);
 
  //average multiplicity (only for other system of Eq. - up to 8th order)
  fAvMultIntFlow8GFC = new TProfile("fAvMultIntFlow8GFC","Average Multiplicity",1,0,1,"s");
  fAvMultIntFlow8GFC->SetXTitle("");
  fAvMultIntFlow8GFC->SetYTitle("");
  fAvMultIntFlow8GFC->SetLabelSize(0.06);
  fAvMultIntFlow8GFC->SetMarkerStyle(25);
  fAvMultIntFlow8GFC->SetLabelOffset(0.01);
  (fAvMultIntFlow8GFC->GetXaxis())->SetBinLabel(1,"Average Multiplicity");
  fHistList->Add(fAvMultIntFlow8GFC);
 
  //average multiplicity (only for other system of Eq. - up to 16th order)
  fAvMultIntFlow16GFC = new TProfile("fAvMultIntFlow16GFC","Average Multiplicity",1,0,1,"s");
  fAvMultIntFlow16GFC->SetXTitle("");
  fAvMultIntFlow16GFC->SetYTitle("");
  fAvMultIntFlow16GFC->SetLabelSize(0.06);
  fAvMultIntFlow16GFC->SetMarkerStyle(25);
  fAvMultIntFlow16GFC->SetLabelOffset(0.01);
  (fAvMultIntFlow16GFC->GetXaxis())->SetBinLabel(1,"Average Multiplicity");
  fHistList->Add(fAvMultIntFlow16GFC);
 }
 
 //avarage of the real part of generating function for differential flow in Pt <Re(D[b][p][q])>
 fDiffFlowPtRPGenFunRe = new TProfile3D("fDiffFlowPtRPGenFunRe","<Re(D[b][p][q])>",fgknBinsPt,(Double_t)(fPtMin/fBinWidthPt),(Double_t)(fPtMax/fBinWidthPt),fgkPmax,0.,(Double_t)fgkPmax,fgkQmax,0.,(Double_t)fgkQmax);
 fDiffFlowPtRPGenFunRe->SetXTitle("b");
 fDiffFlowPtRPGenFunRe->SetYTitle("p");
 fDiffFlowPtRPGenFunRe->SetZTitle("q");
 fDiffFlowPtRPGenFunRe->SetTitleOffset(1.44,"X");
 fDiffFlowPtRPGenFunRe->SetTitleOffset(1.44,"Y");
 fHistList->Add(fDiffFlowPtRPGenFunRe);
 
 //avarage of the imaginary part of generating function for differential flow in Pt <Im(D[b][p][q])>
 fDiffFlowPtRPGenFunIm = new TProfile3D("fDiffFlowPtRPGenFunIm","<Im(D[b][p][q])>",fgknBinsPt,(Double_t)(fPtMin/fBinWidthPt),(Double_t)(fPtMax/fBinWidthPt),fgkPmax,0.,(Double_t)fgkPmax,fgkQmax,0.,(Double_t)fgkQmax);
 fDiffFlowPtRPGenFunIm->SetXTitle("b");
 fDiffFlowPtRPGenFunIm->SetYTitle("p");
 fDiffFlowPtRPGenFunIm->SetZTitle("q");
 fDiffFlowPtRPGenFunIm->SetTitleOffset(1.44,"X");
 fDiffFlowPtRPGenFunIm->SetTitleOffset(1.44,"Y");
 fHistList->Add(fDiffFlowPtRPGenFunIm);
 
 //number of particles per pt bin
 fPtBinRPNoOfParticles = new TProfile("fPtBinRPNoOfParticles","Number of particles per #p_{t} bin",fgknBinsPt,fPtMin,fPtMax);
 fPtBinRPNoOfParticles->SetXTitle("pt [GeV]");
 fPtBinRPNoOfParticles->SetYTitle("");
 fHistList->Add(fPtBinRPNoOfParticles);
 
 //avarage of the real part of generating function for differential flow in Eta <Re(D[b][p][q])>
 fDiffFlowEtaRPGenFunRe = new TProfile3D("fDiffFlowEtaRPGenFunRe","<Re(D[b][p][q])>",fgknBinsEta,(Double_t)(fEtaMin/fBinWidthEta),(Double_t)(fEtaMax/fBinWidthEta),fgkPmax,0.,(Double_t)fgkPmax,fgkQmax,0.,(Double_t)fgkQmax);
 fDiffFlowEtaRPGenFunRe->SetXTitle("b");
 fDiffFlowEtaRPGenFunRe->SetYTitle("p");
 fDiffFlowEtaRPGenFunRe->SetZTitle("q");
 fDiffFlowEtaRPGenFunRe->SetTitleOffset(1.44,"X");
 fDiffFlowEtaRPGenFunRe->SetTitleOffset(1.44,"Y");
 fHistList->Add(fDiffFlowEtaRPGenFunRe);
 
 //avarage of the imaginary part of generating function for differential flow in Eta <Im(D[b][p][q])>
 fDiffFlowEtaRPGenFunIm = new TProfile3D("fDiffFlowEtaRPGenFunIm","<Im(D[b][p][q])>",fgknBinsEta,(Double_t)(fEtaMin/fBinWidthEta),(Double_t)(fEtaMax/fBinWidthEta),fgkPmax,0.,(Double_t)fgkPmax,fgkQmax,0.,(Double_t)fgkQmax);
 fDiffFlowEtaRPGenFunIm->SetXTitle("b");
 fDiffFlowEtaRPGenFunIm->SetYTitle("p");
 fDiffFlowEtaRPGenFunIm->SetZTitle("q");
 fDiffFlowEtaRPGenFunIm->SetTitleOffset(1.44,"X");
 fDiffFlowEtaRPGenFunIm->SetTitleOffset(1.44,"Y");
 fHistList->Add(fDiffFlowEtaRPGenFunIm);
 
 //number of particles per eta bin
 fEtaBinRPNoOfParticles = new TProfile("fEtaBinRPNoOfParticles","Number of particles per #eta bin",fgknBinsEta,fEtaMin,fEtaMax);
 fEtaBinRPNoOfParticles->SetXTitle("#eta");
 fEtaBinRPNoOfParticles->SetYTitle("");
 fHistList->Add(fEtaBinRPNoOfParticles);
 
 //avarage of the real part of generating function for differential flow in Pt <Re(D[b][p][q])>
 fDiffFlowPtPOIGenFunRe = new TProfile3D("fDiffFlowPtPOIGenFunRe","<Re(D[b][p][q])>",fgknBinsPt,(Double_t)(fPtMin/fBinWidthPt),(Double_t)(fPtMax/fBinWidthPt),fgkPmax,0.,(Double_t)fgkPmax,fgkQmax,0.,(Double_t)fgkQmax);
 fDiffFlowPtPOIGenFunRe->SetXTitle("b");
 fDiffFlowPtPOIGenFunRe->SetYTitle("p");
 fDiffFlowPtPOIGenFunRe->SetZTitle("q");
 fDiffFlowPtPOIGenFunRe->SetTitleOffset(1.44,"X");
 fDiffFlowPtPOIGenFunRe->SetTitleOffset(1.44,"Y");
 fHistList->Add(fDiffFlowPtPOIGenFunRe);
 
 //avarage of the imaginary part of generating function for differential flow in Pt <Im(D[b][p][q])>
 fDiffFlowPtPOIGenFunIm = new TProfile3D("fDiffFlowPtPOIGenFunIm","<Im(D[b][p][q])>",fgknBinsPt,(Double_t)(fPtMin/fBinWidthPt),(Double_t)(fPtMax/fBinWidthPt),fgkPmax,0.,(Double_t)fgkPmax,fgkQmax,0.,(Double_t)fgkQmax);
 fDiffFlowPtPOIGenFunIm->SetXTitle("b");
 fDiffFlowPtPOIGenFunIm->SetYTitle("p");
 fDiffFlowPtPOIGenFunIm->SetZTitle("q");
 fDiffFlowPtPOIGenFunIm->SetTitleOffset(1.44,"X");
 fDiffFlowPtPOIGenFunIm->SetTitleOffset(1.44,"Y");
 fHistList->Add(fDiffFlowPtPOIGenFunIm);
 
 //number of particles per pt bin
 fPtBinPOINoOfParticles = new TProfile("fPtBinPOINoOfParticles","Number of particles per #p_{t} bin",fgknBinsPt,fPtMin,fPtMax);
 fPtBinPOINoOfParticles->SetXTitle("pt [GeV]");
 fPtBinPOINoOfParticles->SetYTitle("");
 fHistList->Add(fPtBinPOINoOfParticles);
 
 //avarage of the real part of generating function for differential flow in Eta <Re(D[b][p][q])>
 fDiffFlowEtaPOIGenFunRe = new TProfile3D("fDiffFlowEtaPOIGenFunRe","<Re(D[b][p][q])>",fgknBinsEta,(Double_t)(fEtaMin/fBinWidthEta),(Double_t)(fEtaMax/fBinWidthEta),fgkPmax,0.,(Double_t)fgkPmax,fgkQmax,0.,(Double_t)fgkQmax);
 fDiffFlowEtaPOIGenFunRe->SetXTitle("b");
 fDiffFlowEtaPOIGenFunRe->SetYTitle("p");
 fDiffFlowEtaPOIGenFunRe->SetZTitle("q");
 fDiffFlowEtaPOIGenFunRe->SetTitleOffset(1.44,"X");
 fDiffFlowEtaPOIGenFunRe->SetTitleOffset(1.44,"Y");
 fHistList->Add(fDiffFlowEtaPOIGenFunRe);
 
 //avarage of the imaginary part of generating function for differential flow in Eta <Im(D[b][p][q])>
 fDiffFlowEtaPOIGenFunIm = new TProfile3D("fDiffFlowEtaPOIGenFunIm","<Im(D[b][p][q])>",fgknBinsEta,(Double_t)(fEtaMin/fBinWidthEta),(Double_t)(fEtaMax/fBinWidthEta),fgkPmax,0.,(Double_t)fgkPmax,fgkQmax,0.,(Double_t)fgkQmax);
 fDiffFlowEtaPOIGenFunIm->SetXTitle("b");
 fDiffFlowEtaPOIGenFunIm->SetYTitle("p");
 fDiffFlowEtaPOIGenFunIm->SetZTitle("q");
 fDiffFlowEtaPOIGenFunIm->SetTitleOffset(1.44,"X");
 fDiffFlowEtaPOIGenFunIm->SetTitleOffset(1.44,"Y");
 fHistList->Add(fDiffFlowEtaPOIGenFunIm);
 
 //number of particles per eta bin
 fEtaBinPOINoOfParticles = new TProfile("fEtaBinPOINoOfParticles","Number of particles per #eta bin",fgknBinsEta,fEtaMin,fEtaMax);
 fEtaBinPOINoOfParticles->SetXTitle("#eta");
 fEtaBinPOINoOfParticles->SetYTitle("");
 fHistList->Add(fEtaBinPOINoOfParticles);
 
 /*
 fDiffFlowGenFunRe0 = new TProfile2D("fDiffFlowGenFunRe0","Re(<D[b][0][q]>)",fgknBins,(Double_t)(fPtMin/fBinWidth),(Double_t)(fPtMax/fBinWidth),fgkQmax,0.,(Double_t)fgkQmax);
 fDiffFlowGenFunRe0->SetXTitle("b");
 fDiffFlowGenFunRe0->SetYTitle("q");
 fHistList->Add(fDiffFlowGenFunRe0);

 fDiffFlowGenFunIm0 = new TProfile2D("fDiffFlowGcout<<"HEY M1"<<endl;enFunIm0","Im(<D[b][0][q]>)",fgknBins,(Double_t)(fPtMin/fBinWidth),(Double_t)(fPtMax/fBinWidth),fgkQmax,0.,(Double_t)fgkQmax);
 fDiffFlowGenFunIm0->SetXTitle("b");
 fDiffFlowGenFunIm0->SetYTitle("q");
 fHistList->Add(fDiffFlowGenFunIm0);
 
 fDiffFlowGenFunRe1 = new TProfile2D("fDiffFlowGenFunRe1","Re(<D[b][1][q]>)",fgknBins,(Double_t)(fPtMin/fBinWidth),(Double_t)(fPtMax/fBinWidth),fgkQmax,0.,(Double_t)fgkQmax);
 fDiffFlowGenFunRe1->SetXTitle("b");
 fDiffFlowGenFunRe1->SetYTitle("q");
 fHistList->Add(fDiffFlowGenFunRe1);
 
 fDiffFlowGenFunIm1 = new TProfile2D("fDiffFlowGenFunIm1","Im(<D[b][1][q]>)",fgknBins,(Double_t)(fPtMin/fBinWidth),(Double_t)(fPtMax/fBinWidth),fgkQmax,0.,(Double_t)fgkQmax);
 fDiffFlowGenFunIm1->SetXTitle("b");
 fDiffFlowGenFunIm1->SetYTitle("q");
 fHistList->Add(fDiffFlowGenFunIm1);
 
 fDiffFlowGenFunRe2 = new TProfile2D("fDiffFlowGenFunRe2","Re(<D[b][2][q]>)",fgknBins,(Double_t)(fPtMin/fBinWidth),(Double_t)(fPtMax/fBinWidth),fgkQmax,0.,(Double_t)fgkQmax);
 fDiffFlowGenFunRe2->SetXTitle("b");
 fDiffFlowGenFunRe2->SetYTitle("q");
 fHistList->Add(fDiffFlowGenFunRe2);
 
 fDiffFlowGenFunIm2 = new TProfile2D("fDiffFlowGenFunIm2","Im(<D[b][2][q]>)",fgknBins,(Double_t)(fPtMin/fBinWidth),(Double_t)(fPtMax/fBinWidth),fgkQmax,0.,(Double_t)fgkQmax);
 fDiffFlowGenFunIm2->SetXTitle("b");
 fDiffFlowGenFunIm2->SetYTitle("q");
 fHistList->Add(fDiffFlowGenFunIm2);
 
 fDiffFlowGenFunRe3 = new TProfile2D("fDiffFlowGenFunRe3","Re(<D[b][3][q]>)",fgknBins,(Double_t)(fPtMin/fBinWidth),(Double_t)(fPtMax/fBinWidth),fgkQmax,0.,(Double_t)fgkQmax);
 fDiffFlowGenFunRe3->SetXTitle("b");
 fDiffFlowGenFunRe3->SetYTitle("q");
 fHistList->Add(fDiffFlowGenFunRe3);
 
 fDiffFlowGenFunIm3 = new TProfile2D("fDiffFlowGenFunIm3","Im(<D[b][3][q]>)",fgknBins,(Double_t)(fPtMin/fBinWidth),(Double_t)(fPtMax/fBinWidth),fgkQmax,0.,(Double_t)fgkQmax);
 fDiffFlowGenFunIm3->SetXTitle("b");
 fDiffFlowGenFunIm3->SetYTitle("q");
 fHistList->Add(fDiffFlowGenFunIm3);
 
 fDiffFlowGenFunRe4 = new TProfile2D("fDiffFlowGenFunRe4","Re(<D[b][4][q]>)",fgknBins,(Double_t)(fPtMin/fBinWidth),(Double_t)(fPtMax/fBinWidth),fgkQmax,0.,(Double_t)fgkQmax);
 fDiffFlowGenFunRe4->SetXTitle("b");
 fDiffFlowGenFunRe4->SetYTitle("q");
 fHistList->Add(fDiffFlowGenFunRe4);
 
 fDiffFlowGenFunIm4 = new TProfile2D("fDiffFlowGenFunIm4","Im(<D[b][4][q]>)",fgknBins,(Double_t)(fPtMin/fBinWidth),(Double_t)(fPtMax/fBinWidth),fgkQmax,0.,(Double_t)fgkQmax);
 fDiffFlowGenFunIm4->SetXTitle("b");
 fDiffFlowGenFunIm4->SetYTitle("q");
 fHistList->Add(fDiffFlowGenFunIm4);
 
 fDiffFlowGenFunRe5 = new TProfile2D("fDiffFlowGenFunRe5","Re(<D[b][5][q]>)",fgkQmax,0.,(Double_t)fgkQmax,fgknBins,(Double_t)(fPtMin/fBinWidth),(Double_t)(fPtMax/fBinWidth));
 fDiffFlowGenFunRe5->SetXTitle("b");
 fDiffFlowGenFunRe5->SetYTitle("q");
 fHistList->Add(fDiffFlowGenFunRe5);
 
 fDiffFlowGenFunIm5 = new TProfile2D("fDiffFlowGenFunIm5","Im(<D[b][5][q]>)",fgknBins,(Double_t)(fPtMin/fBinWidth),(Double_t)(fPtMax/fBinWidth),fgkQmax,0.,(Double_t)fgkQmax);
 fDiffFlowGenFunIm5->SetXTitle("b");
 fDiffFlowGenFunIm5->SetYTitle("q");
 fHistList->Add(fDiffFlowGenFunIm5);
 
 fDiffFlowGenFunRe6 = new TProfile2D("fDiffFlowGenFunRe6","Re(<D[b][6][q]>)",fgknBins,(Double_t)(fPtMin/fBinWidth),(Double_t)(fPtMax/fBinWidth),fgkQmax,0.,(Double_t)fgkQmax);
 fDiffFlowGenFunRe6->SetXTitle("b");
 fDiffFlowGenFunRe6->SetYTitle("q");
 fHistList->Add(fDiffFlowGenFunRe6);
 
 fDiffFlowGenFunIm6 = new TProfile2D("fDiffFlowGenFunIm6","Im(<D[b][6][q]>)",fgknBins,(Double_t)(fPtMin/fBinWidth),(Double_t)(fPtMax/fBinWidth),fgkQmax,0.,(Double_t)fgkQmax);
 fDiffFlowGenFunIm6->SetXTitle("b");
 fDiffFlowGenFunIm6->SetYTitle("q");
 fHistList->Add(fDiffFlowGenFunIm6);
 
 fDiffFlowGenFunRe7 = new TProfile2D("fDiffFlowGenFunRe7","Re(<D[b][7][q]>)",fgknBins,(Double_t)(fPtMin/fBinWidth),(Double_t)(fPtMax/fBinWidth),fgkQmax,0.,(Double_t)fgkQmax);
 fDiffFlowGenFunRe7->SetXTitle("b");
 fDiffFlowGenFunRe7->SetYTitle("q");
 fHistList->Add(fDiffFlowGenFunRe7);
 
 fDiffFlowGenFunIm7 = new TProfile2D("fDiffFlowGenFunIm7","Im(<D[b][7][q]>)",fgknBins,(Double_t)(fPtMin/fBinWidth),(Double_t)(fPtMax/fBinWidth),fgkQmax,0.,(Double_t)fgkQmax);
 fDiffFlowGenFunIm7->SetXTitle("b");
 fDiffFlowGenFunIm7->SetYTitle("q");
 fHistList->Add(fDiffFlowGenFunIm7);
 */
 
 //common control histograms
 fCommonHists = new AliFlowCommonHist("AliFlowCommonHistGFC");
 fHistList->Add(fCommonHists);  
 
 //common histograms for final results (2nd order)
 fCommonHistsResults2nd = new AliFlowCommonHistResults("AliFlowCommonHistResults2ndOrderGFC");
 fHistList->Add(fCommonHistsResults2nd); 
 
 //common histograms for final results (4th order)
 fCommonHistsResults4th = new AliFlowCommonHistResults("AliFlowCommonHistResults4thOrderGFC");
 fHistList->Add(fCommonHistsResults4th);
 
 //common histograms for final results (6th order)
 fCommonHistsResults6th = new AliFlowCommonHistResults("AliFlowCommonHistResults6thOrderGFC");
 fHistList->Add(fCommonHistsResults6th); 
 
 //common histograms for final results (8th order)
 fCommonHistsResults8th = new AliFlowCommonHistResults("AliFlowCommonHistResults8thOrderGFC");
 fHistList->Add(fCommonHistsResults8th);
 
 //<w^2>
 fAverageOfSquaredWeight = new TProfile("fAverageOfSquaredWeight","<w^{2}>",1,0,1);
 fAverageOfSquaredWeight->SetLabelSize(0.06);
 fAverageOfSquaredWeight->SetMarkerStyle(25);
 fAverageOfSquaredWeight->SetLabelOffset(0.01);
 (fAverageOfSquaredWeight->GetXaxis())->SetBinLabel(1,"<w^{2}>");
 fHistList->Add(fAverageOfSquaredWeight);
 
 // add list fWeightsList with weights to the main list
 fHistList->Add(fWeightsList); 
}//end of Init()

//================================================================================================================

void AliFlowAnalysisWithCumulants::Make(AliFlowEventSimple* anEvent)
{
 //running over data:
 Int_t nPrim = anEvent->NumberOfTracks(); //total multiplicity
 
 Int_t nEventNSelTracksRP = anEvent->GetEventNSelTracksRP(); //selected multiplicity (particles used for int. flow)
 
 Int_t n = 2; // int flow harmonic (to be improved)
 
 //---------------------------------------------------------------------------------------------------------
 // weights:
 Bool_t useWeights = fUsePhiWeights||fUsePtWeights||fUseEtaWeights;
 
 TH1F *phiWeights = NULL; // histogram with phi weights
 TH1D *ptWeights  = NULL; // histogram with pt weights
 TH1D *etaWeights = NULL; // histogram with eta weights

 Double_t wPhi = 1.; // phi weight
 Double_t wPt  = 1.; // pt weight
 Double_t wEta = 1.; // eta weight
   
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
 //---------------------------------------------------------------------------------------------------------
  
 if(nEventNSelTracksRP>9) //generating function formalism applied here make sense only for selected multiplicity >= 10 
 { 
 //fill the common control histograms
 fCommonHists->FillControlHistograms(anEvent);   
  
 //initializing the generating function G[p][q] for integrated flow 
 Double_t genfunG[fgkPmax][fgkQmax];
 
 for(Int_t p=0;p<fgkPmax;p++)
  {
   for(Int_t q=0;q<fgkQmax;q++)
   {
    genfunG[p][q]=1.;
   }   
  }

 Int_t nSelTracksRP = 0; //cross-checking the selected multiplicity
  
 Double_t dPhi = 0.;
 Double_t dPt  = 0.;
 Double_t dEta = 0.;
 Int_t nBinsPhi = 0;
 
 for(Int_t i=0;i<nPrim;i++)
 {
  fTrack=anEvent->GetTrack(i);
  if(fTrack && fTrack->InRPSelection())
  {
   nSelTracksRP++;
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
    if(fBinWidthPt) 
    {
     wPt = ptWeights->GetBinContent(1+(Int_t)(TMath::Floor((dPt-fPtMin)/fBinWidthPt)));
    }
   }             
   // eta weights:
   if(fUseEtaWeights)
   {    
    if(fBinWidthEta)
    {
     wEta=etaWeights->GetBinContent(1+(Int_t)(TMath::Floor((dEta-fEtaMin)/fBinWidthEta))); 
    }
   }
   // evaluate the generating function:
   for(Int_t p=0;p<fgkPmax;p++)
   {
    for(Int_t q=0;q<fgkQmax;q++)
    {
     genfunG[p][q]*=(1.+wPhi*wPt*wEta*(2.*fR0*sqrt(p+1.)/nEventNSelTracksRP)*cos(fgkFlow*dPhi-2.*q*TMath::Pi()/fgkQmax));
    }
   }
   // calculate <w^2> 
   fAverageOfSquaredWeight->Fill(0.,pow(wPhi*wPt*wEta,2.),1); 
  }
 } // end of for(Int_t i=0;i<nPrim;i++) 
 
 //storing the value of G[p][q] in 2D profile in order to get automatically the avarage <G[p][q]>
 for(Int_t p=0;p<fgkPmax;p++)
 {
  for(Int_t q=0;q<fgkQmax;q++)
  {
   fIntFlowGenFun->Fill((Double_t)p,(Double_t)q,genfunG[p][q],1);
  }
 } 
 
 //storing the selected multiplicity (if fAvMultIntFlow is not filled here then you had wrongly selected the particles used for integrated flow)
 if(nSelTracksRP==nEventNSelTracksRP)
 {
  fAvMultIntFlowGFC->Fill(0.,nSelTracksRP,1);
 }
 
 // calculating Q-vector of event (needed for errors)
 AliFlowVector fQVector;
 fQVector.Set(0.,0.);
 fQVector.SetMult(0);
 fQVector=anEvent->GetQ(1*n,fWeightsList,fUsePhiWeights,fUsePtWeights,fUseEtaWeights); // get the Q vector for this event
 fQVectorComponentsGFC->Fill(0.,fQVector.X(),1);         // in the 1st bin fill Q_x
 fQVectorComponentsGFC->Fill(1.,fQVector.Y(),1);         // in the 2nd bin fill Q_y
 fQVectorComponentsGFC->Fill(2.,pow(fQVector.X(),2.),1); // in the 3rd bin fill (Q_x)^2
 fQVectorComponentsGFC->Fill(3.,pow(fQVector.Y(),2.),1); // in the 4th bin fill (Q_y)^2
 
 //3D profiles for differential flow in pt and eta
 //second loop over data: evaluating the generating function D[b][p][q] for differential flow 
 //remark 0: D[b][p][q] is a complex number => real and imaginary part are calculated separately
 //remark 1: note that bellow G[p][q] is needed, the value of generating function for integrated flow for the CURRENT event
 //remark 2: results are stored in 3D profiles in order to automatically get <Re(D[b][p][q])> and <Im(D[b][p][q])>
 for(Int_t i=0;i<nPrim;i++)
 {
  fTrack=anEvent->GetTrack(i);
  if(fTrack)
  {
   if(fTrack->InRPSelection() && fTrack->InPOISelection())
   {
    fPtBinPOINoOfParticles->Fill(fTrack->Pt(),fTrack->Pt(),1.);
    fEtaBinPOINoOfParticles->Fill(fTrack->Eta(),fTrack->Eta(),1.);
    for(Int_t p=0;p<fgkPmax;p++)
    {
     for(Int_t q=0;q<fgkQmax;q++)
     {
      //real part (Pt)
      fDiffFlowPtPOIGenFunRe->Fill(fTrack->Pt()/fBinWidthPt,(Double_t)p,(Double_t)q,genfunG[p][q]*cos(fgkMltpl*fgkFlow*fTrack->Phi())/(1.+(2.*fR0*sqrt(p+1.)/nSelTracksRP)*cos(fgkFlow*fTrack->Phi()-2.*q*TMath::Pi()/fgkQmax)),1.);
      //imaginary part (Pt)
      fDiffFlowPtPOIGenFunIm->Fill(fTrack->Pt()/fBinWidthPt,(Double_t)p,(Double_t)q,genfunG[p][q]*sin(fgkMltpl*fgkFlow*fTrack->Phi())/(1.+(2.*fR0*sqrt(p+1.)/nSelTracksRP)*cos(fgkFlow*fTrack->Phi()-2.*q*TMath::Pi()/fgkQmax)),1.); 
      //real part (Eta)
      fDiffFlowEtaPOIGenFunRe->Fill(fTrack->Eta()/fBinWidthEta,(Double_t)p,(Double_t)q,genfunG[p][q]*cos(fgkMltpl*fgkFlow*fTrack->Phi())/(1.+(2.*fR0*sqrt(p+1.)/nSelTracksRP)*cos(fgkFlow*fTrack->Phi()-2.*q*TMath::Pi()/fgkQmax)),1.);
      //imaginary part (Eta)
      fDiffFlowEtaPOIGenFunIm->Fill(fTrack->Eta()/fBinWidthEta,(Double_t)p,(Double_t)q,genfunG[p][q]*sin(fgkMltpl*fgkFlow*fTrack->Phi())/(1.+(2.*fR0*sqrt(p+1.)/nSelTracksRP)*cos(fgkFlow*fTrack->Phi()-2.*q*TMath::Pi()/fgkQmax)),1.); 
     }
    }
   } 
   else if(fTrack->InPOISelection() && !(fTrack->InRPSelection()))
   {
    fPtBinPOINoOfParticles->Fill(fTrack->Pt(),fTrack->Pt(),1.);
    fEtaBinPOINoOfParticles->Fill(fTrack->Eta(),fTrack->Eta(),1.);
    for(Int_t p=0;p<fgkPmax;p++)
    {
     for(Int_t q=0;q<fgkQmax;q++)
     {
      //real part (Pt)
      fDiffFlowPtPOIGenFunRe->Fill(fTrack->Pt()/fBinWidthPt,(Double_t)p,(Double_t)q,genfunG[p][q]*cos(fgkMltpl*fgkFlow*fTrack->Phi()),1.);
      //imaginary part (Pt)
      fDiffFlowPtPOIGenFunIm->Fill(fTrack->Pt()/fBinWidthPt,(Double_t)p,(Double_t)q,genfunG[p][q]*sin(fgkMltpl*fgkFlow*fTrack->Phi()),1.); 
      //real part (Eta)
      fDiffFlowEtaPOIGenFunRe->Fill(fTrack->Eta()/fBinWidthEta,(Double_t)p,(Double_t)q,genfunG[p][q]*cos(fgkMltpl*fgkFlow*fTrack->Phi()),1.);
      //imaginary part (Eta)
      fDiffFlowEtaPOIGenFunIm->Fill(fTrack->Eta()/fBinWidthEta,(Double_t)p,(Double_t)q,genfunG[p][q]*sin(fgkMltpl*fgkFlow*fTrack->Phi()),1.);                     
     } 
    }
   }
   //RP:
   if(fTrack->InRPSelection())       
   {
    fPtBinRPNoOfParticles->Fill(fTrack->Pt(),fTrack->Pt(),1.);
    fEtaBinRPNoOfParticles->Fill(fTrack->Eta(),fTrack->Eta(),1.);
    for(Int_t p=0;p<fgkPmax;p++)
    {
     for(Int_t q=0;q<fgkQmax;q++)
     {
      //real part (Pt)
      fDiffFlowPtRPGenFunRe->Fill(fTrack->Pt()/fBinWidthPt,(Double_t)p,(Double_t)q,genfunG[p][q]*cos(fgkMltpl*fgkFlow*fTrack->Phi())/(1.+(2.*fR0*sqrt(p+1.)/nSelTracksRP)*cos(fgkFlow*fTrack->Phi()-2.*q*TMath::Pi()/fgkQmax)),1.);
      //imaginary part (Pt)
      fDiffFlowPtRPGenFunIm->Fill(fTrack->Pt()/fBinWidthPt,(Double_t)p,(Double_t)q,genfunG[p][q]*sin(fgkMltpl*fgkFlow*fTrack->Phi())/(1.+(2.*fR0*sqrt(p+1.)/nSelTracksRP)*cos(fgkFlow*fTrack->Phi()-2.*q*TMath::Pi()/fgkQmax)),1.); 
      //real part (Eta)
      fDiffFlowEtaRPGenFunRe->Fill(fTrack->Eta()/fBinWidthEta,(Double_t)p,(Double_t)q,genfunG[p][q]*cos(fgkMltpl*fgkFlow*fTrack->Phi())/(1.+(2.*fR0*sqrt(p+1.)/nSelTracksRP)*cos(fgkFlow*fTrack->Phi()-2.*q*TMath::Pi()/fgkQmax)),1.);
      //imaginary part (Eta)
      fDiffFlowEtaRPGenFunIm->Fill(fTrack->Eta()/fBinWidthEta,(Double_t)p,(Double_t)q,genfunG[p][q]*sin(fgkMltpl*fgkFlow*fTrack->Phi())/(1.+(2.*fR0*sqrt(p+1.)/nSelTracksRP)*cos(fgkFlow*fTrack->Phi()-2.*q*TMath::Pi()/fgkQmax)),1.); 
     }
    }
   }//end of if(fTrack->InRPSelection())                   
  }//end of if(fTrack)  
 }//ending the second loop over data    
 
 
 
 /*
 //sixteen 2D profiles for differential flow          
 for(Int_t i=0;i<nPrim;i++){
  fTrack=anEvent->GetTrack(i);
  if (fTrack && fTrack->InPOISelection()){
   //for(Int_t p=0;p<fgkPmax;p++){
    for(Int_t q=0;q<fgkQmax;q++){
     //real part
     fDiffFlowGenFunRe0->Fill(fTrack->Pt()/fBinWidth,(Double_t)q,genfunG[0][q]*cos(fgkMltpl*fgkFlow*fTrack->Phi())/(1.+(2.*fR0*sqrt(0.+1.)/nSelTracksRP)*cos(fgkFlow*fTrack->Phi()-2.*q*TMath::Pi()/fgkQmax)),1.);
     //imaginary part
     fDiffFlowGenFunIm0->Fill(fTrack->Pt()/fBinWidth,(Double_t)q,genfunG[0][q]*sin(fgkMltpl*fgkFlow*fTrack->Phi())/(1.+(2.*fR0*sqrt(0.+1.)/nSelTracksRP)*cos(fgkFlow*fTrack->Phi()-2.*q*TMath::Pi()/fgkQmax)),1.); 
     //-----------------------------------------------------------------------
     //real part
     fDiffFlowGenFunRe1->Fill(fTrack->Pt()/fBinWidth,(Double_t)q,genfunG[1][q]*cos(fgkMltpl*fgkFlow*fTrack->Phi())/(1.+(2.*fR0*sqrt(1.+1.)/nSelTracksRP)*cos(fgkFlow*fTrack->Phi()-2.*q*TMath::Pi()/fgkQmax)),1.);
     //imaginary part
     fDiffFlowGenFunIm1->Fill(fTrack->Pt()/fBinWidth,(Double_t)q,genfunG[1][q]*sin(fgkMltpl*fgkFlow*fTrack->Phi())/(1.+(2.*fR0*sqrt(1.+1.)/nSelTracksRP)*cos(fgkFlow*fTrack->Phi()-2.*q*TMath::Pi()/fgkQmax)),1.);             
    //-----------------------------------------------------------------------
     //real part
     fDiffFlowGenFunRe2->Fill(fTrack->Pt()/fBinWidth,(Double_t)q,genfunG[2][q]*cos(fgkMltpl*fgkFlow*fTrack->Phi())/(1.+(2.*fR0*sqrt(2.+1.)/nSelTracksRP)*cos(fgkFlow*fTrack->Phi()-2.*q*TMath::Pi()/fgkQmax)),1.);
     //imaginary part
     fDiffFlowGenFunIm2->Fill(fTrack->Pt()/fBinWidth,(Double_t)q,genfunG[2][q]*sin(fgkMltpl*fgkFlow*fTrack->Phi())/(1.+(2.*fR0*sqrt(2.+1.)/nSelTracksRP)*cos(fgkFlow*fTrack->Phi()-2.*q*TMath::Pi()/fgkQmax)),1.);  
     //-----------------------------------------------------------------------
     //real part
     fDiffFlowGenFunRe3->Fill(fTrack->Pt()/fBinWidth,(Double_t)q,genfunG[3][q]*cos(fgkMltpl*fgkFlow*fTrack->Phi())/(1.+(2.*fR0*sqrt(3.+1.)/nSelTracksRP)*cos(fgkFlow*fTrack->Phi()-2.*q*TMath::Pi()/fgkQmax)),1.);
     //imaginary part
     fDiffFlowGenFunIm3->Fill(fTrack->Pt()/fBinWidth,(Double_t)q,genfunG[3][q]*sin(fgkMltpl*fgkFlow*fTrack->Phi())/(1.+(2.*fR0*sqrt(3.+1.)/nSelTracksRP)*cos(fgkFlow*fTrack->Phi()-2.*q*TMath::Pi()/fgkQmax)),1.);  
     //-----------------------------------------------------------------------
     //real part
     fDiffFlowGenFunRe4->Fill(fTrack->Pt()/fBinWidth,(Double_t)q,genfunG[4][q]*cos(fgkMltpl*fgkFlow*fTrack->Phi())/(1.+(2.*fR0*sqrt(4.+1.)/nSelTracksRP)*cos(fgkFlow*fTrack->Phi()-2.*q*TMath::Pi()/fgkQmax)),1.);
     //imaginary part
     fDiffFlowGenFunIm4->Fill(fTrack->Pt()/fBinWidth,(Double_t)q,genfunG[4][q]*sin(fgkMltpl*fgkFlow*fTrack->Phi())/(1.+(2.*fR0*sqrt(4.+1.)/nSelTracksRP)*cos(fgkFlow*fTrack->Phi()-2.*q*TMath::Pi()/fgkQmax)),1.);  
     //-----------------------------------------------------------------------
     //real part
     fDiffFlowGenFunRe5->Fill(fTrack->Pt()/fBinWidth,(Double_t)q,genfunG[5][q]*cos(fgkMltpl*fgkFlow*fTrack->Phi())/(1.+(2.*fR0*sqrt(5.+1.)/nSelTracksRP)*cos(fgkFlow*fTrack->Phi()-2.*q*TMath::Pi()/fgkQmax)),1.);
     //imaginary part
     fDiffFlowGenFunIm5->Fill(fTrack->Pt()/fBinWidth,(Double_t)q,genfunG[5][q]*sin(fgkMltpl*fgkFlow*fTrack->Phi())/(1.+(2.*fR0*sqrt(5.+1.)/nSelTracksRP)*cos(fgkFlow*fTrack->Phi()-2.*q*TMath::Pi()/fgkQmax)),1.);  
     //-----------------------------------------------------------------------
     //real part
     fDiffFlowGenFunRe6->Fill(fTrack->Pt()/fBinWidth,(Double_t)q,genfunG[6][q]*cos(fgkMltpl*fgkFlow*fTrack->Phi())/(1.+(2.*fR0*sqrt(6.+1.)/nSelTracksRP)*cos(fgkFlow*fTrack->Phi()-2.*q*TMath::Pi()/fgkQmax)),1.);
     //imaginary part
     fDiffFlowGenFunIm6->Fill(fTrack->Pt()/fBinWidth,(Double_t)q,genfunG[6][q]*sin(fgkMltpl*fgkFlow*fTrack->Phi())/(1.+(2.*fR0*sqrt(6.+1.)/nSelTracksRP)*cos(fgkFlow*fTrack->Phi()-2.*q*TMath::Pi()/fgkQmax)),1.);
     //-----------------------------------------------------------------------
     //real part
     fDiffFlowGenFunRe7->Fill(fTrack->Pt()/fBinWidth,(Double_t)q,genfunG[7][q]*cos(fgkMltpl*fgkFlow*fTrack->Phi())/(1.+(2.*fR0*sqrt(7.+1.)/nSelTracksRP)*cos(fgkFlow*fTrack->Phi()-2.*q*TMath::Pi()/fgkQmax)),1.);
     //imaginary part
     fDiffFlowGenFunIm7->Fill(fTrack->Pt()/fBinWidth,(Double_t)q,genfunG[7][q]*sin(fgkMltpl*fgkFlow*fTrack->Phi())/(1.+(2.*fR0*sqrt(7.+1.)/nSelTracksRP)*cos(fgkFlow*fTrack->Phi()-2.*q*TMath::Pi()/fgkQmax)),1.);    
   }
   //}
  }  
 }//ending the second loop over data            
 */
      
 }//end of if(nEventNSelTracksRP>9)                   
 

 
  
   
    
     
      
       
        
         
          
           
 //off the record: numerical equations for cumulants solved up to different highest order  
 if(fOtherEquations)
 {
  //running over data
  Int_t nPrimOE = anEvent->NumberOfTracks();//total multiplicity 
  
  Int_t nEventNSelTracksRPOE = anEvent->GetEventNSelTracksRP();
  
  Double_t genfunG4[fgkPmax4][fgkQmax4];
  Double_t genfunG6[fgkPmax6][fgkQmax6];
  Double_t genfunG8[fgkPmax8][fgkQmax8];
  Double_t genfunG16[fgkPmax16][fgkQmax16];
  for(Int_t p=0;p<fgkPmax16;p++)
  {
   for(Int_t q=0;q<fgkQmax16;q++)
   {
    genfunG16[p][q]=1.;
    if(p<fgkPmax8 && q<fgkQmax8)
    {
     genfunG8[p][q]=1.;
     if(p<fgkPmax6 && q<fgkQmax6)
     {
      genfunG6[p][q]=1.;
      if(p<fgkPmax4 && q<fgkQmax4)
      {
       genfunG4[p][q]=1.;
      }
     }
    } 
   }
  } 
   
  //multiplicities: 
  if(nEventNSelTracksRPOE>15) fAvMultIntFlow16GFC->Fill(0.,nEventNSelTracksRPOE,1);
  if(nEventNSelTracksRPOE>7) fAvMultIntFlow8GFC->Fill(0.,nEventNSelTracksRPOE,1);
  if(nEventNSelTracksRPOE>5) fAvMultIntFlow6GFC->Fill(0.,nEventNSelTracksRPOE,1);
  if(nEventNSelTracksRPOE>3) fAvMultIntFlow4GFC->Fill(0.,nEventNSelTracksRPOE,1);  
  
  //first loop over data: evaluating the generating function G[p][q] for integrated flow 
  for(Int_t i=0;i<nPrimOE;i++)
  {
   fTrack=anEvent->GetTrack(i);
   if(fTrack && fTrack->InRPSelection())
   {
    for(Int_t p=0;p<fgkPmax16;p++)
    {
     for(Int_t q=0;q<fgkQmax16;q++)
     {
      if(nEventNSelTracksRPOE>15)
      {
       genfunG16[p][q]*=(1.+(2.*fR0*sqrt(p+1.)/nEventNSelTracksRPOE)*cos(fgkFlow*fTrack->Phi()-2.*q*TMath::Pi()/fgkQmax16));
      }       
      if(p<fgkPmax8 && q<fgkQmax8)
      {
       if(nEventNSelTracksRPOE>7)
       { 
        genfunG8[p][q]*=(1.+(2.*fR0*sqrt(p+1.)/nEventNSelTracksRPOE)*cos(fgkFlow*fTrack->Phi()-2.*q*TMath::Pi()/fgkQmax8));
       }
       if(p<fgkPmax6 && q<fgkQmax6)
       {
        if(nEventNSelTracksRPOE>5) 
        {
         genfunG6[p][q]*=(1.+(2.*fR0*sqrt(p+1.)/nEventNSelTracksRPOE)*cos(fgkFlow*fTrack->Phi()-2.*q*TMath::Pi()/fgkQmax6));
        }
        if(p<fgkPmax4 && q<fgkQmax4)
        {
         if(nEventNSelTracksRPOE>3)
         {
          genfunG4[p][q]*=(1.+(2.*fR0*sqrt(p+1.)/nEventNSelTracksRPOE)*cos(fgkFlow*fTrack->Phi()-2.*q*TMath::Pi()/fgkQmax4));
         } 
        }
       }
      }
     } 
    }  
   }//end of if(fTrack && fTrack->InRPSelection())
  }//ending the loop over data
 
 //storing the value of G[p][q] in 2D profile in order to get automatically the avarage <G[p][q]>
 for(Int_t p=0;p<fgkPmax16;p++)
 {
  for(Int_t q=0;q<fgkQmax16;q++)
  {
   if(nEventNSelTracksRPOE>15) fIntFlowGenFun16->Fill((Double_t)p,(Double_t)q,genfunG16[p][q],1);
   if(p<fgkPmax8 && q<fgkQmax8)
   {
    if(nEventNSelTracksRPOE>7) fIntFlowGenFun8->Fill((Double_t)p,(Double_t)q,genfunG8[p][q],1);
    if(p<fgkPmax6 && q<fgkQmax6)
    {
     if(nEventNSelTracksRPOE>5) fIntFlowGenFun6->Fill((Double_t)p,(Double_t)q,genfunG6[p][q],1);
     if(p<fgkPmax4 && q<fgkQmax4)
     {
      if(nEventNSelTracksRPOE>3) fIntFlowGenFun4->Fill((Double_t)p,(Double_t)q,genfunG4[p][q],1);
     }
    }
   } 
  }
 }
}//end of if(fOtherEquations)                  
                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
}//end of Make()

//================================================================================================================

void AliFlowAnalysisWithCumulants::Finish()
{
 //calculate the final results
 //AliCumulantsFunctions finalResults(fIntFlowGenFun,NULL,NULL, fIntFlowResults,fDiffFlowResults2,fDiffFlowResults4,fDiffFlowResults6,fDiffFlowResults8,fAvMultIntFlow,fQVectorComponents,  fQDist,fDiffFlowGenFunRe0,fDiffFlowGenFunRe1,fDiffFlowGenFunRe2, fDiffFlowGenFunRe3,fDiffFlowGenFunRe4,fDiffFlowGenFunRe5,fDiffFlowGenFunRe6,fDiffFlowGenFunRe7,fDiffFlowGenFunIm0,fDiffFlowGenFunIm1, fDiffFlowGenFunIm2,fDiffFlowGenFunIm3,fDiffFlowGenFunIm4,fDiffFlowGenFunIm5,fDiffFlowGenFunIm6,fDiffFlowGenFunIm7);

 AliCumulantsFunctions finalResults(fIntFlowGenFun,fIntFlowGenFun4,fIntFlowGenFun6,fIntFlowGenFun8,fIntFlowGenFun16,fAvMultIntFlow4GFC, fAvMultIntFlow6GFC,fAvMultIntFlow8GFC,fAvMultIntFlow16GFC,fDiffFlowPtRPGenFunRe,fDiffFlowPtRPGenFunIm,fPtBinRPNoOfParticles, fDiffFlowEtaRPGenFunRe,fDiffFlowEtaRPGenFunIm,fEtaBinRPNoOfParticles,fDiffFlowPtPOIGenFunRe,fDiffFlowPtPOIGenFunIm,fPtBinPOINoOfParticles, fDiffFlowEtaPOIGenFunRe,fDiffFlowEtaPOIGenFunIm,fEtaBinPOINoOfParticles,fIntFlowResultsGFC,fDiffFlowResults2ndOrderGFC,fDiffFlowResults4thOrderGFC,fDiffFlowResults6thOrderGFC,fDiffFlowResults8thOrderGFC, fAvMultIntFlowGFC,fQVectorComponentsGFC,fAverageOfSquaredWeight,fCommonHistsResults2nd, fCommonHistsResults4th,fCommonHistsResults6th,fCommonHistsResults8th,fCommonHists);
                           
 finalResults.Calculate();  
}

//================================================================================================================

void AliFlowAnalysisWithCumulants::WriteHistograms(TString* outputFileName)
{
 //store the final results in output .root file
 TFile *output = new TFile(outputFileName->Data(),"RECREATE");
 output->WriteObject(fHistList, "cobjGFC","SingleKey");
 delete output;
}

//================================================================================================================

//================================================================================================================

void AliFlowAnalysisWithCumulants::WriteHistograms(TString outputFileName)
{
 //store the final results in output .root file
 TFile *output = new TFile(outputFileName.Data(),"RECREATE");
 output->WriteObject(fHistList, "cobjGFC","SingleKey");
 delete output;
}

//================================================================================================================




