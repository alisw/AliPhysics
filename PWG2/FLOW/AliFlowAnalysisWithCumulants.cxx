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
 fR0(0),
 fPtMax(0),
 fPtMin(0),
 fBinWidth(0),
 fgknBins(0),
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
 fIntFlowGenFun4(NULL),
 fIntFlowGenFun6(NULL),
 fIntFlowGenFun8(NULL),
 fIntFlowGenFun16(NULL),
 fDiffFlowGenFunRe(NULL),
 fDiffFlowGenFunIm(NULL),
 fBinNoOfParticles(NULL),
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
 fOtherEquations(kFALSE)
{
 //constructor 
 fHistList = new TList();
 fR0=AliFlowCumuConstants::fgR0;
 fPtMax=AliFlowCommonConstants::GetPtMax(); 
 fPtMin=AliFlowCommonConstants::GetPtMin();
 fgknBins=AliFlowCommonConstants::GetNbinsPt();
 fOtherEquations=AliFlowCumuConstants::fgOtherEquations;
 if(fgknBins)
 {
  fBinWidth=(fPtMax-fPtMin)/fgknBins;  
 } 
}

AliFlowAnalysisWithCumulants::~AliFlowAnalysisWithCumulants()
{
 //desctructor
 delete fHistList; 
}

//================================================================================================================

void AliFlowAnalysisWithCumulants::CreateOutputObjects()
{
 //various output histograms
 
 //average multiplicity 
 fAvMultIntFlowGFC = new TProfile("fAvMultIntFlowGFC","Average Multiplicity",1,0,1,"s");
 fAvMultIntFlowGFC->SetXTitle("");
 fAvMultIntFlowGFC->SetYTitle("");
 fAvMultIntFlowGFC->SetLabelSize(0.06);
 fAvMultIntFlowGFC->SetMarkerStyle(25);
 fAvMultIntFlowGFC->SetLabelOffset(0.01);
 (fAvMultIntFlowGFC->GetXaxis())->SetBinLabel(1,"Average Multiplicity");
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
 fDiffFlowResults2ndOrderGFC = new TH1D("fDiffFlowResults2ndOrderGFC","v'_2/2{2}",fgknBins,fPtMin,fPtMax);
 fDiffFlowResults2ndOrderGFC->SetXTitle("pt [GeV]");
 fDiffFlowResults2ndOrderGFC->SetYTitle("");
 fHistList->Add(fDiffFlowResults2ndOrderGFC);

 //final results for differential flow v_p/n{4} (by default p=n=2) 
 fDiffFlowResults4thOrderGFC = new TH1D("fDiffFlowResults4thOrderGFC","v'_2/2{4}",fgknBins,fPtMin,fPtMax);
 fDiffFlowResults4thOrderGFC->SetXTitle("pt [GeV]");
 fDiffFlowResults4thOrderGFC->SetYTitle("");
 fHistList->Add(fDiffFlowResults4thOrderGFC);
 
 //final results for differential flow v_p/n{6} (by default p=n=2)  
 fDiffFlowResults6thOrderGFC = new TH1D("fDiffFlowResults6thOrderGFC","v'_2/2{6}",fgknBins,fPtMin,fPtMax);
 fDiffFlowResults6thOrderGFC->SetXTitle("pt [GeV]");
 fDiffFlowResults6thOrderGFC->SetYTitle("");
 fHistList->Add(fDiffFlowResults6thOrderGFC);
 
 //final results for differential flow v_p/n{8} (by default p=n=2)
 fDiffFlowResults8thOrderGFC = new TH1D("fDiffFlowResults8thOrderGFC","v'_2/2{8}",fgknBins,fPtMin,fPtMax);
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
  //avarage of the generating function for integrated flow <G[p][q]> (up to 4th order)
  fIntFlowGenFun4 = new TProfile2D("fIntFlowGenFun4","<G4[p4][q4]>",fgkPmax4,0.,(Double_t)fgkPmax4,fgkQmax4,0.,(Double_t)fgkQmax4);
  fIntFlowGenFun4->SetXTitle("p4");
  fIntFlowGenFun4->SetYTitle("q4");
  fHistList->Add(fIntFlowGenFun4);

  //avarage of the generating function for integrated flow <G[p][q]> (up to 6th order) 
  fIntFlowGenFun6 = new TProfile2D("fIntFlowGenFun6","<G6[p6][q6]>",fgkPmax6,0.,(Double_t)fgkPmax6,fgkQmax6,0.,(Double_t)fgkQmax6);
  fIntFlowGenFun6->SetXTitle("p6");
  fIntFlowGenFun6->SetYTitle("q6");
  fHistList->Add(fIntFlowGenFun6);

  //avarage of the generating function for integrated flow <G[p][q]> (up to 8th order)
  fIntFlowGenFun8 = new TProfile2D("fIntFlowGenFun8","<G8[p8][q8]>",fgkPmax8,0.,(Double_t)fgkPmax8,fgkQmax8,0.,(Double_t)fgkQmax8);
  fIntFlowGenFun8->SetXTitle("p8");
  fIntFlowGenFun8->SetYTitle("q8");
  fHistList->Add(fIntFlowGenFun8);

  //avarage of the generating function for integrated flow <G[p][q]> (up to 16th order)
  fIntFlowGenFun16 = new TProfile2D("fIntFlowGenFun16","<G16[p16][q16]>",fgkPmax16,0.,(Double_t)fgkPmax16,fgkQmax16,0.,(Double_t)fgkQmax16);
  fIntFlowGenFun16->SetXTitle("p16");
  fIntFlowGenFun16->SetYTitle("q16");
  fHistList->Add(fIntFlowGenFun16);
 }
 
 //avarage of the real part of generating function for differential flow <Re(D[b][p][q])>
 fDiffFlowGenFunRe = new TProfile3D("fDiffFlowGenFunRe","<Re(D[b][p][q])>",fgknBins,(Double_t)(fPtMin/fBinWidth),(Double_t)(fPtMax/fBinWidth),fgkPmax,0.,(Double_t)fgkPmax,fgkQmax,0.,(Double_t)fgkQmax);
 fDiffFlowGenFunRe->SetXTitle("b");
 fDiffFlowGenFunRe->SetYTitle("p");
 fDiffFlowGenFunRe->SetZTitle("q");
 fDiffFlowGenFunRe->SetTitleOffset(1.44,"X");
 fDiffFlowGenFunRe->SetTitleOffset(1.44,"Y");
 fHistList->Add(fDiffFlowGenFunRe);
 
 //avarage of the imaginary part of generating function for differential flow <Im(D[b][p][q])>
 fDiffFlowGenFunIm = new TProfile3D("fDiffFlowGenFunIm","<Im(D[b][p][q])>",fgknBins,(Double_t)(fPtMin/fBinWidth),(Double_t)(fPtMax/fBinWidth),fgkPmax,0.,(Double_t)fgkPmax,fgkQmax,0.,(Double_t)fgkQmax);
 fDiffFlowGenFunIm->SetXTitle("b");
 fDiffFlowGenFunIm->SetYTitle("p");
 fDiffFlowGenFunIm->SetZTitle("q");
 fDiffFlowGenFunIm->SetTitleOffset(1.44,"X");
 fDiffFlowGenFunIm->SetTitleOffset(1.44,"Y");
 fHistList->Add(fDiffFlowGenFunIm);
 
 //number of particles per pt bin
 fBinNoOfParticles = new TProfile("fBinNoOfParticles","Number of particles per pt bin",fgknBins,fPtMin,fPtMax);
 fBinNoOfParticles->SetXTitle("pt [GeV]");
 fBinNoOfParticles->SetYTitle("");
 fHistList->Add(fBinNoOfParticles);
 
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
  
}//end of CreateOutputObjects()

//================================================================================================================

void AliFlowAnalysisWithCumulants::Make(AliFlowEventSimple* anEvent)
{
 //running over data
 Int_t nPrim = anEvent->NumberOfTracks();//total multiplicity
  
 //if(nPrim>30)//generating function formalism can be applied only for large multiplicities (to be improved in the future) 
 //{ 
 //fill the common control histograms
 fCommonHists->FillControlHistograms(anEvent);   
  
 //initializing the generating function G[p][q] for integrated flow 
 Double_t G[fgkPmax][fgkQmax];
 
 for(Int_t p=0;p<fgkPmax;p++)
  {
   for(Int_t q=0;q<fgkQmax;q++)
   {
    G[p][q]=1.;
   }   
  }
 
 Int_t nEventNSelTracksIntFlow = anEvent->GetEventNSelTracksIntFlow(); //selected multiplicity (particles used for int. flow)

 Int_t nSelTracksIntFlow = 0; //cross-checking the selected multiplicity
 
 for(Int_t i=0;i<nPrim;i++)
 {
  fTrack=anEvent->GetTrack(i);
  if(fTrack && fTrack->UseForIntegratedFlow())
  {
   nSelTracksIntFlow++;
   for(Int_t p=0;p<fgkPmax;p++)
   {
    for(Int_t q=0;q<fgkQmax;q++)
    {
     G[p][q]*=(1.+(2.*fR0*sqrt(p+1.)/nEventNSelTracksIntFlow)*cos(fgkFlow*fTrack->Phi()-2.*q*TMath::Pi()/fgkQmax));
    }
   } 
  }
 }  
 
 //storing the value of G[p][q] in 2D profile in order to get automatically the avarage <G[p][q]>
 for(Int_t p=0;p<fgkPmax;p++)
 {
  for(Int_t q=0;q<fgkQmax;q++)
  {
   fIntFlowGenFun->Fill((Double_t)p,(Double_t)q,G[p][q],1);
  }
 } 
 
 //storing the selected multiplicity (if fAvMultIntFlow is not filled here then you had wrongly selected the particles used for integrated flow)
 if(nSelTracksIntFlow==nEventNSelTracksIntFlow)
 {
  fAvMultIntFlowGFC->Fill(0.,nSelTracksIntFlow,1);
 }
 
 //calculating Q-vector of event (needed for errors)
 AliFlowVector fQVector;
 fQVector.Set(0.,0.);
 fQVector.SetMult(0);
 fQVector=anEvent->GetQ();                               //get the Q vector for this event
 fQVectorComponentsGFC->Fill(0.,fQVector.X(),1);         //in the 1st bin fill Q_x
 fQVectorComponentsGFC->Fill(1.,fQVector.Y(),1);         //in the 2nd bin fill Q_y
 fQVectorComponentsGFC->Fill(2.,pow(fQVector.X(),2.),1); //in the 3rd bin fill (Q_x)^2
 fQVectorComponentsGFC->Fill(3.,pow(fQVector.Y(),2.),1); //in the 4th bin fill (Q_y)^2
 
 //two 3D profiles for differential flow
 //second loop over data: evaluating the generating function D[b][p][q] for differential flow 
 //remark 0: D[b][p][q] is a complex number => real and imaginary part are calculated separately
 //remark 1: note that bellow G[p][q] is needed, the value of generating function for integrated flow for the CURRENT event
 //remark 2: results are stored in two 3D profiles in order to automatically get <Re(D[b][p][q])> and <Im(D[b][p][q])>
 for(Int_t i=0;i<nPrim;i++)
 {
  fTrack=anEvent->GetTrack(i);
  if(fTrack)
  {
   fBinNoOfParticles->Fill(fTrack->Pt(),fTrack->Pt(),1.);
   if(fTrack->UseForIntegratedFlow() && fTrack->UseForDifferentialFlow())
   {
    for(Int_t p=0;p<fgkPmax;p++)
    {
     for(Int_t q=0;q<fgkQmax;q++)
     {
      //real part
      fDiffFlowGenFunRe->Fill(fTrack->Pt()/fBinWidth,(Double_t)p,(Double_t)q,G[p][q]*cos(fgkMltpl*fgkFlow*fTrack->Phi())/(1.+(2.*fR0*sqrt(p+1.)/nSelTracksIntFlow)*cos(fgkFlow*fTrack->Phi()-2.*q*TMath::Pi()/fgkQmax)),1.);
      //imaginary part
      fDiffFlowGenFunIm->Fill(fTrack->Pt()/fBinWidth,(Double_t)p,(Double_t)q,G[p][q]*sin(fgkMltpl*fgkFlow*fTrack->Phi())/(1.+(2.*fR0*sqrt(p+1.)/nSelTracksIntFlow)*cos(fgkFlow*fTrack->Phi()-2.*q*TMath::Pi()/fgkQmax)),1.);        
     }
    }
   } 
   else if(fTrack->UseForDifferentialFlow() && !(fTrack->UseForIntegratedFlow()))
   {
    for(Int_t p=0;p<fgkPmax;p++)
    {
     for(Int_t q=0;q<fgkQmax;q++)
     {
      //real part
      fDiffFlowGenFunRe->Fill(fTrack->Pt()/fBinWidth,(Double_t)p,(Double_t)q,G[p][q]*cos(fgkMltpl*fgkFlow*fTrack->Phi()),1.);
      //imaginary part
      fDiffFlowGenFunIm->Fill(fTrack->Pt()/fBinWidth,(Double_t)p,(Double_t)q,G[p][q]*sin(fgkMltpl*fgkFlow*fTrack->Phi()),1.);        
     } 
    }
   }      
  }//end of if(fTrack)  
 }//ending the second loop over data    
 
 
 
 /*
 //sixteen 2D profiles for differential flow          
 for(Int_t i=0;i<nPrim;i++){
  fTrack=anEvent->GetTrack(i);
  if (fTrack && fTrack->UseForDifferentialFlow()){
   //for(Int_t p=0;p<fgkPmax;p++){
    for(Int_t q=0;q<fgkQmax;q++){
     //real part
     fDiffFlowGenFunRe0->Fill(fTrack->Pt()/fBinWidth,(Double_t)q,G[0][q]*cos(fgkMltpl*fgkFlow*fTrack->Phi())/(1.+(2.*fR0*sqrt(0.+1.)/nSelTracksIntFlow)*cos(fgkFlow*fTrack->Phi()-2.*q*TMath::Pi()/fgkQmax)),1.);
     //imaginary part
     fDiffFlowGenFunIm0->Fill(fTrack->Pt()/fBinWidth,(Double_t)q,G[0][q]*sin(fgkMltpl*fgkFlow*fTrack->Phi())/(1.+(2.*fR0*sqrt(0.+1.)/nSelTracksIntFlow)*cos(fgkFlow*fTrack->Phi()-2.*q*TMath::Pi()/fgkQmax)),1.); 
     //-----------------------------------------------------------------------
     //real part
     fDiffFlowGenFunRe1->Fill(fTrack->Pt()/fBinWidth,(Double_t)q,G[1][q]*cos(fgkMltpl*fgkFlow*fTrack->Phi())/(1.+(2.*fR0*sqrt(1.+1.)/nSelTracksIntFlow)*cos(fgkFlow*fTrack->Phi()-2.*q*TMath::Pi()/fgkQmax)),1.);
     //imaginary part
     fDiffFlowGenFunIm1->Fill(fTrack->Pt()/fBinWidth,(Double_t)q,G[1][q]*sin(fgkMltpl*fgkFlow*fTrack->Phi())/(1.+(2.*fR0*sqrt(1.+1.)/nSelTracksIntFlow)*cos(fgkFlow*fTrack->Phi()-2.*q*TMath::Pi()/fgkQmax)),1.);             
    //-----------------------------------------------------------------------
     //real part
     fDiffFlowGenFunRe2->Fill(fTrack->Pt()/fBinWidth,(Double_t)q,G[2][q]*cos(fgkMltpl*fgkFlow*fTrack->Phi())/(1.+(2.*fR0*sqrt(2.+1.)/nSelTracksIntFlow)*cos(fgkFlow*fTrack->Phi()-2.*q*TMath::Pi()/fgkQmax)),1.);
     //imaginary part
     fDiffFlowGenFunIm2->Fill(fTrack->Pt()/fBinWidth,(Double_t)q,G[2][q]*sin(fgkMltpl*fgkFlow*fTrack->Phi())/(1.+(2.*fR0*sqrt(2.+1.)/nSelTracksIntFlow)*cos(fgkFlow*fTrack->Phi()-2.*q*TMath::Pi()/fgkQmax)),1.);  
     //-----------------------------------------------------------------------
     //real part
     fDiffFlowGenFunRe3->Fill(fTrack->Pt()/fBinWidth,(Double_t)q,G[3][q]*cos(fgkMltpl*fgkFlow*fTrack->Phi())/(1.+(2.*fR0*sqrt(3.+1.)/nSelTracksIntFlow)*cos(fgkFlow*fTrack->Phi()-2.*q*TMath::Pi()/fgkQmax)),1.);
     //imaginary part
     fDiffFlowGenFunIm3->Fill(fTrack->Pt()/fBinWidth,(Double_t)q,G[3][q]*sin(fgkMltpl*fgkFlow*fTrack->Phi())/(1.+(2.*fR0*sqrt(3.+1.)/nSelTracksIntFlow)*cos(fgkFlow*fTrack->Phi()-2.*q*TMath::Pi()/fgkQmax)),1.);  
     //-----------------------------------------------------------------------
     //real part
     fDiffFlowGenFunRe4->Fill(fTrack->Pt()/fBinWidth,(Double_t)q,G[4][q]*cos(fgkMltpl*fgkFlow*fTrack->Phi())/(1.+(2.*fR0*sqrt(4.+1.)/nSelTracksIntFlow)*cos(fgkFlow*fTrack->Phi()-2.*q*TMath::Pi()/fgkQmax)),1.);
     //imaginary part
     fDiffFlowGenFunIm4->Fill(fTrack->Pt()/fBinWidth,(Double_t)q,G[4][q]*sin(fgkMltpl*fgkFlow*fTrack->Phi())/(1.+(2.*fR0*sqrt(4.+1.)/nSelTracksIntFlow)*cos(fgkFlow*fTrack->Phi()-2.*q*TMath::Pi()/fgkQmax)),1.);  
     //-----------------------------------------------------------------------
     //real part
     fDiffFlowGenFunRe5->Fill(fTrack->Pt()/fBinWidth,(Double_t)q,G[5][q]*cos(fgkMltpl*fgkFlow*fTrack->Phi())/(1.+(2.*fR0*sqrt(5.+1.)/nSelTracksIntFlow)*cos(fgkFlow*fTrack->Phi()-2.*q*TMath::Pi()/fgkQmax)),1.);
     //imaginary part
     fDiffFlowGenFunIm5->Fill(fTrack->Pt()/fBinWidth,(Double_t)q,G[5][q]*sin(fgkMltpl*fgkFlow*fTrack->Phi())/(1.+(2.*fR0*sqrt(5.+1.)/nSelTracksIntFlow)*cos(fgkFlow*fTrack->Phi()-2.*q*TMath::Pi()/fgkQmax)),1.);  
     //-----------------------------------------------------------------------
     //real part
     fDiffFlowGenFunRe6->Fill(fTrack->Pt()/fBinWidth,(Double_t)q,G[6][q]*cos(fgkMltpl*fgkFlow*fTrack->Phi())/(1.+(2.*fR0*sqrt(6.+1.)/nSelTracksIntFlow)*cos(fgkFlow*fTrack->Phi()-2.*q*TMath::Pi()/fgkQmax)),1.);
     //imaginary part
     fDiffFlowGenFunIm6->Fill(fTrack->Pt()/fBinWidth,(Double_t)q,G[6][q]*sin(fgkMltpl*fgkFlow*fTrack->Phi())/(1.+(2.*fR0*sqrt(6.+1.)/nSelTracksIntFlow)*cos(fgkFlow*fTrack->Phi()-2.*q*TMath::Pi()/fgkQmax)),1.);
     //-----------------------------------------------------------------------
     //real part
     fDiffFlowGenFunRe7->Fill(fTrack->Pt()/fBinWidth,(Double_t)q,G[7][q]*cos(fgkMltpl*fgkFlow*fTrack->Phi())/(1.+(2.*fR0*sqrt(7.+1.)/nSelTracksIntFlow)*cos(fgkFlow*fTrack->Phi()-2.*q*TMath::Pi()/fgkQmax)),1.);
     //imaginary part
     fDiffFlowGenFunIm7->Fill(fTrack->Pt()/fBinWidth,(Double_t)q,G[7][q]*sin(fgkMltpl*fgkFlow*fTrack->Phi())/(1.+(2.*fR0*sqrt(7.+1.)/nSelTracksIntFlow)*cos(fgkFlow*fTrack->Phi()-2.*q*TMath::Pi()/fgkQmax)),1.);    
   }
   //}
  }  
 }//ending the second loop over data            
 */
      
//}//end of if(nPrim>30)
                      
 
 //numerical equations for cumulants solved up to different highest order  
 
 
 
 
 
  
    
 if(fOtherEquations)
 {
 
  
  
  Double_t G4[fgkPmax4][fgkQmax4];
  Double_t G6[fgkPmax6][fgkQmax6];
  Double_t G8[fgkPmax8][fgkQmax8];
  Double_t G16[fgkPmax16][fgkQmax16];
  for(Int_t p=0;p<fgkPmax16;p++)
  {
   for(Int_t q=0;q<fgkQmax16;q++)
   {
    G16[p][q]=1.;
    if(p<fgkPmax8 && q<fgkQmax8)
    {
     G8[p][q]=1.;
     if(p<fgkPmax6 && q<fgkQmax6)
     {
      G6[p][q]=1.;
      if(p<fgkPmax4 && q<fgkQmax4)
      {
       G4[p][q]=1.;
      }
     }
    } 
   }
  } 
 
  
 
  //Int_t nEventNSelTracksIntFlowOE = anEvent->GetEventNSelTracksIntFlow();
  
  //first loop over data: evaluating the generating function G[p][q] for integrated flow 
  for(Int_t i=0;i<nPrim;i++)
  {
   fTrack=anEvent->GetTrack(i);
   if(fTrack && fTrack->UseForIntegratedFlow())
   {
    for(Int_t p=0;p<fgkPmax16;p++)
    {
     for(Int_t q=0;q<fgkQmax16;q++)
     {
      G16[p][q]*=(1.+(2.*fR0*sqrt(p+1.)/nEventNSelTracksIntFlow)*cos(fgkFlow*fTrack->Phi()-2.*q*TMath::Pi()/fgkQmax16));     
      if(p<fgkPmax8 && q<fgkQmax8)
      {
       G8[p][q]*=(1.+(2.*fR0*sqrt(p+1.)/nEventNSelTracksIntFlow)*cos(fgkFlow*fTrack->Phi()-2.*q*TMath::Pi()/fgkQmax8));
       if(p<fgkPmax6 && q<fgkQmax6)
       {
        G6[p][q]*=(1.+(2.*fR0*sqrt(p+1.)/nEventNSelTracksIntFlow)*cos(fgkFlow*fTrack->Phi()-2.*q*TMath::Pi()/fgkQmax6));
        if(p<fgkPmax4 && q<fgkQmax4)
        {
         G4[p][q]*=(1.+(2.*fR0*sqrt(p+1.)/nEventNSelTracksIntFlow)*cos(fgkFlow*fTrack->Phi()-2.*q*TMath::Pi()/fgkQmax4));
        }
       }
      }
     } 
    }  
   }//end of if(fTrack && fTrack->UseForIntegratedFlow())
  }//ending the loop over data
 
 //storing the value of G[p][q] in 2D profile in order to get automatically the avarage <G[p][q]>
 for(Int_t p=0;p<fgkPmax16;p++)
 {
  for(Int_t q=0;q<fgkQmax16;q++)
  {
   fIntFlowGenFun16->Fill((Double_t)p,(Double_t)q,G16[p][q],1);
   if(p<fgkPmax8 && q<fgkQmax8)
   {
    fIntFlowGenFun8->Fill((Double_t)p,(Double_t)q,G8[p][q],1);
    if(p<fgkPmax6 && q<fgkQmax6)
    {
     fIntFlowGenFun6->Fill((Double_t)p,(Double_t)q,G6[p][q],1);
     if(p<fgkPmax4 && q<fgkQmax4)
     {
      fIntFlowGenFun4->Fill((Double_t)p,(Double_t)q,G4[p][q],1);
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

 AliCumulantsFunctions finalResults(fIntFlowGenFun,fIntFlowGenFun4,fIntFlowGenFun6,fIntFlowGenFun8,fIntFlowGenFun16,fDiffFlowGenFunRe,fDiffFlowGenFunIm,fBinNoOfParticles, fIntFlowResultsGFC,fDiffFlowResults2ndOrderGFC,fDiffFlowResults4thOrderGFC,fDiffFlowResults6thOrderGFC,fDiffFlowResults8thOrderGFC, fAvMultIntFlowGFC,fQVectorComponentsGFC,fCommonHistsResults2nd, fCommonHistsResults4th,fCommonHistsResults6th,fCommonHistsResults8th);
                           
 finalResults.Calculate();  
}

//================================================================================================================

void AliFlowAnalysisWithCumulants::WriteHistograms(TString* outputFileName)
{
 //store the final results in output .root file
 TFile *output = new TFile(outputFileName->Data(),"RECREATE");
 output->mkdir("cobjGFC","cobjGFC");
 output->cd("cobjGFC");
 fHistList->Write(); 
 delete output;
}

//================================================================================================================




