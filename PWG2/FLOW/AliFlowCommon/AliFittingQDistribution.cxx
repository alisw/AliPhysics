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
 * integrated flow estimate by  *
 *   fitting q-distribution     * 
 *                              *
 * author: Ante Bilandzic       * 
 *          (anteb@nikhef.nl)   *
 *******************************/ 

#define AliFittingQDistribution_cxx

#include "Riostream.h"
#include "AliFlowCommonConstants.h"
#include "AliFlowCommonHist.h"
#include "AliFlowCommonHistResults.h"
#include "TChain.h"
#include "TFile.h"
#include "TList.h" 
#include "TParticle.h"
#include "TProfile.h"
#include "AliFlowEventSimple.h"
#include "AliFlowTrackSimple.h"
#include "AliFittingQDistribution.h"
#include "AliFittingFunctionsForQDistribution.h"

class TH1;
class TGraph;
class TPave;
class TLatex;
class TMarker;
class TObjArray;
class TList;
class TCanvas;
class TSystem;
class TROOT;
class AliFlowVector;
class TVector;

//================================================================================================================

ClassImp(AliFittingQDistribution)

AliFittingQDistribution::AliFittingQDistribution():  
 fTrack(NULL),
 fHistList(NULL),
 fWeightsList(NULL),
 fAvMultIntFlowFQD(NULL),
 fIntFlowResultsFQD(NULL),
 fSigma2(NULL),
 fCommonHists(NULL),
 fCommonHistsResults(NULL),
 fQDistributionFQD(NULL),
 fUsePhiWeights(kFALSE)
{
 //constructor 
 fHistList = new TList();
 fWeightsList = new TList(); 
}

AliFittingQDistribution::~AliFittingQDistribution()
{
 //desctructor
 delete fHistList; 
 delete fWeightsList;
}

//================================================================================================================

void AliFittingQDistribution::Init()
{
 //various output histograms
 
 //avarage multiplicity 
 fAvMultIntFlowFQD = new TProfile("fAvMultIntFlowFQD","Average Multiplicity",1,0,1,"s");
 fAvMultIntFlowFQD->SetXTitle("");
 fAvMultIntFlowFQD->SetYTitle("");
 fAvMultIntFlowFQD->SetLabelSize(0.06);
 fAvMultIntFlowFQD->SetMarkerStyle(25);
 fAvMultIntFlowFQD->SetLabelOffset(0.02);
 (fAvMultIntFlowFQD->GetXaxis())->SetBinLabel(1,"Average Multiplicity");
 fHistList->Add(fAvMultIntFlowFQD);
 
 //final result for integrated flow 
 fIntFlowResultsFQD = new TH1D("fIntFlowResultsFQD","Integrated Flow By Fitting q-distribution",1,0,1);
 fIntFlowResultsFQD->SetXTitle("");
 fIntFlowResultsFQD->SetYTitle("");
 fIntFlowResultsFQD->SetMarkerStyle(25);
 fIntFlowResultsFQD->SetLabelSize(0.06);
 fIntFlowResultsFQD->SetLabelOffset(0.02);
 (fIntFlowResultsFQD->GetXaxis())->SetBinLabel(1,"v_{n}{FQD}");
 fHistList->Add(fIntFlowResultsFQD);
 
 //sigma^2
 fSigma2 = new TH1D("fSigma2","#sigma^{2}",1,0,1);
 fSigma2->SetXTitle("");
 fSigma2->SetYTitle("");
 fSigma2->SetMarkerStyle(25);
 fSigma2->SetLabelSize(0.06);
 fSigma2->SetLabelOffset(0.02);
 (fSigma2->GetXaxis())->SetBinLabel(1,"#sigma^{2}");
 fHistList->Add(fSigma2);
 
 //q-distribution 
 fQDistributionFQD = new TH1D("fQDistributionFQD","q-distribution",100,0,10);
 fQDistributionFQD->SetXTitle("q_{n}=Q_{n}/#sqrt{M}");
 fQDistributionFQD->SetYTitle("Counts");
 fHistList->Add(fQDistributionFQD);
  
 //common control histograms
 fCommonHists = new AliFlowCommonHist("AliFlowCommonHistFQD");
 fHistList->Add(fCommonHists);  
 
 //common histograms for final results (2nd order)
 fCommonHistsResults= new AliFlowCommonHistResults("AliFlowCommonHistResultsFQD");
 fHistList->Add(fCommonHistsResults); 
 
}//end of Init()

//================================================================================================================

void AliFittingQDistribution::Make(AliFlowEventSimple* anEvent)
{
 
 //Int_t nPrim = anEvent->NumberOfTracks();//total multiplicity
  
 Int_t n=2;//harmonic (to be improved)  
   
 //fill the common control histograms
 fCommonHists->FillControlHistograms(anEvent);   

 //calculating Q-vector of event
 AliFlowVector fQVector;
 fQVector.Set(0.,0.);
 fQVector.SetMult(0);
 fQVector=anEvent->GetQ(n,fWeightsList,fUsePhiWeights);                                                                                  
                                                                                                                                                                      
 //multiplicity
 fAvMultIntFlowFQD->Fill(0.,fQVector.GetMult(),1.);
 
 //q = Q/sqrt(M)
 Double_t q=0.;
 
 if(fQVector.GetMult()!=0)
 {
  q = fQVector.Mod()/sqrt(fQVector.GetMult());
  fQDistributionFQD->Fill(q,1.);
 }  
}//end of Make()

//================================================================================================================

void AliFittingQDistribution::Finish()
{
 //calculate the final results
 AliFittingFunctionsForQDistribution finalFitting(fAvMultIntFlowFQD,fQDistributionFQD,fIntFlowResultsFQD,fSigma2,fCommonHistsResults);
         
 finalFitting.Calculate();            
}

//================================================================================================================

void AliFittingQDistribution::WriteHistograms(TString* outputFileName)
{
 //store the final results in output .root file
 TFile *output = new TFile(outputFileName->Data(),"RECREATE");
 output->WriteObject(fHistList, "cobjFQD","SingleKey");
 delete output;
}

//================================================================================================================

void AliFittingQDistribution::WriteHistograms(TString outputFileName)
{
 //store the final results in output .root file
 TFile *output = new TFile(outputFileName.Data(),"RECREATE");
 output->WriteObject(fHistList, "cobjFQD","SingleKey");
 delete output;
}

//================================================================================================================


















