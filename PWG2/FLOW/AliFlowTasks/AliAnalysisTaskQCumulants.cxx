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

/**************************************
 * analysis task for Q-cumulants      * 
 *                                    * 
 * authors: Naomi van der Kolk        *
 *           (kolk@nikhef.nl)         *  
 *          Raimond Snellings         *
 *           (snelling@nikhef.nl)     * 
 *          Ante Bilandzic            *
 *           (anteb@nikhef.nl)        * 
 * ***********************************/
 
#include "Riostream.h"
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TList.h"
#include "TH1.h"
#include "TGraph.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TProfile3D.h"
#include "TBits.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisDataSlot.h"
#include "AliAnalysisDataContainer.h"
#include "AliAnalysisManager.h"

#include "AliFlowEventSimple.h"
#include "AliAnalysisTaskQCumulants.h"
#include "AliFlowAnalysisWithQCumulants.h"
#include "AliFlowCumuConstants.h"
#include "AliFlowCommonConstants.h"
#include "AliFlowCommonHist.h"
#include "AliFlowCommonHistResults.h"

ClassImp(AliAnalysisTaskQCumulants)

//================================================================================================================

AliAnalysisTaskQCumulants::AliAnalysisTaskQCumulants(const char *name, Bool_t useWeights): 
 AliAnalysisTask(name,""), 
 fEvent(NULL),
 fQCA(NULL), // Q-cumulant Analysis (QCA) object
 fListHistos(NULL),
 fUseWeights(useWeights),
 fUsePhiWeights(kFALSE),
 fUsePtWeights(kFALSE),
 fUseEtaWeights(kFALSE),
 fListWeights(NULL)
{
 // constructor
 cout<<"AliAnalysisTaskQCumulants::AliAnalysisTaskQCumulants(const char *name)"<<endl;
 
 // Define input and output slots here
 // Input slot #0 works with a TChain
 DefineInput(0, AliFlowEventSimple::Class());
  
 // Input slot #1 is needed for the weights 
 if(useWeights)
 {
  DefineInput(1, TList::Class());   
 }
        
 // Output slot #0 writes into a TList container
 DefineOutput(0, TList::Class());  
 
}

AliAnalysisTaskQCumulants::AliAnalysisTaskQCumulants(): 
 fEvent(NULL),
 fQCA(NULL),//Q-cumulant Analysis (QCA) object
 fListHistos(NULL),
 fUseWeights(kFALSE),
 fUsePhiWeights(kFALSE),
 fUsePtWeights(kFALSE),
 fUseEtaWeights(kFALSE),
 fListWeights(NULL)
{
 // dummy constructor
 cout<<"AliAnalysisTaskQCumulants::AliAnalysisTaskQCumulants()"<<endl;
}

//================================================================================================================

void AliAnalysisTaskQCumulants::ConnectInputData(Option_t *) 
{
 // connect ESD or AOD (called once)
 cout<<"AliAnalysisTaskQCumulants::ConnectInputData(Option_t *)"<<endl;
}

//================================================================================================================

void AliAnalysisTaskQCumulants::CreateOutputObjects() 
{
 // called at every worker node to initialize
 cout<<"AliAnalysisTaskQCumulants::CreateOutputObjects()"<<endl;

 // analyser
 fQCA = new AliFlowAnalysisWithQCumulants();
 fQCA->Init();
 
 //weights:
 if(fUseWeights)
 {
  //pass the flags to class:
  if(fUsePhiWeights) fQCA->SetUsePhiWeights(fUsePhiWeights);
  if(fUsePtWeights) fQCA->SetUsePtWeights(fUsePtWeights);
  if(fUseEtaWeights) fQCA->SetUseEtaWeights(fUseEtaWeights);
  //get data from input slot #1 which is used for weights:
  if(GetNinputs()==2) 
  {                   
   fListWeights = (TList*)GetInputData(1); 
  }
  //pass the list with weights to class:
  if(fListWeights) fQCA->SetWeightsList(fListWeights);
 }
 
 if(fQCA->GetHistList()) 
 {
  fListHistos = fQCA->GetHistList();
  //fListHistos->Print();
 }
 else 
 {
  Printf(" ERROR: Could not retrieve histogram list (QC, Task::COO)"); 
 }

 //PostData(0,fListHistos);
 
}

//================================================================================================================

void AliAnalysisTaskQCumulants::Exec(Option_t *) 
{
  // main loop (called for each event)
 fEvent = dynamic_cast<AliFlowEventSimple*>(GetInputData(0));

 // Q-cumulants
 if(fEvent) 
 {
  fQCA->Make(fEvent);
 }else 
  {
   cout<<" WARNING: No input data (QC, Task::E) !!!"<<endl;
   cout<<endl;
  }
  
 PostData(0,fListHistos);
}

//================================================================================================================

void AliAnalysisTaskQCumulants::Terminate(Option_t *) 
{
 //accessing the output list which contains the merged 2D and 3D profiles from all worker nodes
 fListHistos = (TList*)GetOutputData(0);
 //fListHistos->Print();

 if(fListHistos)
 {	
  // with or without weights
  TBits *useWeightsBits = dynamic_cast<TBits*>(fListHistos->FindObject("TBits"));
         
  //final results (no-name integrated flow without weights)
  TH1D *intFlowResultsQC = dynamic_cast<TH1D*>((dynamic_cast<TList*>(fListHistos->FindObject("Results")))->FindObject("fIntFlowResultsQC"));
  
  //final results (no-name integrated flow with weights)
  TH1D *intFlowResultsQCW = dynamic_cast<TH1D*>((dynamic_cast<TList*>(fListHistos->FindObject("Results")))->FindObject("fIntFlowResultsQCW"));
  
  //final results (POIs integrated flow without weights)
  TH1D *intFlowResultsPOIQC = dynamic_cast<TH1D*>((dynamic_cast<TList*>(fListHistos->FindObject("Results")))->FindObject("fIntFlowResultsPOIQC"));
  
  //final results (POIs integrated flow with weights)
  TH1D *intFlowResultsPOIQCW = dynamic_cast<TH1D*>((dynamic_cast<TList*>(fListHistos->FindObject("Results")))->FindObject("fIntFlowResultsPOIQCW"));
  
  //final results (RPs integrated flow without weights)
  TH1D *intFlowResultsRPQC = dynamic_cast<TH1D*>((dynamic_cast<TList*>(fListHistos->FindObject("Results")))->FindObject("fIntFlowResultsRPQC"));
  
  //final results (RPs integrated flow with weights)
  TH1D *intFlowResultsRPQCW = dynamic_cast<TH1D*>((dynamic_cast<TList*>(fListHistos->FindObject("Results")))->FindObject("fIntFlowResultsRPQCW"));
  
  //final results (differential flow)
  TH1D *diffFlowResults2ndOrder = dynamic_cast<TH1D*>(fListHistos->FindObject("fDiffFlowResults2ndOrderQC"));
  TH1D *diffFlowResults4thOrder = dynamic_cast<TH1D*>(fListHistos->FindObject("fDiffFlowResults4thOrderQC"));
  
  //final results for covariances (1st bin <2*4>-<2>*<4>, 2nd bin <2*6>-<2>*<6>, ...)
  TH1D *covariances = dynamic_cast<TH1D*>(fListHistos->FindObject("fCovariances"));
  
  //common control histograms (taking into account only the events with 2 and more particles)  
  AliFlowCommonHist *commonHist2nd = dynamic_cast<AliFlowCommonHist*>(fListHistos->FindObject("AliFlowCommonHist2ndOrderQC"));
  
  //common control histograms (taking into account only the events with 4 and more particles)  
  AliFlowCommonHist *commonHist4th = dynamic_cast<AliFlowCommonHist*>(fListHistos->FindObject("AliFlowCommonHist4thOrderQC"));
  
  //common control histograms (taking into account only the events with 6 and more particles)  
  AliFlowCommonHist *commonHist6th = dynamic_cast<AliFlowCommonHist*>(fListHistos->FindObject("AliFlowCommonHist6thOrderQC"));
  
  //common control histograms (taking into account only the events with 8 and more particles)  
  AliFlowCommonHist *commonHist8th = dynamic_cast<AliFlowCommonHist*>(fListHistos->FindObject("AliFlowCommonHist8thOrderQC"));
  
  //common histograms to store the final results for the 2nd order integrated and differential flow
  AliFlowCommonHistResults *commonHistRes2nd = dynamic_cast<AliFlowCommonHistResults*>(fListHistos->FindObject("AliFlowCommonHistResults2ndOrderQC"));
  
  //common histograms to store the final results for the 4th order integrated and differential flow
  AliFlowCommonHistResults *commonHistRes4th = dynamic_cast<AliFlowCommonHistResults*>(fListHistos->FindObject("AliFlowCommonHistResults4thOrderQC"));
  
  //common histograms to store the final results for the 6th order integrated and differential flow
  AliFlowCommonHistResults *commonHistRes6th = dynamic_cast<AliFlowCommonHistResults*>(fListHistos->FindObject("AliFlowCommonHistResults6thOrderQC"));
  
  //common histograms to store the final results for the 8th order integrated and differential flow
  AliFlowCommonHistResults *commonHistRes8th = dynamic_cast<AliFlowCommonHistResults*>(fListHistos->FindObject("AliFlowCommonHistResults8thOrderQC"));
  
  //average selected multiplicity (for int. flow) 
  TProfile *AvMult = dynamic_cast<TProfile*>(fListHistos->FindObject("fAvMultIntFlowQC"));
  
  //xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  //                        !!!! to be removed !!!!
  //profiles containing the Q-vectors from all events 
  TProfile *qvectorForEachEventX = dynamic_cast<TProfile*>(fListHistos->FindObject("fQvectorForEachEventX"));
  TProfile *qvectorForEachEventY = dynamic_cast<TProfile*>(fListHistos->FindObject("fQvectorForEachEventY"));  
  //xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    
  //multi-particle correlations calculated from Q-vectors
  TProfile *qCorrelations = dynamic_cast<TProfile*>(fListHistos->FindObject("fQCorrelations"));
  
  //weighted multi-particle correlations calculated from Q-vectors
  TProfile *qCorrelationsW = dynamic_cast<TProfile*>(fListHistos->FindObject("fQCorrelationsW"));
  
  //average of products: 1st bin: <2*4>, 2nd bin: <2*6>, ...
  TProfile *QProduct = dynamic_cast<TProfile*>(fListHistos->FindObject("fQProduct"));
  
  //average 2- and 4-particle correlations per pt-bin 
  TProfile *binnedPt2p1n1nRP = dynamic_cast<TProfile*>((dynamic_cast<TList*>(fListHistos->FindObject("DifferentialFlow")))->FindObject("f2PerPtBin1n1nRP"));
  TProfile *binnedPt4p1n1n1n1nRP = dynamic_cast<TProfile*>((dynamic_cast<TList*>(fListHistos->FindObject("DifferentialFlow")))->FindObject("f4PerPtBin1n1n1n1nRP"));
  
  //average 2- and 4-particle correlations per eta-bin 
  TProfile *binnedEta2p1n1nRP = dynamic_cast<TProfile*>((dynamic_cast<TList*>(fListHistos->FindObject("DifferentialFlow")))->FindObject("f2PerEtaBin1n1nRP"));
  TProfile *binnedEta4p1n1n1n1nRP = dynamic_cast<TProfile*>((dynamic_cast<TList*>(fListHistos->FindObject("DifferentialFlow")))->FindObject("f4PerEtaBin1n1n1n1nRP"));  
  
  //average 2- and 4-particle correlations per pt-bin 
  TProfile *binnedPt2p1n1nPOI = dynamic_cast<TProfile*>((dynamic_cast<TList*>(fListHistos->FindObject("DifferentialFlow")))->FindObject("f2PerPtBin1n1nPOI"));
  TProfile *binnedPt4p1n1n1n1nPOI = dynamic_cast<TProfile*>((dynamic_cast<TList*>(fListHistos->FindObject("DifferentialFlow")))->FindObject("f4PerPtBin1n1n1n1nPOI"));
  
  //average 2- and 4-particle correlations per eta-bin 
  TProfile *binnedEta2p1n1nPOI = dynamic_cast<TProfile*>((dynamic_cast<TList*>(fListHistos->FindObject("DifferentialFlow")))->FindObject("f2PerEtaBin1n1nPOI"));
  TProfile *binnedEta4p1n1n1n1nPOI = dynamic_cast<TProfile*>((dynamic_cast<TList*>(fListHistos->FindObject("DifferentialFlow")))->FindObject("f4PerEtaBin1n1n1n1nPOI"));  
  
  //average 2- and 4-particle correlations per pt-bin 
  TProfile *binnedWPt2p1n1nPOI = dynamic_cast<TProfile*>((dynamic_cast<TList*>(fListHistos->FindObject("DifferentialFlow")))->FindObject("f2WPerPtBin1n1nPOI"));
  TProfile *binnedWPt4p1n1n1n1nPOI = dynamic_cast<TProfile*>((dynamic_cast<TList*>(fListHistos->FindObject("DifferentialFlow")))->FindObject("f4WPerPtBin1n1n1n1nPOI"));
 
  TProfile *binnedWEta2p1n1nPOI = dynamic_cast<TProfile*>((dynamic_cast<TList*>(fListHistos->FindObject("DifferentialFlow")))->FindObject("f2WPerEtaBin1n1nPOI"));
  TProfile *binnedWEta4p1n1n1n1nPOI = dynamic_cast<TProfile*>((dynamic_cast<TList*>(fListHistos->FindObject("DifferentialFlow")))->FindObject("f4WPerEtaBin1n1n1n1nPOI"));
  
  TProfile *binnedWPt2p1n1nRP = dynamic_cast<TProfile*>((dynamic_cast<TList*>(fListHistos->FindObject("DifferentialFlow")))->FindObject("f2WPerPtBin1n1nRP"));
  TProfile *binnedWPt4p1n1n1n1nRP = dynamic_cast<TProfile*>((dynamic_cast<TList*>(fListHistos->FindObject("DifferentialFlow")))->FindObject("f4WPerPtBin1n1n1n1nRP"));
  
  TProfile *binnedWEta2p1n1nRP = dynamic_cast<TProfile*>((dynamic_cast<TList*>(fListHistos->FindObject("DifferentialFlow")))->FindObject("f2WPerEtaBin1n1nRP"));
  TProfile *binnedWEta4p1n1n1n1nRP = dynamic_cast<TProfile*>((dynamic_cast<TList*>(fListHistos->FindObject("DifferentialFlow")))->FindObject("f4WPerEtaBin1n1n1n1nRP"));
  
  //average values of Q-vector components (1st bin: <Q_x>, 2nd bin: <Q_y>, 3rd bin: <(Q_x)^2>, 4th bin: <(Q_y)^2>) 
  TProfile *QVectorComponents = dynamic_cast<TProfile*>(fListHistos->FindObject("fQvectorComponents"));
  
  // multi-particle correlations calculated with nested loop (needed for int. flow)
  TProfile *directCorrelations = dynamic_cast<TProfile*>(fListHistos->FindObject("fDirectCorrelations"));
  
  // multi-particle correlations calculated with nested loop (needed for weighted int. flow)
  TProfile *directCorrelationsW = dynamic_cast<TProfile*>(fListHistos->FindObject("fDirectCorrelationsW"));
  
  // multi-particle correlations calculated with nested loop (needed for diff. flow)
  TProfile *directCorrelationsDiffFlow = dynamic_cast<TProfile*>(fListHistos->FindObject("fDirectCorrelationsDiffFlow"));
  
  // multi-particle correlations calculated with nested loop (needed for int. flow)
  TProfile *directCorrelationsDiffFlowW = dynamic_cast<TProfile*>(fListHistos->FindObject("fDirectCorrelationsDiffFlowW"));
  
  
  
  
  
  
  // ...............................................................................................................................................
  // non-weighted correlations for each (pt,eta) bin for POIs:
  TProfile2D *twoPtEtaPOI   = dynamic_cast<TProfile2D*>((dynamic_cast<TList*>(fListHistos->FindObject("DifferentialFlow")))->FindObject("f2pPtEtaPOI"));
  TProfile2D *fourPtEtaPOI  = dynamic_cast<TProfile2D*>((dynamic_cast<TList*>(fListHistos->FindObject("DifferentialFlow")))->FindObject("f4pPtEtaPOI"));
  TProfile2D *sixPtEtaPOI   = dynamic_cast<TProfile2D*>((dynamic_cast<TList*>(fListHistos->FindObject("DifferentialFlow")))->FindObject("f6pPtEtaPOI"));
  TProfile2D *eightPtEtaPOI = dynamic_cast<TProfile2D*>((dynamic_cast<TList*>(fListHistos->FindObject("DifferentialFlow")))->FindObject("f8pPtEtaPOI"));
 
  // non-weighted final results for differential flow for each for POIs:
  // 3D (pt,eta)
  TH2D *vn2ndPtEtaPOI = dynamic_cast<TH2D*>((dynamic_cast<TList*>(fListHistos->FindObject("Results")))->FindObject("fvn2ndPtEtaPOI"));  
  TH2D *vn4thPtEtaPOI = dynamic_cast<TH2D*>((dynamic_cast<TList*>(fListHistos->FindObject("Results")))->FindObject("fvn4thPtEtaPOI")); 
  TH2D *vn6thPtEtaPOI = dynamic_cast<TH2D*>((dynamic_cast<TList*>(fListHistos->FindObject("Results")))->FindObject("fvn6thPtEtaPOI")); 
  TH2D *vn8thPtEtaPOI = dynamic_cast<TH2D*>((dynamic_cast<TList*>(fListHistos->FindObject("Results")))->FindObject("fvn8thPtEtaPOI")); 
  // 2D (pt)
  TH1D *vn2ndPtPOI = dynamic_cast<TH1D*>((dynamic_cast<TList*>(fListHistos->FindObject("Results")))->FindObject("fvn2ndPtPOI"));  
  TH1D *vn4thPtPOI = dynamic_cast<TH1D*>((dynamic_cast<TList*>(fListHistos->FindObject("Results")))->FindObject("fvn4thPtPOI")); 
  TH1D *vn6thPtPOI = dynamic_cast<TH1D*>((dynamic_cast<TList*>(fListHistos->FindObject("Results")))->FindObject("fvn6thPtPOI")); 
  TH1D *vn8thPtPOI = dynamic_cast<TH1D*>((dynamic_cast<TList*>(fListHistos->FindObject("Results")))->FindObject("fvn8thPtPOI"));
  // 2D (eta)
  TH1D *vn2ndEtaPOI = dynamic_cast<TH1D*>((dynamic_cast<TList*>(fListHistos->FindObject("Results")))->FindObject("fvn2ndEtaPOI"));  
  TH1D *vn4thEtaPOI = dynamic_cast<TH1D*>((dynamic_cast<TList*>(fListHistos->FindObject("Results")))->FindObject("fvn4thEtaPOI")); 
  TH1D *vn6thEtaPOI = dynamic_cast<TH1D*>((dynamic_cast<TList*>(fListHistos->FindObject("Results")))->FindObject("fvn6thEtaPOI")); 
  TH1D *vn8thEtaPOI = dynamic_cast<TH1D*>((dynamic_cast<TList*>(fListHistos->FindObject("Results")))->FindObject("fvn8thEtaPOI"));
  
  // weighted correlations for each (pt,eta) bin for POIs:
  TProfile2D *twoPtEtaPOIW   = dynamic_cast<TProfile2D*>((dynamic_cast<TList*>(fListHistos->FindObject("DifferentialFlow")))->FindObject("f2pPtEtaPOIW"));
  TProfile2D *fourPtEtaPOIW  = dynamic_cast<TProfile2D*>((dynamic_cast<TList*>(fListHistos->FindObject("DifferentialFlow")))->FindObject("f4pPtEtaPOIW"));
  TProfile2D *sixPtEtaPOIW   = dynamic_cast<TProfile2D*>((dynamic_cast<TList*>(fListHistos->FindObject("DifferentialFlow")))->FindObject("f6pPtEtaPOIW"));
  TProfile2D *eightPtEtaPOIW = dynamic_cast<TProfile2D*>((dynamic_cast<TList*>(fListHistos->FindObject("DifferentialFlow")))->FindObject("f8pPtEtaPOIW"));
  
  // weighted final results for differential flow for each for POIs:
  // 3D (pt,eta)
  TH2D *vn2ndPtEtaPOIW = dynamic_cast<TH2D*>((dynamic_cast<TList*>(fListHistos->FindObject("Results")))->FindObject("fvn2ndPtEtaPOIW"));  
  TH2D *vn4thPtEtaPOIW = dynamic_cast<TH2D*>((dynamic_cast<TList*>(fListHistos->FindObject("Results")))->FindObject("fvn4thPtEtaPOIW")); 
  TH2D *vn6thPtEtaPOIW = dynamic_cast<TH2D*>((dynamic_cast<TList*>(fListHistos->FindObject("Results")))->FindObject("fvn6thPtEtaPOIW")); 
  TH2D *vn8thPtEtaPOIW = dynamic_cast<TH2D*>((dynamic_cast<TList*>(fListHistos->FindObject("Results")))->FindObject("fvn8thPtEtaPOIW")); 
  // 2D (pt)
  TH1D *vn2ndPtPOIW = dynamic_cast<TH1D*>((dynamic_cast<TList*>(fListHistos->FindObject("Results")))->FindObject("fvn2ndPtPOIW"));  
  TH1D *vn4thPtPOIW = dynamic_cast<TH1D*>((dynamic_cast<TList*>(fListHistos->FindObject("Results")))->FindObject("fvn4thPtPOIW")); 
  TH1D *vn6thPtPOIW = dynamic_cast<TH1D*>((dynamic_cast<TList*>(fListHistos->FindObject("Results")))->FindObject("fvn6thPtPOIW")); 
  TH1D *vn8thPtPOIW = dynamic_cast<TH1D*>((dynamic_cast<TList*>(fListHistos->FindObject("Results")))->FindObject("fvn8thPtPOIW"));
  // 2D (eta)
  TH1D *vn2ndEtaPOIW = dynamic_cast<TH1D*>((dynamic_cast<TList*>(fListHistos->FindObject("Results")))->FindObject("fvn2ndEtaPOIW"));  
  TH1D *vn4thEtaPOIW = dynamic_cast<TH1D*>((dynamic_cast<TList*>(fListHistos->FindObject("Results")))->FindObject("fvn4thEtaPOIW")); 
  TH1D *vn6thEtaPOIW = dynamic_cast<TH1D*>((dynamic_cast<TList*>(fListHistos->FindObject("Results")))->FindObject("fvn6thEtaPOIW")); 
  TH1D *vn8thEtaPOIW = dynamic_cast<TH1D*>((dynamic_cast<TList*>(fListHistos->FindObject("Results")))->FindObject("fvn8thEtaPOIW"));
  
  // non-weighted correlations for each (pt,eta) bin for RPs:
  TProfile2D *twoPtEtaRP   = dynamic_cast<TProfile2D*>((dynamic_cast<TList*>(fListHistos->FindObject("DifferentialFlow")))->FindObject("f2pPtEtaRP"));
  TProfile2D *fourPtEtaRP  = dynamic_cast<TProfile2D*>((dynamic_cast<TList*>(fListHistos->FindObject("DifferentialFlow")))->FindObject("f4pPtEtaRP"));
  TProfile2D *sixPtEtaRP   = dynamic_cast<TProfile2D*>((dynamic_cast<TList*>(fListHistos->FindObject("DifferentialFlow")))->FindObject("f6pPtEtaRP"));
  TProfile2D *eightPtEtaRP = dynamic_cast<TProfile2D*>((dynamic_cast<TList*>(fListHistos->FindObject("DifferentialFlow")))->FindObject("f8pPtEtaRP"));
 
  // non-weighted final results for differential flow for RPs:
  // 3D (pt,eta)
  TH2D *vn2ndPtEtaRP = dynamic_cast<TH2D*>((dynamic_cast<TList*>(fListHistos->FindObject("Results")))->FindObject("fvn2ndPtEtaRP"));  
  TH2D *vn4thPtEtaRP = dynamic_cast<TH2D*>((dynamic_cast<TList*>(fListHistos->FindObject("Results")))->FindObject("fvn4thPtEtaRP")); 
  TH2D *vn6thPtEtaRP = dynamic_cast<TH2D*>((dynamic_cast<TList*>(fListHistos->FindObject("Results")))->FindObject("fvn6thPtEtaRP")); 
  TH2D *vn8thPtEtaRP = dynamic_cast<TH2D*>((dynamic_cast<TList*>(fListHistos->FindObject("Results")))->FindObject("fvn8thPtEtaRP")); 
  // 2D (pt)
  TH1D *vn2ndPtRP = dynamic_cast<TH1D*>((dynamic_cast<TList*>(fListHistos->FindObject("Results")))->FindObject("fvn2ndPtRP"));  
  TH1D *vn4thPtRP = dynamic_cast<TH1D*>((dynamic_cast<TList*>(fListHistos->FindObject("Results")))->FindObject("fvn4thPtRP")); 
  TH1D *vn6thPtRP = dynamic_cast<TH1D*>((dynamic_cast<TList*>(fListHistos->FindObject("Results")))->FindObject("fvn6thPtRP")); 
  TH1D *vn8thPtRP = dynamic_cast<TH1D*>((dynamic_cast<TList*>(fListHistos->FindObject("Results")))->FindObject("fvn8thPtRP"));
  // 2D (eta)
  TH1D *vn2ndEtaRP = dynamic_cast<TH1D*>((dynamic_cast<TList*>(fListHistos->FindObject("Results")))->FindObject("fvn2ndEtaRP"));  
  TH1D *vn4thEtaRP = dynamic_cast<TH1D*>((dynamic_cast<TList*>(fListHistos->FindObject("Results")))->FindObject("fvn4thEtaRP")); 
  TH1D *vn6thEtaRP = dynamic_cast<TH1D*>((dynamic_cast<TList*>(fListHistos->FindObject("Results")))->FindObject("fvn6thEtaRP")); 
  TH1D *vn8thEtaRP = dynamic_cast<TH1D*>((dynamic_cast<TList*>(fListHistos->FindObject("Results")))->FindObject("fvn8thEtaRP"));
  
  // weighted correlations for each (pt,eta) bin for RPs:
  TProfile2D *twoPtEtaRPW   = dynamic_cast<TProfile2D*>((dynamic_cast<TList*>(fListHistos->FindObject("DifferentialFlow")))->FindObject("f2pPtEtaRPW"));
  TProfile2D *fourPtEtaRPW  = dynamic_cast<TProfile2D*>((dynamic_cast<TList*>(fListHistos->FindObject("DifferentialFlow")))->FindObject("f4pPtEtaRPW"));
  TProfile2D *sixPtEtaRPW   = dynamic_cast<TProfile2D*>((dynamic_cast<TList*>(fListHistos->FindObject("DifferentialFlow")))->FindObject("f6pPtEtaRPW"));
  TProfile2D *eightPtEtaRPW = dynamic_cast<TProfile2D*>((dynamic_cast<TList*>(fListHistos->FindObject("DifferentialFlow")))->FindObject("f8pPtEtaRPW"));
  
  // weighted final results for differential flow for RPs:
  // 3D (pt,eta)
  TH2D *vn2ndPtEtaRPW = dynamic_cast<TH2D*>((dynamic_cast<TList*>(fListHistos->FindObject("Results")))->FindObject("fvn2ndPtEtaRPW"));  
  TH2D *vn4thPtEtaRPW = dynamic_cast<TH2D*>((dynamic_cast<TList*>(fListHistos->FindObject("Results")))->FindObject("fvn4thPtEtaRPW")); 
  TH2D *vn6thPtEtaRPW = dynamic_cast<TH2D*>((dynamic_cast<TList*>(fListHistos->FindObject("Results")))->FindObject("fvn6thPtEtaRPW")); 
  TH2D *vn8thPtEtaRPW = dynamic_cast<TH2D*>((dynamic_cast<TList*>(fListHistos->FindObject("Results")))->FindObject("fvn8thPtEtaRPW")); 
  // 2D (pt)
  TH1D *vn2ndPtRPW = dynamic_cast<TH1D*>((dynamic_cast<TList*>(fListHistos->FindObject("Results")))->FindObject("fvn2ndPtRPW"));  
  TH1D *vn4thPtRPW = dynamic_cast<TH1D*>((dynamic_cast<TList*>(fListHistos->FindObject("Results")))->FindObject("fvn4thPtRPW")); 
  TH1D *vn6thPtRPW = dynamic_cast<TH1D*>((dynamic_cast<TList*>(fListHistos->FindObject("Results")))->FindObject("fvn6thPtRPW")); 
  TH1D *vn8thPtRPW = dynamic_cast<TH1D*>((dynamic_cast<TList*>(fListHistos->FindObject("Results")))->FindObject("fvn8thPtRPW"));
  // 2D (eta)
  TH1D *vn2ndEtaRPW = dynamic_cast<TH1D*>((dynamic_cast<TList*>(fListHistos->FindObject("Results")))->FindObject("fvn2ndEtaRPW"));  
  TH1D *vn4thEtaRPW = dynamic_cast<TH1D*>((dynamic_cast<TList*>(fListHistos->FindObject("Results")))->FindObject("fvn4thEtaRPW")); 
  TH1D *vn6thEtaRPW = dynamic_cast<TH1D*>((dynamic_cast<TList*>(fListHistos->FindObject("Results")))->FindObject("fvn6thEtaRPW")); 
  TH1D *vn8thEtaRPW = dynamic_cast<TH1D*>((dynamic_cast<TList*>(fListHistos->FindObject("Results")))->FindObject("fvn8thEtaRPW"));
  // ...............................................................................................................................................
  
  
 
  //----------------------------------------------------
 
  fQCA = new AliFlowAnalysisWithQCumulants();  
 
  fQCA->SetUseWeightsBits(useWeightsBits);
  fQCA->SetIntFlowResults(intFlowResultsQC); 
  fQCA->SetIntFlowResultsW(intFlowResultsQCW);
  fQCA->SetIntFlowResultsPOI(intFlowResultsPOIQC); 
  fQCA->SetIntFlowResultsPOIW(intFlowResultsPOIQCW); 
  fQCA->SetIntFlowResultsRP(intFlowResultsRPQC); 
  fQCA->SetIntFlowResultsRPW(intFlowResultsRPQCW); 

  fQCA->SetDiffFlowResults2nd(diffFlowResults2ndOrder);
  fQCA->SetDiffFlowResults4th(diffFlowResults4thOrder); 
  fQCA->SetCovariances(covariances); 
  
  fQCA->SetCommonHists2nd(commonHist2nd); 
  fQCA->SetCommonHists4th(commonHist4th);
  fQCA->SetCommonHists6th(commonHist6th);
  fQCA->SetCommonHists8th(commonHist8th);

  fQCA->SetCommonHistsResults2nd(commonHistRes2nd); 
  fQCA->SetCommonHistsResults4th(commonHistRes4th);
  fQCA->SetCommonHistsResults6th(commonHistRes6th);
  fQCA->SetCommonHistsResults8th(commonHistRes8th);
 
  fQCA->SetAverageMultiplicity(AvMult);
  //xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  //             !!!! to be removed !!!!
  fQCA->SetQvectorForEachEventX(qvectorForEachEventX);
  fQCA->SetQvectorForEachEventY(qvectorForEachEventY);
  //xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  fQCA->SetQCorrelations(qCorrelations);
  fQCA->SetQCorrelationsW(qCorrelationsW);
  fQCA->SetQProduct(QProduct);
  fQCA->SetQVectorComponents(QVectorComponents);
 
  fQCA->SetTwo1n1nPerPtBinRP(binnedPt2p1n1nRP);
  fQCA->SetFour1n1n1n1nPerPtBinRP(binnedPt4p1n1n1n1nRP);
  
  fQCA->SetTwo1n1nPerEtaBinRP(binnedEta2p1n1nRP);
  fQCA->SetFour1n1n1n1nPerEtaBinRP(binnedEta4p1n1n1n1nRP); 
  
  fQCA->SetTwo1n1nPerPtBinPOI(binnedPt2p1n1nPOI);
  fQCA->SetFour1n1n1n1nPerPtBinPOI(binnedPt4p1n1n1n1nPOI);
  
  fQCA->SetTwo1n1nPerEtaBinPOI(binnedEta2p1n1nPOI);
  fQCA->SetFour1n1n1n1nPerEtaBinPOI(binnedEta4p1n1n1n1nPOI);    
 
  fQCA->SetTwo1n1nWPerPtBinPOI(binnedWPt2p1n1nPOI);
  fQCA->SetFour1n1n1n1nWPerPtBinPOI(binnedWPt4p1n1n1n1nPOI);
  
  fQCA->SetTwo1n1nWPerEtaBinPOI(binnedWEta2p1n1nPOI);
  fQCA->SetFour1n1n1n1nWPerEtaBinPOI(binnedWEta4p1n1n1n1nPOI);
  
  fQCA->SetTwo1n1nWPerPtBinRP(binnedWPt2p1n1nRP);
  fQCA->SetFour1n1n1n1nWPerPtBinRP(binnedWPt4p1n1n1n1nRP);
  
  fQCA->SetTwo1n1nWPerEtaBinRP(binnedWEta2p1n1nRP);
  fQCA->SetFour1n1n1n1nWPerEtaBinRP(binnedWEta4p1n1n1n1nRP);
  
  // nested loops results:
  fQCA->SetDirectCorrelations(directCorrelations);
  fQCA->SetDirectCorrelationsW(directCorrelationsW);
  fQCA->SetDirectCorrelationsDiffFlow(directCorrelationsDiffFlow);
  fQCA->SetDirectCorrelationsDiffFlowW(directCorrelationsDiffFlowW);
  
  // non-weighted correlations for each (pt,eta) bin for POIs:
  fQCA->Set2pPtEtaPOI(twoPtEtaPOI);
  fQCA->Set4pPtEtaPOI(fourPtEtaPOI);
  fQCA->Set6pPtEtaPOI(sixPtEtaPOI);
  fQCA->Set8pPtEtaPOI(eightPtEtaPOI);
  
  // non-weighted final results for differential flow for POIs:
  // 3D (pt,eta)
  fQCA->Setvn2ndPtEtaPOI(vn2ndPtEtaPOI);   
  fQCA->Setvn4thPtEtaPOI(vn4thPtEtaPOI);  
  fQCA->Setvn6thPtEtaPOI(vn6thPtEtaPOI);  
  fQCA->Setvn8thPtEtaPOI(vn8thPtEtaPOI);   
  // 2D (pt)
  fQCA->Setvn2ndPtPOI(vn2ndPtPOI);   
  fQCA->Setvn4thPtPOI(vn4thPtPOI);  
  fQCA->Setvn6thPtPOI(vn6thPtPOI);  
  fQCA->Setvn8thPtPOI(vn8thPtPOI);   
  // 2D (eta)
  fQCA->Setvn2ndEtaPOI(vn2ndEtaPOI);   
  fQCA->Setvn4thEtaPOI(vn4thEtaPOI);  
  fQCA->Setvn6thEtaPOI(vn6thEtaPOI);  
  fQCA->Setvn8thEtaPOI(vn8thEtaPOI);   
  
  // weighted correlations for each (pt,eta) bin for POIs:
  fQCA->Set2pPtEtaPOIW(twoPtEtaPOIW);
  fQCA->Set4pPtEtaPOIW(fourPtEtaPOIW);
  fQCA->Set6pPtEtaPOIW(sixPtEtaPOIW);
  fQCA->Set8pPtEtaPOIW(eightPtEtaPOIW);
  
  // weighted final results for differential flow for POIs:
  // 3D (pt,eta)
  fQCA->Setvn2ndPtEtaPOIW(vn2ndPtEtaPOIW);   
  fQCA->Setvn4thPtEtaPOIW(vn4thPtEtaPOIW);  
  fQCA->Setvn6thPtEtaPOIW(vn6thPtEtaPOIW);  
  fQCA->Setvn8thPtEtaPOIW(vn8thPtEtaPOIW); 
  // 2D (pt)
  fQCA->Setvn2ndPtPOIW(vn2ndPtPOIW);   
  fQCA->Setvn4thPtPOIW(vn4thPtPOIW);  
  fQCA->Setvn6thPtPOIW(vn6thPtPOIW);  
  fQCA->Setvn8thPtPOIW(vn8thPtPOIW);   
  // 2D (eta)
  fQCA->Setvn2ndEtaPOIW(vn2ndEtaPOIW);   
  fQCA->Setvn4thEtaPOIW(vn4thEtaPOIW);  
  fQCA->Setvn6thEtaPOIW(vn6thEtaPOIW);  
  fQCA->Setvn8thEtaPOIW(vn8thEtaPOIW);     
  
  // non-weighted correlations for each (pt,eta) bin for RPs:
  fQCA->Set2pPtEtaRP(twoPtEtaRP);
  fQCA->Set4pPtEtaRP(fourPtEtaRP);
  fQCA->Set6pPtEtaRP(sixPtEtaRP);
  fQCA->Set8pPtEtaRP(eightPtEtaRP);
  
  // non-weighted final results for differential flow for RPs:
  // 3D (pt,eta)
  fQCA->Setvn2ndPtEtaRP(vn2ndPtEtaRP);   
  fQCA->Setvn4thPtEtaRP(vn4thPtEtaRP);  
  fQCA->Setvn6thPtEtaRP(vn6thPtEtaRP);  
  fQCA->Setvn8thPtEtaRP(vn8thPtEtaRP);  
  // 2D (pt)
  fQCA->Setvn2ndPtRP(vn2ndPtRP);   
  fQCA->Setvn4thPtRP(vn4thPtRP);  
  fQCA->Setvn6thPtRP(vn6thPtRP);  
  fQCA->Setvn8thPtRP(vn8thPtRP);   
  // 2D (eta)
  fQCA->Setvn2ndEtaRP(vn2ndEtaRP);   
  fQCA->Setvn4thEtaRP(vn4thEtaRP);  
  fQCA->Setvn6thEtaRP(vn6thEtaRP);  
  fQCA->Setvn8thEtaRP(vn8thEtaRP);      
  
  // weighted correlations for each (pt,eta) bin for RPs:
  fQCA->Set2pPtEtaRPW(twoPtEtaRPW);
  fQCA->Set4pPtEtaRPW(fourPtEtaRPW);
  fQCA->Set6pPtEtaRPW(sixPtEtaRPW);
  fQCA->Set8pPtEtaRPW(eightPtEtaRPW);
  
  // weighted final results for differential flow for RPs:
  // 3D (pt,eta)
  fQCA->Setvn2ndPtEtaRPW(vn2ndPtEtaRPW);   
  fQCA->Setvn4thPtEtaRPW(vn4thPtEtaRPW);  
  fQCA->Setvn6thPtEtaRPW(vn6thPtEtaRPW);  
  fQCA->Setvn8thPtEtaRPW(vn8thPtEtaRPW);  
  // 2D (pt)
  fQCA->Setvn2ndPtRPW(vn2ndPtRPW);   
  fQCA->Setvn4thPtRPW(vn4thPtRPW);  
  fQCA->Setvn6thPtRPW(vn6thPtRPW);  
  fQCA->Setvn8thPtRPW(vn8thPtRPW);   
  // 2D (eta)
  fQCA->Setvn2ndEtaRPW(vn2ndEtaRPW);   
  fQCA->Setvn4thEtaRPW(vn4thEtaRPW);  
  fQCA->Setvn6thEtaRPW(vn6thEtaRPW);  
  fQCA->Setvn8thEtaRPW(vn8thEtaRPW);  
  
  fQCA->Finish();
  
  //----------------------------------------------------
 }
 else
 {
  cout<<" WARNING: histogram list pointer is empty (QC, Task::Terminate())"<<endl;
  cout<<endl;
 }
}





















