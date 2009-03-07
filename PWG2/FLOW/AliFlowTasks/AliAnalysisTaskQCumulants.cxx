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
#include "AliQCumulantsFunctions.h"

ClassImp(AliAnalysisTaskQCumulants)

//================================================================================================================

AliAnalysisTaskQCumulants::AliAnalysisTaskQCumulants(const char *name, Bool_t useWeights): 
 AliAnalysisTask(name,""), 
 fEvent(NULL),
 fQCA(NULL), // Q-cumulant Analysis (QCA) object
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
  //final results (integrated flow)
  TH1D *intFlowResults = dynamic_cast<TH1D*>(fListHistos->FindObject("fIntFlowResultsQC"));
  
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
  TProfile *QCorrelations = dynamic_cast<TProfile*>(fListHistos->FindObject("fQCorrelations"));
  
  //weighted multi-particle correlations calculated from Q-vectors
  TProfile *QCorrelationsW = dynamic_cast<TProfile*>(fListHistos->FindObject("fQCorrelationsW"));
  
  //average of products: 1st bin: <2*4>, 2nd bin: <2*6>, ...
  TProfile *QProduct = dynamic_cast<TProfile*>(fListHistos->FindObject("fQProduct"));
  
  //average 2-, 3- and 4-particle correlations per pt-bin 
  TProfile *binnedPt2p1n1nRP = dynamic_cast<TProfile*>((dynamic_cast<TList*>(fListHistos->FindObject("DifferentialFlow")))->FindObject("f2PerPtBin1n1nRP"));
  TProfile *binnedPt2p2n2nRP = dynamic_cast<TProfile*>((dynamic_cast<TList*>(fListHistos->FindObject("DifferentialFlow")))->FindObject("f2PerPtBin2n2nRP"));
  TProfile *binnedPt3p2n1n1nRP = dynamic_cast<TProfile*>((dynamic_cast<TList*>(fListHistos->FindObject("DifferentialFlow")))->FindObject("f3PerPtBin2n1n1nRP"));
  TProfile *binnedPt3p1n1n2nRP = dynamic_cast<TProfile*>((dynamic_cast<TList*>(fListHistos->FindObject("DifferentialFlow")))->FindObject("f3PerPtBin1n1n2nRP"));
  TProfile *binnedPt4p1n1n1n1nRP = dynamic_cast<TProfile*>((dynamic_cast<TList*>(fListHistos->FindObject("DifferentialFlow")))->FindObject("f4PerPtBin1n1n1n1nRP"));
  
  //average 2-, 3- and 4-particle correlations per eta-bin 
  TProfile *binnedEta2p1n1nRP = dynamic_cast<TProfile*>((dynamic_cast<TList*>(fListHistos->FindObject("DifferentialFlow")))->FindObject("f2PerEtaBin1n1nRP"));
  TProfile *binnedEta2p2n2nRP = dynamic_cast<TProfile*>((dynamic_cast<TList*>(fListHistos->FindObject("DifferentialFlow")))->FindObject("f2PerEtaBin2n2nRP"));
  TProfile *binnedEta3p2n1n1nRP =  dynamic_cast<TProfile*>((dynamic_cast<TList*>(fListHistos->FindObject("DifferentialFlow")))->FindObject("f3PerEtaBin2n1n1nRP"));
  TProfile *binnedEta3p1n1n2nRP = dynamic_cast<TProfile*>((dynamic_cast<TList*>(fListHistos->FindObject("DifferentialFlow")))->FindObject("f3PerEtaBin1n1n2nRP"));
  TProfile *binnedEta4p1n1n1n1nRP = dynamic_cast<TProfile*>((dynamic_cast<TList*>(fListHistos->FindObject("DifferentialFlow")))->FindObject("f4PerEtaBin1n1n1n1nRP"));  
  
  //average 2-, 3- and 4-particle correlations per pt-bin 
  TProfile *binnedPt2p1n1nPOI = dynamic_cast<TProfile*>((dynamic_cast<TList*>(fListHistos->FindObject("DifferentialFlow")))->FindObject("f2PerPtBin1n1nPOI"));
  TProfile *binnedPt2p2n2nPOI = dynamic_cast<TProfile*>((dynamic_cast<TList*>(fListHistos->FindObject("DifferentialFlow")))->FindObject("f2PerPtBin2n2nPOI"));
  TProfile *binnedPt3p2n1n1nPOI = dynamic_cast<TProfile*>((dynamic_cast<TList*>(fListHistos->FindObject("DifferentialFlow")))->FindObject("f3PerPtBin2n1n1nPOI"));
  TProfile *binnedPt3p1n1n2nPOI = dynamic_cast<TProfile*>((dynamic_cast<TList*>(fListHistos->FindObject("DifferentialFlow")))->FindObject("f3PerPtBin1n1n2nPOI"));
  TProfile *binnedPt4p1n1n1n1nPOI = dynamic_cast<TProfile*>((dynamic_cast<TList*>(fListHistos->FindObject("DifferentialFlow")))->FindObject("f4PerPtBin1n1n1n1nPOI"));
  
  //average 2-, 3- and 4-particle correlations per eta-bin 
  TProfile *binnedEta2p1n1nPOI = dynamic_cast<TProfile*>((dynamic_cast<TList*>(fListHistos->FindObject("DifferentialFlow")))->FindObject("f2PerEtaBin1n1nPOI"));
  TProfile *binnedEta2p2n2nPOI = dynamic_cast<TProfile*>((dynamic_cast<TList*>(fListHistos->FindObject("DifferentialFlow")))->FindObject("f2PerEtaBin2n2nPOI"));
  TProfile *binnedEta3p2n1n1nPOI =  dynamic_cast<TProfile*>((dynamic_cast<TList*>(fListHistos->FindObject("DifferentialFlow")))->FindObject("f3PerEtaBin2n1n1nPOI"));
  TProfile *binnedEta3p1n1n2nPOI = dynamic_cast<TProfile*>((dynamic_cast<TList*>(fListHistos->FindObject("DifferentialFlow")))->FindObject("f3PerEtaBin1n1n2nPOI"));
  TProfile *binnedEta4p1n1n1n1nPOI = dynamic_cast<TProfile*>((dynamic_cast<TList*>(fListHistos->FindObject("DifferentialFlow")))->FindObject("f4PerEtaBin1n1n1n1nPOI"));  
  
  
  
  //average 2-, 3- and 4-particle correlations per pt-bin 
  TProfile *binnedWPt2p1n1nPOI = dynamic_cast<TProfile*>((dynamic_cast<TList*>(fListHistos->FindObject("DifferentialFlow")))->FindObject("f2WPerPtBin1n1nPOI"));
  TProfile *binnedWPt2p2n2nPOI = dynamic_cast<TProfile*>((dynamic_cast<TList*>(fListHistos->FindObject("DifferentialFlow")))->FindObject("f2WPerPtBin2n2nPOI"));
  TProfile *binnedWPt3p2n1n1nPOI = dynamic_cast<TProfile*>((dynamic_cast<TList*>(fListHistos->FindObject("DifferentialFlow")))->FindObject("f3WPerPtBin2n1n1nPOI"));
  TProfile *binnedWPt3p1n1n2nPOI = dynamic_cast<TProfile*>((dynamic_cast<TList*>(fListHistos->FindObject("DifferentialFlow")))->FindObject("f3WPerPtBin1n1n2nPOI"));
  TProfile *binnedWPt4p1n1n1n1nPOI = dynamic_cast<TProfile*>((dynamic_cast<TList*>(fListHistos->FindObject("DifferentialFlow")))->FindObject("f4WPerPtBin1n1n1n1nPOI"));
  
  cout<<endl;
  cout<<"ab"<<endl;
  cout<<binnedWPt4p1n1n1n1nPOI->GetName()<<endl;
  cout<<endl;
  
  
  
   TProfile *binnedWEta2p1n1nPOI = dynamic_cast<TProfile*>((dynamic_cast<TList*>(fListHistos->FindObject("DifferentialFlow")))->FindObject("f2WPerEtaBin1n1nPOI"));
  
  TProfile *binnedWEta4p1n1n1n1nPOI = dynamic_cast<TProfile*>((dynamic_cast<TList*>(fListHistos->FindObject("DifferentialFlow")))->FindObject("f4WPerEtaBin1n1n1n1nPOI"));
  
  
  TProfile *binnedWPt2p1n1nRP = dynamic_cast<TProfile*>((dynamic_cast<TList*>(fListHistos->FindObject("DifferentialFlow")))->FindObject("f2WPerPtBin1n1nRP"));
  
  TProfile *binnedWPt4p1n1n1n1nRP = dynamic_cast<TProfile*>((dynamic_cast<TList*>(fListHistos->FindObject("DifferentialFlow")))->FindObject("f4WPerPtBin1n1n1n1nRP"));
  
  TProfile *binnedWEta2p1n1nRP = dynamic_cast<TProfile*>((dynamic_cast<TList*>(fListHistos->FindObject("DifferentialFlow")))->FindObject("f2WPerEtaBin1n1nRP"));
  
  TProfile *binnedWEta4p1n1n1n1nRP = dynamic_cast<TProfile*>((dynamic_cast<TList*>(fListHistos->FindObject("DifferentialFlow")))->FindObject("f4WPerEtaBin1n1n1n1nRP"));
  
  
  
  
  
  
    
  //average values of Q-vector components (1st bin: <Q_x>, 2nd bin: <Q_y>, 3rd bin: <(Q_x)^2>, 4th bin: <(Q_y)^2>) 
  TProfile *QVectorComponents = dynamic_cast<TProfile*>(fListHistos->FindObject("fQvectorComponents"));
  
  //multi-particle correlations calculated with nested loop 
  TProfile *DirectCorrelations = dynamic_cast<TProfile*>(fListHistos->FindObject("fDirectCorrelations"));
 
  //----------------------------------------------------
 
  fQCA = new AliFlowAnalysisWithQCumulants();  
 
  fQCA->SetIntFlowResults(intFlowResults); 
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
  fQCA->SetQCorrelations(QCorrelations);
  fQCA->SetQCorrelationsW(QCorrelationsW);
  fQCA->SetQProduct(QProduct);
  fQCA->SetQVectorComponents(QVectorComponents);
 
  fQCA->SetTwo1n1nPerPtBinRP(binnedPt2p1n1nRP);
  fQCA->SetTwo2n2nPerPtBinRP(binnedPt2p2n2nRP);
  fQCA->SetThree2n1n1nPerPtBinRP(binnedPt3p2n1n1nRP);
  fQCA->SetThree1n1n2nPerPtBinRP(binnedPt3p1n1n2nRP);
  fQCA->SetFour1n1n1n1nPerPtBinRP(binnedPt4p1n1n1n1nRP);
  
  fQCA->SetTwo1n1nPerEtaBinRP(binnedEta2p1n1nRP);
  fQCA->SetTwo2n2nPerEtaBinRP(binnedEta2p2n2nRP);
  fQCA->SetThree2n1n1nPerEtaBinRP(binnedEta3p2n1n1nRP);
  fQCA->SetThree1n1n2nPerEtaBinRP(binnedEta3p1n1n2nRP);
  fQCA->SetFour1n1n1n1nPerEtaBinRP(binnedEta4p1n1n1n1nRP); 
  
  fQCA->SetTwo1n1nPerPtBinPOI(binnedPt2p1n1nPOI);
  fQCA->SetTwo2n2nPerPtBinPOI(binnedPt2p2n2nPOI);
  fQCA->SetThree2n1n1nPerPtBinPOI(binnedPt3p2n1n1nPOI);
  fQCA->SetThree1n1n2nPerPtBinPOI(binnedPt3p1n1n2nPOI);
  fQCA->SetFour1n1n1n1nPerPtBinPOI(binnedPt4p1n1n1n1nPOI);
  
  fQCA->SetTwo1n1nPerEtaBinPOI(binnedEta2p1n1nPOI);
  fQCA->SetTwo2n2nPerEtaBinPOI(binnedEta2p2n2nPOI);
  fQCA->SetThree2n1n1nPerEtaBinPOI(binnedEta3p2n1n1nPOI);
  fQCA->SetThree1n1n2nPerEtaBinPOI(binnedEta3p1n1n2nPOI);
  fQCA->SetFour1n1n1n1nPerEtaBinPOI(binnedEta4p1n1n1n1nPOI);    
 
 
 
 
 
  fQCA->SetTwo1n1nWPerPtBinPOI(binnedWPt2p1n1nPOI);
  fQCA->SetTwo2n2nWPerPtBinPOI(binnedWPt2p2n2nPOI);
  fQCA->SetThree2n1n1nWPerPtBinPOI(binnedWPt3p2n1n1nPOI);
  fQCA->SetThree1n1n2nWPerPtBinPOI(binnedWPt3p1n1n2nPOI);
  fQCA->SetFour1n1n1n1nWPerPtBinPOI(binnedWPt4p1n1n1n1nPOI);
  
  fQCA->SetTwo1n1nWPerEtaBinPOI(binnedWEta2p1n1nPOI);
  fQCA->SetFour1n1n1n1nWPerEtaBinPOI(binnedWEta4p1n1n1n1nPOI);
  
  fQCA->SetTwo1n1nWPerPtBinRP(binnedWPt2p1n1nRP);
  fQCA->SetFour1n1n1n1nWPerPtBinRP(binnedWPt4p1n1n1n1nRP);
  
  fQCA->SetTwo1n1nWPerEtaBinRP(binnedWEta2p1n1nRP);
  fQCA->SetFour1n1n1n1nWPerEtaBinRP(binnedWEta4p1n1n1n1nRP);
  
  
  
  
  
  
  
  fQCA->SetDirectCorrelations(DirectCorrelations);
 
  fQCA->Finish();
  
  //----------------------------------------------------
 }
 else
 {
  cout<<" WARNING: histogram list pointer is empty (QC, Task::Terminate())"<<endl;
  cout<<endl;
 }
}





















