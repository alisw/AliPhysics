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
 * analysis task for cumulant method  * 
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
#include "TProfile.h"
#include "TProfile2D.h"
#include "TProfile3D.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisDataSlot.h"
#include "AliAnalysisDataContainer.h"
#include "AliAnalysisManager.h"

#include "AliFlowEventSimple.h"
#include "AliAnalysisTaskCumulants.h"
#include "AliFlowAnalysisWithCumulants.h"
#include "AliFlowCumuConstants.h"
#include "AliFlowCommonConstants.h"
#include "AliFlowCommonHist.h"
#include "AliFlowCommonHistResults.h"
#include "AliCumulantsFunctions.h"

ClassImp(AliAnalysisTaskCumulants)

//================================================================================================================

AliAnalysisTaskCumulants::AliAnalysisTaskCumulants(const char *name, Bool_t useWeights): 
 AliAnalysisTask(name,""), 
 fEvent(NULL),
 fGFCA(NULL), // Generating Function Cumulant (GFCA) analysis object
 fListHistos(NULL),
 fUseWeights(useWeights),
 fUsePhiWeights(kFALSE),
 fUsePtWeights(kFALSE),
 fUseEtaWeights(kFALSE),
 fListWeights(NULL)
{
 // constructor
 cout<<"AliAnalysisTaskCumulants::AliAnalysisTaskCumulants(const char *name)"<<endl;
 
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

AliAnalysisTaskCumulants::AliAnalysisTaskCumulants():
 fEvent(NULL),
 fGFCA(NULL), // Generating Function Cumulant (GFCA) analysis object
 fListHistos(NULL),
 fUseWeights(kFALSE),
 fUsePhiWeights(kFALSE),
 fUsePtWeights(kFALSE),
 fUseEtaWeights(kFALSE),
 fListWeights(NULL)
{
 // dummy constructor
 cout<<"AliAnalysisTaskCumulants::AliAnalysisTaskCumulants()"<<endl;
}

//================================================================================================================

void AliAnalysisTaskCumulants::ConnectInputData(Option_t *) 
{
 // connect ESD or AOD (called once)
 cout<<"AliAnalysisTaskCumulants::ConnectInputData(Option_t *)"<<endl;
}

//================================================================================================================

void AliAnalysisTaskCumulants::CreateOutputObjects() 
{
 // called at every worker node to initialize
 cout<<"AliAnalysisTaskCumulants::CreateOutputObjects()"<<endl;

 // analyser
 fGFCA = new AliFlowAnalysisWithCumulants();
 fGFCA->Init();
 
 //weights:
 if(fUseWeights)
 {
  //pass the flags to class:
  if(fUsePhiWeights) fGFCA->SetUsePhiWeights(fUsePhiWeights);
  if(fUsePtWeights) fGFCA->SetUsePtWeights(fUsePtWeights);
  if(fUseEtaWeights) fGFCA->SetUseEtaWeights(fUseEtaWeights);
  //get data from input slot #1 which is used for weights:
  if(GetNinputs()==2) 
  {                   
   fListWeights = (TList*)GetInputData(1); 
  }
  //pass the list with weights to class:
  if(fListWeights) fGFCA->SetWeightsList(fListWeights);
 }

 if(fGFCA->GetHistList()) 
 {
  fListHistos = fGFCA->GetHistList();
  //fListHistos->Print();
 }
 else
 {
  Printf(" ERROR: Could not retrieve histogram list (GFCA, Task::COO)"); 
 }
}

//================================================================================================================

void AliAnalysisTaskCumulants::Exec(Option_t *) 
{
 // main loop (called for each event)
 fEvent = dynamic_cast<AliFlowEventSimple*>(GetInputData(0));

 // generating function cumulants
 if(fEvent) 
 {
  fGFCA->Make(fEvent);
 }else 
  {
   cout<<" WARNING: No input data (GFCA, Task::E) !!!"<<endl;
   cout<<endl;
  }
  
 PostData(0,fListHistos);
}

//================================================================================================================

void AliAnalysisTaskCumulants::Terminate(Option_t *) 
{  
 //accessing the output list which contains the merged 2D and 3D profiles from all worker nodes
 fListHistos = (TList*)GetOutputData(0);
 //fListHistos->Print();
 
 if(fListHistos)
 {
  //histograms to store the final results
  TH1D *intFlowResults   = dynamic_cast<TH1D*>(fListHistos->FindObject("fIntFlowResultsGFC"));
  TH1D *diffFlowResults2 = dynamic_cast<TH1D*>(fListHistos->FindObject("fDiffFlowResults2ndOrderGFC"));
  TH1D *diffFlowResults4 = dynamic_cast<TH1D*>(fListHistos->FindObject("fDiffFlowResults4thOrderGFC"));
  TH1D *diffFlowResults6 = dynamic_cast<TH1D*>(fListHistos->FindObject("fDiffFlowResults6thOrderGFC"));
  TH1D *diffFlowResults8 = dynamic_cast<TH1D*>(fListHistos->FindObject("fDiffFlowResults8thOrderGFC"));
 	    	    
  //common histograms to store the final results  the integrated and differential flow
  AliFlowCommonHistResults *commonHistRes2nd = dynamic_cast<AliFlowCommonHistResults*>(fListHistos->FindObject("AliFlowCommonHistResults2ndOrderGFC"));
  AliFlowCommonHistResults *commonHistRes4th = dynamic_cast<AliFlowCommonHistResults*>(fListHistos->FindObject("AliFlowCommonHistResults4thOrderGFC"));
  AliFlowCommonHistResults *commonHistRes6th = dynamic_cast<AliFlowCommonHistResults*>(fListHistos->FindObject("AliFlowCommonHistResults6thOrderGFC"));
  AliFlowCommonHistResults *commonHistRes8th = dynamic_cast<AliFlowCommonHistResults*>(fListHistos->FindObject("AliFlowCommonHistResults8thOrderGFC"));
  
  //common control histogram
  AliFlowCommonHist *commonHists = dynamic_cast<AliFlowCommonHist*>(fListHistos->FindObject("AliFlowCommonHistGFC"));
  
  //profiles with average values of generating functions for int. and diff. flow
  TProfile2D *intFlowGenFun    = dynamic_cast<TProfile2D*>(fListHistos->FindObject("fIntFlowGenFun")); 
  
  TProfile2D *intFlowGenFun4   = dynamic_cast<TProfile2D*>(fListHistos->FindObject("fIntFlowGenFun4"));  //only for other system of Eq.
  TProfile2D *intFlowGenFun6   = dynamic_cast<TProfile2D*>(fListHistos->FindObject("fIntFlowGenFun6"));  //only for other system of Eq. 
  TProfile2D *intFlowGenFun8   = dynamic_cast<TProfile2D*>(fListHistos->FindObject("fIntFlowGenFun8"));  //only for other system of Eq.
  TProfile2D *intFlowGenFun16  = dynamic_cast<TProfile2D*>(fListHistos->FindObject("fIntFlowGenFun16")); //only for other system of Eq.  
  
  //RP, Pt:
  TProfile3D *diffFlowPtRPGenFunRe = dynamic_cast<TProfile3D*>(fListHistos->FindObject("fDiffFlowPtRPGenFunRe"));
  TProfile3D *diffFlowPtRPGenFunIm = dynamic_cast<TProfile3D*>(fListHistos->FindObject("fDiffFlowPtRPGenFunIm"));
  TProfile *ptBinRPNoOfParticles = dynamic_cast<TProfile*>(fListHistos->FindObject("fPtBinRPNoOfParticles"));
  
  //RP, Eta:
  TProfile3D *diffFlowEtaRPGenFunRe = dynamic_cast<TProfile3D*>(fListHistos->FindObject("fDiffFlowEtaRPGenFunRe"));
  TProfile3D *diffFlowEtaRPGenFunIm = dynamic_cast<TProfile3D*>(fListHistos->FindObject("fDiffFlowEtaRPGenFunIm"));
  TProfile *etaBinRPNoOfParticles = dynamic_cast<TProfile*>(fListHistos->FindObject("fEtaBinRPNoOfParticles"));
  
  //POI, Pt:
  TProfile3D *diffFlowPtPOIGenFunRe = dynamic_cast<TProfile3D*>(fListHistos->FindObject("fDiffFlowPtPOIGenFunRe"));
  TProfile3D *diffFlowPtPOIGenFunIm = dynamic_cast<TProfile3D*>(fListHistos->FindObject("fDiffFlowPtPOIGenFunIm"));
  TProfile *ptBinPOINoOfParticles = dynamic_cast<TProfile*>(fListHistos->FindObject("fPtBinPOINoOfParticles"));
  
  //POI, Eta:
  TProfile3D *diffFlowEtaPOIGenFunRe = dynamic_cast<TProfile3D*>(fListHistos->FindObject("fDiffFlowEtaPOIGenFunRe"));
  TProfile3D *diffFlowEtaPOIGenFunIm = dynamic_cast<TProfile3D*>(fListHistos->FindObject("fDiffFlowEtaPOIGenFunIm"));
  TProfile *etaBinPOINoOfParticles = dynamic_cast<TProfile*>(fListHistos->FindObject("fEtaBinPOINoOfParticles"));
  
  //average selected multiplicity (for int. flow) 
  TProfile *avMult = dynamic_cast<TProfile*>(fListHistos->FindObject("fAvMultIntFlowGFC"));
  
  TProfile *avMult4  = dynamic_cast<TProfile*>(fListHistos->FindObject("fAvMultIntFlow4GFCA"));  //only for other system of Eq.
  TProfile *avMult6  = dynamic_cast<TProfile*>(fListHistos->FindObject("fAvMultIntFlow6GFCA"));  //only for other system of Eq.
  TProfile *avMult8  = dynamic_cast<TProfile*>(fListHistos->FindObject("fAvMultIntFlow8GFCA"));  //only for other system of Eq.
  TProfile *avMult16 = dynamic_cast<TProfile*>(fListHistos->FindObject("fAvMultIntFlow16GFCA")); //only for other system of Eq.
  
  //average values of Q-vector components (1st bin: <Q_x>, 2nd bin: <Q_y>, 3rd bin: <(Q_x)^2>, 4th bin: <(Q_y)^2>) 
  TProfile *qVectorComponents = dynamic_cast<TProfile*>(fListHistos->FindObject("fQVectorComponentsGFC"));
  
  //<w^2> 
  TProfile *averageOfSquaredWeight = dynamic_cast<TProfile*>(fListHistos->FindObject("fAverageOfSquaredWeight"));
      
  /*
  TProfile2D *diffFlowPtGenFunRe0 = dynamic_cast<TProfile2D*>(fListHistos->FindObject("fdiffFlowPtGenFunRe0"));
  TProfile2D *diffFlowPtGenFunRe1 = dynamic_cast<TProfile2D*>(fListHistos->FindObject("fdiffFlowPtGenFunRe1")); 
  TProfile2D *diffFlowPtGenFunRe2 = dynamic_cast<TProfile2D*>(fListHistos->FindObject("fdiffFlowPtGenFunRe2")); 
  TProfile2D *diffFlowPtGenFunRe3 = dynamic_cast<TProfile2D*>(fListHistos->FindObject("fdiffFlowPtGenFunRe3")); 
  TProfile2D *diffFlowPtGenFunRe4 = dynamic_cast<TProfile2D*>(fListHistos->FindObject("fdiffFlowPtGenFunRe4")); 
  TProfile2D *diffFlowPtGenFunRe5 = dynamic_cast<TProfile2D*>(fListHistos->FindObject("fdiffFlowPtGenFunRe5")); 
  TProfile2D *diffFlowPtGenFunRe6 = dynamic_cast<TProfile2D*>(fListHistos->FindObject("fdiffFlowPtGenFunRe6")); 
  TProfile2D *diffFlowPtGenFunRe7 = dynamic_cast<TProfile2D*>(fListHistos->FindObject("fdiffFlowPtGenFunRe7")); 
  TProfile2D *diffFlowPtGenFunIm0 = dynamic_cast<TProfile2D*>(fListHistos->FindObject("fdiffFlowPtGenFunIm0")); 
  TProfile2D *diffFlowPtGenFunIm1 = dynamic_cast<TProfile2D*>(fListHistos->FindObject("fdiffFlowPtGenFunIm1")); 
  TProfile2D *diffFlowPtGenFunIm2 = dynamic_cast<TProfile2D*>(fListHistos->FindObject("fdiffFlowPtGenFunIm2")); 
  TProfile2D *diffFlowPtGenFunIm3 = dynamic_cast<TProfile2D*>(fListHistos->FindObject("fdiffFlowPtGenFunIm3")); 
  TProfile2D *diffFlowPtGenFunIm4 = dynamic_cast<TProfile2D*>(fListHistos->FindObject("fdiffFlowPtGenFunIm4")); 
  TProfile2D *diffFlowPtGenFunIm5 = dynamic_cast<TProfile2D*>(fListHistos->FindObject("fdiffFlowPtGenFunIm5")); 
  TProfile2D *diffFlowPtGenFunIm6 = dynamic_cast<TProfile2D*>(fListHistos->FindObject("fdiffFlowPtGenFunIm6")); 
  TProfile2D *diffFlowPtGenFunIm7 = dynamic_cast<TProfile2D*>(fListHistos->FindObject("fdiffFlowPtGenFunIm7")); 
  */

  //profile with avarage selected multiplicity for int. flow 
  //TProfile *avMult = dynamic_cast<TProfile*>(fListHistos->FindObject("fAvMultIntFlow"));
  
  //profile with avarage values of Q-vector components (1st bin: <Q_x>, 2nd bin: <Q_y>, 3rd bin: <(Q_x)^2>, 4th bin: <(Q_y)^2>) 
  //TProfile *QVectorComponents = dynamic_cast<TProfile*>(fListHistos->FindObject("fQVectorComponents"));
  
  //q-distribution
  //TH1D *qDist = dynamic_cast<TH1D*>(fListHistos->FindObject("fQDist"));
  
  //AliCumulantsFunctions finalResults(intFlowGenFun,NULL,NULL, intFlowResults,diffFlowResults2,diffFlowResults4,diffFlowResults6,diffFlowResults8,avMult,QVectorComponents,qDist,diffFlowPtGenFunRe0,diffFlowPtGenFunRe1,diffFlowPtGenFunRe2, diffFlowPtGenFunRe3,diffFlowPtGenFunRe4,diffFlowPtGenFunRe5,diffFlowPtGenFunRe6,diffFlowPtGenFunRe7,diffFlowPtGenFunIm0,diffFlowPtGenFunIm1, diffFlowPtGenFunIm2,diffFlowPtGenFunIm3,diffFlowPtGenFunIm4,diffFlowPtGenFunIm5,diffFlowPtGenFunIm6,diffFlowPtGenFunIm7);
  
  //AliCumulantsFunctions finalResults(intFlowGenFun,diffFlowPtGenFunRe,diffFlowPtGenFunIm, intFlowResults,diffFlowResults2,diffFlowResults4,diffFlowResults6,diffFlowResults8,avMult,QVectorComponents,qDist);
         
  //finalResults.Calculate();  
  
  
  
  //----------------------------------------------------
 
  fGFCA = new AliFlowAnalysisWithCumulants();  
 
  fGFCA->SetIntFlowResults(intFlowResults); 
  fGFCA->SetDiffFlowResults2nd(diffFlowResults2);
  fGFCA->SetDiffFlowResults4th(diffFlowResults4);
  fGFCA->SetDiffFlowResults6th(diffFlowResults6);
  fGFCA->SetDiffFlowResults8th(diffFlowResults8); 
  
  fGFCA->SetCommonHistsResults2nd(commonHistRes2nd); 
  fGFCA->SetCommonHistsResults4th(commonHistRes4th);
  fGFCA->SetCommonHistsResults6th(commonHistRes6th);
  fGFCA->SetCommonHistsResults8th(commonHistRes8th);
  
  fGFCA->SetCommonHists(commonHists);
  
  fGFCA->SetIntFlowGenFun(intFlowGenFun);
  
  fGFCA->SetIntFlowGenFun4(intFlowGenFun4);   //only for other system of Eq.
  fGFCA->SetIntFlowGenFun6(intFlowGenFun6);   //only for other system of Eq.
  fGFCA->SetIntFlowGenFun8(intFlowGenFun8);   //only for other system of Eq.
  fGFCA->SetIntFlowGenFun16(intFlowGenFun16); //only for other system of Eq. 
  
  fGFCA->SetDiffFlowPtRPGenFunRe(diffFlowPtRPGenFunRe);
  fGFCA->SetDiffFlowPtRPGenFunIm(diffFlowPtRPGenFunIm);
  fGFCA->SetNumberOfParticlesPerPtBinRP(ptBinRPNoOfParticles);
  
  fGFCA->SetDiffFlowEtaRPGenFunRe(diffFlowEtaRPGenFunRe);
  fGFCA->SetDiffFlowEtaRPGenFunIm(diffFlowEtaRPGenFunIm);
  fGFCA->SetNumberOfParticlesPerEtaBinRP(etaBinRPNoOfParticles);     
  
  fGFCA->SetDiffFlowPtPOIGenFunRe(diffFlowPtPOIGenFunRe);
  fGFCA->SetDiffFlowPtPOIGenFunIm(diffFlowPtPOIGenFunIm);
  fGFCA->SetNumberOfParticlesPerPtBinPOI(ptBinPOINoOfParticles);
  
  fGFCA->SetDiffFlowEtaPOIGenFunRe(diffFlowEtaPOIGenFunRe);
  fGFCA->SetDiffFlowEtaPOIGenFunIm(diffFlowEtaPOIGenFunIm);
  fGFCA->SetNumberOfParticlesPerEtaBinPOI(etaBinPOINoOfParticles);
  
  fGFCA->SetAverageMultiplicity(avMult);
  
  fGFCA->SetAverageMultiplicity4(avMult4);   //only for other system of Eq.
  fGFCA->SetAverageMultiplicity6(avMult6);   //only for other system of Eq.
  fGFCA->SetAverageMultiplicity8(avMult8);   //only for other system of Eq.
  fGFCA->SetAverageMultiplicity16(avMult16); //only for other system of Eq.
  
  fGFCA->SetQVectorComponents(qVectorComponents);
  
  fGFCA->SetAverageOfSquaredWeight(averageOfSquaredWeight);
  
  fGFCA->Finish();
  
  //----------------------------------------------------
 }
 else
 {
  cout<<" WARNING: histogram list pointer is empty (GFC, Task::T)"<<endl;
  cout<<endl;
 }
}





















