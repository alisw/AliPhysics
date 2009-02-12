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

#include "AliESDEvent.h"
#include "AliESDInputHandler.h"

#include "AliAODEvent.h"
#include "AliAODInputHandler.h"

#include "AliMCEventHandler.h"
#include "AliMCEvent.h"

#include "AliCFManager.h"

#include "AliAnalysisTaskCumulants.h"
#include "AliFlowEventSimpleMaker.h"
#include "AliFlowAnalysisWithCumulants.h"
#include "AliFlowCumuConstants.h"
#include "AliFlowCommonConstants.h"
#include "AliFlowCommonHist.h"
#include "AliFlowCommonHistResults.h"
#include "AliCumulantsFunctions.h"

ClassImp(AliAnalysisTaskCumulants)

//================================================================================================================

AliAnalysisTaskCumulants::AliAnalysisTaskCumulants(const char *name, Bool_t on): 
 AliAnalysisTask(name,""), 
 fESD(NULL),
 fAOD(NULL),
 fGFC(NULL),//Generating Function Cumulant (GFC) analysis object
 fEventMaker(NULL),
 fAnalysisType("ESD"), 
 fCFManager1(NULL),
 fCFManager2(NULL),
 fListHistos(NULL),
 fQAInt(NULL),
 fQADiff(NULL),
 fQA(on)
{
//constructor
 cout<<"AliAnalysisTaskCumulants::AliAnalysisTaskCumulants(const char *name)"<<endl;
 
 // Define input and output slots here
 // Input slot #0 works with a TChain
 DefineInput(0, TChain::Class());
  
 // Output slot #0 writes into a TList container
 DefineOutput(0, TList::Class());  
 if(on) 
 {
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class()); 
 }  
}

AliAnalysisTaskCumulants::AliAnalysisTaskCumulants(): 
 fESD(NULL),
 fAOD(NULL), 
 fGFC(NULL),//Generating Function Cumulant (GFC) analysis object
 fEventMaker(NULL),
 fAnalysisType("ESD"),
 fCFManager1(NULL),
 fCFManager2(NULL),
 fListHistos(NULL),  
 fQAInt(NULL),
 fQADiff(NULL),
 fQA(kFALSE)
{
 //dummy constructor
 cout<<"AliAnalysisTaskCumulants::AliAnalysisTaskCumulants()"<<endl;
}

//================================================================================================================

void AliAnalysisTaskCumulants::ConnectInputData(Option_t *) 
{
 //connect ESD or AOD (called once)
 cout<<"AliAnalysisTaskCumulants::ConnectInputData(Option_t *)"<<endl;

 TTree* tree = dynamic_cast<TTree*> (GetInputData(0));
 if (!tree) 
 {
  Printf("ERROR: Could not read chain from input slot 0");
 } 
 else 
 {
 //disable all branches and enable only the needed ones
  if (fAnalysisType == "MC") {
     // we want to process only MC
      tree->SetBranchStatus("*", kFALSE);

      AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());

      if (!esdH) {
	Printf("ERROR: Could not get ESDInputHandler");
      } else {
	fESD = esdH->GetEvent();
      }
    }
    else if (fAnalysisType == "ESD" || fAnalysisType == "ESDMC0" || fAnalysisType == "ESDMC1"  ) {
      tree->SetBranchStatus("*", kFALSE);
      tree->SetBranchStatus("Tracks.*", kTRUE);

      AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());

      if (!esdH) {
	Printf("ERROR: Could not get ESDInputHandler");
      } else
	fESD = esdH->GetEvent();
    }
    else if (fAnalysisType == "AOD") {
      AliAODInputHandler *aodH = dynamic_cast<AliAODInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());

      if (!aodH) {
	Printf("ERROR: Could not get AODInputHandler");
      }
      else {
	fAOD = aodH->GetEvent();
      }
    }
    else {
      Printf("Wrong analysis type: Only ESD, ESDMC0, ESDMC1, AOD and MC types are allowed!");

    }
  }
}

//================================================================================================================

void AliAnalysisTaskCumulants::CreateOutputObjects() 
{
 //called at every worker node to initialize
 cout<<"AliAnalysisTaskCumulants::CreateOutputObjects()"<<endl;

 
 //OpenFile(0);
 

 if(!(fAnalysisType == "AOD" || fAnalysisType == "ESD" || fAnalysisType == "ESDMC0" || fAnalysisType == "ESDMC1" || fAnalysisType == "MC")) 
 {
  cout<<"WRONG ANALYSIS TYPE! only ESD, ESDMC0, ESDMC1, AOD and MC are allowed."<<endl;
  exit(1);
 }
 
 //event maker
 fEventMaker = new AliFlowEventSimpleMaker();
  
 //analyser
 fGFC = new AliFlowAnalysisWithCumulants();
 fGFC->CreateOutputObjects();

 if(fGFC->GetHistList()) 
 {
  fListHistos = fGFC->GetHistList();
  //fListHistos->Print();
 }
 else 
 {
  Printf("ERROR: Could not retrieve histogram list"); 
 }
}

//================================================================================================================

void AliAnalysisTaskCumulants::Exec(Option_t *) 
{
 //main loop (called for each event)
 if (fAnalysisType == "MC") {
    // Process MC truth, therefore we receive the AliAnalysisManager and ask it for the AliMCEventHandler
    // This handler can return the current MC event

    AliMCEventHandler* eventHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
    if (!eventHandler) {
      Printf("ERROR: Could not retrieve MC event handler");
      return;
    }

    AliMCEvent* mcEvent = eventHandler->MCEvent();
    if (!mcEvent) {
      Printf("ERROR: Could not retrieve MC event");
      return;
    }


    fCFManager1->SetEventInfo(mcEvent);
    fCFManager2->SetEventInfo(mcEvent);

    Printf("MC particles: %d", mcEvent->GetNumberOfTracks());

    //cumulant analysis 
    AliFlowEventSimple* fEvent = fEventMaker->FillTracks(mcEvent,fCFManager1,fCFManager2);
    fGFC->Make(fEvent);
    delete fEvent;
  }
  else if (fAnalysisType == "ESD") {
    if (!fESD) {
      Printf("ERROR: fESD not available");
      return;
    }
    Printf("There are %d tracks in this event", fESD->GetNumberOfTracks());
    
    //cumulant analysis 
    AliFlowEventSimple* fEvent = fEventMaker->FillTracks(fESD,fCFManager1,fCFManager2);//cuts
    //AliFlowEventSimple* fEvent = fEventMaker->FillTracks(fESD);
    
    fGFC->Make(fEvent);
    delete fEvent;
  }
  else if (fAnalysisType == "ESDMC0") {
    if (!fESD) {
      Printf("ERROR: fESD not available");
      return;
    }
    Printf("There are %d tracks in this event", fESD->GetNumberOfTracks());
    
    AliMCEventHandler* eventHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
    if (!eventHandler) {
      Printf("ERROR: Could not retrieve MC event handler");
      return;
    }

    AliMCEvent* mcEvent = eventHandler->MCEvent();
    if (!mcEvent) {
      Printf("ERROR: Could not retrieve MC event");
      return;
    }

    fCFManager1->SetEventInfo(mcEvent);
    fCFManager2->SetEventInfo(mcEvent);

    //cumulant analysis 
    AliFlowEventSimple* fEvent=NULL;
    if (fAnalysisType == "ESDMC0") { 
      fEvent = fEventMaker->FillTracks(fESD, mcEvent, fCFManager1, fCFManager2, 0); //0 = kine from ESD, 1 = kine from MC
    } else if (fAnalysisType == "ESDMC1") {
      fEvent = fEventMaker->FillTracks(fESD, mcEvent, fCFManager1, fCFManager2, 1); //0 = kine from ESD, 1 = kine from MC
    }
    fGFC->Make(fEvent);
    delete fEvent;
    //delete mcEvent;
  }
  
  else if (fAnalysisType == "AOD") {
    if (!fAOD) {
      Printf("ERROR: fAOD not available");
      return;
    }
    Printf("There are %d tracks in this event", fAOD->GetNumberOfTracks());

    // analysis 
    //For the moment don't use CF //AliFlowEventSimple* fEvent = fEventMaker->FillTracks(fAOD,fCFManager1,fCFManager2);
    AliFlowEventSimple* fEvent = fEventMaker->FillTracks(fAOD);
    fGFC->Make(fEvent);
    delete fEvent;
  }


  PostData(0,fListHistos); 
  if(fQA) 
  {
   PostData(1,fQAInt);
   PostData(2,fQADiff); 
  }
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
  
  TProfile *avMult4  = dynamic_cast<TProfile*>(fListHistos->FindObject("fAvMultIntFlow4GFC"));  //only for other system of Eq.
  TProfile *avMult6  = dynamic_cast<TProfile*>(fListHistos->FindObject("fAvMultIntFlow6GFC"));  //only for other system of Eq.
  TProfile *avMult8  = dynamic_cast<TProfile*>(fListHistos->FindObject("fAvMultIntFlow8GFC"));  //only for other system of Eq.
  TProfile *avMult16 = dynamic_cast<TProfile*>(fListHistos->FindObject("fAvMultIntFlow16GFC")); //only for other system of Eq.
  
  //average values of Q-vector components (1st bin: <Q_x>, 2nd bin: <Q_y>, 3rd bin: <(Q_x)^2>, 4th bin: <(Q_y)^2>) 
  TProfile *qVectorComponents = dynamic_cast<TProfile*>(fListHistos->FindObject("fQVectorComponentsGFC"));
      
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
 
  fGFC = new AliFlowAnalysisWithCumulants();  
 
  fGFC->SetIntFlowResults(intFlowResults); 
  fGFC->SetDiffFlowResults2nd(diffFlowResults2);
  fGFC->SetDiffFlowResults4th(diffFlowResults4);
  fGFC->SetDiffFlowResults6th(diffFlowResults6);
  fGFC->SetDiffFlowResults8th(diffFlowResults8); 
  
  fGFC->SetCommonHistsResults2nd(commonHistRes2nd); 
  fGFC->SetCommonHistsResults4th(commonHistRes4th);
  fGFC->SetCommonHistsResults6th(commonHistRes6th);
  fGFC->SetCommonHistsResults8th(commonHistRes8th);
  
  fGFC->SetCommonHists(commonHists);
  
  fGFC->SetIntFlowGenFun(intFlowGenFun);
  
  fGFC->SetIntFlowGenFun4(intFlowGenFun4);   //only for other system of Eq.
  fGFC->SetIntFlowGenFun6(intFlowGenFun6);   //only for other system of Eq.
  fGFC->SetIntFlowGenFun8(intFlowGenFun8);   //only for other system of Eq.
  fGFC->SetIntFlowGenFun16(intFlowGenFun16); //only for other system of Eq. 
  
  fGFC->SetDiffFlowPtRPGenFunRe(diffFlowPtRPGenFunRe);
  fGFC->SetDiffFlowPtRPGenFunIm(diffFlowPtRPGenFunIm);
  fGFC->SetNumberOfParticlesPerPtBinRP(ptBinRPNoOfParticles);
  
  fGFC->SetDiffFlowEtaRPGenFunRe(diffFlowEtaRPGenFunRe);
  fGFC->SetDiffFlowEtaRPGenFunIm(diffFlowEtaRPGenFunIm);
  fGFC->SetNumberOfParticlesPerEtaBinRP(etaBinRPNoOfParticles);     
  
  fGFC->SetDiffFlowPtPOIGenFunRe(diffFlowPtPOIGenFunRe);
  fGFC->SetDiffFlowPtPOIGenFunIm(diffFlowPtPOIGenFunIm);
  fGFC->SetNumberOfParticlesPerPtBinPOI(ptBinPOINoOfParticles);
  
  fGFC->SetDiffFlowEtaPOIGenFunRe(diffFlowEtaPOIGenFunRe);
  fGFC->SetDiffFlowEtaPOIGenFunIm(diffFlowEtaPOIGenFunIm);
  fGFC->SetNumberOfParticlesPerEtaBinPOI(etaBinPOINoOfParticles);
  
  fGFC->SetAverageMultiplicity(avMult);
  
  fGFC->SetAverageMultiplicity4(avMult4);   //only for other system of Eq.
  fGFC->SetAverageMultiplicity6(avMult6);   //only for other system of Eq.
  fGFC->SetAverageMultiplicity8(avMult8);   //only for other system of Eq.
  fGFC->SetAverageMultiplicity16(avMult16); //only for other system of Eq.
  
  fGFC->SetQVectorComponents(qVectorComponents);
  
  fGFC->Finish();
  
  //----------------------------------------------------
 }
 else
 {
  cout<<"histogram list pointer is empty"<<endl;
 }
}





















