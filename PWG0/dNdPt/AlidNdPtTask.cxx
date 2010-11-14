/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

#include "iostream"

#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TList.h"
#include "TFile.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDVertex.h"
#include "AliTracker.h"
#include "AliGeomManager.h"

#include "AliESDVZERO.h"
#include "AliMultiplicity.h"

#include "AliESDtrackCuts.h"
#include "AliMCEventHandler.h"
#include "dNdPt/AlidNdPt.h"
#include "dNdPt/AlidNdPtEventCuts.h"
#include "dNdPt/AlidNdPtAcceptanceCuts.h"

#include "dNdPt/AlidNdPtTask.h"

using namespace std;

ClassImp(AlidNdPtTask)

//_____________________________________________________________________________
AlidNdPtTask::AlidNdPtTask(const char *name) 
  : AliAnalysisTaskSE(name)
  , fESD(0)
  , fMC(0)
  , fOutput(0)
  , fPitList(0)
  , fCompList(0)
  , fUseMCInfo(kFALSE)
  , fUseCentrality(0)    // default = off
  , fUseCentralityBin(0)
{
  // Constructor

  // Define input and output slots here
  DefineOutput(1, TList::Class());

  // create the list for comparison objects
  fCompList = new TList;
}

//_____________________________________________________________________________
AlidNdPtTask::~AlidNdPtTask()
{
  if(fOutput) delete fOutput;  fOutput =0; 
  if(fCompList) delete fCompList;  fCompList =0; 
}

//____________________________________________________________________________
Bool_t AlidNdPtTask::Notify()
{
  static Int_t count = 0;
  count++;
  //Printf("Processing %d. file: %s", count, ((TTree*) GetInputData(0))->GetCurrentFile()->GetName());
  TTree *chain = (TChain*)GetInputData(0);
  if(chain)
    Printf("Processing %d. file: %s", count, chain->GetCurrentFile()->GetName());

  /*
  TChain *chain = (TChain*)GetInputData(0);
  AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if (!esdH) {
    Printf("ERROR: Could not get ESDInputHandler");
    return kFALSE;
  } else {
    if(chain)
    Printf("chain->GetCurrentFile()->GetName() %s", chain->GetCurrentFile()->GetName());
  }
  */

return kTRUE;
}

//_____________________________________________________________________________
Bool_t AlidNdPtTask::AddAnalysisObject(AlidNdPt *pObj) 
{
  // add analysis object to the list
  if(pObj == 0) {
    Printf("ERROR: Could not add comparison object");
    return kFALSE;
  }

  // add object to the list
  fCompList->AddLast(pObj);
       
return kTRUE;
}

//_____________________________________________________________________________
void AlidNdPtTask::UserCreateOutputObjects()
{
  // Create histograms
  // Called once

  OpenFile(1, "RECREATE");

  //
  // create output list
  //
  fOutput = new TList;
  fOutput->SetOwner();
  fPitList = fOutput->MakeIterator();

  // add dNdPt analysis objects to the output
  AlidNdPt *pObj=0;
  Int_t count=0;
  TIterator *pitCompList = fCompList->MakeIterator();
  pitCompList->Reset();
  while(( pObj = (AlidNdPt *)pitCompList->Next()) != NULL) {
    fOutput->Add(pObj);
    count++;
  }
  Printf("UserCreateOutputObjects(): Number of output objects: %d \n", count);

  if(fUseCentrality) {
    Printf("Use Centrality - Bin %d", fUseCentralityBin);
  }
}

//_____________________________________________________________________________
void AlidNdPtTask::UserExec(Option_t *) 
{
  //
  // Called for each event
  //

  // ESD event
  fESD = (AliESDEvent*) (InputEvent());
  if (!fESD) {
    Printf("ERROR: ESD event not available");
    return;
  }

  // MC event
  if(fUseMCInfo) {
    fMC = MCEvent();
    if (!fMC) {
      Printf("ERROR: MC event not available");
      return;
    }
  }

  // track cuts from Jochen
  const AliESDVertex* vtxESDTPC = fESD->GetPrimaryVertexTPC();
  if( vtxESDTPC->GetNContributors() < 1 ) {   
    return;
  }

  const AliMultiplicity* multESD = fESD->GetMultiplicity();
  if( vtxESDTPC->GetNContributors() < (-10.+0.25*multESD->GetNumberOfITSClusters(0)) ) {
    return;
  } 
  
  
  // Process analysis
  
  Bool_t process = kTRUE;

  // Check fo centrality
  if (fUseCentrality) {
    if ( CalculateCentralityBin() != fUseCentralityBin ) {
      process = kFALSE;
    }
  }

  if (process) {
    AlidNdPt *pObj = 0;
    fPitList->Reset();
    while((pObj = (AlidNdPt *)fPitList->Next()) != NULL) {
      pObj->Process(fESD,fMC);
    }
  }


  // Post output data.
  PostData(1, fOutput);
}

//_____________________________________________________________________________
void AlidNdPtTask::Terminate(Option_t *) 
{
  // Called one at the end 
  
  // check output data
  fOutput = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutput) {
    Printf("ERROR: AlidNdPtTask::Terminate(): Output data not avaiable GetOutputData(0)==0x0 ..." );
    return;
  }
}

//________________________________________________________________________
Int_t AlidNdPtTask::CalculateCentralityBin(){
  // Get Centrality bin

  Int_t centrality = -1;

  if (fUseCentrality == 0)
    return centrality;

  AliESDVZERO* esdV0 = fESD->GetVZEROData();
  Float_t multV0 = esdV0->GetMTotV0A() + esdV0->GetMTotV0C();
  
  Float_t nClusters[6];

  const AliMultiplicity *mult = fESD->GetMultiplicity();
  for(Int_t ilay = 0; ilay < 6; ilay++){
   nClusters[ilay] = mult->GetNumberOfITSClusters(ilay);
  }
  
  if ( fUseCentrality == 1 ) {
    // -- centrality cuts V0
#if 0 
    // 2010-11-10 - now old cuts
    if (      multV0 >=    0.   && multV0 <=   106.75 ) centrality = 90;
    else if ( multV0 >   106.75 && multV0 <=   277.55 ) centrality = 80;
    else if ( multV0 >   277.55 && multV0 <=   661.85 ) centrality = 70;
    else if ( multV0 >   661.85 && multV0 <=  1345.05 ) centrality = 60;
    else if ( multV0 >  1345.05 && multV0 <=  2412.55 ) centrality = 50;
    else if ( multV0 >  2412.55 && multV0 <=  3907.05 ) centrality = 40;
    else if ( multV0 >  3907.05 && multV0 <=  6042.05 ) centrality = 30;
    else if ( multV0 >  6042.05 && multV0 <=  8902.95 ) centrality = 20;
    else if ( multV0 >  8902.95 && multV0 <= 12788.6  ) centrality = 10;
    else if ( multV0 > 12788.6  && multV0 <= 15222.5  ) centrality = 5;
    else if ( multV0 > 15222.5  && multV0 <= 19449.8  ) centrality = 0;
#else
    // 2010-11-14
    if (      multV0 >=    0.  && multV0 <=   124.5 ) centrality = 90;
    else if ( multV0 >   124.5 && multV0 <=   274.5 ) centrality = 80;
    else if ( multV0 >   274.5 && multV0 <=   574.5 ) centrality = 70;
    else if ( multV0 >   574.5 && multV0 <=  1224.5 ) centrality = 60;
    else if ( multV0 >  1224.5 && multV0 <=  2174.5 ) centrality = 50;
    else if ( multV0 >  2174.5 && multV0 <=  3624.5 ) centrality = 40;
    else if ( multV0 >  3624.5 && multV0 <=  5574.5 ) centrality = 30;
    else if ( multV0 >  5574.5 && multV0 <=  8274.5 ) centrality = 20;
    else if ( multV0 >  8274.5 && multV0 <= 12024.5 ) centrality = 10;
    else if ( multV0 > 12024.5 && multV0 <= 14674.5 ) centrality = 5;
    else if ( multV0 > 14674.5 && multV0 <= 19449.5 ) centrality = 0;
#endif
  }
  else if ( fUseCentrality == 2 ) {
#if 0 
    // 2010-11-10 - now old cuts
    if (      nClusters[1] >=    0.  && nClusters[1] <=    7.18 )  centrality = 100;
    else if ( nClusters[1] >    7.18 && nClusters[1] <=   35.9  )  centrality = 90;
    else if ( nClusters[1] >   35.9  && nClusters[1] <=   93.34 )  centrality = 80;
    else if ( nClusters[1] >   93.34 && nClusters[1] <=  222.58 )  centrality = 70;
    else if ( nClusters[1] >  222.58 && nClusters[1] <=  437.98 )  centrality = 60;
    else if ( nClusters[1] >  437.98 && nClusters[1] <=  768.26 )  centrality = 50;
    else if ( nClusters[1] >  768.26 && nClusters[1] <= 1242.14 )  centrality = 40;
    else if ( nClusters[1] > 1242.14 && nClusters[1] <= 1888.34 )  centrality = 30;
    else if ( nClusters[1] > 1888.34 && nClusters[1] <= 2735.58 )  centrality = 20;
    else if ( nClusters[1] > 2735.58 && nClusters[1] <= 3884.38 )  centrality = 10;
    else if ( nClusters[1] > 3884.38 && nClusters[1] <= 4573.66 )  centrality = 5;
    else if ( nClusters[1] > 4573.66 && nClusters[1] <= 6540.98 )  centrality = 0;
#else
    // 2010-11-14
    if      ( nClusters[1] >     0. && nClusters[1] <=   29.5 )  centrality = 90;
    else if ( nClusters[1] >   29.5 && nClusters[1] <=   69.5 )  centrality = 80;
    else if ( nClusters[1] >   69.5 && nClusters[1] <=  149.5 )  centrality = 70;
    else if ( nClusters[1] >  149.5 && nClusters[1] <=  309.5 )  centrality = 60;
    else if ( nClusters[1] >  309.5 && nClusters[1] <=  589.5 )  centrality = 50;
    else if ( nClusters[1] >  589.5 && nClusters[1] <=  989.5 )  centrality = 40;
    else if ( nClusters[1] >  989.5 && nClusters[1] <= 1569.5 )  centrality = 30;
    else if ( nClusters[1] > 1569.5 && nClusters[1] <= 2369.5 )  centrality = 20;
    else if ( nClusters[1] > 2369.5 && nClusters[1] <= 3509.5 )  centrality = 10;
    else if ( nClusters[1] > 3509.5 && nClusters[1] <= 4349.5 )  centrality = 5;
    else if ( nClusters[1] > 4349.5 && nClusters[1] <= 6540.5 )  centrality = 0;
#endif
  }
  
  return centrality;
}
