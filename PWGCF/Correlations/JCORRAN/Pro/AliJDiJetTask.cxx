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

//______________________________________________________________________________
// Analysis task for di-jet Analysis
// author: R.Diaz, J. Rak,  D.J. Kim
// ALICE Group University of Jyvaskyla 
// Finland 
// Fill the analysis containers for ESD or AOD
// Adapted for AliAnalysisTaskSE and AOD objects  
//////////////////////////////////////////////////////////////////////////////

#include "AliAnalysisTaskSE.h"
#include "AliAODHandler.h"
#include "AliAnalysisManager.h"
#include "AliLog.h"
#include "AliJDiJetTask.h" 
#include "AliJEventHeader.h"
#include "AliJRunHeader.h"
#include "AliJCORRANTask.h"
#include "AliJCard.h"
#include "AliJRunTable.h"
#include "AliAnalysisUtils.h"

//______________________________________________________________________________
AliJDiJetTask::AliJDiJetTask() :   
  AliAnalysisTaskSE("AliJDiJetTaskTask"),
  fJetTask(NULL),
  fJetTaskName(""),
  fJDiJetAnalysis(0x0),
  fOutput(NULL),
  fFirstEvent(kTRUE),
  fAnaUtils(NULL),
  fRunTable(NULL),
  fCard(NULL)
{
  DefineOutput (1, TDirectory::Class());
}

//______________________________________________________________________________
AliJDiJetTask::AliJDiJetTask(const char *name, TString inputformat):
  AliAnalysisTaskSE(name), 
  fJetTask(NULL),
  fJetTaskName(""),
  fJDiJetAnalysis(0x0),
  fOutput(NULL),
  fFirstEvent(kTRUE),
  fAnaUtils(NULL),
  fRunTable(NULL),
  fCard(NULL)
{
  // Constructor
  // AliInfo("---- AliJDiJetTask Constructor ----");

  JUNUSED(inputformat);
  DefineOutput (1, TDirectory::Class());
}

//____________________________________________________________________________
AliJDiJetTask::AliJDiJetTask(const AliJDiJetTask& ap) :
  AliAnalysisTaskSE(ap.GetName()), 
  fJetTask(ap.fJetTask),
  fJetTaskName(ap.fJetTaskName),
  fJDiJetAnalysis( ap.fJDiJetAnalysis ),
  fOutput( ap.fOutput ),
  fFirstEvent( ap.fFirstEvent),
  fAnaUtils(ap.fAnaUtils),
  fRunTable(ap.fRunTable),
  fCard( ap.fCard )
{ 

  //AliInfo("----DEBUG AliJDiJetTask COPY ----");

}

//_____________________________________________________________________________
AliJDiJetTask& AliJDiJetTask::operator = (const AliJDiJetTask& ap)
{
  // assignment operator

  //AliInfo("----DEBUG AliJDiJetTask operator= ----");
  this->~AliJDiJetTask();
  new(this) AliJDiJetTask(ap);
  return *this;
}

//______________________________________________________________________________
AliJDiJetTask::~AliJDiJetTask()
{
  // destructor 

  delete fJDiJetAnalysis;
  delete fAnaUtils;

}

//________________________________________________________________________

void AliJDiJetTask::UserCreateOutputObjects()
{  
  //=== create the jcorran outputs objects
  if(fDebug > 1) printf("AliJDiJetTask::UserCreateOutPutData() \n");

  fAnaUtils = new AliAnalysisUtils();
  fAnaUtils->SetUseOutOfBunchPileUp( kTRUE );

  //=== Get AnalysisManager
  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();

  OpenFile(1);
  fOutput = gDirectory;//->mkdir("JDiHadronCorr");
  fOutput->cd();

  fJDiJetAnalysis = new AliJDiJetAnalysis(fCard);
  fJDiJetAnalysis->UserCreateOutputObjects();
  fCard->WriteCard( gDirectory );

  PostData( 1, fOutput );

  fJetTask = (AliJJetTask*)(man->GetTask( fJetTaskName));

  for( int i=0;i<fJetTask->GetNumberOfJetCollections(); i++ ){
    fJDiJetAnalysis->AddJets( fJetTask->GetAliJJetList( i ), fJetTask->GetTrackOrMCParticle(i) );
    // TODO: fChargedOrFull;
    // TODO: Raidus
  }

  //fJDiJetAnalysis->CreateHistos();


  cout << "Add(fAliRunHeader) in UserCreateObject() ======= " << endl;
}

//______________________________________________________________________________
void AliJDiJetTask::UserExec(Option_t* /*option*/) 
{

  // Processing of one event
  if(fDebug > 5) cout << "------- AliJDiJetTask Exec-------"<<endl;

  // Check Event
  if( fJetTask->GetTaskEntry() != fEntry ) return;
  AliVEvent *event = InputEvent();
  if(!event) return;
  AliAODEvent* aodEvent = dynamic_cast<AliAODEvent*>(event);
  if(!aodEvent) return;
  if( fFirstEvent ) {
    fRunTable = & AliJRunTable::GetSpecialInstance();
    fRunTable->SetRunNumber( aodEvent->GetRunNumber() );
    fFirstEvent = kFALSE;
  }
  if(!IsGoodEvent( event )) return; // zBin is set there

  // Call DiJetAnalysis
  fJDiJetAnalysis->ClearBeforeEvent();
  fJDiJetAnalysis->UserExec();
  PostData(1, fOutput );

  if(fDebug > 5) cout << "\t------- End UserExec "<<endl;
}

//______________________________________________________________________________
void AliJDiJetTask::Init()
{
  // Intialisation of parameters
  AliInfo("Doing initialization") ; 
  //fJDiJetAnalysis->Init();
}

//______________________________________________________________________________
void AliJDiJetTask::Terminate(Option_t *)
{
  cout<<"AliJDiJetTask Analysis DONE !!"<<endl; 

}

//________________________________________________________________________
bool AliJDiJetTask::IsGoodEvent(AliVEvent *event) {

  // TODO pile up test for PP
  if(fRunTable->IsPP() && fAnaUtils->IsPileUpEvent(event)) {
    return kFALSE;
  } else {
    Bool_t triggeredEventMB = kFALSE; //init

    Bool_t triggerkMB = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & ( AliVEvent::kMB );

    if( triggerkMB ){
      triggeredEventMB = kTRUE;  //event triggered as minimum bias
    }
    //--------------------------------------------------------------
    // check reconstructed vertex
    int ncontributors = 0;
    Bool_t goodRecVertex = kFALSE;
    const AliVVertex *vtx = event->GetPrimaryVertex();
    if(vtx){
      ncontributors = vtx->GetNContributors();
      if(ncontributors > 0){
        double zVert = vtx->GetZ();
        if(fCard->VertInZRange(zVert)) {
          goodRecVertex = kTRUE;
        }
      }
    }
    return goodRecVertex;
  }
  //---------------------------------
}
