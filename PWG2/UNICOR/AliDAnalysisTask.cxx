//Author: Dariusz Miskowiec 2007

//=============================================================================
// my analysis task
//=============================================================================
#include "TChain.h"
#include "AliESDInputHandler.h"
#include "AliAnalysisManager.h"
#include "AliDAnalysisTask.h"
#include "AliDAnalGlobal.h"
#include "AliDAnalSingle.h"
#include "AliDAnalCorrel.h"
#include "AliDAnalPtfluc.h"
#include "AliDEventAliceESD.h"

ClassImp(AliDAnalysisTask)

/*****************************************************************************/
AliDAnalysisTask::AliDAnalysisTask() : 
  AliAnalysisTask("dali", ""), 
  fESD(0), fEv0(0), fEv1(0), 
  fDag(0), fAll(0), fPim(0), fPip(0), 
  fCnn(0), fCpp(0), fPtf(0),
  fOutputList(0)
{
  // constructor

  fEv0 = new AliDEventAliceESD();
  DefineInput(0, TChain::Class());
  DefineOutput(0, TList::Class());
}
/*****************************************************************************/
void AliDAnalysisTask::CreateOutputObjects() 
{
  // executed once on each worker 

  fDag = new AliDAnalGlobal("dag");
  fAll = new AliDAnalSingle("all",fEv0->Etamin(),fEv0->Etamax(),0);
  fPim = new AliDAnalSingle("pim",fEv0->Etamin(),fEv0->Etamax(),-211);
  fPip = new AliDAnalSingle("pip",fEv0->Etamin(),fEv0->Etamax(), 211);
  fCnn = new AliDAnalCorrel("cnn",fEv0->Etamin(),fEv0->Etamax(),-211,-211);
  fCpp = new AliDAnalCorrel("cpp",fEv0->Etamin(),fEv0->Etamax(), 211, 211);
  fPtf = new AliDAnalPtfluc("ptf",0,0);
  fOutputList = new TList();
  fOutputList->Add(fDag);
  fOutputList->Add(fAll);
  fOutputList->Add(fPim);
  fOutputList->Add(fPip);
  fOutputList->Add(fCnn);
  fOutputList->Add(fCpp);
  fOutputList->Add(fPtf);
}
/*****************************************************************************/
void AliDAnalysisTask::ConnectInputData(Option_t *)
{
// connect ESD or AOD here
// called on each input data change.

  TTree* tree = dynamic_cast<TTree*> (GetInputData(0));
  if (!tree) {
    Printf("ERROR: Could not read chain from input slot 0");
  } else {
    fEv0->AttachTree(tree);
  } 
  AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*>
        (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());

  if (!esdH) Printf("ERROR: Could not get ESDInputHandler");
  else fESD = esdH->GetEvent();
} 
/*****************************************************************************/
void AliDAnalysisTask::Exec(Option_t */*option*/)
{
  // process one event


  if (!fESD) {
    Printf("ERROR: fESD not available");
    return;
  }
  // If AttachTree aka ReadFromTree above cannot be done then this is 
  // the alternative: shallow copy of the current alice esd event to fEv0
  // memcpy((AliESDEvent*)fEv0, fESD, sizeof(AliESDEvent)); 

  if (!fEv0->Good()) return;
  fDag->Process(fEv0);
  fPim->Process(fEv0);
  fAll->Process(fEv0);
  fPip->Process(fEv0);
  fCnn->Process(0,fEv0,fEv0,0);
  fCnn->Process(2,fEv0,fEv0,TMath::DegToRad()*180);	  
  fCpp->Process(0,fEv0,fEv0,0);
  fCpp->Process(2,fEv0,fEv0,TMath::DegToRad()*180);	  
  fPtf->Process(0,fEv0,fEv0);	  
  PostData(0, fOutputList);
} 
/*****************************************************************************/
void AliDAnalysisTask::Terminate(Option_t *)
{
  // terminate
  printf("terminate\n");
  fDag->Save("unicor-result.root","recreate");
  fAll->Save("unicor-result.root");
  fPim->Save("unicor-result.root");
  fPip->Save("unicor-result.root");
  fCnn->Save("unicor-result.root");
  fCpp->Save("unicor-result.root");
  fPtf->Save("unicor-result.root");
}
/*****************************************************************************/
