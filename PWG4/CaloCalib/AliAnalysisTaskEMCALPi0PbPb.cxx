// $Id$

#include "AliAnalysisTaskEMCALPi0PbPb.h"
#include <TChain.h>
#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TList.h>
#include <TLorentzVector.h>
#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskEMCALClusterizeFast.h"
#include "AliCentrality.h"
#include "AliEMCALGeoUtils.h"
#include "AliESDEvent.h"
#include "AliESDVertex.h"

ClassImp(AliAnalysisTaskEMCALPi0PbPb)

//________________________________________________________________________
AliAnalysisTaskEMCALPi0PbPb::AliAnalysisTaskEMCALPi0PbPb() 
  : AliAnalysisTaskSE(),
    fCentVar(),
    fCentFrom(0),
    fCentTo(100),
    fVtxZMin(-7),
    fVtxZMax(+7),
    fClusName(),
    fOutput(0),
    fEsdEv(0),
    fAodEv(0),
    fEsdClusters(0),
    fEsdCells(0),
    fAodClusters(0),
    fAodCells(0)
{
  // ROOT constructor.
}

//________________________________________________________________________
AliAnalysisTaskEMCALPi0PbPb::AliAnalysisTaskEMCALPi0PbPb(const char *name) 
  : AliAnalysisTaskSE(name),
    fCentVar("V0M"),
    fCentFrom(0),
    fCentTo(100),
    fVtxZMin(-7),
    fVtxZMax(+7),
    fClusName(),
    fOutput(0),
    fEsdEv(0),
    fAodEv(0),
    fEsdClusters(0),
    fEsdCells(0),
    fAodClusters(0),
    fAodCells(0)
{
  // Constructor.

  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  fBranchNames="ESD:AliESDRun.,AliESDHeader.,PrimaryVertex,EMCALCells.,CaloClusters "
               "AOD:header,vertices,emcalCells,caloClusters";
}

//________________________________________________________________________
AliAnalysisTaskEMCALPi0PbPb::~AliAnalysisTaskEMCALPi0PbPb()
{
  // Destructor.

  delete fOutput; fOutput = 0;
}

//________________________________________________________________________
void AliAnalysisTaskEMCALPi0PbPb::UserCreateOutputObjects()
{
  // Create user objects here.

  fOutput = new TList();
  fOutput->SetOwner();

  PostData(1, fOutput); 
}

//________________________________________________________________________
void AliAnalysisTaskEMCALPi0PbPb::UserExec(Option_t *) 
{
  // Called for each event.

  if (!InputEvent())
    return;

  const AliCentrality *centP = InputEvent()->GetCentrality();
  Double_t cent = centP->GetCentralityPercentileUnchecked(fCentVar);
  if (cent<fCentFrom||cent>fCentTo)
    return; //todo bookeeping
  // todo test quality flag

  const AliVVertex *vertex =  InputEvent()->GetPrimaryVertex();
  if (!vertex)
    return;

  if(vertex->GetZ()<fVtxZMin||vertex->GetZ()>fVtxZMax)
    return;
  //todo bookkeeping

  AliAnalysisManager *am = AliAnalysisManager::GetAnalysisManager();
  fEsdEv = dynamic_cast<AliESDEvent*>(InputEvent());
  fAodEv = dynamic_cast<AliAODEvent*>(InputEvent());

  fEsdClusters = 0;
  fEsdCells    = 0;
  fAodClusters = 0;
  fAodCells    = 0;

  // deal with AOD output first
  if (AODEvent() && !fClusName.IsNull()) {
    TList *l = AODEvent()->GetList();
    TClonesArray *clus = 0;
    if (l) {
      clus = dynamic_cast<TClonesArray*>(l->FindObject(fClusName));
      fAodClusters = clus;
    }
  }
//todo need to figure out how to get modified cells

  if (fEsdEv) { // ESD input mode
    if (!fEsdClusters) {
      AliAnalysisTaskEMCALClusterizeFast *cltask = 0;
      TObjArray *ts = am->GetTasks();
      for (Int_t i=0; i<ts->GetEntries(); ++i) {
        cltask = dynamic_cast<AliAnalysisTaskEMCALClusterizeFast*>(ts->At(i));
        if (cltask && cltask->GetClusters()) {
          fEsdClusters = cltask->GetClusters();
          break;
        }
      }
    }
    if (!fEsdClusters) {
      am->LoadBranch("CaloClusters");
      TList *l = fEsdEv->GetList();
      TClonesArray *clus = 0;
      if (l) {
        clus = dynamic_cast<TClonesArray*>(l->FindObject("CaloClusters"));
        fEsdClusters = clus;
      }
    }
    if (!fEsdCells) {
      am->LoadBranch("EMCALCells.");
      fEsdCells = fEsdEv->GetEMCALCells();
    }
  } else if (fAodEv) { // AOD input mode
    if (!fAodClusters) {
      am->LoadBranch("caloClusters");
      TList *l = fAodEv->GetList();
      TClonesArray *clus = 0;
      if (l) {
        clus = dynamic_cast<TClonesArray*>(l->FindObject("caloClusters"));
        fAodClusters = clus;
      }
    }
    if (!fAodCells) {
      am->LoadBranch("emcalCells");
      fAodCells = fAodEv->GetEMCALCells();
    }
  } else {
    AliFatal("Impossible to not have either pointer to ESD or AOD event");
  }
  
#if 0
  AliAnalysisManager *am = AliAnalysisManager::GetAnalysisManager();
  AliAODHandler *oh = dynamic_cast<AliAODHandler*>(am->GetOutputEventHandler());
  if (!oh)
    return;
  AliAODEvent *aod = oh->GetAOD();
  if (!aod)
    return;

  for (Int_t i=0; i<10; ++i) {
    TClonesArray *myclusters = fMyExClusters[i];
    if (!myclusters)
      break;
    myclusters->Clear();

    TClonesArray *arr = dynamic_cast<TClonesArray*>(aod->FindListObject(fMyExClusters[i]->GetName()));
    if (!arr) 
      continue;

    virtual AliVEvent*      {return fInputEvent;}
    virtual AliESDfriend* ESDfriend()   {return fESDfriend; }
    virtual AliAODEvent*  AODEvent()    {return fOutputAOD; }
    virtual TTree*        OutputTree()  {return fTreeA;     }
    virtual AliMCEvent*   MCEvent()     {return fMCEvent;   }


  // Main loop
  // Called for each event
      TRefArray * clusterListESD = new TRefArray();
      event->GetEMCALClusters(clusterListESD); 
      TClonesArray * clusterListAOD = new TClonesArray();
      if(AODEvent()){clusterListAOD = dynamic_cast<TClonesArray*> (AODEvent()->FindListObject("newEMCALClusters"));}

#endif
  PostData(1, fOutput);
}      

//________________________________________________________________________
void AliAnalysisTaskEMCALPi0PbPb::Terminate(Option_t *) 
{
  // Terminate called at the end of analysis.
}

