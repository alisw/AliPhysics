// $Id$

#include "AliEsdSkimTask.h"
#include <TClonesArray.h>
#include <TFile.h>
#include <TTree.h>
#include "AliESDEvent.h"
#include "AliESDtrackCuts.h"
#include "AliMultiplicity.h"

//_________________________________________________________________________________________________
AliEsdSkimTask::AliEsdSkimTask(const char *opt) :
  AliAnalysisTaskSE(opt), fEvent(0), fTree(0), fCuts(0),
  fDoZDC(1), fDoV0(1), fDoT0(1), fDoTPCv(1), fDoSPDv(1), fDoPriv(1),
  fDoEmCs(1), fDoPCs(1), fDoEmT(1), fDoPT(1), fDoTracks(1), fDoMult(1),
  fDoTof(1), fDoPileup(1), fDoClus(1)
{
  // Constructor.

  if (!opt)
    return;

  DefineOutput(1, TTree::Class());
}

//_________________________________________________________________________________________________
void AliEsdSkimTask::UserExec(Option_t */*opt*/) 
{
  // Process event.

  AliESDEvent *esdin = dynamic_cast<AliESDEvent*>(InputEvent());
  if (!esdin)
    return;

  fEvent->Reset();

  TList* objs   = fEvent->GetList();

  AliESDHeader *header = dynamic_cast<AliESDHeader*>(objs->FindObject("AliESDHeader"));
  if (header) {
    *header = *esdin->GetHeader();
  }
  AliESDRun *run = dynamic_cast<AliESDRun*>(objs->FindObject("AliESDRun"));
  if (run) {
    *run = *esdin->GetESDRun();
  }
  AliESDZDC *zdc = dynamic_cast<AliESDZDC*>(objs->FindObject("AliESDZDC"));
  if (zdc) {
    *zdc = *esdin->GetESDZDC();
  }
  AliESDVZERO *v0 = dynamic_cast<AliESDVZERO*>(objs->FindObject("AliESDVZERO"));
  if (v0) {
    *v0 = *esdin->GetVZEROData();
  }
  AliESDTZERO *t0 = dynamic_cast<AliESDTZERO*>(objs->FindObject("AliESDTZERO"));
  if (t0) {
    *t0 = *esdin->GetESDTZERO();
  }
  AliESDVertex *tpcv = dynamic_cast<AliESDVertex*>(objs->FindObject("TPCVertex"));
  if (tpcv) {
    *tpcv = *esdin->GetPrimaryVertexTPC();
  }
  AliESDVertex *spdv = dynamic_cast<AliESDVertex*>(objs->FindObject("SPDVertex"));
  if (spdv) {
    *spdv = *esdin->GetPrimaryVertexSPD();
  }
  AliESDVertex *priv = dynamic_cast<AliESDVertex*>(objs->FindObject("PrimaryVertex"));
  if (priv) {
    *priv = *esdin->GetPrimaryVertexTracks();
  }
  AliESDCaloCells *ecells = dynamic_cast<AliESDCaloCells*>(objs->FindObject("EMCALCells"));
  if (ecells) {
    *ecells = *esdin->GetEMCALCells();
  }
  AliESDCaloCells *pcells = dynamic_cast<AliESDCaloCells*>(objs->FindObject("PHOSCells"));
  if (pcells) {
    *pcells = *esdin->GetPHOSCells();
  }
  AliESDCaloTrigger *etrig = dynamic_cast<AliESDCaloTrigger*>(objs->FindObject("EMCALTrigger"));
  if (etrig) {
    *etrig = *esdin->GetCaloTrigger("EMCAL");
  }
  AliESDCaloTrigger *ptrig = dynamic_cast<AliESDCaloTrigger*>(objs->FindObject("PHOSTrigger"));
  if (ptrig) {
    *ptrig = *esdin->GetCaloTrigger("PHOS");
  }

  AliMultiplicity *mult = dynamic_cast<AliMultiplicity*>(objs->FindObject("AliMultiplicity"));
  if (mult) {
    *mult = *esdin->GetMultiplicity();
  }

  AliTOFHeader *tofh = dynamic_cast<AliTOFHeader*>(objs->FindObject("AliTOFHeader"));
  if (tofh) {
    *tofh = *esdin->GetTOFHeader();
  }
  TClonesArray *spup = dynamic_cast<TClonesArray*>(objs->FindObject("SPDPileupVertices"));
  if (spup) {
    Int_t N = esdin->GetNumberOfPileupVerticesSPD();
    for (Int_t i=0; i<N; ++i) {
      const AliESDVertex *vtx = esdin->GetPileupVertexSPD(i);
      if (vtx)
        fEvent->AddPileupVertexSPD(vtx);
    }
  }
  TClonesArray *tpup = dynamic_cast<TClonesArray*>(objs->FindObject("TrkPileupVertices"));
  if (tpup) {
    Int_t N = esdin->GetNumberOfPileupVerticesTracks();
    for (Int_t i=0; i<N; ++i) {
      const AliESDVertex *vtx = esdin->GetPileupVertexTracks(i);
      if (vtx)
        fEvent->AddPileupVertexTracks(vtx);
    }
  }
  TClonesArray *clus = dynamic_cast<TClonesArray*>(objs->FindObject("CaloClusters"));
  if (clus) {
    Int_t N = esdin->GetNumberOfCaloClusters();
    for (Int_t i=0; i<N; ++i) {
      AliESDCaloCluster *c = esdin->GetCaloCluster(i);
      if (c)
        fEvent->AddCaloCluster(c);
    }
  }
  if (fDoTracks) {
    const Int_t Ntracks = esdin->GetNumberOfTracks();
    Int_t nacc = 0;
    for (Int_t iTracks = 0; iTracks < Ntracks; ++iTracks) {
      AliESDtrack *track = esdin->GetTrack(iTracks);
      if (!track)
        continue;
      if (fCuts) {
        if (!fCuts->IsSelected(track))
          continue;
      }
      fEvent->AddTrack(track);
      ++nacc;
    }
    if (fCuts) 
      printf("selected %d out of %d \n", nacc, Ntracks);
  }
  fTree->Fill();
}

//_________________________________________________________________________________________________
void AliEsdSkimTask::UserCreateOutputObjects() 
{
  // Create output objects.

  fTree = new TTree("esdTree", "Tree with ESD objects");
  fEvent = new AliESDEvent;
  fEvent->AddObject(new AliESDHeader());
  fEvent->AddObject(new AliESDRun());
  if (fDoZDC) 
    fEvent->AddObject(new AliESDZDC());
  if (fDoV0)
    fEvent->AddObject(new AliESDVZERO());
  if (fDoT0)
    fEvent->AddObject(new AliESDTZERO());
  if (fDoTPCv) {
    AliESDVertex *tpcv = new AliESDVertex();
    tpcv->SetName("TPCVertex");
    fEvent->AddObject(tpcv);
  }
  if (fDoSPDv) {
    AliESDVertex *spdv = new AliESDVertex();
    spdv->SetName("SPDVertex");
    fEvent->AddObject(spdv);
  }
  if (fDoPriv) {
    AliESDVertex *priv = new AliESDVertex();
    priv->SetName("PrimaryVertex");
    fEvent->AddObject(priv);
  }
  if (fDoEmCs) {
    fEvent->AddObject(new AliESDCaloCells("EMCALCells","EMCALCells"));
  }
  if (fDoPCs) {
    fEvent->AddObject(new AliESDCaloCells("PHOSCells","PHOSCells"));
  }
  if (fDoEmT) {
    AliESDCaloTrigger *etrig = new AliESDCaloTrigger;
    etrig->SetName("EMCALTrigger");
    fEvent->AddObject(etrig);
  }
  if (fDoPT) {
    AliESDCaloTrigger *ptrig = new AliESDCaloTrigger;
    ptrig->SetName("PHOSTrigger");
    fEvent->AddObject(ptrig);
  }
  if (fDoTracks) {
    TClonesArray *arr = new TClonesArray("AliESDtrack",0);
    arr->SetName("Tracks");
    fEvent->AddObject(arr);
  }
  if (fDoMult) {
    fEvent->AddObject(new AliMultiplicity());
  }
  if (fDoPileup) {
    TClonesArray *arr1 = new TClonesArray("AliESDVertex",0);
    arr1->SetName("SPDPileupVertices");
    fEvent->AddObject(arr1);
    TClonesArray *arr2 = new TClonesArray("AliESDVertex",0);
    arr2->SetName("TPCPileupVertices");
    fEvent->AddObject(arr2);
  }
  if (fDoTof) { 
    fEvent->AddObject(new AliTOFHeader());
  }
  if (fDoClus) {
    TClonesArray *arr = new TClonesArray("AliESDCaloCluster",0);
    arr->SetName("CaloClusters");
    fEvent->AddObject(arr);
  }
  fEvent->GetStdContent();
  fEvent->WriteToTree(fTree);
  fTree->GetUserInfo()->Add(fEvent);
  TFile *file = OpenFile(1);
  fTree->SetDirectory(file);
  fTree->SetAutoFlush(-1024*1024*1024);
  fTree->SetAutoSave(-1024*1024*1024);
  PostData(1,fTree);
}
