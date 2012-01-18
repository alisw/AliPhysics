// $Id$
//
// Task to skim ESD files.
//
//

#include "AliEsdSkimTask.h"
#include <TClonesArray.h>
#include <TFile.h>
#include <TTree.h>
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliESDtrackCuts.h"
#include "AliMultiplicity.h"

//_________________________________________________________________________________________________
AliEsdSkimTask::AliEsdSkimTask(const char *opt) :
  AliAnalysisTaskSE(opt), fEvent(0), fTree(0), fCuts(0),
  fDoZDC(1), fDoV0(1), fDoT0(1), fDoTPCv(1), fDoSPDv(1), fDoPriv(1),
  fDoEmCs(1), fDoPCs(1), fDoEmT(1), fDoPT(1), fDoTracks(1), fDoMult(1),
  fDoTof(1), fDoPileup(1), fDoClus(1), fEmcNames(""), 
  fDoMiniTracks(0), fTracks("Tracks"), fPhosClusOnly(0)
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

  AliAnalysisManager *am = AliAnalysisManager::GetAnalysisManager();

  fEvent->Reset();

  TList* objsin  = esdin->GetList();
  TList* objsout = fEvent->GetList();

  AliESDHeader *header = dynamic_cast<AliESDHeader*>(objsin->FindObject("AliESDHeader"));
  if (header) {
    am->LoadBranch("AliESDHeader.");
    *header = *esdin->GetHeader();
  }
  AliESDRun *run = dynamic_cast<AliESDRun*>(objsin->FindObject("AliESDRun"));
  if (run) {
    am->LoadBranch("AliESDRun.");
    *run = *esdin->GetESDRun();
  }
  AliESDZDC *zdc = dynamic_cast<AliESDZDC*>(objsin->FindObject("AliESDZDC"));
  if (zdc) {
    am->LoadBranch("AliESDZDC.");
    *zdc = *esdin->GetESDZDC();
  }
  AliESDVZERO *v0 = dynamic_cast<AliESDVZERO*>(objsin->FindObject("AliESDVZERO"));
  if (v0) {
    am->LoadBranch("AliESDVZERO.");
    *v0 = *esdin->GetVZEROData();
  }
  AliESDTZERO *t0 = dynamic_cast<AliESDTZERO*>(objsin->FindObject("AliESDTZERO"));
  if (t0) {
    am->LoadBranch("AliESDTZERO.");
    *t0 = *esdin->GetESDTZERO();
  }
  AliESDVertex *tpcv = dynamic_cast<AliESDVertex*>(objsin->FindObject("TPCVertex"));
  if (tpcv) {
    am->LoadBranch("TPCVertex.");
    *tpcv = *esdin->GetPrimaryVertexTPC();
  }
  AliESDVertex *spdv = dynamic_cast<AliESDVertex*>(objsin->FindObject("SPDVertex"));
  if (spdv) {
    am->LoadBranch("SPDVertex.");
    *spdv = *esdin->GetPrimaryVertexSPD();
  }
  AliESDVertex *priv = dynamic_cast<AliESDVertex*>(objsin->FindObject("PrimaryVertex"));
  if (priv) {
    am->LoadBranch("PrimaryVertex.");
    *priv = *esdin->GetPrimaryVertexTracks();
  }
  AliESDCaloCells *ecells = dynamic_cast<AliESDCaloCells*>(objsin->FindObject("EMCALCells"));
  if (ecells) {
    am->LoadBranch("EMCALCells.");
    *ecells = *esdin->GetEMCALCells();
  }
  AliESDCaloCells *pcells = dynamic_cast<AliESDCaloCells*>(objsin->FindObject("PHOSCells"));
  if (pcells) {
    am->LoadBranch("PHOSCells.");
    *pcells = *esdin->GetPHOSCells();
  }
  AliESDCaloTrigger *etrig = dynamic_cast<AliESDCaloTrigger*>(objsin->FindObject("EMCALTrigger"));
  if (etrig) {
    am->LoadBranch("EMCALTrigger.");
    *etrig = *esdin->GetCaloTrigger("EMCAL");
  }
  AliESDCaloTrigger *ptrig = dynamic_cast<AliESDCaloTrigger*>(objsin->FindObject("PHOSTrigger"));
  if (ptrig) {
    am->LoadBranch("PHOSTrigger.");
    *ptrig = *esdin->GetCaloTrigger("PHOS");
  }

  AliMultiplicity *mult = dynamic_cast<AliMultiplicity*>(objsin->FindObject("AliMultiplicity"));
  if (mult) {
    am->LoadBranch("AliMultiplicity.");
    *mult = *esdin->GetMultiplicity();
  }

  AliTOFHeader *tofh = dynamic_cast<AliTOFHeader*>(objsin->FindObject("AliTOFHeader"));
  if (tofh) {
    am->LoadBranch("AliTOFHeader.");
    *tofh = *esdin->GetTOFHeader();
  }
  TClonesArray *spup = dynamic_cast<TClonesArray*>(objsin->FindObject("SPDPileupVertices"));
  if (spup) {
    am->LoadBranch("SPDPileupVertices");
    Int_t N = esdin->GetNumberOfPileupVerticesSPD();
    for (Int_t i=0; i<N; ++i) {
      const AliESDVertex *vtx = esdin->GetPileupVertexSPD(i);
      if (vtx)
        fEvent->AddPileupVertexSPD(vtx);
    }
  }
  TClonesArray *tpup = dynamic_cast<TClonesArray*>(objsin->FindObject("TrkPileupVertices"));
  if (tpup) {
    am->LoadBranch("TrkPileupVertices");
    Int_t N = esdin->GetNumberOfPileupVerticesTracks();
    for (Int_t i=0; i<N; ++i) {
      const AliESDVertex *vtx = esdin->GetPileupVertexTracks(i);
      if (vtx)
        fEvent->AddPileupVertexTracks(vtx);
    }
  }
  TClonesArray *clus = dynamic_cast<TClonesArray*>(objsin->FindObject("CaloClusters"));
  if (clus) {
    am->LoadBranch("");
    Int_t N = esdin->GetNumberOfCaloClusters();
    for (Int_t i=0; i<N; ++i) {
      AliESDCaloCluster *c = esdin->GetCaloCluster(i);
      if (fPhosClusOnly && c->IsEMCAL())
        continue;
      if (c)
        fEvent->AddCaloCluster(c);
    }
  }
  TObjArray *namearr = fEmcNames.Tokenize(";");
  if (namearr) {
    for (Int_t i=0; i<namearr->GetEntries(); ++i) {
      TString cname(namearr->At(i)->GetName());
      if (cname.Length()<=0)
        continue;
      TClonesArray *arrin  = dynamic_cast<TClonesArray*>(objsin->FindObject(cname));
      TClonesArray *arrout = dynamic_cast<TClonesArray*>(objsout->FindObject(cname));
      //AliFatal(Form("Can not find tracks with name %s", fTracks.Data()));
      arrout->Delete();
      const Int_t N = arrin->GetEntries();
      for (Int_t iC=0, nC=0; iC<N; ++iC) {
        AliESDCaloCluster *c = dynamic_cast<AliESDCaloCluster*>(arrin->At(iC));
        if (!c)
          continue;
        new ((*arrout)[nC++]) AliESDCaloCluster(*c);
      }
    }
    delete namearr;
  }
  if (fDoTracks) {
    am->LoadBranch("Tracks");
    TClonesArray *tracks = dynamic_cast<TClonesArray*>(objsin->FindObject(fTracks));
    if (!tracks) {
      AliFatal(Form("Can not find tracks with name %s", fTracks.Data()));
      return;
    }
    const Int_t Ntracks = tracks->GetEntries();
    Int_t nacc = 0;
    for (Int_t iTracks = 0; iTracks < Ntracks; ++iTracks) {
      AliESDtrack *track = dynamic_cast<AliESDtrack*>(tracks->At(iTracks));
      if (!track)
        continue;
      if (fCuts) {
        if (!fCuts->IsSelected(track))
          continue;
      }
      if (fDoMiniTracks) {
        track->MakeMiniESDtrack();
      }
      fEvent->AddTrack(track);
      ++nacc;
    }
    if (fCuts) 
      AliInfo(Form("Selected %d out of %d \n", nacc, Ntracks));
  }
  fTree->Fill();
}

//_________________________________________________________________________________________________
void AliEsdSkimTask::UserCreateOutputObjects() 
{
  // Create output objects.

  fTree = new TTree("esdTree", "Tree with skimmed ESD objects");
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
  TObjArray *namearr = fEmcNames.Tokenize(";");
  if (namearr) {
    for (Int_t i=0; i<namearr->GetEntries(); ++i) {
      TString cname(namearr->At(i)->GetName());
      if (cname.Length()<=0)
        continue;
      TClonesArray *arr = new TClonesArray("AliESDCaloCluster",0);
      arr->SetName(cname);
      fEvent->AddObject(arr);
    }
    delete namearr;
  }
  if (fDoTracks) {
    TClonesArray *arr = new TClonesArray("AliESDtrack",0);
    arr->SetName("Tracks");
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
