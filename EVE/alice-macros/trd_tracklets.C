#include "TEveLine.h"

TEveElementList *trd_tracklets()
{
  AliRunLoader* rl =  AliEveEventManager::AssertRunLoader();
  AliLoader *loader = rl ? rl->GetLoader("TRDLoader") : 0x0;
  AliDataLoader *dl = loader ? loader->GetDataLoader("tracklets") : 0x0;
  if (!dl)
    return;

  dl->Load();
  TTree *trklTree = dl->Tree();

  TBranch *trklBranch = 0x0;

  gEve->DisableRedraw();

  if (trklBranch = trklTree->GetBranch("trkbranch")) {
    TEveElementList* listOfTracklets = new TEveElementList("Online tracklets");
    gEve->AddElement(listOfTracklets);

    UInt_t *leaves = new UInt_t[258];
    trklBranch->SetAddress(leaves);

    for (Int_t iEntry = 0; iEntry < trklBranch->GetEntries(); iEntry++) {
      trklBranch->GetEntry(iEntry);
      for (Int_t iTracklet = 0; iTracklet < 256; iTracklet++) {
        if (leaves[2 + iTracklet] == 0)
          break;
        AliEveTRDTrackletOnline *evetrkl = new AliEveTRDTrackletOnline(new AliTRDtrackletWord(leaves[2 + iTracklet], leaves[0] + leaves[1]));
        gEve->AddElement(evetrkl, listOfTracklets);
      }
    }
    delete [] leaves;
  }

  if (trklBranch = trklTree->GetBranch("mcmtrklbranch")) {
    AliTRDtrackletMCM *trkl = 0x0; //new AliTRDtrackletMCM;
    trklBranch->SetAddress(&trkl);

    TEveElementList* listOfTracklets = new TEveElementList("MCM tracklets");
    gEve->AddElement(listOfTracklets);

    for (Int_t i = 0; i < trklBranch->GetEntries(); i++) {
      trklBranch->GetEntry(i);
      if (!trkl)
	continue;
      gEve->AddElement(new AliEveTRDTrackletOnline(trkl), listOfTracklets);
    }
  }

  gEve->EnableRedraw();
  gEve->Redraw3D();

  return 0x0;
}

