#include "TParticlePDG.h"
#include "TEveLine.h"

TEveElementList *trd_mcmtracklets()
{
  AliRunLoader* rl =  AliEveEventManager::AssertRunLoader();
  rl->LoadKinematics();
  AliLoader *loader = rl->GetLoader("TRDLoader");
  AliDataLoader *dl = loader->GetDataLoader("tracklets");
  dl->Load();
  TTree *trklTree = dl->Tree();
  TBranch *mcmBranch = trklTree->GetBranch("mcmtrklbranch");
  AliTRDtrackletMCM *trkl = new AliTRDtrackletMCM;
  mcmBranch->SetAddress(&trkl);

  gEve->DisableRedraw();
  TEveElementList* listOfTracklets = new TEveElementList("MCM tracklets");
  gEve->AddElement(listOfTracklets);

  AliTRDgeometry *geo = new AliTRDgeometry;
  for (Int_t i = 0; i < mcmBranch->GetEntries(); i++) {
    mcmBranch->GetEntry(i);
    if (!trkl)
	continue;
    AliEveTRDTracklet *evetrkl = new AliEveTRDTracklet(trkl);
    gEve->AddElement(evetrkl, listOfTracklets);

  }

  gEve->EnableRedraw();
  gEve->Redraw3D();

  return listOfTracklets;
}

