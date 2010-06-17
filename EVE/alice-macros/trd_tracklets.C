#ifdef __CINT__
class TEveLine;
#else
#include <TEveManager.h>
#include "TEveLine.h"
#include "TClonesArray.h"
#include <EveBase/AliEveEventManager.h>

#include "AliRunLoader.h"
#include "AliLoader.h"
#include "AliDataLoader.h"
#include "AliTreeLoader.h"
#include "TRD/AliTRDarrayADC.h"
#include "EveDet/AliEveTRDData.h"
#include "TRD/AliTRDtrackletWord.h"
#include "TRD/AliTRDtrackletMCM.h"
#endif

TEveElementList *trd_tracklets()
{
  AliRunLoader* rl =  AliEveEventManager::AssertRunLoader();
  AliLoader *loader = rl ? rl->GetLoader("TRDLoader") : 0x0;

  TTree *trklTree = 0x0;

  AliDataLoader *dl = loader ? loader->GetDataLoader("tracklets") : 0x0;
  if (!dl) {
    printf("No tracklet loader\n");
    return 0x0;
  }

  gEve->DisableRedraw();

  // ----- simulated tracklets -----
  dl->Load();
  trklTree = dl->Tree();
  
  if (trklTree) {
    TBranch *trklBranch = 0x0;
    if ((trklBranch = trklTree->GetBranch("mcmtrklbranch"))) {
      AliTRDtrackletMCM *trkl = 0x0; 
      trklBranch->SetAddress(&trkl);
      
      TEveElementList* listOfTracklets = new TEveElementList("TRD tracklets (sim)");
      gEve->AddElement(listOfTracklets);
      
      for (Int_t i = 0; i < trklBranch->GetEntries(); i++) {
	trklBranch->GetEntry(i);
	if (!trkl)
	  continue;
	gEve->AddElement(new AliEveTRDTrackletOnline(trkl), listOfTracklets);
      }
    }
  }

  // raw tracklets
  AliTreeLoader *tl = (AliTreeLoader*) dl->GetBaseLoader("tracklets-raw");
  if (tl) {
    tl->Load();
    trklTree = tl->Tree();
  }
  else 
    trklTree = 0x0;
  //  trklTree = tl ? tl->Load(), tl->Tree : 0x0;

  if (trklTree) {
    TEveElementList* listOfTracklets = new TEveElementList("TRD tracklets (raw)");
    gEve->AddElement(listOfTracklets);
    
    Int_t hc; 
    TClonesArray *ar = 0x0;
    trklTree->SetBranchAddress("hc", &hc);
    trklTree->SetBranchAddress("trkl", &ar);

    for (Int_t iEntry = 0; iEntry < trklTree->GetEntries(); iEntry++) {
      trklTree->GetEntry(iEntry);
      //      printf("%i tracklets in HC %i\n", ar->GetEntriesFast(), hc);
      for (Int_t iTracklet = 0; iTracklet < ar->GetEntriesFast(); iTracklet++) {
	AliTRDtrackletWord *trklWord = (AliTRDtrackletWord*) (*ar)[iTracklet];
        AliEveTRDTrackletOnline *evetrkl = new AliEveTRDTrackletOnline(new AliTRDtrackletWord(trklWord->GetTrackletWord(), hc));
        gEve->AddElement(evetrkl, listOfTracklets);
      }
    }
  }

  gEve->EnableRedraw();
  gEve->Redraw3D();

  return 0x0;
}

