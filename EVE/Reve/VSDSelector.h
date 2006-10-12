#ifndef REVE_VSDSelector_H
#define REVE_VSDSelector_H

#include "RGBrowser.h"
#include <Reve/VSD.h>

#include <TGTextEntry.h>

namespace Reve {

class RenderElement;
class TrackRnrStyle;
class TrackList;

class VSDSelector : public VSD
{
  VSDSelector(const VSDSelector&);            // Not implemented
  VSDSelector& operator=(const VSDSelector&); // Not implemented

protected:
  TGTextEntry*              mParticleSelection;
  TGCheckButton*            fRecursiveSelect;

  TGTextEntry*              mHitSelection;
  TGTextEntry*              mClusterSelection;
  TGTextEntry*              mRecSelection;

public: 

  VSDSelector(TGCompositeFrame *tFrame);

  virtual void LoadVSD(const Text_t* vsd_file_name,
                       const Text_t* dir_name="Event0");

  void SelectParticles (const Text_t* selection=0);
  void ImportDaughtersRec(RenderElement* parent, TrackList* cont,
			  Int_t first, Int_t last);

  void SelectHits();
  void SelectClusters();
  void SelectRecTracks();

  void SetRecursiveSelection(Bool_t rec){fRecursiveSelect->SetOn(rec,1);}
  //      printf("SetRecursiveSelection is %d on %d \n", rec?1:0,fRecursiveSelect->IsOn()?1:0);}
  Bool_t GetRecursiveSelection(){return fRecursiveSelect->IsOn();}
  ClassDef(VSDSelector, 1);
};

}

#endif
