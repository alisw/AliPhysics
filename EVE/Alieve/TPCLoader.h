// $Header$

#ifndef ALIEVE_TPCLoader_H
#define ALIEVE_TPCLoader_H

#include <Reve/RenderElement.h>

class AliRawReaderRoot;

namespace Alieve {

class TPCData;
class TPCSector2D;
class TPCSector3D;

class TPCLoader : public Reve::RenderElementList
{
  friend class TPCLoaderEditor;

protected:
  TString           fFile;
  Int_t             fEvent;
  Bool_t            fDoubleSR;

  AliRawReaderRoot* fReader;
  TPCData*          fData;

  std::vector<TPCSector2D*> fSec2Ds;
  std::vector<TPCSector3D*> fSec3Ds;

  Bool_t   fSetInitSectorParams;
  Int_t    fInitMinTime;
  Int_t    fInitMaxTime;
  Int_t    fInitThreshold;

public:
  TPCLoader(const Text_t* n="TPCLoader", const Text_t* t=0);
  virtual ~TPCLoader();

  virtual Bool_t CanEditMainColor() { return kFALSE; }

  void SetFile(const Text_t* f) { fFile = f; }
  void SetDoubleSR(Bool_t d)    { fDoubleSR = d; }

  TPCData* GetData() const { return fData; }
  void     SetData(TPCData* d);

  void OpenFile();
  void LoadEvent();
  void NextEvent(Bool_t rewindOnEnd=kTRUE);
  void GotoEvent(Int_t event);

  void UpdateSectors();
  void CreateSectors3D();
  void DeleteSectors3D();

  void SetInitParams(Int_t mint, Int_t maxt, Int_t thr);

  ClassDef(TPCLoader, 1);
}; // endclass TPCLoader

}

#endif
