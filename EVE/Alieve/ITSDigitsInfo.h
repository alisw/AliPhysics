// $Header$

#ifndef ALIEVE_ITSDigitsInfo_H
#define ALIEVE_ITSDigitsInfo_H

#include <Reve/VSD.h>

#include <map>

#include <TObject.h>
#include <TClonesArray.h>
#include <TTree.h>

#include <AliITS.h>
#include <AliITSgeom.h>
#include <AliITSsegmentationSPD.h>
#include <AliITSsegmentationSDD.h>
#include <AliITSsegmentationSSD.h>

static const int NSCALE = 5;

namespace Alieve {

class ITSDigitsInfo : public TObject
{
private:
  void Init();
  Float_t fSPDZCoord[192];

protected:
  Int_t                      fRefCount;

  map<Int_t,  TClonesArray*> fSPDmap;
  map<Int_t,  TClonesArray*> fSDDmap;
  map<Int_t,  TClonesArray*> fSSDmap;

  void        SetITSSegmentation();

public:
  TTree*                   fTree;
  AliITSgeom*              fGeom;
  AliITSsegmentationSPD*   fSegSPD;
  AliITSsegmentationSDD*   fSegSDD;
  AliITSsegmentationSSD*   fSegSSD;

  Int_t        fSPDScaleX[NSCALE];
  Int_t        fSPDScaleZ[NSCALE];
  Int_t        fSDDScaleX[NSCALE];
  Int_t        fSDDScaleZ[NSCALE];
  Int_t        fSSDScale[NSCALE];
    
  ITSDigitsInfo(const Text_t* /*n*/="ITSDigitsInfo", const Text_t* /*t*/=0) :
    TObject()
  { Init(); } 
 
  virtual ~ITSDigitsInfo();

  void SetTree(TTree* tree);
 
  TClonesArray* GetDigits(Int_t moduleID, Int_t detector);

  void GetSPDLocalZ(Int_t j, Float_t& z);


  void IncRefCount() { ++fRefCount; }
  void DecRefCount() { --fRefCount; if(fRefCount <= 0) delete this; }

  virtual void Print(Option_t* opt="") const;

  ClassDef(ITSDigitsInfo, 1);
}; // endclass ITSDigitsInfo

}
#endif
