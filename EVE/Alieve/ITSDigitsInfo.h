// $Header$

#ifndef ALIEVE_ITSDigitsInfo_H
#define ALIEVE_ITSDigitsInfo_H

#include <Reve/Reve.h>

#include <map>

#include <TObject.h>
#include <TClonesArray.h>
#include <TTree.h>

#include <AliITS.h>
#include <AliITSgeom.h>
#include <AliITSsegmentationSPD.h>
#include <AliITSsegmentationSDD.h>
#include <AliITSsegmentationSSD.h>


namespace Alieve {

class ITSDigitsInfo : public TObject, public Reve::ReferenceCount
{
  ITSDigitsInfo(const ITSDigitsInfo&);            // Not implemented
  ITSDigitsInfo& operator=(const ITSDigitsInfo&); // Not implemented

private:
  Float_t fSPDZCoord[192];

protected:
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

  Int_t        fSPDScaleX[5];
  Int_t        fSPDScaleZ[5];
  Int_t        fSDDScaleX[5];
  Int_t        fSDDScaleZ[5];
  Int_t        fSSDScale [5];
    
  ITSDigitsInfo();
  virtual ~ITSDigitsInfo();

  void SetTree(TTree* tree);
 
  TClonesArray* GetDigits(Int_t moduleID, Int_t detector);

  void GetSPDLocalZ(Int_t j, Float_t& z);

  virtual void Print(Option_t* opt="") const;

  ClassDef(ITSDigitsInfo, 1);
}; // endclass ITSDigitsInfo

}
#endif
