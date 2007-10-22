// $Header$

#ifndef ALIEVE_ITSDigitsInfo_H
#define ALIEVE_ITSDigitsInfo_H

#include <Reve/Reve.h>

#include <map>
#include <vector>

#include <TObject.h>
#include <TClonesArray.h>
#include <TTree.h>

#include <AliITS.h>
#include <AliITSgeom.h>
#include <AliITSsegmentationSPD.h>
#include <AliITSsegmentationSDD.h>
#include <AliITSsegmentationSSD.h>

class AliRawReader;

namespace Alieve {
/**************************************************************************/
// ITSModuleSelection
/**************************************************************************/
class ITSModuleSelection
{
public:
  Int_t    fType;
  Int_t    fLayer;
  Float_t  fMinPhi;
  Float_t  fMaxPhi;
  Float_t  fMinTheta;
  Float_t  fMaxTheta; 
  
  ITSModuleSelection();
  virtual ~ITSModuleSelection() {}

  ClassDef(ITSModuleSelection, 1);
};

/**************************************************************************/
// ITSDigitsInfo
/**************************************************************************/
class ITSDigitsInfo : public TObject, public Reve::ReferenceCount
{
  ITSDigitsInfo(const ITSDigitsInfo&);            // Not implemented
  ITSDigitsInfo& operator=(const ITSDigitsInfo&); // Not implemented

private:
  Float_t fSPDZCoord[192];

  void InitInternals();

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

  Int_t                    fSPDMinVal;
  Int_t                    fSSDMinVal;
  Int_t                    fSDDMinVal;
  Int_t                    fSPDMaxVal;
  Int_t                    fSSDMaxVal;
  Int_t                    fSDDMaxVal;

  Int_t                    fSPDHighLim;
  Int_t                    fSDDHighLim;
  Int_t                    fSSDHighLim;

  Int_t                    fSPDScaleX[5];
  Int_t                    fSPDScaleZ[5];
  Int_t                    fSDDScaleX[5];
  Int_t                    fSDDScaleZ[5];
  Int_t                    fSSDScale [5];

  ITSDigitsInfo();
  virtual ~ITSDigitsInfo();

  void SetTree(TTree* tree);
  void ReadRaw(AliRawReader* raw);

  TClonesArray* GetDigits(Int_t moduleID, Int_t detector);

  void GetSPDLocalZ(Int_t j, Float_t& z);

  void GetModuleIDs(ITSModuleSelection* sel, std::vector<UInt_t>& ids);

  virtual void Print(Option_t* opt="") const;

  ClassDef(ITSDigitsInfo, 1);
}; // endclass ITSDigitsInfo

}
#endif
