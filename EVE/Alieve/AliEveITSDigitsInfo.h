// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 * 
 **************************************************************************/

#ifndef ALIEVE_ITSDigitsInfo_H
#define ALIEVE_ITSDigitsInfo_H

#include <TEveUtil.h>

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

/**************************************************************************/
// AliEveITSModuleSelection
/**************************************************************************/
class AliEveITSModuleSelection
{
public:
  Int_t    fType;
  Int_t    fLayer;
  Float_t  fMinPhi;
  Float_t  fMaxPhi;
  Float_t  fMinTheta;
  Float_t  fMaxTheta; 
  
  AliEveITSModuleSelection();
  virtual ~AliEveITSModuleSelection() {}

  ClassDef(AliEveITSModuleSelection, 1);
};

/**************************************************************************/
// AliEveITSDigitsInfo
/**************************************************************************/
class AliEveITSDigitsInfo : public TObject, public TEveRefCnt
{
  AliEveITSDigitsInfo(const AliEveITSDigitsInfo&);            // Not implemented
  AliEveITSDigitsInfo& operator=(const AliEveITSDigitsInfo&); // Not implemented

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

  AliEveITSDigitsInfo();
  virtual ~AliEveITSDigitsInfo();

  void SetTree(TTree* tree);
  void ReadRaw(AliRawReader* raw, Int_t mode);

  TClonesArray* GetDigits(Int_t moduleID, Int_t detector);

  void GetSPDLocalZ(Int_t j, Float_t& z);

  void GetModuleIDs(AliEveITSModuleSelection* sel, std::vector<UInt_t>& ids);

  virtual void Print(Option_t* opt="") const;

  ClassDef(AliEveITSDigitsInfo, 1);
}; // endclass AliEveITSDigitsInfo

#endif
