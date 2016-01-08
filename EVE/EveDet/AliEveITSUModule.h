
/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef ALI_EVE_ITS_U_MODULE_H
#define ALI_EVE_ITS_U_MODULE_H

#include <TEveQuadSet.h>
class AliITSUGeomTGeo;
class AliITSMFTDigitPix;

class AliEveITSUModule : public TEveQuadSet
{
public:

  AliEveITSUModule(const Text_t* n="AliEveITSUModule", const Text_t* t=0);
  AliEveITSUModule(AliITSUGeomTGeo *gm,Int_t id, Int_t layer, Int_t ladder, Int_t detector);
  virtual ~AliEveITSUModule();

  static void InitStatics();

  void SetDigitInQuad(AliITSMFTDigitPix *pDig);
  void SetTrans();
  void DigitSelected(Int_t idx);

  virtual void Print(Option_t* opt="") const;

  //  virtual void LoadQuads() {};
  //  virtual void DigitSelected(Int_t idx) {};

  Int_t GetID() const { return fID; }
  void  SetID(Int_t gid, Bool_t trans=kTRUE);

protected:

  static TEveFrameBox    *fgITSUFrameBox;     // Module frame for ITS Upgrade.
  static TEveFrameBox    *fgITSUFrameBoxDead; // Dead-module frame for ITS Upgrade.
  static TEveRGBAPalette *fgITSUPalette;  // Signal to color mapping for ITS Upgrade.

  Int_t       fID;      // Module id.
  const Int_t fkLayer;  // which layer
  Int_t fkLadder;       // which ladder
  Int_t fkDetector;     // which detector (module within ladder)

  Float_t fDpx;     // Digit size in x.
  Float_t fDpz;     // Digit size in z.

  Bool_t fAtLeastOneDigit;   // is there already a digit put into the geometry?

  static Bool_t fgStaticInitDone; // Flag for static variable initialization.

private:
  AliEveITSUModule(const AliEveITSUModule&);            // Not implemented
  AliEveITSUModule& operator=(const AliEveITSUModule&); // Not implemented

  ClassDef(AliEveITSUModule, 0); // Visualization of an ITS Upgrade module.
};

#endif
