// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/
#ifndef AliEveITSModule_H
#define AliEveITSModule_H

#include <TEveQuadSet.h>

#include <AliEveITSDigitsInfo.h>


class AliEveITSModule : public TEveQuadSet
{
public:
  AliEveITSModule(const Text_t* n="AliEveITSModule", const Text_t* t=0);
  AliEveITSModule(Int_t gid, AliEveITSDigitsInfo* info);
  virtual ~AliEveITSModule();

  static void InitStatics(AliEveITSDigitsInfo* info);

  AliEveITSDigitsInfo* GetDigitsInfo() const { return fInfo; }
  void SetDigitsInfo(AliEveITSDigitsInfo* info);

  Int_t GetSubDetID() const { return fDetID; }

  Int_t GetID() const { return fID; }
  void  SetID(Int_t gid, Bool_t tran=kTRUE);

  virtual void LoadQuads();
  void SetTrans();

  virtual void DigitSelected(Int_t idx);

  virtual void Print(Option_t* opt="") const;

  static TEveFrameBox    *fgSPDFrameBox;     // Module frame for SPD.
  static TEveFrameBox    *fgSPDFrameBoxDead; // Dead-module frame for SPD.
  static TEveFrameBox    *fgSDDFrameBox;     // Module frame for SDD.
  static TEveFrameBox    *fgSDDFrameBoxDead; // Dead-module frame for SPD.
  static TEveFrameBox    *fgSSDFrameBox;     // Module frame for SSD.
  static TEveFrameBox    *fgSSDFrameBoxDead; // Dead-module frame for SPD.

  static TEveRGBAPalette *fgSPDPalette;  // Signal to color mapping for SPD.
  static TEveRGBAPalette *fgSDDPalette;  // Signal to color mapping for SDD.
  static TEveRGBAPalette *fgSSDPalette;  // Signal to color mapping for SSD.

protected:
  AliEveITSDigitsInfo* fInfo; // Source of geometry and data.

  Int_t       fID;      // Module id.
  Int_t       fDetID;   // Detector id (0~SPD, 1~SDD, 2~SSD).

  Int_t       fLayer;   // Layer (0 - 5).
  Int_t       fLadder;  // Ladder.
  Int_t       fDet;     // Detector.

  Float_t     fDx;      // Digit half-size in x.
  Float_t     fDz;      // Digit half-size in z.
  Float_t     fDy;      // Digit half-size in y.

  static Bool_t fgStaticInitDone; // Flag for static variable initialization.

private:
  AliEveITSModule(const AliEveITSModule&);            // Not implemented
  AliEveITSModule& operator=(const AliEveITSModule&); // Not implemented

  ClassDef(AliEveITSModule, 0); // Visualization of an ITS module.
};

#endif
