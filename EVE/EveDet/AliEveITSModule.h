// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/
#ifndef ALIEVE_ITSModule_H
#define ALIEVE_ITSModule_H

#include <TEveQuadSet.h>

#include <EveDet/AliEveITSDigitsInfo.h>


class AliEveITSModule : public TEveQuadSet
{
  AliEveITSModule(const AliEveITSModule&);            // Not implemented
  AliEveITSModule& operator=(const AliEveITSModule&); // Not implemented

protected:
  AliEveITSDigitsInfo* fInfo;

  Int_t       fID;    // Module   id
  Int_t       fDetID; // Detector id (0~SPD, 1~SDD, 2~SSD)

  Int_t       fLayer;
  Int_t       fLadder;
  Int_t       fDet;

  Float_t     fDx;
  Float_t     fDz;
  Float_t     fDy;

  static Bool_t fgStaticInitDone;

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

  static TEveFrameBox* fgSPDFrameBox;
  static TEveFrameBox* fgSDDFrameBox;
  static TEveFrameBox* fgSSDFrameBox;

  static TEveRGBAPalette* fgSPDPalette;
  static TEveRGBAPalette* fgSDDPalette;
  static TEveRGBAPalette* fgSSDPalette;

  ClassDef(AliEveITSModule, 1);
};

#endif
