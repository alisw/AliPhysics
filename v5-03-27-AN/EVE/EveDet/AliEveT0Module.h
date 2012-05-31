// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/
#ifndef AliEveT0Module_H
#define AliEveT0Module_H

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// The main AliEVE drawing module for the T0 detector                   //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <TEveQuadSet.h>

class AliT0digit;
class AliRawReader;
class AliT0RawReader;
class TTree;

class AliEveT0Module : public TEveQuadSet
{
public:
  AliEveT0Module(const Text_t* n="AliEveT0Module", Int_t sigType=0,
		 AliT0digit *digits=0, AliT0RawReader *start=0,
		 Double_t zvertex=0);
  virtual ~AliEveT0Module() {}

  virtual void DigitSelected(Int_t idx);

  void PrintEventInfo(); // *MENU*

  static void LoadRaw(AliRawReader* reader);

  static void MakeModules(AliT0digit *digits);

protected:
  Int_t           fSigType; // 0 ~ ADC, 1 ~ TDC
  AliT0digit     *fDigits;  // Digits.
  AliT0RawReader *fStart;   // Reader.

  Double_t        fZVertex; // Reconstructed vertex position.

private:
  AliEveT0Module(const AliEveT0Module&);
  AliEveT0Module& operator=(const AliEveT0Module&);

  ClassDef(AliEveT0Module, 0); // Representation of a T0 module.
};

#endif
