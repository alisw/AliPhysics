// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/
#ifndef ALIEVE_T0Module_H
#define ALIEVE_T0Module_H

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// The main AliEVE drawing module for the T0 detector                   //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#include <TEveQuadSet.h>
#include <AliT0digit.h>
#include <AliT0RawReader.h>


class AliEveT0Module : public TEveQuadSet
{

  AliEveT0Module(const AliEveT0Module&);
  AliEveT0Module& operator=(const AliEveT0Module&);

public:

  AliEveT0Module(const Text_t* n="AliEveT0Module", Int_t sigType=0, AliT0digit *digits=0,AliT0RawReader *start=0);
  virtual ~AliEveT0Module();

  virtual void DigitSelected(Int_t idx);

  void LoadRaw(TString fileName, Int_t ievt);

  static void MakeModules(AliT0digit *digits);

protected:
  Int_t           fSigType; // 0 ~ ADC, 1 ~ TDC
  AliT0digit     *fDigits;
  AliT0RawReader *fStart;

   ClassDef(AliEveT0Module,1);
};

/*
 class T0ModuleTDC : public AliEveT0Module
 {
 public:
   // constructor

    virtual void QuadSelected(Int_t idx);
 };
*/

#endif
