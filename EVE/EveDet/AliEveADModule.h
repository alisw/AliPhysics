/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/
#ifndef AliEveADModule_H
#define AliEveADModule_H

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// The drawing module for the AD detector                               //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <TEveQuadSet.h>
#include "AliADConst.h"

class AliRawReader;
class AliADRawStream;
class AliESDEvent;

class AliEveADModule : public TEveQuadSet
{
public:
  AliEveADModule(const Text_t* n="AliEveADModule", Bool_t side = kTRUE, Int_t maxCharge = 1023);
  virtual ~AliEveADModule();

  void LoadRaw(AliRawReader *rawReader);

protected:

  AliADRawStream *fStream;      // Raw-stream
  Bool_t          fIsASide;     // A or C side module

private:
  AliEveADModule(const AliEveADModule&);
  AliEveADModule& operator=(const AliEveADModule&);

  ClassDef(AliEveADModule,0) // Representation of a AD module
};

#endif
