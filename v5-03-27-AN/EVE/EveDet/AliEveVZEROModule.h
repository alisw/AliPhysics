/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/
#ifndef AliEveVZEROModule_H
#define AliEveVZEROModule_H

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// The drawing module for the VZERO detector                            //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <TEveQuadSet.h>

class AliRawReader;
class AliVZERORawStream;
class AliESDEvent;

class AliEveVZEROModule : public TEveQuadSet
{
public:
  AliEveVZEROModule(const Text_t* n="AliEveVZEROModule", Bool_t side = kTRUE);
  virtual ~AliEveVZEROModule();

  virtual void DigitSelected(Int_t idx);

  void LoadRaw(AliRawReader *rawReader);

  Int_t GetSampleIndex() const { return fSampleIndex; }
  void  SetSampleIndex(Int_t index);

protected:

  AliVZERORawStream *fStream;      // Raw-stream
  Int_t              fSampleIndex; // Current sample index used
  Bool_t             fIsASide;     // A or C side module

private:
  AliEveVZEROModule(const AliEveVZEROModule&);
  AliEveVZEROModule& operator=(const AliEveVZEROModule&);

  ClassDef(AliEveVZEROModule,0) // Representation of a VZERO module
};

#endif
