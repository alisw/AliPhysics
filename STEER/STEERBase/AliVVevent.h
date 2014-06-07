#ifndef ALIVVEVENT_H
#define ALIVVEVENT_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               *
 * Primary Authors : Sergey Gorbunov, Jochen Thaeder, Chiara Zampolli     */

/*
 * See implementation file for documentation
 */
#include "Rtypes.h"
#include "TString.h"
class AliVVevent;
class AliVVvertex;
class AliVVtrack;
class AliMultiplicity;
class AliVVkink;
class AliVVeventFriend;

class AliVVevent {
 public:
  // --------------------------------------------------------------------------------
  // -- Constructor / Destructors
  AliVVevent() {}
  virtual ~AliVVevent() {} 

  // --------------------------------------------------------------------------------
  virtual void Reset() = 0;

  // --------------------------------------------------------------------------------
  // Access methods

  virtual const AliVVvertex* GetPrimaryVertex() const {return NULL;}
  virtual const AliVVvertex* GetPrimaryVertexSPD() const {return NULL;}
  virtual const AliVVvertex* GetPrimaryVertexTracks() const {return NULL;}
  virtual const AliVVvertex* GetPrimaryVertexTPC() const {return NULL;}
  virtual AliVVtrack* GetTrack(Int_t /*i*/) const {return NULL;}
  virtual AliVVkink* GetKink(Int_t /*i*/) const {return NULL;}
  virtual AliVVtrack* GetV0(Int_t /*i*/) const {return 0;}
  virtual Int_t GetNumberOfTracks() const {return 0;}
  virtual Int_t GetNumberOfV0s() const {return 0;}
  virtual Int_t GetNumberOfKinks() const {return 0;}
  virtual Int_t GetEventNumberInFile() const {return -1;}
  virtual const AliMultiplicity* GetMultiplicity() const {return NULL;} //by default SPDmult
  virtual Int_t GetRunNumber() const {return -1;}
  virtual TString GetFiredTriggerClasses() const {TString string; return string;}
  virtual TObject* FindListObject(const char* /*name*/) const {return NULL;}
  virtual ULong64_t GetTriggerMask() const {return 0;}
  virtual Double_t GetMagneticField() const {return 0;}
  virtual UInt_t GetTimeStamp() const { return 0;}
  virtual UInt_t GetEventSpecie() const { return 0;}
  //virtual AliVVeventFriend* FindFriend() const { return NULL; }

  //  ClassDef(AliVVevent,1)  // base class for AliEvent data

};
#endif
