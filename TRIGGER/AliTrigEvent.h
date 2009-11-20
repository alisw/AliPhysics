#ifndef ALITRIGEVENT_H
#define ALITRIGEVENT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
// Author: Andrei Gheata, 28/07/2009

//==============================================================================
//   AliTrigEvent - Base class for generic information exchanged by a trigger
//                  device. Trigger inputs and outputs are represented and
//                  handled via AliTrigEvent objects. Trigger events are typically
//                  wrappers for the information exchanged on a single I/O slot
//                  or a group of correlated inputs.
//==============================================================================

#ifndef ROOT_TObject
#include "TObject.h"
#endif

class TClass;
class AliTrigEvent : public TObject {

public:
enum ETrigSignalFlags {
  kActive = BIT(14)
}
  
  AliTrigEvent() : TObject {}
  virtual ~AliTrigEvent() {}

  void                      Activate(Bool_t flag);
  Bool_t                    IsActive() const {return TObject::TestBit(kActive);}
  
  // Import data from a source signal. Has to be implemented by derived signals.
  virtual Bool_t            ImportData(AliTrigEvent *source) = 0;
  virtual Bool_t            SetType(const char *classname) {return kTRUE;}
  virtual TClass           *GetType() {return NULL;}
     
  ClassDef(AliTrigEvent,1)  // Generic event embedding data.
};


//==============================================================================
//   AliTrigEventWithMask - Event embedding a bit mask as TBits
// 
//==============================================================================
class TBits;
class AliTrigEventWithMask : public AliTrigEvent {

public:
  AliTrigEventWithMask() : AliTrigEvent(), fValue() {}
  virtual ~AliTrigEventWithMask() {}
  
  virtual Bool_t            ImportData(AliTrigEvent *source);

  const TBits              &GetValue() const {return fValue;}
  void                      SetValue(const TBits &value) {fValue = value;}

private:
  TBits                     fValue;  // Mask value
     
  ClassDef(AliTrigEventWithMask,1)  // Signal embedding a TBits object.
};

//==============================================================================
//   AliTrigEventWithObject - Event embedding a TObject
// 
//==============================================================================

class TClass;
class AliTrigEventWithObject : public AliTrigEvent {

public:
  AliTrigEventWithObject() : AliTrigEvent(), fValue(0) {}
  AliTrigEventWithObject(const char *classname);
  AliTrigEventWithObject(const AliTrigEventWithObject &other);
  virtual ~AliTrigEventWithObject() {}
  
  AliTrigEventWithObject   &operator=(const AliTrigEventWithObject &other);
  virtual Bool_t            ImportData(AliTrigEvent *source);
  virtual TClass           *GetType() const  {return fType;}
  virtual Bool_t            SetType(const char *classname);
  TObject                  *GetValue() const {return fValue;}
  Bool_t                    SetValue(TObject *value);

private:
  TObject                  *fValue;  // Embedded object
  TClass                   *fType;   //! Object type
     
  ClassDef(AliTrigEventWithObject,1)  // Signal embedding an object.
};

#endif
