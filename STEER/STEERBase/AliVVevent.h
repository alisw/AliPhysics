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
class TList;
//class AliVVvertex;
class AliVVtrack;
class AliMultiplicity;
class AliVVfriendEvent;
class AliESDkink;
class TTree;

class AliVVevent{
 public:  
  // --------------------------------------------------------------------------------
  // -- Constructor / Destructor
  // --------------------------------------------------------------------------------

  AliVVevent() {}

 protected:

  virtual ~AliVVevent() {}

 public:
  // --------------------------------------------------------------------------------
  // -- Methods, common for AliVEvent and AliVVevent interface
  // --------------------------------------------------------------------------------

  // Access methods
  
  virtual Double_t  GetMagneticField() const = 0;  
  virtual UInt_t    GetPeriodNumber()  const = 0;
  virtual Int_t     GetRunNumber()     const = 0;
  virtual UInt_t    GetOrbitNumber()   const = 0;
  virtual UShort_t  GetBunchCrossNumber() const = 0;
  virtual ULong64_t GetTriggerMask()   const = 0;
  virtual TString   GetFiredTriggerClasses() const = 0;
  
  virtual Int_t GetNumberOfTracks() const = 0;
  virtual Int_t GetNumberOfV0s() const = 0;
 
  // Other methods

  virtual void Reset() = 0;
  virtual void ConnectTracks() = 0;
  
  // --------------------------------------------------------------------------------
  // -- Methods, specific for AliVVevent (are disabled in AliVEvent interface)
  // --------------------------------------------------------------------------------
  
  virtual UInt_t    GetTimeStamp()   const = 0;
  virtual UInt_t    GetEventSpecie() const = 0;  
  virtual ULong64_t GetTriggerMaskNext50()   const = 0;

  virtual Int_t GetNumberOfKinks() const = 0;
  
  /* 
  virtual const AliVVvertex* GetPrimaryVertex() const {return NULL;}
  virtual const AliVVvertex* GetPrimaryVertexSPD() const {return NULL;}
  virtual const AliVVvertex* GetPrimaryVertexTracks() const {return NULL;}
  virtual const AliVVvertex* GetPrimaryVertexTPC() const {return NULL;}  
  */
  
  virtual const AliVVtrack* GetVVTrack(Int_t /*i*/) const = 0;
  virtual AliESDkink* GetKink(Int_t /*i*/) const = 0;
  //virtual AliVVtrack* GetV0(Int_t /*i*/) const = 0;

 
  // some ESD interfaces, not yet implemented in flat esd (may be not needed, may be need some changes)

  //virtual Int_t GetEventNumberInFile() const = 0;
  //virtual const AliMultiplicity* GetMultiplicity() const = 0; //by default SPDmult  
  //virtual AliVVfriendEvent* FindFriend() const = 0;

  //ClassDef(AliVVevent, 1)

};
#endif
