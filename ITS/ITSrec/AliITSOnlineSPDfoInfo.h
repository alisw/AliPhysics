#ifndef ALIITSONLINESPDFOINFO_H
#define ALIITSONLINESPDFOINFO_H
/* Copyright(c) 2008-2010, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////////////
// Author: A. Mastroserio                                     // 
// This class is used within the detector algorithm framework //
// to collect information on how the scan was arranged.       //
////////////////////////////////////////////////////////////////


#include <TObject.h>
#include <TArrayS.h>
#include <TBits.h>

class AliITSOnlineSPDfoInfo :  public TObject {

 public:
  AliITSOnlineSPDfoInfo();
  virtual ~AliITSOnlineSPDfoInfo();

  virtual void   ClearThis();
  virtual void   AddDACindex(Short_t index);

  // SETTERS
  virtual void SetRunNumber(UInt_t val)  {fRunNumber=val;}
  virtual void SetRouter(UShort_t val)   {fRouter=val;}
  virtual void SetNumTriggers(UInt_t val){fNumTriggers=val;}
  virtual void SetDBversion(Int_t val)   {fDBversion=val;}
  void SetActiveChipsAndHS(UInt_t hs, UInt_t chip) {fActiveChipsAndHS.SetBitNumber(10*hs+chip);}

  // GETTERS
  UInt_t   GetRunNumber() const     {return fRunNumber;}
  UShort_t GetRouter() const        {return fRouter;}
  UInt_t   GetNumTriggers() const   {return fNumTriggers;}
  Int_t    GetDBversion() const     {return fDBversion;}

  UShort_t GetNumDACindex() const   {return fNumDACindex;}
  Short_t  GetDACindex(UShort_t id) const; // returns -1 if ID not present

  TArrayS GetDACIndexArray() const  {return fDACindex;}

  Bool_t IsActiveHS(UInt_t hs) const ;
  Bool_t IsActiveChip(UInt_t hs, UInt_t chip) const;
  TBits  GetActiveChipsAndHS() const {return fActiveChipsAndHS;}

 protected:
  UInt_t   fRunNumber;   // run number
  UShort_t fRouter;      // router number (should be same as eq number)
  UInt_t   fNumTriggers; // number of triggers sent for each scan step
  Int_t    fDBversion;   // global configuration db version
			 
  UShort_t fNumDACindex; // number of DAC indices in TArrayI below
  TArrayS  fDACindex;    // list of DAC indices related to each DAC value
  TBits    fActiveChipsAndHS;

  ClassDef(AliITSOnlineSPDfoInfo,2)
};
    
#endif
