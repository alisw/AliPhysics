#ifndef ALITPCCALPADREGION_H
#define ALITPCCALPADREGION_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TNamed.h>
#include <TString.h>
#include <TObjArray.h>
#include <TIterator.h>

class AliTPCCalPadRegion: public TNamed {
public:
   AliTPCCalPadRegion();
   AliTPCCalPadRegion(const AliTPCCalPadRegion& obj);
   AliTPCCalPadRegion(const char *name, const char *title);
   //AliTPCCalPadRegion(const TString &name, const TString &title) : TNamed(name, title) { }
   virtual ~AliTPCCalPadRegion() { delete fObjects; }
   AliTPCCalPadRegion& operator=(const AliTPCCalPadRegion& rhs);
   
  virtual TObject*   GetObject(UInt_t segment, UInt_t padType);
  virtual void       SetObject(TObject* obj, UInt_t segment, UInt_t padType);
  virtual void       Delete(Option_t* option = "") { if (fObjects) fObjects->Delete(option); }
   virtual TIterator* MakeIterator(Bool_t direction = kIterForward) const { return fObjects->MakeIterator(direction); }
   static  UInt_t     GetNSegments() { return fgkNSegments; }
   static  UInt_t     GetNPadTypes() { return fgkNPadTypes; }
   static  void       GetPadRegionCenterLocal(UInt_t padType, Double_t* xy);
//   static  UInt_t     GetStartRow(UInt_t padType);
//   static  UInt_t     GetEndRow(UInt_t padType);
    
protected:
   virtual Bool_t BoundsOk(const char* where, UInt_t segment, UInt_t padType) const
      { return (segment >= fgkNSegments || padType >= fgkNPadTypes) ? OutOfBoundsError(where, segment, padType) : kTRUE; }
   virtual Bool_t OutOfBoundsError(const char* where, UInt_t segment, UInt_t padType) const
      { Error(where, "Index out of bounds (trying to access segment %d, pad type %d).", segment, padType); return kFALSE; }

   TObjArray* fObjects;     // array containing an object for each pad region

   static const UInt_t fgkNSegments = 36;    // number of TPC sectors, 0-17: A side, 18-35: C side (IROC and OROC are treated as one sector)
   static const UInt_t fgkNPadTypes = 3;     // number of pad types, 0: short pads, 1: medium pads, 2: long pads

   ClassDef(AliTPCCalPadRegion, 1)
};


#endif
