////////////////////////////////////////////////
//  Digitization class for set: TOF           //
//  AliTOFRawSector class                     //
//  Interface                                 //
//  Description                               //
//*-- Authors: Pierella, Seganti, Vicinanza   //
//    (Bologna and Salerno University)        //
////////////////////////////////////////////////

#ifndef ALITOFRAWSECTOR_H
#define ALITOFRAWSECTOR_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


#include "TObject.h"
#include "TClonesArray.h"

//_______________________________________________________
class AliTOFRawSector : public TObject{

 public:
  AliTOFRawSector();
// dtor
  virtual ~AliTOFRawSector();
// copy ctor  (required also by RC10 Coding Convention)
  AliTOFRawSector(const AliTOFRawSector& tofrawsector);
// assignment operator (required also by RC10 Coding Convention)
  AliTOFRawSector& operator = (const AliTOFRawSector& tofrawsector);
  void   WriteSector(); // write a DAQ sector
  void   ReadSector();  // read  a DAQ sector

// getters for AliTOFRawSector object  
  TClonesArray* GetRocData()        const {return fRocData;}
  UInt_t        GetHeader()         const {return fHeader;}
  UInt_t        GetGlobalCheckSum() const {return fGlobalCheckSum;}

// setters for AliTOFRawSector object
  void SetGlobalCS(UInt_t gcs){fGlobalCheckSum=gcs;}
  void SetHeader  (UInt_t hdr){fHeader = hdr;}
  
 protected:
  TClonesArray* fRocData; // pointer to the TClonesArray of Roc Data
  UInt_t        fHeader;    // RawSector header number   
  UInt_t        fGlobalCheckSum; // check flag

  ClassDef(AliTOFRawSector,2) // Container Class for AliTOFRoc objects
};

#endif /* ALITOFRAWSECTOR_H */
