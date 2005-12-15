#ifndef ALITOFRAWSTREAM_H
#define ALITOFRAWSTREAM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////////
//
// This class provides the key-reading for TOF raw data.
//
////////////////////////////////////////////////////////////

#include <TObject.h>

class AliRawReader;


class AliTOFRawStream: public TObject {
  public :

  AliTOFRawStream(AliRawReader* rawReader);
  virtual ~AliTOFRawStream();

  virtual Bool_t   Next();
  
  //void  ResetCounter() {fCounter = -1;}; // v0.01

  Int_t GetDDL()       const {return fDDL;};
  Int_t GetTRM()       const {return fTRM;};
  Int_t GetTDC()       const {return fTDC;};
  Int_t GetChannel()   const {return fTDCChannel;};
  
  Int_t GetSector() const;
  Int_t GetPlate()  const;
  Int_t GetStrip()  const;
  Int_t GetPadZ()   const;
  Int_t GetPadX()   const;
  
  Int_t GetTofBin() const {return fTof;};
  Int_t GetADCbin() const {return fADC;};
    
  enum {kDDLOffset = 0x500};      // offset for DDL numbers
    
  private :
 
  AliTOFRawStream(const AliTOFRawStream& stream);
  AliTOFRawStream& operator = (const AliTOFRawStream& stream);

  AliRawReader*    fRawReader;  // object for reading the raw data

  Int_t            fDDL;        // index of current DDL file
  Int_t            fTRM;        // index of current TRM
  Int_t            fTDC;        // index of current TDC
  Int_t            fTDCChannel; // index of current channel of the TDC
  Int_t            fTof;        // time-of-flight measurement
  Int_t            fADC;        // 'charge' measurement
  Int_t            fErrorFlag;  // error flag
  
  AliTOFGeometry *fTOFGeometry; // pointer to the TOF geometry

  //Int_t            fCounter;    // counter for TOF raw data rows in DDL files // v0.01

  ClassDef(AliTOFRawStream, 1)  // class for reading TOF raw digits
};

#endif
