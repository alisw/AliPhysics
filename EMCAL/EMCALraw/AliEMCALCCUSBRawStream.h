#ifndef ALIEMCALCCUSBRAWSTREAM_H
#define ALIEMCALCCUSBRAWSTREAM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

///////////////////////////////////////////////////////////////////////////////
///
/// This class provides access to CC-USB data in EMCAL test bench raw data.
/// Author: guernane@lpsc.in2p3.fr
///
///////////////////////////////////////////////////////////////////////////////

#include <TObject.h>

class AliRawReader;

class AliEMCALCCUSBRawStream: public TObject {
  public :
    AliEMCALCCUSBRawStream(AliRawReader* rawReader);
    virtual ~AliEMCALCCUSBRawStream() {};

    virtual Bool_t   Next();

    UInt_t           GetTDC(Int_t iTDC) const
    {return fTDC[iTDC];}
    
    UInt_t           GetQDC(Int_t iQDC) const
    {return fQDC[iQDC];}
    UInt_t           GetScalerCCUSB(Int_t iScaler) const 
    {return fScalerCCUSB[iScaler];}
    UInt_t           GetScalerLecroy(Int_t iScaler) const
    {return fScalerLecroy[iScaler];}

  private :
    AliEMCALCCUSBRawStream(const AliEMCALCCUSBRawStream& stream);
    AliEMCALCCUSBRawStream& operator = (const AliEMCALCCUSBRawStream& stream);

    AliRawReader*    fRawReader;     // object for reading the raw data
    UInt_t           fData;          // data read for file
    UInt_t           fHeader;        // bit 15=1 indicates a watchdog buffer
                                     // bit 14=1 indicates a scaler buffer
                                     // bits 0-9 represent the number of events in the buffer
    UInt_t           fOptHeader;     // bits 0-11 represent the number of words in the buffer
    UInt_t           fEventLength;   // event length including terminator words
    UInt_t           fEOBuffer;      // event terminator
    
    static const Int_t fgkNScalerCCUSB  =  2; // number of internal CC-USB scalers
    static const Int_t fgkNScalerLecroy = 12; // number of Lecroy scalers
    static const Int_t fgkNTDC    = 40;       // number of TDC
    static const Int_t fgkNQDC    = 32;       // number of QDC
    
    UInt_t         fTDC[fgkNTDC];                   // TDC channels
    UInt_t         fQDC[fgkNQDC];                   // QDC values         
    UInt_t         fScalerCCUSB[fgkNScalerCCUSB];   // Internal scaler values         
    UInt_t         fScalerLecroy[fgkNScalerLecroy]; // Lecroy scaler values         

    ClassDef(AliEMCALCCUSBRawStream, 0)  // class for reading CC-USB raw digits
};

#endif
