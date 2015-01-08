#ifndef ALIMUONRAWSTREAM_H
#define ALIMUONRAWSTREAM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup raw
/// \class AliMUONRawStream
/// \brief Base class for reading MUON raw digits
///
//  Author: Christian Finck

#include <TObject.h>

class AliRawReader;

class AliMUONRawStream: public TObject {
  public :
    AliMUONRawStream();
    AliMUONRawStream(AliRawReader* rawReader);
    virtual ~AliMUONRawStream();

    /// Initialize iterator
    virtual void First() {return;} // not yet virtual pure, waiting for trigger
    
    /// DDL iterator 
    virtual Bool_t NextDDL() = 0;

    /// Whether the iteration is finished or not
    virtual Bool_t IsDone() const {return kTRUE;} // not yet virtual pure, waiting for trigger

    /// add error message into error logger
    virtual void AddErrorMessage() = 0;

    /// Set object for reading the raw data
    virtual void SetReader(AliRawReader* rawReader) {fRawReader = rawReader;}

    /// Get object for reading the raw data
    virtual AliRawReader* GetReader() {return fRawReader;}

    /// Enable error info logger
    virtual void EnabbleErrorLogger() {fEnableErrorLogger = kTRUE;}

    /// Check if error info logger enable
    virtual Bool_t IsErrorLogger() const {return fEnableErrorLogger;}

    /// swap method for Power PC
    virtual void Swap(UInt_t *buffer, Int_t size) const;


  private :
    /// Not implemented
    AliMUONRawStream(const AliMUONRawStream& stream);
    /// Not implemented
    AliMUONRawStream& operator = (const AliMUONRawStream& stream);

    typedef struct {
     UInt_t fB1:8; ///< first byte word
     UInt_t fB2:8; ///< second byte word
     UInt_t fB3:8; ///< third byte word
     UInt_t fB4:8; ///< fourth byte word
    } RawWord;

    AliRawReader* fRawReader;    //!< object for reading the raw data  
    Bool_t fEnableErrorLogger;   //!< flag to enable the error info logger
    
    ClassDef(AliMUONRawStream, 1)    // base class for reading MUON raw digits
};

#endif
