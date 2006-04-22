#ifndef ALIEMCALRAWSTREAM_H
#define ALIEMCALRAWSTREAM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
///
/// This class provides access to EMCAL digits in raw data.
///
///////////////////////////////////////////////////////////////////////////////

// --- ROOT system ---

// --- AliRoot header files ---
#include "AliAltroRawStream.h"
class AliRawReader;


class AliEMCALRawStream: public AliAltroRawStream {

public :
  
  AliEMCALRawStream(AliRawReader* rawReader);
  virtual ~AliEMCALRawStream();
  
  virtual void             Reset();
  virtual Bool_t           Next();

  Int_t            GetId() const {return fId;};
  Int_t            GetPrevId() const {return fPrevId;};
  Int_t            GetModule() const {return fModule;}
  Int_t            GetPrevModule() const {return fPrevModule;}
  Bool_t           IsNewId() const {return (fId != fPrevId);};
  Bool_t           IsNewModule() const {return (fModule != fPrevModule) || (fId != fPrevId);};
  
protected:
  AliEMCALRawStream(const AliEMCALRawStream& stream);
  AliEMCALRawStream& operator = (const AliEMCALRawStream& stream);

  virtual void ApplyAltroMapping();

  Int_t            fId;         //Id of channel
  Int_t            fPrevId;     //previous id
  Int_t            fModule;     //module containing channel
  Int_t            fPrevModule; //previous module

  ClassDef(AliEMCALRawStream, 0)   // class for reading EMCAL raw digits
    };

#endif
