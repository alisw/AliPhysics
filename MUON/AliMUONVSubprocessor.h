#ifndef ALIMUONVSUBPROCESSOR_H
#define ALIMUONVSUBPROCESSOR_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup shuttle
/// \class AliMUONVSubprocessor
/// \brief Base class for a shuttle sub-task for MUON (either TRK or TRG)
/// 
/// \author Laurent Aphecetche

#ifndef ROOT_TNamed
#  include "TNamed.h"
#endif

class TMap;
class AliMUONPreprocessor;

class AliMUONVSubprocessor : public TNamed
{
public:
  AliMUONVSubprocessor(AliMUONPreprocessor* master,
                       const char* name="", const char* title="");
  virtual ~AliMUONVSubprocessor();
  
  virtual void Initialize(Int_t run, UInt_t startTime, UInt_t endTime);
  virtual UInt_t Process(TMap* dcsAliasMap) = 0;
  
protected:
  AliMUONPreprocessor* Master() const { return fMaster; }
  
  AliMUONVSubprocessor();
  AliMUONVSubprocessor(const AliMUONVSubprocessor&);
  AliMUONVSubprocessor& operator=(const AliMUONVSubprocessor&);
  
private:
  AliMUONPreprocessor* fMaster; ///< Pointer to our master
  
  ClassDef(AliMUONVSubprocessor,1) // Base class of MUON shuttle sub(pre)processors
};

#endif
