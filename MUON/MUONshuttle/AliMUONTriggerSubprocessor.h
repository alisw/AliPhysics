#ifndef ALIMUONTRIGGERSUBPROCESSOR_H
#define ALIMUONTRIGGERSUBPROCESSOR_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup shuttle
/// \class AliMUONTriggerSubprocessor
/// \brief Implementation of AliMUONVSubprocessor for MUON TRK masks
/// 
//  Author Laurent Aphecetche, Subatech

#ifndef ALIMUONVSUBPROCESSOR_H
#  include "AliMUONVSubprocessor.h"
#endif

class AliMUONTriggerLut;
class AliMUONRegionalTriggerConfig;
class AliMUONGlobalCrateConfig;
class AliMUONVCalibParam;
class AliMUONVStore;
class TString;
class TClonesArray;

class AliMUONTriggerSubprocessor : public AliMUONVSubprocessor
{
public:
  AliMUONTriggerSubprocessor(AliMUONPreprocessor* master);
  virtual ~AliMUONTriggerSubprocessor();
  
  Bool_t Initialize(Int_t run, UInt_t startTime, UInt_t endTime);
  UInt_t Process(TMap* dcsAliasMap);
  
private:

  TString GetFileName(const char* fid) const;
  
  /// Not implemented
  AliMUONTriggerSubprocessor(const AliMUONTriggerSubprocessor&);
  /// Not implemented
  AliMUONTriggerSubprocessor& operator=(const AliMUONTriggerSubprocessor&);
  
  Int_t TestFile(const char* baseName, Bool_t shouldBeThere) const;

  void WhichFilesToRead(const char* exportedFiles,
                        Bool_t& globalFile,
                        Bool_t& regionalFile,
                        Bool_t& localFile,
                        Bool_t& lutFile,
			Bool_t& trigScalFile);
  
private:
  AliMUONRegionalTriggerConfig* fRegionalConfig; //!<! regional config
  AliMUONVStore* fLocalMasks; //!<! local masks
  AliMUONGlobalCrateConfig*     fGlobalConfig;   //!<! global config
  AliMUONTriggerLut* fLUT; //!<! look-up table(s)
  TClonesArray* fTrigScalers; //!<! trigger scalers
  Bool_t fRegionalConfigToOCDB; //!<! to store in the OCDB
  Bool_t fLocalMasksToOCDB; //!<! to store in the OCDB
  Bool_t fGlobalConfigToOCDB; //!<! to store in the OCDB
  Bool_t fLUTToOCDB; //!<! to store in the OCDB
  Bool_t fTrigScalersToOCDB; //!<! to store in the OCDB
  
  ClassDef(AliMUONTriggerSubprocessor,4) // A shuttle preprocessor for MUON TRK masks
};

#endif
