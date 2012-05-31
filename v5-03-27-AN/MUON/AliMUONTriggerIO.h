#ifndef ALIMUONTRIGGERIO_H
#define ALIMUONTRIGGERIO_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup calib
/// \class AliMUONTriggerIO
/// \brief Handles read/write of masks and LUT to/from online files
/// 
//  Author Laurent Aphecetche, Subatech

#ifndef ROOT_TObject
#  include <TObject.h>
#endif

#include "AliMpExMap.h"
#include "AliMpRegionalTrigger.h"



#ifndef ROOT_TArrayI
#  include <TArrayI.h>
#endif

class AliMUONTriggerLut;
class AliMUONVStore;
class AliMpExMap;
class AliMpDDL;
class AliMpTriggerCrate;
class AliMpLocalBoard;
class AliMUONGlobalCrateConfig;
class AliMUONRegionalTriggerConfig;


class AliMUONTriggerIO : public TObject
{
public:
  AliMUONTriggerIO();
  AliMUONTriggerIO(const char* regionalFileToRead);
  virtual ~AliMUONTriggerIO();


  Bool_t ReadConfig(const char* localFile,
                   const char* regionalFile,
                   const char* globalFile,
                   AliMUONVStore* localMasks,
                    AliMUONRegionalTriggerConfig* regionalConfig,
                    AliMUONGlobalCrateConfig* globalConfig);

  Bool_t ReadLUT(const char* lutFileToRead, AliMUONTriggerLut& lut);
  
  Bool_t WriteLUT(const AliMUONTriggerLut& lut,
                  const char* lutFileToWrite);
  
  Bool_t WriteConfig(const char* localFile,
		    const char* regionalFile,
		    const char* globalFile,
		    const AliMUONVStore* localMasks,
                    AliMUONRegionalTriggerConfig* regionalConfig,
                    AliMUONGlobalCrateConfig* globalConfig) const;

  
  Int_t  ReadGlobalConfig(const char* globalFile, AliMUONGlobalCrateConfig* globalConfig) const;

  Bool_t WriteGlobalConfig(const char* globalFile, AliMUONGlobalCrateConfig* globalConfig) const;

  Int_t  ReadRegionalConfig(const char* regionalFile, AliMUONRegionalTriggerConfig* regionalConfig);

  Bool_t WriteRegionalConfig(const char* regionalFile, AliMUONRegionalTriggerConfig* regionalConfig) const;

  Int_t  ReadLocalMasks(const char* localFile, AliMUONVStore& localMasks) const;
  
  Bool_t WriteLocalMasks(const char* localFile, const AliMUONVStore& localMasks) const;

  void   ReadLocalLUT(AliMUONTriggerLut& lut, Int_t localBoardId, FILE* flut);
  
  void   WriteLocalLUT(const AliMUONTriggerLut& lut, Int_t localBoardId, 
                       FILE* flut);
                   
  Int_t LocalBoardId(Int_t index) const;
  Int_t LocalBoardId(Int_t ddlId, Int_t crateId, Int_t localId) const;


private:
  
  Bool_t DeCompAddress(UChar_t &ypos, UChar_t &ytri, UChar_t &xdev, UChar_t &xpos, 
                     UShort_t address) const;
    
  void FillLut(AliMUONTriggerLut& lut,
               Int_t icirc, UChar_t istripX, UChar_t idev,  
               Int_t lutLpt[16][2], Int_t lutHpt[16][2]) ;
  
  
  /// Return number of local boards
  Int_t NofLocalBoards() const { return fRegionalTrigger.GetNofLocalBoards(); }
  
  
  
private:
  AliMpRegionalTrigger  fRegionalTrigger; //!< Regional trigger
 
  static const UInt_t  fgkLocalLutSize;  ///< length of the lut for one local board

  
  ClassDef(AliMUONTriggerIO,2) // Read/Write trigger masks and LUT to/from online files
};

#endif
