#ifndef ALIMUONTRIGGERIO_H
#define ALIMUONTRIGGERIO_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup trigger
/// \class AliMUONTriggerIO
/// \brief Handles read/write of masks and LUT to/from online files
/// 
//  Author Laurent Aphecetche, Subatech

#ifndef ROOT_TObject
#  include <TObject.h>
#endif

#include "AliMpExMap.h"
#include "AliMpGlobalCrate.h"
#include "AliMpRegionalTrigger.h"

#ifndef ROOT_TArrayI
#  include <TArrayI.h>
#endif

class AliMUONTriggerLut;
class AliMUONVCalibParam;
class AliMUONVStore;
class AliMpExMap;
class AliMpDDL;
class AliMpTriggerCrate;
class AliMpLocalBoard;

class AliMUONTriggerIO : public TObject
{
public:
  AliMUONTriggerIO();
  AliMUONTriggerIO(const char* regionalFileToRead);
  virtual ~AliMUONTriggerIO();

  Bool_t ReadMasks(const char* localFile,
                   const char* regionalFile,
                   const char* globalFile,
                   AliMUONVStore* localMasks,
                   AliMUONVStore* regionalMasks,
                   AliMUONVCalibParam* globalMasks,
		   Bool_t warn = true);
  
  Bool_t ReadLUT(const char* lutFileToRead, AliMUONTriggerLut& lut);
  
  Bool_t WriteLUT(const AliMUONTriggerLut& lut,
                  const char* lutFileToWrite);
  
  Bool_t WriteMasks(const char* localFile,
		    const char* regionalFile,
		    const char* globalFile,
		    AliMUONVStore* localMasks,
                    AliMUONVStore* regionalMasks,
                    AliMUONVCalibParam* globalMasks) const;
  
  Int_t LocalBoardId(Int_t index) const;

  void UpdateMapping(Bool_t writeFile = true) const;

private:
  
  Bool_t DeCompAddress(UChar_t &ypos, UChar_t &ytri, UChar_t &xdev, UChar_t &xpos, 
                     UShort_t address) const;
    
  void FillLut(AliMUONTriggerLut& lut,
               Int_t icirc, UChar_t istripX, UChar_t idev,  
               Int_t lutLpt[16][2], Int_t lutHpt[16][2]) ;
  
  
  /// Return number of local boards
  Int_t NofLocalBoards() const { return fRegionalTrigger.GetNofLocalBoards(); }
  
  Int_t  ReadGlobal(const char* globalFile, AliMUONVCalibParam* globalMasks);

  Bool_t WriteGlobal(const char* globalFile, AliMUONVCalibParam* globalMasks) const;

  Int_t  ReadRegional(const char* regionalFile, AliMUONVStore* regionalMasks, Bool_t warn = true);

  Bool_t WriteRegional(const char* regionalFile, AliMUONVStore* regionalMasks) const;

  Int_t  ReadLocalMasks(const char* localFile, AliMUONVStore& localMasks) const;
  
  Bool_t WriteLocalMasks(const char* localFile, AliMUONVStore& localMasks) const;

  void   ReadLocalLUT(AliMUONTriggerLut& lut, Int_t localBoardId, FILE* flut);
  
  void   WriteLocalLUT(const AliMUONTriggerLut& lut, Int_t localBoardId, 
                       FILE* flut);
    
  
private:
  AliMpRegionalTrigger  fRegionalTrigger; //!< Regional trigger
  AliMpGlobalCrate      fGlobalCrate;     //!< Global crate object
 
  ClassDef(AliMUONTriggerIO,0) // Read/Write trigger masks and LUT to/from online files
};

#endif
