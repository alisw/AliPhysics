/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

//-----------------------------------------------------------------------
// Utility class to change from online to offline numbering scheme
// for SPD. Most of these methods already exist in AliITSRawStreamSPD.
// Author : A. Mastroserio
//-----------------------------------------------------------------------

#ifndef ALISPDUTILS_H
#define ALISPDUTILS_H

#include "TObject.h"

class AliSPDUtils : public TObject {

 public:
 
  AliSPDUtils(){;}  
   
 
  virtual ~AliSPDUtils();

 
  // module mapping
  static Int_t GetModuleNumber(UInt_t iDDL, UInt_t iModule); 
  static Int_t GetModuleNumber(UInt_t iDDL, UInt_t iHS, UInt_t iChip)  {return GetOfflineModuleFromOnline(iDDL,iHS,iChip);}

  // general coordinate conversions:
  static Bool_t OfflineToOnline(UInt_t module, UInt_t colM, UInt_t RowM, UInt_t& eq, UInt_t& hs, UInt_t& chip, UInt_t& col, UInt_t& row); 
  static Bool_t OnlineToOffline(UInt_t eq, UInt_t hs, UInt_t chip, UInt_t col, UInt_t row, UInt_t& module, UInt_t& colM, UInt_t& rowM);
  static Bool_t GetOfflineFromOfflineChipKey(UInt_t chipkey,UInt_t& module, UInt_t& chip);
 
  // specific coordinate conversions - offline->online
  static UInt_t GetOnlineEqIdFromOffline(UInt_t module);
  static UInt_t GetOnlineHSFromOffline(UInt_t module);
  static UInt_t GetOnlineChipFromOffline(UInt_t module, UInt_t colM);
  static UInt_t GetOnlineColFromOffline(UInt_t module, UInt_t colM); 
  static UInt_t GetOnlineRowFromOffline(UInt_t module, UInt_t rowM);    
  static Bool_t GetOnlineFromOfflineChipKey(UInt_t chipkey,UInt_t& eq, UInt_t& hs, UInt_t& chip);
 
  // specific coordinate conversions - online->offline
  static UInt_t GetOfflineModuleFromOnline(UInt_t eqId, UInt_t hs, UInt_t chip);
  static UInt_t GetOfflineChipKeyFromOnline(UInt_t eqId, UInt_t hs, UInt_t chip); 
  static UInt_t GetOfflineColFromOnline(UInt_t eqId, UInt_t hs, UInt_t chip, UInt_t col);
  static UInt_t GetOfflineRowFromOnline(UInt_t eqId, UInt_t hs, UInt_t chip, UInt_t row);

 
 

  private :

    static const Int_t fgkDDLModuleMap[20][12];  // mapping DDL/module -> module number
  AliSPDUtils& operator= (const AliSPDUtils& c); // dummy
  AliSPDUtils(const AliSPDUtils& c);             // dummy

  ClassDef(AliSPDUtils,0);
};

#endif
