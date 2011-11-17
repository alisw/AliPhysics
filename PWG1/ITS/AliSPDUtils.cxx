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
// Author : A. Mastroserio
//-----------------------------------------------------------------------


#ifndef ALISPDUTILS_CXX
#define ALISPDUTILS_CXX

#include <TString.h>
#include "AliSPDUtils.h"
ClassImp(AliSPDUtils)
  //______________________________________________________________________________
  const Int_t AliSPDUtils::fgkDDLModuleMap[20][12] = {
    { 4, 5, 0, 1, 80, 81, 84, 85, 88, 89, 92, 93},
    {12,13, 8, 9, 96, 97,100,101,104,105,108,109},
    {20,21,16,17,112,113,116,117,120,121,124,125},
    {28,29,24,25,128,129,132,133,136,137,140,141},
    {36,37,32,33,144,145,148,149,152,153,156,157},
    {44,45,40,41,160,161,164,165,168,169,172,173},
    {52,53,48,49,176,177,180,181,184,185,188,189},
    {60,61,56,57,192,193,196,197,200,201,204,205},
    {68,69,64,65,208,209,212,213,216,217,220,221},
    {76,77,72,73,224,225,228,229,232,233,236,237},
    { 7, 6, 3, 2, 83, 82, 87, 86, 91, 90, 95, 94},
    {15,14,11,10, 99, 98,103,102,107,106,111,110},
    {23,22,19,18,115,114,119,118,123,122,127,126},
    {31,30,27,26,131,130,135,134,139,138,143,142},
    {39,38,35,34,147,146,151,150,155,154,159,158},
    {47,46,43,42,163,162,167,166,171,170,175,174},
    {55,54,51,50,179,178,183,182,187,186,191,190},
    {63,62,59,58,195,194,199,198,203,202,207,206},
    {71,70,67,66,211,210,215,214,219,218,223,222},
    {79,78,75,74,227,226,231,230,235,234,239,238}
  };
//___________________________________________________________________________
AliSPDUtils::~AliSPDUtils() {
  //
  //destructor
  //
}
//__________________________________________________________________________
Bool_t AliSPDUtils::OfflineToOnline(UInt_t module, UInt_t colM, UInt_t rowM, UInt_t& eq, UInt_t& hs, UInt_t& chip, UInt_t& col, UInt_t& row) {
  // converts offline coordinates to online
  eq = GetOnlineEqIdFromOffline(module);
  hs = GetOnlineHSFromOffline(module);
  chip = GetOnlineChipFromOffline(module,colM);
  col = GetOnlineColFromOffline(module,colM);
  row = GetOnlineRowFromOffline(module,rowM);
  if (eq>=20 || hs>=6 || chip>=10 || col>=32 || row>=256) return kFALSE;
  else return kTRUE;
}
//__________________________________________________________________________
Bool_t AliSPDUtils::OnlineToOffline(UInt_t eq, UInt_t hs, UInt_t chip, UInt_t col, UInt_t row, UInt_t& module, UInt_t& colM, UInt_t& rowM) {
  // converts online coordinates to offline
  module = GetOfflineModuleFromOnline(eq,hs,chip);
  colM = GetOfflineColFromOnline(eq,hs,chip,col);
  rowM = GetOfflineRowFromOnline(eq,hs,chip,row);
  if (module>=240 || colM>=160 || rowM>=256) return kFALSE;
  else return kTRUE;
}
//__________________________________________________________________________
Bool_t AliSPDUtils::GetOfflineFromOfflineChipKey(UInt_t chipkey,UInt_t& module, UInt_t& chip){
  // converts offline chip key to offline chip coordinates (V. Altini)
  if (chipkey>=1200) {
    TString errMess = Form("%d is not a valid Chip Key number",chipkey);
    printf(" ERROR : %s\n", errMess.Data());
    return 0;
  }

  module = chipkey/5;
  chip=chipkey%20%5;

  return 1;
}
//________________________________________________________________________
UInt_t AliSPDUtils::GetOfflineModuleFromOnline(UInt_t eqId, UInt_t hs, UInt_t chip) {
  // online->offline (module)
  if (eqId<20 && hs<6 && chip<10) return fgkDDLModuleMap[eqId][hs*2+chip/5];
  else return 240;
}
//________________________________________________________________________
UInt_t AliSPDUtils::GetOfflineChipKeyFromOnline(UInt_t eqId, UInt_t hs, UInt_t chip) {
  // online->offline (chip key: 0-1199)
  if (eqId<20 && hs<6 && chip<10) {
    UInt_t module = GetOfflineModuleFromOnline(eqId,hs,chip);
    UInt_t chipInModule = ( chip>4 ? chip-5 : chip );
    if(eqId>9) chipInModule = 4 - chipInModule;  // side C only
    return (module*5 + chipInModule);
  } else return 1200;
}
//__________________________________________________________________________
UInt_t AliSPDUtils::GetOnlineEqIdFromOffline(UInt_t module) {
  // offline->online (eq)
  for (UInt_t eqId=0; eqId<20; eqId++) {
    for (UInt_t iModule=0; iModule<12; iModule++) {
      if (GetModuleNumber(eqId,iModule)==(Int_t)module) return eqId;
    }
  }
  return 20; // error
}
//__________________________________________________________________________
UInt_t AliSPDUtils::GetOnlineHSFromOffline(UInt_t module) {
  // offline->online (hs)
  for (UInt_t eqId=0; eqId<20; eqId++) {
    for (UInt_t iModule=0; iModule<12; iModule++) {
      if (GetModuleNumber(eqId,iModule)==(Int_t)module) return iModule/2;
    }
  }
  return 6; // error
}
//__________________________________________________________________________
UInt_t AliSPDUtils::GetOnlineChipFromOffline(UInt_t module, UInt_t colM) {
  // offline->online (chip)
  for (UInt_t eq=0; eq<20; eq++) {
    for (UInt_t iModule=0; iModule<12; iModule++) {
      if (GetModuleNumber(eq,iModule)==(Int_t)module) {
	if (module<80) {
	  if (eq<10) { // side A
	    return (159-colM)/32 + 5*(iModule%2);
	  }
	  else { // side C
	    return colM/32 + 5*(iModule%2);
	  }
	}
	else if (module<240) {
	  if (eq<10) { // side A
	    return colM/32 + 5*(iModule%2);
	  }
	  else { // side C
	    return (159-colM)/32 + 5*(iModule%2);
	  }
	}
      }
    }
  }
  return 10; // error
}
//__________________________________________________________________________
Int_t AliSPDUtils::GetModuleNumber(UInt_t iDDL, UInt_t iModule) {
  if (iDDL<20 && iModule<12) return fgkDDLModuleMap[iDDL][iModule];
  else return 240;
}
//__________________________________________________________________________
UInt_t AliSPDUtils::GetOnlineColFromOffline(UInt_t module, UInt_t colM) {
  // offline->online (col)
  if (module<80) { // inner layer
    return colM%32;
  }
  else if (module<240) { // outer layer
    return colM%32;
  }
  return 32; // error
}
//__________________________________________________________________________
UInt_t AliSPDUtils::GetOnlineRowFromOffline(UInt_t module, UInt_t rowM) {
  // offline->online (row)
  if (module<80) { // inner layer
    return (255-rowM);
  }
  else if (module<240) { // outer layer
    return (255-rowM);
  }
  return 256; // error
}
//__________________________________________________________________________
UInt_t AliSPDUtils::GetOfflineColFromOnline(UInt_t eqId, UInt_t hs, UInt_t chip, UInt_t col) {
  // online->offline (col)
  if (eqId>=20 || hs>=6 || chip>=10 || col>=32) return 160; // error
  UInt_t offset = 32 * (chip % 5);
  if (hs<2) {
    if (eqId<10) {
      return 159 - (31-col + offset); // inner layer, side A
    }
    else {
      return col + offset; // inner layer, side C
    }
  }
  else {
    if (eqId<10) {
      return (col + offset); // outer layer, side A
    }
    else {
      return 159 - (31-col + offset); // outer layer, side C
    }
  }
}
//__________________________________________________________________________
UInt_t AliSPDUtils::GetOfflineRowFromOnline(UInt_t eqId, UInt_t hs, UInt_t chip, UInt_t row) {
  // online->offline (row)
  if (eqId>=20 || hs>=6 || chip>=10 || row>=256) return 256; // error
  return 255-row;
}
//__________________________________________________________________________
Bool_t AliSPDUtils::GetOnlineFromOfflineChipKey(UInt_t chipkey,UInt_t& eq, UInt_t& hs, UInt_t& chip){
  // online Eq, hs and chip from offline chipkey (V. Altini)
  if (chipkey>=1200) {
    TString errMess = Form("%d is not a valid Chip Key number",chipkey);
    printf("ERROR : %s \n", errMess.Data());
    return 0;
  }

  eq = GetOnlineEqIdFromOffline(chipkey/5);
  hs = GetOnlineHSFromOffline(chipkey/5);
  chip=chipkey%20;
  if(chip>9) chip=19-chip;

  return 1;
}

#endif
