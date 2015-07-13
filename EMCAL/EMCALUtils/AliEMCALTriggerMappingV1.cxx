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

#include "AliEMCALTriggerMapping.h"
#include "AliEMCALTriggerMappingV1.h"
#include "AliEMCALGeometry.h"
#include "AliLog.h"

ClassImp(AliEMCALTriggerMappingV1)

//________________________________________________________________________________________________
AliEMCALTriggerMappingV1::AliEMCALTriggerMappingV1() : AliEMCALTriggerMapping()
{
  // Ctor
  
  SetUniqueID(1);
}
//________________________________________________________________________________________________
AliEMCALTriggerMappingV1::AliEMCALTriggerMappingV1(Int_t ntru, const AliEMCALGeometry* geo) : AliEMCALTriggerMapping(ntru, geo)
{
  // Ctor
  
  SetUniqueID(1);
}

//________________________________________________________________________________________________
Bool_t AliEMCALTriggerMappingV1::GetAbsFastORIndexFromTRU(Int_t iTRU, Int_t iADC, Int_t& id) const
{
  //Trigger mapping method, get  FastOr Index from TRU

  if (iTRU > fNTRU - 1 || iTRU < 0 || iADC > 95 || iADC < 0) {
    AliError("TRU out of range!");
    return kFALSE;
  }

  id  = ( iTRU % 2 ) ? iADC%4 + 4 * (23 - int(iADC/4)) : (3 - iADC%4) + 4 * int(iADC/4);
  id += iTRU * 96;
  return kTRUE;
}

//________________________________________________________________________________________________
Bool_t AliEMCALTriggerMappingV1::GetTRUFromAbsFastORIndex(Int_t id, Int_t& iTRU, Int_t& iADC) const
{
  // Trigger mapping method, get TRU number from FastOr Index

  if (id > fNTRU * 96 - 1 || id < 0) {
    AliError("Fast-OR ID is out of range!");
    return kFALSE;
  }
	
  iTRU = id / 96;
  iADC = id % 96;
  iADC = ( iTRU % 2 ) ? iADC%4 + 4 * (23 - int(iADC/4)) : (3 - iADC%4) + 4 * int(iADC/4);
  return kTRUE;
}

//________________________________________________________________________________________________
Bool_t AliEMCALTriggerMappingV1::GetPositionInTRUFromAbsFastORIndex(Int_t id, Int_t& iTRU, Int_t& iEta, Int_t& iPhi) const
{
  // Trigger mapping method, get position in TRU from FasOr Index
	
  Int_t iADC = -1;	
  if (!GetTRUFromAbsFastORIndex(id, iTRU, iADC)) return kFALSE;
	
  Int_t x = iADC / 4;
  Int_t y = iADC % 4;
  if ( iTRU % 2 ) { // C side 
    iEta = 23 - x;
    iPhi =      y;
  }
  else {            // A side
    iEta =      x;
    iPhi =  3 - y;
  }
  return kTRUE;
}

//________________________________________________________________________________________________
Bool_t AliEMCALTriggerMappingV1::GetPositionInSMFromAbsFastORIndex(Int_t id, Int_t& iSM, Int_t& iEta, Int_t& iPhi) const
{
  //Trigger mapping method, get position in Super Module from FasOr Index

  Int_t iTRU = -1;
  if (!GetPositionInTRUFromAbsFastORIndex(id, iTRU, iEta, iPhi)) return kFALSE;
  if (iTRU % 2) { // C side
    iSM  = 2 * (int(int(iTRU / 2) / 3)) + 1;
  } else {        // A side
    iSM  = 2 * (int(int(iTRU / 2) / 3));
  }
  iPhi += 4 * int((iTRU % 6) / 2);
  return kTRUE;
}

//________________________________________________________________________________________________
Bool_t AliEMCALTriggerMappingV1::GetPositionInEMCALFromAbsFastORIndex(Int_t id, Int_t& iEta, Int_t& iPhi) const
{
  //Trigger mapping method, get position in EMCAL from FastOR index

  Int_t iSM = -1;
  if (GetPositionInSMFromAbsFastORIndex(id, iSM, iEta, iPhi)) {
    if (iSM % 2) iEta += 24; 
    iPhi += 12 * int(iSM / 2);
    return kTRUE;
  }
  return kFALSE;
}

//________________________________________________________________________________________________
Bool_t AliEMCALTriggerMappingV1::GetAbsFastORIndexFromPositionInTRU(Int_t iTRU, Int_t iEta, Int_t iPhi, Int_t& id) const
{
  //Trigger mapping method, get Index if FastOr from Position in TRU
  if (iTRU < 0 || iTRU > fNTRU - 1 || iEta < 0 || iEta > 23 || iPhi < 0 || iPhi > 3) {
    AliError(Form("Out of range! iTRU=%d, iEta=%d, iPhi=%d", iTRU, iEta, iPhi));	
    return kFALSE;
  }
  id =  iPhi  + 4 * iEta + iTRU * 96;
  return kTRUE;
}

//________________________________________________________________________________________________
Bool_t AliEMCALTriggerMappingV1::GetAbsFastORIndexFromPositionInSM(Int_t  iSM, Int_t iEta, Int_t iPhi, Int_t& id) const
{
  // Trigger mapping method, from position in SM Index get FastOR index 

  if (iSM < 0 || iSM >= 12 || iEta < 0 || iEta > 23 || iPhi < 0 || iPhi > 11) {
    AliError(Form("Out of range! iSM=%d, iEta=%d, iPhi=%d", iSM, iEta, iPhi));
    return kFALSE;
  }
  
  Int_t x = iEta;
  Int_t y = iPhi % 4;	
  Int_t iOff = (iSM % 2) ? 1 : 0;
  Int_t iTRU = 2 * int(iPhi / 4) + 6 * int(iSM / 2) + iOff;
  if (GetAbsFastORIndexFromPositionInTRU(iTRU, x, y, id)) {
    return kTRUE;
  }
  return kFALSE;
}

//________________________________________________________________________________________________
Bool_t AliEMCALTriggerMappingV1::GetAbsFastORIndexFromPositionInEMCAL(Int_t iEta, Int_t iPhi, Int_t& id) const
{
  // Trigger mapping method, from position in EMCAL Index get FastOR index 

  if (iEta < 0 || iEta > 47 || iPhi < 0 || iPhi > 63) {
    AliError(Form("Out of range! iEta: %2d iPhi: %2d", iEta, iPhi));
    return kFALSE;
  }

  Int_t s = int(iEta / 24) + 2 * int(iPhi / 12);
  Int_t x = iEta % 24;
  Int_t y = iPhi % 12;
  if (GetAbsFastORIndexFromPositionInSM(s, x, y, id)) {
      return kTRUE;
  }
  return kFALSE;
}

//________________________________________________________________________________________________
Bool_t AliEMCALTriggerMappingV1::GetFastORIndexFromCellIndex(Int_t id, Int_t& idx) const
{
  // Trigger mapping method, from cell index get FastOR index 

  Int_t iSupMod, nModule, nIphi, nIeta, iphim, ietam;
  Bool_t isOK = fGeometry->GetCellIndex( id, iSupMod, nModule, nIphi, nIeta );
  fGeometry->GetModulePhiEtaIndexInSModule( iSupMod, nModule, iphim, ietam );
  
  if (isOK && GetAbsFastORIndexFromPositionInSM(iSupMod, ietam, iphim, idx)) return kTRUE;
  return kFALSE;
}

//________________________________________________________________________________________________
Bool_t AliEMCALTriggerMappingV1::GetCellIndexFromFastORIndex(Int_t id, Int_t idx[4]) const
{
  // Trigger mapping method, from FASTOR index get cell index 

  Int_t iSM=-1, iEta=-1, iPhi=-1;
  if (GetPositionInSMFromAbsFastORIndex(id, iSM, iEta, iPhi)) {
    Int_t ix = 2 * iEta;
    Int_t iy = 2 * iPhi;
    for (Int_t i = 0; i < 2; i++) {
      for (Int_t j = 0; j < 2; j++) {
        idx[2 * i + j] = fGeometry->GetAbsCellIdFromCellIndexes(iSM, iy + i, ix + j);
      }
    }
    return kTRUE;
  }
  return kFALSE;
}

//________________________________________________________________________________________________
Bool_t AliEMCALTriggerMappingV1::GetTRUIndexFromSTUIndex(Int_t id, Int_t& idx, Int_t detector) const
{
  // STU mapping: TRU# 0 and 16 are missing 

   if (id > 31 || id < 0) {
     AliError(Form("TRU index out of range: %d",id));
     return kFALSE;
   }
   
   idx = GetTRUIndexFromSTUIndex(id, detector);
   return kTRUE;
}

//________________________________________________________________________________________________
Int_t AliEMCALTriggerMappingV1::GetTRUIndexFromSTUIndex(Int_t id, Int_t detector) const
{
  // STU mapping: TRU# 0 and 16 are missing 

  if (id > 31 || id < 0) {
    AliError(Form("TRU index out of range: %d",id));
  }
  return (id > 15) ? 2 * (31 - id) : 2 * (15 - id) + 1;
}

//________________________________________________________________________________________________
Bool_t AliEMCALTriggerMappingV1::GetTRUIndexFromOnlineIndex(Int_t id, Int_t& idx) const
{
  //Trigger mapping method, from STU index get TRU index 

   idx = GetTRUIndexFromOnlineIndex(id);
   if (idx > fNTRU - 1 || idx < 0) {
     AliError(Form("TRU index out of range: %d",idx));
     return kFALSE;
   }
   return kTRUE;
}

//________________________________________________________________________________________________
Int_t AliEMCALTriggerMappingV1::GetTRUIndexFromOnlineIndex(Int_t id) const
{
  //Trigger mapping method, from STU index get TRU index 
	
  if (id > fNTRU - 1 || id < 0) {
    AliError(Form("TRU index out of range: %d",id));
  }
  
  if (id == 31) {
    return 31;
  }

  Int_t idx = ((id % 6) < 3) ? 6 * int(id / 6) + 2 * (id % 3) : 6 * int(id / 6) + 2 * (2 - (id % 3)) + 1;

  return idx;
}

///
/// \return TRU  global offline number from:
/// \param hwAdd: hardware address
/// \param ddl number
/// \param sm: uper-module number
/// Used in AliEMCALTriggerRawDigitMaker::Add()
///
Int_t  AliEMCALTriggerMappingV1::GetTRUIndexFromOnlineHwAdd(Int_t hwAdd, Int_t ddl, Int_t sm) const
{  
  // 1/3 SMs
  if ( sm == 10 ) return 30;
  if ( sm == 11 ) return 31;
  
  // Full EMCal
  
  UShort_t iBranch = ( hwAdd >> 11 ) & 0x1;              // 0/1
  
  Int_t iTRU = ( (ddl << 1) | iBranch ) - 1; // 0..2
    
  iTRU  = (sm % 2) ? 2 * (2 - iTRU) + 1 : 2 * iTRU;
    
  iTRU += 6 * int(sm/2);

  if (iTRU > 31 || iTRU < 0) 
  {
    AliError(Form("TRU index out of range: %d",iTRU));
    return -1;
  }
  
  return iTRU;
}


//________________________________________________________________________________________________
Bool_t AliEMCALTriggerMappingV1::GetOnlineIndexFromTRUIndex(Int_t id, Int_t& idx) const
{
  //Trigger mapping method, from STU index get TRU index 
 
    idx = GetOnlineIndexFromTRUIndex(id);
    if (idx > fNTRU-1 || idx < 0)
   {
     AliError(Form("TRU index out of range: %d",idx));
     return kFALSE;
   }
   return kTRUE;
}
//________________________________________________________________________________________________
Int_t AliEMCALTriggerMappingV1::GetOnlineIndexFromTRUIndex(Int_t id) const
{
  //Trigger mapping method, from STU index get TRU index 
	
  if (id > fNTRU-1 || id < 0) 
  {
    AliError(Form("TRU index out of range: %d",id));
  }
  if (id == 31) {
    return 31;
  }

  Int_t idx = (id % 2) ? int((6 - (id % 6)) / 2) + 3 * (2 * int(id / 6) + 1) : 3 * int(id / 6) + int(id / 2);

  return idx;
}

//________________________________________________________________________________________________
Bool_t AliEMCALTriggerMappingV1::GetFastORIndexFromL0Index(Int_t iTRU, Int_t id, Int_t idx[], Int_t size) const
{
  //Trigger mapping method, from L0 index get FastOR index 

  if (size <= 0 || size > 4) {
    AliError("Size not supported!");
    return kFALSE;
  }
		
  Int_t motif[4] = {0, 1, 4, 5};
  switch (size) {
    case 1: // Cosmic trigger
      if (!GetAbsFastORIndexFromTRU(iTRU, id, idx[1])) return kFALSE;
      break;
    case 4: // 4 x 4
      for (Int_t k = 0; k < 4; k++) {
        Int_t iADC = motif[k] + 4 * int(id / 3) + (id % 3);
				
        if (!GetAbsFastORIndexFromTRU(iTRU, iADC, idx[k])) return kFALSE;
      }
      break;
    default:
      break;
  }
	
  return kTRUE;
}



