/**************************************************************************
 * Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id$  */

///////////////////////////////////////////////////////////////////////////////
///
/// This class provides storage container ITS SSD module callibration data
/// used by DA. 
/// 
/// Date: 18/07/2008
///////////////////////////////////////////////////////////////////////////////

#include "AliITSNoiseSSD.h"
#include "AliITSPedestalSSD.h"
#include "AliITSBadChannelsSSD.h"
#include "AliITSModuleDaSSD.h"
#include "TString.h"
#include "AliLog.h"

ClassImp(AliITSModuleDaSSD)


const Int_t   AliITSModuleDaSSD::fgkStripsPerModule   = 1536;   // Number of strips per SSD module
const Int_t   AliITSModuleDaSSD::fgkPNStripsPerModule = 768;    // Number of N/P strips per SSD module
const Int_t   AliITSModuleDaSSD::fgkStripsPerChip     = 128;    // Number of strips per chip HAL25
const UChar_t AliITSModuleDaSSD::fgkMaxAdNumber       = 9;      // MAx SSD FEROM AD number
const UChar_t AliITSModuleDaSSD::fgkMaxAdcNumber      = 13;     // MAx SSD FEROM ADC number
const Int_t   AliITSModuleDaSSD::fgkChipsPerModule    = 12;     // Number of HAL25 chips per SSD module



using namespace std;

//______________________________________________________________________________
AliITSModuleDaSSD::AliITSModuleDaSSD() :
  fEquipId(0),
  fEquipType(0),
  fDdlId(0),
  fAd(0),
  fAdc(0),
  fModuleId(0),
  fNumberOfStrips(0),
  fStrips(NULL),
  fNumberOfChips(0),
  fCm(NULL),
  fCmFerom(NULL),
  fEventsNumber(0)
{
// Default constructor
}


//______________________________________________________________________________
AliITSModuleDaSSD::AliITSModuleDaSSD(const UChar_t ddlID, const UChar_t ad, const UChar_t adc, const UShort_t moduleID) :
  fEquipId(0),
  fEquipType(0),
  fDdlId(ddlID),
  fAd(ad),
  fAdc(adc),
  fModuleId(moduleID),
  fNumberOfStrips(0),
  fStrips(NULL),
  fNumberOfChips(0),
  fCm(NULL),
  fCmFerom(NULL),
  fEventsNumber(0)
{
// Constructor, set module id data
}



//______________________________________________________________________________
AliITSModuleDaSSD::AliITSModuleDaSSD(const Int_t numberofstrips) :
  fEquipId(0),
  fEquipType(0),
  fDdlId(0),
  fAd(0),
  fAdc(0),
  fModuleId(0),
  fNumberOfStrips(0),
  fStrips(NULL),
  fNumberOfChips(0),
  fCm(NULL),
  fCmFerom(NULL),
  fEventsNumber(0)
{
// Constructor, allocates memory for AliITSChannelDaSSD* and TArrayS* array for CM calculated in FEROM
  if (numberofstrips != fgkStripsPerModule) 
    AliWarning(Form("AliITSModuleDaSSD: ALICE ITS SSD Module contains %i strips", fgkStripsPerModule));
  fStrips = new (nothrow) AliITSChannelDaSSD* [numberofstrips];
  if (fStrips) {
     fNumberOfStrips = numberofstrips;
     for (Int_t i = 0; i < numberofstrips; i++) fStrips[i]= NULL;
  } else {
     AliError(Form("AliITSModuleDaSSD: Error allocating memory for %i AliITSChannelDaSSD* objects!", numberofstrips));
     fNumberOfStrips = 0;
     fStrips = NULL;
  }  
  fCmFerom = new (nothrow) TArrayS [fgkChipsPerModule];
  if (!fCmFerom) {
     AliError(Form("AliITSModuleDaSSD: Error allocating memory for %i TArrayS objects!", fgkChipsPerModule));
     fCmFerom = NULL;
  }  
}


//______________________________________________________________________________
AliITSModuleDaSSD::AliITSModuleDaSSD(const Int_t numberofstrips, const Long_t eventsnumber) :
  fEquipId(0),
  fEquipType(0),
  fDdlId(0),
  fAd(0),
  fAdc(0),
  fModuleId(0),
  fNumberOfStrips(0),
  fStrips(NULL),
  fNumberOfChips(0),
  fCm(NULL),
  fCmFerom(NULL),
  fEventsNumber(0)
{
// Constructor, allocates memory for AliITSChannelDaSSD* and events data
  if (numberofstrips != fgkStripsPerModule) 
    AliWarning(Form("AliITSModuleDaSSD: ALICE ITS SSD Module contains %i strips", fgkStripsPerModule));
  fStrips = new (nothrow) AliITSChannelDaSSD* [numberofstrips];
  if (fStrips) {
     fNumberOfStrips = numberofstrips;
     memset(fStrips, 0, numberofstrips * sizeof(AliITSChannelDaSSD*));
     for (Int_t i = 0; i < fNumberOfStrips; i++) {
       fStrips[i] = new AliITSChannelDaSSD(i, eventsnumber);
       if (!fStrips[i]) AliError(Form("AliITSModuleDaSSD: Error allocating memory for AliITSChannelDaSSD %i-th object", i));
     }
  } else {
     AliError(Form("AliITSModuleDaSSD: Error allocating memory for %i AliITSChannelDaSSD* objects!", numberofstrips));
     fNumberOfStrips = 0;
     fStrips = NULL;
  }  
  fCmFerom = new (nothrow) TArrayS [fgkChipsPerModule];
  if (fCmFerom) {
    for (Int_t i = 0; i < fgkChipsPerModule; i++) {
      fCmFerom[i].Set(eventsnumber);
      fCmFerom[i].Reset(0);
    }  
  }
  else {
     AliError(Form("AliITSModuleDaSSD: Error allocating memory for %i TArrayS objects!", fgkChipsPerModule));
     fCmFerom = NULL;
  }  
}



//______________________________________________________________________________
AliITSModuleDaSSD::AliITSModuleDaSSD(const AliITSModuleDaSSD& module) :
  TObject(module),
  fEquipId(module.fEquipId),
  fEquipType(module.fEquipType),
  fDdlId(module.fDdlId),
  fAd(module.fAd),
  fAdc(module.fAdc),
  fModuleId(module.fModuleId),
  fNumberOfStrips(module.fNumberOfStrips),
  fStrips(NULL),
  fNumberOfChips(module.fNumberOfChips),
  fCm(NULL),
  fCmFerom(NULL),
  fEventsNumber(module.fEventsNumber)
{
// copy constructor
  if ((module.fNumberOfStrips > 0) && (module.fStrips)) {
    fStrips = new (nothrow) AliITSChannelDaSSD* [module.fNumberOfStrips];
    if (fStrips) {
      for (Int_t strind = 0; strind < module.fNumberOfStrips; strind++) {
        if (module.fStrips[strind]) {
	  fStrips[strind] = new AliITSChannelDaSSD(*(module.fStrips[strind]));
	  if (!fStrips[strind]) { 
	    AliError("AliITSModuleDaSSD: Error copy constructor");
            for (Int_t i = (strind - 1); i >= 0; i--) delete fStrips[strind];
	    delete [] fStrips;
	    fStrips = NULL;
	    break;
	  }
	} else fStrips[strind] = NULL; 
      }	  
    } else {
       AliError(Form("AliITSModuleDaSSD: Error allocating memory for %i AliITSChannelDaSSD* objects!", module.fNumberOfStrips));
       fNumberOfStrips = 0;
       fStrips = NULL;
    }  
  }
  if (module.fCm) {
    fCm = new (nothrow) TArrayF [module.fNumberOfChips];
    if (fCm) {
      for (Int_t chind = 0; chind < module.fNumberOfChips; chind++) fCm[chind] = module.fCm[chind]; 
    } else {
       AliError(Form("AliITSModuleDaSSD: Error allocating memory for %i TArrayF objects!", module.fNumberOfChips));
       fNumberOfChips = 0;
       fCm = NULL;
    }  
  }
  
  if (module.fCmFerom) {
    fCmFerom = new (nothrow) TArrayS [fgkChipsPerModule];
    if (fCmFerom) {
      for (Int_t chind = 0; chind < fgkChipsPerModule; chind++) fCmFerom[chind] = module.fCmFerom[chind]; 
    } else {
       AliError(Form("AliITSModuleDaSSD: Error allocating memory for %i TArrayS objects!", fgkChipsPerModule));
       fCmFerom = NULL;
    }  
  }
}



//______________________________________________________________________________
AliITSModuleDaSSD& AliITSModuleDaSSD::operator = (const AliITSModuleDaSSD& module)
{
// assignment operator
  if (this == &module)  return *this;  
  TObject::operator=(module);
  if (fStrips) {
    for (Long_t i = 0; i < fNumberOfStrips; i++) if (fStrips[i]) delete fStrips[i];
    delete [] fStrips;
    fStrips = NULL;
  } 
  fEquipId = module.fEquipId;
  fEquipType = module.fEquipType;
  fDdlId = module.fDdlId;
  fAd = module.fAd;
  fAdc = module.fAdc;
  fModuleId = module.fModuleId;
  fStrips = NULL;
  fNumberOfChips = module.fNumberOfChips;
  fCm = NULL;
  fEventsNumber = module.fEventsNumber;
  if ((module.fNumberOfStrips > 0) && (module.fStrips)) {
    fStrips = new (nothrow) AliITSChannelDaSSD* [module.fNumberOfStrips];
    if (fStrips) {
      memset(fStrips, 0, (sizeof(AliITSChannelDaSSD*) * module.fNumberOfStrips));
      for (Int_t strind = 0; strind < module.fNumberOfStrips; strind++) {
        if (module.fStrips[strind]) {
          fStrips[strind] = new AliITSChannelDaSSD(*(module.fStrips[strind]));
          if (!fStrips[strind]) { 
            AliError("AliITSModuleDaSSD: Error copy constructor");
            for (Int_t i = (strind - 1); i >= 0; i--) delete fStrips[strind];
            delete [] fStrips;
            fStrips = NULL;
            fNumberOfStrips = 0;
            break;
          }
        } else fStrips[strind] = NULL; 
      }
      fNumberOfStrips = module.fNumberOfStrips;
    } else {
       AliError(Form("AliITSModuleDaSSD: Error allocating memory for %i AliITSChannelDaSSD* objects!", module.fNumberOfStrips));
       fNumberOfStrips = 0;
       fStrips = NULL;
    }  
  }
  if (fCm) delete [] fCm;
  if (module.fCm) {
    fCm = new (nothrow) TArrayF [module.fNumberOfChips];
    if (fCm) {
      for (Int_t chind = 0; chind < module.fNumberOfChips; chind++) fCm[chind] = module.fCm[chind]; 
    } else {
       AliError(Form("AliITSModuleDaSSD: Error allocating memory for %i TArrayF objects!", module.fNumberOfChips));
       fNumberOfChips = 0;
       fCm = NULL;
    }  
  }  
  if (fCmFerom) delete [] fCmFerom;
  if (module.fCmFerom) {
    fCmFerom = new (nothrow) TArrayS [module.fNumberOfChips];
    if (fCmFerom) {
      for (Int_t chind = 0; chind < fgkChipsPerModule; chind++) fCmFerom[chind] = module.fCmFerom[chind]; 
    } else {
       AliError(Form("AliITSModuleDaSSD: Error allocating memory for %i TArrayS objects!", fgkChipsPerModule));
       fCmFerom = NULL;
    }  
  }  
  return *this;
}
    

    
//______________________________________________________________________________
AliITSModuleDaSSD::~AliITSModuleDaSSD()
{
// Destructor
  if (fStrips)
  {
    for (Long_t i = 0; i < fNumberOfStrips; i++)
    { 
      if (fStrips[i]) delete fStrips[i];
    }
    delete [] fStrips;
  } 
  if (fCm) delete [] fCm;
  if (fCmFerom) delete [] fCmFerom;
}



//______________________________________________________________________________
Bool_t AliITSModuleDaSSD::SetNumberOfStrips(const Int_t numberofstrips)
{
// Allocates memory for AliITSChannelDaSSD*
  if (fStrips) {
    for (Int_t i = 0; i < fNumberOfStrips; i++) if (fStrips[i]) delete fStrips[i];
    delete [] fStrips;
    fStrips = NULL;
  }  
  if (numberofstrips <= 0) {fNumberOfStrips = 0; return kTRUE; } 
  if (numberofstrips != fgkStripsPerModule) 
    AliWarning(Form("AliITSModuleDaSSD: ALICE ITS SSD Module contains %i strips", fgkStripsPerModule));
  fStrips = new (nothrow) AliITSChannelDaSSD* [numberofstrips];
  if (fStrips) {
     fNumberOfStrips = numberofstrips;
     memset(fStrips, 0, sizeof(AliITSChannelDaSSD*) * numberofstrips);
     return kTRUE;
  } else {
     AliError(Form("AliITSModuleDaSSD: Error allocating memory for %i AliITSChannelDaSSD* objects!", numberofstrips));
     fNumberOfStrips = 0;
     fStrips = NULL;
     return kFALSE;
  }  
}


//______________________________________________________________________________
Bool_t AliITSModuleDaSSD::SetNumberOfChips(const Int_t nchips)
{
// Allocate nchips TArrayF objects to save Common Mode
  DeleteCM();
  if (nchips <= 0) {fNumberOfChips = 0; return kTRUE; } 
  if (nchips != fgkChipsPerModule) 
    AliWarning(Form("AliITSModuleDaSSD: ALICE ITS SSD Module contains %i HAL25 chips", fgkChipsPerModule));
  fCm = new (nothrow) TArrayF [nchips];
  if (fCm) {
     fNumberOfChips = nchips;
     return kTRUE;
  } else {
     AliError(Form("AliITSModuleDaSSD: Error allocating memory for %i TArrayF objects!", nchips));
     fNumberOfChips = 0;
     fCm = NULL;
     return kFALSE;
  }  
}

  
//______________________________________________________________________________
Bool_t AliITSModuleDaSSD::SetModuleIdData (const UChar_t ddlID, const UChar_t ad, const UChar_t adc, const Short_t moduleID)
{
// SetModuleIdData
  if (ad > fgkMaxAdNumber) {
    AliWarning(Form("AliITSModuleDaSSD: Wrong AD number: %i", ad));
    return kFALSE;
  }  
  if (adc > fgkMaxAdcNumber || ForbiddenAdcNumber(adc)) {
    AliWarning(Form("AliITSModuleDaSSD: Wrong ADC number: %i", adc));
    return kFALSE;
  }  
  fDdlId = ddlID;
  fAd = ad;
  fAdc = adc;
  fModuleId = moduleID;
  return kTRUE;
}


//______________________________________________________________________________
void AliITSModuleDaSSD::SetModuleFEEId (const UChar_t ddlID, const UChar_t ad, const UChar_t adc)
{
// Set id data of FEE connected to the Module
  fDdlId = ddlID;
  fAd = ad;
  fAdc = adc;
}


//______________________________________________________________________________
void AliITSModuleDaSSD::SetModuleRorcId (const Int_t equipid, const Int_t equiptype)
{
// Set data to access FEROM registres via DDL
  fEquipId = equipid; 
  fEquipType = equiptype;
}


//______________________________________________________________________________
Bool_t AliITSModuleDaSSD::SetEventsNumber(const Long_t eventsnumber)
{
// Allocate the memory for the events data
  Int_t i;
  if (!fStrips) return kFALSE;
  for (i = 0; i < fNumberOfStrips; i++) {
    if (fStrips[i]) { 
      if (!fStrips[i]->SetEvenetsNumber(eventsnumber)) {
        for (Int_t j = 0; j < i; j++) fStrips[j]->DeleteSignal();
        AliError(Form("AliITSModuleDaSSD: Error allocating memory for i% events for module %i, strip %i",eventsnumber, (Int_t)fModuleId, i));
        return kFALSE;
      }
    }
    else { 
      if (!(fStrips[i] = new AliITSChannelDaSSD(i, eventsnumber))) {
        for (Int_t j = 0; j < i; j++) delete fStrips[j];
        delete [] fStrips;
        fNumberOfStrips = 0;
        fStrips = NULL;
        AliError(Form("AliITSModuleDaSSD: Error allocating memory for strip %i of module %i!", (Int_t)fModuleId, i));
        return kFALSE;
      }
    }
  } 
  if (fCmFerom) {
    for (Int_t ie = 0; ie < fgkChipsPerModule; ie++) {
      fCmFerom[ie].Set(eventsnumber);
      fCmFerom[ie].Reset(0);
    }  
  }
  else  AliError("AliITSModuleDaSSD: No memory was allocated for fCmFerom!");
  fEventsNumber = eventsnumber;
  return kTRUE;
}



//______________________________________________________________________________
Bool_t AliITSModuleDaSSD::SetCM (const Float_t cm, const Int_t chipn, const Int_t evn)
{ 
// Set value of CM for a given chip and event 
  if ((!fCm) || (chipn >= fNumberOfChips)) return kFALSE;
  if (evn >= fCm[chipn].GetSize()) return kFALSE;
  else fCm[chipn][evn] = cm;
  return kTRUE;
}



//______________________________________________________________________________
Float_t  AliITSModuleDaSSD::GetCM(const Int_t chipn, const Long_t evn)   const 
{ 
// Get value of CM for a given chip and event 
  if ((!fCm) || (chipn >= fNumberOfChips)) return 0.0f;
  if (evn >= fCm[chipn].GetSize()) return 0.0f;
  else return fCm[chipn][evn];
}



//______________________________________________________________________________
Bool_t  AliITSModuleDaSSD::SetCMFeromEventsNumber(const Long_t eventsnumber)
{
// Allocates memory for the values of CM calculated in Ferom
  if (!fCmFerom) return kFALSE;
  for (Int_t chipind = 0; chipind < fgkChipsPerModule; chipind++) {
      fCmFerom[chipind].Set(eventsnumber);
      fCmFerom[chipind].Reset(0);
  }
  return kTRUE;
}

//______________________________________________________________________________
Bool_t AliITSModuleDaSSD::SetCMFerom (const Short_t cm, const Int_t chipn, const Int_t evn)
{ 
// Set value of FeromCM for a given chip and event 
  if ((!fCmFerom) || (chipn >= fgkChipsPerModule)) return kFALSE;
  if (evn >= fCmFerom[chipn].GetSize()) return kFALSE;
  else fCmFerom[chipn][evn] = cm;
  return kTRUE;
}



//______________________________________________________________________________
Short_t  AliITSModuleDaSSD::GetCMFerom(const Int_t chipn, const Long_t evn)   const 
{ 
// Get value of FeromCM for a given chip and event 
  if ((!fCmFerom) || (chipn >= fgkChipsPerModule)) return 0;
  if (evn >= fCmFerom[chipn].GetSize()) return 0;
  else return fCmFerom[chipn][evn];
}



UChar_t AliITSModuleDaSSD::CheckIfBad(const Int_t stripn) const
{
//Applies the bad channel creteria and set the appropriate flags for returned valie 
  UInt_t          bcflags = 0;
  const UInt_t    WOffsetMask = 0x000003FF;
  if (!fStrips[stripn]) bcflags |= 3;
  else {
    if (fStrips[stripn]->GetNoiseCM() == AliITSChannelDaSSD::GetUndefinedValue()) bcflags |= 8;
    if (fStrips[stripn]->GetNoiseCM() > 20) bcflags |= 8;
    if (fStrips[stripn]->GetNoiseCM() < 1) bcflags |= 16;
    if (fStrips[stripn]->GetPedestal() > ((WOffsetMask >> 1) - 1))  bcflags |= 4;
    else if ((-(fStrips[stripn]->GetPedestal())) > (WOffsetMask >> 1))  bcflags |= 4;
    if (bcflags) bcflags |= 3;
  }
  return bcflags;
}
