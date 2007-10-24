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
///////////////////////////////////////////////////////////////////////////////

#include "AliITSNoiseSSD.h"
#include "AliITSModuleDaSSD.h"

ClassImp(AliITSModuleDaSSD)

using namespace std;

AliITSModuleDaSSD::AliITSModuleDaSSD() :
  fEquipId(0),
  fEquipType(0),
  fDdlId(0),
  fAd(0),
  fAdc(0),
  fModuleId(0),
  fNumberOfStrips(0),
  fStrips(NULL),
  fEventsNumber(0)
{
// Default constructor
}


AliITSModuleDaSSD::AliITSModuleDaSSD(const UChar_t ddlID, const UChar_t ad, const UChar_t adc, const UShort_t moduleID) :
  fEquipId(0),
  fEquipType(0),
  fDdlId(ddlID),
  fAd(ad),
  fAdc(adc),
  fModuleId(moduleID),
  fNumberOfStrips(0),
  fStrips(NULL),
  fEventsNumber(0)
{
// Constructor, set module id data
}



AliITSModuleDaSSD::AliITSModuleDaSSD(const Int_t numberofstrips) :
  fEquipId(0),
  fEquipType(0),
  fDdlId(0),
  fAd(0),
  fAdc(0),
  fModuleId(0),
  fNumberOfStrips(0),
  fStrips(NULL),
  fEventsNumber(0)
{
// Constructor, allocates memory for AliITSChannelDaSSD*
  if (numberofstrips != fgkStripsPerModule) 
    Warning("AliITSModuleDaSSD", "ALICE ITS SSD Module contains %i strips", fgkStripsPerModule);
  fStrips = new (nothrow) AliITSChannelDaSSD* [numberofstrips];
  if (fStrips) {
     fNumberOfStrips = numberofstrips;
     for (Int_t i = 0; i < numberofstrips; i++) fStrips[i]= NULL;
  } else {
     Error("AliITSModuleDaSSD", "Error allocating memory for %i AliITSChannelDaSSD* objects!", numberofstrips);
     fNumberOfStrips = 0;
     fStrips = NULL;
  }  
}


AliITSModuleDaSSD::AliITSModuleDaSSD(const Int_t numberofstrips, const Long_t eventsnumber) :
  fEquipId(0),
  fEquipType(0),
  fDdlId(0),
  fAd(0),
  fAdc(0),
  fModuleId(0),
  fNumberOfStrips(0),
  fStrips(NULL),
  fEventsNumber(0)
{
// Constructor, allocates memory for AliITSChannelDaSSD* and events data
  if (numberofstrips != fgkStripsPerModule) 
    Warning("AliITSModuleDaSSD", "ALICE ITS SSD Module contains %i strips", fgkStripsPerModule);
  fStrips = new (nothrow) AliITSChannelDaSSD* [numberofstrips];
  if (fStrips) {
     fNumberOfStrips = numberofstrips;
     memset(fStrips, 0, numberofstrips * sizeof(AliITSChannelDaSSD*));
     for (Int_t i = 0; i < fNumberOfStrips; i++) {
       fStrips[i] = new AliITSChannelDaSSD(i, eventsnumber);
       if (!fStrips[i]) Error("AliITSModuleDaSSD", "Error allocating memory for AliITSChannelDaSSD %i-th object", i);
     }
  } else {
     Error("AliITSModuleDaSSD", "Error allocating memory for %i AliITSChannelDaSSD* objects!", numberofstrips);
     fNumberOfStrips = 0;
     fStrips = NULL;
  }  
}



AliITSModuleDaSSD::AliITSModuleDaSSD(const AliITSModuleDaSSD& module) :
  TObject(module),
  fEquipId(module.fEquipId),
  fEquipType(module.fEquipType),
  fDdlId(module.fDdlId),
  fAd(module.fAd),
  fAdc(module.fAdc),
  fModuleId(module.fModuleId),
  fNumberOfStrips(module.fNumberOfStrips),
  fStrips(module.fStrips),
  fEventsNumber(module.fEventsNumber)
{
// copy constructor

  Fatal("AliITSModuleDaSSD", "copy constructor not implemented");
}



AliITSModuleDaSSD& AliITSModuleDaSSD::operator = (const AliITSModuleDaSSD& module)
{
// assignment operator

  Fatal("AliITSModuleDaSSD: operator =", "assignment operator not implemented");
  return *this;
}
    

    
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
}


  
Bool_t AliITSModuleDaSSD::SetModuleIdData (const UChar_t ddlID, const UChar_t ad, const UChar_t adc, const UShort_t moduleID)
{
// SetModuleIdData
  if (ad > fgkMaxAdNumber) {
    Warning("AliITSModuleDaSSD", "Wrong AD number: %i", ad);
    return kFALSE;
  }  
  if (adc > fgkMaxAdcNumber || ForbiddenAdcNumber(adc)) {
    Warning("AliITSModuleDaSSD", "Wrong ADC number: %i", adc);
    return kFALSE;
  }  
  fDdlId = ddlID;
  fAd = ad;
  fAdc = adc;
  fModuleId = moduleID;
  return kTRUE;
}



void AliITSModuleDaSSD::SetModuleFEEId (const UChar_t ddlID, const UChar_t ad, const UChar_t adc)
{
// Set id data of FEE connected to the Module
  fDdlId = ddlID;
  fAd = ad;
  fAdc = adc;
}


void AliITSModuleDaSSD::SetModuleRorcId (const Int_t equipid, const Int_t equiptype)
{
// Set data to access FEROM registres via DDL
  fEquipId = equipid; 
  fEquipType = equiptype;
}


Bool_t AliITSModuleDaSSD::SetEventsNumber(const Long_t eventsnumber)
{
// Allocate the memory for the events data
  Int_t i;
  if (!fStrips) return kFALSE;
  for (i = 0; i < fNumberOfStrips; i++) {
    if (fStrips[i])  
      if (!fStrips[i]->SetEvenetsNumber(eventsnumber)) {
        for (Int_t j = 0; j < i; j++) fStrips[j]->DeleteSignal();
        Error("AliITSModuleDaSSD", "Error allocating memory for i% events for module %i, strip %i", 
	                            eventsnumber, (Int_t)fModuleId, i);
        return kFALSE;
      }
    else 
      if (!(fStrips[i] = new AliITSChannelDaSSD(i, eventsnumber))) {
        for (Int_t j = 0; j < i; j++) delete fStrips[j];
        delete [] fStrips;
        fNumberOfStrips = 0;
        fStrips = NULL;
        Error("AliITSModuleDaSSD", "Error allocating memory for strip %i of module %i!", (Int_t)fModuleId, i);
        return kFALSE;
      }
  } 
  return kTRUE;
}



AliITSNoiseSSD* AliITSModuleDaSSD::GetCalibrationSSDModule() const
{
// Creates the AliITSNoiseSSD objects with callibration data
  AliITSNoiseSSD  *mc;
  if (!fStrips) return NULL;
  mc = new AliITSNoiseSSD();
  mc->SetMod(fModuleId);
  mc->SetNNoiseP(fgkPNStripsPerModule);
  mc->SetNNoiseN(fgkPNStripsPerModule);
  for (Int_t i = 0; i < fNumberOfStrips; i++) {
    if (!fStrips[i]) {
      delete mc;
      return NULL;
    }
    if (i < fgkPNStripsPerModule)
          mc->AddNoiseP(i, fStrips[i]->GetNoise());
    else  mc->AddNoiseN((i - fgkPNStripsPerModule), fStrips[i]->GetNoise());                     
  }
  return mc;
}
