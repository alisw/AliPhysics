

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
  if (numberofstrips != fgkStripsPerModule) 
    Warning("AliITSModuleDaSSD", "ALICE ITS SSD Module contains %d strips", fgkStripsPerModule);
  try  {
     fStrips = new AliITSChannelDaSSD* [numberofstrips];
     fNumberOfStrips = numberofstrips;
     for (Int_t i = 0; i < numberofstrips; i++) fStrips[i]= NULL;
  }
  catch (bad_alloc&) {
     Error("AliITSModuleDaSSD", "Error allocating memory for %d AliITSChannelDaSSD objects!", numberofstrips);
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
  if (numberofstrips != fgkStripsPerModule) 
    Warning("AliITSModuleDaSSD", "ALICE ITS SSD Module contains %d strips", fgkStripsPerModule);
  try  {
     fStrips = new AliITSChannelDaSSD* [numberofstrips];
     fNumberOfStrips = numberofstrips;
  }
  catch (bad_alloc&) {
     Error("AliITSModuleDaSSD", "Error allocating memory for %d AliITSChannelDaSSD objects!", numberofstrips);
     fNumberOfStrips = 0;
     fStrips = NULL;
  }
  if (fStrips) {
    Int_t  i;
    try {
       for (i = 0; i < fNumberOfStrips; i++) fStrips[i] = new AliITSChannelDaSSD(i, eventsnumber);
    }  
    catch (bad_alloc&) {
       Error("AliITSModuleDaSSD", "Error allocating memory for %d-th AliITSChannelDaSSD objects!", i);
       for (Int_t j = 0; j < i; j++) delete fStrips[j];
       delete [] fStrips;
       fNumberOfStrips = 0;
       fStrips = NULL;
    }
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
  if (fStrips)
  {
    for (Long_t i = 0; i < fNumberOfStrips; i++)
    { 
      if (fStrips[i]) delete fStrips[i];
//      if (!(i % 100)) cout << "Deleted fStrips[i], i = " << i << endl;
    }
    delete [] fStrips;
  } 
}


  
Bool_t AliITSModuleDaSSD::SetModuleIdData (const UChar_t ddlID, const UChar_t ad, const UChar_t adc, const UShort_t moduleID)
{
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
  fDdlId = ddlID;
  fAd = ad;
  fAdc = adc;
}


void AliITSModuleDaSSD::SetModuleRorcId (const Int_t equipid, const Int_t equiptype)
{
  fEquipId = equipid; 
  fEquipType = equiptype;
}


Bool_t AliITSModuleDaSSD::SetEventsNumber(const Long_t eventsnumber)
{
  Int_t i;
  if (!fStrips) return kFALSE;
  try {
     for (i = 0; i < fNumberOfStrips; i++) {
       if (fStrips[i]) fStrips[i]->SetEvenetsNumber(eventsnumber);
       else fStrips[i] = new AliITSChannelDaSSD(i, eventsnumber);
     } 
  }  
  catch (bad_alloc&) {
     Error("AliITSModuleDaSSD", "Error allocating memory for %d-th AliITSChannelDaSSD objects!", i);
     for (Int_t j = 0; j < i; j++) delete fStrips[j];
     delete [] fStrips;
     fNumberOfStrips = 0;
     fStrips = NULL;
     return kFALSE;
  }
  return kTRUE;
}



AliITSNoiseSSD* AliITSModuleDaSSD::GetCalibrationSSDModule() const
{
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
