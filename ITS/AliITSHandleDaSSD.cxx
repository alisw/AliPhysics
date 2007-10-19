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

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
///
/// This class provides ITS SSD data handling
/// used by DA. 
///
///////////////////////////////////////////////////////////////////////////////

#include <Riostream.h> 
#include "AliITSHandleDaSSD.h"
#include <math.h>
#include <sstream>
#include <string>
#include "event.h"
#include "TFile.h"
#include "TObjArray.h"
#include "AliITSNoiseSSD.h"
#include "AliRawReaderDate.h"

#include "AliITSChannelDaSSD.h"


ClassImp(AliITSHandleDaSSD)

using namespace std;


AliITSHandleDaSSD::AliITSHandleDaSSD() :
  fNumberOfModules(0),
  fModules(NULL),
  fLdcId(0),
  fRunId(0)
{
// Default constructor
}


AliITSHandleDaSSD::AliITSHandleDaSSD(const Int_t numberofmodules) :
  fNumberOfModules(0),
  fModules(NULL),
  fLdcId(0),
  fRunId(0)
{
// Constructor allocates memory for AliITSModuleDaSSD objects
  if (numberofmodules > fgkNumberOfSSDModules) 
      cout << "ALICE ITS SSD contains " <<  fgkNumberOfSSDModules << "modules!"<< endl;
  fModules = new (nothrow) AliITSModuleDaSSD* [numberofmodules];
  if (fModules) {
    fNumberOfModules = numberofmodules;
    memset(fModules, 0, sizeof(AliITSModuleDaSSD*) * numberofmodules);
  } else {
     Error("AliITSHandleDaSSD", "Error allocating memory for %i AliITSModulesDaSSD* objects!", numberofmodules);
     fNumberOfModules = 0;
     fModules = NULL;
  }
}


AliITSHandleDaSSD::AliITSHandleDaSSD(const AliITSHandleDaSSD& ssdadldc) :
  TObject(ssdadldc),
  fNumberOfModules(ssdadldc.fNumberOfModules),
  fModules(ssdadldc.fModules),
  fLdcId(ssdadldc.fLdcId),
  fRunId(ssdadldc.fRunId)
{
  // copy constructor

  Fatal("AliITSHandleDaSSD", "copy constructor not implemented");
}


AliITSHandleDaSSD& AliITSHandleDaSSD::operator = (const AliITSHandleDaSSD& ssdadldc)
{
// assignment operator

  Fatal("operator =", "assignment operator not implemented");
  return *this;
}


AliITSHandleDaSSD::~AliITSHandleDaSSD()
{
// Default destructor 
  if (fModules) 
  {
    for (Int_t i = 0; i < fNumberOfModules; i++)
    { 
      if (fModules[i]) delete fModules[i];
    }
    delete [] fModules;
  }
}



AliITSModuleDaSSD* AliITSHandleDaSSD::GetModule (const UChar_t ddlID, const UChar_t ad, const UChar_t adc) const
{
// Retrieve AliITSModuleDaSSD object from the array
   if (!fModules) return NULL;
   for (Int_t i = 0; i < fNumberOfModules; i++) {
     if (fModules[i]) {
       if (    (fModules[i]->GetDdlId() == ddlID)
            && (fModules[i]->GetAD() == ad)
            && (fModules[i]->GetADC() == adc))
	    return fModules[i];
     }	    
   }
   return NULL;
}



AliITSChannelDaSSD* AliITSHandleDaSSD::GetStrip (const UChar_t ddlID, const UChar_t ad, const UChar_t adc, const UShort_t stripID) const
{
// Retrieve AliITSChannalDaSSD object from the array
  for (Int_t modind = 0; modind < fNumberOfModules; modind++) {
    if (    (fModules[modind]->GetDdlId() == ddlID)
         && (fModules[modind]->GetAD() == ad)
         && (fModules[modind]->GetADC() == adc)  ) 
        {
       return fModules[modind]->GetStrip(stripID);
    }
  }
  Error("AliITSHandleDaSSD", "Error GetStrip (%i, %i, %i, %i), strip not found, returns NULL!", ddlID, ad, adc, stripID);
  return NULL;
}



Bool_t AliITSHandleDaSSD::SetModule(AliITSModuleDaSSD *const module, const Int_t index)
{ 
// Assign array element with AliITSModuleDaSSD object
   if ((index < fNumberOfModules) && (index >= 0)) 
     { 
        if (fModules[index]) delete fModules[index];
	fModules[index] = module;
	return kTRUE;
      }
   else return kFALSE;   		       
}


Bool_t AliITSHandleDaSSD::SetNumberOfModules (const Int_t numberofmodules)
{
// Allocates memory for AliITSModuleDaSSD objects
  if (numberofmodules > fgkNumberOfSSDModules)
    Warning("AliITSHandleDaSSD", "The number of modules %i you use exceeds the maximum %i for ALICE ITS SSD", numberofmodules, fgkNumberOfSSDModules);
  if (fModules) { delete [] fModules; fModules = NULL; }
  fModules = new (nothrow) AliITSModuleDaSSD* [numberofmodules];
  if (fModules) {
    fNumberOfModules = numberofmodules;
    memset(fModules, 0, sizeof(AliITSModuleDaSSD*) * numberofmodules);
    return kTRUE;
  } else {
     Error("AliITSHandleDaSSD", "Error allocating memory for %i AliITSModulesDaSSD* objects!", numberofmodules);
     fNumberOfModules = 0;
     fModules = NULL;
  }
  return kFALSE;
}


Bool_t AliITSHandleDaSSD::ReadCalibrationDataFile (const char* fileName, const Long_t eventsnumber)
{
// Reads raw data from file
  AliRawReaderDate    *rawreaderdate = new AliRawReaderDate(fileName, 0);
  AliITSModuleDaSSD   *module;
  AliITSChannelDaSSD  *strip;
  Long_t            datasize, eventind = 0;
  Int_t             nofstrips, eqbelsize;
  UShort_t          modind;
  long32           *data;
  if (!fModules) {
    Error("AliITSHandleDaSSD", "Error ReadCalibrationDataFile: no structure was allocated for data");
    return kFALSE;
  }
  rawreaderdate->SelectEvents(-1);
  while (rawreaderdate->NextEvent()) {
    if ((rawreaderdate->GetType() != PHYSICS_EVENT) && (rawreaderdate->GetType() != CALIBRATION_EVENT)) continue;
    fLdcId = rawreaderdate->GetLDCId();
    fRunId = rawreaderdate->GetRunNumber(); 
    modind = 0;
    while (rawreaderdate->ReadNextData((UChar_t*&)data)) {
      Int_t     equipid    = rawreaderdate->GetEquipmentId();              //  EquipmentID required to access to rorc
      Int_t     equiptype  = rawreaderdate->GetEquipmentType();            //
      UChar_t   ddlID      = (UChar_t)rawreaderdate->GetDDLID();           // GetDDLID(); index of DDL, ITS SSD: 33-48
      datasize = rawreaderdate->GetDataSize();
      eqbelsize = rawreaderdate->GetEquipmentElementSize();
      if ( (datasize % eqbelsize) || (eqbelsize != sizeof(long32)) ) {
        Error("AliITSHandleDaSSD", "Error ReadCalibrationDataFile: event data size %i is not an integer of equipment data size %i", datasize,
               eqbelsize);
        return kFALSE;
      }
      nofstrips = (Int_t) (datasize / eqbelsize);		 
      for (Int_t strind = 0; strind < nofstrips; strind++) {
        UChar_t   ad      = (UChar_t) (data[strind] >> 28) & 0x0000000F;  // index of AD module     0-9
        UChar_t   adc     = (UChar_t) (data[strind] >> 24) & 0x0000000F;  // index of ADC module    0-5, 8-13
        UShort_t  stripID = (UShort_t)(data[strind] >> 12) & 0x000007FF;  // strip number           0-1535
        Short_t   signal  = (Short_t)(data[strind] & 0x00000FFF);
                  signal  = (signal > AliITSChannelDaSSD::GetUnderflowConst()) ? (signal - 2 * AliITSChannelDaSSD::GetUnderflowConst()) 
		                                                               : signal;
        if (!(module = GetModule(ddlID, ad, adc))) {
          module = new AliITSModuleDaSSD(AliITSModuleDaSSD::GetStripsPerModuleConst());
	  if (!module->SetModuleIdData (ddlID, ad, adc, modind)) return kFALSE;
          module->SetModuleRorcId (equipid, equiptype);
	  fModules[modind++] = module;
        }
        if (!(strip = module->GetStrip(stripID))) {
          strip = new AliITSChannelDaSSD(stripID, eventsnumber);
          module->SetStrip(strip, stripID);
        }
        strip->SetSignal(eventind, signal);
     }
     if (modind) cout << "The memory was allocated for " <<  modind << " modules." << endl;
   }
   if (++eventind > eventsnumber) break;
  }
  delete rawreaderdate;
  return RelocateModules();
}


Bool_t AliITSHandleDaSSD::RelocateModules()
{
// Relocate memory for AliITSModuleDaSSD object array
  Int_t         nm = 0;
  AliITSModuleDaSSD **marray = NULL;
  for (Int_t modind = 0; modind < fNumberOfModules; modind++)
    if (fModules[modind]) nm += 1;
  if (nm == fNumberOfModules) return kTRUE;
  marray = new (nothrow) AliITSModuleDaSSD* [nm];
  if (!marray) {
    Error("AliITSHandleDaSSD", "Error relocating memory for %i AliITSModuleDaSSD* objects!", nm);
    return kFALSE;
  }
  nm = 0;
  for (Int_t modind = 0; modind < fNumberOfModules; modind++)  
    if (fModules[modind]) marray[nm++] = fModules[modind];
  delete [] fModules;
  fModules = marray;
  fNumberOfModules = nm;
  return kTRUE;
}


Bool_t AliITSHandleDaSSD::CalculatePedestal()
{
// Calculates Pedestal
  Float_t     pedestal;
  Short_t    *signal;
  AliITSChannelDaSSD *strip;
  Long_t      ovev, ev;
  if (!fModules) return kFALSE;
  for (Int_t modind = 0; modind < fNumberOfModules; modind++) {
    if (!fModules[modind]) {
      Error("AliITSHandleDaSSD", "Error CalculatePedestal(): No AliITSChannelDaSSD object with index %i is allocated in AliITSModuleDaSSD\n", modind);
      return kFALSE;
    }
    for (Int_t strind = 0; strind < fModules[modind]->GetNumberOfStrips(); strind++) {
      if (!(strip = fModules[modind]->GetStrip(strind))) return kFALSE;
      if (!(signal = strip->GetSignal())) {
        Error("AliITSHandleDaSSD", "Error CalculatePedestal(): there are no events data for module[%i] strip[%i]->GetSignal()", modind, strind);
        return kFALSE;
      }
      pedestal = 0.0f;
      ovev = 0l;
      for (ev = 0; ev < strip->GetEventsNumber(); ev++)
        if (SignalOutOfRange(signal[ev])) ovev += 1;
//	else pedestal += signal[ev];
	else pedestal = ((ev - ovev)) ? pedestal + (signal[ev] - pedestal) / static_cast<Float_t>(ev - ovev + 1) : signal[ev];
	if (ev == ovev) pedestal = AliITSChannelDaSSD::GetUndefinedValue();
//      if ((Long_t n = strip->GetEventsNumber() - ovev)) pedestal /= static_cast<Float_t>(n);
//      else return kFALSE;
      strip->SetPedestal(pedestal);	
    }
  }
  return kTRUE;
}



Bool_t AliITSHandleDaSSD::CalculateNoise()
{
// Calculates Noise
  AliITSChannelDaSSD *strip;
  Short_t    *signal;
  Float_t     noise;
  Long_t      ovev, n;
  if (!fModules) return kFALSE;
  for (Int_t modind = 0; modind < fNumberOfModules; modind++) {
    if (!fModules[modind]) {
      Error("AliITSHandleDaSSD", "Error CalculateNoise(): No AliITSChannelDaSSD object with index %i is allocated in AliITSModuleDaSSD\n", modind);
      return kFALSE;
    }
    for (Int_t strind = 0; strind < fModules[modind]->GetNumberOfStrips(); strind++) {
      if (!(strip = fModules[modind]->GetStrip(strind))) return kFALSE;
      if (!(signal = strip->GetSignal())) {
        Error("AliITSHandleDaSSD", "Error CalculateNoise(): there are no events data for module[%i] strip[%i]->GetSignal()", modind, strind);
        return kFALSE;
      }
      Double_t nsum = 0.0L;
      ovev = 0l;
      for (Long_t ev = 0; ev < strip->GetEventsNumber(); ev++) {
        if (SignalOutOfRange(signal[ev])) ovev += 1;
	else nsum += pow((signal[ev] - strip->GetPedestal()), 2);
//	else nsum = ((ev - ovev)) ? nsum + (pow((signal[ev] - strip->GetPedestal()), 2) - nsum) / static_cast<Double_t>(ev - ovev)
//	                          : pow((signal[ev] - strip->GetPedestal()), 2);
      } 
//      noise =  sqrt(nsum);
      if ((n = strip->GetEventsNumber() - ovev - 1) > 0) noise =  sqrt(nsum / (Float_t)(n));
      else  noise = AliITSChannelDaSSD::GetUndefinedValue();
      strip->SetNoise(noise);
    }
  }
  return kTRUE;
}



Bool_t AliITSHandleDaSSD::CalculateCM(const Int_t modind, const Int_t stripind, Float_t* const cm)
{
// Calculates CM
  AliITSChannelDaSSD *strip = NULL;
  Short_t    *signal;
  Long_t      ovstr, evn, n;
  if ((stripind + AliITSModuleDaSSD::GetStripsPerChip()) > fModules[modind]->GetNumberOfStrips()) return kFALSE;
  if (!(strip = fModules[modind]->GetStrip(stripind))) return kFALSE;
  evn = fModules[modind]->GetStrip(stripind)->GetEventsNumber();
  for (Long_t ev = 0; ev < evn; ev++) {
    Double_t cmsum = 0.0L;
    ovstr = 0l;
    for (Int_t strind = stripind; strind < (stripind + AliITSModuleDaSSD::GetStripsPerChip()); strind++) {
      if (!(strip = fModules[modind]->GetStrip(strind))) return kFALSE;
      if (!(signal = strip->GetSignal())) {
        Error("AliITSHandleDaSSD", "Error CalculateCM: there are no events data for module[%i] strip[%i]->GetSignal()", modind, strind);
        return kFALSE;
      }
      if (SignalOutOfRange(signal[ev])) ovstr += 1;
      else cmsum += (signal[ev] - strip->GetPedestal());
    }
    if ((n = AliITSModuleDaSSD::GetStripsPerChip() - ovstr)) cm[ev] = cmsum / (Float_t)(n);
    else cm[ev] = 0.0;
  } 
  return kTRUE; 
}



Bool_t  AliITSHandleDaSSD::CalculateNoiseCM()
{
// Calculates Noise with CM correction
  Short_t     *signal;
  AliITSChannelDaSSD  *strip = NULL;
  Float_t      noise, *cm = NULL;
  Long_t       ovev, n;
  if (!fModules) return kFALSE;
  for (Int_t modind = 0; modind < fNumberOfModules; modind++) {
    if (!fModules[modind]) {
      Error("AliITSHandleDaSSD", "Error CalculateNoiseCM(): No AliITSChannelDaSSD object with index %i is allocated in AliITSModuleDaSSD", modind);
      return kFALSE;
    }
    for (Int_t strind = 0; strind < fModules[modind]->GetNumberOfStrips(); strind++) {
      if (!(strip = fModules[modind]->GetStrip(strind))) return kFALSE;
      if (!(signal = strip->GetSignal())) {
        Error("AliITSHandleDaSSD", "Error CalculateNoiseCM(): there are no events data for module[%i] strip[%i]->GetSignal()", modind, strind);
        return kFALSE;
      }
      if (!(strind % AliITSModuleDaSSD::GetStripsPerChip())) {
        if (!cm) {
          try { 
	    cm = new Float_t [strip->GetEventsNumber()];  }
          catch (bad_alloc&) {
            Warning("AliITSHandleDaSSD", "Noise calculation with common mode correction failed becouse of memory allocation problems.");
            return kFALSE; 
          }  
        }
// calculate cm;
        if (!CalculateCM(modind, strind, cm)) return kFALSE;
      }
// calculate noise;	
      Double_t nsum = 0.0L;
      ovev = 0l;
      for (Long_t ev = 0; ev < strip->GetEventsNumber(); ev++) {
        if (SignalOutOfRange(signal[ev])) ovev += 1;
	else nsum += pow((signal[ev] - strip->GetPedestal() - cm[ev]), 2);
//	else nsum = ((ev - ovev)) ? nsum + (pow((signal[ev] - strip->GetPedestal() - cm[ev]), 2) - nsum) / static_cast<Double_t>(ev - ovev)
//	                          : pow((signal[ev] - strip->GetPedestal() - cm[ev]), 2);
      } 
//      noise =  sqrt(nsum);
      if ((n = strip->GetEventsNumber() - ovev - 1) > 0) noise =  sqrt(nsum / (Float_t)(n));
      else  noise = AliITSChannelDaSSD::GetUndefinedValue();
      strip->SetNoise(noise);
      }
    }
  if (cm) delete [] cm;
  return kTRUE;
}



TObjArray* AliITSHandleDaSSD::GetCalibrationSSDLDC() const
{
// Fill in the array for OCDB 
  TObjArray *ldcc;
  if (!fModules) return NULL;
  ldcc = new TObjArray(fNumberOfModules, 0);
  for (Int_t i = 0; i < fNumberOfModules; i++) {
    if (!fModules[i]) {
      delete ldcc;
      return NULL;
    }
    ldcc->AddAt(dynamic_cast<TObject*>(fModules[i]->GetCalibrationSSDModule()), i);
  }
  ldcc->Compress();
  return ldcc;
}


Bool_t AliITSHandleDaSSD::SaveCalibrationSSDLDC(string& dafname) const
{
// Save Calibration data locally
  ostringstream   dadatafilename;
  TObjArray      *ldcc;
  if (!fModules) return kFALSE;
  ldcc = new TObjArray(fNumberOfModules, 0);
  for (Int_t i = 0; i < fNumberOfModules; i++) {
    if (!fModules[i]) {
      delete ldcc;
      return kFALSE;
    }
    ldcc->AddAt(dynamic_cast<TObject*>(fModules[i]->GetCalibrationSSDModule()), i);
  }
  ldcc->Compress();
  dadatafilename << dafname << "/ITSSSDda_" << fLdcId << "_" << fRunId << ".root";
  dafname = dadatafilename.str();
  TFile *fileRun = new TFile (dadatafilename.str().data(),"RECREATE");
  if (fileRun->IsZombie()) {
    Error("AliITSHandleDaSSD", "SaveCalibrationSSDLDC() error open file %s", dadatafilename.str().data());
    ldcc->Delete();
    delete fileRun;
    return kFALSE;
  }
  fileRun->WriteTObject(ldcc);
  fileRun->Close();
  ldcc->Delete();
  delete fileRun;
  return kTRUE;
}
