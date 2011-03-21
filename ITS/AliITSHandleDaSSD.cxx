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
//  Author: Oleksandr Borysov
//  Date: 18/07/2008
///////////////////////////////////////////////////////////////////////////////

#include <string>
#include <fstream>
#include <new>
#include <Riostream.h> 
#include "AliITSHandleDaSSD.h"
#include <math.h>
#include <limits.h>
#include "event.h"
#include "TFile.h"
#include "TString.h"
#include "TMath.h"
#include "AliLog.h"
#include "AliITSNoiseSSDv2.h"
#include "AliITSPedestalSSDv2.h"
#include "AliITSBadChannelsSSDv2.h"
#include "AliITSRawStreamSSD.h"
#include "AliRawReaderDate.h"
#include "AliITSRawStreamSSD.h"
#include "AliITSChannelDaSSD.h"


ClassImp(AliITSHandleDaSSD)

using namespace std;


const Int_t    AliITSHandleDaSSD::fgkNumberOfSSDModules = 1698;       // Number of SSD modules in ITS
const Int_t    AliITSHandleDaSSD::fgkNumberOfSSDModulesPerDdl = 108;  // Number of SSD modules in DDL
const Int_t    AliITSHandleDaSSD::fgkNumberOfSSDModulesPerSlot = 12;  // Number of SSD modules in Slot
const Short_t  AliITSHandleDaSSD::fgkMinSSDModuleId = 500;            // Initial SSD modules number
const Int_t    AliITSHandleDaSSD::fgkNumberOfSSDSlotsPerDDL = 9;      // Number of SSD slots per DDL
const Int_t    AliITSHandleDaSSD::fgkNumberOfSSDDDLs = 16;            // Number of SSD modules in Slot
const Float_t  AliITSHandleDaSSD::fgkPedestalThresholdFactor = 3.0;   // Defalt value for fPedestalThresholdFactor 
const Float_t  AliITSHandleDaSSD::fgkCmThresholdFactor = 3.0;         // Defalt value for fCmThresholdFactor

const UInt_t   AliITSHandleDaSSD::fgkZsBitMask      = 0x0000003F;     // Bit mask for FEROM ZS
const UInt_t   AliITSHandleDaSSD::fgkOffSetBitMask  = 0x000003FF;     // Bit mask for FEROM Offset correction
const UInt_t   AliITSHandleDaSSD::fgkBadChannelMask = 0x00000001;     // Mask to suppress the channel from the bad channel list
const Int_t    AliITSHandleDaSSD::fgkAdcPerDBlock   = 6;              // FEROM configuration file constant

//______________________________________________________________________________
AliITSHandleDaSSD::AliITSHandleDaSSD() :
  fRawDataFileName(NULL),
  fNumberOfModules(0),
  fModules(NULL),
  fModIndProcessed(0),
  fModIndRead(0),
  fModIndex(NULL),
  fEqIndex(0),
  fNumberOfEvents(0),
  fBadChannelsList(NULL),
  fDDLModuleMap(NULL),
  fALaddersOff(0),
  fCLaddersOff(0),
  fLdcId(0),
  fRunId(0),
  fPedestalThresholdFactor(fgkPedestalThresholdFactor),
  fCmThresholdFactor(fgkCmThresholdFactor),
  fZsDefault(-1),
  fOffsetDefault(INT_MAX),
  fZsMinimum(2),
  fMergeBCLists(1),
  fZsFactor(3.0)
{
// Default constructor
}


//______________________________________________________________________________
AliITSHandleDaSSD::AliITSHandleDaSSD(Char_t *rdfname) :
  fRawDataFileName(NULL),
  fNumberOfModules(0),
  fModules(NULL),
  fModIndProcessed(0),
  fModIndRead(0),
  fModIndex(NULL),
  fEqIndex(0),
  fNumberOfEvents(0),
  fBadChannelsList(NULL),
  fDDLModuleMap(NULL),
  fALaddersOff(0),
  fCLaddersOff(0),
  fLdcId(0),
  fRunId(0),
  fPedestalThresholdFactor(fgkPedestalThresholdFactor) ,
  fCmThresholdFactor(fgkCmThresholdFactor),
  fZsDefault(-1),
  fOffsetDefault(INT_MAX),
  fZsMinimum(2),
  fMergeBCLists(1),
  fZsFactor(3.0)
{
  if (!Init(rdfname)) AliError("AliITSHandleDaSSD::AliITSHandleDaSSD() initialization error!");
}


//______________________________________________________________________________
AliITSHandleDaSSD::AliITSHandleDaSSD(const AliITSHandleDaSSD& ssdadldc) :
  TObject(ssdadldc),
  fRawDataFileName(ssdadldc.fRawDataFileName),
  fNumberOfModules(ssdadldc.fNumberOfModules),
  fModules(NULL),
  fModIndProcessed(ssdadldc.fModIndProcessed),
  fModIndRead(ssdadldc.fModIndRead),
  fModIndex(NULL),
  fEqIndex(0),
  fNumberOfEvents(ssdadldc.fNumberOfEvents),
  fBadChannelsList(NULL),
  fDDLModuleMap(NULL),
  fALaddersOff(ssdadldc.fALaddersOff),
  fCLaddersOff(ssdadldc.fCLaddersOff),
  fLdcId(ssdadldc.fLdcId),
  fRunId(ssdadldc.fRunId),
  fPedestalThresholdFactor(ssdadldc.fPedestalThresholdFactor),
  fCmThresholdFactor(ssdadldc.fCmThresholdFactor),
  fZsDefault(ssdadldc.fZsDefault),
  fOffsetDefault(ssdadldc.fOffsetDefault),
  fZsMinimum(ssdadldc.fZsMinimum),
  fMergeBCLists(ssdadldc.fMergeBCLists),
  fZsFactor(ssdadldc.fZsFactor)
{
  // copy constructor
  if ((ssdadldc.fNumberOfModules > 0) && (ssdadldc.fModules)) {
    fModules = new (nothrow) AliITSModuleDaSSD* [ssdadldc.fNumberOfModules];
    if (fModules) {
      for (Int_t modind = 0; modind < ssdadldc.fNumberOfModules; modind++) {
        if (ssdadldc.fModules[modind]) {
	  fModules[modind] = new AliITSModuleDaSSD(*(ssdadldc.fModules[modind]));
	  if (!fModules[modind]) { 
	    AliError("AliITSHandleDaSSD: Error copy constructor");
            for (Int_t i = (modind - 1); i >= 0; i--) delete fModules[modind];
	    delete [] fModules;
	    fModules = NULL;
	    break;
	  }
	} else fModules[modind] = NULL; 
      }	  
    } else {
      AliError(Form("AliITSHandleDaSSD: Error allocating memory for %i AliITSModulesDaSSD* objects!", ssdadldc.fNumberOfModules));
      fNumberOfModules = 0;
      fModules = NULL;
    }
  }
  if (ssdadldc.fBadChannelsList) AliWarning("fBadChannelsList is not copied by copy constructor, use other methods to init it!");
  if (ssdadldc.fDDLModuleMap) AliWarning("fDDLModuleMap is not copied by copy constructor, use other methods to init it!");
}


//______________________________________________________________________________
AliITSHandleDaSSD& AliITSHandleDaSSD::operator = (const AliITSHandleDaSSD& ssdadldc)
{
// assignment operator
  if (this == &ssdadldc)  return *this;
  TObject::operator=(ssdadldc);
  if (fModules) {
    for (Int_t i = 0; i < fNumberOfModules; i++) if (fModules[i]) delete fModules[i];
    delete [] fModules; 
    fModules = NULL;
  }
  if (fModIndex) { delete [] fModIndex; fModIndex = NULL; }
  if ((ssdadldc.fNumberOfModules > 0) && (ssdadldc.fModules)) {
    fModules = new (nothrow) AliITSModuleDaSSD* [ssdadldc.fNumberOfModules];
    if (fModules) {
      for (Int_t modind = 0; modind < ssdadldc.fNumberOfModules; modind++) {
        if (ssdadldc.fModules[modind]) {
	      fModules[modind] = new AliITSModuleDaSSD(*(ssdadldc.fModules[modind]));
          if (!fModules[modind]) { 
	        AliError("AliITSHandleDaSSD: Error assignment operator");
            for (Int_t i = (modind - 1); i >= 0; i--) delete fModules[modind];
            delete [] fModules;
            fModules = NULL;
            break;
          }
        } else fModules[modind] = NULL; 
      }	  
    } else {
      AliError(Form("AliITSHandleDaSSD: Error allocating memory for %i AliITSModulesDaSSD* objects!", ssdadldc.fNumberOfModules));
      fNumberOfModules = 0;
      fModules = NULL;
    }
  }
  fRawDataFileName = NULL;
  fModIndProcessed = 0;
  fModIndRead = 0;
  fModIndex = NULL;
  fEqIndex = ssdadldc.fEqIndex;
  fNumberOfEvents = ssdadldc.fNumberOfEvents;
  fLdcId = ssdadldc.fLdcId;
  fRunId = ssdadldc.fRunId;
  fPedestalThresholdFactor = ssdadldc.fPedestalThresholdFactor;
  fCmThresholdFactor = ssdadldc.fCmThresholdFactor;
  fZsDefault = ssdadldc.fZsDefault;
  fOffsetDefault = ssdadldc.fOffsetDefault;
  fZsMinimum = ssdadldc.fZsMinimum;
  fMergeBCLists = ssdadldc.fMergeBCLists;
  fZsFactor = ssdadldc.fZsFactor;
  fALaddersOff = ssdadldc.fALaddersOff;
  fCLaddersOff = ssdadldc.fCLaddersOff;
  fBadChannelsList = NULL;
  fDDLModuleMap = NULL;
  fModIndex = NULL;
  if (ssdadldc.fBadChannelsList) AliWarning("fBadChannelsList is not copied by assignment operator, use other methods to init it!");
  if (ssdadldc.fDDLModuleMap) AliWarning("fDDLModuleMap is not copied by assignment operator, use other methods to init it!");
  return *this;
}


//______________________________________________________________________________
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
  if (fModIndex) delete [] fModIndex;
  if (fBadChannelsList)  delete fBadChannelsList;
  if (fDDLModuleMap) delete [] fDDLModuleMap;
}



//______________________________________________________________________________
void AliITSHandleDaSSD::Reset()
{
// Delete array of AliITSModuleDaSSD* objects. 
  if (fModules) {
    for (Int_t i = 0; i < fNumberOfModules; i++) if (fModules[i]) delete fModules[i];
    delete [] fModules;
    fModules = NULL;
  }
  if (fModIndex) { delete [] fModIndex; fModIndex = NULL; }
/*
  if (fBadChannelsList) {
    delete fBadChannelsList;
    fBadChannelsList = NULL; 
  }    
  if (fDDLModuleMap) { delete [] fDDLModuleMap; fDDLModuleMap = NULL; }
*/
  fALaddersOff.Set(0);
  fCLaddersOff.Set(0);
  fRawDataFileName = NULL;
  fModIndProcessed = fModIndRead = 0;
  fNumberOfEvents = 0;
  fLdcId = fRunId = 0;
  fPedestalThresholdFactor = fgkPedestalThresholdFactor;
  fCmThresholdFactor = fgkCmThresholdFactor;
}



//______________________________________________________________________________
Bool_t AliITSHandleDaSSD::Init(Char_t *rdfname)
{
// Read raw data file and set initial and configuration parameters
  Long_t physeventind = 0, strn = 0, nofstripsev, nofstrips = 0;
  Int_t  nofeqipmentev, nofeqipment = 0, eqn = 0; 
  AliRawReaderDate  *rawreaderdate = NULL;
  UChar_t *data = NULL;
  Long_t datasize = 0, eqbelsize = 1;

  rawreaderdate = new AliRawReaderDate(rdfname, 0);
  if (!(rawreaderdate->GetAttributes() || rawreaderdate->GetEventId())) {
    AliError(Form("AliITSHandleDaSSD: Error reading raw data file %s by RawReaderDate", rdfname));
    MakeZombie();
    return kFALSE;
  }
  if (rawreaderdate->NextEvent()) {
    fRunId = rawreaderdate->GetRunNumber(); 
    rawreaderdate->RewindEvents();
  } else { MakeZombie(); return kFALSE; }
  if (fModules) Reset();
  //rawreaderdate->SelectEvents(-1);
  rawreaderdate->Select("ITSSSD");  
  nofstrips = 0;
  while (rawreaderdate->NextEvent()) {
    if ((rawreaderdate->GetType() != PHYSICS_EVENT) && (rawreaderdate->GetType() != CALIBRATION_EVENT)) continue;
    nofstripsev = 0;
    nofeqipmentev = 0;
    while (rawreaderdate->ReadNextData(data)) {
      fLdcId = rawreaderdate->GetLDCId();
      nofeqipmentev += 1;
      datasize = rawreaderdate->GetDataSize();
      eqbelsize = rawreaderdate->GetEquipmentElementSize();
      if ( datasize % eqbelsize ) {
        AliError(Form("AliITSHandleDaSSD: Error Init(%s): event data size %ld is not an integer of equipment data size %ld", 
	                            rdfname, datasize, eqbelsize));
        MakeZombie();
	    return kFALSE;
      }
      nofstripsev += (Int_t) (datasize / eqbelsize);
    }
    if (physeventind++) {
      if (nofstrips != nofstripsev) AliWarning(Form("AliITSHandleDaSSD: number of strips varies from event to event, ev = %ld, %ld", 
                                                     physeventind, nofstripsev));
      if (nofeqipment != nofeqipmentev) AliWarning("AliITSHandleDaSSD: number of DDLs varies from event to event");
    }
    nofstrips = nofstripsev;
    nofeqipment = nofeqipmentev;
    if (strn < nofstrips)  strn = nofstrips;
    if (eqn < nofeqipment) eqn = nofeqipment;
  }
  delete rawreaderdate;
  if ((physeventind > 0) && (strn > 0))
  {
    fNumberOfEvents = physeventind;
    fRawDataFileName = rdfname;
    fEqIndex.Set(eqn);
    fEqIndex.Reset(-1);
    fModIndex = new (nothrow) Int_t [fgkNumberOfSSDModulesPerDdl * eqn];
    if (fModIndex) 
      for (Int_t i = 0; i < fgkNumberOfSSDModulesPerDdl * eqn; i++) fModIndex[i] = -1; 
    else AliWarning(Form("AliITSHandleDaSSD: Error Init(%s): Index array for %i modules was not created", 
                                     rdfname, fgkNumberOfSSDModulesPerDdl * eqn));
    if (SetNumberOfModules(fgkNumberOfSSDModulesPerDdl * eqn)) {
      TString str = TString::Format("Max number of equipment: %d, max number of channels: %ld\n", eqn, strn);
      DumpInitData(str.Data());
      return kTRUE;
    }  
  }
  MakeZombie();
  return kFALSE;
}



//______________________________________________________________________________
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


Int_t AliITSHandleDaSSD::GetModuleIndex (const UChar_t ddlID, const UChar_t ad, const UChar_t adc) const
{
// Retrieve the position of AliITSModuleDaSSD object in the array
   if (!fModules) return -1;
   for (Int_t i = 0; i < fNumberOfModules; i++) {
     if (fModules[i]) {
       if (    (fModules[i]->GetDdlId() == ddlID)
            && (fModules[i]->GetAD() == ad)
            && (fModules[i]->GetADC() == adc))
	    return i;
     }	    
   }
   return -1;
}



//______________________________________________________________________________
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
  AliError(Form("AliITSHandleDaSSD: Error GetStrip (%i, %i, %i, %i), strip not found, returns NULL!",
                ddlID, ad, adc, stripID));
  return NULL;
}



//______________________________________________________________________________
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



//______________________________________________________________________________
Bool_t AliITSHandleDaSSD::SetNumberOfModules (const Int_t numberofmodules)
{
// Allocates memory for AliITSModuleDaSSD objects
  if (numberofmodules > fgkNumberOfSSDModules)
    AliWarning(Form("AliITSHandleDaSSD: the number of modules %i you use exceeds the maximum %i for ALICE ITS SSD", 
                     numberofmodules, fgkNumberOfSSDModules));
  if (fModules) {
    for (Int_t i = 0; i < fNumberOfModules; i++) if (fModules[i]) delete fModules[i];
    delete [] fModules; 
    fModules = NULL;
  }
  fModules = new (nothrow) AliITSModuleDaSSD* [numberofmodules];
  if (fModules) {
    fNumberOfModules = numberofmodules;
    memset(fModules, 0, sizeof(AliITSModuleDaSSD*) * numberofmodules);
    return kTRUE;
  } else {
     AliError(Form("AliITSHandleDaSSD: Error allocating memory for %i AliITSModulesDaSSD* objects!", numberofmodules));
     fNumberOfModules = 0;
     fModules = NULL;
  }
  return kFALSE;
}



//______________________________________________________________________________
Bool_t AliITSHandleDaSSD::ReadStaticBadChannelsMap(const Char_t *filename)
{
// Reads Static Bad Channels Map from the file
  TFile *bcfile;
  if (!filename) {
    AliWarning("No file name is specified for Static Bad Channels Map!");
    return kFALSE;
  } 
  bcfile = new TFile(filename, "READ");
  if (bcfile->IsZombie()) {
    AliWarning(Form("Error reading file %s with Static Bad Channels Map!", filename));
    return kFALSE;
  }
  bcfile->GetObject("AliITSBadChannelsSSDv2;1", fBadChannelsList);
  if (!fBadChannelsList) {
    AliWarning("Error fBadChannelsList == NULL!");
    bcfile->Close();
    delete bcfile;
    return kFALSE;
  }
  bcfile->Close();
  delete bcfile;
  return kTRUE;
}



Bool_t AliITSHandleDaSSD::ReadDDLModuleMap(const Char_t *filename)
{
// Reads the SSD DDL Map from the file
  ifstream             ddlmfile;
  AliRawReaderDate    *rwr = NULL;
  AliITSRawStreamSSD  *rsm = NULL;
  void                *event = NULL;
  if (fDDLModuleMap) { delete [] fDDLModuleMap; fDDLModuleMap = NULL;}
  fDDLModuleMap = new (nothrow) Int_t [fgkNumberOfSSDDDLs * fgkNumberOfSSDModulesPerDdl];
  if (!fDDLModuleMap) {
    AliWarning("Error allocation memory for DDL Map!");
    return kFALSE;
  }    
  if (!filename) {
    AliWarning("No file name is specified for SSD DDL Map, using the one from AliITSRawStreamSSD!");
    rwr = new AliRawReaderDate(event);
    rsm = new AliITSRawStreamSSD(rwr);
    rsm->Setv11HybridDDLMapping();
    for (Int_t ddli = 0; ddli < fgkNumberOfSSDDDLs; ddli++)
      for (Int_t mi = 0; mi < fgkNumberOfSSDModulesPerDdl; mi++)
        fDDLModuleMap[(ddli * fgkNumberOfSSDModulesPerDdl + mi)] = rsm->GetModuleNumber(ddli, mi);
    if (rsm) delete rsm;
    if (rwr) delete rwr;	
    return kTRUE;
  } 
  ddlmfile.open(filename, ios::in);
  if (!ddlmfile.is_open()) {
    AliWarning(Form("Error reading file %s with SSD DDL Map!", filename));
    if (fDDLModuleMap) { delete [] fDDLModuleMap; fDDLModuleMap = NULL;}
    return kFALSE;
  }
  Int_t ind = 0;
  while((!ddlmfile.eof()) && (ind < (fgkNumberOfSSDDDLs * fgkNumberOfSSDModulesPerDdl))) {
    ddlmfile >> fDDLModuleMap[ind++];
    if (ddlmfile.fail()) AliError(Form("Error extracting number from the DDL map file %s, ind: %d ", filename, ind));
  }
  if (ind != (fgkNumberOfSSDDDLs * fgkNumberOfSSDModulesPerDdl))
    AliWarning(Form("Only %i (< %i) entries were read from DDL Map!", ind, (fgkNumberOfSSDDDLs * fgkNumberOfSSDModulesPerDdl)));
  ddlmfile.close();
  return kTRUE;
}



//______________________________________________________________________________
Int_t AliITSHandleDaSSD::ReadCalibrationDataFile (char* fileName, const Long_t eventsnumber)
{
// Reads raw data from file
  if (!Init(fileName)){
    AliError("AliITSHandleDaSSD: Error ReadCalibrationDataFile");
    return kFALSE;
  }
  fNumberOfEvents = eventsnumber;
  return ReadModuleRawData (fNumberOfModules);  
}



//______________________________________________________________________________
Int_t AliITSHandleDaSSD::ReadModuleRawData (const Int_t modulesnumber)
{
// Reads raw data from file
  AliRawReader        *rawreaderdate = NULL;
  AliITSRawStreamSSD  *stream = NULL;
  AliITSModuleDaSSD   *module;
  AliITSChannelDaSSD  *strip;
  Long_t            eventind = 0;
  Int_t             nofeqipment, eqind;
  Short_t           equipid, prequipid;
  Short_t           modind;
  if (!(rawreaderdate = new AliRawReaderDate(fRawDataFileName, 0))) return 0;
  if (!fModules) {
    AliError("AliITSHandleDaSSD: Error ReadModuleRawData: no structure was allocated for data");
    return 0;
  }
  if (!fDDLModuleMap) if (!ReadDDLModuleMap()) AliWarning("DDL map is not defined, ModuleID will be set to 0!");
  stream = new AliITSRawStreamSSD(rawreaderdate);
  stream->Setv11HybridDDLMapping();
  //rawreaderdate->SelectEvents(-1);
  rawreaderdate->Select("ITSSSD");  
  modind = 0;
  nofeqipment = 0;
  while (rawreaderdate->NextEvent()) {
    if ((rawreaderdate->GetType() != PHYSICS_EVENT) && (rawreaderdate->GetType() != CALIBRATION_EVENT)) continue;
    prequipid = -1;
    eqind = -1;
    while (stream->Next()) {
      equipid = rawreaderdate->GetEquipmentId(); 
      if ((equipid != prequipid)) {
        if ((eqind = GetEqIndex(equipid)) < 0) { fEqIndex.AddAt(equipid, nofeqipment); eqind = nofeqipment++; }
        prequipid = equipid;
      }
      Int_t     equiptype  = rawreaderdate->GetEquipmentType();
      UChar_t   ddlID      = (UChar_t)rawreaderdate->GetDDLID();
      UChar_t   ad      = stream->GetAD();
      UChar_t   adc     = stream->GetADC();
      UShort_t  stripID = stream->GetSideFlag() ? AliITSChannelDaSSD::GetMaxStripIdConst() - stream->GetStrip() : stream->GetStrip();
      Short_t   signal  = stream->GetSignal();

      Int_t indpos = (eqind * fgkNumberOfSSDModulesPerDdl)
                   + ((ad - 1) * fgkNumberOfSSDModulesPerSlot) + (adc < 6 ? adc : (adc - 2));
      Int_t modpos = fModIndex[indpos];
      if (((modpos > 0) && (modpos < fModIndRead)) || ((modpos < 0) && (modind == modulesnumber))) continue;
      if ((modpos < 0) && (modind < modulesnumber)) {
        module = new AliITSModuleDaSSD(AliITSModuleDaSSD::GetStripsPerModuleConst());
        Int_t mddli;
        if (fDDLModuleMap) mddli = RetrieveModuleId(ddlID, ad, adc);
        else  mddli = 0;
        if (!module->SetModuleIdData (ddlID, ad, adc, mddli)) return 0;
        module->SetModuleRorcId (equipid, equiptype);
        modpos = fModIndRead + modind;
        modind += 1;
        fModules[modpos] = module;
        fModIndex[indpos] = modpos;
      } 
      if (stripID < AliITSModuleDaSSD::GetStripsPerModuleConst()) {
        if (!(strip = fModules[modpos]->GetStrip(stripID))) {
          strip = new AliITSChannelDaSSD(stripID, fNumberOfEvents);
          fModules[modpos]->SetStrip(strip, stripID);
        }
        strip->SetSignal(eventind, signal);
	  } else  {
        if (!(fModules[modpos]->GetCMFerom())) {
          fModules[modpos]->AllocateCMFeromArray();
          fModules[modpos]->SetCMFeromEventsNumber(fNumberOfEvents);
		}
	    fModules[modpos]->SetCMFerom(signal, (stripID - AliITSModuleDaSSD::GetStripsPerModuleConst()), eventind);
	  }
    }
    if (++eventind > fNumberOfEvents) break;
  }
  delete stream;
  delete rawreaderdate;
  if (modind) cout << "The memory was allocated for " <<  modind << " modules." << endl;
  fModIndRead += modind;
  if (modind < modulesnumber) RelocateModules();
  return modind;
}



//______________________________________________________________________________
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
    AliError(Form("AliITSHandleDaSSD: Error relocating memory for %i AliITSModuleDaSSD* objects!", nm));
    return kFALSE;
  }
  nm = 0;
  for (Int_t modind = 0; modind < fNumberOfModules; modind++)  
    if (fModules[modind]) marray[nm++] = fModules[modind];
  delete [] fModules;
  fModules = marray;
  fNumberOfModules = fModIndRead = nm;
  return kTRUE;
}



//______________________________________________________________________________
Bool_t AliITSHandleDaSSD::AddFeromCm(const AliITSModuleDaSSD *const module)
{
// Restore the original signal value adding CM calculated and subtracted in ferom
  AliITSChannelDaSSD *strip;
  Short_t            *signal, *cmferom;

  if (!module) return kFALSE;
  if (!module->GetCMFerom()) return kTRUE;
  for (Int_t chipind = 0; chipind < AliITSModuleDaSSD::GetChipsPerModuleConst(); chipind++) {
    if (!(cmferom = module->GetCMFerom(chipind))) {
      AliWarning(Form("AliITSHandleDaSSD: There is no Ferom CM values for chip %i, module %i!", chipind, module->GetModuleId()));
      continue;
    }
    for (Int_t strind = (chipind * AliITSModuleDaSSD::GetStripsPerChip()); 
               strind < ((chipind + 1) * AliITSModuleDaSSD::GetStripsPerChip()); strind++) {
      if (!(strip = module->GetStrip(strind))) continue;
      if (!(signal = strip->GetSignal())) continue;
//      if (strip->GetEventsNumber() != module->GetEventsNumber()) return kFALSE;
      Long_t ovev = 0;
      for (Long_t ev = 0; ev < strip->GetEventsNumber(); ev++) {
        if (SignalOutOfRange(signal[ev]) || SignalOutOfRange(cmferom[ev])) ovev += 1;
        else {
          Short_t signal1 = signal[ev] + cmferom[ev];
          strip->SetSignal(ev, signal1);
	    }  
      }	
    } 
  }  
  return kTRUE;
}



//______________________________________________________________________________
Bool_t AliITSHandleDaSSD::CalculatePedestal(const AliITSModuleDaSSD *const module)
{
// Calculates Pedestal
  AliITSChannelDaSSD *strip;
  Double_t            pedestal, noise;
  Short_t            *signal;
  Long_t              ovev, ev, n;
  if (!module) return kFALSE;
  for (Int_t strind = 0; strind < module->GetNumberOfStrips(); strind++) {
    if (!(strip = module->GetStrip(strind))) continue;
    if (!(signal = strip->GetSignal())) {
      AliError(Form("AliITSHandleDaSSD: Error CalculatePedestal(): there are no events data for module[%i] strip[%i]->GetSignal()",
                     module->GetModuleId(), strind));
      continue;	
    }
//************* pedestal first pass ****************
    pedestal = 0.0L;
    ovev = 0l;
    for (ev = 0; ev < strip->GetEventsNumber(); ev++)
      if (SignalOutOfRange(signal[ev])) ovev += 1;
      else pedestal = ((ev - ovev)) ? pedestal + (signal[ev] - pedestal) / static_cast<Double_t>(ev - ovev + 1) : signal[ev];
    if (ev == ovev) pedestal = AliITSChannelDaSSD::GetUndefinedValue();
    strip->SetPedestal(static_cast<Float_t>(pedestal));	
//************* noise *******************************
    Double_t nsum = 0.0L;
    ovev = 0l;
    for (ev = 0; ev < strip->GetEventsNumber(); ev++) {
      if (SignalOutOfRange(signal[ev])) ovev += 1;
      else nsum += pow((signal[ev] - strip->GetPedestal()), 2);
    } 
    if ((n = strip->GetEventsNumber() - ovev - 1) > 0) noise =  sqrt(nsum / (Double_t)(n));
    else  noise = AliITSChannelDaSSD::GetUndefinedValue();
    strip->SetNoise(static_cast<Float_t>(noise));
//************* pedestal second pass ****************
    pedestal = 0.0L;
    ovev = 0l;
    for (ev = 0; ev < strip->GetEventsNumber(); ev++)
      if (   SignalOutOfRange(signal[ev]) 
          || TMath::Abs(signal[ev] - strip->GetPedestal()) > (fPedestalThresholdFactor * strip->GetNoise())) ovev += 1;
      else pedestal = ((ev - ovev)) ? pedestal + (signal[ev] - pedestal) / static_cast<Double_t>(ev - ovev + 1) : signal[ev];
    if (ev == ovev) pedestal = AliITSChannelDaSSD::GetUndefinedValue();
    strip->SetPedestal(static_cast<Float_t>(pedestal));	
    strip->SetOverflowNumber(ovev);	
  }
  return kTRUE;
}


//______________________________________________________________________________
Bool_t AliITSHandleDaSSD::CalculateNoise(const AliITSModuleDaSSD *const module)
{
// Calculates Noise
  AliITSChannelDaSSD *strip;
  Short_t    *signal;
  Float_t     noise;
  Long_t      ovev, n;
  if (!module) return kFALSE;
  for (Int_t strind = 0; strind < module->GetNumberOfStrips(); strind++) {
    if (!(strip = module->GetStrip(strind))) continue;
    if (!(signal = strip->GetSignal())) {
      strip->SetNoise(AliITSChannelDaSSD::GetUndefinedValue());
      AliError(Form("AliITSHandleDaSSD: Error CalculateNoise(): there are no events data for module[%i] strip[%i]->GetSignal()",
                     module->GetModuleId(), strind));
      continue;
    }
    Double_t nsum = 0.0L;
    ovev = 0l;
    for (Long_t ev = 0; ev < strip->GetEventsNumber(); ev++) {
      if (SignalOutOfRange(signal[ev])) ovev += 1;
      else nsum += pow((signal[ev] - strip->GetPedestal()), 2);
    } 
    if ((n = strip->GetEventsNumber() - ovev - 1) > 0) noise =  static_cast<Float_t>(sqrt(nsum / (Float_t)(n)));
    else  noise = AliITSChannelDaSSD::GetUndefinedValue();
    strip->SetNoise(noise);
  }
  return kTRUE;
}



//______________________________________________________________________________
Bool_t AliITSHandleDaSSD::CalculateNoiseCM(AliITSModuleDaSSD *const module)
{
// Calculates Noise with CM correction
  AliITSChannelDaSSD  *strip = NULL;
  Short_t     *signal;
  Float_t      noise;
  Long_t       ovev, n;
  if (!CalculateCM(module)) { 
    AliError("Error: AliITSHandleDaSSD::CalculateCM() returned kFALSE");
    return kFALSE;
  }  
  for (Int_t strind = 0; strind < module->GetNumberOfStrips(); strind++) {
    if (!(strip = module->GetStrip(strind))) continue; //return kFALSE;
    if (!(signal = strip->GetSignal())) {
      strip->SetNoiseCM(AliITSChannelDaSSD::GetUndefinedValue());
      AliError(Form("AliITSHandleDaSSD: Error CalculateNoiseCM(): there are no events data for module[%i] strip[%i]->GetSignal()", 
                     module->GetModuleId(), strind));
      continue; //return kFALSE;
    }
    Int_t chipind = strind / AliITSModuleDaSSD::GetStripsPerChip();
    Double_t nsum = 0.0L;
    ovev = 0l;
    for (Long_t ev = 0; ev < strip->GetEventsNumber(); ev++) {
      if (SignalOutOfRange(signal[ev])) ovev += 1;
      else nsum += pow((signal[ev] - strip->GetPedestal() - module->GetCM(chipind, ev)), 2);
    } 
    if ((n = strip->GetEventsNumber() - ovev - 1) > 0) noise = static_cast<Float_t>(sqrt(nsum / (Double_t)(n)));
    else  noise = AliITSChannelDaSSD::GetUndefinedValue();
    strip->SetNoiseCM(noise);
  }
  return kTRUE;
}



//______________________________________________________________________________
Bool_t AliITSHandleDaSSD::CalculateCM(AliITSModuleDaSSD *const module)
{
// Calculates CM
  AliITSChannelDaSSD  *strip = NULL;
  Short_t             *signal;
  Long_t               ovstr, n;
  Int_t                stripind;
  module->SetNumberOfChips(AliITSModuleDaSSD::GetChipsPerModuleConst());
  for (Int_t chipind = 0; chipind < module->GetNumberOfChips(); chipind++) {
    stripind = chipind * module->GetStripsPerChip();
    module->GetCM()[chipind].Set(fNumberOfEvents);
    module->GetCM()[chipind].Reset(0.0f);
    for (Long_t ev = 0; ev < fNumberOfEvents; ev++) {
    // calculate firs approximation of CM.
      Double_t cm0 = 0.0L;
      ovstr = 0l;
      for (Int_t strind = stripind; strind < (stripind + AliITSModuleDaSSD::GetStripsPerChip()); strind++) {
        if (!(strip = module->GetStrip(strind))) continue; //return kFALSE; 
        if (!(signal = strip->GetSignal())) {
          AliError(Form("AliITSHandleDaSSD: Error CalculateCM: there are no events data for module[%i] strip[%i]->GetSignal()", 
	                module->GetModuleId(), strind));
          return kFALSE;
        }
        if ((SignalOutOfRange(signal[ev])) || (strip->GetPedestal() == AliITSChannelDaSSD::GetUndefinedValue())) ovstr += 1;
        else cm0 += (signal[ev] - strip->GetPedestal());
      }
      if ((n = AliITSModuleDaSSD::GetStripsPerChip() - ovstr)) cm0 /= (Float_t)(n);
      else { module->SetCM(0.0f, chipind, ev); continue; }
    // calculate avarage (cm - (signal - pedestal)) over the chip
      Double_t cmsigma = 0.0L;
      ovstr = 0l;
      for (Int_t strind = stripind; strind < (stripind + AliITSModuleDaSSD::GetStripsPerChip()); strind++) {
        if (!(strip = module->GetStrip(strind))) continue;
        if (!(signal = strip->GetSignal())) {
          AliError(Form("AliITSHandleDaSSD: Error CalculateCM: there are no events data for module[%i] strip[%i]->GetSignal()\n",
	                 module->GetModuleId(), strind));
          return kFALSE;
        }
        if ((SignalOutOfRange(signal[ev])) || (strip->GetPedestal() == AliITSChannelDaSSD::GetUndefinedValue())) ovstr += 1;
        else cmsigma += pow((cm0 - (signal[ev] - strip->GetPedestal())), 2);
      }
      if ((n = AliITSModuleDaSSD::GetStripsPerChip() - ovstr)) cmsigma = sqrt(cmsigma / (Float_t)(n));
      else { module->SetCM(0.0f, chipind, ev); continue; }
   // calculate cm with threshold
      Double_t cmsum = 0.0L;
      ovstr = 0l;
      for (Int_t strind = stripind; strind < (stripind + AliITSModuleDaSSD::GetStripsPerChip()); strind++) {
        if (!(strip = module->GetStrip(strind))) continue;
        signal = strip->GetSignal();
        if ( (SignalOutOfRange(signal[ev])) || (strip->GetPedestal() == AliITSChannelDaSSD::GetUndefinedValue()) 
	       || (TMath::Abs(cm0 - (signal[ev] - strip->GetPedestal())) > (fCmThresholdFactor * cmsigma)) ) ovstr += 1;
        else cmsum += (signal[ev] - strip->GetPedestal());
      }
      if ((n = AliITSModuleDaSSD::GetStripsPerChip() - ovstr)) cmsum /= (Float_t)(n);
      else cmsum = 0.0L;
      if (!(module->SetCM(cmsum, chipind, ev))) 
        AliError(Form("AliITSHandleDaSSD: Error, module->SetCM(...) returned kFALSE module:chip:event : [%d]:[%d]:[%ld]\n",
	                   module->GetModuleId(), chipind, ev));
    } 
  }
  return kTRUE; 
}


//______________________________________________________________________________
Bool_t AliITSHandleDaSSD::ProcessRawData(const Int_t nmread, const Bool_t usewelford)
{
// Performs calculation of calibration parameters (pedestal, noise, ...) 
  Int_t nm = 0;
  if (nmread <= 0) return kFALSE;
  if (!fModules) return kFALSE;
  while ((nm = ReadModuleRawData(nmread)) > 0) {
    cout << "Processing next " << nm << " modules;" << endl;  
    for (Int_t modind = fModIndProcessed; modind < fModIndRead; modind++) {
      if (!fModules[modind]) {
        AliError(Form("AliITSHandleDaSSD: Error ProcessRawData(): No AliITSModuleDaSSD object with index %i is allocated in AliITSHandleDaSSD\n",
	               modind));
        return kFALSE;
      }
      AddFeromCm(fModules[modind]);
      if (usewelford) {
        CalculatePedNoiseW(fModules[modind]);
        CalculateNoiseCMW(fModules[modind]);
      } else {
        CalculatePedestal(fModules[modind]);
        CalculateNoise(fModules[modind]);
        CalculateNoiseCM(fModules[modind]);
      }
    }
    DeleteSignal();
    DeleteCM();
    DeleteCMFerom();
    fModIndProcessed = fModIndRead;
    cout << fModIndProcessed << " - done" << endl;
    if (nm < nmread ) break;
  }
  return kTRUE;  
}


//______________________________________________________________________________
Short_t AliITSHandleDaSSD::RetrieveModuleId(const UChar_t ddlID, const UChar_t ad, const UChar_t adc) const
{
// Retrieve ModuleId from the DDL map which is defined in AliITSRawStreamSSD class
  if (!fDDLModuleMap) {
    AliError("Error DDLMap is not initialized, return 0!");
	return 0;
  }
  Int_t mddli = ((ad - 1) * fgkNumberOfSSDModulesPerSlot) + (adc < 6 ? adc : (adc - 2));
  if ((ddlID < fgkNumberOfSSDDDLs) && (mddli < fgkNumberOfSSDModulesPerDdl)) {
    mddli = fDDLModuleMap[ddlID * fgkNumberOfSSDModulesPerDdl + mddli];
  }  
  else {
    AliWarning(Form("Module index  = %d or ddlID = %d is out of range 0 is rturned", ddlID, mddli));
    mddli = 0;
  }
  if (mddli > SHRT_MAX) return SHRT_MAX;
  else return (Short_t)mddli;
}



//______________________________________________________________________________
AliITSNoiseSSDv2* AliITSHandleDaSSD::GetCalibrationOCDBNoise() const
{
// Fill in the array for OCDB 
  AliITSNoiseSSDv2         *ldcn = NULL;
  AliITSModuleDaSSD      *module = NULL;
  AliITSChannelDaSSD     *strip = NULL; 
  if (!fModules) return NULL;
  ldcn = new AliITSNoiseSSDv2;
  if (!ldcn) {
    AliError("Error allocation mamory for AliITSBadChannelsSSDv2 object, return NULL!");
    return NULL;
  }
  for (Int_t i = 0; i < fNumberOfModules; i++) {
    if (!(module = fModules[i])) continue;
    if (module->GetModuleId() < fgkMinSSDModuleId) continue;
    for (Int_t strind = 0; strind < module->GetNumberOfStrips(); strind++) {
      if (!(strip = module->GetStrip(strind))) continue;
      Short_t modid = module->GetModuleId() - fgkMinSSDModuleId;
      if (strip->GetStripId() < AliITSModuleDaSSD::GetPNStripsPerModule() ) 
        ldcn->AddNoiseP(modid, strip->GetStripId(), strip->GetNoiseCM());
      else
        ldcn->AddNoiseN(modid, (AliITSChannelDaSSD::GetMaxStripIdConst() - strip->GetStripId()), strip->GetNoiseCM());
    }
  }
  return ldcn;
}


//______________________________________________________________________________
AliITSBadChannelsSSDv2* AliITSHandleDaSSD::GetCalibrationBadChannels() const
{
// Fill in the TObjArray with the list of bad channels 
  AliITSBadChannelsSSDv2   *ldcbc = NULL;
  AliITSModuleDaSSD      *module = NULL;
  AliITSChannelDaSSD     *strip = NULL; 
  if (!fModules) return NULL;
  ldcbc = new AliITSBadChannelsSSDv2;
  if (!ldcbc) {
    AliError("Error allocation mamory for AliITSBadChannelsSSDv2 object, return NULL!");
    return NULL;
  }
  for (Int_t i = 0; i < fNumberOfModules; i++) {
    if (!(module = fModules[i])) continue;
    if (module->GetModuleId() < fgkMinSSDModuleId) continue;
    for (Int_t strind = 0; strind < module->GetNumberOfStrips(); strind++) {
      if (!(strip = module->GetStrip(strind))) continue;
      Short_t modid = module->GetModuleId() - fgkMinSSDModuleId;
      if (strip->GetStripId() < AliITSModuleDaSSD::GetPNStripsPerModule() )
        ldcbc->AddBadChannelP(modid, strip->GetStripId(), EvaluateIfChannelIsBad(module, strip->GetStripId()));
      else 
        ldcbc->AddBadChannelN(modid, (AliITSChannelDaSSD::GetMaxStripIdConst() - strip->GetStripId()), 
                                     EvaluateIfChannelIsBad(module, strip->GetStripId()));
    }
  }
  return ldcbc;
}



//______________________________________________________________________________
Bool_t AliITSHandleDaSSD::SaveCalibrationSSDLDC(Char_t*& dafname)
{
// Save Calibration data locally
  AliITSBadChannelsSSDv2   *ldcbc = NULL;
  AliITSPedestalSSDv2      *ldcp = NULL;
  AliITSNoiseSSDv2         *ldcn = NULL;
  AliITSModuleDaSSD      *module = NULL;
  AliITSChannelDaSSD     *strip = NULL; 
  Char_t         *tmpfname;
  TString         dadatafilename("");
  if (!fModules) return kFALSE;
  ldcn = new AliITSNoiseSSDv2;
  ldcp = new AliITSPedestalSSDv2;
  ldcbc = new AliITSBadChannelsSSDv2;
  if ((!ldcn) || (!ldcp) || (!ldcp)) {
    AliError("Error allocation mamory for calibration objects, return kFALSE!");
    return kFALSE;
  }
  for (Int_t i = 0; i < fNumberOfModules; i++) {
    if (!(module = fModules[i])) continue;
    if (module->GetModuleId() < fgkMinSSDModuleId) continue;
    for (Int_t strind = 0; strind < module->GetNumberOfStrips(); strind++) {
      if (!(strip = module->GetStrip(strind))) continue;
      Short_t modid = module->GetModuleId() - fgkMinSSDModuleId;
      if (strip->GetStripId() < AliITSModuleDaSSD::GetPNStripsPerModule() ) {
        ldcn->AddNoiseP(modid, strip->GetStripId(), strip->GetNoiseCM());
        ldcp->AddPedestalP(modid, strip->GetStripId(), strip->GetPedestal());
        ldcbc->AddBadChannelP(modid, strip->GetStripId(), EvaluateIfChannelIsBad(module, strip->GetStripId()));
      } else {
        ldcn->AddNoiseN(modid, (AliITSChannelDaSSD::GetMaxStripIdConst() - strip->GetStripId()), strip->GetNoiseCM());
        ldcp->AddPedestalN(modid, (AliITSChannelDaSSD::GetMaxStripIdConst() - strip->GetStripId()), strip->GetPedestal()); 
        ldcbc->AddBadChannelN(modid, (AliITSChannelDaSSD::GetMaxStripIdConst() - strip->GetStripId()), 
                                     EvaluateIfChannelIsBad(module, strip->GetStripId()));
      }
    }
  }
  if (dafname) dadatafilename.Form("%s/", dafname);
  dadatafilename += TString::Format("ITSSSDda_%i.root", fLdcId);
  tmpfname = new Char_t[dadatafilename.Length()+1];
  Int_t sz = dadatafilename.Sizeof();
  dafname = strncpy(tmpfname, dadatafilename.Data(),sz);
  TFile *fileRun = new TFile (dadatafilename.Data(),"RECREATE");
  if (fileRun->IsZombie()) {
    AliError(Form("AliITSHandleDaSSD: SaveCalibrationSSDLDC() error open file %s", dadatafilename.Data()));
    ldcn->Delete();
    delete fileRun;
    delete ldcn;
    delete ldcp;
    delete ldcbc;
    return kFALSE;
  }
  fileRun->WriteTObject(ldcn);
  fileRun->WriteTObject(ldcp);
  if (fBadChannelsList) 
    if (fMergeBCLists) {
      MergeBadChannels(ldcbc);
      fileRun->WriteTObject(ldcbc);
    } else fileRun->WriteTObject(fBadChannelsList);
  else fileRun->WriteTObject(ldcbc);
  fileRun->Close();
  delete fileRun;
  delete ldcn;
  delete ldcp;
  delete ldcbc;
  return kTRUE;
}


//______________________________________________________________________________
Int_t AliITSHandleDaSSD::MergeBadChannels(AliITSBadChannelsSSDv2*&  bcl) const
{
// Merges the statick bad channels list with bad channels got upon calibration
  AliITSModuleDaSSD     *module = 0;
  Int_t                  nmpch = 0, nmnch = 0, ngpch = 0, ngnch = 0;
  if (!fBadChannelsList || !bcl) {
    AliWarning("Either fBadChannelsList == NULL or bad_channels_list == NULL, there is nothing to merge!");
	return -1;
  }
  for (Int_t modind = 0; modind < GetNumberOfModules(); modind++) {
    if (!(module = fModules[modind])) continue;
    if (module->GetModuleId() < fgkMinSSDModuleId) continue;
    Short_t modid = module->GetModuleId() - fgkMinSSDModuleId;
    for (Int_t strind = 0; strind < AliITSModuleDaSSD::GetPNStripsPerModule(); strind++) {
      if ( (!(fBadChannelsList->GetBadChannelP(modid, strind) & fgkBadChannelMask)) 
          && (bcl->GetBadChannelP(modid, strind) & fgkBadChannelMask) ) 
        ngpch++;
      if ( (!(fBadChannelsList->GetBadChannelN(modid, strind) & fgkBadChannelMask)) 
          && (bcl->GetBadChannelN(modid, strind) & fgkBadChannelMask) ) 
        ngnch++;
      if ( (!(bcl->GetBadChannelP(modid, strind) & fgkBadChannelMask)) 
          && (fBadChannelsList->GetBadChannelP(modid, strind) & fgkBadChannelMask) ) {
        bcl->AddBadChannelP(modid, strind, fBadChannelsList->GetBadChannelP(modid, strind));
        nmpch++;  
      }    
      if ( (!(bcl->GetBadChannelN(modid, strind) & fgkBadChannelMask)) 
          && (fBadChannelsList->GetBadChannelN(modid, strind) & fgkBadChannelMask) ) {
        bcl->AddBadChannelN(modid, strind, fBadChannelsList->GetBadChannelN(modid, strind));
        nmnch++;  
      }
    }	  
  }
  AliInfo(Form("Static bad, dynamic good: P%d,  N%d", nmpch, nmnch));
  AliInfo(Form("Static good, dynamic bad: P%d,  N%d", ngpch, ngnch));
  return (nmnch + nmpch);
}



//______________________________________________________________________________
Bool_t AliITSHandleDaSSD::DumpModInfo(const Float_t meannosethreshold) const
{
// Dump calibration parameters
  AliITSModuleDaSSD    *mod;
  AliITSChannelDaSSD   *strip;
  if (!fModules) {
    cout << "SSDDALDC::DumpModInfo():  Error, no modules are allocated!" << endl;
    return kFALSE;
  }  
  cout << "Modules with MeanNoise > " << meannosethreshold << endl;
  for (Int_t i = 0; i < fNumberOfModules; i++) {
    if (!(mod = fModules[i])) continue;
    Double_t  maxnoise = 0.0L, meannoise = 0.0L, meanovf = 0.0L;
    Float_t   maxped = 0.0f;
    Int_t     maxstrind = 0, novfstr = 0, maxovf = 0; 
    for (Int_t strind = 0; strind < mod->GetNumberOfStrips(); strind++) {
      if (!(strip = mod->GetStrip(strind))) {novfstr++;  continue; }
      if (strip->GetNoiseCM() >= AliITSChannelDaSSD::GetUndefinedValue() ) {novfstr++;  continue; }
      if (maxnoise < strip->GetNoiseCM()) { 
        maxnoise = strip->GetNoiseCM();
        maxstrind = strind;
      }	
      meannoise = (strind - novfstr) ? meannoise + (strip->GetNoiseCM() - meannoise) / static_cast<Double_t>(strind - novfstr + 1) 
                                    : strip->GetNoiseCM();
      if (TMath::Abs(maxped) < TMath::Abs(strip->GetPedestal())) maxped = strip->GetPedestal();
      meanovf = (strind - novfstr) ? meanovf + (strip->GetOverflowNumber() - meanovf) / static_cast<Double_t>(strind - novfstr + 1) 
                                    : strip->GetOverflowNumber();
      if (strip->GetOverflowNumber() > maxovf) maxovf = strip->GetOverflowNumber();
    } 
    if (meannoise > meannosethreshold)
      cout << "Mod: " << i << ";  DDl: " << (int)mod->GetDdlId() << ";  AD: " << (int)mod->GetAD()  
                           << ";  ADC: " << (int)mod->GetADC() << "; MaxPed = " << maxped
			   << ";  MeanNoise = " << meannoise 
			   << ";  NOfStrips = " << (mod->GetNumberOfStrips() - novfstr) << endl;
	if (maxovf > 10) cout << "Max number of events with overflow :  " << maxovf << ";  mean : " << meanovf << endl;
  }
  return kTRUE;
}



//______________________________________________________________________________
Bool_t AliITSHandleDaSSD::PrintModCalibrationData(const UChar_t ddlID, const UChar_t ad, const UChar_t adc, const Char_t *fname) const
{
// Print Module calibration data whether in file or in cout
   AliITSChannelDaSSD  *strip;
   ofstream             datafile;
   ostream             *outputfile;
   if (!fname) { outputfile = &cout; }
   else {
     datafile.open(fname, ios::out);
     if  (datafile.fail()) return kFALSE;
     outputfile = dynamic_cast<ostream*>(&datafile); 
   }
   *outputfile  << "DDL = " << (int)ddlID << ";  AD = " << (int)ad << ";  ADC = " << (int)adc << endl;
   for (Int_t strind = 0; strind < AliITSModuleDaSSD::GetStripsPerModuleConst(); strind++) {
     if ( (strip = GetStrip(ddlID, ad, adc, strind)) ) {
       *outputfile << "Str = " << strind << ";  ped = " << strip->GetPedestal() 
                   << ";  noise = " << strip->GetNoiseCM() << endl;
     }
     else continue;
   }
  if (datafile.is_open()) datafile.close(); 
  return kTRUE;
}



//______________________________________________________________________________
void AliITSHandleDaSSD::DumpInitData(const Char_t *str) const
{
// Print general information retrieved from raw data file
  cout << "Raw data file: " << fRawDataFileName << endl
       << "LDC: " << (Int_t)fLdcId << "; RunId: " << fRunId << endl
       << "Number of physics events: " << fNumberOfEvents << endl
       << str;
}


//______________________________________________________________________________
Bool_t AliITSHandleDaSSD::AllocateSimulatedModules(const Int_t copymodind)
{
// Used to allocate simulated modules to test the performance  
  AliITSModuleDaSSD *module;
  UChar_t      ad, adc, ddlID;
  ad = adc = ddlID = 0;
  if (!(fModules[copymodind])) return kFALSE;
  for (Int_t modind = 0; modind < fNumberOfModules; modind++) {
    if (!fModules[modind]) {
      module = new AliITSModuleDaSSD(AliITSModuleDaSSD::GetStripsPerModuleConst());
      if ((adc < 5) || ((adc > 7) && (adc < 13)) ) adc += 1;
      else if (adc == 5) adc = 8;
      else if (adc == 13) {
        adc = 0;
        if (ad < 8) ad += 1;
        else {
          ad = 0;
          ddlID +=1;
        }
      }
      if (!module->SetModuleIdData (ddlID, ad, adc, modind)) return kFALSE;
      fModules[modind] = module;
      for (Int_t strind = 0; strind < module->GetNumberOfStrips(); strind++) {
        AliITSChannelDaSSD *cstrip = fModules[copymodind]->GetStrip(strind);
	if(!cstrip)return kFALSE;
        Long_t      eventsnumber = cstrip->GetEventsNumber();
        AliITSChannelDaSSD *strip = new AliITSChannelDaSSD(strind, eventsnumber);
        for (Long_t evind = 0; evind < eventsnumber; evind++) {
          Short_t sign = cstrip->GetSignal(evind);
          if (!strip->SetSignal(evind, sign)) AliError(Form("AliITSHandleDaSSD: Copy events error! mod = %i, str = %i", modind, strind));
        }
	module->SetStrip(strip, strind);
      }
    }
    else {
      ddlID = fModules[modind]->GetDdlId();
      ad    = fModules[modind]->GetAD();
      adc   = fModules[modind]->GetADC();
    }
  }
  for (UShort_t modind = 0; modind < fNumberOfModules; modind++) fModules[modind]->SetModuleId(modind + 1080);  
  return kTRUE;
}


    
//___________________________________________________________________________________________
Bool_t AliITSHandleDaSSD::AdDataPresent(const Int_t ddl, const Int_t ad) const
{
// Check if there are calibration data for given ddl and slot
  for (Int_t modind = 0; modind < fNumberOfModules; modind++)
    if ((GetModule(modind)->GetAD() == ad) && (GetModule(modind)->GetDdlId() == ddl)) return kTRUE;
  return kFALSE;
}

   

//___________________________________________________________________________________________
Bool_t AliITSHandleDaSSD::SaveEqSlotCalibrationData(const Int_t ddl, const Int_t ad, const Char_t *fname) const
{
// Saves calibration files for selected equipment (DDL)
  fstream    feefile;
  Int_t      zsml, offsetml;
  ULong_t    zsth, offset, zsoffset;
  if (!fname) {
    AliError("File name must be specified!");
    return kFALSE;   
  }
  if (!AdDataPresent(ddl, ad)) {
    AliError(Form("Error SaveEqSlotCalibrationData(ddl = %i, ad = %i) no data present", ddl, ad)); 
    return kFALSE;
  }
  feefile.open(fname, ios::out);
  if (!feefile.is_open()) {
	AliError(Form("Can not open the file %s for output!", fname)); 
    return kFALSE;
  }
  for (zsml = 0; fgkZsBitMask >> zsml; zsml++) ;
  for (offsetml = 0; fgkOffSetBitMask >> offsetml; offsetml++) ;
  for (Int_t strind = 0; strind < AliITSModuleDaSSD::GetStripsPerModuleConst(); strind++) {
    for (Int_t adcb = 0; adcb < fgkAdcPerDBlock; adcb++) {
      zsoffset = 0x0;
      for (Int_t j = 0; j < 2; j++) {
        Int_t adc = adcb + j * 8;
        zsth = ZsThreshold(ddl, ad, adc, strind);
        offset = OffsetValue(ddl, ad, adc, strind);
        zsoffset = zsoffset | (((zsth << offsetml) | offset) << ((j == 0) ? (offsetml + zsml) : 0));
      }
      feefile << "0x" << ConvBase(static_cast<unsigned long>(zsoffset), 16) << endl;
    }
  }
  feefile.close();
  return kTRUE;
}


//______________________________________________________________________________
Int_t AliITSHandleDaSSD::ChannelIsBad(const UChar_t ddl, const UChar_t ad, const UChar_t adc, const Int_t strn) const
{
// Check if the channel is bad
  AliITSModuleDaSSD     *module = NULL;
  Int_t                  modn = -1;
  if (fBadChannelsList && fDDLModuleMap) {
    modn = RetrieveModuleId(ddl, ad, adc);
    if (modn < 0) return -1;
    if (modn < fgkMinSSDModuleId) {
      AliWarning(Form("Module ddl/ad/adc: %d/%d/%d has number %d which is wrong for SSD module %d", ddl, ad, adc, strn, modn));
	  return -1;
    }
    Short_t modid = modn - fgkMinSSDModuleId;
    if (strn < AliITSModuleDaSSD::GetPNStripsPerModule()) 
      return (fBadChannelsList->GetBadChannelP(modid, strn)  & fgkBadChannelMask);
    else return (fBadChannelsList->GetBadChannelN(modid, (AliITSChannelDaSSD::GetMaxStripIdConst() - strn)) & fgkBadChannelMask);
  } else {
    AliError("Error ether bad channels list or DDLMap is not initialized or both, EvaluateIfChannelIsBad(module, strip) is used!");
    if ((module = GetModule(ddl, ad, adc))) {
      return (EvaluateIfChannelIsBad(module, strn) & fgkBadChannelMask);
    } else {
      AliWarning(Form("There is no calibration data for ddl = %i,  ad = %i,  adc = %i, 0 is used!", ddl, ad, adc));
      return 0ul;
    }  
	return 0;
  }
}



//______________________________________________________________________________
Int_t AliITSHandleDaSSD::LadderIsOff(const UChar_t ddl, const UChar_t ad, const UChar_t adc) const
{
//Checks if the module with given ddl, ad, adc is on the ladder which is in the list of ladders which are off
  const Int_t nm5 =  500;
  const Int_t nm6 =  1248;
  const Int_t nml5a = 12;
  const Int_t nml5c = 10;
  const Int_t nml6a = 12;
  const Int_t nml6c = 13;
  Int_t               modn, ladder, layer, side;
  AliITSModuleDaSSD  *module;
  if (!(module = GetModule(ddl, ad, adc))) return 0;
  if ((modn = module->GetModuleId()) <= 0)  modn = RetrieveModuleId(ddl, ad, adc);
  if (modn <= 0) return 0;
  layer = modn >= nm6 ? 1 : 0;     // 6 : 5
  ladder = (modn - (layer ? nm6 : nm5)) / (layer ? (nml6a + nml6c) : (nml5a + nml5c));
  if ( ((modn - (layer ? nm6 : nm5)) % (layer ? (nml6a + nml6c) : (nml5a + nml5c))) < (layer ? nml6a : nml5a))
    side = 0;      // A
  else side = 1;   // C
  ladder += (layer ? 600 : 500);
  layer += 5;
  if (side)
    if (fCLaddersOff.GetSize()) {
      for(Int_t i = 0; i < fCLaddersOff.GetSize(); i++) 
        if (fCLaddersOff.At(i) == ladder) return fCLaddersOff.At(i);
      return 0;
    } else return 0;
  else
    if (fALaddersOff.GetSize()) {
      for(Int_t i = 0; i < fALaddersOff.GetSize(); i++) 
        if (fALaddersOff.At(i) == ladder) return fALaddersOff.At(i);
      return 0;
    } else return 0;
  return 0;  
}


        
//______________________________________________________________________________
ULong_t AliITSHandleDaSSD::OffsetValue(const AliITSChannelDaSSD *strip, 
                                       const UChar_t ddl, const UChar_t ad, const UChar_t adc, const Int_t strn) const
{
// Calculate the offset value to be upload to FEROM	
  Int_t  pedint;
  if (fOffsetDefault < INT_MAX) pedint = fOffsetDefault;
  else pedint = TMath::Nint(strip->GetPedestal());
  if (pedint > static_cast<Int_t>((fgkOffSetBitMask >> 1))) {
    if (!ChannelIsBad(ddl, ad, adc, strn) && !((fMergeBCLists) && (EvaluateIfChannelIsBad(GetModule(ddl, ad, adc), strn)))) 
      AliError(Form("Offset %i, channel(ddl/ad/adc/strip) %i/%i/%i/%i  can not be represented with mask 0x%s, Offset = %i",
                   pedint, ddl, ad, adc, strn, ConvBase(fgkOffSetBitMask, 16).c_str(), (fgkOffSetBitMask >> 1)));
    return (fgkOffSetBitMask >> 1);
  }  
  if ((-pedint) > static_cast<Int_t>(((fgkOffSetBitMask + 1) >> 1))) {
    if (!ChannelIsBad(ddl, ad, adc, strn) && !((fMergeBCLists) && (EvaluateIfChannelIsBad(GetModule(ddl, ad, adc), strn)))) 
      AliError(Form("Offset %i, channel(ddl/ad/adc/strip) %i/%i/%i/%i  can not be represented with mask 0x%s, Offset = %i", 
               pedint, ddl, ad, adc, strn, ConvBase(fgkOffSetBitMask, 16).c_str(), 
               ((fgkOffSetBitMask & (~fgkOffSetBitMask >> 1)) - fgkOffSetBitMask - 1)));
    return fgkOffSetBitMask & (~fgkOffSetBitMask >> 1);
  }
  return fgkOffSetBitMask & (pedint >= 0 ? pedint : pedint + fgkOffSetBitMask + 1);
}



//______________________________________________________________________________
ULong_t AliITSHandleDaSSD::OffsetValue(const UChar_t ddl, const UChar_t ad, const UChar_t adc, const Int_t strn) const
{
// Calculate the offset value to be upload to FEROM	
  AliITSChannelDaSSD    *strip = NULL;
  AliITSModuleDaSSD     *module = NULL;
  if ((module = GetModule(ddl, ad, adc))) {
    if ((strip = module->GetStrip(strn))) return OffsetValue(strip, ddl, ad, adc, strn);
    else {
      AliWarning(Form("There is no calibration data for ddl = %i,  ad = %i,  adc = %i,  strip = %i, 0 is used!", ddl, ad, adc, strn));
      return 0ul;
    }
  } else {
    AliWarning(Form("There is no calibration data for ddl = %i,  ad = %i,  adc = %i, 0 is used!", ddl, ad, adc));
    return 0ul;
  }  
}



//______________________________________________________________________________
ULong_t AliITSHandleDaSSD::ZsThreshold(const AliITSChannelDaSSD *strip) const
{ 
// Calculate the value of zero suppression threshold to be upload to FEROM
  ULong_t zs;
  if (fZsDefault < 0) {
    zs = TMath::Nint(fZsFactor * strip->GetNoiseCM());
    if (zs < static_cast<ULong_t>(fZsMinimum)) zs = static_cast<ULong_t>(fZsMinimum);
  }
  else zs = fZsDefault;
  return (zs < fgkZsBitMask) ? (zs & fgkZsBitMask) : fgkZsBitMask;
}


//______________________________________________________________________________
ULong_t AliITSHandleDaSSD::ZsThreshold(const UChar_t ddl, const UChar_t ad, const UChar_t adc, const Int_t strn) const
{
// Calculate the value of zero suppression threshold to be upload to FEROM, account bad channels list
  AliITSChannelDaSSD    *strip = NULL;
  AliITSModuleDaSSD     *module = NULL;
  if (ChannelIsBad(ddl, ad, adc, strn)) return fgkZsBitMask;
  if (LadderIsOff(ddl, ad, adc)) return fgkZsBitMask;
  if (fZsDefault > 0) if (static_cast<ULong_t>(fZsDefault) >= fgkZsBitMask) return fgkZsBitMask;
  if ((module = GetModule(ddl, ad, adc))) {
	if (fMergeBCLists) if (EvaluateIfChannelIsBad(module, strn)) return fgkZsBitMask;
    if ((strip = module->GetStrip(strn)))  return ZsThreshold(strip);
    else {
      AliWarning(Form("There is no calibration data for ddl = %i,  ad = %i,  adc = %i,  strip = %i, 0 is used!", ddl, ad, adc, strn));
      return 0ul;
    }
  } else {
    AliWarning(Form("There is no calibration data for ddl = %i,  ad = %i,  adc = %i, 0 is used!", ddl, ad, adc));
    return 0ul;
  }
}

            
//______________________________________________________________________________
string AliITSHandleDaSSD::ConvBase(const unsigned long value, const long base) const
{
// Converts the unsigned long number into that in another base 
  string digits = "0123456789ABCDEF";
  string result;
  unsigned long v = value;
  if((base < 2) || (base > 16)) {
    result = "Error: base out of range.";
  }
  else {
    int i = 0;
    do {
      result = digits[v % base] + result;
      v /= base;
      i++;
    }
    while((v) || (i<8));
  }
  return result;
}



//______________________________________________________________________________
Int_t AliITSHandleDaSSD::CheckOffChips() const
{
// Check if the chip, module are off
  AliITSChannelDaSSD *strip;
  Int_t       offthreshold;
  Int_t       strnd, chipnd, modnd, stroff, chipoff, modoff;
  offthreshold = TMath::Nint(fZsMinimum/fZsFactor);
  modnd = modoff = 0;
  for (Int_t mi = 0; mi < fNumberOfModules; mi++) {
    if (!fModules[mi]) { modnd++; continue; }
    if (fModules[mi]->GetModuleId() < 0) continue;
    if (LadderIsOff(fModules[mi]->GetDdlId(), fModules[mi]->GetAD(), fModules[mi]->GetADC()) ) continue;
    chipoff = chipnd = 0;
    for (Int_t chipind = 0; chipind < AliITSModuleDaSSD::GetChipsPerModuleConst(); chipind++) {
      strnd = stroff = 0;
      Int_t stripind = chipind * AliITSModuleDaSSD::GetStripsPerChip();
      for (Int_t strind = stripind; strind < (stripind + AliITSModuleDaSSD::GetStripsPerChip()); strind++) {
        if (!(strip = fModules[mi]->GetStrip(strind))) { strnd++;  continue; }
        if (strip->GetNoiseCM() < offthreshold ) stroff++;
      }
      if (strnd == AliITSModuleDaSSD::GetStripsPerChip()) chipnd++;
      else if (stroff == AliITSModuleDaSSD::GetStripsPerChip()) chipoff++;
      else if ((stroff + strnd) == AliITSModuleDaSSD::GetStripsPerChip()) chipoff++;
    }
    if ((!chipoff) && (!chipnd)) continue;
    if (chipnd == AliITSModuleDaSSD::GetChipsPerModuleConst()) {
      AliInfo(Form("Module: (ddl/ad/adc) %i/%i/%i seems to be off and it is not on the ladders which are off!", 
               fModules[mi]->GetDdlId(), fModules[mi]->GetAD(), fModules[mi]->GetADC()));
      modnd++;
    }
    if (chipoff == AliITSModuleDaSSD::GetChipsPerModuleConst()) {
      AliInfo(Form("Module (ddl/ad/adc): %i/%i/%i seems to be off and it is not on the ladders which are off!", 
               fModules[mi]->GetDdlId(), fModules[mi]->GetAD(), fModules[mi]->GetADC()));
      modoff++;
    }
    else if ((chipoff + chipnd) == AliITSModuleDaSSD::GetChipsPerModuleConst()) {
      AliInfo(Form("Module: (ddl/ad/adc): %i/%i/%i seems to be off and it is not on the ladders which are off!", 
               fModules[mi]->GetDdlId(), fModules[mi]->GetAD(), fModules[mi]->GetADC()));
      modoff++;
    }
    else if (chipoff) {
      AliInfo(Form("Module: (ddl/ad/adc): %i/%i/%i has %i chips which are off!", 
               fModules[mi]->GetDdlId(), fModules[mi]->GetAD(), fModules[mi]->GetADC(), chipoff));
      modoff++;
    }
  }
  return (modoff + modnd);
}


//______________________________________________________________________________
Bool_t AliITSHandleDaSSD::CalculatePedNoiseW(const AliITSModuleDaSSD *const module)
{
// Calculates Pedestal and Noise using Welford algorithm
  AliITSChannelDaSSD *strip;
  Double_t            pedestal, noise, p0, s0;
  Short_t            *signal;
  Int_t               ovev, n;
  if (!module) return kFALSE;
  for (Int_t strind = 0; strind < module->GetNumberOfStrips(); strind++) {
    if (!(strip = module->GetStrip(strind))) continue;
    if (!(signal = strip->GetSignal())) {
      AliError(Form("AliITSHandleDaSSD: Error CalculatePedestal(): there are no events data for module[%i] strip[%i]->GetSignal()",
                     module->GetModuleId(), strind));
      continue;	
    }
//************* pedestal and noise first pass ****************
    pedestal = p0 = noise = 0.0L;
    ovev = 0;
    for (Int_t ev = 0; ev < strip->GetEventsNumber(); ev++) {
      if (SignalOutOfRange(signal[ev])) ovev += 1;
      else 
        if (!(ev - ovev)) {
          pedestal = p0 = signal[ev];
          noise = 0.0L;
        } else {
          p0 = pedestal + (signal[ev] - pedestal) / static_cast<Double_t>(ev - ovev + 1);
          s0 = noise + (signal[ev] - pedestal) * (signal[ev] - p0);
          pedestal = p0;
          noise = s0;
        }
    }
    if (strip->GetEventsNumber() == ovev) pedestal = AliITSChannelDaSSD::GetUndefinedValue();
    strip->SetPedestal(static_cast<Float_t>(pedestal));
    if ((n = strip->GetEventsNumber() - ovev - 1) > 0) 
      strip->SetNoise( static_cast<Float_t>(sqrt(noise / static_cast<Double_t>(n))) );
    else {
      strip->SetNoise(AliITSChannelDaSSD::GetUndefinedValue());
      continue;
    }
//************* Second pass excluds event with |p - s|>f*noise *****************
    pedestal = p0 = noise = 0.0L;
    ovev = 0;
    for (Int_t ev = 0; ev < strip->GetEventsNumber(); ev++) {
      if (   SignalOutOfRange(signal[ev]) 
          || TMath::Abs(signal[ev] - strip->GetPedestal()) > (fPedestalThresholdFactor * strip->GetNoise())) ovev += 1;
      else
        if (!(ev - ovev)) {
          pedestal = p0 = signal[ev];
          noise = 0.0L;
        } else {
          p0 = pedestal + (signal[ev] - pedestal) / static_cast<Double_t>(ev - ovev + 1);
          s0 = noise + (signal[ev] - pedestal) * (signal[ev] - p0);
          pedestal = p0;
          noise = s0;
        }
    }      
    if ((n = strip->GetEventsNumber() - ovev - 1) > 0) noise = sqrt(noise / static_cast<Double_t>(n));
    else  noise = AliITSChannelDaSSD::GetUndefinedValue();
    strip->SetNoise(static_cast<Float_t>(noise));
    if (strip->GetEventsNumber() == ovev) pedestal = AliITSChannelDaSSD::GetUndefinedValue();
    strip->SetPedestal(static_cast<Float_t>(pedestal));
    strip->SetOverflowNumber(ovev);	
  }
  return kTRUE;
}


//______________________________________________________________________________
Bool_t AliITSHandleDaSSD::CalculateCMW(AliITSModuleDaSSD *const module)
{
// Calculates CM using Welford algorithm
  AliITSChannelDaSSD  *strip = NULL;
  Short_t             *signal;
  Int_t                ovstr, n;
  Int_t                stripind;
  Double_t             cm0, cm1, cmsigma, cms1;
  module->SetNumberOfChips(AliITSModuleDaSSD::GetChipsPerModuleConst());
  for (Int_t chipind = 0; chipind < module->GetNumberOfChips(); chipind++) {
    stripind = chipind * module->GetStripsPerChip();
    module->GetCM()[chipind].Set(fNumberOfEvents);
    module->GetCM()[chipind].Reset(0.0f);
    for (Long_t ev = 0; ev < fNumberOfEvents; ev++) {
    // calculate firs approximation of CM and SigmaCM.
      cm0 = cm1 = cmsigma = 0.0L;
      ovstr = 0;
      for (Int_t strind = stripind; strind < (stripind + AliITSModuleDaSSD::GetStripsPerChip()); strind++) {
        if (!(strip = module->GetStrip(strind))) { ovstr += 1; continue; } //return kFALSE; 
        if (!(signal = strip->GetSignal())) { ovstr += 1; continue; }  //return kFALSE; 
        if ((SignalOutOfRange(signal[ev])) || (strip->GetPedestal() == AliITSChannelDaSSD::GetUndefinedValue())) ovstr += 1;
        else {
          if (!(strind - stripind - ovstr)) {
          cm0 = cm1 = signal[ev] - strip->GetPedestal();
          cmsigma = 0.0L;
        } else {
          cm1 = cm0 + (signal[ev] - strip->GetPedestal() - cm0) / static_cast<Double_t>(strind - stripind - ovstr + 1);
          cms1 = cmsigma + (signal[ev] - strip->GetPedestal() - cm0) * (signal[ev] - strip->GetPedestal() - cm1);
          cm0 = cm1;
          cmsigma = cms1;
        } }
      }
      if ((n = AliITSModuleDaSSD::GetStripsPerChip() - ovstr - 1) > 0) cmsigma = sqrt(cmsigma / static_cast<Double_t>(n));
      else {
        AliWarning(Form("AliITSHandleDaSSD: Too little number of strips have a signal for module:chip:event : [%d]:[%d]:[%ld]\n",
                   module->GetModuleId(), chipind, ev));
        if (!(module->SetCM(0.0f, chipind, ev))) 
          AliError(Form("AliITSHandleDaSSD: Error, module->SetCM(...) returned kFALSE module:chip:event : [%d]:[%d]:[%ld]\n",
                   module->GetModuleId(), chipind, ev));
        continue;
      }
   // calculate cm with threshold
      Double_t cmsum = 0.0L;
      ovstr = 0;
      for (Int_t strind = stripind; strind < (stripind + AliITSModuleDaSSD::GetStripsPerChip()); strind++) {
        if (!(strip = module->GetStrip(strind))) { ovstr += 1; continue; }
        if (!(signal = strip->GetSignal())) { ovstr += 1; continue; }
        if ( (SignalOutOfRange(signal[ev])) || (strip->GetPedestal() == AliITSChannelDaSSD::GetUndefinedValue()) 
	       || (TMath::Abs(cm0 - (signal[ev] - strip->GetPedestal())) > (fCmThresholdFactor * cmsigma)) ) ovstr += 1;
        else cmsum += (signal[ev] - strip->GetPedestal());
      }
      if ((n = AliITSModuleDaSSD::GetStripsPerChip() - ovstr)) cmsum /= (Double_t)(n);
      else cmsum = 0.0L;
      if (!(module->SetCM(static_cast<Float_t>(cmsum), chipind, ev))) 
        AliError(Form("AliITSHandleDaSSD: Error, module->SetCM(...) returned kFALSE module:chip:event : [%d]:[%d]:[%ld]\n",
                 module->GetModuleId(), chipind, ev));
    } 
  }
  return kTRUE; 
}


//______________________________________________________________________________
Bool_t AliITSHandleDaSSD::CalculateNoiseCMW(AliITSModuleDaSSD *const module)
{
// Calculates Noise with CM correction
  AliITSChannelDaSSD  *strip = NULL;
  Short_t     *signal;
  Int_t        ovev, n;
  if (!CalculateCMW(module)) { 
    AliError("Error: AliITSHandleDaSSD::CalculateCMW() returned kFALSE");
    return kFALSE;
  }  
  for (Int_t strind = 0; strind < module->GetNumberOfStrips(); strind++) {
    if (!(strip = module->GetStrip(strind))) continue; //return kFALSE;
    if (!(signal = strip->GetSignal())) {
      strip->SetNoiseCM(AliITSChannelDaSSD::GetUndefinedValue());
      AliError(Form("AliITSHandleDaSSD: Error CalculateNoiseCMW(): there are no events data for module[%i] strip[%i]->GetSignal()", 
                     module->GetModuleId(), strind));
      continue; //return kFALSE;
    }
//** To get exactly the same set of events as for pedestal and noise calculation **
    Double_t pedestal, noise, p0, s0;
    pedestal = p0 = noise = 0.0L;
    ovev = 0;
    for (Int_t ev = 0; ev < strip->GetEventsNumber(); ev++) {
      if (SignalOutOfRange(signal[ev])) ovev += 1;
      else 
        if (!(ev - ovev)) {
          pedestal = p0 = signal[ev];
          noise = 0.0L;
        } else {
          p0 = pedestal + (signal[ev] - pedestal) / static_cast<Double_t>(ev - ovev + 1);
          s0 = noise + (signal[ev] - pedestal) * (signal[ev] - p0);
          pedestal = p0;
          noise = s0;
        }
    }
    if ((n = strip->GetEventsNumber() - ovev - 1) > 0) noise = sqrt(noise / static_cast<Double_t>(n));
    else  {
	  strip->SetNoiseCM(AliITSChannelDaSSD::GetUndefinedValue());
	  continue;
    }
//** Calculation of CM corrected noise **
    Int_t chipind = strind / AliITSModuleDaSSD::GetStripsPerChip();
    Double_t nsum = 0.0L;
    ovev = 0;
    for (Int_t ev = 0; ev < strip->GetEventsNumber(); ev++) {
      if (   SignalOutOfRange(signal[ev])
          || TMath::Abs(signal[ev] - pedestal) > (fPedestalThresholdFactor * noise)) ovev += 1;
      else nsum += pow((signal[ev] - strip->GetPedestal() - module->GetCM(chipind, ev)), 2);
    } 
    if ((n = strip->GetEventsNumber() - ovev - 1) > 0) noise =  sqrt(nsum / static_cast<Double_t>(n));
    else  noise = AliITSChannelDaSSD::GetUndefinedValue();
    strip->SetNoiseCM(static_cast<Float_t>(noise));
  }
  return kTRUE;
}



//______________________________________________________________________________
UChar_t AliITSHandleDaSSD::EvaluateIfChannelIsBad(const AliITSModuleDaSSD *const module, const Int_t stripn) const
{
//Applies the bad channel creteria and set the appropriate flags for returned value
  AliITSChannelDaSSD  *strip = 0;
  UInt_t               bcflags = 0;
  if (fZsDefault >= 0) { if (static_cast<ULong_t>(fZsDefault) >= fgkZsBitMask) bcflags |= 3; }
  else if (static_cast<ULong_t>(fZsMinimum) >= fgkZsBitMask) bcflags |= 3;
  if (LadderIsOff(module->GetDdlId(), module->GetAD(), module->GetADC()) )  bcflags |= 3;
  
  if (!(strip = module->GetStrip(stripn))) bcflags |= 3;
  else {
    if (strip->GetNoiseCM() == AliITSChannelDaSSD::GetUndefinedValue()) bcflags |= 8;
    if (static_cast<ULong_t>(TMath::Nint(fZsFactor * strip->GetNoiseCM())) >= fgkZsBitMask) bcflags |= 8;
    if (strip->GetNoiseCM() < 1) bcflags |= 16;
    if (strip->GetPedestal() > ((fgkOffSetBitMask >> 1) - 1))  bcflags |= 4;
    else if ((-(strip->GetPedestal())) > (fgkOffSetBitMask >> 1))  bcflags |= 4;
    if (bcflags) bcflags |= 3;
  }
  return bcflags;
}

