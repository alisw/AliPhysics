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

////////////////////////////////////
//
// MUON Raw Data generator in ALICE-MUON
//
// This class v-1:
// * generates raw data for MUON tracker only (for the moment)
// * a simple mapping is used (see below)
// * the bus patch id is calculated but not stored properly in the raw data
//   this has to be changed
// one DDL per 1/2 chamber is created for both cathode.
//
////////////////////////////////////

#include "AliMUONRawData.h"
#include "AliMUONDigit.h"
#include "AliMUONConstants.h"
#include "AliMUONData.h"
#include "AliMUONDDLTracker.h"
#include "AliRun.h" 
#include "AliRunLoader.h"
#include "AliLoader.h"
#include "AliBitPacking.h" 
#include "AliRawDataHeader.h"

const Int_t AliMUONRawData::fgkDefaultPrintLevel = 0;

ClassImp(AliMUONRawData) // Class implementation in ROOT context

//__________________________________________________________________________
AliMUONRawData::AliMUONRawData(AliLoader* loader)
  : TObject()
{
  // Standard Constructor
 
  fDebug           = 0;
  fNCh             = 0;
  fNTrackingCh     = 0;
  fNTriggerCh      = 0;
  fMUONData        = 0;

  fPrintLevel = fgkDefaultPrintLevel;

  // initialize loader's
  fLoader = loader;

  // initialize container
  fMUONData  = new AliMUONData(fLoader,"MUON","MUON");
}

//__________________________________________________________________________
AliMUONRawData::AliMUONRawData()
  : TObject(),
    fNCh(0),
    fNTrackingCh(0),
    fNTriggerCh(0),
    fMUONData(0),
    fPrintLevel(fgkDefaultPrintLevel),
    fDebug(0),
    fLoader(0)
{
  // Default Constructor
}

//_______________________________________________________________________
AliMUONRawData::AliMUONRawData (const AliMUONRawData& rhs)
  : TObject(rhs)
{
// Protected copy constructor

  Fatal("AliMUONRawData", "Not implemented.");
}

//_______________________________________________________________________
AliMUONRawData & 
AliMUONRawData::operator=(const AliMUONRawData& rhs)
{
// Protected assignement operator

  if (this == &rhs) return *this;

  Fatal("operator=", "Not implemented.");
    
  return *this;  
}

//__________________________________________________________________________
AliMUONRawData::~AliMUONRawData(void)
{
  if (fMUONData)
    delete fMUONData;
  return;
}
//____________________________________________________________________
Int_t AliMUONRawData::WriteRawData()
{
 // convert digits of the current event to raw data


  Int_t DDLId;
  Char_t name[20];

  fLoader->LoadDigits("READ");
  fMUONData->SetTreeAddress("D");
  //  printf("je suis dans WriteRawData\n");

  for (Int_t ich = 0; ich < AliMUONConstants::NTrackingCh(); ich++) {

    // open files
    DDLId = ich * 2  + 0x900;
    sprintf(name, "MUON_%d.ddl",DDLId);
    fFile1 = fopen(name,"w");

    DDLId = (ich * 2) + 1 + 0x900 ;
    sprintf(name, "MUON_%d.ddl",DDLId);
    fFile2 = fopen(name,"w");

    WriteDDL(ich);
  
    // reset and close
    fclose(fFile1);
    fclose(fFile2);
    fMUONData->ResetDigits();
  }
  
  return kTRUE;
}
//____________________________________________________________________
Int_t AliMUONRawData::WriteDDL(Int_t iCh)
{

  TClonesArray* muonDigits = 0;

  // DDL event one per half chamber
  AliMUONDDLTracker* eventDDL1;
  AliMUONDDLTracker* eventDDL2;

  // data format
  Char_t parity = 0x4;
  UShort_t manuId;
  UChar_t channelId;
  UShort_t charge;
  Int_t busPatchId;

  UInt_t word;
  Int_t nWord;
  Int_t dummy;

  Int_t offsetX = 0; // offet row
  Int_t offsetY = 0; // offset columns
  Int_t offsetCath = 0; //offset from one cathod to the other
  Int_t maxChannel = 0; // maximum nb of channel in 1/2 chamber
  Int_t id;

  Int_t nDigits;
  AliMUONDigit* digit;

  AliRawDataHeader header;

  eventDDL1 = new AliMUONDDLTracker();
  eventDDL2 = new AliMUONDDLTracker();

  for (Int_t iCath = 0; iCath < 2; iCath++) {

    if (fPrintLevel)
      printf("WriteDDL chamber %d and cathode %d\n", iCh+1, iCath);

    fMUONData->ResetDigits();
    fMUONData->GetCathode(iCath);
    muonDigits = fMUONData->Digits(iCh);

    nDigits = muonDigits->GetEntriesFast();
    if (fPrintLevel)
      printf("ndigits = %d\n",nDigits);

    // open DDL file, on per 1/2 chamber
 
    for (Int_t idig = 0; idig < nDigits; idig++) {

      digit = (AliMUONDigit*) muonDigits->UncheckedAt(idig);

      switch (iCh+1) {
      case 1:
      case 2:
      case 3:
      case 4:
	offsetX = 512;
	offsetY = 256;
	offsetCath = 65536;
	maxChannel = (offsetY * offsetX + 2* offsetY + offsetCath);
	break;
      case 5:
      case 6:
      case 7:
      case 8:
      case 9:
      case 10:
	offsetX = 1024;
	offsetY = 0;
	offsetCath = 65536;
	maxChannel = (256 * offsetX + offsetX + offsetCath);
	break;
      }
      // dummy mapping
      // manu Id directly from a matrix 8*8, same segmentation for B and NB
      // 50 buspatches for 1/2 chamber
 
      id =  (TMath::Abs(digit->PadX()) * offsetX + digit->PadY() + offsetY +
	     offsetCath*iCath);
      busPatchId = id/50;

      Int_t inBusId =  (id - maxChannel/50 * busPatchId);// id channel in buspatch
      Int_t manuPerBus = (maxChannel/(50*64)); // number of manus per buspatch

      // 64 manu cards for one buspatch
      manuId = inBusId/manuPerBus;
      manuId &= 0x7FF; // 11 bits 

      // channel id
      channelId = (inBusId % manuPerBus);
      channelId &= 0x3F; // 6 bits

      // charge
      charge = digit->Signal();
      charge &= 0xFFF;

      if (fPrintLevel)
	printf("id: %d, busPatchId %d, manuId: %d, channelId: %d, padx: %d pady %d, charge %d\n",
	       id, busPatchId, manuId, channelId, digit->PadX(), digit->PadY(), digit->Signal());
      //packing word
      AliBitPacking::PackWord((UInt_t)parity,word,29,31);
      AliBitPacking::PackWord((UInt_t)manuId,word,18,28);
      AliBitPacking::PackWord((UInt_t)channelId,word,12,17);
      AliBitPacking::PackWord((UInt_t)charge,word,0,11);

      // set  DDL Event
      if (digit->PadX() > 0) {
	eventDDL1->SetRawData(word);
	eventDDL1->SetBusPatchId(busPatchId); // information not usable only last bus patch stored
      } else {
	eventDDL2->SetRawData(word);
	eventDDL2->SetBusPatchId(busPatchId);
      }
      if (fPrintLevel) {
	printf("word: 0x%x, ",word);
	printf("manuId back %d, ",eventDDL1->GetManuId(eventDDL1->GetLength()-1));
	printf("channelId back %d, ",eventDDL1->GetChannelId(eventDDL1->GetLength()-1));
	printf("charge back %d\n",eventDDL1->GetCharge(eventDDL1->GetLength()-1));
      }
    }
  }
  // write DDL event
  nWord = eventDDL1->GetLength();
  if (nWord > 0) { // not empty event
    header.fSize = nWord*4 + sizeof(AliRawDataHeader) + 12; // include EoD & length
    fwrite((char*)(&header),sizeof(header),1,fFile1);
    fwrite(eventDDL1->GetAddress(),sizeof(int),nWord+2,fFile1);
    dummy = eventDDL1->GetEoD(); 
    fwrite(&dummy,sizeof(int),1,fFile1);// could be nicer !
  }
  nWord = eventDDL2->GetLength();
  if (nWord > 0) {
    header.fSize = nWord*4 + sizeof(AliRawDataHeader) + 12;
    fwrite((char*)(&header),sizeof(header),1,fFile2);
    fwrite(eventDDL2->GetAddress(),sizeof(int),nWord+2,fFile2);
    dummy = eventDDL2->GetEoD();
    fwrite(&dummy,sizeof(int),1,fFile2);
  }
  delete eventDDL1;
  delete eventDDL2;
  return kTRUE;
}
