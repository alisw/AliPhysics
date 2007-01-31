/**************************************************************************
 * Copyright(c) 1998-2003, ALICE Experiment at CERN, All rights reserved. *
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

////////////////////////////////////////////////////////////
//This class contains the methods to create raw data      //
//for CTP (trigger). The CTP data format is taken as      //
//described in the Trigger TDR - pages 134 and 135.       //
////////////////////////////////////////////////////////////

#include <TObject.h>
#include <TString.h>
#include <Riostream.h>
#include "AliCTPRawData.h"
#include "AliRunLoader.h"
#include "AliCentralTrigger.h"
#include "AliRawDataHeader.h"
#include "AliLog.h"
#include "AliDAQ.h"

ClassImp(AliCTPRawData)
////////////////////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
AliCTPRawData::AliCTPRawData()
{
  // Default constructor
}

//______________________________________________________________________________
AliCTPRawData::AliCTPRawData(const AliCTPRawData &source):
  TObject(source)
{
  // Copy Constructor
}

//______________________________________________________________________________
AliCTPRawData& AliCTPRawData::operator=(const AliCTPRawData &source)
{
  // Assigment operator
  if(this==&source) return *this;
  ((TObject *)this)->operator=(source);
 
 return *this;
}

//______________________________________________________________________________
void AliCTPRawData::RawData()
{
  // This method writes the CTP (trigger)
  // raw data in a DDL file
  ULong64_t l2class = 0;
  UChar_t l2cluster = 0;
  AliInfo("Storing the CTP DDL raw data...");
  AliRunLoader *runloader = AliRunLoader::GetRunLoader();
  if (runloader) {
    if (!runloader->LoadTrigger()) {
      AliCentralTrigger *aCTP = runloader->GetTrigger();
      if (AliDebugLevel() > 0)
	aCTP->Dump();
      // First get the trigger mask
      l2class = aCTP->GetClassMask();
      // Then get the detector cluster to be read out
      l2cluster = aCTP->GetClusterMask();
    }
    else
      AliWarning("No trigger can be loaded! Putting empty trigger class into the CTP raw data !");
  }
  else
    AliError("No run loader is available! Putting empty trigger class into the CTP raw data !");

  AliDebug(1,Form("CTP trigger mask = 0x%llx",l2class));
  AliDebug(1,Form("CTP detector cluster = 0x%x",l2cluster));

  char  fileName[15];
  strcpy(fileName,AliDAQ::DdlFileName("TRG",0));
  AliInfo(Form("Storing CTP raw data in %s",fileName));
  ofstream outfile;         // logical name of the output file 
 #ifndef __DECCXX
    outfile.open(fileName,ios::binary);
#else
    outfile.open(fileName);
#endif

  AliRawDataHeader header;
  // Write a dummy header
  UInt_t dataHeaderPosition=outfile.tellp();
  outfile.write((char*)(&header),sizeof(header));

  // Writing CTP raw data here
  // The format is taken as in
  // pages 134 and 135 of the
  // Trigger TDR

  // If not 0, from where??
  UInt_t bunchCross = 0;
  UInt_t orbitId = 0;
  Bool_t esr = 0; // Enable Segmented Readout flag

  // First 3 words here - Bunch-crossing and orbit numbers
  UInt_t word = 0;
  word |= 0 << 13; // BlockID = 0 in case of CTP readout
  word |= bunchCross & 0xFFF;
  AliDebug(1,Form("CTP word1 = 0x%x",word));
  outfile.write((char*)(&word),sizeof(UInt_t));

  word = 0;
  word |= 0 << 13; // BlockID = 0 in case of CTP readout
  word |= (orbitId >> 12) & 0xFFF;
  AliDebug(1,Form("CTP word2 = 0x%x",word));
  outfile.write((char*)(&word),sizeof(UInt_t));
  word = 0;
  word |= 0 << 13; // BlockID = 0 in case of CTP readout
  word |= orbitId & 0xFFF;
  AliDebug(1,Form("CTP word3 = 0x%x",word));
  outfile.write((char*)(&word),sizeof(UInt_t));

  // Now the 4th word
  word = 0;
  word |= 0 << 13; // BlockID = 0 in case of CTP readout
  word |= ((UInt_t)esr) << 10;
  word |= 0 << 8; // L2SwC - physics trigger
  word |= (l2cluster & 0x3F) << 2; // L2Cluster
  word |= (UInt_t)((l2class >> 48) & 0x3);
  AliDebug(1,Form("CTP word4 = 0x%x",word));
  outfile.write((char*)(&word),sizeof(UInt_t));

  // Then the last 4 words with the trigger classes
  word = 0;
  word |= 0 << 13; // BlockID = 0 in case of CTP readout
  word |= (UInt_t)((l2class >> 36) & 0xFFF);
  AliDebug(1,Form("CTP word5 = 0x%x",word));
  outfile.write((char*)(&word),sizeof(UInt_t));
  word = 0;
  word |= 0 << 13; // BlockID = 0 in case of CTP readout
  word |= (UInt_t)((l2class >> 24) & 0xFFF);
  AliDebug(1,Form("CTP word6 = 0x%x",word));
  outfile.write((char*)(&word),sizeof(UInt_t));
  word = 0;
  word |= 0 << 13; // BlockID = 0 in case of CTP readout
  word |= (UInt_t)((l2class >> 12) & 0xFFF);
  AliDebug(1,Form("CTP word7 = 0x%x",word));
  outfile.write((char*)(&word),sizeof(UInt_t));
  word = 0;
  word |= 0 << 13; // BlockID = 0 in case of CTP readout
  word |= (UInt_t)(l2class & 0xFFF);
  AliDebug(1,Form("CTP word8 = 0x%x",word));
  outfile.write((char*)(&word),sizeof(UInt_t));

  // Write the real data header
  UInt_t currentFilePosition=outfile.tellp();
  outfile.seekp(dataHeaderPosition);
  header.fSize=currentFilePosition-dataHeaderPosition;
  header.SetAttribute(0);  // valid data
  header.SetTriggerClass(l2class);
  outfile.write((char*)(&header),sizeof(header));
  outfile.close();

  return;
}
