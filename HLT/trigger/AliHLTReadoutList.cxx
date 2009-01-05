/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Artur Szostak <artursz@iafrica.com>                   *
 *                  for The ALICE HLT Project.                            *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/// @file   AliHLTReadoutList.cxx
/// @author Artur Szostak <artursz@iafrica.com>
/// @date   19 Nov 2008
/// @brief  Implementation of the AliHLTReadoutList class.
///
/// The AliHLTReadoutList class is used as an interface to the AliHLTEventDDL
/// structure. It makes it easy to manipulate the bits in this structure, which
/// define what DDLs should be readout by DAQ.
/// Several operators are also overloaded which are meant to be used in the trigger
/// menu specification for the AliHLTGlobalTrigger. It allows one to construct
/// expressions for the readout lists, which is necessary to be able to evaluate
/// or compose the final readout list, given multiple input readout lists received
/// from individual components that derive from AliHLTTrigger.

#include "AliHLTReadoutList.h"
#include "AliDAQ.h"
#include "Riostream.h"
#include "TString.h"
#include "TObjString.h"
#include "TObjArray.h"
#include <cassert>

ClassImp(AliHLTReadoutList)


AliHLTReadoutList::AliHLTReadoutList() :
	TObject(),
	fReadoutList()
{
  // Default constructor.
  
  fReadoutList.fCount = gkAliHLTDDLListSize;  // Required by ALICE-INT-2007-015
  memset(fReadoutList.fList, 0x0, sizeof(fReadoutList.fList));
}


AliHLTReadoutList::AliHLTReadoutList(Int_t enabledDetectors) :
	TObject(),
	fReadoutList()
{
  // Constructor to select which detectors to enable for readout.
  // See header file for more details.
  
  fReadoutList.fCount = gkAliHLTDDLListSize;  // Required by ALICE-INT-2007-015
  memset(fReadoutList.fList, 0x0, sizeof(fReadoutList.fList));
  Enable(enabledDetectors);
}


AliHLTReadoutList::AliHLTReadoutList(const char* enabledList) :
	TObject(),
	fReadoutList()
{
  // Constructor to select which detectors and DDLs to enable for readout.
  // See header file for more details.
  
  fReadoutList.fCount = gkAliHLTDDLListSize;  // Required by ALICE-INT-2007-015
  memset(fReadoutList.fList, 0x0, sizeof(fReadoutList.fList));
  
  TString str(enabledList);
  str.ToUpper();
  Int_t enabledDetectors = 0;
  if (str.Contains("ITSSPD")) enabledDetectors |= kITSSPD;
  if (str.Contains("ITSSDD")) enabledDetectors |= kITSSDD;
  if (str.Contains("ITSSSD")) enabledDetectors |= kITSSSD;
  if (str.Contains("TPC")) enabledDetectors |= kTPC;
  if (str.Contains("TRD")) enabledDetectors |= kTRD;
  if (str.Contains("TOF")) enabledDetectors |= kTOF;
  if (str.Contains("HMPID")) enabledDetectors |= kHMPID;
  if (str.Contains("PHOS")) enabledDetectors |= kPHOS;
  if (str.Contains("CPV")) enabledDetectors |= kCPV;
  if (str.Contains("PMD")) enabledDetectors |= kPMD;
  if (str.Contains("MUONTRK")) enabledDetectors |= kMUONTRK;
  if (str.Contains("MUONTRG")) enabledDetectors |= kMUONTRG;
  if (str.Contains("FMD")) enabledDetectors |= kFMD;
  if (str.Contains("T0")) enabledDetectors |= kT0;
  if (str.Contains("V0")) enabledDetectors |= kV0;
  if (str.Contains("ZDC")) enabledDetectors |= kZDC;
  if (str.Contains("ACORDE")) enabledDetectors |= kACORDE;
  if (str.Contains("TRG")) enabledDetectors |= kTRG;
  if (str.Contains("EMCAL")) enabledDetectors |= kEMCAL;
  if (str.Contains("DAQTEST")) enabledDetectors |= kDAQTEST;
  if (str.Contains("HLT")) enabledDetectors |= kHLT;
  if (str.Contains("ALL")) enabledDetectors |= kALLDET;
  Enable(enabledDetectors);
  
  TObjArray* list = str.Tokenize(" ");
  TIter next(list);
  const TObjString* objstr = NULL;
  while ((objstr = dynamic_cast<const TObjString*>(next())) != NULL)
  {
    str = objstr->GetString();
    if (str.IsDigit()) EnableDDLBit(str.Atoi());
  }
  delete list;
}


AliHLTReadoutList::AliHLTReadoutList(const AliHLTEventDDL& list) :
	TObject(),
	fReadoutList()
{
  // Constructor to create readout list from AliHLTEventDDL structure.
  // See header file for more details.
  
  memcpy(&fReadoutList, &list, sizeof(fReadoutList));
}


AliHLTReadoutList::AliHLTReadoutList(const AliHLTReadoutList& list) :
	TObject(list),
	fReadoutList()
{
  // Copy constructor performs a deep copy.
  
  memcpy(&fReadoutList, &list.fReadoutList, sizeof(fReadoutList));
}


AliHLTReadoutList::~AliHLTReadoutList()
{
  // Default destructor.
}


bool AliHLTReadoutList::DecodeDDLID(Int_t ddlId, Int_t& wordIndex, Int_t& bitIndex)
{
  // Decodes the word index and bit index within that word for the readout list structure.
  // See header file for more details.
  
  // The detector number is bits 15..8 of ddlId and DDL number is bits 7..0.
  Int_t detNum = ddlId >> 8;
  Int_t ddlNum = ddlId & 0xFF;
  
  if (detNum < 3)
  {
    wordIndex = detNum;
  }
  else if (detNum == 3)
  {
    wordIndex = detNum + (ddlNum >> 5);
  }
  else if (detNum == 4)
  {
    wordIndex = detNum + 7;
  }
  else if (detNum == 5)
  {
    wordIndex = detNum + 7 + (ddlNum >> 5);
  }
  else if (detNum == 30)
  {
    wordIndex = 29;
  }
  else
  {
    wordIndex = detNum + 9;
  }
  
  if (wordIndex < 0 or gkAliHLTDDLListSize <= wordIndex) return false;
  
  // The bit index within the word indicated by wordIndex.
  bitIndex = (ddlId & 0xFF) % 32;
  return true;
}


Bool_t AliHLTReadoutList::GetDDLBit(Int_t ddlId) const
{
  // Fetches the bit value for a particular DDL in the readout list.
  // See header file for more details.
  
  Int_t wordIndex, bitIndex;
  if (! DecodeDDLID(ddlId, wordIndex, bitIndex)) return kFALSE;
  return ((fReadoutList.fList[wordIndex] >> bitIndex) & 0x1) == 0x1;
}


void AliHLTReadoutList::SetDDLBit(Int_t ddlId, Bool_t state)
{
  // Sets the bit value for a particular DDL in the readout list.
  // See header file for more details.
  
  Int_t wordIndex, bitIndex;
  if (! DecodeDDLID(ddlId, wordIndex, bitIndex)) return;

  // To set, 'OR' word with bit mask
  if ( state )
    fReadoutList.fList[wordIndex] |= (0x00000001 << bitIndex);
  // To unset, 'AND' word with bit mask
  else
    fReadoutList.fList[wordIndex] &= (0xFFFFFFFF ^ (0x00000001 << bitIndex));
}


void AliHLTReadoutList::Enable(Int_t detector)
{
  // Enables all DDLs for a particular detector or detectors.
  // See header file for more details.
  
  if ((detector & kITSSPD) != 0) fReadoutList.fList[0] = 0x000FFFFF;
  if ((detector & kITSSDD) != 0) fReadoutList.fList[1] = 0x00FFFFFF;
  if ((detector & kITSSSD) != 0) fReadoutList.fList[2] = 0x0000FFFF;
  if ((detector & kTPC) != 0)
  {
    fReadoutList.fList[3] = 0xFFFFFFFF;
    fReadoutList.fList[4] = 0xFFFFFFFF;
    fReadoutList.fList[5] = 0xFFFFFFFF;
    fReadoutList.fList[6] = 0xFFFFFFFF;
    fReadoutList.fList[7] = 0xFFFFFFFF;
    fReadoutList.fList[8] = 0xFFFFFFFF;
    fReadoutList.fList[9] = 0x00FFFFFF;
    fReadoutList.fList[10] = 0x00000000;
  }
  if ((detector & kTRD) != 0) fReadoutList.fList[11] = 0x0003FFFF;
  if ((detector & kTOF) != 0)
  {
    fReadoutList.fList[12] = 0xFFFFFFFF;
    fReadoutList.fList[13] = 0xFFFFFFFF;
    fReadoutList.fList[14] = 0x000000FF;
  }
  if ((detector & kHMPID) != 0) fReadoutList.fList[15] = 0x00003FFF;
  if ((detector & kPHOS) != 0) fReadoutList.fList[16] = 0x000FFFFF;
  if ((detector & kCPV) != 0) fReadoutList.fList[17] = 0x000003FF;
  if ((detector & kPMD) != 0) fReadoutList.fList[18] = 0x0000003F;
  if ((detector & kMUONTRK) != 0) fReadoutList.fList[19] = 0x000FFFFF;
  if ((detector & kMUONTRG) != 0) fReadoutList.fList[20] = 0x00000003;
  if ((detector & kFMD) != 0) fReadoutList.fList[21] = 0x00000007;
  if ((detector & kT0) != 0) fReadoutList.fList[22] = 0x00000001;
  if ((detector & kV0) != 0) fReadoutList.fList[23] = 0x00000001;
  if ((detector & kZDC) != 0) fReadoutList.fList[24] = 0x00000001;
  if ((detector & kACORDE) != 0) fReadoutList.fList[25] = 0x00000001;
  if ((detector & kTRG) != 0) fReadoutList.fList[26] = 0x00000001;
  if ((detector & kEMCAL) != 0) fReadoutList.fList[27] = 0x00FFFFFF;
  if ((detector & kDAQTEST) != 0) fReadoutList.fList[28] = 0x00000001;
  if ((detector & kHLT) != 0) fReadoutList.fList[29] = 0x000003FF;
}


void AliHLTReadoutList::Disable(Int_t detector)
{
  // Disables all DDLs for a particular detector or detectors.
  // See header file for more details.
  
  if ((detector & kITSSPD) != 0) fReadoutList.fList[0] = 0x00000000;
  if ((detector & kITSSDD) != 0) fReadoutList.fList[1] = 0x00000000;
  if ((detector & kITSSSD) != 0) fReadoutList.fList[2] = 0x00000000;
  if ((detector & kTPC) != 0)
  {
    fReadoutList.fList[3] = 0x00000000;
    fReadoutList.fList[4] = 0x00000000;
    fReadoutList.fList[5] = 0x00000000;
    fReadoutList.fList[6] = 0x00000000;
    fReadoutList.fList[7] = 0x00000000;
    fReadoutList.fList[8] = 0x00000000;
    fReadoutList.fList[9] = 0x00000000;
    fReadoutList.fList[10] = 0x00000000;
  }
  if ((detector & kTRD) != 0) fReadoutList.fList[11] = 0x00000000;
  if ((detector & kTOF) != 0)
  {
    fReadoutList.fList[12] = 0x00000000;
    fReadoutList.fList[13] = 0x00000000;
    fReadoutList.fList[14] = 0x00000000;
  }
  if ((detector & kHMPID) != 0) fReadoutList.fList[15] = 0x00000000;
  if ((detector & kPHOS) != 0) fReadoutList.fList[16] = 0x00000000;
  if ((detector & kCPV) != 0) fReadoutList.fList[17] = 0x00000000;
  if ((detector & kPMD) != 0) fReadoutList.fList[18] = 0x00000000;
  if ((detector & kMUONTRK) != 0) fReadoutList.fList[19] = 0x00000000;
  if ((detector & kMUONTRG) != 0) fReadoutList.fList[20] = 0x00000000;
  if ((detector & kFMD) != 0) fReadoutList.fList[21] = 0x00000000;
  if ((detector & kT0) != 0) fReadoutList.fList[22] = 0x00000000;
  if ((detector & kV0) != 0) fReadoutList.fList[23] = 0x00000000;
  if ((detector & kZDC) != 0) fReadoutList.fList[24] = 0x00000000;
  if ((detector & kACORDE) != 0) fReadoutList.fList[25] = 0x00000000;
  if ((detector & kTRG) != 0) fReadoutList.fList[26] = 0x00000000;
  if ((detector & kEMCAL) != 0) fReadoutList.fList[27] = 0x00000000;
  if ((detector & kDAQTEST) != 0) fReadoutList.fList[28] = 0x00000000;
  if ((detector & kHLT) != 0) fReadoutList.fList[29] = 0x00000000;
}


bool AliHLTReadoutList::DetectorEnabled(Int_t detector) const
{
  // Checks if a particular detector's DDLs are enabled.
  // See header file for more details.
  
  bool result = true;
  if ((detector & kITSSPD) != 0) result &= fReadoutList.fList[0] == 0x000FFFFF;
  if ((detector & kITSSDD) != 0) result &= fReadoutList.fList[1] == 0x00FFFFFF;
  if ((detector & kITSSSD) != 0) result &= fReadoutList.fList[2] == 0x0000FFFF;
  if ((detector & kTPC) != 0)
  {
    result &= fReadoutList.fList[3] == 0xFFFFFFFF;
    result &= fReadoutList.fList[4] == 0xFFFFFFFF;
    result &= fReadoutList.fList[5] == 0xFFFFFFFF;
    result &= fReadoutList.fList[6] == 0xFFFFFFFF;
    result &= fReadoutList.fList[7] == 0xFFFFFFFF;
    result &= fReadoutList.fList[8] == 0xFFFFFFFF;
    result &= fReadoutList.fList[9] == 0x00FFFFFF;
  }
  if ((detector & kTRD) != 0) result &= fReadoutList.fList[11] == 0x0003FFFF;
  if ((detector & kTOF) != 0)
  {
    result &= fReadoutList.fList[12] == 0xFFFFFFFF;
    result &= fReadoutList.fList[13] == 0xFFFFFFFF;
    result &= fReadoutList.fList[14] == 0x000000FF;
  }
  if ((detector & kHMPID) != 0) result &= fReadoutList.fList[15] == 0x00003FFF;
  if ((detector & kPHOS) != 0) result &= fReadoutList.fList[16] == 0x000FFFFF;
  if ((detector & kCPV) != 0) result &= fReadoutList.fList[17] == 0x000003FF;
  if ((detector & kPMD) != 0) result &= fReadoutList.fList[18] == 0x0000003F;
  if ((detector & kMUONTRK) != 0) result &= fReadoutList.fList[19] == 0x000FFFFF;
  if ((detector & kMUONTRG) != 0) result &= fReadoutList.fList[20] == 0x00000003;
  if ((detector & kFMD) != 0) result &= fReadoutList.fList[21] == 0x00000007;
  if ((detector & kT0) != 0) result &= fReadoutList.fList[22] == 0x00000001;
  if ((detector & kV0) != 0) result &= fReadoutList.fList[23] == 0x00000001;
  if ((detector & kZDC) != 0) result &= fReadoutList.fList[24] == 0x00000001;
  if ((detector & kACORDE) != 0) result &= fReadoutList.fList[25] == 0x00000001;
  if ((detector & kTRG) != 0) result &= fReadoutList.fList[26] == 0x00000001;
  if ((detector & kEMCAL) != 0) result &= fReadoutList.fList[27] == 0x00FFFFFF;
  if ((detector & kDAQTEST) != 0) result &= fReadoutList.fList[28] == 0x00000001;
  if ((detector & kHLT) != 0) result &= fReadoutList.fList[29] == 0x000003FF;
  
  return result;
}


void AliHLTReadoutList::Print(Option_t* /*option*/) const
{
  // Prints the DDLs that will be readout according to this readout list.
  
  cout << "Readout enabled for DDLs:" << endl;
  for (Int_t i = 0; i < AliDAQ::kNDetectors; i++)
  {
    Int_t maxddls = AliDAQ::NumberOfDdls(i);
    cout << AliDAQ::DetectorName(i) << ":";
    bool nonefound = true;
    for (Int_t j = 0; j < maxddls; j++)
    {
      Int_t ddlId = ( ((i == AliDAQ::kNDetectors-1) ? 30 : i) << 8 ) + j;
      if (GetDDLBit(ddlId))
      {
        cout << " " << ddlId;
        nonefound = false;
      }
    }
    if (nonefound) cout << " none";
    cout << endl;
  }
}


AliHLTReadoutList& AliHLTReadoutList::operator = (const AliHLTReadoutList& list)
{
  // Assignment operator performs a deep copy.
  
  TObject::operator = (list);
  if (&list != this)
  {
    memcpy(&fReadoutList, &list.fReadoutList, sizeof(fReadoutList));
  }
  return *this;
}


AliHLTReadoutList& AliHLTReadoutList::operator |= (const AliHLTReadoutList& list)
{
  // This operator performs a bitwise inclusive or operation on all DDL bits.
  // See header file for more details.
  
  assert( fReadoutList.fCount == gkAliHLTDDLListSize );
  for (Int_t i = 0; i < gkAliHLTDDLListSize; i++)
  {
    fReadoutList.fList[i] |= list.fReadoutList.fList[i];
  }
  return *this;
}


AliHLTReadoutList& AliHLTReadoutList::operator ^= (const AliHLTReadoutList& list)
{
  // This operator performs a bitwise exclusive or (xor) operation on all DDL bits.
  // See header file for more details.
  
  assert( fReadoutList.fCount == gkAliHLTDDLListSize );
  for (Int_t i = 0; i < gkAliHLTDDLListSize; i++)
  {
    fReadoutList.fList[i] ^= list.fReadoutList.fList[i];
  }
  return *this;
}


AliHLTReadoutList& AliHLTReadoutList::operator &= (const AliHLTReadoutList& list)
{
  // This operator performs a bitwise and operation on all DDL bits.
  // See header file for more details.
  
  assert( fReadoutList.fCount == gkAliHLTDDLListSize );
  for (Int_t i = 0; i < gkAliHLTDDLListSize; i++)
  {
    fReadoutList.fList[i] &= list.fReadoutList.fList[i];
  }
  return *this;
}


AliHLTReadoutList& AliHLTReadoutList::operator -= (const AliHLTReadoutList& list)
{
  // This operator removes all the DDLs specified in list from this readout list.
  // See header file for more details.
  
  assert( fReadoutList.fCount == gkAliHLTDDLListSize );
  for (Int_t i = 0; i < gkAliHLTDDLListSize; i++)
  {
    // Effectively apply: this = this & (~ (this & list))
    // i.e. this = this & (this ^ list)
    fReadoutList.fList[i] &= fReadoutList.fList[i] ^ list.fReadoutList.fList[i];
  }
  return *this;
}


AliHLTReadoutList AliHLTReadoutList::operator ~ () const
{
  // This operator performs a bitwise ones compliment on all DDL bits.
  // See header file for more details.
  
  AliHLTReadoutList readoutlist;
  readoutlist.fReadoutList.fCount = fReadoutList.fCount;
  readoutlist.fReadoutList.fList[0] = 0x000FFFFF & (~fReadoutList.fList[0]);
  readoutlist.fReadoutList.fList[1] = 0x00FFFFFF & (~fReadoutList.fList[1]);
  readoutlist.fReadoutList.fList[2] = 0x0000FFFF & (~fReadoutList.fList[2]);
  readoutlist.fReadoutList.fList[3] = 0xFFFFFFFF & (~fReadoutList.fList[3]);
  readoutlist.fReadoutList.fList[4] = 0xFFFFFFFF & (~fReadoutList.fList[4]);
  readoutlist.fReadoutList.fList[5] = 0xFFFFFFFF & (~fReadoutList.fList[5]);
  readoutlist.fReadoutList.fList[6] = 0xFFFFFFFF & (~fReadoutList.fList[6]);
  readoutlist.fReadoutList.fList[7] = 0xFFFFFFFF & (~fReadoutList.fList[7]);
  readoutlist.fReadoutList.fList[8] = 0xFFFFFFFF & (~fReadoutList.fList[8]);
  readoutlist.fReadoutList.fList[9] = 0x00FFFFFF & (~fReadoutList.fList[9]);
  readoutlist.fReadoutList.fList[10] = 0x00000000 & (~fReadoutList.fList[10]);
  readoutlist.fReadoutList.fList[11] = 0x0003FFFF & (~fReadoutList.fList[11]);
  readoutlist.fReadoutList.fList[12] = 0xFFFFFFFF & (~fReadoutList.fList[12]);
  readoutlist.fReadoutList.fList[13] = 0xFFFFFFFF & (~fReadoutList.fList[13]);
  readoutlist.fReadoutList.fList[14] = 0x000000FF & (~fReadoutList.fList[14]);
  readoutlist.fReadoutList.fList[15] = 0x00003FFF & (~fReadoutList.fList[15]);
  readoutlist.fReadoutList.fList[16] = 0x000FFFFF & (~fReadoutList.fList[16]);
  readoutlist.fReadoutList.fList[17] = 0x000003FF & (~fReadoutList.fList[17]);
  readoutlist.fReadoutList.fList[18] = 0x0000003F & (~fReadoutList.fList[18]);
  readoutlist.fReadoutList.fList[19] = 0x000FFFFF & (~fReadoutList.fList[19]);
  readoutlist.fReadoutList.fList[20] = 0x00000003 & (~fReadoutList.fList[20]);
  readoutlist.fReadoutList.fList[21] = 0x00000007 & (~fReadoutList.fList[21]);
  readoutlist.fReadoutList.fList[22] = 0x00000001 & (~fReadoutList.fList[22]);
  readoutlist.fReadoutList.fList[23] = 0x00000001 & (~fReadoutList.fList[23]);
  readoutlist.fReadoutList.fList[24] = 0x00000001 & (~fReadoutList.fList[24]);
  readoutlist.fReadoutList.fList[25] = 0x00000001 & (~fReadoutList.fList[25]);
  readoutlist.fReadoutList.fList[26] = 0x00000001 & (~fReadoutList.fList[26]);
  readoutlist.fReadoutList.fList[27] = 0x00FFFFFF & (~fReadoutList.fList[27]);
  readoutlist.fReadoutList.fList[28] = 0x00000001 & (~fReadoutList.fList[28]);
  readoutlist.fReadoutList.fList[29] = 0x000003FF & (~fReadoutList.fList[29]);
  return readoutlist;
}

