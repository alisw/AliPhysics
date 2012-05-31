// $Id$
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
#include "AliHLTDAQ.h"
#include "Riostream.h"
#include "TString.h"
#include "TObjString.h"
#include "TObjArray.h"
#include <cassert>

ClassImp(AliHLTReadoutList)


const char* AliHLTReadoutList::DetectorIdToString(EDetectorId id)
{
  // Converts a detector ID to a user readable string.
  switch (id)
  {
  case kNoDetector: return "kNoDetector";
  case kITSSPD:  return "kITSSPD";
  case kITSSDD:  return "kITSSDD";
  case kITSSSD:  return "kITSSSD";
  case kTPC:     return "kTPC";
  case kTRD:     return "kTRD";
  case kTOF:     return "kTOF";
  case kHMPID:   return "kHMPID";
  case kPHOS:    return "kPHOS";
  case kCPV:     return "kCPV";
  case kPMD:     return "kPMD";
  case kMUONTRK: return "kMUONTRK";
  case kMUONTRG: return "kMUONTRG";
  case kFMD:     return "kFMD";
  case kT0:      return "kT0";
  case kV0:      return "kV0";
  case kZDC:     return "kZDC";
  case kACORDE:  return "kACORDE";
  case kTRG:     return "kTRG";
  case kEMCAL:   return "kEMCAL";
  case kDAQTEST: return "kDAQTEST";
  case kHLT:     return "kHLT";
  case kALLDET:  return "kALLDET";
  default:       return "UNKNOWN!";
  }
}


AliHLTReadoutList::AliHLTReadoutList() :
        TNamed("AliHLTReadoutList", "Readout list object used for manipulating and storing an AliHLTEventDDL structure."),
	fReadoutList()
{
  // Default constructor.
  
  fReadoutList.fCount = gkAliHLTDDLListSize;  // Required by ALICE-INT-2007-015
  memset(fReadoutList.fList, 0x0, sizeof(fReadoutList.fList));
}


AliHLTReadoutList::AliHLTReadoutList(Int_t enabledDetectors) :
        TNamed("AliHLTReadoutList", "Readout list object used for manipulating and storing an AliHLTEventDDL structure."),
	fReadoutList()
{
  // Constructor to select which detectors to enable for readout.
  // See header file for more details.
  
  fReadoutList.fCount = gkAliHLTDDLListSize;  // Required by ALICE-INT-2007-015
  memset(fReadoutList.fList, 0x0, sizeof(fReadoutList.fList));
  Enable(enabledDetectors);
}


AliHLTReadoutList::AliHLTReadoutList(const char* enabledList) :
        TNamed("AliHLTReadoutList", "Readout list object used for manipulating and storing an AliHLTEventDDL structure."),
	fReadoutList()
{
  // Constructor to select which detectors and DDLs to enable for readout.
  // See header file for more details.
  
  fReadoutList.fCount = gkAliHLTDDLListSize;  // Required by ALICE-INT-2007-015
  memset(fReadoutList.fList, 0x0, sizeof(fReadoutList.fList));
  
  TString str(enabledList);
  str.ToUpper();
  Int_t enabledDetectors = 0;
  TObjArray* list = str.Tokenize(" ");
  TIter next(list);
  const TObjString* objstr = NULL;
  while ((objstr = dynamic_cast<const TObjString*>(next())) != NULL)
  {
    str = objstr->GetString();
    if (str.IsDigit()) EnableDDLBit(str.Atoi());
    if (str == "ITSSPD") enabledDetectors |= kITSSPD;
    if (str == "ITSSDD") enabledDetectors |= kITSSDD;
    if (str == "ITSSSD") enabledDetectors |= kITSSSD;
    if (str == "TPC") enabledDetectors |= kTPC;
    if (str == "TRD") enabledDetectors |= kTRD;
    if (str == "TOF") enabledDetectors |= kTOF;
    if (str == "HMPID") enabledDetectors |= kHMPID;
    if (str == "PHOS") enabledDetectors |= kPHOS;
    if (str == "CPV") enabledDetectors |= kCPV;
    if (str == "PMD") enabledDetectors |= kPMD;
    if (str == "MUONTRK") enabledDetectors |= kMUONTRK;
    if (str == "MUONTRG") enabledDetectors |= kMUONTRG;
    if (str == "FMD") enabledDetectors |= kFMD;
    if (str == "T0") enabledDetectors |= kT0;
    if (str == "V0") enabledDetectors |= kV0;
    if (str == "ZDC") enabledDetectors |= kZDC;
    if (str == "ACORDE") enabledDetectors |= kACORDE;
    if (str == "TRG") enabledDetectors |= kTRG;
    if (str == "EMCAL") enabledDetectors |= kEMCAL;
    if (str == "DAQTEST") enabledDetectors |= kDAQTEST;
    if (str == "HLT") enabledDetectors |= kHLT;
    if (str == "ALL") enabledDetectors |= kALLDET;
  }
  delete list;
  Enable(enabledDetectors);
}


AliHLTReadoutList::AliHLTReadoutList(const AliHLTEventDDL& list) :
        TNamed("AliHLTReadoutList", "Readout list object used for manipulating and storing an AliHLTEventDDL structure."),
	fReadoutList()
{
  // Constructor to create readout list from AliHLTEventDDL structure.
  // See header file for more details.
  FillStruct(list);
}


AliHLTReadoutList::AliHLTReadoutList(const AliHLTReadoutList& list) :
	TNamed(list),
	fReadoutList()
{
  // Copy constructor performs a deep copy.
  
  if (list.fReadoutList.fCount == (unsigned)gkAliHLTDDLListSize)
  {
    memcpy(&fReadoutList, &list.fReadoutList, sizeof(fReadoutList));
  }
  else
  {
    FillStruct(list);
  }
}


void AliHLTReadoutList::FillStruct(const AliHLTEventDDL& list)
{
  // Fills internal DDL bits structure.

  fReadoutList.fCount = gkAliHLTDDLListSize;  // Required by ALICE-INT-2007-015
  memset(fReadoutList.fList, 0x0, sizeof(fReadoutList.fList));
  // Handle lists of different sizes. If the size is for a known version
  // of AliHLTEventDDL then handle appropriately, otherwise just copy only
  // the overlapping part of the list.
  if (list.fCount == (unsigned)gkAliHLTDDLListSizeV0)
  {
    memcpy(&fReadoutList.fList[0], &list.fList[0], sizeof(AliHLTUInt32_t)*28);
    memcpy(&fReadoutList.fList[29], &list.fList[28], sizeof(AliHLTUInt32_t)*2);
  }
  else if (list.fCount == (unsigned)gkAliHLTDDLListSizeV1)
  {
    memcpy(&fReadoutList, &list, sizeof(AliHLTEventDDL));
  }
  else
  {
    memcpy(&fReadoutList.fList, &list.fList, (fReadoutList.fCount<list.fCount?fReadoutList.fCount:list.fCount)*sizeof(AliHLTUInt32_t));
  }
}


AliHLTReadoutList::~AliHLTReadoutList()
{
  // Default destructor.
}


bool AliHLTReadoutList::Empty() const
{
  // Returns true if the readout list has no DDLs enabled.

  for (size_t i = 0; i < sizeof(fReadoutList.fList) / sizeof(fReadoutList.fList[0]); i++)
  {
    if (fReadoutList.fList[i] != 0x0) return false;
  }
  return true;
}


void AliHLTReadoutList::Clear(Option_t* /*option*/)
{
  // Resets all the DDL readout bits.
  memset(fReadoutList.fList, 0x0, sizeof(fReadoutList.fList));

#if 1 // ROOT_SVN_REVISION < 9999  //FIXME: after fixed https://savannah.cern.ch/bugs/?69241
  // Check if we need to convert to new format and do so.
  if (fReadoutList.fCount == (unsigned)gkAliHLTDDLListSizeV0)
  {
    fReadoutList.fCount = gkAliHLTDDLListSize;
  }
#endif
}


bool AliHLTReadoutList::DecodeDDLID(Int_t ddlId, Int_t& wordIndex, Int_t& bitIndex)
{
  // Decodes the word index and bit index within that word for the readout list structure.
  // See header file for more details.
  
  // The detector number is bits 15..8 of ddlId and DDL number is bits 7..0.
  Int_t detNum = ddlId >> 8;
  Int_t ddlNum = ddlId & 0xFF;
  
  switch (detNum)
  {
  case 0: // SPD
  case 1: // SDD
  case 2: // SSD
    if (ddlNum >= 32) return false; // only have 1 32-bit word.
    // the 3 ITS detectors have one word each
    wordIndex = detNum;
    break;
  case 3: // TPC
    // the TPC bitfield has in total 8 words
    wordIndex = detNum + (ddlNum >> 5);
    break;
  case 4: // TRD
    if (ddlNum >= 32) return false; // only have 1 32-bit word.
    // the TRD bitfield starts at word 11 (3 words ITS + 8 words TPC)
    wordIndex = 11;
    break;
  case 5: // TOF
    if (ddlNum >= 3*32) return false; // only have 3 32-bit words.
    // TOF has 72 DDLs, the bitfield is 3 words starting at position 12
    wordIndex = 12 + (ddlNum >> 5);
    break;
  case 6: // HMPID
  case 7: // PHOS
  case 8: // CPV
  case 9: // PMD
  case 10: // MUONTRK (MCH)
  case 11: // MUONTRG (MTR)
  case 12: // FMD
  case 13: // T0
  case 14: // V0
  case 15: // ZDC
  case 16: // ACORDE
  case 17: // TRG
    if (ddlNum >= 32) return false; // only have 1 32-bit word.
    // all these detectors fit into one word, the offset is due to
    // TPC and TOF
    wordIndex = detNum + 9;
    break;
  case 18: // EMCAL
    if (ddlNum >= 2*32) return false; // only have 2 32-bit words.
    // 2 words for EMCAL + DCAL
    wordIndex = detNum + 7;
    wordIndex = 27 + (ddlNum >> 5);
    break;
  case 19: // DAQTEST
    if (ddlNum >= 32) return false; // only have 1 32-bit word.
    wordIndex = 29;
    break;
  case 30: // HLT
    if (ddlNum >= 32) return false; // only have 1 32-bit word.
    // the HLT bitfield is in the last word
    wordIndex = 30;
    break;
  default:
    return false;
  }
  
  if (ddlNum >= AliHLTDAQ::NumberOfDdls(detNum == 30 ? 20 : detNum)) return false;
  
  // The bit index within the word indicated by wordIndex.
  bitIndex = ddlNum % 32;
  return true;
}


Bool_t AliHLTReadoutList::GetDDLBit(Int_t ddlId) const
{
  // Fetches the bit value for a particular DDL in the readout list.
  // See header file for more details.
  
  Int_t wordIndex, bitIndex;
  if (! DecodeDDLID(ddlId, wordIndex, bitIndex)) return kFALSE;

#if 1 // ROOT_SVN_REVISION < 9999  //FIXME: after fixed https://savannah.cern.ch/bugs/?69241
  // Check if we need to convert to new format and do so.
  if (fReadoutList.fCount == (unsigned)gkAliHLTDDLListSizeV0)
  {
    if (wordIndex == 27)
    {
      if (bitIndex >= 24) return kFALSE;
    }
    else if (wordIndex == 28)
    {
      return kFALSE;
    }
    else if (wordIndex > 28)
    {
      --wordIndex;
    }
  }
#endif

  return ((fReadoutList.fList[wordIndex] >> bitIndex) & 0x1) == 0x1;
}


void AliHLTReadoutList::SetDDLBit(Int_t ddlId, Bool_t state)
{
  // Sets the bit value for a particular DDL in the readout list.
  // See header file for more details.

#if 1 // ROOT_SVN_REVISION < 9999  //FIXME: after fixed https://savannah.cern.ch/bugs/?69241
  // Check if we need to convert to new format and do so.
  if (fReadoutList.fCount == (unsigned)gkAliHLTDDLListSizeV0)
  {
    AliHLTEventDDL copy = fReadoutList;
    FillStruct(copy);
  }
#endif
  assert(fReadoutList.fCount == gkAliHLTDDLListSize);
  
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

#if 1 // ROOT_SVN_REVISION < 9999  //FIXME: after fixed https://savannah.cern.ch/bugs/?69241
  // Check if we need to convert to new format and do so.
  if (fReadoutList.fCount == (unsigned)gkAliHLTDDLListSizeV0)
  {
    AliHLTEventDDL copy = fReadoutList;
    FillStruct(copy);
  }
#endif
  assert(fReadoutList.fCount == gkAliHLTDDLListSize);
  
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
  if ((detector & kHMPID) != 0) fReadoutList.fList[15] = 0x000FFFFF;
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
  if ((detector & kEMCAL) != 0)
  {
    fReadoutList.fList[27] = 0xFFFFFFFF;
    fReadoutList.fList[28] = 0x00003FFF;
  }
  if ((detector & kDAQTEST) != 0) fReadoutList.fList[29] = 0x00000001;
  if ((detector & kHLT) != 0) fReadoutList.fList[30] = 0x0FFFFFFF;
}


void AliHLTReadoutList::Disable(Int_t detector)
{
  // Disables all DDLs for a particular detector or detectors.
  // See header file for more details.

#if 1 // ROOT_SVN_REVISION < 9999  //FIXME: after fixed https://savannah.cern.ch/bugs/?69241
  // Check if we need to convert to new format and do so.
  if (fReadoutList.fCount == (unsigned)gkAliHLTDDLListSizeV0)
  {
    AliHLTEventDDL copy = fReadoutList;
    FillStruct(copy);
  }
#endif
  assert(fReadoutList.fCount == gkAliHLTDDLListSize);
  
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
  if ((detector & kEMCAL) != 0)
  {
    fReadoutList.fList[27] = 0x00000000;
    fReadoutList.fList[28] = 0x00000000;
  }
  if ((detector & kDAQTEST) != 0) fReadoutList.fList[29] = 0x00000000;
  if ((detector & kHLT) != 0) fReadoutList.fList[30] = 0x00000000;
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
  if ((detector & kHMPID) != 0) result &= fReadoutList.fList[15] == 0x000FFFFF;
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
#if 1 // ROOT_SVN_REVISION < 9999  //FIXME: after fixed https://savannah.cern.ch/bugs/?69241
  if (fReadoutList.fCount == (unsigned)gkAliHLTDDLListSizeV0)
  {
    if ((detector & kEMCAL) != 0) result &= fReadoutList.fList[27] == 0x00FFFFFF;
    if ((detector & kDAQTEST) != 0) result &= fReadoutList.fList[28] == 0x00000001;
    if ((detector & kHLT) != 0) result &= fReadoutList.fList[29] == 0x000003FF;
  }
  else
#endif
  {
    if ((detector & kEMCAL) != 0)
    {
      result &= fReadoutList.fList[27] == 0xFFFFFFFF;
      result &= fReadoutList.fList[28] == 0x00003FFF;
    }
    if ((detector & kDAQTEST) != 0) result &= fReadoutList.fList[29] == 0x00000001;
    if ((detector & kHLT) != 0) result &= fReadoutList.fList[30] == 0x0FFFFFFF;
  }
  
  return result;
}


bool AliHLTReadoutList::DetectorDisabled(Int_t detector) const
{
  // Checks if a particular detector's DDLs are disabled.
  // See header file for more details.
  
  bool result = true;
  if ((detector & kITSSPD) != 0) result &= fReadoutList.fList[0] == 0x00000000;
  if ((detector & kITSSDD) != 0) result &= fReadoutList.fList[1] == 0x00000000;
  if ((detector & kITSSSD) != 0) result &= fReadoutList.fList[2] == 0x00000000;
  if ((detector & kTPC) != 0)
  {
    result &= fReadoutList.fList[3] == 0x00000000;
    result &= fReadoutList.fList[4] == 0x00000000;
    result &= fReadoutList.fList[5] == 0x00000000;
    result &= fReadoutList.fList[6] == 0x00000000;
    result &= fReadoutList.fList[7] == 0x00000000;
    result &= fReadoutList.fList[8] == 0x00000000;
    result &= fReadoutList.fList[9] == 0x00000000;
  }
  if ((detector & kTRD) != 0) result &= fReadoutList.fList[11] == 0x00000000;
  if ((detector & kTOF) != 0)
  {
    result &= fReadoutList.fList[12] == 0x00000000;
    result &= fReadoutList.fList[13] == 0x00000000;
    result &= fReadoutList.fList[14] == 0x00000000;
  }
  if ((detector & kHMPID) != 0) result &= fReadoutList.fList[15] == 0x00000000;
  if ((detector & kPHOS) != 0) result &= fReadoutList.fList[16] == 0x00000000;
  if ((detector & kCPV) != 0) result &= fReadoutList.fList[17] == 0x00000000;
  if ((detector & kPMD) != 0) result &= fReadoutList.fList[18] == 0x00000000;
  if ((detector & kMUONTRK) != 0) result &= fReadoutList.fList[19] == 0x00000000;
  if ((detector & kMUONTRG) != 0) result &= fReadoutList.fList[20] == 0x00000000;
  if ((detector & kFMD) != 0) result &= fReadoutList.fList[21] == 0x00000000;
  if ((detector & kT0) != 0) result &= fReadoutList.fList[22] == 0x00000000;
  if ((detector & kV0) != 0) result &= fReadoutList.fList[23] == 0x00000000;
  if ((detector & kZDC) != 0) result &= fReadoutList.fList[24] == 0x00000000;
  if ((detector & kACORDE) != 0) result &= fReadoutList.fList[25] == 0x00000000;
  if ((detector & kTRG) != 0) result &= fReadoutList.fList[26] == 0x00000000;
#if 1 // ROOT_SVN_REVISION < 9999  //FIXME: after fixed https://savannah.cern.ch/bugs/?69241
  if (fReadoutList.fCount == (unsigned)gkAliHLTDDLListSizeV0)
  {
    if ((detector & kEMCAL) != 0) result &= fReadoutList.fList[27] == 0x00000000;
    if ((detector & kDAQTEST) != 0) result &= fReadoutList.fList[28] == 0x00000000;
    if ((detector & kHLT) != 0) result &= fReadoutList.fList[29] == 0x00000000;
  }
  else
#endif
  {
    if ((detector & kEMCAL) != 0)
    {
      result &= fReadoutList.fList[27] == 0x00000000;
      result &= fReadoutList.fList[28] == 0x00000000;
    }
    if ((detector & kDAQTEST) != 0) result &= fReadoutList.fList[29] == 0x00000000;
    if ((detector & kHLT) != 0) result &= fReadoutList.fList[30] == 0x00000000;
  }
  
  return result;
}


Int_t AliHLTReadoutList::GetFirstWord(EDetectorId detector)
{
  // See header file for more details.
  switch (detector)
  {
  case kITSSPD:  return 0;
  case kITSSDD:  return 1;
  case kITSSSD:  return 2;
  case kTPC:     return 3;
  case kTRD:     return 11;
  case kTOF:     return 12;
  case kHMPID:   return 15;
  case kPHOS:    return 16;
  case kCPV:     return 17;
  case kPMD:     return 18;
  case kMUONTRK: return 19;
  case kMUONTRG: return 20;
  case kFMD:     return 21;
  case kT0:      return 22;
  case kV0:      return 23;
  case kZDC:     return 24;
  case kACORDE:  return 25;
  case kTRG:     return 26;
  case kEMCAL:   return 27;
  case kDAQTEST: return 29;
  case kHLT:     return 30;
  default:       return -1;
  }
}


Int_t AliHLTReadoutList::GetWordCount(EDetectorId detector)
{
  // See header file for more details.
  switch (detector)
  {
  case kITSSPD:  return 1;
  case kITSSDD:  return 1;
  case kITSSSD:  return 1;
  case kTPC:     return 8;
  case kTRD:     return 1;
  case kTOF:     return 3;
  case kHMPID:   return 1;
  case kPHOS:    return 1;
  case kCPV:     return 1;
  case kPMD:     return 1;
  case kMUONTRK: return 1;
  case kMUONTRG: return 1;
  case kFMD:     return 1;
  case kT0:      return 1;
  case kV0:      return 1;
  case kZDC:     return 1;
  case kACORDE:  return 1;
  case kTRG:     return 1;
  case kEMCAL:   return 2;
  case kDAQTEST: return 1;
  case kHLT:     return 1;
  default:       return 0;
  }
}


AliHLTReadoutList::EDetectorId AliHLTReadoutList::GetDetectorFromWord(Int_t wordindex)
{
  // See header file for more details.
  switch (wordindex)
  {
  case 0: return kITSSPD;
  case 1: return kITSSDD;
  case 2: return kITSSSD;
  case 3: return kTPC;
  case 4: return kTPC;
  case 5: return kTPC;
  case 6: return kTPC;
  case 7: return kTPC;
  case 8: return kTPC;
  case 9: return kTPC;
  case 10: return kTPC;
  case 11: return kTRD;
  case 12: return kTOF;
  case 13: return kTOF;
  case 14: return kTOF;
  case 15: return kHMPID;
  case 16: return kPHOS;
  case 17: return kCPV;
  case 18: return kPMD;
  case 19: return kMUONTRK;
  case 20: return kMUONTRG;
  case 21: return kFMD;
  case 22: return kT0;
  case 23: return kV0;
  case 24: return kZDC;
  case 25: return kACORDE;
  case 26: return kTRG;
  case 27: return kEMCAL;
  case 28: return kEMCAL;
  case 29: return kDAQTEST;
  case 30: return kHLT;
  default: return kNoDetector;
  }
}


AliHLTReadoutList::EDetectorId AliHLTReadoutList::GetFirstUsedDetector(EDetectorId startAfter) const
{
  // See header file for more details.

  if (startAfter < kITSSPD and fReadoutList.fList[0] != 0x00000000) return kITSSPD;
  if (startAfter < kITSSDD and fReadoutList.fList[1] != 0x00000000) return kITSSDD;
  if (startAfter < kITSSSD and fReadoutList.fList[2] != 0x00000000) return kITSSSD;
  if (startAfter < kTPC and fReadoutList.fList[3] != 0x00000000) return kTPC;
  if (startAfter < kTPC and fReadoutList.fList[4] != 0x00000000) return kTPC;
  if (startAfter < kTPC and fReadoutList.fList[5] != 0x00000000) return kTPC;
  if (startAfter < kTPC and fReadoutList.fList[6] != 0x00000000) return kTPC;
  if (startAfter < kTPC and fReadoutList.fList[7] != 0x00000000) return kTPC;
  if (startAfter < kTPC and fReadoutList.fList[8] != 0x00000000) return kTPC;
  if (startAfter < kTPC and fReadoutList.fList[9] != 0x00000000) return kTPC;
  if (startAfter < kTPC and fReadoutList.fList[10] != 0x00000000) return kTPC;
  if (startAfter < kTRD and fReadoutList.fList[11] != 0x00000000) return kTRD;
  if (startAfter < kTOF and fReadoutList.fList[12] != 0x00000000) return kTOF;
  if (startAfter < kTOF and fReadoutList.fList[13] != 0x00000000) return kTOF;
  if (startAfter < kTOF and fReadoutList.fList[14] != 0x00000000) return kTOF;
  if (startAfter < kHMPID and fReadoutList.fList[15] != 0x00000000) return kHMPID;
  if (startAfter < kPHOS and fReadoutList.fList[16] != 0x00000000) return kPHOS;
  if (startAfter < kCPV and fReadoutList.fList[17] != 0x00000000) return kCPV;
  if (startAfter < kPMD and fReadoutList.fList[18] != 0x00000000) return kPMD;
  if (startAfter < kMUONTRK and fReadoutList.fList[19] != 0x00000000) return kMUONTRK;
  if (startAfter < kMUONTRG and fReadoutList.fList[20] != 0x00000000) return kMUONTRG;
  if (startAfter < kFMD and fReadoutList.fList[21] != 0x00000000) return kFMD;
  if (startAfter < kT0 and fReadoutList.fList[22] != 0x00000000) return kT0;
  if (startAfter < kV0 and fReadoutList.fList[23] != 0x00000000) return kV0;
  if (startAfter < kZDC and fReadoutList.fList[24] != 0x00000000) return kZDC;
  if (startAfter < kACORDE and fReadoutList.fList[25] != 0x00000000) return kACORDE;
  if (startAfter < kTRG and fReadoutList.fList[26] != 0x00000000) return kTRG;
  if (startAfter < kEMCAL and fReadoutList.fList[27] != 0x00000000) return kEMCAL;
#if 1 // ROOT_SVN_REVISION < 9999  //FIXME: after fixed https://savannah.cern.ch/bugs/?69241
  // Check if we need to convert to new format and do so.
  if (fReadoutList.fCount == (unsigned)gkAliHLTDDLListSizeV0)
  {
    if (startAfter < kDAQTEST and fReadoutList.fList[28] != 0x00000000) return kDAQTEST;
    if (startAfter < kHLT and fReadoutList.fList[29] != 0x00000000) return kHLT;
  }
  else
#endif
  {
    if (startAfter < kEMCAL and fReadoutList.fList[28] != 0x00000000) return kEMCAL;
    if (startAfter < kDAQTEST and fReadoutList.fList[29] != 0x00000000) return kDAQTEST;
    if (startAfter < kHLT and fReadoutList.fList[30] != 0x00000000) return kHLT;
  }
  return kNoDetector;
}


void AliHLTReadoutList::Print(Option_t* /*option*/) const
{
  // Prints the DDLs that will be readout according to this readout list.
  
  cout << "Readout enabled for DDLs:" << endl;
  for (Int_t i = 0; i < AliHLTDAQ::NumberOfDetectors(); i++)
  {
    Int_t maxddls = AliHLTDAQ::NumberOfDdls(i);
    cout << AliHLTDAQ::DetectorName(i) << ":";
    bool nonefound = true;
    for (Int_t j = 0; j < maxddls; j++)
    {
      Int_t ddlId = ( ((i == AliHLTDAQ::NumberOfDetectors()-1) ? 30 : i) << 8 ) + j;
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
    if (list.fReadoutList.fCount == (unsigned)gkAliHLTDDLListSize)
    {
      memcpy(&fReadoutList, &list.fReadoutList, sizeof(fReadoutList));
    }
    else
    {
      FillStruct(list);
    }
  }
  return *this;
}


AliHLTReadoutList& AliHLTReadoutList::operator |= (const AliHLTReadoutList& list)
{
  // This operator performs a bitwise inclusive or operation on all DDL bits.
  // See header file for more details.
  this->OrEq(list);
  return *this;
}

AliHLTReadoutList& AliHLTReadoutList::OrEq(const AliHLTReadoutList& list)
{
  // a bitwise inclusive or operation on all DDL bits.
  // See header file for more details.
  
  assert( fReadoutList.fCount == (unsigned)gkAliHLTDDLListSize );
  assert( fReadoutList.fCount == list.fReadoutList.fCount );
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

  this->XorEq(list);
  return *this;
}

AliHLTReadoutList& AliHLTReadoutList::XorEq(const AliHLTReadoutList& list)
{
  // bitwise exclusive or (xor) operation on all DDL bits.
  // See header file for more details.
  
  assert( fReadoutList.fCount == (unsigned)gkAliHLTDDLListSize );
  assert( fReadoutList.fCount == list.fReadoutList.fCount );
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

  this->AndEq(list);
  return *this;
}

AliHLTReadoutList& AliHLTReadoutList::AndEq(const AliHLTReadoutList& list)
{
  // bitwise and operation on all DDL bits.
  // See header file for more details.

  assert( fReadoutList.fCount == (unsigned)gkAliHLTDDLListSize );
  assert( fReadoutList.fCount == list.fReadoutList.fCount );
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
  
  assert( fReadoutList.fCount == (unsigned)gkAliHLTDDLListSize );
  assert( fReadoutList.fCount == list.fReadoutList.fCount );
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
  readoutlist.fReadoutList.fCount = gkAliHLTDDLListSize;
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
  readoutlist.fReadoutList.fList[10] = 0x00000000;// & (~fReadoutList.fList[10]); // Commented out the end part to suppress coverty warning.
  readoutlist.fReadoutList.fList[11] = 0x0003FFFF & (~fReadoutList.fList[11]);
  readoutlist.fReadoutList.fList[12] = 0xFFFFFFFF & (~fReadoutList.fList[12]);
  readoutlist.fReadoutList.fList[13] = 0xFFFFFFFF & (~fReadoutList.fList[13]);
  readoutlist.fReadoutList.fList[14] = 0x000000FF & (~fReadoutList.fList[14]);
  readoutlist.fReadoutList.fList[15] = 0x000FFFFF & (~fReadoutList.fList[15]);
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
#if 1 // ROOT_SVN_REVISION < 9999  //FIXME: after fixed https://savannah.cern.ch/bugs/?69241
  // Check if we need to convert to new format and do so.
  if (fReadoutList.fCount == (unsigned)gkAliHLTDDLListSizeV0)
  {
    readoutlist.fReadoutList.fList[27] = 0x00FFFFFF & (~fReadoutList.fList[27]);
    readoutlist.fReadoutList.fList[28] = 0x00000000;
    readoutlist.fReadoutList.fList[29] = 0x00000001 & (~fReadoutList.fList[28]);
    readoutlist.fReadoutList.fList[30] = 0x000003FF & (~fReadoutList.fList[29]);
  }
  else
#endif
  {
    readoutlist.fReadoutList.fList[27] = 0xFFFFFFFF & (~fReadoutList.fList[27]);
    readoutlist.fReadoutList.fList[28] = 0x00003FFF & (~fReadoutList.fList[28]);
    readoutlist.fReadoutList.fList[29] = 0x00000001 & (~fReadoutList.fList[29]);
    readoutlist.fReadoutList.fList[30] = 0x0FFFFFFF & (~fReadoutList.fList[30]);
  }
  return readoutlist;
}

#if ROOT_VERSION_CODE < ROOT_VERSION(5,26,0)
void AliHLTReadoutList::Streamer(TBuffer &R__b)
{
   // Stream an object of class AliHLTReadoutList.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(AliHLTReadoutList::Class(),this);
      // Convert old structure to new version if necessary.
      if (fReadoutList.fCount == (unsigned)gkAliHLTDDLListSizeV0)
      {
        fReadoutList.fList[30] = fReadoutList.fList[29];
        fReadoutList.fList[29] = fReadoutList.fList[28];
        fReadoutList.fList[28] = 0x0;
        fReadoutList.fCount = gkAliHLTDDLListSizeV1;
      }
   } else {
      R__b.WriteClassBuffer(AliHLTReadoutList::Class(),this);
   }
}
#endif // ROOT version check.
