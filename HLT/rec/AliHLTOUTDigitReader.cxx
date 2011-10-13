// $Id$

//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
//*                  for The ALICE HLT Project.                            *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************

/** @file   AliHLTOUTDigitReader.cxx
    @author Matthias Richter
    @date   
    @brief  HLTOUT data wrapper for simulated AliRoot HLT digit data.
*/

#include "AliHLTOUTDigitReader.h"
#include "AliRawDataHeader.h"
#include "AliDAQ.h"
#include "TTree.h"
#include "TFile.h"
#include "TArrayC.h"
#include "TBranch.h"
#include "TString.h"
#include <cassert>

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTOUTDigitReader)

AliHLTOUTDigitReader::AliHLTOUTDigitReader(int event, AliHLTEsdManager* pEsdManager, const char* digitFile)
  :
  AliHLTOUTHomerCollection(event, pEsdManager),
  fDigitFileName(digitFile),
  fpDigitFile(NULL),
  fpDigitTree(NULL),
  fMinDDL(-1),
  fMaxDDL(-1),
  fppDigitArrays(NULL),
  fpEquipments(NULL),
  fNofDDLs(0),
  fCurrent(-1)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTOUTDigitReader::~AliHLTOUTDigitReader()
{
  // see header file for class documentation
  CloseTree();

  if (fpDigitFile) {
    fpDigitFile->Close();
    delete fpDigitFile;
    fpDigitFile=NULL;
  }
}

Bool_t AliHLTOUTDigitReader::ReadNextData(UChar_t*& data)
{
  // see header file for class documentation
  if (!fppDigitArrays && (!ReadArrays() || !fppDigitArrays))
    return kFALSE;

  if (fCurrent>=fNofDDLs)
    return kFALSE;

  while (++fCurrent<fNofDDLs) {
    if (fMinDDL>-1 && fMinDDL>fpEquipments[fCurrent]) continue;
    if (fMaxDDL>-1 && fMaxDDL<fpEquipments[fCurrent]) continue;
    if (fppDigitArrays[fCurrent]->GetSize()>(int)sizeof(AliRawDataHeader)) {
      data=reinterpret_cast<UChar_t*>(fppDigitArrays[fCurrent]->GetArray());
      data+=sizeof(AliRawDataHeader);
      return kTRUE;
    }
  }
  return kFALSE;
}

int AliHLTOUTDigitReader::Reset()
{
  // see header file for class documentation
  fCurrent=-1;
  return 0;
}

int AliHLTOUTDigitReader::GetDataSize()
{
  // see header file for class documentation
  if (fCurrent>=0 && fCurrent<fNofDDLs && fppDigitArrays) {
    return fppDigitArrays[fCurrent]->GetSize()-sizeof(AliRawDataHeader);
  }
  return 0;
}

const AliRawDataHeader* AliHLTOUTDigitReader::GetDataHeader()
{
  // see header file for class documentation
  if (fCurrent>=0 && fCurrent<fNofDDLs && fppDigitArrays) {
    return reinterpret_cast<AliRawDataHeader*>(fppDigitArrays[fCurrent]->GetArray());
  }
  return NULL;
}

void AliHLTOUTDigitReader::SelectEquipment(int /*equipmentType*/, int minEquipmentId, int maxEquipmentId)
{
  // see header file for class documentation
  fMinDDL=minEquipmentId;
  fMaxDDL=maxEquipmentId;
}

int AliHLTOUTDigitReader::GetEquipmentId()
{
  // see header file for class documentation
  if (fCurrent>=0 && fCurrent<fNofDDLs && fpEquipments) {
    return fpEquipments[fCurrent];
  }
  return -1;
}

bool AliHLTOUTDigitReader::ReadArrays()
{
  // see header file for class documentation
  bool result=false;

  if (GetCurrentEventNo()<0) {
    HLTWarning("no event selected, no data available");
    return false;
  }

  // 2011-06-06 in order to support AliHLTReconstructor option 'ignore-hltout' for
  // digits, the file name can be set to NULL, nothing done in that case 
  if (!fDigitFileName) return false;

  if (!fpDigitFile) {
    fpDigitFile=new TFile(fDigitFileName);
  }
  if (!fpDigitFile) return false;
  if (fpDigitFile->IsZombie()) return false;

  fpDigitFile->GetObject("rawhltout", fpDigitTree);
  if (!fpDigitTree) {
    return false;
  }

  vector<int> equipments;
  int iTotal=0;
  TIter iter(fpDigitTree->GetListOfBranches());
  while (TBranch* br=(TBranch *)iter.Next()) {
    iTotal++;
    TString bname=br->GetName();
    if (bname.BeginsWith("HLT_") && bname.EndsWith(".ddl")) {
      bname.ReplaceAll("HLT_", "");
      bname.ReplaceAll(".ddl", "");
      if (bname.IsDigit()) {
	int equipment=bname.Atoi();
	int index=equipment-AliDAQ::DdlIDOffset("HLT");
	if (index>=0 && strcmp(br->GetName(), AliDAQ::DdlFileName("HLT", index))==0) {
	  equipments.push_back(equipment);
	} else {
	  HLTWarning("equipment mismatch: skipping id %d (out of range of HLT equipments)", equipment);
	}
      }
    }
  }
  HLTDebug("found %d HLT branche(s) out of %d", equipments.size(), iTotal);

  if (equipments.size()>0) {
    fNofDDLs=equipments.size();
    fpEquipments=new int[fNofDDLs];
    fppDigitArrays=new TArrayC*[fNofDDLs];
    if (fpEquipments && fppDigitArrays) {
      memcpy(fpEquipments, &equipments[0], fNofDDLs*sizeof(int));
      int i=0;
      for (i=0; i<fNofDDLs; i++) {
	fppDigitArrays[i]=NULL;
	fpDigitTree->SetBranchAddress(AliDAQ::DdlFileName("HLT", fpEquipments[i]-AliDAQ::DdlIDOffset("HLT")), &fppDigitArrays[i]);
      }
      fpDigitTree->GetEvent(GetCurrentEventNo());
      for (i=0; i<fNofDDLs; i++) {
	if (fppDigitArrays[i]) {
	  HLTDebug("branch %s: %d byte(s)", AliDAQ::DdlFileName("HLT", fpEquipments[i]-AliDAQ::DdlIDOffset("HLT")), fppDigitArrays[i]);
	}
      }
      result=true;
    } else {
      fNofDDLs=0;
    }
  }
  return result;
}

int AliHLTOUTDigitReader::CloseTree()
{
  // see header file for class documentation
  int iResult=0;
  for (int i=0; i<fNofDDLs; i++) {
    if (fppDigitArrays[i]) delete fppDigitArrays[i];
    fppDigitArrays[i]=NULL;
  }
  if (fppDigitArrays) delete[] fppDigitArrays;
  fppDigitArrays=NULL;
  if (fpEquipments) delete[] fpEquipments;
  fpEquipments=NULL;
  if (fpDigitTree) delete fpDigitTree;
  fpDigitTree=NULL;
  fNofDDLs=0;
  fCurrent=-1;

  return iResult;
}

void AliHLTOUTDigitReader::SetParam(TTree* /*pDigitTree*/, int event)
{
  // see header file for class documentation

  // TODO: here we have the implemented to correct loading of
  // the digit file from the run loader. At time of writing
  // (Jun 2008) the AliLoader for the HLT is missing in the AliRoot
  // framework.
  fEvent=event;
}

void AliHLTOUTDigitReader::SetParam(const char* filename, int event)
{
  // see header file for class documentation

  if (filename && filename[0]!=0) {
    fDigitFileName=filename;
  } else {
    HLTWarning("no valid digit file provided, using default file %s", fDigitFileName.Data());
  }
  fEvent=event;
}
