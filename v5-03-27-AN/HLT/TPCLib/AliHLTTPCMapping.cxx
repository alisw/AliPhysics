// $Id$
/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Kenneth Aamodt <kenneth@ift.uib.no>                   *
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

/** @file   AliHLTTPCMapping.cxx
    @author Kenneth Aamodt
    @date   
    @brief  A mapping class for the TPC.
*/

// see header file for class documentation                                   //
// or                                                                        //
// refer to README to build package                                          //
// or                                                                        //
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt                          //

#if __GNUC__>= 3
using namespace std;
#endif
#include <cassert>
#include "AliHLTTPCMapping.h"
#include "AliHLTTPCTransform.h"

ClassImp(AliHLTTPCMapping)

AliHLTTPCMapping::AliHLTTPCMapping(UInt_t patch)
  :
  fPatch(patch),
  fCurrentRowMapping(NULL),
  fCurrentPadMapping(NULL),
  fNofRows(0),
  fMaxHWAdd(0),
  fRowOffset(0)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

  assert(patch<fgkNofPatches);
  if (patch>=fgkNofPatches) {
    HLTFatal("invalid partition number %d, allowed range [0,%d]", patch, fgkNofPatches-1);
    fPatch=0; // just to avoid boundary overflow
    return;
  }

  // Get the row offset.
  // we have several possibilities to count the rows:
  // - A: absolute number: 0 to 158 over all partitions
  // - B: sectorwise: 0 to 62 for inner, 0 to 95 for outer sector
  // - C: within a partition: the mappping class is designed to return the
  //   mapping within a partition.
  // The mapping files use scheme B. We have to subtract the first row.
  // AliHLTTPCTransform::GetFirstRow returns first row in scheme A.
  fNofRows=AliHLTTPCTransform::GetNRows(patch);
  fRowOffset=AliHLTTPCTransform::GetFirstRow(patch);

  if(!fgMappingIsDone[patch]){
    ReadArray(patch, fgkMappingSize[patch], fgRowMapping[patch], fgPadMapping[patch], fgkMappingHwaSize, fgHwaMapping[patch]);
    fgMappingIsDone[patch]=kTRUE;
  }
  fCurrentRowMapping=fgRowMapping[patch];
  fCurrentPadMapping=fgPadMapping[patch];
  fMaxHWAdd=fgkMappingSize[patch];
}

Bool_t AliHLTTPCMapping::fgMappingIsDone[fgkNofPatches]={
  kFALSE,
  kFALSE,
  kFALSE,
  kFALSE,
  kFALSE,
  kFALSE
};

const UInt_t AliHLTTPCMapping::fgkMappingSize[fgkNofPatches]={
  AliHLTTPCMapping::fgkMapping0Size,
  AliHLTTPCMapping::fgkMapping1Size,
  AliHLTTPCMapping::fgkMapping2Size,
  AliHLTTPCMapping::fgkMapping3Size,
  AliHLTTPCMapping::fgkMapping4Size,
  AliHLTTPCMapping::fgkMapping5Size
};

UInt_t AliHLTTPCMapping::fgRowMapping0[fgkMapping0Size];
UInt_t AliHLTTPCMapping::fgPadMapping0[fgkMapping0Size];
UInt_t AliHLTTPCMapping::fgHwaMapping0[fgkMappingHwaSize];
UInt_t AliHLTTPCMapping::fgRowMapping1[fgkMapping1Size];
UInt_t AliHLTTPCMapping::fgPadMapping1[fgkMapping1Size];
UInt_t AliHLTTPCMapping::fgHwaMapping1[fgkMappingHwaSize];
UInt_t AliHLTTPCMapping::fgRowMapping2[fgkMapping2Size];
UInt_t AliHLTTPCMapping::fgPadMapping2[fgkMapping2Size];
UInt_t AliHLTTPCMapping::fgHwaMapping2[fgkMappingHwaSize];
UInt_t AliHLTTPCMapping::fgRowMapping3[fgkMapping3Size];
UInt_t AliHLTTPCMapping::fgPadMapping3[fgkMapping3Size];
UInt_t AliHLTTPCMapping::fgHwaMapping3[fgkMappingHwaSize];
UInt_t AliHLTTPCMapping::fgRowMapping4[fgkMapping4Size];
UInt_t AliHLTTPCMapping::fgPadMapping4[fgkMapping4Size];
UInt_t AliHLTTPCMapping::fgHwaMapping4[fgkMappingHwaSize];
UInt_t AliHLTTPCMapping::fgRowMapping5[fgkMapping5Size];
UInt_t AliHLTTPCMapping::fgPadMapping5[fgkMapping5Size];
UInt_t AliHLTTPCMapping::fgHwaMapping5[fgkMappingHwaSize];

UInt_t* AliHLTTPCMapping::fgRowMapping[fgkNofPatches]={
  AliHLTTPCMapping::fgRowMapping0,
  AliHLTTPCMapping::fgRowMapping1,
  AliHLTTPCMapping::fgRowMapping2,
  AliHLTTPCMapping::fgRowMapping3,
  AliHLTTPCMapping::fgRowMapping4,
  AliHLTTPCMapping::fgRowMapping5
};

UInt_t* AliHLTTPCMapping::fgPadMapping[fgkNofPatches]={
  AliHLTTPCMapping::fgPadMapping0,
  AliHLTTPCMapping::fgPadMapping1,
  AliHLTTPCMapping::fgPadMapping2,
  AliHLTTPCMapping::fgPadMapping3,
  AliHLTTPCMapping::fgPadMapping4,
  AliHLTTPCMapping::fgPadMapping5
};

UInt_t* AliHLTTPCMapping::fgHwaMapping[fgkNofPatches]={
  AliHLTTPCMapping::fgHwaMapping0,
  AliHLTTPCMapping::fgHwaMapping1,
  AliHLTTPCMapping::fgHwaMapping2,
  AliHLTTPCMapping::fgHwaMapping3,
  AliHLTTPCMapping::fgHwaMapping4,
  AliHLTTPCMapping::fgHwaMapping5
};

AliHLTTPCMapping::~AliHLTTPCMapping(){
  // see header file for class documentation
}

void AliHLTTPCMapping::InitializeMap(UInt_t patch)
{
  // see header file for class documentation
  // method is deprecated, just kept as history for a while

  UInt_t fNRowsToSubtract=0;
  //The row numbers returned by digit readers used are not from zero always
  switch(patch){
  case 0:
    fNRowsToSubtract=0;
    break;
  case 1:
    fNRowsToSubtract=30;
    break;
  case 2:
    fNRowsToSubtract=0;
    break;
  case 3:
    fNRowsToSubtract=27;
    break;
  case 4:
    fNRowsToSubtract=54;
    break;
  case 5:
    fNRowsToSubtract=76;
    break;
  }

  //Making mapping arrays for the given patch
  ifstream inFile;
  TString filename;
  const char* basePath=getenv("ALICE_ROOT");
  if (basePath) {
    filename.Form("%s/TPC/mapping/Patch%d.data", basePath,patch);
  }
  inFile.open(filename.Data());
  if (!inFile) {
    HLTFatal("Unable to open file: %s      This means no mapping is provided.", filename.Data());
  }
  else{
    UInt_t nHWAdd=0;
    if(inFile >> nHWAdd){
      if(inFile >> fMaxHWAdd){
        UInt_t hwAdd=0;
        UInt_t row=0;
        UInt_t pad=0;
	switch(patch){
	case 0:
	  if(fgkMappingSize[patch]<fMaxHWAdd){
	    HLTFatal("Max hardware address exceeded for patch %d, max number is %d, number from mapping file is %d.",patch,fgkMappingSize[patch] ,fMaxHWAdd);
	    break;
	  }
	  while(inFile>>hwAdd && inFile>>row && inFile>>pad){
	    if (hwAdd>fMaxHWAdd) {
	      HLTFatal("hardware address exceeds max hwAddress %d, mapping file %s corrupted?", fMaxHWAdd);
	      break;
	    }
	    fgRowMapping0[hwAdd]=row-fNRowsToSubtract;
	    fgPadMapping0[hwAdd]=pad;
	  }
	  fgMappingIsDone[patch]=kTRUE;
	  break;
	case 1:
	  if(fgkMappingSize[patch]<fMaxHWAdd){
	    HLTFatal("Max hardware address exceeded for patch %d, max number is %d, number from mapping file is %d.",patch,fgkMappingSize[patch] ,fMaxHWAdd);
	    break;
	  }
	  while(inFile>>hwAdd && inFile>>row && inFile>>pad){
	    if (hwAdd>fMaxHWAdd) {
	      HLTFatal("hardware address exceeds max hwAddress %d, mapping file %s corrupted?", fMaxHWAdd);
	      break;
	    }
	    fgRowMapping1[hwAdd]=row-fNRowsToSubtract;
	    fgPadMapping1[hwAdd]=pad;
	  }
	  fgMappingIsDone[patch]=kTRUE;
	  break;
	case 2:
	  if(fgkMappingSize[patch]<fMaxHWAdd){
	    HLTFatal("Max hardware address exceeded for patch %d, max number is %d number from mapping file is %d.",patch,fgkMappingSize[patch] ,fMaxHWAdd);
	    break;
	  }
	  while(inFile>>hwAdd && inFile>>row && inFile>>pad){
	    if (hwAdd>fMaxHWAdd) {
	      HLTFatal("hardware address exceeds max hwAddress %d, mapping file %s corrupted?", fMaxHWAdd);
	      break;
	    }
	    fgRowMapping2[hwAdd]=row-fNRowsToSubtract;
	    fgPadMapping2[hwAdd]=pad;
	  }
	  fgMappingIsDone[patch]=kTRUE;
	  break;
	case 3:
	  if(fgkMappingSize[patch]<fMaxHWAdd){
	    HLTFatal("Max hardware address exceeded for patch %d, max number is %d number from mapping file is %d.",patch,fgkMappingSize[patch] ,fMaxHWAdd);
	    break;
	  }
	  while(inFile>>hwAdd && inFile>>row && inFile>>pad){
	    if (hwAdd>fMaxHWAdd) {
	      HLTFatal("hardware address exceeds max hwAddress %d, mapping file %s corrupted?", fMaxHWAdd);
	      break;
	    }
	    fgRowMapping3[hwAdd]=row-fNRowsToSubtract;
	    fgPadMapping3[hwAdd]=pad;
	  }
	  fgMappingIsDone[patch]=kTRUE;
	  break;
	case 4:
	  if(fgkMappingSize[patch]<fMaxHWAdd){
	    HLTFatal("Max hardware address exceeded for patch %d, max number is %d number from mapping file is %d.",patch,fgkMappingSize[patch] ,fMaxHWAdd);
	    break;
	  }
	  while(inFile>>hwAdd && inFile>>row && inFile>>pad){
	    if (hwAdd>fMaxHWAdd) {
	      HLTFatal("hardware address exceeds max hwAddress %d, mapping file %s corrupted?", fMaxHWAdd);
	      break;
	    }
	    fgRowMapping4[hwAdd]=row-fNRowsToSubtract;
	    fgPadMapping4[(UInt_t)hwAdd]=(UInt_t)pad;
	  }
	  fgMappingIsDone[patch]=kTRUE;
	  break;
	case 5:
	  if(fgkMappingSize[patch]<fMaxHWAdd){
	    HLTFatal("Max hardware address exceeded for patch %d, max number is %d number from mapping file is %d.",patch,fgkMappingSize[patch] ,fMaxHWAdd);
	    break;
	  }
	  while(inFile>>hwAdd && inFile>>row && inFile>>pad){
	    if (hwAdd>fMaxHWAdd) {
	      HLTFatal("hardware address exceeds max hwAddress %d, mapping file %s corrupted?", fMaxHWAdd);
	      break;
	    }
	    fgRowMapping5[hwAdd]=row-fNRowsToSubtract;
	    fgPadMapping5[hwAdd]=pad;
	  }
	  fgMappingIsDone[patch]=kTRUE;
	  break;
	}
      }
    }
  }
  inFile.close();
}

UInt_t AliHLTTPCMapping::GetRow(UInt_t hwadd) const
{
  // see header file for class documentation
  assert(fCurrentRowMapping);
  if (!fCurrentRowMapping) return 0;
  if (hwadd>fMaxHWAdd) return 0;
  return fCurrentRowMapping[hwadd];
}

UInt_t AliHLTTPCMapping::GetPad(UInt_t hwadd) const
{
  // see header file for class documentation
  assert(fCurrentPadMapping);
  if (!fCurrentPadMapping) return 0;
  if (hwadd>fMaxHWAdd) return 0;
  return fCurrentPadMapping[hwadd];
}

UInt_t AliHLTTPCMapping::GetHwAddress(UInt_t row, UInt_t pad) const
{
  // see header file for class documentation
  assert(fPatch<fgkNofPatches);
  assert(fgHwaMapping[fPatch]);
  UInt_t hwaIndex=(row<<fgkMappingHwaRowBitShift) + pad;
  if (hwaIndex<fgkMappingHwaSize) {
    return (fgHwaMapping[fPatch])[hwaIndex];
  }
  return ~((UInt_t)0);
}

Bool_t AliHLTTPCMapping::ReadArray(UInt_t patch, UInt_t arraySize, UInt_t rowArray[], UInt_t padArray[], UInt_t hwaMappingSize, UInt_t hwaArray[]) const
{
  // see header file for class documentation
  //Making mapping arrays for the given patch
  Bool_t result=kTRUE;

  assert(rowArray!=NULL || padArray!=NULL || hwaArray!=NULL);
  if (!rowArray || !padArray || !hwaArray) return kFALSE;

  memset(rowArray, -1, arraySize*sizeof(rowArray[0]));
  memset(padArray, -1, arraySize*sizeof(padArray[0]));
  memset(hwaArray, -1, hwaMappingSize*sizeof(padArray[0]));

  ifstream inFile;
  TString filename;
  const char* basePath=getenv("ALICE_ROOT");
  if (basePath) {
    filename.Form("%s/TPC/mapping/Patch%d.data", basePath,patch);
  }
  inFile.open(filename.Data());
  if (!inFile) {
    HLTFatal("Unable to open file: %s      This means no mapping is provided.", filename.Data());
    return kFALSE;
  }

  UInt_t nHWAdd=0;
  UInt_t maxHWAdd=0;
  if(inFile >> nHWAdd){
    if(inFile >> maxHWAdd){
      if(arraySize<maxHWAdd){
	HLTFatal("Max hardware address exceeded for patch %d, max number is %d, number from mapping file is %d.",patch, arraySize ,fMaxHWAdd);
	result=kFALSE;
      }
      UInt_t hwAdd=0;
      UInt_t row=0;
      UInt_t pad=0;
      while(result && inFile>>hwAdd && inFile>>row && inFile>>pad){
	if (hwAdd>maxHWAdd) {
	  HLTFatal("hardware address exceeds max hwAddress %d, mapping file %s corrupted?", maxHWAdd);
	  result=kFALSE;
	  break;
	}
	Int_t dummy=0;
	// AliHLTTPCTransform::GetFirstRow returns first row in scheme A.
	// We have to transform to scheme B by AliHLTTPCTransform::Slice2Sector.
	Int_t offsetSchemeB=0;
	AliHLTTPCTransform::Slice2Sector(0, fRowOffset, dummy, offsetSchemeB);
	rowArray[hwAdd]=row-offsetSchemeB;
	padArray[hwAdd]=pad;
	UInt_t hwaIndex=(rowArray[hwAdd]<<fgkMappingHwaRowBitShift) + padArray[hwAdd];
	assert(hwaIndex<hwaMappingSize);
	if (hwaIndex<hwaMappingSize) {
	  hwaArray[hwaIndex]=hwAdd;
	}
      }
    }
  }
  inFile.close();
  return result;
}

Bool_t AliHLTTPCMapping::IsValidHWAddress(UInt_t hwadd) const
{
  if (hwadd>fMaxHWAdd){
    return kFALSE;
  }
  else{
    return kTRUE;
  }
}
