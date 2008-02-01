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

ClassImp(AliHLTTPCMapping)

AliHLTTPCMapping::AliHLTTPCMapping(UInt_t patch)
  :
  fNHWAdd(0),
  fMaxHWAdd(0),
  fCurrentRowMapping(NULL),
  fCurrentPadMapping(NULL)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

  switch(patch){
  case 0:
    if(!fMapping0IsDone){
      memset(fgRowMapping0, 0, fgkMapping0Size*sizeof(UInt_t));
      memset(fgPadMapping0, 0, fgkMapping0Size*sizeof(UInt_t));
      InitializeMap(patch);
    }
    fCurrentRowMapping=fgRowMapping0;
    fCurrentPadMapping=fgPadMapping0;
    break;
  case 1:
    if(!fMapping1IsDone){
      memset(fgRowMapping1, 0, fgkMapping1Size*sizeof(UInt_t));
      memset(fgPadMapping1, 0, fgkMapping1Size*sizeof(UInt_t));
      InitializeMap(patch);
    }
    fCurrentRowMapping=fgRowMapping1;
    fCurrentPadMapping=fgPadMapping1;
    break;
  case 2:
    if(!fMapping2IsDone){
      memset(fgRowMapping2, 0, fgkMapping2Size*sizeof(UInt_t));
      memset(fgPadMapping2, 0, fgkMapping2Size*sizeof(UInt_t));
      InitializeMap(patch);
    }
    fCurrentRowMapping=fgRowMapping2;
    fCurrentPadMapping=fgPadMapping2;
    break;
  case 3:
    if(!fMapping3IsDone){
      memset(fgRowMapping3, 0, fgkMapping3Size*sizeof(UInt_t));
      memset(fgPadMapping3, 0, fgkMapping3Size*sizeof(UInt_t));
      InitializeMap(patch);
    }
    fCurrentRowMapping=fgRowMapping3;
    fCurrentPadMapping=fgPadMapping3;
    break;
  case 4:
    if(!fMapping4IsDone){
      memset(fgRowMapping4, 0, fgkMapping4Size*sizeof(UInt_t));
      memset(fgPadMapping4, 0, fgkMapping4Size*sizeof(UInt_t));
      InitializeMap(patch);
    }
    fCurrentRowMapping=fgRowMapping4;
    fCurrentPadMapping=fgPadMapping4;
    break;
  case 5:
    if(!fMapping5IsDone){
      memset(fgRowMapping5, 0, fgkMapping5Size*sizeof(UInt_t));
      memset(fgPadMapping5, 0, fgkMapping5Size*sizeof(UInt_t));
      InitializeMap(patch);
    }
    fCurrentRowMapping=fgRowMapping5;
    fCurrentPadMapping=fgPadMapping5;
    break;
  }
}

Bool_t AliHLTTPCMapping::fMapping0IsDone=kFALSE;
Bool_t AliHLTTPCMapping::fMapping1IsDone=kFALSE;
Bool_t AliHLTTPCMapping::fMapping2IsDone=kFALSE;
Bool_t AliHLTTPCMapping::fMapping3IsDone=kFALSE;
Bool_t AliHLTTPCMapping::fMapping4IsDone=kFALSE;
Bool_t AliHLTTPCMapping::fMapping5IsDone=kFALSE;
UInt_t AliHLTTPCMapping::fgRowMapping0[fgkMapping0Size];
UInt_t AliHLTTPCMapping::fgPadMapping0[fgkMapping0Size];
UInt_t AliHLTTPCMapping::fgRowMapping1[fgkMapping1Size];
UInt_t AliHLTTPCMapping::fgPadMapping1[fgkMapping1Size];
UInt_t AliHLTTPCMapping::fgRowMapping2[fgkMapping2Size];
UInt_t AliHLTTPCMapping::fgPadMapping2[fgkMapping2Size];
UInt_t AliHLTTPCMapping::fgRowMapping3[fgkMapping3Size];
UInt_t AliHLTTPCMapping::fgPadMapping3[fgkMapping3Size];
UInt_t AliHLTTPCMapping::fgRowMapping4[fgkMapping4Size];
UInt_t AliHLTTPCMapping::fgPadMapping4[fgkMapping4Size];
UInt_t AliHLTTPCMapping::fgRowMapping5[fgkMapping5Size];
UInt_t AliHLTTPCMapping::fgPadMapping5[fgkMapping5Size];

AliHLTTPCMapping::~AliHLTTPCMapping(){
  // see header file for class documentation
}

void AliHLTTPCMapping::InitializeMap(UInt_t patch)
{
  // see header file for class documentation

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
    fNHWAdd=0;
    if(inFile >> fNHWAdd){
      if(inFile >> fMaxHWAdd){
        UInt_t hwAdd=0;
        UInt_t row=0;
        UInt_t pad=0;
	switch(patch){
	case 0:
	  if(fgkMapping0Size<fMaxHWAdd){
	    HLTFatal("Max hardware address exceeded for patch %d, max number is %d, number from mapping file is %d.",patch,fgkMapping0Size ,fMaxHWAdd);
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
	  fMapping0IsDone=kTRUE;
	  break;
	case 1:
	  if(fgkMapping1Size<fMaxHWAdd){
	    HLTFatal("Max hardware address exceeded for patch %d, max number is %d, number from mapping file is %d.",patch,fgkMapping1Size ,fMaxHWAdd);
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
	  fMapping1IsDone=kTRUE;
	  break;
	case 2:
	  if(fgkMapping2Size<fMaxHWAdd){
	    HLTFatal("Max hardware address exceeded for patch %d, max number is %d number from mapping file is %d.",patch,fgkMapping2Size ,fMaxHWAdd);
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
	  fMapping2IsDone=kTRUE;
	  break;
	case 3:
	  if(fgkMapping3Size<fMaxHWAdd){
	    HLTFatal("Max hardware address exceeded for patch %d, max number is %d number from mapping file is %d.",patch,fgkMapping3Size ,fMaxHWAdd);
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
	  fMapping3IsDone=kTRUE;
	  break;
	case 4:
	  if(fgkMapping4Size<fMaxHWAdd){
	    HLTFatal("Max hardware address exceeded for patch %d, max number is %d number from mapping file is %d.",patch,fgkMapping4Size ,fMaxHWAdd);
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
	  fMapping4IsDone=kTRUE;
	  break;
	case 5:
	  if(fgkMapping5Size<fMaxHWAdd){
	    HLTFatal("Max hardware address exceeded for patch %d, max number is %d number from mapping file is %d.",patch,fgkMapping5Size ,fMaxHWAdd);
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
	  fMapping5IsDone=kTRUE;
	  break;
	}
      }
    }
  }
  inFile.close();
}

UInt_t AliHLTTPCMapping::GetRow(UInt_t hwadd)
{
  // see header file for class documentation
  assert(fCurrentRowMapping);
  assert(hwadd<=fMaxHWAdd);
  if (!fCurrentRowMapping) return 0;
  if (hwadd>fMaxHWAdd) return 0;
  return fCurrentRowMapping[hwadd];
}

UInt_t AliHLTTPCMapping::GetPad(UInt_t hwadd)
{
  // see header file for class documentation
  assert(fCurrentPadMapping);
  assert(hwadd<=fMaxHWAdd);
  if (!fCurrentPadMapping) return 0;
  if (hwadd>fMaxHWAdd) return 0;
  return fCurrentPadMapping[hwadd];
}
