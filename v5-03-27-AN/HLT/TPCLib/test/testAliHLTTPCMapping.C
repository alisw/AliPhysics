// $Id$

/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
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

/** @file   testAliHLTTPCMapping.C
    @author Matthias Richter
    @date   
    @brief  Test macro/program for the AliHLTTPCMapping class
 */

#ifndef __CINT__
#include "TSystem.h"
#include "AliHLTTPCMapping.h"
#include "AliHLTTPCTransform.h"
#include <ostream>
#include <istream>
#endif //__CINT__

/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
//
// configuration of the test program
//

// printouts or not
const bool bVerbose=true;

/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////

class AliHLTTPCMapping;
bool compareMapping(int patch, AliHLTTPCMapping* mapper);

int testAliHLTTPCMapping()
{
  int iResult=0;
#ifdef __CINT__
  gSystem->Load("libAliHLTUtil.so");
  gSystem->Load("libAliHLTRCU.so");
  gSystem->Load("libAliHLTTPC.so");
#endif
  //AliHLTSystem gHLT;

  const int nofMappers=6;
  AliHLTTPCMapping* mappers[nofMappers];
  AliHLTTPCMapping* mappers2[nofMappers];
  for (int i=0; i<nofMappers; i++) {
    mappers[i]=NULL;
    mappers2[i]=NULL;
  }

  // create mappers
  for (int i=0; i<nofMappers; i++) {
    mappers[i]=new AliHLTTPCMapping(i);
  }

  // check mappers
  for (int i=0; i<nofMappers; i++) {
    if (!compareMapping(i, mappers[i])) {
      iResult=-1;
    }
  }
  cout << "checking: 1st instance:";
  if (iResult<0) 
    cout << " failed" << endl;
  else {
    cout << " ok" << endl;

    // create 2nd instance
    for (int i=0; i<nofMappers; i++) {
      mappers2[i]=new AliHLTTPCMapping(i);
    }

    // check 2nd instance
    for (int i=0; i<nofMappers; i++) {
      if (!compareMapping(i, mappers2[i])) {
	iResult=-1;
      }
    }
    cout << "checking: 2nd instance:";
    if (iResult<0) 
      cout << " failed" << endl;
    else
      cout << " ok" << endl;
  }

  // delete mappers
  for (int i=0; i<nofMappers; i++) {
    if (mappers[i]) delete mappers[i];
    if (mappers2[i]) delete mappers2[i];
  }
  return iResult;
}

bool compareMapping(int patch, AliHLTTPCMapping* mapper)
{
  bool result=true;
  if (!mapper) return false;
  ifstream inFile;
  TString filename;
  const char* basePath=getenv("ALICE_ROOT");
  if (basePath) {
    filename.Form("%s/TPC/mapping/Patch%d.data", basePath,patch);
  }
  inFile.open(filename.Data());
  if (!inFile) {
    cout << "Unable to open file: " << filename << endl;
    return false;
  }

  UInt_t nHWAdd=0;
  UInt_t maxHWAdd=0;
  UInt_t hwAdd=0;
  UInt_t row=0;
  UInt_t pad=0;
  Int_t dummy=0;
  Int_t rowOffset=0;
  result=AliHLTTPCTransform::Slice2Sector(0, AliHLTTPCTransform::GetFirstRow(patch), dummy, rowOffset);
  
  if(inFile >> nHWAdd && inFile >> maxHWAdd) {
    while(result && inFile>>hwAdd && inFile>>row && inFile>>pad){
      row-=rowOffset;
      if (row!=mapper->GetRow(hwAdd) || pad!=mapper->GetPad(hwAdd)) {
	cout << "channel to row/pad mapping mismatch at channel " << hwAdd << ": expected " << row << "/" << pad << "  got " << mapper->GetRow(hwAdd) << "/" << mapper->GetPad(hwAdd) << endl;
	result=false;
	break;
      }
      if (hwAdd!=mapper->GetHwAddress(row, pad)) {
	cout << "row/pad to channel mapping mismatch for " << row << "/" << pad << ": expected channel:" << hwAdd << "  got " << mapper->GetHwAddress(row,pad) << endl;
	result=false;
	break;
      }
    }
  }
  inFile.close();
  return result;  
}

int main(int /*argc*/, const char** /*argv*/)
{
  return testAliHLTTPCMapping();
}
