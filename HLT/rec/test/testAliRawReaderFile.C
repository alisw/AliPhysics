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

/** @file   testAliRawReaderMemory.C
    @author Matthias Richter
    @date   
    @brief  Test macro/program for the AliRawReaderMemory
 */

#ifndef __CINT__
#include "TFile.h"
#include "TDatime.h"
#include "TRandom.h"
#include "TArrayI.h"
#include "TArrayC.h"
#include "TSystem.h"
#include "AliRawDataHeader.h"
#include "AliRawReaderFile.h"
#include "AliDAQ.h"
#include "AliHLTSystem.h"
#include <ostream>
#endif //__CINT__

#ifndef __CINT__
const int sizeofAliRawDataHeader=sizeof(AliRawDataHeader);
#else
// cint does not handle sizeof correctly
const int sizeofAliRawDataHeader=32;
#endif

/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
//
// configuration of the test program
//
const char* tmpdir="/tmp";
const int nofDDLs=2;


/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
int GetRandom(int min, int max);
int FillRandomDDL(TArrayC& target);
int WriteDDL(TArrayC& ddl, int ddlid);
int CheckDDL_ReadNext(AliRawReader* pRawReader, TArrayC& ddl, int ddlid);
int CheckDDL_ReadNextData(AliRawReader* pRawReader, TArrayC& ddl, int ddlid);
int CheckDDL_ReadNextChar(AliRawReader* pRawReader, TArrayC& ddl, int ddlid);

int testAliRawReaderMemory()
{
  int iResult=0;
  TArrayC ddl[nofDDLs];
  int ddlid[nofDDLs]={768, 769};

  for (int i=0; i<nofDDLs; i++) {
    FillRandomDDL(ddl[i]);
    WriteDDL(ddl[i], ddlid[i]);
  }

  TString rawdir=tmpdir;
  if (!rawdir.EndsWith("/")) rawdir+="/";
  AliRawReader* pRawReader=AliRawReader::Create(rawdir);
  if (!pRawReader) {
    cerr << "can not create RawReader" << endl;
    return -1;
  }

  if (!pRawReader->NextEvent()) {
    cerr << "error: getting event failed" << endl;
    return -1;
  }

  // 1st pass: check AliRawReader::ReadNext()
  while (pRawReader->ReadHeader()) {
    int id=pRawReader->GetEquipmentId();
    int i=0;
    for (i=0; i<nofDDLs; i++) {
      if (ddlid[i]==id) break;
    }
    if (i==nofDDLs) {
      cerr << "error: can not find ddl id " << id << endl;
      return -1;
    }

    if ((iResult=CheckDDL_ReadNext(pRawReader, ddl[i], ddlid[i]))<0) return iResult;
  }

  // 2nd pass: check AliRawReader::ReadNextData()
  pRawReader->Reset();
  while (pRawReader->ReadHeader()) {
    int id=pRawReader->GetEquipmentId();
    int i=0;
    for (i=0; i<nofDDLs; i++) {
      if (ddlid[i]==id) break;
    }
    if (i==nofDDLs) {
      cerr << "error: can not find ddl id " << id << endl;
      return -1;
    }

    if ((iResult=CheckDDL_ReadNextData(pRawReader, ddl[i], ddlid[i]))<0) return iResult;
  }

  // 3rd pass: check AliRawReader::ReadNextData()
  pRawReader->Reset();
  while (pRawReader->ReadHeader()) {
    int i=0;
    do {
      int id=pRawReader->GetEquipmentId();
      for (i=0; i<nofDDLs; i++) {
	if (ddlid[i]==id) break;
      }
      if (i==nofDDLs) {
	cerr << "error: can not find ddl id " << id << endl;
	return -1;
      }

      if ((iResult=CheckDDL_ReadNextChar(pRawReader, ddl[i], ddlid[i]))<0) return iResult;
    }
    // note: the ReadHeader is hidden in the ReadNextChar
    while (iResult!=ddlid[i]);
  }

  delete pRawReader;
  return 0;
}

Bool_t seedSet=kFALSE;

/**
 * Get a random number in the given range.
 */
int GetRandom(int min, int max)
{
  if (max-min<2) return min;
  static TRandom rand;
  if (!seedSet) {
    TDatime dt;
    rand.SetSeed(dt.Get());
    seedSet=kTRUE;
  }
  return rand.Integer(max-min);
}

int FillRandomDDL(TArrayC& target)
{
  int size=0;
  do size=GetRandom(100, 10000);
  while (size<100);

  target.Set(size);
  *((Int_t*)target.GetArray())=size;
  for (int i=sizeofAliRawDataHeader; i<size; i++) {
    Int_t data=GetRandom(0,255);
    target.AddAt((UChar_t)data, i);
  }
  return size;
}

int WriteDDL(TArrayC& ddl, int ddlid)
{
  Int_t dummy=0;
  TString filename;
  filename.Form("%s/raw0", tmpdir);
  gSystem->mkdir(filename);
  filename.Form("%s/%s_%d.ddl", filename.Data(), AliDAQ::DetectorNameFromDdlID(ddlid, dummy), ddlid);
  FILE* fp=fopen(filename.Data(), "w");
  if (!fp) return -1;

  fwrite(ddl.GetArray(), 1, ddl.GetSize(), fp);
  fclose(fp);

  return 0;
}

int CheckDDL_ReadNext(AliRawReader* pRawReader, TArrayC& ddl, int ddlid)
{
  int dataSize=pRawReader->GetDataSize();
  if (ddl.GetSize()!=dataSize+sizeofAliRawDataHeader) {
    cerr << "error: size mismatch in ddl " << ddlid << ": " << dataSize+sizeofAliRawDataHeader << " required " << ddl.GetSize() << endl;
    return -1;
  }
  TArrayC buffer(dataSize);
  UChar_t* pTgt=(UChar_t*)buffer.GetArray();
  
  if (!pRawReader->ReadNext(pTgt, buffer.GetSize())) {
    cerr << "error: reading " << buffer.GetSize() << " byte(s) from ReadNext (ddl " << ddlid << ")" << endl;
    return -1;
  }

  cout << "verify ReadNext: ddl " << ddlid << endl;
  if (memcmp(ddl.GetArray()+sizeofAliRawDataHeader, buffer.GetArray(), buffer.GetSize())!=0) {
    cerr << "error: verification of ddl " << ddlid << ")" << endl;
    return -1;
  }
  return 0;
}

int CheckDDL_ReadNextData(AliRawReader* pRawReader, TArrayC& ddl, int ddlid)
{
  int dataSize=pRawReader->GetDataSize();
  if (ddl.GetSize()!=dataSize+sizeofAliRawDataHeader) {
    cerr << "error: size mismatch in ddl " << ddlid << ": " << dataSize+sizeofAliRawDataHeader << " required " << ddl.GetSize() << endl;
    return -1;
  }

  UChar_t* pTgt=NULL;  
  if (!pRawReader->ReadNextData(pTgt) || pTgt==NULL) {
    cerr << "error: reading " << dataSize << " byte(s) from ReadNext (ddl " << ddlid << ")" << endl;
    return -1;
  }

  cout << "verify ReadNextData: ddl " << ddlid << endl;
  if (memcmp(ddl.GetArray()+sizeofAliRawDataHeader, pTgt, dataSize)!=0) {
    cerr << "error: verification of ddl " << ddlid << endl;
    return -1;
  }
  return 0;
}

int CheckDDL_ReadNextChar(AliRawReader* pRawReader, TArrayC& ddl, int ddlid)
{
  int iResult=ddlid;
  int dataSize=pRawReader->GetDataSize();
  if (ddl.GetSize()!=dataSize+sizeofAliRawDataHeader) {
    cerr << "error: size mismatch in ddl " << ddlid << ": " << dataSize+sizeofAliRawDataHeader << " required " << ddl.GetSize() << endl;
    return -1;
  }

  cout << "verify ReadNextChar: ddl " << ddlid << endl;
  static int offset=0;
  int i=sizeofAliRawDataHeader+offset;
  UChar_t data=0;
  Bool_t haveData=kFALSE;
  for (; (haveData=pRawReader->ReadNextChar(data)) && i<ddl.GetSize(); i++) {
    if (data!=(UChar_t)ddl.At(i)) {
      cerr << "error: at position " << i << " of ddl " << ddlid << ": read " << (UShort_t)data << " required " << (UShort_t)(UChar_t)ddl.At(i) << endl;
      return -1;
    }
  }

  if (i<ddl.GetSize()) {
    cerr << "error: reading data of ddl " << ddlid << ": " << i << " of " << ddl.GetSize() << endl;
    return -1;
  }
  if (pRawReader->GetEquipmentId()>=0 && ddlid!=pRawReader->GetEquipmentId()) {
    //cerr << "error: ddl id missmatch, expecting " << ddlid << " got " << pRawReader->GetEquipmentId() << endl;
    // thats not an error condition, RawReader just changes silently to the next DDL
    iResult=pRawReader->GetEquipmentId();
    offset=1;
  } else {
    offset=0;
  }
  if (haveData && iResult==ddlid) {
    cerr << "error: size missmatch in ddl " << ddlid << ": still data available after " << ddl.GetSize() << " byte(s)"<< endl;
    return -1;
  }
  return iResult;
}

int main(int /*argc*/, const char** /*argv*/)
{
  return testAliRawReaderMemory();
}
