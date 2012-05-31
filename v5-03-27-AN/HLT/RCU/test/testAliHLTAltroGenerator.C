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

/** @file   testAliHLTAltroGenerator.C
    @author Matthias Richter
    @date   
    @brief  Test macro/program for the AliHLTAltroGenerator
 */

#ifndef __CINT__
#include "TSystem.h"
#include "AliHLTSystem.h"
#include "AliRawDataHeader.h"
#include "AliAltroDecoder.h"
#include "AliAltroData.h"
#include "AliAltroBunch.h"
#include "AliHLTAltroGenerator.h"
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

// printouts or not
const bool bVerbose=false;

// some defaults
const int maxChannels=1000;
const int maxBunches=50;
const int maxBunchLength=10;
const int maxTimebin=1024;
const int maxSignal=1024;

/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////

int testAliHLTAltroGenerator()
{
  int iResult=0;
#ifdef __CINT__
  string ld_library_path="../.libs:";
  ld_library_path+=gSystem->GetDynamicPath();
  gSystem->SetDynamicPath(ld_library_path.c_str());
  gSystem->Load("libAliHLTRCU.so");
#endif
  AliHLTSystem gHLT;

  AliHLTAltroGenerator g(maxChannels, maxBunches, maxBunchLength, maxTimebin, maxSignal);
  //g.SetDirection(AliHLTAltroGenerator::kForwards);
  if ((iResult=g.Generate())<0) return iResult;

  if (bVerbose) {
    cout << "***************************************************************" << endl;
    cout << "************** Dumping simulated Altro data *******************" << endl;
    g.Print();
    cout << endl;
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  if (bVerbose) {
    g.Rewind();
    cout << "***************************************************************" << endl;
    cout << "********************** reading bunch model  *******************" << endl;
    while (iResult>=0 && g.NextChannel()) {
    cout << "***************************************************************" << endl;
      cout << "channel address: " << g.GetHwAddress() << "    " << g.GetBunchCount() << " bunch(es)" << endl;

      while (iResult>=0 && g.NextBunch()) {
	int bunchLength=g.GetBunchSize();
	cout << "   length " << bunchLength << " start time " << g.GetStartTime() << ":     ";
	const Short_t* pData=g.GetSignals();
	while (bunchLength-->0 && pData) {
	  cout << " " << *pData++;
	}
	cout << "      -> end time " << g.GetEndTime() << endl;
      }
    }
    cout << endl;
  }

  if (bVerbose) {
    g.Rewind();
    int lastChannel=-1;
    int lastTime=-1;
    cout << "***************************************************************" << endl;
    cout << "********************** reading stream model *******************" << endl;
    while (iResult>=0 && g.Next()) {
      if (lastTime>=0 && lastTime!=g.GetStartTime()+1 && lastTime!=g.GetStartTime()-1)
	cout << endl;

      if (lastChannel<0 || lastChannel!=g.GetHwAddress()) {
	cout << "***************************************************************" << endl;
	cout << "channel address: " << g.GetHwAddress() << endl;
      }

      if (lastTime<0 || (lastTime!=g.GetStartTime()+1 && lastTime!=g.GetStartTime()-1))
	cout << " time " << g.GetStartTime() << ":     ";

      cout << " " << g.GetSignal();

      lastChannel=g.GetHwAddress();
      lastTime=g.GetStartTime();
    }
    cout << endl;
    cout << endl;
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  AliRawDataHeader cdh;
  g.SetCDH(&cdh, 32);

  UInt_t trailer=0;
  g.SetRCUTrailer((UChar_t*)&trailer, 4);

  UChar_t* pBuffer=NULL;
  Int_t size=g.GetData(pBuffer);

  /*
  ios::openmode filemode=(ios::openmode)0;
  ofstream rawfile("/tmp/altro-enc.dat", filemode);
  if (rawfile.good()) {
    rawfile.write(reinterpret_cast<const char*>(pBuffer), size);
  }
  rawfile.close();
  */

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  // can not have a static AltroDecoder, crash when function
  // is called. I had a similar problem in the AliHLTAltroChannelSelectorComponent

  if (bVerbose) {
    cout << "***************************************************************" << endl;
    cout << "********************* comparing encoded data ******************" << endl;
    cout << "***************************************************************" << endl;
  }

  g.Rewind();
  AliAltroDecoder* decoder=new AliAltroDecoder;
  if (iResult>=0 && decoder->SetMemory(pBuffer, size)<0) {
    cout << "error setting up decoder " << endl;
    iResult=-1;
  }

  if (iResult>=0 && !decoder->Decode()) {
    cout << "error decoding data" << endl;
    iResult=-1;    
  }

  AliAltroData altrochannel;
  while (iResult>=0 && decoder->NextChannel(&altrochannel)) {
    if (!g.NextChannel()) {
      cout << "error getting next simulated channel" << endl;
      iResult=-1;
      break;
    }
    int hwadd=altrochannel.GetHadd();
    if (hwadd!=g.GetHwAddress()) {
      cout << "channel address missmatch: simulated " << g.GetHwAddress() << " encoded " << hwadd << endl;
      iResult=-1;
      break;
    }

    if (bVerbose) cout << "comparing channel " << hwadd << endl;

    AliAltroBunch altrobunch;
    while (iResult>=0 && altrochannel.NextBunch(&altrobunch)) {
      if (!g.NextBunch()) {
	cout << "error getting bunch in simulated data" <<endl;
	iResult=-1;
	break;
      }
      int bunchLength=altrobunch.GetBunchSize();
      if (bunchLength!=(int)g.GetBunchSize()) {
	cout << "bunch length missmatch: simulated " << g.GetBunchSize() << " encoded " << bunchLength << hex << " (" << bunchLength << ")" << dec << endl;
	iResult=-1;
	break;
      }
      int bunchEndTime=altrobunch.GetEndTimeBin();
	if (bunchEndTime!=(int)g.GetEndTime()) {
	cout << "bunch end time missmatch: simulated " << g.GetEndTime() << " encoded " << bunchEndTime << endl;
	iResult=-1;
	break;
      }
      if (bVerbose) cout << " bunch length " << bunchLength << ", end time " << bunchEndTime << endl;
      const  UInt_t* bunchData=altrobunch.GetData();
      const  Short_t* simData=g.GetSignals();
      for (int bin=0; bin<bunchLength; bin++) {
	if ((Short_t)bunchData[bin]!=simData[bin]) {
	  cout << "data missmatch at bunch position " << bin << " : simulated " << simData[bin] << " encoded " << bunchData[bin] << endl;
	  iResult=-1;
	  break;
	}
      }
    }

  }
  
  delete decoder;

  return iResult;
}

int main(int /*argc*/, const char** /*argv*/)
{
  int iResult=0;
  // this test takes ~20 times longer than the testAliHLTAltroEncoder
  // no clue why, has to be traced with vtune
  int iCount=500;
  for (int i=0; i<iCount; i++) {
    if ((iResult=testAliHLTAltroGenerator())<0) {
      cout << "missmatch in cycle no " << i << endl;
      return iResult;
    }
  }
  cout << "checking: "<< iCount << " encoding cycle(s) successfully tested" << endl;
  return 0;
}
