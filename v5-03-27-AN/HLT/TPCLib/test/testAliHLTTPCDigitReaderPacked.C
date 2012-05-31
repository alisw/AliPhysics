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

/** @file   testAliHLTTPCDigitReaderPacked.C
    @author Matthias Richter
    @date   
    @brief  Test macro/program for the AliHLTTPCDigitReaderPacked
 */

#ifndef __CINT__
#include "TSystem.h"
#include "AliHLTSystem.h"
#include "AliRawDataHeader.h"
#include "AliHLTAltroGenerator.h"
#include "AliHLTTPCDigitReaderPacked.h"
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
const bool bVerbose=true;

// some defaults
const int maxChannels=10;
const int maxBunches=10;
const int maxBunchLength=10;
const int maxTimebin=1024;
const int maxSignal=1024;

/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////

int testAliHLTTPCDigitReaderPacked()
{
  int iResult=0;
#ifdef __CINT__
  gSystem->Load("libAliHLTUtil.so");
  gSystem->Load("libAliHLTRCU.so");
  gSystem->Load("libAliHLTTPC.so");
#endif
  AliHLTSystem gHLT;

  AliHLTAltroGenerator generator(maxChannels, maxBunches, maxBunchLength, maxTimebin, maxSignal);
  //generator.SetDirection(AliHLTAltroGenerator::kForwards);
  if ((iResult=generator.Generate())<0) return iResult;

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  if (bVerbose) {
    cout << "***************************************************************" << endl;
    cout << "************** Dumping simulated Altro data *******************" << endl;
    generator.Print();
    cout << endl;
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  AliRawDataHeader cdh;
  generator.SetCDH(&cdh, 32);

  UInt_t trailer=0;
  generator.SetRCUTrailer((UChar_t*)&trailer, 4);

  UChar_t* pBuffer=NULL;
  Int_t size=generator.GetData(pBuffer);

  int partition=0;
  if (bVerbose) {
    AliHLTTPCDigitReaderPacked decoder;
    decoder.SetOldRCUFormat(true);
    decoder.SetUnsorted(true);
    if ((iResult=decoder.InitBlock(pBuffer, size, partition, 0))>=0) {
      cout << "***************************************************************" << endl;
      cout << "********************** reading bunch model  *******************" << endl;
      while (iResult>=0 && decoder.NextChannel()) {
	cout << "***************************************************************" << endl;
	cout << "channel address: " << decoder.GetAltroBlockHWaddr() << endl;

	while (iResult>=0 && decoder.NextBunch()) {
	  int bunchLength=decoder.GetBunchSize();
	  cout << "   length " << bunchLength << " time " << decoder.GetTime() << ":     ";
	  const UInt_t* pData=decoder.GetSignals();
	  while (bunchLength-->0 && pData) {
	    cout << " " << *pData++;
	  }
	  cout << endl;
	}
      }
      cout << endl;
    }
  }   

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  if (bVerbose) {
    AliHLTTPCDigitReaderPacked decoder;
    decoder.SetOldRCUFormat(true);
    decoder.SetUnsorted(true);
    if ((iResult=decoder.InitBlock(pBuffer, size, partition, 0))>=0) {
      int lastChannel=-1;
      int lastTime=-1;
      cout << "***************************************************************" << endl;
      cout << "********************** reading stream model *******************" << endl;
      while (iResult>=0 && decoder.Next()) {
	if (lastTime>=0 && lastTime!=decoder.GetTime()+1 && lastTime!=decoder.GetTime()-1)
	  cout << endl;
	
	if (lastChannel<0 || lastChannel!=(int)decoder.GetAltroBlockHWaddr()) {
	  cout << "***************************************************************" << endl;
	  cout << "channel address: " << decoder.GetAltroBlockHWaddr() << endl;
	}

	if (lastTime<0 || (lastTime!=decoder.GetTime()+1 && lastTime!=decoder.GetTime()-1))
	  cout << " time " << decoder.GetTime() << ":     ";

	cout << " " << decoder.GetSignal();

	lastChannel=decoder.GetAltroBlockHWaddr();
	lastTime=decoder.GetTime();
      }
      cout << endl;
      cout << endl;
    }    
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  if (bVerbose) {
    cout << "***************************************************************" << endl;
    cout << "********************* comparing encoded data ******************" << endl;
    cout << "***************************************************************" << endl;
  }

  generator.Rewind();
  AliHLTTPCDigitReaderPacked decoder;
  decoder.SetUnsorted(true);
  if (iResult>=0) iResult=decoder.InitBlock(pBuffer, size, partition, 0);
  while (iResult>=0 && decoder.NextChannel()) {
    if (!generator.NextChannel()) {
      cout << "error getting next simulated channel" << endl;
      iResult=-1;
      break;
    }
    int hwadd=decoder.GetAltroBlockHWaddr();
    if (hwadd!=generator.GetHwAddress()) {
      cout << "channel address missmatch: simulated " << generator.GetHwAddress() << " encoded " << hwadd << endl;
      iResult=-1;
      break;
    }

    if (bVerbose) cout << "comparing channel " << hwadd << endl;

    while (iResult>=0 && decoder.NextBunch()) {
      if (!generator.NextBunch()) {
	cout << "error getting bunch in simulated data" <<endl;
	iResult=-1;
	break;
      }
      int bunchLength=decoder.GetBunchSize();
      if (bunchLength!=(int)generator.GetBunchSize()) {
	cout << "bunch length missmatch: simulated " << generator.GetBunchSize() << " encoded " << bunchLength << hex << " (" << bunchLength << ")" << dec << endl;
	iResult=-1;
	break;
      }
      int bunchStartTime=decoder.GetTime();
	if (bunchStartTime!=(int)generator.GetStartTime()) {
	cout << "bunch end time missmatch: simulated " << generator.GetStartTime() << " encoded " << bunchStartTime << endl;
	iResult=-1;
	break;
      }
      if (bVerbose) cout << " bunch length " << bunchLength << ", end time " << bunchStartTime << endl;
      const  UInt_t* bunchData=decoder.GetSignals();
      const  Short_t* simData=generator.GetSignals();
      for (int bin=0; bin<bunchLength; bin++) {
	if ((Short_t)bunchData[bin]!=simData[bin]) {
	  cout << "data missmatch at bunch position " << bin << " : simulated " << simData[bin] << " encoded " << bunchData[bin] << endl;
	  iResult=-1;
	  break;
	}
      }
    }
  }

  return 0;
}

int main(int /*argc*/, const char** /*argv*/)
{
  int iResult=0;
  int iCount=1;
  for (int i=0; i<iCount; i++) {
    if ((iResult=testAliHLTTPCDigitReaderPacked())<0) {
      cout << "missmatch in cycle no " << i << endl;
      return iResult;
    }
  }
  cout << "checking: "<< iCount << " encoding cycle(s) successfully tested" << endl;
  return 0;
}
