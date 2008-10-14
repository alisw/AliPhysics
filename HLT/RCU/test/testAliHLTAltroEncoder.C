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

/** @file   testAliHLTAltroEncoder.C
    @author Matthias Richter
    @date   
    @brief  Test macro/program for the AliHLTAltroEncoder
 */

#ifndef __CINT__
#include "TFile.h"
#include "TDatime.h"
#include "TRandom.h"
#include "TArrayI.h"
#include "TArrayC.h"
#include "TSystem.h"
#include "AliRawDataHeader.h"
#include "AliAltroDecoder.h"
#include "AliAltroData.h"
#include "AliAltroBunch.h"
#include "AliHLTAltroEncoder.h"
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

// printouts or not
const bool bVerbose=false;

// the encoder can be used in two different modes: using
// AddSignal and terminating each channel by explicitely calling
// SetChannel, or the combined version AddChannelSignal
const bool bUseAddChannelSignal=true;

// some defaults
const int maxChannels=100;
const int maxBunches=50;
const int maxTimebin=1024;
const int maxSignal=1024;
const int maxEncodedDataSize=4000000+sizeofAliRawDataHeader;

// file dumps
const char* encDataDumpFile=NULL;//"/tmp/altro-enc.dat";
const char* encodedDataFile="/tmp/altro-encoded.dat";
const char* simulatedDataFile="/tmp/altro-simulated.dat";

/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////

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

/**
 * Print the bunch information of a channel.
 * The function goes recursively to the first bunch of the
 * channel and prints this one first.
 * @return number of 10bit words
 */
int DumpBunchRecursive(int nofBunches, int* pData, int size)
{
  if (!pData || size<2) return -1;
  int pos=size;
  int bunchLength=pData[--pos];
  int bunchEndTime=pData[--pos];
  int read=0;
  if (bunchLength<pos && nofBunches>1) {
    read=DumpBunchRecursive(nofBunches-1, pData, pos-bunchLength);
  }
  cout << "       bunch length " << bunchLength << ", end time " << bunchEndTime << ":";

  for (int time=bunchEndTime; time>bunchEndTime-bunchLength; time--) {	
    if (pos<1) {
      pos=0;
      cerr << "no data for bunch " << endl;
      return -1;
    }
    //cout << " " << pData[--pos];
  }
  cout << endl;
  if (read<0) return read;
  return read+2+bunchLength;
}

/**
 * Write the encoded data to a file.
 */
template<typename T> 
void DumpDataArray(T *array, const char* filename, int memberSize)
{
  if (array->GetSize()==0) return;

  ios::openmode filemode=(ios::openmode)0;
  ofstream rawfile(filename, filemode);
  if (rawfile.good()) {
    rawfile.write(reinterpret_cast<const char*>(array->GetArray()), array->GetSize()*memberSize);
  }
  rawfile.close();
}

/**
 * Print the simulated data.
 */
void PrintSimulatedData(TArrayI& data)
{
  cout << endl << "******  readback *******" << endl << endl;
  int dataPos=data.GetSize();
  while (dataPos>0) {
    int channelAddress=data[--dataPos];
    int nofBunches=data[--dataPos];
    cout << " ------------------ new channel -------------------------" << endl;
    if (nofBunches>0) {
      int read=DumpBunchRecursive(nofBunches, data.GetArray(), dataPos);
      if (read<0) return;
      dataPos-=read;
    }
    cout << " channel " << channelAddress << ":  number of bunches " << nofBunches << endl;
  }
}

int Compare(TArrayI& simData, TArrayC& encData)
{
  if (bVerbose) cout << endl << "******  compare *******" << endl << endl;

  // can not have a static AltroDecoder, crash when function
  // is called. I had a similar problem in the AliHLTAltroChannelSelectorComponent

  int iResult=0;
  AliAltroDecoder* decoder=new AliAltroDecoder;
  if (iResult>=0 && decoder->SetMemory((UChar_t*)encData.GetArray(), (UInt_t)encData.GetSize())<0) {
    cout << "error setting up decoder " << endl;
    iResult=-1;
  }

  if (iResult>=0 && !decoder->Decode()) {
    cout << "error decoding data" << endl;
    iResult=-1;    
  }

  AliAltroData altrochannel;
  int dataPos=simData.GetSize();
  while (iResult>=0 && decoder->NextChannel(&altrochannel)) {
    int hwadd=altrochannel.GetHadd();
    if (dataPos<2) {
      cout << "missmatch in simulated data" << endl;
      iResult=-1;
      break;
    }
    int channelAddress=simData[--dataPos];
    int nofBunches=simData[--dataPos];
    int nofBunchesBackup=nofBunches;
    if (channelAddress!=hwadd) {
      cout << "channel address missmatch: simulated " << channelAddress << " encoded " << hwadd << endl;
      iResult=-1;
      break;
    }

    if (bVerbose) cout << "comparing channel " << channelAddress << ", " << nofBunches << " bunche(s)" << endl;

    AliAltroBunch altrobunch;
    while (iResult>=0 && altrochannel.NextBunch(&altrobunch) && nofBunches-->0) {
      if (dataPos<2) {
	cout << "error reading bunch size and time bin from simulated data" << endl;
	iResult=-1;
	break;
      }
      int bunchLength=simData[--dataPos];
      int bunchEndTime=simData[--dataPos];
      if (bunchLength!=(int)altrobunch.GetBunchSize()) {
	cout << "bunch length missmatch: simulated " << bunchLength << " encoded " << altrobunch.GetBunchSize() << hex << " (" << altrobunch.GetBunchSize() << ")" << dec << endl;
	iResult=-1;
	break;
      }
      if (bunchEndTime!=(int)altrobunch.GetEndTimeBin()) {
	cout << "bunch end time missmatch: simulated " << bunchEndTime << " encoded " << altrobunch.GetEndTimeBin() << endl;
	iResult=-1;
	break;
      }
      if (bVerbose) cout << " bunch length " << bunchLength << ", end time " << bunchEndTime << endl;
      const  UInt_t* bunchData=altrobunch.GetData();
      if (bunchLength>dataPos) {
	cout << "error reading simulated bunch data, required "<< bunchLength << "  available " << dataPos << endl;
	iResult=-1;
	break;
      }
      Int_t* pSim=simData.GetArray();
      for (int bin=0; bin<bunchLength; bin++) {
	if ((int)bunchData[bin]!=pSim[dataPos-bunchLength+bin]) {
	  cout << "data missmatch at bunch position " << bin << " : simulated " << simData[dataPos] << " encoded " << bunchData[bin] << endl;
	  iResult=-1;
	  break;
	}
      }
      dataPos-=bunchLength;
    }
    if (nofBunches>0 && iResult==0) {
      cout << "error getting " << nofBunches << " of " << nofBunchesBackup << " bunches from channel " << channelAddress << endl;
      iResult=-1;
    }

  }

  delete decoder;
  return iResult;
}

void CompareDumpFiles()
{
  TString param=encodedDataFile;
  param+="?filetype=raw";
  TFile encFile(param);
  if (encFile.IsZombie()) {
    cout << "can not open file " << encodedDataFile << endl;
    return;
  }

  TArrayC encData(encFile.GetSize());
  if (encFile.ReadBuffer(encData.GetArray(), encData.GetSize())) {
    cout << "error reading file " << encodedDataFile << endl;
    return;
  }  

  param=simulatedDataFile;
  param+="?filetype=raw";
  TFile simFile(param);
  if (simFile.IsZombie()) {
    cout << "can not open file " << simulatedDataFile << endl;
    return;
  }

  TArrayI simData(simFile.GetSize()/4);
  if (simFile.ReadBuffer((char*)simData.GetArray(), simData.GetSize()*4)) {
    cout << "error reading file " << simulatedDataFile << endl;
    return;
  }

  Compare(simData, encData);
}

int testAliHLTAltroEncoder()
{
#ifdef __CINT__
  string ld_library_path="../.libs:";
  ld_library_path+=gSystem->GetDynamicPath();
  gSystem->SetDynamicPath(ld_library_path.c_str());
  gSystem->Load("libAliHLTRCU.so");
#endif
  AliHLTSystem gHLT;

  int nofChannels=GetRandom(1, maxChannels);
  if (nofChannels==0) nofChannels=1;
  TArrayI simData;
  TArrayC encData(maxEncodedDataSize);
  const int maxAltroDataSize=encData.GetSize()-sizeofAliRawDataHeader;
  int dataPos=0;

  AliHLTAltroEncoder encoder;
  encoder.SetBuffer((AliHLTUInt8_t*)encData.GetArray(), encData.GetSize());

  // set the common data header
  TArrayC dummyCdh(sizeofAliRawDataHeader);
  encoder.SetCDH((AliHLTUInt8_t*)dummyCdh.GetArray(), dummyCdh.GetSize());

  // set a trailer like in the real data format of the v1 RCU format (1 trailer word)
  Int_t trailer=0;
  encoder.SetRCUTrailer((AliHLTUInt8_t*)&trailer, 4);

  if (bVerbose) cout << "number of channels: " << nofChannels << endl;
  int channelAddress=-1;
  int lastChannel=-1;
  int nof10BitWords=0;
  // The nof10BitWords value is aligned to 4. In addition, one
  // 40bit word is needed for the channel trailer
  for (int channel=0; channel<nofChannels && (((nof10BitWords/4)+2)*5)<maxAltroDataSize; channel++) {
    channelAddress=GetRandom(0, 0xfff);
    if (channelAddress==lastChannel) continue;
    int nofBunches=GetRandom(1, maxBunches);
    if (nofBunches==0) continue;
    int totalBunches=0;
    int bunchEndTime=0;
    int lastBunchEndTime=0;
    
    if (bVerbose) cout << " ------------------ new channel -------------------------" << endl;
    int bunch=0;

    // Matthias Oct 2008:
    // Data was simulated in the wrong order: bunches for the higher timebins
    // first. But real data is the other way round: bunches are written in the order
    // of ascending timebins. A new consistency check was added to AliAltroDecoder
    // (AliAltroBunch) in revision 29090. Now the channels are checked for overlapping
    // bunches, the time of a bunch must be smaller than time of the previous one
    // minus its length.
    // Data is now simulated in the right order in order to fullfil this check.
    for (bunch=0; bunch<nofBunches && bunchEndTime<maxTimebin-3; bunch++) {
      while ((bunchEndTime+=GetRandom(0, maxTimebin-bunchEndTime))-lastBunchEndTime<3);
      int bunchLength=GetRandom(0, bunchEndTime-lastBunchEndTime);
      if (bunchLength==0) continue;
      // check if there is enough space for all the signals, end time
      // and bunch length. The value is aligned to 4. In addition, one
      // 40bit word is needed for the channel trailer
      if (((((nof10BitWords+bunchLength+2)/4)+2)*5)>=maxAltroDataSize) break;
      totalBunches++;

      if (bVerbose) cout << "       bunch " << bunch << ", length " << bunchLength << ", end time " << bunchEndTime << ":";

      if (simData.GetSize()<dataPos+bunchLength+2) simData.Set(dataPos+bunchLength+2);
      for (int time=bunchEndTime-bunchLength+1; time<=bunchEndTime; time++) {	
	int signal=GetRandom(0, maxSignal);
	simData[dataPos++]=signal;
	if (bVerbose) cout << " " << signal;
	int added=1;
	if (bUseAddChannelSignal) {
	  added=encoder.AddChannelSignal(signal, time, channelAddress);
	} else {
	  if (encoder.AddSignal(signal, time)<0) {
	    cout << "AddSignal failed" << endl;
	  }
	}
	lastChannel=channelAddress;
	if (added>0) nof10BitWords+=added;
      }
      if (bVerbose) cout << endl;
      simData[dataPos++]=bunchEndTime;
      simData[dataPos++]=bunchLength;
      lastBunchEndTime=bunchEndTime;
      bunchEndTime+=bunchLength;
    }
    if (simData.GetSize()<dataPos+2) simData.Set(dataPos+2);
    if (totalBunches>0) {
      simData[dataPos++]=totalBunches;
      simData[dataPos++]=channelAddress;
    }

    if (!bUseAddChannelSignal && lastChannel>=0) {
      encoder.SetChannel(lastChannel);
    }

    if (bVerbose) cout << " channel " << channelAddress << ":  number of bunches " << totalBunches << endl;

  }

  int dataSize=encoder.SetLength();
  if (dataSize<0) {
    cerr << "error finalizing encoded buffer" << endl;
    return -1;
  }
  int nof40bitWords=encoder.GetTotal40bitWords();
  encData.Set(dataSize);

  if (bVerbose) cout << "simulated data array:" << simData.GetSize() << " , ALTRO block length: " << nof40bitWords << " ALTRO words -> encoded data: " << encData.GetSize() << endl;

  /////////////////////////////////////////////////////////
  // some debugging options

  if (encDataDumpFile) DumpDataArray(&encData, encDataDumpFile, 1);

  simData.Set(dataPos);
  if (bVerbose) PrintSimulatedData(simData);

  /////////////////////////////////////////////////////////
  // read back end compare the encoded data

  if (simData.GetSize()>0 && Compare(simData, encData)<0) {
    // dump files for further debugging
    DumpDataArray(&encData, encodedDataFile, 1);
    DumpDataArray(&simData, simulatedDataFile, 4);
    return -1;
  } else {
    if (bVerbose) cout << "check successfully completed: " << nof40bitWords << " ALTRO words" << endl;
  }

  return 0;
}

int main(int /*argc*/, const char** /*argv*/)
{
//   CompareDumpFiles();
//   return 0;

  int iResult=0;
  int iCount=1000;
  for (int i=0; i<iCount; i++) {
    if ((iResult=testAliHLTAltroEncoder())<0) {
      cout << "missmatch in block no " << i << endl;
      return iResult;
    }
  }
  cout << "checking: " << iCount << " encoding cycle(s) successfully tested" << endl;
  return 0;
}
