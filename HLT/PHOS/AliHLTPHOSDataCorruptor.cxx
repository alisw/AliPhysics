/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Author:  Per Thomas Hille  <perthi@fys.uio.no>                 *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          * 
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include "AliHLTPHOSDataCorruptor.h"
#include "AliHLTPHOSPulseGenerator.h"

#include "TRandom.h"

#include <iostream>

using namespace std;


AliHLTPHOSDataCorruptor::AliHLTPHOSDataCorruptor():fPulseGeneratorPtr(0), fRandomGeneratorPtr(0)
{
  //  cout << " AliHLTPHOSDataCorruptor::AliHLTPHOSDataCorruptor()  creating new datacorruptor" << endl;
  fPulseGeneratorPtr = new AliHLTPHOSPulseGenerator(100, 0, 300, 2, 10);
  fRandomGeneratorPtr = new TRandom();
}


AliHLTPHOSDataCorruptor::AliHLTPHOSDataCorruptor(const AliHLTPHOSDataCorruptor & ):fPulseGeneratorPtr(0), fRandomGeneratorPtr(0)
{

}


AliHLTPHOSDataCorruptor::~AliHLTPHOSDataCorruptor()
{
  delete fPulseGeneratorPtr;
}


/**
 * Takes as input good data from an altro and makes
 * bad data from it by flipping bits in individual samples, adding noise
 * add a double pulse.. etc
 * @param dataArray data that should be corrupted
 * @param N  the number of samples in the array
 */
void
AliHLTPHOSDataCorruptor::MakeCorruptedData(Double_t * /*dataArray*/, int /*N*/)
{
  

}


void
AliHLTPHOSDataCorruptor::MakeCorruptedDataTest(Double_t *dataArray, int N)
{
  //  double testPulse[300];
  //  dataArray[300];
  int* quantisized = new int[N];

  fPulseGeneratorPtr->SetSampleFreq(10);
  fPulseGeneratorPtr->SetAmplitude(100);
  fPulseGeneratorPtr->SetTZero(0);
  fPulseGeneratorPtr->MakePulse(dataArray, N);

  cout <<endl <<endl;
  cout << "AliHLTPHOSDataCorruptor::MakeCorruptedDataTest: printing data array before corruption" <<endl;

  for(int i=0; i< N; i++)
    {
      cout << dataArray[i] <<"\t";
      quantisized[i] = (int)(dataArray[i]);
    }
 
  cout << endl;
  cout << "AliHLTPHOSDataCorruptor::MakeCorruptedDataTest: printing data after quantization" <<endl;  

  for(int i=0; i< N; i++)
    {
      cout << quantisized[i] <<"\t";
    }
  cout << endl;

  int bit = fRandomGeneratorPtr->Integer(10); 

  cout << "AliHLTPHOSDataCorruptor::MakeCorruptedDataTest: printing data flipping sample 20, bit " << bit <<endl; 
  
  FlipBit(&quantisized[10], bit);

  cout << endl;

  for(int i=0; i< N; i++)
    {
      cout << quantisized[i] <<"\t";
    }
  cout << endl;
  cout << endl;
  delete [] quantisized;
}


/**
 *Flipping a single bit in a sample from 0 if it is one, and 1 if it
 *is zero. This emluates a situaion that can occur on the readout bus
 *if the rcu driver is not able to drive the bus because of to low
 *impedance.
 *@param sample The samle the should have the bit(s) flipped
 *@param n      The number of bits to flip, the highest number is 10 (most significant) . 
*/
void 
AliHLTPHOSDataCorruptor::FlipBit(int *sample, int n)
{

  int mask = 1 << n;
  cout << "n = "<< n <<" mask = " << mask << endl;
  cout << "before flip"<< *sample << endl;
  *sample = *sample ^ mask;
  cout <<"after flip" << *sample <<endl;

}
