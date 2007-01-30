/**************************************************************************
 * Copyright(c) 2007, ALICE Experiment at CERN, All rights reserved.      *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

// This class decodes the stream of ALTRO samples to extract
// the PHOS "digits" of current event.
// 
// Typical use case:
//     AliRawReader* rf = new AliRawReaderDate("2006run2211.raw");
//     while (rf->NextEvent()) {
//       AliPHOSRawDecoder dc(rf);
//       dc.SetOldRCUFormat(kTRUE);
//       dc.SubtractPedestals(kTRUE);
//       while ( dc.NextDigit() ) {
//         Int_t module = dc.GetModule();
//         Int_t column = dc.GetColumn();
//         Int_t row = dc.GetRow();
//         Double_t amplitude = dc.GetEnergy();
//         Double_t time = dc.GetTime();
//         Bool_t IsLowGain = dc.IsLowGain();
//            ..........
//       }
//     }

// Author: Boris Polichtchouk

// --- ROOT system ---
#include "TH1.h"

// --- AliRoot header files ---
#include "AliPHOSRawDecoder.h"
#include "AliPHOSPulseGenerator.h"

ClassImp(AliPHOSRawDecoder)

//-----------------------------------------------------------------------------
AliPHOSRawDecoder::AliPHOSRawDecoder():
  fRawReader(0),fCaloStream(0),fPedSubtract(kFALSE),fEnergy(-111),fTime(-111),fModule(-1),fColumn(-1),fRow(-1)
{
  //Default constructor.
}

//-----------------------------------------------------------------------------
AliPHOSRawDecoder::AliPHOSRawDecoder(AliRawReader* rawReader):
  fRawReader(0),fCaloStream(0),fPedSubtract(kFALSE),fEnergy(-111),fTime(-111),fModule(-1),fColumn(-1),fRow(-1)
{
  //Construct a decoder object.
  //Is is user responsibility to provide next raw event 
  //using AliRawReader::NextEvent().

  fRawReader =  rawReader;
  fCaloStream = new AliCaloRawStream(rawReader,"PHOS");
  fCaloStream->SetOldRCUFormat(kFALSE);
}

//-----------------------------------------------------------------------------
AliPHOSRawDecoder::~AliPHOSRawDecoder()
{
  //Destructor.

  if(fCaloStream) delete fCaloStream;
}

//-----------------------------------------------------------------------------
AliPHOSRawDecoder::AliPHOSRawDecoder(const AliPHOSRawDecoder &phosDecoder ):
  fRawReader(phosDecoder.fRawReader),fCaloStream(phosDecoder.fCaloStream),
  fPedSubtract(phosDecoder.fPedSubtract),
  fEnergy(phosDecoder.fEnergy),fTime(phosDecoder.fTime),
  fModule(phosDecoder.fModule),fColumn(phosDecoder.fColumn),
  fRow(phosDecoder.fRow)
{
  //Copy constructor.
}

//-----------------------------------------------------------------------------
AliPHOSRawDecoder& AliPHOSRawDecoder::operator = (const AliPHOSRawDecoder &phosDecode)
{
  if(this != &phosDecode) {
    fRawReader = phosDecode.fRawReader;

    if(fCaloStream) delete fCaloStream;
    fCaloStream = phosDecode.fCaloStream;

    fEnergy = phosDecode.fEnergy;
    fTime = phosDecode.fTime;
    fModule = phosDecode.fModule;
    fColumn = phosDecode.fColumn;
    fRow = phosDecode.fRow;
  }

  return *this;
}

//-----------------------------------------------------------------------------

Bool_t AliPHOSRawDecoder::NextDigit()
{
  //Extract an energy deposited in the crystal,
  //crystal' position (module,column,row),
  //time and gain (high or low).

  AliCaloRawStream* in = fCaloStream;
 
  // Create a shaper pulse object
  AliPHOSPulseGenerator pulse; 

  Bool_t   lowGainFlag = kFALSE ; 
  Int_t    iBin     = 0;

  // Create histogram to store samples
  TH1F hSamples("hSamples","ALTRO samples",in->GetTimeLength(),0,in->GetTimeLength());
  
  while ( in->Next() ) { 

    lowGainFlag = in->IsLowGain();
    // Fill histograms with samples
    hSamples.SetBinContent(in->GetTimeLength()-iBin-1,in->GetSignal());
    iBin++;

    // Fit the full sample
    if(iBin==in->GetTimeLength()) {
      iBin=0;

      // Temporarily we do not fit the sample graph, but
      // take the energy from the graph maximum, and the pedestal 
      // from the 0th point (30 Aug 2006).
      // Time is not evaluated for the moment (12.01.2007). 
      // Take is as a first time bin multiplied by the sample tick time

      fTime = pulse.GetRawFormatTimeTrigger() * in->GetTime();

      fModule = in->GetModule()+1;
      fRow = in->GetRow()   +1;
      fColumn = in->GetColumn()+1;

      fEnergy = hSamples.GetMaximum();       // "digit amplitude"
      if(fPedSubtract) 
	fEnergy-= hSamples.GetBinContent(0); // "pedestal subtraction"

      if(lowGainFlag)
	fEnergy *= pulse.GetRawFormatHighLowGainFactor(); // *16 
      
      return kTRUE;
    }

  } // in.Next()


  return kFALSE;
}
