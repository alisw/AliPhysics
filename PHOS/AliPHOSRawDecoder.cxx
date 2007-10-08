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

/* $Id$ */

// This class decodes the stream of ALTRO samples to extract
// the PHOS "digits" of current event.
// 
// Typical use case:
//     AliRawReader* rf = new AliRawReaderDate("2006run2211.raw");
//     AliPHOSRawDecoder dc(rf);
//     while (rf->NextEvent()) {
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
#include "TArrayI.h"

// --- AliRoot header files ---
#include "AliPHOSRawDecoder.h"
#include "AliPHOSPulseGenerator.h"

ClassImp(AliPHOSRawDecoder)

//-----------------------------------------------------------------------------
AliPHOSRawDecoder::AliPHOSRawDecoder():
  fRawReader(0),fCaloStream(0),fPedSubtract(kFALSE),fEnergy(-111),fTime(-111),fModule(-1),fColumn(-1),fRow(-1),fLowGainFlag(kFALSE),fSamples(0),fPulseGenerator(0)
{
  //Default constructor.
}

//-----------------------------------------------------------------------------
AliPHOSRawDecoder::AliPHOSRawDecoder(AliRawReader* rawReader):
  fRawReader(0),fCaloStream(0),fPedSubtract(kFALSE),fEnergy(-111),fTime(-111),fModule(-1),fColumn(-1),fRow(-1),fLowGainFlag(kFALSE),fSamples(0),fPulseGenerator(0)
{
  //Construct a decoder object.
  //Is is user responsibility to provide next raw event 
  //using AliRawReader::NextEvent().

  fRawReader =  rawReader;
  fCaloStream = new AliCaloRawStream(rawReader,"PHOS");
  fCaloStream->SetOldRCUFormat(kFALSE);
  fSamples = new TArrayI(100);
  fPulseGenerator = new AliPHOSPulseGenerator();
}

//-----------------------------------------------------------------------------
AliPHOSRawDecoder::~AliPHOSRawDecoder()
{
  //Destructor.

  if(fCaloStream) delete fCaloStream;
  if(fSamples) delete fSamples;
  if(fPulseGenerator) delete fPulseGenerator;
}

//-----------------------------------------------------------------------------
AliPHOSRawDecoder::AliPHOSRawDecoder(const AliPHOSRawDecoder &phosDecoder ):
  fRawReader(phosDecoder.fRawReader),fCaloStream(phosDecoder.fCaloStream),
  fPedSubtract(phosDecoder.fPedSubtract),
  fEnergy(phosDecoder.fEnergy),fTime(phosDecoder.fTime),
  fModule(phosDecoder.fModule),fColumn(phosDecoder.fColumn),
  fRow(phosDecoder.fRow),fLowGainFlag(phosDecoder.fLowGainFlag),
  fSamples(phosDecoder.fSamples),
  fPulseGenerator(phosDecoder.fPulseGenerator)
{
  //Copy constructor.
}

//-----------------------------------------------------------------------------
AliPHOSRawDecoder& AliPHOSRawDecoder::operator = (const AliPHOSRawDecoder &phosDecode)
{
  //Assignment operator.

  if(this != &phosDecode) {
    fRawReader = phosDecode.fRawReader;

    if(fCaloStream) delete fCaloStream;
    fCaloStream = phosDecode.fCaloStream;

    fEnergy = phosDecode.fEnergy;
    fTime = phosDecode.fTime;
    fModule = phosDecode.fModule;
    fColumn = phosDecode.fColumn;
    fRow = phosDecode.fRow;
    fLowGainFlag = phosDecode.fLowGainFlag;
    
    if(fSamples) delete fSamples;
    fSamples = phosDecode.fSamples;

    if(fPulseGenerator) delete fPulseGenerator;
    fPulseGenerator = phosDecode.fPulseGenerator;
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
  
  Int_t    iBin     = 0;
  Int_t    mxSmps   = fSamples->GetSize();
  Int_t    tLength  = 0;
  fEnergy = -111;
  Float_t pedMean = 0;
  Int_t   nPed = 0;
  Float_t baseLine = 1.0;
  const Float_t nPreSamples = 10;
  
  fSamples->Reset();
   while ( in->Next() ) { 

     if(!tLength) {
       tLength = in->GetTimeLength();
       if(tLength>mxSmps) {
	 fSamples->Set(tLength);
       }
     }
     
     // Fit the full sample
     if(in->IsNewHWAddress() && iBin>0) {
       
       iBin=0;
       
       // Temporarily we take the energy as a maximum amplitude
       // and the pedestal from the 0th point (30 Aug 2006).
       // Time is not evaluated for the moment (12.01.2007). 
       // Take is as a first time bin multiplied by the sample tick time
       
       if(fPedSubtract) 
	 fEnergy -= (Double_t)(pedMean/nPed); // "pedestal subtraction"
       
       if(fLowGainFlag)
	 fEnergy *= fPulseGenerator->GetRawFormatHighLowGainFactor(); // *16 

       if (fEnergy < baseLine) fEnergy = 0;

       pedMean = 0;
       return kTRUE;
     }

     fLowGainFlag = in->IsLowGain();
     fTime = fPulseGenerator->GetRawFormatTimeTrigger() * in->GetTime();
     fModule = in->GetModule()+1;
     fRow    = in->GetRow()   +1;
     fColumn = in->GetColumn()+1;

     // Fill array with samples
     iBin++;                                                             
     if(tLength-iBin < nPreSamples) {
       pedMean += in->GetSignal();
       nPed++;
     }
     fSamples->AddAt(in->GetSignal(),tLength-iBin);
     if((Double_t)in->GetSignal() > fEnergy) fEnergy = (Double_t)in->GetSignal();
     
   } // in.Next()
   
   return kFALSE;
}
