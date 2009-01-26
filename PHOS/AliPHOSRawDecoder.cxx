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
#include "TMath.h"

// --- AliRoot header files ---
#include "AliPHOSRawDecoder.h"
#include "AliRawReader.h"
#include "AliPHOSCalibData.h"
#include "AliLog.h"

ClassImp(AliPHOSRawDecoder)

//-----------------------------------------------------------------------------
AliPHOSRawDecoder::AliPHOSRawDecoder():
  fRawReader(0),fCaloStream(0),fPedSubtract(kFALSE),fEnergy(-111),fTime(-111),fQuality(0.),fPedestalRMS(0.),
  fAmpOffset(0),fModule(-1),fColumn(-1),fRow(-1),fNewModule(-1),fNewColumn(-1),fNewRow(-1),fNewAmp(0),fNewTime(0), 
  fLowGainFlag(kFALSE),fNewLowGainFlag(kFALSE),fOverflow(kFALSE),fSamples(0),fTimes(0),fCalibData(0)
{
  //Default constructor.
}

//-----------------------------------------------------------------------------
AliPHOSRawDecoder::AliPHOSRawDecoder(AliRawReader* rawReader,  AliAltroMapping **mapping):
  fRawReader(0),fCaloStream(0),fPedSubtract(kFALSE),fEnergy(-111),fTime(-111),fQuality(0.),fPedestalRMS(0.),
  fAmpOffset(0),fModule(-1),fColumn(-1),fRow(-1),fNewModule(-1),fNewColumn(-1),fNewRow(-1),fNewAmp(0),fNewTime(0),
  fLowGainFlag(kFALSE),fNewLowGainFlag(kFALSE),fOverflow(kFALSE),fSamples(0),fTimes(0),fCalibData(0)
{
  //Construct a decoder object.
  //Is is user responsibility to provide next raw event 
  //using AliRawReader::NextEvent().

  fRawReader =  rawReader;
  fCaloStream = new AliCaloRawStream(rawReader,"PHOS",mapping);
  fSamples = new TArrayI(100);
  fTimes = new TArrayI(100);
}

//-----------------------------------------------------------------------------
AliPHOSRawDecoder::~AliPHOSRawDecoder()
{
  //Destructor.

  if(fCaloStream){ delete fCaloStream; fCaloStream=0;}
  if(fSamples){ delete fSamples; fSamples=0 ;}
  if(fTimes){ delete fTimes; fTimes=0 ;}
}

//-----------------------------------------------------------------------------
AliPHOSRawDecoder::AliPHOSRawDecoder(const AliPHOSRawDecoder &phosDecoder ):
  fRawReader(phosDecoder.fRawReader),fCaloStream(phosDecoder.fCaloStream),
  fPedSubtract(phosDecoder.fPedSubtract),
  fEnergy(phosDecoder.fEnergy),fTime(phosDecoder.fTime),fQuality(phosDecoder.fQuality),fPedestalRMS(phosDecoder.fPedestalRMS),
  fAmpOffset(phosDecoder.fAmpOffset),fModule(phosDecoder.fModule),fColumn(phosDecoder.fColumn),
  fRow(phosDecoder.fRow),fNewModule(phosDecoder.fNewModule),fNewColumn(phosDecoder.fNewColumn),
  fNewRow(phosDecoder.fNewRow),fNewAmp(phosDecoder.fNewAmp),fNewTime(phosDecoder.fNewTime),
  fLowGainFlag(phosDecoder.fLowGainFlag),fNewLowGainFlag(phosDecoder.fNewLowGainFlag),
  fOverflow(phosDecoder.fOverflow),fSamples(phosDecoder.fSamples),
  fTimes(phosDecoder.fTimes),fCalibData(phosDecoder.fCalibData) 
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
    fQuality = phosDecode.fQuality ;
    fPedestalRMS = phosDecode.fPedestalRMS ;
    fAmpOffset = phosDecode.fAmpOffset ;
    fModule = phosDecode.fModule;
    fColumn = phosDecode.fColumn;
    fRow = phosDecode.fRow;
    fNewModule = phosDecode.fNewModule;
    fNewColumn = phosDecode.fNewColumn;
    fNewRow = phosDecode.fNewRow;
    fNewAmp = phosDecode.fNewAmp ;
    fNewTime= phosDecode.fNewTime ;
    fLowGainFlag = phosDecode.fLowGainFlag;
    fNewLowGainFlag = phosDecode.fNewLowGainFlag;
    fOverflow = phosDecode.fOverflow ;
    
    if(fSamples) delete fSamples;
    fSamples = phosDecode.fSamples;

    if(fTimes) delete fTimes;
    fTimes = phosDecode.fTimes;
    fCalibData = phosDecode.fCalibData; 
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
  
  fEnergy = -111;
  Float_t pedMean = 0;
  Float_t pedRMS = 0;
  Int_t   nPed = 0;
  Float_t baseLine = 1.0;
  const Int_t kPreSamples = 10;
  
  while ( in->Next() ) { 
    // Evaluate previous sample
    if(in->IsNewHWAddress() && fEnergy!=-111){ //Do not return at first sample
       
       //First remember new sample
       fNewLowGainFlag = in->IsLowGain();
       fNewModule = in->GetModule()+1;
       fNewRow    = in->GetRow()   +1;
       fNewColumn = in->GetColumn()+1;
       fNewAmp = in->GetSignal() ;
       fNewTime=in->GetTime() ;                                                                                                                               
       // We take the energy as a maximum amplitude
       // and the pedestal from the 0th point (30 Aug 2006).
       // Time is not evaluated 
       // Take it as a first time bin multiplied by the sample tick time
       
       if(fPedSubtract) {
	 if (nPed > 0){
           fPedestalRMS=(pedRMS-pedMean*pedMean/nPed)/nPed ;
           if(fPedestalRMS > 0.) 
            fPedestalRMS = TMath::Sqrt(fPedestalRMS) ;
	   fEnergy -= (Double_t)(pedMean/nPed); // pedestal subtraction
         }
	 else
	   return kFALSE;
       }
       else{
         //take pedestals from DB
         Double_t pedestal = (Double_t) fAmpOffset ;
         if(fCalibData){
           Float_t truePed = fCalibData->GetADCpedestalEmc(fModule, fColumn, fRow) ;
           Int_t   altroSettings = fCalibData->GetAltroOffsetEmc(fModule, fColumn, fRow) ;
           pedestal += truePed - altroSettings ;
         }
         else{
           printf("AliPHOSRawDecoder::NextDigit() Can not read data from OCDB \n") ;
         }
         fEnergy-=pedestal ;
       }
       if (fEnergy < baseLine) fEnergy = 0;

       return kTRUE;
     }

     fLowGainFlag = in->IsLowGain();
     fTime = 1;
     fModule = in->GetModule()+1;
     fRow    = in->GetRow()   +1;
     fColumn = in->GetColumn()+1;

    if(fLowGainFlag==fNewLowGainFlag && fModule==fNewModule &&
       fRow==fNewRow && fColumn==fNewColumn ){
       if(fNewAmp>fEnergy)  fEnergy = (Double_t)fNewAmp ;
       fNewModule=-1 ;  //copyed, do not copy more
    } 

     //Calculate pedestal if necessary
     if(fPedSubtract && (in->GetTime() < kPreSamples)) {
       pedMean += in->GetSignal();
       pedRMS+=in->GetSignal()*in->GetSignal() ;
       nPed++;
     }
     if((Double_t)in->GetSignal() > fEnergy)
       fEnergy = (Double_t)in->GetSignal();
     
   } // in.Next()
   
   return kFALSE;
}
