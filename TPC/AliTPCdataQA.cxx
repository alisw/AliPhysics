/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

/*
  February 2008

  The code has been heavily modified so that now the RAW data is
  "expanded" for each sector and stored in a big signal array. Then a
  simple version of the code in AliTPCclustererMI is used to identify
  the local maxima and these are then used for QA. This gives a better
  estimate of the charge (both max and total) and also limits the
  effect of noise.

  Implementation:

  In Update the RAW signals >= 3 ADC channels are stored in the arrays.
  
  There are 3 arrays:
  Float_t** fAllBins       2d array [row][bin(pad, time)] ADC signal
  Int_t**   fAllSigBins    2d array [row][signal#] bin(with signal)
  Int_t*    fAllNSigBins;  1d array [row] Nsignals

  This is done sector by sector.

  When data from a new sector is encountered, the method
  FindLocalMaxima is called on the data from the previous sector, and
  the calibration/data objects are updated with the "cluster"
  info. Finally the arrays are cleared.

  The requirements for a local maxima is:
  Charge in bin is >= 5 ADC channels.
  Charge in bin is larger than all the 8 neighboring bins.
  At least one of the two pad neighbors has a signal.
  At least one of the two time neighbors has a signal.

  Before accessing the data it is expected that the Analyse method is
  called. This normalizes some of the data objects to per event or per
  cluster. 
  If more data is passed to the class after Analyse has been called
  the normalization is reversed and Analyse has to be called again.
*/


// stl includes
#include <iostream>

using namespace std;

//Root includes
#include <TH1F.h>
#include <TH2F.h>
#include <TString.h>
#include <TMath.h>
#include <TF1.h>
#include <TRandom.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TError.h>
#include <TMap.h>
//AliRoot includes
#include "AliRawReader.h"
#include "AliRawReaderRoot.h"
#include "AliRawReaderDate.h"
#include "AliTPCRawStream.h"
#include "AliTPCRawStreamV3.h"
#include "AliTPCCalROC.h"
#include "AliTPCROC.h"
#include "AliMathBase.h"
#include "TTreeStream.h"
#include "AliTPCRawStreamFast.h"

//date
#include "event.h"
#include "AliTPCCalPad.h"
#include "AliTPCPreprocessorOnline.h"

//header file
#include "AliTPCdataQA.h"

ClassImp(AliTPCdataQA)

AliTPCdataQA::AliTPCdataQA() : /*FOLD00*/  
  fFirstTimeBin(60),
  fLastTimeBin(1000),
  fAdcMin(1),
  fAdcMax(100),
  fMapping(NULL),
  fPedestal(0),
  fNoise(0),
  fNLocalMaxima(0),
  fMaxCharge(0),
  fMeanCharge(0),
  fNoThreshold(0),
  fNTimeBins(0),
  fNPads(0),
  fTimePosition(0),
  fOverThreshold10(0),
  fOverThreshold20(0),
  fOverThreshold30(0),
  fHistQVsTimeSideA(0),
  fHistQVsTimeSideC(0),
  fHistQMaxVsTimeSideA(0),
  fHistQMaxVsTimeSideC(0),
  fEventCounter(0),
  fIsAnalysed(kFALSE),
  fAllBins(0),
  fAllSigBins(0),
  fAllNSigBins(0),
  fRowsMax(0),
  fPadsMax(0),
  fTimeBinsMax(0)
{
  //
  // default constructor
  //
}

//_____________________________________________________________________
AliTPCdataQA::AliTPCdataQA(AliRecoParam::EventSpecie_t es) :
fFirstTimeBin(60),
fLastTimeBin(1000),
fAdcMin(1),
fAdcMax(100),
fMapping(NULL),
fPedestal(0),
fNoise(0),
fNLocalMaxima(0),
fMaxCharge(0),
fMeanCharge(0),
fNoThreshold(0),
fNTimeBins(0),
fNPads(0),
fTimePosition(0),
fOverThreshold10(0),
fOverThreshold20(0),
fOverThreshold30(0),
fHistQVsTimeSideA(0),
fHistQVsTimeSideC(0),
fHistQMaxVsTimeSideA(0),
fHistQMaxVsTimeSideC(0),
fEventCounter(0),
fIsAnalysed(kFALSE),
fAllBins(0),
fAllSigBins(0),
fAllNSigBins(0),
fRowsMax(0),
fPadsMax(0),
fTimeBinsMax(0)
{
// ctor creating the histogram
  char *  name = Form("TPCRAW_%s", AliRecoParam::GetEventSpecieName(es)) ; 
  TH1F(name, name,100,0,100) ; 
}

//_____________________________________________________________________
AliTPCdataQA::AliTPCdataQA(const AliTPCdataQA &ped) : /*FOLD00*/
  TH1F(ped),
  fFirstTimeBin(ped.GetFirstTimeBin()),
  fLastTimeBin(ped.GetLastTimeBin()),
  fAdcMin(ped.GetAdcMin()),
  fAdcMax(ped.GetAdcMax()),
  fMapping(NULL),
  fPedestal(0),
  fNoise(0),
  fNLocalMaxima(0),
  fMaxCharge(0),
  fMeanCharge(0),
  fNoThreshold(0),
  fNTimeBins(0),
  fNPads(0),
  fTimePosition(0),
  fOverThreshold10(0),
  fOverThreshold20(0),
  fOverThreshold30(0),
  fHistQVsTimeSideA(0),
  fHistQVsTimeSideC(0),
  fHistQMaxVsTimeSideA(0),
  fHistQMaxVsTimeSideC(0),
  fEventCounter(ped.GetEventCounter()),
  fIsAnalysed(ped.GetIsAnalysed()),
  fAllBins(0),
  fAllSigBins(0),
  fAllNSigBins(0),
  fRowsMax(0),
  fPadsMax(0),
  fTimeBinsMax(0)
{
  //
  // copy constructor
  //
  if(ped.GetNLocalMaxima())
    fNLocalMaxima  = new AliTPCCalPad(*ped.GetNLocalMaxima());
  if(ped.GetMaxCharge())
    fMaxCharge      = new AliTPCCalPad(*ped.GetMaxCharge());
  if(ped.GetMeanCharge())
    fMeanCharge     = new AliTPCCalPad(*ped.GetMeanCharge());
  if(ped.GetNoThreshold())
    fNoThreshold  = new AliTPCCalPad(*ped.GetNoThreshold());
  if(ped.GetNTimeBins())
    fNTimeBins  = new AliTPCCalPad(*ped.GetNTimeBins());
  if(ped.GetNPads())
    fNPads  = new AliTPCCalPad(*ped.GetNPads());
  if(ped.GetTimePosition())
    fTimePosition  = new AliTPCCalPad(*ped.GetTimePosition());
  if(ped.GetOverThreshold10())
    fOverThreshold10  = new AliTPCCalPad(*ped.GetOverThreshold10());
  if(ped.GetOverThreshold20())
    fOverThreshold20  = new AliTPCCalPad(*ped.GetOverThreshold20());
  if(ped.GetOverThreshold30())
    fOverThreshold30  = new AliTPCCalPad(*ped.GetOverThreshold30());
  if(ped.GetHistQVsTimeSideA())
    fHistQVsTimeSideA = new TProfile(*ped.GetHistQVsTimeSideA());
  if(ped.GetHistQVsTimeSideC())
    fHistQVsTimeSideC = new TProfile(*ped.GetHistQVsTimeSideC());
  if(ped.GetHistQMaxVsTimeSideA())
    fHistQMaxVsTimeSideA = new TProfile(*ped.GetHistQMaxVsTimeSideA());
  if(ped.GetHistQMaxVsTimeSideC())
    fHistQMaxVsTimeSideC = new TProfile(*ped.GetHistQMaxVsTimeSideC());
}

//_____________________________________________________________________
AliTPCdataQA::AliTPCdataQA(const TMap *config) : /*FOLD00*/  
  TH1F("TPCRAW","TPCRAW",100,0,100),
  fFirstTimeBin(60),
  fLastTimeBin(1000),
  fAdcMin(1),
  fAdcMax(100),
  fMapping(NULL),
  fPedestal(0),
  fNoise(0),
  fNLocalMaxima(0),
  fMaxCharge(0),
  fMeanCharge(0),
  fNoThreshold(0),
  fNTimeBins(0),
  fNPads(0),
  fTimePosition(0),
  fOverThreshold10(0),
  fOverThreshold20(0),
  fOverThreshold30(0),
  fHistQVsTimeSideA(0),
  fHistQVsTimeSideC(0),
  fHistQMaxVsTimeSideA(0),
  fHistQMaxVsTimeSideC(0),
  fEventCounter(0),
  fIsAnalysed(kFALSE),
  fAllBins(0),
  fAllSigBins(0),
  fAllNSigBins(0),
  fRowsMax(0),
  fPadsMax(0),
  fTimeBinsMax(0)
{
  //
  // default constructor
  //
  if (config->GetValue("FirstTimeBin")) fFirstTimeBin = ((TObjString*)config->GetValue("FirstTimeBin"))->GetString().Atoi();
  if (config->GetValue("LastTimeBin"))  fLastTimeBin = ((TObjString*)config->GetValue("LastTimeBin"))->GetString().Atoi();
  if (config->GetValue("AdcMin"))       fAdcMin = ((TObjString*)config->GetValue("AdcMin"))->GetString().Atoi();
  if (config->GetValue("AdcMax"))       fAdcMax = ((TObjString*)config->GetValue("AdcMax"))->GetString().Atoi();
}

//_____________________________________________________________________
AliTPCdataQA& AliTPCdataQA::operator = (const  AliTPCdataQA &source)
{
  //
  // assignment operator
  //
  if (&source == this) return *this;
  new (this) AliTPCdataQA(source);

  return *this;
}


//_____________________________________________________________________
AliTPCdataQA::~AliTPCdataQA() /*FOLD00*/
{
  //
  // destructor
  //

  // do not delete fMapping, because we do not own it.
  // do not delete fMapping, because we do not own it.
  // do not delete fNoise and fPedestal, because we do not own them.

  delete fNLocalMaxima;
  delete fMaxCharge;
  delete fMeanCharge;
  delete fNoThreshold;
  delete fNTimeBins;
  delete fNPads;
  delete fTimePosition;
  delete fOverThreshold10;
  delete fOverThreshold20;
  delete fOverThreshold30;
  delete fHistQVsTimeSideA;
  delete fHistQVsTimeSideC;
  delete fHistQMaxVsTimeSideA;
  delete fHistQMaxVsTimeSideC;

  for (Int_t iRow = 0; iRow < fRowsMax; iRow++) {
    delete [] fAllBins[iRow];
    delete [] fAllSigBins[iRow];
  }  
  delete [] fAllBins;
  delete [] fAllSigBins;
  delete [] fAllNSigBins;
}
//_____________________________________________________________________
Bool_t AliTPCdataQA::ProcessEvent(AliTPCRawStreamV3 *rawStreamV3)
{
  //
  // Event Processing loop - AliTPCRawStreamV3
  //
  Bool_t withInput = kFALSE;
  Int_t nSignals = 0;
  Int_t lastSector = -1;
  
  while ( rawStreamV3->NextDDL() ){
    while ( rawStreamV3->NextChannel() ){
      Int_t iSector = rawStreamV3->GetSector(); //  current sector
      Int_t iRow    = rawStreamV3->GetRow();    //  current row
      Int_t iPad    = rawStreamV3->GetPad();    //  current pad
      if (iRow<0 || iPad<0) continue;
      // Call local maxima finder if the data is in a new sector
      if(iSector != lastSector) {
        
        if(nSignals>0)
          FindLocalMaxima(lastSector);
        
        CleanArrays();
        lastSector = iSector;
        nSignals = 0;
      }
      
      while ( rawStreamV3->NextBunch() ){
        Int_t  startTbin    = (Int_t)rawStreamV3->GetStartTimeBin();
        Int_t  bunchlength  = (Int_t)rawStreamV3->GetBunchLength();
        const UShort_t *sig = rawStreamV3->GetSignals();
        
        for (Int_t iTimeBin = 0; iTimeBin<bunchlength; iTimeBin++){
          Float_t signal=(Float_t)sig[iTimeBin];
          nSignals += Update(iSector,iRow,iPad,startTbin--,signal);
          withInput = kTRUE;
        }
      }
    }
  }
  
  if (lastSector>=0&&nSignals>0)
    FindLocalMaxima(lastSector);
  
  return withInput;
}

//_____________________________________________________________________
Bool_t AliTPCdataQA::ProcessEvent(AliRawReader *rawReader)
{
  //
  //  Event processing loop - AliRawReader
  //
  AliTPCRawStreamV3 *rawStreamV3 = new AliTPCRawStreamV3(rawReader, (AliAltroMapping**)fMapping);
  Bool_t res=ProcessEvent(rawStreamV3);
  delete rawStreamV3;
  if(res)
    fEventCounter++; // only increment event counter if there is TPC data
  return res;
}

//_____________________________________________________________________
Bool_t AliTPCdataQA::ProcessEventFast(AliTPCRawStreamFast *rawStreamFast)
{
  //
  // Event Processing loop - AliTPCRawStream
  //
  Bool_t withInput = kFALSE;
  Int_t nSignals = 0;
  Int_t lastSector = -1;

  while ( rawStreamFast->NextDDL() ){
    while ( rawStreamFast->NextChannel() ){
      
      Int_t iSector  = rawStreamFast->GetSector(); //  current sector
      Int_t iRow     = rawStreamFast->GetRow();    //  current row
      Int_t iPad     = rawStreamFast->GetPad();    //  current pad
  // Call local maxima finder if the data is in a new sector
      if(iSector != lastSector) {
        
        if(nSignals>0)
          FindLocalMaxima(lastSector);
        
        CleanArrays();
        lastSector = iSector;
        nSignals = 0;
      }
      
      while ( rawStreamFast->NextBunch() ){
        Int_t startTbin = (Int_t)rawStreamFast->GetStartTimeBin();
        Int_t endTbin = (Int_t)rawStreamFast->GetEndTimeBin();
        
        for (Int_t iTimeBin = startTbin; iTimeBin < endTbin; iTimeBin++){
          Float_t signal = rawStreamFast->GetSignals()[iTimeBin-startTbin];
          nSignals += Update(iSector,iRow,iPad,iTimeBin+1,signal);
          withInput = kTRUE;
        }
      }
    }
  }
  
  return withInput;
}
//_____________________________________________________________________
Bool_t AliTPCdataQA::ProcessEventFast(AliRawReader *rawReader)
{
  //
  //  Event processing loop - AliRawReader
  //
  AliTPCRawStreamFast *rawStreamFast = new AliTPCRawStreamFast(rawReader, (AliAltroMapping**)fMapping);
  Bool_t res=ProcessEventFast(rawStreamFast);
  if(res)
    fEventCounter++; // only increment event counter if there is TPC data
                     // otherwise Analyse (called in QA) fails

  delete rawStreamFast;
  return res;
}

//_____________________________________________________________________
Bool_t AliTPCdataQA::ProcessEvent(AliTPCRawStream *rawStream)
{
  //
  // Event Processing loop - AliTPCRawStream
  //

  Bool_t withInput = kFALSE;
  Int_t nSignals = 0;
  Int_t lastSector = -1;

  while (rawStream->Next()) {

    Int_t iSector  = rawStream->GetSector();      //  current ROC
    Int_t iRow     = rawStream->GetRow();         //  current row
    Int_t iPad     = rawStream->GetPad();         //  current pad
    Int_t iTimeBin = rawStream->GetTime();        //  current time bin
    Float_t signal = rawStream->GetSignal();      //  current ADC signal
    
    // Call local maxima finder if the data is in a new sector
    if(iSector != lastSector) {
      
      if(nSignals>0)
        FindLocalMaxima(lastSector);
      
      CleanArrays();
      lastSector = iSector;
      nSignals = 0;
    }
    
    // Sometimes iRow==-1 if there is problems to read the data
    if(iRow>=0) {
      nSignals += Update(iSector,iRow,iPad,iTimeBin,signal);
      withInput = kTRUE;
    }
  }

  if (lastSector>=0&&nSignals>0)
    FindLocalMaxima(lastSector);
  
  return withInput;
}


//_____________________________________________________________________
Bool_t AliTPCdataQA::ProcessEventOld(AliRawReader *rawReader)
{
  //
  //  Event processing loop - AliRawReader
  //

  // if fMapping is NULL the rawstream will crate its own mapping
  AliTPCRawStream rawStream(rawReader, (AliAltroMapping**)fMapping);
  rawReader->Select("TPC");
  Bool_t res =  ProcessEvent(&rawStream);

  if(res)
    fEventCounter++; // only increment event counter if there is TPC data
                     // otherwise Analyse (called in QA) fails
  return res;
}


//_____________________________________________________________________
Bool_t AliTPCdataQA::ProcessEvent(eventHeaderStruct *event)
{
  //
  //  process date event
  //

  AliRawReader *rawReader = new AliRawReaderDate((void*)event);
  Bool_t result=ProcessEvent(rawReader);
  delete rawReader;
  return result;
}



//_____________________________________________________________________
void AliTPCdataQA::DumpToFile(const Char_t *filename, const Char_t *dir, Bool_t append) /*FOLD00*/
{
  //
  //  Write class to file
  //

  TString sDir(dir);
  TString option;

  if ( append )
    option = "update";
  else
    option = "recreate";

  TDirectory *backup = gDirectory;
  TFile f(filename,option.Data());
  f.cd();
  if ( !sDir.IsNull() ){
    f.mkdir(sDir.Data());
    f.cd(sDir);
  }
  this->Write();
  f.Close();

  if ( backup ) backup->cd();
}


//_____________________________________________________________________
Int_t AliTPCdataQA::Update(const Int_t iSector, /*FOLD00*/
			   const Int_t iRow,
			   const Int_t iPad,
			   const Int_t iTimeBin,
			   Float_t signal)
{
  //
  // Signal filling method
  //
  
  //
  // Define the calibration objects the first time Update is called
  // NB! This has to be done first even if the data is rejected by the time 
  // cut to make sure that the objects are available in Analyse
  //
  if (!fNLocalMaxima) fNLocalMaxima = new AliTPCCalPad("NLocalMaxima","NLocalMaxima");
  if (!fMaxCharge) fMaxCharge = new AliTPCCalPad("MaxCharge","MaxCharge");
  if (!fMeanCharge) fMeanCharge = new AliTPCCalPad("MeanCharge","MeanCharge");
  if (!fNoThreshold) fNoThreshold = new AliTPCCalPad("NoThreshold","NoThreshold");
  if (!fNTimeBins) fNTimeBins = new AliTPCCalPad("NTimeBins","NTimeBins");
  if (!fNPads) fNPads = new AliTPCCalPad("NPads","NPads");
  if (!fTimePosition) fTimePosition = new AliTPCCalPad("TimePosition","TimePosition");
  if (!fOverThreshold10) fOverThreshold10 = new AliTPCCalPad("OverThreshold10","OverThreshold10");
  if (!fOverThreshold20) fOverThreshold20 = new AliTPCCalPad("OverThreshold20","OverThreshold20");
  if (!fOverThreshold30) fOverThreshold30 = new AliTPCCalPad("OverThreshold30","OverThreshold30");
  if (!fHistQVsTimeSideA)
    fHistQVsTimeSideA  = new TProfile("hQVsTimeSideA", "Q vs time (side A); Time [Timebin]; Q [ADC ch]", 100, 0, 1000);
  if (!fHistQVsTimeSideC)
    fHistQVsTimeSideC  = new TProfile("hQVsTimeSideC", "Q vs time (side C); Time [Timebin]; Q [ADC ch]", 100, 0, 1000);
  if (!fHistQMaxVsTimeSideA)
    fHistQMaxVsTimeSideA  = new TProfile("hQMaxVsTimeSideA", "Q_{MAX} vs time (side A); Time [Timebin]; Q_{MAX} [ADC ch]", 100, 0, 1000);
  if (!fHistQMaxVsTimeSideC)
    fHistQMaxVsTimeSideC  = new TProfile("hQMaxVsTimeSideC", "Q_{MAX} vs time (side C); Time [Timebin]; Q_{MAX} [ADC ch]", 100, 0, 1000);
  
  // Make the arrays for expanding the data

  if (!fAllBins)
    MakeArrays();

  //
  // If Analyse has been previously called we need now to denormalize the data
  // as more data is coming
  //
  if(fIsAnalysed == kTRUE) {
    
    const Int_t nTimeBins = fLastTimeBin - fFirstTimeBin +1;
    const Float_t denormalization = Float_t(fEventCounter * nTimeBins);
    fNoThreshold->Multiply(denormalization);  
    
    fMeanCharge->Multiply(fNLocalMaxima);
    fMaxCharge->Multiply(fNLocalMaxima);
    fNTimeBins->Multiply(fNLocalMaxima);
    fNPads->Multiply(fNLocalMaxima);
    fTimePosition->Multiply(fNLocalMaxima);
    fIsAnalysed = kFALSE;
  }

  //
  // TimeBin cut
  //
  if (iTimeBin<fFirstTimeBin) return 0;
  if (iTimeBin>fLastTimeBin) return 0;
  
  // if pedestal calibrations are loaded subtract pedestals
  if(fPedestal) {

    Float_t ped = fPedestal->GetCalROC(iSector)->GetValue(iRow, iPad);
    // Don't use data from pads where pedestals are abnormally small or large
    if(ped<10 || ped>90)
      return 0;
    signal -= ped;
  }
  
  // In fNoThreshold we fill all data to estimate the ZS volume
  Float_t count = fNoThreshold->GetCalROC(iSector)->GetValue(iRow, iPad);
  fNoThreshold->GetCalROC(iSector)->SetValue(iRow, iPad,count+1);
  
  // Require at least 3 ADC channels
  if (signal < 3.0)
    return 0;

  // if noise calibrations are loaded require at least 3*sigmaNoise
  if(fNoise) {
    
    Float_t noise = fNoise->GetCalROC(iSector)->GetValue(iRow, iPad);
    
    if(signal < noise*3.0)
      return 0;
  }

  //
  // This signal is ok and we store it in the cluster map
  //

  SetExpandDigit(iRow, iPad, iTimeBin, signal);
  
  return 1; // signal was accepted
}

//_____________________________________________________________________
void AliTPCdataQA::FindLocalMaxima(const Int_t iSector)
{
  //
  // This method is called after the data from each sector has been
  // exapanded into an array
  // Loop over the signals and identify local maxima and fill the
  // calibration objects with the information
  //

  Int_t nLocalMaxima = 0;
  const Int_t maxTimeBin = fTimeBinsMax+4; // Used to step between neighboring pads 
                                           // Because we have tha pad-time data in a 
                                           // 1d array

  for (Int_t iRow = 0; iRow < fRowsMax; iRow++) {

    Float_t* allBins = fAllBins[iRow];
    Int_t* sigBins   = fAllSigBins[iRow];
    const Int_t nSigBins   = fAllNSigBins[iRow];
    
    for (Int_t iSig = 0; iSig < nSigBins; iSig++) {

      Int_t bin  = sigBins[iSig];
      Float_t *b = &allBins[bin];

      //
      // Now we check if this is a local maximum
      //

      Float_t qMax = b[0];

      // First check that the charge is bigger than the threshold
      if (qMax<5) 
	continue;
      
      // Require at least one neighboring pad with signal
      if (b[-maxTimeBin]+b[maxTimeBin]<=0) continue;

      // Require at least one neighboring time bin with signal
      if (b[-1]+b[1]<=0) continue;
      
      //
      // Check that this is a local maximum
      // Note that the checking is done so that if 2 charges has the same
      // qMax then only 1 cluster is generated 
      // (that is why there is BOTH > and >=)
      //
      if (b[-maxTimeBin]   >= qMax) continue;
      if (b[-1  ]          >= qMax) continue; 
      if (b[+maxTimeBin]   > qMax)  continue; 
      if (b[+1  ]          > qMax)  continue; 
      if (b[-maxTimeBin-1] >= qMax) continue;
      if (b[+maxTimeBin-1] >= qMax) continue; 
      if (b[+maxTimeBin+1] > qMax)  continue; 
      if (b[-maxTimeBin+1] >= qMax) continue;
      
      //
      // Now we accept the local maximum and fill the calibration/data objects
      //
      nLocalMaxima++;

      Int_t iPad, iTimeBin;
      GetPadAndTimeBin(bin, iPad, iTimeBin);
      
      Float_t count = fNLocalMaxima->GetCalROC(iSector)->GetValue(iRow, iPad);
      fNLocalMaxima->GetCalROC(iSector)->SetValue(iRow, iPad, count+1);

      count = fTimePosition->GetCalROC(iSector)->GetValue(iRow, iPad);
      fTimePosition->GetCalROC(iSector)->SetValue(iRow, iPad, count+iTimeBin);
      
      Float_t charge = fMaxCharge->GetCalROC(iSector)->GetValue(iRow, iPad);
      fMaxCharge->GetCalROC(iSector)->SetValue(iRow, iPad, charge + qMax);
      
      if(qMax>=10) {
	count = fOverThreshold10->GetCalROC(iSector)->GetValue(iRow, iPad);
	fOverThreshold10->GetCalROC(iSector)->SetValue(iRow, iPad, count+1);
      }
      if(qMax>=20) {
	count = fOverThreshold20->GetCalROC(iSector)->GetValue(iRow, iPad);
	fOverThreshold20->GetCalROC(iSector)->SetValue(iRow, iPad, count+1);
      }
      if(qMax>=30) {
	count = fOverThreshold30->GetCalROC(iSector)->GetValue(iRow, iPad);
	fOverThreshold30->GetCalROC(iSector)->SetValue(iRow, iPad, count+1);
      }

      //
      // Calculate the total charge as the sum over the region:
      //
      //    o o o o o
      //    o i i i o
      //    o i C i o
      //    o i i i o
      //    o o o o o
      //
      // with qmax at the center C.
      //
      // The inner charge (i) we always add, but we only add the outer
      // charge (o) if the neighboring inner bin (i) has a signal.
      //
      Int_t minP = 0, maxP = 0, minT = 0, maxT = 0;
      Float_t qTot = qMax;
      for(Int_t i = -1; i<=1; i++) {
	for(Int_t j = -1; j<=1; j++) {
	  
	  if(i==0 && j==0)
	    continue;
	  
	  Float_t charge1 = GetQ(b, i, j, maxTimeBin, minT, maxT, minP, maxP);
	  qTot += charge1;
	  if(charge1>0) {
	    // see if the next neighbor is also above threshold
	    if(i*j==0) {
	      qTot += GetQ(b, 2*i, 2*j, maxTimeBin, minT, maxT, minP, maxP); 
	    } else {
	      // we are in a diagonal corner
	      qTot += GetQ(b,   i, 2*j, maxTimeBin, minT, maxT, minP, maxP); 
	      qTot += GetQ(b, 2*i,   j, maxTimeBin, minT, maxT, minP, maxP); 
	      qTot += GetQ(b, 2*i, 2*j, maxTimeBin, minT, maxT, minP, maxP); 
	    }
	  }
	}
      }
      
      charge = fMeanCharge->GetCalROC(iSector)->GetValue(iRow, iPad);
      fMeanCharge->GetCalROC(iSector)->SetValue(iRow, iPad, charge + qTot);
      
      count = fNTimeBins->GetCalROC(iSector)->GetValue(iRow, iPad);
      fNTimeBins->GetCalROC(iSector)->SetValue(iRow, iPad, count + maxT-minT+1);
      
      count = fNPads->GetCalROC(iSector)->GetValue(iRow, iPad);
      fNPads->GetCalROC(iSector)->SetValue(iRow, iPad, count + maxP-minP+1);
      
      if((iSector%36)<18) { // A side
	fHistQVsTimeSideA->Fill(iTimeBin, qTot);
	fHistQMaxVsTimeSideA->Fill(iTimeBin, qMax);
      } else {
	fHistQVsTimeSideC->Fill(iTimeBin, qTot);
	fHistQMaxVsTimeSideC->Fill(iTimeBin, qMax);      
      }
    } // end loop over signals
  } // end loop over rows
  
  //  cout << "Number of local maximas found: " << nLocalMaxima << endl;
}

//_____________________________________________________________________
void AliTPCdataQA::Analyse()
{
  //
  //  Calculate calibration constants
  //
  
  cout << "Analyse called" << endl;

  if(fIsAnalysed == kTRUE) {
    
    cout << "No new data since Analyse was called last time" << endl;
    return;
  }

  if(fEventCounter==0) {
    
    cout << "EventCounter == 0, Cannot analyse" << endl;
    return;
  }
  
  Int_t nTimeBins = fLastTimeBin - fFirstTimeBin +1;
  cout << "EventCounter: " << fEventCounter << endl
       << "TimeBins: " << nTimeBins << endl;

  Float_t normalization = 1.0 / Float_t(fEventCounter * nTimeBins);
  fNoThreshold->Multiply(normalization);  
  
  fMeanCharge->Divide(fNLocalMaxima);
  fMaxCharge->Divide(fNLocalMaxima);
  fNTimeBins->Divide(fNLocalMaxima);
  fNPads->Divide(fNLocalMaxima);
  fTimePosition->Divide(fNLocalMaxima);

  fIsAnalysed = kTRUE;
}


//_____________________________________________________________________
void AliTPCdataQA::MakeTree(const char *fname){
  //
  // Export result to the tree -located in the file
  // This file can be analyzed using AliTPCCalibViewer
  // 
  AliTPCPreprocessorOnline preprocesor;

  if (fNLocalMaxima) preprocesor.AddComponent(fNLocalMaxima);
  if (fMaxCharge) preprocesor.AddComponent(fMaxCharge);  
  if (fMeanCharge) preprocesor.AddComponent(fMeanCharge);  
  if (fNoThreshold) preprocesor.AddComponent(fNoThreshold);
  if (fNTimeBins) preprocesor.AddComponent(fNTimeBins);
  if (fNPads) preprocesor.AddComponent(fNPads);
  if (fTimePosition) preprocesor.AddComponent(fTimePosition);
  if (fOverThreshold10) preprocesor.AddComponent(fOverThreshold10);
  if (fOverThreshold20) preprocesor.AddComponent(fOverThreshold20);
  if (fOverThreshold30) preprocesor.AddComponent(fOverThreshold30);

  preprocesor.DumpToFile(fname);  
}


//_____________________________________________________________________
void AliTPCdataQA::MakeArrays(){
  //
  // The arrays for expanding the raw data are defined and 
  // som parameters are intialised
  //
  AliTPCROC * roc = AliTPCROC::Instance();
  //
  // To make the array big enough for all sectors we take 
  // the dimensions from the outer row of an OROC (the last sector)
  //
  fRowsMax     = roc->GetNRows(roc->GetNSector()-1);
  fPadsMax     = roc->GetNPads(roc->GetNSector()-1,fRowsMax-1);
  fTimeBinsMax = fLastTimeBin - fFirstTimeBin +1; 

  //
  // Since we have added 2 pads (TimeBins) before and after the real pads (TimeBins) 
  // to make sure that we can always query the exanded table even when the 
  // max is on the edge
  //

 
  fAllBins = new Float_t*[fRowsMax];
  fAllSigBins = new Int_t*[fRowsMax];
  fAllNSigBins = new Int_t[fRowsMax];

  for (Int_t iRow = 0; iRow < fRowsMax; iRow++) {
    //
    Int_t maxBin = (fTimeBinsMax+4)*(fPadsMax+4);  
    fAllBins[iRow] = new Float_t[maxBin];
    memset(fAllBins[iRow],0,sizeof(Float_t)*maxBin); // set all values to 0
    fAllSigBins[iRow] = new Int_t[maxBin];
    fAllNSigBins[iRow] = 0;
  }
}


//_____________________________________________________________________
void AliTPCdataQA::CleanArrays(){
  //
  //
  //

  for (Int_t iRow = 0; iRow < fRowsMax; iRow++) {

    // To speed up the performance by a factor 2 on cosmic data (and
    // presumably pp data as well) where the ocupancy is low, the
    // memset is only called if there is more than 1000 signals for a
    // row (of the order 1% occupancy)
    if(fAllNSigBins[iRow]<1000) {
      
      Float_t* allBins = fAllBins[iRow];
      Int_t* sigBins   = fAllSigBins[iRow];
      const Int_t nSignals = fAllNSigBins[iRow];
      for(Int_t i = 0; i < nSignals; i++)
	allBins[sigBins[i]]=0;      
    } else {

      Int_t maxBin = (fTimeBinsMax+4)*(fPadsMax+4); 
      memset(fAllBins[iRow],0,sizeof(Float_t)*maxBin);
    }

    fAllNSigBins[iRow]=0;
  }
}

//_____________________________________________________________________
void AliTPCdataQA::GetPadAndTimeBin(Int_t bin, Int_t& iPad, Int_t& iTimeBin){
  //
  // Return pad and timebin for a given bin
  //
  
  //  Int_t bin = iPad*(fTimeBinsMax+4)+iTimeBin;
  iTimeBin  = bin%(fTimeBinsMax+4);
  iPad      = (bin-iTimeBin)/(fTimeBinsMax+4);

  iPad     -= 2;
  iTimeBin -= 2;
  iTimeBin += fFirstTimeBin;

  R__ASSERT(iPad>=0 && iPad<=fPadsMax);
  R__ASSERT(iTimeBin>=fFirstTimeBin && iTimeBin<=fLastTimeBin);
}

//_____________________________________________________________________
void AliTPCdataQA::SetExpandDigit(const Int_t iRow, Int_t iPad, 
				  Int_t iTimeBin, const Float_t signal) 
{
  //
  // 
  //
  R__ASSERT(iRow>=0 && iRow<fRowsMax);
  R__ASSERT(iPad>=0 && iPad<=fPadsMax);
  R__ASSERT(iTimeBin>=fFirstTimeBin && iTimeBin<=fLastTimeBin);
  
  iTimeBin -= fFirstTimeBin;
  iPad     += 2;
  iTimeBin += 2;
  
  Int_t bin = iPad*(fTimeBinsMax+4)+iTimeBin;
  fAllBins[iRow][bin] = signal;
  fAllSigBins[iRow][fAllNSigBins[iRow]] = bin;
  fAllNSigBins[iRow]++;
}

Float_t AliTPCdataQA::GetQ(const Float_t* adcArray, const Int_t time, 
			   const Int_t pad, const Int_t maxTimeBins, 
			   Int_t& timeMin, Int_t& timeMax, 
			   Int_t& padMin,  Int_t& padMax) 
{
  //
  // This methods return the charge in the bin time+pad*maxTimeBins
  // If the charge is above 0 it also updates the padMin, padMax, timeMin
  // and timeMax if necessary
  //
  Float_t charge = adcArray[time + pad*maxTimeBins];
  if(charge > 0) {
    timeMin = TMath::Min(time, timeMin); timeMax = TMath::Max(time, timeMax);
    padMin = TMath::Min(pad, padMin); padMax = TMath::Max(pad, padMax);
  }
  return charge; 
}

//______________________________________________________________________________
void AliTPCdataQA::Streamer(TBuffer &R__b)
{
  // Automatic schema evolution was first used from revision 4
  // Code based on:
  // http://root.cern.ch/root/roottalk/roottalk02/3207.html

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c);
      //we use the automatic algorithm for class version > 3
      if (R__v > 3) {
	AliTPCdataQA::Class()->ReadBuffer(R__b, this, R__v, R__s,
					  R__c);
	return;
      }
      TH1F::Streamer(R__b);
      R__b >> fFirstTimeBin;
      R__b >> fLastTimeBin;
      R__b >> fAdcMin;
      R__b >> fAdcMax;
      R__b >> fNLocalMaxima;
      R__b >> fMaxCharge;
      R__b >> fMeanCharge;
      R__b >> fNoThreshold;
      R__b >> fNTimeBins;
      R__b >> fNPads;
      R__b >> fTimePosition;
      R__b >> fEventCounter;
      R__b >> fIsAnalysed;
      R__b.CheckByteCount(R__s, R__c, AliTPCdataQA::IsA());
   } else {
     AliTPCdataQA::Class()->WriteBuffer(R__b,this);
   }
}
