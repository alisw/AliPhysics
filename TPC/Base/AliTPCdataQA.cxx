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
  July 2011:

  Changes to accomodate updates of general DQM/QA changes to have per trigger
  histograms (for a given event specie).

  AliTPCdataQA has a new flag for only keeping DQM info event by
  event!
  The expert/DA functionality has been kept exactly the same!
  

  June 2010

  This update should solve two problems mainly:
  * The vs event histograms have been limited to a fixed size for the
  DQM. The 500k seemed to be a big size but is no longer so, so we
  need to dynamically expand the range. The non-trivial point is that
  we also have to do it for the copy owned by AliTPCQADataMakerRec.
  * The amoreGui now remembers the x-range of the first visualization so
  the trick of setting the relevant event range as the histogram is
  filled no longer works.

  The fix is a bit crude but avoids creating a new histogram. Instead
  the range is expanded (max events and events per bin is doubled) but
  the number of bins is kept constant! In this way we can change just
  double the max of the X-axis of the hist and rebin the data. The
  same can easily be done for the copy owned by AliTPCQADataMakerRec.

  CAUTION:
  If we change the number of bins we could crash the whole system
  because ROOT does not create space for extra bins! (but we do not do
  this). In that way it is a crude solution.
  The rebinning in the code only works for an even number of bins.

  In addition to the above a bug in the reading of the config file was
  also found and corrected. fAdcMax was set instead of fEventsPerBin.

  Finally cout was changes to AliInfo.

  February 2008

  The code has been heavily modified so that now the RAW data is
  "expanded" for each sector and stored in a big signal array. Then a
  simple version of the code in AliTPCclusterer is used to identify
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


//Root includes
#include <TH1F.h>
#include <TString.h>
#include <TMath.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TError.h>
#include <TMap.h>
#include <TProfile.h>
//AliRoot includes
#include "AliRawReader.h"
#include "AliRawReaderRoot.h"
#include "AliRawReaderDate.h"
#include "AliTPCRawStreamV3.h"
#include "AliTPCCalROC.h"
#include "AliTPCROC.h"
#include "AliMathBase.h"
#include "TTreeStream.h"

//date
#include "event.h"
#include "AliTPCCalPad.h"
#include "AliTPCPreprocessorOnline.h"

//header file
#include "AliTPCdataQA.h"
#include "AliLog.h"


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
  fHistOccupancyVsEvent(0),
  fHistNclustersVsEvent(0),
  fEventCounter(0),
  fIsAnalysed(kFALSE),
  fMaxEvents(500000),           // Max events for event histograms
  fEventsPerBin(1000),          // Events per bin for event histograms
  fSignalCounter(0),            // Signal counter
  fClusterCounter(0),           // Cluster counter
  fAllBins(0),
  fAllSigBins(0),
  fAllNSigBins(0),
  fRowsMax(0),
  fPadsMax(0),
  fTimeBinsMax(0),
  fIsDQM(kFALSE),
  fHistOccVsSector(0x0),
  fHistOcc2dVsSector(0x0),
  fHistQVsSector(0x0),
  fHistQmaxVsSector(0x0),
  fOccVec(0x0),
  fOccMaxVec(0x0),
  fOccVecFine(0x0),
  fOccMaxVecFine(0x0)
{
  //
  // default constructor
  //
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
  fHistOccupancyVsEvent(0),
  fHistNclustersVsEvent(0),
  fEventCounter(ped.GetEventCounter()),
  fIsAnalysed(ped.GetIsAnalysed()),
  fMaxEvents(ped.GetMaxEvents()),
  fEventsPerBin(ped.GetEventsPerBin()),
  fSignalCounter(ped.GetSignalCounter()),
  fClusterCounter(ped.GetClusterCounter()),
  fAllBins(0),
  fAllSigBins(0),
  fAllNSigBins(0),
  fRowsMax(0),
  fPadsMax(0),
  fTimeBinsMax(0),
  fIsDQM(ped.GetIsDQM()),
  fHistOccVsSector(0x0),
  fHistOcc2dVsSector(0x0),
  fHistQVsSector(0x0),
  fHistQmaxVsSector(0x0),
  fOccVec(0x0),
  fOccMaxVec(0x0),
  fOccVecFine(0x0),
  fOccMaxVecFine(0x0)
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
  if(ped.GetHistQVsTimeSideA()) {
    fHistQVsTimeSideA = new TProfile(*ped.GetHistQVsTimeSideA());
    fHistQVsTimeSideA->SetDirectory(0);
  }
  if(ped.GetHistQVsTimeSideC()) {
    fHistQVsTimeSideC = new TProfile(*ped.GetHistQVsTimeSideC());
    fHistQVsTimeSideC->SetDirectory(0);
  }
  if(ped.GetHistQMaxVsTimeSideA()) {
    fHistQMaxVsTimeSideA = new TProfile(*ped.GetHistQMaxVsTimeSideA());
    fHistQMaxVsTimeSideA->SetDirectory(0);
  }
  if(ped.GetHistQMaxVsTimeSideC()) {
    fHistQMaxVsTimeSideC = new TProfile(*ped.GetHistQMaxVsTimeSideC());
    fHistQMaxVsTimeSideC->SetDirectory(0);
  }
  if(ped.GetHistOccupancyVsEventConst()) {
    fHistOccupancyVsEvent  = new TH1F(*ped.GetHistOccupancyVsEventConst());
    fHistOccupancyVsEvent->SetDirectory(0);
  }
  if(ped.GetHistNclustersVsEventConst()) {
    fHistNclustersVsEvent  = new TH1F(*ped.GetHistNclustersVsEventConst());
    fHistNclustersVsEvent->SetDirectory(0);
  }
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
  fHistOccupancyVsEvent(0),
  fHistNclustersVsEvent(0),
  fEventCounter(0),
  fIsAnalysed(kFALSE),
  fMaxEvents(500000),
  fEventsPerBin(1000),
  fSignalCounter(0),
  fClusterCounter(0),
  fAllBins(0),
  fAllSigBins(0),
  fAllNSigBins(0),
  fRowsMax(0),
  fPadsMax(0),
  fTimeBinsMax(0),
  fIsDQM(kFALSE),
  fHistOccVsSector(0x0),
  fHistOcc2dVsSector(0x0),
  fHistQVsSector(0x0),
  fHistQmaxVsSector(0x0),
  fOccVec(0x0),
  fOccMaxVec(0x0),
  fOccVecFine(0x0),
  fOccMaxVecFine(0x0)
{
  //
  // default constructor
  //
  if (config->GetValue("FirstTimeBin")) fFirstTimeBin = ((TObjString*)config->GetValue("FirstTimeBin"))->GetString().Atoi();
  if (config->GetValue("LastTimeBin"))  fLastTimeBin = ((TObjString*)config->GetValue("LastTimeBin"))->GetString().Atoi();
  if (config->GetValue("AdcMin"))       fAdcMin = ((TObjString*)config->GetValue("AdcMin"))->GetString().Atoi();
  if (config->GetValue("AdcMax"))       fAdcMax = ((TObjString*)config->GetValue("AdcMax"))->GetString().Atoi();
  if (config->GetValue("MaxEvents"))    fMaxEvents = ((TObjString*)config->GetValue("MaxEvents"))->GetString().Atoi();
  if (config->GetValue("EventsPerBin")) fEventsPerBin = ((TObjString*)config->GetValue("EventsPerBin"))->GetString().Atoi();
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
  delete fHistOccupancyVsEvent;
  delete fHistNclustersVsEvent;

  // DQM
  delete fHistOccVsSector;
  delete fHistOcc2dVsSector;
  delete fHistQVsSector;
  delete fHistQmaxVsSector;
  delete fOccVec;
  delete fOccMaxVec;
  delete fOccVecFine;
  delete fOccMaxVecFine;
  
  for (Int_t iRow = 0; iRow < fRowsMax; iRow++) {
    delete [] fAllBins[iRow];
    delete [] fAllSigBins[iRow];
  }  
  delete [] fAllBins;
  delete [] fAllSigBins;
  delete [] fAllNSigBins;
}

//_____________________________________________________________________
TH1F* AliTPCdataQA::GetHistOccupancyVsEvent()
{
  //
  // Create Occupancy vs event histogram
  // (we create this histogram differently then the other histograms
  //  because this we want to be able to access and copy
  //  from AliTPCQAMakerRec before it normally would be created)
  //
  if(!fHistOccupancyVsEvent) {

    Int_t nBins = fMaxEvents/fEventsPerBin;
    fHistOccupancyVsEvent = new TH1F("hOccupancyVsEvent", "Occupancy vs event number (~time); Event number; Occupancy", nBins, 0, nBins*fEventsPerBin);
    fHistOccupancyVsEvent->SetDirectory(0);
  }
  
  return fHistOccupancyVsEvent;
}

//_____________________________________________________________________
TH1F* AliTPCdataQA::GetHistNclustersVsEvent()
{
  //
  // Create Nclusters vs event histogram
  // (we create this histogram differently then the other histograms
  //  because this we want to be able to access and copy
  //  from AliTPCQAMakerRec before it normally would be created)
  //
  if(!fHistNclustersVsEvent) {

    Int_t nBins = fMaxEvents/fEventsPerBin;
    fHistNclustersVsEvent = new TH1F("hNclustersVsEvent", "Nclusters vs event number (~time); Event number; Nclusters per event", nBins, 0, nBins*fEventsPerBin);
    fHistNclustersVsEvent->SetDirectory(0);
  }
  
  return fHistNclustersVsEvent;
}

//_____________________________________________________________________
void AliTPCdataQA::UpdateEventHistograms()
{
  // Update histograms that display occupancy and 
  // number of clusters as a function of number of 
  // events
  if (!fHistOccupancyVsEvent)
    GetHistOccupancyVsEvent();
  if (!fHistNclustersVsEvent)
    GetHistNclustersVsEvent();
  
  if(fEventCounter > fMaxEvents) {
    
    // we have to expand the histogram to handle the larger number of
    // events. The way it is done now is to double the range and the
    // number of events per bin (so the number of histogram bins stays
    // constant)
    fEventsPerBin *= 2;
    fMaxEvents *= 2;

    // Change histogram limits
    const Int_t nBins = fHistOccupancyVsEvent->GetXaxis()->GetNbins();
    fHistOccupancyVsEvent->GetXaxis()->Set(nBins, fHistOccupancyVsEvent->GetXaxis()->GetNbins(), fMaxEvents);
    fHistNclustersVsEvent->GetXaxis()->Set(nBins, fHistNclustersVsEvent->GetXaxis()->GetNbins(), fMaxEvents);

    // Rebin the histogram
    for(Int_t bin = 1; bin <= nBins; bin+=2) {

      Int_t newBin = TMath::Nint(Float_t(bin+1)/2.0);
      Float_t newContent = (fHistOccupancyVsEvent->GetBinContent(bin)+
			    fHistOccupancyVsEvent->GetBinContent(bin+1))/2.0;
      fHistOccupancyVsEvent->SetBinContent(newBin, newContent); 

      newContent = (fHistNclustersVsEvent->GetBinContent(bin)+
		    fHistNclustersVsEvent->GetBinContent(bin+1))/2.0;
      fHistNclustersVsEvent->SetBinContent(newBin, newContent); 
    }

    // Set the remaining bins to 0
    Int_t lastHalf = nBins/2 +1;
    for(Int_t bin = lastHalf; bin <= nBins; bin++) {

      fHistOccupancyVsEvent->SetBinContent(bin, 0); 
      fHistNclustersVsEvent->SetBinContent(bin, 0); 
    }

    // In this case we should nut update but wait untill the new
    // number of events per bin is reached!
    return;
  }

  const Int_t bin = TMath::Nint(Float_t(fEventCounter)/fEventsPerBin);

  Float_t averageOccupancy =
    Float_t(fSignalCounter)/fEventsPerBin/(fLastTimeBin - fFirstTimeBin +1.0)
    / 570132.0; // 570,132 is number of pads
  fHistOccupancyVsEvent->SetBinContent(bin, averageOccupancy);
  fSignalCounter = 0;
  
  Float_t averageNclusters =
    Float_t(fClusterCounter)/fEventsPerBin;
  fHistNclustersVsEvent->SetBinContent(bin, averageNclusters);
  fClusterCounter = 0;
}

//_____________________________________________________________________
Bool_t AliTPCdataQA::ProcessEvent(AliTPCRawStreamV3 *const rawStreamV3)
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
      Int_t iPatch  = rawStreamV3->GetPatchIndex(); //  current patch
      Int_t iBranch = rawStreamV3->GetBranch();    //  current branch
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
          nSignals += Update(iSector,iRow,iPad,startTbin--,signal, iPatch, iBranch);
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
Bool_t AliTPCdataQA::ProcessEvent(AliRawReader *const rawReader)
{
  //
  //  Event processing loop - AliRawReader
  //
  AliTPCRawStreamV3 rawStreamV3(rawReader,(AliAltroMapping**)fMapping);
  Bool_t res=ProcessEvent(&rawStreamV3);
  if(res) {
    fEventCounter++; // only increment event counter if there is TPC data

    if(fEventCounter%fEventsPerBin==0) 
      UpdateEventHistograms();
  }
  return res;
}

//_____________________________________________________________________
Bool_t AliTPCdataQA::ProcessEvent(eventHeaderStruct *const event)
{
  //
  //  process date event
  //

  AliRawReaderDate rawReader((void*)event);
  Bool_t result=ProcessEvent(&rawReader);
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
			   Float_t signal,
			   const Int_t iPatch,
			   const Int_t iBranch)
{
  //
  // Signal filling method
  //
  
  //
  // Define the calibration objects the first time Update is called
  // NB! This has to be done first even if the data is rejected by the time 
  // cut to make sure that the objects are available in Analyse
  //
  if(!fIsDQM) {

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
    if (!fHistQVsTimeSideA) {
      fHistQVsTimeSideA  = new TProfile("hQVsTimeSideA", "Q vs time (side A); Time [Timebin]; Q [ADC ch]", 100, 0, 1000);
      fHistQVsTimeSideA->SetDirectory(0);
    }
    if (!fHistQVsTimeSideC) {
      fHistQVsTimeSideC  = new TProfile("hQVsTimeSideC", "Q vs time (side C); Time [Timebin]; Q [ADC ch]", 100, 0, 1000);
      fHistQVsTimeSideC->SetDirectory(0);
  }
    if (!fHistQMaxVsTimeSideA) {
      fHistQMaxVsTimeSideA  = new TProfile("hQMaxVsTimeSideA", "Q_{MAX} vs time (side A); Time [Timebin]; Q_{MAX} [ADC ch]", 100, 0, 1000);
      fHistQMaxVsTimeSideA->SetDirectory(0);
    }
    if (!fHistQMaxVsTimeSideC) {
      fHistQMaxVsTimeSideC  = new TProfile("hQMaxVsTimeSideC", "Q_{MAX} vs time (side C); Time [Timebin]; Q_{MAX} [ADC ch]", 100, 0, 1000);
      fHistQMaxVsTimeSideC->SetDirectory(0);
    }
  } else { // DQM histograms and array
    
    if (!fHistOccVsSector) {
      fHistOccVsSector  = new TProfile("hOccVsSector", "Occupancy vs sector; Sector; Occupancy", 72, 0, 72);
      fHistOccVsSector->SetDirectory(0);

      fHistOcc2dVsSector  = new TProfile2D("hOcc2dVsSector", "Occupancy vs sector and patch; Sector; Patch", 72, 0, 36, 6, 0, 6);
      fHistOcc2dVsSector->SetDirectory(0);

      fHistQVsSector  = new TProfile("hQVsSector", "Q vs sector; Sector; Q [ADC ch]", 72, 0, 72);
      fHistQVsSector->SetDirectory(0);

      fHistQmaxVsSector  = new TProfile("hQmaxVsSector", "Qmax vs sector; Sector; Qmax [ADC ch]", 72, 0, 72);
      fHistQmaxVsSector->SetDirectory(0);

      fOccVec = new TArrayD(72);
      for(Int_t i = 0; i < 72; i++)
	fOccVec->GetArray()[i] = 0;

      fOccMaxVec = new TArrayD(72);
      const Double_t nTimeBins = fLastTimeBin - fFirstTimeBin +1;
      for(Int_t i = 0; i < 72; i++)
	
	if(i<36) // IROCs (5504 pads)
	  fOccMaxVec->GetArray()[i] = nTimeBins*5504;
	else     // OROCs (9984 pads)
	  fOccMaxVec->GetArray()[i] = nTimeBins*9984;

      // 12 branches for each full sector
      const Int_t nBranches = 36*12;
      fOccVecFine = new TArrayD(nBranches);
      for(Int_t i = 0; i < nBranches; i++)
	fOccVecFine->GetArray()[i] = 0;

      // Pads per patch same for all sectors
      Int_t nPads0[6] = {1152, 1536, 1152, 1280, 1280, 1280};
      Int_t nPads1[6] = {1152, 1664, 1152, 1280, 1280, 1280};

      fOccMaxVecFine = new TArrayD(nBranches);
      for(Int_t i = 0; i < nBranches; i++) {
	
	const Int_t fullSector = Int_t(i/12);
	Int_t branch = i - fullSector*12;
	R__ASSERT(branch>=0 && branch<12);
	
	const Int_t patch = Int_t(branch/2);
	branch -= patch*2;
	
	R__ASSERT(branch>=0 && branch<2);
	if(branch == 0)
	  fOccMaxVecFine->GetArray()[i] = nTimeBins*nPads0[patch];
	else     // OROCs (9984 pads)
	  fOccMaxVecFine->GetArray()[i] = nTimeBins*nPads1[patch];
      }
    }
  }
  // Make the arrays for expanding the data

  if (!fAllBins)
    MakeArrays();

  //
  // If Analyse has been previously called we need now to denormalize the data
  // as more data is coming
  //
  if(fIsAnalysed == kTRUE && !fIsDQM) {
    
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
  
  if(fIsDQM) {

    fOccVec->GetArray()[iSector] += 1.0;
    // To change before committing
    if(iPatch>=0 && iBranch>=0 && iPatch<=5 && iBranch <= 1)
      fOccVecFine->GetArray()[(iSector%36)*12+iPatch*2+iBranch] += 1.0;
  } else {
    // In fNoThreshold we fill all data to estimate the ZS volume
    Float_t count = fNoThreshold->GetCalROC(iSector)->GetValue(iRow, iPad);
    fNoThreshold->GetCalROC(iSector)->SetValue(iRow, iPad,count+1);
  }

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

  fSignalCounter++;
  
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
      
      if(!fIsDQM) {
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
      
      if(fIsDQM) {
	fHistQVsSector->Fill(iSector, qTot);
	fHistQmaxVsSector->Fill(iSector, qMax);
      } else {
	Float_t charge = fMeanCharge->GetCalROC(iSector)->GetValue(iRow, iPad);
	fMeanCharge->GetCalROC(iSector)->SetValue(iRow, iPad, charge + qTot);
	
	Float_t count = fNTimeBins->GetCalROC(iSector)->GetValue(iRow, iPad);
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
      }
    } // end loop over signals
  } // end loop over rows
  
  fClusterCounter += nLocalMaxima;
}

//_____________________________________________________________________
void AliTPCdataQA::Analyse()
{
  //
  //  Calculate calibration constants
  //
  
  AliInfo("Analyse called");

  if(fIsDQM == kTRUE) {
    
    AliInfo("DQM flas is set -> No 2d information to analyze");
    return;
  }

  if(fIsAnalysed == kTRUE) {
    
    AliInfo("No new data since Analyse was called last time");
    return;
  }

  if(fEventCounter==0) {
    
      AliInfo("EventCounter == 0, Cannot analyse");
    return;
  }
  
  Int_t nTimeBins = fLastTimeBin - fFirstTimeBin +1;
  AliInfo(Form("EventCounter: %d , TimeBins: %d", fEventCounter, nTimeBins));

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
void AliTPCdataQA::MakeTree(const char *fname) const {
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

//______________________________________________________________________________
Float_t AliTPCdataQA::GetQ(const Float_t* adcArray, const Int_t time, 
			   const Int_t pad, const Int_t maxTimeBins, 
			   Int_t& timeMin, Int_t& timeMax, 
			   Int_t& padMin,  Int_t& padMax) const
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
void AliTPCdataQA::Streamer(TBuffer &xRuub)
{
  // Automatic schema evolution was first used from revision 4
  // Code based on:
  // http://root.cern.ch/root/roottalk/roottalk02/3207.html

   UInt_t xRuus, xRuuc;
   if (xRuub.IsReading()) {
      Version_t xRuuv = xRuub.ReadVersion(&xRuus, &xRuuc);
      //we use the automatic algorithm for class version > 3
      if (xRuuv > 3) {
	AliTPCdataQA::Class()->ReadBuffer(xRuub, this, xRuuv, xRuus,
					  xRuuc);
	return;
      }
      TH1F::Streamer(xRuub);
      xRuub >> fFirstTimeBin;
      xRuub >> fLastTimeBin;
      xRuub >> fAdcMin;
      xRuub >> fAdcMax;
      xRuub >> fNLocalMaxima;
      xRuub >> fMaxCharge;
      xRuub >> fMeanCharge;
      xRuub >> fNoThreshold;
      xRuub >> fNTimeBins;
      xRuub >> fNPads;
      xRuub >> fTimePosition;
      xRuub >> fEventCounter;
      xRuub >> fIsAnalysed;
      xRuub.CheckByteCount(xRuus, xRuuc, AliTPCdataQA::IsA());
   } else {
     AliTPCdataQA::Class()->WriteBuffer(xRuub,this);
   }
}

//____________________________________________________________________________________________
void AliTPCdataQA::FillOccupancyProfile()
{
  // This has to be filled at the end of the loop over data
  if(!fIsDQM) 
    AliInfo("Method only meaningful for DQM");
  
  for(Int_t i = 0; i < 72; i++) {

    fOccVec->GetArray()[i] /= fOccMaxVec->GetArray()[i];
    fHistOccVsSector->Fill(i, fOccVec->GetArray()[i]);
  }

  const Int_t nBranches = 36*12;
  for(Int_t i = 0; i < nBranches; i++) {

    fOccVecFine->GetArray()[i] /= fOccMaxVecFine->GetArray()[i];

    const Int_t fullSector = Int_t(i/12);

    Int_t branch = i - fullSector*12;
    const Int_t patch = Int_t(branch/2);

    branch -= patch*2;

    fHistOcc2dVsSector->Fill(fullSector+0.5*branch+0.1, patch+0.5, fOccVecFine->GetArray()[i]);
  }
}

//____________________________________________________________________________________________
void AliTPCdataQA::ResetProfiles()
{
  if(!fIsDQM) 
    AliInfo("Method only meaningful for DQM");
  
  if(fHistQVsSector)
    fHistQVsSector->Reset();
  if(fHistQmaxVsSector)
    fHistQmaxVsSector->Reset();
  if(fHistOccVsSector)
    fHistOccVsSector->Reset();
  if(fHistOcc2dVsSector)
    fHistOcc2dVsSector->Reset();

  if(fOccVec)
    for(Int_t i = 0; i < 72; i++)
      fOccVec->GetArray()[i] = 0.0;
  if(fOccVecFine)
    for(Int_t i = 0; i < 36*12; i++)
      fOccVecFine->GetArray()[i] = 0.0;
}
