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

///////////////////////////////////////////////////////////////////////////////
// Class AliPHOSCalibHistoProducer accumulating histograms
// with amplitudes per PHOS channel
// It is intended to run at DAQ computers (LDC, GDC, HLT or MOOD)
// and it fills the histograms with amplitudes per channel.
// Usage example see in PHOS/macros/Shuttle/AliPHOSCalibHistoProducer.C
//
// Author: Boris Polichtchouk, 4 October 2006
///////////////////////////////////////////////////////////////////////////////

#include "AliLog.h"
#include "AliPHOSCalibHistoProducer.h"
#include "TH1.h"
#include "TFile.h"
#include "AliRawReader.h"
#include "AliPHOSDigit.h"
#include "AliPHOSRawDecoder.h"
#include "AliPHOSRawDigiProducer.h"
#include "AliPHOSGeometry.h"
#include "TClonesArray.h"

ClassImp(AliPHOSCalibHistoProducer)

//-----------------------------------------------------------------------------
AliPHOSCalibHistoProducer::AliPHOSCalibHistoProducer() : 
  fRawReader(0),fHistoFile(0),fUpdatingRate(100),fIsOldRCUFormat(kFALSE),
  fDigits(0),fEvents(0)
{
  // Constructor: initializes data members
  // Checks existence of histograms which might have been left
  // from the previous runs to continues their filling

  fHistoFile =  new TFile("calibHisto.root","update");
  fDigits = new TClonesArray("AliPHOSDigit",100);

  for(Int_t module=0; module<5; module++) {
    for(Int_t column=0; column<56; column++) {
      for(Int_t row=0; row<64; row++) {
	char hname[128];
	sprintf(hname,"mod%dcol%drow%d",module,column,row);
	TH1F* hist = (TH1F*)fHistoFile->Get(hname);
	if(hist) 
	  fAmpHisto[module][column][row]=hist;
	else
	  fAmpHisto[module][column][row]=0;
      }
    }
  }
}

//-----------------------------------------------------------------------------
AliPHOSCalibHistoProducer::~AliPHOSCalibHistoProducer()
{
  // Destructor
  
  UpdateHistoFile();

  if(fHistoFile) delete fHistoFile;
  if (fDigits) { fDigits->Delete(); delete fDigits; }

}

//-----------------------------------------------------------------------------
AliPHOSCalibHistoProducer::AliPHOSCalibHistoProducer(const AliPHOSCalibHistoProducer &histoproducer) :
  TObject(histoproducer),fRawReader(histoproducer.fRawReader),fHistoFile(histoproducer.fHistoFile),
  fUpdatingRate(histoproducer.fUpdatingRate),fIsOldRCUFormat(histoproducer.fIsOldRCUFormat),
  fDigits(histoproducer.fDigits),fEvents(histoproducer.fEvents)
{
  //Copy constructor.

  for(Int_t module=0; module<5; module++) {
    for(Int_t column=0; column<56; column++) {
      for(Int_t row=0; row<64; row++) {
	char hname[128];
	sprintf(hname,"mod%dcol%drow%d",module,column,row);
	TH1F* hist = (TH1F*)histoproducer.fHistoFile->Get(hname);
	if(hist) 
	  fAmpHisto[module][column][row]= new TH1F(*hist);
	else
	  fAmpHisto[module][column][row]=0;
      }
    }
  }
}

//-----------------------------------------------------------------------------
AliPHOSCalibHistoProducer& AliPHOSCalibHistoProducer::operator= 
(const AliPHOSCalibHistoProducer &histoproducer)
{
  //Assignment operator.

  if(this != &histoproducer) {

    fRawReader = histoproducer.fRawReader;
    fHistoFile = histoproducer.fHistoFile;
    fUpdatingRate = histoproducer.fUpdatingRate;
    fIsOldRCUFormat = histoproducer.fIsOldRCUFormat;
    fDigits = histoproducer.fDigits;
    fEvents = histoproducer.fEvents;

    for(Int_t module=0; module<5; module++) {
      for(Int_t column=0; column<56; column++) {
	for(Int_t row=0; row<64; row++) {
	  if(fAmpHisto[module][column][row]){
	    delete fAmpHisto[module][column][row];
	    fAmpHisto[module][column][row] = histoproducer.fAmpHisto[module][column][row];
	  }
	  else
	  fAmpHisto[module][column][row] = histoproducer.fAmpHisto[module][column][row];
	}
      }
    }


  }

  return *this;
}
//-----------------------------------------------------------------------------
void AliPHOSCalibHistoProducer::Run()
{
  // Reads raw data of current event and fills amplitude histograms
  // The histograms are written to file every fUpdatingRate events

  Int_t relId[4];
  AliPHOSDigit* digit;
  Double_t energy;
  
  AliPHOSGeometry* geo = AliPHOSGeometry::GetInstance();
  if(!geo) geo = AliPHOSGeometry::GetInstance("IHEP");
  if(!geo) AliFatal(Form("Cannot get PHOS geometry"));

  AliPHOSRawDecoder dc(fRawReader);
  if(fIsOldRCUFormat)
    dc.SetOldRCUFormat(kTRUE);
  
  AliPHOSRawDigiProducer producer;  
  producer.MakeDigits(fDigits,&dc);
    
//   printf("% digits.\n",fDigits->GetEntries());
  for(Int_t iDigit=0; iDigit<fDigits->GetEntries(); iDigit++) {
    digit = (AliPHOSDigit*)fDigits->At(iDigit);
    energy = digit->GetEnergy(); // no pedestal subtraction!
//     printf("[absId,E]: [%d,%.3f]\n",digit->GetId(),digit->GetEnergy());
    geo->AbsToRelNumbering(digit->GetId(),relId);	    
    Int_t mod = relId[0]-1;
    Int_t col = relId[3]-1;
    Int_t row = relId[2]-1;
//     printf("\t(mod,col,row): (%d,%d,%d)\n",mod,col,row);

    if(fAmpHisto[mod][col][row]) {
      fAmpHisto[mod][col][row]->Fill(energy);
    }
    else {
      char hname[128];
      sprintf(hname,"mod%dcol%drow%d",mod,col,row);
      fAmpHisto[mod][col][row] = new TH1F(hname,hname,100,0.,1000.);
    }
    
    // update histograms in local file every 100th event
    if(fEvents%fUpdatingRate == 0) {
      AliInfo(Form("Updating histo file, event %d, run %d\n",fEvents,fRawReader->GetRunNumber()));
      UpdateHistoFile();
    } 
  }
  
  fDigits->Delete();
  fEvents++;

//   }

//   UpdateHistoFile();
//   AliInfo(Form("%d events of run %d processed.",iEvent,runNum));
}

//-----------------------------------------------------------------------------
void AliPHOSCalibHistoProducer::UpdateHistoFile()
{
  // Write histograms to file

  if(!fHistoFile) return;
  if(!fHistoFile->IsOpen()) return;

  TH1F* hist=0;

  for(Int_t module=0; module<5; module++) {
    for(Int_t column=0; column<56; column++) {
      for(Int_t row=0; row<64; row++) {
	hist = fAmpHisto[module][column][row]; 
	if(hist) hist->Write(hist->GetName(),TObject::kWriteDelete);
      }
    }
  }

}
