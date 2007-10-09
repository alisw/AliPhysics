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
#include "TH2F.h"
#include "TFile.h"
#include "AliPHOSRawDecoder.h"

ClassImp(AliPHOSCalibHistoProducer)

//-----------------------------------------------------------------------------
AliPHOSCalibHistoProducer::AliPHOSCalibHistoProducer() : 
  fRawDecoder(0),fHistoFile(0),fUpdatingRate(100),fIsOldRCUFormat(kFALSE),
  fEvents(0),fNbins(100),fXlow(0.),fXup(1000.)
{
  // Constructor: initializes data members
  // Checks existence of histograms which might have been left
  // from the previous runs to continues their filling

  fHistoFile =  new TFile("calibHisto.root","update");
  
  for(Int_t module=0; module<5; module++) {
    for(Int_t column=0; column<56; column++) {
      for(Int_t row=0; row<64; row++) {
	char hname[128];
	sprintf(hname,"mod%dcol%drow%d",module,column,row);
	TH1F* hist = (TH1F*)fHistoFile->Get(hname);
	if(hist) 
	  fAmpHisto[module][column][row]=hist;
	else 
	  fAmpHisto[module][column][row] = 0;     
      }
    }
  }
}

//-----------------------------------------------------------------------------            
AliPHOSCalibHistoProducer::AliPHOSCalibHistoProducer(Int_t nbinsx, Double_t xlow, Double_t xup) :
  fRawDecoder(0),fHistoFile(0),fUpdatingRate(100),fIsOldRCUFormat(kFALSE),
  fEvents(0),fNbins(nbinsx),fXlow(xlow),fXup(xup)
{
  // Constructor: initializes data members.
  // Checks existence of histograms which might have been left
  // from the previous runs to continues their filling.
  // In addition sets number of bins, low and upper limits common for all histograms.

  fHistoFile =  new TFile("calibHisto.root","update");

  for(Int_t module=0; module<5; module++) {
    for(Int_t column=0; column<56; column++) {
      for(Int_t row=0; row<64; row++) {
        char hname[128];
        sprintf(hname,"mod%dcol%drow%d",module,column,row);
        TH1F* hist = (TH1F*)fHistoFile->Get(hname);
        if(hist)
          fAmpHisto[module][column][row]=hist;
        else
          fAmpHisto[module][column][row] = 0;
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

}

//-----------------------------------------------------------------------------
AliPHOSCalibHistoProducer::AliPHOSCalibHistoProducer(const AliPHOSCalibHistoProducer &histoproducer) :
  TObject(histoproducer),fRawDecoder(histoproducer.fRawDecoder),fHistoFile(histoproducer.fHistoFile),
  fUpdatingRate(histoproducer.fUpdatingRate),fIsOldRCUFormat(histoproducer.fIsOldRCUFormat),
  fEvents(histoproducer.fEvents),fNbins(histoproducer.fNbins),fXlow(histoproducer.fXlow),fXup(histoproducer.fXup)
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

    fRawDecoder = histoproducer.fRawDecoder;
    fHistoFile = histoproducer.fHistoFile;
    fUpdatingRate = histoproducer.fUpdatingRate;
    fIsOldRCUFormat = histoproducer.fIsOldRCUFormat;
    fEvents = histoproducer.fEvents;
    fEvents = histoproducer.fEvents;
    fNbins = histoproducer.fNbins;
    fXlow = histoproducer.fXlow;
    fXup = histoproducer.fXup;   
    
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

  if(!fRawDecoder) AliFatal("Raw decoder not set!");
  
  Double_t energy;
  Int_t mod,col,row;
  
  if(fIsOldRCUFormat)
    fRawDecoder->SetOldRCUFormat(kTRUE);

  while(fRawDecoder->NextDigit()) {
    
    if(fRawDecoder->IsLowGain()) continue; 

    energy = fRawDecoder->GetEnergy();
    
    mod = fRawDecoder->GetModule()-1;
    col = fRawDecoder->GetColumn()-1;
    row = fRawDecoder->GetRow()-1;
    
    if(fAmpHisto[mod][col][row]) {
      fAmpHisto[mod][col][row]->Fill(energy);
    }
    else {
      char hname[128];
      sprintf(hname,"mod%dcol%drow%d",mod,col,row);
      fAmpHisto[mod][col][row] = new TH1F(hname,hname,fNbins,fXlow,fXup);
      fAmpHisto[mod][col][row]->Fill(energy);
    }
  }
    // update histograms in local file every 100th event
    if(fEvents != 0 && fEvents%fUpdatingRate == 0) {
      AliInfo(Form("Updating histo file, event %d, run %d\n",
		   fEvents,fRawDecoder->GetRawReader()->GetRunNumber()));
      UpdateHistoFile();
    }
    
    //   UpdateHistoFile();
    //   AliInfo(Form("%d events of run %d processed.",iEvent,runNum));
  
  fEvents++;
  
}

//-----------------------------------------------------------------------------
void AliPHOSCalibHistoProducer::UpdateHistoFile()
{
  // Write histograms to file

  if(!fHistoFile) return;
  if(!fHistoFile->IsOpen()) return;

  TH1F* hist=0;
  char hname[128];
  char htitle[128];

  for(Int_t module=0; module<5; module++) {
    sprintf(hname,"hMeanE%d",module);
    sprintf(htitle,"Mean energies in module %d",module);
    TH2F hMeanE(hname,htitle,56,0.,56.,64,0.,64);

    for(Int_t column=0; column<56; column++) {
      for(Int_t row=0; row<64; row++) {
	hist = fAmpHisto[module][column][row]; 
	if(hist) hist->Write(hist->GetName(),TObject::kWriteDelete);
	if(hist) hMeanE.SetBinContent(column,row,hist->GetMean());
      }
    }
    hMeanE.Write(hMeanE.GetName(),TObject::kWriteDelete);
  }

}
