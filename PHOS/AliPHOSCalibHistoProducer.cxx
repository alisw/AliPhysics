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
#include "AliPHOSRawStream.h"

ClassImp(AliPHOSCalibHistoProducer)

//-----------------------------------------------------------------------------
AliPHOSCalibHistoProducer::AliPHOSCalibHistoProducer(AliRawReader* rawReader) : 
  fRawReader(rawReader),fHistoFile(0),fUpdatingRate(100),fIsOldRCUFormat(kFALSE)
{
  // Constructor: initializes input data stream supplied by rawReader
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
	  fAmpHisto[module][column][row]=0;
      }
    }
  }
}

//-----------------------------------------------------------------------------
AliPHOSCalibHistoProducer::AliPHOSCalibHistoProducer() :
  fRawReader(0),fHistoFile(0),fUpdatingRate(0),fIsOldRCUFormat(kFALSE)
{
  // Default constructor
}

//-----------------------------------------------------------------------------
AliPHOSCalibHistoProducer::~AliPHOSCalibHistoProducer()
{
  // Destructor
  if(fHistoFile) delete fHistoFile;
}

//-----------------------------------------------------------------------------
void AliPHOSCalibHistoProducer::Run()
{
  // Reads raw data stream and fills amplitude histograms
  // The histograms are written to file every fUpdatingRate events

  TH1F* gHighGain = 0;
  TH1F* gLowGain = 0;
  Int_t iBin = 0;
  Int_t iEvent = 0;
  Int_t runNum = 0;

  AliPHOSRawStream in(fRawReader);
  if(fIsOldRCUFormat)
    in.SetOldRCUFormat(kTRUE);

  // Read raw data event by event

  while (fRawReader->NextEvent()) {
    runNum = fRawReader->GetRunNumber();

    while ( in.Next() ) { 
       
      if(!gHighGain) gHighGain = new TH1F("gHighGain","High gain",
					  in.GetTimeLength(),0,in.GetTimeLength());
      else
	if(gHighGain->GetNbinsX() != in.GetTimeLength()) {
	  delete gHighGain;
	  gHighGain = new TH1F("gHighGain","High gain",in.GetTimeLength(),0,in.GetTimeLength());
	}

      if(!gLowGain)  gLowGain = new TH1F("gLowGain","Low gain",
					 in.GetTimeLength(),0,in.GetTimeLength());
      else
	if(gLowGain->GetNbinsX() != in.GetTimeLength()) {
	  delete gLowGain;
	  gLowGain = new TH1F("gLowGain","Low gain",in.GetTimeLength(),0,in.GetTimeLength());
	}

      Bool_t lowGainFlag = in.IsLowGain();
      
      if(lowGainFlag) 
	gLowGain->SetBinContent(in.GetTimeLength()-iBin-1,in.GetSignal());
      else {
	gHighGain->SetBinContent(in.GetTimeLength()-iBin-1,in.GetSignal());
      }

      iBin++;

      if(iBin==in.GetTimeLength()) {
	iBin=0;

	Double_t energy;

	if(!lowGainFlag) {
	  energy = gHighGain->GetMaximum(); // no pedestal subtraction!
	}
	else {
	  energy = gLowGain->GetMaximum(); // no pedestal subtraction!
	}
	    
	Int_t mod = in.GetModule();
	Int_t col = in.GetColumn();
	Int_t row = in.GetRow();

	if(fAmpHisto[mod][col][row]) {
	  if(!lowGainFlag) {
	    fAmpHisto[mod][col][row]->Fill(energy);
	  }
	}
	else {
	  char hname[128];
	  sprintf(hname,"mod%dcol%drow%d",mod,col,row);
	  fAmpHisto[mod][col][row] = new TH1F(hname,hname,100,0.,1000.);
	}


      }
    }

    // update histograms in local file every 100th event
    if(iEvent%fUpdatingRate == 0) {
      AliInfo(Form("Updating histo file, event %d, run %d\n",iEvent,runNum));
      UpdateHistoFile();
    } 
    iEvent++;
  }

  UpdateHistoFile(); 
  AliInfo(Form("%d events of run %d processed.",iEvent,runNum));
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
