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
/* History of cvs commits:
 *
 * $Log$ 
 *
*/
///////////////////////////////////////////////////////////////////////////////
// Class AliEMCALCalibHistoProducer accumulating histograms
// with amplitudes per EMCAL channel
// It is intended to run at DAQ computers (LDC, GDC, HLT or MOOD)
// and it fills the histograms with amplitudes per channel.
// Usage example see in EMCAL/macros/Shuttle/AliEMCALCalibHistoProducer.C
//
// Author: Boris Polichtchouk, 4 October 2006
// Adapted for EMCAL by Gustavo Conesa Balbastre, October 2006
///////////////////////////////////////////////////////////////////////////////

#include "AliLog.h"
#include "AliEMCALCalibHistoProducer.h"
#include "TH1.h"
#include "TFile.h"
#include "TProfile.h"
#include "AliRawReader.h"
#include "AliCaloRawStream.h"

ClassImp(AliEMCALCalibHistoProducer)

//-----------------------------------------------------------------------------
AliEMCALCalibHistoProducer::AliEMCALCalibHistoProducer(AliRawReader* rawReader) : 
  fRawReader(rawReader),fHistoFile(0),fHistoFileName("calibHisto.root"),
  fUpdatingRate(100),fIsOldRCUFormat(kFALSE), fNSuperModules(12),  fNCellsEta (48),   
  fNCellsPhi(24),  fNCellsPhiHalfSM(12)
{
  // Constructor

  for(Int_t ism=0; ism<fNSuperModules; ism++) {
    fAmpProf[ism] = 0;
    fSMInstalled[ism]=kTRUE;
    for(Int_t icol=0; icol<fNCellsEta; icol++) 
      for(Int_t irow=0; irow<fNCellsPhi; irow++) 
	  fAmpHisto[ism][icol][irow]=0;
  }

}

//-----------------------------------------------------------------------------
AliEMCALCalibHistoProducer::AliEMCALCalibHistoProducer() :
  fRawReader(0),fHistoFile(0),fUpdatingRate(0),fIsOldRCUFormat(kFALSE),
  fNSuperModules(12), fNCellsEta (48), fNCellsPhi(24), fNCellsPhiHalfSM(12)
{
  // Default constructor
}

//-----------------------------------------------------------------------------
AliEMCALCalibHistoProducer::~AliEMCALCalibHistoProducer()
{
  // Destructor
  if(fHistoFile) {
    fHistoFile->Close();
    delete fHistoFile;
  }
}
//-----------------------------------------------------------------------------
void AliEMCALCalibHistoProducer::Init()
{
  // initializes input data stream supplied by rawReader
  // Checks existence of histograms which might have been left
  // from the previous runs to continue their filling
  fHistoFile =  new TFile(fHistoFileName,"update");
  char hname[128];
  Int_t nRow =  fNCellsPhi ;

  for(Int_t supermodule=0; supermodule<fNSuperModules; supermodule++) {
    //Check installed supermodules
    if(fSMInstalled[supermodule]==kFALSE) continue;
    //Check created profiles
    sprintf(hname,"mod%d",supermodule);
    TProfile* prof = (TProfile*)fHistoFile->Get(hname);
    if(prof)
      fAmpProf[supermodule]=prof;
    
    //Check created histograms
    if(supermodule > 10) nRow = fNCellsPhiHalfSM ; //Supermodules 11 and 12 are half supermodules
    for(Int_t column=0; column<fNCellsEta; column++) {
      for(Int_t row=0; row<nRow; row++) {
	sprintf(hname,"mod%dcol%drow%d",supermodule,column,row);
	TH1F* hist = (TH1F*)fHistoFile->Get(hname);
	if(hist) 
	  fAmpHisto[supermodule][column][row]=hist;
      }
    }
  }
  
}
//-----------------------------------------------------------------------------
void AliEMCALCalibHistoProducer::Run()
{
  // Reads raw data stream and fills amplitude histograms
  // The histograms are written to file every fUpdatingRate events
  //Also fills profiles to study the stability of supermodules during runs.

  Init();
  
//   TH1F* gHighGain = 0;
//   TH1F* gLowGain = 0;
  Int_t iBin = 0;
  Int_t iEvent = 0;
  Int_t runNum = 0;
  Int_t nProfFreq = 1000; //Number of events with which a bin of the TProfile if filled
  Int_t nEvtBins = 1000; //Total number of the profile survey bins.

  AliCaloRawStream in(fRawReader,"EMCAL");
  if(fIsOldRCUFormat)
    in.SetOldRCUFormat(kTRUE);

  // Read raw data event by event

  while (fRawReader->NextEvent()) {
    runNum = fRawReader->GetRunNumber();
    Float_t energy = 0;
     while ( in.Next() ) { 

      if(fSMInstalled[in.GetModule()]==kFALSE) continue;
       
      if (in.GetSignal() > energy) {
	energy = (Double_t) in.GetSignal();
      }    

      iBin++;

      if(iBin==in.GetTimeLength()) {
	iBin=0;
	    
	Int_t mod = in.GetModule();
	Int_t col = in.GetColumn();
	Int_t row = in.GetRow();
	Int_t evtbin = iEvent/nProfFreq;
	char hname[128];

	//Check if histogram/profile already exist, if not create it.
	if(!fAmpHisto[mod][col][row]) {
	  sprintf(hname,"mod%dcol%drow%d",mod,col,row);
	  fAmpHisto[mod][col][row] = new TH1F(hname,hname,1024,-0.5,1023.);
	}
	if(!fAmpProf[mod]) {
	  sprintf(hname,"mod%d",mod);
	  fAmpProf[mod] = new TProfile(hname,hname,nEvtBins,0.,nEvtBins);
	}
	//Fill histogram/profile 
	Bool_t lowGainFlag = in.IsLowGain();
	if(!lowGainFlag) {
	  fAmpHisto[mod][col][row]->Fill(energy);
	  fAmpProf[mod]->Fill(evtbin, energy);
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
void AliEMCALCalibHistoProducer::UpdateHistoFile()
{
  // Write histograms to file

  if(!fHistoFile) return;
  if(!fHistoFile->IsOpen()) return;

  TH1F* hist=0;
  TProfile* prof =0;
 
  Int_t nRow =  fNCellsPhi ;
  for(Int_t supermodule=0; supermodule<fNSuperModules; supermodule++) {
    
    prof = fAmpProf[supermodule]; 
    if(prof) prof->Write(prof->GetName(),TObject::kWriteDelete);
    
    if(supermodule > 10)  nRow = fNCellsPhiHalfSM ; //Supermodules 11 and 12 are half supermodules
    for(Int_t column=0; column<fNCellsEta; column++) {
      for(Int_t row=0; row<nRow; row++) {
	hist = fAmpHisto[supermodule][column][row]; 
	if(hist) hist->Write(hist->GetName(),TObject::kWriteDelete);
      }
    }
  }
  
}
