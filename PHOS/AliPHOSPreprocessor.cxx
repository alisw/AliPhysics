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
// PHOS Preprocessor class. It runs by Shuttle at the end of the run,
// calculates calibration coefficients and dead/bad channels
// to be posted in OCDB
//
// Author: Boris Polichtchouk, 4 October 2006
///////////////////////////////////////////////////////////////////////////////

#include "AliPHOSPreprocessor.h"
#include "AliLog.h"
#include "AliCDBMetaData.h"
#include "AliPHOSEmcCalibData.h"
#include "TFile.h"
#include "TH1.h"
#include "TMap.h"
#include "TRandom.h"
#include "TKey.h"
#include "TList.h"
#include "TObjString.h"
#include "AliPHOSEmcBadChannelsMap.h"

ClassImp(AliPHOSPreprocessor)

//_______________________________________________________________________________________
AliPHOSPreprocessor::AliPHOSPreprocessor() :
AliPreprocessor("PHS",0)
{
  //default constructor
}

//_______________________________________________________________________________________
AliPHOSPreprocessor::AliPHOSPreprocessor(AliShuttleInterface* shuttle):
AliPreprocessor("PHS",shuttle)
{
  // Constructor
}

//_______________________________________________________________________________________
UInt_t AliPHOSPreprocessor::Process(TMap* /*valueSet*/)
{
  // process data retrieved by the Shuttle

  // The fileName with the histograms which have been produced by
  // AliPHOSCalibHistoProducer.
  // It is a responsibility of the SHUTTLE framework to form the fileName
  
  TString runType = GetRunType();
  Log(Form("Run type: %s",runType.Data()));

  if(runType=="LED") {
    Bool_t ledOK = ProcessLEDRun();
    if(ledOK) return 0;
    else
      return 1;
  }

  gRandom->SetSeed(0); //the seed is set to the current  machine clock!
  AliPHOSEmcCalibData calibData;

  TList* list = GetFileSources(kDAQ, "AMPLITUDES");
  if(!list) {
    Log("Sources list not found, exit.");
    return 1;
  }

  TIter iter(list);
  TObjString *source;
  
  while ((source = dynamic_cast<TObjString *> (iter.Next()))) {
    AliInfo(Form("found source %s", source->String().Data()));

    TString fileName = GetFile(kDAQ, "AMPLITUDES", source->GetName());
    AliInfo(Form("Got filename: %s",fileName.Data()));

    TFile f(fileName);

    if(!f.IsOpen()) {
      Log(Form("File %s is not opened, something goes wrong!",fileName.Data()));
      return 1;
    }

    const Int_t nMod=5; // 1:5 modules
    const Int_t nCol=56; //1:56 columns in each module
    const Int_t nRow=64; //1:64 rows in each module

    Double_t coeff;
    char hnam[80];
    TH1F* histo=0;
    
    //Get the reference histogram
    //(author: Gustavo Conesa Balbastre)

    TList * keylist = f.GetListOfKeys();
    Int_t nkeys   = f.GetNkeys();
    Bool_t ok = kFALSE;
    TKey  *key;
    TString refHistoName= "";
    Int_t ikey = 0;
    Int_t counter = 0;
    TH1F* hRef = 0;
    
    //Check if the file contains any histogram
    
    if(nkeys< 2){
      Log(Form("Not enough histograms (%d) for calibration.",nkeys));
      return 1;
    }

    while(!ok){
      ikey = gRandom->Integer(nkeys);
      key = (TKey*)keylist->At(ikey);
      refHistoName = key->GetName();
      hRef = (TH1F*)f.Get(refHistoName);
      counter++;
      // Check if the reference histogram has too little statistics
      if(hRef->GetEntries()>2) ok=kTRUE;
      if(!ok && counter >= nMod*nCol*nRow){
	Log("No histogram with enough statistics for reference.");
	return 1;
      }
    }
    
    Log(Form("reference histogram %s, %.1f entries, mean=%.3f, rms=%.3f.",
	     hRef->GetName(),hRef->GetEntries(),
	     hRef->GetMean(),hRef->GetRMS()));

    Double_t refMean=hRef->GetMean();
    
    // Calculates relative calibration coefficients for all non-zero channels
    
    for(Int_t mod=0; mod<nMod; mod++) {
      for(Int_t col=0; col<nCol; col++) {
	for(Int_t row=0; row<nRow; row++) {
	  sprintf(hnam,"mod%dcol%drow%d",mod,col,row);
	  histo = (TH1F*)f.Get(hnam);
	  //TODO: dead channels exclusion!
	  if(histo) {
	    coeff = histo->GetMean()/refMean;
	    if(coeff>0)
	      calibData.SetADCchannelEmc(mod+1,col+1,row+1,1./coeff);
	    else 
	      calibData.SetADCchannelEmc(mod+1,col+1,row+1,1.);
	    AliInfo(Form("mod %d col %d row %d  coeff %f\n",mod,col,row,coeff));
	  }
	  else
	    calibData.SetADCchannelEmc(mod+1,col+1,row+1,1.); 
	}
      }
    }
    
    f.Close();
  }
  
  //Store EMC calibration data
  
  AliCDBMetaData emcMetaData;
  Bool_t emcOK = Store("Calib", "EmcGainPedestals", &calibData, &emcMetaData);

  if(emcOK) return 0;
  else
    return 1;

}


Bool_t AliPHOSPreprocessor::ProcessLEDRun()
{
  //Process LED run, fill bad channels map.

  AliPHOSEmcBadChannelsMap badMap;

  TList* list = GetFileSources(kDAQ, "LED");
  if(!list) {
    Log("Sources list for LED run not found, exit.");
    return kFALSE;
  }

  TIter iter(list);
  TObjString *source;
  char hnam[80];
  TH1F* histo=0;
  
  while ((source = dynamic_cast<TObjString *> (iter.Next()))) {

    AliInfo(Form("found source %s", source->String().Data()));

    TString fileName = GetFile(kDAQ, "LED", source->GetName());
    AliInfo(Form("Got filename: %s",fileName.Data()));

    TFile f(fileName);

    if(!f.IsOpen()) {
      Log(Form("File %s is not opened, something goes wrong!",fileName.Data()));
      return kFALSE;
    }

    const Int_t nMod=5; // 1:5 modules
    const Int_t nCol=56; //1:56 columns in each module
    const Int_t nRow=64; //1:64 rows in each module

    // Check for dead channels    
    Log(Form("Begin check for dead channels."));

    for(Int_t mod=0; mod<nMod; mod++) {
      for(Int_t col=0; col<nCol; col++) {
	for(Int_t row=0; row<nRow; row++) {
	  sprintf(hnam,"mod%dcol%drow%d",mod,col,row);
	  histo = (TH1F*)f.Get(hnam);
	  if(histo)
	    if (histo->GetMean()<1) {
	      Log(Form("Channel: [%d,%d,%d] seems dead, <E>=%.1f.",mod,col,row,histo->GetMean()));
	      badMap.SetBadChannel(mod,col,row);
	    }
	}
      }
    }

  }

  //Store bad channels map
  AliCDBMetaData badMapMetaData;

  //Bad channels data valid from current run fRun until updated (validityInfinite=kTRUE)
  Bool_t result = Store("Calib", "EmcBadChannels", &badMap, &badMapMetaData, fRun, kTRUE);
  return result;

}
