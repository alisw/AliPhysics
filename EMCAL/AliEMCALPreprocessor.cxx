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
 * Revision 1.2  2006/12/12 17:16:09  gustavo
 * Detector name hardcoded in Preprocesor with new detector name notation (3 letters). New way to take reference histogram to avoid problems in case of low number of entries or no existing histogram. Change return 0 by return 1
 *
 * Revision 1.1  2006/12/07 16:32:16  gustavo
 * First shuttle code, online calibration histograms producer, EMCAL preprocessor
 * 
 *
*/
///////////////////////////////////////////////////////////////////////////////
// EMCAL Preprocessor class. It runs by Shuttle at the end of the run,
// calculates calibration coefficients and dead/bad channels
// to be posted in OCDB
//
// Author: Boris Polichtchouk, 4 October 2006
// Adapted for EMCAL by Gustavo Conesa Balbastre, October 2006
///////////////////////////////////////////////////////////////////////////////

//Root
#include "TFile.h"
#include "TH1.h"
#include "TMap.h"
#include "TRandom.h"
#include "TKey.h"
#include "TList.h"
#include "TString.h"

//AliRoot
#include "AliEMCALPreprocessor.h"
#include "AliLog.h"
#include "AliCDBMetaData.h"
#include "AliEMCALCalibData.h"

ClassImp(AliEMCALPreprocessor)

//_______________________________________________________________________________________
AliEMCALPreprocessor::AliEMCALPreprocessor() :
AliPreprocessor("EMC",0)
{
  //default constructor
}

//_______________________________________________________________________________________
AliEMCALPreprocessor::AliEMCALPreprocessor(AliShuttleInterface* shuttle):
AliPreprocessor("EMC",shuttle)
{
  // Constructor
}

//_______________________________________________________________________________________
UInt_t AliEMCALPreprocessor::Process(TMap* /*valueSet*/)
{
  // process data retrieved by the Shuttle

  // The fileName with the histograms which have been produced by
  // AliEMCALCalibHistoProducer.
  // It is a responsibility of the SHUTTLE framework to form the fileName

  TString  fileName = GetFile(kDAQ, "AMPLITUDES", "GDC");
  Log(Form("Got filename: %s",fileName.Data()));

  TFile f(fileName);

  if(!f.IsOpen()) {
    Log(Form("File %s is not opened, something goes wrong!",fileName.Data()));
    return 0;
  }

  const Int_t nMod=12; // 1:5 modules
  const Int_t nCol=48; //1:56 columns in each module
  Int_t nRow=24; //1:64 rows in each module
  const Int_t nRowHalfSM = 12; //Supermodules 11 and 12 are half supermodules

  Double_t coeff;
  char hnam[80];
  TH1F* histo=0;


  //Get reference histogram
  TList * keylist = f.GetListOfKeys();
  Int_t nkeys   = f.GetNkeys();
  Bool_t ok = kFALSE;
  TKey  *key;
  TString refHistoName= "";
  Int_t ikey = 0;
  Int_t counter = 0;
  TH1F* hRef = new TH1F();
  
  //Check if the file contains any histogram
  if(nkeys< 2){
    Log(Form("Not enough histograms for calibration, nhist = %d",nkeys));
    return 1;
  }
  
  while(!ok){
    ikey = gRandom->Integer(nkeys);
    key = (TKey*)keylist->At(ikey);
    refHistoName = key->GetName();
    hRef = (TH1F*)f.Get(refHistoName);
    counter++;
    // Check if the reference has too little statistics and 
    // if the histogram has the correct name (2 kinds, mod#col#row for 
    // reference here, and mod#, see AliEMCALHistoProducer.
    if(refHistoName.Contains("col") && hRef->GetEntries()>2 && hRef->GetMean()>0) 
      ok=kTRUE;
    if(!ok && counter >= nMod*nCol*nRow+nMod){
      Log("No histogram with enough statistics for reference");
      return 1;
    }
  }
  
  Double_t refMean=hRef->GetMean();

  // Calculates relative calibration coefficients for all non-zero channels
  AliEMCALCalibData calibData;

  for(Int_t mod=0; mod<nMod; mod++) {
    if(mod > 10) nRow = nRowHalfSM ;
    for(Int_t col=0; col<nCol; col++) {
      for(Int_t row=0; row<nRow; row++) {
        sprintf(hnam,"mod%dcol%drow%d",mod,col,row);
        histo = (TH1F*)f.Get(hnam);
	//TODO: dead channels exclusion!
        if(histo && histo->GetMean() > 0) {
	  coeff = histo->GetMean()/refMean;
	  calibData.SetADCchannel(mod+1,col+1,row+1,1./coeff);
	  AliDebug(1,Form("mod %d col %d row %d  coeff %f\n",mod,col,row,coeff));
	}
        else
          calibData.SetADCchannel(mod+1,col+1,row+1,-111); 
      }
    }
  }

  AliCDBMetaData metaData;
  Int_t result = Store("Calib", "Data", &calibData, &metaData);

  f.Close();

  return result;

}
