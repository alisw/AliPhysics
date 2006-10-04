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

ClassImp(AliPHOSPreprocessor)

//_______________________________________________________________________________________
AliPHOSPreprocessor::AliPHOSPreprocessor() :
AliPreprocessor("PHOS",0)
{
  //default constructor
}

//_______________________________________________________________________________________
AliPHOSPreprocessor::AliPHOSPreprocessor(const char* detector, AliShuttleInterface* shuttle):
AliPreprocessor(detector,shuttle)
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

  const char* fileName = GetFile(kDAQ, "AMPLITUDES", "GDC");
  AliInfo(Form("Got filename: %s",fileName));

  TFile f(fileName);

  if(!f.IsOpen()) {
    AliInfo(Form("File %s is not opened, something goes wrong!",fileName));
    return 0;
  }

  const Int_t nMod=5; // 1:5 modules
  const Int_t nCol=56; //1:56 columns in each module
  const Int_t nRow=64; //1:64 rows in each module

  Double_t coeff;
  char hnam[80];
  TH1F* histo=0;

  // Generate name of the reference histogram
  TString refHistoName = "mod";
  refHistoName += gRandom->Integer(nMod)+1;
  refHistoName += "col";
  refHistoName += gRandom->Integer(nCol)+1;
  refHistoName += "row";
  refHistoName += gRandom->Integer(nRow)+1;
  TH1F* hRef = (TH1F*)f.Get(refHistoName);

  // If the reference histogram does not exist or has too little statistics,
  // it is better to give up preprocessing
  if(!hRef || hRef->GetEntries()<2) {
    AliInfo(Form("Cannot get reference histogram %s",refHistoName.Data()));
    return 0;
  }

  AliPHOSEmcCalibData calibData;
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
	  calibData.SetADCchannelEmc(mod+1,col+1,row+1,0.001/coeff);
	  AliInfo(Form("mod %d col %d row %d  coeff %f\n",mod,col,row,coeff));
	}
        else
          calibData.SetADCchannelEmc(mod+1,col+1,row+1,-111); 
      }
    }
  }

  AliCDBMetaData metaData;
  Int_t result = Store("Calib", "Data", &calibData, &metaData);

  f.Close();
  return result;

}
