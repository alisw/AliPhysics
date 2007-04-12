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


//-------------------------------------------------------
//          Implementation of the TPC pedestal and noise calibration
//
//   Origin: Jens Wiechula, Marian Ivanov   J.Wiechula@gsi.de, Marian.Ivanov@cern.ch
// 
// 
//-------------------------------------------------------


/* $Id$ */



//Root includes
#include <TObjArray.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TH2S.h>
#include <TH1S.h>
#include <TString.h>
#include <TVectorF.h>
#include <TMath.h>
#include <TF1.h>
#include <TRandom.h>
#include <TROOT.h>
#include <TDirectory.h>
#include <TSystem.h>
#include <TFile.h>


//AliRoot includes
#include "AliRawReader.h"
#include "AliRawReaderRoot.h"
#include "AliTPCRawStream.h"
#include "AliTPCCalROC.h"
#include "AliTPCCalPad.h"
#include "AliTPCROC.h"
#include "AliTPCCalibPedestal.h"
#include "AliMathBase.h"

#include "TTreeStream.h"



ClassImp(AliTPCCalibPedestal) /*FOLD00*/







AliTPCCalibPedestal::AliTPCCalibPedestal() : /*FOLD00*/
  TObject(),
  fFirstTimeBin(60),
  fLastTimeBin(1000),
  fAdcMin(1),
  fAdcMax(100),
  fROC(AliTPCROC::Instance()),
  fCalRocArrayPedestal(72),
  fCalRocArrayRMS(72),
  fHistoPedestalArray(72),
  fDebugStreamer(0),
  fDebugLevel(0)
{
    //
    // AliTPCSignal default constructor
    //
    //debug stream
    TDirectory *backup = gDirectory;
    fDebugStreamer = new TTreeSRedirector("deb2.root");
    if ( backup ) backup->cd();  //we don't want to be cd'd to the debug streamer    
}

//_____________________________________________________________________
AliTPCCalibPedestal::~AliTPCCalibPedestal() /*FOLD00*/
{
  //
  // destructor
  //
    if ( fDebugStreamer ) delete fDebugStreamer;
    delete fROC;
}




//_____________________________________________________________________
Int_t AliTPCCalibPedestal::Update(const Int_t icsector, /*FOLD00*/
				const Int_t icRow,
				const Int_t icPad,
				const Int_t icTimeBin,
				const Float_t csignal)
{
    //
    // Signal filling methode 
    //
    if ( (icTimeBin>fLastTimeBin) || (icTimeBin<fFirstTimeBin)   ) return 0;

    Int_t iChannel  = fROC->GetRowIndexes(icsector)[icRow]+icPad; //  global pad position in sector

    // dirty and fast filling methode. No boundary checking!!!
    // Attention: the entry counter of the histogram is not increased
    //            this means that e.g. the colz draw option gives an empty plot
    Int_t bin = (iChannel+1)*(fAdcMax-fAdcMin+2)+((Int_t)csignal-fAdcMin+1);

    GetHistoPedestal(icsector,kTRUE)->GetArray()[bin]++;

    return 0;
}
//_____________________________________________________________________
Bool_t AliTPCCalibPedestal::ProcessEvent(AliRawReader *rawReader) /*FOLD00*/
{
  //
  //  simple event processing loop
  //


    AliTPCRawStream input(rawReader);

  rawReader->Select("TPC");

  input.SetOldRCUFormat(1);
  printf("start event processing\n");


  Bool_t withInput = kFALSE;

  while (input.Next()) {

      Int_t isector  = input.GetSector();                       //  current sector
      Int_t iRow     = input.GetRow();                          //  current row
      Int_t iPad     = input.GetPad();                          //  current pad
      Int_t iTimeBin = input.GetTime();                         //  current time bin
      Float_t signal = input.GetSignal();                       //  current ADC signal

      Update(isector,iRow,iPad,iTimeBin,signal);
      withInput = kTRUE;
  }

  printf("end event processing\n");
  if ( fDebugLevel>0 )
      fDebugStreamer->GetFile()->Write();
  return withInput;
}
//_____________________________________________________________________
Bool_t AliTPCCalibPedestal::TestEvent() /*FOLD00*/
{
  //
  //  Test event loop
  //

    gRandom->SetSeed(0);

    for (UInt_t iSec=0; iSec<72; iSec++){
        if (iSec%36>0) continue;
	for (UInt_t iRow=0; iRow < fROC->GetNRows(iSec); iRow++){
	    for (UInt_t iPad=0; iPad < fROC->GetNPads(iSec,iRow); iPad++){
		for (UInt_t iTimeBin=0; iTimeBin<1024; iTimeBin++){
		    Float_t signal=(Int_t)(iRow+5+gRandom->Gaus(0,.7));
		    if ( signal>0 )Update(iSec,iRow,iPad,iTimeBin,signal);
		}
	    }
	}
    }
    return kTRUE;
}
//_____________________________________________________________________
TH2S* AliTPCCalibPedestal::GetHisto(Int_t sector, TObjArray *arr, /*FOLD00*/
				  Int_t nbinsY, Float_t ymin, Float_t ymax,
				  Char_t *type, Bool_t force)
{
    //
    // return pointer to Q histogram
    // if force is true create a new histogram if it doesn't exist allready
    //
    if ( !force || arr->UncheckedAt(sector) )
	return (TH2S*)arr->UncheckedAt(sector);

    // if we are forced and histogram doesn't yes exist create it
    Char_t name[255], title[255];

    sprintf(name,"hCalib%s%.2d",type,sector);
    sprintf(title,"%s calibration histogram sector %.2d",type,sector);

    // new histogram with Q calib information. One value for each pad!
    TH2S* hist = new TH2S(name,title,
			  nbinsY, ymin, ymax,
			  fROC->GetNChannels(sector),0,fROC->GetNChannels(sector)
			 );
    hist->SetDirectory(0);
    arr->AddAt(hist,sector);
    return hist;
}
//_____________________________________________________________________
TH2S* AliTPCCalibPedestal::GetHistoPedestal(Int_t sector, Bool_t force) /*FOLD00*/
{
    //
    // return pointer to T0 histogram
    // if force is true create a new histogram if it doesn't exist allready
    //
    TObjArray *arr = &fHistoPedestalArray;
    return GetHisto(sector, arr, fAdcMax-fAdcMin, fAdcMin, fAdcMax, "Pedestal", force);
}
//_____________________________________________________________________
AliTPCCalROC* AliTPCCalibPedestal::GetCalRoc(Int_t sector, TObjArray* arr, Bool_t force) /*FOLD00*/
{
    //
    // return pointer to ROC Calibration
    // if force is true create a new histogram if it doesn't exist allready
    //
    if ( !force || arr->UncheckedAt(sector) )
	return (AliTPCCalROC*)arr->UncheckedAt(sector);

    // if we are forced and histogram doesn't yes exist create it

    // new AliTPCCalROC for T0 information. One value for each pad!
    AliTPCCalROC *croc = new AliTPCCalROC(sector);
    //init values
    for ( UInt_t iChannel = 0; iChannel<croc->GetNchannels(); iChannel++){
	croc->SetValue(iChannel, 0);
    }
    arr->AddAt(croc,sector);
    return croc;
}
//_____________________________________________________________________
AliTPCCalROC* AliTPCCalibPedestal::GetCalRocPedestal(Int_t sector, Bool_t force) /*FOLD00*/
{
    //
    // return pointer to Carge ROC Calibration
    // if force is true create a new histogram if it doesn't exist allready
    //
    TObjArray *arr = &fCalRocArrayPedestal;
    return GetCalRoc(sector, arr, force);
}
//_____________________________________________________________________
AliTPCCalROC* AliTPCCalibPedestal::GetCalRocRMS(Int_t sector, Bool_t force) /*FOLD00*/
{
    //
    // return pointer to signal width ROC Calibration
    // if force is true create a new histogram if it doesn't exist allready
    //
    TObjArray *arr = &fCalRocArrayRMS;
    return GetCalRoc(sector, arr, force);
}
//_____________________________________________________________________
void AliTPCCalibPedestal::Analyse() /*FOLD00*/
{
    //
    //  Calculate calibration constants
    //

    Int_t nbinsAdc = fAdcMax-fAdcMin;

    TH1F *py = new TH1F("htemp_py","htemp_py", nbinsAdc, fAdcMin, fAdcMax);
    TF1 *gaus = new TF1("fit","gaus");
    TVectorD param(3);

    Float_t *array_py=0;
    Short_t *array_hP=0;

    array_py = py->GetArray();

    for (Int_t iSec=0; iSec<72; iSec++){
	TH2S *hP = GetHistoPedestal(iSec);
        if ( !hP ) continue;

	AliTPCCalROC *rocPedestal = GetCalRocPedestal(iSec,kTRUE);
	AliTPCCalROC *rocRMS      = GetCalRocRMS(iSec,kTRUE);

	array_hP = hP->GetArray();
        UInt_t nChannels = fROC->GetNChannels(iSec);

	for (UInt_t iChannel=0; iChannel<nChannels; iChannel++){

            // set bin content of py in a dirty but fast way
	    for (Int_t iAdc=0;iAdc<nbinsAdc;iAdc++)
		array_py[iAdc+1] = (Float_t)array_hP[(iChannel+1)*(nbinsAdc+2)+(iAdc+1)];

	    gaus->SetParameters(0,0);
	    py->Fit(gaus,"nq");

	    rocPedestal->SetValue(iChannel,gaus->GetParameter(1));
            rocRMS->SetValue(iChannel,gaus->GetParameter(2));

	    //AliMathBase::FitGaus(nbinsAdc,array_hP+(iChannel+1)*(nbinsAdc+2),&param)
	    //rocPedestal->SetValue(iChannel,param[1]);
            //rocRMS->SetValue(iChannel,param[2]);
	}
    }
    delete py;
    delete gaus;
    delete fDebugStreamer;
    fDebugStreamer = 0x0;
}
//_____________________________________________________________________
void AliTPCCalibPedestal::DumpToFile(const Char_t *filename, const Char_t *dir, Bool_t append) /*FOLD00*/
{
    //
    //  Write class to file
    //

    TDirectory *backup = gDirectory;
    TString sDir(dir);
    TString option;

    if ( append )
	option = "update";
    else
        option = "recreate";

    TFile f(filename,option.Data());
    if ( !sDir.IsNull() ){
	f.mkdir(sDir.Data());
	f.cd(sDir);
    }
    gDirectory->WriteTObject(this);
    f.Close();

    if ( backup ) backup->cd();

}
