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


//Root includes
#include <TH1F.h>
#include <TH2F.h>
#include <TString.h>
#include <TMath.h>
#include <TF1.h>
#include <TRandom.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TMap.h>
//AliRoot includes
#include "AliRawReader.h"
#include "AliRawReaderRoot.h"
#include "AliRawReaderDate.h"
#include "AliTPCRawStream.h"
#include "AliTPCCalROC.h"
#include "AliTPCROC.h"
#include "AliMathBase.h"
#include "TTreeStream.h"

//date
#include "event.h"

//header file
#include "AliTPCCalibPedestal.h"


///////////////////////////////////////////////////////////////////////////////////////
//          Implementation of the TPC pedestal and noise calibration
//
//   Origin: Jens Wiechula, Marian Ivanov   J.Wiechula@gsi.de, Marian.Ivanov@cern.ch
// 
// 
// *************************************************************************************
// *                                Class Description                                  *
// *************************************************************************************
//
// Working principle:
// ------------------
// Raw pedestal data is processed by calling one of the ProcessEvent(...) functions
// (see below). These in the end call the Update(...) function, where the data is filled
// into histograms.
//
// For each ROC one TH2F histo (ROC channel vs. ADC channel) is created when
// it is filled for the first time (GetHistoPedestal(ROC,kTRUE)). All histos are stored in the
// TObjArray fHistoPedestalArray.
//
// For a fast filling of the histogram the corresponding bin number of the channel and ADC channel
// is computed by hand and the histogram array is accessed directly via its pointer.
// ATTENTION: Doing so the the entry counter of the histogram is not increased
//            this means that e.g. the colz draw option gives an empty plot unless
//	    calling 'histo->SetEntries(1)' before drawing.
//
// After accumulating the desired statistics the Analyse() function has to be called.
// Whithin this function the pedestal and noise values are calculated for each pad, using
// the fast gaus fit function  AliMathBase::FitGaus(...), and the calibration
// storage classes (AliTPCCalROC) are filled for each ROC.
// The calibration information is stored in the TObjArrays fCalRocArrayPedestal and fCalRocArrayRMS;
//
//
//
// User interface for filling data:
// --------------------------------
//
// To Fill information one of the following functions can be used:
//
// Bool_t ProcessEvent(eventHeaderStruct *event);
//   - process Date event
//   - use AliTPCRawReaderDate and call ProcessEvent(AliRawReader *rawReader)
//
// Bool_t ProcessEvent(AliRawReader *rawReader);
//  - process AliRawReader event
//   - use AliTPCRawStream to loop over data and call ProcessEvent(AliTPCRawStream *rawStream)
//
// Bool_t ProcessEvent(AliTPCRawStream *rawStream);
//   - process event from AliTPCRawStream
//   - call Update function for signal filling
//
// Int_t Update(const Int_t isector, const Int_t iRow, const Int_t
//              iPad,  const Int_t iTimeBin, const Float_t signal);
//   - directly  fill signal information (sector, row, pad, time bin, pad)
//     to the reference histograms
//
// It is also possible to merge two independently taken calibrations using the function
//
// void Merge(AliTPCCalibPedestal *ped)
//   - copy histograms in 'ped' if the do not exist in this instance
//   - Add histograms in 'ped' to the histograms in this instance if the allready exist
//   - After merging call Analyse again!
//
//
//
// -- example: filling data using root raw data:
// void fillPedestal(Char_t *filename)
// {
//    rawReader = new AliRawReaderRoot(fileName);
//    if ( !rawReader ) return;
//    AliTPCCalibPedestal *calib = new AliTPCCalibPedestal;
//    while (rawReader->NextEvent()){
//      calib->ProcessEvent(rawReader);
//    }
//    calib->Analyse();
//    calib->DumpToFile("PedestalData.root");
//    delete rawReader;
//    delete calib;
// }
//
//
// What kind of information is stored and how to retrieve them:
// ------------------------------------------------------------
//
// - Accessing the 'Reference Histograms' (pedestal distribution histograms):
//
// TH2F *GetHistoPedestal(Int_t sector);
//
// - Accessing the calibration storage objects:
//
// AliTPCCalROC *GetCalRocPedestal(Int_t sector);  - for the pedestal values, mean from gaus fit
// AliTPCCalROC *GetCalRocSigma(Int_t sector);     - for the Noise values, sigma from guas fit
// AliTPCCalROC *GetCalRocMean(Int_t sector);  - for the pedestal values, truncated mean
// AliTPCCalROC *GetCalRocRMS(Int_t sector);     - for the Noise values, rms from truncated mean
//
// example for visualisation:
// if the file "PedestalData.root" was created using the above example one could do the following:
//
// TFile filePedestal("PedestalData.root")
// AliTPCCalibPedestal *ped = (AliTPCCalibPedestal*)filePedestal->Get("AliTPCCalibPedestal");
// ped->GetCalRocPedestal(0)->Draw("colz");
// ped->GetCalRocRMS(0)->Draw("colz");
//
// or use the AliTPCCalPad functionality:
// AliTPCCalPad padPedestal(ped->GetCalPadPedestal());
// AliTPCCalPad padNoise(ped->GetCalPadRMS());
// padPedestal->MakeHisto2D()->Draw("colz");  //Draw A-Side Pedestal Information
// padNoise->MakeHisto2D()->Draw("colz");  //Draw A-Side Noise Information
//
/*
 example: fill pedestal with gausschen noise
 AliTPCCalibPedestal ped;
 ped.TestEvent();
 ped.Analyse();
 //Draw output;
 TCanvas* c1 = new TCanvas;
 c1->Divide(1,2);
 c1->cd(1);
 ped.GetHistoPedestal(0)->SetEntries(1); //needed in order for colz to work, reason: fast filling does not increase the entries counter
 ped.GetHistoPedestal(0)->Draw("colz");
 c1->cd(2);
 ped.GetHistoPedestal(36)->SetEntries(1); //needed in order for colz to work, reason: fast filling does not increase the entries counter
 ped.GetHistoPedestal(36)->Draw("colz");
 TCanvas* c2 = new TCanvas;
 c2->Divide(2,2);
 c2->cd(1);
 ped.GetCalRocPedestal(0)->Draw("colz");
 c2->cd(2);
 ped.GetCalRocRMS(0)->Draw("colz");
 c2->cd(3);
 ped.GetCalRocPedestal(36)->Draw("colz");
 c2->cd(4);
 ped.GetCalRocRMS(36)->Draw("colz");
*/
//
// Time dependent pedestals:
//
// If wished there is the possibility to calculate for each channel and time bin
// the mean pedestal [pedestals(t)]. This is done by
//
// 1) setting SetTimeAnalysis(kTRUE),
// 2) processing the data by looping over the events using ProcessEvent(..)
// 3) calling the Analyse() and AnalyseTime(nevents) functions (providing nevents)
// 4) getting the pedestals(t) using   TArrayF **timePed = calibPedestal.GetTimePedestals();
// 5) looking at values using   timePed[row][pad].At(timebin)
//
// This functionality is intended to be used on an LDC bu the detector algorithm
// (TPCPEDESTALda) to generate a data set used for configuration of the pattern
// memory for baseline subtraction in the ALTROs. Later the information should also
// be stored as reference data.
//


ClassImp(AliTPCCalibPedestal)

AliTPCCalibPedestal::AliTPCCalibPedestal() : 
  AliTPCCalibRawBase(),
  fAdcMin(1),
  fAdcMax(100),
  fAnaMeanDown(0.),
  fAnaMeanUp(1.),
  fTimeAnalysis(kFALSE),
  fCalRocArrayPedestal(72),
  fCalRocArraySigma(72),
  fHistoPedestalArray(72),
  fTimeSignal(NULL),
  fCalRocArrayMean(72),
  fCalRocArrayRMS(72)
{
  //
  // default constructor
  //
  SetNameTitle("AliTPCCalibPedestal","AliTPCCalibPedestal");
  fFirstTimeBin=60;
  fLastTimeBin=1000;
}


//_____________________________________________________________________
AliTPCCalibPedestal::AliTPCCalibPedestal(const AliTPCCalibPedestal &ped) : 
  AliTPCCalibRawBase(ped),
  fAdcMin(ped.GetAdcMin()),
  fAdcMax(ped.GetAdcMax()),
  fAnaMeanDown(ped.fAnaMeanDown),
  fAnaMeanUp(ped.fAnaMeanUp),
  fTimeAnalysis(ped.fTimeAnalysis),
  fCalRocArrayPedestal(72),
  fCalRocArraySigma(72),
  fHistoPedestalArray(72),
  fTimeSignal(ped.fTimeSignal),
  fCalRocArrayMean(72),
  fCalRocArrayRMS(72)
{
  //
  // copy constructor
  //
  for (Int_t iSec = 0; iSec < 72; ++iSec){
    const AliTPCCalROC *calPed = (AliTPCCalROC*)ped.fCalRocArrayPedestal.UncheckedAt(iSec);
    const AliTPCCalROC *calRMS = (AliTPCCalROC*)ped.fCalRocArrayRMS.UncheckedAt(iSec);
    const TH2F         *hPed   = (TH2F*)ped.fHistoPedestalArray.UncheckedAt(iSec);
    
    if ( calPed != 0x0 ) fCalRocArrayPedestal.AddAt(new AliTPCCalROC(*calPed), iSec);
    if ( calRMS != 0x0 ) fCalRocArrayRMS.AddAt(new AliTPCCalROC(*calRMS), iSec);
    
    if ( hPed != 0x0 ){
      TH2F *hNew = new TH2F(*hPed);
      hNew->SetDirectory(0);
      fHistoPedestalArray.AddAt(hNew,iSec);
    }
  }
}
AliTPCCalibPedestal::AliTPCCalibPedestal(const TMap *config): 
  AliTPCCalibRawBase(),
  fAdcMin(1),
  fAdcMax(100),
  fAnaMeanDown(0.),
  fAnaMeanUp(1.),
  fTimeAnalysis(kFALSE),
  fCalRocArrayPedestal(72),
  fCalRocArraySigma(72),
  fHistoPedestalArray(72),
  fTimeSignal(NULL),
  fCalRocArrayMean(72),
  fCalRocArrayRMS(72)  
{
 //
 // This constructor uses a TMap for setting some parametes
 //
  SetNameTitle("AliTPCCalibPedestal","AliTPCCalibPedestal");
  fFirstTimeBin=60;
  fLastTimeBin=1000;
  if (config->GetValue("FirstTimeBin")) fFirstTimeBin = ((TObjString*)config->GetValue("FirstTimeBin"))->GetString().Atoi();
  if (config->GetValue("LastTimeBin"))  fLastTimeBin = ((TObjString*)config->GetValue("LastTimeBin"))->GetString().Atoi();
  if (config->GetValue("AdcMin"))       fAdcMin = ((TObjString*)config->GetValue("AdcMin"))->GetString().Atoi();
  if (config->GetValue("AdcMax"))       fAdcMax = ((TObjString*)config->GetValue("AdcMax"))->GetString().Atoi();
  if (config->GetValue("TimeAnalysis")) SetTimeAnalysis(((TObjString*)config->GetValue("TimeAnalysis"))->GetString().Atoi());
} 


//_____________________________________________________________________
AliTPCCalibPedestal& AliTPCCalibPedestal::operator = (const  AliTPCCalibPedestal &source)
{
  //
  // assignment operator
  //
  if (&source == this) return *this;
  new (this) AliTPCCalibPedestal(source);

  return *this;
}


//_____________________________________________________________________
AliTPCCalibPedestal::~AliTPCCalibPedestal() 
{
  //
  // destructor
  //

  fCalRocArrayPedestal.Delete();
  fCalRocArrayRMS.Delete();
  fCalRocArraySigma.Delete();
  fHistoPedestalArray.Delete();

  if ( fTimeSignal ) {
    for (Int_t i = 0; i < 159; i++) {
      delete [] fTimeSignal[i];
      fTimeSignal[i] = 0;
    }
    delete [] fTimeSignal;
    fTimeSignal = 0;
  }

  // do not delete fMapping, because we do not own it.

}


//_____________________________________________________________________
void AliTPCCalibPedestal::SetTimeAnalysis(Bool_t time)
{
  //
  // Use time dependent analysis: Pedestals are analysed as a function
  // of the drift time. There is one mean value generated for each time
  // bin and each channel. It can be used as reference data and for
  // configuration of the ALTRO pattern memory for baseline subtraction.
  //
  // ATTENTION: Use only on LDC in TPCPEDESTALda! On a LDC we get data
  // only from one sector. For the full TPC we would need a lot of
  // memory (36*159*140*1024*4bytes = 3.3GB)!
  //

  fTimeAnalysis = time;

  if ( !fTimeAnalysis ) return;

  // prepare array for one sector (159*140*1024*4bytes = 92MB):
  fTimeSignal = new TArrayF*[159];
  for (Int_t i = 0; i < 159; i++) {  // padrows
    fTimeSignal[i] = new TArrayF[140];
    for (Int_t j = 0; j < 140; j++) {  // pads per row
      fTimeSignal[i][j].Set(1024);
      for (Int_t k = 0; k < 1024; k++) {  // time bins per pad
	fTimeSignal[i][j].AddAt(0., k);
      }
    }      
  }
}


//_____________________________________________________________________
Int_t AliTPCCalibPedestal::Update(const Int_t icsector, 
				  const Int_t icRow,
				  const Int_t icPad,
				  const Int_t icTimeBin,
				  const Float_t csignal)
{
  //
  // Signal filling method
  //
  if (icRow<0) return 0;
  if (icPad<0) return 0;
  if (icTimeBin<0) return 0;
 
  // Time dependent pedestals
  if ( fTimeAnalysis ) {
    if ( icsector < 36 ) // IROC
      fTimeSignal[icRow][icPad].AddAt(fTimeSignal[icRow][icPad].At(icTimeBin)+csignal, icTimeBin);
    else 
      fTimeSignal[icRow+63][icPad].AddAt(fTimeSignal[icRow+63][icPad].At(icTimeBin)+csignal, icTimeBin);
  }
  //return if we are out of the specified time bin or adc range
  if ( (icTimeBin>fLastTimeBin) || (icTimeBin<fFirstTimeBin) ) return 0;
  if ( ((Int_t)csignal>fAdcMax) || ((Int_t)csignal<fAdcMin)  ) return 0;

  Int_t iChannel  = fROC->GetRowIndexes(icsector)[icRow]+icPad; //  global pad position in sector

  // fast filling method
  // Attention: the entry counter of the histogram is not increased
  //            this means that e.g. the colz draw option gives an empty plot
  Int_t bin = (iChannel+1)*(fAdcMax-fAdcMin+2)+((Int_t)csignal-fAdcMin+1);

  GetHistoPedestal(icsector,kTRUE)->GetArray()[bin]++;

  return 0;
}


//_____________________________________________________________________
Bool_t AliTPCCalibPedestal::TestEvent() 
{
  //
  //  Test event loop
  // fill one oroc and one iroc with random gaus
  //

  gRandom->SetSeed(0);

  for (UInt_t iSec=0; iSec<72; ++iSec){
    if (iSec%36>0) continue;
    for (UInt_t iRow=0; iRow < fROC->GetNRows(iSec); ++iRow){
      for (UInt_t iPad=0; iPad < fROC->GetNPads(iSec,iRow); ++iPad){
        for (UInt_t iTimeBin=0; iTimeBin<1024; ++iTimeBin){
          Float_t signal=(Int_t)(iRow+3+gRandom->Gaus(0,.7));
          if ( signal>0 )Update(iSec,iRow,iPad,iTimeBin,signal);
        }
      }
    }
  }
  return kTRUE;
}


//_____________________________________________________________________
TH2F* AliTPCCalibPedestal::GetHisto(Int_t sector, TObjArray *arr, 
				    Int_t nbinsY, Float_t ymin, Float_t ymax,
				    const Char_t *type, Bool_t force)
{
    //
    // return pointer to Q histogram
    // if force is true create a new histogram if it doesn't exist allready
    //
    if ( !force || arr->UncheckedAt(sector) )
      return (TH2F*)arr->UncheckedAt(sector);

    // if we are forced and histogram doesn't yes exist create it
    // new histogram with Q calib information. One value for each pad!
    TH2F* hist = new TH2F(Form("hCalib%s%.2d",type,sector),
                          Form("%s calibration histogram sector %.2d;ADC channel;Channel (pad)",type,sector),
                          nbinsY, ymin, ymax,
                          fROC->GetNChannels(sector),0,fROC->GetNChannels(sector)
                         );
    hist->SetDirectory(0);
    arr->AddAt(hist,sector);
    return hist;
}


//_____________________________________________________________________
TH2F* AliTPCCalibPedestal::GetHistoPedestal(Int_t sector, Bool_t force) 
{
    //
    // return pointer to T0 histogram
    // if force is true create a new histogram if it doesn't exist allready
    //
    TObjArray *arr = &fHistoPedestalArray;
    return GetHisto(sector, arr, fAdcMax-fAdcMin, fAdcMin, fAdcMax, "Pedestal", force);
}


//_____________________________________________________________________
AliTPCCalROC* AliTPCCalibPedestal::GetCalRoc(Int_t sector, TObjArray* arr, Bool_t force) 
{
    //
    // return pointer to ROC Calibration
    // if force is true create a new histogram if it doesn't exist allready
    //
    if ( !force || arr->UncheckedAt(sector) )
	return (AliTPCCalROC*)arr->UncheckedAt(sector);

    // if we are forced and the histogram doesn't yet exist create it

    // new AliTPCCalROC for T0 information. One value for each pad!
    AliTPCCalROC *croc = new AliTPCCalROC(sector);
    arr->AddAt(croc,sector);
    return croc;
}


//_____________________________________________________________________
AliTPCCalROC* AliTPCCalibPedestal::GetCalRocPedestal(Int_t sector, Bool_t force) 
{
    //
    // return pointer to ROC with Pedestal data
    // if force is true create a new histogram if it doesn't exist allready
    //
    TObjArray *arr = &fCalRocArrayPedestal;
    return GetCalRoc(sector, arr, force);
}


//_____________________________________________________________________
AliTPCCalROC* AliTPCCalibPedestal::GetCalRocSigma(Int_t sector, Bool_t force) 
{
    //
    // return pointer to  ROC with signal witdth in sigma
    // if force is true create a new histogram if it doesn't exist allready
    //
    TObjArray *arr = &fCalRocArraySigma;
    return GetCalRoc(sector, arr, force);
}
//_____________________________________________________________________
AliTPCCalROC* AliTPCCalibPedestal::GetCalRocMean(Int_t sector, Bool_t force)
{
  //
    // return pointer to ROC with signal mean information
    // if force is true create a new histogram if it doesn't exist allready
  //
  TObjArray *arr = &fCalRocArrayMean;
  return GetCalRoc(sector, arr, force);
}

//_____________________________________________________________________
AliTPCCalROC* AliTPCCalibPedestal::GetCalRocRMS(Int_t sector, Bool_t force) 
{
  //
    // return pointer to signal width ROC Calibration
    // if force is true create a new histogram if it doesn't exist allready
  //
  TObjArray *arr = &fCalRocArrayRMS;
  return GetCalRoc(sector, arr, force);
}


//_____________________________________________________________________
void AliTPCCalibPedestal::Merge(AliTPCCalibPedestal * const ped)
{
  //
  //  Merge reference histograms of sig to the current AliTPCCalibSignal
  //
  MergeBase(ped);
  // merge histograms
  for (Int_t iSec=0; iSec<72; ++iSec){
    TH2F *hRefPedMerge   = ped->GetHistoPedestal(iSec);
    
    if ( hRefPedMerge ){
      TDirectory *dir = hRefPedMerge->GetDirectory(); hRefPedMerge->SetDirectory(0);
      TH2F *hRefPed   = GetHistoPedestal(iSec);
      if ( hRefPed ) hRefPed->Add(hRefPedMerge);
      else {
        TH2F *hist = new TH2F(*hRefPedMerge);
        hist->SetDirectory(0);
        fHistoPedestalArray.AddAt(hist, iSec);
      }
      hRefPedMerge->SetDirectory(dir);
    }
  }
  
  // merge array
  // ...
  
}

//_____________________________________________________________________
Long64_t AliTPCCalibPedestal::Merge(TCollection * const list)
{
  //
  // Merge all objects of this type in list
  //
  
  Long64_t nmerged=1;
  
  TIter next(list);
  AliTPCCalibPedestal *ce=0;
  TObject *o=0;
  
  while ( (o=next()) ){
    ce=dynamic_cast<AliTPCCalibPedestal*>(o);
    if (ce){
      Merge(ce);
      ++nmerged;
    }
  }
  
  return nmerged;
}

//_____________________________________________________________________
void AliTPCCalibPedestal::Analyse() 
{
  //
  //  Calculate calibration constants
  //

  Int_t nbinsAdc = fAdcMax-fAdcMin;

  TVectorD param(4);
  TMatrixD dummy(3,3);

  TH1F *hChannel=new TH1F("hChannel","hChannel",nbinsAdc,fAdcMin,fAdcMax);
  
  Float_t *arrayhP=0;  

  for (Int_t iSec=0; iSec<72; ++iSec){
    TH2F *hP = GetHistoPedestal(iSec);
    if ( !hP ) continue;

    AliTPCCalROC *rocPedestal = GetCalRocPedestal(iSec,kTRUE);
    AliTPCCalROC *rocSigma    = GetCalRocSigma(iSec,kTRUE);
    AliTPCCalROC *rocMean     = GetCalRocMean(iSec,kTRUE);
    AliTPCCalROC *rocRMS      = GetCalRocRMS(iSec,kTRUE);

    arrayhP = hP->GetArray();
    UInt_t nChannels = fROC->GetNChannels(iSec);

    for (UInt_t iChannel=0; iChannel<nChannels; ++iChannel){
      Int_t offset = (nbinsAdc+2)*(iChannel+1)+1;
      //calculate mean and sigma using a gaus fit
      //Double_t ret =
      AliMathBase::FitGaus(arrayhP+offset,nbinsAdc,fAdcMin,fAdcMax,&param,&dummy);
      // if the fitting failed set noise and pedestal to 0
      // is now done in AliMathBase::FitGaus !
//       if ( ret == -4 ) {
// 	param[1]=0;
// 	param[2]=0;
//       }
      if ( param[1]<fAdcMin || param[1]>fAdcMax ){
        param[1]=0;
        param[2]=0;
      }
      rocPedestal->SetValue(iChannel,param[1]);
      rocSigma->SetValue(iChannel,param[2]);
      //calculate mean and RMS using a truncated means
      hChannel->Set(nbinsAdc+2,arrayhP+offset-1);
      hChannel->SetEntries(param[3]);
      param[1]=0;
      param[2]=0;
      if ( param[3]>0 ) AliMathBase::TruncatedMean(hChannel,&param,fAnaMeanDown,fAnaMeanUp);
      rocMean->SetValue(iChannel,param[1]);
      rocRMS->SetValue(iChannel,param[2]);
    }
  }
  delete hChannel;
}


//_____________________________________________________________________
void AliTPCCalibPedestal::AnalyseTime(Int_t nevents)
{
  //
  // Calculate for each channel and time bin the mean pedestal. This
  // is used on LDC by TPCPEDESTALda to generate data used for configuration
  // of the pattern memory for baseline subtraction in the ALTROs.
  //

  if ( nevents <= 0 ) return;
  if ( fTimeAnalysis ) {
    for (Int_t i = 0; i < 159; i++) {  // padrows
      for (Int_t j = 0; j < 140; j++) {  // pads per row
	for (Int_t k = 0; k < 1024; k++) {  // time bins per pad
	  fTimeSignal[i][j].AddAt(fTimeSignal[i][j].At(k)/(Float_t)nevents, k);
	}
      }
    }
  }
}
