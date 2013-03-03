/**************************************************************************
 * Copyright(c) 2007-08, ALICE Experiment at CERN, All rights reserved. *
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


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Class for Evaluation and Validation of the ALTRO Tail Cancelation Filter  //
// (TCF) parameters out of TPC Raw data                                      //
//                                                                           //
// Author: Stefan Rossegger, Simon Feigl                                     //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliTPCCalibTCF.h"

#include <TObject.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TKey.h>
#include <TStyle.h>
#include <TMinuit.h>
#include <TH1F.h>
#include <TH2F.h>
#include <AliSysInfo.h>

#include <TMath.h>
#include <TNtuple.h>
#include <TEntryList.h>
#include "AliRawReaderRoot.h"
#include "AliRawHLTManager.h"
#include "AliTPCRawStreamV3.h"
#include "AliTPCROC.h"

#include "AliTPCAltroEmulator.h"

#include "AliTPCmapper.h"
#include <fstream>

ClassImp(AliTPCCalibTCF)
  
AliTPCCalibTCF::AliTPCCalibTCF() :
  TNamed(),
  fGateWidth(50),
  fSample(900),
  fPulseLength(400),
  fLowPulseLim(30),
  fUpPulseLim(900),
  fRMSLim(1.0),
  fRatioIntLim(2)

{
  //
  //  AliTPCCalibTCF standard constructor
  //
}
  
//_____________________________________________________________________________
AliTPCCalibTCF::AliTPCCalibTCF(Int_t gateWidth, Int_t sample, Int_t pulseLength, Int_t lowPulseLim, Int_t upPulseLim, Double_t rmsLim, Double_t ratioIntLim) : 
  TNamed(),
  fGateWidth(gateWidth),
  fSample(sample),
  fPulseLength(pulseLength),
  fLowPulseLim(lowPulseLim),
  fUpPulseLim(upPulseLim),
  fRMSLim(rmsLim),
  fRatioIntLim(ratioIntLim)
{
  //
  //  AliTPCCalibTCF constructor with specific (non-standard) thresholds
  //
}
  
//_____________________________________________________________________________
AliTPCCalibTCF::AliTPCCalibTCF(const AliTPCCalibTCF &tcf) : 
  TNamed(tcf),
  fGateWidth(tcf.fGateWidth),
  fSample(tcf.fSample),
  fPulseLength(tcf.fPulseLength),
  fLowPulseLim(tcf.fLowPulseLim),
  fUpPulseLim(tcf.fUpPulseLim),
  fRMSLim(tcf.fRMSLim),
  fRatioIntLim(tcf.fRatioIntLim)
{
  //
  //  AliTPCCalibTCF copy constructor
  //
}


//_____________________________________________________________________________
AliTPCCalibTCF& AliTPCCalibTCF::operator = (const AliTPCCalibTCF &source)
{
  //
  // AliTPCCalibTCF assignment operator
  //
 
  if (&source == this) return *this;
  new (this) AliTPCCalibTCF(source);

  return *this;

}

//_____________________________________________________________________________
AliTPCCalibTCF::~AliTPCCalibTCF()
{
  //
  // AliTPCCalibTCF destructor
  //
}


//_____________________________________________________________________________
void AliTPCCalibTCF::ProcessRawFileV3(const char *nameRawFile, const char *nameFileOut) {
  //
  // New RCU data format!: Standard middle of 2009 
  //
  // Loops over all events within one RawData file and collects proper pulses 
  // (according to given tresholds) per pad
  // Histograms per pad are stored in 'nameFileOut'
  //
  
  AliRawReader *rawReader = AliRawReader::Create(nameRawFile);
  if (!rawReader) {
    printf("Could not create a raw reader for %s\n",nameRawFile);
    return;
  } 

  rawReader->RewindEvents(); // just to make sure
  
  rawReader->Select("TPC");

  if (!rawReader->NextEvent()) {
    printf("no events found in %s\n",nameRawFile);
    return;
  }

  // TPC stream reader 
  AliTPCRawStreamV3 rawStream(rawReader);
  
  Int_t ievent=0;
  do {  
    AliSysInfo::AddStamp(Form("start_event_%d",ievent), ievent,-1,-1);
    printf("Reading next event ... Nr: %d\n",ievent);
    // Start the basic data extraction
    ProcessRawEventV3(rawReader, &rawStream, nameFileOut);
    AliSysInfo::AddStamp(Form("end_event_%d",ievent), ievent,-1,-1);
    ievent++;

  } while (rawReader->NextEvent());

  rawReader->~AliRawReader();
  
}

//_____________________________________________________________________________
void AliTPCCalibTCF::ProcessRawEventV3( AliRawReader *rawReader, AliTPCRawStreamV3 *rawStream, const char *nameFileOut) {
  //
  // New RCU data format!: Standard middle of 2009 
  //
  // Extracts proper pulses (according the given tresholds) within one event
  // and accumulates them into one histogram per pad. All histograms are
  // saved in the file 'nameFileOut'. 
  // The first bins of the histograms contain the following information:
  //   bin 1: Number of accumulated pulses
  //   bin 2;3;4: Sector; Row; Pad; 
  // 
  
  TFile fileOut(nameFileOut,"UPDATE");
  fileOut.cd();  
  
  TH1I *tempHis = new TH1I("tempHis","tempHis",fSample,fGateWidth,fSample+fGateWidth);
  TH1I *tempRMSHis = new TH1I("tempRMSHis","tempRMSHis",2000,0,2000);
  
  // loop over the data in this event

  while (rawStream->NextDDL() ) { 

    Int_t ddl = rawReader->GetDDLID();
    
    while (rawStream->NextChannel() ) {
      
      while (rawStream->NextBunch() ) {

	Int_t t0 = rawStream->GetStartTimeBin();
	Int_t bl = rawStream->GetBunchLength();

	if (bl<fSample+fGateWidth) continue;

	Int_t sector = rawStream->GetSector();
	Int_t row =    rawStream->GetRow();
	Int_t pad =    rawStream->GetPad();

	UShort_t *signals=(UShort_t*)rawStream->GetSignals();
	if (!signals) continue;
	
	// Write to temporary histogramm
	for (Int_t i=0;i<bl;++i) {
	  UShort_t time=t0-i;
	  UShort_t signal=signals[i];
	  if ( (fGateWidth<time) && (time<=fSample+fGateWidth) ) {
	    tempHis->SetBinContent(time-fGateWidth,signal);
	  }
	}
         
	// calculation of the pulse properties and comparison to thresholds settings
	
	Int_t max = (Int_t)tempHis->GetMaximum(FLT_MAX);
	Int_t maxpos =  tempHis->GetMaximumBin();
	
	Int_t first = (Int_t)TMath::Max(maxpos-10, 0);
	Int_t last  = TMath::Min((Int_t)maxpos+fPulseLength-10, fSample+fGateWidth);
	
	// simple baseline substraction ? better one needed ? (pedestalsubstr.?)
	// and RMS calculation with timebins before the pulse and at the end of
	// the signal 
	for (Int_t ipos = 0; ipos<6; ipos++) {
	  // before the pulse
	  tempRMSHis->Fill(tempHis->GetBinContent(first+ipos));
	}
	for (Int_t ipos = 0; ipos<20; ipos++) {
	  // at the end to get rid of pulses with serious baseline fluctuations
	  tempRMSHis->Fill(tempHis->GetBinContent(last-ipos)); 
	}
	
	Double_t baseline = tempRMSHis->GetMean();
	Double_t rms = tempRMSHis->GetRMS();
	tempRMSHis->Reset();
	
	Double_t lowLim = fLowPulseLim+baseline;
	Double_t upLim = fUpPulseLim+baseline;

	// get rid of pulses which contain gate signal and/or too much noise
	// with the help of ratio of integrals
	Double_t intHist = 0;
	Double_t intPulse = 0;
	Double_t binValue;
	for(Int_t ipos=first; ipos<=last; ipos++) {
	  binValue = TMath::Abs(tempHis->GetBinContent(ipos) - baseline);
	  intHist += binValue;
	  if(ipos>=first+5 && ipos<=first+15) {intPulse += binValue;}
	}
        
	// gets rid of high frequency noise:
	// calculating ratio (value one to the right of maximum)/(maximum)
	// has to be >= 0.1; if maximum==0 set ratio to 0.1
	Double_t maxCorr = max - baseline;
	Double_t binRatio = 0.1;
	if(TMath::Abs(maxCorr)>1e-5) {
	  binRatio = (tempHis->GetBinContent(maxpos+1) - baseline) / maxCorr;
	}
	
	// Decision if found pulse is a proper one according to given tresholds
        if (max>lowLim && max<upLim && !((last-first)<fPulseLength) && rms<fRMSLim && (intHist/intPulse)<fRatioIntLim &&intPulse>10&& (binRatio >= 0.1) ) {

	  // 1D histogramm for mean pulse per pad
	  char hname[100];
	  snprintf(hname,100,"sec%drow%dpad%d",sector,row,pad);
	  
	  TH1F *his = (TH1F*)fileOut.Get(hname);
	  
	  if (!his ) { // new entry (pulse in new pad found)
	    
	    his = new TH1F(hname,hname, fPulseLength+5, 0, fPulseLength+5);
	    his->SetBinContent(1,1);        //  pulse counter (1st pulse)
	    his->SetBinContent(2,sector);   //  sector
	    his->SetBinContent(3,row);      //  row
	    his->SetBinContent(4,pad);      //  pad	  
	    
	    for (Int_t ipos=0; ipos<last-first; ipos++){
	      Int_t signal = (Int_t)(tempHis->GetBinContent(ipos+first)-baseline);
	      his->SetBinContent(ipos+5,signal);
	    }
	    his->Write(hname);
	    printf("new  %s: Signal %d at bin %d \n", hname, max-(Int_t)baseline, maxpos+fGateWidth);
	    
	  } else {  // adding pulse to existing histogram (pad already found)
	    
	    his->AddBinContent(1,1); //  pulse counter for each pad
	    for (Int_t ipos=0; ipos<last-first; ipos++){
	      Int_t signal= (Int_t)(tempHis->GetBinContent(ipos+first)-baseline);
	      his->AddBinContent(ipos+5,signal);
	    }
	    printf("adding ...  %s: Signal %d at bin %d \n", hname, max-(Int_t)baseline, maxpos+fGateWidth);
	    his->Write(hname,kOverwrite);
	  }	


	  // 2D histogramm for pulse spread within a DDL (normalized to one)
	  char hname2d[100];
	  snprintf(hname2d,100,"2Dhisto_ddl%d",ddl);
	  TH2F *his2d = (TH2F*)fileOut.Get(hname2d);
	  if (!his2d ) { // new entry (ddl was not seen before)

	    his2d = new TH2F(hname2d,hname2d, fPulseLength, 0., (Double_t)fPulseLength, 50,-0.02,0.02);
	    for (Int_t ipos=0; ipos<last-first; ipos++){
	      Double_t signal = tempHis->GetBinContent(ipos+first)-baseline;
	      if (TMath::Abs(signal/maxCorr)>1e-10)  // zero bins are biased
		his2d->Fill(ipos,signal/maxCorr);
	    }
	    his2d->Write(hname2d);
	    printf("new  %s: \n", hname2d);
	  } else {  // adding pulse to existing histogram 

	    for (Int_t ipos=0; ipos<last-first; ipos++){
	      Double_t signal= tempHis->GetBinContent(ipos+first)-baseline;
	      if (TMath::Abs(signal/maxCorr)>1e-10) // zero bins are biased
		his2d->Fill(ipos,signal/maxCorr);
	    }
	    his2d->Write(hname2d,kOverwrite);
	  }
	  
	  tempHis->Reset();

        } // pulse stored

      } // bunch loop
    }// channel loop
  } // ddl loop
  
  tempHis->~TH1I();
  tempRMSHis->~TH1I();
  printf("Finished to read event ... \n");
  fileOut.Close();

}

//____________________________________________________________________________
void AliTPCCalibTCF::MergeHistoPerSector(const char *nameFileIn) {
  //
  // Merges all histograms within one sector, calculates the TCF parameters
  // of the 'histogram-per-sector' and stores (histo and parameters) into 
  // seperated files ...
  //
  // note: first 4 timebins of a histogram hold specific informations
  //       about number of collected pulses, sector, row and pad
  //
  // 'nameFileIn':  root file produced with Process function which holds
  //                one histogram per pad (sum of signals of proper pulses)
  // 'Sec+nameFileIn': root file with one histogram per sector
  //                   (information of row and pad are set to -1)
  //

  TFile fileIn(nameFileIn,"READ");
  TH1F *hisPad = 0;
  TKey *key = 0;
  TIter next( fileIn.GetListOfKeys() );

  char nameFileOut[100];
  snprintf(nameFileOut,100,"Sec-%s",nameFileIn);

  TFile fileOut(nameFileOut,"RECREATE");
  fileOut.cd();
  
  Int_t nHist = fileIn.GetNkeys();
  Int_t iHist = 0; // histogram counter for merge-status print
  
  while ( (key=(TKey*)next()) ) {

    iHist++;
    TString name(key->GetName());
    if (name.Contains("ddl") ) continue;  // ignore the 2d histogramms per ddl

    hisPad = (TH1F*)fileIn.Get(name.Data()); // copy object to memory

    Int_t pulseLength = hisPad->GetNbinsX() -4; 
    // -4 because first four timebins contain pad specific informations
    Int_t npulse = (Int_t)hisPad->GetBinContent(1);
    Int_t sector = (Int_t)hisPad->GetBinContent(2);
  
    char hname[100];
    snprintf(hname,100,"sector%d",sector);
    TH1F *his = (TH1F*)fileOut.Get(hname);
    
    if (!his ) { // new histogram (new sector)
      his = new TH1F(hname,hname, pulseLength+4, 0, pulseLength+4);
      his->SetBinContent(1,npulse); // pulse counter
      his->SetBinContent(2,sector); // set sector info 
      his->SetBinContent(3,-1); // set to dummy value 
      his->SetBinContent(4,-1); // set to dummy value
      for (Int_t ipos=0; ipos<pulseLength; ipos++){
	his->SetBinContent(ipos+5,hisPad->GetBinContent(ipos+5));
      }
      his->Write(hname);
      printf("found  %s ...\n", hname);
    } else { // add to existing histogram for sector
      his->AddBinContent(1,npulse); // pulse counter      
      for (Int_t ipos=0; ipos<pulseLength; ipos++){
	his->AddBinContent(ipos+5,hisPad->GetBinContent(ipos+5));
      }
      his->Write(hname,kOverwrite); 
    }

    if (iHist%500==0) {
      printf("merging status: \t %d pads out of %d \n",iHist, nHist);
    }
  }

  printf("merging done ...\n");
  fileIn.Close();
  fileOut.Close();


}


//____________________________________________________________________________
void AliTPCCalibTCF::AnalyzeRootFile(const char *nameFileIn, Int_t minNumPulse, Int_t histStart, Int_t histEnd) {
  //
  // This function takes a prepeared root file (accumulated histograms: output
  // of process function) and performs an analysis (fit and equalization) in 
  // order to get the TCF parameters. These are stored in an TNtuple along with 
  // the pad and creation infos. The tuple is written to the output file 
  // "TCFparam+nameFileIn"
  // To reduce the analysis time, the minimum number of accumulated pulses within 
  // one histogram 'minNumPulse' (to perform the analysis on) can be set
  //

  TFile fileIn(nameFileIn,"READ");
  TH1F *hisIn;
  TKey *key;
  TIter next( fileIn.GetListOfKeys() );

  char nameFileOut[100];
  snprintf(nameFileOut,100,"TCF-%s",nameFileIn);
  
  TFile fileOut(nameFileOut,"RECREATE");
  fileOut.cd();

  TNtuple *paramTuple = new TNtuple("TCFparam","TCFparameter","sec:row:pad:npulse:Z0:Z1:Z2:P0:P1:P2");
  
  Int_t nHist = fileIn.GetNkeys(); 
  Int_t iHist = 0;  // counter for print of analysis-status
  
  while ((key = (TKey *) next())) { // loop over histograms
    ++iHist;
    if(iHist < histStart || iHist  > histEnd) {continue;}

    TString name(key->GetName());
    if (name.Contains("ddl") ) continue;  // ignore the 2d histogramms per ddl

    hisIn = (TH1F*)fileIn.Get(key->GetName()); // copy object to memory
  
    Int_t numPulse = (Int_t)hisIn->GetBinContent(1); 
    if ( numPulse >= minNumPulse ) {
      printf("Analyze histogram %d out of %d\n",iHist,nHist);
      Double_t* coefP = new Double_t[3];
      Double_t* coefZ = new Double_t[3];
      for(Int_t i = 0; i < 3; i++){
	coefP[i] = 0;
	coefZ[i] = 0;
      }
      // perform the analysis on the given histogram 
      Int_t fitOk = AnalyzePulse(hisIn, coefZ, coefP);    
      if (fitOk) { // Add found parameters to file 
	Int_t sector = (Int_t)hisIn->GetBinContent(2);
	Int_t row = (Int_t)hisIn->GetBinContent(3);
	Int_t pad = (Int_t)hisIn->GetBinContent(4);
	paramTuple->Fill(sector,row,pad,numPulse,coefZ[0],coefZ[1],coefZ[2],coefP[0],coefP[1],coefP[2]);
      }
      coefP->~Double_t();
      coefZ->~Double_t();
    } else {
      printf("Skip histogram %d out of %d | not enough accumulated pulses\n",iHist,nHist);
    }
    
  }

  fileIn.Close();
  paramTuple->Write();
  fileOut.Close();

}


//____________________________________________________________________________
Int_t AliTPCCalibTCF::AnalyzePulse(TH1F * const hisIn, Double_t *coefZ, Double_t *coefP) {
  //
  // Performs the analysis on one specific pulse (histogram) by means of fitting
  // the pulse and equalization of the pulseheight. The found TCF parameters 
  // are stored in the arrays coefZ and coefP
  //

  Int_t pulseLength = hisIn->GetNbinsX() -4; 
  // -4 because the first four timebins usually contain pad specific informations
  Int_t npulse = (Int_t)hisIn->GetBinContent(1);
  Int_t sector = (Int_t)hisIn->GetBinContent(2);
  Int_t row = (Int_t)hisIn->GetBinContent(3);
  Int_t pad = (Int_t)hisIn->GetBinContent(4);
  
  // write pulseinformation to TNtuple and normalize to 100 ADC (because of 
  // given upper and lower fit parameter limits) in order to pass the pulse
  // to TMinuit

  TNtuple *dataTuple = new TNtuple("ntupleFit","Pulse","timebin:sigNorm:error");  
  Double_t error  = 0.05;
  Double_t max = hisIn->GetMaximum(FLT_MAX);
  for (Int_t ipos=0; ipos<pulseLength; ipos++) {
    Double_t errorz=error;
    if (ipos>100) { errorz = error*100; } // very simple weight: FIXME in case
    Double_t signal = hisIn->GetBinContent(ipos+5);
    Double_t signalNorm = signal/max*100; //pulseheight normaliz. to 100ADC
    dataTuple->Fill(ipos, signalNorm, errorz);
  }
   
  // Call fit function (TMinuit) to get the first 2 PZ Values for the 
  // Tail Cancelation Filter
  Int_t fitOk = FitPulse(dataTuple, coefZ, coefP);
 
  if (fitOk) {
    // calculates the 3rd set (remaining 2 PZ values) in order to restore the
    // original height of the pulse
    Int_t equOk = Equalization(dataTuple, coefZ, coefP);
    if (!equOk) {
      Error("FindFit", "Pulse equalisation procedure failed - pulse abandoned ");
      printf("in Sector %d | Row %d | Pad %d |", sector, row, pad);
      printf(" Npulses: %d \n\n", npulse);
      coefP[2] = 0; coefZ[2] = 0;
      dataTuple->~TNtuple();
      return 0;
    }  
    printf("Calculated TCF parameters for: \n");
    printf("Sector %d | Row %d | Pad %d |", sector, row, pad);
    printf(" Npulses: %d \n", npulse);
    for(Int_t i = 0; i < 3; i++){
      printf("P[%d] = %f     Z[%d] = %f \n",i,coefP[i],i,coefZ[i]);
      if (i==2) { printf("\n"); }
    }
    dataTuple->~TNtuple();
    return 1;
  } else { // fit did not converge
    Error("FindFit", "TCF fit not converged - pulse abandoned ");
    printf("in Sector %d | Row %d | Pad %d |", sector, row, pad);
    printf(" Npulses: %d \n\n", npulse);
    coefP[2] = 0; coefZ[2] = 0;
    dataTuple->~TNtuple();
    return 0;
  }
  
}



//____________________________________________________________________________
void AliTPCCalibTCF::TestTCFonRootFile(const char *nameFileIn, const char *nameFileTCF,  Int_t minNumPulse, Int_t plotFlag, Int_t lowKey, Int_t upKey)
{
  //
  // Performs quality parameters evaluation of the calculated TCF parameters in 
  // the file 'nameFileTCF' for every (accumulated) histogram within the 
  // prepeared root file 'nameFileIn'. 
  // The found quality parameters are stored in an TNtuple which will be saved
  // in a Root file 'Quality-*'. 
  // If the parameter for the given pulse (given pad) was not found, the pulse 
  // is rejected.
  //

  TFile fileIn(nameFileIn,"READ");

  Double_t* coefP = new Double_t[3];
  Double_t* coefZ = new Double_t[3];
  for(Int_t i = 0; i < 3; i++){
    coefP[i] = 0;
    coefZ[i] = 0;
  }

  char nameFileOut[100];
  snprintf(nameFileOut,100,"Quality_%s_AT_%s",nameFileTCF, nameFileIn);
  TFile fileOut(nameFileOut,"RECREATE");

  TNtuple *qualityTuple = new TNtuple("TCFquality","TCF quality Values","sec:row:pad:npulse:heightDev:areaRed:widthRed:undershot:maxUndershot");
 
  TH1F *hisIn;
  TKey *key;
  TIter next( fileIn.GetListOfKeys() );

  Int_t nHist = fileIn.GetNkeys();
  Int_t iHist = 0;
  
  for(Int_t i=0;i<lowKey-1;i++){++iHist; key = (TKey *) next();}
  while ((key = (TKey *) next())) { // loop over saved histograms
    
    //  loading pulse to memory;
    TString name(key->GetName());
    if (name.Contains("ddl") ) continue;  // ignore the 2d histogramms per ddl

    printf("validating pulse %d out of %d\n",++iHist,nHist);
    hisIn = (TH1F*)fileIn.Get(key->GetName()); 
 

    // find the correct TCF parameter according to the his infos (first 4 bins)
    Int_t nPulse = FindCorTCFparam(hisIn, nameFileTCF, coefZ, coefP); 
    if (nPulse>=minNumPulse) {  // doing the TCF quality analysis 
      Double_t *quVal = GetQualityOfTCF(hisIn,coefZ,coefP, plotFlag);
      Int_t sector = (Int_t)hisIn->GetBinContent(2);
      Int_t row = (Int_t)hisIn->GetBinContent(3);
      Int_t pad = (Int_t)hisIn->GetBinContent(4);      
      qualityTuple->Fill(sector,row,pad,nPulse,quVal[0],quVal[1],quVal[2],quVal[3],quVal[4],quVal[5]);
      quVal->~Double_t();
    }
    
    if (iHist>=upKey) {break;}
    
  }

  fileOut.cd();
  qualityTuple->Write();

  coefP->~Double_t();
  coefZ->~Double_t();

  fileOut.Close();
  fileIn.Close();

}



//_____________________________________________________________________________
void AliTPCCalibTCF::TestTCFonRawFile(const char *nameRawFile, const char *nameFileOut, const char *nameFileTCF, Int_t minNumPulse, Int_t plotFlag, bool bUseHLTOUT) {
  //
  // Performs quality parameters evaluation of the calculated TCF parameters in 
  // the file 'nameFileTCF' for every proper pulse (according to given thresholds)
  // within the RAW file 'nameRawFile'. 
  // The found quality parameters are stored in a TNtuple which will be saved
  // in the Root file 'nameFileOut'. If the parameter for the given pulse 
  // (given pad) was not found, the pulse is rejected.
  //

  //
  // Reads a RAW data file, extracts Pulses (according the given tresholds)
  // and test the found TCF parameters on them ...
  // 
  

  // create the data reader
  AliRawReader *rawReader = new AliRawReaderRoot(nameRawFile);
  if (!rawReader) {
    return;
  }

  // create HLT reader for redirection of TPC data from HLTOUT to TPC reconstruction
  AliRawReader *hltReader=AliRawHLTManager::AliRawHLTManager::CreateRawReaderHLT(rawReader, "TPC");

  // now choose the data source
  if (bUseHLTOUT) rawReader=hltReader;

  //  rawReader->Reset();
  rawReader->RewindEvents();

  if (!rawReader->NextEvent()) {
    printf("no events found in %s\n",nameRawFile);
    return;
  }

  Double_t* coefP = new Double_t[3];
  Double_t* coefZ = new Double_t[3];
  for(Int_t i = 0; i < 3; i++){
    coefP[i] = 0;
    coefZ[i] = 0;
  }

  Int_t ievent = 0;
  
  TH1I *tempHis = new TH1I("tempHis","tempHis",fSample+fGateWidth,fGateWidth,fSample+fGateWidth);
  TH1I *tempRMSHis = new TH1I("tempRMSHis","tempRMSHis",2000,0,2000);
  
  TFile fileOut(nameFileOut,"UPDATE"); // Quality Parameters storage
  TNtuple *qualityTuple = (TNtuple*)fileOut.Get("TCFquality");
  if (!qualityTuple) { // no entry in file
    qualityTuple = new TNtuple("TCFquality","TCF quality Values","sec:row:pad:npulse:heightDev:areaRed:widthRed:undershot:maxUndershot:pulseRMS");
  }

  do {

    printf("Reading next event ... Nr:%d\n",ievent);
    AliTPCRawStreamV3 *rawStream = new AliTPCRawStreamV3(rawReader);
    rawReader->Select("TPC");
    ievent++;

    while ( rawStream->NextDDL() ){
      while ( rawStream->NextChannel() ){
        
        const Int_t sector = rawStream->GetSector();
        const Int_t row    = rawStream->GetRow();
        const Int_t pad = rawStream->GetPad();
        
        while ( rawStream->NextBunch() ){
          UInt_t  startTbin    = rawStream->GetStartTimeBin();
          Int_t  bunchlength  = rawStream->GetBunchLength();
          const UShort_t *sig = rawStream->GetSignals();
          for (Int_t iTimeBin = 0; iTimeBin<bunchlength; iTimeBin++){
            const Int_t time = startTbin-iTimeBin;
            Float_t signal=(Float_t)sig[iTimeBin];

            // this pad always gave a useless signal, probably induced by the supply
            // voltage of the gate signal (date:2008-Aug-07)
            if(sector==51 && row==95 && pad==0) {
              continue;
            }

            // only process pulses of pads with correct address
            if(sector<0 || sector+1 > Int_t(AliTPCROC::Instance()->GetNSector())) {
              continue;
            }
            if(row<0 || row+1 > Int_t(AliTPCROC::Instance()->GetNRows(sector))) {
              continue;
            }
            if(pad<0 || pad+1 > Int_t(AliTPCROC::Instance()->GetNPads(sector,row))) {
              continue;
            }

            // still the same pad, save signal to temporary histogram
            if (time<=fSample+fGateWidth && time>fGateWidth) {
              tempHis->SetBinContent(time,signal);
            }
          }
        }
        
        Int_t max = (Int_t)tempHis->GetMaximum(FLT_MAX);
        Int_t maxpos =  tempHis->GetMaximumBin();

        Int_t first = (Int_t)TMath::Max(maxpos-10, 0);
        Int_t last  = TMath::Min((Int_t)maxpos+fPulseLength-10, fSample+fGateWidth);


        // simple baseline substraction ? better one needed ? (pedestalsubstr.?)
        // and RMS calculation with timebins before the pulse and at the end of
        // the signal
        for (Int_t ipos = 0; ipos<6; ipos++) {
          // before the pulse
          tempRMSHis->Fill(tempHis->GetBinContent(first+ipos));
        }
        for (Int_t ipos = 0; ipos<20; ipos++) {
          // at the end to get rid of pulses with serious baseline fluctuations
          tempRMSHis->Fill(tempHis->GetBinContent(last-ipos));
        }
        Double_t baseline = tempRMSHis->GetMean();
        Double_t rms = tempRMSHis->GetRMS();
        tempRMSHis->Reset();

        Double_t lowLim = fLowPulseLim+baseline;
        Double_t upLim = fUpPulseLim+baseline;

        // get rid of pulses which contain gate signal and/or too much noise
        // with the help of ratio of integrals
        Double_t intHist = 0;
        Double_t intPulse = 0;
        Double_t binValue;
        for(Int_t ipos=first; ipos<=last; ipos++) {
          binValue = TMath::Abs(tempHis->GetBinContent(ipos) - baseline);
          intHist += binValue;
          if(ipos>=first+5 && ipos<=first+15) {intPulse += binValue;}
        }

        // gets rid of high frequency noise:
        // calculating ratio (value one to the right of maximum)/(maximum)
        // has to be >= 0.1; if maximum==0 set ratio to 0.1
        Double_t maxCorr = max - baseline;
        Double_t binRatio = 0.1;
        if(TMath::Abs(maxCorr) > 1e-5 ) {
          binRatio = (tempHis->GetBinContent(maxpos+1) - baseline) / maxCorr;
        }

        // Decision if found pulse is a proper one according to given tresholds
        if (max>lowLim && max<upLim && !((last-first)<fPulseLength) && rms<fRMSLim && intHist/intPulse<fRatioIntLim && (binRatio >= 0.1) ){
          // note:
          // assuming that lowLim is higher than the pedestal value!
          char hname[100];
          snprintf(hname,100,"sec%drow%dpad%d",sector,row,pad);
          TH1F *his = new TH1F(hname,hname, fPulseLength+4, 0, fPulseLength+4);
          his->SetBinContent(1,1); //  pulse counter (1st pulse)
          his->SetBinContent(2,sector);  //  sector
          his->SetBinContent(3,row);  //  row
          his->SetBinContent(4,pad);  //  pad

          for (Int_t ipos=0; ipos<last-first; ipos++){
            const Double_t signal = tempHis->GetBinContent(ipos+first)-baseline;
            his->SetBinContent(ipos+5,signal);
          }

          printf("Pulse found in %s: ADC %d at bin %d \n", hname, max, maxpos+fGateWidth);

          // find the correct TCF parameter according to the his infos
          // (first 4 bins)
          Int_t nPulse = FindCorTCFparam(his, nameFileTCF, coefZ, coefP);

          if (nPulse>=minNumPulse) {  // Parameters found - doing the TCF quality analysis
	    Double_t *quVal = GetQualityOfTCF(his,coefZ,coefP, plotFlag);
            qualityTuple->Fill(sector,row,pad,nPulse,quVal[0],quVal[1],quVal[2],quVal[3],quVal[4],quVal[5]);
            quVal->~Double_t();
          }
          his->~TH1F();
        }
        tempHis->Reset();
      }
    }
  
    
    printf("Finished to read event ... \n");   

    delete rawStream;


  } while (rawReader->NextEvent()); // event loop

  printf("Finished to read file - close output file ... \n");
  
  fileOut.cd();
  qualityTuple->Write("TCFquality",kOverwrite);
  fileOut.Close();
  
  tempHis->~TH1I();
  tempRMSHis->~TH1I();

  coefP->~Double_t();
  coefZ->~Double_t();

  rawReader->~AliRawReader();
  
}

//____________________________________________________________________________
TH2F *AliTPCCalibTCF::PlotOccupSummary2Dhist(const char *nameFileIn, Int_t side) {
  //
  // Plots the number of summed pulses per pad on a given TPC side
  // 'nameFileIn': root-file created with the Process function
  //

  TFile fileIn(nameFileIn,"READ");
  TH1F *his;
  TKey *key;
  TIter next(fileIn.GetListOfKeys());

  TH2F * his2D = new TH2F("his2D","his2D", 250,-250,250,250,-250,250);

  AliTPCROC * roc  = AliTPCROC::Instance();

  Int_t nHist=fileIn.GetNkeys();
  if (!nHist) { return 0; }

  Int_t iHist = 0;
  Float_t xyz[3];

  Int_t binx = 0;
  Int_t biny = 0;

  Int_t npulse = 0;
  Int_t sec = 0;
  Int_t row = 0;
  Int_t pad = 0;

  while ((key = (TKey *) next())) { // loop over histograms within the file
    iHist++;
    
    TString name(key->GetName());
    if (name.Contains("ddl") ) continue;  // ignore the 2d histogramms per ddl

    his = (TH1F*)fileIn.Get(key->GetName()); // copy object to memory

    npulse = (Int_t)his->GetBinContent(1);
    sec = (Int_t)his->GetBinContent(2);
    row = (Int_t)his->GetBinContent(3);
    pad = (Int_t)his->GetBinContent(4);

    if ( (side==0) && (sec%36>=18) ) continue;
    if ( (side>0) && (sec%36<18) ) continue;

    if ( (row<0) && (pad<0) ) { // row and pad are equal to -1, then -> summed pulses per sector
      // fill all pad with this values
      for (UInt_t rowi=0; rowi<roc->GetNRows(sec); rowi++) {
        for (UInt_t padi=0; padi<roc->GetNPads(sec,rowi); padi++) {
          roc->GetPositionGlobal(sec,rowi,padi,xyz);
          binx = 1+TMath::Nint((xyz[0]+250.)*0.5);
          biny = 1+TMath::Nint((xyz[1]+250.)*0.5);
          his2D->SetBinContent(binx,biny,npulse);
        }
      }
    } else {
      roc->GetPositionGlobal(sec,row,pad,xyz);
      binx = 1+TMath::Nint((xyz[0]+250.)*0.5);
      biny = 1+TMath::Nint((xyz[1]+250.)*0.5);

      his2D->SetBinContent(binx,biny,npulse);
    }
    if (iHist%100==0){ printf("hist %d out of %d\n",iHist,nHist);}
  }
  his2D->SetXTitle("x (cm)");
  his2D->SetYTitle("y (cm)");
  his2D->SetStats(0);

  his2D->DrawCopy("colz");

  if (!side) {
    gPad->SetTitle("A side");
  } else {
    gPad->SetTitle("C side");
  }

  return his2D;
}


//____________________________________________________________________________
void AliTPCCalibTCF::PlotOccupSummary(const char *nameFile, Int_t side, Int_t nPulseMin) {
  //
  // Plots the number of summed pulses per pad above a given minimum at the 
  // pad position at a given TPC side
  // 'nameFile': root-file created with the Process function
  //

  TFile *file = new TFile(nameFile,"READ");
  TH1F *his;
  TKey *key;
  TIter next( file->GetListOfKeys() );


  char nameFileOut[100];
  snprintf(nameFileOut,100,"Occup-%s",nameFile);
  TFile fileOut(nameFileOut,"RECREATE");
  // fileOut.cd();

  TNtuple *ntuple = new TNtuple("ntuple","ntuple","x:y:z:npulse");
  // ntuple->SetDirectory(0); // force to be memory resistent

  Int_t nHist=file->GetNkeys();
  if (!nHist) { return; }
  Int_t iHist = 0;

  Int_t secWise = 0;

  while ((key = (TKey *) next())) { // loop over histograms within the file
    
    TString name(key->GetName());
    if (name.Contains("ddl") ) continue;  // ignore the 2d histogramms per ddl

    his = (TH1F*)file->Get(key->GetName()); // copy object to memory
    iHist++;
    Int_t npulse = (Int_t)his->GetBinContent(1);
    Int_t sec = (Int_t)his->GetBinContent(2);
    Int_t row = (Int_t)his->GetBinContent(3);
    Int_t pad = (Int_t)his->GetBinContent(4);

    if ( (row<0) && (pad<0) ) { // row and pad are equal to -1, then -> summed pulses per sector
      row = 40; pad = 40;    // set to approx middle row for better plot
      secWise=1;
    }

    Float_t *pos = new Float_t[3];
    // find x,y,z position of the pad
    AliTPCROC::Instance()->GetPositionGlobal(sec,row,pad,pos); 
    if (npulse>=nPulseMin) { 
      ntuple->Fill(pos[0],pos[1],pos[2],npulse);
      if (iHist%100==0){ printf("hist %d out of %d\n",iHist,nHist);}
    }
    pos->~Float_t();
  }

  if (secWise) { // pulse per sector
    ntuple->SetMarkerStyle(8);
    ntuple->SetMarkerSize(4);
  } else {        // pulse per Pad
    ntuple->SetMarkerStyle(7);
  }

  char cSel[100];
  if (!side) {
    snprintf(cSel,100,"z>0&&npulse>=%d",nPulseMin);
    ntuple->Draw("y:x:npulse",cSel,"colz");
  } else {
    snprintf(cSel,100,"z<0&&npulse>=%d",nPulseMin);
    ntuple->Draw("y:x:npulse",cSel,"colz");
  }

  if (!side) {
    gPad->SetTitle("A side");
  } else {
    gPad->SetTitle("C side");
  }


  ntuple->Write();
  fileOut.Close();
  file->Close();
}

//____________________________________________________________________________
void AliTPCCalibTCF::PlotQualitySummary(const char *nameFileQuality, const char *plotSpec, const char *cut, const char *pOpt)
{
  // 
  // This function is an easy interface to load the QualityTuple (produced with
  // the function 'TestOn%File' and plots them according to the plot specifications
  // 'plotSpec' e.g. "widthRed:maxUndershot"
  // One may also set cut and plot options ("cut","pOpt") 
  //
  // The stored quality parameters are ...
  //   sec:row:pad:npulse: ... usual pad info
  //   heightDev ... height deviation in percent
  //   areaRed ... area reduction in percent
  //   widthRed ... width reduction in percent
  //   undershot ... mean undershot after the pulse in ADC
  //   maxUndershot ... maximum of the undershot after the pulse in ADC
  //   pulseRMS ... RMS of the pulse used to calculate the Quality parameters in ADC
  //

  TFile file(nameFileQuality,"READ");
  TNtuple *qualityTuple = (TNtuple*)file.Get("TCFquality");
  //gStyle->SetPalette(1);
  
  TH2F *his2D = new TH2F(plotSpec,nameFileQuality,11,-10,1,25,1,100);
  char plSpec[100];
  snprintf(plSpec,100,"%s>>%s",plotSpec,plotSpec);
  qualityTuple->Draw(plSpec,cut,pOpt);

  gStyle->SetLabelSize(0.03,"X");
  gStyle->SetLabelSize(0.03,"Y");
  gStyle->SetLabelSize(0.03,"Z");
  gStyle->SetLabelOffset(-0.02,"X");
  gStyle->SetLabelOffset(-0.01,"Y");
  gStyle->SetLabelOffset(-0.03,"Z");

  his2D->GetXaxis()->SetTitle("max. undershot [ADC]");
  his2D->GetYaxis()->SetTitle("width Reduction [%]");

  his2D->DrawCopy(pOpt);

  gPad->SetPhi(0.1);gPad->SetTheta(90);
  
  his2D->~TH2F();
  
}

//_____________________________________________________________________________
Int_t AliTPCCalibTCF::FitPulse(TNtuple *dataTuple, Double_t *coefZ, Double_t *coefP) {
  //
  // function to fit one pulse and to calculate the according pole-zero parameters
  //
 
  // initialize TMinuit with a maximum of 8 params
  TMinuit *minuitFit = new TMinuit(8);
  minuitFit->mncler();                    // Reset Minuit's list of paramters
  minuitFit->SetPrintLevel(-1);           // No Printout
  minuitFit->SetFCN(AliTPCCalibTCF::FitFcn); // To set the address of the 
                                           // minimization function  
  minuitFit->SetObjectFit(dataTuple);
  
  Double_t arglist[10];
  Int_t ierflg = 0;
  
  arglist[0] = 1;
  minuitFit->mnexcm("SET ERR", arglist ,1,ierflg);
  
  // Set standard starting values and step sizes for each parameter
  // upper and lower limit (in a reasonable range) are set to improve 
  // the stability of TMinuit
  static Double_t vstart[8] = {125, 4.0, 0.3, 0.5, 5.5, 100,    1, 2.24};
  static Double_t step[8]   = {0.1, 0.1,  0.1, 0.1, 0.1, 0.1,  0.1,  0.1};
  static Double_t min[8]    = {100,  3.,  0.1, 0.2,  3.,  60.,  0.,  2.0};
  static Double_t max[8]    = {200, 20.,   5.,  3., 30., 300., 20., 2.5};
  
  minuitFit->mnparm(0, "A1", vstart[0], step[0], min[0], max[0], ierflg);
  minuitFit->mnparm(1, "A2", vstart[1], step[1], min[1], max[1], ierflg);
  minuitFit->mnparm(2, "A3", vstart[2], step[2], min[2], max[2], ierflg);
  minuitFit->mnparm(3, "T1", vstart[3], step[3], min[3], max[3], ierflg);
  minuitFit->mnparm(4, "T2", vstart[4], step[4], min[4], max[4], ierflg);
  minuitFit->mnparm(5, "T3", vstart[5], step[5], min[5], max[5], ierflg);
  minuitFit->mnparm(6, "T0", vstart[6], step[6], min[6], max[6], ierflg);
  minuitFit->mnparm(7, "TTP", vstart[7], step[7], min[7], max[7],ierflg);
  minuitFit->FixParameter(7); // 2.24 ... out of pulserRun Fit (->IRF)

  // Now ready for minimization step
  arglist[0] = 2000;   // max num of iterations
  arglist[1] = 0.1;    // tolerance

  minuitFit->mnexcm("MIGRAD", arglist ,2,ierflg);
  
  Double_t p1 = 0.0 ;
  minuitFit->mnexcm("SET NOW", &p1 , 0, ierflg) ;  // No Warnings
  
  if (ierflg == 4) { // Fit failed
    for (Int_t i=0;i<3;i++) { 
      coefP[i] = 0; 
      coefZ[i] = 0; 
    }
    minuitFit->~TMinuit();
    return 0;
  } else { // Fit successfull

    // Extract parameters from TMinuit
    Double_t *fitParam = new Double_t[6];
    for (Int_t i=0;i<6;i++) {
      Double_t err = 0;
      Double_t val = 0;
      minuitFit->GetParameter(i,val,err);
      fitParam[i] = val;
    } 
    
    // calculates the first 2 sets (4 PZ values) out of the fitted parameters
    Double_t *valuePZ = ExtractPZValues(fitParam);
   
    // TCF coefficients which are used for the equalisation step (stage)
    // ZERO/POLE Filter
    coefZ[0] = TMath::Exp(-1/valuePZ[2]);
    coefZ[1] = TMath::Exp(-1/valuePZ[3]);
    coefP[0] = TMath::Exp(-1/valuePZ[0]);
    coefP[1] = TMath::Exp(-1/valuePZ[1]);
   
    fitParam->~Double_t();
    valuePZ->~Double_t();
    minuitFit->~TMinuit();

    return 1;

  }

}


//____________________________________________________________________________
void AliTPCCalibTCF::FitFcn(Int_t &/*nPar*/, Double_t */*grad*/, Double_t &f, Double_t * const par, Int_t /*iflag*/)
{
  //
  // Minimization function needed for TMinuit with FitFunction included 
  // Fit function: Sum of three convolution terms (IRF conv. with Exp.)
  //

  // Get Data ...
  TNtuple *dataTuple = (TNtuple *) gMinuit->GetObjectFit();

  //calculate chisquare
  Double_t chisq = 0;
  Double_t delta = 0;
  for (Int_t i=0; i<dataTuple->GetEntries(); i++) { // loop over data points
    dataTuple->GetEntry(i);
    Float_t *p = dataTuple->GetArgs();
    Double_t t = p[0];
    Double_t signal = p[1];   // Normalized signal
    Double_t error = p[2]; 

    // definition and evaluation if the IonTail specific fit function
    Double_t sigFit = 0;
    
    Double_t ttp = par[7];   // signal shaper raising time
    t=t-par[6];              // time adjustment
    
    if (t<0) {
      sigFit = 0;
    } else {
      Double_t f1 = 1/TMath::Power((4-ttp/par[3]),5)*(24*ttp*TMath::Exp(4)*(TMath::Exp(-t/par[3]) - TMath::Exp(-4*t/ttp) * ( 1+t*(4-ttp/par[3])/ttp+TMath::Power(t*(4-ttp/par[3])/ttp,2)/2 + TMath::Power(t*(4-ttp/par[3])/ttp,3)/6 + TMath::Power(t*(4-ttp/par[3])/ttp,4)/24)));
      
      Double_t f2 = 1/TMath::Power((4-ttp/par[4]),5)*(24*ttp*TMath::Exp(4)*(TMath::Exp(-t/par[4]) - TMath::Exp(-4*t/ttp) * ( 1+t*(4-ttp/par[4])/ttp+TMath::Power(t*(4-ttp/par[4])/ttp,2)/2 + TMath::Power(t*(4-ttp/par[4])/ttp,3)/6 + TMath::Power(t*(4-ttp/par[4])/ttp,4)/24)));
      
      Double_t f3 = 1/TMath::Power((4-ttp/par[5]),5)*(24*ttp*TMath::Exp(4)*(TMath::Exp(-t/par[5]) - TMath::Exp(-4*t/ttp) * ( 1+t*(4-ttp/par[5])/ttp+TMath::Power(t*(4-ttp/par[5])/ttp,2)/2 + TMath::Power(t*(4-ttp/par[5])/ttp,3)/6 + TMath::Power(t*(4-ttp/par[5])/ttp,4)/24)));
      
      sigFit = par[0]*f1 + par[1]*f2 +par[2]*f3;
    }

    // chisqu calculation
    delta  = (signal-sigFit)/error;
    chisq += delta*delta;
  }

  f = chisq;

}



//____________________________________________________________________________
Double_t* AliTPCCalibTCF::ExtractPZValues(Double_t *param) {
  //
  // Calculation of Pole and Zero values out of fit parameters
  //

  Double_t vA1, vA2, vA3, vTT1, vTT2, vTT3, vTa, vTb;
  vA1 = 0;  vA2 = 0;  vA3 = 0;
  vTT1 = 0; vTT2 = 0; vTT3 = 0;
  vTa = 0; vTb = 0;
  
  // nasty method of sorting the fit parameters to avoid wrong mapping
  // to the different stages of the TCF filter
  // (e.g. first 2 fit parameters represent the electron signal itself!)

  if ((param[3]-param[4]) <1e-5 ) {param[3]=param[3]+0.0001;} // if equal
  if ((param[5]-param[4]) <1e-5 ) {param[5]=param[5]+0.0001;} // if equal
  
  if ((param[5]>param[4])&&(param[5]>param[3])) {
    if (param[4]>=param[3]) {
      vA1 = param[0];  vA2 = param[1];  vA3 = param[2];
      vTT1 = param[3]; vTT2 = param[4]; vTT3 = param[5];
    } else {
      vA1 = param[1];  vA2 = param[0];  vA3 = param[2];
      vTT1 = param[4]; vTT2 = param[3]; vTT3 = param[5];
    }
  } else if ((param[4]>param[5])&&(param[4]>param[3])) {
    if (param[5]>=param[3]) {
      vA1 = param[0];  vA2 = param[2];  vA3 = param[1];
      vTT1 = param[3]; vTT2 = param[5]; vTT3 = param[4];
    } else {
      vA1 = param[2];  vA2 = param[0];  vA3 = param[1];
      vTT1 = param[5]; vTT2 = param[3]; vTT3 = param[4];
    }
  } else if ((param[3]>param[4])&&(param[3]>param[5])) {
    if (param[5]>=param[4]) {
      vA1 = param[1];  vA2 = param[2];  vA3 = param[0];
      vTT1 = param[4]; vTT2 = param[5]; vTT3 = param[3];
    } else {
      vA1 = param[2];  vA2 = param[1];  vA3 = param[0];
      vTT1 = param[5]; vTT2 = param[4]; vTT3 = param[3];
    }    
  }
  

  // Transformation of fit parameters into PZ values (needed by TCF) 
  Double_t beq = (vA1/vTT2+vA1/vTT3+vA2/vTT1+vA2/vTT3+vA3/vTT1+vA3/vTT2)/(vA1+vA2+vA3);
  Double_t ceq = (vA1/(vTT2*vTT3)+vA2/(vTT1*vTT3)+vA3/(vTT1*vTT2))/(vA1+vA2+vA3);
  
  Double_t  s1 = -beq/2-sqrt((beq*beq-4*ceq)/4);
  Double_t  s2 = -beq/2+sqrt((beq*beq-4*ceq)/4);
  
  if (vTT2<vTT3) {// not necessary but avoids significant undershots in first PZ 
    vTa = -1/s1;
    vTb = -1/s2;
  }else{ 
    vTa = -1/s2;
    vTb = -1/s1;
  }
    
  Double_t *valuePZ = new Double_t[4];
  valuePZ[0]=vTa;
  valuePZ[1]=vTb;
  valuePZ[2]=vTT2;
  valuePZ[3]=vTT3;
      
  return valuePZ;
  
}


//____________________________________________________________________________
Int_t AliTPCCalibTCF::Equalization(TNtuple *dataTuple, Double_t *coefZ, Double_t *coefP) {
  //
  // calculates the 3rd set of TCF parameters (remaining 2 PZ values) in 
  // order to restore the original pulse height and adds them to the passed arrays
  //

  const Int_t kPulseLength = dataTuple->GetEntries();

  if (kPulseLength<2) {
    //    prinft("PulseLength does not make sense\n");
    return 0;
  }

  Double_t *s0 = new Double_t[kPulseLength]; // original pulse
  Double_t *s1 = new Double_t[kPulseLength]; // pulse after 1st PZ filter
  Double_t *s2 = new Double_t[kPulseLength]; // pulse after 2nd PZ filter

  for (Int_t ipos=0; ipos<kPulseLength; ipos++) {
    dataTuple->GetEntry(ipos);
    Float_t *p = dataTuple->GetArgs();
    s0[ipos] = p[1]; 
  }
  
  // non-discret implementation of the first two TCF stages (recursive formula)
  // discrete Altro emulator is not used because of accuracy!
  s1[0] = s0[0]; // 1st PZ filter
  for(Int_t ipos = 1; ipos < kPulseLength ; ipos++){
    s1[ipos] = s0[ipos] + coefP[0]*s1[ipos-1] - coefZ[0]*s0[ipos-1];
  }
  s2[0] = s1[0]; // 2nd PZ filter
  for(Int_t ipos = 1; ipos < kPulseLength ; ipos++){
    s2[ipos] = s1[ipos] + coefP[1]*s2[ipos-1] - coefZ[1]*s1[ipos-1];
  }
  
  // find maximum amplitude and position of original pulse and pulse after 
  // the first two stages of the TCF 
  Int_t s0pos = 0; 
  Double_t s0ampl = s0[0], s2ampl = s2[0]; // start values
  for(Int_t ipos = 1; ipos < kPulseLength; ipos++){
    if (s0[ipos] > s0ampl){
      s0ampl = s0[ipos]; 
      s0pos = ipos;      // should be pos 11 ... check?
    }
    if (s2[ipos] > s2ampl){
      s2ampl = s2[ipos];
    }    
  }
  // calculation of 3rd set ...
  if(s0ampl > s2ampl){
    coefZ[2] = 0;
    coefP[2] = (s0ampl - s2ampl)/s0[s0pos-1];
  } else if (s0ampl < s2ampl) {
    coefP[2] = 0;
    coefZ[2] = (s2ampl - s0ampl)/s0[s0pos-1];
  } else { // same height ? will most likely not happen ?
    printf("No equalization because of identical height\n");
    coefP[2] = 0;
    coefZ[2] = 0;
  }

  delete [] s0;
  delete [] s1;
  delete [] s2;
  
  // if equalization out of range (<0 or >=1) it failed!
  // if ratio of amplitudes of fittet to original pulse < 0.9 it failed!
  if (coefP[2]<0 || coefZ[2]<0 || coefP[2]>=1 || coefZ[2]>=1 || TMath::Abs(s2ampl / s0ampl)<0.9) {
    return 0; 
  } else {
    return 1;
  }
  
}



//____________________________________________________________________________
Int_t AliTPCCalibTCF::FindCorTCFparam(TH1F * const hisIn, const char *nameFileTCF, Double_t *coefZ, Double_t *coefP) {
  //
  // This function searches for the correct TCF parameters to the given
  // histogram 'hisIn' within the file 'nameFileTCF' 
  // If no parameters for this pad (padinfo within the histogram!) where found
  // the function returns 0

  //  Int_t numPulse = (Int_t)hisIn->GetBinContent(1); // number of pulses
  Int_t sector = (Int_t)hisIn->GetBinContent(2);
  Int_t row = (Int_t)hisIn->GetBinContent(3);
  Int_t pad = (Int_t)hisIn->GetBinContent(4);
  Int_t nPulse = 0; 

  //-- searching for calculated TCF parameters for this pad/sector
  TFile fileTCF(nameFileTCF,"READ");
  TNtuple *paramTuple = (TNtuple*)fileTCF.Get("TCFparam");

  // create selection criteria to find the correct TCF params
  char sel[100];   
  if ( paramTuple->GetEntries("row==-1&&pad==-1") ) { 
    // parameters per SECTOR
    snprintf(sel,100,"sec==%d&&row==-1&&pad==-1",sector);
  } else {            
    // parameters per PAD
    snprintf(sel,100,"sec==%d&&row==%d&&pad==%d",sector,row,pad);
  }

  // list should contain just ONE entry! ... otherwise there is a mistake!
  Long64_t entry = paramTuple->Draw(">>list",sel,"entrylist");
  TEntryList *list = (TEntryList*)gDirectory->Get("list");
  
  if (entry) { // TCF set was found for this pad
    Long64_t pos = list->GetEntry(0);
    paramTuple->GetEntry(pos);   // get specific TCF parameters       
    Float_t *p = paramTuple->GetArgs();
    // check ...
    if((sector-p[0])<1e-5) {printf("sector ok ... "); }          
    if((row-p[1])<1e-5) {printf("row ok ... "); }          
    if((pad-p[2])<1e-5) {printf("pad ok ... \n"); }          
    
    // number of averaged pulses used to produce TCF params
    nPulse = (Int_t)p[3]; 
    // TCF parameters
    coefZ[0] = p[4];  coefP[0] = p[7];
    coefZ[1] = p[5];  coefP[1] = p[8];
    coefZ[2] = p[6];  coefP[2] = p[9];
      
  } else { // no specific TCF parameters found for this pad 
    
    printf("  no specific TCF paramaters found for pad in ...\n");
    printf("  Sector %d | Row %d | Pad %d |\n", sector, row, pad);
    nPulse = 0;
    coefZ[0] = 0;  coefP[0] = 0;
    coefZ[1] = 0;  coefP[1] = 0;
    coefZ[2] = 0;  coefP[2] = 0;

  }

  fileTCF.Close();

  return nPulse; // number of averaged pulses for producing the TCF params
  
}


//____________________________________________________________________________
Double_t *AliTPCCalibTCF::GetQualityOfTCF(TH1F *hisIn, Double_t *coefZ, Double_t *coefP, Int_t plotFlag) {
  //
  // This function evaluates the quality parameters of the given TCF parameters
  // tested on the passed pulse (hisIn)
  // The quality parameters are stored in an array. They are ...
  //    height deviation [ADC]
  //    area reduction [percent]
  //    width reduction [percent]
  //    mean undershot [ADC]
  //    maximum of undershot after pulse [ADC]
  //    Pulse RMS [ADC]

  // perform ALTRO emulator
  TNtuple *pulseTuple = ApplyTCFilter(hisIn, coefZ, coefP, plotFlag); 

  printf("calculate quality val. for pulse in ... ");
  printf(" Sector %d | Row %d | Pad %d |\n", (Int_t)hisIn->GetBinContent(2),  (Int_t)hisIn->GetBinContent(3), (Int_t)hisIn->GetBinContent(4));
  
  // Reasonable limit for the calculation of the quality values
  Int_t binLimit = 80; 
  
  // ============== Variable preparation

  // -- height difference in percent of orginal pulse
  Double_t maxSig = pulseTuple->GetMaximum("sig");
  Double_t maxSigTCF = pulseTuple->GetMaximum("sigAfterTCF");      
  // -- area reduction (above zero!)
  Double_t area = 0;
  Double_t areaTCF = 0;    
  // -- width reduction at certain ADC treshold
  // TODO: set treshold at ZS treshold? (3 sigmas of noise?)
  Int_t threshold = 3; // treshold in percent
  Int_t threshADC = (Int_t)(maxSig/100*threshold);  
  Int_t startOfPulse = 0;   Int_t startOfPulseTCF = 0;
  Int_t posOfStart = 0;     Int_t posOfStartTCF = 0;
  Int_t widthFound = 0;     Int_t widthFoundTCF = 0;
  Int_t width = 0;          Int_t widthTCF = 0;
  // -- Calcluation of Undershot (mean of negavive signal after the first 
  // undershot)
  Double_t undershotTCF = 0;  
  Double_t undershotStart = 0;
  // -- Calcluation of Undershot (Sum of negative signal after the pulse)
  Double_t maxUndershot = 0;


  // === loop over timebins to calculate quality parameters
  for (Int_t i=0; i<binLimit; i++) {
   
    // Read signal values
    pulseTuple->GetEntry(i); 
    Float_t *p = pulseTuple->GetArgs();
    Double_t sig = p[1]; 
    Double_t sigTCF = p[2];

    // calculation of area (above zero)
    if (sig>0) {area += sig; }
    if (sigTCF>0) {areaTCF += sigTCF; }
    

    // Search for width at certain ADC treshold 
    // -- original signal
    if (widthFound == 0) {
      if( (sig > threshADC) && (startOfPulse == 0) ){
	startOfPulse = 1;
	posOfStart = i;
      }
      if( (sig <= threshADC) && (startOfPulse == 1) ){
	widthFound = 1;
	width = i - posOfStart + 1;	
      }
    }
    // -- signal after TCF
    if (widthFoundTCF == 0) {
      if( (sigTCF > threshADC) && (startOfPulseTCF == 0) ){
	startOfPulseTCF = 1;
	posOfStartTCF = i;
      }
      if( (sigTCF <= threshADC) && (startOfPulseTCF == 1) ){
	widthFoundTCF = 1;
	widthTCF = i -posOfStartTCF + 1;
      }
      
    }
      
    // finds undershot start
    if  ( (widthFoundTCF==1) && (sigTCF<0) ) {
      undershotStart = 1;
    }

    // Calculation of undershot sum (after pulse)
    if ( widthFoundTCF==1 ) {
      undershotTCF += sigTCF; 
    }

    // Search for maximal undershot (is equal to minimum after the pulse)
    if ( ((undershotStart-1)<1e-7)&&(i<(posOfStartTCF+widthTCF+20)) ) {
      if (maxUndershot>sigTCF) { maxUndershot = sigTCF; }
    }

  }  

  // ==  Calculation of Quality parameters

  // -- height difference in ADC
  Double_t heightDev = maxSigTCF-maxSig; 

  // Area reduction of the pulse in percent
  Double_t areaReduct = 100-areaTCF/area*100; 

  // Width reduction in percent
  Double_t widthReduct = 0;
  if ((widthFound==1)&&(widthFoundTCF==1)) { // in case of not too big IonTail 
    widthReduct = 100-(Double_t)widthTCF/(Double_t)width*100; 
    if (widthReduct<0) { widthReduct = 0;}  
  }

  // Undershot - mean of neg.signals after pulse
  Double_t length = 1;
  if (binLimit-widthTCF-posOfStartTCF) { length = (binLimit-widthTCF-posOfStartTCF);}
  Double_t undershot = undershotTCF/length; 


  // calculation of pulse RMS with timebins before and at the end of the pulse
  TH1I *tempRMSHis = new TH1I("tempRMSHis","tempRMSHis",100,-50,50);
  for (Int_t ipos = 0; ipos<6; ipos++) {
    // before the pulse
    tempRMSHis->Fill(hisIn->GetBinContent(ipos+5));
    // at the end
    tempRMSHis->Fill(hisIn->GetBinContent(hisIn->GetNbinsX()-ipos));
  }
  Double_t pulseRMS = tempRMSHis->GetRMS();
  tempRMSHis->~TH1I();
  
  if (plotFlag) {
    // == Output 
    printf("height deviation [ADC]:\t\t\t %3.1f\n", heightDev);
    printf("area reduction [percent]:\t\t %3.1f\n", areaReduct);
    printf("width reduction [percent]:\t\t %3.1f\n", widthReduct);
    printf("mean undershot [ADC]:\t\t\t %3.1f\n", undershot);
    printf("maximum of undershot after pulse [ADC]: %3.1f\n", maxUndershot);
    printf("RMS of the original (or summed) pulse [ADC]: \t %3.2f\n\n", pulseRMS);

  }

  Double_t *qualityParam = new Double_t[6];
  qualityParam[0] = heightDev;
  qualityParam[1] = areaReduct;
  qualityParam[2] = widthReduct;
  qualityParam[3] = undershot;
  qualityParam[4] = maxUndershot;
  qualityParam[5] = pulseRMS;

  pulseTuple->~TNtuple();

  return qualityParam;
}


//____________________________________________________________________________
TNtuple *AliTPCCalibTCF::ApplyTCFilter(TH1F * const hisIn, Double_t * const coefZ, Double_t * const coefP, Int_t plotFlag) {
  //
  // Applies the given TCF parameters on the given pulse via the ALTRO emulator 
  // class (discret values) and stores both pulses into a returned TNtuple
  //

  Int_t nbins = hisIn->GetNbinsX() -4; 
  // -1 because the first four timebins usually contain pad specific informations  
  Int_t nPulse = (Int_t)hisIn->GetBinContent(1); // Number of summed pulses
  Int_t sector = (Int_t)hisIn->GetBinContent(2);
  Int_t row = (Int_t)hisIn->GetBinContent(3);
  Int_t pad = (Int_t)hisIn->GetBinContent(4);
 
  // redirect histogram values to arrays (discrete for altro emulator)
  Double_t *signalIn = new Double_t[nbins];
  Double_t *signalOut = new Double_t[nbins];
  short *signalInD = new short[nbins]; 
  short *signalOutD = new short[nbins];
  for (Int_t ipos=0;ipos<nbins;ipos++) {
    Double_t signal = hisIn->GetBinContent(ipos+5); // summed signal
    signalIn[ipos]=signal/nPulse;                 // mean signal
    signalInD[ipos]=(short)(TMath::Nint(signalIn[ipos])); //discrete mean signal 
    signalOutD[ipos]=signalInD[ipos];    // will be overwritten by AltroEmulator    
  }

  // transform TCF parameters into ALTRO readable format (Integer)
  Int_t valK[3];
  Int_t valL[3];
  for (Int_t i=0; i<3; i++) {
    valK[i] = (Int_t)(coefP[i]*(TMath::Power(2,16)-1));
    valL[i] = (Int_t)(coefZ[i]*(TMath::Power(2,16)-1));
  }
    
  // discret ALTRO EMULATOR ____________________________
  AliTPCAltroEmulator *altro = new AliTPCAltroEmulator(nbins, signalOutD);
  altro->ConfigAltro(0,1,0,0,0,0); // perform just the TailCancelation
  altro->ConfigTailCancellationFilter(valK[0],valK[1],valK[2],valL[0],valL[1],valL[2]);
  altro->RunEmulation();
  delete altro;
  
  // non-discret implementation of the (recursive formula)
  // discrete Altro emulator is not used because of accuracy!
  Double_t *s1 = new Double_t[1000]; // pulse after 1st PZ filter
  Double_t *s2 = new Double_t[1000]; // pulse after 2nd PZ filter
  s1[0] = signalIn[0]; // 1st PZ filter
  for(Int_t ipos = 1; ipos<nbins; ipos++){
    s1[ipos] = signalIn[ipos] + coefP[0]*s1[ipos-1] - coefZ[0]*signalIn[ipos-1];
  }
  s2[0] = s1[0]; // 2nd PZ filter
  for(Int_t ipos = 1; ipos<nbins; ipos++){
    s2[ipos] = s1[ipos] + coefP[1]*s2[ipos-1] - coefZ[1]*s1[ipos-1];
  }
  signalOut[0] = s2[0]; // 3rd PZ filter
  for(Int_t ipos = 1; ipos<nbins; ipos++){
    signalOut[ipos] = s2[ipos] + coefP[2]*signalOut[ipos-1] - coefZ[2]*s2[ipos-1];
  }
  s1->~Double_t();
  s2->~Double_t();

  // writing pulses to tuple
  TNtuple *pulseTuple = new TNtuple("ntupleTCF","PulseTCF","timebin:sig:sigAfterTCF:sigND:sigNDAfterTCF");
  for (Int_t ipos=0;ipos<nbins;ipos++) {
    pulseTuple->Fill(ipos,signalInD[ipos],signalOutD[ipos],signalIn[ipos],signalOut[ipos]);
  }

  if (plotFlag) {
    char hname[100];
    snprintf(hname,100,"sec%drow%dpad%d",sector,row,pad);
    new TCanvas(hname,hname,600,400);
    //just plotting non-discret pulses | they look pretties in case of mean sig ;-)
    pulseTuple->Draw("sigND:timebin","","L");
    // pulseTuple->Draw("sig:timebin","","Lsame");
    pulseTuple->SetLineColor(3);
    pulseTuple->Draw("sigNDAfterTCF:timebin","","Lsame");
    // pulseTuple->Draw("sigAfterTCF:timebin","","Lsame");
  }
  
  delete [] signalIn;
  delete [] signalOut;
  delete [] signalInD;
  delete [] signalOutD;
 
  return pulseTuple;

}


//____________________________________________________________________________
void AliTPCCalibTCF::PrintPulseThresholds() {
  //
  // Prints the pulse threshold settings
  //

  printf("   %4.0d [ADC] ... expected Gate fluctuation length \n", fGateWidth);
  printf("   %4.0d [ADC] ... expected usefull signal length \n",  fSample);
  printf("   %4.0d [ADC] ... needed pulselength for TC characterisation \n", fPulseLength);
  printf("   %4.0d [ADC] ... lower pulse height limit \n", fLowPulseLim);
  printf("   %4.0d [ADC] ... upper pulse height limit \n", fUpPulseLim);
  printf("   %4.1f [ADC] ... maximal pulse RMS \n", fRMSLim);
  printf("   %4.1f [ADC] ... pulse/tail integral ratio \n", fRatioIntLim);

} 


//____________________________________________________________________________
void AliTPCCalibTCF::MergeHistoPerFile(const char *fileNameIn, const char *fileNameSum, Int_t mode)
{
  // Gets histograms from fileNameIn and adds contents to fileSum
  //
  // If fileSum doesn't exist, fileSum is created
  //   mode = 0, just ONE BIG FILE ('fileSum') will be used
  //   mode = 1, one file per sector ('fileSum-Sec#.root') will be used 
  // mode=1 is much faster, but the additional function 'MergeToOneFile' has to be used in order to  
  // get one big and complete collection file again ...
  //
  // !Make sure not to add the same file more than once!
  
  TFile fileIn(fileNameIn,"READ");
  TH1F *hisIn;                             
  TKey *key;                                          
  TIter next(fileIn.GetListOfKeys());  
  // opens a file, although, it might not be uses (see "mode")
  TFile *fileOut = new TFile(fileNameSum,"UPDATE"); 
  //fileOut.cd();
  
  Int_t nHist=fileIn.GetNkeys();
  Int_t iHist=0;

  Int_t secPrev = -1;
  char fileNameSumSec[100];


  while((key=(TKey*)next())) {
    const char *hisName = key->GetName();

    TString name(key->GetName());
    if (name.Contains("ddl") ) continue;  // ignore the 2d histogramms per ddl

    hisIn=(TH1F*)fileIn.Get(hisName);          
    
    Int_t numPulse=(Int_t)hisIn->GetBinContent(1);
    Int_t sec=(Int_t)hisIn->GetBinContent(2);
    Int_t pulseLength= hisIn->GetNbinsX()-4;    

    // in case of mode 1, store histos in files per sector
    if (sec!=secPrev && mode != 0) {
      if (secPrev>0) { // closing old file
        fileOut->Close();
      }
      // opening new file 
      snprintf(fileNameSumSec,100,"%s-Sec%d.root",fileNameSum,sec);
      fileOut = new TFile(fileNameSumSec,"UPDATE");
      secPrev = sec;
    }

    // search for existing histogram
    TH1F *his=(TH1F*)fileOut->Get(hisName);
    if (iHist%100==0) {
      printf("Histogram %d / %d, %s, Action: ",iHist,nHist,hisName);
      if (!his) {
	printf("NEW\n"); 
      } else {
	printf("ADD\n"); 
      }
    }
    iHist++;
    
    if (!his) {
      his=hisIn;
      his->Write(hisName);
    } else {
      his->AddBinContent(1,numPulse);
      for (Int_t ii=5; ii<pulseLength+5; ii++) {
	his->AddBinContent(ii,hisIn->GetBinContent(ii));
      }
      his->Write(hisName,TObject::kOverwrite);
    }
  }

  printf("closing files (may take a while)...\n");
  fileOut->Close();
  

  fileIn.Close();
  printf("...DONE\n\n");
}


//____________________________________________________________________________
void AliTPCCalibTCF::MergeToOneFile(const char *nameFileSum) {

  // Merges all Sec-files together ...
  // this is an additional functionality for the function MergeHistsPerFile
  // if for example mode=1

  TH1F *hisIn;
  TKey *key;

  // just delete the file entries ...
  TFile fileSumD(nameFileSum,"RECREATE");
  fileSumD.Close();

  char nameFileSumSec[100];

  for (Int_t sec=0; sec<72; sec++) { // loop over all possible filenames

    snprintf(nameFileSumSec,100,"%s-Sec%d.root",nameFileSum,sec);
    TFile *fileSumSec = new TFile(nameFileSumSec,"READ");

    Int_t nHist=fileSumSec->GetNkeys();
    Int_t iHist=0;

    if (nHist) { // file found \ NKeys not empty

      TFile fileSum(nameFileSum,"UPDATE");
      fileSum.cd();

      printf("Sector file %s found\n",nameFileSumSec);
      TIter next(fileSumSec->GetListOfKeys());
      while( (key=(TKey*)next()) ) {
        const char *hisName = key->GetName();
	TString name(hisName);
	if (name.Contains("ddl") ) continue;  // ignore the 2d histogramms per ddl
        hisIn=(TH1F*)fileSumSec->Get(hisName);


        if (iHist%100==0) {
          printf("found histogram %d / %d, %s\n",iHist,nHist,hisName);
        }
        iHist++;

	//        TH1F *his = (TH1F*)hisIn->Clone(hisName);
        hisIn->Write(hisName);

      }
      printf("Saving histograms from sector %d (may take a while) ...",sec);
      fileSum.Close();

    }
    fileSumSec->Close();
  }
  printf("...DONE\n\n");
}


//____________________________________________________________________________
Int_t AliTPCCalibTCF::DumpTCFparamToFilePerPad(const char *nameFileTCFPerPad,const char *nameFileTCFPerSec, const char *nameMappingFile) {
  //
  // Writes TCF parameters per PAD to .data file
  //
  // from now on: "roc" refers to the offline sector numbering
  //              "sector" refers to the 18 sectors per side
  //
  // Gets TCF parameters of single pads from nameFileTCFPerPad and writes them to
  // the file 'tpcTCFparamPAD.data'
  //
  // If there are parameters for a pad missing, then the parameters of the roc,
  // in which the pad is located, are used as the pad parameters. The parameters for
  // the roc are retreived from nameFileTCFPerSec. If there are parameters for
  // a roc missing, then the parameters are set to -1.  

  Float_t k0 = -1, k1 = -1, k2 = -1, l0 = -1, l1 = -1, l2 = -1;
  Int_t roc, row, pad, side, sector, rcu, hwAddr; 
  Int_t entryNum = 0;
  Int_t checksum = 0;
  Int_t tpcPadNum = 557568;
  Int_t validFlag = 1; // 1 if parameters for pad exist, 0 if they are only inherited from the roc

  // get file/tuple with parameters per pad
  TFile fileTCFparam(nameFileTCFPerPad);
  TNtuple *paramTuple = (TNtuple*)fileTCFparam.Get("TCFparam");

  // get mapping file
  // usual location of mapping file: $ALICE_ROOT/TPC/Calib/tpcMapping.root
  TFile *fileMapping = new TFile(nameMappingFile, "read");
  AliTPCmapper *mapping = (AliTPCmapper*) fileMapping->Get("tpcMapping");
  delete fileMapping;

  if (mapping == 0) {
    printf("Failed to get mapping object from %s.  ...\n", nameMappingFile);
    return -1;
  } else {
    printf("Got mapping object from %s\n", nameMappingFile);
  }

  Bool_t *entryID = new Bool_t[7200000]; // helping vector
  for (Int_t ii = 0; ii<7200000; ii++) {
    entryID[ii]=0;
  }

  // creating outputfile
  ofstream fileOut;
  char nameFileOut[255];
  snprintf(nameFileOut,255,"tpcTCFparamPAD.data");
  fileOut.open(nameFileOut);
  // following not used:
  // char headerLine[255];
  // snprintf(headerLine,255,"15\tside\tsector\tRCU\tHWadr\tk0\tk1\tk2\tl0\tl1\tl2\tValidFlag");
  // fileOut << headerLine << std::endl;
  fileOut << "15" << std::endl;
 
  // loop over nameFileTCFPerPad, write parameters into outputfile
  // NOTE: NO SPECIFIC ORDER !!!
  printf("\nstart assigning parameters to pad...\n");  
  for (Int_t iParam = 0; iParam < paramTuple->GetEntries(); iParam++) {
    paramTuple->GetEntry(iParam);
    Float_t *paramArgs = paramTuple->GetArgs();
    roc = Int_t(paramArgs[0]);
    row = Int_t(paramArgs[1]);
    pad = Int_t(paramArgs[2]);
    side = Int_t(mapping->GetSideFromRoc(roc));
    sector = Int_t(mapping->GetSectorFromRoc(roc));
    rcu = Int_t(mapping->GetRcu(roc,row,pad));
    hwAddr = Int_t(mapping->GetHWAddress(roc,row,pad));
    k0 = TMath::Nint(paramArgs[7] * (TMath::Power(2,16) - 1));
    k1 = TMath::Nint(paramArgs[8] * (TMath::Power(2,16) - 1));
    k2 = TMath::Nint(paramArgs[9] * (TMath::Power(2,16) - 1));
    l0 = TMath::Nint(paramArgs[4] * (TMath::Power(2,16) - 1));
    l1 = TMath::Nint(paramArgs[5] * (TMath::Power(2,16) - 1));
    l2 = TMath::Nint(paramArgs[6] * (TMath::Power(2,16) - 1));
    if (entryNum%10000==0) {
      printf("assigned pad %i / %i\n",entryNum,tpcPadNum);
    }
    
    fileOut << entryNum++ << "\t" << side << "\t" << sector << "\t" << rcu << "\t" << hwAddr << "\t";
    fileOut << k0 << "\t" << k1 << "\t" << k2 << "\t" << l0 << "\t" << l1 << "\t" << l2 << "\t" << validFlag << std::endl;
    entryID[roc*100000 + row*1000 + pad] = 1;
  }

  // Wrote all found TCF params per pad into data file
  // NOW FILLING UP THE REST WITH THE PARAMETERS FROM THE ROC MEAN
  
  // get file/tuple with parameters per roc
  TFile fileSecTCFparam(nameFileTCFPerSec);
  TNtuple *paramTupleSec = (TNtuple*)fileSecTCFparam.Get("TCFparam");

  // loop over all pads and get/write parameters for pads which don't have
  // parameters assigned yet
  validFlag = 0; 
  for (roc = 0; roc<72; roc++) {
    side = Int_t(mapping->GetSideFromRoc(roc));
    sector = Int_t(mapping->GetSectorFromRoc(roc));
    for (Int_t iParamSec = 0; iParamSec < paramTupleSec->GetEntries(); iParamSec++) {
      paramTupleSec->GetEntry(iParamSec);
      Float_t *paramArgsSec = paramTupleSec->GetArgs();
      if ((paramArgsSec[0]-roc)<1e-7) { // if roc is found
	k0 = TMath::Nint(paramArgsSec[7] * (TMath::Power(2,16) - 1));
	k1 = TMath::Nint(paramArgsSec[8] * (TMath::Power(2,16) - 1));
	k2 = TMath::Nint(paramArgsSec[9] * (TMath::Power(2,16) - 1));
	l0 = TMath::Nint(paramArgsSec[4] * (TMath::Power(2,16) - 1));
	l1 = TMath::Nint(paramArgsSec[5] * (TMath::Power(2,16) - 1));
	l2 = TMath::Nint(paramArgsSec[6] * (TMath::Power(2,16) - 1));
	break;
      } else {
	k0 = k1 = k2 = l0 = l1 = l2 = -1;
      }
    }
    for (row = 0; row<mapping->GetNpadrows(roc); row++) {
      for (pad = 0; pad<mapping->GetNpads(roc,row); pad++) {
	if (entryID[roc*100000 + row*1000 + pad]==1) {
	  continue;
	}

	entryID[roc*100000 + row*1000 + pad] = 1;
	rcu = Int_t(mapping->GetRcu(roc,row,pad));
	hwAddr = Int_t(mapping->GetHWAddress(roc,row,pad));
	if (entryNum%10000==0) {
	  printf("assigned pad %i / %i\n",entryNum,tpcPadNum);
	}

	fileOut << entryNum++ << "\t" << side << "\t" << sector << "\t" << rcu << "\t" << hwAddr << "\t";
	fileOut << k0 << "\t" << k1 << "\t" << k2 << "\t" << l0 << "\t" << l1 << "\t" << l2 << "\t" << validFlag << std::endl;
      }
    }
  }

  printf("assigned pad %i / %i\ndone assigning\n",entryNum,tpcPadNum);
  
  // check if correct amount of sets of parameters were written
  for (Int_t ii = 0; ii<7200000; ii++) {
    checksum += entryID[ii];
  }
  if (checksum == tpcPadNum) {
    printf("checksum ok, sets of parameters written = %i\n",checksum);
  } else {
    printf("\nCHECKSUM WRONG, sets of parameters written = %i, should be %i\n\n",checksum,tpcPadNum);
  }
  
  // closing & destroying
  fileOut.close();
  fileTCFparam.Close();
  fileSecTCFparam.Close();
  delete [] entryID;
  printf("output written to file: %s\n",nameFileOut);
  return 0;
}



//____________________________________________________________________________
Int_t AliTPCCalibTCF::DumpTCFparamToFilePerSector(const char *nameFileTCFPerSec, const char *nameMappingFile) {
  //
  // Writes TCF parameters per SECTOR (=ROC) to .data file
  //
  // from now on: "roc" refers to the offline sector numbering
  //              "sector" refers to the 18 sectors per side
  //
  // Gets TCF parameters of a roc from nameFileTCFPerSec and writes them to
  // the file 'tpcTCFparamSector.data'
  //
  // If there are parameters for a roc missing, then the parameters are set to -1
  
  Float_t k0 = -1, k1 = -1, k2 = -1, l0 = -1, l1 = -1, l2 = -1;
  Int_t entryNum = 0;
  Int_t validFlag = 0; // 1 if parameters for roc exist
  
  // get file/tuple with parameters per roc
  TFile fileTCFparam(nameFileTCFPerSec);
  TNtuple *paramTupleSec = (TNtuple*)fileTCFparam.Get("TCFparam");
  
  
  // get mapping file
  // usual location of mapping file: $ALICE_ROOT/TPC/Calib/tpcMapping.root
  TFile *fileMapping = new TFile(nameMappingFile, "read");
  AliTPCmapper *mapping = (AliTPCmapper*) fileMapping->Get("tpcMapping");
  delete fileMapping;
  
  if (mapping == 0) {
    printf("Failed to get mapping object from %s.  ...\n", nameMappingFile);
    return -1;
  } else {
    printf("Got mapping object from %s\n", nameMappingFile);
  }
  
  
  // creating outputfile
  
  ofstream fileOut;
  char nameFileOut[255];
  snprintf(nameFileOut,255,"tpcTCFparamSector.data");
  fileOut.open(nameFileOut);
  // following not used:   
  // char headerLine[255];
  // snprintf(headerLine,255,"16\tside\tsector\tRCU\tHWadr\tk0\tk1\tk2\tl0\tl1\tl2\tValidFlag");
  // fileOut << headerLine << std::endl;
  fileOut << "16" << std::endl;
  
  // loop over all rcu's in the TPC (6 per sector)
  printf("\nstart assigning parameters to rcu's...\n");
  for (Int_t side = 0; side<2; side++) {
    for (Int_t sector = 0; sector<18; sector++) {
      for (Int_t rcu = 0; rcu<6; rcu++) {
	
	validFlag = 0;
	Int_t roc = Int_t(mapping->GetRocFromPatch(side, sector, rcu));
	
	// get parameters (through loop search) for rcu from corresponding roc
	for (Int_t iParam = 0; iParam < paramTupleSec->GetEntries(); iParam++) {
	  paramTupleSec->GetEntry(iParam);
	  Float_t *paramArgs = paramTupleSec->GetArgs();
	  if ((paramArgs[0]-roc)<1e-7) { // if roc is found
	    validFlag = 1; 
	    k0 = TMath::Nint(paramArgs[7] * (TMath::Power(2,16) - 1));
	    k1 = TMath::Nint(paramArgs[8] * (TMath::Power(2,16) - 1));
	    k2 = TMath::Nint(paramArgs[9] * (TMath::Power(2,16) - 1));
	    l0 = TMath::Nint(paramArgs[4] * (TMath::Power(2,16) - 1));
	    l1 = TMath::Nint(paramArgs[5] * (TMath::Power(2,16) - 1));
	    l2 = TMath::Nint(paramArgs[6] * (TMath::Power(2,16) - 1));
	    break;
	  }
	}
	if (!validFlag) { // No TCF parameters found for this roc 
	  k0 = k1 = k2 = l0 = l1 = l2 = -1;
	}
	
	fileOut << entryNum++ << "\t" << side << "\t" << sector << "\t" << rcu << "\t" << -1 << "\t";
	fileOut << k0 << "\t" << k1 << "\t" << k2 << "\t" << l0 << "\t" << l1 << "\t" << l2 << "\t" << validFlag << std::endl;
      }
    }
  }

  printf("done assigning\n");
  
  // closing files
  fileOut.close();
  fileTCFparam.Close();
  printf("output written to file: %s\n",nameFileOut);
  return 0;

}
