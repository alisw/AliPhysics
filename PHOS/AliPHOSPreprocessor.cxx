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
#include "AliCDBEntry.h"
#include "AliPHOSEmcCalibData.h"
#include "TFile.h"
#include "TH1.h"
#include "TF1.h"
#include "TH2.h"
#include "TMap.h"
#include "TRandom.h"
#include "TKey.h"
#include "TList.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "TObject.h"
#include "TString.h"
#include "TMath.h"
#include "TAxis.h"
#include "AliPHOSEmcBadChannelsMap.h"

ClassImp(AliPHOSPreprocessor)

  Double_t rgaus(Double_t *x, Double_t *par)
{
  Double_t gaus = par[0] * TMath::Exp( -(x[0]-par[1])*(x[0]-par[1]) /
                                       (2*par[2]*par[2]) );
  return gaus;
}

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

  AddRunType("PHYSICS");
  AddRunType("LED");
}

//_______________________________________________________________________________________
UInt_t AliPHOSPreprocessor::Process(TMap* /*valueSet*/)
{
  // process data retrieved by the Shuttle
  
  TString runType = GetRunType();
  Log(Form("Run type: %s",runType.Data()));

  if(runType=="LED") {
    Bool_t ledOK = ProcessLEDRun();
    Bool_t badmap_OK = FindBadChannelsEmc();
    if(!badmap_OK) Log(Form("WARNING! FindBadChannels() completed with BAD status!"));
    if(ledOK) return 0;
    else
      return 1;
  }
  
  if(runType=="PHYSICS") {
    
    Bool_t calibEmc_OK = CalibrateEmc();
    
    if(calibEmc_OK) return 0;
    else
      return 1;
  }

  Log(Form("Unknown run type %s. Do nothing and return OK.",runType.Data()));
  return 0;

}


Bool_t AliPHOSPreprocessor::ProcessLEDRun()
{
  //Process LED run, update High Gain/Low Gain ratios.

  AliPHOSEmcCalibData calibData;

  TList* list = GetFileSources(kDAQ, "LED");
  if(!list) {
    Log("Sources list for LED run not found, exit.");
    return kFALSE;
  }

  if(!list->GetEntries()) {
    Log("Sources list for LED run is empty, exit.");
    return kFALSE;
  }

  //Retrieve the last EMC calibration object
  const AliPHOSEmcCalibData* clb=0;
  AliCDBEntry* entryCalib = GetFromOCDB("Calib", "EmcGainPedestals");

  if(!entryCalib)
    Log(Form("Cannot find any AliCDBEntry for [Calib, EmcGainPedestals]!"));
  else
    clb = (AliPHOSEmcCalibData*)entryCalib->GetObject();
  
  TIter iter(list);
  TObjString *source;
  
  while ((source = dynamic_cast<TObjString *> (iter.Next()))) {
    
    AliInfo(Form("found source %s", source->String().Data()));

    TString fileName = GetFile(kDAQ, "LED", source->GetName());
    AliInfo(Form("Got filename: %s",fileName.Data()));

    TFile f(fileName);

    if(!f.IsOpen()) {
      Log(Form("File %s is not opened, something goes wrong!",fileName.Data()));
      return kFALSE;
    }
    
    TH1I* fFiredCells = (TH1I*)f.Get("fFiredCells");
    
    if(fFiredCells) {
      const Double_t nFiredCells = fFiredCells->GetMean();
      Log(Form("Number of fired cells per event is %.1f",nFiredCells));
    }
    
    const Int_t nMod=5; // 1:5 modules
    const Int_t nCol=56; //1:56 columns in each module
    const Int_t nRow=64; //1:64 rows in each module
    
    for(Int_t mod=0; mod<nMod; mod++) {
      for(Int_t col=0; col<nCol; col++) {
	for(Int_t row=0; row<nRow; row++) {
	  
	  if(clb) {
	    Float_t  hg2lg = clb->GetHighLowRatioEmc(5-mod,col+1,row+1);
	    Double_t coeff = clb->GetADCchannelEmc(5-mod,col+1,row+1);
	    calibData.SetADCchannelEmc(5-mod,col+1,row+1,coeff);
	    calibData.SetHighLowRatioEmc(5-mod,col+1,row+1,hg2lg);
	  }	
  
	  //High Gain to Low Gain ratio
	  Float_t ratio = HG2LG(mod,row,col,&f);
	  if(ratio != 16.) {
	    calibData.SetHighLowRatioEmc(5-mod,col+1,row+1,ratio);
	    AliInfo(Form("mod %d iX %d iZ %d  ratio %.3f\n",mod,row,col,ratio));
	  }
	}
      }
    }

  } // end of loop over files

  //Store the updated High Gain/Low Gain ratios
  AliCDBMetaData emcMetaData;

  //Data valid from current run until updated (validityInfinite=kTRUE)
  //Bool_t result = Store("Calib","EmcGainPedestals",&calibData,&emcMetaData,0,kTRUE);

  //Store reference data
  Bool_t refOK = StoreReferenceLED(list);
  if(refOK) Log(Form("LED reference data successfully stored."));
  
  //return result;
  return kTRUE;
}

Float_t AliPHOSPreprocessor::HG2LG(Int_t mod, Int_t X, Int_t Z, TFile* f)
{
  //Calculates High gain to Low gain ratio 
  //for crystal at the position (X,Z) in the PHOS module mod.
  
  char hname[128];
  sprintf(hname,"%d_%d_%d",mod,X,Z);

  TH1F* h1 = (TH1F*)f->Get(hname);
  if(!h1) return 16.;
  
  if(h1->GetEntries()<2000.) return 16.;
  
  if(h1->GetMaximum()<10.) h1->Rebin(4);
  if(h1->GetMaximum()<10.) return 16.;
  
  Double_t max = h1->GetBinCenter(h1->GetMaximumBin()); // peak
  Double_t xmin = max - (h1->GetRMS()/3);
  Double_t xmax = max + (h1->GetRMS()/2);
  //       Double_t xmin = max - (h1->GetRMS());
  //       Double_t xmax = max + (h1->GetRMS());

  TF1* gaus1 = new TF1("gaus1",rgaus,xmin,xmax,3);
  gaus1->SetParNames("Constant","Mean","Sigma");
  gaus1->SetParameter("Constant",h1->GetMaximum());
  gaus1->SetParameter("Mean",max);
  gaus1->SetParameter("Sigma",1.);
  gaus1->SetLineColor(kBlue);
  
  Double_t mean_min = h1->GetXaxis()->GetXmin();
  Double_t mean_max = h1->GetXaxis()->GetXmax();
  gaus1->SetParLimits(1,mean_min,mean_max);
  
  h1->Fit(gaus1,"RQ+");
  Double_t hg2lg = gaus1->GetParameter("Mean");
  if( (hg2lg-mean_min<0.001) || (mean_max-hg2lg<0.001)) hg2lg=max;
  
  AliInfo(Form("%s: %.1f entries, mean=%.3f, peak=%.3f, rms= %.3f. HG/LG = %.3f\n",
	  h1->GetTitle(),h1->GetEntries(),h1->GetMean(),max,h1->GetRMS(),hg2lg)); 

  return hg2lg;
  
}

Bool_t AliPHOSPreprocessor::FindBadChannelsEmc()
{
  //Loop over two systems: DAQ and HLT.
  //For each system the same algorithm implemented in DoFindBadChannelsEmc() invokes.

  TList* list=0;
  TString path;
  
  Int_t system[2] = { kDAQ, kHLT };
  const char* sysn[] = { "DAQ","HLT" };
  Bool_t result[2] = { kTRUE, kTRUE };

  for (Int_t i=0; i<2; i++) {
    
    if(system[i] == kHLT) continue;
    
    AliPHOSEmcBadChannelsMap badMap;
    list = GetFileSources(system[i], "BAD_CHANNELS");

    if(!list) {
      Log(Form("%s sources list for BAD_CHANNELS not found!",sysn[i]));
      result[i] = kFALSE;
      continue;
    }

    if(!list->GetEntries()) {
      Log(Form("Got empty sources list. It seems %s DA2 did not produce any files!",sysn[i]));
      result[i] = kFALSE;
      continue;
    }

    result[i] *= DoFindBadChannelsEmc(system[i],list,badMap);

    // Store the bad channels map.
  
    AliCDBMetaData md;
    md.SetResponsible("Boris Polishchuk");

    if(system[i] == kDAQ) 
      path = "Calib";
    else 
      path = "HLT";
  
    // Data valid from current run until being updated (validityInfinite=kTRUE)
    result[i] *= Store(path.Data(), "EmcBadChannels", &badMap, &md, 0, kTRUE);
    
  }
  
  if(result[0] || result[1]) return kTRUE;
  else return kFALSE;
}

Bool_t AliPHOSPreprocessor::DoFindBadChannelsEmc(Int_t system, TList* list, AliPHOSEmcBadChannelsMap& badMap)
{
  //Creates the bad channels map for PHOS EMC.

  // The file fileName contains histograms which have been produced by DA2 detector algorithm.
  // It is a responsibility of the SHUTTLE framework to form the fileName.

  TIter iter(list);
  TObjString *source;
  char hnam[80];
  TH1F* h1=0;

  const Float_t fQualityCut = 1.;
  Int_t nGoods[5] = {0,0,0,0,0};
  
  while ((source = dynamic_cast<TObjString *> (iter.Next()))) {

    AliInfo(Form("found source %s", source->String().Data()));

    TString fileName = GetFile(system, "BAD_CHANNELS", source->GetName());
    AliInfo(Form("Got filename: %s",fileName.Data()));

    TFile f(fileName);

    if(!f.IsOpen()) {
      Log(Form("File %s is not opened, something goes wrong!",fileName.Data()));
      return kFALSE;
    }

    Log(Form("Begin check for bad channels."));

    for(Int_t mod=0; mod<5; mod++) {
      for(Int_t iX=0; iX<64; iX++) {
	for(Int_t iZ=0; iZ<56; iZ++) {
	  
	  sprintf(hnam,"%d_%d_%d_%d",mod,iX,iZ,1); // high gain	
	  h1 = (TH1F*)f.Get(hnam);

	  if(h1) {
	    Double_t mean = h1->GetMean();
	    
	    if(mean)
	      Log(Form("iX=%d iZ=%d gain=%d   mean=%.3f\n",iX,iZ,1,mean));

	    if( mean>0 && mean<fQualityCut ) { 
	      nGoods[mod]++; 
	    }
	    else
	      badMap.SetBadChannel(mod+1,iZ+1,iX+1); //module, col,row
	  }

	}
      }

      if(nGoods[mod])
	Log(Form("Module %d: %d good channels.",mod,nGoods[mod]));
    }

    
  } // end of loop over sources
  
  return kTRUE;
}

Bool_t AliPHOSPreprocessor::CalibrateEmc()
{
  //Loop over two systems: DAQ and HLT.
  //For each system the same algorithm implemented in DoCalibrateEmc() invokes.

  AliPHOSEmcCalibData*   lastCalib=0;
  const AliPHOSEmcBadChannelsMap* badMap=0;
  AliCDBEntry* entryBCM=0;
  AliCDBEntry* entryEmc=0;
  TList* list=0;
  TString path;
  
  Int_t system[2] = { kDAQ, kHLT };
  const char* sysn[] = { "DAQ","HLT" };
  Bool_t result[2] = { kTRUE, kTRUE };

  for (Int_t i=0; i<2; i++) {

    if(system[i] == kHLT) continue;

    AliPHOSEmcCalibData calibData;
    list = GetFileSources(system[i], "AMPLITUDES");
  
    if(!list) {
      Log(Form("%s sources list not found!",sysn[i]));
      result[i] = kFALSE;
      continue;
    }

    if(!list->GetEntries()) {
      Log(Form("Got empty sources list. It seems %s DA1 did not produce any files!",sysn[i]));
      result[i] = kFALSE;
      continue;
    }

    // Retrieve the Bad Channels Map (BCM)

    if(system[i] == kDAQ) 
      path = "Calib";
    else 
      path = "HLT";
  
    entryBCM = GetFromOCDB(path.Data(), "EmcBadChannels");

    if(!entryBCM)
      Log(Form("WARNING!! Cannot find any AliCDBEntry for [%s, EmcBadChannels]!",path.Data()));
    else
      badMap = (AliPHOSEmcBadChannelsMap*)entryBCM->GetObject();

    if(!badMap)
      Log(Form("WARNING!! Nothing for %s in AliCDBEntry. All cells considered GOOD!",sysn[i]));

    // Retrieve  the last EMC calibration object
    
    entryEmc = GetFromOCDB(path.Data(), "EmcGainPedestals");
    
    if(!entryEmc) 
      Log(Form("Cannot find any EmcGainPedestals entry for this run and path %s",path.Data()));
    else
      lastCalib = (AliPHOSEmcCalibData*)entryEmc->GetObject();

    if(lastCalib) 
      result[i] *= DoCalibrateEmc(system[i],list,badMap,*lastCalib);    
    else 
      result[i] *= DoCalibrateEmc(system[i],list,badMap,calibData);
    
    //Store EMC calibration data
    AliCDBMetaData emcMetaData;
    
    // if(lastCalib)
    //   result[i] *= Store(path.Data(), "EmcGainPedestals", lastCalib, &emcMetaData, 0, kFALSE);
    // else
    //   result[i] *= Store(path.Data(), "EmcGainPedestals", &calibData, &emcMetaData, 0, kFALSE);

    //Store reference data
    Bool_t refOK = StoreReferenceEmc(system[i],list);
    if(refOK) Log(Form("Reference data for %s amplitudes successfully stored.",sysn[i]));
    
  }
  
  if(result[0] || result[1]) return kTRUE;
  else return kFALSE;
}

Bool_t AliPHOSPreprocessor::StoreReferenceEmc(Int_t system, TList* list)
{
  //Put 2D calibration histograms (E vs Time) prepared by DAQ/HLT to the reference storage.
  //system is DAQ or HLT, TList is the list of FES sources.

  if(system!=kDAQ) return kFALSE;

  TObjString *source = dynamic_cast<TObjString *> (list->First());
  if(!source) return kFALSE;

  TString fileName = GetFile(system, "AMPLITUDES", source->GetName());

  Bool_t resultRef = StoreReferenceFile(fileName.Data(),"CalibRefPHOS.root");
  return resultRef;

}

Bool_t AliPHOSPreprocessor::StoreReferenceLED(TList* list)
{
  //Put HG/LG histograms to the reference storage.
  
  TObjString *source = dynamic_cast<TObjString *> (list->First());
  if(!source) return kFALSE;
  
  TString fileName = GetFile(kDAQ, "LED", source->GetName());
  
  Bool_t resultRef = StoreReferenceFile(fileName.Data(),"LEDRefPHOS.root");
  return resultRef;
  
}


Bool_t AliPHOSPreprocessor::DoCalibrateEmc(Int_t system, TList* list, const AliPHOSEmcBadChannelsMap* badMap, AliPHOSEmcCalibData& calibData)
{
  // Return kTRUE if OK.
  // I.  Calculates the set of calibration coefficients to equalyze the mean energies deposited at high gain.
  // II. Extracts High_Gain/Low_Gain ratio for each channel.

  // The file fileName contains histograms which have been produced by DA1 detector algorithm.
  // It is a responsibility of the SHUTTLE framework to form the fileName.

  gRandom->SetSeed(0); //the seed is set to the current  machine clock!
  Int_t minEntries=1000; // recalculate calibration coeff. if Nentries > minEntries.

  TIter iter(list);
  TObjString *source;
  
  while ((source = dynamic_cast<TObjString *> (iter.Next()))) {
    AliInfo(Form("found source %s", source->String().Data()));

    TString fileName = GetFile(system, "AMPLITUDES", source->GetName());
    AliInfo(Form("Got filename: %s",fileName.Data()));

    TFile f(fileName);

    if(!f.IsOpen()) {
      Log(Form("File %s is not opened, something goes wrong!",fileName.Data()));
      return kFALSE;
    }
    
    const Int_t nMod=5; // 1:5 modules
    const Int_t nCol=56; //1:56 columns in each module
    const Int_t nRow=64; //1:64 rows in each module

    Double_t coeff;
    char hnam[80];
    TH2F* h2=0;
    TH1D* h1=0;
    
    //Get the reference histogram
    //(author: Gustavo Conesa Balbastre)

    TList * keylist = f.GetListOfKeys();
    Int_t nkeys   = f.GetNkeys();
    Bool_t ok = kFALSE;
    TKey  *key;
    Int_t ikey = 0;
    Int_t counter = 0;
    TH1D* hRef = 0;

    //Check if the file contains any histogram
    
    if(nkeys< 2){
      Log(Form("Not enough histograms (%d) for calibration.",nkeys));
      return 1; // it's not fatal! May be short run..
    }
    
    while(!ok){
      ikey = gRandom->Integer(nkeys);
      key = (TKey*)keylist->At(ikey);
      TObject* obj = f.Get(key->GetName());
      TString cname(obj->ClassName());
      if(cname == "TH2F") {
	h2 = (TH2F*)obj;
	TString htitl = h2->GetTitle();
	if(htitl.Contains("and gain 1")) {
	  hRef = h2->ProjectionX();
	  hRef->GetXaxis()->SetRangeUser(10.,1000.); // to cut off saturation peak and noise
	  // Check if the reference histogram has too little statistics
	  if(hRef->GetMean() && hRef->GetEntries()>minEntries) ok=kTRUE;

	  const TString delim = "_";
	  TString str = hRef->GetName();
	  TObjArray* tks = str.Tokenize(delim);
	  const Int_t md = ((TObjString*)tks->At(0))->GetString().Atoi();
	  const Int_t X   = ((TObjString*)tks->At(1))->GetString().Atoi();
	  const Int_t Z   = ((TObjString*)tks->At(2))->GetString().Atoi();

	  if(badMap) {
	    if(badMap->IsBadChannel(5-md,Z+1,X+1)) {
	      AliInfo(Form("Cell mod=%d col=%d row=%d is bad. Histogram %s rejected.",
			   5-md,Z+1,X+1,hRef->GetName()));
	      ok=kFALSE;
	    }
	  }

	}
      }
      
      counter++;
      
      if(!ok && counter > nkeys){
	Log("No histogram with enough statistics for reference. Exit.");
	return 1; // Not fatal, just wait..
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
	  
	  sprintf(hnam,"%d_%d_%d_1",mod,row,col); // high gain!
	  h2 = (TH2F*)f.Get(hnam);
	  
	  //TODO: dead channels exclusion!
	  if(h2) {
	    h1 = h2->ProjectionX();
	    h1->GetXaxis()->SetRangeUser(10.,1000.); //to cut off saturation peak and noise
	    coeff = h1->GetMean()/refMean;
	    if(coeff>0 && h1->GetEntries()>minEntries) {
	      calibData.SetADCchannelEmc(5-mod,col+1,row+1,0.005/coeff);
	      AliInfo(Form("mod %d col %d row %d  coeff %f\n",mod,col,row,coeff));
	    }
	  }
	}
      }
    }
    
    f.Close();
  }
  
  return 1;
}

     

