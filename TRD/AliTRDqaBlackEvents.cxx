/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * withount fee, provided that the abov copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is without express or implied warranty.                  *
 **************************************************************************/

/* $Id: AliTRDqaBlackEvents.cxx 23387 2008-01-17 17:25:16Z cblume $ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  QA of black events                                                    //
//                                                                        //
//  Author:                                                               //
//    Sylwester Radomski (radomski@physi.uni-heidelberg.de)               //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include "TH1D.h"
#include "TH2D.h"
#include "TH2S.h"
#include "TH3F.h"
#include "TF1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TGraph.h"

#include "AliTRDgeometry.h"
#include "AliTRDrawStreamTB.h"
#include "AliTRDqaBlackEvents.h"

ClassImp(AliTRDqaBlackEvents)

///////////////////////////////////////////////////////////////////////////////////////////////////

AliTRDqaBlackEvents::AliTRDqaBlackEvents() 
  :TObject() 
  ,fnEvents(0)
  ,fCreateFull(0)
  ,fThresh(0)
  ,fCount(0)
  ,fOccupancy(0)
  ,fDetRob(0)
  ,fTBEvent(0)
  ,fRefHistPed(0)
  ,fRefHistNoise(0)
  ,fErrorHC(0)
  ,fErrorMCM(0)
  ,fErrorADC(0)
  ,fErrorSMHC(0)
  ,fErrorSMMCM(0)
  ,fErrorSMADC(0)
  ,fErrorGraphHC(0)
  ,fErrorGraphMCM(0)
  ,fErrorGraphADC(0)
  ,fGraphMCM(0)
  ,fMcmTracks(0)
  ,fMapMCM(0)
  ,fFracMCM(0)
  ,fSMHCped(0)
  ,fSMHCerr(0)
  ,fNoiseTotal(0)
  ,fPP(0)
  ,fMinNoise(0.5)
  ,fMaxNoise(2)
  ,fFitType(0) 
  //  ,fRefFileName("")
{
  //
  // Constructor 
  // to create the histograms call Init()
  //

  strcpy(fRefFileName, "");
}

///////////////////////////////////////////////////////////////////////////////////////////////////

AliTRDqaBlackEvents::AliTRDqaBlackEvents(const AliTRDqaBlackEvents &qa) 
  :TObject(qa) 
  ,fnEvents(0)
  ,fCreateFull(0)
  ,fThresh(0)
  ,fCount(0)
  ,fOccupancy(0)
  ,fDetRob(0)
  ,fTBEvent(0)
  ,fRefHistNoise(0)
  ,fErrorHC(0)
  ,fErrorMCM(0)
  ,fErrorADC(0)
  ,fErrorSMHC(0)
  ,fErrorSMMCM(0)
  ,fErrorSMADC(0)
  ,fErrorGraphHC(0)
  ,fErrorGraphMCM(0)
  ,fErrorGraphADC(0)
  ,fGraphMCM(0)
  ,fMcmTracks(0)
  ,fMapMCM(0)
  ,fFracMCM(0)
  ,fSMHCped(0)
  ,fSMHCerr(0)
  ,fNoiseTotal(0)
  ,fPP(0)
  ,fMinNoise(0.5)
  ,fMaxNoise(2) 
  ,fFitType(0)
   //,fRefFileName("")
{
  //
  // Copy constructor 
  // to create the histograms call Init()
  //
  
  strcpy(fRefFileName, "");
}

///////////////////////////////////////////////////////////////////////////////////////////////////

void AliTRDqaBlackEvents::Init() 
{
  //
  // creates histograms 
  // 

  //TFile *file = new 
  //Info("Init", "Statring");

  fnEvents = 0;

  // histograms for chambers
  for(Int_t det=0; det<kDET; det++) {

    fNPoint[det]  = new TH2D(Form("entries_%d", det), "",  16, -0.5, 15.5, 144, -0.5, 143.5);
    fData[det]    = new TH3F(Form("data_%d", det), "", 16, -0.5, 15.5, 144, -0.5, 143.5, 50, -0.5, 49.5);

    // pedestal noise maps using RMS and Fit
    fChPed[det]   = new TH2D(Form("ped_%d", det), "", 16, -0.5, 15.5, 144, -0.5, 143.5);
    fChNoise[det] = new TH2D(Form("noise_%d", det), "", 16, -0.5, 15.5, 144, -0.5, 143.5);
    
    //fChPed[det]   = new TH2D(Form("ped_%d", det), "", 16, -0.5, 15.5, 144, -0.5, 143.5);
    //fChNoise[det] = new TH2D(Form("noise_%d", det), "", 16, -0.5, 15.5, 144, -0.5, 143.5);    

    // distribution per detector
    fPed[det]     = new TH1D(Form("pedDist_%d", det), ";pedestals (ADC counts)", 100, 5, 15);
    fNoise[det]   = new TH1D(Form("noiseDist_%d", det), ";noise (ADC counts)", 100, 0, 5); 
    fSignal[det]  = new TH1D(Form("signal_%d", det), ";signal (ADC counts)", 100, -0.5, 99.5);
    fChPP[det]    = new TH1D(Form("pp_%d", det), ";pp (ADC)", 200, -0.5, 199.5);

    fnEntriesRM[det] = new TH2D(Form("entriesRM_%d", det), ";ROB,MCM", 8, -0.5, 7.5, 16, -0.5, 15.5);
    
    // histograms after reference subtraction
    fChPedRes[det]    = new TH2D(Form("pedRef_%d", det), "", 16, -0.5, 15.5, 144, -0.5, 143.5);
    fChNoiseRes[det]  = new TH2D(Form("noiseRef_%d", det), "", 16, -0.5, 15.5, 144, -0.5, 143.5);


    // error codes
    fErrorLocMCM[det] = new TH2D(Form("errorLocMCM_%d", det), "", 16, -0.5, 15.5, 144, -0.5, 143.5);
    fErrorLocADC[det] = new TH2D(Form("errorLocADC_%d", det), "", 16, -0.5, 15.5, 144, -0.5, 143.5);    
    fErrorLocHC[det]  = new TH2D(Form("errorLocHC_%d", det), "", 16, -0.5, 15.5, 144, -0.5, 143.5);
  }

  // histogram for each MCM
  for(Int_t i=0; i < kDET * kROB * kMCM; i++)
    fFullCounter[i] = 0;

  // histograms from the whole detector
  fOccupancy = new TH1D("occupancy", "", 20, -0.5, 19.5);
  fDetRob    = new TH2D("DetRob", ";detector;ROB", kDET, -0.5, 539.5, 8, -0.5, 7.5);
  fTBEvent   = new TH2D("tbEvent", ";event ID;time bin", 100, -0.5, 99.5, 30, -0.5, 29.5);

  // errors statistics and location
  fErrorHC  = new TH1D("errorHC", ";error ID;", 7, -0.5, 6.5);
  fErrorMCM = new TH1D("errorMCM", ";error ID;", 7, -0.5, 6.5);
  fErrorADC = new TH1D("errorADC", ";error ID;", 7, -0.5, 6.5);
  
  fErrorSMHC  = new TH1D("errorSM_HC", ";SM id", 18, -0.5, 17.5);
  fErrorSMMCM = new TH1D("errorSM_MCM", ";SM id", 18, -0.5, 17.5);
  fErrorSMADC = new TH1D("errorSM_ADC", ";SM id", 18, -0.5, 17.5);


  fErrorGraphHC  = new TGraph();
  fErrorGraphMCM = new TGraph();
  fErrorGraphADC = new TGraph();

  fGraphMCM = new TGraph();
  
  fMapMCM = new TH2D("mapMCM", ";det;mcm", 540, -0.5, 539.5, kROB*kMCM, -0.5, kROB*kMCM-0.5);
  fFracMCM = new TH1D("fracMCM", ";frequency", 100, 0, 1);
  

  fErrorGraphHC->GetHistogram()->SetTitle("Error HC;event number;fraction with error (%)");
  fErrorGraphMCM->GetHistogram()->SetTitle("Error MCM;event number;fraction with error (%)");
  fErrorGraphADC->GetHistogram()->SetTitle("Error ADC;event number;fraction with error (%)"); 


  fSMHCped = new TH2D("smHcPed", ";super module;half chamber", 18, -0.5, 17.5, 60, -0.5, 59.5);
  //fSMHCerr = 0;

  //Info("Init", "Done");

  // number of ADC channels fired per SM and in total
  for(Int_t sm=0; sm<kSM+1; sm++)
    fNumberADC[sm] = new TGraph();

  //
  fNoiseTotal = new TH1D("noiseTotal", "noise (ADC)", 250, 0, 10);
  fPP = new TH1D("peakPeak", "p-p (ADC)", 200, -0.5, 199.5);

  for(Int_t sm=0; sm<kSM; sm++) {
    fSmNoiseRms[sm] = new TH1D(Form("noiseRms_sm%d", sm), ";noise from RMS (ADC)", 100, 0, 10);
    fSmNoiseFit[sm] = new TH1D(Form("noiseFit_sm%d", sm), ";noise frim Fit (ADC)", 100, 0, 10);
    fSmPP[sm] = new TH1D(Form("peakPeak_sm%d", sm), ";peak-peak (ADC)", 200, -0.5, 199.5); 
  }

  fMcmTracks = new TObjArray();
}


///////////////////////////////////////////////////////////////////////////////////////////////////

void AliTRDqaBlackEvents::Reset() 
{
  //
  // Resets the histograms
  //

  for(Int_t i=0; i<kDET; i++) {
    fData[i]->Reset();
    fChPed[i]->Reset();
    fChNoise[i]->Reset();
  }
}


///////////////////////////////////////////////////////////////////////////////////////////////////

void AliTRDqaBlackEvents::SetRefFile(const char *filename) {
  
  strcpy(fRefFileName, filename);
}

///////////////////////////////////////////////////////////////////////////////////////////////////

void AliTRDqaBlackEvents::ReadRefHists(Int_t det) {
  
  fRefHistPed = 0;
  fRefHistNoise = 0;
  
  TFile *file = 0;
  if (fRefFileName) TFile::Open(fRefFileName);
  if (!file) return;

  fRefHistPed   = (TH2D*)file->Get(Form("ped_%d",det));
  fRefHistNoise = (TH2D*)file->Get(Form("noise_%d", det));

  if (file) file->Close();
}

///////////////////////////////////////////////////////////////////////////////////////////////////

Int_t AliTRDqaBlackEvents::AddEvent(AliTRDrawStreamTB *data) 
{
  //
  // Add an event
  //

  // structure to keep track if particular chanel is used
  Char_t isUsed[kDET][kCOL][kPAD]; 
  for(Int_t i=0; i<kDET; i++)
    for(Int_t j=0; j<kCOL; j++)
      for(Int_t k=0; k<kPAD; k++)
   	isUsed[i][j][k] = 0;
  

  // clear the mcm data
  for(Int_t i=0; i < kDET * kROB * kMCM; i++) {
    if (fFullSignal[i]) fFullSignal[i]->Reset();
    fFullCounter[i] = 0;
  }

  Int_t nb = 0;
  
  Int_t lastdet  = -1;
  Int_t lastside = -1;
  Int_t lastmcm  = -1;

  Int_t rob_last  = -1;
  Int_t mcm_last  = -1;

  Int_t nGoodHC  = 0;
  Int_t nGoodMCM = 0;
  Int_t nGoodADC = 0;

  Int_t nErrorHC  = 0;
  Int_t nErrorMCM = 0;
  Int_t nErrorADC = 0;
  

  //Int_t adc_last  = -1;
  // Int_t sm_01  = -1;
  
  // number of ADCs per SM
  Int_t nADCinSM[kSM+1];
  for(Int_t sm=0; sm<kSM+1; sm++) nADCinSM[sm] = 0;



  while (data->Next()) {

    Int_t sm  = data->GetSM();
    Int_t layer = data->GetLayer();
    Int_t stack = data->GetStack();

    Int_t det = data->GetDet();
    Int_t side = data->GetSide();

    Int_t row = data->GetRow();
    Int_t col = data->GetCol();

    Int_t rob = data->GetROB();
    Int_t mcm = data->GetMCM();
    Int_t adc = data->GetADC();

 
    Int_t *sig = data->GetSignals();
    nb++;

    nADCinSM[sm]++;
    nADCinSM[kSM]++;
    
    // memory coruption protection
    if (det<0 || det>=kDET) continue;
    
    // check errors
   
    // tests
    //fErrorHC->Fill(data->GetHCErrorCode());
    if (data->GetMCMErrorCode() > 0) fErrorLocMCM[det]->Fill(row, col);
    if (data->GetADCErrorCode() > 0) fErrorLocADC[det]->Fill(row, col);    

    // new HC found
    if ((det + side*kDET) != (lastdet + lastside*kDET)) {
      Int_t code = data->GetHCErrorCode();
      // if (code) { 
      fErrorHC->Fill(code);
      
      if (code) fErrorSMHC->Fill(sm);
      if (code) nErrorHC++;
      nGoodHC++;

      //Int_t mask = 1;
      //for(Int_t cc = 0; cc < 3; cc++) {
      //  if (code & mask) fErrorHC->Fill(cc);
      //  cc *= 2;
      //	}
      //}
    }
    lastdet  = det;
    lastside = side;
    
    // new MCM found  
    if (mcm != lastmcm){
      Int_t code = data->GetMCMErrorCode();
      fErrorMCM->Fill(code);

      if (code) fErrorSMMCM->Fill(sm);
      if (code) nErrorMCM++;
      nGoodMCM++;
    }
    lastmcm = mcm;
    
    // new ADC channel found
    Int_t code = data->GetADCErrorCode();
    fErrorADC->Fill(code);
    if (code) fErrorSMADC->Fill(sm);
    if (code) nErrorADC++;
    nGoodADC++;

    // end of error checking
    
    // check the ROBs
    fDetRob->Fill(det, rob, 1./(kMCM*18));
    isUsed[det][row][col]++;

    // check if mcm signal is continuus
    if ((rob_last != rob) || (mcm_last != mcm)) {
      rob_last = rob;
      mcm_last = mcm;
      fnEntriesRM[det]->Fill(rob,mcm);
    }
    
    // number of entries for each channels
    fNPoint[det]->Fill(row, col);
    

    // create a structure for an MCM if needed
    Int_t mcmIndex = det * (kMCM * kROB) + rob * kMCM + mcm;
    if (fCreateFull && !fFullSignal[mcmIndex])
      fFullSignal[mcmIndex] = new TH2S(Form("mcm_%d_%d_%d_%d_%d", sm, stack, layer, rob, mcm), 
				       Form("mcm-%d-%d-%d-%d-%d;ADC;time bin", sm, stack, layer, rob, mcm),
				       21, -0.5, 20.5, 30, -0.5, 29.5);
    

    // loop over Time Bins and fill histograms
    Int_t minV = 1024;
    Int_t maxV = 0;

    for(Int_t k=0; k<kTB; k++) { /// to be corrected

      //if (data->GetADCErrorCode() > 0) continue;

      //if (col == 0 || col == 143)
      //printf("TB: %d %d %d\n", row, col, sig[k]);

      //if (sig[k] < 1) 
      //printf("det = %d rob = %d mcm = %d adc = %d k = %d S = %d\n", det, rob, mcm, adc, k, sig[k]);
      
      fSignal[det]->Fill(sig[k]);
      fData[det]->Fill(row, col, sig[k]);
      
      minV = (minV < sig[k]) ? minV : sig[k];
      maxV = (maxV > sig[k]) ? maxV : sig[k];


      // check if data strange enought
      if (fCreateFull && fFullSignal[mcmIndex]) {
	//if (sm == 17 && )
	//if (det != 29) {
	if (sig[k] > fThresh || sig[k] < 1) fFullCounter[mcmIndex]++;
	//if (sig[k] < 1) fFullCounter[mcmIndex] = 0; // remove austrian flag 
	//}
	fFullSignal[mcmIndex]->Fill(adc, k, sig[k]);
      }
      
      // noisy chamber
      if (det == 29 && col > 7) {
	fTBEvent->Fill(fnEvents, k, sig[k]);
      }
    }
    
    fPP->Fill(maxV-minV);
    fChPP[det]->Fill(maxV-minV);
    fSmPP[sm]->Fill(maxV-minV);
  }
  
  // is the dead-alive status changing during the run
  for(Int_t i=0; i<kDET; i++) {
    for(Int_t j=0; j<kCOL; j++)
      for(Int_t k=0; k<kPAD; k++)
  	fOccupancy->Fill(isUsed[i][j][k]);
  }

  // save interesting histos
  Int_t mcmTrackCandidate = 0;
  for(Int_t i = 0; i < kDET * kROB * kMCM; i++) { 
    if (fFullCounter[i] && fFullSignal[i] && CheckMCM(i) )  {
      
      fMcmTracks->AddLast(fFullSignal[i]->Clone(Form("event_%d_%s", fnEvents, fFullSignal[i]->GetName())));
      mcmTrackCandidate++;
      
      Int_t mcmTrackletDet = i/(kROB * kMCM); 	
      Int_t mcmTrackletMcm = i%(kROB * kMCM);
      fMapMCM->Fill(mcmTrackletDet, mcmTrackletMcm);
    }
  }
  
  fGraphMCM->SetPoint(fnEvents, fnEvents, mcmTrackCandidate);
  printf("Number of MCM track candidates = %d\n", mcmTrackCandidate);
  

  // update fraction of error graphs
  Double_t err;

  err = (nGoodHC > 0)? 100.*nErrorHC/nGoodHC : -1;
  fErrorGraphHC->SetPoint(fnEvents, fnEvents, err);
  
  err = (nGoodMCM > 0)? 100.*nErrorMCM/nGoodMCM : -1;
  fErrorGraphMCM->SetPoint(fnEvents, fnEvents, err);
  
  err = (nGoodADC > 0)? 100.*nErrorADC/nGoodADC : -1;
  fErrorGraphADC->SetPoint(fnEvents, fnEvents, err);

  // number of fired ADC per SM
  for(Int_t sm=0; sm<kSM+1; sm++) 
    fNumberADC[sm]->SetPoint(fnEvents, fnEvents, nADCinSM[sm]);


  fnEvents++;
  return nb;
}

///////////////////////////////////////////////////////////////////////////////////////////////////

void AliTRDqaBlackEvents::Process(const char *filename) 
{
  //
  // Process something
  //
  
  char fn[256];
  strcpy(fn, filename);
  
  //printf("FILENAME = %s (%s)\n", filename, fn);

  Int_t map[kDET];
  
  TH1D *hist = new TH1D("fitSignal", "", 50, -0.5, 49.5);
  TF1 *fit = new TF1("fit", "gaus(0)", 0, 20);
  fit->SetParameters(1e3, 10, 1);
    
  for(Int_t det=0; det<kDET; det++) {

    //printf("processing chamber %d\n", det);   

    map[det] = 0;
    if (fData[det]->GetSum() < 10) continue;
    map[det] = 1;

    // read reference distributions
    ReadRefHists(det);

    for(Int_t row=0; row<fData[det]->GetXaxis()->GetNbins(); row++) {
      for(Int_t pad=0; pad<fData[det]->GetYaxis()->GetNbins(); pad++) {
	
	// project the histogramm
	hist->Reset();
	for(Int_t bb=0; bb<50; bb++) {
	  Int_t dataBin = fData[det]->FindBin(row, pad, bb);
	  Double_t v = fData[det]->GetBinContent(dataBin);
	  hist->SetBinContent(bb+1, v);
	}

	Int_t bin = fChPed[det]->FindBin(row, pad);

	if (hist->GetSum() > 1) {
	  
	  Double_t ped = 0, noise = 0;

	  if (fFitType == 0) {
	    fit->SetParameters(1e3, 10, 1);
	    hist->Fit(fit, "q0", "goff", 0, 20);
	    TF1 *f = hist->GetFunction("fit");
	    ped = TMath::Abs(f->GetParameter(1));
	    noise = TMath::Abs(f->GetParameter(2));
	    fSmNoiseFit[det/30]->Fill(noise);
	  } else {
	    ped = hist->GetMean();
	    noise = hist->GetRMS();
	    fSmNoiseRms[det/30]->Fill(noise);
	    //if (pad == 0)
	    //  printf("data %f %f %f\n", hist->GetSum(), ped, noise);
	  }

	  fChPed[det]->SetBinContent(bin, ped);
	  fChNoise[det]->SetBinContent(bin, noise);
	  fNoiseTotal->Fill(noise);

	  // subtract reference values
	  Double_t refped = 0;
	  Double_t refnoise = 0;
	  
	  if (fRefHistPed)   refped   = fRefHistPed->GetBinContent(bin);
	  if (fRefHistNoise) refnoise = fRefHistPed->GetBinContent(bin);

	  fChPedRes[det]->SetBinContent(bin, ped-refped);
	  fChNoiseRes[det]->SetBinContent(bin, noise-refnoise);
	  
	  fPed[det]->Fill(ped);
	  fNoise[det]->Fill(noise);

	  // fill SM-HC plot
	  Int_t sm = det / 30;
	  Int_t hc = (pad < kPAD/2) ? 2* (det % 30) : 2* (det % 30) + 1;
	  if (ped > 9. && ped < 11) fSMHCped->Fill(sm, hc, 1./1152.); // number of pads in HC

	} else {
	  
	  // not enought data found 
	  fChPed[det]->SetBinContent(bin, 0);
	  fChNoise[det]->SetBinContent(bin, 0);
	  fChPedRes[det]->SetBinContent(bin, 0);
	  fChNoiseRes[det]->SetBinContent(bin, 0);
	}
	
	//delete hist;
      }
    }
  }


  //printf("Number of events = %d\n", fnEvents);

  // normalize number of entries histos
  Int_t max = 0;
  for(Int_t i=0; i<kDET; i++) { 
    if (!map[i]) continue;
    for(Int_t j=0; j<fNPoint[i]->GetXaxis()->GetNbins(); j++) {
      for(Int_t k=0; k<fNPoint[i]->GetYaxis()->GetNbins(); k++) {
	Int_t dataBin = fNPoint[i]->FindBin(j, k);
	Double_t v = fNPoint[i]->GetBinContent(dataBin);
	if (v > max) max = (Int_t)v;
      }
    }
  }
  
  char entriesDistName[100];

  for(Int_t i=0; i<kDET; i++) {
    
    if (!map[i]) continue;
    
    sprintf(entriesDistName, "entriesDist_%d", i);
    fNPointDist[i] = new TH1D(entriesDistName, ";number of events", max+2, -0.5, max+1.5);
    
    for(Int_t j=0; j<fNPoint[i]->GetXaxis()->GetNbins(); j++) {
      for(Int_t k=0; k<fNPoint[i]->GetYaxis()->GetNbins(); k++) {
	Int_t dataBin = fNPoint[i]->FindBin(j, k);
	Double_t v = fNPoint[i]->GetBinContent(dataBin);
	//if (v > fnEvents) printf("N = %d V = %lf\n", fnEvents, v);
	fNPointDist[i]->Fill(v); 
      }
    }
    
    fNPoint[i]->Scale(1./fnEvents);
  }
  

  for(Int_t i=0; i<kDET; i++) {
    fnEntriesRM[i]->SetMaximum(fnEvents * 1.5);
  }

  // save histograms

  //printf("FILENAME 2 = %s (%d)\n", fn, fn);
  TFile *file = new TFile(fn, "recreate");
  for(Int_t det = 0; det < kDET; det++) {
    if (!map[det]) continue; 
    fChPed[det]->Write();
    fChNoise[det]->Write();
    fNPoint[det]->Write();
    fNPointDist[det]->Write();
    fPed[det]->Write();
    fNoise[det]->Write();
    fSignal[det]->Write();
    fnEntriesRM[det]->Write();
    fChPP[det]->Write();

    fChPedRes[det]->Write();
    fChNoiseRes[det]->Write();

    // save error hists
    fErrorLocMCM[det]->SetMinimum(0);
    fErrorLocMCM[det]->SetMaximum(fnEvents);
    fErrorLocMCM[det]->Write();

    fErrorLocADC[det]->SetMinimum(0);
    fErrorLocADC[det]->SetMaximum(fnEvents);
    fErrorLocADC[det]->Write();
  }

  for(Int_t sm=0; sm<kSM; sm++) {
    fSmNoiseRms[sm]->Write();
    fSmNoiseFit[sm]->Write();
    fSmPP[sm]->Write();
  }



  Int_t nMcm = 0;
  for(Int_t i=0; i < kDET * kROB * kMCM; i++) {
    if (fFullSignal[i] && fFullCounter[i] > fCount) {
      fFullSignal[i]->Write();
      nMcm++;
    }
  }

  printf("Number of saved MCMs = %d\n", nMcm);
  
  fMcmTracks->Write();
  printf("Number of tracks = %d\n", fMcmTracks->GetEntries());
  
  // permanently problematic MCMs
  for(Int_t det=0; det<kDET; det++) {
    for(Int_t mcm=0; mcm<kROB*kMCM; mcm++) {
      
      Int_t mRob = mcm / kMCM;
      Int_t mMcm = mcm % kMCM;
      Int_t bin = fMapMCM->FindBin(det, mcm);
      Double_t frac = 1. * fMapMCM->GetBinContent(bin) / fnEvents;	
      fFracMCM->Fill(frac);
      
      if (frac > 0.7) {
	printf("{%d, %d, %d}, \n", det, mRob, mMcm, frac);
      }      
    }
  }



  fOccupancy->Write();
  fDetRob->Write();
  fTBEvent->Write();
  
  // error hists
  fErrorHC->Write();
  fErrorMCM->Write();
  fErrorADC->Write();

  fErrorSMHC->Write();
  fErrorSMMCM->Write();
  fErrorSMADC->Write();  
  
  // write graphs
  fErrorGraphHC->Write("trendErrorHC");
  fErrorGraphMCM->Write("trendErrorMCM");
  fErrorGraphADC->Write("trendErrorADC");
  
  fGraphMCM->Write("trendMCM");

  fMapMCM->SetMaximum(fnEvents);
  fMapMCM->Write();
  fFracMCM->Write();
  
  fSMHCped->Write();

  for(Int_t sm=0; sm<kSM; sm++)
    fNumberADC[sm]->Write(Form("nADCinSM%d",sm));
  
  fNumberADC[kSM]->Write("nADCinEvent");

  fNoiseTotal->Write();
  fPP->Write();

  file->Close();
  delete file;
}

///////////////////////////////////////////////////////////////////////////////////////////////////

Int_t AliTRDqaBlackEvents::CheckMCM(Int_t index) {
  
  return 1;
  
  static Int_t data[21][3] = {
    {1, 0, 1}, 
    {242, 0, 0}, 
    {242, 0, 1}, 
    {242, 0, 2}, 
    {242, 0, 4}, 
    {242, 0, 5}, 
    {242, 0, 6}, 
    {242, 0, 8}, 
    {242, 0, 12}, 
    {251, 7, 7}, 
    {254, 3, 11}, 
    {259, 3, 14}, 
    {260, 1, 9}, 
    {260, 3, 15}, 
    {273, 1, 7}, 
    {273, 1, 15}, 
    {276, 5, 11}, 
    {280, 6, 2}, 
    {299, 6, 4}, 
    {511, 2, 9}, 
    {517, 7, 15}
  };
  
  for(Int_t i=0; i<21; i++) {
    Int_t wIndex = data[i][0] * kROB*kMCM + data[i][1] * kMCM + data[i][2];
    if (index == wIndex) return 0;
  }

  return 1;
}

///////////////////////////////////////////////////////////////////////////////////////////////////







void AliTRDqaBlackEvents::DrawChamber(const char *filename, Int_t det, Int_t w, Int_t h) 
{
  //
  // Draw raport for one chamber: 
  // pedestal map, noise map, distribution of pedestal and noise
  // 
  // input:
  // name of the file with histograms (created with Process())
  // detector Id (0 - 539)
  // 

  // setup global style
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  gStyle->SetPadTopMargin(0.02);
  gStyle->SetPadBottomMargin(0.05);

  TFile *file = new TFile(filename, "READ");

  TCanvas *c = new TCanvas("blackEvents",Form("blackEvents %d",det), w, h);
  c->SetVertical(kFALSE);
  c->Divide(3,1, 0.01, 0.01);
  c->cd(3);
  
  TPad *mPad = (TPad*) gPad;
  mPad->Divide(1,2,0.01,0.01);
  
  c->cd(1);
  TH2D *h2 = (TH2D*)file->Get(Form("ped_%d",det));
  h2->SetMinimum(5);
  h2->SetMaximum(15);
  h2->SetTitle(";Z direction;#phi direction");
  h2->Draw("colz");
  
  c->cd(2);
  h2 = (TH2D*)file->Get(Form("noise_%d",det));
  h2->SetMinimum(fMinNoise);
  h2->SetMaximum(fMaxNoise);
  h2->SetTitle(";Z direction;#phi direction");
  h2->Draw("colz");
  
  mPad->cd(1);
  //gPad->SetLogy();
  TH1D *h1 = (TH1D*)file->Get(Form("pedDist_%d", det));
  h1->Draw();
  
  mPad->cd(2);
  gPad->SetLogy();
  h1 = (TH1D*)file->Get(Form("noiseDist_%d", det));
  h1->Draw();			 
  
  h1->Fit("gaus");
  TF1 *f = h1->GetFunction("gaus");
  const char *tt = Form("#mu = %.2f #sigma = %0.2f ", f->GetParameter(1),f->GetParameter(2));
  TLatex *ll = new TLatex(2, 100, tt);
  ll->SetTextSize(0.06);
  ll->Draw();
}

///////////////////////////////////////////////////////////////////////////////////////////////////

void AliTRDqaBlackEvents::DrawSm(const char *filename, Int_t sm, Int_t w, Int_t h) 
{
  //
  // ????????????
  //
  
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  
  gStyle->SetPadTopMargin(0.02);
  //gStyle->SetPadBottomMargin(0.05);  
  //gStyle->SetPadLeftMargin(0.02);  
  //gStyle->SetPadRightMargin(0.02);

  TFile *file = new TFile(filename, "READ");

  TCanvas *c = new TCanvas("blackEventsSM",Form("blackEvents SM %d",sm), w, h);
  c->SetVertical(kFALSE);
  c->Divide(5, 6, 0.001, 0.01);
  
  for(Int_t i=0; i<30; i++) {
    
    TH2D *h2 = (TH2D*)file->Get(Form("noise_%d",i+30*sm));
    if (!h2) continue;
    h2->SetMinimum(fMinNoise);
    h2->SetMaximum(fMaxNoise);

    // to be replaced by the official calculation
    Int_t stack = i/6;
    Int_t layer = i%6;
    Int_t index = (5-layer)*5 + stack + 1;
    //printf("%d %d %d %d\n", i, stack, layer, index);
    c->cd(index);
    gPad->SetBottomMargin(0.02);
    gPad->SetTopMargin(0.02);

    h2->Draw("col");
  }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
