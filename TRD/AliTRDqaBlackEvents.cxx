/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * withount fee, provided thats the abov copyright notice appears in all   *
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
#include "TF1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TGraph.h"

#include "AliLog.h"
#include "AliRawReader.h"

#include "AliTRDrawStreamOld.h"
#include "AliTRDqaBlackEvents.h"

ClassImp(AliTRDqaBlackEvents)

///////////////////////////////////////////////////////////////////////////////////////////////////

AliTRDqaBlackEvents::AliTRDqaBlackEvents() 
  :TObject() 
  ,fnEvents(0)
  ,fCreateFull(0)
  ,fThresh(0)
  ,fCount(0)
  ,fRefEv(0)
  ,fRefFileName(0x0)
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
{
  //
  // Constructor 
  // to create the histograms call Init()
  //

  for (Int_t i = 0; i < kDET; i++) {
    fPed[i]            = 0x0;
    fNoise[i]          = 0x0;
    fChPP[i]           = 0x0;
    fNPointDist[i]     = 0x0;
    fChPed[i]          = 0x0;
    fChNoise[i]        = 0x0;
    fNPoint[i]         = 0x0;
    fData[i]           = 0x0;
    fSignal[i]         = 0x0;
    fnEntriesRM[i]     = 0x0;
    fnEntriesRMDist[i] = 0x0;
    fChPedRes[i]       = 0x0;
    fChNoiseRes[i]     = 0x0;
    fErrorLocHC[i]     = 0x0;
    fErrorLocMCM[i]    = 0x0;
    fErrorLocADC[i]    = 0x0;
  }
  for (Int_t i = 0; i < 3; i++) {
    fGraphPP[i]        = 0x0;
    fSMLink[i]         = 0x0;
    fGrLink[i]         = 0x0;
    fppThresh[i]       = 0;
    fnPP[i]            = 0;
    fnLink[i]          = 0;
  }
  for (Int_t i = 0; i < 2; i++) {
    fnErrorHC[i]       = 0;
    fnErrorMCM[i]      = 0;
    fnErrorADC[i]      = 0;
  }
  for (Int_t i = 0; i < kSM; i++) {
    fSmNoiseRms[i]     = 0x0;
    fSmNoiseFit[i]     = 0x0;
    fSmPP[i]           = 0x0;
  }
  for (Int_t i = 0; i < kSM+1; i++) {
    fNumberADC[i]      = 0x0;
    fnADCinSM[i]       = 0;
  }
  for (Int_t i = 0; i < 1000; i++) {
    fEvNoDist[i]       = 0;
  }
  for (Int_t i = 0; i < kDET*kROB*kMCM; i++) {
    fFullSignal[i]     = 0x0;
    fFullCounter[i]    = 0;
  }
  //strncpy(fRefFileName,"",256);

}

///////////////////////////////////////////////////////////////////////////////////////////////////

AliTRDqaBlackEvents::AliTRDqaBlackEvents(const AliTRDqaBlackEvents &qa) 
  :TObject(qa) 
  ,fnEvents(0)
  ,fCreateFull(0)
  ,fThresh(0)
  ,fCount(0)
  ,fRefEv(0)
  ,fRefFileName(0x0)
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
{
  //
  // Copy constructor 
  // to create the histograms call Init()
  //
  
  for (Int_t i = 0; i < kDET; i++) {
    fPed[i]            = 0x0;
    fNoise[i]          = 0x0;
    fChPP[i]           = 0x0;
    fNPointDist[i]     = 0x0;
    fChPed[i]          = 0x0;
    fChNoise[i]        = 0x0;
    fNPoint[i]         = 0x0;
    fData[i]           = 0x0;
    fSignal[i]         = 0x0;
    fnEntriesRM[i]     = 0x0;
    fnEntriesRMDist[i] = 0x0;
    fChPedRes[i]       = 0x0;
    fChNoiseRes[i]     = 0x0;
    fErrorLocHC[i]     = 0x0;
    fErrorLocMCM[i]    = 0x0;
    fErrorLocADC[i]    = 0x0;
  }
  for (Int_t i = 0; i < 3; i++) {
    fGraphPP[i]        = 0x0;
    fSMLink[i]         = 0x0;
    fGrLink[i]         = 0x0;
    fppThresh[i]       = 0;
    fnPP[i]            = 0;
    fnLink[i]          = 0;
  }
  for (Int_t i = 0; i < 2; i++) {
    fnErrorHC[i]       = 0;
    fnErrorMCM[i]      = 0;
    fnErrorADC[i]      = 0;
  }
  for (Int_t i = 0; i < kSM; i++) {
    fSmNoiseRms[i]     = 0x0;
    fSmNoiseFit[i]     = 0x0;
    fSmPP[i]           = 0x0;
  }
  for (Int_t i = 0; i < kSM+1; i++) {
    fNumberADC[i]      = 0x0;
    fnADCinSM[i]       = 0;
  }
  for (Int_t i = 0; i < 1000; i++) {
    fEvNoDist[i]       = 0;
  }
  for (Int_t i = 0; i < kDET*kROB*kMCM; i++) {
    fFullSignal[i]     = 0x0;
    fFullCounter[i]    = 0;
  }
  //strncpy(fRefFileName,"",256);

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
    //fData[det]    = new TH3F(Form("data_%d", det), "", 16, -0.5, 15.5, 144, -0.5, 143.5, 50, -0.5, 49.5);

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
  fErrorHC  = new TH1D("errorHC", ";error ID;", 18, -3.5, 14.5);
  fErrorMCM = new TH1D("errorMCM", ";error ID;", 18, -3.5, 14.5);
  fErrorADC = new TH1D("errorADC", ";error ID;", 18, -3.5, 14.5);
  
  fErrorSMHC  = new TH1D("errorSM_HC", ";SM id", 18, -0.5, 17.5);
  fErrorSMMCM = new TH1D("errorSM_MCM", ";SM id", 18, -0.5, 17.5);
  fErrorSMADC = new TH1D("errorSM_ADC", ";SM id", 18, -0.5, 17.5);


  fErrorGraphHC  = new TGraph();
  fErrorGraphMCM = new TGraph();
  fErrorGraphADC = new TGraph();

  fGraphMCM = new TGraph();

  for(Int_t i=0; i<3; i++) {
    fGraphPP[i] = new TGraph();
  }


  fMapMCM = new TH2D("mapMCM", ";det;mcm", 540, -0.5, 539.5, kROB*kMCM, -0.5, kROB*kMCM-0.5);
  fFracMCM = new TH1D("fracMCM", ";frequency", 100, 0, 1);
  

  fErrorGraphHC->GetHistogram()->SetTitle("fraction of events with HC error;event number");
  fErrorGraphMCM->GetHistogram()->SetTitle("fraction of events with MCM error;event number;");
  fErrorGraphADC->GetHistogram()->SetTitle("fraction of events with ADC error;event number;"); 


  fSMHCped = new TH2D("smHcPed", ";super module;half chamber", 18, -0.5, 17.5, 60, -0.5, 59.5);
  
  // link monitor
  const char *linkName[3] = {"smLink", "smBeaf", "smData"};
  const char *linkGrName[3] = {"grSmLink", "grSmBeaf", "grSmData"};
  for(Int_t i=0; i<3; i++) {
    fSMLink[i] = new TH2D(linkName[i], ";super module;link", 18, -0.5, 17.5, 60, -0.5, 59.5);
    fGrLink[i] = new TGraph();
    fGrLink[i]->SetName(linkGrName[i]);
  }

  //fZSsize = new TH1D("zssizeSingle", ";threshold;nADC", 40, -0.5, 39.5);
  

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

  // event number consistancy
  for(Int_t i=0; i<1000; i++) {
    fEvNoDist[i] = new TH1D(Form("mcmEvDist_%d", i), ";#Delta Events", 201, -100.5, 100.5);
  }
  fRefEv = -1;

  fMcmTracks = new TObjArray();

  // clean data direct
  for(Int_t det=0; det<kDET; det++) 
    for(Int_t row=0; row<kROW; row++)
      for(Int_t pad=0; pad<kPAD; pad++)
	for(Int_t ch=0; ch<kCH; ch++) {
	  fDataDirect[det][row][pad][ch] = 0;
	  fSignalDirect[det][ch] = 0;  // overdone
	}
}


///////////////////////////////////////////////////////////////////////////////////////////////////

void AliTRDqaBlackEvents::Reset() 
{
  //
  // Resets the histograms
  //

  for(Int_t i=0; i<kDET; i++) {
    //fData[i]->Reset();
    fChPed[i]->Reset();
    fChNoise[i]->Reset();
  }
}


///////////////////////////////////////////////////////////////////////////////////////////////////

void AliTRDqaBlackEvents::SetRefFile(const char *filename) {
  
  //strncpy(fRefFileName,filename,256);
  fRefFileName = filename;

}

///////////////////////////////////////////////////////////////////////////////////////////////////

void AliTRDqaBlackEvents::ReadRefHists(Int_t det) {
  //  
  // Read the reference histograms
  //

  fRefHistPed = 0;
  fRefHistNoise = 0;
  
  TFile *file = 0x0;
  if (fRefFileName) file = TFile::Open(fRefFileName);
  if (!file) return;

  fRefHistPed   = (TH2D*)file->Get(Form("ped_%d",det));
  fRefHistNoise = (TH2D*)file->Get(Form("noise_%d", det));

  if (file) file->Close();
}

///////////////////////////////////////////////////////////////////////////////////////////////////

void AliTRDqaBlackEvents::StartEvent()
{
  //
  // start an event
  //

  // clear the mcm data
  for(Int_t i=0; i < kDET * kROB * kMCM; i++) {
    if (fFullSignal[i]) fFullSignal[i]->Reset();
    fFullCounter[i] = 0;
  }

  for(Int_t i=0; i<2; i++) {
    fnErrorHC[i] = 0;
    fnErrorMCM[i] = 0;
    fnErrorADC[i] = 0;
  }

 
  Int_t ppThresh[3] = {10, 20, 40};
  for(Int_t i=0; i<3; i++) {
    fppThresh[i] = ppThresh[i];
    fnPP[i] = 0;
    fnLink[i] = 0;
  }

  for(Int_t sm=0; sm<kSM+1; sm++) fnADCinSM[sm] = 0;

  if (fRefEv > 0) fRefEv++;
  fEvNoDist[999]->Reset();  // keep only the last event

}

///////////////////////////////////////////////////////////////////////////////////////////////////

void AliTRDqaBlackEvents::AddBuffer(AliTRDrawStreamOld *data, AliRawReader * const reader) 
{
  
  //printf ("try to read data\n");
  Int_t nextBuff  = data->NextBuffer();
  //printf("done ...\n");

  if (nextBuff == 0) return;
 
  Int_t sm = reader->GetEquipmentId() - 1024;
  //printf("reading SM %d\n", sm);
  AliInfo(Form("reading SM %d", sm));
  
  if (sm < 0 || sm > 17) return;

  // lopp over stacks, links ...

  for (Int_t istack = 0; istack < 5; istack++) {	
    for (Int_t ilink = 0; ilink < 12; ilink++) {
      
      //printf("HC = %d %d\n", istack, ilink);

      Int_t det = sm * 30 + istack * 6 + ilink/2;	
      
      // check if data delivered
      if (!(data->IsLinkActiveInStack(istack, ilink))) continue;
      fSMLink[0]->Fill(sm, istack * 12 + ilink);
      fnLink[0]++;

      // check if beaf-beaf
      if (data->GetLinkMonitorError(istack, ilink)) {
	fSMLink[1]->Fill(sm, istack * 12 + ilink);
	fnLink[1]++;
	continue;
      }
	
      // fill histogram with HC header errors
      Int_t nErrHc = 0;
      Int_t nErrHcTot = 0;
      
      nErrHc = FillBits(fErrorHC, data->GetH0ErrorCode(istack, ilink), 0);
      if (!nErrHc) fErrorHC->Fill(-3);
      nErrHcTot += nErrHc;
      
      nErrHc = FillBits(fErrorHC, data->GetH1ErrorCode(istack, ilink), 2);
      if (!nErrHc) fErrorHC->Fill(-2);
      nErrHcTot += nErrHc;
      
      nErrHc = FillBits(fErrorHC, data->GetHCErrorCode(istack, ilink), 4);
      if (!nErrHc) fErrorHC->Fill(-1);
      nErrHcTot += nErrHc;
      
      // trending
      fnErrorHC[0]++;
      if (nErrHcTot > 0) { 
	fnErrorHC[1]++;
	fErrorSMHC->Fill(sm);
      }
      
      // data integrity protection 

      //if (data->GetHCErrorCode(istack, ilink) > 0) continue;
      if (data->GetH0ErrorCode(istack, ilink) > 0) continue;
      if (data->GetH1ErrorCode(istack, ilink) > 0) continue;
      
      fSMLink[2]->Fill(sm, istack * 12 + ilink);	
      fnLink[2]++;
             
      
      for (Int_t imcm = 0; imcm < data->GetHCMCMmax(istack, ilink); imcm++ ){
	  
	//printf("mcm = %d %d %d\n", istack, ilink, imcm);
	
	// fill MCM error code
	
	Int_t nErrMcm = 0;
	Int_t nErrMcmTot = 0;
	
	nErrMcm = FillBits(fErrorMCM, data->GetMCMhdErrorCode(istack, ilink, imcm), 0);
	if (!nErrMcm) fErrorMCM->Fill(-3);
	nErrMcmTot += nErrMcm;
	
	nErrMcm = FillBits(fErrorMCM, data->GetMCMADCMaskErrorCode(istack, ilink, imcm), 5);
	if (!nErrMcm) fErrorMCM->Fill(-2);
	nErrMcmTot += nErrMcm;
	
	nErrMcm = FillBits(fErrorMCM, data->GetMCMErrorCode(istack, ilink, imcm), 10);
	if (!nErrMcm) fErrorMCM->Fill(-1);
	nErrMcmTot += nErrMcm;			  
	
	// trending
	fnErrorMCM[0]++;
	if (nErrMcmTot > 0) { 
	  fnErrorMCM[1]++;
	  fErrorSMMCM->Fill(sm);
	}
	
	// MCM protection
	if ( (data->GetMCMhdErrorCode(istack,ilink,imcm)) & 2 ) continue;
	//if ((data->GetMCMADCMaskErrorCode(istack,ilink,imcm))) continue;
	//if ((data->GetMCMErrorCode(istack,ilink,imcm))) continue;
	
	Int_t mcmEvent = data->GetEventNumber(istack, ilink, imcm);
	
	// set the reference event number 
	if (fRefEv < 0) {
	  fRefEv = mcmEvent;
	  printf("Reference Event Number = %d (%d %d %d)\n", fRefEv, istack, ilink, imcm);
	}
	
	// fill event distribution
	if (!(fnEvents%10)) {
	  fEvNoDist[fnEvents/10]->Fill(mcmEvent - fRefEv);
	}
	
	fEvNoDist[999]->Fill(mcmEvent - fRefEv);
	
	Int_t mcm = data->GetMCM(istack, ilink, imcm);
	Int_t rob = data->GetROB(istack, ilink, imcm);

	// create a structure for an MCM if needed
	Int_t mcmIndex = det * (kMCM * kROB) + rob * kMCM + mcm;
	if (fCreateFull && !fFullSignal[mcmIndex])
	  fFullSignal[mcmIndex] = 
	    new TH2S(Form("mcm_%d_%d_%d_%d_%d", sm, istack, ilink/2, rob, mcm), 
		     Form("mcm-%d-%d-%d-%d-%d;ADC;time bin", sm, istack, ilink/2, rob, mcm),
		     21, -0.5, 20.5, 30, -0.5, 29.5);
	
	
	//Int_t zsADC[21][40];
	/*
	  for(Int_t ina=0; ina<21; ina++) 
	  for(Int_t th=0; th<40; th++)
	  zsADC[ina][th] = 0;
	*/
	
	// first loop over ADC chanels 
	
	for (Int_t iadc=0; iadc < data->GetADCcount(istack, ilink, imcm); iadc++) {
	  
	  //printf("ADC = %d\n", iadc);
	  

	  // fill ADC error bits
	  Int_t nErrAdc = FillBits(fErrorADC, data->GetADCErrorCode(), 0);
	  if (!nErrAdc) fErrorADC->Fill(-1);
	  
	  fnErrorADC[0]++;
	  if (nErrAdc > 0) {
	    fnErrorADC[1]++;
	    fErrorSMADC->Fill(sm);
	  }
	  
	  // ADC protection
	  if ((data->GetADCErrorCode(istack,ilink,imcm,iadc))) continue;

	  Int_t minV = 1024;
	  Int_t maxV = 0;
	    
	  Int_t *sig = data->GetSignalDirect(istack, ilink, imcm, iadc);

	  //Int_t adc = data->GetADCnumber(istack, ilink, imcm, iadc);
	  Int_t row = data->GetRow(istack, ilink, imcm);
	  Int_t col = data->GetCol(istack, ilink, imcm, iadc);	    	    
	    
	  // loop over Time Bins and fill histograms
	  for(Int_t k=0; k < data->GetNumberOfTimeBins(istack, ilink); k++) { 
	      	      
	    //fSignal[det]->Fill(sig[k]);
	    //fData[det]->Fill(row, col, sig[k]); // slow

	    if ((sig[k] >=0) && (sig[k] < kCH)) {
	      fSignalDirect[det][sig[k]]++;
	      fDataDirect[det][row][col][sig[k]]++; // direct data
	    }

	    // peak-peak
	    minV = (minV < sig[k]) ? minV : sig[k];
	    maxV = (maxV > sig[k]) ? maxV : sig[k];
	    
	    // check for active MCMs
	    if (fCreateFull && fFullSignal[mcmIndex]) {
	      if (sig[k] > fThresh || sig[k] < 0) fFullCounter[mcmIndex]++;
	      //if (sm == 0 && istack == 0 && ilink/2 == 1 && rob == 1 && mcm == 15) fFullCounter[mcmIndex]++; // special
	      //fFullSignal[mcmIndex]->Fill(adc, k, sig[k]); // slow
	    }
	    
	    // zero suppresion tests
	    /*
	      for(Int_t th=0; th<40; th++) {
	      if (sig[k] > th) {
	      zsADC[iadc][th] = 1;
	      if (iadc > 0) zsADC[iadc-1][th] = 1;
	      if (iadc < 10) zsADC[iadc+1][th] = 1;
	      }
	      }
	    */
	    
	  } // tb
	  
	  if (maxV > 0) { 
	    fnADCinSM[sm]++;
	    fnADCinSM[kSM]++;
	  }
	  
	  Int_t adcPP = maxV - minV;
	  //if (adcPP == 100) fFullCounter[mcmIndex] += 10;

	  fPP->Fill(adcPP);
	  fChPP[det]->Fill(adcPP);
	  fSmPP[sm]->Fill(adcPP);
	  
	  for(Int_t i=0; i<3; i++) {
	    if ((adcPP) > fppThresh[i]) fnPP[i]++;
	  }
	  
	    
	} // adc

	  // fill ZS histos
	  /*
	    for(Int_t th=0; th<40; th++) {
	    Int_t nnADC = 0;
	    for(Int_t ins=0; ins<21; ins++) 
	    nnADC += zsADC[ins][th];
	    fZSsize->Fill(th, nnADC);
	    }
	  */

	// fill active MCMs

	if (fCreateFull && fFullSignal[mcmIndex] && (fFullCounter[mcmIndex] > fCount)) {
	  
	  for (Int_t iadc=0; iadc < data->GetADCcount(istack, ilink, imcm); iadc++) {
	    
	    // ADC protection
	    if ((data->GetADCErrorCode(istack,ilink,imcm,iadc))) continue;	  
	    
	    //Int_t row = data->GetRow(istack, ilink, imcm);
	    //Int_t col = data->GetCol(istack, ilink, imcm, iadc);	    	    
	    Int_t adc = data->GetADCnumber(istack, ilink, imcm, iadc);
	    
	    Int_t *sig = data->GetSignalDirect(istack, ilink, imcm, iadc);
	    
	    // loop over Time Bins and fill histograms
	    for(Int_t k=0; k < data->GetNumberOfTimeBins(istack, ilink); k++) { 
	      fFullSignal[mcmIndex]->Fill(adc, k, sig[k]); // slow
	    }
	  } // tb
	}
	
      } // mcm 
    } // link
  } // stack  
  
  // printf("end of loops\n");
} 


///////////////////////////////////////////////////////////////////////////////////////////////////

void AliTRDqaBlackEvents::FinishEvent()
{
  //
  // Processing at the end of the current event
  //

  for(Int_t i=0; i<3; i++) {
    fGraphPP[i]->SetPoint(fnEvents, fnEvents, fnPP[i]);
  }

  // trend of the number of links
  for(Int_t i=0; i<3; i++) {
    fGrLink[i]->SetPoint(fnEvents, fnEvents, fnLink[i]);
  }

  // save interesting histos
  Int_t mcmTrackCandidate = 0;
  for(Int_t i = 0; i < kDET * kROB * kMCM; i++) { 
    if ((fFullCounter[i] > fCount) && fFullSignal[i] && CheckMCM(i) )  {
      
      fMcmTracks->AddLast(fFullSignal[i]->Clone(Form("event_%d_%s", fnEvents, fFullSignal[i]->GetName())));
      mcmTrackCandidate++;
      
      Int_t mcmTrackletDet = i/(kROB * kMCM); 	
      Int_t mcmTrackletMcm = i%(kROB * kMCM);
      fMapMCM->Fill(mcmTrackletDet, mcmTrackletMcm);
    }
  }
  
  fGraphMCM->SetPoint(fnEvents, fnEvents, mcmTrackCandidate);
  AliInfo(Form("Number of MCM track candidates = %d\n", mcmTrackCandidate));
  
  
  // update fraction of error graphs
  Double_t err;
  
  err = (fnErrorHC[0] > 0)? 100.*fnErrorHC[1]/fnErrorHC[0] : -1;
  fErrorGraphHC->SetPoint(fnEvents, fnEvents, err);
  
  err = (fnErrorMCM[0] > 0)? 100.*fnErrorMCM[1]/fnErrorMCM[0] : -1;
  fErrorGraphMCM->SetPoint(fnEvents, fnEvents, err);
  
  err = (fnErrorADC[0] > 0)? 100.*fnErrorADC[1]/fnErrorADC[0] : -1;
  fErrorGraphADC->SetPoint(fnEvents, fnEvents, err);

  // number of fired ADC per SM
  for(Int_t sm=0; sm<kSM+1; sm++) 
    fNumberADC[sm]->SetPoint(fnEvents, fnEvents, fnADCinSM[sm]);

  fnEvents++;
}

///////////////////////////////////////////////////////////////////////////////////////////////////

void AliTRDqaBlackEvents::Process(const char *filename) 
{
  //
  // Process something
  //
  
  //char fn[256];
  //strncpy(fn,filename,256);
  
  //AliInfo(Form("FILENAME = %s (%s)\n", filename, fn));

  Int_t map[kDET];
  
  TH1D *hist = new TH1D("fitSignal", "", 50, -0.5, 49.5);
  TF1 *fit = new TF1("fit", "gaus(0)", 0, 20);
  fit->SetParameters(1e3, 10, 1);
    
  for(Int_t det=0; det<kDET; det++) {

    //AliInfo(Form("processing chamber %d\n", det));   

    map[det] = 0;
    //if (fData[det]->GetSum() < 10) continue;
    //if (fDataDirect[det][10][10][10] < 20) continue;
    //map[det] = 1;


    // rewrite signal-direct
    for(Int_t ch=0; ch<kCH; ch++) {
      fSignal[det]->Fill(ch, fSignalDirect[det][ch]);
    }

    // read reference distributions
    ReadRefHists(det);

    //for(Int_t row=0; row<fData[det]->GetXaxis()->GetNbins(); row++) {
    //for(Int_t pad=0; pad<fData[det]->GetYaxis()->GetNbins(); pad++) {
	
    for(Int_t row=0; row<kROW; row++) {
      for(Int_t pad=0; pad<kPAD; pad++) {

	// project the histogramm
	hist->Reset();
	//for(Int_t bb=0; bb<50; bb++) {
	for(Int_t bb=0; bb<kCH; bb++) {
	  //Int_t dataBin = fData[det]->FindBin(row, pad, bb);
	  //Double_t v = fData[det]->GetBinContent(dataBin);
	  hist->SetBinContent(bb+1, fDataDirect[det][row][pad][bb]);
	}

	Int_t bin = fChPed[det]->FindBin(row, pad);

	if (hist->GetSum() > 1) {
	  
	  map[det] = 1;
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
	    //  AliInfo(Form("data %f %f %f\n", hist->GetSum(), ped, noise));
	  }

	  fChPed[det]->SetBinContent(bin, ped);
	  fChNoise[det]->SetBinContent(bin, noise);
	  fNoiseTotal->Fill(noise);

	  // subtract reference values
	  Double_t refped = 0;
	  Double_t refnoise = 0;
	  
	  if (fRefHistPed)   refped   = fRefHistPed->GetBinContent(bin);
	  if (fRefHistPed)   refnoise = fRefHistPed->GetBinContent(bin);
	  // Original code, should it not be fRefHistNoise->GetBinContent(bin)
	  // instead of fRefHistPed->GetBinContent(bin)  (CBL) ???
	  //if (fRefHistNoise) refnoise = fRefHistPed->GetBinContent(bin);

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


  //AliInfo(Form("Number of events = %d\n", fnEvents));

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
    
    snprintf(entriesDistName,100,"entriesDist_%d",i);
    fNPointDist[i] = new TH1D(entriesDistName, ";number of events", max+2, -0.5, max+1.5);
    
    for(Int_t j=0; j<fNPoint[i]->GetXaxis()->GetNbins(); j++) {
      for(Int_t k=0; k<fNPoint[i]->GetYaxis()->GetNbins(); k++) {
	Int_t dataBin = fNPoint[i]->FindBin(j, k);
	Double_t v = fNPoint[i]->GetBinContent(dataBin);
	//if (v > fnEvents) AliInfo(Form("N = %d V = %lf\n", fnEvents, v));
	fNPointDist[i]->Fill(v); 
      }
    }
    
    fNPoint[i]->Scale(1./fnEvents);
  }
  

  for(Int_t i=0; i<kDET; i++) {
    fnEntriesRM[i]->SetMaximum(fnEvents * 1.5);
  }

  // save histograms

  //AliInfo(Form("FILENAME 2 = %s (%d)\n", fn, fn));
  TFile *file = new TFile(filename, "recreate");
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

  AliInfo(Form("Number of saved MCMs = %d\n", nMcm));
  
  fMcmTracks->Write();
  AliInfo(Form("Number of tracks = %d\n", fMcmTracks->GetEntries()));
  
  // permanently problematic MCMs
  for(Int_t det=0; det<kDET; det++) {
    for(Int_t mcm=0; mcm<kROB*kMCM; mcm++) {
      
      Int_t mRob = mcm / kMCM;
      Int_t mMcm = mcm % kMCM;
      Int_t bin = fMapMCM->FindBin(det, mcm);
      Double_t frac = 1. * fMapMCM->GetBinContent(bin) / fnEvents;	
      fFracMCM->Fill(frac);
      
      if (frac > 0.7) {
	AliInfo(Form("{%d, %d, %d, %f}, \n", det, mRob, mMcm, frac));
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
  
  for(Int_t i=0; i<3; i++) {
    fGraphPP[i]->Write(Form("fracPP_%d", i));
  }
  
  //fZSsize->Scale(1./fnEvents);
  //fZSsize->Write();

  fMapMCM->SetMaximum(fnEvents);
  fMapMCM->Write();
  fFracMCM->Write();
  
  fSMHCped->Write();

  for(Int_t i=0; i<3; i++ ) {
    fSMLink[i]->Write();
    fGrLink[i]->Write();
  }

  for(Int_t sm=0; sm<kSM; sm++)
    fNumberADC[sm]->Write(Form("nADCinSM%d",sm));
  
  fNumberADC[kSM]->Write("nADCinEvent");

  fNoiseTotal->Write();
  fPP->Write();

  for(Int_t i=0; i<1000; i++) {
    if (fEvNoDist[i]->GetSum() > 0) fEvNoDist[i]->Write();
  }


  file->Close();
  delete file;
}

///////////////////////////////////////////////////////////////////////////////////////////////////

Int_t AliTRDqaBlackEvents::CheckMCM(Int_t /*index*/) const {
  //
  // Checks a single MCM
  //  

  return 1;
  
  // static Int_t data[21][3] = {
  //   {1, 0, 1}, 
  //   {242, 0, 0}, 
  //   {242, 0, 1}, 
  //   {242, 0, 2}, 
  //   {242, 0, 4}, 
  //   {242, 0, 5}, 
  //   {242, 0, 6}, 
  //   {242, 0, 8}, 
  //   {242, 0, 12}, 
  //   {251, 7, 7}, 
  //   {254, 3, 11}, 
  //   {259, 3, 14}, 
  //   {260, 1, 9}, 
  //   {260, 3, 15}, 
  //   {273, 1, 7}, 
  //   {273, 1, 15}, 
  //   {276, 5, 11}, 
  //   {280, 6, 2}, 
  //   {299, 6, 4}, 
  //   {511, 2, 9}, 
  //   {517, 7, 15}
  // };
  
  // for(Int_t i=0; i<21; i++) {
  //   Int_t wIndex = data[i][0] * kROB*kMCM + data[i][1] * kMCM + data[i][2];
  //   if (index == wIndex) return 0;
  // }

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
    //AliInfo(Form("%d %d %d %d\n", i, stack, layer, index));
    c->cd(index);
    gPad->SetBottomMargin(0.02);
    gPad->SetTopMargin(0.02);

    h2->Draw("col");
  }
}

///////////////////////////////////////////////////////////////////////////////////////////////////

Int_t AliTRDqaBlackEvents::FillBits(TH1D *hist, Int_t code, Int_t offset) {
  //
  // Fill bits
  //

  Int_t nb = 0;
  UInt_t test = 1;
  for(Int_t i=0; i<8; i++) {
    if (code & test) {
      hist->Fill(i+offset);
      nb++;
    }
    test *= 2;       
  }
  
  return nb;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
