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

///_________________________________________________________________________
///
/// This class constructs Digits out of Hits
///
///

// --- Standard library ---

// --- ROOT system ---
#include <TMath.h>
#include <TTree.h>
#include <TMap.h>
#include <TGeoManager.h>
#include <TGeoPhysicalNode.h>
#include <AliGeomManager.h>
#include <TRandom.h>
#include <TF1.h>
#include <TH1F.h>

// --- AliRoot header files ---
#include "AliRun.h"
#include "AliVZERO.h"
#include "AliVZEROhit.h"
#include "AliRunLoader.h"
#include "AliLoader.h"
#include "AliGRPObject.h"
#include "AliRunDigitizer.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliCDBEntry.h"
#include "AliVZEROCalibData.h"
#include "AliCTPTimeParams.h"
#include "AliVZEROdigit.h"
#include "AliVZERODigitizer.h"

ClassImp(AliVZERODigitizer)

 AliVZERODigitizer::AliVZERODigitizer()
                   :AliDigitizer(),
                    fCalibData(GetCalibData()),
                    fPhotoCathodeEfficiency(0.18),
                    fNdigits(0),
                    fDigits(0),
                    fSignalShape(NULL),
                    fPMResponse(NULL),
                    fSinglePhESpectrum(NULL)
   
{
  // default constructor

//    fNdigits = 0;
//    fDigits  = 0;
//   
//    fPhotoCathodeEfficiency =   0.18;
//    fPMVoltage              =  768.0;
//    fPMGain = TMath::Power((fPMVoltage / 112.5) ,7.04277); 
   
//   fCalibData = GetCalibData();

  fSignalShape = new TF1("VZEROSignalShape",this,&AliVZERODigitizer::SignalShape,0,200,6,"AliVZERODigitizer","SignalShape");
  //  fSignalShape->SetParameters(0,1.57345e1,-4.25603e-1,2.9,6.40982,3.69339e-01);
  //  fSignalShape->SetParameters(1.34330e+00,1.13007e+02,-4.95705e-01,
  //			      3.68911e+00,1.01040e+00, 3.94675e-01);
  fSignalShape->SetParameters(-1.07335e+00,2.16002e+01,-1.26133e-01,
			      1.41619e+00,5.50334e-01,3.86111e-01);

  fPMResponse = new TF1("VZEROPMResponse",this,&AliVZERODigitizer::PMResponse,-kPMRespTime,2.*kPMRespTime,0,"AliVZERODigitizer","PMResponse");
  fSinglePhESpectrum = new TF1("VZEROSinglePhESpectrum",this,&AliVZERODigitizer::SinglePhESpectrum,0,20,0,"AliVZERODigitizer","SinglePhESpectrum");

  // Now get the CTP L0->L1 delay
  AliCDBEntry *entry = AliCDBManager::Instance()->Get("GRP/CTP/CTPtiming");
  if (!entry) AliFatal("CTP timing parameters are not found in OCDB !");
  AliCTPTimeParams *ctpParams = (AliCTPTimeParams*)entry->GetObject();
  Float_t l1Delay = (Float_t)ctpParams->GetDelayL1L0()*25.0;

  AliCDBEntry *entry2 = AliCDBManager::Instance()->Get("VZERO/Calib/TimeDelays");
  if (!entry2) AliFatal("VZERO time delays are not found in OCDB !");
  TH1F *delays = (TH1F*)entry2->GetObject();

  for(Int_t i = 0 ; i < 64; ++i) {

    for(Int_t j = 0; j < kNClocks; ++j) fAdc[i][j] = 0;
    fLeadingTime[i] = fTimeWidth[i] = 0;

    fPmGain[i] = fCalibData->GetGain(i);

    fAdcPedestal[i][0] = fCalibData->GetPedestal(i);
    fAdcSigma[i][0]    = fCalibData->GetSigma(i); 
    fAdcPedestal[i][1] = fCalibData->GetPedestal(i+64);
    fAdcSigma[i][1]    = fCalibData->GetSigma(i+64); 

    Int_t board = AliVZEROCalibData::GetBoardNumber(i);
    fNBins[i] = TMath::Nint(((Float_t)(fCalibData->GetMatchWindow(board)+1)*25.0+
			     (Float_t)kMaxTDCWidth*fCalibData->GetWidthResolution(board))/
			    fCalibData->GetTimeResolution(board));
    fNBinsLT[i] = TMath::Nint(((Float_t)(fCalibData->GetMatchWindow(board)+1)*25.0)/
			      fCalibData->GetTimeResolution(board));
    fBinSize[i] = fCalibData->GetTimeResolution(board);
    fHptdcOffset[i] = (((Float_t)fCalibData->GetTriggerCountOffset(board)-
			(Float_t)fCalibData->GetRollOver(board))*25.0+
		       fCalibData->GetTimeOffset(i)+
		       l1Delay+
		       delays->GetBinContent(i+1)+
		       kV0Offset);

    fTime[i] = new Float_t[fNBins[i]];
    memset(fTime[i],0,fNBins[i]*sizeof(Float_t));
  }

}

//____________________________________________________________________________ 
  AliVZERODigitizer::AliVZERODigitizer(AliRunDigitizer* manager)
                    :AliDigitizer(manager),
		     fCalibData(GetCalibData()),
                     fPhotoCathodeEfficiency(0.18),
		     fNdigits(0),
                     fDigits(0),
		     fSignalShape(NULL),
                     fPMResponse(NULL),
                     fSinglePhESpectrum(NULL)
		                        
{
  // constructor
  // Initialize OCDB and containers used in the digitization
  
  fSignalShape = new TF1("VZEROSignalShape",this,&AliVZERODigitizer::SignalShape,0,200,6,"AliVZERODigitizer","SignalShape");
  //  fSignalShape->SetParameters(0,1.57345e1,-4.25603e-1,2.9,6.40982,3.69339e-01);
  //  fSignalShape->SetParameters(1.34330e+00,1.13007e+02,-4.95705e-01,
  //			      3.68911e+00,1.01040e+00, 3.94675e-01);
  fSignalShape->SetParameters(-1.07335e+00,2.16002e+01,-1.26133e-01,
			      1.41619e+00,5.50334e-01,3.86111e-01);
  fPMResponse = new TF1("VZEROPMResponse",this,&AliVZERODigitizer::PMResponse,-kPMRespTime,2.*kPMRespTime,0,"AliVZERODigitizer","PMResponse");
  fSinglePhESpectrum = new TF1("VZEROSinglePhESpectrum",this,&AliVZERODigitizer::SinglePhESpectrum,0,20,0,"AliVZERODigitizer","SinglePhESpectrum");
  
  // Now get the CTP L0->L1 delay
  AliCDBEntry *entry = AliCDBManager::Instance()->Get("GRP/CTP/CTPtiming");
  if (!entry) AliFatal("CTP timing parameters are not found in OCDB !");
  AliCTPTimeParams *ctpParams = (AliCTPTimeParams*)entry->GetObject();
  Float_t l1Delay = (Float_t)ctpParams->GetDelayL1L0()*25.0;

  AliCDBEntry *entry2 = AliCDBManager::Instance()->Get("VZERO/Calib/TimeDelays");
  if (!entry2) AliFatal("VZERO time delays are not found in OCDB !");
  TH1F *delays = (TH1F*)entry2->GetObject();

  for(Int_t i = 0 ; i < 64; ++i) {

    for(Int_t j = 0; j < kNClocks; ++j) fAdc[i][j] = 0;
    fLeadingTime[i] = fTimeWidth[i] = 0;

    fPmGain[i] = fCalibData->GetGain(i);

    fAdcPedestal[i][0] = fCalibData->GetPedestal(i);
    fAdcSigma[i][0]    = fCalibData->GetSigma(i); 
    fAdcPedestal[i][1] = fCalibData->GetPedestal(i+64);
    fAdcSigma[i][1]    = fCalibData->GetSigma(i+64); 

    Int_t board = AliVZEROCalibData::GetBoardNumber(i);
    fNBins[i] = TMath::Nint(((Float_t)(fCalibData->GetMatchWindow(board)+1)*25.0+
			     (Float_t)kMaxTDCWidth*fCalibData->GetWidthResolution(board))/
			    fCalibData->GetTimeResolution(board));
    fNBinsLT[i] = TMath::Nint(((Float_t)(fCalibData->GetMatchWindow(board)+1)*25.0)/
			      fCalibData->GetTimeResolution(board));
    fBinSize[i] = fCalibData->GetTimeResolution(board);
    fHptdcOffset[i] = (((Float_t)fCalibData->GetTriggerCountOffset(board)-
			(Float_t)fCalibData->GetRollOver(board))*25.0+
		       fCalibData->GetTimeOffset(i)+
		       l1Delay+
		       delays->GetBinContent(i+1)+
		       kV0Offset);

    fTime[i] = new Float_t[fNBins[i]];
    memset(fTime[i],0,fNBins[i]*sizeof(Float_t));
  }
}
           
//____________________________________________________________________________ 
  AliVZERODigitizer::~AliVZERODigitizer()
{
  // destructor
  
  if (fDigits) {
    fDigits->Delete();
    delete fDigits;
    fDigits=0; 
  }

  if (fSignalShape) {
    delete fSignalShape;
    fSignalShape = NULL;
  }
  if (fPMResponse) {
    delete fPMResponse;
    fPMResponse = NULL;
  }
  if (fSinglePhESpectrum) {
    delete fSinglePhESpectrum;
    fSinglePhESpectrum = NULL;
  }

  for(Int_t i = 0 ; i < 64; ++i) {
    if (fTime[i]) delete [] fTime[i];
  }
}

//_____________________________________________________________________________
Bool_t AliVZERODigitizer::Init()
{
  // Initialises the digitizer

  // Initialises the Digit array
  fDigits = new TClonesArray ("AliVZEROdigit", 1000);
  
  return kTRUE;
}

//____________________________________________________________________________
void AliVZERODigitizer::Exec(Option_t* /*option*/) 
{   
  // Creates digits from hits
  fNdigits     =    0;  

  Int_t labels[64][3];
  for(Int_t i = 0 ; i < 64; ++i) {
    memset(fTime[i],0,fNBins[i]*sizeof(Float_t));
    for(Int_t j = 0; j < kNClocks; ++j) fAdc[i][j] = 0;
    fLeadingTime[i] = fTimeWidth[i] = 0;
    labels[i][0] = labels[i][1] = labels[i][2] = -1;
  }
  Float_t integral = fPMResponse->Integral(-kPMRespTime,2.*kPMRespTime);
  Float_t meansPhE = fSinglePhESpectrum->Mean(0,20);

  AliRunLoader* outRunLoader =  AliRunLoader::GetRunLoader(fManager->GetOutputFolderName());    
  if (!outRunLoader) {
    Error("Exec", "Can not get output Run Loader");
    return;
  }
    
  AliLoader* outLoader = outRunLoader->GetLoader("VZEROLoader");

  if (!outLoader) {
    Error("Exec", "Can not get output VZERO Loader");
    return;
  }

  const char* mode = "update";
  if(outRunLoader->GetEventNumber() == 0) mode = "recreate";
  outLoader->LoadDigits(mode);

  if (!outLoader->TreeD()) outLoader->MakeTree("D");
  outLoader->MakeDigitsContainer();
  TTree* treeD  = outLoader->TreeD();
  Int_t bufsize = 16000;
  treeD->Branch("VZERODigit", &fDigits, bufsize); 
  
  for (Int_t iInput = 0; iInput < fManager->GetNinputs(); iInput++) {
     AliRunLoader* runLoader = AliRunLoader::GetRunLoader(fManager->GetInputFolderName(iInput));
     AliLoader* loader = runLoader->GetLoader("VZEROLoader");
     if (!loader) {
       Error("Exec", "Can not get VZERO Loader for input %d", iInput);
       continue;
	 }
      
     if (!runLoader->GetAliRun()) runLoader->LoadgAlice();

     AliVZERO* vzero = (AliVZERO*) runLoader->GetAliRun()->GetDetector("VZERO");
     if (!vzero) {
       Error("Exec", "No VZERO detector for input %d", iInput);
       continue;
	 }
      
     loader->LoadHits();
     TTree* treeH = loader->TreeH();
     if (!treeH) {
       Error("Exec", "Cannot get TreeH for input %d", iInput);
       continue; 
	 }
       
     TClonesArray* hits = vzero->Hits();

     //     Float_t lightYieldCorr[64] = {0.00707,0.00517,0.00520,0.00537,0.00735,0.00537,0.00733,0.00605,0.00778,0.00749,0.00701,0.00755,0.00732,0.00617,0.00669,0.00525,0.00752,0.00820,0.00797,0.01107,0.01080,0.00889,0.00880,0.01712,0.00866,0.00701,0.00811,0.00602,0.00879,0.00821,0.00861,0.01433,0.00061,0.00032,0.00099,0.00061,0.00034,0.00046,0.00031,0.00122,0.00155,0.00091,0.00032,0.00096,0.00120,0.00067,0.00113,0.00060,0.00158,0.00136,0.00340,0.00066,0.00076,0.00119,0.00129,0.00147,0.00137,0.00117,0.00088,0.00164,0.00128,0.00081,0.00121,0.00250};
     Float_t lightYieldCorr[64] = {0.01173,0.00874,0.00878,0.00886,0.01151,0.00925,0.01167,0.00983,0.01181,0.01243,0.01115,0.01220,0.01228,0.01053,0.01021,0.00930,0.01270,0.01411,0.01316,0.01894,0.01923,0.01860,0.01738,0.00305,0.01584,0.01251,0.01344,0.00310,0.01302,0.01266,0.01407,0.00338,0.00089,0.00100,0.00130,0.00081,0.00052,0.01230,0.00059,0.02452,0.02980,0.00137,0.01421,0.00116,0.00141,0.00092,0.02480,0.00096,0.00182,0.00174,0.00218,0.00106,0.00116,0.00160,0.00162,0.03097,0.00194,0.00171,0.00132,0.00239,0.00173,0.00118,0.00163,0.00262};
//  Now makes Digits from hits
     Int_t nTracks = (Int_t) treeH->GetEntries();
     for (Int_t iTrack = 0; iTrack < nTracks; iTrack++) {
         vzero->ResetHits();
         treeH->GetEvent(iTrack);
         Int_t nHits = hits->GetEntriesFast();
         for (Int_t iHit = 0; iHit < nHits; iHit++) {
	   AliVZEROhit* hit = (AliVZEROhit *)hits->UncheckedAt(iHit);
	   Int_t nPhot = hit->Nphot();
	   Int_t cell  = hit->Cell();                          
	   Int_t pmt = Cell2Pmt(cell);
	   Int_t trackLabel = hit->GetTrack();
	   for(Int_t l = 0; l < 3; ++l) {
	     if (labels[pmt][l] < 0) {
	       labels[pmt][l] = trackLabel;
	       break;
	     }
	   }
	   Float_t dt_scintillator = gRandom->Gaus(0,kIntTimeRes);
	   Float_t t = dt_scintillator + 1e9*hit->Tof();
	   if (pmt < 32) t += kV0CDelayCables;
	   t += fHptdcOffset[pmt];
	   Int_t nPhE;
	   Float_t prob = lightYieldCorr[pmt]*fPhotoCathodeEfficiency; // Optical losses included!
	   if (nPhot > 100)
	     nPhE = (Int_t)gRandom->Gaus(prob*Float_t(nPhot)+0.5,
					 sqrt(Float_t(nPhot)*prob*(1.-prob)));
	   else
	     nPhE = gRandom->Binomial(nPhot,prob);
	   Float_t charge = TMath::Qe()*fPmGain[pmt]*fBinSize[pmt]/integral;
	   for (Int_t iPhE = 0; iPhE < nPhE; ++iPhE) {
	     Float_t tPhE = t + fSignalShape->GetRandom(0,fBinSize[pmt]*Float_t(fNBins[pmt]));
	     Float_t gainVar = fSinglePhESpectrum->GetRandom(0,20)/meansPhE;
	     Int_t firstBin = TMath::Max(0,(Int_t)((tPhE-kPMRespTime)/fBinSize[pmt]));
	     Int_t lastBin = TMath::Min(fNBins[pmt]-1,(Int_t)((tPhE+2.*kPMRespTime)/fBinSize[pmt]));
	     for(Int_t iBin = firstBin; iBin <= lastBin; ++iBin) {
	       Float_t tempT = fBinSize[pmt]*(0.5+iBin)-tPhE;
	       fTime[pmt][iBin] += gainVar*charge*fPMResponse->Eval(tempT);
	     }
	   }         // ph.e. loop
         }           // hit loop
     }               // track loop
     loader->UnloadHits();
  }                  // input loop

  Float_t maximum = 0.9*fSignalShape->GetMaximum(0,200); // Not exact, one needs to do this on the convoluted
  Float_t integral2 = fSignalShape->Integral(0,200); // function. Anyway the effect is small <10% on the 2.5 ADC thr
  for (Int_t ipmt = 0; ipmt < 64; ++ipmt) {
    Float_t thr = fCalibData->GetDiscriThr(ipmt)*kChargePerADC*maximum*fBinSize[ipmt]/integral2;
    Bool_t ltFound = kFALSE, ttFound = kFALSE;
    for (Int_t iBin = 0; iBin < fNBins[ipmt]; ++iBin) {
      Float_t t = fBinSize[ipmt]*Float_t(iBin+0.5);
      if (fTime[ipmt][iBin] > thr) {
	if (!ltFound && (iBin < fNBinsLT[ipmt])) {
	  ltFound = kTRUE;
	  fLeadingTime[ipmt] = t;
	}
      }
      else {
	if (ltFound) {
	  if (!ttFound) {
	    ttFound = kTRUE;
	    fTimeWidth[ipmt] = t - fLeadingTime[ipmt];
	  }
	}
      }
      Float_t tadc = t - kClockOffset - fCalibData->GetTimeOffset(ipmt);
      Int_t clock = kNClocks/2 - Int_t(tadc/25.0);
      if (clock >= 0 && clock < kNClocks)
	fAdc[ipmt][clock] += fTime[ipmt][iBin]/kChargePerADC;
    }
    Int_t board = AliVZEROCalibData::GetBoardNumber(ipmt);
    if (ltFound && ttFound) {
      fTimeWidth[ipmt] = fCalibData->GetWidthResolution(board)*
	Float_t(Int_t(fTimeWidth[ipmt]/fCalibData->GetWidthResolution(board)));
      if (fTimeWidth[ipmt] < Float_t(kMinTDCWidth)*fCalibData->GetWidthResolution(board))
	fTimeWidth[ipmt] = Float_t(kMinTDCWidth)*fCalibData->GetWidthResolution(board);
      if (fTimeWidth[ipmt] > Float_t(kMaxTDCWidth)*fCalibData->GetWidthResolution(board))
	fTimeWidth[ipmt] = Float_t(kMaxTDCWidth)*fCalibData->GetWidthResolution(board);
    }
  }

  Int_t evenOrOdd = gRandom->Integer(2);
  for (Int_t j=0; j<64; ++j){
    for (Int_t iClock = 0; iClock < kNClocks; ++iClock) {
      Int_t integrator = (iClock + evenOrOdd) % 2;
      fAdc[j][iClock]  += gRandom->Gaus(fAdcPedestal[j][integrator], fAdcSigma[j][integrator]);
    }
  }
	
  // Now add digits to the digit Tree 

  Short_t *chargeADC = new Short_t[kNClocks];
  for (Int_t i=0; i<64; i++) {      
    Float_t totADC = 0;
    for (Int_t j = 0; j < kNClocks; ++j) {
      Int_t tempadc = Int_t(fAdc[i][j]);
      if (tempadc > 1023) tempadc = 1023;
      chargeADC[j] = tempadc;
      if (j >= 8 && j <= 11) {
	Int_t integrator = (j + evenOrOdd) % 2;
	if ((Float_t(tempadc) - fAdcPedestal[i][integrator]) > (2.*fAdcSigma[i][integrator]))
	  totADC += (Float_t(tempadc) - fAdcPedestal[i][integrator]);
      }
    }
    totADC += fAdcPedestal[i][(10+evenOrOdd)%2];
    AddDigit(i, totADC, fLeadingTime[i], fTimeWidth[i], Bool_t((10+evenOrOdd)%2), chargeADC, labels[i]);
  }
  delete [] chargeADC;

  treeD->Fill();
  outLoader->WriteDigits("OVERWRITE");  
  outLoader->UnloadDigits();     
  ResetDigit();
}

//____________________________________________________________________________
void AliVZERODigitizer::AddDigit(Int_t PMnumber, Float_t adc, Float_t time, Float_t width, Bool_t integrator, Short_t *chargeADC, Int_t *labels) 
 { 
 
// Adds Digit 
 
  TClonesArray &ldigits = *fDigits;  
	 
  new(ldigits[fNdigits++]) AliVZEROdigit(PMnumber,adc,time,width,kFALSE,kFALSE,integrator,chargeADC,labels);
	 
}
//____________________________________________________________________________
void AliVZERODigitizer::ResetDigit()
{

// Clears Digits

  fNdigits = 0;
  if (fDigits) fDigits->Delete();
}

//____________________________________________________________________________
AliVZEROCalibData* AliVZERODigitizer::GetCalibData() const

{
  AliCDBManager *man = AliCDBManager::Instance();

  AliCDBEntry *entry=0;

  entry = man->Get("VZERO/Calib/Data");

//   if(!entry){
//     AliWarning("Load of calibration data from default storage failed!");
//     AliWarning("Calibration data will be loaded from local storage ($ALICE_ROOT)");
//     Int_t runNumber = man->GetRun();
//     entry = man->GetStorage("local://$ALICE_ROOT/OCDB")
//       ->Get("VZERO/Calib/Data",runNumber);
// 	
//   }

  // Retrieval of data in directory VZERO/Calib/Data:


  AliVZEROCalibData *calibdata = 0;

  if (entry) calibdata = (AliVZEROCalibData*) entry->GetObject();
  if (!calibdata)  AliFatal("No calibration data from calibration database !");

  return calibdata;

}

double AliVZERODigitizer::SignalShape(double *x, double *par)
{
  // this function simulates the time
  // of arrival of the photons at the
  // photocathode
  Double_t xx = x[0];
  if (xx <= par[0]) return 0;
  Double_t a = 1./TMath::Power((xx-par[0])/par[1],1./par[2]);
  if (xx <= par[3]) return a;
  Double_t b = 1./TMath::Power((xx-par[3])/par[4],1./par[5]);
  Double_t f = a*b/(a+b);
  AliDebug(100,Form("x=%f func=%f",xx,f));
  return f;
}

double AliVZERODigitizer::PMResponse(double *x, double * /* par */)
{
  // this function describes the
  // PM time response to a single
  // photoelectron
  Double_t xx = x[0]+kPMRespTime;
  return xx*xx*TMath::Exp(-xx*xx/(kPMRespTime*kPMRespTime));
}

double AliVZERODigitizer::SinglePhESpectrum(double *x, double * /* par */)
{
  // this function describes the
  // PM amplitude response to a single
  // photoelectron
  Double_t xx = x[0];
  if (xx < 0) return 0;
  return (TMath::Poisson(xx,kPMNbOfSecElec)+kPMTransparency*TMath::Poisson(xx,1.0));
}

Int_t AliVZERODigitizer::Cell2Pmt(Int_t cell) const
{
  // The method maps the scintillator
  // indexes to the PM ones
  if (cell < 0 || cell >= 80) {
    AliError(Form("Wrong VZERO cell index %d",cell));
    return -1;
  }
  if (cell < 16) return cell;
  if (cell < 48) return 8 + cell/2;
  return cell - 16;
}
