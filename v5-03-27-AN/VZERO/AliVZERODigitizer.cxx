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
#include "AliDigitizationInput.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliCDBEntry.h"
#include "AliVZEROCalibData.h"
#include "AliCTPTimeParams.h"
#include "AliLHCClockPhase.h"
#include "AliVZEROdigit.h"
#include "AliVZERODigitizer.h"
#include "AliVZEROSDigit.h"

ClassImp(AliVZERODigitizer)

 AliVZERODigitizer::AliVZERODigitizer()
                   :AliDigitizer(),
                    fCalibData(GetCalibData()),
                    fPhotoCathodeEfficiency(0.18),
                    fNdigits(0),
                    fDigits(0),
                    fSignalShape(NULL),
                    fPMResponse(NULL),
                    fSinglePhESpectrum(NULL),
		    fEvenOrOdd(kFALSE),
		    fTask(kHits2Digits),
		    fVZERO(NULL)
{
  // default constructor
  // Initialize OCDB and containers used in the digitization

  Init();
}

//____________________________________________________________________________ 
  AliVZERODigitizer::AliVZERODigitizer(AliVZERO *vzero, DigiTask_t task)
                    :AliDigitizer(),
		     fCalibData(GetCalibData()),
                     fPhotoCathodeEfficiency(0.18),
		     fNdigits(0),
                     fDigits(0),
		     fSignalShape(NULL),
                     fPMResponse(NULL),
                     fSinglePhESpectrum(NULL),
		     fEvenOrOdd(kFALSE),
		     fTask(task),
		     fVZERO(vzero)
{
  // constructor
  // Initialize OCDB and containers used in the digitization

  Init();
}
           
//____________________________________________________________________________ 
  AliVZERODigitizer::AliVZERODigitizer(AliDigitizationInput* digInput)
                    :AliDigitizer(digInput),
		     fCalibData(GetCalibData()),
                     fPhotoCathodeEfficiency(0.18),
		     fNdigits(0),
                     fDigits(0),
		     fSignalShape(NULL),
                     fPMResponse(NULL),
                     fSinglePhESpectrum(NULL),
		     fEvenOrOdd(kFALSE),
		     fTask(kHits2Digits),
		     fVZERO(NULL)
{
  // constructor
  // Initialize OCDB and containers used in the digitization

  Init();
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
  // Initialize OCDB and containers used in the digitization

  // check if the digitizer was already initialized
  if (fSignalShape) return kTRUE;

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

  AliCDBEntry *entry1 = AliCDBManager::Instance()->Get("GRP/CTP/TimeAlign");
  if (!entry1) AliFatal("CTP time-alignment is not found in OCDB !");
  AliCTPTimeParams *ctpTimeAlign = (AliCTPTimeParams*)entry1->GetObject();
  l1Delay += ((Float_t)ctpTimeAlign->GetDelayL1L0()*25.0);

  AliCDBEntry *entry2 = AliCDBManager::Instance()->Get("VZERO/Calib/TimeDelays");
  if (!entry2) AliFatal("VZERO time delays are not found in OCDB !");
  TH1F *delays = (TH1F*)entry2->GetObject();

  AliCDBEntry *entry3 = AliCDBManager::Instance()->Get("GRP/Calib/LHCClockPhase");
  if (!entry3) AliFatal("LHC clock-phase shift is not found in OCDB !");
  AliLHCClockPhase *phase = (AliLHCClockPhase*)entry3->GetObject();

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
    fHptdcOffset[i] = (((Float_t)fCalibData->GetRollOver(board)-
			(Float_t)fCalibData->GetTriggerCountOffset(board))*25.0+
		       fCalibData->GetTimeOffset(i)-
		       l1Delay-
		       phase->GetMeanPhase()+
		       delays->GetBinContent(i+1)+
		       kV0Offset);
    fClockOffset[i] = (((Float_t)fCalibData->GetRollOver(board)-
			(Float_t)fCalibData->GetTriggerCountOffset(board))*25.0+
		       fCalibData->GetTimeOffset(i)-
		       l1Delay+
		       kV0Offset);

    fTime[i] = new Float_t[fNBins[i]];
    memset(fTime[i],0,fNBins[i]*sizeof(Float_t));
  }

  return kTRUE;
}

//____________________________________________________________________________
void AliVZERODigitizer::Digitize(Option_t* /*option*/) 
{   
  // Creates digits from hits
  fNdigits = 0;  

  if (fVZERO && !fDigInput) {
    AliLoader *loader = fVZERO->GetLoader();
    if (!loader) {
      AliError("Can not get VZERO Loader via AliVZERO object!");
      return;
    }
    AliRunLoader* runLoader = AliRunLoader::Instance();
    for (Int_t iEvent = 0; iEvent < runLoader->GetNumberOfEvents(); ++iEvent) {
      runLoader->GetEvent(iEvent);
      if (fTask == kHits2Digits) {
	DigitizeHits();
	DigitizeSDigits();
	WriteDigits(loader);
      }
      else {
	DigitizeHits();
	WriteSDigits(loader);
      }
    }
  }
  else if (fDigInput) {
      ReadSDigits();
      DigitizeSDigits();
      AliRunLoader *currentLoader = AliRunLoader::GetRunLoader(fDigInput->GetOutputFolderName());
      AliLoader *loader = currentLoader->GetLoader("VZEROLoader");
      if (!loader) { 
	AliError("Cannot get VZERO Loader via RunDigitizer!");
	return;
      }
      WriteDigits(loader);
  }
  else {
    AliFatal("Invalid digitization task! Exiting!");
  }
}

//____________________________________________________________________________
void AliVZERODigitizer::AddDigit(Int_t pmnumber, Float_t time, Float_t width, Bool_t integrator, Short_t *chargeADC, Int_t *labels) 
 { 
 
// Adds Digit 
 
  TClonesArray &ldigits = *fDigits;  
	 
  new(ldigits[fNdigits++]) AliVZEROdigit(pmnumber,time,width,integrator,chargeADC,labels);
	 
}
//____________________________________________________________________________
void AliVZERODigitizer::AddSDigit(Int_t pmnumber, Int_t nbins, Float_t *charges, Int_t *labels) 
 { 
 
// Adds SDigit 
 
  TClonesArray &ldigits = *fDigits;  
	 
  new(ldigits[fNdigits++]) AliVZEROSDigit(pmnumber,nbins,charges,labels);
	 
}
//____________________________________________________________________________
void AliVZERODigitizer::ResetDigits()
{

// Clears Digits

  fNdigits = 0;
  if (fDigits) fDigits->Clear();
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

void AliVZERODigitizer::DigitizeHits()
{
  // Digitize the hits to the level of
  // SDigits (fTime arrays)

  for(Int_t i = 0 ; i < 64; ++i) {
    memset(fTime[i],0,fNBins[i]*sizeof(Float_t));
    fLabels[i][0] = fLabels[i][1] = fLabels[i][2] = -1;
  }
  Float_t integral = fPMResponse->Integral(-kPMRespTime,2.*kPMRespTime);
  Float_t meansPhE = fSinglePhESpectrum->Mean(0,20);

     AliLoader* loader = fVZERO->GetLoader();
     if (!loader) {
       AliError("Can not get VZERO Loader!");
       return;
     }
     loader->LoadHits();
     TTree* treeH = loader->TreeH();
     if (!treeH) {
       AliError("Cannot get TreeH!");
       return;
     }
     TClonesArray* hits = fVZERO->Hits();

//  Now makes Digits from hits
     Int_t nTracks = (Int_t) treeH->GetEntries();
     for (Int_t iTrack = 0; iTrack < nTracks; iTrack++) {
         fVZERO->ResetHits();
         treeH->GetEvent(iTrack);
         Int_t nHits = hits->GetEntriesFast();
         for (Int_t iHit = 0; iHit < nHits; iHit++) {
	   AliVZEROhit* hit = (AliVZEROhit *)hits->UncheckedAt(iHit);
	   Int_t nPhot = hit->Nphot();
	   Int_t cell  = hit->Cell();                          
	   Int_t pmt = Cell2Pmt(cell);
	   if (pmt < 0) continue;
	   Int_t trackLabel = hit->GetTrack();
	   for(Int_t l = 0; l < 3; ++l) {
	     if (fLabels[pmt][l] < 0) {
	       fLabels[pmt][l] = trackLabel;
	       break;
	     }
	   }
	   Float_t dt_scintillator = gRandom->Gaus(0,kIntTimeRes);
	   Float_t t = dt_scintillator + 1e9*hit->Tof();
	   if (pmt < 32) t += kV0CDelayCables;
	   t += fHptdcOffset[pmt];
	   Int_t nPhE;
	   Float_t prob = fCalibData->GetLightYields(pmt)*fPhotoCathodeEfficiency; // Optical losses included!
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
}


void AliVZERODigitizer::DigitizeSDigits()
{
  // Digitize the fTime arrays (SDigits) to the level of
  // Digits (fAdc arrays)
  for(Int_t i = 0 ; i < 64; ++i) {
    for(Int_t j = 0; j < kNClocks; ++j) fAdc[i][j] = 0;
    fLeadingTime[i] = fTimeWidth[i] = 0;
  }

  Float_t maximum = 0.9*fSignalShape->GetMaximum(0,200); // Not exact, one needs to do this on the convoluted
  Float_t integral2 = fSignalShape->Integral(0,200); // function. Anyway the effect is small <10% on the 2.5 ADC thr
  for (Int_t ipmt = 0; ipmt < 64; ++ipmt) {
    Float_t thr = fCalibData->GetCalibDiscriThr(ipmt,kFALSE)*kChargePerADC*maximum*fBinSize[ipmt]/integral2;
    Bool_t ltFound = kFALSE, ttFound = kFALSE;
    for (Int_t iBin = 0; iBin < fNBins[ipmt]; ++iBin) {
      Float_t t = fBinSize[ipmt]*Float_t(iBin);
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
      Float_t tadc = t - fClockOffset[ipmt];
      Int_t clock = kNClocks/2 - Int_t(tadc/25.0);
      if (clock >= 0 && clock < kNClocks)
	fAdc[ipmt][clock] += fTime[ipmt][iBin]/kChargePerADC;
    }
    AliDebug(1,Form("Channel %d Offset %f Time %f",ipmt,fClockOffset[ipmt],fLeadingTime[ipmt]));
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

  fEvenOrOdd = gRandom->Integer(2);
  for (Int_t j=0; j<64; ++j){
    for (Int_t iClock = 0; iClock < kNClocks; ++iClock) {
      Int_t integrator = (iClock + fEvenOrOdd) % 2;
      AliDebug(1,Form("ADC %d %d %f",j,iClock,fAdc[j][iClock]));
      fAdc[j][iClock]  += gRandom->Gaus(fAdcPedestal[j][integrator], fAdcSigma[j][integrator]);
    }
  }
	
}

void AliVZERODigitizer::WriteDigits(AliLoader *loader)
{
  // Take fAdc arrays filled by the previous
  // method and produce and add digits to the digit Tree

  loader->LoadDigits("UPDATE");

  if (!loader->TreeD()) loader->MakeTree("D");
  loader->MakeDigitsContainer();
  TTree* treeD  = loader->TreeD();
  DigitsArray();
  treeD->Branch("VZERODigit", &fDigits); 
  
  Short_t *chargeADC = new Short_t[kNClocks];
  for (Int_t i=0; i<64; i++) {      
    for (Int_t j = 0; j < kNClocks; ++j) {
      Int_t tempadc = Int_t(fAdc[i][j]);
      if (tempadc > 1023) tempadc = 1023;
      chargeADC[j] = tempadc;
    }
    AddDigit(i, fLeadingTime[i], fTimeWidth[i], Bool_t((10+fEvenOrOdd)%2), chargeADC, fLabels[i]);
  }
  delete [] chargeADC;

  treeD->Fill();
  loader->WriteDigits("OVERWRITE");  
  loader->UnloadDigits();     
  ResetDigits();
}

void AliVZERODigitizer::WriteSDigits(AliLoader *loader)
{
  // Take fTime arrays filled by the previous
  // method and produce and add sdigits to the sdigit Tree

  loader->LoadSDigits("UPDATE");

  if (!loader->TreeS()) loader->MakeTree("S");
  loader->MakeSDigitsContainer();
  TTree* treeS  = loader->TreeS();
  SDigitsArray();
  treeS->Branch("VZEROSDigit", &fDigits); 
  
  for (Int_t ipmt = 0; ipmt < 64; ++ipmt) {
    AddSDigit(ipmt,fNBins[ipmt],fTime[ipmt],fLabels[ipmt]);
  }

  treeS->Fill();
  loader->WriteSDigits("OVERWRITE");  
  loader->UnloadSDigits();     
  ResetDigits();
}

void AliVZERODigitizer::ReadSDigits()
{
  // Read SDigits which are then to precessed
  // in the following method
  for(Int_t i = 0 ; i < 64; ++i) {
    memset(fTime[i],0,fNBins[i]*sizeof(Float_t));
    fLabels[i][0] = fLabels[i][1] = fLabels[i][2] = -1;
  }

  // Loop over input files
  Int_t nFiles= fDigInput->GetNinputs();
  for (Int_t inputFile = 0; inputFile < nFiles; inputFile++) {
    // Get the current loader 
    AliRunLoader* currentLoader = 
      AliRunLoader::GetRunLoader(fDigInput->GetInputFolderName(inputFile));

    AliLoader *loader = currentLoader->GetLoader("VZEROLoader");
    loader->LoadSDigits("READ");
  
    // Get the tree of summable digits
    TTree* sdigitsTree = loader->TreeS();
    if (!sdigitsTree)  {
      AliError("No sdigit tree from digInput");
      continue;
    }

    // Get the branch 
    TBranch* sdigitsBranch = sdigitsTree->GetBranch("VZEROSDigit");
    if (!sdigitsBranch) {
      AliError("Failed to get sdigit branch");
      return;
    }

    // Set the branch address
    TClonesArray *sdigitsArray = NULL;
    sdigitsBranch->SetAddress(&sdigitsArray);

    // Sum contributions from the sdigits
    // Get number of entries in the tree 
    Int_t nentries  = Int_t(sdigitsBranch->GetEntries());
    for (Int_t entry = 0; entry < nentries; ++entry)  {
      sdigitsBranch->GetEntry(entry);
      // Get the number of sdigits 
      Int_t nsdigits = sdigitsArray->GetEntries();
      for (Int_t sdigit = 0; sdigit < nsdigits; sdigit++) {
	AliVZEROSDigit* sDigit = static_cast<AliVZEROSDigit*>(sdigitsArray->UncheckedAt(sdigit));
	Int_t pmNumber = sDigit->PMNumber();
	Int_t nbins = sDigit->GetNBins();
	if (nbins != fNBins[pmNumber]) {
	  AliError(Form("Incompatible number of bins between digitizer (%d) and sdigit (%d) for PM %d! Skipping sdigit!",
			fNBins[pmNumber],nbins,pmNumber));
	  continue;
	}
	// Sum the charges
	Float_t *charges = sDigit->GetCharges();
	for(Int_t iBin = 0; iBin < nbins; ++iBin) fTime[pmNumber][iBin] += charges[iBin];
	// and the labels
	Int_t *labels = sDigit->GetTracks();
	Int_t j = 0;
	for(Int_t i = 0; i < 3; ++i) {
	  if (fLabels[pmNumber][i] < 0) {
	    if (labels[j] < 0) break;
	    fLabels[pmNumber][i] = labels[j];
	    j++;
	  }
	}
      }
    }
    loader->UnloadSDigits();
  }
}

//____________________________________________________________________
TClonesArray*
AliVZERODigitizer::DigitsArray() 
{
  // Initialize digit array if not already and
  // return pointer to it. 
  if (!fDigits) { 
    fDigits = new TClonesArray("AliVZEROdigit", 64);
    fNdigits = 0;
  }
  return fDigits;
}

//____________________________________________________________________
TClonesArray*
AliVZERODigitizer::SDigitsArray() 
{
  // Initialize sdigit array if not already and
  // return pointer to it. 
  if (!fDigits) { 
    fDigits = new TClonesArray("AliVZEROSDigit", 64);
    fNdigits = 0;
  }
  return fDigits;
}
