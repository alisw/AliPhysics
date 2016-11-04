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

/* $Id: AliADDigitizer.cxx  $ */

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
#include <TH2F.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TTree.h>
#include <TSpline.h>

// --- AliRoot header files ---
#include "AliRun.h"
#include "AliDetector.h"
#include "AliAD.h"
#include "AliADhit.h"
#include "AliADConst.h"
#include "AliADRecoParam.h"
#include "AliRunLoader.h"
#include "AliLoader.h"
#include "AliGRPObject.h"
#include "AliDigitizationInput.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliCDBEntry.h"
#include "AliADCalibData.h"
#include "AliCTPTimeParams.h"
#include "AliLHCClockPhase.h"
#include "AliADdigit.h"
#include "AliADDigitizer.h"
#include "AliADSDigit.h"
#include "AliADTriggerSimulator.h"
#include "AliLog.h"
#include "AliHeader.h"
#include "AliGenEventHeader.h"

ClassImp(AliADDigitizer);

//____________________________________________________________________________
AliADDigitizer::AliADDigitizer()
  : AliDigitizer()
  , fCalibData(GetCalibData())
  , fNdigits(0)
  , fDigits(0)
  , fChargeSignalShape(NULL)
  , fTimeSignalShape(NULL)
  , fThresholdShape(NULL)
  , fTailBegin(16)
  , fTailEnd(20)
  , fTS(NULL)
  , fEvenOrOdd(kFALSE)
  , fTask(kHits2Digits)
  , fAD(NULL)
{
  // default constructor
  // Initialize OCDB and containers used in the digitization
  Init();
}

//____________________________________________________________________________
AliADDigitizer::AliADDigitizer(AliAD *AD, DigiTask_t task)
  : AliDigitizer()
  , fCalibData(GetCalibData())
  , fNdigits(0)
  , fDigits(0)
  , fChargeSignalShape(NULL)
  , fTimeSignalShape(NULL)
  , fThresholdShape(NULL)
  , fTailBegin(16)
  , fTailEnd(20)
  , fTS(NULL)
  , fEvenOrOdd(kFALSE)
  , fTask(task)
  , fAD(AD)
{
  // constructor
  // Initialize OCDB and containers used in the digitization
  Init();
}

//____________________________________________________________________________
AliADDigitizer::AliADDigitizer(AliDigitizationInput* digInput)
  : AliDigitizer(digInput)
  , fCalibData(GetCalibData())
  , fNdigits(0)
  , fDigits(0)
  , fChargeSignalShape(NULL)
  , fTimeSignalShape(NULL)
  , fThresholdShape(NULL)
  , fTailBegin(16)
  , fTailEnd(20)
  , fTS(NULL)
  , fEvenOrOdd(kFALSE)
  , fTask(kHits2Digits)
  , fAD(NULL)
{
  // constructor
  // Initialize OCDB and containers used in the digitization
  Init();
}

//____________________________________________________________________________
AliADDigitizer::~AliADDigitizer()
{
  // destructor
  if (fDigits) {
    fDigits->Delete();
    delete fDigits;
    fDigits = NULL;
  }

  if (fTimeSignalShape) {
    delete fTimeSignalShape;
    fTimeSignalShape = NULL;
  }
  if (fChargeSignalShape) {
    delete fChargeSignalShape;
    fChargeSignalShape = NULL;
  }

  for(Int_t ch=0; ch<16; ++ch) {
    if (fTime[ch]) {
      delete [] fTime[ch];
    }
    fTime[ch] = NULL;

    for (Int_t centerBCint=0; centerBCint<2; ++centerBCint) {
      if (fTailVsTotalCharge[ch][centerBCint])
	delete fTailVsTotalCharge[ch][centerBCint];
      fTailVsTotalCharge[ch][centerBCint] = NULL;
    }
  }
}

Bool_t AliADDigitizer::SetupTailVsTotalCharge()
{
  const AliADRecoParam p; // for now tail begin does not depend on the type of recoparam obj
  fTailBegin = p.GetTailBegin();
  fTailEnd   = p.GetTailEnd();

  for (Int_t ch=0; ch<16; ++ch) {
    fTailVsTotalCharge[ch][0] = NULL;
    fTailVsTotalCharge[ch][1] = NULL;
  }
  AliCDBManager *man = AliCDBManager::Instance();
  if (!man)
    AliFatal("!man");

  AliCDBEntry *entry = AliCDBManager::Instance()->Get("AD/Calib/Saturation");
  if (!entry)
    AliFatal("!entry");

  fTS = dynamic_cast<TTree*>(entry->GetObject());
  if (!fTS)
    AliFatal("!fTS");

  Int_t chOffline,chOnline;
  TClonesArray *f_Int[2] = { NULL, NULL };
  Bool_t  doExtrapolation[kADNClocks];
  Float_t extrapolationThresholds[kADNClocks];
  fTS->SetBranchAddress("chOffline", &chOffline);
  fTS->SetBranchAddress("chOnline",  &chOnline);
  fTS->SetBranchAddress("f_Int0",    &f_Int[0]);
  fTS->SetBranchAddress("f_Int1",    &f_Int[1]);
  fTS->SetBranchAddress("doExtrapolation",         &doExtrapolation);
  fTS->SetBranchAddress("extrapolationThresholds", &extrapolationThresholds);

  TF1 *f0=NULL, *f1=NULL;
  for (Int_t ch=0; ch<16; ++ch) {
    f_Int[0]->Clear();
    f_Int[1]->Clear();
    fTS->GetEntry(ch);

    fTailVsTotalCharge[ch][0] = new TGraph;
    fTailVsTotalCharge[ch][1] = new TGraph;

    for (Int_t j=0; j<1000; ++j) {
      const Float_t tail = (j<20 ? 0.5*j : Float_t(j-10));
      Float_t charge0 = tail;
      Float_t charge1 = tail;
      for (Int_t bc=0; bc<fTailBegin; ++bc) {
	if (!doExtrapolation[bc])
	  continue;
	const Bool_t integrator = ((bc%2) == 0); // bc=10 -> integrator=kTRUE
	f0 = dynamic_cast<TF1*>(f_Int[!integrator]->At(bc));
	f1 = dynamic_cast<TF1*>(f_Int[ integrator]->At(bc));

	if (!f0 || !f1) {
	  AliWarning("f0==NULL || f1==NULL");
	  continue;
	}

	charge0 += TMath::Max(0.0, f0->Eval(tail));
	charge1 += TMath::Max(0.0, f1->Eval(tail));
      }
      fTailVsTotalCharge[ch][0]->SetPoint(fTailVsTotalCharge[ch][0]->GetN(), charge0, tail);
      fTailVsTotalCharge[ch][1]->SetPoint(fTailVsTotalCharge[ch][1]->GetN(), charge1, tail);
    }
  }
  fTS->ResetBranchAddresses();
  return kTRUE;
}

//____________________________________________________________________________
Bool_t AliADDigitizer::Init()
{
  // Initialises the digitizer
  // Initialize OCDB and containers used in the digitization

  // check if the digitizer was already initialized
  if (fTimeSignalShape)
    return kTRUE;

  for (Int_t ch=0; ch<16; ++ch) {
    fTailVsTotalCharge[ch][0] = NULL;
    fTailVsTotalCharge[ch][1] = NULL;
  }

  if (!SetupTailVsTotalCharge())
    return kFALSE;

  fTimeSignalShape   = new TF1("ADTimeSignalShape",  this, &AliADDigitizer::TimeSignalShape,  0,200,6, "AliADDigitizer", "TimeSignalShape");
  fChargeSignalShape = new TF1("ADChargeSignalShape",this, &AliADDigitizer::ChargeSignalShape,0,300,3, "AliADDigitizer", "ChargeSignalShape");
  fThresholdShape    = new TF1("ADThresholdShape",   this, &AliADDigitizer::ThresholdShape,   0, 50,1, "AliADDigitizer", "ThresholdShape");

  fTimeSignalShape->SetParameters(-1.07335e+00,2.16002e+01,-1.26133e-01,
                                  1.41619e+00,5.50334e-01,3.86111e-01);

  AliCDBEntry *entry = NULL;
  // Now get the CTP L0->L1 delay
  entry = AliCDBManager::Instance()->Get("GRP/CTP/CTPtiming");
  if (!entry)
    AliFatal("CTP timing parameters are not found in OCDB !");
  AliCTPTimeParams *ctpParams = (AliCTPTimeParams*)entry->GetObject();
  Float_t l1Delay = (Float_t)ctpParams->GetDelayL1L0()*25.0;

  entry = AliCDBManager::Instance()->Get("GRP/CTP/TimeAlign");
  if (!entry)
    AliFatal("CTP time-alignment is not found in OCDB !");
  AliCTPTimeParams *ctpTimeAlign = (AliCTPTimeParams*)entry->GetObject();
  l1Delay += ((Float_t)ctpTimeAlign->GetDelayL1L0()*25.0);

  entry = AliCDBManager::Instance()->Get("AD/Calib/TimeDelays");
  if (!entry)
    AliFatal("AD time delays are not found in OCDB !");
  TH1F *TimeDelays = (TH1F*)entry->GetObject();

  entry = AliCDBManager::Instance()->Get("GRP/Calib/LHCClockPhase");
  if (!entry)
    AliFatal("LHC clock-phase shift is not found in OCDB !");
  AliLHCClockPhase *phase = (AliLHCClockPhase*)entry->GetObject();

  // Get Pulse shape parameters
  entry = AliCDBManager::Instance()->Get("AD/Calib/PulseShapes");
  if (!entry)
    AliFatal("AD pulse shapes are not found in OCDB !");
  TH2F *PulseShapes = (TH2F*)entry->GetObject();

  // Time slewing splines
  GetTimeSlewingSplines();

  // Try to extrapolate the splines
  ExtrapolateSplines();

  for(Int_t i=0; i<16; ++i) {
    fCssOffset[i] = PulseShapes->GetBinContent(i+1,1);
    fCssTau[i]    = PulseShapes->GetBinContent(i+1,2);
    fCssSigma[i]  = PulseShapes->GetBinContent(i+1,3);

    for(Int_t j=0; j<kADNClocks; ++j) fAdc[i][j] = 0;
    fLeadingTime[i] = fTimeWidth[i] = 0;

    fPmGain[i]    = fCalibData->GetGain(i);

    fAdcPedestal[i][0] = fCalibData->GetPedestal(i);
    fAdcSigma[i][0]    = fCalibData->GetSigma(i);
    fAdcPedestal[i][1] = fCalibData->GetPedestal(i+16);
    fAdcSigma[i][1]    = fCalibData->GetSigma(i+16);

    const Int_t board = AliADCalibData::GetBoardNumber(i);
    fNBins[i]   = TMath::Nint(((Float_t)(fCalibData->GetMatchWindow(board)+1)*25.0+
			       (Float_t)kADMaxTDCWidth*fCalibData->GetWidthResolution(board))/
			      fCalibData->GetTimeResolution(board));
    fNBinsLT[i] = TMath::Nint(((Float_t)(fCalibData->GetMatchWindow(board)+1)*25.0)/
                              fCalibData->GetTimeResolution(board));
    fBinSize[i] = fCalibData->GetTimeResolution(board);

    fHptdcOffset[i] = (((Float_t)fCalibData->GetRollOver(board)-
                        (Float_t)fCalibData->GetTriggerCountOffset(board))*25.0
                       +fCalibData->GetTimeOffset(i)
                       -l1Delay
                       -phase->GetMeanPhase()
                       -TimeDelays->GetBinContent(i+1)
                       -kADOffset);

    fClockOffset[i] = (((Float_t)fCalibData->GetRollOver(board)-
                        (Float_t)fCalibData->GetTriggerCountOffset(board))*25.0
                       +fCalibData->GetTimeOffset(i)
                       -l1Delay
                       -kADOffset);

    fTime[i] = new Float_t[fNBins[i]];
    memset(fTime[i],0,fNBins[i]*sizeof(Float_t));

    //std::cout<<"AD: "<<" fNBins = "<<fNBins[i]<<" fNBinsLT = "<<fNBinsLT[i]<<" fHptdcOffset = "<<fHptdcOffset[i]<<" fClockOffset = "<<fClockOffset[i]<<std::endl;
  }

  return kTRUE;
}

//____________________________________________________________________________
void AliADDigitizer::Digitize(Option_t* /*option*/)
{
  // Creates digits from hits
  fNdigits = 0;

  if (fAD && !fDigInput) {
    AliLoader *loader = fAD->GetLoader();
    if (!loader) {
      AliError("Can not get AD Loader via AliAD object!");
      return;
    }
    AliRunLoader* runLoader = AliRunLoader::Instance();

    for (Int_t iEvent=0; iEvent<runLoader->GetNumberOfEvents(); ++iEvent) {
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
    AliLoader *loader = currentLoader->GetLoader("ADLoader");
    if (!loader) {
      AliError("Cannot get AD Loader via RunDigitizer!");
      return;
    }
    WriteDigits(loader);
  }
  else {
    AliFatal("Invalid digitization task! Exiting!");
  }
}

//____________________________________________________________________________
void AliADDigitizer::DigitizeHits()
{
  // Digitize the hits to the level of
  // SDigits (fTime arrays)
  Int_t nTotPhot[16];
  Float_t PMTime[16];
  Float_t PMTimeWeight[16];
  Int_t nPMHits[16];

  for(Int_t i=0; i<16; ++i) {
    memset(fTime[i],0,fNBins[i]*sizeof(Float_t));
    fLabels[i][0] = fLabels[i][1] = fLabels[i][2] = -1;
    nTotPhot[i] = 0;
    PMTime[i] = 10000;
    PMTimeWeight[i] = 0;
    nPMHits[i] = 0;
  }

  AliLoader* loader = fAD->GetLoader();
  if (!loader) {
    AliError("Can not get AD Loader!");
    return;
  }
  loader->LoadHits();
  TTree* treeH = loader->TreeH();
  if (!treeH) {
    AliError("Cannot get TreeH!");
    return;
  }
  TClonesArray* hits = fAD->Hits();

  //Loop over hits
  Int_t nTracks = (Int_t) treeH->GetEntries();
  for(Int_t iTrack=0; iTrack<nTracks; iTrack++) {
    fAD->ResetHits();
    treeH->GetEvent(iTrack);
    Int_t nHits = hits->GetEntriesFast();
    for (Int_t iHit=0; iHit<nHits; iHit++) {
      AliADhit* hit = (AliADhit *)hits->UncheckedAt(iHit);
      Int_t nPhot = hit->GetNphot();
      Int_t pmt  = hit->GetCell();
      if (pmt < 0) continue;
      Int_t trackLabel = hit->GetTrack();
      for(Int_t l=0; l<3; ++l) {
        if (fLabels[pmt][l] < 0) {
          fLabels[pmt][l] = trackLabel;
          break;
        }
      }
      Float_t dt_scintillator = gRandom->Gaus(0,kADIntTimeRes);
      Float_t t = dt_scintillator + hit->GetTof();
      nTotPhot[pmt] += nPhot;
      nPMHits[pmt]++;
      //PMTime[pmt] += t*nPhot*nPhot;
      //PMTimeWeight[pmt] += nPhot*nPhot;
      if (PMTime[pmt]>t) PMTime[pmt] = t;

    }//hit loop
  }//track loop

  //Now makes SDigits from hits
  for(Int_t iPM=0; iPM<16; iPM++) {
    if (nPMHits[iPM]==0 || nTotPhot[iPM]==0) {
      PMTime[iPM] = 0.0;
      continue;
    }
    //PMTime[iPM] = PMTime[iPM]/PMTimeWeight[iPM];
    PMTime[iPM] += fHptdcOffset[iPM];

    fChargeSignalShape->SetParameters(fCssOffset[iPM],fCssTau[iPM],fCssSigma[iPM]);
    Float_t integral = fChargeSignalShape->Integral(0,300);
    //std::cout<<"Integral = "<<integral<<std::endl;

    Float_t charge = nTotPhot[iPM]*fPmGain[iPM]*fBinSize[iPM]/integral;

    Int_t firstBin = TMath::Max(0,(Int_t)((PMTime[iPM])/fBinSize[iPM]));
    Int_t lastBin = fNBins[iPM]-1;
    //std::cout<<"First Bin: "<<firstBin<<std::endl;
    for(Int_t iBin = firstBin; iBin <= lastBin; ++iBin) {
      Float_t tempT = fBinSize[iPM]*(0.5+iBin)-PMTime[iPM];
      if (tempT <= 0) continue;
      fTime[iPM][iBin] += charge*fChargeSignalShape->Eval(tempT);
    }
  }//PM loop
  loader->UnloadHits();
}

//____________________________________________________________________________
void AliADDigitizer::DigitizeSDigits()
{
  // Digitize the fTime arrays (SDigits) to the level of
  // Digits (fAdc arrays)
  Float_t fMCTime[16];
  for(Int_t i=0; i<16; ++i) {
    for(Int_t j=0; j<kADNClocks; ++j) fAdc[i][j] = 0;
    fMCTime[i] = fLeadingTime[i] = fTimeWidth[i] = 0;
  }

  for (Int_t ipmt=0; ipmt<16; ++ipmt) {
    fChargeSignalShape->SetParameters(fCssOffset[ipmt],fCssTau[ipmt],fCssSigma[ipmt]);
    Float_t maximum = 0.9*fChargeSignalShape->GetMaximum(0,300);
    Float_t integral = fChargeSignalShape->Integral(0,300);
    Float_t thr = fCalibData->GetCalibDiscriThr(ipmt)*kADChargePerADC*maximum*fBinSize[ipmt]/integral;
    //Float_t thr = 0;

    Bool_t ltFound = kFALSE, ttFound = kFALSE;
    for (Int_t iBin = 1; iBin<fNBins[ipmt]; ++iBin) {
      Float_t t = fBinSize[ipmt]*Float_t(iBin);
      if (fTime[ipmt][iBin] > 0.0) {
        if (!ltFound && (iBin < fNBinsLT[ipmt])) {
          ltFound = kTRUE;
          fMCTime[ipmt] = t;
          //std::cout<<"Leading Bin: "<<iBin<<std::endl;
          //std::cout<<"Leading TADC: "<<t-fClockOffset[ipmt]<<std::endl;
        }
      }
      if (fTime[ipmt][iBin-1] > thr && fTime[ipmt][iBin] < thr) {
        if (ltFound) {
          if (!ttFound) {
            ttFound = kTRUE;
            fTimeWidth[ipmt] = t - fMCTime[ipmt];
          }
        }
      }
      const Float_t tadc = t - fClockOffset[ipmt];
      const Int_t  clock = kADNClocks/2 + Int_t(tadc/25.0);
      if (clock >= 0 && clock < kADNClocks)
        fAdc[ipmt][clock] += fTime[ipmt][iBin]/kADChargePerADC;
    }
    AliDebug(1,Form("Channel %d Offset %f Time %f",ipmt,fClockOffset[ipmt],fMCTime[ipmt]));
    const Int_t board = AliADCalibData::GetBoardNumber(ipmt);
    if (ltFound && ttFound) {
      fTimeWidth[ipmt] = fCalibData->GetWidthResolution(board)*
        Float_t(Int_t(fTimeWidth[ipmt]/fCalibData->GetWidthResolution(board)));
      if (fTimeWidth[ipmt] < Float_t(kADMinTDCWidth)*fCalibData->GetWidthResolution(board))
        fTimeWidth[ipmt] = Float_t(kADMinTDCWidth)*fCalibData->GetWidthResolution(board);
      if (fTimeWidth[ipmt] > Float_t(kADMaxTDCWidth)*fCalibData->GetWidthResolution(board))
        fTimeWidth[ipmt] = Float_t(kADMaxTDCWidth)*fCalibData->GetWidthResolution(board);
    }
  }

  fEvenOrOdd = gRandom->Integer(2);

  // CM: pulse shape from saturation correction
  AdjustPulseShapeADC();

  for (Int_t j=0; j<16; ++j){
    Float_t adcSignal = 0.0;
    Float_t adcClock = 0.0;
    for (Int_t iClock=0; iClock<kADNClocks; ++iClock) {
      Int_t integrator = (iClock + fEvenOrOdd) % 2;
      AliDebug(1,Form("ADC %d %d %f",j,iClock,fAdc[j][iClock]));
      fAdc[j][iClock]  += gRandom->Gaus(fAdcPedestal[j][integrator], fAdcSigma[j][integrator]);
    }
    for (Int_t iClock=0; iClock<kADNClocks; ++iClock) {
      Int_t integrator = (iClock + fEvenOrOdd) % 2;
      adcClock = (Int_t)fAdc[j][iClock];
      if (fAdc[j][iClock]>1023) adcClock = 1023;
      adcClock -= fAdcPedestal[j][integrator];
      if (adcClock< 4*fAdcSigma[j][integrator]) adcClock = 0;
      adcSignal += adcClock;
    }
    fThresholdShape->SetParameter(0,fCalibData->GetCalibDiscriThr(j));
    if (gRandom->Rndm() > fThresholdShape->Eval(adcSignal)) fMCTime[j] = -1024.0;
    if (fThresholdShape->Eval(adcSignal)<1e-2) fMCTime[j] = -1024.0;
    fLeadingTime[j] = UnCorrectLeadingTime(j,fMCTime[j],adcSignal);

  }
  //Fill BB and BG flags in trigger simulator
  AliADTriggerSimulator * triggerSimulator = new AliADTriggerSimulator();
  triggerSimulator->FillFlags(fBBFlag,fBGFlag,fLeadingTime);
}

void AliADDigitizer::AdjustPulseShapeADC()
{
  TClonesArray *f_Int[2] = { NULL, NULL };
  Float_t extrapolationThresholds[kADNClocks];
  Bool_t  doExtrapolation[kADNClocks];
  fTS->SetBranchAddress("f_Int0", &f_Int[0]);
  fTS->SetBranchAddress("f_Int1", &f_Int[1]);
  fTS->SetBranchAddress("doExtrapolation",         &doExtrapolation);
  fTS->SetBranchAddress("extrapolationThresholds", &extrapolationThresholds);

  // fAdc contains not yet quantized (Float_t) and not yet clipped (at 1024) ADC values without pedestal shift and noise
  for (Int_t ch=0; ch<16; ++ch) {
    f_Int[0]->Clear();
    f_Int[1]->Clear();
    fTS->GetEntry(ch);
    Float_t totalCharge=0.0f;
    for (Int_t bc=0; bc<kADNClocks; ++bc) {
      totalCharge += fAdc[ch][bc];
    }
    const Float_t tail = fTailVsTotalCharge[ch][fEvenOrOdd]->Eval(totalCharge);
    if (tail < 3.0f || !doExtrapolation[kADNClocks/2])
      continue; // do not modify the pulse shape for small pulses or with default saturation object

    Float_t newADC[kADNClocks];
    Float_t chargeLastBC = 0.0; // charge in last extrapolated BC
    for (Int_t bc=0; bc<kADNClocks; ++bc) {
      newADC[bc] = 0.0;
      if (!doExtrapolation[bc])
	continue;
      const Bool_t integrator = ((bc+fEvenOrOdd) % 2);
      TF1 *f = dynamic_cast<TF1*>(f_Int[integrator]->At(bc));
      chargeLastBC = newADC[bc] = f->Eval(tail);
    }
    // we use a linear function for the tail starting at 80% of the last BC before the tail
    // [ (n+1) n(n+1)/2; 1 n] [a;b] = [tail; chargeLastBC]
    const Float_t n = fTailEnd-fTailBegin;
    if (n < 1)
      AliFatalF("n=%f fTailBegin=%d fTailEnd=%d", n, fTailBegin, fTailEnd);
    const Float_t det = n*(n+1)/2;

    // allowed range for chargeLastBC for a>=0 and b>=0
    const Float_t rangeChargeLastBC[2] = { tail/(n+1), 2*tail/(n+1) };
    const Float_t eps=1e-6;
    chargeLastBC = TMath::Min(rangeChargeLastBC[1]-eps, 0.8f*chargeLastBC);
    chargeLastBC = TMath::Max(rangeChargeLastBC[0]+eps,      chargeLastBC);
    const Float_t a = n*(-(n+1)/2*chargeLastBC + tail)/det;
    const Float_t b =   ( (n+1)  *chargeLastBC - tail)/det;
    if (a < 0 || b < 0) {
      AliErrorF("charge redistribution failed for Ch%02d: n=%f a=%e b=%e chargeLastBC=%f tail=%f",
		ch, n, a, b, chargeLastBC, tail);
      continue;
    }
    // fill the tail BCs with linear function
    for (Int_t bc=fTailBegin; bc<=fTailEnd; ++bc) {
      newADC[bc] = a + b*(fTailEnd-bc);
    }

    AliDebugF(5, "OLD: Ch%02d %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f", ch,
	      fAdc[ch][20], fAdc[ch][19], fAdc[ch][18], fAdc[ch][17], fAdc[ch][16], fAdc[ch][15]);
    AliDebugF(5, "NEW: Ch%02d %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f", ch,
	      newADC[20], newADC[19], newADC[18], newADC[17], newADC[16], newADC[15]);
    Float_t newCharge=0.0;
    for (Int_t bc=0; bc<kADNClocks; ++bc) {
      fAdc[ch][bc] = newADC[bc];
      newCharge   += newADC[bc];
    }
    AliDebugF(5, "Ch%02d: totalCharge,newCharge= %f %f  (tail=%f)", ch, totalCharge, newCharge, tail);
    if (TMath::Abs(totalCharge - newCharge) > 2.0f)
      AliWarningF("Ch%02d: difference between totalCharge=%f and newCharge=%f is too large (tail=%f)", ch, totalCharge, newCharge, tail);
  }
  fTS->ResetBranchAddresses();
}

//____________________________________________________________________________
void AliADDigitizer::ReadSDigits()
{
  // Read SDigits which are then to precessed
  // in the following method
  for(Int_t i=0; i<16; ++i) {
    memset(fTime[i],0,fNBins[i]*sizeof(Float_t));
    fLabels[i][0] = fLabels[i][1] = fLabels[i][2] = -1;
  }

  // Loop over input files
  Int_t nFiles= fDigInput->GetNinputs();
  for (Int_t inputFile=0; inputFile<nFiles; inputFile++) {
    // Get the current loader
    AliRunLoader* currentLoader =
      AliRunLoader::GetRunLoader(fDigInput->GetInputFolderName(inputFile));

    AliLoader *loader = currentLoader->GetLoader("ADLoader");
    loader->LoadSDigits("READ");

    // Get the tree of summable digits
    TTree* sdigitsTree = loader->TreeS();
    if (!sdigitsTree)  {
      AliError("No sdigit tree from digInput");
      continue;
    }

    // Get the branch
    TBranch* sdigitsBranch = sdigitsTree->GetBranch("ADSDigit");
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
    for (Int_t entry=0; entry<nentries; ++entry)  {
      sdigitsBranch->GetEntry(entry);
      // Get the number of sdigits
      Int_t nsdigits = sdigitsArray->GetEntries();

      for (Int_t sdigit=0; sdigit<nsdigits; sdigit++) {
        AliADSDigit* sDigit = static_cast<AliADSDigit*>(sdigitsArray->UncheckedAt(sdigit));
        Int_t pmNumber = sDigit->PMNumber();
        Int_t nbins = sDigit->GetNBins();
        if (nbins != fNBins[pmNumber]) {
          AliError(Form("Incompatible number of bins between digitizer (%d) and sdigit (%d) for PM %d! Skipping sdigit!",
                        fNBins[pmNumber],nbins,pmNumber));
          continue;
        }
        // Sum the charges
        Float_t *charges = sDigit->GetCharges();
        for(Int_t iBin=0; iBin<nbins; ++iBin) fTime[pmNumber][iBin] += charges[iBin];
        // and the labels
        Int_t *labels = sDigit->GetTracks();
        Int_t j = 0;
        for(Int_t i=0; i<3; ++i) {
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


//____________________________________________________________________________
void AliADDigitizer::WriteDigits(AliLoader *loader)
{
  // Take fAdc arrays filled by the previous
  // method and produce and add digits to the digit Tree

  loader->LoadDigits("UPDATE");

  if (!loader->TreeD()) loader->MakeTree("D");
  loader->MakeDigitsContainer();
  TTree* treeD  = loader->TreeD();
  DigitsArray();
  treeD->Branch("ADDigit", &fDigits);

  Short_t *chargeADC = new Short_t[kADNClocks];
  for (Int_t i=0; i<16; i++) {
    for (Int_t j=0; j < kADNClocks; ++j) {
      Int_t tempadc = Int_t(fAdc[i][j]);
      if (tempadc > 1023) tempadc = 1023;
      chargeADC[j] = tempadc;
    }
    AddDigit(i, fLeadingTime[i], fTimeWidth[i], Bool_t((10+fEvenOrOdd)%2), chargeADC, fBBFlag[i], fBGFlag[i], fLabels[i]);
  }
  delete [] chargeADC;

  treeD->Fill();
  loader->WriteDigits("OVERWRITE");
  loader->UnloadDigits();
  ResetDigits("C");
}

//____________________________________________________________________________
void AliADDigitizer::WriteSDigits(AliLoader *loader)
{
  // Take fTime arrays filled by the previous
  // method and produce and add sdigits to the sdigit Tree

  loader->LoadSDigits("UPDATE");

  if (!loader->TreeS()) loader->MakeTree("S");
  loader->MakeSDigitsContainer();
  TTree* treeS  = loader->TreeS();
  SDigitsArray();
  treeS->Branch("ADSDigit", &fDigits);
  //fAD->MakeBranchInTree(treeS,"AD",&fDigits,8000,"");

  for (Int_t ipmt=0; ipmt < 16; ++ipmt) {
    AddSDigit(ipmt,fNBins[ipmt],fTime[ipmt],fLabels[ipmt]);
  }

  treeS->Fill();
  loader->WriteSDigits("OVERWRITE");
  loader->UnloadSDigits();
  ResetDigits("C");
}



//____________________________________________________________________________
void AliADDigitizer::AddDigit(Int_t pmnumber, Float_t time, Float_t width, Bool_t integrator, Short_t *chargeADC, Bool_t bbFlag, Bool_t bgFlag, Int_t *labels)
{
  // Adds Digit
  TClonesArray &ldigits = *fDigits;
  new(ldigits[fNdigits++]) AliADdigit(pmnumber,time,width,integrator,chargeADC,bbFlag,bgFlag,labels);

}
//____________________________________________________________________________
void AliADDigitizer::AddSDigit(Int_t pmnumber, Int_t nbins, Float_t *charges, Int_t *labels)
{
  // Adds SDigit
  TClonesArray &ldigits = *fDigits;
  new(ldigits[fNdigits++]) AliADSDigit(pmnumber,nbins,charges,labels);

}
//____________________________________________________________________________
void AliADDigitizer::ResetDigits(Option_t *opt)
{
  // Clears Digits
  fNdigits = 0;
  if (fDigits)
    fDigits->Clear(opt);
}

//____________________________________________________________________________
AliADCalibData* AliADDigitizer::GetCalibData() const
{
  AliCDBManager *man = AliCDBManager::Instance();
  AliCDBEntry *entry = man->Get("AD/Calib/Data");
  if (!entry) {
    AliWarning("Load of calibration data from default storage failed!");
    AliWarning("Calibration data will be loaded from local storage ($ALICE_ROOT)");

    man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
    man->SetRun(1);
    entry = man->Get("AD/Calib/Data");
  }
  // Retrieval of data in directory AD/Calib/Data:
  AliADCalibData *calibdata = NULL;

  if (entry)
    calibdata = (AliADCalibData*) entry->GetObject();
  if (!calibdata)
    AliFatal("No calibration data from calibration database !");

  //calibdata->PrintConfig();
  return calibdata;
}
//____________________________________________________________________________
Float_t AliADDigitizer::UnCorrectLeadingTime(Int_t i, Float_t time, Float_t adc) const
{
  // UnCorrect the MC time
  // for slewing effect and
  // misalignment of the channels
  const Double_t fTOF[4] = {65.2418, 65.1417, 56.6459, 56.7459};

  if (time < 1e-6) return time;
  if (adc < 1) return time;

  // Slewing and offset correction
  Int_t board = AliADCalibData::GetBoardNumber(i);
  //std::cout<<"MC time: "<<time<<std::endl;
  time -= fHptdcOffset[i];
  //std::cout<<"TOF: "<<time<<std::endl;
  time -= fTOF[i/4];
  if(adc<30 && fTimeSlewingExtpol[i]) time += fTimeSlewingExtpol[i]->Eval(TMath::Log10(1/adc))*fCalibData->GetTimeResolution(board);

  else time += fTimeSlewingSpline[i]->Eval(TMath::Log10(1/adc))*fCalibData->GetTimeResolution(board);
  //std::cout<<"Charge: "<<adc<<std::endl;
  //std::cout<<"Leading time: "<<time<<std::endl;

  const Float_t smearedTime = SmearLeadingTime(i,time);
  return smearedTime;
}
//____________________________________________________________________________
Float_t AliADDigitizer::SmearLeadingTime(Int_t i, Float_t time) const
{
  const Int_t runNumber = AliCDBManager::Instance()->GetRun();
  Float_t sigmaADA = 0, sigmaADC = 0;
  if (runNumber < 225753) {
    sigmaADA = 0.45;
    sigmaADC = 0.15;
  }
  if (runNumber > 225753 && runNumber < 226501) {
    sigmaADA = 1.0;
    sigmaADC = 0.15;
  }
  if (runNumber > 226501) {
    sigmaADA = 0.50;
    sigmaADC = 0.15;
  }

  if (i<8) time += gRandom->Gaus(1.25,sigmaADC);
  else     time += gRandom->Gaus(1.05,sigmaADA);

  return time;
}
//_____________________________________________________________________________
void AliADDigitizer::GetTimeSlewingSplines()
{
  AliCDBManager *man = AliCDBManager::Instance();
  AliCDBEntry *entry = man->Get("AD/Calib/TimeSlewing");

  TList *fListSplines = NULL;;
  if (entry) fListSplines = (TList*) entry->GetObject();
  if (!fListSplines)
    AliFatal("No time slewing correction from calibration database !");

  for (Int_t i=0; i<16; i++)
    fTimeSlewingSpline[i] = (TSpline3*)(fListSplines->At(i));
}
//_____________________________________________________________________________
void AliADDigitizer::ExtrapolateSplines()
{
  TH1F *hTimeVsSignal = NULL;
  for (Int_t i=0; i<16; i++) {
    TCanvas *c = new TCanvas("c", " ",0,0,1,1);
    c->cd();
    fTimeSlewingSpline[i]->Paint();
    hTimeVsSignal = fTimeSlewingSpline[i]->GetHistogram();

    TString TimeSlewingFitName = "hTimeSlewingFit";
    TimeSlewingFitName += i;
    fTimeSlewingExtpol[i] = new TF1(TimeSlewingFitName.Data(),"[0]+[1]*TMath::Power(10,-x*[2])",-3,0);
    fTimeSlewingExtpol[i]->SetParameter(0,650);
    fTimeSlewingExtpol[i]->SetParLimits(0,200,3000);
    fTimeSlewingExtpol[i]->SetParameter(1,450);
    fTimeSlewingExtpol[i]->SetParLimits(1,50,1000);
    fTimeSlewingExtpol[i]->SetParameter(2,-0.5);
    fTimeSlewingExtpol[i]->SetParLimits(2,-0.9,-0.05);
    fTimeSlewingExtpol[i]->SetLineColor(kMagenta);
    Int_t fitStatus =  hTimeVsSignal->Fit(TimeSlewingFitName.Data(),"R"," ",-2.5,-1.5);
    if(fitStatus != 0) {
      AliWarning(Form("Extrapolation of spline %d not succesfull",i));
      fTimeSlewingExtpol[i] = 0x0;
    }
    delete c;
  }
}
//____________________________________________________________________________
double AliADDigitizer::ChargeSignalShape(double *x, double *par)
{
  // this function simulates the charge shape
  const Double_t xx = x[0];
  return TMath::Exp(-0.5*TMath::Power(TMath::Log((xx+par[0])/par[1])/par[2],2));
}

//____________________________________________________________________________
double AliADDigitizer::ThresholdShape(double *x, double *par)
{
  // this function simulates the threshold shape
  const Double_t xx = x[0];
  return 1/(1+TMath::Exp(-xx + par[0]));
}

//____________________________________________________________________________
double AliADDigitizer::TimeSignalShape(double *x, double *par)
{
  // this function simulates the time shape

  Double_t xx = x[0];
  if (xx <= par[0]) return 0;
  Double_t a = 1./TMath::Power((xx-par[0])/par[1],1./par[2]);
  if (xx <= par[3]) return a;
  Double_t b = 1./TMath::Power((xx-par[3])/par[4],1./par[5]);
  Double_t f = a*b/(a+b);
  AliDebug(100,Form("x=%f func=%f",xx,f));
  return f;
}
//____________________________________________________________________
TClonesArray* AliADDigitizer::DigitsArray()
{
  // Initialize digit array if not already and
  // return pointer to it.
  if (!fDigits) {
    fDigits = new TClonesArray("AliADdigit", 16);
    fNdigits = 0;
  }
  return fDigits;
}

//____________________________________________________________________
TClonesArray* AliADDigitizer::SDigitsArray()
{
  // Initialize sdigit array if not already and
  // return pointer to it.
  if (!fDigits) {
    fDigits = new TClonesArray("AliADSDigit", 16);
    fNdigits = 0;
  }
  return fDigits;
}
