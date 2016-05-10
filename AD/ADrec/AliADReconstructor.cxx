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

/* $Id: AliADReconstructor.cxx 20956 2007-09-26 14:22:18Z mrodrigu $ */
//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//  Class for AD reconstruction                                         //
//////////////////////////////////////////////////////////////////////////////
#include <algorithm>

#include <TParameter.h>
#include <TGeoManager.h>
#include <TGeoMatrix.h>
#include <TGeoPhysicalNode.h>

#include <TF1.h>

#include <AliGeomManager.h>

#include "AliRawReader.h"
#include "AliGRPObject.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliCDBEntry.h"
#include "AliESDEvent.h"
#include "AliRunInfo.h"
#include "AliCTPTimeParams.h"
#include "AliLHCClockPhase.h"
#include "TSpline.h"

#include "AliADReconstructor.h"
#include "AliADdigit.h"
#include "AliESDAD.h"
#include "AliESDADfriend.h"
#include "AliADConst.h"
#include "AliADCalibData.h"
#include "AliADRawStream.h"
#include "AliADDecision.h"

ClassImp(AliADReconstructor);

//_____________________________________________________________________________
AliADReconstructor::AliADReconstructor()
  : AliReconstructor()
  , fESDAD(NULL)
  , fESD(NULL)
  , fESDADfriend(NULL)
  , fCalibData(NULL)
  , fDigitsArray(NULL)
  , fSaturationCorrection(NULL)
  , fCollisionMode(0)
  , fBeamEnergy(0.)
  , fCorrectForSaturation(kTRUE)
{  
  // Default constructor  

  // Get calibration data
  fCalibData = GetCalibData();

  //Get time slewing
  GetTimeSlewingSplines();

  AliCDBManager *man = AliCDBManager::Instance();
  AliCDBEntry *entry = NULL;

  // Now get the CTP L0->L1 delay
  entry = man->Get("GRP/CTP/CTPtiming");
  if (!entry) AliFatal("CTP timing parameters are not found in OCDB !");
  AliCTPTimeParams *ctpParams = dynamic_cast<AliCTPTimeParams*>(entry->GetObject());
  Float_t l1Delay = (Float_t)ctpParams->GetDelayL1L0()*25.0;

  entry = man->Get("GRP/CTP/TimeAlign");
  if (!entry) AliFatal("CTP time-alignment is not found in OCDB !");
  AliCTPTimeParams *ctpTimeAlign = dynamic_cast<AliCTPTimeParams*>(entry->GetObject());
  l1Delay += ((Float_t)ctpTimeAlign->GetDelayL1L0()*25.0);

  entry = man->Get("AD/Calib/TimeDelays");
  if (!entry) AliFatal("AD time delays are not found in OCDB !");
  TH1F *TimeDelays = dynamic_cast<TH1F*>(entry->GetObject());

  entry = man->Get("GRP/Calib/LHCClockPhase");
  if (!entry) AliFatal("LHC clock-phase shift is not found in OCDB !");
  AliLHCClockPhase *phase = dynamic_cast<AliLHCClockPhase*>(entry->GetObject());

  for (Int_t i=0; i<16; ++i) {
    const Int_t board = AliADCalibData::GetBoardNumber(i); 
    fHptdcOffset[i] = (((Float_t)fCalibData->GetRollOver(board)-
			(Float_t)fCalibData->GetTriggerCountOffset(board))*25.0
		       +fCalibData->GetTimeOffset(i)
		       -l1Delay
		       -phase->GetMeanPhase()
		       -TimeDelays->GetBinContent(i+1)
		       -kADOffset);
  }
   
  const Float_t zADC1 = TMath::Abs(GetZPosition("AD/ADC1"));  // outer 
  const Float_t zADC2 = TMath::Abs(GetZPosition("AD/ADC2"));  // inner
  const Float_t zADA1 = TMath::Abs(GetZPosition("AD/ADA1"));  // inner
  const Float_t zADA2 = TMath::Abs(GetZPosition("AD/ADA2"));  // outer

  // distance in time units from nominal vertex to AD
  fLayerDist[0] = zADC1/TMath::Ccgs()*1e9;
  fLayerDist[1] = zADC2/TMath::Ccgs()*1e9;
  fLayerDist[2] = zADA1/TMath::Ccgs()*1e9;
  fLayerDist[3] = zADA2/TMath::Ccgs()*1e9; 
}
//________________________________________________________________________________
void AliADReconstructor::SetOption(Option_t* opt)
{
  TObjArray *oa = TString(opt).Tokenize(" ");
  for (Int_t i=0, n=oa->GetEntries(); i<n; ++i) {
    const TString s = oa->At(i)->GetName();
    if (s.Contains("SaturationCorrection")) {
      fCorrectForSaturation = !s.Contains("-");
    }
  }
  delete oa;
  AliInfo(Form("opt='%s'  fCorrectForSaturation: %s", opt, fCorrectForSaturation ? "ON" : "OFF"));
}
//________________________________________________________________________________
Double_t AliADReconstructor::GetZPosition(const char* symname)
{
  // Get the global z coordinate of the given AD alignable volume
  //
  Double_t *tr;
  TGeoPNEntry *pne = gGeoManager->GetAlignableEntry(symname);
  if (!pne) {
    AliFatalClass(Form("TGeoPNEntry with symbolic name %s does not exist!", symname));
    return 0;
  }

  TGeoPhysicalNode *pnode = pne->GetPhysicalNode();
  if (NULL != pnode) {
    TGeoHMatrix* hm = pnode->GetMatrix();
    tr = hm->GetTranslation();
  } else {
    const char* path = pne->GetTitle();
    if (!gGeoManager->cd(path)) {
      AliFatalClass(Form("Volume path %s not valid!", path));
      return 0;
    }
    tr = gGeoManager->GetCurrentMatrix()->GetTranslation();
  }
  return tr[2];
}


//_____________________________________________________________________________
AliADReconstructor& AliADReconstructor::operator = (const AliADReconstructor&)
{
  // assignment operator
  Fatal("operator =",  "assignment operator not implemented");
  return *this;
}

//_____________________________________________________________________________
AliADReconstructor::~AliADReconstructor()
{
  // destructor
  delete fESDAD;
  delete fESDADfriend;
  delete fDigitsArray;
}

//_____________________________________________________________________________
void AliADReconstructor::Init()
{
  // initializer
  SetOption(GetOption());

  fESDAD        = new AliESDAD;
  fESDADfriend  = new AliESDADfriend;

  if (fCorrectForSaturation) {
    AliCDBEntry *entry = AliCDBManager::Instance()->Get("AD/Calib/Saturation");
    if (!entry) AliFatal("AD/Calib/Saturation is not found in OCDB !");
    fSaturationCorrection = dynamic_cast<TTree*>(entry->GetObject());
  }
}

UShort_t MakeTriggerFlags(Bool_t UBA, Bool_t UBC,
			  Bool_t UGA, Bool_t UGC) {
  UShort_t flags = 0;
  flags |= (1<< 0)*(UBA && UBC);
  flags |= (1<< 1)*(UBA || UBC);
  flags |= (1<< 2)*(UGA && UBC);
  flags |= (1<< 3)*UGA;
  flags |= (1<< 4)*(UGC && UBA);
  flags |= (1<< 5)*UGC;
  //CTA1 and CTC1
  //CTA1 or CTC1
  //CTA2 and CTC2
  //CTA2 or CTC2 
  //MTA and MTC
  //MTA or MTC
  flags |= (1<<12)*UBA;
  flags |= (1<<13)*UBC;
  flags |= (1<<14)*(UGA || UGC);
  flags |= (1<<15)*((UGA && UBC) || (UGC && UBA));
  return flags;
}

//_____________________________________________________________________________
void AliADReconstructor::ConvertDigits(AliRawReader* rawReader, TTree* digitsTree) const
{
  // converts RAW to digits 
  if (NULL == digitsTree) {
    AliError("No digits tree!");
    return;
  }

  if (NULL == fDigitsArray)
    fDigitsArray = new TClonesArray("AliADdigit", 16);
  digitsTree->Branch("ADDigit", &fDigitsArray);

  rawReader->Reset();
  AliADRawStream rawStream(rawReader);
  if (rawStream.Next()) {   
    Bool_t aBBflag[16];
    Bool_t aBGflag[16];
    
    for (Int_t iChannel=0; iChannel<16; ++iChannel) { // iChannel -> online channel number
      const Int_t offlineCh = kOfflineChannel[iChannel];
      // ADC charge samples
      Short_t chargeADC[kADNClocks] = { 0 };
//       aBBflag[offlineCh] = kFALSE;
//       aBGflag[offlineCh] = kFALSE;
      for (Int_t iClock=0; iClock<kADNClocks; ++iClock) {
	chargeADC[iClock]  = rawStream.GetPedestal(iChannel, iClock);
// 	if (rawStream.GetBBFlag(iChannel, iClock)) aBBflag[offlineCh]=kTRUE;
// 	if (rawStream.GetBGFlag(iChannel, iClock)) aBGflag[offlineCh]=kTRUE;
      }
      // Integrator flag
      const Bool_t integrator = rawStream.GetIntegratorFlag(iChannel, kADNClocks/2);
      const Bool_t BBflag     = rawStream.GetBBFlag(iChannel, kADNClocks/2); 
      const Bool_t BGflag     = rawStream.GetBGFlag(iChannel, kADNClocks/2);
      aBBflag[offlineCh] = BBflag;
      aBGflag[offlineCh] = BGflag;
   
      // HPTDC data (leading time and width)
      const Int_t   board     = AliADCalibData::GetBoardNumber(offlineCh);
      const Float_t time      = rawStream.GetTime(iChannel)  * fCalibData->GetTimeResolution(board);
      const Float_t width     = rawStream.GetWidth(iChannel) * fCalibData->GetWidthResolution(board);
      // Add a digit
      //if(!fCalibData->IsChannelDead(offlineCh)){ Off for the moment
      new ((*fDigitsArray)[fDigitsArray->GetEntriesFast()])
	AliADdigit(offlineCh, time, width,integrator, chargeADC, BBflag, BGflag);
      //}
      
      fESDADfriend->SetBBScalers(offlineCh, rawStream.GetBBScalers(iChannel));
      fESDADfriend->SetBGScalers(offlineCh, rawStream.GetBGScalers(iChannel));
      for (Int_t iEv=0; iEv<kADNClocks; ++iEv) {
	fESDADfriend->SetBBFlag(offlineCh, iEv, rawStream.GetBBFlag(iChannel, iEv));
	fESDADfriend->SetBGFlag(offlineCh, iEv, rawStream.GetBGFlag(iChannel, iEv));
      }
    }
    //BC Unmasked triggers
    Int_t pBBmulADA = 0;
    Int_t pBBmulADC = 0;
    Int_t pBGmulADA = 0;
    Int_t pBGmulADC = 0;
  
    for (Int_t iChannel=0; iChannel<4; ++iChannel) {//Loop over pairs of pads
    //Enable time is used to turn off the coincidence 
      pBBmulADC += ((!fCalibData->GetEnableTiming(iChannel+ 0) || aBBflag[iChannel+ 0]) && (!fCalibData->GetEnableTiming(iChannel+ 4) || aBBflag[iChannel+ 4]));
      pBGmulADC += ((!fCalibData->GetEnableTiming(iChannel+ 0) || aBGflag[iChannel+ 0]) && (!fCalibData->GetEnableTiming(iChannel+ 4) || aBGflag[iChannel+ 4]));
      
      pBBmulADA += ((!fCalibData->GetEnableTiming(iChannel+ 8) || aBBflag[iChannel+ 8]) && (!fCalibData->GetEnableTiming(iChannel+12) || aBBflag[iChannel+12]));
      pBGmulADA += ((!fCalibData->GetEnableTiming(iChannel+ 8) || aBGflag[iChannel+ 8]) && (!fCalibData->GetEnableTiming(iChannel+12) || aBGflag[iChannel+12]));
    }
    
    const Bool_t UBA = (pBBmulADA >= fCalibData->GetBBAThreshold());
    const Bool_t UBC = (pBBmulADC >= fCalibData->GetBBCThreshold());
    const Bool_t UGA = (pBGmulADA >= fCalibData->GetBGAThreshold());
    const Bool_t UGC = (pBGmulADC >= fCalibData->GetBGCThreshold());

    const UShort_t fTrigger = MakeTriggerFlags(UBA, UBC,
					       UGA, UGC);
    
    fESDADfriend->SetTriggerInputs(fTrigger);
    fESDADfriend->SetTriggerInputsMask(rawStream.GetTriggerInputsMask());
  
    for (Int_t iScaler=0; iScaler<AliESDADfriend::kNScalers; ++iScaler) {
      fESDADfriend->SetTriggerScalers(iScaler, rawStream.GetTriggerScalers(iScaler));
    }
    digitsTree->Fill();
  }

  fDigitsArray->Clear();
}

//_____________________________________________________________________________
void AliADReconstructor::FillESD(TTree* digitsTree, TTree* /*clustersTree*/,AliESDEvent* esd) const
{
  printf("Running AD Reconstruction \n");

  if (!digitsTree) {
      AliError("No digits tree!");
      return;
  }

  TBranch* digitBranch = digitsTree->GetBranch("ADDigit");
  digitBranch->SetAddress(&fDigitsArray);

  Float_t   mult[16];
  Float_t    adc[16];
  Float_t   tail[16];
  Float_t   tailComplement[16];
  Float_t   time[16];
  Float_t  width[16];
  Bool_t aBBflag[16];
  Bool_t aBGflag[16];
  Float_t   adcTrigger[16];
  
  for (Int_t i=0; i<16; i++){
    adc[i]            = 0.0f;
    tail[i]           = 0.0f;
    tailComplement[i] = 0.0f;
    mult[i]           = 0.0f;
    time[i]           = kInvalidTime;
    width[i]          = 0.0f;
    aBBflag[i]        = kFALSE;
    aBGflag[i]        = kFALSE;
    adcTrigger[i]     = 0.0f;
  }
  
  TClonesArray *f_Int0 = new TClonesArray;
  TClonesArray *f_Int1 = new TClonesArray;
  Int_t chOffline=0, chOnline=0;
  Float_t extrapolationThresholds[21];
  Bool_t  doExtrapolation[21];
  for (Int_t i=0; i<21; ++i) {
    extrapolationThresholds[i] = -999.0f;
    doExtrapolation[i]         = kFALSE;
  }
  Float_t chargeEqualizationFactor = 1.0f;
  if (fCorrectForSaturation) {
    fSaturationCorrection->SetBranchAddress("f_Int0",                   &f_Int0);
    fSaturationCorrection->SetBranchAddress("f_Int1",                   &f_Int1);
    fSaturationCorrection->SetBranchAddress("chOffline",                &chOffline);
    fSaturationCorrection->SetBranchAddress("chOnline",                 &chOnline);
    fSaturationCorrection->SetBranchAddress("extrapolationThresholds",  &extrapolationThresholds);
    fSaturationCorrection->SetBranchAddress("doExtrapolation",          &doExtrapolation);
    fSaturationCorrection->SetBranchAddress("chargeEqualizationFactor", &chargeEqualizationFactor);
  }
  Bool_t correctedForSaturation = kFALSE;

  for (Long64_t e=0, n=digitsTree->GetEntries(); e<n; ++e) {
    digitsTree->GetEvent(e);
    
    for (Int_t d=0, m=fDigitsArray->GetEntriesFast(); d<m; ++d) {    
      const AliADdigit* digit    = dynamic_cast<const AliADdigit*>(fDigitsArray->At(d));
      const Int_t       pmNumber = digit->PMNumber();
      
      // Pedestal retrieval and suppression
      Bool_t  integrator = digit->Integrator();
      Float_t maxadc     = 0.0f;
      Int_t   imax       = -1;
      Float_t adcPedSub[kADNClocks] = { 0 };
      for (Int_t iClock=0; iClock<kADNClocks; ++iClock) {
	const Short_t charge      = digit->ChargeADC(iClock);
	const Bool_t  iIntegrator = ((iClock%2) == 0 ? integrator : !integrator);

	const Int_t k = pmNumber + 16*iIntegrator;
	adcPedSub[iClock] = static_cast<Float_t>(charge) - fCalibData->GetPedestal(k);
	if(adcPedSub[iClock] <= GetRecoParam()->GetNSigmaPed()*fCalibData->GetSigma(k)) {
	  adcPedSub[iClock] = 0;
	  continue;
	}
	
	if (iClock < GetRecoParam()->GetStartClock() ||
	    iClock > GetRecoParam()->GetEndClock()) continue;

	if (adcPedSub[iClock] > maxadc) {
	  maxadc = adcPedSub[iClock];
	  imax   = iClock;
	}
      }

      // start and end BCs for charge integration
      const Int_t start = TMath::Max( 0, imax - GetRecoParam()->GetNPreClocks());
      const Int_t end   = TMath::Min(20, imax + GetRecoParam()->GetNPostClocks());
      
      // integrated charge without saturation correction
      adc[pmNumber] = 0.0f;
      for (Int_t iClock=start; iClock<=end; ++iClock)
	adc[pmNumber] += adcPedSub[iClock];
      
      // HPTDC leading time and width
      // Correction for slewing and various time delays
      time[pmNumber]    = CorrectLeadingTime(pmNumber, digit->Time(), adc[pmNumber]);
      width[pmNumber]   = digit->Width();
      aBBflag[pmNumber] = digit->GetBBflag();
      aBGflag[pmNumber] = digit->GetBGflag();

      if (imax != -1) {
	const Float_t threshold = 20.0f;
	Bool_t isPileUp = kFALSE;
	for (Int_t bc=13; bc<20 && !isPileUp; ++bc)  
	  isPileUp |= (adcPedSub[bc+1] > adcPedSub[bc] + threshold);
    
	for (Int_t iClock=14; iClock<=20; ++iClock)
	  tail[pmNumber] += adcPedSub[iClock];

	adcTrigger[pmNumber] = adcPedSub[10];

	for (Int_t iClock=end+1; iClock<21; ++iClock)
	  tailComplement[pmNumber] += adcPedSub[iClock];

	if (fCorrectForSaturation) {
	  f_Int0->Clear();
	  f_Int1->Clear();
	  fSaturationCorrection->GetEntry(pmNumber);
	}

	AliDebug(3, Form("pmNumber=%d offlineCh=%d onlineCh=%d isPileUp=%d", pmNumber, chOffline, chOnline, isPileUp));
	adc[pmNumber] = 0.0f;
	for (Int_t iClock=start; iClock<=end; ++iClock) {
	  const Bool_t iIntegrator = ((iClock%2) == 0 ? integrator : !integrator);
	  
	  const TF1* fExtrapolation = (doExtrapolation[iClock]
				       ? static_cast<const TF1*>(iIntegrator
								 ? f_Int1->At(iClock)
								 : f_Int0->At(iClock))
				       : NULL);
	  doExtrapolation[iClock] &= (!isPileUp);
	  correctedForSaturation  |= doExtrapolation[iClock];
	  if (doExtrapolation[iClock] && tail[pmNumber] > extrapolationThresholds[iClock]) {
	    AliDebug(3, Form("extrapolation[%2d] clock=%2d adcBefore=%.1f adcAfter=%.1f",
			     pmNumber, iClock, adcPedSub[iClock], fExtrapolation->Eval(tail[pmNumber])));
	  }

	  adc[pmNumber] += chargeEqualizationFactor*((doExtrapolation[iClock] && tail[pmNumber] > extrapolationThresholds[iClock])
						     ? fExtrapolation->Eval(tail[pmNumber])
						     : adcPedSub[iClock]
						     );
	}
      }

      AliDebug(3, Form("ADreco: GetRecoParam()->GetNPreClocks()=%d, GetRecoParam()->GetNPostClocks()=%d imax=%d maxadc=%.1f adc[%2d]=%.1f tail[%2d]=%.1f tailC=%.1f",
		       GetRecoParam()->GetNPreClocks(),
		       GetRecoParam()->GetNPostClocks(),
		       imax,
		       maxadc,
		       pmNumber,
		       adc[pmNumber],
		       pmNumber,
		       tail[pmNumber],
		       tailComplement[pmNumber]));
      
      // MIP multiplicity
      const Float_t adcPerMIP = const_cast<AliADCalibData*>(fCalibData)->GetADCperMIP(pmNumber);
      mult[pmNumber] = (adcPerMIP != 0 ? adc[pmNumber]/adcPerMIP : 0.0f);
      
      if (adc[pmNumber] > 0) {
	AliDebug(1, Form("PM = %d ADC = %.2f (%.2f) TDC %.2f (%.2f)   Int %d (%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d)    %.2f %.2f   %.2f %.2f    %d %d",
			 pmNumber,  adc[pmNumber],
			 digit->ChargeADC(11)+digit->ChargeADC(10)+digit->ChargeADC(9)+digit->ChargeADC(8)+
			 digit->ChargeADC(7)+digit->ChargeADC(6)+digit->ChargeADC(5)+digit->ChargeADC(4)-
			 4.*fCalibData->GetPedestal(pmNumber)-4.*fCalibData->GetPedestal(pmNumber+16),
			 digit->Time(), time[pmNumber],
			 integrator,
			 digit->ChargeADC(0), digit->ChargeADC(1), digit->ChargeADC(2), digit->ChargeADC(3), digit->ChargeADC(4), digit->ChargeADC(5), digit->ChargeADC(6), digit->ChargeADC(7),
			 digit->ChargeADC(8), digit->ChargeADC(9), digit->ChargeADC(10),
			 digit->ChargeADC(11), digit->ChargeADC(12),
			 digit->ChargeADC(13), digit->ChargeADC(14), digit->ChargeADC(15), digit->ChargeADC(16), digit->ChargeADC(17), digit->ChargeADC(18), digit->ChargeADC(19), digit->ChargeADC(20),
			 fCalibData->GetPedestal(pmNumber), fCalibData->GetSigma(pmNumber),
			 fCalibData->GetPedestal(pmNumber+16), fCalibData->GetSigma(pmNumber+16),
			 aBBflag[pmNumber], aBGflag[pmNumber]));
      };
      
      // Fill ESD friend object
      for (Int_t iEv=0; iEv < kADNClocks; ++iEv) {
	fESDADfriend->SetPedestal      (pmNumber, iEv, static_cast<Float_t>(digit->ChargeADC(iEv)));
	fESDADfriend->SetIntegratorFlag(pmNumber, iEv, ((iEv%2) == 0 ? integrator : !integrator));
      }
      fESDADfriend->SetTime (pmNumber, digit->Time());
      fESDADfriend->SetWidth(pmNumber, digit->Width());      
    } // end of loop over digits
  } // end of loop over events in digits tree

  if (fCorrectForSaturation)
    fSaturationCorrection->ResetBranchAddresses();
  
  fESDAD->SetBit(AliESDAD::kCorrectedLeadingTime,   kTRUE);
  fESDAD->SetBit(AliESDAD::kCorrectedForSaturation, correctedForSaturation);
  fESDAD->SetMultiplicity(mult);
  fESDAD->SetADC(adc);
  fESDAD->SetTime(time);
  fESDAD->SetWidth(width);
  fESDAD->SetBBFlag(aBBflag);
  fESDAD->SetBGFlag(aBGflag);
  fESDAD->SetADCTail(tail);
  fESDAD->SetADCTrigger(adcTrigger);
  
  // Fill BB and BG flags for all channel in 21 clocks (called past-future flags)
  for (Int_t i=0; i<16; ++i) {
    for (Int_t iClock=0; iClock<kADNClocks; ++iClock) {
      fESDAD->SetPFBBFlag(i, iClock, fESDADfriend->GetBBFlag(i, iClock));
      fESDAD->SetPFBGFlag(i, iClock, fESDADfriend->GetBGFlag(i, iClock));
    }
  }
  fESDAD->SetBit(AliESDAD::kPastFutureFlagsFilled, kTRUE);
   
  Int_t	pBBmulADA = 0;
  Int_t	pBBmulADC = 0;
  Int_t	pBGmulADA = 0;
  Int_t	pBGmulADC = 0;
  
  for (Int_t iChannel=0; iChannel<4; ++iChannel) {//Loop over pairs of pads
    //Enable time is used to turn off the coincidence 
    pBBmulADC += ((!fCalibData->GetEnableTiming(iChannel+ 0) || aBBflag[iChannel+ 0]) && (!fCalibData->GetEnableTiming(iChannel+ 4) || aBBflag[iChannel+ 4]));
    pBGmulADC += ((!fCalibData->GetEnableTiming(iChannel+ 0) || aBGflag[iChannel+ 0]) && (!fCalibData->GetEnableTiming(iChannel+ 4) || aBGflag[iChannel+ 4]));

    pBBmulADA += ((!fCalibData->GetEnableTiming(iChannel+ 8) || aBBflag[iChannel+ 8]) && (!fCalibData->GetEnableTiming(iChannel+12) || aBBflag[iChannel+12]));
    pBGmulADA += ((!fCalibData->GetEnableTiming(iChannel+ 8) || aBGflag[iChannel+ 8]) && (!fCalibData->GetEnableTiming(iChannel+12) || aBGflag[iChannel+12]));
  }
	
  const Bool_t UBA = (pBBmulADA >= fCalibData->GetBBAThreshold());
  const Bool_t UBC = (pBBmulADC >= fCalibData->GetBBCThreshold());
  const Bool_t UGA = (pBGmulADA >= fCalibData->GetBGAThreshold());
  const Bool_t UGC = (pBGmulADC >= fCalibData->GetBGCThreshold());

  const UShort_t fTrigger = MakeTriggerFlags(UBA, UBC,
					     UGA, UGC);

  fESDAD->SetTriggerBits(fTrigger);
  fESDAD->SetBit(AliESDAD::kOnlineBitsFilled,kTRUE);

  // now fill the AD decision
  AliADDecision offlineDecision;
  offlineDecision.SetRecoParam(GetRecoParam());
  offlineDecision.FillDecisions(fESDAD);
  
  if (NULL != esd) { 
    AliDebug(1, Form("Writing AD data to ESD tree"));
    esd->SetADData(fESDAD);
    
    AliESDfriend *fr = dynamic_cast<AliESDfriend*>(esd->FindListObject("AliESDfriend"));
    if (NULL != fr) {
      AliDebug(1, Form("Writing AD friend data to ESD tree"));
      fr->SetADfriend(fESDADfriend);
    }
  }
  
  fDigitsArray->Clear();
}

//_____________________________________________________________________________
const AliADCalibData* AliADReconstructor::GetCalibData() const
{
  AliCDBManager *man = AliCDBManager::Instance();
  AliCDBEntry *entry = man->Get("AD/Calib/Data");
  const AliADCalibData *calibdata = NULL;

  if (NULL != entry)
    calibdata = dynamic_cast<const AliADCalibData*>(entry->GetObject());
  if (NULL == calibdata)
    AliFatal("No calibration data from calibration database !");

  return calibdata;
}

//_____________________________________________________________________________
void AliADReconstructor::GetTimeSlewingSplines()
{
  AliCDBManager *man = AliCDBManager::Instance();
  AliCDBEntry *entry = man->Get("AD/Calib/TimeSlewing");  
  TList *fListSplines = NULL;

  if (NULL != entry)
    fListSplines = dynamic_cast<TList*>(entry->GetObject());
  if (NULL == fListSplines)
    AliFatal("No time slewing correction from calibration database !");
  
  for (Int_t i=0; i<16; ++i)
    fTimeSlewingSpline[i] = dynamic_cast<TSpline3*>(fListSplines->At(i));
}

//_____________________________________________________________________________
AliCDBStorage* AliADReconstructor::SetStorage(const char *uri) 
{
  // Sets the storage  
  Bool_t deleteManager = kFALSE;  
  AliCDBManager *manager = AliCDBManager::Instance();
  AliCDBStorage *defstorage = manager->GetDefaultStorage();
  
  if(!defstorage || !(defstorage->Contains("AD"))) { 
     AliWarning("No default storage set or default storage doesn't contain AD!");
     manager->SetDefaultStorage(uri);
     deleteManager = kTRUE;
  }
 
  AliCDBStorage *storage = manager->GetDefaultStorage();
  if(deleteManager) {
     AliCDBManager::Instance()->UnsetDefaultStorage();
     defstorage = 0;   // the storage is killed by AliCDBManager::Instance()->Destroy()
  }

  return storage; 
}

//____________________________________________________________________________
void AliADReconstructor::GetCollisionMode()
{
  // Retrieval of collision mode 
  const TString beamType = GetRunInfo()->GetBeamType();
  if (beamType == AliGRPObject::GetInvalidString()) {
    AliError("AD cannot retrieve beam type");
    return;
  }
  if((beamType.CompareTo("P-P") ==0) ||
     (beamType.CompareTo("p-p") ==0)) {
    fCollisionMode = 0;
  } else if ((beamType.CompareTo("Pb-Pb") == 0)  ||
	     (beamType.CompareTo("A-A") == 0)) {
    fCollisionMode = 1;
  }
    
  fBeamEnergy = GetRunInfo()->GetBeamEnergy();
  if (fBeamEnergy == AliGRPObject::GetInvalidFloat()) {
    AliError("Missing value for the beam energy ! Using 0");
    fBeamEnergy = 0.;
  }
  
  AliDebug(1, Form("\n ++++++ Beam type and collision mode retrieved as %s %d @ %1.3f GeV ++++++\n\n", beamType.Data(), fCollisionMode, fBeamEnergy));
}
//____________________________________________________________________________
Float_t AliADReconstructor::CorrectLeadingTime(Int_t i, Float_t time, Float_t adc) const
{
  // Correct the leading time
  // for slewing effect and
  // misalignment of the channels
  if (time < 1e-6)
    return kInvalidTime;

  // In case of pathological signals
  if (adc < 1.0)
    return time;

  // Slewing and offset correction
  const Int_t board = AliADCalibData::GetBoardNumber(i);
  time -= fTimeSlewingSpline[i]->Eval(TMath::Log10(1/adc))*fCalibData->GetTimeResolution(board);
  time += fLayerDist[i/4];
  
  // Channel alignment and general offset subtraction
  //time -= fHptdcOffset[i];
  return time;
}
