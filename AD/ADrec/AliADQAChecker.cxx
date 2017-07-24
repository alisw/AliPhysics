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


/*
  Checks the quality assurance. Under construction.
  By comparing with reference data

*/

// --- ROOT system ---
#include <TH2.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TPaveText.h>
#include <TLatex.h>
#include <TString.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TSpline.h>

// --- Standard library ---

// --- AliRoot header files ---
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliCDBEntry.h"
#include "AliLog.h"
#include "AliQAv1.h"
#include "AliQAChecker.h"
#include "AliADQAChecker.h"
#include "AliADQADataMakerRec.h"
#include "AliADQAParam.h"

ClassImp(AliADQAChecker);

//__________________________________________________________________
AliADQAChecker::AliADQAChecker()
  : AliQACheckerBase("AD","AD Quality Assurance Data Checker")
  , fQAParam(NULL)
  , fQAParamDefault()
  , fLowEventCut(1000)
  , fORvsANDCut(0.2)
  , fBGvsBBCut(0.2)
{
  fQAParamDefault.SetSatMed(0.1);
  fQAParamDefault.SetSatHigh(0.3);
  fQAParamDefault.SetSatHuge(0.5);
  fQAParamDefault.SetMaxPedDiff(1);
  fQAParamDefault.SetMaxPedWidth(1.5);
  fQAParamDefault.SetChargeChannelZoomMin(0);
  fQAParamDefault.SetChargeChannelZoomMax(50);
  fQAParamDefault.SetTdcTimeMinBBFlag(170);
  fQAParamDefault.SetTdcTimeMaxBBFlag(210);
  fQAParamDefault.SetTdcTimeMinBGFlag( 50);
  fQAParamDefault.SetTdcTimeMaxBGFlag( 90);
  fQAParamDefault.SetChargeTrendMin(0);
  fQAParamDefault.SetChargeTrendMax(1000);
  fQAParamDefault.SetMaxNoTimeRate(10e-3);
  fQAParamDefault.SetMaxNoFlagRate(10e-2);
  fQAParamDefault.SetMaxBBVariation(5);
  fQAParamDefault.SetMaxBGVariation(5);
  fQAParamDefault.SetAsynchronBB(0.5);
  fQAParamDefault.SetAsynchronBG(0.5);
  fQAParam = &fQAParamDefault;
}
//____________________________________________________________________________
AliADQAChecker::~AliADQAChecker()
{
}

//____________________________________________________________________________
const AliADQAParam* AliADQAChecker::GetQAParam() const
{
  AliCDBManager *man = AliCDBManager::Instance();
  AliCDBEntry *entry = man->Get("AD/Calib/QAParam");
  if (!entry) {
    AliWarning("Load of QA param from default storage failed!");
  } else {
    const AliADQAParam *QAParam = dynamic_cast<const AliADQAParam*>(entry->GetObject());
    if (!QAParam)
      AliError("No QA param from calibration database !");
    return QAParam;
  }
  // fall-back:
  return &fQAParamDefault;
}
//__________________________________________________________________
void AliADQAChecker::Check(Double_t *check, AliQAv1::ALITASK_t index, TObjArray **list, const AliDetectorRecoParam * /*recoParam*/)
{
  // Main check function: Depending on the TASK, different checks will be applied
  // Check for missing channels and check on the trigger type for raw data
  // Check for missing disk or rings for esd (to be redone)

  for (Int_t specie=0; specie<AliRecoParam::kNSpecies; ++specie) {
    check[specie] = 1.0;
    if (!AliQAv1::Instance()->IsEventSpecieSet(specie))
      continue;
    if (index == AliQAv1::kRAW) {
      if (AliRecoParam::ConvertIndex(specie) == AliRecoParam::kCalib) {
	check[specie] = CheckPedestals(list[specie]);
      } else {
	check[specie] = CheckRaws(list[specie]);
      }
    } else if (index == AliQAv1::kESD) {
      check[specie] =  CheckEsds(list[specie]);
    }
  }
}

// returns the QA box of a given histogram, if it does not exists it is created
TPaveText* GetQABox(TH1* h, Float_t x1, Float_t y1, Float_t x2, Float_t y2, Option_t *opt) {
  TList *lHistFcn = h->GetListOfFunctions();
  const TString nameOfQABox = TString::Format("QABox_%s", h->GetName());
  TPaveText *qaBox = dynamic_cast<TPaveText*>(lHistFcn->FindObject(nameOfQABox));
  if (!qaBox) {
    qaBox = new TPaveText(x1, y1, x2, y2, opt);
    qaBox->SetName(nameOfQABox);
    lHistFcn->Add(qaBox);
  }
  return qaBox;
}

//_________________________________________________________________
Double_t AliADQAChecker::CheckPedestals(TObjArray* list) const
{
  Double_t test = 1.0;
  if (!list || !list->GetEntries())
    return test;

  const Int_t idx[2] = {
    AliADQADataMakerRec::kPedestalInt0,
    AliADQADataMakerRec::kPedestalInt1
  };
  for (Int_t integrator=0; integrator<2; ++integrator) {
    TH2 *hPed = dynamic_cast<TH2*>(list->At(idx[integrator]));
    if (!hPed) {
      AliWarningF("PedestalInt%d histogram is not found", integrator);
      continue;
    }

    // add the QA box (=TPaveText) if it does not exist
    TPaveText *qaBox = GetQABox(hPed, 0.15, 0.63, 0.85, 0.85, "NDC");

    // set up the QA box
    qaBox->Clear();
    if (!hPed->GetEntries()) {
      qaBox->SetFillColor(kYellow);
      qaBox->AddText("Is this AD pedestal run?");
    } else {
      TString badChannels = "";
      for (Int_t ch=0; ch<16; ch++) {
	TH1* hPedSlice = hPed->ProjectionY("hPedSlice", ch+1, ch+1);
	if(hPedSlice->GetRMS() > fQAParam->GetMaxPedWidth())
	  badChannels += TString::Format("%d, ", ch);
	delete hPedSlice;
      } // next channel
      if (badChannels != "") {
	test = 0.1;
	qaBox->SetFillColor(kRed);
	qaBox->AddText("Noisy pedestal for channel " + badChannels);
      } else {
	qaBox->SetFillColor(kGreen);
	qaBox->AddText("OK");
      }
    }
  } // next integrator
  return test;
}

//_________________________________________________________________
Double_t AliADQAChecker::CheckRaws(TObjArray* list) const
{
  Double_t test = 1.0;
  if (!list || !list->GetEntries()) {
    AliWarning("There are no histograms to be checked");
    return test;
  }

  // (1) charge ADA,ADC
  const Int_t idx[2] = {
    AliADQADataMakerRec::kChargeADA,
    AliADQADataMakerRec::kChargeADC
  };
  const char*   sideName[2]   = { "ADA", "ADC" };
  const Float_t paveTextX1[2] = {  0.43, 0.43 };
  const Float_t paveTextX2[2] = {  0.66, 0.66 };
  const Float_t paveTextY1[2] = {  0.70, 0.55 };
  const Float_t paveTextY2[2] = {  0.83, 0.63 };
  const Int_t   textColor[2]  = { kRed, kBlue };

  Float_t nEvents = 0.0f; // needed later
  for (Int_t side=0; side<2; ++side) {
    TH1 *hCharge = dynamic_cast<TH1*>(list->At(idx[side]));
    if (!hCharge) {
      AliWarningF("Charge%s histogram is not found", sideName[side]);
      continue;
    }
    if (side == 0) // take nEvents from A-side charge histogram
      nEvents = hCharge->Integral(-1,-1);

    TPaveText *qaBox = GetQABox(hCharge, paveTextX1[side], paveTextY1[side], paveTextX2[side], paveTextY2[side], "NDC");
    qaBox->Clear();
    qaBox->SetFillColor(kWhite);
    qaBox->SetTextColor(textColor[side]);
    qaBox->AddText(TString::Format(" Mean charge %s = %4.1f", sideName[side], hCharge->GetMean()));
  }

  // (2) trigger counts
  TH1* hTriggers = dynamic_cast<TH1*>(list->At(AliADQADataMakerRec::kTriggers));
  Int_t nEventsADOR = 0;
  if (!hTriggers) {
    AliWarning("Triggers histogram is not found");
  } else {
    nEventsADOR = TMath::Nint(hTriggers->GetBinContent(6));
    TPaveText *qaBox = GetQABox(hTriggers, 0.56, 0.73, 0.72, 0.85, "NDC");
    qaBox->Clear();
    qaBox->SetFillColor(kWhite);
    qaBox->SetTextColor(kAzure-8);
    qaBox->AddText(TString::Format(" Number of events = %.0f", nEvents));
  }

  // (3) flag and no time
  TH1 *hFlagNoTime = dynamic_cast<TH1*>(list->At(AliADQADataMakerRec::kFlagNoTime));
  if (!hFlagNoTime) {
    AliWarning("FlagNoTime histogram is not found");
  } else {
    TPaveText *qaBox = GetQABox(hFlagNoTime, 0.15, 0.67, 0.85, 0.76, "NDC");

    TH1 *histoRate = dynamic_cast<TH1*>(hFlagNoTime->Clone("histoRateCloned"));
    histoRate->Sumw2();
    if (nEvents)
      histoRate->Scale(1.0f/nEvents);

    TString badChannels = "";
    for (Int_t ch=0; ch<16; ++ch) {
      if (histoRate->GetBinContent(ch+1) > fQAParam->GetMaxNoTimeRate()) {
	test = 0.3;
	badChannels += TString::Format("%d, ", ch);
      }
    }
    delete histoRate;
	qaBox->Clear();
      if (badChannels != "") {
	qaBox->SetFillColor(kRed);
	qaBox->AddText("No time rate too high, ch:" + badChannels);
      } else {
	qaBox->SetFillColor(kGreen);
	qaBox->AddText("No time rate OK");
    }
  }

  // (4) time and no flag
  TH1 *hTimeNoFlag = dynamic_cast<TH1*>(list->At(AliADQADataMakerRec::kTimeNoFlag));
  if (!hTimeNoFlag) {
    AliWarning("TimeNoFlag histogram is not found");
  } else {
    TPaveText *qaBox = GetQABox(hTimeNoFlag, 0.15, 0.56, 0.85, 0.65, "NDC");
    TH1 *histoRate = dynamic_cast<TH1*>(hTimeNoFlag->Clone("histoRateCloned"));
    histoRate->Sumw2();
    if (nEvents)
      histoRate->Scale(1.0f/nEvents);

    TString badChannels = "";
    for (Int_t ch=0; ch<16; ++ch) {
      if(histoRate->GetBinContent(ch+1) > fQAParam->GetMaxNoFlagRate()) {
	test = 0.3;
	badChannels += TString::Format("%d, ", ch);
      }
    }
    delete histoRate;

    qaBox->Clear();
    if (badChannels != "") {
      qaBox->SetFillColor(kRed);
      qaBox->AddText("No flag rate too high, ch: " + badChannels);
    } else {
      qaBox->SetFillColor(kGreen);
      qaBox->AddText("No flag rate OK");
    }
  }

  // (4) BB and BG flag
  const Int_t idxFlag[2] = {
    AliADQADataMakerRec::kNEventsBBFlag,
    AliADQADataMakerRec::kNEventsBGFlag
  };
  const char* flagNames[2] = { "BB", "BG" };
  const Float_t y1[2] = { 0.67, 0.56 };
  const Float_t y2[2] = { 0.76, 0.65 };
  const Float_t limitsFlags[2] = {
    fQAParam->GetMaxBBVariation(),
    fQAParam->GetMaxBGVariation()
  };
  for (Int_t flagType=0; flagType<2; ++flagType) {
    TH1* hFlag = dynamic_cast<TH1*>(list->At(idxFlag[flagType]));
    if (!hFlag) {
      AliWarningF("NEvents%sFlag histogram is not found", flagNames[flagType]);
      continue;
    }
    TPaveText *qaBox = GetQABox(hFlag, 0.15, y1[flagType], 0.85, y2[flagType], "NDC");
    TH1 *histoRate = dynamic_cast<TH1*>(hFlag->Clone("histoRateCloned"));
    histoRate->Sumw2();
    if (nEvents)
      histoRate->Scale(1.0f/nEvents);

    Float_t meanRates[2] = { 0, 0 };
    for (Int_t ch=0; ch<16; ++ch)
      meanRates[ch/8] += histoRate->GetBinContent(1+ch);

    meanRates[0] /= 8;
    meanRates[1] /= 8;

    Bool_t highVar = kFALSE;
    for (Int_t ch=0; ch<16 && !highVar; ++ch) {
      if (TMath::Abs(histoRate->GetBinContent(1+ch) - meanRates[ch/8]) > limitsFlags[flagType]) {
	test    = 0.7;
	highVar = kTRUE;
      }
    }
    delete histoRate;

    qaBox->Clear();
    if (highVar) {
      qaBox->SetFillColor(kYellow);
      qaBox->AddText(TString::Format("%s rate variation too high", flagNames[flagType]));
    } else {
      qaBox->SetFillColor(kGreen);
      qaBox->AddText(TString::Format("%s rate variation OK", flagNames[flagType]));
    }
  } // next flag type

  // (5) charge saturation
  TH2 *hChargeSaturation = dynamic_cast<TH2*>(list->At(AliADQADataMakerRec::kChargeSaturation));
  if (!hChargeSaturation) {
    AliWarning("ChargeSaturation histogram is not found");
  } else {
    TPaveText *qaBox = GetQABox(hChargeSaturation, 0.11, 0.40, 0.89, 0.65, "NDC");
    TString satText = "    ";
    Bool_t medSat  = kFALSE;
    Bool_t highSat = kFALSE;
    Bool_t hugeSat = kFALSE;
    for (Int_t ch=0; ch<16; ++ch) {
      TH1* hChargeSlice = hChargeSaturation->ProjectionY("hChargeSlice", ch+1, ch+1);
      const Double_t saturation = (hChargeSlice->Integral()
				   ? hChargeSlice->Integral(1000,1025)/hChargeSlice->Integral()
				   : 0.0);
      satText += TString::Format("%1.3f   ", saturation);
      if (saturation >= fQAParam->GetSatMed() && saturation < fQAParam->GetSatHigh()) {
	test   = 0.7;
	medSat = kTRUE;
      }
      if (saturation >= fQAParam->GetSatHigh() && saturation < fQAParam->GetSatHuge()) {
	test    = 0.3;
	highSat = kTRUE;
      }
      if (saturation > fQAParam->GetSatHuge()) {
	test    = 0.1;
	hugeSat = kTRUE;
      }
      delete hChargeSlice;
    } // next channel

    qaBox->Clear();
    qaBox->AddText("Saturation");
    qaBox->AddText(satText);
    if (hugeSat) {
      qaBox->SetFillColor(kRed);
    } else if (highSat) {
      qaBox->SetFillColor(kRed);
    } else if (medSat) {
      qaBox->SetFillColor(kOrange);
    } else {
      qaBox->SetFillColor(kGreen);
    }
  }

  // (6) BB flags vs clock
  const Int_t idxFlagVsClock[2] = {
    AliADQADataMakerRec::kBBFlagVsClock,
    AliADQADataMakerRec::kBBFlagVsClock_ADOR
  };
  const char* flagVsClockNames[2] = {
    "BBFlagVsClock",
    "BBFlagVsClock_ADOR"
  };
  const Bool_t isLowStatistics[2] = {
    kFALSE,
    nEventsADOR<50
  };
  for (Int_t flagType=0; flagType<2; ++flagType) {
    TH2 *hFlag = dynamic_cast<TH2*>(list->At(idxFlagVsClock[flagType]));
    if (!hFlag) {
      AliWarningF("%s histogram is not found", flagVsClockNames[flagType]);
      continue;
    }
    TPaveText *qaBox = GetQABox(hFlag, 0.30, 0.15, 0.70, 0.37, "NDC");
    if (isLowStatistics[flagType]) {
      qaBox->Clear();
      qaBox->SetFillColor(kYellow);
      qaBox->AddText("Low statistics");
      continue;
    }
    Bool_t notConfig[2] = { kFALSE, kFALSE };
    Bool_t notSynch[2]  = { kFALSE, kFALSE };
    for (Int_t ch=0; ch<16; ++ch) {
      TH1* hClockSlice = hFlag->ProjectionY("hClockSlice", ch+1, ch+1);
      Double_t center = hClockSlice->GetBinContent(11);
      Double_t flagArray[21] = { 0 };
      for (Int_t iClock = 0; iClock<21; ++iClock)
	flagArray[iClock] = hClockSlice->GetBinContent(iClock+1);

      const Int_t maxClock = TMath::LocMax(21, flagArray);
      if (center == 0) {
	notConfig[ch/8] = kTRUE;
      } else if (maxClock != 10) {
	notSynch[ch/8] = kTRUE;
      }
      delete hClockSlice;
    } // next channel

    qaBox->Clear();
    if (notConfig[0] || notConfig[1]) {
      qaBox->SetFillColor(kViolet);
      qaBox->AddText(notConfig[1]
		     ? "ADA: dead channels!"
		     : "ADA: ok");
      qaBox->AddText(notConfig[0]
		     ? "ADC: dead channels!"
		     : "ADC: ok");
    } else if (notSynch[0] || notSynch[1]) {
      qaBox->SetFillColor(kRed);
      qaBox->AddText(notSynch[1]
		     ? "ADA: misconfigured!"
		     : "ADA: ok");
      qaBox->AddText(notSynch[0]
		     ? "ADC: misconfigured!"
		     : "ADC: ok");
    } else {
      qaBox->SetFillColor(kGreen);
      qaBox->AddText("ADA ok");
      qaBox->AddText("ADC ok");
    }
  } // next flag type

  // (7) pedestal differences
  const Int_t idxPedDiff[2] = {
    AliADQADataMakerRec::kPedestalDiffInt0,
    AliADQADataMakerRec::kPedestalDiffInt1
  };
  for (Int_t integrator=0; integrator<2; ++integrator) {
    TH2* hPedDiff = dynamic_cast<TH2*>(list->At(idxPedDiff[integrator]));
    if (!hPedDiff) {
      AliWarningF("PedestalInt%d histogram is not found", integrator);
      continue;
    }
    TPaveText *qaBox = GetQABox(hPedDiff, 0.15, 0.63, 0.85, 0.85, "NDC");

    TString badChannels = "";
    for (Int_t ch=0; ch<16; ++ch) {
      TH1* hPedestalSlice = hPedDiff->ProjectionY("hPedestalSlice", ch+1, ch+1);
      const Double_t mean = hPedestalSlice->GetMean();
      if (TMath::Abs(mean) > fQAParam->GetMaxPedDiff()) {
	test = 0.3;
	badChannels += TString::Format(" %d,", ch);
      }
      delete hPedestalSlice;
    } // next channel

    qaBox->Clear();
    if (badChannels == "") {
      qaBox->SetFillColor(kGreen);
      qaBox->AddText("OK");
    } else {
      qaBox->SetFillColor(kYellow);
      qaBox->AddText("Unstable pedestal for channel");
      qaBox->AddText(badChannels);
    }
  } // next integrator

  return test;
}

//_________________________________________________________________
Double_t AliADQAChecker::CheckEsds(TObjArray * list) const
{
  const Double_t test = 1.0; // initialisation to OK
  return test;
}

//______________________________________________________________________________
void AliADQAChecker::Init(const AliQAv1::DETECTORINDEX_t det)
{
  // intialises QA and QA checker settings
  fQAParam = GetQAParam();

  AliQAv1::Instance(det);
  Float_t *hiValue  = new Float_t[AliQAv1::kNBIT];
  Float_t *lowValue = new Float_t[AliQAv1::kNBIT];
  lowValue[AliQAv1::kINFO]    =  0.5;
  hiValue[AliQAv1::kINFO]     =  1.0;
  lowValue[AliQAv1::kWARNING] =  0.2;
  hiValue[AliQAv1::kWARNING]  =  0.5;
  lowValue[AliQAv1::kERROR]   =  0.0;
  hiValue[AliQAv1::kERROR]    =  0.2;
  lowValue[AliQAv1::kFATAL]   = -1.0;
  hiValue[AliQAv1::kFATAL]    =  0.0;
  SetHiLo(hiValue, lowValue);
  delete [] hiValue;
  delete [] lowValue;
}

//______________________________________________________________________________
void AliADQAChecker::SetQA(AliQAv1::ALITASK_t index, Double_t * value) const
{
  // sets the QA word according to return value of the Check
  AliQAv1 * qa = AliQAv1::Instance(index);
  for (Int_t specie = 0; specie < AliRecoParam::kNSpecies; specie++) {
    qa->UnSet(AliQAv1::kFATAL, specie);
    qa->UnSet(AliQAv1::kWARNING, specie);
    qa->UnSet(AliQAv1::kERROR, specie);
    qa->UnSet(AliQAv1::kINFO, specie);
    if (! value ) { // No checker is implemented, set all QA to Fatal
      qa->Set(AliQAv1::kFATAL, specie);
    } else {
      if ( value[specie] >= fLowTestValue[AliQAv1::kFATAL] && value[specie] < fUpTestValue[AliQAv1::kFATAL] )
        qa->Set(AliQAv1::kFATAL, specie);
      else if ( value[specie] > fLowTestValue[AliQAv1::kERROR] && value[specie] <= fUpTestValue[AliQAv1::kERROR]  )
        qa->Set(AliQAv1::kERROR, specie);
      else if ( value[specie] > fLowTestValue[AliQAv1::kWARNING] && value[specie] <= fUpTestValue[AliQAv1::kWARNING]  )
        qa->Set(AliQAv1::kWARNING, specie);
      else if ( value[specie] > fLowTestValue[AliQAv1::kINFO] && value[specie] <= fUpTestValue[AliQAv1::kINFO] )
        qa->Set(AliQAv1::kINFO, specie);
    }
  }
}

// helper functions for drawing histograms H1-H3
// (H1)
TLegend* MakeLegend(Double_t x1, Double_t y1, Double_t x2, Double_t y2, Double_t textSize) {
  TLegend *leg = new TLegend(x1, y1, x2, y2);
  leg->SetBit(TObject::kCanDelete); // will be deleted by TCanvas::Clear();
  leg->SetTextFont(42);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetFillColor(0);
  leg->SetMargin(0.25);
  leg->SetTextSize(textSize);
  leg->SetEntrySeparation(0.5);
  return leg;
}
// (H2)
TH1* MakeRatioHistogram(TH1* h1, TH1* h2) { // returns h1 / h2
  TH1* h = dynamic_cast<TH1*>(h1->Clone(Form("Ratio_%s_%s", h1->GetName(), h2->GetName())));
  h->Divide(h2);
  h->SetBit(TObject::kCanDelete); // will be deleted by TCanvas::Clear();
  return h;
}
// (H3)
TH1* MakeScaledHistogram(TH1 *h1, Double_t norm, Bool_t sumw2=kFALSE) { // returns h1 / norm
  TH1 *h = dynamic_cast<TH1*>(h1->Clone(Form("Scaled_%s", h1->GetName())));
  h->Sumw2(sumw2);
  h->SetBit(TObject::kCanDelete); // will be deleted by TCanvas::Clear();
  if (norm)
    h->Scale(1.0/norm);
  return h;
}

TCanvas* AliADQAChecker::CreatePads(TCanvas *c1) const {
  gStyle->SetLabelSize  ( 1.1,  "XYZ");
  gStyle->SetLabelFont  (42,    "XYZ");
  gStyle->SetLabelOffset( 0.01, "XYZ");
  gStyle->SetTitleFont  (42,    "XYZ");
  gStyle->SetTitleOffset( 1.0,  "XYZ");
  gStyle->SetTitleSize  ( 1.1,  "XYZ");

  const UInt_t   nRows       = 11;
  const UInt_t   nLargeRows  =  1;
  const Bool_t   isLarge[11] = { 0,0,1,0,0,0,0,0,0,0,0 };
  const Double_t timesLarger = 1.5;
  Double_t xRow[nRows] = { 0.95, 0,0,0,0,0,0,0,0,0,0 };
  const Double_t xStep = 0.95/(nRows - nLargeRows + timesLarger*nLargeRows);
  for (UInt_t iRow=1; iRow<nRows; ++iRow)
    xRow[iRow] = xRow[iRow-1] - xStep*(isLarge[iRow-1] ? timesLarger : 1.0);

  const Int_t nPads = 12;
  struct PadInfo {
    PadInfo(const char* name, Int_t divX, Double_t x1, Double_t y1, Double_t x2, Double_t y2)
      : fName(name)
      , fDivX(divX)
      , fX1(x1)
      , fY1(y1)
      , fX2(x2)
      , fY2(y2) {}

    TString fName;
    Int_t   fDivX;
    Double_t fX1;
    Double_t fY1;
    Double_t fX2;
    Double_t fY2;
  };

  const struct PadInfo padInfo[12] = {
    PadInfo("Charge",        2, 0.0, xRow[ 1], 1.0, xRow[ 0]),
    PadInfo("ChargeZoom",    6, 0.0, xRow[ 2], 1.0, xRow[ 1]),
    PadInfo("Time",          4, 0.0, xRow[ 3], 1.0, xRow[ 2]),
    PadInfo("TimeRatio",     4, 0.0, xRow[ 4], 1.0, xRow[ 3]),
    PadInfo("MeanTime",      4, 0.0, xRow[ 5], 1.0, xRow[ 4]),
    PadInfo("ClockCfg",      5, 0.0, xRow[ 6], 1.0, xRow[ 5]),
    PadInfo("Pedestal",      2, 0.0, xRow[ 7], 0.5, xRow[ 6]),
    PadInfo("MaxCharge",     1, 0.5, xRow[ 7], 1.0, xRow[ 6]),
    PadInfo("Coincidences",  4, 0.0, xRow[ 8], 1.0, xRow[ 7]),
    PadInfo("Triggers",      1, 0.0, xRow[ 9], 1.0, xRow[ 8]),
    PadInfo("Decisions",     3, 0.0, xRow[10], 1.0, xRow[ 9]),
    PadInfo("ChargeTrend",   1, 0.0, 0.0,      1.0, xRow[10])
  };

  for (Int_t iPad=0; iPad<nPads; ++iPad) {
    TVirtualPad *pad = new TPad(padInfo[iPad].fName, padInfo[iPad].fName,
				padInfo[iPad].fX1, padInfo[iPad].fY1, padInfo[iPad].fX2, padInfo[iPad].fY2);
    if (padInfo[iPad].fDivX != 1)
      pad->Divide(padInfo[iPad].fDivX, 1);
    pad->Draw();
  }

  return c1;
}

TVirtualPad* AliADQAChecker::GetPadByName(TCanvas *c1, const char* name) const {
  if (!c1)
    AliFatal("c1==NULL");
  TVirtualPad *p = dynamic_cast<TVirtualPad*>(c1->FindObject(name));
  if (!p)
    AliFatalF("pad '%s' does not exist", name);
  return p;
}

TSpline3* AliADQAChecker::MakeTimeSlewingSpline(TH2* h) const {
  TH1* hMeanTimeVsCharge = new TH1F("hMeanTimeVsCharge", ";-log_{10}(charge)", 200, -4, 0);
  for (Int_t j=0; j<200; ++j) {
    TH1* hTimeSlice = h->ProjectionY("hTimeSlice", j+1, j+1);
    if (hTimeSlice->GetEntries()<100 && j>100) {
      hMeanTimeVsCharge->SetBinContent(j+1, hMeanTimeVsCharge->GetBinContent(j));
      hMeanTimeVsCharge->SetBinError  (j+1, hMeanTimeVsCharge->GetBinError(j));
    } else {
      hMeanTimeVsCharge->SetBinContent(j+1, hTimeSlice->GetMean());
      if (hTimeSlice->GetEntries())
	hMeanTimeVsCharge->SetBinError(j+1, 1.0/TMath::Power(hTimeSlice->GetEntries(), 2));
    }
    delete hTimeSlice;
  } // next charge bin

  TSpline3 *s = new TSpline3(hMeanTimeVsCharge);
  delete hMeanTimeVsCharge;
  s->SetBit(TObject::kCanDelete);
  s->SetLineWidth(3);
  s->SetLineColor(kPink);
  return s;
}

//____________________________________________________________________________
void AliADQAChecker::MakeImage(TObjArray **list, AliQAv1::TASKINDEX_t task, AliQAv1::MODE_t mode)
{
  //AliInfoF("task=%d mode=%d", task, mode);
  // makes the QA image for sim and rec
  for (Int_t esIndex=0; esIndex<AliRecoParam::kNSpecies; ++esIndex) {
    if (!AliQAv1::Instance(AliQAv1::GetDetIndex(GetName()))->IsEventSpecieSet(AliRecoParam::ConvertIndex(esIndex)) ||
	!list ||
	!list[esIndex] ||
	!list[esIndex]->GetEntries())
      continue;

    TObjArray *listEs = list[esIndex];

    //  SetStats(kFALSE) for all histograms
    for (Int_t i=0, n=listEs->GetEntries(); i<n; ++i) {
      TH1 *h = dynamic_cast<TH1*>(listEs->At(i));
      if (h)
	h->SetStats(kFALSE);
    }

    // make the canvas if it does not exist
    if (!fImage[esIndex]) {
      const TString title = TString::Format("QA_%s_%s_%s", GetName(),
					    AliQAv1::GetTaskName(task).Data(),
					    AliRecoParam::GetEventSpecieName(esIndex));
      fImage[esIndex] = new TCanvas(title, title, 2500, 5500);
    }
    // clear all objects associated with the canvas
    TCanvas *c1 = fImage[esIndex];
    c1->Clear();
    // and (re)create all pads
    CreatePads(c1)->cd();

    // draw top label
    TLatex topText;
    topText.SetTextAlign(23);
    topText.SetTextSize(.038);
    topText.SetTextFont(42);
    topText.SetTextColor(kBlue+3);
    topText.DrawLatexNDC(0.5, 0.99,
			 Form("%s, %s, Run: %d",
			      AliQAv1::GetTaskName(task).Data(),
			      AliRecoParam::GetEventSpecieName(esIndex),
			      AliQAChecker::Instance()->GetRunNumber()));


    Float_t nEvents = 0;

    // (1) Charge pad
    GetPadByName(c1, "Charge")->cd(1)->SetLogy();
    TH1* hChargeADA = dynamic_cast<TH1*>(listEs->At(AliADQADataMakerRec::kChargeADA));
    TH1* hChargeADC = dynamic_cast<TH1*>(listEs->At(AliADQADataMakerRec::kChargeADC));
    if (hChargeADA && hChargeADC) {
      nEvents = hChargeADA->Integral(-1,-1);
      const Float_t maxY = TMath::Max(hChargeADA->GetBinContent(hChargeADA->GetMaximumBin()),
				      hChargeADC->GetBinContent(hChargeADC->GetMaximumBin()));
      const Float_t minY = TMath::Min(hChargeADA->GetBinContent(hChargeADA->GetMinimumBin()),
				      hChargeADC->GetBinContent(hChargeADC->GetMinimumBin()));
      hChargeADA->DrawCopy()->GetYaxis()->SetRangeUser(minY+1 ,2*maxY);
      hChargeADC->Draw("SAME");
      TLegend *leg = MakeLegend(0.70,0.67,0.97,0.82, 0.05);
      leg->AddEntry(hChargeADA, "ADA", "L");
      leg->AddEntry(hChargeADC, "ADC", "L");
      leg->Draw();
    }

    GetPadByName(c1, "Charge")->cd(2)->SetLogz();
    TH1* hChargeBB       = dynamic_cast<TH1*>(listEs->At(AliADQADataMakerRec::kChargeEoIBB));
    TH1 *hChargeBB_zoomY = NULL;
    if (hChargeBB)
      hChargeBB->Draw("COLZ");

    // (1) Charge pad - zoomed
    GetPadByName(c1, "ChargeZoom")->cd(5)->SetLogz();
    if (hChargeBB) {
      hChargeBB_zoomY = hChargeBB->DrawCopy("COLZ");
      hChargeBB_zoomY->GetYaxis()->SetRangeUser(fQAParam->GetChargeChannelZoomMin(), fQAParam->GetChargeChannelZoomMax());
    }
    GetPadByName(c1, "ChargeZoom")->cd(4)->SetLogz();
    TH1 *hChargeAll       = dynamic_cast<TH1*>(listEs->At(AliADQADataMakerRec::kChargeEoI));
    TH1 *hChargeAll_zoomY = NULL;
    if (hChargeAll) {
      hChargeAll_zoomY = hChargeAll->DrawCopy("COLZ");
      hChargeAll_zoomY->GetYaxis()->SetRangeUser(fQAParam->GetChargeChannelZoomMin(), fQAParam->GetChargeChannelZoomMax());
    }
    GetPadByName(c1, "ChargeZoom")->cd(6)->SetLogz();
    if (hChargeBB_zoomY && hChargeAll_zoomY) {
      TH1* hRatioBBToAll = MakeRatioHistogram(hChargeBB_zoomY, hChargeAll_zoomY);
      hRatioBBToAll->SetTitle("Ratio: w_BB_Flag/All");
      hRatioBBToAll->Draw("COLZ"); // will be deleted in c1->Clear()
    }
    for (Int_t iHist=0; iHist<3; ++iHist) {
      GetPadByName(c1, "ChargeZoom")->cd(iHist+1)->SetLogz(iHist == 2);
      TH1 *h =dynamic_cast<TH1*>(listEs->At(AliADQADataMakerRec::kChargeVsClockInt0+iHist));
      if (h)
	h->Draw("COLZ");
    }

    // (3) Time pad
    for (Int_t iHist=0; iHist<4; ++iHist) {
      GetPadByName(c1, "Time")->cd(iHist+1)->SetLogz(iHist == 3);
      TH1 *h = dynamic_cast<TH1*>(listEs->At(AliADQADataMakerRec::kHPTDCTime+iHist));
      if (h)
	h->Draw("COLZ");
    }

    // (4) Time ratio pad
    GetPadByName(c1, "TimeRatio")->cd(1);
    TH1* hFlagNoTime = dynamic_cast<TH1*>(listEs->At(AliADQADataMakerRec::kFlagNoTime));
    TH1* hTimeNoFlag = dynamic_cast<TH1*>(listEs->At(AliADQADataMakerRec::kTimeNoFlag));
    if (hFlagNoTime && hTimeNoFlag) {
      TH1* hFlagNoTimeScaled = MakeScaledHistogram(hFlagNoTime, nEvents);
      TH1* hTimeNoFlagScaled = MakeScaledHistogram(hTimeNoFlag, nEvents);
      hFlagNoTimeScaled->GetYaxis()->SetRangeUser(0, 2*TMath::Max(hFlagNoTimeScaled->GetBinContent(hFlagNoTimeScaled->GetMaximumBin()),
                                                                  hTimeNoFlagScaled->GetBinContent(hTimeNoFlag->GetMaximumBin())));
      hFlagNoTimeScaled->Draw();
      hTimeNoFlagScaled->Draw("SAME");
      TLegend *leg = MakeLegend(0.15,0.78,0.85,0.88, 0.04);
      leg->AddEntry(hFlagNoTimeScaled, "Events with BB/BG flag but no time", "L");
      leg->AddEntry(hTimeNoFlagScaled, "Events with time but no BB/BG flag", "L");
      leg->Draw();
    }
    TH1* hTimeAll = dynamic_cast<TH1*>(listEs->At(AliADQADataMakerRec::kHPTDCTimeRebin));
    if (hTimeAll) {
      for (Int_t iHist=1; iHist<3; ++iHist) {
	GetPadByName(c1, "TimeRatio")->cd(1+iHist)->SetLogz();
	TH1* hTimeSel = dynamic_cast<TH1*>(listEs->At(AliADQADataMakerRec::kHPTDCTimeRebin + iHist));
	if (!hTimeSel)
	  continue;
	TH1* hTimeRatio = MakeRatioHistogram(hTimeSel, hTimeAll);
	if (iHist == 1)
	  hTimeRatio->GetYaxis()->SetRangeUser(fQAParam->GetTdcTimeMinBBFlag(), fQAParam->GetTdcTimeMaxBBFlag());
	if (iHist == 2)
	  hTimeRatio->GetYaxis()->SetRangeUser(fQAParam->GetTdcTimeMinBGFlag(), fQAParam->GetTdcTimeMaxBGFlag());
	hTimeRatio->Draw("COLZ");
      }
    }
    GetPadByName(c1, "TimeRatio")->cd(4);
    TH1 *hTimeBB = dynamic_cast<TH1*>(listEs->At(AliADQADataMakerRec::kNEventsBBFlag));
    TH1 *hTimeBG = dynamic_cast<TH1*>(listEs->At(AliADQADataMakerRec::kNEventsBGFlag));
    if (hTimeBB && hTimeBG) {
      TH1* hTimeBBScaled = MakeScaledHistogram(hTimeBB, nEvents, kTRUE);
      TH1* hTimeBGScaled = MakeScaledHistogram(hTimeBG, nEvents, kTRUE);
      hTimeBBScaled->GetYaxis()->SetRangeUser(0.8*TMath::Min(hTimeBBScaled->GetBinContent(hTimeBBScaled->GetMinimumBin()),
							     hTimeBGScaled->GetBinContent(hTimeBGScaled->GetMinimumBin())),
					      1.8*TMath::Max(hTimeBBScaled->GetBinContent(hTimeBBScaled->GetMaximumBin()),
							     hTimeBGScaled->GetBinContent(hTimeBGScaled->GetMaximumBin())));
      hTimeBBScaled->Draw("E");
      // hTimeBBScaled->Draw("ESAME");
      // hTimeBGScaled->DrawCopy("SAMEHIST");
      hTimeBGScaled->Draw("ESAME");
      TLegend *leg = MakeLegend(0.15,0.78,0.85,0.88, 0.04);
      leg->AddEntry(hTimeBBScaled, "Events with BB flag", "L");
      leg->AddEntry(hTimeBGScaled, "Events with BG flag", "L");
      leg->Draw();
    }

    // (5) Clock Flags pad
    for (Int_t iHist=0; iHist<5; ++iHist) {
      GetPadByName(c1, "ClockCfg")->cd(iHist+1)->SetLogz(iHist > 1);
      TH2* hFlagVsClock = dynamic_cast<TH2*>(listEs->At(AliADQADataMakerRec::kBBFlagVsClock+iHist));
      if (hFlagVsClock)
	hFlagVsClock->Draw("COLZ");
    }

    // (6) Coincidences pad
    for (Int_t iHist=0; iHist<2; ++iHist) {
      GetPadByName(c1, "Coincidences")->cd(iHist+1)->SetLogy();
      TH1* hCoincADC = dynamic_cast<TH1*>(listEs->At(AliADQADataMakerRec::kNBBCoincADC+2*iHist));
      TH1* hCoincADA = dynamic_cast<TH1*>(listEs->At(AliADQADataMakerRec::kNBBCoincADA+2*iHist));
      if (hCoincADC && hCoincADA) {
	hCoincADC->Draw();
	hCoincADA->Draw("SAME");
	TLegend *leg = MakeLegend(0.70,0.67,0.97,0.82, 0.05);
	leg->AddEntry(hCoincADA, "ADA", "L");
	leg->AddEntry(hCoincADC, "ADC", "L");
	leg->Draw();
      }
    }
    for (Int_t iHist=0; iHist<2; ++iHist) {
      GetPadByName(c1, "Coincidences")->cd(iHist+3)->SetLogz();
      TH1 *hBBCoinc = dynamic_cast<TH1*>(listEs->At(AliADQADataMakerRec::kNBBCoincCorr+iHist));
      if (hBBCoinc)
	hBBCoinc->Draw("COLZ");
    }

    // (7) Pedestal monitoring pad
    for (Int_t iHist=0; iHist<2; ++iHist) {
      GetPadByName(c1, "Pedestal")->cd(iHist+1);
      TH1* hPed = (AliRecoParam::ConvertIndex(esIndex) == AliRecoParam::kCalib
		   ? dynamic_cast<TH1*>(listEs->At(AliADQADataMakerRec::kPedestalInt0+iHist))
		   : dynamic_cast<TH1*>(listEs->At(AliADQADataMakerRec::kPedestalDiffInt0+iHist)));
      if (hPed)
	hPed->Draw("COLZ");
    }

    // (8) Saturation monitoring pad
    GetPadByName(c1, "MaxCharge")->cd()->SetLogz();
    TH2* hChargeSaturation = dynamic_cast<TH2*>(listEs->At(AliADQADataMakerRec::kChargeSaturation));
    if (hChargeSaturation) {
      TH2* hChargeSaturationRebinned = hChargeSaturation->RebinY(64, Form("%s_rebinY64", hChargeSaturation->GetName()));
      hChargeSaturationRebinned->SetBit(TObject::kCanDelete);
      hChargeSaturationRebinned->Draw("COLZ");
    }

    // (9) Trigger inputs pad
    GetPadByName(c1, "Triggers")->cd()->SetLogy();
    TH1 *hTrig = dynamic_cast<TH1*>(listEs->At(AliADQADataMakerRec::kTriggers));
    if (hTrig) {
      hTrig->Draw("HIST");
      hTrig->Draw("TEXT0SAME");
    }

    // (10) Mean time pad
    GetPadByName(c1, "MeanTime")->cd(1);
    TH1 *hMeanTimeADA = dynamic_cast<TH1*>(listEs->At(AliADQADataMakerRec::kMeanTimeADA));
    TH1 *hMeanTimeADC = dynamic_cast<TH1*>(listEs->At(AliADQADataMakerRec::kMeanTimeADC));
    if (hMeanTimeADA && hMeanTimeADC) {
      hMeanTimeADA->DrawCopy()->GetYaxis()->SetRangeUser(0, 1.1*TMath::Max(hMeanTimeADA->GetBinContent(hMeanTimeADA->GetMaximumBin()),
									   hMeanTimeADC->GetBinContent(hMeanTimeADC->GetMaximumBin())));
      hMeanTimeADC->Draw("SAME");
      TLegend *leg = MakeLegend(0.70,0.67,0.97,0.82, 0.05);
      leg->AddEntry(hMeanTimeADA, "ADA", "L");
      leg->AddEntry(hMeanTimeADC, "ADC", "L");
      leg->Draw();
    }
    for (Int_t iHist=0; iHist<3; ++iHist) {
      GetPadByName(c1, "MeanTime")->cd(iHist+2)->SetLogz();
      TH1* hMeanTimeDiff = dynamic_cast<TH1*>(listEs->At(AliADQADataMakerRec::kMeanTimeDiff+iHist));
      if (hMeanTimeDiff)
	hMeanTimeDiff->Draw("COLZ");
    }

    // (11) time slewing and Decisions pad
    TH2 *hTimeSlewing[2] = {
      dynamic_cast<TH2*>(listEs->At(AliADQADataMakerRec::kTimeSlewingADA)),
      dynamic_cast<TH2*>(listEs->At(AliADQADataMakerRec::kTimeSlewingADC))
    };
    const char* sideName[2] = { "ADA", "ADC" };
    for (Int_t side=0; side<2; ++side) {
      GetPadByName(c1, "Decisions")->cd(1+2*side)->SetLogz();
      if (!hTimeSlewing[side])
	continue;
      hTimeSlewing[side] = dynamic_cast<TH2*>(hTimeSlewing[side]->DrawCopy("COLZ"));
      hTimeSlewing[side]->GetYaxis()->SetRangeUser(fQAParam->GetTdcTimeMinBBFlag(), fQAParam->GetTdcTimeMaxBBFlag());
      hTimeSlewing[side]->SetTitle(TString::Format("Time slewing %s", sideName[side]));
      MakeTimeSlewingSpline(hTimeSlewing[side])->Draw("PSAME");
    }
    GetPadByName(c1, "Decisions")->cd(2);
    TH1* hDecisions = dynamic_cast<TH1*>(listEs->At(AliADQADataMakerRec::kDecisions));
    if (hDecisions)
      hDecisions->Draw("COLZTEXT");

    // (12) charge trend
    GetPadByName(c1, "ChargeTrend")->cd();
    TH1* hChargeQuantileADA = dynamic_cast<TH1*>(listEs->At(AliADQADataMakerRec::kTrend_TriggerChargeQuantileADA));
    TH1* hChargeQuantileADC = dynamic_cast<TH1*>(listEs->At(AliADQADataMakerRec::kTrend_TriggerChargeQuantileADC));
    if (hChargeQuantileADA && hChargeQuantileADC) {
      hChargeQuantileADA = hChargeQuantileADA->DrawCopy();
      hChargeQuantileADA->GetYaxis()->SetRangeUser(fQAParam->GetChargeTrendMin(), fQAParam->GetChargeTrendMax());
      hChargeQuantileADA->GetYaxis()->SetTitle("Quantile 0.9");
      hChargeQuantileADA->GetXaxis()->SetRange(1, hChargeQuantileADC->GetNbinsX()-1);
      hChargeQuantileADC->Draw("SAME");
      TLegend *leg = MakeLegend(0.70,0.67,0.97,0.82, 0.05);
      leg->AddEntry(hChargeQuantileADA, "ADA", "L");
      leg->AddEntry(hChargeQuantileADC, "ADC", "L");
      leg->Draw();
    }

    c1->SaveAs(Form("QAsummary_%d_%d.png",
		    AliQAChecker::Instance()->GetRunNumber(),
		    esIndex));
    //c1->Print(Form("%s%s%d.%s", AliQAv1::GetImageFileName(), AliQAv1::GetModeName(mode), AliQAChecker::Instance()->GetRunNumber(), AliQAv1::GetImageFileFormat()), "ps");
  } // next species
}
