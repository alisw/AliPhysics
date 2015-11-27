/**************************************************************************
 * Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
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
/**
 * @file AliEmcalTriggerQAAP.cxx
 * @date Nov. 12, 2015
 * @author Salvatore Aiola <salvatore.aiola@cern.ch>, Yale University
 */

#include <THashList.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TProfile.h>
#include <TObjString.h>
#include <TObjArray.h>

#include "AliEmcalTriggerPatchInfoAPV1.h"
#include "AliEmcalTriggerFastORAP.h"
#include "AliLog.h"

#include "AliEmcalTriggerConstantsAP.h"

#include "AliEmcalTriggerQAAP.h"

using namespace EmcalTriggerAP;

/// \cond CLASSIMP
ClassImp(AliEmcalTriggerQAAP)
/// \endcond

const Int_t AliEmcalTriggerQAAP::fgkMaxPatchAmp[6] = {2000, 1500, 1500, 2000, 1500, 2000};
const TString AliEmcalTriggerQAAP::fgkPatchTypes[3] = {"Online", "Recalc", "Offline"};
/**
 * Dummy constructor
 */
AliEmcalTriggerQAAP::AliEmcalTriggerQAAP():
  TNamed(),
  fFastorL0Th(1),
  fFastorL1Th(1),
  fBkgPatchType(kTMEMCalBkg),
  fDebugLevel(0),
  fHistos(0)
{
  for (Int_t i = 0; i < 3; i++) {
    fBkgADCAmpEMCal[i].Set(100);
    fBkgADCAmpDCal[i].Set(100);
    fNBkgPatchesEMCal[i] = 0;
    fNBkgPatchesDCal[i] = 0;
    fEnabledPatchTypes[i] = kFALSE;

    for (Int_t itype = 0; itype < 6; itype++) {
      fMaxPatchEMCal[itype][i] = 0;
      fMaxPatchDCal[itype][i] = 0;
    }
  }
}

/**
 * Constructor
 */
AliEmcalTriggerQAAP::AliEmcalTriggerQAAP(const char* name):
  TNamed(name,name),
  fFastorL0Th(1),
  fFastorL1Th(1),
  fBkgPatchType(kTMEMCalBkg),
  fDebugLevel(0),
  fHistos(0)
{
  for (Int_t i = 0; i < 3; i++) {
    fBkgADCAmpEMCal[i].Set(100);
    fBkgADCAmpDCal[i].Set(100);
    fNBkgPatchesEMCal[i] = 0;
    fNBkgPatchesDCal[i] = 0;
    fEnabledPatchTypes[i] = kFALSE;

    for (Int_t itype = 0; itype < 6; itype++) {
      fMaxPatchEMCal[itype][i] = 0;
      fMaxPatchDCal[itype][i] = 0;
    }
  }
}

/**
 * Destructor
 */
AliEmcalTriggerQAAP::~AliEmcalTriggerQAAP()
{
}

/**
 * Set the patch types to be plotted
 */
void AliEmcalTriggerQAAP::EnablePatchType(PatchTypes_t type, Bool_t e)
{
  fEnabledPatchTypes[type] = e;
}

/**
 * Initialize the class, i.e. allocate histograms.
 */
void AliEmcalTriggerQAAP::Init()
{
  TString hname;
  TString htitle;

  fHistos = new THashList();
  fHistos->SetName(Form("histos%s", GetName()));
  fHistos->SetOwner(kTRUE);

  hname = Form("EMCTRQA_histFastORL0");
  htitle = Form("EMCTRQA_histFastORL0;FastOR abs. ID;entries above L0 threshold");
  CreateTH1(hname, htitle, 5000, 0, 5000);

  hname = Form("EMCTRQA_histFastORL0Mean");
  htitle = Form("EMCTRQA_histFastORL0Mean;FastOR abs. ID;mean ADC counts");
  CreateTProfile(hname, htitle, 5000, 0, 5000);

  hname = Form("EMCTRQA_histLargeAmpFastORL0");
  htitle = Form("EMCTRQA_histLargeAmpFastORL0 (>400);FastOR abs. ID;entries above 400");
  CreateTH1(hname, htitle, 5000, 0, 5000);

  hname = Form("EMCTRQA_histFastORL1");
  htitle = Form("EMCTRQA_histFastORL1;FastOR abs. ID;entries above L1 threshold");
  CreateTH1(hname, htitle, 5000, 0, 5000);

  hname = Form("EMCTRQA_histFastORL1Mean");
  htitle = Form("EMCTRQA_histFastORL1Mean;FastOR abs. ID;mean L1 time sum");
  CreateTProfile(hname, htitle, 5000, 0, 5000);

  hname = Form("EMCTRQA_histLargeAmpFastORL1");
  htitle = Form("EMCTRQA_histLargeAmpFastORL1 (>400);FastOR abs. ID;entries above 400");
  CreateTH1(hname, htitle, 5000, 0, 5000);

  hname = Form("EMCTRQA_histFastORL1AmpVsL0Amp");
  htitle = Form("EMCTRQA_histFastORL1AmpVsL0Amp;L0 amplitude;L1 time sum;entries");
  CreateTH2(hname, htitle, 128, 0, 2048, 128, 0, 2048);

  for (Int_t itype = 0; itype < 3; itype++) {
    if (!fEnabledPatchTypes[itype]) continue;
    hname = Form("EMCTRQA_histEMCalMedianVsDCalMedian%s", fgkPatchTypes[itype].Data());
    htitle = Form("EMCTRQA_histEMCalMedianVsDCalMedian%s;EMCal median;DCal median;entries", fgkPatchTypes[itype].Data());
    CreateTH2(hname, htitle, fgkMaxPatchAmp[fBkgPatchType]/10, 0, fgkMaxPatchAmp[fBkgPatchType], fgkMaxPatchAmp[fBkgPatchType]/10, 0, fgkMaxPatchAmp[fBkgPatchType]);
  }

  for (Int_t itrig = 0; itrig < 6; itrig++) {
    if (kEMCalTriggerNames[itrig].IsNull()) continue;
    for (Int_t itype = 0; itype < 3; itype++) {
      if (!fEnabledPatchTypes[itype]) continue;
      hname = Form("EMCTRQA_histEMCalPatchAmp%s%s", kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
      htitle = Form("EMCTRQA_histEMCalPatchAmp%s%s;amplitude;entries", kEMCalTriggerNames[itype].Data(), fgkPatchTypes[itype].Data());
      CreateTH1(hname, htitle, fgkMaxPatchAmp[itrig]/10, 0, fgkMaxPatchAmp[itrig]);

      hname = Form("EMCTRQA_histDCalPatchAmp%s%s", kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
      htitle = Form("EMCTRQA_histDCalPatchAmp%s%s;amplitude;entries", kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
      CreateTH1(hname, htitle, fgkMaxPatchAmp[itrig]/10, 0, fgkMaxPatchAmp[itrig]);

      hname = Form("EMCTRQA_histEMCalMedianVsDCalMax%s%s", kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
      htitle = Form("EMCTRQA_histEMCalMedianVsDCalMax%s%s;DCal max;EMCal median;entries", kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
      CreateTH2(hname, htitle, fgkMaxPatchAmp[itrig]/10, 0, fgkMaxPatchAmp[itrig], fgkMaxPatchAmp[itrig]/10, 0, fgkMaxPatchAmp[itrig]);

      hname = Form("EMCTRQA_histDCalMedianVsEMCalMax%s%s", kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
      htitle = Form("EMCTRQA_histDCalMedianVsEMCalMax%s%s;EMCal max;DCal median;entries", kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
      CreateTH2(hname, htitle, fgkMaxPatchAmp[itrig]/10, 0, fgkMaxPatchAmp[itrig], fgkMaxPatchAmp[itrig]/10, 0, fgkMaxPatchAmp[itrig]);

      hname = Form("EMCTRQA_histEMCalMedianVsEMCalMax%s%s", kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
      htitle = Form("EMCTRQA_histEMCalMedianVsEMCalMax%s%s;EMCal max;EMCal median;entries", kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
      CreateTH2(hname, htitle, fgkMaxPatchAmp[itrig]/10, 0, fgkMaxPatchAmp[itrig], fgkMaxPatchAmp[itrig]/10, 0, fgkMaxPatchAmp[itrig]);

      hname = Form("EMCTRQA_histDCalMedianVsDCalMax%s%s", kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
      htitle = Form("EMCTRQA_histDCalMedianVsDCalMax%s%s;DCal max;DCal median;entries", kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
      CreateTH2(hname, htitle, fgkMaxPatchAmp[itrig]/10, 0, fgkMaxPatchAmp[itrig], fgkMaxPatchAmp[itrig]/10, 0, fgkMaxPatchAmp[itrig]);

      hname = Form("EMCTRQA_histEMCalMaxVsDCalMax%s%s", kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
      htitle = Form("EMCTRQA_histEMCalMaxVsDCalMax%s%s;EMCal max;DCal max;entries", kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
      CreateTH2(hname, htitle, fgkMaxPatchAmp[itrig]/10, 0, fgkMaxPatchAmp[itrig], fgkMaxPatchAmp[itrig]/10, 0, fgkMaxPatchAmp[itrig]);

      hname = Form("EMCTRQA_histEdgePos%s%s", kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
      htitle = Form("EMCTRQA_histEdgePos%s%s;col;row;entries", kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
      CreateTH2(hname, htitle, 48, 0, 48, 105, 0, 105);

      hname = Form("EMCTRQA_histCMPos%s%s", kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
      htitle = Form("EMCTRQA_histCMPos%s%s;#eta;#phi;entries", kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
      CreateTH2(hname, htitle, 60, -1, 1, 200, 0, TMath::TwoPi());

      hname = Form("EMCTRQA_histGeoPos%s%s", kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
      htitle = Form("EMCTRQA_histGeoPos%s%s;#eta;#phi;entries", kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
      CreateTH2(hname, htitle, 60, -1, 1, 200, 0, TMath::TwoPi());

      hname = Form("EMCTRQA_histLargeAmpEdgePos%s%s", kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
      htitle = Form("EMCTRQA_histLargeAmpEdgePos%s%s (>700);col;row;entries above 700", kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
      CreateTH2(hname, htitle, 48, 0, 48, 105, 0, 105);

      hname = Form("EMCTRQA_histEMCalPatchEnergy%s%s", kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
      htitle = Form("EMCTRQA_histEMCalPatchEnergy%s%s;energy (GeV);entries", kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
      CreateTH1(hname, htitle, 200, 0, 200);

      hname = Form("EMCTRQA_histDCalPatchEnergy%s%s", kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
      htitle = Form("EMCTRQA_histDCalPatchEnergy%s%s;energy (GeV);entries", kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
      CreateTH1(hname, htitle, 200, 0, 200);
    }
  }
}

/**
 * Process a patch, filling relevant histograms.
 * \param patch Pointer to a valid trigger patch
 */
void AliEmcalTriggerQAAP::ProcessPatch(AliEmcalTriggerPatchInfoAPV1* patch)
{
  TString hname;

  Int_t triggerBits[6] = { patch->GetTriggerBitConfig()->GetLevel0Bit(),
      patch->GetTriggerBitConfig()->GetGammaLowBit(),
      patch->GetTriggerBitConfig()->GetGammaHighBit(),
      patch->GetTriggerBitConfig()->GetJetLowBit(),
      patch->GetTriggerBitConfig()->GetJetHighBit(),
      patch->GetTriggerBitConfig()->GetBkgBit()
  };

  Int_t offsets[3] = { 0, AliEmcalTriggerPatchInfoAPV1::kRecalcOffset, AliEmcalTriggerPatchInfoAPV1::kOfflineOffset };
  Int_t amplitudes[3] = { patch->GetADCAmp(),  patch->GetADCAmp(),  patch->GetADCOfflineAmp() };

  for (Int_t itrig = 0; itrig < 6; itrig++) {
    if (kEMCalTriggerNames[itrig].IsNull()) continue;
    for (Int_t itype = 0; itype < 3; itype++) {
      if (!fEnabledPatchTypes[itype]) continue;
      if (!patch->TestTriggerBit(triggerBits[itrig]+offsets[itype])) continue;

      hname = Form("EMCTRQA_histEdgePos%s%s", kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
      FillTH2(hname, patch->GetColStart(), patch->GetRowStart());

      hname = Form("EMCTRQA_histCMPos%s%s", kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
      FillTH2(hname, patch->GetEtaCM(), patch->GetPhiCM());

      hname = Form("EMCTRQA_histGeoPos%s%s", kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
      FillTH2(hname, patch->GetEtaGeo(), patch->GetPhiGeo());

      if (amplitudes[itype] > 700) {
        hname = Form("EMCTRQA_histLargeAmpEdgePos%s%s", kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
        FillTH2(hname, patch->GetColStart(), patch->GetRowStart());
      }

      TString det;

      if (patch->IsEMCal()) {
        det = "EMCal";
        if (fMaxPatchEMCal[itrig][itype] < amplitudes[itype]) fMaxPatchEMCal[itrig][itype] = amplitudes[itype];

        if (itrig == fBkgPatchType) {
          if (fNBkgPatchesEMCal[itype] >= fBkgADCAmpEMCal[itype].GetSize()) {
            fBkgADCAmpEMCal[itype].Set((fNBkgPatchesEMCal[itype]+1)*2);
          }
          fBkgADCAmpEMCal[itype].AddAt(amplitudes[itype], fNBkgPatchesEMCal[itype]);
          fNBkgPatchesEMCal[itype]++;
        }
      }
      else if (patch->IsDCalPHOS()) {
        det = "DCal";
        if (fMaxPatchDCal[itrig][itype] < amplitudes[itype]) fMaxPatchDCal[itrig][itype] = amplitudes[itype];

        if (itrig == fBkgPatchType) {
          if (fNBkgPatchesDCal[itype] >= fBkgADCAmpDCal[itype].GetSize()) {
            fBkgADCAmpDCal[itype].Set((fNBkgPatchesDCal[itype]+1)*2);
          }
          fBkgADCAmpDCal[itype].AddAt(amplitudes[itype], fNBkgPatchesDCal[itype]);
          fNBkgPatchesDCal[itype]++;
        }
      }
      else {
        AliWarning(Form("Patch is not EMCal nor DCal/PHOS (pos: %d, %d)", patch->GetRowStart(), patch->GetColStart()));
      }

      hname = Form("EMCTRQA_hist%sPatchAmp%s%s", det.Data(), kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
      FillTH1(hname, amplitudes[itype]);

      hname = Form("EMCTRQA_hist%sPatchEnergy%s%s", det.Data(), kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
      FillTH1(hname, patch->GetPatchE());

      if (fDebugLevel >= 2) {
        Printf("Type = %s; global pos = (%d, %d); Amp (online) = %d; Amp (offline) = %d; Patch energy = %.3f\n"
            "Position (CM): Eta=%.3f, Phi=%.3f\n"
            "Position (Geo): Eta=%.3f, Phi=%.3f\n",
            kEMCalTriggerNames[itrig].Data(), patch->GetRowStart(), patch->GetColStart(), patch->GetADCAmp(), patch->GetADCOfflineAmp(), patch->GetPatchE(),
            patch->GetEtaCM(), patch->GetPhiCM(),
            patch->GetEtaGeo(), patch->GetPhiGeo());
      }
    }
  }
}

/**
 * Process a FastOR, filling relevant histograms.
 * \param patch Pointer to a valid trigger FastOR
 */
void AliEmcalTriggerQAAP::ProcessFastor(AliEmcalTriggerFastORAP* fastor)
{
  TString hname;

  if (fastor->GetL0Amp() > fFastorL0Th) {
    hname = Form("EMCTRQA_histFastORL0");
    FillTH1(hname, fastor->GetAbsId());

    hname = Form("EMCTRQA_histFastORL0Mean");
    FillTProfile(hname, fastor->GetAbsId(), fastor->GetL0Amp());

    if (fastor->GetL0Amp() > 400) {
      hname = Form("EMCTRQA_histLargeAmpFastORL0");
      FillTH1(hname, fastor->GetAbsId());
    }
  }

  if (fastor->GetL1Amp() > fFastorL1Th) {
    hname = Form("EMCTRQA_histFastORL1");
    FillTH1(hname, fastor->GetAbsId());

    hname = Form("EMCTRQA_histFastORL1Mean");
    FillTProfile(hname, fastor->GetAbsId(), fastor->GetL1Amp());

    if (fastor->GetL1Amp() > 400) {
      hname = Form("EMCTRQA_histLargeAmpFastORL1");
      FillTH1(hname, fastor->GetAbsId());
    }
  }

  if (fastor->GetL1Amp() > fFastorL1Th && fastor->GetL0Amp() > fFastorL0Th) {
    hname = Form("EMCTRQA_histFastORL1AmpVsL0Amp");
    FillTH2(hname, fastor->GetL0Amp(), fastor->GetL1Amp());
  }
}

/**
 * This method should be called at the end of each event.
 */
void AliEmcalTriggerQAAP::EventCompleted()
{
  AliDebug(2, Form("Entering AliEmcalTriggerQAAP::EventCompleted"));

  TString hname;

  Double_t medianEMCal[3] = {0};
  Double_t medianDCal[3] = {0};

  for (Int_t itype = 0; itype < 3; itype++) {
    if (!fEnabledPatchTypes[itype]) continue;

    AliDebug(2, Form("Patch type %s", fgkPatchTypes[itype].Data()));

    medianEMCal[itype] = TMath::Median(fNBkgPatchesEMCal[itype], fBkgADCAmpEMCal[itype].GetArray());
    medianDCal[itype] = TMath::Median(fNBkgPatchesDCal[itype], fBkgADCAmpDCal[itype].GetArray());

    hname = Form("EMCTRQA_histEMCalMedianVsDCalMedian%s", fgkPatchTypes[itype].Data());
    FillTH2(hname, medianEMCal[itype], medianDCal[itype]);

    for (Int_t itrig = 0; itrig < 6; itrig++) {
      if (kEMCalTriggerNames[itrig].IsNull()) continue;
      if (fMaxPatchDCal[itrig][itype] == 0 && fMaxPatchEMCal[itrig][itype] == 0) continue;

      AliDebug(2, Form("Trigger type: %s", kEMCalTriggerNames[itype].Data()));

      hname = Form("EMCTRQA_histEMCalMedianVsDCalMax%s%s", kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
      FillTH2(hname, fMaxPatchDCal[itrig][itype], medianEMCal[itype]);

      hname = Form("EMCTRQA_histDCalMedianVsEMCalMax%s%s", kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
      FillTH2(hname, fMaxPatchEMCal[itrig][itype], medianDCal[itype]);

      hname = Form("EMCTRQA_histDCalMedianVsDCalMax%s%s", kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
      FillTH2(hname, fMaxPatchDCal[itrig][itype], medianDCal[itype]);

      hname = Form("EMCTRQA_histEMCalMedianVsEMCalMax%s%s", kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
      FillTH2(hname, fMaxPatchEMCal[itrig][itype], medianEMCal[itype]);

      hname = Form("EMCTRQA_histEMCalMaxVsDCalMax%s%s", kEMCalTriggerNames[itype].Data(), fgkPatchTypes[itype].Data());
      FillTH2(hname, fMaxPatchEMCal[itrig][itype], fMaxPatchDCal[itrig][itype]);

      fMaxPatchEMCal[itrig][itype] = 0;
      fMaxPatchDCal[itrig][itype] = 0;
    }
    fNBkgPatchesEMCal[itype] = 0;
    fNBkgPatchesDCal[itype] = 0;
    fBkgADCAmpEMCal[itype].Reset();
    fBkgADCAmpDCal[itype].Reset();
  }
}

//______________________________________________________________________________
void AliEmcalTriggerQAAP::CreateTProfile(const char *name, const char *title, int nbins, double xmin, double xmax)
{
  /*
   * Create a new TProfile within the container.
   *
   * @param name: Name of the histogram
   * @param title: Title of the histogram
   * @param nbins: number of bins
   * @param xmin: min. value of the range
   * @param xmax: max. value of the range
   * Raises fatals in case the object is attempted to be duplicated
   */
  if (fHistos->FindObject(name)) {
    Fatal("AliEmcalTriggerQAAP::CreateTProfile", "Object %s already exists", name);
    return;
  }
  TProfile* hist = new TProfile(name, title, nbins, xmin, xmax);
  fHistos->Add(hist);
}

//______________________________________________________________________________
void AliEmcalTriggerQAAP::CreateTH1(const char *name, const char *title, int nbins, double xmin, double xmax)
{
  /*
   * Create a new TH1 within the container.
   *
   * @param name: Name of the histogram
   * @param title: Title of the histogram
   * @param nbins: number of bins
   * @param xmin: min. value of the range
   * @param xmax: max. value of the range
   * Raises fatals in case the object is attempted to be duplicated
   */
  if (fHistos->FindObject(name)) {
    Fatal("AliEmcalTriggerQAAP::CreateTH1", "Object %s already exists", name);
    return;
  }
  TH1* hist = new TH1D(name, title, nbins, xmin, xmax);
  fHistos->Add(hist);
}

//______________________________________________________________________________
void AliEmcalTriggerQAAP::CreateTH2(const char *name, const char *title, int nbinsx, double xmin, double xmax, int nbinsy, double ymin, double ymax)
{
  /*
   * Create a new TH2 within the container.
   *
   * @param name: Name of the histogram
   * @param title: Title of the histogram
   * @param nbinsx: number of bins in x-direction
   * @param xmin: min. value of the range in x-direction
   * @param xmax: max. value of the range in x-direction
   * @param nbinsy: number of bins in y-direction
   * @param ymin: min. value of the range in y-direction
   * @param ymax: max. value of the range in y-direction
   * Raises fatals in case the object is attempted to be duplicated
   */
  if (fHistos->FindObject(name)) {
    Fatal("AliEmcalTriggerQAAP::CreateTH2", "Object %s already exists", name);
    return;
  }
  TH2* hist = new TH2D(name, title, nbinsx, xmin, xmax, nbinsy, ymin, ymax);
  fHistos->Add(hist);
}

//______________________________________________________________________________
void AliEmcalTriggerQAAP::CreateTH3(const char* name, const char* title, int nbinsx, double xmin, double xmax, int nbinsy, double ymin, double ymax, int nbinsz, double zmin, double zmax)
{
  /*
   * Create a new TH3 within the container.
   *
   * @param nbinsx: number of bins in x-direction
   * @param xmin: min. value of the range in x-direction
   * @param xmax: max. value of the range in x-direction
   * @param nbinsy: number of bins in y-direction
   * @param ymin: min. value of the range in y-direction
   * @param ymax: max. value of the range in y-direction
   * @param nbinsz: number of bins in z-direction
   * @param zmin: min. value of the range in z-direction
   * @param zmax: max. value of the range in z-direction
   * Raises fatals in case the object is attempted to be duplicated
   */
  if (fHistos->FindObject(name)) {
    Fatal("AliEmcalTriggerQAAP::CreateTH3", "Object %s already exists", name);
    return;
  }
  TH3* hist = new TH3D(name, title, nbinsx, xmin, xmax, nbinsy, ymin, ymax, nbinsz, zmin, zmax);
  fHistos->Add(hist);
}

//______________________________________________________________________________
void AliEmcalTriggerQAAP::FillTProfile(const char *name, double x, double y, double weight)
{
  /*
   * Fill a profile histogram within the container.
   *
   * @param name: Name of the histogram
   * @param x: x-coordinate
   * @param weight (@default 1): optional weight of the entry
   * Raises fatals in case the histogram is not found
   */
  TProfile* hist = dynamic_cast<TProfile*>(fHistos->FindObject(name));
  if (!hist) {
    Fatal("AliEmcalTriggerQAAP::FillTProfile", "Histogram %s not found", name);
    return;
  }
  hist->Fill(x, y, weight);
}

//______________________________________________________________________________
void AliEmcalTriggerQAAP::FillTH1(const char *name, double x, double weight)
{
  /*
   * Fill a 1D histogram within the container.
   *
   * @param name: Name of the histogram
   * @param x: x-coordinate
   * @param weight (@default 1): optional weight of the entry
   * Raises fatals in case the histogram is not found
   */
  TH1* hist = dynamic_cast<TH1*>(fHistos->FindObject(name));
  if (!hist) {
    Fatal("AliEmcalTriggerQAAP::FillTH1", "Histogram %s not found", name);
    return;
  }
  hist->Fill(x, weight);
}

//______________________________________________________________________________
void AliEmcalTriggerQAAP::FillTH2(const char *name, double x, double y, double weight)
{
  /*
   * Fill a 2D histogram within the container.
   *
   * @param name: Name of the histogram
   * @param x: x-coordinate
   * @param y: y-coordinate
   * @param weight (@default 1): optional weight of the entry
   * Raises fatals in case the histogram is not found
   */
  TH2* hist = dynamic_cast<TH2*>(fHistos->FindObject(name));
  if (!hist) {
    Fatal("AliEmcalTriggerQAAP::FillTH2", "Histogram %s not found", name);
    return;
  }
  hist->Fill(x, y, weight);
}

//______________________________________________________________________________
void AliEmcalTriggerQAAP::FillTH3(const char* name, double x, double y, double z, double weight)
{
  /*
   * Fill a 3D histogram within the container.
   *
   * @param name: Name of the histogram
   * @param x: x-coordinate
   * @param y: y-coordinate
   * @param z: z-coordinate
   * @param weight (@default 1): optional weight of the entry
   * Raises fatals in case the histogram is not found
   */

  TH3* hist = dynamic_cast<TH3*>(fHistos->FindObject(name));
  if (!hist) {
    Fatal("AliEmcalTriggerQAAP::FillTH3", "Histogram %s not found", name);
    return;
  }
  hist->Fill(x, y, z, weight);
}

//______________________________________________________________________________
TObject *AliEmcalTriggerQAAP::FindObject(const char *name) const
{
  /*
   * Find an object inside the container.
   *
   * @param name: Name of the object to find inside the container
   * @return: pointer to the object (NULL if not found)
   */

  return fHistos->FindObject(name);
}

//______________________________________________________________________________
TObject* AliEmcalTriggerQAAP::FindObject(const TObject* obj) const
{
  /*
   * Find and object inside the container.
   *
   * @param obj: the object to find
   * @return: pointer to the object (NULL if not found)
   */
  TString hname(obj->GetName());
  return fHistos->FindObject(hname);
}
