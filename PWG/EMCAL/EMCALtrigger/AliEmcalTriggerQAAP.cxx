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

#include <cstring>

#include <THashList.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TProfile.h>
#include <TObjString.h>
#include <TObjArray.h>

#include "AliEMCALTriggerPatchInfo.h"
#include "AliEMCALTriggerFastOR.h"
#include "AliLog.h"

#include "AliEMCALTriggerConstants.h"
#include "AliEmcalTriggerQAAP.h"


using namespace EMCALTrigger;

/// \cond CLASSIMP
ClassImp(AliEmcalTriggerQAAP)
/// \endcond

const Int_t AliEmcalTriggerQAAP::fgkMaxPatchAmp[6] = {2000, 2000, 2000, 6000, 6000, 5000};
const TString AliEmcalTriggerQAAP::fgkPatchTypes[3] = {"Online", "Recalc", "Offline"};
/**
 * Dummy constructor
 */
AliEmcalTriggerQAAP::AliEmcalTriggerQAAP():
  TNamed(),
  fFastorL0Th(400),
  fFastorL1Th(400),
  fBkgPatchType(kTMEMCalBkg),
  fADCperBin(20),
  fDebugLevel(0),
  fHistos(0)
{
  for (Int_t i = 0; i < 3; i++) {
    fBkgADCAmpEMCal[i].Set(100);
    fBkgADCAmpDCal[i].Set(100);
    fNBkgPatchesEMCal[i] = 0;
    fNBkgPatchesDCal[i] = 0;
    fEnabledPatchTypes[i] = kTRUE;

    for (Int_t itype = 0; itype < 6; itype++) {
      fMaxPatchEMCal[itype][i] = 0;
      fMaxPatchDCal[itype][i] = 0;
    }
  }

  memset(fPatchAreas, 0, sizeof(Int_t)*6);
}

/**
 * Constructor
 */
AliEmcalTriggerQAAP::AliEmcalTriggerQAAP(const char* name):
  TNamed(name,name),
  fFastorL0Th(400),
  fFastorL1Th(400),
  fBkgPatchType(kTMEMCalBkg),
  fADCperBin(20),
  fDebugLevel(0),
  fHistos(0)
{
  for (Int_t i = 0; i < 3; i++) {
    fBkgADCAmpEMCal[i].Set(100);
    fBkgADCAmpDCal[i].Set(100);
    fNBkgPatchesEMCal[i] = 0;
    fNBkgPatchesDCal[i] = 0;
    fEnabledPatchTypes[i] = kTRUE;

    for (Int_t itype = 0; itype < 6; itype++) {
      fMaxPatchEMCal[itype][i] = 0;
      fMaxPatchDCal[itype][i] = 0;
    }
  }

  memset(fPatchAreas, 0, sizeof(Int_t)*6);
}

/**
 * Copy Constructor
 */
AliEmcalTriggerQAAP::AliEmcalTriggerQAAP(const AliEmcalTriggerQAAP& triggerQA) :
  TNamed(triggerQA),
  fFastorL0Th(triggerQA.fFastorL0Th),
  fFastorL1Th(triggerQA.fFastorL1Th),
  fBkgPatchType(triggerQA.fBkgPatchType),
  fADCperBin(triggerQA.fADCperBin),
  fDebugLevel(triggerQA.fDebugLevel),
  fHistos(0)
{
  for (Int_t i = 0; i < 3; i++) {
    fBkgADCAmpEMCal[i].Set(100);
    fBkgADCAmpDCal[i].Set(100);
    fNBkgPatchesEMCal[i] = 0;
    fNBkgPatchesDCal[i] = 0;
    fEnabledPatchTypes[i] = kTRUE;

    for (Int_t itype = 0; itype < 6; itype++) {
      fMaxPatchEMCal[itype][i] = 0;
      fMaxPatchDCal[itype][i] = 0;
    }
  }

  memset(fPatchAreas, 0, sizeof(Int_t)*6);
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
  fHistos->SetName(Form("%s_histos", GetName()));
  fHistos->SetOwner(kTRUE);

  hname = Form("EMCTRQA_histFastORL0");
  htitle = Form("EMCTRQA_histFastORL0;FastOR abs. ID;entries above 0");
  CreateTH1(hname, htitle, 5000, 0, 5000);

  hname = Form("EMCTRQA_histLargeAmpFastORL0");
  htitle = Form("EMCTRQA_histLargeAmpFastORL0 (>%d);FastOR abs. ID;entries above %d", fFastorL0Th, fFastorL0Th);
  CreateTH1(hname, htitle, 5000, 0, 5000);

  hname = Form("EMCTRQA_histFastORL0TimeOk");
  htitle = Form("EMCTRQA_histFastORL0TimeOk;FastOR abs. ID;entries (7 < time < 10)");
  CreateTH1(hname, htitle, 5000, 0, 5000);

  hname = Form("EMCTRQA_histFastORL0Mean");
  htitle = Form("EMCTRQA_histFastORL0Mean;FastOR abs. ID;mean ADC counts");
  CreateTProfile(hname, htitle, 5000, 0, 5000);

  hname = Form("EMCTRQA_histFastORL0MeanTimeOk");
  htitle = Form("EMCTRQA_histFastORL0MeanTimeOk;FastOR abs. ID;mean ADC counts  (7 < time < 10)");
  CreateTProfile(hname, htitle, 5000, 0, 5000);

  hname = Form("EMCTRQA_histFastORL1");
  htitle = Form("EMCTRQA_histFastORL1;FastOR abs. ID;entries above 0");
  CreateTH1(hname, htitle, 5000, 0, 5000);

  hname = Form("EMCTRQA_histLargeAmpFastORL1");
  htitle = Form("EMCTRQA_histLargeAmpFastORL1 (>%d);FastOR abs. ID;entries above %d", fFastorL1Th, fFastorL1Th);
  CreateTH1(hname, htitle, 5000, 0, 5000);

  hname = Form("EMCTRQA_histFastORL1Mean");
  htitle = Form("EMCTRQA_histFastORL1Mean;FastOR abs. ID;mean L1 time sum");
  CreateTProfile(hname, htitle, 5000, 0, 5000);

  hname = Form("EMCTRQA_histFastORL1AmpVsL0Amp");
  htitle = Form("EMCTRQA_histFastORL1AmpVsL0Amp;L0 amplitude;L1 time sum;entries");
  CreateTH2(hname, htitle, 256, 0, 1024, 256, 0, 1024);

  for (Int_t itype = 0; itype < 3; itype++) {
    if (!fEnabledPatchTypes[itype]) continue;
    hname = Form("EMCTRQA_histEMCalMedianVsDCalMedian%s", fgkPatchTypes[itype].Data());
    htitle = Form("EMCTRQA_histEMCalMedianVsDCalMedian%s;EMCal median;DCal median;entries", fgkPatchTypes[itype].Data());
    CreateTH2(hname, htitle, fgkMaxPatchAmp[fBkgPatchType]/fADCperBin/10, 0, fgkMaxPatchAmp[fBkgPatchType]/10, fgkMaxPatchAmp[fBkgPatchType]/fADCperBin/10, 0, fgkMaxPatchAmp[fBkgPatchType]/10);
  }

  const char* det[2] = { "EMCal", "DCal" };

  for (Int_t itrig = 1; itrig < 6; itrig++) {
    if (kEMCalTriggerNames[itrig].IsNull()) continue;
    for (Int_t itype = 0; itype < 3; itype++) {
      if (!fEnabledPatchTypes[itype]) continue;
      for (Int_t idet = 0; idet < 2; idet++) {
        hname = Form("EMCTRQA_hist%sPatchEnergy%s%s", det[idet], kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
        htitle = Form("EMCTRQA_hist%sPatchEnergy%s%s;energy (GeV);entries", det[idet], kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
        CreateTH1(hname, htitle, 200, 0, 200);

        hname = Form("EMCTRQA_hist%sPatchAmp%s%s", det[idet], kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
        htitle = Form("EMCTRQA_hist%sPatchAmp%s%s;amplitude;entries", det[idet], kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
        CreateTH1(hname, htitle, fgkMaxPatchAmp[itrig]/fADCperBin, 0, fgkMaxPatchAmp[itrig]);

        hname = Form("EMCTRQA_hist%sPatchAmpSubtracted%s%s", det[idet], kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
        htitle = Form("EMCTRQA_hist%sPatchAmpSubtracted%s%s;amplitude - background * area;entries", det[idet], kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
        CreateTH1(hname, htitle, fgkMaxPatchAmp[itrig]/fADCperBin, -fgkMaxPatchAmp[itrig]/2, fgkMaxPatchAmp[itrig]/2);

        hname = Form("EMCTRQA_hist%sMaxPatchAmp%s%s", det[idet], kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
        htitle = Form("EMCTRQA_hist%sMaxPatchAmp%s%s;amplitude;entries", det[idet], kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
        CreateTH1(hname, htitle, fgkMaxPatchAmp[itrig]/fADCperBin, 0, fgkMaxPatchAmp[itrig]);

        hname = Form("EMCTRQA_hist%sMaxPatchAmpSubtracted%s%s", det[idet], kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
        htitle = Form("EMCTRQA_hist%sMaxPatchAmpSubtracted%s%s;amplitude - background * area;entries", det[idet], kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
        CreateTH1(hname, htitle, fgkMaxPatchAmp[itrig]/fADCperBin, -fgkMaxPatchAmp[itrig]/2, fgkMaxPatchAmp[itrig]/2);

        for (Int_t jdet = 0; jdet < 2; jdet++) {
          hname = Form("EMCTRQA_hist%sMedianVs%sMax%s%s", det[idet], det[jdet], kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
          htitle = Form("EMCTRQA_hist%sMedianVs%sMax%s%s;%s max;%s median;entries", det[idet], det[jdet], kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data(), det[jdet], det[idet]);
          CreateTH2(hname, htitle, fgkMaxPatchAmp[itrig]/fADCperBin, 0, fgkMaxPatchAmp[itrig], fgkMaxPatchAmp[itrig]/fADCperBin/10, 0, fgkMaxPatchAmp[itrig]/10);
        }
      }

      hname = Form("EMCTRQA_histEMCalMaxVsDCalMax%s%s", kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
      htitle = Form("EMCTRQA_histEMCalMaxVsDCalMax%s%s;EMCal max;DCal max;entries", kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
      CreateTH2(hname, htitle, fgkMaxPatchAmp[itrig]/fADCperBin/4, 0, fgkMaxPatchAmp[itrig], fgkMaxPatchAmp[itrig]/fADCperBin/4, 0, fgkMaxPatchAmp[itrig]);

      hname = Form("EMCTRQA_histEdgePos%s%s", kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
      htitle = Form("EMCTRQA_histEdgePos%s%s;col;row;entries", kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
      CreateTH2(hname, htitle, 48, 0, 48, 105, 0, 105);

      //hname = Form("EMCTRQA_histCMPos%s%s", kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
      //htitle = Form("EMCTRQA_histCMPos%s%s;#eta;#phi;entries", kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
      //CreateTH2(hname, htitle, 60, -1, 1, 150, 0, TMath::TwoPi());

      hname = Form("EMCTRQA_histGeoPos%s%s", kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
      htitle = Form("EMCTRQA_histGeoPos%s%s;#eta;#phi;entries", kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
      CreateTH2(hname, htitle, 60, -1, 1, 150, 0, TMath::TwoPi());

      hname = Form("EMCTRQA_histLargeAmpEdgePos%s%s", kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
      htitle = Form("EMCTRQA_histLargeAmpEdgePos%s%s (>700);col;row;entries above 700", kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
      CreateTH2(hname, htitle, 48, 0, 48, 105, 0, 105);
    }
  }
}

/**
 * Process a patch, filling relevant histograms.
 * \param patch Pointer to a valid trigger patch
 */
void AliEmcalTriggerQAAP::ProcessPatch(AliEMCALTriggerPatchInfo* patch)
{
  TString hname;

  Int_t triggerBits[6] = { patch->GetTriggerBitConfig()->GetLevel0Bit(),
      patch->GetTriggerBitConfig()->GetGammaLowBit(),
      patch->GetTriggerBitConfig()->GetGammaHighBit(),
      patch->GetTriggerBitConfig()->GetJetLowBit(),
      patch->GetTriggerBitConfig()->GetJetHighBit(),
      patch->GetTriggerBitConfig()->GetBkgBit()
  };

  Int_t offsets[3] = { 0, AliEMCALTriggerPatchInfo::kRecalcOffset, AliEMCALTriggerPatchInfo::kOfflineOffset };
  Int_t amplitudes[3] = { patch->GetADCAmp(),  patch->GetADCAmp(),  patch->GetADCOfflineAmp() };
  Double_t bkg[3] = {0};

  for (Int_t itrig = 1; itrig < 6; itrig++) {
    if (kEMCalTriggerNames[itrig].IsNull()) continue;
    if (itrig == fBkgPatchType) continue;

    for (Int_t itype = 0; itype < 3; itype++) {
      if (!fEnabledPatchTypes[itype]) continue;
      if (!patch->TestTriggerBit(triggerBits[itrig]+offsets[itype])) continue;
      fPatchAreas[itrig] = patch->GetPatchSize()*patch->GetPatchSize();

      TString det;

      if (patch->IsEMCal()) {
        det = "EMCal";
        if (fMaxPatchEMCal[itrig][itype] < amplitudes[itype]) fMaxPatchEMCal[itrig][itype] = amplitudes[itype];
        bkg[itype] = fBkgDCal[itype];  // use DCal background for EMCal
      }
      else if (patch->IsDCalPHOS()) {
        det = "DCal";
        if (fMaxPatchDCal[itrig][itype] < amplitudes[itype]) fMaxPatchDCal[itrig][itype] = amplitudes[itype];
        bkg[itype] = fBkgEMCal[itype];  // use EMCal background for DCal
      }
      else {
        AliWarning(Form("Patch is not EMCal nor DCal/PHOS (pos: %d, %d)", patch->GetRowStart(), patch->GetColStart()));
      }

      hname = Form("EMCTRQA_histEdgePos%s%s", kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
      FillTH2(hname, patch->GetColStart(), patch->GetRowStart());

      //hname = Form("EMCTRQA_histCMPos%s%s", kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
      //FillTH2(hname, patch->GetEtaCM(), TVector2::Phi_0_2pi(patch->GetPhiCM()));

      hname = Form("EMCTRQA_histGeoPos%s%s", kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
      FillTH2(hname, patch->GetEtaGeo(), TVector2::Phi_0_2pi(patch->GetPhiGeo()));

      if (amplitudes[itype] > 700) {
        hname = Form("EMCTRQA_histLargeAmpEdgePos%s%s", kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
        FillTH2(hname, patch->GetColStart(), patch->GetRowStart());
      }

      hname = Form("EMCTRQA_hist%sPatchAmp%s%s", det.Data(), kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
      FillTH1(hname, amplitudes[itype]);

      hname = Form("EMCTRQA_hist%sPatchAmpSubtracted%s%s", det.Data(), kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
      FillTH1(hname, Double_t(amplitudes[itype]) - bkg[itype] * fPatchAreas[itrig]);

      hname = Form("EMCTRQA_hist%sPatchEnergy%s%s", det.Data(), kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
      FillTH1(hname, patch->GetPatchE());
    }

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

/**
 * Process a patch, filling relevant histograms.
 * \param patch Pointer to a valid trigger patch
 */
void AliEmcalTriggerQAAP::ProcessBkgPatch(AliEMCALTriggerPatchInfo* patch)
{
  TString hname;

  Int_t offsets[3] = { 0, AliEMCALTriggerPatchInfo::kRecalcOffset, AliEMCALTriggerPatchInfo::kOfflineOffset };
  Int_t amplitudes[3] = { patch->GetADCAmp(),  patch->GetADCAmp(),  patch->GetADCOfflineAmp() };
  Int_t bkgTriggerBit = patch->GetTriggerBitConfig()->GetBkgBit();

  for (Int_t itype = 0; itype < 3; itype++) {
    if (!fEnabledPatchTypes[itype]) continue;
    if (!patch->TestTriggerBit(bkgTriggerBit+offsets[itype])) continue;
    fPatchAreas[bkgTriggerBit] = patch->GetPatchSize()*patch->GetPatchSize();

    TString det;

    if (patch->IsEMCal()) {
      det = "EMCal";
      if (fMaxPatchEMCal[bkgTriggerBit][itype] < amplitudes[itype]) fMaxPatchEMCal[bkgTriggerBit][itype] = amplitudes[itype];

      if (fNBkgPatchesEMCal[itype] >= fBkgADCAmpEMCal[itype].GetSize()) {
        fBkgADCAmpEMCal[itype].Set((fNBkgPatchesEMCal[itype]+1)*2);
      }
      fBkgADCAmpEMCal[itype].AddAt(amplitudes[itype], fNBkgPatchesEMCal[itype]);
      fNBkgPatchesEMCal[itype]++;

    }
    else if (patch->IsDCalPHOS()) {
      det = "DCal";
      if (fMaxPatchDCal[bkgTriggerBit][itype] < amplitudes[itype]) fMaxPatchDCal[bkgTriggerBit][itype] = amplitudes[itype];

      if (fNBkgPatchesDCal[itype] >= fBkgADCAmpDCal[itype].GetSize()) {
        fBkgADCAmpDCal[itype].Set((fNBkgPatchesDCal[itype]+1)*2);
      }
      fBkgADCAmpDCal[itype].AddAt(amplitudes[itype], fNBkgPatchesDCal[itype]);
      fNBkgPatchesDCal[itype]++;
    }
    else {
      AliWarning(Form("Patch is not EMCal nor DCal/PHOS (pos: %d, %d)", patch->GetRowStart(), patch->GetColStart()));
    }

    hname = Form("EMCTRQA_histEdgePos%s%s", kEMCalTriggerNames[bkgTriggerBit].Data(), fgkPatchTypes[itype].Data());
    FillTH2(hname, patch->GetColStart(), patch->GetRowStart());

    //hname = Form("EMCTRQA_histCMPos%s%s", kEMCalTriggerNames[bkgTriggerBit].Data(), fgkPatchTypes[itype].Data());
    //FillTH2(hname, patch->GetEtaCM(), TVector2::Phi_0_2pi(patch->GetPhiCM()));

    hname = Form("EMCTRQA_histGeoPos%s%s", kEMCalTriggerNames[bkgTriggerBit].Data(), fgkPatchTypes[itype].Data());
    FillTH2(hname, patch->GetEtaGeo(), TVector2::Phi_0_2pi(patch->GetPhiGeo()));

    if (amplitudes[itype] > 700) {
      hname = Form("EMCTRQA_histLargeAmpEdgePos%s%s", kEMCalTriggerNames[bkgTriggerBit].Data(), fgkPatchTypes[itype].Data());
      FillTH2(hname, patch->GetColStart(), patch->GetRowStart());
    }

    hname = Form("EMCTRQA_hist%sPatchAmp%s%s", det.Data(), kEMCalTriggerNames[bkgTriggerBit].Data(), fgkPatchTypes[itype].Data());
    FillTH1(hname, amplitudes[itype]);

    hname = Form("EMCTRQA_hist%sPatchEnergy%s%s", det.Data(), kEMCalTriggerNames[bkgTriggerBit].Data(), fgkPatchTypes[itype].Data());
    FillTH1(hname, patch->GetPatchE());
  }

  if (fDebugLevel >= 2) {
    Printf("Type = %s; global pos = (%d, %d); Amp (online) = %d; Amp (offline) = %d; Patch energy = %.3f\n"
        "Position (CM): Eta=%.3f, Phi=%.3f\n"
        "Position (Geo): Eta=%.3f, Phi=%.3f\n",
        kEMCalTriggerNames[bkgTriggerBit].Data(), patch->GetRowStart(), patch->GetColStart(), patch->GetADCAmp(), patch->GetADCOfflineAmp(), patch->GetPatchE(),
        patch->GetEtaCM(), patch->GetPhiCM(),
        patch->GetEtaGeo(), patch->GetPhiGeo());
  }

}

/**
 * Process a FastOR, filling relevant histograms.
 * \param patch Pointer to a valid trigger FastOR
 */
void AliEmcalTriggerQAAP::ProcessFastor(AliEMCALTriggerFastOR* fastor)
{
  TString hname;

  if (fastor->GetL0Amp() > 0) {
    hname = Form("EMCTRQA_histFastORL0");
    FillTH1(hname, fastor->GetAbsId());

    hname = Form("EMCTRQA_histFastORL0Mean");
    FillTProfile(hname, fastor->GetAbsId(), fastor->GetL0Amp());

    if (fastor->GetL0Time() < 10 && fastor->GetL0Time() > 7) {
      hname = Form("EMCTRQA_histFastORL0TimeOk");
      FillTH1(hname, fastor->GetAbsId());

      hname = Form("EMCTRQA_histFastORL0MeanTimeOk");
      FillTProfile(hname, fastor->GetAbsId(), fastor->GetL0Amp());
    }
  }

  if (fastor->GetL0Amp() > fFastorL0Th) {
    hname = Form("EMCTRQA_histLargeAmpFastORL0");
    FillTH1(hname, fastor->GetAbsId());
  }

  if (fastor->GetL1Amp() > 0) {
    hname = Form("EMCTRQA_histFastORL1");
    FillTH1(hname, fastor->GetAbsId());

    hname = Form("EMCTRQA_histFastORL1Mean");
    FillTProfile(hname, fastor->GetAbsId(), fastor->GetL1Amp());
  }

  if (fastor->GetL1Amp() > fFastorL1Th) {
    hname = Form("EMCTRQA_histLargeAmpFastORL1");
    FillTH1(hname, fastor->GetAbsId());
  }

  if (fastor->GetL1Amp() > 0 && fastor->GetL0Amp() > 0) {
    hname = Form("EMCTRQA_histFastORL1AmpVsL0Amp");
    FillTH2(hname, fastor->GetL0Amp(), fastor->GetL1Amp());
  }
}

/**
 * Computes background for the current event
 */
void AliEmcalTriggerQAAP::ComputeBackground()
{
  AliDebug(2, Form("Entering AliEmcalTriggerQAAP::ComputeBackground"));

  for (Int_t itype = 0; itype < 3; itype++) {
    if (!fEnabledPatchTypes[itype]) continue;

    AliDebug(2, Form("Patch type %s", fgkPatchTypes[itype].Data()));

    fMedianEMCal[itype] = TMath::Median(fNBkgPatchesEMCal[itype], fBkgADCAmpEMCal[itype].GetArray());
    fMedianDCal[itype] = TMath::Median(fNBkgPatchesDCal[itype], fBkgADCAmpDCal[itype].GetArray());

    fBkgEMCal[itype] = fMedianEMCal[itype] / fPatchAreas[fBkgPatchType];
    fBkgDCal[itype] = fMedianDCal[itype] / fPatchAreas[fBkgPatchType];
  }
}

/**
 * This method should be called at the end of each event.
 */
void AliEmcalTriggerQAAP::EventCompleted()
{
  AliDebug(2, Form("Entering AliEmcalTriggerQAAP::EventCompleted"));

  TString hname;

  for (Int_t itype = 0; itype < 3; itype++) {
    if (!fEnabledPatchTypes[itype]) continue;

    AliDebug(2, Form("Patch type %s", fgkPatchTypes[itype].Data()));

    hname = Form("EMCTRQA_histEMCalMedianVsDCalMedian%s", fgkPatchTypes[itype].Data());
    FillTH2(hname, fMedianEMCal[itype], fMedianDCal[itype]);

    for (Int_t itrig = 1; itrig < 6; itrig++) {
      if (kEMCalTriggerNames[itrig].IsNull()) continue;
      if (fMaxPatchDCal[itrig][itype] == 0 && fMaxPatchEMCal[itrig][itype] == 0) continue;

      AliDebug(2, Form("Trigger type: %s", kEMCalTriggerNames[itype].Data()));

      hname = Form("EMCTRQA_histEMCalMedianVsDCalMax%s%s", kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
      FillTH2(hname, fMaxPatchDCal[itrig][itype], fMedianEMCal[itype]);

      hname = Form("EMCTRQA_histDCalMedianVsEMCalMax%s%s", kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
      FillTH2(hname, fMaxPatchEMCal[itrig][itype], fMedianDCal[itype]);

      hname = Form("EMCTRQA_histDCalMedianVsDCalMax%s%s", kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
      FillTH2(hname, fMaxPatchDCal[itrig][itype], fMedianDCal[itype]);

      hname = Form("EMCTRQA_histEMCalMedianVsEMCalMax%s%s", kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
      FillTH2(hname, fMaxPatchEMCal[itrig][itype], fMedianEMCal[itype]);

      hname = Form("EMCTRQA_histEMCalMaxVsDCalMax%s%s", kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
      FillTH2(hname, fMaxPatchEMCal[itrig][itype], fMaxPatchDCal[itrig][itype]);

      if (fMaxPatchEMCal[itrig][itype] > 0) {
        hname = Form("EMCTRQA_histEMCalMaxPatchAmpSubtracted%s%s", kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
        FillTH1(hname, Double_t(fMaxPatchEMCal[itrig][itype]) - fBkgDCal[itype] * fPatchAreas[itrig]);

        hname = Form("EMCTRQA_histEMCalMaxPatchAmp%s%s", kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
        FillTH1(hname, fMaxPatchEMCal[itrig][itype]);
      }

      if (fMaxPatchDCal[itrig][itype] > 0) {
        hname = Form("EMCTRQA_histDCalMaxPatchAmpSubtracted%s%s", kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
        FillTH1(hname, Double_t(fMaxPatchDCal[itrig][itype]) - fBkgEMCal[itype] * fPatchAreas[itrig]);

        hname = Form("EMCTRQA_histDCalMaxPatchAmp%s%s", kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
        FillTH1(hname, fMaxPatchDCal[itrig][itype]);
      }

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
