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
 * @file AliEMCALTriggerQA.cxx
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

#include "AliEMCALTriggerPatchInfo.h"
#include "AliEMCALTriggerFastOR.h"
#include "AliLog.h"

#include "AliEMCALTriggerConstants.h"

#include "AliEMCALTriggerQA.h"

using namespace EMCALTrigger;

/// \cond CLASSIMP
ClassImp(AliEMCALTriggerQA)
/// \endcond

const Int_t AliEMCALTriggerQA::fgkMaxPatchAmp[32] = {
    1500, 1500, 1500, 2000, 2000,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    1500, 1500, 2000, 2000,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

/**
 * Dummy constructor
 */
AliEMCALTriggerQA::AliEMCALTriggerQA():
TNamed(),
fFastorL0Th(1),
fFastorL1Th(1),
fBkgPatchType(kEMCalRecalcL1Bkg),
fDebugLevel(0),
fHistos(0)
{
  fBkgADCAmpEMCal[0].Set(100);
  fBkgADCAmpDCal[0].Set(100);

  fBkgADCAmpEMCal[1].Set(100);
  fBkgADCAmpDCal[1].Set(100);

  fNBkgPatchesEMCal[0] = 0;
  fNBkgPatchesDCal[0] = 0;

  fNBkgPatchesEMCal[1] = 0;
  fNBkgPatchesDCal[1] = 0;

  for (int itype = 0; itype < 32; itype++) {
    fMaxPatchEMCal[itype][0] = 0;
    fMaxPatchEMCal[itype][1] = 0;
    fMaxPatchDCal[itype][0] = 0;
    fMaxPatchDCal[itype][1] = 0;
    fTrigThr[itype] = 0;
  }
}

/**
 * Constructor
 */
AliEMCALTriggerQA::AliEMCALTriggerQA(const char* name):
		      TNamed(name,name),
		      fFastorL0Th(1),
		      fFastorL1Th(1),
		      fBkgPatchType(kEMCalRecalcL1Bkg),
		      fDebugLevel(0),
		      fHistos(0)
{
  fBkgADCAmpEMCal[0].Set(100);
  fBkgADCAmpDCal[0].Set(100);

  fBkgADCAmpEMCal[1].Set(100);
  fBkgADCAmpDCal[1].Set(100);

  fNBkgPatchesEMCal[0] = 0;
  fNBkgPatchesDCal[0] = 0;

  fNBkgPatchesEMCal[1] = 0;
  fNBkgPatchesDCal[1] = 0;

  for (int itype = 0; itype < 32; itype++) {
    fMaxPatchEMCal[itype][0] = 0;
    fMaxPatchEMCal[itype][1] = 0;
    fMaxPatchDCal[itype][0] = 0;
    fMaxPatchDCal[itype][1] = 0;
    fTrigThr[itype] = 0;
  }
}

/**
 * Destructor
 */
AliEMCALTriggerQA::~AliEMCALTriggerQA()
{
}

/**
 * Initialize the class, i.e. allocate histograms.
 */
void AliEMCALTriggerQA::Init()
{
  TString hname;
  TString htitle;

  const char *patchtypes[2] = {"Online", "Offline"};

  fHistos = new THashList();
  fHistos->SetName(Form("histos%s", GetName()));
  fHistos->SetOwner(kTRUE);

  hname = Form("EMCTRQA_histFastORL0");
  htitle = Form("EMCTRQA_histFastORL0;FastOR abs. ID;entries above L0 threshold");
  CreateTH1(hname, htitle, 4000, 0, 4000);

  hname = Form("EMCTRQA_histFastORL0Mean");
  htitle = Form("EMCTRQA_histFastORL0Mean;FastOR abs. ID;mean ADC counts");
  CreateTProfile(hname, htitle, 4000, 0, 4000);

  hname = Form("EMCTRQA_histLargeAmpFastORL0");
  htitle = Form("EMCTRQA_histLargeAmpFastORL0 (>400);FastOR abs. ID;entries above 400");
  CreateTH1(hname, htitle, 4000, 0, 4000);

  hname = Form("EMCTRQA_histFastORL1");
  htitle = Form("EMCTRQA_histFastORL1;FastOR abs. ID;entries above L1 threshold");
  CreateTH1(hname, htitle, 4000, 0, 4000);

  hname = Form("EMCTRQA_histFastORL1Mean");
  htitle = Form("EMCTRQA_histFastORL1Mean;FastOR abs. ID;mean L1 time sum");
  CreateTProfile(hname, htitle, 4000, 0, 4000);

  hname = Form("EMCTRQA_histLargeAmpFastORL1");
  htitle = Form("EMCTRQA_histLargeAmpFastORL1 (>400);FastOR abs. ID;entries above 400");
  CreateTH1(hname, htitle, 4000, 0, 4000);

  hname = Form("EMCTRQA_histFastORL1AmpVsL0Amp");
  htitle = Form("EMCTRQA_histFastORL1AmpVsL0Amp;L0 amplitude;L1 time sum;entries");
  CreateTH2(hname, htitle, 128, 0, 2048, 128, 0, 2048);

  hname = Form("EMCTRQA_histEMCalMedianVsDCalMedian%s", patchtypes[0]);
  htitle = Form("EMCTRQA_histEMCalMedianVsDCalMedian%s;EMCal median;DCal median;entries", patchtypes[0]);
  CreateTH2(hname, htitle, fgkMaxPatchAmp[fBkgPatchType]/10, 0, fgkMaxPatchAmp[fBkgPatchType], fgkMaxPatchAmp[fBkgPatchType]/10, 0, fgkMaxPatchAmp[fBkgPatchType]);

  hname = Form("EMCTRQA_histEMCalMedianVsDCalMedian%s", patchtypes[1]);
  htitle = Form("EMCTRQA_histEMCalMedianVsDCalMedian%s;EMCal median;DCal median;entries", patchtypes[1]);
  CreateTH2(hname, htitle, fgkMaxPatchAmp[fBkgPatchType]/10, 0, fgkMaxPatchAmp[fBkgPatchType], fgkMaxPatchAmp[fBkgPatchType]/10, 0, fgkMaxPatchAmp[fBkgPatchType]);

  for (int itype = 0; itype < 32; itype++) {
    if (kEMCalTriggerNames[itype].IsNull()) continue;
    for (const char **patchtype = patchtypes; patchtype < patchtypes + 2; ++patchtype) {
      hname = Form("EMCTRQA_histEMCalPatchAmp%s%s", kEMCalTriggerNames[itype].Data(), *patchtype);
      htitle = Form("EMCTRQA_histEMCalPatchAmp%s%s;amplitude;entries", kEMCalTriggerNames[itype].Data(), *patchtype);
      CreateTH1(hname, htitle, fgkMaxPatchAmp[itype]/10, 0, fgkMaxPatchAmp[itype]);

      hname = Form("EMCTRQA_histDCalPatchAmp%s%s", kEMCalTriggerNames[itype].Data(), *patchtype);
      htitle = Form("EMCTRQA_histDCalPatchAmp%s%s;amplitude;entries", kEMCalTriggerNames[itype].Data(), *patchtype);
      CreateTH1(hname, htitle, fgkMaxPatchAmp[itype]/10, 0, fgkMaxPatchAmp[itype]);

      hname = Form("EMCTRQA_histEMCalMedianVsDCalMax%s%s", kEMCalTriggerNames[itype].Data(), *patchtype);
      htitle = Form("EMCTRQA_histEMCalMedianVsDCalMax%s%s;DCal max;EMCal median;entries", kEMCalTriggerNames[itype].Data(), *patchtype);
      CreateTH2(hname, htitle, fgkMaxPatchAmp[itype]/10, 0, fgkMaxPatchAmp[itype], fgkMaxPatchAmp[itype]/10, 0, fgkMaxPatchAmp[itype]);

      hname = Form("EMCTRQA_histDCalMedianVsEMCalMax%s%s", kEMCalTriggerNames[itype].Data(), *patchtype);
      htitle = Form("EMCTRQA_histDCalMedianVsEMCalMax%s%s;EMCal max;DCal median;entries", kEMCalTriggerNames[itype].Data(), *patchtype);
      CreateTH2(hname, htitle, fgkMaxPatchAmp[itype]/10, 0, fgkMaxPatchAmp[itype], fgkMaxPatchAmp[itype]/10, 0, fgkMaxPatchAmp[itype]);

      hname = Form("EMCTRQA_histEMCalMedianVsEMCalMax%s%s", kEMCalTriggerNames[itype].Data(), *patchtype);
      htitle = Form("EMCTRQA_histEMCalMedianVsEMCalMax%s%s;EMCal max;EMCal median;entries", kEMCalTriggerNames[itype].Data(), *patchtype);
      CreateTH2(hname, htitle, fgkMaxPatchAmp[itype]/10, 0, fgkMaxPatchAmp[itype], fgkMaxPatchAmp[itype]/10, 0, fgkMaxPatchAmp[itype]);

      hname = Form("EMCTRQA_histDCalMedianVsDCalMax%s%s", kEMCalTriggerNames[itype].Data(), *patchtype);
      htitle = Form("EMCTRQA_histDCalMedianVsDCalMax%s%s;DCal max;DCal median;entries", kEMCalTriggerNames[itype].Data(), *patchtype);
      CreateTH2(hname, htitle, fgkMaxPatchAmp[itype]/10, 0, fgkMaxPatchAmp[itype], fgkMaxPatchAmp[itype]/10, 0, fgkMaxPatchAmp[itype]);

      hname = Form("EMCTRQA_histEMCalMaxVsDCalMax%s%s", kEMCalTriggerNames[itype].Data(), *patchtype);
      htitle = Form("EMCTRQA_histEMCalMaxVsDCalMax%s%s;EMCal max;DCal max;entries", kEMCalTriggerNames[itype].Data(), *patchtype);
      CreateTH2(hname, htitle, fgkMaxPatchAmp[itype]/10, 0, fgkMaxPatchAmp[itype], fgkMaxPatchAmp[itype]/10, 0, fgkMaxPatchAmp[itype]);
    }

    hname = Form("EMCTRQA_histEdgePos%s", kEMCalTriggerNames[itype].Data());
    htitle = Form("EMCTRQA_histEdgePos%s;col;row;entries", kEMCalTriggerNames[itype].Data());
    CreateTH2(hname, htitle, 48, 0, 48, 105, 0, 105);

    hname = Form("EMCTRQA_histCMPos%s", kEMCalTriggerNames[itype].Data());
    htitle = Form("EMCTRQA_histCMPos%s;#eta;#phi;entries", kEMCalTriggerNames[itype].Data());
    CreateTH2(hname, htitle, 60, -1, 1, 200, 0, TMath::TwoPi());

    hname = Form("EMCTRQA_histGeoPos%s", kEMCalTriggerNames[itype].Data());
    htitle = Form("EMCTRQA_histGeoPos%s;#eta;#phi;entries", kEMCalTriggerNames[itype].Data());
    CreateTH2(hname, htitle, 60, -1, 1, 200, 0, TMath::TwoPi());

    hname = Form("EMCTRQA_histLargeAmpEdgePos%s", kEMCalTriggerNames[itype].Data());
    htitle = Form("EMCTRQA_histLargeAmpEdgePos%s (>700);col;row;entries above 700", kEMCalTriggerNames[itype].Data());
    CreateTH2(hname, htitle, 48, 0, 48, 105, 0, 105);

    hname = Form("EMCTRQA_histEMCalPatchEnergy%s", kEMCalTriggerNames[itype].Data());
    htitle = Form("EMCTRQA_histEMCalPatchEnergy%s;energy (GeV);entries", kEMCalTriggerNames[itype].Data());
    CreateTH1(hname, htitle, 200, 0, 200);

    hname = Form("EMCTRQA_histDCalPatchEnergy%s", kEMCalTriggerNames[itype].Data());
    htitle = Form("EMCTRQA_histDCalPatchEnergy%s;energy (GeV);entries", kEMCalTriggerNames[itype].Data());
    CreateTH1(hname, htitle, 200, 0, 200);
  }
}

/**
 * Process a patch, filling relevant histograms.
 * \param patch Pointer to a valid trigger patch
 */
void AliEMCALTriggerQA::ProcessPatch(AliEMCALTriggerPatchInfo* patch)
{
  TString hname;

  for (int type = 0; type < 32; type++) {
    if (kEMCalTriggerNames[type].IsNull()) continue;
    if (!patch->TestTriggerBit(type)) continue;

    if (patch->GetADCOfflineAmp() > fTrigThr[type]) {
      hname = Form("EMCTRQA_histEdgePos%s", kEMCalTriggerNames[type].Data());
      FillTH2(hname, patch->GetColStart(), patch->GetRowStart());

      hname = Form("EMCTRQA_histCMPos%s", kEMCalTriggerNames[type].Data());
      FillTH2(hname, patch->GetEtaCM(), patch->GetPhiCM());

      hname = Form("EMCTRQA_histGeoPos%s", kEMCalTriggerNames[type].Data());
      FillTH2(hname, patch->GetEtaGeo(), patch->GetPhiGeo());

      if (patch->GetADCAmp() > 700 || patch->GetADCOfflineAmp() > 700) {
        hname = Form("EMCTRQA_histLargeAmpEdgePos%s", kEMCalTriggerNames[type].Data());
        FillTH2(hname, patch->GetColStart(), patch->GetRowStart());
      }
    }

    TString det;

    Int_t amp = patch->GetADCAmp();
    Int_t offlineamp = patch->GetADCOfflineAmp();

    if (patch->IsEMCal()) {
      det = "EMCal";
      if (fMaxPatchEMCal[type][0] < amp) fMaxPatchEMCal[type][0] = amp;
      if (fMaxPatchEMCal[type][1] < offlineamp) fMaxPatchEMCal[type][1] = offlineamp;

      if (type == fBkgPatchType) {
        if (amp > fTrigThr[fBkgPatchType]) {
          if (fNBkgPatchesEMCal[0] >= fBkgADCAmpEMCal[0].GetSize()) {
            fBkgADCAmpEMCal[0].Set((fNBkgPatchesEMCal[0]+1)*2);
          }
          fBkgADCAmpEMCal[0].AddAt(amp, fNBkgPatchesEMCal[0]);
          fNBkgPatchesEMCal[0]++;
        }

        if (offlineamp > fTrigThr[fBkgPatchType]) {
          if (fNBkgPatchesEMCal[1] >= fBkgADCAmpEMCal[1].GetSize()) {
            fBkgADCAmpEMCal[1].Set((fNBkgPatchesEMCal[1]+1)*2);
          }
          fBkgADCAmpEMCal[1].AddAt(offlineamp, fNBkgPatchesEMCal[1]);
          fNBkgPatchesEMCal[1]++;
        }
      }
    }
    else if (patch->IsDCalPHOS()) {
      det = "DCal";
      if (fMaxPatchDCal[type][0] < amp) fMaxPatchDCal[type][0] = amp;
      if (fMaxPatchDCal[type][1] < offlineamp) fMaxPatchDCal[type][1] = offlineamp;

      if (type == fBkgPatchType) {
        if (amp > fTrigThr[fBkgPatchType]) {
          if (fNBkgPatchesDCal[0] >= fBkgADCAmpDCal[0].GetSize()) {
            fBkgADCAmpDCal[0].Set((fNBkgPatchesDCal[0]+1)*2);
          }
          fBkgADCAmpDCal[0].AddAt(amp, fNBkgPatchesDCal[0]);
          fNBkgPatchesDCal[0]++;
        }

        if (offlineamp > fTrigThr[fBkgPatchType]) {
          if (fNBkgPatchesDCal[1] >= fBkgADCAmpDCal[1].GetSize()) {
            fBkgADCAmpDCal[1].Set((fNBkgPatchesDCal[1]+1)*2);
          }
          fBkgADCAmpDCal[1].AddAt(offlineamp, fNBkgPatchesDCal[1]);
          fNBkgPatchesDCal[1]++;
        }
      }
    }
    else {
      AliWarning(Form("Patch is not EMCal nor DCal/PHOS (pos: %d, %d)", patch->GetRowStart(), patch->GetColStart()));
    }

    if (amp > 0) {
      hname = Form("EMCTRQA_hist%sPatchAmp%sOnline", det.Data(), kEMCalTriggerNames[type].Data());
      FillTH1(hname, amp);
    }

    if (offlineamp > 0) {
      hname = Form("EMCTRQA_hist%sPatchAmp%sOffline", det.Data(), kEMCalTriggerNames[type].Data());
      FillTH1(hname, offlineamp);
    }

    if (patch->GetPatchE() > 0) {
      hname = Form("EMCTRQA_hist%sPatchEnergy%s", det.Data(), kEMCalTriggerNames[type].Data());
      FillTH1(hname, patch->GetPatchE());
    }

    if (fDebugLevel >= 2) {
      Printf("Type = %s; global pos = (%d, %d); Amp (online) = %d; Amp (offline) = %d; Patch energy = %.3f\n"
          "Position (CM): Eta=%.3f, Phi=%.3f\n"
          "Position (Geo): Eta=%.3f, Phi=%.3f\n",
          kEMCalTriggerNames[type].Data(), patch->GetRowStart(), patch->GetColStart(), patch->GetADCAmp(), patch->GetADCOfflineAmp(), patch->GetPatchE(),
          patch->GetEtaCM(), patch->GetPhiCM(),
          patch->GetEtaGeo(), patch->GetPhiGeo());
    }
  }
}

/**
 * Process a FastOR, filling relevant histograms.
 * \param patch Pointer to a valid trigger FastOR
 */
void AliEMCALTriggerQA::ProcessFastor(AliEMCALTriggerFastOR* fastor)
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
void AliEMCALTriggerQA::EventCompleted()
{
  TString hname;

  Double_t medianEMCal[2] = {0};
  Double_t medianDCal[2] = {0};

  for (Int_t i = 0; i < 2; i++) {
    medianEMCal[i] = TMath::Median(fNBkgPatchesEMCal[i], fBkgADCAmpEMCal[i].GetArray());
    medianDCal[i] = TMath::Median(fNBkgPatchesDCal[i], fBkgADCAmpDCal[i].GetArray());
  }

  const char *patchtypes[2] = {"Online", "Offline"};

  for (Int_t i = 0; i < 2; i++) {
    hname = Form("EMCTRQA_histEMCalMedianVsDCalMedian%s", patchtypes[i]);
    FillTH2(hname, medianEMCal[i], medianDCal[i]);

    for (int itype = 0; itype < 32; itype++) {
      if (kEMCalTriggerNames[itype].IsNull()) continue;
      if (fMaxPatchDCal[itype][i] == 0 && fMaxPatchEMCal[itype][i] == 0) continue;

      hname = Form("EMCTRQA_histEMCalMedianVsDCalMax%s%s", kEMCalTriggerNames[itype].Data(), patchtypes[i]);
      FillTH2(hname, medianEMCal[i], fMaxPatchDCal[itype][i]);

      hname = Form("EMCTRQA_histDCalMedianVsEMCalMax%s%s", kEMCalTriggerNames[itype].Data(), patchtypes[i]);
      FillTH2(hname, medianDCal[i], fMaxPatchEMCal[itype][i]);

      hname = Form("EMCTRQA_histDCalMedianVsDCalMax%s%s", kEMCalTriggerNames[itype].Data(), patchtypes[i]);
      FillTH2(hname, medianDCal[i], fMaxPatchDCal[itype][i]);

      hname = Form("EMCTRQA_histEMCalMedianVsEMCalMax%s%s", kEMCalTriggerNames[itype].Data(), patchtypes[i]);
      FillTH2(hname, medianEMCal[i], fMaxPatchEMCal[itype][i]);

      hname = Form("EMCTRQA_histEMCalMaxVsDCalMax%s%s", kEMCalTriggerNames[itype].Data(), patchtypes[i]);
      FillTH2(hname, fMaxPatchEMCal[itype][i], fMaxPatchDCal[itype][i]);

      fMaxPatchEMCal[itype][i] = 0;
      fMaxPatchDCal[itype][i] = 0;
    }
    fNBkgPatchesEMCal[i] = 0;
    fNBkgPatchesDCal[i] = 0;
    fBkgADCAmpEMCal[i].Reset();
    fBkgADCAmpDCal[i].Reset();
  }
}

//______________________________________________________________________________
void AliEMCALTriggerQA::CreateTProfile(const char *name, const char *title, int nbins, double xmin, double xmax)
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
    Fatal("AliEMCALTriggerQA::CreateTProfile", "Object %s already exists", name);
    return;
  }
  TProfile* hist = new TProfile(name, title, nbins, xmin, xmax);
  fHistos->Add(hist);
}

//______________________________________________________________________________
void AliEMCALTriggerQA::CreateTH1(const char *name, const char *title, int nbins, double xmin, double xmax)
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
    Fatal("AliEMCALTriggerQA::CreateTH1", "Object %s already exists", name);
    return;
  }
  TH1* hist = new TH1D(name, title, nbins, xmin, xmax);
  fHistos->Add(hist);
}

//______________________________________________________________________________
void AliEMCALTriggerQA::CreateTH2(const char *name, const char *title, int nbinsx, double xmin, double xmax, int nbinsy, double ymin, double ymax)
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
    Fatal("AliEMCALTriggerQA::CreateTH2", "Object %s already exists", name);
    return;
  }
  TH2* hist = new TH2D(name, title, nbinsx, xmin, xmax, nbinsy, ymin, ymax);
  fHistos->Add(hist);
}

//______________________________________________________________________________
void AliEMCALTriggerQA::CreateTH3(const char* name, const char* title, int nbinsx, double xmin, double xmax, int nbinsy, double ymin, double ymax, int nbinsz, double zmin, double zmax)
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
    Fatal("AliEMCALTriggerQA::CreateTH3", "Object %s already exists", name);
    return;
  }
  TH3* hist = new TH3D(name, title, nbinsx, xmin, xmax, nbinsy, ymin, ymax, nbinsz, zmin, zmax);
  fHistos->Add(hist);
}

//______________________________________________________________________________
void AliEMCALTriggerQA::FillTProfile(const char *name, double x, double y, double weight)
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
    Fatal("AliEMCALTriggerQA::FillTProfile", "Histogram %s not found", name);
    return;
  }
  hist->Fill(x, y, weight);
}

//______________________________________________________________________________
void AliEMCALTriggerQA::FillTH1(const char *name, double x, double weight)
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
    Fatal("AliEMCALTriggerQA::FillTH1", "Histogram %s not found", name);
    return;
  }
  hist->Fill(x, weight);
}

//______________________________________________________________________________
void AliEMCALTriggerQA::FillTH2(const char *name, double x, double y, double weight)
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
    Fatal("AliEMCALTriggerQA::FillTH2", "Histogram %s not found", name);
    return;
  }
  hist->Fill(x, y, weight);
}

//______________________________________________________________________________
void AliEMCALTriggerQA::FillTH3(const char* name, double x, double y, double z, double weight)
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
    Fatal("AliEMCALTriggerQA::FillTH3", "Histogram %s not found", name);
    return;
  }
  hist->Fill(x, y, z, weight);
}

//______________________________________________________________________________
TObject *AliEMCALTriggerQA::FindObject(const char *name) const
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
TObject* AliEMCALTriggerQA::FindObject(const TObject* obj) const
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
