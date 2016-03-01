/**************************************************************************
 * Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
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

#include <cstring>
#include <fstream>

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

#include "AliEmcalTriggerQAPP.h"

using namespace EMCALTrigger;

/// \cond CLASSIMP
ClassImp(AliEmcalTriggerQAPP)
/// \endcond

const Int_t AliEmcalTriggerQAPP::fgkMaxPatchAmp[6] = {3000, 3000, 3000, 6000, 6000, 5000};
const TString AliEmcalTriggerQAPP::fgkPatchTypes[3] = {"Online", "Recalc", "Offline"};

/// Dummy constructor for ROOT I/O
AliEmcalTriggerQAPP::AliEmcalTriggerQAPP():
  TNamed(),
  fFastorL0Th(400),
  fFastorL1Th(400),
  fADCperBin(16),
  fDebugLevel(0),
  fL0MinTime(7),
  fL0MaxTime(10),
  fHistManager()
{
  for (Int_t i = 0; i < 3; i++) {
    fEnabledPatchTypes[i] = kTRUE;

    for (Int_t itype = 0; itype < 6; itype++) {
      fMaxPatchEMCal[itype][i] = 0;
      fMaxPatchDCal[itype][i] = 0;
    }
  }

  fEnabledTriggerTypes[EMCALTrigger::kTMEMCalBkg] = kFALSE;
  fEnabledTriggerTypes[EMCALTrigger::kTMEMCalLevel0] = kTRUE;
  fEnabledTriggerTypes[EMCALTrigger::kTMEMCalGammaL] = kFALSE;
  fEnabledTriggerTypes[EMCALTrigger::kTMEMCalJetL] = kFALSE;
  fEnabledTriggerTypes[EMCALTrigger::kTMEMCalGammaH] = kTRUE;
  fEnabledTriggerTypes[EMCALTrigger::kTMEMCalJetH] = kTRUE;
}

/// Default constructor
///
/// \param name Name of the object
AliEmcalTriggerQAPP::AliEmcalTriggerQAPP(const char* name):
  TNamed(name,name),
  fFastorL0Th(400),
  fFastorL1Th(400),
  fADCperBin(16),
  fDebugLevel(0),
  fL0MinTime(7),
  fL0MaxTime(10),
  fHistManager(name)
{
  for (Int_t i = 0; i < 3; i++) {
    fEnabledPatchTypes[i] = kTRUE;

    for (Int_t itype = 0; itype < 6; itype++) {
      fMaxPatchEMCal[itype][i] = 0;
      fMaxPatchDCal[itype][i] = 0;
    }
  }

  fEnabledTriggerTypes[EMCALTrigger::kTMEMCalBkg] = kFALSE;
  fEnabledTriggerTypes[EMCALTrigger::kTMEMCalLevel0] = kTRUE;
  fEnabledTriggerTypes[EMCALTrigger::kTMEMCalGammaL] = kFALSE;
  fEnabledTriggerTypes[EMCALTrigger::kTMEMCalJetL] = kFALSE;
  fEnabledTriggerTypes[EMCALTrigger::kTMEMCalGammaH] = kTRUE;
  fEnabledTriggerTypes[EMCALTrigger::kTMEMCalJetH] = kTRUE;
}

/// Copy constructor
///
/// \param triggerQA Reference to an object to copy from
AliEmcalTriggerQAPP::AliEmcalTriggerQAPP(const AliEmcalTriggerQAPP& triggerQA) :
  TNamed(triggerQA),
  fFastorL0Th(triggerQA.fFastorL0Th),
  fFastorL1Th(triggerQA.fFastorL1Th),
  fADCperBin(triggerQA.fADCperBin),
  fDebugLevel(triggerQA.fDebugLevel),
  fL0MinTime(triggerQA.fL0MinTime),
  fL0MaxTime(triggerQA.fL0MaxTime),
  fHistManager(triggerQA.GetName())
{
  for (Int_t i = 0; i < 3; i++) {
    fEnabledPatchTypes[i] = kTRUE;

    for (Int_t itype = 0; itype < 6; itype++) {
      fMaxPatchEMCal[itype][i] = 0;
      fMaxPatchDCal[itype][i] = 0;
    }
  }
}

/// Destructor
AliEmcalTriggerQAPP::~AliEmcalTriggerQAPP()
{
}

/// Read the offline bad channel map from a standard stream
///
/// \param stream A reference to a standard stream to read from (can be a file stream)
void AliEmcalTriggerQAPP::ReadOfflineBadChannelFromStream(std::istream& stream)
{
  Short_t absId = 0;

  while (stream.good()) {
    stream >> absId;
    AddOfflineBadChannel(absId);
  }
}


/// Read the offline bad channel map from a text file
///
/// \param fname Path and name of the file
void AliEmcalTriggerQAPP::ReadOfflineBadChannelFromFile(const char* fname)
{
  std::ifstream file(fname);
  ReadOfflineBadChannelFromStream(file);
}

/// Set the patch types to be plotted
///
/// \param type Patch type of which the status is being changed
/// \param e    Either enable or disable
void AliEmcalTriggerQAPP::EnablePatchType(PatchTypes_t type, Bool_t e)
{
  fEnabledPatchTypes[type] = e;
}

 /// Initialize the class, i.e. allocate histograms.
void AliEmcalTriggerQAPP::Init()
{
  TString hname;
  TString htitle;

  hname = Form("EMCTRQA_histFastORL0");
  htitle = Form("EMCTRQA_histFastORL0;FastOR abs. ID;entries above 0");
  fHistManager.CreateTH1(hname, htitle, 5000, 0, 5000);

  hname = Form("EMCTRQA_histLargeAmpFastORL0");
  htitle = Form("EMCTRQA_histLargeAmpFastORL0 (>%d);FastOR abs. ID;entries above %d", fFastorL0Th, fFastorL0Th);
  fHistManager.CreateTH1(hname, htitle, 5000, 0, 5000);

  hname = Form("EMCTRQA_histFastORL0Amp");
  htitle = Form("EMCTRQA_histFastORL0Amp;FastOR abs. ID;ADC counts");
  fHistManager.CreateTH2(hname, htitle, 5000, 0, 5000, 1024, 0, 4096);

  hname = Form("EMCTRQA_histFastORL0Time");
  htitle = Form("EMCTRQA_histFastORL0Time;FastOR abs. ID;time");
  fHistManager.CreateTH2(hname, htitle, 5000, 0, 5000, 21, -1, 20);

  hname = Form("EMCTRQA_histFastORL0AmpVsTime");
  htitle = Form("EMCTRQA_histFastORL0AmpVsTime;time;amplitude");
  fHistManager.CreateTH2(hname, htitle, 21, -1, 20, 1024, 0, 4096);

  hname = Form("EMCTRQA_histFastORL0TimeOk");
  htitle = Form("EMCTRQA_histFastORL0TimeOk;FastOR abs. ID;entries (%d < time < %d)", fL0MinTime, fL0MaxTime);
  fHistManager.CreateTH1(hname, htitle, 5000, 0, 5000);

  hname = Form("EMCTRQA_histFastORL0AmpTimeOk");
  htitle = Form("EMCTRQA_histFastORL0AmpTimeOk;FastOR abs. ID;ADC counts  (%d < time < %d)", fL0MinTime, fL0MaxTime);
  fHistManager.CreateTH2(hname, htitle, 5000, 0, 5000, 1024, 0, 4096);

  hname = Form("EMCTRQA_histFastORL1");
  htitle = Form("EMCTRQA_histFastORL1;FastOR abs. ID;entries above 0");
  fHistManager.CreateTH1(hname, htitle, 5000, 0, 5000);

  hname = Form("EMCTRQA_histLargeAmpFastORL1");
  htitle = Form("EMCTRQA_histLargeAmpFastORL1 (>%d);FastOR abs. ID;entries above %d", fFastorL1Th, fFastorL1Th);
  fHistManager.CreateTH1(hname, htitle, 5000, 0, 5000);

  hname = Form("EMCTRQA_histFastORL1Amp");
  htitle = Form("EMCTRQA_histFastORL1Amp;FastOR abs. ID;L1 time sum");
  fHistManager.CreateTH2(hname, htitle, 5000, 0, 5000, 1024, 0, 4096);

  hname = Form("EMCTRQA_histFastORL1AmpVsL0Amp");
  htitle = Form("EMCTRQA_histFastORL1AmpVsL0Amp;L0 amplitude;L1 time sum;entries");
  fHistManager.CreateTH2(hname, htitle, 256, 0, 1024, 256, 0, 1024);

  hname = Form("EMCTRQA_histCellAmp");
  htitle = Form("EMCTRQA_histCellAmp;cell abs. ID;energy (GeV)");
  fHistManager.CreateTH2(hname, htitle, 20000, 0, 20000, 400, 0, 200);

  hname = Form("EMCTRQA_histCell");
  htitle = Form("EMCTRQA_histCell;cell abs. ID;entries above 0 GeV");
  fHistManager.CreateTH1(hname, htitle, 20000, 0, 20000);

  hname = Form("EMCTRQA_histLargeAmpCell");
  htitle = Form("EMCTRQA_histLargeAmpCell (>%.1f);cell abs. ID;entries above %.1f GeV", fFastorL1Th * EMCALTrigger::kEMCL1ADCtoGeV, fFastorL1Th * EMCALTrigger::kEMCL1ADCtoGeV);
  fHistManager.CreateTH1(hname, htitle, 20000, 0, 20000);

  const char* det[2] = { "EMCal", "DCal" };

  for (Int_t itrig = 0; itrig < 6; itrig++) {
    if (kEMCalTriggerNames[itrig].IsNull() || fEnabledTriggerTypes[itrig] == kFALSE) continue;
    for (Int_t itype = 0; itype < 3; itype++) {
      if (!fEnabledPatchTypes[itype]) continue;
      for (Int_t idet = 0; idet < 2; idet++) {
        hname = Form("EMCTRQA_hist%sPatchAmp%s%s", det[idet], kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
        htitle = Form("EMCTRQA_hist%sPatchAmp%s%s;amplitude;entries", det[idet], kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
        fHistManager.CreateTH1(hname, htitle, fgkMaxPatchAmp[itrig]/fADCperBin, 0, fgkMaxPatchAmp[itrig]);

        hname = Form("EMCTRQA_hist%sMaxPatchAmp%s%s", det[idet], kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
        htitle = Form("EMCTRQA_hist%sMaxPatchAmp%s%s;amplitude;entries", det[idet], kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
        fHistManager.CreateTH1(hname, htitle, fgkMaxPatchAmp[itrig]/fADCperBin, 0, fgkMaxPatchAmp[itrig]);

        for (Int_t itrig2 = itrig+1; itrig2 < 6; itrig2++) {
          if (kEMCalTriggerNames[itrig2].IsNull() || fEnabledTriggerTypes[itrig2] == kFALSE) continue;
          hname = Form("EMCTRQA_hist%s%sMaxVs%sMax%s", det[idet], kEMCalTriggerNames[itrig].Data(), kEMCalTriggerNames[itrig2].Data(), fgkPatchTypes[itype].Data());
          htitle = Form("EMCTRQA_hist%s%sMaxVs%sMax%s;%s max;%s max;entries", det[idet], kEMCalTriggerNames[itrig].Data(), kEMCalTriggerNames[itrig2].Data(), fgkPatchTypes[itype].Data(), kEMCalTriggerNames[itrig2].Data(), kEMCalTriggerNames[itrig].Data());
          fHistManager.CreateTH2(hname, htitle, fgkMaxPatchAmp[itrig2]/fADCperBin, 0, fgkMaxPatchAmp[itrig2], fgkMaxPatchAmp[itrig]/fADCperBin, 0, fgkMaxPatchAmp[itrig]);
        }
      }

      hname = Form("EMCTRQA_histEMCalMaxVsDCalMax%s%s", kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
      htitle = Form("EMCTRQA_histEMCalMaxVsDCalMax%s%s;EMCal max;DCal max;entries", kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
      fHistManager.CreateTH2(hname, htitle, fgkMaxPatchAmp[itrig]/fADCperBin/4, 0, fgkMaxPatchAmp[itrig], fgkMaxPatchAmp[itrig]/fADCperBin/4, 0, fgkMaxPatchAmp[itrig]);

      hname = Form("EMCTRQA_histEdgePos%s%s", kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
      htitle = Form("EMCTRQA_histEdgePos%s%s;col;row;entries", kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
      fHistManager.CreateTH2(hname, htitle, 48, 0, 48, 105, 0, 105);

      hname = Form("EMCTRQA_histGeoPos%s%s", kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
      htitle = Form("EMCTRQA_histGeoPos%s%s;#eta;#phi;entries", kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
      fHistManager.CreateTH2(hname, htitle, 60, -1, 1, 150, 0, TMath::TwoPi());

      hname = Form("EMCTRQA_histLargeAmpEdgePos%s%s", kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
      htitle = Form("EMCTRQA_histLargeAmpEdgePos%s%s (>700);col;row;entries above 700", kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
      fHistManager.CreateTH2(hname, htitle, 48, 0, 48, 105, 0, 105);
    }
  }
}

/**
 * Process a patch, filling relevant histograms.
 * \param patch Pointer to a valid trigger patch
 */
void AliEmcalTriggerQAPP::ProcessPatch(AliEMCALTriggerPatchInfo* patch)
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

  for (Int_t itrig = 0; itrig < 6; itrig++) {
    if (kEMCalTriggerNames[itrig].IsNull() || fEnabledTriggerTypes[itrig] == kFALSE) continue;

    for (Int_t itype = 0; itype < 3; itype++) {
      if (!fEnabledPatchTypes[itype]) continue;
      if (!patch->TestTriggerBit(triggerBits[itrig]+offsets[itype])) continue;

      TString det;

      if (patch->IsEMCal()) {
        det = "EMCal";
        if (fMaxPatchEMCal[itrig][itype] < amplitudes[itype]) fMaxPatchEMCal[itrig][itype] = amplitudes[itype];
      }
      else if (patch->IsDCalPHOS()) {
        det = "DCal";
        if (fMaxPatchDCal[itrig][itype] < amplitudes[itype]) fMaxPatchDCal[itrig][itype] = amplitudes[itype];
      }
      else {
        AliWarning(Form("Patch is not EMCal nor DCal/PHOS (pos: %d, %d)", patch->GetRowStart(), patch->GetColStart()));
      }

      hname = Form("EMCTRQA_histEdgePos%s%s", kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
      fHistManager.FillTH2(hname, patch->GetColStart(), patch->GetRowStart());


      hname = Form("EMCTRQA_histGeoPos%s%s", kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
      fHistManager.FillTH2(hname, patch->GetEtaGeo(), TVector2::Phi_0_2pi(patch->GetPhiGeo()));

      if (amplitudes[itype] > 700) {
        hname = Form("EMCTRQA_histLargeAmpEdgePos%s%s", kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
        fHistManager.FillTH2(hname, patch->GetColStart(), patch->GetRowStart());
      }

      hname = Form("EMCTRQA_hist%sPatchAmp%s%s", det.Data(), kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
      fHistManager.FillTH1(hname, amplitudes[itype]);
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
 * Process a cell, filling relevant histograms.
 * \param cell Const reference to a AliEmcalCellInfo object
 */
void AliEmcalTriggerQAPP::ProcessCell(const AliEmcalCellInfo& cell)
{
  TString hname;

  if (fOfflineBadChannels.find(cell.fAbsId) != fOfflineBadChannels.end()) return;

  if (cell.fEnergy > 0) {
    hname = Form("EMCTRQA_histCell");
    fHistManager.FillTH1(hname, cell.fAbsId);

    hname = Form("EMCTRQA_histCellAmp");
    fHistManager.FillTH2(hname, cell.fAbsId, cell.fEnergy);

    if (cell.fEnergy > fFastorL1Th * EMCALTrigger::kEMCL1ADCtoGeV) {
      hname = Form("EMCTRQA_histLargeAmpCell");
      fHistManager.FillTH1(hname, cell.fAbsId);
    }
  }
}

/**
 * Process a FastOR, filling relevant histograms.
 * \param patch Pointer to a valid trigger FastOR
 */
void AliEmcalTriggerQAPP::ProcessFastor(AliEMCALTriggerFastOR* fastor)
{
  TString hname;

  if (fastor->GetL0Amp() > 0) {
    hname = Form("EMCTRQA_histFastORL0");
    fHistManager.FillTH1(hname, fastor->GetAbsId());

    hname = Form("EMCTRQA_histFastORL0Amp");
    fHistManager.FillTH2(hname, fastor->GetAbsId(), fastor->GetL0Amp());

    hname = Form("EMCTRQA_histFastORL0Time");
    fHistManager.FillTH2(hname, fastor->GetAbsId(), fastor->GetL0Time());

    hname = Form("EMCTRQA_histFastORL0AmpVsTime");
    fHistManager.FillTH2(hname, fastor->GetL0Time(), fastor->GetL0Amp());

    if (fastor->GetL0Time() > fL0MinTime && fastor->GetL0Time() < fL0MaxTime) {
      hname = Form("EMCTRQA_histFastORL0TimeOk");
      fHistManager.FillTH1(hname, fastor->GetAbsId());

      hname = Form("EMCTRQA_histFastORL0AmpTimeOk");
      fHistManager.FillTH2(hname, fastor->GetAbsId(), fastor->GetL0Amp());
    }
  }

  if (fastor->GetL0Amp() > fFastorL0Th) {
    hname = Form("EMCTRQA_histLargeAmpFastORL0");
    fHistManager.FillTH1(hname, fastor->GetAbsId());
  }

  if (fastor->GetL1Amp() > 0) {
    hname = Form("EMCTRQA_histFastORL1");
    fHistManager.FillTH1(hname, fastor->GetAbsId());

    hname = Form("EMCTRQA_histFastORL1Amp");
    fHistManager.FillTH2(hname, fastor->GetAbsId(), fastor->GetL1Amp());
  }

  if (fastor->GetL1Amp() > fFastorL1Th) {
    hname = Form("EMCTRQA_histLargeAmpFastORL1");
    fHistManager.FillTH1(hname, fastor->GetAbsId());
  }

  if (fastor->GetL1Amp() > 0 && fastor->GetL0Amp() > 0) {
    hname = Form("EMCTRQA_histFastORL1AmpVsL0Amp");
    fHistManager.FillTH2(hname, fastor->GetL0Amp(), fastor->GetL1Amp());
  }
}


/**
 * This method should be called at the end of each event.
 */
void AliEmcalTriggerQAPP::EventCompleted()
{
  AliDebug(2, Form("Entering AliEmcalTriggerQAAP::EventCompleted"));

  TString hname;

  for (Int_t itype = 0; itype < 3; itype++) {
    if (!fEnabledPatchTypes[itype]) continue;

    AliDebug(2, Form("Patch type %s", fgkPatchTypes[itype].Data()));

    for (Int_t itrig = 0; itrig < 6; itrig++) {
      if (kEMCalTriggerNames[itrig].IsNull() || fEnabledTriggerTypes[itrig] == kFALSE) continue;

      AliDebug(2, Form("Trigger type: %s", kEMCalTriggerNames[itype].Data()));

      hname = Form("EMCTRQA_histEMCalMaxVsDCalMax%s%s", kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
      fHistManager.FillTH2(hname, fMaxPatchEMCal[itrig][itype], fMaxPatchDCal[itrig][itype]);

      hname = Form("EMCTRQA_histEMCalMaxPatchAmp%s%s", kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
      fHistManager.FillTH1(hname, fMaxPatchEMCal[itrig][itype]);

      hname = Form("EMCTRQA_histDCalMaxPatchAmp%s%s", kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
      fHistManager.FillTH1(hname, fMaxPatchDCal[itrig][itype]);

      for (Int_t itrig2 = itrig+1; itrig2 < 6; itrig2++) {
        if (kEMCalTriggerNames[itrig2].IsNull() || fEnabledTriggerTypes[itrig2] == kFALSE) continue;

        hname = Form("EMCTRQA_histEMCal%sMaxVs%sMax%s", kEMCalTriggerNames[itrig].Data(), kEMCalTriggerNames[itrig2].Data(), fgkPatchTypes[itype].Data());
        fHistManager.FillTH2(hname, fMaxPatchEMCal[itrig2][itype], fMaxPatchEMCal[itrig][itype]);

        hname = Form("EMCTRQA_histDCal%sMaxVs%sMax%s", kEMCalTriggerNames[itrig].Data(), kEMCalTriggerNames[itrig2].Data(), fgkPatchTypes[itype].Data());
        fHistManager.FillTH2(hname, fMaxPatchDCal[itrig2][itype], fMaxPatchDCal[itrig][itype]);
      }

      fMaxPatchEMCal[itrig][itype] = 0;
      fMaxPatchDCal[itrig][itype] = 0;
    }
  }
}
