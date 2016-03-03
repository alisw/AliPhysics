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
#include "AliEMCALGeometry.h"
#include "AliVCaloCells.h"

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
  fOfflineBadChannels(),
  fFastORPedestal(),
  fFastorL0Th(400),
  fFastorL1Th(400),
  fADCperBin(16),
  fDebugLevel(0),
  fL0MinTime(7),
  fL0MaxTime(10),
  fMinCellAmp(0.),
  fMinL0FastORAmp(0),
  fMinL1FastORAmp(0),
  fHistManager(),
  fGeom(0),
  fSumOfflineEMCal(0),
  fSumL0EMCal(0),
  fSumL1EMCal(0),
  fSumOfflineDCal(0),
  fSumL0DCal(0),
  fSumL1DCal(0),
  fNCellEMCal(0),
  fNL0EMCal(0),
  fNL1EMCal(0),
  fNCellDCal(0),
  fNL0DCal(0),
  fNL1DCal(0)
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
  fOfflineBadChannels(),
  fFastORPedestal(5000),
  fFastorL0Th(400),
  fFastorL1Th(400),
  fADCperBin(16),
  fDebugLevel(0),
  fL0MinTime(7),
  fL0MaxTime(10),
  fMinCellAmp(0.),
  fMinL0FastORAmp(0),
  fMinL1FastORAmp(0),
  fHistManager(name),
  fGeom(0),
  fSumOfflineEMCal(0),
  fSumL0EMCal(0),
  fSumL1EMCal(0),
  fSumOfflineDCal(0),
  fSumL0DCal(0),
  fSumL1DCal(0),
  fNCellEMCal(0),
  fNL0EMCal(0),
  fNL1EMCal(0),
  fNCellDCal(0),
  fNL0DCal(0),
  fNL1DCal(0)
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
  fOfflineBadChannels(triggerQA.fOfflineBadChannels),
  fFastORPedestal(triggerQA.fFastORPedestal),
  fFastorL0Th(triggerQA.fFastorL0Th),
  fFastorL1Th(triggerQA.fFastorL1Th),
  fADCperBin(triggerQA.fADCperBin),
  fDebugLevel(triggerQA.fDebugLevel),
  fL0MinTime(triggerQA.fL0MinTime),
  fL0MaxTime(triggerQA.fL0MaxTime),
  fMinCellAmp(triggerQA.fMinCellAmp),
  fMinL0FastORAmp(triggerQA.fMinL0FastORAmp),
  fMinL1FastORAmp(triggerQA.fMinL1FastORAmp),
  fHistManager(triggerQA.GetName()),
  fGeom(0),
  fSumOfflineEMCal(0),
  fSumL0EMCal(0),
  fSumL1EMCal(0),
  fSumOfflineDCal(0),
  fSumL0DCal(0),
  fSumL1DCal(0),
  fNCellEMCal(0),
  fNL0EMCal(0),
  fNL1EMCal(0),
  fNCellDCal(0),
  fNL0DCal(0),
  fNL1DCal(0)
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

/// Set the pedestal value for a FastOR
///
/// \param absId Absolute ID of a FastOR
/// \param ped   Pedestal value
void AliEmcalTriggerQAPP::SetFastORPedestal(Short_t absId, Float_t ped)
{
  if (absId < 0 || absId >= fFastORPedestal.GetSize()) {
    AliWarning(Form("Abs. ID %d out of range (0,5000)", absId));
    return;
  }
  fFastORPedestal[absId] = ped;
}


/// Read the FastOR pedestals from a standard stream
///
/// \param stream A reference to a standard stream to read from (can be a file stream)
void AliEmcalTriggerQAPP::ReadFastORPedestalFromStream(std::istream& stream)
{
  Short_t absId = 0;
  Float_t ped = 0;
  while (stream.good()) {
    stream >> ped;
    SetFastORPedestal(absId, ped);
    absId++;
  }
}

/// Read the FastOR pedestals from a text file
///
/// \param fname Path and name of the file
void AliEmcalTriggerQAPP::ReadFastORPedestalFromFile(const char* fname)
{
  std::ifstream file(fname);
  ReadFastORPedestalFromStream(file);
}

/// Set the patch types to be plotted
///
/// \param type Patch type of which the status is being changed
/// \param e    Either enable or disable
void AliEmcalTriggerQAPP::EnablePatchType(PatchTypes_t type, Bool_t e)
{
  fEnabledPatchTypes[type] = e;
}

/// Set the trigger types to be plotted
///
/// \param type Trigger type of which the status is being changed
/// \param e    Either enable or disable
void AliEmcalTriggerQAPP::EnableTriggerType(EMCalTriggerType_t type, Bool_t e)
{
  fEnabledTriggerTypes[type] = e;
}

/// Actions to be executed only once for the first event
void AliEmcalTriggerQAPP::ExecOnce()
{
  fGeom = AliEMCALGeometry::GetInstance();
  if (!fGeom) {
    AliError("Could not get geometry!");
  }
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

  hname = Form("EMCTRQA_histCellAmpVsFastORL0Amp");
  htitle = Form("EMCTRQA_histCellAmpVsFastORL0Amp;FastOR amplitude;2x2 cell sum energy (GeV)");
  fHistManager.CreateTH2(hname, htitle, 1024, 0, 4096, 400, 0, 200);

  hname = Form("EMCTRQA_histFastORNoOffline");
  htitle = Form("EMCTRQA_histFastORNoOffline;FastOR abs. ID;entries with no offline energy");
  fHistManager.CreateTH1(hname, htitle, 5000, 0, 5000);

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
      fHistManager.CreateTH2(hname, htitle, fgkMaxPatchAmp[itrig]/fADCperBin, 0, fgkMaxPatchAmp[itrig], fgkMaxPatchAmp[itrig]/fADCperBin, 0, fgkMaxPatchAmp[itrig]);

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

  for (Int_t idet = 0; idet < 2; idet++) {
    hname = Form("EMCTRQA_hist%sOfflineSumVsL0Sum", det[idet]);
    htitle = hname + ";Sum of L0 amplitudes;Sum of cell energies (GeV)";
    fHistManager.CreateTH2(hname, htitle, 250, 0, 5000, 400, 0, 200);

    hname = Form("EMCTRQA_hist%sOfflineSumVsL1Sum", det[idet]);
    htitle = hname + ";Sum of L1 amplitudes;Sum of cell energies (GeV)";
    fHistManager.CreateTH2(hname, htitle, 250, 0, 5000, 400, 0, 200);

    hname = Form("EMCTRQA_hist%sL1SumVsL0Sum", det[idet]);
    htitle = hname + ";Sum of L0 amplitudes;Sum of L1 amplitudes";
    fHistManager.CreateTH2(hname, htitle, 250, 0, 5000, 250, 0, 5000);

    hname = Form("EMCTRQA_hist%sNCellVsNL0", det[idet]);
    htitle = hname + ";Number of L0 FastORs;Number of cells";
    fHistManager.CreateTH2(hname, htitle, 250, 0, 500, 250, 0, 2000);

    hname = Form("EMCTRQA_hist%sNCellVsNL1", det[idet]);
    htitle = hname + ";Number of L1 FastORs;Number of cells";
    fHistManager.CreateTH2(hname, htitle, 250, 0, 500, 250, 0, 2000);

    hname = Form("EMCTRQA_hist%sNL1VsNL0", det[idet]);
    htitle = hname + ";Number of L0 FastORs;Number of L1 FastORs";
    fHistManager.CreateTH2(hname, htitle, 250, 0, 500, 250, 0, 500);
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

  if (cell.fEnergy < fMinCellAmp) return;

  hname = Form("EMCTRQA_histCell");
  fHistManager.FillTH1(hname, cell.fAbsId);

  hname = Form("EMCTRQA_histCellAmp");
  fHistManager.FillTH2(hname, cell.fAbsId, cell.fEnergy);

  if (cell.fEnergy > fFastorL1Th * EMCALTrigger::kEMCL1ADCtoGeV) {
    hname = Form("EMCTRQA_histLargeAmpCell");
    fHistManager.FillTH1(hname, cell.fAbsId);
  }

  if (fGeom) {
    Double_t pos[3] = {0};
    fGeom->GetGlobal(cell.fAbsId, pos);
    if (fGeom->IsInDCAL(pos[0], pos[1], pos[2])) {
      fSumOfflineDCal += cell.fEnergy;
      fNCellDCal++;
    }
    else if (fGeom->IsInEMCAL(pos[0], pos[1], pos[2])) {
      fSumOfflineEMCal += cell.fEnergy;
      fNCellEMCal++;
    }
    else {
      AliWarning(Form("Cell with absolute ID %hd was not identified as neither DCal or EMCal!!", cell.fAbsId));
    }
  }
}

/**
 * Process a FastOR, filling relevant histograms.
 * \param patch Pointer to a valid trigger FastOR
 */
void AliEmcalTriggerQAPP::ProcessFastor(AliEMCALTriggerFastOR* fastor, AliVCaloCells* cells)
{
  TString hname;

  Int_t L0amp = fastor->GetL0Amp() - fFastORPedestal[fastor->GetAbsId()];
  Bool_t isDCal = kFALSE;

  if (fGeom) {
    Int_t idx[4] = {-1};
    fGeom->GetCellIndexFromFastORIndex(fastor->GetAbsId(), idx);
    Double_t pos[3] = {0};
    if (idx[0] >= 0) {
      fGeom->GetGlobal(idx[0], pos);
      isDCal = fGeom->IsInDCAL(pos[0], pos[1], pos[2]);

      if (L0amp > fMinL0FastORAmp) {
        if (cells) {
          Double_t offlineAmp = 0;
          for (Int_t i = 0; i < 4; i++) {
            offlineAmp += cells->GetCellAmplitude(idx[i]);
          }
          hname = Form("EMCTRQA_histCellAmpVsFastORL0Amp");
          fHistManager.FillTH2(hname, L0amp, offlineAmp);
          if (offlineAmp == 0) {
            hname = Form("EMCTRQA_histFastORNoOffline");
            fHistManager.FillTH1(hname, fastor->GetAbsId());
          }
        }
      }
    }
  }

  if (L0amp > fMinL0FastORAmp) {
    hname = Form("EMCTRQA_histFastORL0");
    fHistManager.FillTH1(hname, fastor->GetAbsId());

    hname = Form("EMCTRQA_histFastORL0Amp");
    fHistManager.FillTH2(hname, fastor->GetAbsId(), L0amp);

    hname = Form("EMCTRQA_histFastORL0Time");
    fHistManager.FillTH2(hname, fastor->GetAbsId(), fastor->GetL0Time());

    hname = Form("EMCTRQA_histFastORL0AmpVsTime");
    fHistManager.FillTH2(hname, fastor->GetL0Time(), L0amp);

    if (fastor->GetL0Time() > fL0MinTime && fastor->GetL0Time() < fL0MaxTime) {
      hname = Form("EMCTRQA_histFastORL0TimeOk");
      fHistManager.FillTH1(hname, fastor->GetAbsId());

      hname = Form("EMCTRQA_histFastORL0AmpTimeOk");
      fHistManager.FillTH2(hname, fastor->GetAbsId(), L0amp);
    }

    if (isDCal) {
      fSumL0DCal += L0amp;
      fNL0DCal++;
    }
    else {
      fSumL0EMCal += L0amp;
      fNL0EMCal++;
    }
  }

  if (L0amp > fFastorL0Th) {
    hname = Form("EMCTRQA_histLargeAmpFastORL0");
    fHistManager.FillTH1(hname, fastor->GetAbsId());
  }

  if (fastor->GetL1Amp() > fMinL1FastORAmp) {
    hname = Form("EMCTRQA_histFastORL1");
    fHistManager.FillTH1(hname, fastor->GetAbsId());

    hname = Form("EMCTRQA_histFastORL1Amp");
    fHistManager.FillTH2(hname, fastor->GetAbsId(), fastor->GetL1Amp());

    if (isDCal) {
      fSumL1DCal += fastor->GetL1Amp();
      fNL1DCal++;
    }
    else {
      fSumL1EMCal += fastor->GetL1Amp();
      fNL1EMCal++;
    }
  }

  if (fastor->GetL1Amp() > fFastorL1Th) {
    hname = Form("EMCTRQA_histLargeAmpFastORL1");
    fHistManager.FillTH1(hname, fastor->GetAbsId());
  }

  if (fastor->GetL1Amp() > fMinL1FastORAmp || L0amp > fMinL0FastORAmp) {
    hname = Form("EMCTRQA_histFastORL1AmpVsL0Amp");
    fHistManager.FillTH2(hname, L0amp, fastor->GetL1Amp());
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

  fHistManager.FillTH2("EMCTRQA_histDCalOfflineSumVsL0Sum", fSumL0DCal, fSumOfflineDCal);
  fHistManager.FillTH2("EMCTRQA_histDCalOfflineSumVsL1Sum", fSumL1DCal, fSumOfflineDCal);
  fHistManager.FillTH2("EMCTRQA_histDCalL1SumVsL0Sum", fSumL0DCal, fSumL1DCal);
  fHistManager.FillTH2("EMCTRQA_histEMCalOfflineSumVsL0Sum", fSumL0EMCal, fSumOfflineEMCal);
  fHistManager.FillTH2("EMCTRQA_histEMCalOfflineSumVsL1Sum", fSumL1EMCal, fSumOfflineEMCal);
  fHistManager.FillTH2("EMCTRQA_histEMCalL1SumVsL0Sum", fSumL0EMCal, fSumL1EMCal);


  fHistManager.FillTH2("EMCTRQA_histDCalNCellVsNL0", fNL0DCal, fNCellDCal);
  fHistManager.FillTH2("EMCTRQA_histDCalNCellVsNL1", fNL1DCal, fNCellDCal);
  fHistManager.FillTH2("EMCTRQA_histDCalNL1VsNL0", fNL0DCal, fNL1DCal);
  fHistManager.FillTH2("EMCTRQA_histEMCalNCellVsNL0", fNL0EMCal, fNCellEMCal);
  fHistManager.FillTH2("EMCTRQA_histEMCalNCellVsNL1", fNL1EMCal, fNCellEMCal);
  fHistManager.FillTH2("EMCTRQA_histEMCalNL1VsNL0", fNL0EMCal, fNL1EMCal);

  fSumOfflineEMCal = 0;
  fSumL0EMCal = 0;
  fSumL1EMCal = 0;
  fSumOfflineDCal = 0;
  fSumL0DCal = 0;
  fSumL1DCal = 0;

  fNCellEMCal = 0;
  fNL0EMCal = 0;
  fNL1EMCal = 0;
  fNCellDCal = 0;
  fNL0DCal = 0;
  fNL1DCal = 0;
}
