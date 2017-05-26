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

#include <AliEMCALTriggerPatchInfo.h>
#include <AliEMCALTriggerFastOR.h>
#include <AliLog.h>
#include <AliEMCALGeometry.h>
#include <AliVCaloCells.h>
#include <AliEMCALTriggerConstants.h>

#include "AliEMCALTriggerOfflineLightQAPP.h"

/// \cond CLASSIMP
ClassImp(AliEMCALTriggerOfflineLightQAPP);
/// \endcond

/// Dummy constructor for ROOT I/O
AliEMCALTriggerOfflineLightQAPP::AliEMCALTriggerOfflineLightQAPP():
  AliEMCALTriggerQA(),
  fBadChannels(),
  fL0MinTime(7),
  fL0MaxTime(10),
  fMinL0FastORAmp(0),
  fMinL1FastORAmp(0),
  fHistManager(),
  fNTotTRU(0),
  fMaxFORabsId(0)
{
  for (Int_t i = 0; i < 3; i++) {
    for (Int_t itype = 0; itype < 6; itype++) {
      fMaxPatchEMCal[itype][i] = 0;
      fMaxPatchDCal[itype][i] = 0;
    }
  }
}

/// Default constructor
///
/// \param name Name of the object
AliEMCALTriggerOfflineLightQAPP::AliEMCALTriggerOfflineLightQAPP(const char* name):
  AliEMCALTriggerQA(name),
  fBadChannels(),
  fL0MinTime(7),
  fL0MaxTime(10),
  fMinL0FastORAmp(0),
  fMinL1FastORAmp(0),
  fHistManager(name),
  fNTotTRU(0),
  fMaxFORabsId(0)
{
  for (Int_t i = 0; i < 3; i++) {
    for (Int_t itype = 0; itype < 6; itype++) {
      fMaxPatchEMCal[itype][i] = 0;
      fMaxPatchDCal[itype][i] = 0;
    }
  }
}

/// Copy constructor
///
/// \param triggerQA Reference to an object to copy from
AliEMCALTriggerOfflineLightQAPP::AliEMCALTriggerOfflineLightQAPP(const AliEMCALTriggerOfflineLightQAPP& triggerQA) :
  AliEMCALTriggerQA(triggerQA),
  fBadChannels(),
  fL0MinTime(triggerQA.fL0MinTime),
  fL0MaxTime(triggerQA.fL0MaxTime),
  fMinL0FastORAmp(triggerQA.fMinL0FastORAmp),
  fMinL1FastORAmp(triggerQA.fMinL1FastORAmp),
  fHistManager(triggerQA.GetName()),
  fNTotTRU(0),
  fMaxFORabsId(0)
{
  for (Int_t i = 0; i < 3; i++) {
    for (Int_t itype = 0; itype < 6; itype++) {
      fMaxPatchEMCal[itype][i] = 0;
      fMaxPatchDCal[itype][i] = 0;
    }
  }
}

/// Destructor
AliEMCALTriggerOfflineLightQAPP::~AliEMCALTriggerOfflineLightQAPP()
{
}

/// Read the FastOR bad channel map from a standard stream
///
/// \param stream A reference to a standard stream to read from (can be a file stream)
void AliEMCALTriggerOfflineLightQAPP::ReadFastORBadChannelFromStream(std::istream& stream)
{
  Short_t absId = -1;

  while (stream.good()) {
    stream >> absId;
    AddFastORBadChannel(absId);
  }
}

/// Read the FastOR bad channel map from a text file
///
/// \param fname Path and name of the file
void AliEMCALTriggerOfflineLightQAPP::ReadFastORBadChannelFromFile(const char* fname)
{
  std::ifstream file(fname);
  ReadFastORBadChannelFromStream(file);
}

/// Initialize the class, i.e. allocate histograms.
void AliEMCALTriggerOfflineLightQAPP::Init()
{
  TString hname;
  TString htitle;

  // Geometry object not available at this stage, so need to hardcode some numbers
  fNTotTRU = 32 + 20; // there are 32 TRU in the EMCal (10 full SM * 3 TRU + 2 small SM * 1 TRU)
  Int_t ndet = 2;
  fMaxFORabsId = fNTotTRU * 96;  // there are 96 channels in each TRU

  if (fTimeStampBinWidth > 0) fHistManager.CreateHistoGroup("ByTimeStamp");
  fHistManager.CreateHistoGroup("ByTRU");

  hname = "EMCTRQA_histFastORL0Hits";
  htitle = "EMCTRQA_histFastORL0Hits;col;row;entries";
  fHistManager.CreateTH2(hname, htitle, 48, 0, 48, 105, 0, 105);

  hname = "EMCTRQA_histFastORL0AccumulatedAmplitude";
  htitle = "EMCTRQA_histFastORL0AccumulatedAmplitude;col;row;accumulated amplitude";
  fHistManager.CreateTH2(hname, htitle, 48, 0, 48, 105, 0, 105);

  hname = "EMCTRQA_histFastORL0HitsTriggered";
  htitle = "EMCTRQA_histFastORL0HitsTriggered;col;row;entries";
  fHistManager.CreateTH2(hname, htitle, 48, 0, 48, 105, 0, 105);

  hname = "EMCTRQA_histFastORL0AccumulatedAmplitudeTriggered";
  htitle = "EMCTRQA_histFastORL0AccumulatedAmplitudeTriggered;col;row;accumulated amplitude";
  fHistManager.CreateTH2(hname, htitle, 48, 0, 48, 105, 0, 105);

  hname = "EMCTRQA_histFastORL0Time";
  htitle = "EMCTRQA_histFastORL0Time;FastOR abs. ID;time";
  fHistManager.CreateTH2(hname, htitle, fMaxFORabsId, 0, fMaxFORabsId, 21, -1, 20);

  hname = "EMCTRQA_histFastORL0AmpVsTime";
  htitle = "EMCTRQA_histFastORL0AmpVsTime;time;amplitude";
  fHistManager.CreateTH2(hname, htitle, 21, -1, 20, 1024/fADCperBin, 0, 1024);

  hname = "EMCTRQA_histFastORL1Hits";
  htitle = "EMCTRQA_histFastORL1Hits;col;row;entries";
  fHistManager.CreateTH2(hname, htitle, 48, 0, 48, 105, 0, 105);

  hname = "EMCTRQA_histFastORL1AccumulatedAmplitude";
  htitle = "EMCTRQA_histFastORL1AccumulatedAmplitude;col;row;accumulated amplitude";
  fHistManager.CreateTH2(hname, htitle, 48, 0, 48, 105, 0, 105);

  hname = "EMCTRQA_histFastORL1AmpVsL0Amp";
  htitle = "EMCTRQA_histFastORL1AmpVsL0Amp;L0 amplitude;L1 time sum;entries";
  fHistManager.CreateTH2(hname, htitle, 1024/fADCperBin, 0, 1024, 1024/fADCperBin, 0, 1024);

  hname = "EMCTRQA_histFastORL1AmpVsL0AmpTriggered";
  htitle = "EMCTRQA_histFastORL1AmpVsL0AmpTriggered;L0 amplitude;L1 time sum;entries";
  fHistManager.CreateTH2(hname, htitle, 1024/fADCperBin, 0, 1024, 1024/fADCperBin, 0, 1024);

  hname = "EMCTRQA_histCellAmpVsFastORL0Amp";
  htitle = "EMCTRQA_histCellAmpVsFastORL0Amp;FastOR L0 amplitude;2x2 cell sum energy (GeV)";
  fHistManager.CreateTH2(hname, htitle, 1024/fADCperBin, 0, 1024, 1024/fADCperBin, 0, 80);

  hname = "EMCTRQA_histCellAmpVsFastORL0AmpTriggered";
  htitle = "EMCTRQA_histCellAmpVsFastORL0AmpTriggered;FastOR L0 amplitude;2x2 cell sum energy (GeV)";
  fHistManager.CreateTH2(hname, htitle, 1024/fADCperBin, 0, 1024, 1024/fADCperBin, 0, 80);

  hname = "EMCTRQA_histCellAmpVsFastORL1Amp";
  htitle = "EMCTRQA_histCellAmpVsFastORL1Amp;FastOR L1 amplitude;2x2 cell sum energy (GeV)";
  fHistManager.CreateTH2(hname, htitle, 1024/fADCperBin, 0, 1024, 1024/fADCperBin, 0, 80);

  for (Int_t nTRU = 0; nTRU < fNTotTRU; nTRU++) {
    if (nTRU == 34 || nTRU == 35 ||
        nTRU == 40 || nTRU == 41 ||
        nTRU == 46 || nTRU == 47) continue;
    hname = TString::Format("ByTRU/EMCTRQA_histCellAmpVsFastORL0AmpTRU%d", nTRU);
    htitle = TString::Format("ByTRU/EMCTRQA_histCellAmpVsFastORL0Amp%d;FastOR L0 amplitude;2x2 cell sum energy (GeV)", nTRU);
    fHistManager.CreateTH2(hname, htitle, 1024/fADCperBin, 0, 1024, 1024/fADCperBin, 0, 80);

    hname = TString::Format("ByTRU/EMCTRQA_histFastORAmpSTUVsTRU%d", nTRU);
    htitle = TString::Format("ByTRU/EMCTRQA_histFastORAmpSTUVsTRU%d;TRU amplitude;STU amplitude", nTRU);
    fHistManager.CreateTH2(hname, htitle, 1024/fADCperBin, 0, 1024, 1024/fADCperBin, 0, 1024);

    hname = TString::Format("ByTRU/EMCTRQA_histCellAmpVsFastORL1AmpSTU%d", nTRU);
    htitle = TString::Format("ByTRU/EMCTRQA_histCellAmpVsFastORL0Amp%d;FastOR L1 amplitude;2x2 cell sum energy (GeV)", nTRU);
    fHistManager.CreateTH2(hname, htitle, 1024/fADCperBin, 0, 1024, 1024/fADCperBin, 0, 80);
  }

  const char* det[2] = { "EMCal", "DCal" };

  for (Int_t itrig = 0; itrig < fgkNTriggerTypes; itrig++) {

    for (Int_t itype = 0; itype < fgkNPatchTypes; itype++) {
      if (!IsPatchTypeEnabled(itype, itrig))  continue;
      for (Int_t idet = 0; idet < ndet; idet++) {
        hname = TString::Format("EMCTRQA_hist%sPatchAmp%s%s", det[idet], EMCALTrigger::kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
        htitle = TString::Format("EMCTRQA_hist%sPatchAmp%s%s;amplitude;entries", det[idet], EMCALTrigger::kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
        fHistManager.CreateTH1(hname, htitle, fgkMaxPatchAmp[itrig]/fADCperBin, 0, fgkMaxPatchAmp[itrig]);

        hname = TString::Format("EMCTRQA_hist%sMaxPatchAmp%s%s", det[idet], EMCALTrigger::kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
        htitle = TString::Format("EMCTRQA_hist%sMaxPatchAmp%s%s;amplitude;entries", det[idet], EMCALTrigger::kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
        fHistManager.CreateTH1(hname, htitle, fgkMaxPatchAmp[itrig]/fADCperBin, 0, fgkMaxPatchAmp[itrig]);
      }

      hname = TString::Format("EMCTRQA_histMaxEdgePos%s%s", EMCALTrigger::kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
      htitle = TString::Format("EMCTRQA_histMaxEdgePos%s%s;col;row;entries", EMCALTrigger::kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
      fHistManager.CreateTH2(hname, htitle, 48, 0, 48, 105, 0, 105);

      hname = TString::Format("EMCTRQA_histAccAmpEdgePos%s%s", EMCALTrigger::kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
      htitle = TString::Format("EMCTRQA_histAccAmpEdgePos%s%s;col;row;accumulated amplitude", EMCALTrigger::kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
      fHistManager.CreateTH2(hname, htitle, 48, 0, 48, 105, 0, 105);
    }
  }
}

/**
 * Process a patch, filling relevant histograms.
 * \param patch Pointer to a valid trigger patch
 */
void AliEMCALTriggerOfflineLightQAPP::ProcessPatch(const AliEMCALTriggerPatchInfo* patch)
{
  TString hname;

  Bool_t (AliEMCALTriggerPatchInfo::* triggerCheck[fgkNPatchTypes][fgkNTriggerTypes])(void) const = {
      { &AliEMCALTriggerPatchInfo::IsLevel0,
          &AliEMCALTriggerPatchInfo::IsGammaLow,
          &AliEMCALTriggerPatchInfo::IsGammaHigh,
          &AliEMCALTriggerPatchInfo::IsJetLow,
          &AliEMCALTriggerPatchInfo::IsJetHigh,
          &AliEMCALTriggerPatchInfo::IsBkg
      },
      { &AliEMCALTriggerPatchInfo::IsLevel0Recalc,
          &AliEMCALTriggerPatchInfo::IsGammaLowRecalc,
          &AliEMCALTriggerPatchInfo::IsGammaHighRecalc,
          &AliEMCALTriggerPatchInfo::IsJetLowRecalc,
          &AliEMCALTriggerPatchInfo::IsJetHighRecalc,
          &AliEMCALTriggerPatchInfo::IsBkgRecalc
      },
      { &AliEMCALTriggerPatchInfo::IsLevel0Simple,
          &AliEMCALTriggerPatchInfo::IsGammaLowSimple,
          &AliEMCALTriggerPatchInfo::IsGammaHighSimple,
          &AliEMCALTriggerPatchInfo::IsJetLowSimple,
          &AliEMCALTriggerPatchInfo::IsJetHighSimple,
          &AliEMCALTriggerPatchInfo::IsBkgSimple
      }
  };

  TString det;

  if (patch->IsEMCal()) {
    det = "EMCal";
  }
  else if (patch->IsDCalPHOS()) {
    det = "DCal";
  }
  else {
    AliWarning(Form("Patch is not EMCal nor DCal/PHOS (pos: %d, %d)", patch->GetRowStart(), patch->GetColStart()));
  }

  for (Int_t itrig = 0; itrig < fgkNTriggerTypes; itrig++) {
    for (Int_t itype = 0; itype < fgkNPatchTypes; itype++) {
      if (!IsPatchTypeEnabled(itype, itrig))  continue;

      if (!(patch->*(triggerCheck[itype][itrig]))()) continue;

      //if (itrig == EMCALTrigger::kTMEMCalLevel0 && itype == kRecalcPatch) Printf("Recalc amp is %d", GetAmplitude(patch, kRecalcPatch));

      if (patch->IsEMCal()) {
        if (GetAmplitude(fMaxPatchEMCal[itrig][itype], itype) < GetAmplitude(patch, itype)) fMaxPatchEMCal[itrig][itype] = patch;
      }
      else if (patch->IsDCalPHOS()) {
        if (GetAmplitude(fMaxPatchDCal[itrig][itype], itype) < GetAmplitude(patch, itype)) fMaxPatchDCal[itrig][itype] = patch;
      }

      hname = TString::Format("EMCTRQA_histAccAmpEdgePos%s%s", EMCALTrigger::kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
      fHistManager.FillTH2(hname, patch->GetColStart(), patch->GetRowStart(), GetAmplitude(patch, itype));

      hname = Form("EMCTRQA_hist%sPatchAmp%s%s", det.Data(), EMCALTrigger::kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
      fHistManager.FillTH1(hname, GetAmplitude(patch, itype));
    }
    if (fDebugLevel >= 2) {
      Printf("Type = %s; global pos = (%d, %d); Amp (online) = %d; Amp (offline) = %d; Patch energy = %.3f\n"
          "Position (CM): Eta=%.3f, Phi=%.3f\n"
          "Position (Geo): Eta=%.3f, Phi=%.3f\n",
          EMCALTrigger::kEMCalTriggerNames[itrig].Data(), patch->GetRowStart(), patch->GetColStart(), patch->GetADCAmp(), patch->GetADCOfflineAmp(), patch->GetPatchE(),
          patch->GetEtaCM(), patch->GetPhiCM(),
          patch->GetEtaGeo(), patch->GetPhiGeo());
    }
  }
}


/**
 * Process a FastOR, filling relevant histograms.
 * \param patch Pointer to a valid trigger FastOR
 */
void AliEMCALTriggerOfflineLightQAPP::ProcessFastor(const AliEMCALTriggerFastOR* fastor, AliVCaloCells* cells)
{
  TString hname;

  if (fBadChannels.find(fastor->GetAbsId()) != fBadChannels.end()) return;

  Bool_t isDCal = kFALSE;

  Double_t offlineAmp = 0;
  Int_t nTRU = -1;
  Int_t nADC = -1;
  Int_t iSM  = -1;
  Int_t iEta = -1;
  Int_t iPhi = -1;

  if (fGeom) {
    fGeom->GetPositionInSMFromAbsFastORIndex(fastor->GetAbsId(), iSM, iEta, iPhi);
    isDCal = fGeom->IsDCALSM(iSM);

    Int_t idx[4] = {-1};
    fGeom->GetCellIndexFromFastORIndex(fastor->GetAbsId(), idx);
    fGeom->GetTRUFromAbsFastORIndex(fastor->GetAbsId(), nTRU, nADC);

    if (nTRU >= fNTotTRU || nTRU < 0 ||
        nTRU == 34 || nTRU == 35 ||
        nTRU == 40 || nTRU == 41 ||
        nTRU == 46 || nTRU == 47) {
      if (fastor->GetL0Amp() > 0) AliError(Form("FastOR with abs ID %d (amp = %u) and TRU number %d: this TRU does not exist!!", fastor->GetAbsId(), fastor->GetL0Amp(), nTRU));
      return;
    }

    if (idx[0] >= 0) {
      if (fastor->GetL0Amp() > fMinL0FastORAmp) {
        if (cells) {
          for (Int_t i = 0; i < 4; i++) {
            offlineAmp += cells->GetCellAmplitude(idx[i]);
          }
        }
      }
    }
  }

  if (fastor->GetL0Amp() > fMinL0FastORAmp) {
    if (fTimeStampBinWidth > 0) {
      hname = TString::Format("ByTimeStamp/EMCTRQA_histFastORL0_%u_%u", fEventTimeStampBin, fEventTimeStampBin+fTimeStampBinWidth);
      fHistManager.FillTH1(hname, fastor->GetAbsId());
    }

    hname = "EMCTRQA_histFastORL0Hits";
    fHistManager.FillTH2(hname, fastor->GetGlobalCol(), fastor->GetGlobalRow());

    hname = "EMCTRQA_histFastORL0AccumulatedAmplitude";
    fHistManager.FillTH2(hname, fastor->GetGlobalCol(), fastor->GetGlobalRow(), fastor->GetL0Amp());

    hname = "EMCTRQA_histFastORL0Time";
    fHistManager.FillTH2(hname, fastor->GetAbsId(), fastor->GetL0Time());

    hname = "EMCTRQA_histFastORL0AmpVsTime";
    fHistManager.FillTH2(hname, fastor->GetL0Time(), fastor->GetL0Amp());

    hname = "EMCTRQA_histCellAmpVsFastORL0Amp";
    fHistManager.FillTH2(hname, fastor->GetL0Amp(), offlineAmp);

    if (fastor->GetL0Time() > fL0MinTime && fastor->GetL0Time() < fL0MaxTime) {
      hname = "EMCTRQA_histFastORL0HitsTriggered";
      fHistManager.FillTH2(hname, fastor->GetGlobalCol(), fastor->GetGlobalRow());

      hname = "EMCTRQA_histFastORL0AccumulatedAmplitudeTriggered";
      fHistManager.FillTH2(hname, fastor->GetGlobalCol(), fastor->GetGlobalRow(), fastor->GetL0Amp());

      hname = "EMCTRQA_histCellAmpVsFastORL0AmpTriggered";
      fHistManager.FillTH2(hname, fastor->GetL0Amp(), offlineAmp);
    }

    hname = TString::Format("ByTRU/EMCTRQA_histCellAmpVsFastORL0AmpTRU%d",nTRU);
    fHistManager.FillTH2(hname, fastor->GetL0Amp(), offlineAmp);
  }

  if (fastor->GetL1Amp() > fMinL1FastORAmp) {
    hname = "EMCTRQA_histFastORL1Hits";
    fHistManager.FillTH2(hname, fastor->GetGlobalCol(), fastor->GetGlobalRow());

    hname = "EMCTRQA_histFastORL1AccumulatedAmplitude";
    fHistManager.FillTH2(hname, fastor->GetGlobalCol(), fastor->GetGlobalRow(), fastor->GetL1Amp());

    hname = "EMCTRQA_histCellAmpVsFastORL1Amp";
    fHistManager.FillTH2(hname, fastor->GetL1Amp(), offlineAmp);

    hname = TString::Format("ByTRU/EMCTRQA_histCellAmpVsFastORL1AmpSTU%d",nTRU);
    fHistManager.FillTH2(hname, fastor->GetL1Amp(), offlineAmp);
  }

  UInt_t L1ampChopped = (fastor->GetL1Amp() >> 2) * 4;
  if (fastor->GetL0Amp() > fMinL0FastORAmp || L1ampChopped > fMinL1FastORAmp) {
    hname = TString::Format("ByTRU/EMCTRQA_histFastORAmpSTUVsTRU%d",nTRU);
    fHistManager.FillTH2(hname, fastor->GetL0Amp(), L1ampChopped);

    hname = "EMCTRQA_histFastORL1AmpVsL0Amp";
    fHistManager.FillTH2(hname, fastor->GetL0Amp(), L1ampChopped);

    if (fastor->GetL0Time() > fL0MinTime && fastor->GetL0Time() < fL0MaxTime) {
      hname = "EMCTRQA_histFastORL1AmpVsL0AmpTriggered";
      fHistManager.FillTH2(hname, fastor->GetL0Amp(), L1ampChopped);
    }
  }
}


/**
 * This method should be called at the end of each event.
 */
void AliEMCALTriggerOfflineLightQAPP::EventCompleted()
{
  AliDebug(2, Form("Entering AliEmcalTriggerQAAP::EventCompleted"));

  TString hname;

  for (Int_t itrig = 0; itrig < fgkNTriggerTypes; itrig++) {
    AliDebug(2, Form("Trigger type: %s", EMCALTrigger::kEMCalTriggerNames[itrig].Data()));

    for (Int_t itype = 0; itype < fgkNPatchTypes; itype++) {
      if (!IsPatchTypeEnabled(itype, itrig))  continue;

      AliDebug(2, Form("Patch type %s", fgkPatchTypes[itype].Data()));

      if (fMaxPatchEMCal[itrig][itype]) {
        hname = TString::Format("EMCTRQA_histMaxEdgePos%s%s", EMCALTrigger::kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
        fHistManager.FillTH2(hname, fMaxPatchEMCal[itrig][itype]->GetColStart(), fMaxPatchEMCal[itrig][itype]->GetRowStart());
      }

      hname = TString::Format("EMCTRQA_histEMCalMaxPatchAmp%s%s", EMCALTrigger::kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
      fHistManager.FillTH1(hname, GetAmplitude(fMaxPatchEMCal[itrig][itype], itype));

      if (fMaxPatchDCal[itrig][itype]) {
        hname = TString::Format("EMCTRQA_histMaxEdgePos%s%s", EMCALTrigger::kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
        fHistManager.FillTH2(hname, fMaxPatchDCal[itrig][itype]->GetColStart(), fMaxPatchDCal[itrig][itype]->GetRowStart());
      }

      hname = TString::Format("EMCTRQA_histDCalMaxPatchAmp%s%s", EMCALTrigger::kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
      fHistManager.FillTH1(hname, GetAmplitude(fMaxPatchDCal[itrig][itype], itype));
    }
  }

  for (Int_t itype = 0; itype < fgkNPatchTypes; itype++) {
    for (Int_t itrig = 0; itrig < fgkNTriggerTypes; itrig++) {
      fMaxPatchEMCal[itrig][itype] = 0;
      fMaxPatchDCal[itrig][itype] = 0;
    }
  }
}

/**
 * This function should be called every event to set the new time stamp.
 * It sets the time stamp in the internal field, computes the time stamp bin
 * based on fTimeStampBinWidth and creates the "by-time-stamp" histograms.
 * \param timeStamp Time stamp of the event
 */
void AliEMCALTriggerOfflineLightQAPP::EventTimeStamp(UInt_t timeStamp)
{
  AliEMCALTriggerQA::EventTimeStamp(timeStamp);

  if (fTimeStampBinWidth > 0) {
    TString hname = TString::Format("ByTimeStamp/EMCTRQA_histEvents_%u_%u", fEventTimeStampBin, fEventTimeStampBin+fTimeStampBinWidth);
    if (!fHistManager.FindObject(hname)) {
      TString htitle;

      hname = TString::Format("ByTimeStamp/EMCTRQA_histEvents_%u_%u", fEventTimeStampBin, fEventTimeStampBin+fTimeStampBinWidth);
      htitle = TString::Format("EMCTRQA_histEvents;;events");
      TH1* hevents = fHistManager.CreateTH1(hname, htitle, 1, 0, 1);
      hevents->GetXaxis()->SetBinLabel(1, TString::Format("%u <= time stamp < %u", fEventTimeStampBin, fEventTimeStampBin+fTimeStampBinWidth));

      hname = TString::Format("ByTimeStamp/EMCTRQA_histFastORL0_%u_%u", fEventTimeStampBin, fEventTimeStampBin+fTimeStampBinWidth);
      htitle = TString::Format("EMCTRQA_histFastORL0;FastOR abs. ID;entries above 0");
      fHistManager.CreateTH1(hname, htitle, fMaxFORabsId, 0, fMaxFORabsId);
    }

    hname = TString::Format("ByTimeStamp/EMCTRQA_histEvents_%u_%u", fEventTimeStampBin, fEventTimeStampBin+fTimeStampBinWidth);
    fHistManager.FillTH1(hname, 0.);
  }
}
