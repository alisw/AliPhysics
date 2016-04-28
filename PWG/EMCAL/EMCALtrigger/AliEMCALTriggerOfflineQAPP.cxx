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

#include "AliEMCALTriggerOfflineQAPP.h"

/// \cond CLASSIMP
ClassImp(AliEMCALTriggerOfflineQAPP);
/// \endcond

/// Dummy constructor for ROOT I/O
AliEMCALTriggerOfflineQAPP::AliEMCALTriggerOfflineQAPP():
  AliEMCALTriggerQA(),
  fOfflineBadChannels(),
  fBadChannels(),
  fFastORPedestal(),
  fDCalPlots(kTRUE),
  fL0MinTime(7),
  fL0MaxTime(10),
  fMinCellAmp(0.),
  fMinL0FastORAmp(0),
  fMinL1FastORAmp(0),
  fHistManager(),
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
  fNL1DCal(0),
  fNTotTRU(0),
  fMaxFORabsId(0)
{
  for (Int_t i = 0; i < 3; i++) {
    for (Int_t itype = 0; itype < 6; itype++) {
      fMaxPatchEMCal[itype][i] = 0;
      fMaxPatchDCal[itype][i] = 0;
    }
  }

  EnablePatchType(kOnlinePatch, EMCALTrigger::kTMEMCalLevel0, kTRUE);
  EnablePatchType(kOnlinePatch, EMCALTrigger::kTMEMCalGammaL, kTRUE);
  EnablePatchType(kOnlinePatch, EMCALTrigger::kTMEMCalGammaH, kTRUE);
  EnablePatchType(kOnlinePatch, EMCALTrigger::kTMEMCalJetL, kTRUE);
  EnablePatchType(kOnlinePatch, EMCALTrigger::kTMEMCalJetH, kTRUE);

  EnablePatchType(kRecalcPatch, EMCALTrigger::kTMEMCalLevel0, kTRUE);
  EnablePatchType(kRecalcPatch, EMCALTrigger::kTMEMCalGammaH, kTRUE);
  EnablePatchType(kRecalcPatch, EMCALTrigger::kTMEMCalJetH, kTRUE);

  EnablePatchType(kOfflinePatch, EMCALTrigger::kTMEMCalLevel0, kTRUE);
  EnablePatchType(kOfflinePatch, EMCALTrigger::kTMEMCalGammaH, kTRUE);
  EnablePatchType(kOfflinePatch, EMCALTrigger::kTMEMCalJetH, kTRUE);
}

/// Default constructor
///
/// \param name Name of the object
AliEMCALTriggerOfflineQAPP::AliEMCALTriggerOfflineQAPP(const char* name):
  AliEMCALTriggerQA(name),
  fOfflineBadChannels(),
  fBadChannels(),
  fFastORPedestal(5000),
  fDCalPlots(kTRUE),
  fL0MinTime(7),
  fL0MaxTime(10),
  fMinCellAmp(0.),
  fMinL0FastORAmp(0),
  fMinL1FastORAmp(0),
  fHistManager(name),
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
  fNL1DCal(0),
  fNTotTRU(0),
  fMaxFORabsId(0)
{
  for (Int_t i = 0; i < 3; i++) {
    for (Int_t itype = 0; itype < 6; itype++) {
      fMaxPatchEMCal[itype][i] = 0;
      fMaxPatchDCal[itype][i] = 0;
    }
  }

  EnablePatchType(kOnlinePatch, EMCALTrigger::kTMEMCalLevel0, kTRUE);
  EnablePatchType(kOnlinePatch, EMCALTrigger::kTMEMCalGammaL, kTRUE);
  EnablePatchType(kOnlinePatch, EMCALTrigger::kTMEMCalGammaH, kTRUE);
  EnablePatchType(kOnlinePatch, EMCALTrigger::kTMEMCalJetL, kTRUE);
  EnablePatchType(kOnlinePatch, EMCALTrigger::kTMEMCalJetH, kTRUE);

  EnablePatchType(kRecalcPatch, EMCALTrigger::kTMEMCalLevel0, kTRUE);
  EnablePatchType(kRecalcPatch, EMCALTrigger::kTMEMCalGammaH, kTRUE);
  EnablePatchType(kRecalcPatch, EMCALTrigger::kTMEMCalJetH, kTRUE);

  EnablePatchType(kOfflinePatch, EMCALTrigger::kTMEMCalLevel0, kTRUE);
  EnablePatchType(kOfflinePatch, EMCALTrigger::kTMEMCalGammaH, kTRUE);
  EnablePatchType(kOfflinePatch, EMCALTrigger::kTMEMCalJetH, kTRUE);
}

/// Copy constructor
///
/// \param triggerQA Reference to an object to copy from
AliEMCALTriggerOfflineQAPP::AliEMCALTriggerOfflineQAPP(const AliEMCALTriggerOfflineQAPP& triggerQA) :
  AliEMCALTriggerQA(triggerQA),
  fOfflineBadChannels(triggerQA.fOfflineBadChannels),
  fBadChannels(),
  fFastORPedestal(triggerQA.fFastORPedestal),
  fDCalPlots(kTRUE),
  fL0MinTime(triggerQA.fL0MinTime),
  fL0MaxTime(triggerQA.fL0MaxTime),
  fMinCellAmp(triggerQA.fMinCellAmp),
  fMinL0FastORAmp(triggerQA.fMinL0FastORAmp),
  fMinL1FastORAmp(triggerQA.fMinL1FastORAmp),
  fHistManager(triggerQA.GetName()),
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
  fNL1DCal(0),
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
AliEMCALTriggerOfflineQAPP::~AliEMCALTriggerOfflineQAPP()
{
}

/// Read the offline bad channel map from a standard stream
///
/// \param stream A reference to a standard stream to read from (can be a file stream)
void AliEMCALTriggerOfflineQAPP::ReadOfflineBadChannelFromStream(std::istream& stream)
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
void AliEMCALTriggerOfflineQAPP::ReadOfflineBadChannelFromFile(const char* fname)
{
  std::ifstream file(fname);
  ReadOfflineBadChannelFromStream(file);
}

/// Read the FastOR bad channel map from a standard stream
///
/// \param stream A reference to a standard stream to read from (can be a file stream)
void AliEMCALTriggerOfflineQAPP::ReadFastORBadChannelFromStream(std::istream& stream)
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
void AliEMCALTriggerOfflineQAPP::ReadFastORBadChannelFromFile(const char* fname)
{
  std::ifstream file(fname);
  ReadFastORBadChannelFromStream(file);
}

/// Set the pedestal value for a FastOR
///
/// \param absId Absolute ID of a FastOR
/// \param ped   Pedestal value
void AliEMCALTriggerOfflineQAPP::SetFastORPedestal(Short_t absId, Float_t ped)
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
void AliEMCALTriggerOfflineQAPP::ReadFastORPedestalFromStream(std::istream& stream)
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
void AliEMCALTriggerOfflineQAPP::ReadFastORPedestalFromFile(const char* fname)
{
  std::ifstream file(fname);
  ReadFastORPedestalFromStream(file);
}

/// Initialize the class, i.e. allocate histograms.
void AliEMCALTriggerOfflineQAPP::Init()
{
  TString hname;
  TString htitle;

  // Geometry object not available at this stage, so need to hardcode some numbers
  fNTotTRU = 32; // there are 32 TRU in the EMCal (10 full SM * 3 TRU + 2 small SM * 1 TRU)
  Int_t ndet = 1;
  if (fDCalPlots) {
    fNTotTRU += 20; // there are additional 14 TRU in DCal, some skip in between (see below)
    ndet += 1;
  }
  fMaxFORabsId = fNTotTRU * 96;  // there are 96 channels in each TRU

  if (fTimeStampBinWidth > 0) fHistManager.CreateHistoGroup("ByTimeStamp");

  hname = Form("EMCTRQA_histFastORL0");
  htitle = Form("EMCTRQA_histFastORL0;FastOR abs. ID;entries above 0");
  fHistManager.CreateTH1(hname, htitle, fMaxFORabsId, 0, fMaxFORabsId);

  hname = Form("EMCTRQA_histLargeAmpFastORL0");
  htitle = Form("EMCTRQA_histLargeAmpFastORL0 (>%d);FastOR abs. ID;entries above %d", fFastorL0Th, fFastorL0Th);
  fHistManager.CreateTH1(hname, htitle, fMaxFORabsId, 0, fMaxFORabsId);

  hname = Form("EMCTRQA_histFastORL0Amp");
  htitle = Form("EMCTRQA_histFastORL0Amp;FastOR abs. ID;ADC counts");
  fHistManager.CreateTH2(hname, htitle, fMaxFORabsId, 0, fMaxFORabsId, 1024, 0, 4096);

  hname = Form("EMCTRQA_histFastORL0Time");
  htitle = Form("EMCTRQA_histFastORL0Time;FastOR abs. ID;time");
  fHistManager.CreateTH2(hname, htitle, fMaxFORabsId, 0, fMaxFORabsId, 21, -1, 20);

  hname = Form("EMCTRQA_histFastORL0AmpVsTime");
  htitle = Form("EMCTRQA_histFastORL0AmpVsTime;time;amplitude");
  fHistManager.CreateTH2(hname, htitle, 21, -1, 20, 1024, 0, 4096);

  hname = Form("EMCTRQA_histFastORL0TimeOk");
  htitle = Form("EMCTRQA_histFastORL0TimeOk;FastOR abs. ID;entries (%d < time < %d)", fL0MinTime, fL0MaxTime);
  fHistManager.CreateTH1(hname, htitle, fMaxFORabsId, 0, fMaxFORabsId);

  hname = Form("EMCTRQA_histFastORL0AmpTimeOk");
  htitle = Form("EMCTRQA_histFastORL0AmpTimeOk;FastOR abs. ID;ADC counts  (%d < time < %d)", fL0MinTime, fL0MaxTime);
  fHistManager.CreateTH2(hname, htitle, fMaxFORabsId, 0, fMaxFORabsId, 1024, 0, 4096);

  hname = Form("EMCTRQA_histFastORL1");
  htitle = Form("EMCTRQA_histFastORL1;FastOR abs. ID;entries above 0");
  fHistManager.CreateTH1(hname, htitle, fMaxFORabsId, 0, fMaxFORabsId);

  hname = Form("EMCTRQA_histLargeAmpFastORL1");
  htitle = Form("EMCTRQA_histLargeAmpFastORL1 (>%d);FastOR abs. ID;entries above %d", fFastorL1Th, fFastorL1Th);
  fHistManager.CreateTH1(hname, htitle, fMaxFORabsId, 0, fMaxFORabsId);

  hname = Form("EMCTRQA_histFastORL1Amp");
  htitle = Form("EMCTRQA_histFastORL1Amp;FastOR abs. ID;L1 time sum");
  fHistManager.CreateTH2(hname, htitle, fMaxFORabsId, 0, fMaxFORabsId, 1024, 0, 4096);

  hname = Form("EMCTRQA_histFastORL1AmpVsL0Amp");
  htitle = Form("EMCTRQA_histFastORL1AmpVsL0Amp;L0 amplitude;L1 time sum;entries");
  fHistManager.CreateTH2(hname, htitle, 256, 0, 1024, 256, 0, 1024);

  hname = Form("EMCTRQA_histCellAmp");
  htitle = Form("EMCTRQA_histCellAmp;cell abs. ID;energy (GeV)");
  fHistManager.CreateTH2(hname, htitle, fMaxFORabsId*4, 0, fMaxFORabsId*4, 400, 0, 200);

  hname = Form("EMCTRQA_histCell");
  htitle = Form("EMCTRQA_histCell;cell abs. ID;entries above 0 GeV");
  fHistManager.CreateTH1(hname, htitle, fMaxFORabsId*4, 0, fMaxFORabsId*4);

  hname = Form("EMCTRQA_histLargeAmpCell");
  htitle = Form("EMCTRQA_histLargeAmpCell (>%.1f);cell abs. ID;entries above %.1f GeV", fFastorL1Th * EMCALTrigger::kEMCL1ADCtoGeV, fFastorL1Th * EMCALTrigger::kEMCL1ADCtoGeV);
  fHistManager.CreateTH1(hname, htitle, fMaxFORabsId*4, 0, fMaxFORabsId*4);

  hname = Form("EMCTRQA_histCellAmpVsFastORL0Amp");
  htitle = Form("EMCTRQA_histCellAmpVsFastORL0Amp;FastOR L0 amplitude;2x2 cell sum energy (GeV)");
  fHistManager.CreateTH2(hname, htitle, 1024, 0, 4096, 400, 0, 200);

  hname = Form("EMCTRQA_histCellAmpVsFastORL0AmpTriggered");
  htitle = Form("EMCTRQA_histCellAmpVsFastORL0AmpTriggered;FastOR L0 amplitude;2x2 cell sum energy (GeV)");
  fHistManager.CreateTH2(hname, htitle, 1024, 0, 4096, 400, 0, 200);

  hname = Form("EMCTRQA_histCellAmpVsFastORL1Amp");
  htitle = Form("EMCTRQA_histCellAmpVsFastORL1Amp;FastOR L1 amplitude;2x2 cell sum energy (GeV)");
  fHistManager.CreateTH2(hname, htitle, 1024, 0, 4096, 400, 0, 200);

  for (Int_t nTRU = 0; nTRU < fNTotTRU; nTRU++) {
    if (nTRU == 34 || nTRU == 35 ||
        nTRU == 40 || nTRU == 41 ||
        nTRU == 46 || nTRU == 47) continue;
    hname = Form("EMCTRQA_histCellAmpVsFastORL0AmpTRU%d", nTRU);
    htitle = Form("EMCTRQA_histCellAmpVsFastORL0Amp%d;FastOR L0 amplitude;2x2 cell sum energy (GeV)", nTRU);
    fHistManager.CreateTH2(hname, htitle, 1024, 0, 4096, 400, 0, 200);

    hname = Form("EMCTRQA_histCellAmpVsFastORL0AmpTriggeredTRU%d", nTRU);
    htitle = Form("EMCTRQA_histCellAmpVsFastORL0AmpTriggeredTRU%d;FastOR L0 amplitude;2x2 cell sum energy (GeV)", nTRU);
    fHistManager.CreateTH2(hname, htitle, 1024, 0, 4096, 400, 0, 200);
  }

  hname = Form("EMCTRQA_histFastORNoOffline");
  htitle = Form("EMCTRQA_histFastORNoOffline;FastOR abs. ID;entries with no offline energy");
  fHistManager.CreateTH1(hname, htitle, fMaxFORabsId, 0, fMaxFORabsId);

  const char* det[2] = { "EMCal", "DCal" };

  for (Int_t itrig = 0; itrig < fgkNTriggerTypes; itrig++) {

    if (IsPatchTypeEnabled(kOfflinePatch, itrig) && IsPatchTypeEnabled(kRecalcPatch, itrig)) {
      for (Int_t idet = 0; idet < ndet; idet++) {
        hname = Form("EMCTRQA_hist%s%sMaxOfflineVsRecalc", det[idet], EMCALTrigger::kEMCalTriggerNames[itrig].Data());
        htitle = Form("EMCTRQA_hist%s%sMaxOfflineVsRecalc;Recalc;Offline;entries", det[idet], EMCALTrigger::kEMCalTriggerNames[itrig].Data());
        fHistManager.CreateTH2(hname, htitle, fgkMaxPatchAmp[itrig]/fADCperBin, 0, fgkMaxPatchAmp[itrig], fgkMaxPatchAmp[itrig]/fADCperBin, 0, fgkMaxPatchAmp[itrig]);

        hname = Form("EMCTRQA_hist%s%sOfflineVsRecalc", det[idet], EMCALTrigger::kEMCalTriggerNames[itrig].Data());
        htitle = Form("EMCTRQA_hist%s%sOfflineVsRecalc;Recalc;Offline;entries", det[idet], EMCALTrigger::kEMCalTriggerNames[itrig].Data());
        fHistManager.CreateTH2(hname, htitle, fgkMaxPatchAmp[itrig]/fADCperBin, 0, fgkMaxPatchAmp[itrig], fgkMaxPatchAmp[itrig]/fADCperBin, 0, fgkMaxPatchAmp[itrig]);
      }
    }

    for (Int_t itype = 0; itype < fgkNPatchTypes; itype++) {
      if (!IsPatchTypeEnabled(itype, itrig))  continue;
      for (Int_t idet = 0; idet < ndet; idet++) {
        hname = Form("EMCTRQA_hist%sPatchAmp%s%s", det[idet], EMCALTrigger::kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
        htitle = Form("EMCTRQA_hist%sPatchAmp%s%s;amplitude;entries", det[idet], EMCALTrigger::kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
        fHistManager.CreateTH1(hname, htitle, fgkMaxPatchAmp[itrig]/fADCperBin, 0, fgkMaxPatchAmp[itrig]);

        hname = Form("EMCTRQA_hist%sMaxPatchAmp%s%s", det[idet], EMCALTrigger::kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
        htitle = Form("EMCTRQA_hist%sMaxPatchAmp%s%s;amplitude;entries", det[idet], EMCALTrigger::kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
        fHistManager.CreateTH1(hname, htitle, fgkMaxPatchAmp[itrig]/fADCperBin, 0, fgkMaxPatchAmp[itrig]);
      }

      hname = Form("EMCTRQA_histMaxEdgePos%s%s", EMCALTrigger::kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
      htitle = Form("EMCTRQA_histMaxEdgePos%s%s;col;row;entries", EMCALTrigger::kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
      fHistManager.CreateTH2(hname, htitle, 48, 0, 48, 105, 0, 105);

      if (fDCalPlots) {
        hname = Form("EMCTRQA_histEMCal%sMaxVsDCal%sMax%s", EMCALTrigger::kEMCalTriggerNames[itrig].Data(), EMCALTrigger::kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
        htitle = Form("EMCTRQA_histEMCal%sMaxVsDCal%sMax%s;DCal %s max;EMCal %s max;entries", EMCALTrigger::kEMCalTriggerNames[itrig].Data(), EMCALTrigger::kEMCalTriggerNames[itrig].Data(), EMCALTrigger::kEMCalTriggerNames[itrig].Data(), EMCALTrigger::kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
        fHistManager.CreateTH2(hname, htitle, fgkMaxPatchAmp[itrig]/fADCperBin, 0, fgkMaxPatchAmp[itrig], fgkMaxPatchAmp[itrig]/fADCperBin, 0, fgkMaxPatchAmp[itrig]);
      }
    }
  }

  for (Int_t itype = 0; itype < fgkNPatchTypes; itype++) {
    for (Int_t itrig = 0; itrig < fgkNTriggerTypes; itrig++) {
      if (!IsPatchTypeEnabled(itype, itrig))  continue;
      for (Int_t itrig2 = itrig+1; itrig2 < 6; itrig2++) {
        if (!IsPatchTypeEnabled(itype, itrig2))  continue;

        hname = Form("EMCTRQA_histEMCal%sMaxVsEMCal%sMax%s", EMCALTrigger::kEMCalTriggerNames[itrig].Data(), EMCALTrigger::kEMCalTriggerNames[itrig2].Data(), fgkPatchTypes[itype].Data());
        htitle = Form("EMCTRQA_histEMCal%sMaxVsEMCal%sMax%s;EMCal %s max;EMCal %s max;entries", EMCALTrigger::kEMCalTriggerNames[itrig].Data(), EMCALTrigger::kEMCalTriggerNames[itrig2].Data(), fgkPatchTypes[itype].Data(), EMCALTrigger::kEMCalTriggerNames[itrig2].Data(), EMCALTrigger::kEMCalTriggerNames[itrig].Data());
        fHistManager.CreateTH2(hname, htitle, fgkMaxPatchAmp[itrig2]/fADCperBin, 0, fgkMaxPatchAmp[itrig2], fgkMaxPatchAmp[itrig]/fADCperBin, 0, fgkMaxPatchAmp[itrig]);


        if (fDCalPlots) {
          hname = Form("EMCTRQA_histDCal%sMaxVsDCal%sMax%s", EMCALTrigger::kEMCalTriggerNames[itrig].Data(), EMCALTrigger::kEMCalTriggerNames[itrig2].Data(), fgkPatchTypes[itype].Data());
          htitle = Form("EMCTRQA_histDCal%sMaxVsDCal%sMax%s;DCal %s max;DCal %s max;entries", EMCALTrigger::kEMCalTriggerNames[itrig].Data(), EMCALTrigger::kEMCalTriggerNames[itrig2].Data(), fgkPatchTypes[itype].Data(), EMCALTrigger::kEMCalTriggerNames[itrig2].Data(), EMCALTrigger::kEMCalTriggerNames[itrig].Data());
          fHistManager.CreateTH2(hname, htitle, fgkMaxPatchAmp[itrig2]/fADCperBin, 0, fgkMaxPatchAmp[itrig2], fgkMaxPatchAmp[itrig]/fADCperBin, 0, fgkMaxPatchAmp[itrig]);

          hname = Form("EMCTRQA_histEMCal%sMaxVsDCal%sMax%s", EMCALTrigger::kEMCalTriggerNames[itrig].Data(), EMCALTrigger::kEMCalTriggerNames[itrig2].Data(), fgkPatchTypes[itype].Data());
          htitle = Form("EMCTRQA_histEMCal%sMaxVsDCal%sMax%s;DCal %s max;EMCal %s max;entries", EMCALTrigger::kEMCalTriggerNames[itrig].Data(), EMCALTrigger::kEMCalTriggerNames[itrig2].Data(), fgkPatchTypes[itype].Data(), EMCALTrigger::kEMCalTriggerNames[itrig2].Data(), EMCALTrigger::kEMCalTriggerNames[itrig].Data());
          fHistManager.CreateTH2(hname, htitle, fgkMaxPatchAmp[itrig2]/fADCperBin, 0, fgkMaxPatchAmp[itrig2], fgkMaxPatchAmp[itrig]/fADCperBin, 0, fgkMaxPatchAmp[itrig]);

          hname = Form("EMCTRQA_histDCal%sMaxVsEMCal%sMax%s", EMCALTrigger::kEMCalTriggerNames[itrig].Data(), EMCALTrigger::kEMCalTriggerNames[itrig2].Data(), fgkPatchTypes[itype].Data());
          htitle = Form("EMCTRQA_histDCal%sMaxVsEMCal%sMax%s;EMCal %s max;DCal %s max;entries", EMCALTrigger::kEMCalTriggerNames[itrig].Data(), EMCALTrigger::kEMCalTriggerNames[itrig2].Data(), fgkPatchTypes[itype].Data(), EMCALTrigger::kEMCalTriggerNames[itrig2].Data(), EMCALTrigger::kEMCalTriggerNames[itrig].Data());
          fHistManager.CreateTH2(hname, htitle, fgkMaxPatchAmp[itrig2]/fADCperBin, 0, fgkMaxPatchAmp[itrig2], fgkMaxPatchAmp[itrig]/fADCperBin, 0, fgkMaxPatchAmp[itrig]);
        }
      }
    }
  }

  for (Int_t idet = 0; idet < ndet; idet++) {
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
void AliEMCALTriggerOfflineQAPP::ProcessPatch(const AliEMCALTriggerPatchInfo* patch)
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
    if (!fDCalPlots) return;
    det = "DCal";
  }
  else {
    AliWarning(Form("Patch is not EMCal nor DCal/PHOS (pos: %d, %d)", patch->GetRowStart(), patch->GetColStart()));
  }

  for (Int_t itrig = 0; itrig < fgkNTriggerTypes; itrig++) {
    if (IsPatchTypeEnabled(kRecalcPatch, itrig) && IsPatchTypeEnabled(kOfflinePatch, itrig)) {
      hname = Form("EMCTRQA_hist%s%sOfflineVsRecalc", det.Data(), EMCALTrigger::kEMCalTriggerNames[itrig].Data());
      fHistManager.FillTH2(hname, GetAmplitude(patch, kRecalcPatch) , GetAmplitude(patch, kOfflinePatch));
    }

    for (Int_t itype = 0; itype < fgkNPatchTypes; itype++) {
      if (!IsPatchTypeEnabled(itype, itrig))  continue;

      if (!(patch->*(triggerCheck[itype][itrig]))()) continue;

      if (patch->IsEMCal()) {
        if (GetAmplitude(fMaxPatchEMCal[itrig][itype], itype) < GetAmplitude(patch, itype)) fMaxPatchEMCal[itrig][itype] = patch;
      }
      else if (patch->IsDCalPHOS()) {
        if (GetAmplitude(fMaxPatchDCal[itrig][itype], itype) < GetAmplitude(patch, itype)) fMaxPatchDCal[itrig][itype] = patch;
      }

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
 * Process a cell, filling relevant histograms.
 * \param cell Const reference to a AliEmcalCellInfo object
 */
void AliEMCALTriggerOfflineQAPP::ProcessCell(const AliEMCALCellInfo& cell)
{
  TString hname;

  if (fGeom) {
    Double_t pos[3] = {0};
    fGeom->GetGlobal(cell.fAbsId, pos);
    if (fGeom->IsInDCAL(pos[0], pos[1], pos[2])) {
      if (!fDCalPlots) return;
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
}

/**
 * Process a FastOR, filling relevant histograms.
 * \param patch Pointer to a valid trigger FastOR
 */
void AliEMCALTriggerOfflineQAPP::ProcessFastor(const AliEMCALTriggerFastOR* fastor, AliVCaloCells* cells)
{
  TString hname;

  if (fBadChannels.find(fastor->GetAbsId()) != fBadChannels.end()) return;

  Int_t L0amp = fastor->GetL0Amp() - fFastORPedestal[fastor->GetAbsId()];
  Bool_t isDCal = kFALSE;

  Double_t offlineAmp = 0;
  Int_t nTRU = -1;
  Int_t nADC = -1;

  if (fGeom) {
    Int_t idx[4] = {-1};
    fGeom->GetCellIndexFromFastORIndex(fastor->GetAbsId(), idx);
    fGeom->GetTRUFromAbsFastORIndex(fastor->GetAbsId(), nTRU, nADC);
    Double_t pos[3] = {0};
    if (idx[0] >= 0) {
      fGeom->GetGlobal(idx[0], pos);
      isDCal = fGeom->IsInDCAL(pos[0], pos[1], pos[2]);
      if (isDCal && !fDCalPlots) return;

      if (L0amp > fMinL0FastORAmp) {
        if (cells) {
          for (Int_t i = 0; i < 4; i++) {
            offlineAmp += cells->GetCellAmplitude(idx[i]);
          }
        }
      }
    }
  }

  if (L0amp > fMinL0FastORAmp) {
    if (fTimeStampBinWidth > 0) {
      hname = TString::Format("ByTimeStamp/EMCTRQA_histFastORL0_%u_%u", fEventTimeStampBin, fEventTimeStampBin+fTimeStampBinWidth);
      fHistManager.FillTH1(hname, fastor->GetAbsId());
    }

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

      hname = Form("EMCTRQA_histCellAmpVsFastORL0AmpTriggered");
      fHistManager.FillTH2(hname, L0amp, offlineAmp);

      hname = Form("EMCTRQA_histCellAmpVsFastORL0AmpTriggeredTRU%d",nTRU);
      fHistManager.FillTH2(hname, L0amp, offlineAmp);
    }

    hname = Form("EMCTRQA_histCellAmpVsFastORL0Amp");
    fHistManager.FillTH2(hname, L0amp, offlineAmp);

    hname = Form("EMCTRQA_histCellAmpVsFastORL0AmpTRU%d",nTRU);
    fHistManager.FillTH2(hname, L0amp, offlineAmp);

    if (offlineAmp == 0) {
      hname = Form("EMCTRQA_histFastORNoOffline");
      fHistManager.FillTH1(hname, fastor->GetAbsId());
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
    if (fTimeStampBinWidth > 0) {
      hname = TString::Format("ByTimeStamp/EMCTRQA_histLargeAmpFastORL0_%u_%u", fEventTimeStampBin, fEventTimeStampBin+fTimeStampBinWidth);
      fHistManager.FillTH1(hname, fastor->GetAbsId());
    }

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

    hname = Form("EMCTRQA_histCellAmpVsFastORL1Amp");
    fHistManager.FillTH2(hname, fastor->GetL1Amp(), offlineAmp);
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
void AliEMCALTriggerOfflineQAPP::EventCompleted()
{
  AliDebug(2, Form("Entering AliEmcalTriggerQAAP::EventCompleted"));

  TString hname;

  for (Int_t itrig = 0; itrig < fgkNTriggerTypes; itrig++) {
    AliDebug(2, Form("Trigger type: %s", EMCALTrigger::kEMCalTriggerNames[itrig].Data()));

    if (IsPatchTypeEnabled(itrig, kOfflinePatch) && IsPatchTypeEnabled(itrig, kRecalcPatch)) {
      hname = Form("EMCTRQA_histEMCal%sMaxOfflineVsRecalc", EMCALTrigger::kEMCalTriggerNames[itrig].Data());
      fHistManager.FillTH2(hname, GetAmplitude(fMaxPatchEMCal[itrig][1], 1), GetAmplitude(fMaxPatchEMCal[itrig][2], 2));

      if (fDCalPlots) {
        hname = Form("EMCTRQA_histDCal%sMaxOfflineVsRecalc", EMCALTrigger::kEMCalTriggerNames[itrig].Data());
        fHistManager.FillTH2(hname, GetAmplitude(fMaxPatchDCal[itrig][1], 1), GetAmplitude(fMaxPatchDCal[itrig][2], 2));
      }
    }

    for (Int_t itype = 0; itype < fgkNPatchTypes; itype++) {
      if (!IsPatchTypeEnabled(itype, itrig))  continue;

      AliDebug(2, Form("Patch type %s", fgkPatchTypes[itype].Data()));

      if (fMaxPatchEMCal[itrig][itype]) {
        hname = Form("EMCTRQA_histMaxEdgePos%s%s", EMCALTrigger::kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
        fHistManager.FillTH2(hname, fMaxPatchEMCal[itrig][itype]->GetColStart(), fMaxPatchEMCal[itrig][itype]->GetRowStart());
      }

      hname = Form("EMCTRQA_histEMCalMaxPatchAmp%s%s", EMCALTrigger::kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
      fHistManager.FillTH1(hname, GetAmplitude(fMaxPatchEMCal[itrig][itype], itype));

      if (fDCalPlots) {
        if (fMaxPatchDCal[itrig][itype]) {
          hname = Form("EMCTRQA_histMaxEdgePos%s%s", EMCALTrigger::kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
          fHistManager.FillTH2(hname, fMaxPatchDCal[itrig][itype]->GetColStart(), fMaxPatchDCal[itrig][itype]->GetRowStart());
        }

        hname = Form("EMCTRQA_histDCalMaxPatchAmp%s%s", EMCALTrigger::kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
        fHistManager.FillTH1(hname, GetAmplitude(fMaxPatchDCal[itrig][itype], itype));

        hname = Form("EMCTRQA_histEMCal%sMaxVsDCal%sMax%s", EMCALTrigger::kEMCalTriggerNames[itrig].Data(), EMCALTrigger::kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[itype].Data());
        fHistManager.FillTH2(hname, GetAmplitude(fMaxPatchDCal[itrig][itype], itype), GetAmplitude(fMaxPatchEMCal[itrig][itype], itype));
      }
    }
  }

  for (Int_t itype = 0; itype < fgkNPatchTypes; itype++) {
    for (Int_t itrig = 0; itrig < fgkNTriggerTypes; itrig++) {
      if (!IsPatchTypeEnabled(itype, itrig))  continue;

      for (Int_t itrig2 = itrig+1; itrig2 < fgkNTriggerTypes; itrig2++) {
        if (!IsPatchTypeEnabled(itype, itrig2))  continue;

        hname = Form("EMCTRQA_histEMCal%sMaxVsEMCal%sMax%s", EMCALTrigger::kEMCalTriggerNames[itrig].Data(), EMCALTrigger::kEMCalTriggerNames[itrig2].Data(), fgkPatchTypes[itype].Data());
        fHistManager.FillTH2(hname, GetAmplitude(fMaxPatchEMCal[itrig2][itype], itype), GetAmplitude(fMaxPatchEMCal[itrig][itype], itype));

        if (fDCalPlots) {
          hname = Form("EMCTRQA_histDCal%sMaxVsDCal%sMax%s", EMCALTrigger::kEMCalTriggerNames[itrig].Data(), EMCALTrigger::kEMCalTriggerNames[itrig2].Data(), fgkPatchTypes[itype].Data());
          fHistManager.FillTH2(hname, GetAmplitude(fMaxPatchDCal[itrig2][itype], itype), GetAmplitude(fMaxPatchDCal[itrig][itype], itype));

          hname = Form("EMCTRQA_histEMCal%sMaxVsDCal%sMax%s", EMCALTrigger::kEMCalTriggerNames[itrig].Data(), EMCALTrigger::kEMCalTriggerNames[itrig2].Data(), fgkPatchTypes[itype].Data());
          fHistManager.FillTH2(hname, GetAmplitude(fMaxPatchDCal[itrig2][itype], itype), GetAmplitude(fMaxPatchEMCal[itrig][itype], itype));

          hname = Form("EMCTRQA_histDCal%sMaxVsEMCal%sMax%s", EMCALTrigger::kEMCalTriggerNames[itrig].Data(), EMCALTrigger::kEMCalTriggerNames[itrig2].Data(), fgkPatchTypes[itype].Data());
          fHistManager.FillTH2(hname, GetAmplitude(fMaxPatchEMCal[itrig2][itype], itype), GetAmplitude(fMaxPatchDCal[itrig][itype], itype));

        }
      }
    }
  }

  for (Int_t itype = 0; itype < fgkNPatchTypes; itype++) {
    for (Int_t itrig = 0; itrig < fgkNTriggerTypes; itrig++) {
      fMaxPatchEMCal[itrig][itype] = 0;
      fMaxPatchDCal[itrig][itype] = 0;
    }
  }

  fHistManager.FillTH2("EMCTRQA_histEMCalOfflineSumVsL0Sum", fSumL0EMCal, fSumOfflineEMCal);
  fHistManager.FillTH2("EMCTRQA_histEMCalOfflineSumVsL1Sum", fSumL1EMCal, fSumOfflineEMCal);
  fHistManager.FillTH2("EMCTRQA_histEMCalL1SumVsL0Sum", fSumL0EMCal, fSumL1EMCal);
  fHistManager.FillTH2("EMCTRQA_histEMCalNCellVsNL0", fNL0EMCal, fNCellEMCal);
  fHistManager.FillTH2("EMCTRQA_histEMCalNCellVsNL1", fNL1EMCal, fNCellEMCal);
  fHistManager.FillTH2("EMCTRQA_histEMCalNL1VsNL0", fNL0EMCal, fNL1EMCal);

  if (fDCalPlots) {
    fHistManager.FillTH2("EMCTRQA_histDCalOfflineSumVsL0Sum", fSumL0DCal, fSumOfflineDCal);
    fHistManager.FillTH2("EMCTRQA_histDCalOfflineSumVsL1Sum", fSumL1DCal, fSumOfflineDCal);
    fHistManager.FillTH2("EMCTRQA_histDCalL1SumVsL0Sum", fSumL0DCal, fSumL1DCal);
    fHistManager.FillTH2("EMCTRQA_histDCalNCellVsNL0", fNL0DCal, fNCellDCal);
    fHistManager.FillTH2("EMCTRQA_histDCalNCellVsNL1", fNL1DCal, fNCellDCal);
    fHistManager.FillTH2("EMCTRQA_histDCalNL1VsNL0", fNL0DCal, fNL1DCal);
  }

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

/**
 * This function should be called every event to set the new time stamp.
 * It sets the time stamp in the internal field, computes the time stamp bin
 * based on fTimeStampBinWidth and creates the "by-time-stamp" histograms.
 * \param timeStamp Time stamp of the event
 */
void AliEMCALTriggerOfflineQAPP::EventTimeStamp(UInt_t timeStamp)
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

      hname = TString::Format("ByTimeStamp/EMCTRQA_histLargeAmpFastORL0_%u_%u", fEventTimeStampBin, fEventTimeStampBin+fTimeStampBinWidth);
      htitle = TString::Format("EMCTRQA_histLargeAmpFastORL0 (>%d);FastOR abs. ID;entries above %d", fFastorL0Th, fFastorL0Th);
      fHistManager.CreateTH1(hname, htitle, fMaxFORabsId, 0, fMaxFORabsId);
    }

    hname = TString::Format("ByTimeStamp/EMCTRQA_histEvents_%u_%u", fEventTimeStampBin, fEventTimeStampBin+fTimeStampBinWidth);
    fHistManager.FillTH1(hname, 0.);
  }
}
