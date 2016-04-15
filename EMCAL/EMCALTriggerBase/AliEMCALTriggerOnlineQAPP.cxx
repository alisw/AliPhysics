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

#include <TH1F.h>
#include <TH2F.h>

#include "AliEMCALTriggerPatchInfo.h"
#include "AliEMCALTriggerFastOR.h"
#include "AliLog.h"
#include "AliEMCALGeometry.h"
#include "AliVCaloCells.h"
#include "AliEMCALTriggerConstants.h"
#include "AliEMCALTriggerOnlineQAPP.h"

/// \cond CLASSIMP
ClassImp(AliEMCALTriggerOnlineQAPP);
/// \endcond

/// Dummy constructor for ROOT I/O
AliEMCALTriggerOnlineQAPP::AliEMCALTriggerOnlineQAPP():
  AliEMCALTriggerQA(),
  fL0MinTime(7),
  fL0MaxTime(10),
  fMinCellAmp(0.),
  fMinL0FastORAmp(0),
  fMinL1FastORAmp(0),
  fHistograms(),
  fHistEvents(0),
  fHistFastORL0(0),
  fHistFastORL0LargeAmp(0),
  fHistFastORL0Amp(0),
  fHistFastORL0Time(0),
  fHistFastORL1(0),
  fHistFastORL1LargeAmp(0),
  fHistFastORL1Amp(0)
{
  for (Int_t iSM = 0; iSM < fgkSM; iSM++) {
    fHistFastORL0BySM[iSM] = 0;
    fHistFastORL0LargeAmpBySM[iSM] = 0;
    fHistFastORL0AmpBySM[iSM] = 0;
    fHistFEEvsTRUBySM[iSM] = 0;

    fHistFastORL1BySM[iSM] = 0;
    fHistFastORL1LargeAmpBySM[iSM] = 0;
    fHistFastORL1AmpBySM[iSM] = 0;
    fHistFEEvsSTUBySM[iSM] = 0;
  }

  for (Int_t ipatch = 0; ipatch < fgkNPatchTypes; ipatch++) {
    fEnabledPatchTypes[ipatch] = kTRUE;
  }

  for (Int_t itrigger = 0; itrigger < fgkNTriggerTypes; itrigger++) {
    for (Int_t ipatch = 0; ipatch < fgkNPatchTypes; ipatch++) {
      for (Int_t idet = 0; idet < fgkNDet; idet++) {
        fNPatches[idet][itrigger][ipatch] = 0;

        fHistNPatches[idet][itrigger][ipatch] = 0;
        fHistPatchAmp[idet][itrigger][ipatch] = 0;
        fHistMaxPatchAmp[idet][itrigger][ipatch] = 0;
      }
      fMaxPatchEMCal[itrigger][ipatch] = new AliEMCALTriggerPatchInfo;
      fMaxPatchDCal[itrigger][ipatch] = new AliEMCALTriggerPatchInfo;

      fHistMaxEdgePos[itrigger][ipatch] = 0;
      fHistAmpEdgePos[itrigger][ipatch] = 0;
    }
  }

  fEnabledPatchTypes[kOnlinePatch] = kTRUE;
  fEnabledPatchTypes[kRecalcPatch] = kTRUE;
  fEnabledPatchTypes[kOfflinePatch] = kTRUE;
}

/// Default constructor
///
/// \param name Name of the object
AliEMCALTriggerOnlineQAPP::AliEMCALTriggerOnlineQAPP(const char* name):
  AliEMCALTriggerQA(name),
  fL0MinTime(7),
  fL0MaxTime(10),
  fMinCellAmp(0.),
  fMinL0FastORAmp(0),
  fMinL1FastORAmp(0),
  fHistograms(),
  fHistEvents(0),
  fHistFastORL0(0),
  fHistFastORL0LargeAmp(0),
  fHistFastORL0Amp(0),
  fHistFastORL0Time(0),
  fHistFastORL1(0),
  fHistFastORL1LargeAmp(0),
  fHistFastORL1Amp(0)
{
  fHistograms.SetName(name);

  for (Int_t iSM = 0; iSM < fgkSM; iSM++) {
    fHistFastORL0BySM[iSM] = 0;
    fHistFastORL0LargeAmpBySM[iSM] = 0;
    fHistFastORL0AmpBySM[iSM] = 0;
    fHistFEEvsTRUBySM[iSM] = 0;

    fHistFastORL1BySM[iSM] = 0;
    fHistFastORL1LargeAmpBySM[iSM] = 0;
    fHistFastORL1AmpBySM[iSM] = 0;
    fHistFEEvsSTUBySM[iSM] = 0;
  }

  for (Int_t itrigger = 0; itrigger < fgkNTriggerTypes; itrigger++) {
    for (Int_t ipatch = 0; ipatch < fgkNPatchTypes; ipatch++) {
      for (Int_t idet = 0; idet < fgkNDet; idet++) {
        fNPatches[idet][itrigger][ipatch] = 0;

        fHistNPatches[idet][itrigger][ipatch] = 0;
        fHistPatchAmp[idet][itrigger][ipatch] = 0;
        fHistMaxPatchAmp[idet][itrigger][ipatch] = 0;
      }
      fMaxPatchEMCal[itrigger][ipatch] = new AliEMCALTriggerPatchInfo;
      fMaxPatchDCal[itrigger][ipatch] = new AliEMCALTriggerPatchInfo;

      fHistMaxEdgePos[itrigger][ipatch] = 0;
      fHistAmpEdgePos[itrigger][ipatch] = 0;
    }
  }

  fEnabledPatchTypes[kOnlinePatch] = kTRUE;
  fEnabledPatchTypes[kRecalcPatch] = kTRUE;
  fEnabledPatchTypes[kOfflinePatch] = kTRUE;

  fEnabledTriggerTypes[EMCALTrigger::kTMEMCalBkg] = kFALSE;
  fEnabledTriggerTypes[EMCALTrigger::kTMEMCalLevel0] = kTRUE;
  fEnabledTriggerTypes[EMCALTrigger::kTMEMCalGammaL] = kTRUE;
  fEnabledTriggerTypes[EMCALTrigger::kTMEMCalJetL] = kTRUE;
  fEnabledTriggerTypes[EMCALTrigger::kTMEMCalGammaH] = kTRUE;
  fEnabledTriggerTypes[EMCALTrigger::kTMEMCalJetH] = kTRUE;
}

/// Copy constructor
///
/// \param triggerQA Reference to an object to copy from
AliEMCALTriggerOnlineQAPP::AliEMCALTriggerOnlineQAPP(const AliEMCALTriggerOnlineQAPP& triggerQA) :
  AliEMCALTriggerQA(triggerQA),
  fL0MinTime(triggerQA.fL0MinTime),
  fL0MaxTime(triggerQA.fL0MaxTime),
  fMinCellAmp(triggerQA.fMinCellAmp),
  fMinL0FastORAmp(triggerQA.fMinL0FastORAmp),
  fMinL1FastORAmp(triggerQA.fMinL1FastORAmp),
  fHistograms(),
  fHistEvents(0),
  fHistFastORL0(0),
  fHistFastORL0LargeAmp(0),
  fHistFastORL0Amp(0),
  fHistFastORL0Time(0),
  fHistFastORL1(0),
  fHistFastORL1LargeAmp(0),
  fHistFastORL1Amp(0)
{
  fHistograms.SetName(triggerQA.GetName());

  for (Int_t iSM = 0; iSM < fgkSM; iSM++) {
    fHistFastORL0BySM[iSM] = 0;
    fHistFastORL0LargeAmpBySM[iSM] = 0;
    fHistFastORL0AmpBySM[iSM] = 0;
    fHistFEEvsTRUBySM[iSM] = 0;

    fHistFastORL1BySM[iSM] = 0;
    fHistFastORL1LargeAmpBySM[iSM] = 0;
    fHistFastORL1AmpBySM[iSM] = 0;
    fHistFEEvsSTUBySM[iSM] = 0;
  }

  for (Int_t ipatch = 0; ipatch < fgkNPatchTypes; ipatch++) {
    fEnabledPatchTypes[ipatch] = triggerQA.fEnabledPatchTypes[ipatch];
  }

  for (Int_t itrigger = 0; itrigger < fgkNTriggerTypes; itrigger++) {
    fEnabledTriggerTypes[itrigger] = triggerQA.fEnabledTriggerTypes[itrigger];

    for (Int_t ipatch = 0; ipatch < fgkNPatchTypes; ipatch++) {
      for (Int_t idet = 0; idet < fgkNDet; idet++) {
        fNPatches[idet][itrigger][ipatch] = 0;

        fHistNPatches[idet][itrigger][ipatch] = 0;
        fHistPatchAmp[idet][itrigger][ipatch] = 0;
        fHistMaxPatchAmp[idet][itrigger][ipatch] = 0;
      }
      fMaxPatchEMCal[itrigger][ipatch] = new AliEMCALTriggerPatchInfo;
      fMaxPatchDCal[itrigger][ipatch] = new AliEMCALTriggerPatchInfo;

      fHistMaxEdgePos[itrigger][ipatch] = 0;
      fHistAmpEdgePos[itrigger][ipatch] = 0;
    }
  }
}

/// Destructor
AliEMCALTriggerOnlineQAPP::~AliEMCALTriggerOnlineQAPP()
{
  for (Int_t itrigger = 0; itrigger < fgkNTriggerTypes; itrigger++) {
    for (Int_t ipatch = 0; ipatch < fgkNPatchTypes; ipatch++) {
      if (fMaxPatchEMCal[itrigger][ipatch]) delete fMaxPatchEMCal[itrigger][ipatch];
      if (fMaxPatchDCal[itrigger][ipatch]) delete fMaxPatchDCal[itrigger][ipatch];
    }
  }
}

/// Initialize the class, i.e. allocate histograms.
void AliEMCALTriggerOnlineQAPP::Init()
{
  TString hname;
  TString htitle;

  hname = "EMCTRQA_histEvents";
  htitle = "events;;total events";
  fHistEvents = new TH1F(hname, htitle, 1, 0, 1);
  fHistograms.Add(fHistEvents);

  hname = "EMCTRQA_histFastORL0";
  htitle = "L0;FastOR abs. ID;entries above 0";
  fHistFastORL0 = new TH1F(hname, htitle, 5000, 0, 5000);
  fHistograms.Add(fHistFastORL0);

  hname = "EMCTRQA_histFastORL0LargeAmp";
  htitle = TString::Format("L0 (amp>%d);FastOR abs. ID;entries above %d", fFastorL0Th, fFastorL0Th);
  fHistFastORL0LargeAmp = new TH1F(hname, htitle, 5000, 0, 5000);
  fHistograms.Add(fHistFastORL0LargeAmp);

  hname = "EMCTRQA_histFastORL0Amp";
  htitle = "L0 amplitudes;FastOR abs. ID;amplitude";
  fHistFastORL0Amp = new TH2F(hname, htitle, 5000, 0, 5000, 512, 0, 4096);
  fHistograms.Add(fHistFastORL0Amp);

  hname = "EMCTRQA_histFastORL0Time";
  htitle = "L0 trigger time;FastOR abs. ID;L0 trigger time";
  fHistFastORL0Time = new TH2F(hname, htitle, 5000, 0, 5000, 20, 0, 20);
  fHistograms.Add(fHistFastORL0Time);

  hname = "EMCTRQA_histFastORL1";
  htitle = "L1;FastOR abs. ID;entries above 0";
  fHistFastORL1 = new TH1F(hname, htitle, 5000, 0, 5000);
  fHistograms.Add(fHistFastORL1);

  hname = "EMCTRQA_histFastORL1LargeAmp";
  htitle = TString::Format("L1 (amp>%d);FastOR abs. ID;entries above %d", fFastorL1Th, fFastorL1Th);
  fHistFastORL1LargeAmp = new TH1F(hname, htitle, 5000, 0, 5000);
  fHistograms.Add(fHistFastORL1LargeAmp);

  hname = "EMCTRQA_histFastORL1Amp";
  htitle = "L1 amplitudes;FastOR abs. ID;amplitude";
  fHistFastORL1Amp = new TH2F(hname, htitle, 5000, 0, 5000, 512, 0, 4096);
  fHistograms.Add(fHistFastORL1Amp);

  Int_t nSM = fgkSM;
  if (fGeom) nSM = fGeom->GetNumberOfSuperModules();

  for (Int_t iSM = 0; iSM < nSM; iSM++) {
    hname = TString::Format("EMCTRQA_histFastORL0_SM%d", iSM);
    htitle = TString::Format("SM%d L0 (amp>0);eta id;phi id;entries above 0", iSM);
    fHistFastORL0BySM[iSM] = new TH2F(hname, htitle, 24, 0, 24, 12, 0, 12);
    fHistograms.Add(fHistFastORL0BySM[iSM]);

    hname = TString::Format("EMCTRQA_histFastORL0LargeAmp_SM%d", iSM);
    htitle = TString::Format("SM%d L0 (amp>%d);eta id;phi id;entries above %d", iSM, fFastorL0Th, fFastorL0Th);
    fHistFastORL0LargeAmpBySM[iSM] = new TH2F(hname, htitle, 24, 0, 24, 12, 0, 12);
    fHistograms.Add(fHistFastORL0LargeAmpBySM[iSM]);

    hname = TString::Format("EMCTRQA_histFastORL0Amp_SM%d", iSM);
    htitle = TString::Format("SM%d L0 int. amplitude;eta id;phi id;amplitude", iSM);
    fHistFastORL0AmpBySM[iSM] = new TH2F(hname, htitle, 24, 0, 24, 12, 0, 12);
    fHistograms.Add(fHistFastORL0AmpBySM[iSM]);

    hname = TString::Format("EMCTRQA_histFEEvsTRU_SM%d", iSM);
    htitle = TString::Format("SM%d FEE vs TRU;TRU amplitude;FEE energy (GeV)", iSM);
    fHistFEEvsTRUBySM[iSM] = new TH2F(hname, htitle, 128, 0, 2048, 150, 0, 150);
    fHistograms.Add(fHistFEEvsTRUBySM[iSM]);

    hname = TString::Format("EMCTRQA_histFastORL1_SM%d", iSM);
    htitle = TString::Format("SM%d L1 (amp>0);eta id;phi id;entries above 0", iSM);
    fHistFastORL1BySM[iSM] = new TH2F(hname, htitle, 24, 0, 24, 12, 0, 12);
    fHistograms.Add(fHistFastORL1BySM[iSM]);

    hname = TString::Format("EMCTRQA_histFastORL1LargeAmp_SM%d", iSM);
    htitle = TString::Format("SM%d L1 (amp>%d);eta id;phi id;entries above %d", iSM, fFastorL1Th, fFastorL1Th);
    fHistFastORL1LargeAmpBySM[iSM] = new TH2F(hname, htitle, 24, 0, 24, 12, 0, 12);
    fHistograms.Add(fHistFastORL1LargeAmpBySM[iSM]);

    hname = TString::Format("EMCTRQA_histFastORL1Amp_SM%d", iSM);
    htitle = TString::Format("SM%d L1 int. amplitude;eta id;phi id;amplitude", iSM);
    fHistFastORL1AmpBySM[iSM] = new TH2F(hname, htitle, 24, 0, 24, 12, 0, 12);
    fHistograms.Add(fHistFastORL1AmpBySM[iSM]);

    hname = TString::Format("EMCTRQA_histFEEvsSTU_SM%d", iSM);
    htitle = TString::Format("SM%d FEE vs STU;STU amplitude;FEE energy (GeV)", iSM);
    fHistFEEvsSTUBySM[iSM] = new TH2F(hname, htitle, 128, 0, 2048, 150, 0, 150);
    fHistograms.Add(fHistFEEvsSTUBySM[iSM]);
  }

  const char* det[fgkNDet] = { "EMCal", "DCal" };

  for (Int_t itrig = 0; itrig < fgkNTriggerTypes; itrig++) {
    if (EMCALTrigger::kEMCalTriggerNames[itrig].IsNull() || fEnabledTriggerTypes[itrig] == kFALSE) continue;

    for (Int_t ipatch = 0; ipatch < fgkNPatchTypes; ipatch++) {
      if (!fEnabledPatchTypes[ipatch]) continue;
      for (Int_t idet = 0; idet < fgkNDet; idet++) {
        hname = TString::Format("EMCTRQA_hist%sNPatches%s%s", det[idet], EMCALTrigger::kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[ipatch].Data());
        htitle = TString::Format("EMCTRQA_hist%sNPatches%s%s;num. of patches;events", det[idet], EMCALTrigger::kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[ipatch].Data());
        fHistNPatches[idet][itrig][ipatch] = new TH1F(hname, htitle, 100, 0, 5000);
        fHistograms.Add(fHistNPatches[idet][itrig][ipatch]);

        hname = TString::Format("EMCTRQA_hist%sMaxPatchAmp%s%s", det[idet], EMCALTrigger::kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[ipatch].Data());
        htitle = TString::Format("EMCTRQA_hist%sMaxPatchAmp%s%s;amplitude;entries", det[idet], EMCALTrigger::kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[ipatch].Data());
        fHistMaxPatchAmp[idet][itrig][ipatch] = new TH1F(hname, htitle, fgkMaxPatchAmp[itrig]/fADCperBin, 0, fgkMaxPatchAmp[itrig]);
        fHistograms.Add(fHistMaxPatchAmp[idet][itrig][ipatch]);

        hname = TString::Format("EMCTRQA_hist%sPatchAmp%s%s", det[idet], EMCALTrigger::kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[ipatch].Data());
        htitle = TString::Format("EMCTRQA_hist%sPatchAmp%s%s;amplitude;entries", det[idet], EMCALTrigger::kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[ipatch].Data());
        fHistPatchAmp[idet][itrig][ipatch] = new TH1F(hname, htitle, fgkMaxPatchAmp[itrig]/fADCperBin, 0, fgkMaxPatchAmp[itrig]);
        fHistograms.Add(fHistPatchAmp[idet][itrig][ipatch]);
      }

      hname = TString::Format("EMCTRQA_histMaxEdgePos%s%s", EMCALTrigger::kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[ipatch].Data());
      htitle = TString::Format("Edge Position Max %s patch %s;col;row;entries", EMCALTrigger::kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[ipatch].Data());
      fHistMaxEdgePos[itrig][ipatch] = new TH2F(hname, htitle, 48, 0, 48, 105, 0, 105);
      fHistograms.Add(fHistMaxEdgePos[itrig][ipatch]);

      hname = TString::Format("EMCTRQA_histAmpEdgePos%s%s", EMCALTrigger::kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[ipatch].Data());
      htitle = TString::Format("Integrated amplitude %s patch %s;col;row;entries", EMCALTrigger::kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[ipatch].Data());
      fHistAmpEdgePos[itrig][ipatch] = new TH2F(hname, htitle, 48, 0, 48, 105, 0, 105);
      fHistograms.Add(fHistAmpEdgePos[itrig][ipatch]);
    }
  }
}

/**
 * Process a patch, filling relevant histograms.
 * \param patch Pointer to a valid trigger patch
 */
void AliEMCALTriggerOnlineQAPP::ProcessPatch(const AliEMCALTriggerPatchInfo* patch)
{
  Int_t triggerBits[6] = { patch->GetTriggerBitConfig()->GetLevel0Bit(),
      patch->GetTriggerBitConfig()->GetGammaLowBit(),
      patch->GetTriggerBitConfig()->GetGammaHighBit(),
      patch->GetTriggerBitConfig()->GetJetLowBit(),
      patch->GetTriggerBitConfig()->GetJetHighBit(),
      patch->GetTriggerBitConfig()->GetBkgBit()
  };

  Int_t offsets[fgkNPatchTypes]    = { 0, AliEMCALTriggerPatchInfo::kRecalcOffset, AliEMCALTriggerPatchInfo::kOfflineOffset };
  Int_t amplitudes[fgkNPatchTypes] = { patch->GetADCAmp(),  patch->GetADCAmp(),  patch->GetADCOfflineAmp() };

  for (Int_t itrig = 0; itrig < fgkNTriggerTypes; itrig++) {
    if (EMCALTrigger::kEMCalTriggerNames[itrig].IsNull() || fEnabledTriggerTypes[itrig] == kFALSE) continue;

    for (Int_t ipatch = 0; ipatch < fgkNPatchTypes; ipatch++) {
      if (!fEnabledPatchTypes[ipatch]) continue;
      if(!ipatch){
        if (!(patch->TestTriggerBit(triggerBits[itrig]+offsets[ipatch]) ||
            patch->TestTriggerBit(triggerBits[itrig] + patch->GetTriggerBitConfig()->GetTriggerTypesEnd() + offsets[ipatch]))) continue;
      }
      else
        if (!patch->TestTriggerBit(triggerBits[itrig]+offsets[ipatch])) continue;

      Int_t idet = 0;
      if (patch->IsEMCal()) {
        idet = 0;
        if (GetAmplitude(fMaxPatchEMCal[itrig][ipatch], ipatch) < amplitudes[ipatch]) {
          *(fMaxPatchEMCal[itrig][ipatch]) = *patch;
        }
      }
      else if (patch->IsDCalPHOS()) {
        idet = 1;
        if (GetAmplitude(fMaxPatchDCal[itrig][ipatch], ipatch) < amplitudes[ipatch]) *(fMaxPatchDCal[itrig][ipatch]) = *patch;
      }
      else {
        AliWarning(Form("Patch is not EMCal nor DCal/PHOS (pos: %d, %d)", patch->GetRowStart(), patch->GetColStart()));
      }

      fNPatches[idet][itrig][ipatch]++;
      fHistPatchAmp[idet][itrig][ipatch]->Fill(amplitudes[ipatch]);
      fHistAmpEdgePos[itrig][ipatch]->Fill(patch->GetColStart(), patch->GetRowStart(), amplitudes[ipatch]);
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
void AliEMCALTriggerOnlineQAPP::ProcessFastor(const AliEMCALTriggerFastOR* fastor, AliVCaloCells* cells)
{
  UInt_t L0amp = fastor->GetL0Amp();
  UInt_t L1amp = fastor->GetL1Amp();

  if (L0amp > fMinL0FastORAmp) {
    fHistFastORL0->Fill(fastor->GetAbsId());
    if (L0amp > fFastorL0Th) fHistFastORL0LargeAmp->Fill(fastor->GetAbsId());
    fHistFastORL0Amp->Fill(fastor->GetAbsId(), L0amp);
    fHistFastORL0Time->Fill(fastor->GetAbsId(), fastor->GetL0Time());
  }

  if (L1amp > fMinL1FastORAmp) {
    fHistFastORL1->Fill(fastor->GetAbsId());
    if (L1amp > fFastorL1Th) fHistFastORL1LargeAmp->Fill(fastor->GetAbsId());
    fHistFastORL1Amp->Fill(fastor->GetAbsId(), L1amp);
  }

  if (!fGeom) return;
  // After this only instructions that require geometry
  Int_t iSM  = -1;
  Int_t iEta = -1;
  Int_t iPhi = -1;
  fGeom->GetPositionInSMFromAbsFastORIndex(fastor->GetAbsId(), iSM, iEta, iPhi);
  Bool_t isDCal = fGeom->IsDCALSM(iSM);

  if (iSM >=0 && iSM < fgkSM) {
    if (L0amp > fMinL0FastORAmp) {
      fHistFastORL0BySM[iSM]->Fill(iEta, iPhi);
      if (L0amp > fFastorL0Th) fHistFastORL0LargeAmpBySM[iSM]->Fill(iEta, iPhi);
      fHistFastORL0AmpBySM[iSM]->Fill(iEta, iPhi, L0amp);
    }

    if (L1amp > fMinL1FastORAmp) {
      fHistFastORL1BySM[iSM]->Fill(iEta, iPhi);
      if (L1amp > fFastorL1Th) fHistFastORL1LargeAmpBySM[iSM]->Fill(iEta, iPhi);
      fHistFastORL1AmpBySM[iSM]->Fill(iEta, iPhi, L1amp);
    }
  }

  if (!cells) return;
  // After this only instructions that require geometry & cell information
  Double_t offlineAmp = 0;
  Int_t idx[4] = {-1};
  if (fGeom->GetCellIndexFromFastORIndex(fastor->GetAbsId(), idx)) {
    for (Int_t i = 0; i < 4; i++) offlineAmp += cells->GetCellAmplitude(idx[i]);
  }

  if (iSM >=0 && iSM < fgkSM) {
    if (L0amp > fMinL0FastORAmp) fHistFEEvsTRUBySM[iSM]->Fill(L0amp, offlineAmp);
    if (L1amp > fMinL0FastORAmp) fHistFEEvsSTUBySM[iSM]->Fill(L1amp, offlineAmp);
  }
}

/**
 * This method should be called at the end of each event.
 */
void AliEMCALTriggerOnlineQAPP::EventCompleted()
{
  fHistEvents->Fill(0);

  enum {kEMCAL=0,kDCAL=1};
  for (Int_t itrig = 0; itrig < fgkNTriggerTypes; itrig++) {
    if (EMCALTrigger::kEMCalTriggerNames[itrig].IsNull() || fEnabledTriggerTypes[itrig] == kFALSE) continue;

    for (Int_t ipatch = 0; ipatch < fgkNPatchTypes; ipatch++) {
      if (!fEnabledPatchTypes[ipatch]) continue;

      fHistNPatches[kEMCAL][itrig][ipatch]->Fill(fNPatches[kEMCAL][itrig][ipatch]);
      fHistNPatches[kDCAL][itrig][ipatch]->Fill(fNPatches[kDCAL][itrig][ipatch]);

      if (fMaxPatchEMCal[itrig][ipatch]->GetColStart() >= 0) {
        fHistMaxEdgePos[itrig][ipatch]->Fill(fMaxPatchEMCal[itrig][ipatch]->GetColStart(),
            fMaxPatchEMCal[itrig][ipatch]->GetRowStart());
        fHistMaxPatchAmp[kEMCAL][itrig][ipatch]->Fill(GetAmplitude(fMaxPatchEMCal[itrig][ipatch], ipatch));
      }

      if (fMaxPatchDCal[itrig][ipatch]->GetColStart() >= 0) {
        fHistMaxEdgePos[itrig][ipatch]->Fill(fMaxPatchDCal[itrig][ipatch]->GetColStart(),
            fMaxPatchDCal[itrig][ipatch]->GetRowStart());
        fHistMaxPatchAmp[kDCAL][itrig][ipatch]->Fill(GetAmplitude(fMaxPatchDCal[itrig][ipatch], ipatch));
      }

      fNPatches[kEMCAL][itrig][ipatch] = 0;
      fNPatches[kDCAL][itrig][ipatch] = 0;
      fMaxPatchEMCal[itrig][ipatch]->Reset();
      fMaxPatchDCal[itrig][ipatch]->Reset();
    }
  }
}
