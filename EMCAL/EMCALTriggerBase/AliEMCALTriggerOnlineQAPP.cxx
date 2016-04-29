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
  for (Int_t itrigger = 0; itrigger < fgkNTriggerTypes; itrigger++) {
    for (Int_t ipatch = 0; ipatch < fgkNPatchTypes; ipatch++) {
      fMaxPatchEMCal[itrigger][ipatch] = new AliEMCALTriggerPatchInfo;
      fMaxPatchDCal[itrigger][ipatch] = new AliEMCALTriggerPatchInfo;

      fHistMaxEdgePos[itrigger][ipatch] = 0;
      fHistAmpEdgePos[itrigger][ipatch] = 0;
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

  for (Int_t itrigger = 0; itrigger < fgkNTriggerTypes; itrigger++) {
    for (Int_t ipatch = 0; ipatch < fgkNPatchTypes; ipatch++) {
      fMaxPatchEMCal[itrigger][ipatch] = new AliEMCALTriggerPatchInfo;
      fMaxPatchDCal[itrigger][ipatch] = new AliEMCALTriggerPatchInfo;

      fHistMaxEdgePos[itrigger][ipatch] = 0;
      fHistAmpEdgePos[itrigger][ipatch] = 0;
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

  for (Int_t itrigger = 0; itrigger < fgkNTriggerTypes; itrigger++) {

    for (Int_t ipatch = 0; ipatch < fgkNPatchTypes; ipatch++) {
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
  htitle = "L0 amplitudes;col;row;amplitude";
  fHistFastORL0Amp = new TH2F(hname, htitle, 48, 0, 48, 104, 0, 104);
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
  htitle = "L1 amplitudes;col;row;amplitude";
  fHistFastORL1Amp = new TH2F(hname, htitle, 48, 0, 48, 104, 0, 104);
  fHistograms.Add(fHistFastORL1Amp);

  Int_t nSM = fgkSM;
  if (fGeom) nSM = fGeom->GetNumberOfSuperModules();

  const char* det[fgkNDet] = { "EMCal", "DCal" };

  for (Int_t itrig = 0; itrig < fgkNTriggerTypes; itrig++) {
    for (Int_t ipatch = 0; ipatch < fgkNPatchTypes; ipatch++) {
      if (!IsPatchTypeEnabled(ipatch, itrig)) continue;

      hname = TString::Format("EMCTRQA_histMaxEdgePos%s%s", EMCALTrigger::kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[ipatch].Data());
      htitle = TString::Format("Edge Position Max %s patch %s;col;row;entries", EMCALTrigger::kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[ipatch].Data());
      fHistMaxEdgePos[itrig][ipatch] = new TH2F(hname, htitle, 48, 0, 48, 104, 0, 104);
      fHistograms.Add(fHistMaxEdgePos[itrig][ipatch]);

      hname = TString::Format("EMCTRQA_histAmpEdgePos%s%s", EMCALTrigger::kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[ipatch].Data());
      htitle = TString::Format("Integrated amplitude %s patch %s;col;row;entries", EMCALTrigger::kEMCalTriggerNames[itrig].Data(), fgkPatchTypes[ipatch].Data());
      fHistAmpEdgePos[itrig][ipatch] = new TH2F(hname, htitle, 48, 0, 48, 104, 0, 104);
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
  Bool_t (AliEMCALTriggerPatchInfo::* triggerCheck[3][6])(void) const = {
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

  if (fDebugLevel >= 2) {
    Printf("Processing patch -> global pos = (%d, %d); Amp (online) = %d; Amp (offline) = %d; \n",
        patch->GetRowStart(), patch->GetColStart(), patch->GetADCAmp(), patch->GetADCOfflineAmp());
  }

  for (Int_t itrig = 0; itrig < fgkNTriggerTypes; itrig++) {
    for (Int_t ipatch = 0; ipatch < fgkNPatchTypes; ipatch++) {
      if (!IsPatchTypeEnabled(ipatch, itrig)) continue;
      if (!(patch->*(triggerCheck[ipatch][itrig]))()) continue;

      Int_t amp = GetAmplitude(patch, ipatch);
      Int_t idet = 0;
      if (patch->IsEMCal()) {
        idet = 0;
        if (GetAmplitude(fMaxPatchEMCal[itrig][ipatch], ipatch) < amp) {
          *(fMaxPatchEMCal[itrig][ipatch]) = *patch;
        }
      }
      else if (patch->IsDCalPHOS()) {
        idet = 1;
        if (GetAmplitude(fMaxPatchDCal[itrig][ipatch], ipatch) < amp) {
          *(fMaxPatchDCal[itrig][ipatch]) = *patch;
        }
      }
      else {
        AliWarning(Form("Patch is not EMCal nor DCal/PHOS (pos: %d, %d)", patch->GetRowStart(), patch->GetColStart()));
      }

      fHistAmpEdgePos[itrig][ipatch]->Fill(patch->GetColStart(), patch->GetRowStart(), amp);
    }
  }
}

/**
 * Process a FastOR, filling relevant histograms.
 * \param patch Pointer to a valid trigger FastOR
 */
void AliEMCALTriggerOnlineQAPP::ProcessFastor(const AliEMCALTriggerFastOR* fastor, AliVCaloCells* /*cells*/)
{
  UInt_t L0amp = fastor->GetL0Amp();
  UInt_t L1amp = fastor->GetL1Amp();

  if (L0amp > fMinL0FastORAmp) {
    fHistFastORL0->Fill(fastor->GetAbsId());
    if (L0amp > fFastorL0Th) fHistFastORL0LargeAmp->Fill(fastor->GetAbsId());
    fHistFastORL0Amp->Fill(fastor->GetGlobalCol(), fastor->GetGlobalRow(), L0amp);
    fHistFastORL0Time->Fill(fastor->GetAbsId(), fastor->GetL0Time());
  }

  if (L1amp > fMinL1FastORAmp) {
    fHistFastORL1->Fill(fastor->GetAbsId());
    if (L1amp > fFastorL1Th) fHistFastORL1LargeAmp->Fill(fastor->GetAbsId());
    fHistFastORL1Amp->Fill(fastor->GetGlobalCol(), fastor->GetGlobalRow(), L1amp);
  }

  if (!fGeom) return;
  // After this only instructions that require geometry
  Int_t iSM  = -1;
  Int_t iEta = -1;
  Int_t iPhi = -1;
  fGeom->GetPositionInSMFromAbsFastORIndex(fastor->GetAbsId(), iSM, iEta, iPhi);
  Bool_t isDCal = fGeom->IsDCALSM(iSM);
}

/**
 * This method should be called at the end of each event.
 */
void AliEMCALTriggerOnlineQAPP::EventCompleted()
{
  fHistEvents->Fill(0);

  enum {kEMCAL=0,kDCAL=1};
  for (Int_t itrig = 0; itrig < fgkNTriggerTypes; itrig++) {
    for (Int_t ipatch = 0; ipatch < fgkNPatchTypes; ipatch++) {
      if (!IsPatchTypeEnabled(ipatch, itrig)) continue;

      if (fMaxPatchEMCal[itrig][ipatch]->GetColStart() >= 0) {
        fHistMaxEdgePos[itrig][ipatch]->Fill(fMaxPatchEMCal[itrig][ipatch]->GetColStart(),
            fMaxPatchEMCal[itrig][ipatch]->GetRowStart());
        }

      if (fMaxPatchDCal[itrig][ipatch]->GetColStart() >= 0) {
        fHistMaxEdgePos[itrig][ipatch]->Fill(fMaxPatchDCal[itrig][ipatch]->GetColStart(),
            fMaxPatchDCal[itrig][ipatch]->GetRowStart());
      }

      fMaxPatchEMCal[itrig][ipatch]->Reset();
      fMaxPatchDCal[itrig][ipatch]->Reset();
    }
  }
}
