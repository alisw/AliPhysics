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
#include <memory>
#include "THistManager.h"

#include "AliAnalysisTaskEmcalNoiseTriggers.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALTriggerMapping.h"
#include "AliEMCALTriggerPatchInfo.h"
#include "AliEMCALTriggerPatchADCInfoAP.h"
#include "AliInputEventHandler.h"
#include "AliLog.h"
#include "AliVCaloTrigger.h"
#include "AliVEvent.h"

/// \cond CLASSIMP
ClassImp(EMCalTriggerPtAnalysis::AliAnalysisTaskEmcalNoiseTriggers)
/// \endcond

namespace EMCalTriggerPtAnalysis {

const TString AliAnalysisTaskEmcalNoiseTriggers::fgkPatchNames[2] = {"Online", "Recalc"};

AliAnalysisTaskEmcalNoiseTriggers::AliAnalysisTaskEmcalNoiseTriggers():
    AliAnalysisTaskEmcalTriggerBase(),
    fTriggerBits(0),
    fTriggerString("")
{
  SetSelectNoiseEvents(true);
}

AliAnalysisTaskEmcalNoiseTriggers::AliAnalysisTaskEmcalNoiseTriggers(const char *name):
    AliAnalysisTaskEmcalTriggerBase(name),
    fTriggerBits(0),
    fTriggerString("")
{
  SetSelectNoiseEvents(true);
}

void AliAnalysisTaskEmcalNoiseTriggers::CreateUserHistos(){
  const int kMaxCol = 48, kMaxRow = 104, kMaxFastOr = kMaxRow * kMaxCol;
  const TString kMaxTypes[2] = {"All", "Max"};

  // Event level histograms
  fHistos->CreateTH1("hEventCount", "EventCounter", 1, 0.5, 1.5);
  fHistos->CreateTH1("hVertexZ","",  100, -40, 40);

  // Histograms at fastor level
  fHistos->CreateTH2("hFastorL1TimeSums", "FastOR L1 time sum distribution",  kMaxFastOr, -0.5, kMaxFastOr - 0.5, 2049, -0.5, 2048.5);

  // Histograms at trigger patch level
  for(int itype = 0; itype < 2; itype++){
    const TString &mypatchtype = fgkPatchNames[itype];
    fHistos->CreateTH1(Form("hEventSelPatch%s", mypatchtype.Data()), Form("Number of firing %s patches", mypatchtype.Data()), 1001, -0.5, 1000.5);
    for(int imax = 0; imax < 2; imax++){
      fHistos->CreateTH2(Form("hPatchMaxFastorADCvsSumADC%s%s", kMaxTypes[imax].Data(), mypatchtype.Data()), Form("Max FastOR ADC vs. ADC sum for %s %s patches; Sum ADC; Max FastOR ADC", kMaxTypes[imax].Data(), mypatchtype.Data()), 2048, 0., 2048, 2500, 0., 2500);
      fHistos->CreateTH2(Form("hSumADCvsPatchADC%s%s", kMaxTypes[imax].Data(), mypatchtype.Data()), Form("Sum all ADC vs. Patch ADC for %s %s patches; sum ADC; patch ADC", kMaxTypes[imax].Data(), mypatchtype.Data()), 2048, 0., 2048, 2500, 0., 2500);
      fHistos->CreateTH2(Form("hNfastorVsFracMaxFastor%s%s", kMaxTypes[imax].Data(), mypatchtype.Data()), Form("Number of non-zero FastORs vs. ADC fraction of the highest FastOR for %s %s patches; Number of FastORs; Fraction ADC highest Fastor", kMaxTypes[imax].Data(), mypatchtype.Data()), 257, -0.5, 256.5, 100, 0., 1);
    }
  }
}

bool AliAnalysisTaskEmcalNoiseTriggers::IsUserEventSelected(){
  if(!(fInputHandler->IsEventSelected() & fTriggerBits)) return false;
  if(fTriggerString.Length()){
    return fInputEvent->GetFiredTriggerClasses().Contains(fTriggerString);
  }
  return true;
}

bool AliAnalysisTaskEmcalNoiseTriggers::Run(){
  if(!fL1ADC.IsAllocated()){
    int nrows = 64;
    if(fGeom->GetTriggerMappingVersion() == 2) nrows = 104;
    fL1ADC.Allocate(48, nrows);
  }
  PrepareL1FastorADC();

  // Look at fastors inside masked events
  AnalyseFastors();

  // Loop patches
  // In addition mark max patches for the case online and recalc - these
  // patches will be processed later to investigate max patch distributions
  AliEMCALTriggerPatchInfo *maxonline(nullptr), *maxrecalc(nullptr);
  Int_t nselOnline(0), nselRecalc(0), ngoodOnline(0), ngoodRecalc(0);
  for(auto patchit : *(this->fTriggerPatchInfo)){
    AliEMCALTriggerPatchInfo *recpatch = static_cast<AliEMCALTriggerPatchInfo *>(patchit);
    if(fTriggerBits & AliVEvent::kEMCEGA){
      if(fTriggerString.Contains("EG1")){
        if(recpatch->IsGammaHigh()){
          if(AnalyseTriggerPatch(*recpatch, kOnline, false)) nselOnline++;
          if(SelectFiredPatch("EG1", recpatch->GetADCAmp())) ngoodOnline++;
          if(!maxonline || recpatch->GetADCAmp() > maxonline->GetADCAmp()) maxonline = recpatch;
        }
        if(recpatch->IsGammaLowRecalc()){
          if(AnalyseTriggerPatch(*recpatch, kRecalc, false)) nselRecalc++;
          if(SelectFiredPatch("EG1", recpatch->GetADCAmp())) ngoodRecalc++;
          if(!maxrecalc || recpatch->GetADCAmp() > maxrecalc->GetADCAmp()) maxrecalc = recpatch;
        }
      }
      if(fTriggerString.Contains("EG2")){
        if(recpatch->IsGammaLow()){
          if(AnalyseTriggerPatch(*recpatch, kOnline, false)) nselOnline++;
          if(SelectFiredPatch("EG2", recpatch->GetADCAmp())) ngoodOnline++;
          if(!maxonline || recpatch->GetADCAmp() > maxonline->GetADCAmp()) maxonline = recpatch;
        }
        if(recpatch->IsGammaLowRecalc()){
          if(AnalyseTriggerPatch(*recpatch, kRecalc, false)) nselRecalc++;
          if(SelectFiredPatch("EG2", recpatch->GetADCAmp())) ngoodRecalc++;
          if(!maxrecalc || recpatch->GetADCAmp() > maxrecalc->GetADCAmp()) maxrecalc = recpatch;
        }
      }
    } else if(fTriggerBits & AliVEvent::kEMCEJE){
      if(fTriggerString.Contains("EJ1")){
        if(recpatch->IsJetHigh()){
          if(AnalyseTriggerPatch(*recpatch, kOnline, false)) nselOnline++;
          if(SelectFiredPatch("EJ1", recpatch->GetADCAmp())) ngoodOnline++;
          if(!maxonline || recpatch->GetADCAmp() > maxonline->GetADCAmp()) maxonline = recpatch;
        }
        if(recpatch->IsJetLowRecalc()){
          if(AnalyseTriggerPatch(*recpatch, kRecalc, false)) nselRecalc++;
          if(SelectFiredPatch("EJ1", recpatch->GetADCAmp())) ngoodRecalc++;
          if(!maxrecalc || recpatch->GetADCAmp() > maxrecalc->GetADCAmp()) maxrecalc = recpatch;
        }
      }
      if(fTriggerString.Contains("EJ2")){
        if(recpatch->IsJetLow()){
          if(AnalyseTriggerPatch(*recpatch, kOnline, false)) nselOnline++;
          if(SelectFiredPatch("EJ2", recpatch->GetADCAmp())) ngoodOnline++;
          if(!maxonline || recpatch->GetADCAmp() > maxonline->GetADCAmp()) maxonline = recpatch;
        }
        if(recpatch->IsJetLowRecalc()){
          if(AnalyseTriggerPatch(*recpatch, kRecalc, false)) nselRecalc++;
          if(SelectFiredPatch("EJ2", recpatch->GetADCAmp())) ngoodRecalc++;
          if(!maxrecalc || recpatch->GetADCAmp() > maxrecalc->GetADCAmp()) maxrecalc = recpatch;
        }
      }
    }
  }

  // Fill event statistics
  fHistos->FillTH1("hEventSelPatchOnline", nselOnline);
  fHistos->FillTH1("hEventSelPatchRecalc", nselRecalc);
  fHistos->FillTH1("hEventGoodPatchOnline", ngoodOnline);
  fHistos->FillTH1("hEventGoodPatchRecalc", ngoodRecalc);

  if(maxonline) AnalyseTriggerPatch(*maxonline, kOnline, true);
  if(maxrecalc) AnalyseTriggerPatch(*maxrecalc, kRecalc, true);

  return true;
}

void AliAnalysisTaskEmcalNoiseTriggers::PrepareL1FastorADC(){
  fL1ADC.Reset();
  AliVCaloTrigger *emctrigger = fInputEvent->GetCaloTrigger("EMCAL");
  emctrigger->Reset();

  Int_t globCol=-1, globRow=-1, adcAmp=-1;
  while(emctrigger->Next()){
    // get position in global 2x2 tower coordinates
    // A0 left bottom (0,0)
    emctrigger->GetPosition(globCol, globRow);
    emctrigger->GetL1TimeSum(adcAmp);
    AliDebugStream(1) << GetName() << "Fastor at (" << globCol << "," << globRow << ") with ADC " << adcAmp << std::endl;
    if (adcAmp < 0) adcAmp = 0;

    try {
      (fL1ADC)(globCol,globRow) = adcAmp;
    }
    catch (AliEMCALTriggerDataGrid<Int_t>::OutOfBoundsException &e) {
      std::string dirstring = e.GetDirection() == AliEMCALTriggerDataGrid<Int_t>::OutOfBoundsException::kColDir ? "Col" : "Row";
      AliErrorStream() << "Trigger maker task - filling trigger bit grid - index out-of-bounds in " << dirstring << ": " << e.GetIndex() << std::endl;
    }
  }
}

AliEMCALTriggerPatchADCInfoAP *AliAnalysisTaskEmcalNoiseTriggers::MakeFastorADCValuesForPatch(const AliEMCALTriggerPatchInfo &patch ) const {
  AliEMCALTriggerPatchADCInfoAP *adcpatch = new AliEMCALTriggerPatchADCInfoAP(patch.GetPatchSize());
  for(unsigned char icol = 0; icol < patch.GetPatchSize(); icol++){
    for(unsigned char irow = 0; irow < patch.GetPatchSize(); irow++){
      Int_t adc = 0;
      try{
        adc = fL1ADC(icol + patch.GetColStart(), irow + patch.GetRowStart());
      } catch (AliEMCALTriggerDataGrid<Int_t>::OutOfBoundsException &e){
        adc = 0;
      }
      adcpatch->SetADC(adc, icol, irow);
    }
  }
  return adcpatch;
}

void AliAnalysisTaskEmcalNoiseTriggers::UserFillHistosAfterEventSelection(){
  fHistos->FillTH1("hEventCount", 1);
  fHistos->FillTH1("hVertexZ", fVertex[2]);
}

void AliAnalysisTaskEmcalNoiseTriggers::AnalyseFastors(){
  AliVCaloTrigger *emctrigger = fInputEvent->GetCaloTrigger("EMCAL");
  emctrigger->Reset();

  Int_t globCol=-1, globRow=-1, adcAmp=-1, absFastor = -1;
  while(emctrigger->Next()){
    emctrigger->GetPosition(globCol, globRow);
    emctrigger->GetL1TimeSum(adcAmp);
    if(adcAmp < 0) adcAmp = 0;
    fGeom->GetTriggerMapping()->GetAbsFastORIndexFromPositionInEMCAL(globCol, globRow, absFastor);
    fHistos->FillTH2("hFastorL1TimeSums", absFastor, adcAmp);
  }
}

Bool_t AliAnalysisTaskEmcalNoiseTriggers::AnalyseTriggerPatch(const AliEMCALTriggerPatchInfo &recpatch, SelectPatchType_t pt, Bool_t maxpatch){
  std::unique_ptr<AliEMCALTriggerPatchADCInfoAP> adcvalues(MakeFastorADCValuesForPatch(recpatch));

  // cut patch according to ADC sum without masking (select also patches with noise)
  Int_t sumADC = adcvalues->GetSumADC();
  if(pt == kRecalc) {
    if(!SelectFiredPatch(fTriggerString, sumADC)) return kFALSE;
  }

  const TString &mypatchtype =  fgkPatchNames[pt], maxname = maxpatch ? "Max" : "All";
  Int_t maxADC = adcvalues->GetMaxADC(), ncontrib = adcvalues->GetNFastorsContrib();
  Float_t fracMaxFastor = static_cast<Float_t>(maxADC)/static_cast<Float_t>(sumADC);
  fHistos->FillTH2(Form("hPatchMaxFastorADCvsSumADC%s%s", maxname.Data(), mypatchtype.Data()), sumADC, maxADC);
  fHistos->FillTH2(Form("hSumADCvsPatchADC%s%s", maxname.Data(), mypatchtype.Data()), sumADC, recpatch.GetADCAmp());
  fHistos->FillTH2(Form("hNfastorVsFracMaxFastor%s%s", maxname.Data(), mypatchtype.Data()), ncontrib, fracMaxFastor);
  return kTRUE;
}

} /* namespace EMCalTriggerPtAnalysis */
