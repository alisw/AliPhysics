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
#include <TClonesArray.h>
#include <THashList.h>

#include "AliEmcalTriggerPatchInfoAPV1.h"
#include "AliEmcalTriggerQAAP.h"
#include "AliEmcalTriggerFastORAP.h"

#include "AliEmcalTriggerQATask.h"

/// \cond CLASSIMP
ClassImp(AliEmcalTriggerQATask)
/// \endcond

/**
 * Dummy constructor
 */
AliEmcalTriggerQATask::AliEmcalTriggerQATask() : 
  AliAnalysisTaskEmcal("AliEmcalTriggerQATask",kTRUE),
  fTriggerPatchesName("EmcalTriggers"),
  fEMCALTriggerQA(0),
  fBadChannels(),
  fTriggerPatches(0),
  fHistEMCalTriggers(0)
{
}

/**
 * Named constructor.
 * \param name Name of the trigger QA task
 */
AliEmcalTriggerQATask::AliEmcalTriggerQATask(const char *name) :
  AliAnalysisTaskEmcal(name,kTRUE),
  fTriggerPatchesName("EmcalTriggers"),
  fEMCALTriggerQA(0),
  fBadChannels(),
  fTriggerPatches(0),
  fHistEMCalTriggers(0)
{
  // Constructor.
  SetMakeGeneralHistograms(kTRUE);

  TString qaName(Form("%s_AliEmcalTriggerQAAP",name));
  fEMCALTriggerQA = new AliEmcalTriggerQAAP(qaName);
}

/**
 * Destructor
 */
AliEmcalTriggerQATask::~AliEmcalTriggerQATask()
{
  delete fEMCALTriggerQA;
}

/**
 * Init the analysis.
 */
void AliEmcalTriggerQATask::ExecOnce()
{
  AliAnalysisTaskEmcal::ExecOnce();

  if (!fInitialized) return;

  fTriggerPatches = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fTriggerPatchesName));

  TString objname(fTriggerPatches->GetClass()->GetName());
  TClass cls(objname);
  if (!cls.InheritsFrom("AliEmcalTriggerPatchInfoAPV1")) {
    AliError(Form("%s: Objects of type %s in %s are not inherited from AliEmcalTriggerPatchInfoAPV1!",
		  GetName(), cls.GetName(), fTriggerPatchesName.Data()));
    fTriggerPatches = 0;
  }

  if (!fTriggerPatches) {
    fInitialized = kFALSE;
    AliError(Form("%s: Unable to get trigger patch container with name %s. Aborting", GetName(), fTriggerPatchesName.Data()));
    return;
  }
}

/**
 * Sets the calo trigger names used in 2015 data taking
 */
void AliEmcalTriggerQATask::Set2015CaloTriggerNames()
{
  fCaloTriggerNames[kMinBias] = "CINT7-B-NOPF-CENT";
  fCaloTriggerNames[kEMCalL0] = "CEMC7-B-NOPF-CENTNOPMD";
  fCaloTriggerNames[kEMCalL1G1] = "CINT7EG1-B-NOPF-CENTNOPMD";
  fCaloTriggerNames[kEMCalL1G2] = "CINT7EG2-B-NOPF-CENTNOPMD";
  fCaloTriggerNames[kEMCalL1J1] = "CINT7EJ1-B-NOPF-CENTNOPMD";
  fCaloTriggerNames[kEMCalL1J2] = "CINT7EJ2-B-NOPF-CENTNOPMD";
  fCaloTriggerNames[kDCalL0] = "CDMC7-B-NOPF-CENTNOPMD";
  fCaloTriggerNames[kDCalL1G1] = "CINT7DG1-B-NOPF-CENTNOPMD";
  fCaloTriggerNames[kDCalL1G2] = "CINT7DG2-B-NOPF-CENTNOPMD";
  fCaloTriggerNames[kDCalL1J1] = "CINT7DJ1-B-NOPF-CENTNOPMD";
  fCaloTriggerNames[kDCalL1J2] = "CINT7DJ2-B-NOPF-CENTNOPMD";
  fCaloTriggerNames[kPHOSL0] = "CPHI7-B-NOPF-CENTNOPMD";
  fCaloTriggerNames[kPHOSL1H] = "CINT7PHH-B-NOPF-CENTNOPMD";
  fCaloTriggerNames[kPHOSL1M] = "CINT7PHM-B-NOPF-CENTNOPMD";
  fCaloTriggerNames[kPHOSL1L] = "CINT7PHL-B-NOPF-CENTNOPMD";
}

/**
 * Create objects, histograms
 */
void AliEmcalTriggerQATask::UserCreateOutputObjects()
{
  AliAnalysisTaskEmcal::UserCreateOutputObjects();

  fEMCALTriggerQA->SetDebugLevel(DebugLevel());
  fEMCALTriggerQA->Init();

  if (fOutput) {  
    fHistEMCalTriggers = new TH1F("fHistEMCalTriggers","fHistEMCalTriggers; triggers; counts",40,0,40);
#if ROOT_VERSION_CODE < ROOT_VERSION(6,4,2)
    fHistEMCalTriggers->SetBit(TH1::kCanRebin);
#else
    fHistEMCalTriggers->SetCanExtend(TH1::kAllAxes);
#endif
    fHistEMCalTriggers->GetXaxis()->SetBinLabel(1, "EMCal L0");
    fHistEMCalTriggers->GetXaxis()->SetBinLabel(2, "EMCal L1 G1");
    fHistEMCalTriggers->GetXaxis()->SetBinLabel(3, "EMCal L1 G2");
    fHistEMCalTriggers->GetXaxis()->SetBinLabel(4, "EMCal L1 J1");
    fHistEMCalTriggers->GetXaxis()->SetBinLabel(5, "EMCal L1 J2");
    fHistEMCalTriggers->GetXaxis()->SetBinLabel(6, "EMCal L1 Any");
    fHistEMCalTriggers->GetXaxis()->SetBinLabel(7, "EMCal Any");

    fHistEMCalTriggers->GetXaxis()->SetBinLabel(8, "DCal L0");
    fHistEMCalTriggers->GetXaxis()->SetBinLabel(9, "DCal L1 G1");
    fHistEMCalTriggers->GetXaxis()->SetBinLabel(10, "DCal L1 G2");
    fHistEMCalTriggers->GetXaxis()->SetBinLabel(11, "DCal L1 J1");
    fHistEMCalTriggers->GetXaxis()->SetBinLabel(12, "DCal L1 J2");
    fHistEMCalTriggers->GetXaxis()->SetBinLabel(13, "DCal L1 Any");
    fHistEMCalTriggers->GetXaxis()->SetBinLabel(14, "DCal Any");

    fHistEMCalTriggers->GetXaxis()->SetBinLabel(15, "EMCal/DCal L0");
    fHistEMCalTriggers->GetXaxis()->SetBinLabel(16, "EMCal/DCal L1 G1");
    fHistEMCalTriggers->GetXaxis()->SetBinLabel(17, "EMCal/DCal L1 G2");
    fHistEMCalTriggers->GetXaxis()->SetBinLabel(18, "EMCal/DCal L1 J1");
    fHistEMCalTriggers->GetXaxis()->SetBinLabel(19, "EMCal/DCal L1 J2");
    fHistEMCalTriggers->GetXaxis()->SetBinLabel(20, "EMCal/DCal L1 Any");

    fHistEMCalTriggers->GetXaxis()->SetBinLabel(21, "EMCal/DCal Any");

    fHistEMCalTriggers->GetXaxis()->SetBinLabel(22, "PHOS L0");
    fHistEMCalTriggers->GetXaxis()->SetBinLabel(23, "PHOS L1 H");
    fHistEMCalTriggers->GetXaxis()->SetBinLabel(24, "PHOS L1 M");
    fHistEMCalTriggers->GetXaxis()->SetBinLabel(25, "PHOS L1 L");
    fHistEMCalTriggers->GetXaxis()->SetBinLabel(26, "PHOS L1 Any");
    fHistEMCalTriggers->GetXaxis()->SetBinLabel(27, "PHOS Any");

    fHistEMCalTriggers->GetXaxis()->SetBinLabel(28, "CALO L0");
    fHistEMCalTriggers->GetXaxis()->SetBinLabel(29, "CALO L1");
    fHistEMCalTriggers->GetXaxis()->SetBinLabel(30, "CALO Any");

    fHistEMCalTriggers->GetXaxis()->SetBinLabel(31, "MB");

    fHistEMCalTriggers->GetXaxis()->SetBinLabel(32, "MB or CALO");

    fOutput->Add(fHistEMCalTriggers);

    fOutput->Add(fEMCALTriggerQA->GetListOfHistograms());

    PostData(1, fOutput);
  }
}

/**
 * Run analysis.
 * \return Always true.
 */
Bool_t AliEmcalTriggerQATask::Run() 
{
  return kTRUE;
}

UInt_t AliEmcalTriggerQATask::SteerFiredTriggers(const TString& firedTriggersStr) const
{
  UInt_t firedTriggers = 0;

  for (Int_t bit = 0; bit < kLastCaloTrigger; bit++) {
    if (firedTriggersStr.Contains(fCaloTriggerNames[bit])) {
      SETBIT(firedTriggers, bit);
    }
  }

  return firedTriggers;
}


/**
 * Fill QA histograms
 * \return Always true.
 */
Bool_t AliEmcalTriggerQATask::FillHistograms() 
{
  UInt_t firedTriggerBits = SteerFiredTriggers(InputEvent()->GetFiredTriggerClasses());

  if ((firedTriggerBits & kEMCalL0bit) != 0) fHistEMCalTriggers->Fill("EMCal L0", 1);
  if ((firedTriggerBits & kEMCalL1G1bit) != 0) fHistEMCalTriggers->Fill("EMCal L1 G1", 1);
  if ((firedTriggerBits & kEMCalL1G2bit) != 0) fHistEMCalTriggers->Fill("EMCal L1 G2", 1);
  if ((firedTriggerBits & kEMCalL1J1bit) != 0) fHistEMCalTriggers->Fill("EMCal L1 J1", 1);
  if ((firedTriggerBits & kEMCalL1J2bit) != 0) fHistEMCalTriggers->Fill("EMCal L1 J2", 1);
  if ((firedTriggerBits & kEMCalL1Anybit) != 0) fHistEMCalTriggers->Fill("EMCal L1 Any", 1);
  if ((firedTriggerBits & kEMCalAnybit) != 0) fHistEMCalTriggers->Fill("EMCal Any", 1);

  if ((firedTriggerBits & kDCalL0bit) != 0) fHistEMCalTriggers->Fill("DCal L0", 1);
  if ((firedTriggerBits & kDCalL1G1bit) != 0) fHistEMCalTriggers->Fill("DCal L1 G1", 1);
  if ((firedTriggerBits & kDCalL1G2bit) != 0) fHistEMCalTriggers->Fill("DCal L1 G2", 1);
  if ((firedTriggerBits & kDCalL1J1bit) != 0) fHistEMCalTriggers->Fill("DCal L1 J1", 1);
  if ((firedTriggerBits & kDCalL1J2bit) != 0) fHistEMCalTriggers->Fill("DCal L1 J2", 1);
  if ((firedTriggerBits & kDCalL1Anybit) != 0) fHistEMCalTriggers->Fill("DCal L1 Any", 1);
  if ((firedTriggerBits & kDCalAnybit) != 0) fHistEMCalTriggers->Fill("DCal Any", 1);

  if ((firedTriggerBits & kEMCalDCalL0bit) != 0) fHistEMCalTriggers->Fill("EMCal/DCal L0", 1);
  if ((firedTriggerBits & kEMCalDCalL1G1bit) != 0) fHistEMCalTriggers->Fill("EMCal/DCal L1 G1", 1);
  if ((firedTriggerBits & kEMCalDCalL1G2bit) != 0) fHistEMCalTriggers->Fill("EMCal/DCal L1 G2", 1);
  if ((firedTriggerBits & kEMCalDCalL1J1bit) != 0) fHistEMCalTriggers->Fill("EMCal/DCal L1 J1", 1);
  if ((firedTriggerBits & kEMCalDCalL1J2bit) != 0) fHistEMCalTriggers->Fill("EMCal/DCal L1 J2", 1);
  if ((firedTriggerBits & kEMCalDCalL1Anybit) != 0) fHistEMCalTriggers->Fill("EMCal/DCal L1 Any", 1);

  if ((firedTriggerBits & kEMCalDCalAnybit) != 0) fHistEMCalTriggers->Fill("EMCal/DCal Any", 1);

  if ((firedTriggerBits & kPHOSL0bit) != 0) fHistEMCalTriggers->Fill("PHOS L0", 1);
  if ((firedTriggerBits & kPHOSL1Hbit) != 0) fHistEMCalTriggers->Fill("PHOS L1 H", 1);
  if ((firedTriggerBits & kPHOSL1Mbit) != 0) fHistEMCalTriggers->Fill("PHOS L1 M", 1);
  if ((firedTriggerBits & kPHOSL1Lbit) != 0) fHistEMCalTriggers->Fill("PHOS L1 L", 1);
  if ((firedTriggerBits & kPHOSL1Anybit) != 0) fHistEMCalTriggers->Fill("PHOS L1 Any", 1);
  if ((firedTriggerBits & kPHOSAnybit) != 0) fHistEMCalTriggers->Fill("PHOS Any", 1);

  if ((firedTriggerBits & kCALOL0bit) != 0) fHistEMCalTriggers->Fill("CALO L0", 1);
  if ((firedTriggerBits & kCALOL1bit) != 0) fHistEMCalTriggers->Fill("CALO L1", 1);
  if ((firedTriggerBits & kCALOAnybit) != 0) fHistEMCalTriggers->Fill("CALO Any", 1);

  if ((firedTriggerBits & kMinBiasbit) != 0) fHistEMCalTriggers->Fill("MB", 1);

  if ((firedTriggerBits & kCALOMinBias) != 0) fHistEMCalTriggers->Fill("MB or CALO", 1);

  if (fTriggerPatches) {
    Int_t nPatches = fTriggerPatches->GetEntriesFast();

    AliDebug(2, Form("nPatches = %d", nPatches));

    Int_t type = 0;

    for (Int_t i = 0; i < nPatches; i++) {
      AliDebug(2, Form("Processing patch %d", i));

      AliEmcalTriggerPatchInfoAPV1* patch = static_cast<AliEmcalTriggerPatchInfoAPV1*>(fTriggerPatches->At(i));
      if (!patch) continue;

      fEMCALTriggerQA->ProcessPatch(patch);
    }
  }

  if (fCaloTriggers) {
    AliEmcalTriggerFastORAP fastor;
    fCaloTriggers->Reset();
    Int_t globCol = -1, globRow = -1;
    Float_t L0amp = -1;
    Int_t L1amp = -1;
    while (fCaloTriggers->Next()) {
      // get position in global 2x2 tower coordinates
      // A0 left bottom (0,0)
      fCaloTriggers->GetPosition(globCol, globRow);
      // exclude channel completely if it is masked as hot channel
      if (fBadChannels.HasChannel(globCol, globRow)) continue;
      // for some strange reason some ADC amps are initialized in reconstruction
      // as -1, neglect those
      fCaloTriggers->GetL1TimeSum(L1amp);
      if (L1amp < 0) L1amp = 0;
      fCaloTriggers->GetAmplitude(L0amp);
      if (L0amp < 0) L0amp = 0;

      fastor.Initialize(L0amp, L1amp, globRow, globCol, fGeom);

      fEMCALTriggerQA->ProcessFastor(&fastor);
    }
  }
  fEMCALTriggerQA->EventCompleted();

  return kTRUE;
}
