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
#include <THnSparse.h>

#include "AliEMCALTriggerPatchInfo.h"
#include "AliEMCALTriggerFastOR.h"
#include "AliEMCALTriggerConstants.h"
#include "AliEmcalTriggerQATask.h"


using namespace EMCALTrigger;

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
  fADCperBin(20),
  fBkgPatchType(kTMEMCalBkg),
  fBadChannels(),
  fTriggerPatches(0),
  fHistEMCalTriggers(0),
  fHistEventQA(0)
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
  fADCperBin(20),
  fBkgPatchType(kTMEMCalBkg),
  fBadChannels(),
  fTriggerPatches(0),
  fHistEMCalTriggers(0),
  fHistEventQA(0)
{
  // Constructor.
  SetMakeGeneralHistograms(kTRUE);

  fEMCALTriggerQA = new TObjArray((fNcentBins+1)*2);
  fEMCALTriggerQA->SetOwner(kTRUE);

  for (Int_t i = 0; i < fNcentBins; i++) {
    TString qaName(Form("%s_AliEMCALTriggerQA_Cent%d", name, i));
    fEMCALTriggerQA->AddAt(new AliEMCALTriggerQA(qaName), i);
  }
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

  if (fTriggerPatches) {
    TString objname(fTriggerPatches->GetClass()->GetName());
    TClass cls(objname);
    if (!cls.InheritsFrom("AliEMCALTriggerPatchInfo")) {
      AliError(Form("%s: Objects of type %s in %s are not inherited from AliEMCALTriggerPatchInfo!",
          GetName(), cls.GetName(), fTriggerPatchesName.Data()));
      fTriggerPatches = 0;
    }
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

    fHistEMCalTriggers->GetXaxis()->SetBinLabel(33, "EMCal Any and !MB");
    fHistEMCalTriggers->GetXaxis()->SetBinLabel(34, "DCal Any and !MB");
    fHistEMCalTriggers->GetXaxis()->SetBinLabel(35, "EMCal/DCal Any and !MB");

    fOutput->Add(fHistEMCalTriggers);

    for (Int_t i = 0; i < fNcentBins; i++) {
      GetTriggerQA(i)->SetDebugLevel(DebugLevel());
      GetTriggerQA(i)->Init();
      fOutput->Add(GetTriggerQA(i)->GetListOfHistograms());
    }

    Int_t dim = 0;
    TString title[20];
    Int_t nbins[20] = {0};
    Double_t min[20] = {0};
    Double_t max[20] = {0};

    if (fForceBeamType != AliAnalysisTaskEmcal::kpp) {
      title[dim] = "Centrality %";
      nbins[dim] = 101;
      min[dim] = 0;
      max[dim] = 101;
      dim++;
    }

    title[dim] = "DCal median (offline)";
    nbins[dim] = AliEMCALTriggerQA::fgkMaxPatchAmp[fBkgPatchType] / fADCperBin / 10;
    min[dim] = 0;
    max[dim] = AliEMCALTriggerQA::fgkMaxPatchAmp[fBkgPatchType] / 10;
    dim++;

    title[dim] = "EMCal median (offline)";
    nbins[dim] = AliEMCALTriggerQA::fgkMaxPatchAmp[fBkgPatchType] / fADCperBin / 10;
    min[dim] = 0;
    max[dim] = AliEMCALTriggerQA::fgkMaxPatchAmp[fBkgPatchType] / 10;
    dim++;

    title[dim] = "DCal median (recalc)";
    nbins[dim] = AliEMCALTriggerQA::fgkMaxPatchAmp[fBkgPatchType] / fADCperBin / 10;
    min[dim] = 0;
    max[dim] = AliEMCALTriggerQA::fgkMaxPatchAmp[fBkgPatchType] / 10;
    dim++;

    title[dim] = "EMCal median (recalc)";
    nbins[dim] = AliEMCALTriggerQA::fgkMaxPatchAmp[fBkgPatchType] / fADCperBin / 10;
    min[dim] = 0;
    max[dim] = AliEMCALTriggerQA::fgkMaxPatchAmp[fBkgPatchType] / 10;
    dim++;

    fHistEventQA = new THnSparseF("fHistEventQA","fHistEventQA",dim,nbins,min,max);
    for (Int_t i = 0; i < dim; i++)
      fHistEventQA->GetAxis(i)->SetTitle(title[i]);
    fOutput->Add(fHistEventQA);

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

  if ((firedTriggerBits & kMinBiasbit) == 0 &&
      (firedTriggerBits & kEMCalAnybit) != 0) fHistEMCalTriggers->Fill("EMCal Any and !MB", 1);

  if ((firedTriggerBits & kMinBiasbit) == 0 &&
      (firedTriggerBits & kDCalAnybit) != 0) fHistEMCalTriggers->Fill("DCal Any and !MB", 1);

  if ((firedTriggerBits & kMinBiasbit) == 0 &&
      (firedTriggerBits & kEMCalDCalAnybit) != 0) fHistEMCalTriggers->Fill("EMCal/DCal Any and !MB", 1);

  if (fTriggerPatches) {
    Int_t nPatches = fTriggerPatches->GetEntriesFast();

    AliDebug(2, Form("nPatches = %d", nPatches));

    Int_t type = 0;

    for (Int_t i = 0; i < nPatches; i++) {
      AliDebug(2, Form("Processing patch %d", i));

      AliEMCALTriggerPatchInfo* patch = static_cast<AliEMCALTriggerPatchInfo*>(fTriggerPatches->At(i));
      if (!patch) continue;

      GetTriggerQA(fCentBin)->ProcessBkgPatch(patch);
    }

    GetTriggerQA(fCentBin)->ComputeBackground();

    FillEventQA();

    for (Int_t i = 0; i < nPatches; i++) {
      AliDebug(2, Form("Processing patch %d", i));

      AliEMCALTriggerPatchInfo* patch = static_cast<AliEMCALTriggerPatchInfo*>(fTriggerPatches->At(i));
      if (!patch) continue;

      GetTriggerQA(fCentBin)->ProcessPatch(patch);
    }
  }

  if (fCaloTriggers) {
    AliEMCALTriggerFastOR fastor;
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

      Int_t time = 0;
      Int_t nl0times(0);
      fCaloTriggers->GetNL0Times(nl0times);
       if(nl0times) {
         TArrayI l0times(nl0times);
         fCaloTriggers->GetL0Times(l0times.GetArray());
         for(int itime = 0; itime < nl0times; itime++){
           if(l0times[itime] >7 && l0times[itime] < 10){
             time = l0times[itime];
             break;
           }
         }
       }

      fastor.Initialize(L0amp, L1amp, globRow, globCol, time, fGeom);

      GetTriggerQA(fCentBin)->ProcessFastor(&fastor);
    }
  }
  GetTriggerQA(fCentBin)->EventCompleted();

  return kTRUE;
}

/**
 * Fill event QA THnSparse
 */
void AliEmcalTriggerQATask::FillEventQA()
{
  Double_t contents[20] = {0};

  Double_t DCalmedian[3] = {0};
  Double_t EMCalmedian[3] = {0};

  GetTriggerQA(fCentBin)->GetDCalMedian(DCalmedian);
  GetTriggerQA(fCentBin)->GetEMCalMedian(EMCalmedian);

  for (Int_t i = 0; i < fHistEventQA->GetNdimensions(); i++) {
    TString title(fHistEventQA->GetAxis(i)->GetTitle());
    if (title=="Centrality %")
      contents[i] = fCent;
    else if (title=="DCal median (offline)")
      contents[i] = DCalmedian[2];
    else if (title=="EMCal median (offline)")
      contents[i] = EMCalmedian[2];
    else if (title=="DCal median (recalc)")
      contents[i] = DCalmedian[1];
    else if (title=="EMCal median (recalc)")
      contents[i] = EMCalmedian[1];
    else
      AliWarning(Form("Unable to fill dimension %s!",title.Data()));
  }

  fHistEventQA->Fill(contents);
}

/**
 * Set number of ADC per bin in all the trigger QA
 * \param i number of ADC per bin.
 */
void AliEmcalTriggerQATask::SetADCperBin(Int_t n)
{
  fADCperBin = n;

  for (Int_t i = 0; i < fNcentBins; i++) {
    GetTriggerQA(i)->SetADCperBin(n);
  }
}

/**
 * Set background patch type in all the trigger QA
 * \param t background patch type
 */
void AliEmcalTriggerQATask::SetBkgPatchType(Int_t t)
{
  fBkgPatchType = t;

  for (Int_t i = 0; i < fNcentBins; i++) {
    GetTriggerQA(i)->SetBkgPatchType(t);
  }
}

/**
 * Set number of centrality bins and adjust fEMCALTriggerQA array accordingly
 * \param n number of centrality bins
 */
void AliEmcalTriggerQATask::SetNCentBins(Int_t n)
{
  if (n <= 0)  n = 0;

  // delete superflouos items in fEMCALTriggerQA
  for (Int_t i = n; i < fNcentBins; i++) {
    fEMCALTriggerQA->RemoveAt(i);
  }

  for (Int_t i = fNcentBins; i < n; i++) {
    TString qaName(Form("%s_AliEMCALTriggerQA_Cent%d", GetName(), i));
    if (fEMCALTriggerQA->At(0)) {
      AliEMCALTriggerQA* triggerQA = new AliEMCALTriggerQA(*(static_cast<AliEMCALTriggerQA*>(fEMCALTriggerQA->At(0))));
      triggerQA->SetName(qaName);
      triggerQA->SetTitle(qaName);
      fEMCALTriggerQA->AddAt(triggerQA, i);
    }
    else {
      fEMCALTriggerQA->AddAt(new AliEMCALTriggerQA(qaName), i);
    }
  }

  fNcentBins = n;
}
