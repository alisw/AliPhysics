// -*- C++ -*-
// $Id$

/**************************************************************************
 * Author: C. Mayer                                                       *
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

#include <TCutG.h>
#include <TFile.h>
#include <THashList.h>
#include <TH2.h>
#include <TList.h>
#include <TF1.h>
#include <TProfile.h>
#include <TString.h>
#include <TTree.h>

#include "AliLog.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskADCalib.h"

#include "AliESDEvent.h"
#include "AliESDfriend.h"
#include "AliESDAD.h"
#include "AliESDADfriend.h"
#include "AliADCalibData.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliCDBStorage.h"
#include "AliCDBMetaData.h"
#include "AliCDBId.h"

ClassImp(AliAnalysisTaskADCalib);

AliAnalysisTaskADCalib::AliAnalysisTaskADCalib(const char *name)
  : AliAnalysisTaskSE(name)
  , fCalibData(NULL)
  , fList(NULL)
  , fStatus(kOk) {
  fBCRangeTail[0] = 14;
  fBCRangeTail[1] = 20;

  fBCRangeExtrapolation[0] = 10;
  fBCRangeExtrapolation[1] = 13;

  DefineOutput(1, TList::Class());
}

AliAnalysisTaskADCalib::~AliAnalysisTaskADCalib() {
  if (AliAnalysisManager::GetAnalysisManager()->GetAnalysisType() == AliAnalysisManager::kProofAnalysis) {
    delete fList;
    fList = NULL;
  }
}

Bool_t AliAnalysisTaskADCalib::FillHist(TString name, Double_t x, Double_t y) {
  if (NULL == fList) return kFALSE;
  TH2* h = dynamic_cast<TH2*>(fList->FindObject(name.Data()));
  if (NULL == h) return kFALSE;
  h->Fill(x, y);
  return kTRUE;
}

TString AliAnalysisTaskADCalib::GetHistName(Int_t  ch, // offline channel,
					    Int_t  bc, // bc
					    Bool_t integrator) const { // integrator
  return TString::Format("hCh%02d_bc%02d_int%d",
			 ch, bc, integrator);
}
TString AliAnalysisTaskADCalib::GetFcnName(Int_t  ch, // offline channel,
					   Int_t  bc, // bc
					   Bool_t integrator) const { // integrator
  return TString::Format("f_Ch%02d_BC%02d_int%d",
			 ch, bc, integrator);
}
TString AliAnalysisTaskADCalib::GetHistTitle(Int_t  ch, // offline channel,
					     Int_t  bc, // bc
					     Bool_t integrator) const { // integrator
  return TString::Format("chOff=%02d BC=%02d int=%d;charge in tail [%d..%d] (ADC);charge in BC=%02d",
			 ch, bc, integrator,
			 fBCRangeExtrapolation[0],
			 fBCRangeExtrapolation[1],
			 bc);
}

void AliAnalysisTaskADCalib::NotifyRun() {
  // (1) get the AD Calibration OCDB object
  AliCDBManager *man = AliCDBManager::Instance();
  if (NULL == man) {
    AliFatal("CDB manager not found");
    return;
  }

  if (!man->IsDefaultStorageSet())
    man->SetDefaultStorage("raw://");

  man->SetRun(fCurrentRunNumber);

  AliCDBEntry *entry = man->Get("AD/Calib/Data");
  if (NULL == entry) {
    AliFatal("AD/Calib/Data not found");
    return;
  }
  fCalibData = dynamic_cast<AliADCalibData*>(entry->GetObject());
  if (NULL == fCalibData) {
    AliFatal("No calibration data from calibration database");
    return;
  }
}

void AliAnalysisTaskADCalib::UserCreateOutputObjects() {

  fStatus = kOk; // pedantic

  // (1a) set up the output list
  fList = new THashList;
  fList->SetOwner(kTRUE);

  // (1b) populate the list with histograms
  TH2 *h = NULL;
  for (Int_t ch=0; ch<16; ++ch) { // offline channel number
    for (Int_t bc=fBCRangeExtrapolation[0], n=fBCRangeExtrapolation[1]; bc<=n; ++bc) { // bc
      for (Int_t integrator=0; integrator<2; ++integrator) { // integrator flag for bc
	h = new TH2D(GetHistName (ch, bc, integrator),
		     GetHistTitle(ch, bc, integrator),
		     301, -1.0,  601.0,
		     513, -2.0, 4098.0);
	fList->Add(h);
      }
    }
  }
  PostData(1, fList);
}

void AliAnalysisTaskADCalib::UserExec(Option_t* ) {
  AliESDEvent* esdEvent = dynamic_cast<AliESDEvent*>(InputEvent());
  if (NULL == esdEvent) {
    AliError("NULL == esdEvent");
    return;
  }
  AliESDAD* esdAD = esdEvent->GetADData();
  if (NULL == esdAD) {
    AliError("NULL == esdAD");
    return;
  }

  AliESDfriend *esdFriend = esdEvent->FindFriend();  
  if (NULL == esdFriend) {
    AliError("NULL == esdFriend");
    return;
  }
  AliESDADfriend* esdADfriend = esdFriend->GetADfriend();
  if (NULL == esdADfriend) {
    AliError("NULL == esdADfriend");
    return;
  }

  TH2 *h = NULL;
  for (Int_t ch=0; ch<16; ++ch) { // offline channel number
    // (1) compute pedestal subtracted ADC values per BC
    Float_t adcPedSub[21] = { };
    for (Int_t bc=0; bc<21; ++bc) {
      adcPedSub[bc]  = Float_t(esdADfriend->GetPedestal(ch, bc));
      adcPedSub[bc] -= fCalibData->GetPedestal(ch + 16*esdADfriend->GetIntegratorFlag(ch, bc));
    }

    // (2) reject events with secondary peaks in the charge time series
    const Float_t threshold = 20.0f;
    Bool_t isPileUp = kFALSE;
    for (Int_t bc=13; bc<20 && !isPileUp; ++bc)  
      isPileUp |= (adcPedSub[bc+1] > adcPedSub[bc] + threshold);
    
    if (isPileUp)
      continue;

    // (3) compute the charge in the tail
    Float_t tail = 0.0f;
    for (Int_t bc=fBCRangeTail[0], n=fBCRangeTail[1]; bc<=n; ++bc) {
      tail += adcPedSub[bc];
    }

    // (4) fill histograms with charge in BC vs. charge in the tail
    for (Int_t bc=fBCRangeExtrapolation[0], n=fBCRangeExtrapolation[1]; bc<=n; ++bc) {
      const Bool_t  integratorFlag = esdADfriend->GetIntegratorFlag(ch, bc);
      const TString histName = GetHistName(ch, bc, integratorFlag);
      const Bool_t  ok = FillHist(histName, tail, adcPedSub[bc]);
      if (!ok)
	AliError(Form("FillHist failed for %s", histName.Data()));
    }
  }

  PostData(1, fList);
}

void AliAnalysisTaskADCalib::Terminate(Option_t* ) {
  // NOP
}

void AliAnalysisTaskADCalib::ProcessOutput(const Char_t  *filename,
					   AliCDBStorage *cdbStorage,
					   Int_t runNumber) {
  fStatus = kOk;

  if (NULL != fList) {
    delete fList;
    fList = NULL;
  }

  // (1) open file
  TFile *f = TFile::Open(filename);
  if (NULL == f || !f->IsOpen()) {
    AliError(Form("cannot open output file %s", filename));
    fStatus = kInputError;
    return;
  }

  // (2) cd
  if (!f->cd("ADCalib")) {
    fStatus = kDataError;
    return;
  }

  // (3) get list of histograms
  fList = dynamic_cast<TList*>(gDirectory->Get("ADCalibListHist"));
  if (NULL == fList) {
    fStatus = kDataError;
    f->Close();
    return;
  }

  // (4) make the OCDB calibration TTree
  TTree *t = MakeCalibObject();
  if (NULL == t) {
    fStatus = kMeasurementError;
    f->Close();
    return;
  }

  f->Close();

  // (5) update list of histograms
  TFile *fSave = TFile::Open("CalibObjectsQA_AD.root", "RECREATE");
  fList->Write("", TObject::kSingleKey | TObject::kWriteDelete);
  fSave->Write();
  fSave->Close();

  // (6) Creation of OCDB metadata and id
  AliCDBMetaData *md= new AliCDBMetaData();
  md->SetResponsible("C. Mayer");
  md->SetBeamPeriod(0);
  md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
  md->SetComment("AD Saturation");
  AliCDBId id("AD/Calib/Saturation", runNumber, AliCDBRunRange::Infinity());

  // (7) store the AD saturation correction OCDB object
  fStatus = (cdbStorage->Put(t, id, md)
	     ? kOk
	     : kStoreError);
}

Bool_t AliAnalysisTaskADCalib::MakeExtrapolationFit(TH2 *h, TF1 *f, Int_t ch, Int_t bc, Double_t &xMax) {
  // (1) set up the TF1 depending on BC
  switch (bc) {
  case 10:
    f->SetParameters(0, 8);
    f->SetParNames("offset", "slope");
    f->SetParLimits(1, 0.0, 20.0);
    xMax = 120;
    break;
  case 11:
  case 12:
    f->SetParameters(0, 2, 0, 2);
    f->SetParNames("offset", "slope", "p_{0}", "power");
    f->SetParLimits(1, 0.0, 10.0);
    f->SetParLimits(2, 0.0,  1.0);
    f->SetParLimits(3, 1.1,  6.0);
    xMax = 300 + 200*(ch<=8);
    break;
  default:
    return kFALSE;
  }
  f->SetLineStyle(2);

  // (2) fit to the profile
  TProfile *h_pfx = h->ProfileX();
  h_pfx->SetDirectory(NULL);
  h_pfx->Fit(f, "", "", 0, xMax);

  // new xMax is for saturation
  xMax = f->GetX(1024.0, 10.0, 590.0);

  // (3) cut outliers
  // (3a) make up a TCutG
  const Double_t dX =   5.0;
  const Double_t dY = 120.0;
  const Int_t  iMax = Int_t(xMax/dX);
  TString cutName = TString::Format("%s_cutg", h->GetName());
  TCutG *cutg = new TCutG(cutName, 2*(2+iMax)+1);
  fList->Add(cutg);

  Int_t counter=0;
  for (Int_t i=-1; i<=iMax; ++i)
    cutg->SetPoint(counter++, dX*i, f->Eval(dX*i)-dY);
  for (Int_t i=iMax; i>=-1; --i)
    cutg->SetPoint(counter++, dX*i, f->Eval(dX*i)+dY);
  cutg->SetPoint(counter++, -dX, f->Eval(-dX)-dY);

  // (3b) get a new profile with TCutG
  h_pfx = h->ProfileX(Form("%s_pfxCut", h->GetName()), 1, -1, "["+cutName+"]");
  fList->Add(h_pfx);
  h_pfx->SetDirectory(NULL);
  h_pfx->SetLineWidth(2);

  // (4) fit the new profile and update xMax
  h_pfx->Fit(f, "", "", 0, xMax);
  xMax = f->GetX(1024.0, 10.0, 590.0);

  return kTRUE;
}

Int_t AliAnalysisTaskADCalib::GetStatus() const {
  return fStatus;
}

const Int_t gOffline2Online[16] = {
  15,13,11,9,14,12,10,8,
  7,5,3,1,6,4,2,0
};

TTree* AliAnalysisTaskADCalib::MakeCalibObject() {
  if (NULL == fList)
    return NULL;

  // (0) set up the TTree
  TTree *t = new TTree;

  Int_t chOffline=0, chOnline=0;
  t->Branch("chOffline", &chOffline);
  t->Branch("chOnline",  &chOnline);

  TClonesArray f_Int0("TF1", 21);
  TClonesArray f_Int1("TF1", 21);
  t->Branch("f_Int0", &f_Int0, 32000, 0);
  t->Branch("f_Int1", &f_Int1, 32000, 0);
  f_Int0.BypassStreamer();
  f_Int1.BypassStreamer();

  Float_t extrapolationThresholds[21];
  Bool_t  doExtrapolation[21];
  t->Branch("extrapolationThresholds", &extrapolationThresholds, "val[21]/F");
  t->Branch("doExtrapolation",         &doExtrapolation,         "val[21]/O");

  // (1) compute charge quantiles on the tail charge
  Double_t chargeQuantiles[16];
  for (Int_t ch=0; ch<16; ++ch) { // offline channel number
    TH2 *h0 = dynamic_cast<TH2*>(fList->FindObject(GetHistName(ch, 10, 0)));
    Double_t qx = 0.99; // 99% quantile
    TH1* h1 = h0->ProjectionX();
    h1->SetDirectory(NULL);
    h1->GetQuantiles(1, chargeQuantiles+ch, &qx);
    delete h1;
  }
  Double_t meanQuantiles[2] = { 0,0 };
  meanQuantiles[0] = TMath::Mean(8, chargeQuantiles);
  meanQuantiles[1] = TMath::Mean(8, chargeQuantiles+8);

  // (2) compute charge equalization factors
  Float_t chargeEqualizationFactor = 1.0f;
  t->Branch("chargeEqualizationFactor", &chargeEqualizationFactor);
  
  // (3) fill TTree
  for (Int_t ch=0; ch<16; ++ch) { // offline channel number
    f_Int0.Clear();
    f_Int1.Clear();

    chOffline = ch;
    chOnline  = gOffline2Online[ch];

    chargeEqualizationFactor = meanQuantiles[ch/8]/chargeQuantiles[ch];
    Printf("quantile[%2d] = %f f=%f", ch, chargeQuantiles[ch], chargeEqualizationFactor);

    for (Int_t bc=0; bc<21; ++bc) {
      TH2* h0 = dynamic_cast<TH2*>(fList->FindObject(GetHistName(ch, bc, 0)));
      TH2* h1 = dynamic_cast<TH2*>(fList->FindObject(GetHistName(ch, bc, 1)));
      switch (bc) {
      case 10: {
	new (f_Int0[bc]) TF1(GetFcnName(ch, bc, 0), "[0] + [1]*x");
	new (f_Int1[bc]) TF1(GetFcnName(ch, bc, 1), "[0] + [1]*x");
	Bool_t fitOk = kTRUE;
	Double_t thr[2] = { 0, 0 };
	fitOk &= MakeExtrapolationFit(h0, static_cast<TF1*>(f_Int0[bc]), ch, bc, thr[0]);
	fitOk &= MakeExtrapolationFit(h1, static_cast<TF1*>(f_Int1[bc]), ch, bc, thr[1]);
	doExtrapolation[bc] = fitOk;
	extrapolationThresholds[bc] = TMath::Min(thr[0], thr[1]);
	break;
      }
      case 11: {
	new (f_Int0[bc]) TF1(GetFcnName(ch, bc, 0), "[0] + [1]*x + (x>0 ? [2]*x**[3] : 0)");
	new (f_Int1[bc]) TF1(GetFcnName(ch, bc, 1), "[0] + [1]*x + (x>0 ? [2]*x**[3] : 0)");
	Bool_t fitOk = kTRUE;
	Double_t thr[2] = { 0, 0 };
	fitOk &= MakeExtrapolationFit(h0, static_cast<TF1*>(f_Int0[bc]), ch, bc, thr[0]);
	fitOk &= MakeExtrapolationFit(h1, static_cast<TF1*>(f_Int1[bc]), ch, bc, thr[1]);
	doExtrapolation[bc] = fitOk;
	extrapolationThresholds[bc] = TMath::Min(thr[0], thr[1]);
	break;
      }
      default:
	doExtrapolation[bc]         = kFALSE;
	extrapolationThresholds[bc] = -999.9f;
	new (f_Int0[bc]) TF1(GetFcnName(ch, bc, 0), "x");
	new (f_Int1[bc]) TF1(GetFcnName(ch, bc, 1), "x");
      }
    }
    t->Fill();
  } // next channel

  return t;
}
