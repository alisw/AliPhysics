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
#include "ADESDFriendUtils.h"

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
  , fADESDFriendUtils(NULL)
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

    delete fADESDFriendUtils;
    fADESDFriendUtils = NULL;
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
			 fBCRangeTail[0],
			 fBCRangeTail[1],
			 bc);
}

void AliAnalysisTaskADCalib::NotifyRun() {
  if (NULL == fADESDFriendUtils) {
    AliFatal("NULL == fADESDFriendUtils");
    return;
  }
  fADESDFriendUtils->Init(fCurrentRunNumber);
}

void AliAnalysisTaskADCalib::UserCreateOutputObjects() {
  fStatus = kOk; // pedantic

  // (1) set up the output list
  fList = new THashList;
  fList->SetOwner(kTRUE);

  // (2) populate the list with histograms
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

  // (3) set up AD ESD friend helper object
  fADESDFriendUtils = new ADESDFriendUtils;

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

  fADESDFriendUtils->Update(esdADfriend);

  TH2 *h = NULL;
  for (Int_t ch=0; ch<16; ++ch) { // offline channel number
    if (fADESDFriendUtils->IsPileUp(ch))
      continue;
    
    // (3) compute the charge in the tail
    Float_t tail = 0.0f;
    for (Int_t bc=fBCRangeTail[0], n=fBCRangeTail[1]; bc<=n; ++bc)
      tail += fADESDFriendUtils->GetADCPedSub(ch, bc);

    // (4) fill histograms with charge in BC vs. charge in the tail
    for (Int_t bc=fBCRangeExtrapolation[0], n=fBCRangeExtrapolation[1]; bc<=n; ++bc) {
      const Bool_t  integratorFlag = esdADfriend->GetIntegratorFlag(ch, bc);
      const TString histName = GetHistName(ch, bc, integratorFlag);
      const Bool_t  ok = FillHist(histName, tail, fADESDFriendUtils->GetADCPedSub(ch, bc));
      if (!ok)
	AliError(Form("FillHist failed for %s", histName.Data()));
    }
  }

  PostData(1, fList);
}

void AliAnalysisTaskADCalib::Terminate(Option_t* ) {
  // NOP
}

void AliAnalysisTaskADCalib::ProcessOutput(const Char_t  *fileName,
					   AliCDBStorage *cdbStorage,
					   Int_t runNumber) {
  fStatus = kOk;

  // if necessary clean up
  if (NULL != fList) {
    delete fList;
    fList = NULL;
  }

  // (1) open file
  TFile *f = TFile::Open(fileName);
  if (NULL == f || !f->IsOpen()) {
    AliError(Form("cannot open output file %s", fileName));
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
  TTree *tSat = MakeSaturationCalibObject();
  if (NULL == tSat) {
    AliError("Failed to make the AD saturation calibration OCDB object");
    fStatus = kMeasurementError;
    f->Close();
    return;
  }

  // (5) update the gain parameterization
  AliCDBEntry* gainOCDBObject = UpdateGainParameters(runNumber, tSat);
  if (NULL == gainOCDBObject) {
    AliError("Failed to generate the updated AD PM gain OCDB object");
    fStatus = kMeasurementError;
    f->Close();
    return;
  }

  f->Close();

  // (6) update list of histograms
  TString qaFileName = fileName;
  qaFileName.ReplaceAll(".root", "_ADQA.root");
  TFile *fSave = TFile::Open(qaFileName, "RECREATE");
  fList->Write("", TObject::kSingleKey | TObject::kWriteDelete);
  fSave->Write();
  fSave->Close();

  // (7) store OCDB objects
  // (7a) check for cdbStorage
  if (NULL == cdbStorage) {
    AliError("NULL == cdbStorage");
    fStatus = kStoreError;
    return;
  }

  // (7b) Creation of OCDB metadata and id
  AliCDBMetaData *md= new AliCDBMetaData();
  md->SetResponsible("C. Mayer");
  md->SetBeamPeriod(0);
  md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
  md->SetComment("AD Saturation");
  AliCDBId idSat ("AD/Calib/Saturation", runNumber, AliCDBRunRange::Infinity());
  AliCDBId idGain("AD/Calib/PMGains",    runNumber, AliCDBRunRange::Infinity());

  // (7c) put the objects into the OCDB storage
  fStatus = (cdbStorage->Put(tSat, idSat, md) &&
	     cdbStorage->Put(gainOCDBObject->GetObject(),
			     idGain,
			     gainOCDBObject->GetMetaData())
	     ? kOk
	     : kStoreError);
}

Bool_t AliAnalysisTaskADCalib::MakeExtrapolationFit(TH2 *h, TF1 *f, Int_t ch, Int_t bc, Double_t &xMax) {
  if (NULL == h || NULL == f) {
    xMax = -999.9f;
    return kFALSE;
  }

  // (1) set up the TF1 depending on BC
  switch (bc) {
  case 10:
    f->SetParameters(0, 8);
    f->SetParNames("offset", "slope");
    f->SetParLimits(1, 0.0, 40.0);
    break;
  case 11:
  case 12:
    f->SetParameters(0, 2, 0, 2);
    f->SetParNames("offset", "slope", "p_{0}", "power");
    f->SetParLimits(1, 0.0, 10.0);
    f->SetParLimits(2, 0.0,  1.0);
    f->SetParLimits(3, 1.1,  6.0);
    break;
  default:
    return kFALSE;
  }
  f->SetLineStyle(2);

  // (2) fit to the profile
  TProfile *h_pfx = h->ProfileX();
  h_pfx->SetDirectory(NULL);

  xMax = 50.0f;
  for (Int_t i=1, n=h_pfx->GetNbinsX(); i<n; ++i) {
    if (h_pfx->GetBinContent(i) <  10.0)
      continue;
    if (h_pfx->GetBinContent(i) > 750.0 ||
	(xMax != 50.0 && h_pfx->GetBinContent(i) <  10.0))
      break;
    xMax = h_pfx->GetXaxis()->GetBinUpEdge(i);
  }
  AliDebug(3, Form("ch=%02d bc=%2d xMax= %.1f", ch, bc, xMax));
  h_pfx->Fit(f, "WQ0", "", 0, xMax);

  // update xMax based on the fit function
  xMax = f->GetX(1024.0, 10.0, 590.0);

  // (3) cut outliers by adapting a TGutG to the fitted function
  const Int_t nInterations = 6;
  for (Int_t iteration=nInterations-1; iteration>=0; --iteration) {
    // (3a) make up a TCutG based on the last fit
    //      last iteraton: adapt dY to the slope
    const Double_t dX = 5.0;
    const Double_t dY = (iteration > 0 ? 300.0 : 25.0*(1024.0/xMax)*TMath::Sqrt(TMath::Power(1024.0/xMax, -2) + 1.0));

    const Int_t  iMax = Int_t(xMax/dX);
    TString cutName = TString::Format("%s_cutg_%d", h->GetName(), iteration);
    TCutG *cutg = new TCutG(cutName, 2*(2+iMax)+1);
    fList->Add(cutg);
    
    Int_t counter=0;
    for (Int_t i=-1; i<=iMax; ++i)
      cutg->SetPoint(counter++, dX*i, f->Eval(dX*i)-dY);
    for (Int_t i=iMax; i>=-1; --i)
      cutg->SetPoint(counter++, dX*i, f->Eval(dX*i)+dY);
    cutg->SetPoint(counter++, -dX, f->Eval(-dX)-dY);

    // (3b) get a new profile with TCutG
    h_pfx = h->ProfileX(Form("%s_pfxCut_%d", h->GetName(), iteration), 1, -1, "["+cutName+"]");
    fList->Add(h_pfx);
    h_pfx->SetDirectory(NULL);
    h_pfx->SetLineWidth(2);

    // (3c) fit the new profile and update xMax
    h_pfx->Fit(f, "WQ0", "", 0, xMax);
    xMax = f->GetX(1024.0, 10.0, 590.0);
  }
  return kTRUE;
}

Int_t AliAnalysisTaskADCalib::GetStatus() const {
  return fStatus;
}

TTree* AliAnalysisTaskADCalib::MakeSaturationCalibObject() {
  const Int_t gOffline2Online[16] = {
    15, 13, 11,  9, 14, 12, 10,  8,
     7,  5,  3,  1,  6,  4,  2,  0
  };  

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
  const Double_t meanQuantiles[2] = { 
    TMath::Mean(8, chargeQuantiles),
    TMath::Mean(8, chargeQuantiles+8)
  };

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
    AliInfo(Form("quantile[%2d] = %f f=%f", ch, chargeQuantiles[ch], chargeEqualizationFactor));

    for (Int_t bc=0; bc<21; ++bc) {
      TH2* h0 = dynamic_cast<TH2*>(fList->FindObject(GetHistName(ch, bc, 0)));
      TH2* h1 = dynamic_cast<TH2*>(fList->FindObject(GetHistName(ch, bc, 1)));
      switch (bc) {
	if (NULL != h0 && NULL != h1) {
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
	case 11:
	case 12: {
	  new (f_Int0[bc]) TF1(GetFcnName(ch, bc, 0), "[0] + [1]*x + [2]*abs(x)**[3]");
	  new (f_Int1[bc]) TF1(GetFcnName(ch, bc, 1), "[0] + [1]*x + [2]*abs(x)**[3]");
	  Bool_t fitOk = kTRUE;
	  Double_t thr[2] = { 0, 0 };
	  fitOk &= MakeExtrapolationFit(h0, static_cast<TF1*>(f_Int0[bc]), ch, bc, thr[0]);
	  fitOk &= MakeExtrapolationFit(h1, static_cast<TF1*>(f_Int1[bc]), ch, bc, thr[1]);
	  doExtrapolation[bc] = fitOk;
	  extrapolationThresholds[bc] = TMath::Min(thr[0], thr[1]);
	  break;
	}
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

AliCDBEntry* AliAnalysisTaskADCalib::UpdateGainParameters(Int_t runNumber, TTree* tSat)
{
  if (NULL == tSat) {
    AliError("NULL == tSat");
    return NULL;
  }
  // (0) get and set up the OCDB manager
  AliCDBManager *man = AliCDBManager::Instance();
  if (NULL == man) {
    AliFatal("CDB manager not found");
    return NULL;
  }

  if (!man->IsDefaultStorageSet())
    man->SetDefaultStorage("raw://");

  man->SetRun(runNumber);

  // (1a) Get the AD calibration data OCDB object
  AliCDBEntry *entry = man->Get("AD/Calib/Data");
  AliADCalibData* calibData = dynamic_cast<AliADCalibData*>(entry->GetObject());
  if (NULL == calibData) {
    AliFatal("No calibration data from calibration database");
    return NULL;
  }

  // (1b) Get the AD PM gain calibration object (to be updated)
  entry = man->Get("AD/Calib/PMGains");
  if (NULL == entry) {
    AliFatal("AD/Calib/Data not found");
    return NULL;
  }
  TH2* h2Gain = dynamic_cast<TH2*>(entry->GetObject());
  if (NULL == h2Gain) {
    AliFatal("No calibration data from calibration database");
    return NULL;
  }

  // (2a) compute the mean of the charge equalization factors per side
  const Int_t n = tSat->Draw("chargeEqualizationFactor", "", "GOFF");
  if (n != 16) {
    AliError(Form("Invalid Saturation OCDB object: n=%d != 16", n));
    return NULL;
  }
  Double_t *eqFactors = tSat->GetV1();
  const Float_t meanEq[2] = {
    TMath::Mean(8, eqFactors),
    TMath::Mean(8, eqFactors+8)
  };
  AliInfo(Form("meanEq=%.3f %.3f", meanEq[0], meanEq[1]));

  // (2b) extract mip,hv,a,b and compute the mean MIP per side
  Float_t mip[16] = { 0 };
  Float_t hv[16]  = { 0 };
  Float_t a[16]   = { 0 };
  Float_t b[16]   = { 0 };
  for (Int_t ch=0; ch<16; ++ch) {
    a[ch]   = h2Gain->GetBinContent(1+ch, 1);
    b[ch]   = h2Gain->GetBinContent(1+ch, 2);
    hv[ch]  = calibData->GetMeanHV(ch);
    mip[ch] = TMath::Power(hv[ch]/a[ch], b[ch]);
  }
  const Float_t meanMIP[2] = {
    TMath::Mean(8, mip),
    TMath::Mean(8, mip+8)
  };
  AliInfo(Form("meanMIP=%.3f %.3f", meanMIP[0], meanMIP[1]));
  
  // (3) update the a and b parameters in such a way to have the same MIP per side for the given HV
  for (Int_t ch=0; ch<16; ++ch) {
    a[ch] *= TMath::Power(meanEq[ch/8] * meanMIP[ch/8]/mip[ch], -1./b[ch]);
    AliInfo(Form("a=%.2f b=%.2f mip=%.3f mip'=%.3f",
		 a[ch], b[ch], mip[ch], TMath::Power(hv[ch]/a[ch], b[ch])));
    h2Gain->SetBinContent(1+ch, 1, a[ch]);
  }

  return entry;
}
