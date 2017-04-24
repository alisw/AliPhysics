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

#include <TArrayI.h>
#include <TCutG.h>
#include <TFile.h>
#include <THashList.h>
#include <TH2.h>
#include <THnSparse.h>
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
#include "AliADRecoParam.h"

ClassImp(AliAnalysisTaskADCalib);

AliAnalysisTaskADCalib::AliAnalysisTaskADCalib(const char *name)
  : AliAnalysisTaskSE(name)
  , fADESDFriendUtils(NULL)
  , fList(NULL)
  , fStatus(kOk) {
  fBCRangeTail[0] = 14;
  fBCRangeTail[1] = 20;

  fBCRangeExtrapolation[0] =  9;
  fBCRangeExtrapolation[1] = 15;

  fTimeResolution[0] = 25./256.;
  fTimeResolution[1] = 25./256.;

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
Bool_t AliAnalysisTaskADCalib::FillHist(TString name, Double_t x, Double_t y, Double_t z) {
  if (NULL == fList) return kFALSE;
  THnSparse* h = dynamic_cast<THnSparse*>(fList->FindObject(name.Data()));
  if (NULL == h) return kFALSE;
  const Double_t xyz[3] = { x,y,z };
  h->Fill(xyz);
  return kTRUE;
}

TString AliAnalysisTaskADCalib::GetHistName(Int_t  ch,                 // offline channel
					    Int_t  bc,                 // bc
					    Bool_t integrator) const { // integrator
  return TString::Format("hCh%02d_bc%02d_int%d", ch, bc, integrator);
}
TString AliAnalysisTaskADCalib::GetHistName(Int_t  ch,                 // offline channel
					    Bool_t integrator) const { // integrator
  return TString::Format("hCh%02d_int%d", ch, integrator);
}
TString AliAnalysisTaskADCalib::GetFcnName(Int_t  ch,                 // offline channel
					   Int_t  bc,                 // bc
					   Bool_t integrator) const { // integrator
  return TString::Format("f_Ch%02d_BC%02d_int%d", ch, bc, integrator);
}
TString AliAnalysisTaskADCalib::GetHistTitle(Int_t  ch,                 // offline channel
					     Int_t  bc,                 // bc
					     Bool_t integrator) const { // integrator
  return TString::Format("chOff=%02d BC=%02d int=%d;charge in tail [%d..%d] (ADC);charge in BC=%02d",
			 ch, bc, integrator, fBCRangeTail[0], fBCRangeTail[1], bc);
}
TString AliAnalysisTaskADCalib::GetHistTitle(Int_t  ch,                 // offline channel,
					     Bool_t integrator) const { // integrator
  return TString::Format("chOff=%02d int(BC=10)=%d;charge in tail [%d..%d] (ADC);time Ch%02d (ns);time Ch%02d (ns)",
			 ch, integrator, fBCRangeTail[0], fBCRangeTail[1], ch, ch+4);
}

void AliAnalysisTaskADCalib::NotifyRun() {
  // (1) set up ADESDFriendUtils
  if (NULL == fADESDFriendUtils) {
    AliFatal("NULL == fADESDFriendUtils");
    return;
  }
  fADESDFriendUtils->Init(fCurrentRunNumber);

  // get the time resolution from OCDB
  AliCDBManager *man = AliCDBManager::Instance();
  if (NULL == man) {
    AliFatal("CDB manager not found");
    return;
  }
  if (!man->IsDefaultStorageSet())
    man->SetDefaultStorage("raw://");

  man->SetRun(fCurrentRunNumber);

  AliCDBEntry *entry = man->Get("AD/Calib/Data");
  AliADCalibData* calibData = dynamic_cast<AliADCalibData*>(entry->GetObject());
  if (NULL == calibData) {
    AliFatal("No calibration data from calibration database");
    return;
  }

  fTimeResolution[0] = calibData->GetTimeResolution(0);
  fTimeResolution[1] = calibData->GetTimeResolution(1);
  AliInfo(Form("timeResolution: %f %f", fTimeResolution[0], fTimeResolution[1]));

  // get the event specie from AliESDEvent
  AliESDEvent* esdEvent = dynamic_cast<AliESDEvent*>(InputEvent());
  if (NULL == esdEvent) {
    AliFatal("NULL == esdEvent");
    return;
  }
  entry = man->Get("AD/Calib/RecoParam");
  const TObjArray *recoParamArray =  dynamic_cast<const TObjArray*>(entry->GetObject());
  if (NULL == recoParamArray) {
    AliFatal("NULL == recoParamArray");
    return;
  }
  const AliADRecoParam *recoParam = NULL;
  for (Int_t i=0, n=recoParamArray->GetEntries(); i<n; ++i) {
    recoParam = dynamic_cast<const AliADRecoParam*>(recoParamArray->At(i));
    if (NULL == recoParam)
      continue;
    if ((recoParam->GetEventSpecie() & esdEvent->GetEventSpecie()) != 0)
      break;
  }
  if (NULL == recoParam) {
    AliFatal("NULL == recoParam");
    return;
  }
  fBCRangeTail[0] = recoParam->GetTailBegin();
  fBCRangeTail[1] = recoParam->GetTailEnd();
  AliInfo(Form("BCRangeTail: [%2d,%2d]", fBCRangeTail[0], fBCRangeTail[1]));
}

void AliAnalysisTaskADCalib::UserCreateOutputObjects() {
  fStatus = kOk; // pedantic

  // (1) set up the output list
  fList = new THashList;
  fList->SetOwner(kTRUE);

  // (2) populate the list with histograms
  Double_t tMin[2] = {   0.0,   0.0 };
  Double_t tMax[2] = { 400.0, 400.0 };
  for (Int_t ch=0; ch<16; ++ch) { // offline channel number
    const Int_t    side = ch/8;
    const Double_t dt   = fTimeResolution[side];
    tMin[side] = TMath::Nint(tMin[side]/dt/4)*dt*4;
    tMax[side] = TMath::Nint(tMax[side]/dt/4)*dt*4;
    const Int_t tBins = TMath::Nint((tMax[side]-tMin[side])/(dt*4));
    for (Int_t integrator=0; integrator<2; ++integrator) { // integrator flag for bc
      for (Int_t bc=fBCRangeExtrapolation[0], n=fBCRangeExtrapolation[1]; bc<=n; ++bc) { // bc
 	fList->Add(new TH2D(GetHistName (ch, bc, integrator),
 			    GetHistTitle(ch, bc, integrator),
 			    302, -1.0,  601.0,
  			    129, -4.0, 1028.0));
      }
      const Int_t   nBins[3] = { 5000,      tBins, tBins      };
      const Double_t xMin[3] = {   -1, tMin[side], tMin[side] };
      const Double_t xMax[3] = { 9999, tMax[side], tMax[side] };
      fList->Add(new THnSparseI(GetHistName (ch, integrator),
				GetHistTitle(ch, integrator),
				3, nBins, xMin, xMax));
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
  const AliESDAD* esdAD = esdEvent->GetADData();
  if (NULL == esdAD) {
    AliError("NULL == esdAD");
    return;
  }

  AliESDfriend *esdFriend = esdEvent->FindFriend();
  if (NULL == esdFriend) {
    AliError("NULL == esdFriend");
    return;
  }
  const AliESDADfriend* esdADfriend = esdFriend->GetADfriend();
  if (NULL == esdADfriend) {
    AliError("NULL == esdADfriend");
    return;
  }

  fADESDFriendUtils->Update(esdADfriend);

  Float_t tail[16]      = { 0 };      // tail charges
  Float_t time[16]      = { 0 };      // HTPDC time (ns) from esdFriend
  Bool_t  timeValid[16] = { kFALSE }; // true if there is a time measurement and not pileup

  for (Int_t ch=0; ch<16; ++ch) { // offline channel number
    tail[ch]      = 0.0f;
    time[ch]      = 0.0f;
    timeValid[ch] = kFALSE;

    if (fADESDFriendUtils->IsPileUp(ch))
      continue;

    time[ch]      =  esdADfriend->GetTime(ch);
    timeValid[ch] = (esdADfriend->GetWidth(ch) > 1.0f);

    // (3) compute the charge in the tail
    for (Int_t bc=fBCRangeTail[0], n=fBCRangeTail[1]; bc<=n; ++bc)
      tail[ch] += fADESDFriendUtils->GetADCPedSub(ch, bc);

    // (4) fill 2D histograms with charge in BC vs. tail charge
    for (Int_t bc=fBCRangeExtrapolation[0], n=fBCRangeExtrapolation[1]; bc<=n; ++bc) {
      const Bool_t  integratorFlag = esdADfriend->GetIntegratorFlag(ch, bc);
      const TString histName = GetHistName(ch, bc, integratorFlag);
      const Bool_t  ok = FillHist(histName, tail[ch], fADESDFriendUtils->GetADCPedSub(ch, bc));
      if (!ok)
	AliError(Form("FillHist failed for %s", histName.Data()));
    }
  }
  // (5) fill 3D histograms with time2 vs. time1 vs. tail charge
  for (Int_t side=0; side<2; ++side) {
    for (Int_t i=0; i<4; ++i) {
      const Int_t ch[2] = { 8*side + i, 8*side + i + 4 }; // offline channel numbers of adjacent pads
      if (timeValid[ch[0]] && timeValid[ch[1]]) {
	for (Int_t j=0; j<2; ++j) {
	  const Bool_t  integratorFlag = esdADfriend->GetIntegratorFlag(ch[j], 10);
	  const TString histName = GetHistName(ch[j], integratorFlag);
	  const Bool_t  ok = FillHist(histName, tail[ch[j]], time[ch[0]], time[ch[1]]);
	  if (!ok)
	    AliError(Form("FillHist failed for %s", histName.Data()));
	}
      }
    }
  }
  PostData(1, fList);
}

void AliAnalysisTaskADCalib::Terminate(Option_t* ) {
  // NOP
}

void AliAnalysisTaskADCalib::ProcessOutput(const Char_t  *fileName,
					   AliCDBStorage *cdbStorage,
					   Int_t runNumber,
					   Int_t runNumberEnd) {
  fStatus = kOk;

  // if necessary clean up
  if (NULL != fList) {
    delete fList;
    fList = NULL;
  }

  // (0) get and set up the OCDB manager
  AliCDBManager *man = AliCDBManager::Instance();
  if (NULL == man) {
    fStatus = kInputError;
    AliFatal("CDB manager not found");
    return;
  }

  if (!man->IsDefaultStorageSet())
    man->SetDefaultStorage("raw://");

  man->SetRun(runNumber);

  // (1a) Get the AD calibration data OCDB object
  AliCDBEntry *entry = man->Get("AD/Calib/Data");
  AliADCalibData* calibData = dynamic_cast<AliADCalibData*>(entry->GetObject());
  if (NULL == calibData) {
    fStatus = kInputError;
    AliFatal("No calibration data from calibration database");
    return;
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
  TTree *tSat = MakeSaturationCalibObject(calibData);
  if (NULL == tSat) {
    AliError("Failed to make the AD saturation calibration OCDB object");
    fStatus = kMeasurementError;
    f->Close();
    return;
  }

  // (5) update the gain parameterization
  AliCDBEntry* gainOCDBObject = UpdateGainParameters(runNumber, tSat, calibData);
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
  fSave->mkdir("ADCalib");
  fSave->cd("ADCalib");
  fList->Write("ADCalibListHist", TObject::kSingleKey | TObject::kWriteDelete);
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
  AliCDBId idSat ("AD/Calib/Saturation", runNumber, (runNumberEnd == -1 ? runNumber : runNumberEnd));
  AliCDBId idGain("AD/Calib/PMGains",    runNumber, runNumber);

  // (7c) put the objects into the OCDB storage
  fStatus = (cdbStorage->Put(tSat, idSat, md) &&
	     cdbStorage->Put(gainOCDBObject->GetObject(),
			     idGain,
			     gainOCDBObject->GetMetaData())
	     ? kOk
	     : kStoreError);
}

TGraph* AliAnalysisTaskADCalib::MakeGraphSlope(TH2 *h, Double_t &s, const TString& name) const {
  // make up a TGraph containing sum of the histogram
  // along a line through (0,0) with a given slope vs. this slope
  TAxis *ax = h->GetXaxis();
  TAxis *ay = h->GetYaxis();

  const Int_t nx = h->GetNbinsX();
  const Int_t ny = h->GetNbinsY();

  TGraph *g = new TGraph;
  g->SetName(name);
  for (Int_t i=2; i<=nx; ++i) {
    const Double_t x     = ax->GetBinCenter(i);
    const Double_t slope = ay->GetBinCenter(ny)/x;
    Double_t sum=0;
    for (Int_t j=1; j<=ny; ++j) {
      const Double_t y = ay->GetBinCenter(j);
      const Double_t x = y/slope;
      if (x<10 || y<5) continue;
      sum += h->Interpolate(x, y);
    }
    g->SetPoint(g->GetN(), slope, sum);
  }
  for (Int_t i=2; i<=ny; ++i) {
    const Double_t y     = ay->GetBinCenter(i);
    const Double_t slope = y/ax->GetBinCenter(nx);
    Double_t sum=0;
    for (Int_t j=1; j<=nx; ++j) {
      const Double_t x = ax->GetBinCenter(j);
      const Double_t y = slope*x;
      if (x<10 || y<5) continue;
      sum += h->Interpolate(x, y);
    }
    g->SetPoint(g->GetN(), slope, sum);
  }
  TArrayI idx(g->GetN());
  TMath::Sort(g->GetN(), g->GetY(), idx.GetArray());
  s= g->GetX()[idx[0]];
  AliDebugClassF(5, "slope = %f", s);
  return g;
}

TH2* AliAnalysisTaskADCalib::RemoveHorizontalLines(TH2* h, Int_t dx) const {
  // removes horizontal lines which are due to ADC saturation and can bias the fit
  for (Int_t j=2; j<h->GetNbinsY(); ++j) {
    for (Int_t i=dx+1; i<=h->GetNbinsX()-dx; ++i) {
      Int_t sum = 0;
      for (Int_t k=i-dx; k<i+dx; ++k) {
	sum += Bool_t(h->GetBinContent(k,j-1));
	sum += Bool_t(h->GetBinContent(k,j+1));
      }
      if (!sum) {
	for (Int_t k=i-dx; k<i+dx; ++k)
	  h->SetBinContent(k,j, 0);
      }
    }
  }
  return h;
}

Bool_t AliAnalysisTaskADCalib::MakeExtrapolationFit(TH2 *h, TF1 *f, Int_t ch, Int_t bc, Double_t &xMax) {
  if (NULL == h || NULL == f || (ch < 8 && bc == 9)) {
    xMax = -999.9f;
    return kFALSE;
  }

  // (1a) do not fit on (nearly) empty histograms
  if (h->GetEntries() < 10)
    return kFALSE;

  // (1b) initial estimates for the slope, xMax, yMax
  Double_t slope=0;
  fList->Add(MakeGraphSlope(RemoveHorizontalLines(h, 10), // removes activity due to ADC saturation
			    slope, TString::Format("%s_slope", h->GetName())));
  TF1 f0("f0", "[0]*x", 0, 600);
  f0.SetParameter(0, slope);
  xMax = f0.GetX(1024.0, 10.0, 590.0);
  Double_t yMax = f0.Eval(xMax);
  AliDebugF(5, "Ch%02d bc=%2d slope=%.1f", ch, bc, slope);

  // (1c) set up the TF1 depending on the BC
  switch (bc) {
  case  9:
    f->SetParameters(1024.0, slope/1024.0, 0.0, 0.0);
    break;
  case 10:
    f->SetParameters(0, slope);
    f->SetParNames("offset", "slope");
    f->SetParLimits(0, -10, 10);
    f->SetParLimits(1, 0.0, 40.0);
    break;
  case 11:
  case 12:
  case 13:
  case 14:
  case 15:
    f->SetParameters(0, slope, 0.1, 1.5);
    f->SetParNames("offset", "slope", "p_{0}", "power");
    f->SetParLimits(0,-20.0,20.0);
    f->SetParLimits(1, 0.0, 10.0);
    f->SetParLimits(2, 0.0,  1.0);
    f->SetParLimits(3, 1.1,  4.0);
    if (bc >= 13 && slope < 0.6) {
      f->FixParameter(2, 0.0);
      f->FixParameter(3, 1.1);
    }
    break;
  default:
    return kFALSE;
  }
  f->SetLineStyle(2);

  // (2a) initial profile (without TCutG applied)
  TProfile *h_pfx = h->ProfileX(Form("%s_pfxNoCut", h->GetName()), 1, -1);
  fList->Add(h_pfx);
  h_pfx->SetDirectory(NULL);
  h_pfx->SetLineWidth(2);

  // (2b) fit the profile
  f->FixParameter(1, f->GetParameter(1));
  if (bc == 9) {
    Double_t q, p=0.9999;
    h->ProjectionX()->GetQuantiles(1, &q, &p);
    h_pfx->Fit(f, "Q0", "", 0, q);
  } else {
    h_pfx->Fit(f, "Q0", "", 0, xMax);
  }
  f->ReleaseParameter(1);

  // (3) cut outliers+pile-up by adapting iteratively a TGutG to the fitted function
  const Int_t nInterations = 6;
  for (Int_t iteration=nInterations-1; iteration>=0; --iteration) {
    // (3a) make up a TCutG based on the last fit
    //      last iteraton: adapt dY to the slope
    const Double_t dX = 5.0;
    //      a line with distance d to a given linear function with slope m has an offset of d*sqrt(1+m**2)
    const Double_t dY = 15*((bc>=11 || bc==9) ? (iteration+1) : 1)*TMath::Sqrt(1.0+TMath::Power(yMax/xMax, 2));
    const Int_t  iMax = Int_t(0.5+xMax/dX);
    TString cutName = TString::Format("%s_cutg_%d", h->GetName(), iteration);
    TCutG *cutg = new TCutG(cutName, 2*(2+iMax));
    fList->Add(cutg);
    Int_t counter=0;
    for (Int_t i=-1; i<=iMax; ++i)
      cutg->SetPoint(counter++, dX*i, f->Eval(dX*i)-dY);
    for (Int_t i=iMax; i>=-1; --i)
      cutg->SetPoint(counter++, dX*i, f->Eval(dX*i)+dY);

    // (3b) make a new profile with the TCutG
    h_pfx = h->ProfileX(Form("%s_pfxCut_%d", h->GetName(), iteration), 1, -1, "["+cutName+"]");
    fList->Add(h_pfx);
    h_pfx->SetDirectory(NULL);
    h_pfx->SetLineWidth(2);

    // (3c) fit the new profile and update xMax, yMax
    h_pfx->Fit(f, "WQ0", "", 0, xMax);
    xMax = f->GetX(1024.0, 10.0, 590.0);
    yMax = f->Eval(xMax);
  }
  return kTRUE;
}

Int_t AliAnalysisTaskADCalib::GetStatus() const {
  return fStatus;
}

TTree* AliAnalysisTaskADCalib::MakeSaturationCalibObject(AliADCalibData* calibData) {
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
  f_Int0.BypassStreamer(kFALSE);
  f_Int1.BypassStreamer(kFALSE);

  Float_t extrapolationThresholds[21];
  Bool_t  doExtrapolation[21];
  t->Branch("extrapolationThresholds", &extrapolationThresholds, "val[21]/F");
  t->Branch("doExtrapolation",         &doExtrapolation,         "val[21]/O");

  // (1) compute charge quantiles on the tail charge
  Double_t chargeQuantiles[16];
  for (Int_t ch=0; ch<16; ++ch) { // offline channel number
    TH2 *h0Int0 = dynamic_cast<TH2*>(fList->FindObject(GetHistName(ch, 10, 0)));
    TH2 *h0Int1 = dynamic_cast<TH2*>(fList->FindObject(GetHistName(ch, 10, 1)));
    TH2* hSum   = dynamic_cast<TH2*>(h0Int0->Clone("hSum"));
    hSum->Add(h0Int1);
    const Double_t qx = 0.95; // 95% quantile
    TH1* h1 = hSum->ProjectionX();
    h1->GetQuantiles(1, chargeQuantiles+ch, &qx);
    delete h1;
    delete hSum;
  }
  Double_t meanQuantiles[2] = { 0,0 };
  Double_t nGood[2]         = { 0,0 };
  for (Int_t ch=0; ch<16; ++ch) {
    if (calibData->IsChannelDead(ch))
      continue;
    nGood[ch/8]         += 1 ;
    meanQuantiles[ch/8] += chargeQuantiles[ch];
    AliDebug(3, Form("chargeQuantiles[%02d] = %f", ch, chargeQuantiles[ch]));
  };
  if (nGood[0]) meanQuantiles[0] /= nGood[0];
  if (nGood[1]) meanQuantiles[1] /= nGood[1];

  // (2) compute charge equalization factors
  Float_t chargeEqualizationFactor = 1.0f;
  t->Branch("chargeEqualizationFactor", &chargeEqualizationFactor);

  // (3) fill TTree
  for (Int_t ch=0; ch<16; ++ch) { // offline channel number
    f_Int0.Clear("C");
    f_Int1.Clear("C");

    chOffline = ch;
    chOnline  = gOffline2Online[ch];

    chargeEqualizationFactor = (chargeQuantiles[ch] && meanQuantiles[ch/8] ? meanQuantiles[ch/8]/chargeQuantiles[ch] : 1.0 );
    AliInfo(Form("quantile[%2d] = %f f=%f", ch, chargeQuantiles[ch], chargeEqualizationFactor));

    const Double_t largeThr = 1e5;
    for (Int_t bc=0; bc<21; ++bc) {
      switch (bc) {
      case 9: {
	new (f_Int0[bc]) TF1(GetFcnName(ch, bc, 0), "[0]*TMath::TanH([1]*(x-[2]))+[3]", 0, 600);
	new (f_Int1[bc]) TF1(GetFcnName(ch, bc, 1), "[0]*TMath::TanH([1]*(x-[2]))+[3]", 0, 600);
	doExtrapolation[bc] = kTRUE;
	break;
      }
      case 10: {
	new (f_Int0[bc]) TF1(GetFcnName(ch, bc, 0), "[0] + [1]*x");
	new (f_Int1[bc]) TF1(GetFcnName(ch, bc, 1), "[0] + [1]*x");
	doExtrapolation[bc] = kTRUE;
	break;
      }
      case 11:
      case 12:
      case 13:
      case 14:
      case 15: {
	new (f_Int0[bc]) TF1(GetFcnName(ch, bc, 0), "[0] + [1]*x + [2]*abs(x)**[3]", 0, 600);
	new (f_Int1[bc]) TF1(GetFcnName(ch, bc, 1), "[0] + [1]*x + [2]*abs(x)**[3]", 0, 600);
	doExtrapolation[bc] = kTRUE;
	break;
      }
      default:
	doExtrapolation[bc]         = kFALSE;
	extrapolationThresholds[bc] = -999.9f;
	new (f_Int0[bc]) TF1(GetFcnName(ch, bc, 0), "x", 0, 600);
	new (f_Int1[bc]) TF1(GetFcnName(ch, bc, 1), "x", 0, 600);
      }

      if (doExtrapolation[bc]) {
	TH2* h0 = dynamic_cast<TH2*>(fList->FindObject(GetHistName(ch, bc, 0))); // integrator0
	TH2* h1 = dynamic_cast<TH2*>(fList->FindObject(GetHistName(ch, bc, 1))); // integrator1
	Bool_t fitOk[2] = { kTRUE,    kTRUE    };
	Double_t thr[2] = { largeThr, largeThr };
	fitOk[0] &= MakeExtrapolationFit(h0, static_cast<TF1*>(f_Int0[bc]), ch, bc, thr[0]);
	fitOk[1] &= MakeExtrapolationFit(h1, static_cast<TF1*>(f_Int1[bc]), ch, bc, thr[1]);
	doExtrapolation[bc] = (fitOk[0] || fitOk[1]);
	extrapolationThresholds[bc] = ((fitOk[0] && fitOk[1])
				       ? TMath::Min(thr[0], thr[1])
				       : fitOk[0]*thr[0] + fitOk[1]*thr[1]);
	if (!doExtrapolation[bc])
	  extrapolationThresholds[bc] = -999.9f;
      }
    }
    t->Fill();
  } // next channel

  t->ResetBranchAddresses();

  return t;
}

AliCDBEntry* AliAnalysisTaskADCalib::UpdateGainParameters(Int_t runNumber, TTree* tSat, AliADCalibData* calibData)
{
  if (NULL == tSat) {
    AliError("NULL == tSat");
    return NULL;
  }

  AliCDBManager *man = AliCDBManager::Instance();
  if (NULL == man) {
    fStatus = kInputError;
    AliFatal("CDB manager not found");
    return NULL;
  }

  if (!man->IsDefaultStorageSet())
    man->SetDefaultStorage("raw://");

  man->SetRun(runNumber);

  // (1b) Get the AD PM gain calibration object (to be updated)
  AliCDBEntry *entry = man->Get("AD/Calib/PMGains");
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
  Double_t meanEq[2] = { 0,0 };
  Double_t nGood[2]  = { 0,0 };
  for (Int_t ch=0; ch<16; ++ch) {
    if (calibData->IsChannelDead(ch))
      continue;
    nGood[ch/8]  += 1;
    meanEq[ch/8] += eqFactors[ch];
    AliDebug(3, Form("eqFactor[%02d] = %f", ch, eqFactors[ch]));
  }
  if (nGood[0]) meanEq[0] /= nGood[0];
  if (nGood[1]) meanEq[1] /= nGood[1];
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
    (Float_t)TMath::Mean(8, mip),
    (Float_t)TMath::Mean(8, mip+8)
  };
  AliInfo(Form("meanMIP=%.3f %.3f", meanMIP[0], meanMIP[1]));

  // (3) update the a and b parameters in such a way to have the same MIP per side for the given HV
  for (Int_t ch=0; ch<16; ++ch) {
    if (calibData->IsChannelDead(ch))
      continue;
    a[ch] *= TMath::Power(meanEq[ch/8] * meanMIP[ch/8]/mip[ch], -1./b[ch]);
    AliInfo(Form("a=%.2f b=%.2f mip=%.3f mip'=%.3f",
		 a[ch], b[ch], mip[ch], TMath::Power(hv[ch]/a[ch], b[ch])));
    h2Gain->SetBinContent(1+ch, 1, a[ch]);
  }

  return entry;
}
