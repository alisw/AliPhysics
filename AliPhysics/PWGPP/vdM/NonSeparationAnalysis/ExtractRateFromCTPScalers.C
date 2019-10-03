// -*- C++ -*-
// $Id$

#include <TFile.h>
#include <TTree.h>
#include <TGraphErrors.h>
#include <TGraph.h>
#include <TMath.h>
#include <TF1.h>

Double_t ComputeMeanRate(Long64_t nSel, Double_t *dt, Double_t *rs, Double_t &rateErr, Double_t &counts) {
  rateErr = 0;
  counts = 0;
  Double_t deltaT = 0;
  for (Long64_t i=0; i<nSel-1; ++i) {
    counts += rs[i]*dt[i];
    deltaT += dt[i];
    Printf("ComputeMeanRate %5lld %f %f %f", i, rs[i], rs[i]*dt[i], counts);
  }
  Printf("ComputeMeanRate %f %f", counts, deltaT);
  rateErr = TMath::Sqrt(counts)/deltaT;
  return counts/deltaT;
}


struct RateCorr {
  RateCorr(Double_t rA=0.116,
	   Double_t rC=0.266)
    : f("RateCorr", RateCorr::rAD, 0.0, 4.0, 2)
    , ratioA(rA)
    , ratioC(rC) {
    f.SetParameters(ratioA, ratioC);
  }

  Double_t operator()(Double_t y) const {
    Double_t xmin(f.GetXmin());
    Double_t xmax(f.GetXmax());
    Double_t xmed(0), val(0);
    for(Int_t i=0; i<100000; i++){
      xmed=(xmax+xmin)/2.;
      val=f.Eval(xmed);
      if(val>=y)
	xmax=xmed;
      else
	xmin=xmed;
    }
    xmed=(xmax+xmin)/2.;
    return xmed;
  }

  // this function defines the relationship between ratePerBC and mu
  static Double_t rAD(Double_t * x, Double_t * par) {
    const Double_t eMinMu  = TMath::Exp(-x[0]);
    const Double_t eMinMuA = TMath::Exp(-x[0]*par[0]);
    const Double_t eMinMuC = TMath::Exp(-x[0]*par[1]);
    return (1-eMinMu + eMinMu * (1-eMinMuA) * (1-eMinMuC));
  }

  TF1 f;
  const Double_t ratioA;
  const Double_t ratioC;
} ;

struct THeader {
  THeader()
    : tBeg(0)
    , tEnd(0)
    , aqFlag(kTRUE)
    , xSep_mm(0)
    , ySep_mm(0) {}

  Double_t tBeg;    // begin of step (unix seconds)
  Double_t tEnd;    // end of step (unix seconds)
  Bool_t   aqFlag;  // AQ flag (always 1 for this dataset)
  Double_t xSep_mm; // nominal beam separation [mm]
  Double_t ySep_mm; // nominal beam separation [mm]
};

struct TData {
  TData(Int_t nB=1)
    : countsL0b(0)
    , nBunches(nB)
    , rateRaw_Hz(0)
    , rateRawErr_Hz(0)
    , rateBB_Hz(0)
    , rateBBErr_Hz(0)
    , ratePU_Hz(0)
    , ratePUErr_Hz(0)
    , rateDC_Hz(0)
    , rateDCErr_Hz(0)  {}
  Double_t countsL0b;      // L0b counts in the interval
  Int_t    nBunches;       // number of bunches for the trigger class (-B: 27, -I1..8: 1)
  Double_t rateRaw_Hz;     // raw rate 
  Double_t rateRawErr_Hz;  // rate error
  Double_t rateBB_Hz;      // rate with bkgd. subtraction
  Double_t rateBBErr_Hz;   // rate error 
  Double_t ratePU_Hz;      // rate with pile-up correction + bkgd. subtraction
  Double_t ratePUErr_Hz;   // rate error
  Double_t rateDC_Hz;      // rate with intensity decay correction + pile-up correction + bkgd. subtraction 
  Double_t rateDCErr_Hz;   // rate error
};

Double_t FindInGraph(const TGraphErrors *g, Double_t x, Double_t eps, Double_t &errY) {
  const Double_t *gx  = g->GetX();
  const Double_t *gy  = g->GetY();
  const Double_t *gey = g->GetEY();
  for (Int_t i=0, n=g->GetN(); i<n; ++i) {
    if (TMath::Abs(gx[i]-x) < eps) {
      errY = gey[i];
      return gy[i];
    }
  }
  errY = -1;
  return -1;
}

TGraphErrors* DrawRateVsSep(TTree *tRate, TGraph *gSep,
			    TGraph* gCounter,
			    Double_t tBeg, Double_t tEnd,
			    Int_t nB,
			    TGraphErrors *gRateSepNoCorr,
			    Double_t pileupA, Double_t pileupC,
			    TTree *TT,
			    TString dataName,
			    TString scanType, Double_t offset_mm,
			    TGraphErrors *gSatCont,
			    Int_t satContVariation,
			    TF1* fIntDecay) {

  const Double_t *t = gSep->GetX();
  const Double_t *s = gSep->GetY();
  const Double_t *c = gCounter->GetY();
  Long64_t counter = -1;

  const RateCorr fRateCorr(pileupA, pileupC);

  const Bool_t isADAND(kFALSE);//title.Contains("ADAND") || kTRUE);

  TBranch *bHeader = NULL;
  TBranch *bData   = NULL;
  THeader  *header = new THeader;;
  TData    *data   = NULL;
  if (TT) {
    if (!TT->GetBranch("Header")) {
      bHeader = TT->Branch("Header", &header);
    } else {
      bHeader = TT->GetBranch("Header");
      TT->SetBranchAddress("Header", &header);
    }
    if (!TT->GetBranch(dataName)) {
      bData = TT->Branch(dataName, &data);
    } else {
      bData = TT->GetBranch(dataName);
      TT->SetBranchAddress(dataName, &data);
    }
  }

  TGraphErrors *gRateSep = new TGraphErrors;
  Long64_t nSel = 0;
  for (Long64_t i=0, n=gSep->GetN(); i<n; ++i)  {
    if (t[i] >= tBeg && t[i] <= tEnd) {
      header->xSep_mm = 0;
      header->ySep_mm = 0;
      if (c[i] == counter) {
	header->tEnd = t[i];
	Printf("ZZ %.0f %.0f %lld", header->tBeg, header->tEnd, counter);
	// tRate->Scan("rate:rateErr**-2:time", TString::Format("time>=%.0f && time<=%.0f", header->tBeg, header->tEnd));
	nSel = tRate->Draw("rate:rateErr**-2:dt", TString::Format("time>=%.0f && time<=%.0f", header->tBeg+0*3600, header->tEnd+0*3600), "GOFF");
	*data = TData(nB); // reset
	Printf("nSel=%lld", nSel);
	if (nSel > 1) {
	  const Double_t intDecayCorr = (fIntDecay
					 ? fIntDecay->Eval(0.5*(header->tBeg + header->tEnd))
					 : 1.0);
	  Double_t satContErr = 0;
	  Double_t satCont = (gSatCont
			      ? FindInGraph(gSatCont, s[i], 1e-4, satContErr)
			      : 0.0);
	  Printf("SATCONT: %f %e", satCont, satContErr);
	  Double_t satCorr = TMath::Min(1.0, TMath::Max(0.0, 1 - satCont - satContVariation*(satContErr)));

	  Printf("nSel = %lld sep=%f intDecayCorr=%f satCorr=%f", nSel, s[i], intDecayCorr, satCorr);

	  Double_t rateErr = 0;
	  const Double_t rate = ComputeMeanRate(nSel, tRate->GetV3(), tRate->GetV1(), rateErr, data->countsL0b);
	  const Double_t rateCorrOld = -11245.*nB*TMath::Log(1-rate/nB/11245.);
	  const Double_t rateCorr    = nB*11245.*fRateCorr(rate/nB/11245.0);
	  const Double_t rateCorrErr = TMath::Max(TMath::Abs(nB*11245.*fRateCorr(((rate+rateErr)/nB/11245.0)) - rateCorr),
						  TMath::Abs(nB*11245.*fRateCorr(((rate-rateErr)/nB/11245.0)) - rateCorr));

	  data->rateRaw_Hz     = rate;
	  data->rateRawErr_Hz  = rateErr;

	  data->rateBB_Hz      = satCorr * rate;
	  data->rateBBErr_Hz   = satCorr * rateErr;

	  data->ratePU_Hz      = nB*11245.*fRateCorr(data->rateBB_Hz/nB/11245.0);
	  data->ratePUErr_Hz   = TMath::Max(TMath::Abs(nB*11245.*fRateCorr(((data->rateBB_Hz + data->rateBBErr_Hz)/nB/11245.0)) - data->ratePU_Hz),
					    TMath::Abs(nB*11245.*fRateCorr(((data->rateBB_Hz - data->rateBBErr_Hz)/nB/11245.0)) - data->ratePU_Hz));

	  data->rateDC_Hz      = intDecayCorr * data->ratePU_Hz;
	  data->rateDCErr_Hz   = intDecayCorr * data->ratePUErr_Hz;

	  Printf("RateCorr: %f %f | %f %f", rate, rateCorrOld, rateCorr, rateCorrErr);
	  gRateSep->SetPoint     (gRateSep->GetN(), s[i], (isADAND ? rateCorr : rateCorrOld));
	  gRateSep->SetPointError(gRateSep->GetN()-1,
				  0,
				  (isADAND ? rateCorrErr: rateCorrOld/(1e-10+rate)*rateErr));
	  gRateSepNoCorr->SetPoint(gRateSepNoCorr->GetN(), s[i], rate);

	  header->xSep_mm = scanType.Contains("X") ? s[i] : offset_mm;
	  header->ySep_mm = scanType.Contains("Y") ? s[i] : offset_mm;
	}
	if (TT->GetNbranches() == 2)
	  TT->Fill();
	else
	  bData->Fill();
      }
      header->tBeg = t[i];
      counter     = c[i];
    }
  }

  TT->ResetBranchAddresses();

  gRateSep->Fit("gaus");
  TF1* fg = (TF1*)gRateSep->FindObject("gaus");

  TF1* fGaus = new TF1("gaus", "gaus(0)*(1+([3]-[1])*x**2+([4]-[1])*x**4+([5]-[1])*x**6)", -0.4, 0.4);
  fGaus->SetParNames("R^{C} [Hz]", "#mu [mm]", "#sigma_{scan} [mm]", "p_{2} [mm^{-2}]", "p_{4} [mm^{-4}]", "p_{6} [mm^{-6}]", "C [Hz]");
  fGaus->SetParameter(0, fg->GetParameter(0));
  fGaus->SetParameter(1, fg->GetParameter(1));
  fGaus->SetParameter(2, fg->GetParameter(2));
  fGaus->SetParLimits(3, 0, 1e4);
  fGaus->SetParLimits(4, 0, 1e4);
  fGaus->SetParLimits(5, 0, 1e4);
//   fGaus->SetParLimits(6, 0, 100);
  fGaus->SetNpx(10000);

  gRateSep->Fit(fGaus, "Q");
  fGaus->SetNpx(10000);

  return gRateSep;
}

TTree* MakeRateTTree(Long64_t n, Double_t *timeSeconds, Double_t *counter) {
  struct {
    Double_t time;
    Double_t timeErr;
    Double_t rate;
    Double_t rateErr;
    Double_t dt;
  } rateData;
  TTree *tRate = new TTree;
  tRate->Branch("rate", &rateData, "time/D:timeErr:rate:rateErr:dt");

  Double_t dt = 0;
  for (Long64_t i=0; i<n-1; ++i) {
    rateData.dt      = timeSeconds[i+1] - timeSeconds[i];
    rateData.time    = 0.5*(timeSeconds[i]+timeSeconds[i+1]);
    rateData.timeErr = rateData.dt/2;
    rateData.rate    = counter[i+1] - counter[i];
    rateData.rateErr = sqrt(TMath::Max(1., rateData.rate))/rateData.dt;
    rateData.rate   /= rateData.dt;
    Printf("YY %f %f %f %f", rateData.time, rateData.dt, rateData.rate, rateData.rateErr);
    tRate->Fill();
  }
  tRate->ResetBranchAddresses();
  return tRate;
}


void ExtractRateFromCTPScalers(TString classID, TString bunchID, Int_t rn, TString scanInfoFileName, Int_t nBunches, Int_t fillNumber,
			       Double_t pileupA, Double_t pileupC, TString intDecayFileName, TString satContFileName, Int_t satContVariation,
			       TTree *TT) {

  TTree   tScanInfo;
  tScanInfo.ReadFile(scanInfoFileName, "scanType/C:scanBeg/D:scanEnd:offset_mm:yMaxLin:yMaxLog");
  Long64_t nScans   = tScanInfo.GetEntries();

  char     scanType[1024] = { 0 };
  Double_t scanBeg   = 0;
  Double_t scanEnd   = 0;
  Double_t offset_mm = 0;
  Double_t yMaxLin   = 0;
  Double_t yMaxLog   = 0;
  tScanInfo.SetBranchAddress("scanType",  scanType);
  tScanInfo.SetBranchAddress("scanBeg",   &scanBeg);
  tScanInfo.SetBranchAddress("scanEnd",   &scanEnd);
  tScanInfo.SetBranchAddress("offset_mm", &offset_mm);
  tScanInfo.SetBranchAddress("yMaxLin",   &yMaxLin);
  tScanInfo.SetBranchAddress("yMaxLog",   &yMaxLog);

  TFile *fS = TFile::Open(Form("root/scalers_%d.root", rn));
  TTree *TS = (TTree*)gFile->Get("TS");
  const Long64_t n = TS->Draw(TString::Format("fSec+1e-9*fNanoSec:%s_%s.L0b", classID.Data(), bunchID.Data()), "", "GOFF");
  TTree *tRate = MakeRateTTree(n, TS->GetV1(), TS->GetV2());

  tRate->Draw("time:rate:timeErr:rateErr", "", "GOFF");
  TGraphErrors *gRate       = new TGraphErrors(n-1, tRate->GetV1(), tRate->GetV2(), tRate->GetV3(), tRate->GetV4());

  TTree tSep;
  tSep.ReadFile(Form("txt/%d/sep.txt", fillNumber), "ts/D:sep:counter/I:ts1/C:ts2/C");
  Long64_t nSep = tSep.Draw("ts:sep:counter", "", "GOFF");
  TGraph *gSep     = new TGraph(nSep, tSep.GetV1(), tSep.GetV2());
  TGraph *gCounter = new TGraph(nSep, tSep.GetV1(), tSep.GetV3());

  TGraphErrors **gRateVsSep       = new TGraphErrors*[nScans];
  TGraphErrors **gRateVsSepNoCorr = new TGraphErrors*[nScans];
  for (Long64_t i=0; i<nScans; ++i) {
    tScanInfo.GetEntry(i);
    TFile *fID = TFile::Open(intDecayFileName);
    TF1 * fIntDecay = (fID
		       ? (TF1*)fID->Get(Form("fIntCorr_%s_Scan%lld", bunchID.Data(), 1+(i/2)))
		       : NULL);

    TGraphErrors *gSatCont = NULL;
    TFile *fSC = TFile::Open(satContFileName);
    gSatCont = (TGraphErrors*)fSC->Get(Form("gSatScan_%s", scanType));
    Printf("=== %s %.0f %.0f %p %p", scanType, scanBeg, scanEnd, gSatCont, fIntDecay);
    gRateVsSepNoCorr[i] = new TGraphErrors;
    gRateVsSep[i] = DrawRateVsSep(tRate, gSep, gCounter, scanBeg, scanEnd,  nBunches, gRateVsSepNoCorr[i], pileupA, pileupC,
				  TT, classID+"_"+bunchID,
				  TString(scanType), offset_mm, gSatCont, satContVariation, fIntDecay);
    fSC->Close();
    if (fID)
      fID->Close();
  }

  tScanInfo.ResetBranchAddresses();
  tSep.ResetBranchAddresses();

  TFile *fSave = TFile::Open(Form("root/%d/%d_%s-%s.root", fillNumber, rn, classID.Data(), bunchID.Data()), "RECREATE");
  gRate->Write("gRate", TObject::kWriteDelete);
  gSep->Write("gSep", TObject::kWriteDelete);
  gCounter->Write("gCounter", TObject::kWriteDelete);
  for (Long64_t i=0; i<nScans; ++i) {
    gRateVsSep[i]->Write(Form("gRateVsSep_Scan%lld", i));
    gRateVsSepNoCorr[i]->Write(Form("gRateVsSepNoPileupCorr_Scan%lld", i));
  }
  fSave->Write();
  fSave->Close();
  fS->Close();
}
