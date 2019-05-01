#include "LInfo.h"
#include <TCanvas.h>
#include <TDatime.h>
#include <TFile.h>
#include <TGrid.h>
#include <TH2F.h>
#include <TMap.h>
#include <TMath.h>
#include <TNtuple.h>
#include <TProfile2D.h>

void LInfo::Compute()
{
  if (fIsComputed)
    return;

  char id[100];
  char title[100];
  const char *sideStr[] = {"A", "C"};
  const Int_t kNCol     = NCol();
  const Int_t kNRow     = NRow();
  const Int_t kNStrip   = NStrip();

  for (Int_t iSM=0; iSM<kNSM; ++iSM) {
    Int_t isector = iSM/2;
    Int_t iside = iSM%2;
    for (Int_t igain=0; igain<2; ++igain) {
      if (!fhAmpOverMon[iSM][igain]) {
        sprintf(id, "hAmpOverMon%02d%d", iSM, igain);
        sprintf(title, "LED amplitude over LEDMON: SM %d (%1s%d) gain %d", iSM, sideStr[iside], isector, igain);
        fhAmpOverMon[iSM][igain] = new TH2F(id, title, kNCol, -0.5, kNCol-0.5, kNRow, -0.5, kNRow - 0.5);
        fhAmpOverMon[iSM][igain]->SetDirectory(0);
        sprintf(id, "hStripRmsOverMean%02d%d", iSM, igain);
        sprintf(title, "LedMon rms  over mean: SM %d (%1s%d) gain %d", iSM, sideStr[iside], isector, igain);
        fhStripRmsOverMean[iSM][igain] = new TH1F(id, title, kNCol, -0.5, kNCol-0.5);
        fhStripRmsOverMean[iSM][igain]->SetDirectory(0);
        sprintf(id, "hLedRmsOverMean%02d%d", iSM, igain);
        sprintf(title, "Led rms  over mean: SM %d (%1s%d) gain %d", iSM, sideStr[iside], isector, igain);
        fhLedRmsOverMean[iSM][igain] = new TH2F(id, title, kNCol, -0.5, kNCol-0.5, kNRow, -0.5, kNRow-0.5);
        fhLedRmsOverMean[iSM][igain]->SetDirectory(0);
      } else {
        fhAmpOverMon[iSM][igain]->Reset();
        fhStripRmsOverMean[iSM][igain]->Reset();
        fhLedRmsOverMean[iSM][igain]->Reset();
      }
    }
  }

  for (Int_t iSM=0; iSM<kNSM; ++iSM) {
    Int_t isector = iSM/2;
    Int_t iside = iSM%2;
    for (Int_t igain=0; igain<2; ++igain) {
      for (Int_t col=0;col<kNCol;++col) {
        Int_t strip = col / 2;
        Int_t mbin = fhStrip[iSM][igain]->FindBin(strip);
        Double_t ledMonAmp = fhStrip[iSM][igain]->GetBinContent(mbin);
        for (Int_t row=0;row<kNRow;++row) {
          Int_t lbin = fhLed[iSM][igain]->FindBin(col,row);
          Double_t ledAmp = fhLed[iSM][igain]->GetBinContent(lbin);
          if (ledMonAmp!=0) {
            Double_t weightf = ledAmp/ledMonAmp;
            fhAmpOverMon[iSM][igain]->SetBinContent(lbin,weightf);
          }
          Double_t ledRms = fhLedCount[iSM][igain]->GetBinContent(lbin);
          if (ledAmp>0) {
            fhLedRmsOverMean[iSM][igain]->SetBinContent(lbin,ledRms/ledAmp);
          }
        }
      }
      for (Int_t strip=0;strip<kNStrip;++strip) {
        Int_t bin = fhStrip[iSM][igain]->FindBin(strip);
        Double_t mean = fhStrip[iSM][igain]->GetBinContent(bin);
        Double_t rms = fhStripCount[iSM][igain]->GetBinContent(bin);
        if (mean>0)
          fhStripRmsOverMean[iSM][igain]->SetBinContent(bin,rms/mean);
      }
    }
  }
  fIsComputed = 1;
}

void LInfo::CreateHistograms()
{
  // book histograms
  char id[100];
  char title[100];
  const char *sideStr[] = {"A", "C"};
  const Int_t kNCol     = NCol();
  const Int_t kNRow     = NRow();
  const Int_t kNStrip   = NStrip();
  const char *opt = "S";

  for (Int_t iSM=0; iSM<kNSM; ++iSM) {
    Int_t isector = iSM/2;
    Int_t iside = iSM%2;
    for (Int_t igain=0; igain<2; ++igain) {
      sprintf(id, "hStrip%02d%d", iSM, igain);
      sprintf(title, "LEDMon Amplitude: SM %d (%1s%d) gain %d", iSM, sideStr[iside], isector, igain);
      fhStrip[iSM][igain] = new TProfile(id, title, kNStrip, -0.5, kNStrip-0.5, opt);
      fhStrip[iSM][igain]->SetDirectory(0);

      sprintf(id, "hStripCount%02d%d", iSM, igain);
      sprintf(title, "LEDMon Entries: SM %d (%1s%d) gain %d", iSM, sideStr[iside], isector, igain);
      fhStripCount[iSM][igain] = new TProfile(id, title, kNStrip, -0.5, kNStrip-0.5, opt);
      fhStripCount[iSM][igain]->SetDirectory(0);

      sprintf(id, "hStripWeighted%02d%d", iSM, igain);
      sprintf(title, "LEDMon Weighted Amplitude: SM %d (%1s%d) gain %d", iSM, sideStr[iside], isector, igain);
      fhStripWeighted[iSM][igain] = new TProfile(id, title, kNStrip, -0.5, kNStrip-0.5, opt);
      fhStripWeighted[iSM][igain]->SetDirectory(0);

      sprintf(id, "hLed%02d%d", iSM, igain);
      sprintf(title, "Led Amplitude: SM %d (%1s%d) gain %d", iSM, sideStr[iside], isector, igain);
      fhLed[iSM][igain] = new TProfile2D(id, title, kNCol, -0.5, kNCol-0.5, kNRow, -0.5, kNRow - 0.5, opt);
      fhLed[iSM][igain]->SetDirectory(0);

      sprintf(id, "hLedCount%02d%d", iSM, igain);
      sprintf(title, "Led Entries: SM %d (%1s%d) gain %d", iSM, sideStr[iside], isector, igain);
      fhLedCount[iSM][igain] = new TProfile2D(id, title, kNCol, -0.5, kNCol-0.5, kNRow, -0.5, kNRow - 0.5, opt);
      fhLedCount[iSM][igain]->SetDirectory(0);

      sprintf(id, "hLedWeighted%02d%d", iSM, igain);
      sprintf(title, "Led Weighted Amplitude: SM %d (%1s%d) gain %d", iSM, sideStr[iside], isector, igain);
      fhLedWeighted[iSM][igain] = new TProfile2D(id, title, kNCol, -0.5, kNCol-0.5, kNRow, -0.5, kNRow - 0.5, opt);
      fhLedWeighted[iSM][igain]->SetDirectory(0);

      fhAmpOverMon[iSM][igain]       = 0;
      fhLedRmsOverMean[iSM][igain]   = 0;
      fhStripRmsOverMean[iSM][igain] = 0;
    }
  }
}

void LInfo::FillLed(Int_t mod, Int_t gain, Int_t col, Int_t row, Double_t amp, Double_t rms)
{
  if ((amp<0)||(amp>1500))
    return;
  fhLed[mod][gain]->Fill(col, row, amp);
  if (rms>0) {
    fhLedCount[mod][gain]->Fill(col, row, rms);
    fhLedWeighted[mod][gain]->Fill(col, row, amp, 1./TMath::Power(rms,2));
  }
}

TCanvas *LInfo::DrawHist(Int_t which, Int_t gain, const char *opt) const
{
  TString name;
  TString hopt;
  if ((which<1||which>4) && !fIsComputed) {
    cerr << "Execute Linfo::Compute first" << endl;
    return 0;
  }
  TCanvas *c = new TCanvas("c","c",1600,1600);
  c->Divide(2,10);
  c->SetLeftMargin(0.02);
  c->SetRightMargin(0.02);
  c->SetTopMargin(0.02);
  c->SetBottomMargin(0.02);
  for (Int_t sm=0;sm<kNSM;++sm) {
    c->cd(sm+1);
    gPad->SetLeftMargin(0.02);
    gPad->SetRightMargin(0.02);
    gPad->SetTopMargin(0.02);
    gPad->SetBottomMargin(0.02);
    TH1 *h = 0;
    switch (which) {
    case 1:
      h = GetStripHist(sm,gain);
      name = Form("StripMeanHist_%d_%d",gain,fRunNo);
      break;
    case 2:
      h = GetStripRmsHist(sm,gain);
      name = Form("StripRmsHist_%d_%d",gain,fRunNo);
      break;
    case 3:
      h = GetLedHist(sm,gain);
      name = Form("LedMeanHist_%d_%d",gain,fRunNo);
      hopt = "colz";
      break;
    case 4:
      h = GetLedRmsHist(sm,gain);
      name = Form("LedRmsHist_%d_%d",gain,fRunNo);
      hopt = "colz";
      break;
    case 5:
      h = GetLedMonDispHist(sm,gain);
      name = Form("LedMonDispHist_%d_%d",gain,fRunNo);
      break;
    case 6:
      h = GetLedDispHist(sm,gain);
      name = Form("LedDispHist_%d_%d",gain,fRunNo);
      hopt = "colz";
      break;
    default :
      h = GetLedOverMonHist(sm,gain);
      name = Form("LedOverMonHist_%d_%d",gain,fRunNo);
      hopt = "colz";
    }
    if (opt)
      h->Draw(opt);
    else
      h->Draw(hopt);
  }
  c->SetName(name);
  c->SetTitle(name);
  return c;
}


void LInfo::FillStrip(Int_t mod, Int_t gain, Int_t strip, Double_t amp, Double_t rms)
{
  if ((amp<0)||(amp>1500))
    return;
  Int_t spos = strip;
  if (mod%2==1)
    spos = 23-strip;
  fhStrip[mod][gain]->Fill(spos, amp);
  if (rms>0) {
    fhStripCount[mod][gain]->Fill(spos, rms);
    fhStripWeighted[mod][gain]->Fill(spos, amp, 1./TMath::Power(rms,2));
  }
}

Double_t LInfo::FracLeds(Int_t sm, Int_t gain) const
{
  const Int_t kGain=gain;
  Double_t ret=0,all=0;
  for (Int_t iSM=0; iSM<kNSM; ++iSM) {
    if (sm>=0&&iSM!=sm)
      continue;
    Int_t ncols=NCol();
    if (iSM>11 && iSM<18)
      ncols=32;
    Int_t nrows=NRow();
    if (iSM==10||iSM==11||iSM==18||iSM==19)
      nrows=8;
    for (Int_t col=0;col<ncols;++col) {
      for (Int_t row=0;row<nrows;++row) {
        ++all;
        Int_t bin=fhLed[iSM][kGain]->GetBin(col,row);
        if (fhLed[iSM][kGain]->GetBinContent(bin)>0.)
        ++ret;
      }
    }
  }
  if (all>0)
    return ret/all;
  return 0;
}

Double_t LInfo::FracStrips(Int_t sm, Int_t gain) const
{
  const Int_t kstripGain=gain;
  Double_t ret=0,all=0;
  for (Int_t iSM=0; iSM<kNSM; ++iSM) {
    if (sm>=0&&iSM!=sm)
      continue;
    Int_t nstrips=NStrip();
    if (iSM>11 && iSM<18)
      nstrips=32/2;
    for (Int_t strip=1; strip<=nstrips; ++strip) {
      ++all;
      if (fhStrip[iSM][kstripGain]->GetBinContent(strip)>0.)
        ++ret;
    }
  }
  if (all>0)
    return ret/all;
  return 0;
}

void LInfo::Print(Option_t *option) const
{
  cout << "Runno: " << fRunNo << endl;
}
