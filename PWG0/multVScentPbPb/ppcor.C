#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TList.h"
#include "TFile.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TH2F.h"
#include "THnSparse.h"
#include "TLegend.h"
#include "TSystem.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TF1.h"
#include "TLine.h"
#include "TPaletteAxis.h"
#include "TArrayD.h"
#include "TGraphErrors.h"
//
//
#endif


TH2 *hMCAccTr,*hMCAccClsS,*hMCAccCls,*hMCGen,*hzvMC2;
TH2 *hDtAccTr,*hDtAccClsS,*hDtAccCls,*hzvDt2;
TH1 *hzvMC,*hzvDt,*hMCTrue;
//
TH2 *hTrAcceptance,*hTrDtAccNrm,*hTrDtCorr2D;
TH2 *hClAcceptance,*hClDtAccNrm,*hClDtCorr2D;
TH1 *hTrDtCorr1D,*hClDtCorr1D;
//
// for test of MC on MC only
TH2 *hMCGenData = 0;
TH1 *hMCTrueData = 0;

TH1* GetHisto(TList* lst, const char* name, const char* pref="");
void CorrectMult(TH2* hMCAcc, TH2* hMCGn,
		 TH2* hDtAcc, TH1* hzDt,
		 const char* pref,
		 // output
		 TH2 *&hDtAccNrm,    // raw data normalized per Z bin
		 TH2 *&hAcceptance,  // 2D acceptance
		 TH2 *&hDtCorr2D,    // 2D corrected dN/dEta per Z bin (norm up to bin width)
		 TH1 *&hDtCorr1D     // 1D corrected dN/dEta 
		 );


const double minZStat = 50;

void ppcor(const char* flMC, const char* flDt)
{
  TFile *fileMC = TFile::Open(flMC);
  TList* clistMC = (TList*)fileMC->Get("clist");
  TFile *fileDt = TFile::Open(flDt);
  TList* clistDt = (TList*)fileDt->Get("clist");
  //
  hMCAccTr   = (TH2*)GetHisto(clistMC,"b0_TrData_ZvEtaCutT","_mc"); 
  hMCAccClsS = (TH2*)GetHisto(clistMC,"b0_TrData_ZvEtaSPD1","_mc");
  hMCAccCls  = (TH2*)hMCAccClsS->Clone( Form("MCClusSPD_%s","_mc") );
  hMCAccCls->Add(hMCAccTr);
  hMCGen   = (TH2*)GetHisto(clistMC,"b0_zvEtaPrimMC","_mc");
  hzvMC2   = (TH2*)GetHisto(clistMC,"zv","_mc");
  hzvMC    =  hzvMC2->ProjectionX("zvgenMC");
  //
  hMCTrue = hMCGen->ProjectionX("MCTrue");
  if (hzvMC2->Integral()>0) hMCTrue->Scale(1./hzvMC2->Integral());
  hMCTrue->Scale(1./hMCTrue->GetBinWidth(1));
  //
  hDtAccTr   = (TH2*)GetHisto(clistDt,"b0_TrData_ZvEtaCutT","_dt"); 
  hDtAccClsS = (TH2*)GetHisto(clistDt,"b0_TrData_ZvEtaSPD1","_dt"); 
  hDtAccCls  = (TH2*)hDtAccClsS->Clone( Form("DtClusSPD_%s","_dt") );
  hDtAccCls->Add(hDtAccTr);
  hzvDt2   = (TH2*)GetHisto(clistDt,"zv","_dt");
  hzvDt   =  hzvDt2->ProjectionX("zvDt");
  //
  hMCGenData = (TH2*)GetHisto(clistDt,"b0_zvEtaPrimMC","_data");
  if (hMCGenData) {
    hMCTrueData = hMCGenData->ProjectionX("MCTrueData");
    if (hzvDt2->Integral()>0) hMCTrueData->Scale(1./hzvDt2->Integral());
    hMCTrueData->Scale(1./hMCTrueData->GetBinWidth(1));
  }
  //
  CorrectMult(hMCAccTr,hMCGen,
	      hDtAccTr,hzvDt, "trc",
	      hTrDtAccNrm,hTrAcceptance,hTrDtCorr2D,hTrDtCorr1D);
  //
  CorrectMult(hMCAccCls,hMCGen,
	      hDtAccCls,hzvDt, "cls",
	      hClDtAccNrm,hClAcceptance,hClDtCorr2D,hClDtCorr1D);
  //
  //
  hTrDtCorr1D->Draw();  // dndeta from tracklets
  hClDtCorr1D->Draw();  // dndeta from clusters
}


TH1* GetHisto(TList* lst, const char* name, const char* pref)
{
  TH1* hst = (TH1*)lst->FindObject(name);
  if (!hst) return 0;
  TH1* hcl = (TH1*) hst->Clone( Form("%s%ss",hst->GetName(),pref) );
  return hcl;
}

void CorrectMult(TH2* hMCAcc, TH2* hMCGn,
		 TH2* hDtAcc, TH1* hzDt,
		 const char* pref,
		 // output
		 TH2 *&hDtAccNrm,    // raw data normalized per Z bin
		 TH2 *&hAcceptance,  // 2D acceptance
		 TH2 *&hDtCorr2D,    // 2D corrected dN/dEta per Z bin (norm up to bin width)
		 TH1 *&hDtCorr1D     // 1D corrected dN/dEta 
		 )
{
  //
  int neta = hDtAcc->GetXaxis()->GetNbins();
  int nz   = hDtAcc->GetYaxis()->GetNbins();  
  //
  hDtAccNrm = (TH2*) hDtAcc->Clone( Form("rawDataNormZ_%s",pref) );
  hDtAccNrm->Reset();
  for (int iz=1;iz<=nz;iz++) { // normalize data histo per Z bin
    double zv = hDtAcc->GetYaxis()->GetBinCenter(iz);
    int iz0 = hzDt->FindBin(zv);
    double zst  = hzDt->GetBinContent(iz0);
    double zstE = hzDt->GetBinError(iz0);
    if (zst<minZStat) continue;
    for (int ix=1;ix<=neta;ix++) {
      double vl = hDtAcc->GetBinContent(ix,iz);
      double er = hDtAcc->GetBinError(ix,iz);
      if (vl<1e-9 || er<1e-6) continue;
      double rat  = vl/zst;
      double ratE = rat*TMath::Sqrt( er*er/(vl*vl) + (zstE*zstE)/(zst*zst) );
      //printf("%2d %2d  %+e(%+e) -> %+e(%+e)\n",ix,iz,vl,er,rat,ratE);
      hDtAccNrm->SetBinContent(ix,iz,rat);
      hDtAccNrm->SetBinError(ix,iz,ratE);
    }
  }
  //
  hAcceptance = (TH2*)hMCAcc->Clone( Form("Acceptance_%s",pref) );
  hAcceptance->Divide(hMCGn);
  //
  hDtCorr2D = (TH2*)hDtAccNrm->Clone( Form("dNdEtaCor2D_%s",pref) );
  hDtCorr2D->Divide(hAcceptance);
  //
  hDtCorr1D = hDtCorr2D->ProjectionX( Form("dNdEtaCor1D_%s",pref) );
  hDtCorr1D->Reset();
  //
  for (int ix=1;ix<=neta;ix++) { // get weighted average
    double sm=0,sme=0;
    double bw = hDtCorr1D->GetBinWidth(ix);
    for (int iz=1;iz<=nz;iz++) {
      double vl = hDtCorr2D->GetBinContent(ix,iz);
      double er = hDtCorr2D->GetBinError(ix,iz);
      if (vl<1e-6 || er<1e-12) continue;
      sm += vl/(er*er);
      sme += 1./(er*er);
    }
    if (sme<1e-6) continue;
    sm /= sme;
    sme = 1./TMath::Sqrt(sme);
    hDtCorr1D->SetBinContent(ix,sm/bw);
    hDtCorr1D->SetBinError(ix,sme/bw);
  }
  //
}
