/// \file drawMeanPtRaw.C
/// \brief Macro to draw the Output of AliAnalysisTaskMeanPtRaw
///
/// This macro draws the output of AliAnalysisTaskMeanPtRaw. The first argument is the file, which should be used to plot, 
/// the second argument contains additional information, which are put on the plot.
/// This second argument can contain pT cuts in the form lowerPt<pT<higherPt
/// All different information need to be seperated by ; and are then parsed in this macro.
///
/// \author Philipp Luettig <philipp.luettig@cern.ch>, University of Frankfurt, Germany
/// \date May 12th, 2015

#include "TFile.h"
#include "TList.h"
#include "TProfile.h"
#include "TH1.h"
#include "TH2.h"
#include "THnSparse.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TLegend.h"
#include <vector>

Int_t GetMaxBin(std::vector<TProfile> &vec);
Double_t GetMaxPt(std::vector<TProfile> &vec);

void drawMeanPtRaw(const char *fInput = "AnalysisResults.root", const char * cInformation ="")
{
  TFile *f = TFile::Open(fInput,"READ");
  TDirectoryFile *d = (TDirectoryFile*)f->Get("MeanPtRaw");
  TList *l = (TList*)d->Get("MeanPtRaw");
  //l->ls();
  THnSparse *hsMpt = (THnSparse*)l->FindObject("fPtVsMultRaw");
  TH1I *hCutNames = (TH1I*)l->FindObject("fTrackCutName");
  
  Double_t pTmin = 0.15;
  Double_t pTmax = 10;
  
  const Int_t colors[] = {kBlue+1, kRed+1, kBlack, kGreen+3, kMagenta+1, kOrange-1,kCyan+2,kYellow+2};
  
  // parse information given as second argument
  TString sAdditionalInfo(cInformation);
  TObjArray *oAddInfo = sAdditionalInfo.Tokenize(";");
  
  std::vector<Bool_t> bDrawEntry;
  std::vector<Double_t> dVecPtMin;
  std::vector<Double_t> dVecPtMax;
  Int_t nPtCuts = 0;
  
  for(Int_t iEntry = 0; iEntry < oAddInfo->GetEntries(); iEntry++)
  {
	bDrawEntry.push_back(kTRUE);
  }
  
  for(Int_t iEntry = 0; iEntry < oAddInfo->GetEntries(); iEntry++)
  {
	TObjString *osAdd = (TObjString*)oAddInfo->At(iEntry);
	if(osAdd->GetString().Contains("pT"))
	{
	  TObjArray *oPt = osAdd->GetString().Tokenize("<");
	  TObjString *osPtMin = (TObjString*)oPt->At(0);
	  TObjString *osPtMax = (TObjString*)oPt->At(2);
	  
	  dVecPtMin.push_back(atof(osPtMin->GetString().Data()));
	  dVecPtMax.push_back(atof(osPtMax->GetString().Data()));
	  
	  bDrawEntry.at(iEntry) = kFALSE;
	  nPtCuts++;
	}	
  }
  
  std::vector<std::vector<TProfile> > pVecMpt;
  
  Int_t nCutSettings = hsMpt->GetAxis(2)->GetNbins();
  
  for(Int_t iCutSetting = 0; iCutSetting < nCutSettings; iCutSetting++) {
	hsMpt->GetAxis(2)->SetRange(iCutSetting, iCutSetting);
	std::vector<TProfile> pVecTemp;
	for(Int_t iPtCut = 0; iPtCut < nPtCuts; iPtCut++) {
	  hsMpt->GetAxis(1)->SetRangeUser(dVecPtMin.at(iPtCut), dVecPtMax.at(iPtCut));
	  TProfile pMpt(*(TProfile*)hsMpt->Projection(1,0)->ProfileX());
	  pMpt.SetName(Form("pMpt_%i_%i", iCutSetting, iPtCut));
	  pMpt.SetMarkerStyle(20);
	  pMpt.SetMarkerColor(colors[iPtCut%9]);
	  pVecTemp.push_back(pMpt);
	}
	pVecMpt.push_back(pVecTemp);
  }
  
  
  Int_t maxBin = -1;
  Double_t maxPt = -1;
  for(Int_t iCutSetting = 0; iCutSetting < nCutSettings; iCutSetting++) {
	if(GetMaxBin(pVecMpt.at(iCutSetting)) > maxBin)  { maxBin = GetMaxBin(pVecMpt.at(iCutSetting)); }
	if(GetMaxPt(pVecMpt.at(iCutSetting)) > maxPt)  { maxPt = GetMaxPt(pVecMpt.at(iCutSetting)); }
  }
  
  TH1F hDummy("hDummy",";n_{acc} (raw);#LT #it{p}_{T} #GT (raw) (GeV/c)",3500,-0.5,3499.5);
  hDummy.GetXaxis()->SetRangeUser(0, maxBin);
  hDummy.GetYaxis()->SetRangeUser(0, maxPt);
  hDummy.GetYaxis()->SetTitleOffset(1.2);
  hDummy.SetTitleFont(43);
  hDummy.SetTitleSize(18);
  
  TLatex lat;
  lat.SetTextFont(43);
  lat.SetTextSize(16);
  lat.SetTextColor(kBlack);
  
  TLegend leg(0.7, 0.7, 0.9, 0.9);
  leg.SetTextFont(43);
  leg.SetTextSize(15);
  leg.SetTextColor(kBlack);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  for(Int_t iPtCut = 0; iPtCut < nPtCuts; iPtCut++) {
	leg.AddEntry(&pVecMpt.at(0).at(iPtCut), Form("%.2f < #it{p}_{T} < %.2f GeV/c", dVecPtMin.at(iPtCut), dVecPtMax.at(iPtCut)), "P");
  }
  
  
  Double_t dYPositionDefault = 0.38;
  Double_t dDeltaYDefault = 0.04;
  
  gStyle->SetOptStat(0);
  TCanvas c1("c1","c1",nCutSettings*600,600);
  c1.Divide(nCutSettings,1);
  for(Int_t iCutSetting = 0; iCutSetting < nCutSettings; iCutSetting++) {
	c1.cd(iCutSetting+1);
	gPad->SetTopMargin(0.06);
    gPad->SetRightMargin(0.3);

	Double_t dYPosition = dYPositionDefault;
	Double_t dDeltaY = dDeltaYDefault;
	
	hDummy.SetTitle(hCutNames->GetXaxis()->GetBinLabel(iCutSetting+1));
	hDummy.DrawCopy();
	
	for(Int_t iPtCut = 0; iPtCut < nPtCuts; iPtCut++) {
	  pVecMpt.at(iCutSetting).at(iPtCut).Draw("same"); 
	}
	
	for(Int_t iEntry = 0; iEntry < oAddInfo->GetEntries(); iEntry++)
	{
	  if(bDrawEntry.at(iEntry)) {
		TObjString *osAdd = (TObjString*)oAddInfo->At(iEntry);
// 		lat.DrawLatex(10, dYPosition, Form("%s", osAdd->GetString().Data()));
// 		dYPosition -= dDeltaY;
		lat.DrawLatexNDC(0.72, dYPosition, Form("%s", osAdd->GetString().Data()));
		dYPosition -= dDeltaY;
	  }
	}
	
	leg.Draw();
  } // end iCutSetting
  
  
  c1.SaveAs("meanpt_vs_mult_raw.gif");
}

Int_t GetMaxBin(std::vector<TProfile> &vec)
{
  Int_t iMaxBin = -1;
  for(UInt_t iEntry = 0; iEntry < vec.size(); iEntry++) {
	for(Int_t i = 0; i < vec.at(iEntry).GetNbinsX(); i++)
	{
	  if( (vec.at(iEntry).GetBinContent(i)>0) && (vec.at(iEntry).GetBinError(i)>0))
	  {
		iMaxBin = vec.at(iEntry).GetBinCenter(i);
	  }
	}
  }
  return iMaxBin;
}

Double_t GetMaxPt(std::vector<TProfile> &vec)
{
  Double_t dMaxPt = -1;
  for(UInt_t iEntry = 0; iEntry < vec.size(); iEntry++) {
	for(Int_t i = 0; i < vec.at(iEntry).GetNbinsX(); i++)
	{
	  if( ((vec.at(iEntry).GetBinContent(i) + vec.at(iEntry).GetBinError(i)) > dMaxPt) )
	  {
		dMaxPt = (vec.at(iEntry).GetBinContent(i) + vec.at(iEntry).GetBinError(i));
	  }
	}
  }
  return dMaxPt;
}
