#include "TSystem.h"
#include "TGrid.h"

#include "TFile.h"
#include "TString.h"
#include "TMath.h"
#include "TTree.h"

#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"

#include "AliUEHistograms.h"
#include "AliTHn.h"

#include "AliCorrelationAnalysis.h"

AliCorrelationAnalysis::AliCorrelationAnalysis() :
  TObject(),
  fCurrentFileName(""),
  fCacheSameEvent(0x0),
  fCacheMixedEvent(0x0),
  fPtMin(.51),
  fPtMax(49.99),
  fZVtxRange(-1.)
{

}

AliCorrelationAnalysis::~AliCorrelationAnalysis()
{

}

void* AliCorrelationAnalysis::GetUEHistogram(const TString &fileName, TList** listRef, Bool_t mixed, const char* tag)
{
  if ((fCurrentFileName.Length() == 0) ||
      (fCurrentFileName != fileName)) {
    fCurrentFileName = fileName;

    TFile *file = TFile::Open(fCurrentFileName);
    if (!file || file->IsZombie())
      return 0;

    TList *list = (TList*) file->Get("PWG4_LeadingTrackUE/histosLeadingTrackUE");
    if (!list)
      list = (TList*) file->Get(Form("PWG4_PhiCorrelations/histosPhiCorrelations%s", tag));
    if (!list)
      list = (TList*) file->Get("PWG4_PhiCorrelations/histosPhiCorrelations_Syst");

    if (!list)
      return 0;

    if (listRef)
      *listRef = list;

    fCacheMixedEvent = list->FindObject("AliUEHistogramsMixed");
    fCacheSameEvent = list->FindObject("AliUEHistogramsSame");

    if (mixed)
      return fCacheMixedEvent;

    if (list->FindObject("AliUEHistograms"))
      return list->FindObject("AliUEHistograms");

    return fCacheSameEvent;
  } else {
    Printf("GetUEHistogram --> Using cache for %s", fCurrentFileName.Data());

    if (mixed)
      return fCacheMixedEvent;
    else
      return fCacheSameEvent;
  }
}

void AliCorrelationAnalysis::FillParentTHnSparse(TString fileName, Bool_t reduce, const TString &tag)
{
  if (fileName.EndsWith("merge") ||
      fileName.EndsWith("merge_runlist_1") ||
      fileName.EndsWith("merge_runlist_2") ||
      fileName.EndsWith("merge_runlist_3"))
    fileName += "/AnalysisResults.root";

  if (fileName.BeginsWith("alien:"))
    TGrid::Connect("alien:");

  TList* list = 0;

  AliUEHistograms* h = (AliUEHistograms*) GetUEHistogram(fileName, &list, kFALSE, tag);
  Printf("We have %d axes", ((AliTHn*) h->GetUEHist(2)->GetTrackHist(AliUEHist::kToward))->GetNVar());

  if (reduce)
    ((AliTHn*) h->GetUEHist(2)->GetTrackHist(AliUEHist::kToward))->ReduceAxis();
  ((AliTHn*) h->GetUEHist(2)->GetTrackHist(AliUEHist::kToward))->FillParent();
  ((AliTHn*) h->GetUEHist(2)->GetTrackHist(AliUEHist::kToward))->DeleteContainers();

  AliUEHistograms* hMixed = (AliUEHistograms*) GetUEHistogram(fileName, 0, kTRUE, tag);
  if (reduce)
    ((AliTHn*) hMixed->GetUEHist(2)->GetTrackHist(AliUEHist::kToward))->ReduceAxis();
  ((AliTHn*) hMixed->GetUEHist(2)->GetTrackHist(AliUEHist::kToward))->FillParent();
  ((AliTHn*) hMixed->GetUEHist(2)->GetTrackHist(AliUEHist::kToward))->DeleteContainers();

  TString fileNameNew(fileName);

  if (fileName.BeginsWith("alien:"))
    fileNameNew = gSystem->BaseName(fileNameNew);

  fileNameNew.ReplaceAll(".root", "");
  if (reduce)
    fileNameNew += "_.root";
  else
    fileNameNew += "_zvtx.root";

  TFile *file3 = TFile::Open(fileNameNew, "RECREATE");
  file3->mkdir("PWG4_PhiCorrelations");
  file3->cd("PWG4_PhiCorrelations");
  list->Write("histosPhiCorrelations", TObject::kSingleKey);
  file3->Close();
}

void AliCorrelationAnalysis::PlotDeltaPhiEtaGap(const TString &fileNamePbPb, TString fileNamePbPbMix,
						const TString &fileNamepp, const TString &fileNamepp2,
						const TString &outputFile)
{
  if (fileNamePbPbMix.Length() == 0)
    fileNamePbPbMix = fileNamePbPb;

  TFile* file = TFile::Open(outputFile, "RECREATE");
  file->Close();

  Int_t leadingPtOffset = 1;

  Bool_t symmetrizePt = kFALSE;

  // pp HMTF
  Float_t leadingPtArr[] = { 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0 };
  Float_t assocPtArr[] =   { 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0 };
  Int_t maxLeadingPt = sizeof(leadingPtArr) / sizeof(Float_t) - 1;
  Int_t maxAssocPt   = sizeof(assocPtArr)   / sizeof(Float_t) - 1;

  TList* list = 0;
  const char *gTag = "";
  AliUEHistograms *h = (AliUEHistograms*) GetUEHistogram(fileNamePbPb, &list, kFALSE, gTag);
  AliUEHistograms *hMixed = (AliUEHistograms*) GetUEHistogram(fileNamePbPbMix, 0, kTRUE, gTag);

  if (h->GetUEHist(2)->GetTrackHist(AliUEHist::kToward)->GetGrid(6)->GetGrid()->GetNbins() == 0) {
    ((AliTHn*) h->GetUEHist(2)->GetTrackHist(AliUEHist::kToward))->FillParent();
    ((AliTHn*) h->GetUEHist(2)->GetTrackHist(AliUEHist::kToward))->DeleteContainers();
  }

  if (hMixed->GetUEHist(2)->GetTrackHist(AliUEHist::kToward)->GetGrid(6)->GetGrid()->GetNbins() == 0) {
    ((AliTHn*) hMixed->GetUEHist(2)->GetTrackHist(AliUEHist::kToward))->FillParent();
    ((AliTHn*) hMixed->GetUEHist(2)->GetTrackHist(AliUEHist::kToward))->DeleteContainers();
  }

  if (symmetrizePt) {
    h->GetUEHist(2)->SymmetrizepTBins();
    hMixed->GetUEHist(2)->SymmetrizepTBins();
  }

  TList* list2 = 0;
  AliUEHistograms* h2 = 0;
  AliUEHistograms* hMixed2 = 0;
  if (fileNamepp) {
    h2 = (AliUEHistograms*) GetUEHistogram(fileNamepp, &list2);
    hMixed2 = (AliUEHistograms*) GetUEHistogram(fileNamepp, 0, kTRUE);
  }

  TList* list3 = 0;
  AliUEHistograms* h3 = 0;
  AliUEHistograms* hMixed3 = 0;
  if (fileNamepp2) {
    h3 = (AliUEHistograms*) GetUEHistogram(fileNamepp2, &list3);
    hMixed3 = (AliUEHistograms*) GetUEHistogram(fileNamepp2, 0, kTRUE);
  }

  TH2* refMultRaw = (TH2*) list->FindObject("referenceMultiplicity");
  if (refMultRaw) {
    Double_t centrBins[] = { 0., 1., 10., 20., 30., 50., 80., 100. };
    Int_t nCentrBins = sizeof(centrBins) / sizeof(Double_t);
    TH1* refMult = new TH1F("refMult", ";centrality;<Nch>", nCentrBins, centrBins);
    for (Int_t i=0; i<nCentrBins; i++) {
      TH1* proj = refMultRaw->ProjectionY(Form("proj%d", i), refMultRaw->GetXaxis()->FindBin(centrBins[i] + 0.1), refMultRaw->GetXaxis()->FindBin(centrBins[i+1] - 0.1));
      refMult->SetBinContent(refMult->GetXaxis()->FindBin(centrBins[i] + 0.1), proj->GetMean());
      refMult->SetBinError(refMult->GetXaxis()->FindBin(centrBins[i] + 0.1), proj->GetMeanError());
      Printf("Ref multiplicity for centrality %f to %f: %f", centrBins[i], centrBins[i+1], proj->GetMean());
    }
    file = TFile::Open(outputFile, "UPDATE");
    refMult->Write();
    file->Close();
  }

  TTree *tree = (TTree*) list->FindObject("UEAnalysisSettings");

  if (h->GetUEHist(2)->GetTrackEtaCut() == 0) {
    Double_t etaCut = GetEtaCut(tree);
    Printf("Setting eta cut to %f", etaCut);
    h->SetTrackEtaCut(etaCut);
  }

  if (list2) {
    tree = (TTree*) list2->FindObject("UEAnalysisSettings");
    if (tree) {
      Double_t etaCut = GetEtaCut(tree);
      Printf("Setting eta cut to %f", etaCut);
      h2->SetTrackEtaCut(etaCut);
    } else {
      Double_t etaCut = 0.9;
      Printf("WARNING: Setting eta cut to %f without checking", etaCut);
      h2->SetTrackEtaCut(etaCut);
    }
  }

  if (list3) {
    tree = (TTree*) list3->FindObject("UEAnalysisSettings");
    if (tree)	{
      Double_t etaCut = GetEtaCut(tree);
      Printf("Setting eta cut to %f", etaCut);
      h3->SetTrackEtaCut(etaCut);
    } else {
      Double_t etaCut = 0.9;
      Printf("WARNING: Setting eta cut to %f without checking", etaCut);
      h3->SetTrackEtaCut(etaCut);
    }
  }

  Bool_t normalizePerTrigger = kFALSE;
  AliUEHist::CFStep step = (AliUEHist::CFStep) 8;

  for (Int_t i=0; i<maxLeadingPt; i++) {
    for (Int_t j=0; j<maxAssocPt; j++) {

      fPtMin = assocPtArr[j] + 0.01;
      fPtMax = assocPtArr[j+1] - 0.01;

      if(fPtMin >= fPtMax)
	continue;

      SetupRanges(h);
      SetupRanges(hMixed);
      SetupRanges(h2);
      SetupRanges(hMixed2);
      SetupRanges(h3);
      SetupRanges(hMixed3);

      if (assocPtArr[j] >= leadingPtArr[i+leadingPtOffset])
	continue;

      TH1 *hist[8] = { 0x0 };
      Int_t nHist = sizeof(hist) / sizeof(TH1*);

      Bool_t equivMixedBin = kTRUE;
      Bool_t scaleToPairs = kTRUE;
      Int_t histType = 1;

      GetSumOfRatios(h, hMixed, &hist[0],  step,  0,   1, leadingPtArr[i] + 0.01, leadingPtArr[i+leadingPtOffset] - 0.01, normalizePerTrigger);
      GetSumOfRatios(h, hMixed, &hist[1],  step,  0,  10, leadingPtArr[i] + 0.01, leadingPtArr[i+leadingPtOffset] - 0.01, normalizePerTrigger);
      GetSumOfRatios(h, hMixed, &hist[2],  step, 10,  20, leadingPtArr[i] + 0.01, leadingPtArr[i+leadingPtOffset] - 0.01, normalizePerTrigger);
      GetSumOfRatios(h, hMixed, &hist[3],  step, 20,  30, leadingPtArr[i] + 0.01, leadingPtArr[i+leadingPtOffset] - 0.01, normalizePerTrigger);
      GetSumOfRatios(h, hMixed, &hist[4],  step, 30,  50, leadingPtArr[i] + 0.01, leadingPtArr[i+leadingPtOffset] - 0.01, normalizePerTrigger);
      GetSumOfRatios(h, hMixed, &hist[5],  step, 50,  80, leadingPtArr[i] + 0.01, leadingPtArr[i+leadingPtOffset] - 0.01, normalizePerTrigger);
      GetSumOfRatios(h, hMixed, &hist[6],  step, 80, 100, leadingPtArr[i] + 0.01, leadingPtArr[i+leadingPtOffset] - 0.01, normalizePerTrigger);

      if (h2)
	GetSumOfRatios(h2, hMixed2, &hist[7],  step, 0,  -1, leadingPtArr[i] + 0.01, leadingPtArr[i+leadingPtOffset] - 0.01, kTRUE);


      file = TFile::Open(outputFile, "UPDATE");
      for (Int_t iHist = 0; iHist < nHist; ++iHist) {
	if (hist[iHist]) {
	  hist[iHist]->SetName(Form("dphi_%d_%d_%d", i, j, iHist));
	  hist[iHist]->Write();
	}
      }
      file->Close();

      // for (Int_t iHist = 0; iHist < nHist; ++iHist)
      // 	delete hist[i];
    }

    TH1* triggers = h->GetUEHist(2)->GetTriggersAsFunctionOfMultiplicity(step, leadingPtArr[i] + 0.01, leadingPtArr[i+leadingPtOffset] - 0.01);
    triggers->SetName(Form("triggers_%d", i));
    TString str;
    str.Form("%.1f < p_{T,trig} < %.1f", leadingPtArr[i], leadingPtArr[i+leadingPtOffset]);
    triggers->SetTitle(str);

    file = TFile::Open(outputFile, "UPDATE");
    triggers->Write();
    file->Close();
  }

  delete h;
  delete hMixed;
}

Double_t AliCorrelationAnalysis::GetEtaCut(TTree* analysisSettings)
{
  Double_t etaCut = 0;
  if (analysisSettings) {
    analysisSettings->GetBranch("fTrackEtaCut")->SetAddress(&etaCut);
    analysisSettings->GetEntry(0);
  }
  return etaCut;
}

void AliCorrelationAnalysis::SetupRanges(void* obj)
{
  if (!obj)
    return;

  ((AliUEHistograms*) obj)->SetEtaRange(0, 0);
  ((AliUEHistograms*) obj)->SetPtRange(fPtMin, fPtMax);
  ((AliUEHistograms*) obj)->SetCombineMinMax(kTRUE);
  if (fZVtxRange > 0)
    ((AliUEHistograms*) obj)->SetZVtxRange(-fZVtxRange+0.01, fZVtxRange-0.01);
}

void AliCorrelationAnalysis::GetSumOfRatios(void* hVoid, void* hMixedVoid, TH1** hist, AliUEHist::CFStep step,
					    Int_t centralityBegin, Int_t centralityEnd, Float_t ptBegin, Float_t ptEnd,
					    Bool_t normalizePerTrigger, Bool_t useCentralityBinsDirectly)
{
  Printf("GetSumOfRatios | step %d | %d-%d%% | %.1f - %.1f GeV/c | %.1f - %.1f GeV/c", step, centralityBegin, centralityEnd, fPtMin, fPtMax, ptBegin, ptEnd);

  AliUEHistograms *h = (AliUEHistograms*) hVoid;
  AliUEHistograms *hMixed = (AliUEHistograms*) hMixedVoid;

  Int_t centralityBeginBin = 0;
  Int_t centralityEndBin = -1;

  if (!useCentralityBinsDirectly && centralityEnd >= centralityBegin)
    {
      centralityBeginBin = h->GetUEHist(2)->GetEventHist()->GetGrid(step)->GetGrid()->GetAxis(1)->FindBin(0.01 + centralityBegin);
      centralityEndBin = h->GetUEHist(2)->GetEventHist()->GetGrid(step)->GetGrid()->GetAxis(1)->FindBin(-0.01 + centralityEnd);
    }
  else if (useCentralityBinsDirectly)
    {
      centralityBeginBin = centralityBegin;
      centralityEndBin = centralityEnd;
    }

  *hist  = h->GetUEHist(2)->GetSumOfRatios2(hMixed->GetUEHist(2), step, AliUEHist::kToward, ptBegin, ptEnd, centralityBeginBin, centralityEndBin, normalizePerTrigger);

  TString str;
  str.Form("%.1f < p_{T,trig} < %.1f", ptBegin - 0.01, ptEnd + 0.01);

  TString str2;
  str2.Form("%.2f < p_{T,assoc} < %.2f", fPtMin - 0.01, fPtMax + 0.01);

  TString newTitle;
  newTitle.Form("%s - %s - %d-%d", str.Data(), str2.Data(), centralityBegin, centralityEnd);
  if (!useCentralityBinsDirectly)
    newTitle += "%";
  if ((*hist))
    (*hist)->SetTitle(newTitle);
}

void AliCorrelationAnalysis::MergeDPhiFiles(const TString &fileName, const TString &fileName2, const TString &target)
{
  TFile *file = TFile::Open(fileName);
  TFile *file2 = TFile::Open(fileName2);

  TFile *fileTarget = TFile::Open(target, "RECREATE");
  fileTarget->Close();

  Int_t maxLeadingPt = 10;
  Int_t maxAssocPt = 11;

  Int_t nHists = 8;
  for (Int_t i=0; i<maxLeadingPt; i++) {
    TH1* triggers = (TH1*) file->Get(Form("triggers_%d", i));
    if (!triggers)
      continue;

    TH1* triggers2 = (TH1*) ((file2) ? file2->Get(Form("triggers_%d", i)) : 0);
    if (!triggers2)
      Printf("WARNING: trigger %d missing", i);

    for (Int_t j=0; j<maxAssocPt; j++) {
      for (Int_t histId = 0; histId < nHists; histId++) {
	TH2* hist = (TH2*) file->Get(Form("dphi_%d_%d_%d", i, j, histId));
	if (!hist) {
	  TH2* hist2 = (TH2*) ((file2) ? file2->Get(Form("dphi_%d_%d_%d", i, j, histId)) : 0);
	  if (hist2)
	    Printf("WARNING: %d %d %d exists only in file2, not copied!", i, j, histId);
	  continue;
	}

	TString title(hist->GetTitle());
	title.ReplaceAll("%", "");
	TObjArray *tokens = title.Tokenize("-");

	Float_t centralityBegin = ((TObjString*) tokens->At(2))->String().Atoi();
	Float_t centralityEnd = ((TObjString*) tokens->At(3))->String().Atoi();

	Double_t nTriggers = triggers->Integral(triggers->FindBin(centralityBegin + 0.001), triggers->FindBin(centralityEnd - 0.001));
	Double_t nTriggers2 = 0;

	TH2* hist2 = (TH2*) ((file2) ? file2->Get(Form("dphi_%d_%d_%d", i, j, histId)) : 0);
	if (hist2 && triggers2) {
	  hist->Add(hist2);
	  nTriggers2 = triggers2->Integral(triggers2->FindBin(centralityBegin + 0.001), triggers2->FindBin(centralityEnd - 0.001));
	} else {
	  Printf("WARNING: %d %d %d missing", i, j, histId);
	}

	if (nTriggers + nTriggers2 > 0)
	  hist->Scale(1.0 / (nTriggers + nTriggers2));

	Printf("%s %f %f %f %f", hist->GetTitle(), centralityBegin, centralityEnd, nTriggers, nTriggers2);

	fileTarget = TFile::Open(target, "UPDATE");
	hist->Write();
	fileTarget->Close();
      }
    }
  }
}

void AliCorrelationAnalysis::RemoveWing(const TString &fileName, const TString &outputFile)
{
  // remove wing by flattening using the ratio of a flat line to the corr fct at phi = pi +- 1.5 as fct of delta eta

  TFile *file = TFile::Open(fileName);
  TFile *file2 = TFile::Open(outputFile, "RECREATE");
  file2->Close();

  Int_t maxLeadingPt = 10;
  Int_t maxAssocPt = 11;

  TF1* systFunc = new TF1("func", "1.0 + (abs(x) > 0.5) * (abs(x)-0.5) * 9e-4", -2, 2);

  Int_t nHists = 8;
  for (Int_t i=0; i<maxLeadingPt; i++) {
    TH1 *triggers = (TH1*) file->Get(Form("triggers_%d", i));
    if (triggers) {
      file2 = TFile::Open(outputFile, "UPDATE");
      triggers->Write();
      file2->Close();
    }
    for (Int_t j=0; j<maxAssocPt; j++) {
      for (Int_t histId = 0; histId < nHists; histId++) {
	TH2 *hist = (TH2*) file->Get(Form("dphi_%d_%d_%d", i, j, histId));
	if (!hist)
	  continue;

	Float_t width = 1.5;
	TH1* proj = hist->ProjectionY(Form("projx", hist->GetName()), hist->GetXaxis()->FindBin(TMath::Pi() - width), hist->GetXaxis()->FindBin(TMath::Pi()+ width));
	proj->Fit("pol0", "0");
	proj->Divide(proj->GetFunction("pol0"));

	for (Int_t x=1; x<=hist->GetNbinsX(); x++)
	  for (Int_t y=1; y<=hist->GetNbinsY(); y++) {
	    if (proj->GetBinContent(y) <= 0)
	      continue;
	    Double_t divisor = proj->GetBinContent(y);

	    hist->SetBinContent(x, y, hist->GetBinContent(x, y) / divisor);
	    hist->SetBinError(x, y, hist->GetBinError(x, y) / divisor);
	  }

	file2 = TFile::Open(outputFile, "UPDATE");
	hist->Write();
	file2->Close();
      }
    }
  }
}
