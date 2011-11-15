//====================================================================
#include "AliBasedNdetaTask.h"
#include <TMath.h>
#include <TH2D.h>
#include <TH1D.h>
#include <THStack.h>
#include <TList.h>
#include <TParameter.h>
#include <AliAnalysisManager.h>
#include <AliAODEvent.h>
#include <AliAODHandler.h>
#include <AliAODInputHandler.h>
#include "AliForwardUtil.h"
#include "AliAODForwardMult.h"
#include <TFile.h>
#include <TStyle.h>

//____________________________________________________________________
AliBasedNdetaTask::AliBasedNdetaTask()
  : AliAnalysisTaskSE(), 
    fSums(0),		// Container of sums 
    fOutput(0),         // Container of output 
    fVtxMin(0),		// Minimum v_z
    fVtxMax(0),		// Maximum v_z
    fTriggerMask(0),    // Trigger mask 
    fRebin(0),       	// Rebinning factor 
    fCutEdges(false), 
    fSymmetrice(true),
    fCorrEmpty(true), 
    fUseROOTProj(false),
    fTriggerEff(1),
    fShapeCorr(0),
    fListOfCentralities(0),
    fSNNString(0),
    fSysString(0),
    fCent(0),
    fCentAxis(0),
    fNormalizationScheme(kFull), 
    fSchemeString(0), 
    fTriggerString(0),
    fFinalMCCorrFile("")
{
  // 
  // Constructor
  // 
}

//____________________________________________________________________
AliBasedNdetaTask::AliBasedNdetaTask(const char* name)
  : AliAnalysisTaskSE(name), 
    fSums(0),	        // Container of sums 
    fOutput(0),         // Container of output 
    fVtxMin(-10),	// Minimum v_z
    fVtxMax(10),	// Maximum v_z
    fTriggerMask(AliAODForwardMult::kInel), 
    fRebin(5),		// Rebinning factor 
    fCutEdges(false), 
    fSymmetrice(true),
    fCorrEmpty(true), 
    fUseROOTProj(false),
    fTriggerEff(1),
    fShapeCorr(0),
    fListOfCentralities(0),
    fSNNString(0),
    fSysString(0),
    fCent(0),
    fCentAxis(0),
    fNormalizationScheme(kFull), 
    fSchemeString(0),
    fTriggerString(0),
    fFinalMCCorrFile("")
{
  // 
  // Constructor
  // 
  fListOfCentralities = new TObjArray(1);
  
  // Set the normalisation scheme 
  SetNormalizationScheme(kFull);

  // Set the trigger mask
  SetTriggerMask(AliAODForwardMult::kInel);

  // Output slot #1 writes into a TH1 container
  DefineOutput(1, TList::Class()); 
  DefineOutput(2, TList::Class()); 
}

//____________________________________________________________________
AliBasedNdetaTask::AliBasedNdetaTask(const AliBasedNdetaTask& o)
  : AliAnalysisTaskSE(o), 
    fSums(o.fSums),		// TList* - Container of sums 
    fOutput(o.fOutput),         // Container of output 
    fVtxMin(o.fVtxMin),		// Double_t - Minimum v_z
    fVtxMax(o.fVtxMax),		// Double_t - Maximum v_z
    fTriggerMask(o.fTriggerMask),// Int_t - Trigger mask 
    fRebin(o.fRebin),		// Int_t - Rebinning factor 
    fCutEdges(o.fCutEdges),	// Bool_t - Whether to cut edges when rebinning
    fSymmetrice(o.fSymmetrice),
    fCorrEmpty(o.fCorrEmpty), 
    fUseROOTProj(o.fUseROOTProj),
    fTriggerEff(o.fTriggerEff),
    fShapeCorr(o.fShapeCorr),
    fListOfCentralities(o.fListOfCentralities),
    fSNNString(o.fSNNString),
    fSysString(o.fSysString),
    fCent(o.fCent),
    fCentAxis(o.fCentAxis),
    fNormalizationScheme(o.fNormalizationScheme), 
    fSchemeString(o.fSchemeString), 
    fTriggerString(o.fTriggerString),
    fFinalMCCorrFile(o.fFinalMCCorrFile)
{}

//____________________________________________________________________
AliBasedNdetaTask::~AliBasedNdetaTask()
{
  // 
  // Destructor
  // 
  if (fSums) { 
    fSums->Delete();
    delete fSums;
    fSums = 0;
  }
  if (fOutput) { 
    fOutput->Delete();
    delete fOutput;
    fOutput = 0;
  }
}

//________________________________________________________________________
void 
AliBasedNdetaTask::SetCentralityAxis(UShort_t n, Short_t* bins)
{
  if (!fCentAxis) { 
    fCentAxis = new TAxis();
    fCentAxis->SetName("centAxis");
    fCentAxis->SetTitle("Centrality [%]");
  }
  TArrayD dbins(n+1);
  for (UShort_t i = 0; i <= n; i++) 
    dbins[i] = (bins[i] == 100 ? 100.1 : bins[i]);
  fCentAxis->Set(n, dbins.GetArray());
}
    
//________________________________________________________________________
void 
AliBasedNdetaTask::AddCentralityBin(UShort_t at, Short_t low, Short_t high)
{
  // 
  // Add a centrality bin 
  // 
  // Parameters:
  //    low  Low cut
  //    high High cut
  //
  CentralityBin* bin = MakeCentralityBin(GetName(), low, high);
  fListOfCentralities->AddAtAndExpand(bin, at);
}

//________________________________________________________________________
AliBasedNdetaTask::CentralityBin*
AliBasedNdetaTask::MakeCentralityBin(const char* name, 
				     Short_t low, Short_t high) const
{
  // 
  // Make a centrality bin 
  // 
  // Parameters:
  //    name  Name used for histograms
  //    low   Low cut in percent
  //    high  High cut in percent
  // 
  // Return:
  //    A newly created centrality bin 
  //
  return new CentralityBin(name, low, high);
}
//________________________________________________________________________
void 
AliBasedNdetaTask::SetNormalizationScheme(const char* what)
{
  // 
  // Set normalisation scheme 
  // 
  UShort_t    scheme = 0;
  TString     twhat(what);
  twhat.ToUpper();
  TObjString* opt;
  TIter       next(twhat.Tokenize(" ,|"));
  while ((opt = static_cast<TObjString*>(next()))) { 
    TString s(opt->GetString());
    if      (s.IsNull()) continue;
    Bool_t add = true;
    switch (s[0]) { 
    case '-': add = false; // Fall through 
    case '+': s.Remove(0,1);  // Remove character 
    }
    UShort_t bit = 0;
    if      (s.CompareTo("EVENT")     == 0) bit = kEventLevel;
    else if (s.CompareTo("SHAPE")     == 0) bit = kShape;
    else if (s.CompareTo("BACKGROUND")== 0) bit = kBackground;
    else if (s.CompareTo("TRIGGER")   == 0) bit = kTriggerEfficiency;
    else if (s.CompareTo("FULL")      == 0) bit = kFull;
    else if (s.CompareTo("NONE")      == 0) bit = kNone;
    else if (s.CompareTo("ZEROBIN")   == 0) bit = kZeroBin;
    else 
      Warning("SetNormalizationScheme", "Unknown option %s", s.Data());
    if (add) scheme |= bit;
    else     scheme ^= bit;
  }
  SetNormalizationScheme(scheme);
}
//________________________________________________________________________
void 
AliBasedNdetaTask::SetNormalizationScheme(UShort_t scheme) 
{
  fNormalizationScheme = scheme; 
  TString tit = "";
  if (scheme == kFull) tit = "FULL"; 
  else {
    if (scheme & kEventLevel)        tit.Append("EVENT ");
    if (scheme & kShape)             tit.Append("SHAPE ");
    if (scheme & kBackground)        tit.Append("BACKGROUND ");
    if (scheme & kTriggerEfficiency) tit.Append("TRIGGER ");
    if (scheme & kZeroBin)           tit.Append("ZEROBIN ");
  }
  tit = tit.Strip(TString::kBoth);
  if (!fSchemeString) fSchemeString = new TNamed("scheme", "");
  fSchemeString->SetTitle(tit);
  fSchemeString->SetUniqueID(fNormalizationScheme);
}
//________________________________________________________________________
void 
AliBasedNdetaTask::SetTriggerMask(const char* mask)
{
  // 
  // Set the trigger maskl 
  // 
  // Parameters:
  //    mask Trigger mask
  //
  SetTriggerMask(AliAODForwardMult::MakeTriggerMask(mask));
}
//________________________________________________________________________
void 
AliBasedNdetaTask::SetTriggerMask(UShort_t mask) 
{ 
  fTriggerMask = mask; 
  TString tit(AliAODForwardMult::GetTriggerString(mask));
  tit = tit.Strip(TString::kBoth);
  if (!fTriggerString) fTriggerString = new TNamed("trigger", "");
  fTriggerString->SetTitle(tit);
  fTriggerString->SetUniqueID(fTriggerMask);
}

//________________________________________________________________________
void 
AliBasedNdetaTask::SetShapeCorrection(const TH1* c)
{
  // 
  // Set the shape correction (a.k.a., track correction) for selected
  // trigger(s)
  // 
  // Parameters:
  //    h Correction
  //
  if (!c) return;
  fShapeCorr = static_cast<TH1*>(c->Clone());
  fShapeCorr->SetDirectory(0);
}

//________________________________________________________________________
void 
AliBasedNdetaTask::UserCreateOutputObjects()
{
  // 
  // Create output objects.  
  //
  // This is called once per slave process 
  //
  fSums = new TList;
  fSums->SetName(Form("%s_sums", GetName()));
  fSums->SetOwner();

  // Automatically add 'all' centrality bin if nothing has been defined. 
  AddCentralityBin(0, 0, 0);
  if (fCentAxis && fCentAxis->GetNbins() > 0 && fCentAxis->GetXbins()) { 
    const TArrayD* bins = fCentAxis->GetXbins();
    Int_t          nbin = fCentAxis->GetNbins(); 
    for (Int_t i = 0; i < nbin; i++) 
      AddCentralityBin(i+1,  Short_t((*bins)[i]), Short_t((*bins)[i+1]));
  }
  if (fCentAxis) fSums->Add(fCentAxis);


  // Centrality histogram 
  fCent = new TH1D("cent", "Centrality", 100, 0, 100);
  fCent->SetDirectory(0);
  fCent->SetXTitle(0);
  fSums->Add(fCent);

  // Loop over centrality bins 
  TIter next(fListOfCentralities);
  CentralityBin* bin = 0;
  while ((bin = static_cast<CentralityBin*>(next()))) 
    bin->CreateOutputObjects(fSums);

  // Check that we have an AOD input handler 
  AliAnalysisManager* am = AliAnalysisManager::GetAnalysisManager();
  AliAODInputHandler* ah = 
    dynamic_cast<AliAODInputHandler*>(am->GetInputEventHandler());
  if (!ah) AliFatal("No AOD input handler set in analysis manager");

  // Post data for ALL output slots >0 here, to get at least an empty histogram
  PostData(1, fSums); 
}
//____________________________________________________________________
void 
AliBasedNdetaTask::UserExec(Option_t *) 
{
  // 
  // Process a single event 
  // 
  // Parameters:
  //    option Not used
  //
  // Main loop
  AliAODEvent* aod = dynamic_cast<AliAODEvent*>(InputEvent());
  if (!aod) {
    AliError("Cannot get the AOD event");
    return;
  }  
  
  TObject* obj = aod->FindListObject("Forward");
  if (!obj) { 
    AliWarning("No forward object found");
    return;
  }
  AliAODForwardMult* forward = static_cast<AliAODForwardMult*>(obj);
  
  // Fill centrality histogram 
  Float_t cent = forward->GetCentrality();
  fCent->Fill(cent);

  // Get the histogram(s) 
  TH2D* data   = GetHistogram(aod, false);
  TH2D* dataMC = GetHistogram(aod, true);

  Bool_t isZero = ((fNormalizationScheme & kZeroBin) &&
		   !forward->IsTriggerBits(AliAODForwardMult::kNClusterGt0));


  // Loop over centrality bins 
  CentralityBin* allBin = 
    static_cast<CentralityBin*>(fListOfCentralities->At(0));
  allBin->ProcessEvent(forward, fTriggerMask, isZero, 
		       fVtxMin, fVtxMax, data, dataMC);
  
  // Find this centrality bin 
  if (fCentAxis && fCentAxis->GetNbins() > 0) {
    Int_t          icent   = fCentAxis->FindBin(cent);
    CentralityBin* thisBin = 0;
    if (icent >= 1 && icent <= fCentAxis->GetNbins()) 
      thisBin = static_cast<CentralityBin*>(fListOfCentralities->At(icent));
    if (thisBin)
      thisBin->ProcessEvent(forward, fTriggerMask, isZero, fVtxMin, fVtxMax, 
			    data, dataMC);
  }

  // Here, we get the update 
  if (!fSNNString) { 
    UShort_t sNN = forward->GetSNN();
    fSNNString = new TNamed("sNN", "");
    fSNNString->SetTitle(AliForwardUtil::CenterOfMassEnergyString(sNN));
    fSNNString->SetUniqueID(sNN);
    fSums->Add(fSNNString);
  
    UShort_t sys = forward->GetSystem();
    fSysString = new TNamed("sys", "");
    fSysString->SetTitle(AliForwardUtil::CollisionSystemString(sys));
    fSysString->SetUniqueID(sys);
    fSums->Add(fSysString);

    fSums->Add(fSchemeString);
    fSums->Add(fTriggerString);

    // Print();
  }
  
  PostData(1, fSums);
}

//________________________________________________________________________
void 
AliBasedNdetaTask::SetHistogramAttributes(TH1D* h, Int_t colour, Int_t marker,
					  const char* title, const char* ytitle)
{
  // 
  // Set histogram graphical options, etc. 
  // 
  // Parameters:
  //    h       Histogram to modify
  //    colour  Marker color 
  //    marker  Marker style
  //    title   Title of histogram
  //    ytitle  Title on y-axis. 
  //
  h->SetTitle(title);
  h->SetMarkerColor(colour);
  h->SetMarkerStyle(marker);
  h->SetMarkerSize(marker == 29 || marker == 30 ? 1.2 : 1);
  h->SetFillStyle(0);
  h->SetYTitle(ytitle);
  h->GetXaxis()->SetTitleFont(132);
  h->GetXaxis()->SetLabelFont(132);
  h->GetXaxis()->SetNdivisions(10);
  h->GetYaxis()->SetTitleFont(132);
  h->GetYaxis()->SetLabelFont(132);
  h->GetYaxis()->SetNdivisions(10);
  h->GetYaxis()->SetDecimals();
  h->SetStats(0);
}

//________________________________________________________________________
void
AliBasedNdetaTask::ScaleToCoverage(TH2D* copy, const TH1D* norm) 
{
  // Normalize to the acceptance -
  // dndeta->Divide(accNorm);
  for (Int_t i = 1; i <= copy->GetNbinsX(); i++) { 
    Double_t a = norm->GetBinContent(i);
    for (Int_t j = 1; j <= copy->GetNbinsY(); j++) { 
      if (a <= 0) { 
	copy->SetBinContent(i,j,0);
	copy->SetBinError(i,j,0);
	continue;
      }
      Double_t c = copy->GetBinContent(i, j);
      Double_t e = copy->GetBinError(i, j);
      copy->SetBinContent(i, j, c / a);
      copy->SetBinError(i, j, e / a);
    }
  }
}

//________________________________________________________________________
TH1D*
AliBasedNdetaTask::ProjectX(const TH2D* h, 
			    const char* name,
			    Int_t firstbin, 
			    Int_t lastbin, 
			    bool  useRoot,
			    bool  corr,
			    bool  error)
{
  // 
  // Project onto the X axis 
  // 
  // Parameters:
  //    h         2D histogram 
  //    name      New name 
  //    firstbin  First bin to use 
  //    lastbin   Last bin to use
  //    error     Whether to calculate errors
  // 
  // Return:
  //    Newly created histogram or null
  //
  if (!h) return 0;
  if (useRoot) 
    return h->ProjectionX(name, firstbin, lastbin, (error ? "e" : ""));
  
  TAxis* xaxis = h->GetXaxis();
  TAxis* yaxis = h->GetYaxis();
  TH1D*  ret   = new TH1D(name, h->GetTitle(), xaxis->GetNbins(), 
			  xaxis->GetXmin(), xaxis->GetXmax());
  static_cast<const TAttLine*>(h)->Copy(*ret);
  static_cast<const TAttFill*>(h)->Copy(*ret);
  static_cast<const TAttMarker*>(h)->Copy(*ret);
  ret->GetXaxis()->ImportAttributes(xaxis);

  Int_t first = firstbin; 
  Int_t last  = lastbin;
  if (first < 0)                         first = 0;
  else if (first >= yaxis->GetNbins()+1) first = yaxis->GetNbins();
  if      (last  < 0)                    last  = yaxis->GetNbins();
  else if (last  >  yaxis->GetNbins()+1) last  = yaxis->GetNbins();
  if (last-first < 0) { 
    AliWarningGeneral("AliBasedNdetaTask", 
		      Form("Nothing to project [%d,%d]", first, last));
    return 0;
    
  }

  // Loop over X bins 
  // AliInfo(Form("Projecting bins [%d,%d] of %s", first, last, h->GetName()));
  Int_t ybins = (last-first+1);
  for (Int_t xbin = 0; xbin <= xaxis->GetNbins()+1; xbin++) { 
    Double_t content = 0;
    Double_t error2  = 0;
    Int_t    nbins   = 0;
    
    
    for (Int_t ybin = first; ybin <= last; ybin++) { 
      Double_t c1 = h->GetCellContent(xbin, ybin);
      Double_t e1 = h->GetCellError(xbin, ybin);

      // Ignore empty bins 
      if (c1 < 1e-12) continue;
      if (e1 < 1e-12) {
	if (error) continue; 
	e1 = 1;
      }

      content    += c1;
      error2     += e1*e1;
      nbins++;
    } // for (ybin)
    if(content > 0 && nbins > 0) {
      Double_t factor = (corr ? Double_t(ybins) / nbins : 1);
#if 0
      AliWarningGeneral(ret->GetName(), 
			Form("factor @ %d is %d/%d -> %f", 
			     xbin, ybins, nbins, factor));
#endif
      if (error) { 
	// Calculate weighted average
	ret->SetBinContent(xbin, content * factor);
	ret->SetBinError(xbin, factor * TMath::Sqrt(error2));
      }
      else 
	ret->SetBinContent(xbin, factor * content);
    }
  } // for (xbin)
  
  return ret;
}
  
//________________________________________________________________________
void 
AliBasedNdetaTask::Terminate(Option_t *) 
{
  // 
  // Called at end of event processing.. 
  //
  // This is called once in the master 
  // 
  // Parameters:
  //    option Not used 
  //
  // Draw result to screen, or perform fitting, normalizations Called
  // once at the end of the query

  fSums = dynamic_cast<TList*> (GetOutputData(1));
  if(!fSums) {
    AliError("Could not retrieve TList fSums"); 
    return; 
  }
  
  fOutput = new TList;
  fOutput->SetName(Form("%s_result", GetName()));
  fOutput->SetOwner();
  
  fSNNString     = static_cast<TNamed*>(fSums->FindObject("sNN"));
  fSysString     = static_cast<TNamed*>(fSums->FindObject("sys"));
  fCentAxis      = static_cast<TAxis*>(fSums->FindObject("centAxis"));
  fSchemeString  = static_cast<TNamed*>(fSums->FindObject("scheme"));
  fTriggerString = static_cast<TNamed*>(fSums->FindObject("trigger"));

  if(fSysString && fSNNString && 
     fSysString->GetUniqueID() == AliForwardUtil::kPP)
    LoadNormalizationData(fSysString->GetUniqueID(),
			  fSNNString->GetUniqueID());
  // Print before we loop
  Print();

  // Loop over centrality bins 
  TIter next(fListOfCentralities);
  CentralityBin* bin = 0;
  gStyle->SetPalette(1);
  THStack* dndetaStack        = new THStack("dndeta", "dN/d#eta");
  THStack* dndetaStackRebin   = new THStack(Form("dndeta_rebin%02d", fRebin), 
					    "dN_{ch}/d#eta");
  THStack* dndetaMCStack      = new THStack("dndetaMC", "dN_{ch}/d#eta");
  THStack* dndetaMCStackRebin = new THStack(Form("dndetaMC_rebin%02d", fRebin), 
					    "dN_{ch}/d#eta");
  
  TList* mclist = 0;
  TList* truthlist = 0;
  
  if (fFinalMCCorrFile.Contains(".root")) {
    TFile* ftest = TFile::Open(fFinalMCCorrFile.Data());
    if(ftest) {
      mclist    = dynamic_cast<TList*>(ftest->Get(Form("%sResults",GetName())));
      truthlist = dynamic_cast<TList*>(ftest->Get("MCTruthResults"));
    }
    else 
      AliWarning("MC analysis file invalid - no final MC correction possible");
  }
  Int_t style = GetMarker();
  Int_t color = GetColor();
  
  AliInfo(Form("Marker style=%d, color=%d", style, color));
  while ((bin = static_cast<CentralityBin*>(next()))) {
    
    bin->End(fSums, fOutput, fNormalizationScheme, fShapeCorr, fTriggerEff,
	     fSymmetrice, fRebin, fUseROOTProj, fCorrEmpty, fCutEdges, 
	     fTriggerMask, style, color, mclist, truthlist);
    if (fCentAxis && bin->IsAllBin()) continue;
    TH1* dndeta      = bin->GetResult(0, false, "");
    TH1* dndetaSym   = bin->GetResult(0, true,  "");
    TH1* dndetaMC    = bin->GetResult(0, false, "MC");
    TH1* dndetaMCSym = bin->GetResult(0, true,  "MC");
    if (dndeta)      dndetaStack->Add(dndeta);
    if (dndetaSym)   dndetaStack->Add(dndetaSym);
    if (dndetaMC)    dndetaMCStack->Add(dndetaMC);
    if (dndetaMCSym) dndetaMCStack->Add(dndetaMCSym);
    if (fRebin > 1) { 
      dndeta      = bin->GetResult(fRebin, false, "");
      dndetaSym   = bin->GetResult(fRebin, true,  "");
      dndetaMC    = bin->GetResult(fRebin, false, "MC");
      dndetaMCSym = bin->GetResult(fRebin, true,  "MC");
      if (dndeta)      dndetaStackRebin->Add(dndeta);
      if (dndetaSym)   dndetaStackRebin->Add(dndetaSym);
      if (dndetaMC)    dndetaMCStackRebin->Add(dndetaMC);
      if (dndetaMCSym) dndetaMCStackRebin->Add(dndetaMCSym);
    }
  }
  // Output the stack
  fOutput->Add(dndetaStack);

  // If available output rebinned stack 
  if (!dndetaStackRebin->GetHists() || 
      dndetaStackRebin->GetHists()->GetEntries() <= 0) {
    AliWarning("No rebinned histograms found");
    delete dndetaStackRebin;
    dndetaStackRebin = 0;
  }
  if (dndetaStackRebin) fOutput->Add(dndetaStackRebin);

  // If available, output track-ref stack
  if (!dndetaMCStack->GetHists() || 
      dndetaMCStack->GetHists()->GetEntries() <= 0) {
    AliWarning("No MC histograms found");
    delete dndetaMCStack;
    dndetaMCStack = 0;
  }
  if (dndetaMCStack) fOutput->Add(dndetaMCStack);

  // If available, output rebinned track-ref stack
  if (!dndetaMCStackRebin->GetHists() || 
      dndetaMCStackRebin->GetHists()->GetEntries() <= 0) {
    AliWarning("No rebinned MC histograms found");
    delete dndetaMCStackRebin;
    dndetaMCStackRebin = 0;
  }
  if (dndetaMCStackRebin) fOutput->Add(dndetaMCStackRebin);

  // Output collision energy string 
  if (fSNNString) fOutput->Add(fSNNString->Clone());

  // Output collision system string 
  if (fSysString) fOutput->Add(fSysString->Clone());

  // Output centrality axis 
  if (fCentAxis) fOutput->Add(fCentAxis);

  // Output trigger string 
  if (fTriggerString) fOutput->Add(fTriggerString->Clone());
  
  // Normalization string 
  if (fSchemeString) fOutput->Add(fSchemeString->Clone());

  // Output vertex axis 
  TAxis* vtxAxis = new TAxis(1,fVtxMin,fVtxMax);
  vtxAxis->SetName("vtxAxis");
  vtxAxis->SetTitle(Form("v_{z}#in[%+5.1f,%+5.1f]cm", fVtxMin,fVtxMax));
  fOutput->Add(vtxAxis);

  // Output trigger efficiency and shape correction 
  fOutput->Add(new TParameter<Double_t>("triggerEff", fTriggerEff)); 
  if (fShapeCorr) fOutput->Add(fShapeCorr);

  TNamed* options = new TNamed("options","");
  TString str;
  str.Append(Form("Edges %scut, ", fCutEdges ? "" : "not "));
  str.Append(Form("Empty bins %scorrected for, ", fCorrEmpty ? "" : "not "));
  str.Append(Form("TH2::ProjectionX %sused", fUseROOTProj ? "" : "not "));
  options->SetTitle(str);
  fOutput->Add(options);

  PostData(2, fOutput);
}
//________________________________________________________________________
void
AliBasedNdetaTask::LoadNormalizationData(UShort_t sys, UShort_t energy)
{
  // Load the normalisation data for dN/deta for pp INEL, INEL>0, and NSD
  TString type("pp");
  TString snn("900");
  if(energy == 7000) snn.Form("7000");
  if(energy == 2750) snn.Form("2750"); 
  
  if(fShapeCorr &&  (fTriggerEff != 1)) {
    AliInfo("Objects already set for normalization - no action taken"); 
    return; 
  }
  
  TFile* fin = TFile::Open(Form("$ALICE_ROOT/PWG2/FORWARD/corrections/"
				"Normalization/normalizationHists_%s_%s.root",
				type.Data(),snn.Data()));
  if(!fin) {
    AliWarning(Form("no file for normalization of %d/%d", sys, energy));
    return;
  }

  // Shape correction
  if ((fNormalizationScheme & kShape) && !fShapeCorr) {
    TString trigName("All");
    if (fTriggerMask == AliAODForwardMult::kInel || 
	fTriggerMask == AliAODForwardMult::kNClusterGt0) 
      trigName = "Inel";
    else if (fTriggerMask == AliAODForwardMult::kNSD)
      trigName = "NSD";
    else if (fTriggerMask == AliAODForwardMult::kInelGt0)
      trigName = "InelGt0";
    else {
      AliWarning(Form("Normalization for trigger %s not known, using all",
		      AliAODForwardMult::GetTriggerString(fTriggerMask)));
    }
      
    TString shapeCorName(Form("h%sNormalization", trigName.Data()));
    TH2F*   shapeCor = dynamic_cast<TH2F*>(fin->Get(shapeCorName));
    if (shapeCor) SetShapeCorrection(shapeCor);
    else { 
      AliWarning(Form("No shape correction found for %s", trigName.Data()));
    }
  }

  // Trigger efficiency
  TString effName(Form("%sTriggerEff", 
		       fTriggerMask == AliAODForwardMult::kInel ? "inel" :
		       fTriggerMask == AliAODForwardMult::kNSD ? "nsd" :
		       fTriggerMask == AliAODForwardMult::kInelGt0 ?
		       "inelgt0" : "all"));
  TParameter<float>* eff = 0;
  if (fNormalizationScheme & kTriggerEfficiency) 
    eff = static_cast<TParameter<float>*>(fin->Get(effName));
  Double_t trigEff = eff ? eff->GetVal() : 1;
  if (fTriggerEff != 1) SetTriggerEff(trigEff);
  if (fTriggerEff < 0)  fTriggerEff = 1;

  // TEMPORARY FIX
  // Rescale the shape correction by the trigger efficiency 
  if (fShapeCorr) {
    AliWarning(Form("Rescaling shape correction by trigger efficiency: "
		    "1/E_X=1/%f", trigEff));
    fShapeCorr->Scale(1. / trigEff);
  }

  // Print - out
  if (fShapeCorr && fTriggerEff) AliInfo("Loaded objects for normalization.");
}


//________________________________________________________________________
void 
AliBasedNdetaTask::Print(Option_t*) const
{
  // 
  // Print information 
  // 
  std::cout << this->ClassName() << ": " << this->GetName() << "\n"
	    << std::boolalpha 
	    << " Trigger:              " << (fTriggerString ? 
					     fTriggerString->GetTitle() :
					     "none") << "\n"
	    << " Vertex range:         [" << fVtxMin << ":" << fVtxMax << "]\n"
	    << " Rebin factor:         " << fRebin << "\n" 
	    << " Cut edges:            " << fCutEdges << "\n"
	    << " Symmertrice:          " << fSymmetrice << "\n"
	    << " Use TH2::ProjectionX: " << fUseROOTProj << "\n"
	    << " Correct for empty:    " << fCorrEmpty << "\n"
	    << " Normalization scheme: " << (fSchemeString ? 
					     fSchemeString->GetTitle() : 
					     "none") <<"\n"
	    << " Trigger efficiency:   " << fTriggerEff << "\n" 
	    << " Shape correction:     " << (fShapeCorr ? 
					   fShapeCorr->GetName() : 
					   "none") << "\n"
	    << " sqrt(s_NN):           " << (fSNNString ? 
					     fSNNString->GetTitle() : 
					   "unknown") << "\n"
	    << " Collision system:     " << (fSysString ? 
					     fSysString->GetTitle() : 
					   "unknown") << "\n"
	    << " Centrality bins:      " << (fCentAxis ? "" : "none");
  if (fCentAxis) { 
    Int_t           nBins = fCentAxis->GetNbins();
    const Double_t* bins  = fCentAxis->GetXbins()->GetArray();
    for (Int_t i = 0; i <= nBins; i++) 
      std::cout << (i==0 ? "" : "-") << bins[i];
  }
  std::cout << std::noboolalpha << std::endl;
  
}

//________________________________________________________________________
TH1D*
AliBasedNdetaTask::Rebin(const TH1D* h, Int_t rebin, Bool_t cutEdges)
{
  // 
  // Make a copy of the input histogram and rebin that histogram
  // 
  // Parameters:
  //    h  Histogram to rebin
  // 
  // Return:
  //    New (rebinned) histogram
  //
  if (rebin <= 1) return 0;

  Int_t nBins = h->GetNbinsX();
  if(nBins % rebin != 0) {
    AliWarningGeneral("AliBasedNdetaTask", 
		      Form("Rebin factor %d is not a devisor of current number "
			   "of bins %d in the histogram %s", 
			   rebin, nBins, h->GetName()));
    return 0;
  }
    
  // Make a copy 
  TH1D* tmp = static_cast<TH1D*>(h->Clone(Form("%s_rebin%02d", 
					       h->GetName(), rebin)));
  tmp->Rebin(rebin);
  tmp->SetDirectory(0);

  // The new number of bins 
  Int_t nBinsNew = nBins / rebin;
  for(Int_t i = 1;i<= nBinsNew; i++) {
    Double_t content = 0;
    Double_t sumw    = 0;
    Double_t wsum    = 0;
    Int_t    nbins   = 0;
    for(Int_t j = 1; j<=rebin;j++) {
      Int_t    bin = (i-1)*rebin + j;
      Double_t c   =  h->GetBinContent(bin);
      if (c <= 0) continue;
      
      if (cutEdges) {
	if (h->GetBinContent(bin+1)<=0 || 
	    h->GetBinContent(bin-1)<=0) {
#if 0
	  AliWarningGeneral("AliBasedNdetaTask", 
			    Form("removing bin %d=%f of %s (%d=%f,%d=%f)", 
				 bin, c, h->GetName(), 
				 bin+1, h->GetBinContent(bin+1), 
				 bin-1, h->GetBinContent(bin-1)));
#endif
	  continue;
	}	
      }
      Double_t e =  h->GetBinError(bin);
      Double_t w =  1 / (e*e); // 1/c/c
      content    += c;
      sumw       += w;
      wsum       += w * c;
      nbins++;
    }
      
    if(content > 0 && nbins > 0) {
      tmp->SetBinContent(i, wsum / sumw);
      tmp->SetBinError(i,1./TMath::Sqrt(sumw));
    }
  }
  
  return tmp;
}

//__________________________________________________________________
TH1* 
AliBasedNdetaTask::Symmetrice(const TH1* h)
{
  // 
  // Make an extension of @a h to make it symmetric about 0 
  // 
  // Parameters:
  //    h Histogram to symmertrice 
  // 
  // Return:
  //    Symmetric extension of @a h 
  //
  Int_t nBins = h->GetNbinsX();
  TH1* s = static_cast<TH1*>(h->Clone(Form("%s_mirror", h->GetName())));
  s->SetTitle(Form("%s (mirrored)", h->GetTitle()));
  s->Reset();
  s->SetBins(nBins, -h->GetXaxis()->GetXmax(), -h->GetXaxis()->GetXmin());
  s->SetMarkerColor(h->GetMarkerColor());
  s->SetMarkerSize(h->GetMarkerSize());
  s->SetMarkerStyle(FlipHollowStyle(h->GetMarkerStyle()));
  s->SetFillColor(h->GetFillColor());
  s->SetFillStyle(h->GetFillStyle());
  s->SetDirectory(0);

  // Find the first and last bin with data 
  Int_t first = nBins+1;
  Int_t last  = 0;
  for (Int_t i = 1; i <= nBins; i++) { 
    if (h->GetBinContent(i) <= 0) continue; 
    first = TMath::Min(first, i);
    last  = TMath::Max(last,  i);
  }
    
  Double_t xfirst = h->GetBinCenter(first-1);
  Int_t    f1     = h->GetXaxis()->FindBin(-xfirst);
  Int_t    l2     = s->GetXaxis()->FindBin(xfirst);
  for (Int_t i = f1, j=l2; i <= last; i++,j--) { 
    s->SetBinContent(j, h->GetBinContent(i));
    s->SetBinError(j, h->GetBinError(i));
  }
  // Fill in overlap bin 
  s->SetBinContent(l2+1, h->GetBinContent(first));
  s->SetBinError(l2+1, h->GetBinError(first));
  return s;
}

//__________________________________________________________________
Int_t
AliBasedNdetaTask::GetMarkerStyle(UShort_t bits)
{
  Int_t  base   = bits & (0xFE);
  Bool_t hollow = bits & kHollow;
  switch (base) { 
  case kCircle:       return (hollow ? 24 : 20);
  case kSquare:       return (hollow ? 25 : 21);
  case kUpTriangle:   return (hollow ? 26 : 22);
  case kDownTriangle: return (hollow ? 32 : 23);
  case kDiamond:      return (hollow ? 27 : 33); 
  case kCross:        return (hollow ? 28 : 34); 
  case kStar:         return (hollow ? 30 : 29); 
  }
  return 1;
}
//__________________________________________________________________
UShort_t
AliBasedNdetaTask::GetMarkerBits(Int_t style) 
{ 
  UShort_t bits = 0;
  switch (style) { 
  case 24: case 25: case 26: case 27: case 28: case 30: case 32: 
    bits |= kHollow; break;
  }
  switch (style) { 
  case 20: case 24: bits |= kCircle;       break;
  case 21: case 25: bits |= kSquare;       break;
  case 22: case 26: bits |= kUpTriangle;   break;
  case 23: case 32: bits |= kDownTriangle; break;
  case 27: case 33: bits |= kDiamond;      break;
  case 28: case 34: bits |= kCross;        break;
  case 29: case 30: bits |= kStar;         break;
  }
  return bits;
}
//__________________________________________________________________
Int_t
AliBasedNdetaTask::FlipHollowStyle(Int_t style) 
{
  UShort_t bits = GetMarkerBits(style);
  Int_t    ret  = GetMarkerStyle(bits ^ kHollow);
  return ret;
}

//====================================================================
void
AliBasedNdetaTask::Sum::Init(TList* list, const TH2D* data, Int_t col)
{
  TString n(GetHistName(0));
  TString n0(GetHistName(1));
  const char* postfix = GetTitle();

  fSum = static_cast<TH2D*>(data->Clone(n));
  if (postfix) fSum->SetTitle(Form("%s (%s)", data->GetTitle(), postfix));
  fSum->SetDirectory(0);
  fSum->SetMarkerColor(col);
  fSum->SetMarkerStyle(GetMarkerStyle(kCircle|kSolid));
  fSum->Reset();
  list->Add(fSum);

  fSum0 = static_cast<TH2D*>(data->Clone(n0));
  if (postfix) 
    fSum0->SetTitle(Form("%s 0-bin (%s)", data->GetTitle(), postfix));
  else   
    fSum0->SetTitle(Form("%s 0-bin", data->GetTitle()));
  fSum0->SetDirectory(0);
  fSum0->SetMarkerColor(col);
  fSum0->SetMarkerStyle(GetMarkerStyle(kCross|kHollow));
  fSum0->Reset();
  list->Add(fSum0);

  fEvents = new TH1I(GetHistName(2), "Event types", 2, -.5, 1.5);
  fEvents->SetDirectory(0);
  fEvents->GetXaxis()->SetBinLabel(1, "Non-zero");
  fEvents->GetXaxis()->SetBinLabel(2, "Zero");
  list->Add(fEvents);
}

//____________________________________________________________________
TString
AliBasedNdetaTask::Sum::GetHistName(Int_t what) const
{
  TString n(GetName());
  if      (what == 1) n.Append("0");
  else if (what == 2) n.Append("Events");
  const char* postfix = GetTitle();
  if (postfix && postfix[0] != '\0')  n.Append(postfix);
  return n;
}

//____________________________________________________________________
void
AliBasedNdetaTask::Sum::Add(const TH2D* data, Bool_t isZero)
{

  if (isZero) fSum0->Add(data);
  else        fSum->Add(data);
  fEvents->Fill(isZero ? 1 : 0);
}

//____________________________________________________________________
TH2D*
AliBasedNdetaTask::Sum::GetSum(const TList* input, 
			       TList*       output, 
			       Double_t&    ntotal,
			       Double_t     epsilon0, 
			       Double_t     epsilon,
			       Int_t        marker,
			       Bool_t       rootProj, 
			       Bool_t       corrEmpty) const
{
  TH2D* sum      = static_cast<TH2D*>(input->FindObject(GetHistName(0)));
  TH2D* sum0     = static_cast<TH2D*>(input->FindObject(GetHistName(1)));
  TH1I* events   = static_cast<TH1I*>(input->FindObject(GetHistName(2)));
  if (!sum || !sum0 || !events) {
    AliWarning(Form("Failed to find one or more histograms: "
		    "%s (%p) %s (%p) %s (%p)", 
		    GetHistName(0).Data(), sum, 
		    GetHistName(1).Data(), sum0,
		    GetHistName(2).Data(), events));
    return 0;
  }

  TH2D* ret      = static_cast<TH2D*>(sum->Clone(sum->GetName()));
  ret->SetDirectory(0);
  ret->Reset();
  Int_t n        = Int_t(events->GetBinContent(1));
  Int_t n0       = Int_t(events->GetBinContent(2));

  AliInfo(Form("Adding histograms %s and %s with weights %f and %f resp.",
	       sum0->GetName(), sum->GetName(), 1./epsilon, 1./epsilon0));
  // Generate merged histogram 
  ret->Add(sum0, sum, 1. / epsilon0, 1. / epsilon); 
  ntotal = n / epsilon + n0 / epsilon0;

  TList* out = new TList;
  out->SetOwner();
  const char* postfix = GetTitle();
  if (!postfix) postfix = "";
  out->SetName(Form("partial%s", postfix));
  output->Add(out);

  // Now make copies, normalize them, and store in output list 
  TH2D* sumCopy  = static_cast<TH2D*>(sum->Clone("sum"));
  TH2D* sum0Copy = static_cast<TH2D*>(sum0->Clone("sum0"));
  TH2D* retCopy  = static_cast<TH2D*>(ret->Clone("sumAll"));
  sumCopy->SetMarkerStyle(FlipHollowStyle(marker));
  sumCopy->SetDirectory(0);
  sum0Copy->SetMarkerStyle(GetMarkerStyle(GetMarkerBits(marker)+4));
  sum0Copy->SetDirectory(0);
  retCopy->SetMarkerStyle(marker);
  retCopy->SetDirectory(0);

  TH1D* norm    = ProjectX(sum,  "norm",    0, 0, rootProj, corrEmpty, false);
  TH1D* norm0   = ProjectX(sum0, "norm0",   0, 0, rootProj, corrEmpty, false);
  TH1D* normAll = ProjectX(ret,  "normAll", 0, 0, rootProj, corrEmpty, false);
  norm->SetDirectory(0);
  norm0->SetDirectory(0);
  normAll->SetDirectory(0);
  
  ScaleToCoverage(sumCopy, norm);
  ScaleToCoverage(sum0Copy, norm0);
  ScaleToCoverage(retCopy, normAll);

  Int_t nY = sum->GetNbinsY();
  TH1D* sumCopyPx  = ProjectX(sumCopy,  "average",    1, nY,rootProj,corrEmpty);
  TH1D* sum0CopyPx = ProjectX(sum0Copy, "average0",   1, nY,rootProj,corrEmpty);
  TH1D* retCopyPx  = ProjectX(retCopy,  "averageAll", 1, nY,rootProj,corrEmpty);
  sumCopyPx->SetDirectory(0);
  sum0CopyPx->SetDirectory(0);
  retCopyPx->SetDirectory(0);

  // Scale our 1D histograms
  sumCopyPx->Scale(1., "width");
  sum0CopyPx->Scale(1., "width");  
  retCopyPx->Scale(1., "width");  

  AliInfo(Form("Maximum %f,%f changed to %f", sumCopyPx->GetMaximum(), 
	       sum0CopyPx->GetMaximum(), retCopyPx->GetMaximum()));

  // Scale the normalization - they should be 1 at the maximum
  norm->Scale(n > 0   ? 1. / n  : 1);
  norm0->Scale(n0 > 0 ? 1. / n0 : 1);
  normAll->Scale(ntotal > 0 ? 1. / ntotal : 1);

  out->Add(sumCopy);
  out->Add(sum0Copy);
  out->Add(retCopy);
  out->Add(sumCopyPx);
  out->Add(sum0CopyPx);
  out->Add(retCopyPx);
  out->Add(norm);
  out->Add(norm0);
  out->Add(normAll);

  AliInfo(Form("Returning  (1/%f * %s + 1/%f * %s), "
	       "1/%f * %d + 1/%f * %d = %d", 
	       epsilon0, sum0->GetName(), epsilon, sum->GetName(), 
	       epsilon0, n0, epsilon, n, int(ntotal)));
#if 0
  for (Int_t i = 1; i <= ret->GetNbinsX(); i++) { 
    Double_t nc  = sum->GetBinContent(i, 0);
    Double_t nc0 = sum0->GetBinContent(i, 0);
    ret->SetBinContent(i, 0, nc + nc0); // Just count events 
  }
#endif
 
  return ret;
}

//====================================================================
AliBasedNdetaTask::CentralityBin::CentralityBin()
  : TNamed("", ""), 
    fSums(0), 
    fOutput(0),
    fSum(0), 
    fSumMC(0), 
    fTriggers(0), 
    fLow(0), 
    fHigh(0),
    fDoFinalMCCorrection(false)
{
  // 
  // Constructor 
  //
}
//____________________________________________________________________
AliBasedNdetaTask::CentralityBin::CentralityBin(const char* name, 
						Short_t low, Short_t high)
  : TNamed(name, ""), 
    fSums(0), 
    fOutput(0),
    fSum(0), 
    fSumMC(0), 
    fTriggers(0),
    fLow(low), 
    fHigh(high),
    fDoFinalMCCorrection(false)
{
  // 
  // Constructor 
  // 
  // Parameters:
  //    name Name used for histograms (e.g., Forward)
  //    low  Lower centrality cut in percent 
  //    high Upper centrality cut in percent 
  //
  if (low <= 0 && high <= 0) { 
    fLow  = 0; 
    fHigh = 0;
    SetTitle("All centralities");
  }
  else {
    fLow  = low;
    fHigh = high;
    SetTitle(Form("Centrality bin from %3d%% to %3d%%", low, high));
  }
}
//____________________________________________________________________
AliBasedNdetaTask::CentralityBin::CentralityBin(const CentralityBin& o)
  : TNamed(o), 
    fSums(o.fSums), 
    fOutput(o.fOutput),
    fSum(o.fSum), 
    fSumMC(o.fSumMC), 
    fTriggers(o.fTriggers), 
    fLow(o.fLow), 
    fHigh(o.fHigh),
    fDoFinalMCCorrection(o.fDoFinalMCCorrection)
{
  // 
  // Copy constructor 
  // 
  // Parameters:
  //    other Object to copy from 
  //
}
//____________________________________________________________________
AliBasedNdetaTask::CentralityBin::~CentralityBin()
{
  // 
  // Destructor 
  //
  if (fSums) fSums->Delete();
  if (fOutput) fOutput->Delete();
}

//____________________________________________________________________
AliBasedNdetaTask::CentralityBin&
AliBasedNdetaTask::CentralityBin::operator=(const CentralityBin& o)
{
  // 
  // Assignment operator 
  // 
  // Parameters:
  //    other Object to assign from 
  // 
  // Return:
  //    Reference to this 
  //
  if (&o == this) return *this; 
  SetName(o.GetName());
  SetTitle(o.GetTitle());
  fSums      = o.fSums;
  fOutput    = o.fOutput;
  fSum       = o.fSum;
  fSumMC     = o.fSumMC;
  fTriggers  = o.fTriggers;
  fLow       = o.fLow;
  fHigh      = o.fHigh;
  fDoFinalMCCorrection = o.fDoFinalMCCorrection;

  return *this;
}
//____________________________________________________________________
Int_t
AliBasedNdetaTask::CentralityBin::GetColor(Int_t fallback) const 
{
  if (IsAllBin()) return fallback;
  Float_t  fc       = (fLow+double(fHigh-fLow)/2) / 100;
  Int_t    nCol     = gStyle->GetNumberOfColors();
  Int_t    icol     = TMath::Min(nCol-1,int(fc * nCol + .5));
  Int_t    col      = gStyle->GetColorPalette(icol);
  return col;
}
//____________________________________________________________________
const char* 
AliBasedNdetaTask::CentralityBin::GetListName() const
{
  // 
  // Get the list name 
  // 
  // Return:
  //    List Name 
  //
  if (IsAllBin()) return "all"; 
  return Form("cent%03d_%03d", fLow, fHigh);
}
//____________________________________________________________________
void
AliBasedNdetaTask::CentralityBin::CreateOutputObjects(TList* dir)
{
  // 
  // Create output objects 
  // 
  // Parameters:
  //    dir   Parent list
  //
  fSums = new TList;
  fSums->SetName(GetListName());
  fSums->SetOwner();
  dir->Add(fSums);

  fTriggers = AliAODForwardMult::MakeTriggerHistogram("triggers");
  fTriggers->SetDirectory(0);
  fSums->Add(fTriggers);
}
//____________________________________________________________________
void
AliBasedNdetaTask::CentralityBin::CreateSums(const TH2D* data, const TH2D* mc)
{
  // 
  // Create sum histogram 
  // 
  // Parameters:
  //    data  Data histogram to clone 
  //    mc    (optional) MC histogram to clone 
  //
  if (data) {
    fSum = new Sum(GetName(),"");
    fSum->Init(fSums, data, GetColor());
  }

  // If no MC data is given, then do not create MC sum histogram 
  if (!mc) return;

  fSumMC = new Sum(GetName(), "MC");
  fSumMC->Init(fSums, mc, GetColor());
}

//____________________________________________________________________
Bool_t
AliBasedNdetaTask::CentralityBin::CheckEvent(const AliAODForwardMult* forward,
					     Int_t triggerMask,
					     Double_t vzMin, Double_t vzMax)
{
  // 
  // Check the trigger, vertex, and centrality
  // 
  // Parameters:
  //    aod Event input 
  // 
  // Return:
  //    true if the event is to be used 
  //
  if (!forward) return false;

  // We do not check for centrality here - it's already done 
  return forward->CheckEvent(triggerMask, vzMin, vzMax, 0, 0, fTriggers);
}
  
  
//____________________________________________________________________
void
AliBasedNdetaTask::CentralityBin::ProcessEvent(const AliAODForwardMult* forward,
					       Int_t triggerMask, Bool_t isZero,
					       Double_t vzMin, Double_t vzMax,
					       const TH2D* data, const TH2D* mc)
{
  // 
  // Process an event
  // 
  // Parameters:
  //    forward     Forward data (for trigger, vertex, & centrality)
  //    triggerMask Trigger mask 
  //    vzMin       Minimum IP z coordinate
  //    vzMax       Maximum IP z coordinate
  //    data        Data histogram 
  //    mc          MC histogram
  //
  if (!CheckEvent(forward, triggerMask, vzMin, vzMax)) return;
  if (!data) return;
  if (!fSum) CreateSums(data, mc);

  fSum->Add(data, isZero);
  if (mc) fSumMC->Add(mc, isZero);
}

//________________________________________________________________________
Double_t 
AliBasedNdetaTask::CentralityBin::Normalization(const TH1I& t,
						UShort_t    scheme,
						Double_t    trigEff,
						Double_t&   ntotal) const
{
  // 
  // Calculate normalization 
  // 
  // Parameters: 
  //    t       Trigger histogram
  //    scheme  Normaliztion scheme 
  //    trigEff From MC
  //    ntotal  On return, contains the number of events. 
  //
  Double_t nAll        = t.GetBinContent(AliAODForwardMult::kBinAll);
  Double_t nB          = t.GetBinContent(AliAODForwardMult::kBinB);
  Double_t nA          = t.GetBinContent(AliAODForwardMult::kBinA);
  Double_t nC          = t.GetBinContent(AliAODForwardMult::kBinC);
  Double_t nE          = t.GetBinContent(AliAODForwardMult::kBinE);
  Double_t nOffline    = t.GetBinContent(AliAODForwardMult::kBinOffline);
  Double_t nTriggered  = t.GetBinContent(AliAODForwardMult::kWithTrigger);
  Double_t nWithVertex = t.GetBinContent(AliAODForwardMult::kWithVertex);
  Double_t nAccepted   = ntotal; // t.GetBinContent(AliAODForwardMult::kAccepted);
  ntotal               = 0;
  
  if (nTriggered <= 0.1) { 
    AliError("Number of triggered events <= 0");
    return -1;
  }
  if (nWithVertex <= 0.1) { 
    AliError("Number of events with vertex <= 0");
    return -1;
  }
  ntotal          = nAccepted;
  Double_t vtxEff = nWithVertex / nTriggered;
  Double_t scaler = 1;
  Double_t beta   = nA + nC - 2*nE;

  if (scheme & kEventLevel && !(scheme & kZeroBin)) {
    ntotal = nAccepted / vtxEff;
    scaler = vtxEff;
    AliInfo(Form("Calculating event normalisation as\n"
		 " N = N_A * N_T / N_V = %d * %d / %d = %f (%f)",
		 Int_t(nAccepted), Int_t(nTriggered), Int_t(nWithVertex), 
		 ntotal, scaler));
	    
    if (scheme & kBackground) {
      //          1            E_V             E_V
      //   s = --------- = ------------- = ------------ 
      //        1 - beta   1 - beta E_V    1 - beta N_V 
      //       ---  ----       --------        ---- ---
      //       E_V  N_V          N_V           N_V  N_T 
      // 
      //          E_V
      //     = ------------
      //        1 - beta 
      //            ----
      //             N_T 
      // 
      ntotal -= nAccepted * beta / nWithVertex;
      // This one is direct and correct. 
      // scaler = 1. / (1. / vtxEff - beta / nWithVertex);
      // A simpler expresion
      scaler /= (1 - beta / nTriggered); // 0.831631 -> 0.780689
      AliInfo(Form("Calculating event normalisation as\n"
		   " beta = N_a + N_c + 2 N_e = %d + %d - 2 * %d = %d\n"
		   " N = N - N_A * beta / N_V = %f - %d * %d / %d = %f (%f)",
		   Int_t(nA), Int_t(nC), Int_t(nE), Int_t(beta),
		   nAccepted / vtxEff, Int_t(nAccepted), Int_t(beta), 
		   Int_t(nWithVertex), ntotal, scaler));
    }
  }
  if (scheme & kZeroBin) {
    // Calculate as 
    // 
    //  N = N_A + 1/E_X * N_A / N_V (N_T - N_V - beta)
    //    = N_A (1 + 1/E_X (N_T/N_V - 1 - beta / N_V))
    // 
    //  s = N_A/N = 1 / (1 + 1/E_X (N_T/N_V - 1 - beta / N_V))
    //    = N_V / (N_V + 1/E_X (N_T - N_V - beta)) 
    // 
    if (!(scheme & kBackground)) beta = 0;
    ntotal = nAccepted * (1 + 1/trigEff * (nTriggered / nWithVertex - 1 
					 - beta / nWithVertex));
    scaler = nWithVertex / (nWithVertex + 
			    1/trigEff * (nTriggered-nWithVertex-beta));
    AliInfo(Form("Calculating event normalisation as\n"
		 "  beta = N_a + N_c + 2 N_e = %d + %d - 2 * %d = %d\n"
		 "  N = N_A (1 + 1/E_X (N_T/N_V - 1 - beta / N_V)) = "
		 "%d (1 + 1 / %f (%d / %d - 1 - %d / %d)) = %f (%f)",
		 Int_t(nA), Int_t(nC), Int_t(nE), Int_t(beta),
		 Int_t(nAccepted), trigEff, Int_t(nTriggered), 
		 Int_t(nWithVertex), Int_t(beta), Int_t(nWithVertex), 
		 ntotal, scaler));
  }
  if (scheme & kTriggerEfficiency && !(scheme & kZeroBin)) {
    ntotal /= trigEff;
    scaler *= trigEff;
    AliInfo(Form("Correcting for trigger efficiency:\n"
		 " N = 1 / E_X * N = 1 / %f * %d = %f (%f)", 
		 trigEff, Int_t(ntotal), ntotal / trigEff, scaler));
  }

  AliInfo(Form("\n"
	       " Total of        %9d events for %s\n"
	       "  of these       %9d have an offline trigger\n"
	       "  of these N_T = %9d has the selected trigger\n" 
	       "  of these N_V = %9d has a vertex\n" 
	       "  of these N_A = %9d were in the selected range\n"
	       "  Triggers by hardware type:\n"
	       "    N_b   =  %9d\n"
	       "    N_ac  =  %9d (%d+%d)\n"
	       "    N_e   =  %9d\n"
	       "  Vertex efficiency:          %f\n"
	       "  Trigger efficiency:         %f\n"
	       "  Total number of events: N = %f\n"
	       "  Scaler (N_A/N):             %f",
	       Int_t(nAll), GetTitle(), Int_t(nOffline), 
	       Int_t(nTriggered), Int_t(nWithVertex), Int_t(nAccepted),
	       Int_t(nB), Int_t(nA+nC), Int_t(nA), Int_t(nC), Int_t(nE), 
	       vtxEff, trigEff, ntotal, scaler));
  return scaler;
}

//________________________________________________________________________
const char* 
AliBasedNdetaTask::CentralityBin::GetResultName(Int_t rebin,
						Bool_t sym, 
						const char* postfix) const
{
  static TString n;
  n = Form("dndeta%s%s",GetName(), postfix);
  if (rebin > 1) n.Append(Form("_rebin%02d", rebin));
  if (sym)       n.Append("_mirror");
  return n.Data();
}
//________________________________________________________________________
TH1* 
AliBasedNdetaTask::CentralityBin::GetResult(Int_t rebin,
					    Bool_t sym, 
					    const char* postfix) const
{
  if (!fOutput) { 
    AliWarning(Form("No output list defined in %s [%3d,%3d]", GetName(), 
		    fLow, fHigh));
    return 0;
  }
  TString  n = GetResultName(rebin, sym, postfix);
  TObject* o = fOutput->FindObject(n.Data());
  if (!o) { 
    // AliWarning(Form("Object %s not found in output list", n.Data()));
    return 0;
  }
  return static_cast<TH1*>(o);
}

//________________________________________________________________________
void 
AliBasedNdetaTask::CentralityBin::MakeResult(const TH2D* sum,  
					     const char* postfix, 
					     bool        rootProj, 
					     bool        corrEmpty,
					     const TH1*  shapeCorr,
					     Double_t    scaler,
					     bool        symmetrice, 
					     Int_t       rebin, 
					     bool        cutEdges, 
					     Int_t       marker,
					     Int_t       color,
					     TList*      mclist, 
					     TList*      truthlist)
{
  // 
  // Generate the dN/deta result from input 
  // 
  // Parameters: 
  //     sum        Sum of 2D hists 
  //     postfix    Post fix on names
  //     rootProj   Whether to use ROOT TH2::ProjectionX
  //     corrEmpty  Correct for empty bins 
  //     shapeCorr  Shape correction to use 
  //     scaler     Event-level normalization scaler  
  //     symmetrice Whether to make symmetric extensions 
  //     rebin      Whether to rebin
  //     cutEdges   Whether to cut edges when rebinning 
  //
  TH2D* copy    = static_cast<TH2D*>(sum->Clone(Form("d2Ndetadphi%s%s", 
						     GetName(), postfix)));
  TH1D* accNorm = ProjectX(sum, Form("norm%s%s",GetName(), postfix), 0, 0, 
			   rootProj, corrEmpty, false);
  accNorm->SetDirectory(0);

  // ---- Scale by shape correction ----------------------------------
  if (shapeCorr) copy->Divide(shapeCorr);
  else AliInfo("No shape correction specified, or disabled");
  
  // --- Normalize to the coverage -----------------------------------
  ScaleToCoverage(copy, accNorm);

  // --- Event-level normalization -----------------------------------
  copy->Scale(scaler);

  // --- Project on X axis -------------------------------------------
  TH1D* dndeta = ProjectX(copy, Form("dndeta%s%s",GetName(), postfix), 
			  1, sum->GetNbinsY(), rootProj, corrEmpty);
  dndeta->SetDirectory(0);
  // Event-level normalization 
  dndeta->Scale(1., "width");
  copy->Scale(1., "width");
  
  TH1D* dndetaMCCorrection = 0;
  TList* centlist          = 0;
  TH1D* dndetaMCtruth      = 0;
  TList* truthcentlist     = 0;
  
  // Possible final correction to <MC analysis> / <MC truth>
  if(mclist) 
    centlist = static_cast<TList*> (mclist->FindObject(GetListName()));
  if(centlist)
    dndetaMCCorrection = static_cast<TH1D*>(centlist->FindObject(Form("dndeta%s%s",GetName(), postfix)));
  if(truthlist) 
    truthcentlist = static_cast<TList*> (truthlist->FindObject(GetListName()));
  if(truthcentlist)
    dndetaMCtruth =  static_cast<TH1D*> (truthcentlist->FindObject("dndetaTruth"));
  //std::cout<<dndetaMCCorrection<<"  "<<dndetaMCtruth<<std::endl;
  if(dndetaMCCorrection && dndetaMCtruth) {
    AliInfo("Correcting with final MC correction");
    dndetaMCCorrection->Divide(dndetaMCtruth);
    dndeta->Divide(dndetaMCCorrection);
    
    //std::cout<<"histo "<<Form("dndeta%s%s",GetName(), postfix)<<"  "<<GetListName()<<"  "<<dndetaMCCorrection<<std::endl;
    //std::cout<<"truth "<<GetListName()<<"  "<<dndetaMCtruth<<std::endl;
  
  }
  else AliInfo("No final MC correction applied");
  
  // --- Set some histogram attributes -------------------------------
  TString post;
  Int_t rColor = GetColor(color);
  if (postfix && postfix[0] != '\0') post = Form(" (%s)", postfix);
  SetHistogramAttributes(dndeta,  rColor, marker, 
			 Form("ALICE %s%s", GetName(), post.Data()));
  SetHistogramAttributes(accNorm, rColor, marker, 
			 Form("ALICE %s normalisation%s", 
			      GetName(), post.Data())); 

  // --- Make symmetric extensions and rebinnings --------------------
  if (symmetrice)   fOutput->Add(Symmetrice(dndeta));
  fOutput->Add(dndeta);
  fOutput->Add(accNorm);
  fOutput->Add(copy);
  fOutput->Add(Rebin(dndeta, rebin, cutEdges));
  if (symmetrice)   fOutput->Add(Symmetrice(Rebin(dndeta, rebin, cutEdges)));
}  

//________________________________________________________________________
void 
AliBasedNdetaTask::CentralityBin::End(TList*      sums, 
				      TList*      results,
				      UShort_t    scheme,
				      const TH1*  shapeCorr, 
				      Double_t    trigEff,
				      Bool_t      symmetrice,
				      Int_t       rebin, 
				      Bool_t      rootProj,
				      Bool_t      corrEmpty, 
				      Bool_t      cutEdges, 
				      Int_t       triggerMask,
				      Int_t       marker,
				      Int_t       color, 
				      TList*      mclist,
				      TList*      truthlist) 
{
  // 
  // End of processing 
  // 
  // Parameters:
  //    sums        List of sums
  //    results     Output list of results
  //    shapeCorr   Shape correction or nil
  //    trigEff     Trigger efficiency 
  //    symmetrice  Whether to symmetrice the results
  //    rebin       Whether to rebin the results
  //    corrEmpty   Whether to correct for empty bins
  //    cutEdges    Whether to cut edges when rebinning
  //    triggerMask Trigger mask 
  //

  fSums = dynamic_cast<TList*>(sums->FindObject(GetListName()));
  if(!fSums) {
    AliError("Could not retrieve TList fSums"); 
    return; 
  }
  
  fOutput = new TList;
  fOutput->SetName(GetListName());
  fOutput->SetOwner();
  results->Add(fOutput);

  if (!fSum) { 
    AliInfo("This task did not produce any output");
    return;
  }

  fTriggers = static_cast<TH1I*>(fSums->FindObject("triggers"));

  if (!fTriggers) { 
    AliError("Couldn't find histogram 'triggers' in list");
    return;
  }
  if (!fSum) { 
    AliError(Form("No sum object for %s", GetName()));
    return;
  }

  // --- Get normalization scaler ------------------------------------
  Double_t epsilonT  = trigEff;
  Double_t epsilonT0 = trigEff;
  // TEMPORARY FIX
  if (triggerMask == AliAODForwardMult::kNSD) {
    // This is a local change 
    epsilonT = 0.96; 
    AliWarning(Form("Using hard-coded NSD trigger efficiency of %f",epsilonT));
  }
  else if (triggerMask == AliAODForwardMult::kInel) {
    // This is a local change 
    epsilonT = 0.934; 
    AliWarning(Form("Using hard-coded Inel trigger efficiency of %f",epsilonT));
  }
  if (scheme & kZeroBin) { 
    if (triggerMask==AliAODForwardMult::kInel)
      epsilonT0 = 0.785021; // 0.100240;
    else if (triggerMask==AliAODForwardMult::kInelGt0)
      epsilonT0 = 0;
    else if (triggerMask==AliAODForwardMult::kNSD)
      epsilonT0 = .706587;
    epsilonT = 1;
    AliWarning(Form("Using hard-coded NCluster>0 trigger efficiency of %f",
		    epsilonT0));
  }

  // Get our histograms 
  Double_t nSum   = 0;
  TH2D*    sum    = fSum->GetSum(fSums, fOutput, nSum, epsilonT0, 1, 
				 marker, rootProj, corrEmpty);
  Double_t nSumMC = 0;
  TH2D*    sumMC  = 0;
  if (fSumMC) sumMC = fSumMC->GetSum(fSums, fOutput, nSumMC, 
				     epsilonT0, 1, marker,
				     rootProj, corrEmpty);
  if (!sum) { 
    AliError("Failed to get sum from summer - bailing out");
    return;
  }
    
  Double_t ntotal = nSum;
  Double_t scaler = Normalization(*fTriggers, scheme, epsilonT, ntotal);
  if (scaler < 0) { 
    AliError("Failed to calculate normalization - bailing out");
    return;
  }
  fOutput->Add(fTriggers->Clone());

  // --- Make result and store ---------------------------------------
  MakeResult(sum, "", rootProj, corrEmpty, (scheme & kShape) ? shapeCorr : 0,
	     scaler, symmetrice, rebin, cutEdges, marker, color, 
	     mclist, truthlist);

  // --- Process result from TrackRefs -------------------------------
  if (sumMC) 
    MakeResult(sumMC, "MC", rootProj, corrEmpty, 
	       (scheme & kShape) ? shapeCorr : 0,
	       scaler, symmetrice, rebin, cutEdges, 
	       GetMarkerStyle(GetMarkerBits(marker)+4), color, 
	       mclist, truthlist);
  
  // Temporary stuff 
  // if (!IsAllBin()) return;

}

//
// EOF
//
