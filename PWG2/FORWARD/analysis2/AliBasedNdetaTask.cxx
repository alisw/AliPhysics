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
    fTriggerEff(1),
    fShapeCorr(0),
    fListOfCentralities(0),
    fSNNString(0),
    fSysString(0),
    fCent(0),
    fCentAxis(0),
    fNormalizationScheme(kFull)
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
    fTriggerEff(1),
    fShapeCorr(0),
    fListOfCentralities(0),
    fSNNString(0),
    fSysString(0),
    fCent(0),
    fCentAxis(0),
    fNormalizationScheme(kFull)
{
  // 
  // Constructor
  // 
  fListOfCentralities = new TObjArray(1);

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
    fSymmetrice(true),
    fCorrEmpty(true), 
    fTriggerEff(o.fTriggerEff),
    fShapeCorr(o.fShapeCorr),
    fListOfCentralities(o.fListOfCentralities),
    fSNNString(o.fSNNString),
    fSysString(o.fSysString),
    fCent(o.fCent),
    fCentAxis(o.fCentAxis),
    fNormalizationScheme(o.fNormalizationScheme)
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
  AliInfo(Form("Adding centrality bin %p [%3d,%3d] @ %d", bin, low, high, at));
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
    if      (s.CompareTo("EVENT")     == 0) scheme |= kEventLevel;
    else if (s.CompareTo("ALTEVENT")  == 0) scheme |= kAltEventLevel;
    else if (s.CompareTo("SHAPE")     == 0) scheme |= kShape;
    else if (s.CompareTo("BACKGROUND")== 0) scheme |= kBackground;
    else if (s.CompareTo("FULL")      == 0) scheme |= kFull;
    else if (s.CompareTo("ALTFULL")   == 0) scheme |= kAltFull;
    else if (s.CompareTo("NONE")      == 0) scheme |= kNone;
    else 
      Warning("SetNormalizationScheme", "Unknown option %s", s.Data());
  }
  SetNormalizationScheme(scheme);
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
  UShort_t    trgMask = 0;
  TString     trgs(mask);
  trgs.ToUpper();
  TObjString* trg;
  TIter       next(trgs.Tokenize(" ,|"));
  while ((trg = static_cast<TObjString*>(next()))) { 
    TString s(trg->GetString());
    if      (s.IsNull()) continue;
    if      (s.CompareTo("INEL")  == 0) trgMask = AliAODForwardMult::kInel;
    else if (s.CompareTo("INEL>0")== 0) trgMask = AliAODForwardMult::kInelGt0;
    else if (s.CompareTo("NSD")   == 0) trgMask = AliAODForwardMult::kNSD;
    else 
      Warning("SetTriggerMask", "Unknown trigger %s", s.Data());
  }
  if (trgMask == 0) trgMask = 1;
  SetTriggerMask(trgMask);
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
  fSums->Add(fCentAxis);


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

  // Loop over centrality bins 
  CentralityBin* allBin = 
    static_cast<CentralityBin*>(fListOfCentralities->At(0));
  allBin->ProcessEvent(forward, fTriggerMask, fVtxMin, fVtxMax, data, dataMC);
  
  // Find this centrality bin 
  if (fCentAxis && fCentAxis->GetNbins() > 0) {
    Int_t          icent   = fCentAxis->FindBin(cent);
    CentralityBin* thisBin = 0;
    if (icent >= 1 && icent <= fCentAxis->GetNbins()) 
      thisBin = static_cast<CentralityBin*>(fListOfCentralities->At(icent));
    if (thisBin)
      thisBin->ProcessEvent(forward, fTriggerMask, fVtxMin, fVtxMax, 
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
  h->SetMarkerSize(1);
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
TH1D*
AliBasedNdetaTask::ProjectX(const TH2D* h, 
			    const char* name,
			    Int_t firstbin, 
			    Int_t lastbin, 
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
#if USE_ROOT_PROJECT
  return h->ProjectionX(name, firstbin, lastbin, (error ? "e" : ""));
#endif
  
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

  fSNNString = static_cast<TNamed*>(fSums->FindObject("sNN"));
  fSysString = static_cast<TNamed*>(fSums->FindObject("sys"));
  fCentAxis  = static_cast<TAxis*>(fSums->FindObject("centAxis"));

  if(fSysString && fSNNString && 
     fSysString->GetUniqueID() == AliForwardUtil::kPP)
    LoadNormalizationData(fSysString->GetUniqueID(),
			  fSNNString->GetUniqueID());

  // Loop over centrality bins 
  TIter next(fListOfCentralities);
  CentralityBin* bin = 0;
  while ((bin = static_cast<CentralityBin*>(next()))) 
    bin->End(fSums, fOutput, fNormalizationScheme, fShapeCorr, fTriggerEff,
	     fSymmetrice, fRebin, fCorrEmpty, fCutEdges, 
	     fVtxMin, fVtxMax, fTriggerMask, GetColor(), GetMarker());

  // Output collision energy string 
  if (fSNNString) fOutput->Add(fSNNString->Clone());

  // Output collision system string 
  if (fSysString) fOutput->Add(fSysString->Clone());

  // Output centrality axis 
  if (fCentAxis) fOutput->Add(fCentAxis);

  // Output trigger string 
  TNamed* trigString = 
    new TNamed("trigString", AliAODForwardMult::GetTriggerString(fTriggerMask));
  trigString->SetUniqueID(fTriggerMask);
  fOutput->Add(trigString);

  // Output vertex axis 
  TAxis* vtxAxis = new TAxis(1,fVtxMin,fVtxMax);
  vtxAxis->SetName("vtxAxis");
  vtxAxis->SetTitle(Form("v_{z}#in[%+5.1f,%+5.1f]cm", fVtxMin,fVtxMax));
  fOutput->Add(vtxAxis);

  // Output trigger efficiency and shape correction 
  fOutput->Add(new TParameter<Double_t>("triggerEff", fTriggerEff)); 
  if (fShapeCorr) fOutput->Add(fShapeCorr);

  TString tscheme;
  if (fNormalizationScheme & kEventLevel)    tscheme.Append("EVENTLEVEL,");
  if (fNormalizationScheme & kAltEventLevel) tscheme.Append("ALTEVENTLEVEL,");
  if (fNormalizationScheme & kBackground)    tscheme.Append("BACKGROUND,");
  if (fNormalizationScheme & kShape)         tscheme.Append("SHAPE,");
  
  TNamed* scheme = new TNamed("scheme", tscheme.Data()); 
  scheme->SetUniqueID(fNormalizationScheme); 
  fOutput->Add(scheme);

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
  
  if(fShapeCorr && (fTriggerEff != 1)) {
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
  TString shapeCorName(Form("h%sNormalization", 
			    fTriggerMask == AliAODForwardMult::kInel ? "Inel" :
			    fTriggerMask == AliAODForwardMult::kNSD ? "NSD" :
			    fTriggerMask == AliAODForwardMult::kInelGt0 ?
			    "InelGt0" : "All"));
  TH2F* shapeCor = dynamic_cast<TH2F*>(fin->Get(shapeCorName));
  if (shapeCor) SetShapeCorrection(shapeCor);

  // Trigger efficiency
  TString effName(Form("%sTriggerEff", 
		       fTriggerMask == AliAODForwardMult::kInel ? "inel" :
		       fTriggerMask == AliAODForwardMult::kNSD ? "nsd" :
		       fTriggerMask == AliAODForwardMult::kInelGt0 ?
		       "inelgt0" : "all"));
  TParameter<float>* eff = static_cast<TParameter<float>*>(fin->Get(effName));
  if (eff) SetTriggerEff(eff->GetVal());

  // Print - out
  if (shapeCor && eff) AliInfo("Loaded objects for normalization.");
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
  s->SetMarkerStyle(h->GetMarkerStyle()+4);
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

//====================================================================
AliBasedNdetaTask::CentralityBin::CentralityBin()
  : TNamed("", ""), 
    fSums(0), 
    fOutput(0),
    fSum(0), 
    fSumMC(0), 
    fTriggers(0), 
    fLow(0), 
    fHigh(0)
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
    fHigh(high)
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
    fHigh(o.fHigh)
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
  SetName(o.GetName());
  SetTitle(o.GetTitle());
  fSums      = o.fSums;
  fOutput    = o.fOutput;
  fSum       = o.fSum;
  fSumMC     = o.fSumMC;
  fTriggers  = o.fTriggers;
  fLow       = o.fLow;
  fHigh      = o.fHigh;

  return *this;
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
  fSum = static_cast<TH2D*>(data->Clone(GetName()));
  fSum->SetDirectory(0);
  fSum->Reset();
  fSums->Add(fSum);
  
  // If no MC data is given, then do not create MC sum histogram 
  if (!mc) return;

  fSumMC = static_cast<TH2D*>(mc->Clone(Form("%sMC", GetName())));
  fSumMC->SetDirectory(0);
  fSumMC->Reset();
  fSums->Add(fSumMC);
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
					       Int_t triggerMask,
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
  
  fSum->Add(data);
  if (mc) fSumMC->Add(mc);
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
				      Bool_t      corrEmpty, 
				      Bool_t      cutEdges, 
				      Double_t    vzMin, 
				      Double_t    vzMax, 
				      Int_t       triggerMask,
				      Int_t       color, 
				      Int_t       marker) 
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
  //    vzMin       Minimum IP z coordinate
  //    vzMax 	  Maximum IP z coordinate
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

  fSum      = static_cast<TH2D*>(fSums->FindObject(GetName()));
  fSumMC    = static_cast<TH2D*>(fSums->FindObject(Form("%sMC", GetName())));
  fTriggers = static_cast<TH1I*>(fSums->FindObject("triggers"));

  if (!fTriggers) { 
    AliError("Couldn't find histogram 'triggers' in list");
    return;
  }
  if (!fSum) { 
    AliError(Form("Couldn't find histogram '%s' in list", GetName()));
    return;
  }

  TH1I& t = *fTriggers;
  Int_t nAll        = Int_t(t.GetBinContent(AliAODForwardMult::kBinAll));
  Int_t nB          = Int_t(t.GetBinContent(AliAODForwardMult::kBinB));
  Int_t nA          = Int_t(t.GetBinContent(AliAODForwardMult::kBinA));
  Int_t nC          = Int_t(t.GetBinContent(AliAODForwardMult::kBinC));
  Int_t nE          = Int_t(t.GetBinContent(AliAODForwardMult::kBinE));
  Int_t nMB         = Int_t(t.GetBinContent(AliAODForwardMult::kBinInel));
  Int_t nTriggered  = Int_t(t.GetBinContent(AliAODForwardMult::kWithTrigger));
  Int_t nWithVertex = Int_t(t.GetBinContent(AliAODForwardMult::kWithVertex));
  Int_t nAccepted   = Int_t(t.GetBinContent(AliAODForwardMult::kAccepted));
  Int_t nGood       = nB - nA - nC + 2 * nE;

  
  if (nTriggered <= 0.1) { 
    AliError("Number of triggered events <= 0");
    return;
  }
  if (nWithVertex <= 0.1) { 
    AliError("Number of events with vertex <= 0");
    return;
  }
  if (nGood <= 0.1) { 
    AliWarning(Form("Number of good events=%d=%d-%d-%d+2*%d<=0",
		    nGood, nB, nA, nC, nE));
    nGood = nMB;
  }


  // Scaling 
  //  N_A + N_A/N_V (N_T - N_V) = N_A + N_A/N_V*N_T - N_A 
  //                            = N_A/N_V * N_T 
  // where 
  //    N_A = nAccepted 
  //    N_V = nWithVertex 
  //    N_T = nTriggered 
  // so that the normalisation is simply 
  //    N_A / E_V 
  // where 
  //    E_V=N_V/N_T is the vertex efficiency from data 
  // 
  // If we calculate the number of events with no vertex as 
  // (N_T - N_A), we get 
  //  
  //  N_A + N_A/N_V (N_T - N_A) = N_A (1 + 1/N_V (N_T - N_A))
  //                            = N_A (1 + N_T/N_V - N_A/N_V)
  //                            = N_A (1 + 1/E_V - N_A/N_V)
  // 
  // Double_t alpha    = Double_t(nAccepted) / nWithVertex;
  // Double_t vNorm    = nAccepted + alpha*(nTriggered - nWithVertex);
  Double_t ntotal = 1;
  Double_t vtxEff = double(nWithVertex) / double(nTriggered);;
  Bool_t   shape  = (scheme & kShape);
  if     (scheme & kEventLevel) {
    AliInfo(Form("Calculating event normalisation as "
		 "1 / trigEff * nAccepted * nTriggered / nWithVertex "
		 "= 1/%f * %d * %d / %d = %f",
		 trigEff, nAccepted, nTriggered, nWithVertex,
		 1. / trigEff * nAccepted * double(nTriggered) / nWithVertex));
    ntotal = vtxEff * trigEff;
  }
  else if (scheme & kAltEventLevel) {
    AliInfo(Form("Calculating event normalisation as " 
		 "1 / trigEff * nAccepted * (1 + nTriggered / nWithVertex - "
		 "nAccepted / nWithVertex) = "
		 "1/ %f * %d (1 + %d / %d - %d / %d) = %f", 
		 trigEff, nAccepted, nTriggered, nWithVertex, 
		 nAccepted, nWithVertex, 
		 1/trigEff * nAccepted * (1 + double(nTriggered)/nWithVertex - 
					double(nAccepted) / nWithVertex)));
    ntotal = trigEff * 1 / (1 + double(nTriggered)/nWithVertex - 
			    double(nAccepted) / nWithVertex);
  }

  AliInfo(Form("Total of %9d events for %s\n"
	       "                   of these %9d are triggered minimum bias\n"
	       "                   of these %9d has a %s trigger\n" 
	       "                   of these %9d has a vertex\n" 
	       "                   of these %9d were in [%+4.1f,%+4.1f]cm\n"
	       "                   Triggers by type:\n"
	       "                     B   = %9d\n"
	       "                     A|C = %9d (%9d+%-9d)\n"
	       "                     E   = %9d\n"
	       "                   Implies %9d good triggers\n"
	       "                   Vertex efficiency:      %f\n"
	       "                   Trigger efficiency:     %f\n"
	       "                   Total number of events: %f\n",
	       nAll, GetTitle(), nMB, nTriggered, 
	       AliAODForwardMult::GetTriggerString(triggerMask),
	       nWithVertex, nAccepted,
	       vzMin, vzMax, 
	       nB, nA+nC, nA, nC, nE, nGood, vtxEff, trigEff, ntotal));
  
  const char* name = GetName();
  
  // Get acceptance normalisation from underflow bins 
  TH1D* norm = ProjectX(fSum, Form("norm%s",name), 0, 0, corrEmpty, false);
  norm->SetDirectory(0);

  // Scale by shape correction 
  if(shape && shapeCorr) fSum->Divide(shapeCorr);
  else AliInfo("No shape correction specified, or disabled");
  
  // Project on X axis 
  TH1D* dndeta = ProjectX(fSum, Form("dndeta%s",name), 1, fSum->GetNbinsY(),
			  corrEmpty);
  dndeta->SetDirectory(0);

  // Normalize to the acceptance -
  dndeta->Divide(norm);
  
  // Scale by the vertex efficiency 
  dndeta->Scale(ntotal, "width");
  
  SetHistogramAttributes(dndeta, color, marker, Form("ALICE %s", name));
  SetHistogramAttributes(norm,   color, marker, Form("ALICE %s normalisation", 
						     name)); 

  fOutput->Add(fTriggers->Clone());
  if (symmetrice)   fOutput->Add(Symmetrice(dndeta));
  fOutput->Add(dndeta);
  fOutput->Add(norm);
  fOutput->Add(Rebin(dndeta, rebin, cutEdges));
  if (symmetrice)   fOutput->Add(Symmetrice(Rebin(dndeta, rebin, cutEdges)));

  if (fSumMC) { 
    // Get acceptance normalisation from underflow bins 
    TH1D* normMC   = ProjectX(fSumMC,Form("norm%sMC", name), 0, 1, 
			      corrEmpty, false);
    // Project onto eta axis - _ignoring_underflow_bins_!
    TH1D* dndetaMC = ProjectX(fSumMC,Form("dndeta%sMC", name),1,
			      fSum->GetNbinsY(), corrEmpty);
    // Normalize to the acceptance 
    dndetaMC->Divide(normMC);
    // Scale by the vertex efficiency 
    dndetaMC->Scale(ntotal, "width");
    normMC->Scale(1. / nAccepted);

    SetHistogramAttributes(dndetaMC, color+2, marker, 
			   Form("ALICE %s (MC)",name));
    SetHistogramAttributes(normMC,   color+2, marker, 
			   Form("ALICE %s (MC) normalisation",name)); 

    fOutput->Add(dndetaMC);
    if (symmetrice)   fOutput->Add(Symmetrice(dndetaMC));    

    fOutput->Add(normMC);
    fOutput->Add(Rebin(dndetaMC, rebin, cutEdges));

    if (symmetrice)   
      fOutput->Add(Symmetrice(Rebin(dndetaMC, rebin, cutEdges)));
  }

  // Temporary stuff 
  if (!IsAllBin()) return;

  TFile* forward = TFile::Open("forward.root", "READ");
  if (!forward)  { 
    AliWarning(Form("No forward.root file found"));
    return;
  }

  TH1D* shapeCorrProj = 0;
  if (shapeCorr) {
    shapeCorrProj = static_cast<const TH2D*>(shapeCorr)->ProjectionX();
    shapeCorrProj->Scale(1. / shapeCorr->GetNbinsY());
    shapeCorrProj->SetDirectory(0);
    fOutput->Add(shapeCorrProj);
  }

  TList* official = static_cast<TList*>(forward->Get("official"));
  if (official) { 
    TH1F*  histEta  = static_cast<TH1F*>(official->FindObject("fHistEta"));
    if (histEta)  {
      TH1D* oEta = new TH1D("tracks", "", histEta->GetNbinsX(), 
			    histEta->GetXaxis()->GetXmin(), 
			    histEta->GetXaxis()->GetXmax());
      for (Int_t i = 1; i < histEta->GetNbinsX(); i++) {
	oEta->SetBinContent(i, histEta->GetBinContent(i));
	oEta->SetBinError(i, histEta->GetBinError(i));
      }
      if (shapeCorrProj) oEta->Divide(shapeCorrProj);
      oEta->Scale(ntotal/nAccepted, "width");
      oEta->SetDirectory(0);
      oEta->SetMarkerStyle(marker+4);
      oEta->SetMarkerColor(color+5);
      fOutput->Add(oEta);
      fOutput->Add(Rebin(oEta, rebin, false));
    }
    else 
      AliWarning(Form("Couldn't find histogram fHistEta in list %s", 
		      official->GetName()));
  }
  else 
    AliWarning(Form("Couldn't find list 'official' in %s",forward->GetName()));

  TList* tracks = static_cast<TList*>(forward->Get("tracks"));
  if (tracks) { 
    TH1F*  histEta  = static_cast<TH1F*>(tracks->FindObject("fHistEta"));
    if (histEta)  {
      TH1D* oEta = new TH1D("tracks", "", histEta->GetNbinsX(), 
			    histEta->GetXaxis()->GetXmin(), 
			    histEta->GetXaxis()->GetXmax());
      for (Int_t i = 1; i < histEta->GetNbinsX(); i++) {
	oEta->SetBinContent(i, histEta->GetBinContent(i));
	oEta->SetBinError(i, histEta->GetBinError(i));
      }
      if (shapeCorrProj) oEta->Divide(shapeCorrProj);
      oEta->Scale(ntotal/nAccepted, "width");
      oEta->SetDirectory(0);
      oEta->SetMarkerStyle(marker);
      oEta->SetMarkerColor(color+5);
      fOutput->Add(oEta);
      fOutput->Add(Rebin(oEta, rebin, false));
    }
    else 
      AliWarning(Form("Couldn't find histogram fHistEta in list %s", 
		      tracks->GetName()));
  }
  else 
    AliWarning(Form("Couldn't find list 'tracks' in %s",forward->GetName()));

  forward->Close();
}

//
// EOF
//
