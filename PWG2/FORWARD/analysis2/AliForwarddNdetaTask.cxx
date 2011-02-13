//====================================================================
#include "AliForwarddNdetaTask.h"
#include <TMath.h>
#include <TH2D.h>
#include <TH1D.h>
#include <THStack.h>
#include <TList.h>
#include <AliAnalysisManager.h>
#include <AliAODEvent.h>
#include <AliAODHandler.h>
#include <AliAODInputHandler.h>
#include "AliForwardUtil.h"

//____________________________________________________________________
AliForwarddNdetaTask::AliForwarddNdetaTask()
  : AliAnalysisTaskSE(), 
    fSumForward(0),	//  Sum of histograms 
    fSumForwardMC(0),	//  Sum of MC histograms (if any)
    fSumPrimary(0),	//  Sum of primary histograms
    fSumCentral(0),	//  Sum of central histograms
    fCentral(0),	// Cache of central histogram
    fPrimary(0),	// Cache of primary histogram
    fSums(0),		// Container of sums 
    fOutput(0),		// Container of outputs 
    fTriggers(0),	// Histogram of triggers 
    fVtxMin(0),		// Minimum v_z
    fVtxMax(0),		// Maximum v_z
    fTriggerMask(0),    // Trigger mask 
    fRebin(0),       	// Rebinning factor 
    fCutEdges(false),
    fSNNString(0),
    fSysString(0)
{}

//____________________________________________________________________
AliForwarddNdetaTask::AliForwarddNdetaTask(const char* /* name */)
  : AliAnalysisTaskSE("Forward"), 
    fSumForward(0),	//  Sum of histograms 
    fSumForwardMC(0),	//  Sum of MC histograms (if any)
    fSumPrimary(0),	//  Sum of primary histograms
    fSumCentral(0),	//  Sum of central histograms
    fCentral(0),	// Cache of central histogram
    fPrimary(0),	// Cache of primary histogram
    fSums(0),		// Container of sums 
    fOutput(0),		// Container of outputs 
    fTriggers(0),	// Histogram of triggers 
    fVtxMin(-10),	// Minimum v_z
    fVtxMax(10),	// Maximum v_z
    fTriggerMask(AliAODForwardMult::kInel), 
    fRebin(5),		// Rebinning factor 
    fCutEdges(false),
    fSNNString(0),
    fSysString(0)
{
  // Output slot #1 writes into a TH1 container
  DefineOutput(1, TList::Class()); 
  DefineOutput(2, TList::Class()); 
}

//____________________________________________________________________
AliForwarddNdetaTask::AliForwarddNdetaTask(const AliForwarddNdetaTask& o)
  : AliAnalysisTaskSE(o), 
    fSumForward(o.fSumForward),		// TH2D* -  Sum of histograms 
    fSumForwardMC(o.fSumForwardMC),		// TH2D* -  Sum of MC histograms (if any)
    fSumPrimary(o.fSumPrimary),		// TH2D* -  Sum of primary histograms
    fSumCentral(o.fSumCentral),		// TH2D* -  Sum of central histograms
    fCentral(o.fCentral),       //! Cache of central histogram
    fPrimary(o.fPrimary),       //! Cache of primary histogram
    fSums(o.fSums),		// TList* - Container of sums 
    fOutput(o.fOutput),		// TList* - Container of outputs 
    fTriggers(o.fTriggers),		// TH1D* - Histogram of triggers 
    fVtxMin(o.fVtxMin),		// Double_t - Minimum v_z
    fVtxMax(o.fVtxMax),		// Double_t - Maximum v_z
    fTriggerMask(o.fTriggerMask),		// Int_t - Trigger mask 
    fRebin(o.fRebin),		// Int_t - Rebinning factor 
    fCutEdges(o.fCutEdges),		// Bool_t - Whether to cut edges when rebinning
    fSNNString(o.fSNNString),		// TNamed* - 
    fSysString(o.fSysString)		// TNamed* - 
{}

//____________________________________________________________________
AliForwarddNdetaTask::~AliForwarddNdetaTask()
{
  if (fSums) { 
    fSums->Delete();
    delete fSums;
    fSums = 0;
  }
  if (fOutputs) { 
    fOutputs->Delete();
    delete fOutputs;
    fOutputs = 0;
  }
  if (fCentral)    delete fCentral;
  if (fPrimary)    delete fPrimary;
}

//________________________________________________________________________
void 
AliForwarddNdetaTask::SetTriggerMask(const char* mask)
{
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
AliForwarddNdetaTask::UserCreateOutputObjects()
{
  // Create histograms
  // Called once (on the worker node)

  fOutput = new TList;
  fOutput->SetName(Form("%s_result", GetName()));
  fOutput->SetOwner();

  fSums = new TList;
  fSums->SetName(Form("%s_sums", GetName()));
  fSums->SetOwner();


  fTriggers = new TH1D("triggers", "Number of triggers", 
		       kAccepted, 1, kAccepted);
  fTriggers->SetYTitle("# of events");
  fTriggers->GetXaxis()->SetBinLabel(kAll,         "All events");
  fTriggers->GetXaxis()->SetBinLabel(kB,           "w/B trigger");
  fTriggers->GetXaxis()->SetBinLabel(kA,           "w/A trigger");
  fTriggers->GetXaxis()->SetBinLabel(kC,           "w/C trigger");
  fTriggers->GetXaxis()->SetBinLabel(kE,           "w/E trigger");
  fTriggers->GetXaxis()->SetBinLabel(kMB,          "w/Collision trigger");
  fTriggers->GetXaxis()->SetBinLabel(kWithVertex,  "w/Vertex");
  fTriggers->GetXaxis()->SetBinLabel(kWithTrigger, "w/Selected trigger");
  fTriggers->GetXaxis()->SetBinLabel(kAccepted,    "Accepted by cut");
  fTriggers->GetXaxis()->SetNdivisions(kAccepted, false);
  fTriggers->SetFillColor(kRed+1);
  fTriggers->SetFillStyle(3001);
  fTriggers->SetStats(0);
  fSums->Add(fTriggers);

  fSNNString = new TNamed("sNN", "");
  fSysString = new TNamed("sys", "");
  fSums->Add(fSNNString);
  fSums->Add(fSysString);

  // Check that we have an AOD input handler 
  AliAnalysisManager* am = AliAnalysisManager::GetAnalysisManager();
  AliAODInputHandler* ah = 
    dynamic_cast<AliAODInputHandler*>(am->GetInputEventHandler());
  if (!ah) AliFatal("No AOD input handler set in analysis manager");

  // Post data for ALL output slots >0 here, to get at least an empty histogram
  PostData(1, fSums); 
}

//____________________________________________________________________
TH2D*
AliForwarddNdetaTask::CloneHist(TH2D* in, const char* name) 
{
  if (!in) return 0;
  TH2D* ret = static_cast<TH2D*>(in->Clone(name));
  ret->SetDirectory(0);
  ret->Sumw2();
  fSums->Add(ret);

  return ret;
}

//____________________________________________________________________
void 
AliForwarddNdetaTask::UserExec(Option_t *) 
{
  // Main loop
  AliAODEvent* aod = dynamic_cast<AliAODEvent*>(InputEvent());
  if (!aod) {
    AliError("Cannot get the AOD event");
    return;
  }  

  // Get objects from the event structure 
  TObject* oForward   = aod->FindListObject("Forward");
  TObject* oForwardMC = aod->FindListObject("ForwardMC");
  TObject* oPrimary   = aod->FindListObject("primary");
  TObject* oCentral   = aod->FindListObject("Central");

  // We should have a forward object at least 
  if (!oForward) { 
    AliWarning("No Forward object found AOD");
    return;
  }

  // Cast to good types 
  AliAODForwardMult* forward   = static_cast<AliAODForwardMult*>(oForward);
  AliAODForwardMult* forwardMC = static_cast<AliAODForwardMult*>(oForwardMC);
  TH2D*              primary   = static_cast<TH2D*>(oPrimary);
  TH2D*              central   = static_cast<TH2D*>(oCentral);

  static bool first = true;
  if (first) { 
    UShort_t sNN = forward->GetSNN();
    UShort_t sys = forward->GetSystem();
    fSNNString->SetTitle(AliForwardUtil::CenterOfMassEnergyString(sNN));
    fSNNString->SetUniqueID(sNN);
    fSysString->SetTitle(AliForwardUtil::CollisionSystemString(sys));
    fSysString->SetUniqueID(sys);
    first = false;
  }

  // Create our sum histograms 
  if (!fSumForward) fSumForward = CloneHist(&forward->GetHistogram(),"forward");
  if (!fSumForwardMC && forwardMC) 
    fSumForwardMC = CloneHist(&forwardMC->GetHistogram(),"forwardMC");
  if (!fSumPrimary && primary) fSumPrimary = CloneHist(primary, "primary");
  if (!fSumCentral && central) fSumCentral = CloneHist(central, "central");

  // Add contribtion from MC 
  if (primary) fSumPrimary->Add(primary);
  
  // Count event 
  fTriggers->AddBinContent(kAll);
  if (forward->IsTriggerBits(AliAODForwardMult::kB)) 
    fTriggers->AddBinContent(kB);
  if (forward->IsTriggerBits(AliAODForwardMult::kA)) 
    fTriggers->AddBinContent(kA);
  if (forward->IsTriggerBits(AliAODForwardMult::kC)) 
    fTriggers->AddBinContent(kC);
  if (forward->IsTriggerBits(AliAODForwardMult::kE)) 
    fTriggers->AddBinContent(kE);
  if (forward->IsTriggerBits(AliAODForwardMult::kInel)) 
    fTriggers->AddBinContent(kMB);

  // Check if we have an event of interest. 
  if (!forward->IsTriggerBits(fTriggerMask)) return;
  fTriggers->AddBinContent(kWithTrigger);

  // Check that we have a valid vertex
  if (!forward->HasIpZ()) return;
  fTriggers->AddBinContent(kWithVertex);

  // Check that vertex is within cuts 
  if (!forward->InRange(fVtxMin, fVtxMax)) return;
  fTriggers->AddBinContent(kAccepted);

  // Add contribution 
  fSumForward->Add(&(forward->GetHistogram()));
  if (fSumForwardMC) fSumForwardMC->Add(&(forwardMC->GetHistogram()));
  if (fSumPrimary)   fSumPrimary->Add(primary);
  if (fSumCentral)   fSumCentral->Add(central);

  PostData(1, fSums);
}

//________________________________________________________________________
void 
AliForwarddNdetaTask::SetHistogramAttributes(TH1D* h, Int_t colour, Int_t marker,
					  const char* title, const char* ytitle)
{
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
void 
AliForwarddNdetaTask::Terminate(Option_t *) 
{
  // Draw result to screen, or perform fitting, normalizations Called
  // once at the end of the query
        
  fSums = dynamic_cast<TList*> (GetOutputData(1));
  if(!fOutput) {
    AliError("Could not retrieve TList fSums"); 
    return; 
  }
  
  if (!fOutput) { 
    fOutput = new TList;
    fOutput->SetName(Form("%s_result", GetName()));
    fOutput->SetOwner();
  }

  fSumForward     = static_cast<TH2D*>(fSums->FindObject("forward"));
  fSumForwardMC   = static_cast<TH2D*>(fSums->FindObject("forwardMC"));
  fSumPrimary     = static_cast<TH2D*>(fSums->FindObject("primary"));
  fSumCentral     = static_cast<TH2D*>(fSums->FindObject("central"));
  fTriggers       = static_cast<TH1D*>(fSums->FindObject("triggers"));
  fSNNString      = static_cast<TNamed*>(fSums->FindObject("sNN"));
  fSysString      = static_cast<TNamed*>(fSums->FindObject("sys"));

  if (!fTriggers) { 
    AliError("Couldn't find histogram 'triggers' in list");
    return;
  }
  if (!fSumForward) { 
    AliError("Couldn't find histogram 'forward' in list");
    return;
  }

  Int_t nAll        = Int_t(fTriggers->GetBinContent(kAll));
  Int_t nB          = Int_t(fTriggers->GetBinContent(kB));
  Int_t nA          = Int_t(fTriggers->GetBinContent(kA));
  Int_t nC          = Int_t(fTriggers->GetBinContent(kC));
  Int_t nE          = Int_t(fTriggers->GetBinContent(kE));
  Int_t nMB         = Int_t(fTriggers->GetBinContent(kMB));
  Int_t nTriggered  = Int_t(fTriggers->GetBinContent(kWithTrigger));
  Int_t nWithVertex = Int_t(fTriggers->GetBinContent(kWithVertex));
  Int_t nAccepted   = Int_t(fTriggers->GetBinContent(kAccepted));
  Int_t nGood       = nB - nA - nC + 2 * nE;
  Double_t vtxEff   = Double_t(nMB) / nTriggered * Double_t(nAccepted) / nGood;
  AliInfo(Form("Total of %9d events\n"
	       "                   of these %9d are minimum bias\n"
	       "                   of these %9d has a %s trigger\n" 
	       "                   of these %9d has a vertex\n" 
	       "                   of these %9d were in [%+4.1f,%+4.1f]cm\n"
	       "                   Triggers by type:\n"
	       "                     B   = %9d\n"
	       "                     A|C = %9d (%9d+%-9d)\n"
	       "                     E   = %9d\n"
	       "                   Implies %9d good triggers\n"
	       "                   Vertex efficiency: %f",
	       nAll, nMB, nTriggered, 
	       AliAODForwardMult::GetTriggerString(fTriggerMask),
	       nWithVertex, nAccepted,
	       fVtxMin, fVtxMax, 
	       nB, nA+nC, nA, nC, nE, nGood, vtxEff));
  
  // Get acceptance normalisation from underflow bins 
  TH1D* normForward   = fSumForward->ProjectionX("normForward", 0, 1, "");
  // Project onto eta axis - _ignoring_underflow_bins_!
  TH1D* dndetaForward = fSumForward->ProjectionX("dndetaForward", 1, -1, "e");
  // Normalize to the acceptance 
  dndetaForward->Divide(normForward);
  // Scale by the vertex efficiency 
  dndetaForward->Scale(vtxEff, "width");

  SetHistogramAttributes(dndetaForward, kRed+1, 20, "ALICE Forward");
  SetHistogramAttributes(normForward, kRed+1, 20, "ALICE Forward", "Normalisation");

  fOutput->Add(fTriggers->Clone());
  fOutput->Add(fSNNString->Clone());
  fOutput->Add(fSysString->Clone());
  fOutput->Add(dndetaForward);
  fOutput->Add(normForward);
  fOutput->Add(Rebin(dndetaForward));

  if (fSumForwardMC) { 
    // Get acceptance normalisation from underflow bins 
    TH1D* normForwardMC   = fSumForwardMC->ProjectionX("normForwardMC", 0, 1, "");
    // Project onto eta axis - _ignoring_underflow_bins_!
    TH1D* dndetaForwardMC = fSumForwardMC->ProjectionX("dndetaForwardMC", 1, -1, "e");
    // Normalize to the acceptance 
    dndetaForwardMC->Divide(normForwardMC);
    // Scale by the vertex efficiency 
    dndetaForwardMC->Scale(vtxEff, "width");

    SetHistogramAttributes(dndetaForwardMC, kRed+3, 21, "ALICE Forward (MC)");
    SetHistogramAttributes(normForwardMC, kRed+3, 21, "ALICE Forward (MC)", 
			   "Normalisation");

    fOutput->Add(dndetaForwardMC);
    fOutput->Add(normForwardMC);
    fOutput->Add(Rebin(dndetaForwardMC));
  }

  if (fSumPrimary) { 
    TH1D* dndetaTruth = fSumPrimary->ProjectionX("dndetaTruth",-1,-1,"e");
    dndetaTruth->Scale(1./nAll, "width");

    SetHistogramAttributes(dndetaTruth, kGray+3, 22, "Monte-Carlo truth");

    fOutput->Add(dndetaTruth);
    fOutput->Add(Rebin(dndetaTruth));
  }
  if (fSumCentral) { 
    TH1D* dndetaCentral = fSumCentral->ProjectionX("dndetaCentral",-1,-1,"e");
    dndetaCentral->Scale(1./nGood, "width");

    SetHistogramAttributes(dndetaCentral, kGray+3, 22, "ALICE Central - track(let)s");

    dndetaCentral->GetXaxis()->SetRangeUser(-1,1);
    fOutput->Add(dndetaCentral);
    fOutput->Add(Rebin(dndetaCentral));
  }

  TNamed* trigString = 
    new TNamed("trigString", AliAODForwardMult::GetTriggerString(fTriggerMask));
  trigString->SetUniqueID(fTriggerMask);
  fOutput->Add(trigString);

  TAxis* vtxAxis = new TAxis(1,fVtxMin,fVtxMax);
  vtxAxis->SetName("vtxAxis");
  vtxAxis->SetTitle(Form("v_{z}#in[%+5.1f,%+5.1f]cm", fVtxMin,fVtxMax));
  fOutput->Add(vtxAxis);

  // If only there was a away to get sqrt{s_NN} and beam type 
  

  PostData(2, fOutput);
}

//________________________________________________________________________
TH1D*
AliForwarddNdetaTask::Rebin(const TH1D* h) const
{
  if (fRebin <= 1) return 0;

  Int_t nBins = h->GetNbinsX();
  if(nBins % fRebin != 0) {
    Warning("Rebin", "Rebin factor %d is not a devisor of current number "
	    "of bins %d in the histogram %s", fRebin, nBins, h->GetName());
    return 0;
  }
    
  // Make a copy 
  TH1D* tmp = static_cast<TH1D*>(h->Clone(Form("%s_rebin%02d", 
					       h->GetName(), fRebin)));
  tmp->Rebin(fRebin);
  tmp->SetDirectory(0);

  // The new number of bins 
  Int_t nBinsNew = nBins / fRebin;
  for(Int_t i = 1;i<= nBinsNew; i++) {
    Double_t content = 0;
    Double_t sumw    = 0;
    Double_t wsum    = 0;
    Int_t    nbins   = 0;
    for(Int_t j = 1; j<=fRebin;j++) {
      Int_t    bin = (i-1)*fRebin + j;
      Double_t c   =  h->GetBinContent(bin);

      if (c <= 0) continue;

      if (fCutEdges) {
	if (h->GetBinContent(bin+1)<=0 || 
	    h->GetBinContent(bin-1)) {
	  Warning("Rebin", "removing bin %d=%f of %s (%d=%f,%d=%f)", 
		  bin, c, h->GetName(), 
		  bin+1, h->GetBinContent(bin+1), 
		  bin-1, h->GetBinContent(bin-1));
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
