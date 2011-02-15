//====================================================================
#include "AliCentraldNdetaTask.h"
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
#include "AliAODForwardMult.h"
#include <AliAODCentralMult.h>

//____________________________________________________________________
AliCentraldNdetaTask::AliCentraldNdetaTask()
  : AliAnalysisTaskSE(), 
    fSumCentral(0),	//  Sum of histograms 
    fSumCentralMC(0),	//  Sum of MC histograms (if any)
    fSums(0),		// Container of sums 
    fOutput(0),		// Container of outputs 
    fTriggers(0),	// Histogram of triggers 
    fVtxMin(0),		// Minimum v_z
    fVtxMax(0),		// Maximum v_z
    fTriggerMask(0),    // Trigger mask 
    fRebin(0),       	// Rebinning factor 
    fCutEdges(false)
{}

//____________________________________________________________________
AliCentraldNdetaTask::AliCentraldNdetaTask(const char* /* name */)
  : AliAnalysisTaskSE("Central"), 
    fSumCentral(0),	//  Sum of histograms 
    fSumCentralMC(0),	//  Sum of MC histograms (if any)
    fSums(0),		// Container of sums 
    fOutput(0),		// Container of outputs 
    fTriggers(0),	// Histogram of triggers 
    fVtxMin(-10),	// Minimum v_z
    fVtxMax(10),	// Maximum v_z
    fTriggerMask(AliAODForwardMult::kInel), 
    fRebin(5),		// Rebinning factor 
    fCutEdges(false)
{
  // Output slot #1 writes into a TH1 container
  DefineOutput(1, TList::Class()); 
  DefineOutput(2, TList::Class()); 
}

//____________________________________________________________________
AliCentraldNdetaTask::AliCentraldNdetaTask(const AliCentraldNdetaTask& o)
  : AliAnalysisTaskSE(o), 
    fSumCentral(o.fSumCentral),	// TH2D* -  Sum of histograms 
    fSumCentralMC(o.fSumCentralMC),// TH2D* -  Sum of MC histograms (if any)
    fSums(o.fSums),		// TList* - Container of sums 
    fOutput(o.fOutput),		// TList* - Container of outputs 
    fTriggers(o.fTriggers),	// TH1D* - Histogram of triggers 
    fVtxMin(o.fVtxMin),		// Double_t - Minimum v_z
    fVtxMax(o.fVtxMax),		// Double_t - Maximum v_z
    fTriggerMask(o.fTriggerMask),// Int_t - Trigger mask 
    fRebin(o.fRebin),		// Int_t - Rebinning factor 
    fCutEdges(o.fCutEdges)	// Bool_t - Whether to cut edges when rebinning
{}

//____________________________________________________________________
AliCentraldNdetaTask::~AliCentraldNdetaTask()
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
}

//________________________________________________________________________
void 
AliCentraldNdetaTask::SetTriggerMask(const char* mask)
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
AliCentraldNdetaTask::UserCreateOutputObjects()
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
AliCentraldNdetaTask::CloneHist(TH2D* in, const char* name) 
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
AliCentraldNdetaTask::UserExec(Option_t *) 
{
  // Main loop
  AliAODEvent* aod = dynamic_cast<AliAODEvent*>(InputEvent());
  if (!aod) {
    AliError("Cannot get the AOD event");
    return;
  }  

  // Get objects from the event structure 
  TObject* oCentral   = aod->FindListObject("CentralClusters");
  TObject* oCentralMC = aod->FindListObject("CentralClustersMC");
  TObject* oForward   = aod->FindListObject("Forward");

  // We should have a central object at least 
  if (!oCentral) { 
    AliWarning("No Central object found AOD");
    return;
  }

  // Cast to good types 
  AliAODCentralMult* central   = static_cast<AliAODCentralMult*>(oCentral);
  AliAODCentralMult* centralMC = static_cast<AliAODCentralMult*>(oCentralMC);
  AliAODForwardMult* forward   = static_cast<AliAODForwardMult*>(oForward);
  
  // Create our sum histograms 
  if (!fSumCentral) fSumCentral =CloneHist(&central->GetHistogram(),"central");
  if (!fSumCentralMC && centralMC) 
    fSumCentralMC = CloneHist(&centralMC->GetHistogram(),"centralMC");

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
  fSumCentral->Add(&(central->GetHistogram()));
  if (fSumCentralMC) fSumCentralMC->Add(&(centralMC->GetHistogram()));

  PostData(1, fSums);
}

//________________________________________________________________________
void 
AliCentraldNdetaTask::SetHistogramAttributes(TH1D* h, Int_t colour, Int_t marker,
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
AliCentraldNdetaTask::Terminate(Option_t *) 
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

  fSumCentral     = static_cast<TH2D*>(fSums->FindObject("central"));
  fSumCentralMC   = static_cast<TH2D*>(fSums->FindObject("centralMC"));
  fTriggers       = static_cast<TH1D*>(fSums->FindObject("triggers"));

  if (!fTriggers) { 
    AliError("Couldn't find histogram 'triggers' in list");
    return;
  }
  if (!fSumCentral) { 
    AliError("Couldn't find histogram 'central' in list");
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
  TH1D* normCentral   = fSumCentral->ProjectionX("normCentral", 0, 1, "");
  // Project onto eta axis - _ignoring_underflow_bins_!
  TH1D* dndetaCentral = fSumCentral->ProjectionX("dndetaCentral", 1, -1, "e");
  // Normalize to the acceptance 
  dndetaCentral->Divide(normCentral);
  // Scale by the vertex efficiency 
  dndetaCentral->Scale(vtxEff, "width");

  SetHistogramAttributes(dndetaCentral, kRed+1, 20, "ALICE Central");
  SetHistogramAttributes(normCentral, kRed+1, 20, "ALICE Central", 
			 "Normalisation");

  fOutput->Add(fTriggers->Clone());
  fOutput->Add(dndetaCentral);
  fOutput->Add(normCentral);
  fOutput->Add(Rebin(dndetaCentral));

  if (fSumCentralMC) { 
    // Get acceptance normalisation from underflow bins 
    TH1D* normCentralMC   = 
      fSumCentralMC->ProjectionX("normCentralMC", 0, 1, "");
    // Project onto eta axis - _ignoring_underflow_bins_!
    TH1D* dndetaCentralMC = 
      fSumCentralMC->ProjectionX("dndetaCentralMC", 1, -1, "e");
    // Normalize to the acceptance 
    dndetaCentralMC->Divide(normCentralMC);
    // Scale by the vertex efficiency 
    dndetaCentralMC->Scale(vtxEff, "width");

    SetHistogramAttributes(dndetaCentralMC, kRed+3, 21, "ALICE Central (MC)");
    SetHistogramAttributes(normCentralMC, kRed+3, 21, "ALICE Central (MC)", 
			   "Normalisation");

    fOutput->Add(dndetaCentralMC);
    fOutput->Add(normCentralMC);
    fOutput->Add(Rebin(dndetaCentralMC));
  }

  PostData(2, fOutput);
}

//________________________________________________________________________
TH1D*
AliCentraldNdetaTask::Rebin(const TH1D* h) const
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
