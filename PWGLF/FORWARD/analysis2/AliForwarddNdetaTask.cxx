//====================================================================
#include "AliForwarddNdetaTask.h"
#include <TMath.h>
#include <TH2D.h>
#include <TH1D.h>
#include <THStack.h>
#include <TList.h>
#include <TFile.h>
#include <TF1.h>
#include <TSystem.h>
#include <TGraphAsymmErrors.h>
#include <AliAnalysisManager.h>
#include <AliAODEvent.h>
#include <AliAODHandler.h>
#include <AliAODInputHandler.h>
#include "AliForwardUtil.h"
#include "AliAODForwardMult.h"

//____________________________________________________________________
AliForwarddNdetaTask::AliForwarddNdetaTask()
  : AliBasedNdetaTask()
{
  //
  // Constructor 
  // 
  DGUARD(fDebug, 3, "Default CTOR of AliForwarddNdetaTask");
}

//____________________________________________________________________
AliForwarddNdetaTask::AliForwarddNdetaTask(const char* /* name */)
  : AliBasedNdetaTask("Forward")
{
  // 
  // Constructor
  // 
  // Paramters
  //   name    Name of task 
  // SetTitle("FMD");
  DGUARD(fDebug, 3, "Named CTOR of AliForwarddNdetaTask");
}

//____________________________________________________________________
AliForwarddNdetaTask::AliForwarddNdetaTask(const AliForwarddNdetaTask& o)
  : AliBasedNdetaTask(o)
{
  // 
  // Copy constructor
  // 
  DGUARD(fDebug, 3, "Copy CTOR of AliForwarddNdetaTask");
}

//____________________________________________________________________
AliBasedNdetaTask::CentralityBin*
AliForwarddNdetaTask::MakeCentralityBin(const char* name, Float_t l,Float_t h) 
  const 
{
  // 
  // Make a new centrality bin
  // 
  // Parameters:
  //    name   Histogram names
  //    l      Lower cut
  //    h      Upper cut
  // 
  // Return:
  //    Newly allocated object (of our type)
  //
  DGUARD(fDebug, 3,
	 "Make a centrality bin for AliForwarddNdetaTask: %s [%5.1f%%,%5.1f%%]",
	 name, l, h);
  return new AliForwarddNdetaTask::CentralityBin(name, l, h);
}


//____________________________________________________________________
TH2D*
AliForwarddNdetaTask::GetHistogram(const AliAODEvent& aod, Bool_t mc)
{
  // 
  // Retrieve the histogram 
  // 
  // Parameters:
  //    aod AOD event 
  //    mc  Whether to get the MC histogram or not
  // 
  // Return:
  //    Retrieved histogram or null
  //
  // We should have a forward object at least 
  AliAODForwardMult* forward = GetForward(aod, mc, !mc);
  if (!forward) return 0;
  return &(forward->GetHistogram());
}
//____________________________________________________________________
void
AliForwarddNdetaTask::CheckEventData(Double_t vtx, 
				     TH2*     data, 
				     TH2*     dataMC)
{
  // Check if this is satellite
  // if (!fSatelliteVertices) return;

  // If we don't care about satellites, get out
  if (vtx < fIPzAxis.GetXmin() || vtx > fIPzAxis.GetXmax()) return;

  // Check if this is a satelllite 
  Double_t aVtx = TMath::Abs(vtx);
  if (aVtx < 37.5 || aVtx > 400) return;

  DMSG(fDebug,0,"Got satelitte vertex %f", vtx);
  
  TH2* hists[] = { data, dataMC };

  // In satellite vertices FMD2i is cut away manually at this point
  // for certain vertices. It could be done in the ESDs, but as of
  // this writing not for specific vertices.
  // 
  // cholm comment: It would be difficult to setup the filter in the
  // reconstruction pass, but it could perhaps be done in the AOD
  // filtering.
  // 
  // This is what was done for
  // the Pb-Pb paper (arXiv:1304.0347).
  for (Int_t iX = 0; iX<=data->GetNbinsX(); iX++) {
    // Do all checks up front - as soon as we can - branching is
    // expensive!
    Double_t x    = data->GetXaxis()->GetBinCenter(iX);
    Bool_t   zero = false;
    if (((vtx >  60 && vtx <  90) && x < 3) ||
	((vtx > 330 && vtx < 350) && x > -2.5) ||
	((vtx < 100 || vtx > 305) && TMath::Abs(x) < 4.5) || 
	(vtx < 50                 && TMath::Abs(x) < 4.75))
      zero = true;
    if (!zero) continue;
    for (Int_t iH = 0; iH < 2; iH++) {
      if (!hists[iH]) continue;
      // if (iX > hists[iH]->GetNbinsX()+1) continue;
      // Also zero coverage and phi acceptance for this 
      for (Int_t iY = 0; iY<=hists[iH]->GetNbinsY()+1; iY++) {	  
	hists[iH]->SetBinContent(iX, iY, 0);
	hists[iH]->SetBinError(iX, iY, 0);
      }
    }
  }

  if (fCorrEmpty) {
    DMSG(fDebug,1,"Correcting with corrEmpty=true");
    // Now, since we have some dead areas in FMD2i (sectors 16 and
    // 17), we need to remove the corresponding bins from the
    // histogram. However, it is not obvious which bins (in eta) to
    // remove, so remove everything starting from the most negative to
    // the middle of the histogram.
    // 
    // This hack was first introduced by HHD, but was done at the end of
    // the event processing (CentralityBin::MakeResults).  That is,
    // however, not very practical, as we'd like to normalize to the phi
    // acceptance rather than the eta coverage and then correct for
    // empty bins. Since the only way to really update the phi
    // acceptance stored in the overflow bin is on the event level, we
    // should really do it here.
    const Int_t phiBin1 = 17; // Sector 16
    const Int_t phiBin2 = 18; // Sector 17
    for (Int_t iH = 0; iH < 2; iH++) { 
      if (!hists[iH]) continue;
      
      Int_t midX = hists[iH]->GetNbinsX() / 2;
      // Int_t nY   = hists[iH]->GetNbinsY();
      for (Int_t i = 1; i <= midX; i++) { 
	hists[iH]->SetBinContent(i, phiBin1, 0);
	hists[iH]->SetBinContent(i, phiBin2, 0);
	hists[iH]->SetBinError(i, phiBin1, 0);
	hists[iH]->SetBinError(i, phiBin2, 0);
	
	// Here, we should also modify the overflow bin to reflect the
	// new phi acceptance.  First get the old phi acceptance -
	// then multiply this on the number of bins. This gives us -
	// roughly - the number of sectors we had.  Then take out two
	// from that number, and then calculate the new phi
	// Acceptance. Note, if the sectors where already taken out in
	// the AOD production, we _will_ end up with a wrong number,
	// so we should _not_ do that in the AOD production.  This is
	// tricky and may not work at all.  For now, we should rely on
	// the old way of correcting to the eta coverage and
	// correcting for empty bins.
      }
    }
  }
}
//____________________________________________________________________
Bool_t
AliForwarddNdetaTask::LoadEmpirical(const char* prx)
{
  TString path(prx);
  TUrl    empUrl;
  TFile*  empFile = 0;

  {
    AliForwardUtil::SuppressGuard g(5001);  
    if (gSystem->ExpandPathName(path)) {
      // Expand with TString argument return 0 on success, 1 on failure
      return false;
    }
    if (!path.Contains("empirical"))
      path = gSystem->ConcatFileName(path.Data(), "empirical.root");
    empUrl.SetUrl(path);
    if (!empUrl.GetAnchor() || empUrl.GetAnchor()[0] == '\0')
      empUrl.SetAnchor("default");
    empFile = TFile::Open(empUrl.GetUrl());
    if (!empFile) {
      DMSG(fDebug,1,"%s not found", empUrl.GetUrl());
      return false;
    }
  }
  DMSG(fDebug,0,"Got empirical file %s", empUrl.GetUrl());
  
  TString     base(GetName()); base.ReplaceAll("dNdeta", "");
  TString     empAnch = empUrl.GetAnchor();
  TObject*    empObj  = empFile->Get(Form("%s/%s",base.Data(),empAnch.Data()));
  if (!(empObj &&
	(empObj->IsA()->InheritsFrom(TH1::Class()) ||
	 empObj->IsA()->InheritsFrom(TF1::Class()) ||
	 empObj->IsA()->InheritsFrom(TGraphAsymmErrors::Class())))) {
    Warning("LoadEmpirical", "Didn't get Forward/%s from %s",
	    empAnch.Data(), empUrl.GetUrl());
    return false;
  }
  DMSG(fDebug,0,"Got empirical correction %s [%s]", empObj->GetName(),
       empObj->ClassName());
  // Store correction in output list 
  static_cast<TNamed*>(empObj)->SetName("empirical");
  fResults->Add(empObj);

  if (!empObj->IsA()->InheritsFrom(TH1::Class()))
    return true;

  // Release from directory 
  TH1* h  = static_cast<TH1*>(empObj);
  h->SetDirectory(0);

  // Check for IP delta 
  TH1* xy = static_cast<TH2*>(fSums->FindObject("vertexAccXY"));
  Double_t delta = 0;
  if (xy) {
    Double_t       meanIpx = xy->GetMean(1);
    Double_t       meanIpy = xy->GetMean(2);
    const Double_t refX    = -0.004;
    const Double_t refY    = 0.184;
    Double_t       dx   = (meanIpx - refX);
    Double_t       dy   = (meanIpy - refY);
    Info("LoadEmpirical","Shifts (%f-%f)=%f, (%f-%f)=%f",
	 meanIpx, refX, dx,
	 meanIpy, refY, dy);
    if (TMath::Abs(dx) > 1e-3 || TMath::Abs(dy) > 1e-3)
      delta = TMath::Sqrt(dx*dx+dy*dy);
    // Store delta in output list
    fResults->Add(AliForwardUtil::MakeParameter("deltaIP", delta));
  }
  if (delta > 0.2) {
    // Only correct if delta is larger than 2mm (2% correction)
    TF1* f = new TF1("deltaCorr", "1+[2]*([0]+(x<[1])*pow([0]*(x-[1]),2))");
    f->SetParNames("\\delta","\\eta_{0}","a");
    f->SetParameter(0,delta);
    f->SetParameter(1,-2.0);
    f->SetParameter(2,.10); //TMath::Sqrt(2)); // 0.5);
    Info("LoadEmpirical","Appying correction for IP delta=%f", delta);
    for (Int_t i = 1; i <= h->GetNbinsX(); i++) {
      Double_t c   = h->GetBinContent(i);
      if (c < 1e-6) continue;
      
      Double_t e   = h->GetBinError(i);
      Double_t eta = h->GetXaxis()->GetBinCenter(i);
      Double_t cor = f->Eval(eta);
      // Info("", "%5.2f -> %7.4f", eta, cor);
      h->SetBinContent(i, c*cor);
      h->SetBinError(i, e*cor);
    }
    // Adding correction to result list
    fResults->Add(f);
  }
  return true;
}

//____________________________________________________________________
Bool_t
AliForwarddNdetaTask::Finalize()
{
  // See if we can find the empirical correction so that the bins may
  // apply it
  TString oadb(gSystem->ConcatFileName(AliAnalysisManager::GetOADBPath(),
				       "PWGLF/FORWARD/CORRECTIONS/data"));
  const char* dirs[] = {
    "${PWD}",
    "${FWD}",
    oadb.Data(),
    "${OADB_PATH}/PWGLF/FORWARD/EMPIRICAL",
    "${ALICE_PHYSICS}/OADB/PWGLF/FORWARD/EMPIRICAL",
    0
  };
  const char** pdir = dirs;
  Bool_t       ok   = false;
  while (*pdir) {
    const char*  fns[] = { "", "empirical_000138190.root", 0 };
    const char** pfn   = fns;
    const char*  dir   = *pdir;
    pdir++;
    while (*pfn) {
      AliForwardUtil::SuppressGuard g(1000);
      TString path(gSystem->ConcatFileName(dir, *pfn));
      pfn++;
      const char*  ancs[] = { "param", "default", 0 };
      const char** pan    = ancs;
      while (*pan) {
	const char* anch = *pan;
	pan++;
	TString u(path);  u.Append("#"); u.Append(anch);
	if ((ok = LoadEmpirical(u))) break;
      }
      if (ok) break;
    }
    if (ok) break;
  }
  return AliBasedNdetaTask::Finalize();
}

//========================================================================
TH1*
AliForwarddNdetaTask::CentralityBin::EmpiricalCorrection(TList* results)
{
  TList* out = fOutput;

  TString base(GetName()); base.ReplaceAll("dNdeta", "");
  TString hName(Form("dndeta%s",base.Data()));
  TH1* h = static_cast<TH1*>(out->FindObject(hName));
  if (!h) {
    Warning("End", "%s not found in %s",
	    hName.Data(), out->GetName());
    out->ls();
    return 0;
  }

  TObject* o = static_cast<TH1*>(results->FindObject("empirical"));
  if (!o) {
    Info("EmpiricalCorrection", "Empirical not found in %s",
	 results->GetName());
    return 0;
  }

  Info("EmpiricalCorrection", "Correcting %s/%s with %s [%s]",
       out->GetName(), h->GetName(), o->GetName(), o->ClassName());

  // Make a clone 
  h = static_cast<TH1*>(h->Clone(Form("%sEmp", h->GetName())));
  h->SetDirectory(0);
  
  if (o->IsA()->InheritsFrom(TGraphAsymmErrors::Class())) {
    TGraphAsymmErrors* empCorr = static_cast<TGraphAsymmErrors*>(o);
    TAxis* xAxis = h->GetXaxis();
    for (Int_t i = 1; i <= xAxis->GetNbins(); i++) {
      Double_t x = xAxis->GetBinCenter(i);
      Double_t y = h->GetBinContent(i);
      Double_t c = empCorr->Eval(x);
      h->SetBinContent(i, y / c);
    }
  }
  else if (o->IsA()->InheritsFrom(TF1::Class())) {
    TF1* empCorr = static_cast<TF1*>(o);
    h->Divide(empCorr);
  }
  else if (o->IsA()->InheritsFrom(TH1::Class())) {
    TH1* empCorr = static_cast<TH1*>(o);
    h->Divide(empCorr);
  }
  else { 
    Warning("CorrectEmpirical", 
	    "Don't know how to apply a %s as an empirical correction",
	    o->IsA()->GetName());
    delete h;
    return 0;
  }
  // Adding the corrected histogram to output 
  out->Add(h);
  
  return h;
}
//____________________________________________________________________
bool
AliForwarddNdetaTask::CentralityBin::End(TList*      sums, 
					 TList*      results,
					 UShort_t    scheme,
					 Double_t    trigEff,
					 Double_t    trigEff0,
					 Bool_t      rootProj,
					 Bool_t      corrEmpty, 
					 Int_t       triggerMask,
					 Int_t       marker,
					 Int_t       color,
					 TList*      mclist,
					 TList*      truthlist )
{
  DGUARD(fDebug, 1,"In End of %s with corrEmpty=%d, rootProj=%d", 
	 GetName(), corrEmpty, rootProj);
  if (!AliBasedNdetaTask::CentralityBin::End(sums, results, scheme, trigEff, 
					     trigEff0, rootProj, corrEmpty,
					     triggerMask,
					     marker, color, mclist, 
					     truthlist))
    return false;

  TH1* h = EmpiricalCorrection(results);
  Info("End", "Applied empirical correction: %p (%s)",
       h, h ? h->GetName() : "");
  
  if (!IsAllBin()) return true;

  THStack* res = 0;
  {
    if (gSystem->AccessPathName("forward.root")) return true;

    TFile* file = TFile::Open("forward.root", "READ");
    if (!file) return false;
    
    TList* forward = static_cast<TList*>(file->Get("ForwardSums"));
    if (!forward) { 
      AliError("List Forward not found in forward.root");
      return true;
    }
    TList* rings = static_cast<TList*>(forward->FindObject("ringResults"));
    if (!rings) { 
      AliError("List ringResults not found in forward.root");
      return true;
    }
    res = static_cast<THStack*>(rings->FindObject("all"));
    if (!res) { 
      AliError(Form("Stack all not found in %s", rings->GetName()));
      return true;
    }
  }
  if (!fTriggers) { 
    AliError("Triggers histogram not set");
    return false;
  }

  Double_t ntotal   = 0;
  Double_t epsilonT = trigEff;
#if 0
  // TEMPORARY FIX
  if (triggerMask == AliAODForwardMult::kNSD) {
    // This is a local change 
    epsilonT = 0.92; 
    AliWarning(Form("Using hard-coded NSD trigger efficiency of %f",epsilonT));
  }
#endif
  AliInfo("Adding per-ring histograms to output");
  TString text;
  Double_t scaler = Normalization(*fTriggers, scheme, epsilonT, ntotal, &text);
  TIter next(res->GetHists());
  TH1*  hist = 0;
  while ((hist = static_cast<TH1*>(next()))) hist->Scale(scaler);
  res->SetName("dndetaRings");
  fOutput->Add(res);
  fOutput->Add(new TNamed("normCalc", text.Data()));

  return true;
}

//________________________________________________________________________
//
// EOF
//
