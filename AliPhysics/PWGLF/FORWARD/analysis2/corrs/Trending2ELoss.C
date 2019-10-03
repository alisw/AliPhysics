#ifndef __CINT__
# include <TClass.h>
# include <TDirectory.h>
# include <TCollection.h>
# include <TParameter.h>
# include <TH1.h>
# include <TAxis.h>
# include <TError.h>
# include <TGrid.h>
# include <TClonesArray.h>
# include <TFile.h>
# include <TF1.h>
# include <TLine.h>
# include "AliFMDCorrELossFit.h"
# include "AliForwardCorrectionManager.h"
#else
class TClass;
class TCollection;
class TDirectory;
class TH1;
class TAxis;
#endif 

// ===================================================================
Bool_t CheckClass(TObject* o, const TObject* dir,
		  TClass* cls, Bool_t verbose=true)
{
  if (!cls) return true;
  if (!o->IsA()->InheritsFrom(cls)) {
    if (verbose) 
      Warning("CheckClass", "Object \"%s\" from \"%s\" is not a %s but a %s",
	      o->GetName(), dir->GetName(), cls->GetName(), o->ClassName());
    return false;
  }
  return true;
}

// ===================================================================
TObject*
GetObject(TDirectory* dir, const char* name, TClass* cls,
	  Bool_t verbose=true)
{
  if (!dir) {
    if (verbose) Warning("GetObject", "No directory passed");
    return 0;
  }
  TObject* o = dir->Get(name);
  if (!o) {
    if (verbose)  Warning("GetObject", "Object \"%s\" not found in \"%s\"",
			  name, dir->GetName());
    // dir->ls();
    return 0;
  }
  if (!CheckClass(o, dir, cls, verbose)) return 0;
  return o;
}

// ___________________________________________________________________
TObject*
GetObject(const TCollection* dir, const char* name, TClass* cls,
	  Bool_t verbose=true)
{
  if (!dir) {
    if (verbose) Warning("GetObject", "No collection passed");
    return 0;
  }
  TObject* o = dir->FindObject(name);
  if (!o) {
    if (verbose) 
      Warning("GetObject", "Object \"%s\" not found in \"%s\"",
	      name, dir->GetName());
    // dir->ls();
    return 0;
  }
  if (!CheckClass(o, dir, cls,verbose)) return 0;
  return o;
}
// ___________________________________________________________________
TCollection*
GetCollection(TDirectory* dir, const char* name, Bool_t verbose=true)
{
  return static_cast<TCollection*>(GetObject(dir, name, TCollection::Class(),
					     verbose));
}
// ___________________________________________________________________
TCollection*
GetCollection(const TCollection* dir, const char* name, Bool_t verbose=true)
{
  return static_cast<TCollection*>(GetObject(dir, name, TCollection::Class(),
				   verbose));
}
// ___________________________________________________________________
TH1* GetH1(TDirectory* dir, const char* name, Bool_t verbose=true)
{
  return static_cast<TH1*>(GetObject(dir, name, TH1::Class(),verbose));
}
// ___________________________________________________________________
TH1* GetH1(const TCollection* dir, const char* name, Bool_t verbose=true)
{
  return static_cast<TH1*>(GetObject(dir, name, TH1::Class(),verbose));
}
// ___________________________________________________________________
TAxis* GetAxis(TDirectory* dir, const char* name, Bool_t verbose=true)
{
  return static_cast<TAxis*>(GetObject(dir, name, TAxis::Class(),verbose));
}
// ___________________________________________________________________
TAxis* GetAxis(const TCollection* dir, const char* name, Bool_t verbose=true)
{
  return static_cast<TAxis*>(GetObject(dir, name, TAxis::Class(),verbose));
}
// ___________________________________________________________________
Int_t GetInt(TDirectory* dir, const char* name, Int_t def=0,
	     Bool_t verbose=true)
{
  TParameter<int>* p =
    static_cast<TParameter<int>*>(GetObject(dir,name,TParameter<int>::Class(),
					    verbose));
  if (!p) return def;
  return p->GetVal();
}
// ___________________________________________________________________
Int_t GetInt(const TCollection* dir, const char* name,
	     Int_t def=0, Bool_t verbose=true)
{
  TParameter<int>* p =
    static_cast<TParameter<int>*>(GetObject(dir,name,TParameter<int>::Class(),
					    verbose));
  if (!p) return def;
  return p->GetVal();
}
// ___________________________________________________________________
Long_t GetLong(TDirectory* dir, const char* name,
	       Long_t def=0, Bool_t verbose=true)
{
  TParameter<long>* p =
    static_cast<TParameter<long>*>(GetObject(dir,name,
					     TParameter<long>::Class(),
					     verbose));
  if (!p) return def;
  return p->GetVal();
}
// ___________________________________________________________________
Long_t GetLong(const TCollection* dir, const char* name,
	       Long_t def=0, Bool_t verbose=true)
{
  TParameter<long>* p =
    static_cast<TParameter<long>*>(GetObject(dir,name,
					     TParameter<long>::Class(),
					     verbose));
  if (!p) return def;
  return p->GetVal();
}
// ___________________________________________________________________
Double_t GetDouble(TDirectory* dir, const char* name,
		   Double_t def=0, Bool_t verbose=true)
{
  TParameter<double>* p =
    static_cast<TParameter<double>*>(GetObject(dir,name,
					       TParameter<double>::Class(),
					       verbose));
  if (!p) return def;
  return p->GetVal();
}
// ___________________________________________________________________
Double_t GetDouble(const TCollection* dir, const char* name,
		   Double_t def=0, Bool_t verbose=true)
{
  TParameter<double>* p =
    static_cast<TParameter<double>*>(GetObject(dir,name,
					       TParameter<double>::Class(),
					       verbose));
  if (!p) return def;
  return p->GetVal();
}
// ___________________________________________________________________
Bool_t GetBool(TDirectory* dir, const char* name,
	       Bool_t def=false, Bool_t verbose=true)
{
  TParameter<bool>* p =
    static_cast<TParameter<bool>*>(GetObject(dir,name,
					     TParameter<bool>::Class(),
					     verbose));
  if (!p) return def;
  return p->GetVal();
}
// ___________________________________________________________________
Bool_t GetBool(const TCollection* dir, const char* name,
	       Bool_t def=false, Bool_t verbose=true)
{
  TParameter<bool>* p =
    static_cast<TParameter<bool>*>(GetObject(dir,name,
					     TParameter<bool>::Class(),
					     verbose));
  if (!p) return def;
  return p->GetVal();
}

// ===================================================================
Bool_t Fail(const char* msg)
{
  Warning("Trending2ELoss", msg);
  TFile* stamp = TFile::Open("bad.root", "RECREATE");
  TNamed* status = new TNamed("status", "Bad run");
  status->Write();
  stamp->Write();
  stamp->Close();
  return false;
}

// ===================================================================
Bool_t
Trending2ELoss(const char* fileName="trending.root",
	       Double_t minRate=.7,
	       Int_t    maxGap=3)
{
  if (!gROOT->GetClass("AliFMDCorrELossFit")) {
    const char* fwd = "$ALICE_PHYSICS/PWGLF/FORWARD/analysis2";
    gROOT->Macro(Form("%s/scripts/LoadLibs.C", fwd));
  }
  TString fname(fileName);
  if (fname.BeginsWith("alien:")) { 
    TGrid::Connect("alien:");
    if (!gGrid) return Fail(Form("Failed to connect to grid"));    
  }

  TFile* file = TFile::Open(fname, "READ");
  if (!file) return Fail(Form("Failed to open \"%s\"", fileName));  

  TCollection* top = GetCollection(file, "forwardQAResults");
  if (!top) return false;

  TCollection* ins = GetCollection(top, "fmdEventInspector");
  if (!ins) return false;

  UInt_t   run = GetLong(ins, "runNo",     0);
  UShort_t sys = GetInt(ins,  "sys",       0);
  UShort_t sNN = GetInt(ins,  "sNN",       0);
  Short_t  fld = GetInt(ins,  "field",     9999);
  Bool_t   sat = GetBool(ins, "satellite", false);
  Bool_t   mc  = GetBool(ins, "mc",        false);
  TH1*     acc = GetH1(ins, "nEventsAccepted");
  
  if (run <= 0 || sys <= 0 || sNN <= 0 || fld > 100) 
    return Fail(Form("Unknown run (%d) system (%d), energy (%d), or field (%d)",
		     sys, sNN, fld));

  Long_t             minEvents = 10000;
  if      (sys == 1) minEvents = 1000000;
  else if (sys == 3) minEvents = 100000;
  else if (sys == 4) minEvents = 100000;
  if (!acc || acc->GetEntries() < minEvents) 
    return Fail(Form("%09d: Too (%ld) few (<%ld) events for sys=%d",
		     run, acc ? Long_t(acc->GetEntries()) : 0, minEvents, sys));
  
  TCollection* enf = GetCollection(top, "fmdEnergyFitter");
  if (!enf) return false;
  
  UShort_t nps = GetInt(enf, "nParticles", 0);
  Double_t low = GetDouble(enf, "lowCut", 0);
  Double_t mxc = GetDouble(enf, "maxChi2PerNDF", 0);
  Double_t mre = GetDouble(enf, "maxRelPerError", 0);
  Double_t lwt = GetDouble(enf, "minWeight", 0);
  Bool_t   shf = GetBool(enf, "deltaShift", false);
  TH1*     het = GetH1(enf, "etaAxis");
    
  if (nps <= 0) return Fail(Form("%09d: Too (%d) few peaks fitted", run, nps));
  if (low <= 0) return Fail(Form("%09d: Lower bound (%f) too low", run, low));  
  if (mxc <= 0 || mre <= 0 || lwt <= 0)
    return Fail(Form("%09d: Max chi^2/nu (%f), max dp/p (%f), "
		     "least weight (%f) invalid",
		     run, mxc, mre, lwt));
  if (!het) return Fail("No eta axis defined");

  AliFMDCorrELossFit* corr = new AliFMDCorrELossFit();
  corr->SetEtaAxis(*(het->GetXaxis()));
  corr->SetLowCut(low);
  corr->SetBit(AliFMDCorrELossFit::kHasShift, shf);
  
  const TAxis& eta = corr->GetEtaAxis();

  TClonesArray tmp("AliFMDCorrELossFit::ELossFit");
  Int_t   ir = 0;
  TH1* hRate  = new TH1D("rate", "Success rate", 5, .5, 5.5);
  TH1* hTotal = new TH1I("total","Total fits",   5, .5, 5.5);
  TH1* hMax   = new TH1I("max",  "Max gap",      5, .5, 5.5);
  hRate ->SetFillColor(kRed+2);   hRate ->SetFillStyle(3001);
  hTotal->SetFillColor(kGreen+2); hTotal->SetFillStyle(3001);
  hMax  ->SetFillColor(kBlue+2);  hMax  ->SetFillStyle(3001);
  hRate ->SetDirectory(0);   
  hTotal->SetDirectory(0); 
  hMax  ->SetDirectory(0);
  hRate ->GetListOfFunctions()->Add(new TLine(.5,100*minRate,5.5,100*minRate));
  hMax  ->GetListOfFunctions()->Add(new TLine(.5,maxGap,5.5,maxGap));
  hTotal->GetListOfFunctions()->Add(new TLine(.5,15,5.5,15));
  hTotal->GetListOfFunctions()->Add(new TLine(.5,25,5.5,25));

  UInt_t bad = 0;
  for (UShort_t d = 1; d <= 3; d++) {
    UShort_t nq = d/2+1;
    for (UShort_t q = 0; q < nq; q++) {
      Char_t  r   = q == 0 ? 'I' : 'O';
      TString nam = Form("FMD%d%c", d, r);
      TCollection* det = GetCollection(enf, nam);
      if (!det) continue;

      TCollection* eld = GetCollection(det, "elossDists");
      if (!eld) continue;

      Int_t nTotal = 0;
      Int_t nOK    = 0;
      Int_t dist   = 0;
      Int_t max    = 0;
      for (Int_t bin = 1; bin <= eta.GetNbins(); bin++) {
	TH1* dst = GetH1(eld, Form("%s_etabin%03d", nam.Data(), bin), false);
	if (!dst) continue;
	nTotal++;
	
	AliFMDCorrELossFit::ELossFit* fit = 0;
	TList* fcs = dst->GetListOfFunctions();
	TIter  nxt(fcs);
	TF1*   fun = 0;
	Int_t  i   = 0;
	tmp.Clear();
	while ((fun = static_cast<TF1*>(nxt()))) {
	  fit = new (tmp[i++]) AliFMDCorrELossFit::ELossFit(0,*fun);
	  fit->fDet  = d;
	  fit->fRing = r;
	  fit->fBin  = bin;
	  fit->CalculateQuality(mxc, mre, lwt);
	}
	// Sort them 
	tmp.Sort();
	fit = static_cast<AliFMDCorrELossFit::ELossFit*>(tmp.At(i-1));
	if (!fit) {
	  Warning("Trending2ELoss", "No fit found for %s %3d",
		  nam.Data(), bin);
	  continue;
	}
	corr->SetFit(d, r, bin, new AliFMDCorrELossFit::ELossFit(*fit));
	if (fit->GetQuality() >= AliFMDCorrELossFit::kDefaultQuality) {
	  nOK++;
	  max  = TMath::Max(dist,max);
	  dist = 0;
	}
	else {
	  dist++;
	}
      } // for each bin
      Double_t rate = Float_t(nOK)/nTotal;
      Info("Trending2ELoss", "FMD%d%c [%d] %3d/%3d: %5.1f%% (max: %d)",
	   d, r, ir, nOK, nTotal, 100*rate, max);
      if (rate < minRate)           bad |= (1 << ir);
      if (max  > maxGap)            bad |= (1 << ir);
      if (r == 'I' && nTotal < 25)  bad |= (1 << ir);
      else if        (nTotal < 15)  bad |= (1 << ir);
      
      ir++;
      hRate ->GetXaxis()->SetBinLabel(ir, Form("FMD%d%c",d,r));
      hTotal->GetXaxis()->SetBinLabel(ir, Form("FMD%d%c",d,r));
      hMax  ->GetXaxis()->SetBinLabel(ir, Form("FMD%d%c",d,r));
      hRate ->SetBinContent(ir, 100*rate);
      hTotal->SetBinContent(ir, nTotal);
      hMax  ->SetBinContent(ir, max);
    } // for each ring 
  } // for each detector
  // corr->Print("R");

  TFile* diag = TFile::Open("diagnostics.root","RECREATE");
  hRate ->SetMaximum(100); hRate ->SetMinimum(0); hRate ->Write();
  hTotal->SetMaximum(35);  hTotal->SetMinimum(0); hTotal->Write();
  hMax  ->SetMaximum(10);  hTotal->SetMinimum(0); hMax  ->Write();
  diag->Write();
  diag->Close();
  Bool_t testBad = !corr->IsGood(true);
  if (testBad != (bad>0)) {
    Warning("Trending2ELoss", "Mismatch between this (%s) and test (%s)",
	    (bad>0) ? "bad" : "good", testBad ? "bad" : "good");
  }
  if (bad > 0) 
    return Fail(Form("%09d: One or more detectors are bad: 0x%x", run, bad));
  
  AliCorrectionManagerBase& man = AliForwardCorrectionManager::Instance();
  if (!man.Store(corr,run,sys,sNN,fld,mc,sat,"fmd_corrections.root")) 
    return Fail(Form("%09d: Failed to store correction", run));

  return true;
}
//
// EOF
// 

  
