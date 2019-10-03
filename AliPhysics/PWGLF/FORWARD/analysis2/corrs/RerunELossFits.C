/** 
 * Get a collection from a file directory
 * 
 * @param dir   Parent directory 
 * @param name  Name of collection
 * @param verbose Be verbose
 * 
 * @return collection or null
 */
TCollection* GetCollection(TDirectory* dir, const TString& name,
			   Bool_t verbose=true)
{
  if (!dir) { 
    Error("GetCollection", "No parent directory for %s", name.Data());
    return 0;
  }
  TCollection* ret = 0;
  dir->GetObject(name, ret);
  if (!ret) {
    if (verbose) 
      Error("GetCollection", "Couldn't find %s in %s", 
	    name.Data(), dir->GetName());
    return 0;
  }
  return ret;
}
/** 
 * Get an object from a collection.  Optionally, we check that the
 * type of the possibly found object matches the request.
 * 
 * @param parent Parent collection
 * @param name   Name of object 
 * @param cls If specified, check that the found object (if any) is of
 * this class.
 * @param verbose Whether to be verbose 
 * 
 * @return Found object (possibly type-checked) or null
 */
TObject* GetObject(const TCollection* parent, const TString& name, 
		   const TClass* cls=0, Bool_t verbose=true)
{
  if (!parent) {
    if (verbose) 
      Error("GetObject", "No parent collection for %s", name.Data());
    return 0;
  }
  TObject* ret = parent->FindObject(name);
  if (!ret) {
    if (verbose) 
      Error("GetObject", "Couldn't find %s in %s", 
	    name.Data(), parent->GetName());
    return 0;
  }
  if (cls && !ret->IsA()->InheritsFrom(cls)) { 
    if (verbose) 
      Error("GetObject", "%s in %s is a %s, not a %s", name.Data(),
	    parent->GetName(), ret->ClassName(), cls->GetName());
    return 0;
  }
  return ret;
}
void RemoveObject(TCollection* parent, const TString& name)
{
  while (parent->GetEntries() > 0) {
    TObject* o = GetObject(parent, name, 0, false);
    if (!o) return;
    Info("RemoveObject", "Removing %s from %s", name.Data(), parent->GetName());
    TObject* ret = parent->Remove(o);
    if (!ret || ret != o)
      Warning("RemoveObject", "Failed to remove %s from %s (%p %p)",
	      name.Data(), parent->GetName(), o, ret);
    // delete o;
  }
}

/** 
 * Get a collection contained in another collection
 * 
 * @param parent Parent collection
 * @param name   Name of collection to find 
 * @param verbose Be verbose
 * 
 * @return Found collection or null
 */
TCollection* GetCollection(const TCollection* parent, const TString& name,
			   Bool_t verbose=true)
{
  TObject* o = GetObject(parent, name, TCollection::Class(), verbose);
  if (!o) return 0;
  return static_cast<TCollection*>(o);
}
 
void CleanCollection(TCollection* c)
{  
  const char* other[] = { "fmdESDFixer",
			  "fmdSharingFilter",
			  "fmdDensityCalculator",
			  0 };
  const char** ptr = other;
  while (*ptr) { RemoveObject(c, *ptr); ptr++; } 

  TCollection* ef = GetCollection(c, "fmdEnergyFitter");
  const char* stacks[] = { "chi2", "c", "delta", "xi", "sigma", "sigman",
			   "a2", "a3", "a4", "a5", 0 };
  ptr = stacks;
  while (*ptr) { RemoveObject(ef, *ptr); ptr++; }

  const char* dets[] = { "FMD1I", "FMD2I", "FMD2O", "FMD3I", "FMD3O", 0 };
  ptr = dets;
  while (*ptr) {
    TCollection* det = GetCollection(ef, *ptr);
    if (!det) {
      ptr++;
      continue;
    }
    const char* subs[] = { "elossDists", "elossResults", "elossResiduals", 0 };
    const char** sub = subs;
    while (*sub) { RemoveObject(det, *sub); sub++;  }

    // TH1* eloss = static_cast<TH1*>(GetObject(det,"eloss",TH1::Class()));
    // if (eloss) Printf("%s %d entries", *ptr, eloss->GetEntries());
    ptr++;
  }
  // c->ls();
}

/** 
 * Re-run the energy loss fitter on a merged output file 
 * 
 * @param input  File name of merged output file
 * @param output If specified, the file the new results are written
 * to.  If this is not specified, it defaults to the name of the input
 * file with "_rerun" attached to the base name
 *
 * @param forceSet  Forcibly set things 
 * @param input     Input file 
 * @param output    Output file 
 * @param shift     Enable shift 
 * @param flags     0x1: residuals, 0x2: debug
 */
void RerunELossFits(Bool_t forceSet=false, 
		    const TString& input="forward_eloss.root", 
		    Bool_t shift=true,
		    const TString& output="",
		    UShort_t       flags=0x1)
{
  const char* fwd = "$ALICE_PHYSICS/PWGLF/FORWARD/analysis2";
  gROOT->Macro(Form("%s/scripts/LoadLibs.C", fwd));

  TFile*  inFile  = 0;
  TFile*  outFile = 0;
  TString outName(output);
  if (outName.IsNull()) {
    outName = input;
    outName.ReplaceAll(".root", "_rerun.root");
  }
  Bool_t  allOk = false;
  try {
    // --- Open input file ---------------------------------------------
    inFile = TFile::Open(input, "READ");
    if (!inFile)
      throw TString::Format("Failed to open %s", input.Data());

    // --- InFiled input collections --------------------------------------
    TCollection*   inFwdSum = GetCollection(inFile, "ForwardELossSums");
    if (!inFwdSum) inFwdSum = GetCollection(inFile, "forwardQAResults");
    if (!inFwdSum) throw new TString("Cannot proceed without sums");
    inFwdSum->SetName("ForwardELossSums");
    
    TCollection*   inFwdRes = GetCollection(inFile, "ForwardELossResults");
    if (!inFwdRes) inFwdRes = GetCollection(inFile, "forwardQAResults");
    if (!inFwdRes) {
      inFwdRes = 
	static_cast<TCollection*>(inFwdSum->Clone("ForwardELossResults"));
      // throw new TString("Cannot proceed with merged list");
    }
    inFwdRes->SetName("ForwardELossResults");
  
    TCollection* inEFSum = GetCollection(inFwdRes, "fmdEnergyFitter");
    if (!inEFSum) throw new TString("Cannot proceed without accumulated data");
    
    TCollection* inEFRes = GetCollection(inFwdRes, "fmdEnergyFitter");
    if (!inEFRes) throw new TString("Cannot proceed without previous results");

    TCollection* inESSum = GetCollection(inFwdRes, "fmdEventInspector");
    Int_t  sys = 0;
    Long_t nEvAcc = 0;
    if (inESSum) {
      TObject* oSys = GetObject(inESSum, "sys", TParameter<int>::Class());
      if (oSys) sys = (static_cast<TParameter<int>*>(oSys))->GetVal();
      TH1* oEvAcc = static_cast<TH1*>(GetObject(inESSum, "nEventsAccepted",
						TH1::Class()));
      if (oEvAcc) nEvAcc = oEvAcc->GetEntries();
    }
    if (nEvAcc > 0 && sys > 0) {
      Long_t minEvents = 0;
      switch (sys) {
      //case 1: minEvents = 1000000; break; // At least 1M
      case 1: minEvents = 500000; break; // At least 1M
      case 2: minEvents = 10000;   break; // At least 10k
      case 3:                             // Fall-through
      case 4: minEvents = 100000;  break; // At least 100k
      }
      if (nEvAcc < minEvents) {
	Error("RerunELossFits", "Too few events (%ld<%ld) for %d",
	      nEvAcc, minEvents, sys);
	return;
      }
    }
    // --- Open output file --------------------------------------------
    outFile = TFile::Open(outName, "RECREATE");
    if (!outFile) 
      throw TString::Format("Failed to open %s", outName.Data());

    // --- Make our fitter object ------------------------------------
    AliFMDEnergyFitter* fitter = new AliFMDEnergyFitter("energy");
    fitter->SetDoFits(true);
    fitter->SetEnableDeltaShift(shift);
    fitter->Init();
    if (forceSet || !fitter->ReadParameters(inEFSum)) {
      Printf("Forced settings");

      const TAxis* etaAxis = static_cast<TAxis*>(GetObject(inEFSum,"etaAxis"));
      if (!etaAxis) throw new TString("Cannot proceed without eta axis");
      fitter->SetEtaAxis(*etaAxis);
      
      // Set maximum energy loss to consider 
      fitter->SetMaxE(15); 
      // Set number of energy loss bins 
      fitter->SetNEbins(500);
      // Set whether to use increasing bin sizes 
      // fitter->SetUseIncreasingBins(true);
      // Set whether to do fit the energy distributions 
      fitter->SetDoFits(kTRUE);
      // Set whether to make the correction object 
      fitter->SetDoMakeObject(kTRUE);
      // Set the low cut used for energy
      fitter->SetLowCut(0.4);
      // Set the number of bins to subtract from maximum of distributions
      // to get the lower bound of the fit range
      fitter->SetFitRangeBinWidth(4);
      // Set the maximum number of landaus to try to fit (max 5)
      fitter->SetNParticles(5);
      // Set the minimum number of entries in the distribution before
      // trying to fit to the data - 10k seems the least we can do
      fitter->SetMinEntries(10000);
      // fitter->SetMaxChi2PerNDF(10);
      // Enable debug 
    }
    // Set the maximum number of landaus to try to fit (max 5)
    Int_t maxPart = sys == 1 ? 3 : 5;
    fitter->SetNParticles(maxPart);
    fitter->SetDoMakeObject(true);
    fitter->SetMinEntries(10000);
    if (fitter->GetLowCut() < .5) fitter->SetLowCut(0.55);
    if (flags & 0x2)  fitter->SetDebug(3);
    if (flags & 0x1)
      fitter->SetStoreResiduals(AliFMDEnergyFitter::kResidualSquareDifference);
    // fitter->SetRegularizationCut(1e8); // Lower by factor 3
    // Set the number of bins to subtract from maximum of distributions
    // to get the lower bound of the fit range
    // fitter->SetFitRangeBinWidth(2);
    // Skip all of FMD2 and 3 
    // fitter->SetSkips(AliFMDEnergyFitter::kFMD2|AliFMDEnergyFitter::kFMD3);

    // --- Write copy of sum collection to output --------------------
    TCollection* outFwdSum = static_cast<TCollection*>(inFwdSum->Clone());
    outFile->cd();
    CleanCollection(outFwdSum);
    TCollection* outEFSum = GetCollection(outFwdSum, "fmdEnergyFitter");
    if (outEFSum) {
      TObject* np = GetObject(outEFSum,"nParticles",TParameter<int>::Class());
      if (np) (static_cast<TParameter<int>*>(np))->SetVal(maxPart);
    }
    outFwdSum->Write(inFwdSum->GetName(), TObject::kSingleKey);
    outFwdSum->ls();
    
    // --- Now do the fits -------------------------------------------
    fitter->Print();
    TStopwatch timer;
    timer.Start();
    fitter->Fit(static_cast<TList*>(outFwdSum));
    timer.Print();

    // --- Copy full result folder -----------------------------------
    TCollection* outFwdRes = static_cast<TCollection*>(inFwdRes->Clone());
    CleanCollection(outFwdRes);
    // Remove old fits
    while (true) {
      TCollection* outEFRes = GetCollection(outFwdRes,"fmdEnergyFitter",false);
      if (!outEFRes) break;
      outFwdRes->Remove(outEFRes);
    } 
    // outFwdRes->ls();
    // Make our new fit results folder, and add it to results folder
    TCollection* tmp = GetCollection(outFwdSum, "fmdEnergyFitter");
    outEFRes = static_cast<TCollection*>(tmp->Clone());
    outEFRes->Add(new TNamed("refitted", "Refit of the data"));
    outFwdRes->Add(outEFRes);
    // outFwdRes->ls();

    // --- Write out new results folder ------------------------------
    outFile->cd();
    outFwdRes->Write(inFwdRes->GetName(), TObject::kSingleKey);
    Printf("Wrote results to \"%s\" (%s)", outName.Data(), outFile->GetName());
    allOk = true;
  }
  catch (const TString* e) {
    Error("RerunELossFits", e->Data());
  }
  catch (const TString& e) {
    Error("RerunELossFits", e.Data());
  }
  if (inFile)  inFile->Close();
  if (outFile) {
    Printf("Wrote new output to \"%s\"", outName.Data());
    outFile->Close();
  }
    
  if (allOk)
    gROOT->Macro(Form("%s/corrs/DrawCorrELoss.C(false,false,\"%s\")", 
		      fwd, outName.Data()));
}


