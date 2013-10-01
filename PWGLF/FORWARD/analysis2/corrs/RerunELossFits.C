/** 
 * Get a collection from a file directory
 * 
 * @param dir   Parent directory 
 * @param name  Name of collection
 * 
 * @return collection or null
 */
TCollection* GetCollection(TDirectory* dir, const TString& name)
{
  if (!dir) { 
    Error("GetCollection", "No parent directory for %s", name.Data());
    return 0;
  }
  TCollection* ret = 0;
  dir->GetObject(name, ret);
  if (!ret) { 
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
 * 
 * @return Found object (possibly type-checked) or null
 */
TObject* GetObject(const TCollection* parent, const TString& name, 
		   const TClass* cls=0)
{
  if (!parent) { 
    Error("GetObject", "No parent collection for %s", name.Data());
    return 0;
  }
  TObject* ret = parent->FindObject(name);
  if (!ret) {
    Error("GetObject", "Couldn't find %s in %s", 
	  name.Data(), parent->GetName());
    return 0;
  }
  if (cls && !ret->IsA()->InheritsFrom(cls)) { 
    Error("GetObject", "%s in %s is a %s, not a %s", name.Data(),
	  parent->GetName(), ret->ClassName(), cls->GetName());
    return 0;
  }
  return ret;
}
/** 
 * Get a collection contained in another collection
 * 
 * @param parent Parent collection
 * @param name   Name of collection to find 
 * 
 * @return Found collection or null
 */
TCollection* GetCollection(const TCollection* parent, const TString& name)
{
  TObject* o = GetObject(parent, name, TCollection::Class());
  if (!o) return 0;
  return static_cast<TCollection*>(o);
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
 */
void RerunELossFits(Bool_t forceSet=false, 
		    const TString& input="forward_eloss.root", 
		    const TString& output="")
{
  const char* fwd = "$ALICE_ROOT/PWGLF/FORWARD/analysis2";
  gROOT->Macro(Form("%s/scripts/LoadLibs.C", fwd));

  TFile*  inFile  = 0;
  TFile*  outFile = 0;
  TString outName(output);
  if (outName.IsNull()) {
    outName = input;
    outName.ReplaceAll(".root", "_rerun.root");
  }

  try {
    // --- Open input file ---------------------------------------------
    inFile = TFile::Open(input, "READ");
    if (!inFile)
      throw TString::Format("Failed to open %s", input.Data());

    // --- InFiled input collections --------------------------------------
    TCollection* inFwdSum = GetCollection(inFile, "Forward");
    if (!inFwdSum) throw TString("Cannot proceed without sums");

    TCollection* inFwdRes = GetCollection(inFile, "ForwardResults");
    if (!inFwdRes) throw TString("Cannot proceed with merged list");

    TCollection* inEFSum = GetCollection(inFwdRes, "fmdEnergyFitter");
    if (!inEFSum) throw TString("Cannot proceed without accumulated data");
    
    TCollection* inEFRes = GetCollection(inFwdRes, "fmdEnergyFitter");
    if (!inEFRes) throw TString("Cannot proceed without previous results");

    // --- Open output file --------------------------------------------
    outFile = TFile::Open(outName, "RECREATE");
    if (!outFile)
      throw TString::Format("Failed to open %s", outName.Data());

    // --- Write copy of sum collection to output --------------------
    TCollection* outFwdSum = static_cast<TCollection*>(inFwdSum->Clone());
    outFile->cd();
    outFwdSum->Write("Forward", TObject::kSingleKey);
    
    // --- Make our fitter object ------------------------------------
    AliFMDEnergyFitter* fitter = new AliFMDEnergyFitter("energy");
    if (forceSet || !fitter->ReadParameters(inEFSum)) {

      const TAxis* etaAxis = static_cast<TAxis*>(GetObject(inEFSum,"etaAxis"));
      if (!etaAxis) throw TString("Cannot proceed without eta axis");
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
    fitter->SetDebug(1);
    fitter->SetStoreResiduals(AliFMDEnergyFitter::kResidualSquareDifference);
    // fitter->SetRegularizationCut(3e6);
    // Set the number of bins to subtract from maximum of distributions
    // to get the lower bound of the fit range
    // fitter->SetFitRangeBinWidth(2);
    // Skip all of FMD2 and 3 
    // fitter->SetSkips(AliFMDEnergyFitter::kFMD2|AliFMDEnergyFitter::kFMD3);

    // --- Now do the fits -------------------------------------------
    fitter->Print();
    fitter->Fit(static_cast<TList*>(outFwdSum));
    
    // --- Copy full result folder -----------------------------------
    TCollection* outFwdRes = static_cast<TCollection*>(inFwdRes->Clone());
    // Remove old fits 
    TCollection* outEFRes = GetCollection(outFwdRes, "fmdEnergyFitter");
    outFwdRes->Remove(outEFRes);
    // Make our new fit results folder, and add it to results folder
    TCollection* tmp = GetCollection(outFwdSum, "fmdEnergyFitter");
    outEFRes = static_cast<TCollection*>(tmp->Clone());
    outEFRes->Add(new TNamed("refitted", "Refit of the data"));
    outFwdRes->Add(outEFRes);

    // --- Write out new results folder ------------------------------
    outFile->cd();
    outFwdRes->Write("ForwardResults", TObject::kSingleKey);
    Printf("Wrote results to \"%s\" (%s)", outName.Data(), outFile->GetName());
  }
  catch (const TString& e) {
    Error("RerunELossFits", e);
  }
  if (inFile)  inFile->Close();
  if (outFile) {
    Printf("Wrote new output to \"%s\"", outName.Data());
    outFile->Close();
  }
    
  gROOT->Macro(Form("%s/corrs/DrawCorrELoss.C(false,false,\"%s\")", 
		    fwd, outName.Data()));
}


