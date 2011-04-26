//____________________________________________________________________
Bool_t
MakeSecMap(TList* list, Double_t low, Double_t high, 
	   AliFMDCorrSecondaryMap* m)
{
  // --- Get the list ------------------------------------------------
  TString lName(AliForwardMCCorrectionsTask::VtxBin::BinName(low, high));
  TList* vl = static_cast<TList*>(list->FindObject(lName));
  if (!vl) { 
    Error("MakeSecMap", "List %s not found in %s", 
	  lName.Data(), list->GetName());
    return false;
  }
  
  // --- Get primary distribution ------------------------------------
  const char* primaryName = "primary";
  TH2D* primary = static_cast<TH2D*>(vl->FindObject(primaryName));
  if (!primary) { 
    Error("MakeSecMap", "Couldn't not find histogram %s in %s", 
	  primaryName, lName.Data());
    return false;
  }
  TH2D* primaryI = static_cast<TH2D*>(primary->Clone("primaryI"));
  TH2D* primaryO = static_cast<TH2D*>(primary->Clone("primaryO"));
  primaryI->SetDirectory(0);
  primaryO->SetDirectory(0);
  primaryI->RebinY(2);


  // --- Calculate vertex --------------------------------------------
  Double_t vz = (high+low) / 2;

  // --- Loop over rings ---------------------------------------------
  for (UShort_t d = 1; d <= 3; d++) { 
    UShort_t nr = (d == 1 ? 1 : 2);
    for (UShort_t q = 0; q < nr; q++) { 
      Char_t r = (q == 0 ? 'I' : 'O');
      
      const char* ringName = Form("FMD%d%c_cache", d, r);
      TH2D* ring = static_cast<TH2D*>(vl->FindObject(ringName));
      if (!ring) { 
	Error("MakeSecMap", "Didn't find histogram %s in %s", 
	      ringName, vl->GetName());
	vl->ls();
	continue;
      }
      
      TH2D* tmp = static_cast<TH2D*>(ring->Clone("tmp"));
      TString tmpName = tmp->GetName();
      tmpName.ReplaceAll("_cache", "");
      tmp->SetName(tmpName);
      tmp->SetDirectory(0);
      tmp->Divide(q == 0 ? primaryI : primaryO);


      m->SetCorrection(d, r, vz, tmp);
    }
  }
  return true;
}

//____________________________________________________________________
void
MakeCorrSecMap(const char* filename, 
	       const char* sys="pp", 
	       UShort_t    cms=900, 
	       Short_t     field=+5)
{
  // --- Load code ---------------------------------------------------
  // gROOT->Macro("$ALICE_ROOT/PWG2/FORWARD/analysis2/scripts/LoadLibs.C");


  // --- Get the file ------------------------------------------------
  TFile* file = TFile::Open(filename, "READ");
  if (!file) {
    Error("MakeCorrSecMap", "Couldn't open file %s", filename);
    return;
  }

  // --- Get the parent list -----------------------------------------
  const char* forwardName = "ForwardSums";
  TList* forward = static_cast<TList*>(file->Get(forwardName));
  if (!forward) { 
    Error("MakeCorrSecMap", "Couldn't get list %s from %s", 
	  forwardName, filename);
    return;
  }

  // --- Get the vertex axis -----------------------------------------
  const char* vtxName = "vtxAxis";
  TH1* vtxHist = static_cast<TH1*>(forward->FindObject(vtxName));
  if (!vtxHist) {
    Error("MakeCorrSecMap", "Couldn't get histogram %s from %s", 
	  vtxName, forwardName);
    return;
  }
  const TAxis& vtxAxis = *(vtxHist->GetXaxis());

  // --- Get the vertex axis -----------------------------------------
  const char* etaName = "etaAxis";
  TH1* etaHist = static_cast<TH1*>(forward->FindObject(etaName));
  if (!etaHist) {
    Error("MakeCorrSecMap", "Couldn't get histogram %s from %s", 
	  etaName, forwardName);
    return;
  } 
  const TAxis& etaAxis = *(etaHist->GetXaxis());

  // --- Fill correction object --------------------------------------
  AliFMDCorrSecondaryMap* corr = new AliFMDCorrSecondaryMap;
  corr->SetVertexAxis(vtxAxis);
  corr->SetEtaAxis(etaAxis);

  for (Int_t i = 1; i <= vtxAxis.GetNbins(); i++) { 
    Double_t low  = vtxAxis.GetBinLowEdge(i);
    Double_t high = vtxAxis.GetBinUpEdge(i);

    MakeSecMap(forward, low, high, corr);
  }

  // --- Get output filename and open --------------------------------
  UShort_t isys = AliForwardUtil::ParseCollisionSystem(sys);
  AliForwardCorrectionManager& mgr = AliForwardCorrectionManager::Instance();
  TString fname(mgr.GetFileName(AliForwardCorrectionManager::kSecondaryMap,
				isys, cms, field, false));
  TFile* output = TFile::Open(fname.Data(), "RECREATE");
  if (!output) { 
    Warning("Run", "Failed to open output file %s", fname.Data());
    return kFALSE;
  }

  // --- Write to output ---------------------------------------------
  corr->Write(mgr.GetObjectName(AliForwardCorrectionManager::kSecondaryMap));
  output->Write();
  output->Close();
  Info("Run", "File %s created.  It should be copied to %s and stored in SVN",
       fname.Data(),mgr.GetFileDir(AliForwardCorrectionManager::kSecondaryMap));

}
//____________________________________________________________________
//
// EOF
//
   
