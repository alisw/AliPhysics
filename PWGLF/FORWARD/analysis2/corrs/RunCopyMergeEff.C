/** 
 * 
 * 
 * @param d 
 * @param r 
 * 
 * @return 
 * @ingroup pwg2_forward_scripts_corr
 */
Color_t Color(UShort_t d, Char_t r ) const 
{ 
  return ((d == 1 ? kRed : (d == 2 ? kGreen : kBlue))
	  + ((r == 'I' || r == 'i') ? 2 : -2));
}

/** 
 * 
 * 
 * @param sys 
 * @param cms 
 * @param field 
 * @param path 
 * @ingroup pwg2_forward_scripts_corr
 */
void
RunCopyMergeEff(UShort_t sys, UShort_t cms, Short_t field, const Char_t* path=0)
{
  RunCopyMergeEff(sys == 1 ? "pp" : "PbPb", 
		  cms, 
		  field, 
		  path);
}

/** 
 * 
 * @param sys       Collision system 
 * @param cms       Center of mass energy per nucleon in GeV
 * @param field     Magnetic field 
 * @param path      File path
 * 
 * @ingroup pwg2_forward_scripts_corr
 */
void
RunCopyMergeEff(const char* sys, UShort_t cms, 
		Short_t field, const char* path=0)
{
  gROOT->Macro("$ALICE_ROOT/PWG2/FORWARD/analysis2/scripts/LoadLibs.C");
  gSystem->Load("libPWG2forward.so");

  AliFMDAnaParameters* p = AliFMDAnaParameters::Instance();
  p->SetEnergy(Float_t(cms));
  p->SetMagField(Float_t(field));
  p->SetCollisionSystem(sys);
  if (path) {
    p->SetBackgroundPath(path);
    p->SetEnergyPath(path);
    p->SetEventSelectionPath(path);
    p->SetSharingEfficiencyPath(path);
  }
  p->Init(true, AliFMDAnaParameters::kBackgroundCorrection|
	  AliFMDAnaParameters::kSharingEfficiency);
  
  Int_t    nVtx   = p->GetNvtxBins();
  Double_t minVtx = -p->GetVtxCutZ();
  Double_t maxVtx = -p->GetVtxCutZ();
  Int_t    nEta   = p->GetNetaBins();
  Double_t minEta = p->GetEtaMin();
  Double_t maxEta = p->GetEtaMax();
  
  TAxis vtxAxis(nVtx, minVtx, maxVtx);
  AliFMDCorrMergingEfficiency* m = new AliFMDCorrMergingEfficiency;
  m->SetVertexAxis(nVtx, minVtx, maxVtx);
  
  for (UShort_t d = 1; d <= 3; d++) { 
    UShort_t nQ = (d == 1 ? 1 : 2);
    for (UShort_t q = 0; q < nQ; q++) { 
      Char_t r = (q == 0 ? 'I' : 'O');
      
      for (UShort_t b = 1; b <= nVtx; b++) { 
	TH1F* oldmap = p->GetSharingEfficiency(d, r, b-1);
	if (!oldmap) {
	  Warning("RunCopyMergeEff",
		  "Didn't find secondary map correction "
		  "for FMD%d%c, vertex bin %3d", d, r, b);
	  continue;
	}
	
	TH1D* newmap = new TH1D(Form("FMD%d%c_vtxbin%03d", d, r, b), 
				Form("Merging efficiency for FMD%d%c "
				     "in vertex bin %d [%+8.4f,%+8.4f]", 
				     d, r, b, vtxAxis.GetBinLowEdge(b), 
				     vtxAxis.GetBinUpEdge(b)), 
				nEta, minEta, maxEta);
	newmap->SetXTitle("#eta");
	newmap->SetYTitle("dN_{ch,i,incl}/d#eta / #sum_{i} N_{ch,i,FMD}");
	newmap->SetDirectory(0);
	newmap->SetStats(0);
	newmap->Sumw2();
	newmap->Add(oldmap);
	
	Info("RunCopyMergeEff",
	     "Copying %s to %s", oldmap->GetName(), newmap->GetName());
	
	m->SetCorrection(d, r, b, newmap);
      }
    }
  }
  UShort_t isys = AliForwardUtil::ParseCollisionSystem(sys);
  AliForwardCorrectionManager& mgr = AliForwardCorrectionManager::Instance();
  TString fname(mgr.GetFileName(AliForwardCorrectionManager::kMergingEfficiency,
				isys, cms, field, false));
  TFile* output = TFile::Open(fname.Data(), "RECREATE");
  if (!output) { 
    Warning("Run", "Failed to open output file %s", fname.Data());
    return kFALSE;
  }
  m->Write(mgr.GetObjectName(AliForwardCorrectionManager::kMergingEfficiency));
  output->Write();
  output->Close();
  Info("Run", "File %s created.  It should be copied to %s and stored in SVN",
       fname.Data(),
       mgr.GetFileDir(AliForwardCorrectionManager::kMergingEfficiency));
  
}
//
// EOF
//
