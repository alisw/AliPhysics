/** 
 * 
 * 
 * @param d 
 * @param r 
 * 
 * @return 
 * 
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
 * 
 * @ingroup pwg2_forward_scripts_corr
 */
void
RunCopyVtxBias(UShort_t sys, UShort_t cms, Short_t field, const char* path=0)
{
  RunCopyVtxBias(sys == 1 ? "pp" : "PbPb", cms, field, path);
}

/** 
 * 
 * @param sys       Collision system 
 * @param cms       Center of mass energy per nucleon in GeV
 * @param field     Magnetic field 
 * @param path      File path
 * 
 * 
 * @ingroup pwg2_forward_scripts_corr
 */
void
RunCopyVtxBias(const char* sys, UShort_t cms, Short_t field, const char* path=0)
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
	  AliFMDAnaParameters::kEventSelectionEfficiency);
 
  Int_t    nVtx   = p->GetNvtxBins();
  Double_t minVtx = -p->GetVtxCutZ();
  Double_t maxVtx = -p->GetVtxCutZ();
  Int_t    nEta   = p->GetNetaBins();
  Double_t minEta = p->GetEtaMin();
  Double_t maxEta = p->GetEtaMax();

  TAxis vtxAxis(nVtx, minVtx, maxVtx);
  AliFMDCorrVertexBias* obj = new AliFMDCorrVertexBias;
  obj->SetVertexAxis(vtxAxis);

  for (UShort_t b = 1; b <= nVtx; b++) { 
    for (UShort_t q = 0; q < 2; q++) { 
      Char_t r = q == 0 ? 'I' : 'O';

      TH2F* oldcorr = p->GetEventSelectionEfficiency("INEL",b-1,r);
      if (!oldcorr) {
	Warning("RunCopyVtxBias",
		"Didn't find secondary map correction "
		"for ring type %c, vertex bin %3d", r, b);
	continue;
      }
   
      TH2D* newcorr = new TH2D(Form("FMDX%c", r), 
			       Form("Vertex bias correction for %c rings "
				    "in vertex bin %d [%+8.4f,%+8.4f]", 
				    r, b, vtxAxis.GetBinLowEdge(b), 
				    vtxAxis.GetBinUpEdge(b)), 
			       nEta, minEta, maxEta, 
			       oldcorr->GetYaxis()->GetNbins(), 
			       oldcorr->GetYaxis()->GetXmin(), 
			       oldcorr->GetYaxis()->GetXmax());
      newcorr->SetXTitle("#eta");
      newcorr->SetYTitle("#phi [radians]");
      newcorr->SetZTitle("1/N_{t}#sum_{i}^{N_{tv}} N_{ch,i,primary} / "
			 "1/N_{v}#sum_{i}^{N_{v}} N_{ch,i,primary}");
      newcorr->SetDirectory(0);
      newcorr->SetStats(0);
      newcorr->Sumw2();
      newcorr->Add(oldcorr);

      Info("RunCopyVtxBias",
	   "Copying %s to %s", oldcorr->GetName(), newcorr->GetName());

      obj->SetCorrection(r, b, newcorr);
    }
  }
  UShort_t isys = AliForwardUtil::ParseCollisionSystem(sys);
  AliForwardCorrectionManager& mgr = AliForwardCorrectionManager::Instance();
  TString fname(mgr.GetFileName(AliForwardCorrectionManager::kVertexBias,
				isys, cms, field, false));
  TFile* output = TFile::Open(fname.Data(), "RECREATE");
  if (!output) { 
    Warning("Run", "Failed to open output file %s", fname.Data());
    return kFALSE;
  }
  obj->Write(mgr.GetObjectName(AliForwardCorrectionManager::kVertexBias));
  output->Write();
  output->Close();
  Info("Run", "File %s created.  It should be copied to %s and stored in SVN",
       fname.Data(),mgr.GetFileDir(AliForwardCorrectionManager::kVertexBias));

}
//
// EOF
//
