Color_t Color(UShort_t d, Char_t r ) const 
{ 
  return ((d == 1 ? kRed : (d == 2 ? kGreen : kBlue))
	  + ((r == 'I' || r == 'i') ? 2 : -2));
}

void
RunCopySecMap(UShort_t sys, UShort_t cms, Short_t field, const Char_t* path=0)
{
  RunCopySecMap(sys == 1 ? "pp" : "PbPb", 
		cms, 
		field, 
		path);
}
/** 
 * 
 * @param sys       Collision system 
 * @param cms       Center of mass energy per nucleon in GeV
 * @param field     Magnetic field 
 * 
 * @ingroup pwg2_forward_analysis_scripts
 */
void
RunCopySecMap(const char* sys, UShort_t cms, Short_t field,
	      const Char_t* path=0)
{
  gROOT->Macro("$ALICE_ROOT/PWG2/FORWARD/analysis2/scripts/LoadLibs.C");
  gSystem->Load("libPWG2forward.so");

  AliFMDAnaParameters* p = AliFMDAnaParameters::Instance();
  p->SetEnergy(Float_t(cms));
  p->SetMagField(Float_t(field));
  p->SetCollisionSystem(sys);
  if (path) {
    std::cout<<"Setting path to "<<path<<std::endl;
    p->SetBackgroundPath(path);
    p->SetEnergyPath(path);
    p->SetEventSelectionPath(path);
    p->SetSharingEfficiencyPath(path);
  }
  p->Init(true, AliFMDAnaParameters::kBackgroundCorrection);
 
  Int_t    nVtx   = p->GetNvtxBins();
  Double_t minVtx = -p->GetVtxCutZ();
  Double_t maxVtx = p->GetVtxCutZ();
  Int_t    nEta   = p->GetNetaBins();
  Double_t minEta = p->GetEtaMin();
  Double_t maxEta = p->GetEtaMax();

  TAxis vtxAxis(nVtx, minVtx, maxVtx);
  AliFMDCorrSecondaryMap* m = new AliFMDCorrSecondaryMap;
  m->SetVertexAxis(nVtx, minVtx, maxVtx);
  m->SetEtaAxis(nEta,minEta,maxEta);
  AliFMDCorrDoubleHit* dh = new AliFMDCorrDoubleHit;

  for (UShort_t d = 1; d <= 3; d++) { 
    UShort_t nQ = (d == 1 ? 1 : 2);
    for (UShort_t q = 0; q < nQ; q++) { 
      Char_t r = (q == 0 ? 'I' : 'O');
     
      TH1F* old2hit = p->GetDoubleHitCorrection(d, r);
      if (!old2hit) 
	Warning("RunCopySecMap",
		"Didn't find double hit correction for FMD%d%c", d, r);
      else { 
	TH1D* new2hit = new TH1D(Form("FMD%d%c", d, r), 
				 Form("Double hit correction for FMD%d%c",d,r), 
				 nEta, minEta, maxEta);
	new2hit->SetXTitle("#eta");
	new2hit->SetYTitle("#sum_{i} N_{i,strips hit}(#eta)/"
			   "#sum_{i} N_{i,total hits}(#eta)");
	new2hit->SetFillColor(Color(d,r));
	new2hit->SetFillStyle(3001);
	new2hit->SetDirectory(0);
	new2hit->SetStats(0);
	new2hit->Sumw2();
	new2hit->Add(old2hit);

	Info("RunCopySecMap", 
	     "Copying %s to %s", old2hit->GetName(), new2hit->GetName());

	dh->SetCorrection(d, r, new2hit);
      }

      for (UShort_t b = 1; b <= nVtx; b++) { 
	TH2F* oldmap = p->GetBackgroundCorrection(d, r, b-1);
	if (!oldmap) {
	  Warning("RunCopySecMap",
		  "Didn't find secondary map correction "
		  "for FMD%d%c, vertex bin %3d", d, r, b);
	  continue;
	}
       
	TH2D* newmap = new TH2D(Form("FMD%d%c_vtxbin%03d", d, r, b), 
				Form("Secondary map correction for FMD%d%c "
				     "in vertex bin %d [%+8.4f,%+8.4f]", 
				     d, r, b, vtxAxis.GetBinLowEdge(b), 
				     vtxAxis.GetBinUpEdge(b)), 
				nEta, minEta, maxEta, 
				oldmap->GetYaxis()->GetNbins(), 
				oldmap->GetYaxis()->GetXmin(), 
				oldmap->GetYaxis()->GetXmax());
	newmap->SetXTitle("#eta");
	newmap->SetYTitle("#phi [radians]");
	newmap->SetZTitle("#sum_{i} N_{ch,i,primary} / #sum_{i} N_{ch,i,FMD}");
	newmap->SetDirectory(0);
	newmap->SetStats(0);
	newmap->Sumw2();
	newmap->Add(oldmap);

	Info("RunCopySecMap",
	     "Copying %s to %s", oldmap->GetName(), newmap->GetName());

	m->SetCorrection(d, r, b, newmap);
      }
    }
  }
  UShort_t isys = AliForwardUtil::ParseCollisionSystem(sys);
  AliForwardCorrectionManager& mgr = AliForwardCorrectionManager::Instance();
  TString fname(mgr.GetFileName(AliForwardCorrectionManager::kSecondaryMap,
				isys, cms, field, false));
  TFile* output = TFile::Open(fname.Data(), "RECREATE");
  if (!output) { 
    Warning("Run", "Failed to open output file %s", fname.Data());
    return kFALSE;
  }
  m->Write(mgr.GetObjectName(AliForwardCorrectionManager::kSecondaryMap));
  output->Write();
  output->Close();
  Info("Run", "File %s created.  It should be copied to %s and stored in SVN",
       fname.Data(),mgr.GetFileDir(AliForwardCorrectionManager::kSecondaryMap));

  fname = mgr.GetFileName(AliForwardCorrectionManager::kDoubleHit,
			  isys, cms, field, false);
  output = TFile::Open(fname.Data(), "RECREATE");
  if (!output) { 
    Warning("Run", "Failed to open output file %s", fname.Data());
    return kFALSE;
  }
  dh->Write(mgr.GetObjectName(AliForwardCorrectionManager::kDoubleHit));
  output->Write();
  output->Close();
  Info("Run", "File %s created.  It should be copied to %s and stored in SVN",
       fname.Data(),mgr.GetFileDir(AliForwardCorrectionManager::kDoubleHit));
}
//
// EOF
//
