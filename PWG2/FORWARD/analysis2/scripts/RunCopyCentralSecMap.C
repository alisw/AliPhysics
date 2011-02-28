Color_t Color(UShort_t d, Char_t r ) const 
{ 
  return ((d == 1 ? kRed : (d == 2 ? kGreen : kBlue))
	  + ((r == 'I' || r == 'i') ? 2 : -2));
}

/*void
RunCopyCentralSecMap(UShort_t sys, UShort_t cms, Short_t field, const Char_t* path=0)
{
  RunCopyCentralSecMap(sys == 1 ? "pp" : "PbPb", 
		cms, 
		field, 
		path);
		}*/
/** 
 * 
 * @param sys       Collision system 
 * @param cms       Center of mass energy per nucleon in GeV
 * @param field     Magnetic field 
 * 
 * @ingroup pwg2_forward_analysis_scripts
 */
//void
//RunCopyCentralSecMap(const char* sys, UShort_t cms, Short_t field,
		     //	      const Char_t* path=0)
void
RunCopyCentralSecMap(UShort_t sys, UShort_t cms, Short_t field, const Char_t* path=0)
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
  p->Init(true, AliFMDAnaParameters::kBackgroundCorrection);
  
  Int_t    nVtx   = p->GetNvtxBins();
  Double_t minVtx = -p->GetVtxCutZ();
  Double_t maxVtx = p->GetVtxCutZ();
  Int_t    nEta   = p->GetNetaBins();
  Double_t minEta = p->GetEtaMin();
  Double_t maxEta = p->GetEtaMax();

  TAxis vtxAxis(nVtx, minVtx, maxVtx);
  AliCentralCorrSecondaryMap* m = new AliCentralCorrSecondaryMap;
  m->SetVertexAxis(nVtx, minVtx, maxVtx);
  
  AliCentralCorrAcceptance* a = new AliCentralCorrAcceptance;
  a->SetVertexAxis(nVtx, minVtx, maxVtx);
  
  
  //m->SetEtaAxis(nEta,minEta,maxEta);
  //AliFMDCorrDoubleHit* dh = new AliFMDCorrDoubleHit;

  

  for (UShort_t b = 1; b <= nVtx; b++) { 
    TH2F* oldmap = p->GetBackgroundCorrection(0, 'Q', b-1);
    if (!oldmap) {
      Warning("RunCopySecMap",
	      "Didn't find secondary map correction "
	      "for SPD, vertex bin %3d", b);
      continue;
    }
    
	TH2D* newmap = new TH2D(Form("SPD_vtxbin%03d",  b), 
				Form("Secondary map correction for SPD "
				     "in vertex bin %d [%+8.4f,%+8.4f]", 
				     b, vtxAxis.GetBinLowEdge(b), 
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
	
	m->SetCorrection(b, newmap);
	std::cout<<m->GetCorrection(b)<<std::endl;
	
	//Acceptance 
	TH1F* oldacc = p->GetSPDDeadCorrection(b-1);
	if (!oldacc) {
	  Warning("RunCopySecMap",
		  "Didn't find acceptance correction "
		  "for SPD, vertex bin %3d", b);
	  continue;
	}
	
	TH1D* newacc = new TH1D(Form("SPDdead_vtxbin%03d",  b), 
				Form("Acceptance correction for SPD "
				     "in vertex bin %d [%+8.4f,%+8.4f]", 
				     b, vtxAxis.GetBinLowEdge(b), 
				     vtxAxis.GetBinUpEdge(b)), 
				nEta, minEta, maxEta);
	newacc->SetXTitle("#eta");
	newacc->SetYTitle("correction");
	//newacc->SetZTitle("#sum_{i} N_{ch,i,primary} / #sum_{i} N_{ch,i,FMD}");
	newacc->SetDirectory(0);
	newacc->SetStats(0);
	newacc->Sumw2();
	newacc->Add(oldacc);
	
	Info("RunCopySecMap",
	     "Copying %s to %s", oldacc->GetName(), newacc->GetName());
	
	a->SetCorrection(b, newacc);
	std::cout<<a->GetCorrection(b)<<std::endl;
	

  }
  
  AliCentralMultiplicityTask::Manager* task = new AliCentralMultiplicityTask::Manager();

  TFile f(task->GetFullFileName(0.,sys+1,cms,field), "RECREATE");
  m->Write(task->GetSecMapName());
  f.Close();
  TFile f2(task->GetFullFileName(1.,sys+1,cms,field), "RECREATE");
  a->Write(task->GetAcceptanceName());
  f2.Close();

}
//
// EOF
//
