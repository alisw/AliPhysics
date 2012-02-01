//_____________________________________________________________________
/** 
 * 
 * 
 * @param d 
 * @param r 
 * @param vz 
 * @param nDead 
 * 
 * @return 
 *
 * @ingroup pwg2_forward_scripts_corr
 */
TH2D* MakeOneRing(UShort_t d, Char_t r, Double_t vz, Int_t& nDead)
{
  AliFMDGeometry*   geom = AliFMDGeometry::Instance();
  AliFMDParameters* pars = AliFMDParameters::Instance();

  UShort_t nS = (r == 'I' || r == 'i' ?  20 :  40);
  UShort_t nT = (r == 'I' || r == 'i' ? 512 : 256);
  
  // Make our two histograms 
  TH2D* hAll = new TH2D("all","All",200,-4,6,nS,0,2*TMath::Pi());
  hAll->SetXTitle("#eta");
  hAll->SetYTitle("#phi");
  hAll->Sumw2();
  hAll->SetDirectory(0);
  TH2D* hOK  = static_cast<TH2D*>(hAll->Clone());
  hOK->SetDirectory(0);
  
  // Loop over all sectors and strips in this ring 
  Int_t nOK  = 0;
  Int_t nAll = 0;
  for (UShort_t s = 0; s < nS; s++) { 
    for (UShort_t t = 0; t < nT; t++) { 
      // Get eta,phi by quering the geometry (first for (x,y,z), then
      // correcting for the vertex position, and then calculating 
      // (eta, phi))
      Double_t x, y, z;
      geom->Detector2XYZ(d, r, s, t, x, y, z);
      z -= vz;
      Double_t q, eta, phi, theta;
      AliFMDGeometry::XYZ2REtaPhiTheta(x, y, z, q, eta, phi, theta);
      if (phi < 0) phi += 2*TMath::Pi();

      // Check if this is a dead channel or not 
      Bool_t isDead = pars->IsDead(d, r, s, t);
      hAll->Fill(eta, phi);
      nAll++;
      if (!isDead) {
	hOK->Fill(eta, phi);
	nOK++;
      }
      else         nDead++;
    }
  }
  // Divide out the efficiency. 
  hOK->Divide(hOK,hAll,1,1,"B");

  // Clean up
  delete hAll;

  Info("MakeAcceptanceCorrections","Made correction for FMD%d%c at vz=%f - "
       "%d strips out of %d OK", d, r, vz, nOK, nAll);

  // Return result 
  return hOK;
}

//_____________________________________________________________________
/** 
 * 
 * 
 * @param runNo 
 * @param system 
 * @param energy 
 * @param field 
 * @param nVtxBins 
 * @param vtxLow 
 * @param vtxHigh 
 *
 * @ingroup pwg2_forward_scripts_corr
 */
void MakeAcceptanceCorrection(Int_t   runNo=121526, 
			      Int_t   system = 1,
			      Float_t energy = 900,
			      Float_t field  = 5,
			      Int_t   nVtxBins=10, 
			      Float_t vtxLow=-10, 
			      Float_t vtxHigh=10){
  
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libPWG2forward2");
  
  // Float_t delta = (vtxHigh - vtxLow) / (Float_t)nVtxBins;
  
  Bool_t kGridOnline = kTRUE; 
  if(!(TGrid::Connect("alien://",0,0,"t")))
    kGridOnline = kFALSE;
  
  // --- Initialisations ------------------------------------------
  //Set up CDB manager
  Info("MakeAcceptanceCorrections","Setting up OCDB");
  
  AliCDBManager* cdb = AliCDBManager::Instance();
  if(kGridOnline)
    cdb->SetDefaultStorage("alien://Folder=/alice/data/2010/OCDB");
  else
    cdb->SetDefaultStorage("local://$(ALICE_ROOT)/OCDB");
  cdb->SetRun(runNo);
  
  // Get the geometry 
  Info("MakeAcceptanceCorrections","Loading geometry");
  AliGeomManager::LoadGeometry();

  // Get an initialize parameters 
  Info("MakeAcceptanceCorrections","Intialising parameters");
  AliFMDParameters* pars = AliFMDParameters::Instance();
  pars->Init();

  // Get an initialise geometry 
  Info("MakeAcceptanceCorrections","Initialising geomtry");
  AliFMDGeometry* geom = AliFMDGeometry::Instance();
  geom->Init();
  geom->InitTransformations();

  // --- Output object -----------------------------------------------
  // Make our correction object 
  AliFMDCorrAcceptance* corr = new AliFMDCorrAcceptance();
  corr->SetVertexAxis(nVtxBins, vtxLow, vtxHigh);

  // --- Loop over verticies and rings -------------------------------
  Int_t nDead = 0;
  Float_t dV = (vtxHigh - vtxLow) / nVtxBins;
  for (Double_t v = vtxLow+dV/2; v < vtxHigh; v += dV) { 
    for(UShort_t d = 1; d <= 3;d++) { 
      UShort_t nR = (d == 1 ? 1 : 2);
      for (UShort_t q = 0; q < nR; q++) { 
	Char_t   r  = (q == 0 ? 'I' : 'O');

	// Delegate to other function 
	TH2D* ratio = MakeOneRing(d, r, v, nDead);
	if (!ratio) continue;

	// Set the correction 
	corr->SetCorrection(d, r, v, ratio);
      }
    }
  }

  // Write to a file 
#if 0
  TFile* out = TFile::Open(Form("acceptance_%d.root", runNo), "RECREATE");
  corr->Write("acceptance");
  out->Write();
  out->Close();
#else 
  Info("MakeAcceptanceCorrections","Writing to disk");
  AliForwardCorrectionManager& cm = AliForwardCorrectionManager::Instance();
  TString fname = cm.GetFileName(AliForwardCorrectionManager::kAcceptance, 
				 system, energy, field, false);
  TFile* out = TFile::Open(fname.Data(), "RECREATE");
  corr->Write(cm.GetObjectName(AliForwardCorrectionManager::kAcceptance));
  out->Write();
  out->Close();
  Info("MakeAcceptanceCorrections","File %s generated.  "
       "It should be copied to %s", fname.Data(), 
       cm.GetFileDir(AliForwardCorrectionManager::kAcceptance));
  out = TFile::Open(fname.Data(), "READ");
  new TBrowser;
#endif
  
}

//
// EOF
//

