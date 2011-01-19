//_____________________________________________________________________
TH2* MakeOneRing(UShort_t d, Char_t r, Float_t vz, Int_t& nDead)
{
  AliFMDGeometry*   geom = AliFMDGeometry::Instance();
  AliFMDParameters* pars = AliFMDParameters::Instance();

  UShort_t nS = (r == 'I' || r == 'i' ?  20 :  40);
  UShort_t nT = (r == 'I' || r == 'i' ? 512 : 256);
  
  // Make our two histograms 
  TH2* hAll = new TH2D("all","All",200,-4,6,nS,0,2*TMath::Pi());
  hAll->SetXTitle("#eta");
  hAll->SetYTitle("#phi");
  hAll->Sumw2();
  hAll->SetDirectory(0);
  TH2* hOK  = static_cast<TH2*>(hAll->Clone());
  hOK->SetDirectory(0);
  
  // Loop over all sectors and strips in this ring 
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

      // Check if this is a dead channel or not 
      Bool_t isDead = pars->IsDead(d, r, s, t);
      hAll->Fill(eta, phi);
      if (!isDead) hOK->Fill(eta, phi);
      else         nDead++;
    }
  }
  // Divide out the efficiency. 
  hOK->Divide(hOK,hAll,1,1,"B");

  // Clean up
  delete hAll;

  // Return result 
  return hOK;
}

//_____________________________________________________________________
void MakeAcceptanceCorrection(Int_t   runNo, 
			      Int_t   nVtxBins=10, 
			      Float_t vtxLow=-10, 
			      Float_t vtxHigh=10){
  
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libPWG2forward2");
  
  // Float_t delta = (vtxHigh - vtxLow) / (Float_t)nVtxBins;
  
  //TGrid::Connect("alien://",0,0,"t");
  // --- Initialisations ------------------------------------------
  //Set up CDB manager
  AliCDBManager* cdb = AliCDBManager::Instance();
  //cdb->SetDefaultStorage("alien://Folder=/alice/data/2010/OCDB");
  cdb->SetDefaultStorage("local:///home/canute/ALICE/AliRoot/OCDB");
  cdb->SetRun(runNo);
  
  // Get an initialize parameters 
  AliFMDParameters* pars = AliFMDParameters::Instance();
  pars->Init();

  // Get an initialise geometry 
  AliFMDGeometry* geom = AliFMDGeometry::Instance();
  geom->Init();
  geom->InitTransforms();

  // --- Output object -----------------------------------------------
  // Make our correction object 
  AliFMDCorrAcceptance* corr = new AliFMDCorrAcceptance();
  corr->SetVertexAxis(nVtxBins, vtxLow, vtxHigh);

  // --- Loop over verticies and rings -------------------------------
  Int_t nDead = 0;
  Float_t dV = (vtxHigh - vtxLow) / nVtxBins;
  for (Float_t v = vtxLow+dV/2; v < vtxHigh; v += dV) { 
    for(UShort_t d = 1; d <= 3;d++) { 
      UShort_t nR = (d == 1 ? 1 : 2);
      for (UShort_t q = 0; q < nR; q++) { 
	Char_t   r  = (q == 0 ? 'I' : 'O');

	// Delegate to other function 
	TH2* ratio = MakeOneRing(d, r, v, nDead);
	if (!ratio) continue;

	// Set the correction 
	corr->SetCorrection(d, r, v, ratio);
      }
    }
  }

  // Write to a file 
  TFile* out = TFile::Open(Form("acceptance_%d.root", runNo), "RECREATE");
  corr->Write("acceptance");
  out->Write();
  out->Close();
}

//
// EOF
//

