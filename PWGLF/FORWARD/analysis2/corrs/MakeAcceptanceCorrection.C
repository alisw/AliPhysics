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
 * @ingroup pwglf_forward_scripts_corr
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
  Int_t nPhi = hAll->GetNbinsY();
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

      // Special check for FMD2i
      Int_t VA = t/128;
      if(d==2 && r=='I' && VA>1 && (s==16 || s==17)) isDead =true;

      // Find the eta bin number and corresponding overflow bin
      Int_t etaBin = hAll->GetXaxis()->FindBin(eta);
      Int_t ovrBin = hAll->GetBin(etaBin, nPhi+1); 
      
      // Increment all histogram 
      hAll->Fill(eta, phi);
      hAll->AddBinContent(ovrBin);
      nAll++;

      // If not dead, increment OK histogram 
      if (!isDead) {
	hOK->Fill(eta, phi);
	hOK->AddBinContent(ovrBin);
	nOK++;
      }
      else         nDead++;
    }
  }
  // Divide out the efficiency. 
  // Note, that the overflow bins along eta now contains the ratio 
  // nOK/nAll Strips for a given eta bin. 
  hOK->Divide(hOK,hAll,1,1,"B");

  // Clean up
  delete hAll;

  Info("ExtractAcceptances","Made correction for FMD%d%c at vz=%f - "
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
 * @ingroup pwglf_forward_scripts_corr
 */
void ExtractAcceptance(Int_t   runNo=121526, 
		       Int_t   system = 1,
		       Float_t energy = 900,
		       Float_t field  = 5,
		       Int_t   nVtxBins=10, 
		       Float_t vtxLow=-10, 
		       Float_t vtxHigh=10)
{  
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libPWGLFforward2");
  
  // Float_t delta = (vtxHigh - vtxLow) / (Float_t)nVtxBins;
  
  Bool_t kGridOnline = kTRUE; 
  if(!(TGrid::Connect("alien://",0,0,"t")))
    kGridOnline = kFALSE;
  
  // --- Initialisations ------------------------------------------
  //Set up CDB manager
  Info("ExtractAcceptances","Setting up OCDB");
  
  AliCDBManager* cdb = AliCDBManager::Instance();
  if(kGridOnline)
    cdb->SetDefaultStorage("alien://Folder=/alice/data/2012/OCDB");
  else
    cdb->SetDefaultStorage("local://$(ALICE_PHYSICS)/OCDB");
  cdb->SetRun(runNo);
  
  // Get the geometry 
  Info("ExtractAcceptances","Loading geometry");
  AliGeomManager::LoadGeometry();

  // Get an initialize parameters 
  Info("ExtractAcceptances","Intialising parameters");
  AliFMDParameters* pars = AliFMDParameters::Instance();
  pars->Init();

  // Get an initialise geometry 
  Info("ExtractAcceptances","Initialising geomtry");
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
  Info("ExtractAcceptances","Writing to disk");
  AliForwardCorrectionManager& cm = AliForwardCorrectionManager::Instance();
  TString fname = cm.GetFileName(AliForwardCorrectionManager::kAcceptance, 
				 system, energy, field, false);
  TFile* out = TFile::Open(fname.Data(), "RECREATE");
  corr->SetHasOverflow();
  corr->Write(cm.GetObjectName(AliForwardCorrectionManager::kAcceptance));
  out->Write();
  out->Close();

  std::ofstream f("Upload.C");
  if (!f) { 
    Error("ExtractELoss", "Failed to open Upload.C");
    return;
  }
  f << "// Generated by ExtractAcceptance.C\n"
    << "void Upload(const TUrl& url)\n"
    << "{\n"
    << "  if (TString(\"alien\").EqualTo(url.GetProtocol())) {\n"
    << "    if (!TGrid::Connect(\"alien://\")) {\n"
    << "      Error(\"Upload\", \"Failed to connect to AliEn\");\n"
    << "      return;\n"
    << "    }\n"
    << "  }\n\n";

  mgr.SetPrefix("");
  TString fef(mgr.GetFileName(AliForwardCorrectionManager::kAcceptance, 
			      sys, sNN, field, false));
  TString fep(mgr.GetFilePath(AliForwardCorrectionManager::kAcceptance, 
			      sys, sNN, field, false));
  f << "  TString src  = \"" << fef << "\";\n"
    << "  TString dest = \"" << fep << "\";\n"
    << "  TString out; out.Form(\"%s%s\",url.GetUrl(),dest.Data());\n\n"
    << "  TString dir(gSystem->DirName(out));\n"
    << "  if (gSystem->AccessPathName(dir)) {\n"
    << "    if (gSystem->mkdir(dir, true) < 0) {\n"
    << "      Warning(\"Upload\",\"Failed to make directory %s\","
    << "              dir.Data());\n"
    << "      return;\n"
    << "    }\n"
    << "  }\n"
    << "  if (!TFile::Cp(src,out)) \n"
    << "    Warning(\"Upload\",\"Failed to upload %s -> %s\",\n"
    << "            src.Data(), out.Data());\n"
    << "}\n"
    << "// EOF"
    << std::endl;
  f.close();
  
  Info("ExtracAcceptance", 
       "Run generated Upload.C(DEST) script to copy files in place");
}

//
// EOF
//

