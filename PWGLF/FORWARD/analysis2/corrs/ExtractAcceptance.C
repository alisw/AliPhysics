//_____________________________________________________________________
/** 
 * 
 * 
 * @param d      Detector
 * @param r      Ring 
 * @param vz     Z--coordinate of interaction point
 * @param nDead  On returm the number of dead strips
 * @param deadScript  Output stream for dead strips 
 * 
 * @return 
 *
 * @ingroup pwglf_forward_scripts_corr
 */
TH2D* MakeOneRing(UShort_t      d, 
		  Char_t        r, 
		  Double_t      vz, 
		  Int_t&        nDead, 
		  std::ostream* deadScript)
{
  AliFMDGeometry*   geom  = AliFMDGeometry::Instance();
  AliFMDParameters* pars  = AliFMDParameters::Instance();
  Float_t           konst = pars->GetDACPerMIP();

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
  
  if (deadScript)
    *deadScript << "\n  // === FMD" << d << r << " === " << std::endl;
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
      TString reason;
      
      // Check if this is a dead channel or not 
      Bool_t isDead = pars->IsDead(d, r, s, t);
      if (isDead) {
	if (deadScript)
	  Warning("", "FMD%d%c[%2d,%3d] is dead because OCDB says so");	
	reason = "OCDB";
      }
      
      // Pedestal, noise, gain 
      Double_t ped   = pars->GetPedestal     (d, r, s, t);
      Double_t noise = pars->GetPedestalWidth(d, r, s, t);
      Double_t gain  = pars->GetPulseGain    (d, r, s, t);
      if (ped < 10) {
	if (!isDead && deadScript) {
	  reason = "Low pedestal";
	  Warning("", "FMD%d%c[%2d,%3d] is dead because pedestal %f < 10",
		  d, r, s, t, ped);
	}
	isDead = true;
      }
      Float_t corr  = 0;
      if (noise > .5 && gain > .5) corr = noise / (gain * konst);
      if (corr > 0.05 || corr <= 0) {
	if (!isDead && deadScript) {
	  reason = "high noise/low gain";
	  Warning("", "FMD%d%c[%2d,%3d] is dead because %f/(%f*%f)=%f > 0.05 or negative",
		  d, r, s, t, noise, gain, konst, corr);
	}
	isDead = true;
      }
      
      
      

      // Special check for FMD2i - upper part of sectors 16/17 have 
      // have anomalous gains and/or noise - common sources are 
      // power regulartors for bias currents and the like
#if 0
      Int_t VA = t/128;
      if(d==2 && r=='I' && VA>1 && (s==16 || s==17)) {
	if (!isDead && deadScript) {
	  reason = "I say so";
	  Warning("", "FMD%d%c[%2d,%3d] is dead because I say so", d, r, s, t);
	}
	isDead =true;
      }
#endif
      
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
      else {
	nDead++;
	if (deadScript)
	  *deadScript << "  filter->AddDead(" << d << ",'" << r << "'," 
		      << std::setw(2) << s << ','
		      << std::setw(3) << t << "); // "
		      << reason << std::endl;
      }
    }
  }
  // Divide out the efficiency. 
  // Note, that the overflow bins along eta now contains the ratio 
  // nOK/nAll Strips for a given eta bin. 
  hOK->Divide(hOK,hAll,1,1,"B");

  // Invert overflow bin
  for (Int_t etaBin = 1; etaBin <= hOK->GetNbinsX(); etaBin++) { 
    Double_t ovr    = hOK->GetBinContent(etaBin, nPhi+1);
    Double_t novr   = (ovr < 1e-12 ? 0 : 1./ovr);
    hOK->SetBinContent(etaBin, nPhi+1, novr);
#if 0
    if (ovr > 0 && ovr != 1)
      Info("", "Setting overflow bin (%3d,%3d) to 1/%f=%f", etaBin, nPhi+1, 
	   ovr, hOK->GetBinContent(etaBin, nPhi+1));
#endif
  }
  // Clean up
  delete hAll;

  Printf("=== FMD%d%c at vz=%+5.1f - %d/%d=%3d%% OK (w/overflow)", 
	 d, r, vz, nOK, nAll, (100*nOK)/nAll);

  // Return result 
  return hOK;
}

//_____________________________________________________________________
/** 
 * 
 * 
 * @param runNo     Run number
 * @param nVtxBins  Number of @f$IP_{z}@f$ bins
 * @param vtxLow    Least @f$IP_{z}@f$
 * @param vtxHigh   Largest @f$IP_{z}@f$
 *
 * @ingroup pwglf_forward_scripts_corr
 */
void ExtractAcceptance(Int_t   runNo=121526, 
		       Int_t   nVtxBins=10, 
		       Float_t vtxLow=-10, 
		       Float_t vtxHigh=10)
{  
  const char* fwd = "$ALICE_PHYSICS/PWGLF/FORWARD/analysis2";
  gSystem->AddIncludePath(Form("-I%s", fwd));
  gROOT->Macro(Form("%s/scripts/LoadLibs.C", fwd));

  // gSystem->Load("libANALYSIS");
  // gSystem->Load("libANALYSISalice");
  // gSystem->Load("libPWGLFforward2");
  // Float_t delta = (vtxHigh - vtxLow) / (Float_t)nVtxBins;
  
  Bool_t gridOnline = kTRUE; 
  if(!(TGrid::Connect("alien://",0,0,"t")))
    gridOnline = kFALSE;

  // --- Figure out the year --------------------------------------
  UShort_t year = 0;
  if      (runNo <= 99999)  year = 2009;
  else if (runNo <= 139667) year = 2010;
  else if (runNo <= 170718) year = 2011;
  else if (runNo <= 194306) year = 2012;
  else if (runNo <= 197709) year = 2013;
  else if (runNo <= 208364) year = 2014;
  else                      year = 2015;
  if (year <= 0) { 
    Warning("", "Couldn't deduce the year from the run number");
    // return;
  }

  // --- Initialisations ------------------------------------------
  //Set up CDB manager
  Printf("=== Setting up OCDB");
  
  AliCDBManager* cdb = AliCDBManager::Instance();
  if(gridOnline) {
    cdb->SetDefaultStorageFromRun(runNo);
    // cdb->SetDefaultStorage(Form("alien://Folder=/alice/data/%4d/OCDB", year));
  }
  else
    cdb->SetDefaultStorage("local://$(ALICE_PHYSICS)/OCDB");
  cdb->SetRun(runNo);
  
  // Get the geometry 
  Printf("=== Loading geometry");
  AliGeomManager::LoadGeometry();

  // Get an initialize parameters 
  Printf("=== Intialising parameters");
  AliFMDParameters* pars = AliFMDParameters::Instance();
  pars->Init();

  // Get an initialise geometry 
  Printf("=== Initialising geomtry");
  AliFMDGeometry* geom = AliFMDGeometry::Instance();
  geom->Init();
  geom->InitTransformations();

  // --- Get the general run parameters ------------------------------
  AliCDBEntry* grpE = cdb->Get("GRP/GRP/Data");
  if (!grpE) { 
    AliWarningF("No GRP entry found for run %d", runNo);
    return;
  }
  AliGRPObject* grp = static_cast<AliGRPObject*>(grpE->GetObject());
  if (!grp) { 
    AliWarningF("No GRP object found for run %d", runNo);
    return;
  }
  Float_t  beamE = grp->GetBeamEnergy();
  TString  beamT = grp->GetBeamType();
# if 0 
  // This isn't really needed as the acceptance map is indifferent to
  // the field settings.
  Float_t  l3cur = grp->GetL3Current(AliGRPObject::kMean);
  Char_t   l3pol = grp->GetL3Polarity();
  Bool_t   l3lhc = grp->IsPolarityConventionLHC();
  Bool_t   l3uni = grp->IsUniformBMap();
  AliMagF* fldM  = 
    AliMagF::CreateFieldMap(TMath::Abs(l3cur) * (l3pol ? -1:1), 0, 
			    (l3lhc ? 0 : 1), l3uni, beamE, beamT.Data());
  Float_t  l3fld = fldM->SolenoidField();
#endif
  
  UShort_t sys = AliForwardUtil::ParseCollisionSystem(beamT);
  UShort_t sNN = AliForwardUtil::ParseCenterOfMassEnergy(sys, 2 * beamE);
  Short_t  fld = 0; // AliForwardUtil::ParseMagneticField(l3fld);
  Printf("=== Run=%d, year=%d, sys=%d, sNN=%d, fld=%d", 
	 runNo, year, sys, sNN, fld);

  // --- Output object -----------------------------------------------
  // Make our correction object 
  AliFMDCorrAcceptance* corr = new AliFMDCorrAcceptance();
  corr->SetVertexAxis(nVtxBins, vtxLow, vtxHigh);

  // --- Output script -----------------------------------------------
  std::ofstream deadScript("deadstrips.C");
  deadScript << "// Automatically generaeted by ExtractAcceptance.C\n"
	     << "// Add additional dead strips to sharing filter\n"
	     << "// Information taken from OCDB entry for\n" 
	     << "//\n"
	     << "//    run       = " << runNo << "\n"
	     << "//    year      = " << year << "\n"
	     << "//    system    = " << sys << "\n"
	     << "//    sqrt{sNN} = " << sNN << "GeV\n"
	     << "//    L3 field  = " << fld << "kG\n"
	     << "void deadstrips(AliFMDESDFixer* filter)\n"
	     << "{" << std::endl;
  
  // --- Loop over verticies and rings -------------------------------
  Int_t  nDead   = 0;
  Bool_t gotDead = false;
  Float_t dV     = (vtxHigh - vtxLow) / nVtxBins;
  Printf("=== Looping over vertices: %d bins from %+6.2f to %+6.2f "
	 "in steps of %+6.2f", 
	 nVtxBins, vtxLow, vtxHigh, dV);
  for (Double_t v = vtxLow+dV/2; v < vtxHigh; v += dV) { 
    for(UShort_t d = 1; d <= 3;d++) { 
      UShort_t nR = (d == 1 ? 1 : 2);
      for (UShort_t q = 0; q < nR; q++) { 
	Char_t   r  = (q == 0 ? 'I' : 'O');

	// Delegate to other function 
	Int_t nLocal = 0;
	TH2D* ratio = MakeOneRing(d, r, v, nLocal, 
				  !gotDead ? &deadScript : 0);
	nDead += nLocal;
	if (!ratio) {
	  Warning("ExtractAcceptance", "Didn't get correction from "
		  "FMD%d%c @ vz=%+6.2fcm", d, r, v);
	  continue;
	}
	// Printf("v=%+6.2f FMD%d%c, got %d dead strips", v, d, r, nLocal);

	// Set the correction 
	corr->SetCorrection(d, r, v, ratio);
      }
    }
    gotDead = true;
  }
  corr->SetHasOverflow();
  corr->Print();
  // corr->ls();

  deadScript << "}\n"
	     << "// EOF" << std::endl;
  deadScript.close();

  // Write to a file 
  Printf("=== Writing to disk");
  AliForwardCorrectionManager& cm = AliForwardCorrectionManager::Instance();
  if (!cm.Store(corr, runNo, sys, sNN, fld, false, false, 
		"fmd_corrections.root")) { 
    Error("", "Failed to store acceptance correction in local file");
    return;
  }

  std::ofstream f("Upload.C");
  f << "// Generated by ExtractELoss.C\n"
    << "TString MakeDest(const TString& dest, const TString& fname)\n"
    << "{\n"
    << "  TString tmp(dest);\n"
    << "  if (!tmp.IsNull()) {\n"
    << "    if (!tmp.EndsWith(\"/\")) tmp.Append(\"/\");\n"
    << "    tmp.Append(fname);\n"
    << "  }\n"
    << "  return tmp;\n"
    << "}\n\n"
    << "void Upload(const TString& dest=\"\")\n"
    << "{\n"
    << "  gROOT->Macro(\"" << fwd << "/scripts/LoadLibs.C\");\n"
    << "  \n"
    << "  const char* fmdFile = \"fmd_corrections.root\";\n"
    << "  TString fdest = MakeDest(dest, fmdFile);\n"
    << "  \n"
    << "  AliForwardCorrectionManager::Instance().Append(fmdFile, fdest);\n"
    << "}\n"
    << "// EOF\n"
    << std::endl;
  f.close();
  
  Info("ExtracAcceptance", 
       "Run generated Upload.C(DEST) script to copy files in place");
}

//
// EOF
//

  
