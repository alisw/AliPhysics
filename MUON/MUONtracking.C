void MUONtracking (Int_t evNumber1=0, Int_t evNumber2=99, Int_t idres=116, Int_t ireadgeant=1, Int_t ibgr=1, Int_t nev_bgd=10) 
{
  //////////////////////////////////////
  //                                  //
  // ROOT macro for ALICE Dimuon Arm: //
  // Track reconstruction             //
  //                                  //
  //////////////////////////////////////
  //
  // Reconstructs tracks from events in the ROOT file "galice.root".
  // Track reconstruction is performed (argument "ireadgeant")
  // either directly from GEANT hits (tree TH),
  // or from raw clusters (tree TR) constructed from digits (tree TD).
  // Eventually (argument "ibgr"), background GEANT hits
  // are also taken into account from the ROOT file "galice_bgr.root".
  // 
  // Arguments:
  //   evNumber1  = first event number to act on in file "galice.root"
  //   evNumber2  = last event number to act on in file "galice.root"
  //   idres      : used for statistics
  //              = 116 for Upsilon
  //              = 114 for J/Psi
  //   ireadgeant = 1 to reconstruct tracks directly from GEANT hits
  //              = 0 to reconstruct tracks from raw clusters
  //   ibg        : used only if "ireadgeant" = 1
  //              = 0 if no background GEANT hits to be taken into account;
  //              = 1 if  background GEANT hits
  //                  to be taken into account in file "galice_bgr.root"
  //              used only if "ireadgeant" = 1
  //   nev_bgd    : used only if "ireadgeant" = 1 and "ibg" = 1
  //              = number of events in the background file "galice_bgr.root";
  //                successive signal events are mixed
  //                with different background events in file "galice_bgr.root",
  //                  starting from event number 0,
  //                  incrementing it by 1 till event number ("nev_bgd" - 1),
  //                  continuing with event number 0 and so on.
  //                Strictly speaking, "nev_bgd" can be smaller than
  //                the number of events in the background file,
  //                in which case one will only use
  //                the first "nev_bgd" events of the background file.
  //                But it SHOULD NOT BE LARGER THAN
  //                THE ACTUAL NUMBER OF EVENTS IN THE BACKGROUND FILE.
  //
  // Input file(s):
  //   "galice.root" for GEANT hits or raw clusters
  //   "galice_bgr.root" for background GEANT hits
  //                     (used only if "ireadgeant" = 1 and "ibg" = 1)
  //
  // Output file:
  //   "reconst.root" for ROOT ntuples
  //
  //__________________________________________________________________________

// Dynamically link some shared libs                    

   if (gClassTable->GetID("AliRun") < 0) {
      gROOT->LoadMacro("loadlibs.C");
      loadlibs();
   }
// Connect the Root Galice file containing Geometry, Kine and Hits

    TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject("galice.root");
    
    if (!file) {
	printf("\n Creating galice.root \n");
	file = new TFile("galice.root");
    } else {	printf("\n galice.root found in file list");
    }
    
// Get AliRun object from file or create it if not on file
    if (!gAlice) {
	gAlice = (AliRun*)(file->Get("gAlice"));
	if (gAlice) printf("AliRun object found on file\n");
	if (!gAlice) {
	    printf("\n create new gAlice object");
	    gAlice = new AliRun("gAlice","Alice test program");
	}
    }

// seff = efficiency per chamber (ireadgeant=1) 
    Double_t seff  = 1;
    //    Double_t seff  = 1.;
// sb0 = magn. field in dipole, sbl3 = magn. field in L3 
// necessary for trackfinding only.
    Double_t sb0   = 0.7;
    Double_t sbl3  = 0.2;

//  ifit = 0 trackfinding only
//  ifit = 1 trackfinding + fit 
    Int_t ifit   = 1;
// idebug = 0,1,2 print level for reco_muon.F 
    Int_t idebug = 1;

    AliMUON *MUON  = (AliMUON*) gAlice->GetModule("MUON");
    AliMUONTrackReconstructor *Reconstruction = new AliMUONTrackReconstructor();    
    Int_t nparticles = gAlice->GetEvent(evNumber1);
    if (nparticles <= 0) return;

    Reconstruction->Init(seff,sb0,sbl3);

//   Loop over events 
//
    Int_t inev_bgd=0;
    
    for (Int_t nev= evNumber1; nev<= evNumber2; nev++) 
    {
	printf("nev=%d\n",nev);
	if (nev != evNumber1) Int_t nparticles = gAlice->GetEvent(nev);
	if (nev < evNumber1) continue;
	if (nparticles <= 0) return;
	Reconstruction->FinishEvent();
	if (ireadgeant==1 && ibgr==1) {
	  if (inev_bgd==nev_bgd) inev_bgd=0;
	  Reconstruction->Reconst(ifit,idebug,inev_bgd,nev,idres,ireadgeant,"Add","galice_bgr.root");
	  Reconstruction->FinishEvent();
	  inev_bgd++;
	}
	else {
	  //	printf("ireadgeant=%d\n",ireadgeant);
        Reconstruction->Reconst(ifit,idebug,inev_bgd,nev,idres,ireadgeant,"rien1","rien2");
	Reconstruction->FinishEvent();

	}

    } // event loop 
    
    Reconstruction->Close();
}



