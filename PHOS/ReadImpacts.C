ReadImpacts(Int_t nEvents=1,char* file="galice.root")
{

  // Script reads PHOS impacts and prints them
  // Impacts are exacts values of the track coming to EMC, CPV or PPSD
  // and was stored to separate branches of TreeH by AliPHOSvImpacts
  //
  // Yuri Kharlov 4 June 2001

  f = new TFile(file,"readonly");
  (AliRun*)gAlice=(AliRun*)f->Get("gAlice");
  AliPHOSvImpacts *  fPHOS  = (AliPHOSvImpacts *)gAlice->GetDetector("PHOS") ;
  AliPHOSGeometry *  fGeom  = AliPHOSGeometry::GetInstance(fPHOS->GetGeometry()->
                                                           GetName(),
                                                           fPHOS->GetGeometry()->GetTitle() ) ;
  Int_t nPHOSModules = fGeom->GetNModules();
  Int_t nCPVModules  = fGeom->GetNCPVModules();
  Int_t nPPSDModules = fGeom->GetNPPSDModules();

  TBranch * branchEMCimpacts;
  TBranch * branchCPVimpacts;
  TBranch * branchPPSDimpacts;
  TList * fEmcImpacts ;
  TList * fCpvImpacts ;
  TList * fPpsdImpacts ;

  // Loop over events
  for (Int_t iEvent=0; iEvent<nEvents; iEvent++) {
    printf("===========> Event %5d <====================\n",iEvent);
    gAlice->GetEvent(iEvent) ;
    
    // Get branches EMC, CPV and PPSD impacts
    if (! (branchEMCimpacts =gAlice->TreeH()->GetBranch("PHOSEmcImpacts")) )  return 1;
    if (! (branchCPVimpacts =gAlice->TreeH()->GetBranch("PHOSCpvImpacts")) )  return 1;
    if (! (branchPPSDimpacts=gAlice->TreeH()->GetBranch("PHOSPpsdImpacts")) ) return 1;
    
    // Loop over primary tracks
    for (itrack=0; itrack < gAlice->GetNtrack(); itrack++){
      // Set addresses of impacts
      branchEMCimpacts ->SetAddress(&fEmcImpacts) ;
      branchCPVimpacts ->SetAddress(&fCpvImpacts) ;
      branchPPSDimpacts->SetAddress(&fPpsdImpacts) ;
      branchEMCimpacts ->GetEntry(itrack,0);
      branchCPVimpacts ->GetEntry(itrack,0);
      branchPPSDimpacts->GetEntry(itrack,0);
      
      TClonesArray  *impacts;
      AliPHOSImpact *impact;
      Int_t iModule;

      // Do loop over EMC modules
      for (iModule=0; iModule<nPHOSModules; iModule++) {
	impacts = (TClonesArray *)fEmcImpacts->At(iModule);
	// Do loop over impacts in the module
	for (Int_t iImpact=0; iImpact<impacts->GetEntries(); iImpact++) {
	  impact=(AliPHOSImpact*)impacts->At(iImpact);
	  printf("EMC  module %d: ",iModule);
	  impact->Print();
	}
      }

      // Do loop over CPV modules
      for (iModule=0; iModule<nCPVModules; iModule++) {
	impacts = (TClonesArray *)fCpvImpacts->At(iModule);
	// Do loop over impacts in the module
	for (Int_t iImpact=0; iImpact<impacts->GetEntries(); iImpact++) {
	  impact=(AliPHOSImpact*)impacts->At(iImpact);
	  printf("CPV  module %d: ",iModule);
	  impact->Print();
	}
      }

      // Do loop over PPSD modules
      for (iModule=0; iModule<nPPSDModules; iModule++) {
	impacts = (TClonesArray *)fPpsdImpacts->At(iModule);
	// Do loop over impacts in the module
	for (Int_t iImpact=0; iImpact<impacts->GetEntries(); iImpact++) {
	  impact=(AliPHOSImpact*)impacts->At(iImpact);
	  printf("PPSD Module %d: ",iModule);
	  impact->Print();
	}
      }

    }
  }
}
