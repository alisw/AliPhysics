void AliITSPrintGeom(TString hfn="galice.root",Int_t mod=-1){
    // Macro to print out the information kept in the AliITSgeom class, for 
    //all or a specific module

    // Dynamically link some shared libs
    if (gClassTable->GetID("AliRun") < 0) {
        gROOT->LoadMacro("loadlibs.C");
        loadlibs();
    } 
    else {
      if(gAlice){
	delete gAlice->GetRunLoader();
	delete gAlice;
	gAlice=0;
      }
    }

    AliRunLoader* rl = AliRunLoader::Open(hfn.Data());
    if (rl == 0x0){
      cerr<<"AliITSPrintGeom.C : Can not open session RL=NULL"<< endl;
      return;
    }

    Int_t retval = rl->LoadgAlice();
    if (retval)
     {
      cerr<<"AliITSHits2SDigits.C : LoadgAlice returned error"
           << endl;
       return 3;
     }
    //    gAlice=rl->GetAliRun();
    if(!gAlice){
      cerr<<"AliITSPrintGeom.C. AliRun object not found\n";
      return;
    }

    AliITS *ITS = (AliITS*)gAlice->GetDetector("ITS");
    if(!ITS){
        cout << "Error: no ITS found. Aborting"<<endl;
        return;
    } // end if !ITS

    AliITSgeom *gm = ITS->GetITSgeom();
    Int_t mod1 = 0;
    Int_t mod2 = gm->GetIndexMax();
    if(mod>=0){
        mod1 = mod;
        mod2 = mod+1;
    } // end if mod>=0
    AliITSgeomMatrix *gmm = gm->GetGeomMatrix(0);
    Int_t m;
    cout<<endl<<endl<<"====================================\n";
    gmm->PrintComment(&cout); cout << endl;
    cout<<endl<<endl<<"====================================\n";
    for(m=mod1;m<mod2;m++){
	gmm = gm->GetGeomMatrix(m);
	gmm->Print(&cout); cout << endl;
    } // end for m
}
