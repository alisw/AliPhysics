void AliITSPrintGeom(TString hfn="galice.root",Int_t mod=-1){
    // Macro to print out the information kept in the AliITSgeom class, for 
    // all or a specific module.
    // Inputs:
    //   TString hfn    input file name
    //   Int_t   mod    The specific module to print out transformations for.
    //                  if mod<0, does all modules.
    // Output:
    //    none.
    // Return:
    //    none.

    // Dynamically link some shared libs
    if (gClassTable->GetID("AliRun") < 0) {
        gROOT->LoadMacro("loadlibs.C");
        loadlibs();
    }else {
	if(gAlice){
	    delete AliRunLoader::GetRunLoader();
	    delete gAlice;
	    gAlice=0;
	} // end if gAlice
    } // end if aliroot or not.

    AliRunLoader* rl = AliRunLoader::Open(hfn.Data());
    if (rl == 0x0){
	cerr<<"AliITSPrintGeom.C : Can not open session RL=NULL"<< endl;
	return;
    } // end if no loader

    Int_t retval = rl->LoadgAlice();
    if (retval){
	cerr<<"AliITSHits2SDigits.C : LoadgAlice returned error"<< endl;
	return 3;
    } // end if loader error
    if(!gAlice){
	cerr<<"AliITSPrintGeom.C. AliRun object not found\n";
	return;
    } // end if no gAlice error

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
    Int_t m,lay,lad,det;
    Double_t xyz[3],rcyl;
    cout<<endl<<endl<<"====================================\n";
    if(mod<-1){ gmm->PrintComment(&cout); cout << endl;}
    cout<<endl<<endl<<"====================================\n";
    for(m=mod1;m<mod2;m++){
	gm->GetModuleId(m,lay,lad,det);
	gmm = gm->GetGeomMatrix(m);
	gmm->GetTranslation(xyz);
	rcyl = TMath::Sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]);
	if(mod<-1){
	    gmm->Print(&cout);
	    cout << endl;
	}else{
	    cout << "module="<<m<<" "<<lay<<" "<<lad<<" "<<det;
	    cout << " Rcyl="<<rcyl<<" cm"<<endl;
	} // end if
    } // end for m
}
