// Functions used to access input files - shared by several macros
 
//______________________________________________________________________
Bool_t GaliceITSok(){
    // Checks gAlice to see that ITS and the ITS geometry are properly
    // defined. If not return kFALSE and deletes gAlice and set it to zero.

    if(!gAlice){
	return kFALSE;
    } // end if !gAlice
    // gAlice defined check to see if ITS is properly defined.
    AliITS *ITS = (AliITS*)gAlice->GetDetector("ITS"); 
    if(!ITS){ // ITS not defined, delete and reload gAlice
	delete gAlice;
	gAlice = 0;
	return kFALSE;
    } // end if !ITS
    // ITS defined
    if(!(ITS->GetITSgeom())){
	delete gAlice;
	gAlice = 0;
	return kFALSE;
    } // end if !(ITS->GetITSgeom())
    // ITS and ITS geometry properly defined defined.
    return kTRUE;
}
//______________________________________________________________________
AliRunLoader* AccessFile(TString FileName){
  // Function used to open the input file and fetch the AliRun object

  AliRunLoader* rl = AliRunLoader::Open(FileName.Data());
  if (rl == 0x0){
    cerr<<"AccessFile : Can not open session RL=NULL"<< endl;
    return rl;
  }
 
  Int_t retval = rl->LoadgAlice();
  if (retval){
    cerr<<"AccessFile : LoadgAlice returned error"<<endl;
    return rl;
  }
  return rl;
}

