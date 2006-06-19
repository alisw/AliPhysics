//_________________________________________________________________________
// Macro to create an esd chain and process it with a selector.
// The macro takes as an argument the top directory where the 
// ESDs are stored. Then does the following:
// o) Setup of the par file and load the library libESD.so.
// o) Checks the directories that are one level down
// o) If an AliESDs.root is found it adds it to the esd chain.
// o) Then th chain is processed by the selector esdV0.C
//_________________________________________________________________________
void CreateChain(const char* fDataDir = "/home/pchrist/ALICE/PDC06/pp/data") {
  TStopwatch timer;
  timer.Start();
  
  Char_t fWorkDir[256] = gSystem->pwd();
  const char* pararchivename = "ESD";
  //////////////////////////////////////////
  // Libraries required to load
  //////////////////////////////////////////

  //////////////////////////////////////////////////////////////////
  // Setup PAR File
  if (pararchivename) {
    char processline[1024];
    sprintf(processline,".! tar xvzf %s.par",pararchivename);
    gROOT->ProcessLine(processline);
    gSystem->ChangeDirectory(pararchivename);

    // check for BUILD.sh and execute
    if (!gSystem->AccessPathName("PROOF-INF/BUILD.sh")) {
      printf("*******************************\n");
      printf("*** Building PAR archive    ***\n");
      printf("*******************************\n");

      if (gSystem->Exec("PROOF-INF/BUILD.sh")) {
        Error("batchSelector","Cannot Build the PAR Archive! - Abort!");
        return -1;
      }
    }
    // check for SETUP.C and execute
    if (!gSystem->AccessPathName("PROOF-INF/SETUP.C")) {
      printf("*******************************\n");
      printf("*** Setup PAR archive       ***\n");
      printf("*******************************\n");
      gROOT->Macro("PROOF-INF/SETUP.C");
    }
    
    gSystem->ChangeDirectory("../");
  }

  // Create a Chain
  TChain* fESDChain = new TChain("esdTree");
  void *dirp =  gSystem->OpenDirectory(fDataDir);

  Char_t fPath[256];
  const char * name = 0x0;
  const char * dirname = 0x0;
  const char * pattern = "AliESDs"; 

  while(dirname = gSystem->GetDirEntry(dirp)) {
    sprintf(fPath,"%s/%s",fDataDir,dirname);
    TSystemDirectory* fSystemDir = new TSystemDirectory(".", fPath);
    TList* dirList = fSystemDir->GetListOfFiles();
    //loop over the dir entries - files
    for(Int_t i = 0; i < dirList->GetEntries(); i++) {
      TSystemFile* fFile = (TSystemFile*) dirList->At(i);
      if(strstr(fFile->GetName(),pattern)) {
	Char_t fESDFileName[256];
	sprintf(fESDFileName,"%s/%s.root",fPath,pattern);
	fESDChain->Add(fESDFileName);
      }
    }//loop over files
    delete dirList;
    delete fSystemDir;
  }//loop over dirs
  
  gSystem->ChangeDirectory(fWorkDir);
  fESDChain->Process("esdV0.C"); 

  timer.Stop();
  timer.Print();
}
