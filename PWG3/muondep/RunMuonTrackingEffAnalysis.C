// 2008
// Macro for the running of the AliAnalysisTaskMuonTrackingEff
//

void RunMuonTrackingEffAnalysis (Bool_t alien = false,
				 const char * macroFileName = "$ALICE_ROOT/PWG3/muon/MuonTrackingEffAnalysis.C",
				 const char * esdfileName = "AliESDs.root",
				 const char * geometryFileName = "geometry.root",
				 const char * analysisParFile = "ANALYSIS",
				 const char * pwg3ParFile = "PWG3",
				 const Int_t run = 100)
{
    if(alien)
    {
    //Grid connection    
      printf("*** Connect to AliEn ***\n");
      TGrid::Connect("alien://");
      gSystem->Load("libProofPlayer.so");
    }

// //Load relevant libraries:
//     gSystem->Load("libTree.so");
//     gSystem->Load("libGeom.so");
//     gSystem->Load("libSTEERBase.so");

// //  setupPar("MUON");
//     gSystem->Load("libMUONbase.so");
//     gSystem->Load("libMUONgeometry.so");
//     gSystem->Load("libMUONmapping.so");

// //  setupPar("ESD");
//     gSystem->Load("libESD.so");

//  setupPar(analysisParFile);
    gSystem->Load("libANALYSIS.so");
    gSystem->Load("libANALYSISalice.so");
 
//  setupPar(pwg3ParFile);
    gSystem->Load("libPWG3muon.so");

    char macro[1024];
    sprintf(macro,"%s++",macroFileName);
    gROOT->LoadMacro(macro);
//     gROOT->LoadMacro("./BatchMuonTrackingEffAnalysis.C++");
 
    MuonTrackingEffAnalysis(alien,run,esdfileName,geometryFileName);
}



Int_t setupPar(const char* pararchivename)
{
    //Setup PAR File
    if (pararchivename) {
      char processline[1024];
      sprintf(processline,".! tar xvzf %s.par",pararchivename);
      gROOT->ProcessLine(processline);
      const char* ocwd = gSystem->WorkingDirectory();
      gSystem->ChangeDirectory(pararchivename);
      
      // check for BUILD.sh and execute
      if (!gSystem->AccessPathName("PROOF-INF/BUILD.sh")) {
	printf("*******************************\n");
	printf("*** Building PAR archive    ***\n");
	printf("*******************************\n");
      
	if (gSystem->Exec("PROOF-INF/BUILD.sh")) {
	  Error("runAnalysis","Cannot Build the PAR Archive! - Abort!");
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
    return 1;
}
