Int_t setupPar(const char* pararchivename);

void LoadLibraries(Bool_t useParFiles=kFALSE) {
    gSystem->SetIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_ROOT/ANALYSIS -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TPC -I$ALICE_ROOT/CONTAINERS -I$ALICE_ROOT/STEER/STEER -I$ALICE_ROOT/STEER/STEERBase -I$ALICE_ROOT/STEER/CDB -I$ALICE_ROOT/STEER/ESD -I$ALICE_ROOT/STEER/AOD -I$ALICE_ROOT/PWG/Tools -I$ALICE_ROOT/PWGLF/STRANGENESS -I$ALICE_ROOT/TRD -I$ALICE_ROOT/macros -I$ALICE_ROOT/OADB  -g");
	
    gSystem->AddIncludePath("-I$ALICE_ROOT/include");
    gSystem->AddIncludePath("-I$ALICE_ROOT/TOF");
    gSystem->AddIncludePath("-I$ALICE_ROOT/PWG/Tools");
    gSystem->AddIncludePath("-I$ALICE_ROOT/PWGLF/STRANGENESS");
    /* load libraries */
    
    
        gSystem->Load("libTree.so");
	gSystem->Load("libCore.so");
	gSystem->Load("libGeom.so");
	gSystem->Load("libPhysics.so");
	gSystem->Load("libVMC.so");
	gSystem->Load("libMinuit.so");
	gSystem->Load("libSTEERBase.so");
	gSystem->Load("libESD.so");
	gSystem->Load("libAOD.so");
	gSystem->Load("libANALYSIS.so");
	gSystem->Load("libOADB.so");
	gSystem->Load("libANALYSISalice.so");
	gSystem->Load("libCORRFW.so");
	gSystem->Load("libGui.so");
	gSystem->Load("libProof.so");
	gSystem->Load("libXMLParser.so");
	gSystem->Load("libRAWDatabase.so");
	gSystem->Load("libRAWDatarec.so");
	gSystem->Load("libCDB.so");
	gSystem->Load("libSTEER.so");
	gSystem->Load("libTOFbase.so");
	gSystem->Load("libPWGTools.so");
	gSystem->Load("libPWGLFSTRANGENESS.so");

	//gSystem->Load("libTOFrec.so");
    return;
}

//------------------------------------------------------------------------
Int_t setupPar(const char* pararchivename) {
    ///////////////////
    // Setup PAR File//
    ///////////////////
    
    if (pararchivename) {
        char processline[1024];
        TString base = gSystem->BaseName(pararchivename);
        TString dir  = gSystem->DirName(pararchivename);
        TString ocwd = gSystem->WorkingDirectory();
        // Move to dir where the par files are and unpack
        gSystem->ChangeDirectory(dir.Data());
        sprintf(processline,".! tar xvzf %s.par",base.Data());
        gROOT->ProcessLine(processline);
        // Move to par folder
        gSystem->ChangeDirectory(base.Data());
        
        // check for BUILD.sh and execute
        if (!gSystem->AccessPathName("PROOF-INF/BUILD.sh")) {
            printf("*******************************\n");
            printf("*** Building PAR archive    ***\n");
            printf("*******************************\n");
            
            if (gSystem->Exec("PROOF-INF/BUILD.sh")) {
                Error("runAnalysis","Cannot Build the PAR Archive! - Abort!");
                gSystem->ChangeDirectory(ocwd.Data());
                return -1;
            }
        }
        // check for SETUP.C and execute
        if (!gSystem->AccessPathName("PROOF-INF/SETUP.C")) {
            printf("*******************************\n");
            printf("*** Setup PAR archive       ***\n");
            printf("*******************************\n");
            // If dir not empty, set the full include path
            if (dir.Length()) {
                sprintf(processline, ".include %s", pararchivename);
                gROOT->ProcessLine(processline);
            }
            gROOT->Macro("PROOF-INF/SETUP.C");
        }
        gSystem->ChangeDirectory(ocwd.Data());
    }
    return 1;
}
