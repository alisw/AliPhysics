//=============================================================================
// Compilation of all PAR archives for resonance analysis.
// Returns an integer error message if some compilation fails.
//
// Allowed options:
// - "UNTAR"   --> explodes the .par files in the working directory
// - "CLEAR"   --> removes the PAR archive directory and automatically
//                 recreates it by exploding again the .par file
//=============================================================================

Int_t SetupPar(const char* parName, Bool_t untar = kTRUE)
{
//
// Operations to set up the PAR archive
//

    if (!parName) {
        Error("SetupPar", "NULL argument passed - Abort!");
        return -1;
    }

    // if the directory does not exist, the package is
    // un-tarred anyway
    if (!gSystem->OpenDirectory(parName)) untar = kTRUE;

    // unzip + untar (optional)
    if (untar) {
        Char_t processLine[1024];
        sprintf(processLine, ".! tar xvzf ./%s.par", parName);
        gROOT->ProcessLine(processLine);
    }

    // change from working directory to sub-dir containing PAR archive
    const char* ocwd = gSystem->WorkingDirectory();
    gSystem->ChangeDirectory(parName);

    // check for existence of appropriate 'BUILD.sh' script (in PROOF-INF.<parName>)
    // and execute it (non-zero return values mean that errors occurred)
    if (!gSystem->AccessPathName("PROOF-INF/BUILD.sh")) {
        cout << "*** Building PAR archive \"" << parName << "\"" << endl;
        if (gSystem->Exec("PROOF-INF/BUILD.sh")) {
            Error("SetupPar", "Cannot Build the PAR Archive - Abort!");
            return -1;
        }
    }

    // check for existence of appropriate 'SETUP.C' macro (in PROOF-INF.<parName>)
    // and execute it
    if (!gSystem->AccessPathName("PROOF-INF/SETUP.C")) {
        cout << "*** Setting up PAR archive \"" << parName << "\"" << endl;
        gROOT->Macro("PROOF-INF/SETUP.C");
    }

    // return to working directory
    gSystem->ChangeDirectory("../");

    // return 1 if everything is OK
    return 1;
}

Int_t AliRsnLoad(Option_t *option = "BUILD+UNTAR")
{
//
// Main function of this macro
// Sets up all the packages useful for resonance analysis.
// Depending on the support (PROOF / others), the functions are different.
//

    TString opt(option);
    opt.ToUpper();

    // check for the UNTAR option
    Bool_t untar = opt.Contains("UNTAR");

    // check for the CLEAR option
    // if this is found automatically is required the UNTAR
    if (opt.Contains("CLEAR")) {
        gSystem->Exec("rm -rf STEERBase ESD AOD ANALYSIS ANALYSISalice PWG2resonances");
        untar = kTRUE;
    }

    // build all libraries with the defined options
    SetupPar("STEERBase", untar);
    SetupPar("ESD", untar);
    SetupPar("AOD", untar);
    SetupPar("ANALYSIS", untar);
    SetupPar("ANALYSISalice", untar);
    SetupPar("PWG2resonances", untar);
}
