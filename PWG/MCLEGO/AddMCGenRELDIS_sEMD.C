//#include "AliGenerator.h"
//#include "AliGenReadersEMD.h"
//#include "AliGenExtFile.h"

AliGenerator* AddMCGenRELDIS_sEMD(TString reldisFile="epbpb2510nt_na50.root", TString ntupleName="h2034",Int_t firstev=0) {


    AliGenExtFile *gener = new AliGenExtFile(-1);

    AliGenReadersEMD *reader = new AliGenReadersEMD();
    // need to resolve files from grid, make sure we have a valid connection
    if (!gGrid) TGrid::Connect("alien");
    reader->SetFileName(Form("alien:///alice/cern.ch/user/p/pwg_mm/%s", reldisFile.Data())); // put the location of the input file
    reader->SetNtupleName(ntupleName); // name of the tree inside the input file
    reader->SetStartEvent(firstev); // # of event to start with
    reader->TrackOnlyNeutrons(); // include this if you want to track only neutrons

    gener->SetReader(reader); 

    return gener;
}
