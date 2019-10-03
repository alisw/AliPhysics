//#include "AliGenerator.h"
//#include "AliGenReadersEMD.h"
//#include "AliGenExtFile.h"

AliGenerator* AddMCGenRELDIS_EMD(TString reldisFile="mpbpb2510nt_cwn.root", TString ntupleName="h2034",Int_t firstev=0) {


    AliGenExtFile *gener = new AliGenExtFile(-1);

    AliGenReaderEMD *reader = new AliGenReaderEMD();
    // need to resolve files from grid, make sure we have a valid connection
    if (!gGrid) TGrid::Connect("alien");
    reader->SetFileName(Form("alien:///alice/cern.ch/user/p/pwg_mm/%s", reldisFile.Data())); // put the location of the input file
    reader->SetNtupleName(ntupleName); // name of the tree inside the input file
    reader->SetStartEvent(firstev); // # of event to start with

    gener->SetReader(reader); 

    return gener;
}
