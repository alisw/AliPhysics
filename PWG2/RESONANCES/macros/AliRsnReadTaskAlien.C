void AliRsnReadTaskAlien
(
    const char *kCollectionFile="wn.xml",            // XML file containing tags
    Long64_t    nentries=TChain::kBigNumber
)
{
    // connect to grid
    TGrid::Connect("alien://");
    
    // create chain of files to read
    TAlienCollection *myCollection = TAlienCollection::Open(kCollectionFile);
    if (!myCollection) {
        Error("AliRsnReadTaskRL", Form("Cannot create an AliEn collection from %s", kCollectionFile));
        return;
    }
    TChain* analysisChain = new TChain("esdTree");
    myCollection->Reset();
    
    // loop on the entries of the XML input file
    Int_t i = 0;
    while ( myCollection->Next() ) {
        char esdFile[255];
        sprintf(esdFile, "%s", myCollection->GetTURL(""));
        Info("AliRsnReadTaskRL", Form("Adding %s", esdFile));
        analysisChain->Add(esdFile);
        if (++i >= 1) break;
    }
    Info("AliRsnReadTaskRL", Form("CHAIN HAS %d ENTRIES", (Int_t)analysisChain->GetEntries()));
    
    // load read macro
    gROOT->LoadMacro("$(ALICE_ROOT)/PWG2/RESONANCES/macros/AliRsnReadTask.C");
    AliRsnReadTask(analysisChain);
}