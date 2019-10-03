AliJetFinder*  ConfigJetAnalysis()
{
    //
    // Configuration goes here
    // 
    printf("ConfigJetAnalysis() \n");
    AliJetAODReaderHeader *jrh = new AliJetAODReaderHeader();
    jrh->SetComment("AOD Reader");
    jrh->SetPtCut(0.);
    jrh->SetTestFilterMask(1<<0);
    // Define reader and set its header
    AliJetAODReader *er = new AliJetAODReader();
    er->SetReaderHeader(jrh);
   
 
    // Define jet header
    AliDAJetHeader *jh=new AliDAJetHeader();
    jh->SetComment("DA jet code with default parameters");
    jh->SelectJets(kTRUE);
    jh->SetNclust(10);

    // Define jet finder. Set its header and reader
    jetFinder = new AliDAJetFinder();
    jetFinder->SetJetHeader(jh);
    jetFinder->SetJetReader(er);
    //
    return jetFinder;
}
