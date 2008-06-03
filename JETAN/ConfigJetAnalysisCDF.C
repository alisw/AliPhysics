AliJetFinder*  ConfigJetAnalysis()
{
    //
    // Configuration goes here
    // 
    printf("ConfigJetAnalysis() \n");
    
    // Reader Header
    AliJetESDReaderHeader *jrh = new AliJetESDReaderHeader();
    jrh->SetComment("Testing");
    jrh->SetFirstEvent(0);
    jrh->SetLastEvent(40000);
    //jrh->SetPtCut(0.);
    //jrh->SetReadSignalOnly(kFALSE);
    // Define reader and set its header
    
    // ESD Reader
    AliJetESDReader *er = new AliJetESDReader();
    er->SetReaderHeader(jrh);
   
     // Define jet header
    AliCdfJetHeader * jh = new AliCdfJetHeader();
    jh->SetRadius(0.7);
 
    // Define jet finder
    AliCdfJetFinder *jetFinder = new AliCdfJetFinder();
    jetFinder->SetJetHeader(jh);
    jetFinder->SetJetReader(er);
    jetFinder->SetOutputFile("jets.root");
    
    return jetFinder;
}
