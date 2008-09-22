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
    // Define reader and set its header
    // ESD Reader
    AliJetESDReader *er = new AliJetESDReader();
    er->SetReaderHeader(jrh);
   
     // Define jet header
    AliFastJetHeader * jh = new AliFastJetHeader();
    jh->SetRparam(0.7); // setup parameters
 
    // Define jet finder
    AliFastJetFinder *jetFinder = new AliFastJetFinder();
    jetFinder->SetJetHeader(jh);
    jetFinder->SetJetReader(er);
    jetFinder->SetOutputFile("jets.root");
    
    return jetFinder;
}
