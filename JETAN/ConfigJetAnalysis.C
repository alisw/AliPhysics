AliJetFinder*  ConfigJetAnalysis()
{
    //
    // Configuration goes here
    // 
    printf("ConfigJetAnalysis() \n");
    AliJetESDReaderHeader *jrh = new AliJetESDReaderHeader();
    jrh->SetComment("Testing");
    jrh->SetFirstEvent(0);
    jrh->SetLastEvent(1000);
    jrh->SetPtCut(0.);
    jrh->SetReadSignalOnly(kFALSE);
    // Define reader and set its header
    AliJetESDReader *er = new AliJetESDReader();
    er->SetReaderHeader(jrh);
   
 
    // Define jet header
    AliUA1JetHeaderV1 *jh=new AliUA1JetHeaderV1();
    jh->SetComment("UA1 jet code with default parameters");
    jh->BackgMode(0);
    jh->SetRadius(1.0);
    jh->SetEtSeed(2.);
    jh->SetLegoNbinPhi(420.);
    jh->SetLegoNbinEta(120.);
    jh->SetLegoEtaMin(-1.9);
    jh->SetLegoEtaMax(+1.9);  
    jh->SetMinJetEt(5.);
    
    // Define jet finder. Set its header and reader
    jetFinder = new AliUA1JetFinderV1();
    jetFinder->SetJetHeader(jh);
    jetFinder->SetJetReader(er);
    jetFinder->SetPlotMode(kTRUE);
    jetFinder->SetOutputFile("jets.root");
    //
    return jetFinder;
}
