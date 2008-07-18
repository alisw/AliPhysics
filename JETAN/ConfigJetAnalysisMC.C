AliJetFinder*  ConfigJetAnalysis()
{
    //
    // Configuration goes here
    // 
    printf("ConfigJetAnalysis() \n");
    Printf("For MC events");
    AliJetKineReaderHeader *jrh = new AliJetKineReaderHeader();
    jrh->SetComment("MC full Kinematics");
    jrh->SetFastSimTPC(kFALSE);
    jrh->SetFastSimEMCAL(kFALSE);
    jrh->SetPtCut(0.);

    // Define reader and set its header
    AliJetKineReader *er = new AliJetKineReader();
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
    jetFinder->SetOutputFile("jetsMC.root");
    //
    return jetFinder;
}
