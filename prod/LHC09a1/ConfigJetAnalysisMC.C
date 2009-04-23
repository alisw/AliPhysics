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
    jrh->SetFiducialEta(-2.1,2.1);
    jrh->SetPtCut(0.);

    // Define reader and set its header
    AliJetKineReader *er = new AliJetKineReader();
    er->SetReaderHeader(jrh);
   
 
    // Define jet header
    AliUA1JetHeaderV1 *jh=new AliUA1JetHeaderV1();
    jh->SetComment("UA1 jet code with MC parameters");
    jh->BackgMode(0);
    jh->SetRadius(1.0);
    jh->SetEtSeed(4.);
    jh->SetLegoNbinPhi(432.);
    jh->SetLegoNbinEta(274.);
    jh->SetLegoEtaMin(-2);
    jh->SetLegoEtaMax(+2);
    jh->SetJetEtaMin(-1.5);  
    jh->SetJetEtaMax(1.5);  
    jh->SetMinJetEt(5.);
    
    // Define jet finder. Set its header and reader
    jetFinder = new AliUA1JetFinderV1();
    jetFinder->SetJetHeader(jh);
    jetFinder->SetJetReader(er);
    return jetFinder;
}
