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
    jrh->SetFiducialEta(-2.1,2.1); // to take all MC particles default is 0.9

    // Define reader and set its header
    AliJetKineReader *er = new AliJetKineReader();
    er->SetReaderHeader(jrh);
   
 
    // Define jet header
    AliUA1JetHeaderV1 *jh=new AliUA1JetHeaderV1();
    jh->SetComment("UA1 jet code with  parameters similar to PYCELL during sim.");
    jh->BackgMode(0);
    jh->SetRadius(1.0);
    jh->SetEtSeed(4.);
    jh->SetLegoNbinPhi(432);
    jh->SetLegoNbinEta(274);
    jh->SetLegoEtaMin(-2);
    jh->SetLegoEtaMax(+2);  
    jh->SetJetEtaMax(0.5);  
    jh->SetJetEtaMin(-0.5);  
    jh->SetMinJetEt(10.);
    

    // Define jet finder. Set its header and reader
    jetFinder = new AliUA1JetFinderV1();
    jetFinder->SetJetHeader(jh);
    jetFinder->SetJetReader(er);
    return jetFinder;
}
