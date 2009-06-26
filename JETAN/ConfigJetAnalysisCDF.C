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
    jrh->SetFiducialEta (-0.9,0.9) ;
    jrh->SetPtCut(0.5);
    //jrh->SetReadSignalOnly(kFALSE);
    // Define reader and set its header

    // ESD Reader
    AliJetESDReader *er = new AliJetESDReader();
    er->SetReaderHeader(jrh);

     // Define jet header
    AliCdfJetHeader * jh = new AliCdfJetHeader();
    jh->SetRadius(0.7);
      cout << "Radius = " << jh->GetRadius() << endl;

    jh->SetMinPartJet(1);
      cout << "Min Part Jet = " << jh->GetMinPartJet ();

    // Cuts in Pt of Jets
    jh->SetJetPtCut ( 0. );
//    jh->SetJetPtCut ( 2. );
//    jh->SetJetPtCut ( 5. );
//    jh->SetJetPtCut ( 30. );

    jh->SetDebugCDF( true ) ;



    // Define jet finder
    AliCdfJetFinder *jetFinder = new AliCdfJetFinder();
    jetFinder->SetJetHeader(jh);
    jetFinder->SetJetReader(er);
//    jetFinder->SetOutputFile("jets.root");

    return jetFinder;
}

