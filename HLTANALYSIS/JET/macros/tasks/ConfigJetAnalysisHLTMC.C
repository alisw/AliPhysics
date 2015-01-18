/** ---------------------------------------------------------------------
 *  @file   ConfigJetAnalysisHLTMC.C
 *  @author Jochen Thaeder <jochen@thaeder.de>
 *  @brief  Run HLT cone finder in analysis framework, on Kinematics
 *
 *  --------------------------------------------------------------------- 
 */

AliJetFinder*  ConfigJetAnalysis() {
 
  printf("ConfigJetAnalysis() -- HLT \n");
 
  // ---------------------------------------------------------------------
  // -- Defaults
  // ---------------------------------------------------------------------

  TString comment       = "HLT Fast Fixed Seeded Cone finder on MC";
  AliHLTJETBase::JetAlgorithmType_t algorithm = AliHLTJETBase::kFFSCSquareCell;
  
  Bool_t  leading       = kFALSE;
  Float_t coneRadius    =  0.4;
  Float_t trackCutMinPt =  1.0;
  Float_t seedCutMinPt  =  5.0;
  Float_t jetCutMinEt   = 15.0;
  Bool_t  useMC         = kTRUE;

  // -- Jet Track Cuts
  // ---------------------------------------------------------------------
  AliHLTJETTrackCuts *trackCuts = new AliHLTJETTrackCuts();
  trackCuts->SetChargedOnly( kTRUE );
  trackCuts->SetMinPt( trackCutMinPt );
  
  // -- Jet Seed Cuts
  // ---------------------------------------------------------------------
  AliHLTJETConeSeedCuts *seedCuts = new AliHLTJETConeSeedCuts();
  seedCuts->SetMinPt( seedCutMinPt );

  // -- Jet Jet Cuts
  // ---------------------------------------------------------------------
  AliHLTJETJetCuts *jetCuts = new AliHLTJETJetCuts();
  jetCuts->SetMinEt( jetCutMinEt );

  // -- Jet Reader Header
  // ---------------------------------------------------------------------
  AliHLTJETReaderHeader *jetReaderHeader = new AliHLTJETReaderHeader();

  // Set Algorithm 
  jetReaderHeader->SetJetAlgorithm(algorithm);

  // Set prt to track cuts
  jetReaderHeader->SetTrackCuts(trackCuts);
  jetReaderHeader->SetSeedCuts(seedCuts);

  // Set Eta min/max and Phi min/max
  jetReaderHeader->SetFiducialEta( -0.9, 0.9) ;
  jetReaderHeader->SetFiducialPhi(  0.0, TMath::TwoPi() ) ;

  // Set grid binning
  jetReaderHeader->SetGridEtaBinning( 0.05 );
  jetReaderHeader->SetGridPhiBinning( 0.05 );
 
  // Set cone radius
  jetReaderHeader->SetConeRadius(coneRadius);

  // Use Kinematics
  jetReaderHeader->SetUseMC(useMC);

  // -- Jet Reader
  // ---------------------------------------------------------------------
  AliHLTJETReader *jetReader = new AliHLTJETReader();
  jetReader->SetReaderHeader(jetReaderHeader);

  // ---------------------------------------------------------------------
  // -- Jet Container
  // ---------------------------------------------------------------------
  AliHLTJets *jets = new AliHLTJets();
  jets->SetComment(comment);

  // ---------------------------------------------------------------------
  // -- Jet Header
  // ---------------------------------------------------------------------
  AliHLTJETConeHeader *jetHeader = new AliHLTJETConeHeader();
  jetHeader->SetJetCuts(jetCuts);
  jetHeader->SetUseLeading(leading);

  // ---------------------------------------------------------------------
  // -- Jet Finder
  // ---------------------------------------------------------------------
  AliHLTJETConeFinder *jetFinder = new AliHLTJETConeFinder();
  jetFinder->SetJetHeader(jetHeader);
  jetFinder->SetJetReader(jetReader);
  jetFinder->SetOutputJets(jets);

  return jetFinder;
}
