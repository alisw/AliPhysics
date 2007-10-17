
AliAnaGamma*  ConfigGammaAnalysis()
{
    //
    // Configuration goes here
    // 
    printf("ConfigGammaAnalysis() \n");

    //----------------------------------------------------------
    //-----------------------------------------------------------  
    // Define reader , uncomment 1 of the 2 options
    //-----------------------------------------------------------
    //----------------------------------------------------------

    // -----Option 1------ Data, ESDs
    //AliGammaDataReader *reader = new AliGammaDataReader();
    //       AliGammaMCDataReader *reader = new AliGammaMCDataReader(); //Copy of AliGammaDataReader + Montecarlo Information
//     //Set detectors acceptance (in real and montecarlo data)
//     reader->SetCTSEtaCut(1); reader->SetEMCALEtaCut(1); reader->SetPHOSEtaCut(0.2);
//     reader->SetPhiEMCALCut(40*TMath::DegToRad(), 200*TMath::DegToRad());     
//     reader->SetPhiPHOSCut(200*TMath::DegToRad(), 350*TMath::DegToRad()); 
//     //Set minimum pt for particles in analysis  (in real and montecarlo data)
//     reader->SetNeutralPtCut(0.4); reader->SetChargedPtCut(0.4); 
//     //pid of measured particles
    //reader->SetEMCALPIDOn(kTRUE); reader->SetPHOSPIDOn(kTRUE); //No pid, accept all particles 
//     //if previous kTrue
//     reader->SetPHOSPhotonWeight(0.7);    reader->SetPHOSPi0Weight(0.7); 
//     reader->SetEMCALPhotonWeight(0.7);    reader->SetEMCALPi0Weight(0.7);
    //reader->UsePHOSPIDWeightFormula(kTRUE);
    //TFormula * photonF = new TFormula("photonWeight","0.98*(x<40)+ 0.68*(x>=100)+(x>=40 && x<100)*(0.98+x*(6e-3)-x*x*(2e-04)+x*x*x*(1.1e-06))");
    //TFormula * pi0F = new TFormula("pi0Weight","0.98*(x<65)+ 0.915*(x>=100)+(x>=65 && x-x*(1.95e-3)-x*x*(4.31e-05)+x*x*x*(3.61e-07))");
    //reader->SetPHOSPhotonWeightFormula(photonF);
    //reader->SetPHOSPi0WeightFormula(pi0F);

     // -----Option 2------ Kinematics
       AliGammaMCReader *reader = new AliGammaMCReader();
//     //Set detectors acceptance (in real and montecarlo data)
//     reader->SetCTSEtaCut(1); reader->SetEMCALEtaCut(1); reader->SetPHOSEtaCut(0.2);
//     reader->SetPhiEMCALCut(40*TMath::DegToRad(), 200*TMath::DegToRad());     
//     reader->SetPhiPHOSCut(200*TMath::DegToRad(), 350*TMath::DegToRad()); 
//     //Set minimum pt for particles in analysis  (in real and montecarlo data)
//     reader->SetNeutralPtCut(0.4); reader->SetChargedPtCut(0.4); 
     reader->SetDecayPi0Flag(AliGammaMCReader::kGeantDecay) ; //Options
//     //kNoDecay :Do not decay pi0, keep them in the list as they are
//     //kGeantDecay: Look for gamma decayed by GEANT
//     //kDecay: Decay pi0 by hand (geant was not used)
//     //kDecayGamma: Pi0 is decayed by PYTHIA, pi0 is not final check if photons overlapp
//     //parameters to study if decay is overlapped:
//     reader->SetEMCALIPDistance(450.); reader->SetPHOSIPDistance(460.);
//     reader->SetEMCALMinDistance(3.6);    reader->SetPHOSMinDistance(11.); //Miimum overlapp distance
       reader->SetCheckOverlapping(kTRUE);
 
    //----------------------------------------------------------
    //----------------------------------------------------
    //Define analysis algorithms
    //----------------------------------------------------
    //----------------------------------------------------------
    //3 analysis types
    //kPrompt: Find prompt gamma for a fixed cone and pt cut value
    //kIsolationCut: Find prompt gamma for several cones and pt cuts values
    //kCorrelation: Do prompt gamma - something correlation. Something can be
    //      kParton: Gamma-Parton correlation
    //      kHadron: Gamma-Hadron correlation
    //      kJetLeadCone: Gamma-Jet correlation : constructed in cone around leading particle
    //      kJetFinder: Gamma-Jet correlation : Jet reconstructed with standard algorithms. --Still not implemented--
    //One of the first 3 analysis types is selected in the end with  ana->SetAnalysisType(AliAnaGamma::kCorrelation);
    //The 4 correlation analysis are selected when the corresponding class is initialized

    //============================
    //First initialize prompt photon algoritm
    //==============================
    AliAnaGammaDirect *gd = new AliAnaGammaDirect();
    gd->SetMinGammaPt(1.);
    //not used in option kIsolationCut
    //gd->SetConeSize(0.5); gd->SetPtThreshold(0.); gd->SetPtSumThreshold(0.);
    gd->SetICMethod(AliAnaGammaDirect::kPtIC) ;//Options:
          //kNoIC: Accept all photons, no isolation used
          //kPtIC: IC with cut on pT
          //kSumPtIC: IC with cut on pT sum in cone
          //kSeveralIC: Not allowed in kCorrelation analysis
//     //for option kSeveralIC:
//      gd->SetNCones(2); gd->SetNPtThresholds(3);
//      gd->SetConeSizes(0,0.2); gd->SetConeSizes(1,0.3); //gd->SetConeSizes(2,0.4); 
//      gd->SetPtThresholds(0,0); gd->SetPtThresholds(1,1.); gd->SetPtThresholds(2,2);
 
    //============================
    //Second, select the correlation algoritm
    //==============================
    //Uncomment 1 of the 4 options

    //--- Option 1 ---
    //AliAnaGammaParton *gc = new AliAnaGammaParton();
    //No associated setters and getters for the moment.

    //--- Option 2 ---
    AliAnaGammaHadron *gc = new AliAnaGammaHadron();
//     gc->SetDeltaPhiCutRange(1,4); //Correlation in delta phi (gamma-particle), radians
//     gc->SetMinPtHadron(5.);
//     gc->SetJetsOnlyInCTS(kFALSE); // Don't consider particles in opposite calorimeter
    
    //--- Option 3 ---
    //    AliAnaGammaJetLeadCone *gc = new AliAnaGammaJetLeadCone();
//     gc->SetDeltaPhiCutRange(2.8,3.5); //Correlation with leading particle delta phi (gamma-particle), radians
//     gc->SetRatioCutRange(0.1,1.2);
//     gc->SetJetsOnlyInCTS(kFALSE); // Don't consider particles in opposite calorimeter
//     //and many more setters and getters

    //--- Option 4 ---
    //    AliAnaGammaJetFinder *gc = new AliAnaGammaJetFinder();
    //It does nothing for the moment

    //In case of option 1 or 2, we need to select pairs to be candidate to pi0
    //correlated with the opposite prompt photon. Need to play with this class
    AliNeutralMesonSelection *nms = new AliNeutralMesonSelection();
    nms->SetInvMassCutRange(0.1,0.17);
    nms->KeepNeutralMesonSelectionHistos(kTRUE); //keep in file several histograms or not.
    // and other parameters

    //---------------------------------------------------------------------
    // Finally: Set  analysis algorithm and reader
    //---------------------------------------------------------------------
    ana = new AliAnaGamma();
    ana->SetReader(reader);//pointer to reader
    ana->SetAnalysisType(AliAnaGamma::kCorrelation); //set kPrompt, kCorrelation
    ana->SetGammaDirect(gd);//pointer to direct photon algorithm
    ana->SetGammaCorrelation(gc);//pointer to correlation algorithm
    ana->SetNeutralMesonSelection(nms); //pi0 pair selection
    ana->SetCalorimeter("EMCAL"); //Prompt photon calorimeter
   
    //
    return ana ;
}
