/* $Id$ */
/* $Log$ */

//------------------------------------
// Configuration macro example:
//
// Do prompt photon analysis with ESDs, check with kinematics
//
// Author : Gustavo Conesa Balbastre (INFN-LNF)
//------------------------------------

AliAnaGamma*  ConfigGammaAnalysis()
{
    //
    // Configuration goes here
    // 
    printf("ConfigGammaAnalysis() \n");

    //----------------------------------------------------------
    //-----------------------------------------------------------  
    // Define Reader
    //-----------------------------------------------------------
    //----------------------------------------------------------

    AliGammaMCDataReader *reader = new AliGammaMCDataReader();
    // Switch "on" or "off"  detectors used in study
    reader->SwitchOnEMCAL(kTRUE);
    reader->SwitchOnPHOS(kFALSE) ;
    reader->SwitchOnCTS(kTRUE) ;
    // Set detectors acceptance (in real and montecarlo data)
//     reader->SetCTSEtaCut(1); reader->SetEMCALEtaCut(1); reader->SetPHOSEtaCut(0.2);
    reader->SetPhiEMCALCut(60*TMath::DegToRad(), 180*TMath::DegToRad());     
    //   reader->SetPhiPHOSCut(260*TMath::DegToRad(), 320*TMath::DegToRad()); //Example, only modules in  TRD/TOF holes 
    // Set minimum pt for particles in analysis  (in real and montecarlo data)
    //  reader->SetNeutralPtCut(0.5); reader->SetChargedPtCut(0.2); 
    //pid of measured particles
    reader->SetEMCALPIDOn(kFALSE); reader->SetPHOSPIDOn(kTRUE); //No pid, accept all particles 
    //if previous kTrue
    // use selection with simple weights
//     reader->SetPHOSPhotonWeight(0.7);    reader->SetPHOSPi0Weight(0.7); 
//     reader->SetEMCALPhotonWeight(0.7);    reader->SetEMCALPi0Weight(0.7);
    // use more complicated selectionm, particle weight depending on cluster energy
    reader->UsePHOSPIDWeightFormula(kTRUE);
    TFormula * photonF = new TFormula("photonWeight","0.98*(x<40)+ 0.68*(x>=100)+(x>=40 && x<100)*(0.98+x*(6e-3)-x*x*(2e-04)+x*x*x*(1.1e-06))");
    TFormula * pi0F = new TFormula("pi0Weight","0.98*(x<65)+ 0.915*(x>=100)+(x>=65 && x-x*(1.95e-3)-x*x*(4.31e-05)+x*x*x*(3.61e-07))");
    reader->SetPHOSPhotonWeightFormula(photonF);
    reader->SetPHOSPi0WeightFormula(pi0F);

    //parameters to study if decay is overlapped or there was a conversion before the calo:
    reader->SetEMCALIPDistance(450.); reader->SetPHOSIPDistance(460.);
    reader->SetEMCALMinAngle(2.5 * TMath::DegToRad() );    
    reader->SetPHOSMinAngle(0.45 * TMath::DegToRad() ); //Minimum overlapp distance

    //============================
    // Initialize prompt photon algoritm
    //==============================
    AliAnaGammaDirect *gd = new AliAnaGammaDirect();
    gd->SetMinGammaPt(5.);
    gd->SetConeSize(0.5); gd->SetPtThreshold(1.); gd->SetPtSumThreshold(1.);
    gd->SetICMethod(AliAnaGammaDirect::kPtIC) ;//Options:
          //kNoIC: Accept all photons, no isolation used
          //kPtIC: IC with cut on pT
          //kSumPtIC: IC with cut on pT sum in cone

    //---------------------------------------------------------------------
    // Finally: Set  analysis algorithm and reader
    //---------------------------------------------------------------------
    ana = new AliAnaGamma();
    ana->SetReader(reader);//pointer to reader
    ana->SetAnalysisType(AliAnaGamma::kPrompt); //set kPrompt, kCorrelation
    ana->SetGammaDirect(gd);//pointer to direct photon algorithm
    ana->SetCalorimeter("EMCAL"); //Prompt photon calorimeter, PHOS or EMCAL
   
    //
    return ana ;
}
