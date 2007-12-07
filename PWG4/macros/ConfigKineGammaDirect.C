/* $Id$ */
/* $Log$ */

//------------------------------------
// Configuration macro example:
//
// Do prompt photon analysis with kinematics
//
// Author : Gustavo Conesa Balbastre (INFN-LNF)
//------------------------------------
AliAnaGamma*  ConfigGammaAnalysis()
{
    //
    // Configuration goes here
    // 
    printf("ConfigGammaAnalysis() \n");

    //-----------------------------------------------------------  
    // Reader
    //-----------------------------------------------------------

    AliGammaMCReader *reader = new AliGammaMCReader();
    // Switch "on" or "off"  detectors used in study
    reader->SwitchOnEMCAL(kTRUE);
    reader->SwitchOnPHOS(kFALSE) ;
    reader->SwitchOnCTS(kTRUE) ;
    
    //Set detectors acceptance (in real and montecarlo data)
    reader->SetCTSEtaCut(1); reader->SetEMCALEtaCut(0.7); //reader->SetPHOSEtaCut(0.2);
    reader->SetPhiEMCALCut(40*TMath::DegToRad(), 200*TMath::DegToRad());     
    //     reader->SetPhiPHOSCut(200*TMath::DegToRad(), 350*TMath::DegToRad()); 
    
    //Set minimum pt for particles in analysis  (in real and montecarlo data)
    reader->SetNeutralPtCut(0.5); reader->SetChargedPtCut(0.3); 
    
    reader->SetDecayPi0Flag(AliGammaMCReader::kGeantDecay) ; //Options
    //kNoDecay :Do not decay pi0, keep them in the list as they are
    //kGeantDecay: Look for gamma decayed by GEANT
    //kDecay: Decay pi0 by hand (geant was not used)
    //kDecayGamma: Pi0 is decayed by PYTHIA, pi0 is not final check if photons overlapp
    
    //parameters to study if decay is overlapped:
    reader->SetEMCALIPDistance(450.); reader->SetPHOSIPDistance(460.);
    reader->SetEMCALMinAngle(2.5 * TMath::DegToRad() );    
    reader->SetPHOSMinAngle(0.45 * TMath::DegToRad() ); //Minimum overlapp distance
    reader->SetCheckOverlapping(kTRUE);
    
    //============================
    //Prompt photon algoritm
    //==============================
    AliAnaGammaDirect *gd = new AliAnaGammaDirect();
    gd->SetMinGammaPt(5.);
    gd->SetConeSize(0.5); gd->SetPtThreshold(1); gd->SetPtSumThreshold(.);
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
    ana->SetCalorimeter("EMCAL"); //Prompt photon calorimeter
   
    //
    return ana ;
}
