// $Id$

void ConfigureEMCALRecoUtils(AliEMCALRecoUtils* reco,
                             Bool_t  bMC    = kFALSE,
                             Bool_t  bExotic= kTRUE,
                             Bool_t  bNonLin= kFALSE,
                             Bool_t  bRecalE= kTRUE,
                             Bool_t  bBad   = kTRUE,
                             Bool_t  bRecalT= kTRUE)
{  

  // Configure RecoUtils with OADB objects
  
  printf("**** Configure AliEMCALRecoUtils ***\n");
  
  // Exotic cells removal
  
  if(bExotic)
  {
    printf("Remove exotics in EMCAL\n");
    reco->SwitchOnRejectExoticCell() ;
    reco->SwitchOnRejectExoticCluster(); 
    
    reco->SetExoticCellDiffTimeCut(10000);    // Open  
    reco->SetExoticCellFractionCut(0.95);     // 1-Ecross/Ecell > 0.95 -> out
    reco->SetExoticCellMinAmplitudeCut(0.75); // 750 MeV    
  }  
  
  //Recalibration factors
  
  if(bRecalE && ! bMC)
  {
    reco->SwitchOnRecalibration();
  } 

  // Remove EMCAL hot channels 
  
  if(bBad)
  {
    reco->SwitchOnBadChannelsRemoval();
    reco->SwitchOnDistToBadChannelRecalculation();
  }
 
  // *** Time recalibration settings ***
  
  if(bRecalT && ! bMC)
  {
    reco->SwitchOnTimeRecalibration();
  }
    
  // position
    
  reco->SetPositionAlgorithm(AliEMCALRecoUtils::kPosTowerGlobal);   

  // Non linearity
  
  if( bNonLin ) 
  { 
    if(!kSimulation) reco->SetNonLinearityFunction(AliEMCALRecoUtils::kBeamTestCorrected);
    else             reco->SetNonLinearityFunction(AliEMCALRecoUtils::kPi0MC);
  }
}
