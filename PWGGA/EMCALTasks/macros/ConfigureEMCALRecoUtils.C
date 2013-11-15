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
    
    reco->SetExoticCellDiffTimeCut(50);     // If |t cell max - t cell in cross| > 50 do not add its energy 
    reco->SetExoticCellFractionCut(0.97);   // 1-Ecross/Ecell > 0.97 -> out
    reco->SetExoticCellMinAmplitudeCut(4.); // 4 GeV    
  }  
  
  //Recalibration factors
  
  if(bRecalE && ! bMC)
  {
    reco->SwitchOnRecalibration();
    reco->SwitchOnRunDepCorrection();    
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
    if(!bMC)
    {
      printf("xxx SET Non linearity correction kBeamTestCorrected xxx\n");
      reco->SetNonLinearityFunction(AliEMCALRecoUtils::kBeamTestCorrectedv2);
    }
    else
    {       
      printf("xxx SET Non linearity correction kPi0MCv3 xxx\n");
      reco->SetNonLinearityFunction(AliEMCALRecoUtils::kPi0MCv3);
    }
  }
  else 
  {
    printf("xxx DON'T SET Non linearity correction xxx\n");
    reco->SetNonLinearityFunction(AliEMCALRecoUtils::kNoCorrection);
  }
  
}
