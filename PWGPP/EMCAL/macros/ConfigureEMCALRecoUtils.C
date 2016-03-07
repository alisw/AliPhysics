///
/// \file ConfigureEMCALRecoUtils.C
/// \brief Configuration of AliEMCALRecoUtils.
///
/// Example of configuration of AliEMCALRecoUtils. Called in different analysis configuration macros.
/// This class is used to calibrate/correct/accept EMCal clusters.
///
/// The input parameters:
/// \param reco: pointer to object to initialize in this macro.
/// \param bMC: Bool, indicates if data is MC.
/// \param bExotic: Bool, indicates if exotic clusters are removed.
/// \param bNonLin: Bool, indicates if non linearity correction is applied on clusters.
/// \param bRecalE: Bool, indicates if energy recalibration is applied.
/// \param bBad: Bool, indicates if bad channels/clusters are removed.
/// \param bRecalT: Bool, indicates if time is calibrated.
///
/// \author : Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, (LPSC-CNRS)
///
void ConfigureEMCALRecoUtils(AliEMCALRecoUtils* reco,
                             Bool_t  bMC    = kFALSE,
                             Bool_t  bExotic= kTRUE,
                             Bool_t  bNonLin= kFALSE,
                             Bool_t  bRecalE= kTRUE,
                             Bool_t  bBad   = kTRUE,
                             Bool_t  bRecalT= kTRUE)
{
  printf("**** Configure AliEMCALRecoUtils ***\n");
  
  // Exotic cells removal
  
  if(bExotic)
  {
    printf("Remove exotics in EMCAL\n");
    reco->SwitchOnRejectExoticCell() ;
    reco->SwitchOnRejectExoticCluster(); 
    
//  reco->SetExoticCellDiffTimeCut(50);     // If |t cell max - t cell in cross| > 50 do not add its energy, avoid 
    reco->SetExoticCellFractionCut(0.97);   // 1-Ecross/Ecell > 0.97 -> out
    reco->SetExoticCellMinAmplitudeCut(4.); // 4 GeV    
  }  
  
  // Recalibration factors
  
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
    
  // Recalculate position with method
    
  reco->SetPositionAlgorithm(AliEMCALRecoUtils::kPosTowerGlobal);   

  // Non linearity
  
  if( bNonLin ) 
  { 
    if(!bMC)
    {
      printf("xxx SET Non linearity correction kBeamTestCorrected xxx\n");
      reco->SetNonLinearityFunction(AliEMCALRecoUtils::kBeamTestCorrectedv3);
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
