// YS Subatech Mai 2002
// YK Subatech 7 Aug 2002

// EMCAL Reconstruction chain:
// Hits -> SDigits -> Digits -> RecPoints

//Root
#include "TString.h"

//AliRoot
#include "STEER/AliRun.h"
#include "EMCAL/AliEMCALSDigitizer.h"
#include "EMCAL/AliEMCALDigitizer.h"
#include "EMCAL/AliEMCALClusterizerv1.h"

void EMCALHits2SDigits( Bool_t split=kFALSE, TString fileName = "galice.root") {

  // usage : 
  // 1. write SDigits in the same file as Hits --------------- (OK)
  //root [0] .L Reconstruct.C++
  //root [1] SDigits2Digits()
  // 2. write SDigits in a separate file, one per detector, from Hits --------------- (OK)
  //root [0] .L Reconstruct.C++
  //root [1] SDigits2Digits(kTRUE) // SDigits saved in [DET}.SDigits.root (DET=PHOS, EMCAL)

  delete gAlice ; 
  gAlice = 0 ; 
  
  AliEMCALSDigitizer * sde = new AliEMCALSDigitizer(fileName) ; 
  if (split) 
    sde->SetSplitFile() ;
  sde->ExecuteTask("deb") ; 
  
  delete sde ; 

}

//________________________________________________________________________
void EMCALSDigits2Digits( Bool_t split=kFALSE, TString fileName = "galice.root") {
  
 // usage : 
  // 1. write SDigits in the same file as SDigits --------------- (OK)
  //root [0] .L Reconstruct.C++
  //root [1] SDigits2Digits()
  // 2. write SDigits in a separate file, one per detector, from SDigits --------------- (OK)
  //root [0] .L Reconstruct.C++
  //root [1] SDigitsDigits(kTRUE) // Digits saved in [DET}.Digits.root (DET=PHOS, EMCAL)

  delete gAlice ; 
  gAlice = 0 ; 
  
  AliEMCALDigitizer * de = 0 ; 

  if (split) {
    de = new AliEMCALDigitizer("EMCAL.SDigits.root") ;
    de->SetSplitFile() ;
  } else 
    de = new AliEMCALDigitizer(fileName) ; 
  
  de->ExecuteTask("deb") ; 
  
  delete de ; 
}

//________________________________________________________________________
void EMCALDigits2RecPoints( Bool_t split=kFALSE, TString fileName = "galice.root") {
  
 // usage : 
  // 1. write RecPoints in the same file as Digits --------------- OK 
  //root [0] .L Reconstruct.C++
  //root [1] Digits2RecPoints()
  // 2. write RecPoints in a separate file, one per detector, from Digits --------------- OK 
  //root [0] .L Reconstruct.C++
  //root [1] Digits2RecPoints(kTRUE) // RecPoints saved in [DET}.RecPoints.root (DET=PHOS, EMCAL)

  delete gAlice ; 
  gAlice = 0 ; 
 
  AliEMCALClusterizerv1 * ce = 0 ;  

  if (split) {
    ce = new AliEMCALClusterizerv1("EMCAL.Digits.root") ;
    ce->SetSplitFile() ;
  } else 
    ce = new AliEMCALClusterizerv1(fileName) ; 
  
  ce->ExecuteTask("deb") ; 
  
  delete ce ; 
}

//________________________________________________________________________
void EMCALHits2Digits (Bool_t split=kFALSE, TString fileName = "galice.root") {
  // usage : 
  // 1. write (S)Digits in the same file as Hits --------------- (OK)
  //root [0] .L Reconstruct.C++
  //root [1] Hits2Digits()
  // 2. write (S)Digits in a separate file, one per detector, from Hits --------------- (OK)
  //root [0] .L Reconstruct.C++
  //root [1] Hits2Digits(kTRUE) // SDigits saved in [DET}.SDigits.root (DET=PHOS, EMCAL)
                                // Digits  saved in [DET}.Digits.root  (DET=PHOS, EMCAL)

  delete gAlice ; 
  gAlice = 0 ; 
  
  //EMCAL
  AliEMCALSDigitizer * sde = new AliEMCALSDigitizer(fileName) ; 
  if (split) 
    sde->SetSplitFile() ;
  sde->ExecuteTask("deb") ; 
  
  delete sde ; 
  
  AliEMCALDigitizer * de = 0 ; 
  if (split) {
    de = new AliEMCALDigitizer("EMCAL.SDigits.root") ;
    de->SetSplitFile() ;
  } else 
    de = new AliEMCALDigitizer(fileName) ; 
  
  de->ExecuteTask("deb") ; 
  
  delete de ; 

}


 
