// YS Subatech Mai 2002

//Root
#include "TString.h"

//AliRoot
#include "PHOS/AliPHOSSDigitizer.h"
#include "PHOS/AliPHOSDigitizer.h"
#include "EMCAL/AliEMCALSDigitizer.h"
#include "EMCAL/AliEMCALDigitizer.h"

void Hits2SDigits( Bool_t split=kFALSE, TString fileName = "galice.root") {

  // usage : 
  // 1. write SDigits in the same file as Hits --------------- (OK)
  //root [0] .L Reconstruct.C++
  //root [1] SDigits2Digits()
  // 2. write SDigits in a separate file, one per detector, from Hits --------------- (OK)
  //root [0] .L Reconstruct.C++
  //root [1] SDigits2Digits(kTRUE) // SDigits saved in [DET}.SDigits.root (DET=PHOS, EMCAL)

  AliPHOSSDigitizer * sdp = new AliPHOSSDigitizer(fileName) ; 
  if (split) 
    sdp->SetSplitFile() ;
  sdp->ExecuteTask("deb") ; 

  delete sdp ;

  AliEMCALSDigitizer * sde = new AliEMCALSDigitizer(fileName) ; 
  if (split) 
    sde->SetSplitFile() ;
  sde->ExecuteTask("deb") ; 
 
  delete sde ; 

}

//________________________________________________________________________
void SDigits2Digits( Bool_t split=kFALSE, TString fileName = "galice.root") {
  
 // usage : 
  // 1. write SDigits in the same file as Hits --------------- (OK)
  //root [0] .L Reconstruct.C++
  //root [1] Hits2Digits()
  // 2. write SDigits in a separate file, one per detector, from Hits --------------- (OK)
  //root [0] .L Reconstruct.C++
  //root [1] Hits2Digits(kTRUE) // Digits saved in [DET}.Digits.root (DET=PHOS, EMCAL)

  // PHOS
  AliPHOSDigitizer * dp ; 
 
  if (split) {
    dp = new AliPHOSDigitizer("PHOS.SDigits.root") ; 
    dp->SetSplitFile() ; } 
  else 
    dp = new AliPHOSDigitizer(fileName) ; 
  
  dp->ExecuteTask("deb") ; 
  
  delete dp ;

  //EMCAL
  AliEMCALDigitizer * de ; 

  if (split) {
    de = new AliEMCALDigitizer("EMCAL.SDigits.root") ;
    de->SetSplitFile() ;
  } else 
    de = new AliEMCALDigitizer(fileName) ; 
  
  de->ExecuteTask("deb") ; 
  
  delete de ; 
}

//________________________________________________________________________
void Hits2Digits (Bool_t split=kFALSE, TString fileName = "galice.root") {
  // usage : 
  // 1. write (S)Digits in the same file as Hits --------------- (OK)
  //root [0] .L Reconstruct.C++
  //root [1] Hits2Digits()
  // 2. write (S)Digits in a separate file, one per detector, from Hits --------------- (OK)
  //root [0] .L Reconstruct.C++
  //root [1] Hits2Digits(kTRUE) // SDigits saved in [DET}.SDigits.root (DET=PHOS, EMCAL)
                                // Digits  saved in [DET}.Digits.root  (DET=PHOS, EMCAL)

  //PHOS
  AliPHOSSDigitizer * sdp = new AliPHOSSDigitizer(fileName) ; 
  if (split) 
    sdp->SetSplitFile() ;
  sdp->ExecuteTask("deb") ; 

  delete sdp ;

  AliPHOSDigitizer * dp ; 
 
  if (split) {
    dp = new AliPHOSDigitizer("PHOS.SDigits.root") ; 
    dp->SetSplitFile() ; } 
  else 
    dp = new AliPHOSDigitizer(fileName) ; 
  
  dp->ExecuteTask("deb") ; 
  
  delete dp ;

  //EMCAL
  AliEMCALSDigitizer * sde = new AliEMCALSDigitizer(fileName) ; 
  if (split) 
    sde->SetSplitFile() ;
  sde->ExecuteTask("deb") ; 
  
  delete sde ; 
  
  AliEMCALDigitizer * de ; 
  if (split) {
    de = new AliEMCALDigitizer("EMCAL.SDigits.root") ;
    de->SetSplitFile() ;
  } else 
    de = new AliEMCALDigitizer(fileName) ; 
  
  de->ExecuteTask("deb") ; 
  
  delete de ; 
}

 
