#ifndef ALIPHOSCLUSTERIZER_H
#define ALIPHOSCLUSTERIZER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
                            
/* $Id$ */

//_________________________________________________________________________
//  Base class for the clusterization algorithm (pure abstract)
//*-- Author: Yves Schutz (SUBATECH)

// --- ROOT system ---

#include "TObject.h" 
#include "TClonesArray.h"

// --- Standard library ---

// --- AliRoot header files ---

#include "AliPHOSDigit.h"
#include "AliPHOSRecPoint.h"

class AliPHOSClusterizer : public TObject {

public:

  AliPHOSClusterizer() ;          // ctor            
  virtual ~AliPHOSClusterizer() ; // dtor

  virtual Float_t Calibrate(Int_t Amp) = 0 ; 
  virtual Bool_t IsInEmc(AliPHOSDigit * digit)= 0 ;   
  virtual void    GetNumberOfClustersFound(Int_t * numb) = 0 ; 
  virtual void    GetCalibrationParameters(Float_t & A, Float_t &B) = 0 ; 
  virtual Float_t GetEmcClusteringThreshold() = 0 ; 
  virtual Float_t GetEmcEnergyThreshold() = 0 ;  
  virtual Float_t GetLocalMaxCut() = 0 ; 
  virtual Float_t GetLogWeightCut() = 0 ; 
  virtual Float_t GetPpsdClusteringThreshold() = 0 ; 
  virtual Float_t GetPpsdEnergyThreshold() = 0 ; 

  virtual void  MakeClusters(const DigitsList * dl, RecPointsList * emccl, RecPointsList * ppsdl) = 0 ; 
  virtual void PrintParameters() = 0 ;  
  virtual void SetCalibrationParameters(Float_t A, Float_t B) = 0 ; 
  virtual void SetEmcClusteringThreshold(Float_t cluth) = 0 ; 
  virtual void SetEmcEnergyThreshold(Float_t enth) = 0 ;  
  virtual void SetLocalMaxCut(Float_t cut) = 0 ; 
  virtual void SetLogWeightCut(Float_t w) = 0 ; 
  virtual void SetPpsdClusteringThreshold(Float_t cluth) = 0 ; 
  virtual void SetPpsdEnergyThreshold(Float_t enth) = 0 ; 
 
  ClassDef(AliPHOSClusterizer,1)  // Clusterization algorithm class (abstract base class)

} ;

#endif // AliPHOSCLUSTERIZER_H
