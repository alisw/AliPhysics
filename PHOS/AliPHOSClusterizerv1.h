#ifndef ALIPHOSCLUSTERIZERV1_H
#define ALIPHOSCLUSTERIZERV1_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////
//  Clusterizer implementation version 1      //
//  algorithme class                          //
//                                            //
//  Author Yves Schutz     SUBATECH           //
//                                            //  
//                                            //
////////////////////////////////////////////////

// --- ROOT system ---

// --- Standard library ---

// --- AliRoot header files ---

#include "AliPHOSClusterizer.h"
#include "AliPHOSDigit.h" 



class AliPHOSClusterizerv1 : public AliPHOSClusterizer {
  
public:
  
  AliPHOSClusterizerv1() ;             // ctor            
  virtual ~AliPHOSClusterizerv1(){} ;  // dtor
  
  Int_t AreNeighbours(AliPHOSDigit * d1, AliPHOSDigit * d2) ; // Checks if digits are in neighbour cells 
  Float_t Calibrate(Int_t Amp){ return (fA + fB * Amp) ;}     // Tranforms Amp to energy 
  void FillandSort(const DigitsList * dl, TObjArray * tl) ;   // Sorts the list according to increasing id
  Float_t GetLogWeightCut(void){return  fW0 ; }
  Float_t GetLocalMaxCut(void) {return  fLocMaxCut ; }
  virtual void GetNumberOfClustersFound(Int_t * numb) ;   
  Bool_t IsInEmc(AliPHOSDigit * digit) ;                      // Tells if id digit is in EMC
  virtual void MakeClusters(const DigitsList * dl, RecPointsList * emcl, RecPointsList * ppsdl) ; // does the job 
  void PrintParameters() ;  
  void SetCalibrationParameters(Float_t A,Float_t B){ fA = A ; fB = B;} 
  void SetEmcClusteringThreshold(Float_t cluth) { fEmcClusteringThreshold = cluth ; }
  void SetEmcEnergyThreshold(Float_t enth) { fEmcEnergyThreshold = enth ; } 
  void SetLocalMaxCut(Float_t cut) { fLocMaxCut = cut ; }
  void SetLogWeightCut(Float_t w) { fW0 = w ; }
  void SetPpsdClusteringThreshold(Float_t cluth) { fPpsdClusteringThreshold = cluth ; }
  void SetPpsdEnergyThreshold(Float_t enth) { fPpsdEnergyThreshold = enth ; } 
  
private:
  
  Float_t fA ;
  Float_t fB ;
  Float_t fEmcClusteringThreshold ; 
  Float_t fEmcEnergyThreshold ; 
  Float_t fLocMaxCut ;   
  Int_t   fNumberOfEmcClusters ; 
  Int_t   fNumberOfPpsdClusters ; 
  Float_t fPpsdClusteringThreshold ; 
  Float_t fPpsdEnergyThreshold ;  
  Float_t fW0 ;   
  
public: 
  
  ClassDef(AliPHOSClusterizerv1,1)  // Clusterizer implementation , version 1

};

#endif // AliPHOSCLUSTERIZERV1_H
