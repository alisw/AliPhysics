#ifndef ALIPHOSCLUSTERIZERV1_H
#define ALIPHOSCLUSTERIZERV1_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//_________________________________________________________________________
//  Implementation version 1 of the clusterization algorithm                     
//
//*-- Author: Yves Schutz (SUBATECH)

// --- ROOT system ---

// --- Standard library ---

// --- AliRoot header files ---

#include "AliPHOSClusterizer.h"
#include "AliPHOSDigit.h" 
#include "AliPHOSGeometry.h" 



class AliPHOSClusterizerv1 : public AliPHOSClusterizer {
  
public:
  
  AliPHOSClusterizerv1() ;             // ctor            
  virtual ~AliPHOSClusterizerv1(){} ;  // dtor
  
  Int_t AreNeighbours(AliPHOSDigit * d1, AliPHOSDigit * d2) ; // Checks if digits are in neighbour cells 
  Float_t Calibrate(Int_t Amp){ return (fA + fB * Amp) ;}     // Tranforms Amp to energy 
  void FillandSort(const DigitsList * dl, TObjArray * tl) ;   // Sorts the list according to increasing id
  virtual void GetNumberOfClustersFound(Int_t * numb) ; 
 
  virtual void GetCalibrationParameters(Float_t & A, Float_t &B) { A = fA; B = fB; } 
  virtual Float_t GetEmcClusteringThreshold()  { return fEmcClusteringThreshold;}
  virtual Float_t GetEmcEnergyThreshold()      { return fEmcEnergyThreshold; }  
  virtual Float_t GetLocalMaxCut()             { return fLocMaxCut;} 
  virtual Float_t GetLogWeightCut()            { return fW0;}  
  virtual Float_t GetLocalMaxCutCPV()          { return fLocMaxCutCPV;} 
  virtual Float_t GetLogWeightCutCPV()         { return fW0CPV;}  
  virtual Float_t GetPpsdClusteringThreshold() { return fPpsdClusteringThreshold;  } 
  virtual Float_t GetPpsdEnergyThreshold()     { return fPpsdEnergyThreshold;  }
  virtual Float_t GetCpvClusteringThreshold()  { return fCpvClusteringThreshold;  } 
  virtual Float_t GetCpvEnergyThreshold()      { return fCpvEnergyThreshold;  }

  virtual Bool_t IsInEmc (AliPHOSDigit * digit) ;                      // Tells if id digit is in EMC
  virtual Bool_t IsInPpsd(AliPHOSDigit * digit) ;                      // Tells if id digit is in PPSD
  virtual Bool_t IsInCpv (AliPHOSDigit * digit) ;                      // Tells if id digit is in CPV
  virtual void MakeClusters(const DigitsList * dl, 
			    AliPHOSRecPoint::RecPointsList * emcl, 
			    AliPHOSRecPoint::RecPointsList * ppsdl) ; // does the job 
  virtual void PrintParameters() ;  
  virtual void SetCalibrationParameters(Float_t A,Float_t B){ fA = A ; fB = B;} 
  virtual void SetEmcClusteringThreshold(Float_t cluth)  { fEmcClusteringThreshold = cluth ; }
  virtual void SetEmcEnergyThreshold(Float_t enth)       { fEmcEnergyThreshold = enth ; } 
  virtual void SetLocalMaxCut(Float_t cut)               { fLocMaxCut = cut ; }
  virtual void SetLogWeightCut(Float_t w)                { fW0 = w ; }
  virtual void SetLocalMaxCutCPV(Float_t cut)            { fLocMaxCutCPV = cut ; }
  virtual void SetLogWeightCutCPV(Float_t w)             { fW0CPV = w ; }
  virtual void SetPpsdClusteringThreshold(Float_t cluth) { fPpsdClusteringThreshold = cluth ; }
  virtual void SetPpsdEnergyThreshold(Float_t enth)      { fPpsdEnergyThreshold = enth ; } 
  virtual void SetCpvClusteringThreshold(Float_t cluth)  { fCpvClusteringThreshold = cluth ; }
  virtual void SetCpvEnergyThreshold(Float_t enth)       { fCpvEnergyThreshold = enth ; } 
  
private:
  
  Float_t fA ;                       // offset of the energy calibration
  Float_t fB ;                       // gain of the energy calibration
  AliPHOSGeometry * fGeom ;          // pointer to geometry
  Int_t   fNumberOfEmcClusters ;     // number of EMC clusters found 
  Float_t fEmcClusteringThreshold ;  // minimum energy to include a EMC digit in a cluster
  Float_t fEmcEnergyThreshold ;      // minimum energy of EMC digit to be considered
  Int_t   fNumberOfPpsdClusters ;    // number of PPSD clusters found
  Float_t fPpsdClusteringThreshold ; // minimum energy to include a PPSD digit in a cluster
  Float_t fPpsdEnergyThreshold ;     // minimum energy of PPSD digit to be considered
  Int_t   fNumberOfCpvClusters ;     // number of CPV clusters found
  Float_t fCpvClusteringThreshold ;  // minimum energy to include a CPV digit in a cluster
  Float_t fCpvEnergyThreshold ;      // minimum energy of CPV digit to be considered
  Float_t fLocMaxCut ;               // minimum energy difference to distinguish local maxima in a cluster
  Float_t fW0 ;                      // logarithmic weight for the cluster center of gravity calculation
  Float_t fLocMaxCutCPV ;            // minimum energy difference to distinguish local maxima in a CPV cluster
  Float_t fW0CPV ;                   // logarithmic weight for the CPV cluster center of gravity calculation
    
  ClassDef(AliPHOSClusterizerv1,1)   // Clusterizer implementation version 1

};

#endif // AliPHOSCLUSTERIZERV1_H
