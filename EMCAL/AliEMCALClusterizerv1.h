#ifndef ALIEMCALCLUSTERIZERV1_H
#define ALIEMCALCLUSTERIZERV1_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//_________________________________________________________________________
//  Implementation version 1 of the clusterization algorithm                     
//  Performs clusterization (collects neighbouring active cells) and 
//  unfolding of the clusters with several local maxima.  
//  results are stored in TreeR#, branches PHOSEmcRP (EMC recPoints),
//  PHOSCpvRP (CPV RecPoints) and AliPHOSClusterizer
//
//*-- Author: Yves Schutz (SUBATECH)
// Modif: 
//  August 2002 Yves Schutz: clone PHOS as closely as possible and intoduction
//                           of new  IO (à la PHOS)

// --- ROOT system ---

// --- Standard library ---

// --- AliRoot header files ---

#include "AliEMCALClusterizer.h"
class AliEMCALTowerRecPoint ; 
class AliEMCALDigit ;
class AliEMCALDigitizer ;
class AliEMCALGeometry ;


class AliEMCALClusterizerv1 : public AliEMCALClusterizer {
  
public:
  
  AliEMCALClusterizerv1() ;         
  AliEMCALClusterizerv1(const char * headerFile, const char * name = "Default", const Bool_t toSplit=kFALSE);
  virtual ~AliEMCALClusterizerv1()  ;
  
  virtual Int_t   AreNeighbours(AliEMCALDigit * d1, AliEMCALDigit * d2)const ; 
                               // Checks if digits are in neighbour cells 

  virtual Float_t Calibrate(Int_t amp, Int_t where)const ;  // Tranforms Amp to energy 

  virtual void    GetNumberOfClustersFound(int * numb )const{ numb[0] = fNumberOfPREClusters ; 
                                                              numb[1] = fNumberOfECClusters ; 
                                                              numb[2] = fNumberOfHCClusters ; }

  virtual Float_t GetPREClusteringThreshold()const{ return fPREClusteringThreshold;  } 
  virtual Float_t GetECClusteringThreshold()const{ return fECClusteringThreshold;}
  virtual Float_t GetHCClusteringThreshold()const{ return fHCClusteringThreshold;}

  virtual Float_t GetPRELocalMaxCut()const       { return fPRELocMaxCut;} 
  virtual Float_t GetPREShoLogWeight()const      { return fPREW0;}  
  virtual Float_t GetECLocalMaxCut()const        { return fECLocMaxCut;} 
  virtual Float_t GetECLogWeight()const          { return fECW0;}  
  virtual Float_t GetHCLocalMaxCut()const        { return fHCLocMaxCut;} 
  virtual Float_t GetHCLogWeight()const          { return fHCW0;}  

  virtual Float_t GetTimeGate() const            { return fTimeGate ; }
  virtual const char *  GetRecPointsBranch() const{ return GetName() ;}
  virtual const Int_t GetRecPointsInRun() const   {return fRecPointsInRun ;} 

  void    Exec(Option_t *option);                // Does the job

  virtual void Print(Option_t * option)const ;

  virtual void SetECClusteringThreshold(Float_t cluth)  { fECClusteringThreshold = cluth ; }
  virtual void SetECLocalMaxCut(Float_t cut)            { fECLocMaxCut = cut ; }
  virtual void SetECLogWeight(Float_t w)                { fECW0 = w ; }
  virtual void SetHCClusteringThreshold(Float_t cluth)  { fHCClusteringThreshold = cluth ; }
  virtual void SetHCLocalMaxCut(Float_t cut)            { fHCLocMaxCut = cut ; }
  virtual void SetHCLogWeight(Float_t w)                { fHCW0 = w ; }
  virtual void SetTimeGate(Float_t gate)                { fTimeGate = gate ;}
  virtual void SetPREClusteringThreshold(Float_t cluth) { fPREClusteringThreshold = cluth ; }
  virtual void SetPRELocalMaxCut(Float_t cut)           { fPRELocMaxCut = cut ; }
  virtual void SetPRELogWeight(Float_t w)               { fPREW0 = w ; }
  virtual void SetUnfolding(Bool_t toUnfold = kTRUE )      {fToUnfold = toUnfold ;}  
  static Double_t ShowerShape(Double_t r) ; // Shape of EM shower used in unfolding; 
                                            //class member function (not object member function)
  static void UnfoldingChiSquare(Int_t & nPar, Double_t * Grad, Double_t & fret, Double_t * x, Int_t iflag)  ;
                                            // Chi^2 of the fit. Should be static to be passes to MINUIT
  virtual const char * Version() const { return "clu-v1" ; }  

protected:

  void           WriteRecPoints(Int_t event) ;
  virtual void   MakeClusters( ) ;            
  
private:

  const TString BranchName() const ; 
  void    GetCalibrationParameters(void) ;
  
  Bool_t  FindFit(AliEMCALTowerRecPoint * emcRP, AliEMCALDigit ** MaxAt, Float_t * maxAtEnergy, 
		  Int_t NPar, Float_t * FitParametres) const; //Used in UnfoldClusters, calls TMinuit
  void Init() ;
  void InitParameters() ;

  virtual void   MakeUnfolding() ;
  void           UnfoldCluster(AliEMCALTowerRecPoint * iniEmc,Int_t Nmax, 
		       AliEMCALDigit ** maxAt,Float_t * maxAtEnergy ) ; //Unfolds cluster using TMinuit package
  void           PrintRecPoints(Option_t * option) ;

private:

  Bool_t  fDefaultInit;              //! Says if the task was created by defaut ctor (only parameters are initialized)

  Int_t   fNTowers ;                 // number of Towers in EMCAL

  Bool_t  fToUnfold ;                // To perform unfolding 

  Int_t   fNumberOfPREClusters ;     // number of clusters found in PRE section 
  Int_t   fNumberOfECClusters ;      // number of clusters found in EC section
  Int_t   fNumberOfHCClusters ;      // number of clusters found in HC section
  
  //Calibration parameters... to be replaced by database 
  Float_t fADCchannelPRE ;          // width of one ADC channel for PRE section (GeV)
  Float_t fADCpedestalPRE ;         // pedestal of ADC for PRE section (GeV)
  Float_t fADCchannelEC ;           // width of one ADC channel for EC section (GeV)
  Float_t fADCpedestalEC ;          // pedestal of ADC for EC section (GeV) 
  Float_t fADCchannelHC ;           // width of one ADC channel for HC section (GeV)
  Float_t fADCpedestalHC ;          // pedestal of ADC for HC section (GeV) 
 
  Float_t fECClusteringThreshold ;  // minimum energy to include a EC digit in a cluster
  Float_t fHCClusteringThreshold ;  // minimum energy to include a HC digit in a cluster
  Float_t fPREClusteringThreshold ; // minimum energy to include a PRE digit in a cluster
  Float_t fECLocMaxCut ;            // minimum energy difference to distinguish local maxima in a cluster
  Float_t fECW0 ;                   // logarithmic weight for the cluster center of gravity calculation
  Float_t fHCLocMaxCut ;            // minimum energy difference to distinguish local maxima in a cluster
  Float_t fHCW0 ;                   // logarithmic weight for the cluster center of gravity calculation
  Float_t fPRELocMaxCut ;           //  minimum energy difference to distinguish local maxima in a CPV cluster
  Float_t fPREW0 ;                  // logarithmic weight for the CPV cluster center of gravity calculation
  Int_t fRecPointsInRun ;           //! Total number of recpoints in one run
  Float_t fTimeGate ;               // Maximum time difference between the digits in ont EMC cluster
    
  ClassDef(AliEMCALClusterizerv1,2)   // Clusterizer implementation version 1

};

#endif // AliEMCALCLUSTERIZERV1_H
