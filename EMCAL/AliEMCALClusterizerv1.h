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
  AliEMCALClusterizerv1(const char * headerFile, const char * name = "Default");
  virtual ~AliEMCALClusterizerv1()  ;
  
  virtual Int_t   AreNeighbours(AliEMCALDigit * d1, AliEMCALDigit * d2)const ; 
                               // Checks if digits are in neighbour cells 

  const TString BranchName() const ; 
  virtual Float_t Calibrate(Int_t amp, Bool_t inpresho)const ;  // Tranforms Amp to energy 

  virtual void    GetNumberOfClustersFound(int * numb )const{  numb[0] = fNumberOfTowerClusters ; 
                                                               numb[1] = fNumberOfPreShoClusters ; }

  virtual Float_t GetEmcClusteringThreshold()const{ return fTowerClusteringThreshold;}
  virtual Float_t GetEmcLocalMaxCut()const        { return fTowerLocMaxCut;} 
  virtual Float_t GetEmcLogWeight()const          { return fW0;}  
  virtual Float_t GetTimeGate() const             { return fTimeGate ; }
  virtual Float_t GetCpvClusteringThreshold()const{ return fPreShoClusteringThreshold;  } 
  virtual Float_t GetCpvLocalMaxCut()const        { return fPreShoLocMaxCut;} 
  virtual Float_t GetCpvLogWeight()const          { return fW0CPV;}  
  virtual char *  GetRecPointsBranch() const      { return (char*) fRecPointsBranchTitle.Data() ;}
  virtual const Int_t GetRecPointsInRun() const  {return fRecPointsInRun ;} 
  virtual char *  GetDigitsBranch() const         { return (char*) fDigitsBranchTitle.Data() ;}

  void    Exec(Option_t *option);                // Does the job

  virtual void Print(Option_t * option)const ;

  virtual void SetTowerClusteringThreshold(Float_t cluth)  { fTowerClusteringThreshold = cluth ; }
  virtual void SetTowerLocalMaxCut(Float_t cut)            { fTowerLocMaxCut = cut ; }
  virtual void SetTowerLogWeight(Float_t w)                { fW0 = w ; }
  virtual void SetTimeGate(Float_t gate)                   { fTimeGate = gate ;}
  virtual void SetPreShoClusteringThreshold(Float_t cluth) { fPreShoClusteringThreshold = cluth ; }
  virtual void SetPreShoLocalMaxCut(Float_t cut)           { fPreShoLocMaxCut = cut ; }
  virtual void SetPreShoLogWeight(Float_t w)               { fW0CPV = w ; }
  virtual void SetDigitsBranch(const char * title) { fDigitsBranchTitle = title  ;}
  virtual void SetRecPointsBranch(const char *title){fRecPointsBranchTitle = title; }
  virtual void SetUnfolding(Bool_t toUnfold = kTRUE ) {fToUnfold = toUnfold ;}  
  static Double_t ShowerShape(Double_t r) ; // Shape of EM shower used in unfolding; 
                                            //class member function (not object member function)
  static void UnfoldingChiSquare(Int_t & nPar, Double_t * Grad, Double_t & fret, Double_t * x, Int_t iflag)  ;
                                            // Chi^2 of the fit. Should be static to be passes to MINUIT
  virtual const char * Version() const { return "clu-v1" ; }  

protected:

  void           WriteRecPoints(Int_t event) ;
  virtual void   MakeClusters( ) ;            
  virtual Bool_t IsInTower (AliEMCALDigit * digit)const ;     // Tells if id digit is in Tower
  virtual Bool_t IsInPreShower (AliEMCALDigit * digit)const ;     // Tells if id digit is in PreShower

  
private:

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
  TString fHeaderFileName ;          // name of the file which contains gAlice, Tree headers etc.
  TString fDigitsBranchTitle ;       // name of the file, where digits branch is stored
  TString fRecPointsBranchTitle ;    // name of the file, where RecPoints branchs are stored

  Int_t   fNTowers ;                 // number of Towers in EMCAL

  Bool_t  fToUnfold ;                // To perform unfolding 

  Int_t   fNumberOfTowerClusters ;     // number of Tower clusters found 
  Int_t   fNumberOfPreShoClusters ;    // number of PreShower clusters found
 
  Float_t fADCchannelTower ;           // width of one ADC channel for Tower (GeV)
  Float_t fADCpedestalTower ;          // pedestal of ADC for Tower (GeV) 
  Float_t fADCchannelPreSho ;          // width of one ADC channel for Pre Shower (GeV)
  Float_t fADCpedestalPreSho ;         // pedestal of ADC for PreShower (GeV)

  Float_t fTowerClusteringThreshold ;  // minimum energy to include a EMC digit in a cluster
  Float_t fPreShoClusteringThreshold ; // minimum energy to include a CPV digit in a cluster
  Float_t fTowerLocMaxCut ;            // minimum energy difference to distinguish local maxima in a cluster
  Float_t fW0 ;                        // logarithmic weight for the cluster center of gravity calculation
  Float_t fPreShoLocMaxCut ;           //  minimum energy difference to distinguish local maxima in a CPV cluster
  Float_t fW0CPV ;                   // logarithmic weight for the CPV cluster center of gravity calculation
  Int_t fRecPointsInRun ;            //! Total number of recpoints in one run
  Float_t fTimeGate ;                // Maximum time difference between the digits in ont EMC cluster
    
  ClassDef(AliEMCALClusterizerv1,1)   // Clusterizer implementation version 1

};

#endif // AliEMCALCLUSTERIZERV1_H
