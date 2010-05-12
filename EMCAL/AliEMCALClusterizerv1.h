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
class TH1F;
class AliEMCALRecPoint ; 
class AliEMCALDigit ;
class AliEMCALDigitizer ;
class AliEMCALGeometry ;
class AliEMCALCalibData ;
class AliCaloCalibPedestal ;

class AliEMCALClusterizerv1 : public AliEMCALClusterizer {
  
public:
  
  AliEMCALClusterizerv1() ;         
  AliEMCALClusterizerv1(AliEMCALGeometry* geometry);
  AliEMCALClusterizerv1(AliEMCALGeometry* geometry, AliEMCALCalibData * calib, AliCaloCalibPedestal * pedestal);
	
  virtual ~AliEMCALClusterizerv1()  ;

  virtual Int_t   AreNeighbours(AliEMCALDigit * d1, AliEMCALDigit * d2, Bool_t & shared)const ; 
                               // Checks if digits are in neighbour cells 

  virtual Float_t Calibrate(const Float_t amp, const Float_t time, const Int_t cellId) ;  // Tranforms Amp to energy 

  virtual void    GetNumberOfClustersFound(int numb )const { numb = fNumberOfECAClusters ;} 
  virtual Float_t GetECAClusteringThreshold()        const { return fECAClusteringThreshold;}  
  virtual Float_t GetECALocalMaxCut()                const { return fECALocMaxCut;} 
  virtual Float_t GetECALogWeight()                  const { return fECAW0;}
  virtual Float_t GetMinECut()                       const { return fMinECut;}

  virtual Float_t GetTimeCut()                       const { return fTimeCut ; }
  virtual Float_t GetTimeMin()                       const { return fTimeMin ; }
  virtual Float_t GetTimeMax()                       const { return fTimeMax ; }

  virtual void    Digits2Clusters(Option_t *option);                // Does the job

  virtual void    Print(Option_t * option)const ;

  virtual void SetECAClusteringThreshold(Float_t cluth)  { fECAClusteringThreshold = cluth ; }
  virtual void SetMinECut(Float_t mine)                  { fMinECut = mine; }
  virtual void SetECALocalMaxCut(Float_t cut)            { fECALocMaxCut = cut ; }
  virtual void SetECALogWeight(Float_t w)                { fECAW0 = w ; }
  virtual void SetTimeCut(Float_t t)					 { fTimeCut = t ;}
  virtual void SetTimeMin(Float_t t)					 { fTimeMin = t ;}
  virtual void SetTimeMax(Float_t t)					 { fTimeMax = t ;}

  virtual void SetUnfolding(Bool_t toUnfold = kTRUE )    {fToUnfold = toUnfold ;}  
  static Double_t ShowerShape(Double_t x, Double_t y) ; // Shape of EM shower used in unfolding; 
                                            //class member function (not object member function)
  static void UnfoldingChiSquare(Int_t & nPar, Double_t * Grad, Double_t & fret, Double_t * x, Int_t iflag)  ;
                                            // Chi^2 of the fit. Should be static to be passes to MINUIT
  virtual const char * Version() const { return "clu-v1" ; }  

  void   PrintRecoInfo();                        //*MENU*
  void   SetCalibrationParameters(AliEMCALCalibData * calib)   { fCalibData = calib ; }
  void   SetCaloCalibPedestal(AliCaloCalibPedestal  * caloped) { fCaloPed   = caloped ; }

protected:

  virtual void   MakeClusters();            

private:
  AliEMCALClusterizerv1(const AliEMCALClusterizerv1 &); //copy ctor
  AliEMCALClusterizerv1 & operator = (const AliEMCALClusterizerv1 &);

  void    GetCalibrationParameters(void) ;
  void    GetCaloCalibPedestal(void) ;

  Bool_t  FindFit(AliEMCALRecPoint * emcRP, AliEMCALDigit ** MaxAt, const Float_t * maxAtEnergy, 
		  Int_t NPar, Float_t * FitParametres) const; //Used in UnfoldClusters, calls TMinuit
  void Init() ;
  void InitParameters();

  virtual void   MakeUnfolding();
  void           UnfoldCluster(AliEMCALRecPoint * iniEmc, Int_t Nmax, 
			       AliEMCALDigit ** maxAt,
			       Float_t * maxAtEnergy ); //Unfolds cluster using TMinuit package
  void           PrintRecPoints(Option_t * option) ;

private:
  AliEMCALGeometry* fGeom;           //! pointer to geometry for utilities

  Bool_t  fDefaultInit;              //! Says if the task was created by defaut ctor (only parameters are initialized)
  Bool_t  fToUnfold ;                // To perform unfolding 
  Int_t   fNumberOfECAClusters ;     // number of clusters found in EC section

  //Calibration parameters... to be replaced by database 

  AliEMCALCalibData    * fCalibData ;   //! Calibration database if aval
  AliCaloCalibPedestal * fCaloPed   ;   //! Tower status map if aval

  Float_t fADCchannelECA ;          // width of one ADC channel for EC section (GeV)
  Float_t fADCpedestalECA ;         // pedestal of ADC for EC section (GeV) 
 
  Float_t fECAClusteringThreshold ;  // minimum energy to seed a EC digit in a cluster
  Float_t fECALocMaxCut ;            // minimum energy difference to distinguish local maxima in a cluster
  Float_t fECAW0 ;                   // logarithmic weight for the cluster center of gravity calculation
  Float_t fTimeCut ;                 // Maximum time difference between the digits inside EMC cluster
  Float_t fTimeMin ;                 // Minimum time of physical signal in a cell/digiy
  Float_t fTimeMax ;                 // Maximum time of physical signal in a cell/digit
  Float_t fMinECut;                  // Minimum energy for a digit to be a member of a cluster

  ClassDef(AliEMCALClusterizerv1,9)   // Clusterizer implementation version 1

};

#endif // AliEMCALCLUSTERIZERV1_H
