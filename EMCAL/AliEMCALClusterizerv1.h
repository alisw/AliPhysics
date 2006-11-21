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

class AliEMCALClusterizerv1 : public AliEMCALClusterizer {
  
public:
  
  AliEMCALClusterizerv1() ;         
  //cpy ctor required by coding convention
  AliEMCALClusterizerv1(const AliEMCALClusterizerv1& clus);

  AliEMCALClusterizerv1(const TString alirunFileNameFile, const TString eventFolderName = AliConfig::GetDefaultEventFolderName());
  virtual ~AliEMCALClusterizerv1()  ;
  virtual void Browse(TBrowser* b);

  AliEMCALClusterizerv1 & operator = (const AliEMCALClusterizerv1 &) {
    Fatal("operator =", "not implemented") ;
    return *this ;
  }

  virtual Int_t   AreNeighbours(AliEMCALDigit * d1, AliEMCALDigit * d2)const ; 
                               // Checks if digits are in neighbour cells 

  // Checks if digits are in a tower group; for pseudoclusters 
  virtual Int_t   AreInGroup(AliEMCALDigit * d1, AliEMCALDigit * d2)const ; 

  virtual Float_t Calibrate(Int_t amp, Int_t cellId) ;  // Tranforms Amp to energy 

  virtual void    GetNumberOfClustersFound(int numb )const{ numb = fNumberOfECAClusters ;} 
  virtual Float_t GetECAClusteringThreshold()const{ return fECAClusteringThreshold;}  
  virtual Float_t GetECALocalMaxCut()const       { return fECALocMaxCut;} 
  virtual Float_t GetECALogWeight()const         { return fECAW0;}
  virtual Float_t GetMinECut()const              { return fMinECut;}

  virtual Float_t GetTimeGate() const            { return fTimeGate ; }
  virtual const char *  GetRecPointsBranch() const{ return GetName() ;}
  virtual Int_t GetRecPointsInRun() const   {return fRecPointsInRun ;} 

  void    Exec(Option_t *option);                // Does the job

  virtual void Print(Option_t * option)const ;

  virtual void SetECAClusteringThreshold(Float_t cluth)  { fECAClusteringThreshold = cluth ; }
  virtual void SetMinECut(Float_t mine)                  { fMinECut = mine; }
  virtual void SetECALocalMaxCut(Float_t cut)            { fECALocMaxCut = cut ; }
  virtual void SetECALogWeight(Float_t w)                { fECAW0 = w ; }
  virtual void SetTimeGate(Float_t gate)                 { fTimeGate = gate ;}
  virtual void SetUnfolding(Bool_t toUnfold = kTRUE )    {fToUnfold = toUnfold ;}  
  static Double_t ShowerShape(Double_t r) ; // Shape of EM shower used in unfolding; 
                                            //class member function (not object member function)
  static void UnfoldingChiSquare(Int_t & nPar, Double_t * Grad, Double_t & fret, Double_t * x, Int_t iflag)  ;
                                            // Chi^2 of the fit. Should be static to be passes to MINUIT
  void Unload() ; 
  virtual const char * Version() const { return "clu-v1" ; }  
 
  TList* BookHists();
  void   SaveHists(const char *fn="reco.root");  //*MENU*
  void   PrintRecoInfo();                        //*MENU*
  void   DrawLambdasHists();                     //*MENU*
protected:

  void           WriteRecPoints() ;
  virtual void   MakeClusters(char* opt ) ;            
  virtual void   MakeClusters() { Fatal("MakeClusters","not implemented"); }
            
///////////////////// 
   TList  *fHists;   //!
   TH1F* fPointE;    //histogram of point energy
   TH1F* fPointL1;   //histogram of point L1
   TH1F* fPointL2;   //histogram of point L2
   TH1F* fPointDis;  //histogram of point dispersion
   TH1F* fPointMult; //histogram of point multiplicity
   TH1F* fDigitAmp;  //histogram of digit ADC Amplitude
   TH1F* fMaxE;      //histogram of maximum point energy
   TH1F* fMaxL1;     //histogram of largest (first) of eigenvalue of covariance matrix
   TH1F* fMaxL2;     //histogram of smalest (second) of eigenvalue of covariace matrix
   TH1F* fMaxDis;    //histogram of point dispersion
///////////////////////


private:

  const TString BranchName() const ; 
  void    GetCalibrationParameters(void) ;
  
  Bool_t  FindFit(AliEMCALRecPoint * emcRP, AliEMCALDigit ** MaxAt, Float_t * maxAtEnergy, 
		  Int_t NPar, Float_t * FitParametres) const; //Used in UnfoldClusters, calls TMinuit
  void Init() ;
  void InitParameters() ;

  virtual void   MakeUnfolding() const;
  void           UnfoldCluster(AliEMCALRecPoint * /*iniEmc*/, Int_t /*Nmax*/, 
			       AliEMCALDigit ** /*maxAt*/,
			       Float_t * /*maxAtEnergy*/ ) const; //Unfolds cluster using TMinuit package
  void           PrintRecPoints(Option_t * option) ;

private:
  AliEMCALGeometry* fGeom;           //! pointer to geometry for utilities

  Bool_t  fDefaultInit;              //! Says if the task was created by defaut ctor (only parameters are initialized)
  Bool_t  fToUnfold ;                // To perform unfolding 
  Int_t   fNumberOfECAClusters ;     // number of clusters found in EC section

  Int_t   fNTowerInGroup;            // number of towers to group for pseudoclusters

  //Calibration parameters... to be replaced by database 

  AliEMCALCalibData * fCalibData  ;   //! Calibration database if aval
  Float_t fADCchannelECA ;          // width of one ADC channel for EC section (GeV)
  Float_t fADCpedestalECA ;         // pedestal of ADC for EC section (GeV) 
 
  Float_t fECAClusteringThreshold ;  // minimum energy to seed a EC digit in a cluster
  Float_t fECALocMaxCut ;            // minimum energy difference to distinguish local maxima in a cluster
  Float_t fECAW0 ;                   // logarithmic weight for the cluster center of gravity calculation
  Int_t   fRecPointsInRun ;            //! Total number of recpoints in one run
  Float_t fTimeGate ;                // Maximum time difference between the digits in ont EMC cluster
  Float_t fMinECut;                  // Minimum energy for a digit to be a member of a cluster

  ClassDef(AliEMCALClusterizerv1,4)   // Clusterizer implementation version 1

};

#endif // AliEMCALCLUSTERIZERV1_H
