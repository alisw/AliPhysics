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
class AliPHOSEmcRecPoint ; 
class AliPHOSDigit ;
class AliPHOSDigitizer ;
class AliPHOSGeometry ;


class AliPHOSClusterizerv1 : public AliPHOSClusterizer {
  
public:
  
  AliPHOSClusterizerv1() ;             // ctor            
  AliPHOSClusterizerv1(const char * headerFile,const char *digitsBrancheTitle=0);
  virtual ~AliPHOSClusterizerv1(){}    // dtor
  
  Int_t           AreNeighbours(AliPHOSDigit * d1, AliPHOSDigit * d2)const ; 
                               // Checks if digits are in neighbour cells 
  virtual void    GetNumberOfClustersFound(int * numb ){  numb[0] = fNumberOfEmcClusters ; 
                                                          numb[1] = fNumberOfCpvClusters ; }

  virtual Float_t GetEmcClusteringThreshold()const{ return fEmcClusteringThreshold;}
  virtual Float_t GetEmcLocalMaxCut()const        { return fEmcLocMaxCut;} 
  virtual Float_t GetEmcLogWeight()const          { return fW0;}  
  virtual Float_t GetCpvClusteringThreshold()const{ return fCpvClusteringThreshold;  } 
  virtual Float_t GetCpvLocalMaxCut()const        { return fCpvLocMaxCut;} 
  virtual Float_t GetCpvLogWeight()const          { return fW0CPV;}  
  virtual Float_t GetPpsdClusteringThreshold()const{ return fPpsdClusteringThreshold;  } 
  virtual char *  GetRecPointsBranch() const      { return (char*) fRecPointsBranchTitle.Data() ;}
  virtual char *  GetDigitsBranch() const         { return (char*) fDigitsBranchTitle.Data() ;}

  void    Exec(Option_t *option);                // Does the job

  virtual void Print(Option_t * option)const ;

  virtual void SetEmcClusteringThreshold(Float_t cluth)  { fEmcClusteringThreshold = cluth ; }
  virtual void SetEmcLocalMaxCut(Float_t cut)            { fEmcLocMaxCut = cut ; }
  virtual void SetEmcLogWeight(Float_t w)                { fW0 = w ; }
  virtual void SetCpvClusteringThreshold(Float_t cluth)  { fCpvClusteringThreshold = cluth ; }
  virtual void SetCpvLocalMaxCut(Float_t cut)            { fCpvLocMaxCut = cut ; }
  virtual void SetCpvLogWeight(Float_t w)                { fW0CPV = w ; }
  virtual void SetPpsdClusteringThreshold(Float_t cluth) { fPpsdClusteringThreshold = cluth ; }

  virtual void SetDigitsBranch(const char * title) ; 
  virtual void SetRecPointsBranch(const char *title) ;

  virtual void SetUnfolding(Bool_t toUnfold = kTRUE ) {fToUnfold = toUnfold ;}  

  static void UnfoldingChiSquare(Int_t & nPar, Double_t * Grad, Double_t & fret, Double_t * x, Int_t iflag)  ;
                                            // Chi^2 of the fit. Should be static to be passes to MINUIT
  static Double_t ShowerShape(Double_t r) ; // Shape of EM shower used in unfolding; 
                                            //class member function (not object member function)

private:
  virtual Float_t Calibrate(Int_t amp)const {  return (amp-fPedestal)/fSlope ;}  // Tranforms Amp to energy 
  Bool_t  FindFit(AliPHOSEmcRecPoint * emcRP, int * MaxAt, Float_t * maxAtEnergy, 
		  Int_t NPar, Float_t * FitParametres) ; //Used in UnfoldClusters, calls TMinuit
  void Init() ;

  virtual Bool_t IsInEmc (AliPHOSDigit * digit)const ;     // Tells if id digit is in EMC
  virtual Bool_t IsInPpsd(AliPHOSDigit * digit)const ;     // Tells if id digit is in PPSD
  virtual Bool_t IsInCpv (AliPHOSDigit * digit)const ;     // Tells if id digit is in CPV

  virtual void   MakeClusters( ) ;            
  virtual void   MakeUnfolding() ;
  Bool_t         ReadDigits() ;
  void           UnfoldCluster(AliPHOSEmcRecPoint * iniEmc,Int_t Nmax, 
		       int * maxAt,Float_t * maxAtEnergy ) ; //Unfolds cluster using TMinuit package
  void           WriteRecPoints() ;
  void           PrintRecPoints(Option_t * option) ;

private:

  TString fHeaderFileName ;          // name of the file which contains gAlice, Tree headers etc.
  TString fDigitsBranchTitle ;       // name of the file, where digits branch is stored
  TString fRecPointsBranchTitle ;    // name of the file, where RecPoints branchs are stored

  Int_t   fEvent ;                   // Number of event currently processed 
  Bool_t  fToUnfold ;                // To perform unfolding 

  Bool_t  fIsInitialized ;

  AliPHOSGeometry * fGeom ;          // !pointer to PHOS geometry

  AliPHOSDigitizer * fDigitizer ;    // !digitizer which produced Digits we treat

  Int_t   fNumberOfEmcClusters ;     // number of EMC clusters found 
  Int_t   fNumberOfCpvClusters ;     // number of CPV+PPSD clusters found
  TClonesArray * fDigits ;           // ! Initial list of digits
  TObjArray    * fEmcRecPoints ;     // ! Final list of EMC Rec Points
  TObjArray    * fCpvRecPoints ;     // ! Final list of CPV/PPSD recPoints

  Float_t fPedestal ;                // Calibration parameters 
  Float_t fSlope ;                   // read from Digitizer

  Float_t fEmcClusteringThreshold ;  // minimum energy to include a EMC digit in a cluster
  Float_t fPpsdClusteringThreshold ; // minimum energy to include a PPSD digit in a cluster
  Float_t fCpvClusteringThreshold ;  // minimum energy to include a CPV digit in a cluster
  Float_t fEmcLocMaxCut ;            // minimum energy difference to distinguish local maxima in a cluster
  Float_t fW0 ;                      // logarithmic weight for the cluster center of gravity calculation
  Float_t fCpvLocMaxCut ;            // minimum energy difference to distinguish local maxima in a CPV cluster
  Float_t fW0CPV ;                   // logarithmic weight for the CPV cluster center of gravity calculation
    
  ClassDef(AliPHOSClusterizerv1,1)   // Clusterizer implementation version 1

};

#endif // AliPHOSCLUSTERIZERV1_H
