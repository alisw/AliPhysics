#ifndef ALIEMCALUNFOLDING_H
#define ALIEMCALUNFOLDING_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
     
//_________________________________________________________________________
//  Base class for the cluster unfolding algorithm 
//*-- Author: Adam Matyja (SUBATECH)

// --- ROOT system ---
#include "AliLog.h"
#include "TObject.h" 
//class TTree;

// --- Standard library ---

// --- AliRoot header files ---
class AliEMCALGeometry ;
//class AliEMCALCalibData ;
//class AliCaloCalibPedestal ;
class AliEMCALRecPoint ; 
class AliEMCALDigit ;


class AliEMCALUnfolding : public TObject {

public:

  AliEMCALUnfolding() ;        // default ctor
  virtual ~AliEMCALUnfolding() ; // dtorEM
  AliEMCALUnfolding(AliEMCALGeometry* geometry);// constructor
  AliEMCALUnfolding(AliEMCALGeometry* geometry,Float_t ECALocMaxCut,Double_t *SSPars,Double_t *Par5,Double_t *Par6);// constructor

  virtual void Init() ;
  virtual void SetInput(Int_t numberOfECAClusters,TObjArray *recPoints,TClonesArray *digitsArr);

  //setters and getters
  virtual void SetNumberOfECAClusters(Int_t n) { fNumberOfECAClusters = n; }
  virtual Int_t GetNumberOfECAClusters() const { return fNumberOfECAClusters; }
  virtual void SetRecPoints(TObjArray *rec) { fRecPoints = rec; }
  virtual TObjArray * GetRecPoints() const { return fRecPoints; }
  virtual void SetDigitsArr(TClonesArray *digit) { fDigitsArr = digit; }
  virtual TClonesArray * GetDigitsArr() const { return fDigitsArr; }
  virtual void SetECALocalMaxCut(Float_t cut) { fECALocMaxCut = cut ; }
  virtual Float_t GetECALocalMaxCut() const { return fECALocMaxCut; }
  virtual void SetThreshold(Float_t energy) { fThreshold = energy; }
  virtual Float_t GetThreshold() const { return fThreshold; }

  //unfolding main methods
  virtual void   MakeUnfolding();
  static Double_t ShowerShapeV2(Double_t x, Double_t y) ; // Shape of EM shower used in unfolding; 
                                              //class member function (not object member function)
  static void UnfoldingChiSquareV2(Int_t & nPar, Double_t * Grad, Double_t & fret, Double_t * x, Int_t iflag)  ;
                                            // Chi^2 of the fit. Should be static to be passes to MINUIT
  virtual void SetShowerShapeParams(Double_t *pars) ;
  virtual Double_t* GetShowerShapeParams() const { return fSSPars ; }
  virtual void SetPar5(Double_t *pars) ;
  virtual Double_t* GetPar5() const { return fPar5 ; }
  virtual void SetPar6(Double_t *pars) ;
  virtual Double_t* GetPar6() const { return fPar6 ; }

protected:
  Int_t   fNumberOfECAClusters ;     // number of clusters found in EC section
  Float_t fECALocMaxCut ;            // minimum energy difference to distinguish local maxima in a cluster
  Float_t fThreshold ; //minimum energy for cell to be joined to a cluster
  AliEMCALGeometry     * fGeom;       //! pointer to geometry for utilities
  TObjArray    *fRecPoints; // Array with EMCAL clusters
  TClonesArray *fDigitsArr; // Array with EMCAL digits

private:
  AliEMCALUnfolding(const AliEMCALUnfolding &); //copy ctor
  AliEMCALUnfolding & operator = (const AliEMCALUnfolding &);
  
  Bool_t         UnfoldClusterV2(AliEMCALRecPoint * iniEmc, Int_t Nmax, 
				 AliEMCALDigit ** maxAt,
				 Float_t * maxAtEnergy ); //Unfolds cluster using TMinuit package
  Bool_t         UnfoldClusterV2old(AliEMCALRecPoint * iniEmc, Int_t Nmax, 
				    AliEMCALDigit ** maxAt,
				    Float_t * maxAtEnergy ); //Unfolds cluster using TMinuit package
  Bool_t  FindFitV2(AliEMCALRecPoint * emcRP, AliEMCALDigit ** MaxAt, const Float_t * maxAtEnergy, 
		    Int_t NPar, Float_t * FitParametres) const; //Used in UnfoldClusters, calls TMinuit

  static Double_t fSSPars[8];//! Unfolding shower shape parameters
  // function:
  // f(r)=exp(-(p0*r)^p1 * (1/(p2+p3*(p0*r)^p1)+p4/(1+p6*(p0*r)^p5) ) )
  // p0,p1,p2,p3,p4 are fixed
  // params p5 and p6 are phi-dependent and set in ShowerShapeV2
  static Double_t fPar5[3];//! UF SSPar nr 5 = p0 + phi*p1 + phi^2 *p2
  static Double_t fPar6[3];//! UF SSPar nr 6 = p0 + phi*p1 + phi^2 *p2
  static void EvalPar5(Double_t phi);
  static void EvalPar6(Double_t phi);
  static void EvalParsPhiDependence(Int_t absId, AliEMCALGeometry *geom);

  ClassDef(AliEMCALUnfolding,2)  // Unfolding algorithm class 
} ;

#endif // AliEMCALUNFOLDING_H
