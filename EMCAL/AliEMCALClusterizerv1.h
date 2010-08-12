#ifndef ALIEMCALCLUSTERIZERV1_H
#define ALIEMCALCLUSTERIZERV1_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//_________________________________________________________________________
//  Implementation version 1 of the clusterization algorithm                     
//  Performs clusterization (collects neighbouring active cells) and 
//  unfolding of the clusters with several local maxima.  
//  results are stored in TreeR
//
//*-- Author: Yves Schutz (SUBATECH)
//--          Gustavo Conesa (LPSC-Grenoble), move common clusterizer functionalities to mother class
                        

// --- ROOT system ---

// --- Standard library ---

// --- AliRoot header files ---
#include "AliEMCALClusterizer.h"
class AliEMCALRecPoint ; 
class AliEMCALDigit ;

class AliEMCALClusterizerv1 : public AliEMCALClusterizer {
  
public:
  
  AliEMCALClusterizerv1() ;         
  AliEMCALClusterizerv1(AliEMCALGeometry* geometry);
  AliEMCALClusterizerv1(AliEMCALGeometry* geometry, AliEMCALCalibData * calib, AliCaloCalibPedestal * pedestal);
	
  virtual ~AliEMCALClusterizerv1()  ;

  virtual Int_t   AreNeighbours(AliEMCALDigit * d1, AliEMCALDigit * d2, Bool_t & shared)const ; 
                               // Checks if digits are in neighbour cells
  virtual void    Digits2Clusters(Option_t *option);                // Does the job

  static Double_t ShowerShape(Double_t x, Double_t y) ; // Shape of EM shower used in unfolding; 
                                            //class member function (not object member function)
  static void UnfoldingChiSquare(Int_t & nPar, Double_t * Grad, Double_t & fret, Double_t * x, Int_t iflag)  ;
                                            // Chi^2 of the fit. Should be static to be passes to MINUIT
  virtual const char * Version() const { return "clu-v1" ; }  

protected:

  virtual void   MakeClusters();            

private:
  AliEMCALClusterizerv1(const AliEMCALClusterizerv1 &); //copy ctor
  AliEMCALClusterizerv1 & operator = (const AliEMCALClusterizerv1 &);


  Bool_t  FindFit(AliEMCALRecPoint * emcRP, AliEMCALDigit ** MaxAt, const Float_t * maxAtEnergy, 
		  Int_t NPar, Float_t * FitParametres) const; //Used in UnfoldClusters, calls TMinuit

  virtual void   MakeUnfolding();
  void           UnfoldCluster(AliEMCALRecPoint * iniEmc, Int_t Nmax, 
                               AliEMCALDigit ** maxAt,
                               Float_t * maxAtEnergy ); //Unfolds cluster using TMinuit package

private:


  ClassDef(AliEMCALClusterizerv1,10)   // Clusterizer implementation version 1

};

#endif // AliEMCALCLUSTERIZERV1_H
