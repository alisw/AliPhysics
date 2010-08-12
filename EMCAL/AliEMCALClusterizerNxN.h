#ifndef ALIEMCALCLUSTERIZERNXN_H
#define ALIEMCALCLUSTERIZERNXN_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliEMCALClusterizerNxN.h 41181 2010-05-12 13:58:06Z gconesab $ */

//_________________________________________________________________________
// This class derives from AliEMCALClustrerizer but also keeps the API of AliEMCALClusterizerv1
// Algorithm:
// 1. peek the most energetic cell
// 2. assign it as a center of the cluster and add cells surrounding it: 3x3, 5x5...
// 3. remove the cells contributing to the cluster
// 4. start from 1 for the remaining clusters
// 5. cluster splitting (not implemented yet) - use the shape analysis to resolve the energy sharing
// - for high energy clusters check the surrounding of the 3x3 clusters for extra energy 
// (merge 3x3 clusters and resolve the internal energy sharing - case for 2 clusters merged)

#include "AliEMCALClusterizer.h"
class AliEMCALRecPoint ; 
class AliEMCALDigit ;

class AliEMCALClusterizerNxN : public AliEMCALClusterizer {
  
public:
  
  AliEMCALClusterizerNxN() ;         
  AliEMCALClusterizerNxN(AliEMCALGeometry* geometry);
  AliEMCALClusterizerNxN(AliEMCALGeometry* geometry, AliEMCALCalibData * calib, AliCaloCalibPedestal * pedestal);
	
  virtual ~AliEMCALClusterizerNxN()  ;

  virtual Int_t   AreNeighbours(AliEMCALDigit * d1, AliEMCALDigit * d2, Bool_t & shared)const ; 
                               // Checks if digits are in neighbour cells 

   virtual void   Digits2Clusters(Option_t *option);                // Does the job


  static Double_t ShowerShape(Double_t x, Double_t y) ; // Shape of EM shower used in unfolding; 
                                            //class member function (not object member function)
  static void     SetInputCalibrated(Bool_t val);

  virtual const char * Version() const { return "clu-NxN" ; }  
 
  
protected:

  virtual void   MakeClusters();            

  static Bool_t fgkIsInputCalibrated; // to enable reclusterization from ESD cells

private:
  AliEMCALClusterizerNxN(const AliEMCALClusterizerNxN &); //copy ctor
  AliEMCALClusterizerNxN & operator = (const AliEMCALClusterizerNxN &);


  virtual void   MakeUnfolding();
  void           UnfoldCluster(AliEMCALRecPoint * iniEmc, Int_t Nmax, 
                               AliEMCALDigit ** maxAt,
                               Float_t * maxAtEnergy ); //Unfolds cluster using TMinuit package

private:

   ClassDef(AliEMCALClusterizerNxN,1)   // Clusterizer implementation version 1

};

#endif // AliEMCALCLUSTERIZERNXN_H
