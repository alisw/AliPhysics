#ifndef ALIEMCALCLUSTERIZERNXN_H
#define ALIEMCALCLUSTERIZERNXN_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//_________________________________________________________________________
/// \class AliEMCALClusterizerNxN
/// \ingroup EMCALrec
/// \brief Create clusters of maximum size NxM
///
/// This class derives from AliEMCALClusterizer but also keeps the API of AliEMCALClusterizerv1
/// Algorithm:
/// 1. pick the most energetic cell
/// 2. assign it as a center of the cluster and add cells surrounding it: 3x3, 5x5...
/// 3. remove the cells contributing to the cluster
/// 4. start from 1 for the remaining clusters
/// 5. cluster splitting (not implemented yet) - use the shape analysis to resolve the energy sharing
/// - for high energy clusters check the surrounding of the 3x3 clusters for extra energy 
/// (merge 3x3 clusters and resolve the internal energy sharing - case for 2 clusters merged)
/// Use Case:
///  root [0] AliEMCALClusterizerNxN * cl = new AliEMCALClusterizerNxN("galice.root")  
///  Warning in <TDatabasePDG::TDatabasePDG>: object already instantiated
///               //reads gAlice from header file "..."                      
///  root [1] cl->ExecuteTask()  
///               //finds RecPoints in all events stored in galice.root
///  root [2] cl->SetDigitsBranch("digits2") 
///               //sets another title for Digitis (input) branch
///  root [3] cl->SetRecPointsBranch("recp2")  
///               //sets another title four output branches
///  root [4] cl->SetTowerLocalMaxCut(0.03)  
///               //set clusterization parameters
///  root [5] cl->ExecuteTask("deb all time")  
///               //once more finds RecPoints options are 
///               // deb - print number of found rec points
///               // deb all - print number of found RecPoints and some their characteristics 
///               // time - print benchmarking results
///
/// \author Mateusz Ploskon, LBNL
///
//_________________________________________________________________________

#include "AliEMCALClusterizer.h"
class AliEMCALRecPoint ; 
class AliEMCALDigit ;

class AliEMCALClusterizerNxN : public AliEMCALClusterizer
{
public:
  
  AliEMCALClusterizerNxN() ;         
  AliEMCALClusterizerNxN(AliEMCALGeometry* geometry);
  AliEMCALClusterizerNxN(AliEMCALGeometry* geometry, AliEMCALCalibData * calib,
                         AliEMCALCalibTime * calibt, AliCaloCalibPedestal *pedestal);
	
  virtual ~AliEMCALClusterizerNxN()  ;

  virtual Int_t AreNeighbours(AliEMCALDigit * d1, AliEMCALDigit * d2, Bool_t & shared) const ; 
                               // Checks if digits are in neighbour cells 

  virtual void  Digits2Clusters(Option_t *option);                // Does the job

  virtual const char * Version() const { return "clu-NxN" ; }  
  
  void   SetNRowDiff(Int_t nd)   { fNRowDiff = nd     ; }
  void   SetNColDiff(Int_t nd)   { fNColDiff = nd     ; }
  
  Int_t  GetNRowDiff()     const { return fNRowDiff   ; } 
  Int_t  GetNColDiff()     const { return fNColDiff   ; } 
  
  void   SetEnergyGrad(Bool_t b) { fEnergyGrad = b    ; }
  Bool_t GetEnergyGrad()   const { return fEnergyGrad ; }

protected:

  virtual void   MakeClusters();            

private:
  
  AliEMCALClusterizerNxN              (const AliEMCALClusterizerNxN &); //copy ctor
  AliEMCALClusterizerNxN & operator = (const AliEMCALClusterizerNxN &);

  Int_t  fNRowDiff;   ///< How many neighbors to consider along row (phi)
  Int_t  fNColDiff;   ///< How many neighbors to consider along col (eta)
  Bool_t fEnergyGrad; ///< If true only cluster if neighboring cell has less energy

  /// \cond CLASSIMP
  ClassDef(AliEMCALClusterizerNxN,4) ;
  /// \endcond

};

#endif // AliEMCALCLUSTERIZERNXN_H
