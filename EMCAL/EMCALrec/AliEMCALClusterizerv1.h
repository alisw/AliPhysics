#ifndef ALIEMCALCLUSTERIZERV1_H
#define ALIEMCALCLUSTERIZERV1_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//_________________________________________________________________________
/// \class AliEMCALClusterizerv1
/// \ingroup EMCALrec
/// \brief Clusterize neighbour cells, no split, unfolding possible
///
///  Implementation version 1 of the clusterization algorithm                     
///  Performs clusterization (collects neighbouring active cells) and 
///  unfolding of the clusters with several local maxima.  
///  Results are stored in TreeR.
///
/// \author Yves Schutz (SUBATECH)
/// \author Gustavo Conesa (LPSC-Grenoble), move common clusterizer functionalities to mother class
//_________________________________________________________________________
                        
// --- AliRoot header files ---
#include "AliEMCALClusterizer.h"
class AliEMCALRecPoint ; 
class AliEMCALDigit ;

class AliEMCALClusterizerv1 : public AliEMCALClusterizer 
{
  
public:
  
  AliEMCALClusterizerv1() ;         
  AliEMCALClusterizerv1(AliEMCALGeometry* geometry);
  AliEMCALClusterizerv1(AliEMCALGeometry* geometry, AliEMCALCalibData * calib,
                        AliEMCALCalibTime * calibt, AliCaloCalibPedestal *pedestal);
	
  virtual ~AliEMCALClusterizerv1()  ;

  virtual Int_t   AreNeighbours(AliEMCALDigit * d1, AliEMCALDigit * d2, Bool_t & shared)const ; 

  virtual void    Digits2Clusters(Option_t *option);             

  virtual const char * Version() const { return "clu-v1" ; }  

protected:

  virtual void   MakeClusters();            

private:
  
  AliEMCALClusterizerv1              (const AliEMCALClusterizerv1 &); //copy ctor
  AliEMCALClusterizerv1 & operator = (const AliEMCALClusterizerv1 &);

  /// \cond CLASSIMP
  ClassDef(AliEMCALClusterizerv1,10) ;
  /// \endcond

};

#endif // AliEMCALCLUSTERIZERV1_H
