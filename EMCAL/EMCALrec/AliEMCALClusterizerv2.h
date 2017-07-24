#ifndef ALIEMCALCLUSTERIZERV2_H
#define ALIEMCALCLUSTERIZERV2_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//_________________________________________________________________________
/// \class AliEMCALClusterizerv2
/// \ingroup EMCALrec
/// \brief Clusterize neighbour cells, split if several maxima
///
///  Implementation version 2 of the clusterization algorithm                     
///  Performs clusterization (collects neighbouring active cells). 
///  It does not allow having more than one local maxima. 
///  Results are stored in TreeR.
///
/// \author Constantin Loizides (LBNL)
//_________________________________________________________________________

#include "AliEMCALClusterizerv1.h"
class AliEMCALRecPoint; 
class AliEMCALDigit;

class AliEMCALClusterizerv2 : public AliEMCALClusterizerv1
{
  
public:
  
  AliEMCALClusterizerv2() ;         
  AliEMCALClusterizerv2(AliEMCALGeometry* geometry);
  AliEMCALClusterizerv2(AliEMCALGeometry* geometry, AliEMCALCalibData* calib, 
                        AliEMCALCalibTime * calibt, AliCaloCalibPedestal *pedestal);

  virtual ~AliEMCALClusterizerv2();

  virtual       Int_t AreNeighbours(AliEMCALDigit* d1, AliEMCALDigit* d2, Bool_t& shared) const; 
  
  virtual const char *Version() const { return "clu-v2";}  

  void                SetDoEnGradCut(Bool_t b) { fDoEnGradCut = b; }

protected:
  
  virtual void        MakeClusters();            

  Bool_t              fDoEnGradCut; ///< cut on energy gradient

private:

  AliEMCALClusterizerv2              (const AliEMCALClusterizerv2 &); //copy ctor
  AliEMCALClusterizerv2 & operator = (const AliEMCALClusterizerv2 &);

  /// \cond CLASSIMP
  ClassDef(AliEMCALClusterizerv2,1) ;
  /// \endcond

};

#endif // AliEMCALCLUSTERIZERV2_H
