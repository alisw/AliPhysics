#ifndef ALIEMCALCLUSTERIZERFIXEDWINDOW_H
#define ALIEMCALCLUSTERIZERFIXEDWINDOW_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//_________________________________________________________________________
// This class derives from AliEMCALClustrerizer

#include "AliEMCALClusterizer.h"

class AliEMCALRecPoint; 
class AliEMCALDigit;

class AliEMCALClusterizerFixedWindow : public AliEMCALClusterizer {
 public:
  AliEMCALClusterizerFixedWindow() ;         
  AliEMCALClusterizerFixedWindow(AliEMCALGeometry* geometry);
  AliEMCALClusterizerFixedWindow(AliEMCALGeometry* geometry, AliEMCALCalibData * calib, 
                                 AliEMCALCalibTime * calibt, AliCaloCalibPedestal *pedestal);
  virtual ~AliEMCALClusterizerFixedWindow();
	
  virtual void            Digits2Clusters(Option_t *option);
  virtual const char     *Version() const { return "clu-FixedWindow"; }  

  Int_t                   GetNphi ()                                          const { return fNphi;             }
  Int_t                   GetNeta ()                                          const { return fNeta;             }
  Int_t                   GetShiftPhi ()                                      const { return fShiftPhi;         }
  Int_t                   GetShiftEta ()                                      const { return fShiftEta;         }
  Bool_t                  GetTRUshift()                                       const { return fTRUshift;         }
  void                    SetNphi (Int_t n);
  void                    SetNeta (Int_t n);
  void                    SetShiftPhi (Int_t s);
  void                    SetShiftEta (Int_t s);
  void                    SetTRUshift(Bool_t b);
  
protected:
  void MakeClusters(); 
  void ExecOnce(); 
  
  Int_t                   fNphi;                // Fixed window number of cells in phi direction
  Int_t                   fNeta;                // Fixed window number of cells in eta direction
  Int_t                   fShiftPhi;            // Shifting number of cells in phi direction
  Int_t                   fShiftEta;            // Shifting number of cells in eta direction
  Bool_t                  fTRUshift;            // Allows shifting inside a TRU (true) of through the whole calorimeter (false)
  Int_t                   fNEtaDigitsSupMod;    //!Number of digits per SM in eta 
  Int_t                   fNPhiDigitsSupMod;    //!Number of digits per SM in phi
  Int_t                   fNTRUPhi;             //!Number of TRUs in phi
  Int_t                   fNTRUEta;             //!Number of TRUs in eta
  Int_t                   fNEtaDigits;          //!Total number of digits in eta 
  Int_t                   fNPhiDigits;          //!Total number of digits in phi
  Int_t                   fMaxShiftPhi;         //!Max shift index in phi
  Int_t                   fMaxShiftEta;         //!Max shift index in eta
  Int_t                   fNDigitsCluster;      //!Digits per cluster
  Int_t                   fNClusEtaNoShift;     //!Max number of clusters in eta
  Int_t                   fNClusPhiNoShift;     //!Max number of clusters in phi
  Int_t                   fNClusters;           //!fNClusEtaNoShift x fNClusPhiNoShift
  Int_t                   fNTotalClus;          //!Maximum total number of clusters
  AliEMCALDigit        ***fClustersArray;       //!Temporary array that contains clusters
  Int_t                   fInitialized;         //!Initialized clusterizer
	
private:
  AliEMCALClusterizerFixedWindow(const AliEMCALClusterizerFixedWindow &);                 // not implemented
  AliEMCALClusterizerFixedWindow & operator = (const AliEMCALClusterizerFixedWindow &);   // not implemented

  ClassDef(AliEMCALClusterizerFixedWindow,4)   // Clusterizer implementation fixed windows
};
#endif // AliEMCALCLUSTERIZERFIXEDWINDOW_H
