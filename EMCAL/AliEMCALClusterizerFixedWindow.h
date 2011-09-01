#ifndef ALIEMCALCLUSTERIZERFIXEDWINDOW_H
#define ALIEMCALCLUSTERIZERFIXEDWINDOW_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliEMCALClusterizerFixedWindow.h   */

//_________________________________________________________________________
// This class derives from AliEMCALClustrerizer

#include "AliEMCALClusterizer.h"

class AliEMCALRecPoint; 
class AliEMCALDigit;

class AliEMCALClusterizerFixedWindow : public AliEMCALClusterizer {
public:
	AliEMCALClusterizerFixedWindow() ;         
	AliEMCALClusterizerFixedWindow(AliEMCALGeometry* geometry);
	AliEMCALClusterizerFixedWindow(AliEMCALGeometry* geometry, AliEMCALCalibData * calib, AliCaloCalibPedestal * pedestal);
	virtual ~AliEMCALClusterizerFixedWindow();
	
public:
	virtual void            Digits2Clusters(Option_t *option);
	virtual const char     *Version() const { return "clu-FixedWindow"; }  

  Int_t                             GetNphi ()                                          const { return fNphi;             }
	Int_t                             GetNeta ()                                          const { return fNeta;             }
  Int_t                             GetShiftPhi ()                                      const { return fShiftPhi;         }
  Int_t                             GetShiftEta ()                                      const { return fShiftEta;         }
  Bool_t                            GetTRUshift()                                       const { return fTRUshift;         }
	void                              SetNphi (Int_t n);
	void                              SetNeta (Int_t n);
  void                              SetShiftPhi (Int_t s);
  void                              SetShiftEta (Int_t s);
  void                              SetTRUshift(Bool_t b);
  
protected:
	virtual void MakeClusters(); 
  
	Int_t                               fNphi;                // Fixed window number of cells in phi direction
	Int_t                               fNeta;                // Fixed window number of cells in eta direction
  Int_t                               fShiftPhi;            // Shifting number of cells in phi direction
  Int_t                               fShiftEta;            // Shifting number of cells in eta direction
  Bool_t                              fTRUshift;            // Allows shifting inside a TRU (true) of through the whole calorimeter (false)
  AliEMCALDigit                    ***fClustersArray;       //!Temporary array that contains clusters
	
private:
	AliEMCALClusterizerFixedWindow(const AliEMCALClusterizerFixedWindow &);                 // not implemented
	AliEMCALClusterizerFixedWindow & operator = (const AliEMCALClusterizerFixedWindow &);   // not implemented
	

	ClassDef(AliEMCALClusterizerFixedWindow,4)   // Clusterizer implementation version 1
};

#endif // AliEMCALCLUSTERIZERFIXEDWINDOW_H
