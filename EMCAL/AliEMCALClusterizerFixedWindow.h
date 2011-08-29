#ifndef ALIEMCALCLUSTERIZERFIXEDWINDOW_H
#define ALIEMCALCLUSTERIZERFIXEDWINDOW_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliEMCALClusterizerFixedWindows.h   */

//_________________________________________________________________________
// This class derives from AliEMCALClustrerizer

#include "AliEMCALClusterizer.h"

class AliEMCALRecPoint; 
class AliEMCALDigit;
class AliEMCALFixedWindowClusterInfo;

class AliEMCALClusterizerFixedWindow : public AliEMCALClusterizer {
	
public:
	
	AliEMCALClusterizerFixedWindow() ;         
	AliEMCALClusterizerFixedWindow(AliEMCALGeometry* geometry);
	AliEMCALClusterizerFixedWindow(AliEMCALGeometry* geometry, AliEMCALCalibData * calib, AliCaloCalibPedestal * pedestal);
	
	virtual ~AliEMCALClusterizerFixedWindow();
	
	virtual void   Digits2Clusters(Option_t *option);                // Does the job
	
	virtual const char * Version() const { return "clu-FixedWindow" ; }  
	
	void SetnPhi (Int_t n) 
  {
    if (clusters_array)
      AliWarning("Clusterizer already initialized. Unable to change the parameters.");
    else
      nPhi = n;
  }
  
	void SetnEta (Int_t n) 
  {
    if (clusters_array)
      AliWarning("Clusterizer already initialized. Unable to change the parameters.");
    else
      nEta = n;
  }
	
	Int_t GetnPhi () {return nPhi;}
	Int_t GetnEta () {return nEta;}
  
  void SetshiftPhi (Int_t s) 
  {
    if (clusters_array)
      AliWarning("Clusterizer already initialized. Unable to change the parameters.");
    else
      shiftPhi = s;
  }
  
  void SetshiftEta (Int_t s) 
  {
    if (clusters_array)
      AliWarning("Clusterizer already initialized. Unable to change the parameters.");
    else
      shiftEta = s;
  }
  
  Int_t GetshiftPhi () {return shiftPhi;}
  Int_t GetshiftEta () {return shiftEta;}
  
  void SetTRUshift(Bool_t b) 
  {
    if (clusters_array)
      AliWarning("Clusterizer already initialized. Unable to change the parameters.");
    else
      fTRUshift = b;
  }
  
  Bool_t GetTRUshift() {return fTRUshift;}
  
  AliEMCALFixedWindowClusterInfo* GetClustersInfo() {return fClustersInfo;}
  void SetClustersInfo(AliEMCALFixedWindowClusterInfo *ClusInfo) {fClustersInfo = ClusInfo;}

protected:
	
	virtual void   MakeClusters();            
	
private:
	AliEMCALClusterizerFixedWindow(const AliEMCALClusterizerFixedWindow &); //copy ctor
	AliEMCALClusterizerFixedWindow & operator = (const AliEMCALClusterizerFixedWindow &);
	
  // nPhi x nEta clusterizer
	// Those parameter could be changed to get other types of fixed windows.
	Int_t                               nPhi;
	Int_t                               nEta; 
  Int_t                               shiftPhi;
  Int_t                               shiftEta;
  Bool_t                              fTRUshift;
  AliEMCALFixedWindowClusterInfo    *fClustersInfo;
  AliEMCALDigit                    ***clusters_array;
	
	ClassDef(AliEMCALClusterizerFixedWindow,4)   // Clusterizer implementation version 1
};

#endif // AliEMCALCLUSTERIZERFIXEDWINDOW_H
