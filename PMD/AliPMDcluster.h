#ifndef ALIPMDCLUSTER_H
#define ALIPMDCLUSTER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
//-----------------------------------------------------//
//                                                     //
//  Date   : August 05 2003                            //
//                                                     //
//  Store cluster informations for PMD                 //
//                                                     //
//-----------------------------------------------------//
// Author - B.K. Nandi
//
#include "Rtypes.h"
#include "TObject.h"
class TClonesArray;

class AliPMDcluster : public TObject
{
 public:
  AliPMDcluster();
  AliPMDcluster(Int_t idet, Int_t ismn, Float_t *clusdata,
		Int_t *celldataX, Int_t *celldataY,
		Int_t *celltrack, Int_t *cellpid, Float_t *celladc);
  AliPMDcluster(AliPMDcluster *pmdcluster);
  AliPMDcluster (const AliPMDcluster &pmdcluster);  // copy constructor
  AliPMDcluster &operator=(const AliPMDcluster &pmdcluster); // assignment op
  
  virtual ~AliPMDcluster();

  Int_t   GetDetector() const;
  Int_t   GetSMN() const;
  Float_t GetClusX() const;
  Float_t GetClusY() const;
  Float_t GetClusADC() const;
  Float_t GetClusCells() const;
  Float_t GetClusSigmaX() const;
  Float_t GetClusSigmaY() const;
  Int_t   GetClusCellX(Int_t i) const;
  Int_t   GetClusCellY(Int_t i) const;
  Int_t   GetClusCellTrack(Int_t i) const;
  Int_t   GetClusCellPid(Int_t i) const;
  Float_t GetClusCellAdc(Int_t i) const;

 protected:

  Int_t   fDet;               // Detector No (0:PRE, 1:CPV)
  Int_t   fSMN;               // Serial Module No.
  Int_t   fClusCellDataX[19]; // Array containing cell X
  Int_t   fClusCellDataY[19]; // Array containing cell Y
  Int_t   fClusCellTrack[19]; // Array containing cell track no
  Int_t   fClusCellPid[19];   // Array containing cell pid
  Float_t fClusCellAdc[19];   // Array containing cell adc
  Float_t fClusData[6];       // Array containing cluster information

  /*
    fDet         : Det (0:PRE, 1:CPV), fSMN         : SerialModuleNo
    fClusData[0] : Cluster x         , fClusData[1] : Cluster y
    fClusData[2] : Cluster adc       , fClusData[3] : Cluster Cells
    fClusData[4] : Cluster SigmaX    , fClusData[5] : Cluster SigmaY
  */
  
  ClassDef(AliPMDcluster,6) // Keep Cluster information
};

#endif
