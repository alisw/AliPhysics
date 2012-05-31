#ifndef ALIPMDCLUPID_H
#define ALIPMDCLUPID_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
//-----------------------------------------------------//
//                                                     //
//  Date   : March 22 2004                             //
//                                                     //
//  Store cluster informations for PMD                 //
//  after Discrimination                               //
//                                                     //
//-----------------------------------------------------//
// Author - B.K. Nandi
//
#include "Rtypes.h"
#include "TObject.h"
class TClonesArray;

class AliPMDclupid : public TObject
{
 public:
  AliPMDclupid();
  AliPMDclupid(Int_t idet, Int_t ismn, Int_t trno, Int_t trpid,
	       Int_t mstat,Float_t *clusdata);
  AliPMDclupid(AliPMDclupid *pmdclupid);
  AliPMDclupid (const AliPMDclupid &pmdclupid);  // copy constructor
  AliPMDclupid &operator=(const AliPMDclupid &pmdclupid); // assignment op
  
  virtual ~AliPMDclupid();

  Int_t   GetDetector() const;
  Int_t   GetSMN() const;
  Int_t   GetClusTrackNo() const;
  Int_t   GetClusTrackPid() const;
  Int_t   GetClusMatching() const;
  Float_t GetClusX() const;
  Float_t GetClusY() const;
  Float_t GetClusADC() const;
  Float_t GetClusCells() const;
  Float_t GetClusSigmaX() const;
  Float_t GetClusSigmaY() const;
  Float_t GetClusPID() const;

 protected:

  Int_t   fDet;          // Detector No (0:PRE, 1:CPV)
  Int_t   fSMN;          // Serial Module No.
  Int_t   fTrNo;         // Track number from simulation
  Int_t   fTrPid;        // Pid from simulation
  Int_t   fMstatus;      // matching status
  Float_t fClusData[7];  // Array containing clupid information

  /*
    fDet         : Det (0:PRE, 1:CPV), fSMN         : SerialModuleNo
    fClusData[0] : Cluster x         , fClusData[1] : Cluster y
    fClusData[2] : Cluster adc       , fClusData[3] : Cluster Cells
    fClusData[4] : Cluster sigmax    , fClusData[5] : Cluster sigmay
    fClusData[6] : Cluster pid
  */
  
  ClassDef(AliPMDclupid,3) // Keep Cluster information
};

#endif
