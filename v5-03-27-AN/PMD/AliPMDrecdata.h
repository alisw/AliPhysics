#ifndef ALIPMDRECDATA_H
#define ALIPMDRECDATA_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
//-----------------------------------------------------//
//                                                     //
//                                                     //
//  Date   : May 27, 2009                              //
//                                                     //
//  Store reconstructed points, track no and trackpid  //
//  corresponding to the cluster for PMD               //
//                                                     //
//-----------------------------------------------------//

#include "Rtypes.h"
#include "TObject.h"
class TClonesArray;

class AliPMDrecdata : public TObject
{

 public:
  AliPMDrecdata();
  AliPMDrecdata(Int_t idet, Int_t ismn, Int_t trno, Int_t trpid, Float_t *clusdata);
  AliPMDrecdata(AliPMDrecdata *pmdrecdata);
  AliPMDrecdata (const AliPMDrecdata &pmdrecdata);  // copy constructor
  AliPMDrecdata &operator=(const AliPMDrecdata &pmdrecdata); // assignment op
  
  virtual ~AliPMDrecdata();

  Int_t   GetDetector() const;
  Int_t   GetSMNumber() const;
  Int_t   GetClusTrackNo() const;
  Int_t   GetClusTrackPid() const;
  Float_t GetClusX() const;
  Float_t GetClusY() const;
  Float_t GetClusADC() const;
  Float_t GetClusCells() const;
  Float_t GetClusSigmaX() const;
  Float_t GetClusSigmaY() const;
  
 protected:

  Int_t   fDet;          // Detector No (0:PRE, 1:CPV)
  Int_t   fSMN;          // Serial Module No.
  Int_t   fTrackNo;      // Track number of the cluster
  Int_t   fTrackPid;     // Cluster Track pid
  Float_t fClusData[6];  // Array containing cluster information
  /*
    fDet         : Detector Number,  fSMN         : Serial Module Number
    fClusData[0] : Cluster x      ,  fClusData[1] : Cluster y
    fClusData[2] : Cluster adc    ,  fClusData[3] : Cluster Cells
    fClusData[4] : Cluster SigmaX ,  fClusData[5] : Cluster SigmaY
  */
  
  ClassDef(AliPMDrecdata,0) // keep reconstructed points info
};

#endif
