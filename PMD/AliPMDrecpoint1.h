#ifndef ALIPMDRECPOINT1_H
#define ALIPMDRECPOINT1_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
//-----------------------------------------------------//
//                                                     //
//                                                     //
//  Date   : August 05 2003                            //
//                                                     //
//  Store reconstructed points  for PMD                //
//                                                     //
//-----------------------------------------------------//

#include "Rtypes.h"
#include "TObject.h"
class TClonesArray;

class AliPMDrecpoint1 : public TObject
{

 public:
  AliPMDrecpoint1();
  AliPMDrecpoint1(Int_t idet, Int_t ismn, Float_t *clusdata);
  AliPMDrecpoint1(AliPMDrecpoint1 *pmdrecpoint) {*this = *pmdrecpoint;}
  AliPMDrecpoint1 (const AliPMDrecpoint1 &pmdrecpoint);  // copy constructor
  AliPMDrecpoint1 &operator=(const AliPMDrecpoint1 &pmdrecpoint); // assignment op
  
  virtual ~AliPMDrecpoint1();

  Int_t   GetDetector() const;
  Int_t   GetSMNumber() const;
  Float_t GetClusX() const;
  Float_t GetClusY() const;
  Float_t GetClusADC() const;
  Float_t GetClusCells() const;
  Float_t GetClusSigmaX() const;
  Float_t GetClusSigmaY() const;
  
 protected:

  Int_t   fDet;          // Detector No (0:PRE, 1:CPV)
  Int_t   fSMN;          // Serial Module No.
  Float_t fClusData[6];  // Array containing cluster information
  /*
    fDet         : Detector Number,  fSMN         : Serial Module Number
    fClusData[0] : Cluster x      ,  fClusData[1] : Cluster y
    fClusData[2] : Cluster adc    ,  fClusData[3] : Cluster Cells
    fClusData[4] : Cluster SigmaX ,  fClusData[5] : Cluster SigmaY
  */
  
  ClassDef(AliPMDrecpoint1,3) // keep reconstructed points info
};

#endif
