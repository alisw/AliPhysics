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
  AliPMDrecpoint1(Float_t *clusdata);
  AliPMDrecpoint1(AliPMDrecpoint1 *pmdrecpoint) {*this = *pmdrecpoint;}
  AliPMDrecpoint1 (const AliPMDrecpoint1 &pmdrecpoint);  // copy constructor
  AliPMDrecpoint1 &operator=(const AliPMDrecpoint1 &pmdrecpoint); // assignment op
  
  virtual ~AliPMDrecpoint1();

  Float_t GetDetector() const;
  Float_t GetSMNumber() const;
  Float_t GetClusX() const;
  Float_t GetClusY() const;
  Float_t GetClusADC() const;
  Float_t GetClusCells() const;
  Float_t GetClusRadius() const;
  
 protected:

  Float_t fClusData[7];  // Array containing cluster information
  /*
    fClusData[0] : Detector Number,  fClusData[1] : SuperModule Number
    fClusData[2] : Cluster x      ,  fClusData[3] : Cluster y
    fClusData[4] : Cluster adc    ,  fClusData[5] : Cluster Cells
    fClusData[6] : Cluster radius
  */
  
  ClassDef(AliPMDrecpoint1,2) // keep reconstructed points info
};

#endif
