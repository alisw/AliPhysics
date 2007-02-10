#ifndef ALIPMDCLUDATA_H
#define ALIPMDCLUDATA_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
//-----------------------------------------------------//
//                                                     //
//  Date   : February 05 2007                          //
//                                                     //
//  Store cluster informations and used inside         //
//  AliPMDClusteringV1 and AliPMDClusteringV2          //
//  to pass information from one method to another     //
//                                                     //
//-----------------------------------------------------//
// Author -  Ajay Dash
//
#include "Rtypes.h"
#include "TObject.h"
class TClonesArray;

class AliPMDcludata : public TObject
{
 public:
  AliPMDcludata();
  AliPMDcludata( Float_t *clusdata);
  AliPMDcludata(AliPMDcludata *pmdcludata);
  AliPMDcludata (const AliPMDcludata &pmdcludata);  // copy constructor
  AliPMDcludata &operator=(const AliPMDcludata &pmdcludata); // assignment op
  
  virtual ~AliPMDcludata();

  Float_t GetClusX() const;
  Float_t GetClusY() const;
  Float_t GetClusADC() const;
  Float_t GetClusCells() const;
  Float_t GetClusSigmaX() const;
  Float_t GetClusSigmaY() const;

 protected:


  Float_t fClusData[6];       // Array containing cluster information
  /*
    fClusData[0] : Cluster x         , fClusData[1] : Cluster y
    fClusData[2] : Cluster adc       , fClusData[3] : Cluster Cells
    fClusData[4] : Cluster SigmaX    , fClusData[5] : Cluster SigmaY
  */
  
  ClassDef(AliPMDcludata,1) // Keep Cluster information
};

#endif
