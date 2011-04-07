#ifndef ALIITSLAYERUPGRADE_H
#define ALIITSLAYERUPGRADE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
#include <TObject.h>
#include "AliITStrackerMI.h"
#include "AliITSRecPoint.h"
#include "AliITSRecoParam.h"
 class AliITSlayerUpgrade : public TObject {

  public:
    AliITSlayerUpgrade();
    AliITSlayerUpgrade(Double_t p, Double_t z);
    virtual ~AliITSlayerUpgrade();
    Int_t InsertCluster(AliITSRecPoint *c);
    void ResetClusters();
    const AliITSRecPoint *GetNextCluster(Int_t &ci);
    AliITSRecPoint *GetCluster(Int_t i) const { return fClusters[i]; }
    Int_t GetNumberOfClusters() const { return fN; }
  protected:
    AliITSlayerUpgrade(const AliITSlayerUpgrade&);
    AliITSlayerUpgrade &operator=(const AliITSlayerUpgrade &tr);
    Double_t fPhiOffset;        // offset of the first detector in Phi
    Double_t fZOffset;          // offset of the first detector in Z

    AliITSRecPoint *fClusters[AliITSRecoParam::kMaxClusterPerLayer]; // pointers to clusters
    Int_t fNsel;         // numbers of selected clusters 
    Int_t fIndex[AliITSRecoParam::kMaxClusterPerLayer]; // indexes of selected clusters
    Int_t fN;                       // number of clusters


   ClassDef(AliITSlayerUpgrade,1);   

  };
#endif
