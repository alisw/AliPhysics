#ifndef ALIITSTRACKERU_H
#define ALIITSTRACKERU_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//-------------------------------------------------------------------------
//                ITS upgrade tracker base class
//-------------------------------------------------------------------------

#include "AliTracker.h"
#include "AliESDEvent.h"

class TTree;

//-------------------------------------------------------------------------
class AliITStrackerU : public AliTracker {

  public:
 
  AliITStrackerU();
  virtual ~AliITStrackerU();

  virtual Int_t Clusters2Tracks(AliESDEvent *event);
  virtual Int_t PropagateBack(AliESDEvent *event);
  virtual Int_t RefitInward(AliESDEvent *event);
  virtual Int_t LoadClusters(TTree *);
  virtual void UnloadClusters();
  virtual AliCluster *GetCluster(Int_t index) const;
  
  
  private:
      
  AliITStrackerU(const AliITStrackerU&);
  AliITStrackerU &operator=(const AliITStrackerU &tr);
  
  ClassDef(AliITStrackerU,1)   //ITS upgrade tracker
    
};
#endif
