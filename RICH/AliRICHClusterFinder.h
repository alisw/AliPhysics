#ifndef AliRICHClusterFinder_h
#define AliRICHClusterFinder_h

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "TTask.h"

class AliRICH;
class AliHitMap;
class AliRICHcluster;
class AliRICHClusterFinder : public TTask
{
public:
    
  AliRICHClusterFinder(AliRICH *pRICH);
  virtual ~AliRICHClusterFinder()                                          {;}
  
  void   Exec();                                                              //Loop on events and chambers  
  void   FindRawClusters(Int_t iChamber);                                                   //Find initial clusters  
  void   FindLocalMaxima(AliRICHcluster *pCluster);                           //Find local maxima in a cluster
  void   ResolveCluster(AliRICHcluster *pCluster);                            //Try to resolve a cluster with maxima > 2
  void   FormRawCluster(Int_t i, Int_t j, AliRICHcluster *pCluster);          //form a raw cluster
  void   CoG();                                                               //Evaluate the CoG as the best  
  void   WriteRawCluster(AliRICHcluster *pRawCluster);                        //write in the list of the raw clusters  
  AliRICH *Rich() {return fRICH;}
  
  protected:
      
  AliRICH                *fRICH;  
  AliHitMap              *fHitMap;                       //Hit Map with digit positions

  ClassDef(AliRICHClusterFinder,0) //Class for clustering and reconstruction of space points    
};
#endif
