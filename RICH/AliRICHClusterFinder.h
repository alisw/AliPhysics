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
  
  AliRICH *Rich() {return fRICH;}                                             //Pointer to RICH  
  void     Exec();                                                              //Loop on events and chambers  
  void     FindRawClusters(Int_t iChamber);                                     //Find raw clusters  
  void     FormRawCluster(Int_t i, Int_t j, AliRICHcluster &cluster);           //form a raw cluster
  void     FindLocalMaxima(AliRICHcluster &cluster);                            //Find local maxima in a cluster
  void     ResolveCluster(AliRICHcluster  &cluster);                            //Try to resolve a cluster with maxima > 2
  //PH  void     CoG();                                                               //Evaluate the CoG as the best  
  void     WriteRawCluster(AliRICHcluster &cluster);                            //write in the list of the raw clusters  
protected:      
  AliRICH                *fRICH;                         //Pointer to RICH
  AliHitMap              *fHitMap;                       //Hit Map with digit positions  
  ClassDef(AliRICHClusterFinder,0) //Finds raw clusters, trasfers them to resolved clusters through declustering.
};
#endif
