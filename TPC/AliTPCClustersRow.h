#ifndef ALITPCCLUSTERROW_H
#define ALITPCCLUSTERROW_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  Manager class for TPC   clusters                   //
////////////////////////////////////////////////


#include "AliClusters.h"

class TObject;


class AliTPCClustersRow : public AliClusters{
public:
  AliTPCClustersRow();
  virtual TObject *InsertCluster(const TObject *c);

public:
  
  ClassDef(AliTPCClustersRow,1) // Cluster manager 
};  
#endif //ALITPCCLUSTERROW_H







