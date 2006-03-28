#ifndef ALITPCFAST_H
#define ALITPCFAST_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  fast TPC cluster simulation               //
////////////////////////////////////////////////

#include <TObject.h>

class AliRunLoader;
class AliTPCClustersArray;
class AliTPCParam;


class AliTPCFast : public TObject {

public:
  void Hits2Clusters(AliRunLoader* runLoader) const;
  void Hits2ExactClusters(AliRunLoader* runLoader) const;
  void Hits2ExactClustersSector(AliRunLoader* runLoader,
				AliTPCClustersArray* clustersArray,
				Int_t isec) const;

 protected:
  AliTPCParam * fParam;   //!pointer to parameters - NOT OWNER
  ClassDef(AliTPCFast,0)  // fast TPC cluster simulation
};

#endif
