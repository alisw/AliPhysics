#ifndef ALICLUSTERSARRAY_H
#define ALICLUSTERSARRAY_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  Manager class for AliClustersArray        //
////////////////////////////////////////////////

#include "AliSegmentArray.h"

class AliDetectorParam;
class AliClustersFinder;

class AliClustersArray : public AliSegmentArray {
public:
  AliClustersArray();
  virtual   Bool_t Setup(const AliDetectorParam *param);  
  const AliDetectorParam *  GetParam() {return fParam;} 
  AliClustersFinder *  GetFinder() {return fClFinder;}
  virtual Bool_t SetParam(AliDetectorParam * param);
  virtual Bool_t SetFinder(AliClustersFinder * finder);
  Bool_t  SetClusterType(const char *classname );
  TClass * GetClusterType() {return fClusterType;}
protected:  
  AliDetectorParam * fParam;      //pointer to detector parameters 
  AliClustersFinder * fClFinder; //!pointer to cluster finder object
  TClass *fClusterType; //!
  ClassDef(AliClustersArray,1) // Cluster manager
private:
  AliClustersArray(const AliClustersArray& r); //dummy copy constructor
  AliClustersArray &operator=(const AliClustersArray& r);//dummy assignment operator
};
  
#endif
