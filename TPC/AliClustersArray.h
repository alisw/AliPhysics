#ifndef ALICLUSTERSARRAY_H
#define ALICLUSTERSARRAY_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  Manager class for AliClustersArray        //
////////////////////////////////////////////////
class AliDetectorParam;
class AliClustersFinder;

class AliClustersArray : public AliSegmentArray {
public:
  AliClustersArray();
  virtual   Bool_t Setup(AliDetectorParam *param);  
  const AliDetectorParam *  GetParam() {return fParam;} 
  AliClustersFinder *  GetFinder() {return fClFinder;}
  virtual Bool_t SetParam(AliDetectorParam * param);
  virtual Bool_t SetFinder(AliClustersFinder * finder);
  Bool_t  SetClusterType(Text_t *classname );
  TClass * GetClusterType() {return fClusterType;}
protected:  
  AliDetectorParam * fParam;      //pointer to detector parameters 
  AliClustersFinder * fClFinder; //!pointer to cluster finder object
  TClass *fClusterType; //!
  ClassDef(AliClustersArray,1) 
};
  
#endif
