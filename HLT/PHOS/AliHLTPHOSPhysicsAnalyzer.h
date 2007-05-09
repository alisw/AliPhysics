
#ifndef ALIHLTPHOSPHYSICSANALYZER_H
#define ALIHLTPHOSPHYSICSANALYZER_H

/* Copyright(c) 2006, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                          */


#include <TObjArray.h>
#include <TH1F.h>
#include "TH2F.h"
#include "AliHLTPHOSClusterDataStruct.h"

const Float_t CRYSTAL_SIZE = 2.25;

class AliHLTPHOSPhysicsAnalyzer
{
 public:
  AliHLTPHOSPhysicsAnalyzer();
  virtual ~AliHLTPHOSPhysicsAnalyzer();

  AliHLTPHOSPhysicsAnalyzer(const AliHLTPHOSPhysicsAnalyzer & );
  AliHLTPHOSPhysicsAnalyzer & operator = (const AliHLTPHOSPhysicsAnalyzer &) {return *this;}

  void SetHistogram(TH1F* histPtr) {fRootHistPtr = histPtr;}

  void LocalPosition(AliHLTPHOSClusterDataStruct* clusterPtr, Float_t* locPositionPtr);
  void GlobalPosition(AliHLTPHOSClusterDataStruct* clusterPtr, Float_t* positionPtr);
  void GlobalPosition(Float_t* locPositionPtr , Float_t* positionPtr, Int_t module);

  virtual void WriteHistogram(Char_t* fileName = "histogram.root");
  virtual void Analyze(AliHLTPHOSClusterDataStruct* clustersPtr[10000], Int_t nClusters) = 0;

 protected:
  
  TObjArray* fClustersPtr;
  TH1F* fRootHistPtr;

 private:
  Float_t fRotParametersCos[5];
  Float_t fRotParametersSin[5];
  Float_t fPHOSRadius;

  ClassDef(AliHLTPHOSPhysicsAnalyzer,1);
};

#endif
