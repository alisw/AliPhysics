
#ifndef ALIHLTPHOSPHYSICSANALYZER_H
#define ALIHLTPHOSPHYSICSANALYZER_H

/* Copyright(c) 2006, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                          */

//Intended to be a base class for analysis

#include "Rtypes.h"
#include "AliHLTPHOSConstants.h"
using namespace PhosHLTConst;


class TObjArray;
class TH1F;
class AliHLTPHOSClusterDataStruct;

const Float_t kCRYSTAL_SIZE = 2.25;

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
  
  TObjArray* fClustersPtr;                         //! /**<Pointer to the clusters to be analyzed*/
  TH1F* fRootHistPtr;                              //! /**<Pointer to the histograms which is to be filled*/

 private:
  Float_t fRotParametersCos[5];                        /**<Parameters for calculating global position*/
  Float_t fRotParametersSin[5];                        /**<Parameters for calculating global position*/
  Float_t fPHOSRadius;                                 /**<Distance from the IP to the crystals*/

  ClassDef(AliHLTPHOSPhysicsAnalyzer,1);
};

#endif
