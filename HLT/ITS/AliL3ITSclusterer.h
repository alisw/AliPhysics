#ifndef ALIL3ITSCLUSTERER_H
#define ALIL3ITSCLUSTERER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//-------------------------------------------------------------------------
//                          High Level Trigger ITS clusterer
//       This class derives completely from the off-line AliITSclustererV2.
//       The only difference is in the interface of calling it and stoting
//       the clusters's tree.
//      
//           Origin: Cvetan Cheshkov, CERN, Cvetan.Cheshkov@cern.ch 
//-------------------------------------------------------------------------

#include "AliITSclustererV2.h"

class TTree;
class AliITSgeom;
class AliRawReader;

//-------------------------------------------------------------------------
class AliL3ITSclusterer : public AliITSclustererV2 {
public:
  AliL3ITSclusterer():AliITSclustererV2(){fNModule = 0;}
  AliL3ITSclusterer(const AliITSgeom *geom);

  void Digits2Clusters(AliRawReader* rawReader,TTree *cTree);

private:
  Int_t fNModule;             // total number of modules
 
  ClassDef(AliL3ITSclusterer,1)   //HLT ITS clusterer
};

#endif
