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
class AliHLTITSclusterer : public AliITSclustererV2 {
public:
  AliHLTITSclusterer():AliITSclustererV2(){fNModule = 0;}
  AliHLTITSclusterer(const AliITSgeom *geom);

  void Digits2Clusters(AliRawReader* rawReader,TTree *cTree);

private:
  Int_t fNModule;             // total number of modules
 
  ClassDef(AliHLTITSclusterer,1)   //HLT ITS clusterer
};

typedef AliHLTITSclusterer AliL3ITSclusterer; // for backward compatibility

#endif
