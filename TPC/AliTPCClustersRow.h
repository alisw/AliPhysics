#ifndef ALITPCCLUSTERROW_H
#define ALITPCCLUSTERROW_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  Manager class for TPC   clusters                   //
////////////////////////////////////////////////

#include "AliDetector.h"
#include "AliSegmentArray.h"
#include "AliClusters.h"


#include <TClonesArray.h>


class TClonesArray;
class TObjArray;


class AliTPCClustersRow : public AliClusters{
public:
  AliTPCClustersRow();

public:
  
  ClassDef(AliTPCClustersRow,1) 
};  
#endif //ALITPCCLUSTERROW_H
