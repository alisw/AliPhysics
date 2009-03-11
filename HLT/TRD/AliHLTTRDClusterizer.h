#ifndef ALIHLTTRDCLUSTERIZER_H
#define ALIHLTTRDCLUSTERIZER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  HLT TRD cluster finder                                                //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include "AliTRDclusterizer.h"
#include "AliHLTTRDCluster.h"
#include "AliHLTDataTypes.h"

class AliTRDReconstructor;
class TClonesArray;

class AliHLTTRDClusterizer : public AliTRDclusterizer
{
 public:
  AliHLTTRDClusterizer(const AliTRDReconstructor *const rec = 0x0);
  AliHLTTRDClusterizer(const Text_t *const name, const Text_t *const title, const AliTRDReconstructor *const rec = 0x0);
  AliHLTTRDClusterizer(const AliHLTTRDClusterizer& c);
  AliHLTTRDClusterizer& operator=(const AliHLTTRDClusterizer& c);

  void Copy(TObject& c) const;

  void            SetMemBlock(AliHLTUInt8_t* ptr){fMemBlock=ptr;fNoOfClusters=0;}
  AliHLTUInt8_t*  GetMemBlock(){return fMemBlock;}
  UInt_t          GetAddedSize(){return fNoOfClusters*sizeof(AliHLTTRDCluster);}

 protected:
  TClonesArray*   RecPoints(){return 0x0;}  //these are functions in the parents class. must not be used in hlt!
  void SetClustersOwner(Bool_t own){        //if used accidentally it may give an compilation error because are protected,
      if(own){ /*get rid of warning*/}      //but the error can also appear in run time
  } 
                                            
  
  void            AddClusterToArray(AliTRDcluster *cluster);
  
  AliHLTUInt8_t*  fMemBlock;

  ClassDef(AliHLTTRDClusterizer, 0)
};

#endif
