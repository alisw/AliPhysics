//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTTRDCLUSTERIZER_H
#define ALIHLTTRDCLUSTERIZER_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  HLT TRD cluster finder                                                //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include "AliTRDclusterizer.h"
#include "AliTRDReconstructor.h"
#include "AliHLTDataTypes.h"
#include "AliHLTTRDTrackletWordArray.h"

class AliHLTTRDClustersArray;
class AliHLTTRDClusterizer : public AliTRDclusterizer
{
 public:
  AliHLTTRDClusterizer(const AliTRDReconstructor *const rec = 0x0);
  AliHLTTRDClusterizer(const Text_t *const name, const Text_t *const title, const AliTRDReconstructor *const rec = 0x0);
  AliHLTTRDClusterizer(const AliHLTTRDClusterizer& c);
  virtual ~AliHLTTRDClusterizer() {};
  AliHLTTRDClusterizer& operator=(const AliHLTTRDClusterizer& c);

  void            Copy(TObject& c) const;
  void            SetMemBlock(AliHLTUInt8_t* ptr){
    if(fReconstructor->IsProcessingTracklets()){
      fTrMemBlock=ptr; fTrMemCurrPtr=ptr;
      fClMemBlock=ptr+GetTrMemBlockSize();  //if IsProcessingTracklets() is enabled we always reserve a data block of size GetTrMemBlockSize() for the tracklets
    }else{
      fClMemBlock=ptr;
    }
    fNoOfClusters=0;
    fAddedSize=0;
    fLastDet=-1;
    fClusters=NULL;
  }
  AliHLTUInt8_t*  GetClMemBlock(){return fClMemBlock;}
  AliHLTUInt8_t*  GetTrMemBlock(){return fTrMemBlock;}
  UInt_t          GetAddedClSize(){return fAddedSize;}
  UInt_t          GetAddedTrSize(){return (AliHLTUInt8_t*)fTrMemCurrPtr-(AliHLTUInt8_t*)fTrMemBlock;}
  UInt_t          GetTrMemBlockSize(){return 30*(sizeof(AliHLTTRDTrackletWordArray)+512*sizeof(UInt_t));}

 protected:
  void            AddClusterToArray(AliTRDcluster* cluster);
  void            AddTrackletsToArray();

  TClonesArray*   RecPoints(){return 0x0;}       //these are functions in the parents class and must not be used in hlt!
  TClonesArray*   TrackletsArray(){return 0x0;}  //if used accidentally it may give a compilation error because they are protected,
  void  SetClustersOwner(Bool_t /*own*/){}       //but it could be that the error appears only in  run time
  
  AliHLTUInt8_t*  fClMemBlock;
  AliHLTUInt8_t*  fTrMemBlock;
  AliHLTUInt8_t*  fTrMemCurrPtr;
  Int_t           fLastDet;
  AliHLTTRDClustersArray* fClusters;
  AliHLTUInt32_t  fAddedSize;

  ClassDef(AliHLTTRDClusterizer, 1)
};

#endif
