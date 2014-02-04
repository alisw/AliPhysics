#ifndef ALIHLTTPCRAWCLUSTERSDESCRIPTOR_H
#define ALIHLTTPCRAWCLUSTERSDESCRIPTOR_H

#include "Rtypes.h"

/**
 * @struct AliHLTTPCRawClustersDescriptor
 * The class describes properties of the raw clusters
 * @ingroup alihlt_tpc_datastructs
 */
class AliHLTTPCRawClustersDescriptor
{
 public:

  AliHLTTPCRawClustersDescriptor(): fVersion(0), fMergedClustersFlag(1){}
  ~AliHLTTPCRawClustersDescriptor(){}
  AliHLTTPCRawClustersDescriptor(const AliHLTTPCRawClustersDescriptor& other)
    : fVersion(other.fVersion)
    , fMergedClustersFlag(other.fMergedClustersFlag)
  {}

  AliHLTTPCRawClustersDescriptor& operator=(const AliHLTTPCRawClustersDescriptor& other){ 
    if( &other == this ) return *this;
    fVersion = other.fVersion;
    fMergedClustersFlag = other.fMergedClustersFlag;
    return *this;
  }

  Bool_t CheckSize( UInt_t size ) const {
    if( size<sizeof(UInt_t) ) return 0;
    if( fVersion==0 ) return ( size==sizeof(AliHLTTPCRawClustersDescriptor));
    return 0;
  }
  
  UInt_t GetVersion() const { return fVersion; }
  Int_t GetMergedClustersFlag() const { return fMergedClustersFlag; }
  
  void SetMergedClustersFlag( Int_t flag ){ fMergedClustersFlag=flag; }

 private:

  UInt_t fVersion; // version number
  Int_t fMergedClustersFlag; // flag tells if the clusters were merged at the branch borders
};

#endif
