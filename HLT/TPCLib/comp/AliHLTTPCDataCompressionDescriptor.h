#ifndef ALIHLTTPCDATACOMPRESSIONDESCRIPTOR_H
#define ALIHLTTPCDATACOMPRESSIONDESCRIPTOR_H

#include "Rtypes.h"

/**
 * @struct AliHLTTPCDataCompressionDescriptor
 * The class describes properties of the compressed data
 * @ingroup alihlt_tpc_datastructs
 */
class AliHLTTPCDataCompressionDescriptor
{
 public:

  AliHLTTPCDataCompressionDescriptor(): fVersion(0), fMergedClustersFlag(1){}
  ~AliHLTTPCDataCompressionDescriptor(){}
  AliHLTTPCDataCompressionDescriptor(const AliHLTTPCDataCompressionDescriptor& other)
    : fVersion(other.fVersion)
    , fMergedClustersFlag(other.fMergedClustersFlag)
  {}

  AliHLTTPCDataCompressionDescriptor& operator=(const AliHLTTPCDataCompressionDescriptor& other){ 
    if( &other == this ) return *this;
    fVersion = other.fVersion;
    fMergedClustersFlag = other.fMergedClustersFlag;
    return *this;
  }

  Bool_t CheckSize( UInt_t size ) const {
    if( size<sizeof(UInt_t) ) return 0;
    if( fVersion==0 ) return ( size==sizeof(AliHLTTPCDataCompressionDescriptor));
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
