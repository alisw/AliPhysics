// $Id$
#ifndef ALIHLTTPCCLUSTERFLAGSDATA_H
#define ALIHLTTPCCLUSTERFLAGSDATA_H

#include "Rtypes.h"

/**
 * @struct AliHLTTPCClusterFlagsData
 * Primitive data of the flags for a tpc cluster in compressed form
 */

struct AliHLTTPCClusterFlagsData
{
  UShort_t fVersion; // version number of storage format:
                   // 0: no cluster flags stored
                   // 1: Plain flags stored. fNumberOfFlags defines how many bits are used per cluster
                   // 2: RLE encoded storage - TO BE IMPLEMENTED
  UShort_t fNumberOfFlags;      // number of flags (currently by default 2 for splitting in pad and time direction)
  UInt_t fNumberOfClusters;   // number of clusters
  //Access to the following members might be via int32, so make sure it is aligned correctly!
#if defined(__HP_aCC) || defined(__DECCXX) || defined(__SUNPRO_CC)
  unsigned char fData[1]; // flags
#else
  unsigned char fData[0]; // flags 
#endif
};
typedef struct AliHLTTPCClusterFlagsData AliHLTTPCClusterFlagsData;

#endif
