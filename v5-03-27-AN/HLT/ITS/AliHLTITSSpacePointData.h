// $Id$
#ifndef ALIHLTITSSPACEPOINTDATA_H
#define ALIHLTITSSPACEPOINTDATA_H

//#include "AliHLTTPCRootTypes.h"

/**
 * @struct AliHLTITSSpacePointData
 * Primitive data exchange structure for ITS clusters.
 * Together with the AliHLTITSClusterDataFormat this defines
 * the output of the ITS online Cluster Finder.
 *
 * @ingroup alihlt_its_datastructs
 */
struct AliHLTITSSpacePointData{
  Float_t fY;         // Y coordinate in local coordinates  GetY()
  Float_t fZ;         // Z coordinate in local coordinates  GetZ()
  Float_t fSigmaY2;   // error  of the clusters             GetSigmaY2()
  Float_t fSigmaZ2;   // error  of the clusters             GetSigmaZ2()
  Float_t fSigmaYZ;   // error  of the clusters             GetSigmaYZ()
  Float_t fQ;         // Q of cluster (in ADC counts)       GetQ()
  Int_t fNy;          //number of digits in Y direction     GetNy()
  Int_t fNz;          //number of digits in Z direction     GetNz()
  Int_t fLayer;       // layer number                       GetLayer()
  Int_t fTracks[3];   // MC label                           GetLabel(i)
  Int_t fIndex;       // Detector Index                     GetDetectorIndex()
};
typedef struct AliHLTITSSpacePointData AliHLTITSSpacePointData;


#endif /* ALIHLTITSSPACEPOINTDATA_H */
