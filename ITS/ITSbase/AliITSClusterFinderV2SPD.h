#ifndef ALIITSCLUSTERFINDERV2SPD_H
#define ALIITSCLUSTERFINDERV2SPD_H
//--------------------------------------------------------------
//                       ITS clusterer V2 for SPD
//
//   This can be a "wrapping" for the V1 cluster finding classes
//   if compiled with uncommented "#define V1" line 
//   in the AliITSclustererV2.cxx file.
//
//   Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch 
//--------------------------------------------------------------
#include "AliITSClusterFinder.h"

class TClonesArray;
class AliRawReader;
class AliITSRawStream;
class AliITSRawStreamSPD;

class AliITSClusterFinderV2SPD : public AliITSClusterFinder {
public:
  AliITSClusterFinderV2SPD(AliITSDetTypeRec* dettyp);
  virtual ~AliITSClusterFinderV2SPD(){;}
  virtual void FindRawClusters(Int_t mod);
  virtual void RawdataToClusters(AliRawReader* rawReader);
  

 protected:

  void FindClustersSPD(TClonesArray *digits);
  void FindClustersSPD(AliITSRawStreamSPD* input);
  Int_t ClustersSPD(AliBin* bins, TClonesArray* digits,TClonesArray* clusters,Int_t maxBins, Int_t nzbins,Int_t iModule,Bool_t rawdata=kFALSE);

  Int_t fLastSPD1;       //index of the last SPD1 detector
  Int_t fNySPD;          //number of pixels in Y
  Int_t fNzSPD;          //number of pixels in Z
  Float_t fYpitchSPD;    //pixel size in Y
  Float_t fZ1pitchSPD,fZ2pitchSPD;    //pixel sizes in Z
  Float_t fHwSPD;        //half width of the SPD detector
  Float_t fHlSPD;        //half length of the SPD detector
  Float_t fYSPD[260];    //Y-coordinates of pixel centers
  Float_t fZSPD[170];    //Z-coordinates of pixel centers

  ClassDef(AliITSClusterFinderV2SPD,1)  // ITS cluster finder V2 for SPD
};

#endif
