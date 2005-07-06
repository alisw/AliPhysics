#ifndef ALIITSCLUSTERFINDERV2SSD_H
#define ALIITSCLUSTERFINDERV2SSD_H
//--------------------------------------------------------------
//                       ITS clusterer V2 for SSD
//
//   This can be a "wrapping" for the V1 cluster finding classes
//   if compiled with uncommented "#define V1" line 
//   in the AliITSclustererV2.cxx file.
//
//   Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch 
//--------------------------------------------------------------
#include "AliITSClusterFinderV2.h"

class TClonesArray;
class AliRawReader;
class AliITSRawStream;

class AliITSClusterFinderV2SSD : public AliITSClusterFinderV2 {
public:
  AliITSClusterFinderV2SSD(AliITSgeom* geom);
  virtual ~AliITSClusterFinderV2SSD(){;}
  virtual void FindRawClusters(Int_t mod);
  virtual void RawdataToClusters(AliRawReader* rawReader,TClonesArray** clusters);
 protected:

  void FindClustersSSD(TClonesArray *digits);
  void FindClustersSSD(Ali1Dcluster* neg, Int_t nn, 
		       Ali1Dcluster* pos, Int_t np,
		       TClonesArray *clusters=0x0);

  void FindClustersSSD(AliITSRawStream* input,TClonesArray** clusters);

  Int_t fLastSSD1;        //index of the last SSD1 detector   
  Float_t fYpitchSSD;     //strip pitch (cm)
  Float_t fHwSSD;         //half-width of an SSD detector (cm)
  Float_t fHlSSD;         //half-length of an SSD detector (cm)
  Float_t fTanP;          //tangent of the stereo angle on the P side
  Float_t fTanN;          //tangent of the stereo angle on the N side


  ClassDef(AliITSClusterFinderV2SSD,1)  // ITS cluster finder V2 for SDD
};

#endif
