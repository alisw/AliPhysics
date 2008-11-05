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
#include "AliITSClusterFinder.h"
#include "AliITSDetTypeRec.h"

class TClonesArray;
class AliRawReader;
class AliITSRawStream;
class AliITSRawStreamSSD;
class AliITSCalibrationSSD;

class AliITSClusterFinderV2SSD : public AliITSClusterFinder {
public:
  AliITSClusterFinderV2SSD(AliITSDetTypeRec* dettyp);
  virtual ~AliITSClusterFinderV2SSD(){;}
  virtual void FindRawClusters(Int_t mod);
  virtual void RawdataToClusters(AliRawReader* rawReader,TClonesArray** clusters);
 protected:
  AliITSClusterFinderV2SSD(const AliITSClusterFinderV2SSD& cf);
  AliITSClusterFinderV2SSD& operator=(const AliITSClusterFinderV2SSD&  cf );
  void FindClustersSSD(TClonesArray *digits);
  void FindClustersSSD(Ali1Dcluster* neg, Int_t nn, 
		       Ali1Dcluster* pos, Int_t np,
		       TClonesArray *clusters=0x0);

  void FindClustersSSD(AliITSRawStreamSSD* input,TClonesArray** clusters);
  virtual AliITSCalibrationSSD* GetResp(Int_t mod)const{
    return (AliITSCalibrationSSD*) fDetTypeRec->GetCalibrationModel(mod);}

  Int_t fLastSSD1;        //index of the last SSD1 detector   
  static Short_t* fgPairs;       //array used to build positive-negative pairs
  static Int_t    fgPairsSize;    //actual size of pairs array

  ClassDef(AliITSClusterFinderV2SSD,2)  // ITS cluster finder V2 for SDD
};

#endif
