#ifndef ALIITSCLUSTERFINDERV2SDD_H
#define ALIITSCLUSTERFINDERV2SDD_H
//--------------------------------------------------------------
//                       ITS clusterer V2 for SDD
//
//   This can be a "wrapping" for the V1 cluster finding classes
//   if compiled with uncommented "#define V1" line 
//   in the AliITSclustererV2.cxx file.
//
//   Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch 
//--------------------------------------------------------------
#include "AliITSClusterFinderV2.h"
#include "AliITSDetTypeRec.h"

class TClonesArray;
class AliRawReader;
class AliITSRawStream;
class AliITSCalibrationSDD;
class AliITSsegmentationSDD;

class AliITSClusterFinderV2SDD : public AliITSClusterFinderV2 {
public:
  AliITSClusterFinderV2SDD(AliITSDetTypeRec* dettyp);
  virtual ~AliITSClusterFinderV2SDD(){;}
  virtual void FindRawClusters(Int_t mod);
  virtual void RawdataToClusters(AliRawReader* rawReader,TClonesArray** clusters);

 protected:

  void FindClustersSDD(TClonesArray *digits);
  void FindClustersSDD(AliBin* bins[2], Int_t nMaxBin, Int_t nMaxZ,
		       TClonesArray *dig, TClonesArray *clusters=0x0);

  void FindClustersSDD(AliITSRawStream* input,TClonesArray** clusters);
  void CorrectPosition(Float_t &z, Float_t&y);
  virtual AliITSCalibrationSDD* GetResp(Int_t mod)const{
    return (AliITSCalibrationSDD*) fDetTypeRec->GetCalibrationModel(mod);}
  virtual AliITSsegmentationSDD* GetSeg()const{
    return (AliITSsegmentationSDD*)fDetTypeRec->GetSegmentationModel(1);} 


  ClassDef(AliITSClusterFinderV2SDD,4)  // ITS cluster finder V2 for SDD
};

#endif
