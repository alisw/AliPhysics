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
#include "AliITSClusterFinder.h"
#include "AliITSDetTypeRec.h"

class TBits;
class TClonesArray;
class AliRawReader;
class AliITSRawStream;
class AliITSCalibrationSDD;
class AliITSsegmentationSDD;

class AliITSClusterFinderV2SDD : public AliITSClusterFinder {
public:
  AliITSClusterFinderV2SDD(AliITSDetTypeRec* dettyp);
  virtual ~AliITSClusterFinderV2SDD();
  virtual void FindRawClusters(Int_t mod);
  virtual void RawdataToClusters(AliRawReader* rawReader,TClonesArray** clusters);

  enum {kHybridsPerDDL = 24};   // number of hybrids in each DDL 
  enum {kModulesPerDDL = 12};   // number of modules in each DDL 

 protected:
 AliITSClusterFinderV2SDD(const AliITSClusterFinderV2SDD &source); // copy constructor
  // assignment operator
  AliITSClusterFinderV2SDD& operator=(const AliITSClusterFinderV2SDD &source);
  Bool_t NoiseSuppress(Int_t k, Int_t sid, AliBin* bins, AliITSCalibrationSDD* cal) const;
  void FindClustersSDD(TClonesArray *digits);
  void FindClustersSDD(AliBin* bins[2], TBits* anodeFired[2],
		       TClonesArray *dig, TClonesArray *clusters=0x0, Int_t jitter=0);

  void FindClustersSDD(AliITSRawStream* input,TClonesArray** clusters);
  virtual AliITSCalibrationSDD* GetResp(Int_t mod)const{
    return (AliITSCalibrationSDD*) fDetTypeRec->GetCalibrationModel(mod);}
  virtual AliITSsegmentationSDD* GetSeg()const{
    return (AliITSsegmentationSDD*)fDetTypeRec->GetSegmentationModel(1);} 

  Int_t fNAnodes;                   // number of anodes
  Int_t fNTimeBins;                 // number of time bins
  Int_t fNZbins;                    // number of cells along anodes
  Int_t fNXbins;                    // number of cells along time
  AliBin* fDDLBins[kHybridsPerDDL]; // container for digits for 1 DDL

  ClassDef(AliITSClusterFinderV2SDD,5)  // ITS cluster finder V2 for SDD
};

#endif
