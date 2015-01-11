#ifndef ALIITSCLUSTERFINDERSDDFAST_H
#define ALIITSCLUSTERFINDERSDDFAST_H

//----------------------------------------------------------------------
//              ITS clusterer for SDD - fast algorithm
//
//   Origin: Simone Capodicasa, Universita e INFN, capodica@to.infn.it
//----------------------------------------------------------------------

#include "AliITSClusterFinder.h"
#include "AliITSDetTypeRec.h"
#include <vector>

class TBits;
class TClonesArray;
class AliRawReader;
class AliITSRawStream;
class AliITSCalibrationSDD;
class AliITSsegmentationSDD;

class AliITSClusterFinderSDDfast : public AliITSClusterFinder {
 public:
  AliITSClusterFinderSDDfast(AliITSDetTypeRec* dettyp);
  virtual ~AliITSClusterFinderSDDfast();
  virtual void FindRawClusters(Int_t mod);
  virtual void RawdataToClusters(AliRawReader* rawReader);
  void SetPeakSelection(Float_t looseCut=15., Float_t tightCut=30., Float_t maxTime=2000.){
    fCutOnPeakLoose=looseCut;
    fCutOnPeakTight=tightCut;
    fMaxDrTimeForTightCut=maxTime;
  }
  
  enum {kHybridsPerDDL = 24};   // number of hybrids in each DDL
  enum {kModulesPerDDL = 12};   // number of modules in each DDL
  
 protected:
  AliITSClusterFinderSDDfast(const AliITSClusterFinderSDDfast &source); // copy constructor
  // assignment operator
  AliITSClusterFinderSDDfast& operator=(const AliITSClusterFinderSDDfast &source);
  void FindClustersSDD(TClonesArray *digits);
  void FindClustersSDD(std::vector<int>& bins0, std::vector<int>& bins1, const Int_t map0[], const Int_t map1[], TClonesArray *dig, TClonesArray *clusters=0x0, Int_t jitter=0);
  
  void FindClustersSDD(AliITSRawStream* input);
  virtual AliITSCalibrationSDD* GetResp(Int_t mod)const{
    return (AliITSCalibrationSDD*) fDetTypeRec->GetCalibrationModel(mod);}
  virtual AliITSsegmentationSDD* GetSeg()const{
    return (AliITSsegmentationSDD*)fDetTypeRec->GetSegmentationModel(1);}
  
  Int_t fNAnodes;                   // number of anodes
  Int_t fNTimeBins;                 // number of time bins
  Int_t fNZbins;                    // number of cells along anodes
  Int_t fNXbins;                    // number of cells along time
  std::vector<std::vector<int> > fDDLBins; // container for digits for 1 DDL
  Float_t fCutOnPeakLoose;          // loose cut on peak (for all drift times)
  Float_t fCutOnPeakTight;          // tight cut on peak (for small drift times)
  Float_t fMaxDrTimeForTightCut;    // max. drift time for fCutOnPeakTight
  
  ClassDef(AliITSClusterFinderSDDfast,1)  // ITS cluster finder fast for SDD
    };

#endif
