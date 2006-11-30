// @(#) $Id$
// Author: Constantin Loizides <mailto:loizides@ikf.uni-frankfurt.de>
//- Copyright & copy ALICE HLT Group
// See the implementation file for the detailed description

#ifndef AliHLTVHDLClusterFinder_H
#define AliHLTVHDLClusterFinder_H

#include "AliHLTAltroMemHandler.h"

struct VHDLClusterData
{
  UInt_t fTotalCharge; //total charge
  UInt_t fPad;   //mean in pad
  UInt_t fTime;  //mean in time
  UInt_t fPad2;  //for error in XY direction
  UInt_t fTime2; //for error in Z  direction
  UInt_t fMean;  //mean for comparism
  UInt_t fMerge; //number of merges
  UShort_t fRow;     //row of cluster
  UShort_t fLastPad; //last pad on merge
  UInt_t fChargeFalling; //for deconvolution
  UInt_t fLastCharge;    //for deconvolution
};
typedef struct VHDLClusterData VCData;

//size of ring buffer
#define N_mem 2500
//size of cluster list
#define N_clmem 5000


class AliHLTVHDLClusterFinder 
{
 public:
  AliHLTVHDLClusterFinder();
  virtual ~AliHLTVHDLClusterFinder();
  
  void ProcessDigits();

  void SetXYError(Float_t f) {fXYErr=f;}
  void SetZError(Float_t f) {fZErr=f;}
  void SetDeconv(Bool_t f) {fDeconvPad=f; fDeconvTime=f;}
  void SetThreshold(UInt_t i=10) {fThreshold=i;}
  void SetMatchWidth(UInt_t i=4) {fMatch=i;}
  void SetMergeMinimum(UInt_t i=1) {fMinMerge=i;}
  void SetSTDOutput(Bool_t f=kFALSE) {fstdout=f;}  
  void SetCalcErr(Bool_t f=kTRUE) {fcalcerr=f;}
  void SetASCIIInput(FILE *f){fAltromem.SetASCIIInput(f);}

  Int_t GetNumberOfClusters() const {return fNClusters;}
  
 private:
  AliHLTAltroMemHandler fAltromem; //! pointer to the ALTRO memory?
  VCData fSeq; //! C-structure containing the data sequence
  VCData fSeqs[N_clmem]; //! array of C-structures containg the data sequence
  UShort_t fPList[N_mem]; // pad list
  UShort_t fRow,fNRow; // current row and number of rows
  UChar_t  fPad,fNPad; // current pad and number of pads
  UShort_t fRP,fWP,fOP,fEP,fFP; //pointer in ringbuffer
  UShort_t fLast,fFirst;        //free area in memory

  Int_t fTC; //totalcharge
  Int_t fMT; //mean in time
  Int_t fST; //sigma in time
  Int_t fSM; //seq. mean

  Bool_t fDeconvTime; // deconvoluted time
  Bool_t fDeconvPad;  // deconvoluted pad
  Bool_t fstdout;     // output flag
  Bool_t fcalcerr;    // flag to calculate errors
  Float_t fXYErr;     // error in XY
  Float_t fZErr;      // error in Z

  Int_t fMatch;       //match distance
  UInt_t fThreshold;  //threshold for cluster
  UInt_t fMinMerge;   //minimum number of merges for cluster
  Int_t fNClusters;   //number of found clusters

#ifdef VHDLDEBUG
  FILE *fdeb; //! file for debug
#endif

  void Clear();
  void ClearSeq(UShort_t i);
  void FreeSeq(UShort_t i);
  void IncPointer(UShort_t &p,Short_t add=1,UShort_t N=N_mem);
  void IncWPointer();
  void IncRPointer();
  void NextFreeIndex();
  void FlushMemory();
  void PrepareMemory();
  void OutputMemory();
  void CompareSeq();
  void MergeSeq();
  void InsertSeq();
  void MakeSequence();
  void ProcessSequence();
  //void WriteClusters(Int_t n_clusters,ClusterData *list);

  ClassDef(AliHLTVHDLClusterFinder,1)

};

typedef AliHLTVHDLClusterFinder AliL3VHDLClusterFinder; // for backward compatibility

#endif
