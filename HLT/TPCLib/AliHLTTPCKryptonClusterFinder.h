// $Id$

#ifndef AliHLTTPC_KRYPTONCLUSTERFINDER
#define AliHLTTPC_KRYPTONCLUSTERFINDER
/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTTPCKryptonClusterFinder.h
    @author Anders Vestbo, Constantin Loizides, Jochen Thaeder
	    Kenneth Aamodt kenneth.aamodt@student.uib.no
    @date   
    @brief  Cluster Finder for the TPC
*/

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt


#include "AliHLTLogging.h"
#include "AliHLTTPCPad.h"
class AliHLTTPCSpacePointData;
class AliHLTTPCDigitReader;

class AliHLTTPCKryptonClusterFinder : public AliHLTLogging {

 public:
  struct AliClusterData
  {
    UInt_t fTotalCharge;   //tot charge of cluster
    UInt_t fPad;           //pad value
    UInt_t fTime;          //time value
    ULong64_t fPad2;       //for error in XY direction
    ULong64_t fTime2;      //for error in Z  direction
    UInt_t fMean;          //mean in time
    UInt_t fFlags;         //different flags
    UInt_t fChargeFalling; //for deconvolution
    UInt_t fLastCharge;    //for deconvolution
    UInt_t fLastMergedPad; //dont merge twice per pad
    Int_t fRow;             //row value
  };
  typedef struct AliClusterData AliClusterData; //!

  /** standard constructor */
  AliHLTTPCKryptonClusterFinder();
  /** destructor */
  virtual ~AliHLTTPCKryptonClusterFinder();

  void InitSlice(Int_t slice,Int_t patch,Int_t firstrow, Int_t lastrow,Int_t maxpoints);

  void SetReader(AliHLTTPCDigitReader* f){fDigitReader = f;}

  void PrintClusters();

  Int_t GetNumberOfClusters() const {return fNClusters;}

  void ReBunch(const UInt_t * bunchData,Int_t bunchSize);
  void ReadDataUnsorted(void* ptr,unsigned long size);
  void FindNormalClusters();
  void FindKryptonClusters();
  void CheckForCandidateOnPreviousPad(AliHLTTPCClusters* tmpCluster);
  void SetPatch(Int_t patch){fCurrentPatch=patch;}
  void InitializePadArray();
  Int_t DeInitializePadArray();
  Bool_t ComparePads(AliHLTTPCPad *nextPad,AliHLTTPCClusters* candidate,Int_t nextPadToRead);
 private: 
  /** copy constructor prohibited */
  AliHLTTPCKryptonClusterFinder(const AliHLTTPCKryptonClusterFinder&);
  /** assignment operator prohibited */
  AliHLTTPCKryptonClusterFinder& operator=(const AliHLTTPCKryptonClusterFinder&);

  AliHLTTPCSpacePointData *fSpacePointData; //! array of space points
  AliHLTTPCDigitReader *fDigitReader;       //! reader instance

  UChar_t* fPtr;   //! pointer to packed block
  unsigned long fSize; //packed block size

  Int_t fFirstRow;       //first row
  Int_t fLastRow;        //last row
  Int_t fCurrentRow;     //current active row
  Int_t fCurrentSlice;   //current slice
  Int_t fCurrentPatch;   //current patch
  Int_t fMatch;          //size of match
  UInt_t fThreshold;      //Threshold on total charge for krypton cluster

  Int_t fNClusters;      //number of found clusters
  Int_t fMaxNClusters;   //max. number of clusters
  Float_t fXYErr;        //fixed error in XY
  Float_t fZErr;         //fixed error in Z
  
  Bool_t fVectorInitialized;

  typedef vector<AliHLTTPCPad*> AliHLTTPCPadVector;

  vector<AliHLTTPCPadVector> fRowPadVector;                        //! transient

  vector<AliHLTTPCClusters> fClusters;                             //! transient
  
  vector<Int_t> fTimebinsInBunch;                                  //! transient

  vector<Int_t> fIndexOfBunchStart;                                //! transient

  UInt_t* fNumberOfPadsInRow;                                      //! transient
  
  UInt_t fNumberOfRows;                                            //! transient
  
  UInt_t fRowOfFirstCandidate;

  ClassDef(AliHLTTPCKryptonClusterFinder,0) //Fast cluster finder
};
#endif
