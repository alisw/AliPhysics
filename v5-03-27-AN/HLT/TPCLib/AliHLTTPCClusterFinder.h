//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTTPCCLUSTERFINDER_H
#define ALIHLTTPCCLUSTERFINDER_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

//  @file   AliHLTTPCClusterFinder.h
//  @author Anders Vestbo, Constantin Loizides
// 	    Kenneth Aamodt kenneth.aamodt@student.uib.no
//  @brief  HLT Cluster Finder for the TPC
//  @note

#include "AliHLTLogging.h"
#include <vector>
#include "AliHLTTPCTransform.h"
#include "AliHLTTPCDigitData.h"
#include "AliHLTTPCDigitReader.h"
#include "AliTPCRecoParam.h"
#include "AliHLTTPCClusterMCData.h"
#include "AliHLTTPCRawCluster.h"

class AliHLTTPCPad;
struct AliHLTTPCSpacePointData;
class AliHLTTPCClusters;
class AliTPCTransform;
class AliTPCParam;

/**
 * @class AliHLTTPCClusterFinder
 *
 * The current cluster finder for HLT
 * (Based on STAR L3)
 *
 * Basically we have two versions for the cluster finder now.
 * The default version, reads the data pad by pad, and find the
 * clusters as it reads the data. The other version has now been
 * developed to cope with unsorted data. New methods for the unsorted
 * version can  be found at the end of the default one i the source file.
 * Currently the new version is only build to manage zero-suppressed data.
 * More functionality will be added later.
 * 
 * The cluster finder is initialized with the Init function, 
 * providing the slice and patch information to work on. 
 *
 * The input is a provided by the AliHLTTPCDigitReader class,
 * using the init() funktion, and the next() funktion in order 
 * to get the next bin. Either packed or unpacked data can be
 * processed, dependent if one uses AliHLTTPCDigitReaderPacked 
 * class or AliHLTTPCDigitReaderUnpacked class in the 
 * Clusterfinder Component.
 * The resulting space points will be in the
 * array given by the SetOutputArray function.
 * 
 * There are several setters which control the behaviour:
 *
 * - SetXYError(Float_t):   set fixed error in XY direction
 * - SetZError(Float_t):    set fixed error in Z  direction
 *                            (used if errors are not calculated) 
 * - SetDeconv(Bool_t):     switch on/off deconvolution
 * - SetThreshold(UInt_t):  set charge threshold for cluster
 * - SetMatchWidth(UInt_t): set the match distance in 
 *                            time for sequences to be merged 
 * - SetSTDOutput(Bool_t):  switch on/off output about found clusters   
 * - SetCalcErr(Bool_t):    switch on/off calculation of 
 *                          space point errors (or widths in raw system)
 * - SetRawSP(Bool_t):      switch on/off convertion to raw system
 *
 *
 * Example Usage:
 *
 * <pre>
 * AliHLTTPCFileHandler *file = new AliHLTTPCFileHandler();
 * file->SetAliInput(digitfile); //give some input file
 * for(int slice=0; slice<=35; slice++){
 *   for(int patch=0; pat<6; pat++){
 *     file->Init(slice,patch);
 *     UInt_t ndigits=0;
 *     UInt_t maxclusters=100000;
 *     UInt_t pointsize = maxclusters*sizeof(AliHLTTPCSpacePointData);
 *     AliHLTTPCSpacePointData *points = (AliHLTTPCSpacePointData*)memory->Allocate(pointsize);
 *     AliHLTTPCDigitRowData *digits = (AliHLTTPCDigitRowData*)file->AliAltroDigits2Memory(ndigits,event);
 *     AliHLTTPCClusterFinder *cf = new AliHLTTPCClusterFinder();
 *     cf->SetMatchWidth(2);
 *     cf->InitSlice( slice, patch, row[0], row[1], maxPoints );
 *     cf->SetSTDOutput(kTRUE);    //Some output to standard IO
 *     cf->SetRawSP(kFALSE);       //Convert space points to local system
 *     cf->SetThreshold(5);        //Threshold of cluster charge
 *     cf->SetDeconv(kTRUE);       //Deconv in pad and time direction
 *     cf->SetCalcErr(kTRUE);      //Calculate the errors of the spacepoints
 *     cf->SetOutputArray(points); //Move the spacepoints to the array
 *     cf->Read(iter->fPtr, iter->fSize ); //give the data to the cf
 *     cf->ProcessDigits();        //process the rows given by init
 *     Int_t npoints = cf->GetNumberOfClusters();
 *     AliHLTTPCMemHandler *out= new AliHLTTPCMemHandler();
 *     out->SetBinaryOutput(fname);
 *     out->Memory2Binary(npoints,points); //store the spacepoints
 *     out->CloseBinaryOutput();
 *     delete out;
 *     file->free();
 *     delete cf;
 *   }
 * }
 * </pre>
 * @ingroup alihlt_tpc
 */
class AliHLTTPCClusterFinder : public AliHLTLogging {

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
    UInt_t fQMax;           //qmax
  };
  typedef struct AliClusterData AliClusterData; //!

  static Bool_t CompareWeights( const AliHLTTPCClusterMCWeight &mc1,  const AliHLTTPCClusterMCWeight &mc2 ){ return mc1.fWeight > mc2.fWeight; }
  

  /** standard constructor */
  AliHLTTPCClusterFinder();

  /** destructor */
  virtual ~AliHLTTPCClusterFinder();

  /** Initialize the slice */
  void InitSlice(Int_t slice,Int_t patch,Int_t maxpoints);

  /** Initializes the pad array (vector)*/
  void InitializePadArray();

  /** Deinitialize the pad array (vector)*/
  Int_t DeInitializePadArray();

  /** Read the data in unsorted format, storing the clustercandidates */
  void ReadDataUnsorted(void* ptr,unsigned long size);

  /** Read the data in unsorted format, and deconvolute the signals for each pad in time direction */
  void ReadDataUnsortedDeconvoluteTime(void* ptr,unsigned long size);

  /** Loops over all rows finding the clusters */
  void FindClusters();

  /** Compare two neighbouring pads for matching clustercandidates */
  Bool_t ComparePads(AliHLTTPCPad *nextPad,AliHLTTPCClusters* candidate,Int_t nextPadToRead);
  
  /**  Fills the hw address list */
  Int_t FillHWAddressList(AliHLTUInt16_t *hwaddlist, Int_t maxHWAddress);

 /**  Fills the mc info */
  Int_t FillOutputMCInfo(AliHLTTPCClusterMCLabel * outputMCInfo, Int_t maxNumberOfClusterMCInfo);

  Int_t FillOutputRaw(AliHLTTPCRawCluster* rawClusters, unsigned sizeInByte) const;

  /** Set the pointer to the outputbuffer */
  void SetOutputArray(AliHLTTPCSpacePointData *pt);

  /** Returns the number of clusters */
  Int_t GetNumberOfClusters() const {return fNClusters;}

  /** Returns the Ocuppancy limit */
  Float_t GetOccupancyLimit() const {return fOccupancyLimit;}
  
  // setters
  void SetDeconv(Bool_t f) {fDeconvPad=f; fDeconvTime=f;}
  void SetDeconvPad(Bool_t f) {fDeconvPad=f;}
  void SetDeconvTime(Bool_t f) {fDeconvTime=f;}
  void SetUnsorted(Int_t unsorted){fUnsorted=unsorted;}
  void SetPatch(Int_t patch){fCurrentPatch=patch;}
  void SetDoPadSelection(Bool_t input){fDoPadSelection=input;}
  void SetLastTimeBin(Int_t ltb){fLastTimeBin=ltb;}
  void SetFirstTimeBin(Int_t ftb){fFirstTimeBin=ftb;}
  void SetReleaseMemory( Bool_t v ){ fReleaseMemory = v;}
  void UpdateLastTimeBin(){fLastTimeBin=AliHLTTPCTransform::GetNTimeBins();}

//---------------------------------- Under this line the old sorted clusterfinder functions can be found --------------------------------
  void Read(void* ptr,unsigned long size);
  void ProcessDigits();
  void WriteClusters(Int_t n_clusters,AliClusterData *list);
  void WriteClusters(Int_t nclusters,AliHLTTPCClusters *list);
  void PrintClusters();
  void SetXYError(Float_t f) {fXYErr=f;}
  void SetZError(Float_t f) {fZErr=f;}
  void SetThreshold(UInt_t i) {fThreshold=i;}
  void SetOccupancyLimit(Float_t f) {fOccupancyLimit=f;}
  void SetMatchWidth(UInt_t i) {fMatch=i;}
  void SetSTDOutput(Bool_t f=kFALSE) {fStdout=f;}  
  void SetCalcErr(Bool_t f=kTRUE) {fCalcerr=f;}
  void SetFillRawClusters(Bool_t f=kFALSE) {fFillRawClusters=f;}
  void SetReader(AliHLTTPCDigitReader* f){fDigitReader = f;}

  void Set32BitFormat(Bool_t flag){f32BitFormat = flag;}

  void SetDoMC(Bool_t flag){fDoMC = flag;}
  
  vector<AliHLTUInt16_t> fClustersHWAddressVector;  //! transient
  
  typedef vector<AliHLTTPCPad*> AliHLTTPCPadVector;
  
  vector<AliHLTTPCPadVector> fRowPadVector;	 //! transient
  
  void FillMCClusterVector(vector<AliHLTTPCDigitData> *digitData);

  vector<AliHLTTPCClusterMCWeight> GetClusterMCInfo() const {return fClusterMCVector;}

  Bool_t UpdateCalibDB();

 protected: 
  /** copy constructor prohibited */
  AliHLTTPCClusterFinder(const AliHLTTPCClusterFinder&);
  /** assignment operator prohibited */
  AliHLTTPCClusterFinder& operator=(const AliHLTTPCClusterFinder&);

  AliHLTTPCSpacePointData *fSpacePointData; //! array of space points
  AliHLTTPCDigitReader *fDigitReader;       //! reader instance

  UChar_t* fPtr;           //! pointer to packed block
  unsigned long fSize;     //! packed block size
  Bool_t fDeconvTime;      //! deconv in time direction
  Bool_t fDeconvPad;       //! deconv in pad direction
  Bool_t fStdout;          //! have print out in write clusters
  Bool_t fCalcerr;         //! calculate centroid sigmas
  Bool_t fFillRawClusters; //! store centroids in raw system in separate array


  Int_t fFirstRow;       //! first row
  Int_t fLastRow;        //! last row
  Int_t fCurrentRow;     //! current active row
  Int_t fCurrentSlice;   //! current slice
  Int_t fCurrentPatch;   //! current patch
  Int_t fMatch;          //! size of match
  UInt_t fThreshold;     //! threshold for clusters
  Int_t fNClusters;      //! number of found clusters
  Int_t fMaxNClusters;   //! max. number of clusters
  Float_t fXYErr;        //! fixed error in XY
  Float_t fZErr;         //! fixed error in Z
  
  Float_t fOccupancyLimit;    //! Occupancy Limit

  Int_t fUnsorted;            //! enable for processing of unsorted digit data
  Bool_t fVectorInitialized;  //! flag to check if pad vector is initialized
 
  vector<AliHLTTPCClusters> fClusters;                             //! transient

  vector<AliHLTTPCClusterMCLabel> fClustersMCInfo;                           //! transient

  vector<AliHLTTPCDigitData> fMCDigits;                            //! transient

  vector<AliHLTTPCRawCluster> fRawClusters;                        //! transient
  
  UInt_t* fNumberOfPadsInRow;                                      //! transient
  
  UInt_t fNumberOfRows;                                            //! transient
  
  UInt_t fRowOfFirstCandidate;  				   //! transient

  Bool_t fDoPadSelection;					   //! transient

  Int_t fFirstTimeBin;					           //! transient
  
  Int_t fLastTimeBin;					           //! transient

  UInt_t fTotalChargeOfPreviousClusterCandidate;                   //! transient

  Bool_t fChargeOfCandidatesFalling;                               //! transient

  Bool_t f32BitFormat;                                             //! transient

  Bool_t fDoMC;                                                    //! transient

  vector<AliHLTTPCClusterMCWeight> fClusterMCVector;                               //! transient

  AliTPCTransform * fOfflineTransform;                             //! transient

  AliTPCParam   *fOfflineTPCParam;                                 //! transient

  AliTPCRecoParam fOfflineTPCRecoParam;                            //! transient

  Float_t fTimeMeanDiff;                                           //! transient

  Bool_t fReleaseMemory; //! flag to release the memory after each event

#ifdef do_mc
  void GetTrackID(Int_t pad,Int_t time,Int_t *trackID) const;
#endif
  
  ClassDef(AliHLTTPCClusterFinder, 0) //Fast cluster finder
};
#endif
