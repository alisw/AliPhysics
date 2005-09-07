// @(#) $Id$

#ifndef AliHLTTPC_ClustFinderNew
#define AliHLTTPC_ClustFinderNew

#include "AliHLTTPCRootTypes.h"


struct ClusterData
{
  UInt_t fTotalCharge;
  UInt_t fPad;
  UInt_t fTime;
  UInt_t fPad2;     //for error in XY direction
  ULong64_t fTime2; //for error in Z  direction
  UInt_t fMean;
  UInt_t fFlags;
  UInt_t fChargeFalling; //for deconvolution
  UInt_t fLastCharge;    //for deconvolution
  UInt_t fLastMergedPad; //dont merge twice per pad
};
typedef struct ClusterData ClusterData;

class AliHLTTPCDigitRowData;
class AliHLTTPCSpacePointData;

class AliHLTTPCClustFinderNew {
 
 private:
  
  AliHLTTPCDigitRowData *fDigitRowData; //!
  AliHLTTPCSpacePointData *fSpacePointData; //!
  Bool_t fDeconvTime; //deconv in time direction
  Bool_t fDeconvPad;  //deconv in pad direction
  Bool_t fStdout;     //have print out in write clusters
  Bool_t fCalcerr;    //calculate centroid sigmas
  Bool_t fRawSP;      //store centroids in raw system

  UInt_t fNDigitRowData;
  Int_t fFirstRow;
  Int_t fLastRow;
  Int_t fCurrentRow;
  Int_t fCurrentSlice;
  Int_t fCurrentPatch;
  Int_t fMatch;
  UInt_t fThreshold;
  Int_t fNClusters;
  Int_t fMaxNClusters;
  Float_t fXYErr;
  Float_t fZErr;

#ifdef do_mc
  void GetTrackID(Int_t pad,Int_t time,Int_t *trackID);
#endif
  
 public:
  AliHLTTPCClustFinderNew();
  virtual ~AliHLTTPCClustFinderNew();
  
  void Read(UInt_t ndigits,AliHLTTPCDigitRowData *ptr);
  void InitSlice(Int_t slice,Int_t patch,Int_t firstrow, Int_t lastrow,Int_t maxpoints);
  void InitSlice(Int_t slice,Int_t patch,Int_t maxpoints);
  void ProcessDigits();
  void ProcessRow(AliHLTTPCDigitRowData *tempPt);
  void SetOutputArray(AliHLTTPCSpacePointData *pt);
  void WriteClusters(Int_t n_clusters,ClusterData *list);

  void SetXYError(Float_t f) {fXYErr=f;}
  void SetZError(Float_t f) {fZErr=f;}
  void SetDeconv(Bool_t f) {fDeconvPad=f; fDeconvTime=f;}
  void SetThreshold(UInt_t i) {fThreshold=i;}
  void SetMatchWidth(UInt_t i) {fMatch=i;}
  void SetSTDOutput(Bool_t f=kFALSE) {fStdout=f;}  
  void SetCalcErr(Bool_t f=kTRUE) {fCalcerr=f;}
  void SetRawSP(Bool_t f=kFALSE) {fRawSP=f;}
  Int_t GetNumberOfClusters() {return fNClusters;}
  
  ClassDef(AliHLTTPCClustFinderNew,1) //Fast cluster finder

};

#endif


