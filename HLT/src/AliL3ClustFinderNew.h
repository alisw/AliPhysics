// @(#) $Id$

#ifndef AliL3_ClustFinderNew
#define AliL3_ClustFinderNew

#include "AliL3RootTypes.h"


struct AliClusterData
{
  UInt_t fTotalCharge;
  UInt_t fPad;
  UInt_t fTime;
  ULong64_t fPad2;     //for error in XY direction
  ULong64_t fTime2;    //for error in Z  direction
  UInt_t fMean;
  UInt_t fFlags;
  UInt_t fChargeFalling; //for deconvolution
  UInt_t fLastCharge;    //for deconvolution
  UInt_t fLastMergedPad; //dont merge twice per pad
};
typedef struct AliClusterData AliClusterData;

class AliL3DigitRowData;
class AliL3SpacePointData;

class AliL3ClustFinderNew {
 
 private:
  
  AliL3DigitRowData *fDigitRowData; //!
  AliL3SpacePointData *fSpacePointData; //!
  Bool_t fDeconvTime; //deconv in time direction
  Bool_t fDeconvPad;  //deconv in pad direction
  Bool_t fStdout;     //have print out in write clusters
  Bool_t fCalcerr;    //calculate centroid sigmas
  Bool_t fRawSP;      //store centroids in raw system

  UInt_t fNDigitRowData; //ndigts on row
  Int_t fFirstRow;       //first row
  Int_t fLastRow;        //last row
  Int_t fCurrentRow;     //current active row
  Int_t fCurrentSlice;   //current slice
  Int_t fCurrentPatch;   //current patch
  Int_t fMatch;          //size of match
  UInt_t fThreshold;     //threshold for clusters
  Int_t fNClusters;      //number of found clusters
  Int_t fMaxNClusters;   //max. number of clusters
  Float_t fXYErr;        //fixed error in XY
  Float_t fZErr;         //fixed error in Z

#ifdef do_mc
  void GetTrackID(Int_t pad,Int_t time,Int_t *trackID);
#endif
  
 public:
  AliL3ClustFinderNew();
  virtual ~AliL3ClustFinderNew();
  
  void Read(UInt_t ndigits,AliL3DigitRowData *ptr);
  void InitSlice(Int_t slice,Int_t patch,Int_t firstrow, Int_t lastrow,Int_t maxpoints);
  void InitSlice(Int_t slice,Int_t patch,Int_t maxpoints);
  void ProcessDigits();
  void ProcessRow(AliL3DigitRowData *tempPt);
  void SetOutputArray(AliL3SpacePointData *pt);
  void WriteClusters(Int_t n_clusters,AliClusterData *list);

  void SetXYError(Float_t f) {fXYErr=f;}
  void SetZError(Float_t f) {fZErr=f;}
  void SetDeconv(Bool_t f) {fDeconvPad=f; fDeconvTime=f;}
  void SetThreshold(UInt_t i) {fThreshold=i;}
  void SetMatchWidth(UInt_t i) {fMatch=i;}
  void SetSTDOutput(Bool_t f=kFALSE) {fStdout=f;}  
  void SetCalcErr(Bool_t f=kTRUE) {fCalcerr=f;}
  void SetRawSP(Bool_t f=kFALSE) {fRawSP=f;}
  Int_t GetNumberOfClusters() const {return fNClusters;}
  
  ClassDef(AliL3ClustFinderNew,1) //Fast cluster finder

};

#endif


