#ifndef AliL3_ClustFinderNew
#define AliL3_ClustFinderNew

#include "AliL3RootTypes.h"


struct ClusterData
{
  UInt_t fTotalCharge;
  UInt_t fPad;
  UInt_t fTime;
  UInt_t fPad2;  //for error in XY direction
  UInt_t fTime2; //for error in Z  direction
  UInt_t fMean;
  UInt_t fFlags;
  UInt_t fChargeFalling; //for deconvolution
  UInt_t fLastCharge;    //for deconvolution
  UInt_t fLastMergedPad; //dont merge twice per pad
};
typedef struct ClusterData ClusterData;

class AliL3DigitRowData;
class AliL3SpacePointData;

class AliL3ClustFinderNew {
 
 private:
  
  AliL3DigitRowData *fDigitRowData; //!
  AliL3SpacePointData *fSpacePointData; //!
  Bool_t fDeconvTime;
  Bool_t fDeconvPad;
  Bool_t fstdout;
  Bool_t fcalcerr;

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
  AliL3ClustFinderNew();
  virtual ~AliL3ClustFinderNew();
  
  void Read(UInt_t ndigits,AliL3DigitRowData *ptr);
  void InitSlice(Int_t slice,Int_t patch,Int_t firstrow, Int_t lastrow,Int_t maxpoints);
  void InitSlice(Int_t slice,Int_t patch,Int_t maxpoints);
  void ProcessDigits();
  void ProcessRow(AliL3DigitRowData *tempPt);
  void SetOutputArray(AliL3SpacePointData *pt);
  void WriteClusters(Int_t n_clusters,ClusterData *list);

  void SetXYError(Float_t f) {fXYErr=f;}
  void SetZError(Float_t f) {fZErr=f;}
  void SetDeconv(Bool_t f) {fDeconvPad=f; fDeconvTime=f;}
  void SetThreshold(UInt_t i) {fThreshold=i;}
  void SetMatchWidth(UInt_t i) {fMatch=i;}
  void SetSTDOutput(Bool_t f=kFALSE) {fstdout=f;}  
  void SetCalcErr(Bool_t f=kTRUE) {fcalcerr=f;}

  Int_t GetNumberOfClusters() {return fNClusters;}
  
  ClassDef(AliL3ClustFinderNew,1) //Fast cluster finder

};

#endif


