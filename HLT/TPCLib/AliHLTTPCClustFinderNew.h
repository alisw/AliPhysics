// @(#) $Id$
// Original: AliHLTClustFinderNew.h,v 1.13 2004/06/18 10:55:26 loizides 
#ifndef AliHLTTPC_ClustFinderNew
#define AliHLTTPC_ClustFinderNew

class AliHLTTPCDigitRowData;
class AliHLTTPCSpacePointData;

class AliHLTTPCClustFinderNew {

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
  };
  typedef struct AliClusterData AliClusterData; //!

 private:
  AliHLTTPCDigitRowData *fDigitRowData; //!
  AliHLTTPCSpacePointData *fSpacePointData; //!
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
  AliHLTTPCClustFinderNew();
  virtual ~AliHLTTPCClustFinderNew();
  
  void Read(UInt_t ndigits,AliHLTTPCDigitRowData *ptr);
  void InitSlice(Int_t slice,Int_t patch,Int_t firstrow, Int_t lastrow,Int_t maxpoints);
  void InitSlice(Int_t slice,Int_t patch,Int_t maxpoints);
  void ProcessDigits();
  void ProcessRow(AliHLTTPCDigitRowData *tempPt);
  void SetOutputArray(AliHLTTPCSpacePointData *pt);
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
  
  ClassDef(AliHLTTPCClustFinderNew,1) //Fast cluster finder
};
#endif


