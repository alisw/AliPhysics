#ifndef AliL3_ClustFinderNew
#define AliL3_ClustFinderNew

#include "AliL3RootTypes.h"


struct ClusterData
{

  UInt_t fTotalCharge;
  UInt_t fPad;
  UInt_t fTime;
  UInt_t fMean;
  UInt_t fFlags;
  UInt_t fChargeFalling; //for deconvolution
  UInt_t fLastCharge;    //for deconvolution
};
typedef struct ClusterData ClusterData;

class AliL3DigitRowData;
class AliL3Transform;
class AliL3SpacePointData;

class AliL3ClustFinderNew {
 
 private:
  
  AliL3DigitRowData *fDigitRowData; //!
  AliL3Transform *fTransform; //!
  AliL3SpacePointData *fSpacePointData; //!
  Bool_t fDeconvTime;
  Bool_t fDeconvPad;
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
  
 public:
  AliL3ClustFinderNew();
  AliL3ClustFinderNew(AliL3Transform *transform);
  virtual ~AliL3ClustFinderNew();
  
  void Read(UInt_t ndigits,AliL3DigitRowData *ptr);
  void InitSlice(Int_t slice,Int_t patch,Int_t firstrow, Int_t lastrow,Int_t maxpoints);
  void InitSlice(Int_t slice,Int_t patch,Int_t maxpoints);
  void ProcessDigits();
  void ProcessRow(AliL3DigitRowData *tempPt);
  void SetOutputArray(AliL3SpacePointData *pt);
  void WriteClusters(Int_t n_clusters,ClusterData *list);

  void SetXYError(Float_t f) {fXYErr = f;}
  void SetZError(Float_t f) {fZErr = f;}
  void SetTransformer(AliL3Transform *transform) {fTransform = transform;}
  
  Int_t GetNumberOfClusters() {return fNClusters;}
  
  ClassDef(AliL3ClustFinderNew,1)

};

#endif
