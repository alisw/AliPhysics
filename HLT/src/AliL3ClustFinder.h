#ifndef ALIL3_ClustFinder
#define ALIL3_ClustFinder

#include "AliL3RootTypes.h"
class AliL3Transform;

#define FLAG_ONEPAD         1 //cluster had only one pad
#define FLAG_DOUBLE_PAD     2
#define FLAG_DOUBLE_T       1//4
#define PARAM1              1//4
#define MAX_C               112 //Maximum number of clusters on 1 pad.



struct resx 
{
  UInt_t t ;
  UInt_t pad ;
  UInt_t charge ;
  UInt_t flags ;
  Int_t  mean ;
  UInt_t falling ;
  UInt_t scharge ;
  
  Int_t trackID[3];
};
typedef struct resx resx;

class AliL3DigitRowData;
class AliL3SpacePointData;

class AliL3ClustFinder {

 private:
    
  AliL3Transform *fTransform; //!

  UInt_t fNDigitRowData;
  AliL3DigitRowData *fDigitRowData; //!

  
  AliL3SpacePointData *fSpacePointData; //!
  Int_t fNClusters;
  Int_t fMaxNClusters;
  Int_t fFirstRow;
  Int_t fLastRow;
  Int_t fCurrentSlice;
  Int_t fCurrentPatch;
  Int_t fCurrentRow;
  Bool_t fDeconvTime;
  Bool_t fDeconvPad;
  UInt_t fThreshold;

  Float_t fXYErr;
  Float_t fZErr;

 public:
  AliL3ClustFinder();
  AliL3ClustFinder(AliL3Transform *transform);
  virtual ~AliL3ClustFinder();

  void SetTransformer( AliL3Transform *transform );
  void InitSlice(Int_t slice,Int_t patch,Int_t firstrow, Int_t lastrow,Int_t maxpoints);
  void InitSlice(Int_t slice,Int_t patch,Int_t maxpoints);
  void Read(UInt_t ndigits,AliL3DigitRowData *ptr);
  void ProcessDigits();
  void ProcessRow(AliL3DigitRowData *tempPt);
  void SetOutputArray(AliL3SpacePointData *pt);

  void memcpy7(UInt_t *dst, UInt_t *src);
  void mstore(UInt_t *r,UInt_t av,UInt_t pad,UInt_t ch,UInt_t flags,UInt_t mean,Int_t *trackID);
  void WriteClusters(Int_t ncl,resx *r);
  
  void SetXYError(Float_t f) {fXYErr = f;}
  void SetZError(Float_t f) {fZErr = f;}
  
  //Getters
  Int_t GetNumberOfClusters() {return fNClusters;}

  ClassDef(AliL3ClustFinder,1) 
};

#endif
