#ifndef TOPOLOGY_H
#define TOPOLOGY_H
#include "TObject.h"

class TH2F;
class TBits;
class AliITSUClusterPix;

class Topology :public TObject {

 public:
  enum SortMode_t{kFrequency=0, kHashes=1};//fMode
  enum FitIndex_t{kDeltaZmean=0, kDeltaZmeanErr=1, kDeltaXmean=2, kDeltaXmeanErr=3, kDeltaZsigma=4, kDeltaZsigmaErr=5,
		  kDeltaXsigma=6, kDeltaXsigmaErr=7, kChi2z=8, kChi2x=9, kNDFx=10, kNDFz=11, kFitLength=12}; //position in fArrFit

  Topology();
  Topology(const AliITSUClusterPix &cluster);
  Topology(const TBits &top2copy);//UniqueID of the argument must already have been set
  Topology(const Topology &topo);
  virtual ~Topology();

  Bool_t IsEqual(const TObject* obj) const;
  Bool_t IsSortable() const {return kTRUE;}
  Int_t Compare(const TObject* obj) const;
  
  TBits GetPattern() const {return fPattern;};
  UInt_t GetUniqueID() const {return fUID;}
  Int_t GetRowSpan() const {return fRs;}
  Int_t GetColumnSpan() const {return fCs;}
  Int_t GetWordLength() const {return fWordLength;}
  Int_t GetHash() const {return fHash;}
  Float_t GetFreq() const {return fFreq;}
  Int_t GetCounts() const {return fCounts;}
  Int_t GetGroupID() const {return fGroupID;}
  Int_t GetFiredPixels() const {return fFiredPixels;}
  Float_t GetxCOGPix() const {return fxCOGPix;}
  Float_t GetzCOGPix() const {return fzCOGPix;}
  Float_t GetxCOGshift() const {return fxCOGshift;}
  Float_t GetzCOGshift() const {return fzCOGshift;}
  Int_t GetMode() const {return fMode;}
  TH2F* GetHxA() const {return fHxA;}
  TH2F* GetHxB() const {return fHxB;}
  TH2F* GetHzA() const {return fHzA;}
  TH2F* GetHzB() const {return fHzB;}
  Float_t GetFitStuff(Int_t ind) const {return fArrFit[ind];}
  Int_t GetFlag() const {return fFlag;}
  Int_t GetPartialTop() const {return fPartialTop;}
  Int_t GetPattID() const {return fPattID;}

  static UInt_t FuncMurmurHash2(const void * key, Int_t len);
  static void Word2Top(const UChar_t* Word, TBits &top);
  static void Top2Word(const TBits &top, UChar_t* Word);
  static Int_t Top2Int(const TBits &top);
  static std::ostream& printTop(TBits top, std::ostream &out);
  //Word: 1st byte = row span; 2nd = column span; others: pattern.
  //The length must be the minimum possible.

  void SetGroupID(Int_t num){fGroupID=num;}
  void SetFreq(Float_t num){fFreq=num;}
  void SetFitStuff(Float_t value, Int_t ind) {fArrFit[ind]=value;}
  void SetFlag(Int_t num) {fFlag=num;}
  void IncreaseCounts(){fCounts++;}
  void SetPattID(Int_t num) {fPattID=num;}
  static void SetMode(Int_t mode) {fMode=mode;}
  void DeleteHistos();
  void SetHxA(TH2F* ptr) {fHxA=ptr;}
  void SetHzA(TH2F* ptr) {fHzA=ptr;}
  void SetHxB(TH2F* ptr) {fHxB=ptr;}
  void SetHzB(TH2F* ptr) {fHzB=ptr;}

 private:

  TBits fPattern;
  UInt_t fUID; //Unique ID
  Int_t fRs; //row span
  Int_t fCs; //column span
  Int_t fFiredPixels;
  Float_t fxCOGPix;
  Float_t fzCOGPix;
  Float_t fxCOGshift;
  Float_t fzCOGshift;
  Int_t fHash;
  Int_t fWordLength;
  UChar_t* fWord;//[fWordLength]
  Int_t fCounts;
  Float_t fFreq;
  Int_t fGroupID;
  static Int_t fMode; //DEFAULT kHashes
  TH2F* fHxA;
  TH2F* fHxB;
  TH2F* fHzA;
  TH2F* fHzB;
  Float_t fArrFit[kFitLength];
  Int_t fFlag;
  Int_t fPartialTop;
  Int_t fPattID;
 
  ClassDef(Topology,1)
    
};

#endif
