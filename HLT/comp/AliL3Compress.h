#ifndef AliL3_Compress
#define AliL3_Compress

#include "AliL3RootTypes.h"
#include "AliL3Models.h"
#include "AliL3DigitData.h"

class AliL3TrackArray;

class AliL3Compress {
  
 private:
  AliL3TrackArray *fTracks; //!
  AliL3RandomDigitData *fDigits; //!
  AliL3RandomDigitData **fDPt; //!
  Int_t fNDigits;
  Int_t fNUsed;
  Int_t fMaxDigits;
  
  Int_t fNumPadBits;
  Int_t fNumTimeBits;
  Int_t fNumChargeBits;
  Int_t fNumShapeBits;
  Int_t fSlice;
  Int_t fPatch;
  
  void CreateDigitArray(Int_t maxnumber);
  void CreateDigits(Int_t row,Float_t pad,Float_t time,Int_t charge,Float_t ywidth,Float_t zwidth);
  void QSort(AliL3RandomDigitData **a, Int_t first, Int_t last);
  Int_t ComparePoints(Int_t row,UShort_t pad,UShort_t time);
  Int_t CompareDigits(AliL3RandomDigitData *a,AliL3RandomDigitData *b);

 public:
  AliL3Compress();
  AliL3Compress(Int_t slice,Int_t patch,Int_t pad=7,Int_t time=7,Int_t charge=10,Int_t shape=4);
  virtual ~AliL3Compress();
  
  void SetBitNumbers(Int_t pad,Int_t time,Int_t charge,Int_t shape);
  void WriteFile(AliL3TrackArray *tracks,Char_t *filename);
  void ReadFile(Char_t *filename);
  void CompressFile(Char_t *infile,Char_t *outfile);
  void ExpandFile(Char_t *infile,Char_t *outfile);
  void RestoreData(Char_t *uncompfile);
  void WriteRestoredData(Char_t *remainfile,Char_t *restoredfile);
  void PrintDigits();

  AliL3TrackArray *GetTracks() {return fTracks;}
  
  ClassDef(AliL3Compress,1) 

};

inline Int_t  AliL3Compress::ComparePoints(Int_t row,UShort_t pad,UShort_t time)
{
  
  if(fDPt[fNUsed]->fRow != row) return 0;
  
  if(fDPt[fNUsed]->fPad < pad) return 1;
  if(fDPt[fNUsed]->fPad == pad && fDPt[fNUsed]->fTime < time) return 1;
  
  return 0;

}

inline Int_t AliL3Compress::CompareDigits(AliL3RandomDigitData *a,AliL3RandomDigitData *b)
{
  if(a->fRow < b->fRow) return -1;
    
  if(a->fPad==b->fPad && a->fTime == b->fTime && a->fRow == b->fRow) return 0;
  
  if(a->fPad<b->fPad && a->fRow == b->fRow) return -1;
  if(a->fPad==b->fPad && a->fTime<b->fTime && a->fRow == b->fRow) return -1;
  
  return 1;
}

#endif
