#ifndef ALIL3_HOUGH_Eval
#define ALIL3_HOUGH_Eval

#include "AliL3RootTypes.h"

class TClonesArray;
class TObjArray;
class TH2F;
class TH1F;
class AliL3Transform;
class AliL3TrackArray;
class AliL3HoughTransformer;

class AliL3HoughEval : public TObject {
  
 private:
  
  Int_t fSlice;
  AliL3HoughTransformer *fHoughTransformer;
  AliL3Transform *fTransform; //!

 public:
  AliL3HoughEval(); 
  AliL3HoughEval(AliL3HoughTransformer *transformer);
  virtual ~AliL3HoughEval();

  TClonesArray *GetParticles(Char_t *rootfile);
  void LookInsideRoad(AliL3TrackArray *tracks,Bool_t remove=(Bool_t)false,TH2F *hist=0);
  void DisplaySlice(TH2F *hist);
  void CompareMC(Char_t *rootfile,AliL3TrackArray *merged_tracks,Float_t *eta);
  
  void SetTransformer(AliL3HoughTransformer *t) {fHoughTransformer=t;}

  ClassDef(AliL3HoughEval,1)

};

#endif
