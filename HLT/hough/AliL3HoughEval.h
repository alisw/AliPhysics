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
  
  AliL3HoughTransformer *fHoughTransformer;
  AliL3Transform *fTransform; //!
  Int_t fNumOfPadsToLook;
  Int_t fNumOfRowsToMiss;
  Int_t *fMcTrackTable;
  
 public:
  AliL3HoughEval(); 
  AliL3HoughEval(AliL3HoughTransformer *transformer);
  virtual ~AliL3HoughEval();

  void DefineGoodParticles(Char_t *rootfile,Double_t pet);
  Bool_t LookInsideRoad(AliL3HoughTrack *track,Int_t eta_index,Bool_t remove=(Bool_t)kFALSE);
  void LookInsideRawRoad(AliL3TrackArray *tracks,Int_t eta_index,Bool_t remove=(Bool_t)kFALSE);
  void RemoveTrackFromImage(AliL3HoughTrack *track,Int_t eta_index);

  void DisplaySlice(TH2F *hist);
  void CompareMC(Char_t *rootfile,AliL3TrackArray *merged_tracks,Float_t *eta);
  Int_t *GetMcTrackTable() {return fMcTrackTable;}
  
  void SetNumOfPadsToLook(Int_t f) {fNumOfPadsToLook = f;}
  void SetNumOfRowsToMiss(Int_t f) {fNumOfRowsToMiss = f;}
  void SetTransformer(AliL3HoughTransformer *t) {fHoughTransformer=t;}
  


  ClassDef(AliL3HoughEval,1)

};

#endif
