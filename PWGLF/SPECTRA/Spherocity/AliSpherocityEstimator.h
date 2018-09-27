#ifndef ALISPHEROCITYESTIMATOR__H
#define ALISPHEROCITYESTIMATOR__H
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliESDEvent.h"
#include "TNamed.h"

class AliSpherocityEstimator: public TNamed
{
 public:
  AliSpherocityEstimator(); //Default constructor
  AliSpherocityEstimator(const AliSpherocityEstimator &source); //copy consturctor. Implement if found really necessary
  AliSpherocityEstimator &operator=(const AliSpherocityEstimator &source); //operator=. Implement if found really necessary
  ~AliSpherocityEstimator(); //Destructor
  void SetMinMulti(Int_t nTracks) { if(nTracks>0) fMinMulti=nTracks; };
  Int_t GetMinMulti() { return fMinMulti; };
  void SetMinPt(Double_t MinPt) { if(MinPt>=0) fMinPt=MinPt; };
  Double_t GetMinPt() { return fMinPt; };
  Int_t GetTrackMulti() { return fTrackMulti; };
  Double_t GetSpherocity(AliESDEvent *inevent);
  AliESDtrackCuts *GetESDtrackCuts() { return fSphTrackCuts; };
  void SetESDtrackCuts(AliESDtrackCuts *newcuts);
 private:
  AliESDtrackCuts *fSphTrackCuts; //! No need to store
  Int_t fMinMulti;
  Int_t fTrackMulti;
  Double_t fMinPt;
  ClassDef(AliSpherocityEstimator,1);
};
#endif
