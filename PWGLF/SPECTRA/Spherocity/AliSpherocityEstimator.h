#ifndef ALISPHEROCITYESTIMATOR__H
#define ALISPHEROCITYESTIMATOR__H
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliVEvent.h"
#include "AliVTrack.h"
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
  Double_t GetSpherocity(AliVEvent *inevent);
  AliESDtrackCuts *GetESDtrackCuts() { return fSphTrackCuts; };
  void SetESDtrackCuts(AliESDtrackCuts *newcuts);
  void SetAODFilterBit(Int_t newval=128);
 private:
  AliESDtrackCuts *fSphTrackCuts; //! No need to store
  Int_t fAODFilterBit; //! no need to store either
  Int_t fMinMulti;
  Int_t fTrackMulti;
  Double_t fMinPt;
  Bool_t fOnAODs; //! AOD flag
  ClassDef(AliSpherocityEstimator,1);
};
#endif
