#ifndef AliL3_ClusterFitter
#define AliL3_ClusterFitter

#include "AliL3RootTypes.h"
#include "AliL3Modeller.h"

class AliL3ModelTrack;

class AliL3ClusterFitter : public AliL3Modeller {
  
 private:
  Int_t fPadFitRange;
  Int_t fTimeFitRange;
  Int_t fNmaxOverlaps;
  Float_t fChiSqMax;
  
  void FitCluster(AliL3ModelTrack *track);
  void FitClusters(AliL3ModelTrack *track);
  
 public:
  AliL3ClusterFitter();
  virtual ~AliL3ClusterFitter();
  
  void FindClusters();

  void SetFitRange(Int_t p,Int_t t) {fPadFitRange=p; fTimeFitRange=t;}
  void SetNmaxOverlaps(Int_t i) {fNmaxOverlaps=i;}

  ClassDef(AliL3ClusterFitter,1) 

};

#endif
