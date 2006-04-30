// @(#) $Id$

#ifndef ALIL3HOUGHCLUSTERTRANSFORMER_H
#define ALIL3HOUGHCLUSTERTRANSFORMER_H

#include "AliL3RootTypes.h"
#include "AliL3HoughBaseTransformer.h"


class AliL3Histogram;
class AliL3SpacePointData;
class AliL3MemHandler;

class AliL3HoughClusterTransformer : public AliL3HoughBaseTransformer {
  
 public:
  AliL3HoughClusterTransformer(); 
  AliL3HoughClusterTransformer(Int_t slice,Int_t patch,Int_t netasegments);
  virtual ~AliL3HoughClusterTransformer();
  
  void CreateHistograms(Float_t ptmin,Float_t ptmax,Float_t pres,Int_t nybin,Float_t psi) {
    AliL3HoughBaseTransformer::CreateHistograms(ptmin,ptmax,pres,nybin,psi);
  }
  void CreateHistograms(Int_t nxbin,Float_t ptmin,Int_t nybin,Float_t phimin,Float_t phimax);
  void CreateHistograms(Int_t nxbin,Float_t xmin,Float_t xmax,
			Int_t nybin,Float_t ymin,Float_t ymax);
  void FindClusters();
  void Reset();
  void TransformCircle();
  void TransformCircle(Int_t *row_range,Int_t every) {
    AliL3HoughBaseTransformer::TransformCircle(row_range,every);
  }
  void TransformCircleC(Int_t *rowrange,Int_t every);
  void TransformLine(Int_t */*rowrange*/=0,Float_t */*phirange*/=0){};
  
  Int_t GetEtaIndex(Double_t eta) const;
  AliL3Histogram *GetHistogram(Int_t etaindex);
  Double_t GetEta(Int_t etaindex,Int_t slice) const;
  Int_t GetTrackID(Int_t etaindex,Double_t kappa,Double_t psi) const;
  
 private:

  AliL3Histogram **fParamSpace;   //!
  AliL3MemHandler *fMemHandler;   //!
  AliL3SpacePointData *fClusters; //!
  Int_t fNClusters;//Number of clusters
#ifdef do_mc
  AliL3TrackIndex **fTrackID; //!
#endif
  void DeleteHistograms();

  ClassDef(AliL3HoughClusterTransformer,1) //Normal Hough transformation class

};

#endif




