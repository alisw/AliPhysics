// @(#) $Id$

#ifndef ALIL3_HOUGHCLUSTERTRANSFORMER
#define ALIL3_HOUGHCLUSTERTRANSFORMER

#include "AliL3RootTypes.h"
#include "AliL3HoughBaseTransformer.h"


class AliL3Histogram;
class AliL3SpacePointData;
class AliL3MemHandler;

class AliL3HoughClusterTransformer : public AliL3HoughBaseTransformer {
  
 private:

  AliL3Histogram **fParamSpace;   //!
  AliL3MemHandler *fMemHandler;   //!
  AliL3SpacePointData *fClusters; //!
  Int_t fNClusters;
#ifdef do_mc
  TrackIndex **fTrackID; //!
#endif
  void DeleteHistograms();
  
 public:
  AliL3HoughClusterTransformer(); 
  AliL3HoughClusterTransformer(Int_t slice,Int_t patch,Int_t n_eta_segments);
  virtual ~AliL3HoughClusterTransformer();
  
  void CreateHistograms(Int_t nxbin,Float_t ptmin,Int_t nybin,Float_t phimin,Float_t phimax);
  void CreateHistograms(Int_t nxbin,Float_t xmin,Float_t xmax,
			Int_t nybin,Float_t ymin,Float_t ymax);
  void FindClusters();
  void Reset();
  void TransformCircle();
  void TransformCircleC(Int_t *row_range,Int_t every);
  void TransformLine(Int_t */*rowrange*/=0,Float_t */*phirange*/=0){};
  
  Int_t GetEtaIndex(Double_t eta);
  AliL3Histogram *GetHistogram(Int_t eta_index);
  Double_t GetEta(Int_t eta_index,Int_t slice);
  Int_t GetTrackID(Int_t eta_index,Double_t kappa,Double_t psi);

  ClassDef(AliL3HoughClusterTransformer,1) //Normal Hough transformation class

};

#endif




