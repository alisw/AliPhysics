#ifndef ALIL3_HOUGHTRANSFORMER
#define ALIL3_HOUGHTRANSFORMER

#include "AliL3RootTypes.h"
#include "AliL3HoughBaseTransformer.h"

class AliL3Histogram;

class AliL3HoughTransformer : public AliL3HoughBaseTransformer {
  
 private:
  
  AliL3Histogram **fParamSpace; //!
#ifdef do_mc
  TrackIndex **fTrackID; //!
#endif
  Bool_t fDoMC;

  void DeleteHistograms();

 public:
  AliL3HoughTransformer(); 
  AliL3HoughTransformer(Int_t slice,Int_t patch,Int_t n_eta_segments);
  virtual ~AliL3HoughTransformer();
  
  void CreateHistograms(Int_t nxbin,Double_t ptmin,Int_t nybin,Double_t phimin,Double_t phimax);
  void CreateHistograms(Int_t nxbin,Double_t xmin,Double_t xmax,
			Int_t nybin,Double_t ymin,Double_t ymax);
  void Reset();
  void TransformCircle();
  void TransformCircleC(Int_t row_range);
  void TransformLine();

  Int_t GetEtaIndex(Double_t eta);
  AliL3Histogram *GetHistogram(Int_t eta_index);
  Double_t GetEta(Int_t eta_index,Int_t slice);
  Int_t GetTrackID(Int_t eta_index,Double_t kappa,Double_t psi);
  
  //void Init(Int_t slice=0,Int_t patch=0,Int_t n_eta_segments=100);
  
  ClassDef(AliL3HoughTransformer,1) //Normal Hough transformation class

};

#endif




