// @(#) $Id$

#ifndef ALIL3HOUGHTRANSFORMERNEW_H
#define ALIL3HOUGHTRANSFORMERNEW_H

#include "AliHLTRootTypes.h"
#include "AliHLTHoughTransformer.h"

#include <TH3.h>

class AliHLTTrackArray;
class AliHLTHoughTrack;
 
class AliHLTHoughTransformerNew : public AliHLTHoughTransformer {

 public:
  AliHLTHoughTransformerNew(); 
  AliHLTHoughTransformerNew(Int_t slice,Int_t patch,Int_t netasegments);
  virtual ~AliHLTHoughTransformerNew();
  
  void Reset();
  void CreateHistograms(Float_t ptmin,Float_t ptmax,Float_t pres,Int_t nybin,Float_t psi) {
    AliHLTHoughTransformer::CreateHistograms(ptmin,ptmax,pres,nybin,psi);
  }
  void CreateHistograms(Int_t nxbin,Float_t ptmin,Int_t nybin,Float_t phimin,Float_t phimax) {
    AliHLTHoughTransformer::CreateHistograms(nxbin,ptmin,nybin,phimin,phimax);
  }
  void CreateHistograms(Int_t nxbin,Float_t xmin,Float_t xmax,Int_t nybin,Float_t ymin,Float_t ymax) {
    AliHLTHoughTransformer::CreateHistograms(nxbin,xmin,xmax,nybin,ymin,ymax);
  }
  void CreateHistograms(Int_t nxbins,Float_t xlow,Float_t xup,
			Int_t nybins,Float_t ylow,Float_t yup,
			Int_t nzbins,Float_t zlow,Float_t zup);
  void TransformLine(Int_t *rowrange,Float_t *phirange);
  void TransformLineC(Int_t *rowrange,Float_t *phirange);
  
  TH3 *GetHistogram() {return fParamSpace3D;}
  AliHLTHistogram *GetHistogram(Int_t etaindex){
    return AliHLTHoughTransformer::GetHistogram(etaindex);
  }
  
 private:
  
  TH3 *fParamSpace3D;//Histogram containing the hough space
  
  ClassDef(AliHLTHoughTransformerNew,1) //Normal Hough transformation class

};

typedef AliHLTHoughTransformerNew AliL3HoughTransformerNew; // for backward comaptibility

#endif




