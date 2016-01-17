/// \class AliAnalysisC2Utils
/// \brief Utility functions used in the C2 correlation analysis
///
/// \author Christian Bourjau <cbourjau@cern.ch>, University of Copenhagen, Denmark
/// \date Mar 22, 2016

#ifndef AliAnalysisC2Utils_cxx
#define AliAnalysisC2Utils_cxx

#include "TObject.h"
#include "THn.h"

class TAxis;
class AliAODTrack;

class AliAnalysisC2Utils : public TObject {
 public:
  static Double_t Mod(Double_t x, Double_t y);
  static Double_t WrapAngle(Double_t angle, TAxis* ax);
  static Int_t ExcludeOverUnderFlow(THn *h);
  static void CopyAxisFromHist(Int_t from_axis, THn *from_hist, Int_t to_axis, THn *to_hist);
  static void TransformPoints(Double_t *points_in, Double_t *points_out);
  static THnS* CreateTransformedHist(THn *h);
  static void GetDCA(Double_t& DCAtang, Double_t& DCAlong, AliAODTrack* AODt);

 private:
  ClassDef(AliAnalysisC2Utils,1);
};

#endif
