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
  static Double_t Wrap02pi(Double_t angle);
  static Int_t ExcludeOverUnderFlow(THn *h);
  static void CopyAxisFromHist(Int_t from_axis, THn *from_hist, Int_t to_axis, THn *to_hist);
  static void TransformPoints(Double_t *points_in, Double_t *points_out);
  static THnS* CreateTransformedHist(THn *h);
  static void GetDCA(Double_t& DCAtang, Double_t& DCAlong, AliAODTrack* AODt);
  // Create scalar for the two given pt bins
  //
  // Since the analysis requires that pt1 < pt2, it would be a wast of
  // space to maintain both dimensions in the pair histogram. Ie, all
  // bins where pt1 > pt2 would be empty.  Hence, this function
  // computes a unique bin index for the two pt bins in the interval
  // [0, \sum_{i=1}^{pt2} i]
  //
  // This binning requires that the pt1 and pt2 intervals are identical.
  //
  // The formular used is:
  // idx = \sum_{i=1}^{pt2} i + (pt1 - pt2) - 1
  // Returns -1 if either of the bins is 0 (ie underflow bins)
  static Int_t ComputePtPairBin(Int_t pt1Bin, Int_t pt2Bin);

  // Function to check if the current event fits a given trigger (ig AliVEvent::kMB)
  static Bool_t EventFitsTrigger(UInt_t trig);
 private:
  ClassDef(AliAnalysisC2Utils,1);
};

#endif
